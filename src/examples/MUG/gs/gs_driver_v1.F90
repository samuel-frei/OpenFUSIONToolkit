PROGRAM gs_driver
!---Runtime
USE oft_base
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE oft_blag_operators, ONLY: oft_blag_zerob, oft_blag_getmop, oft_blag_project
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: elec_charge, proton_mass, mu0
USE gs_xmhd
USE oft_io, ONLY: hdf5_field_get_sizes, hdf5_read, hdf5_field_exist
USE oft_gs, ONLY: gs_eq, gs_update_bounds, gs_test_bounds
USE oft_gs_util, ONLY: gs_profile_load

IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr, i, j
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_gs_xmhd_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_blag_zerob), TARGET :: blag_zerob ! setting boundary vals to zero
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
TYPE(gs_eq), TARGET :: equil
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 1000
INTEGER(i4) :: rst_freq = 1
INTEGER(i4) :: ndims
INTEGER(i4) :: npoints
integer(i4), allocatable, dimension(:) :: dim_sizes
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: dt = 1.d-9
REAL (r8):: ip_ratio_target = 1.0
REAL (r8):: ip_target = 0.75E6
REAL(r8), allocatable, dimension(:) :: psi_eq, psi_pert, psi_total, eta_reg, areas
REAL (r8):: coords(3)
LOGICAL :: pm=.FALSE.
LOGICAL :: success
CHARACTER(LEN=25) :: filename_eq = 'nsf_eq.h5' !< Name of input file for mesh, fix later for variable length
CHARACTER(LEN=25) :: filename_pert= 'nsf_perturbation.h5' !< Name of input file for mesh, fix later for variable length
CHARACTER(LEN=25) :: tmp_str

!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
! Set region flags for appropriate boundary conditions
! 1 -> plasma
! 2-> vacuum
! 3 -> solid conductor
! 4 -> superconductor
ALLOCATE(mhd_sim%region_flag(mg_mesh%smesh%nreg))
DO j=1, SIZE(mhd_sim%region_flag)
  IF (j==1) mhd_sim%region_flag(j) = 1
  IF (j==2) mhd_sim%region_flag(j) = 2
  if (j==3) mhd_sim%region_flag(j) = 3
  IF (j >=4) mhd_sim%region_flag(j) = 4
END DO
CALL mhd_sim%setup(mg_mesh,order)
!---------------------------------------------------------------------------
! Read equilibrium from file
!---------------------------------------------------------------------------
CALL hdf5_field_get_sizes(TRIM(filename_eq),"tokamaker/PSI",ndims,dim_sizes)
npoints = dim_sizes(1)
ALLOCATE(psi_eq(npoints))
CALL hdf5_read(psi_eq,TRIM(filename_eq),"tokamaker/PSI",success)

CALL hdf5_field_get_sizes(TRIM(filename_pert),"tokamaker/PSI",ndims,dim_sizes)
npoints = dim_sizes(1)
ALLOCATE(psi_pert(npoints))
ALLOCATE(psi_total(npoints))
CALL hdf5_read(psi_pert,TRIM(filename_pert),"tokamaker/PSI",success)
psi_total = psi_eq + 0.1*psi_pert

CALL hdf5_field_get_sizes(TRIM(filename_eq),"region_info/ETA",ndims,dim_sizes)
npoints = dim_sizes(1)
ALLOCATE(eta_reg(npoints))
CALL hdf5_read(eta_reg,TRIM(filename_eq),"region_info/ETA",success)
!---------------------------------------------------------------------------
! Now, need to setup a tokamaker object
!---------------------------------------------------------------------------
CALL equil%setup(ML_oft_blagrange)
CALL equil%init()
CALL equil%psi%restore_local(psi_total)
CALL gs_update_bounds(equil, track_opoint = .TRUE.)
equil%itor_target=ip_target*mu0
equil%ip_ratio_target=ip_ratio_target
equil%pnorm = 2.865970100459969
equil%alam = 6.046117947969407
tmp_str = 'tokamaker_f.prof'
CALL gs_profile_load(tmp_str,equil%I)
tmp_str = 'tokamaker_p.prof'
CALL gs_profile_load(tmp_str,equil%P)
equil%I%plasma_bounds=equil%plasma_bounds
equil%P%plasma_bounds=equil%plasma_bounds
equil%ncoil_regs = 11
ALLOCATE(equil%coil_currs(equil%ncoil_regs))
equil%coil_currs = [-3001104.351425276, -3000075.909961613, -3002206.4249482094, -444076.14401118725, -443334.8458864972, 36268.199840381465, &
36086.31721276968,63457.412912549225, 63709.70927346871, -714203.0889287366, -714583.222811343 ]* mu0
ALLOCATE(areas(equil%ncoil_regs))
areas = [0.0315,0.0585, 0.0315, 0.015625, 0.015625, 0.030625, 0.030625, 0.0225, 0.0225, 0.030625, 0.030625 ]
equil%coil_currs = equil%coil_currs/areas
!---------------------------------------------------------------------------
! Plot initial fields
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(u,v,mop,vec_vals) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)
!---Project psi initial condition onto scalar Lagrange basis
CALL mesh%save_vertex_scalar(psi_total,mhd_sim%xdmf_plot,'psi0')
CALL mhd_sim%u%restore_local(psi_total,1)

!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)

DEALLOCATE(psi_eq)
DEALLOCATE(psi_pert)
DEALLOCATE(psi_total)

!---------------------------------------------------------------------------
! Set simulation settings and run
!---------------------------------------------------------------------------
mhd_sim%dt=dt
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%eq => equil
oft_env%pm=pm

eta_reg = -1.d0
eta_reg(3) = 1.E-6
mhd_sim%eta = eta_reg

CALL mhd_sim%run_simulation()

!---Finalize enviroment
CALL oft_finalize
CONTAINS
END PROGRAM gs_driver