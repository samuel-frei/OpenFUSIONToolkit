PROGRAM gs_driver_full
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
USE oft_io, ONLY: hdf5_field_get_sizes, hdf5_read, hdf5_field_exist
USE oft_gs, ONLY: gs_eq, gs_update_bounds, gs_test_bounds, compute_bcmat
USE oft_gs_util, ONLY: gs_profile_load
USE gs_xmhd
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE fem_base, ONLY: oft_ml_fem_type

IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr, i, j
TYPE(multigrid_mesh) :: mg_mesh, mg_mesh_1
TYPE(oft_ml_fem_type), TARGET :: ML_blagrange_1, ML_blagrange_2
CLASS(oft_scalar_bfem), POINTER :: blagrange_1 => NULL()
CLASS(oft_scalar_bfem), POINTER :: blagrange_2 => NULL()
TYPE(oft_blag_zerob), TARGET :: blag_zerob ! setting boundary vals to zero
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_solver), POINTER :: minv_2 => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: mop_2 => NULL()
CLASS(oft_vector), POINTER :: u,v, u_2, v_2
TYPE(gs_eq), TARGET :: equil
TYPE(oft_gs_xmhd_sim), TARGET :: gs_td
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 5
INTEGER(i4) :: rst_freq = 1
INTEGER(i4) :: ndims, nl_its, l_its, nretry
INTEGER(i4) :: npoints
integer(i4), allocatable, dimension(:) :: dim_sizes
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: dt = 4.346d-4
REAL(r8) :: t = 0.d0
REAL (r8):: ip_ratio_target = 1.0
REAL (r8):: ip_target = 0.75E6
REAL(r8), allocatable, dimension(:) :: psi_eq, psi_pert, psi_total, eta_reg,curr_reg, areas
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
CALL multigrid_construct_surf(mg_mesh_1)
! order = 1
! CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_blagrange_1,minlev=-1)
! IF(.NOT.oft_2D_lagrange_cast(blagrange_1,ML_blagrange_1%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
order = 2
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_blagrange_2,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(blagrange_2,ML_blagrange_2%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
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

!---------------------------------------------------------------------------
! Now, need to setup a tokamaker object
!---------------------------------------------------------------------------
CALL equil%setup(ML_blagrange_2)
equil%region_info%nnonaxi = 0
ALLOCATE(equil%region_info%reg_map(equil%fe_rep%mesh%nreg))
equil%region_info%reg_map=0
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
equil%ncoils = 11
equil%ncoil_regs = 11
ALLOCATE(equil%coil_nturns(equil%mesh%nreg,equil%ncoils))
equil%coil_nturns = 0
DO j=1, equil%ncoils
  equil%coil_nturns(j + 3, j) = 1
END DO
equil%vcontrol_val = 0.d0
equil%coil_vcont = 0.d0
equil%ncond_regs = 1
ALLOCATE(equil%cond_regions(equil%ncond_regs))
equil%cond_regions(1)%id = 3
equil%cond_regions(1)%eta = 1.d-6/mu0
ALLOCATE(equil%coil_regions(equil%ncoil_regs))
ALLOCATE(equil%coil_currs(equil%ncoil_regs))
equil%coil_currs = [-3001104.351425276, -3000075.909961613, -3002206.4249482094, -444076.14401118725, -443334.8458864972, 36268.199840381465, &
36086.31721276968,63457.412912549225, 63709.70927346871, -714203.0889287366, -714583.222811343 ]*mu0
ALLOCATE(areas(equil%ncoil_regs))
areas = [0.0315,0.0585, 0.0315, 0.015625, 0.015625, 0.030625, 0.030625, 0.0225, 0.0225, 0.030625, 0.030625 ]
equil%coil_currs = equil%coil_currs/areas
DO j=1, equil%ncoils
  equil%coil_regions(j)%id = 3 + j
END DO
CALL compute_bcmat(equil)
!---------------------------------------------------------------------------
! Setup time-dependent solver
!---------------------------------------------------------------------------
gs_td%nsteps = 5
gs_td%dt = dt/10.d0
gs_td%lin_tol = 1.d-11
gs_td%nl_tol = 1.d-09
gs_td%eq => equil
gs_td%pm = pm

gs_td%nu = 1.d-3
gs_td%rho = 9806.d0
gs_td%B_0 = 0.d0

ALLOCATE(eta_reg(equil%mesh%nreg))
eta_reg = -1.d0
eta_reg(3) = 0.7d-6/mu0
gs_td%eta = eta_reg
ALLOCATE(curr_reg(equil%mesh%nreg))
curr_reg = -1.d0
DO j=1, equil%mesh%nreg
  IF (j >=4) curr_reg(j) = equil%coil_currs(j-3)
END DO
gs_td%curr = curr_reg
ALLOCATE(gs_td%region_flag(equil%mesh%nreg))
DO j=1, SIZE(gs_td%region_flag)
  IF (j==1) gs_td%region_flag(j) = 5
  IF (j==2) gs_td%region_flag(j) = 2
  if (j==3) gs_td%region_flag(j) = 1
  IF (j >=4) gs_td%region_flag(j) = 4
END DO
CALL gs_td%setup(mg_mesh, mg_mesh_1)
!Set initial values for fields
field_init%mesh=>mg_mesh%smesh
field_init%func=>const_init
NULLIFY(tmp_arr)
CALL equil%psi%get_local(tmp_arr)
CALL gs_td%u%restore_local(tmp_arr,6)


CALL gs_td%run_simulation()
! !---Finalize enviroment
! CALL oft_finalize
CONTAINS

SUBROUTINE const_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = 1.0
END SUBROUTINE const_init

END PROGRAM gs_driver_full