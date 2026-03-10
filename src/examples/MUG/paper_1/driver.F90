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
USE oft_gs, ONLY: gs_eq, gs_update_bounds, gs_test_bounds, compute_bcmat, gs_setup_walls, gs_get_qprof
USE oft_gs_util, ONLY: gs_profile_load
USE gs_xmhd_v8
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
INTEGER(i4) :: nsteps = 1
INTEGER(i4) :: rst_freq = 1
INTEGER(i4) :: ndims, nl_its, l_its, nretry
INTEGER(i4) :: npoints
integer(i4), allocatable, dimension(:) :: dim_sizes
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: dt = 0.04336664911469267 
REAL(r8) :: t = 0.d0
REAL (r8):: ip_ratio_target = 0.205
REAL (r8):: ip_target = 7.87E6
REAL(r8), allocatable, dimension(:) :: psi_eq, psi_pert, psi_total, eta_reg,curr_reg, areas
REAL (r8):: coords(3), psi(1), q(1)
LOGICAL :: pm=.FALSE.
LOGICAL :: success
CHARACTER(LEN=25) :: filename_eq = 'paper_eq_0115.h5' !< Name of input file for mesh, fix later for variable length
CHARACTER(LEN=25) :: filename_pert= 'paper_pert_0115.h5' !< Name of input file for mesh, fix later for variable length
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
psi_total = psi_eq - 0.1*psi_pert
! psi_total = psi_eq

!---------------------------------------------------------------------------
! Now, need to setup a tokamaker object
!---------------------------------------------------------------------------
CALL equil%setup(ML_blagrange_2)
equil%region_info%nnonaxi = 0
ALLOCATE(equil%region_info%reg_map(equil%fe_rep%mesh%nreg))
equil%region_info%reg_map=0
CALL gs_setup_walls(equil)
CALL equil%init()
CALL equil%psi%restore_local(psi_total)
CALL gs_update_bounds(equil, track_opoint = .TRUE.)
write(*,*) equil%plasma_bounds
equil%itor_target=ip_target*mu0
equil%ip_ratio_target=ip_ratio_target
equil%pnorm = 0.28658184156588085
equil%alam = 2.596639717247778
tmp_str = 'tokamaker_f.prof'
CALL gs_profile_load(tmp_str,equil%I)
tmp_str = 'tokamaker_p.prof'
CALL gs_profile_load(tmp_str,equil%P)
equil%I%plasma_bounds=equil%plasma_bounds
equil%P%plasma_bounds=equil%plasma_bounds
equil%ncoils = 7
equil%ncoil_regs = 7
ALLOCATE(equil%coil_nturns(equil%mesh%nreg,equil%ncoils))
equil%coil_nturns = 0
DO j=1, equil%ncoils
  equil%coil_nturns(j + 8, j) = 1
END DO
equil%vcontrol_val = 0.d0
equil%coil_vcont = 0.d0
equil%ncond_regs = 5
ALLOCATE(equil%cond_regions(equil%ncond_regs))
equil%cond_regions(1)%id = 4
equil%cond_regions(1)%eta = 6.9d-7/mu0
equil%cond_regions(2)%id = 5
equil%cond_regions(2)%eta = 7.d-7/mu0
equil%cond_regions(3)%id = 6
equil%cond_regions(3)%eta = 7.d-7/mu0
equil%cond_regions(4)%id = 7
equil%cond_regions(4)%eta = 6.9d-7/mu0
equil%cond_regions(5)%id = 8
equil%cond_regions(5)%eta = 6.9d-7/mu0
ALLOCATE(equil%coil_regions(equil%ncoil_regs))
ALLOCATE(equil%coil_currs(equil%ncoil_regs))
equil%coil_currs = [-10004540.249054534, 8131096.462417697, 8130193.40495448, -2752265.1799410507, -2755148.7923723585, -1429157.477860515,  -1424955.0239338675]*mu0
ALLOCATE(areas(equil%ncoil_regs))
areas = [1.8,0.25, 0.25, 0.25, 0.25, 0.25, 0.25 ]
equil%coil_currs = equil%coil_currs/areas
equil%mode = 0
equil%I%f_offset = 36.d0
DO j=1, equil%ncoils
  equil%coil_regions(j)%id = 8 + j
END DO
CALL compute_bcmat(equil)
!---------------------------------------------------------------------------
! Setup time-dependent solver
!---------------------------------------------------------------------------
gs_td%nsteps = 800
! gs_td%dt = dt/5.d0
gs_td%dt = dt
gs_td%lin_tol = 1.d-11
gs_td%nl_tol = 1.d-9
gs_td%eq => equil
gs_td%pm = pm

gs_td%nu = 1.d-3
gs_td%rho = 9806.d0
gs_td%B_0 = 0.d0


ALLOCATE(eta_reg(equil%mesh%nreg))
ALLOCATE(gs_td%eta_t(equil%mesh%nreg))
ALLOCATE(gs_td%eta_p(equil%mesh%nreg))
! ALLOCATE(gs_td%eta(equil%mesh%nreg))
eta_reg = 1.d-2/mu0
eta_reg(4) = 6.9d-7/mu0
eta_reg(5) = 7.0d-7/mu0
eta_reg(6) = 7.0d-7/mu0
eta_reg(5) = 1.14d-6/mu0
eta_reg(6) = 1.14d-6/mu0
eta_reg(7) = 6.9d-7/mu0
eta_reg(8) = 6.9d-7/mu0
gs_td%eta_t = eta_reg
gs_td%eta_p = eta_reg


ALLOCATE(curr_reg(equil%mesh%nreg))
curr_reg = -1.d0
DO j=1, equil%mesh%nreg
  IF (j >=9) curr_reg(j) = equil%coil_currs(j-8)
END DO
gs_td%curr = curr_reg
ALLOCATE(gs_td%region_flag(equil%mesh%nreg))
DO j=1, SIZE(gs_td%region_flag)
  IF (j==1) gs_td%region_flag(j) = 5
  IF (j==2 .OR. j==3) gs_td%region_flag(j) = 2
  if (j==4) gs_td%region_flag(j) = 3
  if (j==5) gs_td%region_flag(j) = 3
  if (j==6) gs_td%region_flag(j) = 3
  if (j==7) gs_td%region_flag(j) = 3
  if (j==8) gs_td%region_flag(j) = 3
  IF (j >=9) gs_td%region_flag(j) = 4
END DO

gs_td%evolve_F = .TRUE.
CALL gs_td%setup(mg_mesh, mg_mesh_1)
!Set initial values for fields
field_init%mesh=>mg_mesh%smesh
field_init%func=>const_init
NULLIFY(tmp_arr)
CALL equil%psi%get_local(tmp_arr)
CALL gs_td%u%restore_local(tmp_arr,6)
tmp_arr = equil%I%f_offset
CALL gs_td%u%restore_local(tmp_arr,5)

! psi = 0.00d0
! q = 1.d0
! gs_td%dt = .015011d0/50.d0
! DO j=1, 50
!   write(*,*) j
!   write(*,*) gs_td%eq%itor_target 
!   CALL gs_td%add_timestep(gs_td%dt)
!   gs_td%eq%itor_target = gs_td%eq%itor_target - 0.02*ip_target*mu0
! END DO
! equil%Ip_ratio_target = equil%Ip_ratio_target*1.03
! CALL gs_td%add_timestep(gs_td%dt)
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