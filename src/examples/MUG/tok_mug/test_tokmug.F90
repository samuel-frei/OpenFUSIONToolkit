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
USE oft_gs, ONLY: gs_equil, gs_update_bounds, gs_test_bounds, compute_bcmat, gs_setup_walls, gs_get_qprof, gs_factory
USE oft_gs_util, ONLY: gs_profile_load
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE fem_base, ONLY: oft_ml_fem_type
USE oft_blanket_td

IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr, i, j
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_blagrange
CLASS(oft_scalar_bfem), POINTER :: blagrange
TYPE(oft_blag_zerob), TARGET :: blag_zerob ! setting boundary vals to zero
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_solver), POINTER :: minv_2 => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: mop_2 => NULL()
CLASS(oft_vector), POINTER :: u,v, u_2, v_2
TYPE(gs_equil), TARGET :: equil
TYPE(gs_factory), TARGET :: machine
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
REAL(r8), allocatable, dimension(:) :: psi_eq, psi_pert, psi_total, eta_reg,curr_reg, areas, dens_reg, visc_reg
REAL(r8) :: lin_tol = 1.d-11
REAL(r8) :: nl_tol = 1.d-9
REAL (r8):: coords(3), psi(1), q(1)
LOGICAL :: pm=.TRUE.
LOGICAL :: success
CHARACTER(LEN=25) :: filename_eq = 'equilibrium.h5' !< Name of input file for mesh, fix later for variable length
CHARACTER(LEN=25) :: filename_pert= 'paper_pert_0115.h5' !< Name of input file for mesh, fix later for variable length
CHARACTER(LEN=25) :: tmp_str
TYPE(oft_blanket_td_sim):: b_sim 
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)

order = 2
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_blagrange,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(blagrange,ML_blagrange%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
!---------------------------------------------------------------------------
! Read equilibrium from file
!---------------------------------------------------------------------------
CALL hdf5_field_get_sizes(TRIM(filename_eq),"tokamaker/PSI",ndims,dim_sizes)
npoints = dim_sizes(1)
ALLOCATE(psi_eq(npoints))
CALL hdf5_read(psi_eq,TRIM(filename_eq),"tokamaker/PSI",success)

psi_total = psi_eq


!---------------------------------------------------------------------------
! Now, need to setup a tokamaker device
!---------------------------------------------------------------------------
CALL machine%setup(ML_blagrange)
machine%region_info%nnonaxi = 0
ALLOCATE(machine%region_info%reg_map(machine%fe_rep%mesh%nreg))
machine%region_info%reg_map=0
machine%free = .TRUE.

CALL gs_setup_walls(machine)
CALL machine%init()
machine%ncoils = 7
machine%ncoil_regs = 7
ALLOCATE(machine%coil_nturns(machine%fe_rep%mesh%nreg,machine%ncoils))
machine%coil_nturns = 0
DO j=1, machine%ncoils
  machine%coil_nturns(j + 8, j) = 1
END DO

machine%coil_vcont = 0.d0
machine%ncond_regs = 5
ALLOCATE(machine%cond_regions(machine%ncond_regs))
machine%cond_regions(1)%id = 4
machine%cond_regions(1)%eta = 6.9d-7/mu0
machine%cond_regions(2)%id = 5
machine%cond_regions(2)%eta = 1.14d-6/mu0
machine%cond_regions(3)%id = 6
machine%cond_regions(3)%eta = 1.14d-6/mu0
machine%cond_regions(4)%id = 7
machine%cond_regions(4)%eta = 6.9d-7/mu0
machine%cond_regions(5)%id = 8
machine%cond_regions(5)%eta = 6.9d-7/mu0
ALLOCATE(machine%coil_regions(machine%ncoil_regs))
DO j=1, machine%ncoils
  machine%coil_regions(j)%id = 8 + j
END DO
CALL compute_bcmat(machine)

!---------------------------------------------------------------------------
! Now, need to equilibrium object
!---------------------------------------------------------------------------



! CALL equil%psi%restore_local(psi_total)

! CALL gs_update_bounds(equil)


! equil%itor_target=ip_target*mu0
! equil%ip_ratio_target=ip_ratio_target
! equil%p_scale = 0.28658184156588085
! equil%ffp_scale = 2.596639717247778
! tmp_str = 'tokamaker_f.prof'
! CALL gs_profile_load(tmp_str,equil%I)
! tmp_str = 'tokamaker_p.prof'
! CALL gs_profile_load(tmp_str,equil%P)
! write(*,*) 105
! equil%I%plasma_bounds=equil%plasma_bounds
! equil%P%plasma_bounds=equil%plasma_bounds

! write(*,*) 113
! equil%vcontrol_val = 0.d0

! ALLOCATE(equil%coil_currs(machine%ncoil_regs))
! equil%coil_currs = [-10004540.249054534, 8131096.462417697, 8130193.40495448, -2752265.1799410507, -2755148.7923723585, -1429157.477860515,  -1424955.0239338675]*mu0
! ALLOCATE(areas(machine%ncoil_regs))
! areas = [1.8,0.25, 0.25, 0.25, 0.25, 0.25, 0.25 ]
! equil%coil_currs = equil%coil_currs/areas
! equil%mode = 0
! equil%I%f_offset = 36.d0

! write(*,*) 141
! dt = 0.001
! lin_tol = 1.d-11
! nl_tol = 1.d-9
! ALLOCATE(dens_reg(machine%mesh%nreg))
! dens_reg = -1.d0
! dens_reg(5) = 9806.d0
! dens_reg(6) = 9806.d0
! ALLOCATE(visc_reg(machine%mesh%nreg))
! visc_reg = -1.d0
! visc_reg(5) = 1.d-3
! visc_reg(6) = 1.d-3

! equil%device => machine
! write(*,*) 153
! CALL b_sim%setup(mg_mesh, equil, dt, lin_tol, nl_tol, dens_reg, visc_reg)
! write(*,*) 155
! CALL b_sim%step(t, dt, nl_its, l_its, nretry)




END PROGRAM gs_driver_full