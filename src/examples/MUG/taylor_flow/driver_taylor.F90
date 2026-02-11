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
USE fem_utils, ONLY: bfem_map_flag

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
REAL(r8), POINTER :: vec_vals(:)
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 1
INTEGER(i4) :: rst_freq = 1
INTEGER(i4) :: ndims, nl_its, l_its, nretry
INTEGER(i4) :: npoints
integer(i4), allocatable, dimension(:) :: dim_sizes
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: dt = 0.001
REAL(r8) :: t = 0.d0
REAL (r8):: ip_ratio_target = 0.205
REAL (r8):: ip_target = 7.87E6
REAL(r8), allocatable, dimension(:) :: psi_eq, psi_pert, psi_total, eta_reg,curr_reg, areas
LOGICAL, ALLOCATABLE :: vert_flag(:),edge_flag(:), boundary_flag(:)
REAL (r8):: coords(3), psi(1), q(1)
LOGICAL :: pm=.TRUE.
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

order = 2
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_blagrange_2,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(blagrange_2,ML_blagrange_2%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)


!---------------------------------------------------------------------------
! Setup time-dependent solver
!---------------------------------------------------------------------------
gs_td%nsteps = 1
gs_td%dt = dt
gs_td%lin_tol = 1.d-11
gs_td%nl_tol = 1.d-9
gs_td%pm = pm

gs_td%nu = 1.d10
gs_td%rho = 9806.d0
gs_td%B_0 = 0.d0

ALLOCATE(eta_reg(1))
ALLOCATE(gs_td%eta_t(1))
ALLOCATE(gs_td%eta_p(1))
eta_reg = 7.0d-7/mu0
gs_td%eta_t = eta_reg
gs_td%eta_p = eta_reg

ALLOCATE(curr_reg(1))
curr_reg = -1.d0
gs_td%curr = curr_reg

ALLOCATE(gs_td%region_flag(1))
gs_td%region_flag = 1

gs_td%evolve_F = .FALSE.
CALL gs_td%setup(mg_mesh, mg_mesh_1)

gs_td%psi_bc = .TRUE.
gs_td%nlfun%psi_bc = .TRUE.
gs_td%by_bc = .TRUE.
gs_td%nlfun%by_bc = .TRUE.
! gs_td%velx_bc = .TRUE.
! gs_td%nlfun%velx_bc = .TRUE.
! gs_td%vely_bc = .TRUE.
! gs_td%nlfun%vely_bc = .TRUE.
! gs_td%velz_bc = .TRUE.
! gs_td%nlfun%velz_bc = .TRUE.

gs_td%p_bc(1) = .TRUE.
gs_td%nlfun%p_bc(1) = .TRUE.
ALLOCATE(vert_flag(mg_mesh%smesh%np),edge_flag(mg_mesh%smesh%ne))
ALLOCATE(boundary_flag(blagrange_2%ne))
vert_flag=.FALSE.; edge_flag=.FALSE.
DO i=1,mg_mesh%smesh%nbe
  edge_flag(mg_mesh%smesh%lbe(i))=.TRUE.
  vert_flag(mg_mesh%smesh%le(1,mg_mesh%smesh%lbe(i)))=.TRUE.
  vert_flag(mg_mesh%smesh%le(2,mg_mesh%smesh%lbe(i)))=.TRUE.
END DO
CALL bfem_map_flag(blagrange_2,vert_flag,edge_flag,boundary_flag)
WHERE (boundary_flag)
  gs_td%velx_bc = .TRUE.
  gs_td%vely_bc = .TRUE.
  gs_td%velz_bc = .TRUE.
END WHERE


gs_td%nlfun%velx_bc = gs_td%velx_bc
gs_td%nlfun%vely_bc = gs_td%vely_bc
gs_td%nlfun%velz_bc = gs_td%velz_bc



!---Generate mass matrix
NULLIFY(u,v,mop,vec_vals) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_blag_getmop(ML_oft_blagrange_2%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL ML_oft_blagrange_2%vec_create(u)
CALL ML_oft_blagrange_2%vec_create(v)
!Set initial values for fields
field_init%mesh=>mg_mesh%smesh
field_init%func=>utheta_init
CALL oft_blag_project(ML_oft_blagrange_2%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(1.d0)
CALL u%get_local(vec_vals)
CALL gs_td%u%restore_local(vec_vals,3)

CALL gs_td%run_simulation()

! ! psi = 0.00d0
! ! q = 1.d0
! ! gs_td%dt = .015011d0/50.d0
! ! DO j=1, 50
! !   write(*,*) j
! !   write(*,*) gs_td%eq%itor_target 
! !   CALL gs_td%add_timestep(gs_td%dt)
! !   gs_td%eq%itor_target = gs_td%eq%itor_target - 0.02*ip_target*mu0
! ! END DO
! ! equil%Ip_ratio_target = equil%Ip_ratio_target*1.03
! ! CALL gs_td%add_timestep(gs_td%dt)
! CALL gs_td%run_simulation()
! ! !---Finalize enviroment
! CALL oft_finalize
CONTAINS

SUBROUTINE utheta_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
! val = -0.125d0*pt(1)+1.125d0/pt(1)
val = -0.5d0*pt(1)+1.5d0
END SUBROUTINE utheta_init

END PROGRAM gs_driver_full