MODULE oft_blanket_td
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, &
hdf5_field_exist, oft_bin_file, xdmf_plot_file, hdf5_create_file, hdf5_create_group
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_bmesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
USE oft_gauss_quadrature, ONLY: set_quad_1d
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat, oft_vector_ptr, &
  vector_extrapolate, oft_graph, oft_graph_ptr
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lu, ONLY: oft_lusolver
USE oft_la_utils, ONLY: create_matrix, graph_add_dense_blocks, create_identity_graph
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_cast
!
USE fem_base, ONLY: oft_ml_fem_type, fem_common_linkage
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec, bfem_map_flag,  bfem_interp
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE oft_blag_operators, ONLY: oft_blag_vproject,oft_blag_project, oft_blag_getmop, oft_lag_bginterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE oft_gs, ONLY: gs_epsilon, build_dels, build_dels_mug, gs_eq, gs_update_bounds, gs_test_bounds, set_bcmat
USE oft_gs_td, ONLY: oft_tmaker_td_mfop, tMaker_td_mfnk_update
USE oft_mesh_local_util, ONLY: mesh_local_findedge
IMPLICIT NONE
#include "local.h"
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE

TYPE, public :: oft_blanket_td_sim
CLASS(oft_vector), POINTER :: u => NULL() !< current solution vector
CLASS(oft_vector), POINTER :: rhs => NULL() !< Temporary RHS vector
CLASS(oft_vector), POINTER :: tmp => NULL() !< Temporary RHS vector
TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix free operator
TYPE(oft_blanket_td_mfop), POINTER :: nlfun => NULL() ! !< Time-advance operator 
TYPE(oft_lusolver), POINTER :: pre => NULL() !< Preconditioner using jacobian operator
TYPE(oft_nksolver) :: nksolver !< Newton-Krylov solver for time-advance
CLASS(oft_matrix), POINTER :: jacobian => NULL() !< Needs docs
INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
REAL(r8) :: dt = -1.d0 !< Needs docs
REAL(r8) :: t = 0.d0 !< Needs docs
LOGICAL :: pm = .FALSE.

contains
    !> Apply the matrix
    procedure :: setup => setup_blanket_td
    !> Needs Docs
    procedure :: delete => delete_blanket_td
    !> Needs Docs
    procedure :: step => step_blanket_td
END TYPE oft_blanket_td_sim

type, extends(oft_noop_matrix) ::oft_blanket_td_mfop
    real(r8) :: dt = 1.E-3_r8 !< Time step size [s]
    CLASS(oft_matrix), POINTER :: jac_op => NULL() !< Time-advance operator
    TYPE(oft_tmaker_td_mfop), POINTER :: tkmr => NULL() !< TokaMaker time-dependent object
    TYPE(oft_xmhd_2d_sim), POINTER :: mug => NULL() !< MUG time-dependent object
contains
    !> Apply operator
    procedure :: apply_real => nlfun_apply
    !> Needs Docs
    procedure :: delete => delete_mfop
end type oft_blanket_td_mfop

TYPE(oft_blanket_td_sim), POINTER :: current_sim => NULL()
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL()
CLASS(oft_scalar_bfem), POINTER :: lag_rep => NULL()
CLASS(oft_scalar_bfem), POINTER :: lag_rep_p => NULL()

CONTAINS

subroutine setup_blanket_td(self, equil, dt,lin_tol,nl_tol, dens_reg, visc_reg, eta_t_reg, eta_p_reg)
CLASS(oft_blanket_td_sim), INTENT(inout), TARGET :: self
TYPE(gs_eq), TARGET, INTENT(inout) :: equil
REAL(8), INTENT(in) :: dt !< Needs Docs
REAL(8), INTENT(in) :: lin_tol !< Needs Docs
REAL(8), INTENT(in) :: nl_tol !< Needs Docs
REAL, INTENT(in) :: dens_reg(:), visc_reg(:)
REAL, INTENT(in), optional :: eta_t_reg(:), eta_p_reg(:)
TYPE(oft_tmaker_td_mfop), POINTER :: tok_sim
TYPE(oft_xmhd_2d_sim), POINTER :: mhd_sim
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr
CLASS(oft_matrix), pointer, intent(inout) :: vac_op

mesh => equil%mesh
lag_rep=>equil%fe_rep

!------------------------------------------------------------------------------
! Set up TokaMaker and MUG objects
!------------------------------------------------------------------------------
ALLOCATE(self%nlfun)
self%nlfun%dt=dt

CALL self%nlfun%tkmr%setup_gs_td(equil, dt, lin_tol, nl_tol, .FALSE.)

ALLOCATE(mhd_sim%eta(mesh%nreg, 2))
ALLOCATE(mhd_sim%rho(mesh%nreg))
ALLOCATE(mhd_sim%nu(mesh%nreg))
mhd_sim%rho = dens_reg
mhd_sim%nu = visc_reg

IF(PRESENT(eta_t_reg)) THEN
    mhd_sim%eta(:,1) = = eta_t_reg
ELSE
    mhd_sim%eta(:,1) = self%tkmr%eta_reg
END IF
IF(PRESENT(eta_p_reg)) THEN
    mhd_sim%eta(:,2) = = eta_p_reg
ELSE
    mhd_sim%eta(:,2) = self%tkmr%eta_reg
END IF

ALLOCATE(mhd_sim%n_bc(lag_rep%ne))
ALLOCATE(mhd_sim%velx_bc(lag_rep%ne))
ALLOCATE(mhd_sim%vely_bc(lag_rep%ne))
ALLOCATE(mhd_sim%velz_bc(lag_rep%ne))
ALLOCATE(mhd_sim%by_bc(lag_rep%ne))
ALLOCATE(mhd_sim%psi_bc(lag_rep%ne))
mhd_sim%n_bc = (mhd_sim%rho < 0.d0)
mhd_sim%velx_bc = (mhd_sim%rho < 0.d0)
mhd_sim%vely_bc = (mhd_sim%rho < 0.d0)
mhd_sim%velz_bc = (mhd_sim%rho < 0.d0)
mhd_sim%by_bc = (mhd_sim%rho < 0.d0)
mhd_sim%psi_bc = (mhd_sim%rho < 0.d0)

mhd_sim%cyl = .TRUE.
mhd_sim%incomp = .TRUE.

mhd_sim%dt = dt

CALL mhd_sim%setup_from_tok(self%tkmr)
self%nlfun%mug => mhd_sim

lag_rep_p => mhd_sim%fe_rep%fields(1)%fe

!------------------------------------------------------------------------------
! Create Solver fields
!------------------------------------------------------------------------------
ALLOCATE(self%fe_rep) !CHECKBACK
self%fe_rep%nfields=6
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields)) 
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>oft_blagrange_1
self%fe_rep%field_tags(1)='p'
self%fe_rep%fields(2)%fe=>oft_blagrange_2
self%fe_rep%field_tags(2)='velx'
self%fe_rep%fields(3)%fe=>oft_blagrange_2
self%fe_rep%field_tags(3)='vely'
self%fe_rep%fields(4)%fe=>oft_blagrange_2
self%fe_rep%field_tags(4)='velz'
self%fe_rep%fields(5)%fe=>oft_blagrange_2
self%fe_rep%field_tags(5)='by'
self%fe_rep%fields(6)%fe=>oft_blagrange_2
self%fe_rep%field_tags(6)='psi'
CALL self%fe_rep%vec_create(self%u)
call self%fe_rep%vec_create(self%rhs)
call self%fe_rep%vec_create(self%tmp)

!------------------------------------------------------------------------------
! Set initial field values
!------------------------------------------------------------------------------
CALL self%u%set(1000.d0, 1)
CALL self%u%set(0.d0, 2)
CALL self%u%set(0.d0, 3)
CALL self%u%set(0.d0, 4)
CALL self%u%set(equil%I%f_offset, 5)
NULLIFY(tmp_arr)
CALL equil%psi%get_local(tmp_arr)
CALL self%u%restore_local(tmp_arr,6)

!------------------------------------------------------------------------------
! Build Jacobian matrix -> DEAL WITH THIS LATER
!------------------------------------------------------------------------------
! ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
! self%jacobian_block_mask=1
! self%jacobian_block_mask(6,6) = 0 ! Do not populate psi array, fill in later
! CALL fem_mat_create_mod(self%fe_rep, self%nlfun%jac_op, self%jacobian_block_mask)
! CALL self%tkmr%build_vac_op(vac_op)



! Preconditioner should use approximate jacobian
ALLOCATE(self%pre) !CHECKBACK
self%pre%A=>self%nlfun%jac_op 

!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat) 
self%mfmat%f=>self%nlfun
CALL self%rhs%new(self%mfmat%u0)
CALL self%rhs%new(self%mfmat%f0)
CALL self%rhs%new(self%mfmat%tmp)
CALL self%rhs%new(self%mfmat%utyp)


ALLOCATE(self%mf_solver)
self%mfmat%b0=1.d-5
self%mf_solver%A=>self%mfmat
self%mf_solver%its=1000
self%mf_solver%nrits=20
self%mf_solver%atol=self%lin_tol
self%mf_solver%itplot=1
oft_env%pm = self%pm
self%mf_solver%pm=oft_env%pm
self%mf_solver%pre=>self%pre
!------------------------------------------------------------------------------
! Setup Newton Solver
!------------------------------------------------------------------------------
self%nksolver%A=>self%nlfun
self%nksolver%J_inv=>self%mf_solver
self%nksolver%its=20
self%nksolver%atol=self%nl_tol
self%nksolver%rtol=1.d-20 ! Disable relative tolerance
self%nksolver%backtrack=.FALSE.
self%nksolver%J_update=>gs_mfnk_update
self%nksolver%up_freq=1

end subroutine

subroutine apply_rhs(self, a, b)
class(oft_blanket_td_mfop), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector) :: tmp_in, tmp_out!< Result of metric function
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr1, tmp_arr2

self%mug%nlfun%dt = 0.d0
CALL self%mug%nlfun%apply_real(a,b)
NULLIFY(tmp_arr1)
NULLIFY(tmp_arr2)
CALL b%get_local(tmp_arr1, 6)
tmp_arr1 = tmp_arr1/self%dt

CALL a%get_local(tmp_arr2, 6)
CALL lag_rep%vec_create(tmp_in)
CALL lag_rep%vec_create(tmp_out)
tmp_in%set(0.d0)
tmp_out%set(0.d0)
CALL tmp_in%restore_local(tmp_arr2)
CALL self%tkmr%mfop%apply_rhs(tmp_in,tmp_out) !NEED WAY TO IGNORE MHD CELLS HERE
CALL self%tkmr%mfop%gs_eq%zerob_bc%apply(tmp_out)
CALL tmp_out%get_local(tmp_arr2)
tmp_arr1 = tmp_arr1 + tmp_arr2
CALL b%restore_local(tmp_arr1, 6)
end subroutine

subroutine nlfun_apply(self, a, b)
class(oft_blanket_td_mfop), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector) :: tmp_in, tmp_out!< Result of metric function
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr1, tmp_arr2

self%mug%nlfun%dt = self%dt
CALL self%mug%nlfun%apply_real(a,b)
NULLIFY(tmp_arr1)
NULLIFY(tmp_arr2)
CALL b%get_local(tmp_arr1, 6)
tmp_arr1 = tmp_arr1/self%dt

CALL a%get_local(tmp_arr2, 6)
CALL lag_rep%vec_create(tmp_in)
CALL lag_rep%vec_create(tmp_out)
tmp_in%set(0.d0)
tmp_out%set(0.d0)
CALL tmp_in%restore_local(tmp_arr2)
CALL self%tkmr%mfop%apply_mfop(tmp_in,tmp_out) !NEED WAY TO IGNORE MHD CELLS HERE
CALL self%tkmr%mfop%gs_eq%zerob_bc%apply(tmp_out)
CALL tmp_out%get_local(tmp_arr2)
tmp_arr1 = tmp_arr1 + tmp_arr2
CALL b%restore_local(tmp_arr1, 6)
end subroutine

subroutine step_blanket_td(self, time,dt,nl_its,lin_its,nretry)
class(oft_blanket_td_sim), target, intent(inout) :: self !< NL operator object
REAL(8), INTENT(inout) :: time,dt
INTEGER(4), INTENT(out) :: nl_its,lin_its,nretry
INTEGER(4) :: j

! Update plasma time-advance operator
CALL self%nlfun%tkmr%mfop%update()

! Update operators if the timestep has changed
IF(dt/=self%nlfun%dt)THEN
    dt=ABS(dt)
    self%nlfun%dt=dt
    !CALL build_jac_op(self%mfop,self%mfop%vac_op) !NEED TO IMPLEMENT
    CALL self%pre%update(.TRUE.)
END IF

!Build right hand side
CALL self%tmp%add(0.d0,1.d0,self%u)
NULLIFY(tmp_arr)
CALL apply_rhs(self%nlfun,self%u,self%rhs) !FIGURE OUT WHAT TO DO WITH THIS

! Do nonlinear solve
DO j=1,4
  CALL self%nksolver%apply(self%u,self%rhs)
  IF(self%nksolver%cits<0)THEN
    CALL self%u%add(0.d0,1.d0,self%tmp)
    self%nlfun%dt=self%nlfun%dt/2.d0
    CALL build_approx_jacobian(self,self%nlfun%jac_op, self%u)
    CALL self%pre%update(.TRUE.)
    CALL apply_rhs(self%nlfun,self%u,self%rhs)
    CYCLE
  ELSE
    EXIT
  END IF
END DO

time=time+self%nlfun%dt
dt=self%nlfun%dt
nl_its=self%nksolver%nlits
lin_its=self%nksolver%lits
nretry=j-1
IF(j>4)THEN
    nretry=-nretry
ELSE
    self%nlfun%tkmr%mfop%gs_eq%alam=self%mfop%f_scale
    self%nlfun%tmkr%mfop%gs_eq%pnorm=self%mfop%p_scale
END IF
end subroutine

subroutine delete_blanket_td()
class(oft_blanket_td_sim), intent(inout) :: self !< NL operator object
INTEGER(4) :: i
DEBUG_STACK_PUSH

IF(ASSOCIATED(self%nlfun))THEN
    CALL self%nlfun%delete()
    DEALLOCATE(self%nlfun)
END IF
!
IF(ASSOCIATED(self%rhs))THEN
    CALL self%rhs%delete()
    CALL self%tmp%delete()
    DEALLOCATE(self%rhs,self%tmp)
    NULLIFY(self%u)

    CALL self%mfmat%delete()
    DEALLOCATE(self%mfmat)
    !
    CALL self%pre%delete()
    CALL self%mf_solver%delete()
    DEALLOCATE(self%mf_solver,self%vac_pre)
    !
    CALL self%nksolver%delete()
END IF
DEBUG_STACK_POP
end subroutine

subroutine delete_mfop()
class(oft_blanket_td_mfop), intent(inout) :: self !< NL operator object
DEBUG_STACK_PUSH
!
self%dt=-1.d0

CALL self%tkmr%delete()
CALL self%mug%delete() ! NEED TO IMPLEMENT

!
IF(ASSOCIATED(self%jac_op))THEN
    CALL self%jac_op%delete()
    DEALLOCATE(self%jac_op)
END IF
DEBUG_STACK_POP
end subroutine


END MODULE oft_blanket_td