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
USE fem_composite, ONLY: oft_fem_comp_type, fem_graph_create
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec, bfem_map_flag,  bfem_interp
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE oft_blag_operators, ONLY: oft_blag_vproject,oft_blag_project, oft_blag_getmop, oft_lag_bginterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE oft_gs, ONLY: gs_epsilon, build_dels, gs_eq, gs_update_bounds, gs_test_bounds, set_bcmat
USE oft_gs_td, ONLY: oft_tmaker_td_mfop, tMaker_td_mfnk_update, build_vac_op, apply_rhs
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE xmhd_2d, ONLY: oft_xmhd_2d_sim, build_approx_jacobian
IMPLICIT NONE
#include "local.h"
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE

TYPE, public :: oft_blanket_td_sim
REAL(r8) :: lin_tol = 1.d-13 !< absolute tolerance for linear solver
REAL(r8) :: nl_tol = 1.d-11 !< Needs docs
TYPE(oft_tmaker_td_mfop), POINTER :: tkmr => NULL() !< TokaMaker time-dependent object
TYPE(oft_xmhd_2d_sim), POINTER :: mug => NULL() !< MUG time-dependent object
CLASS(oft_vector), POINTER :: u => NULL() !< current solution vector
CLASS(oft_vector), POINTER :: rhs => NULL() !< Temporary RHS vector
CLASS(oft_vector), POINTER :: tmp => NULL() !< Temporary RHS vector
TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix free operator
TYPE(oft_blanket_td_mfop), POINTER :: nlfun => NULL() ! !< Time-advance operator 
TYPE(oft_native_gmres_solver), POINTER :: mf_solver => NULL() !< Outer linear solver
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
    TYPE(oft_blanket_td_sim), pointer :: parent_sim => NULL() !< pointer to parent simulation object for access to parameters
contains
    !> Apply operator
    procedure :: apply_real => nlfun_apply
    !> Needs Docs
    procedure :: delete => delete_mfop
end type oft_blanket_td_mfop

TYPE(oft_blanket_td_sim), POINTER :: current_sim => NULL()
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL()
CLASS(oft_scalar_bfem), POINTER :: lag_rep => NULL()

CONTAINS

subroutine setup_blanket_td(self, mg_mesh, equil, dt,lin_tol,nl_tol, dens_reg, visc_reg, eta_reg, incomp)
CLASS(oft_blanket_td_sim), INTENT(inout), TARGET :: self
CLASS(multigrid_mesh), INTENT(in) :: mg_mesh
TYPE(gs_eq), INTENT(inout), TARGET :: equil
REAL(8), INTENT(in) :: dt !< Needs Docs
REAL(8), INTENT(in) :: lin_tol !< Needs Docs
REAL(8), INTENT(in) :: nl_tol !< Needs Docs
INTEGER(i4) :: i
REAL(8), INTENT(in) :: dens_reg(:), visc_reg(:)
REAL, INTENT(in), optional :: eta_reg(:,:)
LOGICAL, INTENT(in), optional :: incomp
TYPE(oft_tmaker_td_mfop), POINTER :: tok_sim
TYPE(oft_xmhd_2d_sim), POINTER :: mhd_sim
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr
TYPE(oft_graph_ptr), ALLOCATABLE :: graphs(:,:), known_graphs(:)
INTEGER(i4) :: nkgraphs
TYPE(oft_graph), TARGET :: dense_graph
type(oft_1d_int), pointer, dimension(:) :: bc_nodes
integer(i4), allocatable :: dense_flag(:)

mesh => equil%mesh
lag_rep=>equil%fe_rep
!------------------------------------------------------------------------------
! Set up TokaMaker and MUG objects
!------------------------------------------------------------------------------
ALLOCATE(self%tkmr)
CALL self%tkmr%setup(equil)

ALLOCATE(mhd_sim)
ALLOCATE(mhd_sim%eta(mesh%nreg, 2))
ALLOCATE(mhd_sim%m_i(mesh%nreg))
ALLOCATE(mhd_sim%nu(mesh%nreg))
mhd_sim%m_i = dens_reg
mhd_sim%nu = visc_reg
IF(PRESENT(eta_reg)) THEN
    mhd_sim%eta = eta_reg
ELSE
    mhd_sim%eta(:,1) = self%tkmr%eta_reg
    mhd_sim%eta(:,2) = self%tkmr%eta_reg
END IF

ALLOCATE(mhd_sim%n_bc(lag_rep%ne))
ALLOCATE(mhd_sim%velx_bc(lag_rep%ne))
ALLOCATE(mhd_sim%vely_bc(lag_rep%ne))
ALLOCATE(mhd_sim%velz_bc(lag_rep%ne))
ALLOCATE(mhd_sim%by_bc(lag_rep%ne))
ALLOCATE(mhd_sim%psi_bc(lag_rep%ne))
mhd_sim%n_bc = (mhd_sim%m_i < 0.d0)
mhd_sim%velx_bc = (mhd_sim%m_i < 0.d0)
mhd_sim%vely_bc = (mhd_sim%m_i < 0.d0)
mhd_sim%velz_bc = (mhd_sim%m_i < 0.d0)
mhd_sim%by_bc = (mhd_sim%m_i < 0.d0)
mhd_sim%psi_bc = (mhd_sim%m_i < 0.d0)

mhd_sim%cyl_flag = .TRUE.
IF (PRESENT(incomp)) THEN
    mhd_sim%incomp = incomp
ELSE
    mhd_sim%incomp = .TRUE.
END IF

mhd_sim%dt = dt

CALL mhd_sim%setup(mg_mesh, lag_rep%order, lag_rep)
self%mug => mhd_sim


!------------------------------------------------------------------------------
! Create Solver fields
!------------------------------------------------------------------------------
ALLOCATE(self%fe_rep) !CHECKBACK
self%fe_rep%nfields=7
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields)) 
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>lag_rep
self%fe_rep%field_tags(1)='n' ! constant if incompressible
self%fe_rep%fields(2)%fe=>lag_rep
self%fe_rep%field_tags(2)='velx'
self%fe_rep%fields(3)%fe=>lag_rep
self%fe_rep%field_tags(3)='vely'
self%fe_rep%fields(4)%fe=>lag_rep
self%fe_rep%field_tags(4)='velz'
IF (mhd_sim%incomp) THEN
    self%fe_rep%fields(5)%fe=>mhd_sim%fe_rep%fields(5)%fe
    self%fe_rep%field_tags(5)='p'
ELSE
    self%fe_rep%fields(5)%fe=>lag_rep
    self%fe_rep%field_tags(5)='T'
END IF
self%fe_rep%fields(6)%fe=>lag_rep
self%fe_rep%field_tags(6)='psi'
self%fe_rep%fields(7)%fe=>lag_rep
self%fe_rep%field_tags(7)='F'
CALL self%fe_rep%vec_create(self%u)
call self%fe_rep%vec_create(self%rhs)
call self%fe_rep%vec_create(self%tmp)

!------------------------------------------------------------------------------
! Set initial field values
!------------------------------------------------------------------------------
CALL self%u%set(1.d0, 1)
CALL self%u%set(0.d0, 2)
CALL self%u%set(0.d0, 3)
CALL self%u%set(0.d0, 4)
CALL self%u%set(1000.d0, 5)
NULLIFY(tmp_arr)
CALL self%tkmr%gs_eq%psi%get_local(tmp_arr)
CALL self%u%restore_local(tmp_arr,6)
CALL self%u%set(self%tkmr%gs_eq%I%f_offset, 7)

!------------------------------------------------------------------------------
! Setup nl_fun object
!------------------------------------------------------------------------------
ALLOCATE(self%nlfun)
self%nlfun%dt=dt
self%nlfun%parent_sim => self
!------------------------------------------------------------------------------
! Build Jacobian matrix
!------------------------------------------------------------------------------
ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
self%jacobian_block_mask=1
CALL fem_graph_create(self%fe_rep, graphs, known_graphs, nkgraphs, self%jacobian_block_mask)

! Add dense regions to psi block
ALLOCATE(bc_nodes(1))
bc_nodes(1)%n = lag_rep%nbe
bc_nodes(1)%v => lag_rep%lbe
ALLOCATE(dense_flag(lag_rep%ne))
dense_flag = 0
dense_flag(bc_nodes(1)%v) = 1
!---Add dense blocks
CALL graph_add_dense_blocks(graphs(6,6)%g,dense_graph,dense_flag,bc_nodes)
NULLIFY(graphs(6,6)%g%kr,graphs(6,6)%g%lc)
graphs(6,6)%g%nnz=dense_graph%nnz
graphs(6,6)%g%kr=>dense_graph%kr
graphs(6,6)%g%lc=>dense_graph%lc
DEALLOCATE(dense_flag, bc_nodes)
CALL self%fe_rep%mat_create(self%nlfun%jac_op, self%jacobian_block_mask, graphs_in = graphs)
DO i=1,nkgraphs
    DEALLOCATE(known_graphs(i)%g)
END DO
DEALLOCATE(graphs, known_graphs)


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
!self%nksolver%J_update=>gs_mfnk_update
self%nksolver%up_freq=1
end subroutine setup_blanket_td

subroutine step_blanket_td(self,time,dt,nl_its,lin_its,nretry)
CLASS(oft_blanket_td_sim), target, intent(inout) :: self !< NL operator object
REAL(8), INTENT(inout) :: time,dt
INTEGER(4), INTENT(out) :: nl_its,lin_its,nretry
INTEGER(4) :: i,j,k,ierr

current_sim=>self

write(*,*) '275'
CALL apply_rhs_blanket(self%nlfun,self%u,self%rhs)
write(*,*) '277'
! Update time-advance operator
! CALL self%mfop%update()
! Update operators if the timestep has changed
! IF(dt/=self%mfop%dt)THEN
!     dt=ABS(dt)
!     self%mfop%dt=dt
!     CALL build_vac_op(self%mfop,self%mfop%vac_op)
!     IF(ASSOCIATED(self%adv_op))CALL build_jop(self%mfop,self%adv_op,self%psi_sol)
!     CALL self%vac_pre%update(.TRUE.)
! END IF
end subroutine step_blanket_td

subroutine apply_rhs_blanket(self, a, b)
class(oft_blanket_td_mfop), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: tmp_in, tmp_out!< Result of metric function
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr1, tmp_arr2
write(*,*) 296
self%parent_sim%mug%nlfun%dt = 0.d0
CALL self%parent_sim%mug%nlfun%apply_real(a,b)
NULLIFY(tmp_arr1)
NULLIFY(tmp_arr2)
CALL b%get_local(tmp_arr1, 6) 
tmp_arr1 = tmp_arr1/self%dt !Divide by dt so form of psi equation matches tokamaker implementation

CALL a%get_local(tmp_arr2, 6)
CALL lag_rep%vec_create(tmp_in)
CALL lag_rep%vec_create(tmp_out)
CALL tmp_in%set(0.d0)
CALL tmp_out%set(0.d0)
CALL tmp_in%restore_local(tmp_arr2)
CALL apply_rhs(self%parent_sim%tkmr, tmp_in,tmp_out) !NEED WAY TO IGNORE MHD CELLS HERE
CALL self%parent_sim%tkmr%gs_eq%zerob_bc%apply(tmp_out)
CALL tmp_out%get_local(tmp_arr2)
tmp_arr1 = tmp_arr1 + tmp_arr2
CALL b%restore_local(tmp_arr1, 6)
end subroutine apply_rhs_blanket

subroutine nlfun_apply(self, a, b)
class(oft_blanket_td_mfop), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: tmp_in, tmp_out!< Result of metric function
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr1, tmp_arr2

self%parent_sim%mug%nlfun%dt = self%dt
CALL self%parent_sim%mug%nlfun%apply_real(a,b)
NULLIFY(tmp_arr1)
NULLIFY(tmp_arr2)
CALL b%get_local(tmp_arr1, 6)
tmp_arr1 = tmp_arr1/self%dt

CALL a%get_local(tmp_arr2, 6)
CALL lag_rep%vec_create(tmp_in)
CALL lag_rep%vec_create(tmp_out)
CALL tmp_in%set(0.d0)
CALL tmp_out%set(0.d0)
CALL tmp_in%restore_local(tmp_arr2)
CALL self%parent_sim%tkmr%apply_real(tmp_in,tmp_out) !NEED WAY TO IGNORE MHD CELLS IN APPLY MFOP
CALL self%parent_sim%tkmr%gs_eq%zerob_bc%apply(tmp_out)
CALL tmp_out%get_local(tmp_arr2)
tmp_arr1 = tmp_arr1 + tmp_arr2
CALL b%restore_local(tmp_arr1, 6)
end subroutine nlfun_apply

subroutine build_blankettd_jacobian(self, mat, a)
class(oft_blanket_td_sim), intent(inout) :: self
class(oft_matrix), pointer, intent(inout) :: mat
class(oft_matrix), pointer:: vac_op
class(oft_vector), intent(inout) :: a ! Solution for computing Jacobian
CLASS(oft_native_matrix), POINTER :: V => NULL()
INTEGER(4) :: i, n

select type(vac_op => self%tkmr%vac_op)
type is (oft_native_matrix)
    V => vac_op
class default
    call oft_abort("vac_op must be an oft_native_matrix", &
                   "build_blankettd_jacobian", __FILE__)
end select

!Populate MUG and TokaMaker matrices
CALL build_approx_jacobian(self%mug, a)
CALL build_vac_op(self%tkmr,self%tkmr%vac_op)
self%jacobian => self%mug%jacobian

DO i = 1, V%nr
  n = V%kr(i+1) - V%kr(i)
  CALL self%jacobian%add_values( &
    [i], &
    V%lc(V%kr(i):V%kr(i+1)-1), &
    RESHAPE(V%M(V%kr(i):V%kr(i+1)-1), [1,n]), &
    1, n, iblock=6, jblock=6)
END DO

end subroutine build_blankettd_jacobian

! subroutine step_blanket_td(self, time,dt,nl_its,lin_its,nretry)
! class(oft_blanket_td_sim), target, intent(inout) :: self !< NL operator object
! REAL(8), INTENT(inout) :: time,dt
! INTEGER(4), INTENT(out) :: nl_its,lin_its,nretry
! INTEGER(4) :: j

! ! Update plasma time-advance operator
! CALL self%tkmr%mfop%update()

! ! Update operators if the timestep has changed
! IF(dt/=self%nlfun%dt)THEN
!     dt=ABS(dt)
!     self%nlfun%dt=dt
!     !CALL build_jac_op(self%mfop,self%mfop%vac_op) !NEED TO IMPLEMENT
!     CALL self%pre%update(.TRUE.)
! END IF

! !Build right hand side
! CALL self%tmp%add(0.d0,1.d0,self%u)
! NULLIFY(tmp_arr)
! CALL apply_rhs(self%nlfun,self%u,self%rhs) !FIGURE OUT WHAT TO DO WITH THIS

! ! Do nonlinear solve
! DO j=1,4
!   CALL self%nksolver%apply(self%u,self%rhs)
!   IF(self%nksolver%cits<0)THEN
!     CALL self%u%add(0.d0,1.d0,self%tmp)
!     self%nlfun%dt=self%nlfun%dt/2.d0
!     CALL build_approx_jacobian(self,self%nlfun%jac_op, self%u)
!     CALL self%pre%update(.TRUE.)
!     CALL apply_rhs(self%nlfun,self%u,self%rhs)
!     CYCLE
!   ELSE
!     EXIT
!   END IF
! END DO

! time=time+self%nlfun%dt
! dt=self%nlfun%dt
! nl_its=self%nksolver%nlits
! lin_its=self%nksolver%lits
! nretry=j-1
! IF(j>4)THEN
!     nretry=-nretry
! ELSE
!     self%nlfun%tkmr%mfop%gs_eq%alam=self%mfop%f_scale
!     self%nlfun%tmkr%mfop%gs_eq%pnorm=self%mfop%p_scale
! END IF
! end subroutine

subroutine delete_blanket_td(self)
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
    DEALLOCATE(self%mf_solver)
    !
    CALL self%nksolver%delete()
END IF
DEBUG_STACK_POP
end subroutine

subroutine delete_mfop(self)
class(oft_blanket_td_mfop), intent(inout) :: self !< NL operator object
DEBUG_STACK_PUSH
!
self%dt=-1.d0

!
IF(ASSOCIATED(self%jac_op))THEN
    CALL self%jac_op%delete()
    DEALLOCATE(self%jac_op)
END IF
DEBUG_STACK_POP
end subroutine

!IMPLEMENTTTTT
! !---------------------------------------------------------------------------
! !> Update matrix-free Jacobian on all levels with new solution
! !---------------------------------------------------------------------------
! subroutine mfnk_update(uin)
! class(oft_vector), target, intent(inout) :: uin !< Current field
! IF(oft_debug_print(1))write(*,*)'Updating 2D MUG MF-Jacobian'
! CALL current_sim%mf_mat%update(uin)
! END SUBROUTINE mfnk_update
! !---------------------------------------------------------------------------
! !> Update Jacobian matrices on all levels with new solution
! !---------------------------------------------------------------------------
! subroutine update_jacobian(uin)
! class(oft_vector), target, intent(inout) :: uin !< Current solution
! IF(oft_debug_print(1))write(*,*)'Updating 2D MUG approximate Jacobian'
! CALL build_approx_jacobian(current_sim,uin)
! END SUBROUTINE update_jacobian


END MODULE oft_blanket_td