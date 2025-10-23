!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file gs_xmhd.F90
!
!> Solve coupled non-linear grad-shafranov evolution and extended MHD
!---------------------------------------------------------------------------
MODULE gs_xmhd_v3
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, &
  hdf5_field_exist, oft_bin_file, xdmf_plot_file
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_bmesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat, oft_vector_ptr, &
  vector_extrapolate, oft_graph, oft_graph_ptr
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lu, ONLY: oft_lusolver
USE oft_la_utils, ONLY: create_matrix, graph_add_dense_blocks
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_cast
!
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec, bfem_map_flag
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE oft_blag_operators, ONLY: oft_blag_vproject,oft_blag_project, oft_blag_getmop, oft_lag_bginterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE oft_gs, ONLY: gs_epsilon, build_dels, build_dels_mug, gs_eq, gs_update_bounds, gs_test_bounds, set_bcmat
USE oft_gs_td, ONLY: oft_tmaker_td_mfop, tMaker_td_mfnk_update
IMPLICIT NONE
#include "local.h"
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE

!------------------------------------------------------------------------------
!> Simulation object
!------------------------------------------------------------------------------
TYPE, public :: oft_gs_xmhd_sim
  INTEGER(i4) :: nsteps = -1 !< Needs docs
  INTEGER(i4) :: rst_base = 0 !< Needs docs
  INTEGER(i4) :: rst_freq = 1 !< Needs docs
  REAL(r8) :: dt = -1.d0 !< Needs docs
  REAL(r8) :: t = 0.d0 !< Needs docs
  REAL(r8) :: lin_tol = 1.d-13 !< absolute tolerance for linear solver
  REAL(r8) :: nl_tol = 1.d-11 !< Needs docs
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: curr
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< psi BC flag
  LOGICAL :: pm = .FALSE.
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: region_flag
  TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  TYPE(xdmf_plot_file) :: xdmf_plot
  TYPE(gs_eq), POINTER :: eq => NULL() !< Equilibrium object
  CLASS(oft_vector), POINTER :: u => NULL() !< current solution vector
  CLASS(oft_vector), POINTER :: rhs => NULL() !< Temporary RHS vector
  CLASS(oft_vector), POINTER :: psi_tmp => NULL() !< Temporary storage vector
  CLASS(oft_vector), POINTER :: tmp_vec => NULL() !< Temporary storage vector
  TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix free operator
  TYPE(gs_xmhd_nlfun), POINTER :: nlfun => NULL() ! !< Time-advance operator
  TYPE(oft_native_gmres_solver), POINTER :: mf_solver => NULL() !< Outer linear solver
  TYPE(oft_lusolver), POINTER :: vac_pre => NULL() !< Preconditioner using vacuum operator
  TYPE(oft_nksolver) :: nksolver !< Newton-Krylov solver for time-advance
  TYPE(xml_node), POINTER :: xml_root => NULL() !< XML root element
  TYPE(xml_node), POINTER :: xml_pre_def => NULL() !< XML element for preconditioner definition
  contains
  !> Setup
  PROCEDURE :: setup => setup
  !> Run simulation
  PROCEDURE :: run_simulation => run_simulation
  !> Save restart file
  PROCEDURE :: rst_save => rst_save
  !> Load restart file
  PROCEDURE :: rst_load => rst_load
END TYPE oft_gs_xmhd_sim

TYPE(oft_gs_xmhd_sim), POINTER :: current_sim => NULL()

TYPE, extends(oft_noop_matrix) :: gs_xmhd_nlfun
  REAL(r8) :: dt = -1.d0 !< Time step
  REAL(r8) :: f_scale = 1.d0 !< Scale factor for \f$ F*F' \f$ term
  REAL(r8) :: p_scale = 1.d0 !< Scale factor for \f$ P' \f$ term
  REAL(r8) :: diag_vals(2) = 0.d0 !< Needs docs
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: curr

  TYPE(gs_eq), POINTER :: eq => NULL() !< Equilibrium object
  CLASS(oft_matrix), POINTER :: vac_op => NULL() !< Vacuum time-advance operator
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: region_flag
  CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => nlfun_apply
END TYPE gs_xmhd_nlfun
CONTAINS


subroutine setup(self)
class(oft_gs_xmhd_sim), intent(inout) :: self !< NL operator object
integer(i4) :: ierr
integer(i4) :: order = 2
CLASS(oft_native_matrix), POINTER :: A_native
!---Look for XML defintion elements
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"xmhd2d",self%xml_root,ierr)
  IF(ierr==0)THEN
    !---Look for pre node
    CALL xml_get_element(self%xml_root,"pre",self%xml_pre_def,ierr)
    IF(ierr/=0)NULLIFY(self%xml_pre_def)
  ELSE
    NULLIFY(self%xml_root)
  END IF
END IF
#endif

CALL self%xdmf_plot%setup("gs_xmhd")
CALL self%eq%mesh%setup_io(self%xdmf_plot,order)
!------------------------------------------------------------------------------
! Create nonlinear operator and set it up
!------------------------------------------------------------------------------
ALLOCATE(self%nlfun)
self%nlfun%dt=self%dt
self%nlfun%eq => self%eq
self%nlfun%f_scale=self%eq%alam
self%nlfun%p_scale=self%eq%pnorm
ALLOCATE(self%nlfun%eta(self%eq%mesh%nreg))
self%nlfun%eta=self%eta
ALLOCATE(self%nlfun%region_flag(self%eq%mesh%nreg))
self%nlfun%region_flag = self%region_flag

ALLOCATE(self%nlfun%curr(self%eq%mesh%nreg))
self%nlfun%curr = self%curr
CALL build_approx_jacobian(self%nlfun,self%nlfun%vac_op, 'free')

!------------------------------------------------------------------------------
! Create Solver fields
!------------------------------------------------------------------------------

self%u=>self%eq%psi
call self%eq%fe_rep%vec_create(self%rhs)
call self%eq%fe_rep%vec_create(self%psi_tmp)

ALLOCATE(self%vac_pre)
self%vac_pre%A=>self%nlfun%vac_op
!
!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat)
CALL self%mfmat%setup(self%psi_tmp,self%nlfun)
! CALL self%mfmat%utyp%delete()
! DEALLOCATE(self%mfmat%utyp)
ALLOCATE(self%mf_solver)
self%mfmat%b0=1.d-5
self%mf_solver%A=>self%mfmat
self%mf_solver%its=100
self%mf_solver%nrits=20
self%mf_solver%atol=self%lin_tol
self%mf_solver%itplot=1
oft_env%pm = self%pm
self%mf_solver%pm=oft_env%pm

self%mf_solver%pre=>self%vac_pre

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
end subroutine setup

subroutine run_simulation(self)
class(oft_gs_xmhd_sim), target,intent(inout) :: self !< NL operator object
CLASS(oft_native_matrix), POINTER :: A_native
character(LEN=TDIFF_RST_LEN) :: rst_char
real(r8), pointer :: plot_vals(:)
real(r8) :: elapsed_time
integer(i4) :: i,j, io_stat, rst_tmp
type(oft_timer) :: mytimer
current_sim=>self
self%t=0.d0
CALL self%psi_tmp%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(self%u, self%t, self%dt, 'gs_xmhd_'//rst_char//'.rst', 'U')
NULLIFY(plot_vals)
CALL self%xdmf_plot%add_timestep(self%t)
CALL self%psi_tmp%get_local(plot_vals)
CALL self%eq%mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick()
    ! Update time-advance operator
    !CALL self%mfop%update()
    ! Update operators if the timestep has changed
    IF(self%dt/=self%nlfun%dt)THEN
        self%dt=ABS(self%dt)
        self%nlfun%dt=self%dt
        CALL build_approx_jacobian(self%nlfun,self%nlfun%vac_op, 'free')
        CALL self%vac_pre%update(.TRUE.)
    END IF
    !Build right hand side
    CALL self%psi_tmp%add(0.d0,1.d0,self%u)
    CALL apply_rhs(self%nlfun,self%u,self%rhs) !FIGURE OUT WHAT TO DO WITH THIS
    CALL self%eq%zerob_bc%apply(self%rhs)
    ! Do nonlinear solve
    DO j=1,4
        CALL self%nksolver%apply(self%u,self%rhs)
        IF(self%nksolver%cits<0)THEN
            CALL self%u%add(0.d0,1.d0,self%psi_tmp)
            self%nlfun%dt=self%nlfun%dt/2.d0
            CALL build_approx_jacobian(self%nlfun,self%nlfun%vac_op, 'free')
            CALL self%vac_pre%update(.TRUE.)
            CALL apply_rhs(self%nlfun,self%u,self%rhs) !FIGURE OUT WHAT TO DO WITH THIS
            CALL self%eq%zerob_bc%apply(self%rhs)
            CYCLE
        ELSE
            EXIT
        END IF
    END DO
    self%t=self%t+self%nlfun%dt
    self%dt=self%nlfun%dt
    self%eq%alam=self%nlfun%f_scale
    self%eq%pnorm=self%nlfun%p_scale
    IF(MOD(i,self%rst_freq)==0)THEN
        IF(oft_env%head_proc) CALL mytimer%tick
        !---Create restart file
        WRITE(rst_char,104)self%rst_base+i
        READ(rst_char,104,IOSTAT=io_stat)rst_tmp
        IF((io_stat/=0).OR.(rst_tmp/=self%rst_base+i))CALL oft_abort("Step count exceeds format width", "run_simulation", __FILE__)
        CALL self%rst_save(self%u, self%t, self%dt, 'gs_xmhd_'//rst_char//'.rst', 'U')
        IF(oft_env%head_proc)THEN
            elapsed_time=mytimer%tock()
            WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
        END IF
        !---
        CALL self%xdmf_plot%add_timestep(self%t)
        CALL self%u%get_local(plot_vals)
        CALL self%eq%mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
    END IF 
END DO
write(*,*) self%nlfun%f_scale
end subroutine run_simulation

SUBROUTINE nlfun_apply(self, a, b)
class(gs_xmhd_nlfun), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: ptmp
type(oft_quad_type), pointer :: quad
LOGICAL :: curved
INTEGER(i4) :: i,m,jr, k,l
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: p_source, f_source, diag_vals(2)
REAL(r8) ::  psi, dpsi(3), diag(2), coords(3), jac_det, jac_mat(3,4)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals, psi_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads,res_loc
REAL(r8), POINTER, DIMENSION(:) :: psi_weights, psi_res, alam_vals
quad=>self%eq%fe_rep%quad
NULLIFY( psi_weights, psi_res, alam_vals)
!---Get weights from solution vector
CALL a%get_local(psi_weights)
self%eq%psi => a 
CALL gs_update_bounds(self%eq,track_opoint=.TRUE.)
self%eq%I%plasma_bounds=self%eq%plasma_bounds
self%eq%P%plasma_bounds=self%eq%plasma_bounds

CALL b%set(0.d0)
CALL b%get_local(psi_res, 1)
CALL b%get_local(alam_vals, 1)
diag_vals=0.d0
!$omp parallel private(m,jr,curved,coords,cell_dofs,basis_vals,basis_grads, &
!$omp psi_weights_loc,res_loc,jac_mat, &
!$omp jac_det,psi,dpsi) reduction(+:diag_vals)
!Edit for new fields
ALLOCATE(basis_vals(self%eq%fe_rep%nce),basis_grads(3,self%eq%fe_rep%nce))
ALLOCATE(psi_weights_loc(self%eq%fe_rep%nce))
ALLOCATE(cell_dofs(self%eq%fe_rep%nce),res_loc(self%eq%fe_rep%nce,2))
!$omp do schedule(static)

! ---------------------------------------------------------------------------
! PLASMA LOOP
! ---------------------------------------------------------------------------
DO i=1,self%eq%mesh%nc
  IF(self%region_flag(self%eq%mesh%reg(i))/=1)CYCLE
  curved=cell_is_curved(self%eq%mesh,i) ! Straight cell test
  call self%eq%fe_rep%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  psi_weights_loc = psi_weights(cell_dofs)
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call self%eq%mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,self%eq%fe_rep%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(self%eq%fe_rep,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(self%eq%fe_rep,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = self%eq%mesh%log2phys(i,quad%pts(:,m))
    !---Reconstruct values of solution fields
    psi = 0.d0; dpsi=0.d0
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    DO jr=1,self%eq%fe_rep%nce
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
    END DO

    ! new region flag mapping? 1= GS, 2 = vac, 3 = solid conductor, 4 = liquid conductor, 5 = coil?
    IF (gs_test_bounds(self%eq,coords) .AND. psi >self%eq%plasma_bounds(1)) THEN !check that we are in the plasma
        p_source = self%p_scale*self%eq%P%Fp(psi)*coords(1) 
        f_source = self%f_scale*0.5d0* self%eq%I%fp(psi)/ (coords(1) + gs_epsilon)
        diag=diag+[f_source,p_source]*jac_det*quad%wts(m)
        DO jr=1,self%eq%fe_rep%nce
            res_loc(jr,1) = res_loc(jr,1) &
            - self%dt * basis_vals(jr) * p_source * jac_det*quad%wts(m)
            res_loc(jr,2) = res_loc(jr,2) &
            - self%dt * basis_vals(jr) * f_source * jac_det*quad%wts(m)
        END DO
    END IF
  END DO
    !---Add local values to full vector
  DO jr=1,self%eq%fe_rep%nce
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,1)
    alam_vals(cell_dofs(jr)) = alam_vals(cell_dofs(jr)) + res_loc(jr,2)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads, psi_weights_loc, cell_dofs,res_loc)
!$omp end parallel
DO i=1,self%eq%fe_rep%nbe
    psi_res(self%eq%fe_rep%lbe(i))=0.d0
    alam_vals(self%eq%fe_rep%lbe(i))=0.d0
END DO
! RESCALE EQUATIONS --> add some conditions to this?
f_source = self%eq%Itor_target/diag(1)/(1.d0+1.d0/self%eq%Ip_ratio_target)
p_source = self%eq%Itor_target/diag(2)/(self%eq%Ip_ratio_target+1.d0)
psi_res=psi_res*p_source+alam_vals*f_source
self%f_scale=f_source*self%f_scale
self%p_scale=p_source*self%p_scale
!self%eq%alam=f_source*self%f_scale
diag(1)=diag(1)*f_source
!self%eq%pnorm=p_source*self%p_scale
diag(2)=diag(2)*p_source
CALL b%restore_local(psi_res,add=.TRUE.)
CALL b%new(ptmp)
CALL self%vac_op%apply(a,ptmp)
CALL b%add(1.d0,1.d0,ptmp)
CALL ptmp%delete

DEALLOCATE(psi_res,ptmp,alam_vals)
END SUBROUTINE nlfun_apply

SUBROUTINE apply_rhs(self,a,b)
class(gs_xmhd_nlfun), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
LOGICAL :: curved
INTEGER(i4) :: i,m,jr, k,l
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: eta_loc, curr_loc
REAL(r8) ::  psi, dpsi(3), coords(3), jac_det, jac_mat(3,4)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals, psi_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads
REAL(r8), POINTER, DIMENSION(:) :: psi_weights, psi_res, alam_vals, res_loc
quad=>self%eq%fe_rep%quad
NULLIFY( psi_weights, psi_res, alam_vals)
!---Get weights from solution vector
CALL a%get_local(psi_weights)
self%eq%psi => a 
CALL gs_update_bounds(self%eq,track_opoint=.TRUE.)
self%eq%I%plasma_bounds=self%eq%plasma_bounds
self%eq%P%plasma_bounds=self%eq%plasma_bounds

CALL b%set(0.d0)
CALL b%get_local(psi_res)
!$omp parallel private(m,jr,curved,coords,cell_dofs,basis_vals,basis_grads, &
!$omp psi_weights_loc,res_loc,jac_mat, &
!$omp jac_det,psi,dpsi, &
!$omp eta_loc, curr_loc)
!Edit for new fields
ALLOCATE(basis_vals(self%eq%fe_rep%nce),basis_grads(3,self%eq%fe_rep%nce))
ALLOCATE(psi_weights_loc(self%eq%fe_rep%nce))
ALLOCATE(cell_dofs(self%eq%fe_rep%nce),res_loc(self%eq%fe_rep%nce))
!$omp do schedule(static)
DO i=1,self%eq%mesh%nc
  IF(self%region_flag(self%eq%mesh%reg(i))==1)CYCLE
  curved=cell_is_curved(self%eq%mesh,i) ! Straight cell test
  call self%eq%fe_rep%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  psi_weights_loc = psi_weights(cell_dofs)
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call self%eq%mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,self%eq%fe_rep%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(self%eq%fe_rep,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(self%eq%fe_rep,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = self%eq%mesh%log2phys(i,quad%pts(:,m))
    !---Reconstruct values of solution fields
    psi = 0.d0; dpsi=0.d0
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    DO jr=1,self%eq%fe_rep%nce
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
    END DO
    eta_loc = self%eta(self%eq%mesh%reg(i))
    curr_loc = self%curr(self%eq%mesh%reg(i))
    DO jr=1,self%eq%fe_rep%nce
        IF (self%region_flag(self%eq%mesh%reg(i))==3) THEN
            res_loc(jr) = res_loc(jr) &
            + basis_vals(jr)*psi*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
        END IF 
        IF (self%region_flag(self%eq%mesh%reg(i))==4) THEN
            res_loc(jr) = res_loc(jr) &
            + basis_vals(jr)*self%dt*curr_loc*jac_det*quad%wts(m)
        END IF 
    END DO
    ! new region flag mapping? 1= GS, 2 = vac, 3 = solid conductor, 4 = coil, 5 = liquid conductor?
  END DO
    !---Add local values to full vector
  DO jr=1,self%eq%fe_rep%nce
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads, psi_weights_loc, cell_dofs,res_loc)
!$omp end parallel
DO i=1,self%eq%fe_rep%nbe
    psi_res(self%eq%fe_rep%lbe(i))=psi_weights(self%eq%fe_rep%lbe(i))
END DO
CALL b%restore_local(psi_res,add=.TRUE.)
END SUBROUTINE apply_rhs

SUBROUTINE gs_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a
CALL current_sim%mfmat%update(a)
END SUBROUTINE gs_mfnk_update

SUBROUTINE build_approx_jacobian(self, mat, bc)
class (gs_xmhd_nlfun), intent(inout) :: self
class (oft_matrix), pointer, intent(inout) :: mat
type(oft_1d_int), pointer, dimension(:) :: bc_nodes
TYPE(oft_graph_ptr) :: graphs(1,1)
TYPE(oft_graph), TARGET :: graph1,graph2
character(LEN=*), intent(in) :: bc !< Boundary condition
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals, psi_weights_loc
REAL (r8) :: coords(3), eta_loc, jac_det, jac_mat(3,4)
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads, jac_loc(:,:)
CLASS(oft_scalar_bfem), POINTER :: oft_blagrange => NULL()
CLASS(oft_vector), POINTER :: oft_lag_vec
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
integer(i4), allocatable :: dense_flag(:)
integer (i4) :: i, jr, jc, m
type(oft_quad_type), pointer :: quad
LOGICAL :: curved

oft_blagrange => self%eq%fe_rep
quad=>oft_blagrange%quad
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  !---
  graph1%nr=oft_blagrange%ne
  graph1%nrg=oft_blagrange%global%ne
  graph1%nc=oft_blagrange%ne
  graph1%ncg=oft_blagrange%global%ne
  graph1%nnz=oft_blagrange%nee
  graph1%kr=>oft_blagrange%kee
  graph1%lc=>oft_blagrange%lee
  !---Add dense block for boundary
  IF(TRIM(bc)=="free")THEN
    ! CALL gs_mat_create(mat)
    ALLOCATE(bc_nodes(1))
    bc_nodes(1)%n=oft_blagrange%nbe
    bc_nodes(1)%v=>oft_blagrange%lbe
    ALLOCATE(dense_flag(oft_blagrange%ne))
    dense_flag=0
    dense_flag(bc_nodes(1)%v)=1
    !---Add dense blocks
    CALL graph_add_dense_blocks(graph1,graph2,dense_flag,bc_nodes)
    NULLIFY(graph1%kr,graph1%lc)
    graph1%nnz=graph2%nnz
    graph1%kr=>graph2%kr
    graph1%lc=>graph2%lc
    DEALLOCATE(dense_flag)
  END IF
  !---Create matrix
  graphs(1,1)%g=>graph1
  CALL oft_blagrange%vec_create(oft_lag_vec)
  CALL create_matrix(mat,graphs,oft_lag_vec,oft_lag_vec)
  CALL oft_lag_vec%delete
  DEALLOCATE(oft_lag_vec)
  NULLIFY(graphs(1,1)%g)
ELSE
  CALL mat%zero
END IF
!---
!$omp parallel private(m,jr,jc,curved,cell_dofs,basis_vals,basis_grads, &
!$omp  jac_loc,jac_mat,jac_det,eta_loc)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(psi_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce))
ALLOCATE(jac_loc(oft_blagrange%nce,oft_blagrange%nce))
!$omp do schedule(static)ordered
DO i=1,self%eq%mesh%nc
  curved=cell_is_curved(self%eq%mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  !CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
  jac_loc = 0.0
!---------------------------------------------------------------------------
! Quadrature Loop
!---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call self%eq%mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = self%eq%mesh%log2phys(i,quad%pts(:,m))
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    eta_loc = self%eta(self%eq%mesh%reg(i))
    !---Compute local matrix contributions
    DO jr=1,oft_blagrange%nce
      DO jc=1,oft_blagrange%nce
        ! Induction
        jac_loc(jr,jc) = jac_loc(jr,jc) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        IF (self%region_flag(self%eq%mesh%reg(i)) == 3) THEN
            jac_loc(jr,jc) = jac_loc(jr,jc) &
            + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
        END IF
      END DO
    END DO
  END DO
  !---Get local to global DOF mapping
  call oft_blagrange%ncdofs(i,cell_dofs)
    !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_blagrange%nce
        IF(oft_blagrange%be(cell_dofs(jr)))jac_loc(jr,:)=0.d0
      END DO
    CASE("free")
      DO jr=1,oft_blagrange%nce
        IF(oft_blagrange%be(cell_dofs(jr)))jac_loc(jr,:)=0.d0
      END DO
  END SELECT
  !$omp ordered
  call mat%atomic_add_values(cell_dofs,cell_dofs,jac_loc,oft_blagrange%nce,oft_blagrange%nce)
  !$omp end ordered
END DO
deallocate(cell_dofs,basis_vals,basis_grads,jac_loc)
!$omp end parallel
!---Set diagonal entries for dirichlet rows
ALLOCATE(jac_loc(1,1),cell_dofs(1))
SELECT CASE(TRIM(bc))
  CASE("zerob")
    jac_loc(1,1)=1.d0
    DO i=1,oft_blagrange%nbe
      IF(.NOT.oft_blagrange%linkage%leo(i))CYCLE
      cell_dofs=oft_blagrange%lbe(i)
      call mat%add_values(cell_dofs,cell_dofs,jac_loc,1,1)
    END DO
  CASE("free")
    CALL set_bcmat(self%eq,mat)
END SELECT
DEALLOCATE(cell_dofs,jac_loc)

CALL oft_blagrange%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)

END SUBROUTINE build_approx_jacobian

! subroutine build_vac_op(self,mat)
! class(gs_xmhd_nlfun), intent(inout) :: self
! class(oft_matrix), pointer, intent(inout) :: mat
! CALL build_dels(mat,self%eq,'free',self%dt,self%dt)
! end subroutine build_vac_op

!---------------------------------------------------------------------------
!> Save xMHD solution state to a restart file
!---------------------------------------------------------------------------
subroutine rst_save(self,u,t,dt,filename,path)
class(oft_gs_xmhd_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to save
real(r8), intent(in) :: t !< Current solution time
real(r8), intent(in) :: dt !< Current timestep
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
DEBUG_STACK_PUSH
CALL self%eq%fe_rep%vec_save(u,filename,path)
IF(oft_env%head_proc)THEN
  CALL hdf5_write(t,filename,'t')
  CALL hdf5_write(dt,filename,'dt')
END IF
DEBUG_STACK_POP
end subroutine rst_save

!---------------------------------------------------------------------------
!> Load xMHD solution state from a restart file
!---------------------------------------------------------------------------
subroutine rst_load(self,u,filename,path,t,dt)
class(oft_gs_xmhd_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to load
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
real(r8), optional, intent(out) :: t !< Time of loaded solution
real(r8), optional, intent(out) :: dt !< Timestep at loaded time
DEBUG_STACK_PUSH
CALL self%fe_rep%vec_load(u,filename,path)
IF(PRESENT(t))CALL hdf5_read(t,filename,'t')
IF(PRESENT(dt))CALL hdf5_read(dt,filename,'dt')
DEBUG_STACK_POP
end subroutine rst_load

END MODULE gs_xmhd_v3