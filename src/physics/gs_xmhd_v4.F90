!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file gs_xmhd.F90
!
!> Solve coupled non-linear grad-shafranov evolution and extended MHD
!---------------------------------------------------------------------------
MODULE gs_xmhd_v4
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
USE oft_la_utils, ONLY: create_matrix, graph_add_dense_blocks, create_identity_graph
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_cast
!
USE fem_base, ONLY: oft_ml_fem_type, fem_common_linkage
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
  CLASS(oft_vector), POINTER :: tmp => NULL() !< Temporary storage vector
  CLASS(oft_vector), POINTER :: tmp_vec => NULL() !< Temporary storage vector
  TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix free operator
  TYPE(gs_xmhd_nlfun), POINTER :: nlfun => NULL() ! !< Time-advance operator
  TYPE(oft_native_gmres_solver), POINTER :: mf_solver => NULL() !< Outer linear solver
  TYPE(oft_lusolver), POINTER :: vac_pre => NULL() !< Preconditioner using vacuum operator
  TYPE(oft_nksolver) :: nksolver !< Newton-Krylov solver for time-advance
  TYPE(xml_node), POINTER :: xml_root => NULL() !< XML root element
  TYPE(xml_node), POINTER :: xml_pre_def => NULL() !< XML element for preconditioner definition
  CLASS(oft_matrix), POINTER :: jacobian => NULL() !< Needs docs
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

TYPE(oft_gs_xmhd_sim), POINTER :: current_sim => NULL()
CLASS(multigrid_mesh), POINTER :: mg_mesh => NULL()
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL()
TYPE(oft_ml_fem_type), TARGET, PUBLIC :: ML_oft_blagrange
CLASS(oft_scalar_bfem), POINTER :: oft_blagrange => NULL()

CONTAINS

subroutine setup(self, mg_mesh_in)
class(oft_gs_xmhd_sim), intent(inout) :: self !< NL operator object
CLASS(multigrid_mesh), TARGET, intent(in) :: mg_mesh_in
integer(i4) :: ierr
integer(i4) :: order = 2
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr
CLASS(oft_native_matrix), POINTER :: A_native
!------------------------------------------------------------------------------
! Point FEM and mesh to equilibrium quantities
!------------------------------------------------------------------------------
mg_mesh=>mg_mesh_in
mesh=>mg_mesh%smesh
!---Setup FE representation
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building lagrange FE space'
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_oft_blagrange,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(oft_blagrange,ML_oft_blagrange%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)


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
CALL mesh%setup_io(self%xdmf_plot,order)
!------------------------------------------------------------------------------
! Create nonlinear operator and set it up
!------------------------------------------------------------------------------
ALLOCATE(self%nlfun)
self%nlfun%dt=self%dt
self%nlfun%eq => self%eq
self%nlfun%f_scale=self%eq%alam
self%nlfun%p_scale=self%eq%pnorm
ALLOCATE(self%nlfun%eta(mesh%nreg))
self%nlfun%eta=self%eta
ALLOCATE(self%nlfun%region_flag(mesh%nreg))
self%nlfun%region_flag = self%region_flag
ALLOCATE(self%nlfun%curr(mesh%nreg))
self%nlfun%curr = self%curr
!------------------------------------------------------------------------------
! Create Solver fields
!------------------------------------------------------------------------------
ALLOCATE(self%fe_rep)
self%fe_rep%nfields=1
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields))
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>oft_blagrange
self%fe_rep%field_tags(1)='psi'
! self%fe_rep%fields(2)%fe=>oft_blagrange
! self%fe_rep%field_tags(2)='n'
CALL self%fe_rep%vec_create(self%u)
call self%fe_rep%vec_create(self%rhs)
call self%fe_rep%vec_create(self%tmp)

! CALL self%tmp%set(0.d0)
NULLIFY(tmp_arr)
CALL self%u%set(0.d0)
CALL self%eq%psi%get_local(tmp_arr)
CALL self%u%restore_local(tmp_arr,1)
!CALL self%u%restore_local(tmp_arr,2)
!------------------------------------------------------------------------------
! Build Jacobian matrix
!------------------------------------------------------------------------------
ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
self%jacobian_block_mask=1
self%jacobian_block_mask(1,1) = 3
CALL fem_mat_create_mod(self%fe_rep, self%nlfun%vac_op, self%jacobian_block_mask)
self%jacobian_block_mask(1,1) = 1
CALL build_approx_jacobian(self,self%nlfun%vac_op)

ALLOCATE(self%vac_pre)
self%vac_pre%A=>self%nlfun%vac_op
! CALL self%vac_pre%A%apply(self%u, self%tmp)
! CALL self%tmp%get_local(tmp_arr,2)

! CALL self%fe_rep%mat_create(self%jacobian,self%jacobian_block_mask)

! CALL build_approx_jacobian(self%nlfun,self%nlfun%vac_op, 'free')
!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat)
CALL self%mfmat%setup(self%tmp,self%nlfun)
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
character(LEN=TDIFF_RST_LEN) :: rst_char
real(r8), pointer :: plot_vals(:)
real(r8) :: elapsed_time
integer(i4) :: i,j, io_stat, rst_tmp
type(oft_timer) :: mytimer
CLASS(oft_native_matrix), POINTER :: A_native
current_sim=>self
self%t=0.d0
CALL self%tmp%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(self%u, self%t, self%dt, 'gs_xmhd_'//rst_char//'.rst', 'U')
NULLIFY(plot_vals)
CALL self%xdmf_plot%add_timestep(self%t)
CALL self%tmp%get_local(plot_vals,1)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick()
    ! Update time-advance operator
    !CALL self%mfop%update()
    ! Update operators if the timestep has changed
    IF(self%dt/=self%nlfun%dt)THEN
        self%dt=ABS(self%dt)
        self%nlfun%dt=self%dt
        CALL build_approx_jacobian(self,self%nlfun%vac_op)
        CALL self%vac_pre%update(.TRUE.)
    END IF
    !Build right hand side
    CALL self%tmp%add(0.d0,1.d0,self%u)
    CALL apply_rhs(self%nlfun,self%u,self%rhs) !FIGURE OUT WHAT TO DO WITH THIS
    CALL self%eq%zerob_bc%apply(self%rhs)
    ! Do nonlinear solve
    DO j=1,4
        CALL self%nksolver%apply(self%u,self%rhs)
        IF(self%nksolver%cits<0)THEN
            CALL self%u%add(0.d0,1.d0,self%tmp)
            self%nlfun%dt=self%nlfun%dt/2.d0
            CALL build_approx_jacobian(self,self%nlfun%vac_op)
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
        CALL self%u%get_local(plot_vals, 1)
        CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
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
quad=>oft_blagrange%quad
NULLIFY( psi_weights, psi_res, alam_vals)
!---Get weights from solution vector
CALL a%get_local(psi_weights, 1)
CALL self%eq%psi%restore_local(psi_weights)
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
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(psi_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,3))
!$omp do schedule(static)

! ---------------------------------------------------------------------------
! PLASMA LOOP
! ---------------------------------------------------------------------------
DO i=1,mesh%nc
  IF(self%region_flag(mesh%reg(i))/=1)CYCLE
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  psi_weights_loc = psi_weights(cell_dofs)
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,self%eq%fe_rep%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    !---Reconstruct values of solution fields
    psi = 0.d0; dpsi=0.d0
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    DO jr=1,oft_blagrange%nce
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
    END DO
    ! new region flag mapping? 1= GS, 2 = vac, 3 = solid conductor, 4 = liquid conductor, 5 = coil?
    IF (gs_test_bounds(self%eq,coords) .AND. psi >self%eq%plasma_bounds(1)) THEN !check that we are in the plasma
        p_source = self%p_scale*self%eq%P%Fp(psi)*coords(1) 
        f_source = self%f_scale*0.5d0* self%eq%I%fp(psi)/ (coords(1) + gs_epsilon)
        diag=diag+[f_source,p_source]*jac_det*quad%wts(m)
        DO jr=1,oft_blagrange%nce
            res_loc(jr,1) = res_loc(jr,1) &
            - self%dt * basis_vals(jr) * p_source * jac_det*quad%wts(m)
            res_loc(jr,3) = res_loc(jr,3) &
            - self%dt * basis_vals(jr) * f_source * jac_det*quad%wts(m)
        END DO
    END IF
  END DO
    !---Add local values to full vector
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,1)
    alam_vals(cell_dofs(jr)) = alam_vals(cell_dofs(jr)) + res_loc(jr,3)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads, psi_weights_loc, cell_dofs,res_loc)
!$omp end parallel
DO i=1,oft_blagrange%nbe
    psi_res(oft_blagrange%lbe(i))=0.d0
    alam_vals(oft_blagrange%lbe(i))=0.d0
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
CALL b%restore_local(psi_res,1,add=.TRUE.)
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
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads, res_loc
REAL(r8), POINTER, DIMENSION(:) :: psi_weights, psi_res, alam_vals
quad=>self%eq%fe_rep%quad
NULLIFY( psi_weights, psi_res, alam_vals)
!---Get weights from solution vector
CALL a%get_local(psi_weights,1)
CALL self%eq%psi%restore_local(psi_weights)
CALL gs_update_bounds(self%eq,track_opoint=.TRUE.)
self%eq%I%plasma_bounds=self%eq%plasma_bounds
self%eq%P%plasma_bounds=self%eq%plasma_bounds

CALL b%set(0.d0)
CALL b%get_local(psi_res,1)
!$omp parallel private(m,jr,curved,coords,cell_dofs,basis_vals,basis_grads, &
!$omp psi_weights_loc,res_loc,jac_mat, &
!$omp jac_det,psi,dpsi, &
!$omp eta_loc, curr_loc)
!Edit for new fields
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(psi_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce, 2))
!$omp do schedule(static)
DO i=1,mesh%nc
  IF(self%region_flag(mesh%reg(i))==1)CYCLE
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  psi_weights_loc = psi_weights(cell_dofs)
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    !---Reconstruct values of solution fields
    psi = 0.d0; dpsi=0.d0
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    DO jr=1,oft_blagrange%nce
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
    END DO
    eta_loc = self%eta(mesh%reg(i))
    curr_loc = self%curr(mesh%reg(i))
    DO jr=1,oft_blagrange%nce
        IF (self%region_flag(mesh%reg(i))==3) THEN
            res_loc(jr,1) = res_loc(jr,1) &
            + basis_vals(jr)*psi*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
        END IF 
        IF (self%region_flag(mesh%reg(i))==4) THEN
            res_loc(jr,1) = res_loc(jr,1) &
            + basis_vals(jr)*self%dt*curr_loc*jac_det*quad%wts(m)
        END IF 
    END DO
    ! new region flag mapping? 1= GS, 2 = vac, 3 = solid conductor, 4 = coil, 5 = liquid conductor?
  END DO
    !---Add local values to full vector
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,1)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads, psi_weights_loc, cell_dofs,res_loc)
!$omp end parallel
DO i=1,oft_blagrange%nbe
    psi_res(oft_blagrange%lbe(i))=psi_weights(oft_blagrange%lbe(i))
END DO
CALL b%restore_local(psi_res,1,add=.TRUE.)
END SUBROUTINE apply_rhs

SUBROUTINE gs_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a
CALL current_sim%mfmat%update(a)
END SUBROUTINE gs_mfnk_update

SUBROUTINE build_approx_jacobian(self, mat)
class (oft_gs_xmhd_sim), intent(inout) :: self
class (oft_matrix), pointer, intent(inout) :: mat
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals, psi_weights_loc
REAL (r8) :: coords(3), eta_loc, jac_det, jac_mat(3,4)
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
CLASS(oft_vector), POINTER :: oft_lag_vec
INTEGER(i4), ALLOCATABLE, DIMENSION(:), TARGET :: cell_dofs
integer (i4) :: i, jr, jc, m
type(oft_quad_type), pointer :: quad
type(oft_1d_int), allocatable, dimension(:) :: iloc
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
LOGICAL :: curved
quad=>oft_blagrange%quad
CALL mat%zero
!--Setup thread locks
ALLOCATE(tlocks(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
  call omp_init_lock(tlocks(i))
END DO
!---
!$omp parallel private(m,jr,jc,curved,cell_dofs,basis_vals,basis_grads, &
!$omp  jac_loc,jac_mat,jac_det,eta_loc)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(psi_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce))
ALLOCATE(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
ALLOCATE(iloc(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
   iloc(i)%v=>cell_dofs
END DO
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do schedule(static)ordered
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
!---------------------------------------------------------------------------
! Quadrature Loop
!---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    eta_loc = self%eta(mesh%reg(i))
    !---Compute local matrix contributions
    DO jr=1,oft_blagrange%nce
      DO jc=1,oft_blagrange%nce
        ! Induction
        jac_loc(1, 1)%m(jr,jc) = jac_loc(1, 1)%m(jr,jc) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        IF (self%region_flag(self%eq%mesh%reg(i)) == 3) THEN
            jac_loc(1, 1)%m(jr,jc) = jac_loc(1, 1)%m(jr,jc) &
            + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
        END IF
        !jac_loc(2, 2)%m(jr,jc) = jac_loc(1, 1)%m(jr,jc)
        ! jac_loc(2, 2)%m(jr,jc) = jac_loc(2, 2)%m(jr, jc) &
        !   + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)*coords(1) &
        !   + self%dt*1.d0*DOT_PRODUCT(basis_grads(:, jr),basis_grads(:, jc))*jac_det*quad%wts(m)*coords(1)
      END DO
    END DO
  END DO
  !---Get local to global DOF mapping
  call oft_blagrange%ncdofs(i,cell_dofs)
!---Apply bc to local matrix
  DO jr=1,oft_blagrange%nce
    IF(oft_blagrange%be(cell_dofs(jr))) jac_loc(1,1)%m(jr,:)=0.d0
    !IF(oft_blagrange%be(cell_dofs(jr))) jac_loc(2,2)%m(jr,:)=0.d0
  END DO
  CALL self%fe_rep%mat_add_local(mat,jac_loc,iloc,tlocks)
END DO
deallocate(cell_dofs,basis_vals,basis_grads,jac_loc)
!$omp end parallel
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)

! apply free boundary BCs to psi
CALL set_bcmat_mod(self%eq,mat, 1,1)

CALL oft_blagrange%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)

END SUBROUTINE build_approx_jacobian

!------------------------------------------------------------------------------
!> creates matrix with dense blocks for mask = 3 (mask = 0 => nothing, mask = 1 => normal, mask = 2 => identity)
!------------------------------------------------------------------------------
subroutine fem_mat_create_mod(self,new,mask)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
CLASS(oft_matrix), POINTER, INTENT(out) :: new
INTEGER(i4), OPTIONAL, INTENT(in) :: mask(:,:)
INTEGER(i4) :: i,j,k,nknown_graphs
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: mat_mask,graph_ids
CLASS(oft_vector), POINTER :: tmp_vec
TYPE(oft_graph_ptr), ALLOCATABLE :: graphs(:,:),known_graphs(:)
TYPE(oft_graph), TARGET :: dense_graph
type(oft_1d_int), pointer, dimension(:) :: bc_nodes
integer(i4), allocatable :: dense_flag(:)
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Building composite FE matrix'
ALLOCATE(mat_mask(self%nfields,self%nfields))
mat_mask=1
IF(PRESENT(mask))mat_mask=mask
ALLOCATE(graphs(self%nfields,self%nfields))
!---Populate known graphs
ALLOCATE(known_graphs(self%nfields*self%nfields))
ALLOCATE(graph_ids(2,self%nfields*self%nfields))
graph_ids=0
nknown_graphs=0
DO i=1,self%nfields
  DO k=1,nknown_graphs
    IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type,self%fields(i)%fe%type/)))EXIT
  END DO
  IF(k<=nknown_graphs)CYCLE
  !---
  IF(oft_debug_print(3))WRITE(*,'(4X,A,2I4)')'Building graph',i,i
  nknown_graphs=nknown_graphs+1
  graph_ids(:,nknown_graphs)=(/self%fields(i)%fe%type,self%fields(i)%fe%type/)
  ALLOCATE(known_graphs(nknown_graphs)%g)
  known_graphs(nknown_graphs)%g%nr=self%fields(i)%fe%ne
  known_graphs(nknown_graphs)%g%nrg=self%fields(i)%fe%global%ne
  known_graphs(nknown_graphs)%g%nc=self%fields(i)%fe%ne
  known_graphs(nknown_graphs)%g%ncg=self%fields(i)%fe%global%ne
  known_graphs(nknown_graphs)%g%nnz=self%fields(i)%fe%nee
  known_graphs(nknown_graphs)%g%kr=>self%fields(i)%fe%kee
  known_graphs(nknown_graphs)%g%lc=>self%fields(i)%fe%lee
END DO
!---Set graphs
DO i=1,self%nfields
  DO j=1,self%nfields
    IF(mat_mask(i,j)==0)CYCLE
    IF(mat_mask(i,j)==2)THEN
      IF(i/=j)CALL oft_abort('Identity only valid on diagonal.', &
      'fem_mat_create',__FILE__)
      !---Setup identity graph
      CALL self%fields(i)%fe%vec_create(tmp_vec)
      CALL create_identity_graph(graphs(i,j)%g,tmp_vec)
      CALL tmp_vec%delete
      DEALLOCATE(tmp_vec)
      CYCLE
    END IF
    DO k=1,nknown_graphs
      IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type, &
      self%fields(j)%fe%type/)))EXIT
    END DO
    IF(k<=nknown_graphs)THEN
      IF(oft_debug_print(3))WRITE(*,'(4X,A,2I4)')'Using known graph ',i,j
      graphs(i,j)%g=>known_graphs(k)%g
    ELSE
      IF(oft_debug_print(3))WRITE(*,'(4X,A,2I4)')'Building graph ',i,j
      nknown_graphs=nknown_graphs+1
      graph_ids(:,nknown_graphs)=(/self%fields(i)%fe%type, &
      self%fields(j)%fe%type/)
      ALLOCATE(known_graphs(nknown_graphs)%g)
      known_graphs(nknown_graphs)%g%nr=self%fields(i)%fe%ne
      known_graphs(nknown_graphs)%g%nrg=self%fields(i)%fe%global%ne
      known_graphs(nknown_graphs)%g%nc=self%fields(j)%fe%ne
      known_graphs(nknown_graphs)%g%ncg=self%fields(j)%fe%global%ne
      CALL fem_common_linkage(self%fields(i)%fe,self%fields(j)%fe, &
        known_graphs(nknown_graphs)%g%nnz,known_graphs(nknown_graphs)%g%kr, &
        known_graphs(nknown_graphs)%g%lc)
      graphs(i,j)%g=>known_graphs(nknown_graphs)%g
    END IF
    IF(mat_mask(i,j)==3) THEN
      ALLOCATE(bc_nodes(1))
      bc_nodes(1)%n = self%fields(i)%fe%nbe
      bc_nodes(1)%v => self%fields(i)%fe%lbe

      ALLOCATE(dense_flag(self%fields(i)%fe%ne))
      dense_flag = 0
      dense_flag(bc_nodes(1)%v) = 1
      !---Add dense blocks
      CALL graph_add_dense_blocks(graphs(i,j)%g,dense_graph,dense_flag,bc_nodes)
      NULLIFY(graphs(i,j)%g%kr,graphs(i,j)%g%lc)
      graphs(i,j)%g%nnz=dense_graph%nnz
      graphs(i,j)%g%kr=>dense_graph%kr
      graphs(i,j)%g%lc=>dense_graph%lc
      DEALLOCATE(dense_flag)
    END IF
  END DO
END DO
!---
CALL self%vec_create(tmp_vec)
CALL create_matrix(new,graphs,tmp_vec,tmp_vec)
CALL tmp_vec%delete
DO i=1,nknown_graphs
  DEALLOCATE(known_graphs(i)%g)
END DO
DEALLOCATE(graphs,known_graphs,mat_mask,graph_ids,tmp_vec)
DEBUG_STACK_POP
end subroutine fem_mat_create_mod

!------------------------------------------------------------------------------
!> Add boundary condition terms for free-boundary case to matrix
!------------------------------------------------------------------------------
subroutine set_bcmat_mod(self,mat, iblock, jblock)
class(gs_eq), intent(inout) :: self !< G-S object
class(oft_matrix), intent(inout) :: mat !< Matrix object
integer(4), intent(in) :: iblock, jblock
integer(4) :: i,j,i_inds(1),j_inds(1)
real(8) :: one_val(1,1)
!---Add to matrix
! | A_ii A_ib |
! | M_bi M_bb + M*L^-1*M |
one_val=1.d0
DO i=1,self%fe_rep%nbe
  i_inds=self%fe_rep%lbe(i)
  IF(self%fe_flag(self%fe_rep%lbe(i)))THEN
    CALL mat%add_values(i_inds,i_inds,one_val,1,1, iblock, jblock)
  ELSE
    DO j=1,self%bc_nrhs
      j_inds=self%bc_rhs_list(j)
      IF(ABS(self%bc_bmat(i,j))<1.d-20)CYCLE
      CALL mat%add_values(i_inds,j_inds,self%bc_bmat(i:i,j:j),1,1, iblock, jblock)
    END DO
    DO j=1,self%fe_rep%nbe
      j_inds=self%fe_rep%lbe(j)
      CALL mat%add_values(i_inds,j_inds,self%bc_lmat(i:i,j:j),1,1, iblock, jblock)
    END DO
  END IF
END DO
end subroutine set_bcmat_mod

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

END MODULE gs_xmhd_v4