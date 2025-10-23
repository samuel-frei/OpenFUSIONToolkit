!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file gs_xmhd.F90
!
!> Solve coupled non-linear grad-shafranov evolution and extended MHD
!---------------------------------------------------------------------------
MODULE gs_xmhd
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, &
  hdf5_field_exist, oft_bin_file, xdmf_plot_file
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_bmesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat, oft_vector_ptr, &
  vector_extrapolate
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec, bfem_map_flag
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE oft_blag_operators, ONLY: oft_blag_vproject,oft_blag_project, oft_blag_getmop, oft_lag_bginterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE oft_gs, ONLY: gs_epsilon, build_dels, build_dels_mug, gs_eq, gs_update_bounds, gs_test_bounds
IMPLICIT NONE
#include "local.h"
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE

!------------------------------------------------------------------------------
!> Nonlinear function type for psi evolution (add later variables later)
!------------------------------------------------------------------------------
TYPE, extends(oft_noop_matrix) :: gs_xmhd_nlfun
  REAL(r8) :: dt = -1.d0 !< Time step
  REAL(r8) :: f_scale = 1.d0 !< Scale factor for \f$ F*F' \f$ term
  REAL(r8) :: p_scale = 1.d0 !< Scale factor for \f$ P' \f$ term
  REAL(r8) :: diag_vals(1) = 0.d0 !< Needs docs
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta
  
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< psi BC flag

  TYPE(gs_eq), POINTER :: eq => NULL() !< Equilibrium object
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: region_flag => NULL()
  CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => nlfun_apply
END TYPE gs_xmhd_nlfun

!------------------------------------------------------------------------------
!> Simulation object
!------------------------------------------------------------------------------
TYPE, public :: oft_gs_xmhd_sim
  LOGICAL :: mfnk = .TRUE. !< Use matrix free method?
  INTEGER(i4) :: nsteps = -1 !< Needs docs
  INTEGER(i4) :: rst_base = 0 !< Needs docs
  INTEGER(i4) :: rst_freq = 10 !< Needs docs
  REAL(r8) :: dt = -1.d0 !< Needs docs
  REAL(r8) :: t = 0.d0 !< Needs docs
  REAL(r8) :: lin_tol = 1.d-13 !< absolute tolerance for linear solver
  REAL(r8) :: nl_tol = 1.d-11 !< Needs docs
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: region_flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< psi BC flag
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
  TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  TYPE(xdmf_plot_file) :: xdmf_plot
  TYPE(gs_eq), POINTER :: eq => NULL() !< Equilibrium object
  CLASS(oft_vector), POINTER :: u => NULL() !< Needs docs
  CLASS(oft_matrix), POINTER :: jacobian => NULL() !< Needs docs
  TYPE(oft_mf_matrix), POINTER :: mf_mat => NULL() !< Matrix free operator
  TYPE(gs_xmhd_nlfun), POINTER :: nlfun => NULL() !< Needs docs
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
!
CLASS(multigrid_mesh), POINTER :: mg_mesh => NULL()
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL()
TYPE(oft_ml_fem_type), TARGET, PUBLIC :: ML_oft_blagrange
CLASS(oft_scalar_bfem), POINTER :: oft_blagrange => NULL()
CONTAINS
!---------------------------------------------------------------------------
!> RUN ACTUAL SIMULATION
!---------------------------------------------------------------------------
subroutine run_simulation(self)
class(oft_gs_xmhd_sim), target, intent(inout) :: self
type(oft_nksolver) :: nksolver
!---Solver objects
class(oft_vector), pointer :: u,v,up, g
type(oft_native_gmres_solver), target :: solver
type(oft_timer) :: mytimer
!---History file fields
TYPE(oft_bin_file) :: hist_file
integer(i4) :: hist_i4(3)
real(4) :: hist_r4(3)
!---
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
class(oft_vector), pointer :: ux,uy,uz,v_lag
class(oft_matrix), pointer :: lap_mat => NULL()!< Matrix object
type(oft_lag_bginterp) :: grad_psi
TYPE(poss_scalar_bfield) :: field_init
!---Extrapolation fields
integer(i4), parameter :: maxextrap=2
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
!---
character(LEN=TDIFF_RST_LEN) :: rst_char
integer(i4) :: i,j,io_stat,rst_tmp,npre
real(r8) :: psi_avg, elapsed_time
real(r8), pointer :: plot_vals(:),plot_vec(:,:), vals_tmp(:), tmp(:)
current_sim=>self
!---------------------------------------------------------------------------
! Create solver fields
!------------------------------------------------------------------------x---
call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)
call self%fe_rep%vec_create(g)

!---
call oft_blagrange%vec_create(grad_psi%u)
call oft_blagrange%vec_create(ux)
call oft_blagrange%vec_create(uy)
call oft_blagrange%vec_create(uz)
call oft_blagrange%vec_create(v_lag)
NULLIFY(lmop)
call oft_blag_getmop(oft_blagrange,lmop,'none')
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre)
! ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
! DO i=1,maxextrap
!   CALL self%fe_rep%vec_create(extrap_fields(i)%f)
!   extrapt(i)=0.d0
! END DO
nextrap=0
self%t=0.d0
CALL u%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(u, self%t, self%dt, 'gs_xmhd_'//rst_char//'.rst', 'U')
NULLIFY(plot_vals)
ALLOCATE(plot_vec(3,v_lag%n))
CALL self%xdmf_plot%add_timestep(self%t)
CALL self%u%get_local(plot_vals,1)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
!------------------------------------------------------------------------------
! Compute current and plot
!------------------------------------------------------------------------------
CALL build_dels_mug(lap_mat,oft_blagrange)
CALL v_lag%restore_local(plot_vals)
CALL lap_mat%apply(v_lag,ux)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,ux)
CALL v_lag%get_local(plot_vals)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'J')
!------------------------------------------------------------------------------
! Project magnetic field and plot
!------------------------------------------------------------------------------
CALL self%u%get_local(plot_vals,1)
CALL grad_psi%u%restore_local(plot_vals)
CALL grad_psi%setup(oft_blagrange)
CALL oft_blag_vproject(oft_blagrange,grad_psi,ux,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,ux)
CALL ux%add(0.d0,1.d0,v_lag)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uy)
CALL uy%add(0.d0,1.d0,v_lag)
!
CALL uy%get_local(plot_vals)
plot_vec(1,:)=-plot_vals
CALL ux%get_local(plot_vals)
plot_vec(3,:)=plot_vals
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'B')
!------------------------------------------------------------------------------
! Set parameters of nlfun object
!------------------------------------------------------------------------------
ALLOCATE(self%nlfun)
ALLOCATE(self%nlfun%eta(mesh%nreg))
self%nlfun%eq => self%eq
self%nlfun%eta=self%eta
self%nlfun%psi_bc=>self%psi_bc
self%nlfun%p_scale = self%nlfun%eq%pnorm
self%nlfun%f_scale = self%nlfun%eq%alam
IF(ALLOCATED(self%region_flag)) self%nlfun%region_flag => self%region_flag
!---------------------------------------------------------------------------
! Setup linear solver
!---------------------------------------------------------------------------
CALL build_vac_jacobian(self,u)
ALLOCATE(self%mf_mat)
CALL self%mf_mat%setup(u,self%nlfun)
self%mf_mat%b0=1.d-5
solver%A=>self%mf_mat
solver%its=400
solver%atol=self%lin_tol
solver%itplot=1
solver%nrits=20
solver%pm=oft_env%pm
NULLIFY(solver%pre)
IF(ASSOCIATED(self%xml_pre_def))THEN
  CALL create_solver_xml(solver%pre,self%xml_pre_def)
ELSE
  CALL create_diag_pre(solver%pre)
END IF
solver%pre%A=>self%jacobian 
!---------------------------------------------------------------------------
! Setup nonlinear solver
!---------------------------------------------------------------------------
nksolver%A=>self%nlfun
nksolver%J_inv=>solver
nksolver%its=30
nksolver%atol=self%nl_tol
nksolver%rtol=1.d-20 ! Disable relative tolerance
nksolver%J_update=>mfnk_update
nksolver%up_freq=1
!---------------------------------------------------------------------------
! Setup history file
!---------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('oft_xmhd2d.hist', desc="History file for non-linear thermal diffusion run")
  CALL hist_file%add_field('ts',   'i4', desc="Time step index")
  CALL hist_file%add_field('lits', 'i4', desc="Linear iteration count")
  CALL hist_file%add_field('nlits','i4', desc="Non-linear iteration count")
  CALL hist_file%add_field('time', 'r4', desc="Simulation time [s]")
  CALL hist_file%add_field('ti',   'r4', desc="Average ion temperature [eV]")
  CALL hist_file%add_field('te',   'r4', desc="Average electron temperature [eV]")
  CALL hist_file%add_field('stime','r4', desc="Walltime [s]")
  CALL hist_file%write_header
  CALL hist_file%open ! Open history file
END IF

!---------------------------------------------------------------------------
! Begin time stepping
!---------------------------------------------------------------------------
npre=0
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick()
  self%nlfun%dt=self%dt
  CALL self%nlfun%apply_real(u,g)
  CALL g%get_local(plot_vals,1)
  CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'LHS')
  self%nlfun%dt=0.0
  CALL self%nlfun%apply_real(u,v)
  CALL v%get_local(plot_vals,1)
  CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'RHS')
  CALL g%add(1.d0,-1.d0,v)
  CALL g%get_local(plot_vals,1)
  CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'error')
  psi_avg = self%nlfun%diag_vals(1)
  self%nlfun%dt=self%dt
  npre = npre + 1
  IF(MOD(npre,1)==0)THEN
    CALL update_jacobian(u)
    CALL solver%pre%update(.TRUE.)
  END IF
!   DO j=maxextrap,2,-1
!     CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
!     extrapt(j)=extrapt(j-1)
!   END DO
!   IF(nextrap<maxextrap)nextrap=nextrap+1
!   CALL extrap_fields(1)%f%add(0.d0,1.d0,u)
!   extrapt(1)=self%t
!   IF(i>maxextrap)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,self%t+self%dt,u)
  CALL nksolver%apply(u,v)
  CALL self%nlfun%apply_real(u,g)
  CALL g%add(1.d0,-1.d0,v)
  CALL g%get_local(plot_vals,1)
  CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'residual')
  IF(nksolver%cits<0)CALL oft_abort("Nonlinear solve failed","run_simulation",__FILE__)
  !---------------------------------------------------------------------------
  ! Write out initial solution progress
  !---------------------------------------------------------------------------
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    hist_i4=[self%rst_base+i-1,nksolver%lits,nksolver%nlits]
    hist_r4=REAL([self%t,psi_avg,elapsed_time],4)
103 FORMAT(' Timestep',I8,ES14.6,2X,I4,2X,I4,F12.3,ES12.2)
    WRITE(*,103)self%rst_base+i,self%t,nksolver%lits,nksolver%nlits,elapsed_time,self%dt
    IF(oft_debug_print(1))WRITE(*,*)
    CALL hist_file%write(data_i4=hist_i4, data_r4=hist_r4)
  END IF
  !---------------------------------------------------------------------------
  ! Update timestep and save solution
  !---------------------------------------------------------------------------
  self%t=self%t+self%dt
  CALL self%u%add(0.d0,1.d0,u)
  self%nlfun%eq%alam=self%nlfun%f_scale
  self%nlfun%eq%pnorm=self%nlfun%p_scale
  IF(MOD(i,self%rst_freq)==0)THEN
    IF(oft_env%head_proc)CALL mytimer%tick
    !---Create restart file
    WRITE(rst_char,104)self%rst_base+i
    READ(rst_char,104,IOSTAT=io_stat)rst_tmp
    IF((io_stat/=0).OR.(rst_tmp/=self%rst_base+i))CALL oft_abort("Step count exceeds format width", "run_simulation", __FILE__)
    CALL self%rst_save(u, self%t, self%dt, 'gs_xmhd_'//rst_char//'.rst', 'U')
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
      CALL hist_file%flush
    END IF
    !---
    CALL self%xdmf_plot%add_timestep(self%t)
    CALL self%u%get_local(plot_vals,1)
    CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
    !------------------------------------------------------------------------------
    ! Compute current and plot
    !------------------------------------------------------------------------------
    CALL build_dels_mug(lap_mat,oft_blagrange)
    CALL v_lag%restore_local(plot_vals)
    CALL lap_mat%apply(v_lag,ux)
    CALL v_lag%set(0.d0)
    CALL lminv%apply(v_lag,ux)
    CALL v_lag%get_local(plot_vals)
    CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'J')
    !------------------------------------------------------------------------------
    ! Project magnetic field and plot
    !------------------------------------------------------------------------------
    CALL self%u%get_local(plot_vals,1)
    CALL grad_psi%u%restore_local(plot_vals)
    CALL grad_psi%setup(oft_blagrange)
    CALL oft_blag_vproject(oft_blagrange,grad_psi,ux,uy,uz)
    CALL v_lag%set(0.d0)
    CALL lminv%apply(v_lag,ux)
    CALL ux%add(0.d0,1.d0,v_lag)
    CALL v_lag%set(0.d0)
    CALL lminv%apply(v_lag,uy)
    CALL uy%add(0.d0,1.d0,v_lag)
    !
    CALL uy%get_local(plot_vals)
    plot_vec(1,:)=-plot_vals
    CALL ux%get_local(plot_vals)
    plot_vec(3,:)=plot_vals
    plot_vec(1,:) = plot_vec(1,:)
    plot_vec(2,:) = plot_vec(2,:)
    plot_vec(3,:) = plot_vec(3,:)
    CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'B')
  END IF
!   IF(nksolver%lits<4)THEN
!     self%dt=self%dt*2.d0
!     npre=-1
!   ELSE IF(nksolver%lits>150)THEN
!     self%dt=self%dt/2.d0
!     npre=-1
!   END IF
END DO
CALL hist_file%close()
CALL nksolver%delete()
CALL solver%delete()
CALL u%delete()
CALL up%delete()
CALL v%delete()
DEALLOCATE(u,up,v,plot_vals)
end subroutine run_simulation


!---------------------------------------------------------------------------
!> Compute the NL error function, where we are solving F(x) = 0
!!
!! b = F(a)
!---------------------------------------------------------------------------
subroutine nlfun_apply(self,a,b)
class(gs_xmhd_nlfun), intent(inout) :: self !< NL function object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
LOGICAL :: curved
INTEGER(i4) :: i,m,jr, k,l
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: eta_loc, p_source, f_source, diag_vals(1)
REAL(r8) ::  psi, dpsi(3), diag(2), coords(3), jac_det, jac_mat(3,4)
REAL(r8) :: eta(mesh%nreg)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals, psi_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads,res_loc
REAL(r8), POINTER, DIMENSION(:) :: psi_weights, psi_res, alam_vals
quad=>oft_blagrange%quad
NULLIFY( psi_weights, psi_res, alam_vals)
!---Get weights from solution vector
CALL a%get_local(psi_weights,1)
!---
eta = self%eta !< Needs docs
self%eq%psi=>a
CALL gs_update_bounds(self%eq,track_opoint=.TRUE.)
self%f_scale = self%eq%alam
self%p_scale = self%eq%pnorm
self%eq%I%plasma_bounds=self%eq%plasma_bounds
self%eq%P%plasma_bounds=self%eq%plasma_bounds

CALL b%set(0.d0)
CALL b%get_local(psi_res, 1)
CALL b%get_local(alam_vals, 1)
diag_vals=0.d0

!$omp parallel private(m,jr,curved,coords,cell_dofs,basis_vals,basis_grads, &
!$omp psi_weights_loc,res_loc,jac_mat, &
!$omp jac_det,psi,dpsi, &
!$omp eta_loc) reduction(+:diag_vals)
!Edit for new fields
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(psi_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,2))
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

    ! new region flag mapping? 1= GS, 2 = vac, 3 = solid conductor, 4 = liquid conductor, 5 = coil?
    IF (gs_test_bounds(self%eq,coords) .AND. psi >self%eq%plasma_bounds(1)) THEN !check that we are in the plasma
        !IF (coords(1) > 1.25) write(*,*) coords(1)
        p_source = self%p_scale*self%eq%P%Fp(psi)*coords(1) 
        f_source = self%f_scale*0.5d0* self%eq%I%fp(psi)/ (coords(1) + gs_epsilon)
        diag=diag+[f_source,p_source]*jac_det*quad%wts(m)
        DO jr=1,oft_blagrange%nce
            res_loc(jr,1) = res_loc(jr,1) &
            - self%dt * basis_vals(jr) * p_source * jac_det*quad%wts(m)
            res_loc(jr,2) = res_loc(jr,2) &
            - self%dt * basis_vals(jr) * f_source * jac_det*quad%wts(m)
        END DO
    END IF
  END DO
    !---Add local values to full vector
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,1)
    alam_vals(cell_dofs(jr)) = alam_vals(cell_dofs(jr)) + res_loc(jr,2)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads, psi_weights_loc, cell_dofs,res_loc)
!$omp end parallel
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,psi_weights,psi_res,self%psi_bc) !LOOK INTO THIS, DO I NEED TO DO BEFORE RESCALING?
CALL fem_dirichlet_vec(oft_blagrange,psi_weights,alam_vals,self%psi_bc) !LOOK INTO THIS
! RESCALE EQUATIONS --> add some conditions to this?
f_source = self%eq%Itor_target/diag(1)/(1.d0+1.d0/self%eq%Ip_ratio_target)
p_source = self%eq%Itor_target/diag(2)/(self%eq%Ip_ratio_target+1.d0)
psi_res=psi_res*p_source+alam_vals*f_source
self%eq%alam=f_source*self%f_scale
diag(1)=diag(1)*f_source
self%eq%pnorm=p_source*self%p_scale
diag(2)=diag(2)*p_source

!$omp parallel private(m,jr,curved,coords,cell_dofs,basis_vals,basis_grads, &
!$omp psi_weights_loc,res_loc,jac_mat, &
!$omp jac_det,psi,dpsi, &
!$omp eta_loc) reduction(+:diag_vals)
!Edit for new fields
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(psi_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,1))
!$omp do schedule(static)
!---------------------------------------------------------------------------
! VACUUM LOOP
!---------------------------------------------------------------------------
DO i=1,mesh%nc
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

    eta_loc = eta(mesh%reg(i))
    diag_vals = diag_vals + [psi]*jac_det*quad%wts(m)
    ! new region flag mapping? 1= GS, 2 = vac, 3 = solid conductor, 4 = coil?
    DO jr=1,oft_blagrange%nce
      IF (self%region_flag(mesh%reg(i)) == 1 .OR.self%region_flag(mesh%reg(i)) == 2 .OR. self%region_flag(mesh%reg(i)) == 3 ) THEN
        res_loc(jr,1) = res_loc(jr,1) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr), dpsi)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
      END IF

      IF (self%region_flag(mesh%reg(i)) == 3) THEN
        res_loc(jr,1) = res_loc(jr,1) &
        + basis_vals(jr)*mu0*psi*jac_det*quad%wts(m)/(eta_loc*(coords(1) + gs_epsilon))
      END IF
    END DO
  END DO
    !---Add local values to full vector
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,1)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads,psi_weights_loc,cell_dofs,res_loc)
!$omp end parallel
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,psi_weights,psi_res,self%psi_bc)
!---Put results into full vector
CALL b%restore_local(psi_res,1,add=.TRUE.)
! CALL b%get_local(psi_res,1)
self%diag_vals=oft_mpi_sum(diag_vals,1)
!---Cleanup remaining storage
DEALLOCATE( psi_res, psi_weights)
end subroutine nlfun_apply

!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_vac_jacobian(self,a)
class(oft_gs_xmhd_sim), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Solution for computing jacobian
LOGICAL :: curved
INTEGER(i4) :: i,m,jr,jc, k,l
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: k_boltz=elec_charge
REAL(r8) :: diag_vals(1), eta_loc
REAL(r8) :: psi, dpsi(3), jac_mat(3,4), jac_det, coords(3)
REAL (r8) :: eta(mesh%nreg)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals, psi_weights_loc, res_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads
REAL(r8), POINTER, DIMENSION(:) :: psi_weights
type(oft_1d_int), allocatable, dimension(:) :: iloc
class(oft_vector), pointer :: tmp
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
type(oft_quad_type), pointer :: quad
quad=>oft_blagrange%quad
CALL self%jacobian%zero
NULLIFY(psi_weights)
!---Get weights from solution vector
CALL a%get_local(psi_weights,1)
!---
eta = self%eta !< Needs docs

!--Setup thread locks
ALLOCATE(tlocks(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
  call omp_init_lock(tlocks(i))
END DO
!$omp parallel private(m,jr,jc,curved,cell_dofs,basis_vals,basis_grads, psi_weights_loc, &
!$omp  psi,jac_loc,jac_mat,jac_det, dpsi,iloc,eta_loc)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(psi_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce))
ALLOCATE(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
ALLOCATE(iloc(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
   iloc(i)%v=>cell_dofs
END DO
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do schedule(static)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
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
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0

    !---Reconstruct values of solution fields
    psi = 0.d0; dpsi=0.d0
    DO jr=1,oft_blagrange%nce
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
    END DO

    diag_vals = diag_vals + [psi]*jac_det*quad%wts(m)
    eta_loc = eta(mesh%reg(i))
    !---Compute local matrix contributions
    DO jr=1,oft_blagrange%nce
      DO jc=1,oft_blagrange%nce
        ! Induction
        IF (self%region_flag(mesh%reg(i)) == 1 .OR. self%region_flag(mesh%reg(i)) == 2 .OR. self%region_flag(mesh%reg(i)) == 3) THEN
            jac_loc(1, 1)%m(jr,jc) = jac_loc(1, 1)%m(jr,jc) &
            + self%dt*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        END IF
        IF (self%region_flag(mesh%reg(i)) == 3) THEN
            jac_loc(1, 1)%m(jr,jc) = jac_loc(1, 1)%m(jr,jc) &
            + mu0*basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
        END IF
      END DO
    END DO
  END DO
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%psi_bc(cell_dofs),1)
  CALL self%fe_rep%mat_add_local(self%jacobian,jac_loc,iloc,tlocks)
END DO
!---Cleanup thread-local storage
CALL self%fe_rep%mat_destroy_local(jac_loc)
DEALLOCATE(basis_vals,basis_grads,psi_weights_loc, cell_dofs,jac_loc,iloc)
!$omp end parallel
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
IF(oft_debug_print(2))write(*,'(4X,A)')'Setting BCs'
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%psi_bc,1)
!
call self%fe_rep%vec_create(tmp)
call self%jacobian%assemble(tmp)
call tmp%delete
DEALLOCATE(tmp,psi_weights)
end subroutine build_vac_jacobian

!---------------------------------------------------------------------------
!> Setup composite FE representation and ML environment
!---------------------------------------------------------------------------
subroutine setup(self,mg_mesh_in, order)
class(oft_gs_xmhd_sim), intent(inout) :: self
CLASS(multigrid_mesh), TARGET, intent(in) :: mg_mesh_in
integer(i4), intent(in) :: order
integer(i4) :: i,j, ierr,io_unit, cond_ind, coil_ind, type
LOGICAL, ALLOCATABLE :: vert_flag(:),edge_flag(:), boundary_flag(:)
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
mg_mesh=>mg_mesh_in
mesh=>mg_mesh%smesh
IF(ASSOCIATED(self%fe_rep))CALL oft_abort("Setup can only be called once","setup",__FILE__)
IF(ASSOCIATED(oft_blagrange))CALL oft_abort("FE space already built","setup",__FILE__)
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
!---Setup FE representation
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building lagrange FE space'
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_oft_blagrange,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(oft_blagrange,ML_oft_blagrange%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
CALL self%xdmf_plot%setup("xmhd_2d")
CALL mesh%setup_io(self%xdmf_plot,order)
!---Build composite FE definition for solution field
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Creating FE type'
ALLOCATE(self%fe_rep)
self%fe_rep%nfields=1
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields))
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>oft_blagrange
self%fe_rep%field_tags(1)='psi'
!Initialize the eta array
ALLOCATE(self%eta(mesh%nreg))
self%eta = -1.d0
!---Create solution vector
CALL self%fe_rep%vec_create(self%u)
! TODO: Boundary conditions not hard coded
! Apply BCs per region type
IF (ALLOCATED(self%region_flag)) THEN
  ALLOCATE(cell_dofs(oft_blagrange%nce))
  ALLOCATE(self%psi_bc(oft_blagrange%ne)); self%psi_bc=.FALSE.
  IF (SIZE(self%region_flag) /= mesh%nreg) THEN
    CALL oft_abort("Number of region flags does not match number of regions.","setup",__FILE__)
  END IF
  DO i=1, mesh%nc
    type = self%region_flag(mesh%reg(i))
    IF (type == 1 .OR. type == 2 .OR. type == 3 ) THEN
      CALL apply_cond_bcs(self, i, cell_dofs)
    ELSE IF (type == 4) THEN
      CALL apply_supercond_bcs(self, i, cell_dofs)
    ELSE
      CALL oft_abort("Invalid region flag.","setup",__FILE__)
    END IF
  END DO
  !Fix psi on the boundaries
  ALLOCATE(vert_flag(mesh%np),edge_flag(mesh%ne))
  ALLOCATE(boundary_flag(oft_blagrange%ne))
  vert_flag=.FALSE.; edge_flag=.FALSE.
  DO i=1,mesh%nbe
    edge_flag(mesh%lbe(i))=.TRUE.
    vert_flag(mesh%le(1,mesh%lbe(i)))=.TRUE.
    vert_flag(mesh%le(2,mesh%lbe(i)))=.TRUE.
  END DO
  CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,boundary_flag)
  WHERE (self%psi_bc .OR. boundary_flag)
    self%psi_bc = .TRUE.
  END WHERE
END IF
!---Set any BCs that are not yet set
IF(.NOT.ASSOCIATED(self%psi_bc))self%psi_bc=>oft_blagrange%global%gbe


!---Create Jacobian matrix
ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
self%jacobian_block_mask=1
CALL self%fe_rep%mat_create(self%jacobian,self%jacobian_block_mask)
end subroutine setup

!---------------------------------------------------------------------------
!> Update Jacobian matrices on all levels with new solution
!---------------------------------------------------------------------------
subroutine update_jacobian(uin)
class(oft_vector), target, intent(inout) :: uin !< Current solution
IF(oft_debug_print(1))write(*,*)'Updating 2D MUG approximate Jacobian'
CALL build_vac_jacobian(current_sim,uin)
END SUBROUTINE update_jacobian
!---------------------------------------------------------------------------
!> Update Jacobian matrices on all levels with new fields
!---------------------------------------------------------------------------
subroutine mfnk_update(uin)
class(oft_vector), target, intent(inout) :: uin !< Current field
IF(oft_debug_print(1))write(*,*)'Updating 2D MUG MF-Jacobian'
CALL current_sim%mf_mat%update(uin)
END SUBROUTINE mfnk_update



subroutine apply_cond_bcs(self,cell_ind, cell_dofs)
class(oft_gs_xmhd_sim), intent(inout) :: self
INTEGER(i4) , intent(in) :: cell_ind
INTEGER(i4), POINTER, DIMENSION(:), intent(inout) :: cell_dofs
end subroutine apply_cond_bcs

subroutine apply_supercond_bcs(self,cell_ind, cell_dofs)
class(oft_gs_xmhd_sim), intent(inout) :: self
INTEGER(i4) , intent(in) :: cell_ind
INTEGER(i4), POINTER, DIMENSION(:), intent(inout) :: cell_dofs
INTEGER(i4) :: j
call oft_blagrange%ncdofs(cell_ind,cell_dofs) ! Get global index of local DOF
DO j=1, SIZE(cell_dofs)
  self%psi_bc(cell_dofs(j)) = .TRUE. ! prevent psi evolution in superconductor
END DO
end subroutine apply_supercond_bcs

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
CALL self%fe_rep%vec_save(u,filename,path)
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

END MODULE gs_xmhd