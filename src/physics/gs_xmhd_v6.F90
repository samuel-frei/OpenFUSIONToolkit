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
  ! SOLVER PARAMETERS
  INTEGER(i4) :: nsteps = -1 !< Needs docs
  INTEGER(i4) :: rst_base = 0 !< Needs docs
  INTEGER(i4) :: rst_freq = 1 !< Needs docs
  REAL(r8) :: dt = -1.d0 !< Needs docs
  REAL(r8) :: t = 0.d0 !< Needs docs
  REAL(r8) :: lin_tol = 1.d-13 !< absolute tolerance for linear solver
  REAL(r8) :: nl_tol = 1.d-11 !< Needs docs
  LOGICAL :: pm = .FALSE.
  ! PHYSICS PARAMETERS
  REAL(r8) :: chi = -1.d0 !< Needs docs
  REAL(r8) :: nu = -1.d0 !< Needs docs
  REAL(r8) :: gamma = -1.d0
  REAL(r8) :: D_diff = -1.d0
  REAL(r8) :: k_boltz=elec_charge
  REAL(r8) :: m_i=proton_mass
  REAL(r8) :: den_scale = 1.d19 !< Needs docs
  REAL (r8) :: B_0(3) = 0.d0
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: curr
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: region_flag
  TYPE(gs_eq), POINTER :: eq => NULL() !< Equilibrium object
  ! BOUNDARY CONDITIONS
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: n_bc => NULL() !< n BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velx_bc => NULL() !< velx BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: vely_bc => NULL() !< vely BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velz_bc => NULL() !< velz BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: T_bc => NULL() !< T BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: by_bc => NULL() !< By (F) BC flag
  ! SOLVER OBJECTS
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
  TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  TYPE(xdmf_plot_file) :: xdmf_plot
  CLASS(oft_vector), POINTER :: u => NULL() !< current solution vector
  CLASS(oft_vector), POINTER :: rhs => NULL() !< Temporary RHS vector
  CLASS(oft_vector), POINTER :: tmp => NULL() !< Temporary storage vector
  CLASS(oft_vector), POINTER :: tmp_vec => NULL() !< Temporary storage vector
  TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix free operator
  TYPE(gs_xmhd_nlfun), POINTER :: nlfun => NULL() ! !< Time-advance operator
  TYPE(oft_native_gmres_solver), POINTER :: mf_solver => NULL() !< Outer linear solver
  TYPE(oft_lusolver), POINTER :: pre => NULL() !< Preconditioner using jacobian operator
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
  ! SOLVER PARAMETERS
  REAL(r8) :: dt = -1.d0 !< Time step
  ! PHYSICS PARAMETERS
  REAL(r8) :: chi
  REAL(r8) :: nu = -1.d0 !< Needs docs
  REAL(r8) :: gamma = -1.d0
  REAL(r8) :: D_diff = -1.d0
  REAL(r8) :: k_boltz=elec_charge
  REAL(r8) :: m_i=proton_mass
  REAL(r8) :: den_scale = 1.d19 !< Needs docs
  REAL (r8) :: B_0(3) = 0.d0
  TYPE(gs_eq), POINTER :: eq => NULL() !< Equilibrium object
  REAL(r8) :: f_scale = 1.d0 !< Scale factor for \f$ F*F' \f$ term
  REAL(r8) :: p_scale = 1.d0 !< Scale factor for \f$ P' \f$ term
  REAL(r8) :: diag_vals(2) = 0.d0 !< Used to determine f and p scales
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: curr
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: region_flag
  ! BOUNDARY CONDITIONS
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: n_bc => NULL() !< n BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velx_bc => NULL() !< velx BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: vely_bc => NULL() !< vely BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velz_bc => NULL() !< velz BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: T_bc => NULL() !< T BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: by_bc => NULL() !< By (F) BC flag
  ! SOLVER OBJECTS
  CLASS(oft_matrix), POINTER :: vac_op => NULL() !< Vacuum time-advance operator
  CLASS(oft_matrix), POINTER :: jac_op => NULL() !< Vacuum time-advance operator
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
integer(i4) :: ierr, i, type
integer(i4) :: order = 2
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr
CLASS(oft_native_matrix), POINTER :: A_native
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
!------------------------------------------------------------------------------
! Setup mesh and finite element representation
!------------------------------------------------------------------------------
mg_mesh=>mg_mesh_in
mesh=>mg_mesh%smesh
!---Setup FE representation
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building lagrange FE space'
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_oft_blagrange,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(oft_blagrange,ML_oft_blagrange%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
!------------------------------------------------------------------------------
! Setup boundary conditions by region
! 1- extended MHD
! 2 - vacuum
! 3 - solid conductor
! 4 - coil
! 5 - plasma (grad-shafranov)
!------------------------------------------------------------------------------
! Apply BCs per region type
IF (ALLOCATED(self%region_flag)) THEN
  ALLOCATE(cell_dofs(oft_blagrange%nce))
  ALLOCATE(self%n_bc(oft_blagrange%ne)); self%n_bc=.FALSE.
  ALLOCATE(self%velx_bc(oft_blagrange%ne)); self%velx_bc=.FALSE.
  ALLOCATE(self%vely_bc(oft_blagrange%ne)); self%vely_bc=.FALSE.
  ALLOCATE(self%velz_bc(oft_blagrange%ne)); self%velz_bc=.FALSE.
  ALLOCATE(self%T_bc(oft_blagrange%ne)); self%T_bc=.FALSE.
  ALLOCATE(self%by_bc(oft_blagrange%ne)); self%by_bc=.TRUE.  ! FOR NOW WE'RE NOT EVOLVING By (F)
  IF (SIZE(self%region_flag) /= mesh%nreg) THEN
    CALL oft_abort("Number of region flags does not match number of regions.","setup",__FILE__)
  END IF
  DO i=1, mesh%nc
    type = self%region_flag(mesh%reg(i))
    IF (type == 1) THEN
      CALL apply_mhd_bcs(self, i, cell_dofs)
    ELSE IF (type >1 .AND. type < 6) THEN
      CALL apply_bcs(self, i, cell_dofs)
    ELSE
      CALL oft_abort("Invalid region flag.","setup",__FILE__)
    END IF
  END DO
END IF
!------------------------------------------------------------------------------
!---Look for XML defintion elements for preconditioner
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! Setup plotting
!------------------------------------------------------------------------------
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
self%nlfun%chi = self%chi
self%nlfun%nu = self%nu
self%nlfun%gamma = self%gamma
self%nlfun%D_diff = self%D_diff
self%nlfun%k_boltz = self%k_boltz
self%nlfun%m_i = self%m_i
self%nlfun%den_scale = self%den_scale
self%nlfun%B_0 = self%B_0

ALLOCATE(self%nlfun%eta(mesh%nreg))
self%nlfun%eta=self%eta
ALLOCATE(self%nlfun%curr(mesh%nreg))
self%nlfun%curr = self%curr
ALLOCATE(self%nlfun%region_flag(mesh%nreg))
self%nlfun%region_flag = self%region_flag

self%nlfun%n_bc=>self%n_bc
self%nlfun%velx_bc=>self%velx_bc
self%nlfun%vely_bc=>self%vely_bc
self%nlfun%velz_bc=>self%velz_bc
self%nlfun%T_bc=>self%T_bc
self%nlfun%by_bc=>self%by_bc

!------------------------------------------------------------------------------
! Create Solver fields
!------------------------------------------------------------------------------
ALLOCATE(self%fe_rep)
self%fe_rep%nfields=7
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields))
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>oft_blagrange
self%fe_rep%field_tags(1)='n'
self%fe_rep%fields(2)%fe=>oft_blagrange
self%fe_rep%field_tags(2)='velx'
self%fe_rep%fields(3)%fe=>oft_blagrange
self%fe_rep%field_tags(3)='vely'
self%fe_rep%fields(4)%fe=>oft_blagrange
self%fe_rep%field_tags(4)='velz'
self%fe_rep%fields(5)%fe=>oft_blagrange
self%fe_rep%field_tags(5)='T'
self%fe_rep%fields(6)%fe=>oft_blagrange
self%fe_rep%field_tags(6)='psi'
self%fe_rep%fields(7)%fe=>oft_blagrange
self%fe_rep%field_tags(7)='by'
CALL self%fe_rep%vec_create(self%u)
call self%fe_rep%vec_create(self%rhs)
call self%fe_rep%vec_create(self%tmp)
NULLIFY(tmp_arr)
! CALL self%u%set(0.d0)
! CALL self%eq%psi%get_local(tmp_arr)
! CALL self%u%restore_local(tmp_arr,1)
! CALL self%u%restore_local(tmp_arr,2)
! CALL self%u%restore_local(tmp_arr,3)
! CALL self%u%restore_local(tmp_arr,4)
! CALL self%u%restore_local(tmp_arr,5)
! CALL self%u%restore_local(tmp_arr,6)
! CALL self%u%restore_local(tmp_arr,7)
!------------------------------------------------------------------------------
! Build Jacobian matrix
!------------------------------------------------------------------------------
ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
self%jacobian_block_mask=1
self%jacobian_block_mask(6,6) = 3 ! Add dense regions for the boundary in the psi/psi matrix
CALL fem_mat_create_mod(self%fe_rep, self%nlfun%jac_op, self%jacobian_block_mask)
CALL fem_mat_create_mod(self%fe_rep, self%nlfun%vac_op, self%jacobian_block_mask)
self%jacobian_block_mask(6,6) = 1 !Set back to normal for creating local matrices
CALL build_vac_jacobian(self,self%nlfun%vac_op)

! Preconditioner should use approximate jacobian
ALLOCATE(self%pre)
self%pre%A=>self%nlfun%jac_op 
!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat)
!CALL self%mfmat%setup(self%rhs,self%nlfun)
self%mfmat%f=>self%nlfun
CALL self%rhs%new(self%mfmat%u0)
CALL self%rhs%new(self%mfmat%f0)
CALL self%rhs%new(self%mfmat%tmp)
CALL self%rhs%new(self%mfmat%utyp)



ALLOCATE(self%mf_solver)
self%mfmat%b0=1.d-5
self%mf_solver%A=>self%mfmat
self%mf_solver%its=400
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
end subroutine setup

subroutine run_simulation(self)
class(oft_gs_xmhd_sim), target,intent(inout) :: self !< NL operator object
character(LEN=TDIFF_RST_LEN) :: rst_char
real(r8), pointer :: plot_vals(:), tmp_arr(:), plot_vec(:,:)
real(r8) :: elapsed_time
integer(i4) :: i,j, io_stat, rst_tmp
type(oft_timer) :: mytimer
CLASS(oft_native_matrix), POINTER :: A_native
class(oft_vector), pointer :: tmp_vec
current_sim=>self
self%t=0.d0
CALL oft_blagrange%vec_create(tmp_vec)
CALL self%tmp%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(self%u, self%t, self%dt, 'gs_xmhd_'//rst_char//'.rst', 'U')
NULLIFY(plot_vals)
ALLOCATE(plot_vec(3,tmp_vec%n))
CALL self%xdmf_plot%add_timestep(self%t)
CALL self%u%get_local(plot_vals,1)
plot_vals = plot_vals*self%den_scale
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'n')
CALL self%u%get_local(plot_vals,2)
plot_vec(1,:)=plot_vals
CALL self%u%get_local(plot_vals,3)
plot_vec(3,:)=plot_vals
CALL self%u%get_local(plot_vals,4)
plot_vec(2,:)=plot_vals
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'V')
CALL self%u%get_local(plot_vals,5)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'T')
CALL self%u%get_local(plot_vals,6)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
CALL self%u%get_local(plot_vals,7)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'by')
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick()
    ! Update time-advance operator
    CALL build_approx_jacobian(self,self%nlfun%jac_op, self%u)
    CALL self%pre%update(.TRUE.)
    ! Update operators if the timestep has changed
    IF(self%dt/=self%nlfun%dt)THEN
        self%dt=ABS(self%dt)
        self%nlfun%dt=self%dt
    END IF
    !Build right hand side
    CALL self%tmp%add(0.d0,1.d0,self%u)
    NULLIFY(tmp_arr)
    CALL apply_rhs(self%nlfun,self%u,self%rhs) !FIGURE OUT WHAT TO DO WITH THIS
    CALL self%rhs%get_local(tmp_arr,6)
    CALL tmp_vec%restore_local(tmp_arr)
    CALL self%eq%zerob_bc%apply(tmp_vec)
    CALL tmp_vec%get_local(tmp_arr)
    CALL self%rhs%restore_local(tmp_arr,6)
    ! Do nonlinear solve
    DO j=1,4
        CALL self%nksolver%apply(self%u,self%rhs)
        IF(self%nksolver%cits<0)THEN
            CALL self%u%add(0.d0,1.d0,self%tmp)
            self%dt = self%dt/2.d0
            self%nlfun%dt=self%dt
            CALL build_approx_jacobian(self,self%nlfun%jac_op, self%u)
            CALL build_vac_jacobian(self,self%nlfun%vac_op)
            CALL self%pre%update(.TRUE.)
            CALL apply_rhs(self%nlfun,self%u,self%rhs) !FIGURE OUT WHAT TO DO WITH THIS
            CALL self%rhs%get_local(tmp_arr,6)
            CALL tmp_vec%restore_local(tmp_arr)
            CALL self%eq%zerob_bc%apply(tmp_vec)
            CALL tmp_vec%get_local(tmp_arr)
            CALL self%rhs%restore_local(tmp_arr,6)
            CYCLE
        ELSE
            EXIT
        END IF
    END DO
    write(*,*) 'EXIT CODE', j
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
        CALL self%u%get_local(plot_vals,1)
        plot_vals = plot_vals*self%den_scale
        CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'n')
        CALL self%u%get_local(plot_vals,2)
        plot_vec(1,:)=plot_vals
        CALL self%u%get_local(plot_vals,3)
        plot_vec(3,:)=plot_vals
        CALL self%u%get_local(plot_vals,4)
        plot_vec(2,:)=plot_vals
        CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'V')
        CALL self%u%get_local(plot_vals,5)
        CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'T')
        CALL self%u%get_local(plot_vals,6)
        CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
        CALL self%u%get_local(plot_vals,7)
        CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'by')
    END IF 
END DO
write(*,*) self%nlfun%f_scale
end subroutine run_simulation

SUBROUTINE nlfun_apply(self, a, b)
class(gs_xmhd_nlfun), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: ptmp !temporary storage vector
type(oft_quad_type), pointer :: quad
LOGICAL :: curved
INTEGER(i4) :: i,m,jr, k,l !indexing variables for loops
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: eta_loc, curr_loc
REAL(r8) :: chi, nu, D_diff, gamma, k_boltz, m_i, B_0(3) ! physics parameters
REAL(r8) :: p_source, f_source, diag(2) ! Used for scaling P' and FF'
REAL(r8) :: n, dn(3), vel(3), T, dT(3), psi, dpsi(3), by, dby(3), dvel(3,3), div_vel, btmp(3) !reconstructed variables
REAL(r8) :: coords(3), jac_det, jac_mat(3,4), tmp1(3) ! For integration
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals, n_weights_loc, T_weights_loc, psi_weights_loc, by_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads,  vel_weights_loc, res_loc
REAL(r8), POINTER, DIMENSION(:) :: n_weights, T_weights, psi_weights, by_weights
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
REAL(r8), POINTER, DIMENSION(:) :: n_res, velx_res, vely_res, velz_res, T_res, psi_res, pres_vals,alam_vals, by_res, vtmp

quad=>oft_blagrange%quad

NULLIFY( n_weights, n_res, vel_weights, velx_res, vely_res, velz_res, T_weights, T_res,&
 psi_weights, psi_res, pres_vals, alam_vals, by_weights, by_res, vtmp)

!---Get weights from solution vector
 ALLOCATE(vel_weights(3,oft_blagrange%ne))
CALL a%get_local(n_weights, 1)
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(T_weights, 5)
CALL a%get_local(psi_weights, 6)
CALL a%get_local(by_weights, 7)

!--Update equilibrium with current values of psi
CALL self%eq%psi%restore_local(psi_weights)
CALL gs_update_bounds(self%eq,track_opoint=.TRUE.)
self%eq%I%plasma_bounds=self%eq%plasma_bounds
self%eq%P%plasma_bounds=self%eq%plasma_bounds

!--Initialize residuals with zeros
CALL b%set(0.d0)
CALL b%get_local(n_res, 1)
CALL b%get_local(velx_res, 2)
CALL b%get_local(vely_res, 3)
CALL b%get_local(velz_res, 4)
CALL b%get_local(T_res, 5)
CALL b%get_local(psi_res, 6)
CALL b%get_local(alam_vals, 6)
CALL b%get_local(pres_vals, 6)
CALL b%get_local(by_res, 7)

!--Set local physics parameters
chi = self%chi
nu = self%nu
D_diff = self%D_diff
gamma = self%gamma
k_boltz = self%k_boltz
m_i = self%m_i
B_0 = self%B_0

diag = 0.d0
! Declare variables private for OMP
!$omp parallel private(m,jr,curved,coords,cell_dofs,basis_vals,basis_grads, &
!$omp n_weights_loc, vel_weights_loc, T_weights_loc, psi_weights_loc, by_weights_loc, &
!$omp res_loc,jac_mat, jac_det, &
!$omp n, dn, vel, dvel, div_vel,T, dT, psi, dpsi, by, dby) 

! Allocate local variables
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(n_weights_loc(oft_blagrange%nce))
ALLOCATE(vel_weights_loc(3, oft_blagrange%nce))
ALLOCATE(T_weights_loc(oft_blagrange%nce))
ALLOCATE(psi_weights_loc(oft_blagrange%nce))
ALLOCATE(by_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,9)) ! 7 fields + entries to store pres_vals and alam_vals

!$omp do schedule(static)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function

  ! Set local weights
  n_weights_loc = n_weights(cell_dofs)
  vel_weights_loc = vel_weights(:, cell_dofs)
  T_weights_loc = T_weights(cell_dofs)
  psi_weights_loc = psi_weights(cell_dofs)
  by_weights_loc = by_weights(cell_dofs)
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

    ! Switch basis grads from 2D to 3D
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    !---Reconstruct values of solution fields
    n = 0.d0; dn=0.d0
    vel = 0.d0; dvel = 0.d0; div_vel = 0.d0
    T = 0.d0; dT = 0.d0
    psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0
    DO jr=1,oft_blagrange%nce
      n = n + n_weights_loc(jr)*basis_vals(jr)
      vel = vel + vel_weights_loc(:, jr)*basis_vals(jr)
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      T = T + T_weights_loc(jr)*basis_vals(jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      by = by + by_weights_loc(jr)*basis_vals(jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    n = n * self%den_scale
    dn = dn * self%den_scale
    eta_loc = self%eta(mesh%reg(i))

    div_vel = dvel(1,1) +vel(1)/(coords(1)+gs_epsilon) + dvel(3,3)
    btmp = cross_product(dpsi/(coords(1)+gs_epsilon), [0.d0,1.d0,0.d0]) + [0.d0,1.d0,0.d0]*by/(coords(1)+gs_epsilon) + B_0
    IF(self%region_flag(self%eq%mesh%reg(i)) == 1) THEN
      DO jr=1,oft_blagrange%nce
        ! DENSITY
        res_loc(jr,1) = res_loc(jr, 1) &
            + basis_vals(jr)*n*jac_det*quad%wts(m)*coords(1) &
            +self%dt*basis_vals(jr)*DOT_PRODUCT(dn, vel)*jac_det*quad%wts(m)*coords(1) &
            !+self%dt*basis_vals(jr)*n*div_vel*jac_det*quad%wts(m)*coords(1) &
            +self%dt*D_diff*DOT_PRODUCT(dn,basis_grads(:,jr))*jac_det*quad%wts(m)*coords(1)
        ! VELOCITY
        res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          + basis_vals(jr)*vel*jac_det*quad%wts(m)*coords(1) &
          + self%dt*DOT_PRODUCT(btmp,basis_grads(:,jr))*btmp*jac_det*quad%wts(m)*coords(1)/(mu0*m_i*n) & 
          - self%dt*DOT_PRODUCT(btmp,btmp)*basis_grads(:,jr)*jac_det*quad%wts(m)*coords(1)/(2*mu0*m_i*n) & 
          - self%dt*basis_vals(jr)*DOT_PRODUCT(dn,btmp)*btmp*jac_det*quad%wts(m)*coords(1)/(mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(btmp,btmp)*dn*jac_det*quad%wts(m)*coords(1)/(2*mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*dT*jac_det*quad%wts(m)*coords(1)/m_i &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*T*dn*jac_det*quad%wts(m)*coords(1)/(m_i*n)
        DO k=1,3
          res_loc(jr,k+1) = res_loc(jr, k+1) &
            + basis_vals(jr)*self%dt*DOT_PRODUCT(vel,dvel(k,:))*jac_det*quad%wts(m)*coords(1) &
            + nu*self%dt*DOT_PRODUCT(basis_grads(:,jr),dvel(k,:))*jac_det*quad%wts(m)*coords(1)/(m_i*n) &
            - nu*basis_vals(jr)*self%dt*DOT_PRODUCT(dn,dvel(k,:))*jac_det*quad%wts(m)*coords(1)/(m_i*n**2)
        END DO 
        res_loc(jr,2) = res_loc(jr,2) &
          + self%dt*basis_vals(jr)*(btmp(2)**2-btmp(1)**2-btmp(3)**2)*jac_det*quad%wts(m)/(2.d0*mu0*m_i*n) &
          + self%dt*basis_vals(jr)*nu*vel(1)*jac_det*quad%wts(m)/(m_i*n*(coords(1)+gs_epsilon))  &
          - self%dt*basis_vals(jr)*vel(2)**2*jac_det*quad%wts(m)
      
        res_loc(jr,3) = res_loc(jr,3) &
          - self%dt*basis_vals(jr)*btmp(1)*btmp(2)*jac_det*quad%wts(m)/(mu0*m_i*n) & !!nonzero 
          + self%dt*basis_vals(jr)*nu*vel(2)*jac_det*quad%wts(m)/(m_i*n*(coords(1)+gs_epsilon)) &
          + self%dt*basis_vals(jr)*vel(1)*vel(2)*jac_det*quad%wts(m)

        ! TEMPERATURE
        res_loc(jr,5) = res_loc(jr, 5) &
          + basis_vals(jr)*T*jac_det*quad%wts(m)*coords(1)/(gamma-1) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dT)*jac_det*quad%wts(m)*coords(1)/(gamma-1) &
          !+ self%dt*basis_vals(jr)*T*div_vel*jac_det*quad%wts(m)*coords(1) & 
          + self%dt*chi*DOT_PRODUCT(dT, basis_grads(:,jr))*jac_det*quad%wts(m)*coords(1) &
          - self%dt*chi*basis_vals(jr)*DOT_PRODUCT(dn, dT)*jac_det*quad%wts(m)*coords(1)/n
          
        !PSI (Here, I only include the terms that are not included in vac_op)
        IF(self%region_flag(self%eq%mesh%reg(i)) == 1) THEN
          res_loc(jr,6) = res_loc(jr,6) &
            + basis_vals(jr)*DOT_PRODUCT(vel, dpsi)*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon)) &
            + basis_vals(jr)*tmp1(2)*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
        END IF
        !BY (F)
        res_loc(jr,7) = res_loc(jr,7) &
          + basis_vals(jr)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          - self%dt*basis_vals(jr)*tmp1(2)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + self%dt*basis_vals(jr)*dvel(1,1)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + self%dt*basis_vals(jr)*dvel(3,3)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dby)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          - self%dt*basis_vals(jr)*vel(1)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)**2 &
          + self%dt*eta_loc*DOT_PRODUCT(basis_grads(:,jr), dby)*jac_det*quad%wts(m)/(mu0*(coords(1)+gs_epsilon))
      END DO
    END IF

    ! IF WE ARE IN THE PLASMA
    IF (gs_test_bounds(self%eq,coords) .AND. psi >self%eq%plasma_bounds(1)) THEN !check that we are in the plasma
        p_source = self%p_scale*self%eq%P%Fp(psi)*coords(1) 
        f_source = self%f_scale*0.5d0* self%eq%I%fp(psi)/ (coords(1) + gs_epsilon)
        diag=diag+[f_source,p_source]*jac_det*quad%wts(m)
        DO jr=1,oft_blagrange%nce
            res_loc(jr,8) = res_loc(jr,8) &
            - self%dt * basis_vals(jr) * p_source * jac_det*quad%wts(m)
            res_loc(jr,9) = res_loc(jr,9) &
            - self%dt * basis_vals(jr) * f_source * jac_det*quad%wts(m)
        END DO
    END IF
  END DO
    !---Add local values to full vector
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    n_res(cell_dofs(jr)) = n_res(cell_dofs(jr)) + res_loc(jr,1)/self%den_scale
    velx_res(cell_dofs(jr)) = velx_res(cell_dofs(jr)) + res_loc(jr,2)
    vely_res(cell_dofs(jr)) = vely_res(cell_dofs(jr)) + res_loc(jr,3)
    velz_res(cell_dofs(jr)) = velz_res(cell_dofs(jr)) + res_loc(jr,4)
    T_res(cell_dofs(jr)) = T_res(cell_dofs(jr)) + res_loc(jr,5)
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,6)
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr,7)

    pres_vals(cell_dofs(jr)) = pres_vals(cell_dofs(jr)) + res_loc(jr,8)
    alam_vals(cell_dofs(jr)) = alam_vals(cell_dofs(jr)) + res_loc(jr,9)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads, n_weights_loc, vel_weights_loc, T_weights_loc, &
 psi_weights_loc,by_weights_loc, cell_dofs,res_loc)
!$omp end parallel

! Apply BCs
!write(*,*) ALL(self%T_bc)
CALL fem_dirichlet_vec(oft_blagrange,n_weights,n_res,self%n_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(1, :),velx_res,self%velx_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(2, :),vely_res,self%vely_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(3, :),velz_res,self%velz_bc)
CALL fem_dirichlet_vec(oft_blagrange,T_weights,T_res,self%T_bc)
CALL fem_dirichlet_vec(oft_blagrange,by_weights,by_res,self%by_bc)
!write(*,*) 'MIN_VAL ', MINVAL(T_res)
DO i=1,oft_blagrange%nbe
    alam_vals(oft_blagrange%lbe(i))=0.d0
    pres_vals(oft_blagrange%lbe(i))=0.d0
END DO
! RESCALE EQUATIONS --> add some conditions to this?
f_source = self%eq%Itor_target/diag(1)/(1.d0+1.d0/self%eq%Ip_ratio_target)
p_source = self%eq%Itor_target/diag(2)/(self%eq%Ip_ratio_target+1.d0)
psi_res= psi_res + pres_vals*p_source+alam_vals*f_source
self%f_scale=f_source*self%f_scale
self%p_scale=p_source*self%p_scale
diag(1)=diag(1)*f_source
diag(2)=diag(2)*p_source
CALL b%restore_local(n_res,1,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(velx_res,2,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(vely_res,3,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(velz_res,4,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(T_res,5,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(psi_res,6,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(by_res,7,add=.TRUE.)
CALL b%new(ptmp)
CALL self%vac_op%apply(a,ptmp)
CALL b%add(1.d0,1.d0,ptmp)
CALL b%get_local(psi_res, 6)
CALL ptmp%delete
DEALLOCATE(n_res, velx_res, vely_res, velz_res,T_res, psi_res, by_res,pres_vals, alam_vals)
END SUBROUTINE nlfun_apply

SUBROUTINE apply_rhs(self,a,b)
class(gs_xmhd_nlfun), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
LOGICAL :: curved
INTEGER(i4) :: i,m,jr, k,l
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: eta_loc, curr_loc, gamma
REAL(r8) ::  n, dn(3), vel(3), dvel(3,3), T, dT(3), psi, dpsi(3),by, dby(3), coords(3), jac_det, jac_mat(3,4)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals, n_weights_loc, T_weights_loc, psi_weights_loc, by_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads, res_loc, vel_weights_loc
REAL(r8), POINTER, DIMENSION(:) :: n_weights, T_weights, psi_weights, by_weights 
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
REAL(r8), POINTER, DIMENSION(:) :: n_res, velx_res, vely_res, velz_res, T_res, psi_res, by_res, vtmp
quad=>self%eq%fe_rep%quad
NULLIFY( n_weights, vel_weights, T_weights,psi_weights, by_weights, &
n_res, velx_res, vely_res, velz_res, T_res, psi_res, by_res)
!---Get weights from solution vector
ALLOCATE(vel_weights(3,oft_blagrange%ne))
CALL a%get_local(n_weights, 1)
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(T_weights, 5)
CALL a%get_local(psi_weights, 6)
CALL a%get_local(by_weights, 7)
!--Initialize residuals with zeros
CALL b%set(0.d0)
CALL b%get_local(n_res, 1)
CALL b%get_local(velx_res, 2)
CALL b%get_local(vely_res, 3)
CALL b%get_local(velz_res, 4)
CALL b%get_local(T_res, 5)
CALL b%get_local(psi_res, 6)
CALL b%get_local(by_res, 7)

gamma = self%gamma

!$omp parallel private(m,jr,curved,coords,cell_dofs,basis_vals,basis_grads, &
!$omp n_weights_loc, vel_weights_loc, T_weights_loc, psi_weights_loc,by_weights_loc,res_loc,jac_mat, &
!$omp jac_det, n, dn, vel, dvel, T, dT, psi, dpsi, by, dby, &
!$omp eta_loc, curr_loc)
!Allocate local arrays
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(n_weights_loc(oft_blagrange%nce))
ALLOCATE(vel_weights_loc(3, oft_blagrange%nce))
ALLOCATE(T_weights_loc(oft_blagrange%nce))
ALLOCATE(psi_weights_loc(oft_blagrange%nce))
ALLOCATE(by_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce, 7))

!$omp do schedule(static)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function

  ! Set local weights
  n_weights_loc = n_weights(cell_dofs)
  vel_weights_loc = vel_weights(:, cell_dofs)
  T_weights_loc = T_weights(cell_dofs)
  psi_weights_loc = psi_weights(cell_dofs)
  by_weights_loc = by_weights(cell_dofs)
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
    n = 0.d0; dn = 0.d0
    vel = 0.d0; dvel = 0.d0
    T = 0.d0; dT = 0.d0
    psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0

    ! switch from 2D to 3D gradients
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0

    DO jr=1,oft_blagrange%nce
      n = n + n_weights_loc(jr)*basis_vals(jr)
      dn = dn + n_weights_loc(jr)*basis_grads(:,jr)
      vel = vel + vel_weights_loc(:, jr)*basis_vals(jr)
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      T = T + T_weights_loc(jr)*basis_vals(jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      by = by + by_weights_loc(jr)*basis_vals(jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    n = n * self%den_scale
    dn = dn * self%den_scale

    eta_loc = self%eta(mesh%reg(i))
    curr_loc = self%curr(mesh%reg(i))
    DO jr=1,oft_blagrange%nce
      IF(self%region_flag(self%eq%mesh%reg(i)) == 1) THEN
        ! DENSITY
        res_loc(jr,1) = res_loc(jr, 1) &
            + basis_vals(jr)*n*jac_det*quad%wts(m)*coords(1)
        ! VELOCITY
        res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          + basis_vals(jr)*vel*jac_det*quad%wts(m)*coords(1)
        ! TEMPERATURE
        res_loc(jr,5) = res_loc(jr, 5) &
          + basis_vals(jr)*T*jac_det*quad%wts(m)*coords(1)/(gamma-1)
        
        ! B_y (F)
        res_loc(jr,7) = res_loc(jr,7) &
          + basis_vals(jr)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
      END IF
      
      ! PSI
      IF (self%region_flag(mesh%reg(i))==1 .OR. self%region_flag(mesh%reg(i))==3 ) THEN
          res_loc(jr,6) = res_loc(jr,6) &
          + basis_vals(jr)*psi*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
      END IF 
      IF (self%region_flag(mesh%reg(i))==4) THEN
          res_loc(jr,6) = res_loc(jr,6) &
          + basis_vals(jr)*self%dt*curr_loc*jac_det*quad%wts(m)
      END IF 

    END DO
  END DO
    !---Add local values to full vector
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    n_res(cell_dofs(jr)) = n_res(cell_dofs(jr)) + res_loc(jr,1)/self%den_scale
    velx_res(cell_dofs(jr)) = velx_res(cell_dofs(jr)) + res_loc(jr,2)
    vely_res(cell_dofs(jr)) = vely_res(cell_dofs(jr)) + res_loc(jr,3)
    velz_res(cell_dofs(jr)) = velz_res(cell_dofs(jr)) + res_loc(jr,4)
    T_res(cell_dofs(jr)) = T_res(cell_dofs(jr)) + res_loc(jr,5)
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,6)
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr,7)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads, n_weights_loc, vel_weights_loc, T_weights_loc, psi_weights_loc, by_weights_loc, cell_dofs,res_loc)
!$omp end parallel
! SET BOUNDARY CONDITIONS
CALL fem_dirichlet_vec(oft_blagrange,n_weights,n_res,self%n_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(1, :),velx_res,self%velx_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(2, :),vely_res,self%vely_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(3, :),velz_res,self%velz_bc)
CALL fem_dirichlet_vec(oft_blagrange,T_weights,T_res,self%T_bc)
CALL fem_dirichlet_vec(oft_blagrange,by_weights,by_res,self%by_bc)

DO i=1,oft_blagrange%nbe
    psi_res(oft_blagrange%lbe(i))=psi_weights(oft_blagrange%lbe(i))
END DO
CALL b%restore_local(n_res,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velx_res,2,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(vely_res,3,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velz_res,4,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(T_res,5,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(psi_res,6,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(by_res,7,add=.TRUE.)
END SUBROUTINE apply_rhs

SUBROUTINE gs_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a
CALL current_sim%mfmat%update(a)
END SUBROUTINE gs_mfnk_update

SUBROUTINE build_vac_jacobian(self, mat)
class (oft_gs_xmhd_sim), intent(inout) :: self
class (oft_matrix), pointer, intent(inout) :: mat
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals
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
!$omp  jac_loc,jac_mat,jac_det,eta_loc, iloc)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
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
        jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        IF (self%region_flag(self%eq%mesh%reg(i)) == 1 .OR. self%region_flag(self%eq%mesh%reg(i)) == 3) THEN
            jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
            + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
        END IF
      END DO
    END DO
  END DO
  !---Get local to global DOF mapping
  call oft_blagrange%ncdofs(i,cell_dofs)
!---Apply bc to local matrix
  DO jr=1,oft_blagrange%nce
    IF(oft_blagrange%be(cell_dofs(jr))) jac_loc(6,6)%m(jr,:)=0.d0
  END DO
  CALL self%fe_rep%mat_add_local(mat,jac_loc,iloc,tlocks)
END DO
deallocate(cell_dofs,basis_vals,basis_grads,jac_loc, iloc)
!$omp end parallel
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)

! apply free boundary BCs to psi
CALL set_bcmat_mod(self%eq,mat, 6,6)

CALL self%fe_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
END SUBROUTINE build_vac_jacobian

SUBROUTINE build_approx_jacobian(self, mat, a)
class (oft_gs_xmhd_sim), intent(inout) :: self
class (oft_matrix), pointer, intent(inout) :: mat
class(oft_vector), intent(inout) :: a !< Solution for computing jacobian
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals, n_weights_loc, T_weights_loc, psi_weights_loc, by_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads, vel_weights_loc
REAL(r8) :: n, dn(3), vel(3), T, dT(3), psi, dpsi(3), by, dby(3), dvel(3,3), div_vel, btmp(3) !reconstructed variables
REAL (r8) :: coords(3), eta_loc, jac_det, jac_mat(3,4), tmp2(3), tmp3(3)
REAL(r8) :: chi, nu, D_diff, gamma, k_boltz, m_i, B_0(3) ! physics parameters
REAL(r8), POINTER, DIMENSION(:) :: n_weights, T_weights, psi_weights, by_weights
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
REAL(r8), POINTER, DIMENSION(:) ::  vtmp
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
CLASS(oft_vector), POINTER :: oft_lag_vec
INTEGER(i4), ALLOCATABLE, DIMENSION(:), TARGET :: cell_dofs
integer (i4) :: i, jr, jc, m, k, l
type(oft_quad_type), pointer :: quad
type(oft_1d_int), allocatable, dimension(:) :: iloc
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
LOGICAL :: curved
quad=>oft_blagrange%quad
CALL mat%zero

NULLIFY(n_weights,vel_weights, T_weights, &
         psi_weights, by_weights, vtmp)
!---Get weights from solution vector
CALL a%get_local(n_weights,1)
ALLOCATE(vel_weights(3,oft_blagrange%ne))
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(T_weights,5)
CALL a%get_local(psi_weights,6)
CALL a%get_local(by_weights,7)
!--Set local physics parameters
chi = self%chi
nu = self%nu
D_diff = self%D_diff
gamma = self%gamma
k_boltz = self%k_boltz
m_i = self%m_i
B_0 = self%B_0

!--Setup thread locks
ALLOCATE(tlocks(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
  call omp_init_lock(tlocks(i))
END DO

! Declare variables private for OMP
!$omp parallel private(m,jr,jc,curved,cell_dofs,basis_vals,basis_grads, &
!$omp  jac_loc,jac_mat,jac_det,eta_loc, &
!$omp n_weights_loc, vel_weights_loc, T_weights_loc, psi_weights_loc, by_weights_loc, &
!$omp n, dn, vel, dvel, div_vel, T, dT, psi, dpsi, by, dby, iloc)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce))
ALLOCATE(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
ALLOCATE(iloc(self%fe_rep%nfields))
ALLOCATE(n_weights_loc(oft_blagrange%nce),vel_weights_loc(3, oft_blagrange%nce),&
        T_weights_loc(oft_blagrange%nce), psi_weights_loc(oft_blagrange%nce),&
        by_weights_loc(oft_blagrange%nce))
DO i=1,self%fe_rep%nfields
   iloc(i)%v=>cell_dofs
END DO
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do schedule(static)ordered
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix

  ! Set local weights
  n_weights_loc = n_weights(cell_dofs)
  vel_weights_loc = vel_weights(:, cell_dofs)
  T_weights_loc = T_weights(cell_dofs)
  psi_weights_loc = psi_weights(cell_dofs)
  by_weights_loc = by_weights(cell_dofs)
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
    !---Reconstruct values of solution fields
    n = 0.d0; dn = 0.d0; vel = 0.d0; dvel = 0.d0
    T = 0.d0; dT = 0.d0; psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0
    DO jr=1,oft_blagrange%nce
      n = n + n_weights_loc(jr)*basis_vals(jr)
      vel = vel + vel_weights_loc(:, jr)*basis_vals(jr)
      T = T + T_weights_loc(jr)*basis_vals(jr)
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      by = by + by_weights_loc(jr)*basis_vals(jr)
      dn = dn + n_weights_loc(jr)*basis_grads(:,jr)
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    n = n * self%den_scale
    dn = dn * self%den_scale
    div_vel = dvel(1,1) +vel(1)/(coords(1)+gs_epsilon) + dvel(3,3)
    btmp = cross_product(dpsi/coords(1), [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0]/(coords(1)+gs_epsilon) + B_0

    !---Compute local matrix contributions
    ! VACUUM TERMS
    DO jr=1,oft_blagrange%nce
      DO jc=1,oft_blagrange%nce
        ! psi, psi
        jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        IF (self%region_flag(self%eq%mesh%reg(i)) == 1 .OR. self%region_flag(self%eq%mesh%reg(i)) == 3) THEN
            jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
            + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
        END IF
      END DO
    END DO
  ! EVERYTHING ELSE
    IF (self%region_flag(self%eq%mesh%reg(i)) == 1) THEN
      DO jr=1,oft_blagrange%nce
        DO jc=1,oft_blagrange%nce
        ! !n, n
          jac_loc(1, 1)%m(jr,jc) = jac_loc(1, 1)%m(jr, jc) &
          + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)*coords(1) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(basis_grads(:, jc), vel)*jac_det*quad%wts(m)*coords(1) & !delta_n*div(u)
          + self%dt*basis_vals(jr)*basis_vals(jc)*div_vel*jac_det*quad%wts(m)*coords(1) & ! u dot grad(delta_n)
          + self%dt*D_diff*DOT_PRODUCT(basis_grads(:, jr),basis_grads(:, jc))*jac_det*quad%wts(m)*coords(1)
        !n,vel
          DO l=1,3
            jac_loc(1, l+1)%m(jr,jc) = jac_loc(1, l+1)%m(jr, jc) &
            + basis_vals(jr)*self%dt*n*basis_grads(l, jc)*jac_det*quad%wts(m)*coords(1) & ! n*div(delta_u) = n*SUM(basis_grads)?
            + basis_vals(jr)*self%dt*basis_vals(jc)*dn(l)*jac_det*quad%wts(m)*coords(1) ! u dot grad(n)
          END DO
          jac_loc(1, 2)%m(jr,jc) = jac_loc(1, 2)%m(jr, jc) &
          + basis_vals(jr)*self%dt*n*basis_vals(jc)*jac_det*quad%wts(m)
        !vel,n
          DO l=1,3
            jac_loc(l+1,1)%m(jr,jc) = jac_loc(l+1,1)%m(jr,jc) &
            + self%dt*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dT(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n) &
            + self%dt*basis_vals(jr)*2.d0*k_boltz*T*basis_grads(l,jc)*jac_det*quad%wts(m)*coords(1)/(m_i*n) &
            - self%dt*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*n*dT(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n**2.d0) &
            - self%dt*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dn(l)*T*jac_det*quad%wts(m)*coords(1)/(m_i*n**2.d0) &
            - self%dt*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), btmp)*btmp(l)*jac_det*quad%wts(m)*coords(1)/(mu0*m_i*n**2.d0) &
            - self%dt*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), btmp)*btmp(l)*jac_det*quad%wts(m)*coords(1)/(mu0*m_i*n**2.d0) &
            + self%dt*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn, btmp)*btmp(l)*jac_det*quad%wts(m)*coords(1)/(mu0*m_i*n**3.d0) &
            + self%dt*basis_vals(jc)*basis_grads(l,jr)*DOT_PRODUCT(btmp, btmp)*jac_det*quad%wts(m)*coords(1)/(2.d0*mu0*m_i*n**2.d0) &
            + self%dt*basis_vals(jr)*basis_grads(l,jc)*DOT_PRODUCT(btmp, btmp)*jac_det*quad%wts(m)*coords(1)/(2.d0*mu0*m_i*n**2.d0) &
            - self%dt*basis_vals(jr)*basis_vals(jc)*dn(l)*DOT_PRODUCT(btmp, btmp)*jac_det*quad%wts(m)*coords(1)/(mu0*m_i*n**3.d0) &
            - self%dt*nu*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), dvel(l, :))*jac_det*quad%wts(m)*coords(1)/(m_i*n**2) & !-- not sure if indexing on dvel is right here
            - self%dt*nu*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), dvel(l, :))*jac_det*quad%wts(m)*coords(1)/(m_i*n**2) &
            + self%dt*nu*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn, dvel(l, :))*jac_det*quad%wts(m)*coords(1)/(m_i*n**3)
          END DO
          jac_loc(2,1)%m(jr,jc) = jac_loc(2,1)%m(jr,jc) &
          - self%dt*basis_vals(jr)*basis_vals(jc)*(btmp(2)**2-btmp(1)**2-btmp(3)**2)*jac_det*quad%wts(m)/(mu0*m_i*n**2) &
          - self%dt*basis_vals(jr)*basis_vals(jc)*nu*vel(1)*jac_det*quad%wts(m)/(m_i*n**2*(coords(1)+gs_epsilon))
          jac_loc(3,1)%m(jr,jc) = jac_loc(3,1)%m(jr,jc) &
          + self%dt*basis_vals(jr)*basis_vals(jc)*(btmp(1)*btmp(2))*jac_det*quad%wts(m)/(mu0*m_i*n**2) &
          - self%dt*basis_vals(jr)*basis_vals(jc)*nu*vel(2)*jac_det*quad%wts(m)/(m_i*n**2*(coords(1)+gs_epsilon))
  
          ! vel,vel
          DO k=1,3
            jac_loc(k+1,k+1)%m(jr,jc)= jac_loc(k+1,k+1)%m(jr,jc) &
              + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m) &
              + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, basis_grads(:,jc))*jac_det*quad%wts(m)*coords(1) &
              + self%dt*nu*DOT_PRODUCT(basis_grads(:,jr), basis_grads(:,jc))*jac_det*quad%wts(m)*coords(1)/(m_i*n) &
              - self%dt*nu*basis_vals(jr)*DOT_PRODUCT(dn, basis_grads(:,jc))*jac_det*quad%wts(m)*coords(1)/(m_i*n**2)
            DO l=1,3
              jac_loc(k+1,l+1)%m(jr,jc)= jac_loc(k+1,l+1)%m(jr,jc) &
              + self%dt*basis_vals(jr)*basis_vals(jc)*dvel(k, l)*jac_det*quad%wts(m)*coords(1)
            END DO
          END DO
          jac_loc(2,2)%m(jr,jc)= jac_loc(2,2)%m(jr,jc) &
            + self%dt*nu*basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)/(m_i*n*(coords(1)+gs_epsilon))
          jac_loc(2,3)%m(jr,jc)= jac_loc(2,3)%m(jr,jc) &
            -self%dt*basis_vals(jr)*2.d0*vel(2)*basis_vals(jc)*jac_det*quad%wts(m)
          jac_loc(3,2)%m(jr,jc)= jac_loc(3,2)%m(jr,jc) &
            +self%dt*basis_vals(jr)*basis_vals(jc)*vel(2)*jac_det*quad%wts(m)
          jac_loc(3,3)%m(jr,jc)= jac_loc(3,3)%m(jr,jc) &
            + self%dt*basis_vals(jr)*basis_vals(jc)*vel(1)*jac_det*quad%wts(m) &
            + self%dt*nu*basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)/(m_i*n*(coords(1)+gs_epsilon))
          ! vel, T
          DO l=1,3
            jac_loc(l+1,5)%m(jr,jc) = jac_loc(l+1,5)%m(jr,jc) &
            + self%dt*basis_vals(jr)*2*k_boltz*dn(l)*basis_vals(jc)*jac_det*quad%wts(m)*coords(1)/(m_i*n) &
            + self%dt*basis_vals(jr)*2*k_boltz*basis_grads(l,jc)*jac_det*quad%wts(m)*coords(1)/(m_i) 
          END DO
          !vel, psi
          tmp2 = cross_product(basis_grads(:,jc), [0.d0,1.d0/(coords(1)+gs_epsilon),0.d0]) ! this is 'dB'
          DO l=1,3
            jac_loc(l+1,6)%m(jr,jc) = jac_loc(l+1,6)%m(jr,jc) &
            + self%dt*DOT_PRODUCT(basis_grads(:,jr), tmp2)*btmp(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n*mu0) &
            + self%dt*DOT_PRODUCT(basis_grads(:,jr), btmp)*tmp2(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n*mu0) &
            - self%dt*basis_grads(l,jr)*DOT_PRODUCT(tmp2, btmp)*jac_det*quad%wts(m)*coords(1)/(m_i*n*mu0) &
            - self%dt*basis_vals(jr)*DOT_PRODUCT(dn, tmp2)*btmp(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n**2*mu0) &
            - self%dt*basis_vals(jr)*DOT_PRODUCT(dn, btmp)*tmp2(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n**2*mu0) &
            + self%dt*basis_vals(jr)*dn(l)*DOT_PRODUCT(tmp2, btmp)*jac_det*quad%wts(m)*coords(1)/(m_i*n**2*mu0) 
          END DO
          jac_loc(2,6)%m(jr,jc) = jac_loc(2,6)%m(jr,jc) &
          + self%dt*basis_vals(jr)*(btmp(2)*tmp2(2)-btmp(1)*tmp2(1)-btmp(3)*tmp2(3))*jac_det*quad%wts(m)/(m_i*n*mu0)
          jac_loc(3,6)%m(jr,jc) = jac_loc(3,6)%m(jr,jc) &
          - self%dt*basis_vals(jr)*(btmp(1)*tmp2(2)+btmp(2)*tmp2(1))*jac_det*quad%wts(m)/(m_i*n*mu0)
          !vel, by
          tmp2 = [0.d0,basis_vals(jc)/(coords(1)+gs_epsilon),0.d0]! this is 'delta B_y'
          DO l=1,3
            jac_loc(l+1,7)%m(jr,jc) = jac_loc(l+1,7)%m(jr,jc) &
            + self%dt*DOT_PRODUCT(basis_grads(:,jr),tmp2)*btmp(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n*mu0) &
            + self%dt*DOT_PRODUCT(basis_grads(:,jr), btmp)*tmp2(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n*mu0) &
            - self%dt*basis_grads(l,jr)*DOT_PRODUCT(tmp2, btmp)*jac_det*quad%wts(m)*coords(1)/(m_i*n*mu0) &
            - self%dt*basis_vals(jr)*DOT_PRODUCT(dn, tmp2)*btmp(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n**2*mu0) &
            - self%dt*basis_vals(jr)*DOT_PRODUCT(dn, btmp)*tmp2(l)*jac_det*quad%wts(m)*coords(1)/(m_i*n**2*mu0) &
            + self%dt*basis_vals(jr)*dn(l)*DOT_PRODUCT(tmp2, btmp)*jac_det*quad%wts(m)*coords(1)/(m_i*n**2*mu0)
          END DO
          jac_loc(2,7)%m(jr,jc) = jac_loc(2,7)%m(jr,jc) &
          + self%dt*basis_vals(jr)*btmp(2)*tmp2(2)*jac_det*quad%wts(m)/(m_i*n*mu0)
          jac_loc(3,7)%m(jr,jc) = jac_loc(3,7)%m(jr,jc) &
          - self%dt*basis_vals(jr)*btmp(1)*tmp2(2)*jac_det*quad%wts(m)/(m_i*n*mu0)
          !T, n
          jac_loc(5, 1)%m(jr,jc) = jac_loc(5, 1)%m(jr, jc) &  
          - self%dt*chi*basis_vals(jr)*DOT_PRODUCT(basis_grads(:, jc), dT)*jac_det*quad%wts(m)*coords(1)/n & ! grad(delta_n) 
          + self%dt*chi*basis_vals(jr)*basis_vals(jc)*DOT_PRODUCT(dn, dT)*jac_det*quad%wts(m)*coords(1)/(n**2) ! delta_n 
          ! T, vel
          DO l=1,3
            jac_loc(5, l+1)%m(jr,jc) = jac_loc(5, l+1)%m(jr, jc) &  
            + basis_vals(jr) * self%dt*basis_vals(jc)*dT(l)*jac_det*quad%wts(m)*coords(1)/(gamma-1) & ! delta_u dot Delta_T
            + basis_vals(jr) * self%dt*T*basis_grads(l, jc)*jac_det*quad%wts(m)*coords(1) ! div(delta u) = SUM(basis_grads)?
          END DO
          jac_loc(5, 2)%m(jr,jc) = jac_loc(5, 2)%m(jr, jc) &
          + basis_vals(jr)*self%dt*T*basis_vals(jc)*jac_det*quad%wts(m)
          !T, T
          jac_loc(5, 5)%m(jr,jc) = jac_loc(5,5)%m(jr, jc) &  
          + basis_vals(jr) * basis_vals(jc)*jac_det*quad%wts(m)*coords(1)/(gamma-1) & ! delta_T
          + basis_vals(jr) * self%dt*DOT_PRODUCT(vel, basis_grads(:, jc))*jac_det*quad%wts(m)*coords(1)/(gamma-1) & ! nabla(dT)
          + basis_vals(jr) * self%dt*basis_vals(jc)*div_vel*jac_det*quad%wts(m)*coords(1) & ! dT != nabla(dT)
          + self%dt * chi * DOT_PRODUCT(basis_grads(:, jc),basis_grads(:, jr))*jac_det*quad%wts(m)*coords(1) & ! dT_Chi
          - self%dt*basis_vals(jr)*chi*DOT_PRODUCT(dn, basis_grads(:, jc))*jac_det*quad%wts(m)*coords(1)/n 
          ! psi, psi
          jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
            + self%dt*basis_vals(jr)*DOT_PRODUCT(vel,basis_grads(:,jc))*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
          !psi, vel
          DO l=1,3
            tmp2 = [0.d0, 0.d0,0.d0]
            jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
            + self%dt*basis_vals(jr)*basis_vals(jc)*dpsi(l)*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
            tmp2(l) = 1.d0
            tmp3 = cross_product(B_0, tmp2)
            jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
            + self%dt*basis_vals(jr)*basis_vals(jc)*tmp3(2)*jac_det*quad%wts(m)/(eta_loc*(coords(1)+gs_epsilon))
          END DO
          ! !By, By
          ! jac_loc(7, 7)%m(jr,jc) = jac_loc(7, 7)%m(jr,jc) &
          ! + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          ! + basis_vals(jr)*self%dt*dvel(1,1)*basis_vals(jc)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          ! + basis_vals(jr)*self%dt*dvel(1,1)*basis_vals(jc)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          ! + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, basis_grads(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          ! - basis_vals(jr)*self%dt*vel(1)*basis_vals(jc)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)**2 &
          ! + self%dt*eta_loc*DOT_PRODUCT(basis_grads(:,jr), basis_grads(:,jc))*jac_det*quad%wts(m)/(mu0*(coords(1)+gs_epsilon))
          ! !By, vel
          ! jac_loc(7,2)%m(jr,jc) = jac_loc(7,2)%m(jr,jc) &
          ! + basis_vals(jr)*self%dt*basis_vals(jc)*dby(1)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          ! - basis_vals(jr)*self%dt*basis_vals(jc)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)**2 &
          ! + basis_vals(jr)*self%dt*basis_grads(1,jc)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
          ! jac_loc(7,4)%m(jr,jc) = jac_loc(7,4)%m(jr,jc) &
          ! + basis_vals(jr)*self%dt*basis_vals(jc)*dby(3)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          ! + basis_vals(jr)*self%dt*basis_grads(3,jc)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
          ! jac_loc(7,3)%m(jr,jc) = jac_loc(7,3)%m(jr,jc) &
          !  - basis_vals(jr)*self%dt*tmp2(2)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
          ! !By, psi
          ! tmp2 = cross_product(basis_grads(:,jc),dvel(2,:))
          ! jac_loc(7, 6)%m(jr,jc) = jac_loc(7, 6)%m(jr,jc) &
          !  - basis_vals(jr)*self%dt*tmp2(2)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        END DO
      END DO
      ! write(*,*) 'MAX before', MAXVAL(jac_loc(1,1)%m)
      ! write(*,*) 'MIN before', MINVAL(jac_loc(1,1)%m)
    END IF
  END DO
  DO jr = 1,7
    jac_loc(1, jr)%m = jac_loc(1, jr)%m / self%den_scale
    jac_loc(jr, 1)%m = jac_loc(jr, 1)%m * self%den_scale
  END DO
 !write(*,*) 'hi'
 ! IF(MAXVAL(jac_loc(1,2)%m) > 1.d-28) write(*,*) MAXVAL(jac_loc(1,2)%m)
!---Apply bc to local matrix
  DO jr=1,oft_blagrange%nce
    IF(oft_blagrange%be(cell_dofs(jr))) jac_loc(6,6)%m(jr,:)=0.d0
  END DO
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%n_bc(cell_dofs),1)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velx_bc(cell_dofs),2)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%vely_bc(cell_dofs),3)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velz_bc(cell_dofs),4)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%T_bc(cell_dofs),5)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%by_bc(cell_dofs),7)
  CALL self%fe_rep%mat_add_local(mat,jac_loc,iloc,tlocks)
END DO
deallocate(cell_dofs,basis_vals,basis_grads,jac_loc, iloc)
!$omp end parallel
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
! apply free boundary BCs to psi
CALL set_bcmat_mod(self%eq,mat, 6,6)
CALL fem_dirichlet_diag(oft_blagrange,mat,self%n_bc,1)
CALL fem_dirichlet_diag(oft_blagrange,mat,self%velx_bc,2)
CALL fem_dirichlet_diag(oft_blagrange,mat,self%vely_bc,3)
CALL fem_dirichlet_diag(oft_blagrange,mat,self%velz_bc,4)
CALL fem_dirichlet_diag(oft_blagrange,mat,self%T_bc,5)
CALL fem_dirichlet_diag(oft_blagrange,mat,self%by_bc,7)

CALL self%fe_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec, n_weights, vel_weights, T_weights, psi_weights, by_weights)

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
      ALLOCATE(graphs(i,j)%g)
      CALL graph_deep_copy(known_graphs(k)%g, graphs(i,j)%g)
      !graphs(i,j)%g=>known_graphs(k)%g
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
      !graphs(i,j)%g=>known_graphs(nknown_graphs)%g
      ALLOCATE(graphs(i,j)%g)
      CALL graph_deep_copy(known_graphs(k)%g, graphs(i,j)%g)
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
      DEALLOCATE(dense_flag, bc_nodes)
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
!> copy graph rather than using pointers
!------------------------------------------------------------------------------
SUBROUTINE graph_deep_copy(src, dest)
  TYPE(oft_graph), INTENT(in)  :: src
  TYPE(oft_graph), INTENT(out) :: dest
  ! allocate and copy all fields
  dest%nr  = src%nr
  dest%nrg = src%nrg
  dest%nc  = src%nc
  dest%ncg = src%ncg
  dest%nnz = src%nnz
  ALLOCATE(dest%kr(SIZE(src%kr)))
  ALLOCATE(dest%lc(SIZE(src%lc)))
  dest%kr = src%kr
  dest%lc = src%lc
END SUBROUTINE

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

subroutine apply_mhd_bcs(self,cell_ind, cell_dofs)
class(oft_gs_xmhd_sim), intent(inout) :: self
integer(i4) , intent(in) :: cell_ind
INTEGER(i4), POINTER, DIMENSION(:), intent(inout) :: cell_dofs
INTEGER(i4) :: j
call oft_blagrange%ncdofs(cell_ind,cell_dofs) ! Get global index of local DOFs
DO j=1, SIZE(cell_dofs)
  self%by_bc(cell_dofs(j)) = .TRUE. ! prevent psi evolution in superconductor
END DO
end subroutine apply_mhd_bcs

!---------------------------------------------------------------------------
!> Apply boundary conditions for non-extended MHD regions (plasma, coils, solid conductors, vacuum)
!---------------------------------------------------------------------------
subroutine apply_bcs(self,cell_ind, cell_dofs)
class(oft_gs_xmhd_sim), intent(inout) :: self
INTEGER(i4) , intent(in) :: cell_ind
INTEGER(i4), POINTER, DIMENSION(:), intent(inout) :: cell_dofs
INTEGER(i4) :: j
call oft_blagrange%ncdofs(cell_ind,cell_dofs) ! Get global index of local DOFs
DO j=1, SIZE(cell_dofs)
  self%n_bc(cell_dofs(j)) = .TRUE. ! prevent density evolution in solid conductor
  self%velx_bc(cell_dofs(j)) = .TRUE. ! prevent velocity evolution in solid conductor
  self%vely_bc(cell_dofs(j)) = .TRUE. ! prevent velocity evolution in solid conductor
  self%velz_bc(cell_dofs(j)) = .TRUE. ! prevent velocity evolution in solid conductor
  self%T_bc(cell_dofs(j)) = .TRUE. ! prevent temperature evolution in solid conductor
  self%by_bc(cell_dofs(j)) = .TRUE. ! prevent psi evolution in superconductor
END DO
end subroutine apply_bcs

END MODULE gs_xmhd