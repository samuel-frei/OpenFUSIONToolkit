!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file gs_xmhd.F90
!
!> Solve coupled non-linear grad-shafranov evolution and extended MHD
!---------------------------------------------------------------------------
MODULE gs_xmhd_v8
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
  INTEGER(i4) :: lim_ind = 1 !< Needs docs
  LOGICAL :: pm = .FALSE.
  ! PHYSICS PARAMETERS
  REAL(r8) :: nu = -1.d0 !< Needs docs
  REAL(r8) :: k_boltz=elec_charge
  REAL(r8) :: rho=-1.d0
  REAL(r8) :: den_scale = 1.d19 !< Needs docs
  REAL (r8) :: B_0(3) = 0.d0
  REAL (r8) :: lim_vac_int = 0.d0
  REAL (r8) :: tflux_source = 0.d0
  CLASS(bfem_interp), POINTER :: j_source => NULL() !< Interpolator for current source term
  REAL (r8) :: j_source_scale = 1.d0 !< Scale factor for current source term
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta_t
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta_p
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta_node
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: curr
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: region_flag
  LOGICAL :: evolve_F = .TRUE.
  TYPE(gs_eq), POINTER :: eq => NULL() !< Equilibrium object
  ! BOUNDARY CONDITIONS
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: p_bc => NULL() !< n BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velx_bc => NULL() !< velx BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: vely_bc => NULL() !< vely BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velz_bc => NULL() !< velz BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: by_bc => NULL() !< By (F) BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< By (F) BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: plasma_bc => NULL() !< plasma BC flag for F
  ! SOLVER OBJECTS
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
  TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  TYPE(xdmf_plot_file) :: xdmf_plot
  TYPE(xdmf_plot_file) :: xdmf_plot_1
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
    !> Run single timestep
  PROCEDURE :: add_timestep => add_timestep
  !> Save restart file
  PROCEDURE :: rst_save => rst_save
  !> Load restart file
  PROCEDURE :: rst_load => rst_load
END TYPE oft_gs_xmhd_sim

TYPE, extends(oft_noop_matrix) :: gs_xmhd_nlfun
  ! SOLVER PARAMETERS
  REAL(r8) :: dt = -1.d0 !< Time step
  ! PHYSICS PARAMETERS
  REAL(r8) :: nu = -1.d0 !< Needs docs
  REAL(r8) :: k_boltz=elec_charge
  REAL(r8) :: rho= -1.d0
  INTEGER(i4) :: lim_ind = 1 !< Needs docs
  REAL (r8) :: B_0(3) = 0.d0
  REAL (r8) :: lim_vac_int = 0.d0
  REAL (r8) :: tflux_source = 0.d0
  CLASS(bfem_interp), POINTER :: j_source => NULL() !< Interpolator for current source term
  REAL (r8) :: j_source_scale = 1.d0 !< Scale factor for current source term
  TYPE(gs_eq), POINTER :: eq => NULL() !< Equilibrium object
  REAL(r8) :: f_scale = 1.d0 !< Scale factor for \f$ F*F' \f$ term
  REAL(r8) :: p_scale = 1.d0 !< Scale factor for \f$ P' \f$ term
  REAL(r8) :: diag_vals(2) = 0.d0 !< Used to determine f and p scales
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta_t
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta_p
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: eta_node
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: curr
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: region_flag
  LOGICAL :: evolve_F = .TRUE.
  ! BOUNDARY CONDITIONS
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: p_bc => NULL() !< n BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velx_bc => NULL() !< velx BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: vely_bc => NULL() !< vely BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velz_bc => NULL() !< velz BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: T_bc => NULL() !< T BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: by_bc => NULL() !< By (F) BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< By (F) BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: plasma_bc => NULL() !< plasma BC flag for F
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
CLASS(multigrid_mesh), POINTER :: mg_mesh_1 => NULL()
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh_1 => NULL()
TYPE(oft_ml_fem_type), TARGET, PUBLIC :: ML_oft_blagrange_1 ! 1st order finite element basis
TYPE(oft_ml_fem_type), TARGET, PUBLIC :: ML_oft_blagrange_2 ! 2nd order finite element basis
CLASS(oft_scalar_bfem), POINTER :: oft_blagrange_1 => NULL()
CLASS(oft_scalar_bfem), POINTER :: oft_blagrange_2 => NULL()

CONTAINS

subroutine setup(self, mg_mesh_in, mg_mesh_1_in)
class(oft_gs_xmhd_sim), intent(inout), target :: self !< NL operator object
CLASS(multigrid_mesh), TARGET, intent(in) :: mg_mesh_in
CLASS(multigrid_mesh), TARGET, intent(in) :: mg_mesh_1_in ! for use in plotting lower order fields
integer(i4) :: ierr, i, type, order, order_p
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr
CLASS(oft_native_matrix), POINTER :: A_native
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs_1, cell_dofs_2
LOGICAL, ALLOCATABLE :: p_dir_set(:)
!------------------------------------------------------------------------------
! Setup mesh and finite element representation
!------------------------------------------------------------------------------
order = 2
current_sim=>self
mg_mesh=>mg_mesh_in
mesh=>mg_mesh%smesh
mg_mesh_1=> mg_mesh_1_in
mesh_1=> mg_mesh_1%smesh
!---Setup 1st order FE representation
order_p = order-1
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building 1st order lagrange FE space'
CALL oft_lag_setup(mg_mesh,order_p, ML_blag_obj=ML_oft_blagrange_1,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(oft_blagrange_1,ML_oft_blagrange_1%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)

!---Setup 2nd order FE representation
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building 2nd order lagrange FE space'
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_oft_blagrange_2,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(oft_blagrange_2,ML_oft_blagrange_2%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
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
  ALLOCATE(cell_dofs_1(oft_blagrange_1%nce))
  ALLOCATE(cell_dofs_2(oft_blagrange_2%nce))
  ALLOCATE(self%p_bc(oft_blagrange_1%ne)); self%p_bc=.TRUE.
  ALLOCATE(self%velx_bc(oft_blagrange_2%ne)); self%velx_bc=.FALSE.
  ALLOCATE(self%vely_bc(oft_blagrange_2%ne)); self%vely_bc=.TRUE.
  ALLOCATE(self%velz_bc(oft_blagrange_2%ne)); self%velz_bc=.FALSE.
  ALLOCATE(self%by_bc(oft_blagrange_2%ne)); self%by_bc=.FALSE.  ! FOR NOW WE'RE NOT EVOLVING By (F)
  ALLOCATE(self%psi_bc(oft_blagrange_2%ne)); self%psi_bc=.FALSE. 
  ALLOCATE(self%plasma_bc(oft_blagrange_2%ne)); self%plasma_bc=.FALSE.  ! FOR NOW WE'RE NOT EVOLVING By (F)
  IF (SIZE(self%region_flag) /= mesh%nreg) THEN
    CALL oft_abort("Number of region flags does not match number of regions.","setup",__FILE__)
  END IF
  DO i=1, mesh%nc
    type = self%region_flag(mesh%reg(i))
    IF (type == 1) THEN
      CALL apply_mhd_bcs(self, i, cell_dofs_1, cell_dofs_2)
    ELSE IF(type ==5 .OR. type ==6) THEN
      CALL apply_plasma_bcs(self, i, cell_dofs_1, cell_dofs_2)
    ELSE IF (type >1 .AND. type < 5) THEN
      CALL apply_bcs(self, i, cell_dofs_1, cell_dofs_2)
    ELSE
      CALL oft_abort("Invalid region flag.","setup",__FILE__)
    END IF
  END DO
  ! Pin one node per MHD region after all BCs are set
  ALLOCATE(p_dir_set(mesh%nreg))
  p_dir_set = .FALSE.
  DO i=1, mesh%nc
    type = self%region_flag(mesh%reg(i))
    IF (type == 1 .AND. .NOT. p_dir_set(mesh%reg(i))) THEN
      call oft_blagrange_1%ncdofs(i,cell_dofs_1)
      self%p_bc(cell_dofs_1(1)) = .TRUE.
      p_dir_set(mesh%reg(i)) = .TRUE.
      write(*,*) "Pinning node ", cell_dofs_1(1), " in region ", mesh%reg(i)
    END IF
  END DO
END IF
IF (.NOT. self%evolve_F) THEN
  self%by_bc = .TRUE.
END IF 
DEALLOCATE(cell_dofs_1, cell_dofs_2)
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
CALL self%xdmf_plot%setup("gs_xmhd", "all")
CALL mesh%setup_io(self%xdmf_plot,order)
CALL self%xdmf_plot_1%setup("gs_xmhd_1", "pressure")
CALL mesh_1%setup_io(self%xdmf_plot_1,order_p)
!------------------------------------------------------------------------------
! Create nonlinear operator and set it up
!------------------------------------------------------------------------------
ALLOCATE(self%nlfun) !CHECKBACK
self%nlfun%dt=self%dt
self%nlfun%eq => self%eq
IF(ASSOCIATED(self%eq)) THEN
  self%nlfun%f_scale=self%eq%alam
  self%nlfun%p_scale=self%eq%pnorm
END IF
self%nlfun%nu = self%nu
self%nlfun%rho = self%rho
self%nlfun%B_0 = self%B_0
self%nlfun%evolve_F = self%evolve_F
IF (ALLOCATED(self%eta_node)) THEN
  ALLOCATE(self%nlfun%eta_node(mesh%np))
  self%nlfun%eta_node=self%eta_node
END IF
ALLOCATE(self%nlfun%eta_t(mesh%nreg))
ALLOCATE(self%nlfun%eta_p(mesh%nreg))
self%nlfun%eta_t=self%eta_t
self%nlfun%eta_p=self%eta_p
ALLOCATE(self%nlfun%curr(mesh%nreg))
self%nlfun%curr = self%curr
ALLOCATE(self%nlfun%region_flag(mesh%nreg))
self%nlfun%region_flag = self%region_flag
self%nlfun%lim_vac_int = self%lim_vac_int
self%nlfun%lim_ind = self%lim_ind
self%nlfun%tflux_source = self%tflux_source

self%nlfun%p_bc=>self%p_bc
self%nlfun%velx_bc=>self%velx_bc
self%nlfun%vely_bc=>self%vely_bc
self%nlfun%velz_bc=>self%velz_bc
self%nlfun%by_bc=>self%by_bc
self%nlfun%psi_bc=>self%psi_bc
self%nlfun%plasma_bc => self%plasma_bc

self%nlfun%j_source => self%j_source
self%nlfun%j_source_scale = self%j_source_scale
!------------------------------------------------------------------------------
! Create Solver fields
!------------------------------------------------------------------------------
ALLOCATE(self%fe_rep) !CHECKBACK
self%fe_rep%nfields=6
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields)) 
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>oft_blagrange_1
self%fe_rep%field_tags(1)='p'
self%fe_rep%fields(1)%fe%type = 1
self%fe_rep%fields(2)%fe=>oft_blagrange_2
self%fe_rep%field_tags(2)='velx'
self%fe_rep%fields(2)%fe%type = 2
self%fe_rep%fields(3)%fe=>oft_blagrange_2
self%fe_rep%field_tags(3)='vely'
self%fe_rep%fields(3)%fe%type = 2
self%fe_rep%fields(4)%fe=>oft_blagrange_2
self%fe_rep%field_tags(4)='velz'
self%fe_rep%fields(4)%fe%type = 2
self%fe_rep%fields(5)%fe=>oft_blagrange_2
self%fe_rep%field_tags(5)='by'
self%fe_rep%fields(5)%fe%type = 2
self%fe_rep%fields(6)%fe=>oft_blagrange_2
self%fe_rep%field_tags(6)='psi'
self%fe_rep%fields(6)%fe%type = 2
CALL self%fe_rep%vec_create(self%u)
call self%fe_rep%vec_create(self%rhs)
call self%fe_rep%vec_create(self%tmp)
! CALL self%u%set(1000.d0, 1)
CALL self%u%set(1000.d0, 1)
!CALL self%u%set(0.d0, 1)
CALL self%u%set(0.d0, 2)
CALL self%u%set(0.d0, 3)
CALL self%u%set(0.d0, 4)
CALL self%u%set(0.d0, 5)
CALL self%u%set(0.d0, 6)
! CALL self%rst_load(self%u,'gs_xmhd_00039.rst', 'U')
NULLIFY(tmp_arr)
!------------------------------------------------------------------------------
! Build Jacobian matrix
!------------------------------------------------------------------------------
ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
self%jacobian_block_mask=1
IF(ASSOCIATED(self%eq)) THEN
  self%jacobian_block_mask(6,6) = 3 ! Add dense regions for the boundary in the psi/psi matrix
END IF
IF (self%evolve_F) THEN
  self%jacobian_block_mask(5,5) = 4 ! Add dense regions for F matrix
END IF
CALL fem_mat_create_mod(self%fe_rep, self%nlfun%jac_op, self%jacobian_block_mask)
CALL fem_mat_create_mod(self%fe_rep, self%nlfun%vac_op, self%jacobian_block_mask) 
self%jacobian_block_mask(6,6) = 1 !Set back to normal for creating local matrices
self%jacobian_block_mask(5,5) = 1 !Set back to normal for creating local matrices
CALL build_vac_jacobian(self,self%nlfun%vac_op)
! Preconditioner should use approximate jacobian
ALLOCATE(self%pre) !CHECKBACK
self%pre%A=>self%nlfun%jac_op 
!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat) !CHECKBACK
!CALL self%mfmat%setup(self%rhs,self%nlfun)
self%mfmat%f=>self%nlfun
CALL self%rhs%new(self%mfmat%u0)
CALL self%rhs%new(self%mfmat%f0)
CALL self%rhs%new(self%mfmat%tmp)
CALL self%rhs%new(self%mfmat%utyp)


ALLOCATE(self%mf_solver) !CHECKBACK
self%mfmat%b0=1.d-5
self%mf_solver%A=>self%mfmat
self%mf_solver%its=300
self%mf_solver%nrits=40
self%mf_solver%atol=self%lin_tol
self%mf_solver%itplot=1
oft_env%pm = self%pm
self%mf_solver%pm=oft_env%pm
NULLIFY(self%mf_solver%pre)
IF(ASSOCIATED(self%xml_pre_def))THEN
  CALL create_solver_xml(self%mf_solver%pre,self%xml_pre_def)
  self%mf_solver%pre%A=>self%nlfun%jac_op 
ELSE
  self%mf_solver%pre=>self%pre
END IF
!------------------------------------------------------------------------------
! Setup Newton Solver
!------------------------------------------------------------------------------
self%nksolver%A=>self%nlfun
self%nksolver%J_inv=>self%mf_solver
self%nksolver%its=20
self%nksolver%atol=self%nl_tol
self%nksolver%rtol=1.d-20 ! \
self%nksolver%backtrack=.TRUE.
self%nksolver%J_update=>gs_mfnk_update
self%nksolver%up_freq=1
end subroutine setup

subroutine run_simulation(self)
class(oft_gs_xmhd_sim), intent(inout) :: self !< NL operator object
character(LEN=TDIFF_RST_LEN) :: rst_char
real(r8), pointer :: plot_vals(:), tmp_arr(:), plot_vec(:,:)
real(r8) :: elapsed_time
integer(i4) :: i,j, io_stat, rst_tmp
type(oft_timer) :: mytimer
CLASS(oft_native_matrix), POINTER :: A_native
class(oft_vector), pointer :: tmp_vec
type(oft_lag_bginterp) :: grad_psi
class(oft_vector), pointer :: ux,uy,uz,v_lag, xtmp
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
TYPE(poss_scalar_bfield) :: field_init
self%t=0.d0
CALL oft_blagrange_2%vec_create(tmp_vec)
CALL self%tmp%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(self%u, self%t, self%dt, 'gs_xmhd_'//rst_char//'.rst', 'U')

NULLIFY(plot_vals)
ALLOCATE(plot_vec(3,tmp_vec%n))
call oft_blagrange_2%vec_create(grad_psi%u)
call oft_blagrange_2%vec_create(ux)
call oft_blagrange_2%vec_create(uy)
call oft_blagrange_2%vec_create(uz)
call oft_blagrange_2%vec_create(xtmp)
call oft_blagrange_2%vec_create(v_lag)
NULLIFY(lmop)
call oft_blag_getmop(oft_blagrange_2,lmop,'none')
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre)

CALL self%xdmf_plot%add_timestep(self%t)
CALL self%xdmf_plot_1%add_timestep(self%t)
CALL self%u%get_local(plot_vals,1)
CALL mesh_1%save_vertex_scalar(plot_vals,self%xdmf_plot_1,'p')
NULLIFY(plot_vals)
CALL self%u%get_local(plot_vals,2)
plot_vec(1,:)=plot_vals
CALL self%u%get_local(plot_vals,3)
plot_vec(2,:)=plot_vals
CALL self%u%get_local(plot_vals,4)
plot_vec(3,:)=plot_vals
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'V')
CALL grad_psi%u%restore_local(plot_vec(1,:))
CALL grad_psi%setup(oft_blagrange_2)
CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,ux)
CALL ux%add(0.d0,1.d0,v_lag)
CALL grad_psi%u%restore_local(plot_vec(3,:))
CALL grad_psi%setup(oft_blagrange_2)
CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uy)
CALL uy%add(0.d0,1.d0,v_lag)
!CALL uy%add(1.d0, 1.d0, ux) ! dur/dr + duz/dz
field_init%mesh=>mesh
field_init%func=>r_init
CALL oft_blag_project(oft_blagrange_2,field_init,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uz)
CALL uz%add(0.d0,1.d0,v_lag) !r
CALL uz%get_local(plot_vals)
CALL uz%restore_local(plot_vec(1,:)/plot_vals)
!CALL uy%add(1.d0, 1.d0, uz) ! dur/dr + duz/dz + ur/r
CALL uy%get_local(plot_vals)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'div_v')

CALL self%u%get_local(plot_vals,6)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
CALL grad_psi%u%restore_local(plot_vals)
CALL grad_psi%setup(oft_blagrange_2)
CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,ux)
CALL ux%add(0.d0,1.d0,v_lag)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uy)
CALL uy%add(0.d0,1.d0,v_lag)
CALL uy%get_local(plot_vals)
plot_vec(1,:)=-plot_vals
CALL self%u%get_local(plot_vals,5)
plot_vec(2,:)=plot_vals
CALL ux%get_local(plot_vals)
plot_vec(3,:)=plot_vals
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'B')

!plot poloidal current
CALL self%u%get_local(plot_vals,5)
CALL grad_psi%u%restore_local(plot_vals)
CALL grad_psi%setup(oft_blagrange_2)
CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,ux)
CALL ux%add(0.d0,1.d0,v_lag)
CALL ux%get_local(plot_vals)
plot_vec(3,:)=plot_vals
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uy)
CALL uy%add(0.d0,1.d0,v_lag)
CALL uy%get_local(plot_vals)
plot_vec(1,:)=-plot_vals
plot_vec(2,:) = 0.d0
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'J')

DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick()
    ! Update time-advance operator
    CALL build_approx_jacobian(self,self%nlfun%jac_op, self%u)
    CALL self%pre%update(.TRUE.)
    CALL self%mf_solver%pre%update(.TRUE.)
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
    IF(ASSOCIATED(self%eq)) THEN
        CALL self%eq%zerob_bc%apply(tmp_vec)
    END IF
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
            CALL self%mf_solver%pre%update(.TRUE.)
            CALL apply_rhs(self%nlfun,self%u,self%rhs)
            CALL self%rhs%get_local(tmp_arr,6)
            CALL tmp_vec%restore_local(tmp_arr)
            IF(ASSOCIATED(self%eq)) THEN
                CALL self%eq%zerob_bc%apply(tmp_vec)
            END IF
            CALL tmp_vec%get_local(tmp_arr)
            CALL self%rhs%restore_local(tmp_arr,6)
            CYCLE
        ELSE
            EXIT
        END IF
    END DO
    write(*,*) 'timestep = 1: ', i
    write(*,*) 'EXIT CODE', j
    self%t=self%t+self%nlfun%dt
    self%dt=self%nlfun%dt
    IF(ASSOCIATED(self%eq)) THEN
        self%eq%alam=self%nlfun%f_scale
        self%eq%pnorm=self%nlfun%p_scale
    END IF
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
        CALL self%xdmf_plot_1%add_timestep(self%t)
        NULLIFY(plot_vals)
        CALL self%u%get_local(plot_vals,1)
        CALL mesh_1%save_vertex_scalar(plot_vals,self%xdmf_plot_1,'p')
        NULLIFY(plot_vals)
        CALL self%u%get_local(plot_vals,2)
        plot_vec(1,:)=plot_vals
        CALL self%u%get_local(plot_vals,3)
        plot_vec(2,:)=plot_vals
        CALL self%u%get_local(plot_vals,4)
        plot_vec(3,:)=plot_vals
        CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'V')
        CALL grad_psi%u%restore_local(plot_vec(1,:))
        CALL grad_psi%setup(oft_blagrange_2)
        CALL oft_blag_vproject(oft_blagrange_2,grad_psi,xtmp,uy,uz)
        CALL v_lag%set(0.d0)
        CALL lminv%apply(v_lag,xtmp)
        CALL xtmp%add(0.d0,1.d0,v_lag)
        CALL grad_psi%u%restore_local(plot_vec(2,:))
        CALL grad_psi%setup(oft_blagrange_2)
        CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
        CALL v_lag%set(0.d0)
        CALL lminv%apply(v_lag,uy)
        CALL uy%add(0.d0,1.d0,v_lag)
        CALL uy%add(1.d0,1.d0,xtmp)
        field_init%mesh=>mesh
        field_init%func=>r_init
        CALL oft_blag_project(oft_blagrange_2,field_init,uz)
        CALL v_lag%set(0.d0)
        CALL lminv%apply(v_lag,uz)
        CALL uz%add(0.d0,1.d0,v_lag) !r
        CALL uz%get_local(plot_vals)
        CALL uz%restore_local(plot_vec(1,:)/plot_vals)
        CALL uy%add(1.d0, 1.d0, uz) ! dur/dr + duz/dz + ur/r
        CALL uy%get_local(plot_vals)
        CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'div_v')
        CALL self%u%get_local(plot_vals,6)
        CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
        CALL grad_psi%u%restore_local(plot_vals)
        CALL grad_psi%setup(oft_blagrange_2)
        CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
        CALL v_lag%set(0.d0)
        CALL lminv%apply(v_lag,ux)
        CALL ux%add(0.d0,1.d0,v_lag)
        CALL v_lag%set(0.d0)
        CALL lminv%apply(v_lag,uy)
        CALL uy%add(0.d0,1.d0,v_lag)
        CALL uy%get_local(plot_vals)
        plot_vec(1,:)=-plot_vals
        CALL self%u%get_local(plot_vals,5)
        plot_vec(2,:)=plot_vals
        CALL ux%get_local(plot_vals)
        plot_vec(3,:)=plot_vals
        CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'B')

        !plot poloidal current
        CALL self%u%get_local(plot_vals,5)
        CALL grad_psi%u%restore_local(plot_vals)
        CALL grad_psi%setup(oft_blagrange_2)
        CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
        CALL v_lag%set(0.d0)
        CALL lminv%apply(v_lag,ux)
        CALL ux%add(0.d0,1.d0,v_lag)
        CALL ux%get_local(plot_vals)
        plot_vec(3,:)=plot_vals
        CALL v_lag%set(0.d0)
        CALL lminv%apply(v_lag,uy)
        CALL uy%add(0.d0,1.d0,v_lag)
        CALL uy%get_local(plot_vals)
        plot_vec(1,:)=-plot_vals
        plot_vec(2,:) = 0.d0
        CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'J')

    END IF 
  ! self%eq%Ip_ratio_target = self%eq%Ip_ratio_target*1.02
END DO
DEALLOCATE(self%p_bc, self%velx_bc, self%vely_bc, self%velz_bc, self%by_bc, self%plasma_bc, self%psi_bc)
DEALLOCATE(self%nlfun%eta_t, self%nlfun%eta_p, self%nlfun%curr, self%nlfun%region_flag)
DEALLOCATE(self%fe_rep%fields, self%fe_rep%field_tags, self%jacobian_block_mask) 
DEALLOCATE(plot_vec)
CALL self%mf_solver%delete()
CALL self%nksolver%delete()
end subroutine run_simulation

subroutine add_timestep(self, dt, t)
class(oft_gs_xmhd_sim), intent(inout) :: self !<simulation object
real(r8), intent(in) :: dt
real(r8), intent(in), optional :: t
character(LEN=TDIFF_RST_LEN) :: rst_char
real(r8), pointer :: plot_vals(:), tmp_arr(:), plot_vec(:,:)
real(r8) :: elapsed_time
integer(i4) :: i,j, io_stat, rst_tmp
type(oft_timer) :: mytimer
CLASS(oft_native_matrix), POINTER :: A_native
class(oft_vector), pointer :: tmp_vec
type(oft_lag_bginterp) :: grad_psi
class(oft_vector), pointer :: ux,uy,uz,v_lag, xtmp
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
TYPE(poss_scalar_bfield) :: field_init
CALL oft_blagrange_2%vec_create(tmp_vec)
CALL self%tmp%add(0.d0,1.d0,self%u)
self%dt = dt
self%nlfun%dt = self%dt
self%nlfun%tflux_source = self%tflux_source
self%nlfun%region_flag = self%region_flag
self%nlfun%lim_vac_int = self%lim_vac_int
self%nlfun%j_source_scale = self%j_source_scale
self%nlfun%j_source => self%j_source
CALL update_bcs(self)
oft_env%pm = self%pm
self%mf_solver%pm=oft_env%pm
write(*,*) self%t
IF (PRESENT(t)) THEN
  self%t = t 
END IF
IF(oft_env%head_proc)CALL mytimer%tick()
! Update time-advance operator
CALL build_approx_jacobian(self,self%nlfun%jac_op, self%u)
CALL build_vac_jacobian(self,self%nlfun%vac_op)
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
    CALL apply_rhs(self%nlfun,self%u,self%rhs)
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
self%eq%alam=self%nlfun%f_scale
self%eq%pnorm=self%nlfun%p_scale
IF(oft_env%head_proc) CALL mytimer%tick
!---Create restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)self%rst_base
READ(rst_char,104,IOSTAT=io_stat)rst_tmp
IF((io_stat/=0).OR.(rst_tmp/=self%rst_base))CALL oft_abort("Step count exceeds format width", "run_simulation", __FILE__)
CALL self%rst_save(self%u, self%t, self%dt, 'gs_xmhd_'//rst_char//'.rst', 'U')
IF(oft_env%head_proc)THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
END IF
!---
CALL self%xdmf_plot%add_timestep(self%t)
CALL self%xdmf_plot_1%add_timestep(self%t)
NULLIFY(plot_vals)
ALLOCATE(plot_vec(3,tmp_vec%n))
call oft_blagrange_2%vec_create(grad_psi%u)
call oft_blagrange_2%vec_create(ux)
call oft_blagrange_2%vec_create(uy)
call oft_blagrange_2%vec_create(uz)
call oft_blagrange_2%vec_create(xtmp)
call oft_blagrange_2%vec_create(v_lag)
NULLIFY(lmop)
call oft_blag_getmop(oft_blagrange_2,lmop,'none')
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre)
CALL self%u%get_local(plot_vals,1)
CALL mesh_1%save_vertex_scalar(plot_vals,self%xdmf_plot_1,'p')
NULLIFY(plot_vals)
CALL self%u%get_local(plot_vals,2)
plot_vec(1,:)=plot_vals
CALL self%u%get_local(plot_vals,3)
plot_vec(2,:)=plot_vals
CALL self%u%get_local(plot_vals,4)
plot_vec(3,:)=plot_vals
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'V')
CALL grad_psi%u%restore_local(plot_vec(1,:))
CALL grad_psi%setup(oft_blagrange_2)
CALL oft_blag_vproject(oft_blagrange_2,grad_psi,xtmp,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,xtmp)
CALL xtmp%add(0.d0,1.d0,v_lag)
CALL grad_psi%u%restore_local(plot_vec(2,:))
CALL grad_psi%setup(oft_blagrange_2)
CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uy)
CALL uy%add(0.d0,1.d0,v_lag)
CALL uy%add(1.d0,1.d0,xtmp)
field_init%mesh=>mesh
field_init%func=>r_init
CALL oft_blag_project(oft_blagrange_2,field_init,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uz)
CALL uz%add(0.d0,1.d0,v_lag) !r
CALL uz%get_local(plot_vals)
CALL uz%restore_local(plot_vec(1,:)/plot_vals)
CALL uy%add(1.d0, 1.d0, uz) ! dur/dr + duz/dz + ur/r
CALL uy%get_local(plot_vals)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'div_v')
CALL self%u%get_local(plot_vals,6)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
CALL grad_psi%u%restore_local(plot_vals)
CALL grad_psi%setup(oft_blagrange_2)
CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,ux)
CALL ux%add(0.d0,1.d0,v_lag)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uy)
CALL uy%add(0.d0,1.d0,v_lag)
CALL uy%get_local(plot_vals)
plot_vec(1,:)=-plot_vals
CALL self%u%get_local(plot_vals,5)
plot_vec(2,:)=plot_vals
CALL ux%get_local(plot_vals)
plot_vec(3,:)=plot_vals
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'B')
!plot poloidal current
CALL self%u%get_local(plot_vals,5)
CALL grad_psi%u%restore_local(plot_vals)
CALL grad_psi%setup(oft_blagrange_2)
CALL oft_blag_vproject(oft_blagrange_2,grad_psi,ux,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,ux)
CALL ux%add(0.d0,1.d0,v_lag)
CALL ux%get_local(plot_vals)
plot_vec(3,:)=plot_vals
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uy)
CALL uy%add(0.d0,1.d0,v_lag)
CALL uy%get_local(plot_vals)
plot_vec(1,:)=-plot_vals
plot_vec(2,:) = 0.d0
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'J')
! self%t = self%t + self%dt
! write(*,*) self%nlfun%f_scale
! write(*,*) self%eq%diverted
! write(*,*) self%eq%plasma_bounds
! DEALLOCATE(self%p_bc, self%velx_bc, self%vely_bc, self%velz_bc, self%by_bc, self%plasma_bc)
! DEALLOCATE(self%nlfun%eta_t, self%nlfun%eta_p, self%nlfun%curr, self%nlfun%region_flag)
! DEALLOCATE(self%fe_rep%fields, self%fe_rep%field_tags, self%jacobian_block_mask) 
DEALLOCATE(plot_vec)
! CALL self%mf_solver%delete()
! CALL self%nksolver%delete()
end subroutine add_timestep

SUBROUTINE nlfun_apply(self, a, b)
class(gs_xmhd_nlfun), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: ptmp !temporary storage vector
type(oft_quad_type), pointer :: quad
type(oft_quad_type) :: quad_1d
LOGICAL :: curved
INTEGER(i4) :: i,m,jr, k,l,j, v !indexing variables for loops
INTEGER(i4) :: np_lim, cell, ed
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs_1, cell_dofs_2, cell_b_dofs
INTEGER(i4), allocatable :: elist(:,:)
REAL(r8) :: eta_t_loc, eta_p_loc, curr_loc
REAL(r8) :: nu, rho, B_0(3) ! physics parameters
REAL(r8) :: p_source, f_source, diag(2) ! Used for scaling P' and FF'
REAL(r8) :: p, dp(3), vel(3), by, dby(3), psi, dpsi(3), dvel(3,3), div_vel, btmp(3), F0_res !reconstructed variables
REAL(r8) :: coords(3), jac_det, jac_mat(3,4), tmp1(3), pts(2,2), dl(2), dn(3), dl_mag, f(3) ! For integration
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals_1,basis_vals_2,basis_vals, p_weights_loc, by_weights_loc, psi_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads_1, basis_grads_2, basis_grads, vel_weights_loc, res_loc
REAL(r8), POINTER, DIMENSION(:) :: p_weights, psi_weights, by_weights
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
REAL(r8), POINTER, DIMENSION(:) :: p_res, velx_res, vely_res, velz_res, by_res, psi_res, pres_vals,alam_vals, vtmp, by_res_plasma
REAL(r8) :: signed_area, x1, y1, x2, y2

quad=>oft_blagrange_2%quad

NULLIFY( p_weights, p_res, vel_weights, velx_res, vely_res, velz_res, by_weights, by_res, &
 psi_weights, psi_res, pres_vals, alam_vals, vtmp)

!---Get weights from solution vector
ALLOCATE(vel_weights(3,oft_blagrange_2%ne))
CALL a%get_local(p_weights, 1)
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(by_weights, 5)
CALL a%get_local(psi_weights, 6)

!--Update equilibrium with current values of psi
IF(ANY(self%region_flag ==5)) THEN
  CALL self%eq%psi%restore_local(psi_weights)
  CALL gs_update_bounds(self%eq,track_opoint=.TRUE.)
  self%eq%I%plasma_bounds=self%eq%plasma_bounds
  self%eq%P%plasma_bounds=self%eq%plasma_bounds
END IF
!write(*,*) self%eq%lim_point
!--Initialize residuals with zeros
CALL b%set(0.d0)
CALL b%get_local(p_res, 1)
CALL b%get_local(velx_res, 2)
CALL b%get_local(vely_res, 3)
CALL b%get_local(velz_res, 4)
CALL b%get_local(by_res, 5)
CALL b%get_local(psi_res, 6)
CALL b%get_local(alam_vals, 6)
CALL b%get_local(pres_vals, 6)

!--Set local physics parameters
nu = self%nu
rho = self%rho
B_0 = self%B_0
diag = 0.d0
! Declare variables private for OMP
!$omp parallel num_threads(1) private(i, m,j,jr,k,l,curved,coords,cell_dofs_1, cell_dofs_2,basis_vals_1,&
!$omp basis_vals_2,basis_grads_1, basis_grads_2, p_weights_loc, vel_weights_loc,&  
!$omp  psi_weights_loc, by_weights_loc,res_loc,jac_mat, jac_det, &
!$omp p, dp, vel, dvel, div_vel, psi, dpsi, by, dby, eta_t_loc, eta_p_loc, v, btmp, &
!$omp p_source, f_source)reduction(+:diag)

! Allocate local variables
ALLOCATE(basis_vals_1(oft_blagrange_1%nce),basis_grads_1(3,oft_blagrange_1%nce))
ALLOCATE(basis_vals_2(oft_blagrange_2%nce),basis_grads_2(3,oft_blagrange_2%nce))
ALLOCATE(p_weights_loc(oft_blagrange_1%nce))
ALLOCATE(vel_weights_loc(3, oft_blagrange_2%nce))
ALLOCATE(by_weights_loc(oft_blagrange_2%nce))
ALLOCATE(psi_weights_loc(oft_blagrange_2%nce))
ALLOCATE(cell_dofs_1(oft_blagrange_1%nce), cell_dofs_2(oft_blagrange_2%nce))
ALLOCATE(res_loc(oft_blagrange_2%nce,8)) ! 7 fields + entries to store pres_vals and alam_vals
!$omp do schedule(static)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange_1%ncdofs(i,cell_dofs_1) ! Get global index of local DOFs
  call oft_blagrange_2%ncdofs(i,cell_dofs_2) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function

  ! Set local weights
  p_weights_loc = p_weights(cell_dofs_1)
  vel_weights_loc = vel_weights(:, cell_dofs_2)
  by_weights_loc = by_weights(cell_dofs_2)
  psi_weights_loc = psi_weights(cell_dofs_2)
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange_1%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange_1,i,jr,quad%pts(:,m),basis_vals_1(jr))
      CALL oft_blag_geval(oft_blagrange_1,i,jr,quad%pts(:,m),basis_grads_1(:,jr),jac_mat)
    END DO
    DO jr=1,oft_blagrange_2%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange_2,i,jr,quad%pts(:,m),basis_vals_2(jr))
      CALL oft_blag_geval(oft_blagrange_2,i,jr,quad%pts(:,m),basis_grads_2(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    ! Switch basis grads from 2D to 3D
    basis_grads_1(3, :) = basis_grads_1(2,:)
    basis_grads_1(2,:) = 0.d0
    basis_grads_2(3, :) = basis_grads_2(2,:)
    basis_grads_2(2,:) = 0.d0
    !---Reconstruct values of solution fields
    p = 0.d0; dp=0.d0
    vel = 0.d0; dvel = 0.d0; div_vel = 0.d0
    by = 0.d0; dby = 0.d0
    psi = 0.d0; dpsi=0.d0

    DO jr=1, oft_blagrange_1%nce
      p = p + p_weights_loc(jr)*basis_vals_1(jr)
      dp = dp + p_weights_loc(jr)*basis_grads_1(:,jr)
    END DO
    DO jr=1,oft_blagrange_2%nce
      vel = vel + vel_weights_loc(:, jr)*basis_vals_2(jr)
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads_2(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads_2(3, jr)
      by = by + by_weights_loc(jr)*basis_vals_2(jr)
      dby = dby + by_weights_loc(jr)*basis_grads_2(:,jr)
      psi = psi + psi_weights_loc(jr)*basis_vals_2(jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads_2(:,jr)
    END DO

    eta_t_loc = self%eta_t(mesh%reg(i))
    eta_p_loc = self%eta_p(mesh%reg(i))
    IF(ALLOCATED(self%eta_node)) THEN
      eta_t_loc = 0.d0 
      DO v = 1, SIZE(mesh%lc, 1)
        eta_t_loc = eta_t_loc + self%eta_node(mesh%lc(v, i)) * quad%pts(v, m)  
      END DO
      eta_p_loc = eta_t_loc
    END IF
    div_vel = dvel(1,1) +vel(1)/(coords(1)+gs_epsilon) + dvel(3,3)

    btmp = cross_product(dpsi/(coords(1)+gs_epsilon), [0.d0,1.d0,0.d0]) + [0.d0,1.d0,0.d0]*by/(coords(1)+gs_epsilon) + B_0
    !If WE ARE IN AN MHD REGION
    IF(self%region_flag(mesh%reg(i)) == 1) THEN
    ! First, compute pressure residual (incompressibility equation)
      DO jr=1,oft_blagrange_1%nce
        !PRESSURE
        res_loc(jr,1) = res_loc(jr,1) &
        + basis_vals_1(jr)*div_vel*jac_det*quad%wts(m)*coords(1)
      END DO
    ! Now compute the other xMHD residuals
      DO jr=1,oft_blagrange_2%nce
        !VELOCITY
        res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          + basis_vals_2(jr)*vel*jac_det*quad%wts(m)*coords(1) &
          + self%dt*DOT_PRODUCT(btmp,basis_grads_2(:,jr))*btmp*jac_det*quad%wts(m)*coords(1)/(mu0*rho) & !MAGNETIC FORCES
          - self%dt*DOT_PRODUCT(btmp,btmp)*basis_grads_2(:,jr)*jac_det*quad%wts(m)*coords(1)/(2*mu0*rho) & !MAGNETIC FORCES
          - self%dt*p*basis_grads_2(:,jr)*jac_det*quad%wts(m)*coords(1)/rho !PRESSURE FORCE
        DO k=1,3
          res_loc(jr,k+1) = res_loc(jr, k+1) &
            + basis_vals_2(jr)*self%dt*DOT_PRODUCT(vel,dvel(k,:))*jac_det*quad%wts(m)*coords(1) &
            + nu*self%dt*DOT_PRODUCT(basis_grads_2(:,jr),dvel(k,:))*jac_det*quad%wts(m)*coords(1)/rho
        END DO 
        res_loc(jr,2) = res_loc(jr,2) &
          + self%dt*basis_vals_2(jr)*(btmp(2)**2-btmp(1)**2-btmp(3)**2)*jac_det*quad%wts(m)/(2.d0*mu0*rho) &
          + self%dt*basis_vals_2(jr)*nu*vel(1)*jac_det*quad%wts(m)/(rho*(coords(1)+gs_epsilon))  &
          - self%dt*basis_vals_2(jr)*vel(2)**2*jac_det*quad%wts(m) &
          - self%dt*basis_vals_2(jr)*p*jac_det*quad%wts(m)/rho
      
        res_loc(jr,3) = res_loc(jr,3) &
          - self%dt*basis_vals_2(jr)*btmp(1)*btmp(2)*jac_det*quad%wts(m)/(mu0*rho) & !!nonzero 
          + self%dt*basis_vals_2(jr)*nu*vel(2)*jac_det*quad%wts(m)/(rho*(coords(1)+gs_epsilon)) &
          + self%dt*basis_vals_2(jr)*vel(1)*vel(2)*jac_det*quad%wts(m)
        
        ! !BY (F) (ONLY INCLUDE XMHD SPECIFIC TERMS)
        tmp1 = cross_product(dpsi,dvel(2, :))
        res_loc(jr,5) = res_loc(jr,5) &
          - self%dt*basis_vals_2(jr)*tmp1(2)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + self%dt*basis_vals_2(jr)*dvel(1,1)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + self%dt*basis_vals_2(jr)*dvel(3,3)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + self%dt*basis_vals_2(jr)*DOT_PRODUCT(vel, dby)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + self%dt*basis_vals_2(jr)*vel(1)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)**2
        !PSI (Here, I only include the terms that are not included in vac_op)
        res_loc(jr,6) = res_loc(jr,6) &
          + self%dt*basis_vals_2(jr)*DOT_PRODUCT(vel, dpsi)*jac_det*quad%wts(m)/(eta_t_loc*(coords(1)+gs_epsilon)) &
          + self%dt*basis_vals_2(jr)*tmp1(2)*jac_det*quad%wts(m)/(eta_t_loc*(coords(1)+gs_epsilon))
      END DO
    END IF

    ! ADD F DIFFUSION RESIDUAL EVERYWHERE (will overwrite in plasma region)
    DO jr=1,oft_blagrange_2%nce
      res_loc(jr,5) = res_loc(jr,5) &
        + basis_vals_2(jr)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
        + self%dt*eta_p_loc*DOT_PRODUCT(basis_grads_2(:,jr), dby)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
    END DO

    ! IF WE ARE IN THE PLASMA
    IF(self%region_flag(mesh%reg(i)) == 5) THEN
      IF (gs_test_bounds(self%eq,coords) .AND. psi >self%eq%plasma_bounds(1)) THEN !check that we are in the plasma
          p_source = self%p_scale*self%eq%P%Fp(psi)*coords(1) 
          f_source = self%f_scale*0.5d0* self%eq%I%fp(psi)/ (coords(1) + gs_epsilon)
          diag=diag+[f_source,p_source]*jac_det*quad%wts(m)
          DO jr=1,oft_blagrange_2%nce
              res_loc(jr,7) = res_loc(jr,7) &
              - self%dt * basis_vals_2(jr) * p_source * jac_det*quad%wts(m)
              res_loc(jr,8) = res_loc(jr,8) &
              - self%dt * basis_vals_2(jr) * f_source * jac_det*quad%wts(m)
          END DO
      END IF
    END IF
  END DO

    !---Add local values to full vector
  DO jr=1, oft_blagrange_1%nce
   !$omp atomic
    p_res(cell_dofs_1(jr)) = p_res(cell_dofs_1(jr)) + res_loc(jr,1)
  END DO
  DO jr=1,oft_blagrange_2%nce
    !$omp atomic
    velx_res(cell_dofs_2(jr)) = velx_res(cell_dofs_2(jr)) + res_loc(jr,2)
    !$omp atomic
    vely_res(cell_dofs_2(jr)) = vely_res(cell_dofs_2(jr)) + res_loc(jr,3)
    !$omp atomic
    velz_res(cell_dofs_2(jr)) = velz_res(cell_dofs_2(jr)) + res_loc(jr,4)
    !$omp atomic
    by_res(cell_dofs_2(jr)) = by_res(cell_dofs_2(jr)) + res_loc(jr,5)
    !$omp atomic
    psi_res(cell_dofs_2(jr)) = psi_res(cell_dofs_2(jr)) + res_loc(jr,6)
    !$omp atomic
    pres_vals(cell_dofs_2(jr)) = pres_vals(cell_dofs_2(jr)) + res_loc(jr,7)
    !$omp atomic
    alam_vals(cell_dofs_2(jr)) = alam_vals(cell_dofs_2(jr)) + res_loc(jr,8)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals_1, basis_vals_2, basis_grads_1, basis_grads_2, p_weights_loc, vel_weights_loc,  &
 psi_weights_loc,by_weights_loc, cell_dofs_1, cell_dofs_2,res_loc)
!$omp end parallel

 DO i=1,oft_blagrange_2%nbe
    alam_vals(oft_blagrange_2%lbe(i))=0.d0
    pres_vals(oft_blagrange_2%lbe(i))=0.d0
END DO
! RESCALE EQUATIONS --> add some conditions to this?
IF(ANY(self%region_flag ==5)) THEN
  f_source = self%eq%Itor_target/diag(1)/(1.d0+1.d0/self%eq%Ip_ratio_target)
  p_source = self%eq%Itor_target/diag(2)/(self%eq%Ip_ratio_target+1.d0)
  psi_res= psi_res + pres_vals*p_source+alam_vals*f_source
  self%f_scale=f_source*self%f_scale
  self%p_scale=p_source*self%p_scale
  diag(1)=diag(1)*f_source
  diag(2)=diag(2)*p_source
END IF
!Loop to compute F0 residual
F0_res = 0.d0
IF (any(self%region_flag == 5) .OR. any(self%region_flag == 6)) THEN
  ! Declare variables private for OMP
  !$omp parallel num_threads(1) private(i, m,jr,curved,coords,basis_vals_2, psi_weights_loc, cell_dofs_2,&
  !$omp  by_weights_loc, jac_mat, jac_det, psi) reduction(+:F0_res)
  ! Allocate local variables
  ALLOCATE(basis_vals_2(oft_blagrange_2%nce))
  ALLOCATE(psi_weights_loc(oft_blagrange_2%nce))
  ALLOCATE(cell_dofs_2(oft_blagrange_2%nce))
  ALLOCATE(by_weights_loc(oft_blagrange_2%nce))
  !$omp do schedule(static)
  DO i=1,mesh%nc
    curved=cell_is_curved(mesh,i) ! Straight cell test
    call oft_blagrange_2%ncdofs(i,cell_dofs_2) ! Get global index of local DOFs
    by_weights_loc = by_weights(cell_dofs_2)
    ! Set local weights
    psi_weights_loc = psi_weights(cell_dofs_2)
    !---------------------------------------------------------------------------
    ! Quadrature Loop
    !---------------------------------------------------------------------------
    DO m=1,quad%np
      if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
      !---Evaluate value and gradients of basis functions at current point
      DO jr=1,oft_blagrange_2%nce ! Loop over degrees of freedom
        CALL oft_blag_eval(oft_blagrange_2,i,jr,quad%pts(:,m),basis_vals_2(jr))
      END DO
      !--Extract spatial coordinates at current point
      coords = mesh%log2phys(i,quad%pts(:,m))
      !---Reconstruct values of solution fields
      psi = 0.d0
      DO jr=1,oft_blagrange_2%nce
        psi = psi + psi_weights_loc(jr)*basis_vals_2(jr)
      END DO
      ! Add contributions to F0 residual from inside limiter
      IF(self%region_flag(mesh%reg(i)) == 5) THEN
        ! IF WE ARE IN THE PLASMA
        IF (gs_test_bounds(self%eq,coords) .AND. psi >self%eq%plasma_bounds(1)) THEN !check that we are in the plasma
          F0_res = F0_res + (SQRT(self%f_scale*self%eq%I%f(psi) + by_weights(self%lim_ind)**2))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        ELSE
          F0_res = F0_res + by_weights(self%lim_ind)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        END IF
      END IF
    END DO
  END DO
  !---Cleanup thread-local storage
  DEALLOCATE(basis_vals_2,   psi_weights_loc,cell_dofs_2, by_weights_loc)
  !$omp end parallel
  !BOUNDARY INTEGRAL CONTRIBUTION TO F0
  !For each edge, find cell and corresponding local edge index
  signed_area = 0.0d0 ! Initialize before the loop
  ALLOCATE(cell_b_dofs(oft_blagrange_2%nce))
  ALLOCATE(basis_vals(oft_blagrange_2%nce),basis_grads(3,oft_blagrange_2%nce))
  ALLOCATE(by_weights_loc(oft_blagrange_2%nce))
  ALLOCATE(elist(2,SIZE(self%eq%lim_con)))
  DO i = 1, SIZE(self%eq%lim_con)-1
    ! --- NEW: Accumulate Shoelace Area ---
    x1 = mesh%r(1, self%eq%lim_con(i))
    y1 = mesh%r(2, self%eq%lim_con(i))
    x2 = mesh%r(1, self%eq%lim_con(i+1))
    y2 = mesh%r(2, self%eq%lim_con(i+1))
    signed_area = signed_area + (x1 * y2 - x2 * y1)

    j=ABS(mesh_local_findedge(mesh,[self%eq%lim_con(i),self%eq%lim_con(i+1)]))
    IF(self%region_flag(mesh%reg(mesh%lec(mesh%kec(j)))) /= 5 .AND. self%region_flag(mesh%reg(mesh%lec(mesh%kec(j)))) /= 6) THEN
      elist(2,i)=mesh%lec(mesh%kec(j))
    ELSE
      elist(2,i)=mesh%lec(mesh%kec(j)+1)
    END IF
    DO m=1,3
      IF(j==ABS(mesh%lce(m,elist(2,i))))THEN
        elist(1,i)=m
        EXIT
      END IF
    END DO
  END DO
  i = SIZE(self%eq%lim_con)
  j=ABS(mesh_local_findedge(mesh,[self%eq%lim_con(SIZE(self%eq%lim_con)),self%eq%lim_con(1)]))
  IF(self%region_flag(mesh%reg(mesh%lec(mesh%kec(j)))) /= 5) THEN
    elist(2,i)=mesh%lec(mesh%kec(j))
  ELSE
    elist(2,i)=mesh%lec(mesh%kec(j)+1)
  END IF
  DO m=1,3
    IF(j==ABS(mesh%lce(m,elist(2,i))))THEN
      elist(1,i)=m
      EXIT
    END IF
  END DO
  !Setup 1D quadrature
  CALL set_quad_1d(quad_1d,oft_blagrange_2%order+2)
  !Begin integrating
  DO i = 1, SIZE(self%eq%lim_con)
    cell=elist(2,i)
    ed=elist(1,i)
    eta_p_loc = self%eta_p(mesh%reg(cell))
    pts(:,1)=mesh%r(1:2,mesh%lc(mesh%cell_ed(1,ed),cell))
    pts(:,2)=mesh%r(1:2,mesh%lc(mesh%cell_ed(2,ed),cell))
    dl=pts(:,1)-pts(:,2)
    IF(self%eq%lim_con(i)==mesh%lc(mesh%cell_ed(2,ed),cell))dl=-dl
    dl_mag=SQRT(SUM(dl**2))
    dn=[-dl(2),dl(1), 0.d0]
    CALL oft_blagrange_2%ncdofs(cell,cell_b_dofs)
    by_weights_loc = by_weights(cell_b_dofs)
    DO k=1,quad%np
      f = 0.d0
      f(mesh%cell_ed(1,ed))=quad%pts(1,k)
      f(mesh%cell_ed(2,ed))=1.d0 - quad%pts(1,k)
      coords=mesh%log2phys(cell,f)
      CALL mesh%jacobian(cell,f,jac_mat,jac_det)
      DO jr=1,oft_blagrange_2%nce
        CALL oft_blag_eval(oft_blagrange_2,cell,jr,f,basis_vals(jr))
        CALL oft_blag_geval(oft_blagrange_2,cell,jr,f,basis_grads(:,jr),jac_mat)
      END DO
      dby = 0.d0
      DO jr=1,oft_blagrange_2%nce
        dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
      END DO
      F0_res = F0_res  - SIGN(1.0d0, signed_area)*eta_p_loc * self%dt * DOT_PRODUCT(dby, dn)*quad%wts(k)/(coords(1)+gs_epsilon)
    END DO
  END DO
  DEALLOCATE(basis_vals, basis_grads ,cell_b_dofs, by_weights_loc)
  IF (.NOT. any(self%region_flag == 5)) THEN
    F0_res = F0_res + by_weights(self%lim_ind)*self%lim_vac_int
  END IF
END IF
! Apply BCs
CALL fem_dirichlet_vec(oft_blagrange_1,p_weights,p_res,self%p_bc)
CALL fem_dirichlet_vec(oft_blagrange_2,vel_weights(1, :),velx_res,self%velx_bc)
CALL fem_dirichlet_vec(oft_blagrange_2,vel_weights(2, :),vely_res,self%vely_bc)
CALL fem_dirichlet_vec(oft_blagrange_2,vel_weights(3, :),velz_res,self%velz_bc)
CALL fem_dirichlet_vec(oft_blagrange_2,by_weights,by_res,self%by_bc)
IF (.NOT.ASSOCIATED(self%eq)) THEN
  CALL fem_dirichlet_vec(oft_blagrange_2,psi_weights,psi_res,self%psi_bc)
END IF
! For F, want to overwrite residual with current weight - F0 for plasma nodes
IF (self%evolve_F) THEN
  ALLOCATE(by_res_plasma(oft_blagrange_2%ne))
  by_res_plasma = by_weights - by_weights(self%lim_ind)
  CALL fem_dirichlet_vec(oft_blagrange_2,by_res_plasma,by_res,self%plasma_bc)
  ! Overwrite F0 node with correct residual
  by_res(self%lim_ind) = F0_res
  DEALLOCATE(by_res_plasma)
END IF

!PUT IN OUTPUT VECTOR
CALL b%restore_local(p_res,1,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(velx_res,2,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(vely_res,3,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(velz_res,4,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(by_res,5,add=.TRUE., wait = .TRUE.)
CALL b%restore_local(psi_res,6,add=.TRUE.)
CALL b%new(ptmp)
CALL self%vac_op%apply(a,ptmp)
CALL b%add(1.d0,1.d0,ptmp)
CALL b%get_local(psi_res,6)
CALL ptmp%delete
DEALLOCATE(p_res, velx_res, vely_res, velz_res,psi_res, by_res,pres_vals, alam_vals)
DEALLOCATE(vel_weights, p_weights, psi_weights, by_weights)

END SUBROUTINE nlfun_apply

SUBROUTINE apply_rhs(self,a,b)
class(gs_xmhd_nlfun), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
LOGICAL :: curved
INTEGER(i4) :: i,m,jr, k,l, v, fu, ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs_1, cell_dofs_2
REAL(r8) :: eta_t_loc, curr_loc, source_tmp(1)
REAL(r8) ::  p, dp(3), vel(3), dvel(3,3),  psi, dpsi(3),by, dby(3), coords(3), jac_det, jac_mat(3,4), F0_res
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals_1, basis_vals_2, p_weights_loc, psi_weights_loc, by_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads_1, basis_grads_2, res_loc, vel_weights_loc
REAL(r8), POINTER, DIMENSION(:) :: p_weights,  psi_weights, by_weights 
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
REAL(r8), POINTER, DIMENSION(:) :: p_res, velx_res, vely_res, velz_res, psi_res, by_res, vtmp, by_res_plasma

quad=>oft_blagrange_2%quad
NULLIFY( p_weights, vel_weights, psi_weights, by_weights, &
p_res, velx_res, vely_res, velz_res, psi_res, by_res)
!---Get weights from solution vector
ALLOCATE(vel_weights(3,oft_blagrange_2%ne))
CALL a%get_local(p_weights, 1)
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(by_weights, 5)
CALL a%get_local(psi_weights, 6)

!--Initialize residuals with zeros
CALL b%set(0.d0)
CALL b%get_local(p_res, 1)
CALL b%get_local(velx_res, 2)
CALL b%get_local(vely_res, 3)
CALL b%get_local(velz_res, 4)
CALL b%get_local(by_res, 5)
CALL b%get_local(psi_res, 6)
!$omp parallel num_threads(1) private(i,m,jr,curved,coords,cell_dofs_1, cell_dofs_2,basis_vals_1, basis_vals_2,basis_grads_1, basis_grads_2, &
!$omp p_weights_loc, vel_weights_loc, psi_weights_loc, by_weights_loc,res_loc,jac_mat, &
!$omp jac_det, p, dp, vel, dvel, psi, dpsi, by, dby, eta_t_loc, curr_loc, source_tmp, v) reduction(+:F0_res)
!Allocate local arrays
ALLOCATE(basis_vals_1(oft_blagrange_1%nce),basis_grads_1(3,oft_blagrange_1%nce))
ALLOCATE(basis_vals_2(oft_blagrange_2%nce),basis_grads_2(3,oft_blagrange_2%nce))
ALLOCATE(p_weights_loc(oft_blagrange_1%nce))
ALLOCATE(vel_weights_loc(3, oft_blagrange_2%nce))
ALLOCATE(psi_weights_loc(oft_blagrange_2%nce))
ALLOCATE(by_weights_loc(oft_blagrange_2%nce))
ALLOCATE(cell_dofs_1(oft_blagrange_1%nce),cell_dofs_2(oft_blagrange_2%nce), res_loc(oft_blagrange_2%nce, 6))

!$omp do schedule(static)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange_1%ncdofs(i,cell_dofs_1) ! Get global index of local DOFs
  call oft_blagrange_2%ncdofs(i,cell_dofs_2) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function

  ! Set local weights
  p_weights_loc = p_weights(cell_dofs_1)
  vel_weights_loc = vel_weights(:, cell_dofs_2)
  by_weights_loc = by_weights(cell_dofs_2)
  psi_weights_loc = psi_weights(cell_dofs_2)
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange_1%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange_1,i,jr,quad%pts(:,m),basis_vals_1(jr))
      CALL oft_blag_geval(oft_blagrange_1,i,jr,quad%pts(:,m),basis_grads_1(:,jr),jac_mat)
    END DO

    DO jr=1,oft_blagrange_2%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange_2,i,jr,quad%pts(:,m),basis_vals_2(jr))
      CALL oft_blag_geval(oft_blagrange_2,i,jr,quad%pts(:,m),basis_grads_2(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    !---Reconstruct values of solution fields
    p = 0.d0; dp = 0.d0
    vel = 0.d0; dvel = 0.d0
    psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0

    ! switch from 2D to 3D gradients
    basis_grads_1(3, :) = basis_grads_1(2,:)
    basis_grads_1(2,:) = 0.d0
    basis_grads_2(3, :) = basis_grads_2(2,:)
    basis_grads_2(2,:) = 0.d0
    
    !Reconstruct values of solution fields at current point
    DO jr=1,oft_blagrange_1%nce
      p = p + p_weights_loc(jr)*basis_vals_1(jr)
      dp = dp + p_weights_loc(jr)*basis_grads_1(:,jr)
    END DO

    DO jr=1,oft_blagrange_2%nce
      vel = vel + vel_weights_loc(:, jr)*basis_vals_2(jr)
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads_2(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads_2(3, jr)
      psi = psi + psi_weights_loc(jr)*basis_vals_2(jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads_2(:,jr)
      by = by + by_weights_loc(jr)*basis_vals_2(jr)
      dby = dby + by_weights_loc(jr)*basis_grads_2(:,jr)
    END DO

    eta_t_loc = self%eta_t(mesh%reg(i))
    IF(ALLOCATED(self%eta_node)) THEN
      eta_t_loc = 0.d0 
      DO v = 1, SIZE(mesh%lc, 1)
        eta_t_loc = eta_t_loc + self%eta_node(mesh%lc(v, i)) * quad%pts(v, m)  
      END DO
    END IF
    curr_loc = self%curr(mesh%reg(i))


    !No RHS for incompressibility
    DO jr=1,oft_blagrange_2%nce
      IF(self%region_flag(mesh%reg(i)) == 1) THEN
        ! VELOCITY
        res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          + basis_vals_2(jr)*vel*jac_det*quad%wts(m)*coords(1)
      END IF  
      ! PSI
      IF (self%region_flag(mesh%reg(i))==1 .OR. self%region_flag(mesh%reg(i))==3 ) THEN
          res_loc(jr,6) = res_loc(jr,6) &
          + basis_vals_2(jr)*psi*jac_det*quad%wts(m)/(eta_t_loc*(coords(1)+gs_epsilon))
      END IF 
      IF (self%region_flag(mesh%reg(i))==4) THEN
          res_loc(jr,6) = res_loc(jr,6) &
          + basis_vals_2(jr)*self%dt*curr_loc*jac_det*quad%wts(m)
      END IF 
      IF (ASSOCIATED(self%j_source)) THEN
        CALL self%j_source%interp(i,quad%pts(:,m),jac_mat,source_tmp)
        res_loc(jr,6) = res_loc(jr,6) &
          + basis_vals_2(jr)*self%dt*self%j_source_scale*source_tmp(1)*jac_det*quad%wts(m)
      END IF
      ! BY (F)
      ! Add F diffusion residual everywhere
      res_loc(jr,5) = res_loc(jr,5) &
        + basis_vals_2(jr)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
    END DO
    ! Add contributions to F0 residual from inside limiter
    IF(self%region_flag(mesh%reg(i)) == 5) THEN
      ! IF WE ARE IN THE PLASMA
      IF (gs_test_bounds(self%eq,coords) .AND. psi >self%eq%plasma_bounds(1)) THEN !check that we are in the plasma
        F0_res = F0_res + (SQRT(self%eq%alam*self%eq%I%f(psi) + by_weights(self%lim_ind)**2))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
      ELSE
        F0_res = F0_res + by_weights(self%lim_ind)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
      END IF
    END IF
  END DO
    !---Add local values to full vector
  DO jr=1,oft_blagrange_2%nce
    !$omp atomic
    velx_res(cell_dofs_2(jr)) = velx_res(cell_dofs_2(jr)) + res_loc(jr,2)
    !$omp atomic
    vely_res(cell_dofs_2(jr)) = vely_res(cell_dofs_2(jr)) + res_loc(jr,3)
    !$omp atomic
    velz_res(cell_dofs_2(jr)) = velz_res(cell_dofs_2(jr)) + res_loc(jr,4)
    !$omp atomic
    by_res(cell_dofs_2(jr)) = by_res(cell_dofs_2(jr)) + res_loc(jr,5)
    !$omp atomic
    psi_res(cell_dofs_2(jr)) = psi_res(cell_dofs_2(jr)) + res_loc(jr,6)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals_1, basis_vals_2,basis_grads_1, basis_grads_2, cell_dofs_1, cell_dofs_2)
DEALLOCATE(p_weights_loc, vel_weights_loc, psi_weights_loc, by_weights_loc,res_loc)
!$omp end parallel


IF (.NOT. any(self%region_flag == 5)) THEN
  F0_res = F0_res + by_weights(self%lim_ind)*self%lim_vac_int + self%tflux_source
END IF

! SET BOUNDARY CONDITIONS
CALL fem_dirichlet_vec(oft_blagrange_1,p_weights,p_res,self%p_bc)
CALL fem_dirichlet_vec(oft_blagrange_2,vel_weights(1, :),velx_res,self%velx_bc)
CALL fem_dirichlet_vec(oft_blagrange_2,vel_weights(2, :),vely_res,self%vely_bc)
CALL fem_dirichlet_vec(oft_blagrange_2,vel_weights(3, :),velz_res,self%velz_bc)
CALL fem_dirichlet_vec(oft_blagrange_2,by_weights,by_res,self%by_bc)
IF (.NOT.ASSOCIATED(self%eq)) THEN
  CALL fem_dirichlet_vec(oft_blagrange_2,psi_weights,psi_res,self%psi_bc)
END IF
IF (self%evolve_F) THEN
  ALLOCATE(by_res_plasma(oft_blagrange_2%ne))
  ! For F, want RHS = 0 inside plasma
  by_res_plasma = 0.d0
  CALL fem_dirichlet_vec(oft_blagrange_2,by_res_plasma,by_res,self%plasma_bc)
  by_res(self%lim_ind) = F0_res
  DEALLOCATE(by_res_plasma)
END IF

DO i=1,oft_blagrange_2%nbe
    psi_res(oft_blagrange_2%lbe(i))=psi_weights(oft_blagrange_2%lbe(i))
END DO
CALL b%restore_local(p_res,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velx_res,2,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(vely_res,3,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velz_res,4,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(by_res,5,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(psi_res,6,add=.TRUE.)
DEALLOCATE(vel_weights, p_weights, psi_weights, by_weights)
END SUBROUTINE apply_rhs

SUBROUTINE gs_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a
CALL current_sim%mfmat%update(a)
END SUBROUTINE gs_mfnk_update

SUBROUTINE build_vac_jacobian(self, mat)
class (oft_gs_xmhd_sim), intent(inout) :: self
class (oft_matrix), pointer, intent(inout) :: mat
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals_1, basis_vals_2
REAL (r8) :: coords(3), eta_t_loc, jac_det, jac_mat(3,4)
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads_1, basis_grads_2
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
CLASS(oft_vector), POINTER :: oft_lag_vec
INTEGER(i4), ALLOCATABLE, DIMENSION(:), TARGET :: cell_dofs_1, cell_dofs_2
integer (i4) :: i, jr, jc, m, v
type(oft_quad_type), pointer :: quad
type(oft_1d_int), allocatable, dimension(:) :: iloc
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
LOGICAL :: curved
quad=>oft_blagrange_2%quad
CALL mat%zero
!--Setup thread locks
ALLOCATE(tlocks(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
  call omp_init_lock(tlocks(i))
END DO
!---
!$omp parallel num_threads(1) private(i,m,jr,jc,curved,coords, cell_dofs_1, cell_dofs_2,basis_vals_1, basis_vals_2, &
!$omp  basis_grads_1, basis_grads_2, jac_loc,jac_mat,jac_det,eta_t_loc, iloc, v)
ALLOCATE(basis_vals_1(oft_blagrange_1%nce),basis_grads_1(3,oft_blagrange_1%nce))
ALLOCATE(basis_vals_2(oft_blagrange_2%nce),basis_grads_2(3,oft_blagrange_2%nce))
ALLOCATE(cell_dofs_1(oft_blagrange_1%nce), cell_dofs_2(oft_blagrange_2%nce))
ALLOCATE(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
ALLOCATE(iloc(self%fe_rep%nfields))
iloc(1)%v=>cell_dofs_1
DO i=2,self%fe_rep%nfields
   iloc(i)%v=>cell_dofs_2
END DO
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do schedule(static)ordered
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange_1%ncdofs(i,cell_dofs_1) ! Get global index of local DOFs
  call oft_blagrange_2%ncdofs(i,cell_dofs_2) ! Get global index of local DOFs
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
!---------------------------------------------------------------------------
! Quadrature Loop
!---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange_1%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange_1,i,jr,quad%pts(:,m),basis_vals_1(jr))
      CALL oft_blag_geval(oft_blagrange_1,i,jr,quad%pts(:,m),basis_grads_1(:,jr),jac_mat)
    END DO
    DO jr=1,oft_blagrange_2%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange_2,i,jr,quad%pts(:,m),basis_vals_2(jr))
      CALL oft_blag_geval(oft_blagrange_2,i,jr,quad%pts(:,m),basis_grads_2(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    basis_grads_1(3, :) = basis_grads_1(2,:)
    basis_grads_1(2,:) = 0.d0
    basis_grads_2(3, :) = basis_grads_2(2,:)
    basis_grads_2(2,:) = 0.d0

    eta_t_loc = self%eta_t(mesh%reg(i))
    IF(ALLOCATED(self%eta_node)) THEN
      eta_t_loc = 0.d0 
      DO v = 1, SIZE(mesh%lc, 1)
        eta_t_loc = eta_t_loc + self%eta_node(mesh%lc(v, i)) * quad%pts(v, m)  
      END DO
    END IF
    !---Compute local matrix contributions
    DO jr=1,oft_blagrange_2%nce
      DO jc=1,oft_blagrange_2%nce
        ! Induction
        jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
        + self%dt*DOT_PRODUCT(basis_grads_2(:,jr),basis_grads_2(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        IF (self%region_flag(mesh%reg(i)) == 1 .OR. self%region_flag(mesh%reg(i)) == 3) THEN
            jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
            + basis_vals_2(jr)*basis_vals_2(jc)*jac_det*quad%wts(m)/(eta_t_loc*(coords(1)+gs_epsilon))
        END IF
      END DO
    END DO
  END DO
  !---Apply bc to local matrix
  DO jr=1,oft_blagrange_2%nce
    IF(oft_blagrange_2%be(cell_dofs_2(jr))) jac_loc(6,6)%m(jr,:)=0.d0
  END DO
  IF(.NOT.ASSOCIATED(self%eq)) THEN
    CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%psi_bc(cell_dofs_2),6)
  END IF
  CALL self%fe_rep%mat_add_local(mat,jac_loc,iloc,tlocks)
END DO
deallocate(cell_dofs_1, cell_dofs_2,basis_vals_1, basis_vals_2,basis_grads_1, basis_grads_2,jac_loc, iloc)
!$omp end parallel
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)

! apply free boundary BCs to psi
IF(ASSOCIATED(self%eq)) THEN
  CALL set_bcmat_mod(self%eq,mat, 6,6)
END IF

CALL self%fe_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
END SUBROUTINE build_vac_jacobian

SUBROUTINE build_approx_jacobian(self, mat, a)
class (oft_gs_xmhd_sim), intent(inout) :: self
class (oft_matrix), pointer, intent(inout) :: mat
class(oft_vector), intent(inout) :: a !< Solution for computing jacobian
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals_1, basis_vals_2, p_weights_loc, psi_weights_loc, by_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads_1, basis_grads_2, vel_weights_loc
REAL(r8) :: p, dp(3), vel(3), psi, dpsi(3), by, dby(3), dvel(3,3), div_vel, btmp(3) !reconstructed variables
REAL (r8) :: coords(3), eta_t_loc, eta_p_loc,jac_det, jac_mat(3,4), tmp2(3), tmp3(3)
REAL(r8) :: nu,rho, B_0(3), F0_entry, F0_entry_real(1,1) ! physics parameters
REAL(r8), POINTER, DIMENSION(:) :: p_weights, psi_weights, by_weights
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
REAL(r8), POINTER, DIMENSION(:) ::  vtmp
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
CLASS(oft_vector), POINTER :: oft_lag_vec
INTEGER(i4), ALLOCATABLE, DIMENSION(:), TARGET :: cell_dofs_1, cell_dofs_2
integer (i4) :: i, jr, jc, m, k, l, lim_ind_tmp(1), v
type(oft_quad_type), pointer :: quad
type(oft_1d_int), allocatable, dimension(:) :: iloc
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
LOGICAL :: curved
quad=>oft_blagrange_2%quad
CALL mat%zero
NULLIFY(p_weights,vel_weights, &
         psi_weights, by_weights, vtmp)
!---Get weights from solution vector
CALL a%get_local(p_weights,1)
ALLOCATE(vel_weights(3,oft_blagrange_2%ne))
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(by_weights,5)
CALL a%get_local(psi_weights,6)
!--Set local physics parameters
nu = self%nu
rho = self%rho
B_0 = self%B_0

!--Setup thread locks
ALLOCATE(tlocks(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
  call omp_init_lock(tlocks(i))
END DO
F0_entry = 0.d0
! Declare variables private for OMP
!$omp parallel num_threads(1) private(i,m,jr,jc,curved,coords, cell_dofs_1, cell_dofs_2,basis_vals_1, basis_vals_2, &
!$omp basis_grads_1, basis_grads_2, &
!$omp  jac_loc,jac_mat,jac_det,eta_t_loc, eta_p_loc, &
!$omp p_weights_loc, vel_weights_loc, psi_weights_loc, by_weights_loc, &
!$omp p, dp, vel, dvel, div_vel, psi, dpsi, by, dby, iloc, btmp, tmp2, tmp3, v) reduction(+:F0_entry)
ALLOCATE(basis_vals_1(oft_blagrange_1%nce),basis_grads_1(3,oft_blagrange_1%nce))
ALLOCATE(basis_vals_2(oft_blagrange_2%nce),basis_grads_2(3,oft_blagrange_2%nce))
ALLOCATE(cell_dofs_1(oft_blagrange_1%nce), cell_dofs_2(oft_blagrange_2%nce))
ALLOCATE(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
ALLOCATE(iloc(self%fe_rep%nfields))
ALLOCATE(p_weights_loc(oft_blagrange_1%nce),vel_weights_loc(3, oft_blagrange_2%nce),&
        psi_weights_loc(oft_blagrange_2%nce),&
        by_weights_loc(oft_blagrange_2%nce))
iloc(1)%v=> cell_dofs_1
DO i=2,self%fe_rep%nfields
   iloc(i)%v=>cell_dofs_2
END DO
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do schedule(static)ordered
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange_1%ncdofs(i,cell_dofs_1) ! Get global index of local DOFs
  call oft_blagrange_2%ncdofs(i,cell_dofs_2) ! Get global index of local DOFs
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
  ! Set local weights
  p_weights_loc = p_weights(cell_dofs_1)
  vel_weights_loc = vel_weights(:, cell_dofs_2)
  psi_weights_loc = psi_weights(cell_dofs_2)
  by_weights_loc = by_weights(cell_dofs_2)
!---------------------------------------------------------------------------
! Quadrature Loop
!---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange_1%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange_1,i,jr,quad%pts(:,m),basis_vals_1(jr))
      CALL oft_blag_geval(oft_blagrange_1,i,jr,quad%pts(:,m),basis_grads_1(:,jr),jac_mat)
    END DO
    DO jr=1,oft_blagrange_2%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange_2,i,jr,quad%pts(:,m),basis_vals_2(jr))
      CALL oft_blag_geval(oft_blagrange_2,i,jr,quad%pts(:,m),basis_grads_2(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    basis_grads_1(3, :) = basis_grads_1(2,:)
    basis_grads_1(2,:) = 0.d0
    basis_grads_2(3, :) = basis_grads_2(2,:)
    basis_grads_2(2,:) = 0.d0

    eta_t_loc = self%eta_t(mesh%reg(i))
    eta_p_loc = self%eta_p(mesh%reg(i))
    IF(ALLOCATED(self%eta_node)) THEN
      eta_t_loc = 0.d0 
      DO v = 1, SIZE(mesh%lc, 1)
        eta_t_loc = eta_t_loc + self%eta_node(mesh%lc(v, i)) * quad%pts(v, m)  
      END DO
      eta_p_loc = eta_t_loc
    END IF
    !---Reconstruct values of solution fields
    p = 0.d0; dp = 0.d0; vel = 0.d0; dvel = 0.d0
    psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0

    DO jr=1, oft_blagrange_1%nce
      p = p + p_weights_loc(jr)*basis_vals_1(jr)
      dp = dp + p_weights_loc(jr)*basis_grads_1(:,jr)
    END DO

    DO jr=1,oft_blagrange_2%nce
      vel = vel + vel_weights_loc(:, jr)*basis_vals_2(jr)
      psi = psi + psi_weights_loc(jr)*basis_vals_2(jr)
      by = by + by_weights_loc(jr)*basis_vals_2(jr)
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads_2(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads_2(3, jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads_2(:,jr)
      dby = dby + by_weights_loc(jr)*basis_grads_2(:,jr)
    END DO


    div_vel = dvel(1,1) +vel(1)/(coords(1)+gs_epsilon) + dvel(3,3)
    btmp = cross_product(dpsi/coords(1), [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0]/(coords(1)+gs_epsilon) + B_0
    !---Compute local matrix contributions
    ! VACUUM TERMS
    DO jr=1,oft_blagrange_2%nce
      DO jc=1,oft_blagrange_2%nce
        ! psi, psi
        jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
        + self%dt*DOT_PRODUCT(basis_grads_2(:,jr),basis_grads_2(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        IF (self%region_flag(mesh%reg(i)) == 1 .OR. self%region_flag(mesh%reg(i)) == 3) THEN
            jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
            + basis_vals_2(jr)*basis_vals_2(jc)*jac_det*quad%wts(m)/(eta_t_loc*(coords(1)+gs_epsilon))
        END IF
      END DO
    END DO
  ! EVERYTHING ELSE
    IF (self%region_flag(mesh%reg(i)) == 1) THEN
      !p,p
      IF(ASSOCIATED(self%xml_pre_def))THEN
        DO jr=1,oft_blagrange_1%nce
          DO jc=1,oft_blagrange_1%nce
            jac_loc(1,1)%m(jr,jc) = jac_loc(1,1)%m(jr,jc) &
            + self%dt*DOT_PRODUCT(basis_grads_1(:,jr), basis_grads_1(:,jc))*jac_det*quad%wts(m)*coords(1)/rho
          END DO
        END DO
      END IF

      DO jr=1,oft_blagrange_1%nce
        DO jc=1,oft_blagrange_2%nce
          !p, vel
          DO l=1,3
            jac_loc(1,l+1)%m(jr,jc) = jac_loc(1,l+1)%m(jr,jc) &
            + basis_vals_1(jr)*basis_grads_2(l,jc)*jac_det*quad%wts(m)*coords(1)
          END DO
          jac_loc(1,2)%m(jr,jc) = jac_loc(1,2)%m(jr,jc) &
          + basis_vals_1(jr)*basis_vals_2(jc)*jac_det*quad%wts(m)
        END DO
      END DO
      DO jr=1,oft_blagrange_2%nce
        DO jc=1,oft_blagrange_1%nce
          !vel,p
          DO k=1,3
            jac_loc(k+1,1)%m(jr,jc)= jac_loc(k+1,1)%m(jr,jc) &
              -self%dt*basis_vals_1(jc)*basis_grads_2(k,jr)*jac_det*quad%wts(m)*coords(1)/rho
          END DO    
          jac_loc(2,1)%m(jr,jc)= jac_loc(2,1)%m(jr,jc) &  
            -self%dt*basis_vals_1(jc)*basis_vals_2(jr)*jac_det*quad%wts(m)/rho 
        END DO
      END DO
      DO jr=1,oft_blagrange_2%nce
        DO jc=1,oft_blagrange_2%nce
          ! vel,vel
          DO k=1,3
            jac_loc(k+1,k+1)%m(jr,jc)= jac_loc(k+1,k+1)%m(jr,jc) &
              + basis_vals_2(jr)*basis_vals_2(jc)*jac_det*quad%wts(m)*coords(1) &
              + self%dt*basis_vals_2(jr)*DOT_PRODUCT(vel, basis_grads_2(:,jc))*jac_det*quad%wts(m)*coords(1) &
              + self%dt*nu*DOT_PRODUCT(basis_grads_2(:,jr), basis_grads_2(:,jc))*jac_det*quad%wts(m)*coords(1)/rho
            DO l=1,3
              jac_loc(k+1,l+1)%m(jr,jc)= jac_loc(k+1,l+1)%m(jr,jc) &
              + self%dt*basis_vals_2(jr)*basis_vals_2(jc)*dvel(k, l)*jac_det*quad%wts(m)*coords(1)
            END DO
          END DO
          jac_loc(2,2)%m(jr,jc)= jac_loc(2,2)%m(jr,jc) &
            + self%dt*nu*basis_vals_2(jr)*basis_vals_2(jc)*jac_det*quad%wts(m)/(rho*(coords(1)+gs_epsilon))
          jac_loc(2,3)%m(jr,jc)= jac_loc(2,3)%m(jr,jc) &
            -self%dt*basis_vals_2(jr)*2.d0*vel(2)*basis_vals_2(jc)*jac_det*quad%wts(m)
          jac_loc(3,2)%m(jr,jc)= jac_loc(3,2)%m(jr,jc) &
            +self%dt*basis_vals_2(jr)*basis_vals_2(jc)*vel(2)*jac_det*quad%wts(m)
          jac_loc(3,3)%m(jr,jc)= jac_loc(3,3)%m(jr,jc) &
            + self%dt*basis_vals_2(jr)*basis_vals_2(jc)*vel(1)*jac_det*quad%wts(m) &
            + self%dt*nu*basis_vals_2(jr)*basis_vals_2(jc)*jac_det*quad%wts(m)/(rho*(coords(1)+gs_epsilon))
          !vel, by
          tmp2 = 0.d0
          tmp2 = [0.d0,basis_vals_2(jc)/(coords(1)+gs_epsilon),0.d0]! this is 'delta B_y'
          DO l=1,3
            jac_loc(l+1,5)%m(jr,jc) = jac_loc(l+1,5)%m(jr,jc) &
            + self%dt*DOT_PRODUCT(basis_grads_2(:,jr),tmp2)*btmp(l)*jac_det*quad%wts(m)*coords(1)/(rho*mu0) &
            + self%dt*DOT_PRODUCT(basis_grads_2(:,jr), btmp)*tmp2(l)*jac_det*quad%wts(m)*coords(1)/(rho*mu0) &
            - self%dt*basis_grads_2(l,jr)*DOT_PRODUCT(tmp2, btmp)*jac_det*quad%wts(m)*coords(1)/(rho*mu0)
          END DO
          jac_loc(2,5)%m(jr,jc) = jac_loc(2,5)%m(jr,jc) &
          + self%dt*basis_vals_2(jr)*btmp(2)*tmp2(2)*jac_det*quad%wts(m)/(rho*mu0)
          jac_loc(3,5)%m(jr,jc) = jac_loc(3,5)%m(jr,jc) &
          - self%dt*basis_vals_2(jr)*btmp(1)*tmp2(2)*jac_det*quad%wts(m)/(rho*mu0)
          !vel, psi
          tmp2 = cross_product(basis_grads_2(:,jc), [0.d0,1.d0/(coords(1)+gs_epsilon),0.d0]) ! this is 'dB'
          !DO l=1,3
            !jac_loc(l+1,6)%m(jr,jc) = jac_loc(l+1,6)%m(jr,jc) &
            !+ self%dt*DOT_PRODUCT(basis_grads_2(:,jr), tmp2)*btmp(l)*jac_det*quad%wts(m)*coords(1)/(rho*mu0) &
            !+ self%dt*DOT_PRODUCT(basis_grads_2(:,jr), btmp)*tmp2(l)*jac_det*quad%wts(m)*coords(1)/(rho*mu0) &
            !- self%dt*basis_grads_2(l,jr)*DOT_PRODUCT(tmp2, btmp)*jac_det*quad%wts(m)*coords(1)/(rho*mu0)
          !END DO
          ! jac_loc(2,6)%m(jr,jc) = jac_loc(2,6)%m(jr,jc) &
          ! + self%dt*basis_vals_2(jr)*(-btmp(1)*tmp2(1)-btmp(3)*tmp2(3))*jac_det*quad%wts(m)/(rho*mu0)
          ! jac_loc(3,6)%m(jr,jc) = jac_loc(3,6)%m(jr,jc) &
          ! - self%dt*basis_vals_2(jr)*(btmp(2)*tmp2(1))*jac_det*quad%wts(m)/(rho*mu0)
          ! psi, psi
          jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
            + self%dt*basis_vals_2(jr)*DOT_PRODUCT(vel,basis_grads_2(:,jc))*jac_det*quad%wts(m)/(eta_t_loc*(coords(1)+gs_epsilon))
          !psi, vel
          DO l=1,3
            tmp2 = [0.d0, 0.d0,0.d0]
            jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
            + self%dt*basis_vals_2(jr)*basis_vals_2(jc)*dpsi(l)*jac_det*quad%wts(m)/(eta_t_loc*(coords(1)+gs_epsilon))
            tmp2(l) = 1.d0
            tmp3 = cross_product(B_0, tmp2)
            jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
            + self%dt*basis_vals_2(jr)*basis_vals_2(jc)*tmp3(2)*jac_det*quad%wts(m)/(eta_t_loc*(coords(1)+gs_epsilon))
          END DO
          !By, By (everything except diffusion term)
          jac_loc(5, 5)%m(jr,jc) = jac_loc(5, 5)%m(jr,jc) &
          + basis_vals_2(jr)*self%dt*dvel(1,1)*basis_vals_2(jc)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + basis_vals_2(jr)*self%dt*dvel(3,3)*basis_vals_2(jc)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + basis_vals_2(jr)*self%dt*DOT_PRODUCT(vel, basis_grads_2(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)&
          + basis_vals_2(jr)*self%dt*vel(1)*basis_vals_2(jc)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)**2
          !By, vel
          jac_loc(5,2)%m(jr,jc) = jac_loc(5,2)%m(jr,jc) &
          + basis_vals_2(jr)*self%dt*basis_vals_2(jc)*dby(1)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + basis_vals_2(jr)*self%dt*basis_vals_2(jc)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)**2 &
          + basis_vals_2(jr)*self%dt*basis_grads_2(1,jc)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
          jac_loc(5,4)%m(jr,jc) = jac_loc(5,4)%m(jr,jc) &
          + basis_vals_2(jr)*self%dt*basis_vals_2(jc)*dby(3)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + basis_vals_2(jr)*self%dt*basis_grads_2(3,jc)*by*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
          tmp2 = cross_product(dpsi,basis_grads_2(:,jc))
          jac_loc(5,3)%m(jr,jc) = jac_loc(5,3)%m(jr,jc) &
           - basis_vals_2(jr)*self%dt*tmp2(2)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
          !By, psi
          tmp2 = cross_product(basis_grads_2(:,jc),dvel(2,:))
          jac_loc(5, 6)%m(jr,jc) = jac_loc(5, 6)%m(jr,jc) &
           - basis_vals_2(jr)*self%dt*tmp2(2)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
        END DO
      END DO
    END IF
    ! ADD F DIFFUSION TERMS EVERYWHERE
    DO jr=1,oft_blagrange_2%nce
      DO jc=1,oft_blagrange_2%nce
          !By, By
          jac_loc(5, 5)%m(jr,jc) = jac_loc(5, 5)%m(jr,jc) &
          + basis_vals_2(jr)*basis_vals_2(jc)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon) &
          + self%dt*eta_p_loc*DOT_PRODUCT(basis_grads_2(:,jr), basis_grads_2(:,jc))*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
      END DO
    END DO
    !Compute approximate F0/F0 jacobian factor (very approximate for now)
    IF (self%region_flag(mesh%reg(i)) == 5) THEN
      F0_entry = F0_entry + jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
    END IF
  END DO

!---Apply bc to local matrix
  DO jr=1,oft_blagrange_2%nce
    IF(oft_blagrange_2%be(cell_dofs_2(jr))) jac_loc(6,6)%m(jr,:)=0.d0
  END DO
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%p_bc(cell_dofs_1),1)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velx_bc(cell_dofs_2),2)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%vely_bc(cell_dofs_2),3)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velz_bc(cell_dofs_2),4)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%by_bc(cell_dofs_2),5)
  IF(.NOT.ASSOCIATED(self%eq)) THEN
    CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%psi_bc(cell_dofs_2),6)
  END IF
  CALL self%fe_rep%mat_add_local(mat,jac_loc,iloc,tlocks)
END DO
deallocate(cell_dofs_1, cell_dofs_2,basis_vals_1, basis_vals_2,basis_grads_1, basis_grads_2,jac_loc, iloc)
deallocate(p_weights_loc, vel_weights_loc, psi_weights_loc, by_weights_loc)
!$omp end parallel

IF (.NOT. any(self%region_flag== 5)) THEN
  F0_entry = F0_entry + self%lim_vac_int     
END IF
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
! apply free boundary BCs to psi
IF (ASSOCIATED(self%eq)) THEN
  CALL set_bcmat_mod(self%eq,mat, 6,6)
END IF
CALL fem_dirichlet_diag(oft_blagrange_1,mat,self%p_bc,1)
CALL fem_dirichlet_diag(oft_blagrange_2,mat,self%velx_bc,2)
CALL fem_dirichlet_diag(oft_blagrange_2,mat,self%vely_bc,3)
CALL fem_dirichlet_diag(oft_blagrange_2,mat,self%velz_bc,4)
IF(.NOT.ASSOCIATED(self%eq)) THEN
  CALL fem_dirichlet_diag(oft_blagrange_2,mat,self%psi_bc,6)
END IF

IF (self%evolve_F) THEN
  !Set entries for F in plasma
  CALL fem_dirichlet_diag(oft_blagrange_2,mat,self%plasma_bc,5)
  CALL set_f0mat(self, mat, 5, 5)
  lim_ind_tmp = self%lim_ind
  F0_entry_real = F0_entry
  CALL mat%add_values(lim_ind_tmp,lim_ind_tmp,F0_entry_real,1,1, 5, 5)
ELSE
  CALL fem_dirichlet_diag(oft_blagrange_2,mat,self%by_bc,5)
END IF

CALL self%fe_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec, p_weights, vel_weights, psi_weights, by_weights)

END SUBROUTINE build_approx_jacobian

!------------------------------------------------------------------------------
!> creates matrix with dense blocks on boundaries for mask = 3 (mask = 0 => nothing, mask = 1 => normal, mask = 2 => identity)
! mask = 4 -> couple plasma nodes to one boundary node
!------------------------------------------------------------------------------
subroutine fem_mat_create_mod(self,new,mask)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
CLASS(oft_matrix), POINTER, INTENT(out) :: new
INTEGER(i4), OPTIONAL, INTENT(in) :: mask(:,:)
INTEGER(i4) :: i,j,k,nknown_graphs, l
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: mat_mask,graph_ids
CLASS(oft_vector), POINTER :: tmp_vec
TYPE(oft_graph_ptr), ALLOCATABLE :: graphs(:,:),known_graphs(:)
TYPE(oft_graph), TARGET :: dense_graph
type(oft_1d_int), pointer, dimension(:) :: bc_nodes, F0_node
integer(i4), allocatable :: dense_flag(:), plasma_flag(:)
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
    IF (mat_mask(i,j)==4) THEN
      ALLOCATE(F0_node(1))
      F0_node(1)%n = 1
      ALLOCATE(F0_node(1)%v(1))
      F0_node(1)%v(1) = current_sim%lim_ind
      ALLOCATE(plasma_flag(self%fields(i)%fe%ne))
      plasma_flag = 0
      DO l=1, self%fields(i)%fe%ne
        IF (current_sim%plasma_bc(l)) THEN !if in plasma region
          plasma_flag(l) = 1
        END IF
      END DO
      !---Add dense blocks
      CALL graph_add_dense_blocks(graphs(i,j)%g,dense_graph,plasma_flag,F0_node)      
      NULLIFY(graphs(i,j)%g%kr,graphs(i,j)%g%lc)
      graphs(i,j)%g%nnz=dense_graph%nnz
      graphs(i,j)%g%kr=>dense_graph%kr
      graphs(i,j)%g%lc=>dense_graph%lc
      DEALLOCATE(F0_node, plasma_flag)
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

subroutine update_bcs(self)
class(oft_gs_xmhd_sim), intent(inout) :: self
integer(i4) :: i, type
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs_1, cell_dofs_2
LOGICAL, ALLOCATABLE :: p_dir_set(:)
IF (ALLOCATED(self%region_flag)) THEN
  ALLOCATE(cell_dofs_1(oft_blagrange_1%nce))
  ALLOCATE(cell_dofs_2(oft_blagrange_2%nce))
  ALLOCATE(self%p_bc(oft_blagrange_1%ne)); self%p_bc=.TRUE.
  ALLOCATE(self%velx_bc(oft_blagrange_2%ne)); self%velx_bc=.FALSE.
  ALLOCATE(self%vely_bc(oft_blagrange_2%ne)); self%vely_bc=.TRUE.
  ALLOCATE(self%velz_bc(oft_blagrange_2%ne)); self%velz_bc=.FALSE.
  ALLOCATE(self%by_bc(oft_blagrange_2%ne)); self%by_bc=.FALSE.  ! FOR NOW WE'RE NOT EVOLVING By (F)
  ALLOCATE(self%psi_bc(oft_blagrange_2%ne)); self%psi_bc=.FALSE. 
  ALLOCATE(self%plasma_bc(oft_blagrange_2%ne)); self%plasma_bc=.FALSE.  ! FOR NOW WE'RE NOT EVOLVING By (F)
  IF (SIZE(self%region_flag) /= mesh%nreg) THEN
    CALL oft_abort("Number of region flags does not match number of regions.","setup",__FILE__)
  END IF
  DO i=1, mesh%nc
    type = self%region_flag(mesh%reg(i))
    IF (type == 1) THEN
      CALL apply_mhd_bcs(self, i, cell_dofs_1, cell_dofs_2)
    ELSE IF(type ==5 .OR. type ==6) THEN
      CALL apply_plasma_bcs(self, i, cell_dofs_1, cell_dofs_2)
    ELSE IF (type >1 .AND. type < 5) THEN
      CALL apply_bcs(self, i, cell_dofs_1, cell_dofs_2)
    ELSE
      CALL oft_abort("Invalid region flag.","setup",__FILE__)
    END IF
  END DO
  ! Pin one node per MHD region after all BCs are set
  ALLOCATE(p_dir_set(mesh%nreg))
  p_dir_set = .FALSE.
  DO i=1, mesh%nc
    type = self%region_flag(mesh%reg(i))
    IF (type == 1 .AND. .NOT. p_dir_set(mesh%reg(i))) THEN
      call oft_blagrange_1%ncdofs(i,cell_dofs_1)
      self%p_bc(cell_dofs_1(1)) = .TRUE.
      p_dir_set(mesh%reg(i)) = .TRUE.
    END IF
  END DO
END IF
IF (.NOT. self%evolve_F) THEN
  self%by_bc = .TRUE.
END IF 
DEALLOCATE(cell_dofs_1, cell_dofs_2)
self%nlfun%p_bc=>self%p_bc
self%nlfun%velx_bc=>self%velx_bc
self%nlfun%vely_bc=>self%vely_bc
self%nlfun%velz_bc=>self%velz_bc
self%nlfun%by_bc=>self%by_bc
self%nlfun%plasma_bc => self%plasma_bc
end subroutine update_bcs

subroutine apply_mhd_bcs(self,cell_ind, cell_dofs_a, cell_dofs_b)
class(oft_gs_xmhd_sim), intent(inout) :: self
integer(i4) , intent(in) :: cell_ind
INTEGER(i4), POINTER, DIMENSION(:), intent(inout) :: cell_dofs_a, cell_dofs_b
INTEGER(i4) :: j
call oft_blagrange_2%ncdofs(cell_ind,cell_dofs_b) ! Get global index of local DOFs
call oft_blagrange_1%ncdofs(cell_ind,cell_dofs_a) ! Get global index of local DOFs
! DO j=1, SIZE(cell_dofs_b)
!   !self%velx_bc(cell_dofs_b(j)) = .TRUE. ! prevent psi evolution in superconductor
!   !self%vely_bc(cell_dofs_b(j)) = .TRUE. ! prevent psi evolution in superconductor
!   !self%velz_bc(cell_dofs_b(j)) = .TRUE. ! prevent psi evolution in superconductor
! END DO
DO j=1, SIZE(cell_dofs_a)
  self%p_bc(cell_dofs_a(j)) = .FALSE. ! prevent velocity evolution in solid conductor
END DO
end subroutine apply_mhd_bcs

!---------------------------------------------------------------------------
!> Apply boundary conditions for non-extended MHD regions (plasma, coils, solid conductors, vacuum)
!---------------------------------------------------------------------------
subroutine apply_bcs(self,cell_ind, cell_dofs_a, cell_dofs_b)
class(oft_gs_xmhd_sim), intent(inout) :: self
INTEGER(i4) , intent(in) :: cell_ind
INTEGER(i4), POINTER, DIMENSION(:), intent(inout) :: cell_dofs_a, cell_dofs_b
INTEGER(i4) :: j
call oft_blagrange_1%ncdofs(cell_ind,cell_dofs_a) ! Get global index of local DOFs
call oft_blagrange_2%ncdofs(cell_ind,cell_dofs_b) ! Get global index of local DOFs
DO j=1, SIZE(cell_dofs_b)
  self%velx_bc(cell_dofs_b(j)) = .TRUE. ! prevent velocity evolution in solid conductor
  self%vely_bc(cell_dofs_b(j)) = .TRUE. ! prevent velocity evolution in solid conductor
  self%velz_bc(cell_dofs_b(j)) = .TRUE. ! prevent velocity evolution in solid conductor
END DO
! DO j=1, SIZE(cell_dofs_a)
!   self%p_bc(cell_dofs_a(j)) = .TRUE. ! prevent velocity evolution in solid conductor
! END DO
end subroutine apply_bcs


!---------------------------------------------------------------------------
!> Apply boundary conditions for non-extended MHD regions (plasma, coils, solid conductors, vacuum)
!---------------------------------------------------------------------------
subroutine apply_plasma_bcs(self,cell_ind, cell_dofs_a, cell_dofs_b)
class(oft_gs_xmhd_sim), intent(inout) :: self
INTEGER(i4) , intent(in) :: cell_ind
INTEGER(i4), POINTER, DIMENSION(:), intent(inout) :: cell_dofs_a, cell_dofs_b
INTEGER(i4) :: j
call oft_blagrange_1%ncdofs(cell_ind,cell_dofs_a) ! Get global index of local DOFs
call oft_blagrange_2%ncdofs(cell_ind,cell_dofs_b) ! Get global index of local DOFs
DO j=1, SIZE(cell_dofs_b)
  self%velx_bc(cell_dofs_b(j)) = .TRUE. ! prevent velocity evolution in solid conductor
  self%vely_bc(cell_dofs_b(j)) = .TRUE. ! prevent velocity evolution in solid conductor
  self%velz_bc(cell_dofs_b(j)) = .TRUE. ! prevent velocity evolution in solid conductor
  self%by_bc(cell_dofs_b(j)) = .TRUE. ! prevent psi evolution in superconductor
  self%plasma_bc(cell_dofs_b(j)) = .TRUE. ! prevent psi evolution in superconductor
END DO
! DO j=1, SIZE(cell_dofs_a)
!   self%p_bc(cell_dofs_a(j)) = .TRUE. ! prevent velocity evolution in solid conductor
! END DO
self%plasma_bc(self%lim_ind) = .FALSE. ! Turn off Dirichlet conditions for node determining value of F0
end subroutine apply_plasma_bcs

SUBROUTINE r_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = pt(1)
END SUBROUTINE r_init

!------------------------------------------------------------------------------
!> Add F0 constraint couplings
!------------------------------------------------------------------------------
subroutine set_f0mat(self, mat, iblock, jblock)
class(oft_gs_xmhd_sim), intent(inout) :: self
class(oft_matrix), intent(inout) :: mat !< Matrix object
integer(4), intent(in) :: iblock, jblock
integer(4) :: i, i_inds(1),j_inds(1)
real(8) :: one_val(1,1)
!---Add to matrix
! | A_ii A_ib |
! | M_bi M_bb + M*L^-1*M |
one_val=-1.d0
DO i=1, oft_blagrange_2%ne
  IF(self%plasma_bc(i)) THEN
    i_inds = i 
    j_inds = self%lim_ind
    CALL mat%add_values(i_inds,j_inds,one_val,1,1, iblock, jblock)
  END IF
END DO
end subroutine set_f0mat

END MODULE gs_xmhd_v8