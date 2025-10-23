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
USE oft_lu, ONLY: oft_lusolver
!
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec, bfem_map_flag
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE oft_blag_operators, ONLY: oft_blag_vproject,oft_blag_project, oft_blag_getmop, oft_lag_bginterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE oft_gs, ONLY: gs_epsilon, build_dels, build_dels_mug, gs_eq, gs_update_bounds, gs_test_bounds
USE oft_gs_td, ONLY: oft_tmaker_td_mfop, tMaker_td_mfnk_update, build_vac_op, apply_rhs
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
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< psi BC flag
  LOGICAL :: pm = .FALSE.
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
  TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  TYPE(xdmf_plot_file) :: xdmf_plot
  TYPE(gs_eq), POINTER :: eq => NULL() !< Equilibrium object
  CLASS(oft_vector), POINTER :: u => NULL() !< current solution vector
  CLASS(oft_vector), POINTER :: rhs => NULL() !< Temporary RHS vector
  CLASS(oft_vector), POINTER :: psi_tmp => NULL() !< Temporary storage vector
  CLASS(oft_vector), POINTER :: tmp_vec => NULL() !< Temporary storage vector
  TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix free operator
  TYPE(oft_tmaker_td_mfop), POINTER :: mfop => NULL() !< Time-advance operator
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
CONTAINS


subroutine setup(self)
class(oft_gs_xmhd_sim), intent(inout) :: self !< NL operator object
integer(i4) :: ierr
integer(i4) :: order = 2
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
! Create operator
!------------------------------------------------------------------------------
ALLOCATE(self%mfop)
self%mfop%dt=self%dt
CALL self%mfop%setup(self%eq)

!------------------------------------------------------------------------------
! Create Solver fields
!------------------------------------------------------------------------------

self%u=>self%mfop%gs_eq%psi
call self%eq%fe_rep%vec_create(self%rhs)
call self%eq%fe_rep%vec_create(self%psi_tmp)

ALLOCATE(self%vac_pre)
self%vac_pre%A=>self%mfop%vac_op
!
!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat)
CALL self%mfmat%setup(self%psi_tmp,self%mfop)
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
self%nksolver%A=>self%mfop
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
    CALL self%mfop%update()
    ! Update operators if the timestep has changed
    IF(self%dt/=self%mfop%dt)THEN
        self%dt=ABS(self%dt)
        self%mfop%dt=self%dt
        CALL build_vac_op(self%mfop,self%mfop%vac_op)
        CALL self%vac_pre%update(.TRUE.)
    END IF
    !Build right hand side
    CALL self%psi_tmp%add(0.d0,1.d0,self%u)
    CALL apply_rhs(self%mfop,self%u,self%rhs)
    CALL self%mfop%gs_eq%zerob_bc%apply(self%rhs)
    ! Do nonlinear solve
    DO j=1,4
        CALL self%nksolver%apply(self%u,self%rhs)
        IF(self%nksolver%cits<0)THEN
            CALL self%u%add(0.d0,1.d0,self%psi_tmp)
            self%mfop%dt=self%mfop%dt/2.d0
            CALL build_vac_op(self%mfop,self%mfop%vac_op)
            CALL self%vac_pre%update(.TRUE.)
            CALL apply_rhs(self%mfop,self%u,self%rhs)
            CALL self%mfop%gs_eq%zerob_bc%apply(self%rhs)
            CYCLE
        ELSE
            EXIT
        END IF
    END DO
    self%t=self%t+self%mfop%dt
    self%dt=self%mfop%dt
    self%mfop%gs_eq%alam=self%mfop%f_scale
    self%mfop%gs_eq%pnorm=self%mfop%p_scale

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
write(*,*) self%mfop%gs_eq%alam
end subroutine run_simulation

SUBROUTINE gs_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a
CALL current_sim%mfmat%update(a)
END SUBROUTINE gs_mfnk_update

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

END MODULE gs_xmhd