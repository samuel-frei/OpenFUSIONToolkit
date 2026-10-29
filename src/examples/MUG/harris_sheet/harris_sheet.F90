!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: Harris Sheet Reconnection    {#doc_mug2d_harris_ex}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of the \ref xmhd2d "extended MHD" module in OFT. In
!! this example reconnection in a periodic-slab current sheet will be simulated following
!! the setup used to the [GEM Challange](https://doi.org/10.1029/1999JA900449).
!!
!!\section doc_mug2d_harris_ex_code_helper Code Walk Through
!!
!! The code consists of the main run program and subroutines for field initialization.
!!
!! By default the run starts from the analytic Harris sheet plus a cosine flux seed of
!! relative size `delta`. Setting `init_rst` instead starts it from a restart file:
!!
!! | `init_rst_pert` | `linear` | initial state |
!! |---|---|---|
!! | `T` (default) | `F` | Harris equilibrium + `init_rst_scale` x file (no cosine seed) |
!! | `T` (default) | `T` | perturbation `init_rst_scale` x file, evolved about the equilibrium |
!! | `F` | `F` | the file itself, as a full state (e.g. a nonlinear run's own restart) |
!!
!! A perturbation file is what a linear run writes (`xmhd2d_00010.rst`), or an eigenmode or
!! pseudomode written by nonmodal analysis: the same `U_*` field layout either way. Time and
!! restart numbering start from zero in every case.
!!
!!\section doc_mug2d_harris_ex_code_helper Main Code
!!
!! As usual, we start with a few module imports for required functionality.
! START SOURCE
PROGRAM xmhd_tearing
!---Runtime
USE oft_base
USE oft_io, ONLY: oft_file_exist
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE oft_blag_operators, ONLY: oft_blag_zerob, oft_blag_getmop, oft_blag_getlop, oft_blag_project
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE xmhd_2d
IMPLICIT NONE
#include "local.h"
INTEGER(i4) :: io_unit,ierr
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_blag_zerob), TARGET :: blag_zerob ! setting boundary vals to zero
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: mop_bc => NULL() !< Mass matrix exported for nonmodal analysis
CLASS(oft_matrix), POINTER :: lop_bc => NULL() !< Laplacian stiffness matrix exported for nonmodal analysis
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_vector), POINTER :: u_rst => NULL() !< State read from `init_rst`
LOGICAL :: use_rst !< Start from `init_rst` rather than the analytic seed
!!\subsection doc_mug2d_recon_ex_code_vars Local Variables
!! Next we define the local variables needed to initialize our case and
!! run the time-dependent solve and post-processing. Compared to previous
!! examples, we now have specialized initial condition and post-processing
!! objects.
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 100
INTEGER(i4) :: rst_freq = 10
REAL(r8) :: B_0(3)=0.d0  !< Background magnetic field
REAL(r8) :: n0 = 1.d0 !< Density scale
REAL(r8) :: velx0 = 1.d0 !< Velocity scale (x-direction)
REAL(r8) :: vely0 = 1.d0 !< Velocity scale (y-direction)
REAL(r8) :: velz0 = 1.d0 !< Velocity scale (z-direction)
REAL(r8) :: t0 = 1.d0 !< Temperature scale
REAL(r8) :: psi0 = 1.d0 !< Magnetic flux scale
REAL(r8) :: by0 = 1.d0 !< Magnetic field scale (y-direction)
REAL(r8) :: chi=1.d0 !< thermal diffusivity
REAL(r8) :: eta=1.d0 !< electrical resistivity
REAL(r8) :: nu=1.d0 !< viscosity
REAL(r8) :: gamma=2.d0/3.d0 !< adiabatic index
REAL(r8) :: D_diff=1.d0 !< Diffusion coefficient
REAL(r8) :: den_scale=1.d19 !< Density scale factor for normalization of equations (used for post-processing)
REAL(r8) :: den_inf = 2.d-1 !< Background density
REAL(r8) :: lam_n = 2.d0 !< eq density width
REAL(r8) :: lam_b = 2.d0 !< eq psi perturbation width
REAL(r8) :: L_x = 2.d0 !< Box length in x
REAL(r8) :: L_z = 2.d0 !< Box length in z
REAL(r8) :: delta = 1.d-1 !< Relative size of perturbation (<<1)
REAL(r8) :: dt = 1.d-3
LOGICAL :: pm=.FALSE.
LOGICAL :: use_mfnk=.FALSE.
LOGICAL :: linear=.FALSE.
LOGICAL :: plot_run=.FALSE.
CHARACTER(LEN=OFT_PATH_SLEN) :: init_rst = '' !< Restart file holding the initial state ('' = analytic seed)
LOGICAL :: init_rst_pert = .TRUE. !< `init_rst` holds a perturbation about the Harris equilibrium
REAL(r8) :: init_rst_scale = 1.d0 !< Factor applied to a perturbation read from `init_rst`
NAMELIST/harris_sheet_options/order, chi, eta, nu, gamma, D_diff, &
dt, nsteps, rst_freq, use_mfnk, pm, linear, n0, psi0, velx0, vely0, velz0, &
t0, by0, den_scale, delta, lam_n, lam_b, L_x, L_z,plot_run, &
init_rst, init_rst_pert, init_rst_scale
!!\subsection doc_mug_recon_ex_driver_setup OFT Initialization
!! See \ref doc_api_ex1 for a detailed description of calls in this section.
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,harris_sheet_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr/=0)CALL oft_abort('Error reading harris_sheet_options','xmhd_tearing',__FILE__)
use_rst=LEN_TRIM(init_rst)>0
IF(use_rst)THEN
  IF(.NOT.oft_file_exist(TRIM(init_rst)))CALL oft_abort('init_rst file "'//TRIM(init_rst)// &
    '" not found','xmhd_tearing',__FILE__)
  IF(linear.AND.(.NOT.init_rst_pert))CALL oft_abort('init_rst_pert=F (a full state) needs '// &
    'linear=F: a linear run evolves a perturbation','xmhd_tearing',__FILE__)
  IF((.NOT.init_rst_pert).AND.(init_rst_scale/=1.d0))CALL oft_abort('init_rst_scale applies '// &
    'only to a perturbation (init_rst_pert=T)','xmhd_tearing',__FILE__)
END IF
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh,order)
blag_zerob%ML_lag_rep=>ML_oft_blagrange

!---------------------------------------------------------------------------
! Set intial conditions from analytic functions
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(u,v,mop,vec_vals) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop) ! Construct mass matrix with "none" BC
!---Linear runs export the scalar mass and Laplacian stiffness matrices next to
!---the Jacobian (lin_ops.h5). Both zero the same Dirichlet boundary rows.
!---Each goes in its own file because mat_save truncates the file it writes.
IF(linear)THEN
  CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop_bc,ML_oft_blagrange%current_level%global%gbe)
  CALL mop_bc%save('mass_mat.h5','massmat')
  CALL mop_bc%delete
  DEALLOCATE(mop_bc)
  CALL oft_blag_getlop(ML_oft_blagrange%current_level,lop_bc,'zerob')
  CALL lop_bc%save('stiff_mat.h5','stiffmat')
  CALL lop_bc%delete
  DEALLOCATE(lop_bc)
END IF
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)

IF (plot_run) THEN
  write(*,*) "Plotting initial conditions"
  mhd_sim%rst_freq=rst_freq
  mhd_sim%linear=linear
  mhd_sim%den_scale=den_scale ! the plot rescales density by it
  CALL xmhd_2d_plot(mhd_sim)
  CALL oft_finalize
END IF

!---Project n initial condition onto scalar Lagrange basis
field_init%func=>dens_harris
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(n0)
CALL u%get_local(vec_vals)
vec_vals = vec_vals / den_scale
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,1)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,1)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,1)
END IF

!---Set v_x initial condition
CALL u%set(velx0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,2)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,2)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,2)
END IF

!---Set v_y initial condition
CALL u%set(vely0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,3)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,3)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,3)
END IF

!---Set v_z initial condition
CALL u%set(velz0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,4)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,4)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,4)
END IF

!---Set T initial condition
CALL u%set(t0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,5)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,5)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,5)
END IF

!---Project psi initial condition onto scalar Lagrange basis
field_init%func=>psi_harris_eq
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(psi0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,6)
ELSE
    CALL mhd_sim%u%restore_local(vec_vals,6)
END IF
!---Analytic flux seed, unless the perturbation comes from init_rst instead
IF(.NOT.use_rst)THEN
  field_init%func=>psi_harris_pert
  CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  CALL u%scale(psi0)
  IF(linear)THEN
      CALL u%get_local(vec_vals)
      CALL mhd_sim%u%restore_local(vec_vals,6)
  ELSE
      CALL mhd_sim%u%get_local(vec_vals,6)
      CALL v%restore_local(vec_vals)
      CALL u%add(1.d0, 1.d0, v)
      CALL u%get_local(vec_vals)
      CALL mhd_sim%u%restore_local(vec_vals,6)
  END IF
END IF

!---Set By initial condition
CALL u%set(by0)
CALL u%get_local(vec_vals)
IF(linear)THEN
    CALL mhd_sim%u0%restore_local(vec_vals,7)
    vec_vals = 0.d0
     CALL mhd_sim%u%restore_local(vec_vals,7)
ELSE
     CALL mhd_sim%u%restore_local(vec_vals,7)
END IF

!---Initial state from a restart file (see the table at the top of this file)
IF(use_rst)THEN
  CALL mhd_sim%fe_rep%vec_create(u_rst)
  CALL mhd_sim%rst_load(u_rst,TRIM(init_rst),'U')
  IF(.NOT.init_rst_pert)THEN
    CALL mhd_sim%u%add(0.d0,1.d0,u_rst) ! Full state replaces the analytic one
  ELSE IF(linear)THEN
    CALL mhd_sim%u%add(0.d0,init_rst_scale,u_rst) ! Perturbation about u0
  ELSE
    CALL mhd_sim%u%add(1.d0,init_rst_scale,u_rst) ! Equilibrium + perturbation
  END IF
  IF(oft_env%head_proc)WRITE(*,'(2X,3A,L1,A,ES12.4)')'Initial state from "',TRIM(init_rst), &
    '": perturbation=',init_rst_pert,', scale=',init_rst_scale
  CALL u_rst%delete
  DEALLOCATE(u_rst)
END IF

!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)

!---------------------------------------------------------------------------
! Run simulation
!---------------------------------------------------------------------------
mhd_sim%chi=chi
mhd_sim%eta=eta
mhd_sim%nu=nu
mhd_sim%B_0=B_0
mhd_sim%gamma=gamma
mhd_sim%D_diff=D_diff
mhd_sim%dt=dt
mhd_sim%den_scale=den_scale
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%mfnk=use_mfnk
mhd_sim%linear=linear
oft_env%pm=pm

IF(linear)THEN
  CALL mhd_sim%run_lin_simulation()
ELSE
  CALL mhd_sim%run_simulation()
END IF

CALL xmhd_2d_plot(mhd_sim)
!---Finalize enviroment
CALL oft_finalize
CONTAINS

!! To set the initial conditions we define a set of functions
!!
!! The non-uniform initial conditions for this case is given by
!! \f[ n_0 = (cosh(z / \lambda))^{-2} + n_{\infty} \f]
!! \f[ B_0 = tanh(z / \lambda) \hat{x} \f]
!!
!! The perturbation is then given by
!! \f[ \delta B = -\frac{\pi}{L_z} cos(2 \pi x / L_x) sin(\pi z / L_z) \hat{x} + \frac{2 \pi}{L_x} sin(2 \pi x / L_x) cos(\pi z / L_z) \hat{x} \f]

SUBROUTINE dens_harris(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = den_inf + 1/((COSH(pt(2)/lam_n))**2)
END SUBROUTINE dens_harris

SUBROUTINE psi_harris_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = -lam_b*LOG(COSH(pt(2)/lam_b))
END SUBROUTINE psi_harris_eq

SUBROUTINE psi_harris_pert(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = - delta*COS(2*pi*pt(1)/L_x)*COS(pi*pt(2)/L_z)
END SUBROUTINE psi_harris_pert

END PROGRAM xmhd_tearing
