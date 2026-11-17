!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file xmhd_2d_f.F90
!
!> Fortran part of the Python wrapper for the 2D MUG xMHD model (`xmhd_2d`)
!!
!! Exposes the pieces an optimisation or adjoint loop needs without going through
!! files or a fresh executable per evaluation:
!! - the state vector (local layout, with the global index of every entry)
!! - the time-advance residual `F_dt(u) = M u + dt N(u)` and its Jacobian, in CSR
!! - fixed-timestep nonlinear advance with OFT's Newton solver, optionally
!!   recording every state of the trajectory
!! - restart files in the format `xmhd_2d` itself writes
!!
!! Only one simulation can exist per process, because `xmhd_2d` keeps its finite
!! element spaces in module variables.
!!
!! @ingroup doxy_oft_python
!---------------------------------------------------------------------------------
MODULE xmhd_2d_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_ptr, c_f_pointer, &
  c_bool, c_associated
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, vector_extrapolate
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_cast
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre
USE xmhd_2d, ONLY: oft_xmhd_2d_sim, build_approx_jacobian
USE oft_base_f, ONLY: copy_string, copy_string_rev
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------------
!> 2D xMHD wrapper object for the Python API
!---------------------------------------------------------------------------------
TYPE :: xmhd2d_instance
  TYPE(oft_xmhd_2d_sim), POINTER :: sim => NULL() !< Simulation object
  CLASS(oft_vector), POINTER :: u => NULL() !< Work vector (state)
  CLASS(oft_vector), POINTER :: v => NULL() !< Work vector (right-hand side)
  TYPE(oft_native_gmres_solver), POINTER :: solver => NULL() !< Linear solver for Newton steps
  TYPE(oft_nksolver), POINTER :: nksolver => NULL() !< Newton solver
END TYPE xmhd2d_instance
!> Simulation the Jacobian update callback acts on (Newton's callback takes only a vector)
TYPE(oft_xmhd_2d_sim), POINTER :: active_sim => NULL()
CONTAINS
!---------------------------------------------------------------------------------
!> Cast a C pointer back to the Fortran wrapper object
!---------------------------------------------------------------------------------
FUNCTION xmhd2d_ccast(xmhd2d_cptr,obj,error_str) RESULT(success)
TYPE(c_ptr), INTENT(in) :: xmhd2d_cptr !< C pointer to wrapper object
TYPE(xmhd2d_instance), POINTER, INTENT(out) :: obj !< Fortran wrapper object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
LOGICAL :: success
CALL copy_string('',error_str)
IF(.NOT.c_associated(xmhd2d_cptr))THEN
  CALL copy_string('xMHD 2D object not associated',error_str)
  success=.FALSE.
  RETURN
END IF
CALL c_f_pointer(xmhd2d_cptr,obj)
success=.TRUE.
END FUNCTION xmhd2d_ccast
!---------------------------------------------------------------------------------
!> Jacobian update callback for the Newton solver
!---------------------------------------------------------------------------------
SUBROUTINE wrapper_update_jacobian(uin)
CLASS(oft_vector), TARGET, INTENT(inout) :: uin !< Current solution
CALL build_approx_jacobian(active_sim,uin)
END SUBROUTINE wrapper_update_jacobian
!---------------------------------------------------------------------------------
!> Build the simulation on a surface mesh
!!
!! The mesh comes from `oft_setup_smesh`. If `xml_ptr` is associated it becomes the
!! global XML root first, so an `<xmhd2d><pre>` preconditioner definition is used
!! exactly as the `harris_sheet` driver would use it.
!---------------------------------------------------------------------------------
SUBROUTINE xmhd2d_setup(xmhd2d_ptr,mesh_ptr,order,xml_ptr,nlocal,nglobal,xml_pre,error_str) &
  BIND(C,NAME="xmhd2d_setup")
TYPE(c_ptr), INTENT(out) :: xmhd2d_ptr !< Pointer to new wrapper object
TYPE(c_ptr), VALUE, INTENT(in) :: mesh_ptr !< Pointer to `multigrid_mesh`
INTEGER(c_int), VALUE, INTENT(in) :: order !< Lagrange element order
TYPE(c_ptr), VALUE, INTENT(in) :: xml_ptr !< XML root node (may be NULL)
INTEGER(c_int), INTENT(out) :: nlocal !< Local size of the state vector
INTEGER(c_int), INTENT(out) :: nglobal !< Global size of the state vector
LOGICAL(c_bool), INTENT(out) :: xml_pre !< An `<xmhd2d><pre>` preconditioner was found
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(xmhd2d_instance), POINTER :: obj
TYPE(multigrid_mesh), POINTER :: mg_mesh
CALL copy_string('',error_str)
IF(ASSOCIATED(active_sim))THEN
  CALL copy_string('Only one xMHD 2D simulation can exist per process',error_str)
  RETURN
END IF
IF(oft_env%nprocs>1)THEN
  CALL copy_string('The xMHD 2D wrapper is serial only',error_str)
  RETURN
END IF
IF(.NOT.c_associated(mesh_ptr))THEN
  CALL copy_string('Mesh object not associated',error_str)
  RETURN
END IF
CALL c_f_pointer(mesh_ptr,mg_mesh)
IF(c_associated(xml_ptr))CALL c_f_pointer(xml_ptr,oft_env%xml)
ALLOCATE(obj)
ALLOCATE(obj%sim)
CALL obj%sim%setup(mg_mesh,order)
CALL obj%sim%fe_rep%vec_create(obj%u)
CALL obj%sim%fe_rep%vec_create(obj%v)
active_sim=>obj%sim
nlocal=obj%sim%u%n
nglobal=INT(obj%sim%u%ng,4)
xml_pre=obj%sim%xml_pre_def%associated()
xmhd2d_ptr=C_LOC(obj)
END SUBROUTINE xmhd2d_setup
!---------------------------------------------------------------------------------
!> Set physical parameters (uniform over regions) and solver tolerances
!---------------------------------------------------------------------------------
SUBROUTINE xmhd2d_set_params(xmhd2d_ptr,chi,eta,nu,gamma,d_diff,den_scale,b0, &
  lin_tol,nl_tol,error_str) BIND(C,NAME="xmhd2d_set_params")
TYPE(c_ptr), VALUE, INTENT(in) :: xmhd2d_ptr !< Pointer to wrapper object
REAL(c_double), VALUE, INTENT(in) :: chi !< Thermal diffusivity [m^2/s]
REAL(c_double), VALUE, INTENT(in) :: eta !< Magnetic diffusivity [m^2/s]
REAL(c_double), VALUE, INTENT(in) :: nu !< Kinematic viscosity [m^2/s]
REAL(c_double), VALUE, INTENT(in) :: gamma !< Adiabatic index
REAL(c_double), VALUE, INTENT(in) :: d_diff !< Density diffusivity [m^2/s]
REAL(c_double), VALUE, INTENT(in) :: den_scale !< Density normalisation [m^-3]
REAL(c_double), INTENT(in) :: b0(3) !< Background field [T]
REAL(c_double), VALUE, INTENT(in) :: lin_tol !< Linear solver tolerance (<=0 keeps default)
REAL(c_double), VALUE, INTENT(in) :: nl_tol !< Newton tolerance (<=0 keeps default)
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(xmhd2d_instance), POINTER :: obj
IF(.NOT.xmhd2d_ccast(xmhd2d_ptr,obj,error_str))RETURN
obj%sim%chi=chi
obj%sim%eta=eta
obj%sim%nu=nu
obj%sim%gamma=gamma
obj%sim%D_diff=d_diff
obj%sim%den_scale=den_scale
obj%sim%B_0=b0
IF(lin_tol>0.d0)obj%sim%lin_tol=lin_tol
IF(nl_tol>0.d0)obj%sim%nl_tol=nl_tol
END SUBROUTINE xmhd2d_set_params
!---------------------------------------------------------------------------------
!> Global index (zero based) of every local entry of the state vector
!!
!! The same numbering `mat_save` writes as `lg`, so global vectors and matrices
!! line up with the exports `nonmodal` reads.
!---------------------------------------------------------------------------------
SUBROUTINE xmhd2d_get_lg(xmhd2d_ptr,lg_ptr,error_str) BIND(C,NAME="xmhd2d_get_lg")
TYPE(c_ptr), VALUE, INTENT(in) :: xmhd2d_ptr !< Pointer to wrapper object
TYPE(c_ptr), VALUE, INTENT(in) :: lg_ptr !< Output `(nlocal)` int32 array
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(xmhd2d_instance), POINTER :: obj
CLASS(oft_native_matrix), POINTER :: jac
INTEGER(c_int), POINTER, DIMENSION(:) :: lg
INTEGER(i4) :: i,j
IF(.NOT.xmhd2d_ccast(xmhd2d_ptr,obj,error_str))RETURN
IF(.NOT.native_matrix_cast(jac,obj%sim%jacobian))THEN
  CALL copy_string('Jacobian is not a native matrix',error_str)
  RETURN
END IF
CALL c_f_pointer(lg_ptr,lg,[jac%nr])
DO i=1,jac%ni
  DO j=1,jac%i_map(i)%n
    lg(jac%i_map(i)%offset+j)=INT(jac%i_map(i)%offsetg+jac%i_map(i)%lge(j),4)-1
  END DO
END DO
END SUBROUTINE xmhd2d_get_lg
!---------------------------------------------------------------------------------
!> Copy the simulation's current state into `vals`, or `vals` into it
!---------------------------------------------------------------------------------
SUBROUTINE xmhd2d_state(xmhd2d_ptr,vals_ptr,set,error_str) BIND(C,NAME="xmhd2d_state")
TYPE(c_ptr), VALUE, INTENT(in) :: xmhd2d_ptr !< Pointer to wrapper object
TYPE(c_ptr), VALUE, INTENT(in) :: vals_ptr !< `(nlocal)` float64 array
LOGICAL(c_bool), VALUE, INTENT(in) :: set !< Copy `vals` into the state instead
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(xmhd2d_instance), POINTER :: obj
REAL(r8), POINTER, DIMENSION(:) :: vals
IF(.NOT.xmhd2d_ccast(xmhd2d_ptr,obj,error_str))RETURN
CALL c_f_pointer(vals_ptr,vals,[obj%sim%u%n])
IF(set)THEN
  CALL obj%sim%u%restore_local(vals)
ELSE
  CALL obj%sim%u%get_local(vals)
END IF
END SUBROUTINE xmhd2d_state
!---------------------------------------------------------------------------------
!> Read a restart file into `vals`, or write `vals` to one
!---------------------------------------------------------------------------------
SUBROUTINE xmhd2d_restart(xmhd2d_ptr,filename,vals_ptr,t,dt,save,error_str) &
  BIND(C,NAME="xmhd2d_restart")
TYPE(c_ptr), VALUE, INTENT(in) :: xmhd2d_ptr !< Pointer to wrapper object
CHARACTER(KIND=c_char), INTENT(in) :: filename(OFT_PATH_SLEN) !< Restart file
TYPE(c_ptr), VALUE, INTENT(in) :: vals_ptr !< `(nlocal)` float64 array
REAL(c_double), VALUE, INTENT(in) :: t !< Time written with the state (save only)
REAL(c_double), VALUE, INTENT(in) :: dt !< Timestep written with the state (save only)
LOGICAL(c_bool), VALUE, INTENT(in) :: save !< Write instead of read
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(xmhd2d_instance), POINTER :: obj
REAL(r8), POINTER, DIMENSION(:) :: vals
CHARACTER(LEN=OFT_PATH_SLEN) :: fname
IF(.NOT.xmhd2d_ccast(xmhd2d_ptr,obj,error_str))RETURN
CALL copy_string_rev(filename,fname)
CALL c_f_pointer(vals_ptr,vals,[obj%sim%u%n])
IF(save)THEN
  CALL obj%u%restore_local(vals)
  CALL obj%sim%rst_save(obj%u,t,dt,TRIM(fname),'U')
ELSE
  CALL obj%sim%rst_load(obj%u,TRIM(fname),'U')
  CALL obj%u%get_local(vals)
END IF
END SUBROUTINE xmhd2d_restart
!---------------------------------------------------------------------------------
!> Evaluate the time-advance function `F_dt(u)` (`dt=0` gives the mass term)
!!
!! `diag` receives OFT's diagnostics for `u`: magnetic energy [J*2*mu0], kinetic
!! energy [J], integrated density, integrated temperature and area.
!---------------------------------------------------------------------------------
SUBROUTINE xmhd2d_residual(xmhd2d_ptr,u_ptr,dt,res_ptr,diag,error_str) &
  BIND(C,NAME="xmhd2d_residual")
TYPE(c_ptr), VALUE, INTENT(in) :: xmhd2d_ptr !< Pointer to wrapper object
TYPE(c_ptr), VALUE, INTENT(in) :: u_ptr !< `(nlocal)` state
REAL(c_double), VALUE, INTENT(in) :: dt !< Timestep in `F_dt`
TYPE(c_ptr), VALUE, INTENT(in) :: res_ptr !< `(nlocal)` result
REAL(c_double), INTENT(out) :: diag(5) !< Diagnostics of `u`
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(xmhd2d_instance), POINTER :: obj
REAL(r8), POINTER, DIMENSION(:) :: vals
IF(.NOT.xmhd2d_ccast(xmhd2d_ptr,obj,error_str))RETURN
CALL c_f_pointer(u_ptr,vals,[obj%sim%u%n])
CALL obj%u%restore_local(vals)
obj%sim%nlfun%dt=dt
CALL obj%sim%nlfun%apply(obj%u,obj%v)
CALL c_f_pointer(res_ptr,vals,[obj%sim%u%n])
CALL obj%v%get_local(vals)
diag=obj%sim%nlfun%diag_vals
END SUBROUTINE xmhd2d_residual
!---------------------------------------------------------------------------------
!> Assemble the Jacobian of `F_dt` at `u`; returns its local CSR sizes
!!
!! The matrix stays in the simulation object until the next assembly or advance;
!! fetch it with @ref xmhd2d_jacobian_get. Rows are local (periodic seam copies
!! hold partial rows), so the global operator is `P^T J P` with `P` from `lg`.
!---------------------------------------------------------------------------------
SUBROUTINE xmhd2d_jacobian(xmhd2d_ptr,u_ptr,dt,nr,nnz,error_str) BIND(C,NAME="xmhd2d_jacobian")
TYPE(c_ptr), VALUE, INTENT(in) :: xmhd2d_ptr !< Pointer to wrapper object
TYPE(c_ptr), VALUE, INTENT(in) :: u_ptr !< `(nlocal)` state
REAL(c_double), VALUE, INTENT(in) :: dt !< Timestep in `F_dt`
INTEGER(c_int), INTENT(out) :: nr !< Number of local rows
INTEGER(c_int), INTENT(out) :: nnz !< Number of stored entries
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(xmhd2d_instance), POINTER :: obj
CLASS(oft_native_matrix), POINTER :: jac
REAL(r8), POINTER, DIMENSION(:) :: vals
IF(.NOT.xmhd2d_ccast(xmhd2d_ptr,obj,error_str))RETURN
CALL c_f_pointer(u_ptr,vals,[obj%sim%u%n])
CALL obj%u%restore_local(vals)
obj%sim%jac_dt=dt
CALL build_approx_jacobian(obj%sim,obj%u)
IF(.NOT.native_matrix_cast(jac,obj%sim%jacobian))THEN
  CALL copy_string('Jacobian is not a native matrix',error_str)
  RETURN
END IF
nr=jac%nr
nnz=jac%nnz
END SUBROUTINE xmhd2d_jacobian
!---------------------------------------------------------------------------------
!> Copy the last assembled Jacobian out in zero-based CSR form
!---------------------------------------------------------------------------------
SUBROUTINE xmhd2d_jacobian_get(xmhd2d_ptr,kr_ptr,lc_ptr,m_ptr,error_str) &
  BIND(C,NAME="xmhd2d_jacobian_get")
TYPE(c_ptr), VALUE, INTENT(in) :: xmhd2d_ptr !< Pointer to wrapper object
TYPE(c_ptr), VALUE, INTENT(in) :: kr_ptr !< `(nr+1)` int32 row pointers
TYPE(c_ptr), VALUE, INTENT(in) :: lc_ptr !< `(nnz)` int32 column indices
TYPE(c_ptr), VALUE, INTENT(in) :: m_ptr !< `(nnz)` float64 values
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(xmhd2d_instance), POINTER :: obj
CLASS(oft_native_matrix), POINTER :: jac
INTEGER(c_int), POINTER, DIMENSION(:) :: kr,lc
REAL(c_double), POINTER, DIMENSION(:) :: mvals
IF(.NOT.xmhd2d_ccast(xmhd2d_ptr,obj,error_str))RETURN
IF(.NOT.native_matrix_cast(jac,obj%sim%jacobian))THEN
  CALL copy_string('Jacobian is not a native matrix',error_str)
  RETURN
END IF
CALL c_f_pointer(kr_ptr,kr,[jac%nr+1])
CALL c_f_pointer(lc_ptr,lc,[jac%nnz])
CALL c_f_pointer(m_ptr,mvals,[jac%nnz])
kr=jac%kr-1
lc=jac%lc-1
mvals=jac%M
END SUBROUTINE xmhd2d_jacobian_get
!---------------------------------------------------------------------------------
!> Set up the Newton and linear solvers used by @ref xmhd2d_advance
!!
!! Mirrors `run_simulation`: GMRES on the assembled Jacobian, preconditioned as
!! the XML says (diagonal scaling otherwise), Newton refreshing the Jacobian every
!! 4 iterations on top of the per-step refresh.
!---------------------------------------------------------------------------------
SUBROUTINE setup_solvers(obj)
TYPE(xmhd2d_instance), INTENT(inout) :: obj
IF(ASSOCIATED(obj%nksolver))RETURN
ALLOCATE(obj%solver,obj%nksolver)
obj%solver%A=>obj%sim%jacobian
obj%solver%its=obj%sim%ittarget*4
obj%solver%atol=obj%sim%lin_tol
obj%solver%itplot=1
obj%solver%nrits=40
obj%solver%pm=.FALSE.
NULLIFY(obj%solver%pre)
IF(obj%sim%xml_pre_def%associated())THEN
  CALL create_solver_xml(obj%solver%pre,obj%sim%xml_pre_def)
ELSE
  CALL create_diag_pre(obj%solver%pre)
END IF
obj%solver%pre%A=>obj%sim%jacobian
obj%nksolver%A=>obj%sim%nlfun
obj%nksolver%J_inv=>obj%solver
obj%nksolver%its=30
obj%nksolver%atol=obj%sim%nl_tol
obj%nksolver%rtol=1.d-20
obj%nksolver%J_update=>wrapper_update_jacobian
obj%nksolver%up_freq=4
END SUBROUTINE setup_solvers
!---------------------------------------------------------------------------------
!> Advance `u` by `nsteps` backward-Euler steps of fixed size `dt`
!!
!! Each step solves `F_dt(u_{n+1}) = F_0(u_n)` with OFT's Newton solver, starting
!! from a quadratic extrapolation of the previous states as `run_simulation` does.
!! Unlike `run_simulation` the timestep never adapts and nothing is written to
!! disk. If `traj_ptr` is associated it receives all `nsteps+1` states, starting
!! with the initial one. `diag` receives, per state, magnetic and kinetic energy,
!! and `iters` the Newton and linear iteration counts of each step. On a failed
!! Newton solve the advance stops, `nfailed` is the failing step (1-based), and
!! `u` holds the last converged state.
!---------------------------------------------------------------------------------
SUBROUTINE xmhd2d_advance(xmhd2d_ptr,u_ptr,dt,nsteps,traj_ptr,diag_ptr,iters_ptr,verbose,nfailed,error_str) &
  BIND(C,NAME="xmhd2d_advance")
TYPE(c_ptr), VALUE, INTENT(in) :: xmhd2d_ptr !< Pointer to wrapper object
TYPE(c_ptr), VALUE, INTENT(in) :: u_ptr !< `(nlocal)` state, advanced in place
REAL(c_double), VALUE, INTENT(in) :: dt !< Timestep [s]
INTEGER(c_int), VALUE, INTENT(in) :: nsteps !< Number of steps
TYPE(c_ptr), VALUE, INTENT(in) :: traj_ptr !< `(nsteps+1,nlocal)` states (may be NULL)
TYPE(c_ptr), VALUE, INTENT(in) :: diag_ptr !< `(nsteps+1,2)` magnetic and kinetic energy
TYPE(c_ptr), VALUE, INTENT(in) :: iters_ptr !< `(nsteps,2)` Newton and linear iterations
LOGICAL(c_bool), VALUE, INTENT(in) :: verbose !< Print solver progress
INTEGER(c_int), INTENT(out) :: nfailed !< 0 on success, else the failing step
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
TYPE(xmhd2d_instance), POINTER :: obj
REAL(r8), POINTER, DIMENSION(:) :: vals
REAL(r8), POINTER, DIMENSION(:,:) :: traj,diag
INTEGER(c_int), POINTER, DIMENSION(:,:) :: iters
INTEGER(i4), PARAMETER :: maxextrap=2
INTEGER(i4) :: i,j,nextrap,n
REAL(r8) :: t,extrapt(maxextrap)
LOGICAL :: pm_save
REAL(r8), POINTER, DIMENSION(:) :: tmp
TYPE(oft_vector_ptr) :: extrap_fields(maxextrap)
IF(.NOT.xmhd2d_ccast(xmhd2d_ptr,obj,error_str))RETURN
nfailed=0
n=obj%sim%u%n
CALL c_f_pointer(u_ptr,vals,[n])
CALL c_f_pointer(diag_ptr,diag,[2,nsteps+1])
CALL c_f_pointer(iters_ptr,iters,[2,nsteps])
NULLIFY(traj)
IF(c_associated(traj_ptr))CALL c_f_pointer(traj_ptr,traj,[n,nsteps+1])
active_sim=>obj%sim
pm_save=oft_env%pm
oft_env%pm=verbose
obj%sim%dt=dt
obj%sim%jac_dt=dt
CALL setup_solvers(obj)
obj%solver%pm=verbose
DO i=1,maxextrap
  CALL obj%sim%fe_rep%vec_create(extrap_fields(i)%f)
END DO
extrapt=0.d0
nextrap=0
t=0.d0
CALL obj%u%restore_local(vals)
IF(ASSOCIATED(traj))traj(:,1)=vals
DO i=1,nsteps
  !---Right-hand side F_0(u_n): the mass term
  obj%sim%nlfun%dt=0.d0
  CALL obj%sim%nlfun%apply(obj%u,obj%v)
  diag(:,i)=obj%sim%nlfun%diag_vals(1:2)
  obj%sim%nlfun%dt=dt
  CALL build_approx_jacobian(obj%sim,obj%u)
  CALL obj%solver%pre%update(.TRUE.)
  !---Initial guess from the previous states
  DO j=maxextrap,2,-1
    CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
    extrapt(j)=extrapt(j-1)
  END DO
  IF(nextrap<maxextrap)nextrap=nextrap+1
  CALL extrap_fields(1)%f%add(0.d0,1.d0,obj%u)
  extrapt(1)=t
  IF(i>maxextrap)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,t+dt,obj%u)
  CALL obj%nksolver%apply(obj%u,obj%v)
  IF(obj%nksolver%cits<0)THEN
    nfailed=i
    CALL extrap_fields(1)%f%get_local(vals) ! Last converged state
    EXIT
  END IF
  t=t+dt
  iters(:,i)=[obj%nksolver%nlits,obj%nksolver%lits]
  IF(ASSOCIATED(traj))THEN
    tmp=>traj(:,i+1)
    CALL obj%u%get_local(tmp)
  END IF
END DO
IF(nfailed==0)THEN
  !---Diagnostics of the final state
  obj%sim%nlfun%dt=0.d0
  CALL obj%sim%nlfun%apply(obj%u,obj%v)
  diag(:,nsteps+1)=obj%sim%nlfun%diag_vals(1:2)
  CALL obj%u%get_local(vals)
END IF
DO i=1,maxextrap
  CALL extrap_fields(i)%f%delete
  DEALLOCATE(extrap_fields(i)%f)
END DO
oft_env%pm=pm_save
END SUBROUTINE xmhd2d_advance
END MODULE xmhd_2d_f
