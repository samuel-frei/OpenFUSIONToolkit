!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file xmhd_2d.F90
!
!> Solve time-dependent MHD equations with lagrange basis
!---------------------------------------------------------------------------
MODULE xmhd_2d
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, &
    hdf5_create_timestep, hdf5_field_exist, oft_bin_file
USE oft_quadrature
USE oft_mesh_type, ONLY: smesh, cell_is_curved
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat, oft_vector_ptr, &
    vector_extrapolate
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
!
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec
USE oft_lag_basis, ONLY: oft_lag_setup, oft_blagrange, oft_blag_eval, oft_blag_geval
IMPLICIT NONE
#include "local.h"
PRIVATE

TYPE, extends(oft_noop_matrix) :: xmhd_2d_nlfun
  REAL(r8) :: dt = -1.d0 !< Time step
  REAL(r8) :: kappa_e !< Needs docs
  REAL(r8) :: kappa_i !< Needs docs
  REAL(r8) :: tau_eq !< Needs docs
  REAL(r8) :: diag_vals(2) = 0.d0 !< Needs docs
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: Ti_bc => NULL() !< Ti BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: Te_bc => NULL() !< Te BC flag
CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => nlfun_apply
END TYPE xmhd_2d_nlfun

subroutine nlfun_apply(self,a,b)
class(tdiff_nlfun), intent(inout) :: self !< NL function object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
LOGICAL :: curved
INTEGER(i4) :: i,m,jr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: D_dens, temp_gamma, k_boltz, chi_temp,diag_vals(2)  
REAL(r8) :: n,T,u(3), dn(3),dT(3),du(3,3), jac_mat(3,4),jac_det #figure out dimensions on jac_mat
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,n_weights_loc,T_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: u_weights_loc, basis_grads,res_loc
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights,u_weights, n_res,T_res
REAL(r8), POINTER, DIMENSION(:, :) :: u_weights
quad=>oft_blagrange%quad
NULLIFY(n_weights,T_weights,u_weights, n_res,T_res)
!---Get weights from solution vector
CALL a%get_local(n_weights,1)
CALL a%get_local(T_weights,2)
CALL a%get_local(u_weights,3)
!---
D_dens=self%D_dens
temp_gamma=self%temp_gamma
k_boltz=self%k_boltz
chi_temp=self%chi_temp
IF(self%tau_eq>0.d0)THEN
  tau_eq_inv=1.d0/self%tau_eq
ELSE
  tau_eq_inv=0.d0
END IF
!---Zero result and get storage array
CALL b%set(0.d0)
CALL b%get_local(n_res,1)
CALL b%get_local(T_res,2)
CALL b%get_local(u_res,3)
diag_vals=0.d0
!$omp parallel private(m,jr,curved,cell_dofs,basis_vals,basis_grads,n_weights_loc, &
!$omp T_weights_loc,u_weights_loc, res_loc,jac_mat,jac_det, &
!$omp n,T,u, dn,dT, dV) reduction(+:diag_vals)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(n_weights_loc(oft_blagrange%nce),T_weights_loc(oft_blagrange%nce),u_weights_loc(oft_blagrange%nce,3) )
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,2))
!$omp do schedule(static)
DO i=1,smesh%nc
  curved=cell_is_curved(smesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  n_weights_loc = n_weights(cell_dofs)
  T_weights_loc = T_weights(cell_dofs)
  u_weights_loc = u_weights(cell_dofs,3)
!---------------------------------------------------------------------------
! Quadrature Loop
!---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call smesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !---Reconstruct values of solution fields
    n = 0.d0; dn = 0.d0; T = 0.d0; dT = 0.d0; u = 0.d0, du = 0.d0
    DO jr=1,oft_blagrange%nce
      n = n + n_weights_loc(jr)*basis_vals(jr)
      T = T + T_weights_loc(jr)*basis_vals(jr)
      u = u + u_weights_loc(jr)*basis_vals(jr)
      dn = dn + n_weights_loc(jr)*basis_grads(:,jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      du = du + u_weights_loc(:,jr)*basis_grads(:,jr)
    END DO
    diag_vals = diag_vals + [n,T, u]*jac_det*quad%wts(m)
    !---Compute local function contributions
    DO jr=1,oft_blagrange%nce
      res_loc(jr,1) = res_loc(jr,1) &
       + basis_vals(jr)*n*jac_det*quad%wts(m) &
       + self%dt*basis_vals(jr)*DOT_PRODUCT(dn,u) &
       + self%dt*basis_vals(jr)*n*(du(1,1) + du(2,2) + du(3,3)) ##ASK CHRIS HOW TO DO DIVERGENCE
       + self%dt*D_dens*DOT_PRODUCT(basis_grads(:,jr),dn)
      
      res_loc(jr,2) = res_loc(jr,2) &
        + 1/(temp_gamma+1)*basis_vals(jr)*T*jac_det*quad%wts(m) &
        + 1/(temp_gamma+1)*self%dt*basis_vals(jr)*DOT_PRODUCT(u,dT) &
        - self%dt*k_boltz*basis_vals(jr)*T*(du(1,1) + du(2,2) + du(3,3)) &
        + self%dt*chi_temp*DOT_PRODUCT(dT,basis_grads(:,jr)) &
        - self%dt*chi_temp*basis_vals(jr)*DOT_PRODUCT(dT,dn)/n
    END DO
  END DO
  !---Add local values to full vector
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    n_res(cell_dofs(jr)) = n_res(cell_dofs(jr)) + res_loc(jr,1)
    !$omp atomic
    T_res(cell_dofs(jr)) = T_res(cell_dofs(jr)) + res_loc(jr,2)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads,n_weights_loc,T_weights_loc,u_weights_loc,cell_dofs,res_loc)
!$omp end parallel
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,n_weights,n_res,self%n_bc)
CALL fem_dirichlet_vec(oft_blagrange,T_weights,T_res,self%T_bc)
!---Put results into full vector
CALL b%restore_local(n_res,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(T_res,2,add=.TRUE.)
self%diag_vals=oft_mpi_sum(diag_vals,2)
!---Cleanup remaining storage
DEALLOCATE(n_res,T_res,n_weights,T_weights, u_weights)
end subroutine nlfun_apply