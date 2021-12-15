module hf
   use const

   implicit none

   contains
      subroutine do_hartree_fock(sys, int_store)
         use error_handling
         use system
         use read_in_integrals
         use linalg

         type(int_store_t), intent(in) :: int_store
         type(system_t), intent(in) :: sys

         type(mat_t) :: ovlp
         real(p) :: ortmat(sys%nbasis,sys%nbasis)
         real(p) :: fock(sys%nbasis,sys%nbasis)
         real(p) :: coeff(sys%nbasis,sys%nbasis)
         real(p) :: density(sys%nbasis,sys%nbasis)
         real(p) :: tmpmat(sys%nbasis,sys%nbasis)

         ! Initialise the overlap matrix stored in int_store_t into a mat_t object
         call init_mat_t(sys, int_store%ovlp, ovlp)

         ! S^(-1/2) = C * L^(-1/2) * C^(T)
         call orthogonalise(ovlp)
         call build_ortmat(ovlp, ortmat)

         ! Initial guess Fock matrix
         ! F_0 = S^T(-1/2) * H_core * S^(-1/2)
         fock = 0.0_p
         fock = matmul(int_store%core_hamil, ortmat)
         fock = matmul(transpose(ortmat), fock)

         call init_mat_


      end subroutine do_hartree_fock


      subroutine orthogonalise(mat)
         ! For a symmetric matrix A we produce
         ! A = CLC^(-1)
         ! With L and C stored in a mat_t data structure

         use linalg
         use error_handling

         type(mat_t), intent(inout) :: mat

         integer :: ierr

         call eigs(mat, ierr)
         if (ierr /= 0) call error('hf::orthogonalise', 'Orthogonalisation failed!')

      end subroutine orthogonalise

      subroutine init_mat_t(sys, matrix, mat)
         use read_in_integrals
         use linalg, only: mat_t
         use system
        
         type(system_t), intent(in) :: sys
         real(p), intent(in) :: matrix(:,:)
         type(mat_t), intent(out) :: mat

         mat%N = sys%nbasis
         allocate(mat%A(mat%N, mat%N))
         mat%A(:,:) = matrix(:,:)
         allocate(mat%W(sys%nbasis))

      end subroutine init_mat_t

      subroutine build_ortmat(ovlp, ortmat)

         use linalg, only: mat_t

         type(mat_t), intent(in) :: ovlp
         real(p), allocatable, intent(out) :: ortmat(:,:)

         real(p), allocatable :: tmpmat(:,:)

         integer :: i
         
         associate(n=>ovlp%N, diag=>ovlp%W)
            allocate(ortmat(n, n), source = 0.0_p)
            allocate(tmpmat(n, n), source = 0.0_p)
            do i = 1, n
               ortmat(i,i) = 1/sqrt(diag(i))
            end do
         end associate
         
         tmpmat = matmul(ortmat, transpose(ovlp%A))
         ortmat = matmul(ovlp%A, tmpmat)

         deallocate(tmpmat)

      end subroutine build_ortmat
         
end module hf
