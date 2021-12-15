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

         type(mat_t) :: ovlp, fockmat
         real(p) :: ortmat(sys%nbasis,sys%nbasis)
         real(p) :: fock(sys%nbasis,sys%nbasis)
         real(p) :: coeff(sys%nbasis,sys%nbasis)
         real(p) :: density(sys%nbasis,sys%nbasis)
         real(p) :: tmpmat(sys%nbasis,sys%nbasis)

         ! Initialise the overlap matrix stored in int_store_t into a mat_t object
         call init_mat_t(int_store%ovlp, ovlp)

         ! S^(-1/2) = C * L^(-1/2) * C^(T)
         call diagonalise(ovlp)
         call build_ortmat(ovlp, ortmat)

         ! Initial guess Fock matrix
         ! F_0 = S^T(-1/2) * H_core * S^(-1/2)
         fock = 0.0_p
         fock = matmul(int_store%core_hamil, ortmat)
         fock = matmul(transpose(ortmat), fock)

         call init_mat_t(fock, fockmat)
         call diagonalise(fockmat)

         ! Transform coefficient into AO basis
         fockmat%A = matmul(ortmat, fockmat%A)

         ! Build initial density
         ! [TODO]: this is just hard-coded for now for water, URGENT: add geom_read_in
         call build_density(fockmat%A, density, 5)

         print*, density

      end subroutine do_hartree_fock


      subroutine diagonalise(mat)
         ! For a symmetric matrix A we produce
         ! A = C * L * C^T
         ! With L and C stored in a mat_t data structure

         use linalg
         use error_handling

         type(mat_t), intent(inout) :: mat

         integer :: ierr

         call eigs(mat, ierr)
         if (ierr /= 0) call error('hf::orthogonalise', 'Orthogonalisation failed!')

      end subroutine diagonalise

      subroutine init_mat_t(matrix, mat)
         use read_in_integrals
         use linalg, only: mat_t
         use system
        
         real(p), intent(in) :: matrix(:,:)
         type(mat_t), intent(out) :: mat

         associate(N=>size(matrix, dim=1))
            mat%N = N
            allocate(mat%A(N,N))
            mat%A(:,:) = matrix(:,:)
            allocate(mat%W(N))
         end associate

      end subroutine init_mat_t

      subroutine build_ortmat(ovlp, ortmat)

         use linalg, only: mat_t

         type(mat_t), intent(in) :: ovlp
         real(p), intent(out) :: ortmat(:,:)

         real(p) :: tmpmat(size(ortmat, dim=1), size(ortmat, dim=1))

         integer :: i
         
         tmpmat(:,:) = 0.0_p
         ortmat(:,:) = 0.0_p

         associate(diag=>ovlp%W)
            do i = 1, size(ortmat, dim=1)
               ortmat(i,i) = 1/sqrt(diag(i))
            end do
         end associate
         
         tmpmat = matmul(ortmat, transpose(ovlp%A))
         ortmat = matmul(ovlp%A, tmpmat)

      end subroutine build_ortmat

      subroutine build_density(coeff, density, nocc)
         real(p), intent(in) :: coeff(:,:)
         integer, intent(in) :: nocc
         real(p), intent(out) :: density(:,:)

         density = matmul(coeff(:, 1:nocc), transpose(coeff(:, 1:nocc)))

      end subroutine build_density


         
end module hf
