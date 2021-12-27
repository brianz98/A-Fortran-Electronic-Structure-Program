module linalg
   use const

   implicit none

   type mat_t
      ! A structure that holds the original matrix, the eigenvalues and the eigenvalue matrix
      integer :: N = -1 ! Dimension
      real(p), allocatable :: ao(:,:) ! AO basis matrix
      real(p), allocatable :: ao_ort(:,:) ! Orthonormal AO basis
      real(p), allocatable :: A(:,:) ! eigenvector matrix after diagonalisation
      real(p), allocatable :: W(:) ! Eigenvalues
   end type mat_t

   contains
      subroutine eigs(mat, ierr)
         ! Wrapper around LAPACK dsyev, so named to be reminiscent of the NumPy package 
         ! Finds the eigenvalues and eigenvectors of a double-precision symmetric matrix
         ! This further automates finding the optimal workspace by first passing lwork = -1
         ! to LAPACK dsyev, where work(1) contains the optimised workspace dimension,
         ! then we pass that to dsyev again for actual computation of eigen quantities

         type(mat_t), intent(inout) :: mat
         integer, intent(out) :: ierr

         real(p), allocatable :: work(:)
         integer :: lwork, i

         lwork = -1
         do i = 1, 2
            allocate(work(abs(lwork)))
            call dsyev('V', 'U', mat%N, mat%A(:,:), mat%N, mat%W, work, lwork, ierr)
            lwork = nint(work(1))
            deallocate(work)
         end do
      end subroutine eigs

      subroutine linsolve(A, B, ierr)
         real(p), intent(inout) :: A(:,:)
         real(p), intent(inout) :: B(:)
         integer, intent(out) :: ierr

         real(p), allocatable :: work(:)
         real(p), allocatable :: ipiv(:)
         integer :: lwork, i

         lwork = -1
         allocate(ipiv(size(B)))
         do i = 1, 2
            allocate(work(abs(lwork)))
            call dsysv('L', size(B), 1, A, size(B), ipiv, B, size(B), work, lwork, ierr)
            lwork = nint(work(1))
            deallocate(work)
         end do

      end subroutine linsolve

      elemental subroutine zero_mat(matel)
         ! Zero out entries smaller than machine precision

         real(p), intent(inout) :: matel

         if (abs(matel)<depsilon) matel = 0.0_p
      end subroutine zero_mat
      
end module linalg
