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

      subroutine dgemm_wrapper(transA, transB, outer_row, outer_col, inner_dim, A, B, C)
         ! Wraps around dgemm with fewer arguments

         character(1), intent(in) :: transA, transB
         integer, intent(in) :: outer_row, outer_col, inner_dim
         real(p), intent(in) :: A(*), B(*)
         real(p), intent(inout) :: C(*)

         call dgemm(transA, transB, outer_row, outer_col, inner_dim, 1.0_dp, &
                    A, outer_row, B, inner_dim, 0.0_dp, C, outer_row)
      end subroutine dgemm_wrapper

      subroutine tensor_dot_product_2_4(A, B, C)
         ! Computes tensor contractions of the form C_i^a = A_e^m B_mi^ea
         real(p), intent(in) :: A(:,:), B(:,:,:,:)
         real(p), intent(inout) :: C(:,:)
         integer :: cu, cl, au, al, bu, bl
         integer :: cll, clu, cul, cuu, all, alu, aul, auu

         cll = lbound(C, dim=1); clu = ubound(C, dim=1); cul = lbound(C, dim=2); cuu = ubound(C, dim=2)
         all = lbound(A, dim=1); alu = ubound(A, dim=1); aul = lbound(A, dim=2); auu = ubound(A, dim=2)

         !$omp parallel default(none)&
         do cu = cul, cuu
            do cl = cll, clu
               tmp = 0.0_p
               do au = aul, auu
                  do al = all, alu
                     tmp = tmp + A(al, au)*B(au,cl,)

      end subroutine

      elemental subroutine zero_mat(matel)
         ! Zero out entries smaller than machine precision

         real(p), intent(inout) :: matel

         if (abs(matel)<depsilon) matel = 0.0_p
      end subroutine zero_mat
      
end module linalg
