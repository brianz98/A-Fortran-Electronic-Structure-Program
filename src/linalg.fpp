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

      subroutine dgemm_wrapper(transA, transB, outer_row, outer_col, inner_dim, A, B, C, alpha, beta)
         ! Wraps around dgemm with fewer arguments

         character(1), intent(in) :: transA, transB
         integer, intent(in) :: outer_row, outer_col, inner_dim
         real(p), intent(in) :: A(*), B(*)
         real(p), optional, intent(in) :: alpha, beta
         real(p), intent(inout) :: C(*)
         integer :: LDA, LDB
         real(p) :: alpha_loc, beta_loc

         if (transA == 'T') then
            LDA = inner_dim
         else
            LDA = outer_row
         end if

         if (transB == 'T') then
            LDB = outer_col
         else
            LDB = inner_dim
         end if

         alpha_loc = 1.0_p
         if (present(alpha)) alpha_loc = alpha

         beta_loc = 0.0_p
         if (present(beta)) beta_loc = beta
         
         call dgemm(transA, transB, outer_row, outer_col, inner_dim, alpha_loc, &
                    A, LDA, B, LDB, beta_loc, C, outer_row)
      end subroutine dgemm_wrapper

      elemental subroutine zero_mat(matel)
         ! Zero out entries smaller than machine precision

         real(p), intent(inout) :: matel

         if (abs(matel)<depsilon) matel = 0.0_p
      end subroutine zero_mat

      subroutine omp_reshape(out_arr, in_arr, arr_order)
         ! Performs an equivalent operation to the Fortran intrinsic `reshape' but with OMP acceleration.
         ! The `dimension(:,:,:,:)' specification also means array lower bound must be 1 (default Fortran)

         ! The limitation (perhaps of my understanding) of the Fypp preprocessor means that 
         ! we can't achieve the Fortran `reshape' style permutation, and instead we can only offer the more 'common-sense'
         ! version of i,j,k,l (2,3,1,4)-> j,k,i,l
         
         ! [todo] - Checks on arr_order

         real(dp), dimension(:,:,:,:), intent(in) :: in_arr
         character(4), intent(in) :: arr_order
         real(dp), dimension(:,:,:,:), intent(out) :: out_arr

         integer :: i, j, k, l, iu, ju, ku, lu

         iu = ubound(in_arr, dim=1)
         ju = ubound(in_arr, dim=2)
         ku = ubound(in_arr, dim=3)
         lu = ubound(in_arr, dim=4)

#:set ORDER_STRS = [''.join(_) for _ in list(itertools.permutations('1234'))]
#:set INDEX_LOOKUP = ['i','j','k','l']
#!
#:def permute_index(ORD_STR)
$:','.join([INDEX_LOOKUP[int(_)-1] for _ in ORD_STR])
#:enddef permute_index
#!
#:for ORDER_STR in ORDER_STRS
         if (arr_order == '${ORDER_STR}$') then
            !$omp parallel do default(none) shared(in_arr, out_arr, iu, ju, ku, lu)
            do l = 1, lu
               do k = 1, ku
                  do j = 1, ju
                     do i = 1, iu
                        out_arr(${permute_index(ORDER_STR)}$) = in_arr(i,j,k,l)
                     end do
                  end do
               end do
            end do
            !$omp end parallel do
         end if
#:endfor
            
      end subroutine omp_reshape

      
end module linalg
