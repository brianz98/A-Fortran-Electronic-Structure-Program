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

         call init_mat_t(sys, int_store%ovlp, ovlp)

         call orthogonalise(ovlp)

         print*, ovlp%W
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
         use linalg
         use system
        
         type(system_t), intent(in) :: sys
         real(p), intent(in) :: matrix(:,:)
         type(mat_t), intent(out) :: mat

         mat%N = sys%nbasis
         allocate(mat%A(mat%N, mat%N))
         mat%A(:,:) = matrix(:,:)
         allocate(mat%W(sys%nbasis))
      end subroutine init_mat_t

end module hf
