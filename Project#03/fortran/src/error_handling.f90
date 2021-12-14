module error_handling
   implicit none

   contains
      subroutine error(procedure, error_msg)
         use, intrinsic :: iso_fortran_env, only: error_unit, output_unit

         character(*), intent(in) :: procedure, error_msg
         character(3), parameter :: error_str = '999'

         write(error_unit, '(1X, A6)') 'ERROR.'
         write(error_unit, '(1X, A)') 'Programme stops in procedure: '//adjustl(procedure)//'.'
         write(error_unit, '(1X, A)') 'Reason: '//adjustl(error_msg)//'.'
         write(error_unit, '(1X, A10)') 'EXITING...'

         stop error_str
         
         return
      end subroutine error

end module error_handling
