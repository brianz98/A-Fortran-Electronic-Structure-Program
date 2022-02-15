module error_handling
   use, intrinsic :: iso_fortran_env, only : iunit=>output_unit, stderr=>error_unit

   implicit none

   contains
      subroutine error(procedure, error_msg)

         character(*), intent(in) :: procedure, error_msg
         character(3), parameter :: error_str = '999'

         write(stderr, '(1X, A6)') 'ERROR.'
         write(stderr, '(1X, A)') 'Programme stops in procedure: '//adjustl(procedure)//'.'
         write(stderr, '(1X, A)') 'Reason: '//adjustl(error_msg)//'.'
         write(stderr, '(1X, A10)') 'EXITING...'

         stop error_str
         
         return
      end subroutine error

      subroutine check_allocate(array_name, array_size, ierr)
         ! Wrapper around an ierr check, so error messages are more helpful.
         ! In:
         !     array_name: the name of the array being allocated
         !     array_size: the intended size of the array
         !     ierr: error code returned from allocate()
         character(*), intent(in) :: array_name
         integer, intent(in) :: array_size
         integer, intent(in) :: ierr

         if (array_size < 0 .or. ierr /= 0) then
            ! If via some previous calculation array_size is negative we intercept it
            write(stderr, '(1X, A, 1X, A, I0)') 'Error in allocating array', trim(array_name), 'with size', array_size
            write(stderr, '(1X, A, 1X, I0)') 'Error code:', ierr
            call error('check_allocate', 'Allocation error')
         end if
      end subroutine check_allocate

end module error_handling
