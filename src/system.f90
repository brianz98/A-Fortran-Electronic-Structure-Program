module system
   use const
   
   implicit none

   type system_t
      integer :: nbasis = 0
   end type system_t

   type state_t
      integer :: iter = 0
      real(p) :: energy = 0.0_p
      real(p) :: energy_old = 0.0_p
   end type state_t

end module system