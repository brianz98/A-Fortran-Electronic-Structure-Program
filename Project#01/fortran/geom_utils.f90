module geom_utils
   use const
   
   implicit none

   contains
      subroutine get_bond_lengths(natoms, coords, bondlengths)
         integer, intent(in) :: natoms
         real(p), intent(in) :: coords(:,:)
         real(p), allocatable, intent(out) :: bondlengths(:,:)

         integer :: i, j, nelements

         ! Negligible memeory cost to store the full matrix rather than the lower triangular
         nelements = natoms*(natoms-1)/2
         allocate(bondlengths(natoms, natoms), source=0.0_p)

         do i = 1, natoms-1
            do j = i+1, natoms
               ! Diagonals set to 0.0 in allocation
               bondlengths(i, j) = norm2(coords(:,i)-coords(:,j))
               bondlengths(j, i) = bondlengths(i, j)
            end do
         end do
      end subroutine get_bond_lengths


end module geom_utils
