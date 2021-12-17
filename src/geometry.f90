module geometry
    use const

    implicit none

    contains

        subroutine read_geometry_in(sys, int_store)
            ! Read in geometry data supplied in a file
            use system, only: system_t
            use read_in_integrals, only: int_store_t

            type(system_t), intent(inout) :: sys
            type(int_store_t), intent(inout) :: int_store

            integer :: i, iunit, ir
            real(p) :: charge

            ! stdout unit is 6
            iunit = 6

            ! newunit chooses an unused file handle, 'old' means existing file
            open(newunit=ir, file='dat/geom.dat', status='old', form='formatted')
            read(ir, *) sys%natoms
            
            ! Lazy allocation! 
            ! [TODO] make check_allocate subroutine
            allocate(sys%charges(sys%natoms))
            allocate(sys%coords(3,sys%natoms))
            
            do i = 1, sys%natoms
                read(ir, *) charge, sys%coords(:,i)
                sys%charges(i) = int(charge)
            end do
            
            close(ir)

            sys%nel = sum(sys%charges)

            call get_e_nuc(sys, int_store%e_nuc)


        end subroutine read_geometry_in

        subroutine get_bond_lengths(coords, bondlengths)
            real(p), intent(in) :: coords(:,:)
            real(p), allocatable, intent(out) :: bondlengths(:,:)

            integer :: i, j

            ! Negligible memeory cost to store the full matrix rather than the lower triangular
            associate(natoms=>size(coords, dim=1))
                allocate(bondlengths(natoms, natoms), source=0.0_p)

                do i = 1, natoms-1
                    do j = i+1, natoms
                        ! Diagonals set to 0.0 in allocation
                        bondlengths(j, i) = norm2(coords(:,i)-coords(:,j))
                        bondlengths(i, j) = bondlengths(j, i)
                    end do
                end do
            end associate
        end subroutine get_bond_lengths

        subroutine get_e_nuc(sys, e_nuc)
            use system, only: system_t

            type(system_t), intent(in) :: sys
            real(p), intent(out) :: e_nuc
            real(p), allocatable :: bondlengths(:,:)

            integer :: i, j

            call get_bond_lengths(sys%coords, bondlengths)

            e_nuc = 0.0_p
            do j = 2, size(sys%charges)
                do i = 1, j-1
                    e_nuc = e_nuc + sys%charges(i)*sys%charges(j)/bondlengths(i,j)
                end do
            end do

            ! Explicitly deallocate instead of using scoping rules.. is it necessary?
            deallocate(bondlengths)

        end subroutine get_e_nuc


            
end module geometry
