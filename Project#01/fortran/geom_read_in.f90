module geom_read_in
    use const

    implicit none

    contains

        subroutine read_geometry_in(geom_filename, charges, coords)
            ! Read in geometry data supplied in a file
            character(10), intent(in) :: geom_filename
            integer, allocatable, intent(out) :: charges(:)
            real(p), allocatable, intent(out) :: coords(:,:)

            integer :: i, iunit, ir, natoms

            ! stdout unit is 6
            iunit = 6

            ! newunit chooses an unused file handle, 'old' means existing file
            open(newunit=ir, file=geom_filename, status='old', form='formatted')
            read(ir, *) natoms
            
            ! Lazy allocation! 
            ! [TODO] make check_allocate subroutine
            allocate(charges(natoms))
            allocate(coords(3,natoms))
            
            do i = 1, natoms
                read(ir, *) charges(i), coords(:,i)
            end do
            
            close(ir)

        end subroutine read_geometry_in
            
end module geom_read_in
