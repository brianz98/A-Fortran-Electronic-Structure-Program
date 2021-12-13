program main

    use const
    use geom_read_in
    use geom_utils

    implicit none

    integer :: i, iunit
    character(20) :: geometry_filename
    integer, allocatable :: charges(:)
    real(p), allocatable :: coords(:,:), bondlengths(:,:)

    ! write to stdout
    iunit = 6

    write(iunit, *) 'Electronic Structure Programme'

    geometry_filename = 'geom.dat'

    call read_geometry_in(geometry_filename, charges, coords)
    
    write(iunit, *) 'The geometry is as follows'
    do i = 1, size(charges)
        write(iunit, *) charges(i), coords(:,i)
    end do

    call get_bond_lengths(size(charges), coords, bondlengths)

    write(iunit, *) 'The interatomic distances are as follows'
    do i = 1, size(charges)
        write(iunit, *) bondlengths(:,i)
    end do
    
end program main

