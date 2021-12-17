program main

    use const
    use read_in_integrals, only: read_integrals_in, print_sys_info, int_store_t
    use system, only: system_t
    use hf, only: do_hartree_fock
    use geometry, only: read_geometry_in

    implicit none

    type(int_store_t) :: int_store
    type(system_t) :: sys
    integer :: i, iunit

    ! write to stdout
    iunit = 6

    write(iunit, '(1X, 64("="))')
    write(iunit, '(1X, A)') 'A Fortran Electronic Structure Programme (AFESP)'
    write(iunit, '(1X, 64("="))')


    call read_integrals_in(sys, int_store)
    call read_geometry_in(sys, int_store)

    call print_sys_info(sys, int_store)

    call do_hartree_fock(sys, int_store)
    
end program main

