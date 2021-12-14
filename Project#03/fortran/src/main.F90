program main

    use const
    use read_in_integrals
    use system
    use hf

    implicit none

    type(int_store_t) :: int_store
    type(system_t) :: sys
    integer :: i, iunit

    ! write to stdout
    iunit = 6

    write(iunit, *) 'Electronic Structure Programme'

    call read_integrals_in(sys, int_store)
    call print_int_info(sys, int_store)

    call do_hartree_fock(sys, int_store)
    
end program main

