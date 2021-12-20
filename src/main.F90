program main

    use const
    use integrals, only: read_integrals_in, print_sys_info, int_store_t
    use system, only: system_t
    use hf, only: do_hartree_fock
    use geometry, only: read_geometry_in
    use mp2, only: do_mp2
    use ccsd, only: do_ccsd

    implicit none

    type(int_store_t) :: int_store
    type(system_t) :: sys
    integer :: i, iunit
    integer(kind=8) :: t0, t1, count_rate, count_max
    integer :: date_values(8)

    ! write to stdout
    iunit = 6

    write(iunit, '(1X, 64("="))')
    write(iunit, '(1X, A)') 'A Fortran Electronic Structure Programme (AFESP)'
    write(iunit, '(1X, 64("="))')

    call date_and_time(VALUES=date_values)
    write(iunit, '(1X, A, 1X, I2.2, "/", I2.2, "/", I4.4, 1X, A, 1X, I2.2, 2(":", I2.2))') &
                "Started running on", date_values(3:1:-1), "at", date_values(5:7)

    call system_clock(t0, count_rate, count_max)

    call read_integrals_in(sys, int_store)
    call read_geometry_in(sys, int_store)
    call print_sys_info(sys, int_store)
    call system_clock(t1)
    if (t1<t0) t1 = t1+count_max ! Moving time blocks, see system_clock documentation
    write(iunit, '(1X, A, 1X, F7.5, A)') 'Time taken for system initialisation:', real(t1-t0)/count_rate, "s"
    t0=t1

    call do_hartree_fock(sys, int_store)
    call system_clock(t1)
    if (t1<t0) t1 = t1+count_max
    write(iunit, '(1X, A, 1X, F7.5, A)') 'Time taken for Hartree-Fock:', real(t1-t0)/count_rate, "s"
    t0=t1

    call do_mp2(sys, int_store)
    call system_clock(t1)
    if (t1<t0) t1 = t1+count_max
    write(iunit, '(1X, A, 1X, F7.5, A)') 'Time taken for MP2:', real(t1-t0)/count_rate, "s"
    t0=t1

    call do_ccsd(sys, int_store)
    call system_clock(t1)
    if (t1<t0) t1 = t1+count_max
    write(iunit, '(1X, A, 1X, F7.5, A)') 'Time taken for CCSD:', real(t1-t0)/count_rate, "s"
    t0=t1

    write(iunit, '(1X, 64("="))')
    write(iunit, '(1X, A)') 'Final energy breakdown'
    write(iunit, '(1X, A, 1X, F12.7)') 'E_HF:                   ', sys%e_hf
    write(iunit, '(1X, A, 1X, F12.7)') 'E_MP2:                  ', sys%e_mp2
    write(iunit, '(1X, A, 1X, F12.7)') 'E_CCSD:                 ', sys%e_ccsd
    write(iunit, '(1X, A, 1X, F12.7)') 'E_CCSD(T):              ', sys%e_ccsd_t
    write(iunit, '(1X, 40("-"))')
    write(iunit, '(1X, A, 1X, F12.7)') 'Total electronic energy:', sys%e_hf+sys%e_mp2+sys%e_ccsd+sys%e_ccsd_t
    write(iunit, '(1X, A, 1X, F12.7)') 'Total energy:           ', sys%e_hf+sys%e_mp2+sys%e_ccsd+sys%e_ccsd_t + int_store%e_nuc

    call date_and_time(VALUES=date_values)
    write(iunit, '(1X, 64("="))')
    write(iunit, '(1X, A, 1X, I2.2, "/", I2.2, "/", I4.4, 1X, A, 1X, I2.2, 2(":", I2.2))') &
                "Finished running on", date_values(3:1:-1), "at", date_values(5:7)

    
end program main