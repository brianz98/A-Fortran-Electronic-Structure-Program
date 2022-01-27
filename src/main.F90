program main

    use const
    use omp_lib
    use integrals, only: read_integrals_in, print_sys_info, int_store_t
    use system, only: system_t, read_system_in
    use system, only: RHF, UHF, MP2_spinorb, MP2_spatial, CCSD_spinorb, CCSD_spatial, CCSD_T_spinorb, CCSD_T_spatial
    use hf, only: do_rhf, do_uhf
    use geometry, only: read_geometry_in
    use mp2, only: do_mp2_spinorb, do_mp2_spatial
    use ccsd, only: do_ccsd_spinorb, do_ccsd_spatial, do_ccsd_t_spinorb, do_ccsd_t_spatial

    implicit none

    type(int_store_t) :: int_store
    type(system_t) :: sys
    integer :: iunit
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

    call read_system_in(sys)
    call read_integrals_in(sys, int_store)
    call read_geometry_in(sys, int_store)
    call print_sys_info(sys, int_store)

    call system_clock(t1)
    if (t1<t0) t1 = t1+count_max ! Moving time blocks, see system_clock documentation
    write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for system initialisation:', real(t1-t0)/count_rate, "s"
    t0=t1

    
    associate(ct=>sys%calc_type)
    if (any((/ct==MP2_spinorb, ct==CCSD_spinorb, ct==CCSD_T_spinorb/))) then
        ! Currently UHF is not implemented, so we convert spatial HF orbitals into spinorbitals
        call do_rhf(sys, int_store)
        call system_clock(t1)
        if (t1<t0) t1 = t1+count_max
        write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for restricted Hartree-Fock:', real(t1-t0)/count_rate, "s"
        t0=t1

        ! Spinorbital MP2 is not currently implemented
        call do_mp2_spatial(sys, int_store)
        call system_clock(t1)
        if (t1<t0) t1 = t1+count_max
        write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for UMP2:', real(t1-t0)/count_rate, "s"
        t0=t1
        
        if (any((/ct==CCSD_spinorb, ct==CCSD_T_spinorb/))) then
            call do_ccsd_spinorb(sys, int_store)
            call system_clock(t1)
            if (t1<t0) t1 = t1+count_max
            write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for UCCSD:', real(t1-t0)/count_rate, "s"
            t0=t1

            if (ct==CCSD_T_spinorb) then
                call do_ccsd_t_spinorb(sys, int_store)
                call system_clock(t1)
                if (t1<t0) t1 = t1+count_max
                write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for UCCSD(T):', real(t1-t0)/count_rate, "s"
                t0=t1
            end if
        end if
    end if

    if (any((/ct==MP2_spatial, ct==CCSD_spatial, ct==CCSD_T_spatial/))) then
        call do_rhf(sys, int_store)
        call system_clock(t1)
        if (t1<t0) t1 = t1+count_max
        write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for restricted Hartree-Fock:', real(t1-t0)/count_rate, "s"
        t0=t1

        call do_mp2_spatial(sys, int_store)
        call system_clock(t1)
        if (t1<t0) t1 = t1+count_max
        write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for RMP2:', real(t1-t0)/count_rate, "s"
        t0=t1
        
        if (any((/ct==CCSD_spatial, ct==CCSD_T_spatial/))) then
            call do_ccsd_spatial(sys, int_store)
            call system_clock(t1)
            if (t1<t0) t1 = t1+count_max
            write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for RCCSD:', real(t1-t0)/count_rate, "s"
            t0=t1

            if (ct==CCSD_T_spatial) then
                call do_ccsd_t_spatial(sys, int_store)
                call system_clock(t1)
                if (t1<t0) t1 = t1+count_max
                write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for RCCSD(T):', real(t1-t0)/count_rate, "s"
                t0=t1
            end if
        end if
    end if
    end associate

    write(iunit, '(1X, 64("="))')
    write(iunit, '(1X, A)') 'Final energy breakdown'
    write(iunit, '(1X, A, 1X, F12.7)') 'E_HF:                   ', sys%e_hf
    write(iunit, '(1X, A, 1X, F12.7)') 'E_MP2:                  ', sys%e_mp2
    write(iunit, '(1X, A, 1X, F12.7)') 'E_CCSD:                 ', sys%e_ccsd
    write(iunit, '(1X, A, 1X, F12.7)') 'E_CCSD(T):              ', sys%e_ccsd_t
    write(iunit, '(1X, 40("-"))')
    write(iunit, '(1X, A, 1X, F12.7)') 'Total electronic energy:', sys%e_hf+sys%e_ccsd+sys%e_ccsd_t
    write(iunit, '(1X, A, 1X, F12.7)') 'Total energy:           ', sys%e_hf+sys%e_ccsd+sys%e_ccsd_t + int_store%e_nuc

    call date_and_time(VALUES=date_values)
    write(iunit, '(1X, 64("="))')
    write(iunit, '(1X, A, 1X, I2.2, "/", I2.2, "/", I4.4, 1X, A, 1X, I2.2, 2(":", I2.2))') &
                "Finished running on", date_values(3:1:-1), "at", date_values(5:7)
   
end program main
