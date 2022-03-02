program main

    use const
    use omp_lib
    use, intrinsic :: iso_fortran_env, only : iunit=>output_unit
    use integrals, only: read_integrals_in, print_sys_info, int_store_t, int_store_cc_t
    use system, only: system_t, read_system_in
    use system, only: Hartree_Fock, MP_2, CC_SD, CC_SD_T
    use hf, only: do_rhf, do_uhf
    use geometry, only: read_geometry_in
    use mp2, only: do_mp2_spinorb, do_mp2_spatial
    use ccsd, only: do_ccsd_spinorb, do_ccsd_spatial, do_ccsd_t_spinorb, do_ccsd_t_spatial

    implicit none

    type(int_store_t) :: int_store
    type(int_store_cc_t) :: int_store_cc
    type(system_t) :: sys
    integer(kind=8) :: t0, t1, count_rate, count_max
    integer :: date_values(8)
    character(80) :: calcname
    real(p) :: highest_theory_energy

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
    ! The spinorbital formulations, as UHF is not yet available, HF and MP2 are actually spatial, and then converted to
    ! spin-orbital basis for CCSD and CCSD(T)
    if (.not. sys%restricted) then
        ! Currently UHF is not implemented, so we convert spatial HF orbitals into spinorbitals
        call do_rhf(sys, int_store)
        call system_clock(t1)
        if (t1<t0) t1 = t1+count_max
        write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for restricted Hartree-Fock:', real(t1-t0)/count_rate, "s"
        t0=t1

        if (any(ct == [MP_2, CC_SD, CC_SD_T])) then
            ! Spinorbital MP2 is not currently implemented
            call do_mp2_spatial(sys, int_store)
            call system_clock(t1)
            if (t1<t0) t1 = t1+count_max
            write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for restricted MP2:', real(t1-t0)/count_rate, "s"
            t0=t1
            
            if (any(ct == [CC_SD, CC_SD_T])) then
                call do_ccsd_spinorb(sys, int_store, int_store_cc)
                call system_clock(t1)
                if (t1<t0) t1 = t1+count_max
                write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for unrestricted CCSD:', real(t1-t0)/count_rate, "s"
                t0=t1

                if (ct == CC_SD_T) then
                    call do_ccsd_t_spinorb(sys, int_store, int_store_cc)
                    call system_clock(t1)
                    if (t1<t0) t1 = t1+count_max
                    write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for unrestricted CCSD(T):', real(t1-t0)/count_rate, "s"
                    t0=t1
                end if
            end if
        end if
    end if

    if (sys%restricted) then
        call do_rhf(sys, int_store)
        call system_clock(t1)
        if (t1<t0) t1 = t1+count_max
        write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for restricted Hartree-Fock:', real(t1-t0)/count_rate, "s"
        t0=t1

        if (any(ct == [MP_2, CC_SD, CC_SD_T])) then 
            call do_mp2_spatial(sys, int_store)
            call system_clock(t1)
            if (t1<t0) t1 = t1+count_max
            write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for restricted MP2:', real(t1-t0)/count_rate, "s"
            t0=t1
            
            if (any(ct == [CC_SD, CC_SD_T])) then
                call do_ccsd_spatial(sys, int_store, int_store_cc)
                call system_clock(t1)
                if (t1<t0) t1 = t1+count_max
                write(iunit, '(1X, A, 1X, F7.4, A)') 'Time taken for restricted CCSD:', real(t1-t0)/count_rate, "s"
                t0=t1

                if (ct == CC_SD_T) then
                    call do_ccsd_t_spatial(sys, int_store, calcname, int_store_cc)
                    call system_clock(t1)
                    if (t1<t0) t1 = t1+count_max
                    write(iunit,'(1X, A, 1X, F7.4, A)')'Time taken for restricted '//trim(calcname)//':',real(t1-t0)/count_rate,"s"
                    t0=t1
                end if
            end if
        end if
    end if
    

    write(iunit, '(1X, 64("="))')
    write(iunit, '(1X, A)') 'Final energy breakdown'
    write(iunit, '(1X, A, 1X, F15.10)') 'RHF energy:                    ', sys%e_hf + int_store%e_nuc

    if (any(ct == [MP_2, CC_SD, CC_SD_T])) then
    write(iunit, '(1X, A, 1X, F15.10)') 'MP2 correlation energy:        ', sys%e_mp2
    write(iunit, '(1X, A, 1X, F15.10)') 'MP2 energy:                    ', sys%e_mp2 + sys%e_hf + int_store%e_nuc
    if (any(ct == [CC_SD, CC_SD_T])) then
    write(iunit, '(1X, A, 1X, F15.10)') 'CCSD correlation energy:       ', sys%e_ccsd
    write(iunit, '(1X, A, 1X, F15.10)') 'CCSD energy:                   ', sys%e_ccsd + sys%e_hf + int_store%e_nuc
    if (ct == CC_SD_T) then
    write(iunit, '(1X, A, 1X, F15.10)') 'CCSD[T] correlation energy:    ', sys%e_ccsd_t
    write(iunit, '(1X, A, 1X, F15.10)') 'CCSD[T] energy:                ', sys%e_ccsd_t + sys%e_hf + int_store%e_nuc
    if (sys%ccsd_t_paren) then
    write(iunit, '(1X, A, 1X, F15.10)') 'CCSD(T) correlation energy:    ', sys%e_ccsd_tt
    write(iunit, '(1X, A, 1X, F15.10)') 'CCSD(T) energy:                ', sys%e_ccsd_tt + sys%e_hf + int_store%e_nuc
    end if
    if (sys%ccsd_t_renorm .or. sys%ccsd_t_comp_renorm) then
    write(iunit, '(1X, A, 1X, F15.10)') 'R-CCSD[T] correlation energy:  ', sys%e_rccsd_t
    write(iunit, '(1X, A, 1X, F15.10)') 'R-CCSD[T] energy:              ', sys%e_rccsd_t + sys%e_hf + int_store%e_nuc
    if (sys%ccsd_t_paren) then
    write(iunit, '(1X, A, 1X, F15.10)') 'R-CCSD(T) correlation energy:  ', sys%e_rccsd_tt
    write(iunit, '(1X, A, 1X, F15.10)') 'R-CCSD(T) energy:              ', sys%e_rccsd_tt + sys%e_hf + int_store%e_nuc
    end if
    if (sys%ccsd_t_comp_renorm) then
    write(iunit, '(1X, A, 1X, F15.10)') 'CR-CCSD[T] correlation energy: ', sys%e_crccsd_t
    write(iunit, '(1X, A, 1X, F15.10)') 'CR-CCSD[T] energy:             ', sys%e_crccsd_t + sys%e_hf + int_store%e_nuc
    if (sys%ccsd_t_paren) then
    write(iunit, '(1X, A, 1X, F15.10)') 'CR-CCSD(T) correlation energy: ', sys%e_crccsd_tt
    write(iunit, '(1X, A, 1X, F15.10)') 'CR-CCSD(T) energy:             ', sys%e_crccsd_tt + sys%e_hf + int_store%e_nuc
    end if
    end if
    end if
    end if
    end if
    end if
    write(iunit, '(1X, 47("-"))')
    if (any(ct == [CC_SD, CC_SD_T]) .and. sys%restricted) then
    write(iunit, '(1X, A, 1X, F15.10)') 'T1 diagnostic:                 ', sys%t1_diagnostic
    end if
    if (sys%ccsd_t_renorm .or. sys%ccsd_t_comp_renorm) then
    write(iunit, '(1X, A, 1X, F15.10)') 'D[T]:                          ', sys%D_T
    if (sys%ccsd_t_paren) then
    write(iunit, '(1X, A, 1X, F15.10)') 'D(T):                          ', sys%D_TT
    end if
    end if
    write(iunit, '(1X, 47("-"))')
    write(iunit, '(1X, A, 1X, F15.10)') 'Total electronic energy:       ', sys%e_hf+sys%e_highest
    write(iunit, '(1X, A, 1X, F15.10)') 'Nuclear repulsion:             ', int_store%e_nuc
    write(iunit, '(1X, A, 1X, F15.10)') 'Total energy:                  ', sys%e_hf+ sys%e_highest + int_store%e_nuc
    
    end associate

    call date_and_time(VALUES=date_values)
    write(iunit, '(1X, 64("="))')
    write(iunit, '(1X, A, 1X, I2.2, "/", I2.2, "/", I4.4, 1X, A, 1X, I2.2, 2(":", I2.2))') &
                "Finished running on", date_values(3:1:-1), "at", date_values(5:7)
   
end program main
