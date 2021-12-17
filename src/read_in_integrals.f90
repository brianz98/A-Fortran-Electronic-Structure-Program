module read_in_integrals
    use const

    implicit none

    type int_store_t
        real(p) :: e_nuc = 0.0_p
        real(p), allocatable :: ovlp(:,:)
        real(p), allocatable :: ke(:,:)
        real(p), allocatable :: ele_nuc(:,:)
        real(p), allocatable :: core_hamil(:,:)
        ! Too big, we compress with symmetry
        real(p), allocatable :: eri(:)
    end type int_store_t
    
    contains
        subroutine read_integrals_in(sys, int_store)
            ! Read in integral data supplied in separate files
            use, intrinsic :: iso_fortran_env, only: iostat_end
            use system

            type(int_store_t), intent(inout) :: int_store
            type(system_t), intent(inout) :: sys

            character(20) :: enuc_f, ovlp_f, ke_f, ele_nuc_f, eri_f
            integer :: i, j, iunit, ir, natoms, ios, nbasis, ibasis, jbasis, abasis, bbasis
            real(p) :: intgrl

            ! stdout unit is 6
            iunit = 6

            ! I/O status
            ios = 0

            nbasis = 0

            ! Hard-coded filenames
            enuc_f = 'dat/enuc.dat'
            ovlp_f = 'dat/s.dat'
            ke_f = 'dat/t.dat'
            ele_nuc_f = 'dat/v.dat'
            eri_f = 'dat/eri.dat'

            write(iunit, '(1X, 16("-"))')
            write(iunit, '(1X, A)') 'Integral read-in'
            write(iunit, '(1X, 16("-"))')


            !########### Get basis info ###############
            write(iunit, *) 'Getting number of basis functions...'
            open(newunit=ir, file=ovlp_f, status='old', form='formatted')
            ! Pass through the file once to find nbasis first
            do
                read(ir, *, iostat=ios) ibasis, jbasis, intgrl
                ! If EOF reached
                if (ios==iostat_end) exit
                if (ibasis>nbasis .or. jbasis>nbasis) then
                    nbasis = max(ibasis, jbasis)
                end if
            end do

            sys%nbasis = nbasis

            write(iunit, *) 'Allocating integral store...'
            call init_int_store(sys, int_store)
            
            close(ir)

            !############# Overlap matrix ##################
            write(iunit, *) 'Reading overlap matrix...'
            open(newunit=ir, file=ovlp_f, status='old', form='formatted')

            do
                read(ir, *, iostat=ios) ibasis, jbasis, intgrl
                ! If EOF reached
                if (ios==iostat_end) exit
                ! Simple permutational symmetry
                int_store%ovlp(ibasis,jbasis) = intgrl
                int_store%ovlp(jbasis,ibasis) = intgrl
            end do
            close(ir)

            !############ Kinetic energy ############## 
            write(iunit, *) 'Reading kinetic integrals...'
            open(newunit=ir, file=ke_f, status='old', form='formatted')

            do
                read(ir, *, iostat=ios) ibasis, jbasis, intgrl
                ! If EOF reached
                if (ios==iostat_end) exit
                ! Simple permutational symmetry
                int_store%ke(ibasis,jbasis) = intgrl
                int_store%ke(jbasis,ibasis) = intgrl
            end do
            close(ir)

            !########### Nuclear-electron integral ######
            write(iunit, *) 'Reading nuclear-electron integrals...'
            open(newunit=ir, file=ele_nuc_f, status='old', form='formatted')

            do
                read(ir, *, iostat=ios) ibasis, jbasis, intgrl
                ! If EOF reached
                if (ios==iostat_end) exit
                ! Simple permutational symmetry
                int_store%ele_nuc(ibasis,jbasis) = intgrl
                int_store%ele_nuc(jbasis,ibasis) = intgrl
            end do
            close(ir)

            !######### Core Hamiltonian ##############
            write(iunit, *) 'Constructing core Hamiltonian...'
            int_store%core_hamil = int_store%ke + int_store%ele_nuc

            !######### ERI / Two-body integrals ########
            write(iunit, *) 'Reading two-body integrals...'
            open(newunit=ir, file=eri_f, status='old', form='formatted')

            do
                read(ir, *, iostat=ios) ibasis, jbasis, abasis, bbasis, intgrl
                ! If EOF reached
                if (ios==iostat_end) exit
                ! 8-fold permutational symmetry, we use compression to save space
                int_store%eri(eri_ind(ibasis, jbasis, abasis, bbasis)) = intgrl
            end do
            close(ir)

            write(iunit, *) 'Done reading integrals!'

        end subroutine read_integrals_in

        subroutine init_int_store(sys, int_store)
            use system

            type(system_t), intent(in) :: sys
            type(int_store_t), intent(inout) :: int_store

            integer :: neri

            neri = sys%nbasis*(sys%nbasis-1)/2 + sys%nbasis
            neri = neri*(neri-1)/2 + neri

            associate(nb=>sys%nbasis, is=>int_store)
                allocate(is%ovlp(nb,nb), source=0.0_p)
                allocate(is%ke(nb,nb), source=0.0_p)
                allocate(is%ele_nuc(nb,nb), source=0.0_p)
                allocate(is%core_hamil(nb,nb), source=0.0_p)
                allocate(is%eri(neri), source=0.0_p)
            end associate
        end subroutine init_int_store

        elemental function eri_ind(i, j, a, b) result(ind)
            integer, intent(in) :: i, j, a, b
            integer :: ind

            integer :: ij, ab

            if (i >= j) then
                ij = i*(i-1)/2 + j
            else
                ij = j*(j-1)/2 + i
            end if

            if (a >= b) then
                ab = a*(a-1)/2 + b
            else
                ab = b*(b-1)/2 + a
            end if

            if (ij >= ab) then
                ind = ij*(ij-1)/2 + ab
            else
                ind = ab*(ab-1)/2 + ij
            end if
        end function eri_ind

        subroutine print_sys_info(sys, int_store)
            use system

            type(system_t), intent(in) :: sys
            type(int_store_t), intent(in) :: int_store

            integer :: iunit = 6
            integer :: i

            write(iunit, '(1X, 20("-"))')
            write(iunit, '(1X, A)') 'System information'
            write(iunit, '(1X, 20("-"))')
            write(iunit, '(1X, A, 1X, I0)') 'Number of electrons:', sys%nel
            write(iunit, '(1X, A, 1X, I0)') 'Number of basis functions:', sys%nbasis
            write(iunit, '(1X, A, 1X, ES15.8)') 'E_nuc:', int_store%e_nuc

        end subroutine print_sys_info

end module read_in_integrals
