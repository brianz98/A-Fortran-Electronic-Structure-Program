module integrals
    use const

    implicit none

    public
    private :: istart, bignum

    ! These define a module-internal lookup array of triangular numbers, used in eri_ind below
    ! This is to avoid the repeated calculation of i(i+1)/2.
    ! istart stores the flattened index of the first elements of each row of a lower triangular matrix minus one:
    ! 1
    ! 2 3               ==> istart = [0, 1, 3, 6, ...]
    ! 4 5 6
    ! 7 8 9 10
    ! So when we want to look up the array element (4,2) of the triangular array, stored in a flattened array, 
    ! the index of element at (4,2) in the flattened array is just istart(4) + 2 = 8

    ! bignum should be bigger than n(n-1)/2 where n is the number of basis functions, so currently we can accommodate ~100 basis functions
    integer, parameter :: bignum = 5000
    integer :: istart(bignum)

    type int_store_t
        real(p) :: e_nuc = 0.0_p
        real(p), allocatable :: ovlp(:,:)
        real(p), allocatable :: ke(:,:)
        real(p), allocatable :: ele_nuc(:,:)
        real(p), allocatable :: core_hamil(:,:)
        ! Too big, we compress with symmetry
        real(p), allocatable :: eri(:)
        real(p), allocatable :: eri_mo(:)
        real(p), allocatable :: asym_spinorb(:,:,:,:)
        ! needed for CCSD -> CCSD(T), spinorbital formulation
        real(p), dimension(:,:,:,:), allocatable :: vvoo, vovv, ovoo
    end type int_store_t
    
    contains
        subroutine read_integrals_in(sys, int_store)
            ! Read in integral data supplied in separate files
            use, intrinsic :: iso_fortran_env, only: iostat_end
            use system

            type(int_store_t), intent(inout) :: int_store
            type(system_t), intent(inout) :: sys

            character(20) :: enuc_f, ovlp_f, ke_f, ele_nuc_f, eri_f
            integer :: iunit, ir, ios, nbasis
            integer :: ibasis, jbasis, abasis, bbasis, ij_ind, ab_ind
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

            call init_istart()

            do
                read(ir, *, iostat=ios) ibasis, jbasis, abasis, bbasis, intgrl
                ! If EOF reached
                if (ios==iostat_end) exit
                ! 8-fold permutational symmetry, we use compression to save space
                ij_ind = eri_ind(ibasis, jbasis)
                ab_ind = eri_ind(abasis, bbasis)
                int_store%eri(eri_ind(ij_ind, ab_ind)) = intgrl
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

        subroutine init_istart()
            integer :: i

            istart(1) = 0
                do i = 2, bignum
                    istart(i) = istart(i-1) + i-1
                end do
        end subroutine init_istart

        elemental function eri_ind(i, j) result(ind)
            ! This function can be a general version for i,j,a,b, but 
            ! we can also use storage to cut down on unnecessary instructrions:
            ! i.e. by storing intermediate results ij and ab

            integer, intent(in) :: i, j
            integer :: ind

            if (i >= j) then
                ind = istart(i) + j
            else
                ind = istart(j) + i
            end if

        end function eri_ind

        subroutine print_sys_info(sys, int_store)
            use system

            type(system_t), intent(in) :: sys
            type(int_store_t), intent(in) :: int_store

            integer :: iunit = 6

            write(iunit, '(1X, 20("-"))')
            write(iunit, '(1X, A)') 'System information'
            write(iunit, '(1X, 20("-"))')
            write(iunit, '(1X, A, 1X, I0)') 'Number of electrons:', sys%nel
            write(iunit, '(1X, A, 1X, I0)') 'Number of basis functions:', sys%nbasis
            write(iunit, '(1X, A, 1X, ES15.8)') 'E_nuc:', int_store%e_nuc
            write(iunit, '(1X, A, 1X, ES8.2)') 'scf_e_tol:', sys%scf_e_tol
            write(iunit, '(1X, A, 1X, ES8.2)') 'scf_d_tol:', sys%scf_d_tol
            write(iunit, '(1X, A, 1X, ES8.2)') 'ccsd_e_tol:', sys%ccsd_e_tol
            write(iunit, '(1X, A, 1X, ES8.2)') 'ccsd_t_tol:', sys%ccsd_t_tol

        end subroutine print_sys_info

end module integrals
