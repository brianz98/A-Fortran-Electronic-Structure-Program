module system
   use const
   
   implicit none

   enum, bind(c)
      enumerator :: Hartree_Fock, MP_2, CC_SD, CC_SD_T
   end enum

   type system_t
      ! Basic system information
      integer :: natoms = 0
      integer :: nel = 0
      integer :: nbasis = 0
      integer :: nocc = 0
      integer :: nvirt = 0
      integer, allocatable :: charges(:)
      real(p), allocatable :: coords(:,:)

      ! Calculation info
      real(p) :: e_hf = 0.0_p
      real(p) :: e_mp2 = 0.0_p
      real(p) :: e_ccsd = 0.0_p
      real(p) :: e_ccsd_t = 0.0_p

      ! Coefficient matrix of the canonical HF orbitals in the original, nonorthogonal AO basis
      real(p), allocatable :: canon_coeff(:,:)
      real(p), allocatable :: canon_levels(:)
      real(p), allocatable :: canon_levels_spinorb(:)

      ! Tolerances
      real(p) :: scf_e_tol = 1e-6
      real(p) :: scf_d_tol = 1e-6
      integer :: scf_diis_n_errmat = 6
      real(p) :: ccsd_e_tol = 1e-6
      real(p) :: ccsd_t_tol = 1e-6
      integer :: ccsd_diis_n_errmat = 8
      integer :: scf_maxiter = 50
      integer :: ccsd_maxiter = 50

      ! Level of theory (up to)
      logical :: restricted
      integer :: calc_type

      ! There's currently support for 6 spatial CCSD(T) methods:
      ! CCSD[T]/(T), renormalised CCSD[T]/(T), completely renormalised CCSD[T]/(T)
      logical :: ccsd_t_paren = .false.
      logical :: ccsd_t_renorm = .false.
      logical :: ccsd_t_comp_renorm = .false.

      ! Do we additionally write a FCIDUMP file?
      logical :: write_fcidump = .false.

   end type system_t

   type state_t
      integer :: iter = 0
      real(p) :: energy = 0.0_p
      real(p) :: energy_old = 0.0_p
      real(p), allocatable :: density_old(:,:)
      real(p), allocatable :: t2_old(:,:,:,:)
   end type state_t

   contains
       subroutine read_system_in(sys)
           ! Reads in input file called `els.in` in the current working directory specifying the levels of theory and tolerance
           ! In/out:
           !    sys: system under study

           use error_handling, only: error
           type(system_t), intent(inout) :: sys
           real(p) :: scf_e_tol, scf_d_tol, ccsd_e_tol, ccsd_t_tol
           integer :: scf_diis_n_errmat, ccsd_diis_n_errmat, scf_maxiter, ccsd_maxiter
           logical :: write_fcidump
           character(40) :: calc_type

           integer :: ir, ierr

           ! Define namelist
           namelist /elsinput/ calc_type, scf_e_tol, scf_d_tol, scf_diis_n_errmat, ccsd_e_tol, ccsd_t_tol, ccsd_diis_n_errmat, &
                               scf_maxiter, ccsd_maxiter, write_fcidump

           ! Check if file exists
           inquire(file='els.in', iostat=ierr)
           if (ierr /= 0) call error('system::read_system_in', 'input file els.in does not exist')

           ! Read namelist
           open(newunit=ir, file='els.in', action='read')
           read(unit=ir, nml=elsinput, iostat=ierr)

           if (ierr /= 0) call error('system::read_system_in', 'invalid input file format!')

           close(ir)

           sys%scf_e_tol = scf_e_tol; sys%scf_d_tol = scf_d_tol; sys%scf_diis_n_errmat = scf_diis_n_errmat
           sys%ccsd_e_tol = ccsd_e_tol; sys%ccsd_t_tol = ccsd_t_tol; sys%ccsd_diis_n_errmat = ccsd_diis_n_errmat
           sys%scf_maxiter = scf_maxiter; sys%ccsd_maxiter = ccsd_maxiter; sys%write_fcidump = write_fcidump

           select case(trim(calc_type))
               case("RHF")
                  sys%calc_type = Hartree_Fock
                  sys%restricted = .true.
               case("UHF")
                  sys%calc_type = Hartree_Fock
                  sys%restricted = .false.
               case("MP2_spinorb")
                  sys%calc_type = MP_2
                  sys%restricted = .false.
               case("MP2_spatial")
                  sys%calc_type = MP_2
                  sys%restricted = .true.
               case("CCSD_spinorb")
                  sys%calc_type = CC_SD
                  sys%restricted = .false.
               case("CCSD_spatial")
                  sys%calc_type = CC_SD
                  sys%restricted = .true.
               case("CCSD(T)_spinorb")
                  sys%calc_type = CC_SD_T
                  sys%restricted = .false.
               case("CCSD(T)_spatial")
                  sys%calc_type = CC_SD_T
                  sys%restricted = .true.
                  sys%ccsd_t_paren = .true.
               case("CCSD[T]_spatial")
                  sys%calc_type = CC_SD_T
                  sys%restricted = .true.
               case("RCCSD(T)_spatial")
                  sys%calc_type = CC_SD_T
                  sys%restricted = .true.
                  sys%ccsd_t_paren = .true.
                  sys%ccsd_t_renorm = .true.
               case("RCCSD[T]_spatial")
                  sys%calc_type = CC_SD_T
                  sys%restricted = .true.
                  sys%ccsd_t_renorm = .true.
               case("CRCCSD(T)_spatial")
                  sys%calc_type = CC_SD_T
                  sys%restricted = .true.
                  sys%ccsd_t_paren = .true.
                  sys%ccsd_t_comp_renorm = .true.
               case("CRCCSD[T]_spatial")
                  sys%calc_type = CC_SD_T
                  sys%restricted = .true.
                  sys%ccsd_t_comp_renorm = .true.
               case default
                    call error('system::read_system_in', 'Unrecognised calculation type!')
           end select

        end subroutine

end module system
