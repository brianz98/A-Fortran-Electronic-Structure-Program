module system
   use const
   
   implicit none

   type system_t
      ! Basic system information
      integer :: nel = 0
      integer :: nbasis = 0
      integer :: natoms = 0
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

      ! Converged CCSD amplitudes
      real(p), allocatable :: t1(:,:)
      real(p), allocatable :: t2(:,:,:,:)

      ! Tolerances
      real(p) :: scf_e_tol = 1e-6
      real(p) :: scf_d_tol = 1e-6
      real(p) :: ccsd_e_tol = 1e-6
      real(p) :: ccsd_t_tol = 1e-6

   end type system_t

   type state_t
      integer :: iter = 0
      real(p) :: energy = 0.0_p
      real(p) :: energy_old = 0.0_p
      real(p), allocatable :: density_old(:,:)
      real(p), allocatable :: t2_old(:,:,:,:)
   end type state_t

end module system