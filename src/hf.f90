module hf
   use const

   implicit none

   type diis_t
      ! Contains matrices for use in the DIIS procedure
      integer :: n_errmat = 6
      integer :: n_active = 0
      integer :: iter = 0
      real(p), allocatable :: e(:,:,:)
      real(p), allocatable :: F(:,:,:)
      real(p), allocatable :: B(:,:)
      real(p), allocatable :: c(:)
      real(p), allocatable :: rhs(:)
   end type diis_t

   contains
      subroutine do_hartree_fock(sys, int_store)
         use error_handling
         use system
         use integrals
         use linalg

         type(int_store_t), intent(in) :: int_store
         type(system_t), intent(inout) :: sys

         type(mat_t) :: ovlp, fockmat
         real(p) :: ortmat(sys%nbasis,sys%nbasis)
         real(p) :: fock(sys%nbasis,sys%nbasis)
         real(p) :: coeff(sys%nbasis,sys%nbasis)
         real(p) :: density(sys%nbasis,sys%nbasis)
         real(p) :: tmpmat(sys%nbasis,sys%nbasis)

         type(state_t) :: st
         integer :: iter, maxiter, iunit, i, j, ierr
         logical :: conv = .false.
         type(diis_t) :: diis

         iunit = 6
         maxiter = 100
         write(iunit, '(1X, 12("-"))')
         write(iunit, '(1X, A)') 'Hartree-Fock'
         write(iunit, '(1X, 12("-"))')

         allocate(st%density_old(sys%nbasis,sys%nbasis), source=0.0_p)

         ! Build the orthogonalisation matrix S^-1/2
         ! Initialise the overlap matrix stored in int_store_t into a mat_t object
         call init_mat_t(int_store%ovlp, ovlp)
         ovlp%ao_ort = ovlp%ao

         ! S^(-1/2) = C * L^(-1/2) * C^(T)
         call diagonalise(ovlp)
         call build_ortmat(ovlp, ortmat)

         call deallocate_mat_t(ovlp)

         ! Initial guess Fock matrix is H_core
         ! F'_0 = S^T(-1/2) * F * S^(-1/2)
         call init_mat_t(int_store%core_hamil, fockmat)

         call init_diis_t(sys, diis)

         do iter = 1, maxiter
            fockmat%ao_ort = matmul(transpose(ortmat), matmul(fockmat%ao, ortmat))

            call diagonalise(fockmat)

            ! Transform coefficient into nonorthogonal AO basis
            fockmat%A = transpose(matmul(ortmat, fockmat%A))

            call build_density(fockmat%A, density, sys%nel/2)

            call update_scf_energy(density, fockmat%ao, st, int_store, conv)
            if (conv) then
               write(iunit, '(1X, A)') 'Convergence reached within tolerance.'
               write(iunit, '(1X, A, 1X, F15.8)') 'Final SCF Energy (Hartree):', st%energy
               write(iunit, '(1X, A)') 'Orbital energies (Hartree):'
               do i = size(fockmat%W), 1, -1
                  write(iunit, '(1X, I3, 1X, F15.8)') i, fockmat%W(i)
               end do

               ! Copied for return / use in MP2
               sys%e_hf = st%energy
               allocate(sys%canon_coeff, source=fockmat%A)
               allocate(sys%canon_levels, source=fockmat%W)
               exit
            end if
            write(iunit, '(1X, A, 1X, I3, F15.8)') 'Iteration', iter, st%energy

            call build_fock(sys, density, fockmat%ao, int_store)

            ! DIIS update
            call update_diis(diis, fockmat%ao, density, int_store%ovlp)
         end do

         if (.not. conv) then
            write(iunit, '(1X, A)') 'Convergence not reached, please increase maxiter.'
         end if

         call deallocate_diis_t(diis)
         deallocate(st%density_old)

      end subroutine do_hartree_fock

      subroutine update_diis(diis, fock, density, ovlp)
         use linalg, only: linsolve
         use error_handling, only: error

         type(diis_t), intent(inout) :: diis
         real(p), intent(inout) :: fock(:,:)
         real(p), intent(in) :: ovlp(:,:), density(:,:)

         integer :: i, j, ierr

         diis%iter = diis%iter+1
         if (diis%iter > 6) diis%iter = diis%iter - 6
         if (diis%n_active < diis%n_errmat) diis%n_active = diis%n_active+1
         diis%F(diis%iter,:,:) = 0.0_p
         diis%e(diis%iter,:,:) = 0.0_p
         diis%F(diis%iter,:,:) = fock(:,:)
         diis%e(diis%iter,:,:) = matmul(fock, matmul(density,ovlp)) &
                               - matmul(ovlp, matmul(density, fock))

         associate(n=>diis%n_active, nerr=>diis%n_errmat)
         if (n > 1) then
            ! Construct the B matrix
            if (n <= nerr) then
               if (allocated(diis%B)) deallocate(diis%B, diis%c, diis%rhs)
               allocate(diis%B(n+1,n+1), diis%c(n+1), diis%rhs(n+1), source=0.0_p)
            end if
            diis%B(n+1,:) = -1.0_p
            diis%B(n+1,n+1) = 0.0_p
            diis%rhs(n+1) = -1.0_p
            diis%c = diis%rhs
            do i = 1, n
               do j = 1, i
                  diis%B(i,j) = sum(diis%e(i,:,:)*diis%e(j,:,:))
               end do
            end do
            call linsolve(diis%B, diis%c, ierr)
            if (ierr /= 0) call error('hf::update_diis', 'Linear solve failed!')
            fock = 0.0_p
            do i = 1, n
               fock = fock + diis%c(i) * diis%F(i,:,:)
            end do
         end if
         end associate
      end subroutine update_diis

      subroutine init_diis_t(sys, diis)
         use system, only: system_t

         type(system_t), intent(in) :: sys
         type(diis_t), intent(out) :: diis

         associate(n=>sys%nbasis)
            allocate(diis%e(diis%n_errmat,n,n), source=0.0_p)
            allocate(diis%F(diis%n_errmat,n,n), source=0.0_p)
         end associate
      end subroutine init_diis_t

      subroutine deallocate_diis_t(diis)
         type(diis_t), intent(inout) :: diis

         deallocate(diis%e, diis%F, diis%B, diis%c, diis%rhs)
      end subroutine deallocate_diis_t

      subroutine diagonalise(mat)
         ! For a symmetric matrix A we produce
         ! A = C * L * C^T
         ! With L and C stored in a mat_t data structure

         use linalg
         use error_handling

         type(mat_t), intent(inout) :: mat

         integer :: ierr

         ! Overwrite the previous coeff matrix
         mat%A = mat%ao_ort
         call eigs(mat, ierr)
         if (ierr /= 0) call error('hf::diagonalise', 'Diagonalisation failed!')
         !call zero_mat(mat%A)

      end subroutine diagonalise

      subroutine init_mat_t(matrix, mat)
         use integrals
         use linalg, only: mat_t
         use system
        
         real(p), intent(in) :: matrix(:,:)
         type(mat_t), intent(out) :: mat

         associate(N=>size(matrix, dim=1))
            mat%N = N
            allocate(mat%ao(N,N))
            allocate(mat%A(N,N))
            allocate(mat%ao_ort(N,N))
            mat%ao(:,:) = matrix(:,:)
            allocate(mat%W(N))
         end associate

      end subroutine init_mat_t

      subroutine deallocate_mat_t(mat)
         use linalg, only: mat_t

         type(mat_t), intent(inout) :: mat

         mat%N = 0
         deallocate(mat%ao)
         deallocate(mat%ao_ort)
         deallocate(mat%A)
         deallocate(mat%W)
      end subroutine deallocate_mat_t

      subroutine build_ortmat(ovlp, ortmat)

         use linalg, only: mat_t

         type(mat_t), intent(in) :: ovlp
         real(p), intent(out) :: ortmat(:,:)

         real(p) :: tmpmat(size(ortmat, dim=1), size(ortmat, dim=1))

         integer :: i
         
         tmpmat(:,:) = 0.0_p
         ortmat(:,:) = 0.0_p

         associate(diag=>ovlp%W)
            do i = 1, size(ortmat, dim=1)
               ortmat(i,i) = 1/sqrt(diag(i))
            end do
         end associate
         
         tmpmat = matmul(ortmat, transpose(ovlp%A))
         ortmat = matmul(ovlp%A, tmpmat)

      end subroutine build_ortmat

      subroutine build_density(coeff, density, nocc)
         real(p), intent(in) :: coeff(:,:)
         integer, intent(in) :: nocc
         real(p), intent(out) :: density(:,:)

         density = matmul(transpose(coeff(1:nocc, :)), coeff(1:nocc, :))

      end subroutine build_density

      subroutine update_scf_energy(density, fock, st, int_store, conv)
         
         use system, only: state_t
         use integrals, only: int_store_t

         real(p), intent(in) :: density(:,:)
         real(p), intent(in) :: fock(:,:)
         type(int_store_t), intent(in) :: int_store
         type(state_t), intent(inout) :: st
         logical, intent(inout) :: conv

         st%energy_old = st%energy
         st%energy = sum(density * (int_store%core_hamil + fock))

         if (sqrt(sum((density-st%density_old)**2)) < 1e-7 .and. abs(st%energy-st%energy_old) < 1e-7) conv = .true.
         st%density_old = density

      end subroutine update_scf_energy
         
      subroutine build_fock(sys, density, fock, int_store)

         use integrals, only: int_store_t, eri_ind
         use system, only: system_t

         real(p), intent(in) :: density(:,:)
         type(int_store_t), intent(in) :: int_store
         type(system_t), intent(in) :: sys
         real(p), intent(out) :: fock(:,:)

         integer :: i, j, k, l
         integer :: ij, kl, ik, jl

         associate(eri=>int_store%eri)
            !$omp parallel do default(none) &
            !$omp schedule(dynamic, 2) collapse(2) &
            !$omp private(ij, jl, kl, ik) &
            !$omp shared(sys, density, fock, int_store)
            do j = 1, sys%nbasis
               do i = 1, sys%nbasis
                  fock(i,j) = int_store%core_hamil(i,j)
                  ij = eri_ind(i,j)
                  do l = 1, sys%nbasis
                     jl = eri_ind(j,l)
                     do k = 1, sys%nbasis
                        ! [TODO]: This step can be sped up with pre-computed lookup arrays
                        kl = eri_ind(k,l)
                        ik = eri_ind(i,k)
                        fock(i,j) = fock(i,j) + density(k,l) * (2*eri(eri_ind(ij,kl)) - eri(eri_ind(ik,jl)))
                     end do
                  end do
               end do
            end do
            !$omp end parallel do
         end associate

      end subroutine build_fock

end module hf
