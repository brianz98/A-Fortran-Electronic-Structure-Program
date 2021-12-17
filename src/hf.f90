module hf
   use const

   implicit none

   contains
      subroutine do_hartree_fock(sys, int_store)
         use error_handling
         use system
         use read_in_integrals
         use linalg

         type(int_store_t), intent(in) :: int_store
         type(system_t), intent(in) :: sys

         type(mat_t) :: ovlp, fockmat
         real(p) :: ortmat(sys%nbasis,sys%nbasis)
         real(p) :: fock(sys%nbasis,sys%nbasis)
         real(p) :: coeff(sys%nbasis,sys%nbasis)
         real(p) :: density(sys%nbasis,sys%nbasis)
         real(p) :: tmpmat(sys%nbasis,sys%nbasis)

         type(state_t) :: st
         integer :: iter, maxiter, iunit, i
         logical :: conv = .false.

         iunit = 6
         maxiter = 50
         write(iunit, '(1X, 12("-"))')
         write(iunit, '(1X, A)') 'Hartree-Fock'
         write(iunit, '(1X, 12("-"))')

         allocate(st%density(sys%nbasis,sys%nbasis), source=0.0_p)
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

         do iter = 1, maxiter
            fockmat%ao_ort = matmul(transpose(ortmat), matmul(fockmat%ao, ortmat))

            call diagonalise(fockmat)

            ! Transform coefficient into nonorthogonal AO basis
            fockmat%A = transpose(matmul(ortmat, fockmat%A))

            ! [TODO]: this is just hard-coded for now for water, URGENT: add geom_read_in
            call build_density(fockmat%A, density, 5)

            call update_scf_energy(density, fockmat%ao, st, int_store, conv)
            if (conv) then
               write(iunit, '(1X, A)') 'Convergence reached within tolerance.'
               write(iunit, '(1X, A, 1X, F15.8)') 'Final SCF Energy (Hartree):', st%energy
               write(iunit, '(1X, A)') 'Orbital energies (Hartree):'
               do i = size(fockmat%W), 1, -1
                  write(iunit, '(1X, I0, 1X, F15.8)') i, fockmat%W(i)
               end do
               exit
            else
               write(iunit, '(1X, A, 1X, I0, F15.8)') 'Iteration', iter, st%energy
            end if

            call build_fock(sys, density, fockmat%ao, int_store)
         end do

      end subroutine do_hartree_fock

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

      end subroutine diagonalise

      subroutine init_mat_t(matrix, mat)
         use read_in_integrals
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
         use read_in_integrals, only: int_store_t

         real(p), intent(in) :: density(:,:)
         real(p), intent(in) :: fock(:,:)
         type(int_store_t), intent(in) :: int_store
         type(state_t), intent(inout) :: st
         logical, intent(inout) :: conv

         st%density_old = st%density
         st%density = density
         
         st%energy_old = st%energy
         st%energy = sum(density * (int_store%core_hamil + fock))

         if (sqrt(sum((density-st%density_old)**2)) < 1e5*depsilon .and. abs(st%energy-st%energy_old) < 1e5*depsilon) conv = .true.

      end subroutine update_scf_energy
         
      subroutine build_fock(sys, density, fock, int_store)

         use read_in_integrals, only: int_store_t, eri_ind
         use system, only: system_t

         real(p), intent(in) :: density(:,:)
         type(int_store_t), intent(in) :: int_store
         type(system_t), intent(in) :: sys
         real(p), intent(out) :: fock(:,:)

         integer :: i, j, k, l

         associate(eri=>int_store%eri)
            do j = 1, sys%nbasis
               do i = 1, sys%nbasis
                  fock(i,j) = int_store%core_hamil(i,j)
                  do l = 1, sys%nbasis
                     do k = 1, sys%nbasis
                        ! [TODO]: This step can be sped up with pre-computed lookup arrays
                        fock(i,j) = fock(i,j) + density(k,l) * (2*eri(eri_ind(i,j,k,l)) - eri(eri_ind(i,k,j,l)))
                     end do
                  end do
               end do
            end do
         end associate

      end subroutine build_fock

end module hf
