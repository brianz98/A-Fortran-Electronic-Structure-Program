module ccsd

   use const
   use system, only: system_t

   implicit none

   type cc_amp_t
      ! Contains singles and doubles amplitudes
      real(p), allocatable :: t_ia(:,:)
      real(p), allocatable :: t_ijab(:,:,:,:)
   end type cc_amp_t

   type cc_int_t
      ! Intermediate quantities as defined in Stanton et al. (see do_ccsd below)
      ! Curly F tensors (Eqs. 3-5)
      real(dp), allocatable :: F_oo(:,:), F_vv(:,:), F_ov(:,:)
      ! Curly W tensors (Eqs. 6,8), see Appendix for why 7 is avoided
      real(dp), allocatable :: W_oooo(:,:,:,:), W_vvvv(:,:,:,:), W_ovvo(:,:,:,:)
      ! Effective two-particle excitation operators tau and tilde tau (Eqs. 9-10)
      real(dp), allocatable :: tau_tilde(:,:,:,:), tau(:,:,:,:)
      ! Energy denominators (Eqs. 12-13)
      real(dp), allocatable :: D_ia(:,:), D_ijab(:,:,:,:)
      real(dp), allocatable :: tmp_tia(:,:), tmp_tijab(:,:,:,:)
   end type cc_int_t

   type diis_cc_t
      ! Contains information for use in the CCSD DIIS procedure 
      ! Accelerating the convergence of the coupled-cluster approach: The use of the DIIS method
      ! Gustavo E. Scuseria, Timothy J. Lee, Henry F. Schaefer III
      ! (https://doi.org/10.1016/0009-2614(86)80461-4)

      ! These are identical as the HF-SCF ones, other than that we store and use the T matrices, 
      ! instead of Fock matrices for error matrix computation. Also that we don't store the error matrices as that will 
      ! double the memory usage what is a fast step (compared to updating amplitudes) anyway. Instead we use the idx lookup array
      ! to make sure we take the right differences in computing e_i matrices on the fly.
      !
      ! In more detail, as e_i = T_i+1 - T_i, and we overwrite T1 when T8 is the newest T matrix, we must make sure we never 
      ! subtract the newest T matrix with the oldest T matrix, and since doing that is tricky in a loop we keep the idx lookup array

      ! [TODO]: kick the slowest index to the last!
      integer :: n_errmat = 8
      integer :: n_active = 0
      integer :: iter = 0
      real(p), allocatable :: t1(:,:,:)
      real(p), allocatable :: t2(:,:,:,:,:)
      real(p), allocatable :: B(:,:)
      real(p), allocatable :: c(:)
      real(p), allocatable :: rhs(:)
      integer, allocatable :: idx(:,:)
   end type diis_cc_t
   
   contains

      subroutine do_ccsd(sys, int_store)
         ! This is the spinorbital formulation of J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett, J. Chem. Phys. volume 94, pp. 4334-4345 (1991) (https://doi.org/10.1063/1.460620)
         ! [TODO]: closed-shell (spatial) formulation of Hirata, S.; Podeszwa, R.; Tobita, M.; Bartlett, R. J. J. Chem. Phys. 2004, 120 (6), 2581 (https://doi.org/10.1063/1.1637577)

         use integrals, only: int_store_t, eri_ind
         use error_handling, only: check_allocate
         use system, only: state_t

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(inout) :: int_store
         integer :: i, j, a, b, p, q, r, s
         integer :: pr, qr, pa, pb, qa, qb, ra, rb, sa, sb
         real(dp) :: prqs, psqr

         type(state_t) :: st
         type(diis_cc_t) :: diis
         
         type(cc_amp_t) :: cc_amp
         integer :: ierr, iter
         integer, parameter :: iunit = 6
         integer, parameter :: maxiter = 100

         type(cc_int_t) :: cc_int
         logical :: conv

         write(iunit, '(1X, 10("-"))')
         write(iunit, '(1X, A)') 'CCSD'
         write(iunit, '(1X, 10("-"))')

         ! #############################################################################
         ! Form the antisymmetrised spinorbital basis ERIs: <pq||rs> = <pq|rs> - <pq|sr>
         ! where <pq|rs> = (pr|qs) * delta(sigma_p,sigma_r) * delta(sigma_q,sigma_s)
         ! #############################################################################
         write(iunit, '(1X, A)') 'Forming antisymmetrised spinorbital ERIs...'
         allocate(int_store%asym_spinorb(sys%nbasis*2,sys%nbasis*2,sys%nbasis*2,sys%nbasis*2), source=0.0_dp, stat=ierr)
         call check_allocate('int_store%asym_spinorb', sys%nbasis**4*16, ierr)

         call init_diis_cc_t(sys, diis)

         associate(n=>sys%nbasis, asym=>int_store%asym_spinorb, eri=>int_store%eri_mo)
         ! [TODO]: this isn't necessary but easier for debugging, change to spin lookup arrays or something similar
         do p = 1, n
            pa = 2*p-1
            pb = pa+1

            do q = 1, n
               qa = 2*q-1
               qb = qa+1

               do r = 1, n
                  ra = 2*r-1
                  rb = ra+1

                  pr = eri_ind(p,r)
                  qr = eri_ind(q,r)

                  do s = 1, n
                     sa = 2*s-1
                     sb = sa+1
                     prqs = eri(eri_ind(pr,eri_ind(q,s)))
                     psqr = eri(eri_ind(eri_ind(p,s),qr))
                     ! Draw a spin decision tree to visualise
                     asym(pa,qa,ra,sa) = prqs-psqr
                     asym(pb,qb,rb,sb) = prqs-psqr
                     asym(pa,qb,ra,sb) = prqs
                     asym(pb,qa,rb,sa) = prqs
                     asym(pa,qb,rb,sa) = -psqr
                     asym(pb,qa,ra,sb) = -psqr
                  end do
               end do
            end do
         end do
         end associate

         ! ############################################
         ! Initialise intermediate arrays and DIIS data
         ! ############################################
         call init_cc(sys, int_store, cc_amp, cc_int)
         allocate(st%t2_old, source=cc_amp%t_ijab)
         call init_diis_cc_t(sys, diis)

         write(iunit, '(1X, A)') 'Initialisation done, now entering iterative CC solver...'

         ! ###################
         ! Iterative CC solver
         ! ###################
         associate(nbasis=>sys%nbasis, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab, asym=>int_store%asym_spinorb)
         do iter = 1, maxiter
            ! Update intermediate tensors
            call build_tau(sys, cc_amp, cc_int)
            call build_F(sys, cc_int, cc_amp, asym)
            call build_W(sys, cc_int, cc_amp, asym)

            ! CC amplitude equations
            call update_amplitudes(sys, cc_amp, int_store, cc_int)        
            call update_cc_energy(sys, st, int_store, cc_amp, conv)
            if (conv) then
               write(iunit, '(1X, A)') 'Convergence reached within tolerance.'
               write(iunit, '(1X, A, 1X, F15.8)') 'Final CCSD Energy (Hartree):', st%energy

               ! Copied for return 
               sys%e_ccsd = st%energy
               call move_alloc(cc_amp%t_ia, sys%t1)
               call move_alloc(cc_amp%t_ijab, sys%t2)
               exit
            end if
            write(iunit, '(1X, A, 1X, I3, 1X, F15.8)') 'Iteration', iter, st%energy
            call update_diis_cc(diis, cc_amp)
         end do
         end associate

         ! Nonlazy deallocations
         deallocate(st%t2_old)
         call deallocate_cc_int_t(cc_int)
         call deallocate_diis_cc_t(diis)

      end subroutine do_ccsd

      subroutine init_cc(sys, int_store, cc_amp, cc_int)
         ! Initialise many CC related quantities.
         ! In:
         !     int_store: int_store_t object holding integrals.
         ! In/out:
         !     sys: system under study.
         !     cc_amp: CC amplitudes
         !     cc_int: CC intermediates being initialised

         use integrals, only: int_store_t
         use error_handling, only: check_allocate 

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(in) :: int_store
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         integer :: p, q, r, s, i, j, a, b, ierr
         integer, parameter :: iunit = 6

         write(iunit, '(1X, A)') 'Forming energy denominator matrices...'
         ! Forming the energy denominator matrices
         associate(n=>sys%nbasis, nocc=>sys%nocc, e=>sys%canon_levels)
         allocate(cc_int%D_ia(2*nocc, 2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%D_ia', 4*nocc*(n-nocc), ierr)
         allocate(cc_int%D_ijab(2*nocc, 2*nocc, 2*nocc+1:2*n, 2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%D_ijab', 16*(nocc**2*(n-nocc)**2), ierr)
         do i = 1, nocc
            do a = nocc+1, n
               cc_int%D_ia(2*i-1:2*i,2*a-1:2*a) = e(i)-e(a)
               do j = 1, nocc
                  do b = nocc+1, n
                     cc_int%D_ijab(2*i-1:2*i, 2*j-1:2*j, 2*a-1:2*a, 2*b-1:2*b) = e(i)+e(j)-e(a)-e(b)
                  end do
               end do
            end do
         end do

         ! This is only later useful in CCSD(T) but we can just do it here...
         allocate(sys%canon_levels_spinorb(2*n), source=0.0_dp)
         do i = 1, n
            sys%canon_levels_spinorb(2*i-1:2*i) = e(i)
         end do
      
         write(iunit, '(1X, A)') 'Allocating amplitude tensors...'
         ! Initialise t_i^a=0 and t_{ij}^{ab}=MP1 WF
         allocate(cc_amp%t_ijab(2*nocc,2*nocc,2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_amp%t_ijab', 16*(nocc**2*(n-nocc)**2), ierr)
         allocate(cc_int%tmp_tijab, mold=cc_amp%t_ijab, stat=ierr)
         call check_allocate('cc_int%tmp_tijab', 16*(nocc**2*(n-nocc)**2), ierr)
         allocate(cc_amp%t_ia(2*nocc, 2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_amp%t_ia', 4*nocc*(n-nocc), ierr)
         allocate(cc_int%tmp_tia, mold=cc_amp%t_ia, stat=ierr)
         call check_allocate('cc_int%tmp_tia', 4*nocc*(n-nocc), ierr)
         end associate

         write(iunit, '(1X, A)') 'Forming initial amplitude guesses...'
         associate(doubles=>cc_amp%t_ijab, nocc=>sys%nocc, asym=>int_store%asym_spinorb)
         do p = 1, 2*nocc
            do q = 1, 2*nocc
               do r = 2*sys%nocc+1, 2*sys%nbasis
                  do s = 2*sys%nocc+1, 2*sys%nbasis
                     doubles(p,q,r,s) = asym(p,q,r,s)/cc_int%D_ijab(p,q,r,s)
                  end do
               end do
            end do
         end do
         end associate

         write(iunit, '(1X, A)') 'Allocating intermediate tensors...'
         ! Allocate the intermediate tensors
         associate(n=>sys%nbasis, nocc=>sys%nocc)
         allocate(cc_int%F_vv(2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%F_vv', 4*(n-nocc)**2, ierr)
         allocate(cc_int%F_oo(2*nocc,2*nocc), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%F_oo', 4*nocc**2, ierr)
         allocate(cc_int%F_ov(2*nocc,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%F_ov', 4*nocc*(n-nocc), ierr)
         allocate(cc_int%W_oooo(2*nocc,2*nocc,2*nocc,2*nocc), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%W_oooo', 16*nocc**4, ierr)
         allocate(cc_int%W_vvvv(2*nocc+1:2*n,2*nocc+1:2*n,2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%W_vvvv', 16*(n-nocc)**4, ierr)
         allocate(cc_int%W_ovvo(2*nocc,2*nocc+1:2*n,2*nocc+1:2*n,2*nocc), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%W_ovvo', 16*nocc**2*(n-nocc)**2, ierr)
         allocate(cc_int%tau(2*nocc,2*nocc,2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%tau', 16*nocc**2*(n-nocc)**2, ierr)
         allocate(cc_int%tau_tilde(2*nocc,2*nocc,2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%tau_tilde', 16*nocc**2*(n-nocc)**2, ierr)
         end associate

      end subroutine init_cc

      subroutine deallocate_cc_int_t(cc_int)
         ! Nonlazy deallocation of the relatively large CC intermediate arrays.
         ! In/out:
         !     cc_int: cc_int_t object being deallocated.

         type(cc_int_t), intent(inout) :: cc_int

         deallocate(cc_int%F_oo, cc_int%F_vv, cc_int%F_ov)
         deallocate(cc_int%W_oooo, cc_int%W_vvvv, cc_int%W_ovvo)
         deallocate(cc_int%tau, cc_int%tau_tilde)
         deallocate(cc_int%D_ia, cc_int%D_ijab)
         deallocate(cc_int%tmp_tia, cc_int%tmp_tijab)
      end subroutine

      subroutine init_diis_cc_t(sys, diis)
         ! Initialise CCSD DIIS quantities. Some are left unallocated as they need to be resized during the first few iterations.
         ! In:
         !     sys: system under study.
         ! In/out:
         !     diis: diis_cc_t object being initialised.

         use error_handling, only: check_allocate

         type(system_t), intent(in) :: sys
         type(diis_cc_t), intent(out) :: diis

         integer :: ierr

         allocate(diis%idx(diis%n_errmat,2), source=0)

         associate(n=>sys%nbasis, nocc=>sys%nocc)
            allocate(diis%t1(diis%n_errmat,2*nocc,2*nocc+1:2*n), source=0.0_p, stat=ierr)
            call check_allocate('diis%t1', diis%n_errmat*4*nocc*(n-nocc), ierr)
            allocate(diis%t2(diis%n_errmat,2*nocc,2*nocc,2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_p, stat=ierr)
            call check_allocate('diis%t1', diis%n_errmat*16*(nocc*(n-nocc))**2, ierr)
         end associate
      end subroutine init_diis_cc_t

      subroutine update_diis_cc(diis, cc_amp)
         ! Update a DIIS iteration and form the new guess CC amplitudes (both singles and doubles).
         ! In/out:
         !     diis: the diis_cc_t object being updated
         !     cc_amp: CC amplitudes being updated

         use linalg, only: linsolve
         use error_handling, only: error

         type(diis_cc_t), intent(inout) :: diis
         type(cc_amp_t), intent(inout) :: cc_amp

         integer :: i, j, ierr

         diis%iter = diis%iter+1
         if (diis%iter > diis%n_errmat) diis%iter = diis%iter - diis%n_errmat
         if (diis%n_active < diis%n_errmat) diis%n_active = diis%n_active+1
         diis%t1(diis%iter,:,:) = 0.0_p
         diis%t2(diis%iter,:,:,:,:) = 0.0_p
         diis%t1(diis%iter,:,:) = cc_amp%t_ia
         diis%t2(diis%iter,:,:,:,:) = cc_amp%t_ijab

         associate(n=>diis%n_active, nerr=>diis%n_errmat, it=>diis%iter, idx=>diis%idx)
         if (n > 2) then
            ! We need at least three CC iterations to have at least two error matrices

            ! Construct index list (idx), see diis_cc_t documentation for motivation
            ! Concrete example, say n_errmat=8, and we're at iteration 4 (overwriting the 4th T matrices), so the error matrices
            ! should be: 
            ! e_1 = T4-T3 (== 4-3), e_2=3-2, e_3=2-1, e_4=1-8, e_5=8-7, e_6=7-6, e_7=6-5
            ! So our idx array should look like
            ! 4 3 2 1 8 7 6
            ! 3 2 1 8 7 6 5
            ! (note that the code below gives something like this but with the columns swapped around but that's ok)
            do i = 1, n-1
               idx(i,1) = mod(it+i+1, n)+1
               idx(i,2) = mod(it+i, n)+1
            end do
            ! Construct the B matrix
            if (n <= nerr) then
               if (allocated(diis%B)) deallocate(diis%B, diis%c, diis%rhs)
               allocate(diis%B(n+1,n+1), diis%c(n+1), diis%rhs(n+1), source=0.0_p)
            end if
            diis%B(n+1,:) = -1.0_p
            diis%B(n+1,n+1) = 0.0_p
            diis%rhs(n+1) = -1.0_p
            diis%c = diis%rhs
            do i = 1, n-1
               do j = 1, i
                  ! We compute the error matrices on the fly
                  ! [TODO] - shift iter index to the end to give better cache behaviour
                  diis%B(i,j) = sum((diis%t1(idx(i,1),:,:)-diis%t1(idx(i,2),:,:))*(diis%t1(idx(j,1),:,:)-diis%t1(idx(j,2),:,:))) &
                  + sum((diis%t2(idx(i,1),:,:,:,:)-diis%t2(idx(i,2),:,:,:,:))*(diis%t2(idx(j,1),:,:,:,:)-diis%t2(idx(j,2),:,:,:,:)))
               end do
            end do
            call linsolve(diis%B, diis%c, ierr)
            if (ierr /= 0) call error('ccsd::update_diis_cc', 'Linear solve failed!')
            cc_amp%t_ia = 0.0_p
            cc_amp%t_ijab = 0.0_p
            do i = 1, n
               cc_amp%t_ia = cc_amp%t_ia + diis%c(i) * diis%t1(i,:,:)
               cc_amp%t_ijab = cc_amp%t_ijab + diis%c(i) * diis%t2(i,:,:,:,:)
            end do
         end if
         end associate
      end subroutine update_diis_cc

      subroutine deallocate_diis_cc_t(diis)
         ! Nonlazy deallocation of the diis_cc_t object as they take up significant memory
         ! In/out:
         !     diis: the diis_cc_t object being deallocated

         type(diis_cc_t), intent(inout) :: diis

         deallocate(diis%t1, diis%t2, diis%B, diis%c, diis%rhs)
      end subroutine deallocate_diis_cc_t

      subroutine build_tau(sys, cc_amp, cc_int)
         ! Build the effective singles and doubles operators.
         ! In:
         !     sys: system under study.
         !     cc_amp: CC amplitudes.
         ! In/out:
         !     cc_int: CC intermediates with tau and tau_tilde updated

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         integer :: i, j, a, b
         real(p) :: ia, ja, x

         associate(n=>sys%nbasis, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab,&
            tau=>cc_int%tau, tau_tilde=>cc_int%tau_tilde)
            tau = 0.0_p
            tau_tilde = 0.0_p
            !$omp parallel do default(none) &
            !$omp private(i,j,a,b,ia,ja,x) &
            !$omp shared(sys, cc_amp, cc_int) &
            !$omp schedule(dynamic, 2) collapse(2)
            do b = 2*nocc+1, 2*n
               do a = 2*nocc+1, 2*n
                  do j = 1, 2*nocc
                     ja = t_ia(j,a)
                     do i = 1, 2*nocc
                        x = t_ia(i,a)*t_ia(j,b) - t_ia(i,b)*ja
                        tau_tilde(i,j,a,b) = t_ijab(i,j,a,b) + 0.5*x
                        tau(i,j,a,b) = tau_tilde(i,j,a,b) + 0.5*x
                     end do
                  end do
               end do
            end do
            !$omp end parallel do
         end associate

      end subroutine build_tau

      subroutine build_F(sys, cc_int, cc_amp, asym)
         ! Build the two-index F intermediates,
         ! Eqs. 3-5, note as we use HF reference all terms involving the fock matrix vanishes.
         ! In:
         !     sys: System under study.
         !     cc_amp: CC amplitudes
         !     asym: the antisymmetrised two-electron integrals <pq||rs>
         ! In/out:
         !     cc_int: CC intermediates with the F matrices updated.

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         real(p), intent(in) :: asym(:,:,:,:)
         type(cc_int_t), intent(inout) :: cc_int
         
         integer :: a, e, f, m, n, i

         associate(nbasis=>sys%nbasis, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab,&
            F_vv=>cc_int%F_vv, F_oo=>cc_int%F_oo, F_ov=>cc_int%F_ov, tau_tilde=>cc_int%tau_tilde)
         ! F_ae = \sum_{mf} t_m^f * <ma||fe> - 0.5 \sum_{mnf} \tau~_{mn}^{af} <mn||ef>
         F_vv = 0.0_p
         F_oo = 0.0_p
         F_ov = 0.0_p

         !$omp parallel default(none) &
         !$omp private(a, e, f, m, n, i) &
         !$omp shared(sys, cc_int, cc_amp, asym)

         !$omp do schedule(dynamic, 2) collapse(2)
         do a = 2*nocc+1, 2*nbasis
            do e = 2*nocc+1, 2*nbasis
               do m = 1, 2*nocc
                  do f = 2*nocc+1, 2*nbasis
                     F_vv(a,e) = F_vv(a,e) + t_ia(m,f) * asym(m,a,f,e)
                     do n = 1, 2*nocc
                        F_vv(a,e) = F_vv(a,e) - 0.5*tau_tilde(m,n,a,f)*asym(m,n,e,f)
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do

         ! F_mi = \sum_{en} t_n^e <mn||ie> + 0.5 \sum_{nef} \tau~_{in}^{ef} <mn||ef>
         !$omp do schedule(dynamic, 2) collapse(2)
         do m = 1, 2*nocc
            do i = 1, 2*nocc
               do e = 2*nocc+1, 2*nbasis
                  do n = 1, 2*nocc
                     F_oo(m,i) = F_oo(m,i) + t_ia(n,e)*asym(m,n,i,e)
                     do f = 2*nocc+1, 2*nbasis
                        F_oo(m,i) = F_oo(m,i) + 0.5*tau_tilde(i,n,e,f)*asym(m,n,e,f)
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do

         ! F_me = \sum_{nf} t_n^f * <mn||ef>
         !$omp do schedule(dynamic, 2) collapse(2)
         do m = 1, 2*nocc
            do e = 2*nocc+1, 2*nbasis
               do n = 1, 2*nocc
                  do f = 2*nocc+1, 2*nbasis
                     F_ov(m,e) = F_ov(m,e) + t_ia(n,f)*asym(m,n,e,f)
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel
         end associate
      end subroutine build_F

      subroutine build_W(sys, cc_int, cc_amp, asym)
         ! Build the four-index W intermediates, Eqs. 6 and 8.
         ! Note that we use a factor of 0.5 in the innermost line of the W_oooo loop instead of 0.25 per main text of Stanton,
         ! this is to make allowance for when W_vvvv cannot be stored due to memory limitation and has to be computed on the fly.
         ! See appendix of the same paper.
         ! In:
         !     sys: System under study.
         !     cc_amp: CC amplitudes
         !     asym: the antisymmetrised two-electron integrals <pq||rs>
         ! In/out:
         !     cc_int: CC intermediates with the F matrices updated.

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         real(p), intent(in) :: asym(:,:,:,:)
         type(cc_int_t), intent(inout) :: cc_int
         integer :: m, n, i, j, e, f, b, a
         real(p) :: x

         associate(nbasis=>sys%nbasis, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab,&
            W_oooo=>cc_int%W_oooo, W_ovvo=>cc_int%W_ovvo, tau=>cc_int%tau, W_vvvv=>cc_int%W_vvvv)
         W_oooo = 0.0_p
         W_ovvo = 0.0_p
         !$omp parallel default(none) &
         !$omp private(m, n, i, j, e, f, b, x, a) &
         !$omp shared(sys, cc_int, cc_amp, asym)

         !$omp do schedule(dynamic, 2) collapse(4)
         do m = 1, 2*nocc
            do n = 1, 2*nocc
               do i = 1, 2*nocc
                  do j = 1, 2*nocc
                     W_oooo(m,n,i,j) = asym(m,n,i,j)
                     do e = 2*nocc+1, 2*nbasis
                        W_oooo(m,n,i,j) = W_oooo(m,n,i,j) + t_ia(j,e)*asym(m,n,i,e) + t_ia(i,e)*asym(e,j,m,n)
                        do f = 2*nocc+1, 2*nbasis
                           W_oooo(m,n,i,j) = W_oooo(m,n,i,j) - 0.5*t_ijab(i,j,e,f)*asym(f,e,m,n)
                        end do
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do

         ! [TODO]: make a switch to on-the-fly W_vvvv computation, for when memory is limited
         !$omp do schedule(dynamic, 2) collapse(4)
         do f = 2*nocc+1, 2*nbasis
            do e = 2*nocc+1, 2*nbasis
               do b = 2*nocc+1, 2*nbasis
                  do a = 2*nocc+1, 2*nbasis
                     W_vvvv(a,b,e,f) = asym(a,b,e,f)
                     do n = 1, 2*nocc
                        W_vvvv(a,b,e,f) = W_vvvv(a,b,e,f) + t_ia(n,b)*asym(n,a,e,f) - t_ia(n,a)*asym(n,b,e,f)
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(dynamic, 2) collapse(4)
         do m = 1, 2*nocc
            do b = 2*nocc+1, 2*nbasis
               do e = 2*nocc+1, 2*nbasis
                  do j = 1, 2*nocc
                     W_ovvo(m,b,e,j) = asym(m,b,e,j)
                     do f = 2*nocc+1, 2*nbasis
                        W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) - t_ia(j,f)*asym(f,e,m,b)
                     end do

                     do n = 1, 2*nocc
                        x = t_ia(n,b)
                        W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) - x*asym(m,n,e,j)
                        do f = 2*nocc+1, 2*nbasis
                           W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) - (0.5*t_ijab(j,n,f,b) - x*t_ia(j,f)) * asym(f,e,m,n)
                        end do
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel

         end associate
      end subroutine build_W

      subroutine update_amplitudes(sys, cc_amp, int_store, cc_int)
         ! Perform the CC amplitude equations
         ! In:
         !     sys: system under study.
         !     int_store: integral information
         ! In/out:
         !     cc_amp: CC amplitudes being updated
         !     cc_int: CC intermediates used in computing the updates

         use integrals, only: int_store_t

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(in) :: int_store
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         integer :: i, j, a, b, m, n, e, f
         
         associate(tmp_t1=>cc_int%tmp_tia,tmp_t2=>cc_int%tmp_tijab,t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab,asym=>int_store%asym_spinorb, &
            F_oo=>cc_int%F_oo,F_ov=>cc_int%F_ov,F_vv=>cc_int%F_vv,D_ia=>cc_int%D_ia,D_ijab=>cc_int%D_ijab,W_oooo=>cc_int%W_oooo, &
            W_ovvo=>cc_int%W_ovvo,tau=>cc_int%tau,tau_tilde=>cc_int%tau_tilde,nocc=>sys%nocc,nbasis=>sys%nbasis,&
            W_vvvv=>cc_int%W_vvvv)

         ! Update T1
         !$omp parallel default(none) &
         !$omp private(i, j, a, b, m, n, e, f) &
         !$omp shared(sys, cc_amp, int_store, cc_int)

         ! Implied barrier here, better than master which doesn't have implied barrier
         ! tmp_t1 can be shared as each loop iteration will update a unique element of it.
         !$omp single
         tmp_t1 = 0.0_dp
         !$omp end single
         !$omp do schedule(dynamic, 2) collapse(2)
         do i = 1, 2*nocc
            do a = 2*nocc+1, 2*nbasis
               do m = 1, 2*nocc
                  tmp_t1(i,a) = tmp_t1(i,a) - t1(m,a)*F_oo(m,i)
                  do e = 2*nocc+1, 2*nbasis
                     tmp_t1(i,a) = tmp_t1(i,a) + t2(i,m,a,e)*F_ov(m,e)
                     do f = 2*nocc+1, 2*nbasis
                        tmp_t1(i,a) = tmp_t1(i,a) + 0.5*t2(i,m,e,f)*asym(f,e,m,a)
                     end do
                     do n = 1, 2*nocc
                        tmp_t1(i,a) = tmp_t1(i,a) - 0.5*t2(m,n,a,e)*asym(n,m,e,i)
                     end do
                  end do
               end do
               do e = 2*nocc+1, 2*nbasis
                  tmp_t1(i,a) = tmp_t1(i,a) + t1(i,e)*F_vv(a,e)
                  do n = 1, 2*nocc
                     tmp_t1(i,a) = tmp_t1(i,a) - t1(n,e)*asym(n,a,i,e)
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp single
         tmp_t1 = tmp_t1/D_ia
         
         ! Update T2
         tmp_t2 = 0.0_dp
         !$omp end single

         !$omp do schedule(dynamic, 2) collapse(4)
         do i = 1, 2*nocc
            do j = 1, 2*nocc
               do a = 2*nocc+1, 2*nbasis
                  do b = 2*nocc+1, 2*nbasis
                     tmp_t2(i,j,a,b) = asym(i,j,a,b)
                     do e = 2*nocc+1, 2*nbasis
                        tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) + t2(i,j,a,e)*F_vv(b,e) - t2(i,j,b,e)*F_vv(a,e) &
                        + tmp_t1(i,e)*asym(a,b,e,j) - tmp_t1(j,e)*asym(a,b,e,i)
                        do m = 1, 2*nocc
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - 0.5*F_ov(m,e)*(t2(i,j,a,e)*tmp_t1(m,b)-&
                              t2(i,j,b,e)*tmp_t1(m,a))
                        end do
                        do f = 2*nocc+1, 2*nbasis
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) + 0.5*tau(i,j,e,f)*W_vvvv(a,b,e,f)
                        end do
                     end do
                     do m = 1, 2*nocc
                        tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - t2(i,m,a,b)*F_oo(m,j) + t2(j,m,a,b)*F_oo(m,i) &
                        - tmp_t1(m,a)*asym(m,b,i,j) + tmp_t1(m,b)*asym(m,a,i,j)
                        do e = 2*nocc+1, 2*nbasis
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - 0.5*F_ov(m,e)*(t2(i,m,a,b)*tmp_t1(j,e)-&
                              t2(j,m,a,b)*tmp_t1(i,e))
                        end do
                        do n = 1, 2*nocc
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) + 0.5*tau(m,n,a,b)*W_oooo(m,n,i,j)
                        end do
                     end do

                     do m = 1, 2*nocc
                        do e = 2*nocc+1, 2*nbasis
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) &
                           + t2(i,m,a,e)*W_ovvo(m,b,e,j) - tmp_t1(i,e)*tmp_t1(m,a)*asym(e,j,m,b)&
                           - t2(j,m,a,e)*W_ovvo(m,b,e,i) + tmp_t1(j,e)*tmp_t1(m,a)*asym(e,i,m,b)&
                           - t2(i,m,b,e)*W_ovvo(m,a,e,j) + tmp_t1(i,e)*tmp_t1(m,b)*asym(e,j,m,a)&
                           + t2(j,m,b,e)*W_ovvo(m,a,e,i) - tmp_t1(j,e)*tmp_t1(m,b)*asym(e,i,m,a)
                        end do 
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel
         t1 = tmp_t1
         t2 = tmp_t2/D_ijab
         end associate

      end subroutine update_amplitudes

      subroutine update_cc_energy(sys, st, int_store, cc_amp, conv)
         ! Updates the CC energy and also checks for convergence
         ! In:
         !     sys: system under study.
         !     int_store: integral information.
         !     cc_amp: CC amplitudes.
         ! In/out:
         !     st: state_t object storing energies and RMS(T) info for convergence checking
         ! Out:
         !     conv: returns true if converged within tolerance

         use system, only: state_t
         use integrals, only: int_store_t

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(in) :: int_store
         type(cc_amp_t), intent(in) :: cc_amp
         type(state_t), intent(inout) :: st
         logical, intent(out) :: conv

         integer :: i, j, a, b
         real(p) :: ecc, rmst2

         conv = .false.
         st%energy_old = st%energy
         associate(nocc=>sys%nocc, nbasis=>sys%nbasis, asym=>int_store%asym_spinorb, t1=>cc_amp%t_ia, t2=>cc_amp%t_ijab)
            !$omp parallel default(none) &
            !$omp private(i,j,a,b) &
            !$omp shared(sys,st,int_store,cc_amp,ecc,rmst2) 
            ecc = 0.0_p
            rmst2 = 0.0_p
            !$omp do schedule(dynamic, 2) collapse(4) reduction(+:ecc,rmst2)
            do b = 2*nocc+1, 2*nbasis
               do a = 2*nocc+1, 2*nbasis
                  do j = 1, 2*nocc
                     do i = 1, 2*nocc
                        ecc = ecc + 0.25*asym(i,j,a,b)*(t2(i,j,a,b)+2*t1(i,a)*t1(j,b))
                        ! We only do the RMS on T2 amplitudes
                        rmst2 = rmst2 + (t2(i,j,a,b)-st%t2_old(i,j,a,b))**2
                     end do
                  end do
               end do
            end do
            !$omp end do
            !$omp end parallel
            st%energy = ecc
            st%t2_old = t2
            ! [TODO]: stop hard-coding tolerance, or at least unify it with HF
            if (sqrt(rmst2) < sys%ccsd_t_tol .and. abs(st%energy-st%energy_old) < sys%ccsd_e_tol) conv = .true.

         end associate

      end subroutine update_cc_energy

      subroutine do_ccsd_t(sys, int_store)
         ! Spinorbital formuation of CCSD(T)
         ! In:
         !     int_store: integral information.
         ! In/out:
         !     sys: holds converged amplitudes from CCSD

         use integrals, only: int_store_t
         use error_handling, only: check_allocate

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(in) :: int_store

         integer :: i, j, k, a, b, c, f, m, ierr
         real(p), allocatable :: tmp_t3d(:,:,:), tmp_t3c(:,:,:), tmp_t3c_d(:,:,:)
         real(p) :: e_T
         integer, parameter :: iunit = 6

         write(iunit, '(1X, 10("-"))')
         write(iunit, '(1X, A)') 'CCSD(T)'
         write(iunit, '(1X, 10("-"))')

         associate(nbasis=>sys%nbasis, nocc=>sys%nocc, e=>sys%canon_levels_spinorb, t1=>sys%t1, t2=>sys%t2,&
                   asym=>int_store%asym_spinorb)
         !$omp parallel default(none) &
         !$omp private(i, j, k, a, b, c, f, m, tmp_t3d, tmp_t3c, tmp_t3c_d, ierr) &
         !$omp shared(sys, int_store, e_T)
         e_T = 0.0_p
         
         ! Declare heap-allocated arrays inside parallel region
         allocate(tmp_t3d(2*nocc+1:2*nbasis,2*nocc+1:2*nbasis,2*nocc+1:2*nbasis), source=0.0_p, stat=ierr)
         call check_allocate('tmp_t3d', 8*(nbasis-nocc)**3, ierr)
         allocate(tmp_t3c(2*nocc+1:2*nbasis,2*nocc+1:2*nbasis,2*nocc+1:2*nbasis), source=0.0_p, stat=ierr)
         call check_allocate('tmp_t3c', 8*(nbasis-nocc)**3, ierr)
         allocate(tmp_t3c_d(2*nocc+1:2*nbasis,2*nocc+1:2*nbasis,2*nocc+1:2*nbasis), source=0.0_p, stat=ierr)
         call check_allocate('tmp_t3c_d', 8*(nbasis-nocc)**3, ierr)

         ! We use avoid the storage of full triples (six-dimensional array..) and instead use the strategy of 
         ! batched triples storage, denoted W^{ijk}(abc) in (https://doi.org/10.1016/0009-2614(91)87003-T).
         ! We could of course only compute one element in a loop iteration but that will probably result in bad floating point
         ! performance, here for each thread/loop iteration we use the Fortran intrinsic sum, 
         ! instead of all relying on OMP reduction, to hopefully give better floating point performance.
         !$omp do schedule(dynamic, 2) collapse(3) reduction(+:e_T)
         do i = 1, 2*nocc
            do j = 1, 2*nocc
               do k = 1, 2*nocc
                  ! Compute T3 amplitudes                  
                  do a = 2*nocc+1, 2*nbasis
                     do b = 2*nocc+1, 2*nbasis
                        do c = 2*nocc+1, 2*nbasis
                           ! Disonnected T3: D_{ijk}^{abc}*t_{ijk}^{abc}(d) = P(i/jk)P(a/bc)t_i^a*<jk||bc>
                           ! Can be rewritten directly as no sums in the expression
                           tmp_t3d(a,b,c) = - t1(i,a)*asym(c,b,j,k) + t1(j,a)*asym(c,b,i,k) + t1(k,a)*asym(c,b,j,i) &
                                            + t1(i,b)*asym(c,a,j,k) - t1(j,b)*asym(c,a,i,k) - t1(k,b)*asym(c,a,j,i) &
                                            - t1(i,c)*asym(j,k,b,a) + t1(j,c)*asym(i,k,b,a) + t1(k,c)*asym(j,i,b,a)
                           tmp_t3d(a,b,c) = tmp_t3d(a,b,c)/(e(i)+e(j)+e(k)-e(a)-e(b)-e(c))               

                           ! Connected T3: D_{ijk}^{abc}*t_{ijk}^{abc}(c) = P(i/jk)P(a/bc)[\sum_f t_jk^af <fi||bc> - \sum_m t_im^bc <ma||jk>]
                           ! Has to be zeroed as all terms are sums
                           tmp_t3c(a,b,c) = 0.0_p
                           do f = 2*nocc+1, 2*nbasis
                              tmp_t3c(a,b,c) = tmp_t3c(a,b,c) &
                                             + t2(j,k,a,f)*asym(f,i,b,c) - t2(i,k,a,f)*asym(f,j,b,c) - t2(j,i,a,f)*asym(f,k,b,c) &
                                             - t2(j,k,b,f)*asym(f,i,a,c) + t2(i,k,b,f)*asym(f,j,a,c) + t2(j,i,b,f)*asym(f,k,a,c) &
                                             - t2(j,k,c,f)*asym(f,i,b,a) + t2(i,k,c,f)*asym(f,j,b,a) + t2(j,i,c,f)*asym(f,k,b,a)
                           end do
                           do m = 1, 2*nocc
                              tmp_t3c(a,b,c) = tmp_t3c(a,b,c) &
                                             - t2(i,m,b,c)*asym(m,a,j,k) + t2(j,m,b,c)*asym(m,a,i,k) + t2(k,m,b,c)*asym(m,a,j,i) &
                                             + t2(i,m,a,c)*asym(m,b,j,k) - t2(j,m,a,c)*asym(m,b,i,k) - t2(k,m,a,c)*asym(m,b,j,i) &
                                             + t2(i,m,b,a)*asym(m,c,j,k) - t2(j,m,b,a)*asym(m,c,i,k) - t2(k,m,b,a)*asym(m,c,j,i)
                           end do
                           tmp_t3c_d(a,b,c) = tmp_t3c(a,b,c)/(e(i)+e(j)+e(k)-e(a)-e(b)-e(c))
                        end do
                     end do
                  end do
                  
                  ! Calculate contributions
                  e_T = e_T + sum(tmp_t3c*(tmp_t3c_d+tmp_t3d))/36
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel
         sys%e_ccsd_t = e_T
         write(iunit, '(1X, A, 1X, F15.9)') 'CCSD(T) correlation energy (Hartree):', e_T
         end associate
      end subroutine do_ccsd_t
end module ccsd