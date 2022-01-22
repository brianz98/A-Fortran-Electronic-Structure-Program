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
      real(dp), dimension(:,:), allocatable :: F_oo, F_vv, F_ov
      ! Curly W tensors (Eqs. 6,8), see Appendix for why 7 is avoided
      real(dp), dimension(:,:,:,:), allocatable :: W_oooo, W_vvvv, W_ovvo
      ! Effective two-particle excitation operators tau and tilde tau (Eqs. 9-10)
      real(dp), dimension(:,:,:,:), allocatable :: tau_tilde, tau
      ! Energy denominators (Eqs. 12-13)
      real(dp), allocatable :: D_ia(:,:), D_ijab(:,:,:,:)
      real(dp), allocatable :: tmp_tia(:,:), tmp_tijab(:,:,:,:), tmp_tia_scratch(:,:), tmp_tijab_scratch(:,:,:,:)
      
      ! These are the quantities required by the spin-free formulation, left unallocated if using spinorbital formulation
      real(dp), dimension(:,:,:,:), allocatable :: v_div_d, v_oovv, v_vovo, v_vovv, v_oovo, v_oooo, v_vvvv
      real(dp), dimension(:,:), allocatable :: I_vo, I_vv, I_oo, I_oo_p
      real(dp), dimension(:,:,:,:), allocatable :: c_oovv, x_voov
      real(dp), dimension(:,:,:,:), allocatable :: I_oooo, I_ovvv, I_voov, I_vovv_p, I_ooov_p
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

      subroutine do_ccsd_spinorb(sys, int_store)
         ! This is the spinorbital formulation of J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett, J. Chem. Phys. volume 94, pp. 4334-4345 (1991) (https://doi.org/10.1063/1.460620)
         ! [TODO]: closed-shell (spatial) formulation of Hirata, S.; Podeszwa, R.; Tobita, M.; Bartlett, R. J. J. Chem. Phys. 2004, 120 (6), 2581 (https://doi.org/10.1063/1.1637577)

         use integrals, only: int_store_t, eri_ind
         use error_handling, only: error, check_allocate
         use system, only: state_t

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(inout) :: int_store
         integer :: i, j, a, b, p, q, r, s
         integer :: pr, qr, pa, pb, qa, qb, ra, rb, sa, sb
         real(dp) :: prqs, psqr
         real(dp) :: err
         integer(kind=8) :: c_max, c_rate, t0, t1

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
         call system_clock(count=t0, count_rate=c_rate, count_max=c_max)
         allocate(int_store%asym_spinorb(sys%nbasis*2,sys%nbasis*2,sys%nbasis*2,sys%nbasis*2), source=0.0_dp, stat=ierr)
         call check_allocate('int_store%asym_spinorb', sys%nbasis**4*16, ierr)

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
         call system_clock(t1)
         if (t1<t0) t1 = t1+c_max
         write(iunit, '(1X, A, 1X, F8.6, A)') 'Time taken:', real(t1-t0, kind=dp)/c_rate, " s"
         write(iunit, *)
         t0=t1         

         write(iunit, '(1X, A)') 'Checking that the permuational symmetry of the antisymmetrised integrals hold...'
         err = 0.0_dp
         associate(nbasis=>sys%nbasis, asym=>int_store%asym_spinorb)
         do p = 1, 2*nbasis
            do q = 1, p
               do r = 1, p
                  do s = 1, min(r,p)
                     err = err + abs(asym(p,q,r,s) + asym(p,q,s,r)) + abs(asym(p,q,r,s) - asym(r,s,p,q)) &
                               + abs(asym(p,q,r,s) + asym(s,r,p,q)) + abs(asym(p,q,r,s) - asym(s,r,q,p))
                  end do
               end do
            end do
         end do
         end associate

         if (err > depsilon) then
            write(iunit, '(1X, A, 1X, E15.6)') 'Permutational symmetry error:', err
            call error('ccsd::do_ccsd', 'Permutational symmetry of antisymmetrised integrals does not hold')
         end if

         call system_clock(t1)
         if (t1<t0) t1 = t1+c_max
         write(iunit, '(1X, A, 1X, F8.6, A)') 'Time taken:', real(t1-t0, kind=dp)/c_rate, " s"
         write(iunit, *)
         t0=t1

         ! ############################################
         ! Initialise intermediate arrays and DIIS data
         ! ############################################
         write(iunit, '(1X, A)') 'Initialise CC intermediate tensors and DIIS auxilliary arrays...'
         call init_cc(sys, int_store, cc_amp, cc_int, restricted=.false.)
         allocate(st%t2_old, source=cc_amp%t_ijab)
         call init_diis_cc_t(sys, diis)

         call system_clock(t1)
         if (t1<t0) t1 = t1+c_max
         write(iunit, '(1X, A, 1X, F8.6, A)') 'Time taken:', real(t1-t0, kind=dp)/c_rate, " s"
         write(iunit, *)
         t0=t1

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
            call update_diis_cc(diis, cc_amp)

            call system_clock(t1)
            write(iunit, '(1X, A, 1X, I2, 2X, F10.7, 2X, F8.6, 1X, A)') 'Iteration', iter, st%energy,real(t1-t0, kind=dp)/c_rate,'s'
            t0=t1
         end do
         end associate

         ! Nonlazy deallocations
         deallocate(st%t2_old)
         deallocate(cc_int, diis)
         !call deallocate_cc_int_t(cc_int)
         !call deallocate_diis_cc_t(diis)

      end subroutine do_ccsd_spinorb

      subroutine do_ccsd_spatial(sys, int_store)
         ! This is the spin-free (spatial) formulation of P. Piecuch et al., Computer Physics Communications 149 (2002) 71–96, 
         ! https://doi.org/10.1016/S0010-4655(02)00598-2

         use integrals, only: int_store_t, eri_ind
         use error_handling, only: error, check_allocate
         use system, only: state_t

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(inout) :: int_store
         integer :: i, j, a, b, p, q, r, s
         integer :: pr, qr, pa, pb, qa, qb, ra, rb, sa, sb
         real(dp) :: prqs, psqr
         real(dp) :: err
         integer(kind=8) :: c_max, c_rate, t0, t1

         type(state_t) :: st
         type(diis_cc_t) :: diis
         
         type(cc_amp_t) :: cc_amp
         integer :: ierr, iter
         integer, parameter :: iunit = 6
         integer, parameter :: maxiter = 100

         type(cc_int_t) :: cc_int
         logical :: conv

         write(iunit, '(1X, 14("-"))')
         write(iunit, '(1X, A)') 'Spin-free CCSD'
         write(iunit, '(1X, 14("-"))')

         ! ############################################
         ! Initialise intermediate arrays and DIIS data
         ! ############################################
         write(iunit, '(1X, A)') 'Initialise CC intermediate tensors and DIIS auxilliary arrays...'
         call init_cc(sys, int_store, cc_amp, cc_int, restricted=.true.)
         allocate(st%t2_old, source=cc_amp%t_ijab)
         call init_diis_cc_t(sys, diis)

         call system_clock(t1)
         if (t1<t0) t1 = t1+c_max
         write(iunit, '(1X, A, 1X, F8.6, A)') 'Time taken:', real(t1-t0, kind=dp)/c_rate, " s"
         write(iunit, *)
         t0=t1

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
            call update_diis_cc(diis, cc_amp)

            call system_clock(t1)
            write(iunit, '(1X, A, 1X, I2, 2X, F10.7, 2X, F8.6, 1X, A)') 'Iteration', iter, st%energy,real(t1-t0, kind=dp)/c_rate,'s'
            t0=t1
         end do
         end associate

         ! Nonlazy deallocations
         deallocate(st%t2_old)
         deallocate(cc_int, diis)
         !call deallocate_cc_int_t(cc_int)
         !call deallocate_diis_cc_t(diis)


      end subroutine do_ccsd_spatial

      subroutine init_cc(sys, int_store, cc_amp, cc_int, restricted)
         ! Initialise many CC related quantities.
         ! In:
         !     int_store: int_store_t object holding integrals.
         !     restricted: true if we're using the spin-free formulation
         ! In/out:
         !     sys: system under study.
         !     cc_amp: CC amplitudes
         !     cc_int: CC intermediates being initialised

         use integrals, only: int_store_t
         use error_handling, only: check_allocate 
         use integrals, only: eri_ind

         type(int_store_t), intent(in) :: int_store
         logical, intent(in) :: restricted
         type(system_t), intent(inout) :: sys
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         real(dp), allocatable :: temperi

         integer :: p, q, r, s, i, j, a, b, ierr
         integer, parameter :: iunit = 6

         ! Set loop limits
         if (restricted) then
            sys%nou = sys%nocc
            sys%nvl = sys%nou+1
            sys%nvu = sys%nbasis
         else
            sys%nou = 2*sys%nocc
            sys%nvl = 2*nocc+1
            sys%nvu = 2*n
         end if
         sys%nv = sys%nvu-sys%nvl+1
         sys%no_x_nv = sys%nou*sys%nvirt


         write(iunit, '(1X, A)') 'Forming energy denominator matrices...'
         ! Forming the energy denominator matrices D^{abc...}_{ijk...} = e_i + e_j + ... - e_a - e_b - ...
         associate(e=>sys%canon_levels, nocc=>sys%nou, nvl=>sys%nvl, nvu=>sys%nvu, nvirt=>sys%nv, no_x_nv=>sys%no_x_nv, &
            doubles=>cc_amp%t_ijab, asym=>int_store%asym_spinorb)
         allocate(cc_int%D_ia(nocc, nvl:nvu), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%D_ia', no_x_nv, ierr)
         allocate(cc_int%D_ijab(nocc, nocc, nvl:nvu, nvl:nvu), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%D_ijab', no_x_nv**2, ierr)
         if (restricted) then
            do i = 1, nocc
               do a = nocc+1, n
                  cc_int%D_ia(i,a) = e(i)-e(a)
                  do j = 1, nocc
                     do b = nocc+1, n
                        cc_int%D_ijab(i,j,a,b) = e(i)+e(j)-e(a)-e(b)
                     end do
                  end do
               end do
            end do
         else
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
         end if

         ! This is only later useful in CCSD(T) but we can just do it here...
         allocate(sys%canon_levels_spinorb(2*n), source=0.0_dp)
         do i = 1, n
            sys%canon_levels_spinorb(2*i-1:2*i) = e(i)
         end do
         

         write(iunit, '(1X, A)') 'Allocating amplitude tensors...'
         ! Initialise t_i^a=0 and t_{ij}^{ab}=MP1 WF
         allocate(cc_amp%t_ijab(nocc,nocc,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
         call check_allocate('cc_amp%t_ijab', no_x_nv**2), ierr)
         allocate(cc_int%tmp_tijab, cc_int%tmp_ijab_scratch, mold=cc_amp%t_ijab, stat=ierr)
         call check_allocate('cc_int%tmp_tijab', no_x_nv**2), ierr)
         call check_allocate('cc_int%tmp_ijab_scratch', no_x_nv**2), ierr)
         allocate(cc_amp%t_ia(nocc, nvl:nvu), source=0.0_dp, stat=ierr)
         call check_allocate('cc_amp%t_ia', no_x_nv, ierr)
         allocate(cc_int%tmp_tia, cc_int%tmp_tia_scratch, mold=cc_amp%t_ia, stat=ierr)
         call check_allocate('cc_int%tmp_tia', no_x_nv, ierr)
         call check_allocate('cc_int%tmp_tia_scratch', no_x_nv, ierr)

         if (restricted) then
            ! Re-sort MO ERIs into different slices for efficient vectorised operations
            allocate(cc_int%v_div_d(nocc,nocc,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_div_d', no_x_nv**2), ierr)
            allocate(cc_int%v_abij(nocc,nocc,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_oovv', no_x_nv**2), ierr)
            allocate(cc_int%v_aibj(nocc,nocc,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%bjai', no_x_nv**2), ierr)
            allocate(cc_int%v_abci(nvl:nvu,nocc,nvl:nvu,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_vovv', no_x_nv**2), ierr)
            allocate(cc_int%v_aijk(nvl:nvu,nocc,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_oovo', no_x_nv*nvirt**2), ierr)
            allocate(cc_int%v_ijkl(nocc,nocc,nocc,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_oooo', nocc**4), ierr)
            allocate(cc_int%v_abcd(nvl:nvu,nvl:nvu,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_vvvv', nvirt**4), ierr)

            allocate(temperi(nvu,nvu,nvu,nvu), source=0.0_dp, stat=ierr)
            call check_allocate('temperi', nvu**4, ierr)
            
            do s = 1, nvu
               do r = 1, nvu
                  rs = eri_ind(r,s)
                  do q = 1, nvu
                     do p = 1, nvu
                        temperi(p,q,r,s) = int_store%eri_mo(eri_ind(eri_ind(p,q),rs))
                     end do
                  end do
               end do
            end do

            cc_int%v_oovv = temperi(1:nocc,1:nocc,nvl:nvu,nvl:nvu)
            cc_int%v_div_d = cc_int%v_oovv / D_ijab
            cc_int%v_vovo = temperi(nvl:nvu,1:nocc,nvl:nvu,1:nocc)
            cc_int%v_vovv = temperi(nvl:nvu,1:nocc,nvl:nvu,nvl:nvu)
            cc_int%v_oovo = temperi(1:nocc,1:nocc,nvl:nvu,nvl:nvu)
            cc_int%v_oooo = temperi(1:nocc,1:nocc,1:nocc,1:nocc)
            cc_int%v_vvvv = temperi(nvl:nvu,nvl:nvu,nvl:nvu,nvl:nvu)

            deallocate(temperi)

         end if

         write(iunit, '(1X, A)') 'Forming initial amplitude guesses...'
         if (restricted) then
            doubles = v_div_d
         else
            doubles = asym(1:nocc,1:nocc,nvl:nvu,nvl:nvu)/cc_int%D_ijab
         end if

         write(iunit, '(1X, A)') 'Allocating intermediate tensors...'
         ! Allocate the intermediate tensors
         if (restricted) then
            ! Two-index intermediates
            allocate(cc_int%I_vo(nvl:nvu, nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_vo', no_x_nv, ierr)
            allocate(cc_int%I_vv(nvl:nvu, nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_vv', (nvu-nvl+1)**2, ierr)
            allocate(cc_int%I_oo_p(nocc, nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_oo_p', nocc**2, ierr)
            allocate(cc_int%I_oo(nocc, nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_oo', nocc**2, ierr)

            ! Four-index intermediates
            allocate(cc_int%c_oovv(nocc,nocc,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%c_oovv', no_x_nv**2, ierr)
            allocate(cc_int%x_voov(nvl:nvu,nocc,nocc,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%x_voov', no_x_nv**2, ierr)
            allocate(cc_int%I_oooo(nocc,nocc,nocc,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_oooo', nocc**4, ierr)
            allocate(cc_int%I_ovvv(nocc,nvl:nvu,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_ovvv', nocc*nvirt**3, ierr)
            allocate(cc_int%I_voov(nvl:nvu,nocc,nocc,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_voov', no_x_nv**2, ierr)
            allocate(cc_int%I_vovv_p(nvl:nvu,nocc,nocc,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_vovv_p', nvirt*nocc**3, ierr)
            allocate(cc_int%I_ooov_p(nocc,nocc,nocc,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_ooov_p', nvirt*nocc**3, ierr)

         else
            ! Two-index intermediates
            allocate(cc_int%F_vv(nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%F_vv', nvirt**2, ierr)
            allocate(cc_int%F_oo(nocc,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%F_oo', nocc**2, ierr)
            allocate(cc_int%F_ov(nocc,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%F_ov', no_x_nv, ierr)

            ! Four-index intermediates
            allocate(cc_int%W_oooo(nocc,nocc,nocc,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%W_oooo', nocc**4, ierr)
            allocate(cc_int%W_vvvv(nvl:nvu,nvl:nvu,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%W_vvvv', nvirt**4, ierr)
            allocate(cc_int%W_ovvo(nocc,nvl:nvu,nvl:nvu,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%W_ovvo', no_x_nv**2, ierr)

            ! Tau quantities
            allocate(cc_int%tau(nocc,nocc,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%tau', no_x_nv**2, ierr)
            allocate(cc_int%tau_tilde(nocc,nocc,nvl:nvu,nvl:nvu), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%tau_tilde', no_x_nv**2, ierr)
         end if
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
            !$omp schedule(runtime) collapse(2)
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

         !$omp do schedule(runtime) collapse(2)
         do a = 2*nocc+1, 2*nbasis
            do e = 2*nocc+1, 2*nbasis
               do m = 1, 2*nocc
                  do f = 2*nocc+1, 2*nbasis
                     F_vv(a,e) = F_vv(a,e) + t_ia(m,f) * asym(f,e,m,a)
                     do n = 1, 2*nocc
                        F_vv(a,e) = F_vv(a,e) + 0.5*tau_tilde(n,m,f,a)*asym(n,m,e,f)
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do

         ! F_mi = \sum_{en} t_n^e <mn||ie> + 0.5 \sum_{nef} \tau~_{in}^{ef} <mn||ef>
         !$omp do schedule(runtime) collapse(2)
         do m = 1, 2*nocc
            do i = 1, 2*nocc
               do e = 2*nocc+1, 2*nbasis
                  do n = 1, 2*nocc
                     F_oo(m,i) = F_oo(m,i) - t_ia(n,e)*asym(n,m,i,e)
                     do f = 2*nocc+1, 2*nbasis
                        F_oo(m,i) = F_oo(m,i) - 0.5*tau_tilde(n,i,f,e)*asym(f,e,m,n)
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do

         ! F_me = \sum_{nf} t_n^f * <mn||ef>
         !$omp do schedule(runtime) collapse(2)
         do m = 1, 2*nocc
            do e = 2*nocc+1, 2*nbasis
               do n = 1, 2*nocc
                  do f = 2*nocc+1, 2*nbasis
                     F_ov(m,e) = F_ov(m,e) + t_ia(n,f)*asym(f,e,n,m)
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

         !$omp do schedule(runtime) collapse(4)
         do m = 1, 2*nocc
            do n = 1, 2*nocc
               do i = 1, 2*nocc
                  do j = 1, 2*nocc
                     W_oooo(m,n,i,j) = asym(m,n,i,j)
                     do e = 2*nocc+1, 2*nbasis
                        W_oooo(m,n,i,j) = W_oooo(m,n,i,j) + t_ia(j,e)*asym(e,i,n,m) - t_ia(i,e)*asym(e,j,n,m)
                        do f = 2*nocc+1, 2*nbasis
                           W_oooo(m,n,i,j) = W_oooo(m,n,i,j) - 0.5*t_ijab(j,i,f,e)*asym(f,e,m,n)
                        end do
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do

         ! [TODO]: make a switch to on-the-fly W_vvvv computation, for when memory is limited
         !$omp do schedule(runtime) collapse(4)
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

         !$omp do schedule(runtime) collapse(4)
         do j = 1, 2*nocc
            do e = 2*nocc+1, 2*nbasis
               do b = 2*nocc+1, 2*nbasis
                  do m = 1, 2*nocc
                     W_ovvo(m,b,e,j) = asym(m,b,e,j)
                     do f = 2*nocc+1, 2*nbasis
                        W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) - t_ia(j,f)*asym(f,e,m,b)
                     end do

                     do n = 1, 2*nocc
                        x = t_ia(n,b)
                        W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) + x*asym(n,m,e,j)
                        do f = 2*nocc+1, 2*nbasis
                           W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) - (0.5*t_ijab(j,n,f,b) + x*t_ia(j,f)) * asym(f,e,n,m)
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

      subroutine update_amplitudes_restricted(sys, cc_amp, int_store, cc_int)
         ! Perform the CC amplitude equations for the spin-free formulation
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

         real(dp), dimension(:,:,:,:), allocatable :: reshape_tmp1, reshape_tmp2

         associate(tmp_t1=>cc_int%tmp_tia,tmp_t2=>cc_int%tmp_tijab,t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab,asym=>int_store%asym_spinorb, &
            F_oo=>cc_int%F_oo,F_ov=>cc_int%F_ov,F_vv=>cc_int%F_vv,D_ia=>cc_int%D_ia,D_ijab=>cc_int%D_ijab,W_oooo=>cc_int%W_oooo, &
            W_ovvo=>cc_int%W_ovvo,tau=>cc_int%tau,tau_tilde=>cc_int%tau_tilde,nocc=>sys%nocc,nbasis=>sys%nbasis,&
            W_vvvv=>cc_int%W_vvvv)

         ! Update T1
         
         ! t_i^e * I_e^a, a simple matmul
         ! Benchmarking suggests that the Fortran intrinsic matmul beats threaded dgemm, not to mention naive OpenMP. So here we go
         tmp_t1 = matmul(t1, I_vv)

         ! -I'_i^m * t_m^a, another simple matmul
         tmp_t1 = tmp_t1 - matmul(I_oo_p, t1)

         ! I_e^m ( 2 t_mi^ea - t_im^ea ) + t_m^e (2 v_ei^ma - v_ei%am)
         !$omp parallel do default(none)&
         !$omp schedule(static,10) collapse(2)&
         !$omp shared(nvl, nvu, nocc, tmp_t1, I_vo, t2, v_oovv, v_vovo)
         do a = nvl, nvu
            do i = 1, nocc
               tmp_t1(i,a) = tmp_t1(i,a) + 2*sum(I_vo*transpose(t2(:,i,:,a))) - sum(transpose(I_vo)*t2(i,:,:,a))
                           + 2*sum(t1*v_oovv(:,i,:,a)) - sum(transpose(t1), v_vovo(:,i,a,:))
            end do
         end do
         !$omp end parallel

         ! - v_ei^mn (2 t_mn^ea - t_mn^ae) 
         ! This is a clear case where we can benefit from dgemm 
         ! (see my benchmarking script at https://github.com/brianz98/fortran-tensor-benchmarking).
         ! As the number of contracted indices becomes large, the (significant) overhead of dgemm becomes negligible.

         ! v_ei^mn would require a tensor of type v_vooo, but we only have v_oovo, which is permutationally equivalent
         ! v_ei^mn would be accessed by v_oovo(m,i,e,n), which is why we need order=2,1,4,3 to permute it into (i,m,n,e)
         ! We would like to not store the reshaped arrays but gfortran doesn't like big arrays on the stack (ifort has the 
         ! handy -heap-arrays flag that directs them to the heap, but apparently there's no such flag in gfortran..)
         ! Anyways, we'll have to resort to using temporary arrays to ensure portability..
         ! Note that the temporary arrays do not need to be allocated beforehand, which is a relief..
         
         ! We further need a scratch t1 matrix as dgemm 'C' is overwritten
         allocate(tmp_t1_s(nocc,nvirt))
         reshape_tmp1 = reshape(v_oovo,(/nocc,nocc,nocc,nvirt/),order=(/2,1,4,3/))
         reshape_tmp2 = reshape(t2,(/nocc,nocc,nvirt,nvirt/),order=(/1,2,4,3/))
         call dgemm_wrapper('N','N',nocc,nvirt,nocc**2*nvirt,reshape_tmp1,2*t2-reshape_tmp2,tmp_t1_s)
         tmp_t1 = tmp_t1 - tmp_t1_s

         ! + v_ef^ma (2 t_mi^ef - t_im^ef)
         ! Let's massage the second quantity (the stuff in the brackets) into (i,m,e,f) and the first into (m,e,f,a) and voila
         reshape_tmp1 = reshape(v_vovv,(/nocc,nvirt,nvirt,nvirt/),order=(/2,4,3,1/))
         reshape_tmp2 = reshape(t2,(/nocc,nocc,nvirt,nvirt/),order=(/2,1,3,4/))
         call dgemm_wrapper('N','N',nocc,nvirt,nocc*nvirt**2,2*reshape_tmp2-t2,reshape_tmp1,tmp_t1_s)
         tmp_t1 = tmp_t1 + tmp_t1_s

         tmp_t1 = tmp_t1/D_ia
         ! We no longer need the t1 scratch matrix, but now we need a t2 one
         deallocate(tmp_t1_s)

         allocate(tmp_t2_s(nocc,nocc,nvirt,nvirt))
         ! Now update T2
         tmp_t2 = v_oovv

         ! -t_ij^ae*I_e^b, a nice dgemm
         call dgemm_wrapper('N','N',nocc**2*nvirt,nvirt,nvirt,t2,I_vv,tmp_t2_s)
         tmp_t2 = tmp_t2 + tmp_t2_s

         ! t_im^ab*I_j^m, benchmarking shows naive OMP is the fastest
         !$omp parallel do default(none)&
         !$omp schedule(static,10) collapse(2)&
         !$omp shared(nvl, nvu, nocc, tmp_t2, I_oo, t2)
         do b = nvl, nvu
            do a = nvl, nvu
               do j = 1, nocc
                  do i = 1, nocc
                     tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - sum(t2(i,:,a,b)*I_oo(j,:))
                  end do
               end do
            end do
         end do
         !$omp end parallel

         ! 1/2 v_ef^ab c_ij^ef, nice dgemm again
         call dgemm_wrapper('N','N',nocc**2,nvirt**2,nvirt**2,c_oovv,v_vvvv,tmp_t2_s)
         tmp_t2 = tmp_t2 + tmp_t2_s/2

         ! 1/2 c_mn^ab I_ij^mn, nice dgemm
         call dgemm_wrapper('N','N',nocc**2,nvirt**2,nocc**2,I_oooo,c_oovv,tmp_t2_s)
         tmp_t2 = tmp_t2 + tmp_t2_s/2

         ! -t_mj^ae I_ie^mb - I_ie^ma t_mj^ab + (2t_mi^ea - t_im^ea) I_ej^mb, seems hopeless, use OMP
         !$omp parallel do default(none)&
         !$omp schedule(static,10) collapse(4)&
         !$omp shared(nvl, nvu, nocc, tmp_t2, I_oo, t2)&
         !$omp private(tmp)
         do b = nvl, nvu
            do a = nvl, nvu
               do j = 1, nocc
                  do i = 1, nocc
                     tmp = 0.0_dp
                     do e = nvl, nvu
                        do m = 1, nocc
                           tmp = tmp - t2(m,j,a,e)*I_ovoo(i,e,m,b) - I_ovoo(i,e,m,a)*t2(m,j,a,b) &
                                 + 2*(t2(m,i,e,a)-t2(i,m,e,a))*I_voov(e,j,m,b)
                        end do
                     end do
                     tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) + tmp
                  end do
               end do
            end do
         end do
         !$omp end parallel

         ! t_i^e I'_ej^ab - t_m^a I'_ij^mb, first dgemm, second OMP
         call dgemm('N','N',nocc,nocc*nvirt**2,nvirt,t1,I_vovv_p,tmp_t2_s)
         tmp_t2 = tmp_t2 + tmp_t2_s

         !$omp parallel do default(none)&
         !$omp schedule(static,10) collapse(2)&
         !$omp shared(nvl, nvu, nocc, tmp_t2, I_ooov_p, t1)
         do b = nvl, nvu
            do a = nvl, nvu
               do j = 1, nocc
                  do i = 1, nocc
                     tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - sum(t1(:,a)*I_ooov_p(i,j,:,b))
                  end do
               end do
            end do
         end do
         !$omp end parallel


         end associate

      end subroutine update_amplitudes_restricted

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
         !$omp do schedule(runtime) collapse(2)
         do a = 2*nocc+1, 2*nbasis
            do i = 1, 2*nocc
               do m = 1, 2*nocc
                  tmp_t1(i,a) = tmp_t1(i,a) - t1(m,a)*F_oo(m,i)
                  do e = 2*nocc+1, 2*nbasis
                     tmp_t1(i,a) = tmp_t1(i,a) + t2(m,i,e,a)*F_ov(m,e)
                     do f = 2*nocc+1, 2*nbasis
                        tmp_t1(i,a) = tmp_t1(i,a) + 0.5*t2(m,i,f,e)*asym(f,e,m,a)
                     end do
                     do n = 1, 2*nocc
                        tmp_t1(i,a) = tmp_t1(i,a) - 0.5*t2(n,m,e,a)*asym(n,m,e,i)
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

         !$omp do schedule(runtime) collapse(4)
         do b = 2*nocc+1, 2*nbasis
            do a = 2*nocc+1, 2*nbasis
               do j = 1, 2*nocc
                  do i = 1, 2*nocc
                     tmp_t2(i,j,a,b) = asym(i,j,a,b)
                     do e = 2*nocc+1, 2*nbasis
                        tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) + t2(j,i,e,a)*F_vv(b,e) - t2(j,i,e,b)*F_vv(a,e) &
                        + tmp_t1(i,e)*asym(e,j,a,b) - tmp_t1(j,e)*asym(e,i,a,b)
                        do m = 1, 2*nocc
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - 0.5*F_ov(m,e)*(t2(i,j,a,e)*tmp_t1(m,b)-&
                              t2(i,j,b,e)*tmp_t1(m,a))
                        end do
                        do f = 2*nocc+1, 2*nbasis
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) + 0.5*tau(j,i,f,e)*W_vvvv(a,b,e,f)
                        end do
                     end do
                     do m = 1, 2*nocc
                        tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - t2(m,i,b,a)*F_oo(m,j) + t2(m,j,b,a)*F_oo(m,i) &
                        - tmp_t1(m,a)*asym(m,b,i,j) + tmp_t1(m,b)*asym(m,a,i,j)
                        do e = 2*nocc+1, 2*nbasis
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - 0.5*F_ov(m,e)*(t2(i,m,a,b)*tmp_t1(j,e)-&
                              t2(j,m,a,b)*tmp_t1(i,e))
                        end do
                        do n = 1, 2*nocc
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) + 0.5*tau(n,m,b,a)*W_oooo(m,n,i,j)
                        end do
                     end do

                     do m = 1, 2*nocc
                        do e = 2*nocc+1, 2*nbasis
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) &
                           + t2(m,i,e,a)*W_ovvo(m,b,e,j) - tmp_t1(i,e)*tmp_t1(m,a)*asym(e,j,m,b)&
                           - t2(m,j,e,a)*W_ovvo(m,b,e,i) + tmp_t1(j,e)*tmp_t1(m,a)*asym(e,i,m,b)&
                           - t2(m,i,e,b)*W_ovvo(m,a,e,j) + tmp_t1(i,e)*tmp_t1(m,b)*asym(e,j,m,a)&
                           + t2(m,j,e,b)*W_ovvo(m,a,e,i) - tmp_t1(j,e)*tmp_t1(m,b)*asym(e,i,m,a)
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
            !$omp do schedule(runtime) collapse(4) reduction(+:ecc,rmst2)
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

      subroutine do_ccsd_t_spinorb(sys, int_store)
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
         !$omp do schedule(runtime) collapse(3) reduction(+:e_T)
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
                                             + t2(k,j,f,a)*asym(f,i,b,c) - t2(k,i,f,a)*asym(f,j,b,c) - t2(i,j,f,a)*asym(f,k,b,c) &
                                             - t2(k,j,f,b)*asym(f,i,a,c) + t2(k,i,f,b)*asym(f,j,a,c) + t2(i,j,f,b)*asym(f,k,a,c) &
                                             - t2(k,j,f,c)*asym(f,i,b,a) + t2(k,i,f,c)*asym(f,j,b,a) + t2(i,j,f,c)*asym(f,k,b,a)
                           end do
                           do m = 1, 2*nocc
                              tmp_t3c(a,b,c) = tmp_t3c(a,b,c) &
                                             - t2(m,i,c,b)*asym(m,a,j,k) + t2(m,j,c,b)*asym(m,a,i,k) + t2(m,k,c,b)*asym(m,a,j,i) &
                                             + t2(m,i,c,a)*asym(m,b,j,k) - t2(m,j,c,a)*asym(m,b,i,k) - t2(m,k,c,a)*asym(m,b,j,i) &
                                             + t2(m,i,a,b)*asym(m,c,j,k) - t2(m,j,a,b)*asym(m,c,i,k) - t2(m,k,a,b)*asym(m,c,j,i)
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
      end subroutine do_ccsd_t_spinorb

      subroutine do_ccsd_t_spatial(sys, int_store)
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
         !$omp do schedule(runtime) collapse(3) reduction(+:e_T)
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
                                             + t2(k,j,f,a)*asym(f,i,b,c) - t2(k,i,f,a)*asym(f,j,b,c) - t2(i,j,f,a)*asym(f,k,b,c) &
                                             - t2(k,j,f,b)*asym(f,i,a,c) + t2(k,i,f,b)*asym(f,j,a,c) + t2(i,j,f,b)*asym(f,k,a,c) &
                                             - t2(k,j,f,c)*asym(f,i,b,a) + t2(k,i,f,c)*asym(f,j,b,a) + t2(i,j,f,c)*asym(f,k,b,a)
                           end do
                           do m = 1, 2*nocc
                              tmp_t3c(a,b,c) = tmp_t3c(a,b,c) &
                                             - t2(m,i,c,b)*asym(m,a,j,k) + t2(m,j,c,b)*asym(m,a,i,k) + t2(m,k,c,b)*asym(m,a,j,i) &
                                             + t2(m,i,c,a)*asym(m,b,j,k) - t2(m,j,c,a)*asym(m,b,i,k) - t2(m,k,c,a)*asym(m,b,j,i) &
                                             + t2(m,i,a,b)*asym(m,c,j,k) - t2(m,j,a,b)*asym(m,c,i,k) - t2(m,k,a,b)*asym(m,c,j,i)
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
      end subroutine do_ccsd_t_spatial

end module ccsd
