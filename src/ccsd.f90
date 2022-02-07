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

      ! Various slices of asymmetrised integrals
      real(dp), dimension(:,:,:,:), allocatable :: oooo, ovoo, oovo, ooov, oovv, ovov, ovvv, vovv, vvvv, ovvo
   end type cc_int_t

   type diis_cc_t
      ! Contains information for use in the CCSD DIIS procedure 
      ! Accelerating the convergence of the coupled-cluster approach: The use of the DIIS method
      ! Gustavo E. Scuseria, Timothy J. Lee, Henry F. Schaefer III
      ! (https://doi.org/10.1016/0009-2614(86)80461-4)

      ! These are identical as the HF-SCF ones, other than that we store and use the T matrices, 
      ! instead of Fock matrices for error matrix computation. 
      ! The details of the update is slightly subtle and quite easy to get wrong:
      ! Diagrammatically, 
      ! MP1 amp   (T0')-------> T1 ====> T1' 
      !                amp. eqs.   DIIS
      !             T1'-------> T2 ====> T2'
      !             T2'-------> T3 ====> T3'
      ! and so on. 
      ! The error matrices are defined: e_i = T_i - T_{i-1}'
      ! and the linear combination is over the unprimed amplitudes.
      ! This means we **have to** store both the unprimed amplitudes and the error matrices, plus two temporary arrays
      ! to compute the error matrix at each iteration

      logical :: use_diis = .true.
      integer :: n_errmat = 8
      integer :: n_active = 0
      integer :: iter = 0
      real(p), dimension(:,:,:), allocatable :: t1, e1, t1_s(:,:)
      real(p), dimension(:,:,:,:,:), allocatable :: t2, e2, t2_s(:,:,:,:)
      real(p), allocatable :: B(:,:)
      real(p), allocatable :: c(:)
      real(p), allocatable :: rhs(:)
   end type diis_cc_t
   
   contains

      subroutine do_ccsd_spinorb(sys, int_store)
         ! This is the spinorbital formulation of J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett, 
         ! J. Chem. Phys. volume 94, pp. 4334-4345 (1991) (https://doi.org/10.1063/1.460620)

         use integrals, only: int_store_t, eri_ind
         use error_handling, only: error, check_allocate
         use system, only: state_t, CCSD_T_spinorb

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(inout) :: int_store
         integer :: p, q, r, s
         integer :: pr, qr, pa, pb, qa, qb, ra, rb, sa, sb
         real(dp) :: prqs, psqr
         real(dp) :: err
         integer(kind=8) :: c_max, c_rate, t0, t1

         type(state_t) :: st
         type(diis_cc_t) :: diis
         
         type(cc_amp_t) :: cc_amp
         integer :: ierr, iter
         integer, parameter :: iunit = 6
         integer, parameter :: maxiter = 20

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
         associate(nbasis=>sys%nbasis,nocc=>sys%nocc,nvirt=>sys%nvirt,asym=>int_store%asym_spinorb)
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
         
         if (err > depsilon) then
            write(iunit, '(1X, A, 1X, E15.6)') 'Permutational symmetry error:', err
            call error('ccsd::do_ccsd', 'Permutational symmetry of antisymmetrised integrals does not hold')
         end if

         call system_clock(t1)
         if (t1<t0) t1 = t1+c_max
         write(iunit, '(1X, A, 1X, F8.6, A)') 'Time taken:', real(t1-t0, kind=dp)/c_rate, " s"
         write(iunit, *)
         t0=t1

         write(iunit, '(1X, A)') 'Forming slices of antisymmetrised spinorbital ERIs'

         ! We store only the slices needed: note that since we're moving to restricted indices,
         ! the virtuals are 1:nvirt instead of nocc+1:nbasis. Also in system.f90 we take care of the factor of 2 
         ! between restricted/unrestricted bounds

         ! All 4 occupied
         allocate(cc_int%oooo, source=asym(1:nocc,1:nocc,1:nocc,1:nocc))
         ! 3 occupied, 1 virtual, of which there are 3 types: ovoo, oovo, ooov, which are all related by permutational sym
         allocate(cc_int%ooov, source=asym(1:nocc,1:nocc,1:nocc,nocc+1:2*nbasis))
         allocate(cc_int%ovoo, source=asym(1:nocc,nocc+1:2*nbasis,1:nocc,1:nocc))
         allocate(cc_int%oovo, source=asym(1:nocc,1:nocc,nocc+1:2*nbasis,1:nocc))
         ! 2o2v
         allocate(cc_int%oovv, source=asym(1:nocc,1:nocc,nocc+1:2*nbasis,nocc+1:2*nbasis))
         allocate(cc_int%ovvo, source=asym(1:nocc,nocc+1:2*nbasis,nocc+1:2*nbasis,1:nocc))
         ! 1o3v
         allocate(cc_int%ovvv, source=asym(1:nocc,nocc+1:2*nbasis,nocc+1:2*nbasis,nocc+1:2*nbasis))
         allocate(cc_int%vovv, source=asym(nocc+1:2*nbasis,1:nocc,nocc+1:2*nbasis,nocc+1:2*nbasis))
         ! 4v
         allocate(cc_int%vvvv, source=asym(nocc+1:2*nbasis,nocc+1:2*nbasis,nocc+1:2*nbasis,nocc+1:2*nbasis))

         call system_clock(t1)
         if (t1<t0) t1 = t1+c_max
         write(iunit, '(1X, A, 1X, F8.6, A)') 'Time taken:', real(t1-t0, kind=dp)/c_rate, " s"
         write(iunit, *)
         t0=t1

         deallocate(int_store%asym_spinorb)
         end associate

         ! ############################################
         ! Initialise intermediate arrays and DIIS data
         ! ############################################
         write(iunit, '(1X, A)') 'Initialise CC intermediate tensors and DIIS auxilliary arrays...'
         call init_cc(sys, cc_amp, cc_int)
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
         associate(nbasis=>sys%nbasis, nvirt=>sys%nvirt, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab)
         do iter = 1, sys%ccsd_maxiter
            if (diis%use_diis) then
               ! These are for computing the error matrices for DIIS procedure
               diis%t1_s = cc_amp%t_ia
               diis%t2_s = cc_amp%t_ijab
            end if

            ! Update intermediate tensors
            call build_tau(sys, cc_amp, cc_int)
            call build_F(sys, cc_int, cc_amp)
            call build_W(sys, cc_int, cc_amp)

            ! CC amplitude equations
            call update_amplitudes(sys, cc_amp, int_store, cc_int)
            call update_cc_energy(sys, st, cc_int, cc_amp, conv)
            call system_clock(t1)
            write(iunit, '(1X, A, 1X, I2, 2X, F15.12, 2X, F8.6, 1X, A)') 'Iteration',iter,st%energy,real(t1-t0, kind=dp)/c_rate,'s'
            t0=t1
            if (conv) then
               write(iunit, '(1X, A)') 'Convergence reached within tolerance.'
               write(iunit, '(1X, A, 1X, F15.12)') 'Final CCSD Energy (Hartree):', st%energy

               ! Copied for return 
               sys%e_ccsd = st%energy
               call move_alloc(cc_amp%t_ia, sys%t1)
               call move_alloc(cc_amp%t_ijab, sys%t2)
               if (sys%calc_type == CCSD_T_spinorb) then
                  call move_alloc(cc_int%vovv, int_store%vovv)
                  call move_alloc(cc_int%ovoo, int_store%ovoo)
                  allocate(int_store%vvoo(nvirt,nvirt,nocc,nocc))
                  int_store%vvoo = reshape(cc_int%oovv,shape(int_store%vvoo),order=(/3,4,1,2/))
                  deallocate(cc_int%oooo, cc_int%ooov, cc_int%oovo, cc_int%oovv, cc_int%ovvo, cc_int%ovvv, &
                     cc_int%vvvv)
               else
                  deallocate(cc_int%oooo, cc_int%ooov, cc_int%ovoo, cc_int%oovo, cc_int%oovv, cc_int%ovvo, cc_int%ovvv, &
                     cc_int%vovv, cc_int%vvvv)
               end if
               exit
            end if
            call update_diis_cc(diis, cc_amp)
         end do
         end associate

         ! Nonlazy deallocations
         deallocate(st%t2_old)
         call deallocate_cc_int_t(cc_int)
         if (diis%use_diis) call deallocate_diis_cc_t(diis)

      end subroutine do_ccsd_spinorb

      subroutine do_ccsd_spatial(sys, int_store)
         ! This is the spinorbital formulation of J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett, J. Chem. Phys. volume 94, pp. 4334-4345 (1991) (https://doi.org/10.1063/1.460620)
         ! [TODO]: closed-shell (spatial) formulation of Hirata, S.; Podeszwa, R.; Tobita, M.; Bartlett, R. J. J. Chem. Phys. 2004, 120 (6), 2581 (https://doi.org/10.1063/1.1637577)

         use integrals, only: int_store_t, eri_ind
         use error_handling, only: error, check_allocate
         use system, only: state_t

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(inout) :: int_store
         integer :: p, q, r, s
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
         call init_cc(sys, cc_amp, cc_int)
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
            call build_F(sys, cc_int, cc_amp)
            call build_W(sys, cc_int, cc_amp)

            ! CC amplitude equations
            call update_amplitudes(sys, cc_amp, int_store, cc_int)
            call update_cc_energy(sys, st, cc_int, cc_amp, conv)
            
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
         call deallocate_cc_int_t(cc_int)
         call deallocate_diis_cc_t(diis)


      end subroutine do_ccsd_spatial

      subroutine init_cc(sys, cc_amp, cc_int)
         ! Initialise many CC related quantities.
         ! In/out:
         !     sys: system under study.
         !     cc_amp: CC amplitudes
         !     cc_int: CC intermediates being initialised

         use integrals, only: int_store_t
         use error_handling, only: check_allocate 

         type(system_t), intent(inout) :: sys
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         integer :: i, j, a, b, ierr
         integer, parameter :: iunit = 6

         write(iunit, '(1X, A)') 'Forming energy denominator matrices...'
         ! Forming the energy denominator matrices
         ! Note that abcd...ijkl... etc are **restricted** indices, so the lower bounds can be set to one, i.e. default Fortran
         associate(nvirt=>sys%nvirt, nocc=>sys%nocc, e=>sys%canon_levels)
         allocate(cc_int%D_ia(nocc, nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%D_ia', nocc*nvirt, ierr)
         allocate(cc_int%D_ijab(nocc, nocc, nvirt, nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%D_ijab', nocc**2*nvirt**2, ierr)
         do i = 1, nocc/2
            do a = 1, nvirt/2
               cc_int%D_ia(2*i-1:2*i,2*a-1:2*a) = e(i)-e(a+nocc/2)
               do j = 1, nocc/2
                  do b = 1, nvirt/2
                     cc_int%D_ijab(2*i-1:2*i, 2*j-1:2*j, 2*a-1:2*a, 2*b-1:2*b) = e(i)+e(j)-e(a+nocc/2)-e(b+nocc/2)
                  end do
               end do
            end do
         end do

         ! This is only later useful in CCSD(T) but we can just do it here...
         allocate(sys%canon_levels_spinorb(2*sys%nbasis), source=0.0_dp)
         do i = 1, sys%nbasis
            sys%canon_levels_spinorb(2*i-1:2*i) = e(i)
         end do
      
         write(iunit, '(1X, A)') 'Allocating amplitude tensors...'
         ! Initialise t_i^a=0 and t_{ij}^{ab}=MP1 WF
         allocate(cc_amp%t_ijab(nocc,nocc,nvirt,nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_amp%t_ijab', (nocc*nvirt)**2, ierr)
         allocate(cc_int%tmp_tijab, mold=cc_amp%t_ijab, stat=ierr)
         call check_allocate('cc_int%tmp_tijab', (nocc*nvirt)**2, ierr)
         allocate(cc_amp%t_ia(nocc, nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_amp%t_ia', nocc*nvirt, ierr)
         allocate(cc_int%tmp_tia, mold=cc_amp%t_ia, stat=ierr)
         call check_allocate('cc_int%tmp_tia', nocc*nvirt, ierr)

         write(iunit, '(1X, A)') 'Forming initial amplitude guesses...'
         cc_amp%t_ijab = cc_int%oovv / cc_int%D_ijab

         write(iunit, '(1X, A)') 'Allocating stored intermediate tensors...'
         ! Allocate the intermediate tensors
         allocate(cc_int%F_vv(nvirt,nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%F_vv', nvirt**2, ierr)
         allocate(cc_int%F_oo(nocc,nocc), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%F_oo', nocc**2, ierr)
         allocate(cc_int%F_ov(nocc,nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%F_ov', nocc*nvirt, ierr)
         allocate(cc_int%W_oooo(nocc,nocc,nocc,nocc), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%W_oooo', nocc**4, ierr)
         allocate(cc_int%W_vvvv(nvirt,nvirt,nvirt,nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%W_vvvv', nvirt**4, ierr)
         allocate(cc_int%W_ovvo(nocc,nvirt,nvirt,nocc), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%W_ovvo', (nocc*nvirt)**2, ierr)
         allocate(cc_int%tau(nocc,nocc,nvirt,nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%tau', (nocc*nvirt)**2, ierr)
         allocate(cc_int%tau_tilde(nocc,nocc,nvirt,nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%tau_tilde', (nocc*nvirt)**2, ierr)
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

         diis%n_errmat = sys%ccsd_diis_n_errmat

         if (diis%n_errmat < 2) then
            ! Switch off diis
            diis%use_diis = .false.
         else
            diis%iter = 0; diis%n_active = 0

            associate(nvirt=>sys%nvirt, nocc=>sys%nocc)
               allocate(diis%t1(nocc,nvirt,diis%n_errmat), source=0.0_p, stat=ierr)
               call check_allocate('diis%t1', diis%n_errmat*nocc*nvirt, ierr)
               allocate(diis%e1, source=diis%t1, stat=ierr)
               call check_allocate('diis%e1', diis%n_errmat*nocc*nvirt, ierr)
               allocate(diis%t1_s(nocc,nvirt), source=0.0_p, stat=ierr)
               call check_allocate('diis%t1_s', nocc*nvirt, ierr)

               allocate(diis%t2(nocc,nocc,nvirt,nvirt,diis%n_errmat), source=0.0_p, stat=ierr)
               call check_allocate('diis%t2', diis%n_errmat*(nocc*nvirt)**2, ierr)
               allocate(diis%e2(nocc,nocc,nvirt,nvirt,diis%n_errmat), source=diis%t2, stat=ierr)
               call check_allocate('diis%e2', diis%n_errmat*(nocc*nvirt)**2, ierr)
               allocate(diis%t2_s(nocc,nocc,nvirt,nvirt), source=0.0_p, stat=ierr)
               call check_allocate('diis%t2_s', diis%n_errmat*(nocc*nvirt)**2, ierr)
            end associate
         end if
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

         associate(n=>diis%n_active, nerr=>diis%n_errmat, it=>diis%iter)
         if (diis%use_diis) then
            diis%iter = diis%iter+1
            ! A modulo behaviour, but makes sure we never get iter = 0
            if (diis%iter > diis%n_errmat) diis%iter = diis%iter - diis%n_errmat

            ! In the first few iterations we need to make sure we don't accidentally access the empty arrays
            if (diis%n_active < diis%n_errmat) diis%n_active = diis%n_active+1

            ! Store the current amplitudes
            diis%t1(:,:,diis%iter) = cc_amp%t_ia
            diis%t2(:,:,:,:,diis%iter) = cc_amp%t_ijab

            ! Calculate error matrices for the current iteration
            diis%e1(:,:,diis%iter) = cc_amp%t_ia - diis%t1_s
            diis%e2(:,:,:,:,diis%iter) = cc_amp%t_ijab - diis%t2_s

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
                  ! We compute the error matrices on the fly, linsolve/dsysv only need the lower triangle
                  ! Because of how the idx array works, c1 corresponds to the weightage of the newest error vector, and so on.
                  diis%B(i,j) = sum(diis%e1(:,:,i)*diis%e1(:,:,j)) + sum(diis%e2(:,:,:,:,i)*diis%e2(:,:,:,:,j))
               end do
            end do

            call linsolve(diis%B, diis%c, ierr)
            if (ierr /= 0) call error('ccsd::update_diis_cc', 'Linear solve failed!')

            cc_amp%t_ia = 0.0_p
            cc_amp%t_ijab = 0.0_p
            do i = 1, n
               cc_amp%t_ia = cc_amp%t_ia + diis%c(i) * diis%t1(:,:,i)
               cc_amp%t_ijab = cc_amp%t_ijab + diis%c(i) * diis%t2(:,:,:,:,i)
            end do
         end if
         end associate
      end subroutine update_diis_cc

      subroutine deallocate_diis_cc_t(diis)
         ! Nonlazy deallocation of the diis_cc_t object as they take up significant memory
         ! In/out:
         !     diis: the diis_cc_t object being deallocated

         type(diis_cc_t), intent(inout) :: diis

         deallocate(diis%t1, diis%t2, diis%B, diis%c, diis%rhs, diis%e1, diis%e2, diis%t1_s, diis%t2_s)
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
         real(p) :: x

         associate(nvirt=>sys%nvirt, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab,&
            tau=>cc_int%tau, tau_tilde=>cc_int%tau_tilde)
            tau = 0.0_p
            tau_tilde = 0.0_p
            !$omp parallel do default(none) &
            !$omp private(i,j,a,b,x) &
            !$omp shared(sys, cc_amp, cc_int) &
            !$omp schedule(static, 10) collapse(2)
            do b = 1, nvirt
               do a = 1, nvirt
                  do j = 1, nocc
                     do i = 1, nocc
                        x = t_ia(i,a)*t_ia(j,b) - t_ia(i,b)*t_ia(j,a)
                        tau_tilde(i,j,a,b) = t_ijab(i,j,a,b) + 0.5*x
                        tau(i,j,a,b) = tau_tilde(i,j,a,b) + 0.5*x
                     end do
                  end do
               end do
            end do
            !$omp end parallel do
         end associate

      end subroutine build_tau

      subroutine build_F(sys, cc_int, cc_amp)
         ! Build the two-index F intermediates,
         ! Eqs. 3-5, note as we use HF reference all terms involving the fock matrix vanishes.
         ! In:
         !     sys: System under study.
         !     cc_amp: CC amplitudes
         ! In/out:
         !     cc_int: CC intermediates with the F matrices updated.

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int
         
         integer :: a, e, f, m, n, i

         associate(nvirt=>sys%nvirt, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab,&
            F_vv=>cc_int%F_vv, F_oo=>cc_int%F_oo, F_ov=>cc_int%F_ov, tau_tilde=>cc_int%tau_tilde)
         
         F_vv = 0.0_p
         F_oo = 0.0_p
         F_ov = 0.0_p

         !$omp parallel default(none) &
         !$omp private(a, e, f, m, n, i) &
         !$omp shared(sys, cc_int, cc_amp)

         ! F_ae = t_m^f * <ma||fe> - 0.5 tau~_{mn}^{af} <mn||ef>
         !$omp do schedule(static, 10) collapse(2)
         do a = 1, nvirt
            do e = 1, nvirt
               F_vv(a,e) = F_vv(a,e) + sum(t_ia(:,:)*cc_int%ovvv(:,a,:,e)) + 0.5*sum(tau_tilde(:,:,:,a)*cc_int%oovv(:,:,e,:))
            end do
         end do
         !$omp end do

         ! F_mi = t_n^e <mn||ie> + 0.5 tau~_{in}^{ef} <mn||ef>
         !$omp do schedule(static, 10) collapse(2)
         do m = 1, nocc
            do i = 1, nocc
               F_oo(m,i) = F_oo(m,i) + 0.5*sum(tau_tilde(i,:,:,:)*cc_int%oovv(m,:,:,:)) - sum(t_ia(:,:)*cc_int%ooov(:,m,i,:))
            end do
         end do
         !$omp end do

         ! F_me = \sum_{nf} t_n^f * <mn||ef>
         !$omp do schedule(static, 10) collapse(2)
         do e = 1, nvirt
            do m = 1, nocc
               F_ov(m,e) = F_ov(m,e) + sum(t_ia(:,:)*cc_int%oovv(m,:,e,:))
            end do
         end do
         !$omp end do
         !$omp end parallel
         end associate
      end subroutine build_F

      subroutine build_W(sys, cc_int, cc_amp)
         ! Build the four-index W intermediates, Eqs. 6 and 8.
         ! Note that we use a factor of 0.5 in the innermost line of the W_oooo loop instead of 0.25 per main text of Stanton,
         ! this is to make allowance for when W_vvvv cannot be stored due to memory limitation and has to be computed on the fly.
         ! See appendix of the same paper.
         ! In:
         !     sys: System under study.
         !     cc_amp: CC amplitudes
         ! In/out:
         !     cc_int: CC intermediates with the F matrices updated.

         use linalg, only: dgemm_wrapper

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int
         integer :: m, n, j, e, f, b
         real(p), dimension(:,:,:,:), allocatable :: scratch, reshape_scratch

         associate(nbasis=>sys%nbasis, nocc=>sys%nocc, nvirt=>sys%nvirt, t1=>cc_amp%t_ia, t2=>cc_amp%t_ijab,&
            W_oooo=>cc_int%W_oooo, W_ovvo=>cc_int%W_ovvo, tau=>cc_int%tau, W_vvvv=>cc_int%W_vvvv)

         ! Instead of laying it out like (m,n,i,j) we use (i,j,m,n) because then the contraction tau_mn^ab W_mnij can be processed
         ! by a simple dgemm call
         allocate(scratch,reshape_scratch, mold=W_oooo)
         ! #########################################################################################################################
         ! Eq. (6): W_mnij = <mn||ij> + P_(ij) t_j^e <mn||ie> + 1/2 tau_ij^ef <mn||ef>
         ! #########################################################################################################################
         ! + P_(ij) t_j^e <mn||ie>
         call dgemm_wrapper('N','T',nocc**3,nocc,nvirt,cc_int%ooov,t1,scratch)
         reshape_scratch = reshape(scratch, shape(reshape_scratch), order=(/1,2,4,3/))
         W_oooo = cc_int%oooo + scratch - reshape_scratch
         deallocate(reshape_scratch)

         ! 1/2 tau_ijef <mn||ef>
         allocate(reshape_scratch(nvirt,nvirt,nocc,nocc))
         reshape_scratch = reshape(tau, shape(reshape_scratch), order=(/3,4,1,2/))
         call dgemm_wrapper('N','N',nocc**2,nocc**2,nvirt**2,cc_int%oovv,reshape_scratch,W_oooo,0.5_dp,1.0_dp)
         deallocate(reshape_scratch)

         ! We reshape it to W_ijmn to make contractions more amenable to dgemm
         allocate(reshape_scratch, mold=W_oooo)
         reshape_scratch = reshape(W_oooo, shape(reshape_scratch), order=(/3,4,1,2/))
         W_oooo = reshape_scratch
         deallocate(scratch,reshape_scratch)

         ! #########################################################################################################################
         ! Eq. (7): W_abef = <ab||ef> - P_(ab)(t_mb<am||ef>)
         ! #########################################################################################################################
         ! - P_(ab) t_m^b <am||ef> = + P_(ab) (t_b^m)^T <ma||ef>
         allocate(scratch,reshape_scratch, mold=W_vvvv)
         call dgemm_wrapper('T','N',nvirt,nvirt**3,nocc,t1,cc_int%ovvv,scratch)
         reshape_scratch = reshape(scratch,shape(reshape_scratch),order=(/2,1,3,4/))
         ! The reshaped tensor is the unpermuted one!
         W_vvvv = cc_int%vvvv + reshape_scratch - scratch

         ! We reshape it to W_efab to make contractions more amenable to dgemm
         reshape_scratch = reshape(W_vvvv,shape(reshape_scratch),order=(/3,4,1,2/))
         W_vvvv = reshape_scratch
         deallocate(scratch,reshape_scratch)

         ! #########################################################################################################################
         ! Eq. (8): W_mbej = <mb||ej> + t_jf<mb||ef> - t_nb<mn||ej> - (1/2 t_jnfb + t_jf t_nb)<mn||ef>
         ! #########################################################################################################################
         ! <mb||ej> + t_jf<mb||ef>
         call dgemm_wrapper('N','T',nocc*nvirt**2,nocc,nvirt,cc_int%ovvv,t1,W_ovvo)
         W_ovvo = W_ovvo + cc_int%ovvo

         ! - t_nb<mn||ej> = +(t_bn)^T<nm||ej>
         allocate(scratch, reshape_scratch, mold=W_ovvo)
         call dgemm_wrapper('T','N',nvirt,nocc**2*nvirt,nocc,t1,cc_int%oovo,scratch)
         ! We need to reshape it because currently the scratch tensor is b,m,e,j
         reshape_scratch = reshape(scratch,shape(reshape_scratch),order=(/2,1,3,4/))
         W_ovvo = W_ovvo + reshape_scratch
         deallocate(scratch, reshape_scratch)

         ! - (1/2 t_jnfb + t_jf t_nb)<mn||ef>
         !$omp parallel default(none) &
         !$omp shared(cc_amp,cc_int)
         !$omp do schedule(static, 10) collapse(4)
         do j = 1, nocc
            do e = 1, nvirt
               do b = 1, nvirt
                  do m = 1, nocc
                     do n = 1, nocc
                        do f = 1, nvirt
                           ! [todo]: bad cache behaviour (last term), 
                           ! however dgemm isn't that much better here since two reshapes needed.
                           W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) - (0.5*t2(j,n,f,b) + t1(n,b)*t1(j,f)) * cc_int%oovv(m,n,e,f)
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
         use linalg, only: dgemm_wrapper

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(in) :: int_store
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         real(p), dimension(:,:,:,:), allocatable :: tmp_t2_s, reshape_tmp

         integer :: i, j, a, b, m, e

         associate(tmp_t1=>cc_int%tmp_tia,tmp_t2=>cc_int%tmp_tijab,t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab,asym=>int_store%asym_spinorb, &
            F_oo=>cc_int%F_oo,F_ov=>cc_int%F_ov,F_vv=>cc_int%F_vv,D_ia=>cc_int%D_ia,D_ijab=>cc_int%D_ijab,W_oooo=>cc_int%W_oooo, &
            W_ovvo=>cc_int%W_ovvo,tau=>cc_int%tau,tau_tilde=>cc_int%tau_tilde,nocc=>sys%nocc,nvirt=>sys%nvirt,&
            W_vvvv=>cc_int%W_vvvv)

         ! #########################################################################################################################
         ! Update T1
         ! Eq. (1) D_ia t_i^a = t_ie F_ae - t_ma F_mi + t_imae F_me - t_nf <na||if> -1/2 t_imef <ma||ef> - 1/2 t_mnae <nm||ei>
         ! #########################################################################################################################

         ! #### t_ie F_ae
         tmp_t1 = matmul(t1,transpose(F_vv))
         ! #### -t_ma F_mi
         tmp_t1 = tmp_t1 - matmul(transpose(F_oo),t1)

         !$omp parallel default(none) &
         !$omp shared(sys, cc_amp, int_store, cc_int)
         ! tmp_t1 can be shared as each loop iteration will update a unique element of it.
         !$omp do schedule(static, 10) collapse(2)
         do a = 1, nvirt
            do i = 1, nocc
               ! We can't do sum(A + B) because the order of (i,a) doesn't generally match between these
               tmp_t1(i,a) = tmp_t1(i,a) + sum(t1(:,:)*cc_int%ovvo(:,a,:,i)) + sum(t2(:,i,:,a)*F_ov(:,:)) &
                                         + 0.5*(sum(t2(:,i,:,:)*cc_int%ovvv(:,a,:,:)) - sum(t2(:,:,:,a)*cc_int%oovo(:,:,:,i)))
            end do
         end do
         !$omp end do
         !$omp end parallel

         tmp_t1 = tmp_t1/D_ia
         
         ! #########################################################################################################################
         ! Update T2
         ! Eq. (2) D_ijab t_ijab = <ij||ab> + P_ab(t_ijae (F_be - 1/2 t_mb F_me)) - P_(ij)(t_imab (F_mj + 1/2 t_je F_me)) 
         !                         + 1/2 tau_mnab W_mnij + 1/2 tau_ijef W_abef + P_(ij)P_(ab)(t_imae W_mbej - t_ie t_ma <mb||ej>) 
         !                         + P_(ij)(t_ie<ab||ej>) - P_(ab)(t_ma <mb||ij>)
         ! #########################################################################################################################
         ! <ij||ab>
         tmp_t2 = cc_int%oovv
         allocate(tmp_t2_s, reshape_tmp, mold=t2)
         !$omp parallel default(none) &
         !$omp shared(sys, cc_amp, int_store, cc_int, tmp_t2_s)
         !$omp do schedule(static, 10) collapse(4)
         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     tmp_t2_s(i,j,a,b) = sum(t2(:,i,:,a)*W_ovvo(:,b,:,j))
                     do m = 1, nocc
                        do e = 1, nvirt
                           tmp_t2_s(i,j,a,b) = tmp_t2_s(i,j,a,b) - t1(i,e)*t1(m,a)*cc_int%ovvo(m,b,e,j)
                        end do 
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel
         reshape_tmp = reshape(tmp_t2_s,shape(reshape_tmp),order=(/2,1,3,4/))
         tmp_t2 = tmp_t2 + tmp_t2_s - reshape_tmp
         reshape_tmp = reshape(tmp_t2_s,shape(reshape_tmp),order=(/1,2,4,3/))
         tmp_t2 = tmp_t2 - reshape_tmp
         reshape_tmp = reshape(tmp_t2_s,shape(reshape_tmp),order=(/2,1,4,3/))
         tmp_t2 = tmp_t2 + reshape_tmp

         ! P_(ab)(t_ijae F_be)
         call dgemm_wrapper('N','T',nocc**2*nvirt,nvirt,nvirt,t2,F_vv,tmp_t2_s)
         reshape_tmp = reshape(tmp_t2_s,shape(reshape_tmp),order=(/1,2,4,3/))
         tmp_t2 = tmp_t2 + tmp_t2_s - reshape_tmp
         ! -1/2 P_(ab)(t_ijae t_mb F_me)
         call dgemm_wrapper('N','T',nocc**2*nvirt,nvirt,nvirt,t2,matmul(transpose(t1),F_ov),tmp_t2_s)
         reshape_tmp = reshape(tmp_t2_s,shape(reshape_tmp),order=(/1,2,4,3/))
         tmp_t2 = tmp_t2 - (tmp_t2_s - reshape_tmp)/2
         ! -1/2 P_(ij)(t_je F_me t_imab)
         call dgemm_wrapper('N','N',nocc,nocc*nvirt**2,nocc,matmul(t1,transpose(F_ov)),t2,tmp_t2_s)
         reshape_tmp = reshape(tmp_t2_s,shape(reshape_tmp),order=(/2,1,3,4/))
         tmp_t2 = tmp_t2 - (tmp_t2_s - reshape_tmp)/2
         ! P_(ij)(t_ie<ab||ej>) = P_(ij)(t_ie<ej||ab>)
         call dgemm_wrapper('N','N',nocc,nocc*nvirt**2,nvirt,t1,cc_int%vovv,tmp_t2_s)
         reshape_tmp = reshape(tmp_t2_s,shape(reshape_tmp),order=(/2,1,3,4/))
         tmp_t2 = tmp_t2 + tmp_t2_s - reshape_tmp
         ! -P_(ab)(t_ma <mb||ij>) = +P_(ab)(<ij||bm> t_ma)
         call dgemm_wrapper('N','N',nocc**2*nvirt,nvirt,nocc,cc_int%oovo,t1,tmp_t2_s)
         reshape_tmp = reshape(tmp_t2_s,shape(reshape_tmp),order=(/1,2,4,3/))
         ! The reshaped tensor is the unpermuted one, since tmp_t2_s is i,j,b,a
         tmp_t2 = tmp_t2 + reshape_tmp - tmp_t2_s
         ! -P_(ij)(t_imab F_mj)
         call dgemm_wrapper('T','N',nocc,nocc*nvirt**2,nocc,F_oo,t2,tmp_t2_s)
         reshape_tmp = reshape(tmp_t2_s,shape(reshape_tmp),order=(/2,1,3,4/))
         tmp_t2 = tmp_t2 - tmp_t2_s + reshape_tmp
         ! + 1/2 tau_mnab W_mnij, but remember out W_mnij is reshaped as W_ijmn
         call dgemm_wrapper('N','N',nocc**2,nvirt**2,nocc**2,W_oooo,tau,tmp_t2, 0.5_dp, 1.0_dp)
         ! + 1/2 tau_ijef W_abef
         call dgemm_wrapper('N','N',nocc**2,nvirt**2,nvirt**2,tau,W_vvvv,tmp_t2, 0.5_dp, 1.0_dp)
         deallocate(tmp_t2_s, reshape_tmp)

         t1 = tmp_t1
         t2 = tmp_t2/D_ijab

         end associate

      end subroutine update_amplitudes

      subroutine update_cc_energy(sys, st, cc_int, cc_amp, conv)
         ! Updates the CC energy and also checks for convergence
         ! In:
         !     sys: system under study.
         !     cc_int: CC intermediates, of which we're using the cc_int%oovv slice.
         !     cc_amp: CC amplitudes.
         ! In/out:
         !     st: state_t object storing energies and RMS(T) info for convergence checking
         ! Out:
         !     conv: returns true if converged within tolerance

         use system, only: state_t
         use integrals, only: int_store_t

         type(system_t), intent(in) :: sys
         type(cc_int_t), intent(in) :: cc_int
         type(cc_amp_t), intent(in) :: cc_amp
         type(state_t), intent(inout) :: st
         logical, intent(out) :: conv

         integer :: i, j, a, b
         real(p) :: ecc, rmst2

         conv = .false.
         st%energy_old = st%energy
         associate(nocc=>sys%nocc, nvirt=>sys%nvirt, oovv=>cc_int%oovv, t1=>cc_amp%t_ia, t2=>cc_amp%t_ijab)
            ecc = 0.0_p
            rmst2 = 0.0_p
            !$omp parallel default(none) &
            !$omp private(i,j,a,b) &
            !$omp shared(sys,st,cc_amp,ecc,rmst2,cc_int) 
            !$omp do schedule(static, 10) collapse(2) reduction(+:ecc,rmst2)
            do b = 1, nvirt
               do a = 1, nvirt
                  do j = 1, nocc
                     do i = 1, nocc
                        ecc = ecc + 0.25*oovv(i,j,a,b)*(t2(i,j,a,b)+2*t1(i,a)*t1(j,b))
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
         real(p), dimension(:,:,:), allocatable :: tmp_t3d, tmp_t3c, tmp_t3c_d, reshape_tmp1, reshape_tmp2, t2_reshape(:,:,:,:)
         real(p) :: e_T
         integer, parameter :: iunit = 6

         write(iunit, '(1X, 10("-"))')
         write(iunit, '(1X, A)') 'CCSD(T)'
         write(iunit, '(1X, 10("-"))')

         associate(nvirt=>sys%nvirt, nocc=>sys%nocc, e=>sys%canon_levels_spinorb, t1=>sys%t1, t2=>sys%t2,&
                   vvoo=>int_store%vvoo,vovv=>int_store%vovv,ovoo=>int_store%ovoo)

         ! For efficient cache access during tmp_t3c updates
         allocate(t2_reshape(nvirt,nvirt,nocc,nocc), source=0.0_p, stat=ierr)
         call check_allocate('t2_reshape', (nvirt*nocc)**2, ierr)
         ! ijab -> baij
         t2_reshape = reshape(t2, shape(t2_reshape), order=(/4,3,2,1/))

         e_T = 0.0_p

         !$omp parallel default(none) &
         !$omp private(tmp_t3d, tmp_t3c, tmp_t3c_d, ierr, reshape_tmp1, reshape_tmp2) &
         !$omp shared(sys, int_store, t2_reshape) reduction(+:e_T)
         
         ! Declare heap-allocated arrays inside parallel region
         allocate(tmp_t3d(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('tmp_t3d', nvirt**3, ierr)
         allocate(tmp_t3c(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('tmp_t3c', nvirt**3, ierr)
         allocate(tmp_t3c_d(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('tmp_t3c_d', nvirt**3, ierr)
         allocate(reshape_tmp1(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('reshape_tmp1', nvirt**3, ierr)
         allocate(reshape_tmp2(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('reshape_tmp2', nvirt**3, ierr)

         ! We use avoid the storage of full triples (six-dimensional array..) and instead use the strategy of 
         ! batched triples storage, denoted W^{ijk}(abc) in (https://doi.org/10.1016/0009-2614(91)87003-T).
         ! We could of course only compute one element in a loop iteration but that will probably result in bad floating point
         ! performance, here for each thread/loop iteration we use the Fortran intrinsic sum, 
         ! instead of all relying on OMP reduction, to hopefully give better floating point performance.

         ! P(i/jk)P(a/bc) = (1-jik-kji)(1-bac-cba)
         !$omp do schedule(static, 10) collapse(3) 
         do i = 1, nocc
            do j = 1, nocc
               do k = 1, nocc
                  ! Compute T3 amplitudes
                  do a = 1, nvirt
                     do c = 1, nvirt
                        do b = 1, nvirt
                           ! Disonnected T3: D_{ijk}^{abc}*t_{ijk}^{abc}(d) = P(i/jk)P(a/bc)t_i^a*<jk||bc>
                           ! Can be rewritten directly as no sums in the expression
                           tmp_t3d(a,b,c) = (t1(i,a)*vvoo(b,c,j,k) - t1(j,a)*vvoo(b,c,i,k) - t1(k,a)*vvoo(b,c,j,i)) &
                                             /(e(i)+e(j)+e(k)-e(a+nocc)-e(b+nocc)-e(c+nocc))

                           ! Connected T3: D_{ijk}^{abc}*t_{ijk}^{abc}(c) = 
                           !        P(i/jk)P(a/bc)[\sum_f t_jk^af <fi||bc> - \sum_m t_im^bc <ma||jk>]
                           tmp_t3c(a,b,c) = sum(vovv(:,i,b,c)*t2_reshape(:,a,k,j)) - sum(vovv(:,j,b,c)*t2_reshape(:,a,k,i)) &
                                          - sum(vovv(:,k,b,c)*t2_reshape(:,a,i,j)) &
                              - sum(t2(:,i,c,b)*ovoo(:,a,j,k)) + sum(t2(:,j,c,b)*ovoo(:,a,i,k)) + sum(t2(:,k,c,b)*ovoo(:,a,j,i))
                           
                           tmp_t3c_d(a,b,c) = tmp_t3c(a,b,c)/(e(i)+e(j)+e(k)-e(a+nocc)-e(b+nocc)-e(c+nocc))
                        end do
                     end do
                  end do
                  ! Do reshapes
                  reshape_tmp1 = reshape(tmp_t3d, shape(reshape_tmp1), order=(/2,1,3/))
                  reshape_tmp2 = reshape(tmp_t3d, shape(reshape_tmp1), order=(/3,2,1/))
                  tmp_t3d = tmp_t3d - reshape_tmp1 - reshape_tmp2
                  
                  reshape_tmp1 = reshape(tmp_t3c, shape(reshape_tmp1), order=(/2,1,3/))
                  reshape_tmp2 = reshape(tmp_t3c, shape(reshape_tmp1), order=(/3,2,1/))
                  tmp_t3c = tmp_t3c - reshape_tmp1 - reshape_tmp2

                  reshape_tmp1 = reshape(tmp_t3c_d, shape(reshape_tmp1), order=(/2,1,3/))
                  reshape_tmp2 = reshape(tmp_t3c_d, shape(reshape_tmp1), order=(/3,2,1/))
                  tmp_t3c_d = tmp_t3c_d - reshape_tmp1 - reshape_tmp2

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
         !$omp do schedule(static, 10) collapse(3) reduction(+:e_T)
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
