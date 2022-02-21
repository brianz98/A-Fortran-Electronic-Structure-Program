module ccsd

   use const
   use, intrinsic :: iso_fortran_env, only : iunit=>output_unit
   use system, only: system_t
   use linalg, only: omp_reshape

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

      ! These are the quantities required by the spin-free formulation, left unallocated if using spinorbital formulation
      real(dp), dimension(:,:,:,:), allocatable :: v_oovv, v_ovov, v_vvov, v_oovo, v_oooo, v_vvvv
      real(dp), dimension(:,:), allocatable :: I_vo, I_vv, I_oo, I_oo_p
      real(dp), dimension(:,:,:,:), allocatable :: c_oovv, x_voov, asym_t2
      real(dp), dimension(:,:,:,:), allocatable :: I_oooo, I_ovov, I_voov, I_vovv_p, I_ooov_p
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

      subroutine do_ccsd_spinorb(sys, int_store, int_store_cc)
         ! This is the spinorbital formulation of J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett, 
         ! J. Chem. Phys. volume 94, pp. 4334-4345 (1991) (https://doi.org/10.1063/1.460620)

         use integrals, only: int_store_t, int_store_cc_t, eri_ind
         use error_handling, only: error, check_allocate
         use system, only: state_t, CC_SD_T

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(inout) :: int_store
         type(int_store_cc_t), intent(inout) :: int_store_cc
         integer :: p, q, r, s
         integer :: pr, qr, pa, pb, qa, qb, ra, rb, sa, sb
         real(dp) :: prqs, psqr
         real(dp) :: err
         integer(kind=8) :: c_max, c_rate, t0, t1

         type(state_t) :: st
         type(diis_cc_t) :: diis
         
         type(cc_amp_t) :: cc_amp
         integer :: ierr, iter
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
         call init_cc(sys, cc_amp, cc_int, int_store, restricted=.false.)
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
            call update_cc_energy(sys, st, cc_int, cc_amp, conv, restricted=.false.)
            call system_clock(t1)
            write(iunit, '(1X, A, 1X, I2, 2X, F15.12, 2X, F8.6, 1X, A)') 'Iteration',iter,st%energy,real(t1-t0, kind=dp)/c_rate,'s'
            t0=t1
            if (conv) then
               write(iunit, '(1X, A)') 'Convergence reached within tolerance.'
               write(iunit, '(1X, A, 1X, F15.12)') 'Final CCSD Energy (Hartree):', st%energy

               ! Copied for return 
               sys%e_ccsd = st%energy
               call move_alloc(cc_amp%t_ia, int_store_cc%t1)
               call move_alloc(cc_amp%t_ijab, int_store_cc%t2)
               if (sys%calc_type == CC_SD_T) then
                  call move_alloc(cc_int%vovv, int_store_cc%vovv)
                  call move_alloc(cc_int%ovoo, int_store_cc%ovoo)
                  allocate(int_store_cc%vvoo(nvirt,nvirt,nocc,nocc))
                  int_store_cc%vvoo = reshape(cc_int%oovv,shape(int_store_cc%vvoo),order=(/3,4,1,2/))
               end if
               exit
            end if
            call update_diis_cc(diis, cc_amp)
         end do
         end associate

         ! Nonlazy deallocations only for derived types that are intent(out) or intent(inout)
         deallocate(st%t2_old)

      end subroutine do_ccsd_spinorb

      subroutine do_ccsd_spatial(sys, int_store, int_store_cc)

         ! This is the spin-free (spatial) formulation of P. Piecuch et al., Computer Physics Communications 149 (2002) 71–96, 
         ! https://doi.org/10.1016/S0010-4655(02)00598-2

         use integrals, only: int_store_t, int_store_cc_t, eri_ind
         use error_handling, only: error, check_allocate
         use system, only: state_t, CC_SD_T

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(inout) :: int_store
         type(int_store_cc_t), intent(out) :: int_store_cc
         real(dp) :: t1_diagnostic
         integer(kind=8) :: c_max, c_rate, t0, t1

         type(state_t) :: st
         type(diis_cc_t) :: diis
         
         type(cc_amp_t) :: cc_amp
         integer :: iter
         integer, parameter :: maxiter = 20

         type(cc_int_t) :: cc_int
         logical :: conv, debug

         write(iunit, '(1X, 10("-"))')
         write(iunit, '(1X, A)') 'CCSD'
         write(iunit, '(1X, 10("-"))')
         call system_clock(count=t0, count_rate=c_rate, count_max=c_max)

         ! ############################################
         ! Initialise intermediate arrays and DIIS data
         ! ############################################
         write(iunit, '(1X, A)') 'Initialise CC intermediate tensors and DIIS auxilliary arrays...'
         call init_cc(sys, cc_amp, cc_int, int_store, restricted=.true.)
         allocate(st%t2_old, source=cc_amp%t_ijab)
         call init_diis_cc_t(sys, diis)

         call system_clock(t1)
         if (t1<t0) t1 = t1+c_max
         write(iunit, '(1X, A, 1X, F8.6, A)') 'Time taken:', real(t1-t0, kind=dp)/c_rate, " s"
         write(iunit, *)
         t0=t1

         write(iunit, '(1X, A)') 'Initialisation done, now entering iterative CC solver...'

         debug = .false.

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
            if (debug) then
               call update_restricted_intermediates_debug(sys, cc_amp, cc_int)
            else
               call update_restricted_intermediates(sys, cc_amp, cc_int)
            end if

            ! CC amplitude equations
            if (debug) then
               call update_amplitudes_restricted_debug(sys, cc_amp, cc_int)
            else
               call update_amplitudes_restricted(sys, cc_amp, cc_int)
            end if

            call update_cc_energy(sys, st, cc_int, cc_amp, conv, restricted=.true.)
            call system_clock(t1)
            write(iunit, '(1X, A, 1X, I2, 2X, F15.12, 2X, F8.6, 1X, A)') 'Iteration',iter,st%energy,real(t1-t0, kind=dp)/c_rate,'s'
            t0=t1
            if (conv) then
               write(iunit, '(1X, A)') 'Convergence reached within tolerance.'
               write(iunit, '(1X, A, 1X, F15.12)') 'Final CCSD Energy (Hartree):', st%energy
               ! The 'T1 diagnostic' as defined by https://doi.org/10.1007/978-94-011-0193-6_2
               ! For now since we're not doing actual unrestricted calculations this is only implemented for restricted,
               ! as the half of all t1 amp will be zero for 'unrestricted'.
               t1_diagnostic = sqrt(sum(cc_amp%t_ia**2))/sqrt(real(sys%nel,kind=dp))
               write(iunit, '(1X, A, 1X, F8.5)') 'T1 diagnostic:', t1_diagnostic
               if (t1_diagnostic > 0.02_dp) then
                  write(iunit, '(1X, A)') 'Significant multireference character detected, CCSD result might be unreliable!'
               end if
               if (sys%ccsd_t_comp_renorm) then
                  ! Too many of the quantities needed in forming the 'generalised moment' quantity M_ijk^abc(2) is in the
                  ! cc_int derived type, which we want to get automatically deallocated upon existing 'do_ccsd_spatial',
                  ! so we have to build M_ijk^abc(2) here
                  call build_cr_ccsd_t_intermediates(sys, cc_int, cc_amp, int_store_cc)
               end if
               ! Copied for return 
               sys%e_ccsd = st%energy
               call move_alloc(cc_amp%t_ia, int_store_cc%t1)
               call move_alloc(cc_amp%t_ijab, int_store_cc%t2)
               if (sys%calc_type == CC_SD_T) then
                  call move_alloc(cc_int%v_vvov, int_store_cc%v_vvov)
                  call move_alloc(cc_int%v_oovo, int_store_cc%v_oovo)
                  call move_alloc(cc_int%v_oovv, int_store_cc%v_oovv)
               end if
               exit
            end if
            call update_diis_cc(diis, cc_amp)
         end do
         end associate

         ! Nonlazy deallocations only for derived types that are intent(out) or intent(inout)
         deallocate(st%t2_old)

      end subroutine do_ccsd_spatial

      subroutine init_cc(sys, cc_amp, cc_int, int_store, restricted)

         ! Initialise many CC related quantities.
         ! In:
         !     int_store: int_store_t object with MO eri information in it
         !     restricted: true if we're using the spin-free formulation
         ! In/out:
         !     sys: system under study.
         !     cc_amp: CC amplitudes
         !     cc_int: CC intermediates being initialised

         use integrals, only: int_store_t, eri_ind
         use error_handling, only: check_allocate 

         type(int_store_t), intent(in) :: int_store
         logical, intent(in) :: restricted
         type(system_t), intent(inout) :: sys
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         real(dp), allocatable :: temperi(:,:,:,:)
         integer :: i, j, a, b, p, q, r, s, ierr

         write(iunit, '(1X, A)') 'Forming energy denominator matrices...'
         ! Forming the energy denominator matrices
         ! Note that abcd...ijkl... etc are **restricted** indices, so the lower bounds can be set to one, i.e. default Fortran
         associate(nbasis=>sys%nbasis, nvirt=>sys%nvirt, nocc=>sys%nocc, e=>sys%canon_levels)
         allocate(cc_int%D_ia(nocc, nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%D_ia', nocc*nvirt, ierr)
         allocate(cc_int%D_ijab(nocc, nocc, nvirt, nvirt), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%D_ijab', nocc**2*nvirt**2, ierr)
         if (restricted) then
            do i = 1, nocc
               do a = 1, nvirt
                  cc_int%D_ia(i,a) = e(i)-e(a+nocc)
                  do j = 1, nocc
                     do b = 1, nvirt
                        cc_int%D_ijab(i,j,a,b) = e(i)+e(j)-e(a+nocc)-e(b+nocc)
                     end do
                  end do
               end do
            end do
         else
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
         end if

         ! This is only later useful in CCSD(T) but we can just do it here...
         allocate(sys%canon_levels_spinorb(2*sys%nbasis), source=0.0_dp)
         do i = 1, sys%nbasis
            sys%canon_levels_spinorb(2*i-1:2*i) = e(i)
         end do
      
         ! Shared by both spinorbital and spatial formulations
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

         write(iunit, '(1X, A)') 'Forming ERI slices...'
         if (restricted) then
            ! Re-sort MO ERIs into different slices for efficient vectorised operations
            allocate(cc_int%v_oovv(nocc,nocc,nvirt,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_oovv', (nocc*nvirt)**2, ierr)
            allocate(cc_int%v_ovov(nocc,nvirt,nocc,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_ovov', (nocc*nvirt)**2, ierr)
            allocate(cc_int%v_vvov(nvirt,nvirt,nocc,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_vvov', nocc*nvirt**3, ierr)
            allocate(cc_int%v_oovo(nocc,nocc,nvirt,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_oovo', nocc**3*nvirt, ierr)
            allocate(cc_int%v_oooo(nocc,nocc,nocc,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_oooo', nocc**4, ierr)
            allocate(cc_int%v_vvvv(nvirt,nvirt,nvirt,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%v_vvvv', nvirt**4, ierr)

            allocate(temperi(nbasis,nbasis,nbasis,nbasis), source=0.0_dp, stat=ierr)
            call check_allocate('temperi', nbasis**4, ierr)
            
            do s = 1, nbasis
               do r = 1, nbasis
                  do q = 1, nbasis
                     do p = 1, nbasis
                        ! We're going from the chemists' to the physicists' notation
                        temperi(p,q,r,s) = int_store%eri_mo(eri_ind(eri_ind(p,r),eri_ind(q,s)))
                     end do
                  end do
               end do
            end do

            cc_int%v_oovv = temperi(1:nocc,1:nocc,nocc+1:nbasis,nocc+1:nbasis)
            cc_int%v_ovov = temperi(1:nocc,nocc+1:nbasis,1:nocc,nocc+1:nbasis)
            cc_int%v_vvov = temperi(nocc+1:nbasis,nocc+1:nbasis,1:nocc,nocc+1:nbasis)
            cc_int%v_oovo = temperi(1:nocc,1:nocc,nocc+1:nbasis,1:nocc)
            cc_int%v_oooo = temperi(1:nocc,1:nocc,1:nocc,1:nocc)
            cc_int%v_vvvv = temperi(nocc+1:nbasis,nocc+1:nbasis,nocc+1:nbasis,nocc+1:nbasis)

            deallocate(temperi)

         end if

         write(iunit, '(1X, A)') 'Forming initial amplitude guesses...'
         if (restricted) then
            cc_amp%t_ia = 0.0_dp
            cc_amp%t_ijab = cc_int%v_oovv / cc_int%D_ijab
         else
            cc_amp%t_ijab = cc_int%oovv / cc_int%D_ijab
         end if

         write(iunit, '(1X, A)') 'Allocating stored intermediate tensors...'
         if (restricted) then
            ! Two-index intermediates
            allocate(cc_int%I_vo(nvirt, nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_vo', nocc*nvirt, ierr)
            allocate(cc_int%I_vv(nvirt, nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_vv', nvirt**2, ierr)
            allocate(cc_int%I_oo_p(nocc, nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_oo_p', nocc**2, ierr)
            allocate(cc_int%I_oo(nocc, nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_oo', nocc**2, ierr)

            ! Four-index intermediates
            allocate(cc_int%c_oovv(nocc,nocc,nvirt,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%c_oovv', (nocc*nvirt)**2, ierr)
            allocate(cc_int%x_voov(nvirt,nocc,nocc,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%x_voov', (nocc*nvirt)**2, ierr)
            allocate(cc_int%I_oooo(nocc,nocc,nocc,nocc), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_oooo', nocc**4, ierr)
            allocate(cc_int%I_ovov(nocc,nvirt,nocc,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_ovov', (nocc*nvirt)**2, ierr)
            allocate(cc_int%I_voov(nvirt,nocc,nocc,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_voov', (nocc*nvirt)**2, ierr)
            allocate(cc_int%I_vovv_p(nvirt,nocc,nvirt,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_vovv_p', nocc*nvirt**3, ierr)
            allocate(cc_int%I_ooov_p(nocc,nocc,nocc,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%I_ooov_p', nvirt*nocc**3, ierr)
            allocate(cc_int%asym_t2(nocc,nocc,nvirt,nvirt), source=0.0_dp, stat=ierr)
            call check_allocate('cc_int%asym_t2', (nocc*nvirt)**2, ierr)
         else
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
         end if
         end associate

      end subroutine init_cc

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
         use linalg, only: dgemm_wrapper

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         real(p), dimension(:,:,:,:), allocatable :: tmp_4_1, tmp_4_2

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
               do f = 1, nvirt
                  do m = 1, nocc
                     F_vv(a,e) = F_vv(a,e) + t_ia(m,f)*cc_int%ovvv(m,a,f,e)
                  end do
               end do
            end do
         end do
         !$omp end do

         ! F_mi = t_n^e <mn||ie> + 0.5 tau~_{in}^{ef} <mn||ef>
         !$omp do schedule(static, 10) collapse(2)
         do m = 1, nocc
            do i = 1, nocc
               do e = 1, nvirt
                  do n = 1, nocc
                     F_oo(m,i) = F_oo(m,i) - t_ia(n,e)*cc_int%ooov(n,m,i,e)
                  end do
               end do
            end do
         end do
         !$omp end do

         ! F_me = \sum_{nf} t_n^f * <mn||ef>
         !$omp do schedule(static, 10) collapse(2)
         do e = 1, nvirt
            do m = 1, nocc
               do f = 1, nvirt
                  do n = 1, nocc
                     F_ov(m,e) = F_ov(m,e) + t_ia(n,f)*cc_int%oovv(m,n,e,f)
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel

         ! amnf, mnfe -> ae
         allocate(tmp_4_1(nvirt,nocc,nocc,nvirt))
         call omp_reshape(tmp_4_1,tau_tilde,'3124')
         call dgemm_wrapper('N','N',nvirt,nvirt,nocc**2*nvirt,tmp_4_1,cc_int%oovv,F_vv,alpha=0.5_p,beta=1.0_p)
         deallocate(tmp_4_1)

         ! mnef -> nefm
         allocate(tmp_4_1(nocc,nvirt,nvirt,nocc))
         call omp_reshape(tmp_4_1,cc_int%oovv,'2341')
         call dgemm_wrapper('N','N',nocc,nocc,nocc*nvirt**2,tau_tilde,tmp_4_1,F_oo,alpha=0.5_p,beta=1.0_p)
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

         use linalg, only: dgemm_wrapper, omp_reshape

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int
         integer :: m, n, j, e, f, b
         real(p), dimension(:,:,:,:), allocatable :: scratch, reshape_scratch, tmp

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
         call omp_reshape(reshape_scratch, scratch, '1243')
         W_oooo = cc_int%oooo + scratch - reshape_scratch
         deallocate(reshape_scratch)

         ! 1/2 tau_ijef <mn||ef>
         allocate(reshape_scratch(nvirt,nvirt,nocc,nocc))
         call omp_reshape(reshape_scratch, tau, '3412')
         call dgemm_wrapper('N','N',nocc**2,nocc**2,nvirt**2,cc_int%oovv,reshape_scratch,W_oooo,alpha=0.5_dp,beta=1.0_dp)
         deallocate(reshape_scratch)

         ! We reshape it to W_ijmn to make contractions more amenable to dgemm
         allocate(reshape_scratch, mold=W_oooo)
         call omp_reshape(reshape_scratch, W_oooo, '3412')
         W_oooo = reshape_scratch
         deallocate(scratch,reshape_scratch)

         ! #########################################################################################################################
         ! Eq. (7): W_abef = <ab||ef> - P_(ab)(t_mb<am||ef>)
         ! #########################################################################################################################
         ! - P_(ab) t_m^b <am||ef> = + P_(ab) (t_b^m)^T <ma||ef>
         allocate(scratch,reshape_scratch, mold=W_vvvv)
         call dgemm_wrapper('T','N',nvirt,nvirt**3,nocc,t1,cc_int%ovvv,scratch)
         call omp_reshape(reshape_scratch, scratch, '2134')
         ! The reshaped tensor is the unpermuted one!
         W_vvvv = cc_int%vvvv + reshape_scratch - scratch

         ! We reshape it to W_efab to make contractions more amenable to dgemm
         call omp_reshape(reshape_scratch, W_vvvv, '3412')
         W_vvvv = reshape_scratch
         deallocate(scratch,reshape_scratch)

         ! #########################################################################################################################
         ! Eq. (8): W_mbej = <mb||ej> + t_jf<mb||ef> - t_nb<mn||ej> - (1/2 t_jnfb + t_jf t_nb)<mn||ef>
         ! #########################################################################################################################
         ! <mb||ej> + t_jf<mb||ef>
         call dgemm_wrapper('N','T',nocc*nvirt**2,nocc,nvirt,cc_int%ovvv,t1,W_ovvo)
         W_ovvo = W_ovvo + cc_int%ovvo

         ! - t_nb<mn||ej> = +(t_bn)^T<nm||ej>
         allocate(reshape_scratch, mold=W_ovvo)
         allocate(scratch(nvirt,nocc,nvirt,nocc))
         call dgemm_wrapper('T','N',nvirt,nocc**2*nvirt,nocc,t1,cc_int%oovo,scratch)
         ! We need to reshape it because currently the scratch tensor is b,m,e,j
         call omp_reshape(reshape_scratch, scratch, '2134')
         W_ovvo = W_ovvo + reshape_scratch
         deallocate(scratch,reshape_scratch)


         ! - (1/2 t_jnfb + t_jf t_nb)<mn||ef>
         ! - (fnjb) <me||fn>, then reshape
         allocate(scratch(nocc,nvirt,nocc,nvirt))
         allocate(reshape_scratch(nocc,nvirt,nocc,nvirt))
         call omp_reshape(reshape_scratch, cc_int%oovv, '1324')
         !$omp parallel default(none) &
         !$omp shared(cc_amp,cc_int, scratch)
         !$omp do schedule(static, 10) collapse(4)
         do b = 1, nvirt
            do j = 1, nocc
               do f = 1, nvirt
                  do n = 1, nocc
                     scratch(n,f,j,b) = 0.5*t2(j,n,f,b) + t1(j,f)*t1(n,b)
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel
         allocate(tmp(nocc,nvirt,nocc,nvirt))
         call dgemm_wrapper('N','N',nocc*nvirt,nocc*nvirt,nocc*nvirt,reshape_scratch,scratch,tmp)
         deallocate(reshape_scratch, scratch)
         allocate(reshape_scratch(nocc,nvirt,nvirt,nocc))
         call omp_reshape(reshape_scratch, tmp, '1423')
         W_ovvo = W_ovvo - reshape_scratch
         deallocate(reshape_scratch,tmp)

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

         integer :: i, j, a, b, m, e, f, n

         associate(tmp_t1=>cc_int%tmp_tia,tmp_t2=>cc_int%tmp_tijab,t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab,asym=>int_store%asym_spinorb, &
            F_oo=>cc_int%F_oo,F_ov=>cc_int%F_ov,F_vv=>cc_int%F_vv,D_ia=>cc_int%D_ia,D_ijab=>cc_int%D_ijab,W_oooo=>cc_int%W_oooo, &
            W_ovvo=>cc_int%W_ovvo,tau=>cc_int%tau,tau_tilde=>cc_int%tau_tilde,nocc=>sys%nocc,nvirt=>sys%nvirt,&
            W_vvvv=>cc_int%W_vvvv)

         ! #########################################################################################################################
         ! Update T1
         ! Eq. (1) D_ia t_i^a = t_ie F_ae - t_ma F_mi + t_imae F_me - t_nf <na||if> -1/2 t_imef <ma||ef> - 1/2 t_mnae <nm||ei>
         ! #########################################################################################################################

         ! #### t_ie F_ae
         call dgemm_wrapper('N','T',nocc,nvirt,nvirt,t1,F_vv,tmp_t1)
         ! #### -t_ma F_mi
         call dgemm_wrapper('T','N',nocc,nvirt,nocc,F_oo,t1,tmp_t1,alpha=-1.0_p,beta=1.0_p)

         !$omp parallel default(none) &
         !$omp shared(sys, cc_amp, int_store, cc_int)
         ! tmp_t1 can be shared as each loop iteration will update a unique element of it.
         !$omp do schedule(static, 10) collapse(2)
         do a = 1, nvirt
            do i = 1, nocc
               do e = 1, nvirt
                  do m = 1, nocc
                     tmp_t1(i,a) = tmp_t1(i,a) + t1(m,e)*cc_int%ovvo(m,a,e,i) + t2(m,i,e,a)*F_ov(m,e)
                     do f = 1, nvirt
                        tmp_t1(i,a) = tmp_t1(i,a) + 0.5*t2(m,i,f,e)*cc_int%ovvv(m,a,f,e) 
                     end do
                     do n = 1, nocc
                        tmp_t1(i,a) = tmp_t1(i,a) - 0.5*t2(m,n,e,a)*cc_int%oovo(m,n,e,i)
                     end do
                  end do
               end do
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
         tmp_t2_s = 0.0_dp
         !$omp parallel default(none) &
         !$omp shared(sys, cc_amp, int_store, cc_int, tmp_t2_s)
         !$omp do schedule(static, 10) collapse(4)
         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     do e = 1, nvirt
                        do m = 1, nocc
                           tmp_t2_s(i,j,a,b) = tmp_t2_s(i,j,a,b) - t1(i,e)*t1(m,a)*cc_int%ovvo(m,b,e,j) +t2(m,i,e,a)*W_ovvo(m,b,e,j)
                        end do 
                     end do
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel
         call omp_reshape(reshape_tmp, tmp_t2_s, '2134')
         tmp_t2 = tmp_t2 + tmp_t2_s - reshape_tmp
         call omp_reshape(reshape_tmp, tmp_t2_s, '1243')
         tmp_t2 = tmp_t2 - reshape_tmp
         call omp_reshape(reshape_tmp, tmp_t2_s, '2143')
         tmp_t2 = tmp_t2 + reshape_tmp

         ! P_(ab)(t_ijae F_be)
         call dgemm_wrapper('N','T',nocc**2*nvirt,nvirt,nvirt,t2,F_vv,tmp_t2_s)
         call omp_reshape(reshape_tmp, tmp_t2_s, '1243')
         tmp_t2 = tmp_t2 + tmp_t2_s - reshape_tmp
         ! -1/2 P_(ab)(t_ijae t_mb F_me)
         call dgemm_wrapper('N','T',nocc**2*nvirt,nvirt,nvirt,t2,matmul(transpose(t1),F_ov),tmp_t2_s)
         call omp_reshape(reshape_tmp, tmp_t2_s, '1243')
         tmp_t2 = tmp_t2 - (tmp_t2_s - reshape_tmp)/2
         ! -1/2 P_(ij)(t_je F_me t_imab)
         call dgemm_wrapper('N','N',nocc,nocc*nvirt**2,nocc,matmul(t1,transpose(F_ov)),t2,tmp_t2_s)
         call omp_reshape(reshape_tmp, tmp_t2_s, '2134')
         tmp_t2 = tmp_t2 - (tmp_t2_s - reshape_tmp)/2
         ! P_(ij)(t_ie<ab||ej>) = P_(ij)(t_ie<ej||ab>)
         call dgemm_wrapper('N','N',nocc,nocc*nvirt**2,nvirt,t1,cc_int%vovv,tmp_t2_s)
         call omp_reshape(reshape_tmp, tmp_t2_s, '2134')
         tmp_t2 = tmp_t2 + tmp_t2_s - reshape_tmp
         ! -P_(ab)(t_ma <mb||ij>) = +P_(ab)(<ij||bm> t_ma)
         call dgemm_wrapper('N','N',nocc**2*nvirt,nvirt,nocc,cc_int%oovo,t1,tmp_t2_s)
         call omp_reshape(reshape_tmp, tmp_t2_s, '1243')
         ! The reshaped tensor is the unpermuted one, since tmp_t2_s is i,j,b,a
         tmp_t2 = tmp_t2 + reshape_tmp - tmp_t2_s
         ! -P_(ij)(t_imab F_mj)
         call dgemm_wrapper('T','N',nocc,nocc*nvirt**2,nocc,F_oo,t2,tmp_t2_s)
         call omp_reshape(reshape_tmp, tmp_t2_s, '2134')
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

      subroutine update_restricted_intermediates(sys, cc_amp, cc_int)

         ! Build the "recursively generated" intermediates defined in Piecuch et al Table 1.

         use linalg, only: dgemm_wrapper

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         integer :: i, j, a, b, c
         real(dp), allocatable, dimension(:,:,:,:) :: reshape_tmp, scratch

         associate(tmp_t1=>cc_int%tmp_tia,tmp_t2=>cc_int%tmp_tijab,t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab, &
            D_ia=>cc_int%D_ia,D_ijab=>cc_int%D_ijab, nocc=>sys%nocc, nvirt=>sys%nvirt, nbasis=>sys%nbasis, &
            I_vv=>cc_int%I_vv, I_oo=>cc_int%I_oo, I_vo=>cc_int%I_vo, I_oo_p=>cc_int%I_oo_p, &
            I_oooo=>cc_int%I_oooo, I_ovov=>cc_int%I_ovov, I_voov=>cc_int%I_voov, &
            I_vovv_p=>cc_int%I_vovv_p, I_ooov_p=>cc_int%I_ooov_p, &
            v_oovv=>cc_int%v_oovv, v_ovov=>cc_int%v_ovov, v_vvov=>cc_int%v_vvov, v_vvvv=>cc_int%v_vvvv, v_oooo=>cc_int%v_oooo, &
            v_oovo=>cc_int%v_oovo, asym_t2=>cc_int%asym_t2, c_oovv=>cc_int%c_oovv, x_voov=>cc_int%x_voov)

         ! ----------------------------------------------------------------
         ! asym_t2 = 2*t_miea - t_imea
         asym_t2 = reshape(t2, shape(asym_t2), order=(/2,1,3,4/))
         asym_t2 = -asym_t2 + 2*t2

         ! ----------------------------------------------------------------
         ! c_ijab = t_ijab + t_ia t_jb
         !$omp parallel default(none) shared(sys, cc_amp, cc_int) 

         !$omp do collapse(2) schedule(static, 10)
         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     c_oovv(i,j,a,b) = t2(i,j,a,b) + t1(i,a)*t1(j,b)
                  end do
               end do
            end do
         end do
         !$omp end do

         ! ----------------------------------------------------------------
         ! I_ai = (2 v_aeim - v_eaim) t_me
         ! v_aeim = v_oovv(m,i,e,a), v_eaim = v_oovv(m,i,a,e)
         !$omp do collapse(2) schedule(static, 10) 
         do i = 1, nocc
            do a = 1, nvirt
               I_vo(a,i) = sum((2*v_oovv(:,i,:,a) - v_oovv(:,i,a,:))*t1)
            end do 
         end do
         !$omp end do

         ! ----------------------------------------------------------------
         ! I_ba = (2 v_beam - v_bema) t_me - (2 v_ebmn - v_bemn) c_mnea
         !  v_beam = v_vvov(e,b,m,a), v_bema = v_vvov(b,e,m,a)
         !  v_ebmn = v_oovv(m,n,e,b), v_bemn = v_oovv(m,n,b,e)
         !$omp do collapse(2) schedule(static, 10) 
         do a = 1, nvirt
            do b = 1, nvirt
               I_vv(b,a) = sum((2*v_vvov(:,b,:,a) - v_vvov(b,:,:,a))*transpose(t1)) &
                         - sum((2*v_oovv(:,:,:,b) - v_oovv(:,:,b,:))*c_oovv(:,:,:,a))
            end do 
         end do
         !$omp end do

         ! ----------------------------------------------------------------
         ! I_ji' = (2 v_jeim - v_ejim) t_me + (v_efmi - v_efim) t_mjef
         !  v_jeim = v_oovo(m,i,e,j), v_ejim = v_oovo(i,m,e,j)
         !  v_efmi = v_oovv(m,i,e,f), v_efim = v_oovv(m,i,f,e)
         !  to align the indices the second term has to be broken into two sums, and t_mjef = t_jmfe is used
         !$omp do collapse(2) schedule(static, 10)
         do i = 1, nocc
            do j = 1, nocc
               I_oo_p(j,i) = sum((2*v_oovo(:,i,:,j) - v_oovo(i,:,:,j))*t1) &
                         + sum(v_oovv(:,i,:,:)*t2(:,j,:,:)) - sum(v_oovv(:,i,:,:)*t2(j,:,:,:))
            end do 
         end do
         !$omp end do

         !$omp end parallel

         ! ----------------------------------------------------------------
         ! I_ji = I_ji' + I_ei t_je
         ! When it gets big of course we can use dgemm
         ! [todo]: conditional switches?
         I_oo = I_oo_p + matmul(t1, I_vo)

         ! ----------------------------------------------------------------
         ! I_klij = v_klij + v_efij c_klef + P(ik/jl) t_ke v_elij
         !     v_efij <= reshape v_oovv(i,j,e,f) into efij then we have a dgemm
         !     v_elij <= reshape v_oovo(i,l,e,j) into elij
         I_oooo = v_oooo
         allocate(reshape_tmp(nvirt, nvirt, nocc, nocc))
         reshape_tmp = reshape(v_oovv, shape(reshape_tmp), order=(/3,4,1,2/))
         call dgemm_wrapper('N','N',nocc**2,nocc**2,nvirt**2,c_oovv,reshape_tmp,I_oooo,beta=1.0_dp)
         deallocate(reshape_tmp)
         allocate(reshape_tmp(nvirt,nocc,nocc,nocc))
         reshape_tmp = reshape(v_oovo, shape(reshape_tmp), order=(/3,2,1,4/))
         allocate(scratch(nocc,nocc,nocc,nocc))
         call dgemm_wrapper('N','N',nocc,nocc**3,nvirt,t1,reshape_tmp,scratch)
         deallocate(reshape_tmp)
         allocate(reshape_tmp(nocc,nocc,nocc,nocc))
         reshape_tmp = reshape(scratch, shape(reshape_tmp), order=(/2,1,4,3/))
         I_oooo = I_oooo + scratch + reshape_tmp
         deallocate(reshape_tmp, scratch) 

         ! ----------------------------------------------------------------
         ! I_jbia = v_jbia - 1/2 v_ebim c_jmea - v_jbim t_ma + v_ebia t_je
         ! v_jbia = v_ovov
         ! v_ebim c_jmea = v_oovv(m,i,b,e) * c_oovv(m,j,a,e)
         ! v_jbim <= v_oovo(m,i,b,j) reshape to jbim, then dgemm with t_ma
         ! v_ebia = v_vvov, direct dgemm
         
         I_ovov = v_ovov
         
         !$omp parallel default(none) shared(sys, cc_amp, cc_int) 

         !$omp do collapse(2) schedule(static, 10)
         do a = 1, nvirt
            do i = 1, nocc
               do b = 1, nvirt
                  do j = 1, nocc
                     I_ovov(j,b,i,a) = I_ovov(j,b,i,a) - 0.5*sum(v_oovv(:,i,b,:)*c_oovv(:,j,a,:))
                  end do
               end do
            end do 
         end do
         !$omp end do

         !$omp end parallel

         allocate(reshape_tmp(nocc,nvirt,nocc,nocc))
         reshape_tmp = reshape(v_oovo, shape(reshape_tmp), order=(/4,3,2,1/))
         call dgemm_wrapper('N','N',nocc**2*nvirt,nvirt,nocc,reshape_tmp,t1,I_ovov,alpha=-1.0_dp,beta=1.0_dp)
         deallocate(reshape_tmp)
         call dgemm_wrapper('N','N',nocc,nocc*nvirt**2,nvirt,t1,v_vvov,I_ovov,beta=1.0_dp)

         ! ----------------------------------------------------------------
         ! I_bjia = v_bjia + v_beim (t_mjea - 0.5 c_mjae) - 0.5 v_mbie t_jmae + v_beia t_je - v_bjim t_ma
         ! v_bjia = v_oovv(j,i,a,b)
         ! v_beim = v_oovv(m,i,e,b)
         ! v_mbie = v_ovov(m,e,i,b), both contracted indices in front
         ! v_beia = v_vvov(b,e,i,a), can't really help
         ! v_jbim t_ma = v_oovo(m,i,b,j) * t1(m,a) (contracted indices at the front)

         !$omp parallel default(none) shared(sys, cc_amp, cc_int) 

         !$omp do collapse(2) schedule(static, 10)
         do a = 1, nvirt
            do i = 1, nocc
               do b = 1, nvirt
                  do j = 1, nocc
                     I_voov(b,j,i,a) = v_oovv(j,i,a,b) + sum(v_oovv(:,i,:,b)*(t2(:,j,:,a)-0.5*c_oovv(:,j,a,:))) - &
                     0.5*sum(v_ovov(:,:,i,b)*t2(:,j,:,a)) + sum(v_vvov(b,:,i,a)*t1(j,:)) - sum(v_oovo(i,:,b,j)*t1(:,a))
                  end do
               end do
            end do
         end do
         !$omp end do

         ! ----------------------------------------------------------------
         ! I_ciab' = v_ciab - v_ciam t_mb - t_ma v_cimb
         ! v_ciab = v_vovv = v_vvov(b,a,i,c)
         ! v_ciam = v_vovo = v_ovov(m,a,i,c)
         ! v_cimb = v_voov = v_oovv(m,i,c,b)
         !$omp do collapse(2) schedule(static, 10)
         do b = 1, nvirt
            do a = 1, nvirt
               do i = 1, nocc
                  do c = 1, nvirt
                     I_vovv_p(c,i,a,b) = v_vvov(b,a,i,c) - sum(v_ovov(:,a,i,c)*t1(:,b)) - sum(v_oovv(:,i,c,b)*t1(:,a))
                  end do
               end do
            end do
         end do

         ! ----------------------------------------------------------------
         ! x_bjia = v_beia t_je
         ! v_beia = v_vvov

         !$omp do collapse(2) schedule(static, 10)
         do a = 1, nvirt
            do i = 1, nocc
               do j = 1, nocc
                  do b = 1, nvirt
                     x_voov(b,j,i,a) = sum(v_vvov(b,:,i,a)*t1(j,:))
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel

         ! ----------------------------------------------------------------
         ! I_jkia' = v_jkia + v_efia t_jkef + t_je x_ekia
         ! v_jkia = v_oovo(k,j,a,i), reshape / loop might be similar
         ! v_efia = v_vvov, a nice dgemm finally
         I_ooov_p = reshape(v_oovo, shape(I_ooov_p), order=(/2,1,4,3/))
         call dgemm_wrapper('N','N',nocc**2,nocc*nvirt,nvirt**2,t2,v_vvov,I_ooov_p,beta=1.0_p)
         call dgemm_wrapper('N','N',nocc,nocc**2*nvirt,nvirt,t1,x_voov,I_ooov_p,beta=1.0_p)

         end associate

      end subroutine update_restricted_intermediates

      subroutine update_restricted_intermediates_debug(sys, cc_amp, cc_int)

         ! Build the "recursively generated" intermediates defined in Piecuch et al Table 1.

         use linalg, only: dgemm_wrapper

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         integer :: i, j, a, b, c, k ,l

         associate(tmp_t1=>cc_int%tmp_tia,tmp_t2=>cc_int%tmp_tijab,t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab, &
            D_ia=>cc_int%D_ia,D_ijab=>cc_int%D_ijab, nocc=>sys%nocc, nvirt=>sys%nvirt, nbasis=>sys%nbasis, &
            I_vv=>cc_int%I_vv, I_oo=>cc_int%I_oo, I_vo=>cc_int%I_vo, I_oo_p=>cc_int%I_oo_p, &
            I_oooo=>cc_int%I_oooo, I_ovov=>cc_int%I_ovov, I_voov=>cc_int%I_voov, &
            I_vovv_p=>cc_int%I_vovv_p, I_ooov_p=>cc_int%I_ooov_p, &
            v_oovv=>cc_int%v_oovv, v_ovov=>cc_int%v_ovov, v_vvov=>cc_int%v_vvov, v_vvvv=>cc_int%v_vvvv, v_oooo=>cc_int%v_oooo, &
            v_oovo=>cc_int%v_oovo, asym_t2=>cc_int%asym_t2, c_oovv=>cc_int%c_oovv, x_voov=>cc_int%x_voov)

         do i = 1, nocc
            do a = 1, nvirt
               I_vo(a,i) = sum((2*v_oovv(:,i,:,a)-v_oovv(:,i,a,:))*t1(:,:))
            end do
         end do

         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     c_oovv(i,j,a,b) = t2(i,j,a,b) + t1(i,a)*t1(j,b)
                  end do
               end do
            end do
         end do

         do a = 1, nvirt
            do b = 1, nvirt
               I_vv(b,a) = sum((2*v_vvov(:,b,:,a)-v_vvov(b,:,:,a))*transpose(t1)) &
                - sum((2*v_oovv(:,:,:,b)-v_oovv(:,:,b,:))*c_oovv(:,:,:,a))
            end do
         end do

         do i = 1, nocc
            do j = 1, nocc
               I_oo_p(j,i) = sum((2*v_oovo(:,i,:,j)-v_oovo(i,:,:,j))*t1(:,:)) + sum((v_oovv(:,i,:,:)-v_oovv(i,:,:,:))*t2(:,j,:,:))
            end do
         end do

         do i = 1, nocc
            do j = 1, nocc
               I_oo(j,i) = I_oo_p(j,i) + sum(I_vo(:,i)*t1(j,:))
            end do
         end do

         ! No change

         do j = 1, nocc
            do i = 1, nocc
               do l = 1, nocc
                  do k = 1, nocc
                     I_oooo(k,l,i,j) = v_oooo(k,l,i,j) + sum(v_oovv(i,j,:,:)*c_oovv(k,l,:,:)) + sum(t1(k,:)*v_oovo(i,j,:,l)) &
                     + sum(t1(l,:)*v_oovo(j,i,:,k))
                  end do
               end do
            end do
         end do

         ! No change

         do b = 1, nvirt
            do a = 1, nvirt
               do i = 1, nocc
                  do c = 1, nvirt
                     I_vovv_p(c,i,a,b) = v_vvov(b,a,i,c) - sum(v_ovov(:,a,i,c)*t1(:,b)) - sum(t1(:,a)*v_oovv(:,i,c,b))
                  end do
               end do
            end do
         end do

         ! No change

         do a = 1, nvirt
            do i = 1, nocc
               do b = 1, nvirt
                  do j = 1, nocc
                     I_ovov(j,b,i,a) = v_ovov(j,b,i,a) - sum(v_oovv(i,:,:,b)*c_oovv(j,:,:,a))*0.5 - sum(v_oovo(:,i,b,j)*t1(:,a)) &
                     + sum(v_vvov(:,b,i,a)*t1(j,:))
                  end do
               end do
            end do
         end do

         ! No change

         do a = 1, nvirt
            do i = 1, nocc
               do j = 1, nocc
                  do b = 1, nvirt
                     I_voov(b,j,i,a) = v_oovv(i,j,b,a) + sum(v_oovv(i,:,b,:)*(t2(:,j,:,a)-0.5*c_oovv(:,j,a,:))) &
                     - 0.5*sum(v_ovov(:,b,i,:)*t2(j,:,a,:)) + sum(v_vvov(b,:,i,a)*t1(j,:)) - sum(v_oovo(i,:,b,j)*t1(:,a))
                  end do
               end do
            end do
         end do

         ! No change

         do a = 1, nvirt
            do i = 1, nocc
               do j = 1, nocc
                  do b = 1, nvirt
                     x_voov(b,j,i,a) = sum(v_vvov(b,:,i,a)*t1(j,:))
                  end do
               end do
            end do
         end do

         ! No change

         do a = 1, nvirt
            do i = 1, nocc
               do k = 1, nocc
                  do j = 1, nocc
                     I_ooov_p(j,k,i,a) = v_oovo(k,j,a,i) + sum(v_vvov(:,:,i,a)*t2(j,k,:,:)) + sum(t1(j,:)*x_voov(:,k,i,a))
                  end do
               end do
            end do
         end do

         ! No change

         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     asym_t2(i,j,a,b) = 2*t2(i,j,a,b) - t2(j,i,a,b)
                  end do
               end do
            end do
         end do

         end associate

      end subroutine update_restricted_intermediates_debug

      subroutine update_amplitudes_restricted_debug(sys, cc_amp, cc_int)

         ! Perform the CC amplitude equations for the spin-free formulation
         ! In:
         !     sys: system under study.
         !     int_store: integral information
         ! In/out:
         !     cc_amp: CC amplitudes being updated

         use linalg, only: dgemm_wrapper

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         real(dp), dimension(:,:,:,:), allocatable :: tmp_t2_s
         integer :: i, j, a, b, m, e

         associate(tmp_t1=>cc_int%tmp_tia,tmp_t2=>cc_int%tmp_tijab,t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab, &
            D_ia=>cc_int%D_ia,D_ijab=>cc_int%D_ijab, nocc=>sys%nocc, nvirt=>sys%nvirt, nbasis=>sys%nbasis, &
            I_vv=>cc_int%I_vv, I_oo=>cc_int%I_oo, I_vo=>cc_int%I_vo, I_oo_p=>cc_int%I_oo_p, &
            I_oooo=>cc_int%I_oooo, I_ovov=>cc_int%I_ovov, I_voov=>cc_int%I_voov, &
            I_vovv_p=>cc_int%I_vovv_p, I_ooov_p=>cc_int%I_ooov_p, &
            v_oovv=>cc_int%v_oovv, v_ovov=>cc_int%v_ovov, v_vvov=>cc_int%v_vvov, v_vvvv=>cc_int%v_vvvv, v_oooo=>cc_int%v_oooo, &
            v_oovo=>cc_int%v_oovo, asym_t2=>cc_int%asym_t2, c_oovv=>cc_int%c_oovv, x_voov=>cc_int%x_voov)


         do a = 1, nvirt
            do i = 1, nocc
               tmp_t1(i,a) = sum(I_vv(:,a)*t1(i,:)) - sum(I_oo_p(i,:)*t1(:,a)) + sum(transpose(I_vo)*(2*t2(:,i,:,a)-t2(i,:,:,a))) &
                  + sum(t1*(2*v_oovv(:,i,:,a)-v_ovov(:,a,i,:))) - sum(v_oovo(:,:,:,i)*(2*t2(:,:,:,a)-t2(:,:,a,:))) 
               do m = 1, nocc
                  tmp_t1(i,a) = tmp_t1(i,a) + sum(v_vvov(:,:,m,a)*(2*t2(m,i,:,:)-t2(i,m,:,:)))
               end do
            end do
         end do

         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     tmp_t2(i,j,a,b) = sum(t2(i,j,a,:)*I_vv(:,b)) - sum(t2(i,:,a,b)*I_oo(j,:)) + &
                     0.5*sum(v_vvvv(:,:,a,b)*c_oovv(i,j,:,:)) + 0.5*sum(c_oovv(:,:,a,b)*I_oooo(i,j,:,:)) &
                     + sum(t1(i,:)*I_vovv_p(:,j,a,b)) - sum(t1(:,a)*I_ooov_p(i,j,:,b))
                     do e = 1, nvirt
                        do m = 1, nocc
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - t2(m,j,a,e)*I_ovov(i,e,m,b) - I_ovov(i,e,m,a)*t2(m,j,e,b) + &
                           (2*t2(m,i,e,a)-t2(i,m,e,a))*I_voov(e,j,m,b) 
                        end do
                     end do
                  end do
               end do
            end do
         end do

         allocate(tmp_t2_s(nocc,nocc,nvirt,nvirt))

         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     tmp_t2_s(i,j,a,b) = v_oovv(i,j,a,b) + tmp_t2(i,j,a,b) + tmp_t2(j,i,b,a)
                  end do
               end do
            end do
         end do
         tmp_t2 = tmp_t2_s
         deallocate(tmp_t2_s)

         t1 = tmp_t1/D_ia
         t2 = tmp_t2/D_ijab

         ! All no change. Intermediates and amplitude update code is **correct**

         end associate

      end subroutine update_amplitudes_restricted_debug

      subroutine update_amplitudes_restricted(sys, cc_amp, cc_int)

         ! Perform the CC amplitude equations for the spin-free formulation
         ! In:
         !     sys: system under study.
         ! In/out:
         !     cc_amp: CC amplitudes being updated
         !     cc_int: CC intermediates used in computing the updates

         use linalg, only: dgemm_wrapper

         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         real(dp), dimension(:,:,:,:), allocatable :: reshape_tmp1, tmp_t2_s
         real(dp) :: tmp
         integer :: i, j, a, b, m, e, f

         associate(tmp_t1=>cc_int%tmp_tia,tmp_t2=>cc_int%tmp_tijab,t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab, &
            D_ia=>cc_int%D_ia,D_ijab=>cc_int%D_ijab, nocc=>sys%nocc, nvirt=>sys%nvirt, nbasis=>sys%nbasis, &
            I_vv=>cc_int%I_vv, I_oo=>cc_int%I_oo, I_vo=>cc_int%I_vo, I_oo_p=>cc_int%I_oo_p, &
            I_oooo=>cc_int%I_oooo, I_ovov=>cc_int%I_ovov, I_voov=>cc_int%I_voov, &
            I_vovv_p=>cc_int%I_vovv_p, I_ooov_p=>cc_int%I_ooov_p, &
            v_oovv=>cc_int%v_oovv, v_ovov=>cc_int%v_ovov, v_vvov=>cc_int%v_vvov, v_vvvv=>cc_int%v_vvvv, v_oooo=>cc_int%v_oooo, &
            v_oovo=>cc_int%v_oovo, asym_t2=>cc_int%asym_t2, c_oovv=>cc_int%c_oovv, x_voov=>cc_int%x_voov)

         ! ###################################################################################
         ! ############################### T1 Updates ########################################
         ! ###################################################################################

         ! ------------------ Eq. 43, terms 2-3 ----------------------------------------------
         ! t_ie * I_ea -I'_im * t_ma, simple matmuls
         ! Benchmarking suggests that the Fortran intrinsic matmul beats threaded dgemm, not to mention naive OpenMP. So here we go
         tmp_t1 = matmul(t1, I_vv) - matmul(I_oo_p, t1)

         ! ------------------ Eq. 43, terms 4-5 ----------------------------------------------
         ! I_em ( 2 t_miea - t_imea ) + t_me (2 v_eima - v_eiam)
         !  = I_vo^T(m,e)*asym_t2(m,i,e,a) + 2*t1(m,e)*v_oovv(m,i,e,a) - t1(m,e)*v_ovov(m,a,i,e)
         !$omp parallel do default(none)&
         !$omp schedule(static,10) collapse(2)&
         !$omp shared(sys, cc_amp, cc_int)
         do a = 1, nvirt
            do i = 1, nocc
               ! Hopefully the compiler can cache transpose(I_vo) somewhere so we don't have to create a tmp array
               tmp_t1(i,a) = tmp_t1(i,a) + sum(transpose(I_vo)*asym_t2(:,i,:,a)) &
                           + sum(t1*(2*v_oovv(:,i,:,a) - v_ovov(:,a,i,:)))
            end do
         end do
         !$omp end parallel do

         ! ------------------ Eq. 43, term 6 -------------------------------------------------
         ! - v_ei^mn (2 t_mn^ea - t_mn^ae) 
         ! This is a clear case where we can benefit from dgemm 
         ! (see my benchmarking script at https://github.com/brianz98/fortran-tensor-benchmarking).
         ! As the number of contracted indices becomes large, the (significant) overhead of dgemm becomes negligible.
         !
         ! v_ei^mn would require a tensor of type v_vooo, but we only have v_oovo, which is permutationally equivalent
         ! v_ei^mn would be accessed by v_oovo(m,i,e,n), which is why we need order=2,1,4,3 to permute it into (i,m,n,e)
         ! We would like to not store the reshaped arrays but gfortran doesn't like big arrays on the stack (ifort has the 
         ! handy -heap-arrays flag that directs them to the heap, but apparently there's no such flag in gfortran..)
         ! Anyways, we'll have to resort to using temporary arrays to ensure portability.
         
         ! We further need a scratch t1 matrix as dgemm 'C' is overwritten
         allocate(reshape_tmp1(nocc,nocc,nocc,nvirt))
         reshape_tmp1 = reshape(v_oovo,(/nocc,nocc,nocc,nvirt/),order=(/2,1,4,3/))
         call dgemm_wrapper('N','N',nocc,nvirt,nocc**2*nvirt,reshape_tmp1,asym_t2,tmp_t1,alpha=-1.0_p,beta=1.0_p)
         deallocate(reshape_tmp1)

         ! ------------------ Eq. 43, term 7 -------------------------------------------------
         ! + v_ef^ma (2 t_mi^ef - t_im^ef)
         ! v_efma is type v_vvov, which can be accessed via v_vvov(e,f,m,a), and reshape into (m,e,f,a)
         ! second quantity is the 'asymmetrised t2'(m,i,e,f), which we need to reshape into (i,m,e,f)
         ! Benchmarking shows double reshape is slower than direct OMP
         !$omp parallel do default(none)&
         !$omp schedule(static,10) collapse(2)&
         !$omp shared(sys, cc_amp, cc_int) private(tmp)
         do a = 1, nvirt
            do i = 1, nocc
               tmp = 0.0_dp
               do e = 1, nvirt
                  do f = 1, nvirt
                     do m = 1, nocc
                        tmp = tmp + v_vvov(e,f,m,a)*asym_t2(m,i,e,f)
                     end do
                  end do
               end do
               tmp_t1(i,a) = tmp_t1(i,a) + tmp
            end do
         end do
         !$omp end parallel do

         ! ###################################################################################
         ! ############################### T2 Updates ########################################
         ! ###################################################################################         

         ! ------------------ Eq. 44, term 1 -------------------------------------------------
         ! We delay the addition of v_ijab until the end so we can accumulate the entire bracketed quantity in tmp_t2, 
         ! then apply the permutation directly, so there's nothing here

         ! We have a permutation operator P(ia/jb), after v_ij^ab, meaning i becomes j and a becomes b, 
         ! and conveniently what we actually need to do to perform the action of the permutation is to just reshape
         ! tmp_t2 from (i,j,a,b) to (j,i,b,a) and add it back to the original tmp_t2,

         ! ------------------ Eq. 44, term 2 -------------------------------------------------
         ! t_ij^ae*I_e^b, a nice dgemm
         call dgemm_wrapper('N','N',nocc**2*nvirt,nvirt,nvirt,t2,I_vv,tmp_t2)

         ! ------------------ Eq. 44, term 3 -------------------------------------------------
         ! -t_im^ab*I_j^m, benchmarking shows naive OMP is the fastest
         !$omp parallel do default(none)&
         !$omp schedule(static,10) collapse(2)&
         !$omp shared(sys, cc_amp, cc_int)
         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - sum(t2(:,i,b,a)*I_oo(j,:))
                  end do
               end do
            end do
         end do
         !$omp end parallel do

         ! ------------------ Eq. 44, term 4 -------------------------------------------------
         ! 1/2 v_ef^ab c_ij^ef, nice dgemm again
         call dgemm_wrapper('N','N',nocc**2,nvirt**2,nvirt**2,c_oovv,v_vvvv,tmp_t2,alpha=0.5_p,beta=1.0_p)

         ! ------------------ Eq. 44, term 5 -------------------------------------------------
         ! 1/2 c_mn^ab I_ij^mn, nice dgemm
         call dgemm_wrapper('N','N',nocc**2,nvirt**2,nocc**2,I_oooo,c_oovv,tmp_t2,alpha=0.5_p,beta=1.0_p)

         ! ------------------ Eq. 44, terms 6-8 ----------------------------------------------
         ! -t_mj^ae I_ie^mb - I_ie^ma t_mj^eb + (2t_mi^ea - t_im^ea) I_ej^mb, seems hopeless, use OMP
         !$omp parallel do default(none)&
         !$omp schedule(static,10) collapse(3)&
         !$omp shared(sys, cc_amp, cc_int) private(tmp)
         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     tmp = 0.0_dp
                     do e = 1, nvirt
                        do m = 1, nocc
                           tmp = tmp - t2(m,j,a,e)*I_ovov(i,e,m,b) - I_ovov(i,e,m,a)*t2(m,j,e,b) &
                                 + asym_t2(m,i,e,a)*I_voov(e,j,m,b)
                        end do
                     end do
                     tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) + tmp
                  end do
               end do
            end do
         end do
         !$omp end parallel do

         ! ------------------ Eq. 44, terms 9-10 ---------------------------------------------
         ! t_i^e I'_ej^ab - t_m^a I'_ij^mb, first dgemm, second OMP
         call dgemm_wrapper('N','N',nocc,nocc*nvirt**2,nvirt,t1,I_vovv_p,tmp_t2,beta=1.0_dp)

         !$omp parallel do default(none)&
         !$omp schedule(static,10) collapse(2)&
         !$omp shared(sys, cc_amp, cc_int) 
         do b = 1, nvirt
            do a = 1, nvirt
               do j = 1, nocc
                  do i = 1, nocc
                     tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) - sum(t1(:,a)*I_ooov_p(i,j,:,b))
                  end do
               end do
            end do
         end do
         !$omp end parallel do

         ! P(ia/jb) just means we swap i/j and a/b and add it back to the tensor
         ! Again, it would be nice if we can just say tmp_t2 = reshape(tmp_t2,...) but gfortran can't do it
         allocate(tmp_t2_s(nocc,nocc,nvirt,nvirt))
         tmp_t2_s = reshape(tmp_t2, (/nocc,nocc,nvirt,nvirt/),order=(/2,1,4,3/))
         tmp_t2 = tmp_t2 + tmp_t2_s + v_oovv
         deallocate(tmp_t2_s)

         ! We now write the cluster amplitudes into the actual tensors
         ! Benchmarking shows we don't have to worry about doing this without parallelisation
         t1 = tmp_t1/D_ia
         t2 = tmp_t2/D_ijab

         end associate

      end subroutine update_amplitudes_restricted

      subroutine update_cc_energy(sys, st, cc_int, cc_amp, conv, restricted)
         ! Updates the CC energy and also checks for convergence
         ! In:
         !     sys: system under study.
         !     cc_int: CC intermediates, of which we're using the cc_int%oovv slice.
         !     cc_amp: CC amplitudes.
         !     restricted: whether we're performing a restricted calculation
         ! In/out:
         !     st: state_t object storing energies and RMS(T) info for convergence checking
         ! Out:
         !     conv: returns true if converged within tolerance

         use system, only: state_t
         use integrals, only: int_store_t

         type(system_t), intent(in) :: sys
         type(cc_int_t), intent(in) :: cc_int
         type(cc_amp_t), intent(in) :: cc_amp
         logical, intent(in) :: restricted
         type(state_t), intent(inout) :: st
         logical, intent(out) :: conv

         integer :: i, j, a, b
         real(p) :: ecc, rmst2

         conv = .false.
         st%energy_old = st%energy
         associate(v_oovv=>cc_int%v_oovv, c_oovv=>cc_int%c_oovv, nocc=>sys%nocc, nvirt=>sys%nvirt, &
            oovv=>cc_int%oovv, t1=>cc_amp%t_ia, t2=>cc_amp%t_ijab)

         if (restricted) then
            ecc = 0.0_p
            rmst2 = 0.0_p
            !$omp parallel default(none) &
            !$omp shared(sys,st,cc_amp,ecc,rmst2,cc_int) 
            !$omp do schedule(static, 10) collapse(2) reduction(+:ecc,rmst2)
            do b = 1, nvirt
               do a = 1, nvirt
                  do j = 1, nocc
                     do i = 1, nocc
                        ecc = ecc + (2*v_oovv(i,j,a,b)-v_oovv(i,j,b,a))*(t2(i,j,a,b)+t1(i,a)*t1(j,b))
                        ! We only do the RMS on T2 amplitudes
                        rmst2 = rmst2 + (t2(i,j,a,b)-st%t2_old(i,j,a,b))**2
                     end do
                  end do
               end do
            end do
            !$omp end do
            !$omp end parallel
         else
            ecc = 0.0_p
            rmst2 = 0.0_p
            !$omp parallel default(none) &
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
         end if
         st%energy = ecc
         st%t2_old = t2
         if (sqrt(rmst2) < sys%ccsd_t_tol .and. abs(st%energy-st%energy_old) < sys%ccsd_e_tol) conv = .true.

         end associate

      end subroutine update_cc_energy

      subroutine do_ccsd_t_spinorb(sys, int_store, int_store_cc)
         ! Spinorbital formuation of CCSD(T)
         ! In:
         !     int_store: integral information.
         ! In/out:
         !     sys: holds converged amplitudes from CCSD

         use integrals, only: int_store_t, int_store_cc_t
         use error_handling, only: check_allocate

         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(in) :: int_store
         type(int_store_cc_t), intent(in) :: int_store_cc

         integer :: i, j, k, a, b, c, ierr
         real(p), dimension(:,:,:), allocatable :: tmp_t3d, tmp_t3c, tmp_t3c_d, reshape_tmp1, reshape_tmp2, t2_reshape(:,:,:,:)
         real(p) :: e_T

         write(iunit, '(1X, 10("-"))')
         write(iunit, '(1X, A)') 'CCSD(T)'
         write(iunit, '(1X, 10("-"))')

         associate(nvirt=>sys%nvirt, nocc=>sys%nocc, e=>sys%canon_levels_spinorb, t1=>int_store_cc%t1, t2=>int_store_cc%t2,&
                   vvoo=>int_store_cc%vvoo,vovv=>int_store_cc%vovv,ovoo=>int_store_cc%ovoo)

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
         write(iunit, '(1X, A, 1X, F15.9)') 'Unrestricted CCSD(T) correlation energy (Hartree):', e_T
         end associate

      end subroutine do_ccsd_t_spinorb

      subroutine do_ccsd_t_spatial(sys, int_store, calcname, int_store_cc)

         ! Spatial formuation of CCSD(T) according
         ! P. Piecuch et al., Computer Physics Communications 149 (2002) 71–96, 
         ! https://doi.org/10.1016/S0010-4655(02)00598-2
         ! In/out:
         !     sys: holds converged amplitudes from CCSD
         !     int_store: integral information.
         !     int_store_cc: intermediates carried over from the CCSD calculation.
         ! Out:
         !     calcname: the specific version of CCSD(T)/[T] we're doing

         use integrals, only: int_store_t, int_store_cc_t
         use error_handling, only: check_allocate

         type(int_store_cc_t), intent(inout) :: int_store_cc
         type(system_t), intent(inout) :: sys
         type(int_store_t), intent(inout) :: int_store
         character(*) :: calcname

         integer :: i, j, k, a, b, c, ierr
         real(p), dimension(:,:,:), allocatable :: tmp_t3_D, tmp_t3, z3, t_bar, reshape_tmp, tmp_y, z3_bar, tmp_m3
         real(p), dimension(:,:,:,:), allocatable :: t2_reshape, v_vovv, v_ovoo, asym_t2, c_oovv
         real(p) :: e_T, D_T, e_CR

         write(iunit, '(1X, 10("-"))')
         write(iunit, '(1X, A)') 'CCSD(T)'
         write(iunit, '(1X, 10("-"))')

         associate(nvirt=>sys%nvirt, nocc=>sys%nocc, e=>sys%canon_levels, t1=>int_store_cc%t1, t2=>int_store_cc%t2,&
                   v_oovv=>int_store_cc%v_oovv,v_vvov=>int_store_cc%v_vvov,v_oovo=>int_store_cc%v_oovo, &
                   doing_T=>sys%ccsd_t_paren, doing_R=>sys%ccsd_t_renorm, doing_CR=>sys%ccsd_t_comp_renorm, &
                   I_vovv_pp=>int_store_cc%I_vovv_pp, I_ooov_pp=>int_store_cc%I_ooov_pp)

         ! For efficient cache access during tmp_t3c updates
         allocate(t2_reshape(nvirt,nvirt,nocc,nocc), source=0.0_p, stat=ierr)
         call check_allocate('t2_reshape', (nvirt*nocc)**2, ierr)
         ! ijab -> baij
         t2_reshape = reshape(t2, shape(t2_reshape), order=(/4,3,2,1/))

         ! Reshape some integral slices so the contracted indices end up first
         allocate(v_vovv(nvirt,nocc,nvirt,nvirt),source=0.0_dp,stat=ierr)
         call check_allocate('v_vovv', nocc*nvirt**3, ierr)
         v_vovv = reshape(v_vvov,shape(v_vovv),order=(/4,3,2,1/))
         deallocate(int_store_cc%v_vvov)

         allocate(v_ovoo(nocc,nvirt,nocc,nocc),source=0.0_dp,stat=ierr)
         call check_allocate('v_ovoo', nocc**3*nvirt, ierr)
         v_ovoo = reshape(v_oovo,shape(v_ovoo),order=(/4,3,2,1/))
         deallocate(int_store_cc%v_oovo)

         ! ############################ (R)CCSD(T) or [T] #################################
         if (doing_R .or. doing_CR) then
            allocate(asym_t2(nocc,nocc,nvirt,nvirt), source=0.0_p, stat=ierr)
            call check_allocate('asym_t2', (nocc*nvirt)**2, ierr)
            asym_t2 = reshape(t2,shape(asym_t2),order=(/2,1,3,4/))
            asym_t2 = -asym_t2 + 2*t2

            allocate(c_oovv(nocc,nocc,nvirt,nvirt), source=0.0_p, stat=ierr)
            call check_allocate('c_oovv', (nocc*nvirt)**2, ierr)
            c_oovv = t2 
         end if

         e_T = 0.0_p
         D_T = 0.0_p
         e_CR = 0.0_p
         
         !$omp parallel default(none) &
         !$omp private(tmp_t3_D, tmp_t3, z3, z3_bar, tmp_y, t_bar, ierr, reshape_tmp, tmp_m3) &
         !$omp shared(sys, int_store, t2_reshape, v_vovv, v_ovoo, c_oovv, e_T, D_T, e_CR) 

         if (doing_R .or. doing_CR) then
            ! Both need the denominator, which require c_oovv
            !$omp do collapse(2) schedule(static, 10)
            do b = 1, nvirt
               do a = 1, nvirt
                  do j = 1, nocc
                     do i = 1, nocc
                        c_oovv(i,j,a,b) = c_oovv(i,j,a,b) + t1(i,a)*t1(j,b)
                     end do
                  end do
               end do
            end do
            !$omp end do
         end if

         ! Declare heap-allocated arrays inside parallel region
         allocate(tmp_t3_D(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('tmp_t3_D', nvirt**3, ierr)
         allocate(tmp_t3(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('tmp_t3', nvirt**3, ierr)
         allocate(reshape_tmp(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('reshape_tmp', nvirt**3, ierr)

         if (doing_T) then
            allocate(z3(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
            call check_allocate('z3', nvirt**3, ierr)
            allocate(t_bar(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
            call check_allocate('t_bar', nvirt**3, ierr)
         end if
         
         if (doing_R .or. doing_CR) then
            allocate(tmp_y(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
            call check_allocate('tmp_y', nvirt**3, ierr)
            allocate(z3_bar(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
            call check_allocate('z3_bar', nvirt**3, ierr)
         end if

         if (doing_CR) then
            allocate(tmp_m3(nvirt,nvirt,nvirt), source=0.0_p, stat=ierr)
            call check_allocate('tmp_m3', nvirt**3, ierr)
         end if

         ! We use avoid the storage of full triples (six-dimensional array..) and instead use the strategy of 
         ! batched triples storage, denoted W^{ijk}(abc) in (https://doi.org/10.1016/0009-2614(91)87003-T).
         ! We could of course only compute one element in a loop iteration but that will probably result in bad floating point
         ! performance, here for each thread/loop iteration we use the Fortran intrinsic sum, 
         ! instead of all relying on OMP reduction, to hopefully give better floating point performance.

         ! Based on Eqs. 55-56 and 60, we need three arrays, z3 (z_abc^ijk), t3 (t_ijk^abc(2)), and t3_D (t_ijk^abc(2)*D_ijk^abc)
         ! Eq. 52: E(T) = (t_bar_abc^ijk + z_abc^ijk) t_ijk^abc(2) D_ijk^abc
         ! t_bar_abcijk = 4/3 t_ijkabc -2t_ijkacb + 2/3 t_ijkbca
         ! However, as all indices are contracted, the ordering of the indices doesn't matter
         !$omp do schedule(static, 10) collapse(3) reduction(+:e_T, D_T, e_CR)
         do i = 1, nocc
            do j = 1, nocc
               do k = 1, nocc
                  ! Compute batched T3 amplitudes
                  do a = 1, nvirt
                     do b = 1, nvirt
                        do c = 1, nvirt
                           ! Note that in Piecuch et al. there's mention of complex conjugates (x_abc^ijk = x_ijk^abc*), but 
                           ! since everything's strictly real here we can just follow the normal ordering of indices
                           ! i.e., _ijk^abc
                           !
                           ! We would like to ideally use the intrinsic `reshape`, as we could do for the spin-orbital
                           ! formulation of CCSD(T), unfortunately in that case the reshape involved was a P(i/jk)P(a/bc),
                           ! which was natually factorisable into ijk and abc, so we can reshape only the abc slices while 
                           ! holding ijk fixed. However here the transpositions are 'correlated' and aren't factorisable.
                           ! There might be a way to do it but I can't think of one for now. Anyway the cost of reshapes might 
                           ! be similar to the vector operations anyway(?)
                           tmp_t3_D(a,b,c) = sum(t2_reshape(:,a,j,i)*v_vovv(:,k,b,c)) - sum(t2(:,i,b,a)*v_ovoo(:,c,j,k)) + &
                                             sum(t2_reshape(:,b,i,j)*v_vovv(:,k,a,c)) - sum(t2(:,j,a,b)*v_ovoo(:,c,i,k)) + &
                                             sum(t2_reshape(:,c,j,k)*v_vovv(:,i,b,a)) - sum(t2(:,k,b,c)*v_ovoo(:,a,j,i)) + &
                                             sum(t2_reshape(:,a,k,i)*v_vovv(:,j,c,b)) - sum(t2(:,i,c,a)*v_ovoo(:,b,k,j)) + &
                                             sum(t2_reshape(:,b,k,j)*v_vovv(:,i,c,a)) - sum(t2(:,j,c,b)*v_ovoo(:,a,k,i)) + &
                                             sum(t2_reshape(:,c,i,k)*v_vovv(:,j,a,b)) - sum(t2(:,k,a,c)*v_ovoo(:,b,i,j))
                           tmp_t3(a,b,c) = tmp_t3_D(a,b,c)/(e(i)+e(j)+e(k)-e(a+nocc)-e(b+nocc)-e(c+nocc))
                           if (doing_T) then
                              ! The z3 array as defined in eq.60 of Piecuch
                              z3(a,b,c) = (t1(i,a)*v_oovv(j,k,b,c) + t1(j,b)*v_oovv(i,k,a,c) + &
                                           t1(k,c)*v_oovv(i,j,a,b))/(e(i)+e(j)+e(k)-e(a+nocc)-e(b+nocc)-e(c+nocc))
                           end if
                           if (doing_R .or. doing_CR) then
                              ! The 'y' array as defined in eq. 66 of Piecuch
                              tmp_y(a,b,c) = t1(i,a)*t1(j,b)*t1(k,c) + t1(i,a)*t2(j,k,b,c) + &
                                             t1(j,b)*t2(i,k,a,c) + t1(k,c)*t2(i,j,a,b)
                           end if
                           if (doing_CR) then
                              ! The 'M' array defined in eq. 62 of Piecuch
                              tmp_m3(a,b,c) = sum(t2(i,j,a,:)*I_vovv_pp(:,k,b,c)) - sum(t2(:,i,b,a)*I_ooov_pp(j,k,:,c)) + &
                                              sum(t2(j,i,b,:)*I_vovv_pp(:,k,a,c)) - sum(t2(:,j,a,b)*I_ooov_pp(i,k,:,c)) + &
                                              sum(t2(k,j,c,:)*I_vovv_pp(:,i,b,a)) - sum(t2(:,k,b,c)*I_ooov_pp(j,i,:,a)) + &
                                              sum(t2(i,k,a,:)*I_vovv_pp(:,j,c,b)) - sum(t2(:,i,c,a)*I_ooov_pp(k,j,:,b)) + &
                                              sum(t2(j,k,b,:)*I_vovv_pp(:,i,c,a)) - sum(t2(:,j,c,b)*I_ooov_pp(k,i,:,a)) + &
                                              sum(t2(k,i,c,:)*I_vovv_pp(:,j,a,b)) - sum(t2(:,k,a,c)*I_ooov_pp(i,j,:,b))
                           end if
                        end do
                     end do
                  end do
                  ! Do reshapes for t_bar
                  t_bar = 4*tmp_t3/3
                  reshape_tmp = reshape(tmp_t3,shape(t_bar),order=(/1,3,2/))
                  t_bar = t_bar - 2*reshape_tmp
                  reshape_tmp = reshape(tmp_t3,shape(t_bar),order=(/3,1,2/))
                  t_bar = t_bar + 2*reshape_tmp/3

                  ! Summary of the family of (T) corrections:
                  ! E_[T] = t_bar*t3_D
                  ! E_(T) = (z3 + t_bar)*t3_D
                  ! E_R[T] = E[T] / (1 + asym_t2*c_oovv + t_bar*y)
                  ! E_R(T) = E(T) / (1 + asym_t2*c_oovv + t_bar*y + z3_bar*y)
                  ! E_CR[T] = t_bar*m3 / (1 + asym_t2*c_oovv + t_bar*y)
                  ! E_CR(T) = (t_bar + z3_bar)*m3 / (1 + asym_t2*c_oovv + t_bar*y + z3_bar*y)

                  if (doing_T .and. (doing_R .or. doing_CR)) then
                     ! (C)RCCSD(T) involves z3_bar whereas (C)RCCSD[T] does not
                     ! If we're doing (C)RCCSD(T)
                     z3_bar = 4*z3/3
                     reshape_tmp = reshape(z3,shape(z3_bar),order=(/1,3,2/))
                     z3_bar = z3_bar - 2*reshape_tmp
                     reshape_tmp = reshape(z3,shape(z3_bar),order=(/3,1,2/))
                     z3_bar = z3_bar + 2*reshape_tmp/3
                  end if
                  
                  if (.not. doing_CR) then
                     ! for R/CCSD(T)/[T] we always need the base quantity of t_bar*t3_D
                     e_T = e_T + sum(t_bar*tmp_t3_D)
                     if (doing_T) e_T = e_T + sum(z3*tmp_t3_D)
                  else
                     e_CR = e_CR + sum(t_bar*tmp_m3)
                     if (doing_T) e_CR = e_CR + sum(z3_bar*tmp_m3)
                  end if

                  if (doing_R .or. doing_CR) then
                     ! Separately accumulate the denominator
                     D_T = D_T + sum(t_bar*tmp_y)
                     if (doing_T) D_T = D_T + sum(z3_bar*tmp_y)
                  end if

               end do
            end do
         end do
         !$omp end do
         !$omp end parallel
         
         if (doing_R .or. doing_CR) then
            ! These terms are always present in the denominator
            D_T = D_T + 1 + 2*sum(t1**2) + sum(asym_t2*c_oovv)
         end if

         ! Get the name of the calculation right
         calcname = 'CCSD'
         if (sys%ccsd_t_paren) then
            calcname = trim(calcname)//'(T)'
         else
            calcname = trim(calcname)//'[T]'
         end if

         if (sys%ccsd_t_renorm) calcname = 'renormalised '//trim(calcname)
         if (sys%ccsd_t_comp_renorm) calcname = 'completely renormalised '//trim(calcname)

         if ((.not. sys%ccsd_t_renorm) .and. (.not. sys%ccsd_t_comp_renorm)) then
            ! No additional denominator
            write(iunit, '(1X, A, 1X, F15.9)') 'Restricted '//trim(calcname)//' correlation energy (Hartree):', e_T
            sys%e_ccsd_t = e_T
         else if (sys%ccsd_t_renorm) then
            write(iunit, '(1X, A, 1X, F15.9)') 'Restricted '//trim(calcname)//' correlation energy (Hartree):', e_T/D_T
            sys%e_ccsd_t = e_T/D_T
         else if (sys%ccsd_t_comp_renorm) then
            write(iunit, '(1X, A, 1X, F15.9)') 'Restricted '//trim(calcname)//' correlation energy (Hartree):', e_CR/D_T
            sys%e_ccsd_t = e_CR/D_T
         end if
         end associate

      end subroutine do_ccsd_t_spatial

      subroutine build_cr_ccsd_t_intermediates(sys, cc_int, cc_amp, int_store_cc)

         ! Builds the two intermediate arrays needed for the completely renormalised CCSD[T]/(T) approaches as 
         ! defined in P. Piecuch et al., Computer Physics Communications 149 (2002) 71–96, 
         ! https://doi.org/10.1016/S0010-4655(02)00598-2
         ! Specifically, Table 1 defines the I_ciab'' and I_jkia'' quantities and other quantities involved in their formation
         ! In:
         !     cc_amp: holds the converged CC amplitudes.
         !     sys: system information
         ! In/out:
         !     cc_int: holds CC intermediate tensors needed.
         !     int_store_cc: stores the intermediate tensors carried over to CCSD(T).

         use integrals, only: int_store_t, int_store_cc_t
         use error_handling, only: check_allocate

         type(cc_int_t), intent(inout) :: cc_int
         type(cc_amp_t), intent(in) :: cc_amp
         type(system_t), intent(in) :: sys
         type(int_store_cc_t), intent(inout) :: int_store_cc

         real(p), allocatable, dimension(:,:,:,:) :: x_ovov_p, x_voov_p, x_vvvo, x_ovoo, x_ovov_pp, x_voov_pp, x_vvvo_p

         integer :: i, j, k, a, b, c, e, m, ierr
         real(p) :: tmp

         deallocate(cc_int%I_vovv_p, cc_int%I_ooov_p, cc_int%I_oooo,cc_int%I_ovov,cc_int%I_voov, cc_int%I_vv, cc_int%I_oo, &
                    cc_int%I_oo_p, cc_int%c_oovv, cc_int%x_voov, cc_int%D_ia, cc_int%D_ijab)

         ! A bug (feature??) means that if we allocate these arrays within the associate block, we have to invoke the 
         ! full name (including the derived type) later on, otherwise a segfault occurs.
         allocate(int_store_cc%I_vovv_pp(sys%nvirt,sys%nocc,sys%nvirt,sys%nvirt), source=0.0_p, stat=ierr)
         call check_allocate('int_store_cc%I_vovv_pp',sys%nvirt**3*sys%nocc,ierr)
         allocate(int_store_cc%I_ooov_pp(sys%nocc,sys%nocc,sys%nocc,sys%nvirt), source=0.0_p, stat=ierr)
         call check_allocate('int_store_cc%I_ooov_pp',sys%nocc**3*sys%nvirt,ierr)

         associate(t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab, &
            nocc=>sys%nocc, nvirt=>sys%nvirt, nbasis=>sys%nbasis, &
            I_vovv_pp=>int_store_cc%I_vovv_pp, I_ooov_pp=>int_store_cc%I_ooov_pp, &
            I_vo=>cc_int%I_vo, v_oovo=>cc_int%v_oovo, asym_t2=>cc_int%asym_t2, &
            v_oovv=>cc_int%v_oovv, v_ovov=>cc_int%v_ovov, v_vvov=>cc_int%v_vvov, v_vvvv=>cc_int%v_vvvv, v_oooo=>cc_int%v_oooo)

         !I_vovv_pp(c,i,a,b) = v_vvov(b,a,i,c) + v_vvvv(c,e,a,b)*t1(i,e) - x_ovov_p(i,c,m,a)*t1(m,b) - t1(m,a)*x_voov_p(c,i,m,b) &
         !            - I_vo(c,m)*t2(m,i,a,b) + x_vvvo(c,e,a,m)*asym_t2(i,m,b,e) - x_vvvo(e,c,a,m)*t2(m,i,e,b) &
         !            - t2(i,m,e,a)*x_vvvo(e,c,b,m) + t2(m,n,b,a)*x_ovoo(i,c,m,n)
         !We need: x_ovov_p, x_voov_p, x_vvvo, x_ovoo
         !
         !I_ooov_pp(j,k,i,a) = v_oovo(k,j,a,i) - v_oooo(j,k,i,m)*t1(m,a) + x_ovov_pp(j,e,i,a)*t1(k,e) + t1(j,e)*x_voov_pp(e,k,i,a) &
         !                            + x_ovoo(j,e,i,m)*asym_t2(k,m,a,e) - x_ovoo(j,e,m,i)*t2(m,k,e,a) &
         !                            - t2(m,j,a,e)*x_ovoo(k,e,m,i) + t2(k,j,e,f)*x_vvvo(e,f,a,i)
         !We additionally need: x_ovov_pp, x_voov_pp
         !
         !x_ovov_p(j,b,i,a) = v_ovov(j,b,i,a) - 0.5*v_oovo(m,i,b,j)*t1(m,a) + t1(j,e)*x_vvvo_p(b,e,a,i)
         !Turns out we need x_vvvo_p in turn
         !
         !x_vvvo_p(b,c,a,i) = v_vvov(c,b,i,a) - 0.5*t1(m,a)*v_oovv(m,i,b,c)
         !
         !x_voov_p(b,j,i,a) = v_oovv(i,j,b,a) - 0.5*v_oovo(i,m,b,j)*t1(m,a) + x_vvvo_p(e,b,a,i)*t1(j,e)
         !
         !x_vvvo(b,c,a,i) = x_vvvo_p(b,c,a,i) - 0.5*t1(m,a)*v_oovv(m,i,b,c)
         !
         !x_ovoo(k,a,i,j) = v_oovo(j,i,a,k) + t1(k,e)*v_oovv(i,j,e,a)
         !
         !x_ovov_pp(j,b,i,a) = v_ovov(j,b,i,a) - v_oovo(m,i,b,j)*t1(m,a) + 0.5*t1(j,e)*x_vvvo(b,e,a,i)
         !
         !x_voov_pp(b,j,i,a) = v_oovv(i,j,b,a) - v_oovo(i,m,b,j)*t1(m,a) + 0.5*x_vvvo(e,b,a,i)*t1(j,e)

         allocate(x_ovov_p(nocc,nvirt,nocc,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('x_ovov_p',(nocc*nvirt)**2,ierr)
         allocate(x_voov_p(nvirt,nocc,nocc,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('x_voov_p',(nocc*nvirt)**2,ierr)
         allocate(x_vvvo(nvirt,nvirt,nvirt,nocc), source=0.0_p, stat=ierr)
         call check_allocate('x_vvvo',nvirt**3*nocc,ierr)
         allocate(x_ovoo(nocc,nvirt,nocc,nocc), source=0.0_p, stat=ierr)
         call check_allocate('x_ovoo',nocc**3*nvirt,ierr)
         allocate(x_ovov_pp(nocc,nvirt,nocc,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('x_ovov_pp',(nocc*nvirt)**2,ierr)
         allocate(x_voov_pp(nvirt,nocc,nocc,nvirt), source=0.0_p, stat=ierr)
         call check_allocate('x_voov_pp',(nocc*nvirt)**2,ierr)
         allocate(x_vvvo_p(nvirt,nvirt,nvirt,nocc), source=0.0_p, stat=ierr)
         call check_allocate('x_vvvo_p',nvirt**3*nocc,ierr)


         !$omp parallel default(none) &
         !$omp shared(sys, cc_amp, cc_int, int_store_cc, x_ovov_p, x_voov_p, x_vvvo, x_ovoo, x_ovov_pp, x_voov_pp, x_vvvo_p)

         !$omp do schedule(static, 10) collapse(2) 
         do i = 1, nocc
            do a = 1, nvirt
               do c = 1, nvirt
                  do b = 1, nvirt
                     x_vvvo_p(b,c,a,i) = v_vvov(c,b,i,a) - sum(0.5*t1(:,a)*v_oovv(:,i,b,c))
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(static, 10) collapse(2) 
         do a = 1, nvirt
            do i = 1, nocc
               do b = 1, nvirt
                  do j = 1, nocc
                     x_ovov_p(j,b,i,a) = v_ovov(j,b,i,a) - 0.5*sum(v_oovo(:,i,b,j)*t1(:,a)) + sum(t1(j,:)*x_vvvo_p(b,:,a,i))
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(static, 10) collapse(2) 
         do a = 1, nvirt
            do i = 1, nocc
               do j = 1, nocc
                  do b = 1, nvirt
                     x_voov_p(b,j,i,a) = v_oovv(i,j,b,a) - 0.5*sum(v_oovo(i,:,b,j)*t1(:,a)) + sum(x_vvvo_p(:,b,a,i)*t1(j,:))
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(static, 10) collapse(2) 
         do i = 1, nocc
            do a = 1, nvirt
               do c = 1, nvirt
                  do b = 1, nvirt
                     x_vvvo(b,c,a,i) = x_vvvo_p(b,c,a,i) - 0.5*sum(t1(:,a)*v_oovv(:,i,b,c))
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(static, 10) collapse(2) 
         do j = 1, nocc
            do i = 1, nocc
               do a = 1, nvirt
                  do k = 1, nocc
                     x_ovoo(k,a,i,j) = v_oovo(j,i,a,k) + sum(t1(k,:)*v_oovv(i,j,:,a))
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(static, 10) collapse(2) 
         do a = 1, nvirt
            do i = 1, nocc
               do b = 1, nvirt
                  do j = 1, nocc
                     x_ovov_pp(j,b,i,a) = v_ovov(j,b,i,a) - sum(v_oovo(:,i,b,j)*t1(:,a)) + 0.5*sum(t1(j,:)*x_vvvo(b,:,a,i))
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(static, 10) collapse(2) 
         do a = 1, nvirt
            do i = 1, nocc
               do j = 1, nocc
                  do b = 1, nvirt
                     x_voov_pp(b,j,i,a) = v_oovv(i,j,b,a) - sum(v_oovo(i,:,b,j)*t1(:,a)) + 0.5*sum(x_vvvo(:,b,a,i)*t1(j,:))
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(static, 10) collapse(3) private(tmp)
         do b = 1, nvirt
            do a = 1, nvirt
               do i = 1, nocc
                  do c = 1, nvirt
                     I_vovv_pp(c,i,a,b) = v_vvov(b,a,i,c) + sum(v_vvvv(:,c,b,a)*t1(i,:)) - sum(x_ovov_p(i,c,:,a)*t1(:,b)) &
                                    - sum(t1(:,a)*x_voov_p(c,i,:,b)) - sum(I_vo(c,:)*t2(:,i,a,b)) + sum(t2(:,:,b,a)*x_ovoo(i,c,:,:)) 
                     tmp = 0.0_p
                     do m = 1, nocc
                        tmp = tmp + sum(x_vvvo(c,:,a,m)*asym_t2(i,m,b,:)) &
                                             - sum(x_vvvo(:,c,a,m)*t2(m,i,:,b)) - sum(t2(m,i,a,:)*x_vvvo(:,c,b,m))
                     end do
                     I_vovv_pp(c,i,a,b) = I_vovv_pp(c,i,a,b) + tmp
                  end do
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(static, 10) collapse(3) private(tmp)
         do a = 1, nvirt
            do i = 1, nocc
               do k = 1, nocc
                  do j = 1, nocc
                     I_ooov_pp(j,k,i,a) = v_oovo(k,j,a,i) - sum(v_oooo(:,i,k,j)*t1(:,a)) + sum(x_ovov_pp(j,:,i,a)*t1(k,:)) &
                                          + sum(t1(j,:)*x_voov_pp(:,k,i,a)) + sum(t2(k,j,:,:)*x_vvvo(:,:,a,i)) 
                     tmp = 0.0_p
                     do e = 1, nocc
                        tmp = tmp + sum(x_ovoo(j,e,i,:)*asym_t2(:,k,e,a)) &
                                             - sum(x_ovoo(j,e,:,i)*t2(:,k,e,a)) - sum(t2(:,j,a,e)*x_ovoo(k,e,:,i)) 
                     end do
                     I_ooov_pp(j,k,i,a) = I_ooov_pp(j,k,i,a) + tmp
                  end do
               end do
            end do
         end do
         !$omp end do
         !$omp end parallel

         deallocate(x_ovov_p, x_voov_p, x_vvvo, x_ovoo, x_ovov_pp, x_voov_pp, x_vvvo_p)

         end associate

      end subroutine build_cr_ccsd_t_intermediates

end module ccsd
