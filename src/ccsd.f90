module ccsd

   use const
   use system, only: system_t

   implicit none

   type cc_amp_t
      real(p), allocatable :: t_ia(:,:)
      real(p), allocatable :: t_ijab(:,:,:,:)
   end type cc_amp_t

   type cc_int_t
      ! Intermediate quantities
      ! Curly F tensors (Eqs. 3-5)
      real(dp), allocatable :: F_oo(:,:), F_vv(:,:), F_ov(:,:)
      ! Curly W tensors (Eqs. 6,8), see Appendix for why 7 is avoided
      real(dp), allocatable :: W_oooo(:,:,:,:), W_ovvo(:,:,:,:)
      ! Effective two-particle excitation operators tau and tilde tau (Eqs. 9-10)
      real(dp), allocatable :: tau_tilde(:,:,:,:), tau(:,:,:,:)
      ! Energy denominators (Eqs. 12-13)
      real(dp), allocatable :: D_ia(:,:), D_ijab(:,:,:,:)
      real(dp), allocatable :: tmp_tia(:,:), tmp_tijab(:,:,:,:)
   end type cc_int_t
   
   contains

      subroutine do_ccsd(sys, int_store)
         ! This is the spinorbital formulation of J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett, J. Chem. Phys. volume 94, pp. 4334-4345 (1991) (https://doi.org/10.1063/1.460620)
         ! [TODO]: closed-shell (spatial) formulation of Hirata, S.; Podeszwa, R.; Tobita, M.; Bartlett, R. J. J. Chem. Phys. 2004, 120 (6), 2581 (https://doi.org/10.1063/1.1637577)

         use integrals, only: int_store_t, eri_ind
         use error_handling, only: check_allocate
         use system, only: state_t

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(inout) :: int_store
         integer :: i, j, a, b, p, q, r, s
         integer :: pr, qr, pa, pb, qa, qb, ra, rb, sa, sb
         real(dp) :: prqs, psqr

         type(state_t) :: st
         
         type(cc_amp_t) :: cc_amp
         integer :: ierr, iter
         integer, parameter :: iunit = 6
         integer, parameter :: maxiter = 50

         type(cc_int_t) :: cc_int

         write(iunit, '(1X, 10("-"))')
         write(iunit, '(1X, A)') 'CCSD'
         write(iunit, '(1X, 10("-"))')

         ! Form the antisymmetrised spinorbital basis ERIs: <pq||rs> = <pq|rs> - <pq|sr>
         ! where <pq|rs> = (pr|qs) * delta(sigma_p,sigma_r) * delta(sigma_q,sigma_s)
         write(iunit, '(1X, A)') 'Forming antisymmetrised spinorbital ERIs...'
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

         call init_cc(sys, int_store, cc_amp, cc_int)

         associate(nbasis=>sys%nbasis, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab, asym=>int_store%asym_spinorb)
         do iter = 1, maxiter
            ! Update intermediate tensors
            call update_cc_energy(sys, st, int_store, cc_amp)
            call build_tau(sys, cc_amp, cc_int)
            call build_F(sys, cc_int, cc_amp, asym)
            call build_W(sys, cc_int, cc_amp, asym)
            call update_amplitudes(sys, cc_amp, int_store, cc_int)
            print*, st%energy        
         end do
         end associate

      end subroutine do_ccsd

      subroutine init_cc(sys, int_store, cc_amp, cc_int)
         use integrals, only: int_store_t
         use error_handling, only: check_allocate 

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(in) :: int_store
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         integer :: p, q, r, s, i, j, a, b, ierr

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
         allocate(cc_int%W_ovvo(2*nocc,2*nocc+1:2*n,2*nocc+1:2*n,2*nocc), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%W_ovvo', 16*nocc**2*(n-nocc)**2, ierr)
         allocate(cc_int%tau(2*nocc,2*nocc,2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%tau', 16*nocc**2*(n-nocc)**2, ierr)
         allocate(cc_int%tau_tilde(2*nocc,2*nocc,2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_int%tau_tilde', 16*nocc**2*(n-nocc)**2, ierr)
         end associate

      end subroutine init_cc

      subroutine update_amplitudes(sys, cc_amp, int_store, cc_int)
         use integrals, only: int_store_t

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(in) :: int_store
         type(cc_amp_t), intent(inout) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         integer :: i, j, a, b, m, n, e, f
         real(p) :: W_abef
         
         associate(tmp_t1=>cc_int%tmp_tia,tmp_t2=>cc_int%tmp_tijab,t1=>cc_amp%t_ia,t2=>cc_amp%t_ijab,asym=>int_store%asym_spinorb, &
            F_oo=>cc_int%F_oo,F_ov=>cc_int%F_ov,F_vv=>cc_int%F_vv,D_ia=>cc_int%D_ia,D_ijab=>cc_int%D_ijab,W_oooo=>cc_int%W_oooo, &
            W_ovvo=>cc_int%W_ovvo,tau=>cc_int%tau,tau_tilde=>cc_int%tau_tilde,nocc=>sys%nocc,nbasis=>sys%nbasis)

         ! Update T1
         tmp_t1 = 0.0_dp
         do i = 1, 2*nocc
            do a = 2*nocc+1, 2*nbasis
               do m = 1, 2*nocc
                  tmp_t1(i,a) = tmp_t1(i,a) - t1(m,a)*F_oo(m,i)
                  do e = 2*nocc+1, 2*nbasis
                     tmp_t1(i,a) = tmp_t1(i,a) + t2(i,m,a,e)*F_ov(m,e)
                     do f = 2*nocc+1, 2*nbasis
                        tmp_t1(i,a) = tmp_t1(i,a) - 0.5*t2(i,m,e,f)*asym(m,a,e,f)
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
         tmp_t1 = tmp_t1/D_ia

         ! Update T2
         tmp_t2 = 0.0_dp
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
                           W_abef = asym(a,b,e,f)
                           do n = 1, 2*nocc
                              W_abef = W_abef - tmp_t1(n,b)*asym(a,n,e,f) + tmp_t1(n,a)*asym(b,n,e,f)
                           end do
                           tmp_t2(i,j,a,b) = tmp_t2(i,j,a,b) + 0.5*tau(i,j,e,f)*W_abef
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
                           + t2(i,m,a,e)*W_ovvo(m,b,e,j) - tmp_t1(i,e)*tmp_t1(m,a)*asym(m,b,e,j)&
                           - t2(j,m,a,e)*W_ovvo(m,b,e,i) + tmp_t1(j,e)*tmp_t1(m,a)*asym(m,b,e,i)&
                           - t2(i,m,b,e)*W_ovvo(m,a,e,j) + tmp_t1(i,e)*tmp_t1(m,b)*asym(m,a,e,j)&
                           + t2(j,m,b,e)*W_ovvo(m,a,e,i) - tmp_t1(j,e)*tmp_t1(m,b)*asym(m,a,e,i)
                        end do 
                     end do
                  end do
               end do
            end do
         end do
         t1 = tmp_t1
         t2 = tmp_t2/D_ijab
         end associate

      end subroutine update_amplitudes

      subroutine update_cc_energy(sys, st, int_store, cc_amp)
         use system, only: state_t
         use integrals, only: int_store_t

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(in) :: int_store
         type(cc_amp_t), intent(in) :: cc_amp
         type(state_t), intent(inout) :: st

         integer :: i, j, a, b

         st%energy_old = st%energy
         st%energy = 0.0_p
         associate(nocc=>sys%nocc, nbasis=>sys%nbasis, asym=>int_store%asym_spinorb, t1=>cc_amp%t_ia, t2=>cc_amp%t_ijab)
            do i = 1, 2*nocc
               do j = 1, 2*nocc
                  do a = 2*nocc+1, 2*nbasis
                     do b = 2*nocc+1, 2*nbasis
                        st%energy = st%energy + 0.25*asym(i,j,a,b)*(t2(i,j,a,b)+2*t1(i,a)*t1(j,b))
                     end do
                  end do
               end do
            end do
         end associate

      end subroutine update_cc_energy

      subroutine build_tau(sys, cc_amp, cc_int)
         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         type(cc_int_t), intent(inout) :: cc_int

         integer :: i, j, a, b
         real(p) :: ia, ja, x

         associate(n=>sys%nbasis, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab,&
            tau=>cc_int%tau, tau_tilde=>cc_int%tau_tilde)
            tau = 0.0_p
            tau_tilde = 0.0_p
            do i = 1, 2*nocc
               do a = 2*nocc+1, 2*n
                  ia = t_ia(i,a)
                  do j = 1, 2*nocc
                     ja = t_ia(j,a)
                     do b = 2*nocc+1, 2*n
                        x = ia*t_ia(j,b) - t_ia(i,b)*ja
                        tau_tilde(i,j,a,b) = t_ijab(i,j,a,b) + 0.5*x
                        tau(i,j,a,b) = tau_tilde(i,j,a,b) + 0.5*x
                     end do
                  end do
               end do
            end do
         end associate

      end subroutine build_tau

      subroutine build_F(sys, cc_int, cc_amp, asym)
         ! Eqs. 3-5, as we use HF reference all terms involving the fock matrix vanishes
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

         ! F_mi = \sum_{en} t_n^e <mn||ie> + 0.5 \sum_{nef} \tau~_{in}^{ef} <mn||ef>
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

         ! F_me = \sum_{nf} t_n^f * <mn||ef>
         do m = 1, 2*nocc
            do e = 2*nocc+1, 2*nbasis
               do n = 1, 2*nocc
                  do f = 2*nocc+1, 2*nbasis
                     F_ov(m,e) = F_ov(m,e) + t_ia(n,f)*asym(m,n,e,f)
                  end do
               end do
            end do
         end do
         end associate
      end subroutine build_F

      subroutine build_W(sys, cc_int, cc_amp, asym)
         ! Eqs. 6 and 8
         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         real(p), intent(in) :: asym(:,:,:,:)
         type(cc_int_t), intent(inout) :: cc_int
         integer :: m, n, i, j, e, f, b
         real(p) :: x

         associate(nbasis=>sys%nbasis, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab,&
            W_oooo=>cc_int%W_oooo, W_ovvo=>cc_int%W_ovvo, tau=>cc_int%tau)
         W_oooo = 0.0_p
         W_ovvo = 0.0_p
         do m = 1, 2*nocc
            do n = 1, 2*nocc
               do i = 1, 2*nocc
                  do j = 1, 2*nocc
                     W_oooo(m,n,i,j) = asym(m,n,i,j)
                     do e = 2*nocc+1, 2*nbasis
                        W_oooo(m,n,i,j) = W_oooo(m,n,i,j) + t_ia(j,e)*asym(m,n,i,e) - t_ia(i,e)*asym(m,n,j,e)
                        do f = 2*nocc+1, 2*nbasis
                           W_oooo(m,n,i,j) = W_oooo(m,n,i,j) + 0.5*t_ijab(i,j,e,f)*asym(m,n,e,f)
                        end do
                     end do
                  end do
               end do
            end do
         end do

         do m = 1, 2*nocc
            do b = 2*nocc+1, 2*nbasis
               do e = 2*nocc+1, 2*nbasis
                  do j = 1, 2*nocc
                     W_ovvo(m,b,e,j) = asym(m,b,e,j)
                     do f = 2*nocc+1, 2*nbasis
                        W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) + t_ia(j,f)*asym(m,b,e,f)
                     end do

                     do n = 1, 2*nocc
                        x = t_ia(n,b)
                        W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) - x*asym(m,n,e,j)
                        do f = 2*nocc+1, 2*nbasis
                           W_ovvo(m,b,e,j) = W_ovvo(m,b,e,j) - (0.5*t_ijab(j,n,f,b) + x*t_ia(j,f)) * asym(m,n,e,f)
                        end do
                     end do
                  end do
               end do
            end do
         end do

         end associate
      end subroutine build_W

end module ccsd