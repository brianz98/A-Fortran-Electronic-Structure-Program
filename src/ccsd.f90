module ccsd

   use const
   use system, only: system_t

   implicit none

   type cc_amp_t
      real(p), allocatable :: t_ia(:,:)
      real(p), allocatable :: t_ijab(:,:,:,:)
   end type cc_amp_t
   
   contains

      subroutine do_ccsd(sys, int_store)
         ! This is the spinorbital formulation of J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett, J. Chem. Phys. volume 94, pp. 4334-4345 (1991) (https://doi.org/10.1063/1.460620)
         ! [TODO]: closed-shell (spatial) formulation of Hirata, S.; Podeszwa, R.; Tobita, M.; Bartlett, R. J. J. Chem. Phys. 2004, 120 (6), 2581 (https://doi.org/10.1063/1.1637577)

         use integrals, only: int_store_t, eri_ind
         use error_handling, only: check_allocate

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(inout) :: int_store
         integer :: i, j, a, b, p, q, r, s, rasym, sasym
         integer :: pr, qr, pa, pb, qa, qb, ra, rb, sa, sb
         real(dp) :: prqs, psqr, e_mp2
         type(cc_amp_t) :: cc_amp
         integer :: ierr
         integer, parameter :: iunit = 6

         ! Intermediate quantities
         ! Curly F tensors (Eqs. 3-5)
         real(dp), allocatable :: F_oo(:,:), F_vv(:,:), F_ov(:,:)
         ! Curly W tensors (Eqs. 6,8), see Appendix for why 7 is avoided
         real(dp), allocatable :: W_oooo(:,:,:,:), W_ovvo(:,:,:,:)
         ! Effective two-particle excitation operators tau and tilde tau (Eqs. 9-10)
         real(dp), allocatable :: tau_til(:,:,:,:), tau(:,:,:,:)
         ! Energy denominators (Eqs. 12-13)
         real(dp), allocatable :: D_ia(:,:), D_ijab(:,:,:,:)

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

         ! Forming the energy denominator matrices
         associate(n=>sys%nbasis, nocc=>sys%nocc, e=>sys%canon_levels)
         allocate(D_ia(2*nocc, 2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('D_ia', 4*nocc*(n-nocc), ierr)
         allocate(D_ijab(2*nocc, 2*nocc, 2*nocc+1:2*n, 2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('D_ijab', 16*(nocc**2*(n-nocc)**2), ierr)
         do i = 1, nocc
            do a = nocc+1, n
               D_ia(2*i-1:2*i,2*a-1:2*a) = e(i)-e(a)
               do j = 1, nocc
                  do b = nocc+1, n
                     D_ijab(2*i-1:2*i, 2*j-1:2*j, 2*a-1:2*a, 2*b-1:2*b) = e(i)+e(j)-e(a)-e(b)
                  end do
               end do
            end do
         end do

         ! Initialise t_i^a=0 and t_{ij}^{ab}=MP1 WF
         allocate(cc_amp%t_ijab(2*nocc,2*nocc,2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_amp%t_ijab', 16*(nocc**2*(n-nocc)**2), ierr)
         allocate(cc_amp%t_ia(2*nocc, 2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('cc_amp%t_ia', 4*nocc*(n-nocc), ierr)
         end associate

         associate(doubles=>cc_amp%t_ijab, nocc=>sys%nocc, nvirt=>sys%nvirt, asym=>int_store%asym_spinorb)
         e_mp2 = 0.0_dp
         do p = 1, 2*nocc
            do q = 1, 2*nocc
               do r = 2*sys%nocc+1, 2*sys%nbasis
                  do s = 2*sys%nocc+1, 2*sys%nbasis
                     doubles(p,q,r,s) = asym(p,q,r,s)/D_ijab(p,q,r,s)
                     e_mp2 = e_mp2 + doubles(p,q,r,s)*asym(p,q,r,s)
                  end do
               end do
            end do
         end do
         print*, e_mp2/4
         end associate

         ! Allocate the intermediate quantities
         associate(n=>sys%nbasis, nocc=>sys%nocc)
         allocate(F_vv(2*nocc+1:2*n,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('F_vv', 4*(n-nocc)**2, ierr)
         allocate(F_oo(2*nocc,2*nocc), source=0.0_dp, stat=ierr)
         call check_allocate('F_oo', 4*nocc**2, ierr)
         allocate(F_ov(2*nocc,2*nocc+1:2*n), source=0.0_dp, stat=ierr)
         call check_allocate('F_ov', 4*nocc*(n-nocc), ierr)
         end associate

      end subroutine do_ccsd

      subroutine build_tau(sys, cc_amp, tau, tau_tilde)
         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         real(p), intent(out) :: tau(:,:,:,:), tau_tilde(:,:,:,:)

         integer :: i, j, a, b
         real(p) :: ia, ja, x

         associate(n=>sys%nbasis, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab)
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

      subroutine build_F(sys, F_vv, F_oo, F_ov, cc_amp, asym, tau_tilde)
         ! Eqs. 3-5, as we use HF reference all terms involving the fock matrix vanishes
         type(system_t), intent(in) :: sys
         type(cc_amp_t), intent(in) :: cc_amp
         real(p), intent(in) :: asym(:,:,:,:)
         real(p), intent(in) :: tau_tilde(:,:,:,:)
         real(p), intent(out) :: F_vv(:,:), F_oo(:,:), F_ov(:,:)

         integer :: a, e, f, m, n, i

         associate(nbasis=>sys%nbasis, nocc=>sys%nocc, t_ia=>cc_amp%t_ia, t_ijab=>cc_amp%t_ijab)
         ! F_ae = \sum_{mf} t_m^f * <ma||fe> - 0.5 \sum_{mnf} \tau~_{mn}^{af} <mn||ef>
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

end module ccsd