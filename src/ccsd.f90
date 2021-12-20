module ccsd

   use const

   implicit none

   type cc_amp_t
      real(p), allocatable :: diag(:)
      real(p), allocatable :: t_ia(:,:)
      real(p), allocatable :: t_ijab(:,:,:,:)
   end type cc_amp_t
   
   contains

      subroutine do_ccsd(sys, int_store)
         use system, only: system_t
         use integrals, only: int_store_t, eri_ind

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(inout) :: int_store
         integer :: i, p, q, r, s
         integer :: pr, qr, pa, pb, qa, qb, ra, rb, sa, sb
         real(dp) :: a, b, e_mp2
         type(cc_amp_t) :: cc_amp
         integer, parameter :: iunit = 6

         write(iunit, '(1X, 10("-"))')
         write(iunit, '(1X, A)') 'CCSD'
         write(iunit, '(1X, 10("-"))')

         write(iunit, '(1X, A)') 'Performing spatial to spinorbital MO ERI transformation...'

         !if (.false.) then
         ! Permutational symmetry hurts my brain, let's get a naive version working first
         !associate(n=>sys%nbasis, asym=>int_store%asym_spinorb, eri=>int_store%eri_mo)
         ! [TODO]: this isn't necessary but easier for debugging, change to spin lookup arrays or something similar
         !allocate(asym(n*2,n*2,n*2,n*2), source=0.0_dp)
         !do p = 1, n
         !   pa = 2*p-1
         !   do q = 1, p
         !      qa = 2*q-1
         !      do r = 1, p
         !         ra = 2*r-1
         !         pr = eri_ind(p,r)
         !         qr = eri_ind(q,r)
         !         if (r==p) then
         !            s_up = q
         !         else
         !            s_up = r
         !         end if
         !         do s = 1, s_up
         !            sa = 2*s-1
         !            prqs = eri(eri_ind(pr, eri_ind(q,s)))
         !            psqr = eri(eri_ind(eri_ind(p,s), qr))
         !            
         !            ! <aa|aa> nonzero: <pq||rs> = -<pq||sr> = <rs||pq> = -<sr||pq>
         !            asym(pa,qa,ra,sa) = prqs - psqr
         !            asym(pa,qa,sa,ra) = psqr - prqs
         !            asym(ra,sa,pa,qa) = prqs - psqr
         !            asym(sa,ra,pa,qa) = psqr - prqs

                     ! <bb|bb> nonzero
         !            asym(pa+1,qa+1,ra+1,sa+1) = prqs - psqr
         !            asym(pa+1,qa+1,sa+1,ra+1) = psqr - prqs
         !            asym(ra+1,sa+1,pa+1,qa+1) = prqs - psqr
         !            asym(sa+1,ra+1,pa+1,qa+1) = psqr - prqs

                     ! <ab|ab> nonzero: <pq||rs> = <pq|rs>, <pq||sr> = -<pq|rs> etc
         !            asym(pa,qa+1,ra,sa+1) = prqs
         !            asym(pa,qa+1,sa+1,ra) = -prqs
         !            asym(ra,sa+1,pa,qa+1) = prqs
         !            asym(sa+1,ra,pa,qa+1) = -prqs

                     ! <ba|ba> nonzero
         !            asym(pa+1,qa,ra+1,sa) = prqs
         !            asym(pa+1,qa,sa,ra+1) = -prqs
         !            asym(ra+1,sa,pa+1,qa) = prqs
         !            asym(sa,ra+1,pa+1,qa) = -prqs
         !         end do
         !      end do
         !   end do
         !end do
         !end associate
         !end if

         allocate(int_store%asym_spinorb(sys%nbasis*2,sys%nbasis*2,sys%nbasis*2,sys%nbasis*2), source=0.0_dp)
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
                     a = eri(eri_ind(pr,eri_ind(q,s)))
                     b = eri(eri_ind(eri_ind(p,s),qr))
                     ! Draw a spin decision tree to visualise
                     asym(pa,qa,ra,sa) = a-b
                     asym(pb,qb,rb,sb) = a-b
                     asym(pa,qb,ra,sb) = a
                     asym(pb,qa,rb,sa) = a
                     asym(pa,qb,rb,sa) = -b
                     asym(pb,qa,ra,sb) = -b
                  end do
               end do
            end do
         end do
         end associate

         allocate(cc_amp%diag(2*sys%nbasis))
         do i = 1, sys%nbasis
            cc_amp%diag(2*i-1:2*i) = sys%canon_levels(i)
         end do

         allocate(cc_amp%t_ijab(2*sys%nbasis,2*sys%nbasis,2*sys%nbasis,2*sys%nbasis), source=0.0_dp)
         associate(doubles=>cc_amp%t_ijab, nocc=>sys%nocc, n=>sys%nbasis, asym=>int_store%asym_spinorb, e=>cc_amp%diag)
         e_mp2 = 0.0_dp
         do p = 1, 2*nocc
            do q = 1, 2*nocc
               do r = 2*nocc+1, 2*n
                  do s = 2*nocc+1, 2*n
                     doubles(p,q,r,s) = asym(p,q,r,s)/(e(p)+e(q)-e(r)-e(s))
                     e_mp2 = e_mp2 + doubles(p,q,r,s)*asym(p,q,r,s)
                  end do
               end do
            end do
         end do
         print*, e_mp2/4
         end associate



      end subroutine do_ccsd


end module ccsd