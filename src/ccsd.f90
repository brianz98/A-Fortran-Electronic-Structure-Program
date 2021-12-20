module ccsd

   use const

   implicit none

   contains

      subroutine do_ccsd(sys, int_store)
         use system, only: system_t
         use integrals, only: int_store_t, eri_ind

         type(system_t), intent(in) :: sys
         type(int_store_t), intent(inout) :: int_store
         integer :: p, q, r, s
         integer :: pr, qr, pa, pb, qa, qb, ra, rb, sa, sb
         real(dp) :: a, b
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

         associate(n=>sys%nbasis, asym=>int_store%asym_spinorb, eri=>int_store%eri_mo)
         ! [TODO]: this isn't necessary but easier for debugging, change to spin lookup arrays or something similar
         allocate(int_store%asym_spinorb(n*2,n*2,n*2,n*2), source=0.0_dp)
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

      end subroutine do_ccsd


end module ccsd