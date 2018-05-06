


subroutine reference_schroe( r, l )
   use common_variables
   implicit real(8) (a-h,o-z)
   integer, intent(in) :: l
   real(8), intent(in) :: r(1:msh)
   integer, parameter :: jjj = 10
   real(8), parameter :: alp = 1.d0/274.074d0
   real(8) :: dvdx(1:msh)
   real(8) :: abm(1:msh), a(1:msh), b(1:msh)
   real(8) :: aa(1:jjj,1:jjj), bb(1:jjj), aaa(1:jjj)
   real(8) :: y(1:msh), xx(1:msh), f(1:msh)
   character(13) :: ofile

   cll1 = dble( l*(l+1) )

   call fd1 ( 1, icl+60, mesh, v, dx, dvdx )

   do nrep = 1, nref, 1
      e = eref(nrep)

      do i = 1, icl+60
         ri2 = r(i)*r(i)
         abm(i) = 1.d0 + alp*alp*(e - v(i))
           a(i) = alp*alp/abm(i)*dvdx(i)/r(i)
           b(i) = -a(i)/r(i) - cll1/ri2 - abm(i)*(v(i) - e)
      end do
      do i = 1, 5
         do j = 1, 5
            ri2 = r(i)*r(i)
            aa(i,j) = ((l+j+1)*((l+j) + a(i)*r(i))+b(i)*ri2) * r(i)**(l+j-1)
         end do
         bb(i) = - ((l+1)*(l + a(i)*r(i)) + b(i)*ri2) * r(i)**(l-1)
      end do

      call gauss( aa, bb, aaa, 5, jjj )

      do i = 1, 5
          y(i) = r(i)**(l+0.5)
         xx(i) = (l+0.5)*r(i)**(l+0.5)
          f(i) = 0
         do j = 1, 5
             y(i) = y(i) + aaa(j)*r(i)**(l+j+0.5)
            xx(i) = xx(i) + (l+j+0.5)*aaa(j)*r(i)**(l+j+0.5)
         enddo
         ri2 = r(i)*r(i)
         f(i) = -a(i)*r(i)*(xx(i) - y(i)*0.5) + ((l+0.5d0)**2.d0 &
         &    + abm(i)*(v(i) - e)*ri2)*y(i)
      end do

      !--- Adams-Bashforth Adams-Moulton -------------------------------------

      do i = 5, icl + 50
         y(i+1)=y(i)+dx/720*(1901*xx(i)-2774*xx(i-1)+2616*xx(i-2) &
         -1274*xx(i-3)+251*xx(i-4))
         xx(i+1)=xx(i)+dx/720*(1901*f(i)-2774*f(i-1)+2616*f(i-2) &
         -1274*f(i-3)+251*f(i-4))
         f(i+1)=-a(i+1)*r(i+1)*(xx(i+1)-y(i+1)*0.5)+((0.5+l)**2.d0 &
         +abm(i+1)*(v(i+1)-e)*r(i+1)*r(i+1))*y(i+1)

         y(i+1)=y(i)+dx/720*(251*xx(i+1)+646*xx(i)-264*xx(i-1) &
         +106*xx(i-2)-19*xx(i-3))
         xx(i+1)=xx(i)+dx/720*(251*f(i+1)+646*f(i)-264*f(i-1) &
         +106*f(i-2)-19*f(i-3))
         f(i+1)=-a(i+1)*r(i+1)*(xx(i+1)-y(i+1)*0.5)+((0.5+l)**2.d0 &
         +abm(i+1)*(v(i+1)-e)*r(i+1)*r(i+1))*y(i+1)
      end do

      do i = 1, icl+50
!          pae0(i) = sqrt(r(i))*y(i)
         plj(nrep,i) = sqrt(r(i))*y(i)
         pae(nrep,i) = sqrt(r(i))*y(i)
      end do
   end do

   ofile = (ms(l+1)//'_Pae.dat')
   open(9,file=ofile)
   do i = 1 , icl + 50
      write(9,*) r(i), ( plj(j,i), j = 1, nref )
   end do
   write(*,'(4x,2a)') 'Create : ', ofile
   close(9)

   return
end subroutine reference_schroe



subroutine rrkj2( r, l )
!---- Ultrasoft pseudopotential ----------------------------------------
!  input :
!         rcl1  ...... cutoff radius for Pus
!         rcl2  ...... cutoff radius for Pnc
!         nl    ...... reference number
!         elj   ...... reference energy
!         v     ...... all electron potential
!  out put :
!         pae   ...... the radial AE wave functions
!         pus   ...... the US pseudo-wave-functions
!         d2p   ...... the second derivative of Pus
!-----------------------------------------------------------------------
   use common_variables
   implicit real(8) (a-h,o-z)
   integer, parameter :: jjj = 10
   real(8) :: r(1:msh)
   real(8) :: b(1:msh)
   real(8) :: aa(1:jjj,1:jjj)
   real(8) :: fx(1:msh)
   real(8) :: pae0(1:msh), pae1(1:msh), pae2(1:msh), pae3(1:msh), pae4(1:msh)
   real(8) :: pus0(1:msh), pus1(1:msh), pus2(1:msh), pus3(1:msh), pus4(1:msh)
   real(8) :: d2f(1:msh), d2ft(1:msh)
   real(8) :: x0(0:10)                          ! zero point of bess
   real(8) :: f(1:msh), q(1:msh)
   real(8) :: ali(jjj), c2n(jjj), d2n(jjj), c2n2(jjj), d2n2(jjj)
   character(13) :: ofile
!-----------------------------------------------------------------------

   write(*,'(/, 25a)') ' ---- RRKJ2: Construction of ',  &
   &                   'Ultrasoft pseudo-wave-functions ',  &
   &                   ( '-', i = 1, 17 )

   x    = exp(dx)
   cr   = rmax/x**mesh
   rcl  = max(rcl1,rcl2)
   iicl = log(rcl/cr)/dx
   rcl  = r(iicl)

   ius  = log(rcl1/cr)/dx
   rcl1 = r(ius)

   write(*,'(4x,a,i4,a,f11.8)') 'rcut(',ius,') =', r(ius)

   ! ----- search zero point of Sphere bessel function
   delx = acos(-1.d0)/100.d0

   ifrag = 1
   x     = 0.d0
   x0(0) = 0.d0
   eps   = 0.1d-10

   do while ( ifrag < 5 )
      x = x + delx
      if( bess(l,x)*bess(l,x+delx) <= 0.d0) then
          x1 = x
          x2 = x + delx
          do while ( abs(bess(l,x)) > eps )
             x = x1 - bess(l,x1)/(bess(l,x2) - bess(l,x1))*(x2 - x1)
             if( bess(l,x1)*bess(l,x) < 0.d0 ) then
                x2 = x
             else if( bess(l,x2)*bess(l,x) < 0.d0 ) then
                x1 = x
             end if
          end do
          x0(ifrag) = x
          ifrag = ifrag + 1
      end if
   end do
   write(*,'(4x,a,i1,a)') 'zero point of sphere bessel function (l=',l,')'
   write(*,'(6x,a,5f10.6, /)') 'x0 = ', x0(0:4)

   ! ----- determine c_i & d_i  (i = 0, 2, ..., 8)
   qc    = acos(-1.d0)/rcl1
   gamm  = 10.d0/rcl1
   rtil  = 0.1d0*rcl1
   theta = qc*rtil

   sinr  = sin(theta)
   cosr  = cos(theta)
   sinr2 = sinr *sinr
   sinr3 = sinr2*sinr
   sinr4 = sinr3*sinr
   cosr2 = cosr *cosr
   cosr3 = cosr2*cosr
   cosr4 = cosr3*cosr

   qc2 = qc *qc
   qc3 = qc2*qc
   qc4 = qc3*qc

   g0r = exp( gamm*(rtil - rcl1) )
   g1r = gamm * g0r
   g2r = gamm * g1r
   g3r = gamm * g2r
   g4r = gamm * g3r

   e0r = sinr3
   e1r = 3.d0*qc  * sinr2*cosr
   e2r = 3.d0*qc2 * (2.d0*sinr*cosr2 - sinr3)
   e3r = 3.d0*qc3 * (2.d0*cosr3 -  7.d0*sinr2*cosr)
   e4r = 3.d0*qc4 * (7.d0*sinr3 - 20.d0*cosr2*sinr)

   f0r = sinr4
   f1r = 4.d0*qc  * sinr3*cosr
   f2r = 4.d0*qc2 * (3.d0*sinr2*cosr2 - sinr4)
   f3r = 4.d0*qc3 * (6.d0*sinr*cosr3 - 10.d0*sinr3*cosr)
   f4r = 4.d0*qc4 * (6.d0*cosr4 - 48.d0*sinr2*cosr2 + 10.d0*sinr4)

   do j = 1, 5
      aa(1,j) = rtil**( l + 2*j-1 )
   end do
   do i = 2, 5
      do j = 1, 5
         aa(i,j) = dble( l + 2*j-1-i+2 )*aa(i-1,j)/rtil
      end do
   end do

   b(1) = e0r*g0r
   b(2) = e1r*g0r + e0r*g1r
   b(3) = e2r*g0r + 2.d0*e1r*g1r + e0r*g2r
   b(4) = e3r*g0r + 3.d0*e2r*g1r + 3.d0*e1r*g2r + e0r*g3r
   b(5) = e4r*g0r + 4.d0*e3r*g1r + 6.d0*e2r*g2r + 4.d0*e1r*g3r + e0r*g4r

   call gauss( aa, b, c2n, 5, jjj )

   do j = 1, 5
      aa(1,j) = rtil**( l + 2*j-1 )
   end do
   do i = 2, 5
      do j = 1, 5
         aa(i,j) = dble( l + 2*j-1-i+2 )*aa(i-1,j)/rtil
      end do
   end do

   b(1) = f0r*g0r
   b(2) = f1r*g0r + f0r*g1r
   b(3) = f2r*g0r + 2.d0*f1r*g1r + f0r*g2r
   b(4) = f3r*g0r + 3.d0*f2r*g1r + 3.d0*f1r*g2r + f0r*g3r
   b(5) = f4r*g0r + 4.d0*f3r*g1r + 6.d0*f2r*g2r + 4.d0*f1r*g3r + f0r*g4r

   call gauss( aa, b, d2n, 5, jjj )

   ! ----- second derivative of F_lj and F^tilda_lj
   do i = 0, 4, 1
      c2n2(i+1) = dble( (l+2*i+1)*(l+2*i) )*c2n(i+1)
      d2n2(i+1) = dble( (l+2*i+1)*(l+2*i) )*d2n(i+1)
   end do
   do i = 1, icl
      if(      r(i) <= rtil) then
          ri2 = r(i)*r(i)
          d2f(i)  = r(i)**(l-1) * ( c2n2(1) + ( c2n2(2) + ( c2n2(3)  &
          &                     + ( c2n2(4) + c2n2(5)*ri2 )*ri2 )*ri2 )*ri2 )
          d2ft(i) = r(i)**(l-1) * ( d2n2(1) + ( d2n2(2) + ( d2n2(3)  &
          &                     + ( d2n2(4) + d2n2(5)*ri2 )*ri2 )*ri2 )*ri2 )
      else if( r(i) > rtil) then
          theta = qc*r(i)

          sinr  = sin(theta)
          cosr  = cos(theta)
          sinr2 = sinr *sinr
          sinr3 = sinr2*sinr
          sinr4 = sinr3*sinr
          cosr2 = cosr *cosr

          g0r = exp( gamm*(r(i) - rcl1) )
          g1r = gamm * g0r
          g2r = gamm * g1r

          e0r = sinr3
          e1r = 3.d0*qc  * sinr2*cosr
          e2r = 3.d0*qc2 * (2.d0*sinr*cosr2 - sinr3)

          f0r = sinr4
          f1r = 4.d0*qc  * sinr3*cosr
          f2r = 4.d0*qc2 * (3.d0*sinr2*cosr2 - sinr4)

          d2f(i)  = e2r*g0r + 2.d0*e1r*g1r + e0r*g2r
          d2ft(i) = f2r*g0r + 2.d0*f1r*g1r + f0r*g2r
      end if
   end do

   do nrep = 1, nref
      write(*,'(4x,a,f15.10)') 'Reference energy :', eref(nrep)

      do i = 1, icl + 50
         pae0(i) = plj(nrep,i)
      end do

   !--- derivative of P,  finite difference  -----------------------
      call fd1 ( 1, icl+10, mesh, pae0, dx, pae1 )
      call fd2 ( 1, icl+10, mesh, pae0, dx, pae2 )
      call fd3 ( 1, icl+10, mesh, pae0, dx, pae3 )
      call fd4 ( 1, icl+10, mesh, pae0, dx, pae4 )
   !-----------------------------------------------------------------------

      do j = 0, 1
         fx(1) = x0(j) + delx
         fx(2) = fx(1) + delx

         f(1) = 1.d0
         f(2) = 1.d0
         f(3) = 3.d0

         do while( f(1)*f(2) > 0.d0 )
            do i = 1, 2
               f(i) = - dble(l) + fx(i)*bess(l-1,fx(i))/bess(l,fx(i)) &
               &      - 1.d0/pae0(ius)*pae1(ius)
            end do
            fx(1) = fx(1) + delx
            fx(2) = fx(1) + delx
            if( fx(1) == x0(j+1) ) then
                if(      j == 0 ) then
                    x0(j)   = x0(j+1)
                    x0(j+1) = x0(j+2)
                    write(*,'(6x,a)') 'x0(0) <=== x0(1), x0(1) <=== x0(2)'
                else if( j == 1 ) then
                    x0(j)   = x0(j+1)
                    write(*,'(6x,a)') 'x0(1) <=== x0(2)'
               end if
               fx(1) = x0(j) + delx
               fx(2) = fx(1) + delx
            end if
         end do

         fx(1) = fx(1) - delx
         fx(2) = fx(2) - delx

         iter = 0
         do while( abs(f(3)) > eps )
            do i = 1, 2
               f(i) = - dble(l) + fx(i)*bess(l-1,fx(i))/bess(l,fx(i)) &
               &      - 1.d0/pae0(ius)*pae1(ius)
            end do

            fx(3) = fx(1) - f(1)/(f(2) - f(1))*(fx(2) - fx(1))
             f(3) = - dble(l) + fx(3)*bess(l-1,fx(3))/bess(l,fx(3)) &
             &      - 1.d0/pae0(ius)*pae1(ius)

            if(       f(1)*f(3) < 0.d0) then
                fx(2) = fx(3)
            else if ( f(2)*f(3) < 0.d0) then
                fx(1) = fx(3)
            else
               stop 'error2'
            end if
         !       write(*,*)'x0:',x0(0),x0(1)
         !       write(*,*)'fx:',fx(1),fx(2)
         !       write(*,*)'f:',f(1),f(2),f(3)
         !       write(*,*)fx(i)*bess(l-1,fx(i))/bess(l,fx(i))
         !     &               ,1.d0/pae0(ius)*pae1(ius)
         !       write(11,*)fx(i)*bess(l-1,fx(i))/bess(l,fx(i))

         !       write(*,*)bess(l-1,fx(i)),bess(l,fx(i))
            if ( iter>=2000000 .or. abs(fx(1)-fx(2))<1.d-5 ) then
               write(*,'(/, a)') ' ***** error : Did not decide qi.'
               write(*,'(a,2f18.7)')   ' x : ', fx(1:2)
               write(*,'(a,3es18.10)') ' f : ',  f(1:3)
               stop
            end if
            iter = iter + 1
         end do

         q(j+1) = fx(3)/rcl1
      end do
      write(*,'(6x,a,6f10.6, /)') 'qi = ', q(1:2)

   !------------------------------------------------------------------------
   !     Pus,lj(r)
   !------------------------------------------------------------------------
      rcl1_2 = rcl1  *rcl1
      rcl1_3 = rcl1_2*rcl1
      rcl1_4 = rcl1_3*rcl1
      ! ----- condition for first and second derivative of P(r_c)
      b(1) = pae1(ius)/rcl1
      b(2) = (pae2(ius) - pae1(ius))/rcl1_2

      do i = 1, 2
         t = q(i)*rcl1
         aa(1,i) = bess(l,t) + t*dbess(l,t,1)
         aa(2,i) = q(i)*(2.d0*dbess(l,t,1) + t*dbess(l,t,2))
      end do

      call gauss( aa, b, ali, 2, jjj )

!       open(9,file='Pus.dat')
!       do i = 1, ius
!          t = 0.d0
!          do j = 1, 2
!             t = t + ali(j)*r(i)*bess(l,q(j)*r(i))
!          end do
!          write(9,*) r(i), t, pae0(i)
!       end do
!       close(9)

      ! ----- condition for third derivative of P(r_c)
      x = (pae3(ius) - 3.d0*pae2(ius) + 2.d0*pae1(ius))/rcl1_3
      ay = 0.d0
      do i = 1, 2
         t = q(i)*rcl1
         qi2 = q(i)*q(i)
         ay = ay + ali(i)*qi2 * (3.d0*dbess(l,t,2) + t*dbess(l,t,3))
      end do
      ali(3) = (ay - x)/(6.d0*qc3)

      ! ----- condition for fourth derivative of P(r_c)
      x = (pae4(ius) - 6.d0*pae3(ius) + 11.d0*pae2(ius) - 6.d0*pae1(ius))/rcl1_4
      ay = 0.d0
      do i = 1, 2
         t = q(i)*rcl1
         qi3 = q(i)*q(i)
         qi3 = qi3 *q(i)
         ay = ay + ali(i)*qi3 * (4.d0*dbess(l,t,3) + t*dbess(l,t,4))
      end do
      ali(4) = (x - ay + 24.d0*qc3*gamm*ali(3))/(24.d0*qc4)

      write(*,'(6x,a,4es14.6)') 'alpha =', ali(1:4)

      do i = 1, icl + 50
         if(      i < ius) then
             pus0(i) = 0.d0
             do j = 1, 2
                t = q(j)*r(i)
                pus0(i) = pus0(i) + ali(j)*r(i)*bess(l,t)
             end do

             if(     r(i) < rtil ) then
                ri2 = r(i)*r(i)
                pus0(i) = pus0(i)  &
                &       + ali(3) * r(i)**(l+1)*( c2n(1) + ( c2n(2) + ( c2n(3)  &
                &       + ( c2n(4) + c2n(5)*ri2 )*ri2 )*ri2 )*ri2 )
                pus0(i) = pus0(i)  &
                &       + ali(4) * r(i)**(l+1)*( d2n(1) + ( d2n(2) + ( d2n(3)  &
                &       + ( d2n(4) + d2n(5)*ri2 )*ri2 )*ri2 )*ri2 )
            else if( r(i) >= rtil ) then
                theta = qc*r(i)

                sinr  = sin(theta)
                sinr2 = sinr *sinr
                sinr3 = sinr2*sinr
                sinr4 = sinr3*sinr

                g0r = exp( - gamm*(rcl1 - r(i)) )
                e0r = sinr3
                f0r = sinr4

                pus0(i) = pus0(i) + ali(3)*e0r*g0r + ali(4)*f0r*g0r
            end if
         else if( i >= ius ) then
            pus0(i) = pae0(i)
         end if
      end do

      call fd1( 1, ius+6, msh, pus0, dx, pus1 )
      call fd2( 1, ius+6, msh, pus0, dx, pus2 )
      call fd3( 1, ius+6, msh, pus0, dx, pus3 )
      call fd4( 1, ius+6, msh, pus0, dx, pus4 )

      write(*,'(22x,a)') 'P_US              P_AE              diff.'
      write(*,'(6x,a,3es18.10)')  &
      &     '1st derivs. : ', pus1(ius), pae1(ius), abs(pus1(ius) - pae1(ius))
      write(*,'(6x,a,3es18.10)')  &
      &     '2nd derivs. : ', pus2(ius), pae2(ius), abs(pus2(ius) - pae2(ius))
      write(*,'(6x,a,3es18.10)')  &
      &     '3rd derivs. : ', pus3(ius), pae3(ius), abs(pus3(ius) - pae3(ius))
      write(*,'(6x,a,3es18.10, /)')  &
      &     '4th derivs. : ', pus4(ius), pae4(ius), abs(pus4(ius) - pae4(ius))

      ! ----- second derivative of Pus
      do i = 1, icl
         d2p(nrep,i) = 0.d0
         do j = 1, 2
            t = q(j)*r(i)
            dq2 = 2.d0*dbess(l,t,1) + t*dbess(l,t,2)
            d2p(nrep,i) = d2p(nrep,i) + ali(j)*q(j)*dq2
         end do

         d2p(nrep,i) = d2p(nrep,i) + ali(3)*d2f(i) + ali(4)*d2ft(i)
      end do

      do i = 1, icl + 50
         pae(nrep,i) = pae0(i)
         pus(nrep,i) = pus0(i)
      end do
   end do

   ofile = (ms(l+1)//'_Pus.dat')
   open(9,file=ofile)
   do i = 1 , icl + 50
      write(9,*) r(i), ( pus(j,i), j = 1, nref )
   end do
   write(*,'(4x,2a, /)') 'Create : ', ofile
   close(9)

   return
end subroutine rrkj2



   SUBROUTINE rrkj3( r, l )
!---- Ultrasoft pseudopotential ----------------------------------------
!  input :
!         rcl1  ...... cutoff radius for Pus
!         rcl2  ...... cutoff radius for Pnc
!         nl    ...... reference number
!         elj   ...... reference energy
!         v     ...... all electron potential
!         plj   ...... the radial AE wave functions

!  out put :
!         pnc   ...... the NC pseudo-wave-functions

!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   ! parameter     ( iii = 5  )
   parameter     ( jjj  = 10 )
   ! parameter     ( alp  = 1/274.074  )
   parameter     ( irrkj = 3 )
   dimension     r(msh)
   ! dimension     a(msh),b(msh),aa(jjj,jjj),sl(jjj,jjj)
   dimension     b(msh),aa(jjj,jjj),sl(jjj,jjj)
   dimension     y(msh),fx(msh)
   ! dimension     bb(jjj)!,aaa(jjj)!,xx(msh)!,pae5(msh),pae6(msh)
   dimension     pae0(msh),pae1(msh),pae2(msh),pae3(msh),pae4(msh)
   dimension     pnc0(msh),pnc1(msh),pnc2(msh),pnc3(msh),pnc4(msh)
   dimension     x0(0:10),f(msh),q(msh),flj(msh),ftil(msh)
   dimension     ali(jjj),ali3(jjj),c2n(jjj),d2n(jjj)
   CHARACTER     OFILE*13
!-----------------------------------------------------------------------

   write(*,'(/, 25a)') ' ---- RRKJ3: Construction of ',  &
   &                   'Normconserving pseudo-wave-functions ',  &
   &                   ( '-', i = 1, 12 )

   x=exp(dx)
   cr=rmax/x**mesh
   rcl=max(rcl1,rcl2)
   iicl=log(rcl/cr)/dx
   rcl=r(iicl)

   inc=log(rcl2/cr)/dx
   rcl2=r(inc)
   write(*,'(4x,a,i4,a,f11.8)') 'rcut(',inc,') =', r(inc)

   qc=acos(-1.d0)/rcl2
   gamm=10.d0/rcl2
   rtil=0.1d0*rcl2
   the=qc*rtil

!---- c2n d2n ----------------------------------------------------------
   s0=sin(the)
   c0=cos(the)
   s02=s0*s0
   s03=s0*s02
   s04=s0*s03
   c02=c0*c0
   c03=c0*c02
   c04=c0*c03
   qc2=qc*qc
   qc3=qc*qc2
   qc4=qc*qc3
   s1=3*qc*c0*s02
   s2=-3*qc2*s03+6*qc2*c02*s0
   s3=-21*qc3*s02*c0+6*qc3*c03
   s4=-60*qc4*s0*c02+21*qc4*s03
   g0=exp(-gamm*(rcl2-rtil))
   g1=gamm*g0
   g2=gamm*g1
   g3=gamm*g2
   g4=gamm*g3
   b(1)=s03*g0

   b(2)=s1*g0+s03*g1

   b(3)=s2*g0+2*s1*g1+s03*g2

   b(4)=s3*g0+3*s2*g1+3*s1*g2+s03*g3

   b(5)=s4*g0+4*s3*g1+6*s2*g2+4*s1*g3+s03*g4


   do j=0,4
      aa(1,j+1)=rtil**(l+1+2*j)
      aa(2,j+1)=(l+1+2*j)*rtil**(l+2*j)
      aa(3,j+1)=(l+1+2*j)*(l+2*j)*rtil**(l+2*j-1)
      aa(4,j+1)=(l+1+2*j)*(l+2*j)*(l+2*j-1)*rtil**(l+2*j-2)
      aa(5,j+1)=(l+1+2*j)*(l+2*j)*(l+2*j-1)*(l+2*j-2)*rtil**(l+2*j-3)
   end do

!-----------------------------------------------------------------------
   call gauss(aa,b,c2n,5,jjj)
!-----------------------------------------------------------------------

   s1=qc*4*c0*s03
   s2=qc2*(-4*s04+12*c02*s02)
   s3=qc3*(-40*c0*s03+24*c03*s0)
   s4=qc4*(40*s04-192*c02*s02+24*c04)

   b(1)=s04*g0

   b(2)=s1*g0+s04*g1

   b(3)=s2*g0+2*s1*g1+s04*g2

   b(4)=s3*g0+3*s2*g1+3*s1*g2+s04*g3

   b(5)=s4*g0+4*s3*g1+6*s2*g2+4*s1*g3+s04*g4

   do j=0,4
      aa(1,j+1)=rtil**(l+1+2*j)
      aa(2,j+1)=(l+1+2*j)*rtil**(l+2*j)
      aa(3,j+1)=(l+1+2*j)*(l+2*j)*rtil**(l+2*j-1)
      aa(4,j+1)=(l+1+2*j)*(l+2*j)*(l+2*j-1)*rtil**(l+2*j-2)
      aa(5,j+1)=(l+1+2*j)*(l+2*j)*(l+2*j-1)*(l+2*j-2)*rtil**(l+2*j-3)
   end do
!-----------------------------------------------------------------------
   call gauss(aa,b,d2n,5,jjj)
!-----------------------------------------------------------------------


!--- bess = 0 -------------------------------------------------------------
   delx = acos(-1.d0)/100.d0

   ifrag = 1
   x=0.d0
   x0(0)=0.d0
   eps=0.1d-10

   do while ( ifrag < 6 )
      x=x+delx
      if(bess(l,x)*bess(l,x+delx) <= 0.d0) then
         x1=x
         x2=x+delx
         do while (abs(bess(l,x)) > eps)
            x=x1-bess(l,x1)/(bess(l,x2)-bess(l,x1))*(x2-x1)
            if(bess(l,x1)*bess(l,x) < 0.d0) then
               x2=x
            else if( bess(l,x2)*bess(l,x) < 0.d0) then
               x1=x
            end if
         end do
         x0(ifrag)=x
         ifrag=ifrag+1
      end if
   end do

!      do while ( ifrag .lt. irrkj+1 )
!         x=x+delx
!         if(bess(l,x)*bess(l,x+delx) .le. 0.d0) then
!            x0(ifrag)=x
!            ifrag=ifrag+1
!         end if
!      end do
   write(*,'(4x,a,i1,a)') 'zero point of sphere bessel function (l=',l,')'
   write(*,'(6x,a,6f10.6, /)') 'x0 = ', x0(0:5)



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   do nrep=1, nref
   !-----------------------------------------------------------------------
   !      if( eref(nrep) .lt. 0 ) then
      ! e=eref(nrep)
   !      end if
      write(*,'(4x,a,f15.10)') 'Reference energy :', eref(nrep)

      do i=1, icl + 50
         pae0(i) = plj(nrep,i)
      end do

   !--- derivative of P,  finite difference  -----------------------
      call fd1 (1, icl+10, mesh, pae0, dx, pae1)
      call fd2 (1, icl+10, mesh, pae0, dx, pae2)
      call fd3 (1, icl+10, mesh, pae0, dx, pae3)
   !      call fd3n (1, icl+10, mesh, pae0, dx, pae5)
      call fd4 (1, icl+10, mesh, pae0, dx, pae4)
   !      call fd4n (1, icl+10, mesh, pae0, dx, pae6)
   !-----------------------------------------------------------------------

   !--- sphere bessel function --------------------------------------------


      do j=0, irrkj-1
         fx(1)=x0(j)+delx
         fx(2)=fx(1)+delx
      ! x(2)=x0(j+1)-delx

         f(1) = 1.d0
         f(2) = 1.d0
         f(3) = 3.d0
         do while( f(1)*f(2) > 0.d0 )
            do i=1,2
               f(i) = -l+fx(i)*bess(l-1,fx(i))/bess(l,fx(i)) &
               -1.d0/pae0(inc)*pae1(inc)
            end do
            fx(1)=fx(1)+delx
            fx(2)=fx(1)+delx
            if(fx(1)==x0(j+1)) then
            !              stop 'error6'
               if(j==0) then
                  x0(j)   = x0(j+1)
                  x0(j+1) = x0(j+2)
                  x0(j+2) = x0(j+3)
                  x0(j+3) = x0(j+4)
                  write(*,'(6x,a)') 'x0(i) <=== x0(i+1),  i = 0, 3'
               else if(j==1) then
                  x0(j)   = x0(j+1)
                  x0(j+1) = x0(j+2)
                  x0(j+2) = x0(j+3)
                  write(*,'(6x,a)') 'x0(i) <=== x0(i+1),  i = 1, 3'
               else if(j==2) then
                  x0(j)   = x0(j+1)
                  x0(j+1) = x0(j+2)
                  write(*,'(6x,a)') 'x0(1,2) <=== x0(2,3)'
               end if
               fx(1)=x0(j)+delx
               fx(2)=fx(1)+delx
            end if
         end do
         fx(1)=fx(1)-delx
         fx(2)=fx(2)-delx

         do while(abs(f(3)) > eps)
            do i=1,2
               f(i) = -l+fx(i)*bess(l-1,fx(i))/bess(l,fx(i)) &
               -1.d0/pae0(inc)*pae1(inc)
            end do

            fx(3) = fx(1)-f(1)/(f(2)-f(1))*(fx(2)-fx(1))
            f(3) = -l+fx(3)*bess(l-1,fx(3))/bess(l,fx(3)) &
            -1.d0/pae0(inc)*pae1(inc)

            if(f(1)*f(3) < 0.d0) then
               fx(2)=fx(3)
            else if (f(2)*f(3) < 0.d0) then
               fx(1)=fx(3)
            else
               stop'error3'
            end if
         end do

         q(j+1)=fx(3)/rcl2
      enddo
   !------------------------------------------------------------------------
   !     Pnc,lj(r)
   !------------------------------------------------------------------------

   !----- ( 1 ) ------------------------------------------------------------

   !----- The fourth differentiation of P ---------------------------------

   !----- Flj(r), Flj~(r) -------------------------------------------------
      do i = 1, inc+6
         flj(i)  = 0.d0
         ftil(i) = 0.d0
         if( r(i) < rtil ) then
            do j=1, 5
               flj(i)=flj(i)+r(i)**(l+1+2*(j-1))*c2n(j)
               ftil(i)=ftil(i)+r(i)**(l+1+2*(j-1))*d2n(j)
            end do
         else if(r(i) >= rtil) then
            flj(i) = flj(i) +sin(qc*r(i))**3*exp(-gamm*(rcl2-r(i)))
            ftil(i) = ftil(i) +sin(qc*r(i))**4*exp(-gamm*(rcl2-r(i)))
         end if
      end do

      do i = 1, 3
         do k = i, 3
            do j = 1, inc
               qi = q(i)*r(j)
               qk = q(k)*r(j)
               y(j) = r(j)**3*bess(l,qi)*bess(l,qk)
            end do
         !         write(*,*)'y:',y(1),y(2)
            call intgbp( inc, dx, y, t)
            sl(i,k)=t
         end do
         do j = 1, inc
            qi = q(i)*r(j)
            y(j) = r(j)*r(j)*bess(l,qi)*flj(j)
         end do
         call intgbp( inc, dx, y, t)
         sl(i,4)=t

         do j = 1, inc
            qi = q(i)*r(j)
            y(j) = r(j)*r(j)*bess(l,qi)*ftil(j)
         end do
         call intgbp( inc, dx, y, t)
         sl(i,5)=t

      end do

      do j = 1, inc
         y(j) = r(j)*flj(j)*flj(j)
      end do
      call intgbp( inc, dx, y, t)
      sl(4,4)=t

      do j = 1, inc
         y(j) = r(j)*flj(j)*ftil(j)
      end do
      call intgbp( inc, dx, y, t)
      sl(4,5)=t

      do j = 1, inc
         y(j) = r(j)*ftil(j)*ftil(j)
      end do
      call intgbp( inc, dx, y, t)
      sl(5,5)=t


      do i=1, 5
         do j=i+1, 5
            sl(j,i)=sl(i,j)
         enddo
      enddo

   !----- ( 2 ) -------------------------------------------------------------

      ali(3)=0.d0


      b(1)=pae1(inc)/rcl2
      b(2)=(pae2(inc)-pae1(inc))/rcl2/rcl2


   !----- ( 3 ) -------------------------------------------------------------

      do i=1, 2
         t=q(i)*rcl2
         aa(1,i)=bess(l,t)+t*dbess(l,t,1)
         aa(2,i)=q(i)*(2*dbess(l,t,1)+t*dbess(l,t,2))
      end do


   !      write(*,*) 'inc, icl',inc,icl
   !      write(*,*) 'x0',x0(1),x0(2),x0(3)
   !      write(*,*) 'q(i)',q(1),q(2),q(3)
   !      write(*,*) aa(1,1), aa(1,2)
   !      write(*,*) aa(2,1), aa(2,2)
   !      write(*,*) b(1), b(2)

   !-----------------------------------------------------------------------
      call gauss(aa,b,ali,2,jjj)
   !-----------------------------------------------------------------------

      x=(pae3(inc)-3*pae2(inc)+2*pae1(inc))/rcl2**3
      ay=0.d0
      do i=1, 2
         t=q(i)*rcl2
         ay=ay+ali(i)*q(i)**2*(3*dbess(l,t,2)+t*dbess(l,t,3))
      end do
      ali(4)=(ay-x)/(6.d0*qc**3)

      x=(pae4(inc)-6*pae3(inc)+11*pae2(inc)-6*pae1(inc))/rcl2**4
      ay=0.d0
      do i=1, 2
         t=q(i)*rcl2
         ay=ay+ali(i)*q(i)**3*(4*dbess(l,t,3)+t*dbess(l,t,4))
      end do
      ali(5)=(x-ay+24*qc**3*gamm*ali(4))/(24*qc**4)

      do i=1,inc
         y(i)=r(i)*abs(pae0(i)*pae0(i))
      enddo

      call intgbp( inc, dx, y, pint)

      f(1)=0.d0

      do i=1, 5
         do j=1, 5
            f(1)=f(1)+ali(i)*ali(j)*sl(i,j)
         enddo
      enddo
      f(1)=f(1)-pint

   !----- ( 4 ) -------------------------------------------------------------


   !      s=0.d0
   !      do i=1, 5
   !         s=s+2*ali(i)*sl(i,3)
   !      enddo

   !----- ( 5 ) -------------------------------------------------------------

   !      if((f(1) > 0.d0 .and. s > 0.d0 ).or. (f(1) < 0.d0 .and. s < 0.d0))
   !     &                                            then
   !         del=-0.01d0
   !      else
   !         del=0.01d0
   !      end if

      del=-0.00001d0

      ! f1 = f(1)
      iter=1
      do while(f(1)*f(2) > 0 .OR. iter == 1)
         if( iter == 2 .AND. ( abs(f(1)) < abs(f(2)) ) )then
            del = -del
         end if

         if( iter /= 1 )then
            f(1)=f(2)
         end if


      !----- ( 6 ) --------------------------------------------------------------

         ali(3)=ali(3)+del
      !----- ( 7 ) ---------------------------------------------------------------
         t=q(3)*rcl2
         b(1)=pae1(inc)/rcl2-ali(3)*(bess(l,t)+t*dbess(l,t,1))
         b(2)=(pae2(inc)-pae1(inc))/rcl2/rcl2 &
         -ali(3)*(q(3)*(2*dbess(l,t,1)+t*dbess(l,t,2)))

         do i=1, 2
            t=q(i)*rcl2
            aa(1,i)=bess(l,t)+t*dbess(l,t,1)
            aa(2,i)=q(i)*(2*dbess(l,t,1)+t*dbess(l,t,2))
         end do

      !-----------------------------------------------------------------------
         call gauss(aa,b,ali,2,jjj)
      !-----------------------------------------------------------------------

         x=(pae3(inc)-3*pae2(inc)+2*pae1(inc))/rcl2**3
         ay=0.d0
         do i=1, 3
            t=q(i)*rcl2
            ay=ay+ali(i)*q(i)**2*(3*dbess(l,t,2)+t*dbess(l,t,3))
         end do
         ali(4)=(ay-x)/(6.d0*qc**3)

         x=(pae4(inc)-6*pae3(inc)+11*pae2(inc)-6*pae1(inc))/rcl2**4
         ay=0.d0
         do i=1, 3
            t=q(i)*rcl2
            ay=ay+ali(i)*q(i)**3*(4*dbess(l,t,3)+t*dbess(l,t,4))
         end do
         ali(5)=(x-ay+24*qc**3*gamm*ali(4))/(24*qc**4)

         f(2)=0.d0

         do i=1, 5
            do j=1, 5
               f(2)=f(2)+ali(i)*ali(j)*sl(i,j)
            enddo
         enddo
         f(2)=f(2)-pint

      !------ ( 8 ) ----------------------------------------------------------
         iter=iter+1
      !         write(*,*)f(1),f(2), e
            if( iter>=2000000 ) then
               write(*,'(/, a)')       ' ***** error : Did not decide qi.'
               write(*,'(a,2f18.7)')   ' x : ', fx(1:2)
               write(*,'(a,3es18.10)') ' f : ',  f(1:3)
               stop
            end if
      enddo
      write(*,'(6x,a,6f10.6, /)') 'qi = ', q(1:irrkj)

   !------ ( 9 ) ----------------------------------------------------------


      f(3)=1.d0
      iter=1
      ali3(1)=ali(3)-del
      ali3(2)=ali(3)
      do while (abs( f(3) ) > eps)

         if(( f(1)*f(3) < 0.d0 ) .AND. (iter /= 1)  ) then
            ali3(2) = ali(3)
            f(2)   = f(3)
         else if( f(2)*f(3) < 0.d0 .AND. iter /= 1) then
            ali3(1) = ali(3)
            f(1)   = f(3)
         else if (iter /= 1) then
            write(*,*)f(1),f(2),f(3),iter
             stop ' ***** error : Did not decide alpha(3).'
         end if


      !------ ( a ) ----------------------------------------------------------

         ali(3)=ali3(1)-f(1)/(f(2)-f(1))*(ali3(2)-ali3(1))

         t=q(3)*rcl2
         b(1)=pae1(inc)/rcl2-ali(3)*(bess(l,t)+t*dbess(l,t,1))
         b(2)=(pae2(inc)-pae1(inc))/rcl2/rcl2 &
         -ali(3)*q(3)*(2*dbess(l,t,1)+t*dbess(l,t,2))

         do i=1, 2
            t=q(i)*rcl2
            aa(1,i)=bess(l,t)+t*dbess(l,t,1)
            aa(2,i)=q(i)*(2*dbess(l,t,1)+t*dbess(l,t,2))
         end do
      !-----------------------------------------------------------------------
         call gauss(aa,b,ali,2,jjj)
      !-----------------------------------------------------------------------

         x=(pae3(inc)-3*pae2(inc)+2*pae1(inc))/rcl2**3
         ay=0.d0
         do i=1, 3
            t=q(i)*rcl2
            ay=ay+ali(i)*q(i)**2*(3*dbess(l,t,2)+t*dbess(l,t,3))
         end do
         ali(4)=(ay-x)/(6.d0*qc**3)

         x=(pae4(inc)-6*pae3(inc)+11*pae2(inc)-6*pae1(inc))/rcl2**4
         ay=0.d0
         do i=1, 3
            t=q(i)*rcl2
            ay=ay+ali(i)*q(i)**3*(4*dbess(l,t,3)+t*dbess(l,t,4))
         end do
         ali(5)=(x-ay+24*qc**3*gamm*ali(4))/(24*qc**4)


         f(3)=0.d0

         do i=1, 5
            do j=1, 5
               f(3)=f(3)+ali(i)*ali(j)*sl(i,j)
            enddo
         enddo
         f(3)=f(3)-pint
         iter=iter+1

      end do
      write(*,'(6x,a,5es14.6)') 'alpha =', ali(1:5)
   ! eck-----------------------------------
   !      ali(4)=0
   !      ali(5)=0.d0
   !---------------------------------------
      do i=1, icl+50
         if(i < inc) then
            pnc0(i)=0.d0
            do j=1,3
               pnc0(i)=pnc0(i)+ali(j)*r(i)*bess(l,q(j)*r(i))
            end do
            if(r(i) < rtil) then
               do j=1, 5
                  pnc0(i)=pnc0(i)+ali(4)*r(i)**(l+1+2*(j-1))*c2n(j) &
                  +ali(5)*r(i)**(l+1+2*(j-1))*d2n(j)
               end do
            else if(r(i) >= rtil) then
               pnc0(i)=pnc0(i)+ali(4)*sin(qc*r(i))**3 &
               *exp(-gamm*(rcl2-r(i))) &
               + ali(5)*sin(qc*r(i))**4*exp(-gamm*(rcl2-r(i)))
            end if
         else if (i >= inc) then
            pnc0(i) = pae0(i)
         end if
      end do

      call fd1(1,inc+6,msh,pnc0,dx,pnc1)
      call fd2(1,inc+6,msh,pnc0,dx,pnc2)
      call fd3(1,inc+6,msh,pnc0,dx,pnc3)
      call fd4(1,inc+6,msh,pnc0,dx,pnc4)

      do i=1,inc
         y(i)=r(i)*pnc0(i)*pnc0(i)
      end do
      call intgbp( inc, dx, y, t )

      write(*,'(22x,a)') 'P_US              P_AE              diff.'
      write(*,'(6x,a,3es18.10)')  &
      &     '1st derivs. : ', pnc1(inc), pae1(inc), abs(pnc1(inc) - pae1(inc))
      write(*,'(6x,a,3es18.10)')  &
      &     '2nd derivs. : ', pnc2(inc), pae2(inc), abs(pnc2(inc) - pae2(inc))
      write(*,'(6x,a,3es18.10)')  &
      &     '3rd derivs. : ', pnc3(inc), pae3(inc), abs(pnc3(inc) - pae3(inc))
      write(*,'(6x,a,3es18.10)')  &
      &     '4th derivs. : ', pnc4(inc), pae4(inc), abs(pnc4(inc) - pae4(inc))
      write(*,'(6x,a,3es18.10, /)')  &
      &     'norm [0,rc] : ', t, pint, abs( t - pint)
   !         t=0
   !         do j=1,2
   !           t=t+ali(j)*rcl2*bess(l,q(j)*rcl2)
   !         end do
   !      write(*,*)r(i),rcl2
   !      write(*,*)t,pae0(inc)
   !      do i=7, inc
   !      write(9,'(6e18.10)')r(i),pnc0(i), pnc1(i),pnc2(i),pnc3(i),pnc4(i)
   !      end do
      close(9)

   !----------------------------------------------------------------------
      do i=1, icl+50
         pae(nrep,i)=pae0(i)
         pnc(nrep,i)=pnc0(i)
      end do
   !---------------------------------------------------------------------
   end do
   ofile=(ms(l+1)//'_Pnc.dat')
   open(9,file=ofile)
   do i = 1, icl + 50
!       write(9,*)r(i),pnc0(i)
      write(9,*) r(i), ( pnc(j,i), j = 1, nref )
   end do
   write(*,'(4x,2a, /)') 'Create : ', ofile
!----------------------------------------------------------------------



   return
   end SUBROUTINE rrkj3




subroutine ulps ( r, l, iflag )
!-----------------------------------------------------------------------
!  input :
!         vloc  ...... local potential
!         pae   ...... the radial AE wave functions
!         pus   ...... the US pseudo-wave-functions
!         pnc   ...... the NC pseudo-wave-functions
!         iflag ...... 1 --> USPP, else --> PAW
!
!  out put :
!         chi   ...... local functions
!         bele  ...... matrix elements B
!         bas   ...... basis functions
!         qau   ...... augmentation functions
!         dele  ...... matrix elements D for USPP
!         dele2 ...... matrix elements D for PAW
!
!-----------------------------------------------------------------------
   use common_variables
   implicit real(8) (a-h,o-z)
   real(8) :: r(msh)
   real(8) :: y(msh)
   real(8) :: bele(5,5), bret(5,5)
   character(13) :: ofile
   integer :: ipair(mxref*(mxref+1)/2), jpair(mxref*(mxref+1)/2)
   character(1) :: num(0:9) = (/ '0','1','2','3','4','5','6','7','8','9' /)

   write(*,'(30a)') ' ---- Prepare for Construction of ',  &
   &                'pseudo-potentials ',  &
   &                ( '-', i = 1, 26 )

!----- local functions -------------------------------------------------
   cll1 = dble(l*(l+1))
   do i = 1, mesh
      do j = 1, nref
         if( i <= ius ) then
!              chi(j,i) = d2p(j,i) + (-(l*(l+1)/r(i)/r(i))+eref(j)-vloc(i))*pus(j,i)
             ri2 = r(i)*r(i)
             chi(j,i) = d2p(j,i) - cll1/ri2*pus(j,i) + (eref(j) - vloc(i))*pus(j,i)
         else if( i > ius ) then
             chi(j,i) = 0.d0
         end if
      end do
   end do

!----- matrix elements -------------------------------------------------
   do i=1, nref
      do j=1, nref
         do k=1, icl+2
            y(k) = pus(i,k) * chi(j,k) * r(k)
         end do
         call intgbp( icl+2, dx, y , t)
         bele(i,j) = t
      end do
   end do
!----- basis functions -------------------------------------------------
   if(      nref == 1 ) then
      bret(1,1) = 1.d0/bele(1,1)
   else if( nref == 2 ) then
      a = bele(1,1)
      b = bele(1,2)
      c = bele(2,1)
      d = bele(2,2)
      t = a*d - b*c
      bret(1,1) =  d/t
      bret(1,2) = -b/t
      bret(2,1) = -c/t
      bret(2,2) =  a/t
   end if

   do j = 1, nref
      do i = 1, mesh
         bas(j,i) = 0.d0
         do k = 1, nref
            bas(j,i) = bas(j,i) + bret(k,j)*chi(k,i)
         end do
      end do
   end do

   selectif: if( iflag == 1 ) then
!----------------------------------------------------------------------!
!       Ultrasoft pseudopotential                                      !
!----------------------------------------------------------------------!
!----- augmentation functions ------------------------------------------
      do i = 1, mesh
         do j = 1, nref
            do k = 1, nref
               if( i <=icl ) qau(j,k,i) = pnc(j,i)*pnc(k,i) - pus(j,i)*pus(k,i)
               if( i > icl ) qau(j,k,i) = 0.d0
            end do
         end do
      end do

!----- augmentation charges --------------------------------------------
      do j = 1, nref
         do k = 1, nref
            do i = 1, icl+2
               y(i) = qau(j,k,i)*r(i)
            end do
            call intgbp( icl+2, dx, y, t )
            qch(j,k) = t
         end do
      end do

!----- matrix elements -------------------------------------------------
      do j = 1, nref
         do k=1, nref
            dele(j,k) = bele(j,k) + eref(k)*qch(j,k)
         end do
      end do
   else selectif
!----------------------------------------------------------------------!
!       Projector augmented wave method                                !
!----------------------------------------------------------------------!
!----- augmentation functions for Projector augmented wave method -------------
      do i = 1, mesh
         do j = 1, nref
            do k = 1, nref
               if( i <=icl) qau(j,k,i) = pae(j,i)*pae(k,i) - pus(j,i)*pus(k,i)
               if( i > icl) qau(j,k,i) = 0.d0
            end do
         end do
      end do

!----- augmentation charges for Projector augmented wave method ---------------
      do j = 1, nref
         do k = 1, nref
            do i = 1, icl+2
               y(i) = qau(j,k,i)*r(i)
            end do
            call intgbp( icl+2, dx, y, t )
            qch(j,k) = t
            do lll = 0, 2*l, 2
               do i = 1, nrcut+2
                  y(i) = qau(j,k,i) * r(i)**(lll+1)
               end do
               call intgbp( nrcut+2, dx, y, t )
               qcomp(lll,j,k) = t
               do i = 1, mesh
                  ri2 = r(i)*r(i)
                  if( i<=nrcut ) qau2(lll,j,k,i) = qcomp(lll,j,k)*gf(lll,i)*ri2
                  if( i> nrcut ) qau2(lll,j,k,i) = 0.d0
               end do
            end do
!             qch(j,k) = qcomp(0,j,k)
         end do
      end do

!----- matrix elements -------------------------------------------------
      do j = 1, nref
         do k = 1, nref
            dele(j,k) = bele(j,k) + eref(k)*qch(j,k)
         end do
      end do
   end if selectif

! ----- output
   napir = 0
   do i = 1, nref, 1
      do j = i, nref, 1
         napir = napir + 1
         ipair(napir) = i
         jpair(napir) = j
      end do
   end do

   ofile=(ms(l+1)//'_chi.dat')
   open(9,file=ofile)
   ofile=(ms(l+1)//'_beta.dat')
   open(10,file=ofile)

   do i = 1, icl+50, 1
      write( 9,*) r(i), ( chi(j,i), j = 1, nref)
      write(10,*) r(i), ( bas(j,i), j = 1, nref)
   end do

   close(9)
   ofile=(ms(l+1)//'_chi.dat')
   write(*,'(4x,2a)') 'Create : ', ofile
   close(10)
   ofile=(ms(l+1)//'_beta.dat')
   write(*,'(4x,2a)') 'Create : ', ofile

   if( iflag == 1 ) then
       ofile=(ms(l+1)//'_Q.dat')
       open(11,file=ofile)
       do i = 1, icl+50, 1
          write(11,'(20es18.10)') r(i), ( qau(ipair(k),jpair(k),i), k = 1, napir )
       end do
       close(11)
      write(*,'(4x,2a, /)') 'Create : ', ofile
   else
      do lll = 0, 2*l, 2
         ofile=(ms(l+1)//'_Q_L='//num(lll)//'.dat')
         open(11,file=ofile)
         do i = 1, icl+50, 1
            write(11,'(20es18.10)') r(i), ( qau2(lll,ipair(k),jpair(k),i), k = 1, napir )
         end do
         close(11)
         write(*,'(4x,2a)') 'Create : ', ofile
      end do
   end if

   return
end subroutine ulps



! subroutine gee( zatm, n, L, R, e)
subroutine gee( n, L, R, e )
!-----------------------------------------------------------------------
!     Generalization eigenvalue equation
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   parameter     (jjj=10)
   dimension     r(msh), s(msh)
   dimension     y(msh),ula(msh),usm(5,5),beti(5,msh),wla(msh)
   dimension     wsm(5,msh), vla(msh),vsm(5,msh),vla1(msh)
   dimension     vsm1(5,msh)
   dimension     wtla(5),wtsm(5,5)
   dimension     aa(jjj,jjj),bb(jjj),aaa(jjj) , f(msh)
   dimension     pnl(msh),t(msh) ,d(msh)
   CHARACTER     OFILE*13
!-----------------------------------------------------------------------
   write(*,'(/, 50a)') ' ---- Solve the Generalized eigenequations ',  &
   &                   ( '-', i = 1, 35 )
   CLL1  =   DBLE( L * ( L + 1 ) )
   DELK  =   0.01D0
   EMIN  =  - ( ZATM / DBLE( N ) )**2
   EMAX  =  - 0.00001D0
   IF( ( E > EMAX ) .OR. ( E < EMIN ) ) THEN
      E     =  0.5D0 * ( EMAX + EMIN )
   ENDIF
   1 continue
   do i=1, mesh
      ula(i) = (l+0.5d0)**2.d0+(vloc(i)-e)*r(i)*r(i)
      do j=1, nref
         if ( i <= icl)beti(j,i) = r(i)**(1.5d0)*bas(j,i)
         if ( i > icl) beti(j,i) = 0.d0
      end do
   end do

   ! pi = dacos( -1.d0 )
   usm =  0.d0
   do j=1, nref
      do k=1, nref
         usm(j,k) = dele(j,k) - e*qch(j,k)

      !            usm(j,k) = usm(j,k) + dele2(0,j,k)
      !     &                - e*qcomp(0,j,k)

      end do
   end do


   do i=1,5
      do j=1,5
         aa(i,j)=(ula(i)-(l+j+0.5d0)*(l+j+0.5d0))*r(i)**(l+j+0.5d0)
      end do
      bb(i)=-(ula(i)-(l+0.5d0)**2)*r(i)**(l+0.5d0)
   end do

   call gauss(aa,bb,aaa,5,jjj)

   do i=1,5
      wla(i) = r(i)**(l+0.5)                                       !Y
      vla(i) = (l+0.5)*r(i)**(l+0.5)                                !Y'
   !         vla1(i)= (l+0.5)**2*r(i)**(l+0.5)

      do j=1,5
         wla(i)=wla(i)+aaa(j)*r(i)**(l+j+0.5d0)
         vla(i)=vla(i)+(l+j+0.5)*aaa(j)*r(i)**(l+j+0.5d0)
      !            vla1(i)=vla1(i)+(l+j+0.5)**2.d0*aaa(j)*r(i)**(l+j+0.5d0)
      end do

      vla1(i)=ula(i)*wla(i)                    !Y''
   !         write(*,*)'f(i)',f(i),vla1(i)
   end do
!--- integrates outward to classical turning point ---------------------
!--- determine turning point ( KI ) ---
   DO 33 I = 3, MESH - 10
      KM  =  MESH - I
      IF( E * R(KM) - vloc(km) *r(km) - CLL1 / R(KM) > 0.D0 ) GOTO 34
   !      write(*,*)'ki:',r(km), E * R(KM) - vloc(km) *r(km) - CLL1 / R(KM)
   !      IF( E * R(KM) - VR(KM) - CLL1 / R(KM) .GT. 0.D0 ) GOTO 34
   33 END DO
   34 KI  =  KM + 1
   KI  =  max(KI, icl+10)
!      WRITE(*,*) 'ki, e:',ki,e
!--- Adams-Bashforth Adams-Moulton -------------------------------------

   do i = 5, ki + 6
      wla(i+1)=wla(i)+dx/720*(1901*vla(i)-2774*vla(i-1)+2616*vla(i-2) &
      -1274*vla(i-3)+251*vla(i-4))
      vla(i+1)=vla(i)+dx/720*(1901*vla1(i)-2774*vla1(i-1) &
      +2616*vla1(i-2)-1274*vla1(i-3)+251*vla1(i-4))
      vla1(i+1)=ula(i+1)*wla(i+1)

      wla(i+1)=wla(i)+dx/720*(251*vla(i+1)+646*vla(i)-264*vla(i-1) &
      +106*vla(i-2)-19*vla(i-3))
      vla(i+1)=vla(i)+dx/720*(251*vla1(i+1)+646*vla1(i)-264*vla1(i-1) &
      +106*vla1(i-2)-19*vla1(i-3))
      vla1(i+1)=ula(i+1)*wla(i+1)
   enddo


!-----wlj(x)-----------------------------------------------------------
   do k = 1, nref
      do i=1,5
         do j=1,5
            aa(i,j)=(ula(i)-(l+(j-1)+0.5d0)*(l+j-1+0.5d0)) &
            *r(i)**(l+j-1+0.5d0)
         !         write(*,*)'aa:',aa(i,j)
         end do
         bb(i)=beti(k,i)
      !         write(*,*)'b:',bb(i)
      end do

      call gauss(aa,bb,aaa,5,jjj)

      do i=1,5
         wsm(k,i) = r(i)**(l+0.5)*aaa(1)                           !Y
         vsm(k,i) = (l+0.5)*r(i)**(l+0.5)*aaa(1)                   !Y
      !         vsm1(k,i)= (l+0.5)**2*r(i)**(l+0.5)*aaa(1)
         do j=2,5
            wsm(k,i)=wsm(k,i)+aaa(j)*r(i)**(l+j-0.5d0)
            vsm(k,i)=vsm(k,i)+(l+j-0.5)*aaa(j)*r(i)**(l+j-0.5d0)
         !           vsm1(k,i)=vsm1(k,i)+(l+j-0.5)**2.d0*aaa(j)*r(i)**(l+j-0.5d0)
         enddo
         vsm1(k,i)=ula(i)*wsm(k,i)-beti(k,i)                    !Y''
      !         write(*,*)'f(i)',vsm1(k,i)
      !         write(*,*)'a:',aaa(i)

      enddo

   !      stop

   !--- integrates outward to classical turning point ---------------------
   !--- Adams-Bashforth Adams-Moulton -------------------------------------

      do i = 5, ki + 6
         wsm(k,i+1)=wsm(k,i)+dx/720*(1901*vsm(k,i)-2774*vsm(k,i-1) &
         +2616*vsm(k,i-2) &
         -1274*vsm(k,i-3)+251*vsm(k,i-4))
         vsm(k,i+1)=vsm(k,i)+dx/720*(1901*vsm1(k,i)-2774*vsm1(k,i-1) &
         +2616*vsm1(k,i-2) &
         -1274*vsm1(k,i-3)+251*vsm1(k,i-4))
         vsm1(k,i+1)=ula(i+1)*wsm(k,i+1)-beti(k,i+1)

         wsm(k,i+1)=wsm(k,i)+dx/720*(251*vsm(k,i+1)+646*vsm(k,i) &
         -264*vsm(k,i-1) &
         +106*vsm(k,i-2)-19*vsm(k,i-3))
         vsm(k,i+1)=vsm(k,i)+dx/720*(251*vsm1(k,i+1)+646*vsm1(k,i) &
         -264*vsm1(k,i-1) &
         +106*vsm1(k,i-2)-19*vsm1(k,i-3))
         vsm1(k,i+1)=ula(i+1)*wsm(k,i+1)-beti(k,i+1)
      enddo
   enddo

!-----------------------------------------------------------------------
!                     starts inward integration
!                     outer boundary condition is set here.
!-----------------------------------------------------------------------
!--- determine the upper limit of the wave function ( KJ ) ---
   DO 21 I = KI, MESH
      KJ   =  I
   !C    IF( ( E * R(I) - VR(I) ) * R(I) + RINF ) 22, 21, 21
      IF( ( E   - vloc(i) )* R(i) + 2.D0 * RINF ) 22, 21, 21
   21 END DO
   22 CONTINUE
!      kj = mesh
!--- asymptotic form of the wave function ---
   s(KJ)    =  1.D-20                                        !Y

   t(kj)    =  -(0.5+sqrt(-e)*r(kj))*s(kj)                  !xx


   do i = kj-1, kj-4, -1
      s(i) = s(kj)*sqrt(r(kj)/r(i))*exp(sqrt(-e)*(r(kj)-r(i)))
      t(i) = -(0.5d0+sqrt(-e)*r(i))*s(i)
   end do

   do j = 1, 5

      do i = kj, kj-4,-1
         d(i)=ula(i)*s(i)                    !Y''
      end do

      s(kj-1) = s(kj)+(-dx)/720.d0*(251.d0*t(kj)+646.d0*t(kj-1) &
      -264.d0*t(kj-2)+106.d0*t(kj-3)-19.d0*t(kj-4))
      t(kj-1) = t(kj)+(-dx)/720.d0*(251.d0*d(kj)+646.d0*d(kj-1) &
      -264.d0*d(kj-2)+106.d0*d(kj-3)-19.d0*d(kj-4))
      s(kj-2) = s(kj)+(-dx)/90.d0*(29.d0*t(kj)+124*t(kj-1) &
      +24.d0*t(kj-2)+4.d0*t(kj-3)-t(kj-4))
      t(kj-2) = t(kj)+(-dx)/90.d0*(29.d0*d(kj)+124.d0*d(kj-1) &
      +24.d0*d(kj-2)+4.d0*d(kj-3)-d(kj-4))
      s(kj-3) = s(kj)+(-dx)/240.d0*(81.d0*t(kj)+306.d0*t(kj-1) &
      +216.d0*t(kj-2)+126.d0*t(kj-3)-9.d0*t(kj-4))
      t(kj-3) = t(kj)+(-dx)/240.d0*(81.d0*d(kj)+306.d0*d(kj-1) &
      +216.d0*d(kj-2)+126.d0*d(kj-3)-9.d0*d(kj-4))
      s(kj-4) = s(kj)+(-dx)/180.d0*(56.d0*t(kj)+256.d0*t(kj-1) &
      +96.d0*t(kj-2)+256.d0*t(kj-3)+56.d0*t(kj-4))
      t(kj-4) = t(kj)+(-dx)/180.d0*(56.d0*d(kj)+256.d0*d(kj-1) &
      +96.d0*d(kj-2)+256.d0*d(kj-3)+56.d0*d(kj-4))
   end do


!--- Adams-Bashforth Adams-Moulton -------------------------------------

   do i = kj-4, ki - 1, -1
      s(i-1)=s(i)-dx/720.d0*(1901.d0*t(i)-2774.d0*t(i+1) &
      +2616.d0*t(i+2)-1274.d0*t(i+3)+251.d0*t(i+4))
      t(i-1)=t(i)-dx/720.d0*(1901.d0*d(i)-2774.d0*d(i+1) &
      +2616.d0*d(i+2)-1274.d0*d(i+3)+251.d0*d(i+4))
      d(i-1)=ula(i-1)*s(i-1)

      s(i-1)=s(i)-dx/720.d0*(251.d0*t(i-1)+646.d0*t(i) &
      -264.d0*t(i+1)+106.d0*t(i+2)-19.d0*t(i+3))
      t(i-1)=t(i)-dx/720.d0*(251.d0*d(i-1)+646.d0*d(i) &
      -264.d0*d(i+1)+106.d0*d(i+2)-19.d0*d(i+3))
      d(i-1)=ula(i-1)*s(i-1)
   enddo




   do k = 1, nref
      do i = 1, mesh
         if (i >= icl) then
            f(i) = 0.d0

         else if (i <= ki) then
            f(i) = beti(k,i) * wla(i)

         else
         !               f(i) = beti(k,i) * s(i)
            f(i) = beti(k,i) * wla(i)
         end if
      enddo
      call intgbp( mesh, dx, f, tint1)
      wtla(k) = tint1
   !      write(*,*)'check-30:',tint1
   enddo
!    write(*,*)'icl, ki:',icl ,ki

   do k = 1, nref
      do i = 1, nref
      !            do j =1, mesh
      !            if (j >= icl) then
      !               f(j) = 0.d0
      !            else if (j <= ki) then
      !               f(j) = beti(k,j) * wsm(i,j)
      !            else
      !!               f(j) = beti(k,j) * s(j)
      !!               stop 'error5'
      !               f(j) = beti(k,j) * wsm(i,j)
      !            endif
      !            end do
         do j = 1, icl-1
            f(j) = beti(k,j) * wsm(i,j)
         end do
         do j = icl, mesh
            f(j) = 0
         end do
         call intgbp (mesh, dx, f, tint1)

      !      write(*,*)'check-3',tint1
         wtsm(k,i) = tint1
      enddo
   enddo


   do j = 1, nref
      do i = 1, nref
         aa(j,i) = 0.d0
         do k = 1, nref
            aa(j,i) = aa(j,i) + usm(j,k) * wtsm(k,i)
         end do
         if ( i == j ) aa(j,i) = aa(j,i) + 1.d0
      !      write(*,*)'aa',aa(j,i)
      enddo
   enddo

   do j = 1, nref
      bb(j) = 0.d0
      do k = 1, nref
         bb(j) = bb(j) -usm(j,k) * wtla(k)
      enddo
   !      write(*,*)'bb',bb(j)
   enddo

   call gauss(aa,bb,aaa,nref,jjj)
!      write(*,*)'check-2',aaa(1),aaa(2)
   do i = 1, ki+6
      y(i) = wla(i)
      do j =1, nref
         y(i) = y(i) + aaa(j) * wsm(j,i)
      enddo
   enddo

   yki = y(ki)


! eck
!       if(l==2) then
!!      open(9,file='yout.txt')
!      do i = 1 , ki
!         write(81,*)r(i), y(i)!, wla(i),wsm(1,i)!,wsm(2,i)
!      end do
!!      close(9)
!       stop
!       end if

!      stop
!--- LOG - DER at KI ---
   I  =  KI

   dyki = vla(ki)
   do j = 1, nref
      dyki = dyki + aaa(j) * vsm(j,ki)
   enddo
! eck
!      call fd1(ki-1,ki,mesh,y,dx,vla)

!      write(*,*)'wsm',dyki,vla(ki)
!      stop

   RLDM  =  ( dyki / y(i) - 0.5d0 ) / r(i)

!--- LOG - DER at KI ---

   RLDP  =  ( t(i) / s(i) - 0.5d0 ) / r(i)

!--- continuation of the wave function ---
   RATIO  =  YKI / S(KI)
!      WRITE(*,*) 'RATIO :',RATIO,YKI,S(KI)
   DO 5 I = KI, KJ
      Y(I)   =  S(I) * RATIO
   5 END DO
!--- tail ---
   IF( KJ < MESH ) THEN
      YKJR  =  Y(KJ) * SQRT( R(KJ) )
      DO 7 I = KJ, MESH
         X  =  ( R(KJ) - R(I) ) * SQRT( -E )
         IF( ABS( X ) < 80.D0 ) THEN
            W  =  EXP( X )
         ELSE
            W  =  0.D0
         ENDIF
         Y(I)  =  YKJR * W / SQRT( R(I) )
      7 END DO
   ENDIF
!--- normalization constant --------------------------------------------
   do i = 1, mesh
      p(i) = y(i) * sqrt(r(i))
   enddo

!      open(9,file='pi.txt')
!      do i = 1, mesh
!         write(9,*)r(i), p(i)
!      end do
!      close(9)

   do i = 1, mesh
      f(i) = p(i)*p(i)*r(i)
   enddo
   call intgbp( mesh, dx, f, tint1)
   w=tint1

   do j = 1, nref
      do i = 1, mesh
         f(i) = p(i) * bas(j,i) * r(i)
      end do
      call intgbp ( mesh, dx, f, tint1)

      do k = 1, nref
         do i = 1, mesh
            f(i) = p(i) * bas(k,i) * r(i)
         enddo
         call intgbp ( mesh, dx, f, tint2)

         w = w + qch(j,k) * tint1 * tint2

      enddo
   enddo

!    write(*,*)'check0 e:',e

!--- energy correction due to diff. of LOG - DER ---
   DE  =  R(KI) * Y(KI) * Y(KI) * ( RLDM - RLDP ) / W
!        write(*,*)'check-1',de,r(ki),y(ki),rldm,rldp,w

   DL1 =  ABS( EMAX - EMIN ) / ( ABS( E ) + 1.D0 )
   DL  =  ABS( DE / E )
   IF( ( DL > DELK ) .AND. ( DL1 < DELK ) ) GOTO 10
   IF( DL1 < DELL ) GOTO 10
   IF( DE  > 0.D0 ) THEN
      EMIN  =  E
!       write(*,*)'check1 e:',e
   ELSE
      EMAX  =  E
!       write(*,*)'check2 e:',e
   ENDIF
   DEP   =  DE
   E  =  E + DEP
!    write(*,*)'check3 e:',e
   IF( DL  > DELL ) THEN
      IF( E > EMAX ) THEN
         E  =  0.5D0 * ( E - DEP + EMAX )
!          write(*,*)'check4 e:',e
      ENDIF
      IF( E < EMIN ) THEN
         E  =  0.5D0 * ( E - DEP + EMIN )
!          write(*,*)'check5 e:',e
      ENDIF
      GOTO 1
   ENDIF
!--- eigenvalue is conserved -------------------------------------------
   10 CR    =  1.D0 / SQRT( W )
   DO 11 I = 1, MESH
      PNL(I)  =  CR * Y(I) * SQRT( R(I) )
      pusnl(i)=pnl(i)
   11 END DO
   EIG  =  E

   ofile=(ms(l+1)//'_Pus,nl.dat')
   open(9,file=ofile)

   do i = 1, mesh
      write(9,*)r(i), pnl(i)
   end do
   write(*,'(4x,a,f15.10)') 'eig  : ', eig
   close(9)
!--- count the number of nodes ---
   N1  =  0
   DO 35 K = 5, KI
      JK  =  INT( SIGN( 1.D0, Y(K)   ) )
      JK1 =  INT( SIGN( 1.D0, Y(K-1) ) )
      N1  =  N1 + ABS( JK - JK1 )
   35 END DO
   write(*,'(4x,a,i3)')'node : ', n1/2
   RETURN
   end subroutine gee




   SUBROUTINE GFUNC( l, NRCUT, GF, WORK, R, DX, MESH )
!-----------------------------------------------------------------------
!                                                          2002/02/04
!     calculate g functions (GF)

!-----------------------------------------------------------------------
   implicit real*8 ( a-h, o-z )
   DIMENSION  WORK(*), GF(*), R(*)
   logical :: lnbzr
   real*8 :: wm(4,5)
   integer :: ip(4)

   errcl  = 1.d-10
   zerop1 = 0.d0
   delx = acos(-1.d0)/100.d0

   if( l == 0 ) then
      xx1  = -0.5d0*delx
      val1 = bess( 1, xx1 )
   else
      xx1  = 1.d0
      val1 = -dble(l+1)*bess( l, xx1 )/xx1 + bess( l-1, xx1 )
   end if

   lnbzr = .TRUE.

!------- search zero-point -----------------
   outerdo: do
      do
      xx2 = xx1 + delx

      if( l == 0 ) then
         val2 = bess( 1, xx2 )
      else
         val2 = -dble(l+1)*bess( l, xx2 )/xx2 + bess( l-1, xx2 )
      end if

      if( val1*val2 < 0.d0 ) exit
      xx1  = xx2
      val1 = val2
   end do

   do
   xest = xx1 + ( xx2 - xx1 )*val1/( val1 - val2 )

   if( l == 0 ) then
      vest = bess( 1, xest )
   else
      vest = -dble(l+1)*bess(l, xest)/xest + bess(l-1, xest)
   end if

   if( vest*val1 < 0.d0 ) then
      xx2  = xest
      val2 = vest
   else
      xx1  = xest
      val1 = vest
   endif
   if( abs(vest) < errcl ) exit
end do
!-------------------------------------------
   zerop2 = xest
   xx1  = xx2
   val1 = val2
   if( .NOT. lnbzr ) exit outerdo
   lnbzr = .FALSE.
   zerop1 = zerop2
end do outerdo

! heck
!      write(*,*) zerop1, zerop2

   rcomp = R(NRCUT)
   q1 = zerop1/rcomp
   q2 = zerop2/rcomp


!-----determine c_i, d_i, i = 0, 6
   RCL = rcomp
   FQL = acos(-1.d0)/RCL
   gbeta  = 7.d0/RCL
   rhcl = RCL * 0.1d0
   RCLrhcl = RCL - rhcl
   epower = gbeta*(RCLrhcl)
   if( epower < 90.d0 ) then
      sinr = sin(FQL*rhcl)
      cosr = cos(FQL*rhcl)
      sinr2 = sinr *sinr
      sinr3 = sinr2*sinr
      cosr2 = cosr *cosr
      cosr3 = cosr2*cosr

      FQL2 = FQL *FQL
      FQL3 = FQL2*FQL

      g0r   = exp(-epower)
      g1r   = gbeta*g0r
      g2r   = gbeta*g1r
      g3r   = gbeta*g2r

      e0r   = sinr2
      e1r   = 2.d0*FQL*sinr*cosr
      e2r   = 2.d0*FQL2*(cosr2 - sinr2)
      e3r   = 2.d0*FQL3*(-4.d0*cosr*sinr)

      f0r   = sinr3
      f1r   = 3.d0*FQL*sinr2*cosr
      f2r   = 3.d0*FQL2*(2.d0*sinr*cosr2 - sinr3)
      f3r   = 3.d0*FQL3*(2.d0*cosr3 - 7.d0*sinr2*cosr)

      nbx = 4
      ln = l
      do j = 1, nbx
         wm(1,j) = rhcl**( ln + 2*j-2 )
      end do
      do i = 2, nbx
         do j = 1, nbx
            wm(i,j) = dble( ln + 2*j-2-i+2 )*wm(i-1,j)/rhcl
         end do
      end do
      wm(1,nbx+1) = e0r*g0r
      wm(2,nbx+1) = e1r*g0r + e0r*g1r
      wm(3,nbx+1) = e2r*g0r + 2.d0*e1r*g1r + e0r*g2r
      wm(4,nbx+1) = e3r*g0r + 3.d0*e2r*g1r + 3.d0*e1r*g2r + e0r*g3r
      CALL GSSJOR( WM, IP, nbx, 4, 5, INDER )
      IF( INDER /= 0 ) return
      c_0 = wm(1,nbx+1)
      c_2 = wm(2,nbx+1)
      c_4 = wm(3,nbx+1)
      c_6 = wm(4,nbx+1)

      do j = 1, nbx
         wm(1,j) = rhcl**( ln + 2*j-2 )
      end do
      do i = 2, nbx
         do j = 1, nbx
            wm(i,j) = dble( ln + 2*j-2-i+2 )*wm(i-1,j)/rhcl
         end do
      end do
      wm(1,nbx+1) = f0r*g0r
      wm(2,nbx+1) = f1r*g0r + f0r*g1r
      wm(3,nbx+1) = f2r*g0r + 2.d0*f1r*g1r + f0r*g2r
      wm(4,nbx+1) = f3r*g0r + 3.d0*f2r*g1r &
      + 3.d0*f1r*g2r + f0r*g3r
      CALL GSSJOR( WM, IP, nbx, 4, 5, INDER )
      IF( INDER /= 0 ) return
      d_0 = wm(1,nbx+1)
      d_2 = wm(2,nbx+1)
      d_4 = wm(3,nbx+1)
      d_6 = wm(4,nbx+1)

   else
      c_0 = 0.d0
      c_2 = 0.d0
      c_4 = 0.d0
      c_6 = 0.d0
      d_0 = 0.d0
      d_2 = 0.d0
      d_4 = 0.d0
      d_6 = 0.d0
   end if


!----- condition for g_l(rcomp) = 0
   wm(1,1) = bess(l, q1*rcomp)
   wm(1,2) = bess(l, q2*rcomp)
   wm(1,3) = 0.d0
   wm(1,4) = 0.d0
   wm(1,5) = 0.d0

!----- condition for integration of g_l(r) = 1
   DO I = 1, MESH
      WORK(I)  =  0.D0
   end do

   do I = 1, NRCUT
      WORK(I) = bess(l, q1*R(I)) * R(I)**(l+2)
      WORK(I) = WORK(I)*R(I) !integral R(I)
   end do
   CALL INTGBP( NRCUT, DX, WORK, wm(2,1) )

   do I = 1, NRCUT
      WORK(I) = bess(l, q2*R(I)) * R(I)**(l+2)
      WORK(I) = WORK(I)*R(I) !integral R(I)
   end do
   CALL INTGBP( NRCUT, DX, WORK, wm(2,2) )


   do I = 1, NRCUT
      if( r(i) >= rhcl ) then
         epower = gbeta*(RCL-r(i))
         if( epower < 90.d0 ) then
            sinr = sin(FQL*r(i))
            WORK(I) = sinr*sinr * exp(-epower)
         else
            WORK(I) = 0.d0
         end if
      else
         ri2 = r(i)*r(i)
         if( ln == 0 ) then
            riln = 1.d0
         else
            riln = r(i)**(ln)
         end if
         WORK(I) = riln &
         *( c_0 + ( c_2 + ( c_4 + c_6*ri2 )*ri2 )*ri2 )
      end if
      WORK(I) = WORK(I) * R(I)**(l+2)
      WORK(I) = WORK(I)*R(I) !integral R(I)
   end do
   CALL INTGBP( NRCUT, DX, WORK, wm(2,3) )

   do I = 1, NRCUT
      if( r(i) >= rhcl ) then
         epower = gbeta*(RCL-r(i))
         if( epower < 90.d0 ) then
            sinr = sin(FQL*r(i))
            WORK(I) = sinr*sinr*sinr * exp(-epower)
         else
            WORK(I) = 0.d0
         end if
      else
         ri2 = r(i)*r(i)
         if( ln == 0 ) then
            riln = 1.d0
         else
            riln = r(i)**(ln)
         end if
         WORK(I) = riln &
         *( d_0 + ( d_2 + ( d_4 + d_6*ri2 )*ri2 )*ri2 )
      end if
      WORK(I) = WORK(I) * R(I)**(l+2)
      WORK(I) = WORK(I)*R(I) !integral R(I)
   end do
   CALL INTGBP( NRCUT, DX, WORK, wm(2,4) )

   wm(2,5) = 1.d0

!----- condition for second derivative of g_l(rcomp) = 0
   if( abs(q1) < 1.d-10 ) then
      wm(3,1) = 0.d0
   else
      wm(3,1) = q1*q1*dbess(l,q1*rcomp,2)
   end if
   wm(3,2) = q2*q2*dbess(l,q2*rcomp,2)
   wm(3,3) = 2.d0*FQL2
   wm(3,4) = 0.d0
   wm(3,5) = 0.d0

!----- condition for third derivative of g_l(rcomp) = 0
   if( abs(q1) < 1.d-10 ) then
      wm(4,1) = 0.d0
   else
      wm(4,1) = q1*q1*q1*dbess(l,q1*rcomp,3)
   end if
   wm(4,2) = q2*q2*q2*dbess(l,q2*rcomp,3)
   wm(4,3) =  6.d0*FQL2*gbeta
   wm(4,4) = -6.d0*FQL3
   wm(4,5) = 0.d0

   CALL GSSJOR( WM, IP, 4, 4, 5, INDER )
   IF( INDER /= 0 ) return
   AAA1 = wm(1,5)
   AAA2 = wm(2,5)
   AAA3 = wm(3,5)
   AAA4 = wm(4,5)

   do I = 1, NRCUT
      GF(I) = AAA1*bess(l, q1*R(I)) + AAA2*bess(l, q2*R(I))
      if( r(i) >= rhcl ) then
         epower = gbeta*(RCL-r(i))
         if( epower < 90.d0 ) then
            sinr = sin(FQL*r(i))
            GF(I) = GF(I) + sinr*sinr * exp(-epower) &
            *( AAA3 + AAA4*sinr )
         end if
      else
         ri2 = r(i)*r(i)
         if( ln == 0 ) then
            riln = 1.d0
         else
            riln = r(i)**(ln)
         end if
         GF(I) = GF(I) + riln &
         *( AAA3*( c_0 + ( c_2 + ( c_4 + c_6*ri2 )*ri2 )*ri2 ) &
         + AAA4*( d_0 + ( d_2 + ( d_4 + d_6*ri2 )*ri2 )*ri2 ) )
      end if
   end do
   do I = NRCUT+1, mesh
      GF(I) = 0.d0
   end do

! heck
!      do I = 1, NRCUT
!         write(90,*) r(i), GF(I)
!      end do
!      do I = 1, NRCUT
!         WORK(I) = GF(i) * R(I)**(l+2)
!         WORK(I) = WORK(I)*R(I) !integral R(I)
!      end do
!      CALL INTGBP( NRCUT, DX, R, WORK, sum )
!      write(*,*)l, sum
!      if(l == 4)then
!      stop
!      end if

   RETURN
   END SUBROUTINE GFUNC



   subroutine ele1( r,rho, vnoe, iflag )
!-----------------------------------------------------------------------
!     Electronic density for Projector augmented wave method
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   dimension     r(msh), rho(msh)
   dimension     y(msh),f(msh) ,x(5,5)
!-----------------------------------------------------------------------
   do j = 1, nref
      do k = 1, nref
         do i = 1, icl+2
            y(i) = pusnl(i)*bas(j,i)*r(i)
            f(i) = pusnl(i)*bas(k,i)*r(i)
         end do
         call intgbp( icl+2, dx, y, t)
         call intgbp( icl+2, dx, f, s)

         x(j,k)=t*s

      end do
   end do

   do i = 1, mesh
      rho(i) = pusnl(i)*pusnl(i)
      ssg(i) = ssg(i) + vnoe * pusnl(i)*pusnl(i)
      do j = 1, nref
         do k = 1, nref
            if( iflag == 1 ) then
            ! ----- USPP
                rho(i) = rho(i) +      x(j,k)*qau(k,j,i)
                hsg(i) = hsg(i) + vnoe*x(j,k)*qau(k,j,i)
            else
            ! ----- PAW
                rho(i) = rho(i) +      x(j,k)*qau2(0,k,j,i)
                hsg(i) = hsg(i) + vnoe*x(j,k)*qau2(0,k,j,i)
            end if
         end do
      end do
   end do

   return
   end subroutine ele1


   subroutine ele2(r)
!-----------------------------------------------------------------------
!     Electronic density
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   dimension     r(msh)
   dimension     f(msh)
!-----------------------------------------------------------------------
   open(9,file='ele.dat')
   do i = 1, mesh
      write(9,*) r(i),sgv(i),sgc(i),rhous(i)
   end do
   close(9)
   write(*,'(4x,2a)') 'Create : ', 'ele.dat'

   do i = 1, mesh
      f(i) = sgv(i) * r(i)
   end do
   call intgbp (mesh, dx, f, ava)
   do i = 1, mesh
      f(i) = sgc(i) * r(i)
   end do
   call intgbp (mesh, dx, f, aco)
   do i = 1, mesh
      f(i) = rhous(i) * r(i)
   end do
   call intgbp (mesh, dx, f, uva)

   write(*,'(4x,a,3f15.10)') 'valence : ', ava, aco, uva

   return
   end subroutine ele2




