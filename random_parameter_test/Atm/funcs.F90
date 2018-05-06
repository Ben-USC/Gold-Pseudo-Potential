   subroutine gauss(a,b,x,n,nn)
!----------------------------------------------------------------------

!----------------------------------------------------------------------
   implicit real*8 (a-h,o-z)
   real*8 :: a(nn,nn), b(nn), x(nn), p, q
   integer :: i,j,k,n,nn,d

   do k=1, n-1

      d = k
      g=abs(a(k,k))
      do i = k+1, n
         h = abs(a(i,k))
         if(h > g) then
            d = i
            g = h
         end if
      end do

      do j = k, n
         w = a(k,j)
         a(k,j) = a(d,j)
         a(d,j) = w
      end do
      w=b(k)
      b(k)=b(d)
      b(d)=w

      p=a(k,k)
      do i = k+1, n
         q = a(i,k)/p
         do j = k, n
            a(i,j)=a(i,j)-q * a(k,j)
         end do
         b(i) = b(i) - q * b(k)
      end do
   end do

   x(n) =b(n) / a(n,n)

   do k = n-1, 1, -1
      s = b(k)
      do j = k + 1, n
         s = s - a(k,j) * x(j)
      end do
      x(k) = s / a(k,k)
   end do
   return
   end subroutine gauss




   SUBROUTINE GSSJOR( W, IP, N, MM, MM1, INDER )
!-----------------------------------------------------------------------
!        <<<   Gauss-Jurdan elemination method >>>   since  1990/09/07
!-----------------------------------------------------------------------
!  ( input )
!     N        ......  number of variables
!     W(i,j)   ......  coefficients
!     IP       ......  work area
!  ( output )
!     W(i,N1)  ......  solution x_i
!-----------------------------------------------------------------------
   IMPLICIT REAL*8 ( A-H, O-Z )
   DIMENSION  W(MM,MM1),  IP(MM)
   DATA DZERO / 1.0D-30 /

   N1 = N + 1

   DO 10 I = 1, N
      IP(I) = I
   10 END DO

   DO 20 I = 1, N
      PIVOT = 0.0
      DO 21 J = 1, N1
         PIVOT = MAX( PIVOT, ABS(W(I,J)) )
      21 END DO
      IF( PIVOT < DZERO ) THEN
         INDER = 5
         RETURN
      ENDIF
      DO 22 J = 1, N1
         W(I,J) = W(I,J)/PIVOT
      22 END DO
   20 END DO


!--- reducing the matrix to triangular form ----------------------------
   DO 1000 I = 1, N - 1

      PIVOT = 0.0
      DO 1100 J = I, N
         JJ = IP(J)
         IF( ABS(W(JJ,I)) > PIVOT ) THEN
            PIVOT = ABS(W(JJ,I))
            JPV   = J
         ENDIF
      1100 END DO
      IF( PIVOT < DZERO ) THEN
         INDER = 6
         RETURN
      ENDIF
      IF( JPV /= I ) THEN
         IPOL    = IP(JPV)
         IP(JPV) = IP(I)
         IP(I)   = IPOL
      ENDIF

      IPV = IP(I)
      DO 1200 J = I + 1, N
         JJ  = IP(J)
         DMJ = W(JJ,I)/W(IPV,I)
      !            W(JJ,I) = DMJ
         W(JJ,I) = 0.0
         DO 1200 J2 = I + 1, N1
            W(JJ,J2) = W(JJ,J2) - DMJ*W(IPV,J2)
      1200 END DO

   1000 END DO



!--- solutions ---------------------------------------------------------
   DO 1500 I = 1, N
      IV = IP(I)
      DIEL = W(IV,I)
      W(IV,I) = 1.0D0
      DO 1500 J = I + 1, N1
         W(IV,J) = W(IV,J)/DIEL
   1500 END DO

   IPN = IP(N)
   DO 1550 I = 1, N
      IV = IP(I)
      W(IPN,I) = W(IV,N1)
   1550 END DO
   DO 1600 I = N - 1, 1, -1
      IV = IP(I)
      DO 1600 J = I + 1, N
         W(IPN,I) = W(IPN,I) - W(IV,J)*W(IPN,J)
   1600 END DO
   DO 1650 I = 1, N
      W(I,N1) = W(IPN,I)
   1650 END DO

!--- value for the normal return
   INDER = 0


   RETURN
   END SUBROUTINE GSSJOR


   subroutine fd1 (min, max, n, f, dx, dfdx)
! --- ------------------------------------------------------------------
!     finite difference  first derivative
!-----------------------------------------------------------------------
   implicit real*8 (A-H,O-Z)
   dimension      f(n), dfdx(n)

   idf = 0          ! flag dfdx

   do i=min, max

      idf=idf+1

      if (idf == min) then
         dfdx(i) = (-147*f(i)+360*f(i+1)-450*f(i+2) &
         +400*f(i+3)-225*f(i+4) &
         +72*f(i+5)-10*f(i+6))/60/dx
      else if(idf == min+1) then
         dfdx(i) = (-10*f(i-1)-77*f(i) &
         +150*f(i+1)-100*f(i+2) &
         +50*f(i+3)-15*f(i+4)+2*f(i+5))/60/dx
      else if(idf == min+2) then
         dfdx(i) = (2*f(i-2)-24*f(i-1) &
         -35*f(i)+80*f(i+1) &
         -30*f(i+2)+8*f(i+3)-1*f(i+4))/60/dx
      else if (idf == max-2) then
         dfdx(i) = (1*f(i-4)-8*f(i-3) &
         +30*f(i-2)-80*f(i-1) &
         +35*f(i)+24*f(i+1)-2*f(i+2))/60/dx
      else if(idf == max-1) then
         dfdx(i) = (-2*f(i-5)+15*f(i-4) &
         -50*f(i-3)+100*f(i-2) &
         -150*f(i-1)+77*f(i)+10*f(i+1))/60/dx
      else if(idf == max) then
         dfdx(i) = (10*f(i-6)-72*f(i-5) &
         +225*f(i-4)-400*f(i-3) &
         +450*f(i-2)-360*f(i-1)+147*f(i))/60/dx
      else
         dfdx(i) = (-1*f(i-3)+9*f(i-2) &
         -45*f(i-1)+0*f(i) &
         +45*f(i+1)-9*f(i+2)+1*f(i+3))/60/dx
      end if

   end do

   return
   end subroutine fd1

   subroutine fd2 (min, max, n, f, dx, dfdx)
! --- ------------------------------------------------------------------
!     finite difference   second derivative
!-----------------------------------------------------------------------
   implicit real*8 (A-H,O-Z)
   dimension      f(n), dfdx(n)

   idf = 0          ! flag dfdx

   do i=min, max

      idf=idf+1

      if (idf == min) then
         dfdx(i) = (812*f(i)-3132*f(i+1)+5265*f(i+2) &
         -5080*f(i+3)+2970*f(i+4) &
         -972*f(i+5)+137*f(i+6))/180/dx/dx
      else if(idf == min+1) then
         dfdx(i) = (137*f(i-1)-147*f(i) &
         -255*f(i+1)+470*f(i+2) &
         -285*f(i+3)+93*f(i+4)-13*f(i+5))/180/dx/dx
      else if(idf == min+2) then
         dfdx(i) = (-13*f(i-2)+228*f(i-1) &
         -420*f(i)+200*f(i+1) &
         +15*f(i+2)-12*f(i+3)+2*f(i+4))/180/dx/dx
      else if (idf == max-2) then
         dfdx(i) = (2*f(i-4)-12*f(i-3) &
         +15*f(i-2)+200*f(i-1) &
         -420*f(i)+228*f(i+1)-13*f(i+2))/180/dx/dx
      else if(idf == max-1) then
         dfdx(i) = (-13*f(i-5)+93*f(i-4) &
         -285*f(i-3)+470*f(i-2) &
         -255*f(i-1)-147*f(i)+137*f(i+1))/180/dx/dx
      else if(idf == max) then
         dfdx(i) = (137*f(i-6)-972*f(i-5) &
         +2970*f(i-4)-5080*f(i-3) &
         +5265*f(i-2)-3132*f(i-1)+812*f(i))/180/dx/dx
      else
         dfdx(i) = (2*f(i-3)-27*f(i-2) &
         +270*f(i-1)-490*f(i) &
         +270*f(i+1)-27*f(i+2)+2*f(i+3))/180/dx/dx
      end if

   end do

   return
   end subroutine fd2

   subroutine fd3n (min, max, n, f, dx, dfdx)
! --- ------------------------------------------------------------------
!     finite difference   third derivative
!-----------------------------------------------------------------------
   implicit real*8 (A-H,O-Z)
   dimension      f(n), dfdx(n)

   do i=min+6, max-6
   !            dfdx(i) = (-479*f(i-6)+6840*f(i-5)
   !     &                 -46296*f(i-4)+198760*f(i-3)
   !     &                 -603315*f(i-2)+764208*f(i-1)+0*f(i)
   !     &                 -764208*f(i+1)+603315*f(i+2)
   !     &                 -198760*f(i+3)+46296*f(i+4)-6840*f(i+5)
   !     &                 +479*f(i+6))/302400.d0/dx/dx/dx
      dfdx(i) = ( 1*f(i-3) &
      -8*f(i-2)+13*f(i-1)+0*f(i) &
      -13*f(i+1)+8*f(i+2) &
      -1*f(i+3))/8.d0/dx/dx/dx
   end do

   return
   end subroutine fd3n

   subroutine fd3 (min, max, n, f, dx, dfdx)
! --- ------------------------------------------------------------------
!     finite difference   third derivative
!-----------------------------------------------------------------------
   implicit real*8 (A-H,O-Z)
   dimension      f(n), dfdx(n)

   do i=min+6, max-6
      dfdx(i) = (-479*f(i-6)+6840*f(i-5) &
      -46296*f(i-4)+198760*f(i-3) &
      -603315*f(i-2)+764208*f(i-1)+0*f(i) &
      -764208*f(i+1)+603315*f(i+2) &
      -198760*f(i+3)+46296*f(i+4)-6840*f(i+5) &
      +479*f(i+6))/302400.d0/dx/dx/dx
   end do

   return
   end subroutine fd3


   subroutine fd4 (min, max, n, f, dx, dfdx)
! --- ------------------------------------------------------------------
!     finite difference   fourth derivative
!-----------------------------------------------------------------------
   implicit real*8 (A-H,O-Z)
   dimension      f(n), dfdx(n)

   do i=min+6, max-6
      dfdx(i) = (479*f(i-6)-8208*f(i-5) &
      +69444*f(i-4)-397520*f(i-3) &
      +1809945*f(i-2)-4585248*f(i-1)+6222216*f(i) &
      -4585248*f(i+1)+1809945*f(i+2) &
      -397520*f(i+3)+69444*f(i+4)-8208*f(i+5) &
      +479*f(i+6))/453600.d0/dx/dx/dx/dx
   end do

   return
   end subroutine fd4

   subroutine fd4n (min, max, n, f, dx, dfdx)
! --- ------------------------------------------------------------------
!     finite difference   fourth derivative
!-----------------------------------------------------------------------
   implicit real*8 (A-H,O-Z)
   dimension      f(n), dfdx(n)

   do i=min+6, max-6
      dfdx(i) = (-1*f(i-3) &
      +12*f(i-2)-39*f(i-1)+56*f(i) &
      -39*f(i+1)+12*f(i+2) &
      -1*f(i+3)) /6.d0/dx/dx/dx/dx
   end do

   return
   end subroutine fd4n


module param_for_Bessel
!-----------------------------------------------------------------------
! type declaration of coefficients for Bessel functions
!-----------------------------------------------------------------------
   implicit none

   real(8) :: FJ02, FJ04, FJ06, FJ08, FJ010, FJ012
   real(8) :: FJ11, FJ13, FJ15, FJ17, FJ19,  FJ111
   real(8) :: FJ22, FJ24, FJ26, FJ28, FJ210, FJ212
   real(8) :: FJ33, FJ35, FJ37, FJ39, FJ311
   save

end module

   SUBROUTINE CNSTTS
!-----------------------------------------------------------------------
!     some constants
!-----------------------------------------------------------------------
   use param_for_Bessel
   PARAMETER( JISP = 15 )
   IMPLICIT REAL*8 (A-H,O-Z)
   dimension    CIF(0:JISP)

   CN = 0.0
   FN = 1.D0
   CIF(0) = 1.D0
   DO L = 1, JISP
      CN = CN + 1.D0
      FN = FN*CN
      CIF(L) = 1.D0/FN
   END DO
!-----------------------------------------------------------------------
!             coefficients for Bessel functions

   ! FCC2  = - CIF(2)
   ! FCC4  =   CIF(4)
   ! FCC6  = - CIF(6)
   ! FCC8  =   CIF(8)
   ! FCC10 = - CIF(10)
   ! FCC12 =   CIF(12)

   FJ02  = - CIF(3)
   FJ04  =   CIF(5)
   FJ06  = - CIF(7)
   FJ08  =   CIF(9)
   FJ010 = - CIF(11)
   FJ012 =   CIF(13)

   FJ11  =   2.D0*CIF(3)
   FJ13  = - 4.D0*CIF(5)
   FJ15  =   6.D0*CIF(7)
   FJ17  = - 8.D0*CIF(9)
   FJ19  =  10.D0*CIF(11)
   FJ111 = -12.D0*CIF(13)

   FJ22  =   2.D0* 4.D0*CIF(5)
   FJ24  = - 4.D0* 6.D0*CIF(7)
   FJ26  =   6.D0* 8.D0*CIF(9)
   FJ28  = - 8.D0*10.D0*CIF(11)
   FJ210 =  10.D0*12.D0*CIF(13)
!     FJ212 = -12.D0*14.D0*CIF(15)

   FJ33  =   2.D0* 4.D0* 6.D0*CIF(7)
   FJ35  = - 4.D0* 6.D0* 8.D0*CIF(9)
   FJ37  =   6.D0* 8.D0*10.D0*CIF(11)
   FJ39  = - 8.D0*10.D0*12.D0*CIF(13)
   FJ311 =  10.D0*12.D0*14.D0*CIF(15)

   RETURN
   END SUBROUTINE CNSTTS


   DOUBLE PRECISION FUNCTION BSL0( Q )
   use param_for_Bessel
   IMPLICIT REAL*8 (A-H,O-Z)

   IF( ABS(Q) < 7.D-01 ) THEN
      XF12 = Q*Q
      BSL0 = 1.0D0 + XF12*( FJ02 + XF12*( FJ04 + XF12*( FJ06 &
      + XF12*( FJ08 + XF12*FJ010 ) ) ) )
   ELSE
      BSL0 = SIN(Q)/Q
   ENDIF

   RETURN
   END

   DOUBLE PRECISION FUNCTION BSL1( Q )
   use param_for_Bessel
   IMPLICIT REAL*8 (A-H,O-Z)

   IF( ABS(Q) < 7.D-01 ) THEN
      XF12 = Q*Q
      BSL1 = Q*( FJ11 + XF12*( FJ13 + XF12*( FJ15 &
      + XF12*( FJ17 + XF12*FJ19 ) ) ) )
   ELSE
      BSL1 = ( SIN(Q)/Q - COS(Q) )/Q
   ENDIF

   RETURN
   END

   DOUBLE PRECISION FUNCTION BSL2( Q )
   use param_for_Bessel
   IMPLICIT REAL*8 (A-H,O-Z)

   IF( ABS(Q) < 7.D-01 ) THEN
      XF12 = Q*Q
      BSL2 = XF12*( FJ22 + XF12*( FJ24 + XF12*( FJ26 &
      + XF12*( FJ28 + XF12*FJ210 ) ) ) )
   ELSE
      BS0  = SIN(Q)/Q
      BS1  = ( BS0 - COS(Q) )/Q
      BSL2 = 3.D0*BS1/Q - BS0
   ENDIF

   RETURN
   END


   DOUBLE PRECISION FUNCTION BSL3( Q )
   use param_for_Bessel
   IMPLICIT REAL*8 (A-H,O-Z)

   IF( ABS(Q) < 7.D-01 ) THEN
      XF12 = Q*Q
      BSL3 = Q*XF12*( FJ33 + XF12*( FJ35 + XF12*( FJ37 &
      + XF12*( FJ39 + XF12*FJ311 ) ) ) )
   ELSE
      BS0  = SIN(Q)/Q
      BS1  = ( BS0 - COS(Q) )/Q
      BS2  = 3.D0*BS1/Q - BS0
      BSL3 = 5.D0*BS2/Q - BS1
   ENDIF

   RETURN
   END

!       DOUBLE PRECISION FUNCTION sbessl( q, l )
   DOUBLE PRECISION FUNCTION bess( l, q )
!-----------------------------------------------------------------------
!   spherical Bessel functions
!-----------------------------------------------------------------------
   IMPLICIT REAL*8 (A-H,O-Z)
   logical :: lflag
   save lflag
   data lflag / .TRUE. /

   if( lflag ) call cnstts
   lflag = .FALSE.

   if( l == 0 ) then
      bess = BSL0( q )
   else if( l == 1 ) then
      bess = BSL1( q )
   else if( l == 2 ) then
      bess = BSL2( q )
   else if( l == 3 ) then
      bess = BSL3( q )
   else if( l == -1 ) then
      if( q < 1.d-15 ) then
         bess = 0.d0
      else
         bess = cos(q)/q
      endif
   else if( l > 3 ) then
      if( q < 1.d-15 ) then
         bess = 0.d0
      else
         sbm2 = BSL2( q )
         sbm1 = BSL3( q )
         do i = 4, l
            bess = (2.d0*dble(i)-1.d0)*sbm1/q - sbm2
            sbm2 = sbm1
            sbm1 = bess
         enddo
      endif
   else
      write(*,*) ' error in bess :',l
      bess = 0.d0
   !-----Finalize the parallel environment
   !         call end_parallel(ierr)
   !         stop
   endif

   return
   END



real(8) function dbess (l,q,ld)
!-----------------------------------------------------------------------
!   derivatives of spherical Bessel functions
!-----------------------------------------------------------------------
   implicit real*8 (a-h,o-z)

   if( abs(q).lt.1.d-10 ) then
       write(*,*) ' error in dbess: not supported for q.le.0 yet.',q
       stop
   endif

   dbess = bess( l, q )
   if( ld.eq.0 ) return

   b00 = dbess
   b01 = bess( l+1, q )
   dbess = - b01 + dble(l)*b00/q
   if( ld.eq.1 ) return

   b10 = dbess
   b02 = bess( l+2, q )
   b11 = - b02 + dble(l+1)*b01/q
   dbess = - b11 + dble(l)*( b10 - b00/q )/q
   if( ld.eq.2 ) return

   b20 = dbess
   b03 = bess( l+3, q )
   b12 = - b03 + dble(l+2)*b02/q
   b21 = - b12 + dble(l+1)*( b11 - b01/q )/q
   dbess = - b21 + dble(l)*( b20 - 2.d0*( b10 - b00/q )/q )/q
   if( ld.eq.3 ) return

   b30 = dbess
   b04 = bess( l+4, q )
   b13 = - b04 + dble(l+3)*b03/q
   b22 = - b13 + dble(l+2)*( b12 - b02/q )/q
   b31 = - b22 + dble(l+1)*( b21 - 2.d0*( b11 - b01/q )/q )/q
   dbess = - b31  &
   &    + dble(l)*( b30 - 3.d0*( b20 - 2.d0*( b10 - b00/q )/q )/q )/q
   if( ld.eq.4 ) return


   write(*,*) ' error in dbess: not supported for ld.gt.4 yet.',ld
   stop
end function dbess


!    SUBROUTINE INTGBP( MSH, DX, R, BL, CL )
   SUBROUTINE INTGBP( MSH, DX, BL, CL )
!-----------------------------------------------------------------------
!                                                          1992/12/08
!  Simpson's quadrature by using log-scale mesh
!    Caution !!   The integrand BL should be multiplied by R(I)

!    output :  CL ... integration of BL
!-----------------------------------------------------------------------
   IMPLICIT REAL*8 (A-H,O-Z)
   DIMENSION      BL(*)
   DATA DZERO / 1.0D-20 /

   MESH = MSH
   IF( MOD(MSH,2) == 0 ) MESH = MESH - 1
   CL  =  0.D0
   DO 10 I = 2, MESH - 1, 2
      CL  =  CL + 2.D0*BL(I) + BL(I+1)
   10 END DO
   IX = MESH
   CL  =  BL(1) - BL(IX) + 2.D0*CL
   CL  =  CL*DX/3.D0
   IF( MOD(MSH,2) == 0 ) THEN
      CL  =  CL + ( BL(MSH) + BL(IX) ) * DX * .5D0
   ENDIF
!      IF( ABS(BL(1)).GT.DZERO ) THEN
   IF( ABS(BL(1)) > DZERO .AND. BL(1)*BL(2) > 0.d0 ) THEN
      S2 = LOG( BL(2)/BL(1) )/DX
      CL = CL + BL(1)/S2
   ENDIF


   RETURN
   END SUBROUTINE INTGBP


