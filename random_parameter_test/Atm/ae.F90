


   SUBROUTINE ICDENS( LSG, SG, R, ZATM, XION, MESH, JOT )
!-----------------------------------------------------------------------
!         prepare starting charge density
!-----------------------------------------------------------------------
   IMPLICIT REAL*8 (A-H,O-Z)
   CHARACTER(2) ::   DUMMY
   DIMENSION  SG(*),R(*)

   logical :: LSG


   IF( LSG ) THEN
   !-------------  read charge density from data file  --------------------
      100 CONTINUE
      READ(JOT,1100,ERR=110,END=110) DUMMY
      IF( DUMMY /= '03' ) GO TO 100

      READ(JOT,1000,ERR=110,END=110) ( SG(J),J=1,MESH )
      REWIND JOT
      RETURN

      110 CONTINUE
      LSG = .FALSE.
      REWIND JOT
   ENDIF
   1000 FORMAT(4D18.10)
   1100 FORMAT(A2)

!------------   Thomas - Fermi approximation  --------------------------
   PAI    =  ACOS( -1.D0 )
   ALPHA  =  0.3058D0 * ( ZATM - XION )**( 1.D0 / 3.D0 )
   BETA   =  SQRT( 108.D0 / PAI ) * &
   MAX( 0.D0, ZATM - XION - 2.D0 ) * ALPHA
   GAMMA  =  SQRT( 4.D0 + 3.2D0 * XION )
   IF( ZATM - XION < 2.D0 ) THEN
      G2  =  0.5D0 * ( ZATM - XION ) * GAMMA**3
   ELSE
      G2  =  GAMMA**3
   ENDIF
   DO 35 I = 1, MESH
      X   =  ALPHA * R(I)
      XA  =  3.D0 * X
      IF( XA < 80.D0 ) THEN
         EXA  =  EXP( - XA )
      ELSE
         EXA  =  0.D0
      ENDIF
      XB  =  GAMMA * R(I)
      IF( XB < 80.D0 ) THEN
         EXB  =  EXP( - XB )
      ELSE
         EXB  =  0.D0
      ENDIF
      SG(I)  =  BETA * SQRT( X ) * EXA + G2 * R(I) * R(I) * EXB
   35 END DO


   RETURN
   END SUBROUTINE ICDENS



!   SUBROUTINE SCHROE( ZATM, N, L, EIG, R, VR, PNL, Y, S, G, &
!   KI, KJ, J06 )
SUBROUTINE SCHROE( N, L, EIG, R, VR, PNL, Y, S, G, &
&  KI, KJ, J06 )
!-----------------------------------------------------------------------
!    This routine solves radial Schroedinger equation for given (N,L)
!    and potential.
!          Y(I) --- SQRT( R ) * RNL( R ) = PNL( R ) / SQRT( R )
!          R    --- EXP( X )
!                        X: EQUAL MESH ( HX )
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   DIMENSION     R(*), VR(*), PNL(*)
   DIMENSION     Y(*), S(*), G(*)
   parameter     ( iii = 5  )
   parameter     ( jjj  = 10 )
   parameter     ( alp  = 1/274.074  )
   dimension     dvdx(mesh),abm(mesh)
   dimension     a(mesh),b(mesh),aa(jjj,jjj)
   dimension     bb(jjj),aaa(jjj),xx(mesh),t(mesh),f(mesh),d(mesh)
!-----------------------------------------------------------------------
   ! HH    =   DX * DX
   CLL1  =   DBLE( L * ( L + 1 ) )
   DELK  =   0.01D0
   EMIN  =  - ( ZATM / DBLE( N ) )**2
   EMAX  =  - 0.00001D0
   E     =   EIG
   IF( ( E > EMAX ) .OR. ( E < EMIN ) ) THEN
      E     =  0.5D0 * ( EMAX + EMIN )
   ENDIF
!--- coefficients for power series expansion of potential: V0, V1 ------
!         U1  =  ( VR(1) + 2.D0 * ZATM ) / R(1)
!         U2  =  ( VR(2) + 2.D0 * ZATM ) / R(2)
!         V0  =  ( U1 * R(2) - U2 * R(1) ) / ( R(2) - R(1) )
!         V1  =  ( U2 - U1 ) / ( R(2) - R(1) )
!--- (L+1)-th power of R -----------------------------------------------
!         R1L1  =  1.D0
!         R2L1  =  1.D0
!      DO 71 I = 1, L + 1
!         R1L1  =  R1L1 * R(1)
!         R2L1  =  R2L1 * R(2)
!   71 CONTINUE
!-----------------------------------------------------------------------
!           starts outward integration with a power series.
!-----------------------------------------------------------------------
!    1 CONTINUE
!      CA  =  -2.D0 * ZATM / DBLE( 2*L+2 )
!      CB  =  ( -2.D0 * ZATM * CA + ( -E + V0 ) ) / DBLE( 4*L+6 )
!      CC  =  ( -2.D0 * ZATM * CB + ( -E + V0 ) * CA + V1 )
!     &                                                  / DBLE( 6*L+12 )
!      PNL1 = ( 1.D0 + R(1) * ( CA + R(1) * ( CB + R(1) * CC ) ) ) * R1L1
!      PNL2 = ( 1.D0 + R(2) * ( CA + R(2) * ( CB + R(2) * CC ) ) ) * R2L1
!      Y(1)  =  PNL1 / SQRT( R(1) )
!      Y(2)  =  PNL2 / SQRT( R(2) )
   1 CONTINUE

   do i=1, mesh
      v(i)=vr(i)/r(i)
   enddo

!---- finite difference ------------------------------------------------
   call fd1 (1, mesh, mesh, v, dx, dvdx)
!-----------------------------------------------------------------------


   do i=1, mesh
      abm(i)=1+alp*alp*(e-v(i))                            !M(r)
      a(i)=alp*alp/abm(i)*dvdx(i)/r(i)                       !A(r)
      b(i)=-a(i)/r(i)-(l*(l+1)/r(i)/r(i)+abm(i)*(v(i)-e))    !B(r)
   enddo
   do i=1,5
      do j=1,5
         aa(i,j)=((l+j+1)*((l+j)+a(i)*r(i))+b(i)*r(i)*r(i)) &
         *r(i)**(l+j-1)
      end do
      bb(i)=-((l+1)*(l+a(i)*r(i))+b(i)*r(i)*r(i)) &
      *r(i)**(l-1)
   end do

   call gauss(aa,bb,aaa,iii,jjj)

   do i=1,5
      y(i) = r(i)**(l+0.5)                                         !Y
      xx(i) = (l+0.5)*r(i)**(l+0.5)                                !Y'
      do j=1,5
         y(i)=y(i)+aaa(j)*r(i)**(l+j+0.5)
         xx(i)=xx(i)+(l+j+0.5)*aaa(j)*r(i)**(l+j+0.5)
      enddo
      f(i)=-a(i)*r(i)*(xx(i)-y(i)*0.5)+((l+0.5d0)**2.d0 &
      +abm(i)*(v(i)-e)*r(i)*r(i))*y(i)                        !Y''
   enddo
!--- integrates outward to classical turning point ---------------------
!--- determine turning point ( KI ) ---
   DO 33 I = 3, MESH - 10
      KM  =  MESH - I
      IF( E * R(KM) - VR(KM) - CLL1 / R(KM) > 0.D0 ) GOTO 34
   33 END DO
   34 KI  =  KM + 1
!      WRITE(*,*) L,KAISU,E,KI
!      km=1359
!      WRITE(*,*) l, E,  R(KM),  VR(KM),  CLL1
!--- prepare G - function ---
!      DO 72 I = 1, KI + 2
!      G(I) = ( 0.25D0 - E * R(I) * R(I) + VR(I) * R(I) + CLL1 ) / 12.D0
   72 CONTINUE
!--- integartion by Noumerov scheme ---
!      DO 73 I = 3, KI + 2
!      Y(I) = ( 2.D0 * Y(I-1) - Y(I-2)
!     &       + HH * ( 10.D0 * G(I-1) * Y(I-1) + G(I-2) * Y(I-2) ) )
!     &       / ( 1.D0 - HH * G(I) )
   73 CONTINUE

!--- Adams-Bashforth Adams-Moulton -------------------------------------

   do i = 5, ki + 2
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
   enddo


!--- LOG - DER at KI ---
   I  =  KI
!      YP = ( Y(I-2) - 8.D0 * Y(I-1) + 8.D0 * Y(I+1) - Y(I+2) )
!     &     / ( 12.D0 * DX )
!      RLDM  =  ( YP / Y(I) - 0.5D0 ) / R(I)
   RLDM  =  ( xx(i) / Y(I) - 0.5D0 ) / R(I)
!--- count the number of nodes ---
   N1  =  0
   DO 35 K = 5, KI
      JK  =  INT( SIGN( 1.D0, Y(K)   ) )
      JK1 =  INT( SIGN( 1.D0, Y(K-1) ) )
      N1  =  N1 + ABS( JK - JK1 )
   35 END DO
   N1  =  1 + L + N1 / 2
!--- N1 = number of nodes + L + 1.  N1 should equal N ------------------
   IF( N - N1 ) 2, 4, 3
!--- too many nodes ----------------------------------------------------
   2 IF( E < EMAX ) THEN
      EMAX = E
   ENDIF
   E    =  0.5D0 * ( E + EMIN )
   DL1  =  ABS( EMAX - EMIN ) / ( ABS( E ) + 1.D0 )
   IF( DL1 > DELK) GOTO 1
   EMIN =  1.2D0 * EMIN
   WRITE(J06,630) N1, N, L, EMIN, E, EMAX
   630 FORMAT( ' too many nodes: N1, N, L =', 3I3, 2X, 3D16.8 )
   GOTO 1
!--- too few nodes -----------------------------------------------------
   3 IF( E > EMIN ) THEN
      EMIN  =  E
   ENDIF
   E     =  0.5D0 * ( E + EMAX )
   DL1   =  ABS( EMAX - EMIN ) / ( ABS( E ) + 1.D0 )
   IF( DL1 > DELK ) GOTO 1
!C         EMAX  =  MAX( 1.D0, 10.D0 / RMAX**2 )
   EMAX  =  1.D0
   E     =  0.5D0 * ( E + EMAX )
   WRITE(J06,631) N1, N, L, EMIN, E, EMAX
   631 FORMAT( ' too few  nodes: N1, N, L =', 3I3, 2X, 3D16.8 )
   GOTO 1
!--- correct number of nodes -------------------------------------------
   4 YKI  =  Y(KI)
!-----------------------------------------------------------------------
!                     starts inward integration
!                     outer boundary condition is set here.
!-----------------------------------------------------------------------
!--- determine the upper limit of the wave function ( KJ ) ---
   DO 21 I = KI, MESH
      KJ   =  I
   !C    IF( ( E * R(I) - VR(I) ) * R(I) + RINF ) 22, 21, 21
      IF( ( E * R(I) - VR(I) ) * R(I) + 2.D0 * RINF ) 22, 21, 21
   21 END DO
   22 CONTINUE
!--- prepare G - function ---
   DO 74 I = KI - 2, KJ
      G(I) = ( 0.25D0 - E * R(I) * R(I) + VR(I) * R(I) + CLL1 ) / 12.D0
   74 END DO
!--- asymptotic form of the wave function ---
   S(KJ)    =  1.D-20                                              !Y
!      S(KJ-1)  =  S(KJ) * SQRT( R(KJ) / R(KJ-1) )
!     &                  * EXP( ( R(KJ) - R(KJ-1) ) * SQRT( - E ) )

   t(kj)    =  -(0.5+sqrt(-e)*r(kj))*s(kj)                         !xx


   do i = kj-1, kj-4, -1
      s(i) = s(kj)*sqrt(r(kj)/r(i))*exp(sqrt(-e)*(r(kj)-r(i)))
      t(i) = -(0.5d0+sqrt(-e)*r(i))*s(i)
   end do

   do j = 1, 5

      do i = kj, kj-4,-1
         d(i)=-a(i)*r(i)*(t(i)-s(i)*0.5d0)+((l+0.5d0)**2.d0         & !f
         +abm(i)*(v(i)-e)*r(i)*r(i))*s(i)
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

!--- integration by Noumerov scheme ---
!      DO 75 J = KI - 2, KJ - 2
!      I    =  KJ - 2 - ( J - KI + 2 )
!      S(I) = ( 2.D0 * S(I+1) - S(I+2)
!     &       + HH * ( 10.D0 * G(I+1) * S(I+1) + G(I+2) * S(I+2) ) )
!     &       / ( 1.D0 - HH * G(I) )
   75 CONTINUE
!--- Adams-Bashforth Adams-Moulton -------------------------------------

   do i = kj-4, ki - 1, -1
      s(i-1)=s(i)-dx/720.d0*(1901.d0*t(i)-2774.d0*t(i+1) &
      +2616.d0*t(i+2)-1274.d0*t(i+3)+251.d0*t(i+4))
      t(i-1)=t(i)-dx/720.d0*(1901.d0*d(i)-2774.d0*d(i+1) &
      +2616.d0*d(i+2)-1274.d0*d(i+3)+251.d0*d(i+4))
      d(i-1)=-a(i-1)*r(i-1)*(t(i-1)-s(i-1)*0.5d0)+((l+0.5d0)**2.d0 &
      +abm(i-1)*(v(i-1)-e)*r(i-1)*r(i-1))*s(i-1)

      s(i-1)=s(i)-dx/720.d0*(251.d0*t(i-1)+646.d0*t(i) &
      -264.d0*t(i+1)+106.d0*t(i+2)-19.d0*t(i+3))
      t(i-1)=t(i)-dx/720.d0*(251.d0*d(i-1)+646.d0*d(i) &
      -264.d0*d(i+1)+106.d0*d(i+2)-19.d0*d(i+3))
      d(i-1)=-a(i-1)*r(i-1)*(t(i-1)-s(i-1)*0.5d0)+((l+0.5d0)**2.d0 &
      +abm(i-1)*(v(i-1)-e)*r(i-1)*r(i-1))*s(i-1)
   enddo

!--- LOG - DER at KI ---
   I  =  KI
!      SP = ( S(I-2) - 8.D0 * S(I-1) + 8.D0 * S(I+1) - S(I+2) )
!     &     / ( 12.D0 * DX )
!      RLDP  =  ( SP / S(I) - 0.5D0 ) / R(I)
   RLDP  =  ( t(i) / S(I) - 0.5D0 ) / R(I)
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
   W  =  0.D0
   DO 9 I = 2, MESH - 1, 2
      W  =  W + 2.D0 * ( Y(I) * R(I) )**2 + ( Y(I+1) * R(I+1) )**2
   9 END DO
   W  =  ( Y(1) * R(1) )**2 - ( Y(MESH) * R(MESH) )**2 + 2.D0 * W
   W  =  ( Y(1) * R(1) )**2 / DBLE( 2*L+3 ) + W * DX / 3.D0
!--- energy correction due to diff. of LOG - DER ---
   DE  =  R(KI) * Y(KI) * Y(KI) * ( RLDM - RLDP ) / W

   DL1 =  ABS( EMAX - EMIN ) / ( ABS( E ) + 1.D0 )
   DL  =  ABS( DE / E )
   IF( ( DL > DELK ) .AND. ( DL1 < DELK ) ) GOTO 10
   IF( DL1 < DELL ) GOTO 10
   IF( DE  > 0.D0 ) THEN
      EMIN  =  E
   ELSE
      EMAX  =  E
   ENDIF
   DEP   =  DE
   E  =  E + DEP
   IF( DL  > DELL ) THEN
      IF( E > EMAX ) THEN
         E  =  0.5D0 * ( E - DEP + EMAX )
      ENDIF
      IF( E < EMIN ) THEN
         E  =  0.5D0 * ( E - DEP + EMIN )
      ENDIF
      GOTO 1
   ENDIF
!--- eigenvalue is conserved -------------------------------------------
   10 CR    =  1.D0 / SQRT( W )
   DO 11 I = 1, MESH
      PNL(I)  =  CR * Y(I) * SQRT( R(I) )
   11 END DO
   EIG  =  E

!      if (lus.eq.l .and. nus .eq. n) then
!      ofile=(ms(l+1)//'_Pae.dat')
!      open(9,file=ofile)
!         do i=1, mesh
!!            vus(i)=v(i)
!!            p0ae(i)=pnl(i)
!            write(9,*) r(i),pnl(i)
!         end do
!      end if
!      close(9)
   RETURN
   END SUBROUTINE SCHROE






   SUBROUTINE KINE( L, PNL, R, EKIN, A )
!-----------------------------------------------------------------------
!    Kinetic Energy
!  input :
!         PNL  ...... radial wave function
!           A  ...... work area
!  output :
!       EKIN   ......  Kinetic energy
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   DIMENSION     R(*), PNL(*)
   DIMENSION     A(*)

!--- EKIN : kinetic E. ----
   CLL1  =   DBLE( L * ( L + 1 ) )
   GI = ( PNL(2) - PNL(1) ) / ( R(2) - R(1) )
   A(1) = GI*GI*R(1) + CLL1*PNL(1)*PNL(1)/R(1)
   GI = ( PNL(3) - PNL(1) ) / ( R(3) - R(1) )
   A(2) = GI*GI*R(2) + CLL1*PNL(2)*PNL(2)/R(2)
   DO I = 3, MESH - 2
      GI = ( PNL(I-2) - 8.D0*PNL(I-1) + 8.D0*PNL(I+1) - PNL(I+2) ) &
      / ( 12.D0 * DX )

      A(I) = ( GI*GI + CLL1*PNL(I)*PNL(I) )/R(I)
   ENDDO
   A(MESH-1) = 0.D0
   A(MESH)   = 0.D0
   CALL INTGBP( MESH, DX, A, EKIN )

   RETURN
   END SUBROUTINE KINE


subroutine enegyf( vext, vhar, vexc, vcor, eexc, ecor, sg, a, &
sume, epart )
!-----------------------------------------------------------------------
!    Energy
!  input :
!      VEXT(I) ...... external potential
!      VHAR(I) ...... Hartree potential
!      VEXC(I) ...... Exchange potential
!      VCOR(I) ...... Correlation potential
!      EEXC(I) ...... Exchange energy
!      ECOR(I) ...... Correlation energy
!        SG(I) ...... charge density
!           A  ...... work area
!  output :
!       SUME   ...... total energy
!       EPART  ...... total-energy parts
!            1 ......  Kinetic energy
!            2 ......  External energy
!            3 ......  Hartree energy
!            4 ......  Exchange energy
!            5 ......  Correlation energy
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   DIMENSION     SG(*)
   DIMENSION     A(*)
   DIMENSION     EPART(*)
   DIMENSION     VEXT(*), VHAR(*), VEXC(*), VCOR(*), EEXC(*), ECOR(*)


   ! ekin = epart(1)
!--- eext : external potential e. ----
   do i = 1, mesh
      a(i) = vext(i)*sg(i)
   enddo
   call intgbp( mesh, dx, a, eext )

!--- ehar : hartree potential e. ----
   do i = 1, mesh
      a(i) = vhar(i)*sg(i)
   enddo
   call intgbp( mesh, dx, a, ehar )
   ehar = 0.5d0*ehar

!--- eexcv : exchange potential e.    ----
!--- eexce : exchange e.              ----
!--- ecorv : correlation potential e. ----
!--- ecore : correlation e.           ----
   do i = 1, mesh
      a(i) = vexc(i)*sg(i)
   enddo
   call intgbp( mesh, dx, a, eexcv )
   do i = 1, mesh
      a(i) = eexc(i)*sg(i)
   enddo
   call intgbp( mesh, dx, a, eexce )
   do i = 1, mesh
      a(i) = vcor(i)*sg(i)
   enddo
   call intgbp( mesh, dx, a, ecorv )
   do i = 1, mesh
      a(i) = ecor(i)*sg(i)
   enddo
   call intgbp( mesh, dx, a, ecore )

   sume = sume - ehar + eexce + ecore - ( eexcv + ecorv )
   epart(2) = eext
   epart(3) = ehar
   epart(4) = eexce
   epart(5) = ecore

   return
end subroutine enegyf



