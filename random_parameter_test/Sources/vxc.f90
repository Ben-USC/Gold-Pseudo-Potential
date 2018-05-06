
SUBROUTINE VXC( RS, EX, VX, EC, VC )
!-----------------------------------------------------------------------
!     one - electron potential
!     energy units are a.u. ( = Ryd * 0.5 )
!-----------------------------------------------------------------------
   IMPLICIT REAL*8 (A-H,O-Z)
!S-- Slater exchange ---------------------------------------------------
!S    VXC   = -3.D0 * 0.6108870575D0 / RS * 0.5
!-----------------------------------------------------------------------
!--- Gunarsson-Lundqvist exchange correlation potential. ---------------
!      BETA  =  1.D0 + 0.0545D0 * RS * LOG( 1.D0 + 11.4D0/RS )
!      VXC   = - 0.6108870575D0*BETA / RS
!-----------------------------------------------------------------------
   DATA  A1 / -0.4582D0 /
   DATA  B1,B2,B3 / -0.1423D0, 1.0529, 0.3334 /
   DATA  C1,C2,C3,C4 / -0.0480D0, 0.0311, -0.0116, 0.0020 /
!      DATA  BB / 0.0140D0 /
!--- Ceperly and Alder as parametrized by Perdew and Zunger ------------
!       EX    ...... exchange energy
!       VX    ...... exchange potential ( non-relativistic )
!       RCORR ...... relativistic correction fot exchange potential
!       EC    ...... correlation energy
!       VC    ...... correlation potential
!-----------------------------------------------------------------------
   EX    =  A1/RS
   VX    =  4.0D0*EX/3.0D0
!C      BETA  =  BB/RS
!C      YETA  =  SQRT( 1.0D0 + BETA*BETA )
!C      ZETA  =  LOG( BETA + YETA )
!C      RCORR =  -0.5D0 + 1.5D0*ZETA/BETA/YETA
!C      VXR   =  VX * RCORR
   IF( RS >= 1.0D0 ) THEN
      BUNB = 1.0D0 + B2*SQRT(RS) + B3*RS
      EC   = B1/BUNB
      VC   = EC*(1.0D0 + (B2*SQRT(RS)/6.0D0 + B3*RS/3.0D0)/BUNB )
   ELSE
      EC = C1 + C2*LOG(RS) + C3*RS + C4*RS*LOG(RS)
      VC = EC - (C2 + (C3 + C4)*RS + C4*RS*LOG(RS))/3.0D0
   ENDIF
!C      VXC   = VXR + VC
!      VXC   = VX + VC
!-----------------------------------------------------------------------

   RETURN
END SUBROUTINE VXC


subroutine pbe( rho, drho, d2rho, ri, ex, vx, ec, vc, iflag )
!-----------------------------------------------------------------------
! Generalized gradient corrected exchange correlation potential by PBE
!     for spin unpolarized case
!              ^^^^^^^^^^^
! input
!     rho    : density
!     drho   : derivative of density          d rho/ dr
!     d2rho  : second derivative of density  d^2rho/ dr^2
!     iflag  : 1 -> relativistic correction / else -> not
!
!       RCORRE...... relativistic correction for exchange energy
!       RCORR ...... relativistic correction for exchange potential
!
! output
!     ex     : exchange energy
!     vx     : exchange potential
!     ec     : correlation energy
!     vc     : correlation potential
!
!   ***  energy units are hartree ( = Ryd * 0.5 )
!-----------------------------------------------------------------------
   implicit real*8 (a-h,o-z)
   parameter( c1b3 = 1.d0/3.d0 )
   parameter( c4b3 = 4.d0/3.d0 )
   parameter( c7b3 = 7.d0/3.d0 )
   save rs0, ax, bx, um, uk, &
   a0, alp01, bet01, bet02, bet03, bet04, p0, &
   betapw, BB, t0, h1

!--- parameter for relativistic correction
   DATA  BB / 0.0140D0 /

   data rs0 / 0.62035049090d0 /
!--- parameter for exchange energy
   data ax / -0.738558766382022405884230032680836d0 /
   data bx /  0.16162045967d0 /
   data um / 0.2195149727645171d0 /
   data uk / 0.8040d0 /

!--- parameter for correlation energy
   data a0, alp01, bet01, bet02, bet03, bet04, p0 &
   / 0.0621814d0, 0.21370d0, &
   &      7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0, 1.d0 /
   data betapw / 0.06672455060314922d0 /
   data t0 / 0.2519289703424d0 /
   data h1 / 0.03109069086965489503494086371273d0 /


   if( rho < 1.d-30 ) then
      ex = 0.d0
      vx = 0.d0
      ec = 0.d0
      vc = 0.d0
      return
   endif

   rho3 = rho**c1b3
   rs   = rs0/rho3


   adrho = abs(drho)
   adrbro = adrho / rho

!--- exchange terms ----------------------------------------------------
   ex0  = ax * rho3

   s  = bx * adrbro / rho3
   if( s < 1.d-15 ) then
      ex = ex0
      vx = c4b3*ex0
   else
      s2 = s*s

      ul = um/uk
      fsbb = 1.d0 + s2*ul
      fsbb = 1.d0/fsbb
      fst  = - uk*fsbb
      fs   = 1.d0 + uk + fst

   !     --- exchange energy ---
      ex = ex0*fs


      dfsds  = - 2.d0*ul*s * fst*fsbb

      dfxdn  = c4b3*ex0 * ( fs - s*dfsds )
      dfxddn = ax * bx * dfsds

   !=== from here, for olny sherical symmetric case ===
      d2fsds = - 2.d0*ul*( fst + 2.d0*dfsds*s )*fsbb

      drbro  = drho / rho
      d2rbro = d2rho/drho
      dsdr = s*( d2rbro - c4b3 * drbro )

      dfxdr = ax*bx*d2fsds*dsdr

   !     --- exchange potential ---
      vx = dfxdn - ( dfxdr + 2.d0*dfxddn/ri )*sign(1.d0,drho)

   !          write(*,'(10e15.7)') ri, ex0
   !          write(*,'(10e15.7)') ri, dfxdn, -dfxdr, dfxddn
   !          write(*,'(10e15.7)') ri, fs, dfsds, d2fsds
   !          write(*,'(10e15.7)') ri, rho, drho, d2rho
   !          write(*,'(10e15.7)') ri, s, d2rho/drho, 4.d0*drho/(3.d0*rho)

   !===       end  for olny sherical symmetric case ===

   endif

   if( iflag == 1 ) then
   !   --- relativistic corrections ---------------------------------------
      betar =  BB/rs
      YETA  =  SQRT( 1.0D0 + betar*betar )
      ZETA  =  LOG( betar + YETA )
      R2    =  ( YETA - ZETA/betar )/betar
      RCORRE=  1.D0 - 1.5D0*R2*R2
      RCORR =  -0.5D0 + 1.5D0*ZETA/betar/YETA

      ex = ex * RCORRE
      vx = vx * RCORR
   !   --------------------------------------------------------------------
   endif




!--- correlation terms -------------------------------------------------

   rs1h = dsqrt(rs)

   g0bb = rs1h*( bet01 + rs1h*( bet02 + rs1h*( bet03 + rs1h*bet04 )))
   g0b2 = 1.d0/(a0*g0bb)
   if( abs(g0b2) > 1.d-10 ) then
      g0ln = log( 1.d0 + g0b2 )
   else
      g0ln = g0b2
   endif
   g0rs1 = 1.d0 + alp01*rs
   g0rs  = -a0*g0rs1*g0ln

   ecrsz = g0rs
!  --- correlation energy ---
   ec    = ecrsz


   dg0dr1 = 0.5d0*bet01/rs1h + bet02 &
   + rs1h*(1.5d0*bet03 + 2.d0*bet04*rs1h)
   dg0drs = -a0*alp01*g0ln &
   + g0rs1*dg0dr1 / ( ( 1.d0 + g0b2 )*g0bb*g0bb )

   drsdn = - rs * c1b3

   decdn = drsdn * dg0drs


!==== gradient terms ================
   t  = t0*adrbro
   t2 = t*t / rho3
   if( t2 < 1.d-30 ) then
      vc = ec + decdn
   else

      h2 = betapw/h1

      eecrsz = exp(-ecrsz/h1)
      arsz   = h2 / ( eecrsz - 1.d0 )
      arsz2  = arsz*arsz

      t4 = t2*t2

      h0bs =             t2 + arsz *t4
      h0bb = 1.d0 + arsz*t2 + arsz2*t4
      h0b2 = h0bs/h0bb
      if( abs(h0b2) > 1.d-10 ) then
         h0ln = log( 1.d0 + h2*h0b2 )
      else
         h0ln = h2*h0b2
      endif
      h0trsz = h1*h0ln

      ech  = h0trsz

   !  --- correlation energy by gradient terms ---
      ec = ec + ech



      h0bsa =                t4
      h0bba = t2 + 2.d0*arsz*t4
      h0bst = 1.d0 + 2.d0*arsz *t2
      h0bbt = arsz + 2.d0*arsz2*t2

      dh0c1 = h2*h0bs + h0bb
      dh0c  = betapw/dh0c1
      dh0da = dh0c*( h0bsa - h0b2*h0bba )
      dh0dt = dh0c*( h0bst - h0b2*h0bbt )

      dadn = eecrsz*arsz2/betapw * decdn
      dtdn = - c7b3 * t2

      dh0dn = dh0da*dadn + dh0dt*dtdn


      dfcdn = ec + decdn + dh0dn



      dtddn = 2.d0*t2/adrbro

      dh0ddn = dh0dt*dtddn

      dfcddn = dh0ddn


   !=== from here, for olny sherical symmetric case ===

      h0bstt = 2.d0*arsz
      h0bsta = 2.d0*t2
      h0bbtt = 2.d0*arsz2
      h0bbta = 1.d0 + 4.d0*arsz*t2

      dh0c2  = -dh0dt/dh0c1
      dh0dtt = dh0c2*( h2*h0bst + h0bbt ) &
      + dh0c*( h0bstt &
      - ( h0bst*h0bbt + h0bs*h0bbtt - h0b2*h0bbt*h0bbt )/h0bb )
      dh0dat = dh0c2*( h2*h0bsa + h0bba ) &
      + dh0c*( h0bsta &
      - ( h0bsa*h0bbt + h0bs*h0bbta - h0b2*h0bbt*h0bba )/h0bb )

      dtdr   = t2*( 2.d0*d2rbro - c7b3 * drbro )
      dadr   = dadn * drbro

      dh0drt = dh0dtt*dtdr + dh0dat*dadr


      dtdrdn = dtddn*( d2rbro - c4b3 * drbro )

      dfcdr = dh0drt*dtddn + dh0dt*dtdrdn

   !     --- correlation potential ---
      vc = dfcdn - ( dfcdr + 2.d0*dfcddn/ri )*sign(1.d0,drho)

   !      if( iflag.eq.1 ) then
   !      write(181,'(10e20.12)') t2, dfcdn, dh0dn, dh1dn,
   !     &                        dfcdr, 2.d0*dfcddn/ri
   !      write(181,'(20e14.6)') t2,
   !     &    dfcdr, dh0drt, dh1drt, dtddn,  dh0dt, dh1dt, dtdrdn,
   !     &         dh1drs,  dsddn,           dh1ds,  dsdrdn
   !      endif

   !===       end  for olny sherical symmetric case ===

   endif

   return
end subroutine pbe

