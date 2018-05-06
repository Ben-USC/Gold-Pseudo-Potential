module common_variables
!--------------------------------------------------------------------
!          MSH   ...... total number of mesh point
!          MOR   ...... maximum number of electron orbitals
!          MXL   ...... maximum number of angular momentum
!                       quantum number, l
!          MXREF ...... maximum number of reference energy
!          MXVEL ...... maximum number of valence electron
!--------------------------------------------------------------------
   integer, parameter :: msh = 2001
   integer, parameter :: mor = 18
   integer, parameter :: mxl = 4
   integer, parameter :: mxref = 4
   integer, parameter :: mxvel = 4
   real(8) :: zatm                        ! atomic number
   real(8) :: xion                        ! valence ion
   real(8) :: wnlj(1:mor)                 ! number of electrons of each orbitals
   real(8) :: edel                        ! accuracy of eigenvalues
   real(8) :: vdel                        ! accuracy of potential
   real(8) :: phi                         ! mixing ratio
   real(8) :: rmax1                       ! maximum of radial mesh
   logical :: lsg                         ! .true. = read charge density
   logical :: lprn
   logical :: lplot
   logical :: latt                        ! latter correction
   integer :: nshl                        ! No. of orbitals
   integer :: nljc(1:mor)                 ! electron orbitals
   integer :: mesh1                       ! total No. of radial mesh
   integer :: nval                        ! No. of valences
   real(8) :: xh                          ! width of logarithmic mesh
   real(8) :: v0                          ! local potential at r = 0
   real(8) :: rloc                        ! cutoff radius for the local potential
   real(8) :: rus(1:mxref)                ! cutoff radius for USPP
   real(8) :: rnc(1:mxref)                ! cutoff radius for NCPP
   real(8) :: iref(1:mxref)               ! reference energy
   real(8) :: rtra                        ! transferability (?)
   real(8) :: v(1:msh)                    ! all electron potential
   real(8) :: vloc(1:msh)                 ! local potential
   real(8) :: pusnl(1:msh)                ! pseudo-wave-functions by GEE
   real(8) :: plj(1:mxref,1:msh)          ! radial AE wave functions
   real(8) :: rcl                         ! max(rcl1, rcl2)
   real(8) :: rcl1                        ! cutoff radius for Pus
   real(8) :: rcl2                        ! cutoff radius for Pnc
   real(8) :: ref(1:mxref)                ! reference energy
   real(8) :: eref(1:mxref)               ! reference energy
   real(8) :: p(1:msh)                    ! ???
   real(8) :: pae(1:mxref,1:msh)          ! radial AE wave functions
   real(8) :: pnc(1:mxref,1:msh)          ! NC pseudo-wave-functions
   real(8) :: pus(1:mxref,1:msh)          ! US pseudo-wave-functions
   real(8) :: d2p(1:mxref,1:msh)          ! second derivative of Pus
   real(8) :: bas(1:mxref,1:msh)          ! basis functions
   real(8) :: qch(1:mxref,1:mxref)        ! augmentation charge
   real(8) :: dele(1:mxref,1:mxref)       ! matrix elements D for USPP
   real(8) :: qau(1:mxref,1:mxref,1:msh)  ! augmentation functions for USPP
   real(8) :: qau2(0:2*mxl,1:mxref,1:mxref,1:msh)  ! augmentation functions for PAW
   real(8) :: chi(1:mxref,1:msh)          ! local functions
   real(8) :: sgc(1:msh)                  ! all-electron core-electron density
   real(8) :: sgv(1:msh)                  ! all-electron valence-electron density
   real(8) :: rhous(1:msh)                ! USPP valence-electron density
   real(8) :: ssg(1:msh)                  ! soft part of the pseudo-valence-electron density
   real(8) :: hsg(1:msh)                  ! hard part of the pseudo-valence-electron density
   integer :: nus, lus, nref, icl, ius, inc
   integer :: nrcut

! Projector augmented wave
   real(8) :: qcomp(0:2*mxl,mxref,mxref)
   real(8) :: fcomp,rcomp
   real(8) :: gf(0:2*mxl,msh), work1(msh)



   integer :: methodpp
   logical :: lncpp, luspp, lpaw

   real(8), save :: rmax, dx, rinf, dell
   integer, save :: mesh

   character(1) :: ms(1:7) =                                       &
& (/ 'S', 'P', 'D', 'F', 'G', 'H', 'I' /)

   character(2) :: aname(1:103) =                                  &
& (/ 'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',            &
&    'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA','SC','TI',  &
&    'V ','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE',  &
&    'BR','KR','RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PD',  &
&    'AG','CD','IN','SN','SB','TE','I ','XE','CS','BA','LA','CE',  &
&    'PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',  &
&    'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG','TL','PB',  &
&    'BI','PO','AT','RN','FR','RA','AC','TH','PA','U ','NP','PU',  &
&    'AM','CM','BK','CF','ES','FM','MD','NO','LR' /)

end module common_variables


   PROGRAM LDA
!----------------------------------------------------------------------
!                                                          1992/09/13
!                                                  update  1994/07/18
!                                                  by      F. Shimojo
!       ZATM ...... Atomic number
!       XION ...... valence ion
!       NSHL ...... total number of electron orbitals
!       NLJC ...... electron orbitals
!         (  NLJC = 100*n + 10*l + k  ( k : arbitrary )   )
!       WNLJ ...... number of electrons of each orbitals

!       MSH  ...... total number of mesh point
!       MOR  ...... maximum number of electron orbitals
!-----------------------------------------------------------------------
   use common_variables
!       PARAMETER  ( MOR = 18  )
   IMPLICIT REAL*8 (A-H,O-Z)
   CHARACTER    OFILE*13
!    DIMENSION  NLJC(MOR),WNLJ(MOR)
   DATA  JIN/ 01 /, JOT/ 02 /, J06/ 06 /
   DATA  INDER / 0 /

!    logical :: LSG, LPRN, LPLOT, LATT


   call read_inpufile


!----   open file for output  ( file name = 'out_atom' )
   NZN = IDNINT(ZATM)
   nacr = 2
   if( ANAME(NZN)(2:2) == ' ' ) nacr = 1
!      OFILE = ( 'data/out_'//ANAME(NZN)(1:nacr) )
   OFILE = ( 'out_'//ANAME(NZN)(1:nacr) )
   OPEN(JOT,FILE=OFILE)
!-----------------------------------------------------------------------


!     ------------------------------------------------------
!   CALL ATMLDA( ZATM, XION, NSHL, NLJC, WNLJ, EDEL, VDEL, &
!      &  PHI, RMAX1,mesh1, JSG, JPRN, JPLOT, JLAT, INDER, JOT, J06  )
   CALL ATMLDA( INDER, JOT, J06  )
!     ------------------------------------------------------


   WRITE(J06,210) INDER
   210 FORMAT( / ' ----- atomic calculation ended ( code :', &
   I2,' ) -----')


   STOP
   end




!   SUBROUTINE ATMLDA( ZATM, XION, NSHL, NLJC, WNLJ, EDEL, VDEL, &
!     &      PHI, RMAX1, mesh1, JSG, JPRN, JPLOT, JLAT, INDER, JOT, J06 )
   SUBROUTINE ATMLDA( INDER, JOT, J06 )
!-----------------------------------------------------------------------
!     Non-relativistic atomic calculation.
!        ( Local exchange-correlation potential )
!     Atomic calculation by the noumerov scheme with log-scale meshes.
!-----------------------------------------------------------------------
!       R(I)  ...... real mesh  =  RMAX*exp( (I - M)*FDX )
!       VR(I) ...... R(I) * potential V(r)
!       SG(I) ...... 4 * pai * R(I) * R(I) * electron density
!       XE(N) ...... eigenvalue of (n, l)-orbital
!-----------------------------------------------------------------------
   use common_variables
!      PARAMETER  ( MSH = 2001 )
!       PARAMETER  ( MOR = 18  )
   IMPLICIT REAL*8 (A-H,O-Z)
!   DIMENSION     NLJC(MOR), WNLJ(MOR)
   DIMENSION     MN(MOR), ML(MOR)
   DIMENSION     R(MSH), VR(MSH) ,PNL(MSH)
   DIMENSION     XE(MOR), SG(MSH)
   CHARACTER(1) ::   NUM(0:9)
   DIMENSION     VEXT(MSH), VHAR(MSH), VEXC(MSH), VCOR(MSH)
   DIMENSION     EEXC(MSH), ECOR(MSH)
   DIMENSION     EPART(5)

   dimension     drhox(msh),d2rhox(msh),rho(msh)
   dimension     usrho(msh)
!      DIMENSION     SGC(MSH), SGV(MSH)
   DIMENSION     SGNEW(MSH)
   DIMENSION     A(MSH), B(MSH), G(MSH), DA(5), DB(5)
   CHARACTER(20) ::   FNAME
!   LOGICAL ::       LPRN, LCNV, LATT
   LOGICAL ::        LCNV

!       PARAMETER  ( MXL = 2   )
!       PARAMETER  ( MXREF = 2   )
!       PARAMETER  ( MXVEL = 2   )
   real*8 ::   REFET(MXREF,0:MXL)
   integer :: NREFET(0:MXL)
   real*8 ::   RCUT(0:MXL)
   real*8 ::   CUTQ(0:MXL)
   CHARACTER     calwrd*20, award*12
   integer :: ICHK(MXVEL,0:MXL)
   real*8 ::   OCCN(MXVEL,0:MXL)
   real*8 ::   EIGAE(MXVEL,0:MXL)
   integer :: NOAE(MXVEL,0:MXL)
   CHARACTER(4) ::    ECFGRS(7)
   real*8 ::         ECFG(7)
   integer ::        NECFG
   CHARACTER(4) ::    VECGRS(7)
   real*8 ::         VECG(7)
   integer ::        NVECG
   real*8 ::  REFE(MXREF)

   real*8 ::  SGhard(msh)
   real*8 ::  psrl(msh,MXVEL)
   real*8 ::  fkai(msh,MXVEL)
   real*8 ::  psqr(msh,MXVEL)
   real*8 ::  pnlae(msh,MXVEL)
   real*8 ::  ppnlae(msh,MXVEL)
   real*8 ::  lpsrl(msh,MXVEL,0:MXL)
   real*8 ::  lfkai(msh,MXVEL,0:MXL)
   real*8 ::  lpsqr(msh,MXVEL,0:MXL)
   real*8 ::  lpnlae(msh,MXVEL,0:MXL)
   real*8 ::  lppnlae(msh,MXVEL,0:MXL)

! rojector augmented wave
   real*8 ::  WORK(msh)


!   logical :: LSG, LPLOT
!-----------------------------------------------------------------------


!--- Latter correction -------------------------------------------------
   XLATTR  =   1.D0
!--- parameters of real space ------------------------------------------
!         MESH  =  MSH
!         XH    =   100.d0 !73.D0
!         RMAX  =   100.d0 !60.D0

   MESH  =  mesh1
   RMAX  =   rmax1
   RINF  =   75.D0
!          write(*,*) xh,rmax
!--- parameters of convergence -----------------------------------------
!         PHI   =    0.3D0
!         EDEL  =    5.D-6
!         VDEL  =    1.D-4
   NC1   =  300
   TIME  =   60.D0
!-----------------------------------------------------------------------

!--- prepare radial meshes ---------------------------------------------
   DX   =  1.D0 / XH
   D    =  EXP( DX )
   CR   =  RMAX / D ** MESH
   DO 2 I = 1, MESH
      VR(I)  =  0.D0
      CR     =  CR * D
      R(I)   =  CR
   2 END DO
!       write(*,*)'r(1):',r(1)
!--- prepare starting charge density -----------------------------------
   CALL ICDENS( LSG, SG, R, ZATM, XION, MESH, JOT )
!-----------------------------------------------------------------------


!      CALL SEC( T1 )
!           T2  =  T1

   CYCLET  =  0.D0
!            LATT    =  JLAT .EQ. 1
!            LPRN    =  JPRN .EQ. 1


   NZN = IDNINT(ZATM)
   WRITE(J06,200) ANAME(NZN)
   WRITE(JOT,200) ANAME(NZN)
   200 FORMAT( / ' Self-consistent LDA calculation for atom: ', A2 )
   IF( LPRN ) THEN
      WRITE(J06,38) MESH, RMAX, PHI, ZATM, RINF, EDEL, XION, &
   !     &                  XH, VDEL, NC1, TIME, JLAT
      XH, VDEL, NC1, TIME, LATT
      WRITE(JOT,38) MESH, RMAX, PHI, ZATM, RINF, EDEL, XION, &
   !     &                  XH, VDEL, NC1, TIME, JLAT
      XH, VDEL, NC1, TIME, LATT
      38 FORMAT( / ' (initial data):'/ &
   !      &     1H ,9X,'MESH =',I4,7X,  'RMAX =',F7.2,5X,'   PHI =',F11.7/
   !      &     1H ,9X,'   Z =',F6.1,5X,'RINF =',F7.2,5X,'  EDEL =',F11.7/
   !      &     1H ,9X,'XION =',F6.1,5X,'   H =',F7.2,5X,'  VDEL =',F11.7/
   !     &     1H ,9X,' NC1 =',I4,7X,  'TIME =',F7.2,5X,'  JLAT =',I3   /
      &      1H ,9X,'MESH =',I5,6X,  'RMAX =',F7.2,5X,'   PHI =',F11.7/ &
      &      1H ,9X,'   Z =',F7.1,4X,'RINF =',F7.2,5X,'  EDEL =',F11.7/ &
      &      1H ,9X,'XION =',F7.1,4X,'   H =',F7.2,5X,'  VDEL =',F11.7/ &
      &      1H ,9X,' NC1 =',I5,6X,  'TIME =',F7.2,5X,'  LATT =',L3   / &
      / ' (electron configuration):'/ &
      &      1H ,6X,'NLJC',3X,'orbit',4X,'WNLJ',6X,'energy (Ryd.)' )
   ENDIF

!--- electron configuration --------------------------------------------
   SUM  =  0.D0
   DO 31 N = 1, NSHL
      XE(N)  =  0.D0
      MN(N)  =  NLJC(N) / 100
      ML(N)  =  NLJC(N) / 10 - 10 * MN(N)
      IF( LPRN ) THEN
         WRITE(J06,36) N, NLJC(N), MN(N), &
         MS( ML(N)+1 ), WNLJ(N), XE(N)
         WRITE(JOT,36) N, NLJC(N), MN(N), &
         MS( ML(N)+1 ), WNLJ(N), XE(N)
         36 FORMAT( 1H , I3, 3X, I4, I6, A1, 3X, F7.4, F16.7 )
      ENDIF
      SUM  =  SUM + WNLJ(N)
   31 END DO
!--- abnormal return ---------------------------------------------------
   IF( ABS( SUM + XION - ZATM ) > 1.D-8 ) THEN
      WRITE(J06,39) SUM, XION, ZATM
      39 FORMAT( / 'total number of electrons in error:'/ &
      &              10X, ' SUM, XION, Z  =', 2X, 3D16.7 )
      INDER = 2
      RETURN
   ENDIF
!-----------------------------------------------------------------------
   SM4   =  0.5D0
   SM3   =  0.D0
!-----------------------------------------------------------------------

   IF( LPRN ) THEN
      WRITE(J06,21)
      WRITE(JOT,21)
   !    21    FORMAT( / ' (convergence test):'/7X,'CYCLET',5X,'PH',6X,'RLC',
   !      &              5X,'VERR(NERR)',7X,'SUME',5X,'electrons')
      21 FORMAT( / ' (convergence test):'/5X,'CYCLET',4X,'PH',8X,'RLC', &
      &               4X,'VERR      (NERR)',3X,'SUME',7X,'electrons')
   ENDIF




!-----------------------------------------------------------------------
!                    outer iteration starts.
!-----------------------------------------------------------------------
!      PH  =  2.D0 - PHI
   PH =    PHI


   DO 350 JCYCL = 1, NC1
   !--- potential function from charge density ----------------------------
      S2  =  LOG( SG(2) / SG(1) ) / DX
      DO 81 I = 1, 2
         A(I)  =  SG(I) * R(I) / ( S2 + 1.D0 )
         B(I)  =  SG(I) / S2
         DB(I) =  DX * SG(I) / 3.D0
         DA(I) =  DB(I) * R(I)
      81 END DO
      DO 82 I = 3, MESH
         DB(3) =  DX * SG(I) / 3.D0
         DA(3) =  DB(3) * R(I)
         A(I)  =  A(I-2) + DA(3) + 4.D0 * DA(2) + DA(1)
         B(I)  =  B(I-2) + DB(3) + 4.D0 * DB(2) + DB(1)
         DO 82 L = 1, 2
            DA(L) =  DA(L+1)
            DB(L) =  DB(L+1)
      82 END DO
      QE  =  A(MESH)
      BM  =  B(MESH)

   !--- max error in potential ( VERR ) and Latter correction -------------
      RU  =  - 2.D0 * ( XLATTR + XION )
      IF( RU < - 2.D0 * ZATM ) RU  =    0.D0
      VERR  =  0.D0

      pi=acos(-1.d0)
      do i=1,mesh
         rho(i) = SG(i)/( 4.D0 * pi * R(I) * R(I) )
      end do
      call fd1(1,mesh,mesh,rho,dx,drhox)
      call fd2(1,mesh,mesh,rho,dx,d2rhox)


      DO 94 I = 1, MESH
         IF( SG(I) < 1.D-30 ) THEN
            VEXC(I) = 0.0
            VCOR(I) = 0.0
            EEXC(I) = 0.0
            ECOR(I) = 0.0
         ELSE
         !---         LDA
         !             RS  =  ( 3.D0 * R(I) * R(I) / SG(I) )**( 1.D0 / 3.D0 )
         !             CALL VXC( RS, EX, VX, EC, VC )

         !---         GGA
         ! drhox  =  d rho/ dx
         ! d2rhox = d^2rho/ dx^2

            drho   = drhox(i) / R(i)
            d2rho  = ( d2rhox(i) - drhox(i) ) / ( R(i) * R(i) )

            CALL pbe(rho(i), drho, d2rho, r(i), EX, VX, EC, VC, 0)

            VEXC(I) = 2.0D0 * R(I) * VX
            VCOR(I) = 2.0D0 * R(I) * VC
            EEXC(I) = 2.0D0 * R(I) * EX
            ECOR(I) = 2.0D0 * R(I) * EC
         ENDIF

         VEXT(I)  =  - 2.D0 * ZATM
         VHAR(I)  =  + 2.D0 * ( A(I) + R(I) * ( BM - B(I) ) )
         RV  =  VEXT(I) + VHAR(I) + VEXC(I) + VCOR(I)
         IF( LATT ) THEN
            IF( RV > RU ) THEN
               RV  =  RU
            ELSE
               NLC =  I
            ENDIF
         ENDIF
         DV  =  ABS( RV - VR(I) )! / R(I)
      ! eck
      !         sgc(i) = dv
      ! eck
         if( r(i) > 1.d-02 ) then
            IF( DV > VERR ) THEN
               VERR  =  DV
               NERR  =  I
            ENDIF
         end if
         VR(I) =  RV
      94 END DO
   ! eck
   !      write(99,'(i4,10e14.6)') jcycl, (sgc(i*100+1), i=0,9)
   ! eck
      RLC   =  0.D0
      IF( LATT ) THEN
         NLC   =  MIN( NLC + 1, MESH )
         RLC   =  R( NLC )
      ENDIF
   !--- convergence criterion ---------------------------------------------

      LCNV  =  VERR < VDEL

      IF( LCNV ) THEN
         DELL = EDEL
      ELSE
         DELM = 0.01D0
         IF( ABS( SM3 - SM4 ) < 0.1D0   )  DELM = 0.001D0
         IF( ABS( SM3 - SM4 ) < 0.001D0 )  DELM = EDEL
         DELL = MIN( DELM, EDEL * SQRT( VERR / VDEL ) )
      ENDIF
   !-----------------------------------------------------------------------
   !--- eigenvalues and new charge densities ------------------------------
      SUME  =  0.D0
   !H
      IF( LCNV ) THEN
         DO 83 I = 1, MESH
            SGC(I)  =  0.D0
            SGV(I)  =  0.D0
         83 END DO
      ENDIF
   !H
      DO 9  I = 1, MESH
         SGNEW(I)  =  0.D0
      9 END DO
      DO I = 1, 5
         EPART(I)  =  0.D0
      ENDDO
      DO 13 N = 1, NSHL
         NQ  =  MN(N)
         LQ  =  ML(N)
         IF( ABS( XE(N) ) > 1.D-20 ) THEN
            EIG    =  XE(N)
         ELSE
            EIG    =  0.D0
         ENDIF

      !--- solve Schroedinger equation for given (n,l) -----------------------

         ! CALL SCHROE( ZATM, NQ, LQ, EIG, R, VR, PNL, A, B, G, KI, KJ, J06 )
         CALL SCHROE( NQ, LQ, EIG, R, VR, PNL, A, B, G, KI, KJ, J06 )

      !--- kinetic energy
         CALL KINE( LQ, PNL, R, EKIN, A )
         EPART(1) = EPART(1) + WNLJ(N) * EKIN

         XE(N)  =  EIG
         SUME  =  SUME + WNLJ(N) * XE(N)
      !         IF( WNLJ(N) .GT. 1.D-8 ) THEN
      !H

         IF( LCNV ) THEN
            JCORE  =  MOD( NLJC(N), 10 )
            DO 85 I = 1, MESH
               SGI  =  WNLJ(N) * PNL(I) * PNL(I)
               IF( JCORE < 1 ) THEN
                  SGC(I)  =  SGC(I) + SGI
               ELSE
                  SGV(I)  =  SGV(I) + SGI
               ENDIF
            85 END DO


            IF( LPLOT ) THEN
            !--- output PNL ( : wave function ) -------------------------------
               DATA  JPL / 9 /
               DATA  (NUM(L),L=0,9) &
               / '0','1','2','3','4','5','6','7','8','9' /
               nacr = 2
               if( ANAME(NZN)(2:2) == ' ' ) nacr = 1
               FNAME = ( (((( 'data/'//ANAME(NZN)(1:nacr) )//'_' )// &
               NUM(NQ) )//MS(LQ+1) )//'.wvf' )
               WRITE(J06,*) FNAME
               OPEN(JPL,FILE=FNAME)

               WRITE(JPL,2010) 'Wave Function                 ', &
               ANAME(NZN),NQ,MS(LQ+1)
               WRITE(JPL,2001)
               WRITE(JPL,2020) ( NLJC(I), I = 1, NSHL )
               WRITE(JPL,2021) ( WNLJ(I), I = 1, NSHL )
               WRITE(JPL,2030) 'e', 'eigenvalue : ',EIG
               WRITE(JPL,2030) 't', 'turning p. : ',R(KI), KI
               WRITE(JPL,2030) 'r', 'R_inf      : ',R(KJ), KJ

               WRITE(JPL,2090) '#  r [ bohr ]     ','  w.f.            '
               DO 15 I = 1, MESH
                  WRITE(JPL,2000) R(I), PNL(I)
               15 END DO
               ENDFILE(JPL)
               CLOSE(JPL)
            ENDIF
            2010 FORMAT('#',5X,A31,': ',A2, &
            '   ( for nl =',I2,A1,' )'/ &
            '#',5X,'   by non-relativistic atomic calculation.')
            2001 FORMAT('#',5X,'---- valence electron configuration ----')
            2000 FORMAT(3D18.10)
            2090 FORMAT(3A18)
            2020 FORMAT('#',5X,18I5)
            2021 FORMAT('#',5X,18F5.1)
            2030 FORMAT('#',A1,4X,A13,D18.10,I5)
         !------------------------------------------------------------------
         ENDIF
      !H
         DO 12 I = 1, MESH
            SGNEW(I)  =  SGNEW(I) + WNLJ(N) * PNL(I)**2
         12 END DO
      ! eck
      !      S2  =  LOG( SGnew(2) / SGnew(1) ) / DX
      !      DO I = 1, 2
      !         A(I)  =  SGnew(I) * R(I) / ( S2 + 1.D0 )
      !         DB(I) =  DX * SGnew(I) / 3.D0
      !         DA(I) =  DB(I) * R(I)
      !      end do
      !      DO I = 3, MESH
      !         DB(3) =  DX * SGnew(I) / 3.D0
      !         DA(3) =  DB(3) * R(I)
      !         A(I)  =  A(I-2) + DA(3) + 4.D0 * DA(2) + DA(1)
      !         DO L = 1, 2
      !            DA(L) =  DA(L+1)
      !         end do
      !      end do
      !      write(*,*) '# of el. ', A(MESH)
      ! eck
      !         ENDIF
      13 END DO
   !--- energy ------------------------------------------------------------
      CALL ENEGYF( VEXT, VHAR, VEXC, VCOR, EEXC, ECOR, SGNEW, A, &
      SUME, EPART )
      IF( LCNV ) GOTO 300
   !--- charge density for next iteration ---------------------------------
      PH  =  ( PH + PHI ) / 2.D0
      DO 17 I = 1, MESH
         SG(I) =  ( 1.D0 - PH ) * SG(I) + PH * SGNEW(I)
      17 END DO
      SM4   =   SM3
      SM3   =   SUME

   !      CALL SEC( T3 )
   !      CYCLET =  T3 - T2
   !      T2    =   T3
      IF( LPRN ) THEN
         WRITE(J06,62) JCYCL, CYCLET, PH, RLC, VERR, NERR, SUME, QE
         WRITE(JOT,62) JCYCL, CYCLET, PH, RLC, VERR, NERR, SUME, QE
         62 FORMAT( I3, F7.3, ' sec', F9.6, F8.4, E11.4, &
      !      &                      '(', I3, ')', F11.4, F12.7 )
         '(', I4, ')', F12.4, F12.7 )
      ENDIF


   350 END DO
!--- JCYCL - loop ended ------------------------------------------------
!--- abnormal return ---------------------------------------------------
   WRITE(J06,47) NC1
   47 FORMAT( / ' problem not converged after', I4, ' cycles.' )
   INDER  =  1
   RETURN


!-----------------------------------------------------------------------
!                          problem converged.
!-----------------------------------------------------------------------
   300 INDER  =  0
   ELV = 0.136058D+02
   IF( LPRN ) THEN
      WRITE(J06,48) ANAME(NZN), ZATM, SUM, XION
      WRITE(JOT,48) ANAME(NZN), ZATM, SUM, XION
      WRITE(J06,49)
      WRITE(JOT,49)
      48 FORMAT( / ' ---- LDA calculation for ',A2,' ( Z =',F7.2,' ) ', &
      &       35('-')/ &
      &       6X,'Total No. of electrons :',F11.6,'  ( ion =',F10.6,' )' )
      49 FORMAT( / ' (problem converged):'/ &
      &              1H , 6X, 'orbit', 3X, 'electrons', &
      &                   3X, 'energy (Ryd.)',10X,'(eV)' )
      DO N = 1, NSHL
         WRITE(J06,50) N, MN(N), MS( ML(N)+1 ), WNLJ(N), &
         XE(N), XE(N)*ELV
         WRITE(JOT,50) N, MN(N), MS( ML(N)+1 ), WNLJ(N), &
         XE(N), XE(N)*ELV
         50 FORMAT( 1H , I3, I5, A1, 4X, F10.6, F16.7, F14.5 )
      ENDDO
   ENDIF

   ETOT = 0.D0
   DO I = 1, 5
      ETOT = ETOT + EPART(I)
   ENDDO
   WRITE(J06,1100) sume, sume*ELV, &
   ( EPART(I), I = 1, 5 ), &
   ( EPART(I)*ELV, I = 1, 5 )
   WRITE(JOT,1100) sume, sume*ELV, &
   ( EPART(I), I = 1, 5 ), &
   ( EPART(I)*ELV, I = 1, 5 )
   1100 FORMAT(/'  * Total energy = ',F15.6,' ( Ryd. )',F13.4,' (eV)'/ &
   '  * Energy parts ( Ryd. / eV )'/ &
   &        6X,'     Kinetic        External       Hartree        ', &
   'Exchange    Correlation'/ &
   &        4X,5F15.6/4X,5(F13.4,2X)/ &
   ' ',78('-') )

!--- output data  ------------------------------------------------------
   WRITE(JOT,1010)
   WRITE(JOT,1000) ( R(J), J = 1, MESH )
   WRITE(JOT,1020)
   WRITE(JOT,1000) ( VR(J), J = 1, MESH )
   WRITE(JOT,1030)
   WRITE(JOT,1000) ( SG(J), J = 1, MESH )
!-----------------------------------------------------------------------


   if(LATT) RETURN



   J06 = 6
   JPL = 10

!         ZATM = atomic number
!         RCUT(L) = cutoff radius for the USPP wave functions
!         NREFET(L) = the number of reference energies
!         REFET(I,L) = the I-th reference energy
!         CUTQ(L) = cutoff radius for the NCPP wave functions
!         RCUTL = cutoff radius for the local potential

   iter = 1
   do n = 1, NSHL
      if ( MOD( NLJC(N), 10 ) ==1) then
         LQ  =  ML(N)
         eq  =  xe(n)

         RCUT(LQ) = rus(iter)

         NREFET(LQ) = iref(iter)

         do i = 1 , nrefet(LQ)
            if(i ==1) then
               REFET(I,LQ) = eq
            else if(i ==2) then
               refet(i,LQ) = ref(iter)
            end if
         enddo

         CUTQ(LQ) = rnc(iter)

         RCUTL = rloc

         iter = iter +1
      end if
   enddo





   do i = 1, MXVEL
      do ii = 0, MXL
         ICHK(i,ii) = 0
         OCCN(i,ii) = 0.d0
      enddo
   enddo

!--- electron configuration --------------------------------------------
   SUM  =  0.D0
   SUMC =  0.D0
   NECFG = 0
   NVECG = 0
   DO 41 N = 1, NSHL
      MNN  =  NLJC(N) / 100
      MLN  =  NLJC(N) / 10 - 10 * MNN
      SUM  =  SUM + WNLJ(N)

      JCORE  =  MOD( NLJC(N), 10 )
      IF( JCORE <= 0 ) THEN
         SUMC =  SUMC + WNLJ(N)
      ELSE
         KKY = 1
         IF( JCORE == 1 ) THEN
         !    ----------- JCORE .eq. 1 --------------
            IF( ICHK(KKY,MLN) >= 1 ) stop 'error'

            NECFG = NECFG + 1
         !                 ECFGRS(NECFG) = (( NUMBER(MNN)//MS( MLN+1 ) )
            ECFGRS(NECFG) = (( NUM(MNN)//MS( MLN+1 ) ) &
            //'  ' )
            ECFG(NECFG) = WNLJ(N)
         ELSE
         !    ----------- JCORE .eq. 2 --------------
            KKY = KKY + 1
            3000 IF( ICHK(KKY,MLN) == 1 ) THEN
               KKY = KKY + 1
               GO TO 3000
            ENDIF

            NVECG = NVECG + 1
         !                 VECGRS(NVECG) = (( NUMBER(MNN)//MS( MLN+1 ) )
            VECGRS(NVECG) = (( NUM(MNN)//MS( MLN+1 ) ) &
            //'  ' )
            VECG(NVECG) = WNLJ(N)
         !    ---------------------------------------
         ENDIF
         ICHK(KKY,MLN) = 1
         OCCN(KKY,MLN) = WNLJ(N)
         EIGAE(KKY,MLN)= XE(N)
         NOAE(KKY,MLN) = MNN
      ENDIF
   41 END DO
   ZV = ZATM - SUMC




   NZN = IDNINT(ZATM)
   nacr = 2
   if( ANAME(NZN)(2:2) == ' ' ) nacr = 1
   D    =  EXP( DX )
   CR   =  RMAX / D ** MESH

   rcutmx = 0.d0
   DO L = 0, MXL
      IF( ICHK(1,L) > 0 ) rcutmx = max( rcutmx, RCUT(L) )
   end do
   KI = NINT( LOG(rcutmx/CR)/DX )
   KIMAX = MIN( KI + 50, MESH )


   calwrd = 'scalar relativistic '
   award = ' GGA (PBE)  '



   icl = KIMAX - 50


!-----------------------------------------------------------------------
   if( lpaw ) then
!--define of rcomp-- for PAW
       rcomp = rcutmx/fcomp
       do i = 1, MESH
          if( r(i) > rcomp ) then
             nrcut = i - 1
             exit
          end if
       end do
       write(*,'(4x,a,f13.6,a)') "fcomp    =", fcomp
       write(*,'(4x,a,f13.6,a)') "rcutmx   =", rcutmx, &
       ' ...... max(rcut)'
       write(*,'(4x,a,f13.6,a)') "rcomp    =", rcomp, &
       ' ...... fcomp/max(rcut)'
       write(*,'(4x,a,i13,a)')   "NRCUT    =", NRCUT, &
       ' ...... r(NRCUT) = rcomp'
       write(*,'(4x,a,f13.6,a)') "r(NRCUT) =", r(NRCUT)

!--making g functions--
       write( unit=j06, fmt='(/, a,i4)' ) '  * G-function'
       do i = 0, 2*MXL, 2
          ! call GFUNC( i, nrcut, GF(i,1:MESH)
          ! &             , WORK(i:MESH), R, DX, MESH )
          ! call GFUNC( i, nrcut, GF_, WORK1, R, DX, MESH )
          call gfunc( i, nrcut, work1, work, r, dx, mesh )
          gf(i,1:mesh) = work1(1:mesh)
          write( unit=j06, fmt='(4x,a,i3)' ) 'l =', i
       end do
   end if

   call vlocal(r)

   do i = 1 ,mesh
      rhous(i) = 0.d0
      ssg(i) = 0.d0
      hsg(i) = 0.d0
   end do

   do lq = 0, MXL
      do KKY = 1, MXVEL

         if( ichk(kky,lq) == 1 ) then
            NQ  =  noae(kky,lq)
            write(*,'(/, 3a)') '  * ', NUM(NQ), MS(lq+1)
         end if
      !----------------------------------------------------------------------
         if( ichk(kky,lq)==1 .AND. kky==1 ) then
            NQ  =  noae(kky,lq)
            eq  =  eigae(kky,lq)
            vnoe = occn(kky,lq)
            nref = nrefet(lq)
            rcl1 = rcut(lq)
            rcl2 = cutq(lq)

            do i = 1, nref
               if (i==1) then
                  eref(i) = refet(i,lq)
               else
                  eref(i) = refet(i,lq)
                  iter = iter + 1
               end if
            end do
         !     ------------------------------------------------------------------
            if( lncpp ) then
               call reference_schroe(r, lq)
               call rrkj3( r, lq )

               do i = 1, icl+50
                  do j = 1, nref
                     lpsqr(i,j,lq) = pnc(j,i)
                  end do
               end do
            elseif( luspp ) then
               call reference_schroe(r, lq)
               call rrkj2( r, lq )

               do i = 1, icl+50
                  do j = 1, nref
                     lpsrl(i,j,lq) = pus(j,i)
                     lpnlae(i,j,lq) = plj(j,i)
                  end do
               end do
               call rrkj3( r, lq )

               do i = 1, icl+50
                  do j = 1, nref
                     lpsqr(i,j,lq) = pnc(j,i)
                  end do
               end do
            elseif( lpaw ) then
               call reference_schroe(r, lq)
               call rrkj2( r, lq )

               do i = 1, icl+50
                  do j = 1, nref
                     lpsrl(i,j,lq) = pus(j,i)
                     lpnlae(i,j,lq) = plj(j,i)
                  end do
               end do
            end if

         !     ------------------------------------------------------------------
            if( luspp .or. lpaw ) then
               if( luspp ) then
                   iflag = 1
               elseif( lpaw ) then
                   iflag = 2
               end if

               call ulps( r, lq, iflag )

               do i = 1, icl+50
                  do j = 1, nref
                     lfkai(i,j,lq) = chi(j,i)
                  end do
               end do

               ! call gee( zatm, nq, lq, r, eq )
               call gee( nq, lq, r, eq )

               do i = 1, mesh
                  lppnlae(i,kky,lq) = pusnl(i)
               end do

               call ele1( r, usrho, vnoe, iflag )

               do i = 1 ,mesh
                  rhous(i) = rhous(i) + occn(kky,lq) * usrho(i)
               end do

               write(*,'(/, 40a)') ' ---- Estimate plane-wave cutoff energies ',  &
               &                   ( '-', i = 1, 36 )
               call kec( lq, r )

               call afq2( lq, r, iflag )

               write(*,'(/, 50a)') ' ---- Check for transferability ',  &
               &                   ( '-', i = 1, 46 )
               call tran( lq, r )

               call tran1( lq, r )
            end if

         !-----------------------------------------------------------------------
         else if( ichk(kky,lq) == 1 .AND. kky==2) then

            NQ  =  noae(kky,lq)
            eq  =  eigae(kky,lq)
            vnoe = occn(kky,lq)
         !     ------------------------------------------------------------------
            if( luspp .or. lpaw ) then
                if( luspp ) then
                    iflag = 1
                elseif( lpaw ) then
                    iflag = 2
                end if

               ! call gee( zatm, nq, lq, r, eq )
               call gee( nq, lq, r, eq )

                do i = 1, mesh
                   lppnlae(i,kky,lq) = pusnl(i)
                end do

                call ele1( r, usrho, vnoe, iflag )

                do i = 1 ,mesh
                   rhous(i) = rhous(i) + OCCN(kky,lq) * usrho(i)
                end do
            end if

         !     ------------------------------------------------------------------
         end if
      end do
   end do
!     ------------------------------------------------------------------
   if( luspp .or. lpaw ) then
       write(*,'(/, 50a)') ' ---- Estimate electronic density ',  &
       &                   ( '-', i = 1, 44 )
       call ele2( r )
   end if
!     ------------------------------------------------------------------


!---------- output general data ------------------------------------------
   lsomax = -4
   do L = 0, MXL
      IF( ICHK(1,L) > 0 .OR. ICHK(2,L) > 0 ) lmax = L
   enddo
   FNAME = ( ANAME(NZN)(1:nacr)//'_A.wps' )
   WRITE(J06,*) FNAME
   OPEN( JPL, FILE = FNAME )
   write(JPL,'(5x,d18.10,3i7)') ZV, lmax, lsomax, lmax
   write(JPL,'(5x,i7,2d22.14)') MESH, RMAX, DX
   do l = 0, lmax
      jmax = 0
      do jj = 1, MXVEL
         IF( ICHK(jj,l) >= 1 ) jmax = jj
      enddo
      write(JPL,'(5x,2i7)') l, jmax
   enddo
   ENDFILE(JPL)
   CLOSE(JPL)
!------------------------------------------------------------------------




   ldo: do L = 0, lmax
      IF( ICHK(1,L) >= 1 ) then

         MJ2 = 2 * L + 1
         NREFE = NREFET(L)
         DO NI = 1, NREFE
            REFE(NI) = REFET(NI,L)
         enddo

         DO MM = 1, KIMAX
            do I = 1,NREFE
               PNLAE(MM,I) = lpnlae(mm,i,l) !the all-electron wave-function
               if( luspp .or. lpaw ) then
                  PSRL(MM,I) = lpsrl(mm,i,l) !the USPP pseudo-wave-function
                  FKAI(MM,I) = lfkai(mm,i,l) !the local function
               end if
               if( lncpp ) then
                  PSQR(MM,I) = lpsqr(mm,i,l) !the NCPP pseudo-wave-function
               end if
            end do
         end do

         if( luspp .or. lpaw ) then
            do jj = 1, MXVEL
               if( ICHK(jj,L) >= 1 ) then
                  DO I = 1, MESH
                     PPNLAE(I,jj) = lppnlae(i,jj,l) !the normalized USPP pseudo-wave-function obtained by solving the generalized KS equiation
                  end do
               end if
            end do
         end if

         KI = NINT( LOG(RCUT(L)/CR)/DX )
         IF( CUTQ(L) > R(1) ) THEN
            KIIQ  = NINT( LOG(CUTQ(L)/CR)/DX )
         else
            KIIQ  = 0
         endif
         KIQ  = KIIQ
         KIQM = 0

      !---------- output pseudo-wave-functions -------------------------------
if( luspp .or. lpaw ) then
         FNAME = ( ((( ANAME(NZN)(1:nacr)//'_' )//MS(L+1) )//NUM(MJ2) ) &
         //'.wps' )
         WRITE(J06,'(1X,A20)') FNAME
         OPEN( JPL, FILE = FNAME )
         WRITE(JPL,4000) ANAME(NZN),award,calwrd,'#','( for l =',L
         WRITE(JPL,4001)
         WRITE(JPL,2002) ( ECFGRS(I2), I2 = 1, NECFG )
         WRITE(JPL,2003) ( ECFG(I2), I2 = 1, NECFG )
         if( NVECG > 0 ) then
            WRITE(JPL,2007)
            WRITE(JPL,2002) ( VECGRS(I2), I2 = 1, NVECG )
            WRITE(JPL,2003) ( VECG(I2), I2 = 1, NVECG )
         endif
         WRITE(JPL,2004) ( '----------', I2 = 1, NECFG )
         WRITE(JPL,2005) '#C',KI,R(KI),' : R_cut'
         WRITE(JPL,2006) '#L',1.D0/DX, ' : 1 / logarithmic mesh'
         WRITE(JPL,2006) '#Z',ZV,      ' : core charge         '
         WRITE(JPL,4010) '#N',NREFE,   ' : No. of reference energies'
         WRITE(JPL,4020) ( '#E',I,REFE(I),'[Ryd.]', I = 1, NREFE )
         WRITE(JPL,4030) '  r [ bohr ]      ', &
         ( (( ' pseudo-w.f.('//NUM(NI) )//') ' ), &
         (( 'local-w.f.('//NUM(NI) )//')  ' ), &
         (( 'p-w.f.minor('//NUM(NI) )//')' ), &
         NI = 1, NREFE )
         FIHT = 0.d0
         DO 910 MM = 1, KIMAX
            WRITE(JPL,2200) R(MM), &
            ( PSRL(MM,I), FKAI(MM,I), FIHT, I=1,NREFE )
         910 END DO
         ENDFILE(JPL)
         CLOSE(JPL)
end if
      !-----------------------------------------------------------------------


      !---------- output norm-conserving pseudo-wave-functions ---------------
if( lncpp ) then
         FNAME = ( ((( ANAME(NZN)(1:nacr)//'_' )//MS(L+1) )//NUM(MJ2) ) &
         //'.nwp' )
         WRITE(J06,'(1X,A20)') FNAME
         OPEN( JPL, FILE = FNAME )
         WRITE(JPL,4000) ANAME(NZN),award,calwrd,'#','( for l =',L
         WRITE(JPL,4001)
         WRITE(JPL,2002) ( ECFGRS(I2), I2 = 1, NECFG )
         WRITE(JPL,2003) ( ECFG(I2), I2 = 1, NECFG )
         if( NVECG > 0 ) then
            WRITE(JPL,2007)
            WRITE(JPL,2002) ( VECGRS(I2), I2 = 1, NVECG )
            WRITE(JPL,2003) ( VECG(I2), I2 = 1, NVECG )
         endif
         WRITE(JPL,2004) ( '----------', I2 = 1, NECFG )
         WRITE(JPL,2005) '#C',KI,R(KI),' : R_cut'
         RKIQ  = 0.D0
         RKIQM = 0.D0
         IF( KIQ >= 1 ) RKIQ  = R(KIQ)
         IF( KIQM >= 1 ) RKIQM = R(KIQM)
         WRITE(JPL,2005) '#P',KIQ,RKIQ,' : R_cut',' for pse', &
         'udizatio','n of maj','or part '
         WRITE(JPL,2005) '#Q',KIQM,RKIQM,' : R_cut',' for pse', &
         'udizatio','n of min','or part '
         WRITE(JPL,2006) '#L',1.D0/DX, ' : 1 / logarithmic mesh'
         WRITE(JPL,2006) '#Z',ZV,      ' : core charge         '
         WRITE(JPL,4010) '#N',NREFE,   ' : No. of reference energies'
         WRITE(JPL,4020) ( '#E',I,REFE(I),'[Ryd.]', I = 1, NREFE )
         WRITE(JPL,4030) '   << norm-conserv','ing pseudo-wave-fu', &
         'nctions >>        '
         WRITE(JPL,4030) '  r [ bohr ]      ', &
         ( (( 'major part('//NUM(NI) )//') ' ), &
         (( '  minor part('//NUM(NI) )//') ' ), &
         NI = 1, NREFE )
         PSQI = 0.d0
         DO MM = 1, KIMAX
            WRITE(JPL,2200) R(MM), &
            ( PSQR(MM,I), PSQI, I = 1, NREFE )
         enddo
         ENDFILE(JPL)
         CLOSE(JPL)
end if
      !-----------------------------------------------------------------------


      !--- output pseudo-wave-functions obtained by solving ------------------
      !--- the pseudo Schroedinger equations ---------------------------------
if( luspp .or. lpaw ) then
         NQ = L
         LQ = L
         do jj = 1, MXVEL
            NQ = NQ + 1
            if( ICHK(jj,L) >= 1 ) then

               if( NQ == LQ+1 ) then
                  FNAME = ( ((( ANAME(NZN)(1:nacr)//'_' )//MS(L+1) )//NUM(MJ2) ) &
                  //'.pwf' )
               else
                  FNAME = (((((( ANAME(NZN)(1:nacr)//'_' )//MS(L+1) )//NUM(MJ2) ) &
                  //'_' )//NUM(NQ-LQ) )//'.pwf' )
               endif
               WRITE(J06,'(1X,A20)') FNAME
               OPEN( JPL, FILE = FNAME )
               WRITE(JPL,4000) ANAME(NZN),award,calwrd,'#','( for l =',L
               WRITE(JPL,4001)
               WRITE(JPL,2002) ( ECFGRS(I2), I2 = 1, NECFG )
               WRITE(JPL,2003) ( ECFG(I2), I2 = 1, NECFG )
               if( NVECG > 0 ) then
                  WRITE(JPL,2007)
                  WRITE(JPL,2002) ( VECGRS(I2), I2 = 1, NVECG )
                  WRITE(JPL,2003) ( VECG(I2), I2 = 1, NVECG )
               endif
               WRITE(JPL,2004) ( '----------', I2 = 1, NECFG )
               WRITE(JPL,2005) '#C',KI,R(KI),' : R_cut'
               WRITE(JPL,2006) '#L',1.D0/DX, ' : 1 / logarithmic mesh'
               WRITE(JPL,2006) '#Z',ZV,      ' : core charge         '
               WRITE(JPL,2006) '#W',OCCN(jj,L),' : No. of electrons    '
               WRITE(JPL,2006) '#G',EIGAE(jj,L),' : eigenvalues [Ryd.]  '
               WRITE(JPL,4030)
               WRITE(JPL,4030) '     pseudo-wavefu', &
               (( 'nction for j='//NUM(MJ2) )//'/2  ' )
               WRITE(JPL,4030) '    r [ bohr ]    ', &
               '  major part      ','  minor part      '

               FACT = SIGN( 1.D0, PPNLAE(KI,jj) )
               PPMINR = 0.d0
               DO I = 1, MESH
                  WRITE(JPL,2200) R(I), FACT*PPNLAE(I,jj), &
                  FACT*PPMINR
               enddo
               ENDFILE(JPL)
               CLOSE(JPL)

            endif
         enddo
end if
      !-----------------------------------------------------------------------


      !---------- output all electron wave-functions -------------------------
         FNAME = ( ((( ANAME(NZN)(1:nacr)//'_' )//MS(L+1) )//NUM(MJ2) ) &
         //'.ae' )
         WRITE(J06,'(1X,A20)') FNAME
         OPEN( JPL, FILE = FNAME )
         WRITE(JPL,4000) ANAME(NZN),award,calwrd,'#','( for l =',L
         WRITE(JPL,4001)
         WRITE(JPL,2002) ( ECFGRS(I2), I2 = 1, NECFG )
         WRITE(JPL,2003) ( ECFG(I2), I2 = 1, NECFG )
         if( NVECG > 0 ) then
            WRITE(JPL,2007)
            WRITE(JPL,2002) ( VECGRS(I2), I2 = 1, NVECG )
            WRITE(JPL,2003) ( VECG(I2), I2 = 1, NVECG )
         endif
         WRITE(JPL,2004) ( '----------', I2 = 1, NECFG )
         WRITE(JPL,2005) '#C',KI,R(KI),' : R_cut'
         WRITE(JPL,2006) '#L',1.D0/DX, ' : 1 / logarithmic mesh'
         WRITE(JPL,2006) '#Z',ZV,      ' : core charge         '
         WRITE(JPL,4010) '#N',NREFE,   ' : No. of reference energies'
         WRITE(JPL,4020) ( '#E',I,REFE(I),'[Ryd.]', I = 1, NREFE )
         WRITE(JPL,4030) '   << all electron',' wave functions >>'
         WRITE(JPL,4030) '  r [ bohr ]      ', &
         ( (( 'major part('//NUM(NI) )//') ' ), &
         (( '  minor part('//NUM(NI) )//') ' ), &
         NI = 1, NREFE )
         PMINR = 0.d0
         DO 900 MM = 1, KIMAX
            WRITE(JPL,2200) R(MM), &
            ( PNLAE(MM,I), PMINR, I = 1, NREFE )
         900 END DO
         ENDFILE(JPL)
         CLOSE(JPL)
      !-----------------------------------------------------------------------

      end if
   end do ldo


   KLOCAL = NINT( LOG(RCUTL/CR)/DX )
!         DO MM = 1, MESH
!            VLOC(MM) = the local potential
!         end do

!---------- output local potential -------------------------------------
   FNAME = ( ANAME(NZN)(1:nacr)//'_local' )
   WRITE(J06,'(1X,A20)') FNAME
   OPEN( JPL, FILE = FNAME )
   WRITE(JPL,4000) ANAME(NZN),award,calwrd,'#'
   WRITE(JPL,4001)
   WRITE(JPL,2002) ( ECFGRS(I2), I2 = 1, NECFG )
   WRITE(JPL,2003) ( ECFG(I2), I2 = 1, NECFG )
   if( NVECG > 0 ) then
      WRITE(JPL,2007)
      WRITE(JPL,2002) ( VECGRS(I2), I2 = 1, NVECG )
      WRITE(JPL,2003) ( VECG(I2), I2 = 1, NVECG )
   endif
   WRITE(JPL,2004) ( '----------', I2 = 1, NECFG )
   WRITE(JPL,2005) '#C',KLOCAL,R(KLOCAL),' : R_cut','_local-V'
   WRITE(JPL,2006) '#L',1.D0/DX, ' : 1 / logarithmic mesh'
   WRITE(JPL,2006) '#Z',ZV,      ' : core charge         '
   WRITE(JPL,4030) '  r [ bohr ]      ','  local potential '
   DO 950 MM = 1, MESH
      WRITE(JPL,2200) R(MM), VLOC(MM)
   950 END DO
   ENDFILE(JPL)
   CLOSE(JPL)
!-----------------------------------------------------------------------


   DO MM = 1, MESH
   !            SGV(MM)    = the all-electron valence-electron density
   !            SGC(MM)    = the all-electron core-electron density
      SG(MM)     = ssg(mm)*R(mm)  !the soft part of the USPP valence-electron density
      SGhard(MM) = hsg(mm)*R(mm)  !the hard part of the USPP valence-electron density
   ! cc1(MM)   = ncc1(mm)*R(mm)
   ! tc1(MM)   = ntc1(mm)*R(mm)
   end do

!---------- output valence and core charge density ---------------------
! eck
   CALL INTGBP( MESH, DX, SG,     elsoft )
   CALL INTGBP( MESH, DX, SGhard, elhard )
! ALL INTGBP( MESH, DX, R, ncc1, sncc1 )
! ALL INTGBP( MESH, DX, R, ntc1, sntc1 )
   write(*,*) ' No. of valence electrons (soft part):', elsoft
   write(*,*) ' No. of valence electrons (hard part):', elhard
   write(*,*) ' No. of valence electrons (  total  ):', &
   elhard+elsoft
! rite(3776,*) ' No. of ncc1:', sncc1
! rite(3776,*) ' No. of ntc1:', sntc1
! eck
   FNAME = ( ANAME(NZN)(1:nacr)//'_val' )
   WRITE(J06,'(1X,A20)') FNAME
   OPEN( JPL, FILE = FNAME )
   WRITE(JPL,4000) ANAME(NZN),award,calwrd,'#'
   WRITE(JPL,4001)
   WRITE(JPL,2002) ( ECFGRS(I2), I2 = 1, NECFG )
   WRITE(JPL,2003) ( ECFG(I2), I2 = 1, NECFG )
   if( NVECG > 0 ) then
      WRITE(JPL,2007)
      WRITE(JPL,2002) ( VECGRS(I2), I2 = 1, NVECG )
      WRITE(JPL,2003) ( VECG(I2), I2 = 1, NVECG )
   endif
   WRITE(JPL,2004) ( '----------', I2 = 1, NECFG )
   WRITE(JPL,2006) '#L',1.D0/DX, ' : 1 / logarithmic mesh'
   WRITE(JPL,2006) '#Z',ZV,      ' : core charge         '
   WRITE(JPL,4030) '    rho_v ... vale','nce charge density'
   WRITE(JPL,4030) '    rho_c ... core','    charge density'
   WRITE(JPL,4030) '    rho_vs... vale','nce charge density', &
   ' by pseudo wave fu','nctions           '
   WRITE(JPL,2006) '# ',elsoft,       ' : #el (soft part)     '
   WRITE(JPL,2006) '# ',elhard,       ' : #el (hard part)     '
   WRITE(JPL,2006) '# ',elhard+elsoft,' : #el ( total   )     '
   WRITE(JPL,4030) '  r [ bohr ]      ',' 4*pi*r*r* rho_v  ', &
   ' 4*pi*r*r* rho_c  ',' 4*pi*r*r* rho_vs '
   DO 960 MM = 1, MESH
   !            WRITE(JPL,2222) R(MM), SGV(MM), SGC(MM), SG(MM)+SGhard(MM)
      WRITE(JPL,2222) R(MM), SGV(MM), SGC(MM), ssg(MM)+hsg(MM)
   960 END DO
   ENDFILE(JPL)
   CLOSE(JPL)
   2222 FORMAT( 20D23.15 )
!-----------------------------------------------------------------------





   4000 FORMAT('#',5X,'Ultrasoft pseudopotential             : ', &
   A2/'#',5X,9X,'by ',a12,'/ ',a20,'atomic calculation.'/ &
   A1,5X,9X,A9,I2,' )')
   4001 FORMAT('#',5X,'---- valence electron configuration ', &
   'for constructing pp. ----')
   2007 FORMAT('#',5X,'---- valence electron configuration ')
   2002 FORMAT('#',5X,7(1X,A4,'  ',3X))
   2003 FORMAT('#',5X,7(F7.4,3X))
   2004 FORMAT('#',5X,7(A10))
   2005 FORMAT(A2,I5,' (',D16.10,' )  ',10A8)
   2006 FORMAT(A2,D18.10,9X,A23)
   4010 FORMAT(A2,I5,22X,A28)
   4020 FORMAT((A2,I10,D18.10,2X,A6 ))
   4030 FORMAT('#',20A18)
   2200 FORMAT(20D23.15)





   RETURN

   1000 FORMAT(4D18.10)
   1010 FORMAT( / '01  ( real mesh )')
   1020 FORMAT( / '02  ( r * potential V(r) )')
   1030 FORMAT( / '03  ( 4 * pai * r * r * electron density )')

   END SUBROUTINE ATMLDA



