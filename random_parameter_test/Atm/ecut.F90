

   subroutine kec( L, R)
!-----------------------------------------------------------------------
!     Kinetic energy criterion
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   parameter     (nq=2000)
   dimension     pbar(0:nq)
   dimension     r(msh) , y(msh) ,f(msh)
   CHARACTER     OFILE*13
!-----------------------------------------------------------------------

   qmax = sqrt( 200.d0 )
   qdel = qmax / nq
   pi = acos( -1.d0 )
   h = rmax / mesh

!---- Pbar -------------------------------------------------------------

   k = 1
   do j = 1, mesh-1
      rj = h * j
      rk = r(k)
   !            write(14,*)'rj,rk', rj , rk

      do while( ( rj - rk ) >= 0.d0 )
         k = k + 1
         rk = r(k)
      !               write(*,*) rj , rk
         if (k==2002) then
            write(*,*) rj , rk
            stop
         end if
      end do

      if ( rj == rk ) then
         y(j)= pusnl(k)

      else
         a = ( pusnl(k)-pusnl(k-1) )/ ( r(k)-r(k-1) )
         b = pusnl(k) - a*r(k)
         y(j) = a*rj + b
      end if
   !            write(20,*)rj,y(j)
   !            write(21,*)r(j),pusnl(j)

   end do

   y(mesh) = pusnl(mesh)

   do i = 0, nq

      do j = 1, mesh
         the = i * qdel * h*j
      !            if (i==nq) write(22,*)h*j, Y(j)
      !            if (i==nq) write(*,*)the, bess(l,the)
         f(j) = y(j) * bess(l,the) * the
      !            if (i==nq) write(*,*) h*j,bess(l,the),the, f(j)
      ! heck
      !            if( l==1 .and. i==nq ) then
      !               write(95,*) r(j), f(j)
      !            end if
      end do


   ! eck
      if( l==1 .AND. i==nq ) then
         do j = 1 ,mesh
         !               write(13,*) h*j, f(j)
         end do
      end if



   !         call intgbp (mesh, dx, r, f, t)
      t = 0.d0
      do j =2, mesh-1,2
         t = t + 2.d0 * f(j) + f(j+1)
      end do
      t = f(1) - f(mesh) + 2.d0 *  t
      t = t * h / 3.d0

   !         write(*,*) t,y(5),pusnl(5)
      pbar(i) = sqrt( 2.d0 / pi ) * t
   end do

!---- Ekin -------------------------------------------------------------
   ekin = 0.d0

   do i = 1, nq
      ekin = ekin + ( pbar(i)*( i * qdel ) ) ** 2.d0
   end do

   ekin = ekin * qdel

! heck
!      do i = 0, nq
!         write(90+l,*) i*qdel, pbar(i)
!      end do
!---- del Ekin ---------------------------------------------------------
   ofile=(ms(l+1)//'_delE.dat')
   open(9,file=ofile)
   do n = 1, nq
      ecut = ( n * qdel ) **2.d0

      de = 0.d0
      do i = 0, n
         de = de + ( pbar( i ) * ( i * qdel ) ) ** 2.d0
      end do

   !         write(*,*) ecut, de
      de = ekin - qdel * de
      write(9,*) ecut, de


   end do
   write(*,'(4x,2a, /)') 'Create : ', ofile
   close( 9 )

   return
   end subroutine kec

   subroutine afq(l,r)
!-----------------------------------------------------------------------
!     augmentation functions Qbar for Projector augmented wave
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   parameter     (nq=1000)
   dimension     qbar(0:nq)
   dimension     r(msh) , f(msh)
   CHARACTER     OFILE*13
!---- qbar -------------------------------------------------------------
   qmax = sqrt( 300.d0 )
   delq = qmax / nq

   ofile=(ms(l+1)//'_Qbar.dat')
   open(9,file=ofile)

   do lq = 0, 2*l,2
      do i = 0, nq
         do j = 1, mesh
            the = i * delq * r(j)
            f(j) = qau2(lq,1,1,j) * bess(lq,the) * r(j)
         !            f(j) = qau(1,1,j) * bess(lq,the) * r(j)
         end do
         call intgbp( mesh, dx, f, t)
         qbar(i) = (i * delq) **2.d0 * t
         write(9,*)(delq * i)**2,qbar(i)
      end do
   end do

   close(9)

   return
   end subroutine afq

subroutine afq2( l, r, iflag )
!-----------------------------------------------------------------------
!     augmentation functions Qbar
!-----------------------------------------------------------------------
!  input :
!         iflag ...... 1 --> USPP, else --> PAW
!
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   parameter     (nq=1000)
   dimension     qbar(0:mxl,0:nq)
   dimension     r(msh) , f(msh)
   CHARACTER     OFILE*13
!---- qbar -------------------------------------------------------------
   qmax = sqrt( 300.d0 )
   delq = qmax / dble(nq)

   ofile = (ms(l+1)//'_Qbar.dat')
   open(9,file=ofile)

   if( iflag == 1 ) then
!----------------------------------------------------------------------!
!       Ultrasoft pseudopotential                                      !
!----------------------------------------------------------------------!
      do lq = 0, 2*l, 2
         do i = 0, nq
            do j = 1, mesh
               the = dble(i) * delq * r(j)
               f(j) = qau(1,1,j) * bess(lq,the) * r(j)
            end do
            call intgbp( mesh, dx, f, t)
            ecut = delq*dble(i) * delq*dble(i)
            qbar(lq/2,i) = ecut * t
            if( lq == 2*l ) write(9,'(10es18.10)') ecut, ( qbar(lll,i), lll = 0, l )
         end do
      end do
   else
!----------------------------------------------------------------------!
!       Projector augmented wave method                                !
!----------------------------------------------------------------------!
      do lq = 0, 2*l, 2
         do i = 0, nq
            do j = 1, mesh
               the = dble(i) * delq * r(j)
               f(j) = qau2(lq,1,1,j) * bess(lq,the) * r(j)
            end do
            call intgbp( mesh, dx, f, t)
            ecut = delq*dble(i) * delq*dble(i)
            qbar(lq/2,i) = ecut * t
            if( lq == 2*l ) write(9,'(10es18.10)') ecut, ( qbar(lll,i), lll = 0, l )
         end do
      end do
   end if
   close(9)
   write(*,'(4x,2a, /)') 'Create : ', ofile

   return
   end subroutine afq2

