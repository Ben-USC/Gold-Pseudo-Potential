



   subroutine tran(l,r)
!-----------------------------------------------------------------------
!     transferability
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   dimension     r(msh) , f(msh)
   parameter     ( iii = 5  )
   parameter     ( jjj  = 10 )
   parameter     ( alp  = 1/274.074  )
   dimension     dvdx(msh),abm(msh)
   dimension     a(msh),b(msh),aa(jjj,jjj)
   dimension     y(msh)
   dimension     bb(jjj),aaa(jjj),xx(msh)
   CHARACTER     OFILE*13
!-----------------------------------------------------------------------

!      l=lus
   x=exp(dx)
   cr=rmax/x**mesh

   itr=log(rtra/cr)/dx

!      write(*,*)'check'

   ofile=(ms(l+1)//'_chiae.dat')
   open(9,file=ofile)

   !---- finite difference ------------------------------------------------
      call fd1 (1,itr+60,mesh, v, dx, dvdx)
   !-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!       do e = -4.d0, 2.d0, 0.01d0
   do ie = 0, 600
      e = - 4.d0 + 0.01*dble(ie)
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   !           starts outward integration with a power series.
   !-----------------------------------------------------------------------

      do i=1, itr+60
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
   !-----------------------------------------------------------------------
      call gauss(aa,bb,aaa,iii,jjj)
   !-----------------------------------------------------------------------
      do i=1,5
         y(i) = r(i)**(l+0.5)                             !Y
         xx(i) = (l+0.5)*r(i)**(l+0.5)                    !Y'
         f(i) = 0                                         !Y''
         do j=1,5
            y(i)=y(i)+aaa(j)*r(i)**(l+j+0.5)
            xx(i)=xx(i)+(l+j+0.5)*aaa(j)*r(i)**(l+j+0.5)
         enddo
         f(i)=-a(i)*r(i)*(xx(i)-y(i)*0.5)+((l+0.5d0)**2.d0 &
         +abm(i)*(v(i)-e)*r(i)*r(i))*y(i)
      enddo

   !--- Adams-Bashforth Adams-Moulton -------------------------------------

      do i = 5, itr + 50
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

      i = itr
      chiae  =  ( xx(i) / y(i) - 0.5d0 ) / r(i)

      write(9,*)e, chiae

   end do

   close(9)
   write(*,'(4x,2a, /)') 'Create : ', ofile

   return
   end subroutine tran



   subroutine tran1(  L, R)
!-----------------------------------------------------------------------
!     Generalization eigenvalue equation
!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   parameter     (jjj=10)
   dimension     r(msh)
   dimension     y(msh),ula(msh),usm(5,5),beti(5,msh),wla(msh)
   dimension     wsm(5,msh), vla(msh),vsm(5,msh),vla1(msh)
   dimension     vsm1(5,msh)
   dimension     wtla(5),wtsm(5,5)
   dimension     aa(jjj,jjj),bb(jjj),aaa(jjj) , f(msh)
   CHARACTER     OFILE*13
!-----------------------------------------------------------------------
!      l=lus
   x=exp(dx)
   cr=rmax/x**mesh

   itr=log(rtra/cr)/dx

   ofile=(ms(l+1)//'_chil.dat')
   open(9,file=ofile)

!-----------------------------------------------------------------------
!       do e = -4.d0, 2.d0, 0.01d0
   do ie = 0, 600
      e = - 4.d0 + 0.01*dble(ie)
      if( abs(vloc(1)-e) < 1.d-05) cycle
   !-----------------------------------------------------------------------
      1 continue
   !      write(*,*) 'kaisuu',e
      do i=1, mesh
         ula(i) = (l+0.5d0)**2.d0+(vloc(i)-e)*r(i)*r(i)
         do j=1, nref
            if ( i <= icl) beti(j,i) = r(i)**(1.5d0)*bas(j,i)
            if ( i > icl) beti(j,i) = 0.d0
         end do
      end do

      do j=1, nref
         do k=1, nref
            usm(j,k) = dele(j,k) - e*qch(j,k)
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
   !--- Adams-Bashforth Adams-Moulton -------------------------------------

      do i = 5, mesh-1
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


   !-----wljk(x)-----------------------------------------------------------
      do k = 1, nref
         do i=1,5
            do j=1,5
               aa(i,j)=(ula(i)-(l+(j-1)+0.5d0)*(l+j-1+0.5d0)) &
               *r(i)**(l+j-1+0.5d0)
            end do
            bb(i)=beti(k,i)
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
         !         write(*,*)'f(i)',f(i),vsm1(k,i)

         enddo


      !--- integrates outward to classical turning point ---------------------
      !--- Adams-Bashforth Adams-Moulton -------------------------------------

         do i = 5, mesh-1
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


      do k = 1, nref
         do i = 1, mesh
            if (i >= icl) then
               f(i) = 0.d0
            else
               f(i) = beti(k,i) * wla(i)

            end if
         enddo
         call intgbp( mesh, dx, f, tint1)
         wtla(k) = tint1
      enddo
   !      write(*,*)icl ,ki
      do k = 1, nref
         do i = 1, nref
            do j =1, mesh
               if (j >= icl) then
                  f(j) = 0.d0
               else
                  f(j) = beti(k,j) * wsm(i,j)
               endif
            end do
            call intgbp (mesh, dx, f, tint1)

            wtsm(k,i) = tint1
         enddo
      enddo


      do i = 1, nref
         bb(i) = 0.d0
         do j = 1, nref
            aa(i,j) = 0.d0
            do k = 1, nref
               aa(i,j) = aa(i,j) + usm(i,k) * wtsm(k,j)
            end do
            if ( i == j ) aa(i,j) = aa(i,j) + 1.d0
         !      write(*,*)'aa',aa(i,j)
            bb(i) = bb(i) -usm(i,j) * wtla(j)
         enddo
      !      write(*,*)'bb',bb(i)
      enddo

      call gauss(aa,bb,aaa,nref,jjj)

      do i = 1,mesh
         y(i) = wla(i)
         do j =1, nref
            y(i) = y(i) + aaa(j) * wsm(j,i)
         enddo
      enddo



   !--- LOG - DER at KI ---
      I  =  itr

      dyki = vla(i)
      do j = 1, nref
         dyki = dyki + aaa(j) * vsm(j,i)
      enddo

      chil  =  ( dyki / y(i) - 0.5d0 ) / r(i)

      write(9,*) e,chil
   end do
   close(9)
   write(*,'(4x,2a, /)') 'Create : ', ofile

   RETURN
   end subroutine tran1


