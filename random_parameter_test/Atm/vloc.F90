

SUBROUTINE vlocal(r)
!---- Ultrasoft pseudopotential ----------------------------------------
!  in  put :
!         rloc  ...... cut off distance
!         v     ...... all electron potential
!  out put :
!         vloc  ...... local potential

!-----------------------------------------------------------------------
   use common_variables
   IMPLICIT REAL*8 (A-H,O-Z)
   parameter     ( jjj  = 10 )
   dimension     r(msh)
   dimension     b(msh),aa(jjj,jjj)
   dimension     v1(msh),v2(msh),v3(msh),v4(msh)
   dimension     a2n(jjj)
   dimension     dvloc1(msh)
   dimension     dvloc2(msh)
   dimension     dvloc3(msh)
   dimension     dvloc4(msh)
   character(13) :: fname
!-----------------------------------------------------------------------

   x=exp(dx)
   cr=rmax/x**mesh
   nloc=log(rloc/cr)/dx
   rloc=r(nloc)

!      do i=1,msh
!         v(i)=vus(i)
!      end do


!--- derivative of V,  finite difference  -----------------------
   call fd1 (nloc, nloc+1, msh, v, dx, v1)
   call fd2 (nloc, nloc+1, msh, v, dx, v2)
   call fd3 (nloc-6, nloc+7, msh, v, dx, v3)
   call fd4 (nloc-6, nloc+7, msh, v, dx, v4)
!-----------------------------------------------------------------------
   dv0=v(nloc)
   dv1=v1(nloc)/rloc
   dv2=(v2(nloc)-v1(nloc))/rloc/rloc
   dv3=(v3(nloc)-3*v2(nloc)+2*v1(nloc))/rloc**3
   dv4=(v4(nloc)-6*v3(nloc)+11*v2(nloc)-6*v1(nloc))/rloc**4

   do i = 2, 6
      aa(1,i-1)=rloc**(2*i)
      aa(2,i-1)=2*i*rloc**(2*i-1)
      aa(3,i-1)=2*i*(2*i-1)*rloc**(2*i-2)
      aa(4,i-1)=2*i*(2*i-1)*(2*i-2)*rloc**(2*i-3)
      aa(5,i-1)=2*i*(2*i-1)*(2*i-2)*(2*i-3)*rloc**(2*i-4)
   end do

   b(1)=log(dv0/v0)
   b(2)=dv1/dv0
   b(3)=dv2/dv0-(dv1/dv0)**2
   b(4)=dv3/dv0+2*(dv1/dv0)**3-3*(dv1/dv0)*(dv2/dv0)
   b(5)=dv4/dv0-6*(dv1/dv0)**4+12*(dv1/dv0)**2*(dv2/dv0) &
   -4*(dv1/dv0)*(dv3/dv0)-3*(dv2/dv0)**2


!-----------------------------------------------------------------------
   call gauss(aa,b,a2n,5,jjj)
!-----------------------------------------------------------------------


   do i=1, mesh
      if( i <= nloc) then
         vloc(i)=v0
         do j=2, 6
            vloc(i)=vloc(i)*exp(a2n(j-1)*r(i)**(2*j))
         end do
      else
         vloc(i)=v(i)
      end if
   end do


   call fd1(nloc,nloc+1,msh,vloc,dx,dvloc1)
   call fd2(nloc,nloc+1,msh,vloc,dx,dvloc2)
   call fd3(nloc-6,nloc+7,msh,vloc,dx,dvloc3)
   call fd4(nloc-6,nloc+7,msh,vloc,dx,dvloc4)


   write(*,'(/, a)') '  * Local pseudopotential'
   write(*,'(20x,a)') 'V_loc             V_AE              diff.'
   write(*,'(4x,a,3es18.10)')  &
   &     '0th derivs. : ',   vloc(nloc),  v(nloc), abs(  vloc(nloc)- v(nloc))
   write(*,'(4x,a,3es18.10)')  &
   &     '1st derivs. : ', dvloc1(nloc), v1(nloc), abs(dvloc1(nloc)-v1(nloc))
   write(*,'(4x,a,3es18.10)')  &
   &     '2nd derivs. : ', dvloc2(nloc), v2(nloc), abs(dvloc2(nloc)-v2(nloc))
   write(*,'(4x,a,3es18.10)')  &
   &     '3rd derivs. : ', dvloc3(nloc), v3(nloc), abs(dvloc3(nloc)-v3(nloc))
   write(*,'(4x,a,3es18.10, /)')  &
   &     '4th derivs. : ', dvloc4(nloc), v4(nloc), abs(dvloc4(nloc)-v4(nloc))

   fname = 'vloc.dat'
   open(9, file=fname)
   write(9,'(3es18.10)') ( r(i), vloc(i), v(i), i = 1, mesh )
   close(9)
   write(*,'(4x,2a)') 'Create : ', fname
   write(*,'(78a)') ( '-', i = 1, 78 )


   return
   end SUBROUTINE vlocal


