subroutine elementmatrix_pstiff(elem,polyorder,a,nff)
use feminterface, only: pdecoeff, get2Dintegpoints, shapefunction
use femtypes
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer (I4B) :: elem
integer (I4B) :: polyorder, nff
real(DP), pointer :: a(:,:)  ! intent(out)
intent (in) :: elem, polyorder
intent (out) :: nff
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 Institute for Micro Systems Technology,
!                       Hamburg University of Technology.
!
!    This file is part of PolyDE.
!
!    PolyDE is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by the Free Software Foundation; either version 2, or (at your
!    option) any later version.
!
!    PolyDE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!    MA 02110-1301 USA.
!
!    $Revision: 1.6 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!*******************************************************************************
!
!  Stiffness matrix for implicit pressure in step 2
!
!*******************************************************************************
!  K = theta(1)*theta(2)*deltp(inode)*(grad N * grad N)
!
!  K*pdiff = rhs
!
!  pdiff = p(t(n+1))-p(t(n))
!
!  Remark: deltp(inode) is brought to rhs
!*******************************************************************************
!  local variables
!*******************************************************************************
integer(I4B) :: i, j, k, polylo, polyhi, err, npkt
integer(I4B) :: intorder
integer(I4B) :: ii, jj, ivertex, brvmn, zrbpres, ipbc
integer(I4B) :: errcode
integer(I4B), dimension(3), parameter :: inach=(/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg=(/3,1,2/)
real(DP) :: lam(3)
real(DP) :: aw, xs, ys
real(DP) :: theta12
real(DP), allocatable :: weight(:), lambda(:,:), xsi(:), gxsi(:,:)
real(DP), pointer :: b(:)
complex(DPC) :: value
logical, pointer :: dbc(:)
!*******************************************************************************
theta12 = theta(1)*theta(2)  ! appears in pstiff.f90

polyhi=polyorder
polylo=1
!  size of the element matrix
nff=(polyorder+1)*(polyorder+2)/2

allocate(a(nff,nff), b(nff), xsi(nff), gxsi(nff,2) )
a = 0._DP
b = 0._DP

intorder=2*polyorder
!  fetch numerical integration points (Gauss points)
call get2Dintegpoints(intorder, npkt, weight, lambda, err)

do k=1,npkt
!  get shape function and gradients at location of integration points
   call shapefunction(lambda(:,k),xn(e(:,elem)),yn(e(:,elem)),polylo,polyhi,nff,&
                      .true.,xsi,gxsi,errcode)

   aw=theta12*areaf(elem)*weight(k)

   do j=1,nff
      do i=1,nff
         a(j,i) = a(j,i) + aw*(gxsi(j,1)*gxsi(i,1) + gxsi(j,2)*gxsi(i,2))
      end do ! i      
   end do ! j
end do

deallocate(weight,lambda, xsi, gxsi)
!*******************************************************************************
! Dirichlet boundary condition for pressure
!*******************************************************************************
allocate(dbc(nff))
dbc(1:nff)=.false.

do ivertex = 1,3

   if (.NOT. dirichletbc(elem)%bcdof(ivertex)%varmn(1)) cycle
   brvmn = dirichletbc(elem)%bcdof(ivertex)%brmn(1)
   zrbpres = zrbmn(brvmn,1)

! only need dirichlet bc for pressure
   if (zrbpres .EQ. 0) then
      xs= xn(e(ivertex,elem))
      ys= yn(e(ivertex,elem))
      
      call getbcval2D_mn(1,brvmn,xs,ys,value)
            
      lam(ivertex)=1._DP
      lam(inach(ivertex))=0._DP
      lam(ivorg(ivertex))=0._DP

      call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),1,1,nff,.false.,xsi,&
                         gxsi,errcode)
      b(ivertex)  = real(value,DP)/xsi(ivertex)
      dbc(ivertex)= .true.

   end if !zrbpres

   if (npbc .GT. 0) then
      jj = eg(elem,1)%d(ivertex)
      do ipbc = 1,npbc
         if (-(kzi(jj)) .NE. pkp(ipbc)) cycle
         ! want to have only p = 0 at a fixed point.
         if (zrbpbc(ipbc) .NE. 0) then
            write(*,*) 'Warning: fixed value of pressure was not imposed'
            cycle
         end if    
         value = varpbc(npbc)
         ! only p = 0 needs to be imposed.
         b(ivertex) = value   
         dbc(ivertex) = .true.
      end do ! ipbc
   end if  
end do ! ivertex
!*******************************************************************************
do i=1,nff
   if (dbc(i)) then
      a(i,1:nff)= 0._DP
      a(i,i)    = 1._DP
   else
      ii = eg(elem,1)%d(i)
      do j=1,nff
         if (dbc(j)) then
            bglobal(ii) = bglobal(ii) - a(i,j)*b(j)
            a(i,j) = 0._DP
         end if
      end do
   end if
end do
!*******************************************************************************
deallocate(xsi)
end subroutine elementmatrix_pstiff