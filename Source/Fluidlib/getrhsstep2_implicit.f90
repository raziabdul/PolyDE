subroutine getrhsstep2_implicit(elem,polyorder,b,nff)
use feminterface, only: get2Dintegpoints, shapefunction
use femtypes
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer (I4B) :: elem
integer (I4B) :: polyorder, nff
real(DP), pointer :: b(:)   ! intent(out)
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
!    $Revision: 1.5 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!*****************************************************************************************
!
!  RHS vector for implicit pressure in step 2
!
!  13.10.05
!
! to be done: check how to impose dirichlet b.c. when a and b are seperately created
!             a is created in sub. pstiff while b in this sub.
!
! remarks: use with elemmatrix_pstiff, assemblyrhsstep2
!*****************************************************************************************
!  local variables
!*****************************************************************************************
integer(I4B) :: i, j, jj, k, polylo, polyhi, npkt
integer(I4B) :: intorder
integer(I4B) :: errcode, err
real(DP) :: aw
real(DP), allocatable :: weight(:), lambda(:,:), xsi(:), gxsi(:,:)
real(DP), pointer :: elvel(:,:), velterm(:)
!*****************************************************************************************
polyhi=polyorder
polylo=1
!  size of the element matrix
nff=(polyorder+1)*(polyorder+2)/2

allocate(elvel(2,nff), velterm(nff))
allocate(b(nff), xsi(nff), gxsi(nff,2) )

b = 0.0_DP

do j = 1,nff
   jj = eg(elem,1)%d(j)
   elvel(1,j) = theta(1)*unkno(2,jj)
   elvel(2,j) = theta(1)*unkno(3,jj)     
end do ! j

intorder=2*polyorder
!  fetch numerical integration points (Gauss points)
call get2Dintegpoints(intorder, npkt, weight, lambda, err)

do k=1,npkt
!  get shape function and gradients at location of integration points
   call shapefunction(lambda(:,k),xn(e(:,elem)),yn(e(:,elem)),polylo,polyhi,nff,.true.,&
                      xsi,gxsi,errcode)
   aw=areaf(elem)*weight(k)
!  for all shape functions
   do j=1,nff
      velterm(j) = ( gxsi(j,1)*elvel(1,j) + gxsi(j,2)*elvel(2,j) )
      do i=1,nff
         b(i) = b(i) - aw*xsi(i)*velterm(j)
      end do ! i
   end do ! j
end do ! k

deallocate(weight,lambda, gxsi, xsi, velterm, elvel)

end subroutine getrhsstep2_implicit