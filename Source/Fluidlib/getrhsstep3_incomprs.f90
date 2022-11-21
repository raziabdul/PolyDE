subroutine getrhsstep3_incomprs(rhs)
use femtypes
use feminterface
use globalvariables
use fluidvariables
use fluidinterface
implicit none
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
!***********************************************************************************
!
!  Determine RHS of step 3 (velocity correction)
!
!***********************************************************************************
real(DP), pointer :: rhs(:,:) ! intent(out)
!***********************************************************************************
! local variables
!***********************************************************************************
integer(I4B) :: i, ie, ivar, j, jj, k
!integer(I4B) :: iedge, nr, nl, npkt1D
integer(I4B) :: npkt, polylo, polyhi
integer(I4B) :: intorder, polyorder, neldof
integer(I4B) :: err, errcode
integer(I4B), dimension(3), parameter :: inach=(/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg=(/3,1,2/)
real(DP) :: aw, area
real(DP) :: dpdx, dpdy
real(DP), allocatable :: weight(:), lambda(:,:), xsi(:), gxsi(:,:)
real(DP), pointer :: rhsstep3(:,:)
real(DP), pointer :: elprs(:)
!***********************************************************************************

do ie = 1,n

   polyorder = ep(ie,1)
   polyhi = polyorder
   polylo = 1
   neldof = (polyorder+1)*(polyorder+2)/2

   allocate(rhsstep3(2,neldof))
   allocate(elprs(neldof))

   area = areaf(ie)
   do j = 1,neldof
      jj = eg(ie,1)%d(j)
      ! unkno(1,:) => pressure
      elprs(j) = (1.0_DP-theta(2))*unkn1(1,jj) + theta(2)*unkno(1,jj)
   end do

   allocate(xsi(neldof), gxsi(neldof,2) )
   intorder=2*polyorder
   call get2Dintegpoints(intorder, npkt, weight, lambda, err)

   rhsstep3 = 0._DP
   dpdx = 0._DP
   dpdy = 0._DP

   do k = 1,npkt
      call shapefunction(lambda(:,k),xn(e(:,ie)),yn(e(:,ie)),     &
                         polylo,polyhi,neldof,.true.,xsi,gxsi,errcode)
      aw = area*weight(k)
      do j = 1,neldof
         dpdx = gxsi(j,1)*elprs(j)
         dpdy = gxsi(j,2)*elprs(j)
         do i = 1,neldof           
            rhsstep3(1,i) = rhsstep3(1,i) - aw*xsi(i)*dpdx
            rhsstep3(2,i) = rhsstep3(2,i) - aw*xsi(i)*dpdy
         end do ! i
      end do ! j
   end do !k
   deallocate(weight,lambda)
   deallocate(xsi,gxsi)

   do j = 1, neldof
      jj = eg(ie,1)%d(j)
      do ivar = 1,2
         rhs(ivar+1,jj) = rhs(ivar+1,jj) + rhsstep3(ivar,j)
      end do ! ivar
   end do ! j

   deallocate(elprs,rhsstep3)

end do

end subroutine getrhsstep3_incomprs
