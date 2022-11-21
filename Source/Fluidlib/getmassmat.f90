subroutine getmassmat
use femtypes
use feminterface
use globalvariables
use fluidvariables
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
!*******************************************************************************
!
!     Compute lumped mass matrix
!
!*******************************************************************************
! local variables
!*******************************************************************************
integer(I4B) :: ie, j, jj, k
integer(I4B) :: npkt, polylo, polyhi
integer(I4B) :: intorder, polyorder, nff
integer(I4B) :: err, errcode
real(DP) :: aw, areael
real(DP), allocatable :: weight(:), lambda(:,:), xsi(:), gxsi(:,:)
!real(DP), pointer :: massmat(:)
!*******************************************************************************
allocate(dmmat(ndof))
dmmat = 0._DP

do ie = 1, n
   polyorder = ep(ie,1)
   polyhi=polyorder
   polylo=1
   nff = (polyorder+1)*(polyorder+2)/2

   allocate(xsi(nff), gxsi(nff,2) )
   areael = areaf(ie)
   intorder=2*polyorder
!  fetch numerical integration points (Gauss points)
   call get2Dintegpoints(intorder, npkt, weight, lambda, err)

   do k=1,npkt
!  get shape function and gradients at location of integration points
      call shapefunction(lambda(:,k),xn(e(:,ie)),yn(e(:,ie)),     &
                         polylo,polyhi,nff,.true.,xsi,gxsi,errcode)
      aw = areael*weight(k)
!  for all shape functions
      do j = 1,nff
         jj = eg(ie,1)%d(j)
         dmmat(jj) = dmmat(jj) + aw*xsi(j)
      end do ! j
   end do ! k
end do ! ie

deallocate(weight,lambda, xsi, gxsi)

dmmat(:) = 1.0_DP/dmmat(:)

! check mass matrix
!allocate(massmat(ndof))
!massmat = 0._DP
!do ie = 1,n
!   nff = 3
!   areael = areaf(ie)
!   do j = 1,nff
!      jj = eg(ie)%d(j)
!      massmat(jj) = massmat(jj) + areael/3.0_DP
!   end do !j
!end do ! ie
!massmat = 1.0_DP/massmat

end subroutine getmassmat