      pure subroutine get2Dintpolpoints(norder, lambda, npkt, errcode)
      use feminterface , only : get1Dintpolpoints
      use femtypes
      implicit none
      integer (I4B) norder,  npkt, errcode
      real (DP), allocatable :: lambda(:,:)
      intent (in) :: norder
      intent (out) :: errcode, npkt, lambda
!
!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 - 2015 Institute for Micro Systems Technology,
!                              Hamburg University of Technology.
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
!------------------------------------------------------------------------------
!
!    $Revision: 1.2 $
!    $Date: 2015/11/11 17:18:26 $
!    $Author: m_kasper $

!  This program returns Gauss-Lobatto interpolation points in a triangle.
!  It assumes norder for now is like in 1D: 
!        norder = p + 1
!  It first maps from 1D coordinates to 2D barycentric coordinates on Edge 3, 
!  then just builds the points up the triangle, so that the last point is 
!  always (0,0,1)
!
!      L3=1
!       /\
!      /  \
!     /____\
!  L1=1      L2=1
!
!
!  Input:
!            norder   order of the rule (number of points and weights) in the range 2 to 22
!  Output:
!            lambda     array of barycentric coordinates of the interpolation points 
!            errcode  =1001 if an integration routine of this order is not available
!  
!
!  References: See get1Dintpolpoints.f90
!
!  local variables
!
      integer (I4B) i, j, k, idx
      real (DP) :: xtab(norder)

! fetch the Gauss-Lobato interpolation points
      call get1Dintpolpoints(norder, xtab, errcode)

! number of interpolation points is complete for triangle of p-th order
      npkt = (norder)*(norder+1)/2
      allocate(lambda(3,npkt))
! initialize all points as zeros
      lambda(:,:) = 0._DP
! For linear element, there are only 3 points
! just list the points--each lambda is 1 at its own vertex. See above triangle.
      if (norder .eq. 2) then
        lambda(1,1) = 1._DP
        lambda(2,2) = 1._DP
        lambda(3,3) = 1._DP
      else
! build points from the 1D points
        idx=1
        do i=1, norder
          k = 1
          do j=i, norder
            lambda(1,idx) =  (1._DP - xtab(j))/2._DP
            lambda(2,idx) =  (1._DP + xtab(k))/2._DP
            lambda(3,idx) =   1._DP - lambda(2,idx) - lambda(1,idx)
            idx = idx+1
            k = k+1
          end do
        end do
      end if

      return
      end subroutine get2Dintpolpoints
