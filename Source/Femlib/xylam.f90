      pure subroutine xy2lam(xt,yt,elem,lambda,xn,yn,e)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) elem, e(:,:)
      real (DP) xt, yt, lambda(3), xn(:), yn(:)
      intent (in) :: xt, yt, xn, yn, elem, e
      intent (out) :: lambda
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
!    $Revision: 1.4 $
!    $Date: 2010/10/25 11:56:58 $
!    $Author: chvokas $
!
!  conversion of x,y coordinates in element elem to barycentric coordinates lambda
!
!  Input:
!            xt, yt   world coordinates of the point 
!            xn,yn    coordinates of triangle nodes
!            e        element information (nodes of the element)
!            elem     element number with respect to which coordinates
!                     are computed
!  Output:
!            lambda   barycentric coordinates of point x,y 
!  local variables
      real (DP) dx1, dy1, dx2, dy2, fl2
      real (DP) dtx1, dty1, dtx2, dty2, dtx3, dty3
!     
      dx1=yn(e(2,elem))-yn(e(3,elem))
      dx2=yn(e(3,elem))-yn(e(1,elem))
      dy1=xn(e(2,elem))-xn(e(3,elem))
      dy2=xn(e(3,elem))-xn(e(1,elem))
      fl2=dx2*dy1-dx1*dy2
      dtx1=xn(e(1,elem))-xt
      dtx2=xn(e(2,elem))-xt
      dtx3=xn(e(3,elem))-xt
      dty1=yn(e(1,elem))-yt
      dty2=yn(e(2,elem))-yt
      dty3=yn(e(3,elem))-yt
      lambda(1)=(dtx2*dty3-dtx3*dty2)/fl2
      lambda(2)=(dtx3*dty1-dtx1*dty3)/fl2
      lambda(3)=1._DP-lambda(1)-lambda(2)
!      lambda(3)=(dtx1*dty2-dtx2*dty1)/fl2
      return
      end subroutine xy2lam


      pure subroutine lam2xy(lambda,elem,xt,yt,xn,yn,e)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) elem, e(:,:)
      real (DP) xt, yt, lambda(3), xn(:), yn(:)
      intent (in) :: lambda, xn, yn, elem, e
      intent (out) :: xt, yt
!
!  conversion of barycentric coordinates lambda to x,y coordinates in element elem  
!
!  Input:
!            lambda   barycentric coordinates of point x,y 
!            xn,yn    coordinates of triangle nodes
!            e        element information (nodes of the element)
!            elem     element number with respect to which coordinates
!                     are computed
!  Output:
!            xt, yt   world coordinates of the point 
!     
      xt= dot_product( xn(e(:,elem)) , lambda(:) )
      yt= dot_product( yn(e(:,elem)) , lambda(:) )
      return
      end subroutine lam2xy
      
