      subroutine centerofgravity(x,y,n,xctr,yctr)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) x(:), y(:), xctr, yctr
      integer (I4B) n
      intent (in) :: x, y, n
      intent (out) :: xctr, yctr
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
!    $Date: 2006/07/07 21:08:35 $
!    $Author: r_abdul $
!
!  compute the center of gravity of a polygon
!    
!  Input:
!            x, y     vertex coordinates of the polygon
!            n        number of vertices of the polygon
!  Output:
!            xctr     coordinates of the center of gravity
!            yctr
!    
!  local variables      
      real (DP) areasum, sumx, sumy, dxn, dyn, dxo, dyo, twoarea, xs, ys
      integer (I4B) i
!
!  The center of cravity is computed by:
!  1: triangulating the polygon,
!  2: computing the area and center of gravity (midmoint) of each of the triangles 
!  3: taking the (area) weighted sum of the triangle midpoints
!
      areasum=0._DP
      sumx=0._DP
      sumy=0._DP
      dxn=x(2)-x(1)
      dyn=y(2)-y(1)
!  triangulation is done with all triangles hinged on vertex 1, i.e. the triangles (1, i, i+1)
!  this means that some of the triangle may have negative area
!  in total we have n-2 triangles
      do i=2,n-1
        dxo=dxn
        dyo=dyn
        dxn=x(i+1)-x(1)
        dyn=y(i+1)-y(1)
!  twice the triangle area 
        twoarea=dxo*dyn-dyo*dxn
        areasum=areasum+twoarea
!  three times the midpoint of the triange (x1,y1 is added later)
        xs=x(i)+x(i+1)
        ys=y(i)+y(i+1)
!  weighted sum
        sumx=sumx+xs*twoarea
        sumy=sumy+ys*twoarea
      end do
      if (areasum .gt. 0._DP) then
        xctr=sumx/(3._DP*areasum)+x(1)/3._DP
        yctr=sumy/(3._DP*areasum)+y(1)/3._DP
      else
        print*,' Error in centerofgravity'
      end if
      return
      end subroutine centerofgravity
