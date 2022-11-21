      subroutine schnit(x1,y1,x2,y2,x3,y3,x4,y4,t1,t2,ok)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) x1, y1, x2, y2, x3, y3, x4, y4, t1, t2
      logical ok
      intent (in) :: x1, y1, x2, y2, x3, y3, x4, y4
      intent (out) :: t1, t2, ok
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
!    $Date: 2006/07/07 21:08:36 $
!    $Author: r_abdul $
!
!  Computation of the intersection point of two straight lines
!  passing through the points 1,2 and 3,4 respectively
!
!    Intersection point = x = x1 + t1 * (x2 - x1)
!                       = y = y1 + t1 * (y2 - y1)
!                       = x = x3 + t2 * (x4 - x3)
!                       = y = y3 + t2 * (y4 - y3)
!
!   For t1 in [0,1] the intersection is on the line segment 1-2
!   For t2 in [0,1] the intersection is on the line segment 3-4
!
!  local variables
!
      real (DP) dx21, dx31, dx43, dy21, dy31, dy43, del
      dx21=x2-x1
      dy21=y2-y1
      dx43=x4-x3
      dy43=y4-y3
      del=dx21*dy43-dy21*dx43
      if (abs(del).le.tiny(1._DP)) then
!  One of the segments is a point (has no length)
!  or the lines are parallel
!      print*,'***',' intersection in schnit cannot be computed'
        t1=huge(1._DP)
        t2=huge(1._DP)
        ok=.false.
      else
        dx31=x3-x1
        dy31=y3-y1
        t1=(dx31*dy43-dy31*dx43)/del
        t2=(dx31*dy21-dy31*dx21)/del
        ok=.true.
      end if
      return
      end subroutine schnit