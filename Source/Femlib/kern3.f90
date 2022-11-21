      subroutine kern3(pliste,start,pzahl,xn,yn,xk,yk,pz)
      use feminterface, only: leftpolygon
      use femtypes
      implicit none
      integer (I4B) pliste(:), pzahl, pz ,start
      real (DP) xn(:), yn(:), xk(:), yk(:)
      intent (in) :: pliste, start, pzahl, xn, yn
      intent (out) :: xk, yk, pz
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
!    $Date: 2008/07/30 14:01:14 $
!    $Author: m_kasper $
!
!  compute the kernel of a polygon, i.e. the set of points which is visible 
!  from each of the poygon vertices (and edges). 
!  The compexity of this algorithm is O(h^2), which is far from being optimal.
!  Thus the algorithm is only intendet to be use with a small number of vertices.
!
!  Input:
!            pliste   List of vertices of the poygon: pliste(start:start+pzahl)
!            start    Pointer to the first vertex in plist
!            pzahl    Number of vertices of the polygon
!            xn,yn    Coodinates of the vertices
!
!  Output:
!            xk, yk   Coodinates of the kernel vertices
!            pz       Number of vertices of the kernel polygon
!
!  local variables
      integer (I4B) i, pin, pout
      real (DP) ,allocatable :: xi(:), yi(:), xo(:), yo(:)
!
!  the number of vertices of the kernel may be larger than the number of vertices 
!  of the orginal polygon, but is always smaller than the square of this number
      allocate (xi(pzahl**2), yi(pzahl**2), xo(pzahl**2), yo(pzahl**2))
!
      xk(1:pzahl)=xn(pliste(start:start+pzahl-1))
      yk(1:pzahl)=yn(pliste(start:start+pzahl-1))
      pin=pzahl
      xi(1:pzahl)=xk(1:pzahl)
      yi(1:pzahl)=yk(1:pzahl)
!
      do i=1,pzahl
!  remove the part left to the edge i -> i+1
        call leftpolygon(xk(i),yk(i),xk(mod(i,pzahl)+1),                &
     &    yk(mod(i,pzahl)+1),pin,xi,yi,pout,xo,yo)
!
        if (pout .eq.0 .or. i.eq.pzahl) then
          exit
        end if
        pin=pout
        xi(1:pout)=xo(1:pout)
        yi(1:pout)=yo(1:pout)
      end do
!      
      pz=pout
      xk(1:pz)=xo(1:pz)
      yk(1:pz)=yo(1:pz)
      deallocate (xi, yi, xo, yo)
      return
      end subroutine kern3
!
!
!
      subroutine leftpolygon(xa,ya,xb,yb,pin,xi,yi,pout,xo,yo)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) pin, pout
      real (DP) xa, ya, xb, yb, xi(:), yi(:), xo(:), yo(:)
      intent (in) :: xa, ya, xb, yb, pin, xi, yi
      intent (out) :: pout, xo, yo
!    
!  compute the part of the polygon (xi,yi), which is at the left of the ray (xa,ya)->(xb,yb)
!    
!  Input:
!            xa, ya   coordinates of the first point of the line
!            xb, yb   coordinates of the second point of the line
!            pin      number of vertices of the input polygon
!            xi, yi   vertex coordinates of the input polygon
!  Output:
!            pout     number of vertices of the clipped polygon
!            xo, yo   vertex coordinates of the clipped polygon
!    
!  local variables      
      real (DP) f1, f2, ffirst
      integer (I4B) i, k
      logical pleft, acleft, first
!
!  compute the area (double) of the first triangle
!  area is positive if the point 1 is left from the line (A - B)
      f1=(yb-yi(1))*(xa-xb) - (ya-yb)*(xb-xi(1))
      ffirst=f1
      if (f1 .ge. 0._DP) then
!  first point is left, it should be in the output list
        k=1
        xo(1)=xi(1)
        yo(1)=yi(1)
        pleft=.true.
        first=.true.
      else
!  first point is right
        k=0
        pleft=.false.
        first=.false.
      end if
!
      do i=2,pin
!  treat the edge (i-1) -> (i) of the polygon
!  compute the area (double) of the next triangle
!  area is positive if the point i is left from the line (A - B)
        f2=(yb-yi(i))*(xa-xb) - (ya-yb)*(xb-xi(i))
        if (f2 .ge. 0._DP) then
          acleft=.true.
        else
          acleft=.false.
        end if
        if (acleft .and. pleft) then
!  preceding point was left and actual point is left
          k=k+1
          xo(k)=xi(i)
          yo(k)=yi(i)
        else if (acleft .and. (.not. pleft)) then
!  preceding point was right and actual point is left - two new vertices
          k=k+2
!  compute coordinates of the intersection by using the 
!  preceding and actual coordinate weighted with the areas
          xo(k-1)=(xi(i)*f1-xi(i-1)*f2)/(f1-f2)
          yo(k-1)=(yi(i)*f1-yi(i-1)*f2)/(f1-f2)
          xo(k)=xi(i)
          yo(k)=yi(i)
          pleft=.true.
        else if ((.not. acleft) .and. pleft) then
!  preceding point was left and actual point is right
          k=k+1
!  compute coordinates of the intersection by using the 
!  preceding and actual coordinate weighted with the areas
          xo(k)=(xi(i)*f1-xi(i-1)*f2)/(f1-f2)
          yo(k)=(yi(i)*f1-yi(i-1)*f2)/(f1-f2)
          pleft=.false.
!        else
!  preceding point was right and actual point is right:  nothing to do          
        end if
        f1=f2
      end do
!
!  treat the last edge (n) -> (1) of the polygon      
      if ( (first .and. (.not. pleft)) .or.                             &
     &  ((.not. first) .and. pleft) ) then
!  last point was right and first point is left
!  or last point was left and first point is right
        k=k+1
!  compute coordinate of the intersection by using the 
!  preceding and actual coordinate weighted with the areas
        xo(k)=(xi(1)*f1-xi(pin)*ffirst)/(f1-ffirst)
        yo(k)=(yi(1)*f1-yi(pin)*ffirst)/(f1-ffirst)
      end if
!
      pout=k
      return
      end subroutine leftpolygon
