      subroutine encrtest(i,j,x3,y3,encrlimit,e,en,geb,kzi,xn,yn,encroached)
      use feminterface, only: wink1
      use femtypes
      implicit none
      integer (I4B) i, j, e(:,:),  en(:,:), geb(:), kzi(:)
      real (DP) encrlimit, xn(:), yn(:), x3, y3
      logical encroached
      intent (in) :: i, j, x3, y3, encrlimit, e, en, geb, kzi, xn, yn
      intent (out) :: encroached
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
!    $Revision: 1.7 $
!    $Date: 2008/12/22 15:01:48 $
!    $Author: m_kasper $
!
!  determine whether the edge j of element i is a segment i.e. is part of an input branch
!  and whether the point (x3,y3) is encroached on that segment, this is the case if the 
!  angle at which the segment is seen form point (x3,y3) is larger than encrlimit
!
!  Input:  
!            i        element
!            j        vertex
!            x3, y3   node location which it to be tested for encroachment
!          encrlimit  angle at which a point is considered to encroach upon a segment
!            e        element information, nodes of the elements
!            en       neighborhood information, neighbors of the elements
!            geb      assignment of elements to the geometry regions
!            kzi      node branch information
!            zki      branch information, geometry nodes which determine the branches
!            xn, yn   coordinates of the nodes
!  Output:
!            encroached = .true. if the point (x3,y3) is encroached upon the segment n1->n2
!                            
!  local variables
      integer (I4B) n1, n2
      integer (I4B) isucc(3), ipred(3)
      real (DP) wk
      parameter (isucc=(/2,3,1/))  !  successor
      parameter (ipred=(/3,1,2/))  !  predecessor
!
      encroached=.false.
!  it cannot be a segment if both elements are in the same region
      if (en(j,i) .ne. 0) then
        if (geb(i) .eq. geb(en(j,i)) ) return
      end if
      n1=e(isucc(j),i)
      n2=e(ipred(j),i)
!
      wk=wink1(x3,y3,xn(n1),yn(n1),xn(n2),yn(n2))
      if (abs(wk).gt.encrlimit) then
        encroached=.true.
!  segment is encroached 
      end if
      return
      end subroutine encrtest
