      subroutine decroachn(encrlimit,x3,y3,splitted)
      use feminterface, only: lswap, encrtest, splitencrelement
      use globalvariables
      use femtypes
      implicit none
      real (DP) encrlimit, x3, y3
      logical splitted
      intent (in) :: encrlimit, x3, y3
      intent (out) :: splitted
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
!    $Revision: 1.8 $
!    $Date: 2008/12/22 15:03:35 $
!    $Author: m_kasper $
!
!  A new vertex is tested whether it encroaches upon a segment if so, 
!  the segment is splitted
!
!  Input:    encrlimit  angle at which a segment is considered to be encroached
!            x3, y3     node location which it to be tested for encroachment
!  Output:
!            splitted   = .true. if the point (x3,y3) was encroached upon at least one segment 
!                                these segments had been splitted 
!            n          total number of elements
!            p          total number of nodes
!            e          element information, nodes of the elements
!            en         neighborhood information, neighbors of the elements
!            geb        assignment of elements to the geometry regions
!            kzi        node branch information
!            xn, yn     coordinates of the nodes
!            x          solution vector
!
!  local variables
      integer (I4B) nn, i, j, swaps
      integer (I4B), allocatable :: slist(:)
      logical encroached,split
!
      allocate(slist(2*n))
      slist(1:2*n)=0
      splitted=.false.
!  de-encroach segments
      nn=n
      do i=1,nn
        do j=1,3
!  if a neighbor is present, the check will be done later 
          if (en(j,i) .gt. i) cycle
!  is this edge an input segment and the vertex encroached upon ?
          call encrtest(i,j,x3,y3,encrlimit,e,en,geb,kzi,xn,yn,encroached)
          if (encroached) then
!  if so split the segment
            call splitencrelement(i,j,slist,split)
            splitted=splitted.or.split
          end if
        end do
      end do
      if (splitted) then
!  if a segment was splitted restore the delaunay triangulation
        call lswap(n,e,xn,yn,en,geb,0,.false.,swaps,2,slist)
      end if
      deallocate(slist)
      return
      end subroutine decroachn
