      subroutine decroach(encrlimit, meshsiz2)
      use feminterface, only: lswap, encrtest, splitencrelement
      use globalvariables
      use femtypes
      implicit none
      real (DP) encrlimit, meshsiz2(:)
      intent (in) :: encrlimit, meshsiz2
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
!    $Date: 2008/12/22 15:04:11 $
!    $Author: m_kasper $
!
!  Remove vertices encroached upon a segment by subdividing the segment
!  A vertex is call encroached upon a segment if the angle under which the segment 
!  is seen is larger than a certain limit (encrlimit). Such a segment (an edge of 
!  the boundary) is splitted at the mid of the edge, forming one new vertex 
!  and one new element
!  The subroutine repeatedly continues until no more vertex encroaches upon a segment
!
!  Input:    encrlimit  angle at which a segment is considered to be encroached 
!            meshsiz2   square of mesh-size for each of the regions
!  Output:
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
      integer (I4B) i, j, swaps, nn, isucc(3), ipred(3)
      integer (I4B), allocatable :: slist(:)
      logical encroached, newenc, splitted
      real (DP) x3, y3, msize
      parameter (isucc=(/2,3,1/))  !  successor
      parameter (ipred=(/3,1,2/))  !  predecessor
!
      newenc=.true.
      do while (newenc)
        newenc=.false.
!  de-encroach segments
        nn=n
        allocate(slist(3*nn))
        slist=0
        do i=1,nn
          do j=1,3
!  if the element has a neighbor and it is in the same region it cannot be a segment
            if (en(j,i) .ne. 0) then
              if (geb(i) .eq. geb(en(j,i)) ) cycle
            end if
!  if it is an input segment and size is to large : split the segment
            msize=( xn(e(isucc(j),i)) - xn(e(ipred(j),i)) )**2 +        &
     &            ( yn(e(isucc(j),i)) - yn(e(ipred(j),i)) )**2
            if (msize .gt. meshsiz2(geb(i))) then
              call splitencrelement(i,j,slist,splitted)
              newenc=newenc.or.splitted
              exit
            end if
            x3=xn(e(j,i))
            y3=yn(e(j,i))
!  prevent from successive splitting of small input angles
!          do not split if vertex (x3,y3) is on input branch
            if (kzi(e(j,i)) .ne. 0 ) cycle
!  is this edge an input segment and does it encroach on (x3,y3)
            call encrtest(i,j,x3,y3,encrlimit,e,en,geb,kzi,xn,yn,encroached)
            if (encroached) then
!  if so split the segment
              call splitencrelement(i,j,slist,splitted)
              newenc=newenc.or.splitted
              exit
            end if
          end do
        end do
        if (newenc) then
!  if a segment was splitted restore the delaunay triangulation
          call lswap(n,e,xn,yn,en,geb,0,.false.,swaps,5,slist)
        end if
        deallocate(slist)
!  loop until no more point encroaches upon a segment
      end do
      return
      end subroutine decroach
