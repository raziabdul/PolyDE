      subroutine minmesh(bzi,bzip,bzil,zpz,zpp, meshsize)
      use feminterface, only: genkz, regionlist, intbridges
      use feminterface, only: earclipper
      use globalvariables
      use femtypes
      use triangtypes
      implicit none
      integer (I4B) bzi(:), bzip(:), bzil(:), zpz(:)
      real (DP) zpp(:), meshsize(:)
      intent (in) :: bzi, bzip, bzil, zpp, meshsize
      intent (inout) :: zpz
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
!    $Date: 2014/02/17 16:27:26 $
!    $Author: m_kasper $
!
!  Generate a triangle mesh with minimum number of triangles 
!  i.e. no inner nodes are generated
!
!  Input:
!            gbz      number of regions
!            gzz      number of branches
!            gkz      number of key-points
!            bzi      list of branches, including the inner branches
!                     in the sequence of regions
!            bzip     Pointer (index) to the start branch of a region in the  
!                     list bzi (for compact storage)
!            bzil     Pointer (index) onto bzip including all regions 
!                     which are inside a region (multiply connected regions)
!            zki      list of of key-points (geometry nodes) which determine the branches
!                     zki(1,i)  start key-point of the branch
!                     zki(2,i)  end key-point of the branch
!                     zki(3,i)  third key-point of the branch, with
!                               /  = 0  :  straight line 
!                     zki(3,i)  -  > 0  :  node located on an arc              
!                               \  < 0  :  midpoint of the arc
!            zrb      Type of the boundary condition of the branch, 
!                     which can take the following values:
!                       0 -  99: Dirichlet BC, where for 1- 99: the BC has to be set in the file USERBC 
!                     200 - 299: Neumann of general BC
!                     300 - 399: visible contour, no BC
!                     400 - 499: non-visible contour, no BC
!            zpz      number of nodes on the branch (including start-  and endpoint)
!            zpp      control parameter for the distribution of nodes
!                          /  = 0  :  equidistant distribution
!                     zpp  -  > 0  :  nodes are denser at the startpoint
!                          \  < 0  :  nodes are denser at the endpoint
!            xbk,ybk  coordinates of key-points
!
!  Output:
!            n        total number of elements
!            p        total number of nodes
!            e        element information, nodes of the elements
!            en       neighborhood information, neigbors of the elements
!            geb      assignment of elements to the geometry regions
!            kzi      node branch information
!                      kzi=0: inner node
!                      kzi<0: is identical to the keypoint with number (-kzi)
!                      kzi>0: the node is on the branch with the number (kzi)
!            xn, yn   coordinates of the nodes
!            xmin,xmax Extends of the problem,
!            ymin,ymax outermost bounds of the geometry in x- and y-direction
!
!  local variables
      integer (I4B) i, nbds
      integer (I4B), pointer :: branchpt(:), brnodes(:), numnodes(:)
      type (arrpt), pointer :: head(:)
!
!
!  generate all nodes on branches
      call genkz(zpz, zpp, branchpt, brnodes, bzi, bzil, bzip, meshsize)
!  generate a list of the boundary nodes of the region 
!  and (in the case of a multiply connected region) its inner region nodes, i.e. the holes
!
!  check and allocated arrays number of elements roughly twice the node number
      if (associated(geb)) then
        deallocate(geb)
      end if
      if (associated(e)) then
        deallocate(e)
      end if
      allocate(geb(2*p))
      allocate(e(3,2*p))
      n=0
!
      do i=1,gbz
        call regionlist(bzi,bzip,bzil,zpz,branchpt,brnodes,i,           &
     &  head,nbds,numnodes)

!  convert a multiply connected region into a simply conneted one 
!  by introducing bridges to the holes
        if (nbds .gt. 1) then
          call intbridges(head,nbds,numnodes)
        end if
!  triangulate a simply connected region by ear clipping
        call earclipper(head(nbds)%p,numnodes(nbds),i) 
        deallocate(numnodes, head)
      end do
      deallocate(branchpt, brnodes)
      call zeit(' generating minimal mesh')
      return
      end subroutine minmesh
