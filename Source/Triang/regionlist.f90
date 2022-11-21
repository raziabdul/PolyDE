      subroutine regionlist(bzi,bzip,bzil,zpz,branchpt,brnodes,region,  &
     &  head,nbds,numnodes)
      use feminterface, only:
      use globalvariables
      use femtypes
      use triangtypes
      implicit none
      integer (I4B) bzi(:), bzip(:), bzil(:), zpz(:)
      integer (I4B) branchpt(:), brnodes(:), region, nbds
      integer (I4B), pointer :: numnodes(:)
      type (arrpt), pointer :: head(:)
      intent (in) :: bzi, bzip, bzil, zpz, branchpt, brnodes, region
      intent (out) :: nbds
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
!    $Revision: 1.3 $
!    $Date: 2006/07/08 00:18:59 $
!    $Author: r_abdul $
!
!  generate a list of the boundary nodes of the region in the form of a double linked list
!  and (in the case of a multiply connected region) its inner region nodes, i.e. the holes
!
!  Input:
!            bzi      list of branches, including the inner branches
!                     in the sequence of regions
!            bzip     Pointer (index) to the start branch of a region in the  
!                     list bzi (for compact storage)
!            bzil     Pointer (index) onto bzip including all regions 
!                     which are inside a region (multiply connected regions)
!            branchpt index of the first node of branches in brnodes 
!            brnodes  indices the nodes on the branches 
!                     (start- and endpoint of the branch are not stored in this list)
!            region   index of the actual region
!  Output:
!            head     for each boundary a pointer to the starting node
!            nbds     number of boundaries (e.g. =1 for a simply connected region)
!            numnodes number of nodes of the boundaries
!
!  local variables
      integer (I4B) i, j, k, bdy, branch, starti, sumnodes
      type (node), pointer :: nodei, firstnode, lastnode
! 
!  number of boundaries (e.g. =1 for a simply connected region)
      nbds=bzil(region+1)-bzil(region)
!  allocate a vector of pointers 
      allocate (head(nbds),numnodes(nbds))
!
!  for all the boundaries (outer and inner) of this region
      do k=1,nbds
        bdy=bzil(region)+k-1
        numnodes(k)=0
        sumnodes=sum(numnodes(1:k))
        do i=1,bzip(bdy+1)-bzip(bdy)
          branch=bzi(bzip(bdy)+i-1)
!  add the nodes of the branch;    i)   the starting node
!  starting node only needs to be added for the first branch
          if (i.eq.1) then
            allocate(nodei)
            numnodes(k)=numnodes(k)+1
            nodei%name=sumnodes+numnodes(k)
            lastnode=>nodei
            firstnode=>nodei
            head(k)%p=>nodei
            if (branch.gt.0) then
              nodei%nodeindx=zki(1,branch)
            else
              nodei%nodeindx=zki(2,-branch)
            end if
          end if
!  add the nodes of the branch;    ii)  the nodes on the branch
          if (branch.gt.0) then
!  form start- to end-point
            starti=branchpt(branch)
            do j=1,zpz(branch)-2
              allocate(nodei)
              nodei%prev=>lastnode
              nodei%name=sumnodes+numnodes(k)+j
              lastnode%next=>nodei
              lastnode=>nodei
              nodei%nodeindx=brnodes(starti+j-1)
            end do
            numnodes(k)=numnodes(k)+zpz(branch)-2
          else
!  form end- to start-point
            starti=branchpt(-branch+1)-1
            do j=1,zpz(-branch)-2
              allocate(nodei)
              nodei%prev=>lastnode
              nodei%name=sumnodes+numnodes(k)+j
              lastnode%next=>nodei
              lastnode=>nodei
              nodei%nodeindx=brnodes(starti-j+1)
            end do
            numnodes(k)=numnodes(k)+zpz(-branch)-2
          end if
!  add the nodes of the branch;    iii) the endpoint of the branch
          if (i.ne.bzip(bdy+1)-bzip(bdy)) then
            allocate(nodei)
            numnodes(k)=numnodes(k)+1
            nodei%prev=>lastnode
            nodei%name=sumnodes+numnodes(k)
            lastnode%next=>nodei
            lastnode=>nodei
            if (branch.gt.0) then
              nodei%nodeindx=zki(2,branch)
            else
              nodei%nodeindx=zki(1,-branch)
            end if
          else
!  this was the last branch, close the loop, but dont add this node 
            nodei%next=>firstnode
            firstnode%prev=>nodei
          end if
        end do
      end do
      return
      end subroutine regionlist
