      subroutine intbridges(head,nbds,numnodes)
      use feminterface, only: qsortindex, schnit, flaech, wink
      use globalvariables
      use femtypes
      use triangtypes
      implicit none
      integer (I4B) nbds
      integer (I4B) :: numnodes(:)
      type (arrpt), pointer :: head(:)
      intent (in) :: nbds
      intent (inout) :: numnodes
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
!    $Date: 2010/08/27 12:32:57 $
!    $Author: m_kasper $
!
!  convert a multiply connected region into a simply conneted one 
!  by introducing bridges to the holes
!
!  Input:
!            head     for each boundary a pointer to the starting node
!            nbds     number of boundaries (e.g. =1 for a simply connected region)
!            numnodes number of nodes of the boundaries
!
!  local variables
      integer (I4B) :: i, j, k, sumnodes
      integer (I4B), allocatable :: indx(:)
      real (DP) :: xminbd(nbds), xe, ye, xi, yi, xk, yk, xl, yl
      real (DP) :: t1, t2, wmin, w1, w2, w3, w4, anglemin, anglebound
      logical found, ok, new
      type (node), pointer :: enode, inode, knode, lnode, saveprt
      type (node), pointer :: tmpptr, nextptr, newnode1, newnode2
!  some of the ideas in the code can be found in:
!  M. Held: FIST: Fast Industrial-Strength Triangulation of Polygons,
!           Algorithmica 30 (2001) 563--596
!
!  determine leftmost node of the boundaries and 
!  adjust head such that it always points to this node
      do k=1,nbds
        tmpptr=>head(k)%p
        nextptr=>tmpptr%next
        do i=1,numnodes(k)-1
          if (xn(nextptr%nodeindx) .lt. xn(tmpptr%nodeindx)) then
            tmpptr=>nextptr
          end if
          nextptr=>nextptr%next
        end do
        head(k)%p=>tmpptr
        xminbd(k)=xn(tmpptr%nodeindx)
      end do
      if (nbds .eq. 1) return
!  index sort the inner boundaries from left to right
      allocate(indx(nbds-1))
      call qsortindex(-xminbd(1:nbds-1),indx,nbds-1)
!  check whether leftmost node of the outer boundary is left to the 
!  leftmost point of the leftmost inner boundaries
      if (xminbd(nbds) .gt. xminbd(indx(1)) ) then 
        print*,'*** it may be advisable to use a larger number of nodes on arcs'
      end if
!  now introduce bridges
      sumnodes=sum(numnodes(1:nbds))
      do k=1,nbds-1
!  node on the inner boundary from which we look for a bridge
        inode=>head(indx(k))%p
        xi=xn(inode%nodeindx)
        yi=yn(inode%nodeindx)
!  cycle over all outer boundary nodes until a feasible bridge was found
        enode=>head(nbds)%p
        xe=xn(enode%nodeindx)
        ye=yn(enode%nodeindx)
        found=.false.
        wmin=1.e32
!  minimum acceptable angle for a bridge 
        anglebound=3._DP/nbds
        do i=1,numnodes(nbds)
!  a feasible candidate (enode) at the outer boundary must be left of the inner boundary
          if (xe .ge. xi) then
            enode=>enode%next
            xe=xn(enode%nodeindx)
            ye=yn(enode%nodeindx)
            cycle
          end if
!  a feasible candidate (enode) has to have (inode) inside its cone 
          lnode=>enode%prev
          xl=xn(lnode%nodeindx)
          yl=yn(lnode%nodeindx)
          if (flaech(xl,yl,xe,ye,xi,yi) .le. 0._DP) then   !  not needed, known from previous step ??
            enode=>enode%next
            xe=xn(enode%nodeindx)
            ye=yn(enode%nodeindx)
            cycle
          end if
          knode=>enode%next
          xk=xn(knode%nodeindx)
          yk=yn(knode%nodeindx)
          if (flaech(xe,ye,xk,yk,xi,yi) .le. 0._DP) then
            enode=>enode%next
            xe=xn(enode%nodeindx)
            ye=yn(enode%nodeindx)
            cycle
          end if
!  check whether the bridge (inode - enode) intersects with a segment (knode - lnode)
!  of the outer boundary
          lnode=>knode%next
          xl=xn(lnode%nodeindx)
          yl=yn(lnode%nodeindx)
          new=.true.
          do j=1,numnodes(nbds)-2
!  for segments not beeing adjacent to enode
            if (enode%nodeindx.ne.lnode%nodeindx .and.                  &
      &         enode%nodeindx.ne.knode%nodeindx) then
!  compute the intersection point between (enode - inode) and (knode - lnode)
              call schnit(xi,yi,xe,ye,xk,yk,xl,yl,t1,t2,ok)
              if ( (t1.le.1._DP) .and. (t2.le.1._DP) .and.              &
      &            (t1.ge.0._DP) .and. (t2.ge.0._DP) ) then
!  found an intersection; enode is not possible
                new=.false.
                exit
              end if
            end if
            knode=>lnode
            lnode=>knode%next
            xk=xn(knode%nodeindx)
            yk=yn(knode%nodeindx)
            xl=xn(lnode%nodeindx)
            yl=yn(lnode%nodeindx)
          end do
          if (new) then
!  no intersection found, (enode - inode) is possible bridge
            tmpptr=>enode%prev
            w1=wink(xe,ye  ,xi,yi,xn(tmpptr%nodeindx),yn(tmpptr%nodeindx))
            tmpptr=>inode%next
            w2=wink(xi,yi  ,xn(tmpptr%nodeindx),yn(tmpptr%nodeindx),xe,ye)
            tmpptr=>inode%prev
            w3=wink(xi,yi  ,xe,ye,xn(tmpptr%nodeindx),yn(tmpptr%nodeindx))
            tmpptr=>enode%next
            w4=wink(xe,ye  ,xn(tmpptr%nodeindx),yn(tmpptr%nodeindx),xi,yi)
            anglemin=min(w1,w2,w3,w4)
            if (anglemin .gt. anglebound) then
              found=.true.
              saveprt=>enode
              exit
            else
!  the angle is very small, continue and try to find something better
!  ... but use a smaller accepable ange
              anglebound=anglebound/2._DP
              if (found) then
                if (anglemin .gt. wmin) then
                  wmin=anglemin
                  saveprt=>enode
                end if
              else 
                found=.true.
                wmin=anglemin
                saveprt=>enode
              end if
            end if
          end if
          enode=>enode%next
          xe=xn(enode%nodeindx)
          ye=yn(enode%nodeindx)
        end do
        if (.not. found) then
          print*,'***** it was not possible to convert a multiply ',    &
     &            'connected region into a simple one'
          stop
        end if
! use the saved node of best angle
        enode=>saveprt
!  modify the outer polygon, introduce a bridge between enode and inode
        numnodes(nbds)=numnodes(nbds)+numnodes(k)+2
        sumnodes=sumnodes+2
        allocate(newnode1, newnode2)
!  new nodes at the location of  enode and inode
        newnode1%nodeindx=inode%nodeindx
        newnode1%next=>newnode2
        newnode1%prev=>inode%prev
        newnode1%name=sumnodes-1
        newnode2%nodeindx=enode%nodeindx
        newnode2%next=>enode%next
        newnode2%prev=>newnode1
        newnode2%name=sumnodes
!  modify existing links
        tmpptr=>enode%next
        tmpptr%prev=>newnode2
        tmpptr=>inode%prev
        tmpptr%next=>newnode1
        enode%next=>inode
        inode%prev=>enode
!  set drawing color to yellow= 0.00, 1.00, 0.00 (R,G,B)
        call pgsci(20)
        call pgscr(20,  0., 1., 0.)
!  draw bridges
        call pgmove(sngl(xi),sngl(yi))
        call pgdraw(sngl(xn(enode%nodeindx)),sngl(yn(enode%nodeindx)))
      end do
      deallocate(indx)
      return
      end subroutine intbridges
