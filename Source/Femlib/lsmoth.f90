      subroutine lsmoth(kzi,p,xn,yn,e,n,en,itermin,nregen,liste5,       &
     &  fehlkn,zrb,zki,xbk,ybk)
      use feminterface, only: getsetting, nachkn, pktply, zeit, kern3
      use feminterface, only: centerofgravity, reallocate, inkrp
      use femtypes
      implicit none
      integer (I4B) :: e(:,:), kzi(:), n, p, en(:,:), itermin
      integer (I4B) :: liste5(:), fehlkn, zrb(:,:), zki(:,:)
      real (DP) :: xn(:), yn(:), xbk(:), ybk(:)
      logical :: nregen
      intent (in) :: kzi, p, e, n, en, itermin, zrb, zki, xbk, ybk
      intent (out) :: nregen
      intent (inout) :: xn, yn, liste5, fehlkn
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
!    $Date: 2008/12/22 12:46:06 $
!    $Author: m_kasper $
!
!  Subroutine to smooth the triangle mesh by node relocation
!  for Laplace smoothing the nodes are shifted to the centre of gravity of the neighboring elements
!  for kernel smoothing the nodes are shifted to the centre of gravity of the kernel of neighboring elements
!
!  Input:
!            kzi      node branch information
!                     kzi=0: the node is an inner node
!                     kzi<0: the node is identical to the keypoint with number (-kzi)
!                     kzi>0: the node is on the branch with the number (kzi)
!            p        total number of nodes
!            e        element information, nodes of the elements
!            n        total number of elements
!            en       neighborhood information, neighbors of the elements
!            iterat   number of smooth iterations
!            zrb      Type of the boundary condition of the branch, 
!                     which can take the following values:
!                       0 -  99: Dirichlet BC, where for 1- 99: the BC has to be set in the file USERBC
!                     200 - 299: Neumann of general BC
!                     300 - 399: visible contour, no BC
!                     400 - 499: non-visible contour, no BC
!            zki      list of of key-points (geometry nodes) which determine the branches
!                     zki(1,i)  start key-point of the branch
!                     zki(2,i)  end key-point of the branch
!                     zki(3,i)  third key-koint of the branch, with
!                               /  = 0  :  straight line 
!                     zki(3,i)  -  > 0  :  node located on an arc
!                               \  < 0  :  midpoint of the arc
!            xbk, ybk coordinates of the input or geometry node (key-points)
!
!  Output:
!            nregen   =.true. if the mesh has to be regenerated,
!                     i.e. if an element has negative area
!
!  In-/ Output: 
!            xn, yn   coordinates of the nodes
!            liste5   =1 , if the node has to be relocated
!            fehlkn   number of nodes which have to be relocated
!
!  locale variables
      integer (I4B) :: pmax, kanz(p+1), kknot(6*p)
      integer (I4B) :: j, k, pz, iter, itermax, pzahl, kz1
      character (len=30) :: smoothing
      real (DP), pointer :: xk(:),yk(:)
      real (DP) :: rho, xm, ym, rhocp, rho1, alpha, lo, hi
      real (DP) :: xneu, yneu, xi, yi, flae2, ri, xctr, yctr
      logical :: neg, ready
      parameter (itermax=10)
!
      if (p.le.1) goto 998
!
!  get the neighbors of the nodes and whether they are relocated
!  kknot    contains the list neihgbors (compact storage)
!  kanz     contains the number of neighbors for each node
      call nachkn(p,n,kanz,kknot,e,en,kzi,zrb,zki,liste5)
!
!  pmax =maximum number of neighboring nodes to one node
      pmax=25
      allocate(xk(pmax),yk(pmax))
!
      ready=.false.

!  smoothing "Laplace"  = Laplace smoothing
!            "kernel"   = center of gravity of kernel
      call getsetting('SMOOTHING',smoothing)
      iter=0
      do while (.not. ready)
        do j=1,p
          pzahl=kanz(j+1)-kanz(j)
          if (pzahl.le.0) cycle
          if (pzahl**2 .gt. pmax) then
!  increase pmax
            pmax=max(pzahl**2,2*size(xk))
            xk=>reallocate(xk,pmax)
            yk=>reallocate(yk,pmax)
          end if
          kz1=kzi(j)
!
          if (kz1.eq.0) then
            select case(smoothing)
            case('LAPLACE')
!  compute the midpoint (in fact not the center of gravity)
!  i.e. Laplacian smoothing
              xm=sum(xn(kknot(kanz(j):kanz(j+1)-1)))/real(pzahl,DP)
              ym=sum(yn(kknot(kanz(j):kanz(j+1)-1)))/real(pzahl,DP)
            case('KERNEL')         
!  compute the kernel of the neihgboring nodes and it's center 
              call kern3(kknot,kanz(j),pzahl,xn,yn,xk,yk,pz)
              if (pz .gt. 0) then 
                call centerofgravity(xk,yk,pz,xm,ym)
              else
                print*,'empty polygon'
                stop
              end if
!  we have a new valid positioning of the node, correct the list of malpositined nodes
              if (liste5(j) .ne. 0) then
                liste5(j)=0
                fehlkn=fehlkn-1
                write (*,322) j
              end if
!  this should always be a valid position
!  new position of the node
              xn(j)=xm
              yn(j)=ym
              cycle
            case default
              print*,'wrong parameter for SMOOTHING = ',smoothing,      &
     &          ' in setting'
              stop
            end select
!
          else if (kz1.gt.0 .and. liste5(j).ne.0) then
!  if a node lies on branch it can be relocate on the branch
!  by projection of the node to the branch
!
            if (zki(3,kz1).ne.0) then
!  a node on an arc
              if (zki(3,kz1).lt.0) then
!  arc given by startpoint, endpoint and center
                xctr=xbk(-zki(3,kz1))
                yctr=ybk(-zki(3,kz1))
                rho=sqrt((xbk(zki(1,kz1))-xctr)**2+                     &
     &              (ybk(zki(1,kz1))-yctr)**2)
                rhocp=sqrt((xbk(zki(2,kz1))-xctr)**2+                   &
     &              (ybk(zki(2,kz1))-yctr)**2)
                rho=(rho+rhocp)/2._DP
!  rho  = desired radius
!  rho1 = actual radius
              else if (zki(3,kz1) .gt. 0) then
!  arc given by startpoint, endpoint and third point
                call inkrp(xbk(zki(1,kz1)),ybk(zki(1,kz1)),             &
     &            xbk(zki(2,kz1)),ybk(zki(2,kz1)),                      &
     &            xbk(zki(3,kz1)),ybk(zki(3,kz1)),                      &
     &            xi,yi,flae2,ri,xctr,yctr,rho)
              end if
              rho1=sqrt((xn(j)-xctr)**2+(yn(j)-yctr)**2)
              xm=xctr+rho/rho1*(xn(j)-xctr)
              ym=yctr+rho/rho1*(yn(j)-yctr)
            else
!  a node on a staight line or on an arc which will not be relocated
              cycle
            end if
          else
!  nodes on a straight line: do nothing
            cycle
          end if
!  check whether we're done
!  i.e. the node is located inside the polygon of the neighbors
          call pktply(j,pzahl,kanz,kknot,xm,ym,xn,yn,neg)
          if (neg) then
!  we have to find a new valid point
            if (kzi(j) .ne. 0) then
!  for nodes on branches (arcs)
!  xm, ym is the goal push the node into this direction as far as possible
              hi=1._DP
              lo=0._DP
!  use nested dissection
              do k=1,5
                alpha=(lo+hi)/2._DP
                xneu=xm*alpha+xn(j)*(1._DP-alpha)
                yneu=ym*alpha+yn(j)*(1._DP-alpha)
                call pktply(j,pzahl,kanz,kknot,xneu,yneu,xn,yn,neg)
                if (neg) then 
                  hi=alpha
                else 
                  lo=alpha
                end if
              end do
              xm=xm*lo+xn(j)*(1._DP-lo)
              ym=ym*lo+yn(j)*(1._DP-lo)
            else
!  for inner nodes (not on a branch)
!  compute the kernel of the neihgboring nodes and it's center 
              call kern3(kknot,kanz(j),pzahl,xn,yn,xk,yk,pz)
              if (pz .gt. 0) then 
                call centerofgravity(xk,yk,pz,xm,ym)
              else
                print*,'empty polygon'
                cycle
              end if
            end if
          else
!  we have a new valid positioning of the node, correct the list of malpositined nodes
            if (liste5(j) .ne. 0) then
              liste5(j)=0
              fehlkn=fehlkn-1
              write (*,322) j
322           format('   Knoten',i7,' erfolgreich verschoben')
            end if
          end if
!  new position of the node
          xn(j)=xm
          yn(j)=ym
        end do
        iter=iter+1
!  continue until all nodes are correctly positioned
!  or the maximum number of iterations is exceeded
        if ((iter .ge. itermin .and. fehlkn .eq. 0)                     &
     &    .or. iter .ge. itermax) then
          ready=.true.
        end if
      end do
      deallocate(xk,yk)
!
998   if (fehlkn.le.0) nregen=.false.
      call zeit(' Smoothing the Mesh')
      return
      end subroutine lsmoth
!
!
!
      subroutine nachkn(p,n,kanz,kknot,e,en,kzi,zrb,zki,liste5)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: p, n, kanz(:), kknot(:), e(:,:), en(:,:), kzi(:)
      integer (I4B) :: zrb(:,:), zki(:,:), liste5(:)
      intent (in) :: p, n, e, en, kzi, zrb, zki, liste5
      intent (out) :: kanz, kknot
!
!  Determine neighboring nodes of nodes
!
!  Input:
!            p        total number of nodes
!            n        total number of elements
!            e        element information, nodes of the elements
!            en       neighborhood information, neighbors of the elements
!            kzi      Knoten-Zweig-Information; Zweige zu den Knoten
!            kzi      node branch information
!                     kzi=0: the node is an inner node
!                     kzi<0: the node is identical to the keypoint with number (-kzi)
!                     kzi>0: the node is on the branch with the number (kzi)
!            zrb      Type of the boundary condition of the branch, 
!                     which can take the following values:
!                       0 -  99: Dirichlet BC, where for 1- 99: the BC has to be set in the file USERBC
!                     200 - 299: Neumann of general BC
!                     300 - 399: visible contour, no BC
!                     400 - 499: non-visible contour, no BC
!            zki      list of of key-points (geometry nodes) which determine the branches
!                     zki(1,i)  start key-point of the branch
!                     zki(2,i)  end key-point of the branch
!                     zki(3,i)  third key-koint of the branch, with
!                               /  = 0  :  straight line 
!                     zki(3,i)  -  > 0  :  node located on an arc
!                               \  < 0  :  midpoint of the arc
!            liste5   =1 , if the node has to be relocated
!
!  Output:
!            kanz     Vector for compact storage, pointing to the elements in kknot at: kanz(j) .. kanz(j+1)-1
!            kknot    list of neighbors (nodes) for all nodes in compact storage
!                     neighbors are sorted in mathematically positive direction
!
      integer (I4B) :: i, j, jj, node, khanz, copy, pos, start, eneu
      integer (I4B) :: isucc(3), liste(p)
      logical :: done
      parameter (isucc=(/2,3,1/))  !  successor
!
      kanz(1:p+1)=0
!  loop over all elements and for each node determine the number of neighbors
      do i=1,n
!  each inner node has as many neighbors as adjacent elements
!  each boundary node has one neighbors less than adjacent elements
          kanz(e(1:3,i))=kanz(e(1:3,i))+1
      end do
      do i=1,p
        select case (kzi(i))
        case (0)
!  a inner node which will be relocated
          liste(i)=1
        case (1:)
!  nodes on branches
          if (any(zrb(kzi(i),:).lt.300)) then
!  node at a boundary segment
            if (liste5(i).ne.0) then
!  a node at a boundary segment which has to be relocated
              liste(i)=1
!  correct number of neighbors
              kanz(i)=kanz(i)+1
            else
!  a node at a boundary segment which will not be relocated
              kanz(i)=0
              liste(i)=0
            end if
          else
!  nodes on inner branches will be relocated
            liste(i)=1
          end if
        case (:-1)
!  key-point cannot be relocated
          kanz(i)=0
          liste(i)=0
        end select
      end do
!
!  build the sum of neighbors for compact storage
      khanz=kanz(1)
      kanz(1)=1
      do i=2,p+1
        copy=kanz(i)
        kanz(i)=khanz+kanz(i-1)
        khanz=copy
      end do
!
!  collect neighbors in vector kknot
!  sort in mathematically positive orientation
      do i=1,n
        do j=1,3
!  search for the neighbors
          node=e(j,i)
!  if the node had bee treated (liste=2) is should not be considered (liste=0)
          if (liste(node).ne.1) cycle
!  determine starting node for sorting
          start=e(isucc(j),i)
          eneu=en(isucc(j),i)
          kknot(kanz(node))=start
          pos=kanz(node)+1
!
          if (kzi(node).gt.0) then
!  we will not relocate at boundary segments
            if (any(zrb(kzi(node),:).lt.300)) cycle
!  if start is a inner node
            if (kzi(start).eq.0) cycle
            if (kzi(start).gt.0) then
              if (kzi(start).ne.kzi(node)) cycle
            else
              if ((zki(1,kzi(node)).ne.-kzi(start)).and.                &
     &          (zki(2,kzi(node)).ne.-kzi(start))) cycle
            end if
          end if
!  starting with element i we follow the neighboring elements around the node
          done=.false.
          do while(.not.done) 
            do jj=1,3
              if (e(jj,eneu).eq.node) then
                if (e(isucc(jj),eneu).eq.start) then
!  tour is close, we are back at the stating node
                  done=.true.
                  exit
                end if
!  store node in list of neighbors
                kknot(pos)=e(isucc(jj),eneu)
!  choose next element to procede with
                eneu=en(isucc(jj),eneu)
                pos=pos+1
                exit
              end if
            end do
          end do
!  mark the node as done
          liste(node)=2
          if (pos.ne.kanz(node+1)) then
            write(*,1) node
1           format('   ***** Error in lsmoth at node ',i6)
          end if
        end do
      end do
      return
      end subroutine nachkn
!
!
!
      subroutine pktply(j,pzahl,kanz,kknot,xm,ym,xn,yn,neg)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: j, pzahl, kanz(:), kknot(:)
      real (DP)  :: xn(:), yn(:), xm, ym
      logical :: neg
      intent (in) :: j, pzahl, kanz, kknot, xm, ym, xn, yn
      intent (out) :: neg
!
!  Test, whether the point (xm,ym) is visible from each of the vertices (and edges)
!  of the polygon, which is fulfilled if the point is located inside the kernel of the polygon
!
!  Input:
!            j        node index whose neighbors form the polygon
!            pzahl    number of nodes fo the polygon
!            kanz     Aufaddierte Summe der Anzahl der Nachbarn der Knoten
!                     die Anzahl der Nachbarknoten ist pzahl=kanz(j+1)-kanz(j)
!            kknot    Liste der Nachbarknoten fuer alle Knoten 
!            xm, ym   coordinates of the point to test
!            xn, yn   coordinates of the nodes
!
!  Output:
!            neg      .true. if visibility is not fulfilled
!
!  local variables
      integer (I4B) m2, i
      real (DP) dx1, dx2, dy1, dy2
!
      m2=kknot(kanz(j)+pzahl-1)
      dx2=xn(m2)-xm
      dy2=yn(m2)-ym
      do i=0,pzahl-1
!  For each of the triangles formed by (xm,ym) and two consecutive
!  vertices of the polygon, test whether the area is positive
        dx1=dx2
        dy1=dy2
        m2=kknot(kanz(j)+i)
        dx2=xn(m2)-xm
        dy2=yn(m2)-ym
        if ( dy2*dx1-dy1*dx2 .lt. 0._DP) then
          neg=.true.
          return
        end if
      end do
      neg=.false.
      return
      end subroutine pktply
