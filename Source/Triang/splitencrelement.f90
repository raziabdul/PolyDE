      subroutine splitencrelement(i,j,slist,splitted)
      use feminterface, only: neukn, reallocate
      use globalvariables
      use femtypes
      implicit none
      integer (I4B) i, j, slist(:)
      logical splitted
      intent (in) :: i ,j
      intent (inout) :: slist
      intent (out) splitted
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
!    $Date: 2008/12/22 15:00:10 $
!    $Author: m_kasper $
!
!  split an element i into two by subdividing edge j
!            / \                     /|\
!           /   \         =>        / | \
!          /     \                 /  |  \
!         /_ _j_ _\               /_ _|_ _\
!  Input:
!            i        element to be splitted
!            j        edge to be splitted
!            zki      branch information, geometry nodes which determine the branches
!            xbk      coordinates of the input or geometry node (key-points)
!            ybk
!  In-/ Output:
!            n        total number of elements
!            p        total number of nodes
!            e        element information, nodes of the elements
!            en       neighborhood information, neighbors of the elements
!            geb      assignment of elements to the geometry regions
!            kzi      node branch information
!            xn, yn   coordinates of the nodes
!            x        solution vector
!            slist    a list of triangles to be used in the subsequent swaping process
!                     slist(i) is set to 1 if the triangle had been modified or created
!  Output;   splitted =.false. if splitting was not feasible
!
!  local variables
      integer (I4B) nbr, geb2, k, j4, n4, neuk, err, j1, j2, n1, n2, n3
      logical arc
      integer (I4B) isucc(3), ipred(3), jt
      parameter (isucc=(/2,3,1/))  !  successor
      parameter (ipred=(/3,1,2/))  !  predecessor

      if (n+2 .gt. size(geb)) then
        geb=>reallocate(geb,2*n)
        e=>reallocate(e,3,2*n)
        en=>reallocate(en,3,2*n)
      end if
      if (p+1 .gt. size(xn)) then
        xn=>reallocate(xn,2*p)
        yn=>reallocate(yn,2*p)
        kzi=>reallocate(kzi,2*p)
        x=>reallocate(x,2*p)
      end if

      jt=j
      j1=isucc(jt)
      j2=ipred(jt)
      n1=e(j1,i)
      n2=e(j2,i)
      n3=e(jt,i)
      nbr=en(jt,i)
      if (nbr .ne. 0) then
        geb2=geb(nbr)
      else
        geb2=0
      end if
!  determine the node in the neighbor element
      if (nbr.ne.0) then
        do k=1,3
          if (en(k,nbr) .eq. i) then
            j4=k
            n4=e(k,nbr)
            exit
          end if
        end do
      else
        n4=0
      end if
!  generate new node 
      call neukn(n1,n2,n3,n4,neuk,xn,yn,x,kzi,zki,p,xbk,ybk,            &
     &  geb(i),geb2,arc,err)
      if (err .gt. 0) then
!  Splitting of this element-edge is not possible
!  This may happen if the edge lies on an arc and is relatively large
!  compared to the radius
!  It may be useful to split at one of the other edges (kk loop)
        print*,'** an error occured while positioning of a new node'
        splitted=.false.
        return
      end if
!
!  generate new elements
      n=n+1
      slist(i)=1
      slist(n)=1
      e(1,n)=neuk
      e(2,n)=n2
      e(3,n)=n3
      en(1,n)=en(j1,i)
      en(2,n)=i
      if (nbr.ne.0) then
        en(3,n)=n+1
      else
        en(3,n)=0
      end if
      geb(n)=geb(i)
!  modify existing elements
      e(j2,i)=neuk
      en(j1,i)=n
      if (en(1,n).ne.0) then
        do k=1,3
          if (en(k,en(1,n)).eq.i) then
            en(k,en(1,n))=n
            exit
          end if
        end do
      end if
!  generate at the adjacent element
      if (nbr.ne.0) then
        n=n+1
        slist(n)=1
        slist(nbr)=1
        e(1,n)=neuk
        e(2,n)=n4
        e(3,n)=n2
        en(1,n)=en(ipred(j4),nbr)
        en(2,n)=n-1
        en(3,n)=nbr
        geb(n)=geb2
!  modify existing element
        e(isucc(j4),nbr)=neuk
        en(ipred(j4),nbr)=n
        if (en(1,n).ne.0) then
          do k=1,3
            if (en(k,en(1,n)).eq.nbr) then
              en(k,en(1,n))=n
              exit
            end if
          end do
        end if
      end if
      splitted=.true.
      return
      end subroutine splitencrelement
