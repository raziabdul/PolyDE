!
      subroutine lswap(n,e,xn,yn,en,geb,talk,all,szahl,ic,slist)
      use feminterface, only: lswap1
      use femtypes
      implicit none
      integer (I4B) :: n, e(:,:), en(:,:), geb(:), talk, szahl, ic
      integer (I4B), optional, target :: slist(:)
      real (DP) :: xn(:), yn(:)
      logical :: all
      intent (in) :: n, xn, yn, geb, ic, talk, all
      intent (out) :: szahl
      intent (inout) :: e, en, slist
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
!    $Date: 2014/02/17 13:51:37 $
!    $Author: m_kasper $
!
!  Lswap: Swaping of triangle edges for neighboring elements in order to improve mesh quality
!
!  Input: 
!            n        number of elements
!            xn,yn    node coordinates
!            geb      assignment of elements to the geometry regions
!            ic       control parameter for the criterion to choose
!                      1: Babuska, minimize the maximum angle
!                      2: Lawson, maximize the minimum angle
!                      3: shorter diagonal
!                      5: a variant of the circumcircle rule by Lawson (2)
!                          Maximize the minimum angle
!            talk     Statistical output at the screen
!                       2  output at each iteration
!                       1  report at the end of subroutine
!                       0  no output
!            all      edge flipping using a list
!                      .true. : always test all edges
!                      .false.: list oriented, only test flipped edges again
!  Output:
!            szahl    number of edge flips
!  In-/ Output: 
!            e        nodes of elements
!            en       neighbors of elements
!            slist    (optional) swap list with enties 
!                      = 1 if the element should be considered for swapping
!                      = 0 otherwise
!
!  local variables
      integer (I4B) :: i, j, zahl, elnum, ennum, swzahl,zaehl,nextel,maxlop
      integer (I4B), pointer :: liste(:)
      logical :: swappd
!  maximum number of cycles
      parameter (maxlop=200)
!
!  szahl: total number of edge flips, zaehl: counter for cycles
      szahl=0
      zaehl=0
!  set todo-list
!  in the first pass all elements will be tested if slist is not present
      if (present(slist)) then
        liste=>slist
      else
        allocate (liste(n))
        liste(1:n)=1
      end if
!
100   swzahl=0
      do i=1,n
        elnum=i
200     nextel=0
        if (liste(elnum).eq.1) then
!  only test those elements which had been altered during the last cycle
          liste(elnum)=0
!  unmark element for the next cycle, if not (all)
          if (all) liste(elnum)=1
!  for all edges of the element
          do j=1,3
            ennum=en(j,elnum)
!  if the neighbor is marked, we will test later
            if (ennum.eq.0) cycle
            if (ennum.gt.elnum .and. liste(ennum).eq.1) cycle
            if (elnum.eq.ennum) then
              write (*,4100)
4100          format(' *****Fatal Error in Swap:')
              write (*,4101) elnum
4101          format (' the element',i7,' is neighbor to itself')
              write (*,*) e(1,elnum),e(2,elnum),e(3,elnum),             &
     &          en(1,elnum),en(2,elnum),en(3,elnum)
              stop
            end if
!  we only may flip if the neighbor is in the same region
            if (geb(ennum) .eq. geb(elnum)) then
!  search local node number which dose not simultaneously belong
!  to elnum and ennum
              if (en(1,ennum).eq.elnum) then
                zahl=1
              else if (en(2,ennum).eq.elnum) then
                zahl=2
              else 
                zahl=3
              end if
              call lswap1(e,en,xn,yn,elnum,ennum,j,zahl,swappd,ic)
              if (swappd) then
!  swapped, put the elements into the todo-list
                liste(elnum)=1
                liste(ennum)=1
                swzahl=swzahl+1
                nextel=ennum
              end if
            end if
          end do
        end if
        if (nextel.ne.0) then
!  next element is one of the elements just swapped
          elnum=nextel
          goto 200
        end if
      end do
      if (talk.ge.2) write (*,2000) swzahl
2000  format('   ',i7,' swap operations')
!  number of flips in this cycle swzahl
      szahl=szahl+swzahl
      zaehl=zaehl+1
!  still elements in the todo-list?
      if ((swzahl.gt.0).and.(zaehl.le.maxlop)) goto 100
!
      if (talk .ge.1) then
        write (*,2001) szahl
2001    format('  ',i7,' swap operations')
      end if
!
      if (.not. present(slist)) then
        deallocate(liste)
      end if
      return
      end
!
!
!
      subroutine lswap1(e,en,xn,yn,elnum,ennum,z1,z2,swappd,ic)
      use feminterface, only: wink, winkl
      use femtypes
      implicit none
      integer (I4B) :: e(:,:), en(:,:), elnum, ennum, z1, z2, ic
      real (DP) :: xn(:), yn(:)
      logical :: swappd
      intent (in) :: xn, yn, elnum, ennum, z1, z2, ic
      intent (out) :: swappd
      intent (inout) :: e, en
!  Input:
!            xn,yn    node coordinates
!            elnum    first element to test for edge flipping
!            ennum    second element to test for edge flipping
!            z1       local edge of the first element for which to test
!            z2       local edge of the second element for which to test
!            ic       control parameter for the criterion to choose
!                      1: Babuska, minimize the maximum angle
!                      2: Lawson, maximize the minimum angle
!                      3: shorter diagonal
!                      5: a variant of the circumcircle rule by Lawson (2),
!                         maximize the minimum angle
!  Output:
!            swappd  =.true. if the edge was flipped
!  In-/ Output: 
!            e        nodes of elements
!            en       neighbors of elements
!
!  local variables
      real (DP) :: pi, twopi, delta
      parameter (pi=3.14159265358979323846264338327950288419716_DP)
      parameter (twopi=2*pi)
      integer (I4B) :: i1, i2, j1, j2, tmp1, tmp2, isucc(3), ipred(3)
      real (DP) :: w1, w2, w3, w4, w1new, w2new, wmax, wmnew, rl1, rl2
      real (DP) :: dx1, dx2, dx3, dx4, dy1, dy2, dy3, dy4, wlimit
      parameter (isucc=(/2,3,1/))  !  successor
      parameter (ipred=(/3,1,2/))  !  predecessor
!
      swappd=.false.
!  the element  elnum  has global node numbers:   i1,j1,j2
!  the element  ennum  has global node numbers:   i2,j2,j1
      i1=e(z1,elnum)
      i2=e(z2,ennum)
      j1=e(isucc(z1),elnum)
      j2=e(ipred(z1),elnum)
!
      select case(ic)
!  =============================================================
!  criterion of Babuska
!  minimize the maximum angle
      case (1)
        w1=wink(xn(i1),yn(i1),xn(j1),yn(j1),xn(j2),yn(j2))
        w2=wink(xn(i2),yn(i2),xn(j2),yn(j2),xn(j1),yn(j1))
        wmax=max(w1,w2)
        w1new=wink(xn(j1),yn(j1),xn(i2),yn(i2),xn(i1),yn(i1))
!       w2new=wink(xn(j2),yn(j2),xn(i1),yn(i1),xn(i2),yn(i2))
        w2new=twopi-(w1+w2+w1new)
        wmnew=max(w1new,w2new)
!    if (wmnew.lt.wmax-0.02d0) then
        if (wmnew.lt.wmax-0.02e-10_DP) then
          swappd=.true.
          if (w1new.le.0._DP .or. w2new.le.0._DP) then
            print*,' error in swap'
            swappd=.false.
          end if
        end if
!  =============================================================
!  criterion of Lawson
!  maximize the minimum angle
      case (2)
!  for the triangle with nodes 1, 2, 3 check whether the node 4 is inside the circumcircle
!  u=(x2-x4)*(y2-y1)-(x2-x1)*(y2-y4)
!  v=(x3-x4)*(y3-y1)-(x3-x1)*(y3-y4)
!  w=(x2-x4)*(x2-x1)+(y2-y4)*(y2-y1)
!  z=(x3-x4)*(x3-x1)+(y3-y4)*(y3-y1)
!  u*z > v*w ?
        dx2=xn(j1)-xn(i2)
        dx3=xn(j2)-xn(i2)
        dx1=xn(j1)-xn(i1)
        dx4=xn(j2)-xn(i1)
        dy2=yn(j1)-yn(i2)
        dy3=yn(j2)-yn(i2)
        dy1=yn(j1)-yn(i1)
        dy4=yn(j2)-yn(i1)
        if ((dx2*dy1-dx1*dy2)*(dx3*dx4+dy3*dy4) .gt.                    &
     &      (dx3*dy4-dx4*dy3)*(dx2*dx1+dy2*dy1)) then
          swappd=.true.
        end if
!  =============================================================
!  chose the shoter diagonal
      case (3)
        w3=winkl(i1,j1,i2,xn,yn)
!  larger or near to 180 degrees?
        if ((w3.gt.0._DP).and.((pi-abs(w3)).gt.0.05_DP)) then
          w4=winkl(i2,j2,i1,xn,yn)
          if ((w4.gt.0._DP).and.((pi-abs(w4)).gt.0.05_DP)) then
            rl1=(xn(i1)-xn(i2))**2+(yn(i1)-yn(i2))**2
            rl2=(xn(j1)-xn(j2))**2+(yn(j1)-yn(j2))**2
            if ((rl2-rl1).gt.0._DP) then
              swappd=.true.
            end if
          end if
        end if
!
!  =============================================================
      case (4)
!  =============================================================
!  a variant of the circumcircle rule by Lawson
!  Maximize the minimum angle
      case (5)
        w1=wink(xn(i1),yn(i1),xn(j1),yn(j1),xn(j2),yn(j2))
        w2=wink(xn(i2),yn(i2),xn(j2),yn(j2),xn(j1),yn(j1))
        delta=w1+w2-pi
!        wlimit=.01_DP ! 
        wlimit=100._DP*tiny(1._SP)
        if (delta .ge. wlimit) then
          swappd=.true.
        else if (abs(delta) .le. wlimit) then
!  if the four points lie on a circle both alternatives have the same minimum angle
!  instead we minimize the maximum angle in this case
          wmax=max(w1,w2)
          w1new=wink(xn(j1),yn(j1),xn(i2),yn(i2),xn(i1),yn(i1))
          w2new=twopi-(w1+w2+w1new)
          wmnew=max(w1new,w2new)
          if (wmnew.lt.wmax-tiny(1._SP)) then
            swappd=.true.
            if (w1new.le.0._DP .or. w2new.le.0._DP) then
              print*,' error in swap'
              swappd=.false.
            end if
          end if
        end if

      end select
!
      if (swappd) then
!  correct node and neighborhood information of elements
        e(isucc(z1),elnum)=i2
        e(isucc(z2),ennum)=i1
        tmp1=en(ipred(z1),elnum)
        tmp2=en(ipred(z2),ennum)
        en(z1,elnum)=tmp2
        en(z2,ennum)=tmp1
        en(ipred(z1),elnum)=ennum
        en(ipred(z2),ennum)=elnum
        if (tmp2.gt.0) then
          if (en(1,tmp2).eq.ennum) then
            en(1,tmp2)=elnum
          else if (en(2,tmp2).eq.ennum) then 
            en(2,tmp2)=elnum
          else if (en(3,tmp2).eq.ennum) then 
            en(3,tmp2)=elnum
          end if
        end if
        if (tmp1.gt.0) then
          if (en(1,tmp1).eq.elnum) then
            en(1,tmp1)=ennum
          else if (en(2,tmp1).eq.elnum) then
            en(2,tmp1)=ennum
          else if (en(3,tmp1).eq.elnum) then
            en(3,tmp1)=ennum
          end if
        end if
      end if
      return
      end
