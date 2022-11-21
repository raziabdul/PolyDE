      subroutine delaunay_refinement(meshsize)
      use feminterface, only: getsetting, decroach, angles, offcentre,  &
     &    elemnt, zeit, splitencrelement, lswap, flaech, decroachn,     &
     &    reallocate, netbew, triplt, glaett, schnit
      use femtypes
      use globalvariables
      implicit none
      real (DP), pointer:: meshsize(:)
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
!    $Revision: 1.13 $
!    $Date: 2015/06/02 09:22:17 $
!    $Author: m_kasper $
!
!  Delaunay refinement with guaranteed bound on the minimum angle
!  The algorithm is a simplified version of Rupperts algorithm with 
!  some modifications given by Shewchuk
!
!  Input:
!            zki      branch information, geometry nodes which determine the branches
!            xbk      coordinates of the input or geometry node (key-points)
!            ybk
!  In-/Output:
!            n        total number of elements
!            p        total number of nodes
!            e        element information, nodes of the elements
!            en       neighborhood information, neighbors of the elements
!            geb      assignment of elements to the geometry regions
!            kzi      node branch information
!            xn, yn   coordinates of the nodes
!
!  local variables
      integer (I4B) :: i, j, k, nn, ictr, swaps, pass, minang(1)
      integer (I4B) :: isucc(3), ipred(3), maxang(1), br1, br2, stage
      integer (I4B) :: nskinny, fehlkn, numsize
      integer (I4B), allocatable :: slist(:), liste5(:)
      real (DP) ang(3), skinnylimit, skl, encrlimit, betamax, xctr, yctr
      real (DP) summ, f1, f2, f3, del, msize, maxsize, t1,t2
      real (DP), allocatable :: meshsiz2(:)
      logical ok, encroached, newskinny, splitted, nregen, segment
      parameter (isucc=(/2,3,1/))  !  successor
      parameter (ipred=(/3,1,2/))  !  predecessor
!
!  The algorithm starts with an initial triangulation of the domain, applying a
!  successive refinement strategy and using two consecutive steps:
!  1)  The splitting of encroached segments. 
!      A vertex is call encroached upon a segment if the angle under which the segment 
!      is seen is larger than a certain limit (encrlimit). Such a segment (an edge of 
!      the boundary) is splitted at the mid of the edge
!  2)  Subdivision of skinny triangles.
!      If a skinny triangle (having a smallest angle < skinnylimit) is present, it is 
!      removed by insterting a new vertex at its circumcenter and splitting the triangle
!      (in which the circumcenter is located) into three sub-triangles
!  After each modification of the mesh by one of the refinement rules the Delaunay 
!  triangulation is restored by edge swaping with the Lawson criterion
!  References:
!  Ruppert, Jim: "A Delaunay Refinement Algorithm for Quality 2-Dimensional
!    Mesh Generation", Journal of Algorithms, 18, p 548-585, (1994)
!  Shewchuk, Jonathan Richard: "Delaunay Refinement Algorithms for Triangular
!    Mesh Generation, Computational Geometry: Theory and Applcations, 22, p21-74, (2002)
!
!  skinnylimit: bound for the smallest angle in the triangulation
!
!  The algorithms usually works fine for skinnylimit < 30 deg, 
!  for a larger angle than 25 deg the algorithm may not terminate
      call getsetting('SKINNYLIMIT',skl)
!  get maximum mesh size
      call getsetting('MAXMESHSIZE',maxsize)
      allocate (meshsiz2(size(meshsize)))
      do i=1,gbz
!  assign the square of maximum-mesh-size for each of the regions
        meshsiz2(i)=min(meshsize(i),maxsize)**2
      end do
      skinnylimit=skl * pi/180._DP
!  encrlimit:  angle at which a point is considered to encroach upon a segment
      encrlimit=(1._DP-1.e-4_DP)*min(120._DP,180._DP-2._DP*skl) * pi/180._DP
!  the maximum allowed radius edge raion of triangles
      betamax=(1._DP-1.e-3_DP)/(2._DP*sin(skinnylimit))
!
      allocate (x(size(xn)))
      x(1:p)=0._DP
      call lswap(n,e,xn,yn,en,geb,0,.true.,swaps,2)
!
!  wipe the window
      call pgeras
!  draw triangulation
      call triplt(e,en,xn,yn,n,0.80,0.80,0.80)
      newskinny=.true.
      pass=0
!  start by only considering the mesh size criterion
      stage=0
!  check whether a vertex is encroached upon a segment and split the segments
      call decroach(encrlimit,meshsiz2)
!  eleminate skinny triangles
      do while (newskinny .or. stage.eq.0)
        pass=pass+1
        newskinny=.false.
        nskinny=0
        nn=n
        numsize=0
!
        do i=1,nn
!  test for skinny triangles
          call angles(i,e,xn,yn,ang(1),ang(2),ang(3))
!  the largest edge is opposite to the largest angle 
          maxang=maxloc(ang)
          msize=( xn(e(isucc(maxang(1)),i)) - xn(e(ipred(maxang(1)),i)) )**2 + &
     &          ( yn(e(isucc(maxang(1)),i)) - yn(e(ipred(maxang(1)),i)) )**2
          minang=minloc(ang)
!  in stage 0 only consider meshsize
          if ((ang(minang(1)) .gt. skinnylimit .or. stage .eq. 0)        &
     &      .and. msize .le. meshsiz2(geb(i))) cycle
          if (msize .gt. meshsiz2(geb(i))) numsize=numsize+1
!  test whether the node with minimum angle is a key-point
!  and the other two element nodes are located on distinct branches 
!  having the key-point as a common endpoint
!  this is indended to avoid subdivision of small input angles
          if (ang(minang(1)) .lt. skinnylimit) then
            nskinny=nskinny+1
            if (kzi(e(minang(1),i)) .lt. 0 ) then
              br1=kzi(e(isucc(minang(1)),i))
              br2=kzi(e(ipred(minang(1)),i))
              if (br1.gt.0 .and. br2.gt.0 .and. br1.ne.br2) then
                if (zki(1,br1).eq.zki(1,br2) .or.                       &
     &              zki(1,br1).eq.zki(2,br2) .or.                       &
     &              zki(2,br1).eq.zki(1,br2) .or.                       &
     &              zki(2,br1).eq.zki(2,br2)) then
                  cycle
                end if
              end if
            end if
          end if
          if (p+1 .gt. size(xn)) then
            xn=>reallocate(xn,2*p)
            yn=>reallocate(yn,2*p)
            kzi=>reallocate(kzi,2*p)
            x=>reallocate(x,2*p)
          end if
          if (n+3 .gt. size(geb)) then
            geb=>reallocate(geb,2*n)
            e=>reallocate(e,3,2*n)
            en=>reallocate(en,3,2*n)
          end if
!  if the triangle is skinny or minangle is below 50 deg we us the offcenter
          if (ang(minang(1)) .le. pi/3.6_DP ) then
!  the triangel is skinny, compute the circumcenter/offcenter of this triangle
            call offcentre(xn(e(1:3,i)),yn(e(1:3,i)),ang,minang(1),betamax,xctr,yctr)
!  determine the element in which the circumcenter/offcenter is located
            call elemnt(xn,yn,e,n,xctr,yctr,ictr,en,ok)
          else
!  if the minimum angle is large (meshsize requires refinement), we use the mid-point
            xctr=sum(xn(e(1:3,i)))/3._DP
            yctr=sum(yn(e(1:3,i)))/3._DP
            ictr=i
          end if
!  if the element is located adjacent to an input segment 
!  and the line from the host element to the circumcenter
!  intersects the segment => split the segment
!  use the segment opposite to the largest angle
          j=maxang(1)
          if (en(j,i) .eq. 0) then
            segment=.true.
          else
            if (geb(i) .ne. geb(en(j,i)) ) then
              segment=.true.
            else 
              segment=.false.
            end if
          end if
!  test for intersection
          if (segment) then
            call schnit(xn(e(j,i)),yn(e(j,i)),xctr,yctr,                &
     &                  xn(e(isucc(j),i)),yn(e(isucc(j),i)),            &
     &                  xn(e(ipred(j),i)),yn(e(ipred(j),i)),t1,t2,ok)
            if ( (t1.le.1._DP) .and. (t2.le.1._DP) .and.                &
     &           (t1.ge.0._DP) .and. (t2.ge.0._DP) ) then
              allocate(slist(n+2))
              slist(1:n+2)=0
              call splitencrelement(i,j,slist,splitted)
              call lswap(n,e,xn,yn,en,geb,0,.false.,swaps,2,slist)
              deallocate(slist)
              if (splitted) then
                newskinny=.true.
                cycle
              end if
            end if
          end if
!  check whether the circumcenter encroached on a segment of the triangle ictr
          call decroachn(encrlimit,xctr,yctr,encroached)
          if (encroached) cycle
!  check whether the circumcenter is outside the triangulated area or in another region
          if (.not. ok) cycle
          if (geb(i) .ne. geb(ictr)) cycle
!
!  check whether the new point is lying (almost) on one of the edges
!  if the point is (almost) identical with an already existing point 
!  replace it by the center of the element
          newskinny=.true.
          f1=flaech(xn(e(2,ictr)),yn(e(2,ictr)),xn(e(3,ictr)),yn(e(3,ictr)),xctr, yctr)
          f2=flaech(xn(e(3,ictr)),yn(e(3,ictr)),xn(e(1,ictr)),yn(e(1,ictr)),xctr, yctr)
          f3=flaech(xn(e(1,ictr)),yn(e(1,ictr)),xn(e(2,ictr)),yn(e(2,ictr)),xctr, yctr)
          del=1000._DP*spacing(f1+f2+f3)
          if (f1 .lt. del ) then
            if (f2.lt.del .or. f3.lt.del) then
!  this should not happen
              p=p+1
!  split host element instead
              xn(p)=sum(xn(e(1:3,i)))/3._DP
              yn(p)=sum(yn(e(1:3,i)))/3._DP
              ictr=i
            else
!  split at edge 1
              allocate(slist(n+2))
              slist(1:n+2)=0
              call splitencrelement(ictr,1,slist,splitted)
              call lswap(n,e,xn,yn,en,geb,0,.false.,swaps,2,slist)
              deallocate(slist)
              cycle
            end if
          else if (f2 .lt. del ) then
            if (f3.lt.del) then
!  this should not happen
              p=p+1
!  split host element instead
              xn(p)=(xn(e(1,i))+xn(e(2,i))+xn(e(3,i)))/3._DP
              yn(p)=(yn(e(1,i))+yn(e(2,i))+yn(e(3,i)))/3._DP
              ictr=i
            else
!  split at edge 2
              allocate(slist(n+2))
              slist(1:n+2)=0
              call splitencrelement(ictr,2,slist,splitted)
              call lswap(n,e,xn,yn,en,geb,0,.false.,swaps,2,slist)
              deallocate(slist)
              cycle
            end if
          else if (f3 .lt. del ) then
!  split at edge 3
            allocate(slist(n+2))
            slist(1:n+2)=0
            call splitencrelement(ictr,3,slist,splitted)
            call lswap(n,e,xn,yn,en,geb,0,.false.,swaps,2,slist)
            deallocate(slist)
            cycle
          else
!  everything is fine, use the circumcenter 
            p=p+1
            xn(p)=xctr
            yn(p)=yctr
          end if
!  subdivide the element ictr into three elements by introducing a vertex at xctr, yctr
!  and edges between this vertex and the already existing vertices of the triangle
!  generate two new elements and modify the existing element
          kzi(p)=0
!  first element
          e(1,n+1)=e(2,ictr)
          e(2,n+1)=e(3,ictr)
          e(3,n+1)=p
          en(1,n+1)=n+2
          en(2,n+1)=ictr
          en(3,n+1)=en(1,ictr)
!  correct the neigborhood list, identify the edge in the neighbor-element
          if (en(1,ictr).ne.0) then
            do k=1,3
              if (en(k,en(1,ictr)).eq.ictr) then
                en(k,en(1,ictr))=n+1
                exit
              end if
            end do
          end if
!  second element
          geb(n+1)=geb(ictr)
          e(1,n+2)=e(3,ictr)
          e(2,n+2)=e(1,ictr)
          e(3,n+2)=p
          en(1,n+2)=ictr
          en(2,n+2)=n+1
          en(3,n+2)=en(2,ictr)
!  correct the neigborhood list, identify the edge in the neighbor-element
          if (en(2,ictr).ne.0) then
            do k=1,3
              if (en(k,en(2,ictr)).eq.ictr) then
                en(k,en(2,ictr))=n+2
                exit
              end if
            end do
          end if
          geb(n+2)=geb(ictr)
!  modifiy the element
          e(3,ictr)=p
          en(1,ictr)=n+1
          en(2,ictr)=n+2
          n=n+2
!  restore the delaunay triangulation
          allocate(slist(n))
          slist(1:n-2)=0
          slist(ictr)=1
          slist(n)=1
          slist(n-1)=1
          call lswap(n,e,xn,yn,en,geb,0,.false.,swaps,2,slist)
          deallocate(slist)
        end do
        if (stage.le.0) then
!  do mesh smoothing if meshsize is not fulfilled
          allocate (liste5(p))
          nregen=.false.
          fehlkn=0
          liste5(:)=0
          call glaett(n,e,xn,yn,en,geb,kzi,p,nregen,liste5,             &
     &      fehlkn,zrb,zki,xbk,ybk,2)
          deallocate(liste5)
        else
          call netbew(xn,yn,e,summ,n)
        end if
        if ( (stage.eq.0) .and.                                         &
      &  ((numsize .le. 0.05*nn) .or. (newskinny.eqv..false.)) ) then
!  now switch to considere skinny triangles
          stage=1
          newskinny=.true.
        end if
        write(*,100) pass,n
100     format(' pass:',i4,'    elements: ',i7)
        if (stage .eq. 1) then
!  wipe the window
          call pgeras
!  draw triangulation
          call triplt(e,en,xn,yn,n,0.80,0.80,0.80)
        end if
      end do
      print*,"remaining skinny triangles", nskinny
      return
      end subroutine delaunay_refinement


      subroutine offcentre(x,y,alpha,minloc,betamax,xpa,ypa)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) x(:), y(:), alpha(:), betamax, xpa, ypa
      integer(I4B) :: minloc
      intent (in) :: x, y, alpha, betamax, minloc
      intent (out) :: xpa, ypa
!  computation of the circumcircle or the off-center
!
!     x, y       vertex coodinates of the triangle
!     alpha      angles of the triangle
!     minloc     index of the triangle vertex with smallest angle
!     betamax    maximum allowed radius edge ratio
!                radius edge ratio = circumcenter radius / shortest edge
!                shortest edge is opposite to smallest angle
!     xpa,ypa    coordinates of the off-center or circumcenter
!
!  local variables
      real(DP) :: beta, beta2, betanew, delta
      integer(I4B) :: isucc(3), ipred(3), id
      parameter (isucc=(/2,3,1/))  !  successor
      parameter (ipred=(/3,1,2/))  !  predecessor
!  the algorithm computes off-centers as proposed in:
!  Uengoer, Alpen: "Off-centers: A new type of Steiner points for computing size-optimal 
!                 quality-guaranteed Delaunay triangulations", 
!  Proc. Latin American Symposium on Theoretical Informatics (2004), pp. 152--161
!
      id=minloc
!  beta is the actual radius edge ratio of the triangle
      beta = 0.5_DP / sin(alpha(id))
!  beta2 is the radius edge ratio of the new triangle
!        with the shortest edge of the old triangle and the circumcenter
      if (beta .ge. 1._DP) then
!  Case A: smallest angle is below 30 deg
!          shortest edge of the new triangle is shortest edge of the old triangle
        beta2 = beta**2 / sqrt(4._DP*beta**2-1._DP)
!          use circumcenter the new triangle if it is not skinny
!          otherwise use off-center such that new triangle has radius edge ratio = betamax
        betanew=min(beta2,betamax)
      else
!  Case B: smallest angle is larger than 30 deg
!          shortest edge of the new triangle is smaller than that of the old triangle
!          i.e circumradius would be smaller than the smallest edge of the old triangle
!        beta2 = beta    / sqrt(4._DP*beta**2-1._DP)
!  we anyhow use the same construction as in case A
        beta2 = beta**2 / sqrt(4._DP*beta**2-1._DP)
        betanew=min(beta2,betamax)
      end if
      delta = betanew + sqrt(betanew**2-0.25_DP)
!  coordinates of off-center or circumcenter
      xpa = (x(isucc(id))+x(ipred(id)))/2._DP + delta*(y(isucc(id))-y(ipred(id)))
      ypa = (y(isucc(id))+y(ipred(id)))/2._DP - delta*(x(isucc(id))-x(ipred(id)))
      return
      end
