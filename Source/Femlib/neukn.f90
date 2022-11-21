      subroutine neukn(n1,n2,n3,n4,neuk,xn,yn,x,kzi,zki,p,xbk,ybk,      &
     &  geb1,geb2,arc,err)
      use feminterface, only: posarc
      use femtypes
      implicit none
      integer (I4B) n1, n2, n3, n4, neuk, kzi(:), zki(:,:)
      integer (I4B) p,geb1,geb2, err
      real (DP) xn(:),yn(:),xbk(:),ybk(:)
      complex (DPC) x(:)
      logical arc
      intent (in) :: n1, n2, n3, n4, zki, xbk, ybk, geb1, geb2
      intent (out) :: neuk, arc, err
      intent (inout) :: xn, yn, x, kzi, p
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
!    $Date: 2006/07/07 21:08:36 $
!    $Author: r_abdul $
!
!  Input:
!            n1       node number 1 of the edge to be subdivided
!            n2       node number 2 of the edge to be subdivided
!            n3       third node of element 
!            n4       node of the adjacent element (opposite to n3)
!            zki      branch information, geometry nodes which determine the branches
!            xbk      coordinates of the geometry node
!            ybk
!            geb1     geometry region of first element
!            geb2     geometry region of second element
!
!  Output: 
!            neuk     number of the new node
!            arc      =.true. if the new node had been generated on an arc
!                      .false. for a staight line
!            err      error message obtained from subroutine posarc
!
!  In-/ Output: 
!            xn, yn   Coordinates of the nodes
!            x        solution vector
!            kzi      node branch information
!            p        total number of nodes
!
!  local variables
      integer (I4B) neukzi
      real (DP) xneu,yneu,x1,y1,x2,y2, x3, y3, x4, y4
      complex (DPC) aneu, a1, a2
!
!  Unterprogramm Neukn: generiere neuen Knoten
!
!  Standardwerte:
      neukzi=0
!  Knoten liegt auf der Kante zwischen zwei Gebieten?
      if (geb1.ne.geb2) then
!  Beide Gebiete unterschiedlich, also durch Zweig getrennt
!  Betrachte die beiden Eckknoten n1 und n2 des zu
!  verfeinernden Elements
        if (kzi(n1).eq.kzi(n2)) then
!  beide Punkte auf dem selben Zweig? Kein Problem.
          neukzi=kzi(n1)
        else
!  Ist einer der beiden Punkte ein Zweigendpunkt?
          if ((kzi(n1).lt.0).or.(kzi(n2).lt.0)) then
!  Ist es nicht n1?
            if (kzi(n1).gt.0) then
!  Dann muss es n2 sein! Aber: Vertrauen ist gut, Kontrolle ist besser.
              if ((zki(1,kzi(n1)).eq.-kzi(n2)).or.                      &
     &          (zki(2,kzi(n1)).eq.-kzi(n2))) then
                neukzi=kzi(n1)
              end if
            else
!  entsprechendes gilt fuer n2.
              if (kzi(n2).gt.0) then
                if ((zki(1,kzi(n2)).eq.-kzi(n1)).or.                    &
     &            (zki(2,kzi(n2)).eq.-kzi(n1))) then
                  neukzi=kzi(n2)
                end if
              else
!  Beide Knoten koennen nicht ZweigEndpunkte sein, da immer
!  3 Knoten pro Zweig vorausgesetzt werden
              end if
            end if
          end if
        end if
        if (neukzi.eq.0) then
          write (*,765)
765       format('   ***** Fehler in neukn, Zweig zwischen 2 Elementen')
          write (*,766)
766       format('       aus unterschiedlichen Gebieten nicht gefunden')
        end if
      end if
!  Interpoliere die neue Position
      x1=xn(n1)
      y1=yn(n1)
      x2=xn(n2)
      y2=yn(n2)
      a1=x(n1)
      a2=x(n2)
!  Liegt der neue Knoten auf einem kreisfoermigen Zweig?
      if (neukzi.ne.0) then
        if (zki(3,neukzi).ne.0) then
          arc=.true.
!  Dann Sonderpositionierung aufrufen
          x3=xn(n3)
          y3=yn(n3)
          if (n4 .gt. 0) then
!  neighbor is present
            x4=xn(n4)
            y4=yn(n4)
            call posarc(xbk,ybk,zki,neukzi,xneu,yneu,x1,y1,x2,y2,       &
     &        x3,y3,x4,y4,.true.,err)
          else
!  no neighbor of the element
            call posarc(xbk,ybk,zki,neukzi,xneu,yneu,x1,y1,x2,y2,       &
     &        x3,y3,x4,y4,.false.,err)
          end if
          if (err .ne. 0) then
            print*,'*** error in posarc',err
            return
          end if
        else
          arc=.false.
!  Sonst Mitte der Kante waehlen.
          xneu=(x1+x2)/2.0_DP
          yneu=(y1+y2)/2.0_DP
        end if
      else
        arc=.false.
!  Sonst Mitte der Kante waehlen.
        xneu=(x1+x2)/2.0_DP
        yneu=(y1+y2)/2.0_DP
      end if
!  Fuer den Loesungsvektor die einfachstmoegliche Startannahme
      aneu=(a1+a2)/2.0_DP
!
!  Jetzt liegt der neue Knoten mit all seinen Eigenschaften vor!
!  xneu,yneu: neue Positionen
!  aneu:      Schaetzwert fuer Vektorpotential
!  neukzi:    neue Knoten-Zweig-Information
!
!  Jetzt die Programmdaten ergaenzen
      p=p+1
!  Daten fuer den neuen Knoten setzen
      neuk=p
      xn(neuk)=xneu
      yn(neuk)=yneu
      x(neuk)=aneu
      kzi(neuk)=neukzi
      err=0
      return
      end subroutine neukn
!
!
!
      subroutine posarc(xbk,ybk,zki,neukzi,xneu,yneu,x1,y1,x2,y2,       &
     &  x3,y3,x4,y4,neighb,err)
      use feminterface, only: inkrp, snkrgr
      use femtypes
      implicit none
      integer (I4B) zki(:,:),neukzi, err
      real (DP) xbk(:), ybk(:), xneu, yneu
      real (DP) x1, y1, x2, y2, x3, y3, x4, y4
      logical neighb
      intent (in) :: xbk, ybk, zki,neukzi, x1, y1, x2, y2
      intent (in) :: x3, y3, x4, y4, neighb
      intent (out) :: xneu, yneu, err
!
!  Determine the position of a new node on an arc
!
!  Input:
!            xbk      coordinates of the key-points
!            ybk
!            zki      branch information, key-points which determine the branches
!            neukzi   for each of the new nodes the node-edge-information is stored
!                     input branches of the nodes
!            x1,y1    coordinates of node 1 of the element edge
!            x2,y2    coordinates of node 2 of the element edge
!            x3,y3    third node of the element
!            x4,y4    third node of the neighboring element
!            neighb   =.true. if there is a neighbor present opposite to node 3
!
!  Output: 
!            xneu     coordinates of the new nodes, if err <> 0 the coordinates are 
!            yneu     either not computed or probably give no reliable results 
!            err    = 0 if no error occured
!                     1 if no point on the arc can be found which is at the same time 
!                       visible from point 3 and point 4
!                     3 if there is no point on the arc, which is visible from point 3
!                       i.e. the point 3 is located inside the arc
!                     4 if there is no point on the arc, which is visible from point 4
!                       i.e. the point 4 is located inside the arc
!                     9 The new point probably is not a good candidate, since it is very 
!                       near to one of the element vertices, either (x1,y1) or (x2,y2)
!                    10 an error in snkrgr probably either node 1 or 2 are idential with 3 or 4
!                   100 computation would not give reliable results, since the radius to the first
!                       and the second point differ stongly, i.e. (r1^2-r2^2)/(r1^2+r2^2) > eps!                      
!
!  local variables
      real (DP) r1, r2, t1, t2, swp
      real (DP) xi,yi, xa, ya, xb, yb
      real (DP) flae2, ri, rho
      real (DP) pi, zweipi, xctr, yctr
      real (DP) ta,ta2, tb, tb2, tc, tc2, td, td2, tmin, tmax, tnew
      real (DP) tmin1, tmin2, tmax1, tmax2
      real (DP) radius2, wa, wba, zaehl, xnenn
      integer (I4B) err1
      parameter (pi=3.14159265358979323846264338327950288419716_DP)
      parameter (zweipi=2._DP*pi)
!
      xa=xbk(zki(1,neukzi))
      ya=ybk(zki(1,neukzi))
      xb=xbk(zki(2,neukzi))
      yb=ybk(zki(2,neukzi))
!  compute the location of the center
      if (zki(3,neukzi).lt.0) then
!  center is already known
        xctr=xbk(-zki(3,neukzi))
        yctr=ybk(-zki(3,neukzi))
      else
!  output third node on the arc
        call inkrp(x1,y1,x2,y2,xbk(zki(3,neukzi)),ybk(zki(3,neukzi)),   &
     &    xi,yi,flae2,ri,xctr,yctr,rho)
      end if
!
      if (abs((x1-xctr)**2+(y1-yctr)**2-(x2-xctr)**2-(y2-yctr)**2)/     &
     &    ((x1-xctr)**2+(y1-yctr)**2+(x2-xctr)**2+(y2-yctr)**2)         &
     &    .gt. 1.e-4_DP) then
!  radius to the first and the second point differ stongly (r1^2-r2^2)/(r1^2+r2^2)
        xneu=(x1+x2)/2.0_DP
        yneu=(y1+y2)/2.0_DP
        err=100
        return
      end if
!  identify the portion of the arc between the rays 3->1 and 3->2, visible from point 3
!  compute the intersection of the arc with with the line 3->1
      call snkrgr(x3,y3,x1,y1,xa,ya,xb,yb,xctr,yctr,r1,r2,t1,t2,err1)
      if (err1 .gt. 0) then
        err=err1
        return
      end if
      if (abs(r1-1._DP) .lt. abs(r2-1._DP)) then
!  the point 1 is located on the circle, the second intersection can be at an aribitrary distance
        ta=t1
        ta2=t2
      else
        ta=t2
        ta2=t1
      end if
!  compute the intersection of the arc with with the line 3->2
      call snkrgr(x3,y3,x2,y2,xa,ya,xb,yb,xctr,yctr,r1,r2,t1,t2,err1)
      if (err1 .gt. 0) then
        err=err1
        return
      end if
      if (abs(r1-1._DP) .lt. abs(r2-1._DP)) then
!  the point 2 is located on the circle, the second intersection can be at an aribitrary distance
        tb=t1
        tb2=t2
      else
        tb=t2
        tb2=t1
      end if
!  interchange start and end if end is smaller
      if (ta .gt. tb) then
        swp=ta
        ta=tb
        tb=swp
        swp=ta2
        ta2=tb2
        tb2=swp
      end if
!  now determine the visible fraction
      if ((ta2 .gt. ta) .and. (ta2 .lt. tb)) then
        tmin1=ta2
      else
        tmin1=ta
      end if
      if ((tb2 .gt. ta) .and. (tb2 .lt. tb)) then
        tmax1=tb2
      else
        tmax1=tb
      end if
      if (tmax1 .le. tmin1) then
!  if this happens, it means that there is no valid solution 
!  there is no continous fraction of the arc between the rays 3->1 and 3->2
!  the point 4 is located inside the circle'
        xneu=(x1+x2)/2.0_DP
        yneu=(y1+y2)/2.0_DP
        err=3
        return
      end if
!
!  identify the portion of the arc between the rays 4->1 and 4->2, visible from point 4
      if (neighb) then
!  compute the intersection of the arc with with the line 4->1
        call snkrgr(x4,y4,x1,y1,xa,ya,xb,yb,xctr,yctr,r1,r2,t1,t2,err1)
        if (err1 .gt. 0) then
          err=err1
          return
        end if
        if (abs(r1-1._DP) .lt. abs(r2-1._DP)) then
!  the point 3 is located on the circle, the second intersection can be at an aribitrary distance
          tc=t1
          tc2=t2
        else
          tc=t2
          tc2=t1
        end if
!  compute the intersection of the arc with with the line 4->2
        call snkrgr(x4,y4,x2,y2,xa,ya,xb,yb,xctr,yctr,r1,r2,t1,t2,err1)
        if (err1 .gt. 0) then
          err=err1
          return
        end if
        if (abs(r1-1._DP) .lt. abs(r2-1._DP)) then
!  the point 4 is located on the circle, the second intersection can be at an aribitrary distance
          td=t1
          td2=t2
        else
          td=t2
          td2=t1
        end if
!  interchange start and end, if end is smaller
        if (tc .gt. td) then
          swp=tc
          tc=td
          td=swp
          swp=tc2
          tc2=td2
          td2=swp
        end if
!  now determine the visible fraction
        if ((tc2 .gt. tc) .and. (tc2 .lt. td)) then
          tmin2=tc2
        else
          tmin2=tc
        end if
        if ((td2 .gt. tc) .and. (td2 .lt. td)) then
          tmax2=td2
        else
          tmax2=td
        end if
      else
        tmin2=0._DP
        tmax2=1._DP
      end if
      if (tmax2 .le. tmin2) then
!  if this happens, it means that there is no valid solution 
!  there is no continous fraction of the arc between the rays 4->1 and 4->2
!  the point 4 is located inside the circle'
        xneu=(x1+x2)/2.0_DP
        yneu=(y1+y2)/2.0_DP
        err=4
        return
      end if
!
!  now compute the portion of the arc, which is visible from both sides
      tmin=max(tmin1,tmin2)
      tmax=min(tmax1,tmax2)
      if (tmin .le. tmax) then
        tnew=(tmin+tmax)/2._DP
        err=0
      else
!  there is no common visibility are for the arc from both sides
        tnew=(tmin+tmax)/2._DP
        xneu=(x1+x2)/2.0_DP
        yneu=(y1+y2)/2.0_DP
        err=1
        return
      end if
!
      if (abs((tnew-ta)/(tb-ta)-.5_DP).gt.0.45_DP) then
!  the new point on the arc is very near (5%) to one of the element vertices
        err=9
      end if
!
!  radius as the mean value of the two points
      radius2=((xa-xctr)**2+(ya-yctr)**2+                               &
     &         (xb-xctr)**2+(yb-yctr)**2)/2._DP
!
!  angle between the points 1 and 2 seen from the center
      zaehl=(xa-xctr)*(yb-yctr)-(xb-xctr)*(ya-yctr)
      xnenn=(xb-xctr)*(xa-xctr)+(yb-yctr)*(ya-yctr)
      wba=atan2(zaehl,xnenn)
      if (wba .le. 0._DP) wba=wba+zweipi
      wa=atan2(ya-yctr,xa-xctr)
      xneu=xctr+sqrt(radius2)*cos(wa+tnew*wba)
      yneu=yctr+sqrt(radius2)*sin(wa+tnew*wba)
      return
      end subroutine posarc
