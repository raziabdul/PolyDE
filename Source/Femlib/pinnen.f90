      logical function pinnen(bzi,start,ende,zki,xbk,ybk,x,y)
      use feminterface, only: wink, inkrp
      use femtypes
      implicit none
      integer (I4B) start,ende,bzi(:),zki(:,:)
      real (DP) xbk(:),ybk(:),x,y
      intent (in) :: bzi,start,ende,zki,xbk,ybk,x,y
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
!    $Revision: 1.4 $
!    $Date: 2006/07/07 21:08:36 $
!    $Author: r_abdul $
!
!  Input:
!            bzi      Bereich-Zweig-Information Edges of the Domains
!            start    Startlocation (in bzi) of the edges of the actual Domain
!            ende     Endlocation (in bzi) of the edges in the actual Domain
!            zki      Edge Node-Information, Geometry-Nodes of the Edges
!            xbk,ybk  Co-ordinates of the Geometry-Nodes
!            x,y      Co-ordinates of the Point to be tested
!
!  Test whether the Point (x,y) is located inside the Polygon which is given 
!  by the List of Branches bzi
!
!  local variables
!
      real (DP) pi
      parameter (pi=3.14159265358979323846264338327950288419716_DP)
      integer (I4B) i,j,k,dir1
      real (DP) x1,y1,x2,y2,x3,y3,xm,ym,xi,yi,rinnen
      real (DP) sumwi,wi,wink1,flae,radius
!
      sumwi=0._DP
      do j=start,ende
        i=bzi(j)
        k=abs(i)
!  Gerade Berandung
        if (i.gt.0) then
          x1=xbk(zki(1,k))
          y1=ybk(zki(1,k))
          x2=xbk(zki(2,k))
          y2=ybk(zki(2,k))
          dir1=1
        else
          x2=xbk(zki(1,k))
          y2=ybk(zki(1,k))
          x1=xbk(zki(2,k))
          y1=ybk(zki(2,k))
          dir1=-1
        end if
        wi=wink(x,y,x1,y1,x2,y2)
        if (zki(3,k).ne.0) then
!        write (*,*) 'kreisbogen'
          if (zki(3,k).lt.0) then
            xm=xbk(-zki(3,k))
            ym=ybk(-zki(3,k))
            radius=sqrt((xbk(zki(1,k))-xm)**2 + (ybk(zki(1,k))-ym)**2)
          else
            call inkrp(xbk(zki(1,k)),ybk(zki(1,k)),xbk(zki(2,k)),       &
     &        ybk(zki(2,k)),xbk(zki(3,k)),ybk(zki(3,k)),                &
     &        xi,yi,flae,rinnen,xm,ym,radius)
          end if
!  Beim Kreisbogen werden die Koordinaten eines Punkts (x3,y3)
!  vermerkt, der auf der Tangente durch den aktuellen Punkt
!  liegt.
          if (dir1.eq.1) then
            y3=ybk(zki(1,k))+(xbk(zki(1,k))-xm)
            x3=xbk(zki(1,k))-(ybk(zki(1,k))-ym)
          else
            y3=ybk(zki(2,k))+(xbk(zki(2,k))-xm)*(-1._DP)
            x3=xbk(zki(2,k))-(ybk(zki(2,k))-ym)*(-1._DP)
          end if
!  Wenn der Punkt (x,y) innerhalb des Kreisbogens liegt,
!  ist eine besondere Pruefung erforderlich, um den
!  Blickwinkel zu ermitteln
          if (sqrt((x-xm)**2+(y-ym)**2).lt.radius) then
!         write (*,*) 'radiustest positiv'
!  ermittle den Winkel, unter dem die Verbindung weitergeht
            wink1=wink(x1,y1,x3,y3,x2,y2)
!         wink2=wink(x1,y1,x,y,x2,y2)
            if (wink1.le.pi) then
!  Gar nichts
            else
              wi=wi-2._DP*pi
            end if
          else
!  Punkt ausserhalb des Kreisbogens?
!  Dann einfach wie Geradenstueck behanden.
            if (wi.gt.pi) wi=wi-2._DP*pi
          end if
        else
!  Wenn Gerade, dann:
          if (wi.gt.pi) wi=wi-2._DP*pi
        end if
        sumwi=sumwi+wi
      end do
      if (sumwi.gt.pi) then
        pinnen=.true.
      else
        pinnen=.false.
      end if
      return
      end function pinnen
