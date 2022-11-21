! 
      subroutine extend(anzzwg,xpoint,ypoint,zki,                       &
     &                  xmin,xmax,ymin,ymax)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) xpoint(:), ypoint(:), xmin, xmax, ymin, ymax
      integer (I4B) anzzwg, zki(:,:)
      intent (in) anzzwg, xpoint, ypoint, zki 
      intent (out) xmin, xmax, ymin, ymax
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
!    $Revision: 1.6 $
!    $Date: 2014/07/15 13:04:51 $
!    $Author: m_kasper $
!
!  Berechnung der Extends aller Zweige
!
!  Eingabe:
!    anzzwg    Anzahl der Zweige
!    xpoint    x-Koordinaten der Knoten
!    ypoint    y-Koordinaten der Knoten
!    zki       Zweig-Knoten Information
!    zki(1,i)  Anfangspunkt des Zweiges
!    zki(2,i)  Endpunkt des Zweiges
!    zki(3,i)  Dritter Punkt des Zweiges
!              /  = 0  :  Gerade
!    zki(3,i)  -  > 0  :  Punkt auf dem Kreisbogen
!              \  < 0  :  Mittelpunkt des Kreisbogens
!  Ausgabe:
!    xmin
!    xmax      Extends (Maximalabmessungen) aller Zeichnungselemente
!    ymin
!    ymax
!
!  Bugs : es fehlt die Implementierung fuer den Fall zki(3,i) > 0
!
      integer (I4B) i
      real (DP) radius, rad1, rad2, pid180, xmid, ymid, phi1, phi2
      data pid180 /.017453292519943295769_DP/
!
!  Phase 1     Extends der Geraden untersuchen
!
      xmax = maxval( xpoint(:) )
      xmin = minval( xpoint(:) )
      ymax = maxval( ypoint(:) )
      ymin = minval( ypoint(:) )
!
!  Phase 2     Extends der Kreisboegen untersuchen
!
      do i=1,anzzwg
!
!  nur fuer Kreisboegen
!
        if (zki(3,i).lt.0) then
!  Kreisbogen mit Angabe des Mittelpunktes
!  Radius
          rad1 = sqrt( (xpoint(zki(1,i))-xpoint(abs(zki(3,i))))**2 +    &
     &                  (ypoint(zki(1,i))-ypoint(abs(zki(3,i))))**2 )
          rad2 = sqrt( (xpoint(zki(2,i))-xpoint(abs(zki(3,i))))**2 +    &
     &                  (ypoint(zki(2,i))-ypoint(abs(zki(3,i))))**2 )
          radius = ( rad1 + rad2 ) / 2._DP
!  x-Koordinate des Mittelpunkt
          xmid = xpoint(abs(zki(3,i)))
!  y-Koordinate des Mittelpunkt
          ymid = ypoint(abs(zki(3,i)))
!  Liegt der Kreis vollstaendig innerhalb der Extends
          if ( ((xmid+radius) .le. xmax) .and.                          &
     &         ((xmid-radius) .ge. xmin) .and.                          &
     &         ((ymid+radius) .le. ymax) .and.                          &
     &         ((ymid-radius) .ge. ymin) ) cycle
!  Anfangswinkel in Grad
!     antan liefert werte im Bereich von (-pi bis pi)
          phi1 = atan2(ypoint(zki(1,i))-ypoint(abs(zki(3,i))) ,         &
     &                  xpoint(zki(1,i))-xpoint(abs(zki(3,i))) )/pid180
!  Endwinkel in Grad
          phi2 = atan2(ypoint(zki(2,i))-ypoint(abs(zki(3,i))) ,         &
     &                  xpoint(zki(2,i))-xpoint(abs(zki(3,i))) )/pid180
          if (phi2.le.phi1) phi1=phi1-360._DP
!
!    0 Grad im Bereich des Kreisbogens ?
          if ( (phi1 .le.   0._DP) .and. (phi2 .ge.   0._DP) ) then
            xmax = max( xmax , xmid+radius )
          end if
!   90 Grad im Bereich des Kreisbogens ?
          if ( (phi1 .le.  90._DP) .and. (phi2 .ge.  90._DP) ) then
            ymax = max( ymax , ymid+radius )
          end if
!  180 Grad im Bereich des Kreisbogens ?
          if ( (phi1 .le. 180._DP) .and. (phi2 .ge. 180._DP) ) then
            xmin = min( xmin , xmid-radius )
          end if
!  -90 Grad im Bereich des Kreisbogens ?
          if ( (phi1 .le. -90._DP) .and. (phi2 .ge. -90._DP) ) then
            ymin = min( ymin , ymid-radius )
          end if
!
        else if (zki(3,i).gt.0) then
!  Kreisbogen mit Angabe eines dritten Punktes
!          print*,' kann Extends nicht berechnen'
        end if
      end do
      return
      end subroutine extend