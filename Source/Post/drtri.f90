      subroutine drtri(i,e,x,y,z,zcolor,ncolor,palette,proj,            &
     &  xo,yo,zo,transf)
      use feminterface, only : locate, divtri, sptransformpoly
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:), zcolor(:)
      real (DP), optional ::  xo(:), yo(:), zo(:), transf(4,4)
      integer (I4B) i, e(:,:), ncolor, palette(:)
      logical proj
      intent(in):: i, e, x, y, z, xo, yo, zo, zcolor, ncolor, palette,  &
     &  transf, proj
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
!    $Date: 2006/07/07 22:13:33 $
!    $Author: r_abdul $
!
!  Falschfarbendarstellung eines Dreiecks
!  Faebung entsprechend der Hoehenwerte z in den Dreickseckpunkten
!  Interpolation der Farben auf dem Dreieck
!  und optional Projektion zur Darstellung der Oberflaechentopologie
!
!     i        Nummer des Dreiecks
!     e        Elementinformation Knoten der Dreiecke
!     x        x Koordinaten der Dreicksecken
!     y        y Koordinaten der Dreicksecken
!     z        Hoehenwerte in den Dreieckseckpunkte
!     zcolor   Hoehenwerte der Faebstufen 
!     ncolor   Anzahl der Farbwerte
!     palette  Farbwerte der in aufsteigender Reihenfolge
!     proj     .true. fuer einen 3D Plot 
!              .false. fuer einen 2D Plot
!     xo       x Koordinaten der Dreiecksecken im transformierten System
!     yo       y Koordinaten der Dreiecksecken im transformierten System
!     zo       Hoehenwerte der Dreiecksecken im transformierten System
!     transf   Tranformationsmatrix (Homogene Koordinaten)
!
      integer (I4B) imin(1), imax(1), icolor, npkt
      real (DP) xl(3), yl(3), zl(3), normale
      real (DP) zmin, zmax
      integer (I4B) rgb, zff
      real (SP) xpoly(5), ypoly(5), zpoly(5)
      real (SP) xpolyout(5), ypolyout(5), zpolyout(5), cr, cg, cb
      logical done, first, all
      parameter (zff=Z'FF')

!  lokale Kopie der Eckkoordinaten und Hoehenwerte
      xl(:)=x(e(:,i))
      yl(:)=y(e(:,i))
      zl(:)=z(e(:,i))
!  Berechnung der Flaeche im transformierten System (Flaechennormale zur Betrachterrichtung)
      if (proj) then
        normale=-(xo(e(2,i))-xo(e(1,i)))*(yo(e(3,i))-yo(e(1,i)))+  &
     &           (xo(e(3,i))-xo(e(1,i)))*(yo(e(2,i))-yo(e(1,i)))
!       if (normale .le. 0.) cycle
      end if
!  Ecke mit dem kleinsten und groessten z Wert feststellen
      imin=minloc(zl)
      imax=maxloc(zl)
      zmin=zl(imin(1))
      zmax=zl(imax(1))

!  Erste zu verwendende Farbe feststellen
      icolor=locate(zcolor,zmin)
      icolor=max(1,min(icolor,ncolor))
        
      first=.true.
      done=.false.

      do while (.not. done) 
        if (icolor.ge.ncolor) then 
          all=.true.
        else
          all=.false.
        end if
        call divtri(xl,yl,zl,zcolor(min(ncolor,icolor+1)),first,all,done,npkt,xpoly,ypoly,zpoly)
!
!  Farbe festlegen
        rgb    = palette(icolor)
        cr = float(iand(rgb, zff))/255
        cg = float(iand(ishft(rgb,-8), zff))/255
        cb = float(iand(ishft(rgb,-16), zff))/255
        call pgscr(20, cr, cg, cb)
        if (proj) then
!  Transformation des Polygons
          call sptransformpoly(xpoly,ypoly,zpoly,transf,npkt,xpolyout,ypolyout,zpolyout)
!  Zeichnen des Polygons
          call pgpoly(npkt,xpolyout,ypolyout)
        else
          call pgpoly(npkt,xpoly,ypoly)
        end if
        icolor=icolor+1
      end do
      return 
      end subroutine drtri
!
!
!
      subroutine divtri(xl,yl,zl,zval,first,all,done,npkt,xpoly,ypoly,zpoly)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) xl(3), yl(3), zl(3), zval
      real (SP) xpoly(:), ypoly(:), zpoly(:)
      integer (I4B) npkt
      logical done, first, all
      intent(in):: xl, yl, zl, zval, all
      intent(inout):: first
      intent (out) done, npkt, xpoly, ypoly, zpoly
!  Sukzessive Zerlegung eine Dreiecks (Polygons) in Flaechenanteile 
!  zu den zugehoerigen Hoehenwerten.
!  Beim ersten Aufruf wird die Flaeche zwischen dem Eckpunkt mit dem kleinsten 
!  Hoehenwert und der aktuellen Hoehenlinie ausgegeben. 
!  Bei den folgenden Aufrufen wird jeweils ein Polygon zurueckgeliefert, 
!  das die umschlossene Flaeche zwischen dem letzten Aufruf und der 
!  Hoehenlinie zval enthaelt. Das ausgegebene Polygon enthaelt 3-5 Ecken.
!  Das Polygon ist im mathematisch positiven Sinn orientiert
!
!      xl         x-Koordinaten des  Dreiecks
!      yl         x-Koordinaten des  Dreiecks
!      zl         Hoehenwerte in den Ecken des Dreiecks
!      zval       aktuelle Hoehenwert fuer den das Polygon auszugeben ist
!      first      muss beim ersten Aufruf eines Polygons auf true und bei allen 
!                 folgenden Aufrufen fuer das gleiche Ploygon auf false gesetzt sein
!      all        kann auf true gesetzt werden  um den Rest des Dreieck zu bearbeiten
!      done       wird auf true gesetzt, wenn das Dreieck vollstaendig bearbeitet wurde
!      npkt       Anzahl der Ecken des ausgegebenen Polygons
!      xpoly      x-Koordinaten des Polygons
!      ypoly      y-Koordinaten des Polygons
!      zpoly      z-Koordinaten (Hoehenwerte) des Polygons
!                 
      integer (I4B) plinks, prechts, nachfolger(3), vorganger(3),imin(1)
      real (DP) xlinks, ylinks,zlinks, xrechts, yrechts,zrechts,alpha
      logical mark
      save xlinks, ylinks, zlinks, xrechts, yrechts, zrechts
      save prechts, plinks
      parameter (nachfolger=(/2,3,1/))
      parameter (vorganger=(/3,1,2/))
!
      if (first) then
!  Startpunkt auf der Ecke mit dem kleisten Hoehenwert und der niedrigsten Farbe     
        npkt=1
        imin=minloc(zl)
        plinks=imin(1)
        prechts=imin(1)
        done=.false.
        xpoly(npkt)=xl(imin(1))
        ypoly(npkt)=yl(imin(1))
        zpoly(npkt)=zl(imin(1))
        first=.false.
      else 
!  Wenn Dreieck noch nicht fertig bearbeitet naechstes Polygon vorbereiten
        npkt=2
        xpoly(1)=xlinks
        ypoly(1)=ylinks
        zpoly(1)=zlinks
        xpoly(2)=xrechts
        ypoly(2)=yrechts
        zpoly(2)=zrechts
      end if

!  Rest des Dreiecks bearbeiten
        if (all .or. (maxval(zl) .le. zval)) then
          done=.true.
!  Verbleibende Eckpunkte in das  Polygon aufnehmen
          do while ((nachfolger(prechts)) .ne. plinks )
            npkt=npkt+1
            xpoly(npkt)=xl(nachfolger(prechts))
            ypoly(npkt)=yl(nachfolger(prechts))
            zpoly(npkt)=zl(nachfolger(prechts))
            prechts=nachfolger(prechts)
          end do
        else 
!  Rechts weiterverfolgen
          if (zl(nachfolger(prechts)) .le. zval) then
!  naechste Ecke reckts gehoert zum Polygon
            npkt=npkt+1
            xpoly(npkt)=xl(nachfolger(prechts))
            ypoly(npkt)=yl(nachfolger(prechts))
            zpoly(npkt)=zl(nachfolger(prechts))
            prechts=nachfolger(prechts)
          end if
!  Schnittpunkt der Hoehenlinie mit der Dreiecksseite rechts berechnen 
          alpha=(zval-zl(prechts))/(zl(nachfolger(prechts))-zl(prechts))
          xrechts=xl(prechts) + alpha*(xl(nachfolger(prechts))-xl(prechts))
          yrechts=yl(prechts) + alpha*(yl(nachfolger(prechts))-yl(prechts))
          zrechts=zval
          npkt=npkt+1
          xpoly(npkt)=xrechts
          ypoly(npkt)=yrechts
          zpoly(npkt)=zrechts
    
!  Links weiterverfolgen
          if (zl(vorganger(plinks)) .le. zval) then
!  naechste Ecke links gehoert zum Polygon
            xpoly(npkt+2)=xl(vorganger(plinks))
            ypoly(npkt+2)=yl(vorganger(plinks))
            zpoly(npkt+2)=zl(vorganger(plinks))
            plinks=vorganger(plinks)
            mark=.true.
          else 
            mark=.false.
          end if
!  Schnittpunkt der Hoehenlinie mit der Dreiecksseite links berechnen 
          alpha=(zval-zl(plinks))/(zl(vorganger(plinks))-zl(plinks))
          xlinks=xl(plinks) + alpha*(xl(vorganger(plinks))-xl(plinks))
          ylinks=yl(plinks) + alpha*(yl(vorganger(plinks))-yl(plinks))
          zlinks=zval
          xpoly(npkt+1)=xlinks
          ypoly(npkt+1)=ylinks
          zpoly(npkt+1)=zlinks
!
          if (mark) then 
            npkt=npkt+2
          else
            npkt=npkt+1
          end if
        end if
      return
      end subroutine divtri
