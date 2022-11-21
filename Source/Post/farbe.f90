      subroutine farbe(z,x,y,e,ie,palette,ncolor,zcolmin,zcolmax)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:), zcolmin, zcolmax
      integer (I4B) e(:,:), ie,ncolor, palette(:)
      intent(in):: z, x, y, e, ie, palette, ncolor, zcolmin, zcolmax
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
!  Falschfarbendarstellung ueber einem Dreiecksnetz 2D Plot
!  Faebung entsprechend den Hoehenwerte z in den Dreicksmittelpunkten
!
!     z        Hoehenwerte in den Dreieckseckpunkte
!     x        x Koordinaten der Dreicksecken
!     y        y Koordinaten der Dreicksecken
!     e        Elementinformation Knoten der Dreiecke
!     ie       Anzahl der Dreiecke
!     palette  Farbwerte der in aufsteigender Reihenfolge
!     ncolor   Anzahl der Farbwerte
!     zcolmin  Hoehenwert zur niedrigsten Farbe
!     zcolmax  Hoehenwert zur hoechsten Farbe
!
      real (DP) zmitte
      real (SP) xpoly(3), ypoly(3), cr, cg, cb, ocr, ocg, ocb
      integer (I4B) i, rgb, oldindex, icolor, zff
      parameter (zff=Z'FF')
!
      external :: pgqci, pgqcr, pgsfs, pgsci, pgscr, pgpoly

!  Farbstatus retten
      call pgqci(oldindex)
      call pgsci(20)
      call pgqcr(20, ocr, ocg, ocb)
!  Set Filling style
      call pgsfs(1)
!  Schleife ueber alle Elemente
      do i=1,ie
        zmitte=z(i)
        icolor=int((zmitte-zcolmin)/(zcolmax-zcolmin)*ncolor)
        if (icolor.lt.1) icolor=1
        if (icolor.gt.ncolor) icolor=ncolor
!  Farbe festlegen
        rgb    = palette(icolor)
        cr = float(iand(rgb, zff))/255
        cg = float(iand(ishft(rgb,-8), zff))/255
        cb = float(iand(ishft(rgb,-16), zff))/255
        call pgscr(20, cr, cg, cb)

        xpoly(1)=x(e(1,i))
        ypoly(1)=y(e(1,i))
        xpoly(2)=x(e(2,i))
        ypoly(2)=y(e(2,i))
        xpoly(3)=x(e(3,i))
        ypoly(3)=y(e(3,i))
        call pgpoly(3,xpoly,ypoly)
      end do
      call pgsci(oldindex)
      call pgscr(20, ocr, ocg, ocb)
      return
      end subroutine farbe
