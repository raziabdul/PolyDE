      subroutine nfarbe(z,x,y,e,ie,palette,ncolor,zcolmin,zcolmax)
      use feminterface, only : drtri
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:), zcolmin, zcolmax
      integer (I4B) e(:,:), ie, ncolor, palette(:)
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
!  Falschfarbendarstellung ueber einem Dreiecksnetz  2D Plot
!  Faebung entsprechend den Hoehenwerte z in den Dreickseckpunkten
!  Interpolation der Farben auf dem Dreieck
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
      real (DP), allocatable :: zcolor(:)
      real (SP) ocr, ocg, ocb
      integer (I4B) i, oldindex
      logical proj
!
      external :: pgqci, pgqcr, pgsfs, pgsci, pgscr

!  Farbstatus retten
      call pgqci(oldindex)
      call pgsci(20)
      call pgqcr(20, ocr, ocg, ocb)
!  Set Filling style
      call pgsfs(1)

!  Zuordnung der Farbwerte zu den Hoehenwerten
      allocate(zcolor(ncolor))
      do i=1,ncolor
        zcolor(i)=(zcolmax-zcolmin)/dble(ncolor)*dble(i)+zcolmin
      end do
      proj= .false.
!  Schleife ueber alle Elemente

      do i=1,ie
        call drtri(i,e,x,y,z,zcolor,ncolor,palette,proj)
      end do
      deallocate(zcolor)
!  Farbstatus zuruecksetzen
      call pgsci(oldindex)
      call pgscr(20, ocr, ocg, ocb)
      return
      end subroutine nfarbe
