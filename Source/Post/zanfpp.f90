      subroutine zanfpp(plottyp,xmin,xmax,ymin,ymax,h,bildnr,information)
      use feminterface, only: getpostsetting
      use femtypes
      implicit none
      real (DP) xmin,xmax,ymin,ymax,h
      integer (I4B) bildnr
      character (len=*) plottyp
      character (len=*),optional :: information
      intent (in):: plottyp, information
      intent (inout):: xmin,xmax,ymin,ymax
      intent (out):: bildnr,h
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
!    $Revision: 1.17 $
!    $Date: 2010/08/27 12:39:33 $
!    $Author: m_kasper $
!
!  Beginn einer Zeichnung oder einer neuen Teilzeichnung
!      plottyp  Type des Graphik Fensters
!               /WV,    640 x  480   Pixel Windows Grafik Fenster
!               /WS     800 x  600   Pixel Windows Grafik Fenster
!               /WX,   1024 x  768   Pixel Windows Grafik Fenster
!               /WZ    1280 x 1024   Pixel Windows Grafik Fenster
!               /PS    Postscript (monochrome landscape mode, long edge of paper horizontal).
!               /CPS   Postscript (color landscape mode, long edge of paper horizontal).
!               /VPS   Postscript (monochrome portrait mode, short edge of paper horizontal).
!               /VCPS  Postscript  (color portrait mode, short edge of paper horizontal).
!               /HPGL2 Any HPGL-2 device (presently tested only on HP laserjet 3)
!               /NULL  Null Device
!               /HPGL  (Hewlett Packard HPGL plotter, landscape orientation)')
!               /VHPGL (Hewlett Packard HPGL plotter, portrait orientation)')
!
!               /GIF   850 x 680 pixels (translates to 10.0 x  8.0 inch) Graphics Interchange Format (GIF) files.
!               /VGIF  680 x 850 pixels (translates to  8.0 x 10.0 inch) Graphics Interchange Format (GIF) files.
!
!      xmin     x-Koordnate der linken unteren Ecke
!      xmax     x-Koordnate der rechten oberen Ecke
!      ymin     y-Koordnate der linken unteren Ecke
!      ymax     y-Koordnate der rechten oberen Ecke
!      h        Hoehe der Schrift  in normierten Einheiten
!  
!      text     Beliebiger Text, der in die linke obere Ecke
!               der Zeichnung geschrieben wird
!      a,b      Groesse des Rahmens in cm
!      mx,my    Massstabsgfaktoren in x bzw. y Richtung (Dehnungsfaktoren)
!      x0, y0   Koordinaten des Ursprungsauf dem Zeichnungspapier
!
!      nach dem Aufruf stehen die folgenden Groessen zur Verfuegung
!      xmx, xmy aktuelle Massstabsfaktoren in x bzw. y Richtung
!
      integer (I4B) bildnri,istat,pgbeg,pgopen
      real (SP) margin,rxmin,rxmax,rymin,rymax,height
      logical first, nfo
      data bildnri/0/
      data first/.true./
!
!
      external :: grglun, pgbeg, pgbox, pgcurs, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgslw, pgqch, pgqcr
      external :: pgopen, pgpap, pgqwin, pgsvp, pgwnad
!
      if (.not. present(information)) then
        nfo = .false.
      else if (information .eq. 'YES') then
        nfo = .true.
      else
        nfo = .false.
      end if
      if (first) then
        istat = pgbeg(0,plottyp,1,1)
        first=.false.
      else
      istat = pgopen(plottyp)
      end if
      bildnri=bildnri+1
      bildnr=bildnri
!  Set Viewport with margin (0.01 == 1%)
      margin=.08
      if((plottyp .eq. '/VCPS') .or. (plottyp .eq. '/VPS')) then                         
        call pgpap(7.7,1.43)
        if(nfo) then
          call pgsvp(margin, 1.-margin, 0.3, 1.-margin)
        else
          call pgsvp(margin , 1.-margin -0.05, margin, 1.-margin)
        endif
      else
        call pgpap(11.0,0.7)
        if(nfo) then
          call pgsvp(margin , 0.7-margin-0.05, margin, 1.-margin-0.03)
        else
          call pgsvp(margin , 1.-margin-0.05, margin+0.08, 1.-margin-0.03)
        endif
      endif 
!  plot simulation information if needed
!  use equal scaling for x and y and adjust viewport
      rxmin=real(xmin)
      rxmax=real(xmax)
      rymin=real(ymin)
      rymax=real(ymax)
      call pgwnad (rxmin,rxmax,rymin,rymax)
      call pgqwin(rxmin,rxmax,rymin,rymax)     
      xmin=dble(rxmin)
      xmax=dble(rxmax)
      ymin=dble(rymin)
      ymax=dble(rymax)
!  Set default character Font (1=simple, 2=Roman)
      call pgscf(1)
!  Character height 
!  default 1/40th of the height or width of the view surface (whichever is less); 
      call pgqch(height)
      h=dble(height)
      return
      end subroutine zanfpp
