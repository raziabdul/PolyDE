      subroutine triplt(e,en,x,y,ie,r,g,b)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) x(:), y(:)
      real (SP) r, g, b
      integer (I4B) e(:,:), en(:,:), ie
      intent(in)::e, en, x, y, ie, r, g, b
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
!    $Date: 2006/07/07 22:13:33 $
!    $Author: r_abdul $
!
!  Graphische Ausgabe der Triangulierung
!  Eingabe:  
!       e           Elementinformation
!       en          Nachbarn der Elemente
!       x           x-Koordinaten der Knoten
!       y           y-Koordinaten der Knoten
!       ie          Anzahl der Elemente
!       r, g, b     red green and blue value for the potential lines
!  Ausgabe:         graphisch
      integer (I4B) i, l, m, n
      integer (I4B) oldindex
      real (SP) ocr,ocg,ocb
!
!      
      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgqwin, pgsvp, pgwnad


!  Farbstatus retten
      call pgqci(oldindex)
      call pgsci(20)
      call pgqcr(20, ocr, ocg, ocb)
!  Farbe setzen  Gelb= 1.00, 1.00, 0.00 (R,G,B)
      call pgscr(20,  r, g, b)
!  Zeichnen der Elemente
      do i=1,ie
        l=e(1,i)
        m=e(2,i)
        n=e(3,i)
!  Positionieren des Stiftes nur wenn noetig
!  es werden nur die Zweige gezeichnet, die noch
!  nicht in anderen Elementen gezeichnet wurden
        if (en(3,i).lt.i) then
          if (en(1,i).lt.i) then
            call pgmove(real(x(l)),real(y(l)))
            call pgdraw(real(x(m)),real(y(m)))
            call pgdraw(real(x(n)),real(y(n)))
            if (en(2,i).lt.i) then
              call pgdraw(real(x(l)),real(y(l)))
            end if
          else
            call pgmove(real(x(m)),real(y(m)))
            call pgdraw(real(x(l)),real(y(l)))
            if (en(2,i).lt.i) then
              call pgdraw(real(x(n)),real(y(n)))
            end if
          end if
        else
          if (en(1,i).lt.i) then
            call pgmove(real(x(m)),real(y(m)))
            call pgdraw(real(x(n)),real(y(n)))
            if (en(2,i).lt.i) then
              call pgdraw(real(x(l)),real(y(l)))
            end if
          else
            if (en(2,i).lt.i) then
              call pgmove(real(x(l)),real(y(l)))
              call pgdraw(real(x(n)),real(y(n)))
            end if
          end if
        end if
      end do
      call pgsci(oldindex)
      call pgscr(20, ocr, ocg, ocb)
      return
      end subroutine triplt
