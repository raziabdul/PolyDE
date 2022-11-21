      pure real (DP) function wink(x,y,x1,y1,x2,y2)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) x, y, x1, y1, x2, y2
      intent (in) :: x, y, x1, y1, x2, y2
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
!    $Date: 2006/07/07 21:08:37 $
!    $Author: r_abdul $
!
!  Wink gibt den Winkel zwischen zwei Strecken aus
!  Eingaben: 
!     x,y      gemeinsamer Punkt der beiden Geraden
!     x1,y1    Strecke 1 von x,y nach x1,y1
!     x2,y2    Strecke 2 von x,y nach x2,y2
!  Ausgabe: 
!      Wink    Winkel zwischen Strecke 1 und Strecke 2
!              gemessen am Punkt x,y , beginnend an der Strecke 1
!              Wertebereich: 0 bis 2*pi
!
!  local variables
!
      real (DP) x0,y0,xi1,yi1,xi2,yi2,winkel
      real (DP) pi,zweipi
      parameter (pi=3.14159265358979323846264338327950288419716_DP)
      parameter (zweipi=2._DP*pi)
!
      xi1=x1-x
      yi1=y1-y
      xi2=x2-x
      yi2=y2-y
      x0=xi1*xi2+yi1*yi2
      y0=-yi1*xi2+xi1*yi2
!  Koordinaten des Vektors (xi2,yi2) im xi1,yi1 Koordinatensystem
      if ((x0.eq.0._DP) .and. (y0.eq.0._DP)) then
!       write (*,1)
!1      format(' ****wink: Laenge eines Vektors Null')
        wink=-3._DP*pi
      else
        winkel=atan2(y0,x0)
        if (winkel.lt.0._DP) then
          wink=winkel+zweipi
        else
          wink=winkel
        end if
      end if
      return
      end function wink
!
!
!
      pure real (DP) function wink1(x,y,x1,y1,x2,y2)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) x, y, x1, y1, x2, y2
      intent (in) :: x, y, x1, y1, x2, y2
!
!  Wink gibt den Winkel zwischen zwei Strecken aus
!  Eingaben: 
!     x,y      gemeinsamer Punkt der beiden Geraden
!     x1,y1    Strecke 1 von x,y nach x1,y1
!     x2,y2    Strecke 2 von x,y nach x2,y2
!  Ausgabe: 
!      Wink    Winkel zwischen Strecke 1 und Strecke 2
!              gemessen am Punkt x,y , beginnend an der Strecke 1
!              Wertebereich: -pi bis pi
!
!  local variables
!
      real (DP) x0,y0,xi1,yi1,xi2,yi2
!
      xi1=x1-x
      yi1=y1-y
      xi2=x2-x
      yi2=y2-y
      x0=xi1*xi2+yi1*yi2
      y0=-yi1*xi2+xi1*yi2
!  Koordinaten des Vektors (xi2,yi2) im xi1,yi1 Koordinatensystem
      wink1=atan2(y0,x0)
      return
      end function wink1
!
!
!
!
      real (DP) function winkl(la,mi,ne,xn,yn)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) la, mi, ne
      real (DP) xn(:), yn(:)
      intent (in) :: la, mi, ne, xn, yn
!  Eingabe: 
!    la,mi,ne     Knotennummern fuer die der Winkel zu berechnen ist
!    xn,yn        Vektor der x,y Koordinaten der Knoten
!  Ausgabe: 
!    winkl        Winkel im Bogenmass zw. den Schenkeln (la,mi) und (mi,ne)
!                   Ausgabe des Winkels im Bereich ( -pi ... +pi )
!
!  local variables
      real (DP) dxl,dyl,dxn,dyn,zaehl,xnenn
!
!  la (=last), mi (=middle), ne (=next)
      dxl=xn(la)-xn(mi)
      dyl=yn(la)-yn(mi)
      dxn=xn(mi)-xn(ne)
      dyn=yn(mi)-yn(ne)
!  Der Zaehler ist zugleich die doppelte Dreiecksflaeche
      zaehl=dxl*dyn-dxn*dyl
      xnenn=-dxl*dxn-dyl*dyn
      if (zaehl.eq.0._DP) then
        if (xnenn.eq.0._DP) then
          print*,'*** Fehler in winkl, das Dreieck hat keine Flaeche'
          print*,'*** Knotennummer  :  und Koordinaten (x,y) '
          print*,la,' : ',xn(la),yn(la)
          print*,mi,' : ',xn(mi),yn(mi)
          print*,ne,' : ',xn(ne),yn(ne)
      pause
          winkl=0._DP
        else
          winkl=0._DP
        end if
      else
        winkl=atan2(zaehl,xnenn)
      end if
      return
      end function winkl
!
!
!
      subroutine angles(i,e,xn,yn,w1,w2,w3)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) i, e(:,:)
      real (DP) xn(:), yn(:), w1, w2, w3
      intent (in) :: i, e, xn, yn
      intent (out) w1, w2, w3
!  Input: 
!    i            element number
!    e            nodes of elements
!    xn,yn        node coordinates
!  Output: 
!    w1, w2, w3   inner angles of the triangle in the range ( -pi...+pi ) 
!                 wi is the angle at vertex i
!  local variables
      real (DP) dx1, dy1, dx2, dy2, dx3, dy3, zaehl, xnenn2, xnenn3
      real (DP) pi
      parameter (pi=3.14159265358979323846264338327950288419716_DP)
!
      dx3=xn(e(1,i))-xn(e(2,i))
      dy3=yn(e(1,i))-yn(e(2,i))
      dx1=xn(e(2,i))-xn(e(3,i))
      dy1=yn(e(2,i))-yn(e(3,i))
!  twice the area of the triangle
      zaehl=dx3*dy1-dx1*dy3
!
      xnenn2=-dx3*dx1-dy3*dy1
      w2=atan2(zaehl,xnenn2)
      dx2=xn(e(3,i))-xn(e(1,i))
      dy2=yn(e(3,i))-yn(e(1,i))
      xnenn3=-dx1*dx2-dy1*dy2
      w3=atan2(zaehl,xnenn3)
      if (w2 .ge. 0._DP) then
        w1=pi-w2-w3
      else 
        w1=-(5._DP*pi+w2+w3)
      end if
!      w1=atan2(zaehl,-dx2*dx3-dy2*dy3)
      return
      end subroutine angles