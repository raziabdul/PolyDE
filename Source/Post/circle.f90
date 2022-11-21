      subroutine circle(x, y, w1, w2, r1, r2)
      use femtypes
      implicit none
      real(SP) x, y, w1, w2, r1, r2
      intent (in)::  x, y, w1, w2, r1, r2
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
!    $Revision: 1.7 $
!    $Date: 2014/02/18 10:45:08 $
!    $Author: juryanatzki $
!
!  Zeichnen eines Kreisbogens oder einer Spirale
!    x     x-Koordinate des Anfangspunktes
!    y     y-Koordinate des Anfangspunktes
!    w1    Anfangswinkel in Grad
!    w2    Endwinkel in Grad
!    r1    Anfangsradius
!    r2    Endradius
      real (SP) f, ff, da, a1, a2, sda, cda
      real (SP) dr, xx, yy, xm, ym, tt, pgs
      integer (I4B) n,i, npgs
!  pi/180
      parameter (f=0.01745329251994329576923691_SP)
!  the number of segments on 180 degrees
      parameter (npgs=180)
!  Fehlerabfragen
      if(r1 .lt. 0._SP .or.  r2 .lt. 0._SP) return
      if(r1 .eq. 0._SP .and. r2 .eq. 0._SP) return
!  Anzahl der Polygon.Schritte festlegen
      a1 = w1*f
      a2 = w2*f
      if (a1 .eq. a2) a2 = a2 + 6.283185307179586476925287_SP
!  Polygonstuecke auf 180 Grad
      pgs=real(npgs,SP)/180._SP
      n = abs(a2-a1) / f * pgs + 1._SP
!  Anfangspunkt
      call pgmove(x, y)
!  Schleifenkonstante bestimmen
      da = (a2-a1)/n
      sda = sin(da)
      cda = cos(da)
      dr = (r2-r1)/n
      ff = 1._SP + dr/r1
      xx = r1 * cos(a1)
      yy = r1 * sin(a1)
      xm = x - xx
      ym = y - yy
!  Zeichnen der Kurve
      do i=1,n
        tt = ff * (xx*cda - yy*sda)
        yy = ff * (xx*sda + yy*cda)
        xx = tt
        call pgdraw(xx + xm, yy + ym)
        ff = 2._SP - 1._SP/ff
      end do
      return
      end subroutine circle
!
!
!
      subroutine ellipse(xc, yc, w1, w2, a, b, phi0)
      use femtypes
      implicit none
      real(SP) :: xc, yc, w1, w2, a, b, phi0
      intent (in)::  xc, yc, w1, w2, a, b, phi0
!
!  draw an ellipse
!    x0    x-co-ordinate of the centre
!    y0    y-co-ordinate of the centre
!    w1    start angle in degrees
!    w2    end angle in degrees
!    a     semimajor axis
!    b     semiminor axis
!    phi0  angle between the major axis and the x-axis
      real (SP) :: f, a1, a2, xx, yy, phi0r, phi, pgs, da
      real (SP) :: as0, bs0, ac0, bc0
      integer (I4B) :: n,i, npgs
!  pi/180
      parameter (f=0.01745329251994329576923691_SP)
!  the number of segments on 180 degrees
      parameter (npgs=180)
!  error checking
      if(a .lt. 0._SP .or.  b .lt. 0._SP) return
!      tiny=0
      if(a .eq. 0._SP .and. b .eq. 0._SP) return
!  determine number of polygon segments 
      a1 = w1*f
      a2 = w2*f
      if (a1 .eq. a2) a2 = a2 + 6.283185307179586476925287_SP
!  number of polygon segments on an angle of 180 degrees
      pgs = real(npgs,SP)/180._SP
      n = abs(a2-a1) / f * pgs + 1._SP
!  start point
      phi0r = phi0*f
      xx = xc + ( (a+b)*cos(a1+phi0r) + (a-b)*cos(a1-phi0r) ) /2._DP
      yy = yc + ( (a+b)*sin(a1+phi0r) - (a-b)*sin(a1-phi0r) ) /2._DP
      call pgmove(xx, yy)
!  prepare loop constants
      as0 = a * sin(phi0r)
      ac0 = a * cos(phi0r)
      bs0 = b * sin(phi0r)
      bc0 = b * cos(phi0r)
      da = (a2-a1)/n
      phi = phi0r
!  draw the ellipse
      do i=1,n
        phi= phi + da
        xx = xc + ac0 * cos(phi) - bs0 * sin(phi)
        yy = yc + as0 * cos(phi) + bc0 * sin(phi)
        call pgdraw(xx, yy)
      end do
      return
      end subroutine ellipse