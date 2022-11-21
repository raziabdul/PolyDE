      subroutine inkrp(x1,y1,x2,y2,x3,y3,xpi,ypi,flaech2,rho,           &
     &  xpa,ypa,rausen)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) x1,y1,x2,y2,x3,y3,xpi,ypi,flaech2,rho,xpa,ypa,rausen
      intent (in) :: x1,y1,x2,y2,x3,y3
      intent (out) :: xpi,ypi,flaech2,rho,xpa,ypa,rausen
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
!  lokale variablen
      real(DP) dx1,dy1,d1,dx2,dy2,d2,dx3,dy3,d3,umfang,beta
!      real(DP) al2,al3
!  Berechnung von Inkreises und Umkreis und Flaeche
!  des Dreieckes mit den Eckpunkten 1,2,3 
!
!     xpi,ypi    Koordinaten des Inkreismittelpunkt
!     flaech2    doppelte Dreiecksflaeche
!     rho        Radius des Inkreises
!     xpa,ypa    Koordinaten des Umkreimittelpunktes
!     rausen     Umkreisradius
!
!  Der Radius wird negativ ausgeben, wenn die flaeche negativ ist
!
!  Berechnung der Abstaende und Seitenlaengen
      dx1=y2-y3
      dy1=x2-x3
      d1=sqrt(dx1*dx1+dy1*dy1)
      dx2=y3-y1
      dy2=x3-x1
      d2=sqrt(dx2*dx2+dy2*dy2)
      dx3=y1-y2
      dy3=x1-x2
      d3=sqrt(dx3*dx3+dy3*dy3)
      umfang=d1+d2+d3
!  Dreieckskoordinaten des Inkreismittelpunkts
!      al2=d2/umfang
!      al3=d3/umfang
!  "weltkoordinaten" des Innkreimittelpunktes
!      xpi=x1-al2*dy3+al3*dy2
!      ypi=y1-al2*dx3+al3*dx2
      xpi=x1-(d2*dy3+d3*dy2)/umfang
      ypi=y1-(d2*dx3+d3*dx2)/umfang
!  doppelete Dreiecksflaeche
      flaech2=dx3*dy2-dx2*dy3
!  Berechnung des Inkreisradius
      rho=flaech2/umfang
!  Rerechnung des Umkreisradius
      if (flaech2.eq.0._DP) then
        rausen=huge(1._DP)
        xpa=xpi
        ypa=ypi
      else        
        rausen=d1*d2*d3/(2._DP*flaech2)
!  Koordinaten des Umkreismittelpunktes
        beta=(dx3*dx2+dy3*dy2)/flaech2
        xpa=((x2+x3)-beta*dx1)/2._DP
        ypa=((y2+y3)+beta*dy1)/2._DP
      end if
      return
      end
