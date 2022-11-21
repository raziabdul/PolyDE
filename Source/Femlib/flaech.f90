!
      pure real (DP) function flaech(x1,y1,x2,y2,x3,y3)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) x1,y1,x2,y2,x3,y3
      intent (in) :: x1,y1,x2,y2,x3,y3
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
!    $Date: 2010/10/25 11:56:58 $
!    $Author: chvokas $
!
!  Flaech - Berechnung der Flaeche eines Dreiecks
!  Eingaben:
!  xi,yi (real) - Die drei Punkte des Dreiecks
!  Ausgabe:
!  flaech (real) - Flaeche
!
      flaech=((y2-y3)*(x1-x2)-(y1-y2)*(x2-x3))/2._DP
      return
      end function flaech
!
