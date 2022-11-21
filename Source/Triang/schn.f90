      logical function schn(x1,y1,x2,y2,x3,y3,x4,y4)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) x1,y1,x2,y2,x3,y3,x4,y4
      intent (in) :: x1,y1,x2,y2,x3,y3,x4,y4
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
!    $Date: 2006/07/08 00:18:59 $
!    $Author: r_abdul $
!
!  Unterprogramm 'berechne Schnitt':
!  Hier wird die Lage eines eventuellen Schnittpunktes numerisch bestimmt
!  anschliessend wird geprueft, ob der Schnittpunkt ueberhaupt im Bereich
!  der untersuchten Strecken angeordnet ist. Ob dies der Fall ist, wird
!  mit Hilfe der Variablen 'schnitt' an das aufrufende Programm gemeldet.
!
!  Eingabe:
!            x1,y1    Koordinaten des ersten Punktes der ersten Strecke
!            x2,y2    Koordinaten des zweiten Punktes der ersten Strecke
!            x3,y3    Koordinaten des ersten Punktes der zweiten Strecke
!            x4,y4    Koordinaten des zweiten Punktes der zweiten Strecke
!  Ausgabe:
!            schn     =.true. wenn der Schnitt der Geraden auf beiden Strecken liegt
  
!  lokale variablen
      real (DP) deter,mue,nue
      deter = -(x3-x4)*(y2-y1)+(y3-y4)*(x2-x1)
      if (abs(deter).eq.0._DP) then
!  Die Determinante ist zu klein um damit einen Schnittpunkt mit aus-
!  reichender Genauigkeit zu bestimmen. Das ist auch nicht zu tragisch,
!  denn bei einer derartig kleinen Determinante verlaufen die Linien
!  sowieso ganz oder fast parallel. Also: Don't Panic!
        schn = .false.
!  Natuerlich auch kein Schnitt!
      else
!  der Schnittpunkt kann nun berechnet werden!
        mue = ((y3-y4)*(x3-x1)-(x3-x4)*(y3-y1))/deter
        nue = ((y1-y3)*(x1-x2)-(x1-x3)*(y1-y2))/deter
!  nun wird getestet, ob der Schnittpunkt im Bereich der Linien liegt:
        if ( (mue.gt.1._DP) .or. (mue.lt.0._DP) .or.                    &
      &      (nue.gt.1._DP) .or. (nue.lt.0._DP) ) then
          schn=.false.
        else
          schn=.true.
        end if
      end if
      return
      end function schn
