      logical function innen(xk,yk,plist,pz,num1,num2,num3)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) xk(:),yk(:)
      integer (I4B) plist(:),pz
      integer (I4B) num1,num2,num3
      intent (in) :: xk,yk,plist,pz,num1,num2,num3
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
!    $Date: 2006/07/08 00:18:59 $
!    $Author: r_abdul $
!
!  Eingabe:
!            xk,yk    Co-ordinates of the (Mesh-) Nodes
!            plist    List of Node-Numbers
!            pz       Number of Nodes in the List
!            num1     First Node of the Triangle
!            num2     Second Node of the Triangle
!            num3     Third Node of the Triangle
!
! Innen prueft, ob einer der Punkte aus einer Liste innerhalb eines Dreiecks
! liegt
! Eingaben:
! Plist,pz (integer) - Punktliste und -zahl der Punkte
! mumi (integer) - Nummern der Eckknoten des Dreiecks
! Ausgabe:
! Innen: .true. wenn Punkt innerhalb oder auf dem Rand liegt.
!        .false. wenn Punkt ausserhalb
!
!  local variables
      integer (I4B) i,lauf
      real (DP) xi1,yi1,xi2,yi2,xi3,yi3,x1,y1,x2,y2,x3,y3
!
      x1=xk(plist(num1))
      y1=yk(plist(num1))
      x2=xk(plist(num2))
      y2=yk(plist(num2))
      x3=xk(plist(num3))
      y3=yk(plist(num3))
! Bestimme drei Richtungsvektoren der Dreieckskanten:
      xi1=x2-x1
      yi1=y2-y1
      xi2=x3-x2
      yi2=y3-y2
      xi3=x1-x3
      yi3=y1-y3
! Test fuer alle anderen Punkte, ob sie innerhalb
! des Dreiecks liegen:
      innen=.false.
      do lauf=1,pz
        if ((lauf.eq.num1).or.(lauf.eq.num2).or.(lauf.eq.num3)) then
          cycle
        end if
        i=plist(lauf)
        if (((xk(i)-x3)*yi3-(yk(i)-y3)*xi3).le.0._DP) then
          if (((xk(i)-x1)*yi1-(yk(i)-y1)*xi1).le.0._DP) then
            if (((xk(i)-x2)*yi2-(yk(i)-y2)*xi2).le.0._DP) then
              innen=.true.
              exit
            end if
          end if
        end if
      end do
      return
      end function innen
