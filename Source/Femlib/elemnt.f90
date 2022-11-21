! 
      subroutine elemnt(xn,yn,e,n,xt,yt,ielem,en,ok)
      use feminterface, only: elemnd
      use femtypes
      implicit none
      real (DP) xn(:), yn(:), xt, yt
      integer (I4B) e(:,:), n, ielem, en(:,:)
      logical ok
      intent (in) :: xn, yn, e, n, xt, yt, en
      intent (out) :: ielem, ok
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
!    $Date: 2010/12/08 08:52:28 $
!    $Author: chvokas $
!
!  Ueberpruefen in welchem Element ielem der Punkt
!  mit den Koordinaten (xt/yt) liegt
!  Eingabe:  xn,yn  Koordinaten der Knoten
!            e      Elementinformation
!            n      Anzahl der Elemente
!            xt,yt  Koordinaten des zu ueberpruefenden Punktes
!            en     Nachbarn der Elemente
!  Ausgabe:  ielem  Elementnummer in dem der Knoten liegt
!            ok     =false wenn kein Element gefunden wurde
!
!  Es wird zunaechst versucht, ausgehend von zuletzt gefundenen Element,
!  jeweils die Nachbarn in der Richtung des Punktes abzusuchen.
!  Dies fuehrt leider nicht immer zum Ziel, wenn das Gebiet konkav ist
!  bzw. wenn das Punkt ausserhalb des triangulierten Bereichs liegt
!  Wenn das Element mit der Methode nicht gefunden wird, werden alle Elemente abgesucht
      real (DP) x1,y1,x2,y2
      integer (I4B) ialt,i,jend,kend,ineu,nachf(3)
      logical found
      integer (I4B), save :: ilea=1
      parameter (nachf=(/2,3,1/))
!
      ialt=0
!  Element in dem der Punkt beim letzten Aufruf gefunden wurde
      i=ilea
!  Zaehler fuer die besuchten Elemente
!
10    continue
      found=.true.
      x2=xn(e(3,i))-xt
      y2=yn(e(3,i))-yt
      do jend=1,3
!  Anfangsknoten
        x1=x2
        y1=y2
!  Endknoten
        kend=e(jend,i)
        x2=xn(kend)-xt
        y2=yn(kend)-yt
!  Liegt der Punkt rechts von Strahl
!              durch die Punkte xanf/yanf und xend/yend
        if((x1*y2-x2*y1).lt.0._DP) then
!  Element gegenueber dem Strahl janf, jend
          ineu=en(nachf(jend),i)
          if ( (ineu.le.0) .or. (ineu.eq.ialt) ) then
            found=.false.
          else
            ialt=i
            i=ineu
            goto 10
          end if
        end if
      end do
!
      if (found) then
        ielem=i
        ok=.true.
      else
!  Nicht gefunden jetzt hilft nur noch vollstaendiges Absuchen aller Elemente
        call elemnd(xn,yn,e,n,xt,yt,ielem,ok)
        if (.not. ok) ielem=i
      end if
      ilea=ielem
      return
      end subroutine elemnt
!
!
!
      subroutine elemnd(xn,yn,e,n,xt,yt,ielem,ok)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) xn(:), yn(:), xt, yt
      integer (I4B) e(:,:), n, ielem
      logical ok
      intent (in) :: xn, yn, e, n, xt, yt
      intent (out) :: ielem, ok
!  Ueberpruefen in welchem Element ielem der Punkt
!  mit den Koordinaten (xt/yt) liegt
!  Eingabe:  xn,yn  Koordinaten der Knoten
!            e      Elementinformation
!            n      Anzahl der Elemente
!            xt,yt  Koordinaten des zu ueberpruefenden Punktes
!  Ausgabe:  ielem  Elementnummer in dem der Knoten liegt
!            ok     =false wenn kein Element gefunden wurde
      integer (I4B) i
      real (DP) xl(3),yl(3)
!
      ok=.false.
      do i=1,n
        xl(:)=xn(e(:,i))
        if(xt.gt.maxval(xl)) cycle
        if(xt.lt.minval(xl)) cycle
!  Liegt der Punkt im Wertebereich der x-Werte des Elementes?
        yl(:)=yn(e(:,i))
        if(yt.gt.maxval(yl)) cycle
        if(yt.lt.minval(yl)) cycle
!  Liegt der Punkt im Wertebereich der y-Werte des Elementes?
        xl=xl-xt
        yl=yl-yt
!  Liegt der Punkt links vom Strahl durch die Punkte x1/y1 und x2/y2?
!  wenn ja dann ok
        if((xl(1)*yl(2)-xl(2)*yl(1)).lt.0.d0) cycle
!  Liegt der Punkt links vom Strahl durch die Punkte x2/y2 und x3/y3?
        if((xl(2)*yl(3)-xl(3)*yl(2)).lt.0.d0) cycle
!  Liegt der Punkt links vom Strahl durch die Punkte x3/y3 und x1/y1?
        if((xl(3)*yl(1)-xl(1)*yl(3)).lt.0.d0) cycle
        ielem=i
        ok=.true.
        exit
      end do
      return
      end subroutine elemnd
