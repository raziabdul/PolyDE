      subroutine mplist(gkz,bz,bzi,bzip,bzil,zki,zpz,                   &
     &  plist,pz,plp,zahl)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) gkz,bz,bzi(:),bzip(:),bzil(:),zki(:,:)
      integer (I4B) zpz(:),plist(:),pz,plp(:),zahl
      intent (in) :: gkz, bz, bzi, bzip, bzil, zki, zpz
      intent (out) :: plist, plp, pz, zahl
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
!  Input:
!            gkz      Gebiets-Knoten-Zahl
!            bz       Number of the domain
!            bzi      Bereichs-Zweige-Informationen (Liste der Zweige eines Bereichs
!                     bzw, Gebietes)
!            bzip     Pointer auf Anfang eines Gebiets in Bzi-liste 
!                     (kompakte Speicherung)
!            bzil     Pointer auf bzip Anfang eines Gebietes im Vektor einschliesslich
!                     aller darin enhaltenen Gebiete
!            zki      Zweige-Knoten-Informationen (Liste der Knoten zu jedem Zweig)
!                     zki(1,i)  Anfangspunkt des Zweiges
!                     zki(2,i)  Endpunkt des Zweiges
!                     zki(3,i)  Dritter Punkt des Zweiges
!                               /  = 0  :  Gerade
!                     zki(3,i)  -  > 0  :  Punkt auf dem Kreisbogen
!                               \  < 0  :  Mittelpunkt des Kreisbogens
!            zpz      Zahl der Punkte auf dem Zweig (incl. Anfangs- und Endpunkt)
!  Output:
!            plist    List of Nodes of the actual region
!            plp      Point-List Pointer (compact storage)
!            pz       Number of points in the list
!            zahl     Anzahl der Zweiglisten der Berandung
!                     (bei einfach zusammenhaengenden Gebieten = 1)
!
!  Mplist - Mache Punkt Liste
!           erstellt eine sortierte Punktliste des Randes des Gebiets
!           mit der Nummer Bz.
!           Der Rand wird mathematisch positiv umlaufen.
!
!  local variablea
      integer (I4B) i,j,k
      integer (I4B) lauf,rand,zweig
! 
      zahl=bzil(bz+1)-bzil(bz)
!
      pz=1
!  Fuer alle Raender des Gebiets
      do lauf=1,zahl
        plp(lauf)=pz
        rand=bzil(bz)+lauf-1
!  Fuer alle Zweige des Rands
        do i=1,bzip(rand+1)-bzip(rand)
!  Bestimme Nummern der Knoten auf dem Zweig i:
          k=gkz
          zweig=bzi(bzip(rand)+i-1)
          do j=1,abs(zweig)-1
            k=k+zpz(j)-2
          end do
          if (zweig.gt.0) then
!  Wenn Zweignummer positiv
            plist(pz)=zki(1,zweig)
            do j=1,zpz(zweig)-2
              pz=pz+1
              plist(pz)=k+j
            end do
          else
!  Wenn Zweignummer negativ
            plist(pz)=zki(2,-zweig)
            do j=1,zpz(-zweig)-2
              pz=pz+1
              plist(pz)=k+zpz(-zweig)-2-j+1
            end do
          end if
          pz=pz+1
        end do
      end do
      plp(zahl+1)=pz
      pz=pz-1
      return
      end subroutine mplist
