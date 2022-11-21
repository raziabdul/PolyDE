      subroutine genele(ez,ezalt,ele,nb,geb,kdim,bz,                    &
     &  pkt1,pkt2,pkt3,prlist,prz)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) ez,ezalt,ele(:,:),nb(:,:),geb(:),kdim,bz
      integer (I4B) pkt1,pkt2,pkt3,prlist(:,:),prz
      intent(in) :: pkt1, pkt2, pkt3, ezalt, kdim, bz
      intent(out) :: nb, geb
      intent(inout) :: ez, ele, prlist, prz
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
!  Eingabe
!            pkt1     Knotennummer 1 des Elementes
!            pkt2     Knotennummer 2 des Elementes
!            pkt3     Knotennummer 3 des Elementes
!            ezalt    Anzahl der Elemente
!            kdim     Vereinbarte Groesse fuer die Anzahl der Knoten
!            bz       Bereichszahl, Nummer des Gebietes
!  Ausgabe:
!            nb       Nachbarn der Elemente
!            geb      Gebiete der Elemente
!  Ein-/ Ausgabe
!            ez       Anzahl der Elemente
!            ele      Knoten der Elemente
!            prlist   Punkt-Rand-Liste
!                     enthaelt die Verbindungen der Zweigknoten,
!                     die auf inneren Zweigen liegen und somit
!                     einmal in jede Richtung (+,-) durchlaufen
!                     werden muessen (zweck: Generierung der Nachbarn)
!             prlist(1,*) - Anfangsknoten des Abschnitts
!             prlist(2,*) - Endknoten des Abschnitts
!             prlist(3,*) - Nummer des ersten angrenzenden Elements
!             prlist(4,*) - elementbezogene Nummer des gegenueber-
!                           liegenden Knotens
!            prz      
!
!  lokale variablen
      integer (I4B) i,j,merk1,merk2,merk3
!
!  Unterprogramm zur Generierung eines Elements im Gebiet bz
!
      if (ez+1.gt.2*kdim) then
        write (*,1002)
1002    format(' Fehler: zu viele Elemente')
        stop 15
      end if
!  Element zufuegen:
      ez=ez+1
      ele(1,ez)=pkt1
      ele(2,ez)=pkt2
      ele(3,ez)=pkt3
      geb(ez)=bz
      do 1010 i=1,3
        nb(i,ez)=0
1010  continue
!      write (*,*) pkt1,pkt2,pkt3,bz
!  Jetzt suche alle Nachbarn innerhalb des selben Gebiets:
      do 1011 i=ezalt+1,ez-1
        if ((pkt1.eq.ele(2,i)).and.(pkt2.eq.ele(1,i))) then
          nb(3,ez)=i
          nb(3,i)=ez
        end if
        if ((pkt1.eq.ele(3,i)).and.(pkt2.eq.ele(2,i))) then
          nb(3,ez)=i
          nb(1,i)=ez
        end if
        if ((pkt1.eq.ele(1,i)).and.(pkt2.eq.ele(3,i))) then
          nb(3,ez)=i
          nb(2,i)=ez
        end if
        if ((pkt2.eq.ele(2,i)).and.(pkt3.eq.ele(1,i))) then
          nb(1,ez)=i
          nb(3,i)=ez
        end if
        if ((pkt2.eq.ele(3,i)).and.(pkt3.eq.ele(2,i))) then
          nb(1,ez)=i
          nb(1,i)=ez
        end if
        if ((pkt2.eq.ele(1,i)).and.(pkt3.eq.ele(3,i))) then
          nb(1,ez)=i
          nb(2,i)=ez
        end if
        if ((pkt3.eq.ele(2,i)).and.(pkt1.eq.ele(1,i))) then
          nb(2,ez)=i
          nb(3,i)=ez
        end if
        if ((pkt3.eq.ele(3,i)).and.(pkt1.eq.ele(2,i))) then
          nb(2,ez)=i
          nb(1,i)=ez
        end if
        if ((pkt3.eq.ele(1,i)).and.(pkt1.eq.ele(3,i))) then
          nb(2,ez)=i
          nb(2,i)=ez
        end if
1011  continue
!  Jetzt suche auch noch die Nachbarn auf den Randteilen, die aus inneren
!  Zweigen gebildet werden.
      merk1=0
      merk2=0
      merk3=0
      do 3000 i=1,prz
        if (((prlist(2,i).eq.pkt1).and.(prlist(1,i).eq.pkt2)).or.       &
     &    ((prlist(2,i).eq.pkt2).and.(prlist(1,i).eq.pkt1))) then
          if (prlist(3,i).eq.0) then
!  Fuer dieses Randstueck ist noch kein Nachbarelement bekannt
!  also Elementnummer vermerken
            prlist(3,i)=ez
            prlist(4,i)=3
          else
!  Nachbarelement schon vermerkt. Also wird fuer beide Elemente die
!  Nachbarschaftsinformation ergaenzt.
            merk1=i
            nb(3,ez)=prlist(3,i)
            nb(prlist(4,i),prlist(3,i))=ez
          end if
        else
          if (((prlist(2,i).eq.pkt3).and.(prlist(1,i).eq.pkt2)).or.     &
     &      ((prlist(2,i).eq.pkt2).and.(prlist(1,i).eq.pkt3))) then
            if (prlist(3,i).eq.0) then
              prlist(3,i)=ez
              prlist(4,i)=1
            else
              merk2=i
              nb(1,ez)=prlist(3,i)
              nb(prlist(4,i),prlist(3,i))=ez
            end if
          else
            if (((prlist(2,i).eq.pkt3).and.(prlist(1,i).eq.pkt1)).or.   &
     &        ((prlist(2,i).eq.pkt1).and.(prlist(1,i).eq.pkt3))) then
              if (prlist(3,i).eq.0) then
                prlist(3,i)=ez
                prlist(4,i)=2
              else
                merk3=i
                nb(2,ez)=prlist(3,i)
                nb(prlist(4,i),prlist(3,i))=ez
              end if
            end if
          end if
        end if
3000  continue
      if (merk1.gt.0) then
!  Ein inneres Zweigstueck wird von zwei Elementen geteilt,
!  ist also verbraucht
        prz=prz-1
        do 3001 i=merk1,prz
          do 3002 j=1,4
            prlist(j,i)=prlist(j,i+1)
3002      continue
3001    continue
        if (merk2.gt.merk1) merk2=merk2-1
        if (merk3.gt.merk1) merk3=merk3-1
      end if
      if (merk2.gt.0) then
        prz=prz-1
        do 3003 i=merk2,prz
          do 3004 j=1,4
            prlist(j,i)=prlist(j,i+1)
3004      continue
3003    continue
        if (merk3.gt.merk2) merk3=merk3-1
      end if
      if (merk3.gt.0) then
        prz=prz-1
        do 3005 i=merk3,prz
          do 3006 j=1,4
            prlist(j,i)=prlist(j,i+1)
3006      continue
3005    continue
      end if
      return
      end subroutine genele
