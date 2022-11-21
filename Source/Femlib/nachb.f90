      subroutine nachb(e,en,p,n)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) p, e(:,:), en(:,:), n
      intent (in) e, p, n
      intent (out) en
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
!  Datenaufbereitung des Netzes, Erzeugen der Nachbarschaft en(3,*) 
!
!    Eingabe:  
!      e      Elementinformation
!      p      Anzahl der Knoten
!      n      Anzahl der Elemente
!    Ausgabe:
!      en     Nachbarn der Elemente
!
!  Es werden auch offene Elemente beruecksichtigt
!  Der Algorithmus verwendet die Felder ke(p) und kne(3*n) 
!  zum abspeichern temporaerer Daten
!  Der Aufwand waechst linear (!) mit der Knotenanzahl
!
!  Die Ausgabe erfolgt sortiert im mathematisch positiven Sinn 
!
!  local variables
!
      integer (I4B) ke(0:p+1),kne(3*n)
      integer (I4B) i,j,kalt,kk,k,ie1,jj,ie2,l1,k1,k2,l2
      integer (I4B) k3,l,l3
      integer (I4B) inach(3),ivorg(3)
      parameter (inach=(/2,3,1/))
      parameter (ivorg=(/3,1,2/))
!
      kne(:)=0
      ke(:)=0
      en(:,1:n)=0
!
!  Anzahl der Elemente zu den Knoten feststellen
!  offene Elemente haben als dritten (Pseudo-)Knoten die Nummer 0
      do i=1,n
        do j=1,3
          ke(e(j,i))=ke(e(j,i))+1
        end do
      end do
!  innere Knoten haben soviel Nachbarknoten wie adjazente Elemente. Fuer Randknoten ist
!  die Anzahl der adjazenten Elemete um eins geringer als die der Nachbarknoten.
!  Hier wird als zusaetzliches Pseudoelement fuer die Randknoten "0" aufgennommen
!  damit wird spaeter der Nachbar von Randelemente auf "0" gesetzt
!
!  fuer kompaktes Abspeichern vorbereiten
      kalt=ke(1)
      ke(1)=1
      do i=2,p+1
        kk=ke(i)
        ke(i)=ke(i-1)+kalt
        kalt=kk
      end do
!  Elemente zu den Knoten feststellen und kompakt abspeichern
!  in zwei Haeften um den Aufwand fuer das suchen einer freien Position klein zu halten
!  erste Haelfte
      do i=1,n/2
        do j=1,3
          if (e(j,i).eq.0) cycle
!  vorwaerts suchen
          do k=ke(e(j,i)),ke(e(j,i)+1)-1
            if (kne(k).eq.0) then
              kne(k)=i
              exit
            end if
          end do
        end do
      end do
!  zweite Haelfte
      do i=n/2+1,n
        do j=1,3
          if (e(j,i).eq.0) cycle
!  rueckwaerts suchen
          do k=ke(e(j,i)+1)-1,ke(e(j,i)),-1
            if (kne(k).eq.0) then
              kne(k)=i
              exit
            end if
          end do
        end do
      end do
!
!  fuer alle Knoten
!  Nachbarschaft der Elemente feststellen
      do i=1,p
        do j=ke(i),ke(i+1)-1
!  ie1 erstes Element
          ie1=kne(j)
!  suche lokale Knotennummer mit e(k1,ie1) = i
          do jj=1,3
            if (e(jj,ie1).eq.i) then
              k1=jj
              exit
            end if
          end do
          k2=inach(k1)
!
          do k=j+1,ke(i+1)-1
!  ie2 zweites Element
            ie2=kne(k)
!  pruefe ie1 und ie2 auf Nachbarschaft
!  suche lokale Knotennummer mit e(l1,ie2) = i
            do l=1,3
              if (e(l,ie2).eq.i) then
                l1=l
                exit
              end if
            end do
!
!  pruefe zweiten geneinsamen Knoten
            l2=ivorg(l1)
            if (e(k2,ie1).eq.e(l2,ie2)) then
!  Nachbarschaft festgestellt
              l3=inach(l1)
              k3=ivorg(k1)
              en(l3,ie2)=ie1
              en(k3,ie1)=ie2
              exit
            end if
          end do
        end do
      end do
      return
      end subroutine nachb
