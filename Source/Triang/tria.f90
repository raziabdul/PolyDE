      subroutine tria(xk,yk,plist,plp,zahl,prlist,pzr,prz,ele,ez,       &
     &  ezalt,kdim,geb,nb,bz,erfolg)
      use feminterface, only: genele, innen
      use femtypes
      implicit none
      integer (I4B) plist(:),plp(:),zahl,prlist(:,:),pzr,prz
      integer (I4B) ele(:,:),ez,ezalt,kdim,geb(:),nb(:,:),bz
      real (DP) xk(:),yk(:)
      logical erfolg
      intent (in) :: xk, yk, zahl, pzr, ezalt, kdim, bz
      intent (out) :: erfolg
      intent (inout) :: plist, plp, prlist, prz, ele, ez, geb, nb
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
!  Input
!            xk,yk    Co-ordinated of the Nodes
!            zahl     Number of branch lists
!                     (equal to 1 for simply connected regions)
!            pzr      Number of Nodes of this Domain
!            ezalt    former Number of Elements
!            kdim     Declared Size for the Number of Elements
!            bz       Number of the actual Domain
!  Output:
!            erfolg   =.true. if successful
!  In-/ Output:
!            plist    List of Nodes of the actual Domain
!            plp      Point-List Pointer (compact Storage)
!            prlist   Point-Edge-List
!                     contains the connections of the edge-nodes which are located 
!                     on the inner edges and thus have to be used twice with opposite
!                     direction (used for the generation of neighbors)
!              prlist(1,*) - Starting node of the section
!              prlist(2,*) - End node of the section
!              prlist(3,*) - Number of the first neighboring element
!              prlist(4,*) - Number of the opposite node in the Element
!            prz      number of nodes in prlist
!            ele      Nodes of the Elements
!            ez       Number of Elements
!            geb      Domains of the Elements
!            nb       Neighbors of the Elements
!
!  Tria - Aufstellen der Triangulierung
!     Konzept:
!        Die einzelnen Gebiete werden nach der Generierung der Knoten
!        auf den Zweigen durch Randpolygone beschrieben. Da sich alle
!        Polygone in Dreiecke zerlegen lassen, ist eine so erzeugte
!        Triangulierung die groebste Moeglichkeit ein Gebiet zu
!        zerlegen. In Anbetracht der adaptiven Verfeinerung erscheint
!        es sinnvoll, so grob wie moeglich zu triangulieren.
!
!  lokale Variablen
      integer (I4B) pkt,pkt1,pkt2,pkt3,zahl1,zahl2
      integer (I4B) start,iplus0,iplus1,iplus2,ip0,ip1,ip2,lauf,basis,i
      logical zerleg
!
!  Merken der alten Daten:
      do 53124 lauf=1,zahl
        basis=plp(lauf)-1
        zahl2=pzr
        zahl1=plp(lauf+1)-plp(lauf)
        pkt=1
999     continue
        start=pkt
!      if (zahl1.gt.2) then
!  Es ist noch mehr als ein Element uebrig:
        i=start
!
543     continue
        ip0=i+basis
        ip1=mod(i,zahl1)+1+basis
        ip2=mod(i+1,zahl1)+1+basis
        iplus0=plist(ip0)
        iplus1=plist(ip1)
        iplus2=plist(ip2)
        zerleg=.false.
!
!  Wenn aus drei aufeinanderfolgenden Punkten auf dem
!  Rand ein Dreieck mit positiver Flaeche gebildet werden koennte,
        if (((xk(iplus2)-xk(iplus0))*(yk(iplus1)-yk(iplus0))-           &
     &    (yk(iplus2)-yk(iplus0))*(xk(iplus1)-xk(iplus0))).lt.0) then
!  Dann pruefe, ob keine anderen Randpunkte innen liegen.
          if (.not.innen(xk,yk,plist,zahl2,ip0,ip1,ip2)) zerleg=.true.
!  Ein derartiges Dreieck darf abgetrennt werden!
!
        end if
        if (.not.zerleg) then
          i=mod(i,zahl1)+1
          if (i.eq.start) then
            if (zahl.eq.1) then
              write (*,1001)
1001          format(' ***** Fehler, Zerlegung unmoeglich')
!        goto 53124
              stop 15
            else
              goto 53124
            end if
          end if
          goto 543
        end if
!      else
!  Wenn nur noch 3 Knoten uebrig, kann nur noch ein Element verbleiben
!       i=1
!      end if
!  An dieser Stelle zeigt i auf das naechste abzuschneidende Element
!
!  Beginne die Elementgenerierung
!
        pkt=i
        pkt1=plist(pkt+basis)
        pkt2=plist(mod(pkt,zahl1)+1+basis)
        pkt3=plist(mod(pkt+1,zahl1)+1+basis)
        call genele(ez,ezalt,ele,nb,geb,kdim,bz,pkt1,pkt2,pkt3,         &
     &   prlist,prz)
        erfolg=.true.
!  Rand um generiertes Dreieck verkleinern:
        zahl1=zahl1-1
        zahl2=zahl2-1
        do i=mod(pkt,zahl1+1)+1+basis,zahl2
          plist(i)=plist(i+1)
        end do
        do i=lauf+1,zahl+1
          plp(i)=plp(i)-1
        end do
        if (pkt.eq.zahl1+1) pkt=pkt-1
        pkt=pkt-1
        if (pkt.lt.1) pkt=pkt+zahl1
!  Noch genug Punkte fuer ein neues Dreieck uebrig? Dann noch einmal!
        if (zahl1.ge.3) goto 999
53124 continue
      return
      end subroutine tria
