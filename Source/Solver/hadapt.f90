      subroutine hadapt(ext,int,ende,error,ende2,epsgl)
      use feminterface, only: wfein, ffein, zeit, reallocate, glaett
      use feminterface, only: newelements, getsetting, preassemb 
      use feminterface, only: assemblyandsolve
      use globalvariables, only: ep, nnat, n, e, en, geb, kzi, p, x, &
                               & xbk, xn, yn, zki, zrb, ybk, ndof, physics
      use femtypes
      implicit none
      real (DP) :: error, epsgl
      logical :: ende, ende2, ext, int
      intent(in) :: ext, int
      intent(out) :: ende, error, ende2, epsgl
!
!-------------------------------------------------------------------------------
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
!    $Revision: 1.22 $
!    $Date: 2015/04/01 10:57:06 $
!    $Author: juryanatzki $
!
!-------------------------------------------------------------------------------
!
!  Adaptive mesh-refinement:
!  First of all an error estimation is done. Elements with higher residual than
!  adaption_error divided by number of elements are splitted into four subele- 
!  ments. Neighboring elements are splitted into two elements to restore confor-
!  mity if needed.
!
!  Input:
!            matzif   material index with which the region is filled
!            kdim     declared dimension for number of nodes
!            omega    angular frequency
!            zki      branch information, list of indices of key-points (geometry
!                     nodes) which determine the branches
!            kzrb     boundary condition of the key-points, index of the branch
!                     whose boundary condition is valid for the key-point
!            zrb      type of the boundary condition of the branch
!            xbk,ybk  coordinates of the key-points
!
!  Output:
!            ende     =.true. if MAXNODES is reached
!            error    relative error returned by a posteriori error estimation
!            ende2     =.true. if required tolerance is reached
!            epsgl    the error reached by the linear solver
!
!  In-/ Output: 
!            x        solution vector
!            e        node indices of the elements
!            en       neighbour-elements of the elements
!            n        total number of elements
!            xn,yn    coordinates of the nodes
!            p        total number of nodes
!            geb      assignment of elements to the geometry regions
!            kzi      node branch information
!
!  a maximum of n/5+1 elemente is refined.
!
!  local variables 
      integer (I4B) :: i, maxelemt, zahle, zahlk, fehlkn
      integer (I4B) :: adaptsteps, polyold, polyorder, maxnodes
      integer (I4B), allocatable :: ellist(:), liste5(:)
      real (DP) :: adapterror, resgl
      logical :: nregen
!
!  Allocate vector for polynomial degrees of elements.
      if (.not. associated(ep)) then
!  if there is no ep, then there is no solution --> get starting value for
!  polynomial degree from FEMsettings.txt
        allocate (ep(n,nnat))
        call getsetting('POLYORDER',polyorder)
        ep(1:n,1:nnat) = polyorder
        call getsetting('PHYSICS_MODE',physics)
        select case(physics)
          case('FLUIDINCOMPR','FLUIDELMASSCON')
            select case(polyorder)
              case(1,2)
                ep(1:n,1:2) = 2
                ep(1:n,3)   = 1
              case default
                ep(1:n,3)   = polyorder-1
            end select            
          case default
            !No special element order modification
        end select
      else if (associated(ep)) then
!  there is an ep and polyorder in FEMsettings < minimum value of the old vector
!  --> set initial value for polynomial order to minval(ep)
        call getsetting('POLYORDER',polyorder)
        polyold = minval(ep)
        if (polyorder .lt. polyold) then
          ep(1:n,1:nnat) = polyold
          print "(a)","Starting with minimum polynomial order from old solution."
          print "(a,i2,a)","Polynomial order homogenously set to ",polyold,"."
        else if (polyorder .ge. polyold) then
          ep(1:n,1:nnat) = polyorder
          print "(a)","Starting with polynomial order from FEMsettings.txt."
          print "(a,i2,a)","Polynomial order homogenously set to ",polyorder,"."
        end if
!  Special modification for the pressure field in fluidic problems
        call getsetting('PHYSICS_MODE',physics)
        select case(physics)
          case('FLUIDINCOMPR','FLUIDELMASSCON')
            ep(1:n,3)   = max(polyold,polyorder)-1
            print "(a,i2,a)","Polynomial order for fluidic pressure homogenously set to ",ep(1,3),"."
          case default
            !No special element order modification
        end select
      end if
      print *,'First Pre-assembly/solve...'
      call preassemb
!  Assemble and solve for the first time.
      call assemblyandsolve(epsgl,resgl)
!
      call getsetting('MAXNODES',maxnodes)
      if (p .ge. maxnodes) then 
        ende=.true.
        return
      end if
      call getsetting('ADAPTION_ERROR',adapterror)
      call getsetting('ADAPT_STEPS',adaptsteps)
!  START OF ADAPTSTEPS LOOP:
!
      print *, '============================================================'
      print *, '+ Adaptations do to: ', adaptsteps
!
      do i=1,adaptsteps
!
        print*, '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,2) i,n,p
2       format('+ h-adaptation step: ',i3,'   Elements=',i8,'   Nodes=',i8)
!
        maxelemt=n/5+1
!  refinement of elementen with large angles
        allocate(ellist(n))
        call wfein(ellist,zahle,zahlk,maxelemt/2,ende)
        if (.not.ende) then
!  refinement of elementen with large local error indicator
          call ffein(ellist,zahle,zahlk,ende,error,maxelemt,ext,int)
        end if
!
        print "(a,g11.4)", ' + Relative multi-nature a posteriori error: ', error
!
!  no more refinement if relative error is smaller than given error tolerance
        if (error .lt. adapterror)  then
          ende2= .true.
        else
          ende2= .false.
        end if
!
!  the actual refinement starts form here.
!  all elements to be refined are marked in ellist.
!
        if (n+zahle .gt. size(e,2))  e=>reallocate(e,3,n+zahle)
        if (n+zahle .gt. size(en,2)) en=>reallocate(en,3,n+zahle)
        if (n+zahle .gt. size(geb))  geb=>reallocate(geb,n+zahle)
        if (n+zahle .gt. size(ep,1)) ep=>reallocate(ep,n+zahle,nnat)
        if (p+zahlk .gt. size(xn))  xn=>reallocate(xn,p+zahlk)
        if (p+zahlk .gt. size(yn))  yn=>reallocate(yn,p+zahlk)
        if (p+zahlk .gt. size(kzi)) kzi=>reallocate(kzi,p+zahlk)
        if (p+zahlk .gt. size(x))   x=>reallocate(x,p+zahlk)
!
!  do refinement
!  mark in liste5, if movement is required for a node.
        allocate(liste5(zahlk+p))
        liste5(:)=0
!      print*,'###  elements/ nodes:',zahle+n,zahlk+p
        print *, '+ New elements...'
        call newelements(ellist,nregen,fehlkn,liste5)
        deallocate(ellist)
!
        call zeit(' Generating Elements')
        print *, '+ Mesh-smoothing...'
        call glaett(n,e,xn,yn,en,geb,kzi,p,nregen,liste5,fehlkn,zrb,    &
     &    zki,xbk,ybk,2)
        deallocate(liste5)
        print *, '+ Pre-assembly...'
        call preassemb
        print*,'+ Assembly and solve...'
        call assemblyandsolve(epsgl,resgl,i)
        if (ende2) then
          print "(a,g10.3)",  " + Reached stopping criterion with error=", error
          print "(a,i3,a,i3)"," + Stopping on iteration ",i, " of", adaptsteps
          print *, '============================================================'
          return
        end if
        if (p .ge. maxnodes) then
          print*,'** reached maximum number of nodes'
          print *, '============================================================'
          ende=.true.
          return
        end if
!  END OF ADAPTSTEPS LOOP:
      end do
!  Not return by stop criteria nor max number of nodes
      print*, 'Completed :', adaptsteps, ' adaptations.' 
      print*, '============================================================'
!
      return
      end subroutine hadapt
!
!
!
      subroutine ffein(ellist,zahle,zahlk,ende,merror,maxelemt,ext,int)
      use feminterface, only: residual, setze, locate, qsortindex
      use globalvariables
      use femtypes
      implicit none
      integer (I4B) :: ellist(:), zahle, zahlk, maxelemt
      real (DP) :: merror
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int
      intent (out) :: ende, merror 
      intent (inout) :: ellist, zahle, zahlk
!
!  Fehlerauswertung und Verfeinerung der Elemente mit den groessten Fehlern
!
!  Eingabe:
!            n        Anzahl der Elemente
!            p        Anzahl der Knoten bzw. Matrixzeilen
!            kdim     vereinbarte Groesse für die Anzahl der Knoten
!            e        Elementinformation, Knoten der Elemente
!            en       Nacbarschaftsinformation der Elemente
!            geb      Zuordnung der Elemente zu den Gebieten 
!            xn,yn    Knotenkoordinaten
!            kzi      Knoten-Zweig-Information; Zweige zu den Knoten
!            kzrb     Zweige deren Randbedingung fuer den Knoten gueltig sind
!            zrb      Randbedingungen der Zweige
!            x        Loesungsvektor
!            matzif   Materialziffern der Gebiete
!            omega    Kreisfrequenz
!            maxelemt maximum number of new elements
!
!  Ausgabe:
!            ende     =.true. wenn maximale Zahl der Knoten oder Elemente erreicht ist
!            merror   relative multi-nature error returned by a posteriori error estimation
!
!  Eingabe/ Ausgabe: 
!            ellist   Liste der zu verfeinernden Elemente und Verfeinerungsart
!            zahle    Anzahl der neu entstehenden Elemente
!            zahlk    Anzahl der neu entstehenden Knoten
!
      integer (I4B) :: i, zahl2, indx(n), errcode
      real (DP) :: errmax, errmin, sumres(nnat), sumref(nnat), res(n,nnat), error(nnat), estimerror(n)
!
!
      call residual(ext,int,.false.,errcode,res,sumres,sumref)
!
      select case(physics)
        case('FLUIDINCOMPR')
          estimerror=sum(res(:,1:2),dim=2)
        case('FLUIDELMASSCON')
          res(:,3)=0._DP
          estimerror=sum(res(:,:),dim=2)
        case default
          estimerror=sum(res(:,:),dim=2)
      end select
!  Sort residual vector. biggest values first. indx points to entry in res.
      call qsortindex(estimerror,indx,n)
      if (zahle .ge. maxelemt) then
        ende=.true.
        return
      end if
      errmax=estimerror(indx(1))
!  Der Mindestfehler zur Verfeinerung wird so festgelegt, dass immer
!  mindestes n/18+1 Elemente verfeinert werden.
      errmin=min(errmax/10.0_DP,estimerror(indx(n/18+1)))
!  Die Mindest- und Hoechstgrenzen sind Daumenwerte!
!  wieviele der Elemente haben einen Fehler > felmin
      zahl2=locate(estimerror(1:n),errmin,indx(1:maxelemt-zahle))
!
      do i=1,zahl2
!  Elemente kennzeichnen
        call setze(ellist,zahle,zahlk,indx(i))
      end do
!
!  return relative error of the solution in %
      select case(physics)
        case('FLUIDINCOMPR')
          error(1:2) = 100._DP * sqrt(sumres/sumref)
          error(3) = 0._DP
        case('FLUIDELMASSCON')
          error(1:2) = 100._DP * sqrt(sumres/sumref)
          error(3) = 0._DP
          error(4:5) = 100._DP * sqrt(sumres/sumref)
        case default
          error = 100._DP * sqrt(sumres/sumref)
      end select
      merror = maxval(error)
!
!  maximal zahl2 Elemente behandeln
      return
      end subroutine ffein
!
!
!
      subroutine wfein(ellist,zahle,zahlk,maxelemt,ende)
      use feminterface, only: setzeg, angles, qsortindex, locate
      use globalvariables, only: e, xn, yn, n, pi
      use femtypes
      implicit none
      integer (I4B) :: ellist(:), zahle, zahlk, maxelemt
      logical :: ende
      intent (in) :: maxelemt
      intent (out) :: ellist, ende
      intent (inout) :: zahle, zahlk
!  Ausspaehen von Dreiecken mit grossem Winkel:
!  Zum Verfeinern markieren.
!
!  Eingabe:
!            n        Anzahl der Elemente
!            p        Anzahl der Knoten bzw. Matrixzeilen
!            kdim     vereinbarte Groesse für die Anzahl der Knoten
!            e        Elementinformation, Knoten der Elemente
!            en       Nachbarschaftsinformation der Elemente
!            geb      Zuordnung der Elemente zu den Gebieten 
!            maxelemt Obere Grenze für die Anzahl der neuen Elemente
!            xn,yn    Knotenkoordinaten
!
!  Ausgabe:
!            ellist   Liste der zu verfeinernden Elemente und Verfeinerungsart
!            ende     =.true. wenn maximale Zahl der Knoten oder Elemente erreicht ist
!
!  Eingabe/ Ausgabe: 
!            zahle    Anzahl der neu entstehenden Elemente
!            zahlk    Anzahl der neu entstehenden Knoten
!
!  local variables 
!            rliste   Real-Hilfsliste der Winkel 
      real (DP) :: wklmax, wklmin, rliste(n), wk(3)
      integer (I4B) :: i, j, zahl2, node, jmax(1)
      integer (I4B) :: elem, indx(n), liste5(n)
      logical :: done
!  Grenze fuer Winkel die verfeinert werden
      parameter (wklmax=125.0_DP/180.0_DP*pi)
!  Grenze fuer Winkel die nicht verfeinert werden
      parameter (wklmin=25.0_DP/180.0_DP*pi)
!
!  Markierungen loeschen
      ellist(1:n)=0
!
!  fuer alle Elemente die maximalen Winkel bestimmen
      do i=1,n
        call angles(i,e,xn,yn,wk(1),wk(2),wk(3))
!  groessten Winkel im Element feststellen
        jmax=maxloc(wk)
        liste5(i)=jmax(1)
        rliste(i)=wk(jmax(1))
!  zu kleine Winkel nicht verfeinern
        if (wk(1) .lt. wklmin ) then
          ellist(i)=ellist(i)+10
        end if
        if (wk(2) .lt. wklmin ) then
          ellist(i)=ellist(i)+20
        end if
        if (wk(3) .lt. wklmin ) then
          ellist(i)=ellist(i)+40
        end if
      end do
!  Liste sortieren Ergebnis als Indexfeld in indx ablegen
      call qsortindex(rliste,indx,n)
!
!  wieviele der Elemente haben Winkel > winkelmax
      zahl2=locate(rliste(1:n),wklmax,indx(1:maxelemt))
!
      ende = .false.
!
!  Alles Null setzen (Anfangswerte)
!  Zahl der zu verfeinernden Elemente und Knoten
      done=.false.
      do while (.not. done)
        done=.true.
        zahle=0
        zahlk=0
!
        do i=1,zahl2
          elem=indx(i)
          node=liste5(elem)
!
!  nun markieren der Elemente
!
          call setzeg(ellist,zahle,zahlk,elem,node)
!
!  falls zuviele Elemente oder Konten eingefuegt wurden:
!  das ganze nochmal aber mit einer kleineren Grenze fuer die
!  Anzahl der zu verfeinernden Elemente
!
          if (zahle.gt.maxelemt) then
            zahl2=i-1
            ende = .true.
!  Liste der Verfeinerungsart zuruecksetzen
!  nur die Information ueber nicht zu verfeinernde Elemente behalten
            do j=1,n
              if (ellist(j).lt.10) then
                ellist(j)=0
              else
!  Einer-stellen abschneiden
                ellist(j)=int(ellist(j)/10)*10
              end if
            end do
            done=.false.
            exit
          end if
!
!  alle Elemente sind bearbeitet
!
        end do
      end do
      return
      end subroutine wfein
!
!
!
      subroutine setze(ellist,zahle,zahlk,elem)
      use feminterface, only: reallocate
      use globalvariables
      use femtypes
      implicit none
      integer (I4B) :: ellist(:), zahle, zahlk, elem
      intent (in) :: elem
      intent (inout) :: ellist,zahle,zahlk
!  Eingabe:
!            en       Nachbarschaftsinformation der Elemente
!            geb      Zuordnung der Elemente zu den Gebieten 
!            elem     zu verfeinernde Elementnummer
!
!  Ausgabe:
!
!  Eingabe/ Ausgabe: 
!            ellist   Verfeinerungsart der zu verfeinernden Elemente
!            zahle    Anzahl der neu entstehenden Elemente
!            zahlk    Anzahl der neu entstehenden Knoten
!
!  lokale Variablen
      integer (I4B) :: i, j, nbnode, nachb, stack, l3el, l3sgel
      integer (I4B) :: l3nb, l3sgnb, element
!            liste    Hilfsfeld wird als stack benutzt
      integer (I4B), pointer :: liste(:)
!
!  "rote" Verfeinerung eines Dreiecks in vier Dreieck und
!  Verfeinerung der angrenzenden Dreiecke soweit notwendig
!  Die Verfeinerungsart wird nur in der ellist markiert und noch nicht
!  auf dem Netz ausgefuehrt
!
!  liste dient als Stack und zahl als Tiefe des Stacks (pointer)
!  der Stack wird benutzt um noch zu verfeinernde Elemente zu halten
!
!  die Eintraege in ellist koennen die folgenden Werte annehmen:
!     9    Volle Verfeinerung abgeschlossen (rote Verfeinerung)
!     0    Noch keine Verfeinerung vorgesehen
!     1    Gruene Verfeinerung gegenueber Knoten 1
!     2    Gruene Verfeinerung gegenueber Knoten 2
!     3    Gruene Verfeinerung gegenueber Knoten 3
!     8    Rote Verfeinerung, nur temporaer bedingt durch zwei
!          angrenzende Elemente mit roter Verfeinerung
!
!   fuer die Zehnerstellen:
!    10    gegenueber Knoten 1 soll nicht verfeinert werden
!    20    gegenueber Knoten 2 soll nicht verfeinert werden
!    40    gegenueber Knoten 3 soll nicht verfeinert werden
!    30    gegenueber Knoten 1 und 2 soll nicht verfeinert werden
!    50    gegenueber Knoten 1 und 3 soll nicht verfeinert werden
!    60    gegenueber Knoten 2 und 3 soll nicht verfeinert werden
!    70    gegenueber Knoten 1, 2 und 3 soll nicht verfeinert werden
!
      allocate(liste(20))
      liste(1)=elem
      stack =1
!
      do while (stack.gt.0)
!  weiter falls noch etwas auf dem Stack ist
        element=liste(stack)
        stack=stack-1
        l3el=int(ellist(element)/10)
        l3sgel=ellist(element) - l3el*10
        if (l3sgel.eq.9) then
!  Element wird bereits maximal verfeinert? Dann keine
!  weiteren Operationen erforderlich.
          cycle
!  Verfeinere das Element
        else if (l3sgel.eq.0 .or. l3sgel.eq.8) then
!  das Element war bisher noch unverfeinert
          zahle=zahle+3
        else
!  das Element war bereits gruen verfeinert
          zahle=zahle+2
        end if
!  Knotenzahl wird beim besuchen der Nachbarelemente korrigiert
        ellist(element)=l3el*10+9
!
!  Besuchen der Nachbarn
!
        do i=1,3
          nachb=en(i,element)
          if (nachb.le.0) then
!  Zweig liegt auf dem Rand
            if (l3sgel.ne.i) then
!  dann neuen Knoten nur, wenn das Randstueck noch nicht verfeinert wurde
              zahlk=zahlk+1
            end if
            cycle
          end if
!
!  lokale Knotennummer nbnode im Nachbarelement bestimmen
          do j=1,3
            if (en(j,nachb).eq.element) then
              nbnode=j
              exit
            end if
          end do
!  Verfeinerungsart holen
          l3nb=int(ellist(nachb)/10)
          l3sgnb=ellist(nachb) - l3nb*10
!
          if (l3sgnb .eq. 9) then 
!  Nachbar war bereits maximal verfeinert
            cycle
          else if (l3sgnb .eq. 8) then 
!  Nachbar ist bereits fuer maximale Verfeinerung vorgesehen
            zahlk=zahlk+1
            cycle
          else if (l3sgnb .eq. 0) then 
!  Nachbar ist noch unverfeinert
            zahlk=zahlk+1
            if (((nbnode.eq.1) .and.                                    &
     &        ((l3nb.eq.1).or.(l3nb.eq.3).or.(l3nb.eq.5)))              &
     &        .or.((nbnode.eq.2) .and.                                  &
     &        ((l3nb.eq.2).or.(l3nb.eq.3).or.(l3nb.eq.6)))              &
     &        .or.((nbnode.eq.3) .and.                                  &
     &        ((l3nb.eq.4).or.(l3nb.eq.5).or.(l3nb.eq.6)))) then
!  Nachbar ist noch unverfeinert, aber der Zweig soll nicht Gruen verfeinert werden
!  Verfeinere ihn also rot
              ellist(nachb)=l3nb*10+8
              stack=stack+1
              if (stack .ge. size(liste)) then
!  reallocation of stack needed
                liste=>reallocate(liste,2*size(liste))
              end if
              liste(stack)=nachb
            else
!  Verfeinere den Nachbarn gruen
              ellist(nachb)=l3nb*10+nbnode
              zahle=zahle+1
            end if
          else if (l3sgnb .eq. nbnode) then 
!  Nachbar war schon wie gewuenscht verfeinert
          else
!  Nachbar ist bereits einmal verfeinert dann verfeinere ihn rot
            if (en(l3sgnb,nachb).gt.0) then
!  kein neuer Knoten, wenn Nachbar auf dem Rand liegt und der Randzweig bereits verfeinert war 
              zahlk=zahlk+1
            end if
            zahle=zahle-1
            stack=stack+1
            if (stack .ge. size(liste)) then
!  reallocation of stack needed
              liste=>reallocate(liste,2*size(liste))
            end if
            liste(stack)=nachb
            ellist(nachb)=l3nb*10+8
          end if
!
          if (geb(nachb).ne.geb(element)) then
!  Wenn der Nachbar in einem anderen Gebiet liegt, also eine Gebietsgrenze
!  ueberschritten wird, dann soll der Nachbar maximal verfeinert werden.
            stack=stack+1
            if (stack .ge. size(liste)) then
!  reallocation of stack needed
              liste=>reallocate(liste,2*size(liste))
            end if
            liste(stack)=nachb
          end if
!
        end do
!  end Besuchen der Nachbarn
      end do
      deallocate(liste)
      return
      end subroutine setze
!
!
!
      subroutine setzeg(ellist,zahle,zahlk,elem,node)
      use feminterface, only: setze
      use globalvariables
      use femtypes
      implicit none
      integer (I4B) :: ellist(:), zahle, zahlk, elem, node
      intent (in) :: elem, node
      intent (inout) :: ellist, zahle, zahlk
!  Eingabe:
!            en       Nacbarschaftsinformation der Elemente
!            geb      Zuordnung der Elemente zu den Gebieten 
!            elem     zu verfeinernde Elementnummer
!            node     lokale Zweignummer des zu unterteilenden Elementes
!
!  Ausgabe:
!
!  Eingabe/ Ausgabe: 
!            ellist   Liste der zu verfeinernden Elemente und Verfeinerungsart
!            zahle    Anzahl der neu entstehenden Elemente
!            zahlk    Anzahl der neu entstehenden Knoten
!
!  lokale Variablen
      integer (I4B) :: nachb, j, nbnode, l3el, l3nb, l3sgel, l3sgnb
!  Gruene Verfeinerung des Elementes elem gegenueber des Knotens node
!  und Verfeinerung des Nachbarelementes
!  Falls eines der Elemente bereits verfeinert wurde, wird eine
!  rote Verfeinerung durchgefuehrt
!
!  Pruefen ob verfeinert werden darf
!
!  Zehnerstellen abseparieren
!
      l3el=int(ellist(elem)/10)
!  Einerstellen abseparieren
      l3sgel=ellist(elem) - l3el*10
!
      if (((node.eq.1) .and.                                            &
     &  ((l3el.eq.1).or.(l3el.eq.3).or.(l3el.eq.5)))                    &
     &  .or.((node.eq.2) .and.                                          &
     &  ((l3el.eq.2).or.(l3el.eq.3).or.(l3el.eq.6)))                    &
     &  .or.((node.eq.3) .and.                                          &
     &  ((l3el.eq.4).or.(l3el.eq.5).or.(l3el.eq.6)))) then
!  the edge: node is supposed to be refined, 
!  however the edge is as well excludes from
        return
      end if
!
      nachb=en(node,elem)
!  Knoten im Nachbarelement suchen
      if (nachb.gt.0) then
        do j=1,3
          if (en(j,nachb).eq.elem) then
!  Nachbar
            nbnode=j
            exit
          end if
        end do
!  Zehnerstellen abseparieren
        l3nb=int(ellist(nachb)/10)
!  Einerstellen abseparieren
        l3sgnb=ellist(nachb) - l3nb*10
!
        if (((nbnode.eq.1) .and.                                        &
     &    ((l3nb.eq.1).or.(l3nb.eq.3).or.(l3nb.eq.5)))                  &
     &    .or.((nbnode.eq.2) .and.                                      &
     &    ((l3nb.eq.2).or.(l3nb.eq.3).or.(l3nb.eq.6)))                  &
     &    .or.((nbnode.eq.3) .and.                                      &
     &    ((l3nb.eq.4).or.(l3nb.eq.5).or.(l3nb.eq.6)))) then
!  the edge: node (nbnode in the adjacent element) is supposed to be refined, 
!  however the setting of teh adjacent elment forbid refinement
          return
        end if
      else
        nbnode=0
      end if
!
      if (l3sgel.eq.0) then
        if (nachb.gt.0) then
          if (l3sgnb.eq.0) then
!  Fall 1
!  elem unverfeinert, Nachbar existiert und unverfeinert
!  gruene Verfeinerung von elem und des Nachbarn
!         print*,'fall 1',elem
            zahle=zahle+2
            zahlk=zahlk+1
            ellist(elem)=l3el*10+node
            ellist(nachb)=l3nb*10+nbnode
          else
!  Fall 2
!  elem unverfeinert, Nachbar existiert und ist bereits verfeinert
!  rote Verfeinerung des Nachbarn, damit wird elem automatisch gruen
!         print*,'fall 2'
            call setze(ellist,zahle,zahlk,nachb)
          end if
        else
!  Fall 3
!  elem unverfeinert, kein Nachbar existiert
!  gruene Verfeinerung von elem
!         print*,'fall 3'
          zahle=zahle+1
          zahlk=zahlk+1
          ellist(elem)=l3el*10+node
        end if
      else
        if (nachb.gt.0) then
          if (l3sgnb.eq.0) then
!  Fall 4
!  elem beseits verfeinert, nachbar existiert und unverfeinert
!  rote Verfeinerung von elem, damit wird nachb automatisch gruen
!         print*,'fall 4'
            call setze(ellist,zahle,zahlk,elem)
          else
!  Fall 5
!  elem beseits verfeinert, Nachbar existiert und ist bereits verfeinert
!  pruefe ob die beiden elemente schon richtig (doppelgruen) verf. sind
!     wenn ja mache nichts
!     sonst verfeinere einen der beiden rot
            if (l3sgel.eq.node .and. l3sgnb.eq.nbnode) then
!              print*,'fall 5a'
!  mache nichts
            else
!              print*,'fall 5b'
              call setze(ellist,zahle,zahlk,elem)
            end if
          end if
        else
!  Fall 6
!  elem bereits verfeinert, kein Nachbar existiert
!  pruefe ob das bereits die richtige Verfeinerung ist, sonst:
!  rote Verfeinerung von elem
!         print*,'fall 6'
          if (l3sgel.ne.node) then
            call setze(ellist,zahle,zahlk,elem)
          end if
        end if
      end if
      return
      end subroutine setzeg
!
!
!
      subroutine newelements(ellist,nregen,fehlkn,liste5)
      use feminterface, only: angles, nachb, gennode
      use globalvariables, only: e, en, p, n, xn, yn, ep, geb, nnat, pi
      use femtypes
      implicit none
      integer (I4B) :: fehlkn, ellist(:), liste5(:)
      logical :: nregen
      intent (in) ::  ellist
      intent (out) :: nregen, fehlkn, liste5
!
!  Perform the mesh refinement whether red / green or other refinement types is to be done
!  is stored in ellist. As a result of the refinement all the element data are modified
!
!  Input:
!            zki      branch information, geometry nodes which determine the branches
!            xbk      coordinates of the geometry node
!            ybk
!            ellist   list of the elements to refine and the type of refinement
!  Output:
!            nregen   .true. if the mesh needs regeneration due to the occurence of elements 
!                      with a negative area
!            fehlkn   number of nodes which have to be relocated, nodes are leisted in liste5
!            liste5   =1 if the node leads to a negative element area and has to be relocated
!  In-/ Output:
!            n        total number of elements
!            p        total number of nodes
!            e        element information, nodes of the elements
!            en       neighborhood information, neighbors of the elements
!            geb      assignment of elements to the geometry regions
!            xn, yn   Coordinates of the nodes
!            x        solution vector
!            kzi      node branch information
! 
!  local variables      
      integer (I4B) :: i, j, l3sgel, jmax(1), oldn, j1, j2, j3
      integer (I4B) :: isucc(3), ipred(3), nodelist(3,n)
      real (DP) :: zweipi, wk(3)
      parameter (zweipi=2._DP*pi)
      parameter (isucc=(/2,3,1/))  !  successor
      parameter (ipred=(/3,1,2/))  !  predecessor
!
      nregen=.false.
      fehlkn=0
      nodelist=0
!  nodelist stores for each of the element edges the index 
!  of the new node generated on that edge
!  the entry has a negative sign if the node has been positioned
!  on a circular geometry branch (arc)
      oldn=n
!
!  first part: generate the new nodes
!
      do i=1,oldn
        l3sgel=ellist(i) - int(ellist(i)/10)*10
        if (l3sgel.eq.9) then
!  red refinement:
!  generate 3 new nodes and 3 new elements from this element
!
!  check whether the edges are located on a geometry branch
!  and generate the new nodes ...
          do j=1,3
            call gennode(i,j,nodelist,liste5,nregen,fehlkn)
          end do
        else if (l3sgel.gt.0 .and. l3sgel.le.3) then
!  green refinement:
!  generate 1 new node and 1 new element from this element
!
!  check whether the edge is located on a geometry branch
!
!  generate one new node ...
          j=l3sgel
          call gennode(i,j,nodelist,liste5,nregen,fehlkn)
        else if (l3sgel.gt.3 .and. l3sgel.le.6) then
!  yellow refinement:
!  generate 2 new nodes and 2 new elements from this element
!
!  yellow is presently not fully implemented and not tested
          stop
        end if
      end do
!
!  second part: generate the new elements
!
      do i=1,oldn
        l3sgel=ellist(i) - int(ellist(i)/10)*10
        if (l3sgel.eq.9) then
!  red refinement
!  generate 3 new nodes and 3 new elements from this element
!
!  if one of the nodes had been positioned on an arc then choose the secure variant
!  else decide by computing the angles
          if (minval(nodelist(1:3,i)) .lt. 0) then
            jmax=minloc((/nodelist(1,i),nodelist(2,i),nodelist(3,i)/))
            wk(jmax(1))=pi
          else
            call angles(i,e,xn,yn,wk(1),wk(2),wk(3))
            jmax=maxloc(wk)
          end if
!  generate the elements 
          if (wk(jmax(1)) .lt. pi/2._DP) then
!  none of the angles is larger than 90 degrees    
!      /\
!     /\/\
            e(1,n+1)=abs(nodelist(3,i))
            e(2,n+1)=e(2,i)
            e(3,n+1)=abs(nodelist(1,i))
            geb(n+1)=geb(i)
            ep(n+1,1:nnat)=ep(i,1:nnat)
            e(1,n+2)=abs(nodelist(1,i))
            e(2,n+2)=abs(nodelist(2,i))
            e(3,n+2)=abs(nodelist(3,i))
            geb(n+2)=geb(i)
            ep(n+2,1:nnat)=ep(i,1:nnat)
            e(1,n+3)=abs(nodelist(2,i))
            e(2,n+3)=abs(nodelist(1,i))
            e(3,n+3)=e(3,i)
            geb(n+3)=geb(i)
            ep(n+3,1:nnat)=ep(i,1:nnat)
            e(2,i)=abs(nodelist(3,i))
            e(3,i)=abs(nodelist(2,i))
          else 
!  one of the angles is larger than 90 degrees, subdivide this angle
!      /|\
!     /\|/\
            j1=jmax(1)
            j2=isucc(j1)
            j3=ipred(j1)
            e(1,n+1)=abs(nodelist(j2,i))
            e(2,n+1)=e(j1,i)
            e(3,n+1)=abs(nodelist(j1,i))
            geb(n+1)=geb(i)
            ep(n+1,1:nnat)=ep(i,1:nnat)
            e(1,n+2)=abs(nodelist(j1,i))
            e(2,n+2)=e(j1,i)
            e(3,n+2)=abs(nodelist(j3,i))
            geb(n+2)=geb(i)
            ep(n+2,1:nnat)=ep(i,1:nnat)
            e(1,n+3)=abs(nodelist(j1,i))
            e(2,n+3)=abs(nodelist(j3,i))
            e(3,n+3)=e(j2,i)
            geb(n+3)=geb(i)
            ep(n+3,1:nnat)=ep(i,1:nnat)
            e(1,i)=e(j3,i)
            e(2,i)=abs(nodelist(j2,i))
            e(3,i)=abs(nodelist(j1,i))
          end if
          n=n+3
        else if (l3sgel.gt.0 .and. l3sgel.le.3) then
!  green refinement
!  generate 1 new node and 1 new element from this element
!
          j1=l3sgel
          j2=isucc(j1)
          j3=ipred(j1)
          e(j1,n+1)=e(j3,i)
          e(j2,n+1)=e(j1,i)
          e(j3,n+1)=abs(nodelist(j1,i))
          geb(n+1)=geb(i)
          ep(n+1,1:nnat)=ep(i,1:nnat)
          e(j3,i)=abs(nodelist(j1,i))
          n=n+1
        else if (l3sgel.gt.3 .and. l3sgel.le.6) then
!  yellow refinement
!  generate 2 new nodes and 2 new element from this element
!  do not generate a node on the side l3sgel-3
!         /|
!        /_|
!      / \ |
          j1=l3sgel-3
          j2=isucc(j1)
          j3=ipred(j1)
          e(1,n+1)=e(j3,i)
          e(2,n+1)=nodelist(j2,i)
          e(3,n+1)=e(j2,i)
          geb(n+1)=geb(i)
          ep(n+1,1:nnat)=ep(i,1:nnat)
          e(1,n+2)=e(j2,i)
          e(2,n+2)=nodelist(j2,i)
          e(3,n+2)=nodelist(j3,i)
          e(1,i)=e(j1,i)
          e(2,i)=nodelist(j3,i)
          e(3,i)=nodelist(j2,i)
          geb(n+2)=geb(i)
          ep(n+2,1:nnat)=ep(i,1:nnat)
          n=n+2
        end if
!  check whether all elements have positive area
      end do
!  regenerate the neighbors
      call nachb(e,en,p,n)
!
      return
      end subroutine newelements
!
!
!
      subroutine gennode(i,j,nodelist,liste5,nregen,fehlkn)
      use feminterface, only: neukn
      use globalvariables
      use femtypes
      implicit none
      integer (I4B) :: i, j, nodelist(:,:), liste5(:), fehlkn
      logical :: nregen 
      intent (in) :: i, j
      intent (out) :: nodelist, liste5, nregen
      intent (inout) :: fehlkn
!
!  Generate -if not already existing- a new node at the edge j of the element i
!  determine the properties and coordinated of this node
!
!  Input:
!            i        element index to process
!            j        node index to process
!            e        element information, nodes of the elements
!            en       neighborhood information, neighbors of the elements
!            geb      assignment of elements to the geometry regions
!            xbk      coordinates of the geometry node
!            ybk
!            zki      branch information, geometry nodes which determine the branches
!   Output:
!            nodelist
!            liste5   =1 if the node leads to a negative element area and has to be relocated
!            nregen   .true. if the mesh needs regeneration due to the occurence of elements 
!                      with a negative area
!   In-/ Output:
!            p        total number of nodes
!            xn, yn   Coordinates of the nodes
!            x        solution vector
!            kzi      node branch information
!            fehlkn   number of nodes which have to be relocated, nodes are leisted in liste5
 
!              
!            zki      Zweig-Knoten-Information, Eingabeknoten der Eingabezweige
!            xbk,ybk  Koordinaten der Eingabeknoten
!            ez       Anzahl der Elemente vor der Verfeinerung
!
!  local variables
      integer (I4B) :: nbr, n1, n2, n4, geb1, geb2, k, neuk, err
      integer (I4B) :: isucc(3), ipred(3)
      logical :: arc
      parameter (isucc=(/2,3,1/))  !  successor
      parameter (ipred=(/3,1,2/))  !  predecessor
!
      geb1=geb(i)
      nbr=en(j,i)
      if ((i.lt.nbr) .or. (nbr.le.0)) then
!  ... if not already existing
!  otherwise the node will be accessed via the neighboring element
        n1=e(isucc(j),i)
        n2=e(ipred(j),i)
        if (nbr .ne. 0) then
          geb2=geb(nbr)
        else
          geb2=0
        end if
!  determine the node in the neighbor element
        if (nbr.ne.0) then
          do k=1,3
            if (en(k,nbr) .eq. i) then
              n4=k
              exit
            end if
          end do
        else
          n4=0
        end if
        if (nbr.ne.0) then
          call neukn(n1,n2,e(j,i),e(n4,nbr),neuk,xn,yn,x,kzi,zki,p,     &
     &      xbk,ybk,geb1,geb2,arc,err)
        else 
          call neukn(n1,n2,e(j,i),0,neuk,xn,yn,x,kzi,zki,p,             &
     &      xbk,ybk,geb1,geb2,arc,err)
        end if
        if (err .ne. 0) then 
!  an error occured during the positioning of the new node 
!  use the standard positioning and set the flag for the smoothing procedure
!  this will aim to push the node toward the correct position 
          xn(p)=(xn(n1)+xn(n2))/2._DP
          yn(p)=(yn(n1)+yn(n2))/2._DP
          nregen=.true.
          liste5(p)=1
          fehlkn=fehlkn+1
        end if
        if (arc) then
          nodelist(j,i)=-p
        else
          nodelist(j,i)=p
        end if
!  set the node information in the neighbor element as well
        if (nbr.ne.0) then
          if (arc) then 
            nodelist(n4,nbr)=-p
          else
            nodelist(n4,nbr)=p
          end if
        end if
      end if
      return
      end subroutine gennode
