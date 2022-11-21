      subroutine basin(ok,jd,n,p,e,xn,yn,kzi,geb,en)
      use feminterface, only: nachb, getsetting
      use femtypes
      implicit none
      integer (I4B) n, p 
      integer (I4B), pointer :: e(:,:), kzi(:), geb(:), en(:,:)
      real (DP) :: jd
      real (DP), pointer :: xn(:), yn(:)
      logical ok
      intent (out) ok, n, p, jd
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
!    $Revision: 1.8 $
!    $Date: 2014/02/11 16:23:53 $
!    $Author: juryanatzki $
!
!  Einlesen der Basisdaten aus dem File Basis
!  Wenn das Basis-file die Nachbarn nicht enthaelt, werden diese automatisch
!  bestimmt.
!  Ein-/Ausgabe:
!     ok      =.false. wenn ein Fehler aufgetreten ist
!     n       Anzahl der Elemente
!     p       Anzahl der Knoten
!     e       Elementinformation, Knoten der Elemente
!     xn,yn   Knotenkoordinaten
!     kzi     Knoten-Zweig-Information
!            /   = 0  : innerer Knoten
!     kzi    -   < 0  : ist mit Gebietsknoten Nr.(-kzi) identisch
!            \   > 0  : liegt auf Zweig der Nr. (kzi)
!     geb     Gebiete der Elemente
!     en      Nachbarn der Elemente
!     jd
!
! lokale Variablen
      integer (I4B) unitid, ios
      integer (I4B) i
      logical enwri
      character (len=200) path
!
!  File oeffnen

      call getsetting('PROJECTPATH',path)
      inquire(file=path(1:len_trim(path))//'basin',iostat=ios)
      if (ios .ne. 0) then
        print*,'***** File of base data does not exist'
        ok=.false.
        return
      end if
      call grglun(unitid)
      print*, 'Reading the basis file'
     open (unitid,file=path(1:len_trim(path))//'basis',                &
     &  form='unformatted',position='REWIND',action='READ',iostat=ios)
      if (ios .ne. 0) then
        print*,'Error opening file: ',path(1:len_trim(path))//'basis'
        print*,'Error no: ',ios
        ok=.false.
        return
      end if

!
      read (unitid,end=1999,err=1999) jd, n, p, enwri
!
!  wenn enwri true ist sind die Nachbarschaftsinformationen auf dem file
!  fuer alle Elemente die Knotennummern, Nachbarn und die Gebietszugehoerigkeit
!      
      allocate (xn(p), yn(p), kzi(p), e(3,n), en(3,n), geb(n))
      if (enwri) then
       do i=1,n
          read (unitid,end=1999,err=1999) e(:,i),en(:,i),geb(i)
        end do
      else
        do i=1,n
          read (unitid,end=1999,err=1999) e(1:3,i), geb(i)
        end do
!  Feld en der Nachbarn generieren
        call nachb(e,en,p,n)
      end if
      do i=1,p
        read (unitid,end=1999,err=1999) xn(i),yn(i),kzi(i)
      end do
!
!      if (enwri) then
! TO DO: use IOSTAT
!        read (unitid,end=1999,err=1999) e(:,1:n),en(:,1:n),geb(1:n)
!      else
!        read (unitid,end=1999,err=1999) e(:,1:n), geb(1:n)
!  Feld en der Nachbarn generieren
!        call nachb(e,en,p,n)
!      end if
!      read (unitid,end=1999,err=1999) xn(1:p),yn(1:p),kzi(1:p)
!
      close(unitid,iostat=ios)
      if (ios .ne. 0) then
        print*,'Error closing file: ',path(1:len_trim(path))//'basis'
        print*,'Error no: ',ios
        ok=.false.
        return
      end if
!
!  fuer alle Knoten die Koordinaten und die zugehoerige Zweignummer lesen
      ok=.true.
999   return
1999  print*,'***** Error in reading basis'
      stop
      end subroutine basin
