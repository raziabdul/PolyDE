      subroutine basout(ok,n,p,e,xn,yn,kzi,geb,en)
      use feminterface, only: getsetting, timestamp
      use femtypes
      implicit none
      integer (I4B) e(:,:),n,p,kzi(:),geb(:),en(:,:)
      real (DP) xn(:),yn(:)
      logical ok
      intent (in) :: n, p, e, xn, yn, kzi, geb, en
      intent (out) :: ok
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
!    $Revision: 1.6 $
!    $Date: 2006/07/07 21:08:35 $
!    $Author: r_abdul $
!
!  Ausgabe der Basisdaten auf das File Basis
!
!  Input
!     n        Anzahl der Elemente
!     p        Anzahl der Knoten
!     e        Elementinformation
!     xn, yn   Knotenkoordinaten
!     kzi      Knoten-Zweig-Information; Zweige zu den Knoten
!     geb      Gebiete der Elemente
!     en       Nacbarschaftsinformation der Elemente
!
!  Output
!     ok       =.false.  if an I/O Error occured
!
!     enwri    =.true.  fuer alle Elemente die Knoten, Nachbarn und die Gebietsnummer
!              =.false. fuer alle Elemente die Knoten und die Gebietsnummer
!  lokale variablen
!
      real (DP) :: jd
      logical enwri
      integer (I4B) unitid, ios
      integer (I4B) i
      character (len=200) path
!      
      external :: grglun
!
      jd = timestamp()
      enwri=.false.
!
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
      open (unitid,file=path(1:len_trim(path))//'basis',                &
     &  form='unformatted',position='REWIND',action='WRITE',iostat=ios)
      if (ios .ne. 0) then
        print*,'Error opening file: ',path(1:len_trim(path))//'basis'
        print*,'Error no: ',ios
        ok=.false.
        return
      end if
      rewind(unitid)
      write(unitid,err=1999) jd, n, p, enwri
!  wenn enwri true ist werden die Nachbarschaftsinformationen geschrieben
!
      if (enwri) then
        do i=1,n
          write (unitid,err=1999) e(:,i),en(:,i),geb(i)
        end do
      else
        do i=1,n
          write (unitid,err=1999) e(:,i), geb(i)
        end do
      end if
      do i=1,p
        write (unitid,err=1999) xn(i),yn(i),kzi(i)
      end do
!
!      if (enwri) then
!        write (unitid,err=1999) e(:,1:n), en(:,1:n), geb(1:n)
!      else
!        write (unitid,err=1999) e(:,1:n), geb(1:n)
!      end if
!      write (unitid,err=1999) xn(1:p), yn(1:p), kzi(1:p)
!
      close(unitid,iostat=ios)
      if (ios .ne. 0) then
        print*,'Error closing file: ',path(1:len_trim(path))//'basis'
        print*,'Error no: ',ios
        ok=.false.
        return
      end if
      ok=.true.
!
      return
1999  print*,'***** Error in writing basis'
      stop
      end subroutine basout
!
!
!
