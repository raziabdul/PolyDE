!
      subroutine dichk(n,p,e,en,xn,yn,ok,geb,kzi)
      use feminterface, only: zeit
      use femtypes
      implicit none
      integer (I4B) n,p,e(:,:),en(:,:)
      integer (I4B), optional :: geb(:), kzi(:)
      real (DP) xn(:),yn(:)
      logical ok
      intent (in) :: n, p, e, xn, yn, geb, kzi
      intent (out) :: ok
      intent (inout) :: en
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
!    $Date: 2006/07/07 21:08:35 $
!    $Author: r_abdul $
!
!  local variables
      integer (I4B) i, j, k, n1, zahl, inach(3), j2, j3, geb1, geb2, nbr
      real (DP) area
      parameter (inach=(/2,3,1/))
!
! dichk= Data Integrity Check
!
      ok=.true.
      do i=1,n
        do j=1,3
          if ((e(j,i).le.0).or.(e(j,i).gt.p)) then
            write (*,1) i,e(1,i),e(2,i),e(3,i)
1           format(' dichk: Element',i6,' Knoten:',3i6)
            write (*,2)
2           format(' falsche Knotennummer')
            ok=.false.
            pause
          end if
          if (e(j,i).eq.e(inach(j),i)) then
            write (*,1) i,e(1,i),e(2,i),e(3,i)
            write (*,3)
3           format(' Knotennummer doppelt')
            ok=.false.
            pause
          end if
        end do
        do j=1,3
          n1=en(j,i)
          if (n1 .gt. n) then 
            write (*,9) i,n1
9           format(' dichk: Element',i6,' hat Nachbar',i6,'groesser als   &
     &        Elementzahl')
            ok=.false.
          end if
          if (n1.gt.i) then
            zahl=0
            do k=1,3
              if (en(k,n1).eq.i) zahl=k
            end do
            if (zahl.eq.0) then
              write (*,10) i,n1
10            format(' dichk: Element',i6,' Nachbar',i6)
              write (*,11)
11            format(' Fehler in Nachbarliste')
              write (*,12) i,e(1,i),e(2,i),e(3,i),                      &
     &          en(1,i),en(2,i),en(3,i)
              write (*,12) n1,e(1,n1),e(2,n1),e(3,n1),                  &
     &          en(1,n1),en(2,n1),en(3,n1)
12            format(' Element:',i6,' Knoten:',3i6,' Nachbarn:',3i6)
              ok=.false.
            end if
          end if
        end do
        area=((xn(e(2,i))-xn(e(3,i)))*(yn(e(3,i))-yn(e(1,i)))-          &
     &      (xn(e(3,i))-xn(e(1,i)))*(yn(e(2,i))-yn(e(3,i))))/2.0_DP
        if (area.le.0._DP) then
          write (*,1) i,e(1,i),e(2,i),e(3,i)
          write (*,6) en(1,i),en(2,i),en(3,i)
6         format(' Nachbarn: ',3i6)
          write (*,4) area
4         format(' Flaeche negativ oder Null:',e20.10)
          write (*,5)
5         format(' Knoten x y:')
          write (*,*) e(1,i),xn(e(1,i)),yn(e(1,i))
          write (*,*) e(2,i),xn(e(2,i)),yn(e(2,i))
          write (*,*) e(3,i),xn(e(3,i)),yn(e(3,i))
          ok=.false.
        end if
        if (present(geb)) then
          do j=1,3
            j2=inach(j)
            j3=inach(j2)
            geb1=geb(i)
            nbr=en(j,i)
            if (nbr .gt. 0) then
              geb2=geb(nbr)
            else
              geb2=0
            end if
            if (geb1 .eq. geb2) then
              if ((kzi(e(j2,i)) .eq. kzi(e(j3,i))) .and.                &
     &          (kzi(e(j2,i)).ne.0) .and.                               &
     &          (kzi(e(j2,i)).ne.kzi(e(j3,i)))) then
                write(*,21) i,kzi(e(1,i)),kzi(e(2,i)),kzi(e(3,i))
21              format (' Element',i6,' Knoten Zweig Information',      &
     &            ' falsch gesetzt',3i6)
                write(*,6) en(1,i),en(2,i),en(3,i)
                write(*,22) j,geb1,geb2
22              format (' Zweig: ',i6,' Gebiet Element: ',i6,' Gebiet', &
     &            ' des Nachbarn: ',i6)
                write(*,*)
                ok=.false.
                pause
              end if
            else
              if (kzi(e(j2,i)) .eq. 0 .or. kzi(e(j3,i)) .eq. 0) then
                write(*,21) i,kzi(e(1,i)),kzi(e(2,i)),kzi(e(3,i))
                write(*,6) en(1,i),en(2,i),en(3,i)
                write(*,22) j,geb1,geb2
                write(*,*)
                ok=.false.
              end if
            end if
          end do
        end if
      end do
      call zeit(' dichk: data integrity check 1.2')
      return
      end subroutine dichk
!
!
!
