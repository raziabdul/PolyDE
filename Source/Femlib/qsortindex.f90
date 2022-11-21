      subroutine qsortindex_dp(arr,indx,n)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: arr(:)
      integer (I4B) :: indx(:), n
      intent(in) :: arr, n
      intent(out) :: indx
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
!    $Date: 2011/08/16 15:07:03 $
!    $Author: m_kasper $
!
!   Generate a sorted list of indices
!   Input:      arr      array to be sorted (key)
!               n        number of elements in arr
!   Output:     indx     sorted list of indices (largest element comes first)
!                        
      real(DP) :: a
      integer(I4B), parameter :: nn=10
      integer(I4B) k, i, j, indext, jstack, l, r, swp, stacksize, im
      integer(I4B), allocatable :: istack(:)
!
!  Sorting algorithm:   Quick - Sort (adapted from Numerical Recipes)
!
!  sub-arrays of length up to nn are sorted by straight insertion
!
      stacksize=int(2*log(real(n))/log(real(2)))+2
      allocate (istack(stacksize))
      indx=(/(i, i=1,n)/) 
      jstack=0
      l=1
      r=n
      do
        if (r-l .lt. nn) then
!  sort by straight insertion
          do j=r-1,l,-1
            indext=indx(j)
            a=arr(indext)
            im=r+1
            do i=j+1,r
!              if (arr(indx(i)) .ge. a) then ! inverse order
              if (arr(indx(i)) .le. a) then
                im=i
                exit
              end if
              indx(i-1)=indx(i)
            end do
            indx(im-1)=indext
          end do
          if (jstack .eq. 0) then
            deallocate (istack)
            return
          end if
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+r)/2
!  sawp i and  l+1
          swp=indx(k)
          indx(k)=indx(l+1)
          indx(l+1)=swp
!          if (arr(indx(r)) .lt. arr(indx(l))) then ! inverse order
          if (arr(indx(r)) .gt. arr(indx(l))) then
            swp=indx(l)
            indx(l)=indx(r)
            indx(r)=swp
          end if
!          if (arr(indx(r)) .lt. arr(indx(l+1))) then ! inverse order
          if (arr(indx(r)) .gt. arr(indx(l+1))) then
            swp=indx(l+1)
            indx(l+1)=indx(r)
            indx(r)=swp
          end if
!          if (arr(indx(l+1)) .lt. arr(indx(l))) then ! inverse order
          if (arr(indx(l+1)) .gt. arr(indx(l))) then
            swp=indx(l)
            indx(l)=indx(l+1)
            indx(l+1)=swp
          end if
          i=l+1
          j=r
          indext=indx(l+1)
          a=arr(indext)
          do
            do
              i=i+1
!              if (arr(indx(i)) .ge. a) exit ! inverse order
              if (arr(indx(i)) .le. a) exit
            end do
            do
              j=j-1
!              if (arr(indx(j)) .le. a) exit ! inverse order
              if (arr(indx(j)) .ge. a) exit
            end do
            if (j .lt. i) exit
!  sawp i and  j
            swp=indx(i)
            indx(i)=indx(j)
            indx(j)=swp
          end do
          indx(l+1)=indx(j)
          indx(j)=indext
          jstack=jstack+2
          if (r-i+1 .ge. j-l) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          end if
        end if
      end do
      end subroutine qsortindex_dp
!
!
!
      subroutine qsortindex_sp(arr,indx,n)
      use feminterface, only:
      use femtypes
      implicit none
      real (SP) :: arr(:)
      integer (I4B) :: indx(:), n
      intent(in) :: arr, n
      intent(out) :: indx
!
!   Generate a sorted list of indices
!   Input:      arr      array to be sorted (key)
!               n        number of elements in arr
!   Output:     indx     sorted list of indices (largest element comes first)
!                        
      real(SP) :: a
      integer(I4B), parameter :: nn=10
      integer(I4B) k, i, j, indext, jstack, l, r, swp, stacksize, im
      integer(I4B), allocatable :: istack(:)
!
!  Sorting algorithm:   Quick - Sort (adapted from Numerical Recipes)
!
!  sub-arrays of length up to nn are sorted by straight insertion
!
      stacksize=int(2*log(real(n))/log(real(2)))+2
      allocate (istack(stacksize))
      indx=(/(i, i=1,n)/) 
      jstack=0
      l=1
      r=n
      do
        if (r-l .lt. nn) then
!  sort by straight insertion
          do j=r-1,l,-1
            indext=indx(j)
            a=arr(indext)
            im=r+1
            do i=j+1,r
!              if (arr(indx(i)) .ge. a) then ! inverse order
              if (arr(indx(i)) .le. a) then
                im=i
                exit
              end if
              indx(i-1)=indx(i)
            end do
            indx(im-1)=indext
          end do
          if (jstack .eq. 0) then
            deallocate (istack)
            return
          end if
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+r)/2
!  sawp i and  l+1
          swp=indx(k)
          indx(k)=indx(l+1)
          indx(l+1)=swp
!          if (arr(indx(r)) .lt. arr(indx(l))) then ! inverse order
          if (arr(indx(r)) .gt. arr(indx(l))) then
            swp=indx(l)
            indx(l)=indx(r)
            indx(r)=swp
          end if
!          if (arr(indx(r)) .lt. arr(indx(l+1))) then ! inverse order
          if (arr(indx(r)) .gt. arr(indx(l+1))) then
            swp=indx(l+1)
            indx(l+1)=indx(r)
            indx(r)=swp
          end if
!          if (arr(indx(l+1)) .lt. arr(indx(l))) then ! inverse order
          if (arr(indx(l+1)) .gt. arr(indx(l))) then
            swp=indx(l)
            indx(l)=indx(l+1)
            indx(l+1)=swp
          end if
          i=l+1
          j=r
          indext=indx(l+1)
          a=arr(indext)
          do
            do
              i=i+1
!              if (arr(indx(i)) .ge. a) exit ! inverse order
              if (arr(indx(i)) .le. a) exit
            end do
            do
              j=j-1
!              if (arr(indx(j)) .le. a) exit ! inverse order
              if (arr(indx(j)) .ge. a) exit
            end do
            if (j .lt. i) exit
!  sawp i and  j
            swp=indx(i)
            indx(i)=indx(j)
            indx(j)=swp
          end do
          indx(l+1)=indx(j)
          indx(j)=indext
          jstack=jstack+2
          if (r-i+1 .ge. j-l) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          end if
        end if
      end do
      end subroutine qsortindex_sp
!
!
!
      subroutine qsortindex_i4(arr,indx,n)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: arr(:)
      integer (I4B) :: indx(:), n
      intent(in) :: arr, n
      intent(out) :: indx
!
!   Generate a sorted list of indices
!   Input:      arr      array to be sorted (key)
!               n        number of elements in arr
!   Output:     indx     sorted list of indices (largest element comes first)
!                        
      integer(I4B) :: a
      integer(I4B), parameter :: nn=10
      integer(I4B) k, i, j, indext, jstack, l, r, swp, stacksize, im
      integer(I4B), allocatable :: istack(:)
!
!  Sorting algorithm:   Quick - Sort (adapted from Numerical Recipes)
!
!  sub-arrays of length up to nn are sorted by straight insertion
!
      stacksize=int(2*log(real(n))/log(real(2)))+2
      allocate (istack(stacksize))
      indx=(/(i, i=1,n)/) 
      jstack=0
      l=1
      r=n
      do
        if (r-l .lt. nn) then
!  sort by straight insertion
          do j=r-1,l,-1
            indext=indx(j)
            a=arr(indext)
            im=r+1
            do i=j+1,r
!              if (arr(indx(i)) .ge. a) then ! inverse order
              if (arr(indx(i)) .le. a) then
                im=i
                exit
              end if
              indx(i-1)=indx(i)
            end do
            indx(im-1)=indext
          end do
          if (jstack .eq. 0) then
            deallocate (istack)
            return
          end if
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
        else
          k=(l+r)/2
!  sawp i and  l+1
          swp=indx(k)
          indx(k)=indx(l+1)
          indx(l+1)=swp
!          if (arr(indx(r)) .lt. arr(indx(l))) then ! inverse order
          if (arr(indx(r)) .gt. arr(indx(l))) then
            swp=indx(l)
            indx(l)=indx(r)
            indx(r)=swp
          end if
!          if (arr(indx(r)) .lt. arr(indx(l+1))) then ! inverse order
          if (arr(indx(r)) .gt. arr(indx(l+1))) then
            swp=indx(l+1)
            indx(l+1)=indx(r)
            indx(r)=swp
          end if
!          if (arr(indx(l+1)) .lt. arr(indx(l))) then ! inverse order
          if (arr(indx(l+1)) .gt. arr(indx(l))) then
            swp=indx(l)
            indx(l)=indx(l+1)
            indx(l+1)=swp
          end if
          i=l+1
          j=r
          indext=indx(l+1)
          a=arr(indext)
          do
            do
              i=i+1
!              if (arr(indx(i)) .ge. a) exit ! inverse order
              if (arr(indx(i)) .le. a) exit
            end do
            do
              j=j-1
!              if (arr(indx(j)) .le. a) exit ! inverse order
              if (arr(indx(j)) .ge. a) exit
            end do
            if (j .lt. i) exit
!  sawp i and  j
            swp=indx(i)
            indx(i)=indx(j)
            indx(j)=swp
          end do
          indx(l+1)=indx(j)
          indx(j)=indext
          jstack=jstack+2
          if (r-i+1 .ge. j-l) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
          else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
          end if
        end if
      end do
      end subroutine qsortindex_i4
!
!
!
      subroutine fsortn(felst,liste,n)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) n,liste(:)
      real (DP) felst(:)
      intent (in) :: felst,n
      intent (out):: liste
!
!   Erstellen einer sortierten Indexliste
!   Eingabe:    felst    zu sortierendes Feld (Schluessel)
!               n        Groesse von felst
!   Ausgabe:    liste    sortierte Indexliste (groesstes Element zuerst)
!                        
      integer (I4B) i,j,l,im,hold,m,iq,ir,jstack,stacksize
      integer (I4B), allocatable :: stack(:)
      real (DP) a
      real (SP) fm,fa,fc,fmi,fx
      parameter (m=8, fm=7875., fa=211., fc= 1663., fmi=1./fm)
!
!  Sortieralgorithmus:   Quick - Sort (Numerical Recipes)
!
!  Teilfolgen der Laenge m werden durch straight insertion sortiert
!
      stacksize=int(2*log(real(n))/log(real(2)))+1
      allocate (stack(stacksize))
      liste =(/(i, i=1,n)/) 
!
      jstack=0
      l=1
      ir=n
      fx=0._SP
      do
!  sortieren durch straight insertion
      if (ir-l .lt. m) then
        do j=l+1,ir
          hold=liste(j)
          a=felst(hold)
          do i=j-1,1,-1
            if (felst(liste(i)).ge.a) then
              im=i
              goto 12
            end if
            liste(i+1)=liste(i)
          end do
          im=0
12      liste(im+1)=hold
        end do
        if (jstack.eq.0) then
          deallocate (stack)
          return
        end if
        ir=stack(jstack)
        l=stack(jstack-1)
        jstack=jstack-2
      else
        i=l
        j=ir
!  Zufallszahl erzeugen
        fx=mod(fx*fa+fc,fm)
        iq=int(l+(ir-l+1)*(fx*fmi))
        hold=liste(iq)
        a=felst(hold)
        liste(iq)=liste(l)
        do
21        if (j.gt.0) then
            if (a.gt.felst(liste(j))) then
              j=j-1
              goto 21
            end if
          end if
          if (j.le.i) then
            liste(i)=hold
            goto 30
          end if
          liste(i)=liste(j)
          i=i+1
22        if (i.le.n) then
            if (a.lt.felst(liste(i))) then
              i=i+1
              goto 22
            end if
          end if
          if (j.le.i) then
            liste(j)=hold
            i=j
            goto 30
          end if
          liste(j)=liste(i)
          j=j-1
        end do
30      jstack=jstack+2
        if (jstack.gt.stacksize) then
          print*,' Error during sorting , size of stack insufficient'
        end if
        if (ir-i .ge. i-l) then
          stack(jstack)=ir
          stack(jstack-1)=i+1
          ir=i-1
        else
          stack(jstack)=i-1
          stack(jstack-1)=l
          l=i+1
        end if
      end if
      end do
      deallocate (stack)
      end subroutine fsortn