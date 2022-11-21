      subroutine qsort_dp(arr,n)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) arr(:)
      integer (I4B) n
      intent(in) :: n
      intent(inout) :: arr
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
!   Sort an array using Quick - Sort (adapted from Numerical Recipes)
!
!   Input:        n        number of elements in arr
!   In-/ output:  arr      array to be sorted
!                          sorted on output (largest element comes first)
!                        
      real(DP) a,swp
      integer(I4B), parameter :: nn=10
      integer(I4B) k, i, j, jstack, l, r, stacksize, im
      integer(I4B), allocatable :: istack(:)
!
!  sub-arrays of length up to nn are sorted by straight insertion
!
      stacksize=int(2*log(real(n))/log(real(2)))+2
      allocate (istack(stacksize))
      jstack=0
      l=1
      r=n
      do
        if (r-l .lt. nn) then
!  sort by straight insertion
          do j=r-1,l,-1
            a=arr(j)
            im=r+1
            do i=j+1,r
!              if (arr(i) .ge. a) then ! inverse order
              if (arr(i) .le. a) then
                im=i
                exit
              end if
              arr(i-1)=arr(i)
            end do
            arr(im-1)=a
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
          swp=arr(k)
          arr(k)=arr(l+1)
          arr(l+1)=swp
!          if (arr(r) .lt. arr(l)) then ! inverse order
          if (arr(r) .gt. arr(l)) then
            swp=arr(l)
            arr(l)=arr(r)
            arr(r)=swp
          end if
!          if (arr(r) .lt. arr(l+1)) then ! inverse order
          if (arr(r) .gt. arr(l+1)) then
            swp=arr(l+1)
            arr(l+1)=arr(r)
            arr(r)=swp
          end if
!          if (arr(l+1) .lt. arr(l)) then ! inverse order
          if (arr(l+1) .gt. arr(l)) then
            swp=arr(l)
            arr(l)=arr(l+1)
            arr(l+1)=swp
          end if
          i=l+1
          j=r
          a=arr(l+1)
          do
            do
              i=i+1
!              if (arr(i) .ge. a) exit ! inverse order
              if (arr(i) .le. a) exit
            end do
            do
              j=j-1
!              if (arr(j) .le. a) exit ! inverse order
              if (arr(j) .ge. a) exit
            end do
            if (j .lt. i) exit
!  sawp i and  j
            swp=arr(i)
            arr(i)=arr(j)
            arr(j)=swp
          end do
          arr(l+1)=arr(j)
          arr(j)=a
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
      end subroutine qsort_dp
   

