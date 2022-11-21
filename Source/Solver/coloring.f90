      subroutine coloring(indx,offset)
      use feminterface, only: realloc
      use femtypes
      use globalvariables, only: n, e, p
      implicit none
      integer (I4B), allocatable :: indx(:), offset(:)
      intent (out) :: indx, offset
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
!    $Revision: 1.6 $
!    $Date: 2014/02/18 12:15:25 $
!    $Author: juryanatzki $
!
!-------------------------------------------------------------------------------
!
!  coloring of elements such that elements which share a DOF (i.e. vertex) are in different color sets
!  Output:  
!            indx    array of element indices, sorted by color
!            offset  offest(i) to offest(i+1)-1 belong to color i
!
!  local variables:
      integer (I4B) :: numcolors
      logical, allocatable :: vc_color(:,:), usedcol(:)
      integer (I4B), allocatable :: color_dist(:), ec(:)
      integer (I4B) :: i, j, actualcolor
!
      print*,'Entered coloring...'
!
! assumed number of colors; 12-14 usually is a good guess for 2D mesh with nodal elements
      numcolors=14
! element color array
      allocate(ec(n))
      allocate(usedcol(numcolors))
!
! array with vertex colors
      allocate(vc_color(numcolors,p))
      vc_color = .false.
!
      ec(1) = 1
      vc_color(1,e(:,1)) = .true.
! all elements will be colored
      do i = 2,n
        usedcol(:) = vc_color(:,e(1,i)) .or. vc_color(:,e(2,i)) .or. vc_color(:,e(3,i))
! pick lowest color
        actualcolor = numcolors+1
        do j = 1,numcolors
          if (usedcol(j) .eq. .false.) then
            actualcolor = j
            exit
          end if
        end do
!
        if (actualcolor .gt. numcolors) then
! reallocate arrays
          numcolors = numcolors+1
          call realloc(vc_color,numcolors,p)
          vc_color(numcolors,:) = .false.
          call realloc(usedcol,numcolors)
        end if
        ec(i) = actualcolor
        vc_color(actualcolor,e(:,i)) = .true.
!
      end do

      deallocate(vc_color)
      deallocate(usedcol)
!
      allocate(color_dist(numcolors))
!  count number of elements in each of the colors
      color_dist = 0
      do i= 1,n
        color_dist(ec(i)) = color_dist(ec(i)) + 1
      end do
!
!  correct numcolors
      do i=numcolors,1,-1
        if(color_dist(i) .gt. 0) then
          numcolors=i
          exit
        end if
      end do
!
! assign offset 
      allocate(offset(numcolors+1))
      offset(1) = 1
      do i=2,numcolors
        offset(i) = offset(i-1) + color_dist(i-1)
      end do
      offset(numcolors+1) = n+1
!
! some data about the coloring
      print*,'Colored elements  :', sum(color_dist)
      print*,'Colors used       :', numcolors
      print*,'Color distribution:', color_dist(1:numcolors)
!
! sort ec to index array
      allocate(indx(n))
      color_dist(1:numcolors) = offset(1:numcolors)
      do i=1,n
        indx(color_dist(ec(i))) = i
        color_dist(ec(i)) = color_dist(ec(i))+1
      end do
!
      deallocate(color_dist, ec)
! pass  index, offset to caller

      return
      end subroutine coloring