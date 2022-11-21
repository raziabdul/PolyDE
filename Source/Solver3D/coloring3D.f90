      subroutine coloringsc3D(indx,offset)
      use feminterface, only: realloc
      use femtypes
      use globalvariables3d, only: numv, vn, numn
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
!    $Date: 2015/11/10 14:03:58 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
!  coloring of elements such that elements which share a DOF (i.e. a vertex for scalar nodal elements)
!  are in different color sets
!  Output:  
!            indx    array of element indices, sorted by color
!            offset  offest(i) to offest(i+1)-1 belong to color i
!
!  local variables:
      integer (I4B) :: i, j, actualcolor, numcolors
      logical, allocatable :: vc_color(:,:)
      integer (I4B), allocatable :: color_dist(:), ec(:)
!
      print*,'Entered coloring...'
!
! assumed number of colors; 45-55 usually is a good guess for 3D mesh with nodal elements
      numcolors=50
! element color array
      allocate(ec(numv))
!
! array with vertex colors
      allocate(vc_color(numcolors,numn))
      vc_color = .false.
!
      ec(1) = 1
      vc_color(1,vn(:,1)) = .true.
! all elements will be colored
      do i = 2, numv
! pick lowest color
        actualcolor = numcolors+1
        do j = 1,numcolors
          if (all(vc_color(j,vn(:,i)) .eq. .false.)) then
            actualcolor = j
            exit
          end if
        end do
!
        if (actualcolor .gt. numcolors) then
! reallocate array
          numcolors = numcolors*2
          call realloc(vc_color,numcolors,numn)
          vc_color(actualcolor:numcolors,:) = .false.
        end if
        ec(i) = actualcolor
        vc_color(actualcolor,vn(:,i)) = .true.
!
      end do
!
      deallocate(vc_color)

      allocate(color_dist(numcolors))
!  count number of elements in each of the colors
      color_dist = 0
      do i= 1,numv
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
      offset(numcolors+1) = numv+1
!
! some data about the coloring
      print*,'Colored elements  :', sum(color_dist)
      print*,'Colors used       :', numcolors
!      print*,'Color distribution:', color_dist(1:numcolors)
!
! sort ec to index array
      allocate(indx(numv))
      color_dist(1:numcolors)=offset(1:numcolors)
      do i=1,numv
        indx(color_dist(ec(i))) = i
        color_dist(ec(i)) = color_dist(ec(i))+1
      end do
!
      deallocate(color_dist, ec)
! pass  index, offset to caller
      return
      end subroutine coloringsc3D
