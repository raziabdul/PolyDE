      subroutine getes(es)
      use feminterface,       only: destroyarrptr, palloc
      use globalvariables3D,  only: nume, numv, vv, ve
      use femtypes
      implicit none
      type(ARRPTRI),     pointer :: es(:)
      intent(out)                :: es
!
!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 - 2015 Institute for Micro Systems Technology,
!                              Hamburg University of Technology.
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
!------------------------------------------------------------------------------
!
!    $Revision: 1.2 $
!    $Date: 2015/11/11 12:01:24 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
!    Generate the es edge to surface map
! Input:
!    --
! Output:
!    es         surfaces adjacent to an edge
!------------------------------------------------------------------------------
!  internal variables
      integer (I4B)              :: i, j, k, edge, surface
      integer (I4B),   parameter :: f2e(3,4) = reshape((/2,5,6,3,4,6,1,4,5,1,2,3/),(/3,4/))
      logical :: ok
!
      if (associated(es)) ok=destroyarrptr(es)
      allocate(es(nume))
      do edge = 1, nume
        nullify(es(edge)%d)
      end do

! f2e 3x4 array to return edges from face i: 
!
!  2  3  1  1
!  5  4  4  2
!  6  6  5  3

      do i = 1, numv
        do j = 1, 4
!  if vv is negative, it means that face j of volume i is surface element
          if ( vv(j,i) .lt. 0) then
            surface = -vv(j,i)
!  loop over the tree edges of that surface
            do k = 1, 3
!  get the edges of the surface element
              edge = ve( f2e(k,j),i)
              if (.not.associated(es(edge)%d)) then
                call palloc(es(edge)%d,2)
                es(edge)%d(1) = surface
                es(edge)%d(2) = 0
              else if (es(edge)%d(2) .eq. 0) then
                es(edge)%d(2) = surface
              else
                print*,'this should not happen'
              end if
            end do
          end if
        end do
      end do

      end subroutine getes
