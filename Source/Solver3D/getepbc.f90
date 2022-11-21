      subroutine getepbc
      use feminterface, only: destroyarrptr
      use feminterface3D, only: getes
      use femtypes
      use globalvariables3D, only : sbc, edgebc, nnat, nume, bctype
      implicit none
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!    $Revision: 1.10 $
!    $Date: 2014/08/22 11:25:13 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  returns the corresponding BC index for edges on surfaces 
!
!  Output:   edgebc      BC info (global)of the edge for each nature
!
!
!  Input:    sbc         BC of surface elements (global)
!            ve          volume to edge map     (global)
!            vv          volume neighbors info  (global)
!            nume        number of edges
!            numv        number of volume elements (global)
!            nums        number of surface elements (global)
!
! local variables 

      type (ARRPTRI), pointer :: edgesf(:)=>null()
      integer (I4B) :: i, inat
      logical :: ok
!  fetch surfeces at edges
      call getes(edgesf)

      if (associated(edgebc)) deallocate(edgebc)
      allocate( edgebc(nume,nnat) )
      edgebc = -1
! compare bc of each faces-sharing edge
! check for all edges that are connected to surface elements
      do i = 1, nume
! getes returns null d for edges not connected to 
! the surfaces in the first place, so skip this edge
        if ( associated( edgesf(i)%d ) ) then
! essential BC is always lower then general Neumann
          do inat=1,nnat
            if ( bctype(sbc( edgesf(i)%d(1)),inat ) .lt. bctype(sbc( edgesf(i)%d(2)),inat ) ) then
              edgebc(i,inat) = sbc( edgesf(i)%d(1))
! if not, either equal or greater, so get second
            else
              edgebc(i,inat) = sbc( edgesf(i)%d(2))
            end if
          end do
        end if
      end do
      ok = destroyarrptr(edgesf)
      return
      end subroutine getepbc
