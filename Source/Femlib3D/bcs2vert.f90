      subroutine bcs2vert(ok)
      use feminterface, only: getsetting
      use femtypes
      use globalvariables3D, only : sbc, numn, nnat, sn, nums, sfvertbc, bctype
      implicit none
      logical ok
      intent (out) :: ok
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
!    $Revision: 1.7 $
!    $Date: 2015/07/16 12:52:05 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!  INFO:
!  That Routine updates the entries of the vector sfvertbcn(numn)
!  which holds the boundary condition information of vertices
!
!------------------------------------------------------------------------------
!  Output:
!
!      ok        =.false. if an error occurs
!
!  local variables
      integer (I4B) i, inat, vloop



!  read domain and volume element - node map ------------------------------------ changed
!      do i=1,nums
!        read (unitid,*) sbc(i),sn(:,i)
!      end do
      allocate( sfvertbc(numn,nnat) )
      sfvertbc(:,:)= -1
      do i=1,nums
        do inat=1,nnat
          do vloop = 1, 3
            if (sfvertbc(sn(vloop,i),inat) .ge. 0 ) then
              if ( bctype(sbc(i),inat) .lt. bctype(sfvertbc(sn(vloop,i),inat),inat) ) then
                sfvertbc(sn(vloop,i),inat) = sbc(i)
              else
                cycle
              end if
            else
              sfvertbc(sn(vloop,i),inat) = sbc(i)
            end if
          end do
        end do
      end do
      ok = .true.
      return
      end subroutine bcs2vert