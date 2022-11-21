subroutine prepare_h_adapt3D(ev,es,e_len)
      use feminterface,       only: palloc, zeit
      use feminterface3d,     only: reversemap, getes
      use globalvariables3D,  only: ve, nume, nod, en
      use femtypes
      implicit none
      real     (DP),     pointer :: e_len(:)
      type(ARRPTRI),     pointer :: ev(:), es(:)
      intent(out)                :: ev, es, e_len
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
! This Routine: 'prepare_h_adapt3D'
!    Performs a setup of FEM for h-Adaptation procedure:
! Input:
!    --
! Output:
!    ev         volumes adjacent to an edge
!    es         surfaces adjacent to an edge
!    e_len      length of edges
!------------------------------------------------------------------------------
!  internal variables
      integer (I4B)              :: nmin, nmax, edge
!
!  Get ve
      call reversemap(.true.,ve,ev,nmin,nmax)
!  calculate edge lengths
      call palloc(e_len,nume)
      do edge = 1, nume
        e_len(edge) = norm2( nod(:,en(1,edge))-nod(:,en(2,edge)) )
      end do
!  Get es
      call getes(es)

end subroutine prepare_h_adapt3D



subroutine prepare_p_adapt3D
      use feminterface,       only: getsetting
      use globalvariables3D,  only: vp
      use femtypes
      implicit none
!------------------------------------------------------------------------------
!
! This Routine: 'prepare_p_adapt3D'
!    Performs a setup of FEM for p-Adaptation procedure:
! Input:
!    --
! Output:
!    --
!------------------------------------------------------------------------------
!  internal variables
!
!__

end subroutine prepare_p_adapt3D



subroutine prepare_hp_adapt3D(ev,es,e_len)
      use feminterface,       only: palloc, getsetting, zeit
      use feminterface3d,     only: reversemap, getes
      use globalvariables3D,  only: ve, nume, nod, en
      use globalvariables3D,  only: nnat, numv, vp, vp_temp
      use femtypes
      implicit none
      real     (DP),     pointer :: e_len(:)
      type(ARRPTRI),     pointer :: ev(:), es(:)
      intent(out)                :: ev, es, e_len
!------------------------------------------------------------------------------
!
! This Routine: 'prepare_hp_adapt3D'
!    Performs a setup of FEM for hp-Adaptation procedure:
! Input:
!    --
! Output:
!    ev         volumes adjacent to an edge
!    es         surfaces adjacent to an edge
!    e_len      length of edges
!------------------------------------------------------------------------------
!  internal variables
      integer (I4B)              :: nmin, nmax, edge, polyorder, hpfo_order
      character (len=16)         :: hp_type
!
!!  Set Polyorder of Solution:
      call getsetting('HP_ALGORITHM',hp_type)
      call getsetting('POLYORDER',polyorder)
      if (hp_type.eq.'F8R') then
        call palloc(vp_temp,numv,nnat)
        call getsetting('HPF8R_ORDER',hpfo_order)
        vp(:,:)      = hpfo_order
        vp_temp(:,:) = polyorder
      end if
!
!  Get ve
      call reversemap(.true.,ve,ev,nmin,nmax)
!  calculate edge lengths
      call palloc(e_len,nume)
      do edge = 1, nume
        e_len(edge) = norm2( nod(:,en(1,edge))-nod(:,en(2,edge)) )
      end do
!  Get es
      call getes(es)

end subroutine prepare_hp_adapt3D