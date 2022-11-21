subroutine hp_adapt3D(astep,v_res,ev,es,e_len)
      use feminterface,        only: getsetting
      use feminterface3d,      only: hp_alg_F8R

      use femtypes
      implicit none
      
      integer (I4B)          :: astep
      real     (DP)          :: v_res(:,:)
      real     (DP), pointer :: e_len(:)
      type(ARRPTRI), pointer :: ev(:), es(:)
      intent(in) :: astep
      intent(out) :: v_res
      intent(inout) :: ev, es, e_len
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
!    $Revision: 1.5 $
!    $Date: 2015/11/10 14:23:34 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'hp_adapt3D'
!    Performs 
!    Adaptation Routine:
!       - Choose Adaptation Algorithm
!       - Call adaptation Algorithm
! Input:
!    --
!
! Output:
!    --
! 
! 
!------------------------------------------------------------------------------
!  internal variables:
      character (len=20)              :: hp_alg
!__
!  Definitions:
      print "(A1)"        ,'|'
      print "(A79)"       ,' ___________________________HP - A D A P T A T I O N____________________________ '
      print "(A1)"        ,' '
      print "(A44,I2,A28)",'|                                 S T E P -< ',astep,' >-                           |'
!-
      call getsetting('HP_ALGORITHM',hp_alg)
!____
!  Start:
      if (hp_alg == 'TOP5') hp_alg = 'F8R'
      select case (hp_alg)

        case ('F8R')
          call hp_alg_F8R(astep,v_res,ev,es,e_len)

        case ('KEYPOINT')


        case ('M07')


        case ('M04')


        case ('REF_HISTORY')


        case default


      end select
!__
!  Release Memory:

!
!_End.
end subroutine hp_adapt3D
