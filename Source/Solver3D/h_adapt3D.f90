subroutine h_adapt3D(astep,v_res,ev,es,e_len)
      use globalvariables3D,   only: nnat, vp, vn, numv, numn
      use feminterface,        only: getsetting, palloc
      use feminterface3d,      only: element_mark3D, hmesh_refine3D, get_modres, getkpv
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
!    $Revision: 1.6 $
!    $Date: 2015/11/10 14:28:04 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'h_adapt3D'
!    Performs 
!    Adaptation Routine:
!       - Solve the problem
!       - Estimate the error
!       - Mark elements / DOFs
!       - Apply Refinement Algorithms
!       - Set Up adapted DOFs h- and p- Grids
! Input:
!    --
!
! Output:
!    --
! 
! 
!------------------------------------------------------------------------------
!  internal variables
      integer (I4B)               :: i_nat, num_marked, num_refined, v
      integer (I4B),  allocatable :: h_mark(:), mark_confirm(:)
      real     (DP)               :: real_fract, err_bound_h
      real     (DP),  allocatable :: mod_res(:,:), keyPointValues(:,:), sing_val(:,:)
      character (len=16)          :: marking_type_h, h_adapt_crt
!__
!  Allocate Memory:      
      allocate(mod_res(numv,nnat))
      allocate(mark_confirm(numv))
      allocate(h_mark(numv))
!__
!  Definitions
      print "(A1)"        ,'|'
      print "(A79)"       ,' ____________________________H - A D A P T A T I O N____________________________ '
      print "(A44,I2,A28)",'|                                 S T E P -< ',astep,' >-                           |'
      call getsetting('MARKING_TYPE_H',marking_type_h)
      call getsetting('MAX_ERRORBOUND_H',err_bound_h)
      call getsetting('H_ADAPT_CRIT',h_adapt_crt)
!-
      h_mark(:) = 0
      num_marked = 0
      mark_confirm(:) = 1
!-
      real_fract = 0._DP
!____
!  Start:
!  For each Nature:
!  1) Mark elements / DOFs
      select case (h_adapt_crt)

        case('SING_VAL')
          allocate(keyPointValues(numn,nnat))
          allocate(sing_val(numv,nnat))
          do i_nat = 1, nnat
            call get_modres(v_res(:,i_nat), mod_res(:,i_nat),vp(:,i_nat))     
            call getkpv(mod_res(:,i_nat),keyPointValues(:,i_nat))
            do v = 1, numv
              sing_val(v,i_nat) = sum(keyPointValues(vn(:,v),i_nat))/4.
            end do
            call element_mark3D(sing_val(:,i_nat),mark_confirm,h_mark,1,err_bound_h,marking_type_h,num_marked)
          end do
          deallocate(sing_val)
          deallocate(keyPointValues)

        case('RESIDUAL')
          do i_nat = 1, nnat
            call element_mark3D(v_res(:,i_nat),mark_confirm,h_mark,1,err_bound_h,marking_type_h,num_marked)
          end do

      end select

      real_fract = (real(num_marked,DP)/real(numv,DP))*100._DP
      call hmesh_refine3D(h_mark,es,ev,e_len,num_refined)
              print "(A3)" ,'|--'
      print "(A18,I7,A10)" ,'| h-REFINING:     ',num_refined,' elements.'
      print "(A20,F5.2,A3)",'|                   ',real_fract,' %.'
      print "(A79)"        ,' _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  .  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ '
      print "(A1)"         ,'|'
!__
!  Release Memory:
      deallocate(mod_res)
      deallocate(h_mark)
      deallocate(mark_confirm)
!
!_End.
end subroutine h_adapt3D