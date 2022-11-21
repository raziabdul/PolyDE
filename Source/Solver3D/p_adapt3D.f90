subroutine p_adapt3D(astep,v_res)
      use globalvariables3D,   only: nnat, numv
      use feminterface,        only: getsetting, palloc
      use feminterface3d,      only: element_mark3D, pmesh_refine3D
      use femtypes
      implicit none
      integer (I4B)          :: astep
      real     (DP)          :: v_res(:,:)
      intent (in) :: astep
      intent (out) :: v_res
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
!    $Date: 2015/11/10 15:23:38 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'p_adapt3D'
!    Performs 
!    Folowing parts of the adaptation routine:
!       - ...
!       - Mark elements / DOFs
!       - Apply Refinement Algorithms
!       - ...
! Input:
!    --
!
! Output:
!    --
! 
!
!------------------------------------------------------------------------------
!  internal variables
      integer (I4B)               :: i_nat, num_marked, num_refined
      integer (I4B),  allocatable :: p_mark(:), mark_confirm(:)
      real     (DP)               :: err_bound_p, real_fract
      character (len=16)          :: marking_type_p
!__
!  Allocate Memory:
      allocate(p_mark(numv))
      allocate(mark_confirm(numv))
!__
!  Definitions
      print "(A1)"        ,'|'
      print "(A79)"       ,' ____________________________P - A D A P T A T I O N____________________________ '
      print "(A44,I2,A28)",'|                                 S T E P -< ',astep,' >-                           |'
      call getsetting('MAX_ERRORBOUND_P',err_bound_p)
      call getsetting('MARKING_TYPE_P',marking_type_p)
!-
      real_fract = 0._DP
      p_mark(:)  = 0
      mark_confirm(:) = 1
!____
!  Start:
!  For each Nature:
!  1) Mark elements / DOFs
        do i_nat = 1, nnat
!-    Top Error Bound Refinement:
          call element_mark3D(v_res(:,i_nat),mark_confirm,p_mark,1,err_bound_p,marking_type_p,num_marked)
          real_fract = real(num_marked,DP)*100./real(numv,DP)
          print "(A3)"          ,'|--'
          print "(A18,I1)"      ,'| Nature:         ',i_nat
          print "(A18,I7,A10)"  ,'| Marking:        ',num_marked,' elements.'
          print "(A20,F5.2,A3)" ,'|                   ',real_fract,' %.'
!-
          print "(A1)"          ,'|'
          call pmesh_refine3D(p_mark,i_nat,num_refined,ref_vp=.true.)
          print "(A18,I7,A10)"  ,'| Refinement of   ',num_refined,' elements.'
          print "(A1)"          ,'|'
!-
          p_mark(:)  = 0
        end do
      print "(A79)"        ,' _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  .  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ '
      print "(A1)"         ,'|'
!__
!  Release Memory:
      deallocate(p_mark)
      deallocate(mark_confirm)
!
!_End.
end subroutine p_adapt3D