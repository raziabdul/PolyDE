subroutine hp_alg_F8R(astep,v_res,ev,es,e_len)
      use globalvariables3D,   only: nnat, numv, numn, vp, vp_temp, vn
      use feminterface,        only: getsetting, palloc
      use feminterface3d,      only: element_mark3D, hmesh_refine3D, pmesh_refine3D
      use feminterface3d,      only: get_modres, getkpv, get_pconfirm
      use femtypes
      implicit none
      integer (I4B)          :: astep
      real     (DP)          :: v_res(:,:)
      real     (DP), pointer :: e_len(:)
      type(ARRPTRI), pointer :: ev(:), es(:)
      intent(in) :: astep
      intent(inout) :: v_res, ev,es,e_len
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
!    $Revision: 1.7 $
!    $Date: 2015/11/10 16:02:18 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'hp_alg_F8R'
!    Performs 'FIXED ORDER' hp-Adaptation Algorithm:
!       Elements, which belong to the top 5% of residual hight are subdivided
!    Adaptation Routine:
!       - ...
! Input:
!    --
!
! Output:
!    --
! 
! 
!------------------------------------------------------------------------------
!  Internal Variables:
      integer (I4B)               :: i_nat, num_marked, num_refined, adaptsteps, v
      integer (I4B),  allocatable :: mark_confirm(:), p_mark_left(:), max_vp(:), pconfirm(:,:)
      integer (I4B),  allocatable :: p_mark(:), h_mark(:)
      real     (DP)               :: err_bound_h, err_bound_p, real_fract_h, real_fract_p
      real     (DP),  allocatable :: mod_res(:,:), keyPointValues(:,:), sing_val(:,:)
      character (len=16)          :: hp_alg, marking_type_h, marking_type_p, h_adapt_crt
      logical                     :: adapt_vp
!__
!  Definitions:
      print "(A1)"        ,'|'
      print "(A27)"       ,'| ADAPTATION TYPE      : HP'
      print "(A1)"        ,'|'
      call getsetting('MAX_ERRORBOUND_H',err_bound_h)
      call getsetting('MAX_ERRORBOUND_P',err_bound_p)
      call getsetting('MARKING_TYPE_H',marking_type_h)
      call getsetting('MARKING_TYPE_P',marking_type_p)
      call getsetting('H_ADAPT_CRIT',h_adapt_crt)
      call getsetting('ADAPT_STEPS',adaptsteps)
      call getsetting('HP_ALGORITHM',hp_alg)
!__
! Allocate Memory:
      allocate(mark_confirm(numv))
      allocate(h_mark(numv))

!__
! Checks and preparations:
      if (hp_alg == 'TOP5') then
        print "(A39)"       ,'| ALGORITHM  CATEGORY  : FIXED-FRACTION'
        print "(A29)"       ,'| ADAPTATION ALGORITHM : TOP5'
        print "(A1)"        ,'|'
        adapt_vp = .true.
      else
        print "(A42)"       ,'| ALGORITHM  CATEGORY  : FIXED-ORDER (F8R)'
        print "(A28)"       ,'| ADAPTATION ALGORITHM : F8R'
        print "(A1)"        ,'|'
        adapt_vp = .false.
      end if
      mark_confirm(:) = 1
      h_mark(:)       = 0
!____
! Start:
! 1) h-Mark elements / DOFs for each Nature:
      print "(A3)" ,'|--'
      real_fract_h = 0._DP
      select case (h_adapt_crt)
        case('SING_VAL')
          allocate(sing_val(numv,nnat))
          allocate(mod_res(numv,nnat))
          allocate(keyPointValues(numn,nnat))
          do i_nat = 1, nnat
!-        a) Either Modres should be calculated with 'vp' or with 'vp_temp'
            if (adapt_vp) then
              call get_modres(v_res(:,i_nat), mod_res(:,i_nat),vp(:,i_nat))
            else
              call get_modres(v_res(:,i_nat), mod_res(:,i_nat),vp_temp(:,i_nat))
            end if
!-        b) Get singularity values for nodes, map onto elements
            call getkpv(mod_res(:,i_nat),keyPointValues(:,i_nat))
            do v = 1, numv
              sing_val(v,i_nat) = sum(keyPointValues(vn(:,v),i_nat))/4.
            end do

!-        c) Mark elements according their sing_val
            if (i_nat.eq.3) then
              !call element_mark3D(sing_val(:,i_nat),mark_confirm,h_mark,1,err_bound_h,marking_type_h,num_marked,produce_plot='SORTED')
              call element_mark3D(sing_val(:,i_nat),mark_confirm,h_mark,1,err_bound_h,marking_type_h,num_marked)
            else
              call element_mark3D(sing_val(:,i_nat),mark_confirm,h_mark,1,err_bound_h,marking_type_h,num_marked)
            end if
          end do
          deallocate(keyPointValues)
          deallocate(mod_res)

        case('RESIDUAL')
          do i_nat = 1, nnat
            call element_mark3D(v_res(:,i_nat),mark_confirm,h_mark,1,err_bound_h,marking_type_h,num_marked)
          end do

      end select
!__
! 2) Compute 
      allocate(max_vp(nnat))
      do i_nat = 1, nnat
        max_vp(i_nat) = maxval(vp(:,i_nat))
      end do
      if (h_adapt_crt.eq.'SING_VAL') then
        allocate(pconfirm(numv,nnat))
        do i_nat = 1, nnat
          call get_pconfirm(sing_val(:,i_nat),max_vp(i_nat),pconfirm(:,i_nat))
        end do
        deallocate(sing_val)
      end if
      deallocate(max_vp)
!__
! 3) Perform p-refinement on vp_temp, h-refinement is performed in a normal way.
      allocate(p_mark_left(numv))
      p_mark_left(:)  = 0
      allocate(p_mark(numv))
      p_mark(:)       = 0
      real_fract_p = 0._DP
      do v = 1, numv
        if (h_mark(v) .gt. 0) then 
          mark_confirm(v) = 0
          p_mark(v) = p_mark(v) - 1
        end if
      end do
      do i_nat = 1, nnat
!-    Top Error Bound Refinement:
        print "(A40)"         ,'|--------------------------------------|'
        print "(A18,I1)"      ,'| Nature:         ',i_nat
        print "(A1)"          ,'|'
        call element_mark3D(v_res(:,i_nat),mark_confirm,p_mark,1,err_bound_p,marking_type_p,num_marked)
        if (h_adapt_crt .eq. 'SING_VAL') then
          call pmesh_refine3D(p_mark,i_nat,num_refined,ref_vp=adapt_vp,pconfirm=pconfirm)
        else
          call pmesh_refine3D(p_mark,i_nat,num_refined,ref_vp=adapt_vp)
        end if
        print "(A18,I8,A9)"   ,'| Marking:        ',num_marked,' elements'
        real_fract_p = real(num_marked,DP)*100/real(numv,DP)
        print "(A20,F6.3,A2)" ,'|                   ',real_fract_p,' %'
        print "(A18,I8,A9)"   ,'| p-Refinement of ',num_refined,' elements'
        real_fract_p = real(num_refined,DP)*100/real(numv,DP)
        print "(A20,F6.3,A2)" ,'|                   ',real_fract_p,' %'
        print "(A1)"          ,'|'
!__
!  h-Refinement of elements, which have been p-Marked, but not p-Refined
        if (any(p_mark(:) .gt. 0)) then
          do v = 1, numv
            if (p_mark(v) .gt. 0) then
              p_mark_left(v) = 1
              p_mark(v) = 0
            end if
          end do
        end if
      end do
      h_mark(:) = h_mark(:) + p_mark_left(:)
      deallocate(p_mark_left)
      deallocate(mark_confirm)
      deallocate(p_mark)
      if (allocated(pconfirm)) deallocate(pconfirm)
!__
! 4) h-refine
      do v = 1, numv
        if (h_mark(v) .ne. 0) real_fract_h = real_fract_h + 1
      end do
      call hmesh_refine3D(h_mark,es,ev,e_len,num_refined)
      print "(A18,I9,A13)" ,'| h-Refining:     ',num_refined,' new elements'
      print "(A18,I9,A16)" ,'|                 ',int(real_fract_h),' marked elements'
      real_fract_h = real(real_fract_h,DP)*100./real(numv-num_refined,DP)
      print "(A20,F6.3,A2)",'|                   ',real_fract_h,' %'
      print "(A79)"        ,' _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  .  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ '
      print "(A1)"         ,'|'
!__
! 5) In the last adaptation step, set vp = vp_temp
      if ((.not.adapt_vp).and.(astep .eq. adaptsteps)) then
        vp(:,:) = vp_temp(:,:)
      end if

!__
!  Release Memory:
      deallocate(h_mark)
!__
!
!_End.
end subroutine hp_alg_F8R