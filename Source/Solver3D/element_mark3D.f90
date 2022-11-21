subroutine element_mark3D(crt,mark_confirm,element_mark,mark_depth,error_bound,marking_type,num_marked,produce_plot)
      use globalvariables3D,    only: vv
      use feminterface,         only: qsortindex, print_error, getsetting
      use feminterface3d,       only: plot_marking
      use femtypes
      implicit none
      integer (I4B)               :: mark_depth, num_marked
      integer (I4B)               :: element_mark(:)
      integer (I4B)               :: mark_confirm(:)
      real     (DP)               :: error_bound
      real     (DP)               :: crt(:)
      character (len=*), optional :: produce_plot
      character (len=16)          :: marking_type
      intent   (in)               :: mark_depth, error_bound
      intent(inout)               :: element_mark, crt
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
!    $ Revision: 0.1 $
!    $ Date: 2015/01/01 11:11:11 $
!    $ Author: juryanatzki $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'element_mark3D'
! element_mark elements in the mesh for adaptation
!
!-------------------------------------------------------------------------------
! Input:
!           crt (numv, i_nat)   Criterion for element selection (e.g. error estimates, element-shapes-measure, ..)
!         element_mark (numv)   Element markings vector.
!                             - 0: element is not marked
!                             - mark_depth: element is marked for 'mark_depth' refinements (mark_depth > 0)
!                  mark_depth   Number of refinements to be used during marking
!                         acc   Logical Value indicating the type of barrier.
!                             - If True, criteria are accumulated up to the barrier
!                             - If False, all elements below the barrier are excluded
!                 error_bound   Criterion Barrier. Relative Value between 0,xx -> 1. Can not be zero (leads to error)
!                 fixed_fract   Fixed Fraction of amount of elements to mark
!
! Output:
!         element_mark (numv)   Element markings vector.
!                               'element_mark(:)' might already contain a marking pattern. New element_mark will be added onto the old ones.
!                               -  0: element is not marked
!                               -  mark_depth: element is marked for 'mark_depth' refinements (mark_depth > 0)
!                  num_marked   Number of marked elements.
!
!-------------------------------------------------------------------------------
! local variables:
!
      integer (I4B)              :: fixed_size, crt_size, min_mark, i, num_nb, crt_breakloc
      real     (DP)              :: crt_tot, crt_acc, crt_error_bound, fixed_fract
      real     (DP), allocatable :: steep_val(:)
      integer (I4B), allocatable :: mk_idx(:)
      logical                    :: refine_all
!__
! Definitions:
      num_marked      = 0
      crt_acc         = 0
      min_mark        = 1
      refine_all      = .false.
      crt             = abs(crt)
      crt_tot         = sum(crt)
      crt_size        = size(crt,1)
      crt_error_bound = crt_tot*error_bound
      num_nb          = size(vv,1)
!__
! Checks and preparations:
      if (error_bound.gt.1.) then
        call print_error(3,'Can not perform the marking of elements. Given barrier is higher than 100%.','Refining 100% ...')
        refine_all = .true.
      end if
      if (error_bound.eq.0.) then
        call print_error(5,'Can not perform the marking of elements. Given barrier is equal to 0%.')
      end if
      if (error_bound.eq.(1.)) then
        refine_all = .true.
      end if
!-
      if (marking_type .eq. 'FIXED_FRACTION') then

        call getsetting('ADAPT_FRACTION',fixed_fract)

        if (fixed_fract.gt.1.) then
          call print_error(3,'Can not perform the marking of elements.','Given fraction of refinement elements is higher than 100%. Refining 100% ...')
          refine_all = .true.
        end if
        if (fixed_fract.eq.0.) then
          call print_error(5,'Can not perform the marking of elements. Given fraction of refinement elements is equal to 0%.')
        end if
        if ( fixed_fract.eq.(1.) ) then
          refine_all = .true.
        end if
      end if
!- If all elements should be refined:
      if (refine_all) then
        element_mark(:) = element_mark(:) + mark_depth
        num_marked = crt_size
        return
      end if
!__
!
! Allocate Memory
      allocate(mk_idx(crt_size))
!__
! Start:
! 1) Sort elements according to criteria. Elements with highest criterion come first.
      call qsortindex(crt,mk_idx,crt_size)
!__
! 2) Plot Marking Statistics if required
      if (present(produce_plot)) then
        select case (produce_plot)
          case ('SORTED')
            call plot_marking(crt(mk_idx(:)))

          case ('UNSORTED')
            call plot_marking(crt(:))

        end select
      end if
!__
! 2) Start algorithm
      select case (marking_type)
      case ('FIXED_FRACTION')

        fixed_size = nint(crt_size*fixed_fract)
        i = 1
        do while ((i <= fixed_size).and.(i <= crt_size))
          element_mark(mk_idx(i)) = element_mark(mk_idx(i)) + mark_depth
          i = i + 1
        end do

      case('ACCUMULATIVE')

        crt_tot = 0.
        do i = 1, crt_size
          if (mark_confirm(mk_idx(i)).ne.0) crt_tot = crt_tot + crt(mk_idx(i))
        end do
        crt_error_bound = crt_tot*error_bound

        i = 1
        do while ( ((crt_acc <= crt_error_bound) .or. (i < min_mark)) .and. (i <= crt_size))
          if (mark_confirm(mk_idx(i)).ne.0) then
            crt_acc = crt_acc + crt(mk_idx(i))
            element_mark(mk_idx(i)) = element_mark(mk_idx(i)) + mark_depth
          end if
          i = i + 1
        end do

      case('NONACCUMULATIVE')
        crt_error_bound = error_bound*maxval(crt)
        i = 1
        do while ((crt(mk_idx(i)) > crt_error_bound))
          element_mark(mk_idx(i)) = element_mark(mk_idx(i)) + mark_depth
          i = i +1
        end do

      case('STEEPNESS_NACC')
        allocate(steep_val(crt_size-1))

        crt_breakloc = 1
        do i = 1, size(steep_val,1)
          steep_val(i) = crt(mk_idx(i)) - crt(mk_idx(i+1))
          if (steep_val(i).gt.steep_val(crt_breakloc)) crt_breakloc = i
        end do
        element_mark(mk_idx(1:crt_breakloc)) = element_mark(mk_idx(1:crt_breakloc)) + mark_depth
      end select
!__
! 3) Count the number of marked elements:
      do i = 1, crt_size
        if (element_mark(i).gt.0) num_marked = num_marked + 1
      end do
!__
!
! Release Memory:
      deallocate(mk_idx)
      if (allocated(steep_val)) deallocate(steep_val)
!
! End.
!___
 return
end subroutine element_mark3D
