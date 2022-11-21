      function destroyarrptr_i(arrptr)
      use feminterface, only:
      use femtypes
      implicit none
      type(ARRPTRI), pointer :: arrptr(:)
      logical :: destroyarrptr_i
!
!------------------------------------------------------------------------------
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
!    $Revision: 1.10 $
!    $Date: 2015/11/03 16:01:45 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  destroy an array of pointers (derived data type) which contains integer data
!
!------------------------------------------------------------------------------
!  Input:
!     arrptr    array of pointers (integer data)
!
!  Internal variables:
      integer (I4B) :: i, istat
!
      destroyarrptr_i = .true.
!  Deallocate individual entries
      do i = lbound(arrptr,1), ubound(arrptr,1)
        deallocate(arrptr(i)%d, STAT=istat)
        if (istat .ne. 0) destroyarrptr_i = .false.
      end do
!
!  Deallocate whole pointer
      deallocate(arrptr, STAT=istat)
      if (istat .ne. 0) destroyarrptr_i = .false.
!
      end function destroyarrptr_i
!
!
!
      function destroyarrptr_i2(arrptr)
      use feminterface, only:
      use femtypes
      implicit none
      type(ARRPTRI), pointer :: arrptr(:,:)
      logical :: destroyarrptr_i2
!
!------------------------------------------------------------------------------
!
!  destroy an array of pointers (derived data type) which contains integer data.
!
!------------------------------------------------------------------------------
!  Input:
!     arrptr    array of pointers (integer data)
!
!  Internal variables:
      integer (I4B) :: i, j, istat
!
      destroyarrptr_i2 = .true.
!  Deallocate individual entries
      do i = lbound(arrptr,1), ubound(arrptr,1)
        do j = lbound(arrptr,2), ubound(arrptr,2)
          deallocate(arrptr(i,j)%d, STAT=istat)
          if (istat .ne. 0) destroyarrptr_i2 = .false.
        end do
      end do
!
!  Deallocate whole pointer
      deallocate(arrptr, STAT=istat)
      if (istat .ne. 0) destroyarrptr_i2 = .false.
!
      end function destroyarrptr_i2
!
!
!
      function destroyarrptr_r(arrptr)
      use feminterface, only:
      use femtypes
      implicit none
      type(ARRPTRR), pointer :: arrptr(:)
      logical :: destroyarrptr_r
!
!------------------------------------------------------------------------------
!
!  destroy an array of pointers (derived data type) which contains real data.
!
!------------------------------------------------------------------------------
!  Input:
!     arrptr    array of pointers (real data)
!
!  Internal variables:
      integer (I4B) :: i, istat
!
      destroyarrptr_r = .true.
!  Deallocate individual entries
      do i = lbound(arrptr,1), ubound(arrptr,1)
        deallocate(arrptr(i)%d, STAT=istat)
        if (istat .ne. 0) destroyarrptr_r = .false.
      end do
!
!  Deallocate whole pointer
      deallocate(arrptr, STAT=istat)
      if (istat .ne. 0) destroyarrptr_r = .false.
!
      end function destroyarrptr_r
!
!
!
      function destroyarrptr_dpc(arrptr)
      use feminterface, only:
      use femtypes
      implicit none
      type(ARRPTRDPC), pointer :: arrptr(:)
      logical :: destroyarrptr_dpc
!
!------------------------------------------------------------------------------
!
!  destroy an array of pointers (derived data type) which contains complex data
!
!------------------------------------------------------------------------------
!  Input:
!     arrptr    array of pointers (complex data)
!
!  Internal variables:
      integer (I4B) :: i, istat
!
      destroyarrptr_dpc = .true.
!  Deallocate individual entries
      do i = lbound(arrptr,1), ubound(arrptr,1)
        deallocate(arrptr(i)%d, STAT=istat)
        if (istat .ne. 0) destroyarrptr_dpc = .false.
      end do
!
!  Deallocate whole pointer
      deallocate(arrptr, STAT=istat)
      if (istat .ne. 0) destroyarrptr_dpc = .false.
!
      end function destroyarrptr_dpc