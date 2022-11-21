      function reallocate_dp(p,n)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP), pointer :: p(:), reallocate_dp(:)
      integer (I4B) n
      intent(in) :: n
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
!    $Revision: 1.15 $
!    $Date: 2015/11/03 15:59:35 $
!    $Author: m_kasper $
!
!  resize a pointer array: allocate the array of size n and copy the former content
!   input:
!       p      array to resize
!       n      new length of array
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        return
      end if
      if (size(p) .eq. n) then
        reallocate_dp => p
        nullify(p)
      else
        allocate(reallocate_dp(n),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if
        reallocate_dp(1:min(size(p),n))=p(1:min(size(p),n))
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_dp
!
!
!
      function reallocate_dp2(p,n,m)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP), pointer :: p(:,:), reallocate_dp2(:,:)
      integer (I4B)      :: n, m
      intent(in)         :: n, m
!  resize a pointer array: allocate the array of size nxm and copy the former content
!   input:
!       p      two-dimensional array to resize
!       n      new length of array along the first dimension
!       m      new length of array along the second dimension
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        return
      end if
      if ((size(p,1) .eq. n) .and. (size(p,2) .eq. m)) then
        reallocate_dp2 => p
        nullify(p)
      else
        allocate(reallocate_dp2(n,m),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if
        reallocate_dp2(1:min(size(p,1),n),1:min(size(p,2),m))=p(1:min(size(p,1),n),1:min(size(p,2),m))
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_dp2
!
!
!
      function reallocate_sp(p,n)
      use feminterface, only:
      use femtypes
      implicit none
      real (SP), pointer :: p(:), reallocate_sp(:)
      integer (I4B) n
      intent(in) :: n
!  resize a pointer array: allocate the array of size n and copy the former content
!   input:
!       p      array to resize
!       n      new length of array
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        return
      end if
      if (size(p) .eq. n) then
        reallocate_sp => p
        nullify(p)
      else
        allocate(reallocate_sp(n),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if
        reallocate_sp(1:min(size(p),n))=p(1:min(size(p),n))
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_sp
!
!
!
      function reallocate_dpc(p,n)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC), pointer :: p(:), reallocate_dpc(:)
      integer (I4B) n
      intent(in) :: n
!  resize a pointer array: allocate the array of size n and copy the former content
!   input:
!       p      array to resize
!       n      new length of array
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        return
      end if
      if (size(p) .eq. n) then
        reallocate_dpc => p
        nullify(p)
      else
        allocate(reallocate_dpc(n),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if
        reallocate_dpc(1:min(size(p),n))=p(1:min(size(p),n))
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_dpc
!
!
!
      function reallocate_spc(p,n)
      use feminterface, only:
      use femtypes
      implicit none
      complex (SPC), pointer :: p(:), reallocate_spc(:)
      integer (I4B) n
      intent(in) :: n
!  resize a pointer array: allocate the array of size n and copy the former content
!   input:
!       p      array to resize
!       n      new length of array
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        return
      end if
      if (size(p) .eq. n) then
        reallocate_spc => p
        nullify(p)
      else
        allocate(reallocate_spc(n),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if
        reallocate_spc(1:min(size(p),n))=p(1:min(size(p),n))
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_spc
!
!
!
      function reallocate_i(p,n)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B), pointer :: p(:), reallocate_i(:)
      integer (I4B) n
      intent(in) :: n
!  resize a pointer array: allocate the array of size n and copy the former content
!   input:
!       p      array to resize
!       n      new length of array
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        nullify(reallocate_i)
        return
      end if
      if (size(p) .eq. n) then
        reallocate_i => p
        nullify(p)
      else
        allocate(reallocate_i(n),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if
        reallocate_i(1:min(size(p),n))=p(1:min(size(p),n))
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_i
!
!
!
      function reallocate_i2(p,n,m)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B), pointer :: p(:,:), reallocate_i2(:,:)
      integer (I4B) n, m
      intent(in) :: n, m
!  resize a pointer array: allocate the array of size nxm and copy the former content
!   input:
!       p      two-dimensional array to resize
!       n      new length of array along the first dimension
!       m      new length of array along the second dimension
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-asociated pointer array encountered in reallocate'
        return
      end if
      if ((size(p,1) .eq. n) .and. (size(p,2) .eq. m)) then
        reallocate_i2 => p
        nullify(p)
      else
        allocate(reallocate_i2(n,m),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if
        reallocate_i2(1:min(size(p,1),n),1:min(size(p,2),m))=p(1:min(size(p,1),n),1:min(size(p,2),m))
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_i2
!
!
!
      function reallocate_ch(p,n,txtlen)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) n, txtlen
      character (len=txtlen), pointer :: p(:), reallocate_ch(:)
      intent(in) :: n, txtlen
!  resize a pointer array: allocate the array of size n and copy the former content
!   input:
!       p      array to resize
!       n      new length of array
!       txtlen length of string
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        return
      end if
      allocate(reallocate_ch(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      reallocate_ch(1:min(size(p),n))=p(1:min(size(p),n))
      deallocate(p,stat=ierr)
      return
      end function reallocate_ch
!
!
!
      function reallocate_l(p,n)
      use feminterface, only:
      use femtypes
      implicit none
      logical, pointer :: p(:), reallocate_l(:)
      integer (I4B) n
      intent(in) :: n
!  resize a pointer array: allocate the array of size n and copy the former content
!   input:
!       p      array to resize
!       n      new length of array
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        return
      end if
      if (size(p) .eq. n) then
        reallocate_l => p
        nullify(p)
      else
        allocate(reallocate_l(n),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if
        reallocate_l(1:min(size(p),n))=p(1:min(size(p),n))
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_l
!
!
!
      function reallocate_l2(p,n,m)
      use feminterface, only:
      use femtypes
      implicit none
      logical, pointer :: p(:,:), reallocate_l2(:,:)
      integer (I4B) n, m
      intent(in) :: n, m
!  resize a pointer array: allocate the array of size n and copy the former content
!   input:
!       p      two-dimensional array to resize
!       n      new length of array along the first dimension
!       m      new length of array along the second dimension
!
      integer (I4B) ierr
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        return
      end if
      if ((size(p,1) .eq. n) .and. (size(p,2) .eq. m)) then
        reallocate_l2 => p
        nullify(p)
      else
        allocate(reallocate_l2(n,m),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if
        reallocate_l2(1:min(size(p,1),n),1:min(size(p,2),m))=p(1:min(size(p,1),n),1:min(size(p,2),m))
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_l2
!
!
!
      function reallocate_perr(p,n)
      use feminterface, only:
      use femtypes
      implicit none
      type(ARRPTRI), pointer :: p(:), reallocate_perr(:)
      integer (I4B) n
      intent(in) :: n
!  resize a pointer array of pointers: allocate the array of size n and copy the former content
!   input:
!       p      array to resize
!       n      new length of array
!
      integer (I4B) ierr,i
!
      if (.not. associated(p)) then
        print*,'Non-associated pointer array encountered in reallocate'
        return
      end if

      if (size(p) .eq. n) then
        reallocate_perr => p
        nullify(p)
      else
        allocate(reallocate_perr(n),stat=ierr)
        if (ierr .ne. 0) then
          print*,'Error in attempt to reallocate memory'
        end if

        do i = 1, min(size(p),n)
          reallocate_perr(i)%d => p(i)%d
        end do

        if (n.gt.size(p)) then ! Only if the array needs to be enlargened, nullify new entries:
          do i = size(p)+1, n
            nullify(reallocate_perr(i)%d)
          end do
        else if (n.lt.size(p)) then ! Only if the array needs to be reduced, deallocate entries
          do i = n+1, size(p)
            if (associated(p(i)%d)) then
              deallocate(p(i)%d)
            end if
          end do
        end if
        deallocate(p,stat=ierr)
      end if
      return
      end function reallocate_perr
!
!
!
      subroutine realloc_dp(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(DP), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      real(DP), allocatable :: temp(:)

      if (size(array) .eq. n) return
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array),n)) = array(1:min(size(array),n))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_dp
!
!
!
      subroutine realloc_dp_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(DP), allocatable :: arrayin(:)
      real(DP), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      real(DP), allocatable :: temp(:)

      if (size(arrayin) .eq. n) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin),n)) = arrayin(1:min(size(arrayin),n))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_dp_inout
!
!
!
      subroutine realloc_dp2(array, n, m)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      real(DP), allocatable :: array(:,:)
      intent(in) :: n
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      real(DP), allocatable :: temp(:,:)

      if ((size(array,1) .eq. n) .and. (size(array,2) .eq. m)) return
      allocate(temp(n,m),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array,1),n),1:min(size(array,2),m)) = array(1:min(size(array,1),n),1:min(size(array,2),m))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_dp2
!
!
!
      subroutine realloc_dp2_inout(arrayin, n, m, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      real(DP), allocatable :: arrayin(:,:)
      real(DP), allocatable :: arrayout(:,:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      real(DP), allocatable :: temp(:,:)

      if ((size(arrayin,1) .eq. n) .and. (size(arrayin,2) .eq. m)) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n,m),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin,1),n),1:min(size(arrayin,2),m)) = arrayin(1:min(size(arrayin,1),n),1:min(size(arrayin,2),m))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_dp2_inout
!
!
!
      subroutine realloc_sp(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(SP), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      real(SP), allocatable :: temp(:)

      if (size(array) .eq. n) return
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array),n)) = array(1:min(size(array),n))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_sp
!
!
!
      subroutine realloc_sp_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(SP), allocatable :: arrayin(:)
      real(SP), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      real(SP), allocatable :: temp(:)

      if (size(arrayin) .eq. n) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin),n)) = arrayin(1:min(size(arrayin),n))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_sp_inout
!
!
!
      subroutine realloc_dpc(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(DPC), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      complex(DPC), allocatable :: temp(:)

      if (size(array) .eq. n) return
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array),n)) = array(1:min(size(array),n))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_dpc
!
!
!
      subroutine realloc_dpc_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(DPC), allocatable :: arrayin(:)
      complex(DPC), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      complex(DPC), allocatable :: temp(:)

      if (size(arrayin) .eq. n) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin),n)) = arrayin(1:min(size(arrayin),n))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_dpc_inout
!
!
!
      subroutine realloc_spc(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(SPC), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      complex(SPC), allocatable :: temp(:)

      if (size(array) .eq. n) return
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array),n)) = array(1:min(size(array),n))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_spc
!
!
!
      subroutine realloc_spc_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(SPC), allocatable :: arrayin(:)
      complex(SPC), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      complex(SPC), allocatable :: temp(:)

      if (size(arrayin) .eq. n) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin),n)) = arrayin(1:min(size(arrayin),n))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_spc_inout
!
!
!
      subroutine realloc_i(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      integer(I4B), allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      integer(I4B), allocatable :: temp(:)

      if (size(array) .eq. n) return
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array),n)) = array(1:min(size(array),n))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_i
!
!
!
      subroutine realloc_i_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      integer(I4B), allocatable :: arrayin(:)
      integer(I4B), allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      integer(I4B), allocatable :: temp(:)

      if (size(arrayin) .eq. n) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin),n)) = arrayin(1:min(size(arrayin),n))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_i_inout
!
!
!
      subroutine realloc_i2(array, n, m)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      integer(I4B), allocatable :: array(:,:)
      intent(in) :: n
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      integer(I4B), allocatable :: temp(:,:)

      if ((size(array,1) .eq. n) .and. (size(array,2) .eq. m)) return
      allocate(temp(n,m),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array,1),n),1:min(size(array,2),m)) = array(1:min(size(array,1),n),1:min(size(array,2),m))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_i2
!
!
!
      subroutine realloc_i2_inout(arrayin, n, m, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      integer(I4B), allocatable :: arrayin(:,:)
      integer(I4B), allocatable :: arrayout(:,:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      integer(I4B), allocatable :: temp(:,:)

      if ((size(arrayin,1) .eq. n) .and. (size(arrayin,2) .eq. m)) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n,m),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin,1),n),1:min(size(arrayin,2),m)) = arrayin(1:min(size(arrayin,1),n),1:min(size(arrayin,2),m))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_i2_inout
!
!
!
      subroutine realloc_l(array, n)
      use femtypes
      implicit none
      integer(I4B) :: n
      logical, allocatable :: array(:)
      intent(in) :: n
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      logical, allocatable :: temp(:)

      if (size(array) .eq. n) return
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array),n)) = array(1:min(size(array),n))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_l
!
!
!
      subroutine realloc_l_inout(arrayin, n, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n
      logical, allocatable :: arrayin(:)
      logical, allocatable :: arrayout(:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      logical, allocatable :: temp(:)

      if (size(arrayin) .eq. n) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin),n)) = arrayin(1:min(size(arrayin),n))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_l_inout
!
!
!
      subroutine realloc_l2(array, n, m)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      logical, allocatable :: array(:,:)
      intent(in) :: n
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      logical, allocatable :: temp(:,:)

      if ((size(array,1) .eq. n) .and. (size(array,2) .eq. m)) return
      allocate(temp(n,m),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array,1),n),1:min(size(array,2),m)) = array(1:min(size(array,1),n),1:min(size(array,2),m))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_l2
!
!
!
      subroutine realloc_l2_inout(arrayin, n, m, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n, m
      logical, allocatable :: arrayin(:,:)
      logical, allocatable :: arrayout(:,:)
      intent(in) :: n
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      logical, allocatable :: temp(:,:)

      if ((size(arrayin,1) .eq. n) .and. (size(arrayin,2) .eq. m)) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n,m),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin,1),n),1:min(size(arrayin,2),m)) = arrayin(1:min(size(arrayin,1),n),1:min(size(arrayin,2),m))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_l2_inout
!
!
!
      subroutine realloc_ch(array, n, txtlen)
      use femtypes
      implicit none
      integer(I4B) :: n, txtlen
      character (len=txtlen), allocatable :: array(:)
      intent(in) :: n, txtlen
      intent(inout) :: array
!  resize an allocatable array: allocate the array of size n and copy the former content
!   input:
!       array  array to resize
!       n      new length of array
!
!  local variables
      integer(I4B) :: ierr
      character (len=txtlen), allocatable :: temp(:)

      if (size(array) .eq. n) return
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(array),n)) = array(1:min(size(array),n))
      call move_alloc(from=temp, to=array)

      return
      end subroutine realloc_ch
!
!
!
      subroutine realloc_ch_inout(arrayin, n, txtlen, arrayout)
      use femtypes
      implicit none
      integer(I4B) :: n, txtlen
      character (len=txtlen), allocatable :: arrayin(:)
      character (len=txtlen), allocatable :: arrayout(:)
      intent(in) :: n, txtlen
      intent(out) :: arrayout
      intent(inout) :: arrayin
!  resize an allocatable array: allocate the array of size n and copy the former content tp arrayout
!   input:
!       arrayin   input array to resize
!       n         new length of array
!       arrayout  output array of size newsize with values of arrayin
!
!  local variables
      integer(I4B) :: ierr
      character (len=txtlen), allocatable :: temp(:)

      if (size(arrayin) .eq. n) then
        call move_alloc(from=arrayin, to=arrayout)
        return
      end if
      allocate(temp(n),stat=ierr)
      if (ierr .ne. 0) then
        print*,'Error in attempt to reallocate memory'
      end if
      temp(1:min(size(arrayin),n)) = arrayin(1:min(size(arrayin),n))
      deallocate(arrayin)
      call move_alloc(from=temp, to=arrayout)

      return
      end subroutine realloc_ch_inout
