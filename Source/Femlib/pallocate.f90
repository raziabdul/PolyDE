      subroutine palloc_dp(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(DP), pointer:: ptr(:)
      intent(in) :: n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB_d'
        if (size(ptr) .ne. n) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(n))
        end if
      else
        allocate(ptr(n))
      end if

      return
      end


      subroutine palloc_sp(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      real(SP), pointer:: ptr(:)
      intent(in) :: n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB_s'
        if (size(ptr) .ne. n) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(n))
        end if
      else
        allocate(ptr(n))
      end if

      return
      end


      subroutine palloc_dpc(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      complex(DPC), pointer:: ptr(:)
      intent(in) :: n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB_c'
        if (size(ptr) .ne. n) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(n))
        end if
      else
        allocate(ptr(n))
      end if

      return
      end


      subroutine palloc_spc(ptr,n)
!  alloataion of a pointer array includes checking of allocation status and
      use femtypes
      use feminterface, only: reallocate
      implicit none
      integer(I4B) :: n
      complex(SPC), pointer:: ptr(:)
      intent(in) :: n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB_c'
        if (size(ptr) .ne. n) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(n))
        end if
      else
        allocate(ptr(n))
      end if

      return
      end


      subroutine palloc_i(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      integer(I4B), pointer:: ptr(:)
      intent(in) :: n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB_i'
        if (size(ptr) .ne. n) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(n))
        end if
      else
        allocate(ptr(n))
      end if

      return
      end


      subroutine palloc_l(ptr,n)
      use femtypes
      implicit none
      integer(I4B) :: n
      logical, pointer:: ptr(:)
      intent(in) :: n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB_l'
        if (size(ptr) .ne. n) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(n))
        end if
      else
        allocate(ptr(n))
      end if

      return
      end


      subroutine palloc_ptri(ptr,n)
      use femtypes
      use feminterface, only : destroyarrptr
      implicit none
      integer(I4B) :: n
      type(ARRPTRI), pointer:: ptr(:)
      intent(in) :: n
!  allocation of a pointer array including checking of allocation status
!
!  local variables
      integer(I4B) :: i
      logical :: ok
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB_p'
        if (size(ptr) .ne. n) then
!  the array is already allocated
          ok = destroyarrptr(ptr)
          allocate(ptr(n))
        end if
      else
        allocate(ptr(n))
        do i= 1,n
          nullify(ptr(i)%d)
        end do
      end if

      return
      end


      subroutine palloc_dp2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      real(DP), pointer:: ptr(:,:)
      intent(in) :: m, n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB'
        if ((size(ptr,1) .ne. m) .or. (size(ptr,2) .ne. n)) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(m,n))
        end if
      else
        allocate(ptr(m,n))
      end if

      return
      end


      subroutine palloc_sp2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      real(SP), pointer:: ptr(:,:)
      intent(in) :: m, n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB'
        if ((size(ptr,1) .ne. m) .or. (size(ptr,2) .ne. n)) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(m,n))
        end if
      else
        allocate(ptr(m,n))
      end if

      return
      end


      subroutine palloc_dpc2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      complex(DPC), pointer:: ptr(:,:)
      intent(in) :: m, n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB'
        if ((size(ptr,1) .ne. m) .or. (size(ptr,2) .ne. n)) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(m,n))
        end if
      else
        allocate(ptr(m,n))
      end if

      return
      end


      subroutine palloc_spc2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      complex(SPC), pointer:: ptr(:,:)
      intent(in) :: m, n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB'
        if ((size(ptr,1) .ne. m) .or. (size(ptr,2) .ne. n)) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(m,n))
        end if
      else
        allocate(ptr(m,n))
      end if

      return
      end


      subroutine palloc_i2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      integer(I4B), pointer:: ptr(:,:)
      intent(in) :: m, n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB'
        if ((size(ptr,1) .ne. m) .or. (size(ptr,2) .ne. n)) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(m,n))
        end if
      else
        allocate(ptr(m,n))
      end if

      return
      end


      subroutine palloc_l2(ptr,m,n)
      use femtypes
      implicit none
      integer(I4B) :: m, n
      logical, pointer:: ptr(:,:)
      intent(in) :: m, n
!  allocation of a pointer array including checking of allocation status
!
      if (associated(ptr)) then
! for Debugging purpose
write(*,'(F9.3,A)') real(storage_size(ptr)/8*size(ptr))/(1024),'kB'
        if ((size(ptr,1) .ne. m) .or. (size(ptr,2) .ne. n)) then
!  the array is already allocated
          deallocate(ptr)
          allocate(ptr(m,n))
        end if
      else
        allocate(ptr(m,n))
      end if

      return
      end
