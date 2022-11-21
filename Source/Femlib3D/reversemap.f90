      subroutine reversemap(column,array,arrptr,nmin,nmax)
      use feminterface, only: reallocate
      use femtypes
      use globalvariables3D, only:
      implicit none
      integer (I4B) :: nmin, nmax
      integer (I4B) :: array(:,:)
      type (ARRPTRI), pointer :: arrptr(:)
      logical :: column
      intent(in) :: column
      intent(out) :: nmin, nmax
!
!------------------------------------------------------------------------------
!    $Revision: 1.8 $
!    $Date: 2015/11/11 17:21:58 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Subroutine computes reverse map for an array of integer values. It does a 
!  loop over the array and writes the column/row of integer value i into the
!  next free entry in arrptr(i)%d. The absolute value of the integers will we
!  taken for the reverse map.
!
!  e.g. take the array:
!  --------------------
!          / 1 1 1 4 5 \
!  array = | 2 3 4 5 1 |
!          \ 3 4 5 1 2 /
!
!  The arrptr would look like for column = .true.:
!  -----------------------------------------------
!  arrptr(1)%d = ( 1 2 3 4 5 )
!  arrptr(2)%d = ( 1 5 )
!  arrptr(3)%d = ( 1 2 )
!  arrptr(4)%d = ( 2 3 4 )
!  arrptr(5)%d = ( 3 4 5 )
!
!  The arrptr would look like for column = .false.:
!  ------------------------------------------------
!  arrptr(1)%d = ( 1 1 1 3 2 )
!  arrptr(2)%d = ( 2 3 )
!  arrptr(3)%d = ( 3 2 )
!  arrptr(4)%d = ( 3 2 1 )
!  arrptr(5)%d = ( 3 2 1 )

!
!  Additionally nmin = 1 and nmax = 5 will be returned.
!
!  The routine may be used to get a map of nodes --> volume/surface elements
!  from a map of volume/surface elements --> nodes or node --> edge from
!  edge --> node.
!
!------------------------------------------------------------------------------
!  Input:
!    array      array given to subroutine, that contains integer values only
!    column     if .true. arrptr contains the column of integer value i
!               if .false. arrptr contains the row of integer value i
!
!  Output:
!    arrptr     reverse mapped array of pointers
!    nmin       minimum value of array
!    nmax       maximum value of array
!
!  Internal variables:
      integer (I4B) :: i, j, k
      integer (I4B) :: m, n
      integer (I4B) :: initsize
      integer (I4B), allocatable :: length(:)
!
!  Get size of array(m,n)
      m = size(array,1)
      n = size(array,2)
!  Check for the maximum value of the array and allocate array of pointers.
      nmin = minval(abs(array))
      nmax = maxval(abs(array))
      if (associated(arrptr)) then
        print*,'Array Pointer must not be associated (reversemap)'
      end if

      allocate(arrptr(nmin:nmax))
!  estimate for appropriate average length
      initsize = 1.3*m*n/size(arrptr)
      do i = nmin,nmax
        allocate(arrptr(i)%d(initsize))
        arrptr(i)%d(1:initsize) = 0
      end do

      allocate(length(nmin:nmax))
      length = 0

      if (column) then
!  column loop
        do j = 1,n
!  row loop
          do k = 1,m
            i = abs(array(k,j))
            length(i) = length(i)+1
!  reallocate pointer if entry would be written as overflow
            if (length(i) .gt. size(arrptr(i)%d)) then
              arrptr(i)%d => reallocate(arrptr(i)%d,2*length(i))
            end if
            arrptr(i)%d(length(i)) = j
          end do
        end do
      else
!  column loop
        do j = 1,n
!  row loop
          do k = 1,m
            i = abs(array(k,j))
            length(i) = length(i)+1
!  reallocate pointer if entry would be written as overflow
            if (length(i) .gt. size(arrptr(i)%d)) then
              arrptr(i)%d => reallocate(arrptr(i)%d,2*length(i))
            end if
!  write row of value to pointer
            arrptr(i)%d(length(i)) = k
          end do
        end do
      end if

!  resize pointer size, if size is more than length needed by nonzero entries
      do i = nmin,nmax
        if (length(i) .lt. size(arrptr(i)%d)) then
          arrptr(i)%d => reallocate(arrptr(i)%d,length(i))
        end if
!------------------------------------------------------------------------------
!        if (column) then
!          print "(a,i4,a)","node ",i," occurs in column"
!        else
!          print "(a,i4,a)","node ",i," occurs in row"
!        end if
!        print *,arrptr(i)%d
!------------------------------------------------------------------------------
      end do

      deallocate(length)
!
      end subroutine reversemap
