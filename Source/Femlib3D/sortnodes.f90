      subroutine sortnodes
      use feminterface3D, only: compareandswap
      use femtypes
      use globalvariables3D, only: vn, numv, sn, nums
      implicit none
!
!------------------------------------------------------------------------------
!    $Revision: 1.7 $
!    $Date: 2014/08/22 10:40:45 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  For each of the volume/surface elements sort the nodes in ascending order 
!  such that:   vn(1,i)  <  vn(2,i)  <  vn(3,i)  <  vn(4,i)
!  such that:   sn(1,i)  <  sn(2,i)  <  sn(3,i)
!
!
!  local variables:
!      integer (I4B) :: i
!
!!  this loop can be executed in parallel
!      do i=1,numv
!!  a list of 4 items cannot be sorted using fewer than 5 comparisons
!        call compareandswap(vn(1,i),vn(2,i))
!        call compareandswap(vn(3,i),vn(4,i))
!        call compareandswap(vn(1,i),vn(3,i))
!        call compareandswap(vn(2,i),vn(4,i))
!        call compareandswap(vn(2,i),vn(3,i))
!      end do
!
!  sort volume element -> nodes map
!  a list of 4 items cannot be sorted using fewer than 5 comparisons
      call compareandswap(vn(1,1:numv),vn(2,1:numv))
      call compareandswap(vn(3,1:numv),vn(4,1:numv))
      call compareandswap(vn(1,1:numv),vn(3,1:numv))
      call compareandswap(vn(2,1:numv),vn(4,1:numv))
      call compareandswap(vn(2,1:numv),vn(3,1:numv))
!
!  sort surface element -> nodes map
!  a list of 3 items cannot be sorted using fewer than 3 comparisons
      call compareandswap(sn(1,1:nums),sn(2,1:nums))
      call compareandswap(sn(2,1:nums),sn(3,1:nums))
      call compareandswap(sn(1,1:nums),sn(2,1:nums))
!
      end subroutine sortnodes



      elemental subroutine compareandswap(a,b)
      use feminterface3D, only:
      use femtypes
      implicit none
      integer (I4B) :: a, b
      intent (inout) :: a, b
!
!  compare if  a > b   and if so interchange these
!
!  in-/ output:
!    a, b           integers to be compared
!
!  internal variables:
      integer (I4B) :: temp
!
      if (a .gt. b) then
        temp=a
        a=b
        b=temp
      end if 
      return
      end subroutine compareandswap
