subroutine getse(se)
      use globalvariables3D,         only: vv, ve, numv, nums
      use femtypes
      implicit none
      integer  (I4B), allocatable :: se(:,:)
      intent   (out)              :: se
!  Calculate the fe-Set. Assignment between faces and their edges (all global numbers)
!
!  Output:
!            se(3,nums)  List if surfaces <-> edges
! 
!  local variables:
      integer (I4B)              :: f, v
      integer (I4B),   parameter :: f2e(3,4) = reshape((/2,5,6,3,4,6,1,4,5,1,2,3/),(/3,4/))

      allocate(se(3,nums))
      se(:,:) = -1

      do v = 1, numv
        do f = 1, 4
          if (vv(f,v).lt.0) then
            se(1:3,abs(vv(f,v))) = ve(f2e(1:3,f),v)
          end if
        end do
      end do

      return
end subroutine getse