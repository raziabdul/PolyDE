subroutine getfe(fe)
      use globalvariables3D,         only: vf, ve, numf, numv
      use femtypes
      implicit none
      integer  (I4B), allocatable :: fe(:,:)
      intent   (out)              :: fe
!  Calculate the fe-Set. Assignment between faces and their edges (all global numbers)
!
!  Output:
!            fe(3,numf)  List if faces <-> edges
! 
!  local variables:
      integer (I4B)              :: f, v
      integer (I4B),   parameter :: f2e(3,4) = reshape((/2,5,6,3,4,6,1,4,5,1,2,3/),(/3,4/))

      allocate(fe(3,numf))
      fe(:,:) = 0

      do v = 1, numv
        do f = 1, 4
          if (fe(1,vf(f,v)) .eq. 0) then
            fe(1:3,vf(f,v)) = ve(f2e(1:3,f),v)
          end if
        end do
      end do

      return
end subroutine getfe