      subroutine getvf(numf,vf)
      use femtypes
      use globalvariables3D, only: numv, vv
      implicit none
      integer (I4B) :: numf
      integer (I4B), pointer :: vf(:,:)
      intent (out) :: numf, vf
!
!------------------------------------------------------------------------------
!    $Revision: 1.3 $
!    $Date: 2015/11/11 17:15:29 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Subroutine computes maps for volumes -> edges, edges -> nodes, neighbours
!  and number of edges. Number of faces and volume -> faces map can be obtained
!  as optional variables.
!
!------------------------------------------------------------------------------
!  Output:
!     numf      number of faces
!     vf        volume element -> faces    allocate(vf(4,numv))
!
!  Internal variables:
      integer (I4B) :: i, j, k       ! counters
      integer (I4B) :: neighbour
!
!------------------------------------------------------------------------------
!  COMPUTE VF:
!  -----------
!
!  Allocate and initialize vf
      if (associated(vf)) deallocate(vf)
      allocate(vf(4,numv))
      vf = 0
!
      numf = 0
      do i = 1,numv
        do j = 1,4
          neighbour = vv(j,i)
!  we only check neighbours with index smaller than element (or 0 for surface)
          if (neighbour .lt. i) then
!  a new face
            numf = numf + 1
            vf(j,i) = numf
!  if not at the surface ... 
            if (neighbour .gt. 0) then
!  we have a neighbor, visit him
              do k = 1,4
!  find the face in the neighbour element
                if (vv(k,neighbour) .eq. i) then
!  found the element in the list of neighbours of the neighbour
!  assign the face number at the neighbour
                  vf(k,neighbour) = numf
                  exit
                end if
              end do
            end if
          end if
        end do
      end do
!
      end subroutine getvf
