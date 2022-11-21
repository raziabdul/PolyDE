      subroutine bcstovv
      use feminterface3D, only: 
      use femtypes
      use globalvariables3D, only: nums, numv, sn, vn, vv
      implicit none
!
!------------------------------------------------------------------------------
!    $Revision: 1.6 $
!    $Date: 2014/05/23 15:04:33 $
!    $Author: juryanatzki $
!------------------------------------------------------------------------------
!
!  Map the boundary condition index to vv array as a negative value. If an
!  element has no neighbour, the original vv array contains a zero and the face
!  is located on a surface.
!  vv will obtain the negative number of the surface instead of a zero, which
!  is connected to the boundary condition at the same time.
!
!  if (vv(i,elem) .lt. 0) then the boundary condition is bcs(abs(vv(i,elem))).
!
!------------------------------------------------------------------------------
!  local variables:
      integer (I4B) :: i, j, k
      integer (I4B) :: surfnod(3)
      integer (I4B), parameter :: fton(3,4)=reshape((/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/))
!
      do i = 1,numv  ! volume elements
        do j = 1,4  ! faces
!  If no neighbour, we got a surface ==> assign nodes of surface
          if (vv(j,i) .eq. 0) then
            surfnod = (/vn(fton(1,j),i),vn(fton(2,j),i),vn(fton(3,j),i)/)
!  Reorder surface nodes into ascending order.
            call order3(surfnod(1),surfnod(2),surfnod(3))
            do k = 1,nums  ! surface elements
!  Compare surface nodes with sn array entries.
              if ( (surfnod(1) .eq. sn(1,k)) .and. &
                   (surfnod(2) .eq. sn(2,k)) .and. &
                   (surfnod(3) .eq. sn(3,k)) ) then
!  Store negative number of surface element k in vv(j,i).
                vv(j,i) = -k
!                exit
                goto 10
              end if
            end do  ! surface elements
            print*,'Error in BCSTOVV volume',i,' has a surface with nodes:',vn(fton(1,j),i),vn(fton(2,j),i),vn(fton(3,j),i)
            print*,'which was not found in the list of surface elements'
            pause
           stop
10 continue            
!  If there is a neighbour go to next face.
          else
            cycle
          end if
        end do  ! faces
      end do  ! volume elements
!
      end subroutine bcstovv