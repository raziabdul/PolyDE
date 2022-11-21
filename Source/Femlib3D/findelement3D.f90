! 
      subroutine findelement3D(xyzt,nod,vn,vv,numv,ielem,found)
      use feminterface3d, only: xyz2lam, findelement_x3D
      use femtypes
      implicit none
      real (DP) :: xyzt(:), nod(:,:)
      integer (I4B) :: vn(:,:), vv(:,:), numv, ielem
      logical :: found
      intent (in) :: xyzt, nod, vn, vv, numv
      intent (out) :: ielem, found
!
!    $Revision: 1.3 $
!    $Date: 2015/11/11 17:27:12 $
!    $Author: m_kasper $
!
!  Find the element number in which the point xyzt is located
!  Input:    xyzt   coordinates of the point
!            nod    coordinates of mesh nodes         
!            vn     nodes of volume
!            vv     neighbours of elements
!            numv   number of elements
!  Output:   ielem  element number in which the node was found
!                   if found = a nearby element is returned
!            found  = true if the point was found within an element
!                   = false point is outside the triangulated region
!
!  Starting from the element number found by the last call a sequence of neighbour 
!  elements is visited directed toward the point. This however may not succeed if 
!  the meshed region has concave surface or the node is outside. 
!  If not successful (or vv is not present) an exhaustive search is used.
!
!  local variables
      real (DP) :: lambda(4)
      integer (I4B) :: i, imin(1), inew, iold
      integer (I4B), save :: ilea=1
!
!  Element returned by last call
      i=ilea
!
      iold=0
      do
        call xyz2lam(lambda, i, xyzt, nod, vn)
        imin=minloc(lambda)
        if (lambda(imin(1)) .lt. 0._DP) then
          inew=vv(imin(1),i)
!  check to continue with neighbour
!  and prevent from toggling between two elements
          if ( (inew.gt.0) .and. (inew.ne.iold) ) then
            iold=i
            i=inew
          else
            found=.false.
            exit
          end if
        else
!  found
          found=.true.
          exit
        end if
      end do
!
      if (found) then
        ielem=i
     else
!  not found by directed sear now use exhaustive search
        call findelement_x3D(xyzt,nod,vn,numv,ielem,found)
        if (.not. found) ielem=i
      end if
      ilea=ielem
      return
      end subroutine findelement3D
!
!
!
      subroutine findelement_x3D(xyzt,nod,vn,numv,ielem,found)
      use feminterface3d, only: xyz2lam
      use femtypes
      real (DP) :: xyzt(:), nod(:,:)
      integer (I4B) :: vn(:,:), numv, ielem
      logical :: found
      intent (in) :: xyzt, nod, vn, numv
      intent (out) :: ielem, found
!
!  Find the element number in which the point xyzt is located (exhausive search)
!  Input:    xyzt   coordinates of the point
!            nod    coordinates of mesh nodes         
!            vn     nodes of volume
!            numv   number of elements
!  Output:   ielem  element number in which the node was found
!            found  = true if the point was found within an element
!                   = false point is outside the triangulated region
!
!  local variables
      real (DP) :: lambda(4)
      integer (I4B) :: elem
!
      found=.false.
!  Test all elements (until found) whether point xyzt is inside
      do elem = 1,numv
        call xyz2lam(lambda, elem, xyzt, nod, vn)
        if ( all(lambda .ge. 0._DP) ) then
          found=.true.
          ielem = elem
          exit
        end if
      end do
      return
      end subroutine findelement_x3D
