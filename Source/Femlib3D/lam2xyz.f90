      subroutine lam2xyz(lambda, elem, xyzt, nod, vn)
      use femtypes
      implicit none
      integer (I4B) elem, vn(:,:)
      real (DP) xyzt(3), lambda(4), nod(:,:)
      intent (in) :: lambda, nod, elem, vn
      intent (out) :: xyzt

!-----------------------------------------------------------------------------
!    $Revision: 1.7 $
!    $Date: 2015/06/17 11:47:06 $
!    $Author: jimmykamboh $
!-----------------------------------------------------------------------------
!  conversion of barycentric coordinates lambda to x,y,z coordinates in element 
!  volume elem  
!
!  Input:
!            lambda      barycentric coordinates of point x,y,z 
!            nod         coordinates of tetrahedron nodes
!            vn          volume element information (nodes of the element)
!            elem        volume element number with respect to which coordinates
!                        are computed
!  Output:
!            xyzt        world coordinates of the point 
!

      xyzt = matmul( nod(:,vn(:,elem)) , lambda(:) )
      return
      end subroutine lam2xyz


      
      subroutine xyz2lam(lambda, elem, xyzt, nod, vn)
      use femtypes
      implicit none
      integer (I4B) :: elem, vn(:,:)
      real (DP) :: xyzt(3), lambda(4), nod(:,:)
      intent (in) :: xyzt, nod, elem, vn
      intent (out) :: lambda
!
!  conversion of x,y,z coordinates to barycentric coordinates lambda in element 
!  volume elem  
!
!  Input:
!            xyzt        world coordinates of the point 
!            nod         coordinates of tetraehdron nodes
!            vn          volume element information (nodes of the element)
!            elem        volume element number with respect to which coordinates
!                        are computed
!  Output:
!            lambda      barycentric coordinates of point x, y, z
!
!  local variables
      real (DP) :: a(3), b(3), c(3), d(3), s(3), det(4), vol6
!
      a = nod( :,vn(1,elem) ) - xyzt(:)
      b = nod( :,vn(2,elem) ) - xyzt(:)
      c = nod( :,vn(3,elem) ) - xyzt(:)
      d = nod( :,vn(4,elem) ) - xyzt(:)
!
      s(1)=c(2)*d(3) - c(3)*d(2)
      s(2)=c(3)*d(1) - c(1)*d(3)
      s(3)=c(1)*d(2) - c(2)*d(1)
!  det contains the volume of a tetrahedron with faces (1:4) and the point xyzt
      det(1)=-dot_product( b(1:3) , s(1:3) )
      det(2)= dot_product( a(1:3) , s(1:3) )
!
      s(1)=a(2)*b(3) - a(3)*b(2)
      s(2)=a(3)*b(1) - a(1)*b(3)
      s(3)=a(1)*b(2) - a(2)*b(1)
!
      det(3)=-dot_product( d(1:3) , s(1:3) )
      det(4)= dot_product( c(1:3) , s(1:3) )
!  this is six times the tetrahedron volume
      vol6=sum(det(1:4))
      
      lambda=det/vol6
      return
      end subroutine xyz2lam
