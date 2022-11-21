      pure subroutine getgl(vert, gl, n, t1, t2)
      use femtypes
      use feminterface, only: cross_product
      implicit none
      real (DP) :: vert(3,4), gl(3,4)
      real (DP), optional :: n(3,4), t1(3,4), t2(3,4)
      intent (in) :: vert
      intent (out) :: gl, n, t1, t2
!
!    $Revision: 1.11 $
!    $Date: 2014/08/22 10:52:17 $
!    $Author: m_kasper $
!
!  Calculation of the gradients 
!    and optionally compute face unit vectors t1, t2, n
!  possible calls are either:      call getgl (vert,gl)
!                          or:     call getgl (vert,gl,n)
!                          or:     call getgl (vert,gl,n,t1,t2)
!
!  Gradients:  gl(u,i) is the u-component of the gradient of lambda_i
!  u = {x,y,z} = {1,2,3}; i=1..4
!
!  Input:
!          vert       coordinates of the vertices of tetrahedron
!
!  Output:
!          gl         gradients of barycentric coordinates
!          n          is the outward normal for the face, opposite to gl
!          t1         is a vector along the first edge of the face, 
!                     pointing toward the node of highest index
!          t2         is a vector perpendicular to the first edge of the face 
!                     and pointing toward the first vertex of the face
!                     n, t1, t2 form a right hand system such that t1 x t2 = n
!
!  local variables 

      real (DP) :: volume6, del(3,3)
      integer (I4B) :: face
! permutations array  nodes of faces of tetrahedron, e.g. face 1: 2,4,3
      integer (I4B), parameter, dimension (3,4) ::                      &
     &         fnlocal=reshape((/3,2,4, 1,3,4, 2,1,4, 1,2,3/), (/3,4/) )
!
!  for faces 1 to 3 we compute relative coordinates with respect to vertex 4,
!  which is a member of all these faces (but not of face 4)
      del(1:3,1) = vert(1:3,1) - vert(1:3,4)
      del(1:3,2) = vert(1:3,2) - vert(1:3,4)
      del(1:3,3) = vert(1:3,3) - vert(1:3,4)
!  compute gradients (of volume) for faces 1 to 3; face 4 will be obtained later
!      forall( face = 1:3 )
!!  use permutation array in computing product terms
!        gl(:,face)=  cross_product( del(:,fnlocal(2,face)),             &
!     &                              del(:,fnlocal(1,face)) )
!      end forall
      gl(1,1) = del(3,3)*del(2,2)-del(3,2)*del(2,3)
      gl(2,1) = del(1,3)*del(3,2)-del(1,2)*del(3,3)
      gl(3,1) = del(2,3)*del(1,2)-del(2,2)*del(1,3)
      gl(1,2) = del(3,1)*del(2,3)-del(3,3)*del(2,1)
      gl(2,2) = del(1,1)*del(3,3)-del(1,3)*del(3,1)
      gl(3,2) = del(2,1)*del(1,3)-del(2,3)*del(1,1)
      gl(1,3) = del(3,2)*del(2,1)-del(3,1)*del(2,2)
      gl(2,3) = del(1,2)*del(3,1)-del(1,1)*del(3,2)
      gl(3,3) = del(2,2)*del(1,1)-del(2,1)*del(1,2)
!  compute six times tetrahedron volume
      volume6 = dot_product(del(1:3,1),gl(1:3,1))
      gl(:,1:3)=gl(:,1:3)/volume6
!  Sum of gradients is zero 
      gl(:,4) = - ( gl(:,1)+gl(:,2)+gl(:,3) )
!
      if (.not. present(n)) return
!
!  compute face unit vectors n
      forall( face = 1:4 )
        n(:,face) = - gl(:,face) / sqrt(dot_product(gl(:,face),gl(:,face)))
      end forall
!
      if (.not. present(t1)) return
!
!  compute face unit vectors t1
      forall( face = 1:3 )
        t1(:,face)=   del(:,fnlocal(2,face))
        t1(:,face)= - t1(:,face) / sqrt(dot_product(t1(:,face),t1(:,face)))
      end forall
!  seperated treatment of face 4
      t1(:,4) = vert(:,3) - vert(:,2)
      t1(:,4) = t1(:,4) / sqrt(dot_product(t1(:,4),t1(:,4)))
!  compute face unit vectors t2 = n x t1
      forall (face=1:4)
        t2(:,face)=cross_product(n(:,face),t1(:,face))
      end forall

      return
      end subroutine getgl
