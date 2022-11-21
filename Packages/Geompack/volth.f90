function volth ( a, b, c, d )
!
!******************************************************************************
!
!! VOLTH computes the volume of a tetrahedron.
!
!
!  Purpose:
!
!    Compute 6 times volume of tetrahedron.
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, A(1:3), B(1:3), C(1:3), D(1:3) - 4 vertices of tetrahedron.
!
!    Output, VOLTH - 6 * volume of tetrahedron.
!
!  Changes by Razi:
!
!    - removed abs(...) so that volume can be positive/negative to denote the tetrahedron's orientation.
!    - swapped indices so that the function returns positive volume for LEFT-hand-side oriented 
!      tetrahedron, i.e., Face 4, defined by Vertices 1,2,3 is mathematically NEGATIVE viewing from 
!      inside the tet, so Vertex 4 is pointed from Face 4 using the LEFT-hand rule. 
!
  implicit none

  real ( kind = 8 ) a(3)
  real ( kind = 8 ) b(3)
  real ( kind = 8 ) c(3)
  real ( kind = 8 ) d(3)
  real ( kind = 8 ) u(3)
  real ( kind = 8 ) v(3)
  real ( kind = 8 ) volth
  real ( kind = 8 ) w(3)

  u(1:3) = b(1:3) - a(1:3)
  v(1:3) = c(1:3) - a(1:3)
  w(1:3) = d(1:3) - a(1:3)

  volth = u(1)*(v(3)*w(2) - v(2)*w(3)) + u(2)*(v(1)*w(3) - &
    v(3)*w(1)) + u(3)*(v(2)*w(1) - v(1)*w(2))

  return
end
