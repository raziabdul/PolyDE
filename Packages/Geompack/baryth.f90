subroutine baryth ( a, b, c, d, e, alpha, degen )

!******************************************************************************
!
!! BARYTH computes barycentric coordinates of a point in 3D.
!
!  Purpose: 
!
!    Compute barycentric coordinates of 3D point with respect
!    to four vertices of a tetrahedron.
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
!    Input, A(1:3),B(1:3),C(1:3),D(1:3) - 4 vertices of tetrahedron.
!
!    Input, E(1:3) - fifth point for which barycentric coordinates found
!
!    Output, ALPHA(1:4) - scaled barycentric coordinates (if DEGEN = .FALSE.)
!    such that E = (ALPHA(1)*A + ALPHA(2)*B + ALPHA(3)*C +
!    ALPHA(4)*D)/DET where DET = 6 * (volume of tetra ABCD);
!    an ALPHA(I) may be set to 0 after tolerance test to
!    indicate that E is coplanar with a face, so sum of
!    ALPHA(I)/DET may not be 1; if the actual barycentric
!    coordinates rather than just their signs are needed,
!    modify this routine to divide ALPHA(I) by DET.
!
!    Output, logical DEGEN, TRUE iff A,B,C,D are coplanar.
!
!----------------------------------------------------------------------------------
!  REVISED for Polyde3D:
!
!    - original alpha is SCALED, such that alpha(i) = det*lambda(i),
!      i =1,..,4, so to get lambda, simply divide by det
!    - the change makes  alpha = lambda
!    - original tet has different orientation than the proposed tet
!      orientation for Polyde
!    - original baryth has clockwise orientation of A-B-C vertices and points up 
!      to vertex D (left-hand-wise), ours is anti-cw for A-B-C and points up to D 
!      right-hand-wise
!    - the change suits our tet orientation
!
!    $Revision: 1.1 $
!    $Date: 2005/07/25 09:14:05 $
!    $Author: r_abdul $
!
!----------------------------------------------------------------------------------

  implicit none

  real ( kind = 8 ) a(3)
  real ( kind = 8 ) alpha(4)
  real ( kind = 8 ) amax
  real ( kind = 8 ) b(3)
  real ( kind = 8 ) bmax
  real ( kind = 8 ) c(3)
  real ( kind = 8 ) cmax
  real ( kind = 8 ) cp1
  real ( kind = 8 ) cp2
  real ( kind = 8 ) cp3
  real ( kind = 8 ) d(3)
  real ( kind = 8 ) da(3)
  real ( kind = 8 ) db(3)
  real ( kind = 8 ) dc(3)
  real ( kind = 8 ) de(3)
  logical degen
  real ( kind = 8 ) det
  real ( kind = 8 ) dmax
  real ( kind = 8 ) e(3)
  real ( kind = 8 ) ea(3)
  real ( kind = 8 ) eb(3)
  real ( kind = 8 ) ec(3)
  real ( kind = 8 ) emax
  real ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )
  degen = .false.

  da(1:3) = a(1:3) - d(1:3)
  db(1:3) = b(1:3) - d(1:3)
  dc(1:3) = c(1:3) - d(1:3)

  amax = max ( abs(a(1)),abs(a(2)),abs(a(3)))
  bmax = max ( abs(b(1)),abs(b(2)),abs(b(3)))
  cmax = max ( abs(c(1)),abs(c(2)),abs(c(3)))
  dmax = max ( abs(d(1)),abs(d(2)),abs(d(3)))

! change the signs to suit our tet orientation
  cp1 = db(3)*dc(2) - db(2)*dc(3) 
  cp2 = db(1)*dc(3) - db(3)*dc(1) 
  cp3 = db(2)*dc(1) - db(1)*dc(2) 

  det = da(1)*cp1 + da(2)*cp2 + da(3)*cp3

  if (abs(det) <= 0.01d0*tol*max(amax,bmax,cmax,dmax)) then
    degen = .true.
    return
  end if

! keep de as it is
  de(1:3) = e(1:3) - d(1:3)

! change the signs to suit our tet orientation
  ea(1:3) = e(1:3) - a(1:3)
  eb(1:3) = e(1:3) - b(1:3)
  ec(1:3) = e(1:3) - c(1:3)

  alpha(1) = ( de(1)*cp1 + de(2)*cp2 + de(3)*cp3 )/det

! change the signs to suit our tet orientation
  cp1 = da(3)*de(2) - da(2)*de(3)
  cp2 = da(1)*de(3) - da(3)*de(1)
  cp3 = da(2)*de(1) - da(1)*de(2)

  alpha(2) = ( dc(1)*cp1 + dc(2)*cp2 + dc(3)*cp3 )/det
  alpha(3) = ( db(1)*cp1 + db(2)*cp2 + db(3)*cp3 )/det
  alpha(4) = ( ea(1)*(eb(2)*ec(3) - eb(3)*ec(2)) + ea(2)*(eb(3)*ec(1) &
    - eb(1)*ec(3)) + ea(3)*(eb(1)*ec(2) - eb(2)*ec(1)) )/det

  if (det < 0.0d0) then
    alpha(1) = -alpha(1)
    alpha(2) = -alpha(2)
    alpha(4) = -alpha(4)
  else
    alpha(3) = -alpha(3)
  end if

  emax = max(abs(e(1)),abs(e(2)),abs(e(3)))

  if ( abs(alpha(1)) <= tol*max(bmax,cmax,dmax,emax)) then
    alpha(1) = 0.0d0
  end if

  if (abs(alpha(2)) <= tol*max(amax,cmax,dmax,emax)) then
    alpha(2) = 0.0d0
  end if

  if (abs(alpha(3)) <= tol*max(amax,bmax,dmax,emax)) then
    alpha(3) = 0.0d0
  end if

  if (abs(alpha(4)) <= tol*max(amax,bmax,cmax,emax)) then
    alpha(4) = 0.0d0
  end if

  return
end subroutine baryth
