function angle3 ( u, v, rtolsq )

!******************************************************************************
!
!! ANGLE3 computes the size of a plane angle in 3D.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Purpose: 
!
!    Compute angle in range [0,PI] between 3D vectors U and V.
!
!  Parameters:
!
!    Input, U(1:3), V(1:3) - vectors.
!
!    Input, RTOLSQ - relative tolerance used to detect 0 vector based on
!    square of Euclidean length.
!
!    Output, ANGLE3 - angle between 2 vectors in range [0,PI]
!    If U or V is the 0 vector, ANGLE3 = PI is returned.
!
  implicit none

  real ( kind = 8 ) angle3
  real ( kind = 8 ) d_pi
  real ( kind = 8 ) dotp
  real ( kind = 8 ) lu
  real ( kind = 8 ) lv
  real ( kind = 8 ) rtolsq
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) u(3)
  real ( kind = 8 ) v(3)

  tol = 100.0D+00 * epsilon ( tol )

  dotp = dot_product ( u(1:3), v(1:3) )

  lu = dot_product ( u(1:3), u(1:3) )

  lv = dot_product ( v(1:3), v(1:3) )

  if (lu > rtolsq .and. lv > rtolsq) then
    t = dotp / sqrt(lu*lv)
    if (abs(t) > 1.0d0 - tol ) then
      t = sign(1.0d0,t)
    end if
    angle3 = acos ( t )
  else
    angle3 = d_pi ( )
  end if

  return
end
