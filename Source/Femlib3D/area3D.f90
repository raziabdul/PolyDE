      pure function area3d(vert)
      use femtypes
      implicit none
      real (DP) :: area3d, vert(3,3)
      intent(in) :: vert
!
!------------------------------------------------------------------------------
!    $Revision: 1.6 $
!    $Date: 2014/08/22 10:55:23 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  area3d - calculates area of triangle in 3D space
!  area3d always return a positive value independent of the actual orientation of the triangle
!
!  input:
!  vert (real) - 3x3 array of vertices coordinates:
!                vert(i,j), where i=x,y,z; j=V1,V2,V3
!
!  output:
!  area3d (real) - area
!
! Area is calculated by 1/2 * base * height, where
! base = V2-V1
! height = distance of V3 from base, which simplifies to
! Area = | ( V2-V1 ) x ( V1-V3 ) |/2
!
!
      real (DP) :: x12, x31, y12, y31, z12, z31

      x12 = vert(1,2) - vert(1,1)
      x31 = vert(1,1) - vert(1,3)

      y12 = vert(2,2) - vert(2,1)
      y31 = vert(2,1) - vert(2,3)

      z12 = vert(3,2) - vert(3,1)
      z31 = vert(3,1) - vert(3,3)
      
      area3d = sqrt( ( y31 * z12 - z31 * y12 )**2 & 
&                  + ( z31 * x12 - x31 * z12 )**2 &
&                  + ( x31 * y12 - y31 * x12 )**2 ) / 2._DP

      return
      end function area3d
