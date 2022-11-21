      pure FUNCTION tetmeanratio(vert)
      USE femtypes
      USE feminterface3d, ONLY: tetvolume
      IMPLICIT NONE
      REAL   (DP) :: tetmeanratio, vert(3,4)
      INTENT (IN) :: vert
!
!------------------------------------------------------------------------------
!    $Revision: 1.2 $
!    $Date: 2014/04/28 09:54:39 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!
!  input:
!  vert (real) - 3x4 array of vertices coordinates:
!                vert(i,j), where i=x,y,z; j=V1,V2,V3,V4
!
!  output:
!  tetmeanratio (real) - The Mean ratio of a Tetrahedron as an Indicator for it's shape
!                        (GEOMPACK's  eta  quality criterion)
!
! The mean ratio is calculated by (12*(3*volume)**2/3)/(SUM(lengths)**2), based on:
!   A. Liu, B. Joe, On the shape of tetrahedra from bisection, Math. Comp. 63 (1994) 141154.
!   I. Kossaczký, A recursive approach to local mesh refinement in two and three dimensions, J. Comput. Appl. Math. 55 (1994) 275--88
!   A comparative study between some bisection based partitions in 3D, Miguel A. Padrón, José P. Suárezb, Ángel Plaza (2005)
!
! base   : vert(:,1:3)
! height : vert(:,4)
!
      REAL(DP) :: le(6), volume
   !__
   !
   ! STRUCTURE:
   !          / L12 \
   !         |  L13  |
   !         |  L14  |
   ! le(6) = |  L23  |
   !         |  L24  |
   !          \ L34 /
   !     
   !__
   !
   ! 
      volume = abs(tetvolume(vert))
      
      le = (/ (vert(1,2)-vert(1,1))**2._DP + (vert(2,2)-vert(2,1))**2._DP + (vert(3,2)-vert(3,1))**2._DP,  &
     &        (vert(1,3)-vert(1,1))**2._DP + (vert(2,3)-vert(2,1))**2._DP + (vert(3,3)-vert(3,1))**2._DP,  &
     &        (vert(1,4)-vert(1,1))**2._DP + (vert(2,4)-vert(2,1))**2._DP + (vert(3,4)-vert(3,1))**2._DP,  &
     &        (vert(1,3)-vert(1,2))**2._DP + (vert(2,3)-vert(2,2))**2._DP + (vert(3,3)-vert(3,2))**2._DP,  &
     &        (vert(1,4)-vert(1,2))**2._DP + (vert(2,4)-vert(2,2))**2._DP + (vert(3,4)-vert(3,2))**2._DP,  &
     &        (vert(1,4)-vert(1,3))**2._DP + (vert(2,4)-vert(2,3))**2._DP + (vert(3,4)-vert(3,3))**2._DP  /)
      
      tetmeanratio = 12._DP*((3._DP*volume)**(2._DP/3._DP))/(le(1)+le(2)+le(3)+le(4)+le(5)+le(6))
   !__
      RETURN
      END FUNCTION tetmeanratio
