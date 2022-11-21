      FUNCTION tetangles(vert)
      USE femtypes
      USE feminterface3d, ONLY: tetvolume
      IMPLICIT NONE
      REAL   (DP) :: tetangles(4), vert(3,4)
      INTENT (IN) :: vert
!
!------------------------------------------------------------------------------
!    $Revision: 1.3 $
!    $Date: 2014/05/12 11:36:37 $
!    $Author: juryanatzki $
!------------------------------------------------------------------------------
!
!  mintetangle - calculates area of 4 unit sphere projections of unit triangles
!                originated from each of the tetraheder's verticies and choses
!                the smallest. This can be used as a Quality criterion for tetrahedra
!
!  input:
!  vert (real) - 3x4 array of vertices coordinates:
!                vert(i,j), where i=x,y,z; j=V1,V2,V3,V4
!
!  output:
!  mintetangle (real) - Minimum solid Angle of tetrahedron vert(3,4)
!
! The sinus of the mintetangle devided by 2 is calculated by sin(Teta/2) = (12*Volume(vert))/(sqr(~)), where
! base   : vert(:,1:3)
! height : vert(:,4)
! Area   
!
      REAL(DP) :: le(6), volume
      INTEGER (I4B),   PARAMETER :: ind(3,3,4) = RESHAPE((/ 3,1,2,1,2,3,5,4,6, &
      &                                                     1,4,5,4,5,1,2,6,3, &
      &                                                     2,4,2,4,6,6,1,5,3, &
      &                                                     3,5,3,5,6,6,1,4,2 /),(/3,3,4/))
      INTEGER  :: n
      
    volume = ABS(tetvolume(vert))
!
! STRUCTURE:
!          / L12 \
!         |  L13  |
!         |  L14  |
! le(6) = |  L23  |
!         |  L24  |
!          \ L34 /
!
        le = (/ SQRT((vert(1,2)-vert(1,1))**2 + (vert(2,2)-vert(2,1))**2 + (vert(3,2)-vert(3,1))**2),  &
        &       SQRT((vert(1,3)-vert(1,1))**2 + (vert(2,3)-vert(2,1))**2 + (vert(3,3)-vert(3,1))**2),  &
        &       SQRT((vert(1,4)-vert(1,1))**2 + (vert(2,4)-vert(2,1))**2 + (vert(3,4)-vert(3,1))**2),  &
        &       SQRT((vert(1,3)-vert(1,2))**2 + (vert(2,3)-vert(2,2))**2 + (vert(3,3)-vert(3,2))**2),  &
        &       SQRT((vert(1,4)-vert(1,2))**2 + (vert(2,4)-vert(2,2))**2 + (vert(3,4)-vert(3,2))**2),  &
        &       SQRT((vert(1,4)-vert(1,3))**2 + (vert(2,4)-vert(2,3))**2 + (vert(3,4)-vert(3,3))**2)  /)
!
    DO n = 1, 4        
        tetangles(n) = (12*volume)/( SQRT((le(ind(1,1,n))+le(ind(1,2,n))+le(ind(1,3,n)))*(le(ind(1,1,n))+le(ind(1,2,n))-le(ind(1,3,n)))* &
        &                                 (le(ind(2,1,n))+le(ind(2,2,n))+le(ind(2,3,n)))*(le(ind(2,1,n))+le(ind(2,2,n))-le(ind(2,3,n)))* &
        &                                 (le(ind(3,1,n))+le(ind(3,2,n))+le(ind(3,3,n)))*(le(ind(3,1,n))+le(ind(3,2,n))-le(ind(3,3,n))) ))
    END DO
      RETURN
      END FUNCTION tetangles
