      pure function tetvolume( vert )
      use femtypes
      implicit none
      real (DP) :: tetvolume, vert(3,4)
      intent (in) :: vert
!
!    $Revision: 1.4 $
!    $Date: 2014/04/28 09:55:34 $
!    $Author: m_kasper $
!
!  TETVOLUME computes the volume of a tetrahedron.
!
!
!    input
!        vert       vertex coordinats of the terahedron
!
!    output
!        tetvolume  terahedron volume
!
!  The volume has positive sign for a left-handed tetrahedron, 
!  for which the face normal determined by the right-hand rule 
!  following the sequence of vertices (1,2,3) points away form vertex 4.  
!  i.e. vertices (1,2,3) form a counterclockwise sequence 
!  when viewed form the side opposite of vertex 4, 
!  Otherwise volume will be negative
!
      real(DP) :: delx(3), dely(3), delz(3)
!
      delx(1:3) = vert(1,1:3) - vert(1,4)
      dely(1:3) = vert(2,1:3) - vert(2,4)
      delz(1:3) = vert(3,1:3) - vert(3,4)

      tetvolume = delx(1) * ( dely(2)*delz(3) - delz(2)*dely(3) )       &
     &          + delx(2) * ( dely(3)*delz(1) - delz(3)*dely(1) )       &
     &          + delx(3) * ( dely(1)*delz(2) - delz(1)*dely(2) )
      tetvolume = tetvolume / 6._DP
      return
      end function tetvolume
