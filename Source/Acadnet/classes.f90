      subroutine classes
      use feminterface, only: readky
      use femtypes
      implicit none
!
!    $Revision: 1.5 $
!    $Date: 2008/12/22 15:24:29 $
!    $Author: m_kasper $
!
!  Ueberlesen der Section: Classes
      integer (I4B) int,key,line
      real (DP) float
      character (len=1000) string
!
100   call readky(key,string,float,int,line)
      if ((key.eq.0) .and. (string .eq. 'ENDSEC')) return
      goto 100
      end subroutine classes
