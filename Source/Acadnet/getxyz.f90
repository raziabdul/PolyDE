      subroutine getxyz(x,y,z)
      use feminterface, only: readky, unread
      use femtypes
      implicit none
      real (DP) :: x, y, z
      intent (out) :: x, y, z
!
!    $Revision: 1.5 $
!    $Date: 2008/12/22 15:25:14 $
!    $Author: m_kasper $
!
!  read x, y, z coordinates
!
!  Output:
!         x,y,z  Koordinatenwerte
!
!  local variables
      integer (I4B) :: int, key, line
      real (DP) :: float
      character (len=1000) :: string
!
100   call readky(key,string,float,int,line)
      select case (key)
      case(10)
!  x - punkt
        x=float
      case(20)
!  y - punkt
        y=float
      case(30)
!  z - punkt
        z=float
      case(9, 0)
!  ende
        call unread
        return
      case default
        print*,'*ERROR: syntax error'
      end select
      goto 100
      end subroutine getxyz
