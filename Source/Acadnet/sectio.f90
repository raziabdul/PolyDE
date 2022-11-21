      subroutine sectio(sction)
      use feminterface, only: readky
      use femtypes
      implicit none
      integer (I4B) :: sction
      intent (out) :: sction
!
!    $Revision: 1.5 $
!    $Date: 2008/12/22 15:26:59 $
!    $Author: m_kasper $
!
!  read the header of a section and returns its type
!
!  Output:
!             1:   Section Header
!        /    2:   Section Tables
!             3:   Section Blocks
!  sction  -  4:   Section Entities
!             5:   Section Classes
!        \    6:   Section Objects
!             7:   Section Thumbnailimage
!            10:   Section End of File - Ende
!            99:   Unbekannte Section
!
!  local variables
      integer (I4B) :: int, key, lineno
      real (DP) :: float
      character (len=1000) :: string
!
      call readky(key,string,float,int,lineno)
      if (key.ne.0) print*,'*ERROR: syntax error at line ',lineno
      if (string.eq.'EOF')  then
        sction = 10
!  Ende der Zeichnung
        print*,lineno,' lines read'
        return
      end if
      if (string.ne.'SECTION')                                          &
     &  print*,'*ERROR: syntax error in line ',lineno
      call readky(key,string,float,int,lineno)
      if (key.eq.2) then
        select case (string)
        case('HEADER')
          sction=1
        case('TABLES')
          sction=2
        case('BLOCKS')
          sction=3
        case('ENTITIES')
          sction=4
        case('CLASSES')
          sction=5
        case('OBJECTS')
          sction=6
        case('THUMBNAILMAGE')
          sction=7
        case default
          print*,'*ERROR: unknown section found in line:',lineno
          sction=99
        end select
      else
        print*,'*ERROR: syntax error in line ',lineno-1,                      &
     &    ' Beginning of Section'
      end if
      end subroutine sectio