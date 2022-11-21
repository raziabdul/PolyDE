      subroutine laytab(layanz,laytxt,layrb,txtlen)
      use feminterface, only: ftlay, readky
      use femtypes
      implicit none
      integer (I4B) :: layanz, txtlen
      integer (I4B), pointer :: layrb(:,:)
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: txtlen
      intent (out) :: layanz  ! org: layrb, laytxt
      intent (inout) :: layrb, laytxt
!
!    $Revision: 1.6 $
!    $Date: 2008/12/22 15:57:00 $
!    $Author: m_kasper $
!
!  read the layer table
!
!  Input:
!    txtlen    defined length for strings txtgeb
!  Output:
!    layanz    Anzahl der gelesenen Layer
!    laytxt    Text-bezeichnung der Layer
!    layrb     type of boundary condition for branches on that layer (and for each nature)
!              0    = Dirichlet
!              1-99 = User Dirichlet
!              200  = general Neumann
!              300  = innner contour
!              400  = non-visible contour
!              1000 = comment layer
!              1001 = not recognized layer, treated as a comment layer
!
!  local variables
      integer (I4B) :: int, key, lineno
      real (DP) :: float
      character (len=txtlen) :: string
      integer (I4B) :: layer
!
      call readky(key,string,float,int,lineno)
      if (key.eq.70) then
      else
!        print*,'*ERROR: syntax error in line', lineno
      end if
      layanz=0
!
100   call readky(key,string,float,int,lineno)
      select case (key)
      case (0)
        if (string .eq. 'ENDTAB') then 
          return
        else if (string .eq. 'LAYER') then
!  neuer layer
        end if
      case (2)
!  layer name      
        call ftlay(string,layanz,laytxt,layrb,layer,txtlen)
!
!  The following group codes are not evaluated
      case (5)
!  handle
      case (6)
!  linetype name
      case (62)
!  color number
      case (70)
!  flags
      case (100)
!  subclass mask
      case (290)
!  plotting flag
      case (370)
!  lineweight enum value
      case (390)
!  hard pointer ID/ handle of PLotStyleName object
!
!  The following group codes apply to all entity-type graphical object
      case (102)
!    - start/ end of application defined group
      case (330)
!    - soft-pointer ID/ handle to owner directory
      case (360)
!    - hard-owner ID/ handle to owner directory
!
!  The following group codes apply to BLOCK_RECORD symbol table entries
      case (1000)
!    - Xdata string data DesignCenter Data (optional)
      case (1001)
!    - Xdata application name ACAD (optional)
!
!  The following group codes apply to LAYER symbol table entries
      case (347)
!  Hard-pointer ID/handle to Material object
!
      case default
        write(*,1000) ' *WARNING: Group code(ignored):',key,            &
     &        ' unknown entry in the layer table at line:',lineno
1000    format(a,i5,a,i7)
      end select
      goto 100
      end subroutine laytab
