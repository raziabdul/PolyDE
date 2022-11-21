      subroutine vertex(x2,y2,z2,ausbuc,scalfk,layanz,                  &
     &  laytxt,layrb,layer,txtlen)
      use feminterface, only: readky, ftlay, unread
      use femtypes
      implicit none
      integer (I4B) :: layanz, layer, txtlen
      integer (I4B), pointer :: layrb(:,:)
      real (DP) :: x2, y2, z2, ausbuc, scalfk
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: scalfk, txtlen
      intent (out) :: x2, y2, z2, ausbuc, layer
      intent (inout) :: layanz, laytxt, layrb
!
!    $Revision: 1.6 $
!    $Date: 2008/12/22 15:46:17 $
!    $Author: m_kasper $
!
!  Lesen des Entities: VERTEX bis ein neues Entity beginnt und
!  dessen Rueckuebergabe an Polyline
!
!  Input:
!    scalfk    Skalierungsfaktor
!    layanz    Anzahl der Layer
!    laytxt    Text-Bezeichnung der Layer
!    layrb     type of boundary condition for branches on that layer (and for each nature)
!              0    = Dirichlet
!              1-99 = User Dirichlet
!              200  = general Neumann
!              300  = innner contour
!              400  = non-visble contour
!              1000 = comment layer
!              1001 = not recognized layer, treated as a comment layer
!    txtlen    defined length for strings txtgeb
!  Output:
!    x2        x-Koordinaten des Knotens
!    y2        y-Koordinaten des Knotens
!    z2        z-Koordinaten des Knotens
!    ausbuc    Ausbuchtung
!    layer     Layer
!
!  local variables
      integer (I4B) :: int, key, lineno
      real (DP) :: float
      character (len=txtlen) :: string
!
      layer=3
      ausbuc=0._DP
!
100   call readky(key,string,float,int,lineno)
      select case (key)
      case(8)
!    - Layer
!  layer is a text - standard annotation is 0,1,2..
        call ftlay(string,layanz,laytxt,layrb,layer,txtlen)
      case(10)
!  x - location point
        x2=float*scalfk
      case(20)
!  y - location point
        y2=float*scalfk
      case(30)
!  z - location point
        z2=float*scalfk
      case(42)
!  bulge 
        ausbuc=float
      case(0)
!  beginn eines neuen Entities => Abschluss
        call unread
        return
!
!  The following group codes are not evaluated
      case (39)
!  thickness (Object heigh)
      case(40)
!  starting width
      case(41)
!  ending width
      case(50)
!  curve fit tangent direction
      case(70)
!  vertex flags
      case(71)
!  polygon mesh vertex index
      case(72)
!  polygon mesh vertex index
      case(73)
!  polygon mesh vertex index
      case(74)
!  The following group codes apply to all entity-type graphical object
      case (5)
!  handle
      case (6)
!  linetypename
      case (38)
!  entity's elevation
      case (48)
!  linetype scale
      case (60)
!  object visibility 0 = visible
      case (62)
!  color number
      case (67)
!  model or paper space
      case (92)
!  number of bytes in the proxy entity graphics
      case (100)
!  subclass Marker  Acdbxxxx
      case (102)
!  start/ end of application defined group
      case (310)
!  proxy entity graphics data
      case (330)
!  soft-pointer ID/ handle to owner directory
      case (360)
!  hard-owner ID/ handle to owner directory
      case (370)
!  linewidth enum value
      case (410)
!  APP: layout tab name
      case (420)
!  a 24 bit color value
      case (430)
!  color name
      case (440)
!  transparency value
      case default
      end select
      goto 100
      end subroutine vertex
