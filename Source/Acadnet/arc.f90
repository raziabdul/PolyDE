      subroutine arc(accur,fakt1,scalfk,laytxt,layrb,layanz,xpoint,     &
     &  ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
      use feminterface, only: readky, ftlay, putlin, unread
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, txtlen, layanz
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:), layrb(:,:)
      real (DP) :: accur, fakt1, scalfk
      real (DP), pointer :: xpoint(:), ypoint(:), length(:)
      logical :: ok
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: accur, fakt1, scalfk
      intent (inout) :: xpoint, ypoint, zki, lzrb, zpz, length, ok ! Change to inout since old intent gives runtime error to OneAPI
      intent (inout) :: laytxt, layrb, layanz, anzzwg, anzknt
!
!    $Revision: 1.6 $
!    $Date: 2008/12/22 15:55:17 $
!    $Author: m_kasper $
!
!  read and return the entity ARC
!
!  Input:
!    accur     Genauigkeitschranke
!    fakt1     Punktdichte fuer Kreisboegen
!    scalfk    Skalierungsfaktor
!    laytxt    Text-Bezeichnung der Layer
!    layrb     type of boundary condition for branches on that layer (and for each nature)
!              0    = Dirichlet
!              1-99 = User Dirichlet
!              200  = general Neumann
!              300  = innner contour
!              400  = non-visble contour
!              1000 = comment layer
!              1001 = not recognized layer, treated as a comment layer
!    layanz    Anzahl der Layer
!    txtlen    defined length for strings txtgeb
!  Output:
!    xpoint    x-Koordinaten der Knoten
!    ypoint    y-Koordinaten der Knoten
!    zki       Zweig-Knoten Information
!    lzrb      Randbedingungen der Zweige Pointer auf layer
!    zpz       Zahl der Knoten auf dem Zweig
!    length    Laenge der Zweige
!    ok        .true.
!              .false. Wenn das Element nicht in die Liste aufgenommen wurde
!  In-/ Output: 
!    anzknt    ermittelte Anzahl der Knoten
!    anzzwg    ermittelte Anzahl der Zweige
!
!  local variables
      integer (I4B) :: int, key, lineno, layer
      real (DP) :: float, x1, y1, z1, phi1, phi2, r
      character (len=txtlen) :: string
!
      layer=3
100   call readky(key,string,float,int,lineno)
      select case (key)
      case(8)
!    - Layer
!  layer is a text - standard annotation is 0,1,2..
        call ftlay(string,layanz,laytxt,layrb,layer,txtlen)
      case(10)
!  x - center point
        x1=float*scalfk
      case(20)
!  y - center point
        y1=float*scalfk
      case(40)
!  r - radius
        r=float*scalfk
      case(50)
!  start angle
        phi1=float
      case(51)
!  end angle
        phi2=float
      case(30)
!  z - center point
        z1=float*scalfk
      case(0)
!  beginn eines neuen Entities => Abschluss
        if (layrb(layer,1).lt.1000) then
            print *, 'IN ARC. CALL PUTLIN>>>>'

            call putlin(x1,y1,0._DP,0._DP,r,phi1,phi2,.false.,accur,      &
     &      fakt1,layer,ok,xpoint,ypoint,zki,lzrb,zpz,length,           &
     &      anzzwg,anzknt)
        end if
!  zurueckstellen der letzten Eingabe
        call unread
        return
!
!  The following group codes are not evaluated
      case (39)
!  thickness (Object height)
      case (210)
!  x - extrusion direction
      case (220)
!  y - extrusion direction
      case (230)
!  z - extrusion direction
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
      case (390)
!  Hard-pointer ID/handle to the plot style object
      case (410)
!  APP: layout tab name
      case (420)
!  a 24 bit color value
      case (430)
!  color name
      case (440)
!    - transparency value
      case default
        write(*,1000) ' *WARNING: Group code(ignored):',key,' unknown syntax in arc at line:',lineno
1000    format(a,i5,a,i7)
      end select
      goto 100
      end subroutine arc
