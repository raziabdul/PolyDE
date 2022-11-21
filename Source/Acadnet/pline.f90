      subroutine pline(accur,fakt1,scalfk,laytxt,layrb,layanz,xpoint,   &
     &  ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,okall,txtlen)
      use feminterface, only: readky, ftlay, putlin, vertex, unread
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, txtlen, layanz
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:), layrb(:,:)
      real (DP) :: accur, fakt1, scalfk
      real (DP), pointer :: xpoint(:), ypoint(:), length(:)
      logical :: okall
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: accur, fakt1, scalfk, txtlen
      intent (out) :: xpoint, ypoint, zki, lzrb, zpz, length, okall
      intent (inout) :: laytxt, layrb, layanz, anzzwg, anzknt
!
!    $Revision: 1.6 $
!    $Date: 2008/12/22 16:00:59 $
!    $Author: m_kasper $
!
!  read and return the entity POLYLINE until the beninning of a new entity
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
!    okall     .true.   wenn die Polyline vollstaendig aufgenommen werden konnte 
!              .false.  wenn ein Element nicht in die Liste aufgenommen wurde
!  In-/ Output: 
!    anzknt    ermittelte Anzahl der Knoten
!    anzzwg    ermittelte Anzahl der Zweige
!
!  local variables
      integer (I4B) :: int, key, lineno, layer, close
      real (DP) :: float, x1, y1, z1, x2, y2, z2, xanf, yanf, zanf, lenall
      real (DP) :: ausbu1, ausbu2, xc, yc, radius, phi1, phi2, dummy
      real (DP), parameter ::  pid180=.017453292519943295769237_DP
      logical :: ok, firstv
      character (len=txtlen) :: string
!
      layer=3
      okall=.true.
      lenall=0._DP
      firstv=.true.
100   call readky(key,string,float,int,lineno)
      select case (key)
      case(8)
!    - Layer
!  layer is a text - standard annotation is 0,1,2..
        call ftlay(string,layanz,laytxt,layrb,layer,txtlen)
      case(10)
!  x -
        x1=float*scalfk
      case(20)
!  y -
        y1=float*scalfk
      case(30)
!  z -
        z1=float*scalfk
      case(70)
!  polyline flags bit-coded value  1=closed polyline
        close=mod(int,2)
      case(0)
!  es koennen folgen: vertex, seqend oder ein neues entity
        if (string.eq.'VERTEX') then
          call vertex(x2,y2,z2,ausbu2,scalfk,layanz,laytxt,layrb,layer,txtlen)
          if (layrb(layer,1).lt.1000 .and. (.not. firstv) ) then
            if (abs(ausbu1) .lt. tiny(1._SP) ) then
!  Gerade
              call putlin(x1,y1,x2,y2,0._DP,0._DP,0._DP,.true.,accur,   &
     &          fakt1,layer,ok,xpoint,ypoint,zki,lzrb,zpz,              &
     &          length,anzzwg,anzknt)
            else
!  Kreisbogen
              xc=(x1+x2)/2._DP+(y1-y2)*(1._DP-ausbu1**2)/(4._DP*ausbu1)
              yc=(y1+y2)/2._DP+(x2-x1)*(1._DP-ausbu1**2)/(4._DP*ausbu1)
              radius= (sqrt( (x1-xc)**2+(y1-yc)**2) +                   &
     &          sqrt( (x2-xc)**2+(y2-yc)**2) )/2._DP
              phi1=atan2(y1-yc,x1-xc)/pid180
              phi2=atan2(y2-yc,x2-xc)/pid180
              if (ausbu1.lt.0) then
                dummy=phi1
                phi1=phi2
                phi2=dummy
              end if
              call putlin(xc,yc,xc,yc,radius,phi1,phi2,.false.,accur,   &
     &          fakt1,layer,ok,xpoint,ypoint,zki,lzrb,zpz,              &
     &          length,anzzwg,anzknt)
            end if
            if (.not. ok) okall=.false.
            lenall=lenall+length(anzzwg)
          else
            xanf=x2
            yanf=y2
            zanf=z2
!
!  len ???????????????????????????????????????????
          end if
          firstv=.false.
          x1=x2
          y1=y2
          z1=z2
          ausbu1=ausbu2
        else if (string .eq. 'SEQEND') then
          if (layrb(layer,1).lt.1000 .and. close.eq.1 ) then
            x2=xanf
            y2=yanf
            z2=zanf
            if (abs(ausbu1) .lt. tiny(1._SP) ) then
!  Gerade
              call putlin(x1,y1,x2,y2,0._DP,0._DP,0._DP,.true.,accur,   &
     &          fakt1,layer,ok,xpoint,ypoint,zki,lzrb,zpz,              &
     &          length,anzzwg,anzknt)
            else
!  Kreisbogen
              xc=(x1+x2)/2._DP+(y1-y2)*(1._DP-ausbu1**2)/(4._DP*ausbu1)
              yc=(y1+y2)/2._DP+(x2-x1)*(1._DP-ausbu1**2)/(4._DP*ausbu1)
              radius= (sqrt( (x1-xc)**2+(y1-yc)**2) +                   &
     &          sqrt( (x2-xc)**2+(y2-yc)**2) )/2._DP
              phi1=atan2(y1-yc,x1-xc)/pid180
              phi2=atan2(y2-yc,x2-xc)/pid180
              if (ausbu1.lt.0._DP) then
                dummy=phi1
                phi1=phi2
                phi2=dummy
              end if
              call putlin(xc,yc,xc,yc,radius,phi1,phi2,.false.,accur,   &
     &          fakt1,layer,ok,xpoint,ypoint,zki,lzrb,zpz,              &
     &          length,anzzwg,anzknt)
            end if
            if (.not. ok) okall=.false.
            lenall=lenall+length(anzzwg)
          end if
        else
!  beginn eines neuen Entities => Abschluss
!  zurueckstellen der letzten Eingabe
          call unread
          return
        end if
!
!  The following group codes are not evaluated
      case (39)
!  thickness (Object height)
      case(40)
!  default start width
      case(41)
!  default end width
      case(66)
!  obsolet entities follow flag
      case(71)
!  polygon mesh M vertex count
      case(72)
!  polygon mesh N vertex count
      case(73)
!  smooth surface M density
      case(74)
!  smooth surface M density
      case(75)
!  curves and smooth surface type
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
!  transparency value
      case default
        write(*,1000) ' *WARNING: Group code(ignored):',key,            &
     &        ' unknown syntax in pline at line:',lineno
1000    format(a,i5,a,i7)
      end select
      goto 100
      end subroutine pline
