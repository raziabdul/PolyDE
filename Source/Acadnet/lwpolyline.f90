      subroutine lwpolyline(accur,fakt1,scalfk,laytxt,layrb,layanz,     &
     &  xpoint,ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
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
      intent (out) :: xpoint, ypoint, zki, lzrb, zpz, length, ok
      intent (inout) :: laytxt, layrb, layanz, anzzwg, anzknt
!
!    $Revision: 1.7 $
!    $Date: 2008/12/22 15:57:50 $
!    $Author: m_kasper $
!
!  read and return the entity LWPOLYLINE until the beninning of a new entity
!
!  Input:
!    accur     Genauigkeitschranke
!    fakt1     Punktdichte fuer Kreisboegen
!    scalfk    Skalierungsfaktor
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
!    anzknt    ermittelte Anzahl der Knoten
!    anzzwg    ermittelte Anzahl der Zweige
!
!  local variables
      integer (I4B) :: int, key, lineno, layer, nvert, polylflag, iact, i
      real (DP) :: float, x1, x2, y1, y2, xc, yc, radius, phi1, phi2
      real (DP) :: ausb, dummy
      real (DP), allocatable :: x(:), y(:), bulge(:)
      real (DP), parameter ::  pid180=.017453292519943295769237_DP
      character (len=txtlen) :: string
!
      layer=3
100   call readky(key,string,float,int,lineno)
      select case (key)
      case (8)
!    - Layer
!  layer is a text - standard annotation is 0,1,2..
        call ftlay(string,layanz,laytxt,layrb,layer,txtlen)
      case (10)
!  x - vertex coordinate
        iact=iact+1
        x(iact)=float*scalfk
      case (20)
!  y - vertex coordinate
        y(iact)=float*scalfk
      case(42)
!  bulge
        bulge(iact)=float
      case (90)
!  number of vertices
        nvert=int
        allocate(x(nvert),y(nvert),bulge(nvert))
        bulge=0._DP
        iact=0
      case(70)
!  polyline flag  1=closed
        polylflag=int
      case (0)
!  start of a new entity => output
        if (layrb(layer,1).lt.1000) then
          do i=1,nvert-1
            if (abs(bulge(i)) .lt. tiny(1._SP) ) then
!  a staight line 
              call putlin(x(i),y(i),x(i+1),y(i+1),0._DP,0._DP,0._DP,    &
     &          .true.,accur,fakt1,layer,ok,xpoint,ypoint,zki,lzrb,     &
     &          zpz,length,anzzwg,anzknt)
            else
!  an arc
              x1=x(i)
              y1=y(i)
              x2=x(i+1)
              y2=y(i+1)
              ausb=bulge(i)
              xc=(x1+x2)/2._DP+(y1-y2)*(1._DP-ausb**2)/(4._DP*ausb)
              yc=(y1+y2)/2._DP+(x2-x1)*(1._DP-ausb**2)/(4._DP*ausb)
              radius= (sqrt( (x1-xc)**2+(y1-yc)**2) +                   &
     &          sqrt( (x2-xc)**2+(y2-yc)**2) )/2._DP
              phi1=atan2(y1-yc,x1-xc)/pid180
              phi2=atan2(y2-yc,x2-xc)/pid180
              if (ausb.lt.0) then
                dummy=phi1
                phi1=phi2
                phi2=dummy
              end if
              call putlin(xc,yc,xc,yc,radius,phi1,phi2,.false.,accur,   &
     &          fakt1,layer,ok,xpoint,ypoint,zki,lzrb,zpz,              &
     &          length,anzzwg,anzknt)
            end if
          end do
          if (polylflag .eq. 1) then
            if (abs(bulge(i)) .lt. tiny(1._SP) ) then
!  a staight line 
              call putlin(x(nvert),y(nvert),x(1),y(1),0._DP,0._DP,0._DP,&
     &          .true.,accur,fakt1,layer,ok,xpoint,ypoint,zki,lzrb,     &
     &          zpz,length,anzzwg,anzknt)
            else
!  an arc
              x1=x(nvert)
              y1=y(nvert)
              x2=x(1)
              y2=y(1)
              ausb=bulge(nvert)
              xc=(x1+x2)/2._DP+(y1-y2)*(1._DP-ausb**2)/(4._DP*ausb)
              yc=(y1+y2)/2._DP+(x2-x1)*(1._DP-ausb**2)/(4._DP*ausb)
              radius= (sqrt( (x1-xc)**2+(y1-yc)**2) +                   &
     &          sqrt( (x2-xc)**2+(y2-yc)**2) )/2._DP
              phi1=atan2(y1-yc,x1-xc)/pid180
              phi2=atan2(y2-yc,x2-xc)/pid180
              if (ausb.lt.0) then
                dummy=phi1
                phi1=phi2
                phi2=dummy
              end if
              call putlin(xc,yc,xc,yc,radius,phi1,phi2,.false.,accur,   &
     &          fakt1,layer,ok,xpoint,ypoint,zki,lzrb,zpz,              &
     &          length,anzzwg,anzknt)
            end if
          end if
        end if
        deallocate(x,y,bulge)
!  zurueckstellen der letzten Eingabe
        call unread
        return
!
!  The following group codes are not evaluated
      case (39)
!  thickness (Object height)
      case (40)
!  starting width
      case (41)
!  end width
      case(43)
!  constant width
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
     &        ' unknown syntax in lwpolyline at line:',lineno
1000    format(a,i5,a,i7)
      end select
      goto 100
      end subroutine lwpolyline
