      subroutine mtext(accur,scalfk,anzgeb,xgeb,ygeb,txtgeb,layanz,     &
     &  lplayt,laytxt,layrb,txtlen)
      use feminterface, only: ftlay, readky, putnod, unread, reallocate
      use femtypes
      implicit none
      integer (I4B) :: anzgeb, txtlen, layanz
      integer (I4B), pointer :: layrb(:,:), lplayt(:)
      real (DP) :: accur, scalfk
      real (DP), pointer :: xgeb(:), ygeb(:)
      character (len=txtlen), pointer :: txtgeb(:), laytxt(:)
      intent (in) :: accur, scalfk, txtlen
      intent (out) :: lplayt
      intent (inout) :: xgeb, ygeb, txtgeb
      intent (inout) :: anzgeb, layanz, laytxt, layrb
!
!    $Revision: 1.6 $
!    $Date: 2008/12/22 15:58:51 $
!    $Author: m_kasper $
!
!  read and return the entity MTEXT until the beninning of a new entity
!
!  Input:
!    accur     Genauigkeitschranke
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
!    xgeb      Einfuegungspunkt der Texte
!    ygeb
!    txtgeb    Textwert der Texte
!    lplayt    Pointer auf layer fuer die Texte
!    txtlen    defined length for strings txtgeb
!  In-/ Output: 
!    anzgeb     Anzahl der gelesenen Texte 
!
!  local variables
      integer (I4B) :: int, key, lineno, node, layer, txtcont, pos
      real (DP) :: float, x1, y1, z1
      logical :: pfound
      character (len=txtlen) :: string, txtwrt
!
      txtcont=0
      txtwrt=''
100   call readky(key,string,float,int,lineno)
      select case (key)
      case(8)
!    - Layer
!  layer is a text - standard annotation is 0,1,2..
        call ftlay(string,layanz,laytxt,layrb,layer,txtlen)
      case(10)
!  x - first alignment point
        x1=float*scalfk
      case(20)
!  y - first alignment point
        y1=float*scalfk
      case(30)
!  z - first alignment point
        z1=float*scalfk
      case(1)
!    - last chunk of less than 250 characters
        txtwrt(txtcont*250+1:(txtcont+1)*250)=string(1:250)
      case(3)
!    - 250 character chunk 
        txtwrt(txtcont*250+1:max((txtcont+1)*250,txtlen))=string(1:250)
        txtcont=txtcont+1
      case(0)
!  Beginn eines neuen Entities => Abschluss
!      print*,'text: ',x1,y1,txtwrt
!
!  replace new paragraph: "\P" by blanks
        pfound=.true.
        do while (pfound)
          pos=index(txtwrt,'\P')
          if (pos .ne. 0) then
            txtwrt(pos:pos+1)='  '
          else 
            pfound=.false.
          end if
        end do
        if (layrb(layer,1).lt.1000) then
!  reallocate arrays if necessary
          if (anzgeb+1.gt.size(xgeb))                                   &
     &        xgeb=>reallocate(xgeb,2*size(xgeb)+1)
          if (anzgeb+1.gt.size(ygeb))                                   &
     &        ygeb=>reallocate(ygeb,2*size(ygeb)+1)
!
          call putnod(x1,y1,accur,node,xgeb,ygeb,anzgeb)
!  reallocate arrays if necessary
          if (node.gt.size(txtgeb))                                     &
     &        txtgeb=>reallocate(txtgeb,2*size(txtgeb)+1,txtlen)
          txtgeb(node)=txtwrt(1:len_trim(txtwrt))
          if (node.gt.size(lplayt))                                     &
     &        lplayt=>reallocate(lplayt,2*size(lplayt)+1)
          lplayt(node)=layer
        end if
!
!  Zurueckstellen der letzten Eingabe
        call unread
        return
!

!  The following group codes are not evaluated
      case(7)
!  textstyle name
      case(11)
!  x  - value of x-axis direction vector
      case(21)
!  y  - value of x-axis direction vector
      case(31)
!  z  - value of x-axis direction vector
      case (39)
!  thickness (Object height)
      case(40)
!  text height
      case(41)
!  reference rectangle width
      case(42)
!  horizontal width
      case(43)
!  vertical height
      case(44)
!  line spacing factor
      case(45)
!  Fill box scale (optional):
      case(46)
!  Defined annotation height
      case(48)
!  Column width
      case(49)
!  Column gutter
      case(50)
!  rotation angle
      case(63)
!  Background color (if color index number)
      case(71)
!  attachment point
      case(72)
!  drawing direction
      case(73)
!  line spacing style
      case(75)
!  Column type
      case(76)
!  Column count
      case(78)
!  Column Flow Reversed
      case(79)
!  Column Autoheight
      case(90)
!  Background fill setting:
      case (210)
!  x - extrusion direction
      case (220)
!  y - extrusion direction
      case (230)
!  z - extrusion direction
      case(421 : 429)
!  Background color (if RGB color)
      case(431 : 439)
!  Background color (if color name)
      case(441)
!  Transparency of background fill color (not implemented)
!
!  The following group codes apply to all entity-type graphical object
      case (5)
!  handle
      case (6)
!  linetypename
      case (38)
!  entity's elevation
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
     &        ' unknown syntax in mtext at line:',lineno
1000    format(a,i5,a,i7)
      end select
      goto 100
      end subroutine mtext
