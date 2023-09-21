      subroutine entity(accur,fakt1,scalfk,xpoint,ypoint,zki,lzrb,zpz,  &
     &  xgeb,ygeb,txtgeb,length,anzzwg,anzknt,anzgeb,laytxt,layrb,      &
     &  lplayt,layanz,snapof,txtlen)
      use feminterface, only: readky, line, arc, text, mtext
      use feminterface, only: lwpolyline, circ, pline
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, anzgeb, snapof, layanz, txtlen
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:), lplayt(:), layrb(:,:)
      real (DP) :: accur, fakt1, scalfk
      real (DP), pointer :: xgeb(:), ygeb(:), xpoint(:), ypoint(:), length(:)
      character (len=txtlen), pointer :: txtgeb(:), laytxt(:)
      intent (in) :: accur, fakt1, scalfk, txtlen
!      intent (out) :: zki, lzrb, zpz  ! org
      intent (inout) :: zki, lzrb, zpz
!      intent (out) :: anzzwg, anzknt, anzgeb, lplayt, snapof ! org
      intent (inout) :: anzzwg, anzknt, anzgeb, lplayt, snapof ! test

      !      intent (out) :: xpoint, ypoint, length ! org 
      intent (inout) :: laytxt, layrb, layanz, txtgeb, xgeb, ygeb
      intent (inout) :: xpoint, ypoint, length ! test
!
!    $Revision: 1.7 $
!    $Date: 2009/04/29 12:41:13 $
!    $Author: m_kasper $
!
!  read the entity section
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
!    lzrb      Randbedingungen der Zweige, Pointer auf layer
!    zpz       Zahl der Knoten auf dem Zweig
!    xgeb      Einfuegungspunkt der Texte
!    ygeb
!    txtgeb    Textwert der Texte
!    length    Laenge der Zweige
!    anzknt    ermittelte Anzahl der Knoten
!    anzzwg    ermittelte Anzahl der Zweige
!    anzgeb    Anzahl der gelesenen Texte
!    lplayt    Pointer auf layer fuer die Texte
!    snapof    Anzahl der Zeichnungselemente die eliminiert wurden,
!                weil ihre Ausdehnung kleiner als accur ist
!                oder weil Linien mehrfach (mit gleichen Endpunkten)
!                vorhanden waren
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
!
!  local variables
      integer (I4B) :: int, key, lineno
      real (DP) :: float
      logical :: ok
      character (len=txtlen) string
!
      anzzwg=0
      anzknt=0
      anzgeb=0
      snapof=0
!
      call readky(key,string,float,int,lineno)
      if (key.ne.0)                                                     &
     &  print*,'*ERROR syntax error in line ',lineno
!
200   select case (string)
      case('LINE')
!  lines
        call line(accur,fakt1,scalfk,laytxt,layrb,layanz,               &
     &    xpoint,ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
        if ( .not. ok) snapof=snapof+1
      case('ARC')
!  arc
      print *, 'ENTITY. CALL ARC>>>>>'
        call arc(accur,fakt1,scalfk,laytxt,layrb,layanz,                &
     &    xpoint,ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
        if ( .not. ok) snapof=snapof+1
      case('CIRCLE')
!  circle
        call circ(accur,fakt1,scalfk,laytxt,layrb,layanz,               &
     &    xpoint,ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
        if ( .not. ok) snapof=snapof+1
      case('LWPOLYLINE')
!  simple staight polylines
        call lwpolyline(accur,fakt1,scalfk,laytxt,layrb,layanz,         &
     &    xpoint,ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
        if ( .not. ok) snapof=snapof+1
      case('POLYLINE')
!  polyline
        call pline(accur,fakt1,scalfk,laytxt,layrb,layanz,              &
     &    xpoint,ypoint,zki,lzrb,zpz,length,anzzwg,anzknt,ok,txtlen)
        if ( .not. ok) snapof=snapof+1
      case('TEXT')
!  Text
        call text(accur,scalfk,anzgeb,xgeb,ygeb,txtgeb,layanz,          &
     &    lplayt,laytxt,layrb,txtlen)
      case('MTEXT')
!  MText Text of one or more lines
        call mtext(accur,scalfk,anzgeb,xgeb,ygeb,txtgeb,layanz,         &
     &    lplayt,laytxt,layrb,txtlen)
      case('ENDSEC')
        return
!  entities which are ignored
      case('POINT')
      case('VIEWPORT')
!  just ignored
      case default
!  ignore everyting else, however inform
        print*,'*WARNING the entity >> ',string(1:len_trim(string)),    &
     &    ' << in line',lineno,' is not supported (ignored)'
        print*
      end select
!  Ueberlesen der weiteren Eingabe bis zum naechsten Entity
300   call readky(key,string,float,int,lineno)
      if (key .ne. 0)  goto 300
      goto 200
      end subroutine entity
