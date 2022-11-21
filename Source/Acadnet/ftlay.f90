      subroutine ftlay(string,layanz,laytxt,layrb,layer,txtlen)
      use feminterface, only: low2hi, string2number, reallocate
      use femtypes
      implicit none
      integer (I4B) :: layanz, layer, txtlen
      integer (I4B), pointer :: layrb(:,:)
      character (len=txtlen) :: string
      character (len=txtlen), pointer :: laytxt(:)
      intent (in) :: txtlen
      intent (out) :: layer
      intent (inout) :: string, layanz, laytxt, layrb
!
!    $Revision: 1.7 $
!    $Date: 2009/04/29 12:41:59 $
!    $Author: m_kasper $
!
!  search the text: string in the layer-Table, if not found add it to the table
!
!  Input:
!    txtlen    defined length for strings txtgeb
!  In-/ Output: 
!    string    aufzunehmende Layerbezeichnung als String
!    layanz    Anzahl der gelesenen Layer
!    laytxt    Text-Bezeichnung der Layer
!    layrb     type of boundary condition for branches on that layer (and for each nature)
!              0    = Dirichlet             layers: 0, 6-10, DIRIChlet, FIXED
!              1-99 = User Dirichlet        layers: USER number
!              100  = absorbing boundary    layers: ABC 
!              200  = general Neumann       layers: 2, NEUMAnn, FOURIer, ROBINet
!              300  = innner contour        layers: 3, CONTOur, KONTUr, GEOMEtry
!              400  = non-visble contour    layers: 4, INVISable
!              1000 = comment layer         layers: COMMEnt, KOMMEentar, BEMASsung
!              1001 = not recognized layer,
!                     treated as a comment layer
!  Output:
!    layer     Index of the actual layer in the table of layers
!
!  local variables
      integer (I4B) :: i, lyr
!
      call low2hi(string,txtlen)

      do i=1,layanz
        if (string(1:txtlen) .eq. laytxt(i)(1:txtlen)) then
          layer=i
          return
        end if
      end do
!  neuen Layer einlesen
      layanz=layanz+1
      layer=layanz
!  reallocate arrays if neccessary
      if (layanz.gt.size(layrb,1))                                      &
     &    layrb=>reallocate(layrb,2*size(layrb,1)+1,size(layrb,2))
      if (layanz.gt.size(laytxt))                                       &
     &    laytxt=>reallocate(laytxt,2*size(laytxt)+1,txtlen)
!  Name des layers
!  layer ist eine Text - Bezeichnung Standardlayer sind 0,1,2..
      select case (string(1:5))
      case ('DIRIC' )
        layrb(layer,:)=0
      case ('FIXED' )
        layrb(layer,:)=0
      case ('USER' )
        call string2number(string(6:),lyr)
        layrb(layer,:)=lyr
      case ('ABC' )
        layrb(layer,:)=100
      case ('NEUMA' )
        layrb(layer,:)=200
      case ('FOURI' )
        layrb(layer,:)=200
      case ('MIXED' )
        layrb(layer,:)=200
      case ('ROBIN' )
        layrb(layer,:)=200
      case ('KONTU' )
        layrb(layer,:)=300
      case ('CONTO' )
        layrb(layer,:)=300
      case ('GEOME' )
        layrb(layer,:)=300
      case ('INVIS' )
        layrb(layer,:)=400
      case ('KOMME' )
        layrb(layer,:)=1000
      case ('COMME' )
        layrb(layer,:)=1000
      case ('BEMAS' )
        layrb(layer,:)=1000
      case default
!  versuche den Layer-text in eine Zahl umzuwandeln
        read (string,'(i3)',err=477) lyr
        if (lyr.ge.0 .and. lyr.le.10) then
!  layer 6 bis 10 ebenfalls Dirichlet
          if (lyr.ge.6) lyr=0
          layrb(layer,:)=100*lyr
          goto 499
        end if
477     print*
!        print*, '*WARNING unknown layer designation: ',                &
!     &    string(1:len_trim(string))
!        print*, ' ignoring the layer'
!  unknown layer designation set it to dirichlet for the moment
        layrb(layer,:)=0
      end select
499   laytxt(layer)=string(1:txtlen)
      return
      end subroutine ftlay
