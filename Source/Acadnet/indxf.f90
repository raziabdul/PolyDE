      subroutine indxf(anzknt,anzzwg,xpoint,ypoint,anzgeb,xgeb,ygeb,    &
     &  txtgeb,zki,zpz,lzrb,length,nohead,xmin,xmax,ymin,ymax,          &
     &  lytabl,accur,scalfk,lplayt,layanz,laytxt,layrb,snapof,txtlen)
      use feminterface, only: getfile, header, tables, blocks, classes
      use feminterface, only: objects, entity, unsect, reallocate, sectio
      use femtypes
      implicit none
      integer (I4B) :: anzknt, anzzwg, anzgeb, layanz, snapof, txtlen
      integer (I4B), pointer :: zki(:,:), layrb(:,:), lplayt(:), zpz(:)
      integer (I4B), pointer :: lzrb(:)
      real (DP) :: xmin, xmax, ymin, ymax, accur, fakt1, scalfk
      real (DP), pointer :: xgeb(:), ygeb(:), xpoint(:), ypoint(:), length(:)
      logical :: nohead, lytabl
      character (len=txtlen), pointer :: txtgeb(:), laytxt(:)
      intent (in) :: txtlen
!      intent (out) :: anzknt, anzzwg, anzgeb, zki  ! org
!      intent (out) :: zpz, lzrb, nohead, xmin  ! org
!      intent (out) :: xmax, ymin, ymax, lytabl, accur, scalfk, lplayt  ! org
      intent (inout) :: xmax, ymin, ymax, lytabl, accur, scalfk, lplayt  ! test

      intent (out) :: layanz, snapof  ! ORG layrb, laytxt
!      intent (out) :: xpoint, ypoint, length  ! ORG
      intent (inout) :: xpoint, ypoint, length  ! test
      intent (inout) :: anzknt, anzzwg, anzgeb, zki ! test
      intent (inout) :: txtgeb, xgeb, ygeb
      intent (inout) :: zpz, lzrb, nohead, xmin ! test
      !
!    $Revision: 1.10 $
!    $Date: 2008/12/22 15:48:48 $
!    $Author: m_kasper $
!
!  read a DXF-files
!
!  Input:
!              read from the DXF-file
!  Output:
!    anzknt    number of key-points
!    anzzwg    number of branches
!    xpoint    x-coordinates of key-points
!    ypoint    x-coordinates of key-points
!    anzgeb    number of texts read from the DXF-file (number of regions)
!    xgeb      alignment point of text
!    ygeb
!    txtgeb    the string of TEXT and MTEXT entities
!    zki       key-points (geometry nodes) which determine the branches
!               zki(1,i)  start key-point of the branch
!               zki(2,i)  end key-point of the branch
!               zki(3,i)  third key-point of the branch, with
!                         /  = 0  :  straight line 
!               zki(3,i)  -  > 0  :  node located on an arc              
!                         \  < 0  :  midpoint of an arc
!    zpz       number of nodes on the branch (including start- and endpoint)
!    lzrb      for each branch the layer
!    length    length of branches
!    noheaad   =.true. if the DXF-file has no header section
!    xmin
!    xmax      extends of the drawing as determined from the
!    ymin      header section of the DXF file (if present)
!    ymax
!    lytabl    =.true. if the layer table was read
!    accur     accuracy limit
!    scalfk    scaling factor
!    lplayt    layer of the text strings read from the DXF-file
!    layanz    number of layers beeing recognized
!    laytxt    layer name, i.e. designator of layers read from the DXF-file
!    layrb     type of boundary condition for branches on that layer (and for each nature)
!              0    = Dirichlet
!              1-99 = User Dirichlet
!              200  = general Neumann
!              300  = innner contour
!              400  = non-visble contour
!              1000 = comment layer
!              1001 = not recognized layer, treated as a comment layer
!    snapof    number of drawing entities which have been eliminated since there
!              extent is smaller than the drawing accuracy accur or since lines
!              occure repeated with same endpoints
!    txtlen    defined length for strings txtgeb
!
!  local variables
      integer (I4B) :: sction, ios
      real (DP), parameter :: stdacc=1.e-7_DP
      real (DP) :: unitfac
      character (len=4) :: antwrt
!
!  Set the standard values for the accuracy
      accur=stdacc
      nohead=.true.
      lytabl=.false.
      layanz=0
      xmin=0._DP
      xmax=1._DP
      ymin=0._DP
      ymax=1._DP
!
100   call sectio(sction)
!
      select case (sction)
      case(1)
!  Section HEADER
        call header(accur,xmin,xmax,ymin,ymax,unitfac)
        nohead=.false.
      case(2)
!  Section TABLES
        call tables(lytabl,layanz,laytxt,layrb,txtlen)
      case(3)
!  Section BLOCKS
        call blocks
      case(4)
!
!  read ENTITIES i.e. the drawing elements like: lines, circles, texts..
!
!  accuracy in millimeters
        accur=accur * unitfac / 1.e-3_DP
        scalfk=unitfac
        print*
        print*,' do you want to scale the drawing (Y/N)?   <N>'
        print*
        read(*,'(a)') antwrt
        if (antwrt(1:1).eq.'y' .or. antwrt(1:1).eq.'Y' ) then
          print*
          print*,' enter the scaling factor '
          print*
          read(*,*) scalfk
          if ( .not. nohead ) then
            xmax = xmax * scalfk
            xmin = xmin * scalfk
            ymax = ymax * scalfk
            ymin = ymin * scalfk
          end if
          accur = accur * scalfk
          scalfk=scalfk*unitfac
        end if
!
        print*
        print*,' drawing precision is set to be: ',accur,' mm  '
        print*,' do you want to change this value (Y/N)?   <N>'
        print*
        read(*,'(a)') antwrt
        if (antwrt(1:1).eq.'y' .or. antwrt(1:1).eq.'Y' ) then
          print*
          print*,' enter the new precision'
          print*
          read(*,*) accur
        end if
!
!  standard value 4 nodes on a 90 degre arc
        fakt1=4._DP
        print*
        print*,' the number of subdivisions for arcs is set to be ',    &
     &    fakt1,' nodes at 90 degrees'
        print*,' do you want to change this values (Y/N)?   <N>'
        print*
        read(*,'(a)') antwrt
        if (antwrt(1:1).eq.'y' .or. antwrt(1:1).eq.'Y' ) then
          print*
          print*,' enter the new value'
          print*
          read(*,*) fakt1
        end if
!
        call entity(accur,fakt1,scalfk,xpoint,ypoint,zki,lzrb,zpz,      &
     &    xgeb,ygeb,txtgeb,length,anzzwg,anzknt,anzgeb,laytxt,layrb,    &
     &    lplayt,layanz,snapof,txtlen)
      case(5)
!  Section CLASSES
        call classes
      case(6)
!  Section OBJECTS
        call objects
      case(7)
!  Section THUMBNAILIMAGE  ignore
        call unsect
      case(99)
        call unsect
      case(10)
        goto 300
      case default
        print*,'*ERROR: Syntax error in dxf file: unknown section'
        call unsect
      end select
      goto 100
!
!  finish
300   close (38,iostat=ios)
      return
      end subroutine indxf


