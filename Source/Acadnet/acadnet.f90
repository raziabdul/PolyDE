      program acadnet
      use feminterface, only: zanfpp, indxf, sortnt, testac, setzwg
      use feminterface, only: extend, outdxf, struc, getsetting
      use feminterface, only: wnetin, rbtxt, gengeb, elmgeb, matgeb
      use feminterface, only: wpathsetting, zeit
      use globalvariables, only: nnat
      use femtypes
      implicit none
!
!    $Revision: 1.30 $
!    $Date: 2014/02/11 16:55:03 $
!    $Author: juryanatzki $
!
      integer (I4B), parameter :: txtlen=1000
      integer (I4B) :: anzzwg, anzknt, anzgeb, i, j, gbz, layanz
      integer (I4B) :: bildnr, snapof, unitid
      integer (I4B), allocatable:: kzrb(:,:), donetx(:)
      integer (I4B), pointer :: bzil(:), bzi(:), bzip(:), zki(:,:)
      integer (I4B), pointer :: lzrb(:), zpz(:), zrb(:,:)
      integer (I4B), pointer :: lplayt(:), layrb(:,:)
      real (DP) :: xmin, xmax, ymin, ymax, accur
      real (DP) :: h, rand, scalfk, xpl, xpu, ypl, ypu
      real (SP) :: red = 1.0, green = 0.0, blue = 0.0
      real (DP), pointer :: xgeb(:), ygeb(:), xpoint(:), ypoint(:)
      real (DP), pointer :: area1(:), zpp(:), length(:)
      logical :: nohead, lytabl, exist
      complex (DPC), allocatable:: alrb(:,:), btrb(:,:)
      character (len=25) :: system, plotdevice
      character (len=200) :: path
      character (len=25), pointer :: matname(:)
      character (len=txtlen), pointer :: txtgeb(:), laytxt(:)
!
      call acadini
!
!  arrays related to key-point
      allocate (xpoint(200), ypoint(200))
!  arrays related to branches
      allocate (lplayt(200), laytxt(200), zki(3,200), zpz(200), lzrb(200), length(200))
      allocate (layrb(200,nnat))
      allocate(bzi(200))
      allocate(bzip(200))
      allocate(bzil(200))

!  arrays related to regions
      allocate (xgeb(100), ygeb(100), txtgeb(100))
!
!  read the dxf-file
!
      call indxf(anzknt,anzzwg,xpoint,ypoint,anzgeb,xgeb,ygeb,          &
     &  txtgeb,zki,zpz,lzrb,length,nohead,xmin,xmax,ymin,ymax,          &
     &  lytabl,accur,scalfk,lplayt,layanz,laytxt,layrb,snapof,txtlen)
!
!  sorting of key-points and branches from left to right and bottom to top
!
      call sortnt(anzknt,anzzwg,xpoint,ypoint,zki,zpz,lzrb)
!
      write(*,*)
      write(*,*) anzknt,' key-points  ',anzzwg,' branches'
!
!  determine number of nodes and separation on branches
!
      allocate (zpp(anzzwg))
      zpp=0._DP
!      call setzwg(length,anzzwg,anzknt,zki,zpz,zpp)
!
!  recompute Extends
      call extend(anzzwg,xpoint,ypoint,zki,xmin,xmax,ymin,ymax)
      write (*,*)
      write (*,371) xmin,ymin,xmax,ymax
371   format (' Extends of the AUTO-CAD drawing:'                     &
     &  ,2x,'(',g10.3,' , ',g10.3,') to (',g10.3,' , ',g10.3,')')
      write (*,*)
!
!  output drawing as a dxf-file
!
      call outdxf(xpoint,ypoint,zki,lzrb,laytxt,anzzwg,                 &
     &  xgeb,ygeb,txtgeb,(xmax-xmin)/100._DP,anzgeb,                    &
     &  lplayt)
!
!  screen output of the drawing
!
!  margin in percent of total size
      rand=.02_DP
      xpl=xmin-(xmax-xmin)*rand/2._DP
      xpu=xmax+(xmax-xmin)*rand/2._DP
      ypl=ymin-(ymax-ymin)*rand/2._DP
      ypu=ymax+(ymax-ymin)*rand/2._DP
!
      call getsetting('WHATSYSTEM',system)
      select case (system)
      case ('WINDOWS')
        plotdevice='/WX'
      case ('LINUX')
        plotdevice='/XWINDOW'
      case default
        plotdevice='/XWINDOW'
      end select
      call zanfpp(plotdevice,xpl,xpu,ypl,ypu,h,bildnr)
!
!  zrb is set to an arbitrary value
      allocate (zrb(anzzwg,nnat))
      zrb=0
      call struc(anzzwg,zki,zrb,xpoint,ypoint,red,green,blue)
!  write sample for POSTsettings.txt
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
      inquire (file=path(1:len_trim(path))//'POSTsettings.txt',exist=exist)
      if (.not. exist) then
        open (unitid,file=path(1:len_trim(path))//'POSTsettings.txt')
      else
        open (unitid,file=path(1:len_trim(path))//'POSTsettings.bak')
      end if
      write(unitid,*)'!  set diplay device and extents'
      write(unitid,*)'DISPLAYDEVICE = ',plotdevice
      write(unitid,*)'XMIN        = ',xpl
      write(unitid,*)'XMAX        = ',xpu
      write(unitid,*)'YMIN        = ',ypl
      write(unitid,*)'YMAX        = ',ypu
      write(unitid,*)'!-------------------------------------------------'
      write(unitid,*)'!  open window with full field plot and structure'
      write(unitid,*)'INFO        = YES'
      write(unitid,*)'OPEN'
      write(unitid,*)'AMIN        = -1.'
      write(unitid,*)'AMAX        = 1.'
      write(unitid,*)'PHI         = 0.'
      write(unitid,*)'DRAWINGTYPE = FULLPLOT'
      write(unitid,*)'NUMCOLORS   = -128'
      write(unitid,*)'NUMLINES    = -25'
      write(unitid,*)'LINECOLOR   = 000000000'
      write(unitid,*)'FIELDTYPE   = POTENTIAL'
      write(unitid,*)'FIELD'
      write(unitid,*)'LINECOLOR   = 200200200'
      write(unitid,*)'STRUCT'
      write(unitid,*)'LINECOLOR   = 000000000'
      write(unitid,*)'SCALE'
      write(unitid,*)'INFO        = NO'
      write(unitid,*)'!-------------------------------------------------'
      write(unitid,*)'!  open window with full field plot and structure'
      write(unitid,*)'OPEN'
      write(unitid,*)'STRUCT'
      write(unitid,*)'SCALE'
      write(unitid,*)'MOUSE'
      write(unitid,*)'LINECOLOR   = 255000000'
      write(unitid,*)'LINEWIDTH   = 3'
      write(unitid,*)'DIVISION    = 1000'
      write(unitid,*)'OUTPUT      = YES'
      write(unitid,*)'LINEGRAPH'
      close(unitid)
!
!  Output the structure without the regions
      allocate (alrb(layanz,nnat), btrb(layanz,nnat))
      alrb=(0._DP,0._DP)
      btrb=(0._DP,0._DP)
!  set boundary condition of key-points to an arbitrary value
      allocate (kzrb(anzknt,nnat))
      kzrb=1
!
!  write the geometry to Netin.acd (for now without the regions, kzrb is set to an arbitrary value)
      bzi(:) = 0
! Dummy allocate to avoid runtime error in oneAPI 
      allocate(matname(200))
      call wnetin(0,anzzwg,anzknt,bzi,bzip,bzil,zki,kzrb,zpz,zpp,       &
     &  xpoint,ypoint,matname,path(1:len_trim(path))//'netin.acd',      &
     &  lzrb,layrb,alrb,btrb,laytxt,layanz,nnat)
      deallocate(matname)
     !
!  set the boundary conditions of branches according to the text found
      allocate (donetx(anzgeb))
      call rbtxt(txtgeb,anzgeb,alrb,btrb,donetx,lplayt,layrb,nnat)
!
!  test for consistency of the drawing and boundary conditions
      call testac(anzzwg,anzknt,xpoint,ypoint,zki,lzrb,layrb)
!
!  notify about eliminated branches before a potential program abort
      if (snapof .gt. 0) then
        write (*,534) snapof
534     format ('    ',i4,' entities of the drawing had been eliminated')
      end if
!
!  call automatic generation of regions
!  bzi to be allocated in gengeb
      deallocate(bzi) 
      call gengeb(gbz,anzzwg,anzknt,bzi,bzip,bzil,zki,xpoint,ypoint,area1)
      write (*,*)
!
!  assign texts (material names) to the regions
      allocate (matname(gbz))
      call matgeb(gbz,anzgeb,bzi,bzip,bzil,zki,area1,xpoint,ypoint,     &
     &  xgeb,ygeb,txtgeb,matname,donetx)
      write (*,*)
!
!  eliminate holes
      call elmgeb(gbz,bzi,bzip,bzil,lzrb,layrb,matname)
!
!  determine boundary condition of key-points:
!  use BC's of the adjacent branches and the priority rule
      do j=1,nnat
        do i=1,anzzwg
          if ( layrb(lzrb(i),j) .lt. kzrb(zki(1,i),j) ) then
            kzrb(zki(1,i),j)=i
          end if
          if ( layrb(lzrb(i),j) .lt. kzrb(zki(2,i),j) ) then
            kzrb(zki(2,i),j)=i
          end if
        end do
      end do
!
!  write the geometry to Netin.acd (now including the regions and correct kzrb)
      call wnetin(gbz,anzzwg,anzknt,bzi,bzip,bzil,zki,kzrb,zpz,zpp,     &
     &  xpoint,ypoint,matname,path(1:len_trim(path))//'netin.acd',      &
     &  lzrb,layrb,alrb,btrb,laytxt,layanz,nnat)
!
      call zeit(' End')
      call pgend
!
      stop
      end program acadnet
