      program post3D
      use feminterface, only: getsetting, getpostsetting, readpostsetting, timestamp, low2hi
      use feminterface, only: zeit
      use feminterface3D, only: field3D, fieldsc3D, initialize3D, solin, findelement3D
      use feminterface3D, only: vol2aux, xyz2lam, exportvtk, readunv, pointvalue3D, polyordervtk, meshqualityvtk
      use feminterface3D, only: vtk_scalarpointplot, vtk_vectorpointplot, vtk_tensorpointplot
      use feminterface3D, only: vtk_scalarcellplot, vtk_vectorcellplot, exportdata3D
      use feminterface3D, only: exportvtkres, vtk_scalargridplot
      use femtypes
      use globalvariables3D, only : c0, eltype, nnat, nod, numv, omega, pi, vn, vv
      implicit none
!
!------------------------------------------------------------------------------
!
!    $Revision: 1.16 $
!    $Date: 2014/08/27 15:20:53 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
!
!  3D postprocessor
!
      integer (I4B) :: i, div, elem, grid1, grid2, grid3, unitid, inat, dotpos
      real (DP) :: epsgl
      real (DP) :: dist, delta, lambda(4), phi
      real (DP) :: origin(3), p1(3), p2(3), p3(3)
      real (DP) :: startp(3), endp(3), vec(3), unitvec(3)
      real (DP), allocatable :: distval(:), points(:,:), zval(:,:,:)
      complex (DPC), allocatable :: z(:,:)
      complex (DPC), allocatable :: usc(:), alphau(:,:)
      complex (DPC), allocatable :: u(:,:), curlu(:,:)
      character (len=10) :: fieldtype
      character (len=30)  :: callback
      character (len=200) :: path
      character (len=20) :: nature
      character(len=50) :: meshfile
      logical :: ok, gilt, eleinfo, found
      logical :: typ(5)
! For Command Line Argument
      character (len=10) :: argv
      integer numtime
!      integer ii, iargc, nn
!
      external :: grglun
!------------------------------------------------------------------------------
!      nn=iargc()
!      do ii=1,nn
      call getarg(1,argv)
      numtime=0
!      end do
      call initialize3D()
!
      print "(a,g10.3/,a,g10.3)","frequency        : ",omega/(2._DP*pi)&
                              ,"angular frequency: ",omega
      if (omega .gt. 0) then
        print "(a,g10.3)","vacuum wavelength: ",(c0 / omega)*2._DP*pi
      end if
!
!  read mesh from file basis.mesh from PROJECTPATH
      call getsetting('MESHFILE',meshfile)
      dotpos = index(meshfile,'.',BACK=.true.)
      if (dotpos .eq. 0) then 
        print*,'Error: Extension of meshfile: ',meshfile(1:len_trim(meshfile)),' is missing'
      end if
      call low2hi(meshfile,len_trim(meshfile))
      if (meshfile(dotpos+1:dotpos+4) .eq. 'UNV') then
!  read mesh file in UNV format
        call readunv('basis.unv',ok)
!      else if(  .eq. 'VTK) then
!        call vtkin(ok)                  ! to be implemented
      else
!  assume netgen mesh file
        print *, 'Only UNV mesh format accepted'
      end if

      call zeit(' Reading the Mesh')

!
!  read in solution from file
      call solin(ok, gilt, eleinfo, epsgl)
      call zeit(' Reading the Solution')
!
!  stop programm if no valid solution
      if ( (.not. ok) .or. (.not. gilt) .or. (.not. eleinfo) ) stop
!
      call vol2aux
!  loop for drawing commands
      do
        call readpostsetting(callback)
!
        if (numtime .eq. 0) then
          if (argv .eq. 'line') callback='LINEGRAPH'
          if (argv .eq. 'vtk') callback='EXPORT_VTK_UNSTR_ASCII'
          numtime=1
        end if
        select case(callback)
!------------------------------------------------------------------------------
!  leave the program
        case('END')
          print *,"END-command used. Plotting stopped."
          exit
!------------------------------------------------------------------------------
!  End Of File is reached
        case('EOF')
          print *,"End Of File reached"
          exit
!------------------------------------------------------------------------------
!  an error occured during reading the settings in readpostsettings.f90
        case('ERROR')
          exit
!------------------------------------------------------------------------------
!  write field values of given fieldtype to file data.txt
        case ('EXPORTDATA')
          if (.not.gilt) then
            print *,"No valid solution for calculating value"
            cycle
          end if
!
!  get fieldtype / -quantity, number of gridpoints in right- and up-direction and
!  phase from POSTsettings.txt
          call getpostsetting('FIELDTYPE',fieldtype)
          call getpostsetting('PHI',phi)
!  get points to define plane
          call getpostsetting('ORIGINX',origin(1))
          call getpostsetting('ORIGINY',origin(2))
          call getpostsetting('ORIGINZ',origin(3))
          call getpostsetting('P1X',p1(1))
          call getpostsetting('P1Y',p1(2))
          call getpostsetting('P1Z',p1(3))
          call getpostsetting('P2X',p2(1))
          call getpostsetting('P2Y',p2(2))
          call getpostsetting('P2Z',p2(3))
!  get number of grid points
          call getpostsetting('GRID1',grid1)
          call getpostsetting('GRID2',grid2)
!  export grid values to file data_<fieldtype>.txt
          call exportdata3D(fieldtype,phi,origin,p1,p2,grid1,grid2)
!------------------------------------------------------------------------------
!  write field values of given fieldtype in 3D using the format of
!  plot3D.
!  http://www.lerc.nasa.gov/WWW/wind/valid/plot3d.html
!  http:http://people.sc.fsu.edu/~burkardt/data/plot3d/plot3d.html

        case ('EXPORT_VTK_UNSTR_ASCII')
          if (.not.gilt) then
            print *,"No valid solution for calculating value"
            cycle
          end if
!
!  get fieldtype / -quantity, number of gridpoints in right- and up-direction and
!  phase from POSTsettings.txt
          call getpostsetting('FIELDTYPE',fieldtype)
          call getpostsetting('PHI',phi)

!  export  field values to files: 
          call exportvtk(fieldtype,phi)
          call exportvtkres(fieldtype,phi)
!------------------------------------------------------------------------------
!!$        case ('EXPORT3D')
!!$          if (.not.gilt) then
!!$            print *,"No valid solution for calculating value"
!!$            cycle
!!$          end if
!!$!
!!$!  get fieldtype / -quantity, number of gridpoints in right- and up-direction and
!!$!  phase from POSTsettings.txt
!!$          call getpostsetting('FIELDTYPE',fieldtype)
!!$          call getpostsetting('PHI',phi)
!!$!  get points to define plane
!!$          call getpostsetting('ORIGINX',origin(1))
!!$          call getpostsetting('ORIGINY',origin(2))
!!$          call getpostsetting('ORIGINZ',origin(3))
!!$          call getpostsetting('P1X',p1(1))
!!$          call getpostsetting('P1Y',p1(2))
!!$          call getpostsetting('P1Z',p1(3))
!!$          call getpostsetting('P2X',p2(1))
!!$          call getpostsetting('P2Y',p2(2))
!!$          call getpostsetting('P2Z',p2(3))
!!$!  this point gives the "height"
!!$          call getpostsetting('P3X',p3(1))
!!$          call getpostsetting('P3Y',p3(2))
!!$          call getpostsetting('P3Z',p3(3))
!!$!  get number of grid points
!!$          call getpostsetting('GRID1',grid1)
!!$          call getpostsetting('GRID2',grid2)
!!$          call getpostsetting('GRID3',grid3)
!!$
!!$!  export grid values to files: 
!!$          call exportplot3D(fieldtype,phi,origin,p1,p2,p3,grid1,grid2,grid3)
        case('LINEGRAPH')
           print *, 'IN LINEGGGGG'
!  Store solution data for a line in 3D space given by start and end points,
!  division and phase angle in file linegraph3D.txt.
!  get settings for output
          call getpostsetting('STARTX',startp(1))
          call getpostsetting('STARTY',startp(2))
          call getpostsetting('STARTZ',startp(3))
          call getpostsetting('ENDX',endp(1))
          call getpostsetting('ENDY',endp(2))
          call getpostsetting('ENDZ',endp(3))
          call getpostsetting('DIVISION',div)
          call getpostsetting('PHI',phi)
!  compute direction vector from start to end
          vec = endp - startp
!  compute length of vector and unitvector (length 1)
          dist = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
          unitvec = vec/dist
!  compute delta depending on # of divisions on length
          delta = dist/div
!  allocate and initialize arrays
          if ( eltype .eq. 'SCALAR') then
            allocate( z(15,nnat), usc(nnat), alphau(nnat,nnat) )
            allocate( distval(div+1), points(3,div+1), zval(1,div+1,nnat) )
            distval = 0._DP
            points = 0._DP
            zval = 0._DP
!  loop over all points
            do i = 1, div+1
!  compute x,y coordinates for new point depending on # of division
              distval(i) = (i-1)*delta
              points(:,i)= startp + (i-1)*delta*unitvec
              call findelement3D(points(:,i),nod,vn,vv,numv,elem,found)
              if (found) then
                call xyz2lam(lambda, elem, points(:,i), nod, vn)
                typ=(/ .true., .false., .false., .false., .false. /)
                call fieldsc3D(elem,lambda,typ,z,u=usc,alphau=alphau)
                do inat=1,nnat
                  zval(1,i,inat) = real(usc(inat))
!                  zval(:,i) = real(curlu*exp(cmplx(0._DP,phi*pi/180._DP,DPC)),DP)
                end do
              end if
            end do 
          else 
            allocate(u(3,nnat),curlu(3,nnat))
            allocate(distval(div+1), points(3,div+1), zval(3,div+1,nnat))
            distval = 0._DP
            points = 0._DP
            zval = 0._DP
!  loop over all points
            do i = 1, div+1
!  compute x,y coordinates for new point depending on # of division
              distval(i) = (i-1)*delta
              points(:,i)= startp + (i-1)*delta*unitvec
              call findelement3D(points(:,i),nod,vn,vv,numv,elem,found)
              if (found) then
                call xyz2lam(lambda, elem, points(:,i), nod, vn)
                call field3D(elem,lambda,u,curlu)
                do inat=1,nnat
                  zval(:,i,inat) = real(u(:,inat)*exp(cmplx(0._DP,phi*pi/180._DP,DPC)),DP)
!                  zval(:,i) = real(curlu*exp(cmplx(0._DP,phi*pi/180._DP,DPC)),DP)
                end do
              end if
            end do         
          end if 
!  writing to file
          call getsetting('PROJECTPATH',path)
          do inat = 1, nnat
            write(nature,*) inat
            nature=adjustl(nature)
            call grglun(unitid)
!  open a file in the specified project path
            open (unitid,file=path(1:len_trim(path))//'linegraph3D'//trim(nature)//'.txt',  &
     &            form='formatted',position='REWIND',action='WRITE')
!  write the data to the file
            do i = 1, size(distval)
              write (unitid,100) distval(i), zval(:,i,inat)
100           format(4G22.14)
            end do
            close (unitid)
          end do
!         call grglun(unitid)
!         open (unitid,file=path(1:len_trim(path))//'linegraph3D2.txt',    &
!     &   form='formatted',position='REWIND',action='WRITE')
!         do i = 1, size(distval)
!           write (unitid,100) distval(i), zval(:,i,2)
!         end do
!         close (unitid)
          deallocate(distval,points,zval)
        case ('MESHQUALITY')
          call meshqualityvtk()
        case ('POINTVALUE')
          call pointvalue3D
!------------------------------------------------------------------------------
!  draw the polynomial order of the elements
        case ('POLYORDER')
          if (.not.gilt) then
            print *,"**** Error No valid solution for plotting"
            cycle
          end if
          call polyordervtk()
        case ('SCALARPOINTPLOT')
          call vtk_scalarpointplot()
        case ('VECTORPOINTPLOT')
          call vtk_vectorpointplot()
        case ('TENSORPOINTPLOT')
          call vtk_tensorpointplot()
        case ('SCALARCELLPLOT')
          call vtk_scalarcellplot()
        case ('VECTORCELLPLOT')
          call vtk_vectorcellplot()
        case ('SCALARGRIDPLOT')
          call vtk_scalargridplot()
        end select
      end do
!------------------------------------------------------------------------------
!
!  deallocate arrays
      call zeit(' End of Post3D')
      end program post3D
