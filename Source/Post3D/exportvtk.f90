      subroutine exportvtk(fieldtype,phi)
      use feminterface, only: getsetting, getpostsetting
      use feminterface3D, only: field3D, findelement3D, xyz2lam
      use femtypes
      use globalvariables3D, only : nnat, nod, numn, numv, vn, x
      implicit none
      real (DP) :: phi
      character (len=10) :: fieldtype
      intent (in) :: fieldtype, phi
!
!
!  Subroutine to export values for given fieldtype and phase as an
!  unstructured grid output according to the vtk "classic" format.
!  Experimental
!
!  Input:
!     fieldtype          field quantity to compute (e.g. POTENTIAL, EY ...)
!     phi                phase / time instant to compute

!  Internal variables:
      integer (I4B) :: i, ierror, unitid, inat
!!$      real (DP) :: dist1, dist2, vec1(3), vec2(3), unitvec1(3), unitvec2(3)
!!$      real (DP) :: delta(2), lambda(4), point(3)
!!$      real (DP), allocatable :: val1(:), val2(:), zval(:,:)
!!$      character (len=8)  :: dat
      character (len=10) :: tmp
!!$      character (len=23) :: fmtreal
      character (len=200):: path, vtk_file
!!$      complex (DPC) :: u(3), curlu(3)
      logical :: opened
!
!      
      external :: grglun

!  write data to vtk file
      call getsetting('PROJECTPATH',path)
      call getpostsetting('VTK_FILE',vtk_file)
      vtk_file=adjustl(vtk_file)
      inquire(FILE=path(1:len_trim(path))//vtk_file(1:len_trim(vtk_file))//'.vtk', &
     &        OPENED=opened,NUMBER=unitid,IOSTAT=ierror)
      if (.not.opened) then
!  open a file in the specified project path
        call grglun(unitid)
        open (unitid,File=path(1:len_trim(path))//vtk_file(1:len_trim(vtk_file))//'.vtk', &
     &        FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
        if (ierror .ne. 0) then
          print*, 'Cannot open file: ', path(1:len_trim(path))//'Meshquality.vtk'
        end if

!
!  write vtk header
!  make sure # starts in column 1
      write (unitid,'(a)') '# vtk DataFile Version 3.0'
      write (unitid,'(a)') '3D scalar data'
      write (unitid,'(a)') 'ASCII'
      write (unitid,'(a)') 'DATASET UNSTRUCTURED_GRID'
      write (unitid,'(a,i20,a)') 'POINTS', numn, '  double'
!  write coordinates of nodes
      do i=1,numn
         write (unitid,'(3E17.8E3)') nod(1,i), nod(2,i), nod(3,i)
         !write (unitid,'(3E17.8E3)') nod(1,i)+2.*real(x(i)), nod(2,i)+2.*real(x(i+numn)), nod(3,i)+2.*real(x(i+2*numn))
      end do
!  header for elements preceded by a blank line
      write (unitid,'(a,2i20)') 'CELLS', numv, numv*5
!  write volume element connectivity
!  VTK-Indices are 0-offset. Thus the first point is point id 0
      do i=1,numv
         write (unitid,*) '4', vn(:,i)-1
      end do
!  header for element (cell) types preceded by a blank line
!  for tet the vtk code is 10
      write (unitid,'(a,i20)') 'CELL_TYPES', numv
!  write cell type for each volume element
!  10   VTK_TETRA            nodes 0..3 for vertices
!  24   VTK_QUADRATIC_TETRA  4 between vertex 0 and 1
!                            5 between        1     2
!                            6 between        0     2
!                            7 between        0     3
!                            8 between        1     3
!                            9 between        2     3
      do i=1,numv
         write (unitid,*) '10'
      end do
!  header for scalar data for each node preceded by a blank line
      write (unitid,'(a,i20)') 'POINT_DATA', numn

      end if

      do inat=1,nnat

          write (tmp ,'(i2)') inat
      write (unitid,*) 'SCALARS Potential_Nature'//trim(adjustl(tmp)),' double'
      write (unitid,*) 'LOOKUP_TABLE default'

      do i=1+(numn*(inat-1)),numn*(inat)
!  write field value at each node
         write (unitid,'(g13.5)') real(x(i))
      end do
      end do
      
      
      write (unitid,*) 'VECTORS Displacement',' double'
      write (unitid,*) 'LOOKUP_TABLE default'
      do i=1, numn
!  write field value at each node
         write (unitid,'(3g13.5)') real(x(i)), real(x(i+numn)), real(x(i+2*numn))
      end do

!      close (unitid)
!      call grglun(unitid)
!
      end subroutine
