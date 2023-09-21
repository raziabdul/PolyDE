      subroutine polyordervtk()
      use feminterface, only: getsetting, getpostsetting
      use feminterface3D, only:
      use femtypes
      use globalvariables3D, only : nnat, nod, numn, numv, vn, vp
      implicit none
!
!  Subroutine to export Polyoder  into a vtk file (Simple Legacy Formats )
!
!  Local variables:
      integer (I4B) :: i, ierror, unitid, inat
      logical :: opened
      character (len=3):: tmp
      character (len=200):: path, vtk_file
!      
      external :: grglun

!  write data to file data<n>.vtk
      call getsetting('PROJECTPATH',path)
      call getpostsetting('VTK_FILE',vtk_file)
      vtk_file=adjustl(vtk_file)
      inquire(FILE=path(1:len_trim(path))//vtk_file(1:len_trim(vtk_file))//'.vtk', &
     &        OPENED=opened,NUMBER=unitid)
      if (.not.opened) then
!  open a file in the specified project path
        call grglun(unitid)
        open (unitid,File=path(1:len_trim(path))//vtk_file(1:len_trim(vtk_file))//'.vtk', &
     &        FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
        if (ierror .ne. 0) then
          print*, 'Cannot open file: ', path(1:len_trim(path))//'Polyorder.vtk'
        end if
!
!  write vtk header
!  make sure # starts in column 1
        write (unitid,'(a)') '# vtk DataFile Version 3.0'
        write (unitid,'(a)') '3D scalar data'
        write (unitid,'(a)') 'ASCII'
        write (unitid,'(a)') 'DATASET UNSTRUCTURED_GRID'
!
!  POINTS
        write (unitid,'(a,i20,a)') 'POINTS', numn, '  float'
!  write coordinates of nodes
        do i=1,numn
           write (unitid,'(3f15.8)') nod(1,i), nod(2,i), nod(3,i)
        end do
!
!  CELLS
!  header for elements preceded by a blank line
        write (unitid,'(a,2i20)') 'CELLS', numv, numv*5
!  write volume element connectivity
        do i=1,numv
!  element node indices in VTK are 0-based
!  VTK-Indices are 0-offset. Thus the first point is point id 0
           write (unitid,*) '4', vn(:,i)-1
        end do
!
!  CELL TYPES
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
!  CELL_DATA
      write (unitid,'(a,i20)') 'CELL_DATA', numv
      end if
!
!  header for scalar data for each node preceded by a blank line
      do inat=1, nnat
        write (tmp ,'(i3)') inat
        write (unitid,'(a)') 'SCALARS Polyorder_Nature'//adjustl(tmp)//' int'
        write (unitid,*) 'LOOKUP_TABLE default'

!  write field value at each element
        do i=1, numv
          write (unitid,*) vp(i,inat)
        end do

      end do

!      close (unitid)
!
      return
      end subroutine
