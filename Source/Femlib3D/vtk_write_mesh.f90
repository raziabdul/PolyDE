      subroutine vtk_write_mesh(unitid, ellist)
      use femtypes
      use feminterface,       only: getpostsetting
      use globalvariables3D,  only: numn, nod, x, numv, vn
      use postsettings
      implicit none
      integer (I4B) :: unitid
      integer(I4B), optional :: ellist(:)
      intent (in) :: ellist
!
!  Subroutine to export values for given fieldtype and phase as an
!  unstructured grid output according to the vtk "classic" format.
!
!  Input:

!  Internal variables:
      integer (I4B)           :: i
      real (DP)               :: deffac
      character (len=20)      :: shape
      logical                 :: ok
!  make sure # starts in column 1
      write (unitid,'(a)') '# vtk DataFile Version 3.0'
      write (unitid,'(a)') '3D scalar data'
      write (unitid,'(a)') 'ASCII'
      write (unitid,'(a)') 'DATASET UNSTRUCTURED_GRID'
!
!  POINTS
      write (unitid,'(a,i20,a)') 'POINTS', numn, '  double'
!  write coordinates of nodes
      if (fileopened) call getpostsetting('SHAPE',shape)
      if (shape .eq. 'DEFORMED') then
        do i=1,numn
!  TO DO the factor here is a fixed value
          deffac = 2.0_DP
          write (unitid,'(3E17.8E3)') nod(1,i)+deffac*real(x(i)),         &
     &                                nod(2,i)+deffac*real(x(i+numn)),    &
     &                                nod(3,i)+deffac*real(x(i+2*numn))
        end do
      else
        do i=1,numn
          write (unitid,'(3E17.8E3)') nod(1,i), nod(2,i), nod(3,i)
      end do
      end if
!
!  CELLS 'Tetrahedrons'
!  header for elements preceded by a blank line
      if (present(ellist)) then
        write (unitid,'(a,2i20)') 'CELLS', size(ellist), size(ellist)*5
      else
        write (unitid,'(a,2i20)') 'CELLS', numv, numv*5
      end if
!  write volume element connectivity
!  VTK-Indices are 0-offset. Thus the first point is point id 0
      if (present(ellist)) then
        do i=1,size(ellist)
           write (unitid,'(a, 4i)') '4', vn(:,ellist(i))-1
        end do
      else
        do i=1,numv
          write (unitid,'(a, 4i)') '4', vn(:,i)-1
        end do
      end if
!
!  CELL TYPES Tetrahedrons
!  header for element (cell) types preceded by a blank line
!  for tet the vtk code is 10
      if (present(ellist)) then
        write (unitid,'(a,i20)') 'CELL_TYPES', size(ellist)
      else
        write (unitid,'(a,i20)') 'CELL_TYPES', numv
      end if
!  write cell type for each volume element
!  10   VTK_TETRA            nodes 0..3 for vertices
!  24   VTK_QUADRATIC_TETRA  4 between vertex 0 and 1
!                            5 between        1     2
!                            6 between        0     2
!                            7 between        0     3
!                            8 between        1     3
!                            9 between        2     3
      if (present(ellist)) then
        do i=1,size(ellist)
          write (unitid,'(a)') '10'
        end do
      else
        do i=1,numv
          write (unitid,'(a)') '10'
        end do
      end if
!
      return
      end
