      subroutine vtk_write_cells(unitid, cell_type, ellist)
      use femtypes
      use feminterface,       only: getpostsetting, destroyarrptr
      use feminterface3d,     only: reversemap
      use globalvariables3D,  only: numn, numf, nod, x, numv, vn, vf
      implicit none
      integer (I4B) :: unitid, cell_type
      integer(I4B), optional :: ellist(:)
      intent (in) :: ellist, cell_type
!
!-------------------------------------------------------------------------------
!
! This Routine: 'vtk_write_cells'
! Subroutine to write Cells into a Vtk File.
!
!-------------------------------------------------------------------------------
! Input:
!        unitid          unitid of the file to be written in
!     cell_type          The Type if Cell which is written into the file.
!                        2 Types are implemented so far:
!                        - (5)  Triangles
!                        -(10)  Tetrahedrons (Linear, 1-st Order)
!        ellist          Optional input variable. Marks elements which are to be extracted.
!                        If ellist is absent, all elements are written
! Output:
!        valid vtk export, file is not closed.
!
!-------------------------------------------------------------------------------
! local variables:
!
      integer (I4B)           :: i, f, n, size_flist, fn(3), nmin, nmax
      real (DP)               :: deffac
      character (len=20)      :: shape
      type(ARRPTRI),  pointer :: fv(:)=>null()
      logical                 :: ok
!__
! START:
!  write vtk header
!  make sure # starts in column 1
      write (unitid,'(a)') '# vtk DataFile Version 3.0'
      write (unitid,'(a)') '3D scalar data'
      write (unitid,'(a)') 'ASCII'
      write (unitid,'(a)') 'DATASET UNSTRUCTURED_GRID'
!
!  POINTS
      write (unitid,'(a,i20,a)') 'POINTS', numn, '  double'
!  write coordinates of nodes
      call getpostsetting('SHAPE',shape)
      if (shape .eq. 'DEFORMED') then
        do i=1,numn
!  TO DO the factor here is a fixed value
          deffac = 2.0_DP
          write (unitid,'(3E17.8E3)') nod(1,i)+deffac*real(x(i)),         &
     &                              nod(2,i)+deffac*real(x(i+numn)),    &
     &                              nod(3,i)+deffac*real(x(i+2*numn))
        end do
      else
        do i=1,numn
          write (unitid,'(3E17.8E3)') nod(1,i), nod(2,i), nod(3,i)
      end do
      end if
!__
! Write Cells according to the necessary cell_type   
      select case (cell_type)
!-
        case (5) ! CELLS 'Triangles'
!- 1) Get Auxilarity Data: fv
          call reversemap(.true.,vf,fv,nmin,nmax)
          if (present(ellist)) then
          size_flist = 0
            do f = 1, numf
              if (any(ellist(:).eq.fv(f)%d(1)).or.any(ellist(:).eq.fv(f)%d(2))) then
                size_flist = size_flist + 1
              end if
            end do
          end if
!- 2) Write 'Triangles'
          ! header for elements preceded by a blank line
          if (present(ellist)) then
            write (unitid,'(a,2i20)') 'CELLS', size_flist, size_flist*4
          else
            write (unitid,'(a,2i20)') 'CELLS', numf, numf*4
          end if
         ! write volume element connectivity
         ! VTK-Indices are 0-offset. Thus the first point is point id 0
          if (present(ellist)) then
!- If 'ellist' present
          ! Only proceed, if that face belongs to one of the allowed elements 
            if (any(ellist(:).eq.fv(f)%d(1)).or.any(ellist(:).eq.fv(f)%d(2))) then
              do f=1,numf
              ! Determine the local number of the face f within the volume/element/tetrahedron fv(f)%d(1)
                do n = 1, 4
                  if (vf(n,fv(f)%d(1)).eq.f) exit
                end do
              ! That number. in the same time, marks the node, which does not belong to the face
                fn(1:3) = 0
                do i = 1, 4
                  if (i.ne.n) then
                  ! Node for the face found, look for where to store:
                    if (fn(1).eq.0) then
                      fn(1) = vn(i,fv(f)%d(1))
                    else if (fn(2).eq.0)then
                      fn(2) = vn(i,fv(f)%d(1))
                    else
                      fn(3) = vn(i,fv(f)%d(1))
                    end if
                  end if
                end do
                write (unitid,'(a, 3i)') '3', fn(:)-1
              end do
            end if
          else
!- If All elements are to be considered
            do f=1,numf
            ! Determine the local number of the face f within the volume/element/tetrahedron fv(f)%d(1)
              do n = 1, 4
                if (vf(n,fv(f)%d(1)).eq.f) exit
              end do
            ! That number. in the same time, marks the node, which does not belong to the face
              fn(1:3) = 0
              do i = 1, 4
                if (i.ne.n) then
                ! Node for the face found, look for where to store:
                  if (fn(1).eq.0) then
                    fn(1) = vn(i,fv(f)%d(1))
                  else if (fn(2).eq.0)then
                    fn(2) = vn(i,fv(f)%d(1))
                  else
                    fn(3) = vn(i,fv(f)%d(1))
                  end if
                end if
              end do
              write (unitid,'(a, 3i)') '3', fn(:)-1
            end do
          end if
!- 3) Write Celltypes for 'Triangles'
          ! CELL TYPES 'Triangles'
          ! header for element (cell) types preceded by a blank line
          if (present(ellist)) then
            write (unitid,'(a,i20)') 'CELL_TYPES', size_flist
          else
            write (unitid,'(a,i20)') 'CELL_TYPES', numf
          end if
          ! write cell type for each volume element
          !  5   VTK_TRIANGLE            nodes 0..2 for vertices
          ! 22   VTK_QUADRATIC_TRIANGLE  3 between vertex 0 and 1
          !                              4 between        1     2
          !                              5 between        0     2
          if (present(ellist)) then
            do i=1,size(ellist)
              write (unitid,'(a)') '5'
            end do
          else
            do i=1,numf
              write (unitid,'(a)') '5'
            end do
          end if
!__
! 4) Release Memory:
          ok = destroyarrptr(fv)
!-
        case (10) !  CELLS 'Tetrahedrons'
!- 1) Write 'Tetrahedrons'
!     header for elements preceded by a blank line
          if (present(ellist)) then
            write (unitid,'(a,2i20)') 'CELLS', size(ellist), size(ellist)*5
          else
            write (unitid,'(a,2i20)') 'CELLS', numv, numv*5
          end if
          ! write volume element connectivity
          ! VTK-Indices are 0-offset. Thus the first point is point id 0
          if (present(ellist)) then
            do i=1,size(ellist)
               write (unitid,'(a, 4i)') '4', vn(:,ellist(i))-1
            end do
          else
            do i=1,numv
              write (unitid,'(a, 4i)') '4', vn(:,i)-1
            end do
          end if
!- 2) Write Celltypes for 'Tetrahedrons'
!     CELL TYPES Tetrahedrons
!     header for element (cell) types preceded by a blank line
!     for tet the vtk code is 10
          if (present(ellist)) then
            write (unitid,'(a,i20)') 'CELL_TYPES', size(ellist)
          else
            write (unitid,'(a,i20)') 'CELL_TYPES', numv
          end if
!     write cell type for each volume element
!     10   VTK_TETRA            nodes 0..3 for vertices
!     24   VTK_QUADRATIC_TETRA  4 between vertex 0 and 1
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
!-
        case default
          print*,'***WARNING: No CELLS are written. CELL_TYPE',cell_type,' not known.'
!-
      end select
!-
!__END.
!
      return
    end
