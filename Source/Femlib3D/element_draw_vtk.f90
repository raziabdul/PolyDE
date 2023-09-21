      subroutine element_draw_vtk(ellist)
      use feminterface, only: getsetting
      use feminterface3D, only: vtk_write_mesh
      use femtypes
      implicit none
      integer(I4B) :: ellist(:)
      intent (in) :: ellist
!
!  Subroutine to export Polyoder  into a vtk file (Simple Legacy Formats )
!
!  Local variables:
      integer (I4B) :: ierror, unitid
      character (len=200):: path
!
      external    ::    grglun
!

!  write data to file elements.vtk
      call grglun(unitid)
      call getsetting('PROJECTPATH',path)
!  open a file in the specified project path
      open (unitid,FILE=path(1:len_trim(path))//'elements.vtk',        &
            FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
      if (ierror .ne. 0) then
        print*, 'Cannot open file: ', path(1:len_trim(path))//'elements.vtk'
      end if
!
      call vtk_write_mesh(unitid,ellist)
!
      close (unitid)
!
      return
      end subroutine
