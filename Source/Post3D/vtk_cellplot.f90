      subroutine vtk_scalarcellplot()
      use feminterface, only: getsetting, getpostsetting
      use feminterface3D, only: fieldquantity3d, vtk_write_mesh
      use femtypes
      use globalvariables3D, only : numv
      implicit none
!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 Institute for Micro Systems Technology,
!                       Hamburg University of Technology.
!
!    This file is part of PolyDE.
!
!    PolyDE is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by the Free Software Foundation; either version 2, or (at your
!    option) any later version.
!
!    PolyDE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!    MA 02110-1301 USA.
!
!------------------------------------------------------------------------------
!    $Revision: 1.3 $
!    $Date: 2014/08/22 11:07:36 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!
!
!  Subroutine to export values for given fieldtype and phase as an
!  unstructured grid output according to the vtk "classic" format.
!
!  Input:
!     fieldtype          field quantity to compute (e.g. POTENTIAL, EY ...)
!     phi                phase / time instant to compute

!  Internal variables:
      integer (I4B) :: ierror, unitid, elem
      real (DP) :: phi, lambda(4)
      complex (DPC) :: zs
      character (len=200):: path, vtk_file
      character (len=20) :: fieldtype
      character (len=50) :: descriptor
      character (len=12) :: unit
      logical :: opened, ok

!
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
          print*, 'Cannot open file: ', path(1:len_trim(path))//vtk_file(1:len_trim(vtk_file))//'.vtk'
        end if
!
        call vtk_write_mesh(unitid)
!
!  header for scalar data for each cell preceded by a blank line
        write (unitid,'(a,i20)') 'CELL_DATA', numv
      end if

      call getpostsetting('PHI',phi)
      call getpostsetting('FIELDTYPE',fieldtype)

      write (unitid,'(a,a)') 'SCALARS '//trim(adjustl(fieldtype)),' double'
      write (unitid,'(a)') 'LOOKUP_TABLE default'
!  write field value at each cell
      do elem=1,numv
!  center point
!  get barycentric coordinates lambda of point in the element
        lambda = (/0.25_DP,0.25_DP,0.25_DP,0.25_DP/)
        call fieldquantity3d(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
        if (.not.ok) then
          print*, ' failed to evalute fieldquantity'
        end if
        write (unitid,'(g13.5)') real(zs)
      end do
!
      end subroutine



      subroutine vtk_vectorcellplot()
      use feminterface, only: getsetting, getpostsetting
      use feminterface3D, only: fieldquantity3d, vtk_write_mesh
      use femtypes
      use globalvariables3D, only : numv
      implicit none
!
!
!  Subroutine to export values for given fieldtype and phase as an
!  unstructured grid output according to the vtk "classic" format.
!
!  Input:
!     fieldtype          field quantity to compute (e.g. POTENTIAL, EY ...)
!     phi                phase / time instant to compute

!  Internal variables:
      integer (I4B) :: ierror, unitid, elem
      real (DP) :: phi, lambda(4)
      complex (DPC) :: zs(3)
      character (len=200):: path, vtk_file
      character (len=20) :: fieldtype
      character (len=50) :: descriptor
      character (len=12) :: unit
      logical :: opened, ok

!
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
          print*, 'Cannot open file: ', path(1:len_trim(path))//vtk_file(1:len_trim(vtk_file))//'.vtk'
        end if
!
        call vtk_write_mesh(unitid)
!
!  header for scalar data for each cell preceded by a blank line
        write (unitid,'(a,i20)') 'CELL_DATA', numv
      end if

      call getpostsetting('PHI',phi)
      call getpostsetting('FIELDTYPE',fieldtype)

      write (unitid,'(a,a)') 'VECTORS '//trim(adjustl(fieldtype)),' double'
!  write field value at each cell
      do elem=1,numv
!  center point
!  get barycentric coordinates lambda of point in the element
        lambda = (/0.25_DP,0.25_DP,0.25_DP,0.25_DP/)
        call fieldquantity3d(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
        if (.not.ok) then
          print*, ' failed to evalute fieldquantity'
        end if
        write (unitid,'(3g13.5)') real(zs(1)),real(zs(2)), real(zs(3))
      end do
!
      end subroutine
