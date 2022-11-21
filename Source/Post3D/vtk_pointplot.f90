      subroutine vtk_scalarpointplot()
      use feminterface, only: getsetting, getpostsetting
      use feminterface3D, only: fieldquantity3d, findelement3D, xyz2lam, vtk_write_mesh
      use femtypes
      use globalvariables3D, only : nod, numn , numv, vn, vv
      implicit none
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
!    $Revision: 1.3 $
!    $Date: 2014/11/04 15:27:01 $
!    $Author: juryanatzki $
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
      integer (I4B) :: i, ierror, unitid, elem
      real (DP) :: phi, point(3), lambda(4)
      complex (DPC) :: zs
      character (len=200):: path, vtk_file
      character (len=20) :: fieldtype
      character (len=50) :: descriptor
      character (len=12) :: unit
      logical :: opened, found, ok

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
!  header for scalar data for each node preceded by a blank line
        write (unitid,'(a,i20)') 'POINT_DATA', numn
      end if

      call getpostsetting('PHI',phi)
      call getpostsetting('FIELDTYPE',fieldtype)

      write (unitid,'(a,a)') 'SCALARS '//trim(adjustl(fieldtype)),' double'
      write (unitid,'(a)') 'LOOKUP_TABLE default'
!  write field value at each node
      do i=1,numn
        point = nod(:,i)
        call findelement3D(point,nod,vn,vv,numv,elem,found)
        if (.not.found) then
          print *,"Point is not within the computational domain"
        else
!  get barycentric coordinates lambda of point in the element
          call xyz2lam(lambda, elem, point, nod, vn)
        end if
        call fieldquantity3d(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
        if (.not.ok) then
          print*, ' failed to evalute fieldquantity'
        end if
        write (unitid,'(E17.8E3)') real(zs)
      end do
!
      end subroutine



      subroutine vtk_vectorpointplot()
      use feminterface, only: getsetting, getpostsetting
      use feminterface3D, only: fieldquantity3d, findelement3D, xyz2lam, vtk_write_mesh
      use femtypes
      use globalvariables3D, only : nod, numn, numv, vn, vv
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
      integer (I4B) :: i, ierror, unitid, elem
      real (DP) :: phi, point(3), lambda(4)
      complex (DPC) :: zs(3)
      character (len=200):: path, vtk_file
      character (len=20) :: fieldtype
      character (len=50) :: descriptor
      character (len=12) :: unit
      logical :: opened, found, ok

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
!  header for scalar data for each node preceded by a blank line
        write (unitid,'(a,i20)') 'POINT_DATA', numn
      end if

      call getpostsetting('PHI',phi)
      call getpostsetting('FIELDTYPE',fieldtype)

      write (unitid,'(a,a)') 'VECTORS '//trim(adjustl(fieldtype)),' double'
!  write field value at each node
!
      do i=1,numn
        point = nod(:,i)
        call findelement3D(point,nod,vn,vv,numv,elem,found)
        if (.not.found) then
          print *,"Point is not within the computational domain"
        else
!  get barycentric coordinates lambda of point in the element
          call xyz2lam(lambda, elem, point, nod, vn)
        end if
        call fieldquantity3d(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
        if (.not.ok) then
          print*, ' failed to evalute fieldquantity'
        end if
        write (unitid,'(3g13.5)') real(zs(1)),real(zs(2)), real(zs(3))
      end do
!
      end subroutine



      subroutine vtk_tensorpointplot()
      use feminterface, only: getsetting, getpostsetting
      use feminterface3D, only: fieldquantity3d, findelement3D, xyz2lam, vtk_write_mesh
      use femtypes
      use globalvariables3D, only : nod, numn, numv, vn, vv
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
      integer (I4B) :: i, ierror, unitid, elem
      real (DP) :: phi, point(3), lambda(4)
      complex (DPC) :: zs(3,3)
      character (len=200):: path, vtk_file
      character (len=20) :: fieldtype
      character (len=50) :: descriptor
      character (len=12) :: unit
      logical :: opened, found, ok

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
!  header for scalar data for each node preceded by a blank line
        write (unitid,'(a,i20)') 'POINT_DATA', numn
      end if

      call getpostsetting('PHI',phi)
      call getpostsetting('FIELDTYPE',fieldtype)

      write (unitid,'(a,a)') 'TENSORS '//trim(adjustl(fieldtype)),' double'
!  write field value at each node
!
      do i=1,numn
        point = nod(:,i)
        call findelement3D(point,nod,vn,vv,numv,elem,found)
        if (.not.found) then
          print *,"Point is not within the computational domain"
        else
!  get barycentric coordinates lambda of point in the element
          call xyz2lam(lambda, elem, point, nod, vn)
        end if
        call fieldquantity3d(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
        if (.not.ok) then
          print*, ' failed to evalute fieldquantity'
        end if
        write (unitid,'(3g13.5)') real(zs(1,1)),real(zs(1,2)), real(zs(1,3))
        write (unitid,'(3g13.5)') real(zs(2,1)),real(zs(2,2)), real(zs(2,3))
        write (unitid,'(3g13.5)') real(zs(3,1)),real(zs(3,2)), real(zs(3,3))
      end do
!
      end subroutine



