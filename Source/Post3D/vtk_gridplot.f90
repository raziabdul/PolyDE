      subroutine vtk_scalargridplot()
      use feminterface, only: getsetting, getpostsetting
      use feminterface3D, only: fieldquantity3d, findelement3D, xyz2lam
      use femtypes
      use globalvariables3D, only : nod, numv, vn, vv
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
!    $Revision: 1.1 $
!    $Date: 2014/10/30 05:02:37 $
!    $Author: raziabdul $
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
      integer (I4B) :: i, ierror, unitid, elem, nx, ny, nz, j, k
      real (DP) :: xmin, xmax, ymin, ymax, zmin, zmax, startz
      real (DP) :: phi, point(3), lambda(4)
      real (DP) , allocatable :: xgrid(:), ygrid(:), zgrid(:)
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
!  write vtk header
!  make sure # starts in column 1
        write (unitid,'(a)') '# vtk DataFile Version 3.0'
        write (unitid,'(a)') 'Plot on rectilinear grid'
        write (unitid,'(a)') 'ASCII'
        write (unitid,'(a)') 'DATASET RECTILINEAR_GRID'
!  Grid lines
        call getpostsetting('NX',nx)
        nx = max(nx,2)
        call getpostsetting('NY',ny)
        ny = max(ny,2)
        call getpostsetting('NZ',nz)
        nz=1
        write (unitid,'(a,3i4)') 'DIMENSIONS ', nx, ny, nz
!  fetch extends
        xmin = minval(nod(1,:))
        xmax = maxval(nod(1,:))
        ymin = minval(nod(2,:))
        ymax = maxval(nod(2,:))
        zmin = minval(nod(3,:))
        zmax = maxval(nod(3,:))
!
        allocate(xgrid(nx), ygrid(ny), zgrid(nz))
        write (unitid,'(a,i4,a)') 'X_COORDINATES ', nx,' float'
        do i=1,nx
          xgrid(i) = xmin + (xmax-xmin)/real(nx-1,DP)*real(i-1,DP)
        end do
        write (unitid,'(g13.6)') xgrid
        write (unitid,'(a,i4,a)') 'Y_COORDINATES ', ny,' float'
!
        do i=1,nx
          ygrid(i) = ymin + (ymax-ymin)/real(ny-1,DP)*real(i-1,DP)
        end do
        write (unitid,'(g13.6)') ygrid
!
        write (unitid,'(a,i4,a)') 'Z_COORDINATES ', nz,' float'
        call getpostsetting('STARTZ',startz)
        zgrid(1) = startz
        write (unitid,'(g13.6)') zgrid
!

!  header for scalar data for each node preceded by a blank line
        write (unitid,'(a,i20)') 'POINT_DATA', nx*ny*nz
      end if

      call getpostsetting('PHI',phi)
      call getpostsetting('FIELDTYPE',fieldtype)

      write (unitid,'(a,a)') 'SCALARS '//trim(adjustl(fieldtype)),' double'
      write (unitid,'(a)') 'LOOKUP_TABLE default'
!  write field value at each node
      do k=1,nz
        do j=1,ny
          do i=1,nx
            point = (/ xgrid(i), ygrid(j), zgrid(k) /)
            call findelement3D(point,nod,vn,vv,numv,elem,found)
            if (.not.found) then
!              print *,"Point is not within the computational domain"
              zs = 0._DP
            else
!  get barycentric coordinates lambda of point in the element
              call xyz2lam(lambda, elem, point, nod, vn)
              call fieldquantity3d(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
              if (.not.ok) then
                print*, ' failed to evalute fieldquantity'
              end if
            end if
            write (unitid,'(g13.5)') real(zs)
          end do
        end do
      end do
!
      end subroutine

