      subroutine vtk_varmatplot
      use femtypes
      use feminterface,       only: fetchmatparameters, getsetting
      use feminterface3D,     only: findelement3D, xyz2lam, vtk_write_mesh
      use globalvariables3D, only : nod, numn, vn, vv, numv, dommat, dom
      use matconstants,       only: numparam
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
      integer (I4B)           :: i, ierror, unitid, elem
      real (DP)               :: point(3)
      real( DP), allocatable  :: list(:)
      character (len=200)     :: path, vtk_file
      logical                 :: opened, found, ok
      
      external :: grglun

      allocate (list(numparam))
!
!  write data to vtk file
      call getsetting('PROJECTPATH',path)
      vtk_file = 'USERMAT_PROPERTY'
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

      write (unitid,'(a)') 'SCALARS Matparameter double'
      write (unitid,'(a)') 'LOOKUP_TABLE default'
!  write field value at each node
      do i=1,numn
        point = nod(:,i)
        call findelement3D(point,nod,vn,vv,numv,elem,found)
        call fetchmatparameters(list, dommat(dom(elem)), point)
        write (unitid,'(E17.8E3)') real(list(5)-list(4))
      end do
!
      close(unitid)
      deallocate(list)
      end subroutine
      
      
subroutine export2matlab(datatype,val)
    use feminterface,        only: print_error, getsetting, putsetting
    use globalvariables3D,   only: numn, nnat, x, nod
    use femtypes
    implicit none
    real      (DP), optional    :: val(:,:)
    character (len=*)           :: datatype
!
!
!  Subroutine to export values for the solution to matlab
!
!  Input:
!     fieldtype          field quantity to compute (e.g. POTENTIAL, EY ...)
!     phi                phase / time instant to compute
!
!  Internal variables:
    integer   (I4B)             :: max_fp_i, max_and_i, unitid, n
    real      (DP)              :: fp_tol, geo_fac, fp_omega
    character (len=20)          :: solver
    character (len=50)          :: meshfile
    character (len=200)         :: path
!
    external    ::    grglun
!__
! Checks and preparations:
  call getsetting('PROJECTPATH',path)
  call getsetting('NONLINSOLVERTYPE',solver)
  call getsetting('NONLIN_TOLERANCE',fp_tol)
  call getsetting('MAX_ANDERSON_ITER',max_and_i)
  call getsetting('MAX_FIXPOINT_ITER',max_fp_i)
  call getsetting('FPRELAX_PARAMETER', fp_omega)
  call getsetting('MESHFILE',meshfile)
  call getsetting('GEOMETRY_FACTOR',geo_fac)
  
  select case (datatype)
    case ('HEADERS')
      ! Residual Header
      call grglun(unitid)
      open(unitid,file=path(1:len_trim(path))//'Residual.txt', form='formatted',position='REWIND',action='WRITE')
        write(unitid,*) 'Meshfile: ', meshfile
        write(unitid,*) 'Nonlinear Solver Type: ', solver
        write(unitid,*) 'Number of natures: ', nnat
        write(unitid,*) 'Geometry Factor ', geo_fac
        write(unitid,*) 'Nonlinear Tolerance' , fp_tol
        write(unitid,*) 'Fixed-Point Relaxation:', fp_omega
        write(unitid,*) '--------------------------------------------'   
        write(unitid,*) 'Iteration Step', 'FP Residual' 
        write(unitid,*) '--------------------------------------------'   
        write(unitid,*) 'DATA'
      close(unitid)
      
    case ('NODES')
      call grglun(unitid)
      open(unitid,file=path(1:len_trim(path))//'NodeCoordinates.txt', form='formatted',position='REWIND',action='WRITE')
        do n = 1,numn
          write(unitid,*) nod(1,n),nod(2,n),nod(3,n)
        end do
      close(unitid)
    
    case ('RESIDUAL')
      ! val(1,1) = Number of the current fp - iteration
      ! val(1,2) = Residual of the current fp - iteration
      if (.not. present(val)) call print_error(1, 'Writing Residuum to Matlab was called, but values were not given!')
      call grglun(unitid)
      open (unitid,file=path(1:len_trim(path))//'Residual.txt', form='formatted',position='APPEND',action='WRITE')
        write(unitid,*) val(1,1), val(1,2)
      close (unitid)
    
    case ('SOLUTION')
      call grglun(unitid)
      open (unitid,file=path(1:len_trim(path))//'LiveSolution.txt', form='formatted',position='REWIND',action='WRITE')    
      if (nnat .eq. 1) then 
        write(unitid,*) 'NATURE 1'
        do n = 1, numn
          write(unitid,*) real(x(n))
        end do
      elseif (nnat .eq. 2) then
        write(unitid,*) 'NATURE 1                 ', 'NATURE 2               '
        do n = 1, numn
          write(unitid,*) real(x(n)), real(x(numn+n))
        end do
      elseif (nnat .eq. 3) then
        write(unitid,*) 'NATURE 1                 ', 'NATURE 2               ', 'NATURE 3'
        do n = 1, numn
          write(unitid,*) real(x(n)), real(x(numn+n)), real(x(2*numn+n))
        end do
      end if
      close (unitid)

    case default
  
  end select
  

    
end subroutine