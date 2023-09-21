      subroutine louttrans(gilt, eleinfo, n, ndof, omega, epsgl, nnat, ep,     &
     &           eg, x, fem_accuracy)
      use feminterface, only: getsetting, timestamp
      use femtypes
      implicit none
      integer (I4B) :: n
      integer (I4B) :: ndof, nnat
      integer (I4B), optional :: ep(:,:)
      complex (DPC), optional :: x(:)
      type (ARRPTRI), optional :: eg(:,:)
      real (DP) :: omega, epsgl, fem_accuracy
      logical gilt, eleinfo
      intent (in) :: gilt, eleinfo, n, ndof, nnat, omega, epsgl, fem_accuracy
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
!    $Date: 2011/08/09 09:55:53 $
!    $Author: chvokas $
!
!  write the solution to the file solution
!
!  Input
!        gilt         =.true. if the solution is valid
!        eleinfo      =.true. if ep and eg should be written
!        n            number of elements
!        ndof         number of unknowns (number of dof)
!        omega        angular frequency
!        epsgl        accuracy of the solution in the solution of the linear system
!        nnat         number of natures (multiphysics)
!        ep           polynomial degree of elements
!        eg           degrees of freedom (global) for the elements
!        x            solution vector
!        fem_accuracy actual estimated global error of the FEM solution
!
!  local variables
      integer (I4B) unitid, i, j, neldof, inat
      real (DP) jd
      character (len=200) path
!
      external    ::    grglun
!
!  read in timestamp
      jd = timestamp()
!  fuer alle Knoten den Funktionswert
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
!
!  open the file loes in the specified path
      open (unitid,file=path(1:len_trim(path))//'solution',             &
     &  form='unformatted',position='REWIND',action='WRITE')
!
!  write the data to the file
!      write (unitid,err=1000) gilt, omega, epsgl, x(1:ndof)
      write (unitid,err=1000) jd, gilt, eleinfo, n, ndof, omega, epsgl, &    
     &                        fem_accuracy
!
      if (eleinfo) then
        write (unitid,err=1000) nnat
        do inat=1,nnat
          do j=1,n
            write (unitid,err=1000) ep(j,inat)
!  number of Degree Of Freedom (dof) of this element
            neldof=(ep(j,inat)+1)*(ep(j,inat)+2)/2
            do i=1,neldof
               write (unitid,err=1000) eg(j,inat)%d(i)
            end do
          end do
        end do
      end if
!
      if (gilt .and. present(x)) then
         do j=1,ndof
            write (unitid,err=1000) x(j)
         end do
      else
         do j=1,ndof
            write (unitid,err=1000) cmplx(0._DP,0._DP,DPC)
         end do
      end if
!
      close (unitid,err=1000)
      return
1000  print*,' an error occured while writing the solution vector'
      stop
      end subroutine louttrans