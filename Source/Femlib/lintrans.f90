      subroutine lintrans(ok, jdmesh, jdsolu, gilt, eleinfo, n, ndof, omega, &
     &               epsgl, nnat, ep, eg, x, fem_accuracy)
      use feminterface, only: getsetting, juld2cal
      use femtypes
      implicit none
      integer (I4B) :: ndof, n, nnat
      integer (I4B), pointer :: ep(:,:)
      complex (DPC), pointer :: x(:)
      type (ARRPTRI), pointer :: eg(:,:)
      real (DP) :: omega, epsgl, jdmesh, jdsolu, fem_accuracy
      logical ok, gilt, eleinfo
      intent (in) :: jdmesh, n, nnat
      intent (out) :: ok, jdsolu, gilt, eleinfo, ndof, omega
      intent (out) :: epsgl, fem_accuracy
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
!  read in the solution vector from file solution
!
!  Input:
!        n            number of elements of the actual mesh
!  Output:
!        ok           =.false. if an error occured
!        jdsolu       timestamp of solution
!        gilt         =.true. if the solution is valid
!        eleinfo      =.true. if ep and eg has been written
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
!
      integer (I4B) unitid, i, j, neldof, nstored, inat, nfnat
      integer (I4B) year, mon, day, hour, minute, second
      character(len=3) dayofweek, month
      character (len=30) datstr
      character (len=200) path
      logical exist
!
      call getsetting('PROJECTPATH',path)
      inquire (file=path(1:len_trim(path))//'solution',exist=exist)
      if (.not. exist) then
        print*,'***** solution data file not present'
        print*,'***** call ignored'
        ok=.false.
        return
      end if
      call grglun(unitid)
!
!  open the file under the specified path
      open (unitid,file=path(1:len_trim(path))//'solution',             &
     &  form='unformatted',status='old',position='REWIND',action='READ')
!
!  read in the data
      read (unitid,end=1000,err=1000) jdsolu, gilt, eleinfo, nstored,   &
     &                                ndof, omega, epsgl, fem_accuracy
!
      read (unitid,end=1000,err=1000) nfnat
!
      if (nstored.gt.0) then
        if ((nstored .eq. n) .and. (jdsolu .ge. jdmesh) .and.           &
     &                             (nnat.eq.nfnat)) then
          if (eleinfo) then
            if (associated(ep)) then
              deallocate(ep)
            end if
            if (associated(eg)) then
              do inat=1,nnat
                do j=1,n
                  deallocate(eg(j,inat)%d)
                end do
              end do
              deallocate(eg)
            end if
            allocate(ep(n,nnat))
            allocate(eg(n,nnat))
            do inat=1,nnat
              do j=1,n
                read (unitid,end=1000,err=1000) ep(j,inat)
!  number of Degree Of Freedom (dof) of this element
                neldof=(ep(j,inat)+1)*(ep(j,inat)+2)/2
!  allocate array of global dof of the element
                allocate(eg(j,inat)%d(neldof))
                do i=1,neldof
                  read (unitid,end=1000,err=1000) eg(j,inat)%d(i)
                end do
              end do
            end do
          end if
!
          if (associated(x)) then
            deallocate(x)
          end if
          allocate(x(ndof))
!  read in the solution (stored for each DOF)
          do j=1,ndof
            read (unitid,end=1000,err=1000) x(j)
          end do
!
        else
          print*
          write (*,100) '**** solution contains vector for a mesh of: ',nstored,' elements and ',nfnat,' natures'
100       format( A, I8, A, I3,A)
          write (*,100) '**** actual mesh size has: ',n,' elements and ',nnat,' natures'
          print*,'**** ignoring the solution'
          call juld2cal(jdmesh,year,mon,day,month,dayofweek,hour,minute,second)
          write(datstr,'(A3,'', '',I2.2,''. '',A3,''. '',I4,          &
     &      '' at '',I2,'':'',I2.2,'':''I2.2)')                       &
     &      dayofweek,day,month,year,hour,minute,second
          print*,'date of mesh is:  ',datstr
          call juld2cal(jdsolu,year,mon,day,month,dayofweek,hour,minute,second)
          write(datstr,'(A3,'', '',I2.2,''. '',A3,''. '',I4,          &
     &      '' at '',I2,'':'',I2.2,'':''I2.2)')                       &
     &      dayofweek,day,month,year,hour,minute,second
          print*,'date of solution is:  ',datstr
          print*
          gilt=.false.
        end if
      end if
      close (unitid,err=1000)
!
      ok=.true.
      return
1000  print*,' an error occured while reading the solution file'
      stop
      end subroutine lintrans
