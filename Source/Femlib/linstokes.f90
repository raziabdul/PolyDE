      subroutine linstokes(ok, jdmesh, jdsolu, gilt, eleinfo, n, ndof, omega, &
     &  epsgl, ep, eg, xfluid, fem_accuracy)
      use feminterface, only: getsetting, juld2cal
      use femtypes
      implicit none
      integer (I4B) ndof, n !, npdof
      integer (I4B), pointer :: ep(:,:)
      complex (DPC), pointer :: xfluid(:,:)
      type (ARRPTRI), pointer :: eg(:,:)
      real (DP) omega, epsgl, jdmesh, jdsolu, fem_accuracy
      logical ok, gilt, eleinfo
      intent (in) :: jdmesh, n
      intent (out) :: ok, jdsolu, gilt, eleinfo, ndof, omega !, npdof
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
!    $Revision: 1.4 $
!    $Date: 2008/08/20 12:47:58 $
!    $Author: m_kasper $
!
!*******************************************************************************
!
!  read in the solution vector from file solution of STOKES solver
!
!*******************************************************************************
! remarks: Since the mixed method is used, the orders of shape functions of
!          velocity, temperature (if solved) and pressure are different in the
!          way that ep(velocity) = ep(temperature) = ep(pressure) + 1
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
!        ep           polynomial degree of elements
!        eg           degrees of freedom (global) for the elements
!        x            solution vector
!        fem_accuracy actual estimated global error of the FEM solution
!*******************************************************************************
!  local variables
!*******************************************************************************
      integer (I4B) unitid, i, j, neldof, nstored, npdof
      integer (I4B) year, mon, day, hour, minute, second
      integer (I4B) unitid2, unitid3
      character(len=3) dayofweek, month
      character (len=30) datstr
      character (len=200) path
      logical exist, exist_t, exist_p
!      
      external :: grglun
!
      call getsetting('PROJECTPATH',path)
      inquire (file=path(1:len_trim(path))//'solution',exist=exist)
      if (.not. exist) then
        print*,'***** solution data file not present'
        print*,'***** call ignored'
        ok=.false.
        return
      end if

      inquire (file=path(1:len_trim(path))//'solution_t',exist=exist_t)
      inquire (file=path(1:len_trim(path))//'solution_p',exist=exist_p)

      if (exist_t) then
         call grglun(unitid2)
         open (unitid2,file=path(1:len_trim(path))//'solution_t',                 &
     &         form='unformatted',status='old',position='REWIND',action='READ')
      end if

      if (exist_p) then
         call grglun(unitid3)
         open (unitid3,file=path(1:len_trim(path))//'solution_p',                 &
     &         form='unformatted',status='old',position='REWIND',action='READ')
      end if

      call grglun(unitid)
!
!  open the file under the specified path
      open (unitid,file=path(1:len_trim(path))//'solution',                 &
     &  form='unformatted',status='old',position='REWIND',action='READ')
!
!  read in the data
      read (unitid,end=1000,err=1000) jdsolu, gilt, eleinfo, nstored, ndof, &
                                      omega, epsgl, fem_accuracy
      read (unitid3,end=1000,err=1000) npdof

! temporary
      ndof = npdof

!
      if (nstored.gt.0) then
         if ((nstored .eq. n) .and. (jdsolu .ge. jdmesh)) then
            if (eleinfo) then
               if (associated(ep)) then
                  deallocate(ep)
               end if
               if (associated(eg)) then
                  do j=1,n
                     deallocate(eg(j,1)%d)
                  end do
                  deallocate(eg)
               end if
               allocate(ep(n,1))
               allocate(eg(n,1))
               do j=1,n
                  read (unitid,end=1000,err=1000) ep(j,1)
!  number of Degree Of Freedom (dof) of this element
                  neldof=(ep(j,1)+1)*(ep(j,1)+2)/2
!  allocate array of global dof of the element
                  allocate(eg(j,1)%d(neldof))
                  do i=1,neldof
                     read (unitid,end=1000,err=1000) eg(j,1)%d(i)
                  end do
               end do
            end if
!
            if (associated(xfluid)) then
               deallocate(xfluid)
            end if

            if (exist_p .and. exist_t) then
               allocate(xfluid(4,ndof))
            else if (exist_p) then
               allocate(xfluid(3,ndof))
            else
               allocate(xfluid(2,ndof))
            end if
!  read in the solution (stored for each DOF)
            do j=1,ndof
               read (unitid,end=1000,err=1000) (xfluid(i,j),i=1,2)
            end do
            
            if (exist_t) then
               do j=1,ndof
                  read (unitid2,end=1000,err=1000) xfluid(4,j)
               end do
            end if

            if (exist_p) then
               do j=1,npdof
                  read (unitid3,end=1000,err=1000) xfluid(3,j)
               end do
            end if
!
         else
            print*
            print*,'**** solution contains vector for a mesh of: ',nstored,' elements'
            print*,'**** actual mesh size is: ',n,' elements; ignoring the solution'
            call juld2cal(jdmesh,year,mon,day,month,dayofweek,hour,minute,second)
            write(datstr,'(A3,'', '',I2.2,''. '',A3,''. '',I4,          &
     &        '' at '',I2,'':'',I2.2,'':''I2.2)')                       &
     &        dayofweek,day,month,year,hour,minute,second
            print*,'date of mesh is:  ',datstr
            call juld2cal(jdsolu,year,mon,day,month,dayofweek,hour,minute,second)
            write(datstr,'(A3,'', '',I2.2,''. '',A3,''. '',I4,          &
     &        '' at '',I2,'':'',I2.2,'':''I2.2)')                       &
     &        dayofweek,day,month,year,hour,minute,second
            print*,'date of solution is:  ',datstr
            print*
         end if
      end if
      close (unitid,err=1000)
      if (exist_t) close (unitid2,err=1000)
      if (exist_p) close (unitid3,err=1000)
!
      ok=.true.
      return
1000  print*,' an error occured while reading the solution file'
      stop
      end subroutine linstokes
