subroutine lout_fluid(gilt, eleinfo, eps)
use femtypes
use feminterface, only: getsetting, timestamp
use globalvariables
use fluidvariables
use fluidinterface
implicit none
real(DP) :: eps
logical gilt, eleinfo
intent (in) :: gilt, eleinfo, eps
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
!    $Revision: 1.8 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!**********************************************************************************
!  write the solution to the file solution
!**********************************************************************************
!  Input
!            gilt     =.true. if the solution is valid
!            eleinfo  =.true. if ep and eg should be written
!            eps      = stopping criteria for steady-state solution
!**********************************************************************************
!  local variables
!**********************************************************************************
integer(I4B) :: unitid, i, j, neldof
integer(I4B) :: unit2,unit4
real(DP) :: jd
character(len=200) :: path
!      
external :: grglun

!**********************************************************************************
!  read in timestamp
jd = timestamp()
!  fuer alle Knoten den Funktionswert
call getsetting('PROJECTPATH',path)
!**********************************************************************************
call grglun(unitid)
!  open the file in the specified path
open (unitid,file=path(1:len_trim(path))//'solution', &
     &  form='unformatted',position='REWIND',action='WRITE')
!**********************************************************************************
if (optcomprs .OR. optener) then
   call grglun(unit4)
!  open the file in the specified path
   open (unit4,file=path(1:len_trim(path))//'solution_t', &
     &  form='unformatted',position='REWIND',action='WRITE')
end if
!**********************************************************************************
if (optstrm) then
   call grglun(unit2)
   open (unit2,file=path(1:len_trim(path))//'solution_str', &
        &  form='unformatted',position='REWIND',action='WRITE')
   open (50,file=path(1:len_trim(path))//'stream.out', &
     &  form='formatted',position='REWIND',action='WRITE')
end if
! for Matlab plot
open (20,file=path(1:len_trim(path))//'solution.out', &
     &  form='formatted',position='REWIND',action='WRITE')
!**********************************************************************************
!  write the data to the file
!  jdsolu, gilt, eleinfo, nstored,ndof, omega, epsgl, fem_accuracy
write(unitid,err=1000) jd, gilt, eleinfo, n, ndof, omega, eps, 0._DP
write(*,*) jd, gilt, eleinfo, n, ndof, omega, eps, 0._DP

! write the number of step which has been run so far.
write(20,*) istep

if (eleinfo) then
! for fluid variables
   do j=1,n
      write (unitid,err=1000) ep(j,1)
   !  number of Degree Of Freedom (dof) of this element
      neldof=(ep(j,1)+1)*(ep(j,1)+2)/2
      do i=1,neldof
         write (unitid,err=1000) eg(j,1)%d(i)
      end do
   end do
end if

if (gilt) then
   do j=1,ndof
   ! write p,u,v
      write (unitid,err=1000) cmplx(unkno(1,j),0._DP,DPC), cmplx(unkno(2,j),0._DP,DPC), & 
&                             cmplx(unkno(3,j),0._DP,DPC)
   end do !j

   if (optcomprs .OR. optener) then
      do j = 1,ndof
         ! write temperature .bsw. internal energy
         write (unit4,err=1000)  cmplx(unkno(4,j),0._DP,DPC)
         
         ! write output for matlab
         write (20,100) j, (unkno(i,j), i=1,nvar)         
      end do
      100 format(i4,4E20.10)
   else
      ! write output for matlab only for p, u, v
      do j = 1,ndof
         write (20,200) j, (unkno(i,j), i=1,nvar)
      end do
      200 format(i4,3E20.10)
   end if
! write stream function
   if (optstrm) then
      do j = 1,ndof
         ! write to solution_str
         write (unit2,err=1000)  cmplx(phi(j),0._DP,DPC)
         ! write to stream.out
         write(50,*) j, phi(j)
      end do !j 
   end if
   close(50)
else
   do j=1,ndof
      write (unitid,err=1000) cmplx(0._DP,0._DP,DPC)
      if (optstrm) write (unit2,err=1000)  cmplx(0._DP,0._DP,DPC)
      if (optcomprs .OR. optener) write (unit4,err=1000)  cmplx(0._DP,0._DP,DPC)
   end do
end if
close(unitid,err=1000)
if (optstrm) close(unit2,err=1000)
if (optcomprs .OR. optener) close(unit4,err=1000)

if (optener) then
   write(20,110) 'Average Nusselt Number = ', nuaver
   110 format(1x,a25,es17.10)
end if

close(20)
return
1000  print*,' an error occured while writing the solution vector'
stop

end subroutine lout_fluid
