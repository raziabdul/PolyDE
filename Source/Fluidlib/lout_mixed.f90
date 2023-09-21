subroutine lout_mixed(gilt, eleinfo, eps, npdof)
use femtypes
use feminterface, only: getsetting, timestamp
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer(I4B) :: npdof
real(DP) :: eps
logical :: gilt, eleinfo
intent (in) :: gilt, eleinfo, npdof, eps
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
!            eps      = stopping criteria for linear solver.
!            npdof    = number of pressure dofs
!**********************************************************************************
!  local variables
!**********************************************************************************
integer(I4B) :: unitid, i, j, neldof, pdof
integer(I4B) :: unit2, unit3, unit4
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
!  open the file loes in the specified path
open (unitid,file=path(1:len_trim(path))//'solution', &
     &  form='unformatted',position='REWIND',action='WRITE')
!**********************************************************************************
! solution of velocity in x-direction is set to be default for output plot
call grglun(unit2)  
open (unit2,file=path(1:len_trim(path))//'solution_p', &
     &  form='unformatted',position='REWIND',action='WRITE')
!**********************************************************************************
! temporary removed to leave STOKES solver only for fluid simulation
!if (optmagn .OR. optener) then
!   call grglun(unit3)
!open (unit3,file=path(1:len_trim(path))//'solution_t', &
!     &  form='unformatted',position='REWIND',action='WRITE')
!end if
!**********************************************************************************
if (optstrm) then
   call grglun(unit4)
   open (unit4,file=path(1:len_trim(path))//'solution_str', &
        &  form='unformatted',position='REWIND',action='WRITE')
   open (50,file=path(1:len_trim(path))//'stream.out', &
     &  form='formatted',position='REWIND',action='WRITE')
end if
!**********************************************************************************
open (20,file=path(1:len_trim(path))//'solution.out', &
     &  form='formatted',position='REWIND',action='WRITE')
!**********************************************************************************
!  write the data to the file
write (unitid,err=1000) jd, gilt, eleinfo, n, ndof, omega, eps, 0._DP
write (unit2,err=1000)  npdof
write (20,*) npdof

if (eleinfo) then
   do j=1,n
      write (unitid,err=1000) ep(j,1)-1
!  number of Degree Of Freedom (dof) of this element
      !neldof=(ep(j)+1)*(ep(j)+2)/2
      neldof = ep(j,1)*(ep(j,1)+1)/2
      do i=1,neldof
         !if (i.le.pdof) write (unit2,err=1000) eg(j)%d(i)
         write (unitid,err=1000)  eg(j,1)%d(i)
      end do
   end do
end if

if (gilt) then
   ! write solution to file for 4 variables (u,v,p,T)
   if (optmagn .or. optener) then
      do j=1,npdof
         ! u,v
         write (unitid,err=1000) cmplx(unkno(2,j),0._DP,DPC), cmplx(unkno(3,j),0._DP,DPC)
         ! write pressure solution
         write (unit2,err=1000)  cmplx(unkno(1,j),0._DP,DPC)  
         ! write temperature solution
         !write (unit3,err=1000) cmplx(unkno(4,j),0._DP,DPC)
         ! matlab   
         write (20,301) j, (unkno(i,j), i=1,nvar)
         301 format(i4,4E17.8)
      end do

      ! write only velocity or temperature
      do j=npdof+1,ndof
         ! write velocity solution
         write (unitid,err=1000) cmplx(unkno(2:3,j),0._DP,DPC)
         ! write temperature solution
         !write (unit3,err=1000) cmplx(unkno(4,j),0._DP,DPC)
         ! matlab
         write (20,301) j, 0.0_DP, (unkno(i,j), i=2,nvar)
      end do
   else
   ! write solution to file for 3 variables (without coupling with energy equation)
      do j=1,npdof
      ! u,v
         write (unitid,err=1000) cmplx(unkno(2,j),0._DP,DPC), cmplx(unkno(3,j),0._DP,DPC)
      ! write pressure solution
         write (unit2,err=1000)  cmplx(unkno(1,j),0._DP,DPC)
      ! matlab
         write (20,300) j, (unkno(i,j), i=1,nvar)
         300 format(i4,3E17.8)
      end do
        
      ! write only velocity
      do j=npdof+1,ndof
      ! write velocity solution
         write (unitid,err=1000) cmplx(unkno(2:3,j),0._DP,DPC)
      ! matlab
         write (20,300) j, 0.0_DP, (unkno(i,j), i=2,nvar)
      end do
   end if

! write stream function
   if (optstrm) then
      do j = 1,ndof
         ! write to solution_str
         write(unit4,err=1000)  cmplx(phi(j),0._DP,DPC)
         ! write to stream.out
         write(50,*) j, phi(j)
      end do !j 
      close(unit4,err=1000)
      close(50)
   end if

else
   do j=1,ndof
      write (unitid,err=1000)  cmplx(0._DP,0._DP,DPC), cmplx(0._DP,0._DP,DPC)
      !write (unit3,err=1000)  cmplx(0._DP,0._DP,DPC)
      if (j .le. npdof) write (unit2,err=1000) cmplx(0._DP,0._DP,DPC)
   end do
end if

close(unitid,err=1000)
close(unit2,err=1000)
!if (optmagn .or. optener) close(unit3,err=1000)
close(20)

return
1000  print*,' an error occured while writing the solution vector'
stop

end subroutine lout_mixed