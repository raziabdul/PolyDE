subroutine solver_incomprs
use feminterface, only: getsetting
use femtypes
use fluidvariables
use fluidinterface
use globalvariables
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
!    $Revision: 1.7 $
!    $Date: 2010/08/27 12:56:07 $
!    $Author: m_kasper $
!
!***********************************************************************
!
!     Solver for incompressible flow problem using CBS algorithms
!
!***********************************************************************
! 2 Options for step 2
! - Semi-implicit method
! - Explicit method with artificial compressibility
!***********************************************************************
! local variables
integer(I4B) :: itime, nstep
integer(I4B) :: unitid, i, k, ios
integer (I4B), pointer :: ia(:), ja(:)
real(DP), pointer :: lower(:), diag(:) !upper(:)
real(DP), pointer :: rhs(:,:)
real(DP), pointer :: residue(:)
real(DP) :: eps
character(len=200) :: path
logical :: converge, implc
!
external    ::    grglun
!***********************************************************************
! unkno(nvar,ndof) will be set in main sub. "solver_fluid"
allocate(unkn1(nvar,ndof)) 
allocate(deltp(ndof),delte(n))
allocate(residue(nvar))
! assign implc from theta2
implc = theta(2) .GT. 0    

! implicit solver
if (implc) then
   ! calculate stiffness matrix of pressure (only once)
   call pstiff(lower,diag,ia,ja)
end if

! convergence history file
call getsetting('PROJECTPATH',path)
call grglun(unitid)
open(unitid,file=path(1:len_trim(path))//'convergence.out', &
     form='formatted',action='write',position='rewind',iostat=ios)

! stopping criteria for steady-state
call getsetting('EPS',eps)

! time iterations starts
nstep = istep + 1
istep = ntime + istep

timeiter: do itime = nstep, istep

! set unkn1(nvar,ndof) to the solution from (itime-1)-th iteration
   unkn1(:,:) = unkno(:,:)    

! calculation of intermediate velocities
   call step1_incomprs(rhs)

! impost b.c.   
   call boundary_incomprs(0)

! calculation of the pressure
   if (implc) then
      call step2_implicit(lower,diag,ia,ja)
   else
      call step2_incomprs(rhs)
      call boundary_incomprs(2)
   end if

! calculation of velocity
   call step3_incomprs(rhs)

! impose Dirichlet boundary conditions
   call boundary_incomprs(0)

   if (optener) then
! calculation of temperature
      call step4_incomprs(rhs)

! impose Dirichlet boundary conditions
      call boundary_incomprs(1)
   end if

! convergence check
   call checkconverge(eps,residue,converge)

   if (iwrite .GT. 0) then
      !iwrite = ntime/10
      k = mod(itime,iwrite)
      if (k .EQ. 0) then
         write(*,*) 'Pass iteration: ', itime
         write(*,*) (residue(i),i=1,nvar) 
         if (optener) then
            write(unitid,100) itime, (residue(i),i=1,nvar)
            100 format(i5,4(TR1,e14.6))
         else
            write(unitid,200) itime, (residue(i),i=1,nvar)
            200 format(i5,3(TR1,e14.6))
         end if
         ! output at a certain time-step
!         call lout_fluid(.true.,.true., eps)         
      end if
   end if

   if (converge) then
      ! obtain stead-state solution      
      write(*,*) 'Steady-state solution obtained after',itime,' iterations'
      write(unitid,*) 'Steady-state solution obtained after',itime,' iterations'
      write(unitid,*) (residue(i), i=1,nvar)
      istep = itime
      exit timeiter
   end if
end do timeiter ! itime
close(unitid)
deallocate(unkn1)

if (optener) then
   call nusselt
end if

if (optstrm) then
! compute streamline function
   call streamfunction(ia,ja)
end if
deallocate(ia,ja)

end subroutine solver_incomprs
