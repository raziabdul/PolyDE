subroutine solver_comprs(eps)
use feminterface, only: getsetting
use femtypes
use globalvariables
use fluidvariables
use fluidinterface
implicit none
real(DP), intent(in) :: eps
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
!    $Revision: 1.5 $
!    $Date: 2010/08/27 12:53:59 $
!    $Author: m_kasper $
!
!***********************************************************************
!
!     Solver for thermal flow problem using CBS algorithms
!
!***********************************************************************
! local variables
integer(I4B) :: idof, i, ios, k, unitid
integer(I4B) :: itime,intime, nstep
real(DP), pointer :: rhs2(:,:), rhs1(:,:)
real(DP), pointer :: residue(:)
character(len=200) :: path
logical :: converge
!***********************************************************************
! convergence history file
call getsetting('PROJECTPATH',path)
call grglun(unitid)
open(unitid,file=path(1:len_trim(path))//'convergence.out', &
     form='formatted',action='write',position='rewind',iostat=ios)

! unkno(nvar,ndof) will be set in main sub. "solver_fluid"  
allocate(unkn1(nvar,ndof), pres(ndof), pres1(ndof))
allocate(tt1(ndof), sound(ndof))
allocate(residue(nvar))

! convert variables to conservative form
call transform(.true.)

! time iterations starts

nstep = istep + 1
istep = ntime + istep

timeiter: do itime = nstep, istep

   intime = itime - nstep + 1
! set unkn1(nvar,ndof) to the solution from (itime-1)-th iteration
   unkn1(:,:) = unkno(:,:)    

! calculation of intermediate velocities
   call step1(itime, intime, rhs2) !optional , rhs1)

! calculation of the pressure (explicit solver only)
   call step2(rhs2) !optional ,rhs1)

! calculation of velocity and imposing b.c.s
   call step3(rhs2)

! calculation of energy
   call step4(rhs2)

! impose Dirichlet boundary conditions
   call boundary

! artificial viscosity (Lapidus)
   call artivisco
   call boundary

! compute pressure
   call getpres(pres1)

! variables smoothing scheme for viscous flow
   if (csmoo .GT. 0._DP) then
      call rsmooth
      call boundary  ! impost b.c. again
      call getpres(pres1)   
   end if
! update T and pressure
   do idof = 1,ndof
      tt1(idof) = gamma*(unkno(4,idof)-0.5_DP*(unkno(2,idof)**2+unkno(3,idof)**2)/unkno(1,idof)) &
                  /unkno(1,idof)
      pres(idof)= pres1(idof)
   end do

! convergence check 
      call checkconverge(eps,residue,converge)

   if (iwrite .GT. 0) then
      !iwrite = ntime/10
      k = mod(itime,iwrite)
      if (k .EQ. 0) then
         write(*,*) 'Pass iteration: ', itime
         write(*,*) (residue(i),i=1,nvar) 
         write(unitid,*) itime, (residue(i),i=1,nvar)
         100 format(i5,3(TR2,e14.7))
      end if
   end if

   if (converge) then
      ! obtain stead-state solution
      write(unitid,*) 'Steady-state solution obtained after',itime,' iterations'
      exit timeiter
   end if
end do timeiter ! itime
close(unitid)
deallocate(unkn1, pres1)

! compute drag force
!call getfdrag

! compute streamline
!call getstreamline

call transform(.false.)

end subroutine solver_comprs