subroutine solver_fluid
use feminterface, only: getsetting
use femtypes
use globalvariables
use fluidvariables
use fluidinterface, only: lout_fluid
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
!    $Revision: 1.5 $
!    $Date: 2006/07/07 22:34:04 $
!    $Author: r_abdul $
!
!*****************************************************************************
!
!     Solver for thermal flow problem using CBS algorithms
!
!*****************************************************************************
! Important Remarks:
!  1. Read Instruction at \\Mstprt\Polyde\Fluid Solver\
!     Instruction for solvers.doc
!*****************************************************************************
! local variables
!*****************************************************************************
real(DP) :: eps
!*****************************************************************************
! input parameters of fluid, ex. theta, timestep, niter
call inputfluid

!if (optmagn) then
!   call magneticinput
!end if

! no of primary variables
if (.NOT. optcomprs) then
   if (optener) then
      nvar = 4
   else
      nvar = 3
   end if
else
   nvar = 4
end if

! read multivariable boundary conditions
call readnetin_mn

call getsetting('EPS',eps)

allocate(unkno(nvar,ndof))
unkno = 0._DP

! element length, boundary-located element sides.
call geteleminfo
call getmassmat

! setting preliminary value
call presettings

! solver for compressible or incompressible flow
if (optcomprs) then
!   call solver_comprs(eps)
   write(*,*) 'Compressible solver not mounted yet!!!'
else
   if (optstokes) then
      ! Stokes solver for equal-order (see important remark in elementmatrix_stokes.f90   )
      write(*,*) 'solver for Stokes flow'
      call solver_stokes
   else
      write(*,*) 'solver for incompressible flow'
      call solver_incomprs
   end if
end if

! write solution to seperate output file for each variables
call lout_fluid(.true.,.true., eps)

end subroutine solver_fluid

