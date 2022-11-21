subroutine step1_incomprs(rhs)
use femtypes
use globalvariables
use feminterface
use fluidvariables
use fluidinterface
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
!    $Revision: 1.4 $
!    $Date: 2006/07/07 22:34:04 $
!    $Author: r_abdul $
!
!###########################################################################################
!
!     Calculate Intermediate Velocity
!
!###########################################################################################
real(DP), pointer :: rhs(:,:)
!*******************************************************************************************
! local variables
!*******************************************************************************************
integer(I4B) :: idof
!real(DP) :: invmmat
!*******************************************************************************************
allocate(rhs(nvar,ndof))
rhs = 0._DP

! compute timestep
call timestep_incomprs

! compute rhs (advection and diffusion terms)
call getrhsstep1_incomprs(rhs)

! update the solution
do idof = 1,ndof
   unkno(2,idof)= unkno(2,idof) + dmmat(idof)*deltp(idof)*rhs(2,idof)
   unkno(3,idof)= unkno(3,idof) + dmmat(idof)*deltp(idof)*rhs(3,idof)
end do !idof

!*******************************************************************************************
end subroutine step1_incomprs
