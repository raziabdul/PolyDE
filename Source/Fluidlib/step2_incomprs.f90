subroutine step2_incomprs(rhs)
use femtypes
use fluidvariables
use fluidinterface
use globalvariables
implicit none
real(DP), pointer :: rhs(:,:) ! intent(inout)
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
!##########################################################################################
!
!     Solve pressure equation in explicit form
!
!##########################################################################################
! local variables
!******************************************************************************************
integer(I4B) :: idof
real(DP) :: invmmat
real(DP), pointer :: beta(:)
!******************************************************************************************
! compute artificial compressibility
call artificialcomprs(beta)

! compute rhs
call getrhsstep2_incomprs(rhs)

! and multiply by inversed mass
do idof = 1,ndof
   invmmat = dmmat(idof)
   rhs(1,idof) = invmmat*deltp(idof)*(beta(idof)**2)*rhs(1,idof)
end do !idof
!******************************************************************************************
! update the solution: p at (n+1)
!
! p_(n+1) - p_(n) = beta^2*delta_t_ext*M^(-1)*RHS 
!******************************************************************************************
do idof = 1,ndof
   unkno(1,idof)= unkno(1,idof) + rhs(1,idof)
end do !idof
deallocate(beta)
end subroutine step2_incomprs
