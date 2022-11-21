subroutine step2_implicit(lower,diag,ia,ja)
use femtypes
use feminterface
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer(I4B), pointer :: ia(:), ja(:)
real(DP), pointer :: lower(:), diag(:)
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
!    $Date: 2015/04/01 10:53:05 $
!    $Author: juryanatzki $
!
!*****************************************************************************
!
!  Step 2: Pressure calculation by implicit solver
!
!*****************************************************************************
real(DP), pointer :: rhs(:)
real(DP), pointer :: lower1(:), diag1(:)
real(DP) :: eps, epsgl, resgl
character(len=25) :: solver
logical :: symmetric
!*****************************************************************************
allocate(lower1(1:ia(ndof)), diag1(ndof))

! assembly element matrix to the global matrix
call assemblyrhsstep2(rhs)

call getsetting('LINSOLVER_ERROR',eps)
call getsetting('LINSOLVERTYPE',solver)

call getsetting('PHYSICS_MODE',physics)
if (physics .eq. 'FLUID') symmetric =.true.

lower1 = lower
diag1  = diag

select case(solver)

case ('SSORCG')
   call solve2(lower1,lower1,diag1,rhs,unkno(1,:),ndof,eps,ia,ja,epsgl,resgl,symmetric)

case default
   print*,'Step 2 implicit incompressible solver uses only SSORCG'
   stop

end select
deallocate(lower1,diag1,rhs)

end subroutine step2_implicit