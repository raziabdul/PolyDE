subroutine checkconverge(eps,err,converge)
use femtypes
use feminterface
use fluidvariables
use fluidinterface, only:
use globalvariables
implicit none
real(DP), intent(in)  :: eps
real(DP), pointer :: err(:)
logical, intent(out)  :: converge
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
!    $Revision: 1.3 $
!    $Date: 2006/07/07 22:34:04 $
!    $Author: r_abdul $
!
!*******************************************************************************
!
!  Check Convergence for Steady-state Solution
!
!*******************************************************************************
! local variables
!*******************************************************************************
integer(I4B) :: i, ivar
real(DP) :: diff(nvar)
logical :: crit(nvar)
!*******************************************************************************
allocate(err(nvar))
err = 0._DP

do ivar = 1,nvar
   do i = 1,ndof
      diff(ivar) = unkno(ivar,i)-unkn1(ivar,i)
      err(ivar) = err(ivar) + (diff(ivar)*diff(ivar))
   end do
   err(ivar) = dsqrt(err(ivar))
   crit(ivar) = (err(ivar) .LE. eps)
end do

! check convergence for all variables to be less than "eps".
if (optcomprs .OR. optener) then
   converge = (crit(1) .AND. crit(2) .AND. crit(3) .AND. crit(4))
else
   converge = (crit(1) .AND. crit(2) .AND. crit(3))
end if

end subroutine checkconverge