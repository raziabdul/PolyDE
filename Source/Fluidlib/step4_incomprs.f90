subroutine step4_incomprs(rhs)
use femtypes
use fluidvariables
use fluidinterface
use globalvariables
implicit none
real(DP), pointer :: rhs(:,:)
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
!################################################################
!
!     Temperature Calculation
!
!################################################################
! local variables
!****************************************************************
integer(I4B) :: i
!****************************************************************
rhs(4,:) = 0.0_DP

call getrhsstep4_incomprs(rhs)

do i = 1,ndof
   rhs(4,i) = rhs(4,i)*dmmat(i)*deltp(i)
   unkno(4,i) = unkno(4,i) + rhs(4,i)
end do

end subroutine step4_incomprs