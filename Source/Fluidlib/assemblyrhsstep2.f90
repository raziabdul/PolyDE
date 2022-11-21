subroutine assemblyrhsstep2(rhs)
use femtypes
use feminterface
use fluidinterface
use fluidvariables
use globalvariables
implicit none
real(DP), pointer :: rhs(:)
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
!    $Revision: 1.6 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!*******************************************************************************
!
!  Element by element right-hand-side assembly in compact storage  
!
!remark: modified from assembly.f90 in solver
!*******************************************************************************
! local variables
!*******************************************************************************
integer(I4B) :: i,k,kn
integer(I4B) :: nff, vzk
real(DP), pointer :: bi(:)
!*******************************************************************************
! lower, upper, diag, ia, ja are defined in subr. pstiff seperately.

allocate(rhs(ndof))
rhs = 0._DP

!  loop over all elements (it should be possible to execute this loop in parallel)
do i=1,n

   call getrhsstep2_implicit(i,ep(i,1),bi,nff)

   do k=1,nff
!  row index k: local;  kn: global
      kn=abs(eg(i,1)%d(k))

      if (kn .eq. 0) cycle
      vzk=sign(1,eg(i,1)%d(k))

      rhs(kn) = rhs(kn)+bi(k)*vzk
   end do ! k
   deallocate(bi)
end do ! i

rhs(:) = rhs(:) - bglobal(:)
!rhs(:) = rhs(:)/deltp(:)
do i = 1,ndof
   rhs(i) = rhs(i)/(deltp(i)*invrho)
end do

end subroutine assemblyrhsstep2