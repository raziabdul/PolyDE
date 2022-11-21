subroutine solver_stokes
use femtypes
use feminterface, only: getsetting, solve2, umfsolver
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
!    $Revision: 1.5 $
!    $Date: 2015/04/01 10:52:02 $
!    $Author: juryanatzki $
!
!*****************************************************************************
!
!     Solver for Stokes flow with equal-order 
!
!*****************************************************************************
! local variables
!*****************************************************************************
integer(I4B) :: ntdof, nudof, npdof, ndof2
integer(I4B), pointer :: ia(:), ja(:)
real(DP) :: epsgl, resgl
real(DP) :: eps
real(DP), pointer :: lower(:), upper(:), diag(:), rhs(:), acsr(:)
real(DP), pointer :: sol(:)
character(len=25) :: solver, csroption
logical :: csr, jacobi, matvar
!*****************************************************************************
! Get options from FEMsettings.txt
call getsetting('PHYSICS_MODE',physics)
call getsetting('CSRFORMAT',csroption)

! Option for Compact Storage
if (csroption .eq. 'CSR') then
   csr=.true.
! Default is the lower CSR
else
   csr=.false.   
end if

! select which method will be used to solve Stokes flow
! physics == STOKES -> mixed method
! physics == FLUID -> equal-order method with penalty values

! mixed method
if (physics .eq.'STOKES') then
   optmixed = .true.

   call preassemb_mixed(nudof,npdof)

! number of variables, 3 = u,v and p
   nvar = 3
! read multivariable boundary conditions
   call readnetin_mn
   ndof2 = nudof*2
   ntdof = nudof*2 + npdof

! csr is set to be .true. first
!   csr = .true.
   jacobi=.false.
   matvar=.false.
! if this call is used, there is no need to use the lines after this call.
!   call assembly_stokes(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr,matvar,ntdof)

!*********** use setupia ****************
   if (csr) then
      allocate(ia(ntdof+1))
   else 
      allocate(ia(ntdof))
   end if
   call setupia_temp(ia,csr,ntdof)
   if (csr) then
      allocate(ja(1:ia(ntdof+1)-1))
      allocate(acsr(1:ia(ntdof+1)-1),rhs(ntdof))
   else
      allocate(ja(1:ia(ntdof)))
      allocate(lower(1:ia(ntdof)),upper(1:ia(ntdof)),diag(ntdof),rhs(ntdof))
   end if

! assemble element matrices 
   call assembly_nostep(lower,upper,diag,acsr,rhs,ia,ja,csr,npdof)

!****************************************

! equal-order method
else if (physics .eq. 'FLUID') then
   optmixed = .false.

   ndof2 = ndof*2
   ntdof = ndof*3

   allocate(deltp(ndof),delte(n))
   call timestep_incomprs

   ! unsymmetric matrix
   if (optstokes) then
      csr = .true.
   end if ! optstokes

   if (csr) then
      allocate(ia(ntdof+1))
   else 
      allocate(ia(ntdof))
   end if
   ! setup ia(:)   
   call setupia_stokes(ia,csr)

   if (csr) then
      allocate(ja(1:ia(ntdof+1)-1))
      allocate(acsr(1:ia(ntdof+1)-1),rhs(ntdof))
   else
      allocate(ja(1:ia(ntdof)))
      allocate(lower(1:ia(ntdof)),upper(1:ia(ntdof)),diag(ntdof),rhs(ntdof))
   end if            

   ! assemble element matrices 
   call assembly_nostep(lower,upper,diag,acsr,rhs,ia,ja,csr,npdof)   

end if ! if physics

! solution process
allocate(sol(ntdof))
sol = 0._DP

call getsetting('LINSOLVERTYPE',solver)

if (csr) then
   select case(solver)

   case default
      print*,' no such solver: ',solver,' use SSORCG'
   end select
else
   call getsetting('LINSOLVER_ERROR',eps)
   select case(solver)
   case ('SSORCG')
      call solve2(lower,upper,diag,rhs,sol,ntdof,eps,ia,ja,epsgl,resgl,.false.)

   case default
      print*,'uses only SSORCG'
   end select
   deallocate( ia, ja, lower, upper, diag, rhs)
end if

if (.not. associated(unkno)) then
   allocate(unkno(4,ndof))
end if
! transform to unkno for lout.f90
unkno(2,1:ndof) = sol(1:ndof)
unkno(3,1:ndof) = sol(ndof+1:ndof2)
unkno(1,1:npdof) = sol(ndof2+1:ntdof)

deallocate(sol)

! For mixed method, DOFs of velocity and pressure are not identical 
if (optmixed) then
   call lout_mixed(.true., .true., epsgl, npdof)
end if ! optmixed

end subroutine solver_stokes