subroutine timestep_incomprs
use femtypes
use fluidvariables
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
!    $Revision: 1.6 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!##########################################################################################
!
!  Compute Time-Step according to critical time step
!
!##########################################################################################
! for incompressible flow
! *** calculates the critical time step for all the elements and nodes
!******************************************************************************************
integer(I4B) :: i, ii
integer(I4B) :: ip,ie, nff
real(DP) :: alent, aloti, alotd, dm, uvall
real(DP) :: convts
real(DP), pointer :: uvel(:), vvel(:), uv(:)
!real(DP), pointer :: beta(:)
logical :: implc
!******************************************************************************************
! fixed time step
if(optts.le.-1) then
   deltp = dtfix
   delte = dtfix
   return
endif

! local time stepping for stokes equation
if (optstokes) then
! - From experience, delta t for stationary Stokes flow is not required to be calculated!!!.
! - use optts = -1 to set global time step = dtfix.
! - the stability of solution is obtained when delte is set to be around from 1e-7 to 1e-8 
   deltp = 1.0e6_DP
   do ie = 1, n
      nff = (ep(ie,1)+1)*(ep(ie,1)+2)/2
      alent = elh(ie)
! viscous limit only
      aloti = 0.5*(alent**2)/invre
      delte(ie) = csafm*aloti
      do i=1,nff
         ii = eg(ie,1)%d(i)
         deltp(ii) = min(deltp(ii), aloti) !, convts)
      end do ! i
   end do ! ie
   deltp = csafm*deltp
   return
end if

! local time stepping for incompressible flow

implc = (theta(2) .GT. 0)
deltp = 1.0e6_DP
delte = 1.0e6_DP
uvall = 0.0_DP
convts = 1.0e6_DP

do ie = 1, n

   nff = (ep(ie,1)+1)*(ep(ie,1)+2)/2
   alent = elh(ie)

! viscous limit only
   aloti = 0.5_DP*(alent**2)/invre
   
   if (optener) then
      ! diffusion limit   
      alotd = 0.5_DP*(alent**2)/thdiff 
      aloti = min(aloti,alotd)
   end if

!   allocate(uvel(nff), vvel(nff), uv(nff))
!   do i=1,nff
!      ii = eg(ie)%d(i)
!      uvel(i) = unkno(2,ii)
!      vvel(i) = unkno(3,ii)
!      uv(i) = dsqrt(uvel(i)**2 + vvel(i)**2)
!      uvall = max(uv(i),uvall)
!   end do
!   deallocate(uvel,vvel,uv)
!   if (uvall .GE. 1.0e-7_DP) then
!      convts = alent/uvall   
!   end if

   if (implc) then
! for semi-implicit scheme             
      do i=1,nff
         ii = eg(ie,1)%d(i)
         deltp(ii) = min(deltp(ii), aloti, convts)
         delte(ie) = min(delte(ie), deltp(ii))
      end do
   
   else
! for explicit scheme
      
      !call artificialcomprs(beta) ! need to modify for computing convection limit.
      do i=1,nff
         ii = eg(ie,1)%d(i)
         deltp(ii) = min(deltp(ii), aloti, convts)
         delte(ie) = min(delte(ie), deltp(ii))
      end do      
   end if
  
enddo !ie 

dm = 5.0e3_DP
do ip = 1,ndof
   dm = min(deltp(ip),dm)
enddo !ip

if (optts.eq.0) then
! for global time stepping due to stability condition
! minimum time step for all DOFs.
   deltp = csafm*dm
   delte = csafm*dm

! note:
! delte = delta t critical = h^2*Re/2
! P.Nithiarasu, O.C. Zeinkiewicz, "On stabilization of the CBS
! algorithm: Internal and external time steps
! Int. J. Numer. Meth. Engng. 2000,v.48,875,880
!   delte = csafm*dm
!   delte = 2.0_DP*csafm*dm

else if (optts == 1) then
! local time stepping
   deltp = csafm*dm
   delte = csafm*delte  
end if

end subroutine timestep_incomprs


!##########################################################################################
!                   COMPUTE ARTIFICIAL COMPRESSIBILITY
!##########################################################################################
subroutine artificialcomprs(beta)
use femtypes
use fluidvariables
use globalvariables
implicit none
real(DP), pointer :: beta(:)
!******************************************************************************************
!
!  Compute rhs of pressure equation
!
!******************************************************************************************
! local variables
!******************************************************************************************
integer(I4B) :: i
integer(I4B) :: nff
real(DP) :: sigma
!real(DP) :: uconv, udiff, utherm
real(DP) :: conv, diff, therm
!real(DP), pointer :: uvel(:), vvel(:)
!******************************************************************************************
allocate(beta(ndof))

! sigma is selected between 0.1 to 0.5
sigma = 0.5_DP

nff = 3
!allocate(uvel(nff),vvel(nff)) 

do i = 1,ndof

! convection limit
   conv = sqrt(unkno(2,i)**2+unkno(3,i)**2)  ! sqrt(unkn1(2,jj)**2+unkn1(3,jj)**2)
   !uconv = max(uconv,uvsum)     ! if beta is a single value for every node

! diffusion limit   
   diff  = 2*invre/alen(i)
   !udiff = max(udiff,diff)      ! if beta is a single value for every node

   if (optener) then
      therm = invre/(alen(i)*pr)
      !utherm = max(utherm,therm)   ! if beta is a single value for every node
      beta(i) = max(sigma,conv,diff,therm)
   else
      beta(i) = max(sigma,conv,diff)
   end if
end do !i

! if beta is only a single value for the whole domain
!beta = max(sigma,uconv,udiff,utherm)

end subroutine artificialcomprs