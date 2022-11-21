subroutine nusselt
use femtypes
use feminterface, only: reallocate, get1Dintegpoints, shapefunction, get2Dintegpoints
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
!    $Revision: 1.5 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!*********************************************************************************
!
!     Rate of Heat Transfer (Nusselt Number (Nu))
!  
!*********************************************************************************
! Non-dimensional:
!  
!   Nu  =  - d(T)/dn           - at the surface where constant wall temperature
!
!   Nu  =  1/(Tw - Tf)         - at the surface where heat flux is assigned
!
! Dimensional:
!
!   Nu  =  h*L/k = L*(-d(T)/dn)/(Tw-Tf) - at the surface where constant wall temperature
!
!   Nu  =  q*L/(k*(Tw - Tf))   - at the surface where heat flux is assigned
!
!     h  - heat transfer coefficient
!     L  - characteristic length
!     k  - thermal conductivity of fluid
!     Tw - wall temperature
!     Tf - fluid temperature
!     q  - heat flux
!*********************************************************************************
! local variables
!*********************************************************************************
integer(I4B) :: bredge
integer(I4B) :: i, ib, ie, iedge, ii, inu, j, k
integer(I4B) :: intorder
integer(I4B) :: err, npkt       ! area integral
integer(I4B) :: errcode, npkt1D, nr, nl     ! line integral
integer(I4B) :: nunum, nff
integer(I4B) :: polylo, polyhi, polyorder
integer(I4B), pointer :: nubranch(:)
integer(I4B), dimension(3), parameter :: inach=(/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg=(/3,1,2/)
real(DP) :: eleng, lengtot
real(DP) :: dtdn, dtdx, dtdy
real(DP) :: tvalue
real(DP) :: lam(3)
real(DP), allocatable :: xtab(:), weight(:)
real(DP), pointer :: tt(:), nulocal(:,:)
real(DP), pointer :: xsi(:), gxsi(:,:)
real(DP), pointer :: lambda(:,:)
logical :: tdbc
!*********************************************************************************
allocate(nubranch(nbound),nulocal(2,nbound))
nunum = 0
nulocal = 0.0_DP

do ib = 1,nbound
   ie  = bcside(1,ib)
   iedge = bcside(2,ib)
   tdbc = dirichletbc(ie)%bcedge(iedge)%varmn(4)
   if (.NOT. tdbc) cycle
   bredge = branchef(iedge,ie)
   tvalue = alrbmn(bredge,4)
! to find the heat transfer rate at warm surface
   if (optdim) then
      if (tvalue .LE. tfree) cycle
   else
      if (tvalue .LT. 1.0_DP) cycle
   end if
   nunum = nunum + 1
   nubranch(nunum) = ib

   polyorder = ep(ie,1)
   polyhi = polyorder
   polylo = 1
   nff=(polyorder+1)*(polyorder+2)/2
   allocate(tt(nff))

   do i = 1,nff
      ii = eg(ie,1)%d(i)
      tt(i) = unkno(4,ii)
   end do !i

   allocate(xsi(nff), gxsi(nff,2) )
   dtdx = 0._DP
   dtdy = 0._DP

   intorder=2*polyorder
!   call get2Dintegpoints(intorder, npkt, weight, lambda, err)
!********
   call get1Dintegpoints(intorder, xtab, weight, npkt1D, errcode)
!
   nr=inach(iedge)
   nl=ivorg(iedge)
   lam(iedge)=0._DP
!  now integrate with the boundary conditions of the input branch: branche(iedge)
   do k=1, npkt1D
!  get shape function at location of integration points
      lam(nr)=(xtab(k)+1._DP)/2._DP
      lam(nl)=1._DP-lam(nr)
      call shapefunction(lam,xn(e(:,ie)),yn(e(:,ie)),         &
     &        polylo,polyhi,nff,.true.,xsi,gxsi,errcode)

!*******   
!   do k=1,npkt
!   !  get shape function and gradients at location of integration points
!      call shapefunction(lambda(:,k),xn(e(:,ie)),yn(e(:,ie)),     &
!     &    polylo,polyhi,nff,.true.,xsi,gxsi,errcode)

      do j=1,nff
         dtdx = dtdx + weight(k)*gxsi(j,1)*tt(j)
         dtdy = dtdy + weight(k)*gxsi(j,2)*tt(j)
      end do ! j

   end do !k
   deallocate(xtab, weight)
   deallocate(xsi, gxsi)
   deallocate(tt)
   dtdn = dtdx*nside(1,ib) + dtdy*nside(2,ib)

   do i = 1,2
      nulocal(i,nunum) = dtdn
   end do
   
end do !ib
nubranch => reallocate(nubranch,nunum)

nuaver = 0.0_DP
lengtot = 0.0_DP

do inu = 1,nunum
   ib = nubranch(inu)
   eleng = nside(3,ib)
   nuaver = nuaver + 0.5_DP*(nulocal(1,inu)+nulocal(2,inu))*eleng
   lengtot = lengtot + eleng
end do !inu
nuaver = nuaver/lengtot

if (optdim) then
   nuaver = nuaver*thcond/tmpdiff
   write(*,100) 'average heat transfer coefficient (W/m2K) = ', nuaver
   100 format(1x,a45,es17.10)
else
   write(*,101) 'average Nusselt number = ', nuaver
   101 format(1x,a25,es17.10)
end if

deallocate(nubranch,nulocal)
end subroutine nusselt
