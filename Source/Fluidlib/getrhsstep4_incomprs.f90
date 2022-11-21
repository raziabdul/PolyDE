subroutine getrhsstep4_incomprs(rhs)
use femtypes
use globalvariables
use fluidvariables
use feminterface
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
!    $Revision: 1.6 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!***********************************************************************************
!
!  Determine RHS of step 4 (energy equation) - incompressible
!
!***********************************************************************************
real(DP), pointer :: rhs(:,:) ! intent(out)
!***********************************************************************************
! local variables
!***********************************************************************************
integer(I4B) :: i, ii, ie, j, jj, k, m
integer(I4B) :: iedge, nr, nl, npkt1D
integer(I4B) :: npkt, polylo, polyhi
integer(I4B) :: intorder, polyorder, neldof
integer(I4B) :: err, errcode, branch, zrbtemp
integer(I4B), dimension(3), parameter :: inach=(/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg=(/3,1,2/)
real(DP) :: aw, area, advtemp, const, awcon, delteby2
real(DP) :: sumu, sumv, usum, vsum, umean, vmean
real(DP) :: dtdx, dtdy, tadv
real(DP) :: length2, lw, lwcon, lam(3), pp, qq
real(DP), allocatable :: xtab(:), weight(:), lambda(:,:)
real(DP), pointer :: uvel(:), vvel(:), tt(:)
real(DP), pointer :: advterm(:), diffterm(:), stabterm(:), diffside(:)
real(DP), pointer :: xsi(:), gxsi(:,:)
!***********************************************************************************
const = thdiff
! non-dimensional is set in presetting.f90
!const = invre/pr

do ie = 1,n

   polyorder = ep(ie,1)
   polyhi = polyorder
   polylo = 1
   neldof = (polyorder+1)*(polyorder+2)/2

   allocate(uvel(neldof), vvel(neldof), tt(neldof))
   allocate(advterm(neldof), diffterm(neldof), stabterm(neldof), diffside(neldof))

   area = areaf(ie)
   delteby2 = delte(ie)*0.5_DP

! initialization
   usum = 0._DP
   vsum = 0._DP
   do i = 1,neldof
      ii = eg(ie,1)%d(i)
      uvel(i) = unkno(2,ii)
      vvel(i) = unkno(3,ii)
      tt(i)   = unkno(4,ii)
      usum = usum + uvel(i)
      vsum = vsum + vvel(i)
!      do ivar = 1,2
!         ! define local (elemental) velcities
!         elvel(ivar,i) = unkno(ivar+1,ii)
!      end do !ivar
   end do !i

   umean = usum/real(neldof,DP)
   vmean = vsum/real(neldof,DP)

   allocate(xsi(neldof), gxsi(neldof,2) )
   intorder = 2*polyorder
   call get2Dintegpoints(intorder, npkt, weight, lambda, err)

   advterm  = 0._DP
   stabterm = 0._DP
   diffterm = 0._DP

   do k = 1,npkt

      call shapefunction(lambda(:,k),xn(e(:,ie)),yn(e(:,ie)),     &
                         polylo,polyhi,neldof,.true.,xsi,gxsi,errcode)
      aw = area*weight(k)
      awcon = aw*const

      ! advective velocity
      sumu = 0._DP
      sumv = 0._DP
      do m = 1,neldof
         sumu = sumu + xsi(m)*uvel(m)
         sumv = sumv + xsi(m)*vvel(m)
      end do

      do j = 1,neldof

         ! advective terms
         advtemp = aw*(sumu*gxsi(j,1) + sumv*gxsi(j,2))
         tadv = advtemp*tt(j)
!         adv = aw*(umean*gxsi(j,1) + vmean*gxsi(j,2))
         dtdx = awcon*gxsi(j,1)*tt(j)
         dtdy = awcon*gxsi(j,2)*tt(j)

         do i = 1,neldof        
            advterm(i) = advterm(i) - xsi(i)*tadv
            stabterm(i) = stabterm(i) - delteby2*tadv*(gxsi(i,1)*umean + gxsi(i,2)*vmean)
            diffterm(i) = diffterm(i) - gxsi(i,1)*dtdx - gxsi(i,2)*dtdy
         end do ! i
      end do ! j

   end do !k
   deallocate(weight,lambda, gxsi)
!***********************************************************************************
!  boundary side contribution
!***********************************************************************************
   diffside = 0._DP
   if (genbc(ie)) then

      call get1Dintegpoints(intorder, xtab, weight, npkt1D, errcode)

      do iedge=1,3

         if (branchef(iedge,ie) .eq. 0) cycle
         if (.NOT. generalbc(ie)%bcedge(iedge)%varmn(4)) cycle
         branch  = branchef(iedge,ie)
         zrbtemp = zrbmn(branch,4)
         if (zrbtemp.lt.200 .or. zrbtemp.ge.300) cycle

         pp = real(alrbmn(branch,4),DP)
         qq = real(btrbmn(branch,4),DP)

         nr=inach(iedge)
         nl=ivorg(iedge)
         length2 = lengthby2(iedge,ie)

         lam(iedge)=0._DP
!  now integrate with the boundary conditions of the input branch: branche(iedge)
         do k=1, npkt1D
!  get shape function at location of integration points
            lam(nr)=(xtab(k)+1._DP)/2._DP
            lam(nl)=1._DP-lam(nr)
            call shapefunction(lam,xn(e(:,ie)),yn(e(:,ie)),polylo,polyhi, &
                               neldof,.false.,xsi,gxsi,errcode)
            lw = length2*weight(k)
            lwcon = lw*const

            do j=1,neldof
               do i=1,neldof
                  diffside(i) = diffside(i) + lw*qq*xsi(i)*xsi(j)*tt(i) &
                                            + lw*pp*xsi(i)
               end do ! i
            end do ! j
         end do ! k
      end do ! iedge
      deallocate(xtab, weight)
   end if
   deallocate(uvel,vvel,tt)
   deallocate(xsi)
!***********************************************************************************
! combine all terms
!***********************************************************************************
   do j = 1, neldof
      jj = eg(ie,1)%d(j)
      rhs(4,jj) = rhs(4,jj) + advterm(j) + stabterm(j) + diffterm(j) + diffside(j)       
   end do ! j
   deallocate(advterm, diffterm, stabterm, diffside)
!***********************************************************************************
end do ! ie

end subroutine getrhsstep4_incomprs