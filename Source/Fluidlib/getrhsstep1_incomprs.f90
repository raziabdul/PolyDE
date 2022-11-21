subroutine getrhsstep1_incomprs(rhs)
use femtypes
use feminterface
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
!    $Revision: 1.7 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!*******************************************************************************
! get right-hand side step 1 for incompressible
!
! delte = local time-step
!*******************************************************************************
integer(I4B) :: i,j,k,ii,jj, m
integer(I4B) :: ie, ivar
integer(I4B) :: intorder, polyhi, polylo, polyorder, neldof
integer(I4B) :: iedge, nr, nl, branch, zrbtemp
integer(I4B) :: npkt1D, npkt, err, errcode
integer(I4B), dimension(3), parameter :: inach=(/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg=(/3,1,2/)
real(DP) :: uadv, vadv, advtemp
real(DP) :: usum, vsum, umean, vmean
real(DP) :: sumu, sumv
real(DP) :: aw, area, delteby2, awre
real(DP) :: lam(3), lw, length2
real(DP) :: dudx, dudy, dvdx, dvdy
real(DP) :: qq, pp
real(DP) :: const2, const3, tempterm
real(DP) :: fmagterm, pressterm, magterm
real(DP) :: tsum, psum, hsum, tmean, pmean, hmean
real(DP), allocatable :: xtab(:), weight(:), lambda(:,:)
real(DP), pointer :: uvel(:), vvel(:), tt(:)
real(DP), pointer :: hmag(:), press(:), fmag(:)
real(DP), pointer :: xsi(:), gxsi(:,:)
real(DP), pointer :: rhsstep1(:,:), diffside(:,:), elvel(:,:)
real(DP), pointer :: bodyforce(:,:)
!*******************************************************************************
if (optener) then
   if (optdim) then   
      if (optnat) then
      ! find terms of const2 
         const2 = gravity*thexpn
         const3 = -gravity*thexpn*tfree
      else
         const2 = 0.0_DP
         const3 = 0.0_DP
      end if
   else
      if (optnat) then
         const2 = pr*ra
      else 
         const2 = ri
      end if 
   end if
end if

do ie = 1,n

   polyorder = ep(ie,1)
   polyhi = polyorder
   polylo = 1
   neldof = (polyorder+1)*(polyorder+2)/2

   allocate(rhsstep1(2,neldof), diffside(2,neldof), elvel(2,neldof), bodyforce(2,neldof))
   allocate(uvel(neldof),vvel(neldof))
   bodyforce = 0._DP
   
   if (optener) then
      allocate(tt(neldof))
      do i =1,neldof
         ii = eg(ie,1)%d(i)
         tt(i) = unkno(4,ii)
      end do !i
   end if

   area = areaf(ie)
   delteby2 = delte(ie)*0.5_DP
! initialization
   usum = 0._DP
   vsum = 0._DP
   do i = 1,neldof
      ii = eg(ie,1)%d(i)
      uvel(i) = unkno(2,ii)
      vvel(i) = unkno(3,ii)
      usum = usum + uvel(i)
      vsum = vsum + vvel(i)
      
      do ivar = 1,2
         ! define local (elemental) velcities
         elvel(ivar,i) = unkno(ivar+1,ii)
      end do !ivar
   end do !i

   umean = usum/real(neldof,DP)
   vmean = vsum/real(neldof,DP) 

! compute related parameter for the external magnetic force
   if (optmagn) then
      psum = 0._DP
      tsum = 0._DP
      hsum = 0._DP
      allocate(tt(neldof),press(neldof),hmag(neldof),fmag(neldof))
      do i =1,neldof
         ii = eg(ie,1)%d(i)
         press(i) = unkno(1,ii)
!         tt(i) = unkno(4,ii)*delt+tfree
         tt(i) = unkno(4,ii)
         hmag(i) = magn(ii)
         psum = psum + press(i)
         tsum = tsum + tt(i)
         hsum = hsum + hmag(i)
! local definition for paramagnetic force
!         fmag(i) = magconst*press(i)*(hmag(i)**2)/(tt(i)**2)
      end do !i
      pmean = psum/real(neldof,DP)
      tmean = tsum/real(neldof,DP)
      hmean = hsum/real(neldof,DP)
      do i =1,neldof
! with T square (1/T**2)
! first implementation
! all other variables are assumed to use an average value for each element
         fmag(i) = magconst*pmean*hmean*hmag(i)/(tmean**2)
!         fmag(i) = magconst*pmean*(hmag(i)**2)/(tmean**2)
! modified to have only temperature dependence of (1/T)
!         fmag(i) = magconst*hmean*hmag(i)/tmean
!         fmag(i) = magconst*(hmag(i)**2)/tmean
      end do !i
      deallocate(tt,press,hmag)
   end if

   allocate(xsi(neldof), gxsi(neldof,2) )
   intorder=2*polyorder
   call get2Dintegpoints(intorder, npkt, weight, lambda, err)

! initialization
   uadv = 0._DP
   vadv = 0._DP
   rhsstep1 = 0._DP

   do k = 1,npkt
      call shapefunction(lambda(:,k),xn(e(:,ie)),yn(e(:,ie)),     &
                         polylo,polyhi,neldof,.true.,xsi,gxsi,errcode)
      aw = area*weight(k)
      awre = aw*invre  
      
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
         uadv = advtemp*uvel(j)
         vadv = advtemp*vvel(j)

         ! diffusive terms (invre = 1/Re)
         dudx = awre*uvel(j)*gxsi(j,1)
         dudy = awre*uvel(j)*gxsi(j,2)
         dvdx = awre*vvel(j)*gxsi(j,1)
         dvdy = awre*vvel(j)*gxsi(j,2)
         
         do i = 1,neldof

            rhsstep1(1,i) = rhsstep1(1,i) - xsi(i)*uadv &
                                          - delteby2*uadv*(gxsi(i,1)*umean + gxsi(i,2)*vmean) &
                                          - gxsi(i,1)*dudx - gxsi(i,2)*dudy
                                          !- delteby2*uadv*(gxsi(i,1)*uvel(i)+gxsi(i,2)*vvel(i)) &
            rhsstep1(2,i) = rhsstep1(2,i) - xsi(i)*vadv &
                                          - delteby2*vadv*(gxsi(i,1)*umean + gxsi(i,2)*vmean) &
                                          - gxsi(i,1)*dvdx - gxsi(i,2)*dvdy
                                          !- delteby2*vadv*(gxsi(i,1)*uvel(i)+gxsi(i,2)*vvel(i)) &
         end do ! i
! this is just an idea to store all b.c. to an array in which it can be re-used and
! might be suitable for such a multistep method in the fluid problems
! must define boundary side check first
!         lwre = weight(k)*nside(3,ib)/2.0_DP  
!         normal = nside(1,ib)*gxsi(j,1) + nside(2,ib)*gxsi(j,2)
!         diffside(1,j) = diffside(1,j) + lwre*normal*uvel(j)
!         diffside(2,j) = diffside(2,j) + lwre*normal*vvel(j)
      end do ! j

      ! external force
      if (optener) then
         if (optdim) then
            if (optnat) then
               do j = 1,neldof
                  tempterm = const2*aw*xsi(j)
                  do i = 1,neldof
                     bodyforce(2,i) = bodyforce(2,i) + tempterm*xsi(i)*tt(i) + const3
                  end do ! i
               end do ! j
            end if
         else
            do j = 1,neldof
               tempterm = const2*aw*xsi(j)
               do i = 1,neldof
                  bodyforce(2,i) = bodyforce(2,i) + tempterm*xsi(i)*tt(i)
               end do ! i
            end do ! j
         end if      
      end if

      if (optmagn) then
         do j = 1,neldof
            fmagterm = aw*xsi(j)
!            pressterm = xsi(j)*press(j)
!            tempterm = xsi(j)*tt(j)
!            tempterm = 0.0_DP
!            do k = 1,neldof
!               tempterm = tempterm + xsi(j)*tt(j)*xsi(k)*tt(k)
!            end do 
!            magterm = xsi(j)*hmag(j)
!           do i =1,neldof
!            bodyforce(1,i)=bodyforce(1,i) + fmagterm*pressterm*gxsi(i,1)*hmag(i)*magterm/tempterm
!            bodyforce(2,i)=bodyforce(2,i) + fmagterm*pressterm*gxsi(i,2)*hmag(i)*magterm/tempterm

! first implementation
            bodyforce(1,i) = bodyforce(1,i) + fmagterm*gxsi(i,1)*fmag(i)
            bodyforce(2,i) = bodyforce(2,i) + fmagterm*gxsi(i,2)*fmag(i)
!            end do !i
         end do ! j
      end if
   end do !k
   deallocate(weight,lambda, gxsi)   
   if (associated(tt)) deallocate(tt)
   if (associated(press)) deallocate(press)
   if (associated(hmag)) deallocate(hmag)
!*******************************************************************************
! contribution of boundary side
!*******************************************************************************
   diffside = 0._DP
   if (genbc(ie)) then

      call get1Dintegpoints(intorder, xtab, weight, npkt1D, errcode)

      do iedge = 1,3
         if (branchef(iedge,ie) .EQ. 0) cycle

         nr = inach(iedge)
         nl = ivorg(iedge)
         length2 = lengthby2(iedge,ie)
         lam(iedge) = 0._DP

         do ivar = 2,3
            if (.NOT. generalbc(ie)%bcedge(iedge)%varmn(ivar)) cycle
            branch = branchef(iedge,ie)

            zrbtemp = zrbmn(branch,ivar)

            if (zrbtemp .lt. 200 .and. zrbtemp .ge. 300) cycle

            pp = real(alrbmn(branch,ivar),DP)
            qq = real(btrbmn(branch,ivar),DP)

            do k=1, npkt1D
!  get shape function at location of integration points
               lam(nr)=(xtab(k)+1._DP)/2._DP
               lam(nl)=1._DP-lam(nr)
               call shapefunction(lam,xn(e(:,ie)),yn(e(:,ie)),polylo,polyhi, &
                                  neldof,.false.,xsi,gxsi,errcode)

               lw = length2*weight(k)*invre

               do j=1,neldof
                  do i=1,neldof
                     diffside(ivar-1,i) = diffside(ivar-1,i) + lw*qq*xsi(j)*xsi(i)*elvel(ivar-1,i)
                  end do ! i
                  diffside(ivar-1,j) = diffside(ivar-1,j) + lw*pp*xsi(j)
               end do ! j
            end do ! k
         end do !ivar
      end do ! iedge
      deallocate(xtab,weight)
   end if
   deallocate(xsi)
   deallocate(uvel,vvel,elvel)
   if (associated(tt)) deallocate(tt)
!*******************************************************************************
! sum of all rhs terms to create rhs
!*******************************************************************************
   do j = 1, neldof
      jj = eg(ie,1)%d(j)
      rhs(2,jj) = rhs(2,jj) + rhsstep1(1,j) + diffside(1,j) + bodyforce(1,j)
      rhs(3,jj) = rhs(3,jj) + rhsstep1(2,j) + diffside(2,j) + bodyforce(2,j)
   end do ! j
   deallocate(rhsstep1,diffside,bodyforce)

end do ! ie

end subroutine getrhsstep1_incomprs

