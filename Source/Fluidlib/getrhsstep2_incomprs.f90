subroutine getrhsstep2_incomprs(rhs)
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
!    $Revision: 1.5 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!***********************************************************************************
!
!  Compute rhs of pressure equation
!
! remark : this routine is similar to getrhsstep2.f90 for compressible flow
! it is created only to test explicit form of incompressible flow.
!***********************************************************************************
! local variables
!***********************************************************************************
integer(I4B) :: i,j,k,jj
integer(I4B) :: elem, ivar
integer(I4B) :: intorder, polyhi, polylo, polyorder, nff
integer(I4B) :: iedge, nr, nl, branch, zrbtemp
integer(I4B) :: npkt1D, npkt, err, errcode
integer(I4B), dimension(3), parameter :: inach=(/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg=(/3,1,2/)
real(DP) :: aw
real(DP) :: pp, qq, lam(3), lw, length2
real(DP) :: velterm, presdiffx, presdiffy
real(DP), allocatable :: xtab(:), weight(:), lambda(:,:)
real(DP), pointer :: xsi(:), gxsi(:,:)
real(DP), pointer :: rhsvel(:), rhspres(:), diffside(:)
real(DP), pointer :: elpres(:), elvel(:,:)
!***********************************************************************************
polyorder = 1
polyhi=polyorder
polylo=1
!  size of the element matrix
nff=(polyorder+1)*(polyorder+2)/2

allocate( rhsvel(nff), rhspres(nff), diffside(nff) )
allocate( elpres(nff), elvel(2,nff) )

do elem=1,n

   do j = 1,nff
      jj = eg(elem,1)%d(j)
      do ivar = 2,3
         elvel(ivar-1,j) = theta(1)*unkno(ivar,jj)
      end do ! ivar
      elpres(j) = unkno(1,jj)
   end do ! j

   allocate(xsi(nff), gxsi(nff,2) )
   intorder=2*polyorder
   !  fetch numerical integration points (Gauss points)
   call get2Dintegpoints(intorder, npkt, weight, lambda, err)

   rhsvel = 0._DP
   rhspres = 0._DP

   do k=1,npkt
   !  get shape function and gradients at location of integration points
      call shapefunction(lambda(:,k),xn(e(:,elem)),yn(e(:,elem)),polylo,polyhi,nff,.true.,&
                         xsi,gxsi,errcode)
      aw=areaf(elem)*weight(k)

   !  for all shape functions
      do j=1,nff     
         velterm   = aw*xsi(j)
         presdiffx = aw*gxsi(j,1)*elpres(j)
         presdiffy = aw*gxsi(j,2)*elpres(j)
         
         do i=1,nff
            rhsvel(i) =  rhsvel(i) - velterm*( gxsi(i,1)*elvel(1,i) + gxsi(i,2)*elvel(2,i) )
            rhspres(i) = rhspres(i) - ( gxsi(i,1)*presdiffx + gxsi(i,2)*presdiffy )
!            rhspres(i) = rhspres(i) - aw*(gxsi(j,1)*gxsi(i,1) + gxsi(j,2)*gxsi(i,2))*elpres(i)
         end do !i
      end do !j      
   end do !k
   deallocate(weight,lambda)
!***********************************************************************************************
! boundary integral of pressure term
!***********************************************************************************************
!   allocate(diffside(nff)) move to the top of program for equal order
   diffside = 0._DP
   if (genbc(elem)) then

      call get1Dintegpoints(intorder, xtab, weight, npkt1D, errcode)

      do iedge=1,3       
         if (branchef(iedge,elem) .EQ. 0) cycle
        
         nr=inach(iedge)
         nl=ivorg(iedge)
         length2 = lengthby2(iedge,elem)

         lam(iedge)=0._DP
!  now integrate with the boundary conditions of the input branch
         
         if (.NOT. generalbc(elem)%bcedge(iedge)%varmn(1)) cycle
         branch = branchef(iedge,elem)

         zrbtemp = zrbmn(branch,1)

         if (zrbtemp .lt. 200 .and. zrbtemp .ge. 300) cycle 
         
         pp = real(alrbmn(branch,1),DP)
         qq = real(btrbmn(branch,1),DP)
         
         do k=1,npkt1D
!  get shape function at location of integration points
            lam(nr)=(xtab(k)+1._DP)/2._DP
            lam(nl)=1._DP-lam(nr)
            call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),polylo,polyhi, &
                            nff,.false.,xsi,gxsi,errcode)
            lw = length2*weight(k)

            do j=1,nff          
               do i=1,nff          
                  diffside(i) = diffside(i) + lw*qq*xsi(i)*xsi(j)*elpres(i)
               end do ! i
               diffside(j) = diffside(j) + lw*pp*xsi(j)
            end do ! j
         end do ! k
      end do ! iedge
      deallocate(xtab,weight)
   end if
   deallocate(xsi,gxsi)
!**************************************************************************************
! sum of all rhs terms to create rhs
!**************************************************************************************      
   do j = 1,nff
      jj = eg(elem,1)%d(j)
      rhs(1,jj) = rhs(1,jj) + rhspres(j) + rhsvel(j) + diffside(j)
   end do ! j  
end do !elem


end subroutine getrhsstep2_incomprs
!###################################################################################

!###################################################################################
