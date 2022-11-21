subroutine geteleminfo
use femtypes
use feminterface
use globalvariables
use fluidvariables
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
!*******************************************************************************
!  Calculate element area and dof-wise 
!*******************************************************************************
!
! alen = vector containing element length
!
!*******************************************************************************
! local variables
!*******************************************************************************
integer(I4B) :: i, ie, ii, ied1, ied2, iedge, ivar, ivertex, ipbc
integer(I4B) :: j, k, pp
integer(I4B) :: nr, nl, zweigright, zweigleft
integer(I4B) :: zweigr, zweigl, kziver 
integer(I4B) :: brvmn, bredge, zrbtemp, zrbvertex
integer(I4B) :: kziverright,kziverleft
integer(I4B) :: polyorder, nff, polyhi, polylo
integer(I4B) :: branchv(3)
integer(I4B) :: maxbcsize, neighbor
integer(I4B), pointer :: branchemn(:,:), branchvmn(:,:)
integer(I4B), dimension(3), parameter :: inach = (/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg = (/3,1,2/)
real(DP) :: dxl, dyl, xed, yed, elen
real(DP) :: shelleng, bleng
real(DP) :: vlen(3)
real(DP) :: opplen, delx, dely
!*******************************************************************************
allocate(alen(ndof),areaf(n))
! initialize alen
alen = 1.0e6_DP

! temporary variables
allocate(elh(n))     ! shortest element length of each element
elh = 1.0e6_DP   
maxbcsize = min(ndof,n)
allocate(bcside(5,maxbcsize))

nbound = 0

allocate(generalbc(n),dirichletbc(n))
allocate(genbc(n),edbc(n))
allocate(branchef(3,n), lengthby2(3,n))
! local pointers
allocate(branchvmn(3,nvar),branchemn(3,nvar))
genbc = .false.
edbc  = .false.
lengthby2 = 0._DP


do ie = 1, n

! area of the element
   areaf(ie) = ( (yn(e(3,ie))-yn(e(1,ie)))*(xn(e(2,ie))-xn(e(3,ie)))-                             &
                 (yn(e(2,ie))-yn(e(3,ie)))*(xn(e(3,ie))-xn(e(1,ie))) )/2._DP 

   polyorder = ep(ie,1)
   polyhi=polyorder
   polylo=1
   nff = (polyorder+1)*(polyorder+2)/2
   
   pp = polylo
   shelleng = 1.0e6_DP
   if (polyorder .EQ. 1) then
      do i = 1, nff
         ii = eg(ie,1)%d(i)
         ! nodal element-length = distance from vertex to mid of the opposite edge 
         ied1 = inach(i)
         ied2 = ivorg(i)      
         !xed = 0.5_DP*( xn(e(ied1,ie)) + xn(e(ied2,ie)) )
         !yed = 0.5_DP*( yn(e(ied1,ie)) + yn(e(ied2,ie)) )
         !dxl = xn(e(i,ie)) - xed
         !dyl = yn(e(i,ie)) - yed
         !elen = dsqrt(dxl*dxl + dyl*dyl)

         ! alternative nodal element-length = 2*area/length of opposite edge
         ! the results from both ways are the same.
         ! or can be understood as 1/( sqrt(dNdx^2 + dNdy^2) )

         delx = ( xn(e(ied1,ie)) - xn(e(ied2,ie)) )
         dely = ( yn(e(ied1,ie)) - yn(e(ied2,ie)) )
         opplen = dsqrt(delx**2 + dely**2)
         elen = 2.0_DP*areaf(ie)/opplen
! shortest element lenght for each element         
         elh(ie) = min(elh(ie),elen)

! shortest characteristic (element)length for all dofs
         alen(ii) = min(alen(ii), elen)
      end do ! i
   else
      do i = 1, nff
         ii = eg(ie,1)%d(i)
         if (i .LE. 3) then ! if i<=3
         ! nodal element-length = distance from vertex to mid of the opposite edge 
            ied1 = inach(i)
            ied2 = ivorg(i)      
            xed = 0.5_DP*( xn(e(ied1,ie)) + xn(e(ied2,ie)) )
            yed = 0.5_DP*( yn(e(ied1,ie)) + yn(e(ied2,ie)) )
            dxl = xn(e(i,ie)) - xed
            dyl = yn(e(i,ie)) - yed
            elen = dsqrt(dxl*dxl + dyl*dyl)

         ! alternative nodal element-length = 2*area/length of opposite edge
            !delx = ( xn(e(ied1,ie)) - xn(e(ied2,ie)) )
            !dely = ( yn(e(ied1,ie)) - yn(e(ied2,ie)) )
            !opplen = dsqrt(delx**2 + dely**2)
            !elen = 2.0_DP*areaf(ie)/opplen

            ! shortest characteristic element length of a vertex from all adjacent elements
            alen(ii) = min(alen(ii), elen)
            ! element length of vertex at ie-th element
            vlen(i) = elen
            shelleng = min(shelleng,elen)
         else
            do ! do
               if (pp .GT. polyhi) exit
               if (pp .gt. 1) then
                  k = pp*(pp+1)/2
               else
                  k = 0
               end if
            ! find the shortest length of vertices first
               j = i - k
               if (j .LE. 3) then
               ! edge dof
                  alen(ii) = min(alen(ii),vlen(j)) 
               else
               ! if inner dof
                  alen(ii) = shelleng                  
               end if            
               pp = pp + 1
            end do ! do
         end if ! if i<=3
      end do ! i
   end if   
!*******************************************************************************
! set boundary information for each element on boundaries
!*******************************************************************************
!  determine the branch of vertices
   
   do ivertex=1,3
      allocate(dirichletbc(ie)%bcdof(ivertex)%varmn(nvar))
      allocate(dirichletbc(ie)%bcdof(ivertex)%brmn(nvar))
      allocate(generalbc(ie)%bcdof(ivertex)%varmn(nvar))
      generalbc(ie)%bcdof(ivertex)%varmn(:) = .false.
      dirichletbc(ie)%bcdof(ivertex)%varmn(:) =.false.
      kziver = kzi(e(ivertex,ie))
      if (kziver .lt. 0) then
         ! if corresponding vertex is a keypoint
         branchv(ivertex) = kzrb(-kziver,1)
      else
         branchv(ivertex) = kziver
      end if      
      ! store branch on which the vertex locates
      dirichletbc(ie)%bcdof(ivertex)%branch = branchv(ivertex)
      generalbc(ie)%bcdof(ivertex)%branch = branchv(ivertex)

      do ivar = 1,nvar
         if (kziver .lt. 0) then
            ! if corresponding vertex is a keypoint
            branchvmn(ivertex,ivar) = kzrbmn(-kziver,ivar)
         else
            branchvmn(ivertex,ivar) = kziver
         end if
         brvmn = branchvmn(ivertex,ivar)                 
         dirichletbc(ie)%bcdof(ivertex)%brmn(ivar) = brvmn
         if (brvmn .ne. 0) then
            zrbvertex = zrbmn(brvmn,ivar)

            if (ivar .EQ. 1) then
               if (zrbvertex .EQ. 0) then
               ! only pressure 
                  dirichletbc(ie)%bcdof(ivertex)%varmn(ivar) = .true.
                  edbc(ie) = .true.
               end if

               if (npbc .GT. 0) then
                  pbc : do ipbc = 1,npbc
                     if (kziver .EQ. -pkp(ipbc)) then
                        if (zrbpbc(ipbc) .EQ. 0) then
                           dirichletbc(ie)%bcdof(ivertex)%varmn(ivar) = .true.
                           edbc(ie) = .true.
                        end if
                        cycle pbc
                     end if
                  end do pbc ! ipbc
               end if
               cycle
            end if

            if (zrbvertex.ge. 0 .and. zrbvertex.lt. 100) then            
               dirichletbc(ie)%bcdof(ivertex)%varmn(ivar) = .true.
               edbc(ie) = .true.
            end if
         end if
      end do ! ivar
   end do ! vertex

!  determine the branch of edges

   do iedge=1,3      
      allocate(generalbc(ie)%bcedge(iedge)%varmn(nvar))
      allocate(dirichletbc(ie)%bcedge(iedge)%varmn(nvar))
      generalbc(ie)%bcedge(iedge)%varmn(:) = .false.
      generalbc(ie)%bcedge(iedge)%branch = 0
      dirichletbc(ie)%bcedge(iedge)%varmn(:) = .false.
      dirichletbc(ie)%bcedge(iedge)%branch = 0
!      leng = .false.

!  test whether the edge is a boundary edge
      if (en(iedge,ie) .ne. 0) THEN
         branchemn(iedge,:)=0
         branchef(iedge,ie)=0
      else
!  identify the boundary condition of the edge from those of vertices
         nr=inach(iedge)
         nl=ivorg(iedge)

         lengthby2(iedge,ie) = sqrt((xn(e(nr,ie))-xn(e(nl,ie)))**2 +       &
                                    (yn(e(nr,ie))-yn(e(nl,ie)))**2 )/2._DP 

! set boundary side for fluid solver
         zweigr=branchv(nr)
         zweigl=branchv(nl)
         kziverright = kzi(e(nr,ie))
         kziverleft = kzi(e(nl,ie))
         
         if (kziverright .LT. 0 .OR. kziverleft .LT. 0) then
         
            if (kziverright .LT. 0) then
               branchef(iedge,ie) = zweigl
            else if (kziverleft .LT. 0) then
               branchef(iedge,ie) = zweigr
            end if       
         else
            if (zrb(zweigr,1) .le. zrb(zweigl,1)) then           
               branchef(iedge,ie)=zweigl
            else
               branchef(iedge,ie)=zweigr
            end if
         end if
         bredge = branchef(iedge,ie)
         ! store dirichlet branch on which the edge locates         
         dirichletbc(ie)%bcedge(iedge)%branch = bredge
         generalbc(ie)%bcedge(iedge)%branch = bredge
         
         do ivar = 1,nvar ! loop over different variables
            zweigright=branchvmn(nr,ivar)
            zweigleft=branchvmn(nl,ivar)

!  set the branch according to priority rule
            if (zrbmn(zweigright,ivar) .le. zrbmn(zweigleft,ivar)) then
               branchemn(iedge,ivar)=zweigleft
            else
               branchemn(iedge,ivar)=zweigright
            end if          
            if (branchemn(iedge,ivar) .EQ. 0) cycle
            zrbtemp = zrbmn(branchemn(iedge,ivar),ivar)

            if (ivar .eq. 1) then
            ! for pressure. dirichlet bc is only imposed by a value.
               if (zrbtemp .eq. 0) then
                  dirichletbc(ie)%bcedge(iedge)%varmn(ivar) = .true.
                  edbc(ie) = .true.
               else if (zrbtemp .ge. 200 .and. zrbtemp .lt. 300) then
                  generalbc(ie)%bcedge(iedge)%varmn(ivar) = .true.
                  if (.not. genbc(ie)) genbc(ie) = .true.
               end if
               cycle
            end if

            if (zrbtemp .ge. 0 .and. zrbtemp .lt. 100) then
               dirichletbc(ie)%bcedge(iedge)%varmn(ivar) = .true.
               edbc(ie) = .true.
               
            else if (zrbtemp .ge. 200 .and. zrbtemp .lt. 300) then
               generalbc(ie)%bcedge(iedge)%varmn(ivar) = .true.
               if (.not. genbc(ie)) genbc(ie) = .true.
            end if
         end do ! ivar
      end if ! check if edge on boundary
   end do ! iedge
! test only for incompressible flow
  if (.NOT. optcomprs) then
! set up boundary array (bcside(4,nbound),1 = elem no., 2 = node1, 3 = node2, 4 = bc code)
   do iedge = 1,3
      neighbor = en(iedge,ie)
      if (neighbor .NE. 0) cycle
      nbound = nbound + 1
      if (maxbcsize .LT. nbound) then
         write(*,*) 'Error: Size of array BCSIDE(4,:) is exceeded.'
         write(*,*) '(in geteleminfo.90)'
         STOP
      end if
      bcside(1,nbound) = ie
      bcside(2,nbound) = iedge
      bcside(3,nbound) = e(inach(iedge),ie)
      bcside(4,nbound) = e(ivorg(iedge),ie)
      bcside(5,nbound) = zrb(branchef(iedge,ie),1)
   end do ! iedge
  end if ! optcomprs
end do ! ie

! reallocate size of array "bcside" from maxbcsize to nbound
bcside => reallocate(bcside,5,nbound)

if (.NOT. associated(bcside)) then
   write(*,*) 'Error: Array pointer "bcside" cannot reallocated.'
   STOP
end if
deallocate(branchemn, branchvmn)

allocate(nside(3,nbound))
do i = 1,nbound
   ie = bcside(1,i)
   ied1 = bcside(3,i)
   ied2 = bcside(4,i)
   dxl = xn(ied2)-xn(ied1)
   dyl = yn(ied2)-yn(ied1)
   bleng = dsqrt(dxl*dxl+dyl*dyl)
   nside(1,i) = dyl/bleng
   nside(2,i) = -dxl/bleng
   nside(3,i) = bleng
end do !i

end subroutine geteleminfo