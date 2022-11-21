subroutine elementmatrix_stokesmix(elem,polyorder,jacobi,full,matvar,KE,FE,nff,errcode)
use femtypes
use feminterface, only : get2Dintegpoints, shapefunction, get1Dintegpoints, & 
&                        get1Dintpolpoints, lusolver, pdecoeff
use fluidinterface
use fluidvariables
use globalvariables
implicit none
integer (I4B) :: elem
integer (I4B) :: polyorder, nff, errcode
real(DP), pointer :: KE(:,:), FE(:)
logical :: jacobi, full, matvar
intent (in) :: elem, polyorder, jacobi, full, matvar
intent (out) :: nff
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
!    $Revision: 1.11 $
!    $Date: 2014/07/15 13:02:35 $
!    $Author: m_kasper $
!
!****************************************************************************************
!
!  Compute Element Matrix for Stokes Flow using Mixed Method.
!
!****************************************************************************************
! remarks  : modified for Polyde
!
! - add line integral of pressure term in velocity equations
! - exchange recursions of loop over variables (ivar) and over edges (iedge)
!****************************************************************************************
!  Computation of the element matrix for a single element
!  using numerical integration
!
!  Input:    elem     index of the element
!            x        actual solution vector
!            xn,yn    coordinates of triangle nodes
!            e        element information (nodes of the element)
!            matzif   material numbers of the input regions
!            geb      region index of the elements 
!            polyorder   polynomial degree of the element
!            full     =.true. compute the matrix up to degree (always .true.)
!                     =.false. only compute the sub-matrix of degree
!            matvar   =.true. if the material coefficients are varying across the element
!                     =.false. the material is assumed to be constant 
!
!  Output:   KE        element matrix  or  Jacobian matrix
!            FE        right hand side or residual
!            neldof      size of the element matrix and rhs-vector
!                     i.e. number of form functions
!            errcode  =1000 triangle area is zero or negative
!****************************************************************************************
! multi-nature variables
!****************************************************************************************
! zrbmn(gzz,nvar)   = b.c. of egdes corresponding to variables
! gzz               = total number of branches
! nvar              = number of variables
! alrbmn(gzz,nvar)  = alpha b.c. values in du/dn = alpha + beta*u 
! btrbmn(gzz,nvar)  = beta b.c. values in du/dn = alpha + beta*u
!########################################################################################
! if different polynomial orders of elements are needed,
! need to change the number of parameter passing through the subroutine.
!########################################################################################
!   KE(neldof,neldof)   =   system matrix
!   FE(neldof)          =   system load vector (rhs vector)
!
!       |  K1   0   K2 |             | F1 |
!       |              |             |    |
!  KE = |  0    K1  K3 |         FE =| F2 |
!       |              |             |    |
!       | K2T  K3T   0 |             | 0  |
!
!   K1(udof,udof), K2(udof,pdof), K3(udof,pdof), F1(udof), F2(udof)
!
!****************************************************************************************
!  local variables
!****************************************************************************************
integer(I4B) :: a,b
integer(I4B) :: branchemn(3,3), branchvmn(3,3) ! branche(3) branchv(3) 
integer(I4B) :: i,j,k,ii,jj
integer(I4B) :: polylo, polyhi
integer(I4B) :: err, npkt, ln
integer(I4B) :: intpolpoints, intorder, npkt1D
integer(I4B) :: nr, nl, nr_mn, nl_mn
integer(I4B) :: iedge,ivertex
integer(I4B) :: ivar, ivert
integer(I4B) :: dummy, degree_mn
integer(I4B) :: zweigright, zweigleft
integer(I4B), dimension(3), parameter :: inach = (/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg = (/3,1,2/)
!integer(I4B) :: m, mn
integer(I4B) :: udof,udof2,pdof,neldof,ipbc

real(DP) :: area, aw
real(DP) :: xs, ys
real(DP) :: lengby2(3), lam(3)
real(DP) :: nuaw(2,2), alphaaw, betaaw(2), f1aw
real(DP) :: pp, qq, q1, q2, q3, q4, q5
real(DP) :: lineintegral, pvalue
real(DP), pointer :: K1(:,:), K2(:,:), K3(:,:), K2T(:,:), K3T(:,:)
real(DP), pointer :: FU(:), FV(:)
real(DP), pointer :: xsi(:), gxsi(:,:)
real(DP), allocatable :: weight(:), lambda(:,:), xtab(:)
real(DP), allocatable :: mat(:,:), rhs(:)

complex(DPC) :: value
complex(DPC) :: nu(2,2)
complex(DPC) :: alpha, beta(2), f1, f2(2)
complex (DPC) :: nupde(2,2,1,1), gammapde(2,1,1)
complex (DPC) :: alphapde(1,1), betapde(2,1,1), f1pde(1), f2pde(2,1)


logical, pointer :: ldirichletbc(:), lgeneralbc(:) ! multi-nature
logical, pointer :: dbc(:)
logical :: lgenbc
!****************************************************************************************
area=( (yn(e(3,elem))-yn(e(1,elem)))*                             &
       (xn(e(2,elem))-xn(e(3,elem)))-                             &
       (yn(e(2,elem))-yn(e(3,elem)))*                             &
       (xn(e(3,elem))-xn(e(1,elem))) )/2._DP
if (area .le. tiny(1._DP) ) then
   area = 0._DP
   errcode=1000
   return
end if

if (full) then
   polyhi=polyorder
   polylo=1
!  size of the element matrix
   nff=(polyorder+1)*(polyorder+2)/2
else
   polyhi=polyorder
   polylo=polyorder
!  size of the element matrix
   if (polyorder .gt.1) then
      nff=polyorder+1
   else
      nff=3
   end if
end if

udof = nff
udof2 = 2*udof
pdof = polyorder*(polyorder+1)/2
neldof = udof2 + pdof

allocate(KE(neldof,neldof),FE(neldof))
allocate(K1(udof,udof),K2(udof,pdof),K3(udof,pdof))
allocate(K2T(pdof,udof),K3T(pdof,udof))
allocate(FU(udof),FV(udof))

K1 = 0._DP
K2 = 0._DP
K3 = 0._DP
FU = 0._DP
FV = 0._DP

KE = 0._DP
FE = 0._DP

allocate(xsi(udof),gxsi(udof,2))
!****************************************************************************************
!  fetch material coefficient at the center of gravity
!****************************************************************************************
if (.not.matvar) then
!  if the material is constant across the element
   xs= ( xn(e(1,elem)) + xn(e(2,elem)) + xn(e(3,elem)) ) / 3._DP
   ys= ( yn(e(1,elem)) + yn(e(2,elem)) + yn(e(3,elem)) ) / 3._DP
   call pdecoeff(elem,xs,ys,nupde,gammapde,alphapde,betapde,f1pde,f2pde)
   nu=nupde(:,:,1,1)
   alpha=alphapde(1,1)
   beta=betapde(:,1,1)
   f1=f1pde(1)
   f2=f2pde(:,1)
end if

!  determine the order of numerical integration
intorder=2*polyorder
!****************************************************************************************
!  fetch numerical integration points (Gauss points)
!****************************************************************************************
call get2Dintegpoints(intorder, npkt, weight, lambda, err)

do k=1,npkt

!  if the material is vaying across the element
   if (matvar) then
!  compute the cartesian coordinates and fetch the material coefficients
      xs= dot_product( xn(e(:,elem)) , lambda(:,k) )
      ys= dot_product( yn(e(:,elem)) , lambda(:,k) )
      call pdecoeff(elem,xs,ys,nupde,gammapde,alphapde,betapde,f1pde,f2pde)
      nu=nupde(:,:,1,1)
      alpha=alphapde(1,1)
      beta=betapde(:,1,1)
      f1=f1pde(1)
      f2=f2pde(:,1)
   end if

!  get shape function and gradients at location of integration points
   call shapefunction(lambda(:,k),xn(e(:,elem)),yn(e(:,elem)),     &
         polylo,polyhi,nff,.true.,xsi,gxsi,errcode)

   aw = area*weight(k)
   nuaw = real(nu,DP)*aw
   alphaaw = real(alpha,DP)*aw
   betaaw = real(beta,DP)*aw
   f1aw = real(f1,DP)*aw
   q4 = real(f2(1),DP)*aw
   q5 = real(f2(2),DP)*aw

!  for all shape functions
   do i=1,udof

      q1 = nuaw(1,1)*gxsi(i,1)+ nuaw(1,2)*gxsi(i,2)
      q2 = nuaw(2,1)*gxsi(i,1)+ nuaw(2,2)*gxsi(i,2)
      q3 = betaaw(1)*gxsi(i,1)+ betaaw(2)*gxsi(i,2) + alphaaw*xsi(i)

      do j=1,udof
!  integrate the function:
!  (nu*grad(xsi_i)*grad(xsi_j)+(beta*grad(xsi_i)*xsi_j+alpha*(xsi_i)*(xsi_j)

      K1(j,i)=K1(j,i) + q1*gxsi(j,1) + q2*gxsi(j,2) + q3*xsi(j)

         if (j .LE. pdof) then
            K2(i,j) = K2(i,j) - aw*xsi(j)*gxsi(i,1)
            K3(i,j) = K3(i,j) - aw*xsi(j)*gxsi(i,2)
         end if

      end do
!  right hand side contributions, integrate the function:
!             f_1*xsi_i+(nu*f_2)*grad(xsi_i)
!      FE(i)=FE(i)+f1aw*xsi(i) + q4*gxsi(i,1) + q5*gxsi(i,2)
      FU(i) = FU(i) + f1aw*xsi(i) + q4*gxsi(i,1) + q5*gxsi(i,2)
      FV(i) = FV(i) + f1aw*xsi(i) + q4*gxsi(i,1) + q5*gxsi(i,2)
   end do
end do
deallocate(weight,lambda)
deallocate(gxsi)

K2T = transpose(K2)
K3T = transpose(K3)

KE(1:udof,1:udof) = K1(:,:)
KE(1:udof,udof2+1:neldof) = K2(:,:)

KE(udof+1:udof2,udof+1:udof2) = K1(:,:)
KE(udof+1:udof2,udof2+1:neldof) = K3(:,:)

KE(udof2+1:neldof,1:udof) = K2T(:,:)
KE(udof2+1:neldof,udof+1:udof2) = K3T(:,:)

FE(1:udof) = FU(:)
FE(udof+1:udof2) = FV(:)

deallocate(K1,K2,K3,K2T,K3T,FU,FV)
!****************************************************************************************
!  contributions of the boundary conditions: prepare
!****************************************************************************************
!  determine the branch of vertices
!****************************************************************************************
allocate(ldirichletbc(3)) 
ldirichletbc(:)=.false.

do ivertex=1,3
   do ivar = 1,3  ! loop over variables, u,v,p
      if (kzi(e(ivertex,elem)) .lt. 0) then
!         branchv(ivertex)=kzrb(-kzi(e(ivertex,elem)))
         branchvmn(ivertex,ivar) = kzrbmn(-kzi(e(ivertex,elem)),ivar)
      else
         branchvmn(ivertex,ivar) = kzi(e(ivertex,elem))
      end if
    
      if (branchvmn(ivertex,ivar) .ne. 0) then
         if (zrbmn(branchvmn(ivertex,ivar),ivar).ge. 0 .and. & 
            zrbmn(branchvmn(ivertex,ivar),ivar).lt. 100) then
            ldirichletbc(ivar) = .true.
         end if
      end if
   end do   ! ivar
end do   ! vertex
!****************************************************************************************
!  determine the branch of edges
!****************************************************************************************
allocate(lgeneralbc(3))
lgeneralbc(:)=.false.
lgenbc = .FALSE.

do iedge=1,3

! test whether the edge is a boundary edge
   if (en(iedge,elem) .ne. 0) then
      branchemn(iedge,:)=0
   else
!  identify the boundary condition of the edge from those of vertices
      nr=inach(iedge)
      nl=ivorg(iedge)

      lengby2(iedge)=sqrt((xn(e(nr,elem))-xn(e(nl,elem)))**2 + &
                           (yn(e(nr,elem))-yn(e(nl,elem)))**2 )/2._DP

      lgenbc = .TRUE.

      do ivar = 1,3  ! loop variables
         zweigright=branchvmn(nr,ivar)
         zweigleft=branchvmn(nl,ivar)

!  set the branch according to priority rule
         if (zrbmn(zweigright,ivar) .le. zrbmn(zweigleft,ivar)) then
            branchemn(iedge,ivar)=zweigleft
         else
            branchemn(iedge,ivar)=zweigright
         end if

         if (zrbmn(branchemn(iedge,ivar),ivar).ge.200 .and. & 
            zrbmn(branchemn(iedge,ivar),ivar).lt.300) then
            lgeneralbc(ivar)=.true.
         end if
      end do ! ivar
   end if
end do
!****************************************************************************************
!  contributions of the boundary conditions: general boundary condition
!****************************************************************************************
!  determine the number of Gauss-Lobatto points such that a polynomial of 
!  order=2*polyorder is integarate exactly
if (lgenbc) then

   call get1Dintegpoints(intorder, xtab, weight, npkt1D, errcode)

   do iedge=1,3
      do ivar = 1,3 ! multi-nature
         if (.NOT. lgeneralbc(ivar)) cycle
!  check whether it is a general or Neumann BC
         if (branchemn(iedge,ivar) .eq. 0) cycle

         if (zrbmn(branchemn(iedge,ivar),ivar).lt.200 .or. &
            zrbmn(branchemn(iedge,ivar),ivar).ge.300) cycle

!  we assume that the coefficients of the general BC are constant along the edge
         pp=real(alrbmn(branchemn(iedge,ivar),ivar),DP)
         qq=real(btrbmn(branchemn(iedge,ivar),ivar),DP)

         nr=inach(iedge)
         nl=ivorg(iedge)
         lam(iedge)=0._DP

!  now integrate with the boundary conditions of the input branch: branche(iedge)
         do k=1,npkt1D

!  get shape function at location of integration points
            lam(nr)=(xtab(k)+1._DP)/2._DP
            lam(nl)=1._DP-lam(nr)
            call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),     &
                  polylo,polyhi,nff,.false.,xsi,gxsi,errcode)
!  for all shape functions
            select case (ivar)
               case(2)
                  a = 1
                  b = udof
               case(3)
                  a = udof+1
                  b = udof2
               case(1)
                  a = udof2+1
                  b = neldof
            end select
            do i=a,b
               do j=a,b
                  select case (ivar)
                     case(2)
                        ii = i
                        jj = j
                     case(3)
                        ii = i-udof
                        jj = j-udof
                     case(1)
                        ii = i-udof2
                        jj = j-udof2
                  end select
!  integrate the function: (q*xsi_i*xsi_j)
                  KE(j,i)=KE(j,i)-lengby2(iedge)*qq*weight(k)*xsi(ii)*xsi(jj)
               end do ! j
!  right hand side contributions, integrate the function: (p*xsi_i)
               FE(i)=FE(i) + lengby2(iedge)*pp*weight(k)*xsi(ii) 
!  line integral of pressure on the rhs of velocity equations               
               if (ivar == 1) then
                  lineintegral = lengby2(iedge)*pp*weight(k)*xsi(ii)
                  FE(i-udof)  = FE(i-udof) + lineintegral
                  FE(i-udof2) = FE(i-udof2) + lineintegral                        
               end if
            end do ! i
         end do ! k
      end do ! multi-nature
   end do ! loop over edges
   deallocate(xtab, weight)
end if ! if lgenbc
!****************************************************************************************
!  now modify the element matrix at Dirichlet BC
!****************************************************************************************
allocate(dbc(neldof))
dbc(1:neldof)=.false.

do ivar = 1,3

   if (.NOT. ldirichletbc(ivar)) cycle ! if ldirichletbc = .true.

   do ivertex = 1,3

      if (branchvmn(ivertex,ivar) .eq. 0) cycle

      if (zrbmn(branchvmn(ivertex,ivar),ivar).ge.0 .and. &
          zrbmn(branchvmn(ivertex,ivar),ivar).lt.100) then

         xs= xn(e(ivertex,elem))
         ys= yn(e(ivertex,elem))

! dirichlet_mn is the same as dirichlet, only value = alrbmn(branchvmn,ivar)
         call getbcval2D_mn(ivar,branchvmn(ivertex,ivar),xs,ys,value)
            
         lam(ivertex)=1._DP
         lam(inach(ivertex))=0._DP
         lam(ivorg(ivertex))=0._DP
            
         call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),     &
                            polylo,polyhi,nff,.false.,xsi,gxsi,errcode)

         select case (ivar) ! to put related values to the correct local dof
         case(2)
            ivert = ivertex
         case(3)
            ivert = ivertex+udof
         case(1)
            ivert = ivertex+udof2
         end select

         FE(ivert) = real(value,DP)/xsi(ivertex)
         dbc(ivert) = .true.
      end if

      if (ivar .eq. 1 .AND. npbc .GT. 0) then
         jj = eg(elem,1)%d(ivertex)
         do ipbc = 1,npbc
            if (-(kzi(jj)) .NE. pkp(ipbc)) cycle
            ! want to have only p = 0 at a fixed point.
            if (zrbpbc(ipbc) .NE. 0) then
               write(*,*) 'Warning: fixed value of pressure was not imposed'
               cycle
            end if      
            pvalue = varpbc(npbc)
            ! only p = 0 needs to be imposed.
            FE(ivertex+udof2) = pvalue   
            dbc(ivertex+udof2) = .true.
         end do ! ipbc
      end if ! fixed pressure

   end do   ! ivertex

! to check if polynomial order of pressure is greater than 2
   if (ivar == 1) then
      degree_mn = polyorder-1
   else
      degree_mn = polyorder
   end if
   if (degree_mn .ge. 2) then
!  for the edge dof
!  set the number of Gauss-Lobatto points to polyorder+1
      if (ivar == 1) then ! for pressure
         intpolpoints=polyorder
         degree_mn = polyorder-1
      else ! for velocity
         intpolpoints = polyorder+1
         degree_mn = polyorder
      end if

      allocate(xtab(intpolpoints))
      call get1Dintpolpoints(intpolpoints, xtab, errcode)  
      
      allocate(mat(intpolpoints-2,intpolpoints-2), rhs(intpolpoints-2))

      do iedge=1,3

         if (branchemn(iedge,ivar) .eq. 0) cycle

         if (zrbmn(branchemn(iedge,ivar),ivar).lt.0 .or. &
             zrbmn(branchemn(iedge,ivar),ivar).ge. 100) cycle

            nr=inach(iedge)
            nl=ivorg(iedge)
            lam(iedge)=0._DP

         select case (ivar)
         case(2)
            nr_mn = nr
            nl_mn = nl
         case(3)
            nr_mn = nr+udof
            nl_mn = nl+udof
         case(1)
            nr_mn = nr+udof2
            nl_mn = nl+udof2
         end select
            
         do i=2,intpolpoints-1
!  evaluate BC at Gauss-Lobatto points
            lam(nr)=(xtab(i)+1._DP)/2._DP
            lam(nl)=1._DP-lam(nr)

            call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),     &
                               polylo,polyhi,nff,.false.,xsi,gxsi,errcode)

            do j = 2,degree_mn
               ln=j*(j+1)/2+inach(iedge)
               mat(i-1,j-1)=xsi(ln)
            end do
!  setup right hand side 
            xs= dot_product(xn(e(:,elem)),lam(:))
            ys= dot_product(yn(e(:,elem)),lam(:))

            call getbcval2D_mn(ivar,branchemn(iedge,ivar),xs,ys, value)
            rhs(i-1) = real(value,DP)
!  reduce by the contribution of vertex dof
            rhs(i-1) = rhs(i-1)-xsi(nr)*FE(nr_mn)-xsi(nl)*FE(nl_mn)
         end do
!  solve
         call lusolver(mat,rhs)

         do j=2,degree_mn
            select case (ivar)
            case(2)
               dummy = 0
            case(3)
               dummy = udof
            case(1)
               dummy = udof2
            end select

            ln=j*(j+1)/2+inach(iedge)+dummy

            FE(ln)=rhs(j-1)
            dbc(ln)=.true.
         end do
      end do
      deallocate(mat,rhs)
   end if ! if (degree_mn .GE.2)

!  modify the element matrix to enforce essential (Dirichlet) BC
   do i=1,neldof
      if (dbc(i)) then
         KE(i,1:neldof)=0._DP
         KE(i,i)=1._DP
      else 
         do j=1,neldof
            if (dbc(j)) then
               FE(i)=FE(i)-KE(i,j)*FE(j)
               KE(i,j)=0._DP
            end if
         end do
      end if
   end do
end do ! ivar

deallocate(dbc)
deallocate(ldirichletbc,lgeneralbc) 
deallocate(xsi)
!****************************************************************************************
end subroutine elementmatrix_stokesmix