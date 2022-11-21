subroutine elementmatrix_stokes(elem,KE,FE,matvar)
use femtypes
use feminterface, only: get2Dintegpoints, shapefunction, get1Dintegpoints, & 
&                       get1Dintpolpoints, lusolver
use globalvariables
use fluidvariables
use fluidinterface
implicit none
integer(I4B), intent(in) :: elem
real(DP), pointer :: KE(:,:), FE(:)
logical, intent(in) :: matvar
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
!    $Revision: 1.8 $
!    $Date: 2014/07/15 13:02:35 $
!    $Author: m_kasper $
!
!****************************************************************************************
!
!  Calculate Element Matrix for Stokes Flow with Equal-order method
!
!****************************************************************************************
! remarks   : IMPORTANT PLEASE READ BEFORE USE
!
! 1. From experience, delta t for stationary Stokes flow is not required to be 
!    calculated!!!.
! 2. Use optts = -1 to set global time step = dtfix.
!    (set optts in the inputfluid.dat in the project folder.)
! 3. The stability of solution is obtained when delte is set to be around from 
!    1.0e-7 to 1.0e-8 (similar to the penalty method)
!****************************************************************************************
!  Computation of the element matrix for a single element
!  using numerical integration
!
!  Input:    elem     index of the element
!            x        actual solution vector
!            xn,yn    coordinates of triangle nodes
!            e        element information (nodes of the element)
!            matvar   =.true. if the material coefficients are varying across the element
!                     =.false. the material is assumed to be constant 
!
!  Output:   KE       element matrix  or  Jacobian matrix
!            FE       right hand side or residual
!            neldof   size of the element matrix and rhs-vector
!****************************************************************************************
! multi-nature variables
!****************************************************************************************
! zrbmn(gzz,nvar)    = b.c. of egdes corresponding to variables
! gzz                = total number of branches
! nvar               = number of variables
! alrbmn(gzz,nvar)   = alpha b.c. values in du/dn = alpha + beta*u 
! btrbmn(gzz,nvar)   = beta b.c. values in du/dn = alpha + beta*u
!########################################################################################
!   KE(neldof,neldof)  =   system matrix
!   FE(neldof)         =   system load vector (rhs vector)
!
!        |  K1  0   K2T |            | F1 |
!        |              |            |    |
!   KE = |  0   K1  K3T |        FE =| F2 |
!        |              |            |    |
!        |  K2  K3  KP  |            | 0  |
!
!   K1(udof,udof), K2(udof,pdof), K3(udof,pdof), KP(pdof,pdof), F1(udof), F2(udof)
!****************************************************************************************
!  local variables
!****************************************************************************************
integer(I4B) :: a,b
integer(I4B) :: i,j,k,ii,jj
integer(I4B) :: intorder, polylo, polyhi
integer(I4B) :: err, npkt, ln, npkt1D
integer(I4B) :: nr, nl, nr_mn, nl_mn
integer(I4B) :: iedge, ivertex
integer(I4B) :: ivar, ivert, ipbc
integer(I4B) :: dummy
integer(I4B) :: polyorder, errcode, intpolpoints
integer(I4B) :: branchv, branch, zrbtemp
integer(I4B) :: udof, udof2, pdof, neldof, nff
integer(I4B), dimension(3), parameter :: inach = (/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg = (/3,1,2/)
!
real(DP) :: area, aw, length2
real(DP) :: xs, ys
real(DP) :: lam(3)
real(DP) :: pp, qq
real(DP) :: dttheta, pvalue
real(DP), pointer :: K1(:,:), K2(:,:), K3(:,:), K2T(:,:), K3T(:,:), KP(:,:)
real(DP), pointer :: FU(:), FV(:), MINV(:,:)
real(DP), pointer :: xsi(:), gxsi(:,:)
real(DP), allocatable :: xtab(:), lambda(:,:), weight(:)
real(DP), allocatable :: mat(:,:), rhs(:)
complex(DPC) :: value
logical, pointer :: dbc(:)
!****************************************************************************************
dttheta = delte(elem)*theta(1)                     
area = areaf(elem)

polyorder = ep(elem,1)
polylo = 1
polyhi = polyorder
nff    = (polyorder+1)*(polyorder+2)/2
udof   = nff
udof2  = udof*2
pdof   = nff
neldof = udof2 + pdof

allocate(KE(neldof,neldof),FE(neldof))
allocate(K1(udof,udof),K2(udof,pdof),K3(udof,pdof), KP(pdof,pdof))
allocate(K2T(pdof,udof),K3T(pdof,udof))
allocate(FU(udof),FV(udof))

K1 = 0._DP
K2 = 0._DP
K3 = 0._DP
KP = 0._DP

FU = 0._DP
FV = 0._DP

KE = 0._DP
FE = 0._DP

!****************************************************************************************
!  determine the order of numerical integration
!****************************************************************************************
intorder=2*polyorder
allocate(xsi(nff), gxsi(nff,2))
!****************************************************************************************
!  fetch numerical integration points (Gauss points)
!****************************************************************************************
call get2Dintegpoints(intorder, npkt, weight, lambda, err)
if (err .ne. 0) errcode=err

do k=1,npkt

   call shapefunction(lambda(:,k),xn(e(:,elem)),yn(e(:,elem)),     &
     &    polylo,polyhi,nff,.true.,xsi,gxsi,errcode)

   aw=area*weight(k)

   do i=1,udof
      do j=1,udof
         K1(j,i) = K1(j,i) + aw*( gxsi(i,1)*gxsi(j,1) + gxsi(i,2)*gxsi(j,2) )
         K2(j,i) = K2(j,i) - aw*xsi(j)*gxsi(i,1)
         K3(j,i) = K3(j,i) - aw*xsi(j)*gxsi(i,2)
      end do
!  right hand side contributions, integrate the function:
!             f_1*xsi_i+(nu*f_2)*grad(xsi_i)
!      FE(i)=FE(i)+f1aw*xsi(i) + q4*gxsi(i,1) + q5*gxsi(i,2)
!      FU(i) = FU(i) + f1aw*xsi(i)
!      FV(i) = FV(i) + f1aw*xsi(i)
   end do
end do
deallocate(weight,lambda)
deallocate(gxsi)

K2T = transpose(K2)
K3T = transpose(K3)

allocate(MINV(nff,nff))
MINV = 0._DP
do i = 1,nff
   ii = eg(elem,1)%d(i)
   MINV(i,i) = deltp(ii)*dmmat(ii)
end do

KP = matmul(K2,matmul(MINV,K2T)) + matmul(K3,matmul(MINV,K3T)) - K1
KP = dttheta*KP

! different KP, however, no difference can be observed.
!KP = matmul(K2,matmul(MINV,K2T)) + matmul(K3,matmul(MINV,K3T)) - tint*K1
!KP = theta(1)*KP

K1 = invre*K1

KE(1:udof,1:udof) = K1(:,:)
KE(1:udof,udof2+1:neldof) = K2T(:,:)

KE(udof+1:udof2,udof+1:udof2) = K1(:,:)
KE(udof+1:udof2,udof2+1:neldof) = K3T(:,:)
KE(udof2+1:neldof,1:udof) = K2(:,:)
KE(udof2+1:neldof,udof+1:udof2) = K3(:,:)
KE(udof2+1:neldof,udof2+1:neldof) = KP(:,:)

!FE(1:udof) = FU(:)
!FE(udof+1:udof2) = FV(:)

deallocate(K1,K2,K3,KP,K2T,K3T,FU,FV,MINV)
!****************************************************************************************
!  contributions of the boundary conditions: general boundary condition
!****************************************************************************************
!  determine the number of Gauss-Lobatto points such that a polynomial of 
!  order=2*degree is integarate exactly
if (genbc(elem)) then
   call get1Dintegpoints(intorder, xtab, weight, npkt1D, errcode)

   do iedge=1,3       
      if (branchef(iedge,elem) .EQ. 0) cycle
        
      nr=inach(iedge)
      nl=ivorg(iedge)
      length2 = lengthby2(iedge,elem)
      lam(iedge)=0._DP

      do ivar = 1,3
         if (.NOT. generalbc(elem)%bcedge(iedge)%varmn(ivar)) cycle
         branch = branchef(iedge,elem)

         zrbtemp = zrbmn(branch,ivar)

         if (zrbtemp .lt. 200 .and. zrbtemp .ge. 300) cycle 

!  we assume that the coefficients of the general BC are constant along the edge
         pp = real(alrbmn(branch,ivar),DP)
         qq = real(btrbmn(branch,ivar),DP)

!  now integrate with the boundary conditions of the input branch: branche(iedge)
          do k=1, npkt1D
!  get shape function at location of integration points
            lam(nr)=(xtab(k)+1._DP)/2._DP
            lam(nl)=1._DP-lam(nr)
            call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),         &
            &                  polylo,polyhi,nff,.false.,xsi,gxsi,errcode)
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
                  KE(j,i) = KE(j,i) - length2*qq*weight(k)*xsi(ii)*xsi(jj)
               end do
!  right hand side contributions, integrate the function: (p*xsi_i)
               FE(i) = FE(i) + length2*pp*weight(k)*xsi(ii)
!  line integral of pressure on the rhs of velocity equations               
               if (ivar == 1) then                  
                  FE(i-udof)  = FE(i-udof) + length2*pp*weight(k)*xsi(ii)
                  FE(i-udof2) = FE(i-udof2) + length2*pp*weight(k)*xsi(ii)                           
               end if
            end do
         end do
      end do ! ivar
   end do ! iedge
   deallocate(xtab, weight)
end if ! if genbc
!****************************************************************************************
!  now modify the element matrix at Dirichlet BC
!****************************************************************************************
allocate(dbc(neldof))
dbc(1:neldof)=.false.

if (edbc(elem)) then

   do ivar = 1,3
      do ivertex=1,3
         if (.NOT. dirichletbc(elem)%bcdof(ivertex)%varmn(ivar)) cycle
         
         branchv = dirichletbc(elem)%bcdof(ivertex)%brmn(ivar)

         if (zrbmn(branchv,ivar).ge.0 .and. zrbmn(branchv,ivar).lt.100) then
            xs = xn(e(ivertex,elem))
            ys = yn(e(ivertex,elem))

! dirichlet_mn is the same as dirichlet, only value = alrbmn(branchvmn,ivar)
            call getbcval2D_mn(ivar,branchv,xs,ys,value)
            
            lam(ivertex) = 1._DP
            lam(inach(ivertex)) = 0._DP
            lam(ivorg(ivertex)) = 0._DP
            
            call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),         &
     &        1,1,nff,.false.,xsi,gxsi,errcode)

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
! check a fixed pressure
         if (ivar .EQ. 1) then
            ivert = ivertex+udof2
            if (npbc .GT. 0) then
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
                  FE(ivert) = pvalue   
                  dbc(ivert) = .true.
               end do ! ipbc
            end if         
         end if !(ivar == 1)

      end do ! ivertex

      if (polyorder .ge. 2) then
!  for the edge dof
!  set the number of Gauss-Lobatto points to degree+1
         intpolpoints = polyorder + 1

         allocate(xtab(intpolpoints))
         call get1Dintpolpoints(intpolpoints, xtab, errcode)

         allocate(mat(intpolpoints-2,intpolpoints-2), rhs(intpolpoints-2))

         do iedge=1,3

            if (.NOT. dirichletbc(elem)%bcedge(iedge)%varmn(ivar)) cycle            
 
            branch  = branchef(iedge,elem)
            zrbtemp = zrbmn(branch,ivar)

            if ( zrbtemp .lt.0 .or. zrbtemp.ge.100 ) cycle

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
              call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),       &
     &              polylo,polyhi,nff,.false.,xsi,gxsi,errcode)

               do j = 2,polyorder             
                  ln=j*(j+1)/2+inach(iedge)
                  mat(i-1,j-1)=xsi(ln)
               end do
!  setup right hand side
               xs= dot_product(xn(e(:,elem)),lam(:))
               ys= dot_product(yn(e(:,elem)),lam(:))
!...........
               call getbcval2D_mn(ivar,branch,xs,ys,value)
               rhs(i-1) = real(value,DP)

!  reduce by the contribution of vertex dof
               rhs(i-1) = rhs(i-1)-xsi(nr)*FE(nr_mn)-xsi(nl)*FE(nl_mn)
            end do
!  solve
            call lusolver(mat,rhs)

            do j=2,polyorder            
!  determine the local index of this dof
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
      end if ! if (polyoder .GE.2)
   end do ! ivar
end if ! edbc
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
      end do !j
   end if
end do ! i

deallocate(dbc)
deallocate(xsi)
!*****************************************************************************************
end subroutine elementmatrix_stokes
