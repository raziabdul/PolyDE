      subroutine elementmatrix(elem, polyorder, jacobi, full, matvar, a, b, nff, nffsum, errcode)
      use feminterface, only: pdecoeff, get2Dintegpoints, get1Dintegpoints
      use feminterface, only: lusolver, shapefunction, get1Dintpolpoints
      use feminterface, only: getbcval2D, lam2xy
      use femtypes
      use globalvariables, only: x, xn, yn, e, matzif, geb, kzi, kzrb, en, nnat, zrb
      implicit none
      integer (I4B) :: elem
      integer (I4B) :: errcode
      integer (I4B) :: polyorder(:)
      integer (I4B), dimension(:), allocatable :: nff, nffsum
      complex(DPC), pointer :: a(:,:), b(:)
      logical :: jacobi, full, matvar
      intent (in) :: elem, polyorder, jacobi, full, matvar
    !  intent (out) :: nff, nffsum  ! org
      intent (inout) :: nff, nffsum ! test
      !
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
!    $Revision: 1.27 $
!    $Date: 2014/08/22 11:13:46 $
!    $Author: m_kasper $
!
!
!  Computation of the element matrix for a single element
!  using numerical integration
!
!  Input:    elem     index of the element
!            x        actual solution vector
!            xn,yn    coordinates of triangle nodes
!            e        element information (nodes of the element)
!            matzif   material numbers of the input regions
!            geb      region index of the elements 
!            polyorder   polynomial order of the element
!                     / = .true. for the computation of the jacobian
!            jacobi
!                     \ = .false. for the linear case
!            full     = .true.  compute the matrix up to polyorder
!                     = .false. only compute the sub-matrix of polyorder
!            matvar   =.true.   if the material coefficients are varying across the element
!                     =.false.  the material is assumed to be constant
!
!  Output:   a        element matrix  or  Jacobian matrix
!            b        right hand side or residual
!            nff      size of the element matrix and rhs-vector
!                     i.e. number of form functions
!            nffsum   number to add when we want to consider entries in
!                     the submatrix blocks with index inat or jnat higher than 1
!            errcode  =1000 triangle area is zero or negative
!
!  local variables
!
      integer (I4B) i, j, k, err, npkt, ln
      integer (I4B), allocatable :: polylo(:), polyhi(:)
      integer (I4B) iedge, ivertex, zweigright,  zweigleft 
      integer (I4B) intpolpoints, intorder, npkt1D
      integer (I4B) branchvd(3,nnat), brancheg(3,nnat)
      integer (I4B) inach(3), ivorg(3), nr, nl
      integer (I4B) inat, jnat
      integer (I4B) startnode, endnode
      real(DP) area, xs, ys, lengthby2(3), lam(3), aw
      real (DP) :: startpoint(2), endpoint(2), length(3), vec(2), nvec(2)
      real (DP), allocatable :: xsi(:), gxsi(:,:), xtab(:), weight1d(:)
      real (DP), allocatable :: lambda2d(:,:), weight2d(:)
!  pde coefficients
      complex (DPC) nu(2,2,nnat,nnat), gamma(2,nnat,nnat), alpha(nnat,nnat)
      complex (DPC) beta(2,nnat,nnat), f1(nnat), f2(2,nnat)
!  weighted pde coefficients
      complex (DPC) nuaw(2,2), gammaaw(2), alphaaw
      complex (DPC) betaaw(2), f1aw
      complex (DPC) pp(nnat), qq(nnat,nnat), gammanormal, value, q1, q2, q3, q4, q5
      complex (DPC), allocatable :: mat(:,:), rhs(:)
      logical dirichletbc, generalbc
      logical, pointer :: dbc(:)
      parameter (inach=(/2,3,1/))
      parameter (ivorg=(/3,1,2/))
!  there are two possible variants in numerical integration
!  i )  either using the same integration rule for all entries in the matrix,
!       which offers the advantage of evaluating all shape functions at the same locations 
!       thus reducing the work required in the evaluation of shape function
!  ii)  using different integration formulae for the different matrix entries
!       this has the advantage, that the accuracy degree of the integration formula 
!       can be adapted to the requirements and order of the polynomial under investigation
!  unfortunately there is no (real) possibility to profit from both, since this would
!  require embedded (eventually adaptive) Gaussian integration formulae over the element.
!
!  here we use the first approach
!
      errcode=0
!
! deallocate first. Avoids error in oneAPI
      deallocate(nff)
      deallocate(nffsum)

!  Triangle area calculation
      area=( (yn(e(3,elem))-yn(e(1,elem)))*                             &
     &       (xn(e(2,elem))-xn(e(3,elem)))-                             &
     &       (yn(e(2,elem))-yn(e(3,elem)))*                             &
     &       (xn(e(3,elem))-xn(e(1,elem))) )/2._DP                      
      if (area .le. tiny(1._DP) ) then
        area = 0._DP
        errcode=1000
        return
      end if
!
!  polyhi and polylo for every nature
      allocate(polyhi(nnat), polylo(nnat))
      allocate(nff(1:nnat))
!  size of the element matrix
      do inat=1,nnat
        if (full) then
          polyhi(inat)=polyorder(inat)
          polylo(inat)=1
          nff(inat)=(polyorder(inat)+1)*(polyorder(inat)+2)/2
        else
          polyhi(inat)=polyorder(inat)
          polylo(inat)=polyorder(inat)
          if (polyorder(inat) .gt.1) then
            nff(inat)=polyorder(inat)+1
          else
            nff(inat)=3
          end if
        end if
      end do
!
!  nffsum definition
!  example: if nff=(/3,6,3/) then
!  nffsum(1:4)=(/0,3,9,12/)
      allocate (nffsum(1:nnat+1))
      nffsum(1)=0
      do inat=1,nnat
        nffsum(inat+1)=sum(nff(1:inat))
      end do
!  Allocate required space and initialise elementmatrix a and RHS b
      allocate( a(sum(nff),sum(nff)), b(sum(nff)), xsi(maxval(nff)), gxsi(maxval(nff),2),stat=errcode )
      a(:,:)=(0._DP,0._DP)
      b(:)=(0._DP,0._DP)
!
!  fetch material coefficient at the center of gravity
      if (.not.matvar) then
!  if the material is constant across the element
        xs= ( xn(e(1,elem)) + xn(e(2,elem)) + xn(e(3,elem)) ) / 3._DP
        ys= ( yn(e(1,elem)) + yn(e(2,elem)) + yn(e(3,elem)) ) / 3._DP
        call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2)
      end if
!
!  determine the order of numerical integration. 
!  the weak form has the from (xi_i*xi_j) such that maximal polynomial degree is 2*polyorder
!  (it is assumed that material coefficients do not vary strongly across the element)
!
      intorder=2*(maxval(polyorder))
!  fetch numerical integration points (Gauss points)
      call get2Dintegpoints(intorder, npkt, weight2d, lambda2d, err)
      if (err .ne. 0) errcode=err
!
      do k=1,npkt
!  if the material is varying across the element
        if (matvar) then
!  compute the cartesian coordinates and fetch the material coefficients
          call lam2xy(lambda2d(:,k),elem,xs,ys,xn,yn,e)
          call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2)
        end if
!  get shape function and gradients at location of integration points
!
!  polylo and polyhi replace by 1 and polyorder
        call shapefunction(lambda2d(:,k),xn(e(:,elem)),yn(e(:,elem)),     &
     &    minval(polylo),maxval(polyhi),maxval(nff),.true.,xsi,gxsi,errcode)
        aw=area*weight2d(k)
!  loop over the natures in the direction of increasing row number
        do inat=1,nnat
!  loop over the natures in the direction of increasing column number
          do jnat=1,nnat
!  Multiply material coefficients for this nature with the area times the weight
            nuaw=nu(:,:,inat,jnat)*aw
            gammaaw=gamma(:,inat,jnat)*aw
            alphaaw=alpha(inat,jnat)*aw
            betaaw=beta(:,inat,jnat)*aw
            f1aw=f1(inat)*aw
            q4=f2(1,inat)*aw
            q5=f2(2,inat)*aw
!
!  for all shape functions
!  Since every submatrix block of the elementmatrix is symmetric we can
!  outer-loop over the column index j and inner-loop over the row index i
!  this allows for some optimisation by moving common factors out of the inner loop
            do j=1,nff(jnat)
!  nu is a rank-2 tensor and its entries are multipled with
!  the gradient of xsi in x and the gradient of xsi in y
!  gamma is a vector
              q1=nuaw(1,1)*gxsi(j,1)+nuaw(1,2)*gxsi(j,2)+gammaaw(1)*xsi(j)
              q2=nuaw(2,1)*gxsi(j,1)+nuaw(2,2)*gxsi(j,2)+gammaaw(2)*xsi(j)
!  beta is a vector, alpha is a scalar
              q3=betaaw(1)*gxsi(j,1)+betaaw(2)*gxsi(j,2) +  alphaaw*xsi(j)
!
              do i=1,nff(inat)
!
!  integrate the function
!     i: row (degree of freedom);    j: column (test function)
!  (nu*grad(xsi_j)+gamma*xsi_j)*grad(xsi_i)+(beta*grad(xsi_j))*xsi_i+alpha*xsi_j*xsi_i
!
                a(i+nffsum(inat),j+nffsum(jnat)) = a(i+nffsum(inat),j+nffsum(jnat))     &
           &    + q1*gxsi(i,1) + q2*gxsi(i,2) + q3*xsi(i)
!
!  contributions to the jacobian matrix if nonlinear
!  ## to be implemented
!            if (jacobi) then
!              a(i,j)=a(i,j)+weight(k)*(whatjacobi)
!            end if
              end do
!  right hand side contributions, integrate the function:
!             f_1*xsi_i-f_2*grad(xsi_i)
!
              if (inat .eq. jnat) then
                b(j+nffsum(jnat)) = b(j+nffsum(jnat)) + f1aw*xsi(j) - q4*gxsi(j,1) - q5*gxsi(j,2)
              end if
!
            end do
          end do
        end do
      end do
      deallocate(weight2d,lambda2d)
      deallocate(gxsi)
!
!---------------------------------------------------------------------------
!  contributions of the boundary conditions: prepare
!---------------------------------------------------------------------------
!
!  determine the branch of vertices
!  examine the vertices of the element in all natures
!  to see if we need to consider dirichlet boundary conditions for
!  this element
      dirichletbc=.false.
      do inat=1,nnat
        do ivertex=1,3
          if (kzi(e(ivertex,elem)) .lt. 0) then
            branchvd(ivertex,inat)=kzrb(-kzi(e(ivertex,elem)),inat)
          else
            branchvd(ivertex,inat)=kzi(e(ivertex,elem))
          end if
          if (branchvd(ivertex,inat) .ne.  0) then
            if (zrb(branchvd(ivertex,inat),inat).ge. 0 .and.                         &
       &        zrb(branchvd(ivertex,inat),inat).lt. 100) then
              dirichletbc=.true.
            end if
          end if
        end do
      end do
!
! initialize default to brancheg to internal (=0)
! release build puts garbage into it. Logic is missing for case?
      do inat=1,nnat
        do iedge=1,3
          brancheg(iedge,inat)=0
        end do
      end do
!
!  determine the branch of edges
      generalbc=.false.
      do inat=1,nnat
        do iedge=1,3
!  Check to see if this edge is a boundary edge
!  The possibility of internal boundaries is also included 
!  if both side of the edge are in dfferent regions
!  
!  is there a neighbor?
          if (en(iedge,elem) .gt. 0) then
!  if element and its neighbor are in the same region it is not a boundary branch
            if (geb(elem) .eq. geb(en(iedge,elem)) ) then
              brancheg(iedge,inat)=0
              cycle
            end if
          end if
          nr=inach(iedge)
          nl=ivorg(iedge)
!  none of the nodes (vertices) of this edge is inner
!  identify the boundary condition of the edge from those of vertices
!  check whether one of the two vertices of the edge is a key-point
!  select the boundary condition of the non-keypoint vertex for the edge
!  since the edge belongs to this branch
          if (kzi(e(nr,elem)) .lt. 0) then
            zweigleft=branchvd(nl,inat)
            brancheg(iedge,inat)=zweigleft
          else
            zweigright=branchvd(nr,inat)
            brancheg(iedge,inat)=zweigright
          end if
!  if a general or neumann BC is present set the flag an determine the halflength
          if (zrb(brancheg(iedge,inat),inat).ge.200 .and.               &
     &        zrb(brancheg(iedge,inat),inat).lt.300) then
             generalbc=.true.
             lengthby2(iedge)=sqrt(                                     &
     &          (xn(e(nr,elem))-xn(e(nl,elem)))**2 +                    &
     &          (yn(e(nr,elem))-yn(e(nl,elem)))**2 ) / 2._DP
          end if
        end do
      end do
!
!  contributions of the boundary conditions: general or Neumann BC
!
      if (generalbc) then
!  fetch the Gauss-Legendre points such that a polynomial of
!  order of intorder is integrated exactly.
        call get1Dintegpoints(intorder, xtab, weight1d, npkt1D, errcode)
!  Edge lengths
        length(1:3) = sqrt((xn(e(inach(1:3),elem))-xn(e(ivorg(1:3),elem)))**2 + &
        &                  (yn(e(inach(1:3),elem))-yn(e(ivorg(1:3),elem)))**2 )
!  gxsi not used, but need to allocate to avoid error with oneAPI
        allocate( gxsi(maxval(nff),2) )
        do iedge=1,3

!  we assume that the coefficients of the general BC are constant along the edge
!
          nr=inach(iedge)
          nl=ivorg(iedge)
          lam(iedge)=0._DP
          startnode = e(nr,elem)
          endnode   = e(nl,elem)
          startpoint = (/xn(startnode),yn(startnode)/)
          endpoint   = (/xn(endnode) ,yn(endnode) /)
          vec = (endpoint - startpoint) / length(iedge)
          nvec = (/vec(2),-vec(1)/)
!
!  now integrate with the boundary conditions of the input branch: brancheg(iedge,inat)
          do k=1, npkt1D
!  get shape function at location of integration points
            lam(nr)=(xtab(k)+1._DP)/2._DP
            lam(nl)=1._DP-lam(nr)
            call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),         &
         &  minval(polylo),maxval(polyhi),maxval(nff),.false.,xsi,gxsi,errcode)
!
            call lam2xy(lam,elem,xs,ys,xn,yn,e)

            do inat=1,nnat
!  check whether it is a general or Neumann BC
              if (brancheg(iedge,inat) .eq. 0) cycle
              if (zrb(brancheg(iedge,inat),inat).lt.200 .or.             &
     &        zrb(brancheg(iedge,inat),inat).ge.300) cycle

              select case (zrb(brancheg(iedge,inat),inat))
                case (250)
                  do i=1,nff(inat)
              
                    do jnat=1,nnat
                      gammanormal = sum(gamma(:,inat,jnat)*nvec)
!  for all shape functions
!  we may cycle here if the shape function xsi_i is zero at this edge
                      do j=1,nff(jnat)

!  integrate the function: (q*xsi_i*xsi_j)
!  we may cycle here if the shape function xsi_j is zero at this edge
                       a(i+nffsum(inat),j+nffsum(jnat))=a(i+nffsum(inat),j+nffsum(jnat)) - &
                    &  lengthby2(iedge) * gammanormal * weight1d(k) * xsi(i)*xsi(j)

                      end do
                    end do
!  right hand side contributions, integrate the function: (p*xsi_i)
                    b(i+nffsum(inat))=b(i+nffsum(inat)) + lengthby2(iedge)*pp(inat)*weight1d(k) * xsi(i)
                  end do
                case default
                  call getbcval2D(brancheg(iedge,inat),inat,xs,ys,pp(inat),qq(inat,:),elem)
                  do i=1,nff(inat)

                    do jnat=1,nnat
!  for all shape functions
!  we may cycle here if the shape function xsi_i is zero at this edge
                      do j=1,nff(jnat)

!  integrate the function: (q*xsi_i*xsi_j)
!  we may cycle here if the shape function xsi_j is zero at this edge
                       a(i+nffsum(inat),j+nffsum(jnat))=a(i+nffsum(inat),j+nffsum(jnat)) - &
                    &  lengthby2(iedge)*qq(inat,jnat)*weight1d(k) * xsi(i)*xsi(j)

                      end do
                    end do
!  right hand side contributions, integrate the function: (p*xsi_i)
                    b(i+nffsum(inat))=b(i+nffsum(inat)) + lengthby2(iedge)*pp(inat)*weight1d(k) * xsi(i)
                  end do
              end select

            end do ! inat
          end do ! bpkt1D
        end do ! edge
        deallocate(gxsi)
        deallocate(xtab, weight1d)
      end if
!
!  now modify the element matrix at Dirichlet BC
!
      if (dirichletbc) then

        allocate(dbc(nffsum(nnat+1)))
        dbc=.false.
        do inat=1,nnat
!  for the vertices
          do ivertex=1,3
            if (branchvd(ivertex,inat) .eq. 0) cycle
            if (zrb(branchvd(ivertex,inat),inat).ge.0 .and.       &
       &        zrb(branchvd(ivertex,inat),inat).lt.100) then
              xs= xn(e(ivertex,elem))
              ys= yn(e(ivertex,elem))
              call getbcval2D(branchvd(ivertex,inat),inat,xs,ys,value)

              lam(ivertex)=1._DP
              lam(inach(ivertex))=0._DP
              lam(ivorg(ivertex))=0._DP
! Not used. Avoids runtime error in oneAPI
              allocate( gxsi(maxval(nff),2) )
              call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),         &
       &      1,1,nff(inat),.false.,xsi,gxsi,errcode)
              deallocate(gxsi)
              b(ivertex+nffsum(inat))=value/xsi(ivertex)
              dbc(ivertex+nffsum(inat))=.true.
            end if
          end do
          if (polyorder(inat) .ge. 2) then
!  for the edge dof set the number of Gauss-Lobatto points to polyorder+1
            intpolpoints=polyorder(inat)+1
            allocate(xtab(intpolpoints))
            call get1Dintpolpoints(intpolpoints, xtab, errcode)
!
            allocate(mat(intpolpoints-2,intpolpoints-2), rhs(intpolpoints-2))
            do iedge=1,3
              if (brancheg(iedge,inat) .eq. 0) cycle
              if (zrb(brancheg(iedge,inat),inat).lt.0 .or.                           &
      &           zrb(brancheg(iedge,inat),inat).ge. 100) cycle
              nr=inach(iedge)
              nl=ivorg(iedge)
              lam(iedge)=0._DP
              do i=2,intpolpoints-1
!  evaluate shapefunction at Gauss-Lobatto points
                lam(nr)=(xtab(i)+1._DP)/2._DP
                lam(nl)=1._DP-lam(nr)
! Allocate but not used. Avoids oneAPI runtime error
                allocate( gxsi(maxval(nff),2) )
                call shapefunction(lam,xn(e(:,elem)),yn(e(:,elem)),       &
       &          polylo(inat),polyhi(inat),nff(inat),.false.,xsi,gxsi,errcode)
                deallocate(gxsi)
                do j=2,polyorder(inat)
!  determine the local index of this dof
!  4, 7, 11, 16, 22, .. for edge 3
!  5, 8, 12, 17, 23, .. for edge 1
!  6, 9, 13, 18, 24, .. for edge 2
!  2  3   4   5   6     for these polyorders
                  ln=j*(j+1)/2+inach(iedge)
                  mat(i-1,j-1)=xsi(ln)
                end do
!
!  setup right hand side 
                call lam2xy(lam,elem,xs,ys,xn,yn,e)
!  evaluate BC
                call getbcval2D(brancheg(iedge,inat),inat,xs,ys, rhs(i-1))

!  reduce by the contribution of vertex DOF
                rhs(i-1)=rhs(i-1) - xsi(nr)*b(nr+nffsum(inat)) - &
                &   xsi(nl)*b(nl+nffsum(inat))
              end do
!  solve
              call lusolver(mat,rhs)
              do j=2,polyorder(inat)
!  determine the local index of this dof
                ln=j*(j+1)/2+inach(iedge)
                b(ln+nffsum(inat))=rhs(j-1)
                dbc(ln+nffsum(inat))=.true.
              end do

            end do
!
            deallocate(mat, rhs, xtab)
          end if
        end do
!
!  modify the element matrix to enforce essential (Dirichlet) BC
        do i=1,sum(nff)
          if (dbc(i)) then
            a(i,1:sum(nff))=0._DP
            a(i,i)=1._DP
          else 
            do j=1,sum(nff)
              if (dbc(j)) then
                b(i)=b(i)-a(i,j)*b(j)
                a(i,j)=0._DP
              end if
            end do
          end if
        end do
        deallocate (dbc)

      end if

      deallocate(xsi)
      return
      end subroutine elementmatrix
