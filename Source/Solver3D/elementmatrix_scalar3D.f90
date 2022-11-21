      subroutine elementmatrix_scalar3D(elem, jacobi, full, matvar, a, b, nff, nffsum, errcode)
      use feminterface3D, only: apply_dirichlet_bc, apply_neumann_bc, tetvolume, pdecoeff3D
      use feminterface3D, only: get3Dintegpoints, shapefunctionsc3D, lam2xyz
      use femtypes
      use globalvariables3D, only : nod, vn, vp, polymaxsc, sfvertbc, edgebc, vv, ve, nnat
      implicit none
      integer (I4B) elem
      integer (I4B), allocatable :: nff(:), nffsum(:)
      integer (I4B) errcode
      complex (DPC), allocatable :: a(:,:), b(:)
      logical jacobi, full, matvar
      intent (in) :: elem, jacobi, full, matvar
      intent (out) :: a, b, nff, nffsum, errcode
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!    $Revision: 1.15 $
!    $Date: 2015/11/10 13:19:27 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Computation of the element matrix for a single element
!  using numerical integration
!
!  Input:    elem     index of the element
!            vpelem   polynomial order of the element
!                     / = .true. for the computation of the jacobian
!            jacobi
!                     \ = .false. for the linear case
!            full     =.true. compute the matrix up to vpelem
!                     = .false. only compute the sub-matrix of vpelem
!            matvar   =.true. if the material coefficients are varying across the element
!                     =.false. the material is assumed to be constant 
!
!            vn       volume element information (nodes of the element) (global)
!            nod      coordinates of tet vertices (global)
!                     nod(1,:) => x-coordinates of nodes
!                     nod(2,:) => y-coordinates of nodes
!                     nod(3,:) => z-coordinates of nodes
!
!
!  Output:   a        element matrix  or  Jacobian matrix
!            b        right hand side or residual
!            nff      size of the element matrix and rhs-vector
!                     i.e. number of form functions
!            nffsum
!            errcode  =1000 element volume is zero
!                      1001 if an integration routine of this or a higher order is not available
!                      1010 the required element order is higher than polymax (available shape functions)
!                      2001 polylo is greater than polyhi
!                      2002 polynomial degree cannot be higher than polymax
!
!  local variables
!
      integer (I4B) :: i, j, k, err, npkt, sumnff, inat, jnat
      integer (I4B) :: integorder, vpelem(nnat), maxorder, istat
      integer (I4B), allocatable :: polylo(:), polyhi(:)
      real (DP) :: vw, volume, point(3), vert(3,4)
      real (DP), allocatable :: xsi(:), gxsi(:,:)
      real (DP), allocatable :: weight(:), lambda(:,:)
      complex (DPC) :: nu(3,3,nnat,nnat), alpha(nnat,nnat), beta(3,nnat,nnat), f1(nnat), f2(3,nnat)
      complex (DPC) :: nuvw(3,3), betavw(3), gammavw(3), alphavw
      complex (DPC) :: qnu(3), qf2(3), qal, qf1
      complex (DPC), allocatable :: da(:,:)
      complex (DPC) :: gamma(3,nnat,nnat)
      logical :: isboundary, isvertexbc, isedgebc, isfacebc
!
!------------------------------------------------------------------------------------
!
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
!  Cartesian coordinates of the element vertices
      vert(1:3,1:4) = nod( 1:3,vn(1:4,elem) )
      volume = abs(tetvolume(vert))
      if (volume .le. tiny(1._DP) ) then
        volume = 0._DP
        errcode=1000
        return
      end if

!  polyhi and polylo for every nature
      allocate(polyhi(nnat), polylo(nnat), nff(nnat), nffsum(1:nnat+1), STAT=istat)
      vpelem(1:nnat)=vp(elem,1:nnat)
!
      do inat=1,nnat
        if ( full ) then
          polyhi(inat) = vpelem(inat)
          polylo(inat) = 1
!  Size of the element matrix. This would depend on which shape functions 
!  formula we use. For now we fix for Scalar element:
          nff(inat) = ( vpelem(inat) + 1 ) * ( vpelem(inat) + 2 ) * ( vpelem(inat) + 3 ) / 6
        else
          polyhi(inat) = vpelem(inat)
          polylo(inat) = vpelem(inat)
!  size of the element matrix
          if ( vpelem(inat) .gt. 1 ) then
            nff(inat) = ( vpelem(inat) + 1 )*( vpelem(inat) + 2 )/2
          else
            nff(inat) = 4
          end if
        end if
      end do
!  nffsum definition
!  example: if nff=(/4,10,4/) then
!  nffsum(1:4)=(/0,4,14,18/)
      nffsum(1)=0
      do inat=1,nnat
        nffsum(inat+1)=nffsum(inat)+nff(inat)
      end do
!
      maxorder=maxval(vpelem)
      if (maxorder .gt. polymaxsc) then
        errcode = 1010
        return
      end if
!
!----------------------------------------------------------------------------
!
      sumnff=nffsum(nnat+1)
!  Allocate required space and initialise elementmatrix a and RHS b
      allocate( a(sumnff,sumnff), da( sumnff, sumnff ), b(sumnff), STAT=istat )
      allocate( xsi(maxval(nff)), gxsi(3,maxval(nff)), STAT=istat )
!
      a(:,:)=(0._DP,0._DP)
      b(:)=(0._DP,0._DP)
!
!  fetch material coefficient at the center of gravity
      if ( .not. matvar ) then
!  if the material is constant across the element, 
!  get the centroid coordinates
        point(:) = (vert(:,1)+vert(:,2)+vert(:,3)+vert(:,4) ) /4._DP
        call pdecoeff3D( elem, point, nu, alpha, beta, f1, f2, gamma )
      end if
!
!  Determine the order of numerical integration. 
!  The weak form has at most from (xi_i*xi_j) 2*vpelem order, if we take that the 
!  coefficients do not vary much in the element; integorder of twice the vpelem of 
!  the shape functions suffice
      integorder = 2*maxorder

!  fetch numerical integration points (Gauss points)
      call get3Dintegpoints( integorder, npkt, weight, lambda, err )
      if (err .ne. 0) errcode=err
!
! start integration
int:  do k = 1, npkt
!
!  if the material is varying across the element
        if (matvar) then
!  compute the cartesian coordinates and fetch the material coefficients
          call lam2xyz( lambda(:,k), elem, point, nod, vn )
          call pdecoeff3D( elem, point, nu, alpha, beta, f1, f2, gamma)
        end if
!
!  get shape function and its gradient at location of integration point k
        call shapefunctionsc3D(lambda(:,k), vert, minval(polylo), maxval(polyhi), &
     &                         maxval(nff), .true., xsi, gxsi, errcode)
        vw = volume*weight(k)
!  loop over the natures in the direction of increasing row number
        do inat=1,nnat
!  loop over the natures in the direction of increasing column number
          do jnat=1,nnat
!  Multiply material coefficients for this nature with the volume times the weight
            nuvw = nu(:,:,inat,jnat)*vw
            alphavw = alpha(inat,jnat)*vw
            betavw = beta(:,inat,jnat)*vw
            gammavw=gamma(:,inat,jnat)*vw
!  dot product of the source vector 
            qf1 = f1(inat)*vw
!  dot product of the gradient of source vector 
            qf2 = f2(:,inat)*vw
!
!  for all shape functions
!  Since every submatrix block of the elementmatrix is symmetric we can
!  outer-loop over the column index j and inner-loop over the row index i
!  this allows for some optimisation by moving common factors out of the inner loop
!
            do j = 1,nff(jnat)
!  nu is a rank-2 tensor and its entries are multipled with
!  the gradient of xsi in x and the gradient of xsi in y
!  gamma is a vector
!  dot product of the gradient of the shape function
              qnu = matmul( nuvw, gxsi(:,j) )+ gammavw(:)*xsi(j)
!  dot product of the shape function plus the beta term above so that alpha and beta terms
!  are factored out
!  beta is a vector, alpha is a scalar
!  product of the gradient of the shape function
              qal = dot_product( gxsi(:,j), betavw ) + alphavw*xsi(j)
!
              do i=1,nff(inat)
!  integrate the function
!     i: row (degree of freedom);    j: column (test function)
!
!  note that beta is already inside qal(:)
! 
!  (nu*grad(xsi_j)+gamma*xsi_j)*grad(xsi_i)+(beta*grad(xsi_j))*xsi_i+alpha*xsi_j*xsi_i
!   ????????????????????????????????????????????    ! test sign change   ?????????????  TO DO
!  (nu*grad(xsi_j))*grad(xsi_i) - (beta * grad(xsi_j))*xsi_i - alpha*xsi_j*xsi_i
                da(i+nffsum(inat),j+nffsum(jnat)) = sum( qnu*gxsi(:,i) ) + qal*xsi(i)
              end do
!
!  right hand side contributions, integrate the function:
!             f_1*xsi_i - f_2*grad(xsi_i)
!
!             b(j) = b(j) + sum( qf1*xsi(:,j) ) - sum( qf2*gxsi(:,j) )    ! original
              if (inat .eq. jnat) then 
                b(j+nffsum(jnat)) = b(j+nffsum(jnat)) + qf1*xsi(j) - sum( qf2*gxsi(:,j) )
!   ????????????????????????????????????????????    ! test sign change   ?????????????  TO DO
!                b(j+nffsum(jnat)) = b(j+nffsum(jnat)) - qf1*xsi(j) - sum( qf2*gxsi(:,j) )
              end if
            end do
          end do
        end do
        a = a + da
      end do int
      deallocate( xsi, gxsi, da, STAT=istat)
      deallocate(weight, lambda, STAT=istat)
!---------------------------------------------------------------------------
!  contributions of the boundary conditions: prepare
!---------------------------------------------------------------------------
!
!  check if any of the faces or edges are on the boundary
!
      isvertexbc=any(sfvertbc(vn(:,elem),:) .gt. 0)
      isedgebc=any(edgebc(ve(:,elem), : ) .gt. 0)
      isfacebc=any(vv(:,elem) .lt. 0)
      isboundary=isedgebc .or. isfacebc .or. isvertexbc
!
!---------------------------------------------------------------------------
!        contributions of the boundary conditions: general or Neumann BC
!---------------------------------------------------------------------------
!
      if (isfacebc) then
        call apply_neumann_bc(elem, vert, nff, nffsum, polylo, polyhi, a, b)
      end if
!---------------------------------------------------------------------------
!                    contributions of Dirichlet BC
!---------------------------------------------------------------------------
!
      if ( isboundary ) then
        call apply_dirichlet_bc(elem, vert, nff, nffsum, polylo, polyhi, a, b)
      end if
!

      deallocate( polylo, polyhi, STAT=istat)
      return
      end subroutine elementmatrix_scalar3D



      subroutine apply_neumann_bc(elem, vert, nff, nffsum, polylo, polyhi, a, b)
      use feminterface, only: get2Dintegpoints
      use feminterface3D, only: area3D, shapefunctionsc3D, lam2xyz, getbcval
      use femtypes
      use globalvariables3D, only : nnat, vv, bctype, sbc, vn, nod, vp
      implicit none
      integer (I4B) :: elem, nff(:), nffsum(:), polylo(:), polyhi(:)
      real (DP) :: vert(3,4)
      complex (DPC) :: a(:,:), b(:)
      intent (in) :: elem, vert, nff, nffsum, polylo, polyhi
      intent (inout) :: a, b
!
!   Calculate contribution of general Boundary Conditions to element matrix and right hand side
!   by integrating over the boundary face
!
!  local variables
!
      integer (I4B) :: f, i, j, k, err, npkt2D, inat, jnat
      integer (I4B) :: istat, errcode, integorder, sfnod(3)
      real (DP) :: aw, farea, fvert(3,3), xyzt(3)
      real (DP), allocatable :: xsi(:), gxsi(:,:), lambda23(:,:)
      real (DP), allocatable :: weight2d(:), lambda2d(:,:)
      complex (DPC) :: pvalues, qvalues(nnat)
!
      integer (I4B), parameter :: f2n(3,4)=reshape((/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/))
!                                       f2n : face -> nodes
!                                              2 1 1 1
!                                              3 3 2 2
!                                              4 4 4 3
!
!------------------------------------------------------------------------------------
!
!  None of the surfaces has higher order than the element order
      allocate( xsi(maxval(nff)))
faces:  do f = 1, 4
!
        if ( vv(f,elem) .ge. 0) cycle faces
! skip if not general or Neumann (in non of the natures)
        if ( all(bctype(sbc(abs(vv(f,elem))),:) .lt. 200 .or.          &
     &           bctype(sbc(abs(vv(f,elem))),:) .ge. 300 )) cycle faces
! get the face nodes
        sfnod =  (/ vn( f2n(1,f), elem ),                              &
     &              vn( f2n(2,f), elem ),                              &
     &              vn( f2n(3,f), elem ) /)
!  get area of this face, fvert(axis,node)
        fvert(1:3,1:3) = nod( 1:3, sfnod )
        farea = area3D( fvert )
!
        do inat = 1, nnat
!  general Neumann BC is from 200 to 299
!  skip and cycle to the next nature
          if ( (bctype(sbc(abs(vv(f,elem))),inat) .lt. 200 ) .or.      &
     &         (bctype(sbc(abs(vv(f,elem))),inat) .ge. 300 )) cycle
!
!  for a boundary face, the polynomial order of the face equals the order of the volume element
!        fp(vf(f,elem),inat) = vp(elem,inat)
!  so polyhi is the same as for volume integration
!  however, here we have the integration loop inside the nature loop
!           we thus choose the integration order to depend on nature index
          integorder = 2*vp(elem,inat)
!
!  fetch integration points on this triangle face
          call get2Dintegpoints( integorder, npkt2D, weight2d, lambda2d, err )
!  match lambda2d to lambda
          allocate( lambda23(4,npkt2D), STAT=istat )
          lambda23(f,:) = 0._DP
          do i = 1,3
            lambda23( f2n(i,f),: ) = lambda2d( i,: ) 
          end do
!
!---------------------------------------------------------------------------
!                Integration of the shape functions and the RHS
!---------------------------------------------------------------------------
integ:    do k = 1, npkt2D
            aw = farea * weight2d(k)
!
!  get shape functions at this location (valid for all natures)
            call shapefunctionsc3D(lambda23(:,k), vert, minval(polylo), maxval(polyhi), &
     &                             maxval(nff), .false., xsi, gxsi, errcode)
! fetch the coefficients of the general bc.
            call lam2xyz( lambda23(:,k), elem, xyzt, nod, vn )
!---------------------------------------------------------------------------
! for scalar problem, pvalues and qvalues are scalars
!---------------------------------------------------------------------------
            call getbcval( sbc( abs(vv(f,elem))), inat , xyzt, pvalues, qvalues )
!---------------------------------------------------------------------------
! NOTE: the basis funtions on other faces are zero on this face, but we 
!       include them nevertheless for simplicity.
!       might modify, if it speeds up the program. 
!---------------------------------------------------------------------------
            do i = 1, nff(inat)
              if (xsi(i) .eq. 0._DP) cycle
              do jnat=1,nnat
                do j = 1, nff(jnat)
!  integrate the function: ( q * xsi(j) * xsi(i) )
                  a(i+nffsum(inat),j+nffsum(jnat)) = a(i+nffsum(inat),j+nffsum(jnat)) - aw*qvalues(jnat)*xsi(i)*xsi(j)
!  TO DO for the moment we assume  qvalues to only have one index, in general it should have two: qvalues(inat,jnat)
                end do
              end do
!---------------------------------------------------------------------------
! right hand side contributions, integrate the function: (p * xsi_i)
! for scalar problem, pvalues is scalar-valued
!---------------------------------------------------------------------------
              b(i+nffsum(inat)) = b(i+nffsum(inat)) +  aw*xsi(i)*pvalues
            end do
          end do integ
          deallocate(lambda2d, lambda23, weight2d, STAT=istat)
        end do  ! nature loop
      end do faces
      deallocate( xsi,STAT=istat)
      return
      end subroutine



      subroutine apply_dirichlet_bc(elem, vert, nff, nffsum, polylo, polyhi, a, b)
      use feminterface, only: get1Dintpolpoints, lusolver, get2Dintegpoints
      use feminterface3D, only: shapefunctionsc3D, lam2xyz
      use feminterface3D, only: getbcval, area3D
      use femtypes
      use globalvariables3D, only : sfvertbc, bctype, vv, vn, ve, ep, edgebc, nod, vp, sbc, nnat
      implicit none
      integer (I4B) :: elem, nff(:), nffsum(:), polylo(:), polyhi(:)
      real (DP) :: vert(3,4)
      complex (DPC) :: a(:,:), b(:)
      intent (in) :: elem, vert, nff, nffsum, polylo, polyhi
      intent (inout) :: a, b
!
!  local variables
!
      integer (I4B) :: f, i, j, k, err, npkt2D, ivertex, sumnff, maxorder
      integer (I4B) :: edgenum, lnj, nffunix, inat, errcode
      integer (I4B) :: integorder, intpolpoints, sfnod(3), nr, nl, istat
      integer (I4B) :: nffborder, nffface, p, startlni, stoplni
      integer (I4B), allocatable :: startedge(:), startface(:), startvol(:)
      integer (I4B), allocatable :: numfacedof(:), numvoldof(:)
      integer (I4B), allocatable :: borderidx(:), faceidx(:)
      real (DP) :: aw, awxsi, farea, fvert(3,3), xyzt(3), lam(4)
      real (DP), allocatable :: xsi(:), gxsi(:,:), xtab(:), lambda23(:,:)
      real (DP), allocatable :: weight2d(:), lambda2d(:,:)
      complex (DPC) :: pvalues
      complex (DPC), allocatable :: mat(:,:), rhs(:)
      logical isvertexbc, isedgebc, isfacebc
      logical, allocatable :: dbc(:)
!
      integer (I4B), parameter :: f2n(3,4)=reshape((/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/))
      integer (I4B), parameter :: f2e(3,4)=reshape((/2,5,6,3,4,6,1,4,5,1,2,3/),(/3,4/))
      integer (I4B), parameter :: e2lambda(4,6)=reshape((/1,2,3,4,2,3,4,1,1,3,4,2, &
     &                                                    1,4,2,3,2,4,1,3,3,4,1,2/),(/4,6/))
!  f2e : face -> edges                  f2n : face -> nodes
!           2 3 1 1                            2 1 1 1
!           5 4 4 2                            3 3 2 2
!           6 6 5 3                            4 4 4 3
!
!  e2lambda : edge -> lambda. The first 2 are for lambda's on edge i, next 2 to be 0 
!           1 2 1 1 2 3 
!           2 3 3 4 4 4
!           3 4 4 2 1 1
!           4 1 2 3 3 2
!
!------------------------------------------------------------------------------------
!
      maxorder=maxval(polyhi)
      allocate( numfacedof(maxorder), numvoldof(maxorder), STAT=istat )
!
!  number of DOF's of a face and vol at element's order p (not including
!  p-1, p-2, ... )
      numfacedof(1) = 0
      numvoldof(1)  = 0
      do p = 2,maxorder
        numfacedof(p) = p-2
        numvoldof(p)  = (p-2)*(p-3)/2
      end do
!
!----------------------------------------------------------------------------
!                   setup for accessing local dof indices
!----------------------------------------------------------------------------
!
      allocate( startedge(maxorder), startface(maxorder), startvol(maxorder), STAT=istat )
!  create starting indices for local dofs' indexing
!  these are the indices +1 at which edge, face or vol dofs appear at degree p
      startedge(1) = 0
      startface(1) = 0
      startvol(1)  = 4
! sorting of DOFs (shape functions) is:
! number    1..4 | 5..10|11..16|17..20|21..26|27..34| 35 |36..41|42..53|54..56|...
! type     vertex| edge | edge | face | edge | face |vol | edge | face |  vol |...
! order       1  |   2  |   3  |   3  |   4  |   4  |  4 |   5  |   5  |   5  |...
      do p = 2,maxorder
        startedge(p) = startvol(p-1) + numvoldof(p-1)
        startface(p) = startedge(p)  + 6
        startvol(p)  = startface(p)  + numfacedof(p)*4
      end do
!----------------------------------------------------------------------------
      allocate( xsi(maxval(nff)))
      sumnff=nffsum(nnat+1)
      allocate( dbc(sumnff), STAT=istat )
!  dbc(i) will become .true. if the DOF i belongs to a Dirichlet BC
      dbc = .false.
!
!---------------------------------------------------------------------------
! do Dirichlet BC at vertices
!---------------------------------------------------------------------------
!  visit the edges to find which vertices that are attached to the boundary edges
      isvertexbc=any(sfvertbc(vn(:,elem),:) .gt. 0)
! To avoid error with oneAPI. dummy size
      allocate(gxsi(3, 10))
      if (isvertexbc) then
        do inat=1,nnat
          do ivertex = 1, 4
            if (sfvertbc(vn(ivertex,elem),inat) .gt. 0) then
!  Dirichlet BC is from 0 to 99
              if (bctype(sfvertbc(vn(ivertex,elem),inat),inat) .ge. 100) cycle
              lam(:)=0._DP
              lam(ivertex)=1._DP
!  get vertex (order 1) shape functions
              call shapefunctionsc3D(lam, vert, 1, 1, 4, .false., xsi, gxsi, errcode)
              xyzt=vert(:,ivertex)
              call getbcval(sfvertbc(vn(ivertex,elem),inat), inat, xyzt, pvalues)
              b( ivertex + nffsum(inat) ) = pvalues / xsi(ivertex)
              dbc( ivertex + nffsum(inat) ) = .true.
            end if
          end do
        end do ! end inat
      end if
      deallocate(gxsi)

!---------------------------------------------------------------------------
! do Dirichlet BC at edges
!---------------------------------------------------------------------------
      isedgebc=any(edgebc(ve(:,elem), : ) .gt. 0)
      if (isedgebc) then
        do inat=1,nnat
edge_d:   do k = 1, 6
            edgenum = ve(k,elem)
!  Resume for specifying Dirichlet on edges if element order is more than 1
            if (ep( edgenum, inat ) .le. 1) cycle
            if (edgebc(edgenum, inat ) .le. 0) cycle
!  skip if not Dirichlet
!  Dirichlet BC is from 0 to 99
            if (bctype(edgebc(edgenum,inat),inat) .ge. 100) cycle
!  get interpolation points for this edge
!  For a scalar tetrahedron, we follow the implementation
!  as in scalar triangle for the edge BC
!  order p has p+1 Gauss Lobatto interpolation points
            intpolpoints = ep( edgenum, inat ) +1
            allocate(xtab( intpolpoints), STAT=istat)
! Dummy allocate to avoid runtime error in oneAPI
            allocate( gxsi(3,10) )
! the no. of points for BC is reduced by the 2 vertex dofs
            allocate( mat( intpolpoints-2, intpolpoints-2 ), rhs( intpolpoints-2 ), STAT=istat )
!  fetch interpolation points along this edge, given by Gauss Lobatto with range [-1,1]
            call get1Dintpolpoints( intpolpoints, xtab, errcode )
!
!  do up to G-L points, skipping the 2 vertex points
            nr= e2lambda(1,k)
            nl= e2lambda(2,k)
interpol:   do i = 2, intpolpoints-1
!  convert the 1D xtab to lambdas along this edge
              lam( nl ) = 0.5_DP + xtab(i)/2._DP
              lam( nr ) = 0.5_DP - xtab(i)/2._DP
              lam( e2lambda(3,k) ) = 0._DP
              lam( e2lambda(4,k) ) = 0._DP
!
              call shapefunctionsc3D( lam, vert, polylo(inat), polyhi(inat),   &
     &                                nff(inat),.false., xsi, gxsi, errcode )
!
!  do up to edge order, skipping the 2 vertex points
              do j = 2, ep( edgenum, inat )
!  local dof-number on edge k
                lnj = startedge(j) + k
                mat(i-1,j-1) = xsi(lnj)
              end do
!
!  setup RHS, get BC value of the edge
              call lam2xyz( lam, elem, xyzt, nod, vn )
              call getbcval( edgebc( edgenum, inat ), inat, xyzt, pvalues ) 
              rhs(i-1) = pvalues - xsi(nr)*b(nr+nffsum(inat)) - xsi(nl)*b(nl+nffsum(inat))
            end do interpol
            deallocate(gxsi)

!  solve for x
            call lusolver( mat, rhs )
!  assign to b.
!  do up to the edge polynomial degree
            do j = 2, ep( edgenum, inat )
!  update lnj with the approriate "jump" for the next local dof on the edge
              lnj = startedge(j) + k
              b(lnj+nffsum(inat)) = rhs(j-1)
              dbc(lnj+nffsum(inat)) = .true.
            end do
            deallocate( xtab, mat, rhs, STAT=istat )
          end do edge_d
        end do ! end inat
      end if
!
!---------------------------------------------------------------------------
!         Fulfiling  Dirichlet BC for face 
!---------------------------------------------------------------------------
!
!  integrate shape functions to solve for the dofs at (iiner) face, 
!  the order of shape functions on this face equals the order of the volume
!
!  None of the surfaces has higher order than the element order, so 
!  only do for element order greater than 1
      isfacebc=any(vv(:,elem) .lt. 0)
      if ( isfacebc ) then
        do inat=1,nnat
          if ( vp(elem,inat) .le. 2 ) cycle
faces:    do f = 1, 4
            if ( vv(f,elem) .gt. 0) cycle
!  skip if not Dirichlet
!  Dirichlet BC is from 0 to 99
            if ( bctype(sbc( abs(vv(f,elem))),inat ) .ge. 100 ) cycle
!
            sfnod =  (/ vn( f2n(1,f), elem ),                           &
     &                  vn( f2n(2,f), elem ),                           &
     &                  vn( f2n(3,f), elem ) /)
!  get area of this face, fvert(axis,node)
            fvert(1:3,1:3) = nod( 1:3, sfnod )
            farea = area3D( fvert )
!  for a boundary face, the polynomial order of the face equals the order of the volume element
!        fp(vf(f,elem),inat) = vp(elem,inat)
!  so polyhi is the same as for volume integration
!  however, here we have the integration loop inside the nature loop
            integorder = 2*vp(elem,inat)
!  fetch integration points on this triangle face
            call get2Dintegpoints( integorder, npkt2D, weight2d, lambda2d, err )
!  match lambda2d to lambda
            allocate( lambda23(4,npkt2D), STAT=istat )
            lambda23(f,:) = 0._DP
            do i = 1,3
              lambda23( f2n(i,f),: ) = lambda2d( i,: )
            end do
!---------------------------------------------------------------------------
!                Setup for integration 
!---------------------------------------------------------------------------
!
!  number nffborder of vertex+edge functions equals the order ep of the edges
!  access the edge order from its global edge name with the local face-edge map f2e
            nffborder = ep( ve( f2e(1,f), elem ), inat ) +              &
     &                  ep( ve( f2e(2,f), elem ), inat ) +              &
     &                  ep( ve( f2e(3,f), elem ), inat )
!  number nffface of (inner) face functions
            nffface = sum( numfacedof(1:vp(elem,inat) ) )
!
!  build local indexing by copying the local indices of shape functions on this face
!  -  for the face border (vertices, edges)
!  -  and for face (inner face shape functions)
            allocate( borderidx( nffborder ), faceidx(nffface), STAT=istat )
!    copy vertex dofs first
            borderidx(1:3) = f2n(:,f)
!    copy edge dofs, assuming each edge has its own order
            nffunix = 3
            do i = 1, 3
              do p = 2, ep( ve(f2e(i,f),elem), inat )
                nffunix = nffunix + 1
                borderidx(nffunix) = startedge(p) + f2e(i,f)
              end do
            end do
!    copy the face indices
            nffunix = 1
            do p = 3, vp(elem,inat)
              startlni = startface(p) + numfacedof(p)*(f - 1) + 1
              stoplni = startlni + numfacedof(p) - 1
              faceidx( nffunix: nffunix + stoplni - startlni) = (/(i, i=startlni,stoplni)/)
              nffunix = nffunix + stoplni - startlni + 1
            end do
!
            allocate( mat(nffface, nffface), rhs(nffface), STAT=istat )
            mat = 0._DP
            rhs = 0._DP
!
!---------------------------------------------------------------------------
!                Integration of the shape functions and the RHS
!                  for the submatrix corresponding to the face
!---------------------------------------------------------------------------
!
integrate:  do k = 1, npkt2D
              aw = farea * weight2d(k)
!
!  get shape functions at this location.
!  for a boundary face, the order of face functions equals the order of the volume element,
!  polylo and polyhi are used for this nature
!  Dummy allocate to avoid runtime error in oneAPI
              allocate( gxsi(3,10) )
              call shapefunctionsc3D( lambda23(:,k), vert, polylo(inat), polyhi(inat),   &
     &                                nff(inat),.false., xsi, gxsi, errcode )
              deallocate(gxsi)
     !
!  calculate the scalar product  <xsi(i) , xsi(j)>  by integration
              do j = 1, nffface
                awxsi = aw * xsi(faceidx(j))
                do i = 1, nffface
                  mat(i,j) = mat(i,j)+ xsi(faceidx(i))*awxsi
                end do
!
!  integrate rhs on this face. 
                call lam2xyz( lambda23(:,k), elem, xyzt, nod, vn )
                call getbcval( sbc( abs(vv(f,elem))), inat, xyzt, pvalues )
!
!  calculate:  <pval , xsi(j)> - sum(<B(i)*xsi(i), xsi(j)>, over i)
                rhs(j) = rhs(j) + pvalues*awxsi
                do i = 1, nffborder
                  rhs(j) = rhs(j) - b(borderidx(i)+nffsum(inat))*xsi(borderidx(i))*awxsi
                end do
              end do
            end do integrate
!  solve for x (x is returned to rhs vector)
            call lusolver(mat, rhs)
!  fix the face dofs only, the vertex and edge dofs are fixed already
            do i = 1, nffface
              b( faceidx(i)+nffsum(inat) ) = rhs(i)
              dbc( faceidx(i)+nffsum(inat) ) = .true.
            end do
!
            deallocate( lambda23, lambda2d, weight2d, mat, rhs, borderidx, faceidx, STAT=istat )
          end do faces
        end do
      end if
      deallocate (startedge, startface, startvol, numfacedof, numvoldof)
      deallocate (xsi)
!
!-------------------------------------------------------------------
!    modify the element matrix to enforce essential (Dirichlet) BC
!-------------------------------------------------------------------
      do i = 1, sumnff
ifdbc:  if ( dbc(i) ) then
!  this is the entry of submatrix A_DD (ref: Matrix Assembly.doc)
!  for Dirichlet dof, let all the columns zero first, then 
!  set the diagonal to unity, so A_DD = I.
!  the subvector b_D is already fixed above by assigning b with the 
!  correct local indexing: b(lni)
          a(i,1:sumnff) = 0._DP
          a(i,i) = 1._DP
        else
!  this entry has general Neumann, so nullify submatrix A_DN and A_ND;
!  and modify b_N 
          do j = 1, sumnff
            if ( dbc(j) ) then
              b(i) = b(i) - a(i,j)*b(j)
              a(i,j) = 0._DP
            end if
          end do
        end if ifdbc
      end do
      deallocate(dbc, STAT=istat)
      return
      end


