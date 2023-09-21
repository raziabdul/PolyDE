      subroutine elementmatrix3D(elem, vpelem, jacobi, full, matvar, a, b, nff, errcode)
      use feminterface3D, only: get3Dintegpoints
      use feminterface3D, only: pdecoeff3D, shapefunctionv3D, lam2xyz, xyz2lam
      use feminterface3D, only: tetvolume, bcstovv, area3D, getbcval_vec
      use feminterface3D, only: getgl
      use feminterface, only: lusolver, get1Dintpolpoints, get2Dintegpoints, cross_product
      use femtypes
      use globalvariables3D, only: sbc, en, edgebc, nod, polymaxv, ve, vn, vp, vv, bctype, ep
      implicit none
      integer (I4B) elem
      integer (I4B) vpelem, nff, errcode
      complex (DPC), allocatable :: a(:,:), b(:)
      logical jacobi, full, matvar
      intent (in) :: elem, vpelem, jacobi, full, matvar
      intent (out) :: nff, a, b
!
!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 - 2015 Institute for Micro Systems Technology,
!                              Hamburg University of Technology.
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
!
!    $Revision: 1.35 $
!    $Date: 2015/10/27 10:22:09 $
!    $Author: m_kasper $
!
!
!  Computation of the element matrix for a single element
!  using numerical integration
!
!  NOTES: This is the version using matmul which works best to date... [Razi]
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
!            errcode  =1000 triangle area is zero or negative
!
!  local variables
!
      complex (DPC), allocatable :: mat(:,:), rhs(:), da(:,:)
      complex (DPC) :: nu(3,3), alpha(3,3), beta(3), f1(3), f2(3), awxsi(3)
      complex (DPC) :: nuvw(3,3), alphavw(3,3), betavw(3), pvalues(3), qvalues(3,3)
      complex (DPC) :: qbe(3), qnu(3), qal(3), qf1(3), qf2(3) 

      real (DP), allocatable :: weight3d(:), lambda3d(:,:)
      real (DP), allocatable :: xtab(:), weight2d(:), lambda2d(:,:)
      real (DP), allocatable :: xsi(:,:), cxsi(:,:), lambda(:,:)
      real (DP) :: aw, vw
      real (DP) :: farea, volume, point(3), vert(3,4), fvert(3,3), xyzt(3), lam(4)
      real (DP) :: gl(3,4), nvec(3,4), t1(3,4), t2(3,4), tvec(3)

      logical, allocatable :: dbc(:)
      logical :: isboundary

      integer (I4B), allocatable :: startedge(:), startface(:), startvol(:)
      integer (I4B), allocatable :: lnarray(:) 
      integer (I4B) :: f, i, j, k, polylo, polyhi, err, npkt, npkt2D
      integer (I4B) :: edgenum, lnj, nffunix
      integer (I4B) :: integorder, intpolpoints, sfnod(3)
      integer (I4B) :: nffe, nffun, p, startlni, stoplni 

      integer (I4B), parameter :: f2n(3,4)=reshape((/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/))
      integer (I4B), parameter :: f2e(3,4)=reshape((/2,5,6,3,4,6,1,4,5,1,2,3/),(/3,4/))
      integer (I4B), parameter :: e2lambda(4,6)=reshape((/1,2,3,4,2,3,4,1,1,3,4,2, &
     &                                                    1,4,2,3,2,4,1,3,3,4,1,2/),(/4,6/))
!  number of DOF's of a face and vol at element's order p (not including
!  p-1, p-2, ... )
      integer (I4B), parameter :: numfacedof(polymaxv) = (/ ( 2*(p-1), p = 1,polymaxv) /)
      integer (I4B), parameter :: numvoldof(polymaxv)  = (/ ( 3*(p-1)*(p-2)/2, p = 1,polymaxv ) /)
!  f23 : face -> edges                  f2n : face -> nodes      
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

!  There are two possible variants in numerical integration, here we use:
!       the same integration rule for all entries in the matrix,
!       which offers the advantage of evaluating all formfunctions at 
!       the same locations thus reducing the work required in the evaluation of formfunction
!
      errcode=0
!
!----------------------------------------------------------------------------
!                   setup for accessing local dof indices
!----------------------------------------------------------------------------
      allocate( startedge(polymaxv), startface(polymaxv),               &
     &          startvol(polymaxv) )

!  create starting indices for local dofs' indexing
      startedge(1) = 0
      startface(1) = -1
      startvol(1)  = 6

      do p = 2,polymaxv
        startedge(p) = startvol(p-1) + numvoldof(p-1)
        startface(p) = startedge(p)  + 6
        startvol(p)  = startface(p)  + numfacedof(p)*4
      end do
      nff=0
!----------------------------------------------------------------------------
!  Cartesian coordinates of the element vertices
!  this line is auto-parallelized
      vert(1:3,1:4) = nod( 1:3,vn(1:4,elem) )
      volume = abs(tetvolume(vert))
      if (volume .le. tiny(1._DP) ) then
        volume = 0._DP
        errcode=1000
        return
      end if
!
      if ( full ) then
        polyhi = vpelem
        polylo = 1
!  Size of the element matrix. This would depend on which shape functions 
!  formula we use. For now we fix for Nedelec element:
        nff = vpelem * ( vpelem + 2 ) * ( vpelem + 3 ) / 2
      else
        polyhi = vpelem
        polylo = vpelem
!  size of the element matrix
        if ( vpelem .gt. 1 ) then
          nff = ( vpelem + 2 )*( 3*vpelem + 1 )/2
        else
          nff = 6
        end if
      end if
!
      allocate( a(nff,nff), b(nff), xsi(3,nff), cxsi(3,nff), da( nff, nff ) )
      a(:,:)=(0._DP,0._DP)
      b(:)=(0._DP,0._DP)
!
!  fetch material coefficient at the center of gravity
      if ( .not. matvar ) then
!  if the material is constant across the element, 
!  get the centroid coordinates
!  the following 3 points are auto-parallelized
        point(1) = sum( vert(1,:) )
        point(2) = sum( vert(2,:) )
        point(3) = sum( vert(3,:) )
        point = point/4._DP
        call pdecoeff3D( elem, point, nu, alpha, beta, f1, f2 )
      end if
!
!  Determine the order of numerical integration. 
!  The weak form has at most from (xi_i*xi_j) 2*vpelem order, if we take that the 
!  coefficients do not vary much in the element; integorder of twice the vpelem of 
!  the shape functions should suffice
      integorder = 2*vpelem

!  fetch numerical integration points (Gauss points)
      call get3Dintegpoints( integorder, npkt, weight3d, lambda3d, err )
      if (err .ne. 0) errcode=err
!
!  NOTES about integration:
!   - remember here that xsi(:,j) is  a vector x,y, OR z, so 
!     xsi and cxsi are treated in the same way
!   - matmul and dot_product seems to give convergence problem for the 
!     SSORCG. operations are reverted back to scalar ones as in 2D

!  fetch normal and tangential vectors of the faces
      call getgl( vert, gl, nvec, t1, t2 )

! start integration. See NOTES above about vector and matrix operations
! NOTE: this loop makes the whole program slow
int:  do k = 1, npkt

!  if the material is varying across the element
        if (matvar) then
!  compute the cartesian coordinates and fetch the material coefficients
          call lam2xyz( lambda3d(:,k), elem, point, nod, vn )
          call pdecoeff3D( elem, point, nu, alpha, beta, f1, f2 )
        end if
!
!  get shape function and its curls at location of integration point k
!        call shapefunctionv3D(lambda3d(:,k), vert, polylo, polyhi, nff, .true., xsi, cxsi, errcode)
        call shapefunctionv3D(lambda3d(:,k), vert, polylo, polyhi, .true., xsi, cxsi, errcode)

        !
        vw = volume*weight3d(k)
        nuvw = nu*vw
        alphavw = alpha*vw
        betavw = beta*vw
!  dot product of the source vector 
        qf1 = f1*vw
!  dot product of the curl of source vector 
        qf2 = f2*vw
!
!  for all shape functions
! to vectorize from here ..
        do j = 1,nff
!  dot product of the curl of the shape function
          qnu = matmul( nuvw, cxsi(:,j) )         ! unvectorizable: data type not supported??
!  cross product of the curl of the shape function
          qbe = cross_product( cxsi(:,j), betavw )! unvectorizable: data type not supported
!  dot product of the shape function plus the beta term above so that alpha and beta terms
!  are factored out
          qal = matmul( alphavw, xsi(:,j) ) + qbe ! unvectorizable: data type not supported
          do i=1,nff
!  integrate the function
!     i: row (degree of freedom);    j: column (test function)
!
!  note that beta is already inside qal(:)
! 
!  (nu*curl(xsi_j))*curl(xsi_i) - (beta x curl(xsi_j))*xsi_i - alpha*xsi_j*xsi_i
!
!  REMARK: This operation is very expensive, takes >50% of matrix assembly time
            da(i,j) = sum( qnu*cxsi(:,i) ) - sum( qal*xsi(:,i) )     ! test sign change--right per manfred
          end do
          
!  right hand side contributions, integrate the function:
!             f_1*xsi_i - f_2*curl(xsi_i)
!
!          b(j) = b(j) + sum( qf1*xsi(:,j) ) - sum( qf2*cxsi(:,j) )    ! original
          b(j) = b(j) - sum( qf1*xsi(:,j) ) - sum( qf2*cxsi(:,j) )    ! change sign--looks correct  [Razi]
       end do
        a = a + da
      end do int

      deallocate(weight3d, lambda3d, da)

!---------------------------------------------------------------------------
!  contributions of the boundary conditions: prepare
!---------------------------------------------------------------------------

!  check if any of the faces or edges are on the boundary
!
      isboundary=any(edgebc(ve(:,elem ),1) .gt. 0) .or. any(vv(:,elem) .lt. 0)
!

!---------------------------------------------------------------------------
!        contributions of the boundary conditions: general of Neumann BC
!---------------------------------------------------------------------------
bc_n: if (isboundary) then
!  None of the surfaces has higher order than the element order, so 
!  only do for element order greater than 1
        if ( vp(elem,1) > 1 ) then
faceloop1:  do f = 1, 4
            if ( vv(f,elem) < 0) then
! skip if not general of Neumann
! there's only one nature for now 
              if ( bctype(sbc( abs(vv(f,elem))),1 ) .lt. 200 ) cycle
              sfnod =  (/ vn( f2n(1,f), elem ),                        &
     &                    vn( f2n(2,f), elem ),                         &
     &                    vn( f2n(3,f), elem ) /)

!  get area of this face, fvert(axis,node)
              fvert(1:3,1:3) = nod( 1:3, sfnod )
              farea = area3D( fvert )
              integorder = 2*vp(elem,1)
!  fetch integration points on this triangle face
              call get2Dintegpoints( integorder, npkt2D, weight2d, lambda2d, err )
!  match lambda2d to lambda
              allocate( lambda(4,npkt2D) )
              lambda(f,:) = 0._DP
              do i = 1,3
                lambda( f2n(i,f),: ) = lambda2d( i,: ) 
              end do

!---------------------------------------------------------------------------
!                Integration of the shape functions and the RHS
!
! NOTE: the basis funtions on other faces are zero on this face, but we 
!       include them nevertheless for simplicity.  might modify if it speeds 
!       up the program. 
!---------------------------------------------------------------------------
!
integ:        do k = 1, npkt2D
                aw = farea * weight2d(k)
!
!  get shape functions at this location.
!  for a boundary face, the order of face functions is the order of the volume element,
!  so polylo and polyhi don't change from element integration above
                call shapefunctionv3D( lambda(:,k), vert, polylo, polyhi, .false.,  &
     &                                 xsi, cxsi, errcode )
!
! fetch the coefficients of the general bc.
                call lam2xyz( lambda(:,k), elem, xyzt, nod, vn )
! there's only one nature for now so
                call getbcval_vec( sbc( abs(vv(f,elem))), 1 , xyzt, pvalues, qvalues )
! for the pvalues, we are only concerned with the tangential components on this face 
                pvalues = cross_product( nvec(:,f), pvalues )
! for the qvalues, *it seems* that we are also only concerned with the tangential components on this face
! but this is unclear with the ABC, since as in 2D we want the normal to be specified
! should the qval tensor be modified?
! test this first:  
                do i = 1, nff
                  awxsi = aw * cross_product( nvec(:,f), matmul( qvalues, xsi(:,i)) )
!                  if ( bctype(sbc( abs(vv(f,elem)) ) .gt. 200 ) print *, qvalues
                  do j = 1, nff
!  integrate the function: ( q*xsi_i x xsi_j )_n
!  which is equivalent to: xsi_j * ( n x q*xsi_i )
!  NOTE: avoid complex argument of A in dot_product(A,B) since the complex conjugate  
!  of A is used 
                    a(i,j) = a(i,j) - dot_product( xsi(:,j), awxsi  )
!testing:                      a(i,j) = a(i,j) - dot_product( xsi(:,j), aimag(awxsi)  )
                  end do
!
!  right hand side contributions, integrate the function: (p x xsi_i)_n
!  which is equivalent to: xsi_i * ( n x p )
!  pvalues here are already the tangential components 
                  b(i) = b(i) +  aw*dot_product( xsi(:,i), pvalues ) 
                end do
              end do integ
              deallocate(lambda2d, weight2d, lambda)
            end if
          end do faceloop1
        end if
      end if bc_n

!---------------------------------------------------------------------------
!                    contributions of Dirichlet BC
!---------------------------------------------------------------------------

!  do Dirichlet BC on edges
bc:   if (isboundary) then
        allocate( dbc(nff) )
        dbc = .false.
!  visit the edges
edge:   do k = 1,6
!  negative edgebc i of volume elem means it has a no BC
          edgenum = ve( k, elem )
! dirichlet numbers run from 0 to 199 
has_dbc:  if ( edgebc( edgenum, 1 ) .gt. 0 ) then 
! dirichlet numbers run from 0 to 199, so skip if greater (only a single nature)
            if ( bctype(edgebc( edgenum, 1 ),1) .ge. 200) cycle
!  get interpolation points for this edge
            intpolpoints = ep( ve( k, elem ), 1 )
            allocate( xtab( intpolpoints ), mat( intpolpoints, intpolpoints ),  &
     &                rhs( intpolpoints ) )
!  get tangential direction along edge
!              edgelength = sqrt( sum( ( nod( 1:3, en(2,edgenum) ) - nod( 1:3, en(1,edgenum) ) )**2 ) )
!  the direction vector over its length gives a unit vector of this edge
            tvec(1:3) = nod( 1:3, en(2,edgenum) ) - nod( 1:3, en(1,edgenum) )
!  fetch interpolation points along this edge, given by Gauss Lobatto with range [-1,1]
            call get1Dintpolpoints( intpolpoints, xtab, errcode )
interpolate:do i = 1, intpolpoints
!  convert the 1D xtab to lambdas along this edge
              lam( e2lambda(2,k) ) = ( xtab(i) + 1._DP )/2._DP 
              lam( e2lambda(1,k) ) = 1._DP - lam( e2lambda(2,k) )
              lam( e2lambda(3,k) ) = 0._DP
              lam( e2lambda(4,k) ) = 0._DP
!
              call shapefunctionv3D(lam, vert, polylo, polyhi, .false., xsi, cxsi, errcode)
!  For Nedelec's element:
!    #  of functions on an edge = the order of the edge
!  this is then equal the number of interpolation points, so here the indexing is done 
!  with the edge order j
!
!  do up to edge's order
              do j = 1, ep( edgenum , 1 )
!  local dof-number on edge k
                lnj = startedge(j) + k
                mat(i,j) = dot_product( xsi(:,lnj), tvec )
              end do

!  setup RHS, get BC value of the edge (only a single nature)
              call lam2xyz( lam, elem, xyzt, nod, vn )
              call getbcval_vec( edgebc( edgenum, 1 ), 1, xyzt, pvalues )  ! Dirichlet BC only
              rhs(i) = dot_product( tvec, pvalues )
            end do interpolate
!  solve for x
            call lusolver( mat, rhs )
!  assign to b.
!  do up to the edge polynomial degree (which is the same as the number
!  of interpolation points for the Nedelec elements)
            do j = 1, ep( edgenum, 1 )
!  update lnj with the approriate "jump" for the next local dof on the edge
              lnj = startedge(j) + k
              b(lnj) = rhs(j)
              dbc(lnj) = .true.
            end do
            deallocate( xtab, mat, rhs )
          end if has_dbc
        end do edge

!---------------------------------------------------------------------------
!         Fulfiling  Dirichlet BC for face 
!---------------------------------------------------------------------------

!  integrate shape functions to solve for the dofs
!  for a boundary face, the order of shape functions 
!  on this face are the order of the volume

!  None of the surfaces has higher order than the element order, so 
!  only do for element order greater than 1
        if ( vp(elem,1) > 1 ) then
faceloop: do f = 1, 4
            if ( vv(f,elem) < 0) then
! skip if not Dirichlet
! there's only one nature for now 
            if ( bctype(sbc( abs(vv(f,elem))),1 ) .ge. 200 ) cycle
              sfnod =  (/ vn( f2n(1,f), elem ),                         &
     &                    vn( f2n(2,f), elem ),                         &
     &                    vn( f2n(3,f), elem ) /)

!  get area of this face, fvert(axis,node)
              fvert(1:3,1:3) = nod( 1:3, sfnod )

              farea = area3D( fvert )
              integorder = 2*vp(elem,1)
!  fetch integration points on this triangle face
              call get2Dintegpoints( integorder, npkt2D, weight2d, lambda2d, err )
!  match lambda2d to lambda
              allocate( lambda(4,npkt2D) )
              lambda(f,:) = 0._DP
              do i = 1,3
                lambda( f2n(i,f),: ) = lambda2d( i,: ) 
              end do
!---------------------------------------------------------------------------
!                Setup for integration 
!---------------------------------------------------------------------------
!  
!  number of edge functions are the order of the edges. access the edge order 
!  from its global edge name with the local face-edge map f2e
              nffe = ep( ve( f2e(1,f), elem ), 1 ) +                    &
     &               ep( ve( f2e(2,f), elem ), 1 ) +                    &
     &               ep( ve( f2e(3,f), elem ), 1 )
              nffun = sum( numfacedof(1:vp(elem,1) ) ) + nffe

!  build local indexing by copying the local indices of edges and face functions
!  on this face 
              nffunix = 0
              allocate( lnarray( nffun ) )
!    copy edge dofs first, assuming each edge has its own order
              do i = 1, 3
                do p = 1, ep( ve(f2e(i,f),elem), 1 )
                  nffunix = nffunix + 1
                  lnarray(nffunix) = startedge(p) + f2e(i,f)
                end do
              end do
!    copy the face
              nffunix = nffunix + 1
              do p = 2, vp(elem,1)
                startlni = startface(p) + numfacedof(p)*(f - 1) + 1
                stoplni = startlni + numfacedof(p) - 1
                lnarray( nffunix: nffunix + stoplni - startlni) = (/startlni:stoplni/)
                nffunix = nffunix + stoplni - startlni + 1
              end do
!             nffun = numfacedof( vp(elem) )
              allocate( mat(nffun, nffun), rhs(nffun) )
              mat = 0._DP
              rhs = 0._DP

!---------------------------------------------------------------------------
!                Integration of the shape functions and the RHS
!                  for the submatrix corresponding to the face  
!---------------------------------------------------------------------------
!
integrate:    do k = 1, npkt2D
                aw = farea * weight2d(k)
!
!  get shape functions at this location.
!  for a boundary face, the order of face functions is the order of the volume element,
!  so polylo and polyhi don't change from element integration above
                call shapefunctionv3D( lambda(:,k), vert, polylo, polyhi, .false.,  &
&                                      xsi, cxsi, errcode )
!
!  get local dof of the face function on face f
                do j = 1, nffun
!  take the tangential components of xsi to this face by taking the cross
!  product of the face's normal with xsi
                  awxsi = aw * cross_product( nvec(:,f), xsi(1:3,lnarray(j)) )
                  do i = 1, nffun
                    mat(i,j) = mat(i,j)                                                      &
     &                       + dot_product( cross_product(nvec(:,f),xsi(1:3,lnarray(i))), awxsi )
                  end do
!
!  integrate rhs on this face. 
                  call lam2xyz( lambda(:,k), elem, xyzt, nod, vn )
! there's only one nature for now 
                  call getbcval_vec( sbc( abs(vv(f,elem))), 1 , xyzt, pvalues )
!  use sum of vector*vector to avoid complex conjugate of the 1st arg in dot_product
                  rhs(j) = rhs(j) +  sum( cross_product( nvec(:,f), pvalues )*awxsi )
                end do
              end do integrate

              deallocate (xsi, cxsi)

!  solve for x (x is returned to rhs vector)
              call lusolver(mat, rhs)
!  fix the face dofs only, the edge dofs are fixed already
              do i = nffe + 1, nffun
                b( lnarray(i) ) = rhs(i)
                dbc( lnarray(i) ) = .true.
              end do
              deallocate(lambda2d, weight2d, lambda, mat, rhs, lnarray )
            end if ! end if BC
          end do faceloop
        end if ! end if vp > 1
      end if bc

!-------------------------------------------------------------------
!    modify the element matrix to enforce essential (Dirichlet) BC
!-------------------------------------------------------------------

      if ( isboundary ) then
        do i = 1, nff
ifdbc:    if ( dbc(i) ) then
!  this is the entry of submatrix A_DD (ref: Matrix Assembly.doc)
!  for Dirichlet dof, let all the columns zero first, then 
!  set the diagonal to unity, so A_DD = I.
!  the subvector b_D is already fixed above by assigning b with the 
!  correct local indexing: b(lni)
            a(i,1:nff) = 0._DP
            a(i,i) = 1._DP
          else
!  this entry has general Neumann, so nullify submatrix A_DN and A_ND;
!  and modify b_N 
            do j = 1, nff
              if ( dbc(j) ) then
                b(i) = b(i) - a(i,j)*b(j)
                a(i,j) = 0._DP
              end if
            end do
          end if ifdbc
        end do
        deallocate (dbc)
      end if
      return

      end subroutine elementmatrix3D
