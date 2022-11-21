      subroutine residual3D(res,sumref)
      use feminterface, only: get2Dintegpoints, invertmat3, zeit, cross_product
      use feminterface3D, only: area3D, bcstovv, field3D, getgl,&
                                getbcval_vec, lam2xyz, pdecoeff3D
      use femtypes
      use globalvariables3D, only: sbc, eps0, mu0, nod, numv, omega, vn, vp, vv, bctype
      implicit none
      real(DP), pointer :: res(:)
      real(DP), intent(out) :: sumref
!
!------------------------------------------------------------------------------
!    $Revision: 1.14 $
!    $Date: 2015/11/10 15:39:24 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!  Calculate the residual for each element based on the jump gradient over the
!  element's faces.
!
!  TO DO:
!  The error estimate is not proven to be correct. Must be done based on
!  analytically known problem.
!  
!  Output:
!       res     local error indicators for each volume element as a vector
!       sumref  reference value in squared energy norm for solution
!
!  Internal variables:
      integer (I4B) :: i, k, elem, face     ! loop counters
      integer (I4B) :: nb, nbface, bcnum
      integer (I4B) :: intorder, npkt, errcode      ! get2dintegpoints
      integer (I4B), parameter :: fnodes(3,4)=reshape((/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/))
      real (DP) :: area_f, reftemp, reftempnb, restempn, restempt
      real (DP) :: lamfield(4), vert(3,3), xyzt(3)
      real (DP) :: gl(3,4), n(3,4), t1(3,4), t2(3,4)          ! vectors on face
      real (DP), allocatable :: weight2d(:), lambda2d(:,:)    ! get2dintegpoints
      real (DP) :: k0   ! free space wave number
      real (DP), allocatable :: resn(:), rest(:), reference(:,:)
      complex (DPC) :: u1(3), curlu1(3), u2(3), curlu2(3), jumpn, jumpt
      complex (DPC) :: nu1(3,3), alpha1(3,3), beta1(3), f11(3), f21(3)    ! pdecoeffs
      complex (DPC) :: nu2(3,3), alpha2(3,3), beta2(3), f12(3), f22(3)    ! pdecoeffs
      complex (DPC) :: nu1inverted(3,3), nu2inverted(3,3)  ! inverse matrices
      complex (DPC) :: p2(3), q(3,3)    ! Neumann bc
      logical :: neumann
!------------------------------------------------------------------------------
!
!  Initializations:
      k0 = omega * sqrt(mu0*eps0)
      allocate(resn(numv),rest(numv),reference(4,numv))
      resn = 0._DP
      rest = 0._DP
      reference = 0._DP
!
!  Assign list of bc, if not previously done.
      if (minval(vv) .eq. 0) call bcstovv
!
elems:do elem = 1,numv    ! loop over all volume elements
  faces:do face = 1,4    ! loop over 4 faces of element
!  Determine neighbour of elem at local face.
          nb = vv(face,elem)
!  Cycle to next face, if neighbour was treated before.
          if (nb .lt. elem) cycle faces
!
!  Check for neighbourhood information and assign integration order, bc and
!  local number of neighbour's face.
          if (nb .lt. 0) then
! For now only have single nature, so vp(elem,1)
            intorder = 2*vp(elem,1)   ! integration order
            bcnum = bctype(sbc(abs(vv(face,elem))),1)
!  Check for Neumann boundary condition on face.
            if ((bcnum .ge. 200) .and. (bcnum .le. 299) ) then
              neumann = .true.
            else
              neumann = .false.
            end if
          else if (nb .gt. 0) then
            intorder = 2*min(vp(elem,1),vp(nb,1))     ! integration order
!  Determine local face of neighbour.
            do i=1,4    ! loop over 4 faces of neighbour
              if(vv(i,nb) .eq. elem) then
                nbface = i
                exit
              end if
            end do
          end if
!
!  Fetch numerical integration points.
          call get2Dintegpoints(intorder,npkt,weight2d,lambda2d,errcode)
!
!  Get normal vector for face.
          call getgl(nod( 1:3,vn(1:4,elem) ), gl, n, t1, t2)
!
!  Clear restemp for element face.
          restempn = 0._DP
          restempt = 0._DP
!  Clear reftemp for element face.
          reftemp   = 0._DP
          reftempnb = 0._DP
!
    integ:do k = 1, npkt    ! numerical integration loop
!  Compute the barycentric coordinates for intpoint on face.
            lamfield(face) = 0._DP
            lamfield(fnodes(:,face)) = lambda2d(:,k)
!  Determine world coordinates from lamfield and get pde coefficients for elem.
            call lam2xyz(lamfield,elem,xyzt,nod,vn)
            call pdecoeff3D(elem,xyzt,nu1,alpha1,beta1,f11,f21)
            call invertmat3(nu1,nu1inverted)
!  Get potential u and curl(u) from field3D for elem
            call field3D(elem,lamfield,u1,curlu1)
!
!  Calculate reference value in energy norm for elem
            reftemp = abs( sum(matmul(nu1,curlu1)*curlu1) + sum(matmul(nu1,curlu1)*u1) )
!
!  Face has a neighbour.
            if (nb .gt. 0) then
!  Compute the barycentric coordinates for intpoint on nbface.
              lamfield(nbface) = 0._DP
              lamfield(fnodes(:,nbface)) = lambda2d(:,k)
!  Determine world coordinates from lamfield and get pde coefficients for nb.
              call lam2xyz(lamfield,nb,xyzt,nod,vn)
              call pdecoeff3D(nb,xyzt,nu2,alpha2,beta2,f12,f22)
              call invertmat3(nu2,nu2inverted)
!  Get potential u and curl(u) from field3D for nb.
              call field3D(nb,lamfield,u2,curlu2)
!
!  Calculate reference value in energy norm for nb
              reftempnb = abs( sum(matmul(nu2,curlu2)*curlu2) + sum(matmul(nu2,curlu2)*u2) )
!  Calculate jumps
              jumpn = sum( (matmul(alpha2,u2) - matmul(alpha1,u1))*n(:,face) )
              jumpt = sqrt(sum(cross_product(cmplx(n(:,face),0._DP,DPC),((matmul(nu2,curlu2)-f22) - &
                                            (matmul(nu1,curlu1)-f21)))**2))
!  Face has a Neumann BC.
            else if ((nb .lt. 0) .and. neumann) then
              call getbcval_vec(sbc(abs(vv(face,elem))), 1, xyzt, p2, q)    ! only a single nature
              jumpn = 0._DP
              jumpt = sqrt(sum(( (p2 + matmul(q,u1)) - &
                                  cross_product(cmplx(n(:,face),0._DP,DPC),(matmul(nu1,curlu1)-f21)) )**2))
!  Face has a Dirichlet BC.
!  TO BE DONE
            else
              restempn = 0._DP
              restempt = 0._DP
            end if
!  Compute reference value and residual for neighbouring elements in energy norm.
            reftemp   = reftemp   + weight2d(k) * reftemp
            reftempnb = reftempnb + weight2d(k) * reftempnb
            restempn  = restempn  + weight2d(k) * abs(jumpn)**2 * &
                                    0.5_DP*maxval(abs(nu1inverted + nu2inverted))
            restempt  = restempt  + weight2d(k) * abs(jumpt)**2 * &
                                    0.5_DP*maxval(abs(nu1inverted + nu2inverted))
          end do integ   ! end of numerical integration loop
          deallocate(weight2d,lambda2d)
!
!  Determine vertex coordinates of face and compute area.
          vert(:,1) = nod(:,vn(fnodes(1,face),elem))
          vert(:,2) = nod(:,vn(fnodes(2,face),elem))
          vert(:,3) = nod(:,vn(fnodes(3,face),elem))
          area_f = area3d(vert)
!
          if (nb .gt. 0) then
!  Divide residual by two, if neighbour present.
!  has to be multiplied with scaling factor area_f / p^2 and area_f (integration)
            resn(elem) = resn(elem) + 0.5_DP * area_f**2 / vp(elem,1)**2 * restempn 
            resn(nb) = resn(nb) + 0.5_DP * area_f**2 / vp(nb,1)**2 * restempn 
            rest(elem) = rest(elem) + 0.5_DP * area_f**2 / vp(elem,1)**2 * restempt 
            rest(nb) = rest(nb) + 0.5_DP * area_f**2 / vp(nb,1)**2 * restempt
!  assign reference value to elements
            reference(face,elem) = area_f * reftemp
            reference(nbface,nb) = area_f * reftempnb
          else
!  We have a boundary condition --> one element get whole residual.
!  has to be multiplied with scaling factor area_f / p^2 and area_f (integration)
            resn(elem) =  area_f**2 / vp(elem,1)**2 * restempn
            rest(elem) =  area_f**2 / vp(elem,1)**2 * restempt
          end if
        end do faces   ! end of face loop
      end do elems   ! end of element loop
!
!  compute squared residual in energy norm
      if (associated(res)) deallocate (res)
      allocate(res(numv))
      res = resn + rest
      deallocate(resn, rest)
!  compute squared reference value for whole solution
      sumref = sum(reference)
      deallocate(reference)
!
!  plot the residual of the solution as text
      print "(a,g10.3/,a,g10.3)","global error indicator : ",sqrt(sum(res)) &
     &                          ,"   relative error in % : ",100._DP*sqrt(sum(res)/sumref)

!
      call zeit('calculating the error')
!
      end subroutine
