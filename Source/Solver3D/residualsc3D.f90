      subroutine residualsc3D(ext,int,matvar,errcode,res,sumref,resextout)
      use feminterface3D, only: get3Dintegpoints, getbcval, fieldsc3D, pdecoeff3D, lam2xyz, getgl, area3D, tetvolume
      use feminterface,   only: zeit, get2Dintegpoints
      use femtypes
      use globalvariables3D, only : sbc, nnat, nod, numv, physics, vn, vp, vv, bctype
      use matconstants, only:
      implicit none
      integer (I4B) errcode
      real(DP):: res(:,:), sumref(:)
      real(DP), pointer, optional :: resextout(:,:,:)
      logical :: ext, int, matvar
      intent (in) :: ext, int, matvar 
      intent (out) :: errcode, res, sumref, resextout
!
!-------------------------------------------------------------------------------
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
!    $Date: 2015/11/10 15:40:57 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
!
!  Calculate the residual for each element. The residual is given in squared L2
!  norm. The output is a vector of squared local error indicators and a sum for
!  the quadratic natural energy norm over all elements as a reference value.
!  Note that sumref cannot be computed properly if input variable int is set to
!  .false.
!
!  The calculation of the residual is based on two papers:
!
!  R. Verfuerth, "Robust a posteriori error estimates for stationary convection-
!  diffusion equations," unpublished (2004)
!
!  J.M. Melenk, and B.I. Wohlmuth, "On residual-based-a-posteriori error estimation
!  in hp-FEM," Adv. Comput. Math. 15, 311-331 (2001).
!
!  The residual is currently calculated after Melenk/Wohlmuth. The energy norm
!  computation is done following Verfuerth
!
!  The implementation is based on the 2D residuals evaluation
!
!------------------------------------------------------------------------------
!
!  Input:
!            ext      =.true. if boundary error indicator should be calculated
!            int      =.true. if interior error indicator should be calculated
!                     (NB: IF .FALSE. NO CALCULATION OF RELATIVE ERROR POSSIBLE)
!            matvar   =.true. if the material coefficients are varying across the element
!  Output:
!            res      quadratic local error indicators for each element
!                     for each nature as a matrix
!            sumref   reference value as scale basis for relative error
!                     (sum of |nu*grad(u)|^2 over all elements)
!  Errorcodes:
!            errcode  3001 = not used atm
!
!  Internal variables:
! 
!  intorder      order of a polynomial to be integrated
!
!
!  internal variables:
      integer (I4B) :: elem, i, j, i2, k, k2, inat, jnat, l, epmax
      integer (I4B) :: intorder, npkt, npkt2D, nb
      integer (I4B) :: errcodeint, errcode1d = 0, errcode2d = 0
      integer (I4B), parameter :: f2n(3,4)=reshape((/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/))
      integer (I4B), parameter :: e2n(2,6)=reshape((/1,2,2,3,1,3,1,4,2,4,3,4/),(/2,6/))
      real (DP) :: reftemp, restemp, sumres(nnat), lam2(4), vol_e, vol_nb
      real (DP) :: length(6), nvec(3), norvec(3,4)
      real (DP) :: vert(3,4), vert_nb(3,4), point(3), gl(3,4), xyzt(3), vert2(3,3), area_f(4)
      real (DP) :: fem_accuracy
      real (DP), allocatable :: resint(:,:)
      real (DP), allocatable :: evnumin(:,:), reference(:,:,:)
      real (DP), allocatable :: weight3d(:), lambda3d(:,:)
      real (DP), allocatable :: weight2d(:), lambda2d(:,:), lambdaface(:,:)
      real (DP), pointer :: resext(:,:,:)=>null()
      complex (DPC) :: z(15,nnat), pval, qval(nnat), vectmp(3)
      complex (DPC) :: nu(3,3,nnat,nnat), f1(nnat), f2(3,nnat)
      complex (DPC) :: alpha(nnat,nnat), beta(3,nnat,nnat), gamma(3,nnat,nnat)
!  Variables to be obtained from the field subroutine
      complex (DPC) :: u(nnat), alphau(nnat,nnat), gammau(3,nnat,nnat)
      complex (DPC) :: dgradu(nnat,nnat), f(nnat), g(3,nnat), gradu(3,nnat)
      complex (DPC) :: betagradu(nnat,nnat), gammagradu(nnat,nnat), nugradu(3,nnat,nnat)
!  Variables to be obtained from the field subroutine for the neighbouring element (ext. residual)
      complex (DPC) :: un(nnat), gradun(3,nnat), nugradun(3,nnat,nnat)
      complex (DPC) :: gammaun(3,nnat,nnat), gn(3,nnat)
      logical :: typ(5), neumann(nnat)
      character (len=30) :: inature

!
!  if properly selected this should not happen
      if ((.not. ext) .and. (.not. int)) then
        print "(a)"," **** No type of error estimator chosen!"
        stop
      end if
!  allocation/ initialization
!      allocate(normres(numv,nnat))
      allocate(reference(numv,nnat,nnat), evnumin(numv,nnat))
      reference = 0._DP
      if (ext) then
        if (present(resextout)) then
          if (associated(resextout)) then
            resext => resextout
          else
            allocate(resext(4,numv,nnat))
            resextout => resext
          end if
        else
         allocate(resext(4,numv,nnat))
        end if
        resext = 0._DP
      end if
      if (int) then
        allocate(resint(numv,nnat))
        resint = 0._DP
      end if
      errcode = 0

!$OMP PARALLEL DO default (none)                                        &
!$OMP shared( nod, vn, numv, nnat, evnumin )                            &       ! global variables
!$OMP private( elem, vert, point, nu, alpha, beta, f1, f2, gamma )
      do elem = 1,numv
! get vertices coordinate
        vert(1:3,1:4) = nod( 1:3,vn(1:4,elem) )
!
! get mid-point
        point(1) = sum( vert(1,:) )
        point(2) = sum( vert(2,:) )
        point(3) = sum( vert(3,:) )
        point=point/4._DP
!  get nu and compute minimum eigenvalue
        call pdecoeff3D(elem,point,nu,alpha,beta,f1,f2,gamma)
!  minumum eigenvalue of diffusion tensor at midpoint of element
        do inat = 1,nnat
! TO DO here we should calculate the minimum egenvalue of the tensor nu
          evnumin(elem,inat)=min(abs( nu(1,1,inat,inat)) , abs(nu(2,2,inat,inat) ), abs(nu(3,3,inat,inat)))
        end do
      end do
!$OMP END PARALLEL DO
!
! physics, matvar, errcode, res, sumref
!  LOOP OVER ALL ELEMENTS:
!$OMP PARALLEL DO default (none)                                        &
!$OMP shared( sbc, bctype, nnat, nod, numv, vn, vp, vv)                 &       ! global variables
!$OMP shared( errcode2d, int, ext, f2n, e2n, errcode1d)                 &
!$OMP shared( res, resint, resext, reference, evnumin)                  &
!$OMP private(elem, vert, vol_e, i, j, i2, vert2, typ, neumann)         &
!$OMP private(z, u, alphau, gradu, betagradu, gammagradu, nugradu)      &
!$OMP private(dgradu, f, g, un, nugradun, gammaun, gn, gradun, nb, gl)  &
!$OMP private(nvec, norvec, vert_nb, epmax, vol_nb, pval, qval, xyzt)   &
!$OMP private(intorder, npkt, npkt2D, weight3d, lambda3d, weight2d)     &
!$OMP private(lambda2d, lambdaface, lam2, errcodeint, length, area_f)   &
!$OMP private(vectmp, restemp, gammau, reftemp, inat, jnat, k, k2, l)
      do elem = 1,numv
! get vertices coordinate
        vert(1:3,1:4) = nod( 1:3,vn(1:4,elem) )
        vol_e = abs(tetvolume(vert))
!  calculate the area of the faces of the element
        do i=1,4
          vert2(:,1)=nod(:,vn(f2n(1,i),elem))
          vert2(:,2)=nod(:,vn(f2n(2,i),elem))
          vert2(:,3)=nod(:,vn(f2n(3,i),elem))
          area_f(i)=abs(area3d(vert2))
        end do
!  calculate length of edges
        do l = 1, 6
          length(l) = sqrt( ( vert(1,e2n(1,l)) - vert(1,e2n(2,l)) )**2  +  &
     &                      ( vert(2,e2n(1,l)) - vert(2,e2n(2,l)) )**2  +  &
     &                      ( vert(3,e2n(1,l)) - vert(3,e2n(2,l)) )**2 )
        end do
!
!  determine the order of numerical integration
        intorder = 2*maxval(vp(elem,1:nnat))
!  fetch numerical integration points (Gauss points)
        call get3Dintegpoints(intorder,npkt,weight3d,lambda3d,errcodeint)
!  check if an error occured while doing integration
        if (errcodeint .ge. 1001) errcode2d = errcodeint
!
!  LOOP OVER ALL INTPOINTS:
        do k=1,npkt
!  the field routine should return the values for nu*dgrad(u), beta*grad(u), gamma*grad(u), alpha*u and f
!  TO DO   actually the routine does not make use of constant material values (matvar=false)
          if (int) then
            typ = (/.true.,.true.,.true.,.true.,.false./)
            call fieldsc3D(elem,lambda3d(:,k),typ,z,u=u,alphau=alphau,gradu=gradu,          &
     &      betagradu=betagradu,gammagradu=gammagradu,nugradu=nugradu,dgradu=dgradu,f=f,g=g)
          else
            typ = (/.true.,.true.,.false.,.false.,.false./)
            call fieldsc3D(elem,lambda3d(:,k),typ,z,u=u,alphau=alphau,gradu=gradu,          &
     &      betagradu=betagradu,gammagradu=gammagradu,nugradu=nugradu,f=f,g=g)
          end if 
!  LOOP OVER ALL NATURES:
          do inat = 1,nnat
            if (int) then

!------------------------------------------------------------------------------!
!  calculate the interior residual for the element after the formula:          !
!    => resint = f1 + nu*dgrad(u) + gamma*grad(u) - beta*grad(u) - alpha*u     !
!  (f2 is supposed to be constant --> div(f2) vanishes)                        !
!------------------------------------------------------------------------------!
              restemp = abs(  f(inat) + sum( (dgradu(inat,:)) + & !*(nu(1,1,inat,:)) +   &  !
     &          gammagradu(inat,:) - betagradu(inat,:) - alphau(inat,:)) )
!  calculate interior residual for elem by summing up over integration points
!  (L2 norm)
              resint(elem,inat) = resint(elem,inat) + weight3d(k) * restemp**2
            end if
!  calculate squared natural energy norm by |(nu*grad(u),grad(u))|+|(alpha*u,u)|
            do jnat = 1,nnat
              reftemp = abs( alphau(inat,jnat)*u(inat) )                &
     &                + abs( sum(nugradu(1:3,inat,jnat)* gradu(1:3,inat)))
              reference(elem,inat,jnat) = reference(elem,inat,jnat) + weight3d(k) * reftemp
            end do
!  END OF NATURE LOOP:
          end do
        end do
!  END OF LOOP OVER ALL INTPOINTS:
        deallocate(weight3d,lambda3d)
!  multiply by tetrahedron volume to finish numerical integration
        reference(elem,:,:) = reference(elem,:,:) * vol_e
        if (int) then
          resint(elem,:) = resint(elem,:) * vol_e
!  calculate the error estimate in the energy norm by multiplying with a
!  factor hT^2/nu_min (Verfuerth 2004) and the polynomial degree (Melenk/Wohlmuth)
!          resint(elem,:) = maxval(area_f)**2 /vp(elem,:)**2/evnumin(elem,:) * resint(elem,:)
          resint(elem,:) = maxval(length)**2 /vp(elem,:)**2/evnumin(elem,:) * resint(elem,:)
        end if
!
!
        if (ext) then
!------------------------------------------------------------------------------!
!  calculate the boundary residual for the element after the formular:         !
!    ==> resext =   (nu*grad(u1)|n + gamma*u1 + f2_1)                          !
!                 - (nu*grad(u2)|n + gamma*u2 + f2_2)                          !
!------------------------------------------------------------------------------!
!  LOOP OVER ALL FOUR FACES:
          call getgl(vert, gl, norvec)
          do i = 1,4
            nb = vv(i,elem)
            if (nb .gt. elem) cycle
!  OPERATIONS FOR GEOMETRY:
            nvec(:)=norvec(:,i)

!
!  OPERATIONS FOR NUMERICAL INTEGRATION:
!  get integration points for integrating square of polynomial p --> 2*p
            if (nb .le. 0) then
              intorder = 2*maxval(vp(elem,:))
            else
              epmax = max(maxval(vp(elem,:)),maxval(vp(nb,:)))
              intorder = 2*epmax
              vert_nb(1:3,1:4) = nod( 1:3,vn(1:4,nb) )
              vol_nb = abs(tetvolume(vert_nb))
            end if
!  fetch numerical integration points (Gauss-Legendre points)
            call get2Dintegpoints(intorder, npkt2D, weight2d, lambda2d, errcodeint)
            allocate(lambdaface(4,npkt2D))
            lambdaface(:,:)=0._DP
            do j=1,3
              lambdaface(f2n(j,i),:)=lambda2d(j,:)
            end do

!  check if an error occured while fetching integration points
            if (errcodeint .ge. 1001) errcode1d = errcodeint
!
!  find the corresponding face in the neighboring element
            if (nb .gt. 0) then
              do k2 = 1,4
                if (vv(k2,nb) .eq. elem) then
                  i2 = k2
                  exit
                end if
              end do
            end if
!
!  DETERMINE IF THE FACE IS A SURFACE ELEMENT AND LOOK FOR BC:
            neumann(:)=.false.
            if (nb .lt. 0) then
              do inat=1,nnat
                if ( ( bctype(sbc(abs(vv(i,elem))),inat) .lt. 201) .or. ( bctype(sbc(abs(vv(i,elem))),inat) .ge. 300 )) then
                  neumann(inat) = .false.
                else
                  neumann(inat) = .true.
                end if
              end do
            end if
!
!  DO LOOP FOR INTEGRATION OVER ALL INTPOINTS ALONG ELEMENT EDGE:
            do k = 1, npkt2D
!  set typ for fetching grad(u), nu*grad(u), gamma(u) and f2 (or g) for both elements
!
              typ = (/.true.,.true.,.false.,.false.,.true./)
              call lam2xyz( lambdaface(:,k), elem, xyzt, nod, vn )
              call fieldsc3D(elem,lambdaface(:,k),typ,z,u=u,nugradu=nugradu,gammau=gammau,g=g,gradu=gradu)
!
              do inat = 1,nnat
!
!  DO NUMERICAL INTEGRATION ON SURFACE:
!  calculate solution and residual in L2 norm for intpoint k
!
                if (neumann(inat)) then
!  the branch has a Neumann BC
!  ( nu*grad(u) + gamma*u + g ) * n - (p+sum(q*u))

                  call getbcval( sbc(abs(vv(i,elem))), inat, xyzt, pval, qval)
                  vectmp = sum( nugradu(:,inat,:) + gammau(:,inat,:)    &
     &                      , dim=2)+ g(:,inat)
                  restemp = abs( sum(vectmp * nvec(:))                  &
     &                     - pval - sum(qval(:)*u(:)) )
        !          if (inat .eq. 1 ) print*,'restemp = ', restemp
!  calculate interior residual for elem by summing up over integration points (L2 norm)
                  resext(i,elem,inat) = resext(i,elem,inat) + weight2d(k) * restemp**2
                else if (nb .gt. 0) then
                  typ = (/.true.,.true.,.false.,.false.,.true./)
                  lam2(:)=0._DP
                  do j=1,3
                    lam2(f2n(j,i2))=lambda2d(j,k)
                  end do
                  call fieldsc3D(nb,lam2,typ,z,u=un,nugradu=nugradun,gammau=gammaun,g=gn,gradu=gradun)
!
                  vectmp = sum( nugradu(:,inat,:) -nugradun(:,inat,:)   &
     &                     + gammau(:,inat,:) - gammaun(:,inat,:)       &
     &                     , dim=2) + g(:,inat) - gn(:,inat)
                  restemp = abs( sum(vectmp * nvec(:)) )
!  calculate interior residual for elem by summing up over integration points (L2 norm)
                  resext(i,elem,inat) = resext(i,elem,inat) + weight2d(k) * restemp**2
                end if
!  END OF NATURE LOOP:
              end do
!  END OF NUMERICAL INTEGRATION LOOP:
            end do
            deallocate(lambda2d,weight2d,lambdaface)
!  multiply by edge length/2 to finish numerical integration
            resext(i,elem,:)=resext(i,elem,:)*area_f(i)*0.5_DP
!
!  calculate the error estimate in the energy norm by multiplying with a
!  factor h/nu_min (Verfuerth 2004) and the polynomial degree (Melenk/Wohlmuth)
            do inat=1,nnat
              if (nb .gt. 0) then
!  distribute the residual between elem and neighbor by comparing volume,
                resext(i2,nb,inat)  = vol_nb/(vol_e + vol_nb)        &
     &                  * area_f(i) / (2._DP*min(vp(elem,inat),vp(nb,inat)))    &
     &                  / max(evnumin(elem,inat),evnumin(nb,inat))      &
     &                  * resext(i,elem,inat)
                resext(i,elem,inat) = vol_e /(vol_e + vol_nb)        & 
     &                  * area_f(i) / (2._DP*min(vp(elem,inat),vp(nb,inat)))    &
     &                  / max(evnumin(elem,inat),evnumin(nb,inat))      &
     &                  * resext(i,elem,inat) 
              else
!  we have a boundary --> one element get whole residual
                resext(i,elem,inat) = area_f(i) / vp(elem,inat)         &
     &                  / evnumin(elem,inat) * resext(i,elem,inat)
              end if
            end do
!
!  END OF FACE LOOP:
          end do
        end if
!
!  END OF ELEMENT LOOP:
      end do
!$OMP END PARALLEL DO
!
!  COMPUTE THE WHOLE RESIDUAL
      if (int .and. (.not. ext)) then
        res = resint(:,:)
      else if (ext .and. (.not. int)) then
!  sum over edges
        res = sum(resext(:,:,:),dim=1)
      else if (int .and. ext) then
        res = resint(:,:) + sum(resext(:,:,:),dim=1)
      end if

      do inat = 1,nnat
        write(inature,*) inat
        inature=adjustl(inature)

        print*,'-------------Nature ',inat,' ----------------------'
!        print*,'sum resint'//trim(inature),' = ', sum(resint(:,inat))
        print*,'sum resext'//trim(inature),' = ', sum(resext(:,:,inat))
        sumres(inat) = sum(res(:,inat))
        sumref(inat) = sum(reference(:,inat,:))
        if(sumref(inat) .ne. 0.) then
          fem_accuracy = sqrt(sumres(inat)/sumref(inat))
        else
          fem_accuracy = 0._DP
        end if
        print*,'sumref'//trim(inature),' = ',sumref(inat)
        print*,"global error indicator"//trim(inature),":",sqrt(sumres(inat)) 
        print*,"relative error"//trim(inature)," in %:",100._DP*fem_accuracy
      end do
!      print "(a,g10.3)
!
!  Normalise the residual of individual elements in all natures with
!  respect to the energies in other natures and the coupling energies
!  Use temporary array normres to do the calculations and rewrite result
!  in res(elem,inat)
!
   !   do inat = 1,nnat
   !     normres(1:n,inat) = res(1:n,inat)
   !     do jnat = 1,nnat
   !       if (inat .eq. jnat) cycle
   !       where (abs(reference(1:n,inat,jnat)) .gt. tiny(1._DP))
   !         normres(1:n,inat) = normres(1:n,inat) &
   !    &              + (reference(1:n,inat,jnat)/reference(1:n,jnat,jnat)) &
   !    &              * res(1:n,jnat)
   !       end where
   !     end do
   !     normres(1:n,inat) = normres(1:n,inat)/sum(reference(:,inat,inat))
   !   end do
   !   res = normres
!
   !   deallocate (normres)
      deallocate (reference)
      deallocate (evnumin)
!
      if ( .not. present(resextout) .and. associated(resext) ) deallocate(resext)
      if (allocated(resint)) deallocate(resint)
!
      call zeit('calculating the error')
      return
      end subroutine residualsc3D
 
