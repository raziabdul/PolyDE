      subroutine residual(ext,int,matvar,errcode,res,sumres,sumref)
      use feminterface, only: get1Dintegpoints, get2Dintegpoints, getbcval2D, &
                            & getsetting, field, flaech, zeit , pdecoeff, lam2xy
      use femtypes
      use globalvariables, only: ep, nnat, fem_accuracy, n, e, xn, yn, en, &
                               &  kzi, zrb, physics
      implicit none
      integer (I4B) errcode
      real(DP) :: res(:,:), sumres(:), sumref(:)
      logical :: ext, int, matvar
      intent (in) :: ext, int, matvar
      intent (out) :: errcode, res, sumres, sumref
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
!    $Revision: 1.72 $
!    $Date: 2015/04/01 10:55:20 $
!    $Author: juryanatzki $
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
      integer (I4B) :: elem, i, i2, k, k2, inat, start, ende
      integer (I4B) :: rb, bindx, epmax, nb, isucc(3), ipred(3)
      integer (I4B) :: intorder, npkt, npkt1D
      integer (I4B) :: errcodeint, errcode1d = 0, errcode2d = 0
      real (DP) :: endp(2), startp(2), length(3), vec(2), nvec(2)
      real (DP) :: area_e, area_nb, lam1(3), lam2(3), lambdamid(3)
      real (DP) :: reftemp, restemp, xs, ys
      real (DP), allocatable :: resext(:,:,:), resint(:,:), normres(:,:)
      real (DP), allocatable :: evnumin(:,:), reference(:,:)
      real (DP), allocatable :: xtab(:), weight1d(:), weight2d(:), lambda2d(:,:)
      complex (DPC) :: z(15,nnat), pval, qval(nnat), vectmp(2)
      complex (DPC) :: nu(2,2,nnat,nnat), f1(nnat), f2(2,nnat)
      complex (DPC) :: alpha(nnat,nnat), beta(2,nnat,nnat), gamma(2,nnat,nnat)
!  Variables to be obtained from the field subroutine
      complex (DPC) :: u(nnat), alphau(nnat), betagradu(nnat), dgradu(nnat)
      complex (DPC) :: f(nnat),  gammagradu(nnat)
      complex (DPC) :: gradu(2,nnat), nugradu(2,nnat), g(2,nnat), gammau(2,nnat)
!  Variables to be obtained from the field subroutine for the neighbouring element (ext. residual)
      complex (DPC) :: nugradun(2,nnat), gammaun(2,nnat), gn(2,nnat)
      logical :: typ(5), neumann(nnat)
      character (len=16) :: adapttype
!
!  if properly selected this should not happen
      if ((.not. ext) .and. (.not. int)) then
        print "(a)"," **** No type of error estimator chosen!"
        stop
      end if
!  allocation/ initialization
      allocate(reference(n,nnat))
      allocate(evnumin(n,nnat))
      reference = 0._DP
      evnumin = 0._DP
!
      if (ext) then
        allocate(resext(3,n,nnat))
        resext = 0._DP
      end if
      if (int) then
        allocate(resint(n,nnat))
        resint = 0._DP
      end if
      errcode = 0
!
      isucc=(/2,3,1/)  !  successor
      ipred=(/3,1,2/)  !  predecessor
      lambdamid=(/1._DP/3._DP , 1._DP/3._DP , 1._DP/3._DP/)
!
!  Get the physics mode
      call getsetting('PHYSICS_MODE',physics)
!

      select case(physics)
        case('FLUIDINCOMPR','FLUIDELMASSCON')
          do elem = 1,n
            call lam2xy(lambdamid,elem,xs,ys,xn,yn,e)
            call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2)
            do inat = 1,nnat
               evnumin(elem,inat)=abs( nu(1,1,inat,inat) + nu(2,2,inat,inat)-  &
                  &        sqrt( (nu(1,1,inat,inat)-nu(2,2,inat,inat))**2 +    &
                  &        4._DP*nu(1,2,inat,inat)**2 ) )/2._DP
            end do
!  This value doesn't make sense but is used only for programm stability
!  evnumin for pressure is normally 0
            evnumin(elem,3) = 1._DP
          end do
        case default
          do elem = 1,n
            call lam2xy(lambdamid,elem,xs,ys,xn,yn,e)
            call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2)
            do inat = 1,nnat
               evnumin(elem,inat)=abs( nu(1,1,inat,inat) + nu(2,2,inat,inat)-  &
                  &        sqrt( (nu(1,1,inat,inat)-nu(2,2,inat,inat))**2 +    &
                  &        4._DP*nu(1,2,inat,inat)**2 ) )/2._DP
            end do
          end do
      end select
!

!  LOOP OVER ALL ELEMENTS:
      do elem = 1,n
!  calculate the area of the element
        area_e = flaech(xn(e(1,elem)),yn(e(1,elem)),xn(e(2,elem)),yn(e(2,elem)),&
     &                  xn(e(3,elem)),yn(e(3,elem)))
!  edge lengths
        length(1:3) = sqrt((xn(e(isucc(1:3),elem))-xn(e(ipred(1:3),elem)))**2 + &
     &                     (yn(e(isucc(1:3),elem))-yn(e(ipred(1:3),elem)))**2 )
!
!  get nu and compute minimum eigenvalue
!  minumum eigenvalue of diffusion tensor at midpoint of element
!
!  determine the order of numerical integration
        intorder = 2*maxval(ep(elem,:))
!  fetch numerical integration points (Gauss points)
        call get2Dintegpoints(intorder,npkt,weight2d,lambda2d,errcodeint)
!  check if an error occured while doing 2d integration
        if (errcodeint .ge. 1001) errcode2d = errcodeint
!
        typ = (/.true.,.true.,.true.,.true.,.false./)
!  LOOP OVER ALL INTPOINTS:
        do k=1,npkt
!  the field routine should return the values for nu*dgrad(u), beta*grad(u), gamma*grad(u), alpha*u and f
!  Actually the subroutine makes use of constant material values (matvar=false)
          call field(elem,lambda2d(:,k),typ,z,u=u,alphau=alphau,gradu=gradu,&
     &      betagradu=betagradu,gammagradu=gammagradu,nugradu=nugradu,dgradu=dgradu,f=f)
!  LOOP OVER ALL NATURES:
          do inat = 1,nnat
            if (int) then
!------------------------------------------------------------------------------!
!  calculate the interior residual for the element after the formula:          !
!    => resint = f1 + dgrad(u) + gamma*grad(u) - beta*grad(u) - alpha*u        !
!  (f2 is supposed to be constant --> div(f2) vanishes)                        !
!------------------------------------------------------------------------------!
              restemp = abs( f(inat) + dgradu(inat) +                   &
     &          gammagradu(inat) - betagradu(inat) - alphau(inat) )
!  calculate interior residual for elem by summing up over integration points
!  (L2 norm)
              resint(elem,inat) = resint(elem,inat) + weight2d(k) * restemp**2
            end if
!  calculate squared natural energy norm by |(nu*grad(u),grad(u))|+|(alpha*u,u)|
            reftemp = abs( alphau(inat)*u(inat) )                &
     &              + abs( sum(nugradu(1:2,inat)* gradu(1:2,inat)))
            reference(elem,inat) = reference(elem,inat) + weight2d(k) * reftemp
!  END OF NATURE LOOP:
          end do
        end do
        deallocate(weight2d,lambda2d)
!  multiply by triangle area to finish numerical integration
        reference(elem,:) = reference(elem,:) * area_e
        if (int) then
          resint(elem,:) = resint(elem,:) * area_e
!  calculate the error estimate in the energy norm by multiplying with a
!  factor hT^2/nu_min (Verfuerth 2004) and the polynomial degree (Melenk/Wohlmuth)
          resint(elem,:) = maxval(length)**2 /ep(elem,:)**2/evnumin(elem,:) * resint(elem,:)
        end if
!
!
        if (ext) then
!------------------------------------------------------------------------------!
!  calculate the boundary residual for the element after the formular:         !
!    ==> resext =   (nu*grad(u1)|n + gamma*u1 + f2_1)                          !
!                 - (nu*grad(u2)|n + gamma*u2 + f2_2)                          !
!------------------------------------------------------------------------------!
!  LOOP OVER ALL THREE EDGES:
          do i = 1,3
            nb = en(i,elem)
            if (nb .gt. elem) cycle
!  OPERATIONS FOR GEOMETRY:
!  (determination normal vector of edge and of area of neighbor)
!  set start- and endpoints for direction vector (along element edge)
            start  = e(isucc(i),elem)
            ende   = e(ipred(i),elem)
!  set start- and endpoints for direction vector (along element edge)
            startp = (/xn(start),yn(start)/)
            endp   = (/xn(ende) ,yn(ende) /)
!  calculate direction vector
            vec = (endp - startp) / length(i)
!  assign values to outward normal unit vector (perpendicular to element edge)
            nvec = (/vec(2),-vec(1)/)
!  calculate the area of the neighboring element
            area_nb = 0._DP
            if(nb .gt. 0) then
              area_nb = flaech(xn(e(1,nb)),yn(e(1,nb)),xn(e(2,nb)),yn(e(2,nb)),&
     &                         xn(e(3,nb)),yn(e(3,nb)))
            end if
!
!  OPERATIONS FOR NUMERICAL INTEGRATION:
!  get integration points for integrating square of polynomial p --> 2*p
            if (nb .eq. 0) then
              intorder = 2*maxval(ep(elem,:))
            else
              epmax = max(maxval(ep(elem,:)),maxval(ep(nb,:)))
              intorder = 2*epmax
            end if
!  fetch numerical integration points (Gauss-Legendre points)
            call get1Dintegpoints(intorder, xtab, weight1d, npkt1D, errcodeint)
!  check if an error occured while doing 1D integration
            if (errcodeint .ge. 1001) errcode1d = errcodeint
!
!  find the corresponding edge in the neighboring element
            i2 = 0  ! Initialize
            if (nb .gt. 0) then
              do k2 = 1,3
                if (en(k2,nb) .eq. elem) then
                  i2 = k2
                  exit
                end if
              end do
            end if
!
!  DETERMINE IF THE EDGE IS A BRANCH AND LOOK FOR BC:
            if (nb .eq. 0) then
!  get the index of the branch on which the element edge is located
              if (kzi(start) .gt. 0) then
                bindx = kzi(start)
              else if (kzi(start) .lt. 0) then
!  the point is a keypoint
                if (kzi(ende) .gt. 0) then
                  bindx = kzi(ende)
                else if (kzi(ende) .lt. 0) then
!  the point is a keypoint, should not happen both nodes are key-points
                  print*,'*** both nodes are key-points (residual.f90)'
                else
!  one node is an inner node, cannot be a branch                
                  print*,'*** one node is inner node --> cannot be a branch (residual.f90)'
                end if
              else
!  one node is an inner node, cannot be a branch                
                print*,'*** one node is inner node --> cannot be a branch (residual.f90)'
              end if
!
              do inat = 1,nnat
!  set the boundary condition (BC)
                rb = zrb(bindx,inat)
!  check if BC is a Neumann
                if ((rb .ge. 200) .and. (rb .le. 299)) then
                  neumann(inat) = .true.
                else
                  neumann(inat) = .false.
                end if
              end do
            else
              neumann = .false.
            end if
!
!  DO LOOP FOR INTEGRATION OVER ALL INTPOINTS ALONG ELEMENT EDGE:
            do k = 1, npkt1D
!  set typ for fetching grad(u), nu*grad(u), gamma(u) and f2 (or g) for both elements
!
              typ = (/.true.,.true.,.false.,.false.,.true./)
!
!  output from get1dintpoints reaches from -1 to 1. we want to have it from 0 to 1
!  compute the barycentric coordinates for elem
!  lambda(i) has to be 0 to move along edge i
              lam1(i) = 0._DP
              lam1(isucc(i)) = (xtab(k) + 1._DP)/2._DP
              lam1(ipred(i)) = 1._DP - lam1(isucc(i))
              call field(elem,lam1,typ,z,nugradu=nugradu,gammau=gammau,g=g)
!  other wise the branch has a Dirichlet BC
              do inat = 1,nnat
!
!  DO NUMERICAL INTEGRATION ALONG EDGE:
!  calculate solution and residual in L2 norm for intpoint k
!
                if (neumann(inat)) then
!  the branch has a Neumann BC
!  ( nu*grad(u) + gamma*u + g ) * n - (p+sum(q*u))
                  call getbcval2D( bindx, inat, xs, ys, pval, qval)
                  vectmp = nugradu(:,inat) + gammau(:,inat) + g(:,inat)
                  restemp = abs( sum(vectmp * nvec(:))                  &
     &                     - pval - sum(qval(:)*u(:)) )
!  calculate interior residual for elem by summing up over integration points (L2 norm)
                  resext(i,elem,inat) = resext(i,elem,inat) + weight1d(k) * restemp**2
                else if (nb .gt. 0) then
!  the vertex is an inner vertex (no BC)
!  compute the barycentric coordinates for the neighbor
!  lambda(i2) has to be 0 to move along edge i2
                  lam2(i2) = 0._DP
                  lam2(ipred(i2)) = lam1(isucc(i))
                  lam2(isucc(i2)) = lam1(ipred(i))
                  typ = (/.true.,.true.,.false.,.false.,.true./)
                  call field(nb,lam2,typ,z,nugradu=nugradun,gammau=gammaun,g=gn)
!  [ nu*grad(u) + gamma*u + g ] * n
                  vectmp = nugradu(:,inat) -nugradun(:,inat)            &
     &                   + gammau(:,inat) - gammaun(:,inat)             &
     &                   + g(:,inat) - gn(:,inat)
                  restemp = abs( sum(vectmp * nvec(:)) )
!  calculate interior residual for elem by summing up over integration points (L2 norm)
                  resext(i,elem,inat) = resext(i,elem,inat) + weight1d(k) * restemp**2
                end if
!  END OF NATURE LOOP:
              end do
!  END OF NUMERICAL INTEGRATION LOOP:
            end do
            deallocate(xtab,weight1d)
!  multiply by edge length/2 to finish numerical integration
            resext(i,elem,:)=resext(i,elem,:)*length(i)*0.5_DP
!
!  calculate the error estimate in the energy norm by multiplying with a
!  factor h/nu_min (Verfuerth 2004) and the polynomial degree (Melenk/Wohlmuth)
            do inat=1,nnat
              if (nb .gt. 0) then
!  distribute the residual between elem and neighbor by comparing area,
                resext(i2,nb,inat)  = area_nb/(area_e + area_nb)        &
     &                  * length(i) / min(ep(elem,inat),ep(nb,inat))    &
     &                  / max(evnumin(elem,inat),evnumin(nb,inat))      &
     &                  * resext(i,elem,inat)
                resext(i,elem,inat) = area_e /(area_e + area_nb)        & 
     &                  * length(i) / min(ep(elem,inat),ep(nb,inat))    &
     &                  / max(evnumin(elem,inat),evnumin(nb,inat))      &
     &                  * resext(i,elem,inat) 
              else
!  we have a boundary --> one element get whole residual
                resext(i,elem,inat) = length(i) / ep(elem,inat)         &
     &                  / evnumin(elem,inat) * resext(i,elem,inat) 
              end if
            end do
!
!  END OF EDGE LOOP:
          end do
        end if
!
!  END OF ELEMENT LOOP:
      end do  !elem=1,n
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
!  print values for global error indicator and relative error in energy norm
      select case(physics)
        case('FLUIDINCOMPR','FLUIDELMASSCON')
          do inat = 1,nnat
            if (inat .eq. 3) then
              print "(a,i3,a/a,i3,a)","Nature : ",3,"   global error indicator : NOT GIVEN FOR PRESSURE!" &
                                     ,"Nature : ",3,"      relative error in % : NOT GIVEN FOR PRESSURE!"
              cycle
            end if
            sumres(inat) = sum(res(:,inat))
            sumref(inat) = sum(reference(:,inat))
            fem_accuracy = sqrt(sumres(inat)/sumref(inat))
            print "(a,i3,a,g13.6/a,i3,a,g13.6)","Nature : ",inat,"   global error indicator : ",sqrt(sumres(inat)) &
                                               ,"Nature : ",inat,"      relative error in % : ",100._DP*fem_accuracy
          end do
        case default
          do inat = 1,nnat
            sumres(inat) = sum(res(:,inat))
            sumref(inat) = sum(reference(:,inat))
            fem_accuracy = sqrt(sumres(inat)/sumref(inat))
            print "(a,i3,a,g13.6/a,i3,a,g13.6)","Nature : ",inat,"   global error indicator : ",sqrt(sumres(inat)) &
                                               ,"Nature : ",inat,"      relative error in % : ",100._DP*fem_accuracy
          end do
      end select
!
!  For h-adaptation normalise the residual of individual elements in all
!  natures with respect to the energies in other natures and the coupling energies
!  Use temporary array normres to do the calculations and rewrite result
!  in res(elem,inat)
      call getsetting('ADAPTION_TYPE',adapttype)
      adapttype = adapttype(1:len_trim(adapttype))
      if (adapttype .eq. 'H_ADAPT' .or. adapttype .eq. 'HP_ADAPT') then
        allocate(normres(n,nnat))
        normres = 0._DP
        select case(physics)
          case('FLUIDINCOMPR','FLUIDELMASSCON')
            do inat = 1,nnat
              if (inat .eq. 3) cycle
              normres(1:n,inat) = res(1:n,inat)
!              do jnat = 1,nnat
!                if ((inat .eq. jnat) .or. (jnat .eq. 3)) cycle
!                where (abs(reference(1:n,inat,jnat)) .gt. tiny(1._DP))
!                  normres(1:n,inat) = normres(1:n,inat) &
!                    & + (reference(1:n,inat,jnat)/reference(1:n,jnat,jnat)) &
!                    & * res(1:n,jnat)
!                end where
!              end do
              if (abs(sum(reference(:,inat))) .gt. tiny(1._SP)) then
                normres(1:n,inat) = normres(1:n,inat)/sum(reference(:,inat))
              end if
            end do
          case default
            do inat = 1,nnat
              normres(1:n,inat) = res(1:n,inat)
!              do jnat = 1,nnat
!                if (inat .eq. jnat) cycle
!                where (abs(reference(1:n,inat,jnat)) .gt. tiny(1._DP))
!                  normres(1:n,inat) = normres(1:n,inat) &
!                    & + (reference(1:n,inat,jnat)/reference(1:n,jnat,jnat)) &
!                    & * res(1:n,jnat)
!                end where
!              end do
              if (abs(sum(reference(:,inat))) .gt. tiny(1._SP)) then
                normres(1:n,inat) = normres(1:n,inat)/sum(reference(:,inat))
              end if
            end do
        end select
        res = normres
        deallocate (normres)
      end if
!
      deallocate (reference)
      deallocate (evnumin)
!
      if (allocated(resext)) deallocate(resext)
      if (allocated(resint)) deallocate(resint)
!

      call zeit('calculating the error estimator')
      return
      end subroutine
