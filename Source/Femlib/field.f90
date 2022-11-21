      subroutine field(elem,lambda,typ,z,soln,u,alphau,gammau,gradu,betagradu,gammagradu,nugradu,dgradu,f,g)
      use feminterface, only: shapefunction, pdecoeff, dgshape, lam2xy, &
                              invertmat
      use femtypes
      use globalvariables, only: e, eg, x, ndof, xn, yn, nnat, ep
      implicit none
      integer (I4B) elem
      real (DP) lambda(3)
      logical typ(5)
      complex (DPC) z(:,:)      !  second index for nature
      complex (DPC), optional :: soln(:)
!  scalars (indices are nature)
      complex (DPC), optional :: u(:), alphau(:), betagradu(:), gammagradu(:),  dgradu(:), f(:)
!  vectors (first index 1,2 represents the spatial components x,y), second index is nature
      complex (DPC), optional :: gammau(:,:), gradu(:,:), nugradu(:,:), g(:,:) 
      intent (in) ::  elem, lambda, typ, soln
      intent (out) ::  z, u, alphau, gammau, gradu, betagradu, gammagradu, nugradu, dgradu, f, g
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
!    $Revision: 1.48 $
!    $Date: 2015/06/02 09:33:58 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
!
!  Evaluate the potential, the gradient, the dgradient and products with the
!  corresponding pde coefficients.
!  (Remark: The material is assumed to be constant across the element.)
!  The subroutine will return all values connected with a set type in logical
!  vector typ.
!    - if typ(1) is .true. everything connected with u is returned
!      (u, alpha*u, gamma*u and second term for determining energy density)
!
!    - if typ(2) is .true. everything connected with grad(u) is returned
!      (grad(u), beta*grad(u), gamma*grad(u), nu*grad(u) and the curl))
!
!    - if typ(3) is .true. the 2nd gradient dgrad(u) will be returned
!      (it already contains nu like programmed in dgshape.f90 subroutine)
!
!    - if typ(4) is .true. the scalar source term f1 (also called f) is returned
!
!    - if typ(5) is .true. the sourcevector f2 (also called g) is returned
!
!  Input:
!     elem      element number
!     lambda    vector of barycentric coordinates
!     typ       case selector integer
!
!  Output:
!     z         value for typ (all dof values summed up)
!               (this variable is a vector and should be allocated by the
!               calling subroutine)
!     z( 1)    = u                  u                   scalar
!     z( 2)    = alphau             alpha*u             scalar
!     z( 3)    = ualphau            u*conj(alpha*u)     scalar
!     z( 4, 5) = gradu              grad(u)             vector
!     z( 6)    = betagradu          beta*grad(u)        scalar
!     z( 7, 8) = nugradu            nu*grad(u)          vector
!     z( 9)    = dgradu             div(nu*grad(u))     scalar
!     z(10)    = f                  f                   scalar
!     z(11,12) = g                  g                   vector
!     z( 13)    = gammagradu        gamma*grad(u)       scalar
!     z( 14,15) = gammau            gamma*u             vector
!
!  local variables:
      integer (I4B) nff, nffmax, dof, errcode, i, inat, jnat, polyorder
      real (DP) :: xs, ys
      real (DP), allocatable ::  xsi(:), gxsi(:,:)
      complex (DPC) :: nu(2,2,nnat,nnat), gamma(2,nnat,nnat), beta(2,nnat,nnat)
      complex (DPC) :: alpha(nnat,nnat), em(nnat,nnat), f1(nnat), f2(2,nnat)
      complex (DPC) :: temp(2)
      complex (DPC), allocatable :: dgxsi(:)
!  for passing to shapefunction and dgshape:
      real (DP) :: xn1(3), yn1(3)
      complex (DPC) :: nu1(2,2)
!
!
      nffmax=0
      do inat=1,nnat
        nff=size(eg(elem,inat)%d)
        if (nff .gt. nffmax) then
          nffmax=nff
        end if
      end do
      
      z=0._DPC
      if (present (u)) then
        u=0._DPC
      end if
      if (present (alphau)) then
        alphau=0._DPC
      end if
      if (present (gammau)) then
        gammau=0._DPC
      end if
      if (present (gradu)) then
        gradu=0._DPC
      end if
      if (present (betagradu)) then
        betagradu=0._DPC
      end if
      if (present (gammagradu)) then
        gammagradu=0._DPC
      end if
      if (present (nugradu)) then
        nugradu=0._DPC
      end if
      if (present (dgradu)) then
        dgradu=0._DPC
      end if
      if (present (f)) then
        f=0._DPC
      end if
      if (present (g)) then
        g=0._DPC
      end if
!
!  if 'soln' is present, 'x' must be re-initialized by only one thread, one
!  single time during the call of 'field' for the current element
!  soln is only being passed by Post/fieldquantity.f90 
!
      if (present(soln)) then
        if (associated(x)) then
          deallocate(x)
        end if
        allocate(x(ndof))
        x = soln
      end if
!
!  get pde coefficients at the integration points
       call lam2xy(lambda,elem,xs,ys,xn,yn,e)
       call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2,em)
!
!
      polyorder=maxval(ep(elem,:))
      if (typ(2)) then
        allocate(xsi(nffmax))
        allocate(gxsi(nffmax,2))
        xn1 = xn(e(:,elem))
        yn1 = yn(e(:,elem))
        call shapefunction(lambda,xn1,yn1,1,polyorder,nffmax,&
                           .true.,xsi,gxsi,errcode)
      else if (typ(1)) then
        allocate(xsi(nffmax))
! to avoid runtime error with oneAPI
        allocate(gxsi(nffmax,2))
        xn1 = xn(e(:,elem))
        yn1 = yn(e(:,elem))
        call shapefunction(lambda,xn1,yn1,1,polyorder,nffmax,&
                           .false.,xsi,gxsi,errcode)
        deallocate(gxsi)                     
      end if
      if (typ(3)) then
        allocate(dgxsi(nffmax))
      end if
!  if construct for choosing output
!  u, alpha*u and 2nd term for determining the energy density will be returned
      if(typ(1)) then
        do inat=1,nnat
          nff=size(eg(elem,inat)%d)
          do i=1,nff
            dof=eg(elem,inat)%d(i)
            select case(dof)
            case (1:)
              z(1,inat) = z(1,inat) + xsi(i)*x(dof)
            case (:-1)
              z(1,inat) = z(1,inat) - xsi(i)*x(-dof)
            end select
          end do
        end do
        if (present (u)) then
          u(1:nnat) = z(1,1:nnat)
        end if

!  z(2,inat) = sum(alphau(inat,:))
        z(2,:) = matmul(alpha,z(1,:))
        if (present (alphau)) then
          alphau(:) = z(2,:)
        end if

!  gamma*u
        z(14,:) = matmul(gamma(1,:,:),z(1,:))
        z(15,:) = matmul(gamma(2,:,:),z(1,:))
        if (present (gammau)) then
          gammau(1:2,1:nnat) = z(14:15,1:nnat)
        end if

        z(3,:) = z(1,:)*conjg(z(2,:))
      end if
!
!  grad(u) and beta*grad(u), nu*grad(u) and curl will be returned
!  if there should be a vector, the first entry is in x- and the
!  second in y-direction
      if(typ(2)) then
        do inat=1,nnat
          nff=size(eg(elem,inat)%d)
          do i=1,nff
            dof=eg(elem,inat)%d(i)
            select case(dof)
            case (1:)
              z(4,inat) = z(4,inat) + gxsi(i,1)*x(dof)
              z(5,inat) = z(5,inat) + gxsi(i,2)*x(dof)
            case (:-1)
              z(4,inat) = z(4,inat) - gxsi(i,1)*x(-dof)
              z(5,inat) = z(5,inat) - gxsi(i,2)*x(-dof)
            end select
          end do
!  gradu(:,inat) is a vector so it has x and y components  denoted by indices 1 and 2
          if (present (gradu)) then
            gradu(1,1:nnat) = z(4,1:nnat)
            gradu(2,1:nnat) = z(5,1:nnat)
          end if
        end do
        z( 6,:) = 0._DPC
        z( 7,:) = 0._DPC
        z( 8,:) = 0._DPC
        z(13,:) = 0._DPC
        do inat=1,nnat
          do jnat=1,nnat
            z( 6,inat) = z( 6,inat) + sum(beta(1:2,inat,jnat)*z(4:5,jnat))
            z(13,inat) = z(13,inat) + sum(gamma(1:2,inat,jnat)*z(4:5,jnat))
            temp = matmul(nu(:,:,inat,jnat),z(4:5,jnat))
            z(7,inat) = z(7,inat) + temp(1)
            z(8,inat) = z(8,inat) + temp(2)
          end do
        end do
!  nugradu(:,jnat) is a vector so it has x and y components  denoted by indices 1 and 2
          if (present (nugradu)) then
            nugradu(1,1:nnat) = z(7,1:nnat)
            nugradu(2,1:nnat) = z(8,1:nnat)
          end if
!  betagradu(jnat) is a scalar
        if (present (betagradu)) then
          betagradu(1:nnat) = z(6,1:nnat)
        end if
!  gammagradu(jnat) is a scalar
        if (present (gammagradu)) then
          gammagradu(1:nnat) = z(13,1:nnat)
        end if


      end if
!
!  dgrad(u) will be returned (already contains nu!)
      if(typ(3)) then
        do inat=1,nnat
          do jnat=1,nnat
            polyorder=ep(elem,jnat)
            nff=size(eg(elem,jnat)%d)
            nu1 = nu(:,:,inat,jnat)
            xn1 = xn(e(:,elem))
            yn1 = yn(e(:,elem))
            call dgshape(lambda,nu1,xn1,yn1,1,polyorder,nff,dgxsi,errcode)
            do i=1,nff
              dof=eg(elem,jnat)%d(i)
              select case(dof)
              case (1:)
                z(9,inat) = z(9,inat) + dgxsi(i)*x(dof)
              case (:-1)
                z(9,inat) = z(9,inat) - dgxsi(i)*x(-dof)
              end select
            end do
          end do
        end do
!  dgradu(jnat) is a scalar
        if (present (dgradu)) then
          dgradu(:) = z(9,:)
        end if
      end if
!
!  f1 will be returned
      if(typ(4)) then
        z(10,:) = f1(:)
        if (present (f)) then
          f(1:nnat) = f1(1:nnat)
        end if
      end if
!
!  f2 will be returned
      if(typ(5)) then
        z(11,:) = f2(1,:)
        z(12,:) = f2(2,:)
        if (present (f)) then
          g(1,1:nnat) = f2(1,1:nnat)
          g(2,1:nnat) = f2(2,1:nnat)
        end if
      end if
!
      if (allocated(xsi)) deallocate(xsi)
      if (allocated(gxsi)) deallocate(gxsi)
      if (allocated(dgxsi)) deallocate(dgxsi)
      return
      end subroutine field


      pure subroutine fieldu(elem, lambda, u, gradu)
      use feminterface, only: shapefunction
      use femtypes
      use globalvariables, only: e, eg, x, ndof, xn, yn, nnat, ep
      implicit none
      integer (I4B) elem
      real (DP) lambda(3)
      complex (DPC) :: u(:)
      complex (DPC), optional :: gradu(:,:)
      intent (in) ::  elem, lambda
      intent (out) ::  u, gradu
!------------------------------------------------------------------------------
!  minimalistic version of the field routine
!  no call to pdecoeff
!------------------------------------------------------------------------------
!  local variables:
      integer (I4B) nff, nffmax, dof, errcode, i, inat, polyorder
      real (DP), allocatable ::  xsi(:), gxsi(:,:)
!  for passing to shapefunction and dgshape:
      real (DP) :: xn1(3), yn1(3)
!
      nffmax=0
      do inat=1,nnat
        nff=size(eg(elem,inat)%d)
        if (nff .gt. nffmax) then
          nffmax=nff
        end if
      end do
!  Initialize u
      u=0._DPC
!
      polyorder=maxval(ep(elem,:))
!
      xn1 = xn(e(:,elem))
      yn1 = yn(e(:,elem))
      if (present(gradu)) then
        allocate(xsi(nffmax), gxsi(nffmax,2))
        call shapefunction(lambda,xn1,yn1,1,polyorder,nffmax,&
                           .true.,xsi,gxsi,errcode)
      else
        allocate(xsi(nffmax))
        call shapefunction(lambda,xn1,yn1,1,polyorder,nffmax,&
                           .false.,xsi,gxsi,errcode)
      end if

!  u will be returned
      do inat=1,nnat
        nff=size(eg(elem,inat)%d)
        do i=1,nff
          dof=eg(elem,inat)%d(i)
          select case(dof)
          case (1:)
            u(inat) = u(inat) +  xsi(i)*x(dof)
          case (:-1)
            u(inat) = u(inat) - xsi(i)*x(-dof)
          end select
        end do
      end do
      deallocate(xsi)
!
!  grad(u) will be returned
      if (present(gradu)) then
        gradu = 0._DPC
        do inat=1,nnat
          nff=size(eg(elem,inat)%d)
          do i=1,nff
            dof=eg(elem,inat)%d(i)
            select case(dof)
            case (1:)
              gradu(1,inat) = gradu(1,inat) + gxsi(i,1)*x(dof)
              gradu(2,inat) = gradu(2,inat) + gxsi(i,2)*x(dof)
            case (:-1)
              gradu(1,inat) = gradu(1,inat) - gxsi(i,1)*x(-dof)
              gradu(2,inat) = gradu(2,inat) - gxsi(i,2)*x(-dof)
            end select
          end do
        end do
        deallocate(gxsi)
      end if

      return
      end subroutine fieldu