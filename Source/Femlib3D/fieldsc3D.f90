      subroutine fieldsc3D(elem,lambda,typ,z,soln,u,alphau,gammau,gradu,betagradu,gammagradu,nugradu,dgradu,f,g)
      use feminterface3D, only: shapefunctionsc3D, pdecoeff3D, lam2xyz, dgshapesc3D, getgl
      use femtypes
      use globalvariables3D, only: nnat, nod, numdof, vgdof, vn, vp, x
      implicit none
      integer (I4B) elem
      real (DP) lambda(4)
      logical typ(5)
      complex (DPC) z(:,:)
      complex (DPC), optional :: u(:), alphau(:,:), gammau(:,:,:), dgradu(:,:), f(:), g(:,:)
      complex (DPC), optional :: betagradu(:,:), gammagradu(:,:), nugradu(:,:,:)
      complex (DPC), optional :: gradu(:,:)
      complex (DPC), optional :: soln(:)
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
!------------------------------------------------------------------------------
!    $Revision: 1.11 $
!    $Date: 2015/11/11 17:27:57 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Evaluate the potential, the gradient, the dgradient and products with the
!  corresponding pde coefficients.
!  (Remark: The material is assumed to be constant across the element.)
!  The subroutine will return all values connected with a set type in logical
!  vector typ.
!    - if typ(1) is .true. everything connected with u is returned
!      (u, alpha*u, gamma*u and second term for determining energy density)
!    - if typ(2) is .true. everything connected with grad(u) is returned
!      (grad(u), beta*grad(u), gamma*grad(u), nu*grad(u) and the curl))
!    - if typ(3) is .true. the 2nd gradient dgrad(u) will be returned
!      (it already contains nu like programmed in dgshape.f90 subroutine)
!    - if typ(4) is .true. the scalar source term f1 (also called f) is returned
!    - if typ(5) is .true. the sourcevector f2 (also called g) is returned
!
!  Input:
!     elem      element number
!     lambda    vector of barycentric coordinates
!     typ       case selector
!
!  Output:
!     z         value for typ (all dof values summed up)
!               (this variable is a vector and should be allocated by the
!               calling subroutine)
!     u         size is u(nnat)
!     alphau    size is alphau(nnat,nnat)
!     gammau    size is gammau(2,nnat,nnat)
!     gradu     size is gradu(2,nnat)
!     betagradu size is betagradu(nnat,nnat)
!     gammagradu size is gammagradu(nnat,nnat)
!     nugradu   size is nugradu(2,nnat,nnat)
!     dgradu    size is dgradu(nnat,nnat)
!     f         size is f(nnat)
!     g         size is g(2,nnat)
!     z( 1)       = u                  u                   scalar
!     z( 2)       = alphau             alpha*u             scalar
!     z( 3)       = ualphau            u*conj(alpha*u)     scalar
!     z( 4, 5, 6) = gradu              grad(u)             vector
!     z( 7)       = betagradu          beta*grad(u)        scalar
!     z( 8, 9,10) = nugradu            nu*grad(u)          vector
!     z(11)       = dgradu             div(nu*grad(u))     scalar           !  TO DO dgxsi is not available
!     z(12)       = f                  f                   scalar
!     z(13,14,15) = g                  g                   vector

!
!  local variables:
      integer (I4B) nff, nffmax, dof, errcode, i, inat, jnat, polyorder
      real (DP) :: xyzt(3), gl(3,4)
      real (DP), allocatable ::  xsi(:), gxsi(:,:)
      complex (DPC) :: nu(3,3,nnat,nnat), gamma(3,nnat,nnat), alpha(nnat,nnat)
      complex (DPC) :: beta(3,nnat,nnat), f1(nnat), f2(3,nnat)
      complex (DPC) :: temp(3)
      complex (DPC), allocatable :: dgxsi(:)
!  for passing to shapefunction and dgshape:
      real (DP) :: vert(3,4)
      complex (DPC) :: nu1(3,3)
!
!
      nffmax=0
      do inat=1,nnat
        nff=size(vgdof(elem,inat)%d)
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

      if (present(soln)) then
        if (associated(x)) then
          deallocate(x)
        end if
        allocate(x(numdof))
        x = soln
      end if
!
!  get pde coefficients at the integration points
       call lam2xyz(lambda,elem,xyzt,nod,vn)
       call pdecoeff3D(elem,xyzt,nu,alpha,beta,f1,f2,gamma)
!
      polyorder=maxval(vp(elem,:))
      vert(1:3,1:4) = nod(1:3,vn(1:4,elem))
      if (typ(2)) then
        allocate(xsi(nffmax))
        allocate(gxsi(3,nffmax))
        call shapefunctionsc3D(lambda,vert,1,polyorder,nffmax,&
                           .true.,xsi,gxsi,errcode)
      else if (typ(1)) then
        allocate(xsi(nffmax))
! Dummy allocate to avoid runtime error in oneAPI
        allocate( gxsi(3,nffmax) )
        call shapefunctionsc3D(lambda,vert,1,polyorder,nffmax,&
                           .false.,xsi,gxsi,errcode)
        deallocate(gxsi)
      end if
      if (typ(3)) then
!  TO DO dgxsi is not available
        allocate(dgxsi(nffmax))
      end if
!
!  if construct for choosing output
!  u, alpha*u and 2nd term for determining the energy density will be returned
      if(typ(1)) then
        do inat=1,nnat
          nff=size(vgdof(elem,inat)%d)
          do i=1,nff
            dof=vgdof(elem,inat)%d(i)
            if (dof .eq. 0) cycle
            z(1,inat) = z(1,inat) + xsi(i)*x(dof)
          end do
        end do
        if (present (u)) then
          u(1:nnat) = z(1,1:nnat)
        end if
!
        z(2,:) = matmul(alpha,z(1,:))
        do inat=1,nnat
          do jnat=1,nnat
            if (present (alphau)) then
              alphau(inat,jnat) = alpha(inat,jnat)*u(jnat)
            end if
            if (present (gammau)) then
              gammau(1,inat,jnat) = gamma(1,inat,jnat)*u(jnat)
              gammau(2,inat,jnat) = gamma(2,inat,jnat)*u(jnat)
              gammau(3,inat,jnat) = gamma(3,inat,jnat)*u(jnat)
            end if
          end do
        end do
        z(3,:) = z(1,:)*conjg(z(2,:))
      end if
!
!  grad(u) and beta*grad(u), nu*grad(u) and curl will be returned
!  if there should be a vector, the first entry is in x- and the
!  second in y-direction
      if(typ(2)) then
        do inat=1,nnat
          nff=size(vgdof(elem,inat)%d)
          do i=1,nff
            dof=vgdof(elem,inat)%d(i)
            if (dof .eq. 0) cycle
            z(4,inat) = z(4,inat) + gxsi(1,i)*x(dof)
            z(5,inat) = z(5,inat) + gxsi(2,i)*x(dof)
            z(6,inat) = z(6,inat) + gxsi(3,i)*x(dof)
!  gradu(:,inat) is a vector so it has x and y components denoted by indices 1 and 2
          end do
        end do
        if (present (gradu)) then
          gradu(1,1:nnat) = z(4,1:nnat)
          gradu(2,1:nnat) = z(5,1:nnat)
          gradu(3,1:nnat) = z(6,1:nnat)
        end if

        do inat=1,nnat
          do jnat=1,nnat
            z(7,inat) = z(7,inat) + sum(beta(1:3,inat,jnat)*z(4:6,jnat))
!  betagradu(inat,jnat) is a scalar
            if (present (betagradu)) then
!              betagradu(inat,jnat) = sum(beta(1:3,inat,jnat)*gradu(1:3,jnat))
              betagradu(inat,jnat) = sum(beta(1:3,inat,jnat)*z(4:6,jnat))
            end if
!  gammagradu(inat,jnat) is a scalar
            if (present (gammagradu)) then
!              gammagradu(inat,jnat) = sum(gamma(1:3,inat,jnat)*gradu(1:3,jnat))
              gammagradu(inat,jnat) = sum(gamma(1:3,inat,jnat)*z(4:6,jnat))
            end if
!  z(4:6) is gradu
            temp = matmul(nu(:,:,inat,jnat),z(4:6,jnat))
            z( 8,inat) = z( 8,inat) + temp(1)
            z( 9,inat) = z( 9,inat) + temp(2)
            z(10,inat) = z(10,inat) + temp(3)
            if (present (nugradu)) then
!              temp = matmul(nu(:,:,inat,jnat),gradu(1:3,jnat))
!  nugradu(:,inat,jnat) is a vector so it has x and y and z components  denoted by indices 1 and 2 and 3
              nugradu(1,inat,jnat) = temp(1)
              nugradu(2,inat,jnat) = temp(2)
              nugradu(3,inat,jnat) = temp(3)
            end if
          end do
        end do
      end if
!
!  dgrad(u) will be returned (already contains nu!)
      if(typ(3)) then
        do inat=1,nnat
          do jnat=1,nnat
            polyorder=vp(elem,jnat)
            nff=size(vgdof(elem,jnat)%d)
            nu1 = nu(:,:,inat,jnat)
            call getgl(vert,gl)
!  TO DO dgxsi is not available
            call dgshapesc3D(lambda, gl, nu1, 1, polyorder, dgxsi, errcode)
            do i=1,nff
              dof=vgdof(elem,jnat)%d(i)
              if (dof .eq. 0) cycle
              z(11,inat) = z(11,inat) + dgxsi(i)*x(dof)
!  dgradu(inat,jnat) is a scalar
              if (present (dgradu)) then
                dgradu(inat,jnat) = dgradu(inat,jnat) + dgxsi(i)*x(dof)
              end if
            end do
          end do
        end do
      end if
!
!  f1 will be returned
      if(typ(4)) then
        z(12,:) = f1(:)
        if (present (f)) then
          do inat=1,nnat
            f(inat) = f1(inat)
          end do
        end if
      end if
!
!  f2 will be returned
      if(typ(5)) then
        z(13,:) = f2(1,:)
        z(14,:) = f2(2,:)
        z(15,:) = f2(3,:)
        if (present (f)) then
          do inat=1,nnat
            g(1,inat) = f2(1,inat)
            g(2,inat) = f2(2,inat)
            g(3,inat) = f2(3,inat)
          end do
        end if
      end if
!
      if (allocated(xsi)) deallocate(xsi)
      if (allocated(gxsi)) deallocate(gxsi)
      if (allocated(dgxsi)) deallocate(dgxsi)
      return
      end subroutine fieldsc3D



      subroutine fieldsc3D_simple(elem,lambda,u,gu)
      use feminterface3D, only: shapefunctionsc3D
      use femtypes
      use globalvariables3D, only: nnat, nod, numdof, vgdof, vn, vp, x
      implicit none
      integer (I4B)           :: elem
      real (DP)               :: lambda(4)
      complex (DPC)           :: u(:)
      complex (DPC), optional :: gu(:,:)
      intent (in)             :: elem, lambda
      intent (out)            :: u
!------------------------------------------------------------------------------
!  local variables:
      integer (I4B) nff, nffmax, dof, errcode, i, inat, polyorder
      real (DP), allocatable ::  xsi(:), gxsi(:,:)
!  for passing to shapefunction and dgshape:
      real (DP) :: vert(3,4)
!
      nffmax=0
      do inat=1,nnat
        nff=size(vgdof(elem,inat)%d)
        if (nff .gt. nffmax) then
          nffmax=nff
        end if
      end do

      polyorder=maxval(vp(elem,:))
      vert(1:3,1:4) = nod(1:3,vn(1:4,elem))
      
      
      if (present(gu)) then
        
        allocate(xsi(nffmax),gxsi(3,nffmax))
        call shapefunctionsc3D(lambda,vert,1,polyorder,nffmax,&
                         .true.,xsi,gxsi,errcode)
        u=0.
        gu=0.
        do inat=1,nnat
          nff=size(vgdof(elem,inat)%d)
          do i=1,nff
            dof=vgdof(elem,inat)%d(i)
            if (dof .eq. 0) cycle
            u(inat) = u(inat) + xsi(i)*x(dof)
            gu(1,inat) = gu(1,inat) + gxsi(1,i)*x(dof)
            gu(2,inat) = gu(2,inat) + gxsi(2,i)*x(dof)
            gu(3,inat) = gu(3,inat) + gxsi(3,i)*x(dof)
          end do
        end do
        if (allocated(xsi)) deallocate(xsi)
        if (allocated(gxsi)) deallocate(gxsi)
      
      else
        
        allocate(xsi(nffmax))
! Dummy allocate to avoid runtime error in oneAPI
        allocate( gxsi(3,nffmax) )
        call shapefunctionsc3D(lambda,vert,1,polyorder,nffmax,&
                         .false.,xsi,gxsi,errcode)

        u=0.
        do inat=1,nnat
          nff=size(vgdof(elem,inat)%d)
          do i=1,nff
            dof=vgdof(elem,inat)%d(i)
            if (dof .eq. 0) cycle
            u(inat) = u(inat) + xsi(i)*x(dof)
          end do
        end do
        if (allocated(xsi)) deallocate(xsi)
        deallocate(gxsi)
      end if
!      if (elem==81) print *, 'INFIELD E81 ',real(u(1)),'X121 ',real(x(121))
      return
    end subroutine fieldsc3D_simple
      
      subroutine fieldsc3D_simple_aux(elem,lambda,u,gu)
      use feminterface3D, only: shapefunctionsc3D
      use femtypes
      use globalvariables3D, only: nnat_aux, nod, numdof_aux, vgdof_aux, vn, vp_aux, x_aux
      implicit none
      integer (I4B)           :: elem
      real (DP)               :: lambda(4)
      complex (DPC)           :: u(:)
      complex (DPC), optional :: gu(:,:)
      intent (in)             :: elem, lambda
      intent (out)            :: u
!------------------------------------------------------------------------------
!  local variables:
      integer (I4B) nff, nffmax, dof, errcode, i, inat, polyorder
      real (DP), allocatable ::  xsi(:), gxsi(:,:)
!  for passing to shapefunction and dgshape:
      real (DP) :: vert(3,4)
!
      nffmax=0
      do inat=1,nnat_aux
        nff=size(vgdof_aux(elem,inat)%d)
        if (nff .gt. nffmax) then
          nffmax=nff
        end if
      end do

      polyorder=maxval(vp_aux(elem,:))
      vert(1:3,1:4) = nod(1:3,vn(1:4,elem))
      
      
      if (present(gu)) then
        
        allocate(xsi(nffmax),gxsi(3,nffmax))
        call shapefunctionsc3D(lambda,vert,1,polyorder,nffmax,&
                         .true.,xsi,gxsi,errcode)
        u=0.
        gu=0.
        do inat=1,nnat_aux
          nff=size(vgdof_aux(elem,inat)%d)
          do i=1,nff
            dof=vgdof_aux(elem,inat)%d(i)
            if (dof .eq. 0) cycle
            u(inat) = u(inat) + xsi(i)*x_aux(dof)
            gu(1,inat) = gu(1,inat) + gxsi(1,i)*x_aux(dof)
            gu(2,inat) = gu(2,inat) + gxsi(2,i)*x_aux(dof)
            gu(3,inat) = gu(3,inat) + gxsi(3,i)*x_aux(dof)
          end do
        end do
        if (allocated(xsi)) deallocate(xsi)
        if (allocated(gxsi)) deallocate(gxsi)
      
      else
        
        allocate(xsi(nffmax))
! Dummy allocate to avoid runtime error in oneAPI
        allocate( gxsi(3,nffmax) )
        call shapefunctionsc3D(lambda,vert,1,polyorder,nffmax,&
                         .false.,xsi,gxsi,errcode)

        u=0.
        do inat=1,nnat_aux
          nff=size(vgdof_aux(elem,inat)%d)
          do i=1,nff
            dof=vgdof_aux(elem,inat)%d(i)
            if (dof .eq. 0) cycle
            u(inat) = u(inat) + xsi(i)*x_aux(dof)
          end do
        end do
        if (allocated(xsi)) deallocate(xsi)
        deallocate(gxsi)
      end if


      

      return
      end subroutine fieldsc3D_simple_aux
      
      
      
