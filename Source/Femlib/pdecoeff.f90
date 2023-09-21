      subroutine pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2,em)
      use feminterface, only: fetchmatparameters, calctensor, invertmat3, &
                          & multtensor, transftensor, trans_matrix, transten_r2, xy2lam, fieldu, getsetting
                    
      use femtypes
      use globalvariables, only: pi, eps0, mu0, omega, gravc, Elch, Kboltz, Avogadro, physics, geb, matzif, xn, yn, e, nnat
      use matconstants, only: numparam
      use mpmodule, only: matmodify, sourcemodify, multiphysics, trmat
      implicit none
      integer (I4B) elem
      real(DP) xs, ys
      complex (DPC) :: nu(:,:,:,:), gamma(:,:,:)
      complex (DPC) :: beta(:,:,:), alpha(:,:), f1(:), f2(:,:)
      complex (DPC), optional :: em(:,:)
      intent (in) :: elem, xs, ys
      intent (out) :: nu, gamma, beta, alpha, f1, f2, em
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
!    $Revision: 1.47 $
!    $Date: 2015/06/02 09:29:18 $
!    $Author: m_kasper $
!
!  deliver the material coefficients in the form of the general PDE:
!  
!        -div(nu*grad(u)) - div(gamma*u) + (beta*grad(u)) + alpha*u = f1 + div(f2)
!  
!  Input:
!            elem     number of the element
!            x        actual solution vector [not used yet]
!            xn,yn    coordinates of triangle nodes [not used yet]
!            e        element information (nodes of the element) [not used yet]
!            matzif   material name of the input regions
!            geb      region index of the elements 
!            xs, ys   location at which the material coefficients have to be evaluated
!                     [not used yet]
!  Output:
!            nu       rank 2  tensor of the diffusion term (nu(1,1)=nuxx, nu(1,2)=nuxy, ...)
!            gamma    rank 1 tensor
!            beta     vector of the convective term beta(1)=betax
!            alpha    constant
!            em       constant for the second derivative
!            f1       right hand side constant
!            f2       right hand side source vector
!
!  local variables
      integer (I4B) matindex, pml
      integer (I4B) i, j, k, l
      real(DP), allocatable :: list(:)
      real (DP) :: angle, phase, losst, scalar, xys(2)
      real (DP) :: jmx, jmy, gam(2,2), v11, v22
      real (DP) :: Lame1, Lame2
      real (DP) :: TransfTens(3,3)
!  3D Tensors
      real (DP) :: CTensNoSym3D(3,3,3,3), ETensNoSym3D(3,3,3), epsrTensNoSym3D(3,3)
      real (DP) :: LamTensNoSym3D(3,3), ATensNoSym3D(3,3), PTensNoSym3D(3)
      real (DP) :: MTensNoSym3D(3,3)
!  2D Tensors
      real (DP) :: CTensNoSym(2,2,2,2)
      real (DP) :: ATensNoSym(2,2)
      real (DP) :: MTensNoSym(2,2)
      !  For non-linear problems with material parameters that depend on the actual solution
      real (DP) lambda(3)
      complex (DPC) u(nnat), gradu(2,nnat)
!
      complex (DPC) :: var11, var22, zcomp, source1, source2, matten(3,3)
      logical :: invert
      character (len=6) :: nonlinear
      real (DP), parameter :: rad = pi/180._DP

!------------------------------------------------------------------------------
!  new stuff for 3D anisotropic materials
      real (DP) :: vpa1(3), vpa2(3), vpml(3)
      real (DP) :: trmat2d(3,3) = 0._DP
      real (DP) :: delta_eps(3), delta_mu(3), losst_eps(3), losst_mu(3)
      complex (DPC) :: epsr(3,3), mur(3,3), source(3)
      complex (DPC) :: pmlten(3,3) = 0._DPC
!  Used for piezopyroelectricity at the moment
      complex (DPC) :: rho
      complex (DPC) :: pmldir(2) = 0._DPC, pmlprod = 0._DPC
!------------------------------------------------------------------------------


!
!  determine material index for given element in region geb
      matindex=matzif(geb(elem))

!  read parameters from internal list and assign them to local variables
      allocate (list(numparam))
      xys=(/xs,ys/)
      call fetchmatparameters(list,matindex,xys,elem)
!  possible choices for physics can be:
!  ACOUSTIC, FLOW_PROFILE, HEATTR, MAGNETICVP, TEWAVE and TMWAVE.
      if ((physics .eq. 'TEWAVE') .or. (physics .eq. 'TMWAVE')) then
!  assign PML values
        pml  = int(list(1))
        vpml = (/list(2),list(3),list(4)/)
!  assign principal axes (relative to simulation coordinate system)
        vpa1 = (/list(5),list(6),list(7)/)
        vpa2 = (/list(8),list(9),list(10)/)
!  compute transformation matrix from crystal to simulation coordinate system
        trmat = trans_matrix(vpa1,vpa2)
!  assign relative permittivity values (in crystal coordinate system)
        losst_eps = (/list(14),list(15),list(16)/)
        delta_eps = -1._DP * atan(losst_eps)
        epsr      = 0._DPC
        epsr(1,1) = list(11) * exp(cmplx(0._DP, delta_eps(1), DPC))
        epsr(2,2) = list(12) * exp(cmplx(0._DP, delta_eps(2), DPC))
        epsr(3,3) = list(13) * exp(cmplx(0._DP, delta_eps(3), DPC))
!  modify for multiphysics in crystal coordinate system
        if (multiphysics) then
          call matmodify(list,xs,ys,epsr)
        end if
!  assign relative permeability values (in crystal coordinate system)
        losst_mu = (/list(20),list(21),list(22)/)
        delta_mu = -1._DP * atan(losst_mu)
        mur      = 0._DPC
        mur(1,1) = list(17) * exp(cmplx(0._DP, delta_mu(1), DPC))
        mur(2,2) = list(18) * exp(cmplx(0._DP, delta_mu(2), DPC))
        mur(3,3) = list(19) * exp(cmplx(0._DP, delta_mu(3), DPC))
!  execute transformation of material rank 2 tensors from crystal to simulaion
!  coordinate system
        epsr  = transten_r2(trmat,epsr)
        mur   = transten_r2(trmat,mur)
!  modify for pml
        if (pml .eq. 1) then
!  assign pml values for attenuation in x-direction (pml is set to loss tangent of 1)
          pmlten(1,1) = exp(cmplx(0._DP, atan(1._DP), DPC))
          pmlten(2,2) = exp(cmplx(0._DP,-atan(1._DP), DPC))
          pmlten(3,3) = exp(cmplx(0._DP,-atan(1._DP), DPC))
!  compute angle between x-axis and pml-vector
          angle = acos(dot_product(vpml,(/1,0,0/)) / sqrt(vpml(1)**2 + vpml(2)**2 + vpml(3)**2))
!  compute 2d transformation matrix for pml and ...
          trmat2d(1,:) = (/ cos(angle),sin(angle),0._DP/)
          trmat2d(2,:) = (/-sin(angle),cos(angle),0._DP/)
          trmat2d(3,:) = (/ 0._DP,     0._DP,     1._DP/)
!  execute transformation/rotation of pml tensor
          pmlten = transten_r2(trmat2d,pmlten)
!  modify material tensors by elementwise multiplication with pml tensor
          epsr = epsr * pmlten
          mur  = mur  * pmlten
        end if
!  assign current density values
        source = (/list(23),list(24),list(25)/)
        phase  = list(26)
      end if
!
!
!  select case construct to assign pde coefficients
      select case (physics)
      case ('ACOUSTIC')
        pml    = int(list(1))
        angle  = rad*list(7)
        losst  = list(4)
!  check for PML
        if (pml .eq. 1) then
          var11 = (1._DP/list(2))*exp(cmplx(0._DP, atan(losst),DPC))
          var22 = (1._DP/list(2))*exp(cmplx(0._DP,-atan(losst),DPC))
          zcomp = (1._DP/list(2))*exp(cmplx(0._DP,-atan(losst),DPC))
        else
          var11 = cmplx(1._DP/list(2),0._DP,DPC)
          var22 = cmplx(1._DP/list(2),0._DP,DPC)
          zcomp = cmplx(1._DP/list(2),0._DP,DPC)
        end if
!  calculate tensor and assign all other values
        invert = .false.
        call calctensor (var11, var22, angle, invert, nu(:,:,1,1))
        alpha(1,1) = -omega**2*zcomp/list(3)**2
        beta(:,1,1) = 0._DPC
        gamma(:,1,1) = 0._DPC
!  SETTINGS AND SOURCE TERM NOT VERIFIED!!!
        f1(1)      = 0._DPC
        f2(1,1)   = -cmplx(list(5)/list(2),0._DP,DPC)
        f2(2,1)   = -cmplx(list(6)/list(2),0._DP,DPC)
        deallocate (list)
!
!
      case ('ELECTROSTATICS')
!  var11 and var22 contain eps_r in crystal axes
        var11 = list(1)
        var22 = list(2)
        angle = rad*list(3)
!  calculate material tensor for diffusion term and assign result to nu
        invert = .false.
        call calctensor (var11, var22, angle, invert, nu(:,:,1,1))
!  assign scalar source term as charge density divided by eps_0
        f1(1) = cmplx(list(4)/eps0,0._DP,DPC)
!  assign alpha, beta and f2
        alpha(1,1) = 0._DPC
        beta(:,1,1) = 0._DPC
        gamma(:,1,1) = 0._DPC
        f2(:,1) = 0._DPC
        deallocate(list)       
!
!
      case ('FLOW_PROFILE')
        nu(1,1,1,1)  = cmplx(list(1)/list(2),0._DP,DPC)
        nu(2,2,1,1)  = nu(1,1,1,1)
        nu(1,2,1,1)  = 0._DPC
        nu(2,1,1,1)  = nu(1,2,1,1)
        beta(:,1,1)  = 0._DPC
        gamma(:,1,1) = 0._DPC
        alpha(1,1)   = cmplx(0._DP,omega,DPC)
        f1(1)        = cmplx(-list(3)/list(2),0._DP,DPC)
        f2(:,1)      = 0._DPC
        deallocate(list)
!
!
      case ('FLUID')
!       list(1) = mu or eta
!       list(2) = rho
!       If non-dimensional is used, nu = 1/Re and rho = 1
!     
        nu(1,1,1,1)  = cmplx(list(1)/list(2),0._DP,DPC)  
        nu(2,2,1,1)  = nu(1,1,1,1)
        nu(1,2,1,1)  = 0._DPC
        nu(2,1,1,1)  = nu(1,2,1,1)
        beta(:,1,1)  = 0._DPC
        gamma(:,1,1) = 0._DPC
        alpha(1,1)   = cmplx(0._DP,omega,DPC)
        f1(1)        = cmplx(list(3)/list(2),0._DP,DPC)
        f2(1,1)      = cmplx(list(4),0._DP,DPC)
        f2(2,1)      = cmplx(list(5),0._DP,DPC)
        deallocate(list)
!
!
      case ('FLUIDINCOMPR')
!       list(1) = mu
!       list(2) = rho
!
        call xy2lam(xs,ys,elem,lambda,xn,yn,e)
        call fieldu(elem,lambda,u)
!
        nu(:,:,:,:)  = 0._DPC
        nu(1,1,1,1)  = cmplx(list(1),0._DP,DPC)
        nu(2,2,1,1)  = cmplx(list(1),0._DP,DPC)
        nu(1,1,2,2)  = cmplx(list(1),0._DP,DPC)
        nu(2,2,2,2)  = cmplx(list(1),0._DP,DPC)
!
        beta(:,:,:)  = 0._DPC
        beta(1,1,1)  = list(2)*u(1)
        beta(2,1,1)  = list(2)*u(2)
        beta(1,1,3)  = cmplx(1._DP,0._DP,DPC)
        beta(1,2,2)  = list(2)*u(1)
        beta(2,2,2)  = list(2)*u(2)
        beta(2,2,3)  = cmplx(1._DP,0._DP,DPC)
        beta(1,3,1)  = cmplx(1._DP,0._DP,DPC)
        beta(2,3,2)  = cmplx(1._DP,0._DP,DPC)
!
        gamma(:,:,:) = 0._DPC
        alpha(:,:)   = 0._DPC
        f1(:)        = 0._DPC
        f2(:,:)      = 0._DPC
        deallocate(list)
!
!
      case ('FLUIDELINCOMPR')
!       list(1) = mu
!       list(2) = rho
!       list(3) = epsr11
!       list(4) = epsr22
!
        call xy2lam(xs,ys,elem,lambda,xn,yn,e)
        call fieldu(elem,lambda,u)
!
        nu(:,:,:,:)  = 0._DPC
        ! v_x <- v_x
        nu(1,1,1,1)  = cmplx(list(1),0._DP,DPC)
        nu(2,2,1,1)  = cmplx(list(1),0._DP,DPC)
        ! v_y <- v_y
        nu(1,1,2,2)  = cmplx(list(1),0._DP,DPC)
        nu(2,2,2,2)  = cmplx(list(1),0._DP,DPC)
        ! El. Pot. <- El. Pot.
        nu(1,1,4,4)  = cmplx(eps0*list(3),0._DP,DPC)
        nu(2,2,4,4)  = cmplx(eps0*list(4),0._DP,DPC)
!
        beta(:,:,:)  = 0._DPC
        ! v_x <- v_x
        beta(1,1,1)  = list(2)*u(1)
        beta(2,1,1)  = list(2)*u(2)
        ! v_x <- P
        beta(1,1,3)  = cmplx(1._DP,0._DP,DPC)
        ! v_x <- El. Pot.
        beta(1,1,4)  = cmplx(1._DP,0._DP,DPC)
        ! v_y <- v_y
        beta(1,2,2)  = list(2)*u(1)
        beta(2,2,2)  = list(2)*u(2)
        ! v_y <- P
        beta(2,2,3)  = cmplx(1._DP,0._DP,DPC)
        ! v_y <- El. Pot.
        beta(2,2,4)  = cmplx(1._DP,0._DP,DPC)
        ! P <- v_x
        beta(1,3,1)  = cmplx(1._DP,0._DP,DPC)
        ! P <- v_y
        beta(2,3,2)  = cmplx(1._DP,0._DP,DPC)
!
        gamma(:,:,:) = 0._DPC
        alpha(:,:)   = 0._DPC
        f1(:)        = 0._DPC
        f2(:,:)      = 0._DPC
        deallocate(list)
!
!
      case ('FLUIDELMASSCON')
!       list(1) = mu
!       list(2) = rho
!       list(7) = zeta
!       list(8) = T_F
!
        call xy2lam(xs,ys,elem,lambda,xn,yn,e)
        call fieldu(elem,lambda,u)
!
        nu(:,:,:,:)  = 0._DPC
        nu(1,1,1,1)  = cmplx(list(1),0._DP,DPC)
        nu(2,2,1,1)  = cmplx(list(1),0._DP,DPC)
        nu(1,1,2,2)  = cmplx(list(1),0._DP,DPC)
        nu(2,2,2,2)  = cmplx(list(1),0._DP,DPC)
        nu(1,1,4,4)  = cmplx(eps0*list(3),0._DP,DPC)
        nu(2,2,4,4)  = cmplx(eps0*list(4),0._DP,DPC)
        nu(1,1,5,4)  = u(5) * (list(5)*list(7)*Elch) / (Kboltz*list(8))
        nu(2,2,5,4)  = u(5) * (list(6)*list(7)*Elch) / (Kboltz*list(8))
        nu(1,1,5,5)  = cmplx(list(5),0._DP,DPC)
        nu(2,2,5,5)  = cmplx(list(6),0._DP,DPC)
!
        beta(:,:,:)  = 0._DPC
        beta(1,1,1)  = list(2)*u(1)
        beta(2,1,1)  = list(2)*u(2)
        beta(1,1,3)  = cmplx(1._DP,0._DP,DPC)
        beta(1,1,4)  = u(5) * Avogadro*list(7)*Elch
        beta(1,2,2)  = list(2)*u(1)
        beta(2,2,2)  = list(2)*u(2)
        beta(2,2,3)  = cmplx(1._DP,0._DP,DPC)
        beta(2,2,4)  = u(5) * Avogadro*list(7)*Elch
        beta(1,3,1)  = cmplx(1._DP,0._DP,DPC)
        beta(2,3,2)  = cmplx(1._DP,0._DP,DPC)
!
        gamma(:,:,:) = 0._DPC
        gamma(1,5,5) = -u(1)
        gamma(2,5,5) = -u(2)
!
        alpha(:,:)   = 0._DPC
        alpha(1,1)   = cmplx(0._DP,list(2)*omega,DPC)
        alpha(2,2)   = cmplx(0._DP,list(2)*omega,DPC)
        alpha(5,5)   = cmplx(0._DP,omega,DPC)
!
        f1(:)        = 0._DPC
        f2(:,:)      = 0._DPC
        deallocate(list)
!
!
      case ('FLUIDELMASSTEP')
!       list(1) = mu
!       list(2) = rho
!       list(7) = zeta
!       list(8) = lam11
!       list(9) = lam22
!       list(10)= c_p
!       list(11)= p
!
        call xy2lam(xs,ys,elem,lambda,xn,yn,e)
        call fieldu(elem,lambda,u)
!
        nu(:,:,:,:)  = 0._DPC
        nu(1,1,1,1)  = cmplx(list(1),0._DP,DPC)
        nu(2,2,1,1)  = cmplx(list(1),0._DP,DPC)
        nu(1,1,2,2)  = cmplx(list(1),0._DP,DPC)
        nu(2,2,2,2)  = cmplx(list(1),0._DP,DPC)
        nu(1,1,4,4)  = cmplx(eps0*list(3),0._DP,DPC)
        nu(2,2,4,4)  = cmplx(eps0*list(4),0._DP,DPC)
        nu(1,1,5,4)  = cmplx(abs(u(5))*(list(5)*list(7)*Elch)/(Kboltz*list(8)),0._DP,DPC)
        nu(2,2,5,4)  = cmplx(abs(u(5))*(list(6)*list(7)*Elch)/(Kboltz*list(8)),0._DP,DPC)
        nu(1,1,5,5)  = cmplx(list(5),0._DP,DPC)
        nu(2,2,5,5)  = cmplx(list(6),0._DP,DPC)
        nu(1,1,6,6)  = cmplx(list(8)/(list(2)*list(10)),0._DP,DPC)
        nu(2,2,6,6)  = cmplx(list(9)/(list(2)*list(10)),0._DP,DPC)
!
        beta(:,:,:)  = 0._DPC
        beta(1,1,1)  = cmplx(list(2)*abs(u(1)),0._DP,DPC)
        beta(2,1,1)  = cmplx(list(2)*abs(u(2)),0._DP,DPC)
        beta(1,2,2)  = cmplx(list(2)*abs(u(1)),0._DP,DPC)
        beta(2,2,2)  = cmplx(list(2)*abs(u(2)),0._DP,DPC)
        beta(1,1,3)  = cmplx(1._DP,0._DP,DPC)
        beta(2,2,3)  = cmplx(1._DP,0._DP,DPC)
        beta(1,3,1)  = cmplx(1._DP,0._DP,DPC)
        beta(2,3,2)  = cmplx(1._DP,0._DP,DPC)
        beta(1,6,6)  = cmplx(abs(u(1)),0._DP,DPC)
        beta(2,6,6)  = cmplx(abs(u(2)),0._DP,DPC)
!
        gamma(:,:,:) = 0._DPC
        gamma(1,5,5) = cmplx(-abs(u(1)),0._DP,DPC)
        gamma(2,5,5) = cmplx(-abs(u(2)),0._DP,DPC)
!
        alpha(:,:)   = 0._DPC
        alpha(1,1)   = cmplx(0._DP,list(2)*omega,DPC)
        alpha(2,2)   = cmplx(0._DP,list(2)*omega,DPC)
        f1(:)        = 0._DPC
        f1(6)        = cmplx(list(11),0._DP,DPC)
        f2(:,:)      = 0._DPC
        deallocate(list)
!
!
      case ('HEATTR')
        var11  = list(1)
        var22  = list(2)
        angle  = rad*list(3)
!  calculate tensor for diffusion term and assign result to nu
        invert = .false.
        call calctensor (var11, var22, angle, invert, nu(:,:,1,1))
        zcomp = (list(4)+list(5))/list(6)
!  calculate complex value for power density
        phase = rad*list(10)
!  modify for electro-thermal heat source
        scalar = list(9)
        if (multiphysics) then
          call sourcemodify(xs,ys,scalar)
        end if
        source1 = scalar * exp(cmplx(0._DP,phase,DPC))
        source2 = list(11)
        alpha(1,1) = zcomp+cmplx(0._DP,omega*list(7)*list(8),DPC)
        beta(1,1,1) = list(7) * list(8) * list(12)
        beta(2,1,1) = list(7) * list(8) * list(13)
        gamma(:,1,1) = 0._DPC
        f1(1) = source1+zcomp*source2
        f2(:,1) = 0._DPC
        deallocate(list)
!
!
      case ('TRANHEATTR')
!  Transient Heat Transfer physics mode
!  nnat = 1
!
!  nu is of the form nu(2,2,nnat,nnat)
        nu(1,1,1,1)  = cmplx(list(1),0._DP,DPC)
        nu(1,2,1,1)  = 0._DPC
        nu(2,1,1,1)  = 0._DPC
        nu(2,2,1,1)  = cmplx(list(2),0._DP,DPC)
!  gamma is of the form gamma(2,nnat,nnat)
        gamma(:,1,1) = 0._DPC
!  beta is of the form beta(2,nnat,nnat)
        beta(:,1,1) = 0._DPC
!  alpha is of the form alpha(nnat,nnat)
        alpha(1,1) = cmplx(list(4)*list(5),0._DP,DPC)
!  em is of the form em(nnat,nnat)
        em(1,1) = 0._DPC
!  f1 is of the form f1(nnat)
        f1(1) = list(6)
!  f2 is of the form f2(2,nnat)
        f2(:,1) = 0._DPC
        deallocate(list)
!
!
      case ('LAPLACE')
        nu(1,1,1,1)  = cmplx(1._DP,0._DP,DPC)
        nu(1,2,1,1)  = 0._DPC
        nu(2,1,1,1)  = nu(1,2,1,1)
        nu(2,2,1,1)  = cmplx(1._DP,0._DP,DPC)
        alpha(1,1)   = 0._DPC
        beta(:,1,1)  = 0._DPC
        gamma(:,1,1) = 0._DPC
        f1(1)        = 0._DPC
        f2(:,1)      = 0._DPC
!
!
      case ('MAGNETICVP')
!  compute inverse of permeability tensor
        v11 = mu0 * list(1)
        v22 = mu0 * list(2)
        angle = rad * list(3)
        invert = .true.
        call calctensor(v11, v22, angle, invert, gam)
        nu(1,1,1,1)  = gam(2,2)
        nu(1,2,1,1)  =-gam(2,1)
        nu(2,1,1,1)  =-gam(1,2)
        nu(2,2,1,1)  = gam(1,1)
!  beta = kappa * velocity
        beta(1,1,1)  = list(9) * list(5)
        beta(2,1,1)  = list(9) * list(6)
        gamma(:,1,1) = 0._DPC
!  alpha = j * omega * kappa
        alpha(1,1)   = cmplx(0._DP, omega*list(9), DPC)
!  f1 = J_z
        phase = rad*list(8)
        f1(1) = list(7) * exp(cmplx(0._DP,phase,DPC))
!  f2  contribution of permanent magnets
        jmx= list(4) * cos(angle)
        jmy= list(4) * sin(angle)
        f2(1,1) =  gam(2,1) * jmx + gam(2,2) * jmy
        f2(2,1) = -gam(1,1) * jmx - gam(1,2) * jmy
        deallocate(list)
!
!
      case ('STOKES')
!       list(1) = nu
!       list(2) = rho
!       If non-dimensional is used, nu = 1/Re and rho = 1
!     
        nu(1,1,1,1)  = cmplx(list(1)/list(2),0._DP,DPC)  
        nu(2,2,1,1)  = nu(1,1,1,1)
        nu(1,2,1,1)  = 0._DPC
        nu(2,1,1,1)  = nu(1,2,1,1)
        beta(:,1,1)  = 0._DPC
        gamma(:,1,1) = 0._DPC
        alpha(1,1)   = cmplx(0._DP,omega,DPC)
        f1(1)        = cmplx(list(3)/list(2),0._DP,DPC)
        f2(1,1)      = cmplx(list(4)/list(2),0._DP,DPC)
        f2(2,1)      = cmplx(list(5)/list(2),0._DP,DPC)
        deallocate(list)
!
!
      case ('STATCURRENT')
!  var11 and var22 contain electric conductivity kappa in crystal axes
        var11 = list(1)
        var22 = list(2)
        angle = rad*list(3)
!  calculate material tensor for diffusion term and assign result to nu
        invert = .false.
        call calctensor (var11, var22, angle, invert, nu(:,:,1,1))
        nu(:,:,1,1)  = nu(:,:,1,1)*list(4)
!  assign scalar source term d*Q
        f1(1)        = cmplx(list(4)*list(5),0._DP,DPC)
!  assign vectorial source term (external current density)
        f2(1,1)      = cmplx(-list(4)*list(6),0._DP,DPC)
        f2(2,1)      = cmplx(-list(4)*list(7),0._DP,DPC)
!  assign alpha, beta and f2
        alpha(1,1)   = 0._DPC
        beta(:,1,1)  = 0._DPC
        gamma(:,1,1) = 0._DPC
        deallocate(list)       
!
!
      case ('TEWAVE')
!  compute inverted permeability tensor
        call invertmat3(mur,matten)
!  assign nu components due to differences between rot|rot and div|grad pde type
        nu(1,1,1,1)  =  matten(2,2)
        nu(1,2,1,1)  = -matten(2,1)
        nu(2,1,1,1)  = -matten(1,2)
        nu(2,2,1,1)  =  matten(1,1)
!  assign alpha and beta
        alpha(1,1)   = -(omega**2)*mu0*eps0*epsr(3,3)
        beta(:,1,1)  = 0._DPC
        gamma(:,1,1) = 0._DPC
!  compute complex values for current density f1 and assign f2
        f1(1)        = -source(3) * exp(cmplx(0._DP,phase,DPC)) * cmplx(0._DP,omega*mu0,DPC)
        f2(:,1)      = 0._DPC
        deallocate(list)
!
!
      case ('TMWAVE')
!  compute inverted permeability tensor
        call invertmat3(epsr,matten)
!  assign nu components due to differences between rot|rot and div|grad pde type
        nu(1,1,1,1)  =  matten(2,2)
        nu(1,2,1,1)  = -matten(2,1)
        nu(2,1,1,1)  = -matten(1,2)
        nu(2,2,1,1)  =  matten(1,1)
!  assign alpha, beta and f1
        alpha(1,1)   = -(omega**2)*eps0*mu0*mur(3,3)
        beta(:,1,1)  = 0._DPC
        gamma(:,1,1) = 0._DPC
        f1(1) = 0._DPC
!  calculate complex values for current density f2
        f2(1,1)      =  source(2) * exp(cmplx(0._DP,phase,DPC))
        f2(2,1)      = -source(1) * exp(cmplx(0._DP,phase,DPC))
        f2(:,1)      = matmul(nu(:,:,1,1),f2(:,1))
        deallocate(list)
!
!
      case ('THERMOELECTRIC')
!  nnat=2 in the thermoelectric problem
!  the thermal is considered as the first nature
!  the electric is considered as the second nature
!
!        var11  = list(1)
!        var22  = list(2)
!        angle  = rad*list(3)
!
!  nu is of the form nu(2,2,nnat,nnat)
!  thermal nu = lamda/Tref
        nu(1,1,1,1) = cmplx(list(1)/list(23),0._DP,DPC)
        nu(1,2,1,1) = 0._DPC
        nu(2,1,1,1) = 0._DPC
        nu(2,2,1,1) = cmplx(list(2)/list(23),0._DP,DPC)
!  thermal-electric nu = kappa * Seebeck coef
        nu(1,1,1,2) = cmplx(list(14)*list(22),0._DP,DPC)
        nu(1,2,1,2) = 0._DPC
        nu(2,1,1,2) = 0._DPC
        nu(2,2,1,2) = cmplx(list(15)*list(22),0._DP,DPC)
!  electric-thermal nu = kappa * Seebeck coef
        nu(1,1,2,1) = cmplx(list(14)*list(22),0._DP,DPC)
        nu(1,2,2,1) = 0._DPC
        nu(2,1,2,1) = 0._DPC
        nu(2,2,2,1) = cmplx(list(15)*list(22),0._DP,DPC)
!  electric nu
        nu(1,1,2,2) = cmplx(list(14),0._DP,DPC)
        nu(1,2,2,2) = 0._DPC
        nu(2,1,2,2) = 0._DPC
        nu(2,2,2,2) = cmplx(list(15),0._DP,DPC)

!  calculate tensor for diffusion term and assign result to nu
!        invert = .false.
!        call calctensor (var11, var22, angle, invert, nu(:,:,1,1))
!        call calctensor (var11, var22, angle, invert, nu(:,:,1,2))
!        call calctensor (var11, var22, angle, invert, nu(:,:,2,1))
!        call calctensor (var11, var22, angle, invert, nu(:,:,2,2))

!   (hu+ho)/d      z-component of convective heat transfer
        zcomp = (list(4)+list(5))/list(6)
!  calculate complex value for power density
        phase = rad*list(10)
!  modify for electro-thermal heat source
        scalar = list(9)
        if (multiphysics) then
          call sourcemodify(xs,ys,scalar)
        end if
        source1 = scalar * exp(cmplx(0._DP,phase,DPC))

!  electroheating, only if nonlinear
        call getsetting('NONLINEAR',nonlinear)
        if (nonlinear .eq. 'YES') then
          call fieldu(elem,lambda,u,gradu)
!  calculate complex value for power density: J*conjg(E)
          source1 = source1 + dot_product(matmul(nu(:,:,2,2),gradu(1:2,2)) , gradu(1:2,2))
!  kappa = kappa*(1-alpha*(T-Tref))
          nu(1,1,2,2) = nu(1,1,2,2) / (1._DP + list(24)*(u(1)-list(23)) )
          nu(2,2,2,2) = nu(2,2,2,2) / (1._DP + list(24)*(u(1)-list(23)) )
        end if

        source2 = list(11)
!  alpha is of the form alpha(nnat,nnat)
        alpha(1,1) = zcomp
        alpha(1,2) = 0._DPC
        alpha(2,1) = 0._DPC
        alpha(2,2) = 0._DPC
!  beta is of the form beta(2,nnat,nnat)
        beta(1,1,1) = 0._DPC
        beta(2,1,1) = 0._DPC
        beta(1,1,2) = 0._DPC
        beta(2,1,2) = 0._DPC
        beta(1,2,1) = 0._DPC
        beta(2,2,1) = 0._DPC
        beta(1,2,2) = 0._DPC
        beta(2,2,2) = 0._DPC
!  gamma is of the form gamma(2,nnat,nnat)
        gamma(:,:,:) = 0._DPC
!  f1 is of the form f1(nnat)
        f1(1) = source1 + zcomp*source2
        f1(2) =  cmplx(list(17)*list(18),0._DP,DPC)
!  f2 is of the form f2(2,nnat)
        f2(1,1) = 0._DPC
        f2(2,1) = 0._DPC
        f2(1,2) = cmplx(-list(17)*list(19),0._DP,DPC)
        f2(2,2) = cmplx(-list(17)*list(20),0._DP,DPC)
        deallocate(list)
!
      case ('PLANE STRAIN')
!  nnat=2 in the 2D Plane Strain problem
!  in the first nature we consider the displacement solution in x
!  in the second nature we consider the displacement solution in y

!  First calculate Lame coefficients 1 and 2 (lamda and mu) from
!  given Young's modulus (denoted here with E)
!  and Poisson's ratio (denoted here with v)
!  The formulae are: Lame1 = E*v/((1+v)*(1-2v))
!                    Lame2 = E/(2*(1+v))
        Lame1 = ( list(2)*list(3) ) / ( (1._DP+list(3))*(1._DP-2._DP*list(3)) )
        Lame2 = ( list(2) ) / ( 2._DP+(2._DP*list(3)) )

!  nu is of the form nu(2,2,nnat,nnat)
!  submatrix 1-1, x <- x
        nu(1,1,1,1) = cmplx(Lame1+2._DP*Lame2,0._DP,DPC)
        nu(1,2,1,1) = 0._DPC
        nu(2,1,1,1) = 0._DPC
        nu(2,2,1,1) = cmplx(Lame2,0._DP,DPC)
!  submatrix 1-2, x <- y
        nu(1,1,1,2) = 0._DPC
        nu(1,2,1,2) = cmplx(Lame1,0._DP,DPC)
        nu(2,1,1,2) = cmplx(Lame2,0._DP,DPC)
        nu(2,2,1,2) = 0._DPC
!  submatrix 2-1, y <- x
        nu(1,1,2,1) = 0._DPC
        nu(1,2,2,1) = cmplx(Lame2,0._DP,DPC)
        nu(2,1,2,1) = cmplx(Lame1,0._DP,DPC)
        nu(2,2,2,1) = 0._DPC
!  submatrix 2-2, y <- y
        nu(1,1,2,2) = cmplx(Lame2,0._DP,DPC)
        nu(1,2,2,2) = 0._DPC
        nu(2,1,2,2) = 0._DPC
        nu(2,2,2,2) = cmplx(Lame1+2._DP*Lame2,0._DP,DPC)

! for the case of an anisotropic material
        if ( list(14) .eq. 1._DP ) then
!  submatrix 1-1, x <- x
          nu(1,1,1,1) = list(5)
          nu(1,2,1,1) = list(7)
          nu(2,1,1,1) = list(11)
          nu(2,2,1,1) = list(13)
!  submatrix 1-2, x <- y
          nu(1,1,1,2) = list(7)
          nu(1,2,1,2) = list(6)
          nu(2,1,1,2) = list(13)
          nu(2,2,1,2) = list(12)
!  submatrix 2-1, y <- x
          nu(1,1,2,1) = list(11)
          nu(1,2,2,1) = list(13)
          nu(2,1,2,1) = list(8)
          nu(2,2,2,1) = list(10)
!  submatrix 2-2, y <- y
          nu(1,1,2,2) = list(13)
          nu(1,2,2,2) = list(12)
          nu(2,1,2,2) = list(10)
          nu(2,2,2,2) = list(9)
        end if
!
!  gamma is of the form gamma(2,nnat,nnat)
        gamma(:,:,:) = 0._DPC
!
!  alpha is of the form alpha(nnat,nnat)
        alpha(1,1) = cmplx(-((omega**2)*list(1)),0._DP,DPC)
        alpha(1,2) = 0._DPC
        alpha(2,1) = 0._DPC
        alpha(2,2) = cmplx(-((omega**2)*list(1)),0._DP,DPC)
!  beta is of the form beta(2,nnat,nnat)
        beta(:,:,:) = 0._DPC
!  f1 is of the form f1(nnat)
        f1(1) = cmplx(list(15),0._DP,DPC)
        f1(2) = cmplx(list(16),0._DP,DPC)
!  f2 is of the form f2(2,nnat)
        f2(1,1) = 0._DPC
        f2(2,1) = 0._DPC
        f2(1,2) = 0._DPC
        f2(2,2) = 0._DPC
        deallocate(list)
!
      case ('TRANPLSTRAIN')
!  Transient Plane Strain physics mode
!  nnat = 2
        Lame1 = ( list(2)*list(3) ) / ( (1._DP+list(3))*(1._DP-2._DP*list(3)) )
        Lame2 = ( list(2) ) / ( 2._DP+(2._DP*list(3)) )
!  nu is of the form nu(2,2,nnat,nnat)
!  submatrix 1-1, x <- x
        nu(1,1,1,1) = cmplx(Lame1+2._DP*Lame2,0._DP,DPC)
        nu(1,2,1,1) = 0._DPC
        nu(2,1,1,1) = 0._DPC
        nu(2,2,1,1) = cmplx(Lame2,0._DP,DPC)
!  submatrix 1-2, x <- y
        nu(1,1,1,2) = 0._DPC
        nu(1,2,1,2) = cmplx(Lame1,0._DP,DPC)
        nu(2,1,1,2) = cmplx(Lame2,0._DP,DPC)
        nu(2,2,1,2) = 0._DPC
!  submatrix 2-1, y <- x
        nu(1,1,2,1) = 0._DPC
        nu(1,2,2,1) = cmplx(Lame2,0._DP,DPC)
        nu(2,1,2,1) = cmplx(Lame1,0._DP,DPC)
        nu(2,2,2,1) = 0._DPC
!  submatrix 2-2, y <- y
        nu(1,1,2,2) = cmplx(Lame2,0._DP,DPC)
        nu(1,2,2,2) = 0._DPC
        nu(2,1,2,2) = 0._DPC
        nu(2,2,2,2) = cmplx(Lame1+2._DP*Lame2,0._DP,DPC)
!  gamma is of the form gamma(2,nnat,nnat)
        gamma(:,:,:) = 0._DPC
!  beta is of the form beta(2,nnat,nnat)
        beta(:,:,:) = 0._DPC
!  alpha is of the form alpha(nnat,nnat)
        alpha(:,:) = 0._DPC
!  f1 is of the form f1(nnat)
        f1(1) = list(5)
        f1(2) = list(6)
!  f2 is of the form f2(2,nnat)
        f2(:,:) = 0._DPC
        deallocate(list)
!
      case ('PLANE STRESS')
!  nnat=2 in the 2D Plane Stress problem
!  in the first nature we consider the attempt to solve in x
!  in the second nature we consider the attempt to solve in y

!  First calculate Lame coefficients 1 and 2 (lamda and mu) from
!  given Young's modulus (denoted here with E)
!  and Poisson's ratio (denoted here with v)
!  The formulae are: Lame1 = E*v/((1+v)*(1-2v))
!                    Lame2 = E/(2*(1+v))

        Lame1 = ( list(2)*list(3) ) / ( (1._DP+list(3))*(1._DP-2._DP*list(3)) )
        Lame2 = ( list(2) ) / ( 2._DP+(2._DP*list(3)) )
        
!  nu is of the form nu(2,2,nnat,nnat)

!  submatrix 1-1, x <- x
        nu(1,1,1,1) = cmplx(list(2)/(1._DP-(list(3)*list(3))),0._DP,DPC)
        nu(1,2,1,1) = 0._DPC
        nu(2,1,1,1) = 0._DPC
        nu(2,2,1,1) = cmplx(Lame2,0._DP,DPC)

!  submatrix 1-2, x <- y
        nu(1,1,1,2) = 0._DPC
        nu(1,2,1,2) = cmplx((list(2)*list(3))/(1._DP-(list(3)*list(3))),0._DP,DPC)
        nu(2,1,1,2) = cmplx(Lame2,0._DP,DPC)
        nu(2,2,1,2) = 0._DPC
        
!  submatrix 2-1, y <- x
        nu(1,1,2,1) = 0._DPC
        nu(1,2,2,1) = cmplx(Lame2,0._DP,DPC)
        nu(2,1,2,1) = cmplx((list(2)*list(3))/(1._DP-(list(3)*list(3))),0._DP,DPC)
        nu(2,2,2,1) = 0._DPC
        
!  submatrix 2-2, y <- y
        nu(1,1,2,2) = cmplx(Lame2,0._DP,DPC)
        nu(1,2,2,2) = 0._DPC
        nu(2,1,2,2) = 0._DPC
        nu(2,2,2,2) = cmplx(list(2)/(1._DP-(list(3)*list(3))),0._DP,DPC)
        
! for the case of an anisotropic material
        if ( list(14) .eq. 1._DP ) then
!  submatrix 1-1, x <- x
          nu(1,1,1,1) = cmplx(list(5),0._DP,DPC)
          nu(1,2,1,1) = cmplx(list(7),0._DP,DPC)
          nu(2,1,1,1) = cmplx(list(11),0._DP,DPC)
          nu(2,2,1,1) = cmplx(list(13),0._DP,DPC)
!  submatrix 1-2, x <- y
          nu(1,1,1,2) = cmplx(list(7),0._DP,DPC)
          nu(1,2,1,2) = cmplx(list(6),0._DP,DPC)
          nu(2,1,1,2) = cmplx(list(13),0._DP,DPC)
          nu(2,2,1,2) = cmplx(list(12),0._DP,DPC)
!  submatrix 2-1, y <- x
          nu(1,1,2,1) = cmplx(list(11),0._DP,DPC)
          nu(1,2,2,1) = cmplx(list(13),0._DP,DPC)
          nu(2,1,2,1) = cmplx(list(8),0._DP,DPC)
          nu(2,2,2,1) = cmplx(list(10),0._DP,DPC)
!  submatrix 2-2, y <- y
          nu(1,1,2,2) = cmplx(list(13),0._DP,DPC)
          nu(1,2,2,2) = cmplx(list(12),0._DP,DPC)
          nu(2,1,2,2) = cmplx(list(10),0._DP,DPC)
          nu(2,2,2,2) = cmplx(list(9),0._DP,DPC)
        end if
        
!  alpha is of the form alpha(nnat,nnat)
        alpha(1,1)  = 0._DPC
        alpha(1,2)  = 0._DPC
        alpha(2,1)  = 0._DPC
        alpha(2,2)  = 0._DPC
!  submatrix 1-1, x <- x
        alpha(1,1) = cmplx(-((omega**2)*list(1)),0._DP,DPC)
!  submatrix 2-2, y <- y
        alpha(2,2) = cmplx(-((omega**2)*list(1)),0._DP,DPC)

!  beta is of the form beta(2,nnat,nnat)
        beta(:,:,:) = 0._DPC
!  gamma
        gamma(:,:,:) = 0._DPC
!  f1 is of the form f1(nnat)
        f1(1)       = list(15)
        f1(2)       = list(16)
!  f2 is of the form f2(2,nnat)
        f2(1,1)     = 0._DPC
        f2(2,1)     = 0._DPC
        f2(1,2)     = 0._DPC
        f2(2,2)     = 0._DPC
        deallocate(list)
!
      case ('THERMOELASTIC')
!  nnat=3 in the 2D Thermoelastic problem
!  in the first nature we consider the displacement in x
!  in the second nature we consider the displacement in y
!  in the third nature we consider the temperature
!
!  First calculate Lame coefficients 1 and 2 (lamda and mu) from
!  given Young's modulus (denoted here with E)
!  and Poisson's ratio (denoted here with v)
!  The formulae are: Lame1 = E*v/((1+v)*(1-2v))
!                    Lame2 = E/(2*(1+v))
        Lame1 = ( list(2)*list(3) ) / ( (1._DP+list(3))*(1._DP-2._DP*list(3)) )
        Lame2 = ( list(2) ) / ( 2._DP+(2._DP*list(3)) )
!
!  Calculate the Stress-Temperature tensor
!  M: Stress-Temperature coefficients, A: Thermal Expansion coefficients
!  C: Rank-4 Stiffness tensor from Mechanics
!
!  The indices are as:
!  M(i,j) = C(i,j,k,l)*A(k,l)
!
!                |C1111 C1112 | C1211 C1212|
!  |M11 | M12|   |C1121 C1122 | C1221 C1222| |A11 | A12|
!  |----|----| = |------------|------------|*|----|----|
!  |M21 | M22|   |C2111 C2112 | C2211 C2212| |A21 | A22|
!                |C2121 C2122 | C2221 C2222|
!
!  M11 = C1111*A11 + C1112*A12 + C1121*A21 + C1122*A22
!  M12 = C1211*A11 + C1212*A12 + C1221*A21 + C1222*A22
!  M21 = C2111*A11 + C2112*A12 + C2121*A21 + C2122*A22
!  M22 = C2211*A11 + C2212*A12 + C2221*A21 + C2222*A22
!
!  Due to the fact that the Rank-4 Stiffness tensor is symmetric
!  we can write it as follows:
!  We replace indices 11->1, 22->2, 12->3, 21->3 in pairs of 2
!  Then we write the Rank-4 tensor as a Rank-2 tensor
!
!  |C1111 C1112 | C1211 C1212|   |C11 C13 | C31 C33|
!  |C1121 C1122 | C1221 C1222|   |C13 C12 | C33 C32|   |C11 C12 C13|
!  |------------|------------| = |--------|--------| = |C21 C22 C23|
!  |C2111 C2112 | C2211 C2212|   |C31 C33 | C21 C23|   |C31 C32 C33|
!  |C2121 C2122 | C2221 C2222|   |C33 C32 | C23 C22|
!
!  We can now rewrite the calculation of the Stress-Temperature coefficients
!
!  M11 = C11*A11 + C13*A12 + C13*A21 + C12*A22
!  M12 = C31*A11 + C33*A12 + C33*A21 + C32*A22
!  M21 = C31*A11 + C33*A12 + C33*A21 + C32*A22
!  M22 = C21*A11 + C23*A12 + C23*A21 + C22*A22
!
        ATensNoSym = reshape((/ list(15), list(16), list(17), list(18) /) &
                          &  ,  (/2,2/), (/0._DP/), (/2,1/))
        if (list(41) .eq. 1._DP) then
          CTensNoSym = reshape((/ list(5) , list(7) , list(7) , list(6)     &
                            &  ,  list(11), list(13), list(13), list(12)    &
                            &  ,  list(11), list(13), list(13), list(12)    &
                            &  ,  list(8) , list(10), list(10), list(9)  /) &
                            &  ,  (/2,2,2,2/), (/0._DP/), (/4,3,2,1/))
        else
          CTensNoSym = reshape((/ list(25), list(26), list(27), list(28)    &
                            &  ,  list(29), list(30), list(31), list(32)    &
                            &  ,  list(33), list(34), list(35), list(36)    &
                            &  ,  list(37), list(38), list(39), list(40) /) &
                            &  ,  (/2,2,2,2/), (/0._DP/), (/4,3,2,1/))
        end if
        
        call multtensor(CTensNoSym, ATensNoSym, MTensNoSym)
!
!  nu is of the form nu(2,2,nnat,nnat)
!  submatrix 1-1, x <- x
        nu(1,1,1,1) = cmplx(list(2)/(1._DP-(list(3)*list(3))),0._DP,DPC)
        nu(1,2,1,1) = 0._DPC
        nu(2,1,1,1) = 0._DPC
        nu(2,2,1,1) = cmplx(Lame2,0._DP,DPC)
!  submatrix 1-2, x <- y
        nu(1,1,1,2) = 0._DPC
        nu(1,2,1,2) = cmplx((list(2)*list(3))/(1._DP-(list(3)*list(3))),0._DP,DPC)
        nu(2,1,1,2) = cmplx(Lame2,0._DP,DPC)
        nu(2,2,1,2) = 0._DPC
!  submatrix 1-3, x <- T
        nu(:,:,1,3) = 0._DPC
!  submatrix 2-1, y <- x
        nu(1,1,2,1) = 0._DPC
        nu(1,2,2,1) = cmplx(Lame2,0._DP,DPC)
        nu(2,1,2,1) = cmplx((list(2)*list(3))/(1._DP-(list(3)*list(3))),0._DP,DPC)
        nu(2,2,2,1) = 0._DPC
!  submatrix 2-2, y <- y
        nu(1,1,2,2) = cmplx(Lame2,0._DP,DPC)
        nu(1,2,2,2) = 0._DPC
        nu(2,1,2,2) = 0._DPC
        nu(2,2,2,2) = cmplx(list(2)/(1._DP-(list(3)*list(3))),0._DP,DPC)
!  submatrix 2-3, y <- T
        nu(:,:,2,3) = 0._DPC
!  submatrix 3-1, T <- x
        nu(:,:,3,1) = 0._DPC
!  submatrix 3-2, T <- y
        nu(:,:,3,2) = 0._DPC
!  submatrix 3-3, T <- T
        nu(1,1,3,3) = cmplx(list(19),0._DP,DPC)
        nu(1,2,3,3) = 0._DPC
        nu(2,1,3,3) = 0._DPC
        nu(2,2,3,3) = cmplx(list(20),0._DP,DPC)

!  gamma is of the form gamma(2,nnat,nnat)
!  submatrix 1-1, x <- x
        gamma(:,1,1) = 0._DPC
!  submatrix 1-2, x <- y
        gamma(:,1,2) = 0._DPC
!  submatrix 1-3, x <- T
        gamma(1,1,3) = (-MTensNoSym(1,1))
        gamma(2,1,3) = (-MTensNoSym(1,2))
!  submatrix 2-1, y <- x
        gamma(:,2,1) = 0._DPC
!  submatrix 2-2, y <- y
        gamma(:,2,2) = 0._DPC
!  submatrix 2-3, y <- T
        gamma(1,2,3) = (-MTensNoSym(2,1))
        gamma(2,2,3) = (-MTensNoSym(2,2))
!  submatrix 3-1, T <- x
        gamma(:,3,1) = 0._DPC
!  submatrix 3-2, T <- y
        gamma(:,3,2) = 0._DPC
!  submatrix 3-3, T <- T
        gamma(:,3,3) = 0._DPC

!  beta is of the form beta(2,nnat,nnat)
!  submatrix 1-1, x <- x
        beta(:,1,1) = 0._DPC
!  submatrix 1-2, x <- y
        beta(:,1,2) = 0._DPC
!  submatrix 1-3, x <- T
        beta(:,1,3) = 0._DPC
!  submatrix 2-1, y <- x
        beta(:,2,1) = 0._DPC
!  submatrix 2-2, y <- y
        beta(:,2,2) = 0._DPC
!  submatrix 2-3, y <- T
        beta(:,2,3) = 0._DPC
!  submatrix 3-1, T <- x
        beta(1,3,1) = cmplx(0._DP,1._DP,DPC)*(-omega)*list(23)*MTensNoSym(1,1)
        beta(2,3,1) = cmplx(0._DP,1._DP,DPC)*(-omega)*list(23)*MTensNoSym(1,2)
!  submatrix 3-2, T <- y
        beta(1,3,2) = cmplx(0._DP,1._DP,DPC)*(-omega)*list(23)*MTensNoSym(2,1)
        beta(2,3,2) = cmplx(0._DP,1._DP,DPC)*(-omega)*list(23)*MTensNoSym(2,2)
!  submatrix 3-3, T <- T
        beta(:,3,3) = 0._DPC
!
!  alpha is of the form alpha(nnat,nnat)
!  submatrix 1-1, x <- x
        alpha(1,1) = cmplx(-((omega**2)*list(1)),0._DP,DPC) !0._DPC
        alpha(1,2) = 0._DPC
        alpha(1,3) = 0._DPC
        alpha(2,1) = 0._DPC
!  submatrix 2-2, y <- y
        alpha(2,2) = cmplx(-((omega**2)*list(1)),0._DP,DPC) !0._DPC
        alpha(2,3) = 0._DPC
        alpha(3,1) = 0._DPC
        alpha(3,2) = 0._DPC
!  submatrix 3-3, T <- T
        alpha(3,3) = cmplx(0._DP,(omega*list(1)*list(22)),DPC) !0._DPC
!
!  f1 is of the form f1(nnat)
!  x
        f1(1) = 0._DPC
!  y
        f1(2) = 0._DPC
!  T
        f1(3) = 0._DPC

!  f2 is of the form f2(2,nnat)
!  x
        f2(1,1) = 0._DPC
        f2(2,1) = 0._DPC
!  y
        f2(1,2) = 0._DPC
        f2(2,2) = 0._DPC
!  T
        f2(1,3) = 0._DPC
        f2(2,3) = 0._DPC
        
        deallocate(list)
!
      case ('GRAVITOELASTIC')
!  nnat=3 in the 2D Gravitoelastic problem
!  in the first nature we consider the displacement in x
!  in the second nature we consider the displacement in y
!  in the third nature we consider the gravitational field
!
!  First calculate Lame coefficients 1 and 2 (lamda and mu) from
!  given Young's modulus (denoted here with E)
!  and Poisson's ratio (denoted here with v)
!  The formulae are: Lame1 = E*v/((1+v)*(1-2v))
!                    Lame2 = E/(2*(1+v))
        Lame1 = ( list(2)*list(3) ) / ( (1._DP+list(3))*(1._DP-2._DP*list(3)) )
        Lame2 = ( list(2) ) / ( 2._DP+(2._DP*list(3)) )
!
!  Calculate the Stress-Gravity tensor
!        CTensNoSym = reshape((/ Lame1+2._DP*Lame2 , 0._DP , 0._DP , Lame1     &
!                          &  ,  0._DP, Lame2, Lame2, 0._DP    &
!                          &  ,  0._DP, Lame2, Lame2, 0._DP    &
!                          &  ,  Lame1 , 0._DP, 0._DP, Lame1+2._DP*Lame2  /) &
!                          &  ,  (/2,2,2,2/), (/0._DP/), (/4,3,2,1/))
!        GTensNoSym = reshape((/  list(17)       ,0._DP             &
!                          &  ,   0.5_DP*list(17),0._DP             &
!                          &  ,   0._DP          ,0.5_DP*list(17)   &
!                          &  ,   0._DP          ,list(17)   &  /)  &
!                          &  ,    (/2,2,2/), (/0._DP/), (/3,2,1/))
!        call multtensor(CTensNoSym, GTensNoSym, LTensNoSym)
!
!  nu is of the form nu(2,2,nnat,nnat)
!  submatrix 1-1, x <- x
        nu(1,1,1,1) = cmplx(Lame1+2._DP*Lame2,0._DP,DPC)
        nu(1,2,1,1) = 0._DPC
        nu(2,1,1,1) = 0._DPC
        nu(2,2,1,1) = cmplx(Lame2,0._DP,DPC)
!  submatrix 1-2, x <- y
        nu(1,1,1,2) = 0._DPC
        nu(1,2,1,2) = cmplx(Lame1,0._DP,DPC)
        nu(2,1,1,2) = cmplx(Lame2,0._DP,DPC)
        nu(2,2,1,2) = 0._DPC
!  submatrix 1-3, x <- gravity
        nu(:,:,1,3) = 0._DPC
!  submatrix 2-1, y <- x
        nu(1,1,2,1) = 0._DPC
        nu(1,2,2,1) = cmplx(Lame2,0._DP,DPC)
        nu(2,1,2,1) = cmplx(Lame1,0._DP,DPC)
        nu(2,2,2,1) = 0._DPC
!  submatrix 2-2, y <- y
        nu(1,1,2,2) = cmplx(Lame2,0._DP,DPC)
        nu(1,2,2,2) = 0._DPC
        nu(2,1,2,2) = 0._DPC
        nu(2,2,2,2) = cmplx(Lame1+2._DP*Lame2,0._DP,DPC)
!  submatrix 2-3, y <- gravity
        nu(:,:,2,3) = 0._DPC
!  submatrix 3-1, gravity <- x
        nu(:,:,3,1) = 0._DPC
!  submatrix 3-2, gravity <- y
        nu(:,:,3,2) = 0._DPC
!  submatrix 3-3, gravity <- gravity
        nu(1,1,3,3) = cmplx(1._DP,0._DP,DPC)
        nu(1,2,3,3) = 0._DPC
        nu(2,1,3,3) = 0._DPC
        nu(2,2,3,3) = cmplx(1._DP,0._DP,DPC)

!  gamma is of the form gamma(2,nnat,nnat)
        gamma(:,:,:) = 0._DPC
!  submatrix 1-1, x <- x
!  submatrix 1-2, x <- y
!  submatrix 1-3, x <- gravity
!  submatrix 2-1, y <- x
!  submatrix 2-2, y <- y
!  submatrix 2-3, y <- gravity
!  submatrix 3-1, gravity <- x
!  submatrix 3-2, gravity <- y
!  submatrix 3-3, gravity <- gravity

!  beta is of the form beta(2,nnat,nnat)
        beta(:,:,:) = 0._DPC
!  submatrix 1-1, x <- x
!  submatrix 1-2, x <- y
!  submatrix 1-3, x <- gravity
        beta(1,1,3) = cmplx(list(1),0._DP,DPC)
!  submatrix 2-1, y <- x
!  submatrix 2-2, y <- y
!  submatrix 2-3, y <- gravity
        beta(2,2,3) = cmplx(list(1),0._DP,DPC)
!  submatrix 3-1, gravity <- x
        beta(1,3,1) = cmplx(-(4._DP*pi*gravc*list(1)),0._DP,DPC)
!  submatrix 3-2, gravity <- y
        beta(2,3,2) = cmplx(-(4._DP*pi*gravc*list(1)),0._DP,DPC)
!  submatrix 3-3, gravity <- gravity
!
!  alpha is of the form alpha(nnat,nnat)
        alpha(:,:) = 0._DPC
!
!  f1 is of the form f1(nnat)
!  x
        f1(1) = cmplx(list(15),0._DP,DPC)
!  y
        f1(2) = cmplx(list(16),0._DP,DPC)
!  gravity
        f1(3) = cmplx(list(17),0._DP,DPC)

!  f2 is of the form f2(2,nnat)
!  x
        f2(1,1) = 0._DPC
        f2(2,1) = 0._DPC
!  y
        f2(1,2) = 0._DPC
        f2(2,2) = 0._DPC
!  gravity
        f2(1,3) = 0._DPC
        f2(2,3) = 0._DPC
        
        deallocate(list)
!
      case ('PIEZOELECTRIC')
!  nnat=3 in the 2D Piezoelectric problem
!  in the first nature we consider the displacement in x
!  in the second nature we consider the displacement in y
!  in the third nature we consider the electric potential
!  (electrostatics)
!
        CTensNoSym3D = reshape((/list(20),list(25),list(24),list(25),list(21),list(23),list(24),list(23),list(22)   &
                          &  ,   list(25),list(40),list(39),list(40),list(30),list(37),list(39),list(37),list(34)   &
                          &  ,   list(24),list(39),list(38),list(39),list(29),list(36),list(38),list(36),list(33)   &
                          &  ,   list(25),list(40),list(39),list(40),list(30),list(37),list(39),list(37),list(34)   &
                          &  ,   list(21),list(30),list(29),list(30),list(26),list(28),list(29),list(28),list(27)   &
                          &  ,   list(23),list(37),list(36),list(37),list(28),list(35),list(36),list(35),list(32)   &
                          &  ,   list(24),list(39),list(38),list(39),list(29),list(36),list(38),list(36),list(33)   &
                          &  ,   list(23),list(37),list(36),list(37),list(28),list(35),list(36),list(35),list(32)   &
                          &  ,   list(22),list(34),list(33),list(34),list(27),list(32),list(33),list(32),list(31)/) &
                          &  ,    (/3,3,3,3/), (/0._DP/), (/4,3,2,1/))
        ETensNoSym3D = reshape((/list(41),list(42),list(43)   &
                          &  ,   list(56),list(57),list(58)   &
                          &  ,   list(53),list(54),list(55)   &
                          &  ,   list(56),list(57),list(58)   &
                          &  ,   list(44),list(45),list(46)   &
                          &  ,   list(50),list(51),list(52)   &
                          &  ,   list(53),list(54),list(55)   &
                          &  ,   list(50),list(51),list(52)   &
                          &  ,   list(47),list(48),list(49)/) &
                          &  ,    (/3,3,3/), (/0._DP/), (/3,2,1/))
        epsrTensNoSym3D = reshape((/eps0*list(17),  0._DP ,  0._DP    &
                          &  ,        0._DP ,eps0*list(18),  0._DP    &
                          &  ,        0._DP ,  0._DP ,eps0*list(19)/) &
                          &  ,        (/3,3/), (/0._DP/), (/2,1/))
        TransfTens = reshape((/list( 8),list( 9),list(10)   &
                          &  , list(11),list(12),list(13)   &
                          &  , list(14),list(15),list(16)/) &
                          &  ,  (/3,3/), (/0._DP/), (/2,1/))

        call transftensor(CTensNoSym3D, TransfTens)
        call transftensor(ETensNoSym3D, TransfTens)
        call transftensor(epsrTensNoSym3D, TransfTens)
!
!  The indices for the multiplication are as:
!  E(i,j,m) = C(i,j,k,l)*D(k,l,m)
        !call multtensor(CTensNoSym3D(1:2,1:2,1:2,1:2), DTensNoSym3D(1:2,1:2,1:2), ETensNoSym3D(1:2,1:2,1:2))
!        call multtensor(DTensNoSym3D(1:2,1:2,1:2), CTensNoSym3D(1:2,1:2,1:2,1:2), ETensNoSym3D(1:2,1:2,1:2))
!
!  The transpose of E can be used for the Electro-Elastic tensor
!  which has the opposite effect of the Elasto-Electric tensor
!
!  E_transpose(m,i,j) = D(k,l,m)*C(i,j,k,l)
!                ---      ---      ---,---
!                Sym      Sym      Sym,Sym
!
!                |E_t11 E_t12 E_t13|   |E11 E21 E31|
!  E_transpose = |-----------------| = |-----------|
!                |E_t21 E_t22 E_t23|   |E12 E22 E32|
!
!  nu is of the form nu(2,2,nnat,nnat)
!  submatrix 1-1, x <- x
        nu(1,1,1,1) = cmplx(CTensNoSym3D(1,1,1,1),0._DP,DPC)
        nu(1,2,1,1) = cmplx(CTensNoSym3D(1,1,1,2),0._DP,DPC)
        nu(2,1,1,1) = cmplx(CTensNoSym3D(1,2,1,1),0._DP,DPC)
        nu(2,2,1,1) = cmplx(CTensNoSym3D(1,2,1,2),0._DP,DPC)
!  submatrix 1-2, x <- y
        nu(1,1,1,2) = cmplx(CTensNoSym3D(1,1,2,1),0._DP,DPC)
        nu(1,2,1,2) = cmplx(CTensNoSym3D(1,1,2,2),0._DP,DPC)
        nu(2,1,1,2) = cmplx(CTensNoSym3D(1,2,2,1),0._DP,DPC)
        nu(2,2,1,2) = cmplx(CTensNoSym3D(1,2,2,2),0._DP,DPC)
!  submatrix 1-3, x <- el. potential
        nu(1,1,1,3) = cmplx(ETensNoSym3D(1,1,1),0._DP,DPC)
        nu(1,2,1,3) = cmplx(ETensNoSym3D(1,1,2),0._DP,DPC)
        nu(2,1,1,3) = cmplx(ETensNoSym3D(1,2,1),0._DP,DPC)
        nu(2,2,1,3) = cmplx(ETensNoSym3D(1,2,2),0._DP,DPC)
!  submatrix 2-1, y <- x
        nu(1,1,2,1) = cmplx(CTensNoSym3D(2,1,1,1),0._DP,DPC)
        nu(1,2,2,1) = cmplx(CTensNoSym3D(2,1,1,2),0._DP,DPC)
        nu(2,1,2,1) = cmplx(CTensNoSym3D(2,2,1,1),0._DP,DPC)
        nu(2,2,2,1) = cmplx(CTensNoSym3D(2,2,1,2),0._DP,DPC)
!  submatrix 2-2, y <- y
        nu(1,1,2,2) = cmplx(CTensNoSym3D(2,1,2,1),0._DP,DPC)
        nu(1,2,2,2) = cmplx(CTensNoSym3D(2,1,2,2),0._DP,DPC)
        nu(2,1,2,2) = cmplx(CTensNoSym3D(2,2,2,1),0._DP,DPC)
        nu(2,2,2,2) = cmplx(CTensNoSym3D(2,2,2,2),0._DP,DPC)
!  submatrix 2-3, y <- el. potential
        nu(1,1,2,3) = cmplx(ETensNoSym3D(2,1,1),0._DP,DPC)
        nu(1,2,2,3) = cmplx(ETensNoSym3D(2,1,2),0._DP,DPC)
        nu(2,1,2,3) = cmplx(ETensNoSym3D(2,2,1),0._DP,DPC)
        nu(2,2,2,3) = cmplx(ETensNoSym3D(2,2,2),0._DP,DPC)
!  submatrix 3-1, el. potential <- x
        nu(1,1,3,1) = cmplx(ETensNoSym3D(1,1,1),0._DP,DPC)
        nu(1,2,3,1) = cmplx(ETensNoSym3D(1,2,1),0._DP,DPC)
        nu(2,1,3,1) = cmplx(ETensNoSym3D(1,1,2),0._DP,DPC)
        nu(2,2,3,1) = cmplx(ETensNoSym3D(1,2,2),0._DP,DPC)
!  submatrix 3-2, el. potential <- y
        nu(1,1,3,2) = cmplx(ETensNoSym3D(2,1,1),0._DP,DPC)
        nu(1,2,3,2) = cmplx(ETensNoSym3D(2,2,1),0._DP,DPC)
        nu(2,1,3,2) = cmplx(ETensNoSym3D(2,1,2),0._DP,DPC)
        nu(2,2,3,2) = cmplx(ETensNoSym3D(2,2,2),0._DP,DPC)
!  submatrix 3-3, el. potential <- el. potential
        nu(1,1,3,3) = cmplx(-epsrTensNoSym3D(1,1),0._DP,DPC)
        nu(1,2,3,3) = cmplx(-epsrTensNoSym3D(1,2),0._DP,DPC)
        nu(2,1,3,3) = cmplx(-epsrTensNoSym3D(2,1),0._DP,DPC)
        nu(2,2,3,3) = cmplx(-epsrTensNoSym3D(2,2),0._DP,DPC)

!  gamma is of the form gamma(2,nnat,nnat)
        gamma(:,:,:) = 0._DPC

!  beta is of the form beta(2,nnat,nnat)
        beta(:,:,:) = 0._DPC
!
!  alpha is of the form alpha(nnat,nnat)
        alpha(:,:) = 0._DPC
!  submatrix 1-1, x <- x
        alpha(1,1) = cmplx(-((omega**2)*list(1)),0._DP,DPC)
!  submatrix 2-2, y <- y
        alpha(2,2) = cmplx(-((omega**2)*list(1)),0._DP,DPC)
!
!  f1 is of the form f1(nnat)
!  x
        f1(1) = 0._DPC
!  y
        f1(2) = 0._DPC
!  el. potential
        f1(3) = 0._DPC

!  f2 is of the form f2(2,nnat)
!  x
        f2(1,1) = 0._DPC
        f2(2,1) = 0._DPC
!  y
        f2(1,2) = 0._DPC
        f2(2,2) = 0._DPC
!  el. potential
        f2(1,3) = 0._DPC
        f2(2,3) = 0._DPC

        deallocate(list)
!
      case ('PIEZOPYROELEC')
!  nnat=4 in the 2D Piezo-pyro-electric problem
!  in the first nature we consider the displacement in x
!  in the second nature we consider the displacement in y
!  in the third nature we consider the electric potential
!  in the fourth nature we consider the temperature
!
        TransfTens = reshape((/list( 7),list( 8),list( 9)   &
                          &  , list(10),list(11),list(12)   &
                          &  , list(13),list(14),list(15)/) &
                          &  ,  (/3,3/), (/0._DP/), (/2,1/))
        LamTensNoSym3D = reshape((/list(16),  0._DP ,  0._DP     &
                          &  ,       0._DP ,list(17),  0._DP     &
                          &  ,       0._DP ,  0._DP ,list(18)/)  &
                          &  ,       (/3,3/), (/0._DP/), (/2,1/))
        epsrTensNoSym3D = reshape((/eps0*list(19),  0._DP ,  0._DP    &
                          &  ,        0._DP ,eps0*list(20),  0._DP    &
                          &  ,        0._DP ,  0._DP ,eps0*list(21)/) &
                          &  ,        (/3,3/), (/0._DP/), (/2,1/))
        CTensNoSym3D = reshape((/list(22),list(27),list(26),list(27),list(23),list(25),list(26),list(25),list(24)   &
                          &  ,   list(27),list(42),list(41),list(42),list(32),list(39),list(41),list(39),list(36)   &
                          &  ,   list(26),list(41),list(40),list(41),list(31),list(38),list(40),list(38),list(35)   &
                          &  ,   list(27),list(42),list(41),list(42),list(32),list(39),list(41),list(39),list(36)   &
                          &  ,   list(23),list(32),list(31),list(32),list(28),list(30),list(31),list(30),list(29)   &
                          &  ,   list(25),list(39),list(38),list(39),list(30),list(37),list(38),list(37),list(34)   &
                          &  ,   list(26),list(41),list(40),list(41),list(31),list(38),list(40),list(38),list(35)   &
                          &  ,   list(25),list(39),list(38),list(39),list(30),list(37),list(38),list(37),list(34)   &
                          &  ,   list(24),list(36),list(35),list(36),list(29),list(34),list(35),list(34),list(33)/) &
                          &  ,    (/3,3,3,3/), (/0._DP/), (/4,3,2,1/))
        ETensNoSym3D = reshape((/list(43),list(44),list(45)   &
                          &  ,   list(58),list(59),list(60)   &
                          &  ,   list(55),list(56),list(57)   &
                          &  ,   list(58),list(59),list(60)   &
                          &  ,   list(46),list(47),list(48)   &
                          &  ,   list(52),list(53),list(54)   &
                          &  ,   list(55),list(56),list(57)   &
                          &  ,   list(52),list(53),list(54)   &
                          &  ,   list(49),list(50),list(51)/) &
                          &  ,    (/3,3,3/), (/0._DP/), (/3,2,1/))
        ATensNoSym3D = reshape((/list(61),list(62),list(63)   &
                          &  ,   list(64),list(65),list(66)   &
                          &  ,   list(67),list(68),list(69)/) &
                          &  ,    (/3,3/), (/0._DP/), (/2,1/))
        PTensNoSym3D = (/list(70)   &
                  &  ,   list(71)   &
                  &  ,   list(72)   &
                  &      /)

        call transftensor(LamTensNoSym3D, TransfTens)
        call transftensor(epsrTensNoSym3D, TransfTens)
        call transftensor(CTensNoSym3D, TransfTens)
        call transftensor(ETensNoSym3D, TransfTens)
        call transftensor(ATensNoSym3D, TransfTens)
        call transftensor(PTensNoSym3D, TransfTens)

        call multtensor(CTensNoSym3D(1:2,1:2,1:2,1:2), ATensNoSym3D(1:2,1:2), MTensNoSym3D(1:2,1:2))
!
!  Modify coefficients for PML layer
        rho = cmplx(list(1),0._DP,DPC)
        pml  = int(list(73))
        if (pml .eq. 1) then
          vpml = (/list(74),list(75),list(76)/)
          pmldir(1) = cmplx(1._DP,0._DP,DPC) + (vpml(1) * cmplx(list(77),-list(78),DPC))
          pmldir(2) = cmplx(1._DP,0._DP,DPC) + (vpml(2) * cmplx(list(77),-list(78),DPC))
          pmlprod = pmldir(1) * pmldir(2)
!  Update stiffness matrix coefficients
          do j=1,2
            do l=1,2
              CTensNoSym3D(:,j,:,l) = CTensNoSym3D(:,j,:,l) * pmlprod / (pmldir(j)*pmldir(l))
            end do
          end do
!  Update piezoelectric coefficients
          do j=1,2
            do k=1,2
              ETensNoSym3D(:,j,k) = ETensNoSym3D(:,j,k) * pmlprod / (pmldir(j)*pmldir(k))
            end do
          end do
!  Update thermal conductivity coefficients
          do i=1,2
            do j=1,2
              LamTensNoSym3D(i,j) = LamTensNoSym3D(i,j) * pmlprod / (pmldir(i)*pmldir(j))
            end do
          end do
!  Update permittivity coefficients
          do i=1,2
            do j=1,2
              epsrTensNoSym3D(i,j) = epsrTensNoSym3D(i,j) * pmlprod / (pmldir(i)*pmldir(j))
            end do
          end do
!  Update density
          rho = rho*pmlprod
        end if
!
!  nu is of the form nu(2,2,nnat,nnat)
!  submatrix 1-1, x <- x
        nu(1,1,1,1) = cmplx(CTensNoSym3D(1,1,1,1),0._DP,DPC)
        nu(1,2,1,1) = cmplx(CTensNoSym3D(1,1,1,2),0._DP,DPC)
        nu(2,1,1,1) = cmplx(CTensNoSym3D(1,2,1,1),0._DP,DPC)
        nu(2,2,1,1) = cmplx(CTensNoSym3D(1,2,1,2),0._DP,DPC)
!  submatrix 1-2, x <- y
        nu(1,1,1,2) = cmplx(CTensNoSym3D(1,1,2,1),0._DP,DPC)
        nu(1,2,1,2) = cmplx(CTensNoSym3D(1,1,2,2),0._DP,DPC)
        nu(2,1,1,2) = cmplx(CTensNoSym3D(1,2,2,1),0._DP,DPC)
        nu(2,2,1,2) = cmplx(CTensNoSym3D(1,2,2,2),0._DP,DPC)
!  submatrix 1-3, x <- el. potential
        nu(1,1,1,3) = cmplx(ETensNoSym3D(1,1,1),0._DP,DPC)
        nu(1,2,1,3) = cmplx(ETensNoSym3D(1,1,2),0._DP,DPC)
        nu(2,1,1,3) = cmplx(ETensNoSym3D(1,2,1),0._DP,DPC)
        nu(2,2,1,3) = cmplx(ETensNoSym3D(1,2,2),0._DP,DPC)
!  submatrix 1-4, x <- temperature
        nu(:,:,1,4) = 0._DPC
!  submatrix 2-1, y <- x
        nu(1,1,2,1) = cmplx(CTensNoSym3D(2,1,1,1),0._DP,DPC)
        nu(1,2,2,1) = cmplx(CTensNoSym3D(2,1,1,2),0._DP,DPC)
        nu(2,1,2,1) = cmplx(CTensNoSym3D(2,2,1,1),0._DP,DPC)
        nu(2,2,2,1) = cmplx(CTensNoSym3D(2,2,1,2),0._DP,DPC)
!  submatrix 2-2, y <- y
        nu(1,1,2,2) = cmplx(CTensNoSym3D(2,1,2,1),0._DP,DPC)
        nu(1,2,2,2) = cmplx(CTensNoSym3D(2,1,2,2),0._DP,DPC)
        nu(2,1,2,2) = cmplx(CTensNoSym3D(2,2,2,1),0._DP,DPC)
        nu(2,2,2,2) = cmplx(CTensNoSym3D(2,2,2,2),0._DP,DPC)
!  submatrix 2-3, y <- el. potential
        nu(1,1,2,3) = cmplx(ETensNoSym3D(2,1,1),0._DP,DPC)
        nu(1,2,2,3) = cmplx(ETensNoSym3D(2,1,2),0._DP,DPC)
        nu(2,1,2,3) = cmplx(ETensNoSym3D(2,2,1),0._DP,DPC)
        nu(2,2,2,3) = cmplx(ETensNoSym3D(2,2,2),0._DP,DPC)
!  submatrix 2-4, y <- temperature
        nu(:,:,2,4) = 0._DPC
!  submatrix 3-1, el. potential <- x
        nu(1,1,3,1) = cmplx(ETensNoSym3D(1,1,1),0._DP,DPC)
        nu(1,2,3,1) = cmplx(ETensNoSym3D(1,2,1),0._DP,DPC)
        nu(2,1,3,1) = cmplx(ETensNoSym3D(1,1,2),0._DP,DPC)
        nu(2,2,3,1) = cmplx(ETensNoSym3D(1,2,2),0._DP,DPC)
!  submatrix 3-2, el. potential <- y
        nu(1,1,3,2) = cmplx(ETensNoSym3D(2,1,1),0._DP,DPC)
        nu(1,2,3,2) = cmplx(ETensNoSym3D(2,2,1),0._DP,DPC)
        nu(2,1,3,2) = cmplx(ETensNoSym3D(2,1,2),0._DP,DPC)
        nu(2,2,3,2) = cmplx(ETensNoSym3D(2,2,2),0._DP,DPC)
!  submatrix 3-3, el. potential <- el. potential
        nu(1,1,3,3) = cmplx(-epsrTensNoSym3D(1,1),0._DP,DPC)
        nu(1,2,3,3) = cmplx(-epsrTensNoSym3D(1,2),0._DP,DPC)
        nu(2,1,3,3) = cmplx(-epsrTensNoSym3D(2,1),0._DP,DPC)
        nu(2,2,3,3) = cmplx(-epsrTensNoSym3D(2,2),0._DP,DPC)
!  submatrix 3-4, el. potential <- temperature
        nu(:,:,3,4) = 0._DPC
!  submatrix 4-1, temperature <- x
        nu(:,:,4,1) = 0._DPC
!  submatrix 4-2, temperature <- y
        nu(:,:,4,2) = 0._DPC
!  submatrix 4-3, temperature <- el. potential
        nu(:,:,4,3) = 0._DPC
!  submatrix 4-4, temperature <- temperature
        nu(1,1,4,4) = cmplx(LamTensNoSym3D(1,1),0._DP,DPC)
        nu(1,2,4,4) = cmplx(LamTensNoSym3D(1,2),0._DP,DPC)
        nu(2,1,4,4) = cmplx(LamTensNoSym3D(2,1),0._DP,DPC)
        nu(2,2,4,4) = cmplx(LamTensNoSym3D(2,2),0._DP,DPC)

!  gamma is of the form gamma(2,nnat,nnat)
        gamma(:,:,:) = 0._DPC
!  submatrix 1-4, x <- temperature
        gamma(1,1,4) = cmplx(-MTensNoSym3D(1,1),0._DP,DPC)
        gamma(2,1,4) = cmplx(-MTensNoSym3D(1,2),0._DP,DPC)
!  submatrix 2-4, y <- temperature
        gamma(1,2,4) = cmplx(-MTensNoSym3D(2,1),0._DP,DPC)
        gamma(2,2,4) = cmplx(-MTensNoSym3D(2,2),0._DP,DPC)
!  submatrix 3-4, el. potential <- temperature
        gamma(1,3,4) = cmplx(PTensNoSym3D(1),0._DP,DPC)
        gamma(2,3,4) = cmplx(PTensNoSym3D(2),0._DP,DPC)
        
!  beta is of the form beta(2,nnat,nnat)
        beta(:,:,:) = 0._DPC
!  submatrix 4-1, temperature <- x
        beta(1,4,1) = cmplx(0._DP,-omega*MTensNoSym3D(1,1)*list(6),DPC)
        beta(2,4,1) = cmplx(0._DP,-omega*MTensNoSym3D(1,2)*list(6),DPC)
!  submatrix 4-2, temperature <- y
        beta(1,4,2) = cmplx(0._DP,-omega*MTensNoSym3D(2,1)*list(6),DPC)
        beta(2,4,2) = cmplx(0._DP,-omega*MTensNoSym3D(2,2)*list(6),DPC)
!  submatrix 4-3, temperature <- el. potential
        beta(1,4,3) = cmplx(0._DP,omega*PTensNoSym3D(1)*list(6),DPC)
        beta(2,4,3) = cmplx(0._DP,omega*PTensNoSym3D(2)*list(6),DPC)

!  alpha is of the form alpha(nnat,nnat)
        alpha(:,:) = 0._DPC
!  submatrix 1-1, x <- x
        alpha(1,1) = cmplx(-(omega**2),0._DP,DPC) * rho
!  submatrix 2-2, y <- y
        alpha(2,2) = cmplx(-(omega**2),0._DP,DPC) * rho
!  submatrix 4-4, temperature <- temperature
        alpha(4,4) = cmplx(0._DP,omega*list(4),DPC) * rho

!  f1 is of the form f1(nnat)
!  x
        f1(1) = 0._DPC
!  y
        f1(2) = 0._DPC
!  el. potential
        f1(3) = 0._DPC
!  temperature
        f1(4) = cmplx(list(5),0._DP,DPC)

!  f2 is of the form f2(2,nnat)
!  x
        f2(1,1) = 0._DPC
        f2(2,1) = 0._DPC
!  y
        f2(1,2) = 0._DPC
        f2(2,2) = 0._DPC
!  el. potential
        f2(1,3) = 0._DPC
        f2(2,3) = 0._DPC
!  temperature
        f2(1,4) = 0._DPC
        f2(2,4) = 0._DPC

        deallocate(list)
!
      end select
!
      return
      end subroutine pdecoeff
!
!
!
      pure subroutine calctensor_c(var11, var22, angle, invert, ten)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: angle
      complex (DPC) :: var11, var22, ten(2,2)
      logical :: invert
      intent (in) :: angle, var11, var22, invert
      intent (out) :: ten
!
!  Programm calculates tensor with non-diagonal entries from x- and y-value
!  with given principal axis angle. Tensor can be inverted after calculation.
!
!  Input:
!            angle    pricipal axis angle
!            var11    x-value of material parameter
!            var22    y-value of material parameter
!  Output:
!            ten      tensor containing calculated parameters
!                     (have to be assigned to nu at the correct positions)
!
!  local variables
      real (DP) :: si, co, s2, c2
      complex (DPC) :: vartemp11, vartemp22
!
!            si       sine of angle
!            co       cosine of angle
!            s2       sin(angle)*sin(angle)
!            c2       cos(angle)*cos(angle)
!            vartemp  variables used for temporary storage of var11 and var22
!
!  invert tensor, if invert is .true.
      if (invert) then
        vartemp11 = 1._DP / var11
        vartemp22 = 1._DP / var22
      else
        vartemp11 = var11
        vartemp22 = var22
      end if
!
      si = sin(angle)
      co = cos(angle)
      if (si .lt. co) then
        s2 = si*si
        c2 = 1._DP-s2
      else
        c2 = co*co
        s2 = 1._DP-c2
      end if
!
      ten(1,2) = si*co*(vartemp11 - vartemp22)
      ten(2,1) = ten(1,2)
      ten(1,1) = c2*vartemp11 + s2*vartemp22
      ten(2,2) = s2*vartemp11 + c2*vartemp22
!
      return
      end subroutine calctensor_c
      
      
      pure subroutine calctensor_d(var11, var22, angle, invert, ten)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: angle
      real (DP) :: var11, var22, ten(2,2)
      logical :: invert
      intent (in) :: angle, var11, var22, invert
      intent (out) :: ten
!
!  Programm calculates tensor with non-diagonal entries from x- and y-value
!  with given principal axis angle. Tensor can be inverted after calculation.
!
!  Input:
!            angle    pricipal axis angle
!            var11    x-value of material parameter
!            var22    y-value of material parameter
!  Output:
!            ten      tensor containing calculated parameters
!                     (have to be assigned to nu at the correct positions)
!
!  local variables
      real (DP) :: si, co, s2, c2
      real (DP) :: vartemp11, vartemp22
!
!            si       sine of angle
!            co       cosine of angle
!            s2       sin(angle)*sin(angle)
!            c2       cos(angle)*cos(angle)
!            vartemp  variables used for temporary storage of var11 and var22
!
!  invert tensor, if invert is .true.
      if (invert) then
        vartemp11 = 1._DP / var11
        vartemp22 = 1._DP / var22
      else
        vartemp11 = var11
        vartemp22 = var22
      end if
!
      si = sin(angle)
      co = cos(angle)
      if (si .lt. co) then
        s2 = si*si
        c2 = 1._DP-s2
      else
        c2 = co*co
        s2 = 1._DP-c2
      end if
!
      ten(1,2) = si*co*(vartemp11 - vartemp22)
      ten(2,1) = ten(1,2)
      ten(1,1) = c2*vartemp11 + s2*vartemp22
      ten(2,2) = s2*vartemp11 + c2*vartemp22
!
      return
      end subroutine calctensor_d


      pure subroutine transftensor_3Drank1(tensor1, transmat)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: tensor1(3)
      real (DP) :: transmat(3,3)
      intent (in) :: transmat
      intent (inout) :: tensor1
!
!  Program transforms a 3D rank-1 tensor
!  according to the transformation matrix transmat
!
!  Input:
!            tensor1  the rank-1 tensor
!            transmat the transformation matrix
!
!  Output:
!            tensor1  the transformed rank-1 tensor
!
!  local variables
      integer (I4B) m, i
      real (DP) :: temp(1)
      real (DP) :: TempTens(3)
!
      TempTens = 0._DP
      do m=1,3
            temp(1) = 0._DP
            do i=1,3
                temp(1) = temp(1) + transmat(m,i) * tensor1(i)
            end do
            TempTens(m) = temp(1)
      end do
      tensor1 = TempTens
!
      return
      end subroutine transftensor_3Drank1


      pure subroutine transftensor_3Drank2(tensor1, transmat)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: tensor1(3,3)
      real (DP) :: transmat(3,3)
      intent (in) :: transmat
      intent (inout) :: tensor1
!
!  Program transforms a 3D rank-2 tensor
!  according to the transformation matrix transmat
!
!  Input:
!            tensor1  the rank-2 tensor
!            transmat the transformation matrix
!
!  Output:
!            tensor1  the transformed rank-2 tensor
!
!  local variables
      integer (I4B) m, n, i, j
      real (DP) :: temp(2)
      real (DP) :: TempTens(3,3)
!
      TempTens = 0._DP
      do m=1,3
        do n=1,3
            temp(1) = 0._DP
            do i=1,3
              temp(2) = 0._DP
              do j=1,3
                temp(2) = temp(2) + transmat(n,j) * tensor1(i,j)
              end do
              temp(1) = temp(1) + transmat(m,i) * temp(2)
            end do
            TempTens(m,n) = temp(1)
        end do
      end do
      tensor1 = TempTens
!
      return
      end subroutine transftensor_3Drank2


      pure subroutine transftensor_3Drank3(tensor1, transmat)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: tensor1(3,3,3)
      real (DP) :: transmat(3,3)
      intent (in) :: transmat
      intent (inout) :: tensor1
!
!  Program transforms a 3D rank-3 tensor
!  according to the transformation matrix transmat
!
!  Input:
!            tensor1  the rank-3 tensor
!            transmat the transformation matrix
!
!  Output:
!            tensor1  the transformed rank-3 tensor
!
!  local variables
      integer (I4B) m, n, o, i, j, k
      real (DP) :: temp(3)
      real (DP) :: TempTens(3,3,3)
!
      TempTens = 0._DP
      do m=1,3
        do n=1,3
          do o=1,3
            temp(1) = 0._DP
            do i=1,3
              temp(2) = 0._DP
              do j=1,3
                temp(3) = 0._DP
                do k=1,3
                  temp(3) = temp(3) + transmat(o,k) * tensor1(i,j,k)
                end do
                temp(2) = temp(2) + transmat(n,j) * temp(3)
              end do
              temp(1) = temp(1) + transmat(m,i) * temp(2)
            end do
            TempTens(m,n,o) = temp(1)
          end do
        end do
      end do
      tensor1 = TempTens
!
      return
      end subroutine transftensor_3Drank3


      pure subroutine transftensor_3Drank4(tensor1, transmat)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: tensor1(3,3,3,3)
      real (DP) :: transmat(3,3)
      intent (in) :: transmat
      intent (inout) :: tensor1
!
!  Program transforms a 3D rank-4 tensor
!  according to the transformation matrix transmat
!
!  Input:
!            tensor1  the rank-4 tensor
!            transmat the transformation matrix
!
!  Output:
!            tensor1  the transformed rank-4 tensor
!
!  local variables
      integer (I4B) m, n, o, p, i, j, k, l
      real (DP) :: temp(4)
      real (DP) :: TempTens(3,3,3,3)
!
      TempTens = 0._DP
      do m=1,3
        do n=1,3
          do o=1,3
            do p=1,3
              temp(1) = 0._DP
              do i=1,3
                temp(2) = 0._DP
                do j=1,3
                  temp(3) = 0._DP
                  do k=1,3
                    temp(4) = 0._DP
                    do l=1,3
                      temp(4) = temp(4) + transmat(p,l) * tensor1(i,j,k,l)
                    end do
                    temp(3) = temp(3) + transmat(o,k) * temp(4)
                  end do
                  temp(2) = temp(2) + transmat(n,j) * temp(3)
                end do
                temp(1) = temp(1) + transmat(m,i) * temp(2)
              end do
              TempTens(m,n,o,p) = temp(1)
            end do
          end do
        end do
      end do
      tensor1 = TempTens
!
      return
      end subroutine transftensor_3Drank4


      pure subroutine multtensor_ECD(tensor1, tensor2, ten)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: tensor1(2,2,2,2), tensor2(2,2,2), ten(2,2,2)
      intent (in) :: tensor1, tensor2
      intent (out) :: ten
!
!  Program multiplies 2 tensors of rank-4 and -3
!  and produces a rank-3 tensor as in:
!  rank-3 = rank-4 * rank-3 (2 common indices)
!
!  Input:
!            tensor1  the rank-4 tensor
!            tensor2  the rank-3 tensor
!
!  Output:
!            ten      the resulting rank-3 tensor
!
!  local variables
      integer (I4B) i, j, k, l, m
!
      ten = 0._DP
      do i=1,2
        do j=1,2
          do k=1,2
            do l=1,2
              do m=1,2

                ten(i,j,m) = ten(i,j,m) + tensor1(i,j,k,l)*tensor2(k,l,m)

              end do
            end do
          end do
        end do
      end do
!
      return
      end subroutine multtensor_ECD


      pure subroutine multtensor_EDC(tensor1, tensor2, ten)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: tensor1(2,2,2), tensor2(2,2,2,2), ten(2,2,2)
      intent (in) :: tensor1, tensor2
      intent (out) :: ten
!
!  Program multiplies 2 tensors of rank-4 and -3
!  and produces a rank-3 tensor as in:
!  rank-3 = rank-4 * rank-3 (2 common indices)
!
!  Input:
!            tensor1  the rank-4 tensor
!            tensor2  the rank-3 tensor
!
!  Output:
!            ten      the resulting rank-3 tensor
!
!  local variables
      integer (I4B) i, j, k, l, m
!
      ten = 0._DP
      do i=1,2
        do j=1,2
          do k=1,2
            do l=1,2
              do m=1,2

                ten(m,k,l) = ten(m,k,l) + tensor1(i,j,m)*tensor2(i,j,k,l)

              end do
            end do
          end do
        end do
      end do
!
      return
      end subroutine multtensor_EDC


      pure subroutine multtensor_MCA(tensor1, tensor2, ten)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: tensor1(2,2,2,2), tensor2(2,2), ten(2,2)
      intent (in) :: tensor1, tensor2
      intent (out) :: ten
!
!  Program multiplies 2 tensors of rank-4 and -2
!  and produces a rank-2 tensor as in:
!  rank-2 = rank-4 * rank-2 (2 common indices)
!
!  Input:
!            tensor1  the rank-4 tensor
!            tensor2  the rank-2 tensor
!
!  Output:
!            ten      the resulting rank-2 tensor
!
!  local variables
      integer (I4B) i, j, k, l
!
      ten = 0._DP
      do i=1,2
        do j=1,2
          do k=1,2
            do l=1,2

              ten(i,j) = ten(i,j) + tensor1(i,j,k,l)*tensor2(k,l)

            end do
          end do
        end do
      end do
!
      return
      end subroutine multtensor_MCA