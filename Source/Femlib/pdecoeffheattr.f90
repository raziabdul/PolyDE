      subroutine pdecoeffheattr(elem,xs,ys,nu,gamma,alpha,beta,f1,f2)
      use feminterface, only: calctensor, multtensor, transftensor, fetchmatparameters, getsetting, invertmat, &
                              invertmat3, trans_matrix, transten_r2
      use femtypes
      use globalvariables
      use matconstants
      use mpmodule
      implicit none
      integer (I4B) elem
      real(DP) xs, ys
      complex (DPC) :: nu(:,:,:,:), gamma(:,:,:)
      complex (DPC) :: alpha(:,:), beta(:,:,:), f1(:), f2(:,:)
      intent (in) :: elem, xs, ys
      intent (out) :: nu, gamma, alpha, beta, f1, f2
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
!    $Revision: 1.2 $
!    $Date: 2014/07/03 12:53:56 $
!    $Author: m_kasper $
!
!  deliver the material coefficients in the form of the general PDE:
!  
!        -div(nu*grad(u)) -div(gamma*u) + (beta*grad(u)) + alpha*u = f1 + div(f2)
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
!            f1       right hand side constant
!            f2       right hand side source vector
!
!  local variables
      integer (I4B) matindex
      real(DP), allocatable :: list(:)
      real (DP) :: angle, phase, scalar, xys(2)
!  2D Tensors
      complex (DPC) :: var11, var22, zcomp, source1, source2
      logical :: invert
      real (DP), parameter :: rad = pi/180._DP
!
!
!  determine material index for given element in region geb
      matindex=matzif(geb(elem))

!  read parameters from internal list and assign them to local variables
      allocate (list(numparam))
      xys=(/xs,ys/)
      call fetchmatparameters(list,matindex,xys)

!  select case construct to assign pde coefficients
      select case (physics)
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
      case ('THERMOELECTRIC')
!  nnat=2 in the thermoelectric problem
!  the thermal is considered as the first nature
!  the electric is considered as the second nature
!
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

!  gamma is of the form gamma(2,nnat,nnat)
        gamma(:,:,:) = 0._DPC
!  beta is of the form beta(2,nnat,nnat)
        beta(1,1,1) = 0._DPC
        beta(2,1,1) = 0._DPC
        beta(1,1,2) = 0._DPC
        beta(2,1,2) = 0._DPC
        beta(1,2,1) = 0._DPC
        beta(2,2,1) = 0._DPC
        beta(1,2,2) = 0._DPC
        beta(2,2,2) = 0._DPC
!  alpha is of the form alpha(nnat,nnat)
        alpha(1,1) = 0._DPC
        alpha(1,2) = 0._DPC
        alpha(2,1) = 0._DPC
        alpha(2,2) = 0._DPC
!  f1 is of the form f1(nnat)
        f1(1) = source1
        f1(2) =  cmplx(list(6)*list(18),0._DP,DPC)
!  f2 is of the form f2(2,nnat)
        f2(1,1) = 0._DPC
        f2(2,1) = 0._DPC
        f2(1,2) = cmplx(-list(6)*list(19),0._DP,DPC)
        f2(2,2) = cmplx(-list(6)*list(20),0._DP,DPC)
        deallocate(list)
!
      end select
!
      return
      end subroutine pdecoeffheattr