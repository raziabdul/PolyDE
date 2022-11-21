      subroutine pdecoefftrans(elem,xs,ys,alpha)
      use feminterface, only: calctensor, multtensor, transftensor, fetchmatparameters, getsetting, invertmat, &
                              invertmat3, trans_matrix, transten_r2
      use femtypes
      use globalvariables
      use matconstants
      use mpmodule
      implicit none
      integer (I4B) elem
      real(DP) xs, ys
      complex (DPC) :: alpha(:,:)
      intent (in) :: elem, xs, ys
      intent (out) :: alpha
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
!    $Date: 2014/07/03 12:53:31 $
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
!            alpha    coefficient used for the construction of the mass matrix
!
!  local variables
      integer (I4B) matindex, pml
      integer (I4B) i, j, k, l
      real(DP), allocatable :: list(:)
      real (DP) :: angle, phase, losst, scalar, xys(2)
      real (DP) :: jmx, jmy, gam(2,2), v11, v22
      real (DP) :: Lame1, Lame2
      real (DP) :: Tref
      real (DP) :: TransfTens(3,3)
!  3D Tensors
      real (DP) :: CTensNoSym3D(3,3,3,3), ETensNoSym3D(3,3,3), epsrTensNoSym3D(3,3)
      real (DP) :: LamTensNoSym3D(3,3), ATensNoSym3D(3,3), PTensNoSym3D(3)
      real (DP) :: MTensNoSym3D(3,3), DTensNoSym3D(3,3,3)
!  2D Tensors
      real (DP) :: CTensNoSym(2,2,2,2), DTensNoSym(2,2,2), epsrTensNoSym(2,2)
      real (DP) :: ATensNoSym(2,2), PTensNoSym(2)
      real (DP) :: MTensNoSym(2,2), ETensNoSym(2,2,2), EtrTensNoSym(2,2,2)
      !  For gravitational expansion
      real (DP) :: GTensNoSym(2,2,2)
      !  Stress-gravity tensor
      real (DP) :: LTensNoSym(2,2,2)
      complex (DPC) :: var11, var22, zcomp, source1, source2, matten(3,3)
      logical :: invert
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
      call fetchmatparameters(list,matindex,xys)
!
!  Different physics modes
      select case (physics)
      case ('PLANE STRAIN')

!  alpha is of the form alpha(nnat,nnat)
        alpha(1,1) = cmplx(list(1),0._DP,DPC)
        alpha(1,2) = 0._DPC
        alpha(2,1) = 0._DPC
        alpha(2,2) = cmplx(list(1),0._DP,DPC)

        deallocate (list)
!        
      end select
!
      return
      end subroutine pdecoefftrans