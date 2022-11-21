      subroutine userbc(branch,bcnum,xs,ys,nature,pval,qval,elem)
      use feminterface, only: xy2lam, fieldu
      use femtypes
      use globalvariables, only: mu0, eps0, omega, physics, pi, &
                               & zki,xbk, ybk, userparam, xn, yn, e, nnat
      implicit none
      integer (I4B) branch, bcnum, nature
      integer (I4B) , optional :: elem
      real (DP) xs, ys
      complex(DPC) pval
      complex(DPC), optional :: qval(:)
      intent (in) :: branch, bcnum, xs, ys, nature, elem
      intent (out) :: pval, qval
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
!    $Revision: 1.25 $
!    $Date: 2014/07/03 13:21:04 $
!    $Author: m_kasper $
!
!
!  User programming interface for non-constant Dirichlet boundary conditions
!  Return value (of the dirichlet bc) for coordinates xs, ys on a branch. 
!  Case selector bcnum (1..99) indicates the user function according to the
!  layer designation USER<bcnum> in the Autocad drawing (layer of the geometry).
!
!  Input:
!        branch      the branch for which the BC is sought 
!        bcnum       case selector from drawing layer text, i.e. the boundary condition
!        xs, ys      coordinates for determining dirichlet
!        nature      index of the nature (multiphysics)
!  Output:
!    in the case of dirichlet BC:   bcnum (1..99)
!        pval        value for the scalar potential   u=pval
!    in the case of general BC:     bcnum (201..299)
!        pval
!        qval        (nu*grad u + f2)*n = pval + qval*u
!
!  local variables:
      real (DP) s, startp(2), endp(2)
      real (DP) dir(2), length, dirnormal(2)
!-----------------
      real (DP) w, beta, gammap, gamman, k0, kappa, neff, nwg, np, nn
!-----------------
      real (DP) lambda(3)
      complex (DPC) u(nnat)
!-----------------
!
!        s           local coordinate [0..1] along the branch
!        startp      starting point of the branch (startp(1) = x, startp(2) = y)
!        endp        end point of the branch (endp(1) = x, endp(2) = y)
!        dir         direction vector from startp to endp
!        length      distance between startp and endp
!        dirnormal   normalized direction vector dir
!
!  define starting and end points of the branch.
      if (zki(1,branch) .gt. zki(2,branch)) then
!  has to be evaluated for hand-written geometry
        startp(1) = xbk(zki(2,branch))
        startp(2) = ybk(zki(2,branch))
        endp(1) = xbk(zki(1,branch))
        endp(2) = ybk(zki(1,branch))
      else
!  normal sorting order if geometry written by acadnet
        startp(1) = xbk(zki(1,branch))
        startp(2) = ybk(zki(1,branch))
        endp(1) = xbk(zki(2,branch))
        endp(2) = ybk(zki(2,branch))
      end if
!
!  calculate direction vector
      dir = endp - startp
      length = sqrt((dir(1))**2+(dir(2))**2)
      dirnormal = dir / length
!
!  check if branches are lines or arcs.
      if (abs(zki(3,branch)) == 0) then
        !  check to see if line is horizontal, vertical or diagonal
        if (abs(dir(1)) .lt. 10*tiny(0._DP)) then
          s = (ys - startp(2)) / (endp(2) - startp(2))  !  line is vertical
        else
          s = (xs - startp(1)) / (endp(1) - startp(1))  !  line is horizontal
        end if
      else
!       s = r    arc is present --> polar coordinates necessary  [not yet implemented]
      end if
!
!----------------------------------------------------------------------------------------
!  Common settings for optics boundary conditions
!  ----------------------------------------------
!  Set width of waveguide
      w = userparam(1)

!  Set effective index and refractive indices of waveguide, material in positive and
!  negative direction respectively
      neff = userparam(2)
      nn   = userparam(3)
      nwg  = userparam(4)
      np   = userparam(5)

      if (nwg .lt. neff) then
!$omp critical (single_print_stop)
        print *,"***** Effective index greater than waveguide index. Check USERPARAM settings."
        stop
!$omp end critical (single_print_stop)
      end if

!  Compute vacuum wavenumber, propagation and attenuation constants and transverse
!  wavenumber
      k0 = omega*sqrt(mu0*eps0)
      beta = k0*neff
      gammap = sqrt(beta**2-(k0*np)**2)
      gamman = sqrt(beta**2-(k0*nn)**2)
      kappa = sqrt((k0*nwg)**2-beta**2)
!----------------------------------------------------------------------------------------
!
!
!  Select function for given number in drawing layer. The functions are given in depen- 
!  dence of the local coordinate s along to the branch. If there are any new functions,
!  users please add them as a new case.
!
      select case (bcnum)
      case (1)
!  half sine boundary function
        pval = sin(pi*s)

      case (2)
!  sine boundary function
        pval = sin(2*pi*s)

      case (3)
!  3/2 sine boundary function
        pval = sin(3*pi*s)

      case (4)
!  Gauss error distribution curve with Amplitude
        pval = exp(-(8._DP * s - 4._DP)**2)

      case (5)
!  exponential function
        pval = exp(s)

      case (89)    ! Marc
!  sine boundary function for a Heat Transfer test-problem
        pval = sin(pi*s/2.0_DP)

      case (90)   ! Bodin
!  sine boundary function for a Heat Transfer test-problem
        pval = sin(pi*xs/2.0_DP)

      case (91)   ! Pouseuille flow b.c.
!  parabolic velocity profile 
        pval = 0.01_DP*ys*(20.0_DP-ys)

      case (92)   ! Backward-facing step
!  parabolic velocity profile
        pval = 24.0_DP*ys*(0.5_DP-ys)

      case (201)  ! first order ABC for air
        pval = 0._DPC
        qval = cmplx(0._DP,-k0,DPC)

!----------------------------------------------------------------------------------------
!  definition of excitation in dependence of x-axis
      case (202)  ! evanescent field with xs < -w
        select case (physics)
        case ('TEWAVE')
          pval = cmplx(0._DP,beta*(cos(kappa*w)+gammap/kappa*sin(kappa*w))*exp(gamman*(xs+w)),DPC)
          qval = cmplx(0._DP,-beta,DPC)
        case ('TMWAVE')
          pval = cmplx(0._DP,beta/nn**2*(cos(kappa*w)+gammap/kappa*nwg**2/np**2*sin(kappa*w))*exp(gamman*(xs+w)),DPC)
          qval = cmplx(0._DP,-beta/nn**2,DPC)
        case default
          pval = 0._DPC
          qval = 0._DPC
        end select

      case (203)  ! excitation of mode with -w < xs < 0
        select case (physics)
        case ('TEWAVE')
          pval = cmplx(0._DP,beta*(cos(kappa*xs)-gammap/kappa*sin(kappa*xs)),DPC)
          qval = cmplx(0._DP,-beta,DPC)
        case ('TMWAVE')
          pval = cmplx(0._DP,beta/nwg**2*(cos(kappa*xs)-gammap/kappa*nwg**2/np**2*sin(kappa*xs)),DPC)
          qval = cmplx(0._DP,-beta/nwg**2,DPC)
        case default
          pval = 0._DPC
          qval = 0._DPC
        end select

      case (204)  ! evanescent field with xs > 0
        select case (physics)
        case ('TEWAVE')
          pval = cmplx(0._DP,beta*exp(-gammap*xs),DPC)
          qval = cmplx(0._DP,-beta,DPC)
        case ('TMWAVE')
          pval = cmplx(0._DP,beta/np**2*exp(-gammap*xs),DPC)
          qval = cmplx(0._DP,-beta/np**2,DPC)
        case default
          pval = 0._DPC
          qval = 0._DPC
        end select
!----------------------------------------------------------------------------------------
!  definition of excitation in dependence of y-axis
      case (206)  ! evanescent field with ys < -w
        select case (physics)
        case ('TEWAVE')
          pval = cmplx(0._DP,beta*(cos(kappa*w)+gammap/kappa*sin(kappa*w))*exp(gamman*(ys+w)),DPC)
          qval = cmplx(0._DP,-beta,DPC)
        case ('TMWAVE')
          pval = cmplx(0._DP,beta/nn**2*(cos(kappa*w)+gammap/kappa*nwg**2/np**2*sin(kappa*w))*exp(gamman*(ys+w)),DPC)
          qval = cmplx(0._DP,-beta/nn**2,DPC)
        case default
          pval = 0._DPC
          qval = 0._DPC
        end select

      case (207)  ! excitation of mode with -w < ys < 0
        select case (physics)
        case ('TEWAVE')
          pval = cmplx(0._DP,beta*(cos(kappa*ys)-gammap/kappa*sin(kappa*ys)),DPC)
          qval = cmplx(0._DP,-beta,DPC)
        case ('TMWAVE')
          pval = cmplx(0._DP,beta/nwg**2*(cos(kappa*ys)-gammap/kappa*nwg**2/np**2*sin(kappa*ys)),DPC)
          qval = cmplx(0._DP,-beta/nwg**2,DPC)
        case default
          pval = 0._DPC
          qval = 0._DPC
        end select

      case (208)  ! evanescent field with ys > 0
        select case (physics)
        case ('TEWAVE')
          pval = cmplx(0._DP,beta*exp(-gammap*ys),DPC)
          qval = cmplx(0._DP,-beta,DPC)
        case ('TMWAVE')
          pval = cmplx(0._DP,beta/np**2*exp(-gammap*ys),DPC)
          qval = cmplx(0._DP,-beta/np**2,DPC)
        case default
          pval = 0._DPC
          qval = 0._DPC
        end select
!----------------------------------------------------------------------------------------

      case (210)  ! only for test purpose
        pval = sin(pi*s)
        qval = 1._DP

!----------------------------------------------------------------------------------------

      case (250)  ! Convective flux boundary condition for mass concentration field
        select case (physics)
        case ('FLUIDELMASSCON')
          call xy2lam(xs,ys,elem,lambda,xn,yn,e)
          call fieldu(elem,lambda,u)
          pval = 0._DPC
          qval = 0._DPC
          qval(5) = u(1)
        case default
          pval = 0._DPC
          qval = 0._DPC
        end select

!----------------------------------------------------------------------------------------

      case default
!  if user layer exists, but no number is given, set value to 0
        pval = 0._DPC
      end select
!
      return
      end subroutine userbc
