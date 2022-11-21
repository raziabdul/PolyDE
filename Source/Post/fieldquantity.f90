      subroutine fieldquantity_scalar(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use feminterface, only: field, lam2xy, pdecoeff, invertmat, fetchmatparameters, principlevalues
      use femtypes
      use globalvariables, only: xn, yn, e, physics, pi, eps0, matzif, geb, mu0, omega, xfluid, nnat
      use matconstants, only: numparam
      implicit none
      real (DP) :: lambda(3), phi
      complex (DPC) :: zs
      character (len=*) :: fieldtype, unit
      character (len=*) :: descriptor
      integer (I4B) :: elem
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
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
!    $Revision: 1.31 $
!    $Date: 2015/06/02 10:39:48 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
!
!  The routine gives the required field quantity with descriptor and unit
!  depending on the physics mode used for a point given in barycentric coordinates.
!  To obtain the time value of the field quantity (for plotting), the real part
!  of zs has got to be taken.
!
!------------------------------------------------------------------------------
!  Input:
!            elem        number of the element
!            fieldtype   field quantity to obtain
!            lambda      barycentric coordinates for point
!            phi         phase angle (in degrees) for harmonic quantities
!  Output:
!            zs          value of field quantity at lambda in elem
!            descriptor  description of field quantity (e.g. electric field strength)
!            unit        unit of field quantity
!            ok          =.false. if fieldtype is not available for the physics mode
!
!  local variables
      real (DP) :: xs, ys, sigma1, sigma2, xys(2)
      real(DP), allocatable :: list(:)
      integer (I4B) :: chkdim, matindex, nature
      complex (DPC) :: comp1, comp2
      complex (DPC) :: zval(15,nnat), cfak
      complex (DPC) :: nu(2,2,nnat,nnat), gamma(2,nnat,nnat), beta(2,nnat,nnat)
      complex (DPC) :: alpha(nnat,nnat), em(nnat,nnat), f1(nnat), f2(2,nnat)
      complex (DPC) :: temp(2)
      complex (DPC) :: inverted(2,2), source(2), tempten(2,2)
      logical :: typ(5)
!------------------------------------------------------------------------------
!
!  initializations
      ok=.true.
      descriptor = ''
      unit = ''
!
!  get world coordinates
      call lam2xy(lambda,elem,xs,ys,xn,yn,e)
      xys=(/xs,ys/)
!  get pde coefficients
      call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2,em)
      select case (physics)
!
!------------------------------------------------------------------------------
!
      case ('ACOUSTIC')
        nature=1
        select case (fieldtype)
        case ('POTENTIAL')
          descriptor = 'velocity potential'
          unit = 'm2/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('VX')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(4,nature)*cfak
          typ = .false.
        case ('VY')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(5,nature)*cfak
          typ = .false.
        case ('V')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = sqrt((zval(4,nature)*cfak)**2 + (zval(5,nature)*cfak)**2)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('ELECTROSTATICS')
        nature=1
        select case (fieldtype)
        case ('POTENTIAL')
          descriptor = 'electric'
          unit = 'V'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs=zval(1,nature)
          typ = .false.
        case ('EX')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(4,nature)
          typ = .false.
        case ('EY')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(5,nature)
          typ = .false.
        case ('E')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = sqrt(zval(4,nature)**2 + zval(5,nature)**2)
          typ = .false.
        case ('DX')
          descriptor = 'electric flux density'
          unit = 'As/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(7,nature)*eps0
          typ = .false.
        case ('DY')
          descriptor = 'electric flux density'
          unit = 'As/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(8,nature)*eps0
          typ = .false.
        case ('D')
          descriptor = 'electric flux density'
          unit = 'As/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = sqrt((zval(7,nature)*eps0)**2 + (zval(8,nature)*eps0)**2)
          typ = .false.
        case ('WE')
          descriptor = 'electric energy density'
          unit = 'J/m3'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = 0.5_DP*(zval(4,nature)*zval(7,nature)*eps0 + zval(5,nature)*zval(8,nature)*eps0)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLOW_PROFILE')
        nature=1
        select case (fieldtype)
        case ('VZ')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDINCOMPR')
        select case (fieldtype)
        case ('UX')
          nature=1
          descriptor = 'velocity in x-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('UY')
          nature=2
          descriptor = 'velocity in y-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('P')
          nature=3
          descriptor = 'pressure'
          unit = 'N/m2'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDELINCOMPR')
        select case (fieldtype)
        case ('UX')
          nature=1
          descriptor = 'velocity in x-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('UY')
          nature=2
          descriptor = 'velocity in y-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('P')
          nature=3
          descriptor = 'pressure'
          unit = 'N/m2'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('POTENTIAL')
          nature=4
          descriptor = 'electric'
          unit = 'V'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs=zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDELMASSCON')
        select case (fieldtype)
        case ('UX')
          nature=1
          descriptor = 'velocity in x-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('UY')
          nature=2
          descriptor = 'velocity in y-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('P')
          nature=3
          descriptor = 'pressure'
          unit = 'N/m2'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('POTENTIAL')
          nature=4
          descriptor = 'electric'
          unit = 'V'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('MASSCON')
          nature=5
          descriptor = 'mass concentration'
          unit = 'mol/(m3*s)'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDELMASSTEP')
        select case (fieldtype)
        case ('UX')
          nature=1
          descriptor = 'velocity in x-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('UY')
          nature=2
          descriptor = 'velocity in y-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('P')
          nature=3
          descriptor = 'pressure'
          unit = 'N/m2'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('POTENTIAL')
          nature=4
          descriptor = 'electric'
          unit = 'V'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs=zval(1,nature)*cfak
          typ = .false.
        case ('MASSCON')
          nature=5
          descriptor = 'mass concentration'
          unit = 'mol/(m3*s)'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs=zval(1,nature)
          typ = .false.
        case ('T')
          nature=6
          descriptor = 'temperature'
          unit = 'K'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs=zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUID')
!  determine material index for given element in region geb
        matindex=matzif(geb(elem))
!  read parameters from internal list and assign them to local variables
        allocate (list(numparam))
        call fetchmatparameters(list,matindex,xys)
        nature=1
        select case (fieldtype)
        case ('UX')
          descriptor = 'velocity in x-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(2,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('UY')
          descriptor = 'velocity in y-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(3,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('USUM')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval,xfluid(2,:))
          comp1 = (zval(1,nature)*cfak)**2
          call field(elem,lambda,typ,zval,xfluid(3,:))
          comp2 = (zval(1,nature)*cfak)**2
          zs = sqrt(comp1+comp2)
          typ = .false.
        case ('P')
          descriptor = 'pressure'
          unit = 'N/m2'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(1,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('T')
          chkdim = size(xfluid,dim=1)
          if (chkdim .le. 3) then
            write(*,*) 'The solution of energy equation is unavailable'
            stop
          end if
          descriptor = 'temperature'
          unit = 'K'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(4,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('STREAMLINE')
          descriptor = 'streamline'
          unit = ''
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(5,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('STRESSXX')
        ! for only dimensional variables
          descriptor = 'stress in x-component'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval,xfluid(2,:))
          comp1 = zval(4,nature)*cfak
          call field(elem,lambda,typ,zval,xfluid(3,:))
          comp2 = zval(5,nature)*cfak
          zs = list(1)*(2.0_DP/3.0_DP)*(2.0_DP*comp1 - comp2)
          typ = .false.
        case ('STRESSYY')
        ! for only dimensional variables
          descriptor = 'stress in y-component'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval,xfluid(2,:))
          comp1 = zval(4,nature)*cfak
          call field(elem,lambda,typ,zval,xfluid(3,:))
          comp2 = zval(5,nature)*cfak
          zs = list(1)*(2.0_DP/3.0_DP)*(2.0_DP*comp2 - comp1)
          typ = .false.
        case ('SHEARSTRESS')
          descriptor = 'shear stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval,xfluid(2,:))
          ! eta*du/dy
          comp1 = zval(5,nature)*cfak
          call field(elem,lambda,typ,zval,xfluid(3,:))
          ! eta*dv/dx
          comp2 = zval(4,nature)*cfak
          zs = list(1)*(comp1 + comp2)
          typ = .false.
        case default
          ok=.false.
        end select
        deallocate(list)
!
!------------------------------------------------------------------------------
!
      case ('HEATTR')
        nature=1
        select case (fieldtype)
        case ('T')
          descriptor = 'temperature'
          unit = 'K'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('QX')
          descriptor = 'heat flux'
          unit = 'W/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = -zval(7,nature)*cfak
          typ = .false.
        case ('QY')
          descriptor = 'heat flux'
          unit = 'W/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = -zval(8,nature)*cfak
          typ = .false.
        case ('Q')
          descriptor = 'heat flux'
          unit = 'W/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = sqrt((zval(7,nature)*cfak)**2 + (zval(8,nature)*cfak)**2)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('TRANHEATTR')
        nature=1
        select case (fieldtype)
        case ('T')
          descriptor = 'temperature'
          unit = 'K'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          print *,"***** fieldtype ",fieldtype," unknown for ",physics
          stop
        end select
!
!------------------------------------------------------------------------------
!
      case ('LAPLACE')
        nature=1
        select case (fieldtype)
        case ('POTENTIAL')
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = zval(1,nature)
          typ = .false.
        case ('GRADUX')
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(4,nature)
          typ = .false.
        case ('GRADUY')
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(5,nature)
          typ = .false.
        case ('GRADU')
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = sqrt(zval(4,nature)**2 + zval(5,nature)**2)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('MAGNETICVP')
        nature=1
!  get pde coefficients at the integration points
        call lam2xy(lambda,elem,xs,ys,xn,yn,e)
        call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2)
        inverted(1,1) =  nu(2,2,1,1)
        inverted(1,2) = -nu(2,1,1,1)
        inverted(2,1) = -nu(1,2,1,1)
        inverted(2,2) =  nu(1,1,1,1)
        select case (fieldtype)
        case ('AZ')
          descriptor = 'magnetic vector potential'
          unit = 'Vs/m'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('BX')
          descriptor = 'magnetic flux density'
          unit = 'T'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          zs = zval(5,nature)*cfak
          typ = .false.
        case ('BY')
          descriptor = 'magnetic flux density'
          unit = 'T'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          zs = -zval(4,nature)*cfak
          typ = .false.
        case ('B')
          descriptor = 'magnetic flux density'
          unit = 'T'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          zs = sqrt((zval(5,nature)*cfak)**2 + (-zval(4,nature)*cfak)**2)
          typ = .false.
        case ('HX')
          descriptor = 'magnetic field strength'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
          zs = temp(1)*cfak
          typ = .false.
        case ('HY')
          descriptor = 'magnetic field strength'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
          zs = temp(2)*cfak
          typ = .false.
        case ('H')
          descriptor = 'magnetic field strength'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
          zs = sqrt((temp(1)*cfak)**2 + (temp(2)*cfak)**2)
          typ = .false.
        case ('WM')
          descriptor = 'magnetic energy density'
          unit = 'J/m3'

          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
! B
          temp = (/zval(5,nature),-zval(4,nature)/)
! H*conjg(B)
          temp = matmul(inverted,temp)*conjg(temp)
          zs = 0.25_DP*sum(real(temp,DP))

        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('STOKES')
        nature=1
        select case (fieldtype)
        case ('UX')
          descriptor = 'velocity in x-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(1,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('UY')
          descriptor = 'velocity in y-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(2,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('USUM')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval,xfluid(1,:))
          comp1 = (zval(1,nature)*cfak)**2
          call field(elem,lambda,typ,zval,xfluid(2,:))
          comp2 = (zval(1,nature)*cfak)**2
          zs = sqrt(comp1+comp2)
          typ = .false.
        case ('P')
          descriptor = 'pressure'
          unit = 'N/m2'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(3,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('T')
          descriptor = 'temperature'
          unit = 'K'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(4,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('STATCURRENT')
        nature=1
        select case (fieldtype)
        case ('POTENTIAL')
          descriptor = 'electric'
          unit = 'V'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs=zval(1,nature)
          typ = .false.
        case ('EX')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(4,nature)
          typ = .false.
        case ('EY')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(5,nature)
          typ = .false.
        case ('E')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = sqrt(zval(4,nature)**2 + zval(5,nature)**2)
          typ = .false.
        case ('JX')
          descriptor = 'electric current density'
          unit = 'A/m2'
!  determine material index for given element in region geb
          matindex=matzif(geb(elem))
!  read parameters from internal list and assign them to local variables
          allocate (list(numparam))
          call fetchmatparameters(list,matindex,xys)
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(7,nature)/list(4)
          deallocate(list)
          typ = .false.
        case ('JY')
          descriptor = 'electric current density'
          unit = 'A/m2'
!  determine material index for given element in region geb
          matindex=matzif(geb(elem))
!  read parameters from internal list and assign them to local variables
          allocate (list(numparam))
          call fetchmatparameters(list,matindex,xys)
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(8,nature)/list(4)
          deallocate(list)
          typ = .false.
        case ('J')
          descriptor = 'electric current density'
          unit = 'A/m2'
!  determine material index for given element in region geb
          matindex=matzif(geb(elem))
!  read parameters from internal list and assign them to local variables
          allocate (list(numparam))
          call fetchmatparameters(list,matindex,xys)
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = sqrt((zval(7,nature)/list(4))**2 + (zval(8,nature)/list(4))**2)
          deallocate(list)
          typ = .false.
        case ('IX')
          descriptor = 'electric current per meter'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(7,nature)
          typ = .false.
        case ('IY')
          descriptor = 'electric current per meter'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = -zval(8,nature)
          typ = .false.
        case ('I')
          descriptor = 'electric current per meter'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = sqrt(zval(7,nature)**2 + zval(8,nature)**2)
          typ = .false.
        case ('P')
          descriptor = 'electric power density'
          unit = 'W/m3'
          typ = (/.false.,.true.,.false.,.false.,.false./)
!  determine material index for given element in region geb
          matindex=matzif(geb(elem))
!  read parameters from internal list and assign them to local variables
          allocate (list(numparam))
          call fetchmatparameters(list,matindex,xys)
          call field(elem,lambda,typ,zval)
          zs = zval(4,nature)*zval(7,nature)/list(4) + zval(5,nature)*zval(8,nature)/list(4)
          deallocate(list)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('TEWAVE')
        nature=1
        inverted(1,1) =  nu(2,2,1,1)
        inverted(1,2) = -nu(2,1,1,1)
        inverted(2,1) = -nu(1,2,1,1)
        inverted(2,2) =  nu(1,1,1,1)
        inverted = inverted/mu0
        select case (fieldtype)
        case ('EZ')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('DZ')
          descriptor = 'electric flux density'
          unit = 'As/m2'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = -zval(2,nature)/(mu0*omega**2)*cfak
          typ = .false.
        case ('WE')
          descriptor = 'electric energy density'
          unit = 'J/m3'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = 0.25_DP*zval(1,nature)*conjg(-zval(2,nature)/(mu0*omega**2))
          typ = .false.
        case ('HX')
          descriptor = 'magnetic field strength'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
            zs = temp(1)*cfak/cmplx(0._DP,-omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('HY')
          descriptor = 'magnetic field strength'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
            zs = temp(2)*cfak/cmplx(0._DP,-omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('H')
          descriptor = 'magnetic field strength'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
            zs = sqrt((temp(1)*cfak)**2 + (temp(2)*cfak)**2)/omega
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('BX')
          descriptor = 'magnetic flux density'
          unit = 'T'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          zs = zval(5,nature)*cfak/cmplx(0._DP,-omega,DPC)
          typ = .false.
        case ('BY')
          descriptor = 'magnetic flux density'
          unit = 'T'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          zs = -zval(4,nature)*cfak/cmplx(0._DP,-omega,DPC)
          typ = .false.
        case ('B')
          descriptor = 'magnetic flux density'
          unit = 'T'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          zs = sqrt((zval(5,nature)*cfak)**2+(zval(4,nature)*cfak)**2)/omega
          typ = .false.
        case ('WM')
          descriptor = 'magnetic energy density'
          unit = 'J/m3'
        case ('SX')
          descriptor = 'Poynting vector'
          unit = 'VA/m2'
          typ = (/.true.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
            zs = -0.5_DP*zval(1,nature)*conjg(temp(2)/cmplx(0._DP,-omega,DPC))
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('SY')
          descriptor = 'Poynting vector'
          unit = 'VA/m2'
          typ = (/.true.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
            zs = 0.5_DP*zval(1,nature)*conjg(temp(1)/cmplx(0._DP,-omega,DPC))
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('S')
          descriptor = 'Poynting vector'
          unit = 'VA/m2'
          typ = (/.true.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
            zs = sqrt(abs(-0.5_DP*(zval(1,nature)*conjg(temp(2)/cmplx(0._DP,-omega,DPC))))**2 + &
                      abs(0.5_DP*(zval(1,nature)*conjg(temp(1)/cmplx(0._DP,-omega,DPC))))**2)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('TMWAVE')
        nature=1
        inverted(1,1) =  nu(2,2,1,1)
        inverted(1,2) = -nu(2,1,1,1)
        inverted(2,1) = -nu(1,2,1,1)
        inverted(2,2) =  nu(1,1,1,1)
        inverted = inverted/eps0
        select case (fieldtype)
        case ('HZ')
          descriptor = 'magnetic field strength'
          unit = 'A/m'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('BZ')
          descriptor = 'magnetic flux density'
          unit = 'T'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = -zval(2,nature)/(eps0*omega**2)*cfak
          typ = .false.
        case ('WM')
          descriptor = 'magnetic energy density'
          unit = 'J/m3'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs = 0.25_DP*conjg(zval(1,nature))*(-zval(2,nature)/(eps0*omega**2))
          typ = .false.
        case ('EX')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
!  get pde coefficients at the integration points
            call lam2xy(lambda,elem,xs,ys,xn,yn,e)
            call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2)
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/)-source)
            zs = temp(1)*cfak/cmplx(0._DP,omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('EY')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/)-source)
            zs = temp(2)*cfak/cmplx(0._DP,omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('E')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/)-source)
            zs = sqrt((temp(1)*cfak/cmplx(0._DP,omega,DPC))**2 + (temp(2)*cfak/cmplx(0._DP,omega,DPC))**2)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('DX')
          descriptor = 'electric flux density'
          unit = 'As/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = (/zval(5,nature),-zval(4,nature)/)-source
            zs = temp(1)*cfak/cmplx(0._DP,omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('DY')
          descriptor = 'electric flux density'
          unit = 'As/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = (/zval(5,nature),-zval(4,nature)/)-source
            zs = temp(2)*cfak/cmplx(0._DP,omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('D')
          descriptor = 'electric flux density'
          unit = 'As/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = (/zval(5,nature),-zval(4,nature)/)-source
            zs = sqrt((temp(1)*cfak/cmplx(0._DP,omega,DPC))**2 + (temp(2)*cfak/cmplx(0._DP,omega,DPC))**2)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('WE')
          descriptor = 'electric energy density'
          unit = 'J/m3'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
! D
            temp = ((/zval(5,nature),-zval(4,nature)/)-source)/(-omega)
! E*conjg(D)
            temp = matmul(inverted,temp)*conjg(temp)
            zs = 0.25_DP*sum(real(temp,DP))
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('SX')
          descriptor = 'Poynting vector'
          unit = 'VA/m2'
          typ = (/.true.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/)-source)
            zs = 0.5_DP*conjg(zval(1,nature))*temp(2)/cmplx(0._DP,omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('SY')
          descriptor = 'Poynting vector'
          unit = 'VA/m2'
          typ = (/.true.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/)-source)
            zs = -0.5_DP*conjg(zval(1,nature))*temp(1)/cmplx(0._DP,omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('S')
          descriptor = 'Poynting vector'
          unit = 'VA/m2'
          typ = (/.true.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/)-source)
            zs=sqrt((0.5_DP*(conjg(zval(1,nature))*(temp(2)/cmplx(0._DP,omega,DPC))))**2 + &
                    (-0.5_DP*(conjg(zval(1,nature))*(temp(1)/cmplx(0._DP,omega,DPC))))**2)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('THERMOELECTRIC')
        select case (fieldtype)
        case ('T')
          nature=1
          descriptor = 'temperature'
          unit = 'K'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('POTENTIAL')
          nature=2
          descriptor = 'electric'
          unit = 'V'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('QX')
          nature=1
          descriptor = 'heat flux'
          unit = 'W/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = -zval(7,nature)*cfak
          typ = .false.
        case ('QY')
          nature=1
          descriptor = 'heat flux'
          unit = 'W/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = -zval(8,nature)*cfak
          typ = .false.
        case ('JX')
          nature=2
          descriptor = 'electric current density'
          unit = 'A/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = -zval(7,nature)*cfak
          typ = .false.
        case ('JY')
          nature=2
          descriptor = 'electric current density'
          unit = 'A/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = -zval(8,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PLANE STRAIN')
        select case (fieldtype)
        case ('DEL-X')
          nature=1
          descriptor = 'x-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('DEL-Y')
          nature=2
          descriptor = 'y-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PLANE STRESS')
        select case (fieldtype)
        case ('DEL-X')
          nature=1
          descriptor = 'x-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('DEL-Y')
          nature=2
          descriptor = 'y-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('VAN MISES STRESS','VON MISES STRESS')
          descriptor = 'von Mises stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call principlevalues(real(zval(7,1)*cfak), real(zval(8,2)*cfak),  &
     &             0.5_DP*real((zval(7,2)+zval(8,1))*cfak),sigma1,sigma2)
          zs = sqrt(sigma1**2 + sigma2**2 - sigma1*sigma2)
          typ = .false.
        case ('TRESCA STRESS')
          descriptor = 'Trestca stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call principlevalues(real(zval(7,1)*cfak), real(zval(8,2)*cfak),  &
     &             0.5_DP*real((zval(7,2)+zval(8,1))*cfak),sigma1,sigma2)
          zs = max(abs(sigma1),abs(sigma2),abs(sigma1-sigma2))
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('THERMOELASTIC')
        select case (fieldtype)
        case ('DEL-X')
          nature=1
          descriptor = 'x-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('DEL-Y')
          nature=2
          descriptor = 'y-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
         case ('T')
          nature=3
          descriptor = 'temperature'
          unit = 'K'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('GRAVITOELASTIC')
        select case (fieldtype)
        case ('DEL-X')
          nature=1
          descriptor = 'x-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('DEL-Y')
          nature=2
          descriptor = 'y-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
         case ('GPOTENTIAL')
          nature=3
          descriptor = 'gravity potential'
          unit = 'm2/s2'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PIEZOELECTRIC')
        select case (fieldtype)
        case ('DEL-X')
          nature=1
          descriptor = 'x-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('DEL-Y')
          nature=2
          descriptor = 'y-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
         case ('POTENTIAL')
          nature=3
          descriptor = 'electric'
          unit = 'V'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PIEZOPYROELEC')
        select case (fieldtype)
        case ('DEL-X')
          nature=1
          descriptor = 'x-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('DEL-Y')
          nature=2
          descriptor = 'y-displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('POTENTIAL')
          nature=3
          descriptor = 'electric'
          unit = 'V'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case ('T')
          nature=4
          descriptor = 'temperature'
          unit = 'K'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case default
        print *,"***** unknown physics mode:",physics
        stop
      end select
!
      end subroutine fieldquantity_scalar
      
      
      
      subroutine fieldquantity_vector(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use feminterface, only: field, lam2xy, pdecoeff, invertmat, fetchmatparameters
      use femtypes
      use globalvariables
      use matconstants
      implicit none
      real (DP) :: lambda(3), phi
      complex (DPC) :: zs(:)
      character (len=*) :: fieldtype, unit
      character (len=*) :: descriptor
      integer (I4B) :: elem
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs,  ok
!  Input:
!            elem        number of the element
!            fieldtype   field quantity to obtain
!            lambda      barycentric coordinates for point
!            phi         phase angle (in degrees) for harmonic quantities
!  Output:
!            zs          value of field quantity at lambda in elem
!            descriptor  description of field quantity (e.g. electric field strength)
!            unit        unit of field quantity
!            ok          =.false. if fieldtype is not available for the physics mode
!
!  local variables
      real (DP) :: xs, ys, xys(2)
      real(DP), allocatable :: list(:)
      integer (I4B) :: matindex, nature
      complex (DPC) :: zval(14,nnat), cfak
      complex (DPC) :: nu(2,2,nnat,nnat), gamma(2,nnat,nnat), alpha(nnat,nnat)
      complex (DPC) :: beta(2,nnat,nnat), f1(nnat), f2(2,nnat)
      complex (DPC) :: temp(2)
      complex (DPC) :: inverted(2,2), source(2), tempten(2,2)
      logical :: typ(5)
!------------------------------------------------------------------------------
!
!  initializations
      ok=.true.
      descriptor = ''
      unit = ''
!
!  get world coordinates
      call lam2xy(lambda,elem,xs,ys,xn,yn,e)
      xys=(/xs,ys/)

      select case (physics)
!
!------------------------------------------------------------------------------
!
      case ('ACOUSTIC')
        nature=1
        select case (fieldtype)
        case ('V')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = zval(4,nature)*cfak
          zs(2) = zval(5,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('ELECTROSTATICS')
        nature=1
        select case (fieldtype)
        case ('E')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs(1) = -zval(4,nature)
          zs(2) = -zval(5,nature)
          typ = .false.
        case ('D')
          descriptor = 'electric flux density'
          unit = 'As/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs(1) = -zval(7,nature)*eps0
          zs(2) = -zval(8,nature)*eps0
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLOW_PROFILE')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDINCOMPR')
        select case (fieldtype)
        case ('USUM')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=2
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDELINCOMPR')
        select case (fieldtype)
        case ('USUM')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=2
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDELMASSCON')
        select case (fieldtype)
        case ('USUM')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=2
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDELMASSTEP')
        select case (fieldtype)
        case ('USUM')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=2
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUID')
        nature=1
        select case (fieldtype)
        case ('U')
          descriptor = 'velocity in x-direction'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(2,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = zval(1,nature)*cfak
          call field(elem,lambda,typ,zval,xfluid(3,:))
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('HEATTR')
        nature=1
        select case (fieldtype)
        case ('QX')
        case ('Q')
          descriptor = 'heat flux'
          unit = 'W/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = -zval(7,nature)*cfak
          zs(2) = -zval(8,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('LAPLACE')
        nature=1
        select case (fieldtype)
        case ('GRADU')
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs(1) = -zval(4,nature)
          zs(2) = -zval(5,nature)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('MAGNETICVP')
        nature=1
!  get pde coefficients at the integration points
        call lam2xy(lambda,elem,xs,ys,xn,yn,e)
        call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2)
        inverted(1,1) =  nu(2,2,1,1)
        inverted(1,2) = -nu(2,1,1,1)
        inverted(2,1) = -nu(1,2,1,1)
        inverted(2,2) =  nu(1,1,1,1)
        select case (fieldtype)
        case ('B')
          descriptor = 'magnetic flux density'
          unit = 'T'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          zs(1) = zval(5,nature)*cfak
          zs(2) = -zval(4,nature)*cfak
          typ = .false.
        case ('H')
          descriptor = 'magnetic field strength'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
          zs(1) = temp(1)*cfak
          zs(2) = temp(2)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('STOKES')
        nature=1
        select case (fieldtype)
        case ('U')
          descriptor = 'velocity'
          unit = 'm/s'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval,xfluid(1,:))
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = zval(1,nature)*cfak
          call field(elem,lambda,typ,zval,xfluid(2,:))
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('STATCURRENT')
        nature=1
        select case (fieldtype)
        case ('E')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs(1) = -zval(4,nature)
          zs(2) = -zval(5,nature)
          typ = .false.
        case ('J')
          descriptor = 'electric current density'
          unit = 'A/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
!  determine material index for given element in region geb
          matindex=matzif(geb(elem))
!  read parameters from internal list and assign them to local variables
          allocate (list(numparam))
          call fetchmatparameters(list,matindex,xys)
          call field(elem,lambda,typ,zval)
          zs(1) = -zval(7,nature)/list(4)
          zs(2) = -zval(8,nature)/list(4)
          deallocate(list)
          typ = .false.
        case ('I')
          descriptor = 'electric current per meter'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          zs(1) = -zval(7,nature)
          zs(2) = -zval(8,nature)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('TEWAVE')
        nature=1
        inverted(1,1) =  nu(2,2,1,1)
        inverted(1,2) = -nu(2,1,1,1)
        inverted(2,1) = -nu(1,2,1,1)
        inverted(2,2) =  nu(1,1,1,1)
        inverted = inverted/mu0
        select case (fieldtype)
        case ('H')
          descriptor = 'magnetic field strength'
          unit = 'A/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
            zs(1) = temp(1)*cfak/cmplx(0._DP,-omega,DPC)
            zs(2) = temp(2)*cfak/cmplx(0._DP,-omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('B')
          descriptor = 'magnetic flux density'
          unit = 'T'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          zs(1) = zval(5,nature)*cfak/cmplx(0._DP,-omega,DPC)
          zs(1) = -zval(4,nature)*cfak/cmplx(0._DP,-omega,DPC)
          typ = .false.
        case ('S')
          descriptor = 'Poynting vector'
          unit = 'VA/m2'
          typ = (/.true.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/))
            zs(1) = -0.5_DP*zval(1,nature)*conjg(temp(2)/cmplx(0._DP,-omega,DPC))
            zs(2) = 0.5_DP*zval(1,nature)*conjg(temp(1)/cmplx(0._DP,-omega,DPC))
          else
            zs = 0._DPC
          end if
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('TMWAVE')
        nature=1
        inverted(1,1) =  nu(2,2,1,1)
        inverted(1,2) = -nu(2,1,1,1)
        inverted(2,1) = -nu(1,2,1,1)
        inverted(2,2) =  nu(1,1,1,1)
        inverted = inverted/eps0
        select case (fieldtype)
        case ('E')
          descriptor = 'electric field strength'
          unit = 'V/m'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
!  get pde coefficients at the integration points
            call lam2xy(lambda,elem,xs,ys,xn,yn,e)
            call pdecoeff(elem,xs,ys,nu,gamma,alpha,beta,f1,f2)
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/)-source)
            zs(1) = temp(1)*cfak/cmplx(0._DP,omega,DPC)
            zs(2) = temp(2)*cfak/cmplx(0._DP,omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('D')
          descriptor = 'electric flux density'
          unit = 'As/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = (/zval(5,nature),-zval(4,nature)/)-source
            zs(1) = temp(1)*cfak/cmplx(0._DP,omega,DPC)
            zs(2) = temp(2)*cfak/cmplx(0._DP,omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case ('S')
          descriptor = 'Poynting vector'
          unit = 'VA/m2'
          typ = (/.true.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          if (omega .gt. 0._DP) then
            call invertmat(nu(:,:,1,1),tempten)
            temp = matmul(tempten,f2(:,1))
            source = (/-temp(2),temp(1)/)
            temp = matmul(inverted,(/zval(5,nature),-zval(4,nature)/)-source)
            zs(1) = 0.5_DP*conjg(zval(1,nature))*temp(2)/cmplx(0._DP,omega,DPC)
            zs(2) = -0.5_DP*conjg(zval(1,nature))*temp(1)/cmplx(0._DP,omega,DPC)
          else
            zs = 0._DPC
          end if
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('THERMOELECTRIC')
        select case (fieldtype)
        case ('Q')
          nature=1
          descriptor = 'heat flux'
          unit = 'W/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = -zval(7,nature)*cfak
          zs(2) = -zval(8,nature)*cfak
          typ = .false.
        case ('J')
          nature=2
          descriptor = 'electric current density'
          unit = 'A/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = -zval(7,nature)*cfak
          zs(2) = -zval(8,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PLANE STRAIN')
        select case (fieldtype)
        case ('DISPLACEMENT')
          descriptor = 'displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=2
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PLANE STRESS')
        select case (fieldtype)
        case ('DISPLACEMENT')
          descriptor = 'displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=2
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('THERMOELASTIC')
        select case (fieldtype)
        case ('DISPLACEMENT')
          descriptor = 'displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=1
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('GRAVITOELASTIC')
        select case (fieldtype)
        case ('DISPLACEMENT')
          descriptor = 'displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=2
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PIEZOELECTRIC')
        select case (fieldtype)
        case ('DISPLACEMENT')
          descriptor = 'displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=2
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PIEZOPYROELEC')
        select case (fieldtype)
        case ('DISPLACEMENT')
          descriptor = 'displacement'
          unit = 'm'
          typ = (/.true.,.false.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs(1) = zval(1,nature)*cfak
          nature=2
          zs(2) = zval(1,nature)*cfak
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case default
        print *,"***** unknown physics mode:",physics
      end select
!
      end subroutine fieldquantity_vector



      subroutine fieldquantity_tensor(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use feminterface, only: field, lam2xy, pdecoeff, invertmat, fetchmatparameters
      use femtypes
      use globalvariables
      use matconstants
      implicit none
      real (DP) :: lambda(3), phi
      complex (DPC) :: zs(:,:)
      character (len=*) :: fieldtype, unit
      character (len=*) :: descriptor
      integer (I4B) :: elem
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
!  Input:
!            elem        number of the element
!            fieldtype   field quantity to obtain
!            lambda      barycentric coordinates for point
!            phi         phase angle (in degrees) for harmonic quantities
!  Output:
!            zs          value of field quantity at lambda in elem
!            descriptor  description of field quantity (e.g. electric field strength)
!            unit        unit of field quantity
!            ok          =.false. if fieldtype is not available for the physics mode
!
!  local variables
      real (DP) :: xs, ys, xys(2)
      real(DP), allocatable :: list(:)
      integer (I4B) :: matindex, nature
      complex (DPC) :: comp1, comp2
      complex (DPC) :: zval(14,nnat), cfak
      logical :: typ(5)
!------------------------------------------------------------------------------
!
!  initializations
      ok=.true.
      descriptor = ''
      unit = ''
!
!  get world coordinates
      call lam2xy(lambda,elem,xs,ys,xn,yn,e)
      xys=(/xs,ys/)

      select case (physics)
!
!------------------------------------------------------------------------------
!
      case ('ACOUSTIC')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('ELECTROSTATICS')
        nature=1
        select case (fieldtype)
        case ('MAXWELL STRESS')
          descriptor = 'electric stress tensor'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1,1) = 0.5_DP*eps0*(zval(7,nature)*zval(4,nature) - zval(8,nature)*zval(5,nature))
          zs(2,2) = -zs(1,1)
          zs(1,2) = eps0*zval(7,nature)*zval(5,nature)
          zs(2,1) = eps0*zval(8,nature)*zval(4,nature)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLOW_PROFILE')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDINCOMPR')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDELINCOMPR')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDELMASSCON')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUIDELMASSTEP')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('FLUID')
        nature=1
!  determine material index for given element in region geb
        matindex=matzif(geb(elem))
!  read parameters from internal list and assign them to local variables
        allocate (list(numparam))
        call fetchmatparameters(list,matindex,xys)
        select case (fieldtype)
        case ('STRESS')
!  for only dimensional variables
          descriptor = 'stress in x-component'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          call field(elem,lambda,typ,zval,xfluid(2,:))
          comp1 = zval(4,nature)*cfak
!  eta*du/dy
          call field(elem,lambda,typ,zval,xfluid(3,:))
          comp2 = zval(5,nature)*cfak
!  eta*dv/dx
          zs(1,1) = list(1)*(2.0_DP/3.0_DP)*(2.0_DP*comp1 - comp2)
          zs(2,2) = list(1)*(2.0_DP/3.0_DP)*(2.0_DP*comp2 - comp1)
          zs(1,2) = list(1)*(comp1 + comp2)
          zs(2,1) = list(1)*(comp1 + comp2)
          typ = .false.
        case default
          ok=.false.
        end select
        deallocate(list)
!
!------------------------------------------------------------------------------
!
      case ('HEATTR')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('LAPLACE')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('MAGNETICVP')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('STOKES')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('STATCURRENT')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('TEWAVE')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('TMWAVE')
        nature=1
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('THERMOELECTRIC')
        select case (fieldtype)
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PLANE STRAIN')
        select case (fieldtype)
        case ('STRESS')
          descriptor = 'stress tensor'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1,1) = zval(7,1)*cfak
          zs(2,2) = zval(8,2)*cfak
          zs(1,2) = 0.5_DP * (zval(7,2) + zval(8,1))*cfak
          zs(2,1) = zs(1,2)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PLANE STRESS')
        select case (fieldtype)
        case ('STRESS')
          descriptor = 'stress tensor'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1,1) = zval(7,1)*cfak
          zs(2,2) = zval(8,2)*cfak
          zs(1,2) = 0.5_DP * (zval(7,2) + zval(8,1))*cfak
          zs(2,1) = zs(1,2)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('THERMOELASTIC')
        select case (fieldtype)
        case ('STRESS')
          descriptor = 'stress tensor'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1,1) = zval(7,1)*cfak
          zs(2,2) = zval(8,2)*cfak
          zs(1,2) = 0.5_DP * (zval(7,2) + zval(8,1))*cfak
          zs(2,1) = zs(1,2)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('GRAVITOELASTIC')
        select case (fieldtype)
        case ('STRESS')
          descriptor = 'stress tensor'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1,1) = zval(7,1)*cfak
          zs(2,2) = zval(8,2)*cfak
          zs(1,2) = 0.5_DP * (zval(7,2) + zval(8,1))*cfak
          zs(2,1) = zs(1,2)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PIEZOELECTRIC')
        select case (fieldtype)
        case ('STRESS')
          descriptor = 'stress tensor'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1,1) = zval(7,1)*cfak
          zs(2,2) = zval(8,2)*cfak
          zs(1,2) = 0.5_DP * (zval(7,2) + zval(8,1))*cfak
          zs(2,1) = zs(1,2)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case ('PIEZOPYROELEC')
        select case (fieldtype)
        case ('STRESS')
          descriptor = 'stress tensor'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call field(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1,1) = zval(7,1)*cfak
          zs(2,2) = zval(8,2)*cfak
          zs(1,2) = 0.5_DP * (zval(7,2) + zval(8,1))*cfak
          zs(2,1) = zs(1,2)
          typ = .false.
        case default
          ok=.false.
        end select
!
!------------------------------------------------------------------------------
!
      case default
        print *,"***** unknown physics mode:",physics
      end select
!
      end subroutine fieldquantity_tensor
!
!
!
      subroutine principlevalues(vx,vy,vxy,lambda1,lambda2,angle)
      use femtypes
      implicit none
      real (DP) :: vx, vy, vxy, lambda1, lambda2
      real (DP), optional :: angle
      intent (in) :: vx, vy, vxy
      intent (out) :: lambda1, lambda2, angle
!----------------------------------------------------------------------
!  Compute principle values (eigenvalues) and first princple axis of a 2x2 tensor
!----------------------------------------------------------------------
!
!  Input:
!     vx, vy, vxy    tensor components in x and y direction and shear component vxy
!  Output:
!     lambda1        first principle value
!     lambda2        second princple value (abs(lambda2) < abs(lambda1))
!     angle          angle (with respect to the x-axis) of the first principle axis
!
!  local variables
      real (DP) :: b, c, q
!
!  get principle values and orientationn (angle) of first principle axis
      b = -vx-vy
      c = vx*vy-vxy**2
      if (b .gt. 0) then
        q = - 0.5_DP * (b + sqrt((vx-vy)**2+4._DP*vxy**2))
      else
        q = - 0.5_DP * (b - sqrt((vx-vy)**2+4._DP*vxy**2))
      end if
      lambda1= q
      lambda2= (vx*vy-vxy**2)/q
      if (present(angle)) then
        angle=atan2(lambda1-vx,vxy)
      end if
      return
      end 
