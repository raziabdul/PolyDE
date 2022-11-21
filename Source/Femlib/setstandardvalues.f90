      subroutine setstandardvalues
      use feminterface, only: low2hi
      use femtypes
      use matconstants, only: numparam, parameternames, maxmat, param, parameterdefault
      use globalvariables, only: physics
      implicit none
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
!    $Revision: 1.39 $
!    $Date: 2014/07/01 15:41:45 $
!    $Author: m_kasper $
!
!  set the default values for the material parameters 
!
      integer(I4B) i, j
!
      select case(physics)
!
!
      case ('FLUID')
        numparam=5
!
        allocate(parameternames(numparam),parameterdefault(numparam))
!  dynamic viscosity in kg/m.s (mu or eta)
        parameternames(1)='eta'
        parameterdefault(1)= 1.73e-5_DP   ! air at 298 K, 101325 Pa
!  material density in kg/m3
        parameternames(2)='rho'
        parameterdefault(2)= 1.229_DP    ! water at 298 K, 101325 Pa
!  thermal conductivity
        parameternames(3)='lam11'
        parameterdefault(3)= 0.0267_DP
!  specific heat capacitance
        parameternames(4)='c_p'
        parameterdefault(4)= 1000._DP
!  coefficient of thermal expansion of fluid
        parameternames(5)='thexpn'
        parameterdefault(5)= 0._DP
!
!
      case ('STOKES')
        numparam=9
!
        allocate(parameternames(numparam),parameterdefault(numparam))
!  dynamic viscosity (mu or eta)
        parameternames(1)='eta'
        parameterdefault(1)= 1.73e-5_DP       ! air at 15 °C
!  material density in kg/m3
        parameternames(2)='rho'
        parameterdefault(2)= 1.229_DP         ! air at 15 °C
!  constant source term
        parameternames(3)='source'
        parameterdefault(3)= 0._DP
!  external forces x direction
        parameternames(4)='f_x'
        parameterdefault(4)=0._DP
!  external forces y direction
        parameternames(5)='f_y'
        parameterdefault(5)=0._DP
!  thermal conductivity in the first principle axis
        parameternames(6)='lam11'
        parameterdefault(6)=0.0267_DP
!  thermal conductivity in the second principle axis
        parameternames(7)='lam22'
        parameterdefault(7)=0.0267_DP
!  angle between the x-axis and the first principle axis
        parameternames(8)='angle'
        parameterdefault(8)=0._DP
!  specific heat capacitance
        parameternames(9)='c_p'
        parameterdefault(9)=1000._DP
!
!
      end select
!
!  convert parameter names to UPPERCASE
      do i=1,numparam
        call low2hi(parameternames(i),len_trim(parameternames(i)))
      end do
!  allocate the array of parameters and copy the standard values to it
      do i=0,maxmat
        allocate(param(i)%d(numparam))
        do j=1,numparam
          param(i)%d(j)=parameterdefault(j)
        end do
      end do
      deallocate(parameterdefault)
      return
      end subroutine setstandardvalues