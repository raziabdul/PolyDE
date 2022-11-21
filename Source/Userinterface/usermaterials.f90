      pure subroutine usermaterials(matname,found,xys,names,values,numnames,elem)
      use femtypes
      use feminterface, only: xy2lam, fieldu
      use globalvariables, only: xn, yn, e, nnat
      implicit none
      real (DP), optional :: xys(:), values(:)
      integer (I4B),optional :: elem
      integer (I4B), optional :: numnames
      character (len=*) :: matname
      character (len=*), optional :: names(:)
      logical :: found
      intent (in) :: elem, matname, xys
      intent (out) :: names, values, found, numnames
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
!    $Revision: 1.7 $
!    $Date: 2014/07/03 13:22:02 $
!    $Author: m_kasper $
!
!
!  User programming interface for material parameters which cannot be obtained
!  by use of material files (inhomogenous, spatial or parameter dependent materials)
!
!  Input:
!        matname     the actual material, for which parameters should be returned
!        xs, ys      coordinates for which the material values are to be returned
!        elem        element number
!  Output:
!        found       = .true. if the material is defined in this file
!        names       list of names of material parameter to be given
!        values      list of values for the corresponding materialparameters
!        numnames    number of name/ values pairs beeing returned
!
!  Notes:
!  -  in the case of a inhomogenous (spatial dependency) materials use matvar=.true.
!  -  material names (matname) of usermaterials are case sensitive
!     (unlike filenames of materialfiles in Windows)
!  -  the possible material parameters are listed in the documetation (PDEcoeff)
!     and defined in the subroutine setstandardvalues
!  -  the list of possible material parameters depends on the physics mode
!  -  parameter names are not case sensitve (use upper case)
!  -  if a parameter is not defined below, the standard value is used
!  -  you should avoid to give more than the required parameters
!
!
!  local variables:
      real (DP) lambda(3), xs, ys
      complex (DPC) u(nnat)
!
      select case (matname)
!######################################  add materials here below #####################################

      case ('my-material')
        found=.true.
        if (.not.present(names)) return
!  permeability
        names(1)='velx'
        values(1)=-sqrt(xys(1)**2+xys(2)**2)*sin(atan2(xys(2),xys(1)))
        names(2)='vely'
        values(2)=sqrt(xys(1)**2+xys(2)**2)*cos(atan2(xys(2),xys(1)))
        names(3)='kappa'
        values(3)=25.e6
        numnames=3
!
      case ('Si-non-lin')
        found=.true.
        if (.not.present(names)) return
!
!  Get the current temperature from the field routine
        xs = xys(1)
        ys = xys(2)
        call xy2lam(xs,ys,elem,lambda,xn,yn,e)
        call fieldu(elem,lambda,u)
!
!  thermal conductivity in the first principle axis
        names(1)='lam11'
        values(1)=( (163.) * abs(u(1)) )+100.
!  thermal conductivity in the second principle axis
        names(2)='lam22'
        values(2)=( (163.) * abs(u(1)) )+100.
!  thickness of the plate in m
        names(3)='d'
        values(3)=1.
!  material density in kg/m3
        names(4)='rho'
        values(4)=2330.
!  specific heat capacitance
        names(5)='c_p'
        values(5)=703.
!  power density
        names(6)='p'
        values(6)=0.
!  phase of power density
        names(7)='phase'
        values(7)=0.
!
        numnames=7
!#####################################################################################################
      case default
!  if the material is not define in the file
        found=.false.
        return
      end select
!
      return
      end subroutine usermaterials