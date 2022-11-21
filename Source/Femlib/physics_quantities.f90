module physics_quantities
!  global variabled of the module
      use femtypes
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
!    $Revision: 1.1 $
!    $Date: 2014/07/03 12:57:18 $
!    $Author: m_kasper $
!
!  module variables
!    dimensions       no of spatial dimensions, i.e. '2' or '3'
!
!    eltypes          element type for each nature 
!                     can be   1 = Lagrange,   
!                              2 = Nédélec or   
!                             -1 = Lagrange with element order lowerered by 1 (e.g. pressure in Fluid flow)
!
!    pot_names         name of the potential quantity of the nature                                            (e.g. Phi for Electrostatics)
!
!    gx_names          names of the quantity which equal the negative Gradient of nature variable, x-component (e.g. Ex  for Electrostatics)
!    gy_names                                                                                    , y-component (e.g. Ey  for Electrostatics)
!    gz_names                                                                                    , z-component (e.g. Ez  for Electrostatics)
!    g_names                                                    Gradient magnitude ( sqrt(gx^2+gy^2) )         (e.g. E   for Electrostatics)
!
!    fx_names          names of the quantity which Flux of nature variable: -nu*grad(u)|x        , x-component (e.g. Dx  for Electrostatics)
!    fy_names                                                               -nu*grad(u)|y        , y-component (e.g. Dy  for Electrostatics)
!    fz_names                                                               -nu*grad(u)|z        , z-component (e.g. Dz  for Electrostatics)
!    f_names                                              Flux magnitude ( sqrt(fx^2+fy^2) )                   (e.g. D   for Electrostatics)
!
!    pot_descriptor    description of the potential quantity of the nature                    (e.g. "electric Potential" for Electrostatics)
!
!    gx_descriptor     description of the negative Gradient of nature variable, x-component (e.g. "electric field strength" for Electrostatics)
!    gy_descriptor                                                            , y-component
!    gz_descriptor                                                            , z-component
!    g_descriptor                                                               magnitude
!
!    fx_descriptor     description of the Flux quantity of nature variable,   x-component (e.g. "electric flux density"  for Electrostatics)
!    fy_descriptor                                                          , y-component
!    fz_descriptor                                                          , z-component
!    f_descriptor                                                             magnitude
!
!    pot_unit          unit of the potential quantity of the nature                                            (e.g. "V" for Electrostatics)
!
!    gx_unit           unit of the negative Gradient of nature variable, x-component                         (e.g. "V/m" for Electrostatics)
!    gy_unit                                                            , y-component
!    gz_unit                                                            , z-component
!    g_unit                                                               magnitude
!
!    fx_unit           unit of the Flux quantity of nature variable, x-component                         (e.g. "As/m^2"  for Electrostatics)
!    fy_unit                                                       , y-component
!    fz_unit                                                       , z-component
!    f_unit                                                          magnitude

!
      save
      integer(I4B) :: dimensions
      integer(I4B), allocatable :: eltypes(:)
      character (len=20), allocatable :: pot_names(:), gx_names(:), gy_names(:), gz_names(:),                        &
     &                     g_names(:), fx_names(:), fy_names(:), fz_names(:), f_names(:)
      character (len=50), allocatable :: pot_descriptor(:), gx_descriptor(:), gy_descriptor(:), gz_descriptor(:),    &
     &                     g_descriptor(:), fx_descriptor(:), fy_descriptor(:), fz_descriptor(:), f_descriptor(:)
      character (len=12), allocatable :: pot_unit(:), gx_unit(:), gy_unit(:), gz_unit(:),                            &
     &                     g_unit(:), fx_unit(:), fy_unit(:), fz_unit(:), f_unit(:)
end module physics_quantities
