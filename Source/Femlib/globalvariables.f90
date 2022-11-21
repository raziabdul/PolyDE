      module globalvariables
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
!    $Revision: 1.19 $
!    $Date: 2014/02/11 16:23:53 $
!    $Author: juryanatzki $
!
!  Geometry
!  Key-points
      integer (I4B) gkz
      integer (I4B), pointer :: kzrb(:,:)
      real (DP), pointer :: xbk(:), ybk(:)
!  Branches
      integer (I4B) gzz
      integer (I4B), pointer :: zrb(:,:), zki(:,:)
      complex(DPC), pointer :: alrb(:,:), btrb(:,:)
!  Regions
      integer (I4B) gbz
      integer (I4B), pointer :: matzif(:)
!  Auxiliary Information
      real (DP) xmin, xmax, ymin, ymax
!  Mesh Information
!  Elements
      integer (I4B) n
      integer (I4B), pointer :: e(:,:), en(:,:), geb(:)
!  Nodes
      integer (I4B) p
      integer (I4B), pointer :: kzi(:)
      real (DP), pointer :: xn(:),yn(:)
!  Solution
!  actual estimated global error of the FEM solution
      real (DP) :: fem_accuracy
!  Degrees of freedom
      integer (I4B) ndof
!  Polynomial degree of a 2D element, for every nature
      integer (I4B), pointer :: ep(:,:)
      complex (DPC), pointer :: x(:), dxdt(:), d2xdt2(:)
      complex (DPC), pointer :: xfluid(:,:)
      type (ARRPTRI), pointer :: eg(:,:)
!  Number of natures
!  nnat stands for the number of natures of multiphysics problems (coupled problems)
      integer (I4B) nnat
!
!  auxillary
      real (DP) omega, userparam(10)
      character (len=20) physics
      character (len=10) :: whatsystem='Windows'
!
!  constants
      integer (I4B), parameter :: polymax = 20
      real (DP), parameter :: pi = 3.141592653589793238462643_DP
      real (DP), parameter :: mu0 = 0.1256637061435917295385057e-5_DP  ! free space permeability
      real (DP), parameter :: c0 = 299792458._DP                       ! speed of light
      real (DP), parameter :: eps0 = .8854187817620389850536565e-11_DP ! free space permittivity
      real (DP), parameter :: gravc = 6.67428e-11_DP                   ! gravitational constant
      real (DP), parameter :: Elch = 1.60217646e-19_DP                 ! elementary charge (in Coulombs)
      real (DP), parameter :: Kboltz = 8.617343e-5_DP                  ! Boltzmann's constant (in eV/K)
      real (DP), parameter :: Avogadro = 6.02214e23_DP                 ! The number of molecules in a mole (in 1/mol)
!
!
      end module globalvariables
