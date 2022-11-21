module matconstants
!  global variabled of the module
      use femtypes
      integer (I4B) maxmat
      parameter (maxmat=30)
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
!    $Date: 2014/07/01 16:00:08 $
!    $Author: m_kasper $
!
!  module variables
!    maxmat       maximum number of material (size of arrays)
!    ipos         actual number of materials which have been found and stored up to now
!    matnames     names to the materials (array of size maxmat of strings), e.g. copper, polyimid, iron, ...
!    numparam     number of physical parameters 
!    parametername  names of physical parameters (array of size numparam of strings)  e.g.  epsr, poissonratio, tan_delta, ...
!    param        material parameters as a list of real values param(i)%d(j) with i being the index of the material
!
!
      save
      integer (I4B) :: ipos=-1, numparam
      character (len=25) :: matnames(maxmat)
      character (len=20), allocatable :: parameternames(:)
      character (len=12), allocatable :: parameterunit(:)
      real (DP), allocatable :: parameterdefault(:)
      logical :: usermat(maxmat)
      type (ARRPTRR) :: param(0:maxmat)
end module matconstants
