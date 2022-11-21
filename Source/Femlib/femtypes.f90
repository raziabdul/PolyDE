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
!    $Revision: 1.6 $
!    $Date: 2007/10/25 14:20:42 $
!    $Author: m_kasper $
!
module femtypes
       implicit none
       integer, parameter :: I4B = selected_int_kind(9)
       integer, parameter :: I2B = selected_int_kind(4)
       integer, parameter :: I1B = selected_int_kind(2)
       integer, parameter :: SP = kind(1.0)
       integer, parameter :: DP = kind(1.0d0)
       integer, parameter :: QP = kind(1.0q0)
       integer, parameter :: SPC = kind((1.0,1.0))
       integer, parameter :: DPC = kind((1.0d0,1.0d0))
       integer, parameter :: LGT = kind(.true.)
       type :: ARRPTRI
!  array of pointers
         integer(I4B), dimension(:), pointer :: d
       end type
       type :: ARRPTRR
!  array of pointers
         real(DP), dimension(:), pointer :: d
       end type
       type :: ARRPTRDPC
!  array of pointers
         complex(DPC), dimension(:), pointer :: d
       end type

      type :: ARRALLOI
!  array of allocatable
        integer(I4B), dimension(:), allocatable :: d
      end type

end module femtypes