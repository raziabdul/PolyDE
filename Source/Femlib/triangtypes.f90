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
!    $Revision: 1.3 $
!    $Date: 2006/07/07 21:08:37 $
!    $Author: r_abdul $
!
module triangtypes
      use femtypes
      implicit none
!  node of a double linked list
      type :: node
        integer (I4B) :: nodeindx    ! the index of the node
        real (DP) :: angle           ! inner angle 
        logical :: ear               ! is .true. if the node is an ear
        type (node), pointer :: next ! pointer to the next node
        type (node), pointer :: prev ! pointer to the prevoius node
        integer (I4B) :: name        ! a unique integer as an identifyer
      end type node
!  array of pointers to nodes
      type :: arrpt
        type (node), pointer :: p
      end type
end module triangtypes
