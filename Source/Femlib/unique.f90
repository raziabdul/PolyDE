! unique:
!
!This function receives a pointer to a array, removes multiple values from that array and
!returns the pointer to the shortened array
function unique(array_in)
    use femtypes
    use feminterface, only: reallocate
    implicit none
    integer (I4B), pointer :: array_in(:)
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
!    $Date: 2014/02/11 16:23:53 $
!    $Author: juryanatzki $
!
!  resize an array: allocate the array of size n and copy the former content
!   input:
!       array_in   incomming array, which may contain redundant values
!   output:
!       unique     pointer to the array, which misses manifold values from array_in. The last of the multiple entries survive.
!
    integer (I4B), pointer :: unique(:)
    integer nr_vals, size_array, idx_value
    
    nr_vals = 0
    size_array = size(array_in)
        allocate(unique(size_array))
        unique = array_in
        do idx_value = 1, size_array
            if ( any( array_in((idx_value+1):size_array) .eq. array_in(idx_value) ) ) then
            array_in(idx_value) = 0
            else
            nr_vals = nr_vals + 1
            unique(nr_vals) = array_in(idx_value)
            end if
        end do
        if (size_array .gt. nr_vals)  unique => reallocate(unique,nr_vals)
    deallocate (array_in)
    return
end function