function append_dp(array,values)
    use femtypes
    implicit none
    real (DP), pointer :: array(:,:), values(:,:), append_dp(:,:)
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
!    $Date: 2013/07/24 13:41:25 $
!    $Author: juryanatzki $
!
!  resize an array: allocate the array of size n and copy the former content
!   input:
!       array      first array
!       values     seccond array, which will be appended to the first array
!   output:
!       append     an array, which has one common dimension and one that is a sum of the uncommon dimensions of incomming arrays
!
!    

    integer (I4B) :: num_arrow, num_arcol, num_valrow, num_valcol
    
    
     num_arrow = SIZE(array(:,1))
     num_arcol = SIZE(array(1,:))
    num_valrow = SIZE(values(:,1))
    num_valcol = SIZE(values(1,:))
    
    if (num_arrow .eq. num_valrow) then
        allocate(append_dp(num_arrow,num_arcol+num_valcol))
        append_dp(1:num_arrow,1:num_arcol) = array(1:num_arrow,1:num_arcol)
        append_dp(1:num_arrow,num_arcol+1:num_arcol+num_valcol) = values(1:num_valrow,1:num_valcol)
    end if
    
    if (num_arcol .eq. num_valcol) then ! <|--------- TO BE TESTED -!!
        allocate(append_dp(num_arrow+num_valrow,num_arcol))
        append_dp(1:num_arrow,1:num_arcol) = array(1:num_arrow,1:num_arcol)
        append_dp(num_arrow+1:num_arrow+num_valrow,1:num_arcol) = values(1:num_valrow,1:num_valcol)
    end if
    deallocate(array,values)
    return
end function append_dp



function append_int(array,values)
    use femtypes
    implicit none
    integer (I4B), pointer :: array(:,:), values(:,:), append_int(:,:)
!
!  resize an array: allocate the array of size n and copy the former content
!   input:
!       array      first array
!       values     seccond array, which will be appended to the first array
!   output:
!       append     an array, which has one common dimension and one that is a sum of the uncommon dimensions of incomming arrays
!
!

    integer (I4B) :: num_arrow, num_arcol, num_valrow, num_valcol
    
    
     num_arrow = SIZE(array(:,1))
     num_arcol = SIZE(array(1,:))
    num_valrow = SIZE(values(:,1))
    num_valcol = SIZE(values(1,:))
    
    if (num_arrow .eq. num_valrow) then
        allocate(append_int(num_arrow,num_arcol+num_valcol))
        append_int(1:num_arrow,1:num_arcol) = array(1:num_arrow,1:num_arcol)
        append_int(1:num_arrow,num_arcol+1:num_arcol+num_valcol) = values(1:num_valrow,1:num_valcol)
    end if
    
    if (num_arcol .eq. num_valcol) then ! <|--------- TO BE TESTED -!!
        allocate(append_int(num_arrow+num_valrow,num_arcol))
        append_int(1:num_arrow,1:num_arcol) = array(1:num_arrow,1:num_arcol)
        append_int(num_arrow+1:num_arrow+num_valrow,1:num_arcol) = values(1:num_valrow,1:num_valcol)
    end if
    deallocate(array,values)
    return
end function append_int


function append_vec_int(array,values)
    use femtypes
    implicit none
    integer (I4B), pointer :: array(:), values(:), append_vec_int(:)
!
!  resize an array: allocate the array of size n and copy the former content
!   input:
!       array      first array
!       values     seccond array, which will be appended to the first array
!   output:
!       append     an array, which has one common dimension and one that is a sum of the uncommon dimensions of incomming arrays
!
!

    integer (I4B) :: num_arrow, num_valrow
    
    num_arrow = SIZE(array(:))
    num_valrow = SIZE(values(:))
    
      allocate(append_vec_int(num_arrow+num_valrow))
      append_vec_int(1:num_arrow) = array(1:num_arrow)
      append_vec_int(num_arrow+1:num_arrow+num_valrow) = values(1:num_valrow)
   
   deallocate(array,values)
    return
end function append_vec_int


function append_vec_dp(array,values)
    use femtypes
    implicit none
    real (DP), pointer :: array(:), values(:), append_vec_dp(:)
!
!  resize an array: allocate the array of size n and copy the former content
!   input:
!       array      first array
!       values     seccond array, which will be appended to the first array
!   output:
!       append     an array, which has one common dimension and one that is a sum of the uncommon dimensions of incomming arrays
!
!

    integer (I4B) :: num_arrow, num_arcol, num_valrow, num_valcol
    
    num_arrow = SIZE(array(:))
    num_valrow = SIZE(values(:))
    
      allocate(append_vec_dp(num_arrow+num_valrow))
      append_vec_dp(1:num_arrow) = array(1:num_arrow)
      append_vec_dp(num_arrow+1:num_arrow+num_valrow) = values(1:num_valrow)
   
   deallocate(array,values)
    return
end function append_vec_dp