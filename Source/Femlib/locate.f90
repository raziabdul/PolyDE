      function locate_DP(xx,x)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) xx(:),x
      integer(I4B) locate_DP
      intent(in) :: x,xx
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
!    $Revision: 1.4 $
!    $Date: 2006/07/07 21:08:36 $
!    $Author: r_abdul $
!
!  Search an ordered table by bisection
!       xx       ordered table either monotonic increasing or decreasing
!       x        value to locate in table
!       locate   index in table such that x is between x(locate) and x(locate+1)
!                locate = 0 or locate = n is returned to indcate that x is out of range
      integer(I4B) n,jl,jm,ju
      logical ascnd
!
      n=size(xx)
!  ascnd is true if ascending order of table 
      ascnd = (xx(n) .gt. xx(1))
!  initialize lower and upper limits
      jl=0
      ju=n+1
!  repeat  until condition is satisfied
      do
        if (ju-jl .le. 1) exit
!  compute midpoint
        jm=(ju+jl)/2
        if (ascnd .eqv. (x .gt. xx(jm))) then
!  replace lower limit
          jl=jm
        else
!  replace upper limit
          ju=jm
        end if
      end do
!  set output 
      locate_DP=jl
      return
      end function locate_DP



      function locate_idx_DP(xx,x,indx)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) xx(:),x
      integer(I4B) locate_idx_DP, indx(:)
      intent(in) :: x,xx,indx
!  Search an ordered table by bisection
!  Input:
!       xx       table of values
!       x        value to locate in table
!       indx     sorted list of indices, either monotonic increasing or decreasing
!  Output:
!       locate   index in table such that x is between x(locate) and x(locate+1)
!                locate = 0 or locate = n is returned to indcate that x is out of range
      integer(I4B) n,jl,jm,ju
      logical ascnd
!
      n=size(indx)
!  ascnd is true if ascending order of table 
      ascnd = (xx(indx(n)) .gt. xx(indx(1)))
!  initialize lower and upper limits
      jl=0
      ju=n+1
!  repeat  until condition is satisfied
      do
        if (ju-jl .le. 1) exit
!  compute midpoint
        jm=(ju+jl)/2
        if (ascnd .eqv. (x .gt. xx(indx(jm)))) then
!  replace lower limit
          jl=jm
        else
!  replace upper limit
          ju=jm
        end if
      end do
!  set output 
      locate_idx_DP=jl
      return
      end function locate_idx_DP
