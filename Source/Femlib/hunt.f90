      subroutine hunt(xx,x,jlo)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: jlo
      real (DP) :: x, xx(:)
      intent (in) :: x, xx
      intent (inout) :: jlo
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
!    $Date: 2006/07/19 10:54:28 $
!    $Author: m_schober $
!
!
!  This routine was adopted from Numerical Recipes in FORTAN 90:
!  Given an array xx(1:n), and given a value x, returns a value jlo such that x
!  is between xx(jlo) and xx(jlo+1). xx must be monotonic, either increasing or
!  decreasing. jlo = 0 or jlo = n is returned to indicate, that x is out of
!  range. jlo on input is taken as the initial guess for jlo on output.
!
!  Input:
!    x              value to search for / to be bracketed
!    xx             list of ordered double precision values
!
!  In-/Output:
!    jlo            initial value for the search as input
!                   x is found at a position between xx(jlo) and xx(jlo+1)
!
!  Internal variables:
      integer (I4B) :: n, inc, jhi, jm
      logical :: ascnd = .false.
!
!  determine size of array xx
      n = size(xx)
!  check if list is ascending or descending
      if (xx(n) .ge. xx(1)) ascnd = .true.
!
      if ((jlo .le. 0) .or. (jlo .gt. n)) then
!  usual bisecting if jlo not useful
        jlo = 0
        jhi = n+1
      else
!  set increment for hunt
        inc = 1
        if ((x .ge. xx(jlo)) .eqv. ascnd) then
!  HUNT UP
          do
            jhi = jlo + inc
            if (jhi .gt. n) then
!  hunting finished
              jhi = n + 1
              exit
            else
              if ((x .lt. xx(jhi)) .eqv. ascnd) exit
!  hunting not finished => double increment
              jlo = jhi
              inc = 2*inc
            end if
          end do
        else
!  HUNT DOWN
          jhi = jlo
          do
            jlo = jhi - inc
            if (jlo .lt. 1) then
!  hunting finished
              jlo = 0
              exit
            else
              if ((x .ge. xx(jlo)) .eqv. ascnd) exit
!  hunting not finished => double increment
              jhi = jlo
              inc = 2*inc
            end if
          end do
        end if
      end if
!  hunting finished, value bracketed. begin bisecting
      do
        if (jhi-jlo .le. 1) then
          if (x .eq. xx(n)) jlo = n - 1
          if (x .eq. xx(1)) jlo = 1
          exit
        else
          jm = floor((jhi + jlo)/2._DP)
          if ((x .ge. xx(jm)) .eqv. ascnd) then
            jlo = jm
          else
            jhi = jm
          end if
        end if
      end do
!
      end subroutine hunt