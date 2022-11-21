      subroutine energyflow(xa,ya,xe,ye,flux)
      use feminterface, only: elemnt, fieldquantity, getpostsetting, xy2lam
      use femtypes
      use globalvariables
      implicit none
      real (DP) :: xa, ya, xe, ye
      real (DP) :: flux
      intent (in) :: xa, ya, xe, ye
      intent (out) :: flux
!
!------------------------------------------------------------------------------
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
!    $Revision: 1.12 $
!    $Date: 2014/07/15 13:06:21 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Determine the energy flow through a port which is given by a line between
!  (startx,starty) and (endx,endy). The energy is computed by:
!
!    Re{int[S.n] da} = 1/2 Re{int[(ExH*).n] da}
!
!  The subroutine calc_field provides the Re{Sx} and Re{Sy}.
!
!------------------------------------------------------------------------------
!  Input:
!     xa, ya    start coordinates of the line
!     xe, ye    end coordinates of the line
!
!  Output:
!     flux      Energy flow
!
!  local variables:
      integer (I4B) :: elem, div, i
      real (DP) :: dirvec(2), endp(2), nvec(2), startp(2), unitvec(2)
      real (DP) :: delta, dist
      real (DP) :: lambda(3), phi
      real (DP), allocatable :: points(:,:), distval(:)
      complex (DPC) :: poynting(2)
      character (len=10) :: unit
      character (len=50) :: descriptor
      logical :: ok, ok1
!------------------------------------------------------------------------------
!
!  assign start and end points
      startp = (/xa,ya/)
      endp = (/xe,ye/)
!
!  compute direction vector from start to end
      dirvec = endp - startp
!  compute length of vector and unitvector (length 1)
      dist = sqrt(dirvec(1)**2 + dirvec(2)**2)
      unitvec = dirvec/dist
!
!  compute normal unit vector (perpendicular to direction vector)
      nvec = (/unitvec(2),-unitvec(1)/)
!
!  compute the segment length
      call getpostsetting('DIVISION',div)
!  compute delta depending on # of divisions on length
      delta = dist/div
!
!  allocate and initialize distance value and points
      allocate (distval(div+1), points(2,div+1))
      distval = 0._DP
      points = 0._DP
!
!  get phase shift for calc_field routine
      call getpostsetting('PHI',phi)
!
!  initialize flux
      flux = 0._DP
!
!  numerical integration loop for div+1 points
      do i = 1,div+1
!  compute x,y coordinates for new point depending on # of division
        distval(i) = (i-1)*delta
        points(:,i)= startp + (i-1)*delta*unitvec
!  get element for point
        call elemnt(xn,yn,e,n,points(1,i),points(2,i),elem,en,ok)
        if (ok) then
          call xy2lam(points(1,i),points(2,i),elem,lambda,xn,yn,e)
!  get Poynting vector (nature=1)
          call fieldquantity(elem,'SX',lambda,phi,poynting(1),descriptor,unit,ok1)
          if (.not.ok1) then
            print *,"***** Poynting vector not available for ",trim(physics)
            stop
          end if
          call fieldquantity(elem,'SY',lambda,phi,poynting(2),descriptor,unit,ok1)
!  compute energy flow normal to unitvec
          flux = flux + dot_product(real(poynting,DP),nvec)
        else
          print *," **** Error in energyflow:"
          print "(a,g10.3,a,g10.3,a)"," point (",points(1,i),";",points(2,i),") not in an element"
          cycle
        end if
      end do
      deallocate(distval, points)
!
!  multiply by length between start and end point
      flux = flux*dist/(div+1)
      print *,"Energy flow through port in W: ",flux
!
      return
      end subroutine energyflow