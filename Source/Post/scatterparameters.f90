      subroutine scatterparameters(xa,ya,xe,ye,swr,rc,tc)
      use feminterface, only: fetchmatparameters, elemnt, xy2lam, field
      use femtypes
      use globalvariables
      use matconstants
      implicit none
      real (DP) :: xa, xe, ya, ye
      real (DP) :: swr, rc, tc
      intent (in) :: xa, ya, xe, ye
      intent (out) :: swr, rc, tc
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
!    $Revision: 1.8 $
!    $Date: 2014/07/15 13:06:21 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Determine the max and min value of the field and calculate SWR, Reflection
!  coefficient, Transmission coefficient
!
!  Input:
!     xa, ya    start coordinates of the line
!     xe, ye    end coordinates of the line
!
!  Output:
!     swr       Standing Wave Ratio
!     rc        Reflection coefficient
!     tc        Transmission coefficient
!
!  local variables:
      integer (I4B) :: i, j, ielem, nx, ny
      integer (I4B) :: matindex, matindex1, matindex2
      real (DP) :: mag, t_x, t_y, e_x, e_y, F_max, F_min
      real (DP) :: wavelen, delta, lambda(3)
      real (DP) :: mid_x, mid_y, mid_xl, mid_yl, mid_xr, mid_yr
      real (DP) :: x11, y11, x12, y12, x21, y21, x22, y22, d_x, d_y
      real (DP), allocatable :: xp(:,:), yp(:,:), swr_y(:)
      real (DP) :: sum_swr, mean_swr, sum_sqdiff, var_swr
      complex (DPC) :: z(15,1)
      logical :: ok, typ(5)
!
!  assign the number of intended points in x, y direction
      nx = 100
      ny = 50
!  compute the tangent of the intended line
      mag = sqrt((xe-xa)**2 + (ye-ya)**2)
      t_x = (ye-ya)/mag
      t_y = (xa-xe)/mag
!  compute the unit vector
      e_x = (xe-xa)/mag
      e_y = (ye-ya)/mag
!
!  find the wavelength (lambda) & the distance delta (lambda/2)
      wavelen = 2*pi*c0/omega
      delta = wavelen/2
!
!  find the material index at the middle of the intended line
      mid_x = (xa + xe)/2
      mid_y = (ya + ye)/2
      call elemnt(xn,yn,e,n,mid_x,mid_y,ielem,en,ok)
      if (ok) then
        matindex = matzif(geb(ielem))
      end if
!
!  find the material index at a distance delta on the left side
!  of the intended line
      mid_xl = mid_x - delta*t_x
      mid_yl = mid_y - delta*t_y
      call elemnt(xn,yn,e,n,mid_xl,mid_yl,ielem,en,ok)
!  if there is an element exists at the specified location
      if (ok) then
        matindex1 = matzif(geb(ielem))
!  define the left side coordinates of the section
!  if the material index is same
        if (matindex1 .eq. matindex) then
          x12 = xa - delta*t_x
          y12 = ya - delta*t_y
          x22 = xe - delta*t_x
          y22 = ye - delta*t_y
!  if the material index is not same
        else if (matindex1 .ne. matindex) then
!  shift to 2*delta distance on the right side
10        mid_xr = mid_x + 2*delta*t_x
          mid_yr = mid_y + 2*delta*t_y
          call elemnt(xn,yn,e,n,mid_xr,mid_yr,ielem,en,ok)
          if (ok) then
            matindex2 = matzif(geb(ielem))
!  if the material is same
            if (matindex2 .eq. matindex) then
              print *, " Intended section crosses the interface/boundary."
              print *, " Left side of the intended line is neglected."
!  assign left side coordinates as coordinates of intended line
              x12 = xa
              y12 = ya
              x22 = xe
              y22 = ye
!  assign right side coordinates at 2*delta distance
              x11 = xa + 2*delta*t_x
              y11 = ya + 2*delta*t_y
              x21 = xe + 2*delta*t_x
              y21 = ye + 2*delta*t_y
              goto 30           ! skip to calcualate SWR
            else
              print *, " Could not define a section."
              goto 40
            end if
          else if (.not. ok) then
            print *, " Could not define a section."
            goto 40
          end if
        end if
!  if there is no element exsists at the specified location
      else if (.not. ok) then
        goto 10
      end if
!      
!  find the material index at a distance delta on the right side
!  of the intended line
      mid_xr = mid_x + delta*t_x
      mid_yr = mid_y + delta*t_y
      call elemnt(xn,yn,e,n,mid_xr,mid_yr,ielem,en,ok)
!  if there is an element exists at the specified location
      if (ok) then
        matindex2 = matzif(geb(ielem))
!  define the right side coordinates of the section
!  if the material index is same
        if (matindex2 .eq. matindex) then
          x11 = xa + delta*t_x
          y11 = ya + delta*t_y
          x21 = xe + delta*t_x
          y21 = ye + delta*t_y
!  if the material index is not same
        else if (matindex2 .ne. matindex) then
!  shift to 2*delta distance on the left side
20        mid_xl = mid_x - 2*delta*t_x
          mid_yl = mid_y - 2*delta*t_y
          call elemnt(xn,yn,e,n,mid_xl,mid_yl,ielem,en,ok)
          if (ok) then
            matindex1 = matzif(geb(ielem))
!  if the material is same
            if (matindex1 .eq. matindex) then
              print *, " Intended section crosses the interface/boundary."
              print *, " Right side of the intended line is neglected."
!  assign left side coordinates at 2*delta distance
              x12 = xa - 2*delta*t_x
              y12 = ya - 2*delta*t_y
              x22 = xe - 2*delta*t_x
              y22 = ye - 2*delta*t_y
!  assign right side coordinates as coordinates of intended line
              x11 = xa
              y11 = ya
              x21 = xe
              y21 = ye
              goto 30         ! skip to calcualate SWR
            else
              print *, " Could not define a section."
              goto 40
            end if
          else if (.not. ok) then
            print *, " Could not define a section."
            goto 40
          end if
        end if
!  if there is no element exsists at the specified location
      else if (.not. ok) then
        goto 20
      end if   
!
30    typ = (/.true.,.false.,.false.,.false.,.false./)
!  initialize the sum of Standing Wave Ratio (SWR)
      sum_swr = 0._DP
      allocate (xp(ny,nx), yp(ny,nx), swr_y(ny))
!  
!  incremental distances in vertical & horizontal axis
      d_x = (x11-x12)/nx
      d_y = (y22-y12)/ny
!  step forward in vertical axis
      do i=1,ny
!  initialize maximum & minimum of field
        F_min = 10._DP
        F_max = 0._DP
!
!  step forward in horizontal axis
        do j=1,nx
          yp(i,j) = y12 + (i-0.5_DP)*d_y*e_y + (j-0.5_DP)*d_x*t_y
          xp(i,j) = x12 + (i-0.5_DP)*d_y*e_x + (j-0.5_DP)*d_x*t_x
!
          call elemnt(xn,yn,e,n,xp(i,j),yp(i,j),ielem,en,ok)
          if (ok) then
!  convert world coordinates to barycentric coordinates
            call xy2lam(xp(i,j),yp(i,j),ielem,lambda,xn,yn,e)
!  get field components
            call field(ielem,lambda,typ,z)
!  check & assign to F_max or F_min by comparision
!  with previous values
!  only one nature for TE- and TMWAVE
            F_max = max(abs(z(1,1)),F_max)
            F_min = min(abs(z(1,1)),F_min)
          end if
        end do
!  calculate the SWR & sum it
        swr_y(i) = F_max/F_min
        sum_swr = sum_swr + swr_y(i)
      end do
      deallocate(xp,yp)
!
!  compute the mean SWR
      mean_swr = sum_swr/ny
!  compute the sum of squared differences of SWR with the mean SWR
      sum_sqdiff = 0._DP
      do i=1,ny
        sum_sqdiff = sum_sqdiff + (swr_y(i)-mean_swr)**2
      end do
      deallocate(swr_y)
!
      var_swr = sum_sqdiff/ny       ! variance of SWR
      swr = mean_swr
!
!  compute the reflection and transmission coefficient if SWR is not zero
      if (swr .ne. 0) then
        rc = abs((swr - 1._DP)/(swr + 1._DP))
        tc = 1._DP - rc
!  assign the reflection and transmission coefficient is zero, if SWR = 0
      else
40      swr = 0
        rc = 0
        tc = 0
      end if
!
      return
      end subroutine scatterparameters
