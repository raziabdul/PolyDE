      subroutine pointvalue3D
      use femtypes
      use feminterface, only: getpostsetting
      use feminterface3D, only: fieldquantity3d, findelement3D, xyz2lam
      use globalvariables3D, only: nnat, nod, vn, vv, numv, physics
      implicit none
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
!    $Date: 2014/07/30 12:45:31 $
!    $Author: m_kasper $
!
!
!
!  Evaluate Filed value at a point
!
!  Input:
!     fieldtype          field quantity to compute (e.g. POTENTIAL, EY ...)
!     phi                phase / time instant to compute

!  Internal variables:
      integer (I4B) :: elem
      real (DP) :: phi, point(3), lambda(4), startx, starty,  startz
      complex(DPC) :: zs
      logical found, ok
      character (len=20) :: fieldtype
      character (len=50) :: descriptor
      character (len=12) :: unit

!
      call getpostsetting('FIELDTYPE',fieldtype)
      call getpostsetting('STARTX', startx)
      call getpostsetting('STARTY', starty)
      call getpostsetting('STARTZ', startz)
      call getpostsetting('PHI',phi)
      point=(/startx, starty, startz/)
      call findelement3D(point,nod,vn,vv,numv,elem,found)
      if (.not.found) then
        print *,"Point is not within the computational domain"
      else

!  get barycentric coordinates lambda of point in the element
        call xyz2lam(lambda, elem, point, nod, vn)
        call fieldquantity3d(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
        if (ok) then
          print '(a,es10.3,a,es10.3,a)','field value of '//trim(fieldtype)//' at point (',startx,',',starty,') is:'
          print *,'real part     : ',real(zs,DP)
          print *,'imaginary part: ',aimag(zs)
          print *,'absolute value: ',abs(zs)
        else
          print *,"***** fieldtype: ",trim(fieldtype)," unknown for ",trim(physics)
          print*, ' failed to evalute fieldquantity'
        end if

      end if
!
      end subroutine
