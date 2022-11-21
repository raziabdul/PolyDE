      subroutine dptransformpoly(xl,yl,zl,a,n,xo,yo,zo)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) a(4,4), xl(:), yl(:), zl(:), xo(:), yo(:), zo(:)
      integer n
      intent (in) xl, yl, zl, a, n
      intent (out) xo, yo, zo
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
!    $Date: 2006/07/07 22:13:33 $
!    $Author: r_abdul $
!
!  apply homogenous coordinate transformation to a set of points
!   a   Tranformation matrix
!   xl,yl,zl    input coordinates
!   xo,yo,zo    output coordinates
      integer i
!
      do i=1,n
        xo(i)=a(1,1)*xl(i)+a(1,2)*yl(i)+a(1,3)*zl(i)+a(1,4)
        yo(i)=a(2,1)*xl(i)+a(2,2)*yl(i)+a(2,3)*zl(i)+a(2,4)
        zo(i)=a(3,1)*xl(i)+a(3,2)*yl(i)+a(3,3)*zl(i)+a(3,4)
      end do
      return
      end subroutine dptransformpoly



      subroutine sptransformpoly(xl,yl,zl,a,n,xo,yo,zo)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) a(4,4)
      real (SP) xl(:), yl(:), zl(:), xo(:), yo(:), zo(:)
      integer n
      intent (in) xl, yl, zl, a, n
      intent (out) xo, yo, zo
!  apply homogenous coordinate transformation to a set of points
!   a   Tranformation matrix
!   xl,yl,zl    input coordinates
!   xo,yo,zo    output coordinates
      integer i
!
      do i=1,n
        xo(i)=a(1,1)*xl(i)+a(1,2)*yl(i)+a(1,3)*zl(i)+a(1,4)
        yo(i)=a(2,1)*xl(i)+a(2,2)*yl(i)+a(2,3)*zl(i)+a(2,4)
        zo(i)=a(3,1)*xl(i)+a(3,2)*yl(i)+a(3,3)*zl(i)+a(3,4)
      end do
      return
      end subroutine sptransformpoly



      subroutine scalehom(a,sx,sy,sz)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) a(4,4), sx, sy, sz
      intent (in) sx, sy, sz
      intent (inout) a
!  scale in homogeneous coordinates
      real (DP) b(4,4)
!
      b=0._DP
      b(1,1)=sx
      b(2,2)=sy
      b(3,3)=sz
      b(4,4)=1._DP
      a=matmul(b,a)
      return
      end subroutine scalehom



      subroutine transhom(a,tx,ty,tz)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) a(4,4), tx, ty, tz
      intent (in) tx, ty, tz
      intent (inout) a
!  Translation in homogeneous coordinates
      real (DP) b(4,4)
!
      b=0._DP
      b(1,4)=tx
      b(2,4)=ty
      b(3,4)=tz
      b(4,4)=1._DP
      b(1,1)=1._DP
      b(2,2)=1._DP
      b(3,3)=1._DP
      a=matmul(b,a)
      return
      end subroutine transhom



      subroutine rotatehom(a,theta,axis)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) a(4,4), theta
      integer axis
      intent (in) theta, axis
      intent (inout) a
!  Rotate about axis 1,2,3 in homogeneous coordinates
      real (DP) b(4,4)
!
      b=0._DP
      select case (axis)
      case (1)
! rotate about x-axis
        b(1,1)=1._DP
        b(2,2)= cos(theta)
        b(2,3)=-sin(theta)
        b(3,2)= sin(theta)
        b(3,3)= cos(theta)
      case (2)
! rotate about y-axis
        b(2,2)=1._DP
        b(1,1)= cos(theta)
        b(1,3)= sin(theta)
        b(3,1)=-sin(theta)
        b(3,3)= cos(theta)
      case (3)
! rotate about z-axis
        b(3,3)=1._DP
        b(1,1)= cos(theta)
        b(1,2)=-sin(theta)
        b(2,1)= sin(theta)
        b(2,2)= cos(theta)
      end select
      b(4,4)=1._DP
      a=matmul(b,a)
      return
      end subroutine rotatehom
