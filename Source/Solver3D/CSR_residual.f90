      subroutine CSR_residual(a,b,x,ia,ja,n,res)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) a(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: a, b, x, ia, ja, n
      intent (out) :: res

!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 - 2015 Institute for Micro Systems Technology,
!                              Hamburg University of Technology.
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
!------------------------------------------------------------------------------
!
!    $Revision: 1.2 $
!    $Date: 2015/06/23 13:53:22 $
!    $Author: jimmykamboh $
!
!------------------------------------------------------------------------------

!  local variables
      integer (I4B) i, j, k
!
!  Comutation of residual:  res = b - A*x
!
!    a       matrix (compact storage)
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    res     residual vector
!
      res(:)=b(:)
      do i=1,n
        do j=ia(i),ia(i+1)-1
!  k = column
          k=ja(j)
          res(i)=res(i)-a(j)*x(k)
        end do
      end do
      return
      end subroutine CSR_residual
!
!
!
      subroutine LCSR_residual(lower,upper,diag,b,x,ia,ja,n,res)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) lower(:), upper(:), diag(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
!  local variables
      integer (I4B) i, j, k
!
!  Comutation of residual:  res = b - A*x
!
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    res     residual vector
!
      res(1)=b(1)-diag(1)*x(1)
      do i=2,n
        res(i)=b(i)-diag(i)*x(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          res(i)=res(i)-lower(j)*x(k)
          res(k)=res(k)-upper(j)*x(i)
        end do
      end do
      return
      end subroutine LCSR_residual