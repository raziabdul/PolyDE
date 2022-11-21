      pure subroutine invertmat_dpc(matrix,inverted)
      use feminterface, only:
      use femtypes
      implicit none
      complex(DPC), intent(IN) :: matrix(2,2)
      complex(DPC), intent(OUT) :: inverted(2,2)
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
!    $Revision: 1.6 $
!    $Date: 2010/10/25 11:56:58 $
!    $Author: chvokas $
!------------------------------------------------------------------------------
!
!  Invert a quadratic matrix A to obtain A^(-1).
!      / a b \               1    /  d -b \
!  A = |     |   A^(-1) = ------- |       |
!      \ c d /            ad - bc \ -c  a /
!
!------------------------------------------------------------------------------
!  Input:
!     matrix    quadratic (2 x 2) matrix A
!
!  Output:
!     inverted  inverse matrix A^(-1)
!
!  Internal variables:
      complex(DPC) :: detmatrix
!------------------------------------------------------------------------------
!
!
!  Calculate determinant of matrix
      detmatrix = matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
!
!  Assign values for inverted matrix
      inverted(1,1) =  matrix(2,2)
      inverted(1,2) = -matrix(1,2)
      inverted(2,1) = -matrix(2,1)
      inverted(2,2) =  matrix(1,1)
!
!  calculate inverted matrix
      inverted = inverted / detmatrix
!
      end subroutine invertmat_dpc
!
!
!
      pure subroutine invertmat_r(matrix,inverted)
      use feminterface, only:
      use femtypes
      implicit none
      real(DP), intent(IN) :: matrix(2,2)
      real(DP), intent(OUT) :: inverted(2,2)
!
!------------------------------------------------------------------------------
!
!  Invert a quadratic matrix A to obtain A^(-1).
!      / a b \               1    /  d -b \
!  A = |     |   A^(-1) = ------- |       |
!      \ c d /            ad - bc \ -c  a /
!
!------------------------------------------------------------------------------
!  Input:
!     matrix    quadratic (2 x 2) matrix A
!
!  Output:
!     inverted  inverse matrix A^(-1)
!
!  Internal variables:
      real(DP) :: detmatrix
!------------------------------------------------------------------------------
!
!
!  Calculate determinant of matrix
      detmatrix = matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
!
!  Assign values for inverted matrix
      inverted(1,1) =  matrix(2,2)
      inverted(1,2) = -matrix(1,2)
      inverted(2,1) = -matrix(2,1)
      inverted(2,2) =  matrix(1,1)
!
!  calculate inverted matrix
      inverted = inverted / detmatrix
!
      end subroutine invertmat_r
!
!
!
      pure subroutine invertmat3_dpc(matrix,inverted)
      use feminterface, only:
      use femtypes
      implicit none
      complex(DPC), intent(IN) :: matrix(3,3)
      complex(DPC), intent(OUT) :: inverted(3,3)
!
!------------------------------------------------------------------------------
!
!  Invert a 3 x 3 matrix A to obtain A^(-1).
!  A^(-1)=1/(det A) * (Ad A)^T
!
!------------------------------------------------------------------------------
!  Input:
!     matrix    3 x 3 matrix A
!
!  Output:
!     inverted  inverse matrix A^(-1)
!
!  Internal variables:
      complex(DPC) :: detmatrix, admatrix(3,3)
!------------------------------------------------------------------------------
!
!
!  Calculate determinant of matrix by Sarrus law
      detmatrix =   matrix(1,1)*matrix(2,2)*matrix(3,3) &
                  + matrix(1,2)*matrix(2,3)*matrix(3,1) &
                  + matrix(1,3)*matrix(2,1)*matrix(3,2) &
                  - matrix(3,1)*matrix(2,2)*matrix(1,3) &
                  - matrix(3,2)*matrix(2,3)*matrix(1,1) &
                  - matrix(3,3)*matrix(2,1)*matrix(1,2)
!
!  Assign values for Ad A
      admatrix(1,1) =   matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
      admatrix(1,2) = -(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
      admatrix(1,3) =   matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
      admatrix(2,1) = -(matrix(1,2)*matrix(3,3)-matrix(3,2)*matrix(1,3))
      admatrix(2,2) =   matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
      admatrix(2,3) = -(matrix(1,1)*matrix(3,2)-matrix(3,1)*matrix(1,2))
      admatrix(3,1) =   matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
      admatrix(3,2) = -(matrix(1,1)*matrix(2,3)-matrix(2,1)*matrix(1,3))
      admatrix(3,3) =   matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
!
!  calculate inverted matrix
      inverted = transpose(admatrix) / detmatrix
!
      end subroutine invertmat3_dpc
!
!
!
      pure subroutine invertmat3_r(matrix,inverted)
      use feminterface, only:
      use femtypes
      implicit none
      real(DP), intent(IN) :: matrix(3,3)
      real(DP), intent(OUT) :: inverted(3,3)
!
!------------------------------------------------------------------------------
!
!  Invert a 3 x 3 matrix A to obtain A^(-1).
!  A^(-1)=1/(det A) * (Ad A)^T
!
!------------------------------------------------------------------------------
!  Input:
!     matrix    3 x 3 matrix A
!
!  Output:
!     inverted  inverse matrix A^(-1)
!
!  Internal variables:
      real(DP) :: detmatrix, admatrix(3,3)
!------------------------------------------------------------------------------
!
!
!  Calculate determinant of matrix by Sarrus law
      detmatrix =   matrix(1,1)*matrix(2,2)*matrix(3,3) &
                  + matrix(1,2)*matrix(2,3)*matrix(3,1) &
                  + matrix(1,3)*matrix(2,1)*matrix(3,2) &
                  - matrix(3,1)*matrix(2,2)*matrix(1,3) &
                  - matrix(3,2)*matrix(2,3)*matrix(1,1) &
                  - matrix(3,3)*matrix(2,1)*matrix(1,2)
!
!  Assign values for Ad A
      admatrix(1,1) =   matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
      admatrix(1,2) = -(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
      admatrix(1,3) =   matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
      admatrix(2,1) = -(matrix(1,2)*matrix(3,3)-matrix(3,2)*matrix(1,3))
      admatrix(2,2) =   matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
      admatrix(2,3) = -(matrix(1,1)*matrix(3,2)-matrix(3,1)*matrix(1,2))
      admatrix(3,1) =   matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
      admatrix(3,2) = -(matrix(1,1)*matrix(2,3)-matrix(2,1)*matrix(1,3))
      admatrix(3,3) =   matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
!
!  calculate inverted matrix
      inverted = transpose(admatrix) / detmatrix
!
      end subroutine invertmat3_r
