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

! NOTE: ALL ROUTINES AND FUNCTIONS IN THIS FILE ARE DEPRECATED AND SHOULDN'T BE USED ANYMORE
!
! Please use the routines defined in transformtensor.f90 instead: transformtensor, mat2ten, ten2mat
!
! Function trans_matrix has a new implementation in transformtensor.f90 and therefore commented here

!******************************************************************************
!  GENERATION OF TRANSFORMATION MATRIX ('OLD' PRINCIPAL TO 'NEW' REFERENCE)
!******************************************************************************
!
!      function trans_matrix(x1,x2)
!      use femtypes
!      implicit none
!      real (DP) :: x1(3), x2(3), trans_matrix(3,3)
!      intent (in) :: x1, x2
!!
!!
!!  Compute transformation matrix for given 'old' vectors x1 and x2 with respect
!!  to 'new' unit vectors x,y and z.
!!
!!  Input:
!!     x1,x2         vectors of first and second mutually orthogonal direction
!!
!!  Internal variables:
!      real (DP) :: absxn(3), x3(3), err
!      integer (I4B), parameter :: x(3) = (/1, 0, 0/), y(3) = (/0, 1, 0/), & 
!                                  z(3) = (/0, 0, 1/)
!!
!!
!!  compute absolute values / length of x1 and x2 vector
!      absxn(1) = sqrt(x1(1)**2 + x1(2)**2 + x1(3)**2)
!      absxn(2) = sqrt(x2(1)**2 + x2(2)**2 + x2(3)**2)
!
!!  check given vectors
!      err = dot_product(x1,x2) / (absxn(1)*absxn(2))
!      if (err .ne. 0._DP) then
!!  'Critical Pure'
!        print '(a)', ' '
!        print '(a)', 'TRANS_MATRIX - Fatal error!'
!        print '(a)', '  Given vectors are not orthogonal.'
!        stop
!      end if
!
!!  compute third direction vector
!      x3(1)= x1(2)*x2(3)-x1(3)*x2(2)
!      x3(2)= x1(3)*x2(1)-x1(1)*x2(3)
!      x3(3)= x1(1)*x2(2)-x1(2)*x2(1)
!
!!  compute absolute values / length of x3 vector
!      absxn(3) = sqrt(x3(1)**2 + x3(2)**2 + x3(3)**2)
!
!!  compute matrix of direction cosines / transformation matrix by
!      trans_matrix(1,1) = dot_product(x,x1) / absxn(1)
!      trans_matrix(1,2) = dot_product(x,x2) / absxn(2)
!      trans_matrix(1,3) = dot_product(x,x3) / absxn(3)
!      trans_matrix(2,1) = dot_product(y,x1) / absxn(1)
!      trans_matrix(2,2) = dot_product(y,x2) / absxn(2)
!      trans_matrix(2,3) = dot_product(y,x3) / absxn(3)
!      trans_matrix(3,1) = dot_product(z,x1) / absxn(1)
!      trans_matrix(3,2) = dot_product(z,x2) / absxn(2)
!      trans_matrix(3,3) = dot_product(z,x3) / absxn(3)
!!
!      end function trans_matrix
!
!
!
!******************************************************************************
!  TRANSFORMATION OF TENSORS (RANK 1 ... 4)
!******************************************************************************
!
      pure function transten_r1_dp(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3), b(3), transten_r1_dp(3)
      intent (in) :: a, b
!
!
!  Compute transformation of a rank 1 real tensor (vector)
!
!  Input:
!     a             transformation matrix
!     b             vector to transform
!
!
!  Compute transformed vector by matrix multiplication
      transten_r1_dp = matmul(a,b)
!
      end function transten_r1_dp
!------------------------------------------------------------------------------
!
!
!
      pure function transten_r2_dp(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3), b(3,3), transten_r2_dp(3,3)
      intent (in) :: a, b
!
!
!  Compute transformation of a rank 2 real tensor
!
!  Input:
!     a             transformation matrix
!     b             rank 2 tensor to transform
!
!  Internal variables:
      integer (I4B) :: i, j, k, l
      real (DP) :: sum1
!
!
!  Compute transformed rank 2 tensor
      do i = 1,3
        do j = 1,3
          sum1 = 0._DP
          do k = 1,3
            do l = 1,3
              sum1 = sum1 + a(i,k) * a(j,l) * b(k,l)
            end do
          end do
          transten_r2_dp(i,j) = sum1
        end do
      end do
!
      end function transten_r2_dp
!------------------------------------------------------------------------------
!
!
!
      pure function transten_r2_dpc(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3)
      complex (DPC) :: b(3,3), transten_r2_dpc(3,3)
      intent (in) :: a, b
!
!
!  Compute transformation of a rank 2 complex tensor
!
!  Input:
!     a             transformation matrix
!     b             rank 2 tensor to transform
!
!  Internal variables:
      integer (I4B) :: i, j, k, l
      complex (DPC) :: sum1
!
!
!  Compute transformed rank 2 tensor
      do i = 1,3
        do j = 1,3
          sum1 = 0._DP
          do k = 1,3
            do l = 1,3
              sum1 = sum1 + a(i,k) * a(j,l) * b(k,l)
            end do
          end do
          transten_r2_dpc(i,j) = sum1
        end do
      end do
!
      end function transten_r2_dpc
!------------------------------------------------------------------------------
!
!
!
      pure function transten_r3_dp(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3), b(3,3,3), transten_r3_dp(3,3,3)
      intent (in) :: a, b
!
!
!  Compute transformation of a rank 3 real tensor
!
!  Input:
!     a             transformation matrix
!     b             rank 3 tensor to transform
!
!  Internal variables:
      integer (I4B) :: i, j, k, l, m, n
      real (DP) :: sum1
!
!
!  Compute transformed rank 3 tensor
      do i = 1,3
        do j = 1,3
          do k = 1,3
            sum1 = 0._DP
            do l = 1,3
              do m = 1,3
                do n = 1,3
                  sum1 = sum1 + a(i,l) * a(j,m) * a(k,n) * b(l,m,n)
                end do
              end do
            end do
            transten_r3_dp(i,j,k) = sum1
          end do
        end do
      end do
!
      end function transten_r3_dp
!------------------------------------------------------------------------------
!
!
!
      pure function transten_r3_dpc(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3)
      complex (DPC) :: b(3,3,3), transten_r3_dpc(3,3,3)
      intent (in) :: a, b
!
!
!  Compute transformation of a rank 3 complex tensor
!
!  Input:
!     a             transformation matrix
!     b             rank 3 tensor to transform
!
!  Internal variables:
      integer (I4B) :: i, j, k, l, m, n
      complex (DPC) :: sum1
!
!
!  Compute transformed rank 3 tensor
      do i = 1,3
        do j = 1,3
          do k = 1,3
            sum1 = 0._DP
            do l = 1,3
              do m = 1,3
                do n = 1,3
                  sum1 = sum1 + a(i,l) * a(j,m) * a(k,n) * b(l,m,n)
                end do
              end do
            end do
            transten_r3_dpc(i,j,k) = sum1
          end do
        end do
      end do
!
      end function transten_r3_dpc
!------------------------------------------------------------------------------
!
!
!
      pure function transten_r4_dp(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3), b(3,3,3,3), transten_r4_dp(3,3,3,3)
      intent (in) :: a, b
!
!
!  Compute transformation of a rank 4 real tensor
!
!  Input:
!     a             transformation matrix
!     b             rank 4 tensor to transform
!
!  Internal variables:
      integer (I4B) :: i, j, k, l, m, n, o, p
      real (DP) :: sum1
!
!
!  Compute transformed rank 4 tensor
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              sum1 = 0._DP
              do m = 1,3
                do n = 1,3
                  do o = 1,3
                    do p = 1,3
                      sum1 = sum1 + a(i,m) * a(j,n) * a(k,o) * a(l,p) * b(m,n,o,p)
                    end do
                  end do
                end do
              end do
              transten_r4_dp(i,j,k,l) = sum1
            end do
          end do
        end do
      end do
!
      end function transten_r4_dp
!------------------------------------------------------------------------------
!
!
!
      pure function transten_r4_dpc(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3,3)
      complex (DPC) :: b(3,3,3,3), transten_r4_dpc(3,3,3,3)
      intent (in) :: a, b
!
!
!  Compute transformation of a rank 4 complex tensor
!
!  Input:
!     a             transformation matrix
!     b             rank 4 tensor to transform
!
!  Internal variables:
      integer (I4B) :: i, j, k, l, m, n, o, p
      complex (DPC) :: sum1
!
!
!  Compute transformed rank 4 tensor
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              sum1 = 0._DP
              do m = 1,3
                do n = 1,3
                  do o = 1,3
                    do p = 1,3
                      sum1 = sum1 + a(i,m) * a(j,n) * a(k,o) * a(l,p) * b(m,n,o,p)
                    end do
                  end do
                end do
              end do
              transten_r4_dpc(i,j,k,l) = sum1
            end do
          end do
        end do
      end do
!
      end function transten_r4_dpc
!
!
!
!******************************************************************************
!  MATRIX TO TENSOR NOTATION (RANK 2 ... 4)
!
!  (NOTICE THAT THIS PROGRAM PART IS ONLY VALID FOR USUAL NOTATION. THERE IS
!   NO FACTOR OF 0.5 TAKEN INTO ACCOUNT)
!******************************************************************************
!
      pure function mat2ten_r2_dp(mat)
      use femtypes
      implicit none
      real (DP) :: mat(6), mat2ten_r2_dp(3,3)
      intent (in) :: mat
!
!
!  Transform a vector (matrix notation) to a 3x3 tensor (tensor notation).
!  CAUTION: We assume, that the tensor is symmetric!
!
!
      mat2ten_r2_dp(1,:) = (/mat(1), mat(6), mat(5)/)
      mat2ten_r2_dp(2,:) = (/mat(6), mat(2), mat(4)/)
      mat2ten_r2_dp(3,:) = (/mat(5), mat(4), mat(3)/)
!
      end function mat2ten_r2_dp
!------------------------------------------------------------------------------
!
!
!
      pure function mat2ten_r2_dpc(mat)
      use femtypes
      implicit none
      complex (DPC) :: mat(6), mat2ten_r2_dpc(3,3)
      intent (in) :: mat
!
!
!  Transform a vector (matrix notation) to a 3x3 tensor (tensor notation).
!  CAUTION: We assume, that the tensor is symmetric!
!
!
      mat2ten_r2_dpc(1,:) = (/mat(1), mat(6), mat(5)/)
      mat2ten_r2_dpc(2,:) = (/mat(6), mat(2), mat(4)/)
      mat2ten_r2_dpc(3,:) = (/mat(5), mat(4), mat(3)/)
!
      end function mat2ten_r2_dpc
!------------------------------------------------------------------------------
!
!
!
      pure function mat2ten_r3_dp(mat)
      use femtypes
      implicit none
      real (DP) :: mat(6,3), mat2ten_r3_dp(3,3,3)
      intent (in) :: mat
!
!
!  Transform 6,3 matrix (matrix notation) to a rank 3 tensor (tensor notation).
!
!  Internal variables:
      integer (I4B) :: i
!
!
      do i = 1,3
        mat2ten_r3_dp(i,1,:) = (/mat(1,i), mat(6,i), mat(5,i)/)
        mat2ten_r3_dp(i,2,:) = (/mat(6,i), mat(2,i), mat(4,i)/)
        mat2ten_r3_dp(i,3,:) = (/mat(5,i), mat(4,i), mat(3,i)/)
      end do
!
      end function mat2ten_r3_dp
!------------------------------------------------------------------------------
!
!
!
      pure function mat2ten_r3_dpc(mat)
      use femtypes
      implicit none
      complex (DPC) :: mat(6,3), mat2ten_r3_dpc(3,3,3)
      intent (in) :: mat
!
!
!  Transform 6,3 matrix (matrix notation) to a rank 3 tensor (tensor notation).
!
!  Internal variables:
      integer (I4B) :: i
!
!
      do i = 1,3
        mat2ten_r3_dpc(i,1,:) = (/mat(1,i), mat(6,i), mat(5,i)/)
        mat2ten_r3_dpc(i,2,:) = (/mat(6,i), mat(2,i), mat(4,i)/)
        mat2ten_r3_dpc(i,3,:) = (/mat(5,i), mat(4,i), mat(3,i)/)
      end do
!
      end function mat2ten_r3_dpc
!------------------------------------------------------------------------------
!
!
!
      pure function mat2ten_r4_dp_old(mat)
      use femtypes
      implicit none
      real (DP) :: mat(6,6), mat2ten_r4_dp_old(3,3,3,3)
      intent (in) :: mat
!
!
!  Transform a 6x6 matrix (matrix notation) to a rank 4 tensor (tensor notation).
!  CAUTION: We assume, that the tensor is symmetric!
!
!  Internal variables:
      integer (I4B) :: i, j
      integer (I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
      do i = 1,6
        do j = 1,6
          mat2ten_r4_dp_old(a(i),b(i),a(j),b(j)) = mat(i,j)
          mat2ten_r4_dp_old(a(i),b(i),b(j),a(j)) = mat(i,j)
          mat2ten_r4_dp_old(b(i),a(i),b(j),a(j)) = mat(i,j)
          mat2ten_r4_dp_old(b(i),a(i),a(j),b(j)) = mat(i,j)
        end do
      end do
!
      end function mat2ten_r4_dp_old
!------------------------------------------------------------------------------
!
!
!
      pure function mat2ten_r4_dpc_old(mat)
      use femtypes
      implicit none
      complex (DPC) :: mat(6,6), mat2ten_r4_dpc_old(3,3,3,3)
      intent (in) :: mat
!
!
!  Transform a 6x6 matrix (matrix notation) to a rank 4 tensor (tensor notation).
!  CAUTION: We assume, that the tensor is symmetric!
!
!  Internal variables:
      integer (I4B) :: i, j
      integer (I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
      do i = 1,6
        do j = 1,6
          mat2ten_r4_dpc_old(a(i),b(i),a(j),b(j)) = mat(i,j)
          mat2ten_r4_dpc_old(a(i),b(i),b(j),a(j)) = mat(i,j)
          mat2ten_r4_dpc_old(b(i),a(i),b(j),a(j)) = mat(i,j)
          mat2ten_r4_dpc_old(b(i),a(i),a(j),b(j)) = mat(i,j)
        end do
      end do
!
      end function mat2ten_r4_dpc_old
!
!
!
!******************************************************************************
!  TENSOR TO MATRIX NOTATION
!
!  (NOTICE THAT THIS PROGRAM PART IS ONLY VALID FOR USUAL NOTATION. THERE IS
!   NO FACTOR OF 2 TAKEN INTO ACCOUNT)
!******************************************************************************
!
      function ten2mat_r2_dp(ten)
      use femtypes
      implicit none
      real (DP) :: ten(3,3), ten2mat_r2_dp(6)
      intent (in) :: ten
!
!
!  Transform a rank 2 tensor (tensor notation) to a vector (matrix notation).
!  CAUTION: We assume, that the tensor is symmetric!
!
!  Internal variables:
      real (DP) :: sum1
!
!
!  Check for symmetry (difference of symmetric entries must be 0)
      sum1 = ten(1,2)-ten(2,1) + ten(1,3)-ten(3,1) + ten(2,3)-ten(3,2)
      if (sum1 .ne. 0._DP) then
!  "Critical Pure"
        print '(a)', ' '
        print '(a)', 'TEN2MAT - Fatal error!'
        print '(a)', '  Tensor not symmetric.'
        stop
      end if
!
      ten2mat_r2_dp(1)=ten(1,1)
      ten2mat_r2_dp(2)=ten(2,2)
      ten2mat_r2_dp(3)=ten(3,3)
      ten2mat_r2_dp(4)=ten(3,2)
      ten2mat_r2_dp(5)=ten(1,3)
      ten2mat_r2_dp(6)=ten(1,2)
!
      end function ten2mat_r2_dp
!------------------------------------------------------------------------------
!
!
!
      function ten2mat_r2_dpc(ten)
      use femtypes
      implicit none
      complex (DPC) :: ten(3,3), ten2mat_r2_dpc(6)
      intent (in) :: ten
!
!
!  Transform a rank 2 tensor (tensor notation) to a vector (matrix notation).
!  CAUTION: We assume, that the tensor is symmetric!
!
!  Internal variables:
      complex (DPC) :: sum1
!
!
!  Check for symmetry (difference of symmetric entries must be 0)
      sum1 = ten(1,2)-ten(2,1) + ten(1,3)-ten(3,1) + ten(2,3)-ten(3,2)
      if (sum1 .ne. 0._DPC) then
!  "Critical Pure"
        print '(a)', ' '
        print '(a)', 'TEN2MAT - Fatal error!'
        print '(a)', '  Tensor not symmetric.'
        stop
      end if
!
      ten2mat_r2_dpc(1)=ten(1,1)
      ten2mat_r2_dpc(2)=ten(2,2)
      ten2mat_r2_dpc(3)=ten(3,3)
      ten2mat_r2_dpc(4)=ten(3,2)
      ten2mat_r2_dpc(5)=ten(1,3)
      ten2mat_r2_dpc(6)=ten(1,2)
!
      end function ten2mat_r2_dpc
!------------------------------------------------------------------------------
!
!
!
      function ten2mat_r3_dp(ten)
      use femtypes
      implicit none
      real (DP) :: ten(3,3,3), ten2mat_r3_dp(6,3)
      intent (in) :: ten
!
!
!  Transform a rank 3 tensor (tensor notation) to a 6,3 matrix (matrix notation).
!
!  Internal variables:
      integer (I4B) :: i, j
      integer (I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
      do i = 1,6
        do j = 1,3
          ten2mat_r3_dp(i,j) = ten(j,a(i),b(i))
          ten2mat_r3_dp(i,j) = ten(j,b(i),a(i))
        end do
      end do
!
      end function ten2mat_r3_dp
!------------------------------------------------------------------------------
!
!
!
      function ten2mat_r3_dpc(ten)
      use femtypes
      implicit none
      complex (DPC) :: ten(3,3,3), ten2mat_r3_dpc(6,3)
      intent (in) :: ten
!
!
!  Transform a rank 3 tensor (tensor notation) to a 6,3 matrix (matrix notation).
!
!  Internal variables:
      integer (I4B) :: i, j
      integer (I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
      do i = 1,6
        do j = 1,3
          ten2mat_r3_dpc(i,j) = ten(j,a(i),b(i))
          ten2mat_r3_dpc(i,j) = ten(j,b(i),a(i))
        end do
      end do
!
      end function ten2mat_r3_dpc
!------------------------------------------------------------------------------
!
!
!
      function ten2mat_r4_dp_old(ten)
      use femtypes
      implicit none
      real (DP) :: ten(3,3,3,3), ten2mat_r4_dp_old(6,6)
      intent (in) :: ten
!
!
!  Transform a rank 4 tensor (tensor notation) to a 6x6 matrix (matrix notation).
!
!  Internal variables:
      integer (I4B) :: i, j
      integer (I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
       do i = 1,6
         do j = 1,6
           ten2mat_r4_dp_old(i,j) = ten(a(i),b(i),a(j),b(j))
         end do
       end do
!
      end function ten2mat_r4_dp_old
!------------------------------------------------------------------------------
!
!
!
      function ten2mat_r4_dpc_old(ten)
      use femtypes
      implicit none
      complex (DPC) :: ten(3,3,3,3), ten2mat_r4_dpc_old(6,6)
      intent (in) :: ten
!
!
!  Transform a rank 4 tensor (tensor notation) to a 6x6 matrix (matrix notation).
!
!  Internal variables:
      integer (I4B) :: i, j
      integer (I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
       do i = 1,6
         do j = 1,6
           ten2mat_r4_dpc_old(i,j) = ten(a(i),b(i),a(j),b(j))
         end do
       end do
!
      end function ten2mat_r4_dpc_old
!------------------------------------------------------------------------------
