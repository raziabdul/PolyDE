      pure subroutine transformtensor4_DP(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DP) :: a(3,3)
      real(DP) :: intensor(3,3,3,3), outtensor(3,3,3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
!  subroutine to transform a 3D rank-4 tensor such that 
!
!        T'_i,j,k,l = a_i,m * a_j,n * a_k,o * a_l,p * T_m,n,o,p
!
!  with a transformation matrix a
!
!  Input:
!        intensor     rank-4 tensor
!        a            transformation matrix
!
!  Output:
!        outtensor    transformed rank-4 tensor
!
!  local variables
      integer(I4B) :: i, j, k, l, m, n, o
      real(DP) :: ten(3,3,3,3)

      forall (o=1:3, n=1:3, m=1:3, l=1:3)
        ten(l,m,n,o) = a(l,1)*intensor(m,n,o,1) + a(l,2)*intensor(m,n,o,2) + a(l,3)*intensor(m,n,o,3)
      end forall

      forall (n=1:3, m=1:3, l=1:3, k=1:3)
        outtensor(k,l,m,n) = a(k,1)*ten(l,m,n,1) + a(k,2)*ten(l,m,n,2) + a(k,3)*ten(l,m,n,3)
      end forall

      forall (m=1:3, l=1:3, k=1:3, j=1:3)
        ten(j,k,l,m) = a(j,1)*outtensor(k,l,m,1) + a(j,2)*outtensor(k,l,m,2) + a(j,3)*outtensor(k,l,m,3)
      end forall

      forall (l=1:3, k=1:3, j=1:3, i=1:3)
        outtensor(i,j,k,l) = a(i,1)*ten(j,k,l,1) + a(i,2)*ten(j,k,l,2) + a(i,3)*ten(j,k,l,3)
      end forall

      return
      end



      pure subroutine transformtensor3_DP(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DP) :: a(3,3)
      real(DP) :: intensor(3,3,3), outtensor(3,3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
!  subroutine to transform a 3D rank-3 tensor such that 
!
!        T'_i,j,k = a_i,l * a_j,m * a_k,n * T_l,m,n
!
!  with a transformation matrix a
!
!  Input:
!        intensor     rank-3 tensor
!        a            transformation matrix
!
!  Output:
!        outtensor    transformed rank-3 tensor
!
!  local variables
      integer(I4B) :: i, j, k, l, m
      real(DP) :: ten(3,3,3)

      forall (m=1:3, l=1:3, k=1:3)
        outtensor(k,l,m) = a(k,1)*intensor(l,m,1) + a(k,2)*intensor(l,m,2) + a(k,3)*intensor(l,m,3)
      end forall

      forall (l=1:3, k=1:3, j=1:3)
        ten(j,k,l) = a(j,1)*outtensor(k,l,1) + a(j,2)*outtensor(k,l,2) + a(j,3)*outtensor(k,l,3)
      end forall

      forall (k=1:3, j=1:3, i=1:3)
        outtensor(i,j,k) = a(i,1)*ten(j,k,1) + a(i,2)*ten(j,k,2) + a(i,3)*ten(j,k,3)
      end forall

      return
      end



      pure subroutine transformtensor2_DP(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DP) :: a(3,3)
      real(DP) :: intensor(3,3), outtensor(3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
!  subroutine to transform a 3D rank-2 tensor such that 
!
!        T'_i,j = a_i,k * a_j,l * T_k,l
!
!  with a transformation matrix a
!
!  Input:
!        intensor     rank-2 tensor
!        a            transformation matrix
!
!  Output:
!        outtensor    transformed rank-2 tensor
!
!  local variables
      integer(I4B) :: i, j, k
      real(DP) :: ten(3,3)

      forall (k=1:3, j=1:3)
        ten(j,k) = a(j,1)*intensor(k,1) + a(j,2)*intensor(k,2) + a(j,3)*intensor(k,3)
      end forall

      forall (j=1:3, i=1:3)
        outtensor(i,j) = a(i,1)*ten(j,1) + a(i,2)*ten(j,2) + a(i,3)*ten(j,3)
      end forall

      return
      end



      pure subroutine transformtensor1_DP(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DP) :: a(3,3)
      real(DP) :: intensor(3), outtensor(3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
!  subroutine to transform a 3D rank-1 tensor (vector) such that 
!
!        T'_i = a_i,j * T_j
!
!  with a transformation matrix a
!
!  Input:
!        intensor     rank-2 tensor
!        a            transformation matrix
!
!  Output:
!        outtensor    transformed rank-2 tensor
!
!  local variables
      integer(I4B) :: i

      do i=1,3
        outtensor(i) = a(i,1)*intensor(1) + a(i,2)*intensor(2) + a(i,3)*intensor(3)
      end do

      return
      end



      pure subroutine transformtensor4_DPC(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DPC) :: a(3,3)
      complex(DPC) :: intensor(3,3,3,3), outtensor(3,3,3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
!  subroutine to transform a 3D rank-4 tensor such that 
!
!        T'_i,j,k,l = a_i,m * a_j,n * a_k,o * a_l,p * T_m,n,o,p
!
!  with a transformation matrix a
!
!  Input:
!        intensor     rank-4 tensor
!        a            transformation matrix
!
!  Output:
!        outtensor    transformed rank-4 tensor
!
!  local variables
      integer(I4B) :: i, j, k, l, m, n, o
      complex(DPC) :: ten(3,3,3,3)

      forall (o=1:3, n=1:3, m=1:3, l=1:3)
        ten(l,m,n,o) = a(l,1)*intensor(m,n,o,1) + a(l,2)*intensor(m,n,o,2) + a(l,3)*intensor(m,n,o,3)
      end forall

      forall (n=1:3, m=1:3, l=1:3, k=1:3)
        outtensor(k,l,m,n) = a(k,1)*ten(l,m,n,1) + a(k,2)*ten(l,m,n,2) + a(k,3)*ten(l,m,n,3)
      end forall

      forall (m=1:3, l=1:3, k=1:3, j=1:3)
        ten(j,k,l,m) = a(j,1)*outtensor(k,l,m,1) + a(j,2)*outtensor(k,l,m,2) + a(j,3)*outtensor(k,l,m,3)
      end forall

      forall (l=1:3, k=1:3, j=1:3, i=1:3)
        outtensor(i,j,k,l) = a(i,1)*ten(j,k,l,1) + a(i,2)*ten(j,k,l,2) + a(i,3)*ten(j,k,l,3)
      end forall

      return
      end



      pure subroutine transformtensor3_DPC(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DPC) :: a(3,3)
      complex(DPC) :: intensor(3,3,3), outtensor(3,3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
!  subroutine to transform a 3D rank-3 tensor such that 
!
!        T'_i,j,k = a_i,l * a_j,m * a_k,n * T_l,m,n
!
!  with a transformation matrix a
!
!  Input:
!        intensor     rank-3 tensor
!        a            transformation matrix
!
!  Output:
!        outtensor    transformed rank-3 tensor
!
!  local variables
      integer(I4B) :: i, j, k, l, m
      complex(DPC) :: ten(3,3,3)

      forall (m=1:3, l=1:3, k=1:3)
        outtensor(k,l,m) = a(k,1)*intensor(l,m,1) + a(k,2)*intensor(l,m,2) + a(k,3)*intensor(l,m,3)
      end forall

      forall (l=1:3, k=1:3, j=1:3)
        ten(j,k,l) = a(j,1)*outtensor(k,l,1) + a(j,2)*outtensor(k,l,2) + a(j,3)*outtensor(k,l,3)
      end forall

      forall (k=1:3, j=1:3, i=1:3)
        outtensor(i,j,k) = a(i,1)*ten(j,k,1) + a(i,2)*ten(j,k,2) + a(i,3)*ten(j,k,3)
      end forall

      return
      end



      pure subroutine transformtensor2_DPC(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DPC) :: a(3,3)
      complex(DPC) :: intensor(3,3), outtensor(3,3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
!  subroutine to transform a 3D rank-2 tensor such that 
!
!        T'_i,j = a_i,k * a_j,l * T_k,l
!
!  with a transformation matrix a
!
!  Input:
!        intensor     rank-2 tensor
!        a            transformation matrix
!
!  Output:
!        outtensor    transformed rank-2 tensor
!
!  local variables
      integer(I4B) :: i, j, k
      complex(DPC) :: ten(3,3)

      forall (k=1:3, j=1:3)
        ten(j,k) = a(j,1)*intensor(k,1) + a(j,2)*intensor(k,2) + a(j,3)*intensor(k,3)
      end forall

      forall (j=1:3, i=1:3)
        outtensor(i,j) = a(i,1)*ten(j,1) + a(i,2)*ten(j,2) + a(i,3)*ten(j,3)
      end forall

      return
      end



      pure subroutine transformtensor1_DPC(a, intensor, outtensor)
      use femtypes
      implicit none
      real(DPC) :: a(3,3)
      complex(DPC) :: intensor(3), outtensor(3)
      intent (in) :: a, intensor
      intent (out) :: outtensor
!  subroutine to transform a 3D rank-1 tensor (vector) such that 
!
!        T'_i = a_i,j * T_j
!
!  with a transformation matrix a
!
!  Input:
!        intensor     rank-2 tensor
!        a            transformation matrix
!
!  Output:
!        outtensor    transformed rank-2 tensor
!
!  local variables
      integer(I4B) :: i

      do i=1,3
        outtensor(i) = a(i,1)*intensor(1) + a(i,2)*intensor(2) + a(i,3)*intensor(3)
      end do

      return
      end



      pure subroutine mat2ten_r4_DP(matrix, tensor, factors)
      use femtypes
      implicit none
      real(DP) :: matrix(6,6), tensor(3,3,3,3), factors(3)
      intent (in) :: matrix, factors
      intent (out) :: tensor
!
!  rank-4 tensor conversion from two-indices matrix notation to 4-indices tensor notation
!  i.e. convert the 6x6 matrix (matrix notation) to rank 4 tensor (tensor notation)
!  we use Voigt convention for indices 1,1 -> 1;   2,2 -> 2;   3,3 -> 3;   2,3 -> 4,   1,3 -> 5;   1,2 -> 6
!
!  depending on the definition, factors 2 or 4 have to be introduced during conversion
!       stiffness tensor             factors are (1., 1., 1.)
!       compliance tensor            factors are (1., 2., 4.)
!       piezoresistivity tensor      factors are (1., 1., 2.)
!
!  CAUTION: the tensor must be symmetric
!
!  Internal variables:
      integer(I4B) :: i, j
      integer(I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
      do j = 1,3
        do i = 1,3
          tensor(a(i),b(i),a(j),b(j)) = matrix(i,j) / factors(1)
        end do
      end do

      do j = 1,3
        do i = 4,6
          tensor(a(i),b(i),a(j),b(j)) = matrix(i,j) / factors(2)
          tensor(b(i),a(i),a(j),b(j)) = matrix(i,j) / factors(2)
        end do
      end do

      do j = 4,6
        do i = 1,3
          tensor(a(i),b(i),a(j),b(j)) = matrix(i,j) / factors(2)
          tensor(a(i),b(i),b(j),a(j)) = matrix(i,j) / factors(2)
        end do
      end do

      do j = 4,6
        do i = 4,6
          tensor(a(i),b(i),a(j),b(j)) = matrix(i,j) / factors(3)
          tensor(a(i),b(i),b(j),a(j)) = matrix(i,j) / factors(3)
          tensor(b(i),a(i),b(j),a(j)) = matrix(i,j) / factors(3)
          tensor(b(i),a(i),a(j),b(j)) = matrix(i,j) / factors(3)
        end do
      end do

      end subroutine mat2ten_r4_DP



      pure subroutine mat2ten_r4_DPC(matrix, tensor, factors)
      use femtypes
      implicit none
      real(DP) :: factors(3)
      complex(DPC) :: matrix(6,6), tensor(3,3,3,3)
      intent (in) :: matrix, factors
      intent (out) :: tensor
!
!  rank-4 tensor conversion from two-indices matrix notation to 4-indices tensor notation
!  i.e. convert the 6x6 matrix (matrix notation) to rank 4 tensor (tensor notation)
!  we use Voigt convention for indices 1,1 -> 1;   2,2 -> 2;   3,3 -> 3;   2,3 -> 4,   1,3 -> 5;   1,2 -> 6
!
!  depending on the definition, factors 2 or 4 have to be introduced during conversion
!       stiffness tensor             factors are (1., 1., 1.)
!       compliance tensor            factors are (1., 2., 4.)
!       piezoresistivity tensor      factors are (1., 1., 2.)
!
!  CAUTION: the tensor must be symmetric
!
!  Internal variables:
      integer(I4B) :: i, j
      integer(I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
      do j = 1,3
        do i = 1,3
          tensor(a(i),b(i),a(j),b(j)) = matrix(i,j) / factors(1)
        end do
      end do

      do j = 1,3
        do i = 4,6
          tensor(a(i),b(i),a(j),b(j)) = matrix(i,j) / factors(2)
          tensor(b(i),a(i),a(j),b(j)) = matrix(i,j) / factors(2)
        end do
      end do

      do j = 4,6
        do i = 1,3
          tensor(a(i),b(i),a(j),b(j)) = matrix(i,j) / factors(2)
          tensor(a(i),b(i),b(j),a(j)) = matrix(i,j) / factors(2)
        end do
      end do

      do j = 4,6
        do i = 4,6
          tensor(a(i),b(i),a(j),b(j)) = matrix(i,j) / factors(3)
          tensor(a(i),b(i),b(j),a(j)) = matrix(i,j) / factors(3)
          tensor(b(i),a(i),b(j),a(j)) = matrix(i,j) / factors(3)
          tensor(b(i),a(i),a(j),b(j)) = matrix(i,j) / factors(3)
        end do
      end do

      end subroutine mat2ten_r4_DPC



      pure subroutine ten2mat_r4_DP(matrix, tensor, factors)
      use femtypes
      implicit none
      real(DP) :: matrix(6,6), tensor(3,3,3,3), factors(3)
      intent (in) :: tensor, factors
      intent (out) :: matrix
!
!  rank-4 tensor conversion from 4-indices tensor notation to two-indices matrix notation
!  i.e. convert rank 4 tensor (tensor notation) to the 6x6 matrix (matrix notation)
!  we use Voigt convention for indices 1,1 -> 1;   2,2 -> 2;   3,3 -> 3;   2,3 -> 4,   1,3 -> 5;   1,2 -> 6
!
!  depending on the definition, factors 2 or 4 have to be introduced during conversion
!       stiffness tensor             factors are (1., 1., 1.)
!       compliance tensor            factors are (1., 2., 4.)
!       piezoresistivity tensor      factors are (1., 1., 2.)
!
!  CAUTION: the tensor must be symmetric
!
!  Internal variables:
      integer(I4B) :: i, j
      integer(I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
      do j = 1,3
        do i = 1,3
          matrix(i,j) = tensor(a(i),b(i),a(j),b(j)) * factors(1)
        end do
      end do

      do j = 1,3
        do i = 4,6
          matrix(i,j) = tensor(a(i),b(i),a(j),b(j)) * factors(2)
        end do
      end do

      do j = 4,6
        do i = 1,3
          matrix(i,j) = tensor(a(i),b(i),a(j),b(j)) * factors(2)
        end do
      end do

      do j = 4,6
        do i = 4,6
          matrix(i,j) = tensor(a(i),b(i),a(j),b(j)) * factors(3)
        end do
      end do

      end subroutine ten2mat_r4_DP



      pure subroutine ten2mat_r4_DPC(matrix, tensor, factors)
      use femtypes
      implicit none
      real(DPC) :: factors(3)
      complex(DPC) :: matrix(6,6), tensor(3,3,3,3)
      intent (in) :: tensor, factors
      intent (out) :: matrix
!
!  rank-4 tensor conversion from 4-indices tensor notation to two-indices matrix notation
!  i.e. convert rank 4 tensor (tensor notation) to the 6x6 matrix (matrix notation)
!  we use Voigt convention for indices 1,1 -> 1;   2,2 -> 2;   3,3 -> 3;   2,3 -> 4,   1,3 -> 5;   1,2 -> 6
!
!  depending on the definition, factors 2 or 4 have to be introduced during conversion
!       stiffness tensor             factors are (1., 1., 1.)
!       compliance tensor            factors are (1., 2., 4.)
!       piezoresistivity tensor      factors are (1., 1., 2.)
!
!  CAUTION: the tensor must be symmetric
!
!  Internal variables:
      integer(I4B) :: i, j
      integer(I4B), parameter :: a(6)=(/1,2,3,2,3,1/), b(6)=(/1,2,3,3,1,2/)
!
!
      do j = 1,3
        do i = 1,3
          matrix(i,j) = tensor(a(i),b(i),a(j),b(j)) * factors(1)
        end do
      end do

      do j = 1,3
        do i = 4,6
          matrix(i,j) = tensor(a(i),b(i),a(j),b(j)) * factors(2)
        end do
      end do

      do j = 4,6
        do i = 1,3
          matrix(i,j) = tensor(a(i),b(i),a(j),b(j)) * factors(2)
        end do
      end do

      do j = 4,6
        do i = 4,6
          matrix(i,j) = tensor(a(i),b(i),a(j),b(j)) * factors(3)
        end do
      end do

      end subroutine ten2mat_r4_DPC


      function trans_matrix(x1,x2)
      use femtypes
      use feminterface, only : cross_product
      implicit none
      real (DP) :: x1(3), x2(3), trans_matrix(3,3)
      intent (in) :: x1, x2
!  x1  vector of the first crystal axis 
!  x2  vector of the second crystal axis
!
!  Comments:
!      x1, x2 are to be given in given in the coordinate system (x,y,z)
!      with crystal axes system we here mean the (orthogonal) system of the original tensor data
!
!      A vector V is transfoms as:
!        V'  = a * V
!      A tensor T then transforms as:
!        T'_i,j,k,l = a_i,m * a_j,n * a_k,o * a_l,p * T_m,n,o,p
!
!      For back-transformation, the transposed of transformation matrix a is to be used
!        V  = tanspose(a) * V'
!
!  Internal variables:
      real (DP) :: h1(3), h2(3), h3(3), h1norm, h3norm
!
!  renormalize first vector
      h1norm = sqrt(sum(x1**2))
      h1 = x1/h1norm
!  calculate third axis by the corss product between first and second axis
      h3 = cross_product(h1,x2)
      if (abs(h3(1))+abs(h3(2))+abs(h3(3)) .le. tiny(1._DP)) then
        print '(a)', ' '
        print '(a)', 'TRANS_MATRIX - Fatal error!'
        print '(a)', '  Given vectors are not collinear'
        pause
      end if
!  renormalize third vector
      h3norm = sqrt(sum(h3**2))
      h3 = h3/h3norm
!  now h1, h3 are orthonormal
!  reorthogonalize second vector
      h2 = cross_product(h3,h1)
!
      trans_matrix(1,1) = h1(1)
      trans_matrix(2,1) = h1(2)
      trans_matrix(3,1) = h1(3)
      trans_matrix(1,2) = h2(1)
      trans_matrix(2,2) = h2(2)
      trans_matrix(3,2) = h2(3)
      trans_matrix(1,3) = h3(1)
      trans_matrix(2,3) = h3(2)
      trans_matrix(3,3) = h3(3)

      return
      end

