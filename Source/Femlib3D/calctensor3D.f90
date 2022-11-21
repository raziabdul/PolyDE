      subroutine calctensor3D(tensor,vect1,vect2)
      use feminterface, only: cross_product
      use femtypes
      implicit none
      real (DP) :: vect1(3), vect2(3)
      complex (DPC) :: tensor(3,3)
      intent(in) :: vect1, vect2
      intent(inout) :: tensor
!
!------------------------------------------------------------------------------
!    $Revision: 1.4 $
!    $Date: 2014/07/15 13:00:01 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!* general 3D coordinate transformation of rank two tensor *
!
!  In this program, the coordinate transformation of rank 2 tensor has been done
!  from principal  crystal coordinate system to reference coordinate system. The
!  input arguments have to be vectors in the direction of first and second prin-
!  cipal crystal axes in terms of reference coordinate system  
!
!  input:
!        vect1    first principal crystal axis in terms of reference coordinate system 
!        vect2    second principal crystal axis in terms of reference coordinate system
!
!  in-/output:
!        tensor   rank two tensor to be coordinate-transformed 

! internal variables 
      real (DP) :: i, sum, sum1
      real (DP) :: vectn1(3), vectn2(3), vectn3(3), tr_mat(3,3)
      complex (DPC) :: r2ten(3,3), r2tran(3,3), sumc
      integer (I4B) :: ix1, ix2, i1, i2
!------------------------------------------------------------------------------
!
      sum = 0._DP
      sum1 = 0._DP
!  copy tensor to internal variable r2ten
      r2ten = tensor
!
      do i = 1,3
        sum  = sum  + (vect1(i))**2._DP
        sum1 = sum1 + (vect2(i))**2._DP
      end do
!
!  vectn1 : x' directed unit vector
!  vectn2 : y' directed unit vector
      sum  = sqrt(sum)
      sum1 = sqrt(sum1)
      vectn1 = vect1/sum 
      vectn2 = vect2/sum1
!
!  cross product tan_vect1  =  nor_vect X tan_vect
      vectn3 = cross_product(vectn1,vectn2)
! 
!  Coordinate Transformation Matrix 
      tr_mat(1,:) = vectn1 
      tr_mat(2,:) = vectn2 
      tr_mat(3,:) = vectn3 
!
!  coordinate transformation from X',Y',Z'-> X,Y,Z
      tr_mat = transpose(tr_mat)
      do ix1 = 1,3
        do ix2 = 1,3
          sumc = 0 
          do i1 = 1,3 
            do i2 = 1,3
              sumc = sumc + tr_mat(ix1,i1)*tr_mat(ix2,i2)*r2ten(i1,i2)
            end do 
          end do 
          r2tran(ix1,ix2) = sumc ! with respect to the reference coordinate system (x,y,z) 
        end do 
      end do 
!
!  copy transformed tensor r2tran to tensor
      tensor = r2tran
!
      end subroutine calctensor3D