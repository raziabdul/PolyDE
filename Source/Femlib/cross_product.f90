      pure function cross_product_i(a,b)
      use femtypes
      implicit none
      integer (I4B) :: a(3), b(3), cross_product_i(3)
      intent (in) :: a, b
!
!    $Revision: 1.5 $
!    $Date: 2014/08/22 10:42:05 $
!    $Author: m_kasper $
!
!  Compute cross produt of two integer vectors a, b
!
!  local variables 

      cross_product_i(1)= a(2)*b(3)-a(3)*b(2)
      cross_product_i(2)= a(3)*b(1)-a(1)*b(3)
      cross_product_i(3)= a(1)*b(2)-a(2)*b(1)

      return
      end function cross_product_i
!
!
!
      pure function cross_product_r(a,b)
      use femtypes
      implicit none
      real (DP) :: a(3), b(3), cross_product_r(3)
      intent (in) :: a, b
!
!  Compute cross product of two real vectors a, b
!
!  local variables 

      cross_product_r(1)= a(2)*b(3)-a(3)*b(2)
      cross_product_r(2)= a(3)*b(1)-a(1)*b(3)
      cross_product_r(3)= a(1)*b(2)-a(2)*b(1)

      return
      end function cross_product_r
!
!
!
      pure function cross_product_dpc(a,b)
      use femtypes
      implicit none
      complex (DPC) :: a(3), b(3), cross_product_dpc(3)
      intent (in) :: a, b
!
!  Compute cross produt of two complex vectors a, b.
!
!  local variables 

      cross_product_dpc(1)= a(2)*b(3)-a(3)*b(2)
      cross_product_dpc(2)= a(3)*b(1)-a(1)*b(3)
      cross_product_dpc(3)= a(1)*b(2)-a(2)*b(1)

      return
      end function cross_product_dpc
!
!
!
      pure function cross_product_rdpc(a,b)
      use femtypes
      implicit none
      complex (DPC) :: b(3), cross_product_rdpc(3)
      real (DP)     :: a(3)
      intent (in) :: a, b
!
!  Compute cross produt of two complex vectors a, b.
!
!  local variables 

      cross_product_rdpc(1)= a(2)*b(3)-a(3)*b(2)
      cross_product_rdpc(2)= a(3)*b(1)-a(1)*b(3)
      cross_product_rdpc(3)= a(1)*b(2)-a(2)*b(1)

      return
      end function cross_product_rdpc
!
!
!
      pure function cross_product_dpcr(a,b)
      use femtypes
      implicit none
      complex (DPC) :: a(3), cross_product_dpcr(3)
      real (DP)     :: b(3)
      intent (in) :: a, b
!
!
!  Compute cross product of complex vectors a and real vector b.
!
!  local variables 

      cross_product_dpcr(1)= a(2)*b(3)-a(3)*b(2)
      cross_product_dpcr(2)= a(3)*b(1)-a(1)*b(3)
      cross_product_dpcr(3)= a(1)*b(2)-a(2)*b(1)

      return
      end function cross_product_dpcr
