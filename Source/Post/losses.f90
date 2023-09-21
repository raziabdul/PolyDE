      function point_loss(ielem,lambda)
      use feminterface, only: field, fieldquantity
      use femtypes
      use globalvariables
      use matconstants
      implicit none
      integer (I4B) :: ielem
      real (DP) :: lambda(3), point_loss
      intent (in) :: ielem, lambda
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
!    $Revision: 1.15 $
!    $Date: 2014/07/15 13:06:21 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Calculates the Electromagnetic loss at a given point
!
!  Input:
!     ielem         element number
!     lambda        barycentric coordinates
!
!  Output:
!     point_loss    electromagnetic loss
!
!  local variables:
      complex (DPC) :: z(15,1), temp, tempvec(2)
      logical :: typ(5), ok1
      character (len=10) :: unit
      character (len=50) :: descriptor
!------------------------------------------------------------------------------
!      
!  get the field and the curl of field
      typ = (/.true.,.true.,.false.,.true.,.false./)
      call field(ielem,lambda,typ,z)
!
!  calculate energies depending on physics mode
      select case (physics)
      case ('TEWAVE')
        call fieldquantity(ielem,'HX        ',lambda,0._DP,tempvec(1),descriptor,unit,ok1)
        call fieldquantity(ielem,'HY        ',lambda,0._DP,tempvec(2),descriptor,unit,ok1)
        point_loss = -0.5_DP * real( dot_product(tempvec,(/z(5,1),-z(4,1)/)) - &
                     z(1,1)*conjg(cmplx(0._DP,1._DP,DPC)*(-z(2,1) + z(10,1)))/(omega*mu0) ,DP)
!
      case ('TMWAVE')
        call fieldquantity(ielem,'WM        ',lambda,0._DP,temp,descriptor,unit,ok1)
        call fieldquantity(ielem,'EX        ',lambda,0._DP,tempvec(1),descriptor,unit,ok1)
        call fieldquantity(ielem,'EY        ',lambda,0._DP,tempvec(2),descriptor,unit,ok1)
        point_loss = -0.5_DP * real( cmplx(0._DP,-4*omega,DPC)*temp - &
                     dot_product((/z(5,1),-z(4,1)/),tempvec) ,DP)
      end select
!
      return
      end function point_loss
!
!
!
      subroutine total_loss(region,loss)
      use feminterface, only: flaech, get2Dintegpoints, point_loss, low2hi
      use femtypes
      use globalvariables
      use matconstants
      implicit none
      character (len=14), intent(in) :: region
      real (DP), intent (out) :: loss
!
!----------------------------------------------------------------------------
!
!  Calculates the total loss by integrating the loss for each element
!  in the region
!
!  Output:
!     total_loss      total loss in the region
! 
!  local variables:
      integer (I4B) :: elem, intorder, npkt, errcodeint, k, matindex, pos
      real (DP) :: area_e
      real (DP), allocatable :: elemloss(:)
      real (DP), allocatable :: weight(:), lambda(:,:)
      character (len=25) :: mat(maxmat), reg
!----------------------------------------------------------------------------
!      
!  initialize losses
      allocate(elemloss(n))
      elemloss = 0._DP
      loss = 0._DP
!  copy 14 string to 25 string
      reg = region
!  copy matnames to mat and convert to upper case
      mat = matnames
      do k =1,maxmat
        call low2hi(mat(k),25)
        if (reg .eq. mat(k)) then
          pos = k
        end if
      end do
!
!  LOOP OVER ALL ELEMENTS:
      do elem = 1,n
!  go to next element if region does not match element's material
        matindex = matzif(geb(elem))
        if (pos .ne. matindex) cycle
!
!  determine the order of numerical integration (nature=1)
        intorder = 2*ep(elem,1)
!  fetch numerical integration points (Gauss points)
        call get2Dintegpoints(intorder,npkt,weight,lambda,errcodeint)
!
!
!  LOOP OVER ALL INTPOINTS:
        do k = 1,npkt
!  compute the loss for each element by summing up over integration points
          elemloss(elem) = elemloss(elem) + weight(k)*point_loss(elem,lambda(:,k))
        end do
!
!  calculate the area of the element
        area_e = flaech(xn(e(1,elem)),yn(e(1,elem)),xn(e(2,elem)),yn(e(2,elem)),&
                        xn(e(3,elem)),yn(e(3,elem)))
!
        elemloss(elem) = area_e * elemloss(elem)
        deallocate(weight,lambda)
      end do
!  calculate total loss for the region by summing up over all elements;
      loss = sum(elemloss)
!
      deallocate(elemloss)
      return
      end subroutine total_loss