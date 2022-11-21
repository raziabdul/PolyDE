      subroutine lusolver_c(a,b)
      use femtypes
      use feminterface, only: ludcmp, lubksb
      complex(DPC) :: a(:,:), b(:)
      intent(inout) :: a, b
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
!    $Revision: 1.7 $
!    $Date: 2014/08/27 15:22:21 $
!    $Author: m_kasper $
!
!  Solving of a (small and dense) equation system by LU decomposition
!  adapted form numerical recipes
!
      integer(I4B) indx(size(a,1))
      call ludcmp(a,indx)
      call lubksb(a,indx,b)
      end subroutine lusolver_c
!
!
!
      subroutine ludcmp_c(a,indx)
      use femtypes
      use feminterface, only: swap, outerprod
      implicit none
      complex(DPC) a(:,:)
      integer(I4B) indx(:)
      intent(out) :: indx
      intent(inout) :: a
!  This routine replaces the matrix a by the LU decomposition
!  of a rowwise permutation of itself. 
!  indx is an output vector of length N that records the row permutation effected 
!  by the partial pivoting;
!  This routine is used in combination with lubksb to solve linear equations 
!  or invert a matrix.
      real(DP) vv(size(a,1)) 
!  vv stores the implicit scaling of each row.
      real(DP), parameter :: tinny=tiny(1._DP)
      integer(I4B) j, n, imax, imaxl(1)
!
      n=size(a,1)
!  Loop over rows to get the implicit scaling
      vv=maxval(abs(a),dim=2)
      if (any(vv .eq. 0.0)) then
          print*,'singular matrix in ludcmp'
      end if
      vv=1.0_DP/vv
      do j=1,n
!  Find the pivot row
        imaxl=maxloc(vv(j:n)*abs(a(j:n,j)))
        imax=(j-1)+imaxl(1)
        if (j .ne. imax) then
          call swap(a(imax,:),a(j,:))
          vv(imax)=vv(j)
        end if
        indx(j)=imax
        if (abs(a(j,j)) .eq. 0._DP) then
          a(j,j)=tinny
        end if
!  If the pivot element is zero the matrix is singular (at least to the 
!  precision of the algorithm).For some applications on singular matrices,
!  it is desirable to substitute TINY for zero.
        a(j+1:n,j)=a(j+1:n,j)/a(j,j)
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
      end do
      return
      end subroutine ludcmp_c
!
!
!
      pure subroutine lubksb_c(a,indx,b)
      use femtypes
      use feminterface, only:
      implicit none
      complex(DPC) a(:,:), b(:)
      integer(I4B) indx(:)
      intent(in) :: a, indx
      intent(inout) :: b
!  Solves the set linear equations AX = B. 
!  a stores the LU decomposition determined by the routine ludcmp. 
!  indx is input as the permutation vector returned by ludcmp
!  b is input as the right-hand-side vector B and returns with the solution 
!  vector X
      integer(I4B) :: i,n,ii,ll
      complex(DPC) :: summ
!
      n=size(a,1)
      ii=0
      do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (ii .ne. 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
        else if (summ .ne. 0.0) then
          ii=i
        end if
        b(i)=summ
      end do
!  Now we do the backsubstitution
      do i=n,1,-1
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
      end do
      return
      end subroutine lubksb_c
!
!
!
      pure function outerprod_c(a,b)
      use femtypes
      implicit none
      complex(DPC) a(:), b(:)
      complex(DPC) outerprod_c(size(a),size(b))
      intent(in) :: a, b
!
      outerprod_c = spread(a,dim=2,ncopies=size(b)) *                     &
     &   spread(b,dim=1,ncopies=size(a))
      return
      end function outerprod_c
!
!
!
      pure subroutine swap_c(a,b)
      use femtypes
      implicit none
      complex(DPC) a(:), b(:)
      intent(inout) :: a, b
!
      complex(DPC) dum(size(a))
      dum=a
      a=b
      b=dum
      return
      end subroutine swap_c
!
!
!
      subroutine lusolver_d(a,b)
      use femtypes
      use feminterface, only: ludcmp, lubksb
      implicit none
      real(DP) :: a(:,:), b(:) 
      intent(inout) :: a, b
!*************************************************************************
!  Solving of a (small and dense) equation system by LU decomposition
!  adapted form numerical recipes
!*************************************************************************
      integer(I4B) :: indx(size(a,1))
      call ludcmp(a,indx)
      call lubksb(a,indx,b)
      end subroutine lusolver_d
!
!
!
      subroutine ludcmp_d(a,indx)
      use femtypes
      use feminterface, only: swap, outerprod
      implicit none
      real(DP) :: a(:,:)
      integer(I4B) indx(:)
      intent(out) :: indx
      intent(inout) :: a
!*************************************************************************
!  This routine replaces the matrix a by the LU decomposition
!  of a rowwise permutation of itself. 
!  indx is an output vector of length N that records the row permutation effected 
!  by the partial pivoting;
!  This routine is used in combination with lubksb to solve linear equations 
!  or invert a matrix.
!*************************************************************************
! local variables
!*************************************************************************
      real(DP) :: vv(size(a,1)) 
!*************************************************************************
!  vv stores the implicit scaling of each row.
      real(DP), parameter :: tinny=tiny(1._DP)
      integer(I4B) j, n, imax, imaxl(1)
!*************************************************************************
      n=size(a,1)
!  Loop over rows to get the implicit scaling
      vv=maxval(abs(a),dim=2)
      if (any(vv .eq. 0.0)) then 
          print*,'singular matrix in ludcmp'
      end if
      vv=1.0_DP/vv
      do j=1,n
!  Find the pivot row
         imaxl=maxloc(vv(j:n)*abs(a(j:n,j)))
         imax=(j-1)+imaxl(1)
       
         if (j .ne. imax) then
            call swap(a(imax,:),a(j,:))
            vv(imax)=vv(j)        
         end if
        
         indx(j)=imax        
         if (abs(a(j,j)) .eq. 0._DP) then
            a(j,j)=tinny
         end if
!  if the pivot element is zero the matrix is singular (at least to the 
!  femtypes of the algorithm).For some applications on singular matrices,
!  it is desirable to substitute tiny for zero.
         a(j+1:n,j)=a(j+1:n,j)/a(j,j)
         a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
      end do
      return
      end subroutine ludcmp_d
!
!
!
      pure subroutine lubksb_d(a,indx,b)
      use femtypes
      use feminterface, only:
      implicit none      
      real(DP) :: a(:,:), b(:)
      integer(I4B) :: indx(:)
      intent(in) :: a, indx
      intent(inout) :: b
!*************************************************************************
!  Solves the set linear equations AX = B. 
!  a stores the LU decomposition determined by the routine ludcmp. 
!  indx is input as the permutation vector returned by ludcmp
!  b is input as the right-hand-side vector B and returns with the solution 
!  vector X
!*************************************************************************
      integer(I4B) :: i,n,ii,ll
      real(DP) :: summ
!*************************************************************************
      n=size(a,1)
      ii=0
      do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (ii .ne. 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
        else if (summ .ne. 0.0) then
          ii=i
        end if
        b(i)=summ
      end do
!  Now we do the backsubstitution
      do i=n,1,-1
        b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
      end do
      return
      end subroutine lubksb_d
!
!
!
      pure function outerprod_d(a,b)
      use femtypes
      implicit none
      real(DP) :: a(:), b(:)
      real(DP) :: outerprod_d(size(a),size(b))
      intent(in) :: a, b

      outerprod_d = SPREAD(a,dim=2,ncopies=size(b))*SPREAD(b,dim=1,ncopies=size(a))

      return
      end function outerprod_d
!
!
!
      pure subroutine swap_d(a,b)
      use femtypes
      implicit none
      real(DP) :: a(:), b(:)
      intent(inout) :: a, b
!*************************************************************************
! local variables
!*************************************************************************
      real(DP) :: dum(size(a))
!*************************************************************************
      dum = a
      a = b
      b = dum
      return
      end subroutine swap_d
