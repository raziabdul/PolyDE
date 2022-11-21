subroutine streamfunction(ia,ja)
use femtypes
use feminterface
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer(I4B), pointer :: ia(:), ja(:)
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
!    $Revision: 1.6 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!***********************************************************************************************
!
!        STREAM FUNCTION
!
!***********************************************************************************************
real(DP), pointer :: lower(:), upper(:), diag(:), rhs(:)
real(DP) :: eps, epsgl, resgl
character(len=20) :: matrixtype
logical :: symmetric
!***********************************************************************************************
allocate(lower(1:ia(ndof)),diag(ndof),rhs(ndof),phi(ndof))

call getsetting('MATRIXTYPE',matrixtype)
if (matrixtype .eq. 'SYMMETRIC') then
   symmetric=.true.
else
   symmetric=.false.
end if

! assembly
call assembly_streamfn(lower,upper,diag,rhs,ia,ja,symmetric)

call getsetting('LINSOLVER_ERROR',eps)
! solve
call solve2(lower,lower,diag,rhs,phi,ndof,eps,ia,ja,epsgl,resgl,symmetric)

deallocate(lower,diag,rhs)

end subroutine streamfunction
!###############################################################################################
!###############################################################################################
!###############################################################################################




subroutine assembly_streamfn(lower,upper,diag,rhs,ia,ja,symmetric)
use femtypes
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer(I4B), pointer :: ia(:), ja(:)
real(DP), pointer :: lower(:), upper(:), diag(:), rhs(:)
logical :: symmetric
!***********************************************************************************************
!
! ASSEMBLY FOR STREAM FUNCTION
!
!***********************************************************************************************
integer(I4B) :: i, k, l, kn, knn, ln, lnn, k1, l1, mb
integer(I4B) :: nff, vzk, vzl
real(DP), pointer :: ai(:,:), bi(:)
!***********************************************************************************************
allocate(lower(1:ia(ndof)), upper(1:ia(ndof)), diag(ndof), rhs(ndof))
!  clear matrix, right-hand-side
lower= 0._DP
upper= 0._DP
diag = 0._DP
rhs  = 0._DP

do i=1,n
!  compute element matrix
   call elementmatrix_streamfn(i,ep(i,1),ai,bi,nff)
   
   do k=1,nff
!  row index k: local;  kn: global
      kn=abs(eg(i,1)%d(k))
      
!  if this dof is not present
      if (kn .eq. 0) cycle
      vzk=sign(1,eg(i,1)%d(k))
!  update rhs
      rhs(kn)=rhs(kn)+bi(k)*vzk
!  diagonal
      diag(kn)=diag(kn)+ai(k,k)
!  column index l: local;  ln: global
loop4:do l=1,k-1
         ln=abs(eg(i,1)%d(l))
         if (ln .eq. 0) cycle
         vzl=sign(1,eg(i,1)%d(l))

         if (kn .lt. ln) then
!  swap row and column -  column index must be smaller than row index
            lnn=kn
            knn=ln
            k1=k
            l1=l
         else
            lnn=ln
            knn=kn
            k1=l
            l1=k
         end if

         do mb=ia(knn-1)+1,ia(knn)
         !  search for the the column in ja in the list
            if (ja(mb) .eq. lnn) then
            !  found
               lower(mb)=lower(mb)+ai(l1,k1)*vzk*vzl
!  if the matrix is known to be symmetric we do not need to assign upper matrix
               if (.not.symmetric) upper(mb)=upper(mb)+ai(k1,l1)*vzk*vzl
               cycle loop4
            end if
         end do
         print*,'***Error in Matrix Assembly',' kn: ',kn,' ln: ',ln
      end do loop4
   end do
   deallocate(ai, bi)
end do
end subroutine assembly_streamfn
!###############################################################################################
!###############################################################################################
!###############################################################################################




subroutine elementmatrix_streamfn(elem,polyorder,a,b,nff)
use femtypes
use feminterface
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer (I4B) :: elem
integer (I4B) :: polyorder, nff
real(DP), pointer :: a(:,:), b(:)
intent (in) :: elem, polyorder
intent (out) :: nff
!***********************************************************************************************
!     elementmatrix of stream function
!***********************************************************************************************
integer(I4B) :: i,j,k,ii,jj,ivertex, npkt
integer(I4B) :: ivar
integer(I4B) :: intorder, polyhi, polylo
integer(I4B) :: err, errcode
integer(I4B), dimension(3), parameter :: inach=(/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg=(/3,1,2/)
real(DP) :: aw
real(DP), allocatable :: weight(:), lambda(:,:), xsi(:), gxsi(:,:)
real(DP), pointer :: elvel(:,:), velterm(:)
logical, pointer :: dbc(:)
!***********************************************************************************************
polyhi    = polyorder
polylo    = 1
!  size of the element matrix
nff=(polyorder+1)*(polyorder+2)/2

allocate(elvel(2,nff), velterm(nff) )
allocate(a(nff,nff), b(nff), xsi(nff), gxsi(nff,2) )
a = 0.0_DP
b = 0.0_DP
velterm = 0.0_DP

do j = 1,nff
   jj = eg(elem,1)%d(j)
   do ivar = 2,3
      elvel(ivar-1,j) = theta(1)*unkno(ivar,jj)
   end do ! ivar
end do ! j

intorder=2*polyorder

!  fetch numerical integration points (Gauss points)
call get2Dintegpoints(intorder, npkt, weight, lambda, err)

do k=1,npkt
!  get shape function and gradients at location of integration points
   call shapefunction(lambda(:,k),xn(e(:,elem)),yn(e(:,elem)),polylo,polyhi, &
                      nff,.true.,xsi,gxsi,errcode)

   aw=areaf(elem)*weight(k)

   do j=1,nff

      velterm(j) = ( gxsi(j,2)*elvel(1,j) - gxsi(j,1)*elvel(2,j) )

      do i=1,nff
         a(j,i) = a(j,i) + aw*(gxsi(j,1)*gxsi(i,1) + gxsi(j,2)*gxsi(i,2))
         b(i) = b(i) - aw*xsi(i)*velterm(j)
      end do
   end do
end do
deallocate(weight, lambda, gxsi, xsi)
deallocate(elvel,velterm)
!*****************************************************************************************
!  fix starting point of stream function to be at keypoint number 1
!*****************************************************************************************
allocate(dbc(nff))
dbc(1:nff) =.false.

do ivertex = 1,3
   ii = eg(elem,1)%d(ivertex)
   ! keypoint 1 is set to be a starting point for the streamline
   if (-(kzi(ii)) .NE. 1) cycle
   b(ivertex)  = 0.0_DP
   dbc(ivertex)= .true.
end do ! ivertex
!*****************************************************************************************
do i=1,nff
   if (dbc(i)) then
      a(i,1:nff)=0._DP
      a(i,i)=1._DP
   else 
      do j=1,nff
         if (dbc(j)) then
            b(i)=b(i)-a(i,j)*b(j)
            a(i,j)=0._DP
         end if
      end do
   end if
end do
deallocate(dbc)

end subroutine elementmatrix_streamfn