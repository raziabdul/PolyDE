subroutine pstiff(lower,diag,ia,ja)
use femtypes
use feminterface, only: reallocate, destroyarrptr, getsetting
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer (I4B), pointer :: ia(:), ja(:)
real(DP), pointer :: lower(:), diag(:)
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
!    $Date: 2015/04/01 10:51:16 $
!    $Author: juryanatzki $
!
!*****************************************************************************************
!
!  Compute ony stiffness matrix for step 2
!
!remark: modified from assembly.f90 in solver
!*****************************************************************************************
! local variables
!*****************************************************************************************
integer(I4B) :: i, j, k, l, kn, knn, ln, lnn, k1, l1, mb, acsize
integer(I4B) :: nff, vzk, vzl
integer(I4B), parameter :: isize=20
real(DP), pointer :: ai(:,:)
real(DP), pointer :: upper(:)
type (ARRPTRI), pointer :: tja(:)
!character(len=25) :: matrixtype
logical :: lstat, symmetric
!*****************************************************************************************
call getsetting('PHYSICS_MODE',physics)
if (physics .eq. 'FLUID') symmetric =.true.


!call getsetting('MATRIXTYPE',matrixtype)
!if (matrixtype .eq. 'SYMMETRIC') then
!   symmetric=.true.
!else
!   symmetric=.false.
!end if

!  initialize temporary storage (vector of pointers)
allocate(tja(ndof))

allocate(ia(ndof))
ia=0
do i=1,ndof
   allocate(tja(i).d(isize/2))
   tja(i)%d=0
end do

!  loop over all elements
do i=1,n
!  compute element matrix
   do k=1,size(eg(i,1)%d)
!  row index k: local;  kn: global
      kn=abs(eg(i,1)%d(k))
!  if this dof is not present
      if (kn .eq. 0) cycle
!  column index l: local;  ln: global
loop2:do l=1,k-1
         ln=abs(eg(i,1)%d(l))
         if (ln .eq. 0) cycle

         if (kn .lt. ln) then
      !  swap row and column -  column index must be smaller than row index
            lnn=kn
            knn=ln
         else
            lnn=ln
            knn=kn
         end if
      !  the actual size to hold this row
         acsize=size(tja(knn)%d)
         do mb=1,acsize
!  search for the the column in ja, if this dof was processed previously it is in the list
!  otherwise it is stored at the first free position (entry = 0)
            if (tja(knn)%d(mb) .eq. 0) then
            !  not found ==> add new entry
               tja(knn)%d(mb)=lnn
               ia(knn)=ia(knn)+1
               cycle loop2
            else if (tja(knn)%d(mb) .eq. lnn) then
            !  found 
               cycle loop2
            end if
         end do
!  not found no more space ==> resize
         tja(knn)%d=>reallocate(tja(knn)%d,2*acsize)
         tja(knn)%d(acsize+1)=lnn
         tja(knn)%d(acsize+2:2*acsize)=0
         ia(knn)=ia(knn)+1
         acsize=2*acsize
      end do loop2
   end do
end do
!  now we have the number of row entries stored in ia 
!  sum up to get ia in the CSR or LCSR format
do i = 2, ndof
   ia(i) = ia(i) + ia(i-1)
end do

!  generate ja (compact) from temporary tja (vector of pointers)
allocate(ja(1:ia(ndof)))
do i=2,ndof
   do j=1,ia(i)-ia(i-1)
      ja(ia(i-1)+j)=tja(i)%d(j)
   end do
end do

! destroy tja and free memory
lstat=destroyarrptr(tja)

! array for update global rhs
allocate(bglobal(ndof))
bglobal = 0._DP

! start assembly 
if (symmetric) then 
   allocate(lower(1:ia(ndof)), diag(ndof))
else
   allocate(lower(1:ia(ndof)), upper(1:ia(ndof)), diag(ndof))
!  clear matrix, right-hand-side
   upper=0._DP
end if
lower=0._DP
diag=0._DP

!  loop over all elements (it should be possible to execute this loop in parallel)
do i=1,n
!  compute element matrix
   call elementmatrix_pstiff(i,ep(i,1),ai,nff)  ! bi is computed in getrhsstep2_implicit
   
   do k=1,nff
!  row index k: local;  kn: global
      kn=abs(eg(i,1)%d(k))
      
!  if this dof is not present
      if (kn .eq. 0) cycle
      vzk=sign(1,eg(i,1)%d(k))

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
   deallocate(ai)
end do

end subroutine pstiff
