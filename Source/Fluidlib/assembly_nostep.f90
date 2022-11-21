subroutine assembly_nostep(lower,upper,diag,acsr,rhs,ia,ja,csr,npdof)
use femtypes
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer (I4B), pointer :: ia(:), ja(:)
integer(I4B), optional :: npdof
real(DP), pointer :: lower(:), upper(:), diag(:)
real(DP), pointer :: rhs(:), acsr(:)
logical, intent(in) :: csr
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
!    $Date: 2014/07/15 13:02:35 $
!    $Author: m_kasper $
!
!***************************************************************************************
!
!  Assembly routine for Stokes single step (equal-order primitive variables)
!
!***************************************************************************************
!  local variables
!***************************************************************************************
integer(I4B) :: i
integer(I4B) :: k,kn
integer(I4B) :: kk,kkn
integer(I4B) :: kkk,kkkn
integer(I4B) :: l,ln, ll, lll
integer(I4B) :: mb
integer(I4B) :: udof, udof2, pdof, nudof, nudof2, neldof, ntdof, nff !, npdof
integer(I4B) :: temp, errcode
real :: vzk, vzl
real(DP), POINTER :: KE(:,:), FE(:)
logical :: matvar
!***************************************************************************************
matvar = .false.

! use for only equal order of shape function for the whole domain.(NO_ADAPT, H_ADAPT)
nff = (ep(1,1)+1)*(ep(1,1)+2)/2
! elemental
udof = nff

if (optmixed) then
   if (ep(1,1) .LT. 2) then
      write(*,*) 'Polyorder must be at least 2 for Mixed Method'
      stop
   end if
   pdof = ep(1,1)*(ep(1,1)+1)/2
   ! must be change to the correct one
   !npdof = 0
else
   pdof = udof
   npdof = ndof
end if
udof2 = udof*2
neldof = udof2+pdof
! global
nudof = ndof
nudof2 = nudof*2
ntdof = nudof2+npdof

if (csr) then
!***************************************************************************************
!					FULL MATRIX
!***************************************************************************************
acsr = 0._DP
rhs = 0._DP
ja = 0

do i = 1,n ! loop over element
   
   if (optmixed) then
      call elementmatrix_stokesmix(i,ep(i,1),.false.,.true.,.false.,KE,FE,temp,errcode)
   else
      call elementmatrix_stokes(i,KE,FE,matvar)
   end if

   do k = 1,udof  ! loop over row

      kk = k + udof
      kn = abs(eg(i,1)%d(k))
      kkn = kn + nudof

      if (kn == 0) cycle
      vzk = sign(1,eg(i,1)%d(k))

      rhs(kn) = rhs(kn) + FE(k)*vzk
      rhs(kkn) = rhs(kkn) + FE(kk)*vzk
      if (k <= pdof) then
         kkk = k + udof2
         kkkn = kn + nudof2
         rhs(kkkn) = rhs(kkkn) + FE(kkk)*vzk
      endif
!******************************************************************
!   at k-th row for U dof -> KG( 1:nudof , 1:nudof )
!******************************************************************
ucolk: do l = 1,udof
         if (k == l) then
            acsr(ia(kn)) = acsr(ia(kn)) + KE(k,k)
            ja(ia(kn)) = kn
         else
            ln = abs(eg(i,1)%d(l))
            if (ln == 0) cycle
            vzl = sign(1,eg(i,1)%d(l))

            do mb = ia(kn)+1,ia(kn+1)-1
               if (ja(mb) == 0) then
                  ja(mb) = ln
                  acsr(mb) = KE(k,l)*vzk*vzl
                  cycle ucolk
               else if (ja(mb) == ln) then
                  acsr(mb) = acsr(mb) + KE(k,l)*vzk*vzl
                  cycle ucolk
               end if
            end do
            write(*,*)'***** ERROR in MATRIX ASSEMBLY *****'
            write(*,*)'AT ELEMENT:',i
            write(*,*)'AT GLOBAL ROW:',kn,'COLUMN:',ln
         end if
      end do ucolk
!******************************************************************
!   at k-th row for V dof
!******************************************************************
! KG( 1:nudof , nudof+1:2*nudof ) = 0.0 
!******************************************************************
!   at k-th row for P dof -> KG( 1:nudof , 2*nudof+1:ntdof ) 
!******************************************************************
pcolk: do l = udof2+1,neldof
         if (k == l) then
            acsr(ia(kn)) = acsr(ia(kn)) + KE(k,k)
            ja(ia(kn)) = kn
         else
            lll = l - udof2
            ln = abs(eg(i,1)%d(lll)) + nudof2
            if (ln == 0) cycle
            vzl = sign(1,eg(i,1)%d(lll))

            do mb = ia(kn)+1,ia(kn+1)-1
               if (ja(mb) == 0) then
                  ja(mb) = ln
                  acsr(mb) = KE(k,l)*vzk*vzl
                  cycle pcolk
               else if (ja(mb) == ln) then
                  acsr(mb) = acsr(mb) + KE(k,l)*vzk*vzl
                  cycle pcolk
               end if
            end do
            write(*,*)'***** ERROR in MATRIX ASSEMBLY *****'
            write(*,*)'AT ELEMENT:',i
            write(*,*)'AT GLOBAL ROW:',kn,'COLUMN:',ln
         end if
      end do pcolk
!******************************************************************
!   at kk-th row for U dof
!******************************************************************
! KG(nudof+1:2*nudof , 1:nudof) = 0.0 
!******************************************************************
!   at kk-th row for V dof -> KG(nudof+1:2*nudof , nudof+1:2*nudof)
!******************************************************************
vcolkk: do l = udof+1,udof2
         if (kk == l) then
            acsr(ia(kkn)) = acsr(ia(kkn)) + KE(kk,kk)
            ja(ia(kkn)) = kkn
         else
            ll = l - udof
            ln = abs(eg(i,1)%d(ll)) + nudof
            if (ln == 0) cycle
            vzl = sign(1,eg(i,1)%d(ll))

            do mb = ia(kkn)+1,ia(kkn+1)-1
               if (ja(mb) == 0) then
                  ja(mb) = ln
                  acsr(mb) = KE(kk,l)*vzk*vzl
                  cycle vcolkk
               else if (ja(mb) == ln) then
                  acsr(mb) = acsr(mb) + KE(kk,l)*vzk*vzl
                  cycle vcolkk
               end if
            end do
            write(*,*)'***** ERROR in MATRIX ASSEMBLY *****'
            write(*,*)'AT ELEMENT:',i
            write(*,*)'AT GLOBAL ROW:',kkn,'COLUMN:',ln
         end if
      end do vcolkk
!******************************************************************
!   at kk-th row for P dof -> KG(nudof+1:2*nudof , 2*nudof+1:ntdof)
!******************************************************************
pcolkk: do l = udof2+1,neldof
         if (kk == l) then
            acsr(ia(kkn)) = acsr(ia(kkn)) + KE(kk,kk)
            ja(ia(kkn)) = kkn
         else
            lll = l - udof2
            ln = abs(eg(i,1)%d(lll)) + nudof2
            if (ln == 0) cycle
            vzl = sign(1,eg(i,1)%d(lll))

            do mb = ia(kkn)+1,ia(kkn+1)-1
               if (ja(mb) == 0) then
                  ja(mb) = ln
                  acsr(mb) = KE(kk,l)*vzk*vzl
                  cycle pcolkk
               else if (ja(mb) == ln) then
                  acsr(mb) = acsr(mb) + KE(kk,l)*vzk*vzl
                  cycle pcolkk
               end if
            end do
            write(*,*)'***** ERROR in MATRIX ASSEMBLY *****'
            write(*,*)'AT ELEMENT:',i
            write(*,*)'AT GLOBAL ROW:',kkn,'COLUMN:',ln
         end if
      end do pcolkk
!******************************************************************
!   at kkk-th row for U dof -> KG(2*nudof+1:ntdof , 1:nudof)
!******************************************************************
    if (k <= pdof) then

ucolkkk:do l = 1,udof

         if (kkk == l) then
            acsr(ia(kkkn)) = acsr(ia(kkkn)) + KE(kkk,kkk)
            ja(ia(kkkn)) = kkkn
         else
            ln = abs(eg(i,1)%d(l))
            if (ln == 0) cycle
            vzl = sign(1,eg(i,1)%d(l))

            do mb = ia(kkkn)+1,ia(kkkn+1)-1
               if (ja(mb) == 0) then
                  ja(mb) = ln
                  acsr(mb) = KE(kkk,l)*vzk*vzl
                  cycle ucolkkk
               else if (ja(mb) == ln) then
                  acsr(mb) = acsr(mb) + KE(kkk,l)*vzk*vzl
                  cycle ucolkkk
               end if
            end do
            write(*,*)'***** ERROR in MATRIX ASSEMBLY *****'
            write(*,*)'AT ELEMENT:',i
            write(*,*)'AT GLOBAL ROW:',kkkn,'COLUMN:',ln
         end if
      end do ucolkkk
!******************************************************************
!   at kkk-th row for V dof -> KG(2*nudof+1:ntdof , nudof+1:2*nudof)
!******************************************************************
vcolkkk:do l = udof+1,udof2

         if (kkk == l) then
            acsr(ia(kkkn)) = acsr(ia(kkkn)) + KE(kkk,kkk)
            ja(ia(kkkn)) = kkkn
         else
            ll = l - udof
            ln = abs(eg(i,1)%d(ll)) + nudof
            if (ln == 0) cycle
            vzl = sign(1,eg(i,1)%d(ll))

            do mb = ia(kkkn)+1,ia(kkkn+1)-1
               if (ja(mb) == 0) then
                  ja(mb) = ln
                  acsr(mb) = KE(kkk,l)*vzk*vzl
                  cycle vcolkkk
               else if (ja(mb) == ln) then
                  acsr(mb) = acsr(mb) + KE(kkk,l)*vzk*vzl
                  cycle vcolkkk
               end if
            end do
            write(*,*)'***** ERROR in MATRIX ASSEMBLY *****'
            write(*,*)'AT ELEMENT:',i
            write(*,*)'AT GLOBAL ROW:',kkkn,'COLUMN:',ln
         end if
      end do vcolkkk
!******************************************************************
!   at kkk-th row for P dof
!******************************************************************
pcolkkk:do l = udof2+1,neldof

         if (kkk == l) then
            acsr(ia(kkkn)) = acsr(ia(kkkn)) + KE(kkk,kkk)
            ja(ia(kkkn)) = kkkn
         else
            if (optmixed) cycle
            lll = l - udof2
            ln = abs(eg(i,1)%d(lll)) + nudof2
            if (ln == 0) cycle
            vzl = sign(1,eg(i,1)%d(lll))
            do mb = ia(kkkn)+1,ia(kkkn+1)-1
               if (ja(mb) == 0) then
                  ja(mb) = ln
                  acsr(mb) = KE(kkk,l)*vzk*vzl
                  cycle pcolkkk
               else if (ja(mb) == ln) then
                  acsr(mb) = acsr(mb) + KE(kkk,l)*vzk*vzl
                  cycle pcolkkk
               end if
            end do
            write(*,*)'***** ERROR in MATRIX ASSEMBLY *****'
            write(*,*)'AT ELEMENT:',i
            write(*,*)'AT GLOBAL ROW:',kkkn,'COLUMN:',ln
         end if
      end do pcolkkk

    end if ! comparing k <= pdof
!******************************************************************
   end do ! k, loop over row
   deallocate(KE,FE)
end do ! i, loop over element

end if ! if csr

end subroutine assembly_nostep
