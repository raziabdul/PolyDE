subroutine setupia_temp(ia,csr,ntdof)
use femtypes
use feminterface
use fluidinterface
use globalvariables
implicit none
integer(I4B), intent(OUT) :: ia(:)
integer(I4B), intent(in) :: ntdof
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
!    $Revision: 1.6 $
!    $Date: 2014/07/15 13:02:35 $
!    $Author: m_kasper $
!
!#######################################################################################
!
! Generate the vector ia (or n1l) before assembly needed for the compact matrix storage
!
!#######################################################################################
! modified: 04.10.2004
! status  : finished for unsymmetric matrix
! remarks : use for Mixed method (STOKES solver)
!***************************************************************************************
! local variables
!***************************************************************************************
integer(I4B), dimension(3), parameter :: inach = (/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg = (/3,1,2/)
integer(I4B) :: i, j, k
integer(I4B) :: jj, kk, kkk, k1, k2
integer(I4B) :: edja, edjb, edka, edkb
integer(I4B) :: commonedge
integer(I4B) :: degree1, deg
integer(I4B) :: nff, udof, udof2, neldof, pdof
logical, pointer :: conode(:)
!***************************************************************************************
if (csr) then
   ia = 1
else
   ia = 0
end if

! use for only equal order of shape function for the whole domain.(NO_ADAPT, H_ADAPT)
nff = (ep(1,1)+1)*(ep(1,1)+2)/2
! elemental
udof = nff
pdof = (ep(1,1)+1)*ep(1,1)/2
udof2 = udof*2
neldof = udof2+pdof

allocate(conode(ntdof))
conode(:) = .FALSE.

do i = 1,n
!this is for only j=1,udof. there will be j = udof+1,2*udof and j = 2*udof+1,neldof
!####################################################################################
!   ROW loop for U
!####################################################################################
   do j = 1,udof

      k1 = abs(eg(i,1)%d(j))
      if (k1 == 0) cycle

      if (j .le. 3) then
         edja = inach(j)
         edjb = ivorg(j)
      else
         edja = -10
         edjb = -10
         do deg = 2,ep(i,1)
            if (j-(deg+1)*(deg+2)/2 .le. 0) then
               degree1 = deg
               exit
            end if
         end do
         if (j-(degree1)*(degree1+1)/2 .le.3) then
            edja = ivorg(j-(degree1)*(degree1+1)/2)
         end if
      end if
!*********************************************************************************
! column loop for U, corresponding to ROW Udof
!*********************************************************************************
      do k = j+1,udof

         k2 = abs(eg(i,1)%d(k))
         if (k2 .EQ. 0) cycle

         if (k .le. 3) then
            edka = inach(k)
            edkb = ivorg(k)
         else
            edka = -20
            edkb = -20
            do deg = 2,ep(i,1)
               if (k-(deg+1)*(deg+2)/2 .le. 0) then
                  degree1 = deg
                  exit
               end if
            end do
            if (k-(degree1)*(degree1+1)/2 .le.3) then
               edka = ivorg(k-(degree1)*(degree1+1)/2)
            end if
         end if

! determine if both dof j,k are on the same edge
         commonedge = 0

         if (edja .EQ. edka .OR. edja .EQ. edkb) then
            commonedge = edja
         elseif (edjb .EQ. edka .OR. edjb .EQ. edkb) then
            commonedge = edjb
         end if

         if (commonedge .GT. 0) then
            if (en(commonedge,i) .GT. i) then
               cycle
            end if
         end if

         if (csr) then
            ia(k1) = ia(k1)+1
            ia(k2) = ia(k2)+1
         else
            if (k2 > k1) then
               ia(k2) = ia(k2)+1
            else
               ia(k1) = ia(k1)+1
            end if
         end if
      end do ! endfor k loop U (column)
!*********************************************************************************
! column loop for P
!*********************************************************************************
      do k = udof2+1,neldof

         k2 = abs(eg(i,1)%d(k))
         if (k2 == 0) cycle
         kkk = k-udof2

         if (kkk .LE. 3) then ! note1:from here k is changed to kkk
            edka = inach(kkk)
            edkb = ivorg(kkk)
         else
            edka = -20
            edkb = -20
            do deg = 2,ep(i,1)
               if (kkk-(deg+1)*(deg+2)/2 .le. 0) then
                  degree1 = deg
                  exit
               end if
            end do

            if (kkk-(degree1)*(degree1+1)/2 .le.3) then
               edka = ivorg(kkk-(degree1)*(degree1+1)/2)
            end if
         end if ! note1:change of k to kkk ends here

! determine if both dof j,k are on the same edge
         commonedge = 0

         if (edja .EQ. edka .OR. edja .EQ. edkb) then
            commonedge = edja
         elseif (edjb .EQ. edka .OR. edjb .EQ. edkb) then
            commonedge = edjb
         end if

         if (commonedge .GT. 0) then
            if (en(commonedge,i) .GT. i) then
               cycle
            end if
         end if

         if (edja .EQ. edka .AND. edjb .EQ. edkb) then
            if (.NOT. conode(k1) ) then
               conode(k1) = .TRUE.
            else
               cycle
            end if
         end if

         if (csr) then
            ia(k1) = ia(k1)+1
            ia(k2) = ia(k2)+1
         else
            if (k2 > k1) then
               ia(k2) = ia(k2)+1
            else
               ia(k1) = ia(k1)+1
            end if
         end if
      end do ! endfor k loop P (column)
   end do ! endfor j loop for 1,udof
!####################################################################################
!   ROW loop for V
!####################################################################################
   do j = udof+1,udof2

      k1 = abs(eg(i,1)%d(j))
      if (k1 == 0) cycle
      jj = j-udof

      if (jj .le. 3) then ! note3:from here j is changed to jj
         edja = inach(jj)
         edjb = ivorg(jj)
      else
         edja = -10
         edjb = -10
         do deg = 2,ep(i,1)
            if (jj-(deg+1)*(deg+2)/2 .le. 0) then
               degree1 = deg
               exit
            end if
         end do
         if (jj-(degree1)*(degree1+1)/2 .le.3) then
            edja = ivorg(jj-(degree1)*(degree1+1)/2)
         end if
      end if ! note3:change of j to jj ends here
!*********************************************************************************
!   column loop for V 
!*********************************************************************************
      do k = j+1,udof2

         k2 = abs(eg(i,1)%d(k))
         if (k2 == 0) cycle
         kk = k-udof
         if (kk .le. 3) then ! note2:from here k is changed to kk
            edka = inach(kk)
            edkb = ivorg(kk)
         else
            edka = -20
            edkb = -20
            do deg = 2,ep(i,1)
               if (kk-(deg+1)*(deg+2)/2 .le. 0) then
                  degree1 = deg
                  exit
               end if
            end do
            if (kk-(degree1)*(degree1+1)/2 .le.3) then
               edka = ivorg(kk-(degree1)*(degree1+1)/2)
            end if
         end if ! note2:change of k to kk ends here

! determine if both dof j,k are on the same edge
         commonedge = 0

         if (edja .EQ. edka .OR. edja .EQ. edkb) then
            commonedge = edja
         elseif (edjb .EQ. edka .OR. edjb .EQ. edkb) then
            commonedge = edjb
         end if

         if (commonedge .GT. 0) then
            if (en(commonedge,i) .GT. i) then
               cycle
            end if
         end if
         
         if (csr) then
            ia(k1) = ia(k1)+1
            ia(k2) = ia(k2)+1
         else
            if (k2 > k1) then
               ia(k2) = ia(k2)+1
            else
               ia(k1) = ia(k1)+1
            end if
         end if
      end do ! endfor k loop V
!*********************************************************************************
!   column loop for P corresponding to Vdof
!*********************************************************************************
      do k = udof2+1,neldof

         k2 = abs(eg(i,1)%d(k))
         if (k2 == 0) cycle
         kkk = k-udof2
         if (kkk .le. 3) then ! note1:from here k is changed to kkk
            edka = inach(kkk)
            edkb = ivorg(kkk)
         else
            edka = -20
            edkb = -20

            do deg = 2,ep(i,1)
               if (kkk-(deg+1)*(deg+2)/2 .le. 0) then
                  degree1 = deg
                  exit
               end if
            end do

            if (kkk-(degree1)*(degree1+1)/2 .le.3) then
               edka = ivorg(kkk-(degree1)*(degree1+1)/2)
            end if
         end if ! note1:change of k to kkk ends here

! determine if both dof j,k are on the same edge
         commonedge = 0

         if (edja .EQ. edka .OR. edja .EQ. edkb) then
            commonedge = edja
         elseif (edjb .EQ. edka .OR. edjb .EQ. edkb) then
            commonedge = edjb
         end if

         if (commonedge .GT. 0) then
            if (en(commonedge,i) .GT. i) then
               cycle
            end if
         end if

         if (edja .EQ. edka .AND. edjb .EQ. edkb) then
            if (.NOT. conode(k1) ) then
               conode(k1) = .TRUE.
            else
               cycle
            end if
         end if

         if (csr) then
            ia(k1) = ia(k1)+1
            ia(k2) = ia(k2)+1
         else
            if (k2 > k1) then
               ia(k2) = ia(k2)+1
            else
               ia(k1) = ia(k1)+1
            end if
         end if
      end do ! endfor k loop P
   end do ! endfor j loop for udof+1,2*udof
!####################################################################################
!   ROW loop for P
!####################################################################################
! loop for P is not needed because the off-diagnal entries are always zero.
! Moreover, the diagonal entries have been defined at the beginning of programm where
! all ia(:) for full matrix is set to be 1
!************************************************************************************
end do ! end for element loop
deallocate(conode)
close(22)
!************************************************************************************
! put ia(:) in a desired form
!************************************************************************************
if (csr) then
   do i = ntdof,1,-1
      ia(i+1) = ia(i)
   end do
   ia(1) = 1
   do i = 2,ntdof+1
      ia(i) = ia(i)+ia(i-1)
   end do
else
   do i = 2,ntdof
      ia(i) = ia(i)+ia(i-1)
   end do
end if
!************************************************************************************
return
end subroutine setupia_temp