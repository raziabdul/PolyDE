subroutine setupia_stokes(ia,csr)
use femtypes
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer(I4B), intent(out) :: ia(:)
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
! remarks :
!   5. The invalid entries for unsymmeric matrix (csr = .true.) is corrected.
!      It is because of the pressure node, located on the same position of the node,
!      considered in j-loop (for k1).
!   4. No. of invalid entries in UMFPACK is dramatically reduced. However, there are
!      still some left. (27.08.04)
!   3. In UMFPACK some error indicating invalid entries appeared and are ignored during
!      solution step.
!   2. Array ia is modified in csrsetupia2 and run without any error in assembly step.
!   1. Need some implementation for U,V,P -> done
!#######################################################################################
!  Input:
!       en      =   neighborhood information, neighbors of the elements
!       n       =   total number of elements
!       nnode   =   total number of nodes
!       ndof    =   total number of degrees of freedom
!       ep      =   polynomial degree for each of the elements
!       eg      =   degrees of freedom (global) for the elements
!                   packed into an arry of pointers   eg(i)%d(j)
!       csr     =   ia will return the number of non-zero entries in each row of the 
!                   -  full matrix  if                  .true. 
!                   -  lower triangular matrix if       .false.
!  Output:
!       ia      =   a vector for the compact storage of the sparse matrix,
!                   it contains the accumulated sum of the number of non-zero 
!                   elements in matrix rows of the full or lower triangular matrix
!
!  The routine computes the compact storage information for two common storage schemes:
!  the compressed sparse row format: CSR  if csr=.true. or 
!  the compressed lower sparse row format if csr=.false.
!  The CSR format is very common in sparse linear algebra packages however, the 
!  second format needs less memory 
!  
!                 / -2  1  0  4 \
!                 |  8 -8  0  0 |
!  The Matrix     |  0  0 -6  6 |
!                 \  3  0  5 -5 /
!  
!  is stored in the compressed sparse row format: CSR (sometimes termed Harwell-Boeing 
!  compact storage scheme ) by using three vectors 
!  The vector or non-zero matrix entries:             A =(-2, 4, 1, 8, -8, 6, -6, 5, -5, 3)
!  A vector of indices of the entries in A:           ja =( 1, 4, 2, 1,  2, 3,  4, 3,  4, 1)
!  The vector of pointers to the first entry of rows: ia =( 1, 4, 6, 8, 11)
!  
!  The compressed lower sparse row format take advantage of the fact that in FEM the 
!  occupation pattern is symmetric i.e. if a(i,j)=0 then a(j,i)=0. This format uses 
!  three vectors to hold the non-zero matrix entries and two index vectors
!  The lower-triangular matrix                         L = ( 8, 5, 3)
!  The transposed upper triangular matrix              U = ( 1, 6, 4)
!  The diagonal:                                       D = (-2, -8, -6, -5)
!  A vector of indices of the entries in L:          ja = ( 1, 3, 1)
!  The vector of pointers to the last entry of rows: ia = ( 0, 1, 1, 3)
!  In this format the vector n2l is significatly smaller
!
!  Note that n1l has different interpretation and length for the two storage schemes
!**********************************************************************************
! local variables
!**********************************************************************************
! inach = Nachfolger von i
! ivorg = Vorgaenger von i
integer(I4B), dimension(3), parameter :: inach = (/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg = (/3,1,2/)

integer(I4B) :: i, j, k
integer(I4B) :: jj, kk, kkk, k1, k2
integer(I4B) :: edja, edjb, edka, edkb
integer(I4B) :: commonedge
integer(I4B) :: degree1, deg
integer(I4B) :: nff, udof, udof2, neldof, nudof, nudof2, ntdof
logical, pointer :: conode(:)
!**********************************************************************************
if (csr) then
   ia = 1
else
   ia = 0
endif

! use for only equal order of shape function for the whole domain.(NO_ADAPT, H_ADAPT)
nff = (ep(1,1)+1)*(ep(1,1)+2)/2
! elemental DOF
udof = nff
!pdof = udof
udof2 = udof*2
neldof = udof2+udof
! global DOF
nudof = ndof
nudof2 = nudof*2
ntdof = ndof*3

allocate(conode(ntdof))
conode(:) = .FALSE.

do i = 1,n
!this is for only j=1,udof. there will be j = udof+1,2*udof and j = 2*udof+1,neldof
!####################################################################################
!  ROW loop for U
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
            endif
         enddo
         if (j-(degree1)*(degree1+1)/2 .le.3) then
            edja = ivorg(j-(degree1)*(degree1+1)/2)
         endif
      endif
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
               endif
            enddo
            if (k-(degree1)*(degree1+1)/2 .le.3) then
               edka = ivorg(k-(degree1)*(degree1+1)/2)
            endif
         endif

! determine if both dof j,k are on the same edge
         commonedge = 0

         if (edja .EQ. edka .OR. edja .EQ. edkb) then
            commonedge = edja
         elseif (edjb .EQ. edka .OR. edjb .EQ. edkb) then
            commonedge = edjb
         endif

         if (commonedge .GT. 0) then
            if (en(commonedge,i) .GT. i) then
               cycle
            endif
         endif

         if (csr) then
            ia(k1) = ia(k1)+1
            ia(k2) = ia(k2)+1
         else
            if (k2 > k1) then
               ia(k2) = ia(k2)+1
            else
               ia(k1) = ia(k1)+1
            endif
         endif
      enddo ! endfor k loop U (column)
!*********************************************************************************
! column loop for P
!*********************************************************************************
      do k = udof2+1,neldof

         kkk = k-udof2
         k2 = abs(eg(i,1)%d(kkk)) + nudof2
         if (k2 == 0) cycle

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
               endif
            enddo

            if (kkk-(degree1)*(degree1+1)/2 .le.3) then
               edka = ivorg(kkk-(degree1)*(degree1+1)/2)
            endif
         endif ! note1:change of k to kkk ends here

! determine if both dof j,k are on the same edge
         commonedge = 0

         if (edja .EQ. edka .OR. edja .EQ. edkb) then
            commonedge = edja
         elseif (edjb .EQ. edka .OR. edjb .EQ. edkb) then
            commonedge = edjb
         endif

         if (commonedge .GT. 0) then
            if (en(commonedge,i) .GT. i) then
               cycle
            endif
         endif

         if (edja .EQ. edka .AND. edjb .EQ. edkb) then
            if (.NOT. conode(k1) ) then
               conode(k1) = .TRUE.
            else
               cycle
            endif
         endif

         if (csr) then
            ia(k1) = ia(k1)+1
            ia(k2) = ia(k2)+1
         else
            if (k2 > k1) then
               ia(k2) = ia(k2)+1
            else
               ia(k1) = ia(k1)+1
            endif
         endif
      enddo ! endfor k loop P (column)
   enddo ! endfor j loop for 1,udof
!####################################################################################
!   ROW loop for V
!####################################################################################
   do j = udof+1,udof2

      jj = j-udof
      k1 = abs(eg(i,1)%d(jj)) + nudof
      if (k1 == 0) cycle

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
            endif
         enddo
         if (jj-(degree1)*(degree1+1)/2 .le.3) then
            edja = ivorg(jj-(degree1)*(degree1+1)/2)
         endif
      endif ! note3:change of j to jj ends here
!*********************************************************************************
!   column loop for V 
!*********************************************************************************
      do k = j+1,udof2

         kk = k-udof
         k2 = abs(eg(i,1)%d(kk)) + nudof
         if (k2 == 0) cycle

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
               endif
            enddo
            if (kk-(degree1)*(degree1+1)/2 .le.3) then
               edka = ivorg(kk-(degree1)*(degree1+1)/2)
            endif
         endif ! note2:change of k to kk ends here

! determine if both dof j,k are on the same edge
         commonedge = 0

         if (edja .EQ. edka .OR. edja .EQ. edkb) then
            commonedge = edja
         elseif (edjb .EQ. edka .OR. edjb .EQ. edkb) then
            commonedge = edjb
         endif

         if (commonedge .GT. 0) then
            if (en(commonedge,i) .GT. i) then
               cycle
            endif
         endif
         
         if (csr) then
            ia(k1) = ia(k1)+1
            ia(k2) = ia(k2)+1
         else
            if (k2 > k1) then
               ia(k2) = ia(k2)+1
            else
               ia(k1) = ia(k1)+1
            endif
         endif
      enddo ! endfor k loop V
!*********************************************************************************
!   column loop for P corresponding to Vdof
!*********************************************************************************
      do k = udof2+1,neldof

         kkk = k-udof2
         k2 = abs(eg(i,1)%d(kkk)) + nudof2
         if (k2 == 0) cycle

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
               endif
            enddo

            if (kkk-(degree1)*(degree1+1)/2 .le.3) then
               edka = ivorg(kkk-(degree1)*(degree1+1)/2)
            endif
         endif ! note1:change of k to kkk ends here

! determine if both dof j,k are on the same edge
         commonedge = 0

         if (edja .EQ. edka .OR. edja .EQ. edkb) then
            commonedge = edja
         elseif (edjb .EQ. edka .OR. edjb .EQ. edkb) then
            commonedge = edjb
         endif

         if (commonedge .GT. 0) then
            if (en(commonedge,i) .GT. i) then
               cycle
            endif
         endif

         if (edja .EQ. edka .AND. edjb .EQ. edkb) then
            if (.NOT. conode(k1) ) then
               conode(k1) = .TRUE.
            else
               cycle
            endif
         endif

         if (csr) then
            ia(k1) = ia(k1)+1
            ia(k2) = ia(k2)+1
         else
            if (k2 > k1) then
               ia(k2) = ia(k2)+1
            else
               ia(k1) = ia(k1)+1
            endif
         endif
      enddo ! endfor k loop P
   enddo ! endfor j loop for udof+1,2*udof
!####################################################################################
!   ROW loop for P
!####################################################################################
   do j = udof2+1,neldof

      jj = j-udof2
      k1 = abs(eg(i,1)%d(jj)) + nudof2      
      if (k1 == 0) cycle

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
            endif
         enddo
         if (jj-(degree1)*(degree1+1)/2 .le.3) then
            edja = ivorg(jj-(degree1)*(degree1+1)/2)
         endif
      endif ! note3:change of j to jj ends here
!*********************************************************************************
!   column loop for P corresponding to Pdof
!*********************************************************************************
      do k = j+1,neldof
         kkk = k-udof2
         k2 = abs(eg(i,1)%d(kkk)) + nudof2
         if (k2 == 0) cycle

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
               endif
            enddo

            if (kkk-(degree1)*(degree1+1)/2 .le.3) then
               edka = ivorg(kkk-(degree1)*(degree1+1)/2)
            endif
         endif ! note1:change of k to kkk ends here

! determine if both dof j,k are on the same edge
         commonedge = 0

         if (edja .EQ. edka .OR. edja .EQ. edkb) then
            commonedge = edja
         elseif (edjb .EQ. edka .OR. edjb .EQ. edkb) then
            commonedge = edjb
         endif

         if (commonedge .GT. 0) then
            if (en(commonedge,i) .GT. i) then
               cycle
            endif
         endif

         if (csr) then
            ia(k1) = ia(k1)+1
            ia(k2) = ia(k2)+1
         else
            if (k2 > k1) then
               ia(k2) = ia(k2)+1
            else
               ia(k1) = ia(k1)+1
            endif
         endif
      enddo ! endfor k loop P
   enddo ! endfor j loop for udof2+1, neldof

enddo ! endfor element loop

deallocate(conode)
!************************************************************************************
! put ia(:) in a desired form
!************************************************************************************
if (csr) then
   do i = ntdof,1,-1
      ia(i+1) = ia(i)
   enddo
   ia(1) = 1
   do i = 2,ntdof+1
      ia(i) = ia(i)+ia(i-1)
   enddo
else
   do i = 2,ntdof
      ia(i) = ia(i)+ia(i-1)
   enddo
endif
!************************************************************************************
return
endsubroutine setupia_stokes