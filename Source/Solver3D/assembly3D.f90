      subroutine assembly3D(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr, &
     &                      matvar,symmetric)
      use feminterface, only: zeit, reallocate, destroyarrptr, palloc
      use feminterface3D, only: elementmatrix3D, elementmatrix_scalar3D, coloringsc3D
      use femtypes
      use globalvariables3D, only : eltype, nnat, numdof, numv, vgdof, vp, nod
      implicit none
      integer (I4B), pointer :: ia(:), ja(:)
      complex (DPC), pointer :: lower(:), upper(:), diag(:), rhs(:), acsr(:)
      logical jacobi, csr, matvar, symmetric
      intent (in) ::  jacobi, csr, matvar, symmetric
!
!-------------------------------------------------------------------------------
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
!    $Revision: 1.22 $
!    $Date: 2015/11/10 15:36:14 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'assembly3D'
!  Performs element by element FEM matrix and right-hand-side assembly in compact storage
!__
!  The routine computes the compact storage information for two common storage schemes:
!  the compressed sparse row format:        CSR  if csr=.true. or 
!  the compressed lower sparse row format: LCSR  if csr=.false.
!
!  The CSR format is very common in sparse linear algebra packages however, the 
!  second format needs less memory 
!  
!  
!                 / -2  1  0  4 \
!                 |  8 -8  0  0 |
!  The Matrix     |  0  0 -6  6 |
!                 \  3  0  5 -5 /
!  
!  is stored in the compressed sparse row format: CSR by using three vectors
!  The vector or non-zero matrix entries:           acsr=(-2, 4, 1, 8, -8, 6, -6, 5, -5, 3)
!  A vector of indices of the entries in A:           ja=( 1, 4, 2, 1,  2, 3,  4, 3,  4, 1)
!  The vector of pointers to the first entry of rows: ia=( 1, 4, 6, 8, 11)
!
!  Actually we use a slightly modified CSR. The difference to the standard one lies in the
!  accessibility of the diagonal entries. They can be directly accessed by the values given
!  in vector ia in this version. So the diagonal entries for each row are always in first
!  place (referred to by ia).
!  example: the diagonal matrix entry no. A(3,3) is addressed by ia(3). So it can be found
!           at acsr(6) since 6 is the value given in ia(3).
!
!  The vector or non-zero matrix entries:           acsr=(-2, 4, 1, -8, 8, -6, 6, -5, 3, 5)
!  A vector of indices of the entries in A:           ja=( 1, 4, 2,  2, 1,  3, 4,  4, 1, 3)
!  The vector of pointers to the first entry of rows: ia=( 1, 4, 6, 8, 11)
!
!
!  The compressed lower sparse row format (LCSR) takes advantage of the fact that in FEM the 
!  occupation pattern is symmetric i.e. if A(i,j)=0 then A(j,i)=0. This format uses 
!  three vectors to hold the non-zero matrix entries and two index vectors
!  The lower-triangular matrix                         L=( 8, 5, 3)
!  The transposed upper triangular matrix              U=( 1, 6, 4)
!  The diagonal:                                       D=(-2, -8, -6, -5)
!  A vector of column indices of the entries in L:     ja=( 1, 3, 1)
!  The vector of pointers to the last entry of rows:   ia=( 0, 1, 1, 3)
!  In this format the vector n2l is significantly smaller (roughly half the size of CSR)
!
!  Note that ia has different interpretation and length for the two storage schemes!
!
!-------------------------------------------------------------------------------
!
!  Input:  
!                     / = .true. for the computation of the jacobian
!            jacobi
!                     \ = .false. for the linear case
!            numv     total number of elements
!            vgdof    degrees of freedom (global) for the elements for all natures
!                     packed into an array of pointers   vgdof(i,inat)%d(j)
!            ep       polynomial degree for each of the elements and for all natures
!            numdof   total number of degrees of freedom
!            csr      ia contains the number of non-zero entries in each row of the 
!                       -  full matrix  if                  .true. 
!                       -  lower triangular matrix if       .false.
!            matvar   =.true. if the material coefficients are varying across the element
!                     =.false. the material is assumed to be constant 
!            symmetric=.true. matrix is assumed to be symmetric, only the 
!                             lower triangular matrix and diagonal will be assembled
!                             This only works for LCSR-format, 
!                             CSR always assembles the full matrix 
!
!  Output:  
!            lower    lower triangular matrix
!            upper    upper triangular matrix
!            diag     diagonal matrix
!            acsr     vector (compact storage) of non-zero matrix entries
!            rhs      right-hand-side
!            ia       is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            ja       vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix   if csr = .false.
!                       - vector of pointers to the first entry of rows if csr = .true.
!  In-/ Output:  
!            x        actual solution
!
!-------------------------------------------------------------------------------
! I N F O :
!  - 'iloop2' and 'iloop4' : These indixes take care of the multinature LCSR assembly for multinature symmetric matrices
!  - 'inat'                : Index of nature
!  - 'jnat'                : Nested nature index
!  - 'isize'               : This parameter determines the initial size to hold one matrix row
!                            it should be chosen a bit larger than the expected mean number of row entries 
!                            the minimum value (not recommended) is 2
!                            small values cause excessive reallocation and possibly memory fragmentation
!                            large values may (temporarily) allocate much memory
!  - 'indx' and 'offset'   : These pointers are allocated in coloring
!
! N O T E S :
!    On parallel processing:
!                            Observe the OpenMP conditional compilation marker '!$' 
!                            -> Current code is using coloring.
!                            -> Lines must be commented out to remove the coloring loop and variables.
!                            -> Critical sections must be reintroduced if coloring is not used. 
!
!-------------------------------------------------------------------------------
!  Local variables:
      integer (I4B) :: e, i, j, k, l, kn, knn, ln, lnn, k1, l1, mb, acsize
      integer (I4B) iloop2, iloop4, errcode, vzk, vzl
      integer (I4B) count
      integer (I4B) inat, jnat
      integer (I4B), parameter :: isize=28
      integer (I4B), allocatable :: nff(:), nffsum(:)
    ! for single nature nedelec elements only:
      integer (I4B) nffone, vpone
      integer (I4B), allocatable :: indx(:), offset(:)
      integer (I4B) color, cstart, cend, elem, ncolor
      complex (DPC), allocatable :: ai(:,:), bi(:)
      type (ARRPTRI), pointer :: tja(:)=>null()
      logical :: lstat
!__
! Definitions:
   !print "(A79)"        ,' _______________________________________________________________________________ '
   !print "(A46)"        ,'|                           STARTING  ASSEMBLY'
   !print "(A2)"         ,'| '
!__
! Checks and preparations:
!__
!
! Allocate Memory:
!
! Initialize temporary storage (vector of pointers)
! tja(ndof) is a vector where every (say) subvector holds the column indices
! of the matrix entries for one row
      call palloc(tja,numdof)
!__
! Start:
!__
! 1) Perform the Coloring:
!$    call coloringsc3D(indx,offset)     ! <-- comment out to skip the coloring
!__
!  
! 2) Case selection of Matrix storage format: CSR / LCSR
! a) allocation of 'isize' spaces for every 'tja' subvector
!    'isize'/2 in case of LCSR
      if (csr) then
        call palloc(ia,numdof+1)
        ia=1
        do i=1,numdof
!-      a)
          call palloc(tja(i).d,isize) !for modified csr, entries for diagonal are in first place
          tja(i)%d(1)=i
          tja(i)%d(2:isize)=0
        end do
      else
        call palloc(ia,numdof)
        ia=0
        do i=1,numdof
!-      a)
          call palloc(tja(i).d,isize/2)
          tja(i)%d=0
        end do
      end if
!__
!  
! 3) Determine number of row entries and column indices
!    Start color loop if parallel mode
      if (csr) then
        ncolor=1                       !serial
!$      ncolor=size(offset)-1
!$      do color=1,ncolor              ! <-- comment out to skip the coloring
          cstart=1                     !serial
          cend=numv                    !serial
!$        cstart=offset(color)         ! <-- comment out to skip the coloring
!$        cend  =offset(color+1)-1     ! <-- comment out to skip the coloring
!
!$omp parallel do default(none) &
!$omp shared( nnat, tja, ia, vgdof, cstart, cend, indx)                 &
!$omp private( ln, kn, i, inat, k, jnat, l, mb, elem, acsize)
!__
!  
! 4) Work through all DOFs:
!  a) Loop over all elements (loop in parallel, range given by color/nocolor loop)
!  b) Loop over all natures
!  c) 'k' iterates until size of the element matrix (number of rows)
!     row index 'k' is the local number, 'kn' is the global number
!     'kn' = 0 if it's DOF is not present
!  d) 'l'  is the local column index
!     'ln' is the global column index
!     loop over all natures once more, inat is the nature index in the direction of
!     increasing row number whereas jnat in the direction of increasing column number
!     there are nnat**2 submatrices in the global matrix
!  e) This is for modified csr, entries for diagonal are in first place
!     As global diagonal entries we consider only the local diagonal entries of the
!     say diagonal matrix blocks, in the multinature global matrix we have
!     purely nature parts and coupled nature parts
!     the coupled nature parts are located off the global diagonal
!  f) Search for the column in ja, if this dof was processed previously it is in the list
!     otherwise it is stored at the first free position (entry = 0)
!     if this DOF was not found ==> add new entry
!  g) If no more space ==> resize
!     if there is no more space in every tja subvector to hold the column indices
!     of the entries of one row in the global matrix then allocate double the space
!     for storage to this tja subvector
        do elem=cstart,cend
!-      a)
          i=elem                       !serial
!$        i=indx(elem)                 ! <-- comment out to skip the coloring
!-      b)
          do inat=1,nnat
!-      c)
            do k=1,size(vgdof(i,inat)%d)
              kn=abs(vgdof(i,inat)%d(k))
              if (kn .eq. 0) cycle
              acsize=size(tja(kn)%d) ! the actual size to hold this row
!-      d)
              do jnat=1,nnat
loop1:          do l=1,size(vgdof(i,jnat)%d)
!-      e)
                  if (( k .ne. l) .or. (inat .ne. jnat) ) then 
                    ln=abs(vgdof(i,jnat)%d(l))
                    if (ln .eq. 0) cycle
!-
                    do mb=2,acsize
!-      f)
                      if (tja(kn)%d(mb) .eq. 0) then ! DOF not found
                        tja(kn)%d(mb)=ln
                        ia(kn)=ia(kn)+1
                        cycle loop1
                      else if (tja(kn)%d(mb) .eq. ln) then ! DOF found
                        cycle loop1
                      end if
                    end do
!-      g)
                    tja(kn)%d=>reallocate(tja(kn)%d,2*acsize)
                    tja(kn)%d(acsize+1)=ln
                    tja(kn)%d(acsize+2:2*acsize)=0
                    ia(kn)=ia(kn)+1
                    acsize=2*acsize
                  end if
                end do loop1
              end do
            end do
          end do
        end do
!$omp end parallel do
!$      end do                         ! end color loop <-- comment out to skip the coloring
!__
!
      else ! NOT CSR - then it's LCSR
!__
!  
! 5) Determine number of row entries and column indices
!    Start color loop if parallel mode
        ncolor=1                       !serial
!$      ncolor=size(offset)-1
!
! Start color loop
!$      do color=1,ncolor              ! <-- comment out to skip the coloring
          cstart=1                     !serial
          cend=numv                    !serial
!$        cstart=offset(color)         ! <-- comment out to skip the coloring
!$        cend  =offset(color+1)-1     ! <-- comment out to skip the coloring
!
!$omp parallel do default(none) &
!$omp shared( nnat, tja, ia, vgdof, cstart, cend, indx)                 &
!$omp private( ln, kn, i, inat, k, jnat, l, mb, elem, acsize)           &
!$omp private( iloop2, lnn, knn )
!__
!  
! 6) Work through all DOFs:
!
!  a) Loop over all elements (loop in parallel, range given by color/nocolor loop)
!-
!  b) Loop over all natures
!-
!  c) 'k' iterates until size of the element matrix (number of rows)
!     row index 'k' is the local number, 'kn' is the global number
!     'kn' = 0 if it's DOF is not present
!-
!  d) 'l'  is the local column index
!     'ln' is the global column index
!     loop over all natures once more, inat is the nature index in the direction of
!     increasing row number whereas jnat in the direction of increasing column number
!     there are nnat**2 submatrices in the global matrix
!       Additional information:
!     For loop2 we count up to the value iloop2
!     this is because we want to consider the local diagonal entries only for the
!     submatrix blocks where index inat is smaller than 'jnat'
!     this submatrix block lies above the diagonal, the local diagonal entries of the block
!     lying below the global diagonal will be taken care of by swaping the indices
!     of the symmetric entries which lie above the global diagonal
!-
!  e) If the global row number is smaller than the global column number, then
!     swap row and column - column index must be smaller than row index
!-
!  f) search for the the column in ja, if this dof was processed previously it is in the list
!     otherwise it is stored at the first free position (entry = 0)
!     if this DOF was not found ==> add new entry
        do elem=cstart,cend
!-      a)
          i=elem                       !serial
!$        i=indx(elem)                 ! <-- comment out to skip the coloring
!       b)
          do inat=1,nnat
!-      c)
            do k=1,size(vgdof(i,inat)%d)
              kn=abs(vgdof(i,inat)%d(k))
              if (kn .eq. 0) cycle
!-      d)
              do jnat=1,nnat
                if (inat .eq. jnat) then 
                  iloop2=k-1
                else if (inat .gt. jnat) then
                  iloop2=size(vgdof(i,jnat)%d)
                else
                  iloop2=0
                end if
loop2:          do l=1,iloop2
!-      e)
                  ln=abs(vgdof(i,jnat)%d(l))
                  if (ln .eq. 0) cycle
                  if (kn .lt. ln) then ! The global row number is smaller than the global column number
                    lnn=kn
                    knn=ln
                  else
                    lnn=ln
                    knn=kn
                  end if
                  acsize=size(tja(knn)%d) ! The actual size to hold this row
                  do mb=1,acsize
!-      f)
                    if (tja(knn)%d(mb) .eq. 0) then ! not found ==> add new entry
                      tja(knn)%d(mb)=lnn
                      ia(knn)=ia(knn)+1
                      cycle loop2
                    else if (tja(knn)%d(mb) .eq. lnn) then ! found 
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
          end do
        end do
!$omp end parallel do
!
!$      end do                                                                 ! end color loop <-- comment out to skip the coloring
      end if
!__
!
! 7) Now we have the number of row entries stored in ia,
!    sum up to get ia in the CSR or LCSR format
!    If CSR-format, generate ja (compact) from temporary tja (vector of pointers)
      if (csr) then
        do i = numdof, 1, -1
          ia(i+1) = ia(i)
        end do
        ia(1)=1
        do i = 2, numdof+1
          ia(i) = ia(i) + ia(i-1)
        end do
      else
        do i = 2, numdof
          ia(i) = ia(i) + ia(i-1)
        end do
      end if
!-
      if (csr) then ! case of CSR if csr = .true.
        call palloc(ja,ia(numdof+1)-1)
        do i=1,numdof
          do j=1,ia(i+1)-ia(i)
            ja(ia(i)-1+j)=tja(i)%d(j)
          end do
        end do
      else          ! case of LCSR if csr = .false.
        call palloc(ja,ia(numdof))
        do i=2,numdof
          do j=1,ia(i)-ia(i-1)
            ja(ia(i-1)+j)=tja(i)%d(j)
          end do
        end do
      end if
    ! Destroy tja and free memory
      lstat=destroyarrptr(tja)
!__
! 8) Write out Information about NDOFs
1     format('|  Number of degrees of freedom : ',i10)
2     format('|  Number of matrix entries : ',i12)
!-
      print "(A3)" ,'|--'
        write (*,1) numdof
      print "(A3)" ,'|--'
      if (csr) then
        print "(A3)" ,'|--'
        write (*,2) ia(numdof+1)-1
        print "(A3)" ,'|--'
      else
        print "(A3)" ,'|--'
        write (*,2) 2*ia(numdof)+numdof
        print "(A3)" ,'|--'
      end if
!__
!  
! 9) Now compact storage vectors ia, ja are known, we start with assembly
!  a) Clear matrix, right-hand-side
!  b) Nedelec elements are not available for multinature purpose
!  c) Coloring again, this time for Assembly
      if (csr) then
        if (symmetric) print*,'** MATRIXTYPE symmetric is not available with CSR-format use LCSR'
!-    a)
        call palloc(acsr,ia(numdof+1)-1)
! BUG in palloc? for temodule01, cannot use p=1
        call palloc(rhs,numdof)
        acsr=(0._DP,0._DP)
        rhs=(0._DP,0._DP)
!-    b)
        if ( eltype .ne. 'SCALAR') then
          allocate(nffsum(1))
          nffsum(1)=0
        end if
!-    c)
        ncolor=1                       !serial
!$      ncolor=size(offset)-1
!$      do color=1,ncolor              ! <-- comment out to skip the coloring
          cstart=1                     !serial
          cend=numv                    !serial
!$        cstart=offset(color)         ! <-- comment out to skip the coloring
!$        cend  =offset(color+1)-1     ! <-- comment out to skip the coloring
!
!$omp parallel do default(none) &
!$omp schedule (guided,1)       &
!$omp shared( vp, jacobi, matvar, nnat, acsr, rhs, ia, ja, vgdof)       &
!$omp shared( cstart, cend, indx, eltype)                               &
!$omp private( ai, bi, nff, nffsum, errcode, vzk, vzl, ln, kn )         &
!$omp private( i, inat, k, jnat, l, mb, elem, nffone, vpone)
!
!___  _____________________________________S T A R T__OF__A S S E M B L Y_________
!
! 10) Compute element matrixes and assemble global matrix:
!     Work through all elements (loop in parallel, range given by color/nocolor loop)
!-
!  a) Only single nature for vector bases, so vp(i,1)
!-
!  b) Loop over all natures
!     'k' iterates until size of the element matrix (number of rows)
!     'k' is the row index (local number)
!     'kn' is the row index (global number)
!     'kn' = 0 if it's DOF is not present
!     If this DOF is not present, uprade RHS (Right Hand Side)
!-
!  c) 'l'  is the local column index
!     'ln' is the global column index
!     loop over all natures once more, inat is the nature index in the direction of
!     increasing row number whereas jnat in the direction of increasing column number
!     there are nnat**2 submatrices in the global matrix
!-
!  d) This is for modified csr. so entries for diagonal are in first place.
!     set matrix entries in global diagonal
!     nffsum is used to jump to the part of the multinature element matrix that corresponds to
!     the block that we are currently looking at, if nffsum was not used then
!     we would only consider the purely nature 1 part of the element matrix
!-
!  e) Search for the column in ja
!     if this DOF was not found ==> add new entry
        do elem=cstart,cend
           i=elem                      !serial
!$        i=indx(elem)                 ! <-- comment out to skip the coloring
!
!-      a)
          if ( eltype .eq. 'SCALAR') then
            call elementmatrix_scalar3D(i,jacobi,.true.,matvar,ai,bi,nff,nffsum,errcode)
          else  ! vector bases
            vpone=vp(i,1)
            call elementmatrix3D(i,vpone, jacobi,.true.,matvar,ai,bi,nffone,errcode)
          end if
          do e = 1, size(ai,2)
            if (real(ai(e,e)).le.0) then
              print*,'ASSEMBLY ERROR: element: ', elem
              pause
            end if
          end do
!-      b)
          do inat=1,nnat
            do k=1,size(vgdof(i,inat)%d)  !nff
              kn=abs(vgdof(i,inat)%d(k))
              if (kn .eq. 0) cycle
              vzk=sign(1,vgdof(i,inat)%d(k))
              rhs(kn)=rhs(kn)+bi(k+nffsum(inat))*vzk
!-      c)
              do jnat=1,nnat
loop3:          do l=1,size(vgdof(i,jnat)%d)  !nff
!-      d)
                  if ( (k .eq. l) .and. (inat .eq. jnat) ) then 
                    acsr(ia(kn))=acsr(ia(kn))+ai(k+nffsum(inat),k+nffsum(inat))
                  else 
                    ln=abs(vgdof(i,jnat)%d(l))
                    if (ln .eq. 0) cycle
                    vzl=sign(1,vgdof(i,jnat)%d(l))
                    do mb=ia(kn)+1,ia(kn+1)-1
!-      e)
                      if (ja(mb) .eq. ln) then ! found the right column in ja
                        acsr(mb)=acsr(mb)+ai(k+nffsum(inat),l+nffsum(jnat))*vzk*vzl
                        cycle loop3
                      end if
                    end do
!$omp critical (single_print)
                    print*,'***Error in Matrix Assembly',' kn: ',kn,' ln: ',ln
!$omp end critical (single_print)
                  end if
                end do loop3
              end do
            end do
          end do
          deallocate(ai, bi)
          deallocate(nff)      !allocated in elementmatrix
          deallocate(nffsum)
!  End of elements loop
        end do
!$omp end parallel do
!
!$      end do                         ! end color loop <-- comment out to skip the coloring
!
!___  _______________________________________A S S E M B L Y_____F I N I S H E D____
!
! 11) Test the matrix
!     The matrix must have positive diagonal elements and (if symmetric) non-negative matrix entries
!     Except for boundary dofs, the sum of off-diagonal elements shoud equal the (negative) diagonal element
!     The matrix should be a M-matrix
        do i=1,numdof
          if (real(acsr(ia(i))) .le. 0.) then
            print*,'*** Systemmatrix is not a M-Matrix: may result in stability problems'
            !print*,'Processing value: ', real(acsr(ia(i))),' at position: ', i
            exit
          end if
        end do
!__
!  
      else ! NOT CSR - then it's LCSR
!__
!  
! 12) Allocate LCSR Memory
!     Clear matrix, right-hand-side
!     Nedelec elements are not available for multinature purpose
        if (symmetric) then 
          call palloc(lower,ia(numdof))
          call palloc(diag,numdof)
          call palloc(rhs,numdof)
        else
          call palloc(lower,ia(numdof))
          call palloc(upper,ia(numdof))
          call palloc(diag,numdof)
          call palloc(rhs,numdof)
          upper=(0._DP,0._DP)
        end if
        lower=(0._DP,0._DP)
        diag=(0._DP,0._DP)
        rhs=(0._DP,0._DP)
        if ( eltype .ne. 'SCALAR') then
          allocate(nffsum(1))
          nffsum(1)=0
        end if
!__
!  
! 13) Coloring again, this time for Assembly
        ncolor=1                       !serial
!$      ncolor=size(offset)-1
!$      do color=1,ncolor              ! <-- comment out to skip the coloring
          cstart=1                     !serial
          cend=numv                    !serial
!$        cstart=offset(color)         ! <-- comment out to skip the coloring
!$        cend  =offset(color+1)-1     ! <-- comment out to skip the coloring
        ncolor=1                !serial
!
!$omp parallel do default(none) &
!$omp schedule (guided,1)       &
!$omp shared(vp, jacobi, matvar, nnat, acsr, rhs, ia, ja, vgdof)        &
!$omp shared(diag, lower, upper, symmetric, cstart, cend, indx, eltype) &
!$omp private( ai, bi, nff, nffsum, errcode, vzk, vzl, ln, kn )         &
!$omp private( i, inat, k, jnat, l, mb, elem, nffone, vpone)            &
!$omp private( iloop4, lnn, knn, l1, k1)
!
!___  _____________________________________S T A R T__OF__A S S E M B L Y_________
!
! 14) Compute element matrixes and assemble global matrix:
!     Work through all elements (loop in parallel, range given by color/nocolor loop)
!-
!  a) Compute element matrix
!-
!  b) Only single nature for vector bases, so vp(i,1)
        do elem=cstart,cend
          i=elem            !serial
!$        i=indx(elem)                                                         ! <-- comment out to skip the coloring
!
!-      a)
          if ( eltype .eq. 'SCALAR') then
            call elementmatrix_scalar3D(i,jacobi,.true.,matvar,ai,bi,nff,nffsum,errcode)
          else
!-      b)
            vpone=vp(i,1)
            call elementmatrix3D(i,vpone, jacobi,.true.,matvar,ai,bi,nffone,errcode)
          end if
!  loop over all natures
          do inat=1,nnat
            do k=1,nff(inat)
!  row index k: local;  kn: global
              kn=abs(vgdof(i,inat)%d(k))
!  if this dof is not present
              if (kn .eq. 0) cycle
              vzk=sign(1,vgdof(i,inat)%d(k))
!  update rhs
              rhs(kn)=rhs(kn)+bi(k+nffsum(inat))*vzk
!  update diagonal
              diag(kn)=diag(kn)+ai(k+nffsum(inat),k+nffsum(inat))
              do jnat=1,nnat
!  consider local diagonal entries only for the matrix blocks lying above
!  the global diagonal
                if (inat .eq. jnat) then 
                  iloop4=k-1
                else if (inat .gt. jnat) then
                  iloop4=nff(jnat)
                else
                  iloop4=0
                end if
loop4:          do l=1,iloop4
!  column index l: local;  ln: global
                  ln=abs(vgdof(i,jnat)%d(l))
                  if (ln .eq. 0) cycle
                  vzl=sign(1,vgdof(i,jnat)%d(l))
!
                  if (kn .lt. ln) then
!  swap row and column -  column index must be smaller than row index
                    lnn=kn
                    knn=ln
!  l1 and k1 are used for the indices of the multinature element matrix so that
!  we pick the correct value from it to contribute to the corresponding entry in
!  the global matrix
                    l1=k+nffsum(inat)
                    k1=l+nffsum(jnat)
                  else
                    lnn=ln
                    knn=kn
                    l1=l+nffsum(jnat)
                    k1=k+nffsum(inat)
                  end if
!
                  do mb=ia(knn-1)+1,ia(knn)
!  search for the the column in ja in the list
                    if (ja(mb) .eq. lnn) then
!  found
                      lower(mb)=lower(mb)+ai(k1,l1)*vzk*vzl
!  if the matrix is known to be symmetric we do not need to assign upper matrix
                      if (.not.symmetric) upper(mb)=upper(mb)+ai(l1,k1)*vzk*vzl
                      cycle loop4
                    end if
                  end do
!$omp critical (single_print)
                  print*,'***Error in Matrix Assembly',' kn: ',kn,' ln: ',ln,' inat: ',inat,' jnat: ',jnat
!$omp end critical (single_print)
                end do loop4
             end do ! jnat
          end do ! k local dofs
       end do ! inat
          deallocate(ai, bi)
          deallocate(nff) ! arrays allocated in elementmatrix
          deallocate(nffsum)
!  End of elements loop
        end do
!        print *, 'RHSUPDATE', real(rhs(121))

!$omp end parallel do
!
!$      end do                                                                 ! end color loop <-- comment out to skip the coloring
!
!  ASSEMBLY DONE
!
!  Test the matrix
        count=0
        do i=1,numdof
!  The matrix must have positive diagonal elements and (if symmetric) non-negative matrix entries
!
!  Except for boundary dofs, the sum of off-diagonal elements shoud equal the (negative) diagonal element
!  The matrix should be a M-matrix
          if (real(diag(i)) .le. 0.) then
            count = count + 1
            if (count .eq. 1) print*,'*** Systemmatrix is not a M-Matrix: may result in stability problems'
            exit
          end if
        end do
      end if
!
!
      if (allocated(indx)) deallocate(indx,stat=i)
      if (allocated(offset)) deallocate(offset,stat=i)
!
      call zeit(' Matrix Assembly')
!
      return
      end subroutine assembly3D
