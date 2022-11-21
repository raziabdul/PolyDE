      subroutine assemblytrans(lower,upper,diag,acsrs,acsrd,acsrm,rhs,ia,ja, &
      &                        jacobi,csr,matvar,symmetric)
      use feminterface, only: elementmatrixtrans, zeit, reallocate, destroyarrptr, coloring
      use femtypes
      use globalvariables, only: n, ndof, nnat, ep, eg
      implicit none
      integer (I4B), pointer :: ia(:), ja(:)
      complex (DPC), pointer :: lower(:), upper(:), diag(:), acsrs(:), acsrd(:), acsrm(:), rhs(:)
      logical jacobi, csr, matvar, symmetric
      intent (in) ::  jacobi, csr, matvar, symmetric
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
!    $Revision: 1.10 $
!    $Date: 2014/05/13 13:17:46 $
!    $Author: m_kasper $
!
!  Element by element FEM matrix and right-hand-side assembly in compact storage
!
!  Input:  
!                     / = .true. for the computation of the jacobian
!            jacobi
!                     \ = .false. for the linear case
!            n        total number of elements
!            eg       degrees of freedom (global) for the elements for all natures
!                     packed into an array of pointers   eg(i,inat)%d(j)
!            ep       polynomial degree for each of the elements and for all natures
!            ndof     total number of degrees of freedom
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
!            acsrs    vector (compact storage) of non-zero stiffness matrix entries
!            acsrd    vector (compact storage) of non-zero damping   matrix entries
!            acsrm    vector (compact storage) of non-zero mass      matrix entries
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
!  A vector of indices of the entries in A:           ja=( 1, 4, 2, 1,  2, 4,  3, 3,  4, 1)
!  The vector of pointers to the first entry of rows: ia=( 1, 4, 6, 8, 11)
!
!  Actually we use a slightly modified CSR. The difference to the standard one lies in the
!  accessibility of the diagonal entries. They can be directly accessed by the values given
!  in vector ia in this version. So the diagonal entries for each row are always in first
!  place (referred to by ia).
!  example: the diagonal matrix entry no. A(3,3) is addressed by ia(3). So it can be found
!           at acsr(6) - 6 is the value given in ia(3).
!
!  The vector or non-zero matrix entries:           acsr=(-2, 4, 1, -8, 8, -6, 6, -5, 3, 5)
!  A vector of indices of the entries in A:           ja=( 1, 4, 2,  2, 1,  3, 4,  4, 1, 3)
!  The vector of pointers to the first entry of rows: ia=( 1, 4, 6, 8, 11)
!
!
!  The compressed lower sparse row format (LCSR) take advantage of the fact that in FEM the 
!  occupation pattern is symmetric i.e. if A(i,j)=0 then A(j,i)=0. This format uses 
!  three vectors to hold the non-zero matrix entries and two index vectors
!  The lower-triangular matrix                         L=( 8, 5, 3)
!  The transposed upper triangular matrix              U=( 1, 6, 4)
!  The diagonal:                                       D=(-2, -8, -6, -5)
!  A vector of indices of the entries in L:           ja=( 1, 3, 1)
!  The vector of pointers to the last entry of rows:  ia=( 0, 1, 1, 3)
!  In this format the vector n2l is significatly smaller (roughly half the size of CSR)
!
!  Note that ia has different interpretation and length for the two storage schemes!
!
!  local variables
!  indices iloop2 and iloop4 take care of the multinature LCSR assembly since taking advantage of
!  multinature symmetric matrices is tricky
      integer (I4B) :: i, j, k, l, kn, knn, ln, lnn, k1, l1, mb, acsize, iloop2, iloop4
      integer (I4B) errcode, vzk, vzl
      integer (I4B), allocatable :: nff(:), nffsum(:)
      complex (DPC), pointer :: si(:,:), di(:,:), mi(:,:), bi(:)
      type (ARRPTRI), pointer :: tja(:)
      logical :: lstat
!  inat stands for index of nature, jnat for a nested nature index
      integer (I4B) inat, jnat
      
!  parameter isize determines the initial size to hold one matrix row
!  it should be chosen a bit larger than the expected mean number of row entries 
!  the minimum value (not recommended) is 2 
!  small values cause excessive reallocation and possibly memory fragmentation
!  large values may (temporarily) allocate much memory
      integer (I4B), parameter :: isize=20
!
      integer (I4B), allocatable :: indx(:), offset(:)
      integer (I4B) color, cstart, cend, elem, ncolor
!
!  initialize temporary storage (vector of pointers)
!  tja(ndof) is a vector where every (say) subvector holds the column indices
!  of the matrix entries for one row
      allocate(tja(ndof))
      if (csr) then
!  case of CSR format
        allocate(ia(ndof+1))
        ia=1
        do i=1,ndof
!  begin with an allocation of isize spaces for every tja subvector
          allocate(tja(i)%d(isize))
!  this is for modified csr, entries for diagonal are in first place
          tja(i)%d(1)=i
          tja(i)%d(2:isize)=0
        end do
      else
!  case of LCSR format if csr = .false.
        allocate(ia(ndof))
        ia=0
        do i=1,ndof
          allocate(tja(i)%d(isize/2))
          tja(i)%d=0
        end do
      end if
!
!  determine number of row entries and column indices 
      if (csr) then
!  case of CSR format
!  loop over all elements
        do i=1,n
!  loop over all natures
          do inat=1,nnat
!  size (number of rows) of the element matrix
            do k=1,size(eg(i,inat)%d)
!  row index k: local;  kn: global
              kn=abs(eg(i,inat)%d(k))
!  if this dof is not present
              if (kn .eq. 0) cycle
!  the actual size to hold this row
              acsize=size(tja(kn)%d)
!  column index l: local;  ln: global
!  loop over all natures once more, inat is the nature index in the direction of
!  increasing row number whereas jnat in the direction of increasing column number
!  there are nnat**2 submatrices in the global matrix
              do jnat=1,nnat
loop1:          do l=1,size(eg(i,jnat)%d)
!  this is for modified csr, entries for diagonal are in first place
!  As global diagonal entries we consider only the local diagonal entries of the
!  say diagonal matrix blocks, in the multinature global matrix we have
!  purely nature parts and coupled nature parts
!  the coupled nature parts are located off the global diagonal
                  if ( (k .ne. l) .or. (inat .ne. jnat) ) then
                    ln=abs(eg(i,jnat)%d(l))
                    if (ln .eq. 0) cycle
!
                    do mb=2,acsize
!  search for the column in ja, if this dof was processed previously it is in the list
!  otherwise it is stored at the first free position (entry = 0)
                      if (tja(kn)%d(mb) .eq. 0) then
!  not found ==> add new entry
                        tja(kn)%d(mb)=ln
                        ia(kn)=ia(kn)+1
                        cycle loop1
                      else if (tja(kn)%d(mb) .eq. ln) then
!  found
                        cycle loop1
                      end if
                    end do
!  not found no more space ==> resize
!  if there is no more space in every tja subvector to hold the column indices
!  of the entries of one row in the global matrix then allocate double the space
!  for storage to this tja subvector
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
      else
!  case of LCSR if csr = .false.
!  loop over all elements
        do i=1,n
!  loop over all natures
          do inat=1,nnat
!  size (number of rows) of the element matrix
            do k=1,size(eg(i,inat)%d)
!  row index k: local;  kn: global
              kn=abs(eg(i,inat)%d(k))
!  if this dof is not present
              if (kn .eq. 0) cycle
!  column index l: local;  ln: global
              do jnat=1,nnat
!  for loop2 located a few lines below in the code we count up to the value iloop2
!  this is because we want to consider the local diagonal entries only for the
!  submatrix blocks where index inat is smaller than jnat
!  this submatrix block lies above the diagonal, the local diagonal entries of the block
!  lying below the global diagonal will be taken care of by swaping the indices
!  of the symmetric entries which lie above the global diagonal
                if (inat .eq. jnat) then
                  iloop2=k-1
                else  if (inat .gt. jnat) then
                  iloop2=size(eg(i,jnat)%d)
                else
                  iloop2=0
                end if

loop2:          do l=1,iloop2
                  ln=abs(eg(i,jnat)%d(l))
                  if (ln .eq. 0) cycle
!  if global row number is smaller than the global column number then
                  if (kn .lt. ln) then
!  swap row and column - column index must be smaller than row index
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
          end do
        end do
      end if
!  now we have the number of row entries stored in ia 
!  sum up to get ia in the CSR or LCSR format
      if (csr) then
        do i = ndof, 1, -1
          ia(i+1) = ia(i)
        end do
        ia(1)=1
        do i = 2, ndof+1
          ia(i) = ia(i) + ia(i-1)
        end do
      else
        do i = 2, ndof
          ia(i) = ia(i) + ia(i-1)
        end do
      end if
!  generate ja (compact) from temporary tja (vector of pointers)
      if (csr) then
!  case of CSR if csr = .true.
        allocate(ja(1:ia(ndof+1)-1))
        do i=1,ndof
          do j=1,ia(i+1)-ia(i)
            ja(ia(i)-1+j)=tja(i)%d(j)
          end do
        end do
      else
!  case of LCSR if csr = .false.
        allocate(ja(1:ia(ndof)))
        do i=2,ndof
          do j=1,ia(i)-ia(i-1)
            ja(ia(i-1)+j)=tja(i)%d(j)
          end do
        end do
      end if
! destroy tja and free memory
      lstat=destroyarrptr(tja)
!
      write (*,1) ndof
1     format('   Number of degrees of freedom : ',i10)
      if (csr) then
        write (*,2) ia(ndof+1)-1
      else 
        write (*,2) 2*ia(ndof)+ndof
      end if
2     format('   Number of matrix entries : ',i12)
!
! NOTES of parallel processing:
!
! Observe the OpenMP conditional compilation marker '!$' 
! -> Current code is using coloring.
! -> Lines must be commented out to remove the coloring loop and variables.
! -> Critical sections must be reintroduced if coloring is not used. 
!
! pointers indx and offset are allocated in coloring
!
!
!  now compact storage vectors ia, ja are known, we start with assembly
      if (csr) then
        if (symmetric) print*,'** MATRIXTYPE symmetric is not ', &
        &                     'available with CSR-format use LCSR'
!  case of CSR format
!  clear matrix, right-hand-side
        allocate(acsrs(1:ia(ndof+1)-1), acsrd(1:ia(ndof+1)-1), acsrm(1:ia(ndof+1)-1), rhs(ndof))
        acsrs=(0._DP,0._DP)
        acsrd=(0._DP,0._DP)
        acsrm=(0._DP,0._DP)
        rhs=(0._DP,0._DP)
!
        ncolor=1                !serial
!
! Start color loop
          cstart=1              !serial
          cend=n                !serial
!
!  loop over all elements (loop in parallel, range given by color/nocolor loop)
        do elem=cstart,cend
!
          i=elem            !serial
!
!  compute element matrix
!  element matrix to be computed for all natures, argument ep(i,:)
          call elementmatrixtrans(i,ep(i,:),jacobi,.true.,matvar,si,di,mi,bi,nff,nffsum,errcode)
!  loop over all natures
          do inat=1,nnat
            do k=1,size(eg(i,inat)%d)
!  row index k: local;  kn: global
              kn=abs(eg(i,inat)%d(k))
!  if this dof is not present
              if (kn .eq. 0) cycle
              vzk=sign(1,eg(i,inat)%d(k))
!  update rhs
              rhs(kn)=rhs(kn)+bi(k+nffsum(inat))*vzk
!  loop over all natures once more
              do jnat=1,nnat
!  column index l: local;  ln: global
loop3:          do l=1,size(eg(i,jnat)%d)
!  this is for modified csr. so entries for diagonal are in first place.
!  set matrix entries in global diagonal
!  nffsum is used to jump to the part of the multinature element matrix that corresponds to
!  the block that we are currently looking at, if nffsum was not used then
!  we would only consider the purely nature 1 part of the element matrix
                  if ( (k .eq. l) .and. (inat .eq. jnat) ) then 
                    acsrs(ia(kn))=acsrs(ia(kn))+si(k+nffsum(inat),k+nffsum(inat))
                    acsrd(ia(kn))=acsrd(ia(kn))+di(k+nffsum(inat),k+nffsum(inat))
                    acsrm(ia(kn))=acsrm(ia(kn))+mi(k+nffsum(inat),k+nffsum(inat))
                  else
                    ln=abs(eg(i,jnat)%d(l))
                    if (ln .eq. 0) cycle
                    vzl=sign(1,eg(i,jnat)%d(l))
!
                    do mb=ia(kn)+1,ia(kn+1)-1
!  search for the column in ja in the list
                      if (ja(mb) .eq. ln) then
!  found
                        acsrs(mb)=acsrs(mb)+si(k+nffsum(inat),l+nffsum(jnat))*vzk*vzl
                        acsrd(mb)=acsrd(mb)+di(k+nffsum(inat),l+nffsum(jnat))*vzk*vzl
                        acsrm(mb)=acsrm(mb)+mi(k+nffsum(inat),l+nffsum(jnat))*vzk*vzl
                        cycle loop3
                      end if
                    end do
                    print*,'***Error in Matrix Assembly',' kn: ',kn,' ln: ',ln
                  end if
                end do loop3
              end do
            end do
          end do
          deallocate(si, di, mi, bi)
          if (allocated(nff)) deallocate(nff)      !allocated in elementmatrix
          if (allocated(nffsum)) deallocate(nffsum)
!  End of elements loop
        end do
!
!
!  Test the matrix
        do i=1,ndof
!  The matrix must have positive diagonal elements and (if symmetric) non-negative matrix entries
!
!  Except for boundary dofs, the sum of off-diagonal elements shoud equal the (negative) diagonal element
!  The matrix should be a M-matrix
          if (real(acsrs(ia(i))) .le. 0.) then
            print*,'*** Systemmatrix is not a M-Matrix: may result in stability problems'
            exit
          end if
        end do
!
      else
!  case of LCSR if csr = .false.
        if (symmetric) then
          allocate(lower(1:ia(ndof)), diag(ndof), rhs(ndof))
        else
          allocate(lower(1:ia(ndof)), upper(1:ia(ndof)), diag(ndof), rhs(ndof))
!  clear matrix, right-hand-side
          upper=(0._DP,0._DP)
        end if
        lower=(0._DP,0._DP)
        diag=(0._DP,0._DP)
        rhs=(0._DP,0._DP)
!
        ncolor=1                !serial
!
! Start color loop
          cstart=1              !serial
          cend=n                !serial
!
!  loop over all elements (it should be possible to execute this loop in parallel)
        do elem=cstart,cend
          i=elem            !serial

!  compute element matrix
          call elementmatrixtrans(i,ep(i,:),jacobi,.true.,matvar,si,di,mi,bi,nff,nffsum,errcode)
!  loop over all natures
          do inat=1,nnat
            do k=1,nff(inat)
!  row index k: local;  kn: global
              kn=abs(eg(i,inat)%d(k))
!  if this dof is not present
              if (kn .eq. 0) cycle
              vzk=sign(1,eg(i,inat)%d(k))
!  update rhs
              rhs(kn)=rhs(kn)+bi(k+nffsum(inat))*vzk
!  update diagonal
              diag(kn)=diag(kn)+si(k+nffsum(inat),k+nffsum(inat))
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
                  ln=abs(eg(i,jnat)%d(l))
                  if (ln .eq. 0) cycle
                  vzl=sign(1,eg(i,jnat)%d(l))
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
                      lower(mb)=lower(mb)+si(k1,l1)*vzk*vzl
!  if the matrix is known to be symmetric we do not need to assign upper matrix
                      if (.not.symmetric) then
                          upper(mb)=upper(mb)+si(l1,k1)*vzk*vzl
                      endif
                      cycle loop4
                    end if
                  end do
                  print*,'***Error in Matrix Assembly',' kn: ',kn,' ln: ',ln
                end do loop4
              end do
            end do
          end do
          deallocate(si, di, mi, bi)
          if (allocated(nff)) deallocate(nff) !allocated in elementmatrix
          if (allocated(nffsum)) deallocate(nffsum)
!  End of elements loop
        end do
!
!
!  Test the matrix
        do i=1,ndof
!  The matrix must have positive diagonal elements and (if symmetric) non-negative matrix entries
!
!  Except for boundary dofs, the sum of off-diagonal elements shoud equal the (negative) diagonal element
!  The matrix should be a M-matrix
          if (real(diag(i)) .le. 0.) then
            print*,'*** Systemmatrix is not a M-Matrix: may result in stability problems'
            exit
          end if
        end do
      end if
!
      if (allocated(indx)) deallocate(indx)
      if (allocated(offset)) deallocate(offset)
!
      call zeit(' Matrix Assembly')
!
      return
      end subroutine assemblytrans
