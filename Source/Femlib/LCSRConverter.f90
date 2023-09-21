      subroutine sort_asc_order_DP(a, ia, ja)
      use femtypes
      use feminterface, only: qsortindex
      implicit none
      real (DP), pointer:: a(:)
      integer (I4B), pointer :: ia(:), ja(:)
      intent(in) :: ia
      intent (inout) :: a, ja
!  sort a, ja such that
!       column indices (if CSR format) or
!       row indices (if CSC format) 
!  appear in ascending order

!  local variables
      integer(I4B) :: i, j, n, num, maxlen
      integer(I4B), allocatable :: indx(:), jatemp(:)
      real(DP), allocatable :: atemp(:)

      n = size(ia) - 1
      maxlen = 0
      do i = 1, n
        maxlen = max(maxlen, ia(i+1) - ia(i))
      end do
      allocate(indx(maxlen), atemp(maxlen), jatemp(maxlen))

      do i = 1, n
        num = ia(i+1) - ia(i)
        jatemp(1:num) = ja(ia(i):ia(i)+num-1)
        atemp(1:num)  = a(ia(i):ia(i)+num-1)
        call qsortindex(jatemp,indx,num)
        do j = 1, num
          ja(ia(i)+num-j) = jatemp(indx(j))
          a (ia(i)+num-j) = atemp (indx(j))
        end do
      end do

      return
      end

      
      subroutine sort_asc_order_DPC(a, ia, ja)
      use femtypes
      use feminterface, only: qsortindex
      implicit none
      complex (DPC), pointer:: a(:)
      integer (I4B), pointer :: ia(:), ja(:)
      intent(in) :: ia
      intent (inout) :: a, ja
!  sort a, ja such that
!       column indices (if CSR format) or
!       row indices (if CSC format) 
!  appear in ascending order

!  local variables
      integer(I4B) :: i, j, n, num, maxlen
      integer(I4B), allocatable :: indx(:), jatemp(:)
      complex(DPC), allocatable :: atemp(:)

      n = size(ia) - 1
      maxlen = 0
      do i = 1, n
        maxlen = max(maxlen, ia(i+1) - ia(i))
      end do
      allocate(indx(maxlen), atemp(maxlen), jatemp(maxlen))

      do i = 1, n
        num = ia(i+1) - ia(i)
        jatemp(1:num) = ja(ia(i):ia(i)+num-1)
        atemp(1:num)  = a(ia(i):ia(i)+num-1)
        call qsortindex(jatemp,indx,num)
        do j = 1, num
          ja(ia(i)+num-j) = jatemp(indx(j))
          a (ia(i)+num-j) = atemp (indx(j))
        end do
      end do

      return
      end



!  Convertes from LCSR to COO
      subroutine lcsr2coo_DP(diag, lower, upper, symm, lia, lja, ao, ir, jc)
      use femtypes
      implicit none
      real (DP), pointer:: diag(:), lower(:), upper(:), ao(:)
      integer (I4B), pointer :: lia(:), lja(:), ir(:), jc(:)
      logical symm
      intent(in) :: symm
      intent (out) :: ao, ir, jc
      intent (inout) :: diag, lower, upper, lia, lja
!  Convert LCSR sparse matrix format to COO format

!  Input (LCSR)
!            diag     vector of diagonal matrix enties
!            lower    vector of matrix enties in lower triangle matrix
!            upper    vector of matrix enties in lower triangle matrix
!            sym      = .true. if matrix is symmetric
!            lia      is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            lja      vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix
!
!  Output
!            acsr     COO matrix elements
!            ia       COO row indices of matrix entries
!            ja       COO column indices of matrix entries

!            ao       matrix in coo format
!            ir       row indices
!            jc       column indices

!  local variables
      integer (I4B) :: row, k, n, lnnz, nnz

      n = size(diag)                                        ! number of rows
      lnnz = size(lower)                                    ! number of non-zero entries in lower (or upper) matrix
      nnz = n + 2*lnnz                                      ! number of non-zero entries
      allocate(ao(nnz))

!  copy matrix elements
      ao(1:n)               = diag(1:n)
      ao(n+1:n+lnnz)        = lower(1:lnnz)
      ao(n+lnnz+1:n+2*lnnz) = upper(1:lnnz)
      if (symm) then
        deallocate (diag, lower)
      else
        deallocate (diag, lower, upper)
      end if
      nullify (diag, lower, upper)

!  copy indices
!    diagonal
      allocate(ir(nnz), jc(nnz))
      do row = 1, n
        ir(row) = row
        jc(row) = row
      end do

      do row = 2, n
        do k = lia(row-1)+1, lia(row)
!    lower triangle
          ir(n+k) = row
          jc(n+k) = lja(k)
!    upper triangle (row / col are swapped)
          ir(n+lnnz+k) = lja(k)
          jc(n+lnnz+k) = row
        end do
      end do

      deallocate (lia, lja)
      nullify (lia, lja)

      return
      end subroutine lcsr2coo_DP



      subroutine lcsr2coo_SP(diag, lower, upper, symm, lia, lja, ao, ir, jc)
      use femtypes
      implicit none
      real (SP), pointer:: diag(:), lower(:), upper(:), ao(:)
      integer (I4B), pointer :: lia(:), lja(:), ir(:), jc(:)
      logical symm
      intent(in) :: symm
      intent (out) :: ao, ir, jc
      intent (inout) :: diag, lower, upper, lia, lja
!  Convert LCSR sparse matrix format to COO format

!  Input (LCSR)
!            diag     vector of diagonal matrix enties
!            lower    vector of matrix enties in lower triangle matrix
!            upper    vector of matrix enties in lower triangle matrix
!            sym      = .true. if matrix is symmetric
!            lia      is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            lja      vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix
!
!  Output
!            acsr     COO matrix elements
!            ia       COO row indices of matrix entries
!            ja       COO column indices of matrix entries

!            ao       matrix in coo format
!            ir       row indices
!            jc       column indices

!  local variables
      integer (I4B) :: row, k, n, lnnz, nnz

      n = size(diag)                                        ! number of rows
      lnnz = size(lower)                                    ! number of non-zero entries in lower (or upper) matrix
      nnz = n + 2*lnnz                                      ! number of non-zero entries
      allocate(ao(nnz))

!  copy matrix elements
      ao(1:n)               = diag(1:n)
      ao(n+1:n+lnnz)        = lower(1:lnnz)
      ao(n+lnnz+1:n+2*lnnz) = upper(1:lnnz)
      if (symm) then
        deallocate (diag, lower)
      else
        deallocate (diag, lower, upper)
      end if
      nullify (diag, lower, upper)

!  copy indices
!    diagonal
      allocate(ir(nnz), jc(nnz))
      do row = 1, n
        ir(row) = row
        jc(row) = row
      end do

      do row = 2, n
        do k = lia(row-1)+1, lia(row)
!    lower triangle
          ir(n+k) = row
          jc(n+k) = lja(k)
!    upper triangle (row / col are swapped)
          ir(n+lnnz+k) = lja(k)
          jc(n+lnnz+k) = row
        end do
      end do

      deallocate (lia, lja)
      nullify (lia, lja)

      return
      end subroutine lcsr2coo_SP



      subroutine lcsr2coo_DPC(diag, lower, upper, symm, lia, lja, ao, ir, jc)
      use femtypes
      implicit none
      complex (DPC), pointer:: diag(:), lower(:), upper(:), ao(:)
      integer (I4B), pointer :: lia(:), lja(:), ir(:), jc(:)
      logical symm
      intent(in) :: symm
      intent (out) :: ao, ir, jc
      intent (inout) :: diag, lower, upper, lia, lja
!  Convert LCSR sparse matrix format to COO format

!  Input (LCSR)
!            diag     vector of diagonal matrix enties
!            lower    vector of matrix enties in lower triangle matrix
!            upper    vector of matrix enties in lower triangle matrix
!            sym      = .true. if matrix is symmetric
!            lia      is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            lja      vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix
!
!  Output
!            acsr     COO matrix elements
!            ia       COO row indices of matrix entries
!            ja       COO column indices of matrix entries

!            ao       matrix in coo format
!            ir       row indices
!            jc       column indices

!  local variables
      integer (I4B) :: row, k, n, lnnz, nnz

      n = size(diag)                                        ! number of rows
      lnnz = size(lower)                                    ! number of non-zero entries in lower (or upper) matrix
      nnz = n + 2*lnnz                                      ! number of non-zero entries
      allocate(ao(nnz))

!  copy matrix elements
      ao(1:n)               = diag(1:n)
      ao(n+1:n+lnnz)        = lower(1:lnnz)
      ao(n+lnnz+1:n+2*lnnz) = upper(1:lnnz)
      if (symm) then
        deallocate (diag, lower)
      else
        deallocate (diag, lower, upper)
      end if
      nullify (diag, lower, upper)

!  copy indices
!    diagonal
      allocate(ir(nnz), jc(nnz))
      do row = 1, n
        ir(row) = row
        jc(row) = row
      end do

      do row = 2, n
        do k = lia(row-1)+1, lia(row)
!    lower triangle
          ir(n+k) = row
          jc(n+k) = lja(k)
!    upper triangle (row / col are swapped)
          ir(n+lnnz+k) = lja(k)
          jc(n+lnnz+k) = row
        end do
      end do

      deallocate (lia, lja)
      nullify (lia, lja)

      return
      end subroutine lcsr2coo_DPC


!  Convertes from LCSR to CSR
      subroutine lcsr2csr_DP(diag, lower, upper, symm, lia, lja, acsr, ia, ja)
      use femtypes
      implicit none
      real(DP), pointer :: diag(:), lower(:), upper(:), acsr(:)
      integer(I4B), pointer :: lia(:), lja(:), ja(:), ia(:)
      logical symm
      intent(in) :: symm
      intent(out) :: acsr, ia, ja
      intent(inout) :: diag, lower, upper, lia, lja
!  Convert LCSR sparse matrix format to CSR format

!  Input (LCSR)
!            diag     vector of diagonal matrix enties
!            lower    vector of matrix enties in lower triangle matrix
!            upper    vector of matrix enties in lower triangle matrix
!            sym      = .true. if matrix is symmetric
!            lia      is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            lja      vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix
!
!  Output
!            acsr     CSR matrix elements
!            ia       CSR accumulated sum of the number of non-zero elements in matrix
!            ja       CSR column index of matrix entries

!  local variables
      integer(I4B) :: k, n, lnnz, nnz, az, row, col
      integer(I4B), allocatable :: lastfree(:)

!  Build ia
      print *, 'LCSR2CSR DP BUGGY?'

      n = size(diag)                                        ! number of rows
      lnnz = size(lower)                                    ! number of non-zero entries in lower (or upper) matrix
      nnz = n + 2*lnnz                                      ! number of non-zero entries

      allocate (ia(n+1))
      ia = 0
!  the loop passes array lja and counts matrix elements of value i
!  this will give the number of entries in each row resulting from upper
      do k = 1, lnnz
        ia(lja(k)+1) = ia(lja(k)+1) + 1
      end do
!  build ia from (per row) entiers in lower upper and 1 (for diag)
!  (az equals the sum of entries in upper)
      az = 0
      ia(1) = 1
      do row=2,n
        az = az + ia(row)
        ia(row) = row + lia(row-1) + az
      end do
      ia(n+1) = nnz + 1

      allocate (acsr(nnz), ja(nnz))


!  Copy diagonal elements (in first position)
      do row = 1, n
        ja(ia(row)) = row
        acsr(ia(row)) = diag(row)
      end do

      deallocate (diag)
      nullify(diag)

!  copy entries from lower
      do row = 2, n
        do k = 1, lia(row)-lia(row-1)
          col = lja(k+lia(row-1))
          acsr(ia(row)+k) = lower(k+lia(row-1))
          ja(ia(row)+k) = col
        end do
      end do

      if (.not. symm) then
        deallocate (lower)
        nullify(lower)
      else
        upper => lower
      end if

      allocate (lastfree(n-1))

      do row = 1, n-1
        lastfree(row) = ia(row+1) - 1
      end do
    

!  copy entries from upper
      do col = 2, n
        do k = lia(col-1)+1, lia(col)
          row = lja(k)
          acsr(lastfree(row)) = upper(k)
          ja(lastfree(row)) = col
          lastfree(row) = lastfree(row) - 1 
        end do
      end do

      deallocate (lastfree)
      if (symm) then
        deallocate (lower)
        nullify(lower)
      else
        deallocate (upper)
      end if
      nullify(upper)
      deallocate (lia, lja)
      nullify(lia, lja)
      return
      end subroutine lcsr2csr_DP



      subroutine lcsr2csr_DPC(diag, lower, upper, symm, lia, lja, acsr, ia, ja)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag(:), lower(:), upper(:), acsr(:)
      integer(I4B), pointer :: lia(:), lja(:), ja(:), ia(:)
      logical symm
      intent(in) :: symm
      intent(out) :: acsr, ia, ja
      intent(inout) :: diag, lower, upper, lia, lja
!  Convert LCSR sparse matrix format to CSR format

!  Input (LCSR)
!            diag     vector of diagonal matrix enties
!            lower    vector of matrix enties in lower triangle matrix
!            upper    vector of matrix enties in lower triangle matrix
!            sym      = .true. if matrix is symmetric
!            lia      is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            lja      vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix
!
!  Output
!            acsr     CSR matrix elements
!            ia       CSR accumulated sum of the number of non-zero elements in matrix
!            ja       CSR column index of matrix entries

!  local variables
      integer(I4B) :: k, n, lnnz, nnz, az, row, col
      integer(I4B), allocatable :: lastfree(:)

!  Build ia
      n = size(diag)                                        ! number of rows
      lnnz = size(lower)                                    ! number of non-zero entries in lower (or upper) matrix
      nnz = n + 2*lnnz                                      ! number of non-zero entries

      allocate (ia(n+1))
      ia = 0
!  the loop passes array lja and counts matrix elements of value i
!  this will give the number of entries in each row resulting from upper
      do k = 1, lnnz
        ia(lja(k)+1) = ia(lja(k)+1) + 1
      end do

!  build ia from (per row) entiers in lower upper and 1 (for diag)
!  (az equals the sum of entries in upper)

      az = 0
      ia(1) = 1
      do row=2,n
        az = az + ia(row)
        ia(row) = row + lia(row-1) + az
      end do
      ia(n+1) = nnz + 1

      allocate (acsr(nnz), ja(nnz))

!  Copy diagonal elements (in first position)
      do row = 1, n
        ja(ia(row)) = row
        acsr(ia(row)) = diag(row)
      end do

      deallocate (diag)
      nullify(diag)

!  copy entries from lower
      do row = 2, n
        do k = 1, lia(row)-lia(row-1)
          col = lja(k+lia(row-1))
          acsr(ia(row)+k) = lower(k+lia(row-1))
          ja(ia(row)+k) = col
        end do
      end do

      if (.not. symm) then
        deallocate (lower)
        nullify(lower)
      else
        upper => lower
      end if

      allocate (lastfree(n-1))

      do row = 1, n-1
        lastfree(row) = ia(row+1) - 1
      end do

!  copy entries from upper
      do col = 2, n
        do k = lia(col-1)+1, lia(col)
          row = lja(k)
          acsr(lastfree(row)) = upper(k)
          ja(lastfree(row)) = col
          lastfree(row) = lastfree(row) - 1 
        end do
      end do

      deallocate (lastfree)
      if (symm) then
        deallocate (lower)
        nullify(lower)
      else
        deallocate (upper)
      end if
      nullify(upper)
      deallocate (lia, lja)
      nullify(lia, lja)

      return
      end subroutine lcsr2csr_DPC


!  Convertes from LCSR to CSC
      subroutine lcsr2csc_DP(diag, lower, upper, symm, lia, lja, acsr, ia, ja)
      use femtypes
      implicit none
      real(DP), pointer :: diag(:), lower(:), upper(:), acsr(:)
      integer(I4B), pointer :: lia(:), lja(:), ja(:), ia(:)
      logical symm
      intent(in) :: symm
      intent(out) :: acsr, ia, ja
      intent(inout) :: diag, lower, upper, lia, lja
!  Convert LCSR sparse matrix format to CSC format

!  Input (LCSR)
!            diag     vector of diagonal matrix enties
!            lower    vector of matrix enties in lower triangle matrix
!            upper    vector of matrix enties in lower triangle matrix
!            sym      = .true. if matrix is symmetric
!            lia      is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            lja      vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix
!
!  Output
!            acsr     CSR matrix elements
!            ia       CSR accumulated sum of the number of non-zero elements in matrix
!            ja       CSR column index of matrix entries

!  local variables
      integer(I4B) :: k, n, lnnz, nnz, az, row, col
      integer(I4B), allocatable :: lastfree(:)

!  Build ia

      n = size(diag)                                        ! number of rows
      lnnz = size(lower)                                    ! number of non-zero entries in lower (or upper) matrix
      nnz = n + 2*lnnz                                      ! number of non-zero entries

      allocate (ia(n+1))
      ia = 0
!  the loop passes array lja and counts matrix elements of value i
!  this will give the number of entries in each row resulting from upper
      do k = 1, lnnz
        ia(lja(k)+1) = ia(lja(k)+1) + 1
      end do

!  build ia from (per row) entiers in lower upper and 1 (for diag)
!  (az equals the sum of entries in upper)

      az = 0
      ia(1) = 1
      do row=2,n
        az = az + ia(row)
        ia(row) = row + lia(row-1) + az
      end do
      ia(n+1) = nnz + 1

      allocate (acsr(nnz), ja(nnz))

!  Copy diagonal elements (in first position)
      do row = 1, n
        ja(ia(row)) = row
        acsr(ia(row)) = diag(row)
      end do

      deallocate (diag)
      nullify(diag)

!  copy entries from upper
      do row = 2, n
        do k = 1, lia(row)-lia(row-1)
          col = lja(k+lia(row-1))
          acsr(ia(row)+k) = upper(k+lia(row-1))
          ja(ia(row)+k) = col
        end do
      end do

      if (.not. symm) then
        deallocate (upper)
        nullify(upper)
      end if

      allocate (lastfree(n-1))

      do row = 1, n-1
        lastfree(row) = ia(row+1) - 1
      end do

!  copy entries from lower
      do col = 2, n
        do k = lia(col-1)+1, lia(col)
          row = lja(k)
          acsr(lastfree(row)) = lower(k)
          ja(lastfree(row)) = col
          lastfree(row) = lastfree(row) - 1 
        end do
      end do

      deallocate (lastfree)
      deallocate (lower)
      nullify(lower)
      deallocate (lia, lja)
      nullify(lia, lja)

      return
      end subroutine lcsr2csc_DP



      subroutine lcsr2csc_DPC(diag, lower, upper, symm, lia, lja, acsr, ia, ja)
      use femtypes
      implicit none
      complex(DPC), pointer :: diag(:), lower(:), upper(:), acsr(:)
      integer(I4B), pointer :: lia(:), lja(:), ja(:), ia(:)
      logical symm
      intent(in) :: symm
      intent(out) :: acsr, ia, ja
      intent(inout) :: diag, lower, upper, lia, lja
!  Convert LCSR sparse matrix format to CSC format

!  Input (LCSR)
!            diag     vector of diagonal matrix enties
!            lower    vector of matrix enties in lower triangle matrix
!            upper    vector of matrix enties in lower triangle matrix
!            sym      = .true. if matrix is symmetric
!            lia      is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            lja      vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix
!
!  Output
!            acsr     CSR matrix elements
!            ia       CSR accumulated sum of the number of non-zero elements in matrix
!            ja       CSR column index of matrix entries

!  local variables
      integer(I4B) :: k, n, lnnz, nnz, az, row, col
      integer(I4B), allocatable :: lastfree(:)

!  Build ia

      n = size(diag)                                        ! number of rows
      lnnz = size(lower)                                    ! number of non-zero entries in lower (or upper) matrix
      nnz = n + 2*lnnz                                      ! number of non-zero entries

      allocate (ia(n+1))
      ia = 0
!  the loop passes array lja and counts matrix elements of value i
!  this will give the number of entries in each row resulting from upper
      do k = 1, lnnz
        ia(lja(k)+1) = ia(lja(k)+1) + 1
      end do

!  build ia from (per row) entiers in lower upper and 1 (for diag)
!  (az equals the sum of entries in upper)

      az = 0
      ia(1) = 1
      do row=2,n
        az = az + ia(row)
        ia(row) = row + lia(row-1) + az
      end do
      ia(n+1) = nnz + 1

      allocate (acsr(nnz), ja(nnz))

!  Copy diagonal elements (in first position)
      do row = 1, n
        ja(ia(row)) = row
        acsr(ia(row)) = diag(row)
      end do

      deallocate (diag)
      nullify(diag)

!  copy entries from upper
      do row = 2, n
        do k = 1, lia(row)-lia(row-1)
          col = lja(k+lia(row-1))
          acsr(ia(row)+k) = upper(k+lia(row-1))
          ja(ia(row)+k) = col
        end do
      end do

      if (.not. symm) then
        deallocate (upper)
        nullify(upper)
      end if

      allocate (lastfree(n-1))

      do row = 1, n-1
        lastfree(row) = ia(row+1) - 1
      end do

!  copy entries from lower
      do col = 2, n
        do k = lia(col-1)+1, lia(col)
          row = lja(k)
          acsr(lastfree(row)) = lower(k)
          ja(lastfree(row)) = col
          lastfree(row) = lastfree(row) - 1 
        end do
      end do

      deallocate (lastfree)
      deallocate (lower)
      nullify(lower)
      deallocate (lia, lja)
      nullify(lia, lja)

      return
      end subroutine lcsr2csc_DPC



!  Convertes from LCSR to LCSR (change of data type)
      subroutine lcsr2lcsr_ZZ(diag_in, lower_in, upper_in, rhs_in, symm, diag_out, lower_out, upper_out, rhs_out)
      use femtypes

      implicit none
      complex(DPC), pointer :: diag_in(:), lower_in(:), upper_in(:), rhs_in(:)
      complex(DPC), pointer :: diag_out(:), lower_out(:), upper_out(:), rhs_out(:)
      logical symm
      intent(in) :: symm
      intent(out) :: diag_out, lower_out, upper_out, rhs_out
      intent(inout) :: diag_in, lower_in, upper_in, rhs_in

      diag_out => diag_in
      rhs_out => rhs_in
      lower_out => lower_in
      if (symm) then
        upper_out => lower_in
      else
        upper_out => upper_in
      end if

      return
      end


      subroutine lcsr2lcsr_ZD(diag_in, lower_in, upper_in, rhs_in, symm, diag_out, lower_out, upper_out, rhs_out)
      use femtypes

      implicit none
      complex(DPC), pointer :: diag_in(:), lower_in(:), upper_in(:), rhs_in(:)
      real(DP), pointer :: diag_out(:), lower_out(:), upper_out(:), rhs_out(:)
      logical symm
      intent(in) :: symm
      intent(out) :: diag_out, lower_out, upper_out, rhs_out
      intent(inout) :: diag_in, lower_in, upper_in, rhs_in
!  local variables
      integer(I4B) :: n, lnnz, ierr
!
      n = size(diag_in)                                        ! number of rows
      lnnz = size(lower_in)                                    ! number of non-zero entries in lower (or upper) matrix

      allocate (diag_out(n))
      diag_out = real(diag_in)
      deallocate (diag_in,stat=ierr)

      allocate (rhs_out(n))
      rhs_out = real(rhs_in)
      deallocate (rhs_in,stat=ierr)

      allocate (lower_out(lnnz))
      lower_out = real(lower_in)
      deallocate (lower_in,stat=ierr)

      if (symm) then
        upper_out => lower_out
      else
        allocate (upper_out(lnnz))
        upper_out = real(upper_in)
        deallocate (upper_in,stat=ierr)
      end if

      return
      end



      subroutine lcsr2lcsr_ZC(diag_in, lower_in, upper_in, rhs_in, symm, diag_out, lower_out, upper_out, rhs_out)
      use femtypes

      implicit none
      complex(DPC), pointer :: diag_in(:), lower_in(:), upper_in(:), rhs_in(:)
      complex(SPC), pointer :: diag_out(:), lower_out(:), upper_out(:), rhs_out(:)
      logical symm
      intent(in) :: symm
      intent(out) :: diag_out, lower_out, upper_out, rhs_out
      intent(inout) :: diag_in, lower_in, upper_in, rhs_in

      integer(I4B) :: n, lnnz

      n = size(diag_in)                                        ! number of rows
      lnnz = size(lower_in)                                    ! number of non-zero entries in lower (or upper) matrix

      allocate (diag_out(n))
      diag_out = real(diag_in)
      deallocate (diag_in)

      allocate (rhs_out(n))
      rhs_out = real(rhs_in)
      deallocate (rhs_in)

      allocate (lower_out(lnnz))
      lower_out = real(lower_in)
      deallocate (lower_in)

      if (symm) then
        upper_out => lower_out
      else
        allocate (upper_out(lnnz))
        upper_out = real(upper_in)
        deallocate (upper_in)
      end if

      return
      end



      subroutine lcsr2lcsr_ZS(diag_in, lower_in, upper_in, rhs_in, symm, diag_out, lower_out, upper_out, rhs_out)
      use femtypes

      implicit none
      complex(DPC), pointer :: diag_in(:), lower_in(:), upper_in(:), rhs_in(:)
      real(SPC), pointer :: diag_out(:), lower_out(:), upper_out(:), rhs_out(:)
      logical symm
      intent(in) :: symm
      intent(out) :: diag_out, lower_out, upper_out, rhs_out
      intent(inout) :: diag_in, lower_in, upper_in, rhs_in

      integer(I4B) :: n, lnnz

      n = size(diag_in)                                        ! number of rows
      lnnz = size(lower_in)                                    ! number of non-zero entries in lower (or upper) matrix

      allocate (diag_out(n))
      diag_out = real(diag_in)
      deallocate (diag_in)

      allocate (rhs_out(n))
      rhs_out = real(rhs_in)
      deallocate (rhs_in)

      allocate (lower_out(lnnz))
      lower_out = real(lower_in)
      deallocate (lower_in)

      if (symm) then
        upper_out => lower_out
      else
        allocate (upper_out(lnnz))
        upper_out = real(upper_in)
        deallocate (upper_in)
      end if

      return
      end


