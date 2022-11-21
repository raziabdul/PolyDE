      subroutine matrhsin(n,n1l,n2l,csr,rhs,ascii)
      use feminterface, only: getsetting
      use femtypes
      implicit none
      integer (I4B) n
      complex (DPC), allocatable :: rhs(:), csr(:)
      integer (I4B), allocatable ::  n1l(:), n2l(:)
      logical ascii
      intent (out) n, n1l, n2l, csr, rhs
      intent (in) ascii
!
!    $Revision: 1.1 $
!    $Date: 2004/07/22 11:32:30 $
!    $Author: r_abdul $
!
!  Reads the system matrix to the file "matrhs"
!
!  Input
!            ascii    =.true. if the to read in ascii format
!            n        number of unknowns (number of dof)
!            csr     vector of non-zero matrix entries
!            rhs      right-hand-side
!            n1l      is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            n2l      vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix   if csr = .false.
!                       - vector of pointers to the first entry of rows if csr = .true.
!
!  local variables
      integer (I4B) unitid, i,jj, ios
      character (len=200) path
      character (len=8) nlab, n1lab, n2lab, csrlab, blab

!   Get the project's path
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
!
      if (ascii) then
!  read as ascii
!  open the file matrhs in the specified path
         open (unitid,file=path(1:len_trim(path))//'matrhs.txt',              &
     &     form='formatted',position='REWIND',action='READ',iostat=ios)
         if (ios .ne. 0) then
            print*,'Error opening file: ',path(1:len_trim(path))//'matrhs.txt'
            print*,'Error no: ',ios
            return
         end if
         read(unitid,*) nlab
         read(unitid,*) n
! allocate n1l:
         allocate(n1l(n+1))
         read(unitid,*) n1lab
         read(unitid,*) n1l(1:n+1)
!  get size of n2l and csr
         jj=n1l(n+1)-1
! allocate the rest:
         allocate(n2l(jj), csr(jj), rhs(n))
!
         read(unitid,*) n2lab
         read(unitid,*) (n2l(i),i=1,jj)
         read(unitid,*) csrlab
         read(unitid,*) (csr(i),i=1,jj)
         read(unitid,*) blab
         read(unitid,*) (rhs(i),i=1,n)
         close (unitid)
      else 
!  read as binary
!  open the file matrhs in the specified path
         open (unitid,file=path(1:len_trim(path))//'matrhs.bin',              &
     &     form='unformatted',position='REWIND',action='READ',iostat=ios)
         if (ios .ne. 0) then
            print*,'Error opening file: ',path(1:len_trim(path))//'matrhs.bin'
            print*,'Error no: ',ios
            return
         end if
         read(unitid) n
! allocate n1l:
         allocate(n1l(n+1))
         read(unitid) n1l(1:n+1)
!  get size of n2l and csr
         jj=n1l(n+1)-1
! allocate the rest:
         allocate(n2l(jj), csr(jj), rhs(n))
!
         read(unitid) (n2l(i),i=1,jj)
         read(unitid) (csr(i),i=1,jj)
         read(unitid) (rhs(i),i=1,n)
         close (unitid)
      end if
      return
      end subroutine matrhsin
