      subroutine matrhsout(n,ia,ja,csr,rhs1,ascii)
      use feminterface, only: getsetting
      use femtypes
      implicit none
      integer (I4B) n
      complex (DPC) rhs1(:), csr(:)
      integer (I4B) ia(:), ja(:)
      logical ascii
      intent (in) n, ia, ja, csr, rhs1, ascii
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
!    $Revision: 1.8 $
!    $Date: 2015/06/10 11:32:52 $
!    $Author: m_kasper $
!

!  Writes the system matrix of CSR format to the file "matrhs"
!
!  Input
!            ascii    =.true. if the to write in ascii format
!            n        number of unknowns (number of dof)
!            acsr     vector of non-zero matrix entries
!            rhs      right-hand-side
!            ia       is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            ja       vector for the compact storage of the matrix
!                       - vector of pointers to the first entry of rows for CSR
!
!  local variables
      integer (I4B) unitid, i,jj, ios
      character (len=200) path
!
      external    ::    grglun

!   Get the project's path
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
!
      if (ascii) then
!  write as ascii
!  open the file matrhs in the specified path
        open (unitid,file=path(1:len_trim(path))//'matrhs.txt',                 &
     &     form='formatted',position='REWIND',action='WRITE',iostat=ios)
        if (ios .ne. 0) then
          print*,'Error opening file: ',path(1:len_trim(path))//'matrhs.txt'
          print*,'Error no: ',ios
          return
        end if
        write(unitid,*) 'NDOF'
        write(unitid,*) n
        write(unitid,*) 'IA'
        write(unitid,*) ia(1:n+1)
        write(unitid,*) 'JA'
!  Get size of ja and csr
        jj=ia(n+1)-1
        write(unitid,*) (ja(i),i=1,jj)
        write(unitid,*) 'CSR'
        write(unitid,*) (csr(i),i=1,jj)
        write(unitid,*) 'RHS'
        write(unitid,*) (rhs1(i),i=1,n)
        close (unitid)
      else 
!        write as binary
!  open the file matrhs in the specified path
        open (unitid,file=path(1:len_trim(path))//'matrhs.txt',                 &
     &     form='unformatted',position='REWIND',action='WRITE',iostat=ios)
        if (ios .ne. 0) then
          print*,'Error opening file: ',path(1:len_trim(path))//'matrhs.txt'
          print*,'Error no: ',ios
          return
        end if
        write(unitid) n
        write(unitid) ia(1:n+1)
!  Get size of ja and csr
        jj=ia(n+1)-1
        write(unitid) (ja(i),i=1,jj)
        write(unitid) (csr(i),i=1,jj)
        write(unitid) (rhs1(i),i=1,n)
        close (unitid)
      end if
!1000  print*,'*** An error occurred while writing the system matrix and the RHS ***'
      return
      end subroutine matrhsout
!
!
!
      subroutine lcsrmatrhsout(n,ia,ja,diag,lower,upper,rhs1,ascii)
      use feminterface, only: getsetting
      use femtypes
      implicit none
      integer (I4B) n
      complex (DPC) rhs1(:), diag(:) ,lower(:) ,upper(:)
      integer (I4B) ia(:), ja(:)
      logical ascii
      intent (in) n, ia, ja, diag, lower, upper, rhs1, ascii
!
!  Writes the system matrix of LCSR format to the file "matrhs"
!
!  Input
!            ascii    =.true. if the to write in ascii format
!            n        number of unknowns (number of dof)
!            diag     vector of diagonal matrix enties
!            lower    vector of matrix enties in lower triangle matrix
!            upper    vector of matrix enties in lower triangle matrix
!            rhs      right-hand-side
!            ia       is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!            ja       vector for the compact storage of the matrix
!                       - column index of the lower triangular matrix for LCSR
!
!  local variables
      integer (I4B) unitid, i,jj, ios
      character (len=200) path
!
      external    ::    grglun

!   Get the project's path
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
!
      if (ascii) then
!  write as ascii
!  open the file matrhs in the specified path
        open (unitid,file=path(1:len_trim(path))//'matrhs.txt',                 &
     &     form='formatted',position='REWIND',action='WRITE',iostat=ios)
        if (ios .ne. 0) then
          print*,'Error opening file: ',path(1:len_trim(path))//'matrhs.txt'
          print*,'Error no: ',ios
          return
        end if
        write(unitid,*) 'NDOF'
        write(unitid,*) n
        write(unitid,*) 'IA'
        write(unitid,*) ia(1:n)
        write(unitid,*) 'JA'
!  Get size of ja and lower/upper
        jj=ia(n)
        write(unitid,*) (ja(i),i=1,jj)
        write(unitid,*) 'DIAG'
        write(unitid,*) (diag(i),i=1,n)
        write(unitid,*) 'LOWER'
        write(unitid,*) (lower(i),i=1,jj)
        write(unitid,*) 'UPPER'
        write(unitid,*) (upper(i),i=1,jj)
        write(unitid,*) 'RHS'
        write(unitid,*) (rhs1(i),i=1,n)
        close (unitid)
      else 
!        write as binary
!  open the file matrhs in the specified path
        open (unitid,file=path(1:len_trim(path))//'matrhs.txt',                 &
     &     form='unformatted',position='REWIND',action='WRITE',iostat=ios)
        if (ios .ne. 0) then
          print*,'Error opening file: ',path(1:len_trim(path))//'matrhs.txt'
          print*,'Error no: ',ios
          return
        end if
        write(unitid) n
        write(unitid) ia(1:n)
!  Get size of ja and lower/upper
        jj=ia(n)
        write(unitid) (ja(i),i=1,jj)
        write(unitid,*) (diag(i),i=1,n)
        write(unitid,*) (lower(i),i=1,jj)
        write(unitid,*) (upper(i),i=1,jj)
        write(unitid) (rhs1(i),i=1,n)
        close (unitid)
      end if
!1000  print*,'*** An error occurred while writing the system matrix and the RHS ***'
      return
      end subroutine lcsrmatrhsout
