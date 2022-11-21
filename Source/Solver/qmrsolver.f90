      subroutine qmrsolver_DPC(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use feminterface, only: scal, unscal, zeit, axb, atxb
      use feminterface, only: itrhstout, getsetting, hi2low
      use femtypes
      implicit none
      complex (DPC) lower(:),upper(:),diag(:),b(:),x(:)
      integer (I4B) n,ia(:),ja(:)
      real (DP) eps,epsgl,resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
      external zscpl, zucpl
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
!    $Revision: 1.9 $
!    $Date: 2014/03/07 14:04:13 $
!    $Author: m_kasper $
!
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    n       number of equations
!    eps     accuracy to achieve
!    ia      compact storage information
!    ja      compact storage information
!    epsgl   reached error
!    resgl   residual
!    symm    =.true. matrix is assumed to be symmetric such that
!                    lower and upper are identical and share 
!                    idential memory location
!
!  local variables
!
      integer info(4), m, maxpq, maxvw, mvec, ndim, nlen, nlim
      integer colx, colb, i
      integer, allocatable :: idx(:,:), iwk(:,:)
      real (DP) bscal
      real (DP), allocatable :: resvec(:)
      character (len=3) wrthist
      complex(DPC), allocatable :: vecs(:,:), ZWK(:,:)
      double precision, allocatable :: dwk(:,:)
      double precision :: norms(2), tol      
!
! Variables for zscpl:
!
! ndim : size of vecs, input
! nlen : actual matrix size, input
! nlim : maximum iterations number, input/output
! maxpq: largest size look-ahead block for the PQ sequence, input/output 
! maxvw: largest size look-ahead block for the VW sequence, input/output
! m : size of iwk, dwk, zwk, input
! mvec : maximum vector to build, input/output
! norms : estimates of norm of A. input/output
! zwk :  double complex work array , output
! dwk : double work array, output
! idx : indices of all wrapped variables, output
! iwk: indices for, e.g. sequence of indices of regular vectors, output
! vecs: main vector, input/output
! tol: convergence level, input/output
! info: status flag, input/output
!
      call scal(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
!
!  dimensioned size of array vecs
      ndim=n
!  no of degrees of freedom
      nlen=n
!  max no of iterations
      nlim=int(6*sqrt(real(n))+700./real(n))
      nlim=n
!  size look-ahead block PQ
      maxpq=3
!  size look-ahead block VM
      maxvw=3
!  dimensioned size iwk, dwk zwk
      m=maxpq+maxvw+2
!
      mvec=max(maxpq,maxvw)+2
!
      norms=1.d0
      allocate (zwk(m,8*m+7), iwk(m,13))
      if (symm) then
!  symmetric version
        allocate (dwk(m,7))
        allocate (idx(4,nlim+2))
        allocate (vecs(ndim,3*mvec+2))
      else
!  unsymmetric version
        allocate (dwk(m,11))
        allocate (idx(6,nlim+2))
        allocate (vecs(ndim,5*mvec+3))
      end if
      tol=eps
      info(1)=000600
      info(2)=0
!
      do i=1,n
        vecs(i,2)=b(i)
      end do
!
! initialize residual vector
      allocate (resvec(nlim))
      do

        if (symm) then
!  symmetric version
          call zscpl (NDIM,NLEN,NLIM,MAXPQ,MAXVW,M,MVEC,NORMS,ZWK,        &
     &      DWK,IDX,IWK,VECS,TOL,INFO, resvec)
!  unsymmetric version
        else
          call zucpl (NDIM,NLEN,NLIM,MAXPQ,MAXVW,M,MVEC,NORMS,ZWK,        &
     &      DWK,IDX,IWK,VECS,TOL,INFO, resvec)
        end if
!  reverse communication
        colx   = info(3)
        colb   = info(4)
        select case(info(2))
        case(1)
!  multiply by matrix
          call axb(lower,upper,diag,vecs(1:n,colb),vecs(1:n,colx),ia,ja,n)
        case(2)
!  multiply by transposed matrix
          call atxb(lower,upper,diag,vecs(1:n,colb),vecs(1:n,colx),ia,ja,n)
        case default
          exit
        end select

      end do

      if (info(1) .ne. 0) then
! 0   normal termination, algorithm converged,
! 1   invalid reverse-communication call,
! 2   invalid inputs, at least one of the input parameters has a value out of its valid range,
! 4   algorithm reached the limit of look-ahead Lanczos step without converging,
! 8   the last block in one of the sequences cannot be closed, nor is it possible to force closure,
! 16  an A-invariant subspace has been found,
! 32  an A^T-invariant subspace has been found,
! 48  both invariant subspaces have been found,
! 64  insufficient memory allocated in VECS,
! 128 could not rebuild a regular vector.
        print*,'ierr :',info(1)
      else
!  check convergence  ?????????????????????
        print*,'norms :',norms(1),norms(2)
      end if
!  write residual to file if requested
      call getsetting('WRITE ITERATION',wrthist)
      call hi2low(wrthist,3)
      if (wrthist .eq. 'yes') then
         call itrhstout(2*ia(n)+n, nlim, resvec, n)
      end if

      deallocate (zwk, dwk, iwk)
      deallocate (idx,resvec)
!
      do i=1,n
        x(i)=vecs(i,1)
      end do
      call unscal(lower,upper,diag,b,x,ia,ja,n,resgl,bscal,vecs(1:n,2))
      epsgl = resgl
      deallocate (vecs)
!
      return
      end subroutine qmrsolver_DPC



      subroutine qmrsolver_DP(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use feminterface, only: scal, unscal, zeit, axb, atxb
      use feminterface, only: itrhstout, getsetting, hi2low
      use femtypes
      implicit none
      real (DP) lower(:),upper(:),diag(:),b(:),x(:)
      integer (I4B) n,ia(:),ja(:)
      real (DP) eps,epsgl,resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
      external dscpl, ducpl
!
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    n       number of equations
!    eps     accuracy to achieve
!    ia      compact storage information
!    ja      compact storage information
!    epsgl   reached error
!    resgl   residual
!    symm    =.true. matrix is assumed to be symmetric such that
!                    lower and upper are identical and share 
!                    idential memory location
!
!  local variables
!
      integer info(4), m, maxpq, maxvw, mvec, ndim, nlen, nlim
      integer colx, colb, i
      integer, allocatable :: idx(:,:), iwk(:,:)
      real (DP) bscal
      real (DP), allocatable :: resvec(:)
      character (len=3) wrthist
      real(DP), allocatable :: vecs(:,:)
      double precision, allocatable :: dwk(:,:)
      double precision :: norms(2), tol      
!
! Variables for zscpl:
!
! ndim : size of vecs, input
! nlen : actual matrix size, input
! nlim : maximum iterations number, input/output
! maxpq: largest size look-ahead block for the PQ sequence, input/output 
! maxvw: largest size look-ahead block for the VW sequence, input/output
! m : size of iwk, dwk, zwk, input
! mvec : maximum vector to build, input/output
! norms : estimates of norm of A. input/output
! zwk :  double complex work array , output
! dwk : double work array, output
! idx : indices of all wrapped variables, output
! iwk: indices for, e.g. sequence of indices of regular vectors, output
! vecs: main vector, input/output
! tol: convergence level, input/output
! info: status flag, input/output
!
      call scal(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
!
!  dimensioned size of array vecs
      ndim=n
!  no of degrees of freedom
      nlen=n
!  max no of iterations
      nlim=int(6*sqrt(real(n))+700./real(n))
      nlim=n
!  size look-ahead block PQ
      maxpq=3
!  size look-ahead block VM
      maxvw=3
!  dimensioned size iwk, dwk zwk
      m=maxpq+maxvw+2
!
      mvec=max(maxpq,maxvw)+2
!
      norms=1.d0

      allocate (iwk(m,13))
      if (symm) then
!  symmetric version
        allocate (dwk(m,8*m+14))
        allocate (idx(4,nlim+2))
        allocate (vecs(ndim,3*mvec+2))
      else
!  unsymmetric version
        allocate (dwk(m,8*m+18))
        allocate (idx(6,nlim+2))
        allocate (vecs(ndim,5*mvec+3))
      end if
      tol=eps
      info(1)=000600
      info(2)=0
!
      do i=1,n
        vecs(i,2)=b(i)
      end do
!
! initialize residual vector
      allocate (resvec(nlim))
      do

        if (symm) then
!  symmetric version
          call dscpl (NDIM,NLEN,NLIM,MAXPQ,MAXVW,M,MVEC,NORMS,         &
     &      DWK,IDX,IWK,VECS,TOL,INFO, resvec)
        else
!  unsymmetric version
          call ducpl (NDIM,NLEN,NLIM,MAXPQ,MAXVW,M,MVEC,NORMS,         &
     &      DWK,IDX,IWK,VECS,TOL,INFO, resvec)
        end if
!  reverse communication
        colx   = info(3)
        colb   = info(4)
        select case(info(2))
        case(1)
!  multiply by matrix
          call axb(lower,upper,diag,vecs(1:n,colb),vecs(1:n,colx),ia,ja,n)
        case(2)
!  multiply by transposed matrix
          call atxb(lower,upper,diag,vecs(1:n,colb),vecs(1:n,colx),ia,ja,n)
        case default
          exit
        end select

      end do

      if (info(1) .ne. 0) then
! 0   normal termination, algorithm converged,
! 1   invalid reverse-communication call,
! 2   invalid inputs, at least one of the input parameters has a value out of its valid range,
! 4   algorithm reached the limit of look-ahead Lanczos step without converging,
! 8   the last block in one of the sequences cannot be closed, nor is it possible to force closure,
! 16  an A-invariant subspace has been found,
! 32  an A^T-invariant subspace has been found,
! 48  both invariant subspaces have been found,
! 64  insufficient memory allocated in VECS,
! 128 could not rebuild a regular vector.
        print*,'ierr :',info(1)
      else
!  check convergence  ?????????????????????
        print*,'norms :',norms(1),norms(2)
      end if
!  write residual to file if requested
      call getsetting('WRITE ITERATION',wrthist)
      call hi2low(wrthist,3)
      if (wrthist .eq. 'yes') then
         call itrhstout(2*ia(n)+n, nlim, resvec, n)
      end if

      deallocate (dwk, iwk)
      deallocate (idx,resvec)
!
      do i=1,n
        x(i)=vecs(i,1)
      end do
      call unscal(lower,upper,diag,b,x,ia,ja,n,resgl,bscal,vecs(1:n,2))
      epsgl = resgl
      deallocate (vecs)
!
      return
      end subroutine qmrsolver_DP



      subroutine axb_DPC(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) lower(:),upper(:),diag(:)
      complex (DPC) x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i,j,k
!
!        Matrix  x  Vector  = Vector
!
!          A     *    X    =    B
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      b(1)=x(1)
      do i=2,n
        b(i)=x(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          b(i)=b(i)+lower(j)*x(k)
          b(k)=b(k)+upper(j)*x(i)
        end do
      end do
      return
      end subroutine axb_DPC



      subroutine atxb_DPC(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) lower(:),upper(:),diag(:)
      complex (DPC) x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i,j,k
!
!        transposed Matrix  x  Vector  = Vector
!
!          A^T     *    X    =    B
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      b(1)=x(1)
      do i=2,n
        b(i)=x(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          b(i)=b(i)+upper(j)*x(k)
          b(k)=b(k)+lower(j)*x(i)
        end do
      end do
      return
      end subroutine atxb_DPC



      subroutine axb_DP(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) lower(:),upper(:),diag(:)
      real (DP) x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i,j,k
!
!        Matrix  x  Vector  = Vector
!
!          A     *    X    =    B
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      b(1)=x(1)
      do i=2,n
        b(i)=x(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          b(i)=b(i)+lower(j)*x(k)
          b(k)=b(k)+upper(j)*x(i)
        end do
      end do
      return
      end subroutine axb_DP



      subroutine atxb_DP(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) lower(:),upper(:),diag(:)
      real (DP) x(:),b(:)
      integer (I4B) ia(:),ja(:),n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i,j,k
!
!        transposed Matrix  x  Vector  = Vector
!
!          A^T     *    X    =    B
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      b(1)=x(1)
      do i=2,n
        b(i)=x(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          b(i)=b(i)+upper(j)*x(k)
          b(k)=b(k)+lower(j)*x(i)
        end do
      end do
      return
      end subroutine atxb_DP
!
