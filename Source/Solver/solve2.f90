      subroutine solve2_c(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use feminterface, only: scal, resid, bsor, umx, fsor, unscal
      use femtypes
      implicit none
      complex (SPC) diag(:), b(:), x(:)
      complex (SPC), pointer :: lower(:), upper(:)
      integer (I4B) n, ia(:), ja(:)
      real (SP) eps, epsgl, resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
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
!    $Date: 2015/10/02 16:31:18 $
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
!                    identical memory location
!
!  local variables
      real (SP) om1, resnrm, xnrm, fehler, scnrm2, bscal
      complex (SPC), allocatable :: gam(:), p(:), r(:), q(:), h(:)
      complex (SPC) lambda, delta0, delta1, beta, alpha, cdotu
      integer (I4B) iter, itmax, itmin, i
      logical converge
      external scnrm2, cdotu
!
!              *   S S O R - C G - M E T H O D    *
!         (Minimum-Residual for non-symmetric matrices)
!
!  prepare
      if (symm) upper=>lower
      converge=.false.
      allocate(gam(n),p(n),r(n),q(n),h(n))
      om1=1._SP
!  minimum and maximum number of iterations
      itmin=max(5,int( (real(n)/10.)**(1./3.) ))
      itmax=int(10*sqrt(real(n))+700./real(n))
!  matrix scaling to obtain diagonal elements to be = 1.
      resnrm=scnrm2(n,b,1)
      if (resnrm.lt.tiny(1._SP)) then
        print*,'** right hand side is zero'
        x(1:n)=(0._SP,0._SP)
        return
      end if
      call scal(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
!  residual computation
      call resid(lower,upper,diag,b,x,ia,ja,n,r)
      fehler=scnrm2(n,r,1)/resnrm
!  preconditioning with  SSOR-matrix
      call bsor(upper,ia,ja,r,r,n,om1)
      call umx(lower,ia,ja,x,x,n,om1)
!
      p=r
      call fsor(lower,ia,ja,r,gam,n,om1)
      q=r-gam
!    q=r-om1*gam      for over-relaxation parameter om1 not equal to  1.
!
      call bsor(upper,ia,ja,q,q,n,om1)
      q=q+gam
!    q=(q+gam)/om1    for over-relaxation parameter om1 not equal to  1.
!

!  iteration loop
      do iter=1,itmax
        lambda=(0.0_SP,0.0_SP)
        delta0=(0.0_SP,0.0_SP)
        do i = 1,n
          lambda = lambda + q(i)*q(i)
          delta0 = delta0 + q(i)*r(i)
        end do
        if (abs(lambda).lt.tiny(1._SP)) then
          itmax=iter
          converge=.false.
          exit
        end if
        alpha=cmplx(cmplx(delta0,kind=DPC)/cmplx(lambda,kind=DPC))
!  new solution
        x(1:n)=x(1:n)+alpha*p
        r=r-alpha*q
!  preconditioning with  SSOR-matrix
        call fsor(lower,ia,ja,r,gam,n,om1)
        h=r-gam
!        h=r-om1*gam    for over-relaxation parameter om1 not equal to  1.
        call bsor(upper,ia,ja,h,h,n,om1)
        delta1=(0.0_SP,0.0_SP)
        do i=1,n 
          h(i)=h(i)+gam(i)
!        h=(h+gam)/om1  for over-relaxation parameter om1 not equal to  1.
          delta1=delta1+h(i)*q(i)
        end do
        beta=-cmplx(cmplx(delta1,kind=DPC)/cmplx(lambda,kind=DPC))
!  test for convergence
        if (abs(alpha).lt.1.e-3_SP .and.                                 &
     &      abs(beta+(1._SP,0._SP)).lt.1.e-3_SP) then
          print*,'*** the linear solver does not converge'
          itmax=iter
          converge=.false.
          exit
        end if
!  L2 norm of solution vector
        if ((mod(iter,itmin).eq.0).or.(iter.eq.1)) then
          xnrm=scnrm2(n,x,1)
        end if
!    Print*,alpha,beta,lambda,delta0,delta1
!  normalized error
        fehler=sqrt(abs(delta1))/xnrm
!  print iteration info
        if (mod(iter,25).eq.0) write(*,111) iter,fehler
111     format (i7,'-th iteration   SSOR error = ',g13.6)
!  test for stopping criteria
        if ( (iter.ge.itmin .and. fehler.lt.eps) .or.                   &
     &    (fehler.lt.1.e3_SP*tiny(1._SP)) )  then
          epsgl=fehler
          itmax=iter
          converge=.true.
          exit
        end if
        q=h+beta*q
        p=r+beta*p
      end do
!
      if (converge) then
        write (*,222) iter,fehler
222     format ('  reached stopping criterion after ',i4,' iterations ',      &
     &  'SSOR error = ',g13.6)
      else
        write (*,333) itmax,fehler
333     format ('  ** did not reach stopping criterion in ',i4,' iterations ',&
     &  'SSOR error = ',g13.6)
        epsgl = fehler
      end if
300   call fsor(lower,ia,ja,x,x,n,om1)
      call unscal(lower,upper,diag,b,x,ia,ja,n,resgl,bscal,h)
      deallocate(gam,p,r,q,h)

      return
      end subroutine solve2_c
!
!
!
      subroutine mp_c(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only:
      use femtypes
      implicit none
      complex (SPC) lower(:), upper(:), diag(:), x(:), b(:)
      integer (I4B) ia(:), ja(:), n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i, j, k
!
!        Matrix  x  Vector  = Vector
!
!          A     x    X    =    B
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
      end subroutine mp_c
!
!
!
      subroutine resid_c(lower,upper,diag,b,x,ia,ja,n,res)
      use feminterface, only:
      use femtypes
      implicit none
      complex (SPC) lower(:), upper(:), diag(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
!  local variables
      integer (I4B) i, j, k
!
!  Comutation of residual:  res = b - A*x
!
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    res     residual vector
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      res(1)=b(1)-x(1)
      do i=2,n
        res(i)=b(i)-x(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          res(i)=res(i)-lower(j)*x(k)
          res(k)=res(k)-upper(j)*x(i)
        end do
      end do
      return
      end subroutine resid_c
!
!
!
      subroutine umx_c(lower,ia,ja,v1,e,n,om1)
      use feminterface, only:
      use femtypes
      implicit none
      complex (SPC) lower(:), v1(:), e(:)
      real (SP) om1
      integer (I4B) ia(:), ja(:), n
      intent (in) :: lower, ia, ja, v1, n, om1
      intent (out) :: e
!  local variables
      integer (I4B) i, j
      complex (SPC) ak
!
!        Matrix  x  Vector  = Vector
!        (I+om1*lower)  x    V1    =    E
!
!        Matrix  = Identity + om1 * lower
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      do i=n,2,-1
        ak=(0._SP,0._SP)
        do j=ia(i-1)+1,ia(i)
          ak=ak+lower(j)*v1(ja(j))
        end do
        e(i)=v1(i)+ak
!    e(i)=v1(i)+ak*om1  for over-relaxation parameter om1 not equal to  1.
      end do
      e(1)=v1(1)
      return
      end subroutine umx_c
!
!
!
      subroutine scal_c(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
      use feminterface, only:
      use femtypes
      implicit none
      complex (SPC) upper(:), lower(:), diag(:), x(:), b(:)
      real (SP) bscal
      logical symm
      integer (I4B) ia(:), ja(:), n
      intent (in) :: ia, ja, n, symm
      intent (out) :: bscal
      intent (inout) :: lower, upper, diag, b, x
!   Scaling of lower and upper triangular matrix, the solution vector X
!   and the right hand side  B with 1./sqrt(diag(i)*diag(j))
!   --> condition of matrix improves by equilibration (equivalent to jacobi preconditioning)
!   --> preserves symmetry of matrix if existent
!   scaling factors are stored in the diagonal of the matrix
!
!  Input
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    symm    =.true. matrix is assumed to be symmetric such that
!                    lower and upper are identical and share 
!                    identical memory location
!  Output 
!    lower   lower triangular matrix (scaled)
!    upper   upper triangular matrix (scaled)
!    diag    diagonal of the matrix  (scaling factors)
!    b       right hand side vector (scaled)
!    x       solution vector (scaled)
!    bscal   normalization factor of right hand side
!
!  local variables
      real (SP) scnrm2
      complex (SPC) sc
      integer (I4B) i, j, stai, endi
      external scnrm2
!
!  compute scaling factor of rows/ columns
      diag(1:n)=1._SP/sqrt(diag(1:n))
!
!  scaling of matrix - rows
      stai=1
      do i=1,n
        sc=diag(i)
        endi=ia(i)
        lower(stai:endi)=lower(stai:endi)*sc
        if (.not. symm) upper(stai:endi)=upper(stai:endi)*sc
        stai=ia(i)+1
      end do
!
!  scaling of matrix - columns
      do j=1,ia(n)
        sc=diag(ja(j))
        lower(j)=lower(j)*sc
        if (.not. symm) upper(j)=upper(j)*sc
      end do
!
!  scaling of right hand side and solution vector
      b(1:n)=b(1:n)*diag(1:n)
!  additional scaling of right hand side such that the norm becomes unity
      bscal=scnrm2(n,b,1)
      b(1:n)=b(1:n)/bscal
      x(1:n)=x(1:n)/(diag(1:n)*bscal)
!
      return
      end subroutine scal_c
!
!
!
      subroutine unscal_c(lower,upper,diag,b,x,ia,ja,n,resnrm,bscal,h)
      use feminterface, only: resid
      use femtypes
      implicit none
      complex (SPC) lower(:), upper(:), diag(:), x(:), b(:), h(:)
      real (SP) resnrm, bscal
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: lower, upper, diag, ia, ja, n, bscal
      intent (out) :: resnrm, h
      intent (inout) ::  b, x
!  back scaling of solution vector and right hand side
!  and computation of residual
!
!  Input
!    lower   lower triangular matrix - not used
!    upper   upper triangular matrix - not used
!    diag    scaling factors
!    b       right hand side vector (scaled)
!    x       solution vector (scaled)
!    ia      compact storage information - not used
!    ja      compact storage information - not used
!    n       number of equations
!    bscal   normalization factor of right hand side
!    h       temporary vector of lenght n
!  Output 
!    b       right hand side vector
!    x       solution vector
!    resnrm  L2 norm of residal
!
!  local variables
      real (SP) scnrm2
      external scnrm2
!
      call resid(lower,upper,diag,b,x,ia,ja,n,h)
      resnrm=scnrm2(n,h,1)
!
      x(1:n)=x(1:n)*diag(1:n)*bscal
      b(1:n)=b(1:n)/diag(1:n)*bscal
      return
      end subroutine unscal_c
!
!
!
      subroutine fsor_c(lower,ia,ja,b,x,n,om1)
      use feminterface, only:
      use femtypes
      implicit none
      complex (SPC) lower(:), x(:), b(:)
      real (SP) om1
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: lower, ia, ja, b, n, om1
      intent (out) :: x
!  Solve triangular system (forward SOR)
!  with  C = ( I + om1 * lower )
!        I      identity matrix
!        lower   lower triangular matrix (compact storage)
!
!  C is the forward iteration matrix of  SSOR-method
!
!  the matrix must be scaled
!
!  local variables
      integer (I4B) i, j1, j2
!
      x(1)=b(1)
      j1=ia(1)+1
      do i=2,n
        j1=ia(i-1)+1
        j2=ia(i)
!  lower triangular matrix only        j < i
        x(i)=b(i)-sum(lower(j1:j2)*x(ja(j1:j2)))
!    x(i)=b(i)-om1*sum(lower(j1:j2),x(ja(j1:j2)))  for over-relaxation parameter om1 not equal to  1.
      end do
      return
      end subroutine fsor_c
!
!
!
      subroutine bsor_c(upper,ia,ja,b,x,n,om1)
      use feminterface, only:
      use femtypes
      implicit none
      complex (SPC) upper(:), b(:), x(:)
      real (SP) om1
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, ia, ja, n, om1
      intent (out) :: x
      intent (inout) :: b
!  Solve triangular system (backward SOR)       (x=b!)
!  with  C = ( I + om1 * upper )
!        I      identity matrix
!        upper  upper triangular matrix (compact storage)
!
!  C is the backward iteration matrix of  SSOR-method
!
!  the matrix must be scaled
!
!  !!!!!!!!! this subroutine overwrites the right hand side b !!!!!!!!!!!!
!  !!!!!!!!!                x is not used                     !!!!!!!!!!!!
!
!  local variables
      complex (SPC) ak
      integer (I4B) i, j, k
!
      do j=n,2,-1
        ak=b(j)
!       ak=om1*b(j)  for over-relaxation parameter om1 not equal to  1.
        do k=ia(j),ia(j-1)+1,-1
!  upper triangular matrix only           j > i=ja(k)
          i=ja(k)
          b(i)=b(i)-ak*upper(k)
        end do
      end do
      return
      end subroutine bsor_c
!
!
!
      subroutine solve2_z(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use feminterface, only: scal, resid, bsor, umx, fsor, unscal
      use femtypes
      implicit none
      complex (DPC) diag(:), b(:), x(:)
      complex (DPC), pointer :: lower(:), upper(:)
      integer (I4B) n, ia(:), ja(:)
      real (DP) eps, epsgl, resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
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
!                    identical memory location
!
!  local variables
      real (DP) om1, resnrm, xnrm, fehler, dznrm2, bscal
      complex (DPC), allocatable :: gam(:), p(:), r(:), q(:), h(:)
      complex (DPC) lambda, delta0, delta1, beta, alpha, zdotu
      integer (I4B) iter, itmax, itmin, i
      logical converge
      external zdotu, dznrm2
!
!              *   S S O R - C G - M E T H O D    *
!         (Minimum-Residual for non-symmetric matrices)
!
!  prepare
      if (symm) upper=>lower
      converge=.false.
      allocate(gam(n),p(n),r(n),q(n),h(n))
      om1=1._DP
!  minimum and maximum number of iterations
      itmin=max(5,int( (real(n)/10.)**(1./3.) ))
      itmax=int(10*sqrt(real(n))+700./real(n))
!  matrix scaling to obtain diagonal elements to be = 1.
      resnrm=dznrm2(n,b,1)
      if (resnrm.lt.tiny(1._DP)) then
        print*,'** right hand side is zero'
        x(1:n)=(0._DP,0._DP)
        return
      end if
      call scal(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
!  residual computation
      call resid(lower,upper,diag,b,x,ia,ja,n,r)
      fehler=dznrm2(n,r,1)/resnrm
!  preconditioning with  SSOR-matrix
      call bsor(upper,ia,ja,r,r,n,om1)
      call umx(lower,ia,ja,x,x,n,om1)
!
      p=r
      call fsor(lower,ia,ja,r,gam,n,om1)
      q=r-gam
!    q=r-om1*gam      for over-relaxation parameter om1 not equal to  1.
!
      call bsor(upper,ia,ja,q,q,n,om1)
      q=q+gam
!    q=(q+gam)/om1    for over-relaxation parameter om1 not equal to  1.
!

!  iteration loop
      do iter=1,itmax
        lambda=(0.0_DP,0.0_DP)
        delta0=(0.0_DP,0.0_DP)
        do i = 1,n
          lambda = lambda + q(i)*q(i)
          delta0 = delta0 + q(i)*r(i)
        end do
        if (abs(lambda).lt.tiny(1._DP)) then
          itmax=iter
          converge=.false.
          exit
        end if
        alpha=delta0/lambda
!  new solution
        x(1:n)=x(1:n)+alpha*p
        r=r-alpha*q
!  preconditioning with  SSOR-matrix
        call fsor(lower,ia,ja,r,gam,n,om1)
        h=r-gam
!        h=r-om1*gam    for over-relaxation parameter om1 not equal to  1.
        call bsor(upper,ia,ja,h,h,n,om1)
        delta1=(0.0_DP,0.0_DP)
        do i=1,n 
          h(i)=h(i)+gam(i)
!        h=(h+gam)/om1  for over-relaxation parameter om1 not equal to  1.
          delta1=delta1+h(i)*q(i)
        end do
        beta=-delta1/lambda
!  test for convergence
        if (abs(alpha).lt.1.e-3_DP.and.                                 &
     &      abs(beta+(1._DP,0._DP)).lt.1.e-3_DP) then
          print*,'*** the linear solver does not converge'
          itmax=iter
          converge=.false.
          exit
        end if
!  L2 norm of solution vector
        if ((mod(iter,itmin).eq.0).or.(iter.eq.1)) then
          xnrm=dznrm2(n,x,1)
        end if
!    Print*,alpha,beta,lambda,delta0,delta1
!  normalized error
        fehler=sqrt(abs(delta1))/xnrm
!  print iteration info
        if (mod(iter,25).eq.0) write(*,111) iter,fehler
111     format (i7,'-th iteration   SSOR error = ',g13.6)
!  test for stopping criteria
        if ( (iter.ge.itmin .and. fehler.lt.eps) .or.                   &
     &    (fehler.lt.1.e3_DP*tiny(1._DP)) )  then
          epsgl=fehler
          itmax=iter
          converge=.true.
          exit
        end if
        q=h+beta*q
        p=r+beta*p
      end do
!
      if (converge) then
        write (*,222) iter,fehler
222     format ('  reached stopping criterion after ',i4,' iterations ',      &
     &  'SSOR error = ',g13.6)
      else
        write (*,333) itmax,fehler
333     format ('  ** did not reach stopping criterion in ',i4,' iterations ',&
     &  'SSOR error = ',g13.6)
        epsgl = fehler
      end if
300   call fsor(lower,ia,ja,x,x,n,om1)
      call unscal(lower,upper,diag,b,x,ia,ja,n,resgl,bscal,h)
      deallocate(gam,p,r,q,h)

      return
      end subroutine solve2_z
!
!
!
      subroutine mp_z(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) lower(:), upper(:), diag(:), x(:), b(:)
      integer (I4B) ia(:), ja(:), n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i, j, k
!
!        Matrix  x  Vector  = Vector
!
!          A     x    X    =    B
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
      end subroutine mp_z
!
!
!
      subroutine resid_z(lower,upper,diag,b,x,ia,ja,n,res)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) lower(:), upper(:), diag(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
!  local variables
      integer (I4B) i, j, k
!
!  Comutation of residual:  res = b - A*x
!
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    res     residual vector
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      res(1)=b(1)-x(1)
      do i=2,n
        res(i)=b(i)-x(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          res(i)=res(i)-lower(j)*x(k)
          res(k)=res(k)-upper(j)*x(i)
        end do
      end do
      return
      end subroutine resid_z
!
!
!
      subroutine umx_z(lower,ia,ja,v1,e,n,om1)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) lower(:), v1(:), e(:)
      real (DP) om1
      integer (I4B) ia(:), ja(:), n
      intent (in) :: lower, ia, ja, v1, n, om1
      intent (out) :: e
!  local variables
      integer (I4B) i, j
      complex (DPC) ak
!
!        Matrix  x  Vector  = Vector
!        (I+om1*lower)  x    V1    =    E
!
!        Matrix  = Identity + om1 * lower
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      do i=n,2,-1
        ak=(0._DP,0._DP)
        do j=ia(i-1)+1,ia(i)
          ak=ak+lower(j)*v1(ja(j))
        end do
        e(i)=v1(i)+ak
!    e(i)=v1(i)+ak*om1  for over-relaxation parameter om1 not equal to  1.
      end do
      e(1)=v1(1)
      return
      end subroutine umx_z
!
!
!
      subroutine scal_z(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) upper(:), lower(:), diag(:), x(:), b(:)
      real (DP) bscal
      integer (I4B) ia(:), ja(:), n
      logical symm
      intent (in) :: ia, ja, n, symm
      intent (out) :: bscal
      intent (inout) :: lower, upper, diag, b, x
!   Scaling of lower and upper triangular matrix, the solution vector X
!   and the right hand side  B with 1./sqrt(diag(i)*diag(j))
!   --> condition of matrix improves by equilibration (equivalent to jacobi preconditioning)
!   --> preserves symmetry of matrix if existent
!   scaling factors are stored in the diagonal of the matrix
!
!  Input
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    symm    =.true. matrix is assumed to be symmetric such that
!                    lower and upper are identical and share 
!                    identical memory location
!  Output 
!    lower   lower triangular matrix (scaled)
!    upper   upper triangular matrix (scaled)
!    diag    diagonal of the matrix  (scaling factors)
!    b       right hand side vector (scaled)
!    x       solution vector (scaled)
!    bscal   normalization factor of right hand side
!
!  local variables
      real (DP) dznrm2
      complex (DPC) sc
      integer (I4B) i, j, stai, endi
      external dznrm2
!
!  compute scaling factor of rows/ columns
      diag(1:n)=1._DP/sqrt(diag(1:n))
!
!  scaling of matrix - rows
      stai=1
      do i=1,n
        sc=diag(i)
        endi=ia(i)
        lower(stai:endi)=lower(stai:endi)*sc
        if (.not. symm) upper(stai:endi)=upper(stai:endi)*sc
        stai=ia(i)+1
      end do
!
!  scaling of matrix - columns
      do j=1,ia(n)
        sc=diag(ja(j))
        lower(j)=lower(j)*sc
        if (.not. symm) upper(j)=upper(j)*sc
      end do
!
!  scaling of right hand side and solution vector
      b(1:n)=b(1:n)*diag(1:n)
!  additional scaling of right hand side such that the norm becomes unity
      bscal=dznrm2(n,b,1)
      b(1:n)=b(1:n)/bscal
      x(1:n)=x(1:n)/(diag(1:n)*bscal)
!
      return
      end subroutine scal_z
!
!
!
      subroutine unscal_z(lower,upper,diag,b,x,ia,ja,n,resnrm,bscal,h)
      use feminterface, only: resid
      use femtypes
      implicit none
      complex (DPC) lower(:), upper(:), diag(:), x(:), b(:), h(:)
      real (DP) resnrm, bscal
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: lower, upper, diag, ia, ja, n, bscal
      intent (out) :: resnrm, h
      intent (inout) ::  b, x
!  back scaling of solution vector and right hand side
!  and computation of residual
!
!  Input
!    lower   lower triangular matrix - not used
!    upper   upper triangular matrix - not used
!    diag    scaling factors
!    b       right hand side vector (scaled)
!    x       solution vector (scaled)
!    ia      compact storage information - not used
!    ja      compact storage information - not used
!    n       number of equations
!    bscal   normalization factor of right hand side
!    h       temporary vector of lenght n
!  Output 
!    b       right hand side vector
!    x       solution vector
!    resnrm  L2 norm of residal
!
!  local variables
      real (DP) dznrm2
      external dznrm2
!
      call resid(lower,upper,diag,b,x,ia,ja,n,h)
      resnrm=dznrm2(n,h,1)
!
      x(1:n)=x(1:n)*diag(1:n)*bscal
      b(1:n)=b(1:n)/diag(1:n)*bscal
      return
      end subroutine unscal_z
!
!
!
      subroutine fsor_z(lower,ia,ja,b,x,n,om1)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) lower(:), x(:), b(:)
      real (DP) om1
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: lower, ia, ja, b, n, om1
      intent (out) :: x
!  Solve triangular system (forward SOR)
!  with  C = ( I + om1 * lower )
!        I      identity matrix
!        lower   lower triangular matrix (compact storage)
!
!  C is the forward iteration matrix of  SSOR-method
!
!  the matrix must be scaled
!
!  local variables
      integer (I4B) i, j1, j2
!
      x(1)=b(1)
      j1=ia(1)+1
      do i=2,n
        j1=ia(i-1)+1
        j2=ia(i)
!  lower triangular matrix only        j < i
        x(i)=b(i)-sum(lower(j1:j2)*x(ja(j1:j2)))
!    x(i)=b(i)-om1*sum(lower(j1:j2),x(ja(j1:j2)))  for over-relaxation parameter om1 not equal to  1.
      end do
      return
      end subroutine fsor_z
!
!
!
      subroutine bsor_z(upper,ia,ja,b,x,n,om1)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) upper(:), b(:), x(:)
      real (DP) om1
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, ia, ja, n, om1
      intent (out) :: x
      intent (inout) :: b
!  Solve triangular system (backward SOR)       (x=b!)
!  with  C = ( I + om1 * upper )
!        I      identity matrix
!        upper  upper triangular matrix (compact storage)
!
!  C is the backward iteration matrix of  SSOR-method
!
!  the matrix must be scaled
!
!  !!!!!!!!! this subroutine overwrites the right hand side b !!!!!!!!!!!!
!  !!!!!!!!!                x is not used                     !!!!!!!!!!!!
!
!  local variables
      complex (DPC) ak
      integer (I4B) i, j, k
!
      do j=n,2,-1
        ak=b(j)
!       ak=om1*b(j)  for over-relaxation parameter om1 not equal to  1.
        do k=ia(j),ia(j-1)+1,-1
!  upper triangular matrix only           j > i=ja(k)
          i=ja(k)
          b(i)=b(i)-ak*upper(k)
        end do
      end do
      return
      end subroutine bsor_z
!
!
!
      subroutine solve2_d(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use feminterface, only: scal, resid, bsor, umx, fsor, unscal
      use femtypes
      implicit none
      real (DP) diag(:), b(:), x(:)
      real (DP), pointer :: lower(:), upper(:)
      integer (I4B) n, ia(:), ja(:)
      real (DP) eps, epsgl, resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
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
!                    identical memory location
!
!  local variables
      real (DP), allocatable :: gam(:), p(:), r(:), q(:), h(:)
      real (DP) om1, resnrm, xnrm, fehler, bscal
      real (DP) lambda, delta0, delta1, beta, alpha
      integer (I4B) iter, itmax, itmin, i
      logical converge
!
!              *   S S O R - C G - M E T H O D    *
!         (Minimum-Residual for non-symmetric matrices)
!
!  prepare
      if (symm) upper=>lower
      converge=.false.
      allocate(gam(n),p(n),r(n),q(n),h(n))
      om1=1._DP
!  minimum and maximum number of iterations
      itmin=max(5,int( (real(n)/10.)**(1./3.) ))
      itmax=int(10*sqrt(real(n))+700./real(n))
!  matrix scaling to obtain diagonal elements to be = 1.
      resnrm=norm2(b(1:n))
      if (resnrm.eq.0.) then
        print*,'** right hand side is zero'
          x(1:n)=0._DP
        return
      end if
      call scal(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
!  residual computation
      call resid(lower,upper,diag,b,x,ia,ja,n,r)
      fehler=norm2(r(1:n))/resnrm
!  preconditioning with  SSOR-matrix
      call bsor(upper,ia,ja,r,r,n,om1)
      call umx(lower,ia,ja,x,x,n,om1)
!
      p=r
      call fsor(lower,ia,ja,r,gam,n,om1)
      q=r-gam
!    q=r-om1*gam      for over-relaxation parameter om1 not equal to  1.
!
      call bsor(upper,ia,ja,q,q,n,om1)
      q=q+gam
!    q=(q+gam)/om1    for over-relaxation parameter om1 not equal to  1.
!

!  iteration loop
      do iter=1,itmax
        lambda=0.0_DP
        delta0=0.0_DP
        do i = 1,n
          lambda = lambda + q(i)*q(i)
          delta0 = delta0 + q(i)*r(i)
        end do
        if (abs(lambda).lt.10._DP*tiny(1._DP)) then
          itmax=iter
          converge=.false.
          exit
        end if
        alpha=delta0/lambda
!  new solution
        x(1:n)=x(1:n)+alpha*p
        r=r-alpha*q
!  preconditioning with  SSOR-matrix
        call fsor(lower,ia,ja,r,gam,n,om1)
        h=r-gam
!        h=r-om1*gam    for over-relaxation parameter om1 not equal to  1.
        call bsor(upper,ia,ja,h,h,n,om1)
        delta1=0.0_DP
        do i=1,n 
          h(i)=h(i)+gam(i)
!        h=(h+gam)/om1  for over-relaxation parameter om1 not equal to  1.
          delta1=delta1+h(i)*q(i)
        end do
        beta=-delta1/lambda
!  test for convergence
        if (abs(alpha).lt.1.e-3_DP.and.abs(beta+1._DP).lt.1.e-3_DP) then
          print*,'*** the linear solver does not converge'
          itmax=iter
          converge=.false.
          exit
        end if
!  L2 norm of solution vector
        if ((mod(iter,itmin).eq.0).or.(iter.eq.1)) then
          xnrm= norm2(x(1:n))
        end if
!    Print*,alpha,beta,lambda,delta0,delta1
!  normalized error
        fehler=sqrt(abs(delta1))/xnrm
!  print iteration info
        if (mod(iter,25).eq.0) write(*,111) iter,fehler
111     format (i7,'-th iteration   SSOR error = ',g13.6)
!  test for stopping criteria
        if ( (iter.ge.itmin .and. fehler.lt.eps) .or.                   &
     &    (fehler.lt.1.e-30_DP) )  then
          epsgl=real(fehler)
          itmax=iter
          converge=.true.
          exit
        end if
        q=h+beta*q
        p=r+beta*p
      end do
!
      if (converge) then
        write (*,222) iter,fehler
222     format ('  reached stopping criterion after ',i4,' iterations ',      &
     &  'SSOR error = ',g13.6)
      else
        write (*,333) itmax,fehler
333     format ('  ** did not reach stopping criterion in ',i4,' iterations ',&
     &  'SSOR error = ',g13.6)
        epsgl = fehler
      end if
!
      call fsor(lower,ia,ja,x,x,n,om1)
      call unscal(lower,upper,diag,b,x,ia,ja,n,resgl,bscal,h)
      deallocate(gam,p,r,q,h)

      return
      end subroutine solve2_d
!
!
!
      subroutine mp_d(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only: 
      use femtypes
      implicit none
      real (DP) lower(:), upper(:), diag(:), x(:), b(:)
      integer (I4B) ia(:), ja(:), n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i, j, k
!
!        Matrix  x  Vector  = Vector
!
!          A     x    X    =    B
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
      end subroutine mp_d
!
!
!
      subroutine resid_d(lower,upper,diag,b,x,ia,ja,n,res)
      use feminterface, only: 
      use femtypes
      implicit none
      real (DP) lower(:), upper(:), diag(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
!  local variables
      integer (I4B) i, j, k
!
!  Comutation of residual:  res = b - A*x
!
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    res     residual vector
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      res(1)=b(1)-x(1)
      do i=2,n
        res(i)=b(i)-x(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          res(i)=res(i)-lower(j)*x(k)
          res(k)=res(k)-upper(j)*x(i)
        end do
      end do
      return
      end subroutine resid_d
!
!
!
      subroutine umx_d(lower,ia,ja,v1,e,n,om1)
      use feminterface, only: 
      use femtypes
      implicit none
      integer (I4B) ia(:), ja(:), n
      real (DP) lower(:), v1(:), e(:), om1
      intent (in) :: lower, ia, ja, v1, n, om1
      intent (out) :: e
!  local variables
      integer (I4B) i, j
      real (DP) ak
!
!        Matrix  x  Vector  = Vector
!        (I+om1*lower)  x    V1    =    E
!
!        Matrix  = Identity + om1 * lower
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      do i=n,2,-1
        ak=0._DP
        do j=ia(i-1)+1,ia(i)
          ak=ak+lower(j)*v1(ja(j))
        end do
        e(i)=v1(i)+ak
!    e(i)=v1(i)+ak*om1  for over-relaxation parameter om1 not equal to  1.
      end do
      e(1)=v1(1)
      return
      end subroutine umx_d
!
!
!
      subroutine scal_d(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
      use feminterface, only: 
      use femtypes
      implicit none
      real (DP) upper(:), lower(:), diag(:), x(:), b(:), bscal
      integer (I4B) ia(:), ja(:), n
      logical symm
      intent (in) :: ia, ja, n, symm
      intent (out) :: bscal
      intent (inout) :: lower, upper, diag, b, x
!   Scaling of lower and upper triangular matrix, the solution vector X
!   and the right hand side  B with 1./sqrt(diag(i)*diag(j))
!   --> condition of matrix improves by equilibration (equivalent to jacobi preconditioning)
!   --> preserves symmetry of matrix if existent
!   scaling factors are stored in the diagonal of the matrix
!
!  Input
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    symm    =.true. matrix is assumed to be symmetric such that
!                    lower and upper are identical and share 
!                    idential memory location
!  Output 
!    lower   lower triangular matrix (scaled)
!    upper   upper triangular matrix (scaled)
!    diag    diagonal of the matrix  (scaling factors)
!    b       right hand side vector (scaled)
!    x       solution vector (scaled)
!    bscal   normalization factor of right hand side
!
!  local variables
      real (DP) sc
      integer (I4B) i, j, stai, endi
!
!  compute scaling factor of rows/ columns
      diag(1:n)=1./sqrt(diag(1:n))
!
!  scaling of matrix - rows
      stai=1
      do i=1,n
        sc=diag(i)
        endi=ia(i)
        lower(stai:endi)=lower(stai:endi)*sc
        if (.not. symm) upper(stai:endi)=upper(stai:endi)*sc
        stai=endi+1
      end do
!
!  scaling of matrix - columns
      do j=1,ia(n)
        sc=diag(ja(j))
        lower(j)=lower(j)*sc
        if (.not. symm) upper(j)=upper(j)*sc
      end do
!
!  scaling of right hand side and solution vector
      b(1:n)=b(1:n)*diag(1:n)
!  additional scaling of right hand side such that the norm becomes unity
      bscal=norm2(b(1:n))
      b(1:n)=b(1:n)/bscal
      x(1:n)=x(1:n)/(diag(1:n)*bscal)
!
      return
      end subroutine scal_d
!
!
!
      subroutine unscal_d(lower,upper,diag,b,x,ia,ja,n,resnrm,bscal,h)
      use feminterface, only: resid
      use femtypes
      implicit none
      real (DP) lower(:), upper(:), diag(:), x(:), b(:), bscal
      real (DP) h(:)
      real (DP) resnrm
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: lower, upper, diag, ia, ja, n, bscal
      intent (out) :: resnrm, h
      intent (inout) ::  b, x
!  back scaling of solution vector and right hand side
!  and computation of residual
!
!  Input
!    lower   lower triangular matrix - not used
!    upper   upper triangular matrix - not used
!    diag    scaling factors
!    b       right hand side vector (scaled)
!    x       solution vector (scaled)
!    ia      compact storage information - not used
!    ja      compact storage information - not used
!    n       number of equations
!    bscal   normalization factor of right hand side
!    h       temporary vector of lenght n
!  Output 
!    b       right hand side vector
!    x       solution vector
!    resnrm  L2 norm of residal
!
!  local variables
!
      call resid(lower,upper,diag,b,x,ia,ja,n,h)
      resnrm=norm2(h(1:n))
!
      x(1:n)=x(1:n)*diag(1:n)*bscal
      b(1:n)=b(1:n)/diag(1:n)*bscal
      return
      end subroutine unscal_d
!
!
!
      subroutine fsor_d(lower,ia,ja,b,x,n,om1)
      use feminterface, only: 
      use femtypes
      implicit none
      real (DP) lower(:), x(:), b(:), om1
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: lower, ia, ja, b, n, om1
      intent (out) :: x
!  Solve triangular system (forward SOR)
!  with  C = ( I + om1 * lower )
!        I      identity matrix
!        lower   lower triangular matrix (compact storage)
!
!  C is the forward iteration matrix of  SSOR-method
!
!  the matrix must be scaled
!
!  local variables
      integer (I4B) i, j1, j2
!
      x(1)=b(1)
      do i=2,n
        j1=ia(i-1)+1
        j2=ia(i)
!  lower triangular matrix only        j < i
        x(i)=b(i)-dot_product(lower(j1:j2),x(ja(j1:j2)))
!    x(i)=b(i)-om1*dot_product(lower(j1:j2),x(ja(j1:j2)))  for over-relaxation parameter om1 not equal to  1.
      end do
      return
      end subroutine fsor_d
!
!
!
      subroutine bsor_d(upper,ia,ja,b,x,n,om1)
      use feminterface, only: 
      use femtypes
      implicit none
      real (DP) upper(:), b(:), x(:), om1
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, ia, ja, n, om1
      intent (out) :: x
      intent (inout) :: b
!  Solve triangular system (backward SOR)       (x=b!)
!  with  C = ( I + om1 * upper )
!        I      identity matrix
!        upper  upper triangular matrix (compact storage)
!
!  C is the backward iteration matrix of  SSOR-method
!
!  the matrix must be scaled
!
!  !!!!!!!!! this subroutine overwrites the right hand side b !!!!!!!!!!!!
!  !!!!!!!!!                x is not used                     !!!!!!!!!!!!
!
!  local variables
      real (DP) ak
      integer (I4B) i, j, k
!
      do j=n,2,-1
        ak=b(j)
!       ak=om1*b(j)  for over-relaxation parameter om1 not equal to  1.
        do k=ia(j),ia(j-1)+1,-1
!  upper triangular matrix only           j > i=ja(k)
          i=ja(k)
          b(i)=b(i)-ak*upper(k)
        end do
      end do
      return
      end subroutine bsor_d
!
!
!
      subroutine solve2_s(lower,upper,diag,b,x,n,eps,ia,ja,epsgl,resgl,symm)
      use feminterface, only: scal, resid, bsor, umx, fsor, unscal
      use femtypes
      implicit none
      real (SP) diag(:), b(:), x(:)
      real (SP), pointer :: lower(:), upper(:)
      integer (I4B) n, ia(:), ja(:)
      real (SP) eps, epsgl, resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: lower, upper, diag, b, x
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
!                    identical memory location
!
!  local variables
      real (SP), allocatable :: gam(:), p(:), r(:), q(:), h(:)
      real (SP) om1, resnrm, xnrm, fehler, bscal
      real (SP) lambda, delta0, delta1, beta, alpha
      integer (I4B) iter, itmax, itmin, i
      logical converge
!
!              *   S S O R - C G - M E T H O D    *
!         (Minimum-Residual for non-symmetric matrices)
!
!  prepare
      if (symm) upper=>lower
      converge=.false.
      allocate(gam(n),p(n),r(n),q(n),h(n))
      om1=1._SP
!  minimum and maximum number of iterations
      itmin=max(5,int( (real(n)/10.)**(1./3.) ))
      itmax=int(10*sqrt(real(n))+700./real(n))
!  matrix scaling to obtain diagonal elements to be = 1.
      resnrm=norm2(b(1:n))
      if (resnrm.eq.0.) then
        print*,'** right hand side is zero'
          x(1:n)=0._SP
        return
      end if
      call scal(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
!  residual computation
      call resid(lower,upper,diag,b,x,ia,ja,n,r)
      fehler=norm2(r(1:n))/resnrm
!  preconditioning with  SSOR-matrix
      call bsor(upper,ia,ja,r,r,n,om1)
      call umx(lower,ia,ja,x,x,n,om1)
!
      p=r
      call fsor(lower,ia,ja,r,gam,n,om1)
      q=r-gam
!    q=r-om1*gam      for over-relaxation parameter om1 not equal to  1.
!
      call bsor(upper,ia,ja,q,q,n,om1)
      q=q+gam
!    q=(q+gam)/om1    for over-relaxation parameter om1 not equal to  1.
!

!  iteration loop
      do iter=1,itmax
        lambda=0.0_SP
        delta0=0.0_SP
        do i = 1,n
          lambda = lambda + q(i)*q(i)
          delta0 = delta0 + q(i)*r(i)
        end do
        if (abs(lambda).lt.10._SP*tiny(1._SP)) then
          itmax=iter
          converge=.false.
          exit
        end if
        alpha=delta0/lambda
!  new solution
        x(1:n)=x(1:n)+alpha*p
        r=r-alpha*q
!  preconditioning with  SSOR-matrix
        call fsor(lower,ia,ja,r,gam,n,om1)
        h=r-gam
!        h=r-om1*gam    for over-relaxation parameter om1 not equal to  1.
        call bsor(upper,ia,ja,h,h,n,om1)
        delta1=0.0_SP
        do i=1,n 
          h(i)=h(i)+gam(i)
!        h=(h+gam)/om1  for over-relaxation parameter om1 not equal to  1.
          delta1=delta1+h(i)*q(i)
        end do
        beta=-delta1/lambda
!  test for convergence
        if (abs(alpha).lt.1.e-3_SP.and.abs(beta+1._SP).lt.1.e-3_SP) then
          print*,'*** the linear solver does not converge'
          itmax=iter
          converge=.false.
          exit
        end if
!  L2 norm of solution vector
        if ((mod(iter,itmin).eq.0).or.(iter.eq.1)) then
          xnrm=norm2(x(1:n))
        end if
!    Print*,alpha,beta,lambda,delta0,delta1
!  normalized error
        fehler=sqrt(abs(delta1))/xnrm
!  print iteration info
        if (mod(iter,25).eq.0) write(*,111) iter,fehler
111     format (i7,'-th iteration   SSOR error = ',g13.6)
!  test for stopping criteria
        if ( (iter.ge.itmin .and. fehler.lt.eps) .or.                   &
     &    (fehler.lt.1.e-30_SP) )  then
          epsgl=real(fehler)
          itmax=iter
          converge=.true.
          exit
        end if
        q=h+beta*q
        p=r+beta*p
      end do
!
      if (converge) then
        write (*,222) iter,fehler
222     format ('  reached stopping criterion after ',i4,' iterations ',      &
     &  'SSOR error = ',g13.6)
      else
        write (*,333) itmax,fehler
333     format ('  ** did not reach stopping criterion in ',i4,' iterations ',&
     &  'SSOR error = ',g13.6)
        epsgl = fehler
      end if
!
      call fsor(lower,ia,ja,x,x,n,om1)
      call unscal(lower,upper,diag,b,x,ia,ja,n,resgl,bscal,h)
      deallocate(gam,p,r,q,h)

      return
      end subroutine solve2_s
!
!
!
      subroutine mp_s(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only: 
      use femtypes
      implicit none
      real (SP) lower(:), upper(:), diag(:), x(:), b(:)
      integer (I4B) ia(:), ja(:), n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i, j, k
!
!        Matrix  x  Vector  = Vector
!
!          A     x    X    =    B
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
      end subroutine mp_s
!
!
!
      subroutine resid_s(lower,upper,diag,b,x,ia,ja,n,res)
      use feminterface, only: 
      use femtypes
      implicit none
      real (SP) lower(:), upper(:), diag(:), x(:), b(:), res(:)
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, lower, diag, b, x, ia, ja, n
      intent (out) :: res
!  local variables
      integer (I4B) i, j, k
!
!  Comutation of residual:  res = b - A*x
!
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    res     residual vector
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      res(1)=b(1)-x(1)
      do i=2,n
        res(i)=b(i)-x(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          res(i)=res(i)-lower(j)*x(k)
          res(k)=res(k)-upper(j)*x(i)
        end do
      end do
      return
      end subroutine resid_s
!
!
!
      subroutine umx_s(lower,ia,ja,v1,e,n,om1)
      use feminterface, only: 
      use femtypes
      implicit none
      integer (I4B) ia(:), ja(:), n
      real (SP) lower(:), v1(:), e(:), om1
      intent (in) :: lower, ia, ja, v1, n, om1
      intent (out) :: e
!  local variables
      integer (I4B) i, j
      real (SP) ak
!
!        Matrix  x  Vector  = Vector
!        (I+om1*lower)  x    V1    =    E
!
!        Matrix  = Identity + om1 * lower
!
!  Here we assume matrix to be scaled, the diagonal elements to be = 1.
      do i=n,2,-1
        ak=0._SP
        do j=ia(i-1)+1,ia(i)
          ak=ak+lower(j)*v1(ja(j))
        end do
        e(i)=v1(i)+ak
!    e(i)=v1(i)+ak*om1  for over-relaxation parameter om1 not equal to  1.
      end do
      e(1)=v1(1)
      return
      end subroutine umx_s
!
!
!
      subroutine scal_s(lower,upper,diag,b,x,ia,ja,n,bscal,symm)
      use feminterface, only: 
      use femtypes
      implicit none
      real (SP) upper(:), lower(:), diag(:), x(:), b(:), bscal
      integer (I4B) ia(:), ja(:), n
      logical symm
      intent (in) :: ia, ja, n, symm
      intent (out) :: bscal
      intent (inout) :: lower, upper, diag, b, x
!   Scaling of lower and upper triangular matrix, the solution vector X
!   and the right hand side  B with 1./sqrt(diag(i)*diag(j))
!   --> condition of matrix improves by equilibration (equivalent to jacobi preconditioning)
!   --> preserves symmetry of matrix if existent
!   scaling factors are stored in the diagonal of the matrix
!
!  Input
!    lower   lower triangular matrix (compact storage)
!    upper   upper triangular matrix (compact storage)
!    diag    diagonal of the matrix
!    b       right hand side vector
!    x       solution vector
!    ia      compact storage information
!    ja      compact storage information
!    n       number of equations
!    symm    =.true. matrix is assumed to be symmetric such that
!                    lower and upper are identical and share 
!                    idential memory location
!  Output 
!    lower   lower triangular matrix (scaled)
!    upper   upper triangular matrix (scaled)
!    diag    diagonal of the matrix  (scaling factors)
!    b       right hand side vector (scaled)
!    x       solution vector (scaled)
!    bscal   normalization factor of right hand side
!
!  local variables
      real (SP) sc
      integer (I4B) i, j, stai, endi
!
!  compute scaling factor of rows/ columns
      diag(1:n)=1./sqrt(diag(1:n))
!
!  scaling of matrix - rows
      stai=1
      do i=1,n
        sc=diag(i)
        endi=ia(i)
        lower(stai:endi)=lower(stai:endi)*sc
        if (.not. symm) upper(stai:endi)=upper(stai:endi)*sc
        stai=endi+1
      end do
!
!  scaling of matrix - columns
      do j=1,ia(n)
        sc=diag(ja(j))
        lower(j)=lower(j)*sc
        if (.not. symm) upper(j)=upper(j)*sc
      end do
!
!  scaling of right hand side and solution vector
      b(1:n)=b(1:n)*diag(1:n)
!  additional scaling of right hand side such that the norm becomes unity
      bscal=norm2(b(1:n))
      b(1:n)=b(1:n)/bscal
      x(1:n)=x(1:n)/(diag(1:n)*bscal)
!
      return
      end subroutine scal_s
!
!
!
      subroutine unscal_s(lower,upper,diag,b,x,ia,ja,n,resnrm,bscal,h)
      use feminterface, only: resid
      use femtypes
      implicit none
      real (SP) lower(:), upper(:), diag(:), x(:), b(:), bscal
      real (SP) h(:)
      real (SP) resnrm
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: lower, upper, diag, ia, ja, n, bscal
      intent (out) :: resnrm, h
      intent (inout) ::  b, x
!  back scaling of solution vector and right hand side
!  and computation of residual
!
!  Input
!    lower   lower triangular matrix - not used
!    upper   upper triangular matrix - not used
!    diag    scaling factors
!    b       right hand side vector (scaled)
!    x       solution vector (scaled)
!    ia      compact storage information - not used
!    ja      compact storage information - not used
!    n       number of equations
!    bscal   normalization factor of right hand side
!    h       temporary vector of lenght n
!  Output 
!    b       right hand side vector
!    x       solution vector
!    resnrm  L2 norm of residal
!
!  local variables
!
      call resid(lower,upper,diag,b,x,ia,ja,n,h)
      resnrm=norm2(h(1:n))
!
      x(1:n)=x(1:n)*diag(1:n)*bscal
      b(1:n)=b(1:n)/diag(1:n)*bscal
      return
      end subroutine unscal_s
!
!
!
      subroutine fsor_s(lower,ia,ja,b,x,n,om1)
      use feminterface, only: 
      use femtypes
      implicit none
      real (SP) lower(:), x(:), b(:), om1
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: lower, ia, ja, b, n, om1
      intent (out) :: x
!  Solve triangular system (forward SOR)
!  with  C = ( I + om1 * lower )
!        I      identity matrix
!        lower   lower triangular matrix (compact storage)
!
!  C is the forward iteration matrix of  SSOR-method
!
!  the matrix must be scaled
!
!  local variables
      integer (I4B) i, j1, j2
!
      x(1)=b(1)
      j1=ia(1)+1
      do i=2,n
        j1=ia(i-1)+1
        j2=ia(i)
!  lower triangular matrix only        j < i
        x(i)=b(i)-dot_product(lower(j1:j2),x(ja(j1:j2)))
!    x(i)=b(i)-om1*dot_product(lower(j1:j2),x(ja(j1:j2)))  for over-relaxation parameter om1 not equal to  1.
      end do
      return
      end subroutine fsor_s
!
!
!
      subroutine bsor_s(upper,ia,ja,b,x,n,om1)
      use feminterface, only: 
      use femtypes
      implicit none
      real (SP) upper(:), b(:), x(:), om1
      integer (I4B) n, ia(:), ja(:)
      intent (in) :: upper, ia, ja, n, om1
      intent (out) :: x
      intent (inout) :: b
!  Solve triangular system (backward SOR)       (x=b!)
!  with  C = ( I + om1 * upper )
!        I      identity matrix
!        upper  upper triangular matrix (compact storage)
!
!  C is the backward iteration matrix of  SSOR-method
!
!  the matrix must be scaled
!
!  !!!!!!!!! this subroutine overwrites the right hand side b !!!!!!!!!!!!
!  !!!!!!!!!                x is not used                     !!!!!!!!!!!!
!
!  local variables
      real (SP) ak
      integer (I4B) i, j, k
!
      do j=n,2,-1
        ak=b(j)
!       ak=om1*b(j)  for over-relaxation parameter om1 not equal to  1.
        do k=ia(j),ia(j-1)+1,-1
!  upper triangular matrix only           j > i=ja(k)
          i=ja(k)
          b(i)=b(i)-ak*upper(k)
        end do
      end do
      return
      end subroutine bsor_s
!
!
!
      subroutine mxv_unscaled_c(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only:
      use femtypes
      implicit none
      complex (SPC) lower(:), upper(:), diag(:), x(:), b(:)
      integer (I4B) ia(:), ja(:), n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i, j, k
!
!        Matrix  x  Vector  = Vector
!
!          A     x    X    =    B
!
!  Here we assume matrix not to be scaled
      b(1)=x(1)*diag(1)
      do i=2,n
        b(i)=x(i)*diag(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          b(i)=b(i)+lower(j)*x(k)
          b(k)=b(k)+upper(j)*x(i)
        end do
      end do
      return
      end subroutine mxv_unscaled_c
!
!
!
      subroutine mxv_unscaled_z(lower,upper,diag,b,x,ia,ja,n)
      use feminterface, only:
      use femtypes
      implicit none
      complex (DPC) lower(:), upper(:), diag(:), x(:), b(:)
      integer (I4B) ia(:), ja(:), n
      intent (in) :: lower, upper, diag, x, ia, ja, n
      intent (out) :: b
!  local variables
      integer (I4B) i, j, k
!
!        Matrix  x  Vector  = Vector
!
!          A     x    X    =    B
!
!  Here we assume matrix not to be scaled
      b(1)=x(1)*diag(1)
      do i=2,n
        b(i)=x(i)*diag(i)
        do j=ia(i-1)+1,ia(i)
!  k = column ,  k < i alway holds
          k=ja(j)
          b(i)=b(i)+lower(j)*x(k)
          b(k)=b(k)+upper(j)*x(i)
        end do
      end do
      return
      end subroutine mxv_unscaled_z
