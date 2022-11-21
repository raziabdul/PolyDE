      subroutine zqmrsolver_DPC(a,b,xzr,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      complex (DPC) a(:),b(:),xzr(:)
      integer (I4B) n,ia(:),ja(:)
      real (DP) eps,epsgl,resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: a, b, xzr
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
!    $Revision: 1.5 $
!    $Date: 2014/03/07 14:05:34 $
!    $Author: m_kasper $
!
!    a       matrix
!    b       right hand side vector
!    xzr     solution vector
!    n       number of equations
!    eps     accuracy to achieve
!    ia      compact storage information
!    ja      compact storage information
!    epsgl   reached error
!    resgl   residual
!    symm    =.true. if a is symmetric
!
!  local variables
!
      integer info(4), m, maxpq, maxvw, mvec, ndim, nlen, nlim
      integer i
      integer, allocatable :: idx(:,:), iwk(:,:)
      complex(DPC), allocatable :: vecs(:,:), ZWK(:,:)
      double precision, allocatable :: dwk(:,:)
      double precision :: norms(2), tol      
!
      EXTERNAL ZAXPBY
      EXTERNAL zscpl, ZSPAXB, ZSPATX, ZSPU2S, ZSSSST, ZSSS1I, ZSS2ST
      EXTERNAL ZUCPL, ZUIL1I, ZUIL1T, ZUIL2I, ZUIL2T, ZUILST
      EXTERNAL ZUSS1I, ZUSS1T, ZUSS2I, ZUSS2T, ZUSSST,ZUS2ST, ZUS2SI

!
!**********************************************************************
!
! Parameters for the data format and the preconditioners.
!
      INTEGER NZMAX, NZLMAX, NZUMAX
!
! Parameters for some of the solvers.
!
!
! Miscellaneous parameters.
!
      COMPLEX(DPC) ZONE, ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0))
!
! Variables for the data format and the preconditioners.
!
      logical, allocatable :: jt(:)
      CHARACTER(len=3) :: TYPE
      integer, allocatable :: ji(:), il(:), jl(:), iu(:), ju(:), ida(:), resvec(:)

      INTEGER NCOL, NROW, PRECON
      complex(DPC), allocatable :: dr(:), l(:), u(:), duinv(:), ainv(:)
      double precision, allocatable :: dn(:), ds(:)
      DOUBLE PRECISION OMEGA
!
!
! Variables common to all solvers.
! All solvers use: NDIM, NLEN, NLIM, VECS, TOL, and INFO.
!
!
! Variables specific to only some of the solvers.
!
!
! Variables used by reverse communication.
!
      INTEGER CB, CX, IERR, REVCOM
!
! Local driver variables.
!
      INTEGER ALG, MAXXPQ, MAXXVW, PRE
      complex(DPC), allocatable :: xtmp(:)
!
      nzmax=size(ja)
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

      nrow=n
      ncol=n
! matris type C/R: compex or real U/S: symmetric or unsymmetric A: assembled
      type='CUA'



      NLEN = NROW
! Choices of algorithm       alg=1 : CPL
!                                2 : CPX
!                                3 : QBG
!                                4 : QMR
!                                5 : QMX
!                                6 : TFX
      alg=1
!
! Get the convergence tolerance.
!
      tol=eps
!
!
! Solver-specific initialization.
!
      IF (ALG.EQ. 1) THEN
        norms(1)=1.d0
        norms(2)=1.d0
!
! Set the maximum block and storage information.
!
        MAXXPQ = MAXPQ
        MAXXVW = MAXVW
        MVEC = MAXPQ + MAXVW
      ELSE IF (ALG.EQ. 4) THEN
        norms(1)=1.d0
!
! Set the maximum block and storage information.
!
        MAXXVW = MAXVW
      END IF
!
! Select a preconditioner
! Choices of preconditioner  precon=0 : no prec'
!                                   1 : Left ILUT'
!                                   2 : Right ILUT'
!                                   3 : Two-sided ILUT'
!                                   4 : Left SSOR'
!                                   5 : Right SSOR'
!                                   6 : Two-sided SSOR'

        if (symm) then
          precon=6
        else
          precon=4
        end if


      IF ((PRECON.EQ.1).OR.(PRECON.EQ.2).OR.(PRECON.EQ.3)) THEN
        PRE = 1
      ELSE IF ((PRECON.EQ.4).OR.(PRECON.EQ.5).OR.(PRECON.EQ.6)) THEN
        PRE = 2
      ELSE
        PRECON = 0
        pre=0
      END IF
!
! Initialize the preconditioner.
!
      IF (PRE.EQ.1) THEN
!  ILUT Preconditioner
!   array allocation:
!  L = the matrix elements (output).
!  IL = the array of row pointers (output).
!  JL = the array of column indices (output).
!  U = the matrix elements for U (output).
!  IU = the array of row pointers for U (output).
!  JU = the array of column indices for U (output).
!  JT = logical work array of size NROW (output).
!  JI = work array of size NROW (output).
!  DN = work array of size NROW (output).
!  DR = work array of size NROW (output).
!  DS = work array of size NROW (output).
!  DUINV = the array of reciprocals of the diagonal of U (output).
!  NZLMAX = the maximum size of L and JL (input).
!  NZUMAX = the maximum size of U and JU (input).
!  INFO = error checking flag, non-zero if not enough memory was
!  allocated for L or U (input).
        allocate (jt(nrow), ji(nrow), dn(nrow), dr(nrow), ds(nrow))
        nzlmax=4*int( (ia(n+1)-n)*1.1 )
        nzumax=nzlmax
        allocate (l(nzlmax), il(nrow+1), jl(nzlmax), u(nzumax), iu(nrow+1), ju(nzumax), duinv(nrow))
        CALL ZUILST (NROW,A,IA,JA,L,IL,JL,U,IU,JU,JT,JI,DN,DR,DS,       &
     &    DUINV,NZLMAX,NZUMAX,IERR)
        deallocate (jt, ji, dn, dr, ds)
      ELSE IF (PRE.EQ.2) THEN
! SSOR preconditioner
! NROW = the number of rows in the matrix (input).
! A = the matrix elements for A (input).
! IA = the array of row pointers for A (input).
! JA = the array of column indices for A (input).
! AINV = the reciprocals of the diagonal elements of A (output).
! IDA = the array of pointers to the diagonal elements (output).
! OMEGA = the SSOR parameter omega
! PRECON = the preconditioner type flag (input/output).
        allocate (ainv(nrow), ida(nrow))
        omega = 1._DP
        if (symm) then
! Convert to modified CSR format.
          CALL ZSPU2S (NROW,A,IA,JA,ainv)
          CALL ZSSSST (NLEN,A,IA,JA,AINV,OMEGA,PRECON)
        else
          CALL ZUSSST (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON)
        end if
      END IF
!
! Compute the modified right hand side.
!
      allocate (xtmp(nrow))
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
      allocate (resvec(nlim))

      do i=1,n
        vecs(i,2)=b(i)
      end do
!
! Compute the modified right hand side.
!
      CALL ZSPAXB (NLEN,A,IA,JA,TYPE,XZR,XTMP)
      CALL ZAXPBY (NLEN,VECS(1,2),ZONE,VECS(1,2),-ZONE,XTMP)
      IF (PRE.EQ.1) THEN
        CALL ZUIL1I (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,VECS(1,2))
      ELSE IF (PRE.EQ.2) THEN
        if (symm) then
          CALL ZSSS1I (NLEN,A,IA,JA,AINV,OMEGA,PRECON,VECS(1,2))
        else
          CALL ZUSS1I (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,VECS(1,2))
        end if
      END IF
!
! Set up call to linear systems solver.
! Compute true residual norms, generate second starting vector.
!
      INFO(2) = 0
      INFO(1) = 000600
!
! Call the solver.
!

 80   IF (ALG.EQ. 1) THEN
        if (symm) then
!  symmetric version
          call zscpl (NDIM,NLEN,NLIM,MAXXPQ,MAXXVW,M,MVEC,NORMS,ZWK,        &
     &      DWK,IDX,IWK,VECS,TOL,INFO, resvec)
!  unsymmetric version
        else
          call zucpl (NDIM,NLEN,NLIM,MAXXPQ,MAXXVW,M,MVEC,NORMS,ZWK,        &
     &      DWK,IDX,IWK,VECS,TOL,INFO, resvec)
        end if
      ELSE IF (ALG.EQ. 2) THEN
!        CALL ZUCPX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      ELSE IF (ALG.EQ. 3) THEN
!        CALL ZUQBG (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      ELSE IF (ALG.EQ. 4) THEN
!        CALL ZUQMR (NDIM,NLEN,NLIM,MAXXVW,M,NORMS(1),ZWK,DWK,IDX,IWK,   &
!     &    VECS,TOL,INFO)
      ELSE IF (ALG.EQ. 5) THEN
!        CALL ZUQMX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      ELSE IF (ALG.EQ. 6) THEN
!        CALL ZUTFX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      END IF
!
! Perform matrix-vector multiplications when needed.
!
      IERR = INFO(1)
      REVCOM = INFO(2)
      CX = INFO(3)
      CB = INFO(4)
!
! Multiply VECS(1,CX) with the preconditioned matrix.
!
      IF (REVCOM.EQ.1) THEN
        IF (PRE.EQ.0) THEN
          CALL ZSPAXB (NROW,A,IA,JA,TYPE,VECS(1,CX),VECS(1,CB))
        ELSE IF (PRE.EQ.1) THEN
          CALL ZAXPBY (NLEN,XTMP,ZONE,VECS(1,CX),ZZERO,XTMP)
          CALL ZUIL2I (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,XTMP)
          CALL ZSPAXB (NROW,A,IA,JA,TYPE,XTMP,VECS(1,CB))
          CALL ZUIL1I (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,VECS(1,CB))
        ELSE IF ((PRECON.EQ.4).OR.(PRECON.EQ.5)) THEN
          CALL ZAXPBY (NLEN,XTMP,ZONE,VECS(1,CX),ZZERO,XTMP)
          CALL ZUSS2I (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,XTMP)
          CALL ZSPAXB (NROW,A,IA,JA,TYPE,XTMP,VECS(1,CB))
          CALL ZUSS1I (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,VECS(1,CB))
        ELSE IF (PRECON.EQ.6) THEN
          CALL ZAXPBY (NLEN,VECS(1,CB),ZONE,VECS(1,CX),ZZERO,XTMP)
          if (symm) then
            CALL ZSS2SI (NROW,A,IA,JA,AINV,OMEGA,PRECON,                  &
     &      VECS(1,CB),XTMP)
          else
            CALL ZUS2SI (NROW,A,IA,JA,AINV,IDA,OMEGA,PRECON,              &
     &      VECS(1,CB),XTMP)
          end if
        END IF
        GO TO 80
!
! Multiply VECS(1,CX) with the preconditioned transpose.
!
      ELSE IF (REVCOM.EQ.2) THEN
        IF (PRE.EQ.0) THEN
          CALL ZSPATX (NLEN,A,IA,JA,TYPE,VECS(1,CX),VECS(1,CB))
        ELSE IF (PRE.EQ.1) THEN
          CALL ZAXPBY (NLEN,XTMP,ZONE,VECS(1,CX),ZZERO,XTMP)
          CALL ZUIL1T (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,XTMP)
          CALL ZSPATX (NLEN,A,IA,JA,TYPE,XTMP,VECS(1,CB))
          CALL ZUIL2T (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,VECS(1,CB))
        ELSE IF ((PRECON.EQ.4).OR.(PRECON.EQ.5)) THEN
          CALL ZAXPBY (NLEN,XTMP,ZONE,VECS(1,CX),ZZERO,XTMP)
          CALL ZUSS1T (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,XTMP)
          CALL ZSPATX (NLEN,A,IA,JA,TYPE,XTMP,VECS(1,CB))
          CALL ZUSS2T (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,VECS(1,CB))
        ELSE IF (PRECON.EQ.6) THEN
          CALL ZAXPBY (NLEN,VECS(1,CB),ZONE,VECS(1,CX),ZZERO,XTMP)
          if (symm) then
            CALL ZSS2ST (NROW,A,IA,JA,AINV,OMEGA,PRECON,                &
     &      VECS(1,CB),XTMP)
          else
          CALL ZUS2ST (NROW,A,IA,JA,AINV,IDA,OMEGA,PRECON,              &
     &      VECS(1,CB),XTMP)
          end if
        END IF
        GO TO 80
      END IF
!
! Check why the solver stopped (this could be more compact).
!
      IF (IERR.EQ. 0) THEN
        WRITE (6,'(A32)') 'The residual norm has converged.'
        GO TO 90
      ELSE IF (IERR.EQ. 1) THEN
        WRITE (6,'(A35)') 'Invalid reverse communication call.'
        GO TO 90
      ELSE IF (IERR.EQ. 2) THEN
        WRITE (6,'(A27)') 'Invalid inputs encountered.'
        GO TO 90
      ELSE IF (IERR.EQ. 4) THEN
        WRITE (6,'(A31)') 'The algorithm did not converge.'
        GO TO 90
      END IF
      IF (ALG.EQ. 1) THEN
        IF (IERR.LT.0) THEN
          WRITE (6,'(A40,I5)')                                          &
     &      'Error encountered in the ZSVDC routine: ',-IERR
          GO TO 90
        ELSE IF (IERR.EQ. 8) THEN
          WRITE (6,'(A35)') 'The last block could not be closed.'
          GO TO 90
        ELSE IF (IERR.EQ. 16) THEN
          WRITE (6,'(A39)') 'An A-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 32) THEN
          WRITE (6,'(A41)') 'An A^T-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 48) THEN
          WRITE (6,'(A41)') 'Both invariant subspaces have been found.'
           GO TO 90
        ELSE IF (IERR.EQ. 64) THEN
          WRITE (6,'(A30)') 'Insufficient memory allocated.'
          GO TO 90
        ELSE IF (IERR.EQ. 128) THEN
          WRITE (6,'(A32)') 'Cannot convert regular to inner.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 2) THEN
        IF (IERR.EQ. 8) THEN
          WRITE (6,'(A25)') 'The algorithm broke down.'
          GO TO 90
        ELSE IF (IERR.EQ. 16) THEN
          WRITE (6,'(A39)') 'An A-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 32) THEN
          WRITE (6,'(A41)') 'An A^T-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 48) THEN
          WRITE (6,'(A41)') 'Both invariant subspaces have been found.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 3) THEN
        IF (IERR.EQ. 8) THEN
          WRITE (6,'(A25)') 'The algorithm broke down.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 4) THEN
        IF (IERR.LT.0) THEN
          WRITE (6,'(A40,I5)')                                          &
     &     'Error encountered in the ZSVDC routine: ',-IERR
          GO TO 90
        ELSE IF (IERR.EQ. 8) THEN
          WRITE (6,'(A35)') 'The last block could not be closed.'
          GO TO 90
        ELSE IF (IERR.EQ. 16) THEN
          WRITE (6,'(A39)') 'An A-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 32) THEN
          WRITE (6,'(A41)') 'An A^T-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 48) THEN
          WRITE (6,'(A41)') 'Both invariant subspaces have been found.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 5) THEN
        IF (IERR.EQ. 8) THEN
          WRITE (6,'(A25)') 'The algorithm broke down.'
          GO TO 90
        ELSE IF (IERR.EQ. 16) THEN
          WRITE (6,'(A39)') 'An A-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 32) THEN
          WRITE (6,'(A41)') 'An A^T-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 48) THEN
          WRITE (6,'(A41)') 'Both invariant subspaces have been found.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 6) THEN
        IF (IERR.EQ. 8) THEN
          WRITE (6,'(A25)') 'The algorithm broke down.'
          GO TO 90
        END IF
      END IF
      WRITE (6,'(A19,I5)') 'Unknown INFO code: ', IERR
!
! Compute the unpreconditioned solution.
!
 90   IF (PRE.EQ.1) THEN
        CALL ZUIL2I (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,VECS(1,1))
      ELSE IF (PRE.EQ.2) THEN
        if (symm) then
          CALL ZSSS2I (NLEN,A,IA,JA,AINV,OMEGA,PRECON,VECS(1,1))
        else
          CALL ZUSS2I (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,VECS(1,1))
        end if
      END IF
      CALL ZAXPBY (NLEN,XZR,ZONE,XZR,ZONE,VECS(1,1))

      deallocate (zwk, dwk, iwk)
      deallocate (idx)
      deallocate (vecs)
      deallocate (xtmp)
      epsgl = tol
      resgl = tol
!
      return
      end subroutine zqmrsolver_DPC



      subroutine zqmrsolver_DP(a,b,xzr,n,eps,ia,ja,epsgl,resgl,symm)
      use femtypes
      implicit none
      real (DP) a(:),b(:),xzr(:)
      integer (I4B) n,ia(:),ja(:)
      real (DP) eps,epsgl,resgl
      logical symm
      intent (in) :: n, eps, ia, ja, symm
      intent (out) :: epsgl, resgl
      intent (inout) :: a, b, xzr
!
!    a       matrix
!    b       right hand side vector
!    xzr     solution vector
!    n       number of equations
!    eps     accuracy to achieve
!    ia      compact storage information
!    ja      compact storage information
!    epsgl   reached error
!    resgl   residual
!    symm    =.true. if a is symmetric
!
!  local variables
!
      integer info(4), m, maxpq, maxvw, mvec, ndim, nlen, nlim
      integer i
      integer, allocatable :: idx(:,:), iwk(:,:)
      double precision, allocatable :: vecs(:,:)
      double precision, allocatable :: dwk(:,:)
      double precision :: norms(2), tol      
!
      EXTERNAL dAXPBY
      EXTERNAL dscpl, dSPAXB, dSPATX, dSPU2S, dSSSST, dSSS1I, dSS2ST
      EXTERNAL dUCPL, dUIL1I, dUIL1T, dUIL2I, dUIL2T, dUILST
      EXTERNAL dUSS1I, dUSS1T, dUSS2I, dUSS2T, dUSSST, dUS2ST, dUS2SI
!
!**********************************************************************
!
! Parameters for the data format and the preconditioners.
!
      INTEGER NZMAX, NZLMAX, NZUMAX
!
! Parameters for some of the solvers.
!
!
! Miscellaneous parameters.
!
      double precision dONE, dZERO
      PARAMETER (dONE = 1.0D0, dZERO = 0.0D0)
!
! Variables for the data format and the preconditioners.
!
      logical, allocatable :: jt(:)
      CHARACTER(len=3) :: TYPE
      integer, allocatable :: ji(:), il(:), jl(:), iu(:), ju(:), ida(:), resvec(:)

      INTEGER NCOL, NROW, PRECON
      double precision, allocatable :: dr(:), l(:), u(:), duinv(:), ainv(:)
      double precision, allocatable :: dn(:), ds(:)
      DOUBLE PRECISION OMEGA
!
!
! Variables common to all solvers.
! All solvers use: NDIM, NLEN, NLIM, VECS, TOL, and INFO.
!
!
! Variables specific to only some of the solvers.
!
!
! Variables used by reverse communication.
!
      INTEGER CB, CX, IERR, REVCOM
!
! Local driver variables.
!
      INTEGER ALG, MAXXPQ, MAXXVW, PRE
      double precision, allocatable :: xtmp(:)
!
      nzmax=size(ja)
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

      nrow=n
      ncol=n
! matris type C/R: compex or real U/S: symmetric or unsymmetric A: assembled
      type='CUA'



      NLEN = NROW
! Choices of algorithm       alg=1 : CPL
!                                2 : CPX
!                                3 : QBG
!                                4 : QMR
!                                5 : QMX
!                                6 : TFX
      alg=1
!
! Get the convergence tolerance.
!
      tol=eps
!
!
! Solver-specific initialization.
!
      IF (ALG.EQ. 1) THEN
        norms(1)=1.d0
        norms(2)=1.d0
!
! Set the maximum block and storage information.
!
        MAXXPQ = MAXPQ
        MAXXVW = MAXVW
        MVEC = MAXPQ + MAXVW
      ELSE IF (ALG.EQ. 4) THEN
        norms(1)=1.d0
!
! Set the maximum block and storage information.
!
        MAXXVW = MAXVW
      END IF
!
! Select a preconditioner
! Choices of preconditioner  precon=0 : no prec'
!                                   1 : Left ILUT'
!                                   2 : Right ILUT'
!                                   3 : Two-sided ILUT'
!                                   4 : Left SSOR'
!                                   5 : Right SSOR'
!                                   6 : Two-sided SSOR'

        if (symm) then
          precon=6
        else
          precon=4
        end if


      IF ((PRECON.EQ.1).OR.(PRECON.EQ.2).OR.(PRECON.EQ.3)) THEN
        PRE = 1
      ELSE IF ((PRECON.EQ.4).OR.(PRECON.EQ.5).OR.(PRECON.EQ.6)) THEN
        PRE = 2
      ELSE
        PRECON = 0
        pre=0
      END IF
!
! Initialize the preconditioner.
!
      IF (PRE.EQ.1) THEN
!  ILUT Preconditioner
!   array allocation:
!  L = the matrix elements (output).
!  IL = the array of row pointers (output).
!  JL = the array of column indices (output).
!  U = the matrix elements for U (output).
!  IU = the array of row pointers for U (output).
!  JU = the array of column indices for U (output).
!  JT = logical work array of size NROW (output).
!  JI = work array of size NROW (output).
!  DN = work array of size NROW (output).
!  DR = work array of size NROW (output).
!  DS = work array of size NROW (output).
!  DUINV = the array of reciprocals of the diagonal of U (output).
!  NZLMAX = the maximum size of L and JL (input).
!  NZUMAX = the maximum size of U and JU (input).
!  INFO = error checking flag, non-zero if not enough memory was
!  allocated for L or U (input).
        allocate (jt(nrow), ji(nrow), dn(nrow), dr(nrow), ds(nrow))
        nzlmax=4*int( (ia(n+1)-n)*1.1 )
        nzumax=nzlmax
        allocate (l(nzlmax), il(nrow+1), jl(nzlmax), u(nzumax), iu(nrow+1), ju(nzumax), duinv(nrow))
        CALL dUILST (NROW,A,IA,JA,L,IL,JL,U,IU,JU,JT,JI,DN,DR,DS,       &
     &    DUINV,NZLMAX,NZUMAX,IERR)
        deallocate (jt, ji, dn, dr, ds)
      ELSE IF (PRE.EQ.2) THEN
! SSOR preconditioner
! NROW = the number of rows in the matrix (input).
! A = the matrix elements for A (input).
! IA = the array of row pointers for A (input).
! JA = the array of column indices for A (input).
! AINV = the reciprocals of the diagonal elements of A (output).
! IDA = the array of pointers to the diagonal elements (output).
! OMEGA = the SSOR parameter omega
! PRECON = the preconditioner type flag (input/output).
        allocate (ainv(nrow), ida(nrow))
        omega = 1._DP
        if (symm) then
! Convert to modified CSR format.
          CALL dSPU2S (NROW,A,IA,JA,ainv)
          CALL dSSSST (NLEN,A,IA,JA,AINV,OMEGA,PRECON)
        else
          CALL dUSSST (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON)
        end if
      END IF
!
! Compute the modified right hand side.
!
      allocate (xtmp(nrow))
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
      allocate (resvec(nlim))

      do i=1,n
        vecs(i,2)=b(i)
      end do
!
! Compute the modified right hand side.
!
      CALL dSPAXB (NLEN,A,IA,JA,TYPE,XZR,XTMP)
      CALL dAXPBY (NLEN,VECS(1,2),dONE,VECS(1,2),-dONE,XTMP)
      IF (PRE.EQ.1) THEN
        CALL dUIL1I (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,VECS(1,2))
      ELSE IF (PRE.EQ.2) THEN
        if (symm) then
          CALL dSSS1I (NLEN,A,IA,JA,AINV,OMEGA,PRECON,VECS(1,2))
        else
          CALL dUSS1I (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,VECS(1,2))
        end if
      END IF
!
! Set up call to linear systems solver.
! Compute true residual norms, generate second starting vector.
!
      INFO(2) = 0
      INFO(1) = 000600
!
! Call the solver.
!

 80   IF (ALG.EQ. 1) THEN
        if (symm) then
!  symmetric version
          call dscpl (NDIM,NLEN,NLIM,MAXXPQ,MAXXVW,M,MVEC,NORMS,        &
     &      DWK,IDX,IWK,VECS,TOL,INFO, resvec)
!  unsymmetric version
        else
          call ducpl (NDIM,NLEN,NLIM,MAXXPQ,MAXXVW,M,MVEC,NORMS,        &
     &      DWK,IDX,IWK,VECS,TOL,INFO, resvec)
        end if
      ELSE IF (ALG.EQ. 2) THEN
!        CALL dUCPX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      ELSE IF (ALG.EQ. 3) THEN
!        CALL dUQBG (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      ELSE IF (ALG.EQ. 4) THEN
!        CALL dUQMR (NDIM,NLEN,NLIM,MAXXVW,M,NORMS(1),DWK,IDX,IWK,   &
!     &    VECS,TOL,INFO)
      ELSE IF (ALG.EQ. 5) THEN
!        CALL dUQMX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      ELSE IF (ALG.EQ. 6) THEN
!        CALL dUTFX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      END IF
!
! Perform matrix-vector multiplications when needed.
!
      IERR = INFO(1)
      REVCOM = INFO(2)
      CX = INFO(3)
      CB = INFO(4)
!
! Multiply VECS(1,CX) with the preconditioned matrix.
!
      IF (REVCOM.EQ.1) THEN
        IF (PRE.EQ.0) THEN
          CALL dSPAXB (NROW,A,IA,JA,TYPE,VECS(1,CX),VECS(1,CB))
        ELSE IF (PRE.EQ.1) THEN
          CALL dAXPBY (NLEN,XTMP,dONE,VECS(1,CX),dZERO,XTMP)
          CALL dUIL2I (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,XTMP)
          CALL dSPAXB (NROW,A,IA,JA,TYPE,XTMP,VECS(1,CB))
          CALL dUIL1I (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,VECS(1,CB))
        ELSE IF ((PRECON.EQ.4).OR.(PRECON.EQ.5)) THEN
          CALL dAXPBY (NLEN,XTMP,dONE,VECS(1,CX),dZERO,XTMP)
          CALL dUSS2I (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,XTMP)
          CALL dSPAXB (NROW,A,IA,JA,TYPE,XTMP,VECS(1,CB))
          CALL dUSS1I (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,VECS(1,CB))
        ELSE IF (PRECON.EQ.6) THEN
          CALL dAXPBY (NLEN,VECS(1,CB),dONE,VECS(1,CX),dZERO,XTMP)
          if (symm) then
            CALL dSS2SI (NROW,A,IA,JA,AINV,OMEGA,PRECON,                  &
     &      VECS(1,CB),XTMP)
          else
            CALL dUS2SI (NROW,A,IA,JA,AINV,IDA,OMEGA,PRECON,              &
     &      VECS(1,CB),XTMP)
          end if
        END IF
        GO TO 80
!
! Multiply VECS(1,CX) with the preconditioned transpose.
!
      ELSE IF (REVCOM.EQ.2) THEN
        IF (PRE.EQ.0) THEN
          CALL dSPATX (NLEN,A,IA,JA,TYPE,VECS(1,CX),VECS(1,CB))
        ELSE IF (PRE.EQ.1) THEN
          CALL dAXPBY (NLEN,XTMP,dONE,VECS(1,CX),dZERO,XTMP)
          CALL dUIL1T (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,XTMP)
          CALL dSPATX (NLEN,A,IA,JA,TYPE,XTMP,VECS(1,CB))
          CALL dUIL2T (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,VECS(1,CB))
        ELSE IF ((PRECON.EQ.4).OR.(PRECON.EQ.5)) THEN
          CALL dAXPBY (NLEN,XTMP,dONE,VECS(1,CX),dZERO,XTMP)
          CALL dUSS1T (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,XTMP)
          CALL dSPATX (NLEN,A,IA,JA,TYPE,XTMP,VECS(1,CB))
          CALL dUSS2T (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,VECS(1,CB))
        ELSE IF (PRECON.EQ.6) THEN
          CALL dAXPBY (NLEN,VECS(1,CB),dONE,VECS(1,CX),dZERO,XTMP)
          if (symm) then
            CALL dSS2ST (NROW,A,IA,JA,AINV,OMEGA,PRECON,                &
     &      VECS(1,CB),XTMP)
          else
          CALL dUS2ST (NROW,A,IA,JA,AINV,IDA,OMEGA,PRECON,              &
     &      VECS(1,CB),XTMP)
          end if
        END IF
        GO TO 80
      END IF
!
! Check why the solver stopped (this could be more compact).
!
      IF (IERR.EQ. 0) THEN
        WRITE (6,'(A32)') 'The residual norm has converged.'
        GO TO 90
      ELSE IF (IERR.EQ. 1) THEN
        WRITE (6,'(A35)') 'Invalid reverse communication call.'
        GO TO 90
      ELSE IF (IERR.EQ. 2) THEN
        WRITE (6,'(A27)') 'Invalid inputs encountered.'
        GO TO 90
      ELSE IF (IERR.EQ. 4) THEN
        WRITE (6,'(A31)') 'The algorithm did not converge.'
        GO TO 90
      END IF
      IF (ALG.EQ. 1) THEN
        IF (IERR.LT.0) THEN
          WRITE (6,'(A40,I5)')                                          &
     &      'Error encountered in the ZSVDC routine: ',-IERR
          GO TO 90
        ELSE IF (IERR.EQ. 8) THEN
          WRITE (6,'(A35)') 'The last block could not be closed.'
          GO TO 90
        ELSE IF (IERR.EQ. 16) THEN
          WRITE (6,'(A39)') 'An A-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 32) THEN
          WRITE (6,'(A41)') 'An A^T-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 48) THEN
          WRITE (6,'(A41)') 'Both invariant subspaces have been found.'
           GO TO 90
        ELSE IF (IERR.EQ. 64) THEN
          WRITE (6,'(A30)') 'Insufficient memory allocated.'
          GO TO 90
        ELSE IF (IERR.EQ. 128) THEN
          WRITE (6,'(A32)') 'Cannot convert regular to inner.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 2) THEN
        IF (IERR.EQ. 8) THEN
          WRITE (6,'(A25)') 'The algorithm broke down.'
          GO TO 90
        ELSE IF (IERR.EQ. 16) THEN
          WRITE (6,'(A39)') 'An A-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 32) THEN
          WRITE (6,'(A41)') 'An A^T-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 48) THEN
          WRITE (6,'(A41)') 'Both invariant subspaces have been found.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 3) THEN
        IF (IERR.EQ. 8) THEN
          WRITE (6,'(A25)') 'The algorithm broke down.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 4) THEN
        IF (IERR.LT.0) THEN
          WRITE (6,'(A40,I5)')                                          &
     &     'Error encountered in the ZSVDC routine: ',-IERR
          GO TO 90
        ELSE IF (IERR.EQ. 8) THEN
          WRITE (6,'(A35)') 'The last block could not be closed.'
          GO TO 90
        ELSE IF (IERR.EQ. 16) THEN
          WRITE (6,'(A39)') 'An A-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 32) THEN
          WRITE (6,'(A41)') 'An A^T-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 48) THEN
          WRITE (6,'(A41)') 'Both invariant subspaces have been found.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 5) THEN
        IF (IERR.EQ. 8) THEN
          WRITE (6,'(A25)') 'The algorithm broke down.'
          GO TO 90
        ELSE IF (IERR.EQ. 16) THEN
          WRITE (6,'(A39)') 'An A-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 32) THEN
          WRITE (6,'(A41)') 'An A^T-invariant subspace has been found.'
          GO TO 90
        ELSE IF (IERR.EQ. 48) THEN
          WRITE (6,'(A41)') 'Both invariant subspaces have been found.'
          GO TO 90
        END IF
      ELSE IF (ALG.EQ. 6) THEN
        IF (IERR.EQ. 8) THEN
          WRITE (6,'(A25)') 'The algorithm broke down.'
          GO TO 90
        END IF
      END IF
      WRITE (6,'(A19,I5)') 'Unknown INFO code: ', IERR
!
! Compute the unpreconditioned solution.
!
 90   IF (PRE.EQ.1) THEN
        CALL dUIL2I (NROW,L,IL,JL,U,IU,JU,DUINV,PRECON,VECS(1,1))
      ELSE IF (PRE.EQ.2) THEN
        if (symm) then
          CALL dSSS2I (NLEN,A,IA,JA,AINV,OMEGA,PRECON,VECS(1,1))
        else
          CALL dUSS2I (NLEN,A,IA,JA,AINV,IDA,OMEGA,PRECON,VECS(1,1))
        end if
      END IF
      CALL dAXPBY (NLEN,XZR,dONE,XZR,dONE,VECS(1,1))

      deallocate (dwk, iwk)
      deallocate (idx)
      deallocate (vecs)
      deallocate (xtmp)
      epsgl = tol
      resgl = tol
!
      return
      end subroutine zqmrsolver_DP


