      
      subroutine petscsolver(a, b, x, n, eps, ia, ja)
#include <petsc/finclude/petsc.h>
      use feminterface, only: getsetting,itrhstout,zeit
      use femtypes
      use petsc
      implicit none
            

!----------------------------------------------------------------
!                   Macro definitions
!----------------------------------------------------------------
!
!  Macros to make clearer the process of setting values in vectors and
!  getting values from vectors.
!
!   - The element xx_a(ib) is element ib+1 in the vector x
!   - Here we add 1 to the base array index to facilitate the use of
!     conventional Fortran 1-based array indexing.
!
!   This may be obsolete#define xx_a(ib)  xx_v(xx_i + (ib))
!  #define xx_a(ib)
      complex (DP), pointer :: a(:)
      complex (DP) b(:), x(:)
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP) eps
      integer (I4B)  n
      intent (in) :: n, eps
      intent (inout) :: b, x

!    a      matrix (csr format)
!    b       right hand side
!    n       order of matrix
!    x      solution vector
!    eps     accuracy (to do)
!    ia      vector of pointers to the first entry of rows
!    ja      vector of column indices of the entries in a

!========================================================
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
!    $Date: 2006/07/08 02:28:45 $
!    $Author: r_abdul $
!=======================================================
!   solver with petsc interface 
!-------------------------------------------------------

! local variables
      integer (I4B)  i, j, nz, nnz, nzix
      integer (I4B), allocatable :: nnzrow(:)


      !
! Petsc variables types
! IMPORTANT:
!     xx_v   a pointer to retrieve the Vec type solution vector
!               and then store in x
!               Its type is determined by which library it is compiled
!               with, eg, g gives real, g_complex gives complex
!
!=======================================================
      KSP              ksp
      PC               pc
      Mat              petA 
      Vec              petx, petb, D, u
      double precision norm, rnorm
      PetscOffset      xx_i
      PetscViewer      view
      PetscReal        zero
      PetscBool       issymm, isH, iszero, flg
      PetscScalar      sum, atmp
      PetscScalar, pointer :: xx_v(:)
      double precision   prvec(50000)
      PetscInt         its, itrlen
      character (len=4) writeItr
      PCType           pcname
      PCSide           side
      PetscRandom      rctx   
      PetscMPIInt      rank, size  
      PetscErrorCode  ierr 
      
      external mymonitor
      
      print *, 'INIT PETSC'

      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      print *, 'DONE INIT'
!      call MPI_Comm_size(PETSC_COMM_SELF, size, ierr)
!      call MPI_Comm_rank(PETSC_COMM_SELF, rank, ierr)
      print *, 'DONE MPI'
      
    
      allocate(nnzrow(n))
      print *, 'DONE ALLOCATE'
! get number of columns per row first
      do i=1,n
         nnzrow(i)=ia(i+1)-ia(i)
      enddo

! size of csr/ nnzeros:
      nnz=ia(n+1)-1
      nz=maxval(nnzrow)
!
! set A to be symmetric, supposed to be faster ??
      call MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,nz,nnzrow,petA,ierr)
print *, 'DONE MatCreateSeqAIJ'
      call MatSetOption(petA, MAT_SYMMETRIC, PETSC_TRUE, ierr)
!
      call zeit('start of PETSc')
      nzix=1
      do i=1,n
       do j=ia(i),ia(i+1)-1
        atmp = a(nzix)
        call MatSetValue(petA,i-1,ja(nzix)-1,atmp,INSERT_VALUES,ierr)
        nzix=nzix+1
       end do
      end do
      call zeit('Building PETSc Array')

      call MatAssemblyBegin(petA, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(petA, MAT_FINAL_ASSEMBLY, ierr)

! check symmetry :
      zero = 1e-16
      call MatIsSymmetric(petA, zero, issymm, ierr)
      if (issymm == 1 ) then
        print *, 'PETSc: Matrix A is symmetric'
      else
         print *, 'PETSc WARNING: Matrix A is NOT symmetric!'
      end if

! create vector and set the type
      call VecCreate(PETSC_COMM_WORLD, petb, ierr)
      call VecSetSizes(petb ,PETSC_DECIDE, n, ierr)
      call VecSetType(petb, VECSEQ, ierr)

! fill in vector
! TO DO: may be you can use VecGetArray as down below, seems to work for x
!
print*, 'Start PETSc RHS'
      do i=1,n
         call VecSetValue(petb, i-1, b(i), INSERT_VALUES, ierr);
      end do
      call zeit('Building PETSc RHS')

! finalize
print*, 'FINALIZE PETB'
      call VecAssemblyBegin(petb, ierr)
      call VecAssemblyEnd(petb, ierr)

      print*, 'CREATE PETX'
! create solution vector x by duplicating b
      call VecDuplicate(petb, petx, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!         Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Create linear solver context
!  _world and _SELF seem no difference
call KSPCreate(PETSC_COMM_SELF, ksp, ierr)
print*, 'KSPCreate DONE>>>>>>>>>>>>>>>>>>>>'

!  Set operators. Here the matrix that defines the linear system
!  also serves as the preconditioning matrix.
      call KSPSetOperators(ksp,petA,petA,ierr)
print*, 'KSPSetOperators DONE>>>>>>>>>>>>>>>>>>>>'

!  If we want to fix the solver type:
      call KSPSetType(ksp,KSPCG,ierr)
!  This MUST come after KSPSetFromOptions
      call KSPCGSetType(ksp,KSP_CG_SYMMETRIC, ierr )
!  Set runtime options, e.g.,
!      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      call KSPSetFromOptions(ksp,ierr)
print*, 'KSPSetOption DONE>>>>>>>>>>>>>>>>>>>>'


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                      Solve the linear system
!
! if not specified, Petsc uses GMRES with ILU(0) as preprocessor for single
! processor codes
! 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! set residual history vector
      call KSPSetResidualHistory(ksp, prvec, 50000, PETSC_TRUE, ierr)
      print *, 'Start KSPSolver ...'
      call KSPSolve(ksp, petb, petx, ierr)
      call zeit('PETSc matrix solver')

! check to make sure about the initial guess
      call KSPGetInitialGuessNonzero(ksp,iszero,ierr)
      if (iszero == 1 ) then
        print *, 'PETSc: x0 is not zero', iszero
      else 
        print *, 'PETSc : x0 is zero', iszero
     end if

! check to see which side is the preconditioner applied
!      call KSPGetPreconditionerSide(ksp,side,ierr)
!      print *, 'PC SIDE'
!      call PetscObjectView(side,PETSC_VIEWER_STDOUT_WORLD, ierr)

! get residual history vector
      call KSPGetResidualHistory(ksp, itrlen, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                     Check solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! print values
!      call VecView(petx, PETSC_VIEWER_STDOUT_WORLD, ierr)

!  Get Vec petx , and store in our defined format: real (DP)
!  Get a pointer to vector data.
!    - For default PETSc vectors, VecGetArray() returns a pointer to
!      the data array.  Otherwise, the routine is implementation dependent.
      print *, 'READ Petsc X ...'
      call VecGetArrayF90(petx,xx_v,ierr)
      x = xx_v

!  MUST call VecRestoreArrayF90(..) when no longer need access to the array.
       call VecRestoreArrayF90(petx,xx_v,ierr)

!------------------------------------------------------------------
!  Check the error
      call KSPGetIterationNumber(ksp, its, ierr)
! get the last (approximate preconditioned) residual norm that has been computed
      call KSPGetResidualNorm(ksp, rnorm, ierr)

! print some values
      print *,'==============================================='
      print *, 'iterations          :', its
      print *, 'final residual norm :', rnorm
      print *, 'itrlen value        :', itrlen
      print *, 'last entry in prvec :', prvec(itrlen)
      print *,'==============================================='


! write iterations to file; nnzero is needed by itrhstout
      call getsetting('WRITE ITERATION', writeItr)
      if (writeItr == 'YES') call itrhstout(nnz,its+1,prvec(1:its+1),n)

!  Free work space AKA deallocate
      deallocate(nnzrow)
      call KSPDestroy(ksp, ierr)
      call VecDestroy(petx, ierr)
      call VecDestroy(petb, ierr)
      call MatDestroy(petA, ierr)


      end subroutine petscsolver


 
