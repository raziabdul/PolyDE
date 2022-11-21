#undef __FUNCT__
#define __FUNCT__ "petscsolver"

      subroutine petscsolver(a, b, x, n, eps, ia, ja)
      use feminterface, only: itrhstout
      use femtypes
      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscvec.h90"

      complex (DPC), pointer :: a(:)
      complex (DPC) b(:), x(:)
      integer (I4B)  n 
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP), allocatable :: itrvec(:)
      real (DP) eps
      intent (in) :: n, eps
      intent (inout) :: b, x


!    a      matrix (csr format)
!    b       right hand side
!    x      solution vector
!    n       order of matrix
!    eps     accuracy (to do)
!    ia      vector of pointers to the first entry of rows
!    ja      vector of column indices of the entries in a

!=====================================================!
!    $Revision: 1.1 $                                 !
!    $Date: 2004/08/04 23:06:02 $                     !
!    $Author: r_abdul $                               !
!=====================================================!
!   solver with petsc interface 
!-------------------------------------------------------

! local variables
      integer (I4B), allocatable :: nnzrow(:)
      integer (I4B)  i, j, nnzmax, nnzero
!
! Petsc variables types
      KSP              ksp
      PC               pc
      Mat              petA 
      Vec              petx, petb
      PetscScalar,     pointer ::  xff(:)
      double precision norm, rnorm
      PetscViewer      view
      double precision   prvec(1000)
      integer          ierr, rank, flg, size, its, itrlen

!
! read the default ~/.petscrc for options when executing
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
! not useful for single processor, but do anyway
      call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
! 
! get number of columns with nonzeros per row first for
! memory allocation
      allocate(nnzrow(n))
      do i=1,n
         nnzrow(i)=ia(i+1)-ia(i)
      end do
      nnzmax=maxval(nnzrow)
      print *, nnzmax

! allocate memory to Petsc for Mat petA
! PROBLEM: this call doesn't work
      call MatSeqAIJSetPreallocation(petA, 0, nnzrow, ierr)
! Alternative call is also problematic
!      call MatSeqAIJSetPreallocation(petA, nnzmax, PETSC_NULL_INTEGER, ierr)
!      call MatSeqAIJSetPreallocation(petA, 25, PETSC_NULL_INTEGER, ierr)
!      call MatSeqAIJSetPreallocation(petA, nnzmax, 0, ierr)
!      call MatSeqAIJSetPreallocation(petA, nnzmax, nnzrow, ierr)
!      deallocate(nnzrow)

! create petA with own arrays
! ia and ja must suit 0 index in A, so minus 1 for each entry
      call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, n, n, ia-1, &
           ja-1, a, petA, ierr)
! finalize
      call MatAssemblyBegin(petA, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(petA, MAT_FINAL_ASSEMBLY, ierr)
! display pattern
!      call MatView(A, PETSC_VIEWER_DRAW_SELF, ierr)

! print values
!      call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)

! create petb with own arrays
      call VecCreateSeqWithArray(PETSC_COMM_WORLD, n, b, petb, ierr)
!
! finalize
      call VecAssemblyBegin(petb, ierr)
      call VecAssemblyEnd(petb, ierr)

! create solution vector x by duplicating b
      call VecDuplicate(petb, petx, ierr)

! print values
!      call VecView(b, PETSC_VIEWER_STDOUT_WORLD, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!         Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Create linear solver context
      call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)

!  Set operators. Here the matrix that defines the linear system
!  also serves as the preconditioning matrix.
      call KSPSetOperators(ksp, petA, petA, DIFFERENT_NONZERO_PATTERN, ierr)

!  Set runtime options, e.g.,
!      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      call KSPSetFromOptions(ksp,ierr)

! Set Relative tolerance
!   Currently doesn't work, use switch for now
!      call KSPSetTolerances(ksp,1.0e-7,PETSC_DEFAULT_DOUBLE_PRECISION, &
!      &     PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                      Solve the linear system
!
! if not specified, Petsc uses GMRES with ILU(0) as preprocessor for single
! processor codes
! 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Get residuals history
! PROBLEMS:
! the functions KSPSetResidualHistory and KSPGetResidualHistory
! behave really strange:
! You can print ONLY either prvec, or itrlen
! 
! set residual history vector
      call KSPSetResidualHistory(ksp, prvec, itrlen, PETSC_TRUE, ierr)

      call KSPSetRhs(ksp, petb, ierr)
      call KSPSetSolution(ksp, petx, ierr)
      call KSPSolve(ksp, ierr)

! get residual history vector
      call KSPGetResidualHistory(ksp, itrlen, ierr)


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                     Check solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! print values
!      call VecView(petx, PETSC_VIEWER_STDOUT_WORLD, ierr)

!  Get Vec petx , and store in our defined format: complex (DPC)
!  Get a pointer to vector data.
      call VecGetArrayF90(petx, xff, ierr)
          
!  store in actual x to pass to off
       do i=1,n
          x(i)=xff(i)
       end do

! for testing
!       do i=1,n
!          print *, x(i)
!       end do

!  MUST call VecRestoreArray(..) when no longer need access to the array.
       call VecRestoreArrayf90(petx, xff, ierr)

!------------------------------------------------------------------
!  Check the error

! get total iterations
      call KSPGetIterationNumber(ksp, its, ierr)
! get the last (approximate preconditioned) residual norm that has been computed
      call KSPGetResidualNorm(ksp, rnorm, ierr)
!
! print some values
      print *,'==============================================='
      print *, 'iterations :', its
      print *, 'final residual norm :', rnorm
!----------------------------------------------------------------
! BUG? if you print these 2 lines, you itrlen &  prvec become 0
! Need to free the arrays? how ?
!      print *, 'from na : ', itrlen
!      print *, 'from prvec : ', prvec
!---------------------------------------------------------------
      print *,'==============================================='

! nnzero is needed by itrhstout
      nnzero=ia(n+1)-1
! write iterations to file
      call itrhstout(nnzero, its+1, prvec)

!  Free work space AKA deallocate
      call KSPDestroy(ksp, ierr)
      call VecDestroy(petx, ierr)
      call VecDestroy(petb, ierr)
      call MatDestroy(petA, ierr)      
      call PetscFinalize(ierr)
      end subroutine petscsolver

