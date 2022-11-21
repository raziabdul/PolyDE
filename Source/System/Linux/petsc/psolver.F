      program psolver
      use feminterface, only: matrhsin
      use femtypes
      implicit none
!
!    $Revision: 1.3 $
!    $Date: 2004/08/02 18:56:44 $
!    $Author: r_abdul $
!
!  reads binary file matrhs.bin then solves with petsc 

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscvec.h90"


      integer (I4B)  n, i, j, nz, nnz, nzix, unitid
      integer (I4B), pointer :: ia(:), ja(:)
      complex (DPC), pointer :: csr(:), rhs(:)
      logical ascii
!
! Petsc variables types
      KSP              ksp
      PC               pc
      Mat              A 
      Vec              x, b
      double precision norm
      PetscScalar,     pointer :: xff(:)
!      PetscScalar      none
      PetscViewer      view
      integer          ierr, rank, flg, dim, size, its
!
! read the default ~/.petscrc for options when executing
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

! not useful for single processor, but do anyway
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
!
! get values form matrhs.bin file
      ascii=.false.
      call matrhsin(n,ia,ja,csr,rhs,ascii)

! size of csr/ nnzeros:
      nnz=ia(n+1)-1
!
! create A with own arrays
! ia and ja must suit 0 index in A, so minus 1 for each entry
      call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, n, n, ia-1, ja-1, csr, A, ierr)

! finalize
      call MatAssemblyBegin(A , MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(A , MAT_FINAL_ASSEMBLY, ierr)

!      call PetscDrawOpen(PETSC_COMM_SELF, PETSC_NULL_CHARACTER, PETSC_NULL_CHARACTER, &
!& PETSC_DECIDE, PETSC_DECIDE, PETSC_DRAW_FULL_SIZE, PETSC_DRAW_FULL_SIZE, PetscDraw* draw)

! display pattern
!      call MatView(A, PETSC_VIEWER_DRAW_SELF, ierr)

! print values
!      call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)

! create vector with MPI even for 1-processor, since it's simpler
      call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, n, b, ierr)

! fill in vector
! TO DO: may be you can use VecGetArray as down below, seems to work for x
!
      do i=1,n
         call VecSetValue(b, i-1, rhs(i), INSERT_VALUES, ierr);
      end do

! finalize
      call VecAssemblyBegin(b, ierr)
      call VecAssemblyEnd(b, ierr)

! create solution vector x by duplicating b
      call VecDuplicate(b, x, ierr)

! print values
!      call VecView(b, PETSC_VIEWER_STDOUT_WORLD, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!         Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Create linear solver context
      call KSPCreate(PETSC_COMM_WORLD, ksp, ierr)

!  Set operators. Here the matrix that defines the linear system
!  also serves as the preconditioning matrix.
      call KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)

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
      call KSPSetRhs(ksp,b,ierr)
      call KSPSetSolution(ksp,x,ierr)
      call KSPSolve(ksp,ierr)


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                     Check solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! print values
!      call VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)

! save x values to matlab .mat file
!  NOT WORKING: seems not supported by the libraries yet?
!  may need to change compile configuration
!    call PetscViewerMatlabOpen(PETSC_COMM_WORLD,"xp1.mat",PETSC_FILE_CREATE,mviewer,ierr)
!
! save as Matlab ascii instead
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD, 'xpet.m', view, ierr)
      call PetscViewerSetFormat(view,PETSC_VIEWER_ASCII_MATLAB,ierr)
      call VecView(x, view, ierr)
      call PetscViewerDestroy(view,ierr)


!  Get Vec x , and store in our defined format complex (DPC)
       call VecGetArrayF90(x, xff, ierr)
          
!  MUST call VecRestoreArray(..) when no longer need access to the array.
       call VecRestoreArrayf90(x, xff, ierr)


!------------------------------------------------------------------
!  Check the norm of solution (not really usefull)
      call VecNorm(x,NORM_2,norm,ierr)
      call KSPGetIterationNumber(ksp,its,ierr)
      
      print *,'Norm of solution ',norm,'iterations ',its
      
!  Free work space AKA deallocate

      call KSPDestroy(ksp,ierr)
      call VecDestroy(x,ierr)
      call VecDestroy(b,ierr)
      call MatDestroy(A,ierr)      
      call PetscFinalize(ierr)
      stop
      end program psolver

