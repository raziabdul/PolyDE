      subroutine petscassembly(n1l,jacobi,matvar)
      use feminterface, only: elementmatrix, zeit
      use femtypes
      use globalvariables
      implicit none
      integer (I4B) n1l(:)
      logical jacobi, matvar
      intent (in) ::  n1l, jacobi, matvar

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"

!
!    $Revision: 1.1 $
!    $Date: 2004/11/10 11:44:36 $
!    $Author: r_abdul $
!
!  Element by element FEM matrix and right-hand-side assembly by PETSc
!
!  Input:  
!            n1l      is a vector for the compact storage of the sparse matrix,
!                     it contains the accumulated sum of the number of non-zero 
!                     elements in matrix rows of the full or lower triangular matrix
!                     / = .true. for the computation of the jacobian
!            jacobi
!                     \ = .false. for the linear case
!            n        total number of elements
!            xn, yn   Coordinates of the nodes
!            e        element information, nodes of the elements
!            eg       degrees of freedom (global) for the elements
!                     packed into an array of pointers   eg(i)%d(j)
!            ep       polynomial degree for each of the elements
!            ndof     total number of degrees of freedom
!            csr      n1l contains the number of non-zero entries in each row of the 
!                       -  full matrix  if                  .true. 
!                       -  lower triangular matrix if       .false.
!            matvar   =.true. if the material coefficients are varying across the element
!                     =.false. the material is assumed to be constant 
!            matzif   material numbers of the input regions
!            geb      assignment of elements to the geometry regions
!  Output:  
!
!  In-/ Output:  
!            x        actual solution
!
!  local variables

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                   PETSc Variable declarations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Variables:
!     x, b, u  - approx solution, right-hand-side
!     A        - matrix that defines linear system
!
!  Note that vectors are declared as PETSc "Vec" objects.  These vectors
!  are mathematical objects that contain more than just an array of
!  double precision numbers. I.e., vectors in PETSc are not just
!        double precision x(*).
!  However, local vector data can be easily accessed via VecGetArray().
!  See the Fortran section of the PETSc users manual for details.
!

      PetscInt       i, j, II, JJ, ml, nl, its, Istart, Iend, ione
      PetscErrorCode ierr
      PetscMPIInt     rank,size
      PetscTruth  flg
      PetscScalar v,one,neg_one
      Vec         x, b
      Mat         A


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 PETSC initialization
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      one  = 1.0
      neg_one = -1.0
      ione    = 1
      call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!      Compute the matrix and right-hand-side vector that define
!      the linear system, Ax = b.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Create parallel matrix, specifying only its global dimensions.
!
!      call MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,  &
!           ndof, ndof, A, ierr)
! the above must be followed by this
!      call MatSetFromOptions(A,ierr)

! better alternative, is with:
      ml=2
      nl=2
      call MatCreateMPIAIJ(PETSC_COMM_WORLD, ml, ml, ml*nl, ml*nl, 2, &
           PETSC_NULL_INTEGER, 5, PETSC_NULL_INTEGER, A, ierr)

      call MatGetOwnershipRange(A, Istart, Iend, ierr)

      print *, 'petscassembly OK', rank, size
      print *, ndof, Istart, Iend

      return
      end subroutine petscassembly
