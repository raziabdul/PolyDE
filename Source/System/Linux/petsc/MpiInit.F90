      subroutine MPI(ierr)
      implicit none
      
! #include "include/petsc/finclude/petsc.h90"      
#include "include/petsc/finclude/petsc.h"
      
      integer :: ierr 
      
      call MPI_init(ierr)
      
      end subroutine MPI
