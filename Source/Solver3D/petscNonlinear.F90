      subroutine petscNonlinear(xn)
#include <petsc/finclude/petsc.h>
      use feminterface,        only: lcsr2lcsr, lcsrzeroremover, lcsr2csr, itrhstout
      use feminterface3D,      only:  getbcval, solve_linear3D, assembly3D
      use globalvariables3D,   only: numdof, nnat, numv, vgdof, sbc, bctype, pvalue, nod
      use globalvariables3D,   only: x, vn, vv 
      use femtypes
      use petsc
      implicit none
      complex (DPC),      pointer :: xn(:)
      intent (inout)              :: xn

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
!   Interface to PETSc nonlinear solvers 
!----------------------------------------------------
! local variables:
!    a       matrix (csr format)
!    b       right hand side
!    n       order of matrix
!    xn      solution vector
!    eps     accuracy (to do)
!    ia      vector of pointers to the first entry of rows
!    ja      vector of column indices of the entries in a

      real (DP)              :: epsgl, resgl
      real (DP), pointer     :: lower_DP(:)=>null(), upper_DP(:)=>null(), diag_DP(:)=>null(), rhs_DP(:)=>null()
      real (DP), pointer     :: a_aux_DP(:)=>null(), xn2(:)
      complex (DPC), pointer :: res_iter(:)
      complex (DPC), pointer :: a(:)
      complex (DPC), pointer :: lower(:)=>null(), upper(:)=>null(), diag(:)=>null(), rhs(:)=>null()
      complex (DPC), pointer :: acsr(:)=>null()
      complex (DPC)          :: pval
      integer (I4B),    parameter :: f2n(3,4)=reshape((/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/))
      integer (I4B), pointer :: ia(:), ja(:)
      integer (I4B), pointer :: ia_aux(:), ja_aux(:)
      integer (I4B), allocatable :: nnzrow(:)
      integer (I4B)              :: i, j, nz, nnz, nzix, f, elem, surf, bc_index, inat, nd
      integer (I4B)              :: iter, iterMax = 3
      real (DP)              :: xyzs(3)
      logical                    :: jacobi, csr, matvar, symmetric
!
! local Petsc variables:
! IMPORTANT:
!     xx_v   a pointer to retrieve the Vec type solution vector
!               and then store in x
!               Its type is determined by which library it is compiled
!               with, eg, g gives real, g_complex gives complex
!
!=======================================================
      SNES             snes
      SNESLineSearch   linesearch
      KSP              ksp
      PC               pc
      Mat              petA 
      Vec              petx, petb, petr, petfa
      PetscViewer      view
      PetscReal        zero, norm, negone
      PetscBool       issymm, isH, iszero, flg
      PetscScalar      sum,  pfive, atmp
      PetscScalar, pointer :: xx_v(:), xt_v(:)
      double precision   prvec(50000)
      PetscInt         its, itrlen, itksp, n
      character (len=4) writeItr
      PCType           pcname
      PCSide           side
      PetscMPIInt      rank, size  
      PetscErrorCode  ierr 
      
      external petFormFunction
     
!  
      negone = -1.0
! The first part is to be identical with other nonlinear solvers in Polyde.
! Solve assuming it's a linear system to get the initial x for starting the 
! nonlinear solver.
! For now the first linear solver uses Polyde default, not PETSc
! The global variable x is stored, so we copy x to petx
! xn is irrelevant as input for this case
      call solve_linear3D( epsgl, resgl, xn, .false., res_iter)
! From here on we use PETSc routine mostly
!
! This version call snes to solve one iteration, then loop again with updated formfunction
! until converge
!

      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      if (ierr .ne. 0) then
        print*,'Unable to initialize PETSc'
        stop
      endif
! Set Vec petx and copy x values to petx
! create vector and set the type
      call VecCreate(PETSC_COMM_WORLD, petx, ierr)
      call VecSetSizes(petx ,PETSC_DECIDE, numdof, ierr)
      call VecSetType(petx, VECSEQ, ierr)

! fill in vector
      print*, 'Start PETSc petx'
      do i=1, numdof
         call VecSetValue(petx, i-1, x(i), INSERT_VALUES, ierr);
      end do
! petx copied confirmed
!      call VecView(petx, PETSC_VIEWER_STDOUT_WORLD, ierr)
 
     call zeit('Building PETSc petx')

!-------------------------------------------------------------------------------
! Start building the system matrix petA
!-------------------------------------------------------------------------------
!
! Jacobian is irrelevant here
      jacobi = .false.
! matvar should be read by setting. For now set as false
      matvar = .false.
! assume for now all matrices are nonsymmetric. To be read from settings
      symmetric = .false.
! Default now is lcsr only 
      csr = .false. 

      call assembly3D(lower, upper, diag, acsr, rhs, ia, ja, jacobi, csr,       &
           &                matvar, symmetric)


! assembly3D by default provides in LCSR and complex datatype. Need to first convert to
! real datatype, then to CSR
      call lcsr2lcsr(diag, lower, upper, rhs, symmetric, diag_DP, lower_DP, upper_DP, rhs_DP)
      call lcsr2csr(diag_DP, lower_DP, upper_DP, symmetric, ia, ja, a_aux_DP, ia_aux, ja_aux)  
!
! We can now convert to PETSc data structure petA
      allocate(nnzrow(numdof))
! get number of columns per row first
      do i=1,numdof
         nnzrow(i) = ia_aux(i+1) - ia_aux(i)
      enddo
!
! size of csr/ nnzeros:
      nnz = ia_aux(numdof+1)-1
      nz=maxval(nnzrow)
!
      call MatCreateSeqAIJ(PETSC_COMM_SELF, numdof, numdof, nz, nnzrow, petA, ierr)
!
      call zeit('start of PETSc')
      nzix=1
      do i=1,numdof
       do j=ia_aux(i),ia_aux(i+1)-1
        atmp = a_aux_DP(nzix)
        call MatSetValue(petA, i-1, ja_aux(nzix)-1, atmp, INSERT_VALUES, ierr)
        nzix=nzix+1
       end do
      end do
      call zeit('Building PETSc Array')

      call MatAssemblyBegin(petA, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(petA, MAT_FINAL_ASSEMBLY, ierr)
!      call MatView(petA, PETSC_VIEWER_STDOUT_WORLD, ierr)
!     
! check symmetry :
      zero = 1e-16
      call MatIsSymmetric(petA, zero, issymm, ierr)
      if (issymm == 1 ) then
        print *, 'PETSc: Matrix A is symmetric'
      else
         print *, 'PETSc WARNING: Matrix A is NOT symmetric!'
      end if

! create vector and set the type for petb
      call VecCreate(PETSC_COMM_WORLD, petb, ierr)
      call VecSetSizes(petb ,PETSC_DECIDE, numdof, ierr)
      call VecSetType(petb, VECSEQ, ierr)

! fill in vector
!
      print*, 'Start PETSc RHS'
      do i=1,numdof
         call VecSetValue(petb, i-1, rhs_DP(i), INSERT_VALUES, ierr);
      end do
      call zeit('Building PETSc RHS')

! finalize
      print*, 'FINALIZE PETB'
      call VecAssemblyBegin(petb, ierr)
      call VecAssemblyEnd(petb, ierr)
      call VecDuplicate(petb, petr, ierr) 
! print values
!     print *, 'PETB RHS:'
!      call VecView(petb, PETSC_VIEWER_STDOUT_SELF, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!         Create the nonlinear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Initialize with petx with linear solution
      call KSPCreate(PETSC_COMM_SELF, ksp, ierr)
      call KSPSetOperators(ksp, petA, petA, ierr)
!      KSPSetOperators as above
!      call KSPSetFromOptions(ksp,ierr)
      call KSPSolve(ksp, petb, petx, ierr)
! petx should now have the linear solution values.
print *, '=========== DONE KSP ==========='
print *, 'CHECK RESIDUAL...'
      call VecDuplicate(petb, petfa, ierr) ! create petfa vec context 
      call MatMult(petA, petx, petfa, ierr) ! Returns (2)
      call VecAXPY(petfa, negone, petb, ierr)
      call VecNorm(petfa, NORM_2, norm, ierr)
      print *, 'PETFA NORM  ', norm

!  Create nonlinear solver context
!  _world and _SELF seem no difference
      call SNESCreate(PETSC_COMM_SELF,snes, ierr)
      call SNESGetKSP(snes,ksp,ierr)
!
!  Create a linear solver routine computing solution       !     call SNESGetKSP(snes,ksp,ierr)      
!      call SNESSetKSP(snes, ksp, ierr)
      call KSPGetPC(ksp,pc,ierr)
      call KSPSetType(ksp,KSPMINRES,ierr)
!      call KSPSetOperators(ksp,petA,petA,ierr)
      call KSPSetResidualHistory(ksp, prvec, 50000, PETSC_TRUE, ierr)

!!$      print *, 'IN SNES BEFORE ITERATION, X BC:'
!!$      do i=1,2
!!$         print *, real(x(i))
!!$      enddo
!!$      print *, 'IN SNES BEFORE ITERATION, X NOT BC:'
!!$      do i=61,62
!!$         print *, real(x(i))
!!$      enddo

iter = 1
do while (iter .le. 20)
!do while (iter .le. iterMax)
      print *, 'IN SNES ITERATION, ITER:', iter
!      print *, 'IN SNES ITERATION, PETX:'
!      call VecView(petx, PETSC_VIEWER_STDOUT_WORLD, ierr)
!
! Set form function F = A*x - b
! test passing petb
      call SNESSetFunction(snes, petr, petFormFunction, petb, ierr)
!
! We don't use Newton based solver yet, so we need not the Jacobian. 
! We try first with the Quasi-Newton method with BFGS variant.
! Q-N method approximates the Jacobian from PetA. Testing:
! SetJacobian does not make a difference. So remove for now
!      call SNESSetJacobian(snes, petA, petA, PETSC_NULL_FUNCTION, 0, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                      Solve the nonlinear system
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! initial guess
      print *, 'Start SNESSolver ...'
      call SNESSetSolution(snes, petx, ierr)
      call SNESGetLineSearch(snes, linesearch, ierr)
      call SNESLineSearchSetType(linesearch, 'cp', ierr)
! SNESQN SNESGMRES now ok but does not converge
! SNESNCG diverge
! SNESNGS, SNESANDERSON, SNESRICHARDSON, 
! From doc: QN solves for the solution of F(x) = b
! Ref: P Brune, et al., "Composing Scalable Nonlinear Algebraic Solvers", 
!      SIAM Review, 57(4), 2015
      call SNESSetType(snes, SNESQN, ierr) 
      call SNESSetFromOptions(snes, ierr)

! TO DO: Set nonlinear iteration to ONE
! From options:
!	-snes_max_it <max_it>	- maximum number of iterations
!       -snes_max_funcs <max_funcs>	- maximum number of function evaluations
! Confirmed to work:
! -snes_max_it 1 -snes_max_funcs 1
!
!      call SNESSetIterationNumber(snes, 1, ierr) ! doesn't seem to affect
!
! For F(X) = A(x)*x - b, no need for petb, as in ex1f1.F90
      call SNESSolve(snes, PETSC_NULL_VEC, petx, ierr)
!
! Both seems to be the same. Still large error in solution
! to remove readme status. not working
!      call VecLockReadPop(petx, ierr)

      call zeit('PETSc matrix solver')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                     Check solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! print values
!      call VecView(petx, PETSC_VIEWER_STDOUT_WORLD, ierr)

! Update petx to include BC values
      do elem = 1, numv
        do f = 1, 4
          if (vv(f,elem).lt.0) then
            surf = abs(vv(f,elem))
            bc_index = sbc(surf)
            do inat = 1, nnat
              if (bctype(bc_index,inat).lt.200) then
                ! TODO: also consider dof on edges and on faces
                do nd = 1, 3 ! Apply bc to node-dof
                  pval = 0._DPC
                  xyzs(1:3) = nod(1:3,vn(f2n(nd,f),elem))
                  call getbcval( sbc( abs(vv(f,elem))), inat, xyzs, pval ) 
                  call VecSetValue(petx, vgdof(elem,inat)%d(f2n(nd,f)) - 1, pval, INSERT_VALUES, ierr)
                end do
              end if
            end do
          end if
        end do ! face
      end do ! elem

! Update petb also
! To call assembly3D?
      call assembly3D(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr,       &
           &                matvar,symmetric)
      if (symmetric) upper => lower
! clean up unused arrays
!        deallocate( ep, edgebc, sfvertbc )
!        nullify( ep, edgebc, sfvertbc )
! assembly3D by default provides in LCSR and complex datatype. Need to first convert to
! real datatype, then to CSR
      call lcsr2lcsr(diag, lower, upper, rhs, symmetric, diag_DP, lower_DP, upper_DP, rhs_DP)
      call lcsr2csr(diag_DP, lower_DP, upper_DP, symmetric, ia, ja, a_aux_DP, ia_aux, ja_aux)  
!
! We can now convert to PETSc data structure petb
! fill in vector
!
      print*, 'Start PETSc RHS'
      do i=1,numdof
!        if ( i == 121 ) print *, 'IN SNES ITERATION, RHS121 NOTBC:', iter-1, real(rhs_DP(121)) 
        call VecSetValue(petb, i-1, rhs_DP(i), INSERT_VALUES, ierr);
      end do
      call zeit('Building PETSc RHS')

! CHECK NORM OF RHS AGAIN
      call MatMult(petA, petx, petfa, ierr) ! Returns (2)
      call VecAXPY(petfa, negone, petb, ierr)
      call VecNorm(petfa, NORM_2, norm, ierr)
      print *, 'PETFA NORM  ', norm

! finalize
      print*, 'FINALIZE PETB'
      call VecAssemblyBegin(petb, ierr)
      call VecAssemblyEnd(petb, ierr)

! print values
!      print *, 'IN SNES: AFTER BC '
!      call VecView(petx, PETSC_VIEWER_STDOUT_WORLD, ierr)
      
      iter = iter + 1
      print *, 'ITERR', iter
      if (iter .eq. iterMax) then
         print *, 'Max. no. of nonlinear iterations reached, solution did not converge'
      endif
!
      call VecGetArrayF90(petx,xx_v,ierr)
      xn = xx_v
! Update the global variable x
      x = xn
!
!      print *, 'IN SNES ITERATION, RHS121 NOTBC:', iter-1, real(rhs_DP(121))
      print *, 'IN SNES ITERATION, X121 NOTBC:', iter-1, real(xn(121))
!      print *, 'IN SNES ITERATION, X088 BC:', iter-1, real(xn(88))
!      do i=61,62
!         print *, real(xn(i))
!      enddo

end do ! while iter

!  Get Vec petx , and store in our defined format: real (DP)
!  Get a pointer to vector data.
!    - For default PETSc vectors, VecGetArray() returns a pointer to
!      the data array.  Otherwise, the routine is implementation dependent.
      print *, 'READ Petsc X ...'
      call VecGetArrayF90(petx,xx_v,ierr)
      xn = xx_v

!  MUST call VecRestoreArrayF90(..) when no longer need access to the array.
       call VecRestoreArrayF90(petx,xx_v,ierr)

!      call VecGetArrayF90(petx,xx_v,ierr)
!      xn2 = xx_v

!  MUST call VecRestoreArrayF90(..) when no longer need access to the array.
!       call VecRestoreArrayF90(petx,xx_v,ierr)

!------------------------------------------------------------------
!  Check the error
      call SNESGetLinearSolveIterations(snes,itksp, ierr) 
      call SNESGetIterationNumber(snes,its,ierr)

! print some values
      print *,'==============================================='
      print *, 'PETSc Nonlinear iterations          :', its
      print *,'==============================================='


! write iterations to file; nnzero is needed by itrhstout
!      call getsetting('WRITE ITERATION', writeItr)
      writeItr = 'NO'
      if (writeItr == 'YES') call itrhstout(nnz, its+1, prvec(1:its+1), numdof)

!  Free work space AKA deallocate
      deallocate(nnzrow)
      call VecDestroy(petx, ierr)
      call VecDestroy(petb, ierr)
      call MatDestroy(petA, ierr)
      call SNESDestroy(snes, ierr)


      end subroutine petscNonlinear
      
!==============================================================
      subroutine petFormFunction(snes, petx, petf, petb, ierr)
#include <petsc/finclude/petsc.h>
      use feminterface3D, only: assembly3D, getbcval
      use feminterface, only: lcsr2lcsr, lcsrzeroremover, lcsr2csr
      use globalvariables3D, only: x, nod, numv, vgdof, vn, vv, sbc, nnat, bctype
      use femtypes
      use petscsnes
      implicit none
      PetscScalar, pointer :: xtmp_v(:)
      SNES     snes
      Vec      petx, petf, petb
      PetscErrorCode ierr
!      integer dummy(*)
! ------------------------------------------------------------------------
!
!  FormFunction - Evaluates nonlinear function, F(x) = A(x)*x - b(x).
!
!  Input Parameters:
!  snes - the SNES context
!  petx - input solution vector
!  dummy - optional user-defined context (not used here)
!
!  Output Parameter:
!  petf - function vector
!  As of now: petf is A(x)*x-b
! Problem: Formfunction must be in the form
!  SNESFunction(SNES snes,Vec x,Vec f,void *ctx);
!  I need more for petA, petb      
!  See ex2.c on how to read arrays and vectors from the calling function
!
! Local Polyde variables
      integer (I4B), pointer :: ia(:), ja(:)
      integer (I4B), pointer :: ia_aux(:), ja_aux(:)
      integer (I4B), allocatable :: nnzrow(:)
      integer (I4B),    parameter:: f2n(3,4)=reshape( (/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/) ) 
     integer (I4B)          ::  i, j, nz, nnz, nzix, f, elem, surf, bc_index, inat, nd
      complex (DPC), pointer :: lower(:)=>null(), upper(:)=>null(), diag(:)=>null(), rhs(:)=>null()
      complex (DPC), pointer :: acsr(:)=>null()
      complex (DPC)          :: pval
      real (DP), pointer     :: lower_DP(:)=>null(), upper_DP(:)=>null(), diag_DP(:)=>null(), rhs_DP(:)=>null()
      real (DP), pointer     :: a_aux_DP(:)=>null(), xtmp(:)
      real (DP)              :: xyzs(3)
      logical                :: jacobi, csr, matvar, symmetric
!
! Local petsc variables
      Mat         petA
      Vec         petxt, petxu
      PetscScalar, pointer :: xx_v(:)
      PetscScalar atmp, neg
      PetscReal   negone, norm
      PetscInt    n
      PetscViewer viewer
!
! The F function is in terms of A(x) and b(x). So we must assemble this for every nonlinear iteration 
! This is done with DM in PETSc. We cannot do it this way yet. So call Polyde's assembly routine.
!
      neg=-1.0
      negone = -1.0
! Jacobian is irrelevant here
      jacobi=.false.
! matvar should be read by setting. For now set as false
      matvar=.false.
! assume for now all matrices are nonsymmetric. To be read from settings too
      symmetric=.false.
      csr = .false.
!
! We also have confirmed that petx changes for every call to formfunction
      print *, 'FORMFUNCTION: READ Petsc X ...'
!      call VecDuplicate(petx, petxt, ierr)
!      call VecCopy(petx, petxt, ierr)
!
! After a nonlinear iteration, the current solution vector may have offset the BC data.
! Restore the BC data to petx but by way of temp var petxt.
! --->    Removed, since petx cannot be changed. See previous version
!========================================================================
! checking if xx_v changes when petx changes:
!      x = xx_v ! can't copy. got seg fault
      print *, 'IN FORMFUNCTION BEFORE ASSEMBLY'
!      do i=1,5
!         print *, real(x(i))
!      enddo
!      call VecRestoreArrayF90(petxt, xx_v, ierr)

! Assembly must be done with the new x
! TO DO
! to check where in the code the new x changes the A and b
      call assembly3D(lower, upper, diag, acsr, rhs, ia, ja, jacobi, csr,       &
           &                matvar, symmetric)

! assembly3D by default provides in LCSR and complex datatype. Need to first convert to
! real datatype, then to CSR
      call lcsr2lcsr(diag, lower, upper, rhs, symmetric, diag_DP, lower_DP, upper_DP, rhs_DP)
!      call lcsrzeroremover(diag_DP, lower_DP, upper_DP, ia, ja, symmetric, 2, 1.e-4_DP, rhs_DP)
      call lcsr2csr(diag_DP, lower_DP, upper_DP, symmetric, ia, ja, a_aux_DP, ia_aux, ja_aux)  

! We can now convert to PETSc data structure again
!
! Get n from the size of petx
      call VecGetSize(petx, n, ierr)
      allocate(nnzrow(n))
! get number of columns per row first
      do i=1,n
         nnzrow(i)=ia_aux(i+1)-ia_aux(i)
      enddo

! size of csr/ nnzeros:
      nnz=ia_aux(n+1)-1
      nz=maxval(nnzrow)
!
! set A to be symmetric, supposed to be faster ??
      call MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,nz,nnzrow,petA,ierr)
      call MatSetOption(petA,MAT_SYMMETRIC,PETSC_TRUE,ierr)
!
      nzix=1
      do i=1,n
       do j=ia_aux(i),ia_aux(i+1)-1
!        atmp = 0.0 ! test zero A
        atmp = a_aux_DP(nzix)
        call MatSetValue(petA, i-1, ja_aux(nzix)-1, atmp, INSERT_VALUES, ierr)
        nzix=nzix+1
       end do
      end do
      call zeit('Building PETSc nonlinear array')
!
      call MatAssemblyBegin(petA, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(petA, MAT_FINAL_ASSEMBLY, ierr)

! Compute function for F(x):
! Unclear if  F(x) = A(x)*x - b (1) or F(x) = A(x)*x (2)
! It seems from tutorial ex1.c F(x) should follow (1)
! MatMult: (out) petf = petA*petx:

      call MatMult(petA, petx, petf, ierr) ! Returns (2)
!
! If uncomment, returns 
! VecAXPY: petf = negone*petb  +  petf
! Ongoing problem: the RHS changes because of solution too, currently
      call VecAXPY(petf, negone, petb, ierr)

! print values

!      call VecView(petf, PETSC_VIEWER_STDOUT_WORLD, ierr)
!      call VecNorm(petf, NORM_2, norm, ierr)
!      print *, 'PETF NORM  ', norm
      return
      end

