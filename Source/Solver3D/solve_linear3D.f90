      subroutine solve_linear3D(epsgl, resgl, X_0, only_assembly, res_iter)
      use feminterface, only: reallocate, zeit, matrhsout, getsetting,    &
     &  print_error, palloc, linear_sparse_solve
      use feminterface3D, only: assembly3D, preassemb3D, CSR_residual, export2matlab, LCSR_residual
      use femtypes
      use globalvariables3D, only: ep, edgebc, numdof, numv, ve, x, sfvertbc, vn, nod, numn, nnat, vv
      implicit none
      complex (DPC),  pointer, optional :: X_0(:), res_iter(:)
      real (DP)                         :: epsgl, resgl
      logical, optional :: only_assembly
      intent (in)                       :: X_0, only_assembly
      intent(out)                       :: epsgl, resgl, res_iter

!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!    $Revision: 1.4 $
!    $Date: 2015/06/23 13:31:19 $
!    $Author: jimmykamboh $
!------------------------------------------------------------------------------
!
!  prepare assembly, assemble and solve the system matrix
!
!   Input:
!         X_0   :   Initial Solution Vector (used in nonlinear case)
!
!  local variables
      integer (I4B)          :: ndofold
      integer (I4B), pointer :: ia(:)=>null(), ja(:)=>null()
      complex (DPC), pointer :: lower(:)=>null(), upper(:)=>null(), diag(:)=>null(), rhs(:)=>null()
      complex (DPC), pointer :: acsr(:)=>null()
      complex (DPC), pointer :: residual(:)=>null()
      character (len=25)     :: wrtoption, csroption, ascoption, solver, matrixtype, material_var_opt, nonlinear, export_matlab
      logical jacobi, csr, matvar, ascii, symmetric

!__
! Checks and preparations:
      call getsetting('CSRFORMAT',csroption)
      if (csroption .eq. 'CSR') then
        print*,'CSRFORMAT must be LCSR, CSR is no more allowed'
        pause
        stop
      else
        csr = .false.
      end if

      if (present(X_0)) then
        if (size(X_0).ne.numdof) then
          call print_error(5,'Can not solve Problem, dimension mismatch of the initial Solution')
        end if
      end if
!-
      call getsetting('MATRIXTYPE',matrixtype)
      if (matrixtype .eq. 'SYMMETRIC') then
        symmetric=.true.
      else
        symmetric=.false.
      end if
!-
      call getsetting('WRITEMATRIX',wrtoption)
!-
      call getsetting('LINSOLVERTYPE',solver)

      call getsetting('NONLINEAR',nonlinear)
!-
      jacobi=.false.
      call getsetting('MATER_VARIATION',material_var_opt)
      if (material_var_opt .eq. 'YES') then
        matvar=.true.
      else
        matvar=.false.
      end if
!-
      if (.not. associated(x)) then
        call palloc(x,numdof)
      else
        if (size(x) < numdof) then
          ndofold = size(x)
          x=>reallocate(x,numdof)
          if (present(X_0)) then
            x(ndofold+1:numdof) = X_0(ndofold+1:numdof)
          else
            x(ndofold+1:numdof) = 0._DPC
          end if
        else if (size(x) > numdof) then
          deallocate(x)
          call palloc(x,numdof)
        end if
      end if
!- Write X_0 to x ( for the non-linear case)
      if (present(X_0)) then
        x = X_0
      else
        x = 0._DPC
      end if
!-
      call getsetting('EXPORT2MATLAB',export_matlab)
      if (export_matlab.eq.'YES') call export2matlab('SOLUTION')
!__
! Start:
      !if (nonlinear .ne. 'YES') then ! #### EDITED
      ! ! if (.not. associated(ep)) then
      !  ! 1) Perform Pre-assembly of Element Matrices
      !    print*,'Pre - Assembly'
      !    call preassemb3D
      !    call zeit('')
      ! ! end if
      !else 
      !  print*, 'Skipped Pre-Assembly'
      !end if ! ### EDITED
!__ 

!__
! 3) Assemble the global system of equations
!    and clean up unused global variables before solving


      call assembly3D(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr,       &
           &                matvar,symmetric)
      if (symmetric) upper => lower

      if (nonlinear .ne. 'YES') then ! #### EDITED
        deallocate( ep, edgebc, sfvertbc )
        nullify( ep, edgebc, sfvertbc )
      end if

!################## Edit
!  TO Do this should only be called if we have nonlinear iteration
      if (csr) then
        call palloc(residual,numdof)
        call CSR_residual(acsr,rhs,x,ia,ja,numdof,residual)
      else
        call palloc(residual,numdof)
        call LCSR_residual(lower,upper,diag,rhs,x,ia,ja,numdof,residual)
      end if
!################## End Edit
      if (present(only_assembly)) then
        if (only_assembly .eq. .true.) then
          if (csr .eq. .true.) then
            deallocate (acsr,rhs,ia,ja)
          else
! LCSR format
            if (symmetric) then
              deallocate (lower,diag,rhs,ia,ja)
            else
              deallocate (lower,upper,diag,rhs,ia,ja)
            end if
          end if
          return
        end if
      end if
! 'BUG RIGHT BELOW +++++'
! residual should always be present
! print *, 'RES_ITER',present(res_iter) ! True, so res_iter is present
! but compiler says not allocated
! so try:
      if (present(res_iter)) allocate(res_iter(numdof))

!  TO DO this should only be called if we have nonlinear iteration
      if (.not.present(res_iter)) then
        deallocate(residual)
      else
        deallocate(res_iter)
        res_iter=>residual
        nullify(residual)
      end if
!################## Edit
      !diag = diag + 0
      !rhs = rhs + 0*x
!################## End Edit
!__
! 4) Write system matrix and RHS to file if requested
      if (wrtoption .eq. 'YES') then
        call getsetting('MATRIXOUTPUT',ascoption)
        if (ascoption .eq. 'ASCII') then
          ascii=.true.
        else
          ascii=.false.
        end if
        call matrhsout(numdof,ia,ja,acsr,rhs,ascii)
        call zeit('Writing matrix to file...')
      end if
!__
! 5) Solve global FEM: find x for Ax=b
      call linear_sparse_solve(numdof, lower, upper, diag, rhs, x, ia, ja, epsgl, resgl)
!-
      call zeit('solving the matrix')
      return
!
!_End.
      end subroutine solve_linear3D
