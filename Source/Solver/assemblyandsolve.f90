      subroutine assemblyandsolve(epsgl,resgl,adaptstep)
      use feminterface, only: assembly, reallocate, nonlinearsolver,  &
      &   zeit, matrhsout, getsetting, csrmmout, linear_sparse_solve
      use femtypes
      use globalvariables, only: x, ndof, p, fem_accuracy
      use mpmodule, only: multiphysics
      implicit none
      integer (I4B), optional :: adaptstep
      real (DP)   :: epsgl, resgl
      intent(in) :: adaptstep
      intent(out) :: epsgl, resgl
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
!    $Revision: 1.24 $
!    $Date: 2014/07/01 16:15:41 $
!    $Author: m_kasper $
!
!  Input:
!            adaptstep   taken during adaptation, needed to write matrix to file
!
!  Output:
!            epsl
!            resgl
!
!  local variables
      integer (I4B) :: ndofold
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP) eps, accfac
      complex (DPC), pointer :: lower(:), upper(:), diag(:), rhs(:), acsr(:)
      character (len=25) wrtoption, csroption, ascoption, solver, matrixtype, nonlinear
      logical jacobi, csr, matvar, ascii, symmetric
!
!  check if a solution is available as starting vector
      if (.not. associated(x)) then
        allocate(x(ndof))
        x = 0._DPC
      else
!  check if length(x) is ndof. if not reallocate to ndof
        if (size(x) < ndof) then
          ndofold = size(x)
          x=>reallocate(x,ndof)
          x(ndofold+1:ndof) = 0._DPC
        else if (size(x) > ndof) then
          deallocate(x)
          allocate(x(ndof))
          x = 0._DPC
        end if
      end if
!
      call getsetting('NONLINEAR',nonlinear)
      select case (nonlinear)
      case('YES')
!  Do not do the assembly here but in the nonlinearsolver routine
        call nonlinearsolver(acsr,rhs,lower,upper,diag,ia,ja,epsgl,resgl)
!
!  Default case is linear
      case default
!  Get options from FEMsettings.txt
        call getsetting('CSRFORMAT',csroption)
!
!  Option for Compact Storage
        if (csroption .eq. 'CSR') then
          csr=.true.
!  Default is the lower CSR
        else
          csr=.false.
        end if
!
        jacobi=.false.
!  material is non homogeneous for multiphysics -> set matvar to .true.
        if (multiphysics) then
          matvar=.true.
        else
          matvar=.false.
        end if
!
        call getsetting('MATRIXTYPE',matrixtype)
        if (matrixtype .eq. 'SYMMETRIC') then
          symmetric=.true.
        else
          symmetric=.false.
        end if
!
        call assembly(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr, &
     &                matvar,symmetric)
!  Write system matrix and RHS to file if requested
        call getsetting('WRITEMATRIX',wrtoption)
        if (wrtoption .eq. 'YES') then
          call getsetting('MATRIXOUTPUT',ascoption)
          ascoption = ascoption(1:len_trim(ascoption))
!
          select case (ascoption)
          case('ASCII')
            ascii=.true.
            call matrhsout(ndof,ia,ja,acsr,rhs,ascii)
          case('MATRIXMARKET')
            print*,' Doing MM matrix '
            if (present (adaptstep)) then
              print*, 'adap ',adaptstep
            endif
            print*, 'nndof',ndof
            call csrmmout(ia,ja,acsr,adaptstep) !no RHS yet
          case default
!  Default is binary
            ascii=.false.
            call matrhsout(ndof,ia,ja,acsr,rhs,ascii)
          end select
!
          call zeit('Writing matrix to file...')
        end if
!
        call getsetting('LINSOLVER_ERROR',eps)
!  for LINSOLVER_ERROR given as zero or negative
        if (eps .le. 0._DP) then
          accfac=sqrt(real(p,DP))*10._DP
!  determine linear solver error for the next solution step
!  dependent on the global solution error (dicretization error)
          eps=min(1.e-2_DP,fem_accuracy/accfac)
        end if
        call getsetting('LINSOLVERTYPE',solver)
!
        if (csr) then
          print*,'CSR matrix format is obsolet use LCSR'
          stop
        else
          call linear_sparse_solve(ndof, lower, upper, diag, rhs, x, ia, ja, epsgl, resgl)
        end if
!
        call zeit('Solving the matrix')
!  End of the select nonlinear structure
      end select
!
      return
      end subroutine assemblyandsolve