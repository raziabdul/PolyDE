      subroutine assemblyandsolvetrans(epsgl,resgl)
      use feminterface, only: assemblytrans, transientsolver, reallocate,  &
      &  zeit, matrhsout, getsetting
      use femtypes
      use globalvariables, only: x, dxdt, d2xdt2, ndof
      use mpmodule
      implicit none
      real (DP)   :: epsgl, resgl
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
!    $Revision: 1.5 $
!    $Date: 2011/06/23 12:45:48 $
!    $Author: chvokas $
!
!  local variables
      integer (I4B) :: stepcount, stepsnumber, returncount
      integer (I4B) :: ndofold
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP) eps, accfac
      real (DP) :: starttime, timestep, endtime, currenttime
      real (DP) :: trangtol = 0._DP, tranltollo = 0._DP, tranltolup = 0._DP
      complex (DPC), pointer :: acsr(:), acsrs(:), acsrd(:), acsrm(:), rhs(:), lower(:), upper(:), diag(:)
      complex (DPC), pointer :: d2xdt2_at_tminusdt(:)
      character (len=25) wrtoption, csroption, ascoption, solver, matrixtype, nonlinear
      character (len=16) :: transsolver, nonlinsolver
      logical jacobi, csr, matvar, ascii, symmetric
      logical timestepreject
!
! Get options from FEMsettings.txt
      call getsetting('CSRFORMAT',csroption)
!
! Option for Compact Storage
      if (csroption .eq. 'CSR') then
        csr=.true.
! Default is the lower CSR
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

!  initialize first order time derivative of x
      if (.not. associated(dxdt)) then
        allocate(dxdt(ndof))
        dxdt = 0._DPC
      end if
!  initialize second order time derivative of x
      if (.not. associated(d2xdt2)) then
        allocate(d2xdt2(ndof))
        d2xdt2 = 0._DPC
      end if
!
      call getsetting('MATRIXTYPE',matrixtype)
      if (matrixtype .eq. 'SYMMETRIC') then
        symmetric=.true.
      else
        symmetric=.false.
      end if
!
      call getsetting('NONLINEAR',nonlinear)
      call getsetting('NONLINEARSOLVER',nonlinsolver)
      select case (nonlinsolver)
        case ('NOJACOBIAN','MOD_NEWRAPH','ARC_LENGTH','ARC_LENGTH_LINS')
          if (nonlinear .ne. 'YES') then
            print*,' **************! Please look in FEMsettings.txt !**************'
            print*,' For this nonlinear method: ',nonlinsolver,' set NONLINEAR = YES'
            return
          end if
        case default
!  In the case of transient linear problem or nonlinear without iteration
!  (updates matrices only once for every time step) do the assembly for the first time here
          call assemblytrans(lower,upper,diag,acsrs,acsrd,acsrm,rhs,ia,ja,jacobi,csr,  &
          &                  matvar,symmetric)
      end select
!
!  Get start and end times for time integration
      call getsetting('START_TIME',starttime)
      call getsetting('INIT_TIME_STEP',timestep)
      call getsetting('END_TIME',endtime)
!
!  Get time integration tolerances
      call getsetting('TRANSIENTSOLVER',transsolver)
      select case (transsolver)
        case ('NEWMARK')
          allocate(d2xdt2_at_tminusdt(size(d2xdt2)))
          d2xdt2_at_tminusdt = 0._DPC
          call getsetting('TRAN_GTOL',trangtol)
          call getsetting('TRAN_LTOL_LO',tranltollo)
          call getsetting('TRAN_LTOL_UP',tranltolup)
        case default
          ! Do nothing
      end select
!
      stepcount = 0
      returncount = 0
      timestepreject = .false.
      currenttime = starttime+timestep
!  Time integration loop
      do while (currenttime .le. endtime)
        stepcount = stepcount + 1
        call transientsolver(acsr,acsrs,acsrd,acsrm,rhs,ia,ja,d2xdt2_at_tminusdt,stepcount,returncount, &
        &    trangtol,tranltollo,tranltolup,currenttime,timestep,timestepreject,epsgl,resgl)
!  returncount checks to see whether we should exit the time integration loop
        if (returncount .eq. 2) exit
!  Reassemble global matrices in the presence of non-linearities
        if (nonlinear .eq. 'YES') then
          select case (nonlinsolver)
            case ('NOJACOBIAN','MOD_NEWRAPH','ARC_LENGTH','ARC_LENGTH_LINS')
              ! Do nothing, the effective stiffness is updated during the nonlinear iteration
            case default
              call assemblytrans(lower,upper,diag,acsrs,acsrd,acsrm,rhs,ia,ja,jacobi,csr,  &
              &                  matvar,symmetric)
          end select
        end if
      end do
!
!  Deallocate d2xdt2_at_tminusdt in case of Newmark method
      if (transsolver .eq. 'NEWMARK') then
        deallocate(d2xdt2_at_tminusdt)
      end if
!
      call zeit('Solving the matrix')
      return
      end subroutine assemblyandsolvetrans
