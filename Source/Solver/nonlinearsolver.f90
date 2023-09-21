      subroutine nonlinearsolver_newmark(acsr,x_at_t,csrvec,ia,ja,newmark_params,epsgl,resgl)
      use feminterface, only: getsetting, pardiso_solver, assemblytrans, amux
      use femtypes
      use globalvariables, only: x, dxdt, d2xdt2, ndof
      implicit none
      complex (DPC), pointer :: acsr(:), x_at_t(:), csrvec(:)
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP) :: newmark_params(8),epsgl, resgl
      intent (in) :: newmark_params
      intent (out) :: epsgl, resgl
!
!-------------------------------------------------------------------------------
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
!    $Revision: 1.6 $
!    $Date: 2014/06/25 08:28:32 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
!  local variables:
      integer (I4B) :: nonlinstepscount, nonlinitersteps
      real (DP) :: eps, nonlinsteperror, nonlinerror
      character (len=16) :: nonlinsolver, solver
      character (len=25) :: matrixtype, csroption
      complex (DPC), pointer :: acsrs(:), acsrd(:), acsrm(:), rhs(:), lower(:), upper(:), diag(:)
      complex (DPC), pointer :: csrvec2(:)
      logical jacobi, csr, matvar, symmetric
!  For the computation of the energy norm
      complex (DPC), pointer :: x_current(:), Ax(:)

      print*,'nonlinearsolver_newmark is not (no more) available'
!
!! Get options from FEMsettings.txt
!      call getsetting('MATRIXTYPE',matrixtype)
!      if (matrixtype .eq. 'SYMMETRIC') then
!        symmetric=.true.
!      else
!        symmetric=.false.
!      end if
!      call getsetting('CSRFORMAT',csroption)
!      if (csroption .eq. 'CSR') then
!        csr=.true.
!      else
!        csr=.false.
!      end if
!      jacobi=.false.
!      matvar=.false.
!!
!      call getsetting('NONLINEARSOLVER',nonlinsolver)
!      call getsetting('NONLIN_ERROR',nonlinerror)
!      nonlinsteperror = nonlinerror + 0.1_DP
!      nonlinstepscount = 0
!!
!      select case (nonlinsolver)
!!
!        case ('NOJACOBIAN')
!!  Nonlinear iteration without Jacobian computation
!          call getsetting('NONLINITERSTEPS',nonlinitersteps)
!          allocate(x_current(size(x)))
!          x_current = 0._DPC
!          allocate(Ax(size(x)))
!          Ax = 0._DPC
!!  Nonlinear iteration loop
!          do while ((nonlinsteperror .gt. nonlinerror) .and. (nonlinitersteps .gt. nonlinstepscount))
!!  Assemble global matrices
!            call assemblytrans(lower,upper,diag,acsrs,acsrd,acsrm,rhs,ia,ja,jacobi,csr,  &
!            &                  matvar,symmetric)
!!  Calculate effective stiffness matrix
!!  newmark_params = (/ a0,a1,a2,a3,a4,a5,a6,a7 /)
!            allocate(acsr(size(acsrs)))
!            acsr = (0._DP,0._DP)
!            acsr = acsrs + newmark_params(1)*acsrm + newmark_params(2)*acsrd
!!  Calculate effective loads for this time step
!            allocate(csrvec(size(rhs)),csrvec2(size(rhs)))
!            csrvec = 0._DPC
!            csrvec2 = 0._DPC
!            call amux(ndof,(newmark_params(2)*x_at_t + newmark_params(5)*dxdt + newmark_params(6)*d2xdt2),csrvec,acsrd,ja,ia)
!            call amux(ndof,(newmark_params(1)*x_at_t + newmark_params(3)*dxdt + newmark_params(4)*d2xdt2),csrvec2,acsrm,ja,ia)
!            csrvec = rhs + csrvec + csrvec2
!            deallocate(csrvec2)
!!
!!  Calculate non-linear iteration error
!            if (nonlinstepscount .gt. 0) then
!              call amux(ndof,(x - x_current),Ax,acsr,ja,ia)
!              nonlinsteperror = abs(real(dot_product(x,Ax)))
!            end if
!!  Save current x
!            x_current = x
!!
!!  Solve for the next x
!            call getsetting('LINSOLVERTYPE',solver)
!            select case(solver)
!            case ('UMF')
!              call umfsolver(acsr,csrvec,x,ndof,eps,ia,ja,epsgl,resgl)
!            case ('UMFC')
!              call umfcsolver(acsr,csrvec,x,ndof,ia,ja)
!              deallocate(ia,ja)
!            case ('PARDISO')
!              call pardiso_solver(acsr,csrvec,x,ndof,ia,ja)
!              deallocate(ia,ja)
!              deallocate(acsr)
!            case default
!              print*,' no such solver: ',solver,' use UMF or UMFC or PARDISO'
!            end select
!!
!            nonlinstepscount = nonlinstepscount + 1
!!  End nonlinear iteration loop
!          end do
!          deallocate(x_current,Ax)
!!
!        case ('MOD_NEWRAPH')
!
!
!        case ('ARC_LENGTH')
!
!
!        case ('ARC_LENGTH_LINS')
!
!
!        case default
!!  Default case is without any nonlinearity iteration, simply solve the system
!          call getsetting('LINSOLVERTYPE',solver)
!          select case(solver)
!          case ('UMF')
!            call umfsolver(acsr,csrvec,x,ndof,eps,ia,ja,epsgl,resgl)
!          case ('UMFC')
!            call umfcsolver(acsr,csrvec,x,ndof,ia,ja)
!            deallocate(ia,ja)
!          case ('PARDISO')
!            call pardiso_solver(acsr,csrvec,x,ndof,ia,ja)
!            deallocate(ia,ja)
!            deallocate(acsr)
!          case default
!            print*,' no such solver: ',solver,' use UMF or UMFC or PARDISO'
!          end select
!!  End select (nonlinsolver)
!      end select
!
      return
      end subroutine nonlinearsolver_newmark


      subroutine nonlinearsolver_stat(acsr,rhs,lower,upper,diag,ia,ja,epsgl,resgl)
      use feminterface, only: getsetting, umfsolver, umfcsolver, assembly, amux, pardiso_solver
      use femtypes
      use globalvariables, only: x, ndof
      implicit none
      complex (DPC), pointer :: acsr(:), rhs(:), lower(:), upper(:), diag(:)
      integer (I4B), pointer :: ia(:), ja(:)
      real (DP) :: epsgl, resgl
      intent (out) :: epsgl, resgl
!
!-------------------------------------------------------------------------------
!
!  local variables:
      integer (I4B) :: nonlinstepscount, nonlinitersteps
      real (DP) :: eps, xEnorm, DifxEnorm, nonlinresEnorm, nonlinerror
      character (len=16) :: nonlinsolver, solver
      character (len=25) :: matrixtype, csroption
      logical jacobi, csr, matvar, symmetric
!  For the computation of the energy norm
      complex (DPC), pointer :: x_current(:), Ax(:), ADifx(:) ! nonlinres(:)
!
! Get options from FEMsettings.txt
      print*,'nonlinearsolver_stat is not (no more) available'

!      call getsetting('MATRIXTYPE',matrixtype)
!      if (matrixtype .eq. 'SYMMETRIC') then
!        symmetric=.true.
!      else
!        symmetric=.false.
!      end if
!      call getsetting('CSRFORMAT',csroption)
!      if (csroption .eq. 'CSR') then
!        csr=.true.
!      else
!        csr=.false.
!      end if
!      jacobi=.false.
!      matvar=.false.
!!
!      call getsetting('NONLINEARSOLVER',nonlinsolver)
!      call getsetting('NONLIN_ERROR',nonlinerror)
!      xEnorm = 1._DP
!      DifxEnorm = xEnorm + 0.1_DP
!      nonlinstepscount = 0
!!
!      select case (nonlinsolver)
!!
!        case ('NOJACOBIAN')
!!  Nonlinear iteration without Jacobian computation
!          call getsetting('NONLINITERSTEPS',nonlinitersteps)
!          allocate(x_current(size(x)))
!          x_current = 0._DPC
!          allocate(ADifx(size(x)))
!          allocate(Ax(size(x)))
!          ! allocate(nonlinres(size(x)))
!          ADifx = 0._DPC
!          Ax = 0._DPC
!          ! nonlinres = 0._DPC
!!  Nonlinear iteration loop
!          do while ((DifxEnorm .gt. (0.001_DP*xEnorm)) .and. (nonlinitersteps .gt. nonlinstepscount))
!!  Assemble global matrices
!            call assembly(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr, &
!            &             matvar,symmetric)
!!  Check if the solution approaches the steady-state
!            if (nonlinstepscount .gt. 0) then
!              call amux(ndof,(x - x_current),ADifx,acsr,ja,ia)
!              DifxEnorm = abs(real(dot_product((x - x_current),ADifx)))
!              call amux(ndof,x,Ax,acsr,ja,ia)
!              xEnorm = abs(real(dot_product(x,Ax)))
!              ! nonlinres = Ax - rhs
!            end if
!!  Save current x
!            x_current = x
!!
!!  Solve for the next x
!            call getsetting('LINSOLVERTYPE',solver)
!            select case(solver)
!            case ('UMF')
!              call umfsolver(acsr,rhs,x,ndof,eps,ia,ja,epsgl,resgl)
!            case ('UMFC')
!              call umfcsolver(acsr,rhs,x,ndof,ia,ja)
!              deallocate(ia,ja)
!            case ('PARDISO')
!              call pardiso_solver(acsr,rhs,x,ndof,ia,ja)
!              deallocate(ia,ja)
!              deallocate(acsr)
!            case default
!              print*,' no such solver: ',solver,' use UMF or UMFC or PARDISO'
!            end select
!!
!            nonlinstepscount = nonlinstepscount + 1
!!  Print statements
!            print*,"Nonliniterstep ===>",nonlinstepscount
!            print*,"DifxEnorm      ===>",DifxEnorm
!            print*,"xEnorm         ===>",xEnorm
!!  End nonlinear iteration loop
!          end do
!          deallocate(x_current,Ax,ADifx) ! nonlinres
!
!
!        case ('MOD_NEWRAPH')
!
!
!        case ('ARC_LENGTH')
!
!
!        case ('ARC_LENGTH_LINS')
!
!
!        case default
!!  If nonlinearsolver is not correctly defined simply assemble and solve the system once
!          call assembly(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr, &
!          &             matvar,symmetric)
!!
!          call getsetting('LINSOLVERTYPE',solver)
!          select case(solver)
!          case ('UMF')
!            call umfsolver(acsr,rhs,x,ndof,eps,ia,ja,epsgl,resgl)
!          case ('UMFC')
!            call umfcsolver(acsr,rhs,x,ndof,ia,ja)
!            deallocate(ia,ja)
!          case ('PARDISO')
!            call pardiso_solver(acsr,rhs,x,ndof,ia,ja)
!            deallocate(ia,ja)
!            deallocate(acsr)
!          case default
!            print*,' no such solver: ',solver,' use UMF or UMFC or PARDISO'
!          end select
!!  End select (nonlinsolver)
!      end select
!!
      return
      end subroutine nonlinearsolver_stat