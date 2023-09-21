      subroutine linear_sparse_solve(numdof, lower_DPC, upper_DPC, diag_DPC, rhs_DPC, x, ia, ja, epsgl, resgl)
      use feminterface, only: lcsr2lcsr, lcsrzeroremover, lcsr2csr, lcsr2coo, lcsr2csc, &
     &  solve2, pardiso_solver, umfsolver, sort_asc_order, getsetting, petscsolver

! ignore umfcsolver
      use femtypes
      implicit none   

      integer(I4B) :: numdof
      integer(I4B), pointer :: ia(:), ja(:)
      real(DP) :: epsgl, resgl
      complex(DPC), pointer :: lower_DPC(:), upper_DPC(:), diag_DPC(:), rhs_DPC(:), x(:)
      intent(in) :: numdof
      intent(out) :: epsgl, resgl
      intent(inout) :: lower_DPC, upper_DPC, diag_DPC, rhs_DPC, ia, ja, x

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
!  solve the sparse linear equation system
!
!   Input:
!         matrix on input is in LCSR format
!
!    numdof      number of equations
!    lower_DPC   lower triangular matrix (compact storage)
!    upper_DPC   upper triangular matrix (compact storage)
!    diag_DPC    diagonal of the matrix
!    rhs_DPC     right hand side vector
!    x           solution vector
!    ia          compact storage information
!    ja          compact storage information
!
!    eps         accuracy to achieve (settings)
!
!   Output:
!    epsgl   reached error
!    resgl   residual


!  local variables
      real(DP) :: eps
      real(SP) :: epsgl_SP, resgl_SP
      logical :: symmetric
      character (len=25) solver, matrixtype, csroption
      character (len=3) :: datatype
      integer(I4B), pointer :: ia_aux(:)=>null(), ja_aux(:)=>null()
      real(DP), pointer :: a_aux_DP(:)=>null()
      real(SP), pointer :: a_aux_SP(:)=>null()
      real(DP), pointer :: diag_DP(:)=>null(), lower_DP(:)=>null(), upper_DP(:)=>null(), rhs_DP(:)=>null(), x_DP(:)=>null()
      real(SP), pointer :: diag_SP(:)=>null(), lower_SP(:)=>null(), upper_SP(:)=>null(), rhs_SP(:)=>null(), x_SP(:)=>null()
!      complex(SPC) , pointer :: diag_SPC(:)=>null(), lower_SPC(:)=>null(), upper_SPC(:)=>null(), rhs_SPC(:)=>null(), x_SPC(:)=>null()
      complex(DPC), pointer :: a_aux_DPC(:)=>null(), x_DPC(:)=>null()


      call getsetting('CSRFORMAT',csroption)
      if (csroption .eq. 'CSR') then
        print*,'CSRFORMAT must be LCSR, CSR is no more allowed'
        pause
        stop
      end if
!-

!   TO DO=========  we may use analysis type from "physics_quantities"  ====================================
!  analysis type can be:
!  SS = steady state,   TH = time-harmonic   or   TR = transient
!  it is SS or TR we have a real valued matrix equation system
!   TO DO=========  we may use analysis type from "physics_quantities"  ====================================
      call getsetting('LINSYS_DATA_TYPE',datatype)

      call getsetting('LINSOLVER_ERROR',eps)

      call getsetting('MATRIXTYPE',matrixtype)
      if (matrixtype .eq. 'SYMMETRIC') then
        symmetric=.true.
      else
        symmetric=.false.
      end if

      select case (datatype)
      case('DP')
        allocate(x_DP(size(x)))
        x_DP = x
        deallocate(x)
!  copy to real valued matrix equation system
        call lcsr2lcsr(diag_DPC,lower_DPC,upper_DPC,rhs_DPC,symmetric,diag_DP,lower_DP,upper_DP,rhs_DP)
!  remove zero and near to zero entries form the matrix
!   TO DO=========  properties should be obtained from the settings file  ====================================
        call lcsrzeroremover(diag_DP, lower_DP, upper_DP, ia, ja, symmetric, 2, 1.e-4_DP, rhs_DP)
      case('SP')
        allocate(x_SP(size(x)))
        x_SP = x
        deallocate(x)
!  copy to real valued matrix equation system
        call lcsr2lcsr(diag_DPC,lower_DPC,upper_DPC,rhs_DPC,symmetric,diag_SP,lower_SP,upper_SP,rhs_SP)
!  remove zero and near to zero entries form the matrix
!   TO DO=========  properties should be obtained from the settings file  ====================================
        call lcsrzeroremover(diag_SP, lower_SP, upper_SP, ia, ja, symmetric, 2, 1.e-4_DP, rhs_SP)
!      case('SPC')
      case default              ! default is double precision complex
!  remove zero and near to zero entries form the matrix
!   TO DO=========  properties should be obtained from the settings file  ====================================
        x_DPC => x
!   TO DO=========  LCSRREMOVER for the complex case has to implemented  ====================================
!        call lcsrzeroremover(diag_DPC, lower_DPC, upper_DPC, ia, ja, symmetric, 2, 1.e-4_DP, rhs_DPC)
      end select

      call getsetting('LINSOLVERTYPE',solver)

      !  convert storage format if needed
      select case(solver)
      case ('UMF')
!  UMF  requires data to be in COO format
        select case (datatype)
        case('DP')
          call lcsr2coo(diag_DP, lower_DP, upper_DP, symmetric, ia, ja, a_aux_DP, ia_aux, ja_aux)
        case('SP')
          call lcsr2coo(diag_SP, lower_SP, upper_SP, symmetric, ia, ja, a_aux_SP, ia_aux, ja_aux)
!        case('SPC')
        case default              ! default is double precision complex
          call lcsr2coo(diag_DPC, lower_DPC, upper_DPC, symmetric, ia, ja, a_aux_DPC, ia_aux, ja_aux)
        end select

     case ('UMFC')
!  UMFC  requires data to be in CSC format
        select case (datatype)
        case('DP')
          call lcsr2csc(diag_DP, lower_DP, upper_DP, symmetric, ia, ja, a_aux_DP, ia_aux, ja_aux)
!        case('SP')
!        case('SPC')
        case default              ! default is double precision complex
          call lcsr2csc(diag_DPC, lower_DPC, upper_DPC, symmetric, ia, ja, a_aux_DPC, ia_aux, ja_aux)
        end select

      case ('PARDISO')
! PARDISO  requires CSR
        select case (datatype)
        case('DP')
          call lcsr2csr(diag_DP, lower_DP, upper_DP, symmetric, ia, ja, a_aux_DP, ia_aux, ja_aux)
        case default              ! default is double precision complex
          call lcsr2csr(diag_DPC, lower_DPC, upper_DPC, symmetric, ia, ja, a_aux_DPC, ia_aux, ja_aux)
        end select
! receive a_aux_XX, ia_aux, ja_aux

      case ('PETSC')
         call lcsr2csr(diag_DP, lower_DP, upper_DP, symmetric, ia, ja, a_aux_DP, ia_aux, ja_aux)
!        receive a_aux_DP, ia_aux, ja_aux
     end select  ! solver

!  call the solver
      call getsetting('LINSOLVER_ERROR',eps)
      select case(solver)
      case ('SSORCG')
        select case (datatype)
        case('DP')
          call solve2(lower_DP,upper_DP,diag_DP,rhs_DP,x_DP,numdof,eps,ia,ja,             &
     &                epsgl,resgl,symmetric)
          if (symmetric) then
            deallocate( ia, ja, lower_DP, diag_DP, rhs_DP)
            nullify( ia, ja, lower_DP, diag_DP, rhs_DP)
          else 
            deallocate( ia, ja, lower_DP, upper_DP, diag_DP, rhs_DP)
            nullify( ia, ja, lower_DP, upper_DP, diag_DP, rhs_DP)
          end if
        case('SP')
          call solve2(lower_SP,upper_SP,diag_SP,rhs_SP,x_SP,numdof,real(eps,SP),ia,ja,    &
     &                epsgl_SP,resgl_SP,symmetric)
          epsgl = epsgl_SP
          resgl = resgl_SP
          if (symmetric) then
            deallocate( ia, ja, lower_SP, diag_SP, rhs_SP)
            nullify( ia, ja, lower_SP, diag_SP, rhs_SP)
          else 
            deallocate( ia, ja, lower_SP, upper_SP, diag_SP, rhs_SP)
            nullify( ia, ja, lower_SP, upper_SP, diag_SP, rhs_SP)
          end if
        case default              ! default is double precision complex
          call solve2(lower_DPC,upper_DPC,diag_DPC,rhs_DPC,x_DPC,numdof,eps,ia,ja,        &
     &                epsgl,resgl,symmetric)
          if (symmetric) then
            deallocate( ia, ja, lower_DPC, diag_DPC, rhs_DPC)
            nullify( ia, ja, lower_DPC, diag_DPC, rhs_DPC)
          else
            deallocate( ia, ja, lower_DPC, upper_DPC, diag_DPC, rhs_DPC)
            nullify( ia, ja, lower_DPC, upper_DPC, diag_DPC, rhs_DPC)
          end if
        end select ! case SSORCG
!        
      case ('PETSC')
        select case (datatype)
        case('DP')
          print *, 'calling petsc solvers library ...'
!          call petscsolver(a_aux_DP,rhs_DP,x_DP,numdof,eps,ia_aux,ja_aux)  ! original
! TEMPORARY, TO DO ACTUAL resgl
          resgl = 0.
          epsgl = 0.
          if (symmetric) then
            deallocate( a_aux_DP, ia_aux, ja_aux, rhs_DP )
            nullify( a_aux_DP, ia_aux, ja_aux, rhs_DP )
          else 
            deallocate( a_aux_DP, ia_aux, ja_aux, rhs_DP )
            nullify( a_aux_DP, ia_aux, ja_aux, rhs_DP )
          end if
        case('SP')
          print *, 'SP TO IMPLEMENT'
!        case default              ! default is double precision real
          print *, 'calling petsc solvers library ...'
!        call petscsolver(a_aux_DP,rhs_DP,x_DP,numdof,eps,ia_aux,ja_aux)  ! original
! TEMPORARY, TO DO ACTUAL resgl
        resgl = 0.
        epsgl = 0.
          if (symmetric) then
            deallocate( a_aux_DP, ia_aux, ja_aux, rhs_DP )
            nullify( a_aux_DP, ia_aux, ja_aux, rhs_DP )
          else 
           deallocate( a_aux_DP, ia_aux, ja_aux, rhs_DP )
           nullify( a_aux_DP, ia_aux, ja_aux, rhs_DP )
         end if
        end select ! case
!a_aux_DP, ia_aux, ja_aux

      case ('PARDISO')
        select case (datatype)
        case('DP')
          call pardiso_solver( a_aux_DP, rhs_DP, x_DP, numdof, ia_aux, ja_aux )
        case default              ! default is double precision complex
          call pardiso_solver( a_aux_DPC, rhs_DPC, x_DPC, numdof, ia_aux, ja_aux )
        end select
      ! to pass a_aux_DPC, ia_aux, ja_aux
  
      case default
        print*,' no such solver: ', solver,' use PARDISO'
      end select

!  copy back to complex
      select case (datatype)
      case('DP')
        allocate (x(size(x_DP)))
        x = x_DP
      case('SP')
        allocate (x(size(x_SP)))
        x = x_SP
        deallocate(x_SP)
!      case('SPC')
      case default
!  DPC is default
        x => x_DPC
      end select

      return
      end

