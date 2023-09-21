subroutine solve_nonlinear3D()
      use feminterface,        only: print_error, getsetting
      use feminterface3D,      only: solve_andersonr, solve_fixedpt, solve_diag_jacobian 
      use feminterface3D,      only: solve_diag_jacobianid, solve_barbor, solve_barborv, petscNonlinear
      use feminterface3D,      only: solve_andersone, get_starting_solution
      use globalvariables3D,   only: numdof, x
      use femtypes
      implicit none
!
!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 - 2015 Institute for Micro Systems Technology,
!                              Hamburg University of Technology.
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
!
!    $Revision: 1.5 $
!    $Date: 2015/11/10 14:29:35 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
! This subroutine takes initial guess x_0 to compute the solution for nonlinear
! system and returns the final solution. For solution process, it calls the 
! selected nonlinear solver routine for the desired solution process. Routines are
! arranged in alphabatical order below and some routines can possibly be discarded.
! Routines are quite independent of eachother, therefore, changes in one routine
! does not effect another routine.

!-------------------------------------------------------------------------------
! local variables:
!
      complex (DPC), pointer     :: x_0(:)
      character (len=20)         :: nonlinsolver, nonlinear
     
   
! Checks and preparations:
      call getsetting('NONLINSOLVERTYPE',nonlinsolver)
      call getsetting('NONLINEAR',nonlinear)
      

      allocate(x_0(numdof))

      x_0(1:numdof) = 0._DPC
      
      !call get_starting_solution(x_0)
   
      select case (nonlinsolver)
      case ('FIXED_POINT')      !Fixed Point Solver
      print *, 'calling fixed point'
         call solve_fixedpt(x_0)
    
      case ('ANDERSON_EXTRAP')  !for non-generalized Anderson extrapolation solver
         call solve_andersone(x_0)
         
      case ('ANDERSON_RELAX')   !Anderson relaxation solver
         call solve_andersonr(x_0)
         
      case ('BARBOR')           !Barzilai-Borwein solver
         call solve_barbor(x_0)
         
      case ('BARBOR_VAR')       !a variant of Barzilai-Borwein solver with vector parameter in step update
         call solve_barborv(x_0)
         
      case ('DIAG_JACOBIAN')    !Diagonal approx. for Jacobian as in Newton method (based on lin. system residual)
         call solve_diag_jacobian(x_0)
         
      case ('DIAG_JACOBIAN_ID') !Diagonal approx. for Jacobian as in Newton method based upon iteration defect
         call solve_diag_jacobianid(x_0)
         
      case ('PETSCNONLINEAR') !
!         call petscNonlinear(x_0)

      case default
         print*, 'No nonlinear solver chosen, default nonlinear solver is FIXED POINT'
         call solve_fixedpt(x_0)
      end select
       

    x = x_0(:)
    deallocate(x_0)
    return   
end subroutine solve_nonlinear3D




!----------------------------------------------------------------------------------------------!
!--------------------------------------------subroutine Anderson extrapolation solver-----------------------------!
!----------------------------------------------------------------------------------------------!

subroutine solve_andersone(y_n)
      use feminterface,        only: print_error
      use feminterface,        only: getsetting
      use feminterface3D,      only: solve_linear3D, export2matlab, solve_fixedpt
      use globalvariables3D,   only: numdof, x
      use femsettings,         only: readstd
      use femtypes
      implicit none
      complex (DPC),      pointer :: y_n(:)
      intent (inout) :: y_n

!
!--------------------------------
! local variables:
!
      integer (I4B)               :: max_fp_i, i , xi_1, yi_1, xi_2, yi_2, and_i, max_and_i
      integer (I4B)               :: shd, sh_idx
      real    (DP),   allocatable :: val(:,:)
      real     (DP)               :: epsgl,resgl , theeta, ande_tol, and_res, anditer_err, fp_omega
      real     (DP),  allocatable :: diff_r(:), r(:)
      complex (DPC),      pointer :: res_iter(:)
      complex (DPC),      pointer :: x_n(:,:)
      complex (DPC),      pointer :: x_n_ptr(:)
      character (len=20)          :: res_type, export_matlab
      logical                     :: exp_matlb

      exp_matlb = .false.
!__
!
   call getsetting('MAX_FIXPOINT_ITER',max_fp_i)
   call getsetting('MAX_ANDERSON_ITER',max_and_i)
   call getsetting('NONLIN_TOLERANCE',ande_tol)
   call getsetting('FPRELAX_PARAMETER',fp_omega)
   call getsetting('PLOT_RESIDUALTYPE',res_type)
   call getsetting('EXPORT2MATLAB',export_matlab)
   if (export_matlab.eq.'YES') exp_matlb = .true.
   if (exp_matlb) call export2matlab('HEADERS')
   if (exp_matlb) call export2matlab('NODES')
     shd = 4
   
     allocate(diff_r(numdof), res_iter(numdof), r(numdof))
     allocate(x_n(numdof,shd))
     allocate(val(1,2))

     
     if (max_fp_i.lt.2) then  
       call print_error(3,'No. of fixed pt. iterations must be atleast 2')
     end if
   
     x_n(:,1) = y_n(:)
     
     call solve_fixedpt(y_n)
     
     max_fp_i = 4
   
     do i = 1, max_fp_i
        call solve_linear3D(epsgl,resgl,y_n,.false.,res_iter) !this will indeed perform fixed pt. step with input y_n
        y_n(:) = fp_omega*x + (1-fp_omega)*y_n(:)
        if (i.eq.max_fp_i) then
        x_n(:,3) = y_n(:)
        else if ((i.lt.max_fp_i).and.(mod(i,2).eq.1)) then
        x_n(:,2) = y_n(:)
        else if ((i.lt.max_fp_i).and.(mod(i,2).eq.0)) then 
        x_n(:,1) = y_n(:)
        end if
     end do
   
      
      if (mod(max_fp_i,2).eq.0) then
      x_n(:,3:4) = x_n (:,2:3)
      else if (mod(max_fp_i,2).eq.1) then
      x_n(:,4) = x_n(:,3)
      x_n(:,3) = x_n(:,1)
      x_n(:,1) = x_n(:,2)
      x_n(:,2) = x_n(:,3)     
      end if 
      
     print "(A18,I7,A10)" ,'Initial FP steps: ', max_fp_i 
      
     sh_idx = 4 
     and_i = 1
     
     call getsetting('MAX_FIXPOINT_ITER',max_fp_i)
     
     max_fp_i = max_fp_i + 4

     do while (and_i.le.max_and_i)
       
      readstd   = .false.
      call getsetting('FPRELAX_PARAMETER',fp_omega)
     
     
      print "(A18,I7)"    ,'Anderson extrapolation step no.: ', and_i 
      print "(A18,G12.3)" ,'           Relaxation Parameter: ', fp_omega 
     
      
       xi_1 = mod(sh_idx+(shd-4),shd)+1
       yi_1 = mod(sh_idx+(shd-3),shd)+1
       xi_2 = mod(sh_idx+(shd-2),shd)+1
       yi_2 = sh_idx

       diff_r(:) = x_n(:,yi_2) - x_n(:,xi_2) - x_n(:,yi_1) + x_n(:,xi_1) !r_n - r_n-1, difference of last residuals
       r(:) = x_n(:,yi_2) - x_n(:,xi_2)

!__
! 2) Compute theeta - optimal relaxation parameter:
       theeta = dot_product(r(:),diff_r(:))/dot_product(diff_r(:),diff_r(:))
!__
! 3) Update solution
       sh_idx = mod(sh_idx,shd) + 1
       x_n(:,sh_idx) = x_n(:,yi_2) + theeta*(x_n(:,yi_1) - x_n(:,yi_2))
   
       x_n_ptr=>x_n(:,sh_idx)
!__
! 4) Use fixed point iteration to obtain next solution:
       call solve_linear3D(epsgl,resgl,x_n_ptr,.false.,res_iter)
!__
! 5) Compute iteration error:
       y_n(:) = x_n(:,sh_idx)

       sh_idx = mod(sh_idx,shd) + 1
       x_n(:,sh_idx) = fp_omega*x + (1-fp_omega)*y_n(:)
       
       xi_1 = mod(sh_idx+(shd-4),shd)+1
       yi_1 = mod(sh_idx+(shd-3),shd)+1
       xi_2 = mod(sh_idx+(shd-2),shd)+1
       yi_2 = sh_idx
       
       anditer_err    = sqrt(dot_product(x_n(:,yi_2) - x_n(:,xi_2),x_n(:,yi_2) - x_n(:,xi_2)))
       and_res = sqrt(dot_product(res_iter,res_iter))
       
       print "(a30, G12.3)",'Anderson relax. solution residual in L2 norm: ', and_res
       print "(a30, G12.3)",'Anderson relax. iteration defect in L2 norm: ', anditer_err

        if (exp_matlb) then
          val(1,1) = and_i + max_fp_i
          if ( res_type .eq.'SOLUTION_RESIDUAL') val(1,2) = and_res
          if ( res_type .eq.'ITERATION_RESIDUAL') val(1,2) = anditer_err
          
          call export2matlab('RESIDUAL',val)
        end if
   
       if ((and_res.lt.ande_tol).and.(anditer_err.lt.ande_tol)) then
          y_n(:) = x_n(:,sh_idx)
          deallocate(x_n)
          return
       end if
   
       and_i = and_i + 1

  end do
  
   y_n(:) = x_n(:,sh_idx)

! Release Memory:
   deallocate(diff_r,r,x_n,res_iter,x_n_ptr,val)
   
    
 return
end subroutine solve_andersone

!----------------------------------------------------------------------------------------------!
!--------------------------------------------subroutine Anderson relaxation solver------------------------------!
!----------------------------------------------------------------------------------------------!

subroutine solve_andersonr(y_n)
      use feminterface,        only: print_error
      use feminterface,        only: getsetting
      use feminterface3D,      only: solve_linear3D, export2matlab, solve_fixedpt
      use globalvariables3D,   only: numdof, x
      use femsettings,         only: readstd
      use femtypes
      implicit none
      complex (DPC),      pointer :: y_n(:)
      intent (inout) :: y_n

!-------------------------------------------------------------------------------
! local variables:
!
      integer (I4B)               :: max_fp_i, i , xi_1, yi_1, xi_2, yi_2, and_i, max_and_i
      integer (I4B)               :: shd, sh_idx
      real    (DP),   allocatable :: val(:,:)
      real     (DP)               :: epsgl,resgl , eta, andr_tol, and_res, anditer_err , fp_omega
      real     (DP),  allocatable :: diff_r(:), diff_y(:)
      complex (DPC),      pointer :: res_iter(:)
      complex (DPC),      pointer :: x_n(:,:)
      complex (DPC),      pointer :: x_n_ptr(:)
      character (len=20)          :: res_type, export_matlab
      logical                     :: exp_matlb

      exp_matlb = .false.
!__
!
   call getsetting('MAX_FIXPOINT_ITER',max_fp_i)
   call getsetting('MAX_ANDERSON_ITER',max_and_i)
   call getsetting('NONLIN_TOLERANCE',andr_tol)
   call getsetting('FPRELAX_PARAMETER',fp_omega)
   call getsetting('PLOT_RESIDUALTYPE',res_type)
   call getsetting('EXPORT2MATLAB',export_matlab)
   if (export_matlab.eq.'YES') exp_matlb = .true.
   if (exp_matlb) call export2matlab('HEADERS')
   if (exp_matlb) call export2matlab('NODES')
     shd = 4
   
     allocate(diff_r(numdof), res_iter(numdof), diff_y(numdof))
     allocate(x_n(numdof,shd))
     allocate(val(1,2))

     
     if (max_fp_i.lt.2) then  
       call print_error(3,'No. of fixed pt. iterations must be atleast 2')
     end if

   
     x_n(:,1) = y_n(:)
     
     call solve_fixedpt(y_n)
     
     max_fp_i = 4
     
     do i = 1, max_fp_i
        call solve_linear3D(epsgl,resgl,y_n,.false.,res_iter) !this will indeed perform fixed pt. step with input y_n
        y_n(:) = fp_omega*x + (1-fp_omega)*y_n(:)
        if (i.eq.max_fp_i) then
        x_n(:,3) = y_n(:)
        else if ((i.lt.max_fp_i).and.(mod(i,2).eq.1)) then
        x_n(:,2) = y_n(:)
        else if ((i.lt.max_fp_i).and.(mod(i,2).eq.0)) then 
        x_n(:,1) = y_n(:)
        end if
     end do
   
      
      if (mod(max_fp_i,2).eq.0) then
      x_n(:,3:4) = x_n (:,2:3)
      else if (mod(max_fp_i,2).eq.1) then
      x_n(:,4) = x_n(:,3)
      x_n(:,3) = x_n(:,1)
      x_n(:,1) = x_n(:,2)
      x_n(:,2) = x_n(:,3)     
      end if 
      
      print "(A18,I7,A10)" ,'Initial FP steps: ', max_fp_i 
      
     sh_idx = 4 
     and_i = 1
     
     call getsetting('MAX_FIXPOINT_ITER',max_fp_i)
     
     max_fp_i = max_fp_i + 4

     do while (and_i.le.max_and_i)
     
      readstd   = .false.
      call getsetting('FPRELAX_PARAMETER',fp_omega)
     
     
      print "(A18,I7)"    ,'   Anderson Relaxation step no.: ', and_i 
      print "(A18,G12.3)" ,'           Relaxation Parameter: ', fp_omega 
           
       xi_1 = mod(sh_idx+(shd-4),shd)+1
       yi_1 = mod(sh_idx+(shd-3),shd)+1
       xi_2 = mod(sh_idx+(shd-2),shd)+1
       yi_2 = sh_idx

       diff_r(:) = x_n(:,yi_2) - x_n(:,xi_2) - x_n(:,yi_1) + x_n(:,xi_1) !r_n - r_n-1, difference of last residuals
       diff_y(:) = x_n(:,yi_2) - x_n(:,yi_1)

!__
! 2) Compute Eta - optimal relaxation parameter:
       eta = -dot_product(diff_y(:),diff_r(:))/dot_product(diff_r(:),diff_r(:))
!__
! 3) Update solution
       sh_idx = mod(sh_idx,shd) + 1
       x_n(:,sh_idx) = x_n(:,yi_2) + eta*(x_n(:,yi_2) - x_n(:,xi_2))
   
       x_n_ptr=>x_n(:,sh_idx)
!__
! 4) Use fixed point iteration to obtain next solution:
       call solve_linear3D(epsgl,resgl,x_n_ptr,.false.,res_iter)
!__
! 5) Compute iteration error:
       y_n(:) = x_n(:,sh_idx)
       
       sh_idx = mod(sh_idx,shd) + 1
       x_n(:,sh_idx) = fp_omega*x + (1-fp_omega)*y_n(:)
       
       xi_1 = mod(sh_idx+(shd-4),shd)+1
       yi_1 = mod(sh_idx+(shd-3),shd)+1
       xi_2 = mod(sh_idx+(shd-2),shd)+1
       yi_2 = sh_idx
       
       anditer_err    = sqrt(dot_product(x_n(:,yi_2) - x_n(:,xi_2),x_n(:,yi_2) - x_n(:,xi_2)))
       and_res        = sqrt(dot_product(res_iter,res_iter))
       
       print "(a30, G12.3)",'Anderson relax. solution residual in L2 norm: ', and_res
       print "(a30, G12.3)",'Anderson relax. iteration defect in L2 norm: ', anditer_err

        if (exp_matlb) then
          val(1,1) = and_i + max_fp_i
          if ( res_type .eq.'SOLUTION_RESIDUAL') val(1,2) = and_res
          if ( res_type .eq.'ITERATION_RESIDUAL') val(1,2) = anditer_err
          
          call export2matlab('RESIDUAL',val)
        end if
   
       if ((and_res.lt.andr_tol).and.(anditer_err.lt.andr_tol)) then
          y_n(:) = x_n(:,sh_idx)
          deallocate(x_n)
          return
       end if
   
       and_i = and_i + 1
       
       if (and_i.eq.max_and_i) then
       print *, 'Max. no. of nonlinear iterations reached, solution did not converge'
       endif
       
  end do
  
   y_n(:) = x_n(:,sh_idx)

! Release Memory:
   deallocate(diff_r,diff_y,x_n,res_iter,x_n_ptr,val)
   
  return
end subroutine solve_andersonr

!----------------------------------------------------------------------------------------------!
!--------------------------------------------subroutine BB solver------------------------------!
!----------------------------------------------------------------------------------------------!

 subroutine solve_barbor(x_n)
      use feminterface3D,      only: solve_linear3D, export2matlab
      use feminterface,        only: getsetting
      use globalvariables3D,   only: numdof, x
      use femtypes
      implicit none
      complex (DPC),      pointer :: x_n(:)
      intent(inout)               :: x_n

!
!------------------------------------
! local variables:
!
      integer (I4B)               :: bb_i, max_bb_i
      real    (DP),   allocatable :: val(:,:)
      real     (DP)               :: epsgl, resgl, bbiter_err, bb_res, alpha, djac , old_bb_res, bb_tol
      complex (DPC),      pointer :: res_iter(:)
      complex (DPC),  allocatable :: delx(:), psi(:,:)
      character (len=20)          :: res_type, export_matlab
      logical                     :: exp_matlb

      exp_matlb = .false.
   
   call getsetting('NONLINITERSTEPS',max_bb_i)
   call getsetting('NONLIN_TOLERANCE',bb_tol)
   call getsetting('PLOT_RESIDUALTYPE',res_type)
   call getsetting('EXPORT2MATLAB',export_matlab)
   if (export_matlab.eq.'YES') exp_matlb = .true.
   if (exp_matlb) call export2matlab('HEADERS')
   if (exp_matlb) call export2matlab('NODES')
!__
! Allocate Memory:
   allocate(psi(numdof,2))
   allocate(delx(numdof), res_iter(numdof))
   allocate(val(1,2))
!__
! Start:
! 1) Solve the linear system using initial guess x_0 to obtain x_n+1, initial guess taken from solve_nonlinear3D
   
   psi = 0._DP  
   djac = 1._DP
   bb_res = 1.e30
   
   
   call solve_linear3D(epsgl,resgl,x_n,.false.,res_iter)
   
      
     bb_i = 1
     alpha = 1._DP
     psi(:,1) = x_n(:) - x   
   
   do while (bb_i.le.max_bb_i)
   
      print "(A18,I7,A10)" ,'BB step no.: ', bb_i 
           
     delx(:) =  - psi(:,1) / djac
      
     
     x_n(:) = x_n(:) + alpha*delx(:)
     
     call solve_linear3D(epsgl,resgl,x_n,.false.,res_iter) 
     
     psi (:,2) = x_n(:) - x   
         
     old_bb_res = bb_res
     bb_res = sqrt(dot_product(res_iter,res_iter))
     bbiter_err = sqrt(dot_product(psi(:,2),psi(:,2)))    
     print "(a30, G12.3)",'BB L2 Iteration Residual: ', bb_res
     print "(a30, G12.3)",'BB L2 Iteration Defect: ', bbiter_err

      if (exp_matlb) then
        val(1,1) = bb_i
        if ( res_type .eq.'SOLUTION_RESIDUAL') val(1,2) = bb_res
        if ( res_type .eq.'ITERATION_RESIDUAL') val(1,2) = bbiter_err
          
        call export2matlab('RESIDUAL',val)
      end if

     if (bb_res .le. old_bb_res) then
       alpha = min(2.,alpha*1.15)
     else
       alpha = max(1.e-3,alpha/2.)
     end if

     
! 4) Check for convergence, break if reached
     if ((bb_res.lt.bb_tol).and.(bbiter_err.lt.bb_tol)) then
       deallocate(psi,delx,res_iter)
       return
     end if
     
!     if (mod(bb_i,2).eq.1) then
!     
     djac = dot_product(delx,delx)/dot_product(psi(:,2)-psi(:,1),delx)
!         else 
!        djac = dot_product(delx,psi(:,2)-psi(:,1))/dot_product(psi(:,2)-psi(:,1),psi(:,2)-psi(:,1))
!     end if
          
     psi (:,1) = psi (:,2) 
     
     bb_i = bb_i + 1
     
     if (bb_i.eq.max_bb_i) then
       print *, 'Max. no. of nonlinear iterations reached, solution did not converge'
     endif
     
   end do

   deallocate(psi, delx,res_iter, val)

 return
end subroutine solve_barbor

!----------------------------------------------------------------------------------------------!
!--------------------------------------------subroutine variant of BB solver------------------------------!
!----------------------------------------------------------------------------------------------!

subroutine solve_barborv(x_n)
      use feminterface3D,      only: solve_linear3D, export2matlab
      use feminterface,        only: getsetting
      use globalvariables3D,   only: numdof, x
      use femtypes
      implicit none
      complex (DPC),      pointer :: x_n(:)
      intent(inout)               :: x_n


!
!-------------------------------------------------------
! local variables:
!
      integer (I4B)               :: bb_i, max_bb_i, i
      real    (DP),   allocatable :: val(:,:)
      real     (DP)               :: epsgl, resgl, bbiter_err, bb_res, alpha , old_bb_res, bb_tol
      complex (DPC),      pointer :: res_iter(:)
      complex (DPC),  allocatable :: delx(:), psi(:,:), djac(:)
      character (len=20)          :: res_type, export_matlab
      logical                     :: exp_matlb

      exp_matlb = .false.
!__
! Definitions:
   
   
   call getsetting('NONLINITERSTEPS',max_bb_i)
   call getsetting('NONLIN_TOLERANCE',bb_tol)
   call getsetting('PLOT_RESIDUALTYPE',res_type)
   call getsetting('EXPORT2MATLAB',export_matlab)
   if (export_matlab.eq.'YES') exp_matlb = .true.
   if (exp_matlb) call export2matlab('HEADERS')
   if (exp_matlb) call export2matlab('NODES')
! Allocate Memory:
   allocate(psi(numdof,2))
   allocate(delx(numdof), res_iter(numdof),djac(numdof))
   allocate(val(1,2))
   
   psi = 0._DP  
   djac = 1._DP
   bb_res = 1.e30
      
   call solve_linear3D(epsgl,resgl,x_n,.false.,res_iter)
   
     bb_i = 1
     alpha = 1._DP
     psi(:,1) = x_n(:) - x   
   
   do while (bb_i.le.max_bb_i)
   
      print "(A18,I7,A10)" ,'BB (variant)step no.: ', bb_i 
           
     delx(:) =  - psi(:,1) / djac(:)
          
     x_n(:) = x_n(:) + alpha*delx(:)
     
     call solve_linear3D(epsgl,resgl,x_n,.false.,res_iter) 
     
     psi (:,2) = x_n(:) - x   
         
     old_bb_res = bb_res
     bb_res = sqrt(dot_product(res_iter,res_iter))
     bbiter_err = sqrt(dot_product(psi(:,2),psi(:,2)))    
     print "(a30, G12.3)",'BB L2 Iteration Residual: ', bb_res
     print "(a30, G12.3)",'BB L2 Iteration Defect: ', bbiter_err  

      if (exp_matlb) then
        val(1,1) = bb_i
        if ( res_type .eq.'SOLUTION_RESIDUAL') val(1,2) = bb_res
        if ( res_type .eq.'ITERATION_RESIDUAL') val(1,2) = bbiter_err
          
        call export2matlab('RESIDUAL',val)
      end if

     if (bb_res .le. old_bb_res) then
       alpha = min(2.,alpha*1.15)
     else
       alpha = max(1.e-3,alpha/2.)
     end if

! 4) Check for convergence, break if reached
     if ((bb_res.lt.bb_tol).and.(bbiter_err.lt.bb_tol)) then
       deallocate(psi,delx,res_iter,djac)
       return
     end if
               
    do i = 1,numdof
      if (abs(delx(i)) .ge. 1.e3*abs(psi(i,2)-psi(i,1))) then
        djac(i) = 1.e3
      else if (abs(delx(i)) .le. 1.e-3*abs(psi(i,2)-psi(i,1))) then
        djac(i) = 1.e-3
      else
        djac(i) = delx(i)/(psi(i,2)-psi(i,1))
      end if
    end do
     
    psi (:,1) = psi (:,2) 
     
    bb_i = bb_i + 1
    
    if (bb_i.eq.max_bb_i) then
      print *, 'Max. no. of nonlinear iterations reached, solution did not converge'
    endif
     
   end do

   deallocate(psi,delx,res_iter,djac, val)

 return
end subroutine solve_barborv

!----------------------------------------------------------------------------------------------!
!--------------------------------------------subroutine Diagonal Jac. Approx. solver------------------------------!
!----------------------------------------------------------------------------------------------!

 subroutine solve_diag_jacobian(x_n)
      use feminterface3D,      only: solve_linear3D, export2matlab
      use feminterface,        only: getsetting
      use globalvariables3D,   only: numdof, x
      use femtypes
      implicit none
      complex (DPC),      pointer :: x_n(:)
      intent(inout)               :: x_n

         
!
!-----------------------------------------------
! local variables:
!
      integer (I4B)               :: dj_i, max_dj_i, i
      real    (DP),   allocatable :: val(:,:)
      real     (DP)               :: epsgl, resgl, djiter_err, dj_res, alpha , dj_tol
      complex (DPC),      pointer :: djac(:), psi(:), delx(:)
      character (len=20)          :: res_type, export_matlab
      logical                     :: exp_matlb

      exp_matlb = .false.
!__
! Definitions:
   
   call getsetting('NONLINITERSTEPS',max_dj_i)
   call getsetting('NONLIN_TOLERANCE',dj_tol)
   call getsetting('PLOT_RESIDUALTYPE',res_type)
   call getsetting('EXPORT2MATLAB',export_matlab)
   if (export_matlab.eq.'YES') exp_matlb = .true.
   if (exp_matlb) call export2matlab('HEADERS')
   if (exp_matlb) call export2matlab('NODES')
! Allocate Memory:
   allocate(djac(numdof), psi(numdof), delx(numdof))
   allocate(val(1,2))
!__
! Start:
! 1) Solve the linear system using initial guess x_0 to obtain x_n+1, initial guess taken from solve_nonlinear3D
  
    
   djac = 1._DP 
   
      call solve_linear3D(epsgl,resgl,x_n,.true.,psi)
   
      
     dj_i = 1
     alpha = 1._DP
   
   do while (dj_i.le.max_dj_i)
   
      print "(A18,I7,A10)" ,'Diag. Jacobian approx. step no.: ', dj_i 
           
     do i=1,numdof
       if (abs(djac(i)) .ge. 100*tiny(1._SP)) then
         delx(i) =  - psi(i) / djac(i)
       else
         delx(i) = 0._DP
       end if
     end do
              
    
     x_n(:) = x_n(:) + alpha*delx(:)
     
     call solve_linear3D(epsgl,resgl,x_n,.true.,psi) 
         
     dj_res = sqrt(dot_product(psi,psi))
     djiter_err = sqrt(dot_product(delx,delx))  
     
     print "(a30, G12.3)",'DJ L2 Iteration Residual: ', dj_res
     print "(a30, G12.3)",'DJ L2 Iteration Defect: ', djiter_err  
      if (exp_matlb) then
        val(1,1) = dj_i
        if ( res_type .eq.'SOLUTION_RESIDUAL') val(1,2) = dj_res
        if ( res_type .eq.'ITERATION_RESIDUAL') val(1,2) = djiter_err
          
        call export2matlab('RESIDUAL',val)
      end if
     
! 4) Check for convergence, break if reached
     if ((dj_res.lt.dj_tol)) then !.and.(djiter_err.lt.dj_tol)) then
       deallocate(djac, psi, delx)
       return
     end if
     
     djac = djac + (dot_product(delx,psi)/(djiter_err**4))*(delx*delx)
     !djac = djac + (dot_product(delx,psi)/(dj_res**4))*(psi*psi)
        
     dj_i = dj_i + 1
     
     
     if (dj_i.eq.max_dj_i) then
       print *, 'Max. no. of nonlinear iterations reached, solution did not converge'
     endif
     
   end do

   deallocate(djac, psi, delx, val)
!

!___
 return
end subroutine solve_diag_jacobian

!----------------------------------------------------------------------------------------------!
!--------------------------------------------subroutine Diagonal Jac. Approx. using iteration defect in update-------------------!
!----------------------------------------------------------------------------------------------!

 subroutine solve_diag_jacobianid(x_n)
      use feminterface3D,      only: solve_linear3D, export2matlab
      use feminterface,        only: getsetting
      use globalvariables3D,   only: numdof, x
      use femtypes
      implicit none
      complex (DPC),      pointer :: x_n(:)
      intent(inout)               :: x_n

         
!
!--------------------------------------
! local variables:
!
      integer (I4B)               :: dj_i, max_dj_i, i
      real    (DP),   allocatable :: val(:,:)
      real     (DP)               :: epsgl, resgl, djiter_err, dj_res, alpha, dj_tol 
      complex (DPC),      pointer :: djac(:), psi(:), delx(:), res_iter(:)
      character (len=20)          :: res_type, export_matlab
      logical                     :: exp_matlb

      exp_matlb = .false.

   
   call getsetting('NONLINITERSTEPS',max_dj_i)
   call getsetting('NONLIN_TOLERANCE',dj_tol)
   call getsetting('PLOT_RESIDUALTYPE',res_type)
   call getsetting('EXPORT2MATLAB',export_matlab)
   if (export_matlab.eq.'YES') exp_matlb = .true.

! Allocate Memory:
   allocate(djac(numdof), psi(numdof), delx(numdof), res_iter(numdof))
   allocate(val(1,2))
   
   djac = 1._DP 
  
   call solve_linear3D(epsgl,resgl,x_n,.false.,res_iter)
        
     dj_i = 1
     alpha = 1._DP
     psi = x_n(:) - x

     if (exp_matlb) call export2matlab('HEADERS')
     if (exp_matlb) call export2matlab('NODES')
   
   do while (dj_i.le.max_dj_i)
     
      print "(A18,I7,A10)" ,'Diag. Jacobian(id) approx. step no.: ', dj_i 
          
     alpha = min(1.,alpha*1.2)
     
     do i=1,numdof
       if (abs(djac(i)) .ge. 100*tiny(1._SP)) then
         delx(i) =  - psi(i) / djac(i)
       else
        ! delx(i) = 0._DP
       end if
     end do
         
     x_n(:) = x_n(:) + alpha*delx(1:numdof)
     
     call solve_linear3D(epsgl,resgl,x_n,.false.,res_iter) 
     
     psi = x_n(:) - x   
         
     !dj_res = sqrt(dot_product(res_iter,res_iter))
     !djiter_err = sqrt(dot_product(psi,psi)) 
     dj_res = sqrt(dot_product(res_iter,res_iter))
     djiter_err = sqrt(dot_product(delx,delx)) 
     print "(a30, G12.3)",'DJID L2 Iteration Residual: ', dj_res
     print "(a30, G12.3)",'DJID L2 Iteration Defect: ', djiter_err
      if (exp_matlb) then
        val(1,1) = dj_i
        if ( res_type .eq.'SOLUTION_RESIDUAL') val(1,2) = dj_res
        if ( res_type .eq.'ITERATION_RESIDUAL') val(1,2) = djiter_err
          
        call export2matlab('RESIDUAL',val)
      end if
       
! 4) Check for convergence, break if reached
     if ((dj_res.lt.dj_tol)) then !.and.(djiter_err.lt.dj_tol)) then
       deallocate(djac, psi, delx,res_iter)
       return
     end if
     
     !djac = djac + (dot_product(delx,psi)/(djiter_err**4))*(delx*delx)
     djac = djac + (dot_product(delx,psi)/sum(delx**4))*(delx*delx)
             
     dj_i = dj_i + 1
     
     if (dj_i.eq.max_dj_i) then
       print *, 'Max. no. of nonlinear iterations reached, solution did not converge'
     endif
     
   end do

   deallocate(djac, psi, delx, res_iter, val)

 return
end subroutine solve_diag_jacobianid

!----------------------------------------------------------------------------------------------!
!--------------------------------------------subroutine Fixed point solver------------------------------!
!----------------------------------------------------------------------------------------------!

subroutine solve_fixedpt(x_n)
      use feminterface3D,      only: solve_linear3D, export2matlab, getbcval
      use feminterface,        only: getsetting
      use globalvariables3D,   only: numdof, x, nnat, vv, numv, vgdof, sbc, bctype, pvalue, nod, vn
      use femsettings,         only: readstd
      use femtypes
      implicit none
      complex (DPC),      pointer :: x_n(:)
      intent (inout)              :: x_n

!----------------------------------------------------
! local variables:

      real    (DP)                :: epsgl, resgl, fp_tol, fp_omega, fpiter_l1err, fp_l1res, fp_res , fpiter_err, xyzs(3)
      real    (DP),   allocatable :: val(:,:)
      integer (I4B)               :: max_fp_i, fp_i, elem, f, inat, surf, bc_index, n
      integer (I4B),    parameter :: f2n(3,4)=reshape((/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/))
      complex (DPC),  pointer     :: res_iter(:)
      complex (DPC)               :: pval
      character (len=20)          :: res_type, export_matlab
      logical                     :: exp_matlb

      exp_matlb = .false.
      
   
   call getsetting('MAX_FIXPOINT_ITER',max_fp_i)
   call getsetting('NONLIN_TOLERANCE',fp_tol)
   readstd   = .false.
   call getsetting('FPRELAX_PARAMETER',fp_omega)
   readstd   = .true.
   call getsetting('PLOT_RESIDUALTYPE',res_type)
   call getsetting('EXPORT2MATLAB',export_matlab)
   
   if (export_matlab.eq.'YES') exp_matlb = .true.
! Allocate Memory:
   allocate(res_iter(numdof))
   allocate(val(1,2))
!__
      fp_i = 1
      if (exp_matlb) call export2matlab('HEADERS')
      if (exp_matlb) call export2matlab('NODES')
! Start:
! 1) Solve the linear system using initial guess x_0 to obtain x_n+1
   
       do while (fp_i.le.max_fp_i)
          readstd   = .false.
          call getsetting('FPRELAX_PARAMETER',fp_omega)
          readstd   = .true.
         
          print "(A1)"         ,' '
          print "(A79)"        ,' =============================================================================== '
          print "(A46,I)"      ,'|                    FIXED POINT ITERATION No.', fp_i
          print "(A46,G12.3)"  ,'|                    FIXED POINT  RELAXATION: ', fp_omega
          print "(A79)"        ,' ------------------------------------------------------------------------------- '

         call solve_linear3D(epsgl,resgl,x_n,.false.,res_iter)
!__
! 2) Compute the difference between current and former solutions (fp_res , fpiter_err)
         fp_res     = norm2(real(res_iter))
         fp_l1res   = sum(abs(res_iter))
         fpiter_err = sqrt(dot_product(x(:)-x_n(:),x(:)-x_n(:)))
         fpiter_l1err = sum(abs(x(:)-x_n(:)))
         
         
          print "(a30, G12.3)",'FP L1 Iteration Residual: ', fp_l1res
          print "(a30, G12.3)",'FP L2 Iteration Residual: ', fp_res
          print "(a30, G12.3)",'FP L1 Iteration Defect: ', fpiter_l1err
          print "(a30, G12.3)",'FP L2 Iteration Defect: ', fpiter_err
          
          if (exp_matlb) then
            val(1,1) = fp_i
            if ( res_type .eq.'SOLUTION_RESIDUAL') val(1,2) = fp_res
            if ( res_type .eq.'ITERATION_RESIDUAL') val(1,2) = fpiter_err
          
            call export2matlab('RESIDUAL',val)
            !call export2matlab('SOLUTION')
          end if
!__
! 3) Store the current solution vector
         x_n(:) = fp_omega*x + (1-fp_omega)*x_n(:)


         do elem = 1, numv
           do f = 1, 4
             if (vv(f,elem).lt.0) then
               surf = abs(vv(f,elem))
               bc_index = sbc(surf)
               do inat = 1, nnat
                 if (bctype(bc_index,inat).lt.200) then
                   ! TODO: also consider dof on edges and on faces
                   do n = 1, 3 ! Apply bc to node-dof
                     pval = 0._DPC
                     xyzs(1:3) = nod(1:3,vn(f2n(n,f),elem))
                     call getbcval( sbc( abs(vv(f,elem))), inat, xyzs, pval ) 
                     x_n(vgdof(elem,inat)%d(f2n(n,f))) = pval
                   end do
                 end if
               end do
             end if 
           end do 
         end do

         if ((fp_res.lt.fp_tol).and.(fpiter_err.lt.fp_tol)) then
           deallocate(res_iter)
           return
         end if
        
         fp_i = fp_i + 1
         
         if (fp_i.eq.max_fp_i) then
           print *, 'Max. no. of nonlinear iterations reached, solution did not converge'
         endif
       
       end do 
      
    deallocate(res_iter, val)
 return
  end subroutine solve_fixedpt
  
  
subroutine get_starting_solution (x_0)
      use femtypes
      use feminterface,      only: fetchmatparameters
      use feminterface3D,    only: findelement3D
      use globalvariables3D, only: Elch, Kboltz, eps0, mu0, omega, pi, physics, pvalue, nod, vn, vv, numn, numv, dommat, dom, zmax, zmin
      use matconstants
      implicit none

      complex (DPC), pointer    :: x_0(:)
      intent (inout)            :: x_0
!
!-------------------------------------------------------------------------------
! This subroutine generates an initial guess x_0 to help the nonlinear
! system solver get started with a close initial guess.
!
!-------------------------------------------------------------------------------
! local variables:
!
      real (DP)                 :: xyzs(3), T_ref, N_D, N_A, n_i, randnum(2*numn)
      integer (i4B)             :: n, matindex, elem
      real(DP),     allocatable :: list(:)
      logical                   :: found
!__
! Definitions
! So far only initial guess for nature 1 for the semiconductor equations:

!__
! Select Physics:
    select case (physics)
      case('3DSEMICONDUCTOR', '3DMECHANICSSEMICOND')
        do n = 1, numn
            xyzs(:) = nod(:,n)
            allocate(list(numparam))
            call findelement3D(xyzs,nod,vn,vv,numv,elem,found)
            call fetchmatparameters(list,dommat(dom(elem)),xyzs,elem)
          
            N_A   = list(4)
            N_D   = list(5)
            T_ref = list(6)
            n_i   = list(9)
            
            x_0(n)        = ((Kboltz*T_ref)/Elch)*asinh((N_D-N_A)/(2*n_i))

            !x_0(n)        = (-0.15/50E-6)*xyzs(2)-0.15 + ((Kboltz*T_ref)/Elch)*asinh((N_D-N_A)/(2*n_i))
            !x_0(numn+n)   = ((-0.15/50E-6)*xyzs(2)-0.15)*0.3
            !x_0(2*numn+n) = (0.3/20E-6)*xyzs(3)-0.3
            deallocate(list)
          
        end do
        
        case default
        
      end select
end subroutine
