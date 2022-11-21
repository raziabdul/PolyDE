subroutine solve_adapt3D()
      use globalvariables3D,  only: vp, nnat, numv, numf, vf, vp, vn, nod, numn, numdof
      use feminterface,       only: getsetting, destroyarrptr, print_error
      use feminterface3d,     only: writeunv, writeng, writedata, preassemb3D
      use feminterface3d,     only: solve_linear3D, solve_nonlinear3D, vol2aux, bcstovv
      use feminterface3d,     only: h_adapt3D, p_adapt3D, hp_adapt3D
      use feminterface3d,     only: prepare_h_adapt3D, prepare_p_adapt3D, prepare_hp_adapt3D, meshquality3d
      use feminterface3d,     only: estimate_error3D, vtk_varmatplot
      use feminterface3d,     only: getkpv, writeoutnodevalues, writeoutelementvalues, get_modres
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
!    $Date: 2015/11/11 17:02:36 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'solve_adapt3D'
!    This is a core managing routine of the solver for the case of solution by adaptation.
!    It estimates the accuracy of a solution and performs
!    refinements of the computational domain according to settings provided by the user.
!    Three refinement mechanisms are implied right now, that is:
!    1) h-Adaptation
!    2) p-Adaptation
!    3) hp-Adaptation
!
! Error Estimation formulas (to be removed):
!    error  = 100._DP*sqrt(sum(res)/sumref)
!    maxres = (adapterror/(100._DP*numv))**2 * sumref
!    maxres(inat) = (adapterror/(100._DP*numv))**2 * sumref(inat)
!
!    Adaptation Loop calls Adaptation routines.
!    Adaptation Routines (n times):
!       - Solve the problem
!       - Estimate the error
!       - Mark elements / DOFs
!       - Apply Refinement Algorithms
!       - Set Up adapted DOFs h- and p- Grids 
!
!------------------------------------------------------------------------------
!  internal variables
      integer (I4B)              :: i, j, astep, adaptsteps
      integer (I4B)              :: num_prints
      real     (DP)              :: adapterror, epsgl, resgl, vrnorm_max
      real     (DP), allocatable :: plotdata_xy(:,:,:)
!      real     (DP), allocatable :: mod_res(:,:),keyPointValues(:,:), real_vp(:)
      real     (DP),     pointer :: e_len(:)=>null()
      real     (DP), allocatable :: v_res(:,:), v_res_norm(:)
!      real     (DP),     pointer :: res_old(:,:)=>null()
!      real     (DP), allocatable :: sing_val(:,:)=> null()
      type(ARRPTRI),     pointer :: ev(:)=>null(), es(:)=>null(), fv(:)=>null()
      character (len=16)         :: adapt_type, mqplot, filenameFormat, fileName, nonlinear, matvar
      logical                    :: ok
    
!__
! Definitions:
!      allocate (vp(numv,nnat))
      call getsetting('ADAPT_STEPS',adaptsteps)
      call getsetting('ADAPTION_ERROR',adapterror)
      call getsetting('MESHQUALITYPLOT',mqplot)
      call getsetting('MATER_VARIATION',matvar)
      call getsetting('NONLINEAR',nonlinear)
!-
      call getsetting('ADAPTION_TYPE',adapt_type)
      adapt_type = adapt_type(1:len_trim(adapt_type))
!-
      mqplot = mqplot(1:len_trim(mqplot))
      allocate (plotdata_xy(adaptsteps+1,nnat,2))

!      allocate(res_old(numv,nnat))
!____
! Start:
!__
! 1) Get Auxilarity:
      call vol2aux
!      call vol2aux( numf, vf)
!
      call bcstovv
!-
      if (mqplot.eq.'YES') call meshquality3d(vn,numv,nod)
      if (matvar.eq.'YES') then
        call vtk_varmatplot
        print*,' '
        print*,'*** User-defined material is used.'
        print*,'*** Please check the material properties plot in the file "USERMAT_PROPERTY.vtk" before proceeding.'
        print*,' '
        pause
      end if
!__
! 2) Solve FEM to get a Solution
      call preassemb3D

      if (nonlinear.eq.'YES') then
        call solve_nonlinear3D()
      else
        call solve_linear3D(epsgl,resgl)
      end if
!___________________________________________S T A R T__O F__A D A P T A T I O N____
      astep = 1
      do while (astep.le.adaptsteps)
!__
! 3) Estimate the Error:
        allocate(v_res(numv,nnat),v_res_norm(nnat))
        call estimate_error3D(v_res,v_res_norm)
!        if (astep.eq.1) res_old = v_res
!        if (astep.eq.2) then
!          call writedata(0,'Residuals_N1.dat','XYPLOT','Convergence Data','NDOF','Rel. Error', res_old(:,1) ,v_res(:,1))
!          call writedata(0,'Residuals_N2.dat','XYPLOT','Convergence Data','NDOF','Rel. Error', res_old(:,2) ,v_res(:,2))
!          call writedata(0,'Residuals_N3.dat','XYPLOT','Convergence Data','NDOF','Rel. Error', res_old(:,3) ,v_res(:,3))
!          sum_err1 = sum(res_old)
!          sum_err2 = sum(v_res)
!        end if
        vrnorm_max = maxval(v_res_norm(:))
        do i = 1, nnat
          plotdata_xy(astep,i,1) = numdof
          plotdata_xy(astep,i,2) = v_res_norm(i)
        end do
!- IF Error estimate of the worse element is lower that error bound, exit loop
        if (vrnorm_max.lt.adapterror) then
          deallocate(v_res,v_res_norm)
          print "(a2)" ,'--'
          print "(a31)",'STOP: Reached requires accuracy.'
          exit
        end if

        select case (adapt_type)
          case ('NO_ADAPT')

          case ('H_ADAPT')
            call prepare_h_adapt3D(ev,es,e_len)
            call h_adapt3D(astep,v_res,ev,es,e_len)

          case ('P_ADAPT')
            call prepare_p_adapt3D
            call p_adapt3D(astep,v_res)

          case ('HP_ADAPT')
            call prepare_hp_adapt3D(ev,es,e_len)
            call hp_adapt3D(astep,v_res,ev,es,e_len)

          case default
            call print_error(5,'No such adaptation algorithm:',adapt_type)

        end select

        deallocate(e_len)
        ok=destroyarrptr(ev)
        ok=destroyarrptr(es)
        deallocate(v_res,v_res_norm)
!-
        if ((mqplot.eq.'YES').and.(adapt_type.ne.'P_ADAPT')) call meshquality3d(vn,numv,nod)
!-
! 4) Solve FEM to get final Solution
        call preassemb3D
        if (nonlinear.eq.'YES') then
          call solve_nonlinear3D()
        else
          call solve_linear3D(epsgl,resgl)
        end if
        astep = astep + 1
      end do
!_______________________________________________E N D__O F__A D A P T A T I O N____
!__
! 5) Estimate the final error:
      allocate(v_res(numv,nnat),v_res_norm(nnat))
      call estimate_error3D(v_res,v_res_norm)
      vrnorm_max = maxval(v_res_norm(:))
      do i = 1, nnat
        plotdata_xy(adaptsteps+1,i,1) = numdof
        plotdata_xy(adaptsteps+1,i,2) = v_res_norm(i)
      end do
      print "(A79)"       ,' _______________________________________________________________________________ '
!__
! 6) Inform the user about stop without reaching the error bound:
      print "(A2)"          ,'--'
      print "(A48)"         ,'STOP: Maximum number of adaptation steps reached.'
      print "(A53,F5.2,A27)",'WARNING: Error estimate of worst element error is by: ', vrnorm_max/adapterror,' above required error bound!'
!__
! 7) Print out Error History on the screen
      print*,'Error Data:'
      print *,'============================================================================'
      print *,'   AdaptStep   |   Nature   |        NDOF         |        Error'
      print *,'----------------------------------------------------------------------------'
      num_prints = astep
      do j = 1, nnat
        do i = 1, num_prints
          print *, (i-1), j, '    ',plotdata_xy(i,j,1), plotdata_xy(i,j,2)
        end do
        if (j.ne.nnat) print *,' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
      end do
      print *,'============================================================================'
!__
! 8) Write Table to file Convergence.* in PROJECTPATH
      write(filenameFormat,'(A6,I1,A4)') '(A11,I', ceiling(log10(real(nnat))),',A4)'
      do i = 1, nnat
      write(fileName,filenameFormat) 'Convergence', i, '.dat'
        call writedata(0,fileName,'XYPLOT','Convergence Data','NDOF','Rel. Error',plotdata_xy(1:num_prints,i,1),plotdata_xy(1:num_prints,i,2))
      end do
!-
      call writeunv(ok)
      call writeng(ok)
!__
! Release Memory
      deallocate(v_res,v_res_norm)
      ok = destroyarrptr(fv)
      deallocate(plotdata_xy)
!      deallocate(res_old)
!      deallocate(keyPointValues, mod_res, sing_val, real_vp)
!
!_End.
      return
end subroutine solve_adapt3D
