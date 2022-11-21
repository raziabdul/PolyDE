      subroutine transientsolver(acsr,acsrs,acsrd,acsrm,rhs,ia,ja,d2xdt2_at_tminusdt,stepcount,returncount, &
      &          trangtol,tranltollo,tranltolup,currenttime,timestep,timestepreject,epsgl,resgl)
      use feminterface, only: zeit, getsetting, umfsolver, umfcsolver, pardiso_solver, nonlinearsolver, &
                              calctimestep_newmark, calctimestep_runge, amux
      use femtypes
      implicit none
      complex (DPC), pointer :: acsr(:), acsrs(:), acsrd(:), acsrm(:), rhs(:)
      complex (DPC), pointer :: d2xdt2_at_tminusdt(:)
      integer (I4B), pointer :: ia(:), ja(:)
      integer (I4B) :: stepcount, returncount
      real (DP) :: trangtol, tranltollo, tranltolup, currenttime, timestep
      real (DP) :: epsgl, resgl
      logical :: timestepreject
      intent (in) :: trangtol
      intent (out) :: epsgl, resgl
      intent (inout) :: stepcount, returncount, currenttime, timestep, timestepreject
      intent (inout) ::  tranltollo, tranltolup
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
!    $Revision: 1.5 $
!    $Date: 2014/02/18 12:28:27 $
!    $Author: juryanatzki $
!
!-------------------------------------------------------------------------------
!
!  local variables:
      character (len=16) :: transsolver
      real (DP) :: newmark_params(8)
! To check dummy variables
      acsr = 0.
      acsrs = 0.
      acsrd = 0.
      acsrm = 0.
      rhs = 0.
      ia = 0
      ja = 0
      d2xdt2_at_tminusdt = 0.
      stepcount = 0
      returncount = 0
      tranltollo = 0.
      tranltolup = 0.
      currenttime = 0
      timestep = 0.
      timestepreject = 0.
      epsgl = 0.
      resgl = 0.

      call getsetting('TRANSIENTSOLVER',transsolver)
      print*,'TRANSIENTSOLVER is not (no more) available'
!
!      call getsetting('NONLINEAR',nonlinear)
!      call getsetting('NONLINEARSOLVER',nonlinsolver)
!
!      select case (transsolver)
!
!        case ('NEWMARK')
!
!          call getsetting('DELTA_NEWMARK',delta_newmark)
!          call getsetting('ALPHA_NEWMARK',alpha_newmark)
!
!!  Calculate integration constants for the Newmark method
!          a0 = 1._DP/(alpha_newmark*(timestep**2))
!          a1 = delta_newmark/(alpha_newmark*timestep)
!          a2 = 1._DP/(alpha_newmark*timestep)
!          a3 = (1._DP/(2._DP*alpha_newmark))-1._DP
!          a4 = (delta_newmark/alpha_newmark) - 1._DP
!          a5 = (timestep/2._DP)*((delta_newmark/alpha_newmark)-2._DP)
!          a6 = timestep*(1._DP-delta_newmark)
!          a7 = delta_newmark*timestep
!
!          select case(nonlinsolver)
!            case ('NOJACOBIAN','MOD_NEWRAPH','ARC_LENGTH','ARC_LENGTH_LINS')
!              ! Do nothing, the effective stiffness is composed during the nonlinear iteration
!!  Case default cares for both linear and nonlinear cases without nonlinearity iteration in every time step
!            case default
!!  Calculate effective stiffness matrix
!              allocate(acsr(size(acsrs)))
!              acsr = (0._DP,0._DP)
!              acsr = acsrs + a0*acsrm + a1*acsrd
!!  Calculate effective loads for this time step
!              allocate(csrvec(size(rhs)),csrvec2(size(rhs)))
!              csrvec = 0._DPC
!              csrvec2 = 0._DPC
!              call amux(ndof,(a1*x + a4*dxdt + a5*d2xdt2),csrvec,acsrd,ja,ia)
!              call amux(ndof,(a0*x + a2*dxdt + a3*d2xdt2),csrvec2,acsrm,ja,ia)
!              csrvec = rhs + csrvec + csrvec2
!              deallocate(csrvec2)
!          end select
!!
!!  Copy x in x_at_t before we solve for the next x
!          allocate(x_at_t(size(x)))
!          x_at_t = x
!!  Solve for x
!!  Call nonlinearsolver routine which does nonlinear iteration
!          if (nonlinear .eq. 'YES') then
!            newmark_params = (/ a0,a1,a2,a3,a4,a5,a6,a7 /)
!            call nonlinearsolver(acsr,x_at_t,csrvec,ia,ja,newmark_params,epsgl,resgl)
!!  If no iteration is required do not call the nonlinear solver, solve
!!  the system once for this time step and continue
!          else
!!  Copy ia and ja in ia_local and ja_local because the UMFSOLVER deallocates these pointers
!            allocate(ia_local(size(ia)))
!            allocate(ja_local(size(ja)))
!            ia_local = ia
!            ja_local = ja
!!  Solve
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
!!  Recopy the ia and ja vectors
!            allocate(ia(size(ia_local)))
!            allocate(ja(size(ja_local)))
!            ia = ia_local
!            ja = ja_local
!            deallocate(ia_local,ja_local)
!          end if
!!
!!  Calculate d2xdt2_at_tplusdt and dxdt_at_tplusdt
!          allocate(dxdt_at_tplusdt(size(dxdt)))
!          allocate(d2xdt2_at_tplusdt(size(d2xdt2)))
!          d2xdt2_at_tplusdt = a0*(x-x_at_t) - a2*dxdt - a3*d2xdt2
!          dxdt_at_tplusdt = dxdt + a6*d2xdt2 + a7*d2xdt2_at_tplusdt
!!
!!  At this point we have the d2xdt2_at_tminusdt vector from the previous iteration,
!!  the d2xdt2 which has not yet been updated and the d2xdt2_at_tplusdt which has just been calculated
!!  Calculate the next time step. d2xdt2 does not need to be passed because it is a global variable
!          call calctimestep_newmark(d2xdt2_at_tminusdt,d2xdt2_at_tplusdt,timestep,nexttimestep, &
!          &    currenttime,alpha_newmark,tranltollo,tranltolup,timeglobalerror)
!!
!!  Copy current d2xdt2 at d2xdt2_at_tminusdt so that in the next iteration this is
!!  considered as the state, 2 time steps before and is used in the calctimestep routine
!          d2xdt2_at_tminusdt = d2xdt2
!!
!          timestep = nexttimestep
!!  Update d2xdt2 and dxdt so that they can be used in the next iteration
!          d2xdt2 = d2xdt2_at_tplusdt
!          dxdt = dxdt_at_tplusdt
!!  Deallocate tplusdt time derivatives
!          deallocate(d2xdt2_at_tplusdt,dxdt_at_tplusdt)
!          deallocate(x_at_t,csrvec)
!!  Update current time
!!  Update time step
!          currenttime = currenttime+timestep
!          call getsetting('END_TIME',endtime)
!          if (currenttime .ge. endtime) then
!            returncount = returncount+1
!            if (returncount .eq. 2) return
!            timestep = timestep - abs(currenttime-endtime)
!            currenttime = endtime
!          end if
!!  Check for the global tolerance
!!  If exceeded start computation from the beginning
!          if (timeglobalerror .gt. trangtol) then
!            call getsetting('START_TIME',currenttime)
!            call getsetting('INIT_TIME_STEP',timestep)
!            tranltollo = tranltollo * ( (trangtol/timeglobalerror)**1.5_DP )
!            tranltolup = tranltolup * ( (trangtol/timeglobalerror)**1.5_DP )
!
!            print*,"tranltollo    ",tranltollo
!            print*,"tranltolup    ",tranltolup
!
!            stepcount = 0
!            returncount = 0
!!  Reset the solution vectors
!            x = 0._DPC
!            dxdt = 0._DPC
!            d2xdt2 = 0._DPC
!          end if
!!
!          print*,"stepcount     ",stepcount
!          print*,"nexttimestep  ",nexttimestep
!          print*,"currenttime   ",currenttime
!!
!        case ('WILSON_THETA')
!
!
!        case ('HOUBOLT')
!
!
!        case ('CENTRALDIFF')
!
!
!        case ('RUNGE2NDORDER')
!
!          allocate(x_at_t(size(x)))
!          x_at_t = x
!
!          call getsetting('LINSOLVERTYPE',solver)
!
!          allocate(acsr(size(acsrs)))
!          acsr = (0._DP,0._DP)
!          acsr = acsrs + (2._DP/timestep)*acsrd
!
!          allocate(csrvec(size(rhs)))
!          csrvec = 0._DPC
!          call amux(ndof,x_at_t,csrvec,acsrd,ja,ia)
!          csrvec = rhs + (2._DP/timestep)*csrvec
!
!          allocate(ia_local(size(ia)))
!          allocate(ja_local(size(ja)))
!          ia_local = ia
!          ja_local = ja
!
!          select case(solver)
!          case ('UMF')
!            call umfsolver(acsr,csrvec,x,ndof,eps,ia,ja,epsgl,resgl)
!          case ('UMFC')
!            call umfcsolver(acsr,csrvec,x,ndof,ia,ja)
!          case default
!            print*,' no such solver: ',solver,' use UMF or UMFC'
!          end select
!
!          allocate(k2(size(x)))
!          k2 = 2*(x - x_at_t)
!
!          allocate(acsr(size(acsrs)))
!          acsr = (0._DP,0._DP)
!          acsr = acsrs + (1._DP/timestep)*acsrd
!
!          allocate(ia(size(ia_local)))
!          allocate(ja(size(ja_local)))
!          ia = ia_local
!          ja = ja_local
!
!          csrvec = 0._DPC
!          call amux(ndof,x_at_t,csrvec,acsrd,ja,ia)
!          csrvec = rhs + (1._DP/timestep)*csrvec
!          select case(solver)
!          case ('UMF')
!            call umfsolver(acsr,csrvec,x,ndof,eps,ia,ja,epsgl,resgl)
!          case ('UMFC')
!            call umfcsolver(acsr,csrvec,x,ndof,ia,ja)
!          case default
!            print*,' no such solver: ',solver,' use UMF or UMFC'
!          end select
!
!          allocate(k1(size(x)))
!          k1 = x - x_at_t
!
!          allocate(ia(size(ia_local)))
!          allocate(ja(size(ja_local)))
!          ia = ia_local
!          ja = ja_local
!
!          deallocate(ia_local,ja_local)
!
!!  Calculate error in order to determine next time step
!          call calctimestep_runge(k1,k2,timestep,nexttimestep,timestepreject)
!
!          if (timestepreject) then
!            x = x_at_t
!            currenttime = currenttime-timestep+nexttimestep
!          else
!            currenttime = currenttime + nexttimestep
!            call getsetting('END_TIME',endtime)
!            if (currenttime .ge. endtime) then
!              returncount = returncount+1
!              if (returncount .eq. 2) return
!              nexttimestep = nexttimestep - abs(currenttime-endtime)
!              currenttime = endtime
!            end if
!          end if
!          timestep = nexttimestep
!          
!          deallocate(x_at_t,csrvec)
!          
!        case default
!
!      end select
!
      return
      end subroutine transientsolver
!
!
      subroutine calctimestep_newmark(d2xdt2_at_tminusdt,d2xdt2_at_tplusdt,timestep,nexttimestep, &
      &          currenttime,alpha_newmark,tranltollo,tranltolup,timeglobalerror)
      use feminterface, only: getsetting
      use femtypes, only: DP, DPC
      use globalvariables, only: d2xdt2
      implicit none
      complex (DPC), pointer :: d2xdt2_at_tminusdt(:), d2xdt2_at_tplusdt(:)
      real (DP) :: timestep, nexttimestep, currenttime, alpha_newmark
      real (DP) :: tranltollo, tranltolup, timeglobalerror
      intent (in) :: timestep, currenttime, alpha_newmark, tranltollo, tranltolup
      intent (out) :: nexttimestep, timeglobalerror
!
!  local variables:
      real (DP) :: timesteperror = 0._DP
!
      nexttimestep = 0._DP
      timeglobalerror = 0._DP
!
!  timesteperror calculated in Euclidean-norm
      timesteperror = sqrt ( &
                    &   sum ( ( (timestep**2._DP)*((1._DP/6._DP)-alpha_newmark)*(d2xdt2-d2xdt2_at_tplusdt) )**2._DP &
                    &       ) &
                    &      )
!
      timeglobalerror = (currenttime/timestep) * timesteperror
!
!      timesteperror = sqrt ( sum ( (timestep**2._DP) * ((1._DP/6._DP)-alpha_newmark) * (d2xdt2-d2xdt2_at_tplusdt) &
!                    &            ) **2._DP &
!                    &      )
!      timesteperror =  sqrt( sum( (((timestep**2._DP)/24._DP) * &
!                    & (d2xdt2_at_tminusdt + ((2._DP-(24._DP*alpha_newmark))*d2xdt2) + (((24._DP*alpha_newmark)-3._DP)*d2xdt2_at_tplusdt)) &
!                    & )**2._DP ) )
!      timesteperror = timesteperror / (sqrt(sum(x**2._DP)))
!
      if (timesteperror .lt. tranltollo) then
        nexttimestep = timestep * ( (tranltollo/timesteperror)**(1._DP/3._DP) )
      else if (timesteperror .gt. tranltolup) then
        nexttimestep = timestep * ( (tranltolup/timesteperror)**(1._DP/3._DP) )
      else
        nexttimestep = timestep
      end if
!
      print*,"timesteperror  ",timesteperror
      print*,"timeglobalerror",timeglobalerror
!
      return
      end subroutine calctimestep_newmark
!
!
      subroutine calctimestep_runge(k1,k2,timestep,nexttimestep,timestepreject)
      use feminterface, only: getsetting
      use femtypes
      use globalvariables
      implicit none
      real (DP)   :: timestep, nexttimestep
      complex (DPC), pointer :: k1(:), k2(:)
      logical :: timestepreject
      intent(in) :: timestep
      intent(out) :: nexttimestep
      intent(inout) :: timestepreject

!  local variables:
      real (DP) :: timesteperror = 0._DP, tranreltoler = 0._DP
      
      timesteperror = maxval(abs(k2-k1)/(3._DP*abs(x)))
      call getsetting('TRAN_REL_TOLER',tranreltoler)

      if (timesteperror .le. tranreltoler) then
        nexttimestep = timestep*min( 1.5_DP,SQRT(tranreltoler/(1.2_DP*timesteperror)) )
        timestepreject = .false.
      else
        nexttimestep = timestep*max( 0.1_DP,SQRT(tranreltoler/(1.2_DP*timesteperror)) )
        timestepreject = .true.
      end if

      deallocate(k1,k2)

      return
      end subroutine calctimestep_runge
