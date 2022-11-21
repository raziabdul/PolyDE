      subroutine adaptation(epsgl)
      use feminterface, only: assemblyandsolve, assemblyandsolvetrans, basout, getsetting, &
     &                         hadapt, padapt, preassemb, hpadapt, residual
      use femtypes
      use globalvariables, only: physics, ep, nnat, n, xn, yn, p, e, en, &
                               & kzi, geb
      use fluidinterface, only:
      implicit none
      real (DP), intent(out) :: epsgl
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
!    $Revision: 1.28 $
!    $Date: 2015/04/01 10:56:13 $
!    $Author: juryanatzki $
!
!-------------------------------------------------------------------------------
!  This program will do the adaptation if chosen by setting parameters ADAPTION_
!  TYPE in FEMsettings.txt. First the kind of error estimator is chosen.
!
!  Output:
!            epsgl      reached error by linear solver
!
!  internal variables:
      integer (I4B) ::  polyorder, errcode
      real (DP) :: error, resgl, sumres(nnat), sumref(nnat)
      real (DP), allocatable :: res(:,:)
      logical :: ende, ende2, ext, int, ok
      character (len=16) :: adapttype, estimator, algorithm
!
!  select how to calculate the residual
      call getsetting('ERROR_ESTIMATOR',estimator)
      estimator = estimator(1:len_trim(estimator))
      select case (estimator)
        case ('BOUNDARYRESIDUAL')
          ext = .true.
          int = .false.
        case ('INTERIORRESIDUAL')
          ext = .false.
          int = .true.
        case ('FULLRESIDUAL')
          ext = .true.
          int = .true.
        case default
          print "(a)","No specific error estimator chosen. Using default (full residual)."
          ext = .true.
          int = .true.
      end select
!
!  do the adaptation
      call getsetting('ADAPTION_TYPE',adapttype)
      adapttype = adapttype(1:len_trim(adapttype))
      select case (adapttype)
!
!
!  reset the polyorder, assemble and solve
        case ('SET_POLYORDER')
          call getsetting('PHYSICS_MODE',physics)
!
          if (.not. associated(ep)) then
            allocate (ep(n,nnat))
          end if
!  Get and set starting value for polynomial degree from FEMsettings.txt
          call getsetting('POLYORDER',polyorder)
          ep(1:n,1:nnat)=polyorder
          select case(physics)
            case('FLUIDINCOMPR','FLUIDELMASSCON')
              select case(polyorder)
                case(1,2)
                  ep(1:n,1:2) = 2
                  ep(1:n,3)   = 1
                case default
                  ep(1:n,3)   = polyorder-1
              end select
              call preassemb
            case('STOKES')
              !Do not call preassemb
            case default
              call preassemb
          end select
!
          select case(physics)
!  Fluid no-adaption
            case('STOKES')
              call solver_stokes
              return
            case('FLUID')
              call solver_fluid
!
!  Transient without adaptation
            case('TRANHEATTR')
              call assemblyandsolvetrans(epsgl,resgl)
!
!  Default, assemble and solve for the first time.
            case default
              call assemblyandsolve(epsgl,resgl)
          end select
!
!  only assemble and solve
        case ('NO_ADAPT')
          call getsetting('PHYSICS_MODE',physics)
!
          if (.not. associated(ep)) then
            allocate (ep(n,nnat))
!  Get and set starting value for polynomial degree from FEMsettings.txt
            call getsetting('POLYORDER',polyorder)
            ep(1:n,1:nnat)=polyorder
            select case(physics)
              case('FLUIDINCOMPR','FLUIDELMASSCON')
                select case(polyorder)
                  case(1,2)
                    ep(1:n,1:2) = 2
                    ep(1:n,3)   = 1
                  case default
                    ep(1:n,3)   = polyorder-1
                end select
                call preassemb
              case('STOKES')
                !Do not call preassemb
              case default
                call preassemb
            end select
          end if
!
          select case(physics)
!  Fluid no-adaption
            case('STOKES')
              call solver_stokes
              return
            case('FLUID')
              call solver_fluid
!
!  Transient without adaptation
            case('TRANHEATTR')
              call assemblyandsolvetrans(epsgl,resgl)
!
!  Default, assemble and solve for the first time.
            case default
              call assemblyandsolve(epsgl,resgl)
          end select
!
!  do the p-adaptation
        case ('P_ADAPT')
          call getsetting('PHYSICS_MODE',physics)
!
!  Allocate vector for polynomial degrees of elements.
          if (.not. associated(ep)) then
            allocate (ep(n,nnat))
!  Get and set starting value for polynomial degree from FEMsettings.txt
            call getsetting('POLYORDER',polyorder)
            ep(1:n,1:nnat)=polyorder
            select case(physics)
              case('FLUIDINCOMPR','FLUIDELMASSCON')
                select case(polyorder)
                  case(1,2)
                    ep(1:n,1:2) = 2
                    ep(1:n,3)   = 1
                  case default
                    ep(1:n,3)   = polyorder-1
                end select
              case default
                !No special element order modification
            end select
            call preassemb
!  Assemble,solve and calculate residual for the first time.
            call assemblyandsolve(epsgl,resgl)
          end if
          call padapt(ext,int,epsgl)
!
!  do the h-adaptation
        case ('H_ADAPT')
          print*,'/----------------------------------------------------------\'
          print*,'+ Doing H-Adaptation...'
          call hadapt(ext,int,ende,error,ende2,epsgl)
!  write new mesh to file and calculate error
          print*, '+ Writing new mesh to file...'
          call basout(ok,n,p,e,xn,yn,kzi,geb,en)
          allocate(res(n,nnat))
          print*, '+ Computing residual...'
          call residual(ext,int,.false.,errcode,res,sumres,sumref)
          deallocate(res)
          print*,'+ H-Adaptation complete.'
          print*,'\----------------------------------------------------------/'
!
!  do the hp-adaptation
        case ('HP_ADAPT')
          print*,'/----------------------------------------------------------\'  
          print*,'+ Doing HP-Adaptation...'
          call getsetting('HP_ALGORITHM',algorithm)
          algorithm = algorithm(1:len_trim(algorithm))
          select case (algorithm)
            case ('KEYPOINT', 'TOP5', 'KP_PHASELAG')
              call hpadapt(ext,int,ende,error,ende2,epsgl)
            case ('HEUVELINE')
              print*,'hp-adaptation: HEUVELINE not available'
              stop
            case ('MELENK')
              print*,'hp-adaptation: MELENK not available'
              stop
            case default
              call hpadapt(ext,int,ende,error,ende2,epsgl)
          end select
!  write new mesh to file and calculate error
          print*, '+ Writing new mesh to file...'  
          call basout(ok,n,p,e,xn,yn,kzi,geb,en)
          allocate(res(n,nnat))
          print*, '+ Computing residual...'
          call residual(ext,int,.false.,errcode,res,sumres,sumref)
          deallocate(res)
          print*,'+ HP-Adaptation complete.'
          print*,'\----------------------------------------------------------/'
      end select
!
      end subroutine adaptation
