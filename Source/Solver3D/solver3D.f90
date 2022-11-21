      program solver3D
      use feminterface, only: getsetting, zeit, low2hi
      use feminterface3d, only: initialize3D, readng, solout, readunv, solin_aux
      use feminterface3d, only: solve_adapt3D
      use femtypes
      use globalvariables3D, only : c0, eltype, nnat, numv, omega, pi, vp
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
!    $Revision: 1.24 $
!    $Date: 2015/04/02 10:14:53 $
!    $Author: juryanatzki $
!
!------------------------------------------------------------------------------
!
! This Routine: 'solver3D'
!    Main Program
!
! Input:
!    optional console arguments
!
! Output:
!    Solved FEM Problem
!           OR
!    Error
!
!
!------------------------------------------------------------------------------
!  Internal variables:
      integer (I4B)     :: k, dotpos, polyorder
      character(len=50) :: meshfile
      logical           :: ok
      
      ! We need these for solin_aux
      real (DP)         :: epsgl
      logical           :: gilt, eleinfo
      character(len=3)  :: use_aux_sol
      !__
!__
!_1) removed_
! 2) Initialize the solver
      call initialize3D
      
!____________________________________________________D E L E T E_____
!  print frequency and vacuum wavelength
!      print "(a,g8.3/,a,g8.3)","frequency        : ",omega/(2._DP*pi)&
!                              ,"angular frequency: ",omega
!      if (omega .gt. 0) then
!        print "(a,g8.3)","vacuum wavelength: ",(c0 / omega)*2._DP*pi
!      end if
!____________________________________________________________________
!__
! 3) Read mesh from MESHFILE
      call getsetting('MESHFILE',meshfile)
      dotpos = index(meshfile,'.',BACK=.true.)
      if (dotpos .eq. 0) then 
        print*,'Error: Extension of meshfile: ',meshfile(1:len_trim(meshfile)),' is missing'
      end if
      call low2hi(meshfile,len_trim(meshfile))
      if (meshfile(dotpos+1:dotpos+4) .eq. 'UNV') then
        call readunv(meshfile,ok)
!      else if(  .eq. 'VTK) then
!        call vtkin(ok)                  ! to be implemented
      else
        call readng(meshfile, ok)
      end if
!      do i=1,numv
!        print*,'reading vn(:,',i,') = ',vn(:,i)
!      end do
      if (.not. ok) stop
      call zeit('reading mesh')
      
      !__ TODO: Move this somewhere else (e.g. initialize3d)
      call getsetting('USE_AUX_SOLUTION', use_aux_sol)
      if (use_aux_sol .eq. 'YES') call solin_aux(ok, gilt, eleinfo, epsgl)
      !__

!____________________________________________________ TODO _____
!  We here should read a starting solution from a pervious run, 
!  by a call to solin which then returns
!  the solution vector x,
!  the polynomial orders v,
!  and the assignment of dofs vgdof
!  if such a previous solution is available and matches to the actal mesh
!  otherwise we initialize vp=polyorder (and the solution vector)
!____________________________________________________ TODO _____

!  allocate array of volume element polynomial degrees
!  get and set starting value for polynomial degree from FEMsettings.txt.
      allocate (vp(numv,nnat))
      call getsetting('POLYORDER',polyorder)
      vp(1:numv,:) = polyorder

!  get element type
      call getsetting('ELEMENT_TYPE',eltype)
      call low2hi(eltype,16)
!____________________________________________________________________
!__
! 4) Solve FEM
      call solve_adapt3D
!__
! 5) Write solution to file 'solution' in PROJECTPATH
      call solout(.true., .true., 1.e-8_DP)
!      pause
!
!_End.
      end program solver3D
