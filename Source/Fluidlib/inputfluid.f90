subroutine inputfluid
use femtypes
use feminterface, only: getsetting
use fluidvariables
use globalvariables
implicit none
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
!    $Revision: 1.4 $
!    $Date: 2006/07/07 22:34:04 $
!    $Author: r_abdul $
!
!********************************************************************************
!
!  Read Control Parameters for Fluid Computation
!
!********************************************************************************
! local variables
!********************************************************************************
integer :: opttemp
integer(I4B) :: i, ios !j, idof
integer(I4B) :: unitid !, unit2
character(len=100):: text
character(len=200):: path
logical :: opt1
!
external    ::    grglun
!********************************************************************************
call getsetting('PROJECTPATH',path)
call grglun(unitid)
open (unitid,file=path(1:len_trim(path))//'inputfluid.dat',      &
     &  form='formatted',status='old',action='READ',iostat=ios)

! reading compressibility option: 0 = incompressible, 1 = compressible
read(unitid,*) text
read(unitid,*) opttemp

if (opttemp == 1) then
   optcomprs = .TRUE.
   read(unitid,*) optstrm
else
   optcomprs = .FALSE.
   read(unitid,*) text
   read(unitid,*) optener, optstokes, optstrm, optdim
end if

! reading initial-value option
! 0 = read all assumed data
! 1 = read existing data from file
! -1 = assume all rassum, uassum, vassum = 1.0
read(unitid,*) text
read(unitid,*) iopt

read(unitid,*) text
read(unitid,*) rassum, uassum, vassum

if (optcomprs) then
!*******************************************************************************************
! read information for compressible
!*******************************************************************************************
! reading gamma= Cp/Cv and artificial diffusion constant(0.5) 
read(unitid,*) text
read(unitid,*) (conin(i),i=1,2)

! reading free stream values of rho, u, v, p and Mach number
! Note: Pinf is anyway calculated from others. Pinf = cinf(1)*(gamma-1)*Tinf/gamma
read(unitid,*) text
read(unitid,*) (cinf(i),i=1,5)

! reading time-step related parameters
! no. of TS to run      - ntime
! no. of TS ran so far  - istep
! non-dim time elasped  - timt 
! flag for ts calculation(1=local, 0=global, or -1=prescribed TS) - optts
! no. of mass iterations for consistent mass matrix calculation - niter
! no. of steps to write residues - iwrite
! fixed time step value(only if optts <= -1) - dtfix
read(unitid,*) text
read(unitid,*) ntime, istep, timt, optts, niter, iwrite, dtfix  ! compressible

! safety factor
read(unitid,*) text
read(unitid,*) csafm

! constant for residual smoothing and no. of smoothing loops
read(unitid,*) text
read(unitid,*) csmoo, nsmoo

! relaxation parameters for velocities and pressure (theta1 and theta2) 
! Note: compressible   - theta1 = 0.5, theta2 = 0.0
read(unitid,*) text
read(unitid,*) (theta(i),i=1,2)

! Reynolds number, Prandtl number, ratio of stagnation temp. to free stream temp.,
! free stream temperature(in Kelvin) (in Rankine see Nithairasu)
read(unitid,*) text
read(unitid,*) re, pr, ratio, tfree
!*******************************************************************************************
else
!*******************************************************************************************
! read information for incompressible
!*******************************************************************************************

! reading free stream values of p, u, v, T 
read(unitid,*) text
read(unitid,*) (cinf(i),i=1,4)

! reading time-step related parameters
! ntime  = no. of TS to run
! optts  = flag for timestep calculation(1=local, 0=global, or -1=prescribed TS),
! iwrite = no. of steps to write residues, 
! dtfix  = fixed time step value (used only when optts <= -1)
read(unitid,*) text
read(unitid,*) ntime, optts, iwrite, dtfix

! safety factor (used only when optts = 1)
read(unitid,*) text
read(unitid,*) csafm

! relaxation parameters for velocities and pressure (theta1 and theta2) 
! Note: incompressible explicit - theta1 = 0.5 or 1.0 , theta2 = 0.0
!       incompressible implicit - theta1 = 1.0, theta2 = 1.0
read(unitid,*) text
read(unitid,*) (theta(i),i=1,2)

! Reynolds number, Prandtl number, Rayleigh number, Richadson number (1/Froude number)
read(unitid,*) text
read(unitid,*) re, pr, ra, ri

read(unitid,*) text
read(unitid,*) optnat

read(unitid,*) text
read(unitid,*) optmagn

opt1 = (optdim .AND. optener)
if (optmagn .OR. opt1) then
   read(unitid,*) text
   read(unitid,*) twall, lchar
end if

!*******************************************************************************************
end if

close(unitid)
end subroutine inputfluid
