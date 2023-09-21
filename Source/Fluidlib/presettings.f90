subroutine presettings
use femtypes
use feminterface
use fluidvariables
use globalvariables
use matconstants
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
!    $Revision: 1.5 $
!    $Date: 2006/07/07 22:34:04 $
!    $Author: r_abdul $
!
!******************************************************************************************
!
!     Settings global variables and related control parameters
!
!******************************************************************************************
! local variable
!******************************************************************************************
integer(I4B) :: i,j,ip
integer(I4B) :: ios, unitid
real(DP) :: assume
character(len=200):: path
!
external    ::    grglun
!******************************************************************************************

if (optcomprs) then
!******************************************************************************************
!   preliminary set-up for compressible 
!******************************************************************************************
! c = sqrt(u^2+v^2)/Ma
cc = dsqrt(cinf(2)**2 + cinf(3)**2)/cinf(5)

! free stream press = c^2*rho/gamma
pinf = (cc**2)*cinf(1)/gamma

! free stream temperature = c^2/(gamma - 1)
tinf = cc**2/( conin(1)-1.0_DP )

! option for initial values
! 0 = read all assumed data
! 1 = read existing data from file
! -1 = assume all rassum, uassum, vassum = 1.0
if(iopt.eq.0) then
   do ip = 1, ndof
      unkno(1,ip) = rassum
      unkno(2,ip) = uassum
      unkno(3,ip) = vassum
      unkno(4,ip) = tinf/conin(1) + 0.5_DP*( uassum**2 + vassum**2 )
      tt1(ip)     = tinf
      sound(ip)   = cc
      !pres(ip)    = (sound(ip)**2)*unkno(1,ip)/conin(1)
      pres(ip)    = (cc**2)*rassum/conin(1)
      pres1(ip)   = pres(ip)
   enddo !ip
! without loop ip 
!   unkno(1,:) = rassum
!   unkno(2,:) = uassum
!   unkno(3,:) = vassum
!   unkno(4,:) = tinf/conin(1) + 0.5_DP*( uassum**2 + vassum**2 )
!   tt1(:)     = tinf
!   sound(:)   = cc
!   pres(:)    = (cc**2)*rassum/conin(1)
!   pres1(:)   = pres(:)
else if (iopt .EQ. 1) then
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
      open (unitid,file=path(1:len_trim(path))//'solution.out',  &
            form='formatted',status='old',action='READ',iostat=ios)
      read(unitid,*) istep
      do i = 1,ndof
         read(unitid,*) ip, (unkno(j,i),j=1,nvar)
      end do ! i
      close(unitid)     
      do ip = 1, ndof
         tt1(ip)   = conin(1)*( unkno(4,ip) - 0.5_DP*(unkno(2,ip)**2+ unkno(3,ip)**2) )
         pres(ip)  = ( conin(1)-1 )*unkno(1,ip)*tt1(ip)/conin(1)
         pres1(ip) = pres(ip)
      enddo !ip
else if (iopt .EQ. -1) then
   assume     = 1.0_DP
   unkno(1,:) = assume
   unkno(2,:) = assume
   unkno(3,:) = assume
   if (optener) unkno(4,:) = assume
end if ! iopt
               
invre = 1.0_DP/re
gamma = conin(1)
gam1 = gamma - 1.0_DP
su = gam1*(cinf(5)**2)
!sucon = 198.6_DP/(su*tfree)  ! Tfree in Rankine
sucon = 110.4_DP/(su*tfree)   ! Tfree in Kelvin
suth = (1.0_DP/su + sucon)*(su**1.5)
!******************************************************************************************
else
!******************************************************************************************
!   preliminary set-up for incompressible 
!******************************************************************************************
! option for initial values
! 0 = read all assumed data
! 1 = read existing data from file
! -1 = assume all rassum, uassum, vassum = 1.0
if(iopt.eq.0) then
   unkno(1,:) = cinf(1)             ! set pressure as the free stream value
   unkno(2,:) = uassum              ! set u as the given value
   unkno(3,:) = vassum              ! set v as the given value
   if (optener) unkno(4,:) = rassum ! set temperature as the given value 
else if (iopt .EQ. 1) then
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
      open (unitid,file=path(1:len_trim(path))//'solution.out',  &
            form='formatted',status='old',action='READ',iostat=ios)
      read(unitid,*) istep
      do i = 1,ndof
         read(unitid,*) ip, (unkno(j,i),j=1,nvar)
      end do ! i
      close(unitid)     
else if (iopt .EQ. -1) then
   assume     = 1.0_DP
   unkno(1,:) = assume
   unkno(2,:) = assume
   unkno(3,:) = assume
   if (optener) unkno(4,:) = assume
end if ! iopt


if (optdim) then
! visco = mu or eta (param(matindex)%d(numparam) for fluid)
! see details in setstandardvalues.f90
   visco = param(1)%d(1)
   invrho = 1.0_DP/param(1)%d(2)
   invre = visco*invrho
   thcond = param(1)%d(3)
   if (optener) then
   ! thermal diffusivity = thermal conductivity/(density*specific heat)
      thdiff = thcond*invrho/param(1)%d(4)
   ! free-stream temperature
      tfree = cinf(4)
   ! the total temperature difference between hot surface and fluid.
      tmpdiff = twall-tfree
      tfilm = (twall+tfree)/2
      thexpn = 1/tfilm 
   end if
end if

if (optener) then
! coupled with energy equation   
   if (.NOT. optdim) then
   ! non-dimensional
      invrho = 1.0_DP     
      if (optnat) then
      ! for natural convection
         invre = pr
      else
      ! forced convection
         invre = 1.0_DP/re
      end if
      thdiff = invre/pr
   end if
else
! only fluid equations
   if (.NOT. optdim) then
      invre = 1.0_DP/re
      invrho = 1.0_DP
   end if
end if

if (optmagn) then
! external magnetic force acts as the bouyancy force in natural convection
! temporary setting for oxygen
   curie = 1.2923e-5_DP

! only order of T
!   magconst = -0.5_DP*mu0*curie/32e-3_DP
! T square
   magconst = -0.5_DP*invrho*mu0*curie/gasconst

   if (optdim) then
      if (cinf(4) .le. 1.0_DP) then
         write(*,*) 'check value of free-stream temperature again'
      end if      
   else
    ! dimensionless
      invre = pr
      invrho = 1.0_DP
    ! check terms again and compute the correct one
      magconst = -0.5_DP*lchar*mu0*curie/gasconst
    ! free-stream temperature
      tfree = cinf(4)
    ! the total temperature difference between hot surface and fluid.
      tmpdiff = twall-tfree
   end if
end if
!******************************************************************************************
end if

end subroutine presettings