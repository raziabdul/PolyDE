subroutine readnetin_mn
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
!    $Revision: 1.6 $
!    $Date: 2014/07/15 13:02:35 $
!    $Author: m_kasper $
!
!**********************************************************************************
!
!  Read netinmod.acd
!
!**********************************************************************************
! status    : check again for the read data if there is any conflict with 
!             "solver" of Polyde
! modified  : 23.08.2005
! remarks   : taken from F90_stokethermal, input_geom_mn which is adapted 
!             from lese.f90 in Triang
!           - read geometrical data from "netinmod.acd" instead of "netin.acd"
!      TO DO implement multinature boundary conditions
!**********************************************************************************
!   zrb	            Randbedingungen der Zweige
!                     0 -  99: Dirichlet
!                   100      : Mixed Boundary Conditions
!                   200 - 299: Neumannn
!                   300 - 399: Kontur
!                   400 - 499: innen
!   alrb, btrb  (realle) Parameter der allgemeinen Randbedingungen (Boundary values)
!   nkb                  No. of branches which a keypoint locates on
!**********************************************************************************
!   local variables
!**********************************************************************************
integer(I4B), pointer :: zrbbc(:,:), bccode(:)
integer(I4B), pointer :: krb(:)
integer(I4B) :: ios,i,j,k,ivar, ipbc
integer(I4B) :: dummy, unitid
real(DP), pointer :: alrbre(:,:), btrbre(:,:), alrbim(:,:), btrbim(:,:)
real(DP) :: varptemp
character(len=200) :: path
!
external    ::    grglun
!**********************************************************************************
call getsetting('PROJECTPATH',path)
call grglun(unitid)
open(unitid,file=path(1:len_trim(path))//'netinmod.acd', &
     form='formatted',status='old',action='READ',iostat=ios)
!**********************************************************************************
! read boundary conditions for multinature in case zrb(i) > 100
!**********************************************************************************
read(unitid,*) nbc   ! number of different multinature bc

allocate(bccode(nbc),zrbbc(nbc,nvar))
allocate(alrbre(nbc,nvar),alrbim(nbc,nvar),btrbre(nbc,nvar),btrbim(nbc,nvar))
allocate(zrbmn(gzz,nvar),alrbmn(gzz,nvar),btrbmn(gzz,nvar))
120 format(4(f10.5))

do i = 1,nbc
   read(unitid,*)  bccode(i), zrbbc(i,:)
   do k = 1,nvar
      read(unitid,*) alrbre(i,k),alrbim(i,k),btrbre(i,k),btrbim(i,k)
!     write(12,120) alrbre(i,k),alrbim(i,k),btrbre(i,k),btrbim(i,k)
   end do
end do

branch: do i = 1,gzz ! loop over branches
   do j = 1,nbc
      if (zrb(i,1) == bccode(j)) then
         zrbmn(i,:) = zrbbc(j,:)
         alrbmn(i,:) = cmplx(alrbre(j,:), alrbim(j,:),DP)
         btrbmn(i,:) = cmplx(btrbre(j,:), btrbim(j,:),DP)
         cycle branch
      else
         zrbmn(i,:) = zrb(i,1)
      end if
   end do ! j
end do branch ! i

deallocate(bccode,zrbbc) 
deallocate(alrbre,alrbim,btrbre,btrbim)
!**********************************************************************************
! set b.c. to keypoints (kp is set for the min BC no.)
!**********************************************************************************
allocate(kzrbmn(gkz,nvar),krb(gkz))
kzrbmn = 0

do ivar = 1,nvar  ! loop over variables
   krb = 999999
   do i = 1, gzz
      if (zrbmn(i,ivar) < krb(zki(1,i))) then
         krb(zki(1,i)) = zrbmn(i,ivar) ! BC of starting keypoint of branches
         kzrbmn(zki(1,i),ivar) = i     ! tell which branch the keypoint is on.
      end if

      if (zrbmn(i,ivar) < krb(zki(2,i))) then
         krb(zki(2,i)) = zrbmn(i,ivar)
         kzrbmn(zki(2,i),ivar) = i
      end if
   end do ! branch i
end do ! ivar
deallocate(krb)
!**********************************************************************************
! read a fixed b.c. on node
!**********************************************************************************
read(unitid,*) npbc
if (npbc .GT. 0) then
   allocate(pkp(npbc),zrbpbc(npbc),varpbc(npbc))
   do ipbc = 1,npbc
      read(unitid,*) pkp(ipbc), zrbpbc(ipbc)
      dummy = zrbpbc(ipbc)
      read(unitid,*) varpbc(ipbc), varptemp
   end do ! ipbc
end if
!**********************************************************************************
close(unitid)
write(*,*) 'finished reading netinmod.acd information'
return
end subroutine readnetin_mn