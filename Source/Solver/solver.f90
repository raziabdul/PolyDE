      program solver
      use feminterface, only: getsetting, initialize, meshinput
      use feminterface, only: zeit, lin, lout
      use femtypes
      use globalvariables, only: x, fem_accuracy, physics, omega, pi, &
                               & c0, n, ndof, nnat
      use mpmodule, only: initdata
!  Release global allocatables and pointers
      use globalvariables, only: ep, eg, e, en, geb, kzi, xn, yn
      use globalvariables, only: kzrb, xbk, ybk, zki, zrb, alrb, btrb, matzif
      use matconstants, only: parameternames, param, maxmat

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
!    $Revision: 1.41 $
!    $Date: 2015/04/01 11:03:47 $
!    $Author: juryanatzki $
!
      real (DP) jdmesh, jdsolu
      real (DP) error, epsgl, resgl, resnl, omegasol
      logical ok, reell, gilt, eleinfo, istlin, konfom, zyl
      integer inat, i

      call zeit('+ solver_start')
!
      call initialize()
!  check for use of multiphysics and read data if yes
      if ((physics .eq. 'HEATTR') .or. (physics .eq. 'TEWAVE') .or. (physics .eq. 'TMWAVE')) then
        call initdata()
      end if
!
      print "(a,g10.3/,a,g10.3)","frequency        : ",omega/(2._DP*pi)   &
     &                        ,"angular frequency: ",omega
      if (omega .gt. 0) then
        print "(a,g10.3)","vacuum wavelength: ",(c0 / omega)*2._DP*pi
      end if
!
      print*,''
      print*,'+ Reading mesh input...'
      call meshinput(ok,jdmesh,reell,gilt,istlin,konfom,zyl,epsgl,resgl,&
     &               resnl,error)
!  set error to 100% initially
      fem_accuracy=100._DP
!
      print*,'+ Reading solution vector from existing solution file...'
      call lin(ok, jdmesh, jdsolu, gilt, eleinfo, n, ndof, omegasol,    &
     &         epsgl, nnat, ep, eg, x, fem_accuracy)
!
      call zeit('reading mesh input')
!
      print*,'+ Executing adaptation algorithm...'
      call adaptation(epsgl)
!
      call getsetting('PHYSICS_MODE',physics)
!
      if (physics .NE.'FLUID' .AND. physics .NE.'STOKES' ) then
!  write solution to file
        print*, '+ Writing solution to file...'
        call lout(.true., .true., n, ndof, omega, epsgl, nnat, ep,     &
     &              eg, x, fem_accuracy)
      end if
!
!
      call zeit('+ solver_end')
!
      print*,''
      print*,'==> Solver completed execution. <=='
!
!  need to release all the global pointers still associated
      if (associated(e))   deallocate(e)
      if (associated(en))  deallocate(en)
      if (associated(geb)) deallocate(geb)
      if (associated(ep))  deallocate(ep)
      if (associated(kzi)) deallocate(kzi)
      if (associated(xn))  deallocate(xn)
      if (associated(yn))  deallocate(yn)

      if (associated(x))  deallocate(x)

      if (associated(kzrb))   deallocate(kzrb)
      if (associated(xbk))    deallocate(xbk)
      if (associated(ybk))    deallocate(ybk)
      if (associated(zki))    deallocate(zki)
      if (associated(zrb))    deallocate(zrb)
      if (associated(alrb))   deallocate(alrb)
      if (associated(btrb))   deallocate(btrb)
      if (associated(matzif)) deallocate(matzif)

      if (allocated(parameternames)) deallocate(parameternames)
!
!  Check whether the pointer is associated with a target
      if ( associated(eg) ) then

        do inat = 1, nnat
          do i = 1, size(eg, DIM=1)
            deallocate(eg(i,inat)%d)
          end do
        end do
        
        deallocate(eg) ! and nullify
      end if
!
!  allocate the array of parameters and copy the standard values to it
      do i=0,maxmat
        deallocate( param(i)%d )
      end do
!
!
      stop
      end