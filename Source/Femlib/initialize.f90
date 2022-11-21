      subroutine initialize()
      use femtypes
      use feminterface, only: zanfpp, zeit, getsetting, readsetting
      use feminterface, only: setprocesspriority, physicsinfo
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
!    $Revision: 1.15 $
!    $Date: 2015/04/01 10:49:28 $
!    $Author: juryanatzki $
!
!  perform initial settings
!
      integer (I4B) ::  bildnr, priority
      real (DP) :: freq, lam0, h
      character (len=200) path
!      real (DP) :: xmin, xmax, ymin, ymax, h
!
!  read in system info (Windows or Linux) and store it under the global variable
!  whatsystem
      call getsetting('WHATSYSTEM',whatsystem)
!
!  initialize the graphics window
!  Driver is NULL for initialization purpose only
      xmin=0._DP
      xmax=1._DP
      ymin=0._DP
      ymax=1._DP
      call zanfpp('/NULL',xmin,xmax,ymin,ymax,h,bildnr)
      call zeit(' ')
!
      priority=1
      call setprocesspriority(priority)
!
!  display project path
      call getsetting('PROJECTPATH',path)
      print*, 'actual project-path is :   ',path(1:len_trim(path))
!  get the number of natures
      call getsetting('PHYSICS_MODE',physics)
      call physicsinfo(physics,nnat,2)
!
!  call settings for wavelength and frequency
      call getsetting('WAVELENGTH',lam0)
      call getsetting('FREQUENCY',freq)
!  if one value is 0, it is not used. frequency is stronger than wavelength
      if (freq .ge. 0._DP)  then
        omega = 2._DP*pi*freq
      else if (lam0 .gt. 0._DP) then
        omega = (c0 / lam0)*2._DP*pi
      else
        omega = 0._DP
      end if
!
!  read user parameters from file and store as global variable vector userparam
      call getsetting('USERPARAM1',userparam(1))
      call getsetting('USERPARAM2',userparam(2))
      call getsetting('USERPARAM3',userparam(3))
      call getsetting('USERPARAM4',userparam(4))
      call getsetting('USERPARAM5',userparam(5))
      call getsetting('USERPARAM6',userparam(6))
      call getsetting('USERPARAM7',userparam(7))
      call getsetting('USERPARAM8',userparam(8))
      call getsetting('USERPARAM9',userparam(9))
      call getsetting('USERPARAM10',userparam(10))
!
      return
      end subroutine initialize
