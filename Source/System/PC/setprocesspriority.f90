      subroutine setprocesspriority(priority)
      use femtypes
!  following libraries are for Compaq Visual Fortran 6
      use dflib
      use dfwin
!  following library is for Intel Visual Fortran 8
!      use ifwin
      implicit none
      integer (I4B) priority
      intent (in) :: priority
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
!    $Date: 2006/07/08 00:00:25 $
!    $Author: r_abdul $
!
!  set the priority class of the main process
!  Input:
!                    / 1: IDLE_PRIORITY_CLASS
!            priority- 2: NORMAL_PRIORITY_CLASS
!                    \ 3: HIGH_PRIORITY_CLASS    
!                      4 :REALTIME_PRIORITY_CLASS
!  ATTENTION: don't use priority 3 or 4, the highest priority will exclusively execute this 
!             process, even the mouse would be unresponsive
!
      integer (4) pocesshandle
      logical res
!
      pocesshandle=GetCurrentProcess()
      select case (priority)
      case (1)
        res=SetPriorityClass(pocesshandle,IDLE_PRIORITY_CLASS)
      case (2)
        res=SetPriorityClass(pocesshandle,NORMAL_PRIORITY_CLASS)
      case (3)
        res=SetPriorityClass(pocesshandle,HIGH_PRIORITY_CLASS)
      case (4)
        res=SetPriorityClass(pocesshandle,REALTIME_PRIORITY_CLASS)
      case default
        res=SetPriorityClass(pocesshandle,IDLE_PRIORITY_CLASS)
      end select
!      priority=GetPriorityClass(pocesshandle)
      end subroutine setprocesspriority
