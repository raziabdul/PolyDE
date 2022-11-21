      subroutine zeit(text)
      use feminterface, only: second
      use femtypes
      implicit none
      character (len=*) text
      intent (in) :: text
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
!    $Revision: 1.8 $
!    $Date: 2011/08/16 15:12:26 $
!    $Author: m_kasper $
!  Output of total CPU-time  and time elapsed since last call
!  Input:  Text   Output text (comment)
!
!  lokal variables
      real(DP) ecputime, exetime, sumtime, eexetime
      real(DP) , save :: oldtime, starttime, oldexetime
      integer(4) :: v(8)
      character(5) zone
      character(8) date
      character(10) time
      logical , save ::  first=.true.
!
      call cpu_time(sumtime)
      call date_and_time(date,time,zone,v)
      if (first ) then
        oldtime=sumtime
        starttime=(((v(3)-1)*24+v(5))*60+v(6))*60+v(7)+v(8)/1000._DP
        oldexetime=0._DP
        first=.false.
      else if (sumtime .gt. 0._DP) then
!  elapsed CPUtime since last call
        ecputime=sumtime-oldtime
!  total EXEtime
        exetime=(((v(3)-1)*24+v(5))*60+v(6))*60+v(7)+v(8)/1000._DP-starttime
!  elapsed EXEtime since last call
        eexetime=exetime-oldexetime
        write (*,3) adjustl(text), ecputime, eexetime, sumtime, exetime
3       format(a,f7.2'/'f7.2,' sec',t54,'total CPU/EXE time:',f8.2' / 'f8.2,' sec')
        oldtime=sumtime
        oldexetime=exetime
      end if
      return
      end