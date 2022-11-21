      subroutine juld2cal(jd,year,mon,day,month,dayofweek,hour,minute,second)
      use femtypes
      use feminterface, only:
      implicit none
      real (DP) jd
      integer (I4B) year, mon, day, hour, minute, second
      character(len=*) dayofweek, month
      intent (in) jd
      intent (out) year, mon, day, month, dayofweek, hour, minute, second
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
!    $Revision: 1.3 $
!    $Date: 2006/07/07 21:08:36 $
!    $Author: r_abdul $
!
!  convert julian day to calendar day
!  Input:
!            jd           julian day (including fraction for hour, minutes, seconds)
!  Output
!            year         year
!            mon          month 1..12
!            day          day 1..31
!            month        three character month string
!            dayofweek    three character day of week string
!            hour         hour  0.23
!            minute       minute  0..59
!            second       second  0..59
!   
!  local variables
      integer (I4B) l, n, ndw
      real (DP) fract, d
!
      l=int(jd)+68569
      n=4*l/146097
      l=l-(146097*n+3.)/4
      year=4000.*(l+1)/1461001
      l=l-1461*year/4.+31
      mon=80.*l/2447
!  day
      day=l-2447*mon/80
      l=mon/11
!  month
      mon=mon+2-12*l
!  year
      year=100*(n-49)+year+l
      select case (mon)
      case (1)
        month='Jan'
      case (2)
        month='Feb'
      case (3)
        month='Mar'
      case (4)
        month='Apr'
      case (5)
        month='May'
      case (6)
        month='Jun'
      case (7)
        month='Jul'
      case (8)
        month='Aug'
      case (9)
        month='Sep'
      case (10)
        month='Oct'
      case (11)
        month='Nov'
      case (12)
        month='Dec'
      end select
      d=(int(jd)+1.)/7.
      ndw=int((d-int(d))*7+0.5)
      select case (ndw)
      case(0)
        dayofweek='Sun'
      case(1)
        dayofweek='Mon'
      case(2)
        dayofweek='Tue'
      case(3)
        dayofweek='Wed'
      case(4)
        dayofweek='Thu'
      case(5)
        dayofweek='Fri'
      case(6)
        dayofweek='Sat'
      end select
!  compute hour
      hour = int( (jd -int(jd)) * 24.0)
      fract = ( jd -int(jd) - real( hour) / 24.0) * 24.0
!  compute minute
      minute = int( fract * 60.0)
      fract = ( fract - real( minute) / 60.0) * 60.0
!  compute second
      second = int( fract * 60.0)
      return
      end subroutine juld2cal



      subroutine cal2juld(jd,year,mon,day,hour,minute,second)
      use femtypes
      use feminterface, only:
      implicit none
      real (DP) jd
      integer (I4B) year, mon, day, hour, minute, second
      intent (in) year, mon, day, hour, minute, second
      intent (out) jd
!  convert calendar day to julian day
!  Input:
!            year         year
!            mon          month 1..12
!            day          day 1..31
!            hour         hour  0.23
!            minute       minute  0..59
!            second       second  0..59
!  Output
!            jd           julian day (including fraction for hour, minutes, seconds)
!   
!  local variables
      integer (I4B) j2, ju
! 
      j2 = int( (mon - 14)/12 );
      ju = day - 32075 + int(1461 * ( year + 4800 + j2 ) / 4 );
      ju = ju + int( 367 * (mon - 2 - j2 * 12) / 12);
      ju = ju - int(3 * int( (year + 4900 + j2) / 100) / 4);
      jd = real(ju,DP)+(real(hour*3600+minute*60+second,DP)+0.5)/86400._DP
      return
      end subroutine cal2juld



      real (DP) function timestamp()
      use femtypes
      use feminterface, only: cal2juld
      implicit none
!   
!  local variables
      integer dt(8)
      character (len=12) clock(3)
!
      call date_and_time(clock(1),clock(2),clock(3),dt)
      call cal2juld(timestamp,dt(1),dt(2),dt(3),dt(5),dt(6),dt(7))
      return
      end function timestamp

