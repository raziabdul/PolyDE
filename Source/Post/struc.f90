      subroutine struc(gzz,zki,zrb,xbk,ybk,r,g,b)
      use feminterface, only: inkrp,circle
      use femtypes
      integer (I4B) gzz,zki(:,:),zrb(:,:)
      real (DP) xbk(:),ybk(:)
      real (SP) r, g, b
      intent (in):: gzz,zki,zrb,xbk,ybk,r,g,b
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
!    $Revision: 1.7 $
!    $Date: 2008/12/22 13:16:51 $
!    $Author: m_kasper $
!
!  draw the structure (branches)
!
      real (DP) x1,y1,x2,y2,xmitte,ymitte,phi1,phi2,flae2
      real (DP) pi,radius,xinnen,yinnen,rinnen
      integer (I4B) i,oldindex
      real (SP) ocr,ocg,ocb
      parameter (pi=3.141592653589793_DP)

!
      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgqwin, pgsvp, pgwnad

      !  save color status
      call pgqci(oldindex)
      call pgsci(20)
      call pgqcr(20, ocr, ocg, ocb)
!  set color  e.g. yellow= 1.00, 1.00, 0.00 (R,G,B)
      call pgscr(20,  r, g, b)
!  graphical output of stucture
      do i=1,gzz
        if (zrb(i,1).lt.400) then
          x1=xbk(zki(1,i))
          y1=ybk(zki(1,i))
          x2=xbk(zki(2,i))
          y2=ybk(zki(2,i))
          if (zki(3,i).eq.0) then
            call pgmove(sngl(x1),sngl(y1))
            call pgdraw(sngl(x2),sngl(y2))
          else
            if (zki(3,i).lt.0) then
              xmitte=xbk(-zki(3,i))
              ymitte=ybk(-zki(3,i))
              radius=sqrt((x1-xmitte)**2+(y1-ymitte)**2)
            else
              call inkrp(x1,y1,x2,y2,xbk(zki(3,i)),ybk(zki(3,i)),       &
     &          xinnen,yinnen,flae2,rinnen,xmitte,ymitte,radius)
            end if
            phi1=atan2(y1-ymitte,x1-xmitte)
            phi2=atan2(y2-ymitte,x2-xmitte)
            if (phi2.lt.phi1) phi2=phi2+2._DP*pi
            call circle(sngl(x1),sngl(y1),sngl(phi1/pi*180._DP),        &
     &            sngl(phi2/pi*180._DP),sngl(radius),sngl(radius))
          end if
        end if
      end do
      call pgsci(oldindex)
      call pgscr(20, ocr, ocg, ocb)
      return
      end subroutine struc
