      subroutine snkrgr(xa,ya,xb,yb,xc,yc,xd,yd,xm,ym,rl1,rl2,t1,t2,err)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) xa,ya,xb,yb,xc,yc,xd,yd,xm,ym,rl1,rl2,t1,t2
      integer (I4B) err
      intent (in) :: xa, ya, xb, yb, xc, yc, xd, yd, xm, ym
      intent (out) :: rl1, rl2, t1, t2, err
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
!    $Date: 2006/07/07 21:08:37 $
!    $Author: r_abdul $
!
!  Computation of the intersection point between a straight line
!  passing through the points (xa,ya)-(xb,yb) and an arc in mathematical
!  positive orientation from point (xc,yc) to the point (xd,yd)
!  having the centre point at (xm,ym)
!  
!    intersection point: = x = xa + rl * (xb - xa)   (line)
!                        = y = ya + rl * (yb - ya)
!    for the arc:
!                        x = xm  + R * cos( alpha + t * (beta-alpha) )
!                        y = ym  + R * sin( alpha + t * (beta-alpha) )
!
!   There are none, one or two intersections rl1, rl2, t1, t2
!   (the tangent is a special case of two idential inresection points)
!   
!   In the case, that no intersection is present, the return value is set to rl=huge(1._DP)
!   t than results in the nearest point between arc and line
!
!   For  rl in [0,1] the intersection lies on the line A-B
!   For  t  in [0,1] the intersection lies on the arc C-D
!
!   err    =  0   if no erroro occured
!            10   the straight line from A to B is a point 
!
!  lokal variables
!
      real (DP) pi, zweipi
      real (DP) dx1m,dy1m,dx2m,dy2m
      real (DP) dxba,dyba,dxam,dyam,dxcm,dycm,dxdm,dydm
      real (DP) r2,r2k,dsk,zaehl,xnenn, a, bhalbe, c, q
      real (DP) alpha, gamma1, gamma2, tmax
      parameter (pi=3.14159265358979323846264338327950288419716_DP)
      parameter (zweipi=2._DP*pi)
!
      dxba=xb-xa
      dyba=yb-ya
      dxam=xa-xm
      dyam=ya-ym
      dxcm=xc-xm
      dycm=yc-ym
      dxdm=xd-xm
      dydm=yd-ym
!  radius of the circle
      r2=dxcm*dxcm+dycm*dycm
!  again (check)
      r2k=dxdm*dxdm+dydm*dydm
!     if (abs(r2-r2k)/(r2+r2k).gt.01) then
!       print*,'*** startpoit and endpoint are not on the same radius'
!     end if
!  take the mean value of radii computed with start- and endpoint
      r2=(r2+r2k)*0.5_DP
!  length  of the staight line
      a=dxba**2+dyba**2
      if (a .le. 100.*tiny(1._DP)) then
!  this does not give sense, the line A-B is a point
        err=1
        return
      end if
      bhalbe=(dxam*dxba+dyam*dyba)
      c=dxam**2+dyam**2-r2
!  discriminant (a measure fo the distance between line and circle)
      dsk=bhalbe**2-a*c
      if (dsk.lt.0._DP) then
!  no intersection between line and arc
        rl1=-bhalbe/a
        rl2=rl1
      else
!  solving for rl on the line (quadratic equation)
!  we allways have rl1 > rl2
!  distinguish between bhalbe geater or smaller zero for better numerical stability
!  see numerical recipes
        if (bhalbe .ge. 0._DP) then 
          q=-(bhalbe+sqrt(dsk))
          rl1=q/a
          rl2=c/q
        else 
          q=-(bhalbe-sqrt(dsk))
          rl1=c/q
          rl2=q/a
        end if
      end if
!  angle between  point D and C seen from the centre
      zaehl=dxcm*dydm-dxdm*dycm
      xnenn=dxdm*dxcm+dydm*dycm
      alpha=atan2(zaehl,xnenn)
      if (alpha.le.0._DP) alpha=alpha+zweipi
!  relative intersection co-ordinates
      dx1m=dxam+rl1*dxba
      dy1m=dyam+rl1*dyba
!  angle between the intersection point 1 and C
      zaehl=dxcm*dy1m-dx1m*dycm
      if (abs(zaehl) .le. spacing(max(dxcm*dy1m,dx1m*dycm))) then
        gamma1=0._DP
      else
        xnenn=dx1m*dxcm+dy1m*dycm
        gamma1=atan2(zaehl,xnenn)
        if (gamma1.lt.0._DP) gamma1=gamma1+zweipi
      end if
      dx2m=dxam+rl2*dxba
      dy2m=dyam+rl2*dyba
!  angle between the intersection point 2 and C
      zaehl=dxcm*dy2m-dx2m*dycm
      if (abs(zaehl) .le. spacing(max(dxcm*dy2m,dx2m*dycm))) then
        gamma2=0._DP
      else
        xnenn=dx2m*dxcm+dy2m*dycm
        gamma2=atan2(zaehl,xnenn)
        if (gamma2.lt.0._DP) gamma2=gamma2+zweipi
      end if
      t1=gamma1/alpha
      t2=gamma2/alpha
      tmax=zweipi/alpha
      if (t1.ge.tmax) then
        t1=t1-tmax
      end if
      if (t2.ge.tmax) then
        t2=t2-tmax
      end if
      if (dsk .lt. 0._DP) then
        rl1=huge(1._DP)
        rl2=huge(1._DP)
      end if
!  correct special cases (rounding error)
      if (xa.eq.xc .and. ya.eq.yc) then
        if (abs(rl1) .lt. abs(rl2)) then
          rl1=0._DP
          t1=0._DP
        else
          rl2=0._DP
          t2=0._DP
        end if
      else if (xa.eq.xd .and. ya.eq.yd) then
        if (abs(rl1) .lt. abs(rl2)) then
          rl1=0._DP
          t1=1._DP
        else
          rl2=0._DP
          t2=1._DP
        end if
      end if
      if (xb.eq.xc .and. yb.eq.yc) then
        if (abs(rl1-1._DP) .lt. abs(rl2-1._DP)) then
          rl1=1._DP
          t1=0._DP
        else
          rl2=1._DP
          t2=0._DP
        end if
      else if (xb.eq.xd .and. yb.eq.yd) then
        if (abs(rl1-1._DP) .lt. abs(rl2-1._DP)) then
          rl1=1._DP
          t1=1._DP
        else
          rl2=1._DP
          t2=1._DP
        end if
      end if
      err=0
      return
      end subroutine snkrgr