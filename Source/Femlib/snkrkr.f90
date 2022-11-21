      subroutine snkrkr(xa,ya,xb,yb,xc,yc,xd,yd,xm1,ym1,xm2,ym2,        &
     &  t1u,t1v,t2u,t2v,xpu,ypu,xpv,ypv,err)
      use feminterface, only: wink, snkrgr
      use femtypes
      implicit none
      real (DP) xa, ya, xb, yb, xc, yc, xd, yd, xm1, ym1, xm2, ym2
      real (DP) t1u, t1v, t2u, t2v, xpu, ypu, xpv, ypv
      integer (I4B) err
      intent (in) :: xa, ya, xb, yb, xc, yc, xd, yd, xm1, ym1, xm2, ym2
      intent (out) :: t1u, t1v, t2u, t2v, xpu, ypu, xpv, ypv
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
!    $Date: 2006/07/07 21:08:37 $
!    $Author: r_abdul $
!
!  Berechnung des Schnittpunktes von einem mathematisch positiv
!  orientierten Kreisbogen vom Punkt (xa,ya) zum Punkt (xb,yb)
!  mit dem Mittelpunkt (xm1,ym1) und einem weitern mathematisch 
!  positiv orientierten Kreisbogen vom Punkt (xc,yc) zum 
!  Punkt (xd,yd) mit dem Mittelpunkt (xm2,ym2)
!
!    fuer den Kreisbogen gilt die Gleichung:
!                   x = xm  + R * cos( alpha+t*(beta-alpha) )
!                   y = ym  + R * sin( alpha+t*(beta-alpha) )
!
!   es existieren entweder keiner oder zwei Schnittpunkte t1u, t1v, t2u, t2v
!    (die Tangente ist ein Sonderfall mit zwei gleichen Schittpunkten)
!   fuer den Fall, dass kein Schnittpunkt existiert, wird t1u=t1v=t2u=t2v=1.e32 gesetzt
!   fuer t1 aus [0,1] liegt der Schnittpunkt auf dem Kreisbogen A-B
!   fuer t2 aus [0,1] liegt der Schnittpunkt auf dem Kreisbogen C-D
!   Die Koordinaten der beiden Schnittpunkte sind: (xpu,ypu) und (xpv,ypv)
!
!   err    =  0   if no erroro occured
!            10   an error snkrgr the result has no meaning
!
!  local variables
!
      real (DP) xe,ye,xf,yf, r1, r2, rlu, rlv
      real (DP) dx12,dy12,d,dxa1,dya1,dxb1,dyb1,dxc2,dyc2,dxd2,dyd2
      real (DP) dx1m,dx2m,dy1m,dy2m,gamma1,gamma2,alpha,xnenn,zaehl
      real (DP) r1h2,r1h2k,r2h2,r2h2k,w1anf,w2anf,w1end,w2end
      real (DP) pi, zweipi
      parameter (pi=3.14159265358979323846264338327950288419716_DP)
      parameter (zweipi=2._DP*pi)
!
      dx12=xm1-xm2
      dy12=ym1-ym2
! Abstand der beiden Kreismittelpunkte
      d=sqrt(dx12*dx12+dy12*dy12)
!
      dxa1=xa-xm1
      dya1=ya-ym1
!  Radius des Kreises 1
      r1h2=dxa1*dxa1+dya1*dya1
!  nochmals Radius zur Kontrolle
      dxb1=xb-xm1
      dyb1=yb-ym1
      r1h2k=dxb1*dxb1+dyb1*dyb1
!  Mitteln der beiden Resultate fuer Anfangs- bzw. Endpunkt
      r1h2=(r1h2+r1h2k)*0.5_DP
      r1=sqrt(r1h2)
!
      dxc2=xc-xm2
      dyc2=yc-ym2
!  Radius des Kreises 2
      r2h2=dxc2*dxc2+dyc2*dyc2
!  nochmals Radius zur Kontrolle
      dxd2=xd-xm2
      dyd2=yd-ym2
      r2h2k=dxd2*dxd2+dyd2*dyd2
!  Mitteln der beiden Resultate fuer Anfangs- bzw. Endpunkt
      r2h2=(r2h2+r2h2k)*0.5_DP
      r2=sqrt(r2h2)
!
      if ( (r1+r2 .le. d ) .or. (abs(r2-r1) .gt. d ) ) then
!  kein Schittpunkt moeglich
        t1u=huge(1._DP)
        t1v=huge(1._DP)
        t2u=huge(1._DP)
        t2v=huge(1._DP)
        err=0
        return
      end if
!
      if ((d.eq.0._DP) .and. (r1.eq.r2) )then
!  Sonderfall: zwei gleiche Kreismittelpunkte und Radien
        w1anf=wink(xm1,ym1,r1,0.d0,xa,ya)
        w1end=wink(xm1,ym1,r1,0.d0,xb,yb)
        w2anf=wink(xm2,ym2,r2,0.d0,xc,yc)
        w2end=wink(xm2,ym2,r2,0.d0,xd,yd)
        if ( ((w1anf.gt.w2anf) .and. (w1anf.lt.w2end)) .or.             &
     &       ((w2end.gt.w1anf) .and. (w2end.lt.w1end)) ) then
          t1u=.5_DP
          t2u=.5_DP
          t1v=.5_DP
          t2v=.5_DP
          xpu=xm1+r1*cos((w1anf+w2end)/2._DP)
          ypu=ym1+r1*sin((w1anf+w2end)/2._DP)
          xpv=xpu
          ypv=ypu
        else if ( ((w1end.gt.w2anf) .and. (w1end.lt.w2end)) .or.          &
     &            ((w2anf.gt.w1anf) .and. (w2anf.lt.w1end)) ) then
          t1u=.5_DP
          t2u=.5_DP
          t1v=.5_DP
          t2v=.5_DP
          xpu=xm1+r1*cos((w1end+w2anf)/2._DP)
          ypu=ym1+r1*sin((w1end+w2anf)/2._DP)
          xpv=xpu
          ypv=ypu
        else
          t1u=huge(1._DP)
          t1v=huge(1._DP)
          t2u=huge(1._DP)
          t2v=huge(1._DP)
        end if
        err=0
        return
      end if
!
!  Festlegung einer Geraden durch die Punkte  E  und  F
!
      if ( dx12 .eq. 0._DP) then
        xe=xm1-r1
        xf=xm1+r1
        ye=( (r2h2-r1h2)/dy12 + ym1+ym2 )*.5_DP
        yf=ye
!
      else if ( dy12 .eq. 0._DP) then
        ye=ym1-r1
        yf=ym1+r1
        xe=( (r2h2-r1h2)/dx12 + xm1+xm2 )*.5_DP
        xf=xe
      else
!
        ye=ym1-r1
        yf=ym1+r1
        xe=((r2h2-r1h2)+(xm1**2-xm2**2)+(ym1**2-ym2**2))/(2._DP*dx12)   &
     &    -(ym1-ym2)/dx12*ye
        xf=((r2h2-r1h2)+(xm1**2-xm2**2)+(ym1**2-ym2**2))/(2._DP*dx12)   &
     &    -(ym1-ym2)/dx12*yf
      end if
!  Schnittpunkte auf den Kreisbogen berechnen
      call snkrgr(xe,ye,xf,yf,xa,ya,xb,yb,xm1,ym1,rlu,rlv,t1u,t1v,err)
      if (err .gt. 0) return
      xpu=xe+rlu*(xf-xe)
      ypu=ye+rlu*(yf-ye)
      xpv=xe+rlv*(xf-xe)
      ypv=ye+rlv*(yf-ye)
!  Winkel zwischen Punkt D und C  vom Mittelpunkt aus gesehen
      zaehl=dxc2*dyd2-dxd2*dyc2
      xnenn=dxd2*dxc2+dyd2*dyc2
      alpha=atan2(zaehl,xnenn)
      if (alpha.le.0.d0) alpha=alpha+zweipi
!  Relative Schnittpunktskoordinaten 1. Schnittpunkt
      dx1m=xpu-xm2
      dy1m=ypu-ym2
!  Winkel zwischen dem Schnittpunkt 1 und dem Punkt C berechnen
      zaehl=dxc2*dy1m-dx1m*dyc2
      xnenn=dx1m*dxc2+dy1m*dyc2
      gamma1=atan2(zaehl,xnenn)
      if (gamma1.le.0._DP) gamma1=gamma1+zweipi
!  Relative Schnittpunktskoordinaten 2. Schnittpunkt
      dx2m=xpv-xm2
      dy2m=ypv-ym2
!  Winkel zwischen dem 2. Schnittpunkt und dem Punkt C berechnen
      zaehl=dxc2*dy2m-dx2m*dyc2
      xnenn=dx2m*dxc2+dy2m*dyc2
      gamma2=atan2(zaehl,xnenn)
      if (gamma2.le.0._DP) gamma2=gamma2+zweipi
      t2u=gamma1/alpha
      t2v=gamma2/alpha
      return
      end subroutine snkrkr