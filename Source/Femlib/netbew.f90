      subroutine netbew(xn,yn,e,rk1,n)
      use feminterface, only: inkrp, zeit, angles
      use femtypes
      implicit none
      real (DP) xn(:),yn(:), rk1
      integer (I4B) e(:,:), n
      intent (in) :: xn, yn, e,  n
      intent (out) :: rk1
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
!    $Date: 2006/07/07 21:08:36 $
!    $Author: r_abdul $
!
!  Programm zum Bewerten der Triangulierung
!  Es werden alle Elemente untersucht und die Abweichungen der Winkel
!  von 60 Grad Winkeln (euklidisch) aufsummiert
!
!  Eingabe:
!     xn,yn     Koordinaten der Knoten (Anfangswerte)
!     e       Elementinformation
!     n       Anzahl der Elemente
!  Ausgabe:  
!     rk1     mittlere quadratisch Abweichung von 60 Grad-Winkeln
!
!  lokale Variablen
!
      integer (I4B) i
      real (DP) pi,rk2,wmin,wmax,w1,w2,w3,xmp,ymp
      real (DP) xra,yra,rausen,flaech,rho,pid3
      parameter (pi=3.14159265358979323846264338327950288419716_DP)
      parameter (pid3 = pi/3._DP)
!
!      write (*,1)
!1     format('  Netzbewertung')
      rk1=0._DP
      rk2=0._DP
      wmin=7._DP
      wmax=0._DP
!  fuer alle Elemente      I
      do i=1,n
        call angles(i,e,xn,yn,w1,w2,w3)
        wmax=max(wmax,w1,w2,w3)
        wmin=min(wmin,w1,w2,w3)
        w1=w1-pid3
        w2=w2-pid3
        w3=w3-pid3
!  Abweichung vom Standardwinkel (60 Grad) berechnen
        rk1=rk1+w1*w1+w2*w2+w3*w3
!  Innkreisradius rho berechnen
        call inkrp(xn(e(1,i)),yn(e(1,i)),xn(e(2,i)),yn(e(2,i)),         &
     &    xn(e(3,i)),yn(e(3,i)),xmp,ymp,flaech,rho,xra,yra,rausen)
!        rk2=max(rk2,sqrt(abs(flaech))/rho)
        if (rho.ne.0._DP) then
          rk2=max(rk2,rausen/rho)
        else
          rk2=huge(1._DP)
        endif
      end do
      write(*,1555) wmax/pi*180._DP
1555  format('   Maximum Angle: ',g12.6)
      write(*,1556) wmin/pi*180._DP
1556  format('   Minimum Angle: ',g12.6)
!
!   rk1 = L2 Norm (winkel  -  60 Grad)      ueber alle Winkel
!
!  die Zahl rk1 sollte fuer gute Netze moeglichst klein sein
!
      rk1=sqrt(rk1)/dble(3*n)
!
!  rk2= maximales Verhaeltnis von Umkreisradius/Inkreisradius
!
!  Die Zahl rk2 sollte fuer gute Netze moeglichst nahe bei 1. liegen
!
      rk2=rk2/2._DP
      write (*,987) rk1,rk2
987   format('   rk1: ',g12.6,'       rk2: ',1g12.6)
!      call zeit(' Mesh Evaluation')
      return
      end subroutine netbew