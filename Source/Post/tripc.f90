      subroutine tripc(zr,zc,nc,x,y,e,en,ie)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) x(:), y(:), zr(:), zc(:)
      integer (I4B) e(:,:), en(:,:), ie, nc
      intent(in):: zr, zc, nc, x, y, e, en, ie
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
!    $Date: 2006/07/07 22:13:33 $
!    $Author: r_abdul $
!
!  Zeichnen von Hoehenlinien ueber einer Triangulierung
!     zr    Vektor der Hoehenwerte in den Knoten
!     zc    Hoehenerte der zu zeichnenden Hoehenlinien (aufsteigend)
!     nc    Anzahl der Hoehenlinien
!     x     x-Koordinaten der Knoten
!     y     y-Koordinaten der Knoten
!  :  e     Elementinformation
!     en    Nachbarn der Elemente
!     ie    Anzahl der Elemente
      real (DP) h, x1, x2, x3, y1, y2, y3, a1, a2, a3, dh1, dh2, dh3
      real (DP) dd, amin
      real (SP) xe, ye, xx(2), yy(2), xsc, ysc
      integer (I4B) i, j, k, iz, izwg, ifound, iscnd, inext, inxscn
      integer (I4B) klo, khi, iedge(2), done(ie)
      logical second
!
      external :: pgdraw, pgmove
!  Inizialisieren des Bearbeitungs-Feldes
      do i=1,ie
!  done   wird auf den entsprechenden Zaehler der Hoehenlinie
!         gesetzt, wenn das Element bereits abgearbeitet ist
        a1=zr(e(1,i))
        a2=zr(e(2,i))
        a3=zr(e(3,i))
        amin=min(a1,a2,a3)
!  suchen des kleinsten im Element zu beruecksichtigenden Hoehenwertes
        if (amin.le.zc(1)) then
          done(i)=0
          cycle
        end if
        if (amin.gt.zc(nc)) then
          done(i)=nc
          cycle
        end if
        klo=1
        khi=nc
10      continue
        if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if (zc(k).gt.amin) then
            khi=k
          else
            klo=k
          end if
          goto 10
        end if
        if(zc(khi).lt.max(a1,a2,a3)) then
          done(i)=klo
        else
!  h > amax in diesem Element kann keine Hoehenlinie liegen
          done(i)=nc
        end if
      end do
!
!  Schleife ueber die Hoehenwerte
      do j=1,nc
!  h   aktueller Hoehenwert
        h=zc(j)
        ifound=0
250     continue
!  Suchen eines Elementes, das den Hoehenwert enthaelt
!  Schleife ueber die Elemente
        ifound=0
        do i=ifound+1,ie
          if (done(i).lt.j) then
            a1=zr(e(1,i))
            a2=zr(e(2,i))
            a3=zr(e(3,i))
            if (h .lt. max(a1,a2,a3)) then
              ifound=i
              exit
            else
!  h > amax in diesem Element kann keine weitere Hoehenlinie liegen
              done(i)=nc
            end if
          end if
        end do
!  kein weiteres Element mit diesem Hoehenwert gefunden => naechste Hoehenlinie
        if (ifound .eq. 0) cycle
!  Verfolgen der Hoehenlinie
        i=ifound
        x1=x(e(1,i))
        x2=x(e(2,i))
        x3=x(e(3,i))
        y1=y(e(1,i))
        y2=y(e(2,i))
        y3=y(e(3,i))
        iz=0
!  Schnittpunkte ausrechnen
        if((a1.lt.h).and.(a2.gt.h) .or. (a1.gt.h).and.(a2.lt.h)) then
          iz=iz+1
          iedge(iz)=3
          dh1=a1-h
          dd=dh1/(a2-a1)
          xx(iz)=real(x1+(x1-x2)*dd)
          yy(iz)=real(y1+(y1-y2)*dd)
        end if
        if((a2.lt.h).and.(a3.gt.h) .or. (a2.gt.h).and.(a3.lt.h)) then
          iz=iz+1
          iedge(iz)=1
          dh2=a2-h
          dd=dh2/(a3-a2)
          xx(iz)=real(x2+(x2-x3)*dd)
          yy(iz)=real(y2+(y2-y3)*dd)
        end if
!
        if(iz.ne.2) then
!
          if((a1.lt.h).and.(a3.gt.h) .or. (a1.gt.h).and.(a3.lt.h)) then
            iz=iz+1
            iedge(iz)=2
            dh3=a3-h
            dd=dh3/(a1-a3)
            xx(iz)=real(x3+(x3-x1)*dd)
            yy(iz)=real(y3+(y3-y1)*dd)
          end if
!  Fall keine zwei Scnittpunkte des Hoehenwertes mit den Dreiecksseiten gefunden wurden
          if (iz.ne.2) then
            print*,'*** Fehler in Tripc (Fliesspunktgenauigkeit unzureichend)'
            done(i)=j
            goto 250
          end if
        end if
!  wenn moeglich am Rand beginnen
        if (en(iedge(1),i) .le. 0) then
          inext=en(iedge(2),i)
!  positionieren und zeichnen
          call pgmove(xx(1),yy(1))
          call pgdraw(xx(2),yy(2))
          done(i)=j
          second=.false.
        else
          inext=en(iedge(1),i)
!  positionieren und zeichnen
          call pgmove(xx(2),yy(2))
          call pgdraw(xx(1),yy(1))
          done(i)=j
          if (en(iedge(2),i) .le. 0) then
            second=.false.
          else
!  Element und naechsten Nachbarn merken
            second=.true.
            iscnd=i
            inxscn=en(iedge(2),i)
            xsc=xx(2)
            ysc=yy(2)
          end if
        end if
!
!  verfolgen der Hoehenlinie im Element inext
500     if (inext .le. 0) goto 700
        if (done(inext).ge.j) goto 700
        do k=1,3
          if (en(k,inext).eq.i) then
            izwg=k
            exit
          end if
        end do
        i=inext
        x1=x(e(1,i))
        y1=y(e(1,i))
        x2=x(e(2,i))
        y2=y(e(2,i))
        x3=x(e(3,i))
        y3=y(e(3,i))
        a1=zr(e(1,i))
        a2=zr(e(2,i))
        a3=zr(e(3,i))
!  Schnittpunkte ausrechnen
        if (izwg.eq.3) goto 610
        if((a1.lt.h).and.(a2.gt.h) .or. (a1.gt.h).and.(a2.lt.h)) then
          dh1=a1-h
          dd=dh1/(a2-a1)
          xe=real(x1+(x1-x2)*dd)
          ye=real(y1+(y1-y2)*dd)
          inext=en(3,I)
          goto 630
        end if
610     if (izwg.eq.1) goto 620
        if((a2.lt.h).and.(a3.gt.h) .or. (a2.gt.h).and.(a3.lt.h)) then
          dh2=a2-h
          dd=dh2/(a3-a2)
          xe=real(x2+(x2-x3)*dd)
          ye=real(y2+(y2-y3)*dd)
          inext=en(1,I)
          goto 630
        end if
620     if (izwg.eq.2) then
          print*,'*** Fehler in Tripc (Fliesspunktgenauigkeit unzureichend)'
          done(i)=j
          goto 700
        end if
        if((a1.lt.h).and.(a3.gt.h) .or. (a1.gt.h).and.(a3.lt.h)) then
          dh3=a3-h
          dd=dh3/(a1-a3)
          xe=real(x3+(x3-x1)*dd)
          ye=real(y3+(y3-y1)*dd)
          inext=en(2,I)
        end if
!  zeichnen ohne Positionieren
630     call pgdraw(xe,ye)
        done(I)=J
        goto 500
!
700     if (second) then
          second=.false.
          inext=inxscn
          i=iscnd
!  nur Positionieren
          call pgmove(xsc,ysc)
          goto 500
        end if
!
        goto 250
      end do
      return
      end subroutine tripc
