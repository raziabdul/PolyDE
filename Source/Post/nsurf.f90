      subroutine nsurf(z,x,y,p,e,ie,palette,ncolor,zcolmin,zcolmax,     &
     &     zscale,color,hlin,cross,delcross,theta,phi)
      use feminterface, only : transhom, scalehom, rotatehom, fsortn,   &
     &  dptransformpoly, drtri, hline, hlinexy
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:), zcolmin, zcolmax, zscale, delcross 
      real (DP) phi, theta
      integer (I4B) p, e(:,:), ie, ncolor, palette(:)
      logical color, hlin, cross
      intent(in):: z, x, y, p, e, ie, palette, ncolor, zcolmin,         &
     &   zcolmax, zscale, color, hlin, cross, delcross, theta, phi
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
!  Falschfarbendarstellung ueber einem Dreiecksnetz  3D Plot
!  Faebung entsprechend den Hoehenwerte z in den Dreickseckpunkten
!  Interpolation der Farben auf dem Dreieck
!  Projektion zur Darstellung der Oberflaechentopologie
!
!     z        Hoehenwerte in den Dreieckseckpunkte
!     x        x Koordinaten der Dreicksecken
!     y        y Koordinaten der Dreicksecken
!     p        Anzahl der Netzknoten
!     e        Elementinformation Knoten der Dreiecke
!     ie       Anzahl der Dreiecke
!     palette  Farbwerte der in aufsteigender Reihenfolge
!     ncolor   Anzahl der Farbwerte
!     zcolmin  Hoehenwert zur niedrigsten Farbe
!     zcolmax  Hoehenwert zur hoechsten Farbe
!     zscale   Skalierungsfaktor der Hoehenwerte fuer die perspektivische Darstellung
!              zscale > 1 fuerht zu Ueberhoehung
!     color    bei .true. wird entspechend der Farbpalette eingefaerbt
!     hlin     bei .true. werden die Hoehenlinien mit der Abstufung aus zcolor gezeichnet
!     cross    bei .true. Gitterdarstellung, Linien konstanter in x und y Werte
!     delcross Abstand der Gitterlinien in x und y Richtung
!     theta    Winkel zwischen der z-Achse und dem Beobachter in Winkelgrad
!               =0 Grad Aufsicht
!     phi      Drehwinkel der x-y Ebene um die z-Achse

      real (DP) xkmin, xkmax, ykmin, ykmax, zkmin, zkmax, transf(4,4)
      real (DP) xo(p), yo(p), zo(p)
      real (DP), allocatable :: zcolor(:), zsort(:)
      real (SP) ocr, ocg, ocb, xpoly(3), ypoly(3)
      integer (I4B) i, oldindex, ii, iminx(1), imaxx(1), iminy(1)
      integer (I4B) imaxy(1), iminz(1), imaxz(1), liste(ie)
      logical proj
!
      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgsfs, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgqcs, pgqwin, pgqvp, pgsvp, pgswin, pgwnad
      external :: pgrect, pgptxt, pgpoly

      iminx=minloc(x)
      imaxx=maxloc(x)
      iminy=minloc(y)
      imaxy=maxloc(y)
      iminz=minloc(z)
      imaxz=maxloc(z)
      xkmin=x(iminx(1))
      xkmax=x(imaxx(1))
      ykmin=x(iminy(1))
      ykmax=x(imaxy(1))
      zkmin=x(iminz(1))
      zkmax=x(imaxz(1))
!  Transformationsmatrix initialisieren
      transf=0._DP
      transf(1,1)=1._DP
      transf(2,2)=1._DP
      transf(3,3)=1._DP
      transf(4,4)=1._DP
!  Translation in den Ursprung
!      call transhom(transf,-(xkmin+xkmax)/2._DP,-(ykmin+ykmax)/2._DP,-(zkmin+zkmax)/2._DP)
!  Skalieren der z-Skala
      call scalehom(transf,1._DP,1._DP,zscale)
!  Rotation um die z-Achse
!      call rotatehom(transf,3.141596_DP/4._DP,3)
      call rotatehom(transf,phi * 3.141596_DP/180._DP,3)
!  Kippen
!      call rotatehom(transf,3.141596_DP/1.6_DP,2)
      call rotatehom(transf,-theta * 3.141596_DP/180._DP,1)
!
!  Durchfuehren der Transformation
      call dptransformpoly(x,y,z,transf,p,xo,yo,zo)
      proj=.true.

!  Erstelle Indexliste fuer depth sort
      allocate (zsort(ie))
      zsort= zo(e(1,:))+zo(e(2,:))+zo(e(3,:))
      call fsortn(zsort,liste,ie)
      deallocate (zsort)
!  Farbstatus retten
      call pgqci(oldindex)
      call pgsci(20)
      call pgqcr(20, ocr, ocg, ocb)
!  Set Filling style
      call pgsfs(1)

!  Zuordnung der Farbwerte zu den Hoehenwerten
      allocate(zcolor(ncolor))
      do i=1,ncolor
        zcolor(i)=(zcolmax-zcolmin)/dble(ncolor)*dble(i)+zcolmin
      end do
!  Schleife ueber alle Elemente

      do ii=1,ie
        i=liste(ie+1-ii)
!  Falschfarbendarstellung
        if (color) then
          call drtri(i,e,x,y,z,zcolor,ncolor,palette,proj,xo,yo,zo,transf)
        else 
!  Weiss ausmalen
          if (proj) then
            xpoly(:)=xo(e(:,i))
            ypoly(:)=yo(e(:,i))
          else
            xpoly(:)=x(e(:,i))
            ypoly(:)=y(e(:,i))
          end if
          call pgscr(20, 1._SP, 1._SP, 1._SP)
          call pgpoly(3,xpoly,ypoly)
        end if
!  Hoehenlinien
        if (hlin) then
          call pgsci(1)
          call hline(i,e,x,y,z,zcolor,ncolor,proj,transf)
          call pgsci(20)
        end if
!  Gitternetz
        if (cross) then
          call pgsci(1)
          call hlinexy(i,e,x,y,z,0._DP,delcross,proj,1,transf)
          call hlinexy(i,e,x,y,z,0._DP,delcross,proj,2,transf)
          call pgsci(20)
        end if
      end do
!
      deallocate (zcolor)
!  Farbstatus zuruecksetzen
      call pgsci(oldindex)
      call pgscr(20, ocr, ocg, ocb)
      return
      end subroutine nsurf
