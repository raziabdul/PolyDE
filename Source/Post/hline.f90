      subroutine hline(i,e,x,y,z,zcolor,ncolor,proj,transf)
      use feminterface, only : locate, sptransformpoly
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:) ,zcolor(:)
      real (DP), optional ::  transf(4,4)
      integer (I4B) i, e(:,:), ncolor
      logical proj
      intent(in):: i, e, x, y, z, zcolor, ncolor, transf, proj
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
!  Hoehenliniendarstellung eines Dreiecks
!  und optional Projektion zur Darstellung der Oberflaechentopologie
!
!     i        Nummer des Dreiecks
!     e        Elementinformation Knoten der Dreiecke
!     x        x Koordinaten der Dreicksecken
!     y        y Koordinaten der Dreicksecken
!     z        Hoehenwerte in den Dreieckseckpunkte
!     zcolor   Hoehenwerte der Stufen 
!     ncolor   Anzahl der Stufen
!     proj     .true. fuer einen 3D Plot 
!              .false. fuer einen 2D Plot
!     transf   Tranformationsmatrix (Homogene Koordinaten)
!
!  Zeichnen von Hoehenlinien auf einem Dreieck
!
      real (DP) xl(3), yl(3), zl(3), zmin, zmax, zval, dz(3), dd
      real (SP) xx(2), yy(2), zz(2), xout(2), yout(2), zout(2)
      integer (I4B) j, icolor, iz
!
      external :: pgqci, pgqcr, pgsfs, pgsci, pgscr, pgpoly, pgmove, pgdraw

!  lokale Kopie der Eckkoordinaten und Hoehenwerte
      xl=x(e(:,i))
      yl=y(e(:,i))
      zl=z(e(:,i))
!
      zmax=maxval(zl)
      zmin=minval(zl)

!  Erster zu verwendender Hoehenwert feststellen
      icolor=locate(zcolor,zmin)
      if (icolor.ge.ncolor) return
      if (icolor.lt.1) icolor=1
        
!  Schleife ueber die Hoehenwerte
      do j=icolor,ncolor
        zval=zcolor(j)
        if(zval.gt.zmax) exit
        iz=0
        dz(:)=zl(:)-zval
!  Schnittpunkte ausrechnen
        if( dz(1)*dz(2) .lt. 0._DP) then
          iz=iz+1
          dd=dz(1)/(zl(2)-zl(1))
          xx(iz)=sngl(xl(1)+(xl(1)-xl(2))*dd)
          yy(iz)=sngl(yl(1)+(yl(1)-yl(2))*dd)
          zz(iz)=sngl(zl(1)+(zl(1)-zl(2))*dd)
        end if
        if( dz(2)*dz(3) .lt. 0._DP) then                    
          iz=iz+1
          dd=dz(2)/(zl(3)-zl(2))
          xx(iz)=sngl(xl(2)+(xl(2)-xl(3))*dd)
          yy(iz)=sngl(yl(2)+(yl(2)-yl(3))*dd)
          zz(iz)=sngl(zl(2)+(zl(2)-zl(3))*dd)
        end if
        if (iz .eq. 0) cycle
        if (iz .ne. 2) then
          iz=iz+1
          dd=dz(3)/(zl(1)-zl(3))
          xx(iz)=sngl(xl(3)+(xl(3)-xl(1))*dd)
          yy(iz)=sngl(yl(3)+(yl(3)-yl(1))*dd)
          zz(iz)=sngl(zl(3)+(zl(3)-zl(1))*dd)
        end if
        if (proj) then
!  Transformation der Hoehenlinie
          call sptransformpoly(xx,yy,zz,transf,2,xout,yout,zout)
!  Positionieren und Zeichnen
          call pgmove(xout(1),yout(1))
          call pgdraw(xout(2),yout(2))
        else
!  Positionieren und Zeichnen
          call pgmove(xx(1),yy(1))
          call pgdraw(xx(2),yy(2))
        end if
      end do
      return
      end subroutine hline
