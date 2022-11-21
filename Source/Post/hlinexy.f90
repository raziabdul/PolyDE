      subroutine hlinexy(i,e,x,y,z,z0,delta,proj,direction,transf)
      use feminterface, only : sptransformpoly
      use femtypes
      implicit none
      real (DP) x(:), y(:), z(:), z0, delta
      real (DP), optional ::  transf(4,4)
      integer (I4B) i, e(:,:), direction
      logical proj
      intent(in):: i, e, x, y, z, z0, delta, transf, proj
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
!     direction 1=x 2=y 3=z
!     transf   Tranformationsmatrix (Homogene Koordinaten)
!
!  Zeichnen von Hoehenlinien auf einem Dreieck
!
      real (SP) xx(2), yy(2), zz(2), xout(2), yout(2), zout(2)
      real (DP) xl(3), yl(3), zl(3), zmin, zmax, zval, dz(3), dd
      integer (I4B) iz
!
!  lokale Kopie der Eckkoordinaten und Hoehenwerte
      xl=x(e(:,i))
      yl=y(e(:,i))
      zl=z(e(:,i))
!
      select case (direction)
      case (1)
        zmax=maxval(xl)
        zmin=minval(xl)
      case (2)
        zmax=maxval(yl)
        zmin=minval(yl)
      case (3)
        zmax=maxval(zl)
        zmin=minval(zl)
      end select

!  Erster zu verwendender Hoehenwert feststellen
       zval=floor((zmin-z0)/delta) * delta
!  Schleife ueber die Hoehenwerte
      do 
        zval=zval+delta
        if(zval.gt.zmax) exit
        iz=0
        select case (direction)
        case (1)
          dz(:)=xl(:)-zval
        case (2)
          dz(:)=yl(:)-zval
        case (3)
          dz(:)=zl(:)-zval
        end select
!  Schnittpunkte ausrechnen
        if( dz(1)*dz(2) .lt. 0._DP) then
          iz=iz+1
          dd=dz(1)/(dz(2)-dz(1))
          xx(iz)=sngl(xl(1)+(xl(1)-xl(2))*dd)
          yy(iz)=sngl(yl(1)+(yl(1)-yl(2))*dd)
          zz(iz)=sngl(zl(1)+(zl(1)-zl(2))*dd)
        end if
        if( dz(2)*dz(3) .lt. 0._DP) then                    
          iz=iz+1
          dd=dz(2)/(dz(3)-dz(2))
          xx(iz)=sngl(xl(2)+(xl(2)-xl(3))*dd)
          yy(iz)=sngl(yl(2)+(yl(2)-yl(3))*dd)
          zz(iz)=sngl(zl(2)+(zl(2)-zl(3))*dd)
        end if
        if (iz .eq. 0) cycle
        if (iz .ne. 2) then
          iz=iz+1
          dd=dz(3)/(dz(1)-dz(3))
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
      end subroutine hlinexy
