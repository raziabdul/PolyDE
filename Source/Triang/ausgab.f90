      subroutine ausgab(ele,xk,yk,nb,ez,kz,kzi,zrb,krb,geb,konfom,zyl)
      use feminterface, only: basout !, lout
      use femtypes
      implicit none
      integer (I4B) ele(:,:),nb(:,:),ez,kz,kzi(:),zrb(:,:),krb(:,:), geb(:)
      real (DP) xk(:),yk(:)
      logical konfom,zyl
      intent (in) :: nb, ez, kz, zrb, krb, geb, konfom, zyl
      intent (inout) :: ele, xk, yk, kzi
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
!    $Revision: 1.9 $
!    $Date: 2008/08/20 13:02:13 $
!    $Author: m_kasper $
!
!  kti-Knotentabelle indirekt  (Hilfsfeld)
!
!  local variables
      integer (I4B) irb0,irb1,irb2,irb3,irb4
      integer (I4B) irb0z,irb1z,irb2z,irb3z,irb4z,i0z
      integer (I4B) i,i99,kzneu,j, kti(2,kz), kz1(kz)
      real (DP) xtmp(kz),ytmp(kz)
      logical ok
!
!  Ermittele, wieviele Knoten welche Randbedingung haben
      i0z=0
      irb0z=0
      irb1z=0
      irb2z=0
      irb3z=0
      irb4z=0
      i99=0
      do i=1,kz
        if (kzi(i).eq.0) then
          irb4z=irb4z+1
        else
          if (kzi(i).gt.0) then
            select case (zrb(kzi(i),1))
            case (0:99)
              i0z=i0z+1
            case (100:199)
              irb0z=irb0z+1
            case (200:299)
              irb1z=irb1z+1
            case (300:399)
              irb2z=irb2z+1
            case (400:499)
              irb3z=irb3z+1
            end select
          else
            select case (krb(-kzi(i),1))
            case (0:99)
              i0z=i0z+1
            case (100:199)
              irb0z=irb0z+1
            case (200:299)
              irb1z=irb1z+1
            case (300:399)
              irb2z=irb2z+1
            case (400:499)
              irb3z=irb3z+1
            case default
              i99=i99+1
            end select
          end if
        end if
      end do
!  Bestimme die irb Werte
      irb4=irb4z+1
      irb3=irb4+irb3z
      irb2=irb3+irb2z
      irb1=irb2+irb1z
      irb0=irb1+irb0z
      i0z=1
      irb4z=irb4
      irb3z=irb3
      irb2z=irb2
      irb1z=irb1
      irb0z=irb0
      kzneu=kz-i99
!  Sortiere die Knoten nach ihren Randbedingungen.
!  Verwende hierzu die Knoten-Tabelle-indirekt
!
!  Schreibe fuer jeden Knoten die neue Position
!  in die Tabelle kti(1,knotennummer)
!
!  Vermerke fuer die neue Position die alte Nummer in der
!  Tabelle kti(2,neue_nummer)
!
      do i=1,kz
        if (kzi(i).eq.0) then
          kti(1,i)=i0z
          kti(2,i0z)=i
          i0z=i0z+1
        else
          if (kzi(i).gt.0) then
            select case (zrb(kzi(i),1))
            case (0:99)
              kti(1,i)=irb0z
              kti(2,irb0z)=i
              irb0z=irb0z+1
            case (100:199)
              kti(1,i)=irb1z
              kti(2,irb1z)=i
              irb1z=irb1z+1
            case (200:299)
              kti(1,i)=irb2z
              kti(2,irb2z)=i
              irb2z=irb2z+1
            case (300:399)
              kti(1,i)=irb3z
              kti(2,irb3z)=i
              irb3z=irb3z+1
            case (400:499)
              kti(1,i)=irb4z
              kti(2,irb4z)=i
              irb4z=irb4z+1
            end select
          else
            select case (krb(-kzi(i),1))
            case (0:99)
              kti(1,i)=irb0z
              kti(2,irb0z)=i
              irb0z=irb0z+1
            case (100:199)
              kti(1,i)=irb1z
              kti(2,irb1z)=i
              irb1z=irb1z+1
            case (200:299)
              kti(1,i)=irb2z
              kti(2,irb2z)=i
              irb2z=irb2z+1
            case (300:399)
              kti(1,i)=irb3z
              kti(2,irb3z)=i
              irb3z=irb3z+1
            case (400:499)
              kti(1,i)=irb4z
              kti(2,irb4z)=i
              irb4z=irb4z+1
            end select
          end if
        end if
      end do
      do i=1,ez
        do j=1,3
!  Ersetze alte Knotennummer durch neue Nummer
          ele(j,i)=kti(1,ele(j,i))
        end do
      end do
      do i=1,kzneu
!  Setze die alten Koordinaten an die neue Listenposition
!  Speichere die Daten in den Tabellen (x,y,kz1) zwischen
        xtmp(i)=xk(kti(2,i))
        ytmp(i)=yk(kti(2,i))
        kz1(i)=kzi(kti(2,i))
      end do
      do i=1,kzneu
!  Verwende wieder die richtigen Tabellen
        xk(i)=xtmp(i)
        yk(i)=ytmp(i)
        kzi(i)=kz1(i)
      end do
!  Gebe das Netz aus.
      call basout(ok,ez,kzneu,ele,xk,yk,kzi,geb,nb)
!  Erzeuge einen beliebigen, ungueltigen Loesungsvektor,
!  Damit beim Einlesen keine Fehler auftreten.
!      call lout(.false., .false., 0, 0, 0._DP, 0._DP, 0._DP)
      return
      end subroutine ausgab
