      subroutine readnetin(gbz,gzz,gkz,zki,kzrb,zrb,xbk,ybk,matzif,ok,        &
     &                     alrb,btrb,nnat)
      use feminterface, only: getsetting, readmat, reallocate
      use femtypes
      implicit none
      integer (I4B) gbz, gzz, gkz, nnat
      integer (I4B), pointer ::  kzrb(:,:), zrb(:,:), matzif(:),  zki(:,:)
      real (DP), pointer ::  xbk(:), ybk(:)
      complex (DPC), pointer :: alrb(:,:), btrb(:,:)
      logical ok
      intent (in ) :: nnat
      intent (out) :: gbz, gzz, gkz, zki, kzrb, zrb, xbk, ybk, matzif
      intent (out) :: ok, alrb, btrb
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
!    $Revision: 1.8 $
!    $Date: 2010/08/27 13:00:04 $
!    $Author: m_kasper $
!
!  Einlesen der Topologieinformation, Gebiete, Zweige, Knoten und der Randbedingungen
!
!  Ausgabe:
!     gbz         Gebiets-Bereiche-Zahl (Anzahl der Gebiete)
!     gzz         Gebiets-Zweige-Zahl (Anzahl der Zweige)
!     gkz         Gebiets-Knoten-Zahl (Anzahl der Knoten)
!     bzi         Bereichs-Zweige-Informationen (Liste der Zweige eines Bereichs bzw, Gebietes)
!     bzip        Pointer auf Anfang eines Gebiets in bzi-Liste (kompakte Speicherung)
!     bzil        Pointer auf bzip Anfang eines Gebietes im Vektor einschliesslich aller darin enhaltenen Gebiete
!     kzi         Knoten Zweig Information
!                      /   = 0  : innerer Knoten
!                 kzi  -   < 0  : ist mit Gebietsknoten Nr.(-kzi) identisch
!                      \   > 0  : liegt auf Zweig der Nr. (kzi)
!     kzrb        Zweige deren Randbedingung fuer den Knoten gueltig sind
!     zrb         Randbedingungen der Zweige
!                     0 -  99: Dirichlet
!                   200 - 299: Neumannn
!                   300 - 399: Kontur
!                   400 - 499: innen
!     xbk,ybk     Koordinaten der Bereichs-Knoten
!     matzif      Materialkennziffer fuer jedes Gebiet
!     ok          =.false. wenn ein Fehler aufgetreten ist
!     alrb, btrb  (komplexe) Parameter der allgemeinen Randbedingungen
!
!  lokale variablen
      real (DP) zpp
      integer (I4B) i,j,zpz,zahl,point1,point2,unitid, ios, alstat
      integer (I4B) , pointer :: bzi(:), bzip(:), bzil(:)
      logical innen, found
      character (len=50), pointer :: matname(:)
      character (len=200) path
!
      external    ::    grglun
!
      ok=.true.
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
      open (unitid,file=path(1:len_trim(path))//'netinf',               &
     &  form='unformatted',status='old',action='READ',iostat=ios)
      if (ios .ne. 0) then
        print*, 'Error on opening file netinf'
        print*, 'IO Error No: ',ios
        ok=.false.
        return
      end if
      rewind unitid
!
!  prepare reading key-points
      read (unitid,iostat=ios) gkz
      if (ios .ne. 0) ok=.false.
      write (*,100) gkz
100   format('  key-points:',i6)
      allocate(kzrb(gkz,nnat), xbk(gkz), ybk(gkz),stat=alstat)
      if (alstat .ne. 0) then
        ok=.false.
        print*,'Error in array allocation'
        return 
      end if
!
!  read the keypoints
      do i=1,gkz
        read (unitid) xbk(i),ybk(i),(kzrb(i,j),j=1,nnat)
      end do
!
!  prepare reading branches
      read (unitid,iostat=ios) gzz
      if (ios .ne. 0) ok=.false.
      write (*,101) gzz
101   format('  branches:  ',i6)
      allocate(zki(3,gzz),zrb(gzz,nnat),alrb(gzz,nnat),btrb(gzz,nnat),bzi(2*gzz),      &
     &  stat=alstat)
      if (alstat .ne. 0) then
        ok=.false.
        print*,'Error in array allocation'
        return 
      end if
!
!  read the branches
      do i=1,gzz
        read (unitid) zki(:,i),zpz,zpp,(zrb(i,j),alrb(i,j),btrb(i,j),j=1,nnat)
      end do
!     
!  prepare reading regions
      read (unitid,iostat=ios) gbz
      if (ios .ne. 0) ok=.false.
      write (*,102) gbz
102   format('  regions:   ',i6)
      allocate(bzil(gbz+1),bzip(2*gbz+1),matzif(gbz),matname(gbz),      &
     &  stat=alstat)
      if (alstat .ne. 0) then
        ok=.false.
        print*,'Error in array allocation'
        return 
      end if
!
!  read the regions
      point1=1
      point2=1
      do i=1,gbz
        bzil(i)=point1
        innen=.true.
        do while (innen)
          if (point1.gt.size(bzip)) then
            bzip=>reallocate(bzip,2*size(bzip))
          end if
          bzip(point1)=point2
          read (unitid,iostat=ios) zahl,(bzi(point2-1+j),j=1,zahl),matname(i)
          point1=point1+1
          point2=point2+zahl
          if (matname(i).ne.'-1') innen=.false.
        end do
        call readmat(matname(i),found,matzif(i))
      end do
      if (ios .ne. 0) ok=.false.
      bzil(gbz+1)=point1
      if (point1.gt.size(bzip)) then
        bzip=>reallocate(bzip,size(bzip)+1)
      end if
      bzip(point1)=point2
!
      close (unitid)
      deallocate (matname, bzi, bzip, bzil)
      return
      end subroutine readnetin
! 
