      subroutine wnetin(gbz,gzz,gkz,bzi,bzip,bzil,zki,kzrb,zpz,zpp,     &
     &  xbk,ybk,matname,filenm,lzrb,layrb,alrb,btrb,laytxt,layanz,nnat)
      use feminterface, only: getsetting
      use femtypes
      implicit none
      integer (I4B) :: gbz, gzz, gkz, nnat, layanz
      integer (I4B) :: bzi(:), bzip(:), bzil(:), zki(:,:), kzrb(:,:)
      integer (I4B) :: zpz(:), lzrb(:), layrb(:,:)
      real (DP) :: xbk(:), ybk(:), zpp(:)
      character (len=*) :: filenm, matname(:), laytxt(:)
      complex (DPC) :: alrb(:,:), btrb(:,:)
! ORG
      !      intent (in) :: gbz, gzz, gkz, bzi, bzip, bzil, zki, kzrb, zpz, zpp
!!      intent (in) :: xbk, ybk, matname, filenm, lzrb, layrb, alrb, btrb
 !     intent (in) :: laytxt, layanz, nnat
! TEST:
      intent (inout) :: gbz, gzz, gkz, bzi, bzip, bzil, zki, kzrb, zpz, zpp
      intent (inout) :: xbk, ybk, matname, filenm, lzrb, layrb, alrb, btrb
      intent (inout) :: laytxt, layanz, nnat

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
!    $Revision: 1.11 $
!    $Date: 2009/03/23 12:02:00 $
!    $Author: chvokas $
!
!  Schreibe die Problemstellung
!  a) Formatted zum Ergaenzen der Materialien bei automatischer
!               Gebietsgenerierung
!  b) Unformatted zur weiteren Verwendung in Fem
!
!  local variables
      integer (I4B) :: i, j, k
      integer (I4B) :: unitid1, unitid2, ios
      character (len=200) :: path,matz
      character (len=13) :: frm
      character (len=1000) :: fmtstr1, fmtstr2, frmstA
!
      call grglun(unitid1)
!  Format of file filenm (i.e. netin.acd)
!    gkz                                                             number of key-points
!    xbk,ybk                                                         coordinates of key-points
!    ...       for each key-point
!    gzz                                                             number of branches
!    zki(1,.),zki(2,.)zki(3,.),zpz,zpp,laytxt                        start-point, end-point, third point,number of nodes 
!    ...       for each branch                                       .. distribution, name of boundary condition
!    layanz                                                          number of boundary conditions
!    (zrb,alrb,btrb; i=1..nnat), laytxt                              (kind, parameters al, bt) for each nature, name 
!    ...       for each boundary condition
!    gbz                                                             number of regions
!    zwzahl,(bzi(point2-1+j),j=1,zwzahl),matname(i)                  number of branches, branches, material
!    ...       for each region
      open(unitid1,file=filenm,                                         &
     &  form='formatted',position='REWIND',action='WRITE',iostat=ios)
      if (ios .ne. 0) then
        print*,'**** Error opening file: ',filenm
        print*,'**** Error no: ',ios
        return
      end if
      call getsetting('PROJECTPATH',path)
      call grglun(unitid2)
!  Format of file netinf
!    gkz                                                             number of key-points
!    xbk,ybk,kzrb                                                    coordinates of key-points
!    ...       for each key-point
!    gzz                                                             number of branches
!    zki(1,.),zki(2,.)zki(3,.),zpz,zpp,(zrb,alrb,btrb; i=1..nnat)
!    ...       for each branch
!    gbz                                                             number of regions
!    zwzahl,(bzi(point2-1+j),j=1,zwzahl),matname(i)                  number of branches, branches, material
!    ...       for each region
!  boundary conditions are stored for key-points and branches
      open(unitid2,file=path(1:len_trim(path))//'netinf',               &
     &  form='unformatted',position='REWIND',action='WRITE',iostat=ios)
      if (ios .ne. 0) then
        print*,'**** Error opening file: ',path(1:len_trim(path))//'netinf'
        print*,'**** Error no: ',ios
        return
      end if
!
!  write the key-points
      write (unitid1,333) gkz
333   format(i5,46x,'// Key-points')
      write (unitid2) gkz
      do i=1,gkz
        write (unitid1,444) xbk(i),ybk(i),i
444     format(g25.17, ' ', g25.17, '              //', i6)
!
        write (unitid2) xbk(i),ybk(i),(kzrb(i,j),j=1,nnat)
      end do
!
!  write the branches
      write (unitid1,555) gzz
555   format(i5,46x,'// Branches')
      write (unitid2) gzz
      do i=1,gzz
        frmstA='(3i7, i5, f7.3, i6, ''            //'', i6)'
        write (unitid1,frmstA) zki(1,i),zki(2,i),zki(3,i),zpz(i),zpp(i),lzrb(i),i
!
        write (unitid2) zki(:,i),zpz(i),zpp(i),                         &
     &        (layrb(lzrb(i),j),alrb(lzrb(i),j),btrb(lzrb(i),j),j=1,nnat)
      end do
!
!  write the boundary conditions
      write (unitid1,666) layanz
666   format(i5,46x,'// Boundary Conditions')
      do i=1,layanz
        fmtstr1='( '
        fmtstr2=' (i5,''  ('',g13.5,'','',g13.5,'')  '','' &
        & ('',g13.5,'','',g13.5,'')  '',''    ''), a)'
        write(frmstA,'(a,i3,a)') trim(fmtstr1),nnat,trim(fmtstr2)
        write (unitid1,frmstA) (layrb(i,j),alrb(i,j),btrb(i,j),j=1,nnat),trim(laytxt(i))
      end do
!
!  write the regions
      write (unitid1,777) gbz
777   format(i5,46x,'// Regions')
      write (unitid2) gbz
      do k=1,gbz
        do i=bzil(k),bzil(k+1)-1
          if (i .eq. bzil(k+1)-1) then
            matz=matname(k)
          else
            matz='-1'
          end if
          write(frm, '(''('',I4,''I7,A,A)'')') bzip(i+1)-bzip(i)+1
          write (unitid1,FMT=frm) bzip(i+1)-bzip(i),                    &
     &      (bzi(j),j=bzip(i),bzip(i+1)-1),'   ',matz(1:len_trim(matz))
!
          write (unitid2) bzip(i+1)-bzip(i),                            &
     &      (bzi(j),j=bzip(i),bzip(i+1)-1),matz
        end do
      end do
      close (unitid1)
      close (unitid2)
      return
      end subroutine wnetin
