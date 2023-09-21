      subroutine lese(bzi,bzip,bzil,zpz,zpp,matname,krb,meshsize,       &
     &                lzrb,layanz,layrb,laytxt,lalrb,lbtrb)
      use feminterface, only: getsetting, gengeb, zeit, reallocate, low2hi, string2number
      use globalvariables, only: gkz, xbk, ybk, nnat, kzrb, gzz, zki, alrb, btrb, zrb, gbz
      use femtypes
      implicit none
      integer (I4B) :: layanz
      integer (I4B), pointer :: bzi(:), bzip(:), bzil(:), krb(:,:)
      integer (I4B), pointer :: zpz(:), lzrb(:), layrb(:,:)
      real (DP), pointer :: zpp(:), meshsize(:)
      complex (DPC), pointer :: lalrb(:,:), lbtrb(:,:)
      character (len=*), pointer :: matname(:), laytxt(:)
      intent (out) :: bzi, bzip, bzil, zpz, zpp, matname, krb, meshsize
      intent (out) :: lzrb, layanz, layrb, laytxt, lalrb, lbtrb
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
!    $Revision: 1.18 $
!    $Date: 2014/06/24 13:59:09 $
!    $Author: m_kasper $
!
!  Output:
!            gbz      number of regions
!            gzz      number of branches
!            gkz      number of key-points
!            bzi      list of branches, including the inner branches
!                     in the sequence of regions
!            bzip     Pointer (index) to the start branch of a region in the
!                     list bzi (for compact storage)
!            bzil     Pointer (index) onto bzip including all regions
!                     which are inside a region (multiply connected regions)
!            zki      list of of key-points (geometry nodes) which determine the branches
!                     zki(1,i)  start key-point of the branch
!                     zki(2,i)  end key-point of the branch
!                     zki(3,i)  third key-koint of the branch, with
!                               /  = 0  :  straight line
!                     zki(3,i)  -  > 0  :  node located on an arc
!                               \  < 0  :  midpoint of the arc
!            zrb      Type of the boundary condition of the branch,
!                     which can take the following values:
!                       0 -  99: Dirichlet BC, where for 1- 99: the BC has to be set in the file USERBC
!                     200 - 299: Neumann of general BC
!                     300 - 399: visible contour, no BC
!                     400 - 499: non-visible contour, no BC
!            zpz      number of nodes on the branch (including start-  and endpoint)
!            zpp      control parameter for the distribution of nodes
!                          /  = 0  :  equidistant distribution
!                     zpp  -  > 0  :  nodes are denser at the startpoint
!                          \  < 0  :  nodes are denser at the endpoint
!            xbk,ybk  coordinates of key-points
!            krb      boundary condition of key-points
!            kzrb     branch which is valid for the key-point
!            alrb     parameters of the Dirichlet of general BC
!            btrb     parameters of the Dirichlet of general BC
!
!  local variables
      integer (I4B) :: i, j, unitid, inat, point1, point2
      integer (I4B) :: zwzahl, ios, pos
      real (DP) :: ms
      real (DP), pointer:: area1(:)
      character (len=200) :: path
      character (len=10000) :: str
      logical :: continuation
!      
      external :: grglun
!
!  read the geometry information from the file: PROJECTPATH//netin.acd
!
!  Format of input file netin.acd
!    gkz                                                             number of key-points
!    xbk,ybk                                                         coordinates of key-points
!    ...       for each key-point
!    gzz                                                             number of branches
!    zki(1,.),zki(2,.)zki(3,.),zpz,zpp,lzrb                          start-point, end-point, third point,number of nodes 
!    ...       for each branch                                       .. distribution, index of boundary condition
!    layanz                                                          number of boundary conditions
!    (zrb,alrb,btrb; i=1..nnat), laytxt                              (kind, parameters al, bt) for each nature, name 
!    ...       for each boundary condition
!    gbz                                                             number of regions
!    zwzahl,(bzi(point2-1+j),j=1,zwzahl),matname(i)                  number of branches, branches, material
!    ...       for each region
!
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
      open (unitid,file=path(1:len_trim(path))//'netin.acd',            &
     &   form='formatted',position='REWIND',action='READ',iostat=ios)
      if (ios .ne. 0) then
        print*,'**** Error opening file: ',path(1:len_trim(path))//'netin.acd'
        print*,'**** Error no: ',ios
        return
      end if
!
!  read the key-points
      read (unitid,*) gkz
      allocate(xbk(gkz), ybk(gkz), krb(gkz,nnat), kzrb(gkz,nnat))
      i=0
      do while (i .lt. gkz)
        read(unitid,fmt='(a)',iostat=ios) str
        if (ios .ne. 0) then
          print*,'Error in reading NETIN.ACD:',str
          stop
        end if
!  ignore comment lines
!  at each line, everything following either ! or # or/* or // is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
        if (index(str,'//') .ne. 0) str=str(1:index(str,'//')-1)
!  allow for blank or comment lines
        if (len_trim(str) .eq. 0) cycle
        i=i+1
        read (str,*) xbk(i),ybk(i)
!  default values for boundary condition on key-points
        krb(i,:)=999999
        kzrb(i,:)=0
      end do
!
!  read the branches
      read (unitid,*) gzz
      allocate (zki(3,gzz), zpz(gzz), zpp(gzz), lzrb(gzz))
      i=0
      do while (i .lt. gzz)
        read(unitid,fmt='(a)',iostat=ios) str
        if (ios .ne. 0) then
          print*,'Error in reading NETIN.ACD:',str
          stop
        end if
!  ignore comment lines
!  at each line, everything following either ! or # or/* or // is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
        if (index(str,'//') .ne. 0) str=str(1:index(str,'//')-1)
!  allow for blank or comment lines
        if (len_trim(str) .eq. 0) cycle
        i=i+1
        read (str,*) zki(1,i),zki(2,i),zki(3,i),zpz(i),zpp(i),lzrb(i)
!  correct to minimum number of nodes =3 if necessary
        if (zpz(i) .lt. 3) zpz(i)=3
        if ( (zki(1,i) .gt. gkz) .or. (zki(1,i) .lt. 1)) then
          write (*,3000) i,zki(1,i)
3000      format(' error: key-point 1 of branch ',i5,' incorrect:',i5)
          stop 2
        end if
        if ( (zki(2,i) .gt. gkz) .or. (zki(2,i) .lt. 1) ) then
          write (*,3001) i,zki(2,i)
3001      format(' error: key-point 2 of branch ',i5,' incorrect:',i5)
          stop 2
        end if
      end do
!
!  read the boundary conditions
      read (unitid,*) layanz
      allocate (lalrb(layanz,nnat), lbtrb(layanz,nnat),                 &
     &         layrb(layanz,nnat), laytxt(layanz))
      i=0
      do while (i .lt. layanz)
        read(unitid,fmt='(a)',iostat=ios) str
        if (ios .ne. 0) then
          print*,'Error in reading NETIN.ACD:',str
          stop
        end if
!  ignore comment lines
!  at each line, everything following either ! or # or/* or // is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
        if (index(str,'//') .ne. 0) str=str(1:index(str,'//')-1)
!  allow for blank or comment lines
        if (len_trim(str) .eq. 0) cycle
        i=i+1
        read (str,*,iostat=ios) (layrb(i,j),lalrb(i,j),lbtrb(i,j),j=1,nnat),laytxt(i)
        if (ios .ne. 0) then
!  try once more (no comment), if error
          read (str,*,iostat=ios) (layrb(i,j),lalrb(i,j),lbtrb(i,j),j=1,nnat)
          laytxt(i)=''
        end if
!  format is restrictive:  
!     iii         ( xx.xx , xx.xx )  ...  ( xx.xx , xx.xx )      iii         ( xx.xx , xx.xx )  ...  ( xx.xx , xx.xx )
!   BC (eg 200)   complex value ( , ) needed  ...
        if (ios .ne. 0) then
          print*, '***** Error reading: ',trim(path(1:len_trim(path))//'netin.acd'),' while reading:'
          write (*,'(A,I5,A)'), ' ***** Error reading line:',i,' of section boundary conditions:'
          print*, trim(str)
        end if
      end do
!
!  set zrb, alrb, btrb by use of lzrb, lalrb, lbtrb
      allocate (alrb(gzz,nnat),btrb(gzz,nnat), zrb(gzz,nnat))
      do i=1, gzz
        zrb(i,1:nnat)=layrb(lzrb(i),1:nnat)
        alrb(i,1:nnat)=lalrb(lzrb(i),1:nnat)
        btrb(i,1:nnat)=lbtrb(lzrb(i),1:nnat)
      end do
!  set the boundary condition for key-points
      do i=1, gzz
!  first part
        do inat=1,nnat
          if (zrb(i,inat) .lt. krb(zki(1,i),inat)) then
            krb(zki(1,i),inat)=zrb(i,inat)
            kzrb(zki(1,i),inat)=i
          end if
        end do
!  second part
        do inat=1,nnat
          if (zrb(i,inat) .lt. krb(zki(2,i),inat)) then
            krb(zki(2,i),inat)=zrb(i,inat)
            kzrb(zki(2,i),inat)=i
          end if
        end do
      end do
!
!  read the regions
      read (unitid,*) gbz
!  read the number of regions
      if (gbz .gt. 0) then
!  larger than zero?  then read regions from input file
        allocate(bzi(2*gzz), bzip(2*gbz+1), bzil(gbz+1), matname(gbz), meshsize(gbz))
        point1=1
        point2=1
!  read the regions and and store in the compact format
        continuation=.false.
        i=0
        do while (i .lt. gbz .or. continuation)
          read(unitid,fmt='(a)',iostat=ios) str
          if (ios .ne. 0) then
            print*,'Error in reading NETIN.ACD:',str
            stop
          end if
!  ignore comment lines
!  at each line, everything following either ! or # or/* or // is ignored
          if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
          if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
          if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
          if (index(str,'//') .ne. 0) str=str(1:index(str,'//')-1)
!  allow for blank or comment lines
          if (len_trim(str) .eq. 0) cycle
          if (.not.continuation) then
            i=i+1
            bzil(i)=point1
          end if
!5421      continue
          if (point1 .gt. size(bzip)) then
            bzip=>reallocate(bzip,2*size(bzip))
          end if
          bzip(point1)=point2
!
          read (str,*) zwzahl,(bzi(point2-1+j),j=1,zwzahl),matname(i)
!  find postion of matname
          pos=index(str,matname(i)(1:len_trim(matname(i))),back=.true.)
          pos=pos+len_trim(matname(i))
          str=adjustl(str(pos:))
          call low2hi(str,len_trim(str))
!  find 'meshsize' if given
!  the command meshsize=<value>  is allowed to follow the material name
          pos=index(str,'MESHSIZE')
          if (pos .ne. 0) then
            pos=pos+len_trim('MESHSIZE')
!  position of assignment sign =
            pos=pos+index(str(pos:),'=')
            call string2number(str(pos:),ms)
          else
!  'meshsize' not given
            ms=huge(1._DP)
          end if

          point1=point1+1
          point2=point2+zwzahl
          if (matname(i) .eq. '-1') then
            continuation=.true.
          else
            continuation=.false.
            meshsize(i)=ms
          end if
        end do
        bzil(gbz+1)=point1
        if (point1 .gt. size(bzip)) then
          bzip=>reallocate(bzip,size(bzip)+1)
        end if
        bzip(point1)=point2
!  A region is described by the list of the surrounding branches.
!  The orientation of a branch is reversed by assigning a negative sign
!  to the branch number.
!
!  In the case of a simply connected region, the branches are oriented
!  counter clockwise. In this case there is only an outer boundary.
!  The format of a simply connected region region is:
!
!       number of branches,  list of branches,  Materialname of the region
!
!  A multiply connected region additionally contains a preceding list
!  of the inner bounderies. Each inner boundary is the list of branches
!  in the clockwise orientation. For distinction between inner and outer
!  boundaries, inner regions end with a "-1" instead of the material name.
!  The format of a multiply connected region is:
!
!       number of branches,  list of branches,  -1   (inner boundaries)
!            ...                   ...                ...
!       number of branches,  list of branches,  Materialname of the region
!
!  The geometry below has two regions with in total 8 branches.
!             /\            If we assume that branches are allways oriented
!           /    \          from left to right or from bottom to top we get:
!        2/  _5_   \1
!       /  4|   |6   \        4     4  5 -6 -7            -1
!     /     |_7_|      \      5     3  9  8 -1 -2         triangle
!   /__._______._________\    4     7  6 -5 -4            square
!     3     9       8
!
!  Regions are stored in a compact format using three arrays: bzil, bzip, bzi
!    bzi      is a continous list of branches, including the inner branches.
!             The sequence is the same as in the input file (netin.acd)
!    bzip     is a pointer (index) to the first branch in the list bzi of a
!             region. For multiply connected regions inner boundary branches
!             precede the outer boundary.
!    bzil     for each region is a pointer (index) to the first
!             (possibly inner) region in the array bzip
!
!  Compact storage of the geometry below results in:
!             /\
!           /    \          bzi:   4  5 -6 -7  3  9  8 -1 -2  7  6 -5 -4
!        2/  _5_   \1              |           |              |            |
!       /  4|   |6   \      bzip:  1           5              9           13
!     /     |_7_|      \           |                          |            |
!   /__._______._________\  bzil:  1                          3            4
!     3     9       8
!
      else
!  If the number of regions is zero, we try to automatically retrieve the
!  regions by connecting branches
!  However, all regions will be filled with material NULL
        call gengeb(gbz,gzz,gkz,bzi,bzip,bzil,zki,xbk,ybk,area1)
        deallocate(area1)
        allocate(matname(gbz))
        do i=1,gbz
          matname(i)='NULL'
        end do
      end if
      close(unitid)
      call zeit(' Reading Input and Generating Regions')
      return
      end subroutine lese
