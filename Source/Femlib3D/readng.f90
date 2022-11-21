      subroutine readng(meshfile, ok)
      use feminterface, only: getsetting, strtok, readmat
      use feminterface3d, only: sortnodes, readbc
      use femtypes
      use globalvariables3D, only : sbc, dom, nnat, nod, numdom, dommat, numn, nums, numv
      use globalvariables3D, only : sn, vn, numbc, bctype, pvalue, qvalue, bcnames, domnames
      use matconstants, only : matnames
      implicit none
      character(len=*) :: meshfile
      logical ok
      intent (in) :: meshfile
      intent (out) :: ok
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!    $Revision: 1.15 $
!    $Date: 2015/03/11 16:59:48 $
!    $Author: juryanatzki $
!------------------------------------------------------------------------------
!
!  Read mesh data from NETGEN mesh file. Neutral format is the only mesh input
!  format supported atm. The file containing the mesh has got to be named 
!  basis.mesh. The information contained in these formats is:
!
!  Neutral format:
!   1. Number of nodes (numn)
!   2. List of coordintes (nod)
!   3. Number of volume elements (numv)
!   4. domain (1 entry) and nodes (4 entries) for volume elements (dom, vn)
!   5. Number of surface elements (nums)
!   6. bc at surface (1 entry) and nodes (3 entries) for surface elements (bcs, sn)
!
!  Volume mesh format (NETGEN's internal mesh format) is going to follow. This
!  format contains more information about neighboring domains for surface ele-
!  ments. It seems to contain information about the edges as well to include
!  additional degrees of freedom for isoparametric elements.
!
!  Input:
!      meshfile     name of the mesh.file (netgen-file)
!  Output:
!      ok          =.false. if an error occurs
!------------------------------------------------------------------------------
!  Output:
!    Regions
!      numdom    number of domains / regions
!      dommat    number for material in domain / region
!    Elements
!      numv      number of volume elements
!      vn        volume element -> nodes
!      dom       assignment of element to domain / region
!      nums      number of surface elements
!      sn        surface element -> nodes
!    Nodes
!      numn      number of nodes
!      nod       coordinates of nodes
!    Material
!      matnames  names of materials (material-files)
!    Boundary Conditions
!      sbc       boundary condition(index) of surface elements
!      bctype    type of biunary condition (for each nature)
!      pvalue    p-value of boundary condition
!      qvalue    q-value of boundary condition
!
!  local variables
      integer (I4B) :: i, j, unitid, ios, alstat, minbc, nummaterials, domindex, matindex
      real (DP) :: geoscale
      logical :: found, exist
      logical, allocatable :: isbcpresent(:)
      character (len=200) :: path, str
      character (len=40) :: domnm, matnm
      character (len=40), allocatable :: tmpmatname(:)

!
!  get project path from user's environment variables
!      hid_face=0

      ok = .false.
!  get project path from user's environment variables
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
!  open mesh.neutral file in project path for formatted reading (ASCII)
      inquire (file=path(1:len_trim(path))//meshfile(1:len_trim(meshfile)),     &
     &   exist=exist,iostat=ios)
      if (ios .ne. 0 .or. exist.eq..false.) then
        print*,'***** Mesh-file:',path(1:len_trim(path))//meshfile(1:len_trim(meshfile)),' does not exist'
        return
      end if
      open (unitid,file=path(1:len_trim(path))//meshfile(1:len_trim(meshfile)),  &
     &   form='formatted',position='REWIND',action='READ',iostat=ios)
      if (ios .ne. 0) then
        print*, '**** Error on opening file: ',                         &
     &           path(1:len_trim(path))//meshfile(1:len_trim(meshfile))
        print*, '**** IO Error No: ',ios
        return
      end if

!
!  read number of nodes
      read (unitid,*) numn
      print "(a,i8)","  nodes: ", numn
!  allocate and initialize pointers for coordinates
      allocate(nod(3,numn), stat=alstat)
      nod=0._DP
      if (alstat .ne. 0) then
        ok=.false.
        print*,' *** Error in array allocation'
        return 
      end if
!
!  read the coordinates
      do i=1,numn
        read (unitid,*) nod(:,i)
      end do
      call getsetting('GEOMETRY_FACTOR',geoscale)
      if (geoscale .ne. 1._DP) then
        nod = nod * geoscale
      end if
!
!  read number of volume elements
      read (unitid,*) numv
      print "(a,i8)","  volume elements: ", numv
!  allocate and initialize pointers for domain and volume element -> node map
      allocate(dom(numv),vn(4,numv), stat=alstat)
      dom=0; vn=0
      if (alstat .ne. 0) then
        ok=.false.
        print*,' *** Error in array allocation'
        return 
      end if
!
!  read domain and volume element - node map
      do i=1,numv
        read (unitid,*) dom(i),vn(:,i)
!   print*,'reading vn(:,',i,') = ',vn(:,i)
      end do
!
!  read number of surface elements
      read (unitid,*) nums
      print "(a,i8)","  surface elements: ", nums
!  allocate and initialize pointers for domain and surface element -> node map
      allocate(sn(3,nums), sbc(nums), stat=alstat)
      if (alstat .ne. 0) then
        ok=.false.
        print*,' *** Error in array allocation'
        return 
      end if
!
!  read surface elements
      do i=1,nums
        read(unitid,*) sbc(i), sn(:,i)
      end do
      close (unitid)

!  boundary condition lables must run form 1 to numbc
      minbc = minval(sbc(:))
      if (minbc .le. 0) then
        print*,'Boundary condition',minbc,'found in mesh file'
        print*,'Boundary conditions/ tags must be positive'
        return
      end if
      numbc = maxval(sbc(:))
      allocate (isbcpresent(numbc))
      isbcpresent(:) = .false.
      do i=1, numbc
        if (any(sbc(:) .eq. i)) isbcpresent(i) = .true.
      end do
      if ( .not.all(isbcpresent(:)) ) then
        print*,'Boundary conditions/ tags must in Netgen file must be labled as 1, 2, 3, ...'
      end if
      deallocate(isbcpresent)


      allocate (bctype(numbc,nnat), pvalue(numbc,nnat), qvalue(numbc,nnat),bcnames(numbc))
!  read and set Boundary Condidions on surface vertices according to priority,
      call readbc(bcnames, numbc, nnat, .false., bctype, pvalue, qvalue, ok)

!
      ok=.true.
!
!  read the material names from materials.txt and read material data
      inquire (file=path(1:len_trim(path))//'materials.txt',            &
     &        exist=exist,iostat=ios)
      if (ios .ne. 0 .or. exist.eq..false.) then

!  the file 'materials.txt' does not exist, 
        print*,'File materials.txt is missing'
        stop

      else

!  the file 'materials.txt' exists
!  we interprete the names found in the UNV file to represent domain-names and
!  have to map them to material-names
        call grglun(unitid)
!  open materials.txt file in project path for formatted reading (ASCII)
        open (unitid,file=path(1:len_trim(path))//'materials.txt',      &
     &        form='formatted',position='REWIND',action='READ',iostat=ios)
        if (ios .ne. 0) then
          print*, '**** Error on opening file: ',path(1:len_trim(path))//'materials.txt'
          print*, '**** IO Error No: ',ios
          ok=.false.
          return
        end if

!  the number of domains is the maximum value of dom(:)
        numdom = maxval(dom)
        allocate(dommat(numdom),domnames(numdom))
        allocate(tmpmatname(numdom))
        nummaterials = 0
        domindex = 0
        do
!  read the material names from materials.txt and read material data
!  one line for each domain: first entry domain-name; second entry is material with which the domain is filled
!     domain and material are seperated by = or :
!  e.g.   Membrane = Polymide
!         Conductor = Copper
          read(unitid,fmt='(a)',iostat=ios) str
          if (ios .lt.0) then
            exit
          else if (ios .gt. 0) then
            print*,'**** Error while reading materials.txt'
            ok = .false.
            stop
          end if

!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
          if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
          if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
          if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  allow for blank or comment lines
          if (len_trim(str) .eq.0) cycle
          call strtok(str,     '=,;:>'//char(9)//char(13), domnm)     ! char(9)= horizontal tab; char(13)= CR
          call strtok(char(0), '=,;:>'//char(9)//char(13), matnm)

          if (ios .lt.0) then
            exit
          else if (ios .gt. 0) then
            print*,'**** Error while reading materials.txt'
            ok = .false.
            stop
          end if
!  assign
          domnm = domnm(index(domnm,':')+1:)
          call low2hi(domnm,len_trim(domnm))
          domnm = adjustl(trim(domnm))
          matnm = matnm(index(matnm,':')+1:)
!          call low2hi(matnm,len_trim(matnm))
          matnm = adjustl(trim(matnm))

!  set domain index
          domindex = domindex +1
          domnames(domindex) = domnm

!  find material index
          found = .false.
          do j = 1, nummaterials
            if (matnm .eq. tmpmatname(j)) then
              found = .true.
              matindex = j
              exit
            end if
          end do
          if (found) then
            dommat(domindex) = matindex
          else
            nummaterials = nummaterials + 1
            tmpmatname(nummaterials) = matnm
            call readmat(matnm,exist,dommat(nummaterials))
            if (.not. exist) then
              print*,'the material file ',matnm(1:len_trim(matnm)),' does not exist'
              stop
            end if
          end if
        end do

        close (unitid)
        matnames(1:nummaterials) = tmpmatname(1:nummaterials)
        deallocate(tmpmatname)

      end if
!  end of reading the material file

      ok=.true.
!
!  Sort the nodes in array vn.
      call sortnodes
!
      return
      end subroutine readng
