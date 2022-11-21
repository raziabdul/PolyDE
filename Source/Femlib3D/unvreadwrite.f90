      subroutine writeunv(ok)
      use feminterface, only: getsetting
      use femtypes
      use globalvariables3D, only : dom, nnat, nod, numn, nums, numv,   &
     &                              sn, vn, numdom, dommat, numbc, sbc, &
     &                              bcnames, domnames
      use matconstants, only : matnames
      implicit none
      logical ok
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
!    $Revision: 1.14 $
!    $Date: 2015/11/11 17:25:54 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Write mesh data as unv I-DEAS mesh file. The file containing the mesh 
!  will be named basis.unv. The information contained in this format is:
!
!  Neutral format:
!   1. Number of nodes (numn)
!   2. List of coordintes (nod)
!   3. Number of volume elements (numv)
!   4. domain (1 entry) and nodes (4 entries) for volume elements (dom, vn)
!   5. Number of surface elements (nums)
!   6. bc at surface (1 entry) and nodes (3 entries) for surface elements (sbc, sn)
!
!------------------------------------------------------------------------------
!  Input:
!    Regions
!      numdom      number of domains / regions
!      dommat      number for material in domain / region
!    Elements
!      numv        number of volume elements
!      nums        number of surface elements
!      sn          surface element -> nodes
!      vn          volume element -> nodes
!      dom         assignment of element to domain / region
!    Nodes
!      nod         world coordinates of nodes
!      numn        number of nodes
!
!  Output:
!      ok          =.false. if an error occurs
!
!  local variables
      integer (I4B) :: i, j, unitid, ios, firstindex
      integer (I4B), allocatable:: nsingroup(:), nvingroup(:)
      real (DP) :: geoscale
      logical :: second
      character (len=200) :: path
      character (len=6) :: delimiter='    -1'
!
!  get project path from user's environment variables
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
!  open mesh.neutral file in project path for formatted reading (ASCII)
      open (unitid,file=path(1:len_trim(path))//'basis.unv',            &
     &   form='formatted',position='REWIND',action='WRITE',iostat=ios)
      if (ios .ne. 0) then
        print*, '**** Error on opening file: ',path(1:len_trim(path))//'basis.unv'
        print*, '**** IO Error No: ',ios
        ok=.false.
        return
      end if
      rewind unitid
!
!  write nodes
      write (unitid,fmt='(a6)') delimiter
      write (unitid,fmt='(a4)') '2411'
      print "(a,i8,a)","  writing ", numn," nodes"
!  write coordinates
      call getsetting('GEOMETRY_FACTOR',geoscale)
      if (geoscale .ne. 1._DP) then
        do i=1,numn
!  node label,  export coordinate system number,  displacement coordinate system number,  color
          write (unitid,fmt='(4i10)') i,1,1,1
!  node coordinates in the part coordinate system
          write (unitid,fmt='(1p3e25.16)') nod(1,i)/geoscale, nod(2,i)/geoscale, nod(3,i)/geoscale
        end do
      else
        do i=1,numn
!  node label,  export coordinate system number,  displacement coordinate system number,  color
          write (unitid,fmt='(4i10)') i,1,1,1
!  node coordinates in the part coordinate system
          write (unitid,fmt='(1p3e25.16)') nod(1,i), nod(2,i), nod(3,i)
        end do
      end if

      write (unitid,fmt='(a6)') delimiter
!
!  write elements
      write (unitid,fmt='(a6)') delimiter
      write (unitid,fmt='(a4)') '2412'
      print "(a,i8,a)","  writing ", nums," surface elements"
!  write surface elements
      do i=1,nums
!  fe descriptor id,  physical property table number,  material property table number,  color,  number of nodes on element
!       fe descriptor id =41 for triangles
        write (unitid,fmt='(6i10)') i,41,1,1,1,3
!  node labels defining element
        write (unitid,fmt='(8i10)') sn(1,i), sn(2,i), sn(3,i)
      end do

!  write volume elements
      print "(a,i8,a)","  writing ", numv," volume elements"
      do i=1,numv
!  fe descriptor id,  physical property table number,  material property table number,  color,  number of nodes on element
!       fe descriptor id =111 for tetrahedra
        write (unitid,fmt='(6i10)') i+nums,111,1,1,1,4
!  node labels defining element
        write (unitid,fmt='(8i10)') vn(1,i), vn(2,i), vn(3,i), vn(4,i)
      end do
      write (unitid,fmt='(a6)') delimiter

!  write elements
      write (unitid,fmt='(a6)') delimiter
!  write Permanent Groups (groups of boundary conditions and elements in domains)
      write (unitid,fmt='(a4)') '2467'

!  assign the surface elements to groups

!  number of different groups of surface elements(i.e. number of different boundary conditions)
      allocate (nsingroup(numbc))
      nsingroup(1:numbc)=0

      do i=1, nums
        do j=1, numbc
          if (sbc(i) .eq. j) then
            nsingroup(sbc(i)) = nsingroup(sbc(i)) + 1
            exit
          end if
        end do

      end do

!  now write the groups
!   surface goupes
      do i= 1, numbc
! 1 to 4: group number,  active constraint set no. for group,  active restraint set no. for group,  active load set no. for group
! 5 to 8: active dof set no. for group,  active temperature set no. for group,  active contact set no. for group,  number of entities in group
        write (unitid,fmt='(8i10)') i, 0, 0, 0, 0, 0, 0, nsingroup(i)
        write (unitid,fmt=*) bcnames(i)
        second = .false.
        do j=1, nums
          if (sbc(j) .eq. i) then
            if (second) then
              write (unitid,fmt='(8i10)') 8, firstindex, 0, 0, 8, j, 0, 0
              second = .false.
            else
              firstindex = j
              second = .true.
            end if
          end if
        end do
        if (second) then
          write (unitid,fmt='(4i10)') 8, firstindex, 0, 0
          second = .false.
        end if
      end do

!  count number of elements in each domain
!  number of entities in this group
      allocate (nvingroup(numdom))
      nvingroup(1:numdom)=0
      do i=1, numv
        do j=1, numdom
          if (dom(i) .eq. j) then 
            nvingroup(dom(i)) = nvingroup(dom(i)) + 1
            exit
          end if
        end do
      end do

!   volume element groups
      do i= 1, numdom
        write (unitid,fmt='(8i10)') i+numbc, 0, 0, 0, 0, 0, 0, nvingroup(i)
!  TO DO the matnames are not properly assigned by the domain indices
!        should be the domain-names (domnames) instead of material-names (matnames)
!        this would require to store them as globalvariables
        write (unitid,fmt=*) domnames(i), matnames(i), dommat(i)
        second = .false.
        do j=1, numv
          if (dom(j) .eq. i) then
            if (second) then
              write (unitid,fmt='(8i10)') 8, firstindex+nums, 0, 0, 8, j+nums, 0, 0
              second = .false.
            else
              firstindex = j
              second = .true.
            end if
          end if
        end do
        if (second) then
          write (unitid,fmt='(4i10)') 8, firstindex+nums, 0, 0
          second = .false.
        end if
      end do

      write (unitid,fmt='(a6)') delimiter
!
      close (unitid)
      deallocate(nsingroup, nvingroup)
      ok=.true.
      return
      end subroutine writeunv



      subroutine readunv(meshfile,ok)
      use feminterface, only: getsetting, reallocate, low2hi, strtok, readmat, realloc
      use feminterface3d, only: readbc, sortnodes
      use femtypes
      use globalvariables3D, only : numn, nod , dom, nums, numv, sn, sbc,    &
     &                              vn, numdom, dommat, nnat, numbc, bctype, &
     &                              bcnames, pvalue, qvalue, domnames
      use matconstants, only : matnames, maxmat
      implicit none
      character(len=*) :: meshfile
      logical ok
      intent (in) :: meshfile
      intent (out) :: ok
!
!------------------------------------------------------------------------------
!
!  Read mesh data in the  unv-format (I-DEAS mesh file), 
!  the file of boundary conditions: 'boundary conditions.txt',
!  the material specification: 'materials.txt' and 
!  the material files '... .mat'.
!
!  The unv-file is expected to only have Dataset numbers
!      2411 (nodes),   2412 (elements) and   2467 (groups); 
!  everything else will be ignorded.
!
!------------------------------------------------------------------------------
!  Input:
!      meshfile     name of the unv.file
!  Output:
!      ok          =.false. if an error occurs
!
!  local variables
      integer (I4B) :: i, j, unitid, ios, datasetnumber, linenumber
      integer (I4B) :: iecs, idcs, icol, sizeel, sizeen, sizegr, a1, a2, a3, a4, a5, a6
      integer (I4B) :: elnum, feid, phyp, matpro, eltype, type1, type2, k, l, nsingroup, lastelement, ns, nv
      integer (I4B) :: numgroups, nummaterials, domindex, matindex, oldsize, n_bc, n_dom, length
      integer (I4B), allocatable :: tag_type(:), bcface(:), matdomain(:)
      integer (I4B), allocatable:: entity_tag(:), elementtype(:), vertex(:,:)
      real (DP) :: geoscale
      logical :: exist, found
      character (len=200) :: path, str
      character (len=6) :: char6
      character (len=40) :: domnm, matnm
      character (len=40), allocatable :: tmpmatname(:)
      character (len=40), allocatable :: groupnames(:)

!  Variables, defined in globalvariables3d
!      numn                 number of nodes
!      numv                 number of volume elements
!      nums                 number of surfaces (elements)
!      numdom...............number of domains (with UNV equal to number of materials)
!      numbc                Number of boundary conditions
!      sbc(:)               BC-index of a surface element               allocate:  sbc(nums)
!      nod(3,:)             coordinates of nodes                                   nod(3,numn)
!      vn(4,:)              vertices(nodes) of elements                            vn(4,numv)
!      sn(3,:)              nodes of surface element                               sn(3,nums)
!      dom(:)               assignment of element to domain / region               dom(numv)
!      dommat(:)            Material index of a domain                             dommat(numdom)
!      bcnames(:)           Text-specifier of boundary conditions                  bcnames(numbc)
!      bctype(:,:)          BC-type e.g. 0, 100, 200 ,...                          bctype(numbc,nnat)
!      pvalue(:,:)          p-value of boundary condition                          pvalue(numbc,nnat)
!      qvalue(:,:)          q-value of boundary condition                          qvalue(numbc,nnat)
!  from module: matconstants
!      matnames(:)          Material names (Text) of each domain                   matnames(numdom)


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
      linenumber = 0
!  initialize, allocate arrays
      sizeen  = 100000
      sizeel  = 500000
      allocate(nod(3,sizeen))
      allocate(vertex(4,sizeel))
      allocate(elementtype(sizeel))
      elementtype(:) = -1
      numgroups = 0

!  first line must be a delimiter
      read (unitid,fmt='(a6)',iostat=ios) char6
      linenumber = linenumber + 1
      if (ios .ne. 0 .or. char6(5:6) .ne. '-1') then
        print*, '**** Error reading delimiter in unv-file'
        print*, '**** IO Error No: ',ios,' in line:',linenumber
        ok=.false.
        return
      end if

!  start to read unv-file
      do

!  get dataset number
        read (unitid,fmt=*,iostat=ios) datasetnumber
        linenumber = linenumber + 1
        if (ios .ne. 0 .or. datasetnumber .eq. -1) then
          print*, '**** Error reading delimiter in unv-file'
          print*, '**** IO Error No: ',ios,' in line:',linenumber
          ok=.false.
          return
        end if

        select case(datasetnumber)
!  The file is expected to only have Dataset numbers  2411 (nodes), 2412 (elements), 2467 (groups)
!  everything else will be ignorded
        case (2411)
!  read nodes
          print "(a)","  reading nodes"
!  read coordinates
          numn=0
!  node label,  export coordinate system number,  displacement coordinate system number,  color
          do
            read (unitid,fmt='(4i10)',iostat=ios) i,iecs,idcs,icol
            linenumber = linenumber + 1
            if (i.eq.-1 .or. ios .ne. 0) then
!  end of dataset
              exit
            end if
!  reallocate if necessary
            numn=max(numn,i)
            if (numn .gt.sizeen) then
              sizeen = 2*sizeen
              nod => reallocate(nod,3,sizeen)
            end if
!  node coordinates in the part coordinate system
            read (unitid,fmt='(1p3e25.16)') nod(1,i), nod(2,i), nod(3,i)
            linenumber = linenumber + 1
          end do
          call getsetting('GEOMETRY_FACTOR',geoscale)
          if (geoscale .ne. 1._DP) then
            nod = nod * geoscale
          end if

        case (2412)
!  read elements
          print "(a,i6,a)","  reading elements"
          nums  = 0
          numv  = 0
          lastelement = 0

          do
!  fe descriptor id,  physical property table number,  material property table number,  color,  number of nodes on element
            read (unitid,fmt='(6i10)',iostat=ios) elnum, feid, phyp, matpro, icol, eltype
            linenumber = linenumber + 1
            if (elnum.le.0 .or. ios .ne. 0) then
!  end of dataset
              exit
            end if

            if (elnum .gt. sizeel) then
              oldsize = sizeel
              sizeel = max(elnum,2*sizeel)
              call realloc(vertex,4,sizeel)
              call realloc(elementtype,sizeel)
              elementtype(oldsize+1:sizeel) = -1
            end if
            lastelement = max(lastelement,elnum)
!  store
            select case (eltype)
            case(2)
!  a line element (beam neglect)
              read (unitid,fmt=*)
              read (unitid,fmt=*)
              linenumber = linenumber + 2
              elementtype(elnum) = 2

            case(3)
!  a triangle (surface element)
!  fe descriptor id =41 for triangles
              if (feid .ne. 41) then
                print*,'error reading unv-file',' in line:',linenumber
                print*,'expected to receive triangle element of type 41, but received ',feid
              end if
!  node labels defining element
              read (unitid,fmt='(8i10)') vertex(1,elnum), vertex(2,elnum), vertex(3,elnum)
              linenumber = linenumber + 1
              elementtype(elnum) = 3
              nums = nums +1
            case(4)
!  a tetrahedron
!  fe descriptor id =111 for tetrahedra
              if (feid .ne. 111) then
                print*,'error reading unv-file',' in line:',linenumber
                print*,'expected to receive tretrahedron element of type 111, but received ',feid
              end if
!  node labels defining element
              read (unitid,fmt='(8i10)') vertex(1,elnum), vertex(2,elnum), vertex(3,elnum), vertex(4,elnum)
              linenumber = linenumber + 1
              elementtype(elnum) = 4
              numv = numv +1

            case default
              print*,'error reading unv-file',' in line:',linenumber
              print*,'error: unknown element type',eltype
              pause
              stop
            end select
          end do

        case (2467)
!  read Permanent Groups (Groups of boundary conditions and elements in domains)
          sizegr=500
          allocate (entity_tag(sizeel), groupnames(sizegr))
          entity_tag = -1
          do
!! 1 to 4: group number,  active constraint set no. for group,  active restraint set no. for group,  active load set no. for group
!! 5 to 8: active dof set no. for group,  active temperature set no. for group,  active contact set no. for group,  number of entities in group
            read (unitid,fmt='(8i10)') i, a1, a2, a3, a4, a5, a6, nsingroup
            linenumber = linenumber + 1
            if (i.eq.-1 .or. ios .ne. 0) then
              exit
            end if
            numgroups=max(numgroups,i)
            if(numgroups .gt. sizegr) then
              call realloc(groupnames,2*sizegr,40)
              sizegr = 2*sizegr
            end if

            read (unitid,fmt='(A40)') groupnames(i)
            groupnames(i) = adjustl(groupnames(i))
            call low2hi(groupnames(i),len_trim(groupnames(i)))
            linenumber = linenumber + 1
!  now read the group
            do j=1, nsingroup/2
              read (unitid,fmt='(8i10)') type1, k, a1, a2, type2, l, a3, a4
              linenumber = linenumber + 1
              if (max(k,l) .gt. sizeel) then
                call realloc(entity_tag,max(2*sizeel,k,l))
                entity_tag(sizeel+1:max(2*sizeel,k,l)) = -1
                sizeel = max(2*sizeel,k,l)
              end if
              if (type1 .eq. 8) then
                if (entity_tag(k) .ne. -1) then
                  print*,'error reading unv-file',' in line:',linenumber
                  print*,'an element has multiple tags ',len_trim(groupnames(i)), ' and ',len_trim(groupnames(entity_tag(k)))
                end if
                if (entity_tag(l) .ne. -1) then
                  print*,'error reading unv-file',' in line:',linenumber
                  print*,'an element has multiple tags ',len_trim(groupnames(i)), ' and ',len_trim(groupnames(entity_tag(l)))
                end if
                entity_tag(k) = i
                entity_tag(l) = i
              else
!  unnknown type
                print*,'error reading unv-file',' in line:',linenumber
                print*,'unexpected type',type1
              end if

            end do
            if (mod(nsingroup,2) .eq. 1) then
              read (unitid,fmt='(8i10)') type1, k, a1, a2
              linenumber = linenumber + 1
              if (k .gt. sizeel) then
                call realloc(entity_tag,sizeel+1)
                entity_tag(sizeel+1:sizeel+1) = -1
                sizeel = sizeel+1
              end if
              if (entity_tag(k) .ne. -1) then
                print*,'error reading unv-file',' in line:',linenumber
                print*,'an element has multiple tags ',len_trim(groupnames(i)), ' and ',len_trim(groupnames(entity_tag(k)))
              end if
              entity_tag(k) = i
            end if

          end do

        case default
!  some other datasetnumber => ignore
          do
!  advance to next delimiter (i.e. ignore this section)
            read (unitid,fmt='(a6)',iostat=ios) char6
            linenumber = linenumber + 1
            if (ios .lt. 0) then
!  End of File
              print*,'error reading unv-file',' in line:',linenumber
              print*,'unexpected end of file reading ',path(1:len_trim(path))//meshfile(1:len_trim(meshfile))
              print*,'number of nodes:',numn
              print*,'number of elements:',numv
              print*,'number of surfaces',nums
              close (unitid)
              return
            else if (ios .gt. 0) then
!  Error
              print*,'error reading unv-file',' in line:',linenumber
              print*,'error reading ',path(1:len_trim(path))//meshfile(1:len_trim(meshfile)),' : ',char6
              pause
              stop
            end if
!  Check for delimiter
            if (char6(5:6) .eq. '-1') then
!  found delimiter
              exit
            end if
          end do
        end select

!  find next delimiter
        read (unitid,fmt='(a6)',iostat=ios) char6
        linenumber = linenumber + 1
        if (ios .lt. 0) then
!  End of File, terminates as expected -> return
          print*,'number of nodes:   ',numn
          print*,'number of elements:',numv
          print*,'number of surfaces:',nums
          close (unitid)
          exit
          return
        else if (char6(5:6) .ne. '-1') then
          print*, '**** Error reading delimiter in unv-file in line:',linenumber
          print*, '**** IO Error No: ',ios
          close (unitid)
          ok=.false.
          return
        end if

      end do
!  end of reading unv-file

      nod => reallocate(nod,3,numn)
!
!  print extends
      print*, 'Extends of Geometry'
      write(*,'(a,2g13.4)') 'x: ',minval(nod(1,:)),maxval(nod(1,:))
      write(*,'(a,2g13.4)') 'y: ',minval(nod(2,:)),maxval(nod(2,:))
      write(*,'(a,2g13.4)') 'z: ',minval(nod(3,:)),maxval(nod(3,:))

!   TO DO=========  check for inner surfaces having a tag ??  ====================================
      allocate(tag_type(numgroups))

!  check whether the tag is a BC or a domain/ material 
      do i=1, lastelement
        select case (elementtype(i))
        case(3)
!  mark tag to be a bc
          if (entity_tag(i) .le. 0) then
!  The surface-element  i  has been found in the list of surface elements (triangle)
!  However, this element does not appear in a Group-list
     !       print*,'remark, interprete surface element without tag to be an interface'
     !       print 100,'element: ',i,'x =',sum(nod(1,vertex(1:3,i)))/3.,  &
     !&                              'y =',sum(nod(2,vertex(1:3,i)))/3.,  &
     !&                              'z =',sum(nod(3,vertex(1:3,i)))/3.
100         format (A,I7,3(' ;',3X,A,1X,G10.3))
            ok = .false.
            entity_tag(i)=300       ! a inner surface, i.e. an interface
          else if(adjustl(groupnames(entity_tag(i))) .eq. 'INTERFACE') then
            tag_type(entity_tag(i))=300
            entity_tag(i)=300
          else
            tag_type(entity_tag(i))=3
          end if
        case(4)
!  mark tag to be a domain/ material
          if (entity_tag(i) .le. 0) then
!  The volume-element  i  has been found in the list of volume elements (tetrahedra)
!  However, this element does not appear in a Group-list
            print*,'error, a volume element has no tag (Material)'
            print 100,'element: ',i,'x =',sum(nod(1,vertex(1:4,i)))/4.,  &
     &                              'y =',sum(nod(2,vertex(1:4,i)))/4.,  &
     &                              'z =',sum(nod(3,vertex(1:4,i)))/4.
            ok = .false.
            return
          else
            tag_type(entity_tag(i))=4
          end if
        end select
      end do

!  find numdom, numbc
      numbc = 0
      numdom = 0
      do i=1,numgroups
        if (tag_type(i) .eq. 4) then
          numdom = numdom + 1
        else if(tag_type(i) .eq. 300) then
! numbc does not change
        else
          numbc = numbc + 1
        end if
      end do
      if (numdom .gt. maxmat) then
        print*,'error, the number of domains/ materials is too large'
      end if


!  generate the list of BC's and Domains
!   TO DO=========  allocate appropriately (with given size / not a fixed size)      allocate (matnames(numdom/nummat))
!
!  check wether a groupname 'INTERFACE' exists as a BCname otherwise append this
      found = .false.
      do i=1,numgroups
        if (groupnames(i) .eq. 'INTERFACE') then
          found = .true.
          k = i
          exit
        end if
      end do
      if (.not. found) then
        numbc = numbc + 1
      else
!
        if (k .ne. numbc) then
          print*, 'it seems that the name ''INTERFACE'' was used for a Boundary Condition'
          print*, 'this is not allowed ''INTERFACE'' a reserved name'
        end if
      end if
      allocate(bcface(numgroups), matdomain(numgroups), bcnames(numbc), domnames(numdom))
      bcnames(numbc) = 'INTERFACE'

!  assign bcnames, domnames indices to bcface, matdomain
      n_bc = 0
      n_dom = 0
      do i=1,numgroups
        if (tag_type(i) .eq. 3) then
          n_bc = n_bc + 1
          bcnames(n_bc) = groupnames(i)
          bcface(i) = n_bc
        else if (tag_type(i) .eq. 4) then
          n_dom = n_dom + 1
          domnames(n_dom) = groupnames(i)
          matdomain(i) = n_dom
        end if
      end do
      deallocate(tag_type)

      print*,'Boundary conditions found in mesh file:'
      do i=1,numbc-1
        print*,i,' : ',bcnames(i)
        bcnames(i)=bcnames(i)(index(bcnames(i),':')+1:)
        call low2hi(bcnames(i),len_trim(bcnames(i)))
      end do

      print*,'Materials/ domaines  found in mesh file:'
      do i=1,numdom
        print*,i,' : ',domnames(i)
        domnames(i)=domnames(i)(index(domnames(i),':')+1:)
        call low2hi(domnames(i),len_trim(domnames(i)))
      end do

!  sbc   is the boundary condition of a face
      allocate (dom(numv), sbc(nums), dommat(numdom))
      allocate (vn(4,numv), sn(3,nums))
!  relabeling elements
      ns = 0
      nv = 0
      do i=1, lastelement
        select case (elementtype(i))
        case(3)
!  relabel surface-elements
          ns = ns +1
          sn(1,ns) = vertex(1,i)
          sn(2,ns) = vertex(2,i)
          sn(3,ns) = vertex(3,i)
          if (entity_tag(i) .eq. 300) then
!  in the case of interface
            sbc(ns) = numbc
          else
            sbc(ns) = bcface(entity_tag(i))
          end if
        case(4)
!  relabel volume-elements
          nv = nv +1
          vn(1,nv) = vertex(1,i)
          vn(2,nv) = vertex(2,i)
          vn(3,nv) = vertex(3,i)
          vn(4,nv) = vertex(4,i)
          dom(nv) = matdomain(entity_tag(i))
        end select
      end do


      deallocate (groupnames, entity_tag, elementtype, vertex)
      deallocate (matdomain, bcface)

      ok =.true.
      allocate (bctype(numbc,nnat), pvalue(numbc,nnat), qvalue(numbc,nnat))
      bctype(numbc,:) = 300
      pvalue(numbc,:) = cmplx(0._DP,0._DP)
      qvalue(numbc,:) = cmplx(0._DP,0._DP)
!  read and set Boundary Condidions on surface vertices according to priority,
      call readbc(bcnames, numbc-1, nnat, .true., bctype, pvalue, qvalue, ok)

!  read the file: materials.txt
      inquire (file=path(1:len_trim(path))//'materials.txt',            &
     &        exist=exist,iostat=ios)
      if (ios .ne. 0 .or. exist.eq..false.) then

!  the file 'materials.txt' does not exist, 
!  we assume the domain-names found in the UNV file to represent material-names
        matnames(1:numdom)=domnames(1:numdom)
        dommat = (/(i,i=1,numdom)/)

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
!
!  read the material names from materials.txt and read material data
!  one line for each domain: first entry domain-name; second entry is material with which the domain is filled
!     domain and material are seperated by = or :
!  e.g.   Membrane = Polymide
!         Conductor = Copper
        allocate(tmpmatname(numdom))
        nummaterials = 0
        do
!  one line for each domain: first entry domain-name second entry is material with which the domain is filled
          read(unitid,fmt='(a)',iostat=ios) str
          if (ios .lt.0) then
            exit
          else if (ios .gt. 0) then
            print*,'**** Error while reading materials.txt'
            ok = .false.
            pause
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
            pause
            stop
          end if
!  assign
          domnm = domnm(index(domnm,':')+1:)
          call low2hi(domnm,len_trim(domnm))
          domnm = adjustl(trim(domnm))
          matnm = matnm(index(matnm,':')+1:)
!          call low2hi(matnm,len_trim(matnm))
          matnm = adjustl(trim(matnm))

!  find domain index
          found = .false.
          do j = 1, numdom
            length=min(len_trim(domnm),len_trim(domnames(j)))
            if (domnm(1:length) .eq. domnames(j)(1:length)) then
              found = .true.
              domindex = j
              exit
            end if
          end do
          if ( .not. found) then
            print*,'the file materials.txt lists the domain: ',         &
     &            domnm(1:len_trim(domnm)),' which was not found in the UNV file'
            cycle
          end if

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
            call readmat(matnm,exist,dommat(domindex))
            if (.not. exist) then
              print*,'the material file ',matnm(1:len_trim(matnm)),' does not exist'
              pause
              stop
            end if
          end if
        end do

        close (unitid)
        matnames(1:nummaterials) = tmpmatname(1:nummaterials)
        deallocate(tmpmatname)

      end if
!  end of reading the material file
!
!  Sort the nodes in array vn
      call sortnodes

      return
      end subroutine readunv



      subroutine readbc(bcnames, numbc, nnat, hasnames, bctype, pvalue, qvalue, ok)
      use feminterface, only: strtok, getsetting, low2hi
      use json_module
      use femtypes
      implicit none
      integer(I4B) :: numbc, nnat, bctype(:,:)
      complex(DPC) :: pvalue(:,:), qvalue(:,:)
      character (len=*)              :: bcnames(:)
      logical :: ok, hasnames
      intent (in) :: numbc, nnat, hasnames
      intent (out) :: bctype, pvalue, qvalue, ok
      intent (inout) :: bcnames
!
!  Read the file of boundary conditions: 'boundary conditions.txt'
!
!  Each line defines a boundary condition which consits of the bc-descriptors for each nature and the bcname:
!      bc-descriptor1, bc-descriptor2, ... , bcname
!
!  For each nature, boundary condition (i.e. the bc-descriptor) consits of:
!          bctype, 'pvalue'= pvalue, 'qvalue'= pvalue,
!  which are sparated by a delimiter. pavlue and qvalue are complex in general
!
!  The imaginary part may be omitted, but qvalue must be present even if not needed (e.q. in the case of a Dirichlet BC)
!    Expamples of valid descriptors are
!               0, pvalue = 1, 0. qvalue= 1.,1
!        or   200, pvalue = 0     qvalue= 1.
!
!  the bcname is a text which for beeing recognized has to be the same as in the mesh file
!
!      hasnames           = .false. if bcnames are not provided from the mesh file (like with netgen)
!
!  local variables
      integer (I4B)                  :: nature, j, k, l, ios, linenumber, unitid, bcvalue(nnat), bcindex
      real (DP)                      :: pvalue_r(nnat), pvalue_i(nnat), qvalue_r(nnat), qvalue_i(nnat)
      character (len=1000)           :: str
      character (len=80)             :: token, bc_token, bcnam, bc_filename
      character (len=20)             :: delim
      character (len=200)            :: path
      character (len=:), allocatable :: settingsfile, chrctr, error_msg
      logical                        :: imaginary_part, set(numbc), found, status_ok
      type(json_file)                :: json
!
      call getsetting('PROJECTPATH',path)
      call getsetting('BC_FILENAME',bc_filename)
!_
! 
!   TO DO:  clean-up here
!      if (path(len_trim(path)-4:len_trim(path)) == 'json') then
!        call json%load_file(filename = path(1:len_trim(path)))
!        if (json_failed()) then
!          call json_check_for_errors(status_ok, error_msg)
!          print*,'***ERROR reading json:', error_msg
!          pause
!          stop
!        else
!          call json%get('projectPathInfo.projectFile',settingsfile,found)    
!        end if
!        call json%destroy()
!
!        return
!      end if
!__
!  
      ok=.true.
      call grglun(unitid)
!  open the file: boundary conditions.txt in project path for formatted reading (ASCII)
      open (unitid,file=path(1:len_trim(path))//bc_filename(1:len_trim(bc_filename))//'.txt',            &
     &      form='formatted',position='REWIND',action='READ',iostat=ios)
      if (ios .ne. 0) then
        print*, '**** Error opening file: ',path(1:len_trim(path))//bc_filename(1:len_trim(bc_filename))//'.txt'
        print*, '**** IO Error No: ',ios
        ok=.false.
        return
      end if
!
!  Delimiter for tokenizing a string
      delim=' =,;:"{}()!@#$%^&*'
      linenumber = 0
      bcindex = 0
      set = .false.
      do
        read(unitid,fmt='(a)',iostat=ios) str
        if (ios .lt. 0) then
!  reached end of file
          exit
        else if (ios .gt. 0) then 
          print*, '**** Error opening file: ',path(1:len_trim(path))//bc_filename(1:len_trim(bc_filename))//'.txt'
          print*,' in line',linenumber,' : ','"'//str(1:77)//'"'
          ok=.false.
          exit
        end if
        linenumber=linenumber+1
!
!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  allow for blank or comment lines
        if (len_trim(str) .eq.0) cycle
!
        bcindex = bcindex + 1
        if (.not. hasnames .and. bcindex.gt.numbc) exit
!  allow leading blanks
        str=adjustl(str)
!  we do not distinguish between lower and upper case
        call low2hi(str,len_trim(str))

!  read boundary conditions for all natures
        do nature=1,nnat
          if (nature .eq. 1) then
            call strtok(str, delim, token)
          else
            call strtok(char(0), delim, token)
          end if
          if (token(1:1) .ne. char(0)) then
!  first token is a integer number, i.e. the BC-type
            read(token,*,IOSTAT=ios) bcvalue(nature)
            if (ios .ne. 0) then
              write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
              print*,trim(str)
              ok=.false.
              exit
            end if
          else
            write(*,'(a,i3,a)') '*ERRROR: reading line:',linenumber,' of file ''boundary conditions.txt'''
            print*,trim(str)
          end if
!  get next token
          call strtok(char(0), delim, token)
!  second token is a float number, i.e. pvalue (real part)
          if (token(1:1) .ne. 'P' ) then
            write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
            write(*,'(a,i3,a)') ' expecting  ''P=...'' for nature',nature,' but received:'
            print*,trim(str)
            ok=.false.
            exit
          end if
!  get the next token   -> a number, real part of P
          call strtok(char(0), delim, token)
          if (token(1:1) .ne. char(0)) then
!  convert the string:"token" to a number
            read(token,*,IOSTAT=ios) pvalue_r(nature)
            if (ios .ne. 0) then
              write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
              print*,trim(str)
              write(*,'(2a)') '*ERRROR: expecting a number, conversion error for: ',token
              ok=.false.
              exit
            end if
          else
            write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
            print*,trim(str)
            ok=.false.
            exit
          end if

!  get the next token   -> either a float number (imaginary part of pvalue)
!                          or a text 'Q=...'  if the imaginary part is missing
          call strtok(char(0), delim, token)
          if (token(1:1) .ne. char(0)) then
!  convert the string:"token" to a number
            read(token,*,IOSTAT=ios) pvalue_i(nature)
            if (ios .eq. 0) then
!  get the next token   -> a text starting with 'Q'
              call strtok(char(0), delim, token)
              imaginary_part=.true.
            else
!  it is not a number, imaginary part of P is missing
              pvalue_i(nature)=0._DP
              imaginary_part=.false.
            end if
          else
!  we are at the end, Q=part and the BC-name are missing
            write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
            print*,trim(str)
            write(*,'(2a)') '*ERRROR: expecting a number, conversion error for: ',token
            ok=.false.
            exit
          end if
!
          if (token(1:1) .ne. 'Q' ) then
            write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
            write(*,'(a,i3,a)') ' expecting:  ''Q=...'' for nature',j,' but received:'
            print*,trim(str)
            ok=.false.
            exit
          end if
!  get the next token   -> a number, real part of Q
          call strtok(char(0), delim, token)
          if (token(1:1) .ne. char(0)) then
!  convert the string:"token" to a number
            read(token,*,IOSTAT=ios) qvalue_r(nature)
            if (ios .ne. 0) then
              write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
              print*,trim(str)
              write(*,'(2a)') '*ERRROR: expecting a number, conversion error for: ',token
              ok=.false.
              exit
            end if
          else
            write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
            print*,trim(str)
            ok=.false.
            exit
          end if
          
!  get the next token   -> either a float number (imaginary part of pvalue)
!                          or if no imaginary part is present the next nature 
!                          or a text if this was the last nature
          if (imaginary_part) then
            call strtok(char(0), delim, token)
            if (token(1:1) .ne. char(0)) then
!  convert the string:"token" to a number
              read(token,*,IOSTAT=ios) qvalue_i(nature)
              if (ios .ne. 0) then
                write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
                print*,trim(str)
                write(*,'(2a)') '*ERRROR: expecting a number, conversion error for: ',token
                ok=.false.
                exit
              end if
            else
!  we are at the end, at least the BC-name is missing
              write(*,'(a,i3,a)') '*ERRROR: in line:',linenumber,' of file ''boundary conditions.txt'''
              print*,trim(str)
              ok=.false.
              exit
            end if
          else
!  it is not a number, imaginary part of P is missing
            qvalue_i(nature)=0._DP
          end if

        end do
!  end reading the natures 
!
!now read last token, the bc-name
!  
        call strtok(char(0),',;', token)

! if delimiter was comma or semicolon, it may happen that the token was empty
        if (len_trim(token) .eq. 0) call strtok(char(0),',;', token)
        bc_token=adjustl(token(index(token,':')+1:))
        if (.not. hasnames) bcnames(bcindex) = bc_token
!   check bc_name and assign bcvalue palue_r, pvalue_i, qvalue_r, qvalue_i
        l=0
        do k=1, numbc
          bcnam=bcnames(k)(index(bcnames(k),':')+1:)
          call low2hi(bcnam,len_trim(bcnam))
          bcnam=adjustl(bcnam)
          if (bc_token .eq. bcnam) then
            if(set(k)) then
              print*,'repeated assignment of bounadry condition: ',bcnam(1:len_trim(bcnam)),' in line:',linenumber
            end if
            bctype(k,:)=bcvalue(:)
            pvalue(k,:)=cmplx(pvalue_r(:),pvalue_i(:),DPC)
            qvalue(k,:)=cmplx(qvalue_r(:),qvalue_i(:),DPC)
            set(k)=.true.
            l=k
          end if
        end do
        if (l .eq. 0) then
          print*,'the boundary condition ',bc_token(1:len_trim(bc_token)),' cannot be assigned'
        end if
! end of processing this line, proceed to read next line
      end do
      if (any(set .ne. .true.)) then
        print*,'the file  ''boundary conditions.txt'' does not have a definition for:'
        do j=1,numbc
          if (set(j) .ne. .true.) then
            print*,bcnames(j)(1:len_trim(bcnames(j)))
            bctype(j,:)=0
            pvalue(j,:)=cmplx(0._DP,0._DP,DPC)
            qvalue(j,:)=cmplx(0._DP,0._DP,DPC)
            ok=.false.
          end if
        end do
      end if

      close (unitid)
      if (ok .eq. .false.) then 
       pause
       stop
      end if
      return
      end subroutine readbc



! ======================= Dataset Number 2411 =======================================
!
!  DOCUMENTATION OF I-DEAS (UNV-FILE) DATASET NUMBERS
!
! ======================= Dataset Number 2411 =======================================
!Record 1:        FORMAT(4I10)
!                 Field 1       -- node label
!                 Field 2       -- export coordinate system number
!                 Field 3       -- displacement coordinate system number
!                 Field 4       -- color
!Record 2:        FORMAT(1P3D25.16)
!                 Fields 1-3    -- node coordinates in the part coordinate system
!
!Records 1 and 2 are repeated for each node in the model.
!
! ======================= Dataset Number 2412 =======================================
!Record 1:        FORMAT(6I10)
!                 Field 1       -- element label
!                 Field 2       -- fe descriptor id
!                 Field 3       -- physical property table number
!                 Field 4       -- material property table number
!                 Field 5       -- color
!                 Field 6       -- number of nodes on element
!
!Record 2:  *** FOR NON-BEAM ELEMENTS ***
!                 FORMAT(8I10)
!                 Fields 1-n    -- node labels defining element
! 
!Record 2:  *** FOR BEAM ELEMENTS ONLY ***
!                 FORMAT(3I10)
!                 Field 1       -- beam orientation node number
!                 Field 2       -- beam fore-end cross section number
!                 Field 3       -- beam  aft-end cross section number
!
!Record 3:  *** FOR BEAM ELEMENTS ONLY ***
!                 FORMAT(8I10)
!                 Fields 1-n    -- node labels defining element
!
!Records 1 and 2 are repeated for each non-beam element in the model.
!Records 1 - 3 are repeated for each beam element in the model.
!FE Descriptor Id definitions
!   11  Rod
!   21  Linear beam
!   22  Tapered beam
!   23  Curved beam
!   24  Parabolic beam
!   31  Straight pipe
!   32  Curved pipe
!   41  Plane Stress Linear Triangle
!   42  Plane Stress Parabolic Triangle
!   43  Plane Stress Cubic Triangle
!   44  Plane Stress Linear Quadrilateral
!   45  Plane Stress Parabolic Quadrilateral
!   46  Plane Strain Cubic Quadrilateral
!   51  Plane Strain Linear Triangle
!   52  Plane Strain Parabolic Triangle
!   53  Plane Strain Cubic Triangle
!   54  Plane Strain Linear Quadrilateral
!   55  Plane Strain Parabolic Quadrilateral
!   56  Plane Strain Cubic Quadrilateral
!   61  Plate Linear Triangle
!   62  Plate Parabolic Triangle
!   63  Plate Cubic Triangle
!   64  Plate Linear Quadrilateral
!   65  Plate Parabolic Quadrilateral
!   66  Plate Cubic Quadrilateral
!   71  Membrane Linear Quadrilateral
!   72  Membrane Parabolic Triangle
!   73  Membrane Cubic Triangle
!   74  Membrane Linear Triangle
!   75  Membrane Parabolic Quadrilateral
!   76  Membrane Cubic Quadrilateral
!   81  Axisymetric Solid Linear Triangle
!   82  Axisymetric Solid Parabolic Triangle
!   84  Axisymetric Solid Linear Quadrilateral
!   85  Axisymetric Solid Parabolic Quadrilateral
!   91  Thin Shell Linear Triangle
!   92  Thin Shell Parabolic Triangle
!   93  Thin Shell Cubic Triangle
!   94  Thin Shell Linear Quadrilateral
!   95  Thin Shell Parabolic Quadrilateral
!   96  Thin Shell Cubic Quadrilateral
!   101 Thick Shell Linear Wedge
!   102 Thick Shell Parabolic Wedge
!   103 Thick Shell Cubic Wedge
!   104 Thick Shell Linear Brick
!   105 Thick Shell Parabolic Brick
!   106 Thick Shell Cubic Brick
!   111 Solid Linear Tetrahedron
!   112 Solid Linear Wedge
!   113 Solid Parabolic Wedge
!   114 Solid Cubic Wedge
!   115 Solid Linear Brick
!   116 Solid Parabolic Brick
!   117 Solid Cubic Brick
!   118 Solid Parabolic Tetrahedron
!   121 Rigid Bar
!   122 Rigid Element
!   136 Node To Node Translational Spring
!   137 Node To Node Rotational Spring
!   138 Node To Ground Translational Spring
!   139 Node To Ground Rotational Spring
!   141 Node To Node Damper
!   142 Node To Gound Damper
!   151 Node To Node Gap
!   152 Node To Ground Gap
!   161 Lumped Mass
!   171 Axisymetric Linear Shell
!   172 Axisymetric Parabolic Shell
!   181 Constraint
!   191 Plastic Cold Runner
!   192 Plastic Hot Runner
!   193 Plastic Water Line
!   194 Plastic Fountain
!   195 Plastic Baffle
!   196 Plastic Rod Heater
!   201 Linear node-to-node interface
!   202 Linear edge-to-edge interface
!   203 Parabolic edge-to-edge interface
!   204 Linear face-to-face interface
!   208 Parabolic face-to-face interface
!   212 Linear axisymmetric interface
!   213 Parabolic axisymmetric interface
!   221 Linear rigid surface
!   222 Parabolic rigin surface
!   231 Axisymetric linear rigid surface
!   232 Axisymentric parabolic rigid surface
!
! ======================= Dataset Number 2467 =======================================
!Record 1:        FORMAT(8I10)
!                 Field 1       -- group number
!                 Field 2       -- active constraint set no. for group
!                 Field 3       -- active restraint set no. for group
!                 Field 4       -- active load set no. for group
!                 Field 5       -- active dof set no. for group
!                 Field 6       -- active temperature set no. for group
!                 Field 7       -- active contact set no. for group
!                 Field 8       -- number of entities in group
!
!Record 2:        FORMAT(20A2)
!                 Field 1       -- group name
!
!Record 3-N:      FORMAT(8I10)
!                 Field 1       -- entity type code
!                 Field 2       -- entity tag
!                 Field 3       -- entity node leaf id.
!                 Field 4       -- entity component/ ham id.
!                 Field 5       -- entity type code
!                 Field 6       -- entity tag
!                 Field 7       -- entity node leaf id.
!                 Field 8       -- entity component/ ham id.
!
!Repeat record 3 for all entities as defined by record 1, field 8.
!Records 1 thru n are repeated for each group in the model.
!Entity node leaf id. and the component/ ham id. are zero for all
!entities except "reference point", "reference point series"
!and "coordinate system".
!
!          Permanent group entity type codes
!
!    Entity Type Code        Entity Description
!
!           1                coordinate system
!           2                data surface thickness
!           3                force on point
!           4                force on edge
!           5                traction on face
!           6                pressure on face
!           7                nodes
!           8                finite elements
!           9                dof sets, dof entities
!          10                constraint sets, coupled dofs
!          11                constraint sets, mpc equations
!          12                restraint sets, nodal displacements
!          13                restraint sets, nodal temperatures
!          14                load sets, nodal forces
!          15                load sets, nodal temperatures
!          16                load sets, nodal heat sources/sinks
!          17                load sets, face pressures
!          18                load sets, edge pressures
!          19                load sets, face heat fluxes
!          20                load sets, edge heat fluxes
!          21                load sets, face heat convections
!          22                load sets, edge heat convections
!          23                load sets, face heat radiations
!          24                load sets, edge heat radiations
!          25                load sets, element heat generations
!          26                load sets, beam temperatures
!          27                trace lines
!          28                beam force
!          29                beam distributed load
!          30                data surface
!          31                data curve
!          32                displacement on point (restraint)
!          33                displacement on edge (restraint)
!          34                displacement on surface (restraint)
!          35                temperature on point (restraint)
!          36                temperature on edge (restraint) 
!          37                temperature on face (restraint)
!          38                temperature on point (temperature)
!          39                temperature on edge (temperature)
!          40                temperature on face (temperature)
!          41                heat source on point
!          42                heat flux on edge
!          43                convection on edge
!          44                radiation on edge
!          45                heat flux on face
!          46                convection on face
!          47                radiation on face
!          48                geometry contact region
!          49                fe contact region
!          50                contact pair
!          51                kinematic dof on point
!          52                kinematic dof on edge
!          53                kinematic dof on face
!          54                element definition
!          55                anchor node
!          56                edge dependancy mesh definition
!          57                fem point connector
!          58                fem area connector
!          59                vertex
!          60                edge
!          61                face
!          62                region
!          63                wireframe connector
!          64                wireframe curve
!          65                wireframe section
!          66                wireframe region
!          67                reference point
!          68                reference point series
!          69                centerpoint

