      subroutine vol2aux(numf,vf)
      use feminterface, only: reallocate, destroyarrptr, realloc
      use feminterface3D, only: reversemap, getvf
      use femtypes
      use globalvariables3D, only: en, nume, numn, numv, ve, vn, vv
      implicit none
      integer (I4B), intent(out), optional :: numf
      integer (I4B), pointer, optional :: vf(:,:)
!
!------------------------------------------------------------------------------
!    $Revision: 1.10 $
!    $Date: 2015/11/11 17:24:04 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Subroutine computes maps for volumes -> edges, edges -> nodes, neighbours
!  and number of edges. Number of faces and volume -> faces map can be obtained
!  as optional variables.
!
!------------------------------------------------------------------------------
!  Output:
!     nume      number of edges (global)
!     numf      number of faces
!     en        edge -> nodes              allocate(en(2,nume)) (global)
!     ve        volume element -> edges    allocate(ve(6,numv)) (global)
!     vf        volume element -> faces    allocate(vf(4,numv))
!     vv        neighbour information      allocate(vv(4,numv))
!               vv(i,j) is the element adjacent to j at element face i
!               (opposite to the vertex i). If element j has no neighbour at
!               face i,  vv(i,j)=0.
!
!  Internal variables:
      integer (I4B) :: i, j, k, l       ! counters
      integer (I4B) :: elem, elem2      ! element number
      integer (I4B) :: enfree
      integer (I4B) :: le, le2, ln      ! local edge/node
      integer (I4B) :: numsecnode, secnodecand
      integer (I4B) :: nmin, nmax
      integer (I4B), allocatable :: secnode(:,:)
      logical :: destroyed
      type (ARRPTRI), pointer :: ev(:)=>null(), nv(:)=>null()
!  tonode(i,j) represents the three local nodes j=1..3 that are connected to node i
      integer (I4B), parameter :: tonode(4,3)=reshape((/2,1,1,1,3,3,2,2,4,4,4,3/),(/4,3/))
!  edgenum(i,j) represents the local edge numbers for the edge, that connects
!  node i with the node j=tonode(i,:)
      integer (I4B), parameter :: edgenum(4,3)=reshape((/1,1,3,4,3,2,2,5,4,5,6,6/),(/4,3/))
!  opsn(i,le) are the two nodes i=(1,2) opposite to local edge le
      integer (I4B), parameter :: opsn(2,6)=reshape((/3,4,1,4,2,4,2,3,1,3,1,2/),(/2,6/))
!  fen(i,le) is the face defined by local edge le and opposite node i=opsn(:,le)
      integer (I4B), parameter :: fen(2,6)=reshape((/4,3,4,1,4,2,3,2,3,1,2,1/),(/2,6/))
!
!  get the reverse map node -> volumes (nv) from vn
      call reversemap(.true.,vn,nv,nmin,nmax)
!
!  Allocate and initialize arrays for en, fn, ve and vf
      if (associated(en)) deallocate(en)
      if (associated(ve)) deallocate(ve)
      allocate (en(2,2*numv+4), ve(6,numv))
      en = 0
      ve = 0
!  Set next free entry in en and fn
      enfree = 1
!     
!------------------------------------------------------------------------------
!  COMPUTE EN, VE, and NUME:
!  -------------------------
!
      do i = 1,numn  ! global nodes
        allocate(secnode(size(nv(i)%d)/2 + 4,2))
        numsecnode = 0
!  Loop over volume elements to which the node belongs
        do j = 1,size(nv(i)%d)
          elem = nv(i)%d(j)
!  Loop over all 4 nodes of elem to find local position of node i
          do k = 1,4
            if (vn(k,elem) .eq. i) then
              ln = k
              exit
            end if
          end do
!  Determine three edges connected to local node ln
secnodes: do k = 1,3
            secnodecand = vn(tonode(ln,k),elem)
!  If secnodecand is smaller than actual global node, it was treated before ==> cycle
            if(secnodecand .lt. i) cycle
!  Check if secnodecand was written before, if yes write ve map ==> cycle
            do l = 1,numsecnode
              if (secnodecand .eq. secnode(l,1)) then
                ve(abs(edgenum(ln,k)),elem) = sign(secnode(l,2),edgenum(ln,k))
                cycle secnodes
              end if
            end do
!  Resize en, if more space than size needed
            if(enfree .gt. size(en(1,:))) then
              en => reallocate(en,2,2*enfree)
            end if
!  Write the nodes to the next enfree entry in en
            en(:,enfree)=(/i,secnodecand/)
!            ve(edgenum(ln,k),elem) = enfree
!  If edges require sign for orientation later use following line instead
            ve(abs(edgenum(ln,k)),elem) = sign(enfree,edgenum(ln,k))
            numsecnode = numsecnode + 1
            if (numsecnode .gt. size(secnode(:,1))) then
              call realloc(secnode,2*numsecnode,2)
            end if
            secnode(numsecnode,1) = secnodecand
            secnode(numsecnode,2) = enfree
            enfree = enfree + 1
          end do secnodes
        end do  ! volume elements
        deallocate(secnode)
      end do  ! global nodes
!
!  destroy array of pointers nv
      destroyed = destroyarrptr(nv)
!
!  Set number of edges and resize en array
      nume = enfree - 1
      en => reallocate(en,2,nume)
!
!------------------------------------------------------------------------------
!  Print data for edges and nodes
!      print *,"volume -> edge map"
!      do i = 1,numv
!        print "(i4,a,6i5)",i,":  ",ve(:,i)
!      end do
!      print *,"edge -> node map"
!      do i = 1,nume
!        print "(i4,a,2i4)",i,":  ",en(:,i)
!      end do
!------------------------------------------------------------------------------
!  COMPUTE VV:
!  -----------
!
!  Allocate and initialize neighbour array vv
      if (associated(vv)) deallocate(vv)
      allocate(vv(4,numv))
      vv = 0
!
!  get the reverse map edge -> volumes (ev) from ve
      call reversemap(.true.,ve,ev,nmin,nmax)
!     
!  Loop over all global edges
      do i = 1,nume
!  Loop over volume elements the edge belongs to
        do j = 1,size(ev(i)%d)
          elem = ev(i)%d(j)
!  Loop over all edges of elem to find local position le of edge i
          do k = 1,size(ve(:,elem))
            if (ve(k,elem) .eq. i) then
              le = k
              exit
            end if
          end do
          do k = 1,size(ev(i)%d)
            elem2 = ev(i)%d(k)
            if (elem .eq. elem2) cycle
!  Loop over all edges of elem2 to find local position le2 of edge i
            do l = 1,size(ve(:,elem2))
              if (ve(l,elem2) .eq. i) then
                le2 = l
                exit
              end if
            end do
!  If elem and elem2 share the same edge and same opposite global node, they
!  are neighbours.
            if ( vn(opsn(1,le),elem) .eq. vn(opsn(1,le2),elem2) ) then
              vv(fen(1,le),elem) = elem2
              vv(fen(1,le2),elem2) = elem
            else if ( vn(opsn(1,le),elem) .eq. vn(opsn(2,le2),elem2) ) then
              vv(fen(1,le),elem) = elem2
              vv(fen(2,le2),elem2) = elem
            else if ( vn(opsn(2,le),elem) .eq. vn(opsn(1,le2),elem2) ) then
              vv(fen(2,le),elem) = elem2
              vv(fen(1,le2),elem2) = elem
            else if ( vn(opsn(2,le),elem) .eq. vn(opsn(2,le2),elem2) ) then
              vv(fen(2,le),elem) = elem2
              vv(fen(2,le2),elem2) = elem
            end if
          end do
        end do  ! volume elements
      end do  ! global edges
!
!  destroy array of pointers ev
      destroyed = destroyarrptr(ev)
!
!------------------------------------------------------------------------------
!  Print data for neighbours
!      print *,"volume -> volume map (neighbours)"
!      do i = 1,numv
!        print "(i4,a,4i5)",i,":  ",vv(:,i)
!      end do
!
!------------------------------------------------------------------------------
!  COMPUTE VF:
!  -----------
!
      if (present(numf) .or. present(vf)) then
        call getvf(numf,vf)
      end if
!
!
      end subroutine vol2aux
