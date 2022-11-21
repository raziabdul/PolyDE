      subroutine hmesh_refine3D(v_mark,es,ev,e_len,num_refined)
      use feminterface,       only: reallocate, qsortindex
      use feminterface3D,     only: BisectTET, swap23
      use globalvariables3D,  only: nod, vn, ve, vv, vp, vp_temp, en, sn, sbc, dom, numv, numn, nums, nume, nnat
      use femtypes
      implicit none
      integer (I4B)               :: num_refined
      integer (I4B),  allocatable :: v_mark(:)
      real     (DP),      pointer :: e_len(:)
      type(ARRPTRI),      pointer :: ev(:), es(:)
      intent(out)                 :: num_refined
      intent(inout)               :: es, ev, e_len, v_mark
!
!-------------------------------------------------------------------------------
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
!    $Revision: 1.7 $
!    $Date: 2015/11/11 15:41:26 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
! Note:
!        es      Pointer Array which contains surfaces, which are sharing a specific edge or 0,0 in case of an edge not being on the surface
!
! Adaptive 3D mesh-refinement:
!
! Changed:
!           nod      \
!           dom      |
!           vn       |
!           vp       |
!           ve       |
!           vv        >          Resized
!           sn       |
!           sbc      |
!           en       |
!           ev       |
!           es       /
!
! EXPLANATIONS of surface handling according to the Polyeder-Law in 3D
!   IF an Edge lies ON THE SURFACE:
!
!   n3'=n3+n3      <=new number of tetraheda
!   n0'=n0+1       <=new number of nodes
!   n1'=n1+n3+2    <=new number of egdes (*)
!   n1s'=n1s+3     <=new number of the edges on surface
!   n2'=n2+(2*n3+1)<=new number of all new faces (*)
!   n2s'=n2s+2     <=new number of surfaces
!   n2f'=n2f+n3    <=new number of inner faces
!
!   IF an Edge DOES NOT lie on the surface
!   Polyeder-Law-3D for CONVEX Elements:
!
!   n3'=n3+n3      <=new number of tetraeda
!   n0'=n0+1       <=new number of nodes
!   n1'=n1+n3+1    <=new number of egdes (*)
!   n2'=n2+2*n3    <=new number of all new faces (*)
!
!-------------------------------------------------------------------------------
! local variables
!
      integer (I4B)              :: v, vol, mark, max_depth, numv_start, numv_end
      integer (I4B), allocatable :: v_mark_idx(:)
      real     (DP), allocatable :: e_len_max(:)
!__
!
! Definitions:
      numv_start = numv
!_____ 
! Start: 
! 1) Initialize:
!-
      e_len(nume+1:size(e_len)) = 0
      v_mark(numv+1:size(v_mark)) = 0
!__
!  
! 2) Iterate through all elements as long as there are still marked elements left
!    (some elements might have been marked for multiple bisections)
!    1. Sort the Elements in descending order according to the longest edge occuring in them,
!    2. Bisect marked Tetrahedrons starting with those having the longest edges.
!       'max_depth'  gives the highest number of bisections given on an element (depth of refinement)
      max_depth = maxval(v_mark(:))
      do while (max_depth.gt.0)
        allocate(v_mark_idx(numv))
        allocate(e_len_max(numv))
        do v = 1, numv
          e_len_max(v) = maxval(e_len(ve(1:6,v)))
        end do
        call qsortindex(e_len_max,v_mark_idx,numv)
        deallocate(e_len_max)
        do v = 1, numv
          vol = v_mark_idx(v)
          mark = v_mark(vol)
          if ( mark.gt.0 ) call BisectTET(vol,es,ev,e_len,v_mark)
        end do
        deallocate(v_mark_idx)
        max_depth = maxval(v_mark(:))
      end do
      numv_end = numv
      num_refined = numv_end - numv_start
!__
!
! 3) Resize the Mesh to actual size:
      nod => reallocate(nod,3,numn)
      dom => reallocate(dom,numv)
      vn  => reallocate(vn,4,numv)
      vp  => reallocate(vp,numv,nnat)
      if (associated(vp_temp)) vp_temp => reallocate(vp_temp,numv,nnat)
      ve  => reallocate(ve,6,numv)
      vv  => reallocate(vv,4,numv)
      sn  => reallocate(sn,3,nums)
      sbc => reallocate(sbc,nums)
      en  => reallocate(en,2,nume)
      ev  => reallocate(ev,nume)
      es  => reallocate(es,nume)
      e_len => reallocate(e_len,nume)
!
! End.
!___
!      call swap23()
      return
end subroutine hmesh_refine3D



subroutine BisectTET(t_0,es,ev,e_len,v_mark)
      use globalvariables3d,   only: ve
      use feminterface,        only: reallocate
      use feminterface3D,      only: getCriticalPath, BisectHull
      use femtypes
      implicit none
      integer (I4B)               :: t_0
      integer (I4B),  allocatable :: v_mark(:)
      real     (DP),      pointer :: e_len(:)
      type(ARRPTRI),      pointer :: ev(:), es(:)
      intent   (in)               :: t_0
      intent(inout)               :: ev, es, e_len, v_mark
!
!-------------------------------------------------------------------------------
!
! This Routine: 'BisectTET'
!  Performs the mesh refinement
!  Terms - Critical Path (CRITP)
!    and - Longest Edge Propagation Path (LEPP)
!  are user vice versa
!
!-------------------------------------------------------------------------------
! Input:
!           t_0        Tetrahedron to be splitted
! Output:
!           refined and conforming mesh
!
! Local variables:
      integer  (I4B)                    :: i, k, e_0, e_L, subset_size, vol, vol_L, l_6_pos(1)
      integer  (I4B)                    :: last
      real      (DP)                    :: l_6(6)
      integer  (I4B),           pointer :: CP(:)=>null(), n_0(:)=>null(), CP_add(:)=>null(), n_0_add(:)=>null()
!
!__
! Start:
! 1) Get Critical Path:
      call getCriticalPath(t_0,CP,n_0,ev,e_len)
      last = size(CP)
!__
!  
! 2) Work through all edges of the Critical Path and bisect hulls:
!    e_0: The edge subject to bisection.
!    n_0: Global node number of "upright node" of e_0
!         Elements/ Faces/ Edges connected to this node maintain their global names/ numbers
!
!    Before bisection, we check weather e_0 really represents the longest edge of the bisection hull.
!    Check in each surrounding tetrahedron for it's longest edge
!    a) In case of e_0 = 0, cycle because this edge has already been bisected in the process of bisection of former hulls
!    b) Find the longest edge of the Hull of e_0, name it 'e_L'
!    c) If e_L turned out to be e_0, bisect the hull
!    d) If e_L is not e_0, compute the critical path for e_L, name it 'CP_add'
!    e) Iterate through CP_add to check weather it contains edges which are shorter
!       than any edge occuring in CP, if so, cut the CP_add vector to prevent insertion
!       of shorter edges than already contained in CP, mark insertion position in CP_add, name it 'subset_size'
!    f) Insert CP_add into CP, continue loop from the last entry upwards
      do while(CP(last).ne.CP(1))
!- 
        e_0 = CP(last)
!- a)
        if (e_0.eq.0) then
          CP  => reallocate(CP,(size(CP)-1))
          n_0 => reallocate(n_0,(size(n_0)-1))
          last = size(CP)
          cycle
        end if
!- b)
        do i = 1 , size(ev(CP(last))%d)
          vol = ev(CP(last))%d(i)
          l_6(:) = e_len(ve(1:6,vol))
          l_6_pos = maxloc(l_6)
          e_L = ve(l_6_pos(1),vol)
          if (e_len(e_0).lt.e_len(e_L)) then
            e_0 = e_L
            vol_L = vol
          end if
        end do
!- c)
        if (e_0.eq.CP(last)) then
          call BisectHull(CP(last),es,ev,e_len,v_mark,n_0(last))
          CP  => reallocate(CP,(size(CP)-1))
          n_0 => reallocate(n_0,(size(n_0)-1))
        else
!- d)
          call getCriticalPath(vol_L,CP_add,n_0_add,ev,e_len)
          subset_size = size(CP_add)
!- e)
          do i = 1, size(CP_add)
            do k = 1, size(CP)
              if(CP(k).eq.CP_add(size(CP_add)-i+1)) then
                CP(k) = 0
                subset_size = i
              end if
            end do
          end do
!- f)
          CP  => reallocate(CP,(size(CP)+subset_size))
          n_0 => reallocate(n_0,(size(n_0)+subset_size))
          CP(size(CP)-subset_size+1:size(CP)) = CP_add(1:subset_size)
          n_0(size(n_0)-subset_size+1:size(n_0)) = n_0_add(1:subset_size)
          deallocate(CP_add,n_0_add)
          nullify(CP_add,n_0_add)
        end if
        last = size(CP)
      end do
!__
! 3) Bisect last edge, no more longer edges to take care of:
      call BisectHull(CP(1),es,ev,e_len,v_mark,n_0(1))
!__
!
! Release Memory:
      if (associated(n_0)) deallocate(n_0)
      if (associated(CP))  deallocate(CP)
!
! End.
!___
      return
end subroutine BisectTET


subroutine getCriticalPath(t_0,CP,n_0,ev,e_len)
      use feminterface,              only: reallocate, palloc
      use feminterface3D,            only: SortEdges
      use globalvariables3D,         only: ve, en, vn
      use femtypes
      implicit none
      integer  (I4B)                    :: t_0
      integer  (I4B),           pointer :: CP(:), n_0(:)
      real     (DP),            pointer :: e_len(:)
      type(ARRPTRI),            pointer :: ev(:)
      intent    (in)                    :: t_0
      intent   (out)                    :: CP, n_0
      intent (inout)                    :: ev, e_len
!
!-------------------------------------------------------------------------------
!
! This Routine: 'getCriticalPath'
! Computes the critical path of edges-subdivisions needed to perform in order to split Tetrahedron t_0 at it's longest edge.
!
!-------------------------------------------------------------------------------
! Input:
!           t_0                  Global Tetrahedron Number. Tetrahedron to be splitted.
! Output:
!           CP                   Edges of Critical path: Starting from the bottom, all edges have to be bisected in order to bisect CP(1)
!           n_0                  Global Node Number. Created Tetrahedrons at the side of n_0 will be assigned ...
!                                ... To the name of the Father tetrahedrons t_0 if n_0 is present
!           numof_newelems       Numbers of tetrahedron splitted in process of bissection of an edge from CP(edges)
!           surfedge             Logical representation of edge-surface location
!
!-------------------------------------------------------------------------------
!
! EXPLANATIONS of Mesh Extension according to the Polyeder-Law in 3D:
!
! 1)IF an Edge lies ON THE SURFACE:
!   
!   n3'=n3+n3      <=new number of tetraheda
!   n0'=n0+1       <=new number of nodes
!   n1'=n1+n3+2    <=new number of egdes (*)
!   n1s'=n1s+3     <=new number of the edges on surface
!   n2'=n2+(2*n3+1)<=new number of all new faces (*)
!   n2s'=n2s+2     <=new number of surfaces
!   n2f'=n2f+n3    <=new number of inner faces
!   
! 2)IF an Edge DOES NOT lie on the surface (Polyeder-Law-3D for KONVEX Elements:)
!   
!   n3'=n3+n3      <=new number of tetraeda
!   n0'=n0+1       <=new number of nodes
!   n1'=n1+n3+1    <=new number of egdes (*)
!   n2'=n2+2*n3    <=new number of all new faces (*)
!-------------------------------------------------------------------------------
! local variables
!
      integer (I4B)              :: i, v
      integer (I4B)              :: vol, hullsize, e_0, e_L, e_L_loc, e_0_loc
      integer (I4B)              :: CP_size_old, CP_size_new, nume_nested
      integer (I4B),     pointer :: CP_temp(:)=>null(), n_0_temp(:)=>null()
      integer (I4B),   parameter :: e2oe(6) = (/6,4,5,2,3,1/)
      integer (I4B),   parameter :: e2cf(2,6) = reshape((/3,4,1,4,2,4,3,2,1,3,1,2/),(/2,6/))
      integer (I4B),   parameter :: t_02t_p(6,6) = reshape( (/  -1,2,1,1,2,0,  &
                                                             &  2,-1,3,0,2,3,  &
                                                             &  1,3,-1,1,0,3,  &
                                                             &  1,0,1,-1,4,4,  &
                                                             &  2,2,0,4,-1,4,  &
                                                             &  0,3,3,4,4,-1   /),(/6,6/) )

!__
! Checks and preparations:
! Determine the longest Edge of t_0, name it 'e_0'
      e_0 = ve(1,t_0)
      do i = 2, 6
        if (e_len(ve(i,t_0)).gt.e_len(e_0)) then
          e_0 = ve(i,t_0)
        end if
      end do
!__
!
! Allocate Memory:
      call palloc(CP,1)
      call palloc(CP_temp,1)
      call palloc(n_0,1)
      call palloc(n_0_temp,1)
      CP(1) = e_0
      CP_temp = CP
      n_0(1) = en(1,e_0)
      n_0_temp = n_0
!_____
! Start:
!
! 1) Get Critical Path
      do while (e_0.gt.0)
        hullsize = size(ev(e_0)%d)
        CP_size_old = size(CP)
        do v = 1, hullsize
          vol = ev(e_0)%d(v)
       ! Determine the longest Edge of 'vol'
          e_L = ve(1,vol)
          do i = 2, 6
            if (e_len(ve(i,vol)).gt.e_len(e_L)) then
              e_L = ve(i,vol)
           end if
          end do
          if (e_L.ne.e_0) then ! A tetrahedron has been found which has an edge longer that e_0 it's longest edge will be added to the CP if it is not in yet
        ! Determine local numbers of e_0 and e_L:
            do e_L_loc = 1, 6
              if (ve(e_L_loc,vol).eq.e_L) exit
            end do
            do e_0_loc = 1, 6
              if (ve(e_0_loc,vol).eq.e_0) exit
            end do
            if (.not.any(CP(:).eq.e_L)) then
              if (.not.any(CP_temp(:).eq.e_L)) then
                CP => reallocate(CP,size(CP)+1)
                CP(size(CP)) = e_L
                n_0 => reallocate(n_0,size(CP))
                if (t_02t_p(e_0_loc,e_L_loc).gt.0) then ! e_0 and e_L are neighboring edges, bisection of the element must be performed fulfilling condition for n_0
                  n_0(size(CP)) = vn(t_02t_p(e_0_loc,e_L_loc),vol)
                else if (t_02t_p(e_0_loc,e_L_loc).eq.0) then ! e_0 and e_L are oposite edges, bisection of the element can be performed without fulfilling condition for n_0
                  n_0(size(CP)) = en(1,e_L)
                else ! e_0 and e_L are equal, subdivision of e_0 is sufficient
                  n_0(size(CP)) = 0
                end if
              end if
            end if
          end if
        end do
        CP_size_new = size(CP)
        nume_nested = CP_size_new - CP_size_old
        if (nume_nested.gt.0) then
          CP_temp => reallocate(CP_temp,size(CP_temp)+nume_nested)
          CP_temp(size(CP_temp)-nume_nested+1:size(CP_temp)) = CP(CP_size_old+1:CP_size_new)
          n_0_temp => reallocate(n_0_temp,size(CP_temp))
          n_0_temp(size(CP_temp)-nume_nested+1:size(CP_temp)) = n_0(CP_size_old+1:CP_size_new)
        end if
        if (size(CP).gt.1) then
          CP(1:size(CP)-1) = CP(2:size(CP))
          CP => reallocate(CP,size(CP)-1)
        else
          CP(1) = 0
        end if
        e_0 = CP(1)
      end do
!__
!
! 2) Sort the edges of the Critical Path according to their length (longest last):
      if (size(CP).gt.1) call SortEdges(CP_temp,n_0_temp)
!__
!
! Release Memory:
!
      deallocate(CP)
      deallocate(n_0)
      n_0 => n_0_temp
      CP => CP_temp
      nullify(CP_temp)
      nullify(n_0_temp)
!
! End.
!___
      return
end subroutine getCriticalPath


subroutine BisectHull(e_0,es,ev,e_len,v_mark,n_0_in)
      use feminterface,         only: reallocate, realloc
      use feminterface3d,       only: EdgeLength, updateaux
      use globalvariables3D,    only: nod, vn, ve, vv, vp, vp_temp, en, sn, sbc, numv, numn, nums, dom, nnat
      use matconstants
      use femtypes
      implicit none
      integer  (I4B)               :: e_0
      integer  (I4B),  allocatable :: v_mark(:)
      real      (DP),      pointer :: e_len(:)
      type (ARRPTRI),      pointer :: ev(:), es(:)
      integer  (I4B),     optional :: n_0_in
      intent    (in)               :: e_0, n_0_in
      intent (inout)               :: ev, es, e_len, v_mark
!
!-------------------------------------------------------------------------------
!
! This Routine: 'BisectHull'
! Refinement of 1 Edge and surrounding tetrahedra plus conformity
!
!-------------------------------------------------------------------------------
!                      
! Input:
!           e_0        Global Edge Number. Edge to be splitted.
!           n_0        Global Node Number. Created Tetrahedrons at the side of n_0 will be assigned to the name of the Father tetrahedrons t_0 if n_0 is present
!           es         Global numbers of Surfaces at the edge (if given, othervise 0,0)
!           ev         Global numbers of volumes at the edge
! Output:
!           Refined and conforming mesh (updateaux assures mesh conformity)
!           
!
!-------------------------------------------------------------------------------
!
! Local variables:
!
      integer (I4B)              :: m, n, s, v
      integer (I4B)              :: numv_add
      integer (I4B)              :: vol, n_t, t_t, n_0, nb_t, nb_0
      integer (I4B)              :: surf
      integer (I4B), allocatable :: vols_inchange(:), surfs_inchange(:), et_nodes(:), vn0nt(:,:)
      integer (I4B),   parameter :: e2cf(2,6) = reshape((/3,4,1,4,2,4,3,2,1,3,1,2/),(/2,6/))
      real (DP) :: xyz(3)
! _____
! Start:
!
      numv_add = size(ev(e_0)%d)
!__
! 1) Determine n_0 and n_t:
      if (present(n_0_in).and.(n_0_in.gt.0)) then
        n_0 = n_0_in
      else
        n_0 = en(1,e_0)
      end if
      do m = 1, 2 ! n_0 and n_t are the two end points of the mother edge.
        if (en(m,e_0).ne.n_0) n_t = en(m,e_0)
      end do
   ! Determine local numbers on n_0 and n_t and store them: vn0nt(1:3,vol) = (/ vol(Elem.Numb.) , n_0(Local) , n_t(Local) /)
      allocate(vn0nt(3,numv_add))
      do v = 1, numv_add
        vn0nt(1,v) = ev(e_0)%d(v)
        do n = 1, 4
          if (vn(n,vn0nt(1,v)).eq.n_0) vn0nt(2,v) = n
          if (vn(n,vn0nt(1,v)).eq.n_t) vn0nt(3,v) = n
        end do
      end do
!__
!
! 2) Create New Node:
!
!    Insert it on e_0. The Position for the new Node can be chosen in following possible manners:
!
! - Version 1: only nodes of edge surrounding
!      xyz=xyz - (nod(1:3,en(1,e_0))+nod(1:3,en(2,e_0)))*size(ev(e_0)%d)
!      xyz=xyz/( 2 * size( ev(e_0)%d) )
!
! - Version 2: Sum of all nodes
!      xyz=xyz/( 4 * size( ev(e_0)%d) )
!
! - Version 3: average between version 1 and version 2
!      xyz=xyz - (nod(1:3,en(1,e_0))+nod(1:3,en(2,e_0)))/2*size(ev(e_0)%d)
!      xyz=xyz/( 3 * size( ev(e_0)%d) )
!__
!   ALTERNATIVELY (Less desirable): insertion at midpoint of e_0:
!      nod(1:3,numn) = (/(nod(1,en(1,e_0)) + nod(1,en(2,e_0))) ,        &
!     &                 (nod(2,en(1,e_0)) + nod(2,en(2,e_0))) ,         &
!     &                 (nod(3,en(1,e_0)) + nod(3,en(2,e_0))) /)*0.5_DP
!__
!   USED NOW: Version 3 - Average between version 1 and version 2
! a) Start of hull's mass center calculation
! b) Projection onto e_0
! c) Insertion of projected Node
      numn = numn + 1
      if (size(nod,2).lt.numn) nod => reallocate(nod,3,numn*2)
!-  a)
      xyz=0._DP
      do n=1,4
        do v=1,size(ev(e_0)%d)
          xyz=xyz+nod( 1:3 , vn(n , ev(e_0)%d(v) ) )
        end do
      end do
      xyz=(1._DP+0.25_DP)*xyz - (2._DP*0.25_DP)*(nod(1:3,en(1,e_0))+nod(1:3,en(2,e_0)))*size(ev(e_0)%d)
      xyz=xyz/( 4 * size( ev(e_0)%d) )
!-  b)
      nod(1:3,numn) = (/(nod(1,en(2,e_0)) - nod(1,en(1,e_0))) ,         &
     &                  (nod(2,en(2,e_0)) - nod(2,en(1,e_0))) ,         &
     &                  (nod(3,en(2,e_0)) - nod(3,en(1,e_0))) /)
      xyz(1:3) = (/(xyz(1) - nod(1,en(1,e_0))) ,                        &
     &            (xyz(2) - nod(2,en(1,e_0))) ,                         &
     &            (xyz(3) - nod(3,en(1,e_0))) /)
!-  c)
      nod(1:3,numn)=dot_product(xyz,nod(1:3,numn))/dot_product(nod(1:3,numn),nod(1:3,numn))*nod(1:3,numn)+nod(1:3,en(1,e_0))
!__
! 3) Adjust 'en' for e_0:
!    'en' is treated separately to assure that the right edge assignment for e_0 is made.
!    Other edges will be inserted in the process of 'updateaux'
      allocate(et_nodes(2))
      if (any(en(:,e_0).eq.n_0)) then
        do n = 1,2
          if (en(n,e_0).eq.n_0) exit
        end do
        if (n.eq.1) en(1:2,e_0) = (/ en(1,e_0) , numn /)
        if (n.eq.2) en(1:2,e_0) = (/ en(2,e_0) , numn /)
        et_nodes(1:2) = (/ n_t , numn /)
        e_len(e_0) = EdgeLength(en(1,e_0),en(2,e_0))
      else
        print*, '**** ERROR: Error in 3D-h-Adaptation (Tetrahedron-Bisection of Edge ', e_0 ,') - Node selected as n_0 is not part of the Edge to be bisected:'
        print*, 'Nodes of the Edge ', e_0, ':', en(:,e_0)
        print*, 'n_0 :', n_0
        pause
        stop
      end if
!__
!
! 4) Adjust vv:
!    'vv' is treated here only for elements outside the Hull.
!    Treatment for elements inside the Hull is done in 'updateaux'
!    Resize the arrays if needed
      if (size(vv,2).lt.numv+numv_add) then
        dom => reallocate(dom,numv*2)
        vn  => reallocate(vn,4,numv*2)
        vp  => reallocate(vp,numv*2,nnat)
        if (associated(vp_temp)) vp_temp  => reallocate(vp_temp,numv*2,nnat)
        ve  => reallocate(ve,6,numv*2)
        vv  => reallocate(vv,4,numv*2)
        call realloc(v_mark,numv*2)
        v_mark(numv+1:size(v_mark)) = 0
      end if
!__
!
      do v = 1, numv_add
        vol = vn0nt(1,v)
        t_t = numv+v
        nb_t = vv(vn0nt(2,v),vol)
        nb_0 = vv(vn0nt(3,v),vol)
        if (nb_t.gt.0) then ! Only update Neighbour information if there is one
          do n = 1,4
            if (vv(n,nb_t).eq.vol) vv(n,nb_t) = t_t
          end do
        end if
       ! Since the new node will carry the highest number, we can be sure, that the face oposite to it will be of local number 4 (for t_0 and t_t)
        vv(4,vol) = nb_0
        vv(4,t_t) = nb_t
        if (n.gt.5) then
          print*, '**** ERROR: Error in 3D-h-Adaptation (Tetrahedron-Bisection of Tetrahedron ', vol ,') - Fail of Update of neighbour:', nb_t
          pause
          stop
        end if
      end do
!__
!
! 5) Adjust vn:
      allocate(vols_inchange(2*numv_add))
      do v = 1, numv_add
        vol = vn0nt(1,v)
        if (any(vn(:,vol).eq.n_0)) then
          numv = numv + 1
          dom(numv)  = dom(vol)
          vp(numv,:) = vp(vol,:)
          if (associated(vp_temp)) vp_temp(numv,:) = vp_temp(vol,:)
          v_mark(vol)  = max(0,(v_mark(vol)-1))
          v_mark(numv) = v_mark(vol)
!-
          vn(1:4,numv) = vn(1:4,vol)
          vols_inchange(2*v-1) = vol
          vols_inchange(2*v)   = numv

          if (vn0nt(3,v).eq.1)  vn(1:4,vol) = (/ vn(2,vol) , vn(3,vol) , vn(4,vol) , numn /)
          if (vn0nt(3,v).eq.2)  vn(1:4,vol) = (/ vn(1,vol) , vn(3,vol) , vn(4,vol) , numn /)
          if (vn0nt(3,v).eq.3)  vn(1:4,vol) = (/ vn(1,vol) , vn(2,vol) , vn(4,vol) , numn /)
          if (vn0nt(3,v).eq.4)  vn(1:4,vol) = (/ vn(1,vol) , vn(2,vol) , vn(3,vol) , numn /)

          if (vn0nt(2,v).eq.1)  vn(1:4,numv) = (/ vn(2,numv) , vn(3,numv) , vn(4,numv) , numn /)
          if (vn0nt(2,v).eq.2)  vn(1:4,numv) = (/ vn(1,numv) , vn(3,numv) , vn(4,numv) , numn /)
          if (vn0nt(2,v).eq.3)  vn(1:4,numv) = (/ vn(1,numv) , vn(2,numv) , vn(4,numv) , numn /)
          if (vn0nt(2,v).eq.4)  vn(1:4,numv) = (/ vn(1,numv) , vn(2,numv) , vn(3,numv) , numn /)
        else
          print*, '**** ERROR: Error in 3D-h-Adaptation (Tetrahedron-Bisection of Tetrahedron ', vol ,') - Node selected as n_0 is not part of Tetrahedron to be bisected:'
          print*, 'Nodes of Tetrahedron ', vol, ':', vn(:,vol)
          print*, 'n_0 :', n_0
          pause
          stop
        end if
      end do
!__
!
! 6) Adjust sn:
      if (associated(es(e_0)%d)) then
        if (size(sn,2).lt.nums+2) then
          sn  => reallocate(sn,3,nums*2)
          sbc => reallocate(sbc,nums*2)
        end if
        allocate(surfs_inchange(4))
        do s = 1,2
          surf = es(e_0)%d(s)
          if (any(sn(:,surf).eq.n_0)) then
            nums = nums + 1
            sbc(nums) = sbc(surf)
            sn(1:3,nums) = sn(1:3,surf)

            surfs_inchange(2*s-1) = surf
            surfs_inchange(2*s)   = nums

            do n = 1,3
              if (sn(n,surf).eq.n_0) exit
            end do
            do m = 1,3
              if (sn(m,surf).eq.n_t) exit
            end do

            if (m.eq.1)  sn(1:3,surf) = (/ sn(2,surf) , sn(3,surf) , numn /)
            if (m.eq.2)  sn(1:3,surf) = (/ sn(1,surf) , sn(3,surf) , numn /)
            if (m.eq.3)  sn(1:3,surf) = (/ sn(1,surf) , sn(2,surf) , numn /)

            if (n.eq.1)  sn(1:3,nums) = (/ sn(2,nums) , sn(3,nums) , numn /)
            if (n.eq.2)  sn(1:3,nums) = (/ sn(1,nums) , sn(3,nums) , numn /)
            if (n.eq.3)  sn(1:3,nums) = (/ sn(1,nums) , sn(2,nums) , numn /)
          else
            print*, '**** ERROR: Error in 3D-h-Adaptation (Tetrahedron-Bisection of Surface ', surf ,') - Node selected as n_0 is not part of the Surface to be bisected:'
            print*, 'Nodes of Surface ', surf, ':', sn(:,surf)
            print*, 'n_0 :', n_0
            pause
            stop
          end if
        end do
      end if
!__
!  
! 7) Conform the mesh
!
      if (associated(es(e_0)%d)) then
        call updateaux(vols_inchange,e_0,es,ev,e_len,et_nodes,surfs_inchange)
      else
        call updateaux(vols_inchange,e_0,es,ev,e_len,et_nodes)
      end if 
!__
!
! Release Memory:
      if (associated(es(e_0)%d)) then
        deallocate(surfs_inchange)
      end if
      deallocate(vols_inchange)
      deallocate(et_nodes)
      deallocate(vn0nt)
!__
!
      return
end subroutine BisectHull



subroutine updateaux(vols,e_0,es,ev,e_len,et_nodes,surfs)
      use feminterface,              only: reallocate, palloc
      use feminterface3d,            only: EdgeLength
      use globalvariables3D,         only: ve, vn, en, nume, en, sn, vv
      use femtypes
      implicit none
      integer  (I4B)                    :: e_0
      integer  (I4B)                    :: vols(:), et_nodes(:)
      integer  (I4B)         , optional :: surfs(:)
      real      (DP),           pointer :: e_len(:)
      type (ARRPTRI),           pointer :: ev(:), es(:)
      intent    (in)                    :: vols, e_0, et_nodes, surfs
      intent (inout)                    :: es, ev, e_len
!
!-------------------------------------------------------------------------------
!
! This Routine: 'updateaux'
! Conform the mesh:
! Pass to the routine the volumes involved, the surfaces, ...
! ... the old edge e_0 as well as the new edge e_t (the endpoints of which are stored in 'et_nodes').
! vv will be updated around e_0 and e_t
!
!-------------------------------------------------------------------------------
!
!  Input:
!            vols     : Elements(Volumes) which have been worked on in the process of Hull-Bissection (all t_0 and t_t)
!            e_0      : The number of the common edge of all elements 'vols'
!            es       : Global numbers of Surfaces at the edge (if given, othervise 0,0)
!            ev       : Global numbers of volumes at the edge
!            et_nodes : Endpoints of the edge e_t which will be created and which suppose to accure as continuition of edge e_0
!            surfs    : If the hulledge is a surface-edge, the global numbers of surfaces are delivered here which are involved into subdivision
!  Output:
!            conform mesh
!
!  Update:
!            ve
!            en
!            vv
!
!-------------------------------------------------------------------------------
!
!  local variables:
      integer (I4B)              :: e, f, i, k, m, n, s, v
      integer (I4B)              :: edg, e1, e2, ed, vol, vol_2, n_0, n_t, t_0, t_t, e_t
      integer (I4B)              :: hullsize
      integer (I4B)              :: en_candi(2), edges_intern(2)
      integer (I4B),     pointer :: edges_involved(:)=>null()
      !  e_perm(i,e) are the two local nodes 'i' (where i = 1,2) which give the endpoints of an local edge 'e'
      integer (I4B),   parameter :: e_perm(2,6) = reshape((/1,2,2,3,1,3,1,4,2,4,3,4/),(/2,6/))
      !  opsn(i,le) are the two nodes i=(1,2) opposite to local edge le
      integer (I4B),   parameter :: opsn(2,6) = reshape((/3,4,1,4,2,4,2,3,1,3,1,2/),(/2,6/))
      !  fen(i,le) is the face defined by local edge le and opposite node i=opsn(:,le)
      integer (I4B),   parameter :: fen(2,6) = reshape((/4,3,4,1,4,2,3,2,3,1,2,1/),(/2,6/))
      !  fsn(3,f) is the node-set defined by local face 'f'
      integer (I4B),   parameter :: fsn(3,4) = reshape((/2,3,4,1,3,4,1,2,4,1,2,3/),(/3,4/))
      logical                    :: edge_found

      edges_intern(:) = 0
      en_candi(:) = 0
      hullsize = size(vols)/2
!
!______
! Start:
!
! 1) Collect information about which edges are involved into hull subdivision so far, non-existing edges will be added later
!
      call palloc(edges_involved,6)
      edges_involved = ve(1:6,vols(1))
      do v = 2, hullsize
        do e = 1, 6
          if(.not.any(edges_involved(:).eq.ve(e,vols(2*v-1)))) then
            edges_involved => reallocate(edges_involved,size(edges_involved)+1)
            edges_involved(size(edges_involved)) = ve(e,vols(2*v-1))
          end if
        end do
      end do
!__
!
! 2) en, ve: Add non-existing edges, restore confirmity of en and ve:
!
      do v = 1, size(vols)
        vol = vols(v)
        do e = 1,6
        ! Get two points, where an edge update is possibly necessary
          en_candi(1:2) = (/ vn(e_perm(1,e),vol) , vn(e_perm(2,e),vol) /)
          edge_found = .false.
          do i = 1, size(edges_involved)
          ! If the edge does already exist, use it
            if ( (en(1,edges_involved(i)).eq.en_candi(1)).and.(en(2,edges_involved(i)).eq.en_candi(2)) ) then
              ve(e,vol) = edges_involved(i)
              edge_found = .true.
              exit
            end if
          end do
          ! If the edge between the two points is not existing yet, create a new edge, store it in the list and assign it to ve(e,vol)
          if (.not.edge_found) then
            nume = nume + 1
            if (size(en,2).lt.nume) then
              en  => reallocate(en,2,nume*2)
              ev  => reallocate(ev,nume*2)
              es  => reallocate(es,nume*2)
              e_len => reallocate(e_len,nume*2)
              e_len(nume+1:size(e_len)) = 0
            end if
            if ( (et_nodes(1).eq.en_candi(1)).and.(et_nodes(2).eq.en_candi(2)) ) then
              ! Here we have found the edge e_t, which represents the continuation of e_0
              ! We need it's global number for vv - update.
              ! Store: edges_intern(1) = e_0, edges_intern(2) = e_t
              edges_intern(1:2) = (/ e_0 , nume /)
            end if

            en(1:2,nume) = en_candi(1:2)
            ve(e,vol) = nume
            edges_involved => reallocate(edges_involved,size(edges_involved)+1)
            edges_involved(size(edges_involved)) = nume

            e_len(nume) = EdgeLength(en(1,nume),en(2,nume))
          end if
        end do
      end do
!__
!
! 3) ev: Restore Confirmity of ev
!
! Be n = last inserted node and node 4 of all tetrahedrons inside Hull (due to it's high number) --> n = et_nodes(2)
!    m = former node n_t    -----------------------------------------------------------------------> m = et_nodes(1)
! Edges Types:
! Type 1 :  contains n, contains m : 1 Edge   has this type. Edge e_t has been newly created, it is shared by all newly created elements
! Type 2 :  contains n,   not m    : 2 Edges have this type. They are connected to the node 'n' and are shared by 4 elements, two old ones, two new ones
! Type 3 :    not n   , contains m : 2 Edges have this type. They are shared by neighbours of t_t and new elements (t_t)
! Type 4 :    not n   ,   not m    : 1 Edge   has this type. It is oposite to e_t and element t_t needs to be added to it's list beside t_0
      e_t = 0
      do v = 1, hullsize
        vol = vols(2*v)
        do e = 1,6 ! This volume needs to be part of 6 ev
          edg = ve(e,vol)
          if ((en(1,edg).eq.et_nodes(1)).and.(en(2,edg).eq.et_nodes(2))) then           ! Type 1
            e_t = edg
            if (associated(ev(edg)%d)) cycle
            call palloc(ev(edg)%d,hullsize)
            do i = 1, hullsize
              ev(edg)%d(i) = vols(2*i)
            end do
            ! If e_0 is not a surfedge, e_t is not a surfedge either:
            nullify(es(edg)%d)
          else if (.not.any(en(:,edg).eq.et_nodes(1)).and.((any(en(:,edg).eq.et_nodes(2))))) then      ! Type 2
            if (associated(ev(edg)%d)) then ! ev for this edge is already allocated, two entries are filled in
              if (ev(edg)%d(3).eq.0) then
                ev(edg)%d(3)   = vol
                ev(edg)%d(4)   = vols(2*v-1)
                ! If entries 3 and 4 are written for this edge, it can not be an surface edge, it has 4 surrounding volumes:
                nullify(es(edg)%d)
              end if
            else ! ev needs to be allocated, allocate and place the first two entries
              call palloc(ev(edg)%d,4)
              ev(edg)%d(1)   = vol
              ev(edg)%d(2)   = vols(2*v-1)
              ev(edg)%d(3:4) = 0
            end if
          else if (any(en(:,edg).eq.et_nodes(1)).and.(.not.(any(en(:,edg).eq.et_nodes(2))))) then      ! Type 3
            do i = 1, size(ev(edg)%d) ! Replace t_0 by t_t (update)
              if (ev(edg)%d(i).eq.vols(2*v-1)) ev(edg)%d(i) = vol
            end do
          else if (.not.any(en(:,edg).eq.et_nodes(1)).and.(.not.(any(en(:,edg).eq.et_nodes(2))))) then      ! Type 4
            ev(edg)%d => reallocate(ev(edg)%d,size(ev(edg)%d)+1) ! Add t_t (update)
            ev(edg)%d(size(ev(edg)%d)) = vol
          else
            print*, 'ERROR in ev-Update - Edge Type undefined'
          end if
        end do
      end do
!__
!
! 4) es:                    S U R F A C E  L O C A T I O N
! Be n = last inserted node and node 4 of all tetrahedrons inside Hull (due to it's high number) --> n = et_nodes(2)
!    m = former node n_t    -----------------------------------------------------------------------> m = et_nodes(1)
! Edges Types:
! Type 1 : contains n, contains m : 1 Edge   has this type. Edge e_t has been created, it is shared by all created elements
! Type 2 : contains n,   not m    : 2 Edges have this type. They are connected to the node 'n' and are shared by 2 surfaces, old one and new one
! Type 3 :   not n   , contains m : 2 Edges have this type. They are connected to the node 'm' and are shared by neighbour surfaces and new ones (s_0 needs to be replaced by new surface numbers here)
      if (present(surfs).and.associated(es(e_0)%d)) then ! if(surfedge), doublecheck
        ! Type 1
        if (.not.associated(es(e_t)%d)) then
          call palloc(es(e_t)%d,2)
          es(e_t)%d(1:2) = (/ surfs(2) , surfs(4) /)
        end if
        do e = 1, size(edges_involved)
          ! Type 2
          if (((en(1,edges_involved(e)).eq.sn(1,surfs(2))) .and.    &
     &         (en(2,edges_involved(e)).eq.sn(3,surfs(2)))).and.    &
     &         (edges_involved(e).ne.e_t)) then ! Edge, connected to node 'n', but not e_t is found (check for sn(1->3,surf), surf(2) <-- inserted first)
            call palloc(es(edges_involved(e))%d,2)
            es(edges_involved(e))%d(1:2) = (/ surfs(1) , surfs(2) /)
            if (size(ev(edges_involved(e))%d).gt.2) then
              if (ev(edges_involved(e))%d(3).eq.0) ev(edges_involved(e))%d => reallocate(ev(edges_involved(e))%d,2)
            end if
          end if
          if (((en(1,edges_involved(e)).eq.sn(2,surfs(2))) .and.    &
     &         (en(2,edges_involved(e)).eq.sn(3,surfs(2)))).and.    &
     &         (edges_involved(e).ne.e_t)) then ! Edge, connected to node 'n', but not e_t is found (check for sn(2->3,surf), surf(2))
            call palloc(es(edges_involved(e))%d,2)
            es(edges_involved(e))%d(1:2) = (/ surfs(1) , surfs(2) /)
            if (size(ev(edges_involved(e))%d).gt.2) then
              if (ev(edges_involved(e))%d(3).eq.0) ev(edges_involved(e))%d => reallocate(ev(edges_involved(e))%d,2)
            end if
          end if
          if (((en(1,edges_involved(e)).eq.sn(1,surfs(4))) .and.    &
     &         (en(2,edges_involved(e)).eq.sn(3,surfs(4)))).and.    &
     &         (edges_involved(e).ne.e_t)) then ! Edge, connected to node 'n', but not e_t is found (check for sn(1->3,surf), surf(4))
            call palloc(es(edges_involved(e))%d,2)
            es(edges_involved(e))%d(1:2) = (/ surfs(3) , surfs(4) /)
            if (size(ev(edges_involved(e))%d).gt.2) then
              if (ev(edges_involved(e))%d(3).eq.0) ev(edges_involved(e))%d => reallocate(ev(edges_involved(e))%d,2)
            end if
          end if
          if (((en(1,edges_involved(e)).eq.sn(2,surfs(4))) .and.    &
     &         (en(2,edges_involved(e)).eq.sn(3,surfs(4)))).and.    &
     &         (edges_involved(e).ne.e_t)) then ! Edge, connected to node 'n', but not e_t is found (check for sn(2->3,surf), surf(4))
            call palloc(es(edges_involved(e))%d,2)
            es(edges_involved(e))%d(1:2) = (/ surfs(3) , surfs(4) /)
            if (size(ev(edges_involved(e))%d).gt.2) then
              if (ev(edges_involved(e))%d(3).eq.0) ev(edges_involved(e))%d => reallocate(ev(edges_involved(e))%d,2)
            end if
          end if
          ! Type 3
          if ((en(1,edges_involved(e)).eq.sn(1,surfs(2))).and.(en(2,edges_involved(e)).eq.sn(2,surfs(2)))) then ! Edge, represented by the nodes 1 and 2 (all except new inserted node) of new surface is found sn
            do s = 1,2
              if (es(edges_involved(e))%d(s).eq.surfs(1)) es(edges_involved(e))%d(s) = surfs(2) ! Update from s_0 to s_t for neighbour surfaces.
            end do
          end if
          if ((en(1,edges_involved(e)).eq.sn(1,surfs(4))).and.(en(2,edges_involved(e)).eq.sn(2,surfs(4)))) then ! Edge, represented by the nodes 1 and 2 (all except new inserted node) of new surface is found sn
            do s = 1,2
              if (es(edges_involved(e))%d(s).eq.surfs(3)) es(edges_involved(e))%d(s) = surfs(4) ! Update from s_0 to s_t for neighbour surfaces.
            end do
          end if
        end do
      end if
!__
!
! 5) vv:
      ! The old and the new Elements are neighbours to each other:
      do v = 1, hullsize
      n_0 = en(1,e_0)
      n_t = et_nodes(1)
      t_0 = vols(2*v-1) ! Old Element
      t_t = vols(2*v)   ! New element
      ! In t_t, t_0 will be the neighbour at the local n_t - face
      ! In t_0, t_t will be the neighbour at the local n_0 - face
        do n = 1, 4 ! n = n_0 in t_0
          if (vn(n,t_0).eq.n_0) exit
        end do
        do m = 1, 4 ! m = n_t in t_t
          if (vn(m,t_t).eq.n_t) exit
        end do

        vv(n,t_0) = t_t
        vv(m,t_t) = t_0
      end do

      ! Internal sidewise neighbours resolution:
      do e = 1, 2
        ed = edges_intern(e)
        do v = 1, size(ev(ed)%d)
          vol = ev(ed)%d(v)
          ! Loop over all edges of vol to find local position e of edge ed
          do e1 = 1,6
            if (ve(e1,vol).eq.ed) exit
          end do
          do k = 1,size(ev(ed)%d)
            vol_2 = ev(ed)%d(k)
            if (vol .eq. vol_2) cycle
            do e2 = 1,6
              if (ve(e2,vol_2).eq.ed) exit
            end do
            ! If vol and vol_2 share the same edge and same opposite global node, they
            ! are neighbours.
            if ( vn(opsn(1,e1),vol) .eq. vn(opsn(1,e2),vol_2) ) then
              vv(fen(1,e1),vol) = vol_2
              vv(fen(1,e2),vol_2) = vol
            else if ( vn(opsn(1,e1),vol) .eq. vn(opsn(2,e2),vol_2) ) then
              vv(fen(1,e1),vol) = vol_2
              vv(fen(2,e2),vol_2) = vol
            else if ( vn(opsn(2,e1),vol) .eq. vn(opsn(1,e2),vol_2) ) then
              vv(fen(2,e1),vol) = vol_2
              vv(fen(1,e2),vol_2) = vol
            else if ( vn(opsn(2,e1),vol) .eq. vn(opsn(2,e2),vol_2) ) then
              vv(fen(2,e1),vol) = vol_2
              vv(fen(2,e2),vol_2) = vol
            end if
          end do
        end do
      end do
      ! In case of surfedge:
      if (present(surfs)) then
        do v = 1, size(vols)
          vol = vols(v)
          do s = 1, size(surfs)
            do f = 1,4 ! If all three points of a face of vn coincide with the three points of the surface, the surface will be assigned to vv:
              if ( (vn(fsn(1,f),vol).eq.sn(1,surfs(s))).and.    &
     &             (vn(fsn(2,f),vol).eq.sn(2,surfs(s))).and.    &
     &             (vn(fsn(3,f),vol).eq.sn(3,surfs(s))) ) vv(f,vol) = - surfs(s)
            end do
          end do
        end do
      end if
!_End.
!  
      deallocate(edges_involved)
      return
end subroutine updateaux


function EdgeLength( n_0 , n_t )
      use globalvariables3D,   only: nod
      use femtypes
      implicit none
      integer (I4B)               :: n_0, n_t
      real     (DP)               :: EdgeLength
      intent   (in)               :: n_0, n_t
!
!  EdgeLength calculates the length of an edge located between global node numbers n_0 and n_t
!
!
!    input
!        n_0,n_t    End nodes of an edge
!
!    output
!        EdgeLength length of the edge calculated as the distance between it's two endpoints
!
!      EdgeLength = sqrt( ( nod(1,n_0)-nod(1,n_t) )**2 +   &
!     &                   ( nod(2,n_0)-nod(2,n_t) )**2 +   &
!     &                   ( nod(3,n_0)-nod(3,n_t) )**2 )
      EdgeLength = norm2( nod(1:3,n_0)-nod(1:3,n_t) )
      return
end function EdgeLength


subroutine SortEdges( edges , n_0 )
      use feminterface,        only: qsortindex, palloc
      use globalvariables3D,   only: en, nod
      use femtypes
      implicit none
      integer (I4B), pointer, optional :: n_0(:)
      integer (I4B),           pointer :: edges(:)
      intent  (inout)                  :: edges, n_0
!
!  LongestEdge determines for the current mesh, which of the given edges is the longest one
!
!
!    input
!        edges      List of global edge numbers
!
!    output
!        SortEdges  List of global edge numbers sorted ascending according to their lengths
!
!  local variables:
!                  
      integer (I4B)              :: edges_size
      integer (I4B), allocatable :: e_index(:)
      integer (I4B),    pointer  :: edges_temp(:)=>null(), n_0_temp(:)=>null()
      real    (DP) , allocatable :: e_len(:)

      edges_size = size(edges)
      allocate(e_len(edges_size))
      allocate(e_index(edges_size))


      e_len(1:edges_size) = - sqrt( ( nod(1,en(1,edges(:)))-nod(1,en(2,edges(:))) )**2 + &
     &                              ( nod(2,en(1,edges(:)))-nod(2,en(2,edges(:))) )**2 + &
     &                              ( nod(3,en(1,edges(:)))-nod(3,en(2,edges(:))) )**2 )

      call qsortindex(e_len,e_index,edges_size)
!__
!
! 1) edges:
      call palloc(edges_temp,edges_size)
      edges_temp = edges(e_index(:))
      deallocate(edges)
      edges => edges_temp
      nullify(edges_temp)
!__
!
! 2) n_0:
      if (present(n_0)) then
        call palloc(n_0_temp,edges_size)
        n_0_temp   = n_0(e_index(:))
        deallocate(n_0)
        n_0 => n_0_temp
        nullify(n_0_temp)
      end if
!__
!
! Release Memory:
      deallocate(e_index)
      deallocate(e_len)
!__
      return
end subroutine SortEdges