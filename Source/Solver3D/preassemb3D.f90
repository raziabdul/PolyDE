subroutine preassemb3D()
      use feminterface3D, only: getepbc, getvf, bcs2vert, nedelecdof, scalardof
      use femtypes
      use globalvariables3D, only : eltype, numdof, numf, vf
      implicit none
!
!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 - 2015 Institute for Micro Systems Technology,
!                              Hamburg University of Technology.
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
!    $Revision: 1.35 $
!    $Date: 2015/11/10 15:38:22 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'preassemb3D'
!  Generates degrees of freedom (global numbers) for all elements
!  depending on the polynomial degree
!
!-------------------------------------------------------------------------------
!  
!  Input:
!            eltype   type of element used: NEDELEC, SCALAR 
!                     (future: multi-natured, etc.) (global)
!            vn       element information (nodes of the element) (global)
!            vv       neighbour information      allocate(vv(4,numv))
!                     vv(i,j) is the element adjacent to j at element face i
!                     (opposite to the vertex i). If element j has no neighbour at
!                     face i,  vv(i,j)=0. (global)
!            numn     number of nodes (global)
!            numv     number of volume elements (global)
!            vp       polynomial degree of volume element    allocate(vp(numv)) (global)
!
!  Output:
!            vgdof    degrees of freedom (global) for the elements
!                     packed into an array of pointers   vgdof(i)%d(j) (global)
!            numdof   total number of degrees of freedom (global)
!
!-------------------------------------------------------------------------------
!
!  Local variables:
      logical :: ok 
!__
! Definitions:
      !print "(A2)"         ,'| '
      !print "(A79)"        ,' _______________________________________________________________________________ '
      !print "(A49)"        ,'|                           STARTING  PREASSEMBLY'
      !print "(A2)"         ,'| '
!__
! Checks and preparations:
   
!__
!
! Allocate Memory:
   
!__
! Start:
! 1) Get global edges numbers, en & fn is never used so far, so discard immediately:
      call getvf(numf, vf)
!-
      if (eltype .eq. 'NEDELEC') then
        print "(A43)"   ,'| Using Nedelec-Type Element (Vector Field)'
        call nedelecdof(numf, vf)
      elseif (eltype .eq. 'SCALAR') then
        print *, 'Using scalar-type element'
        call scalardof(numf, vf)
        call bcs2vert(ok)
      else
        print *, ' Element type not supported. Program Terminated.'
        return
      end if
!__
! 2) Assign edge info for BC specifications:
      call getepbc
              print "(A3)" ,'|--'
      print "(A29,I)"      ,'| DEGREES OF FREEDOM (NDOF): ', numdof
!
! End.
!___
      return
end subroutine preassemb3D


subroutine nedelecdof(numf, vf)
      use feminterface,            only : destroyarrptr
      use globalvariables3D,       only : ep, fp, nnat, numdof, nume, numv, polymaxv, vgdof, ve, vp
      use femtypes
      implicit none
      integer (I4B), pointer           :: vf(:,:)
      integer (I4B), intent (in)       :: numf
!
!-------------------------------------------------------------------------------
!
! This Routine: 'nedelecdof'
!  Generates degrees of freedom (global numbers) for all elements
!  depending on the polynomial degree according to Nedelec elements
!  The number of the DOF for the volume element follows Nedelec formula  
!  for a polynomial order: 
!                          xxx-xxx-xxx
!
!-------------------------------------------------------------------------------
!  
!  Input:
!            numv     total number of volume elements (global)
!            vn       element information (nodes of the element) (global)
!            vp       polynomial degree for each of the volume elements (global)
!
! USED:      oldvp    if adaption the former vp vector
!                     if oldvp is not present a new eg vector is generated
!                     otherwise eg is updated according to the vp vector
!  Output:
!            vgdof    degrees of freedom (global) for the elements
!                     packed into an array of pointers   eg(i)%d(j)
!            numdof     total number of degrees of freedom
!
!-------------------------------------------------------------------------------
!
!  Local variables:
      type (ARRPTRI), pointer :: eg(:,:)=>null(), fg(:,:)=>null()
      integer (I4B) :: p, i, j, k, maxp, startd, stopd,inat
      integer (I4B) :: dof, edge, face
      integer (I4B), allocatable :: startedge(:), startface(:), startvol(:)
      integer (I4B), allocatable :: numfacedof(:), numvoldof(:) 
      integer (I4B) :: sumface(0:polymaxv) 
      logical :: isnull, ok
!__
! Checks and preparations:
   
!__
!
! Allocate Memory:
      allocate(numfacedof(polymaxv))
      allocate(numvoldof(polymaxv))
      allocate(startedge(polymaxv))
      allocate(startface(polymaxv))
      allocate(startvol(polymaxv))
!__
! Definitions:
! a) Calculate the number of DOF's of a face and vol at element's order p (not including p-1, p-2, ... )
! b) Create starting indices for d in vgdof
! c) Total number of face DOF's of order p
!-  a)
      numfacedof = (/ (2*(p-1), p = 1,polymaxv) /)
      numvoldof  = (/ ( 3*(p-1)*(p-2)/2, p = 1,polymaxv ) /)
!-  b)
      startedge(1) = 0
      startface(1) = -1
      startvol(1)  = 6
      do p = 2,polymaxv
        startedge(p) = startvol(p-1) + numvoldof(p-1)
        startface(p) = startedge(p)  + 6
        startvol(p)  = startface(p)  + numfacedof(p)*4
      end do
!-  c)
      sumface(0) = 0
      do p = 1,polymaxv
        sumface(p) = sum(numfacedof(1:p))
      end do
!__
! Start:
! 1) Compute 'ep' : polynomial order of all edges
      allocate(ep(nume,nnat))
      ep = polymaxv
      do inat=1,nnat
        do i = 1,numv
          do j = 1,6
            edge = abs(ve(j,i))
            ep(edge,inat)  = min(ep(edge,inat),vp(i,inat))
          end do
        end do
      end do
!__
! 2) Compute 'fp' : polynomial order of all faces
      allocate(fp(numf,nnat))
      fp = polymaxv
      do inat=1,nnat
        do i = 1,numv
          do j = 1,4
            face = vf(j,i)
            fp(face,inat) = min( fp(face,inat), vp(i,inat) )
          end do
        end do
      end do
!__
! 3) Compute edges dofs:
      allocate(eg(nume,nnat))
      dof = 0
      do inat=1,nnat
        do i = 1,nume
          allocate(eg(i,inat)%d(ep(i,inat)))
          eg(i,inat)%d = (/ ( j+dof, j = 1,ep(i,inat) ) /)
          dof = dof+ep(i,inat)
        end do
      end do
!__
! 4) Compute faces dofs:
      allocate(fg(numf,inat))
      do inat=1,nnat
        do i = 1,numf
          allocate(fg(i,inat)%d(sumface(fp(i,inat))))
          fg(i,inat)%d = (/ ( ( j+sumface(k-1)+dof, j=1,numfacedof(k) ) , k=2,fp(i,inat) ) /)
          dof = dof+sumface(fp(i,inat))
        end do
      end do
!__
! 5) Assign dof to elements, i.e, copy the dof from above
!    vol: loop over the polynomial orders
!  order: loop over the edges
!      j: loop over the faces
!-
      if (associated(vgdof)) then
        ok = destroyarrptr(vgdof)
      end if
      nullify(vgdof)
      allocate(vgdof(numv,nnat))
!-
      do inat=1,nnat
vol:    do i = 1,numv
          maxp = min(vp(i,inat), polymaxv)
          allocate( vgdof(i,inat)%d( maxp*(maxp+2)*(maxp+3)/2) )
order:    do p = 1,maxp
!-
            do j = 1,6
              edge = abs(ve(j,i))
              if (ep(edge,inat) .ge. p) then
                vgdof(i,inat)%d(startedge(p)+j) = eg(edge,inat)%d(p)
              else
                vgdof(i,inat)%d(startedge(p)+j) = 0
              end if
            end do
!-
            do j = 1,4
              face = vf(j,i)
              startd = startface(p) + numfacedof(p)*(j - 1) + 1
              stopd = startd + numfacedof(p) - 1
              if (fp(face,inat) .ge. p) then
                vgdof(i,inat)%d( startd:stopd ) = fg(face,inat)%d(sumface(p-1)+1:sumface(p-1)+numfacedof(p))
              else
                vgdof(i,inat)%d( startd:stopd ) = 0
              end if
            end do
!- inner dofs:
            vgdof(i,inat)%d(startvol(p)+1:startvol(p)+numvoldof(p)) = (/ (dof+j , j = 1 , numvoldof(p)) /)
            dof = dof + numvoldof(p)
          end do order
        end do vol
      end do
!__
! 6) Compute the total number of dofs:
      numdof = dof
!__
!
! Release Memory:
      isnull = destroyarrptr(eg)
      isnull = destroyarrptr(fg)
      deallocate(numfacedof)
      deallocate(numvoldof)
      deallocate(startedge)
      deallocate(startface)
      deallocate(startvol)
      deallocate(fp)
      nullify(fp)
!
! End.
!___
 return
end subroutine nedelecdof


subroutine scalardof(numf, vf)
      use feminterface,            only : destroyarrptr
      use globalvariables3D,       only : ep, fp, nnat, numdof, nume, numn, numv, polymaxsc, vgdof, ve, vn, vp
      use femtypes
      implicit none
      integer (I4B)                    :: numf
      integer (I4B),           pointer :: vf(:,:)
      intent (in) numf, vf
!
!-------------------------------------------------------------------------------
!
! This Routine: 'scalardof'
!  Generates degrees of freedom (global numbers) for all elements
!  depending on the polynomial degree according to scalar elements.
!  The number of the DOF for the scalar tetrahedron element, dependence of polynomial order p, is: 
!
!   on vertices: m0=1
!      on edges: m1=p-1
!      on faces: m2=(p-1)(p-2)/2
!  at the inner: m2=(p-1)(p-2)(p-3)/6
!  
!-------------------------------------------------------------------------------
!  
!  Input:
!            numv       total number of volume elements (global)
!            vn         element information (nodes of the element) (global)
!            vp         polynomial degree for each of the volume elements (global)
!
! USED:      oldvp      if adaption the former vp vector
!                       if oldvp is not present a new eg vector is generated
!                       otherwise eg is updated according to the vp vector
!  Output:
!            vgdof      degrees of freedom (global) for the elements
!                       packed into an array of pointers   eg(i)%d(j)
!            numdof     total number of degrees of freedom
!
!-------------------------------------------------------------------------------
!
!  Local variables:
      type (ARRALLOI), allocatable :: fg(:,:), eg(:,:)
      integer (I4B) :: inat, elementorder
      integer (I4B) :: p, i, j, k, maxp, startd, stopd
      integer (I4B) :: dof, edge, face, ndof0,  ndof1,  ndof2,  ndof3
      integer (I4B), allocatable :: startedge(:), startface(:), startvol(:)
      integer (I4B), allocatable :: numfacedof(:), numvoldof(:) 
      integer (I4B) :: sumface(0:polymaxsc) 
      logical :: ok
!__
! Checks and preparations:

!__
!
! Allocate Memory:
      allocate(numfacedof(polymaxsc))
      allocate(numvoldof(polymaxsc))
      allocate(startedge(polymaxsc))
      allocate(startface(polymaxsc))
      allocate(startvol(polymaxsc))
!__
! Definitions:
!
!   I N F O R M A T I O N:
! sorting of DOFs (shape functions) is:
! number    1..4 | 5..10|11..16|17..20|21..26|27..34| 35 |36..41|42..53|54..56|...
! type     vertex| edge | edge | face | edge | face |vol | edge | face |  vol |...
! order       1  |   2  |   3  |   3  |   4  |   4  |  4 |   5  |   5  |   5  |...
!   TO DO: change element DOF sequence (for static condensation) to be:
! number    1..4 | 5..10|11..16|17..20|21..26|27..34|35..40|41..52|... xx |+1..+4|...
! type     vertex| edge | edge | face | edge | face | edge | face |...vol |  vol |...
! order       1  |   2  |   3  |   3  |   4  |   4  |   5  |   5  |...  4 |   5  |...
!
!    Calculations:
! a) Calculate the number of DOF's of a face and vol at element's order p (not including p-1, p-2, ... )
! b) Create starting indices for local dofs' indexing, these are the indices +1 at which edge, face or vol dofs appear at degree p
! c) Total number of face DOF's of order p
!-  a)
      numfacedof(1) = 0
      numvoldof(1)  = 0
      do p = 2,polymaxsc
        numfacedof(p) = p-2
        numvoldof(p)  = (p-2)*(p-3)/2
      end do
!-  b)
      startedge(1) = 0
      startface(1) = 0
      startvol(1)  = 4
      do p = 2,polymaxsc
        startedge(p) = startvol(p-1) + numvoldof(p-1) 
        startface(p) = startedge(p)  + 6
        startvol(p)  = startface(p)  + numfacedof(p)*4
      end do
!-  c)
      sumface(0) = 0
      do p = 1,polymaxsc
        sumface(p) = sum(numfacedof(1:p))
      end do
!__
! Start:
! 1) Compute 'ep' : polynomial order of all edges
!    Use minimum rule to assign polynomial order to all edges and faces
!-
!$omp parallel sections private(inat, i, j, elementorder)
!$omp section
      allocate( ep(nume,nnat) )
      ep = polymaxsc
      do inat=1,nnat
        do i = 1,numv
          elementorder = vp(i,inat)
          do j = 1,6
            edge = abs(ve(j,i))
            ep(edge,inat)  = min( ep(edge,inat),elementorder )
          end do
        end do
      end do
!__
! 2) Compute 'fp' : polynomial order of all faces
!-
!$omp section
      allocate( fp(numf,nnat) )
      fp = polymaxsc
      do inat=1,nnat
        do i = 1,numv
          elementorder = vp(i,inat)
          do j = 1,4
            face = vf(j,i)
            fp(face,inat) = min( fp(face,inat),elementorder )
          end do
        end do
      end do
!$omp end parallel sections
!__
! 3) Compute edges dofs:
!    'dof' counts the number of DOF assigned so far
!  a) The vertex DOF will be copied from the labeling of mesh nodes
!  b) The DOF will be labeled in the sequence: vertexdofs , edgedofs, facedofs, voldofs
!
!-  a)
      dof = numn*nnat
      ndof0 = dof
      print*,'| Number of vertex DOFs:',ndof0
!
!-  b) Num. of edge dof = p-1
      allocate(eg(nume,nnat))
      do inat=1,nnat
        do i = 1,nume
          allocate(eg(i,inat)%d(ep(i,inat)-1))
          eg(i,inat)%d = (/ ( j+dof, j = 1,ep(i,inat)-1 ) /)
          dof = dof+ep(i,inat)-1
        end do
      end do
      ndof1 = dof-ndof0
      print*,'| Number of edge DOFs  :',ndof1


      if (associated(vgdof)) then
        ok = destroyarrptr(vgdof)
      end if
      allocate( vgdof(numv,nnat))

! 4) Assign dof to elements, i.e, copy the dof from above
!    a) loop over vertices
!    b) loop over the edges
!$OMP PARALLEL DO default (none)                                        &
!$OMP shared(nnat, numv, vgdof, vp)                                     &
!$OMP shared(numn, ve, vn, ep, startedge, eg)                           &
!$OMP private(inat, i, maxp, edge, j, p)
      do inat=1,nnat
        do i = 1,numv
          maxp = min(vp(i,inat), polymaxsc)
          allocate ( vgdof(i,inat)%d( (maxp+1)*(maxp+2)*(maxp+3)/6) )

!-      a)
          do j = 1,4
            vgdof(i,inat)%d(j) =  vn(j,i)+((inat-1)*numn)
          end do
!-      b)
          do p = 2,maxp
            do j = 1,6
              edge = abs(ve(j,i))
              if ( (ep(edge,inat) .ge. p) .and. (ep(edge,inat) .gt. 1) ) then
                vgdof(i,inat)%d(startedge(p)+j) = eg(edge,inat)%d(p-1)
              else
                vgdof(i,inat)%d(startedge(p)+j) = 0
              end if
            end do
          end do
!
        end do
      end do
!$OMP END PARALLEL DO

      deallocate(eg)

!__
! 5) Compute faces dofs:
      allocate(fg(numf,nnat))
      do inat=1,nnat
        do i = 1,numf
          allocate(fg(i,inat)%d(sumface(fp(i,inat))))
          fg(i,inat)%d = (/ ( ( j+sumface(k-1)+dof, j=1,numfacedof(k) ) , k=3,fp(i,inat) ) /)
          dof = dof+sumface(fp(i,inat))
        end do
      end do
      ndof2 = dof-ndof0-ndof1
      print*,'| Number of face DOFs  :',ndof2
!__
! 6) Assign dof to elements, i.e, copy the dof from above
!    c) loop over the faces
!-
!$OMP PARALLEL DO default (none)                                        &
!$OMP shared(nnat, numv, vgdof, vp)                                     &
!$OMP shared(fp, startface, vf, numfacedof, fg, sumface)                &
!$OMP private(inat, i, maxp, j, p, face, startd, stopd)
      do inat=1,nnat
        do i = 1,numv
          maxp = min(vp(i,inat), polymaxsc)
!
!-      c)
          do p = 3,maxp
            do j = 1,4
              face = vf(j,i)
              startd = startface(p) + numfacedof(p)*(j - 1) + 1
              stopd = startd + numfacedof(p) - 1
              if ((fp(face,inat) .ge. p) .and. (fp(face,inat) .gt. 2) ) then
                vgdof(i,inat)%d(startd:stopd) = fg(face,inat)%d(sumface(p-1)+1:sumface(p-1)+numfacedof(p))
              else
                vgdof(i,inat)%d(startd:stopd) = 0
              end if
            end do
          end do
!- c)
        end do
      end do
!$OMP END PARALLEL DO
!__
! 7) Compute inner dofs:
      do inat=1,nnat
        do i = 1,numv
          maxp = min(vp(i,inat), polymaxsc)
          do p = 4,maxp
            vgdof(i,inat)%d(startvol(p)+1:startvol(p)+numvoldof(p)) = (/ (dof+j, j = 1, numvoldof(p)) /)
            dof = dof + numvoldof(p)
          end do
        end do
      end do
!__
! 8) Compute total number of dofs:
      numdof = dof
      ndof3 = numdof-ndof0-ndof1-ndof2
      print*,'| Number of inner DOFs :',ndof3
!__
!
! Release Memory:
      deallocate(fg)
      deallocate(numfacedof)
      deallocate(numvoldof)
      deallocate(startedge)
      deallocate(startface)
      deallocate(startvol)
      deallocate(fp)
      nullify(fp)
!
! End.
!___
 return
end subroutine scalardof

