      subroutine earclipper(head,numnodes,region)
      use feminterface, only: reallocate, wink, ninnen
      use globalvariables
      use femtypes
      use triangtypes
      implicit none
      integer (I4B) region, numnodes
      type (node), pointer :: head
      intent (in) :: region
      intent (inout) :: numnodes
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
!    $Date: 2008/08/06 12:41:14 $
!    $Author: m_kasper $
!
!  triangulate a simply connected region by ear clipping
!
!  Input:
!            region   index of the actual region
!            numnodes number of nodes of the boundaries
!            head     pointer to the starting node
!
!  local variables
      integer (I4B) i, j, k, ntria
      real (DP) xmi, ymi, xla, yla, xne, yne
      real (DP) xll, yll, xnn, ynn, earangle, reflexangle, minangle
      logical ear, done
      type (node), pointer :: minode, lanode, nenode, knode, tmpnode

!  preparation: compute the inner angles at the nodes
      minode=>head
      minangle=huge(1._DP)
      do i=1,numnodes
        xmi=xn(minode%nodeindx)
        ymi=yn(minode%nodeindx)
        lanode=>minode%prev
        xla=xn(lanode%nodeindx)
        yla=yn(lanode%nodeindx)
        nenode=>minode%next
        xne=xn(nenode%nodeindx)
        yne=yn(nenode%nodeindx)
        minode%angle=wink(xmi,ymi,xne,yne,xla,yla)
!  for starting earclipping find the node of minimum angle
        if (minode%angle .lt. minangle) then 
          minangle=minode%angle
          head=>minode
        end if
        minode=>minode%next
      end do
      write(*,100) region, minangle*180._DP/pi
100   format('Minimum angle of region ',i6,' is :',f8.3,'  degrees')
!  check whether arrays are allocated or need reallocation
      if (associated(geb)) then
        if (size(geb) .lt. n+numnodes-2) then
          e  =>reallocate(e,3,2*(n+numnodes))
          geb=>reallocate(geb,2*(n+numnodes))
        end if
      else 
        allocate(e(3,gbz*numnodes))
        allocate(geb(gbz*numnodes))
      end if
!  the maximum angle which potentially can be an ear
      earangle=pi-1.e1_DP*spacing(pi)
!  angles > 180 are reflex
      reflexangle=pi
!  a polygon of n vertices is triangulated with n-2 triangles
      ntria=numnodes-2
      done=.false.
      do i=1,ntria
!  three consecutive vertices of a polygon V(i-1), V(i), V(i+1) form an ear if no 
!  polygon vertex is inside the triangle  V(i-1), V(i), V(i+1)
!  a 2D polygon has at least 2 ears
        minode=>head
!  find an ear
        if (numnodes .gt. 3) then
          do j=1,numnodes
!  beeing an ear requires the inner angle to be convex (smaller than Pi-eps)
            if (minode%angle .le. earangle) then
!  check whether a (non-adjacent) vertex is inside the ear candidate
              xmi=xn(minode%nodeindx)
              ymi=yn(minode%nodeindx)
              lanode=>minode%prev
              xla=xn(lanode%nodeindx)
              yla=yn(lanode%nodeindx)
              nenode=>minode%next
              xne=xn(nenode%nodeindx)
              yne=yn(nenode%nodeindx)
              knode=>nenode%next
              ear=.true.
              do k=2,numnodes-2
!  it is sufficient to check reflex vertices (angle > Pi)
!  if it is not an ear there is at least one reflex vertex inside the ear candidate
                if (knode%angle .ge. reflexangle) then
                  if ( ninnen( xmi,ymi, xne,yne, xla,yla,               &
     &                 xn(knode%nodeindx),yn(knode%nodeindx) ) ) then
                    ear=.false.
                    exit
                  end if
                end if
                knode=>knode%next
              end do
              if (ear) exit
            end if
!  test the next ear candidate
            minode=>minode%next
          end do
        else
!  the last triangle, i.e numnodes = 3
          lanode=>minode%prev
          nenode=>minode%next
          ear=.true.
          done=.true.
        end if
!  generate a new element and remove the ear from the polygon
        if (ear) then
          n=n+1
          numnodes=numnodes-1
          geb(n)=region
          e(1,n)=minode%nodeindx
          e(2,n)=nenode%nodeindx
          e(3,n)=lanode%nodeindx
          if (.not. done) then
            lanode%next=>nenode
            nenode%prev=>lanode
!  update angles
            tmpnode=>lanode%prev
            xll=xn(tmpnode%nodeindx)
            yll=yn(tmpnode%nodeindx)
            lanode%angle=wink(xn(e(3,n)),yn(e(3,n)),xn(e(2,n)),yn(e(2,n)),xll,yll)
            tmpnode=>nenode%next
            xnn=xn(tmpnode%nodeindx)
            ynn=yn(tmpnode%nodeindx)
            nenode%angle=wink(xn(e(2,n)),yn(e(2,n)),xnn,ynn,xn(e(3,n)),yn(e(3,n)))
            if (nenode%angle .lt. lanode%angle) then
              head=>nenode
            else
              head=>lanode
            end if
            deallocate(minode)
          else 
            deallocate(minode)
            deallocate(lanode)
            deallocate(nenode)
          end if
        else
          print*,'**** was not able to find an ear'
          stop
        end if
      end do
      return
      end subroutine earclipper
!
!
!
      elemental logical function ninnen(x1,y1,x2,y2,x3,y3,xk,yk)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: x1, y1, x2, y2, x3, y3, xk, yk
      intent (in) :: x1, y1, x2, y2, x3, y3, xk, yk
!
!  test whether the point xk, yk is inside a triangle
!
!  Input:
!            x1, y1   coordinates of the first triangle vertex
!            x2, y2   coordinates of the second triangle vertex
!            x3, y3   coordinates of the third triangle vertex
!            xk, yk   coordinates of the point to test
!  Output:
!            ninnen   =.true. if node k is located inside the traingle of nodes 1, 2, 3 
!
!  local variables
!
      ninnen=.false.
      if ((xk-x3)*(y1-y3) .le. (yk-y3)*(x1-x3)) then
        if ((xk-x1)*(y2-y1) .le. (yk-y1)*(x2-x1)) then
          if ((xk-x2)*(y3-y2) .le. (yk-y2)*(x3-x2)) then
            ninnen=.true.
          end if
        end if
      end if
      return
      end function ninnen
