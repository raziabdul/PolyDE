subroutine preassemb_mixed(nudof,npdof)
use femtypes
use feminterface, only: reallocate
use fluidvariables
use fluidinterface
use globalvariables
implicit none
integer(I4B), intent(out) :: nudof,npdof
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
!    $Revision: 1.6 $
!    $Date: 2014/07/15 13:02:35 $
!    $Author: m_kasper $
!
!*********************************************************************************
!  Generation of the degrees of freedom (global numbers) for all elements
!  depending on the polynomial degree
!*********************************************************************************
!  Input:
!       e       element information (nodes of the element)
!       en      neighborhood information, neighbors of the elements
!       nel     total number of elements
!       nnode   total number of nodes
!       ep      polynomial degree for each of the elements
!
!  Output:
!       eg      degrees of freedom (global) for the elements
!               packed into an arry of pointers   eg(i)%d(j)
!       ndof    total number of degrees of freedom
!*********************************************************************************
!      neldof = (ep(i)+1)*(ep(i)+2)/2                             ! normal element
!      neldof = ((ep(i)+1)*(ep(i)+2))+(ep(i)*(ep(i)+1))/2         ! fluid element
!*********************************************************************************
!  local variables
!*********************************************************************************
integer(I4B) :: i, j, k, ln, from, upto
integer(I4B) :: nextdof, neighbor, edgedegree, jnachb
integer(I4B) :: pressdof, ip, degree
integer(I4B) :: udof, udof2, pdof, neldof
integer, parameter, dimension(3) :: inach = (/2,3,1/)
!*********************************************************************************
nextdof = p+1
allocate(eg(n,1))

! use when adaption is activated.
degree = 0 

do i = 1,n
   udof = (ep(i,1)+1)*(ep(i,1)+2)/2
   udof2 = udof*2
   pdof = ep(i,1)*(ep(i,1)+1)/2
   neldof = 2*udof + pdof

   allocate(eg(i,1)%d(neldof)) ! neldof is a global variable
!  assign the vertex dof
   eg(i,1)%d(1) = e(1,i)
   eg(i,1)%d(2) = e(2,i)
   eg(i,1)%d(3) = e(3,i)
   
   degree = max(ep(i,1),degree)
end do

if (degree .GE. 2) pressdof = nextdof

! the node numbering is set for all elements order by order (order = polynomial order).
do ip = 2,degree  ! polydeg loop
   do i = 1,n

      if (ip .gt. ep(i,1)) then
         do j = 1,3
            neighbor = en(j,i)
            
            do k = ep(i,1)+1,ip 
              ln=k*(k+1)/2+inach(j)
              eg(i,1)%d(ln)=0
            end do
            if (neighbor .gt. 0) then
              do k = ep(neighbor,1)+1, ip
                ln=k*(k+1)/2+inach(jnachb)
                eg(neighbor,1)%d(ln)=0
              end do
            end if
         end do
         cycle
      end if

      do j = 1,3  ! edge
         neighbor = en(j,i)
         if ( i .gt. neighbor) then ! if i > neighbor
!  nodes on the edges had not been fixed yet
!  the degree of the edge depends on the neighbor, choose the smaller
            if (neighbor .gt. 0) then
               edgedegree = MIN(ep(i,1),ep(neighbor,1))
            ! search the edge in the neighboring element
               do k = 1,3
                  if (en(k,neighbor) .eq. i) then
                     jnachb = k
                     exit
                  end if
               end do
            else
            ! there is no neighbor
               edgedegree = ep(i,1)
            end if

            k = ip
            ln = k*(k+1)/2+inach(j)
            eg(i,1)%d(ln) = nextdof

            if (neighbor .gt. 0) then
               ln = k*(k+1)/2+inach(jnachb)
         !  assign the same dof in the neighboring element
         !  however, with negative sign if degree is odd
               eg(neighbor,1)%d(ln) = nextdof*(-1)**k

            end if

            nextdof = nextdof+1

            if (ip .LE. (degree-1)) then
               pressdof = pressdof+1
            end if

         end if  ! if i > neighbor
      end do ! edge

      if (ip .ge. 3) then
!  assign the face (inner) dof
!        do k = 3,ep(i)
         k = ip
         from = k*(k+1)/2+4
         upto = (k+1)*(k+2)/2
         do ln = from, upto
            eg(i,1)%d(ln) = nextdof
            nextdof = nextdof+1

            if (k .LE. (ep(i,1)-1)) then
               pressdof = pressdof+1
            end if
         end do
!        end do
      end if
   end do ! element loop
end do ! polydeg loop

ndof = nextdof-1
pressdof = pressdof-1

do i = 1,n
   do j = udof+1,udof2
      if (eg(i,1)%d(j-udof) .GT. 0) then
         eg(i,1)%d(j) = eg(i,1)%d(j-udof)+ndof
      else if (eg(i,1)%d(j-udof) .LT. 0) then
         eg(i,1)%d(j) = eg(i,1)%d(j-udof)-ndof
      end if
   end do

   do k = udof2+1,neldof
      if (eg(i,1)%d(k-udof2) .GT. 0) then
         eg(i,1)%d(k) = eg(i,1)%d(k-udof2)+(2*ndof)
      else if (eg(i,1)%d(k-udof2) .LT. 0) then
         eg(i,1)%d(k) = eg(i,1)%d(k-udof2)-(2*ndof)
      end if
   end do
end do
nudof = ndof            ! number of total dof for only U
npdof = pressdof        ! number of total dof for only P
return
end subroutine preassemb_mixed