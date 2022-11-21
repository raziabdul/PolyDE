    subroutine preassemb(oldep)
    use feminterface, only: reallocate
    use femtypes
    use globalvariables, only: n, p, ndof, en, e, nnat, eg, ep
    implicit none
    integer (I4B), optional, pointer:: oldep(:,:)
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
!    $Date: 2008/08/20 12:58:55 $
!    $Author: m_kasper $
!
!  Generation of the degrees of freedom (global numbers) for all elements
!  depending on the polynomial degree
!  
!  Input:
!            e        element information (nodes of the element)
!            en       neighborhood information, neighbors of the elements
!            n        total number of elements
!            p        total number of nodes
!            ep       polynomial degree for each of the elements
!            oldep    if adaption the former ep vector
!                     if oldep is not present a new eg vector is generated
!                     otherwise eg is updated according to the ep vector
!  Output:
!            eg       degrees of freedom (global) for the elements and natures
!                     packed into an array of pointers   eg(i,inat)%d(j)
!            ndof     total number of degrees of freedom
!
!  Local variables
!  inach for nachfolger, the successor
!  ln stands for local node
    integer (I4B)  inach(3), i, j, k, ln, from, upto, nextdof
!  inat stands for index of nature
    integer (I4B) inat
    integer (I4B)  neldof, neighbor, edgedegree, jnachb, olddegree
    logical adaption
!  DATA statements can initialize array element values      
    data inach/ 2, 3, 1/

!  Check whether the parameter oldep is present
!  if present adaption is .true.
      adaption=present(oldep)

!  If adaption is .false. or in other words oldep is not present
!  then .not. adaption is .true. and do the following
      if (.not. adaption) then

!  Check whether we have an old eg vector
        !Check whether the pointer is associated with a target
        if ( associated(eg) ) then
        
          do inat = 1, nnat
            do i = 1, size(eg, DIM=1)
              deallocate(eg(i,inat)%d)
            end do
          end do
          
          deallocate(eg)
          !Disassociate a pointer from a target (nullify)         
          nullify(eg)
          
        end if

!  Initialization: we do not have an eg vector
        allocate(eg(n,nnat))

      end if
        
!  Now by using the variable nnat the following part of the code in this
!  file can be put inside a loop and repeated a number of times defined
!  by the number of physical natures. Care should be taken to modify some
!  variables so that they can accommodate the elements of more than 1
!  physical natures
!  variable eg modified, assume eg to be of type eg(i,inat)
!  for eg first dimension gives the element and the second the nature
!  variable ep modified, assume ep to be of type ep(i,inat)
!  variable oldep modified, assume oldep to be of type oldep(i,inat)
!  variable e not to be modified (geometrical), assume e to be of type e(3,i)
!  for e first dimension gives the node, the second the element
!  variable en not to be modified, assume en to be of type en(3,i)
!  for en first dimension gives the edge, the second the element

!  The following part of the code is concerned with the assignment
!  of the degrees of freedom for the elements
!  Initialise the ndof (number of degrees of freedom) variable to 0
    ndof=0
!  loop over the natures
    do inat = 1, nnat
      
      if (.not. adaption) then
        
!  For all elements assign first order dofs
        do i = 1, n
          
!  Number of Degrees Of Freedom (dof) of this element
!  Example check: if ep(i)=2 then neldof=(2+1)*(2+2)/2=12/2=6
          neldof=(ep(i,inat)+1)*(ep(i,inat)+2)/2
          
!  Allocate array of global dof of the element
!  Every element can have a different number of degrees of freedom
          allocate(eg(i,inat)%d(neldof))
          
!  Assign the vertex dof
          eg(i,inat)%d(1)=e(1,i)+ndof
          eg(i,inat)%d(2)=e(2,i)+ndof
          eg(i,inat)%d(3)=e(3,i)+ndof
          
        end do

        nextdof=ndof+p+1
        
!  Else we have adaption: we have an old ep and eg vector        
      else
        do i = 1, n

          if (ep(i,inat) .gt. oldep(i,inat)) then
!  Number of Degree Of Freedom (dof) of this element
!  The same formula as before
            neldof=(ep(i,inat)+1)*(ep(i,inat)+2)/2
!  The reallocate function can be found in source file reallocate.f90
            eg(i,inat)%d=>reallocate(eg(i,inat)%d,neldof)
            
          else if (ep(i,inat) .lt. oldep(i,inat)) then
            print*,'**** it is not allowed to decrease polynomial order in adaption'
          end if
        end do
        nextdof=ndof+p+1
      end if


!  For all elements assign edge and face dofs
      do i = 1, n
!  Assign the edge dof, visit the edges
        do j = 1, 3
!  The neighbor at the jth edge of the ith element
          neighbor=en(j,i)
!  Only treat this edge if the number of the element to which the edge belongs
!  is larger than the number of the neighboring element at this edge
          if ( i .gt. neighbor) then
!  Nodes on the edges have not been fixed yet
!  If a neighbor exists at this edge of the element then the neighbor
!  has a value greater than 0
            if (neighbor .gt. 0) then
            
!  The degree of the edge depends on the neighbor, choose the smaller
              edgedegree=min(ep(i,inat),ep(neighbor,inat))
              
              if (.not. adaption) then
                olddegree=1
              else
                olddegree=min(oldep(i,inat),oldep(neighbor,inat))
              end if
              
!  Search the 3 edges in the neighboring element and find the edge of
!  the neighbor which is also an edge of the current element i
              do k=1,3
                if (en(k,neighbor) .eq. i) then
                  jnachb=k
                  exit
                end if
              end do

!  Otherwise there is no neighbor              
            else
            
              edgedegree=ep(i,inat)
              
              if (.not. adaption) then
                olddegree=1
              else
                olddegree=oldep(i,inat)
              end if
              
            end if !(neighbor .gt. 0)
            
            do k=olddegree+1,edgedegree
!  Determine the local index of this dof
!  The local index is given by ln
!  Then use ln as the index in the vector eg(i,inat)%d(ln)
!  for the element and assign to it a global dof
!  4, 7, 11, 16, 22, .. for edge 3
!  5, 8, 12, 17, 23, .. for edge 1
!  6, 9, 13, 18, 24, .. for edge 2
!  2  3   4   5   6     for these degrees
              ln=k*(k+1)/2+inach(j)
              eg(i,inat)%d(ln)=nextdof
              if (neighbor .gt. 0) then
                ln=k*(k+1)/2+inach(jnachb)
!  Assign the same dof in the neighboring element
!  However, with negative sign if degree is odd
                eg(neighbor,inat)%d(ln)=nextdof*((-1)**k)
              end if
              nextdof=nextdof+1
            end do
!  All the others are nullified
            do k=edgedegree+1,ep(i,inat)
              ln=k*(k+1)/2+inach(j)
              eg(i,inat)%d(ln)=0
            end do
            
            if (neighbor .gt. 0) then
            
              do k=edgedegree+1,ep(neighbor,inat)
                ln=k*(k+1)/2+inach(jnachb)
                eg(neighbor,inat)%d(ln)=0
              end do
              
            end if
          end if
        end do

!  Assign the face (inner) dof
        if (.not. adaption) then
          olddegree=2
        else
          olddegree=oldep(i,inat)
        end if
        
        do k=olddegree+1,ep(i,inat)
          from=k*(k+1)/2+4
          upto=(k+1)*(k+2)/2
          do ln = from, upto
            eg(i,inat)%d(ln)=nextdof
            nextdof=nextdof+1
          end do
        end do
        
      end do
      ndof=nextdof-1
      
    end do

    return
    end subroutine preassemb
