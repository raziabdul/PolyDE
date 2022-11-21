      subroutine getbcval2D( branch, nature, xs, ys, pval, qval, elem )
      use feminterface, only: userbc
      use femtypes
      use globalvariables, only: zrb, alrb, btrb
      implicit none
      integer (I4B) :: branch, nature
      integer (I4B) , optional :: elem
      real (DP) :: xs, ys
      complex(DPC) :: pval
      complex(DPC), optional :: qval(:)
      intent (in)   :: branch, nature, xs, ys, elem
      intent (out)  :: pval, qval
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
!    $Date: 2011/08/09 09:55:53 $
!    $Author: chvokas $
!
!  return parameters of  dirichlet or gerenral boundary condition 
!  for given the BC-number at coordinate (xs,ys)
!
!  Input:
!        branch      branch number from geometry
!        xs, ys      coordinates for determining dirichlet
!  Output:
!     in the case of dirichlet BC
!        pval        value for the scalar potential   u=pval
!     in the case of general BC
!        pval
!        qval        (nu*grad u + f2)*n = pval + qval*u
!
!  local variables
      integer (I4B) bcnum
      bcnum=zrb(branch,nature)
!  If there is a Dirichlet- or General-BC then assign values directly. 
!  Call userbc if bcnum is 1 to 99 or 201-299. Print error in other cases.

!  zero bcnum refers homogenous Dirichlet BC, 200 refers to homogenous Neumann
      if ( bcnum .eq. 0) then
        pval = alrb(branch,nature)
      else if ( bcnum .eq. 200) then
        pval = alrb(branch,nature)
        qval = 0._DPC
        qval(nature) = btrb(branch,nature)
      else if (bcnum .ge.   1 .and. bcnum .lt. 100) then
        call userbc(branch,bcnum,xs,ys,nature,pval,qval)
      else if (bcnum .ge. 201 .and. bcnum .lt. 300) then
        call userbc(branch,bcnum,xs,ys,nature,pval,qval,elem)
      else
        print *, "**** No boundary condition defined on this branch"
      end if
!
      return
      end subroutine getbcval2D
