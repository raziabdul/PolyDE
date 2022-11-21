      subroutine glaett(n,e,xn,yn,en,geb,kzi,p,nregen,liste5,           &
     &  fehlkn,zrb,zki,xbk,ybk,loops)
      use feminterface, only: lsmoth, netbew, lswap, inkrp
      use femtypes
      implicit none
      integer (I4B) :: n, e(:,:), en(:,:), geb(:), kzi(:), p
      integer (I4B) :: liste5(:), fehlkn, zrb(:,:), zki(:,:), loops
      real (DP) :: xn(:), yn(:), xbk(:), ybk(:)
      logical :: nregen
      intent (in) :: n, geb, kzi, p, zrb, zki, xbk, ybk, loops
      intent (inout) :: e, xn, yn, en, nregen, liste5, fehlkn
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
!    $Date: 2008/12/22 12:43:49 $
!    $Author: m_kasper $
!
!  Mesh smoothing by use of swaping triangle edges and node relocation (Laplace smoothing)
!
!  Input:
!            n        total number of elements
!            geb      assignment of elements to the geometry regions
!            kzi      node branch information
!                     kzi=0: the node is an inner node
!                     kzi<0: the node is identical to the keypoint with number (-kzi)
!                     kzi>0: the node is on the branch with the number (kzi)
!            p        total number of nodes
!            zrb      Type of the boundary condition of the branch, 
!                     which can take the following values:
!                       0 -  99: Dirichlet BC, where for 1- 99: the BC has to be set in the file USERBC
!                     200 - 299: Neumann of general BC
!                     300 - 399: visible contour, no BC
!                     400 - 499: non-visible contour, no BC
!            zki      list of of key-points (geometry nodes) which determine the branches
!                     zki(1,i)  start key-point of the branch
!                     zki(2,i)  end key-point of the branch
!                     zki(3,i)  third key-koint of the branch, with
!                               /  = 0  :  straight line 
!                     zki(3,i)  -  > 0  :  node located on an arc
!                               \  < 0  :  midpoint of the arc
!            xbk, ybk coordinates of the input or geometry node (key-points)
!            loops    max. number of swap & smooth iterations
!
!  In-/ Output: 
!            e        element information, nodes of the elements
!            xn, yn   coordinates of the nodes
!            en       neighborhood information, neighbors of the elements
!            nregen   =.true. if the mesh has to be regenerated,
!                     i.e. if an element has negative area
!            liste5   =1 , if the node has to be relocated
!            fehlkn   number of nodes which have to be relocated
!
!  local variables
      integer (I4B) :: szahl, zaehl, kz1, i, maxit
      real (DP) :: sum, rho1, rho, rhocp, xi, yi, flae2, ri, xctr, yctr
      logical :: merk, smooth
      parameter (maxit=10)
!
!  check if branche nodes on arcs are at the correct radius
!  otherwise mark them for relocation
      do i=1,p
        if (kzi(i) .gt. 0) then
          if (liste5(i) .eq. 0) then
            kz1=kzi(i)
!  only on arcs
            if (zki(3,kz1) .ne. 0) then
              if (zki(3,kz1) .lt. 0) then
!  arc given by startpoint, endpoint and center
                xctr=xbk(-zki(3,kz1))
                yctr=ybk(-zki(3,kz1))
                rho=sqrt((xbk(zki(1,kz1))-xctr)**2+                     &
     &                   (ybk(zki(1,kz1))-yctr)**2)
                rhocp=sqrt((xbk(zki(2,kz1))-xctr)**2+                   &
     &                     (ybk(zki(2,kz1))-yctr)**2)
                rho=(rho+rhocp)/2._DP
              else if (zki(3,kz1) .gt. 0) then
!  arc given by startpoint, endpoint and third point
                call inkrp(xbk(zki(1,kz1)),ybk(zki(1,kz1)),             &
     &            xbk(zki(2,kz1)),ybk(zki(2,kz1)),                      &
     &            xbk(zki(3,kz1)),ybk(zki(3,kz1)),                      &
     &            xi,yi,flae2,ri,xctr,yctr,rho)
              end if
!  rho  = desired radius
!  rho1 = actual radius
              rho1=sqrt((xn(i)-xctr)**2+(yn(i)-yctr)**2)
              if(abs(rho1-rho)/(rho1+rho) .gt. 1.e-5_DP) then
!  the difference between actual and desired radius is too large
!  mark this node for relocation
                liste5(i)=1
                fehlkn=fehlkn+1
              end if
            end if
          end if
        end if
      end do
!             
      zaehl=1
!  start with the mesh smooting
      if (nregen) then
        call lsmoth(kzi,p,xn,yn,e,n,en,4,nregen,liste5,fehlkn,          &
     &    zrb,zki,xbk,ybk)
        call lswap(n,e,xn,yn,en,geb,1,.false.,szahl,2)
        merk=.true.
      else
        call lswap(n,e,xn,yn,en,geb,1,.false.,szahl,2)
        merk=.false.
      end if
!
      smooth=.true. 
!  do smoothing and swaping alternately until:
!     there is no more node which requires relocation
!  or a maximum number of iteration (maxit) has been performed
      do while (smooth)
        call lsmoth(kzi,p,xn,yn,e,n,en,4,nregen,liste5,fehlkn,          &
     &    zrb,zki,xbk,ybk)
        zaehl=zaehl+1
        if ((zaehl.lt.loops).or.(nregen.and.(zaehl.lt.maxit))) then
          smooth=.true.
        else if ((.not.nregen).and.merk) then
!  guarantee that after relocation of nodes on further pass of smoothing is done
          merk=.false.
          zaehl=1
          smooth=.true.
        else 
          smooth=.false. 
        end if
        if (smooth) call lswap(n,e,xn,yn,en,geb,1,.false.,szahl,2)
        if (szahl .le. 1) smooth=.false.
      end do
!
      if (nregen) then
        write (*,4322) fehlkn
4322    format('**** glaett: was not able to fully relocate the node:',i6)
      end if
      call netbew(xn,yn,e,sum,n)
      return
      end subroutine glaett
