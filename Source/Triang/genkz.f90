      subroutine genkz(zpz, zpp, branchpt, brnodes, bzi, bzil, bzip, meshsize)
      use feminterface, only: inkrp, zeit, reallocate, getsetting
      use globalvariables
      use femtypes
      implicit none
      integer (I4B) :: zpz(:), bzi(:), bzil(:), bzip(:)
      real (DP) :: meshsize(:)
      integer (I4B), pointer :: branchpt(:), brnodes(:)
      real (DP) :: zpp(:)
      intent (in) :: zpp, bzi, bzil, bzip, meshsize
      intent (inout) :: zpz
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
!    $Date: 2008/07/30 14:36:44 $
!    $Author: m_kasper $
!
!
!  generate the nodes on all branches
!
!  Input:
!            gzz      number of branches
!            gkz      number of key-points
!            zki      list of of key-points (geometry nodes) which determine the branches
!                     zki(1,i)  start key-point of the branch
!                     zki(2,i)  end key-point of the branch
!                     zki(3,i)  third key-koint of the branch, with
!                               /  = 0  :  straight line 
!                     zki(3,i)  -  > 0  :  node located on an arc              
!                               \  < 0  :  midpoint of the arc
!            zpz      number of nodes on the branch (including start-  and endpoint)
!            zpp      control parameter for the distribution of nodes
!                          /  = 0  :  equidistant distribution
!                     zpp  -  > 0  :  nodes are denser at the startpoint
!                          \  < 0  :  nodes are denser at the endpoint
!            zrb      Type of the boundary condition of the branch, 
!                     which can take the following values:
!                       0 -  99: Dirichlet BC, where for 1- 99: the BC has to be set in the file USERBC 
!                     200 - 299: Neumann of general BC
!                     300 - 399: visible contour, no BC
!                     400 - 499: non-visible contour, no BC
!            xbk,ybk  coordinates of key-points
!
!  Output:
!            xn, yn   coordinates of the nodes
!            kzi      node branch information
!                      kzi=0: inner node
!                      kzi<0: is identical to the keypoint with number (-kzi)
!                      kzi>0: the node is on the branch with the number (kzi)
!            p        total number of nodes
!            branchpt index of the first node of branches in brnodes 
!            brnodes  indices the nodes on the branches 
!                     (start- and endpoint of the branch are not stored in this list)
!
!  local variables
      integer (I4B) :: i, j, pold, npkt
      integer (I4B) :: nbds, k, brc, bdy
      real (DP) :: maxsize, brlength, msiz
      real (DP) :: xinnen, yinnen, rinnen
      real (DP) :: rm1, rm, zppn, x1, y1, x2, y2, dx, dy, ri, h, radius
      real (DP) :: xcenter, ycenter, area2, phi1, phi2, rho1, rho2
!
!  re-determine number of nodes on branches such that maximum meshsize if fulfilled
!  for all regions 
      call getsetting('MAXMESHSIZE',maxsize)
      do i=1,gbz
        nbds=bzil(i+1)-bzil(i)
!  for all the boundaries (outer and inner) of this region
        do k=1,nbds
          bdy=bzil(i)+k-1
!  for all branches
          do j=1,bzip(bdy+1)-bzip(bdy)
            brc=abs(bzi(bzip(bdy)+j-1))
!  determine length of branch
            x1=xbk(zki(1,brc))
            y1=ybk(zki(1,brc))
            x2=xbk(zki(2,brc))
            y2=ybk(zki(2,brc))
            if (zki(3,brc).eq.0) then
!  a line
              brlength=sqrt( (x1-x2)**2 + (y1-y2)**2 )
            else 
              if (zki(3,brc).lt.0) then
!  an arc, third point is center point
                xcenter=xbk(-zki(3,brc))
                ycenter=ybk(-zki(3,brc))
                radius=sqrt((x1-xcenter)**2+(y1-ycenter)**2)
              else
!  an arc, third point is a point at the arc, compute center point
                call inkrp(x1,y1,x2,y2,xbk(zki(3,brc)),ybk(zki(3,brc)), &
     &            xinnen,yinnen,area2,rinnen,xcenter,ycenter,radius)
              end if
!  start- and endpoint of angle
              phi1=atan2(y1-ycenter,x1-xcenter)
              phi2=atan2(y2-ycenter,x2-xcenter)
              if (phi2.lt.phi1) phi2=phi2+2._DP*pi
              brlength=radius * (phi2-phi1)
            end if
!  set number of nodes on branch according to maximum mesh-size
            msiz=min(meshsize(i),maxsize)
            zpz(brc)=max(3,zpz(brc),ceiling(1.1_DP*brlength/msiz)+1)
          end do
        end do
      end do
!
!  determine the number of points
      npkt=sum(zpz(1:gzz))-2*gzz+gkz
!  allocate arrays
      allocate(xn(npkt), yn(npkt), kzi(npkt))
!
!  copy the key-points to the node list
      do i=1,gkz
        xn(i)=xbk(i)
        yn(i)=ybk(i)
        kzi(i)=-i
      end do
      p=gkz
      allocate(branchpt(gzz+1),brnodes(npkt-gkz))
!  loop over the branches
      do i=1,gzz
        zppn=abs(zpp(i))+1._DP
        rm1=dble(zpz(i)-1)
        rm=dble(zpz(i))
!  startpoint
        x1=xbk(zki(1,i))
        y1=ybk(zki(1,i))
!  endpoint
        x2=xbk(zki(2,i))
        y2=ybk(zki(2,i))
        branchpt(i)=p+1-gkz
        pold=p
        p=p+zpz(i)-2
        kzi(pold+1:p)=i
        if (zki(3,i).eq.0) then
!  generate nodes on a staight line
          dx=x2-x1
          dy=y2-y1
          do j=2,zpz(i)-1
            ri=dble(j)
!  formula for the distribution of nodes along the branch
            h=(ri-1._DP)*(((ri-1._DP)/rm1)**zppn-1._DP)/(zppn*(ri-rm))
            if (zpp(i).ge.0) then
              xn(pold-1+j)=x1+dx*h
              yn(pold-1+j)=y1+dy*h
              brnodes(pold-1+j-gkz)=pold-1+j
            else
              xn(p+2-j)=x2-dx*h
              yn(p+2-j)=y2-dy*h
              brnodes(p+2-j-gkz)=p+2-j
            end if
          end do
        else
!  generate nodes on an arc
          if (zki(3,i).lt.0) then
!  third point is center point
            xcenter=xbk(-zki(3,i))
            ycenter=ybk(-zki(3,i))
          else
!  third point is a point at the arc, compute center point
            call inkrp(x1,y1,x2,y2,xbk(zki(3,i)),ybk(zki(3,i)),         &
     &        xinnen,yinnen,area2,rinnen,xcenter,ycenter,radius)
          end if
!  start- and endpoint radius and angle
          rho1=sqrt((x1-xcenter)**2+(y1-ycenter)**2)
          rho2=sqrt((x2-xcenter)**2+(y2-ycenter)**2)
          phi1=atan2(y1-ycenter,x1-xcenter)
          phi2=atan2(y2-ycenter,x2-xcenter)
!  assure that endpoit has a larger angle
          if (phi2.lt.phi1) phi2=phi2+2._DP*pi
          dx=phi2-phi1
          do j=2,zpz(i)-1
            ri=dble(j)
!  formula for the angular distribution of nodes along the branch
            h=(ri-1._DP)*(((ri-1._DP)/rm1)**zppn-1._DP)/(zppn*(ri-rm))
            if (zpp(i).ge.0) then
              radius=rho1+(rho2-rho1)*h
              xn(pold-1+j)=xcenter+radius*cos(phi1+dx*h)
              yn(pold-1+j)=ycenter+radius*sin(phi1+dx*h)
              brnodes(pold-1+j-gkz)=pold-1+j
            else
              radius=rho2-(rho2-rho1)*h
              xn(p+2-j)=xcenter+radius*cos(phi2-dx*h)
              yn(p+2-j)=ycenter+radius*sin(phi2-dx*h)
              brnodes(p+2-j-gkz)=p+2-j
            end if
          end do
        end if
      end do
      branchpt(gzz+1)=p+1-gkz
      write (*,2001) p
2001  format(' generated',i7,'  nodes on branches')
      call zeit( 'Generating Nodes')
      return
      end subroutine genkz
