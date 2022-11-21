      subroutine sorttr(gbz,gzz,bzi,bzip,bzil,zki,zrb,ok)
      use feminterface, only: zeit
      use femtypes
      implicit none
      integer (I4B) :: gbz,gzz
      integer (I4B) :: bzi(:), bzip(:), bzil(:), zki(:,:), zrb(:,:)
      logical :: ok
      intent (in) :: gbz, gzz, bzip, zki, zrb
      intent (out) :: ok
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
!    $Date: 2008/12/22 15:06:23 $
!    $Author: m_kasper $
!
!  Input:
!            gbz      number of regions
!            gzz      number of branches
!            bzip     Pointer (index) to the start branch of a region in the
!                     list bzi (for compact storage)
!            bzil     Pointer (index) onto bzip including all regions 
!                     which are inside a region (multiply connected regions)
!            zki      list of of key-points (geometry nodes) which determine the branches
!                     zki(1,i)  start key-point of the branch
!                     zki(2,i)  end key-point of the branch
!                     zki(3,i)  third key-koint of the branch, with
!                               /  = 0  :  straight line 
!                     zki(3,i)  -  > 0  :  node located on an arc
!                               \  < 0  :  midpoint of the arc
!            zrb      Type of the boundary condition of the branch, 
!                     which can take the following values:
!                       0 -  99: Dirichlet BC, where for 1- 99: the BC has to be set in the file USERBC
!                     200 - 299: Neumann of general BC
!                     300 - 399: visible contour, no BC
!                     400 - 499: non-visible contour, no BC
!
!  Output:
!            ok       =.true. if no error occured
!
!  In-/ Output:
!            bzi      list of braches which form the regions (compact storage)
!
!  Sort all branches according to their direction
!  (intended for the case that the input is not fully correct)
!
!  local variables
      integer (I4B) :: i, ii, j, lpunkt, apunkt, zweig
      integer (I4B) :: liste(gzz)
!
      liste(:)=0
      ok=.true.
!  for all regions 
      do ii=1,gbz
!  for all branches of the region
        do i=bzil(ii),bzil(ii+1)-1
          zweig=bzi(bzip(i))
          if (zweig.gt.0) then
            apunkt=zki(1,zweig)
            lpunkt=zki(2,zweig)
          else
            apunkt=zki(2,-zweig)
            lpunkt=zki(1,-zweig)
          end if
!  count how often the branch is used
          liste(abs(zweig))=liste(abs(zweig))+1
          do j=2,bzip(i+1)-bzip(i)
            zweig=bzi(bzip(i)+j-1)
            if (zweig.gt.0) then
              if (lpunkt.eq.zki(2,zweig)) then
                lpunkt=zki(1,zweig)
!  correct the direction (invert)
                bzi(bzip(i)+j-1)=-zweig
              else
                if (.not.(lpunkt.eq.zki(1,zweig))) then
!  neither start- nor endpoint match
                  write (*,1000) ii,bzi(zweig)
1000              format(' ERROR in region: ',i4,' branch',i5,' is incorrect')
                  ok=.false.
                  return
                else
                  lpunkt=zki(2,zweig)
                end if
              end if
            else
              if (lpunkt.eq.zki(1,-zweig)) then
                lpunkt=zki(2,-zweig)
!  correct the direction (invert)
                bzi(bzip(i)+j-1)=-zweig
              else
                if (.not.(lpunkt.eq.zki(2,-zweig))) then
!  neither start- nor endpoint match
                  write (*,1000) i,zweig
                  ok=.false.
                  return
                else
                  lpunkt=zki(1,-zweig)
                end if
              end if
            end if
!  count how often the branch is used
            liste(abs(zweig))=liste(abs(zweig))+1
          end do
          if (lpunkt.ne.apunkt) then
!  the last branch does not close the region
            write (*,1001) ii
1001        format(' ERROR: region ',i4,' is not closed')
            ok=.false.
            return
          end if
        end do
      end do
!  now check the correct use of branches
      do i=1,gzz
        if (liste(i).eq.0) then
!  each branch should be used at least once
          write (*,1002) i
1002      format(' ERROR: branch ',i6,' is not used')
          ok=.false.
        end if
        if ((liste(i).eq.1).and.(any(zrb(i,:).ge.300))) then
!  inner branches should be used twice
          write (*,1003) i
1003      format(' ERROR: branch ',i6,                                  &
     &      ' is labeled to be inner but is used only once')
          ok=.false.
        end if
        if ((liste(i).eq.2).and.(any(zrb(i,:).lt.300))) then
!  boundary branches should be used only once 
          if (all(zrb(i,:).lt.200)) then
!  supply warning in the case of Dirichlet BC
            write (*,1004) i
1004        format(' WARNING: branch ',i6,' having a Dirichlet',        &
     &      'boundary condition is used more than once')
          else
!  supply error message in the case of Neumann or general BC
            write (*,1005) i
1005        format(' ERROR: branch ',i6,' having a Neumann or general ',&
     &      'boundary condition is used more than once')
            ok=.false.
          end if
        end if
        if (liste(i).gt.2) then
!  none of the branches should be used more than twice
          write (*,1006) i
1006      format(' ERROR: branch ',i6,' is used more than twice')
          ok=.false.
        end if
      end do
      call zeit(' Sorting')
      return
      end subroutine sorttr
