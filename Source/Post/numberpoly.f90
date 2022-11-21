      subroutine numberpoly(dch,xdmin,xdmax,ydmin,ydmax,nature)
      use feminterface, only: flaech, schnit
      use femtypes
      use globalvariables
      implicit none
      integer (I4B) :: nature
      real (DP) :: dch, xdmin, xdmax, ydmin, ydmax
      intent (in) :: dch, xdmin, xdmax, ydmin, ydmax, nature
!
!------------------------------------------------------------------------------
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
!    $Date: 2008/08/20 12:56:06 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
!
!  Numberpoly draws a number using pgplot routines for the element's polynomial
!  degree. The size of the font is chosen in relation to the element's incircle
!  which touches the three element edges. See wikipedia entry for incircle geo-
!  metry and used terms (http://en.wikipedia.org/wiki/Incircle).
!
!  Input:
!    dch        default character height (minimum of geometry's height or width
!               divided by 40)
!
!  NOTE: pgptxt in Linux can only take real argumements in single precision
!         seems the same for Windows too? -- Razi
!
!  Internal variables:
!
      integer (I4B) :: elem, i, anfang, ende
      integer (I4B) :: inach(3), ivorg(3)
      real (SP) :: fontsize, point(2)
      real (DP) :: angle(3), area, inradius, length(3), t1, t2
      real (DP) :: endp(2), startp(2), sides(3)
      real (DP) :: vec(2,3), icvec(2,3), tempangle
      real (DP) :: from1(2), from2(2), to1(2), to2(2)
      character (len=2) :: polychar(20)
      logical :: ok
      parameter (inach=(/2,3,1/))
      parameter (ivorg=(/3,1,2/))
      parameter (polychar=(/'1','2','3','4','5','6','7','8','9','10','11','12','13','14',&
                            '15','16','17','18','19','20'/))
!
      do elem = 1,n
!  calculate the element's area
        area = flaech(xn(e(1,elem)),yn(e(1,elem)),xn(e(2,elem)),yn(e(2,elem)),&
                        xn(e(3,elem)),yn(e(3,elem)))
        do i = 1,3
!  set start- and endpoints for direction vector (along element edge)
          anfang = e(inach(i),elem)
          ende   = e(ivorg(i),elem)
!  set start- and endpoints for direction vector (along element edge)
          startp = (/xn(anfang),yn(anfang)/)
          endp = (/xn(ende),yn(ende)/)
!  compute length of side i
          vec(:,i) = endp - startp
          sides(i) = sqrt(vec(1,i)**2+vec(2,i)**2)
          vec(:,i) = vec(:,i) / sides(i)
        end do
!
!  compute inradius
        inradius = 2*area / sum(sides)
!  set fontsize and call pgsch
        fontsize = inradius / dch
        call pgsch(fontsize)
!  compute incenter or Gergonne point
        do i = 1,3
!  determine angle at vertex 1
          angle(i) = acos( ( vec(1,ivorg(i))*(-vec(1,inach(i))) + vec(2,ivorg(i))*&
                           (-vec(2,inach(i))) ) )
!  determine length of bisecting line
          length(i) = 2._DP*sides(ivorg(i))*sides(inach(i))*cos(angle(i)/2._DP) / &
                      (sides(ivorg(i)) + sides(inach(i)))
!  determine vector that points towards the incenter
          tempangle = atan2(vec(2,ivorg(i)),vec(1,ivorg(i)))
          icvec(1,i) = length(i)*cos(tempangle+angle(i)/2._DP)
          icvec(2,i) = length(i)*sin(tempangle+angle(i)/2._DP)
        end do
!  determine coordinates of incenter
        from1=(/xn(e(1,elem)),yn(e(1,elem))/)
        from2=(/xn(e(2,elem)),yn(e(2,elem))/)
        to1=(/xn(e(1,elem))+icvec(1,1),yn(e(1,elem))+icvec(2,1)/)
        to2=(/xn(e(2,elem))+icvec(1,2),yn(e(2,elem))+icvec(2,2)/)
!
        call schnit(from1(1),from1(2),to1(1),to1(2),from2(1),from2(2),to2(1),to2(2),t1,t2,ok)
!  draw the text
        if (ok) then
         point(1) = from1(1)+t1*(to1(1)-from1(1))
         point(2) = from1(2)+t1*(to1(2)-from1(2))-inradius*0.5_SP
         if ( (point(1) .lt. xdmin) .or. (point(1) .gt. xdmax) .or. &
              (point(2) .lt. ydmin) .or. (point(2) .gt. ydmax) ) then
           cycle
         else
           call pgptxt(point(1),point(2),0.0,0.5,polychar(ep(elem,nature)))
         end if
        else
          print*,"An error occured for computing the incenter of element:",elem
        end if
      end do
!
      end subroutine numberpoly
