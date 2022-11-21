      subroutine fpotentialp(nc,nl,phi,amax,amin,colorfill,hlin,        &
     &                       r,g,b,fieldtype)
      use feminterface, only: farbpalette, nfarbep, nsurf, getpostsetting
      use femtypes
      implicit none
      integer (I4B) :: nc, nl
      real (SP) :: r, g, b
      real (DP) :: phi, amax, amin
      logical :: colorfill, hlin
      character (len=*) :: fieldtype
      intent (in):: nc, nl, phi, amax, amin, colorfill, hlin
      intent (in):: r, g, b, fieldtype
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
!    $Revision: 1.23 $
!    $Date: 2014/06/24 12:58:35 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
!
!  Coloring of an element (triangle) in dependence of a potential value
!       nc          number of colors
!       nl          number of lines for equipotential lines
!       phi         phase angle in degrees (t*omega)
!       amin, amax  minimum and maximum potential value
!       colorfill   =.true, for color fill plot
!       hlin        =.true. if equipotential should be plotted
!       r, g, b     red green and blue value for the potential lines
!       fieldtype   defines which kind of field quantity should be plotted

      real (DP) zcolmin, zcolmax, zlinmin, zlinmax, amaxl, aminl
      integer (I4B), allocatable :: palette(:)
      integer (I4B) ncolor, nlin
      character (len=30) :: colorscale
!
      amaxl=amax
      aminl=amin
      write (*,765) amaxl
765   format('   maximum potential to display ',g13.6)
      write (*,766) aminl
766   format('   minimum potential to display ',g13.6)
!
      ncolor=abs(nc)
      allocate (palette(ncolor))
!  set color palette e.g. gray scale (GS) or color scale (CS)
      call getpostsetting('PALETTE',colorscale)
      zcolmin=aminl
      zcolmax=amaxl
      call farbpalette(palette,ncolor,colorscale)
!
!  equipotential scale
      nlin=abs(nl)
      zlinmin=aminl
      zlinmax=amaxl
!
      call nfarbep(palette,colorfill,ncolor,zcolmin,zcolmax,            &
     &  hlin,nlin,zlinmin,zlinmax,r,g,b,phi,fieldtype)
!      call nsurf(bb,x,y,p,e,ie,palette,ncolor,zcolmin,zcolmax,          &
!     &  10._DP,.true.,.false.,.true.,0.001_DP,50._DP,30._DP)
      deallocate (palette)
      end subroutine fpotentialp