      subroutine farbpalette(palette,ncolor,colorscheme)
      use feminterface, only : hsb_to_rgb
      use femtypes
      implicit none
      integer (I4B) ncolor
      integer (I4B) palette(:)
      character (len=*) :: colorscheme
      intent (in) ::ncolor, colorscheme
      intent (out):: palette
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
!    $Date: 2012/01/03 15:44:36 $
!    $Author: m_kasper $
!
!  Determine a color palette
!  Input:
!    ncolor      Number of colors
!    color       is  .false.  if a gray scale is required
!  Output:
!    palette     Color value as a RGB-value having each 8 bit for B, G, R
!
!  local variables
      real (DP) hue, saturation, brightness, r, g, b
      integer (I4B) i,rgb
!
      select case(colorscheme)
      case ('CS','STANDARD','COLORSCALE','COLOR')
        saturation=.9_DP
        brightness=1.0_DP
        do i=0,ncolor-1
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(ncolor-i)=rgb
        end do
      case ('DEEP PERCEPTION')
        do i=0,ncolor-1
          r = 255._DP * real(i)/real(ncolor-1)
!          r = 255._DP * (1._DP - real(i)/real(ncolor-1))
          g = 255._DP * real(i)/real(ncolor-1)
!          g = 255._DP * (1._DP - real(i)/real(ncolor-1))
!          g = 255
!          b = 255._DP * real(i)/real(ncolor-1)
          b = 255._DP * (1._DP - real(i)/real(ncolor-1))
!b=0
          rgb=nint(r) + 2**8*nint(g) + 2**16*nint(b)
! magenta 101 yellow 110 cyan 011
          palette(i+1)=rgb
        end do
      case ('REVERSE COLOR')
        saturation=.9_DP
        brightness=1.0_DP
        do i=1,ncolor
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(i)=rgb
        end do
      case ('PALE COLOR')
        saturation=.6_DP
        brightness=1.0_DP
        do i=0,ncolor-1
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(ncolor-i)=rgb
        end do
      case ('REVERSE PALE COLOR')
        saturation=.6_DP
        brightness=1.0_DP
        do i=1,ncolor
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(i)=rgb
        end do
      case ('PASTEL COLOR')
        saturation=.3_DP
        brightness=1.0_DP
        do i=0,ncolor-1
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(ncolor-i)=rgb
        end do
      case ('REVERSE PASTEL COLOR')
        saturation=.3_DP
        brightness=1.0_DP
        do i=1,ncolor
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(i)=rgb
        end do
      case ('MILD COLOR')
        saturation=1._DP
        brightness=1.0_DP
        do i=0,ncolor-1
          hue=30._DP+150._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(ncolor-i)=rgb
        end do
      case ('REVERSE MILD COLOR')
        saturation=1._DP
        brightness=1.0_DP
        do i=1,ncolor
          hue=30._DP+150._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(i)=rgb
        end do
      case ('SCREAMING COLOR')
        saturation=1._DP
        brightness=1.0_DP
        do i=1,ncolor
          hue=240._DP+120._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(i)=rgb
        end do
      case ('REVERSE SCREAMING COLOR')
        saturation=1._DP
        brightness=1.0_DP
        do i=0,ncolor-1
          hue=240._DP+120._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(ncolor-i)=rgb
        end do
      case ('DULL COLOR')
        saturation=.9_DP
        brightness=0.7_DP
        do i=0,ncolor-1
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(ncolor-i)=rgb
        end do
      case ('REVERSE DULL COLOR')
        saturation=.9_DP
        brightness=0.7_DP
        do i=1,ncolor
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(i)=rgb
        end do
      case ('DEEP COLOR')
        saturation=1.0_DP
        brightness=0.85_DP
        do i=0,ncolor-1
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(ncolor-i)=rgb
        end do
      case ('REVERSE DEEP COLOR')
        saturation=1.0_DP
        brightness=0.85_DP
        do i=1,ncolor
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(i)=rgb
        end do
      case ('GS','GRAYSCALE','GREYSCALE')
        saturation=0._DP
        hue=0._DP
        do i=0,ncolor-1
          brightness=real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(ncolor-i)=rgb
        end do
      case ('REVERSE GS','REVERSE GRAYSCALE','REVERSE GREYSCALE')
        saturation=0._DP
        hue=0._DP
        do i=1,ncolor
          brightness=real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(i)=rgb
        end do
      case default  ! same as standard
        saturation=.9_DP
        brightness=1.0_DP
        do i=0,ncolor-1
          hue=240._DP*real(i)/real(ncolor-1)
          call hsb_to_rgb(hue,saturation,brightness,rgb)
          palette(ncolor-i)=rgb
        end do
      end select
      return
      end subroutine farbpalette
