      subroutine convertcolorstring(red_f, green_f, blue_f)
      use feminterface, only : getpostsetting, string2number
      use femtypes
      implicit none
      real (SP) :: red_f, green_f, blue_f
      intent (out) :: red_f, green_f, blue_f
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
!    $Revision: 1.4 $
!    $Date: 2006/07/07 22:13:33 $
!    $Author: r_abdul $
!
!  Program to read a colorstring and convert it to float. So it can
!  be used for drawing the mesh and structure in different colors
!
!  local variables
      character (len=9) :: colorstring
      character (len=3) :: red, green, blue
!
!  fetch colorstring from POSTsettings.txt
      call getpostsetting('LINECOLOR',colorstring)
!
!  store different parts of the colorstring under red, green and blue
      red = colorstring(1:3)
      green = colorstring(4:6)
      blue = colorstring(7:9)
!
      call string2number(red, red_f)
      call string2number(green, green_f)
      call string2number(blue, blue_f)
      red_f = red_f / 255
      green_f = green_f / 255
      blue_f = blue_f / 255
!
      end subroutine convertcolorstring
