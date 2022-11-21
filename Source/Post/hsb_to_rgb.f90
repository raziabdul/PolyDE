      subroutine hsb_to_rgb(hue,saturation,brightness,rgb)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) hue, saturation, brightness
      integer (I4B) rgb
      intent (in) :: hue, saturation, brightness
      intent (out) :: rgb
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
!  Eingabe
!    hue:         Farbwert im Bereich von 0. - 360. Grad
!    saturation:  Saettigung im Bereich von 0. - 1.
!    brightness:  Helligkeit im Bereich von 0. - 1.
!  Ausgabe
!    rgb:         integerwert jeweil 8 bit fuer B, G, R
!
!  Algorithmus nach Foley; van Dam : " Computer Graphics " (HSV_To_RGB)
      real (DP) red, green, blue, t, p, q, fracpart, h
      integer (I4B) intpart
      integer (I4B) zff, z100, z10000
      parameter (zff=Z'FF', z100=Z'100', z10000=Z'10000')
!
      if (saturation .le. 0._DP) then
        red=brightness
        green=brightness
        blue=brightness
      else 
        if (hue .ge. 360.0_DP) then
          h=0._DP
        else 
          h=hue/60._DP
        end if
        intpart=floor(h)
        fracpart=h-intpart
        p=brightness*(1._DP - saturation)
        q=brightness*(1._DP - saturation * fracpart)
        t=brightness*(1._DP - saturation * (1._DP-fracpart) )
        select case (intpart)
        case (0) 
          red=brightness
          green=t
          blue=p
        case (1)
          red=q
          green=brightness
          blue=p
        case (2)
          red=p
          green=brightness
          blue=t
        case (3)
          red=p
          green=q
          blue=brightness
        case (4)
          red=t
          green=p
          blue=brightness
        case (5)
          red=brightness
          green=p
          blue=q
        end select
      end if
!  RGB Wert im Microsoft Format umrechenen
      rgb=nint(red*zff)+z100*nint(green*zff)+z10000*nint(blue*zff)
      return
      end subroutine hsb_to_rgb
