!
      subroutine copyr(text1,text2)
      use feminterface, only: revis
      use femtypes
      implicit none
      character (len=*) text1,text2
      intent (in) :: text1, text2
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
!    $Date: 2006/07/07 21:08:35 $
!    $Author: r_abdul $
!
! Dieses Unterprogramm gibt einen allgemeinen Copyright Text aus
!
      write (*,1)
1     format(' ')
! Leerzeile
      write (*,2) text1
2     format('  ',a,':')
      write (*,3) text2
3     format('  ',a)
      write (*,1)
      call revis()
      return
      end
!
!
!
      subroutine revis()
      use feminterface, only:
      use femtypes
      implicit none
!
      write (*,1000)
1000  format('  Interner Revisionsstand: 3.16 10/02 ')
      end subroutine
