subroutine getfile(file,path)
!------------------------------------------------
!             This is for Linux
!-----------------------------------------------
!     Author: $Author: r_abdul $ 
!       Date: $Date: 2006/07/08 00:56:21 $
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
!   Revision: $Revision: 1.13 $
!-----------------------------------------------
! The same name as in getfile.f90 for PC
!
! 
  use feminterface, only:
  use femtypes
  use ifport
  implicit none
  character (len=*) file, path
  character (256) buf
  integer (I4B) res, lenf
  intent (out) file, path

  read (*,*) file
  lenf =len_trim(file)
! get full path file name
  res = fullpathqq(file, buf)
  if (res-lenf .gt. 0) then
     path=buf(1:res-lenf)
  else
     path=''
  end if

end subroutine getfile

