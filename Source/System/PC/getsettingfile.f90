      subroutine getsettingfile(settingfile)
      use dfport
!      use ifport
      implicit none
      character (len=*) :: settingfile
      intent (out) :: settingfile
!
!     PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 Institute for Micro Systems Technology,
!                       Hamburg University of Technology.
!
!    This file is part of PolyDE.
!
!    PolyDE is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation; either version 2, or (at your
!    option) any later version.
!
!    PolyDE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!    MA 02110-1301 USA.
!
!    $Revision: 1.10 $
!    $Date: 2015/04/01 09:55:30 $
!    $Author: juryanatzki $
      
!  Subroutine to read the PolydePath from the environment settings and return
!  the settingfile to the calling program
!
!  local variables:
      character (len=200) :: polydepath
!
      call getenv("PolydePath",polydepath)
      !print *,"PolydePath specified in the user's environment is: ",polydepath(1:len_trim(polydepath))
      !settingfile = polydepath(1:len_trim(polydepath))//"FEMSettings.txt"
      settingfile = polydepath(1:len_trim(polydepath))//"FEMProject.json"
      end subroutine getsettingfile
