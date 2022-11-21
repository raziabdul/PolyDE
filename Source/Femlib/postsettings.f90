module postsettings
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
!    $Revision: 1.39 $
!    $Date: 2014/08/27 15:23:49 $
!    $Author: m_kasper $
!
!  Module postsettings hold global variables used in the routines:
!  getpostsetting, putpostsetting, readpostsetting
!
      use femtypes
!
!  _c  parameters of character type
!  _i  parameters of integer type
!  _f  prameters of float type
!
!  number of the parameters of each type
      integer (I4B), parameter :: n_c=13, n_i=12, n_f=22
!  parameter values (initial settings)
      character (len=200)       :: param_c(n_c)=(/'GRID',               &   !  1
     &                                            '/WX',                &   !  2
     &                                            'FULLPLOT',           &   !  3
     &                                            'POTENTIAL',          &   !  4
     &                                            'NO',                 &   !  5
     &                                            '000000000',          &   !  6
     &                                            'NO_LOG',             &   !  7
     &                                            'NO',                 &   !  8
     &                                            'CS',                 &   !  9
     &                                            'NO',                 &   ! 10
     &                                            '',                   &   ! 11
     &                                            'UNDEFORMED',         &   ! 12
     &                                            'DATA'/)                  ! 13
!
      integer (I4B)             :: param_i(n_i)=(/100,                  &   !  1
     &                                            20,                   &   !  2
     &                                            20,                   &   !  3
     &                                            20,                   &   !  4
     &                                            20,                   &   !  5
     &                                             1,                   &   !  6
     &                                             1,                   &   !  7
     &                                            -256,                 &   !  8
     &                                            -30,                  &   !  9
     &                                             20,                  &   ! 10
     &                                             20,                  &   ! 11
     &                                              1/)                     ! 12
!
      real (DP)                 :: param_f(n_f)=(/-1._DP,               &   !  1
     &                                            1._DP,                &   !  2
     &                                            0._DP,                &   !  3
     &                                            0._DP,                &   !  4
     &                                            0._DP,                &   !  5
     &                                            0._DP,                &   !  6
     &                                            0._DP,                &   !  7
     &                                            0._DP,                &   !  8
     &                                            1._DP,                &   !  9
     &                                            0._DP,                &   ! 10
     &                                            0._DP,                &   ! 11
     &                                            0._DP,                &   ! 12
     &                                            1._DP,                &   ! 13
     &                                            0._DP,                &   ! 14
     &                                            0._DP,                &   ! 15
     &                                            0._DP,                &   ! 16
     &                                            0._DP,                &   ! 17
     &                                            0._DP,                &   ! 18
     &                                            0._DP,                &   ! 19
     &                                            1._DP,                &   ! 20
     &                                            0._DP,                &   ! 21
     &                                            1._DP/)                   ! 22
!
!  Names of the parameters
      character (len=14), parameter :: name_c(n_c)=(/'DATATYPE',        &   !  1
     &                                               'DISPLAYDEVICE',   &   !  2
     &                                               'DRAWINGTYPE',     &   !  3
     &                                               'FIELDTYPE',       &   !  4
     &                                               'INFO',            &   !  5
     &                                               'LINECOLOR',       &   !  6
     &                                               'LOGSCALE',        &   !  7
     &                                               'OUTPUT',          &   !  8
     &                                               'PALETTE',         &   !  9
     &                                               'PNUMBER',         &   ! 10
     &                                               'REGION',          &   ! 11
     &                                               'SHAPE',           &   ! 12
     &                                               'VTK_FILE'/)           ! 13
!
      character (len=14), parameter :: name_i(n_i)=(/'DIVISION',        &   !  1
     &                                               'GRID1',           &   !  2
     &                                               'GRID2',           &   !  3
     &                                               'GRIDX',           &   !  4
     &                                               'GRIDY',           &   !  5
     &                                               'LINEWIDTH',       &   !  6
     &                                               'NATURE',          &   !  7
     &                                               'NUMCOLORS',       &   !  8
     &                                               'NUMLINES',        &   !  9
     &                                               'NX',              &   ! 10
     &                                               'NY',              &   ! 11
     &                                               'NZ'/)                 ! 12
!
      character (len=14), parameter :: name_f(n_f)=(/'AMIN',            &   !  1
     &                                               'AMAX',            &   !  2
     &                                               'ENDX',            &   !  3
     &                                               'ENDY',            &   !  4
     &                                               'ENDZ',            &   !  5
     &                                               'ORIGINX',         &   !  6
     &                                               'ORIGINY',         &   !  7
     &                                               'ORIGINZ',         &   !  8
     &                                               'P1X',             &   !  9
     &                                               'P1Y',             &   ! 10
     &                                               'P1Z',             &   ! 11
     &                                               'P2X',             &   ! 12
     &                                               'P2Y',             &   ! 13
     &                                               'P2Z',             &   ! 14
     &                                               'PHI',             &   ! 15
     &                                               'STARTX',          &   ! 16
     &                                               'STARTY',          &   ! 17
     &                                               'STARTZ',          &   ! 18
     &                                               'XMIN',            &   ! 19
     &                                               'XMAX',            &   ! 20
     &                                               'YMIN',            &   ! 21
     &                                               'YMAX'/)               ! 22
!
!  fileopened will be set to true, when the file was opened 
      logical :: fileopened=.false.

end module postsettings



      subroutine getpostsetting_c(key,value)
!
!  Return the string value which corresponds to the parameter key
!  in the POSTsetting.txt file
!
!  Specific routine to the genereic routine getsetting
!  for the case of a string argument.
!
!    key     one of the Parameter-keywords
!    value   string corresponding to the Parameter setting
!  
      use feminterface, only:
      use postsettings
      use femtypes
      implicit none
      character (len=*) :: key
      character (len=*) :: value
      intent (in) :: key
      intent (out) :: value
!  local variables
      integer (I4B) :: i
!
      if (.not. fileopened) then
        print*,'Error in getpostsetting_c: file postsettings not opened'
        print*,'key:', key
        return
      end if
      do i=1,n_c
        if (key .eq. name_c(i)) then
          value=param_c(i)
          return
        end if
      end do 
!  key not found 
      print*,' Error (ignored) in Parameter setting: "',key,            &
     &  '" unknown parameter'
      return
      end subroutine getpostsetting_c
!
!
!
      subroutine getpostsetting_i(key,value)
!
!  Return the integer value which corresponds to the parameter key
!  in the POSTsetting.txt file
!
!  Specific routine to the genereic routine getsetting
!  for the case of a integer argument.
!
!    key     one of the Parameter-keywords
!    value   integer corresponding to the Parameter setting
!  
      use feminterface, only:
      use postsettings 
      use femtypes
      implicit none
      character (len=*) :: key
      integer (I4B) :: value
      intent (in) :: key
      intent (out) :: value
!  local variables
      integer (I4B) :: i
!
      if (.not. fileopened) then
        print*,'Error in getpostsetting_i: file postsettings not opened'
        return
      end if
      do i=1,n_i
        if (key .eq. name_i(i)) then
          value=param_i(i)
          return
        end if
      end do 
!  key not found 
      print*,' Error (ignored) in Parameter setting: "',key,            &
     &  '" unknown parameter'
      return
      end subroutine getpostsetting_i
!
!
!
      subroutine getpostsetting_f(key,value)
!
!  Return the float value which corresponds to the parameter key
!  in the POSTsetting.txt file
!
!  Specific routine to the genereic routine getsetting
!  for the case of a float argument.
!
!    key     one of the Parameter-keywords
!    value   float corresponding to the Parameter setting
!  
      use feminterface, only:
      use postsettings 
      use femtypes
      implicit none
      character (len=*) :: key
      real (DP) :: value
      intent (in) :: key
      intent (out) :: value
!  local variables
      integer (I4B) :: i
!
      if (.not. fileopened) then
        print*,'Error in getpostsetting_f: file postsettings not opened'
        return
      end if
      do i=1,n_f
        if (key .eq. name_f(i)) then
          value=param_f(i)
          return
        end if
      end do 
!  key not found 
      print*,' Error (ignored) in Parameter setting: "',key,            &
     &  '" unknown parameter'
      return
      end subroutine getpostsetting_f
!
!
!
      subroutine putpostsetting_c(key,value)
!
!  Set the string value which corresponds to the parameter key
!  in the POSTsetting.txt file
!
!  Specific routine to the genereic routine putsetting
!  for the case of a string argument.
!
!    key     one of the Parameter-keywords
!    value   string corresponding to the Parameter setting
!  
      use feminterface, only:
      use postsettings
      use femtypes
      implicit none
      character (len=*) :: key
      character (len=*) :: value
      intent (in) :: key
      intent (in) :: value
!  local variables
      integer (I4B) :: i
!
      do i=1,n_c
        if (key .eq. name_c(i)) then
          if (len_trim(value) .gt. len(param_c)) then
            print*,'Error in parameter storing, string to long:',value
          end if
          param_c(i)=value
          return
        end if
      end do 
!  key not found 
      print*,' Error (ignored) in Parameter setting: "',key,            &
     &  '" unknown parameter'
      return
      end subroutine putpostsetting_c
!
!
!
      subroutine putpostsetting_i(key,value)
!
!  Set the integer value which corresponds to the parameter key
!  in the POSTsetting.txt file
!
!  Specific routine to the genereic routine putsetting
!  for the case of a integer argument.
!
!    key     one of the Parameter-keywords
!    value   integer corresponding to the Parameter setting
!  
      use feminterface, only:
      use postsettings 
      use femtypes
      implicit none
      character (len=*) :: key
      integer (I4B) :: value
      intent (in) :: key
      intent (in) :: value
!  local variables
      integer (I4B) :: i
!
      do i=1,n_i
        if (key .eq. name_i(i)) then
          param_i(i)=value
          return
        end if
      end do 
!  key not found 
      print*,' Error (ignored) in Parameter setting: "',key,            &
     &  '" unknown parameter'
      return
      end subroutine putpostsetting_i
!
!
!
      subroutine putpostsetting_f(key,value)
!
!  Set the float value which corresponds to the parameter key
!  in the POSTsetting.txt file
!
!  Specific routine to the genereic routine putsetting
!  for the case of a float argument.
!
!    key     one of the Parameter-keywords
!    value   float corresponding to the Parameter setting
!  
      use feminterface, only:
      use postsettings 
      use femtypes
      implicit none
      character (len=*) :: key
      real (DP) :: value
      intent (in) :: key
      intent (in) :: value
!  local variables
      integer (I4B) :: i
!
      do i=1,n_f
        if (key .eq. name_f(i)) then
          param_f(i)=value
          return
        end if
      end do 
!  key not found 
      print*,' Error (ignored) in Parameter setting: "',key,            &
     &  '" unknown parameter'
      return
      end subroutine putpostsetting_f
!
!
!
      subroutine readpostsetting(callback)
!
!  Read the POSTsetting.txt file of parameters
!
!  Kown bugs: the file is opened from the actual directory
!
      use feminterface, only: low2hi, string2number, getsetting
      use postsettings 
      use femtypes
      implicit none
      character (len=*):: callback
!  local variables
      integer (I4B) :: i, pos, epos, ios, rpos, length
      integer (I4B), save :: unitid
      character (len=200) :: str, name
      logical :: exist,search
      character (len=200), save :: postsettingfile
!
!  look for file
      if (.not.fileopened) then
        call getsetting('PROJECTPATH',postsettingfile)
        postsettingfile = postsettingfile(1:len_trim(postsettingfile))//'POSTsettings.txt'
!
        call grglun(unitid)
        inquire (file=postsettingfile,exist=exist)
        if (.not. exist) then
           print*,'Error in readpostsetting: file ',postsettingfile(1:len_trim(postsettingfile)),&
                 ' not present.'
          return
        end if
        open (unitid,file=postsettingfile,form='FORMATTED',iostat=ios)
        if (ios .ne. 0) then
          print*,'Error in opening POSTsettings.txt'
          return
        end if
        fileopened=.true.
      end if
!
!  Line by line loop over all entries in file
      do
        read(unitid,fmt='(a)',iostat=ios) str
        if (ios .lt. 0) then
          callback = 'EOF'
          close(unitid)
          return
        else if (ios .gt. 0) then 
          callback = 'ERROR'
          print*,'Error in readpostsetting: string not read successfully from '&
                 ,postsettingfile(1:len_trim(postsettingfile))
          close(unitid)
          return
        end if
!
!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  replace tabs (char(9)) by a single blank
        do i=1,len_trim(str)
          if (str(i:i) .eq. char(9)) str(i:i)=' '
        end do
!  allow for blank or comment lines
        length=len_trim(str)
        if (length .eq.0) cycle
!
        print*,str(1:length)
        str=adjustl(str)
!  Loop over all parameters of string type
        epos=index(str,'=')
        call low2hi(str,length)
        search=.true.
        i=1
        do while (search .and. i.le.n_c)
          name=name_c(i)
          pos=index(str,name(1:len_trim(name)))
          rpos=index(name(1:len_trim(name)),str(1:len_trim(str(1:epos-1))))
          if (pos .gt. 0 .and. rpos .gt. 0) then
            if (epos .gt. pos) then
              str=str(epos+1:)
              param_c(i)=adjustl(str)
              search=.false.
            end if
          end if
        i=i+1
        end do
!  Loop over all parameters of integer type
        i=1
        do while (search .and. i.le.n_i)
          name=name_i(i)
          pos=index(str,name(1:len_trim(name)))
          rpos=index(name(1:len_trim(name)),str(1:len_trim(str(1:epos-1))))
          if (pos .gt. 0 .and. rpos .gt. 0) then
            if (epos .gt. pos) then
!  Convert string to integer
              call string2number(str(epos+1:),param_i(i))
              search=.false.
            end if
          end if
        i=i+1
        end do
!  Loop over all parameters of float type
        i=1
        do while (search .and. i.le.n_f)
          name=name_f(i)
          pos=index(str,name(1:len_trim(name)))
          rpos=index(name(1:len_trim(name)),str(1:len_trim(str(1:epos-1))))
          if (pos .gt. 0 .and. rpos .gt. 0) then
            if (epos .gt. pos) then
!  Convert string to float
              call string2number(str(epos+1:),param_f(i))
              search=.false.
            end if
          end if
        i=i+1
        end do
!  Unknown Parameter found in file
        if (search) then
          print*
          callback = str
          return
        end if
      end do
!
      close(unitid)
      return
      end subroutine readpostsetting
