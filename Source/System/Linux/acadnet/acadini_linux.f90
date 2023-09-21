      subroutine acadini
      use feminterface, only: getfile, getsetting, putsetting, wpathsetting
      use feminterface, only: readsetting, physicsinfo, zanfpp
      use feminterface, only: setprocesspriority, zeit
      use globalvariables, only: nnat
      use femtypes
! use posix library
      use ifposix
! use ifort
      use ifport

      implicit none
!
!  Open DXF-file and do some setting
!  This is linux version of acadini.f90
!
!  local variables
      integer (I4B) :: ios, bildnr, priority, res, dim = 2, status, strlen
      real (DP) :: xmin, xmax, ymin, ymax, h
      logical exist
      character (len=201) :: path, filena, buf, errmsg
      character (len=25) physics
!
!--------------------------------------------------------
! Linux-centric lines using intel compiler function
! to get argument from linux/unix shell
!--------------------------------------------------------
     call get_command_argument(1, filena, strlen, status, errmsg) 
! get full path of the file
      res = fullpathqq(filena, buf)
      if (res-strlen > 0) then
         path=buf(1:res-strlen)
      else
         path=''
      end if
! 
10   if (status .ne. 0) then
        print*
        print*,' enter the DXF-file containing the structure'
        print*,' Filename:'
        print*
        call getfile(filena, path)
     end if
!  get the filename and open file
      inquire (file=path(1:len_trim(path))//filena,exist=exist)
      if ( .not. exist) then
        print*,' the file: ',path(1:len_trim(path))//filena,' does not exist'
        print*,' please try again'
        goto 10
      end if
!  set the projectpath and read the setting file from the new project path 
      call putsetting('PROJECTPATH',path)
      call getsetting('PHYSICS_MODE',physics)
!  get the number of natures
      call physicsinfo(physics,nnat,dim)
!  open the DXF-file
      open (38,file=path(1:len_trim(path))//filena,action='READ',iostat=ios)
      if (ios .ne. 0) then
        print*,'*ERROR opening file: ',path//filena
        print*,' Error no: ',ios
      end if
!  initialize the graphics window
!  Driver is NULL for initialization purpose only
      xmin=0._DP
      xmax=1._DP
      ymin=0._DP
      ymax=1._DP
      call zanfpp('/NULL',xmin,xmax,ymin,ymax,h,bildnr)
      call zeit(' ')
!
      priority=1
      call setprocesspriority(priority)
!
      return
      end
