      subroutine acadini
      use feminterface, only: getfile, getsetting, putsetting, wpathsetting
      use feminterface, only: readsetting, physicsinfo, zanfpp
      use feminterface, only: setprocesspriority
      use globalvariables, only: nnat
      use femtypes
      implicit none
!
!  Open DXF-file and do some setting
!
!  local variables
      integer (I4B) :: ios, bildnr, priority
      real (DP) :: xmin, xmax, ymin, ymax, h
      logical exist
      character (len=201) :: path, filena
      character (len=25) physics, system
!
      print*
      print*,' enter the DXF-file containing the structure'
      print*,' Filename:'
10    print*
!  get the filename and open file
      call getfile(filena, path)
      inquire (file=path(1:len_trim(path))//filena,exist=exist)
      if ( .not. exist) then
        print*,' the file: ',path(1:len_trim(path))//filena,' does not exist'
        print*,' please try again'
        goto 10
      end if
!  set the projectpath and read the setting file from the new project path 
      call getsetting('WHATSYSTEM',system)
      call putsetting('PROJECTPATH',path)
      if (system.eq.'WINDOWS') then
        call wpathsetting
      end if
      call readsetting
      call getsetting('PHYSICS_MODE',physics)
!  get the number of natures
      call physicsinfo( physics, nnat, 2 )
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
