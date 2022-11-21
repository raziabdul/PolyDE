module femsettings
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
!    $Revision: 1.55 $
!    $Date: 2015/10/27 14:06:48 $
!    $Author: juryanatzki $
!
!
!  Module femsettings hold global variables used in the routines:
!  getsetting, putsetting, readsetting, writsetting
!
      use femtypes
!
!  _c  parameters of character type
!  _i  parameters of integer type
!  _f  prameters of float type
!
!  number of the parameters of each type
      integer (I4B), parameter :: n_c=28, n_i=9, n_f=31
!  parameter values (initial settings)
      character (len=200)       :: param_c(n_c)=(/                      &
     &                                            'HP_ADAPT',           &  !  ADAPTION_TYPE
     &                                            'LCSR',               &  !  CSRFORMAT
     &                                            'NEDELEC',            &  !  ELEMENT_TYPE
     &                                            'FULLRESIDUAL',       &  !  ERROR_ESTIMATOR
     &                                            'KEYPOINT',           &  !  HP_ALGORITHM
     &                                            'SSORCG',             &  !  LINSOLVERTYPE
     &                                            'BINARY',             &  !  MATRIXOUTPUT
     &                                            'NONSYMMETRIC',       &  !  MATRIXTYPE
     &                                            'BASIS.UNV',          &  !  MESHFILE
     &                                            'YES',                &  !  MESHSMOOTHING
     &                                            'TO',                 &  !  MP_COUPLING
     &                                            'NO',                 &  !  MULTIPHYSICS
     &                                            'NO',                 &  !  NONLINEAR
     &                                            'FIXED_POINT',        &  !  NONLINSOLVERTYPE
     &                                            'LAPLACE',            &  !  PHYSICS_MODE
     &                                            '??????',             &  !  PROJECTPATH
     &                                            'LAPLACE',            &  !  SMOOTHING
     &                                            'NEWMARK',            &  !  TRANSIENTSOLVER
     &                                            'WINDOWS',            &  !  WHATSYSTEM
     &                                            'NO',                 &  !  WRITE_ITERATION
     &                                            'NO',                 &  !  WRITEFIGURE
     &                                            'NO',                 &  !  WRITEMATRIX
     &                                            'NO',                 &  !  WRITE_SIMSTATS
     &                                            'YES',                &  !  MESHQUALITYPLOT
     &                                            'NONACCUMULATIVE',    &  !  MARKING_TYPE_H
     &                                            'ACCUMULATIVE',       &  !  MARKING_TYPE_P
     &                                            'SING_VAL',           &  !  ADAPT_CRITERION_H
     &                                            'NO'                  &  !  MATER_VARIATION
     &                                               /)
      integer (I4B)             :: param_i(n_i)=(/                      &
     &                                            25,                   &  !  ADAPTATION_STEPS
     &                                            100000000,            &  !  maximum execution time (runtime) in seconds
     &                                            100000000,            &  !  maximum number of degrees of freedom
     &                                            100000000,            &  !  maximum number of nodes (triangle vetrices)
     &                                            1000,                 &  !  NONLINITERSTEPS
     &                                            2,                    &  !  initial poynomial order
     &                                            1,                    &  !  DELTA_POLYDEG
     &                                            1,                    &  !  HPF8R_ORDER
     &                                            1000                  &  !  MAX_ANDERSON_ITER
     &                                               /)
      real (DP)                 :: param_f(n_f)=(/                      &
     &                                            1._DP,                &  !   1  ADAPTION_ERROR
     &                                            0.25_DP,              &  !   2  ALPHA_NEWMARK
     &                                            0.5_DP,               &  !   3  DELTA_NEWMARK
     &                                            100._DP,              &  !   4  END_TIME
     &                                            1.e-7_DP,             &  !   5  EPS
     &                                            1._DP,                &  !   6  FPRELAX_PARAMETER
     &                                            -9999._DP,            &  !   7  FREQUENCY
     &                                            1.0_DP,               &  !   8  GEOMETRY_FACTOR
     &                                            0.001_DP,             &  !   9  INIT_TIME_STEP
     &                                            1.e-9_DP,             &  !  10  LINSOLVER_ERROR
     &                                            1.e32_DP,             &  !  11  MAXMESHSIZE
     &                                            1.e-5_DP,             &  !  12  NONLIN_ERROR
     &                                            25._DP,               &  !  13  SKINNYLIMIT
     &                                            0._DP,                &  !  14  START_TIME
     &                                            0.5_DP,               &  !  15  TIME_STEP
     &                                            1.0_DP,               &  !  16  TRAN_GTOL
     &                                            0.2_DP,               &  !  17  TRAN_LTOL_LO
     &                                            0.01_DP,              &  !  18  TRAN_LTOL_UP
     &                                            0.01_DP,              &  !  19  TRAN_REL_TOLER
     &                                            0.5_DP,               &  !  20  ADAPT_FRACTION
     &                                            0.1_DP,               &  !  21  MAX_ERRORBOUND_H
     &                                            0.5_DP,               &  !  22  MAX_ERRORBOUND_P
     &                                            1.E-3_DP,             &  !  23  NONLIN_TOLERANCE
     &                                            0._DP,                &  !  24  USERPARAM4
     &                                            0._DP,                &  !  25  USERPARAM5
     &                                            0._DP,                &  !  26  USERPARAM6
     &                                            0._DP,                &  !  27  USERPARAM7
     &                                            0._DP,                &  !  28  USERPARAM8
     &                                            0._DP,                &  !  29  USERPARAM9
     &                                            0._DP,                &  !  30  USERPARAM10
     &                                            -9999._DP             &  !  31  WAVELENGTH
     &                                               /)
!
!  Names of the parameters,                         alphabetic odering
      character (len=24), parameter :: name_c(n_c)=(/                   &
     &                                               'ADAPTION_TYPE',   &
     &                                               'CSRFORMAT',       &
     &                                               'ELEMENT_TYPE',    &
     &                                               'ERROR_ESTIMATOR', &
     &                                               'HP_ALGORITHM',    &
     &                                               'LINSOLVERTYPE',   &
     &                                               'MATRIXOUTPUT',    &
     &                                               'MATRIXTYPE',      &
     &                                               'MESHFILE',        &
     &                                               'MESHSMOOTHING',   &
     &                                               'MP_COUPLING',     &
     &                                               'MULTIPHYSICS',    &
     &                                               'NONLINEAR',       &
     &                                               'NONLINSOLVERTYPE',&
     &                                               'PHYSICS_MODE',    &
     &                                               'PROJECTPATH',     &
     &                                               'SMOOTHING',       &
     &                                               'TRANSIENTSOLVER', &
     &                                               'WHATSYSTEM',      &
     &                                               'WRITE ITERATION', &
     &                                               'WRITEFIGURE',     &
     &                                               'WRITE_SIMSTATS',  &
     &                                               'WRITEMATRIX',     &
     &                                               'MESHQUALITYPLOT', &
     &                                               'MARKING_TYPE_H',  &
     &                                               'MARKING_TYPE_P',  &
     &                                               'H_ADAPT_CRIT',    &
     &                                               'MATER_VARIATION'  &
     &                                               /)
      character (len=24), parameter :: name_i(n_i)=(/                   &
     &                                               'ADAPT_STEPS',     &
     &                                               'MAXEXETIME',      &
     &                                               'MAXDOF',          &
     &                                               'MAXNODES',        &
     &                                               'NONLINITERSTEPS', &
     &                                               'POLYORDER',       &
     &                                               'DELTA_POLYDEG',   &
     &                                               'HPF8R_ORDER',     &
     &                                               'MAX_ANDERSON_ITER'&
     &                                               /)
      character (len=24), parameter :: name_f(n_f)=(/                   &
     &                                             'ADAPTION_ERROR',    &  !   1
     &                                             'ALPHA_NEWMARK',     &  !   2
     &                                             'DELTA_NEWMARK',     &  !   3
     &                                             'END_TIME',          &  !   4
     &                                             'EPS',               &  !   5
     &                                             'FPRELAX_PARAMETER', &  !   6
     &                                             'FREQUENCY',         &  !   7
     &                                             'GEOMETRY_FACTOR',   &  !   8
     &                                             'INIT_TIME_STEP',    &  !   9
     &                                             'LINSOLVER_ERROR',   &  !  10
     &                                             'MAXMESHSIZE',       &  !  11
     &                                             'NONLIN_ERROR',      &  !  12
     &                                             'SKINNYLIMIT',       &  !  13
     &                                             'START_TIME',        &  !  14
     &                                             'TIME_STEP',         &  !  15
     &                                             'TRAN_GTOL',         &  !  16
     &                                             'TRAN_LTOL_LO',      &  !  17
     &                                             'TRAN_LTOL_UP',      &  !  18
     &                                             'TRAN_REL_TOLER',    &  !  19
     &                                             'ADAPT_FRACTION',    &  !  20
     &                                             'MAX_ERRORBOUND_H',  &  !  21
     &                                             'MAX_ERRORBOUND_P',  &  !  22
     &                                             'NONLIN_TOLERANCE',  &  !  23
     &                                             'USERPARAM4',        &  !  24
     &                                             'USERPARAM5',        &  !  25
     &                                             'USERPARAM6',        &  !  26
     &                                             'USERPARAM7',        &  !  27
     &                                             'USERPARAM8',        &  !  28
     &                                             'USERPARAM9',        &  !  29
     &                                             'USERPARAM10',       &  !  30
     &                                             'WAVELENGTH'         &  !  31
     &                                               /)
!
!  readstd will be set to true, when the file was read
      logical :: readstd=.false.
!
end module femsettings


      recursive subroutine getsetting_c(key,value)
!
!  Return the string value which corresponds to the parameter key
!  in the FEMsetting.txt file
!
!  Specific routine to the generic routine getsetting
!  for the case of a string argument.
!
!    key     one of the Parameter-keywords
!    value   string corresponding to the Parameter setting
!  
      use feminterface, only: readsetting
      use femsettings,  only: n_c, readstd, name_c, param_c
      use femtypes,     only: I4B
      implicit none
      character (len=*) :: key
      character (len=*) :: value
      intent (in) :: key
      intent (out) :: value
!  local variables
      integer (I4B) :: i
!
      if (.not. readstd) call readsetting
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
      end subroutine getsetting_c
!
!
!
      recursive subroutine getsetting_i(key,value)
!
!  Return the integer value which corresponds to the parameter key
!  in the FEMsetting.txt file
!
!  Specific routine to the genereic routine getsetting
!  for the case of a integer argument.
!
!    key     one of the Parameter-keywords
!    value   integer corresponding to the Parameter setting
!  
      use feminterface, only: readsetting
      use femsettings 
      use femtypes
      implicit none
      character (len=*) :: key
      integer (I4B) :: value
      intent (in) :: key
      intent (out) :: value
!  local variables
      integer (I4B) :: i
!
      if (.not. readstd) call readsetting
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
      end subroutine getsetting_i
!
!
!
      recursive subroutine getsetting_f(key,value)
!
!  Return the float value which corresponds to the parameter key
!  in the FEMsetting.txt file
!
!  Specific routine to the generic routine getsetting
!  for the case of a float argument.
!
!    key     one of the Parameter-keywords
!    value   float corresponding to the Parameter setting
!  
      use feminterface, only: readsetting
      use femsettings 
      use femtypes
      implicit none
      character (len=*) :: key
      real (DP) :: value
      intent (in) :: key
      intent (out) :: value
!  local variables
      integer (I4B) :: i
!
      if (.not. readstd) call readsetting
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
      end subroutine getsetting_f
!
!
!
      subroutine putsetting_c(key,value)
!
!  Set the string value which corresponds to the parameter key
!  in the FEMsetting.txt file
!
!  Specific routine to the generic routine putsetting
!  for the case of a string argument.
!
!    key     one of the Parameter-keywords
!    value   string corresponding to the Parameter setting
!  
      use feminterface, only: readsetting
      use femsettings
      use femtypes
      implicit none
      character (len=*) :: key
      character (len=*) :: value
      intent (in) :: key
      intent (in) :: value
!  local variables
      integer (I4B) :: i
!
      if (.not. readstd) call readsetting
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
      end subroutine putsetting_c
!
!
!
      subroutine putsetting_i(key,value)
!
!  Set the integer value which corresponds to the parameter key
!  in the FEMsetting.txt file
!
!  Specific routine to the genereic routine putsetting
!  for the case of a integer argument.
!
!    key     one of the Parameter-keywords
!    value   integer corresponding to the Parameter setting
!  
      use feminterface, only: readsetting
      use femsettings 
      use femtypes
      implicit none
      character (len=*) :: key
      integer (I4B) :: value
      intent (in) :: key
      intent (in) :: value
!  local variables
      integer (I4B) :: i
!
      if (.not. readstd) call readsetting
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
      end subroutine putsetting_i
!
!
!
      subroutine putsetting_f(key,value)
!
!  Set the float value which corresponds to the parameter key
!  in the FEMsetting.txt file
!
!  Specific routine to the genereic routine putsetting
!  for the case of a float argument.
!
!    key     one of the Parameter-keywords
!    value   float corresponding to the Parameter setting
!  
      use feminterface, only: readsetting
      use femsettings 
      use femtypes
      implicit none
      character (len=*) :: key
      real (DP) :: value
      intent (in) :: key
      intent (in) :: value
!  local variables
      integer (I4B) :: i
!
      if (.not. readstd) call readsetting
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
      end subroutine putsetting_f
!
!
!
  subroutine readsetting
!
!  Read settings from the Projectfile.json
!  Alternatively start reading Femsettings.txt
!
      use feminterface, only: low2hi, string2number, getsetting, getsettingfile
      use femsettings, only: param_c, param_i, param_f, n_c, n_f, n_i, readstd,&
                           & name_c, name_i, name_f
!      use json_module  >>>> Disable for now  15.11.15
      use femtypes

      implicit none
!  local variables
      integer (I4B)                  :: i, intgr
      real    (DP)                   :: rl
      character (len=200)            :: path
      character (len=:), allocatable :: projectname, settingsfile, chrctr, error_msg
!      type(json_file)    :: json
      logical            :: found, status_ok, json_present
      readstd=.true.
      json_present = .false.  ! was .true.
!__
! 1) Get Projectpath, Load Project File
    call getsettingfile(path)
!!$!- Get Project Path
!!$    call json%load_file(filename = path(1:len_trim(path)))
!!$    if (json_failed()) then
!!$      call json_check_for_errors(status_ok, error_msg)
!!$      print*,'***ERROR reading json:', error_msg
!!$      json_present = .false.
!!$    else
!!$      call json%get('projectName',projectname,found)
!!$        if (.not.found) print*,'***ERROR reading json: /"projectName/" not found', error_msg
!!$      call json%get('projectPathInfo.'//projectname,settingsfile,found)
!!$        if (.not.found) print*,'***ERROR reading json: /"projectPathInfo.projectFile/" not found', error_msg
!!$    end if
!!$    call json%destroy()
!!$!- Load Project File
!!$    call json%load_file(filename = settingsfile )
!!$    if (json_failed()) then
!!$      call json_check_for_errors(status_ok, error_msg)
!!$      print*,'***ERROR reading json:', error_msg
!!$      json_present = .false.
!!$    else
!!$!__
!!$! 2) Load Character-type Settings
!!$      
!!$      do i = 1, size(name_c)
!!$        call json%get('technicalInfo.FEMSettings.'//name_c(i), chrctr, found)
!!$        if (found) param_c(i) = chrctr(1:len_trim(chrctr))
!!$      end do
!!$!__
!!$! 3) Load Integer-type Settings
!!$
!!$      do i = 1, size(name_i)
!!$        call json%get('technicalInfo.FEMSettings.'//name_i(i), param_i(i), found)
!!$      end do
!!$!__
!!$! 4) Load Float-type Settings
!!$
!!$      do i = 1, size(name_f)
!!$        call json%get('technicalInfo.FEMSettings.'//name_f(i), param_f(i), found)
!!$      end do
!!$!-
!!$    end if
!!$    call json%destroy()
!__
! 5) If no json file was found, read FEMSettings.txt
    if (json_present .eqv. .false.) then
      call read_femsettingsfile
    end if
    
    return
  end subroutine readsetting
!
!
      subroutine read_femsettingsfile
!
!  Read the FEMsetting.txt file of parameters
!
      use feminterface, only: low2hi, string2number, getsetting, getsettingfile
      use femsettings, only: param_c, param_i, param_f, n_c, n_f, n_i, readstd,&
                           & name_c, name_i, name_f
      use femtypes

      implicit none
!  local variables
      integer (I4B) :: i, unitid, pos, epos, ios, rpos, length
      logical :: exist, search
      character (len=200) :: str, name, path
      character (len=200) :: settingfile
!
      readstd=.true.

!  look for file
      call getsettingfile(settingfile)
      settingfile = settingfile(1:len_trim(settingfile))


      call grglun(unitid)
      inquire (file=settingfile,exist=exist)
      if (.not. exist) then
        return
      end if
      open (unitid,file=settingfile,form='FORMATTED',iostat=ios)
      if (ios .ne. 0) then
        print*,'Error in opening FEMsettings.txt'
      end if
!
!  Line by line loop over all entries in file
      do
        read(unitid,fmt='(a)',iostat=ios) str
        if (ios .lt. 0) then
          print*
          close(unitid)
          exit
        else if (ios .gt. 0) then 
          print*,'Error in reading FEMsettings.txt :',str
        end if
!
!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  allow for blank or comment lines
        length=len_trim(str)
        if (length .eq.0) cycle
!
        print*,str(1:length)
        str=adjustl(str)
        epos=index(str,'=')
        call low2hi(str,epos)
!  Except for the path, we do not distinguish between lower and upper case
        if (str(1:epos-1).ne.'PROJECTPATH') call low2hi(str,length)
!
!  Loop over all parameters of string type
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
          print*,'Unknown parameter found'
          print*,'Error in reading FEMsettings.txt :',str
          pause
        end if
      end do
!
!  read once more the femsettingfile, now in the project directory
      call grglun(unitid)
      call getsetting('PROJECTPATH',path)
      inquire (file=path(1:len_trim(path))//'FEMsettings.txt',exist=exist)
      if (.not. exist) then
        return
      end if
      open (unitid,file=path(1:len_trim(path))//'FEMsettings.txt',      &
     &  form='FORMATTED',iostat=ios)
      if (ios .ne. 0) then
        print*,'Error in opening FEMsettings.txt'
      end if
!
!  Line by line loop over all entries in file
      do
        read(unitid,fmt='(a)',iostat=ios) str
        if (ios .lt. 0) then
          print*
          close(unitid)
          exit
        else if (ios .gt. 0) then 
          print*,'Error in reading FEMsettings.txt :',str
        end if
!
!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  allow for blank or comment lines
        length=len_trim(str)
        if (length .eq.0) cycle
!
        print*,str(1:length)
        str=adjustl(str)
        epos=index(str,'=')
        call low2hi(str,epos)
!  except for the path, we do not distinguish between lower and upper case
        if (str(1:epos-1).ne.'PROJECTPATH') call low2hi(str,length)
!
!  Loop over all parameters of string type
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
          print*,'Unknown parameter found'
          print*,'Error in reading FEMsettings.txt :',str
          pause
        end if
      end do
      return
      end subroutine read_femsettingsfile

!
      subroutine wpathsetting
!
!  Replace the project path in the FEMsetting.txt file with the actual one
!
      use feminterface, only: low2hi, getsettingfile
      use femsettings 
      use femtypes
      implicit none
!  local variables
      integer (I4B) :: i, pathpos
      integer (I4B) :: unitid, pos, epos, ios, unitidr
      character (len=250) :: str
      character (len=200) :: settingfile
      logical :: exist
!
!  look for file
      call getsettingfile(settingfile)
      settingfile = settingfile(1:len_trim(settingfile))

      call grglun(unitid)
      inquire (file=settingfile,exist=exist)
      if (.not. exist) then
        call writsetting
      end if
      open (unitid,file=settingfile,form='FORMATTED',iostat=ios)
      if (ios .ne. 0) then
        print*,'Error in opening FEMsettings.txt'
      end if
      call grglun(unitidr)
      open (unitidr,file=settingfile//'1',form='FORMATTED',iostat=ios)
      if (ios .ne. 0) then
        print*,'Error opening tmp file'
      end if
!  copy the settings file 
      do
        read(unitid,fmt='(a)',iostat=ios) str
        if (ios .lt. 0) then
          close(unitid,status='DELETE')
          close(unitidr)
          exit
        else if (ios .gt. 0) then 
          print*,'Error in reading FEMsettings.txt :',str
        end if
        write(unitidr,fmt='(a)') str(1:len_trim(str))        
      end do
      call grglun(unitid)
      open (unitid,file=settingfile//'1',form='FORMATTED',iostat=ios)
      call grglun(unitidr)
      open (unitidr,file=settingfile,form='FORMATTED',iostat=ios)
!
!  Line by line loop over all entries in file
      do
        read(unitid,fmt='(a)',iostat=ios) str
        if (ios .lt. 0) then
!  return if the attempt to read a line failed
          close(unitid,status='DELETE')
          close(unitidr)
          return
        else if (ios .gt. 0) then 
          print*,'Error in reading FEMsettings.txt :',str
        end if
        write(unitidr,fmt='(a)') str(1:len_trim(str))        
!
!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  allow for blank or comment lines
        if (len_trim(str) .eq.0) cycle
!
!  Look for the projectpath and replace it with the actual one
        epos=index(str,'=')
        call low2hi(str,epos)
        pos=index(str,'PROJECTPATH')
        if (pos .gt. 0 .and. epos .gt. pos) then
          backspace(unitidr)
!  Look for position of projectpath in name list and store position.
!  Used for sorting purposes only.
          do i = 1,n_c
            if (name_c(i) .eq. 'PROJECTPATH') then
              pathpos = i
              exit
            end if
          end do
          str=adjustl(param_c(pathpos))
          write(unitidr,fmt='(a)') 'PROJECTPATH = '//str(1:len_trim(str))
        end if
      end do
!
      return
      end subroutine wpathsetting
!
!
!
      subroutine writsetting
!
!  Write the FEMsetting.txt file of parameters
      use feminterface, only: getsettingfile
      use femsettings 
      use femtypes
      implicit none
!  local variables
      integer (I4B) :: i, unitid
      character (len=200) :: str,strp
      character (len=200) :: settingfile      
      logical :: opened
!
      call getsettingfile(settingfile)
      settingfile = settingfile(1:len_trim(settingfile))
      inquire (file=settingfile,opened=opened)
      if (opened) then
        inquire (file=settingfile,number=unitid)
        close(unitid)
      end if
      call grglun(unitid)
      open (unitid,file=settingfile,form='FORMATTED')
!
!  Write parameters of string type
      do i=1,n_c
        str=name_c(i)
        str=adjustl(str)
        strp=adjustl(param_c(i))
        write(unitid,*) str(1:len_trim(str)),' = ',strp(1:len_trim(strp))
      end do 
!
!  Write parameters of integer type
      do i=1,n_i
        str=name_i(i)
        str=adjustl(str)
        write(unitid,*) str(1:len_trim(str)),' = ',param_i(i)
      end do
!
!  Write parameters of float type 
      do i=1,n_f
        str=name_f(i)
        str=adjustl(str)
        write(unitid,*) str(1:len_trim(str)),' = ',param_f(i)
      end do 

      close(unitid)
      
      return
      end subroutine writsetting
