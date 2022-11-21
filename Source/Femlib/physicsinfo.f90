      subroutine physicsinfo(physics,nnat,dimensions)
      use feminterface, only: readphysics
      use femtypes
      implicit none
      integer (I4B) :: nnat, dimensions
      character (len=*) :: physics
      intent (in) :: physics, dimensions
      intent (out) :: nnat
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
!    $Revision: 1.11 $
!    $Date: 2014/09/11 14:10:57 $
!    $Author: m_kasper $
!
!
!  the routine is intended to return properties related to a physics mode
!  actually it only returns the number of natures nnat for both 2d and 3d cases
!  in future it may be used to return other properties of the physics mode
!  e.g. symmetric(equation system), real/complex, ...
!
!  local variables
      logical :: ok

      call readphysics(physics, dimensions, nnat, ok)

      if (.not. ok) then
        stop
      end if

      return
      end subroutine physicsinfo


      subroutine readphysics(physics, dim, nnat, ok)
      use femtypes
      use feminterface, only: strtok, low2hi, string2number, tab2blank
      use matconstants, only: numparam, parameternames, parameterunit, parameterdefault
      use physics_quantities
      implicit none
      integer (I4B) :: dim, nnat
      character (len=*) :: physics
      logical :: ok
      intent (in) :: physics, dim
      intent (out) :: nnat, ok
!
!  readphysics - read the physics file and store in module matconstants
!
!  input:
!  physics          the actual physics mode
!  dim              no of spatial dimensions, i.e. '2' or '3'
!
!  output:
!  nnat             number of natures
!  ok               =.true. if completed successfully
!
!

!  local variables
      integer (I4B) :: unitid, ios, linenumber, i, inat, k
      real (DP) :: value
      character (len=200) ::  polydepath, str
      character (len=50) ::  token, token2
      logical :: exist, found, fine

!
      dimensions = dim
      ok = .true.
      call getenv("PolydePath",polydepath)

!  look for file
      inquire (file=polydepath(1:len_trim(polydepath))//'Physics/PhysicsModes.txt', &
     &         exist=exist,iostat=ios)
      if (ios .ne. 0 .or. exist.eq..false.) then
        ok = .false.
        print*,'The file ''PhysicsModes.txt'' in the directory: ',polydepath(1:len_trim(polydepath))//'Physics',' is missing'
        return
      end if
      call grglun(unitid)
      open (unitid,file=polydepath(1:len_trim(polydepath))//'Physics/PhysicsModes.txt', &
     &      form='FORMATTED',position='REWIND',action='READ',iostat=ios)
      if (ios .ne. 0) then
        ok = .false.
        print*,'Error in opening PhysicsModes.txt'
      end if

!  line by line read PhysicsModes.txt and find "Structure = dimensions"
      linenumber = 0
      found = .false.
      do
        read(unitid,fmt='(a)',iostat=ios) str
        linenumber = linenumber + 1
        if (ios .lt. 0) then
          exit
        else if (ios .gt. 0) then
          print*,'**** Error at line:',linenumber,' while reading PhysicsModes.txt'
          ok = .false.
          return
        end if
!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  allow for blank or comment lines
        call tab2blank(str)
        if (len_trim(str) .eq. 0) cycle
!  get first token
        call strtok(str,     '=;:', token)     ! char(9)= horizontal tab; char(13)= CR
        call low2hi(token,len_trim(token))
        token = adjustl(trim(token))
        if (token .ne. 'STRUCTURE') cycle
!  get second token
        call strtok(char(0), '=,;:', token)
        call low2hi(token,len_trim(token))
        token = adjustl(trim(token))
        if ((dimensions .eq. 2 .and. token .eq. '2D') .or.              &
            (dimensions .eq. 3 .and. token .eq. '3D')) then
          found = .true.
          exit
        end if
      end do ! read line for "Structure"

      if (.not.found) then
        ok =.false.
        print*,'Error in file PhysicsModes.txt, section = 2D or 3D not found'
        return
      end if

! line by line read PhysicsModes.txt and find "PHYSICSMODE = physics"
      found = .false.
      do
        read(unitid,fmt='(a)',iostat=ios) str
        linenumber = linenumber + 1
        if (ios .lt. 0) then
          exit
        else if (ios .gt. 0) then
          print*,'**** Error at line:',linenumber,' while reading PhysicsModes.txt'
          ok = .false.
          return
        end if
!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  allow for blank or comment lines
        call tab2blank(str)
        if (len_trim(str) .eq. 0) cycle
!  get first token
        call strtok(str,     '=;:', token)     ! char(9)= horizontal tab; char(13)= CR
        call low2hi(token,len_trim(token))
        token = adjustl(trim(token))
        if (token .ne. 'PHYSICSMODE') cycle
!  get second token
        call strtok(char(0), '=,;:', token)
        call low2hi(token,len_trim(token))
        token = adjustl(trim(token))
        if (token .eq. physics) then
          found = .true.
          exit
        end if
      end do

      if (.not.found) then
        ok =.false.
        print*,'Error in reading the file PhysicsModes.txt'
        print*,'the PHYSICSMODE: ', physics(1:len_trim(physics)), ' was not found'
        return
      end if


!  read properties of the Physicsmode
      nnat = 0
      numparam = 0
      do
        read(unitid,fmt='(a)',iostat=ios) str
        linenumber = linenumber + 1
        if (ios .lt. 0) then
          exit
        else if (ios .gt. 0) then
          print*,'**** Error at line:',linenumber,' while reading ''PhysicsModes.txt'''
          ok = .false.
          return
        end if
!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
        if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
        if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
        if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  allow for blank or comment lines
        call tab2blank(str)
        if (len_trim(str) .eq. 0) cycle
!  get first token
!print*,linenumber, trim(str)
        call strtok(str,     '=(:', token)     ! char(9)= horizontal tab; char(13)= CR
        call low2hi(token,len_trim(token))
        token = adjustl(trim(token))

        select case(token)

        case ('NNAT')
          call strtok(char(0), '=,;:', token2)
          call string2number(token2,nnat)
          allocate (eltypes(nnat))
          allocate (pot_names(nnat), pot_descriptor(nnat), pot_unit(nnat),   &
     &              gx_names(nnat),  gx_descriptor(nnat),  gx_unit(nnat),    &
     &              gy_names(nnat),  gy_descriptor(nnat),  gy_unit(nnat),    &
     &              gz_names(nnat),  gz_descriptor(nnat),  gz_unit(nnat),    &
     &              g_names(nnat),   g_descriptor(nnat),   g_unit(nnat),     &
     &              fx_names(nnat),  fx_descriptor(nnat),  fx_unit(nnat),    &
     &              fy_names(nnat),  fy_descriptor(nnat),  fy_unit(nnat),    &
     &              fz_names(nnat),  fz_descriptor(nnat),  fz_unit(nnat),    &
     &              f_names(nnat),   f_descriptor(nnat),   f_unit(nnat)  )
          pot_names(nnat) = 'unknown'
          gx_names(nnat)  = 'unknown'
          gy_names(nnat)  = 'unknown'
          gz_names(nnat)  = 'unknown'
          g_names(nnat)   = 'unknown'
          fx_names(nnat)  = 'unknown'
          fy_names(nnat)  = 'unknown'
          fz_names(nnat)  = 'unknown'
          f_names(nnat)   = 'unknown'


        case ('ANALYSIS TYPE')


        case ('NONLINEAR')


        case ('ELTYPE')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
            call strtok(char(0), '=,;:', token2)
            call string2number(token2, k)
            eltypes(i) = k
          end if


        case ('NUMPARAM')
          call strtok(char(0), '=,;:', token2)
          call string2number(token2,numparam)
          allocate (parameternames(numparam), parameterdefault(numparam), parameterunit(numparam))


        case ('PARAMETER')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. numparam) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of parameters'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            parameternames(i) = token2
!     Default Value
            call strtok(char(0), ';,', token2)
            call string2number(token2, value, fine)
            if (.not.fine) then
              print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
              print*,'Conversion error for:',token2
              ok = .false.
            else
              parameterdefault(i) = value
            end if
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            parameterunit(i) = token2
          end if


        case ('FIELDQUANTITY')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            pot_names(i) = token2
!     Descriptor
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            pot_descriptor(i) = token2
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            pot_unit(i) = token2
          end if


        case ('GX')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            gx_names(i) = token2
!     Descriptor
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            gx_descriptor(i) = token2
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            gx_unit(i) = token2
          end if


        case ('GY')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            gy_names(i) = token2
!     Descriptor
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            gy_descriptor(i) = token2
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            gy_unit(i) = token2
          end if


        case ('GZ')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            gz_names(i) = token2
!     Descriptor
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            gz_descriptor(i) = token2
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            gz_unit(i) = token2
          end if


        case ('G')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            g_names(i) = token2
!     Descriptor
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            g_descriptor(i) = token2
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            g_unit(i) = token2
          end if


        case ('FX')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            fx_names(i) = token2
!     Descriptor
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            fx_descriptor(i) = token2
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            fx_unit(i) = token2
          end if


        case ('FY')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            fy_names(i) = token2
!     Descriptor
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            fy_descriptor(i) = token2
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            fy_unit(i) = token2
          end if


        case ('FZ')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            fz_names(i) = token2
!     Descriptor
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            fz_descriptor(i) = token2
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            fz_unit(i) = token2
          end if


        case ('F')
          call strtok(char(0), '=,;:', token2)
          token2 = token2(index(token2,'(')+1:index(token2,')')-1)
          call string2number(token2, i)
          if (i .gt. nnat) then 
            print*,'Error while reading ''PhysicsModes.txt'' at line',linenumber
            print*,'Index must not be larger than the number of natures'
            ok = .false.
          else
!     Name
            call strtok(char(0), '=;,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            f_names(i) = token2
!     Descriptor
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            f_descriptor(i) = token2
!     Unit
            call strtok(char(0), ';,', token2)
            call low2hi(token2,len_trim(token2))
            token2 = adjustl(trim(token2))
            f_unit(i) = token2
          end if

        case ('PHYSICSMODE','STRUCTURE') 
!  a further physics mode starts
          exit


        case default
          print*,'Error in reading the file PhysicsModes.txt'
          print*,'for the PHYSICSMODE: ', physics(1:len_trim(physics)), ' at line', linenumber
          print*,'the parameter: ''', token(1:len_trim(token)), ''' is not allowed'
        end select
      end do
      close (unitid)

      return
      end subroutine
