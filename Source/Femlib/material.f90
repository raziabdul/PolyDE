      subroutine readmat(matname,found,matindex)
      use femtypes
      use feminterface,       only: getsetting, low2hi, string2number
      use feminterface,       only: setstandardvalues, putmatparam, usermaterials
      use feminterface3D,     only: usermaterials3D
      use matconstants,       only: ipos, matnames, maxmat, numparam, parameternames, usermat
      use varmat_handler,     only: readmatprops
      use physics_quantities, only: dimensions
      implicit none
      character (len=*) matname
      integer (I4B) matindex
      logical found
      intent (out) :: found, matindex
      intent (in) :: matname
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
!    $Revision: 1.20 $
!    $Date: 2014/07/03 12:50:52 $
!    $Author: m_kasper $
!
!  Read the material file and store it internally
!
!  Input
!    matname      the name of the material (string)
!  Output
!    found        is .true. if the material file was found on the disc
!    matindex     the material number under which the material properties are accessible
!
!  local variables
!
      integer (I4B)       :: i, unitid, pos, epos, ios, iparam, linenumber, rpos
      character (len=100) :: str
      character (len=200) :: path, name
      character (len=24)  :: matvar
      character (len=50)  :: matvar_file, matvar_format
      logical             :: ok, foundparam, ex, firstcall=.true.
!      
      external :: grglun
!__
! Checks and Preparations:

       call getsetting('PROJECTPATH',path)
       call getsetting('MATER_VARIATION',matvar)
       call getsetting('AMV_FILENAME',matvar_file)
       call getsetting('AMV_FILEFORMAT',matvar_format)
!-
      if (firstcall) then 
!  first time that we want to store a material 
!  allocate arrays, initialize names of parameters and set standard values
        call setstandardvalues
        firstcall=.false.
        ipos=0
      end if
!
! special case:  air
      if (matname.eq.'0' .or. matname.eq.'nothing') then
        matindex=0
        return
      end if
!
!  Check whether material was already stored previously
!
      do i = 1,ipos
        if (matname.eq.matnames(i)) then
!  Material was found in the internally stored list
          found=.true.
          matindex=i
          return
        end if
      end do
!
!  Material was not found in internal list, now read from disc
!
     inquire(file=path(1:len_trim(path))//matname(1:len_trim(matname))//'.material',EXIST=ex)
      if (ex) then
        found=.true.
!  the file exists
        ipos=ipos+1
        if (ipos.gt.maxmat) then
          write(*,1000) maxmat
1000      format('**** actual number of materials exceeds the maximum number of ',i2,' matarials')
          stop
! TO DO #######  need to implement reallocation of arrays
        else
!  index at which the material data are stored
          matindex=ipos
          usermat(ipos)=.false.
          matnames(ipos)=matname
        end if
!  now read the file
        call grglun(unitid)
        open (unitid,file=path(1:len_trim(path))//                      &
     &    matname(1:len_trim(matname))//'.material',                    &
     &    form='FORMATTED', status='OLD', iostat=ios, position='REWIND')
        if (ios .ne. 0) then
          print*,'Error reading material file: ',matname(1:len_trim(matname))//'.material'
        end if
!
!  
!  Line by line loop over all entries in file
        linenumber=0
        do
          read(unitid,fmt='(a)',iostat=ios) str
          if (ios .lt. 0) then
!  reached end of file
            close(unitid)
            return
          else if (ios .gt. 0) then 
            print*,'Error reading material file: ',matname(1:len_trim(matname))//'.material'
            print*,' in line',linenumber,' : ','"'//str(1:77)//'"'
          end if
          linenumber=linenumber+1
!
!  ignore comment lines
!  at each line, everything following either ! or # or/* is ignored
          if (index(str,'!') .ne. 0) str=str(1:index(str,'!')-1)
          if (index(str,'#') .ne. 0) str=str(1:index(str,'#')-1)
          if (index(str,'/*') .ne. 0) str=str(1:index(str,'/*')-1)
!  allow for blank or comment lines
          if (len_trim(str) .eq.0) cycle
!
!  allow leading blanks
          str=adjustl(str)
          epos=index(str,'=')
          if (epos .eq. 0 ) then
!  we were expecting parameter = ... , but did not find the assignment sign
            print*,'Error reading material file: ',matname(1:len_trim(matname))//'.material'
            print*,' in line',linenumber,' : ','"'//str(1:77)//'"'
          else            
!  we do not distinguish between lower and upper case
            call low2hi(str,epos)
            foundparam=.false.
!  Loop over the list of physical parameters which are needed for this problem type
            do i=1,numparam
              name=parameternames(i)
              pos=index(str(1:epos),name(1:len_trim(name)))
              rpos=index(name(1:len_trim(name)),str(1:len_trim(str(1:epos-1))))
              if (pos .gt. 0 .and. rpos .gt. 0) then
                if (epos .gt. pos) then
                  str=str(epos+1:)
                  str=adjustl(str)
                  foundparam=.true.
                  iparam=i
                  exit
                end if
              end if
            end do
          end if          
          if (foundparam) then 
!  interpret the assignment
            call putmatparam(iparam,str)
!          else 
!            print*,'Warning, required material parameter: ',name(1:len_trim(name)),&
!     &             ' not found in material file: ',matname(1:len_trim(matname))//'.material'
          end if          
        end do
!
        close(unitid)
      else 
!
!  Material was neither found in internal list nor on disc, try the userinterface
!
        found=.false.
        if (matvar.eq.'YES') then
            if (matvar_file.ne.'NONAME') then
                call readmatprops(matvar_file(1:len_trim(matvar_file)),matvar_format(1:len_trim(matvar_format)), ok)
            end if
        end if
        
!
        if (dimensions .eq. 2) then
          call usermaterials(matname,found)
        else
          call usermaterials3D(matname,found)
        end if
        ipos=ipos+1
        if (ipos.gt.maxmat) then
          write(*,1000) maxmat
          stop
! TO DO #######  need to implement reallocation of arrays
        else
!  index at which the material data are stored
          usermat(ipos)=.true.
          matindex=ipos
          matnames(ipos)=matname
!  copy standartvalues??
        end if

        if (.not. found) then 
          write (*,2000) matname(1:len_trim(matname)),path(1:len_trim(path))
2000      format(' ***** Error, Material: ',A,'.material',        &
     &       ' the Material-file does not exist in directory ',A)
          stop
        end if
      end if
      
      return
      end subroutine readmat
!
!
!
      subroutine putmatparam(iparam,str)
      use femtypes
      use feminterface, only: strtok, string2number
      use matconstants, only: ipos, param
      implicit none
      integer(I4B) iparam
      character (len=*) str
      intent (in) :: iparam,str
!  interpret the string str and assign the parameter to internal arrays
!  Input
!    iparam       actual parameter index
!    str          the string to interpret
!
!  local variables
      real (DP) value
      character (len=3) delim
      character (len=100) token
!
!  delimiter
      delim=' ,;'
!  read a scalar
      call strtok(str, delim, token)
      call string2number(token,value)
      param(ipos)%d(iparam)=value
      return
      end subroutine putmatparam
!
!
!
      pure subroutine fetchmatparameters(list,matindex,xyzs,elem)
      use femtypes
      use feminterface, only: usermaterials, low2hi
      use feminterface3D, only: usermaterials3D
      use matconstants, only: numparam, usermat, parameternames, matnames, param
      use physics_quantities, only: dimensions
      implicit none
      integer (I4B) matindex
      integer (I4B), optional :: elem
      real(DP) list(:), xyzs(:)
      intent (in) :: matindex, xyzs, elem
      intent (out) :: list
!  retrieve the material parameters of this material
!  Input
!    matindex     actual parameter index
!           = 0   the standard values (material = 0 or nothing)
!          <> 0   a material which had been read form a material file, or a user material
!    xyzs         (world-)coordinates for which the material values are to be returned
!    elem         element number (optional, if use in usermaterial)
!  Output
!    list         list of material parameters as defined in the subroutine setstandardvalues
!
!  local variables
      integer (I4B) i, j, numnames
      real(DP) values(numparam)
      character(len=20) names(numparam)
      logical found
!
      list(1:numparam)=param(matindex)%d(1:numparam)
!  in the case of usermaterials 
      if (matindex .ne. 0) then
        if (usermat(matindex)) then
          if (dimensions .eq. 2) then
            call usermaterials(matnames(matindex),found,xyzs,names,values,numnames,elem)
          else
            call usermaterials3D(matnames(matindex),found,xyzs,names,values,numnames,elem)
          end if
!  translate names to list entries
          do i=1,numnames
            call low2hi(names(i),100)
            do j=1,numparam
              if (names(i) .eq. parameternames(j)) then 
                list(j)=values(i)
                exit     
              end if
            end do
          end do
        end if
      end if
      return
      end subroutine fetchmatparameters
