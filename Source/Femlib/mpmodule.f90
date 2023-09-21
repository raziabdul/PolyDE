      module mpmodule
      use femtypes
      implicit none
      public
!
!-------------------------------------------------------------------------------
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
!    $Revision: 1.37 $
!    $Date: 2010/10/25 11:56:58 $
!    $Author: chvokas $
!
!------------------------------------------------------------------------------
!  This module contains subroutines used for multiphysics in PolyDE. All sub-
!  routines only use real variables atm.
!------------------------------------------------------------------------------
!
!  MULTIPHYSICS VARIABLES
!
!  Data storage:
!  -------------
!  It is assumed, that the datatype of individual directions for T(emperature),
!  E(lectric) or M(echanical) field are the same. If data is given as a list,
!  last dimension of z is 1 (zM(:,:,1) or zT(:,1)). Otherwise, arrays for
!  electric and mechanical data are three-dimensional.
      character (len=10), save :: datatypeE, datatypeM, datatypeP, datatypeT
      real (DP), allocatable, save :: xE(:), yE(:), zE(:,:,:)
      real (DP), allocatable, save :: xM(:), yM(:), zM(:,:,:)
      real (DP), allocatable, save :: xP(:), yP(:), zP(:,:)
      real (DP), allocatable, save :: xT(:), yT(:), zT(:,:)
!
!  Transformation variables:
!  -------------------------
      real (DP), save :: trmat(3,3)
!
!  Multiphysics switch and kind of coupling modification:
!  ------------------------------------------------------
      logical, save :: multiphysics, eomod, momod, tomod, etmod



      contains
!
!
!
      subroutine bilintpol(xval,yval,zval,xp,yp,zp,errcode)
      use feminterface, only: hunt
      use femtypes
      implicit none
      integer (I4B) :: errcode
      real (DP) :: xp, yp, zp
      real (DP) :: xval(:), yval(:), zval(:,:)
      intent (in)  :: xval, yval, zval, xp, yp
      intent (out) :: zp
!
!
!  This subroutine performs a bininear interpolation to find the value for a
!  given point (xp,yp) on a regular grid given by coordinate vectors xval and
!  yval. Values are contained in array zval. The interpolated value is returned
!  as zp.
!
!  The algorithm first searches for lower and upper coordinates to determine
!  grid cell for point (xp,yp). Then it assigns local coordinates t and u, both
!  ranging from 0 to 1 and computes the field value zp by bilinear interpolation.
!
!  (xl,yu) *-------------* (xu,yu)   ^
!          |             |           |
!          |  *(xp,yp)   |           |
!          |             |           u
!          |             |           |
!  (xl,yl) *-------------* (xu,yl)   -
!
!          |----- t ----->
!
!
!  Formula is given in:
!  Numerical Recipies in FORTRAN 77: The Art of Scientific Computing,
!  Chapter 3.6, www.nr.com
!
!
!  Input:
!     xval, yval         coordinates of grid points
!     zval               field values of grid points
!     xp, yp             coordinates of interpolation point

!  Output:
!     errcode            if errcode 10 --> (xp,yp) not in grid
!     zp                 field value of interpolation point

!  Internal variables:
      integer (I4B) :: ixl, ixu, iyl, iyu
      integer (I4B), save :: xlo, ylo
      real (DP) :: t, u
!
!
!  cycle over grid points to find grid cell for (xp,yp)
!  find lower x-coord entry
      call hunt(xval,xp,xlo)
      if (((xlo .eq. 0) .or. (xlo .eq. size(xval))) &
          .and. (xp .ne. xval(size(xval)))) then
        errcode = 10
      else
        ixl = xlo
        errcode = 0
      end if
!  find lower y-coord entry
      if (errcode .eq. 0) then
        call hunt(yval,yp,ylo)
        if (((ylo .eq. 0) .or. (ylo .eq. size(yval))) &
            .and. (yp .ne. yval(size(yval)))) then
          errcode = 10
        else
          iyl = ylo
          errcode = 0
        end if
      end if
!
      if (errcode .eq. 0) then
!  assign upper entry and compute local variable (ranging from 0 to 1)
        if (xp .eq. xval(size(xval))) then
          ixu = ixl
          t = 0
        else
          ixu = ixl + 1
          t = (xp - xval(ixl)) / (xval(ixu) - xval(ixl))
        end if
!
        if (yp .eq. yval(size(yval))) then
          iyu = iyl
          u = 0
        else
          iyu = iyl + 1
          u = (yp - yval(iyl)) / (yval(iyu) - yval(iyl))
        end if

!  compute value for zp
        zp = (1._DP - t)*(1._DP - u)*zval(ixl,iyl) + t*(1._DP - u)*zval(ixu,iyl) + &
             t*u*zval(ixu,iyu) + (1._DP - t)*u*zval(ixl,iyu)
      else
!  set zp to 0, if point not in grid
        zp = 0._DP
      end if
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      subroutine exportdata(datatype,fieldtype,phi,xplmin,xplmax,       &
     &                      yplmin,yplmax,gridx,gridy)
      use feminterface, only: elemnt, fieldquantity, get2Dintegpoints,  &
                              getsetting, lam2xy, reallocate, xy2lam
      use femtypes
      use globalvariables
      implicit none
      integer (I4B) :: gridx, gridy
      real (DP) :: phi, xplmin, xplmax, yplmin, yplmax
      character (len=*) :: datatype, fieldtype
      intent (in) :: datatype, fieldtype, phi
      intent (in) :: xplmin, xplmax, yplmin, yplmax, gridx, gridy
!
!
!  Subroutine to export values for given fieldtype and phase as:
!  a) scattered data given in format xval yval zval (datataype = 'SCATTERED')
!  b) grid data on a regular grid of gridx * gridy points  (datataype = 'GRID')
!
!  to the file data<n>.txt, where <n> is an integer. The file contains the data-
!  type, fieldtype, origin of data and:
!  a) number of grid points in x- and y-direction (gridx and gridy), x-coords
!     (xval) and y-coords (yval) of points and matrix of grid values (zval).
!  b) number of total entries (num) and xval, yval, zval list of values
!
!  Comments in file additionally include accuracy, date and time of solution.
!
!
!  Input:
!     datatype           'SCATTERED' or 'GRID' for corresponding type
!     fieldtype          field quantity to compute (e.g. POTENTIAL, EY ...)
!     phi                phase / time instant to compute
!     xplmin, xplmax     extents in x-direction for grid
!     yplmin, yplmax     extents in y-direction for grid
!     gridx, gridy       number of grid points in x- and y-direction

!  Internal variables:
      integer (I4B) :: elem, i, ierror, j, num, unitid
!      integer (I4B) :: errcode, npkt, plength
      real (DP) :: delta(2), lambda(3)
      real (DP), allocatable :: xval(:), yval(:), zval(:,:)
!      real (DP), pointer :: weig(:), lamint(:,:)
      real (DP), pointer :: xtmp(:), ytmp(:), ztmp(:)
      character (len=8)  :: dat
      character (len=10) :: tim, tmp, tmp1, unit
      character (len=23) :: fmtreal
      character (len=50) :: descriptor
      character (len=200):: path
      complex (DPC) :: fieldval
      logical :: ok, ok1
!      
      external :: grglun
!
      select case (datatype)
      case ('GRID')
!  compute distance between neighbouring grid points
        delta(1) = (xplmax - xplmin)/(gridx-1)
        delta(2) = (yplmax - yplmin)/(gridy-1)
!  allocate coordinate vectors and solution value matrix
        allocate(xval(gridx), yval(gridy), zval(gridx,gridy))
        xval = (/(xplmin + (i-1)*delta(1) ,i=1, gridx)/)
        yval = (/(yplmin + (i-1)*delta(2) ,i=1, gridy)/)
!  loop over all grid points to compute fieldvalue
        do i = 1, size(xval)
          do j = 1, size(yval)
!  get element for point
            call elemnt(xn,yn,e,n,xval(i),yval(j),elem,en,ok)
!  if point is in an element, compute value
            if (ok) then
              call xy2lam(xval(i),yval(j),elem,lambda,xn,yn,e)
              call fieldquantity(elem,fieldtype,lambda,phi,fieldval,    &
     &                           descriptor,unit,ok1)
              if (.not.ok1) then
                print *,"***** fieldtype: ",trim(fieldtype)," unknown for ",trim(physics)
                stop
              end if
              zval(i,j) = real(fieldval,DP)
!  if point is not in an element, set value to 0
            else
              zval(i,j) = 0._DP
            end if
          end do
        end do
!
      case ('SCATTERED')
        allocate(xtmp(p), ytmp(p), ztmp(p))
        xtmp = xn
        ytmp = yn
!  loop over all nodes and compute field values
        do i = 1, p
!  get element for point
          call elemnt(xn,yn,e,n,xtmp(i),ytmp(i),elem,en,ok)
!  if point is in an element, compute value
          if (ok) then
            call xy2lam(xtmp(i),ytmp(i),elem,lambda,xn,yn,e)
            call fieldquantity(elem,fieldtype,lambda,phi,fieldval,      &
     &                         descriptor,unit,ok1)
            if (.not.ok1) then
              print *,"***** fieldtype: ",trim(fieldtype)," unknown for ",trim(physics)
              stop
            end if
            ztmp(i) = real(fieldval,DP)
!  if point is not in an element, this may not happen
          else
            ztmp(i) = 0._DP
            print '(a)',"**** Error in exportdata. Point not in an element. This may not happen!"
          end if
        end do
!! NEXT LOOP CORRECTLY WRITES VALUES, BUT IDTANG IN BIVAR.F90 TELLS, THERE ARE
!! IDENTICAL POINTS DURING INTERPOLATION. SO LOOP IS COMMENTED.
!  loop over all elements to get additional field values at 2d integration points
!        do elem = 1,n
!          call get2Dintegpoints(ep(elem),npkt,weig,lamint,errcode)
!  length/size of original pointer
!          plength = size(xtmp)
!  reallocate arrays to make place for new entries
!          xtmp => reallocate(xtmp,plength+npkt)
!          ytmp => reallocate(ytmp,plength+npkt)
!          ztmp => reallocate(ztmp,plength+npkt)
!  loop over all integration points depending on polynomial degree of element
!          do i = 1, npkt
!  position of first free entry in pointer
!            j = plength + i
!  get xtmp and ytmp coordinate from integ point and compute field value
!            call lam2xy(lamint(:,i),elem,xtmp(j),ytmp(j),xn,yn,e)
!            call fieldquantity(elem,fieldtype,lamint(:,i),phi,fieldval,descriptor,unit)
!            ztmp(j) = real(fieldval,DP)
!          end do
!          deallocate(weig,lamint)
!        end do
!  allocate arrays and store values from temporary pointers to xval, yval and zval
        num = size(xtmp)
        allocate(xval(num), yval(num), zval(num,1))
        xval = xtmp
        yval = ytmp
        zval(:,1) = ztmp
        deallocate(xtmp,ytmp,ztmp)
!
      case default
! print error message if datatype unknown
        print '(a)',"***** Datatype "//trim(datatype)//" unknown."
        print '(a)',"Cannot compute values. Check POSTsettings.txt"
      end select

!  write data to file data<n>.txt
      call grglun(unitid)
      call getsetting('PROJECTPATH',path)
!  open a file in the specified project path
      open (unitid,FILE=path(1:len_trim(path))//'data_'//trim(fieldtype)//'.txt',&
     &      STATUS='NEW',FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
      if (ierror .ne. 0) then
        do i = 1,99
          write (tmp ,'(i2)') i
          open (unitid,FILE=path(1:len_trim(path))//'data_'//trim(fieldtype)//trim(adjustl(tmp))//'.txt',&
     &          STATUS='NEW',FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
          if (ierror .ne. 0) then
            cycle
          else
            print '(a)',"Writing output to "//path(1:len_trim(path))//'data_'//trim(fieldtype)//trim(adjustl(tmp))//'.txt'
            exit
          end if
        end do
      else
        print '(a)',"Writing output to "//path(1:len_trim(path))//'data_'//trim(fieldtype)//'.txt'
      end if
!
!  write datatype of data contained in file
      write (unitid,'(a,/2x,a)') '% Datatype',datatype

!  write fieldtype of data contained in file
      write (unitid,'(a,/2x,a)') '% Fieldtype',fieldtype

!  write location of origin (always (0,0) for PolyDE)
      write (unitid,'(a)') '% Origin of data'
      write (unitid,'(2(1x,es15.8))') 0._DP, 0._DP

      select case (datatype)
      case ('GRID')
!  write number of grid points in x- and y-direction
        write (unitid,'(a)') '% Number of grid points in x- and y-direction'
        write (tmp ,'(i10)') size(xval)
        write (tmp1,'(i10)') size(yval)
        write (unitid,'(2x,a)') trim(adjustl(tmp))//' '//trim(adjustl(tmp1))

!  write header for grid of x- and y-values
        write (unitid,'(a)') '% Grid of '//trim(adjustl(tmp)) //        &
     &                       ' x- and '//trim(adjustl(tmp1)) //         &
     &                       ' y-values'
!  write x-values
        fmtreal = '('//trim(adjustl(tmp))//'(1x,es15.8)'//')'
        write (unitid,trim(fmtreal)) xval
!  write y-values
        fmtreal = '('//trim(adjustl(tmp1))//'(1x,es15.8)'//')'
        write (unitid,trim(fmtreal)) yval
!  write date and time of solution
        call date_and_time(dat,tim)
        write (tmp,'(a,g10.3)') '',100._DP*fem_accuracy
        write (unitid,'(a)') '% Data ('//trim(fieldtype)//')' //        &
     &        ', Relative Error in % ='//tmp //                         &
     &        ', Date = '//dat(7:8)//'.'//dat(5:6)//'.'//dat(1:4) //    &
     &        ', Time = '//tim(1:2)//':'//tim(3:4)
!  write field values as a gridx * gridy matrix
        write (tmp,'(i10)') size(zval(:,1))
        fmtreal = '('//trim(adjustl(tmp))//'(1x,es15.8)'//')'
        do i = 1, size(zval(1,:))
          write (unitid,trim(fmtreal)) zval(:,i)
        end do
!
      case ('SCATTERED')
!  write number of list entries
        write (unitid,'(a)') '% Number of list entries'
        write (tmp ,'(i10)') num
        write (unitid,'(2x,a)') trim(adjustl(tmp))
!  write date and time of solution
        call date_and_time(dat,tim)
        write (tmp,'(a,g10.3)') '',100._DP*fem_accuracy
        write (unitid,'(a)') '% Data ('//trim(fieldtype)//')' //        &
     &         ', Relative Error in % ='//tmp //                        &
     &         ', Date = '//dat(7:8)//'.'//dat(5:6)//'.'//dat(1:4) //   &
     &         ', Time = '//tim(1:2)//':'//tim(3:4)
!  write list of values (xval yval zval)
        do i = 1, num
          write (unitid,'(3(1x,g23.16))') xval(i),yval(i),zval(i,1)
        end do
!
      case default
! print error message if datatype unknown
        print '(a)',"***** Datatype "//trim(datatype)//" unknown."
        print '(a)',"No values exported. Check POSTsettings.txt"
      end select

      close (unitid)
      call grglun(unitid)

      deallocate (xval,yval,zval)
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      subroutine importdata(file,datatype,fieldtype,xval,yval,zval,errcode)
      use feminterface, only: getsetting, low2hi
      use femtypes
      implicit none
      integer (I4B) :: errcode
      real (DP), allocatable :: xval(:), yval(:), zval(:,:)
      character (len=*) :: datatype, fieldtype
      character (len=*) :: file
      intent (in)  :: file
      intent (out) :: datatype, fieldtype, xval, yval, zval
!
!
!  Subroutine to import field values as:
!  a) scattered data given in format xval yval zval (datataype = 'SCATTERED')
!  b) grid data on a regular grid of gridx * gridy points  (datataype = 'GRID')
!
!  from the file data.txt. The file contains the fieldtype and:
!  a) number of grid points in x- and y-direction (gridx and gridy), x-coords
!     (xval) and y-coords (yval) of points and matrix of grid values (zval).
!  b) number of total entries (num) and xval, yval, zval list of values
!
!  Comments in file are read and must occur at correct places!
!
!
!  Output:
!     datatype           'SCATTERED' or 'GRID' for corresponding type
!     errcode            20 if an error ocurred while reading
!     fieldtype          field quantity to compute (e.g. POTENTIAL, EY ...)
!     xval, yval         coordinates of grid points / list
!     zval               field values of grid points /list

!  Internal variables:
      integer (I4B) :: gridx, gridy, i, ierror=0, num, unitid
      real (DP) :: origin(2)
      character (len=200):: path, comment
!      
      external :: grglun
!
!  initialize errcode
      errcode = 0

!  read data from file
      call grglun(unitid)
      call getsetting('PROJECTPATH',path)
!  open a file in the specified project path
      open (unitid,FILE=path(1:len_trim(path))//trim(file),STATUS='OLD',&
     &      FORM='formatted',POSITION='REWIND',ACTION='READ',IOSTAT=ierror)
      if (ierror .ne. 0) then
        print '(a)',"***** File "//path(1:len_trim(path))//trim(file)//" does not exist."
        errcode = 20
        return
      end if

!  read datatype
      read (unitid,*) comment
      read (unitid,*) datatype
      datatype = trim(adjustl(datatype))
      call low2hi(datatype,len_trim(datatype))

!  read fieldtype of data contained in file
      read (unitid,*) comment
      read (unitid,*) fieldtype
      fieldtype = trim(adjustl(fieldtype))
      call low2hi(fieldtype,len_trim(fieldtype))

!  read location of origin (always (0,0) for PolyDE)
      read (unitid,*) comment
      read (unitid,*) origin(:)

!  select case for datatype contained in file
      select case (datatype)
      case ('GRID')
!  read number of grid points in x- and y-direction
        read (unitid,*) comment
        read (unitid,*) gridx, gridy
!  allocate arrays according to number of points
        allocate (xval(gridx), yval(gridy), zval(gridx,gridy))
!  read header for grid of x- and y-values
        read (unitid,*) comment
!  read x-values
        read (unitid,*) xval
        xval = xval + origin(1)
!  read y-values
        read (unitid,*) yval
        yval = yval + origin(2)
!  read comment for solution
        read (unitid,*) comment
!  read field values as a gridx * gridy matrix
        do i = 1, gridy
          read (unitid,*) zval(:,i)
        end do
!
      case ('SCATTERED')
!  read number of list entries
        read (unitid,*) comment
        read (unitid,*) num
!  allocate arrays according to number of list entries
        allocate (xval(num), yval(num), zval(num,1))
!  read comment for solution
        read (unitid,*) comment
!  read field values as list of num entries
        do i = 1, num
          read (unitid,*) xval(i),yval(i),zval(i,1)
        end do
!  modify list entries by origin
        xval = xval + origin(1)
        yval = yval + origin(2)
!
      case default
! return error code 20 if datatype unknown
        print '(a)',"***** Datatype "//trim(datatype)//" unknown."
        print '(a)',"Check file "//path(1:len_trim(path))//trim(file)
        errcode = 20
      end select

      close (unitid)
      call grglun(unitid)
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      subroutine initdata()
      use feminterface, only: getsetting
      use femtypes
      implicit none
!
!
!  Subroutine initdata will read required data depending on selected multi-
!  physics mode. Data is stored in logical variables and pointer arrays made
!  available by mpmodule.
!
!  ONLY TO COUPLING WORKS up to now

!  Internal variables:
      character (len=10) :: str
!
!
!  get multiphysics setting from FEMsettings.txt
      call getsetting('MULTIPHYSICS',str)
      if (str .eq. 'YES') then
!  set multiphysics variable
        multiphysics = .true.
!  get kind of multiphysics coupling (kind of coupling must be in alphabetical
!  order with O(ptical) being the last letter)
        call getsetting('MP_COUPLING',str)
!  begin case structure to set variables and read data
        select case (str)
!  E(lectro)-O(ptical) coupling only
          case ('EO')
            eomod = .true.
            etmod = .false.
            momod = .false.
            tomod = .false.
            call loadelectricfield()
!  E(lectro)-T(hermal) coupling only
          case ('ET')
            eomod = .false.
            etmod = .true.
            momod = .false.
            tomod = .false.
            call loadpowerdensity()
!  M(echano)-O(ptical) coupling only (elasto-optic or photo-elatic)
          case ('MO')
            eomod = .false.
            etmod = .false.
            momod = .true.
            tomod = .false.
            call loadstrain()
!  T(hermo)-O(ptical) coupling only
          case ('TO')
            eomod = .false.
            etmod = .false.
            momod = .false.
            tomod = .true.
            call loadtemperature()
!  E(lectro)-M(echano)-O(ptical) coupling
          case ('EMO')
            eomod = .true.
            etmod = .false.
            momod = .true.
            tomod = .false.
            call loadelectricfield()
            call loadstrain()
!  E(lectro)-T(hermo)-O(ptical) coupling
          case ('ETO')
            eomod = .true.
            etmod = .false.
            momod = .false.
            tomod = .true.
            call loadelectricfield()
            call loadtemperature()
!  E(lectro)-M(echano)-T(hermo)-O(ptical) coupling
          case ('EMTO')
            eomod = .true.
            etmod = .false.
            momod = .true.
            tomod = .true.
            call loadelectricfield()
            call loadstrain()
            call loadtemperature()
!  M(echano)-T(hermo)-O(ptical) coupling
          case ('MTO')
            eomod = .false.
            etmod = .false.
            momod = .true.
            tomod = .true.
            call loadstrain()
            call loadtemperature()
        end select
!  if all data were not read successfully, set multiphysics back to .false.
        if ((.not. eomod) .and. (.not. etmod) .and. (.not. momod) .and. (.not. tomod)) then
          print *,"Warning: Multiphysics will be ignored!"
          multiphysics = .false.
        end if
      else
!  set multiphysics variable to .false. for other strings than 'YES'
        multiphysics = .false.
      end if
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      subroutine loadelectricfield()
      use feminterface, only: 
      use femtypes
      implicit none
!
!
!  This subroutine successively loads electric field strength data from files
!  data_EX.txt, data_EY.txt and data_EZ.txt and stores it in arrays xE, yE and
!  zE. The program returns by setting eomod to .false. if something went
!  wrong reading the first to files or contained fieldtype is false. If there is
!  no data for individual components, these components will be set to 0._DP. We
!  assume, that all data has the same datatype and field values are given at the
!  same coordinates!
!
!  Internal variables:
      integer (I4B) :: errcode
      real (DP), allocatable :: xval(:), yval(:), zval(:,:)
      character (len=10) :: fieldtype
      character (len=16) :: file
      logical :: err1st = .false. , err2nd = .false. , err3rd = .false.
!      
!
!  File name for electric field data in x-direction MUST be data_EX.txt
      file = 'data_EX.txt     '
      call importdata(file,datatypeE,fieldtype,xval,yval,zval,errcode)
!  set err1st = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: EX will be ignored!"
        err1st = .true.
      else if (fieldtype .ne. 'EX') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with electro-optic coupling."
        print '(a)',"Warning: EX will be ignored!"
        err1st = .true.
      else
!  allocate common arrays xE, yE and zE depending on datatypeE and copy imported
!  values xval, yval and zval to them
        select case (datatypeE)
!  CAUTION: WE ASSUME, THAT ALL DATA (X,Y,Z-DIRECTION) WAS WRITTEN WITH THE SAME
!           DATATYPE AND THE SAME COORDINATES!
        case ('GRID')
          allocate( xE(size(xval)), yE(size(yval)), zE(3,size(xval),size(yval)) )
        case ('SCATTERED')
          allocate( xE(size(xval)), yE(size(yval)), zE(3,size(xval),1) )
        end select
!  initialize arrays (fill with zero entries)
        xE = 0._DP
        yE = 0._DP
        zE = 0._DP
!  everything is fine, so assign coordinates and electric field strength with
!  respect to x-direction
        xE(:) = xval
        yE(:) = yval
        zE(1,:,:) = zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  NOW LOAD ELECTRIC FIELD WITH RESPECT TO Y-DIRECTION
      file = 'data_EY.txt     '
      call importdata(file,datatypeE,fieldtype,xval,yval,zval,errcode)
!  set err2nd = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: EY will be ignored!"
        err2nd = .true.
      else if (fieldtype .ne. 'EY') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with electro-optic coupling."
        print '(a)',"Warning: EY will be ignored!"
        err2nd = .true.
      else
!  allocate common arrays xE, yE and zE depending on datatypeE and copy imported
!  values xval, yval and zval to them
        if (.not. allocated(xE)) then
          select case (datatypeE)
!  CAUTION: WE ASSUME, THAT ALL DATA (X,Y,Z-DIRECTION) WAS WRITTEN WITH THE SAME
!           DATATYPE AND THE SAME COORDINATES!
          case ('GRID')
            allocate( xE(size(xval)), yE(size(yval)), zE(3,size(xval),size(yval)) )
          case ('SCATTERED')
            allocate( xE(size(xval)), yE(size(yval)), zE(3,size(xval),1) )
          end select
!  initialize arrays (fill with zero entries)
          xE = 0._DP
          yE = 0._DP
          zE = 0._DP
!  everything is fine, so assign coordinates
          xE(:) = xval
          yE(:) = yval
        end if
!  assign Ey field strength
        zE(2,:,:) = zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  NOW LOAD ELECTRIC FIELD WITH RESPECT TO Z-DIRECTION
      file = 'data_EZ.txt     '
      call importdata(file,datatypeE,fieldtype,xval,yval,zval,errcode)
!  set err3rd = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: EZ will be ignored!"
        err3rd = .true.
      else if (fieldtype .ne. 'EZ') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with electro-optic coupling."
        print '(a)',"Warning: EZ will be ignored!"
        err3rd = .true.
      else
!  allocate common arrays xE, yE and zE depending on datatypeE and copy imported
!  values xval, yval and zval to them
        if (.not. allocated(xE)) then
          select case (datatypeE)
!  CAUTION: WE ASSUME, THAT ALL DATA (X,Y,Z-DIRECTION) WAS WRITTEN WITH THE SAME
!           DATATYPE AND THE SAME COORDINATES!
          case ('GRID')
            allocate( xE(size(xval)), yE(size(yval)), zE(3,size(xval),size(yval)) )
          case ('SCATTERED')
            allocate( xE(size(xval)), yE(size(yval)), zE(3,size(xval),1) )
          end select
!  initialize arrays (fill with zero entries)
          xE = 0._DP
          yE = 0._DP
          zE = 0._DP
!  everything is fine, so assign coordinates
          xE(:) = xval
          yE(:) = yval
        end if
!  assign Ez field strength
        zE(3,:,:) = zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  set electro-optic modification to false, if all data was not read successfully
      if (err1st .and. err2nd .and. err3rd) then
        eomod = .false.
      else
        print '(a)',"Electric field data successfully loaded."
      end if
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      subroutine loadstrain()
      use feminterface, only: 
      use femtypes
      implicit none
!
!
!  This subroutine successively loads strain field components of a symmetric
!  strain tensor from files data_Sxx.txt, data_Syy.txt, data_Szz.txt,
!  data_Syz.txt, data_Sxz.txt and data_Sxy.txt and stores it in arrays xM, yM
!  and zM. Due to internal storage in matrix/engineering notation, a factor of
!  two will be applied to entries zM(4:5,:,:). The program returns by setting
!  momod to .false. if something went wrong reading all files or contained
!  fieldtype is false for all files. If there is no data for individual
!  components, this component is set to 0._DP. We assume, that all data has the
!  same datatype and field values are given at the same coordinates!
!
!  Internal variables:
      integer (I4B) :: errcode
      real (DP), allocatable :: xval(:), yval(:), zval(:,:)
      character (len=10) :: fieldtype
      character (len=16) :: file
      logical :: err1st = .false. , err2nd = .false. , err3rd = .false., &
                 err4th = .false. , err5th = .false. , err6th = .false.
!      
!
!  File name for strain data in x-direction MUST be data_SX.txt
      file = 'data_SX.txt     '
      call importdata(file,datatypeM,fieldtype,xval,yval,zval,errcode)
!  set err1st = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: SX will be ignored!"
        err1st = .true.
      else if (fieldtype .ne. 'SX') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with elasto-optic coupling."
        print '(a)',"Warning: SX will be ignored!"
        err1st = .true.
      else
!  allocate common arrays xM, yM and zM depending on datatypeM and copy imported
!  values xval, yval and zval to them
        select case (datatypeE)
!  CAUTION: WE ASSUME, THAT ALL DATA (X,Y,Z-DIRECTION) WAS WRITTEN WITH THE SAME
!           DATATYPE AND THE SAME COORDINATES!
        case ('GRID')
          allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),size(yval)) )
        case ('SCATTERED')
          allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),1) )
        end select
!  initialize arrays (fill with zero entries)
        xM = 0._DP
        yM = 0._DP
        zM = 0._DP
!  everything is fine, so assign coordinates and strain with respect to x
        xM(:) = xval
        yM(:) = yval
        zM(1,:,:) = zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  NOW LOAD STRAIN WITH RESPECT TO Y-DIRECTION
      file = 'data_SY.txt     '
      call importdata(file,datatypeM,fieldtype,xval,yval,zval,errcode)
!  set err2nd = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: SY will be ignored!"
        err2nd = .true.
      else if (fieldtype .ne. 'SY') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with elasto-optic coupling."
        print '(a)',"Warning: SY will be ignored!"
        err2nd = .true.
      else
!  allocate common arrays xM, yM and zM depending on datatypeE and copy imported
!  values xval, yval and zval to them
        if (.not. allocated(xM)) then
          select case (datatypeM)
!  CAUTION: WE ASSUME, THAT ALL DATA (X,Y,Z-DIRECTION) WAS WRITTEN WITH THE SAME
!           DATATYPE AND THE SAME COORDINATES!
          case ('GRID')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),size(yval)) )
          case ('SCATTERED')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),1) )
          end select
!  initialize arrays (fill with zero entries)
          xM = 0._DP
          yM = 0._DP
          zM = 0._DP
!  everything is fine, so assign coordinates
          xM(:) = xval
          yM(:) = yval
        end if
!  assign strain with respect to y
        zM(2,:,:) = zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  NOW LOAD STRAIN WITH RESPECT TO Z-DIRECTION
      file = 'data_SZ.txt     '
      call importdata(file,datatypeM,fieldtype,xval,yval,zval,errcode)
!  set err3rd = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: SZ will be ignored!"
        err3rd = .true.
      else if (fieldtype .ne. 'SZ') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with elasto-optic coupling."
        print '(a)',"Warning: SZ will be ignored!"
        err3rd = .true.
      else
!  allocate common arrays xM, yM and zM depending on datatypeM and copy imported
!  values xval, yval and zval to them
        if (.not. allocated(xM)) then
          select case (datatypeE)
!  CAUTION: WE ASSUME, THAT ALL DATA (X,Y,Z-DIRECTION) WAS WRITTEN WITH THE SAME
!           DATATYPE AND THE SAME COORDINATES!
          case ('GRID')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),size(yval)) )
          case ('SCATTERED')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),1) )
          end select
!  initialize arrays (fill with zero entries)
          xM = 0._DP
          yM = 0._DP
          zM = 0._DP
!  everything is fine, so assign coordinates
          xM(:) = xval
          yM(:) = yval
        end if
!  assign strain with respect to z
        zM(3,:,:) = zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  NOW LOAD STRAIN WITH RESPECT TO YZ-DIRECTION
      file = 'data_SYZ.txt    '
      call importdata(file,datatypeM,fieldtype,xval,yval,zval,errcode)
!  set err4th = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: SYZ will be ignored!"
        err4th = .true.
      else if (fieldtype .ne. 'SYZ') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with elasto-optic coupling."
        print '(a)',"Warning: SYZ will be ignored!"
        err4th = .true.
      else
!  allocate common arrays xM, yM and zM depending on datatypeM and copy imported
!  values xval, yval and zval to them
        if (.not. allocated(xM)) then
          select case (datatypeE)
!  CAUTION: WE ASSUME, THAT ALL DATA (X,Y,Z-DIRECTION) WAS WRITTEN WITH THE SAME
!           DATATYPE AND THE SAME COORDINATES!
          case ('GRID')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),size(yval)) )
          case ('SCATTERED')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),1) )
          end select
!  initialize arrays (fill with zero entries)
          xM = 0._DP
          yM = 0._DP
          zM = 0._DP
!  everything is fine, so assign coordinates
          xM(:) = xval
          yM(:) = yval
        end if
!  assign strain with respect to yz (notice the factor of 2 for tensor to matrix notation)
        zM(4,:,:) = 2._DP*zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  NOW LOAD STRAIN WITH RESPECT TO XZ-DIRECTION
      file = 'data_SXZ.txt    '
      call importdata(file,datatypeM,fieldtype,xval,yval,zval,errcode)
!  set err5th = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: SXZ will be ignored!"
        err5th = .true.
      else if (fieldtype .ne. 'SXZ') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with elasto-optic coupling."
        print '(a)',"Warning: SXZ will be ignored!"
        err5th = .true.
      else
!  allocate common arrays xM, yM and zM depending on datatypeM and copy imported
!  values xval, yval and zval to them
        if (.not. allocated(xM)) then
          select case (datatypeE)
!  CAUTION: WE ASSUME, THAT ALL DATA (X,Y,Z-DIRECTION) WAS WRITTEN WITH THE SAME
!           DATATYPE AND THE SAME COORDINATES!
          case ('GRID')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),size(yval)) )
          case ('SCATTERED')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),1) )
          end select
!  initialize arrays (fill with zero entries)
          xM = 0._DP
          yM = 0._DP
          zM = 0._DP
!  everything is fine, so assign coordinates
          xM(:) = xval
          yM(:) = yval
        end if
!  assign strain with respect to xz (notice the factor of 2 for tensor to matrix notation)
        zM(5,:,:) = 2._DP*zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  NOW LOAD STRAIN WITH RESPECT TO XY-DIRECTION
      file = 'data_SXY.txt    '
      call importdata(file,datatypeM,fieldtype,xval,yval,zval,errcode)
!  set err6th = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: SXY will be ignored!"
        err6th = .true.
      else if (fieldtype .ne. 'SXY') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with elasto-optic coupling."
        print '(a)',"Warning: SXY will be ignored!"
        err6th = .true.
      else
!  allocate common arrays xM, yM and zM depending on datatypeM and copy imported
!  values xval, yval and zval to them
        if (.not. allocated(xM)) then
          select case (datatypeE)
!  CAUTION: WE ASSUME, THAT ALL DATA (X,Y,Z-DIRECTION) WAS WRITTEN WITH THE SAME
!           DATATYPE AND THE SAME COORDINATES!
          case ('GRID')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),size(yval)) )
          case ('SCATTERED')
            allocate( xM(size(xval)), yM(size(yval)), zM(6,size(xval),1) )
          end select
!  initialize arrays (fill with zero entries)
          xM = 0._DP
          yM = 0._DP
          zM = 0._DP
!  everything is fine, so assign coordinates
          xM(:) = xval
          yM(:) = yval
        end if
!  assign strain with respect to xy (notice the factor of 2 for tensor to matrix notation)
        zM(6,:,:) = 2._DP*zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  set elasto-optic modification to false, if all data was not read successfully
      if (err1st .and. err2nd .and. err3rd .and. &
          err4th .and. err5th .and. err6th) then
        momod = .false.
      else
        print '(a)',"Strain data successfully loaded."
      end if
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      subroutine loadpowerdensity()
      use feminterface, only: 
      use femtypes
      implicit none
!
!
!  This subroutine loads power density data from file data_P.txt and stores it
!  in arrays xP, yP and zP.
!
!  Internal variables:
      integer (I4B) :: errcode
      real (DP), allocatable :: xval(:), yval(:), zval(:,:)
      character (len=10) :: fieldtype
      character (len=16) :: file
      logical :: err1st = .false.
!      
!
!  File name for power density field data MUST be data_P.txt
      file = 'data_P.txt      '
      call importdata(file,datatypeP,fieldtype,xval,yval,zval,errcode)
!  set err1st = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: power density data will be ignored!"
        err1st = .true.
      else if (fieldtype .ne. 'P') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with electro-thermal coupling."
        print '(a)',"Warning: power density data will be ignored!"
        err1st = .true.
      else
!  allocate common arrays xP, yP and zP depending on datatypeP and copy imported
!  values xval, yval and zval to them
        select case (datatypeP)
        case ('GRID')
          allocate( xP(size(xval)), yP(size(yval)), zP(size(xval),size(yval)) )
        case ('SCATTERED')
          allocate( xP(size(xval)), yP(size(yval)), zP(size(xval),1) )
        end select
        xP = xval
        yP = yval
        zP = zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  set tomod to false, if temperature data was not read successfully
      if (err1st) then
        etmod = .false.
      else
        print '(a)',"Power density data successfully loaded."
      end if
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      subroutine loadtemperature()
      use feminterface, only: 
      use femtypes
      implicit none
!
!
!  This subroutine loads temperature data from file data_T.txt and stores it
!  in arrays xT, yT and zT.
!
!  Internal variables:
      integer (I4B) :: errcode
      real (DP), allocatable :: xval(:), yval(:), zval(:,:)
      character (len=10) :: fieldtype
      character (len=16) :: file
      logical :: err1st = .false.
!      
!
!  File name for thermal field data MUST be data_T.txt
      file = 'data_T.txt      '
      call importdata(file,datatypeT,fieldtype,xval,yval,zval,errcode)
!  set err1st = .true. if something went wrong with data import or data contained
!               in file does not match required data
      if (errcode .eq. 20) then
        print '(a)',"Warning: temperature data will be ignored!"
        err1st = .true.
      else if (fieldtype .ne. 'T') then
        print '(a)',"Fieldtype in "//trim(file)//" is "//trim(fieldtype)//" with thermo-optic coupling."
        print '(a)',"Warning: temperature data will be ignored!"
        err1st = .true.
      else
!  allocate common arrays xT, yT and zT depending on datatypeT and copy imported
!  values xval, yval and zval to them
        select case (datatypeT)
        case ('GRID')
          allocate( xT(size(xval)), yT(size(yval)), zT(size(xval),size(yval)) )
        case ('SCATTERED')
          allocate( xT(size(xval)), yT(size(yval)), zT(size(xval),1) )
        end select
        xT = xval
        yT = yval
        zT = zval
!  deallocate imported arrays (data is stored in common arrays, now)
        deallocate(xval,yval,zval)
      end if
!
!  set tomod to false, if temperature data was not read successfully
      if (err1st) then
        tomod = .false.
      else
        print '(a)',"Temperature data successfully loaded."
      end if
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      function locateparam(str)
      use feminterface, only: 
      use femtypes
      use matconstants
      implicit none
      integer (I4B) :: locateparam
      character (len=*) :: str
!
!
!  Locate position of str in parameternames array.
!
!  internal variables
      integer (I4B) :: i
!
!
      do i = 1, numparam
        if (str .eq. parameternames(i)) then
          locateparam = i
          exit
        end if
      end do
!
      end function
!------------------------------------------------------------------------------
!
!
!
      subroutine matmodify(list,xs,ys,matten)
      use feminterface, only: idbvip, invertmat3, mat2ten_r2,           &
     &                        mat2ten_r3, mat2ten_r4, ten2mat_r2,       &
     &                        ten2mat_r3, ten2mat_r4, transten_r1,      &
     &                        transten_r2, transten_r3, transten_r4
      use femtypes
      use globalvariables, only:
      implicit none
      real(DP) :: xs, ys, list(:)
      complex (DPC) :: matten(3,3)
      intent (in) :: list,xs,ys
      intent (inout) :: matten
!
!
!  Subroutine matmodify modifies a given 3x3 material tensor depending on kind
!  of multiphysics coupling. It is assumed, that necessary data from previous
!  simulations is already stored in common arrays in mpmodule. The modification
!  will be done for all couplings required in successive order.
!  Field values will be interpolated from common arrays depending on the data- 
!  type ('GRID' or 'SCATTERED'). Necessary constants or coefficients will be
!  fetched from list in material files (e.g. thermo-optical coefficient). Then
!  individual entries or whole array is modified.
!
!  ONLY TO COUPLING WORKS up to now
!
!
!  Input:
!     list               contains data from material files (refer to setstandard-
!                        values.f90 for individual position of entries)
!     xs, ys             coordinate of point with material modification
!
!  Output:
!     matten             3x3 material tensor which contains complex relative per-
!                        mittivity at point (xs,ys)
!
!  Internal variables:
      integer (I4B) :: i, j, pos, errcode, nip
      real (DP), allocatable :: xstmp(:), ystmp(:), zs(:)
      real (DP) :: deltaB_Pockels(6), deltaB_Kerr(6), deltaB_Elasto(6)
      real (DP) :: ten_r2(3,3)
      real (DP) :: tocoeff, matrix_Pockels(6,3), matrix_Kerr(6,6)
      real (DP) :: matrix_Piezo(3,6), matrix_Elasto(6,6), matrix_ThEx(6)
      real (DP) :: Efield(3), Efield2(6)
      real (DP) :: Strain_Ext(6), Strain_Piezo(6), Strain_Temp(6)
      complex (DPC) :: B(6), invten(3,3)
      character (len=20) :: str
      logical, save :: eomessage, momessage, tomessage
      logical, save :: eointpolstarted, mointpolstarted, tointpolstarted
!
!
!
!  check, if there is any modification to do
      if (eomod .or. momod .or. tomod) then
!  initialize interpolation point
        nip = 1
        allocate (xstmp(nip),ystmp(nip),zs(nip))
        xstmp = xs
        ystmp = ys
!
!  initialize data
        deltaB_Pockels = 0._DP
        deltaB_Kerr    = 0._DP
        deltaB_Elasto  = 0._DP
        Efield         = 0._DP
        Efield2        = 0._DP
        Strain_Ext     = 0._DP
        Strain_Piezo   = 0._DP
        Strain_Temp    = 0._DP
!
!  GET COEFFICIENTS
!------------------------------------------------------------------------------
!  locate position of thermo-optic coefficient
        str = 'DEPS_BY_DT'
        pos = locateparam(str)
!  assign thermo-optic coefficient
        tocoeff = list(pos)
!
!  locate position of first linear electro-optic coefficient r11 (matrix notation)
        str = 'R11'
        pos = locateparam(str)
!  assign all linear (Pockels) electro-optic coefficients r (matrix notation)
!  it is assumed, that there are only 18 independent linear electrooptic constants
        matrix_Pockels(1,:) = (/list(pos),    list(pos+1),  list(pos+2)/)
        matrix_Pockels(2,:) = (/list(pos+3),  list(pos+4),  list(pos+5)/)
        matrix_Pockels(3,:) = (/list(pos+6),  list(pos+7),  list(pos+8)/)
        matrix_Pockels(4,:) = (/list(pos+9),  list(pos+10), list(pos+11)/)
        matrix_Pockels(5,:) = (/list(pos+12), list(pos+13), list(pos+14)/)
        matrix_Pockels(6,:) = (/list(pos+15), list(pos+16), list(pos+17)/)
!
!  locate position of first elasto-optic coefficient p11 (matrix notation)
        str = 'K11'
        pos = locateparam(str)
!  assign all quadratic (Kerr) electro-optic coefficients k (matrix notation)
!  it is assumed, that there are 36 independent quadratic electrooptic constants
        matrix_Kerr(1,:) = (/list(pos),    list(pos+1),  list(pos+2),  list(pos+3),  list(pos+4),  list(pos+5)/)
        matrix_Kerr(2,:) = (/list(pos+6),  list(pos+7),  list(pos+8),  list(pos+9),  list(pos+10), list(pos+11)/)
        matrix_Kerr(3,:) = (/list(pos+12), list(pos+13), list(pos+14), list(pos+15), list(pos+16), list(pos+17)/)
        matrix_Kerr(4,:) = (/list(pos+18), list(pos+19), list(pos+20), list(pos+21), list(pos+22), list(pos+23)/)
        matrix_Kerr(5,:) = (/list(pos+24), list(pos+25), list(pos+26), list(pos+27), list(pos+28), list(pos+29)/)
        matrix_Kerr(6,:) = (/list(pos+30), list(pos+31), list(pos+32), list(pos+33), list(pos+34), list(pos+35)/)
!
!  locate position of first elasto-optic coefficient p11 (matrix notation)
        str = 'P11'
        pos = locateparam(str)
!  assign all elasto-optic coefficients p (matrix notation)
!  it is assumed, that there are 36 independent elasto-optic constants
        matrix_Elasto(1,:) = (/list(pos),    list(pos+1),  list(pos+2),  list(pos+3),  list(pos+4),  list(pos+5)/)
        matrix_Elasto(2,:) = (/list(pos+6),  list(pos+7),  list(pos+8),  list(pos+9),  list(pos+10), list(pos+11)/)
        matrix_Elasto(3,:) = (/list(pos+12), list(pos+13), list(pos+14), list(pos+15), list(pos+16), list(pos+17)/)
        matrix_Elasto(4,:) = (/list(pos+18), list(pos+19), list(pos+20), list(pos+21), list(pos+22), list(pos+23)/)
        matrix_Elasto(5,:) = (/list(pos+24), list(pos+25), list(pos+26), list(pos+27), list(pos+28), list(pos+29)/)
        matrix_Elasto(6,:) = (/list(pos+30), list(pos+31), list(pos+32), list(pos+33), list(pos+34), list(pos+35)/)
!
!  locate position of first linear electro-optic coefficient r11 (matrix notation)
        str = 'D11'
        pos = locateparam(str)
!  assign all piezo-electric coefficients d (matrix notation)
!  it is assumed, that there are only 18 independent piezo-electric constants
        matrix_Piezo(1,:) = (/list(pos),    list(pos+1),  list(pos+2),  list(pos+3),  list(pos+4),  list(pos+5)/)
        matrix_Piezo(2,:) = (/list(pos+6),  list(pos+7),  list(pos+8),  list(pos+9),  list(pos+10), list(pos+11)/)
        matrix_Piezo(3,:) = (/list(pos+12), list(pos+13), list(pos+14), list(pos+15), list(pos+16), list(pos+17)/)
!
!  locate position of first thermal expansion coefficient alpha1 (matrix notation)
        str = 'ALPHA1'
        pos = locateparam(str)
!  assign all thermal expansion coefficients
!  it is assumed, that there are 6 independent thermal expansion coefficients
        matrix_ThEx = (/list(pos), list(pos+1), list(pos+2), list(pos+3), list(pos+4), list(pos+5)/)
!
!
!  GET FIELD DATA
!------------------------------------------------------------------------------
!  get thermal field data and compute new material tensor
        if ( tomod .and. (xs .ge. minval(xT)) .and. (xs .le. maxval(xT)) .and. &
                         (ys .ge. minval(yT)) .and. (ys .le. maxval(yT)) ) then
          if (.not. tomessage) then
            print "(a)","INFO: Thermo-optic material modification."
            tomessage = .true.
          end if
          select case (datatypeT)
          case ('GRID')
            call bilintpol(xT,yT,zT,xs,ys,zs(1),errcode)
          case ('SCATTERED')
            if (.not. tointpolstarted) then
              call idbvip (1,size(xT),xT,yT,zT(:,1),nip,xstmp,ystmp,zs)
              tointpolstarted = .true.
            else
              call idbvip (2,size(xT),xT,yT,zT(:,1),nip,xstmp,ystmp,zs)
            end if
          end select
!  modify all entries in the permittivity tensor
          do i = 1,3
            do j = 1,3
!  we assume, that if former permittivity entry was zero, new entry may deviate
!  from zero after modification
!
!  compute modified permittivity (thermo-optic effect is assumed to only influence real part)
!  REMARK: due to this, PMLs might not work perfectly any more
              matten(i,j) = cmplx( (sqrt(real(matten(i,j),DP)) + sqrt(tocoeff)*zs(1))**2, aimag(matten(i,j)), DPC )
            end do
          end do
!
!  Compute Strain due to thermal expansion and change in impermeability
          Strain_Temp    = zs(1)*matrix_ThEx
          deltaB_Elasto  = matmul(matrix_Elasto,Strain_Temp)
        end if
!
!  INVERT MATERIAL TENSOR FOR MODIFICATION OF RELATIVE DIELECTRIC IMPERMEABILITY
!------------------------------------------------------------------------------
!  vector B = 1 / eps_r (matrix notation)
        call invertmat3(matten,invten)
        B = ten2mat_r2(invten)
!
!------------------------------------------------------------------------------
!  get electric field data and compute deltaB due to electro-optic effects
        if ( eomod .and. (xs .ge. minval(xE)) .and. (xs .le. maxval(xE)) .and. &
                       (ys .ge. minval(yE)) .and. (ys .le. maxval(yE)) ) then
!  store E-field in vector (vector(1) = Ex, vector(2) = Ey, vector(3) = Ez)
          if (.not. eomessage) then
            print "(a)","INFO: Material modification due to electro-optic effect."
            eomessage = .true.
          end if
          select case (datatypeE)
          case ('GRID')
!  allocate temporary array for subroutine blintpol
            call bilintpol(xE,yE,zE(1,:,:),xs,ys,Efield(1),errcode)
            call bilintpol(xE,yE,zE(2,:,:),xs,ys,Efield(2),errcode)
            call bilintpol(xE,yE,zE(3,:,:),xs,ys,Efield(3),errcode)
          case ('SCATTERED')
            if (.not. eointpolstarted) then
              call idbvip (1,size(xE),xE,yE,zE(1,:,1),nip,xstmp,ystmp,Efield(1))
              call idbvip (3,size(xE),xE,yE,zE(2,:,1),nip,xstmp,ystmp,Efield(2))
              call idbvip (3,size(xE),xE,yE,zE(3,:,1),nip,xstmp,ystmp,Efield(3))
              eointpolstarted = .true.
            else
              call idbvip (2,size(xE),xE,yE,zE(1,:,1),nip,xstmp,ystmp,Efield(1))
              call idbvip (3,size(xE),xE,yE,zE(2,:,1),nip,xstmp,ystmp,Efield(2))
              call idbvip (3,size(xE),xE,yE,zE(3,:,1),nip,xstmp,ystmp,Efield(3))
            end if
          end select
!  execute rank 1 tensor transformation from simulation to crystal coordinate system
          Efield = transten_r1(transpose(trmat),Efield)
!  compute Efield2 for Kerr effect
          Efield2(1) = Efield(1)**2
          Efield2(2) = Efield(2)**2
          Efield2(3) = Efield(3)**2
          Efield2(4) = Efield(2)*Efield(3)
          Efield2(5) = Efield(3)*Efield(1)
          Efield2(6) = Efield(1)*Efield(2)
!
!  Compute change in impermeability due to Pockels and Kerr effect
          deltaB_Pockels = matmul(matrix_Pockels,Efield)
          deltaB_Kerr    = matmul(matrix_Kerr,Efield2)
!
!  Compute Strain due to converse piezo-electric effect and change in impermeability
          Strain_Piezo   = matmul(transpose(matrix_Piezo),Efield)
          deltaB_Elasto  = deltaB_Elasto + matmul(matrix_Elasto,Strain_Piezo)
        end if
!
!  get strain field data and compute deltaB due to elasto-optic effects
        if ( momod .and. (xs .ge. minval(xM)) .and. (xs .le. maxval(xM)) .and. &
                         (ys .ge. minval(yM)) .and. (ys .le. maxval(yM)) ) then
!  store S(train)-field in zs vector (zs(1) = Sxx, zs(2) = Syy, zs(3) = Szz
!                                     zs(4) = Syz, zs(5) = Sxz, zs(6) = Sxy)
          if (.not. momessage) then
            print "(a)","INFO: Material modification due to elasto-optic effect."
            momessage = .true.
          end if
          select case (datatypeM)
          case ('GRID')
!  allocate temporary array for subroutine blintpol
            call bilintpol(xM,yM,zM(1,:,:),xs,ys,Strain_Ext(1),errcode)
            call bilintpol(xM,yM,zM(2,:,:),xs,ys,Strain_Ext(2),errcode)
            call bilintpol(xM,yM,zM(3,:,:),xs,ys,Strain_Ext(3),errcode)
            call bilintpol(xM,yM,zM(4,:,:),xs,ys,Strain_Ext(4),errcode)
            call bilintpol(xM,yM,zM(5,:,:),xs,ys,Strain_Ext(5),errcode)
            call bilintpol(xM,yM,zM(6,:,:),xs,ys,Strain_Ext(6),errcode)
          case ('SCATTERED')
            if (.not. mointpolstarted) then
              call idbvip (1,size(xM),xM,yM,zM(1,:,1),nip,xstmp,ystmp,Strain_Ext(1))
              call idbvip (3,size(xM),xM,yM,zM(2,:,1),nip,xstmp,ystmp,Strain_Ext(2))
              call idbvip (3,size(xM),xM,yM,zM(3,:,1),nip,xstmp,ystmp,Strain_Ext(3))
              call idbvip (3,size(xM),xM,yM,zM(4,:,1),nip,xstmp,ystmp,Strain_Ext(4))
              call idbvip (3,size(xM),xM,yM,zM(5,:,1),nip,xstmp,ystmp,Strain_Ext(5))
              call idbvip (3,size(xM),xM,yM,zM(6,:,1),nip,xstmp,ystmp,Strain_Ext(6))
              mointpolstarted = .true.
            else
              call idbvip (2,size(xE),xM,yM,zM(1,:,1),nip,xstmp,ystmp,Strain_Ext(1))
              call idbvip (3,size(xE),xM,yM,zM(2,:,1),nip,xstmp,ystmp,Strain_Ext(2))
              call idbvip (3,size(xE),xM,yM,zM(3,:,1),nip,xstmp,ystmp,Strain_Ext(3))
              call idbvip (3,size(xE),xM,yM,zM(4,:,1),nip,xstmp,ystmp,Strain_Ext(4))
              call idbvip (3,size(xE),xM,yM,zM(5,:,1),nip,xstmp,ystmp,Strain_Ext(5))
              call idbvip (3,size(xE),xM,yM,zM(6,:,1),nip,xstmp,ystmp,Strain_Ext(6))
            end if
          end select
!  matrix notation to rank 2 tensor notation
!  NOTICE: S_ij = Sm,     for m being 1, 2, 3
!          S_ij = 0.5*Sm, for m being 4, 5, 6
          ten_r2(1,:) = (/Strain_Ext(1),0.5_DP*Strain_Ext(6),0.5_DP*Strain_Ext(5)/)
          ten_r2(2,:) = (/0.5_DP*Strain_Ext(6),Strain_Ext(2),0.5_DP*Strain_Ext(4)/)
          ten_r2(3,:) = (/0.5_DP*Strain_Ext(5),0.5_DP*Strain_Ext(4),Strain_Ext(3)/)
!  execute rank 2 tensor transformation from simulation to crystal coordinate system
          ten_r2 = transten_r2(transpose(trmat),ten_r2)
!  rank 2 tensor notation to 1x6 matrix notation
          Strain_Ext = (/ten_r2(1,1),ten_r2(2,2),ten_r2(3,3),2._DP*ten_r2(2,3),2._DP*ten_r2(1,3),2._DP*ten_r2(1,2)/)
!
!  Compute change in impermeability due elasto-optic effect
          deltaB_Elasto = deltaB_Elasto + matmul(matrix_Elasto,Strain_Ext)
        end if
!
!  MODIFY ALL ENTRIES IN THE IMPERMEABILITY VECTOR
!------------------------------------------------------------------------------
!  we assume, that if former impermeability entry was zero, new entry may deviate
!  from zero after modification
        B = cmplx( real(B,DP) + deltaB_Pockels + deltaB_Kerr + deltaB_Elasto, aimag(B), DPC )
!  B matrix notation to rank 2 tensor notation
        invten = mat2ten_r2(B)
        call invertmat3(invten,matten)
!
!------------------------------------------------------------------------------
!  END OF MODIFICATION
!
        deallocate (xstmp,ystmp,zs)
      end if
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      subroutine sourcemodify(xs,ys,scalar)
      use feminterface, only: idbvip
      use femtypes
      use globalvariables, only:
      implicit none
      real(DP) :: xs, ys, scalar
      intent (in) :: xs,ys
      intent (inout) :: scalar
!
!
!  Subroutine sourcemodify modifies the power density as source for the heat
!  transfer equation. We assume uni-directional Joule heating with power den-
!  sity from heating current exported to file data_P.txt in the project dir.
!
!  ONLY TO COUPLING WORKS up to now
!
!
!  Input:
!     list               contains data from material files (refer to setstandard-
!                        values.f90 for individual position of entries)
!     xs, ys             coordinate of point with material modification
!
!  Output:
!     scalar             modified value for power density in W/m^3
!
!  Internal variables:
      integer (I4B) :: errcode, nip
      real (DP), allocatable :: xstmp(:), ystmp(:), zs(:)
      logical, save :: etmessage, etintpolstarted
!
!
!
!  initialize interpolation point
      nip = 1
      allocate (xstmp(nip),ystmp(nip),zs(nip))
      xstmp = xs
      ystmp = ys
!
!  GET FIELD DATA
!------------------------------------------------------------------------------
!  get power density field data and compute new source term
      if ( etmod .and. (xs .ge. minval(xP)) .and. (xs .le. maxval(xP)) .and. &
                       (ys .ge. minval(yP)) .and. (ys .le. maxval(yP)) ) then
        if (.not. etmessage) then
          print "(a)","INFO: Electro-thermal material modification."
          etmessage = .true.
        end if
        select case (datatypeP)
        case ('GRID')
          call bilintpol(xP,yP,zP,xs,ys,zs(1),errcode)
        case ('SCATTERED')
          if (.not. etintpolstarted) then
            call idbvip (1,size(xP),xP,yP,zP(:,1),nip,xstmp,ystmp,zs)
            etintpolstarted = .true.
          else
            call idbvip (2,size(xP),xP,yP,zP(:,1),nip,xstmp,ystmp,zs)
          end if
        end select
!  modify source term
        scalar = scalar + zs(1)
!------------------------------------------------------------------------------
!  END OF MODIFICATION
!
        deallocate (xstmp,ystmp,zs)
      end if
!
      end subroutine
!------------------------------------------------------------------------------
!
!
!
      end module mpmodule
