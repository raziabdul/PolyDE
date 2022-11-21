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
!    $Revision:  $
!    $Date:  $
!    $Author:  $
!

module varmat_handler
    use femtypes
    implicit none

    real(DP), allocatable   :: material_map(:,:)
    integer(I4B)            :: nrecords
    logical                 :: varmat_file_read = .false.

    contains
    
    
!    Read file containing spatial material property information (e.g. doping intensity)
    
! For now one-dimensional data in 'CIS Format'
!
! TODO:
!   Better format handling + error handling
!   Add support for comments, delimiters etc.
!   Support 3D data and multiple values per point
!   Store data in a structure instead of array

subroutine varmat_deallocate()
    if (allocated(material_map)) deallocate(material_map)
end subroutine varmat_deallocate
!
!
!
!------------------
!  readmatprops
!------------------
!
!
subroutine readmatprops(filename, fileformat, ok)
    
    use femtypes
    use feminterface,      only :getsetting
    
    implicit none

    character(len=*)                  ::      filename, fileformat
    logical                           ::      ok
    
    intent(in)          ::      filename
    intent(out)         ::      ok
!------------------------------------------------------------------------------
! This Routine: 'readmatprops'
!  Reads a file providing material properties on sampled positions
!
!------------------------------------------------------------------------------
    ! Local variables
    integer(I4B)        ::      ncols, a, ios, unitid, i
    real(DP)            ::      geoscale
    character(len=200)  ::      path
    character(len=10)   ::      name1, name2

    
    call getsetting('PROJECTPATH',path)
    call getsetting('GEOMETRY_FACTOR',geoscale)
    ok = .true.
    
    call grglun(unitid)
    
    ! 1. Open file and handle possible error
    open (unitid,file=path(1:len_trim(path))//filename,            &
     &      form='formatted',position='REWIND',action='READ',iostat=ios)
    
    if (ios .ne. 0) then
        print*, '**** Error opening file: ',path(1:len_trim(path))//filename
        print*, '**** IO Error No: ',ios
        ok=.false.
        return
    end if
    
	! 2. Read Header
	! 			First line is the title which we ignore
	! 			Second line has number of records, number of coloumns and and another unknown number (a)
	!			Fourth line is name of first col
	!			Fifth line is name of second col
    
    if (fileformat .eq. 'SILVACO_DOT') then
	    read(unitid,*)
	    read(unitid, *) nrecords, ncols, a ! TODO: Change this when we know the exact format specification
        read(unitid, *) name1
	    read(unitid, *) name2
    
        allocate(material_map(nrecords,ncols))
    
        ! 3. Read data: for now we just take the 1D-coordinate from first col and the value from second col.
        !    Data is then saved in a nrecords*2 matrix material_map
	    do i=1,nrecords
		    if (ios .lt. 0) then
			    ! EOF
			    exit
		    end if
		
		    read(unitid, *, iostat=ios) material_map(i,1), material_map(i,2)
            material_map(i,1) = material_map(i,1)*geoscale
            material_map(i,2) = material_map(i,2)*10.E6_DP
        end do
        
        varmat_file_read=.true.
        close(unitid)
    else
        print*, '**** Unsupported file format: ', fileformat
        ok=.false.
        return
    end if
    
    
    
    
    !! __Write output file for testing purposes
    !call grglun(unitid)
    !open (unitid,file=path(1:len_trim(path))//'interp_test.txt', form='formatted',position='REWIND',action='WRITE',iostat=ios,status='REPLACE',recl=150)
    !
    !do n=1,numn
    !    write(unitid,*) nod(1,n), nod(2,n), nod(3,n), mat_prop(n)
    !end do
    !
    !call writeoutnodevalues(mat_prop, 1) ! VTK TEST
    !
    !close(unitid)
    !!__
   
       
    if (ok .eqv. .false.) then 
        pause
        stop
    end if
    
end subroutine readmatprops
    
!------------------
!  interpolate_1d
!------------------
! One dimensional linear interpolation
! 
! Input:
!   mapping     the given mapping (:,2) (e.g. z->val)
!   pos         the requested position for which the value will be interpolated
!
! Output:
! interp_val    the interpolated value
subroutine interpolate_1d (mapping, pos, interp_val)
    use femtypes
    implicit none
    
    real(DP), allocatable        ::  mapping(:,:)
    real(DP)                    ::  pos, interp_val
    
    intent(in)  ::  mapping, pos
    intent(out) ::  interp_val
    
    ! Local variables:
    integer(I4B)    ::  closest_point_index, interp_index1, interp_index2
    real(DP)        ::  closest_point, closest_value, y_0, y_1, v_0, v_1
        
    if (pos .gt. maxval(mapping(:,1),1)) then
    ! Add extrapolation...? for now just take last value
        interp_val = mapping(maxloc(mapping(:,1),1) ,2)
        return
    end if
    
        
        closest_point_index =   minloc(abs(mapping(:,1)-pos),1)
        closest_point       =   mapping(closest_point_index,1)
        closest_value       =   mapping(closest_point_index,2)
        
        ! Decide between which two points to interpolate
        if (closest_point .lt. pos) then
            interp_index1 = closest_point_index
            interp_index2 = closest_point_index + 1
        else
            interp_index1 = closest_point_index - 1
            interp_index2 = closest_point_index         
        end if
        
        ! Check if we are on either end of mapping
        if (interp_index1 .eq. 0) then
            y_0 = 0
            y_1 = mapping(interp_index2,1)
            v_0 = 0
            v_1 = mapping(interp_index2,2)
        elseif (interp_index2 .ge. size(mapping(:,1))) then
            y_0 = mapping(interp_index1,1)
            y_1 = mapping(size(mapping(:,1)),1)
            v_0 = mapping(interp_index1,2)
            v_1 = mapping(size(mapping(:,2)),2)
        else
            y_0 = mapping(interp_index1,1)
            y_1 = mapping(interp_index2,1)
            v_0 = mapping(interp_index1,2)
            v_1 = mapping(interp_index2,2)
        end if
        
        ! Do the actual interpolation
        if (y_1 .eq. y_0) then
              interp_val = closest_value ! Avoid dividing by zero
        else
            interp_val = v_0 + ((v_1-v_0)/(y_1-y_0))*(pos - y_0)    
        end if
end subroutine interpolate_1d
!
!------------------
!  extrap_gauss
!------------------
! A normal (gauss) distribution function (semi bell-curve) for extrapolating doping data 
subroutine extrap_gauss(maxval, rel_pos, Gw, extrap_val)
    use femtypes
    use feminterface, only : getsetting
    implicit none
    
    real(DP)    :: maxval, rel_pos, Gw, extrap_val
    
    intent(in)  :: maxval, rel_pos, Gw
    intent(out) :: extrap_val
!
! Local variables:
!
    real(DP)    :: geoscale    
    
    call getsetting('GEOMETRY_FACTOR',geoscale)

    extrap_val = maxval * exp(-((rel_pos/geoscale)**2)/(Gw/geoscale))
    
end subroutine extrap_gauss
!
!
!------------------
!  getvarmatprop
!------------------
!
!
subroutine getvarmatprop(matname,paramname,xyzs,mat_value,set)
    use femtypes
    use feminterface, only: getsetting, print_error, number2string
    
    implicit none
    
    character (len=*)   ::  matname, paramname
    real(DP)            ::  xyzs(:), mat_value
    logical             ::  set
    intent(in)          ::  matname, xyzs
    intent(out)         ::  mat_value, set
!-------------------------------------------------------------------------------
!
! This Routine: 'getvarmatprop'
!    This is a core managing routine, which obtains parameters of variational materials.
!
!------------------------------------------------------------------------------
! Local variables:
    real(DP)            :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, MARGIN
    real(DP)            :: amvcentr_x, amvcentr_y, amvcentr_z
    real(DP)            :: amvsidel_x, amvsidel_y, amvsidel_z
    real(DP)            :: x_coord, y_coord, z_coord, interp_val, x_dist, y_dist, z_dist, dist
    real(DP)            :: geoscale
    character(len=24)   :: amvshape, str
    logical             :: LESS_XMIN, GRTR_XMAX, LESS_YMIN, GRTR_YMAX, LESS_ZMIN, GRTR_ZMAX
    logical             :: CALC_MARGIN, INSIDE_AMV
    logical             :: ok
!__
! Checks and Definitions:

    set          = .false.
    CALC_MARGIN  = .false.
    INSIDE_AMV  = .true.
    
    call getsetting('GEOMETRY_FACTOR',geoscale)
    
    call getsetting('AMV_SHAPE',amvshape)
    
    call getsetting('AMV_MARGIN_SIZE',MARGIN)
    
    if (MARGIN .gt. 0._DP) MARGIN = MARGIN*geoscale
    
!- Coordinates of the point
    x_coord  = xyzs(1)
    y_coord  = xyzs(2)
    z_coord  = xyzs(3)
    
    
!- Define, whether a calculation inside the margin is necessary. Depending on the shape of AMV
!  And calculate the peak value for the margin interpolation function
    select case (amvshape)
      
      case('CUBOID')
        LESS_XMIN    = .false.
        GRTR_XMAX    = .false.
        LESS_YMIN    = .false.
        GRTR_YMAX    = .false.
        LESS_ZMIN    = .false.
        GRTR_ZMAX    = .false.

        call getsetting('AMV_CENTER_X',amvcentr_x)
        amvcentr_x = amvcentr_x*geoscale
    
        call getsetting('AMV_CENTER_Y',amvcentr_y)
        amvcentr_y = amvcentr_y*geoscale
    
        call getsetting('AMV_CENTER_Z',amvcentr_z)
        amvcentr_z = amvcentr_z*geoscale
    
        call getsetting('AMV_SIDE_X',amvsidel_x)
        amvsidel_x = amvsidel_x*geoscale
    
        call getsetting('AMV_SIDE_Y',amvsidel_y)
        amvsidel_y = amvsidel_y*geoscale
    
        call getsetting('AMV_SIDE_Z',amvsidel_z)
        amvsidel_z = amvsidel_z*geoscale
!__
! 1a) Find out the limits of the bounded amv region without considering the optional margin
        XMIN = amvcentr_x - amvsidel_x/2._DP
        XMAX = amvcentr_x + amvsidel_x/2._DP
        YMIN = amvcentr_y - amvsidel_y/2._DP
        YMAX = amvcentr_y + amvsidel_y/2._DP
        ZMIN = amvcentr_z - amvsidel_z/2._DP
        ZMAX = amvcentr_z + amvsidel_z/2._DP
!__
! 1b) Check, whether the point lies outside the mater-variational region, set logicals according to how the position coordinates exceed the bounded amv region
!     This might be usefull for some materials.
        if (x_coord.lt.XMIN) LESS_XMIN = .true.
        if (x_coord.gt.XMAX) GRTR_XMAX = .true.
        if (y_coord.lt.YMIN) LESS_YMIN = .true.
        if (y_coord.gt.YMAX) GRTR_YMAX = .true.
        if (z_coord.lt.ZMIN) LESS_ZMIN = .true.
        if (z_coord.gt.ZMAX) GRTR_ZMAX = .true.
        
        if (LESS_XMIN .or. GRTR_XMAX .or. LESS_YMIN .or. GRTR_YMAX .or. LESS_ZMIN .or. GRTR_ZMAX) then ! point lies outside the bounded mater variation area - compute an attenuation at the margin
          if (MARGIN .gt. 0._DP) CALC_MARGIN = .true.
          INSIDE_AMV = .false.
        end if

      case('SPHERE')
        ! TODO: Checks for location inside the sphere or outside need to be implemented
        
      case default
        call print_error(5,'Unknown mater Variation region shape: '//amvshape)
        
    end select

!__
! 4) Compute interpolation value for the Point. If outside the AMV region, provide the peak value for the margin interpolation Gauss function
    select case (matname)
      case('CIS-pSi-Resistor', 'CIS-pSi-Resistor2')
        if (.not.varmat_file_read) call print_error(5,'Sampled material data was not read during initialization. Please check all settings.')
        select case (paramname)
          case('N_A')
              call interpolate_1d(material_map, z_coord, interp_val)
              
          case default
            call print_error(5,'No specifications for material parameter are found for parameter: '//paramname)
            
        end select
            

      case default
        call print_error(5,'No specifications for the material are found for: '//matname)
        
    end select
        
!__
! 3) Margin interpolation
    if ( CALC_MARGIN ) then
      x_dist = 0._DP
      y_dist = 0._DP
      z_dist = 0._DP
      if ( (x_coord.gt.XMAX) .or. (x_coord.lt.XMIN) ) x_dist = MIN(abs(XMAX - x_coord), abs(XMIN - x_coord)) 
      if ( (y_coord.gt.YMAX) .or. (y_coord.lt.YMIN) ) y_dist = MIN(abs(YMAX - y_coord), abs(YMIN - y_coord))
      if ( (z_coord.gt.ZMAX) .or. (z_coord.lt.ZMIN) ) z_dist = MIN(abs(ZMAX - z_coord), abs(ZMIN - z_coord))
      
      dist = sqrt(x_dist**2 + y_dist**2 + z_dist**2)
      
      call extrap_gauss(interp_val, dist, MARGIN, mat_value)
      set = .true.
    else
      mat_value = interp_val
      set = .true.
    end if
 
end subroutine 


end module varmat_handler