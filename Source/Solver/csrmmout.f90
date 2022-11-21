    
!subroutine to output a csr matrix to MatrixMarket format
    
    subroutine csrmmout(ia,ja,csr,adaptstep)
    use feminterface, only: mmwrite, getsetting
    use globalvariables, only: ndof
    use femtypes
    implicit none
    integer (I4B), pointer :: ia(:) ,ja(:)
    complex (DPC), pointer :: csr(:)
    integer (I4B), optional :: adaptstep
    intent (in) :: ia, ja, csr, adaptstep
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
!    $Revision: 1.1 $
!    $Date: 2010/10/25 11:57:22 $
!    $Author: chvokas $
!
!-------------------------------------------------------------------------------
!
!  local variables:
    integer col_cnt, i, row_i_entries, j

!Adapted from mmwrite     
    ! Matrix Market inputs
    integer (I4B) ounit, ios
    character rep*10 
    character field*7 
    character symm*19
    character ofile*32
    integer nrows,ncols,nnz

    integer indx(size(ja)) ! row indices (automatic shape from ja)
    !integer jndx(10) ! column indices (assume shape from ja)
    integer , allocatable :: ival(:) 
    double precision , allocatable :: rval(:) 
    !complex (DPC), allocatable :: cval(:)     
    character (len=32) filename
    character (len=10) adap
    character (len=16) :: adapttype

    nnz = ia(ndof+1)-1    ! number of nonzeros 
    nrows = ndof          ! number of rows
    ncols = ndof          ! number of columns
    
    ! find explicit array of row indices 
    col_cnt=1 !column counter
    do i=1,nrows
        row_i_entries=ia(i+1)-ia(i)
        do j=1,row_i_entries
            ! store row indices
            indx(col_cnt)=i
            col_cnt = col_cnt+1
        enddo 
    enddo
    
!  check which adaptation is running
    call getsetting('ADAPTION_TYPE',adapttype)
    adapttype = adapttype(1:len_trim(adapttype))
    
! create filenames like, A_.mtx, A_h[1,n].mtx, A_p[1,n].mtx, A_hp[1,n].mtx ...
    filename = 'A_'
    if (present (adaptstep)) then
        print*, '----------= adaptation =-------', adaptstep
        write (adap,'(i2)') adaptstep ! <- convert integer to string char
        !print*,'ad=',adap(1:len_trim(adap))
    else 
        adap='0'  !first pre-solution or no adaptation
    end if
    select case (adapttype)
        case ('NO_ADAPT')
                filename = filename(1:len_trim(filename))//'no_adapt'//adjustl(adap)
        case ('H_ADAPT')
                filename = filename(1:len_trim(filename))//'h'//adjustl(adap)
        case ('P_ADAPT')
                filename = filename(1:len_trim(filename))//'p'//adjustl(adap)
        case ('HP_ADAPT')
                filename = filename(1:len_trim(filename))//'hp'//adjustl(adap)
    end select    
    
    !print*,filename  
    
    ! append extension
    filename = filename(1:len_trim(filename))//'.mtx'

    ! matrix attributes (will be kept on the .mtx file)
    rep = 'coordinate'
    !field = 'complex'
    field = 'pattern'  ! no values just the structure
    symm = 'general'
    !ofile = 'A.mtx'
    ofile = filename
    ! get an available unit
    call grglun(ounit)
    
    ! open a file and check  
    open(ounit,file=ofile,position='REWIND',action='WRITE',iostat=ios)
    if (ios .ne. 0) then
        print*, 'Problem opening Matrix Market file ',ofile
        return
    end if
   
    print *,'+ Writing matrix into Matrix Market format.'
    print *,'+ Writing header and data to file: ',ofile
    ! mmwrite was changed to take crs directly, ival and rval are unused.
    call mmwrite(ounit,rep,field,symm,nrows,ncols,nnz,indx,ja,ival,rval,csr)
    
    ! release file
    close (ounit)

    return 
    end subroutine csrmmout
