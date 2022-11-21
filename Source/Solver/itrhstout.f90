      subroutine itrhstout(nnzero,iter,resvec,ndof)
      use feminterface, only: getsetting, timestamp, hi2low
      use femtypes
      implicit none
      real (DP) resvec(:)
      integer (I4B) :: nnzero, iter, ndof
      intent (in) :: nnzero, iter, resvec, ndof
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
!    $Revision: 1.7 $
!    $Date: 2006/07/07 23:45:11 $
!    $Author: r_abdul $
!
!  Inputs
!
!   ndof     number of degrees of freedom
!   nnzero   number of nonzeros (estimated)
!   iter     number of iterations
!   resvec   iterarion history
!
!
!  writes out iteration history vector to a file 
!  file name format is: yymmdd_hhMM_pp.vvv, where
!    yy  ... year
!    mm  ... month
!    dd  ... date
!    hh  ... hour
!    MM  ... minute
!    pp  ... polynomial order
!    vvv ... solver
!
!  Format inside file:
!
! 1    juliandate
! 2   polynomial order
! 3   number of elements
! 4   number of DOFs
! 5   number of nonzeros
! 6   number of iterations
! 7   [residuals vector]
!
!  local variables
      real (DP) jd
      integer (I4B) i, unitid, ios, porder
      character (len=200) path
      character (len=3) solver
      character (len=8) date
      character (len=10) time      
      character (len=12) pp
      character (len=40) filename

!   Get the project's path
      call getsetting('PROJECTPATH',path)
      call getsetting('LINSOLVERTYPE',solver)
      call getsetting('POLYORDER',porder)

      call hi2low(solver,3)

      call grglun(unitid)
!
!  read in timestamp
      jd = timestamp()

! get date info 
      call date_and_time(date,time)

! convert integer porder to character:
      write (pp,*) porder
      pp=trim(adjustl(pp))
      if (porder .lt. 10) then
         pp='0'//pp
      end if

! construct filename
      filename=date(3:8)//'_'//time(1:4)//'_'//'_'//pp(1:2)//'.'//solver

!  open the file iteration.hst in the specified path
         open (unitid,file=path(1:len_trim(path))//filename,         &
     &     form='formatted',position='REWIND',action='WRITE',iostat=ios)
         if (ios .ne. 0) then
            print*,'Error opening file: ',path(1:len_trim(path))//'ih.'//solver
            print*,'Error no: ',ios
            return
         end if

!  write time stamp to file
         write(unitid,*) jd
!  write polynomial order
         write(unitid,*) porder
!  write number of DOFS
         write(unitid,*) ndof
!  write number of NNZEROS (estimated)
         write(unitid,*) nnzero
!  write total iterations
         write(unitid,*) iter
!  write residual
         do i=1,iter
            write(unitid,*) resvec(i)
         end do

         close (unitid)

       end subroutine itrhstout
