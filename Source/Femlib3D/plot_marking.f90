      subroutine plot_marking(plot_input)
      use globalvariables3D,   only: numv
      use feminterface,        only: print_error
      use feminterface3D,      only: writedata
      use femtypes
      implicit none
      real (DP)                   :: plot_input(:)
      intent(in)                  :: plot_input
!
!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 - 2015 Institute for Micro Systems Technology,
!                              Hamburg University of Technology.
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
!------------------------------------------------------------------------------
!    $Revision: 1.1 $
!    $Date: 2015/06/17 17:04:59 $
!    $Author: juryanatzki $
!------------------------------------------------------------------------------
!
! This Routine: 'plot_marking'
!    
! Input:
!   criterion(numv,nnat)   Marking criterion values for all elements in all natures
!
! Output:
!    --
! 
! //NCHAR(int)
!------------------------------------------------------------------------------
!  Internal Variables:
          integer (I4B),    parameter :: sections=25, height=20, leftmarg=8
          integer (I4B)               :: i, j, crit_size, cell_size, occurence(0:sections)
              real (DP),  allocatable :: criterion(:), data_out(:,:)
              real (DP)               :: min_q, max_q, average_criterion, stddev, mincriterion, maxcriterion
       character(len=16)               :: string
      character(len=78)               :: str
!                      horizontal     vertical       cross          left cross     right cross    fill           half-fill
!  extented Ascii (Codepage 437)
      character (1) :: ch=achar(196), cv=achar(179), cp=achar(197), cl=achar(195), cr=achar(180), cf=achar(219), cx=achar(220)
!  standard Ascii version
!      character (1) :: ch='-',        cv='|',        cp='+',        cl='+',        cr='+',        cf='X'      , cx='x'
!__
! Allocate Memory:
      crit_size = size(plot_input)
      cell_size = INT(crit_size/25)
      
      allocate(criterion(crit_size))
      allocate(data_out(crit_size,2))
!__
!  Definitions:
      
              criterion = plot_input / REAL(maxval(plot_input))
      average_criterion = sum(criterion)/crit_size
                 stddev = sqrt( dot_product(criterion,criterion)/crit_size - average_criterion**2)
           mincriterion = minval(criterion)
           maxcriterion = maxval(criterion)
      
      min_q =  0._DP
      max_q =  1._DP
      
      occurence = 0
!__
! 1) 
      do i = 1, crit_size
        j = nint((criterion(i) - min_q)/(max_q - min_q)*sections)
        occurence(j) = occurence(j) + 1
        
        data_out(i,1:2) = (/ REAL(i) , REAL(plot_input(i)) /)
      end do
      
      str(1:leftmarg)='                      '
      str(leftmarg-1:)= ' 0% '//repeat(ch,5)//' 20% '//repeat(ch,5)         &
     &               //' 40% '//repeat(ch,5)//' 60% '//repeat(ch,5)         &
     &               //' 80% '//repeat(ch,5)//' 100%'
      str(leftmarg+55:leftmarg+55+12)='Marking'
      print*,str

      do i= height, 1, -1
        select case (i)
        case (15)
          str(leftmarg- 5:leftmarg- 2)='75%'
          str(leftmarg+ 1:leftmarg+ 9)=repeat(ch,9)
          str(leftmarg+11:leftmarg+19)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+21:leftmarg+29)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+31:leftmarg+39)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+41:leftmarg+49)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg   :leftmarg   )=cl
          str(leftmarg+10:leftmarg+10)=cp
          str(leftmarg+20:leftmarg+20)=cp
          str(leftmarg+30:leftmarg+30)=cp
          str(leftmarg+40:leftmarg+40)=cp
          str(leftmarg+50:leftmarg+50)=cr
        case (10)
          str(leftmarg- 5:leftmarg- 2)='50%'
          str(leftmarg+ 1:leftmarg+ 9)=repeat(ch,9)
          str(leftmarg+11:leftmarg+19)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+21:leftmarg+29)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+31:leftmarg+39)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+41:leftmarg+49)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg   :leftmarg   )=cl
          str(leftmarg+10:leftmarg+10)=cp
          str(leftmarg+20:leftmarg+20)=cp
          str(leftmarg+30:leftmarg+30)=cp
          str(leftmarg+40:leftmarg+40)=cp
          str(leftmarg+50:leftmarg+50)=cr
        case ( 5)
          str(leftmarg- 5:leftmarg- 2)='25%'
          str(leftmarg+ 1:leftmarg+ 9)=repeat(ch,9)
          str(leftmarg+11:leftmarg+19)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+21:leftmarg+29)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+31:leftmarg+39)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+41:leftmarg+49)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg   :leftmarg   )=cl
          str(leftmarg+10:leftmarg+10)=cp
          str(leftmarg+20:leftmarg+20)=cp
          str(leftmarg+30:leftmarg+30)=cp
          str(leftmarg+40:leftmarg+40)=cp
          str(leftmarg+50:leftmarg+50)=cr
        case default
          str(leftmarg:)='                                                  '
          str(leftmarg- 6:leftmarg- 3)='   '
          str(leftmarg- 1:leftmarg+ 2)=' '//cv//' '
          str(leftmarg+10:leftmarg+10)=cv
          str(leftmarg+20:leftmarg+20)=cv
          str(leftmarg+30:leftmarg+30)=cv
          str(leftmarg+40:leftmarg+40)=cv
          str(leftmarg+50:leftmarg+50)=cv
        end select

!  display some statistic at right margin
        select case (i)
        case (height)
            str(leftmarg+55:leftmarg+55+12)='Criterion'

        case (height- 4)
          str(leftmarg+55:leftmarg+55+12)='aver. Criterion'
        case (height- 5)
          write(str(leftmarg+55:leftmarg+55+12),'(F9.7)') average_criterion
        case (height- 8)
          str(leftmarg+55:leftmarg+55+12)='st. deviation'
        case (height- 9)
          write(str(leftmarg+55:leftmarg+55+12),'(F9.7)') stddev
        case (height-12)
          str(leftmarg+55:leftmarg+55+12)='min Criterion'
        case (height-13)
          write(str(leftmarg+55:leftmarg+55+12),'(F9.7)') minval(plot_input)
        case (height-16)
          str(leftmarg+55:leftmarg+55+12)='max Criterion'
        case (height-17)
          write(str(leftmarg+55:leftmarg+55+12),'(F12.1)') maxval(plot_input)
        case default
          str(leftmarg+55:leftmarg+55+12)='             '
        end select

!__
! 3) Draw Histogram:
!        do j=0, sections
!        !  fill: lower bound is 75% of level, i.e. i-0.25
!          if ( real(occurence(j))/real(crit_size)*100. .ge. (real(i)-0.25)) then 
!            str(2*j+leftmarg : 2*j+leftmarg)=cf
!        !  half-fill: lower bound is 25% of level, i.e. i-0.75
!          else if ( real(occurence(j))/real(crit_size)*100. .ge. (real(i)-0.75)) then 
!            str(2*j+leftmarg : 2*j+leftmarg)=cx
!          end if
!        end do
!        print*,str
!__
! 3) Draw XY-Plot:
        do j=0, sections
        !  fill: lower bound is 75% of level, i.e. i-0.25
          if ( criterion(j*cell_size+1)*100. .ge. (real(i*4)-0.25)) then 
            str(2*j+leftmarg : 2*j+leftmarg)=cf
        !  half-fill: lower bound is 25% of level, i.e. i-0.75
          else if ( criterion(j*cell_size+1)*100. .ge. (real(i*4)-0.75)) then 
            str(2*j+leftmarg : 2*j+leftmarg)=cx
          end if
        end do
        print*,str
      
      end do

      str(1:leftmarg)='                      '
      str(leftmarg-1:)= ' 0% '//repeat(ch,5)//' 20% '//repeat(ch,5)         &
     &               //' 40% '//repeat(ch,5)//' 60% '//repeat(ch,5)         &
     &               //' 80% '//repeat(ch,5)//' 100%'
      print*,str
!      print*, occurence, sum(occurence)
!__
!
      call writedata(0,'Elemarking.txt','XYPLOT','Convergence Data','NDOF','Rel. Error',data_out(:,1),data_out(:,2))
!__
!  Release Memory:
      deallocate(criterion)
      deallocate(data_out)
!__
!
      return
!_End.
end subroutine plot_marking
