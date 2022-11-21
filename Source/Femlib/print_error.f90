      subroutine print_error(error_level,msg,warning)
      use femtypes
      implicit none
      integer (I4B)                   :: error_level
      character (len=*),     optional :: warning
      character (len=*)               :: msg
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
!
!    $Revision: 1.2 $
!    $Date: 2015/11/11 17:05:59 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'print_error'
! Prints out an Error on the screen
!
!-------------------------------------------------------------------------------
! Input:       
!          error_level   Level of error. possible levels:
!                        1) LOW PRIORITY  (Message is displayed, no routine interrupt)
!                        2) SIGNIFICANT   (Message is displayed, no routine interrupt)
!                        3) HIGH PRIORITY (program pauses and waits for user acknowledgement)
!                        4) CRUTIAL       (program pauses and waits for user acknowledgement, user gets additional warning)
!                        5) FATAL         (program will be terminated)
!                  msg   Error message
!              warning   Optional warning
!
! Output:      
!                        Screenprint of the error message
!
!-------------------------------------------------------------------------------
! local variables
   character (len=260)   :: preamble, msg_type, level_type
!__
! Definitions:
   msg_type = 'ERROR: '
   select case (error_level)
     case(1)
       level_type = '* LOW PRIORITY '
     case(2)
       level_type = '** SIGNIFICANT '
     case(3)
       level_type = '*** HIGH PRIORITY '
     case(4)
       level_type = '**** CRUTIAL '
     case(5)
       level_type = '***** FATAL '
     case default
       level_type = '* UNKNOWN PRIORITY '
   end select
!-
! Add Warning if present, use new line (achar(10))
   if (present(warning)) then
     msg = msg(1:len_trim(msg))//achar(10)//warning(1:len_trim(warning))
   end if
   preamble = level_type(1:len_trim(level_type))//' '//msg_type(1:len_trim(msg_type))//' '
!_____
! Start:
   select case (error_level)
     case(1)
       print*,preamble(1:len_trim(preamble))//msg(1:len_trim(msg))
       
     case(2)
       print*,preamble(1:len_trim(preamble))//msg(1:len_trim(msg))
       
     case(3)
       print*,preamble(1:len_trim(preamble))//msg(1:len_trim(msg))
       pause
       
     case(4)
       print*,preamble(1:len_trim(preamble))//msg(1:len_trim(msg))
       print "(A2)" ,'--'
       print "(A33)",'Simulation might loose stability.'
       pause
       
     case(5)
       print*,preamble(1:len_trim(preamble))//msg(1:len_trim(msg))
       print "(A2)" ,'--'
       print "(A23)",'Program will be stoped.'
       pause
       stop
       
     case default
       print*,preamble(1:len_trim(preamble))//msg
       pause

   end select
!
! End.
!__
 return
end subroutine print_error