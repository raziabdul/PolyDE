subroutine get_pconfirm(sing_val,max_vp,pconfirm)
      use globalvariables3D,    only: numv, polymaxsc
      use femtypes
      implicit none
      real     (DP)                :: sing_val(:)
      integer (I4B)                :: max_vp, pconfirm(:)
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
!    $ Revision: 0.1 $
!    $ Date: 2015/02/12 11:11:11 $
!    $ Author: juryanatzki $
!
!-------------------------------------------------------------------------------
!
! This Routine: 'get_pconfirm'
!
! ...
!
! 
!
!-------------------------------------------------------------------------------
! Input:       
!                    sing_val            ...
!                      max_vp            ...
!
! Output:   
!         pconfirm(numv,nnat)            ...
!
!-------------------------------------------------------------------------------
! local variables:
!
   integer (I4B)                :: c, v, num_classes
   real     (DP),   allocatable :: classes(:)
   real     (DP)                :: max_sing_val
!__
! Definitions:
   pconfirm(:)  = 0
   num_classes  = max_vp + 1
   max_sing_val = maxval(sing_val(:))
   allocate(classes(num_classes))
   
   do c = 1, num_classes
     classes(c) = CEILING(REAL(c)/REAL(num_classes)*max_sing_val)
   end do
!__
! Start:
! 1) 
   do v = 1, numv
     do c = 1, num_classes
       if (sing_val(v).le.classes(c)) then
         pconfirm(v) = num_classes + 1 - c
         exit
       end if
     end do
     if (pconfirm(v).eq.0) print*,'WARNING: pconfirm of element ', v,'has an error!'
   end do
!__
!  Release Memory:
   deallocate(classes)
! End.
!___
 return
end subroutine get_pconfirm