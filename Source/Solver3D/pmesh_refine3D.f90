      subroutine pmesh_refine3D(v_mark,nat,num_refined,ref_vp,pconfirm)
      use feminterface,        only: getsetting
      use globalvariables3d,   only: numv, polymaxsc, vp, vp_temp, vv
      use femtypes
      implicit none
      integer (I4B), optional              :: pconfirm(:,:)
      integer (I4B)                        :: nat, num_refined
      integer (I4B)                        :: v_mark(:)
      logical                              :: ref_vp
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
!    $Revision: 1.4 $
!    $Date: 2015/11/11 11:49:18 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!  
!  This Routine: 'padapt3D'
!  Adaptive 3D polydegree refinement. The worst maxraise % of all
!  elements should be raised (or lowered) in polynomial degree by a given factor 'v_mark' if
!  their residual is worse than maxres. Minraise % are always raised in polynomial order by 'v_mark'.
!
!-------------------------------------------------------------------------------
! local variables
!
   integer (I4B)              :: f, v, new_deg, nb, delta_p
   integer (I4B)              :: mean_poly, lowest_poly, highest_poly
   logical                    :: stopdoing
   logical      , allocatable :: todolist(:)
!__
!  
! Definitions:
   call getsetting('DELTA_POLYDEG',delta_p)
!-
   num_refined = 0
!-
   allocate(todolist(numv))
!_____
! Start:
!    If 'ref_vp' is true, perform p-Adaptation on vp, otherwise on vp_temp
   if(ref_vp) then
!__
! 1) Adjust polydegrees according to the marking of elements
     do v = 1, numv
       if ( v_mark(v).ne.0 ) then
         new_deg = vp(v,nat) + v_mark(v)
         if (new_deg .gt. polymaxsc) then
         ! No p-Adapt performed if new_deg is above polymax
           vp(v,nat) = polymaxsc
         else if (new_deg .lt. 2) then
         ! No p-Adapt performed if new_deg is below 1
           vp(v,nat) = 1
           v_mark(v) = 0
         else
         ! P-Adapt element
           if (present(pconfirm)) then
             if (new_deg.le.pconfirm(v,nat)) then
               vp(v,nat) = new_deg
             else if (vp(v,nat).gt.pconfirm(v,nat)) then
               vp(v,nat) = pconfirm(v,nat)
             end if
           else
             vp(v,nat) = new_deg
           end if
         ! Make visible, that p_add has been used by setting it to '0'
           if (v_mark(v) .gt. 0) num_refined = num_refined + 1
           v_mark(v) = 0
         end if
       end if
     end do
!__
! 2) Look if polynomial order of neighbors differ from the element's by more than 'delta_p'.
!    Change lower polyorder in corresponding element if necessary.
     do
       stopdoing = .false.
       todolist  = .false.
       do v = 1,numv
         do f = 1,4
           nb = vv(f,v)
!  Cycle if there is no neighbor or neighbor is smaller than element.
           if ((nb .lt. 0) .or. (nb .lt. v)) cycle
           if (vp(nb,nat) .lt. (vp(v,nat) - delta_p)) then
             todolist(nb) = .true.
           else if (vp(v,nat) .lt. (vp(nb,nat) - delta_p)) then
             todolist(v) = .true.
           end if
         end do
       end do
       stopdoing = .true.
!  Modify (raise) element's polynomial order to reduce the difference to the neighbour element if necessary.
       do v = 1,numv
         if (todolist(v)) then
           vp(v,nat) = vp(v,nat) + delta_p
           stopdoing = .false.
         end if
       end do
       if (stopdoing) exit
     end do
     
     mean_poly    = sum(vp(:,nat))/numv
     lowest_poly  = minval(vp(:,nat))
     highest_poly = maxval(vp(:,nat))
     print "(A16,I2)","|    Mean poly: ", mean_poly
     print "(A16,I2)","|  Lowest poly: ", lowest_poly
     print "(A16,I2)","| Highest poly: ", highest_poly
     print "(A1)"    ,'|'
   else
!__
! 1) Adjust polydegrees according to the marking of elements
         do v = 1, numv
           if ( v_mark(v).ne.0 ) then
             new_deg = vp_temp(v,nat) + v_mark(v)
             if (new_deg .gt. polymaxsc) then
             ! No p-Adapt performed if new_deg is above polymax
               vp_temp(v,nat) = polymaxsc
               v_mark(v) = 0
             else if (new_deg .lt. 2) then
             ! No p-Adapt performed if new_deg is below 1
               vp_temp(v,nat) = 1
               v_mark(v) = 0
             else
             ! P-Adapt element
               if (present(pconfirm)) then
                 if (new_deg.le.pconfirm(v,nat)) then
                   vp_temp(v,nat) = new_deg
                 else
                   v_mark(v) = 0
                 end if
               else
                 vp_temp(v,nat) = new_deg
               end if
             ! Make visible, that p_add has been used by setting it to '0'
               v_mark(v) = 0
               num_refined = num_refined + 1
             end if
           end if
         end do
!__
! 2) Look if polynomial order of neighbors differ from the element's by more than 'delta_p'.
!    Change lower polyorder in corresponding element if necessary.
         do
           stopdoing = .false.
           todolist  = .false.
           do v = 1,numv
             do f = 1,4
               nb = vv(f,v)
!  Cycle if there is no neighbor or neighbor is smaller than element.
               if ((nb .lt. 0) .or. (nb .lt. v)) cycle
               if (vp_temp(nb,nat) .lt. (vp_temp(v,nat) - delta_p)) then
                 todolist(nb) = .true.
               else if (vp_temp(v,nat) .lt. (vp_temp(nb,nat) - delta_p)) then
                 todolist(v) = .true.
               end if
             end do
           end do
           stopdoing = .true.
!  Modify (raise) element's polynomial order to reduce the difference to the neighbour element if necessary.
           do v = 1,numv
             if (todolist(v)) then
               vp_temp(v,nat) = vp_temp(v,nat) + delta_p
               stopdoing = .false.
             end if
           end do
           if (stopdoing) exit
         end do
         mean_poly    = sum(vp_temp(:,nat))/numv
         lowest_poly  = minval(vp_temp(:,nat))
         highest_poly = maxval(vp_temp(:,nat))
         print "(A16,I2)","|    Mean poly: ", mean_poly
         print "(A16,I2)","|  Lowest poly: ", lowest_poly
         print "(A16,I2)","| Highest poly: ", highest_poly
         print "(A1)"    ,'|'
   end if
!__
! Check for negative p_mark values and remove them before ending the routine
     do v = 1, numv
       if (v_mark(v).lt.0) v_mark(v) = 0
     end do
! Release Memory:
   deallocate(todolist)
!__
!  
! End.
!___
 return

end subroutine pmesh_refine3D

