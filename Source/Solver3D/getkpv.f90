subroutine getkpv(mod_res,keyPointValues)
      use globalvariables3D,    only: numv, numn, vv, vn, en, nod, sn, numdom
      use feminterface3d,       only: writedata, writeoutelementvalues, writeoutnodevalues, reversemap
      use feminterface,         only: qsortindex, destroyarrptr
      use femtypes
      implicit none
      
      real     (DP)             :: mod_res(:), keyPointValues(:)
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
! This Routine: 'getkpv'
! 
! Computation of a singularity measure called 'keypoint vaue' for each node in the mesh
!
!
! 
!
!-------------------------------------------------------------------------------
! Input:       
!           v_res(numv)            Error estimate for each Element
!
! Output:   keyPointValues(numn)   
!         
!
!-------------------------------------------------------------------------------
! local variables:
!
    integer (I4B)           :: n_0, n, v, v_nb, num_nb
    real (DP)               :: kp
!__
! Definitions:
    keyPointValues(:) = 0
    num_nb = size(vv,1)
    kp = 0
    
!__
! Start
  
!  do n = 1, num_ns
!    if(size(ns(n)%d) >= 1)then
!      v_surr = size(nv(n)%d)
!      do v=1, v_surr
!        vol = nv(n)%d(v)
!        do n_nb = 1, num_nb
!          v_nb = vv(n_nb,vol)
!          if(v_nb > 0) then 
!            do n_t = 1, num_nb
!              if(n /= vn(1,v_nb) .and. n /= vn(2,v_nb) .and. n /= vn(3,v_nb) .and. n /= vn(4,v_nb)) then
!                kp = mod_res(vol) - mod_res(v_nb)   
!                if(kp > keyPointValues(n)) keyPointValues(n) = kp
!                exit
!              end if
!            end do
!          end if           
!        end do
!      end do
!    end if
!  end do

  do v = 1, numv
    do n = 1, num_nb
      n_0  = vn(n,v)
      v_nb = vv(n,v)
      if(v_nb > 0) then 
        kp = mod_res(v) - mod_res(v_nb)
        if(kp > keyPointValues(n_0)) keyPointValues(n_0) = kp
      end if           
    end do
  end do
! End.
!___
 return
end subroutine getkpv