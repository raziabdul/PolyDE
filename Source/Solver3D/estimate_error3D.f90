subroutine estimate_error3D(v_res,v_res_norm)
      use globalvariables3D,  only: nnat, numv
      use feminterface3d,     only: residualsc3D
      use femtypes
      implicit none
      real (DP) :: v_res(:,:), v_res_norm(:)
      intent(out) :: v_res, v_res_norm
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!
! This Routine: 'estimate_error3D'
!    Pops the error estimate for each nature for:
!    Volume, Face/Surface, Edge and Node
!
! Input:
! subroutine estimate_error3D(v_err,f_err,s_err,e_err,n_err)
!    v_err      : Container for estimates per element per nature
!    f_err      : Container for estimates per face per nature
!    s_err      : Container for estimates per surface per nature
!    e_err      : Container for estimates per edge per nature
!    n_err      : Container for estimates per node per nature
!
! Output:
!    v_err      : Residual per element per nature
!    f_err      : Residual per face per nature
!    s_err      : Residual per surface per nature
!    e_err      : Residual per edge per nature
!    n_err      : Residual per node per nature
!    v_res_norm : Residual Norm per Nature
!    vrnorm_max : Highest Residual Norm among all Natures
!
!
!------------------------------------------------------------------------------
!  Internal variables:
      integer (I4B)              :: i, errcode
      real     (DP),     pointer :: resext(:,:,:)=>null()
      real     (DP), allocatable :: sumref(:)
!__
! Allocate Memory:
      allocate(sumref(nnat))
!__
! Estimate the Error of current Solution:
      call residualsc3D(.true.,.false.,.false.,errcode,v_res,sumref,resext)
      do i = 1, nnat
        if(sumref(i) .ne. 0.) then
          v_res_norm(i) = sqrt(sum(abs(v_res(:,i)))/sumref(i))
        else
          v_res_norm(i) = 0._DP
        end if
      end do
      !vrnorm_max = maxval(v_res_norm(:))
!__
! Release Memory:
      deallocate(sumref)
      deallocate(resext)
!
!_End.
end subroutine estimate_error3D
