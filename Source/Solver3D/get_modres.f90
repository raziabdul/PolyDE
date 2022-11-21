subroutine get_modres(v_res,mod_res,vol_poly)
      use globalvariables3D,    only: numv, vv , vn, nod, sn, polymaxsc, dom
      use feminterface3d,       only: reversemap, tetvolume
      use feminterface,         only: destroyarrptr
      use femtypes
      implicit none

      real     (DP)              :: v_res(:), mod_res(:)
      integer  (I4B)             :: vol_poly(:)
      intent(in) :: v_res, vol_poly
      intent(out) :: mod_res
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
! This Routine: 'get_modres'
!
! Preprocess the residual values in such a way that they become more suitable for singularity indication
!
! 
!
!-------------------------------------------------------------------------------
! Input:       
!           v_res(numv, nnat)            Error estimate for each Element and each Nature
!
! Output:   mod_res(numv, nnat)          Modified error estimate for each Element and each Nature
!         
!
!-------------------------------------------------------------------------------
! local variables:
!
      integer (I4B)               :: n, f, v, i, v_nb, num_nb, num_ns, nmin, nmax, poly_mid
      real (DP), allocatable      :: element_vol(:), element_polyfac(:)
      integer (I4B), allocatable  :: element_fac(:), element_val(:)
      real (DP)                   :: vert(3,4), rising_edge_poly, falling_edge_poly
      type(ARRPTRI), pointer      :: ns(:)=>null()
      logical                     :: ok

!__
! Definitions:
      call reversemap(.true.,sn,ns,nmin,nmax)

      allocate(element_vol(numv), element_fac(numv), element_val(numv), element_polyfac(numv))

      num_nb = size(vv,1)
      num_ns = size(ns)
      poly_mid = anint(maxval(vol_poly)/2.)
      element_fac(:) = 7
      element_val(:) = 0
      rising_edge_poly  =   1*(real(maxval(vol_poly),DP)/real(polymaxsc,DP))
      falling_edge_poly = 100*(real(maxval(vol_poly),DP)/real(polymaxsc,DP))

!__
! Start

! 1) Calculate the volume of the elements
      do v=1, numv
        vert(1:3,1:4) = nod(1:3,vn(1:4,v))
        element_vol(v) = abs(tetvolume(vert))
      end do  

! 2) Weight-Calculation
      do v=1, numv
        do f = 1, num_nb
          v_nb = vv(f,v)
          if (v_nb .gt. 0) then
            if (dom(v) .eq. dom(v_nb)) then
              element_fac(v) = element_fac(v) - 1
            end if
          end if
        end do
        if (element_fac(v) .eq. 3)then
          do n = 1, size(ns)
            do i = 1, num_nb
              if ((vn(i,v) .eq. n) .and. (size(ns(n)%d) .ge. 1))then
                element_val(v) = element_val(v) + 1
              end if
            end do
          end do
          element_fac(v) = element_val(v)
        end if
        if (element_fac(v) .gt. 0) element_fac(v) = (2**(element_fac(v)-1))
      end do

! 3) Calculate the polynomial degree factor
      do v = 1, numv
        element_polyfac(v) = 1._DP/(exp(-(real(vol_poly(v)-poly_mid,DP)/real( rising_edge_poly,DP)))     &
     &                       +      exp(  real(vol_poly(v)-poly_mid,DP)/real(falling_edge_poly,DP)))

! 4) Calculate the new residual
        mod_res(v) = (v_res(v)*(element_fac(v)/32._DP)*element_polyfac(v))/element_vol(v)
      end do


      deallocate(element_vol, element_fac,element_val, element_polyfac)
      ok = destroyarrptr(ns)
! End.
!___
      return
end subroutine get_modres