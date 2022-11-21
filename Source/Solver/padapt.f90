      subroutine padapt(ext,int,epsgl)
      use feminterface, only: preassemb, assemblyandsolve, getsetting, &
                              qsortindex, residual, zeit
      use femtypes
      use globalvariables, only: n, ndof, ep, polymax, en, nnat, physics
      implicit none
      logical, intent(in) :: ext, int
      real (DP), intent(out) :: epsgl
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
!    $Revision: 1.21 $
!    $Date: 2015/04/01 11:02:10 $
!    $Author: juryanatzki $
!
!-------------------------------------------------------------------------------
!
!  This subroutine does the adaptation (only p-adaptation atm). It assigns a starting
!  value for the polynomial degree, which can be set in the FEMsettings.txt. Then
!  it calls the assemblyandsolve subroutine to compute a solution for the chosen
!  polynomial degree. After this is done it computes the residual for each element
!  and compares it with a given value (from FEMsettings). The worst maxraise % of all
!  elements should be raised in polynomial degree by one, if their residual is worse
!  than maxres. minraise % are always raised in polynomial order by one.
!
!  Output:
!            epsgl      reached error by linear solver
!
!  internal variables:
      integer (I4B) :: i, j, k, m, inat
      integer (I4B) :: adaptsteps, deltap, errcode, nb, ndofold
!  mdeltap is the multinature difference for the polynomial order
!  of the same element
!  mmaxp (multi-maxp) is used for the maximum poly order
!  of one element in all natures
      integer (I4B) :: mdeltap, mmaxp
      integer (I4B), pointer :: oldep(:,:)
      integer (I4B), allocatable :: indx(:)
      real (DP) :: adapterror, error(nnat), minraise, maxraise, maxres(nnat), resgl, sumres(nnat), sumref(nnat)
!  merror is the multinature error, the maxval of error(nnat)
      real (DP) :: merror
      real(DP), allocatable :: res(:,:)
      logical :: stopdoing
      logical, allocatable  :: todolist(:)
!  store a poly-order histogram      
      integer (I4B) :: poly_histogram(polymax,nnat)
!
!  Allocate residual vector
      allocate(res(n,nnat))
      call residual(ext,int,.false.,errcode,res,sumres,sumref)
!  return relative error of the solution in %
!  sumres is the sum of the non-normalized element residuals
      call getsetting('PHYSICS_MODE',physics)
      select case(physics)
        case('FLUIDINCOMPR')
          error(1:2) = 100._DP * sqrt(sumres/sumref)
          error(3) = 0._DP
        case('FLUIDELMASSCON')
          error(1:2) = 100._DP * sqrt(sumres/sumref)
          error(3) = 0._DP
          error(4:5) = 100._DP * sqrt(sumres/sumref)
        case default
          error = 100._DP * sqrt(sumres/sumref)
      end select
      merror = maxval(error)
!
!  Get maximum value for residual and number of adaptation steps from settings and
!  compare with maximum residual obtained from solution. Start p-adaptation if
!  this is higher and number of adaptation steps is greater than 0.
      call getsetting('ADAPT_STEPS',adaptsteps)
!  For multi-nature problems the same adaptation error is used in all natures
      call getsetting('ADAPTION_ERROR',adapterror)

      if ((merror .gt. adapterror) .and. (adaptsteps .gt. 0)) then
        print "(a)","Residual bound not reached --> Starting p-adaptation."
        allocate (indx(n),oldep(n,nnat),todolist(n))
        indx = 0
!  deltap is needed for setting the maximum difference of polynomial degree
!  between neighboring elements during the neighbor's polynomial check.
!  meaningful values are whole numbers between 1 and 9
        deltap = 1
!  the multinature maximum polynomial order difference of one element
        mdeltap = 2
!  values for minimum and maximum numbers of elements raised
        minraise = 0.05_DP
        maxraise = 0.3_DP
!
!  loop for maximum amount of adaptations
!  For multi-nature problems the same number of adaptation steps is used in all natures
        do i = 1, adaptsteps
!  value for largest local error indicator maxres (ADAPTION_ERROR divided by
!  number of elements
          ndofold = ndof
          print "(a,i3)","p-adaptation step: ",i
          oldep = ep
!  Loop over Natures
          do inat = 1, nnat
!  Skip adaptation for pressure in fluidic problems
            select case(physics)
              case('FLUIDINCOMPR','FLUIDELMASSCON')
                select case(inat)
                  case (3)
                    cycle
                  case default
                    ! Continue, do not cycle
                end select
              case default
                ! Continue, do not cycle
            end select
!
            print "(a,i3)","Nature: ",inat
!  Sort residual vector. biggest values first. indx points to entry in res.
            call qsortindex(res(:,inat),indx,n)
            maxres(inat) = (adapterror/(100._DP*n))**2 * sumref(inat)
!
!  Always raise polyorder for minraise % worst elements.
            do j = 1, ceiling(minraise*n + 1)
              if (ep(indx(j),inat) .eq. polymax) then
                cycle
              else
                ep(indx(j),inat) = ep(indx(j),inat)+1
              end if
            end do
!  Raise polyorder for maxraise % worst elements, if residual worse than maxres.
            do j = ceiling(minraise*n + 2), ceiling(maxraise*n + 1)
              if (res(indx(j),inat) .gt. maxres(inat)) then
                if (ep(indx(j),inat) .eq. polymax) then
                  cycle
                else
                  ep(indx(j),inat) = ep(indx(j),inat)+1
                end if
              end if
            end do
!
!  Look if polynomial order of neighbors differ from the element's more than 1.
!  Change lower polyorder in corresponding element if necessary.
!  loop over all elements
            do
              stopdoing = .false.
              todolist  = .false.
              do k = 1,n
                do m = 1,3
                  nb = en(m,k)
!  Cycle if there is no neighbor or neighbor is smaller than element.
                  if ((nb .eq. 0) .or. (nb .gt. k)) cycle
                  if (abs(ep(nb,inat)-ep(k,inat)) .gt. deltap) then
                    if (ep(nb,inat)-ep(k,inat) .lt. (-deltap)) then
                      todolist(nb) = .true.
                    else if (ep(k,inat)-ep(nb,inat) .lt. (-deltap)) then
                      todolist(k) = .true.
                    end if
                  end if
                end do
              end do
              stopdoing = .true.
!  Modify element's polynomial order if necessary.
              do k = 1,n
                if (todolist(k)) then
                  ep(k,inat) = ep(k,inat) + deltap
                  stopdoing = .false.
                end if
              end do
              if (stopdoing) exit
            end do
!  END OF NATURES LOOP:
          end do
!
!  Look if polynomial orders of the same element in
!  different natures differ more than 2.
!  Change lower polyorder in corresponding nature if necessary.
!  loop over all elements
          select case(physics)
            case('FLUIDINCOMPR','FLUIDELMASSCON')
              do k = 1,n
                ep(k,3) = minval(ep(k,1:2))-1
              end do
            case default
              do k = 1,n
                mmaxp = maxval(ep(k,:))
                do inat = 1,nnat
                  ep(k,inat) = max(ep(k,inat),(mmaxp-mdeltap))
                end do
              end do
          end select
!
          print "(a,i2,a,i2,a,i2)","mean poly: ",(sum(ep(:,:))/(n*nnat))," lowest poly: "&
                                  ,minval(ep(:,:))," highest poly: ",maxval(ep(:,:))
!
!  Make histogram of polyorder x order count
          poly_histogram = 0
          do inat = 1,nnat
              do j = 1,n
                  poly_histogram(ep(j,inat),inat)=poly_histogram(ep(j,inat),inat)+1
              end do
          end do
        
          print*,'+ Polynomial order distribution (from order 1 to 20):'
          write(*,'(5i10)')((((poly_histogram(j*5+k,inat)),k=1,5),j=0,3),inat=1,nnat)                                  
!
!  Solve and calculate residual for new polynomial orders
          call preassemb
!          call preassemb(oldep)
          call assemblyandsolve(epsgl,resgl,i)
          call residual(ext,int,.false.,errcode,res,sumres,sumref)
!  return relative error of the solution in %
          select case(physics)
            case('FLUIDINCOMPR')
              error(1:2) = 100._DP * sqrt(sumres/sumref)
              error(3) = 0._DP
            case('FLUIDELMASSCON')
              error(1:2) = 100._DP * sqrt(sumres/sumref)
              error(3) = 0._DP
              error(4:5) = 100._DP * sqrt(sumres/sumref)
            case default
              error = 100._DP * sqrt(sumres/sumref)
          end select
          merror = maxval(error)
!  Exit adaptation loop if maxres is reached for every element
          if (merror .le. adapterror) then
            print "(a,i3,a)","Reached stopping criterion for p-adaptation after ",i," steps."
            exit
!  Exit adaptation loop if there are no new degrees of freedom
          else if (ndof .eq. ndofold) then
            print "(a,i3,a)","Maximum amount of p-refinement reached after ",i," steps."
            exit
          end if
!  END OF ADAPTSTEPS LOOP:
        end do
!
        deallocate (indx,oldep,todolist)
        call zeit('doing p-adaptation')
        if (i .ge. adaptsteps) then
          print "(a)","Maximum number of adaptation steps used."
          select case(physics)
!  Print statements to be reviewed
            case('FLUIDINCOMPR')
              print "(a,g10.3,a,g10.3)","desired maximum residual: ",maxval(maxres(1:2))," reached: ",maxval(res(:,1:2))
            case('FLUIDELMASSCON')
              print "(a,g10.3,a,g10.3)","desired maximum residual: ",maxval(maxres(1:2))," reached: ",maxval(res(:,1:2))
            case default
              print "(a,g10.3,a,g10.3)","desired maximum residual: ",maxval(maxres)," reached: ",maxval(res)
          end select
        end if
      else
        print "(a)","Already reached residual bound. No p-adaptation necessary."
      end if
!
      deallocate(res)
      end subroutine padapt
