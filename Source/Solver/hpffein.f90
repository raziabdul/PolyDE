      subroutine ffeinhp(ellist,zahle,zahlk,ende,merror,maxelemt,ext,int)
      use feminterface, only: getsetting, residual, setze, qsortindex
      use globalvariables, only: ep, n, nnat, gzz, polymax, e, zki, kzi, physics
      use femtypes
      implicit none
      integer (I4B) ellist(:), zahle, zahlk, maxelemt
      real (DP) merror
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int
      intent (out) :: ende, merror
      intent (inout) :: ellist, zahle, zahlk
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
!    $Revision: 1.12 $
!    $Date: 2015/04/01 11:01:30 $
!    $Author: juryanatzki $
!
!-------------------------------------------------------------------------------
!
!  Fehlerauswertung und Verfeinerung der Elemente mit den groessten Fehlern
!
!  Eingabe:
!            n        Anzahl der Elemente
!            p        Anzahl der Knoten bzw. Matrixzeilen
!            kdim     vereinbarte Groesse für die Anzahl der Knoten
!            e        Elementinformation, Knoten der Elemente
!            en       Nacbarschaftsinformation der Elemente
!            geb      Zuordnung der Elemente zu den Gebieten
!            xn,yn    Knotenkoordinaten
!            kzi      Knoten-Zweig-Information; Zweige zu den Knoten
!            kzrb     Zweige deren Randbedingung fuer den Knoten gueltig sind
!            zrb      Randbedingungen der Zweige
!            x        Loesungsvektor
!            matzif   Materialziffern der Gebiete
!            omega    Kreisfrequenz
!            maxelemt maximum number of new elements
!
!  Ausgabe:
!            ende     =.true. wenn maximale Zahl der Knoten oder Elemente erreicht ist
!            merror    relative multinature error from a posteriori error estimation
!
!  Eingabe/ Ausgabe: 
!            ellist   Liste der zu verfeinernden Elemente und Verfeinerungsart
!            zahle    Anzahl der neu entstehenden Elemente
!            zahlk    Anzahl der neu entstehenden Knoten
!
      integer (I4B) :: i, j, inat, indx(n), keypoint, errcode
      real (DP) res(n,nnat), sumres(nnat), sumref(nnat), error(nnat)
      real (DP) :: adapterror, minraise, maxraise, maxres(nnat), sigma
      logical :: hrefine
!
!
      call residual(ext,int,.false.,errcode,res,sumres,sumref)
!
      if (zahle .ge. maxelemt) then
        ende=.true.
        return
      end if
!  return relative error of the solution in %
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
!  values for minimum and maximum numbers of elements raised. These values may
!  not be 0!
      maxraise = 0.5_DP
      minraise = maxraise*0.1_DP
!
!  Sigma is a parameter to adjust the error bound for each element. One element
!  should have a maximum residual that is less or equal to sigma % of maxres.
      sigma   = 0.75_DP
!
      call getsetting('ADAPTION_ERROR',adapterror)
!  value for largest local error indicator maxres (ADAPTION_ERROR divided by
!  number of elements, maxres(nnat)
      maxres = (adapterror/(100._DP*n))**2 * sumref
!
      do inat=1,nnat
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
!  Sort residual vector. biggest values first. indx points to entry in res.
        call qsortindex(res(:,inat),indx,n)
!  do minraise elements
        do j = 1, ceiling(minraise*n)
!  initialize hrefine variable --> no h-refinement if hrefine = .false.
          hrefine = .false.
          if (res(indx(j),inat) .le. sigma*maxres(inat)) then
            exit
          else if ( (kzi(e(1,indx(j))).lt.0) .or. (kzi(e(2,indx(j))).lt.0)    &
               .or. (kzi(e(3,indx(j))).lt.0) .or. (ep(indx(j),inat) .eq. polymax)  ) then
!  if element has a keypoint or maximum polydegree is reached --> mark element
!  for h-refinement
!
!  only elements with polymax or keypoints connected to straight lines should be refined
!  search for keypoint index of node
          keypoint = 0 !initialize safe
            do i=1,3
              if (kzi(e(1,indx(j))).lt.0) keypoint = abs(kzi(e(1,indx(j))))
              if (kzi(e(2,indx(j))).lt.0) keypoint = abs(kzi(e(2,indx(j))))
              if (kzi(e(3,indx(j))).lt.0) keypoint = abs(kzi(e(3,indx(j))))
            end do
!  search branch connected to keypoint as start- or end-point of the branch.
!  set hrefine to .true. if connected to a straight line
            do i=1,gzz
              if ( ((zki(1,i) .eq. keypoint) .or. (zki(2,i) .eq. keypoint)) .and. (zki(3,i) .eq. 0)) then
                hrefine = .true.
                exit
              end if
            end do
!  if element already reached polymax, always do h-refinement
            if (ep(indx(j),inat) .eq. polymax) hrefine = .true.
!  mark element for h-refinement and decrease polynomial order if hrefine is .true.
!  do p-refinement if hrefine is .false.
            if (hrefine) then
              call setze(ellist,zahle,zahlk,indx(j))
              ep(indx(j),inat) = ceiling(ep(indx(j),inat)/2._DP)
            else
              ep(indx(j),inat) = ep(indx(j),inat)+1
            end if
          else if (ep(indx(j),inat) .lt. polymax) then
!  if element has polydegree smaller than the maximum --> do p-refinement
            ep(indx(j),inat) = ep(indx(j),inat)+1
          end if
        end do
!  Raise polyorder for maxraise % worst elements.
        do j = ceiling(minraise*n + 1), ceiling(maxraise*n)
          if (res(indx(j),inat) .gt. sigma*maxres(inat)) then
!  If error bigger than given tolerance, do p-refinement for element
            if (ep(indx(j),inat) .eq. polymax) then
!  If maximum polydegree --> do h-refinement
              call setze(ellist,zahle,zahlk,indx(j))
            else
!  If not maximum polydegree --> do p-refinement
              ep(indx(j),inat) = ep(indx(j),inat)+1
            end if
          else
!  If error smaller or equal to given tolereance --> exit refinement
            exit
          end if
        end do
      end do
      return
      end subroutine ffeinhp



      subroutine ffeinhp_top5(ellist,zahle,zahlk,ende,merror,maxelemt,ext,int)
      use feminterface, only: getsetting, residual, setze, qsortindex
      use globalvariables, only: ep, n, nnat, polymax, physics
      use femtypes
      implicit none
      integer (I4B) ellist(:), zahle, zahlk, maxelemt
      real (DP) merror
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int
      intent (out) :: ende, merror 
      intent (inout) :: ellist, zahle, zahlk
!
!  Fehlerauswertung und Verfeinerung der Elemente mit den groessten Fehlern
!
!  Eingabe:
!            n        Anzahl der Elemente
!            p        Anzahl der Knoten bzw. Matrixzeilen
!            kdim     vereinbarte Groesse für die Anzahl der Knoten
!            e        Elementinformation, Knoten der Elemente
!            en       Nacbarschaftsinformation der Elemente
!            geb      Zuordnung der Elemente zu den Gebieten 
!            xn,yn    Knotenkoordinaten
!            kzi      Knoten-Zweig-Information; Zweige zu den Knoten
!            kzrb     Zweige deren Randbedingung fuer den Knoten gueltig sind
!            zrb      Randbedingungen der Zweige
!            x        Loesungsvektor
!            matzif   Materialziffern der Gebiete
!            omega    Kreisfrequenz
!            maxelemt maximum number of new elements
!
!  Ausgabe:
!            ende     =.true. wenn maximale Zahl der Knoten oder Elemente erreicht ist
!            error    relative error from a posteriori error estimation
!
!  Eingabe/ Ausgabe: 
!            ellist   Liste der zu verfeinernden Elemente und Verfeinerungsart
!            zahle    Anzahl der neu entstehenden Elemente
!            zahlk    Anzahl der neu entstehenden Knoten
!
      integer (I4B) :: j, inat, indx(n), errcode
      real (DP) res(n,nnat), sumres(nnat), sumref(nnat), estimerror(n), error(nnat)
      real (DP) :: adapterror, minraise, maxraise, maxres(nnat), sigma
!
      call residual(ext,int,.false.,errcode,res,sumres,sumref)
!
!  Sort residual vector. biggest values first. indx points to entry in res.
      if (zahle .ge. maxelemt) then
        ende=.true.
        return
      end if
!  return relative error of the solution in %
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
!  values for minimum and maximum numbers of elements raised. These values may
!  not be 0!
      maxraise = 0.5_DP
      minraise = maxraise*0.1_DP
!
!  Sigma is a parameter to adjust the error bound for each element. One element
!  should have a maximum residual that is less or equal to sigma % of maxres.
      sigma   = 0.75_DP
!
      call getsetting('ADAPTION_ERROR',adapterror)
!  value for largest local error indicator maxres (ADAPTION_ERROR divided by
!  number of elements
      maxres = (adapterror/(100._DP*n))**2 * sumref
!
      do inat=1,nnat
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
!  Sort residual vector. biggest values first. indx points to entry in res.
        call qsortindex(res(:,inat),indx,n)
!  do minraise elements
        do j = 1, ceiling(minraise*n)
          if (res(indx(j),inat) .le. sigma*maxres(inat)) then
            exit
          else
            call setze(ellist,zahle,zahlk,indx(j))
            if (ep(indx(j),inat) .ge. 2) then
              ep(indx(j),inat) = ep(indx(j),inat) - 1
            end if
          end if
        end do
!  Raise polyorder for maxraise % worst elements.
        do j = ceiling(minraise*n + 1), ceiling(maxraise*n)
          if (res(indx(j),inat) .gt. sigma*maxres(inat)) then
!  If error bigger than given tolerance, do p-refinement for element
            if (ep(indx(j),inat) .eq. polymax) then
!  If maximum polydegree --> do h-refinement
              call setze(ellist,zahle,zahlk,indx(j))
            else
!  If not maximum polydegree --> do p-refinement
              ep(indx(j),inat) = ep(indx(j),inat)+1
            end if
          else
!  If error smaller or equal to given tolereance --> exit refinement
            exit
          end if
        end do
      end do
      return
      end subroutine ffeinhp_top5



      subroutine ffeinhp_KPphas(ellist,zahle,zahlk,ende,merror,maxelemt,ext,int)
      use feminterface, only: getsetting, residual, setze, qsortindex, pdecoeff
      use globalvariables, only: e, ep, kzi, polymax, xn, yn, n, nnat, physics
      use femtypes
      implicit none
      integer (I4B) ellist(:), zahle, zahlk, maxelemt
      real (DP) merror
      logical :: ende, ext, int
      intent (in) :: maxelemt, ext, int
      intent (out) :: ende, merror 
      intent (inout) :: ellist, zahle, zahlk
!
!  Error evaluation and refinement of mesh in elements with large error
!  This version uses h-refinement for singularities at keypoints, 
!  to reduce phase-lag and numerical stabilibity of the linear solver
!
!  Input:
!            maxelemt maximum number of new elements
!
!  Output:
!            ende     =.true. if the maximum number of elements was reached
!            error    relative error from a posteriori error estimation
!
!  In-/ Output: 
!            ellist   liste of elements to be refined and refinement type
!            zahle    number of new elements 
!            zahlk    number of new mesh nodes
!
      integer (I4B) :: j, inat, indx(n), errcode, nlarge, nmedium
! local copy of ep is needed for non-linear problems due to "call pdecoeff"
      integer (I4B), allocatable :: eplocal(:,:)
      real (DP) :: res(n,nnat), sumres(nnat), sumref(nnat),  error(nnat)
      real (DP) :: sumsqrtres(nnat), psumres, msize2, k2, kh2
      real (DP) :: xs, ys, xelem(3), yelem(3)
      complex (DPC) :: nu(2,2,nnat,nnat), gamma(2,nnat,nnat), alpha(nnat,nnat)
      complex (DPC) :: beta(2,nnat,nnat), f1(nnat), f2(2,nnat)
! algorithm parameters
! large_errcont   error fraction of elements in class (A)               usually around 0.25
      real (DP), parameter :: large_errcont=0.25_DP
! medium_errcont  error fraction of elements in class (B)               usually in [0.5...0.7]
      real (DP), parameter :: medium_errcont=0.6_DP
! marigial_errcont  error fraction of elements in class (B)             usually in [0...0.02]
      real (DP), parameter :: marigial_errcont=0.01_DP
! kappa  a safety factor for the beginning of superexponential zone     usually 3 > kappa > 1
      real (DP), parameter :: kappa=2._DP
! khmax  maximum of wavenumber x meshsize for numerical stability of CG usually around 2.
      real (DP), parameter :: khmax=1.8_DP

!  the algorithm distinguishes between four classes (A) - (D)
!  (A)   elements with LARGE error contribution
!  (B)   elements with MEDIUM error contribution
!  (C)   elements with SMALL error contribution
!  (D)   elements with MARIGINAL error contribution
!
!  we use 
!  -  h-refinement on class (A) if the element has a keypoint
!                  intends to control error at singularities
!  -  h-refinement on classes (A) (B) (C) and (D) if mesh size is to large 
!                  intends to control phase-lag and stability of CG matrix solvers
!  -  h-refinement on classes (A) and (B) if maximum p-level is reached (p=pmax)
!  -  p-refinement on classes (A) and (B) if possible (p<pmax)
!                  intends to reduce error at exponential rate
!  -  p-UNrefinement on class (D) if possible (p>1)
!                  intends to reduce number of degrees of freedom
!
!  classe (A) - (D) are treated in three seperated subsequent loops
!
!
!  compute residuals (error estimator)
      call residual(ext,int,.false.,errcode,res,sumres,sumref)
!  return relative error of the solution in %
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
!  element-wise residual are the squared 
      res=sqrt(res)
      sumsqrtres = sum(res(:,:),dim=1)
!  end if maximum nuber of elements is reached
      if (zahle .ge. maxelemt) then
        ende=.true.
        return
      end if
!
      allocate(eplocal(n,nnat))
      eplocal = ep
!
!  Loop over the natures
      do inat=1,nnat
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
!  Sort residual vector. biggest values first. indx points to entry in array res
        call qsortindex(res(:,inat),indx,n)
!
!  class (A)
!  for elements with large error contribution,    either h- or p-adaption
        nlarge=n
!  psumres hold the partial sum of residuals 
        psumres=0._DP
!
        do j = 1, n
!  the sum of element errors in this class must not be larger than 
!  the fraction  large_errcont  of the total error 
          psumres=psumres+res(indx(j),inat)
          if (psumres .gt. large_errcont*sumsqrtres(inat)) then
            nlarge=j
            exit
          end if
!  if element has a keypoint or maximum polydegree is reached 
!  --> mark element for h-refinement
          if ( (kzi(e(1,indx(j))).lt.0) .or. (kzi(e(2,indx(j))).lt.0)    &
               .or. (kzi(e(3,indx(j))).lt.0) ) then
            call setze(ellist,zahle,zahlk,indx(j))
            eplocal(indx(j),inat) = ceiling(eplocal(indx(j),inat)/2._DP)
          else if (eplocal(indx(j),inat) .eq. polymax) then 
            call setze(ellist,zahle,zahlk,indx(j))
          else
!  determine whether mesh size is to large
!  compute square of element mesh size
            xelem(1:3)=xn(e(1:3,indx(j)))
            yelem(1:3)=yn(e(1:3,indx(j)))
            msize2=max( (xelem(1)-xelem(2))**2+(yelem(1)-yelem(2))**2,    &
            &           (xelem(2)-xelem(3))**2+(yelem(2)-yelem(3))**2,    &
            &           (xelem(3)-xelem(1))**2+(yelem(3)-yelem(1))**2 )
!  get coefficients of differential equation at the element
            xs=sum(xelem(1:3))/3._DP
            ys=sum(yelem(1:3))/3._DP
            call pdecoeff(indx(j),xs,ys,nu,gamma,alpha,beta,f1,f2)
!  it is critical if we have the type of a wave equation where alpha/nu is negative!
            k2=max( real(-alpha(inat,inat)/nu(1,1,inat,inat)) , real(-alpha(inat,inat)/nu(2,2,inat,inat)) )
            kh2=k2*msize2
!  the first part of the if statement checks for phase-lag condition kh < (2p+1)/kappa
!  while the second part (independent on polynomial order) improves stability of CG solvers
            if ( kh2 .ge. ((2*eplocal(indx(j),inat)+1)/kappa)**2  .or.              &
            &    kh2 .ge. khmax**2 ) then
!  h-refinement
              call setze(ellist,zahle,zahlk,indx(j))
              eplocal(indx(j),inat) = ceiling(eplocal(indx(j),inat)/2._DP)
            else
!  h-refinement is not required  --> do p-refinement
              eplocal(indx(j),inat) = eplocal(indx(j),inat)+1
            end if
          end if
        end do
!
!  class (B)
!  for elements with medium error contribution,    prefers  p-adaption
        do j = nlarge, n
!  the sum of element errors in this class must not be larger than 
!  the fraction  medium_errcont  of the total error 
          psumres=psumres+res(indx(j),inat)
          if (psumres .gt. medium_errcont*sumsqrtres(inat)) then
            nmedium=j
            exit
          end if
!  if maximum polydegree is reached 
!  --> mark element for h-refinement
          if ( eplocal(indx(j),inat) .eq. polymax ) then
            call setze(ellist,zahle,zahlk,indx(j))
          else
!  determine whether mesh size is to large
!  compute square of element mesh size
            xelem(1:3)=xn(e(1:3,indx(j)))
            yelem(1:3)=yn(e(1:3,indx(j)))
            msize2=max( (xelem(1)-xelem(2))**2+(yelem(1)-yelem(2))**2,    &
            &           (xelem(2)-xelem(3))**2+(yelem(2)-yelem(3))**2,    &
            &           (xelem(3)-xelem(1))**2+(yelem(3)-yelem(1))**2 )
!  get coefficients of differential equation at the element
            xs=sum(xelem(1:3))/3._DP
            ys=sum(yelem(1:3))/3._DP
            call pdecoeff(indx(j),xs,ys,nu,gamma,alpha,beta,f1,f2)
!  it is critical if we have the type of a wave equation where alpha/nu is negative!
            k2=max( real(-alpha(inat,inat)/nu(1,1,inat,inat)) , real(-alpha(inat,inat)/nu(2,2,inat,inat)) )
            kh2=k2*msize2
!  the first part of the if statement checks for phase-lag condition kh < (2p+1)/kappa
!  while the second part (independent on polynomial order) improves stability of CG solvers
            if ( kh2 .ge. ((2*eplocal(indx(j),inat)+1)/kappa)**2  .or.              &
            &    kh2 .ge. khmax**2 ) then
!  h-refinement
              call setze(ellist,zahle,zahlk,indx(j))
              eplocal(indx(j),inat) = ceiling(eplocal(indx(j),inat)/2._DP)
            else
!  h-refinement is not required  --> do p-refinement
              eplocal(indx(j),inat) = eplocal(indx(j),inat)+1
            end if
          end if
        end do
!
!  class (C) and (D)
!  decrease degree p for elements which marginally contribute to error
        do j = nmedium, n
          psumres=psumres+res(indx(j),inat)
!  determine whether mesh size is to large
!  compute square of element mesh size
          xelem(1:3)=xn(e(1:3,indx(j)))
          yelem(1:3)=yn(e(1:3,indx(j)))
          msize2=max( (xelem(1)-xelem(2))**2+(yelem(1)-yelem(2))**2,      &
          &           (xelem(2)-xelem(3))**2+(yelem(2)-yelem(3))**2,      &
          &           (xelem(3)-xelem(1))**2+(yelem(3)-yelem(1))**2 )
!  get coefficients of differential equation at the element
          xs=sum(xelem(1:3))/3._DP
          ys=sum(yelem(1:3))/3._DP
          call pdecoeff(indx(j),xs,ys,nu,gamma,alpha,beta,f1,f2)
!  it is critical if we have the type of a wave equation where alpha/nu is negative!
          k2=max( real(-alpha(inat,inat)/nu(1,1,inat,inat)) , real(-alpha(inat,inat)/nu(2,2,inat,inat)) )
          kh2=k2*msize2
!  here we check (independent on polynomial order) for stability of CG solvers
          if ( kh2 .ge. khmax**2 ) then
!  h-refinement
            call setze(ellist,zahle,zahlk,indx(j))
            eplocal(indx(j),inat) = ceiling(eplocal(indx(j),inat)/2._DP)
          else if (eplocal(indx(j),inat) .gt. 1  .and.                              &
          &      psumres .gt. (1._DP-marigial_errcont)*sumsqrtres(inat)) then
!  the sum of element errors in this class must not be larger than 
!  the fraction  marigial_errcont  of the total error 
!  h-refinement is not required  --> do p-UNrefinement
             eplocal(indx(j),inat) = eplocal(indx(j),inat)-1
          end if
        end do
!  END OF NATURES LOOP:
      end do
!
      ep = eplocal
      deallocate(eplocal)
!
      return
      end subroutine ffeinhp_KPphas
