      subroutine hpadapt(ext,int,ende,error,ende2,epsgl)
      use feminterface, only: wfein, ffeinhp, ffeinhp_KPphas, ffeinhp_top5
      use feminterface, only: getsetting, zeit, reallocate, glaett
      use feminterface, only: preassemb, assemblyandsolve
      use feminterface, only: newelements, transsol
      use globalvariables, only: ep, n, nnat, p, e, en, xn, yn, eg, &
                               & geb, kzi, x, zrb, zki, xbk, ybk, ndof, polymax, physics
      use femtypes
      implicit none
      real (DP) :: error, epsgl
      logical :: ende, ende2, ext, int
      intent(in) :: ext, int
      intent (out) :: ende, error, ende2, epsgl
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
!    $Revision: 1.31 $
!    $Date: 2015/04/01 10:58:13 $
!    $Author: juryanatzki $
!
!-------------------------------------------------------------------------------
!
!  Adaptive mesh- and polynomial refinement:
!  First of all an error estimation is done. Elements with higher residual than
!  adaption_error divided by number of elements are splitted into four subele- 
!  ments. Neighboring elements are splitted into two elements to restore confor-
!  mity if needed.
!
!  Input:
!            matzif   material index with which the region is filled
!            kdim     declared dimension for number of nodes
!            omega    angular frequency
!            zki      branch information, list of indices of key-points (geometry
!                     nodes) which determine the branches
!            kzrb     boundary condition of the key-points, index of the branch
!                     whose boundary condition is valid for the key-point
!            zrb      type of the boundary condition of the branch
!            xbk,ybk  coordinates of the key-points
!
!  Output:
!            ende     =.true. if MAXNODES is reached
!            error    relative error returned by a posteriori error estimation
!            ende2     =.true. if required tolerance is reached
!            epsgl    the error reached by the linear solver
!
!  In-/ Output: 
!            x        solution vector
!            e        node indices of the elements
!            en       neighbour-elements of the elements
!            n        total number of elements
!            xn,yn    coordinates of the nodes
!            p        total number of nodes
!            geb      assignment of elements to the geometry regions
!            kzi      node branch information
!
!  a maximum of n/5+1 elemente is refined.
!
!  local variables 
      integer (I4B) :: i, maxelemt, zahle, zahlk, fehlkn, inat
      integer (I4B) :: adaptsteps, polyold, polyorder, maxnodes
      integer (I4B), allocatable :: ellist(:), liste5(:)
      real (DP) :: adapterror, resgl
      logical :: nregen
      integer (I4B) :: j, k, m, nb, deltap
      logical :: stopdoing
      logical, allocatable  :: todolist(:)
      character (len=16) :: solver, algorithm
!  old mesh
      integer (I4B) :: nold, pold
      integer (I4B), pointer :: eold(:,:), enold(:,:), gebold(:)
      integer (I4B), pointer :: epold(:,:), kziold(:)
      real (DP), pointer :: xnold(:), ynold(:)
      complex (DPC), pointer :: xold(:)
      type (ARRPTRI), pointer :: egold(:,:)
!
!  Allocate vector for polynomial degrees of elements.
      if (.not. associated(ep)) then
!  if there is no ep, then there is no solution --> get starting value for
!  polynomial degree from FEMsettings.txt
        allocate (ep(n,nnat))
        call getsetting('POLYORDER',polyorder)
        ep(1:n,1:nnat) = polyorder
        call getsetting('PHYSICS_MODE',physics)
        select case(physics)
          case('FLUIDINCOMPR','FLUIDELMASSCON')
            select case(polyorder)
              case(1,2)
                ep(1:n,1:2) = 2
                ep(1:n,3)   = 1
              case default
                ep(1:n,3)   = polyorder-1
            end select            
          case default
            !No special element order modification
        end select
      else if (associated(ep)) then
!  there is an ep and polyorder in FEMsettings < mean value of the old vector
!  --> set initial value for polynomial order to minval(ep) 
        call getsetting('POLYORDER',polyorder)
        polyold = sum(ep)/(n*nnat)
        if (polyorder .gt. polyold) then
          ep(1:n,1:nnat) = polyorder
          print "(a)","Starting with polynomial order from FEMsettings.txt."
          print "(a,i2,a)","Polynomial order homogenously set to ",polyorder,"."
        end if
!  Special modification for the pressure field in fluidic problems
        call getsetting('PHYSICS_MODE',physics)
        select case(physics)
          case('FLUIDINCOMPR','FLUIDELMASSCON')
            if (polyorder .gt. polyold) then
              ep(1:n,3)   = polyorder-1
              print "(a,i2,a)","Polynomial order for fluidic pressure homogenously set to ",ep(1,3),"."
            end if
          case default
            !No special element order modification
        end select
      end if
!
      print *,'First Pre-assembly/solve...'
!
      call preassemb
!  Assemble and solve for the first time.
      call assemblyandsolve(epsgl,resgl)
!        
      call getsetting('MAXNODES',maxnodes)
      if (p .ge. maxnodes) then 
        print*,'** reached maximum number of nodes'
        ende=.true.
        return
      end if
      call getsetting('LINSOLVERTYPE',solver)
      solver = solver(1:len_trim(solver))
      call getsetting('ADAPTION_ERROR',adapterror)
      call getsetting('ADAPT_STEPS',adaptsteps)
      call getsetting('HP_ALGORITHM',algorithm)
      error=huge(1._DP)
!
      print *, '============================================================'
      print *, '+ Adaptations do to: ', adaptsteps
!
      do i=1,adaptsteps
        print *, '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,2) i,n,p
2       format('+ hp-adaptation step: ',i3,'   Elements=',i8,'   Nodes=',i8)
        maxelemt=n/5+1
!  refine elements having large inner angles
        allocate(ellist(n))
        call wfein(ellist,zahle,zahlk,maxelemt/2,ende)    
        allocate(epold(n,nnat))
        epold(1:n,1:nnat)=ep(1:n,1:nnat)
        if (.not.ende) then
!  h-adaption (subdivision of elements having large error)
          select case (algorithm)
            case ('KEYPOINT')
              print *, '+ keypoint'
              call ffeinhp(ellist,zahle,zahlk,ende,error,maxelemt,ext,int)
             case ('KP_PHASELAG')
              print *, '+ kp_phaselag'
              call ffeinhp_KPphas(ellist,zahle,zahlk,ende,error,maxelemt,ext,int)
            case ('TOP5')
              print *, '+ top5'
              call ffeinhp_top5(ellist,zahle,zahlk,ende,error,maxelemt,ext,int)
            case default
              print*, '*** unknown HP_ALGORITHM =',algorithm, 'using KEYPOINT'
              call ffeinhp(ellist,zahle,zahlk,ende,error,maxelemt,ext,int)
          end select
        end if
!
        print "(a,g11.4)", ' + Relative multi-nature a posteriori error: ', error
!
!  no more refinement if relative error is smaller or equal to given error tolerance
        if (error .le. adapterror)  then
          ende2= .true.
        else
          ende2= .false.
        end if
!
!  up to here we did tag elements for h-refinement now really modify the mesh
!  elements to be refined are stored in ellist
!
!  save the old mesh
        eold=>e
        enold=>en
        gebold=>geb
        egold=>eg
        xnold=>xn
        ynold=>yn
        kziold=>kzi
        xold=>x
        nold=n
        pold=p
        if (n+zahle .gt. size(ep,1))   ep=>reallocate(ep,n+zahle,nnat)
        allocate(e(3,n+zahle))
        allocate(en(3,n+zahle))
        allocate(geb(n+zahle))
        allocate(xn(p+zahlk))
        allocate(yn(p+zahlk))
        allocate(kzi(p+zahlk))
        allocate(x(p+zahlk))
        e(:,1:n)=eold(:,1:n)
        en(:,1:n)=enold(:,1:n)
        geb(1:n)=gebold(1:n)
        nullify(eg)
        xn(1:p)=xnold(1:p)
        yn(1:p)=ynold(1:p)
        kzi(1:p)=kziold(1:p)
        x(1:p)=xold(1:p)
!  now do the refinement
!  liste5 marks the nodes requiring node delocation in the smoothing process
!  may happen on narrow arcs
        allocate(liste5(zahlk+p))
        liste5(:)=0
        print *, '+ New elements...'
        call newelements(ellist,nregen,fehlkn,liste5)
        deallocate(ellist)
!
        call zeit(' Generating Elements')
!  mesh smoothing
        print *, '+ Mesh-smoothing...'
        call glaett(n,e,xn,yn,en,geb,kzi,p,nregen,liste5,fehlkn,zrb,    &
     &    zki,xbk,ybk,2)
        deallocate(liste5)
!
!  Look if polynomial order of neighbors differ from the element's more than 1.
!  Change lower polyorder in corresponding element if necessary.
!  loop over all elements
        deltap = 1
        allocate (todolist(n))
!  Loop over the natures
        do inat = 1,nnat
!  Skip element order change for pressure in fluidic problems
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
!  End loop of Natures
        end do
        deallocate (todolist)
!
!  Arrange special element order for fluidic pressure
!  loop over all elements
        select case(physics)
          case('FLUIDINCOMPR','FLUIDELMASSCON')
            do k = 1,n
              select case(minval(ep(k,1:2)))
                case (2)
                  ep(k,3) = 1
                case (1)
                  ep(k,1:2) = 2
                  ep(k,3) = 1
                case default
                  ep(k,3) = minval(ep(k,1:2))-1
              end select
            end do
          case default
            ! No special element order modification
        end select
!
        print "(a,i2,a,i2,a,i2)","mean poly: ",(sum(ep(1:n,1:nnat))/(n*nnat)), &
     &                           " lowest poly: ",minval(ep(1:n,1:nnat)),      &
     &                           " highest poly: ",maxval(ep(1:n,1:nnat))
!  solve
!
        print *, '+ Pre-assembly...' 
!
        call preassemb
        x=>reallocate(x,ndof)
        x(p+zahlk:ndof)=0._DP
!  transfer of old solution to the new starting vector
        if ((solver .ne. 'UMF') .and. (solver .ne. 'UMFC') .and. (solver .ne. 'PARDISO')) then
          call transsol(xnold,ynold,eold,enold,epold,egold,xold,nold)
          call zeit(' Solution transfer')
        end if
!  kill the old mesh
        deallocate(eold, enold, gebold, epold, xnold, ynold, kziold, xold)
        do inat = 1, nnat
          do j = 1, size(egold(:,inat))
            deallocate(egold(j,inat)%d)
          end do
        end do
        deallocate(egold)
!
        print*,'+ Assembly and solve...'
        call assemblyandsolve(epsgl,resgl,i)
        if (ende2) then
          print "(a,g10.3)",  " + Reached stopping criterion with error=", error
          print "(a,i3,a,i3)"," + Stopping on iteration ",i, " of", adaptsteps
          print*, '============================================================'
          return
        end if
        if (p .ge. maxnodes) then
          print*,'** reached maximum number of nodes'
          print*, '============================================================'
          ende=.true.
          return
        end if
      end do
!      
!  Not return by stop criteria nor max number of nodes
      print*, 'Completed :', adaptsteps, ' adaptations.' 
      print*, '============================================================'
!      
      return
      end subroutine hpadapt
!
!
!
      subroutine transsol(xnold,ynold,eold,enold,epold,egold,xold,nold)
      use feminterface, only: elemnt, shapefunction, xy2lam, lusolver, get2Dintegpoints
      use globalvariables, only: n, e, xn, yn, nnat, ep, eg, x
      use femtypes
      implicit none
      integer (I4B) :: nold
      integer (I4B) :: eold(:,:), enold(:,:)
      integer (I4B) :: epold(:,:)
      real (DP) :: xnold(:), ynold(:)
      complex (DPC) :: xold(:)
      type (ARRPTRI) :: egold(:,:)
      intent (in) :: xnold, ynold, eold, enold, epold, egold, xold, nold
!
!  transfer of old solution to the new starting vector
!  local variables
      integer (I4B) :: i, j, k, l, nff, dof, ielem, errcode
      integer (I4B) :: intorder, npkt, elem, nffold, inat
      real (DP) :: lambda(3)
      real (DP) :: area, aw, alphaaw, xs, ys, q3
      complex (DPC), pointer :: a(:,:), b(:)
      complex (DPC) :: f1, f1aw 
      real (DP), allocatable :: weight(:), lam(:,:)
      real (DP), allocatable :: xsi(:), gxsi(:,:), xsiold(:)
      logical :: ok
!
!
!  for all elements
!  avoids error with oneAPI. dummy size
      allocate( gxsi(10,2) )
      do elem=1,n
        area=( (yn(e(3,elem))-yn(e(1,elem)))*                           &
     &         (xn(e(2,elem))-xn(e(3,elem)))-                           &
     &         (yn(e(2,elem))-yn(e(3,elem)))*                           &
     &         (xn(e(3,elem))-xn(e(1,elem))) )/2._DP
!
        do inat=1,nnat
!  size of the element matrix
          nff=(ep(elem,inat)+1)*(ep(elem,inat)+2)/2
          allocate(a(nff,nff), b(nff), xsi(nff))
          a(:,:)=(0._DP,0._DP)
          b(:)=(0._DP,0._DP)
!  determine the order of numerical integration
          intorder=2*ep(elem,inat)
!  fetch numerical integration points (Gauss points) lam
          call get2Dintegpoints(intorder,npkt,weight,lam,errcode)
!
          do k=1,npkt
!  get shape function at location of integration points
            call shapefunction(lam(:,k),xn(e(:,elem)),yn(e(:,elem)),    &
     &      1,ep(elem,inat),nff,.false.,xsi,gxsi,errcode)
!
!  compute f1, the old mesh value at this point

!  world coordinates in new mesh
            xs= dot_product( xn(e(:,elem)) , lam(:,k) )
            ys= dot_product( yn(e(:,elem)) , lam(:,k) )
            call elemnt(xnold,ynold,eold,nold,xs,ys,ielem,enold,ok)
!  get barycentric coordinates lambda of point in the old mesh            
            call xy2lam(xs,ys,ielem,lambda,xnold,ynold,eold)
            nffold=size(egold(ielem,inat)%d)
            allocate(xsiold(nffold))
            call shapefunction(lambda,xnold(eold(:,ielem)),ynold(eold(:,ielem)),&
              1,epold(ielem,inat),nffold,.false.,xsiold,gxsi,errcode)
            f1=0._DPC
            do l=1,nffold
              dof=egold(ielem,inat)%d(l)
              select case(dof)
              case (1:)
                f1 = f1 + xsiold(l)*xold(dof)
              case (:-1)
                f1 = f1 - xsiold(l)*xold(-dof)
              end select
            end do
            deallocate(xsiold)
            aw=area*weight(k)
            alphaaw=aw
            f1aw=f1*aw
!  for all shape functions
            do i=1,nff
              q3=alphaaw*xsi(i)
              do j=1,nff
!  integrate the function:
                a(i,j)=a(i,j) + q3*xsi(j)
              end do
!  right hand side contributions, integrate the function:
              b(i)=b(i)+f1aw*xsi(i)
            end do
          end do
          deallocate(weight,lam,xsi)
!  solve
          call lusolver(a,b)
!  set        
          do l=1,nff
            if (eg(elem,inat)%d(l).ne.0) then
              if (eg(elem,inat)%d(l).gt.0) then
                x(eg(elem,inat)%d(l))=b(l)
              else 
                x(-(eg(elem,inat)%d(l)))=-b(l)
              end if
            end if
          end do
          deallocate(a, b)
        end do
      end do
      deallocate( gxsi)

      return
      end subroutine transsol
