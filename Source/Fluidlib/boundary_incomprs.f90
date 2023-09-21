subroutine boundary_incomprs(optvar)
use femtypes
use fluidvariables
use feminterface, only : get1Dintpolpoints
use fluidinterface, only: getbcval2D_mn, solidbc1
use globalvariables
implicit none
integer(I4B), optional, intent(in):: optvar  ! 0 = u,v, 1 = T
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
!    $Revision: 1.6 $
!    $Date: 2008/08/20 12:53:46 $
!    $Author: m_kasper $
!
!***************************************************************************************
!
!  Imposing Dirichlet Boundary Conditions - Incompressible
!
!***************************************************************************************
! need to check variables to be set with ivar: 
!  p = var1
!  u = var2
!  v = var3
!  T = var4
!***************************************************************************************
! dirichletbc(n)%bcedge(3)%varmn(nvar) - tell if the corresponding variable number 
!                                        on edge is prescribed                      
!               %bcdof(3)%varmn(nvar)  - tell if the corresponding variable number
!                                        on vertex is prescribed
!***************************************************************************************
! local variable
!***************************************************************************************
integer(I4B) :: brvertex, brvmn
integer(I4B) :: i, ie, ii, ivertex, ivar, avar, bvar
integer(I4B) :: bredge, iedge, intpolpoints, j, ln, nr, nl
integer(I4B) :: nff, polyorder, errcode
integer(I4B), dimension(3), parameter :: inach = (/2,3,1/)
integer(I4B), dimension(3), parameter :: ivorg = (/3,1,2/)
real(DP) :: xs, ys
real(DP) :: lam(3)
real (DP), pointer :: xtab(:)
complex(DP) :: value ! complex in Polyde
logical :: dbclog, implc
!
external    ::    inoutflowbc1
!***************************************************************************************
implc = (theta(2) .GT. 0)

do ie = 1, n
   if (.NOT. edbc(ie)) cycle
   polyorder = ep(ie,1)
   nff = (polyorder+1)*(polyorder+2)/2

   do ivertex = 1,3
      ii = eg(ie,1)%d(ivertex)
      brvertex = dirichletbc(ie)%bcdof(ivertex)%branch
      if (brvertex .EQ. 0) cycle

      ! check userbc for special type of dirichelt b.c.
      ! imposed values of differenct variables will be calculated.

      if (zrb(brvertex,1) .EQ. 97) then
         ! inflow or outflow b.c.: all prescribed u,v,p and T will be calculated
         ! from the given free stream values.(cinf(:))
         call inoutflowbc1(ii)            

      else if (zrb(brvertex,1) .EQ. 98) then
         ! no-slip b.c.
         call solidbc1(ii,optvar,zrb(brvertex,1))

      else if (zrb(brvertex,1) .EQ. 99) then
         ! no-slip b.c.
         call solidbc1(ii,optvar,zrb(brvertex,1))

      else
         select case (optvar)
         case(0)
            avar = 2
            bvar = 3
         case(1)
            avar = 4
            bvar = 4
         case(2)
            avar = 1
            bvar = 1
         end select

         do ivar = avar,bvar
!            if (ivar .EQ. 1 .AND. implc) cycle
            dbclog = dirichletbc(ie)%bcdof(ivertex)%varmn(ivar)
            if (dbclog) then

               xs= xn(e(ivertex,ie))
               ys= yn(e(ivertex,ie))

               brvmn = dirichletbc(ie)%bcdof(ivertex)%brmn(ivar)

               call getbcval2D_mn(ivar,brvmn,xs,ys,value)

               unkno(ivar,ii) = real(value,DP)   

            end if
         end do ! ivar
      end if ! zrb     
   end do ! ivertex   
   ! polyorder greater than 1
   if (polyorder .GE. 2) then
      intpolpoints=polyorder+1
      allocate(xtab(intpolpoints))
      call get1Dintpolpoints(intpolpoints, xtab, errcode)

      do iedge = 1,3        
         bredge = branchef(iedge,ie)
         if (bredge .eq. 0) cycle
         nr=inach(iedge)
         nl=ivorg(iedge)
         lam(iedge)=0._DP
         do i=2,intpolpoints-1
!  evaluate BC at Gauss-Lobatto points
            lam(nr)=(xtab(i)+1._DP)/2._DP
            lam(nl)=1._DP-lam(nr)

            do j = 2, polyorder
               ln = j*(j+1)/2+inach(iedge)
               ii = eg(ie,1)%d(ln)            
               if (zrb(bredge,1) .EQ. 97) then
               ! inflow or outflow b.c.
                  call inoutflowbc1(ii)            
               else if (zrb(bredge,1) .EQ. 98) then
                  ! no-slip b.c.
                  call solidbc1(ii,optvar,zrb(bredge,1))
               else if (zrb(bredge,1) .EQ. 99) then
                  ! slip b.c.
                  call solidbc1(ii,optvar,zrb(bredge,1))               
               else
                  select case (optvar)
                  case(0)
                     avar = 2
                     bvar = 3
                  case(1)
                     avar = 4
                     bvar = 4
                  case(2)
                     avar = 1
                     bvar = 1
                  end select

                  do ivar = avar,bvar                  
                     dbclog = dirichletbc(ie)%bcedge(iedge)%varmn(ivar)
                     if (.NOT. dbclog) cycle

                     xs= dot_product( xn(e(:,ie)) , lam(:) )
                     ys= dot_product( yn(e(:,ie)) , lam(:) )
                     call getbcval2D_mn(ivar,bredge,xs,ys,value)
                  
                     unkno(ivar,ii) = real(value,DP)
                  end do ! ivar
               end if ! zrb(bredge)
            end do ! j
         end do ! i
      end do ! iedge
      deallocate(xtab)   
   end if ! polyorder
end do ! ie

end subroutine boundary_incomprs
!################################################################################
subroutine inoutflowbc1(glbnum)
use femtypes
use fluidvariables
implicit none
integer(I4B), intent(in) :: glbnum
! local variables
!real(DP) :: ener, ener1, ener2, twall
!*************************************************************************************
! see Zienkiewicz 2000, p. 67: FEM in Fluid
!*************************************************************************************
! velocity
unkno(2,glbnum) = cinf(2)  ! usually u = 1.0_DP
unkno(3,glbnum) = cinf(3)  !         v = 0.0_DP

end subroutine inoutflowbc1
!#####################################################################################
subroutine solidbc1(glbnum,opt,zrbbc)
use femtypes
use fluidvariables
implicit none
integer(I4B), intent(in) :: glbnum
integer(I4B), optional, intent(in) :: opt, zrbbc
! local variables
!real(DP) :: ener, ener1, ener2, twall
!*************************************************************************************
! see Zienkiewicz 2000, p. 67: FEM in Fluid
!*************************************************************************************
! no-slip condition & no temperature-jump condition

select case(opt)
case(0)
! velocity for fixed wall
   unkno(2,glbnum) = 0._DP
   unkno(3,glbnum) = 0._DP

case(1)
! temperature
   if (zrbbc .EQ. 98) then
      unkno(4,glbnum) = 1.0_DP
   else if (zrbbc .EQ. 99) then
      unkno(4,glbnum) = 0.0_DP
   end if

case default
   write(*,*) 'Error: Invalid variable option in boundary_incomprs'
   STOP

end select

end subroutine solidbc1