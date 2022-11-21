module fluidinterface
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
!    $Revision: 1.7 $
!    $Date: 2006/07/07 22:34:04 $
!    $Author: r_abdul $
!
!****************************************************************************************
! This module is seperated into 4 parts
! 1. common subroutines
! 2. incompressible-related subroutines
! 3. stokes-related subroutines
! 4. compressible-related subroutines
!****************************************************************************************
!  COMMON
!****************************************************************************************
interface
   subroutine assembly_streamfn(lower,upper,diag,rhs,ia,ja,symmetric)
   use femtypes
   implicit none
   integer(I4B), pointer :: ia(:), ja(:)
   real(DP), pointer :: lower(:), upper(:), diag(:), rhs(:)
   logical, intent(in) :: symmetric
   end subroutine
end interface

interface
   subroutine checkconverge(eps,err,converge)
   use femtypes
   implicit none
   real(DP), intent(in) :: eps
   real(DP), pointer :: err(:)
   logical, intent(out) :: converge
   end subroutine
end interface

interface
   subroutine elementmatrix_streamfn(elem,polyorder,a,b,nff)
   use femtypes
   implicit none
   integer (I4B) :: elem
   integer (I4B) :: polyorder, nff
   real(DP), pointer :: a(:,:), b(:)
   intent (in) :: elem, polyorder
   intent (out) :: nff
   end subroutine
end interface

interface
      subroutine getbcval2D_mn(ivar, branch, xs, ys, pval, qval )
      use femtypes
      implicit none
      integer (I4B) :: branch, ivar
      real (DP) :: xs, ys
      complex(DPC) :: pval
      complex(DPC), optional :: qval
      intent (in)   :: branch, xs, ys, ivar
      intent (out)  :: pval, qval
      end subroutine getbcval2D_mn
end interface

interface
   subroutine lout_fluid(gilt, eleinfo, eps)
   use femtypes
   implicit none
   real(DP) :: eps
   logical  :: gilt, eleinfo
   intent(in) :: gilt, eleinfo, eps
   end subroutine
end interface

interface
   subroutine streamfunction(ia,ja)
   use femtypes
   implicit none
   integer(I4B), pointer :: ia(:), ja(:)
   end subroutine
end interface
!****************************************************************************************
!  INCOMPRESSIBLE FLOW
!****************************************************************************************
interface
   subroutine artificialcomprs(beta)
   use femtypes
   implicit none
   real(DP), pointer :: beta(:)
   end subroutine
end interface

interface
   subroutine assemblyrhsstep2(rhs)
   use femtypes
   implicit none
   real(DP), pointer :: rhs(:)
   end subroutine
end interface

interface
   subroutine boundary_incomprs(optvar)
   use femtypes
   implicit none
   integer(I4B), optional, intent(in) :: optvar
   end subroutine
end interface

interface
   subroutine elementmatrix_pstiff(elem,polyorder,a,nff)
   use femtypes
   implicit none
   integer (I4B) :: elem
   integer (I4B) :: polyorder, nff
   real(DP), pointer :: a(:,:)  ! intent(out)
   intent (in) :: elem, polyorder
   intent (out) :: nff
   end subroutine
end interface

interface
   subroutine getrhsstep1_incomprs(rhs)
   use femtypes
   implicit none
   real(DP), pointer :: rhs(:,:) 
   end subroutine
end interface

interface
   subroutine getrhsstep2_implicit(elem,polyorder,b,nff)
   use femtypes
   implicit none
   integer (I4B) :: elem
   integer (I4B) :: polyorder, nff
   real(DP), pointer :: b(:)   ! intent(out)
   intent (in) :: elem, polyorder
   intent (out) :: nff
   end subroutine
end interface

interface
   subroutine getrhsstep2_incomprs(rhs)
   use femtypes
   implicit none
   real(DP), pointer :: rhs(:,:)
   end subroutine
end interface

interface
   subroutine getrhsstep3_incomprs(rhs)
   use femtypes
   implicit none
   real(DP), pointer :: rhs(:,:) 
   end subroutine
end interface

interface
   subroutine getrhsstep4_incomprs(rhs)
   use femtypes
   implicit none
   real(DP), pointer :: rhs(:,:) 
   end subroutine
end interface

interface
   subroutine nusselt
   use femtypes
   implicit none
   end subroutine
end interface

interface
   subroutine pstiff(lower,diag,ia,ja)
   use femtypes
   implicit none
   integer (I4B), pointer :: ia(:), ja(:)
   real(DP), pointer :: lower(:), diag(:)
   end subroutine
end interface

interface   ! Calculate intermediate velocity
   subroutine step1_incomprs(rhs)
   use femtypes
   implicit none
   real(DP), pointer :: rhs(:,:)
   end subroutine
end interface

interface
   subroutine step2_implicit(lower,diag,ia,ja)
   use femtypes
   implicit none
   integer(I4B), pointer :: ia(:), ja(:)
   real(DP), pointer :: lower(:), diag(:)
   end subroutine
end interface

interface   ! Calculate pressure explicitly
   subroutine step2_incomprs(rhs)
   use femtypes
   implicit none
   real(DP), pointer :: rhs(:,:)
   end subroutine
end interface

interface   !  Velocity correction
   subroutine step3_incomprs(rhs)
   use femtypes
   implicit none
   real(DP), pointer :: rhs(:,:) ! intent(out)
   end subroutine
end interface

interface
   subroutine step4_incomprs(rhs)
   use femtypes
   implicit none
   real(DP), pointer :: rhs(:,:)
   end subroutine
end interface
!****************************************************************************************
!  STOKES FLOW
!****************************************************************************************
interface
   subroutine assembly_nostep(lower,upper,diag,acsr,rhs,ia,ja,csr,npdof)
   use femtypes
   implicit none
   integer (I4B), pointer :: ia(:), ja(:)
   integer(I4B), optional :: npdof
   real(DP), pointer :: lower(:), upper(:), diag(:)
   real(DP), pointer :: rhs(:), acsr(:)
   logical, intent(in) :: csr
   end subroutine
end interface

!Intentionally left out. Waiting for the correct implementation
!interface
!   subroutine assembly_stokes(lower,upper,diag,acsr,rhs,ia,ja,jacobi,csr,matvar,ntdof)
!   use femtypes
!   implicit none
!   integer(I4B) :: ntdof
!   integer(I4B), pointer :: ia(:), ja(:)
!   real(DP), pointer :: lower(:), upper(:), diag(:), rhs(:), acsr(:)
!   logical :: jacobi, csr, matvar
!   intent(in) :: jacobi, csr, matvar, ntdof
!   end subroutine
!end interface

interface
   subroutine elementmatrix_stokes(elem,KE,FE,matvar)
   use femtypes
   implicit none
   integer(I4B), intent(in) :: elem
   real(DP), pointer :: KE(:,:), FE(:)
   logical, intent(in) :: matvar
   end subroutine
end interface

interface
   subroutine elementmatrix_stokesmix(elem,polyorder,jacobi,full,matvar,KE,FE,nff,errcode)
   use femtypes
   implicit none
   integer (I4B) :: elem
   integer (I4B) :: polyorder, nff, errcode
   real(DP), pointer :: KE(:,:), FE(:)
   logical :: jacobi, full, matvar
   intent (in) :: elem, polyorder, jacobi, full, matvar
   intent (out) :: nff
   end subroutine
end interface

interface
   subroutine lout_mixed(gilt, eleinfo, eps, npdof)
   use femtypes
   implicit none
   integer(I4B) :: npdof
   real(DP) :: eps
   logical :: gilt, eleinfo
   intent (in) :: gilt, eleinfo, npdof, eps
   end subroutine
end interface

interface
   subroutine preassemb_mixed(nudof,npdof)
   use femtypes
   implicit none
   integer(I4B), intent(out) :: nudof, npdof
   end subroutine
end interface

interface
   subroutine setupia_stokes(ia,csr)
   use femtypes
   implicit none
   integer(I4B), intent(out) :: ia(:)
   logical, intent(in) :: csr
   end subroutine
end interface

interface
   subroutine setupia_temp(ia,csr,ntdof)
   use femtypes
   implicit none
   integer(I4B), intent(OUT) :: ia(:)
   integer(I4B), intent(in) :: ntdof
   logical, intent(in) :: csr
   end subroutine
end interface

!****************************************************************************************
!  COMPRESSIBLE FLOW
!****************************************************************************************
interface ! Compute the advective term of the momentum equation
   subroutine advectstep1(rhs2)
   use femtypes
   implicit none
   real(DP), pointer :: rhs2(:,:)
   end subroutine
end interface

interface !  Determine RHS of step 4 only for convective term
   subroutine enerconv(rhs2)
   use femtypes
   real(DP), pointer :: rhs2(:,:) ! intent(inout)
   end subroutine
end interface

interface !  Calculates pressure from velocity and energy.
   subroutine getpres(pre)
   use femtypes
   implicit none
   real(DP), pointer :: pre(:) ! intent(out)
   end subroutine
end interface

interface   ! Calculate mass iterations
   subroutine getrh0(rhs2,ia1,ia2)
   use femtypes
   implicit none
   integer(I4B), intent(in) :: ia1,ia2
   real(DP), pointer :: rhs2(:,:) ! intent(inout)
   end subroutine
end interface

interface !  Determine RHS of step 1 (intermediate velocity)
   subroutine getrhsstep1(rhs2)
   use femtypes
   real(DP), pointer :: rhs2(:,:) ! intent(inout)
   end subroutine
end interface

interface !  Determine RHS of step 2 (pressure or density equation)
   subroutine getrhsstep2(rhs2)
   use femtypes
   implicit none
   real(DP), pointer :: rhs2(:,:) ! intent(inout)
   end subroutine
end interface

interface !  Determine RHS of step 3 (momentum equation)
   subroutine getrhsstep3(rhs2)
   use femtypes
   implicit none
   real(DP), pointer :: rhs2(:,:) ! intent(inout)
   end subroutine
end interface

interface !  Determine RHS of step 4 (energy equation)
   subroutine getrhsstep4(rhs2)
   use femtypes
   real(DP), pointer :: rhs2(:,:) ! intent(inout)
   end subroutine
end interface

interface
   subroutine solver_comprs(eps)
   use femtypes
   implicit none
   real(DP), intent(in) :: eps
   end subroutine
end interface

interface   ! Calculate intermediate velocity
   subroutine step1(itime, intime, rhs2, rhs1)
   use femtypes
   implicit none
   real(DP), pointer :: rhs2(:,:)
   real(DP), pointer, optional :: rhs1(:,:)
   integer(I4B) :: itime, intime
   intent(in) :: itime, intime
   !intent(out) :: rhs2, rhs1
   end subroutine
end interface

interface   ! Calculate pressure or density
   subroutine step2(rhs2,rhs1)
   use femtypes
   implicit none
   real(DP), pointer :: rhs2(:,:)
   real(DP), pointer, optional :: rhs1(:,:)
   end subroutine
end interface

interface   !  Velocity correction
   subroutine step3(rhs2)
   use femtypes
   implicit none
   real(DP), pointer :: rhs2(:,:) ! intent(out)
   end subroutine
end interface

interface
   subroutine step4(rhs2,rhs1)
   use femtypes
   implicit none
   real(DP), pointer :: rhs2(:,:)
   real(DP), pointer, optional :: rhs1(:,:)
   end subroutine
end interface

interface
   subroutine transform(iflag)
!   use femtypes
   implicit none
   logical :: iflag
   intent(in) :: iflag
   end subroutine
end interface

interface
      subroutine solidbc1(glbnum,opt,zrbbc)
      use femtypes
      implicit none
      integer(I4B), intent(in) :: glbnum
      integer(I4B), optional, intent(in) :: opt, zrbbc
      end subroutine
end interface

end module fluidinterface