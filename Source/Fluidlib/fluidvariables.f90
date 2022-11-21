module fluidvariables
use femtypes
use globalvariables
implicit none
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
!    $Date: 2006/07/07 22:34:04 $
!    $Author: r_abdul $
!
!*******************************************************************************
! Module for fluid global variables
!*******************************************************************************
! In compressible or incompressible flow for 2D, there are always at least
! three variables.
!
!   incompressible : u, v and p
!   compressible : rho, rho*u, rho*v, rho*E
!    
! where pressure (p) in case of compressible flow can be found using, for
! example, the ideal gas law, p = rho*R*T
!*******************************************************************************
! starting option: 0 = assume all parameters, 1 = starting unknown values are given
integer(I4B) :: iopt

! no. of step run so far (usually start at 1)
integer(I4B) :: istep

! no. of step to write the residues
integer(I4B) :: iwrite

! no. of mass iterations for consistent mass matrix computation
integer(I4B) :: niter

! no of step to be run.
integer(I4B) :: ntime 

! no of unknowns or variables
integer(I4B) :: nvar

! no. of neighbouring nodes, neigh.node number corresponding to the node 
integer(I4B), pointer :: ncmax(:), matcon(:,:)

! time step options
integer(I4B) :: optts !ilots

! inverse of Re and Rho (density)
real(DP) :: invre, invrho

! speed of sound
! cc = sqrt( (gamma -1)/Tinf )
real(DP) :: cc

! safety factor
real(DP) :: csafm

! assumed values of rho, velocities.
real(DP) :: rassum, uassum, vassum

! guessed values or free stream values for pressure and temperature
real(DP) :: pinf, tinf

! constant for residual smoothing and no. of smoothing loops
real(DP) :: csmoo, nsmoo

! time step for nodes and elements
real(DP) :: dtfix
real(DP), pointer :: deltp(:), delte(:)

! gamma, gamma-1
real(DP) :: gamma, gam1

! ratio of stagnation temp. to free stream temp.
real(DP) :: ratio

! Reynolds number (Re), Prandtl number (Pr), Rayleigh number(Ra), Richarson number(Ri)
real(DP) :: re,pr,ra,ri

! Average Nusselt number
real(DP) :: nuaver

! coefficient of Suntherland's law for temperature-dependent viscosity
real(DP) :: su
real(DP) :: sucon, suth

! free stream temperature (???in Rankine for compute temperature-dependent density)
real(DP) :: tfree

! non.dim time elapsed
real(DP) :: timt

! free stream values of rho, u1, u2, p and Mach
! note : pinf is anyway calculated from other parameters
real(DP) :: cinf(5)

! gamma(Cp/Cv), artificial diffustion constant
real(DP) :: conin(2)

! relaxation parameters for velocities and pressure
real(DP) :: theta(2)

! pressure, pressure at (n-1)-th: pres(ndof), pres1(ndof)
real(DP), pointer :: pres(:), pres1(:)

! local speed of sound: sound(ndof)
real(DP), pointer :: sound(:)

! temp. at (n-1)-th: tt1(ndof)
real(DP), pointer :: tt1(:)

! element properties: alen(n), areaf(n)
real(DP), pointer :: alen(:)
real(DP), pointer :: areaf(:)

! inverse of lumped mass matrix: dmmat(ndof)
real(DP), pointer :: dmmat(:)

! unkno(nvar,ndof) : unknown at n-th iteration
real(DP), pointer :: unkno(:,:)

! unkn1(nvar,ndof) : unknown at (n-1)-th iteration
real(DP), pointer :: unkn1(:,:)

! unkel(nvar,neldof) : unknown at element
real(DP), pointer :: unkel(:,:)

! stream function(phi)
real(DP), pointer :: phi(:)

! magnetic field: magn(ndof)
real(DP), pointer :: magn(:)

real(DP) :: curie, magconst
real(DP) :: twall, tmpdiff    ! wall temperature (hot), reference temperature diff.
real(DP) :: tfilm    ! reference film temperature = (twall+tfree)/2
real(DP) :: lchar    ! characteristic length
real(DP) :: visco    ! viscosity
real(DP) :: thcond   ! thermal conductivity(k)
real(DP) :: thdiff   ! thermal diffusivity = k/rho*c_p
real(DP) :: thexpn   ! coefficient of thermal expansion (ideal gas = 1/tfilm)
real(DP), parameter :: gasconst = 8.3144_DP
real(DP), parameter :: gravity = 9.80665_DP 
! compressible options, F = incompressible, T = compressible
logical :: optcomprs, optener, optstokes, optdim
logical :: optmixed, optstrm, optnat, optmagn
!*******************************************************************************
! element information 
!*******************************************************************************
! no. of multinature b.c.
integer(I4B) :: nbc

! element branch info for fluid (other physics do not need to store globally)
integer(I4B), pointer :: branchef(:,:)    ! branchef(nedge,n)

! netinmod.acd
integer(I4B), pointer :: kzrbmn(:,:)      ! BCs of keypoints in multinature; kzrbmn(gkz,nvar)
integer(I4B), pointer :: zrbmn(:,:)       ! BCs of branches in multinature; zrbmn(gzz,nvar)

! edge lengths of element
real(DP), pointer :: lengthby2(:,:)       ! length/2 (nedge,nel)

complex(DPC), pointer :: alrbmn(:,:)      ! alpha for bc; alrbmn(gzz,nvar)
complex(DPC), pointer :: btrbmn(:,:)      ! beta for bc; btrbmn(gzz,nvar)

! boundary-condition related variables
type :: multinature   
   integer(I4B) :: branch                       ! On which branch the node is.
   integer(I4B), dimension(:), pointer :: brmn  ! Branch of keypoint for each variables.
   logical, dimension(:), pointer :: varmn      ! Variables 
end type multinature

type :: bctype
   type(multinature), dimension(3) :: bcedge
   type(multinature), dimension(3) :: bcdof
end type bctype

type(bctype), pointer :: generalbc(:)
type(bctype), pointer :: dirichletbc(:)
logical, pointer :: genbc(:)
logical, pointer :: edbc(:)
!*******************************************************************************
! boundary array
!*******************************************************************************
integer(I4B) :: nbound                 ! no. of boundary side of all elements
integer(I4B), pointer :: bcside(:,:)   ! bcside(4,nbound) 
                                       ! 1 = elem no.,
                                       ! 2 = edge,
                                       ! 3 = node1,
                                       ! 4 = node2, 
                                       ! 5 = bc code
real(DP), pointer :: nside(:,:)        ! nside(3,nbound)
                                       ! 1 = nx, 2 = ny, 3 = leng 
!*******************************************************************************
! fixed dirictlet b.c. on node
!*******************************************************************************
integer(I4B) :: npbc                   ! no. of dirichlet pressure
integer(I4B), pointer :: pkp(:)        ! keypoint to be fixed for pressure
integer(I4B), pointer :: zrbpbc(:)     ! bc code
real(DP), pointer :: xpbc(:), ypbc(:)  ! coordinate to be fixed
real(DP), pointer :: varpbc(:)
real(DP), pointer :: bglobal(:)
!*******************************************************************************
! temporary variables
!*******************************************************************************
! in geteleminfo and timestep_incomprs
real(DP), pointer :: elh(:)   ! elh(n) - shortest element length of each element
!*******************************************************************************
end module fluidvariables