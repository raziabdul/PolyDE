      subroutine pdecoeff3D_vector(elem, xyzs, nu, alpha, beta, f1, f2)
      use feminterface, only: getsetting, fetchmatparameters
      use feminterface3d, only: calctensor3D
      use femtypes
      use globalvariables3D, only: dom, eps0, mu0, dommat, omega, physics, pi, Kboltz
      use matconstants
      implicit none
      integer (I4B) elem
      real(DP) :: xyzs(3)
      complex (DPC) :: nu(3,3), alpha(3,3), beta(3)
      complex (DPC) :: f1(3), f2(3)
      intent (in) :: elem, xyzs
      intent (out) :: nu, alpha, beta, f1, f2
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
!    $Revision: 1.27 $
!    $Date: 2015/11/11 17:32:57 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------
!  For vector-valued basis functions
!--------------------------------------------------------------------------------------
!
!  deliver the material coefficients in the form of the general PDE
!  with:
!  (1) vector potential (default):
!  
!        -curl(nu*curl(U)) + (beta*curl(U)) + alpha*U = f1 + curl(f2)
!
!
!  Input:
!            elem     number of the element
!            x        actual solution vector [not used yet] (GLOBAL)
!            e        element information (nodes of the element) [not used yet] (GLOBAL)
!            matzif   material name of the input regions (GLOBAL)
!            geb      region index of the elements  (GLOBAL)
!            xyzs     location at which the material coefficients have to be evaluated
!                     this is usually the centroid of the volume element
!                     [not used yet]
!  Output:
!            nu       rank 2 tensor of the diffusion term:      { nu_xx  nu_xy  nu_xz }
!                                                          nu = | nu_yx  nu_yy  nu_yz |
!                                                               { nu_zx  nu_zy  nu_zz }
!            beta     vector of the convective term beta(1)=betax
!            alpha    rank 2 tensor, like nu 
!            f1       right hand side vector
!            f2       right hand side source vector
!
!  local variables
      integer (I4B) matindex, pml
      real(DP), allocatable :: list(:)
      real (DP) :: phase, delta1, delta2, k0
      real (DP) :: xs, ys, zs, p1vect(3), p2vect(3)
      real (DP), parameter :: rad = pi/180._DP
!
!  determine material index for given element in region
      matindex=dommat(dom(elem))
!  read parameters from internal list and assign them to local variables
      allocate (list(numparam))
      xs = xyzs(1)
      ys = xyzs(2)
      zs = xyzs(3)
      call fetchmatparameters(list, matindex, xyzs)
!
!  possible choices for physics can be:
!  TEWAVE, TMWAVE
      select case (physics)
!
      case ('TEWAVE')
!  NO PML YET
        pml    = int(list(1))
!  if delta>0, we have losses
        delta1 = atan(list(2))
        delta2 = atan(list(3))
!  set phase for source
        phase  = rad*list(13)
!  set first and second principal axis direction
        p1vect = (/list(14),list(15),list(16)/)
        p2vect = (/list(17),list(18),list(19)/)
!
!  calculate tensor for diffusion term and assign result to nu

!  assign nu/permeability tensor (assume diagonal tensor)
        nu = 0._DPC
        nu(1,1) = 1._DP / list(4) * exp(cmplx(0._DP,delta1,DPC))
        nu(2,2) = 1._DP / list(5) * exp(cmplx(0._DP,delta1,DPC))
        nu(3,3) = 1._DP / list(6) * exp(cmplx(0._DP,delta1,DPC))
!  transform nu material tensor from principal axes system to x,y,z system
        call calctensor3D(nu,p1vect,p2vect)
!
!  assign alpha/permittivity tensor (assume diagonal tensor)
        alpha = 0._DPC
        alpha(1,1) = list(7) * exp(cmplx(0._DP,-delta2,DPC))
        alpha(2,2) = list(8) * exp(cmplx(0._DP,-delta2,DPC))
        alpha(3,3) = list(9) * exp(cmplx(0._DP,-delta2,DPC))
!  transform alpha material tensor from principal axes system to x,y,z system
        call calctensor3D(alpha,p1vect,p2vect)
!
        k0 = omega*sqrt(eps0*mu0)   ! free space wave number
        alpha = (k0**2)*alpha
!
!  we do not have a convective term => assign 0
        beta = 0._DPC
!
!  calculate complex values for current density f1
        f1(1) = list(10) * exp(cmplx(0._DP,phase,DPC))
        f1(2) = list(11) * exp(cmplx(0._DP,phase,DPC))
        f1(3) = list(12) * exp(cmplx(0._DP,phase,DPC))
!
        f1 = cmplx(0._DP,-omega*mu0,DPC)*f1
!  f2 is 0 for TEWAVE
        f2 = 0._DPC
!
!
      case ('TMWAVE')
!  NO PML YET
        pml    = int(list(1))
!  if delta>0, we have losses
        delta1 = atan(list(2))
        delta2 = atan(list(3))
!  set phase for source
        phase  = rad*list(13)
!  set first and second principal axis direction
        p1vect = (/list(14),list(15),list(16)/)
        p2vect = (/list(17),list(18),list(19)/)
!
!  calculate tensor for diffusion term and assign result to nu

!  assign nu/permittivity tensor (assume diagonal tensor)
        nu = 0._DPC
        nu(1,1) = 1._DP / list(7) * exp(cmplx(0._DP,delta2,DPC))
        nu(2,2) = 1._DP / list(8) * exp(cmplx(0._DP,delta2,DPC))
        nu(3,3) = 1._DP / list(9) * exp(cmplx(0._DP,delta2,DPC))
!  transform nu material tensor from principal axes system to x,y,z system
        call calctensor3D(nu,p1vect,p2vect)
!
!  assign alpha/permeability tensor (assume diagonal tensor)
        alpha = 0._DPC
        alpha(1,1) = list(4) * exp(cmplx(0._DP,-delta1,DPC))
        alpha(2,2) = list(5) * exp(cmplx(0._DP,-delta1,DPC))
        alpha(3,3) = list(6) * exp(cmplx(0._DP,-delta1,DPC))
!  transform alpha material tensor from principal axes system to x,y,z system
        call calctensor3D(alpha,p1vect,p2vect)
!
        k0 = omega*sqrt(eps0*mu0)   ! free space wave number
        alpha = (k0**2)*alpha
!
!  we do not have a convective term => assign 0
        beta = 0._DPC
!
!  f1 is 0 for TMWAVE
        f1 = 0._DPC
!  calculate complex values for current density f2
        f2(1) = list(10) * exp(cmplx(0._DP,phase,DPC))
        f2(2) = list(11) * exp(cmplx(0._DP,phase,DPC))
        f2(3) = list(12) * exp(cmplx(0._DP,phase,DPC))
!
        f2 = -matmul(nu,f2)
!
!
      case default

!  This default case is a normal wave equation [ -curl(curl(u)) + omega**2*u = 0 ]
        nu = 0._DPC
        nu(1,1) = cmplx(1._DP,0._DP,DPC)
        nu(2,2) = cmplx(1._DP,0._DP,DPC)
        nu(3,3) = cmplx(1._DP,0._DP,DPC)
        k0 = omega
        alpha = k0**2 * nu
        beta = 0._DPC
        f1 = 0._DPC
        f2 = 0._DPC
      end select
!
      deallocate(list)
!
      return
      end subroutine pdecoeff3D_vector


!--------------------------------------------------------------------------------------
!  For scalar-valued basis functions
!--------------------------------------------------------------------------------------
  
subroutine pdecoeff3D_scalar(elem, xyzs, nu, alpha, beta, f1, f2, gamma)
      use feminterface,      only: getsetting, fetchmatparameters, multtensor, trans_matrix, transformtensor
      use feminterface3d,    only: calctensor3D, xyz2lam, fieldsc3D_simple, fieldsc3D_simple_aux
      use femtypes
      use globalvariables3D, only: dom, dommat, omega, physics, pi, eps0, nod, vn, nnat, Elch, Kboltz, nnat_aux
      use matconstants
      implicit none
      integer (I4B) elem
      real(DP) :: xyzs(3)
      complex (DPC)             :: nu(:,:,:,:), alpha(:,:), beta(:,:,:), nu_mech(3,3,nnat_aux,nnat_aux)
      complex (DPC)             :: f1(:), f2(:,:)
      complex (DPC), optional   :: gamma(:,:,:)
      intent (in)               :: elem, xyzs
      intent (out)              :: nu, alpha, beta, f1, f2, gamma
!
!
!  (2) scalar potential:
!
!        -grad(nu*grad(U)) + (beta*grad(U)) + alpha*U = f1 + grad(f2)
!
!  Input:
!            elem     number of the element
!            x        actual solution vector [not used yet] (GLOBAL)
!            e        element information (nodes of the element) [not used yet] (GLOBAL)
!            matzif   material name of the input regions (GLOBAL)
!            geb      region index of the elements  (GLOBAL)
!            xyzs    location at which the material coefficients have to be evaluated
!                     this is usually the centroid of the volume element
!                     [not used yet]
!  Output:
!            nu       rank 2 tensor of the diffusion term:      { nu_xx  nu_xy  nu_xz }
!                                                          nu = | nu_yx  nu_yy  nu_yz |
!                                                               { nu_zx  nu_zy  nu_zz }
!            beta     vector of the convective term beta(1)=betax
!            alpha    rank 2 tensor, like nu 
!            f1       right hand side vector
!            f2       right hand side source vector
!
!  local variables
      integer (I4B)          :: matindex, pml, inat, jnat
      real( DP), allocatable :: list(:)
      real (DP)              :: lambda(4)
      complex (DPC)          :: z(15,nnat), u(nnat), gu(3,nnat),temp,ptnl, u_aux(nnat_aux), gu_aux(3,nnat_aux), stress_temp(nnat_aux), stress_tensor(3,3)
      real (DP)              :: phase, delta1, delta2, k0 ,var11, scalar, source1, source2
      real (DP)              :: xs, ys, zs, p1vect(3), p2vect(3), epsr(3)      
      real (DP)              :: Lame1, Lame2
      real (DP)              :: GenRate, N_A, N_D, X_N, X_P, tau_n, tau_p, n_i, mu_n, mu_p, T_ref, X_Psi, V_T, X_phi_n, X_phi_p, Phi_Bn, Phi_Bp
      real (DP)              :: mu_n_x, mu_n_y, mu_n_z, mu_p_x, mu_p_y, mu_p_z
      real (DP),   parameter :: rad = pi/180._DP
      real (DP)              :: CTensor(3,3,3,3), PITensor(3,3,3,3), CTensor_rot(3,3,3,3), PITensor_rot(3,3,3,3), delta_tensor(3,3)
      real (DP)              :: fpa(3), spa(3), t_matrix(3,3)
      real (DP)              :: pelt11, pelt22, pelt33      
      logical                :: typ(5)
!
!  determine material index for given element in region
      matindex=dommat(dom(elem))
!  read parameters from internal list and assign them to local variables
      allocate (list(numparam))
      xs = xyzs(1)
      ys = xyzs(2)
      zs = xyzs(3)
      call fetchmatparameters(list, matindex, xyzs)
!      print*,'size_of_list',size(list)
!
!  possible choices for physics can be:
!  TEWAVE, TMWAVE
      select case (physics)
!
case ('STATCURRENT')
        nu = 0._DPC
        nu(1,1,1,1) = cmplx(list(1),0._DP,DPC)
        nu(2,2,1,1) = cmplx(list(2),0._DP,DPC)
        nu(3,3,1,1) = cmplx(list(3),0._DP,DPC)
        f1 = 0._DPC
        f2 = 0._DPC
        alpha = 0._DPC
        beta = 0._DPC
        gamma = 0._DPC
!
case ('3DHEATTR')
!  This is the only scalar-valued problem for now [ -grad(grad(u)) + omega**2* u = 0 ]
        nu = 0._DPC
! list(1:3) contains thermal conductivity in all 3 principal axes as read from material files
!
        nu(1,1,1,1) = cmplx(list(1),0._DP,DPC)   !
        nu(2,2,1,1) = cmplx(list(2),0._DP,DPC)   !
        nu(3,3,1,1) = cmplx(list(3),0._DP,DPC)   !
        gamma = 0._DPC
        !k0 = 0._DP
        alpha = 0.0_DPC
        beta = 0._DPC
        f1 = list(7)
        f2 = 0._DPC

case ('3DNONLINHEATTR')
!  This is the only scalar-valued problem for now [ -grad(grad(u)) + omega**2* u = 0 ]

        !typ = (/ .true.,.false.,.false.,.false.,.false.  /)
        call xyz2lam(lambda, elem, xyzs, nod, vn)
        call fieldsc3D_simple(elem, lambda, u)
        nu = 0._DPC
! list(1:3) contains thermal conductivity in all 3 principal axes as read from material files
!
        nu(1,1,1,1) = (u(1)**2)*cmplx(list(1),0._DP,DPC)   !
        nu(2,2,1,1) = (u(1)**2)*cmplx(list(2),0._DP,DPC)   !
        nu(3,3,1,1) = (u(1)**2)*cmplx(list(3),0._DP,DPC)   !
        gamma = 0._DPC
        !k0 = 0._DP
        alpha = 0.0_DPC
        beta = 0._DPC
        f1 = list(7)
        f2 = 0._DPC

!
!########################################################
!########################################################
!####                                                ####
!####                 MULTIPHYSIC MODES              ####
!####                                                ####
!########################################################
!########################################################
!
case ('3DFLUIDINCOMPR')
!  nnat=4 in the Incompressible Fluid Problem in 3D
!  This is the only scalar-valued problem for now [ -grad(grad(u)) + omega**2* u = 0 ]
        nu = 0._DPC
! list(1:2) contains thermal conductivity in all 3 principal axes as read from material files
!
        nu(1,1,1,1) = cmplx(list(1),0._DP,DPC)   !
        nu(2,2,1,1) = cmplx(list(2),0._DP,DPC)   !
        nu(3,3,1,1) = cmplx(list(3),0._DP,DPC)   !
        gamma = 0._DPC
        k0 = 0._DP
        alpha = 0.0_DPC
        beta = 0._DPC
        f1 = list(10)
        f2 = 0._DPC
!
case ('3DTHERMOELECTRIC')
!  nnat=2 in the thermoelectric problem
!  the thermal is considered as the first nature
!  the electric is considered as the second nature
      
!        var11  = list(1)
!        var22  = list(2)
!        var33  = list(3)
!        angle  = rad*list(4)
! Source for the Thermal and Electric Equations : ELECTRO-THERMAL ANALYSIS OF PELTIER COOLING USING FEM
! D. ENESCU1, E.O. VÎRJOGHE2, M. IONEL1, M.F. STAN2 

! Thermal Equation -div.((K+alpha*alpha*T*sigma)gradT+ peltier_coeff*sigma*grad phi)= div.(q) 
! nu11 = K+alpha*alpha*T*sigma ; nu12 = peltier_coeff*sigma ; nu21 = alpha*sigma ; nu22 = sigma ; 
! alpha = 0 ; beta = 0 ; gamma = 0 ; f1 = q ; f2 = 0 ; u = T

! Electric Equation -div.(sigma*grad phi + sigma*alpha*grad T)= div.(J) 
! nu11 = K+alpha*alpha*T*sigma ; nu12 = peltier_coeff*sigma ; nu21 = alpha*sigma ; nu22 = sigma ; 
! alpha = 0 ; beta = 0 ; gamma = 0 ; f1 = J ; f2 = 0 ; u = phi

!  nu is of the form nu(3,3,nnat,nnat)
!  thermal nu = lamda/Tref

! To calculate the Peltier's constant

        typ = (/ .true.,.false.,.false.,.false.,.false.  /)
        call xyz2lam(lambda, elem, xyzs, nod, vn)
        call fieldsc3D_simple(elem, lambda, u, gu)
        
        temp   = u(1)  ! u is of nature 1 here (temp) 
        ptnl   = u(2)  ! u is of nature 2 here (potential)
        
        
        pelt11 = list(9)*temp
        pelt22 = list(10)*temp
        pelt33 = list(11)*temp

        nu(:,:,:,:)=0._DPC
 
        nu(1,1,1,1) = cmplx(list(3) + list(9)*list(6)*pelt11,0._DP,DPC)           ! kappa11
        nu(2,2,1,1) = cmplx(list(4) + list(10)*list(7)*pelt22,0._DP,DPC)          ! kappa22
        nu(3,3,1,1) = cmplx(list(5) + list(11)*list(8)*pelt33,0._DP,DPC)          ! kappa33
        
        nu(1,1,2,2) = cmplx(list(6),0._DP,DPC)                                    ! sigma11
        nu(2,2,2,2) = cmplx(list(7),0._DP,DPC)                                    ! sigma22
        nu(3,3,2,2) = cmplx(list(8),0._DP,DPC)                                    ! sigma33
        
        nu(1,1,1,2) = cmplx(list(6)*pelt11,0._DP,DPC)                             ! pelt11
        nu(2,2,1,2) = cmplx(list(7)*pelt22,0._DP,DPC)                             ! pelt22
        nu(3,3,1,2) = cmplx(list(8)*pelt33,0._DP,DPC)                             ! pelt33
        
        nu(1,1,2,1) = cmplx(list(6)*list(9),0._DP,DPC)                            ! seeb11
        nu(2,2,2,1) = cmplx(list(7)*list(10),0._DP,DPC)                           ! seeb22 
        nu(3,3,2,1) = cmplx(list(8)*list(11),0._DP,DPC)                           ! seeb33

!  alpha is of the form alpha(3,3,nnat,nnat)
        alpha = 0._DPC
         
!  beta is of the form beta(3,nnat)
        beta = 0._DPC         
  
!  gamma is of the form gamma(3,nnat,nnat)
        gamma(:,:,:) = 0._DPC         
        
!  f1 is of the form f1(nnat)
        f1(1) = list(6)*(gu(1,1)*list(9) + gu(1,2))*(gu(1,2)) + list(7)*(gu(2,1)*list(10) + gu(2,2))*(gu(2,2)) + list(8)*(gu(3,1)*list(11) + gu(3,2))*(gu(3,2))
        f1(2) = 0._DPC
!  f2 is of the form f2(nnat)        
        f2 = 0._DPC 
    

      case ('3DELECTROSTATICS')

!  This case is for a Gaussian Equation (Electrostatics)
        nu = 0._DPC
        nu(1,1,1,1) = cmplx(list(1),0._DP,DPC) ! epsr11
        nu(2,2,1,1) = cmplx(list(2),0._DP,DPC) ! epsr22
        nu(3,3,1,1) = cmplx(list(3),0._DP,DPC) ! epsr33
        alpha = 0._DPC
        beta = 0._DPC
        gamma = 0._DPC
        f1 = list(4)/eps0                       ! rho/eps0 (source/charge)
        f2 = 0._DPC
        
      case ('3DPOISSON')
        ! nnat = 1
        ! nat = 1 : Electric Potential (phi)
        ! Poisson Eq:     -div( eps/e0 * grad(phi) ) = N_charge
        ! phi      = el. potential (unknown)
        ! N_charge = number of elementary charges 1/m^3 causing potential gradient

        call xyz2lam(lambda, elem, xyzs, nod, vn)        
        
        call fieldsc3D_simple(elem, lambda, u)
        X_Psi   = u(1)

        epsr(1) = list(1)
        epsr(2) = list(2) 
        epsr(3) = list(3)        
        N_A     = list(4)
        N_D     = list(5)        
        T_ref   = list(6)
        mu_n    = list(7)
        mu_p    = list(8)
        n_i     = list(9)
        tau_n   = list(10)
        tau_p   = list(11)  

        V_T     = KBoltz*T_ref/Elch 

        ! nu is of the form nu(3,3,nnat,nnat)
        nu(:,:,:,:) = 0._DPC
        ! nu_11
        nu(1,1,1,1) = cmplx(eps0*epsr(1),0._DP,DPC)
        nu(2,2,1,1) = cmplx(eps0*epsr(2),0._DP,DPC)
        nu(3,3,1,1) = cmplx(eps0*epsr(3),0._DP,DPC)      
        !
        ! alpha is of the form alpha(nnat,nnat)
        alpha(:,:) = 0._DPC
        beta(:,:,:) = 0._DPC
        gamma(:,:,:) = 0._DPC

        f1(1) = cmplx(Elch*(-2*n_i*sinh(X_Psi/V_T) + (N_D - N_A)), 0._DP,DPC) 
        f2(:,:) = 0._DPC


      case ('3DSEMICONDUCTOR')
        ! 3D Semiconductor Mode using the equivalent-Fermi-voltage formulation:
        ! nnat = 3
        ! nat = 1 : Electric Potential (Psi)
        ! nat = 2 : Quasi-fermi-level electrons (phi_n)
        ! nat = 3 : Quasi-fermi-level holes (phi_p)       
        ! Poisson Eq:     -div( eps * grad(Psi*) ) = q* (n_i*exp(-Psi/V_T)*exp(phi_p/V_T) - n_i*exp(Psi/V_T)*exp(-phi_n/V_T) + N_D - N_A)
        ! Contiuity Eqs:  -div( mu_n * n_i * exp( (Psi-phi_n)/V_T)* grad(phi_n) ) = G
        !                  div( mu_p * n_i * exp(-(Psi-phi_p)/V_T)* grad(phi_p) ) = G

        epsr(1) = list(1)
        epsr(2) = list(2)
        epsr(3) = list(3)
        N_A     = list(4)
        N_D     = list(5)
        T_ref   = list(6)
        mu_n    = list(7)
        mu_p    = list(8)
        n_i     = list(9)
        tau_n   = list(10)
        tau_p   = list(11)
        
        call xyz2lam(lambda, elem, xyzs, nod, vn)
        call fieldsc3D_simple(elem, lambda, u)
        
        X_Psi    = u(1)
        X_phi_n  = u(2)
        X_phi_p  = u(3)
        


        V_T     = KBoltz*T_ref/Elch
        
        !X_N = n_i * exp((Phi_Bn)/V_T)
        !X_P = n_i * exp((Phi_Bp)/V_T)
        !
        !if ( abs(X_Psi).gt.1._DP ) then
        !  X_Psi = ((Kboltz*T_ref)/Elch)*asinh((N_D-N_A)/(2*n_i))
        !end if
        
        
        Phi_Bn = X_Psi-X_phi_n
        Phi_Bp = X_phi_p-X_Psi
        
        
        if ( abs(Phi_Bn).gt.11._DP ) then
          X_N = n_i * (1._DP + (Phi_Bn)/V_T + ((Phi_Bn)/V_T)**2/2._DP)
        else
          X_N = n_i * exp((Phi_Bn)/V_T)
        end if
        
        if ( abs(Phi_Bp).gt.11._DP ) then
          X_P = n_i * (1._DP + (Phi_Bp)/V_T + ((Phi_Bp)/V_T)**2/2._DP)
        else
          X_P = n_i * exp((Phi_Bp)/V_T)
        end if
 
        
        ! Calculate Generation/Recombination rate term
        !GenRate = 0._DP
        if (abs(X_phi_p-X_phi_n) .lt. 10*epsilon(X_phi_n)) then
           Genrate = 0._DP         
        else
           GenRate = (X_N*X_P - n_i**2) / (tau_n*(X_P+n_i)+tau_p*(X_N+n_i))
        end if
         
        ! nu is of the form nu(3,3,nnat,nnat)
        nu(:,:,:,:) = 0._DPC
        ! nu_11
        nu(1,1,1,1) = cmplx(eps0*epsr(1),0._DP,DPC)
        nu(2,2,1,1) = cmplx(eps0*epsr(2),0._DP,DPC)
        nu(3,3,1,1) = cmplx(eps0*epsr(3),0._DP,DPC)
        
        ! nu_22
        nu(1,1,2,2) = cmplx(mu_n * X_N,0._DP,DPC)
        nu(2,2,2,2) = cmplx(mu_n * X_N,0._DP,DPC)
        nu(3,3,2,2) = cmplx(mu_n * X_N,0._DP,DPC)
        
        ! nu_33 ! Continuity Eq multiplied by (-1)
        nu(1,1,3,3) = cmplx(mu_p * X_P,0._DP,DPC)
        nu(2,2,3,3) = cmplx(mu_p * X_P,0._DP,DPC)
        nu(3,3,3,3) = cmplx(mu_p * X_P,0._DP,DPC)

        if (any(real(nu(:,:,:,:)).lt.0)) then
          print*,nu
          print*,'Element', elem
        end if
        
        ! alpha is of the form alpha(nnat,nnat)
        alpha(:,:) = 0._DPC        
        beta(:,:,:) = 0._DPC
        gamma(:,:,:) = 0._DPC
        
        f1(1) = cmplx(Elch*(X_P - X_N + N_D - N_A), 0._DP,DPC)
        f1(2) = cmplx(GenRate, 0._DP,DPC)
        f1(3) = cmplx(-GenRate, 0._DP,DPC) ! Continuity Eq multiplied by (-1)
        
        f2(:,:) = 0._DPC

      case ('3DSOLIDMECHANICS')
!  nnat=3 in the 3D Strain problem
!  in the first nature we consider the displacement solution in x
!  in the second nature we consider the displacement solution in y
!  in the third nature we consider the displacement solution in z
!  First calculate Lame coefficients 1 and 2 (lamda and mu) from
!  given Young's modulus (denoted here with E)
!  and Poisson's ratio (denoted here with v)
!  The formulae are: Lame1 = E*v/((1+v)*(1-2v))
!                    Lame2 = E/(2*(1+v))
!-----------------------------------------------------------------
! List of Parameters:
!
! Info: for anisotropic materials the stiffness matrix coefficients
!  |sigma_xx|   |C11 C12 C13 C14 C15 C16|   |epsilon_xx|
!  |sigma_yy|   |C21 C22 C23 C24 C25 C26|   |epsilon_yy|
!  |sigma_zz|   |C31 C32 C33 C34 C35 C36|   |epsilon_zz|
!  | tau_yz | = |C41 C42 C43 C44 C45 C46| * | gamma_yz |
!  | tau_xz |   |C51 C52 C53 C54 C55 C56|   | gamma_xz |
!  | tau_xy |   |C61 C62 C63 C64 C65 C66|   | gamma_xy |
!
! For the indication of an anisotropic material
! default value is 0 so this means that the material is isotropic by default
!
!   list(1)     = 'rho', material density in kg/m3
!   list(2)     = 'Young', Young's modulus E in Pascal
!   list(3)     = 'Poisson', Poisson's ratio nu
!   list(4-6)   = x,y,z - component of the first principle axis
!   list(7-9)   = x,y,z - component of the seccond principle axis
!   list(10-30) = stiffness matrix coefficients
!   list(31)    = indication of an anisotropic material
!   list(32-34) = Body (Volume) force in x,y,z - direction
!              E * nu /(1+nu)/(1-2nu)
        Lame1 = ( list(2)*list(3) ) / ( (1._DP + list(3))*(1._DP - 2._DP*list(3)) )
!              E /(1+nu)/2
        Lame2 = list(2) / ( 1._DP+list(3) ) /2._DP

! for the case of an isotropic material
        if ( list(31) .eq. 0._DP ) then
!  nu is of the form nu(3,3,nnat,nnat)
!  submatrix 1-1, x <- x
! cmplx(Lame1+2._DP*Lame2,0._DP,DPC)
! cmplx(Lame1,0._DP,DPC)
! cmplx(Lame2,0._DP,DPC)
          nu(:,:,:,:) = 0._DP
          nu(1,1,1,1) = cmplx(Lame1+2._DP*Lame2,0._DP,DPC)
          nu(2,2,1,1) = cmplx(Lame2,0._DP,DPC)
          nu(3,3,1,1) = cmplx(Lame2,0._DP,DPC)

          nu(1,2,1,2) = cmplx(Lame1,0._DP,DPC) !Lam
          nu(2,1,1,2) = cmplx(Lame2,0._DP,DPC) !mu

          nu(1,3,1,3) = cmplx(Lame1,0._DP,DPC) !Lam
          nu(3,1,1,3) = cmplx(Lame2,0._DP,DPC) !mu

          nu(1,2,2,1) = cmplx(Lame2,0._DP,DPC) !mu
          nu(2,1,2,1) = cmplx(Lame1,0._DP,DPC) !Lam

          nu(1,1,2,2) = cmplx(Lame2,0._DP,DPC)
          nu(2,2,2,2) = cmplx(Lame1+2._DP*Lame2,0._DP,DPC)
          nu(3,3,2,2) = cmplx(Lame2,0._DP,DPC)

          nu(2,3,2,3) = cmplx(Lame1,0._DP,DPC) !Lam
          nu(3,2,2,3) = cmplx(Lame2,0._DP,DPC) !mu

          nu(1,3,3,1) = cmplx(Lame2,0._DP,DPC) !mu
          nu(3,1,3,1) = cmplx(Lame1,0._DP,DPC) !Lam

          nu(2,3,3,2) = cmplx(Lame2,0._DP,DPC) !mu
          nu(3,2,3,2) = cmplx(Lame1,0._DP,DPC) !Lam

          nu(1,1,3,3) = cmplx(Lame2,0._DP,DPC)
          nu(2,2,3,3) = cmplx(Lame2,0._DP,DPC)
          nu(3,3,3,3) = cmplx(Lame1+2._DP*Lame2,0._DP,DPC)
        else ! Non-Isotropic Material
          
          ! Get prinicpal axes
          fpa    = (/ list(4), list(5), list(6) /)    ! fpa_x, fpa_y, fpa_z
          spa    = (/ list(7), list(8), list(9) /)    ! spa_x, spa_y, spa_z
          
          ! Get stiffness tensor
          CTensor(:,:,:,:) = 0._DP
          CTensor(1,1,1,1) = list(10)    ! C11
          CTensor(1,1,2,2) = list(11)    ! C12
          CTensor(1,1,3,3) = list(12)    ! C13
          CTensor(1,1,2,3) = list(13)    ! C14
          CTensor(1,1,3,1) = list(14)    ! C15
          CTensor(1,1,1,2) = list(15)    ! C16
          CTensor(2,2,2,2) = list(16)    ! C22
          CTensor(2,2,3,3) = list(17)    ! C23
          CTensor(2,2,2,3) = list(18)    ! C24
          CTensor(2,2,3,1) = list(19)    ! C25
          CTensor(2,2,1,2) = list(20)    ! C26
          CTensor(3,3,3,3) = list(21)    ! C33
          CTensor(3,3,2,3) = list(22)    ! C34
          CTensor(3,3,3,1) = list(23)    ! C35
          CTensor(3,3,1,1) = list(24)    ! C36
          CTensor(2,3,2,3) = list(25)    ! C44
          CTensor(2,3,3,1) = list(26)    ! C45
          CTensor(2,3,1,2) = list(27)    ! C46
          CTensor(3,1,3,1) = list(28)    ! C55
          CTensor(3,1,1,2) = list(29)    ! C56
          CTensor(1,2,1,2) = list(30)    ! C66
          
          ! Infer other half of tensor by symmetry:
          CTensor(2,2,1,1) = CTensor(1,1,2,2)
          CTensor(3,3,1,1) = CTensor(1,1,3,3)
          CTensor(2,3,1,1) = CTensor(1,1,2,3)
          CTensor(3,3,2,2) = CTensor(2,2,3,3)
          CTensor(3,1,1,1) = CTensor(1,1,3,1)
          CTensor(2,3,2,2) = CTensor(2,2,2,3)
          CTensor(1,2,1,1) = CTensor(1,1,1,2)
          CTensor(3,1,2,2) = CTensor(2,2,3,1)
          CTensor(2,3,3,3) = CTensor(3,3,2,3)
          CTensor(1,2,2,2) = CTensor(2,2,1,2)
          CTensor(3,1,3,3) = CTensor(3,3,3,1)
          CTensor(1,2,3,3) = CTensor(3,3,1,2)
          CTensor(3,1,2,3) = CTensor(2,3,3,1)
          CTensor(1,2,2,3) = CTensor(2,3,1,2)
          CTensor(1,2,3,1) = CTensor(3,1,1,2)
          
          ! Transform stiffness tensor
          t_matrix = trans_matrix(fpa, spa)
          call transformtensor(t_matrix, CTensor, CTensor_rot)      
          
          !  nu is of the form nu(3,3,nnat,nnat)
          nu(:,:,:,:) = 0._DPC
          nu(1,1,1,1) = cmplx(CTensor_rot(1,1,1,1), 0._DP,DPC) ! C'_1111
          nu(1,2,1,1) = cmplx(CTensor_rot(1,1,1,2), 0._DP,DPC) ! C'_1112
          nu(1,3,1,1) = cmplx(CTensor_rot(1,1,3,1), 0._DP,DPC) ! C'_1131
          nu(2,1,1,1) = cmplx(CTensor_rot(1,1,1,2), 0._DP,DPC) ! C'_1112
          nu(2,2,1,1) = cmplx(CTensor_rot(1,2,1,2), 0._DP,DPC) ! C'_1212
          nu(2,3,1,1) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu(3,1,1,1) = cmplx(CTensor_rot(1,1,3,1), 0._DP,DPC) ! C'_1131
          nu(3,2,1,1) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu(3,3,1,1) = cmplx(CTensor_rot(3,1,3,1), 0._DP,DPC) ! C'_3131
          
          nu(1,1,1,2) = cmplx(CTensor_rot(1,1,1,2), 0._DP,DPC) ! C'_1112
          nu(1,2,1,2) = cmplx(CTensor_rot(1,1,2,2), 0._DP,DPC) ! C'_1122
          nu(1,3,1,2) = cmplx(CTensor_rot(1,1,2,3), 0._DP,DPC) ! C'_1123
          nu(2,1,1,2) = cmplx(CTensor_rot(1,2,1,2), 0._DP,DPC) ! C'_1212
          nu(2,2,1,2) = cmplx(CTensor_rot(2,2,1,2), 0._DP,DPC) ! C'_2212
          nu(2,3,1,2) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu(3,1,1,2) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu(3,2,1,2) = cmplx(CTensor_rot(2,2,3,1), 0._DP,DPC) ! C'_2231
          nu(3,3,1,2) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          
          nu(1,1,1,3) = cmplx(CTensor_rot(1,1,3,1), 0._DP,DPC) ! C'_1131
          nu(1,2,1,3) = cmplx(CTensor_rot(1,1,2,3), 0._DP,DPC) ! C'_1123
          nu(1,3,1,3) = cmplx(CTensor_rot(1,1,3,3), 0._DP,DPC) ! C'_1133
          nu(2,1,1,3) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu(2,2,1,3) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu(2,3,1,3) = cmplx(CTensor_rot(3,3,1,2), 0._DP,DPC) ! C'_3312
          nu(3,1,1,3) = cmplx(CTensor_rot(3,1,3,1), 0._DP,DPC) ! C'_3131
          nu(3,2,1,3) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu(3,3,1,3) = cmplx(CTensor_rot(3,3,3,1), 0._DP,DPC) ! C'_3331
          
          nu(1,1,2,1) = cmplx(CTensor_rot(1,1,1,2), 0._DP,DPC) ! C'_1112
          nu(1,2,2,1) = cmplx(CTensor_rot(1,2,1,2), 0._DP,DPC) ! C'_1212
          nu(1,3,2,1) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu(2,1,2,1) = cmplx(CTensor_rot(1,1,2,2), 0._DP,DPC) ! C'_1122
          nu(2,2,2,1) = cmplx(CTensor_rot(2,2,1,2), 0._DP,DPC) ! C'_2212
          nu(2,3,2,1) = cmplx(CTensor_rot(2,2,3,1), 0._DP,DPC) ! C'_2231
          nu(3,1,2,1) = cmplx(CTensor_rot(1,1,2,3), 0._DP,DPC) ! C'_1123
          nu(3,2,2,1) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu(3,3,2,1) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
           
          nu(1,1,2,2) = cmplx(CTensor_rot(1,2,1,2), 0._DP,DPC) ! C'_1212
          nu(1,2,2,2) = cmplx(CTensor_rot(2,2,1,2), 0._DP,DPC) ! C'_2212
          nu(1,3,2,2) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu(2,1,2,2) = cmplx(CTensor_rot(2,2,1,2), 0._DP,DPC) ! C'_2212
          nu(2,2,2,2) = cmplx(CTensor_rot(2,2,2,2), 0._DP,DPC) ! C'_2222
          nu(2,3,2,2) = cmplx(CTensor_rot(2,2,2,3), 0._DP,DPC) ! C'_2223
          nu(3,1,2,2) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu(3,2,2,2) = cmplx(CTensor_rot(2,2,2,3), 0._DP,DPC) ! C'_2223
          nu(3,3,2,2) = cmplx(CTensor_rot(2,3,2,3), 0._DP,DPC) ! C'_2323
           
          nu(1,1,2,3) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu(1,2,2,3) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu(1,3,2,3) = cmplx(CTensor_rot(3,3,1,2), 0._DP,DPC) ! C'_3312
          nu(2,1,2,3) = cmplx(CTensor_rot(2,2,3,1), 0._DP,DPC) ! C'_2231
          nu(2,2,2,3) = cmplx(CTensor_rot(2,2,2,3), 0._DP,DPC) ! C'_2223
          nu(2,3,2,3) = cmplx(CTensor_rot(2,2,3,3), 0._DP,DPC) ! C'_2233
          nu(3,1,2,3) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu(3,2,2,3) = cmplx(CTensor_rot(2,3,2,3), 0._DP,DPC) ! C'_2323
          nu(3,3,2,3) = cmplx(CTensor_rot(3,3,2,3), 0._DP,DPC) ! C'_3323
          
          nu(1,1,3,1) = cmplx(CTensor_rot(1,1,3,1), 0._DP,DPC) ! C'_1131
          nu(1,2,3,1) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu(1,3,3,1) = cmplx(CTensor_rot(3,1,3,1), 0._DP,DPC) ! C'_3131
          nu(2,1,3,1) = cmplx(CTensor_rot(1,1,2,3), 0._DP,DPC) ! C'_1123
          nu(2,2,3,1) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu(2,3,3,1) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu(3,1,3,1) = cmplx(CTensor_rot(1,1,3,3), 0._DP,DPC) ! C'_1133
          nu(3,2,3,1) = cmplx(CTensor_rot(3,3,1,2), 0._DP,DPC) ! C'_3312
          nu(3,3,3,1) = cmplx(CTensor_rot(3,3,3,1), 0._DP,DPC) ! C'_3331
          
          nu(1,1,3,2) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu(1,2,3,2) = cmplx(CTensor_rot(2,2,3,1), 0._DP,DPC) ! C'_2231
          nu(1,3,3,2) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu(2,1,3,2) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu(2,2,3,2) = cmplx(CTensor_rot(2,2,2,3), 0._DP,DPC) ! C'_2223
          nu(2,3,3,2) = cmplx(CTensor_rot(2,3,2,3), 0._DP,DPC) ! C'_2323
          nu(3,1,3,2) = cmplx(CTensor_rot(3,3,1,2), 0._DP,DPC) ! C'_3312
          nu(3,2,3,2) = cmplx(CTensor_rot(2,2,3,3), 0._DP,DPC) ! C'_2233
          nu(3,3,3,2) = cmplx(CTensor_rot(3,3,2,3), 0._DP,DPC) ! C'_3323
           
          nu(1,1,3,3) = cmplx(CTensor_rot(3,1,3,1), 0._DP,DPC) ! C'_3131
          nu(1,2,3,3) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu(1,3,3,3) = cmplx(CTensor_rot(3,3,3,1), 0._DP,DPC) ! C'_3331
          nu(2,1,3,3) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu(2,2,3,3) = cmplx(CTensor_rot(2,3,2,3), 0._DP,DPC) ! C'_2323
          nu(2,3,3,3) = cmplx(CTensor_rot(3,3,2,3), 0._DP,DPC) ! C'_3323
          nu(3,1,3,3) = cmplx(CTensor_rot(3,3,3,1), 0._DP,DPC) ! C'_3331
          nu(3,2,3,3) = cmplx(CTensor_rot(3,3,2,3), 0._DP,DPC) ! C'_3323
          nu(3,3,3,3) = cmplx(CTensor_rot(3,3,3,3), 0._DP,DPC) ! C'_3333
        end if
!
!  gamma is of the form gamma(3,nnat,nnat)
        gamma(:,:,:) = 0._DPC
!
!  alpha is of the form alpha(nnat,nnat)
        alpha(1,1) = cmplx(-((omega**2)*list(1)),0._DP,DPC)
        alpha(1,2) = 0._DPC
        alpha(1,3) = 0._DPC
        alpha(2,1) = 0._DPC
        alpha(2,2) = cmplx(-((omega**2)*list(1)),0._DP,DPC)
        alpha(2,3) = 0._DPC
        alpha(3,1) = 0._DPC
        alpha(3,2) = 0._DPC
        alpha(3,3) = cmplx(-((omega**2)*list(1)),0._DP,DPC)
!  beta is of the form beta(3,nnat,nnat)
        beta(:,:,:) = 0._DPC
!  f1 is of the form f1(nnat)
        f1(1) = cmplx(list(32),0._DP,DPC)
        f1(2) = cmplx(list(33),0._DP,DPC)
        f1(3) = cmplx(list(34),0._DP,DPC)
!  f2 is of the form f2(3,nnat)
        f2(:,:) = 0._DPC

        
      case ('3DMECHANICSSEMICOND')
        ! 3D Semiconductor Mode using the equivalent-Fermi-voltage formulation,
        ! taking into account an existing mechanical solution (from 3DSOLIDMECHANICS)
        ! nnat = 3
        ! nat = 1 : Electric Potential (Psi)
        ! nat = 2 : Quasi-fermi-level electrons (phi_n)
        ! nat = 3 : Quasi-fermi-level holes (phi_p)       
        ! Poisson Eq:     -div( eps * grad(Psi*) ) = q* (n_i*exp(-Psi/V_T)*exp(phi_p/V_T) - n_i*exp(Psi/V_T)*exp(-phi_n/V_T) + N_D - N_A)
        ! Contiuity Eqs:  -div( mu_n * n_i * exp( (Psi-phi_n)/V_T)* grad(phi_n) ) = G
        !                  div( mu_p * n_i * exp(-(Psi-phi_p)/V_T)* grad(phi_p) ) = G

        
        ! initialize tensors
        CTensor(:,:,:,:)    = 0._DP
        PITensor(:,:,:,:)   = 0._DP
        stress_tensor(:,:)  = 0._DPC
        delta_tensor(:,:)      = 0._DP
        
        
        epsr(1) = list(1)
        epsr(2) = list(2)
        epsr(3) = list(3)
        N_A     = list(4)
        N_D     = list(5)
        T_ref   = list(6)
        mu_n    = list(7)
        mu_p    = list(8)
        n_i     = list(9)
        tau_n   = list(10)
        tau_p   = list(11)
        CTensor(1,1,1,1) = list(12)    ! C11
        CTensor(1,1,2,2) = list(13)    ! C12
        CTensor(1,1,3,3) = list(14)    ! C13
        CTensor(1,1,2,3) = list(15)    ! C14
        CTensor(1,1,3,1) = list(16)    ! C15
        CTensor(1,1,1,2) = list(17)    ! C16
        CTensor(2,2,2,2) = list(18)    ! C22
        CTensor(2,2,3,3) = list(19)    ! C23
        CTensor(2,2,2,3) = list(20)    ! C24
        CTensor(2,2,3,1) = list(21)    ! C25
        CTensor(2,2,1,2) = list(22)    ! C26
        CTensor(3,3,3,3) = list(23)    ! C33
        CTensor(3,3,2,3) = list(24)    ! C34
        CTensor(3,3,3,1) = list(25)    ! C35
        CTensor(3,3,1,1) = list(26)    ! C36
        CTensor(2,3,2,3) = list(27)    ! C44
        CTensor(2,3,3,1) = list(28)    ! C45
        CTensor(2,3,1,2) = list(29)    ! C46
        CTensor(3,1,3,1) = list(30)    ! C55
        CTensor(3,1,1,2) = list(31)    ! C56
        CTensor(1,2,1,2) = list(32)    ! C66
        
        ! Infer lower half of matrix by symmetry:
        CTensor(2,2,1,1) = CTensor(1,1,2,2)
        CTensor(3,3,1,1) = CTensor(1,1,3,3)
        CTensor(2,3,1,1) = CTensor(1,1,2,3)
        CTensor(3,3,2,2) = CTensor(2,2,3,3)
        CTensor(3,1,1,1) = CTensor(1,1,3,1)
        CTensor(2,3,2,2) = CTensor(2,2,2,3)
        CTensor(1,2,1,1) = CTensor(1,1,1,2)
        CTensor(3,1,2,2) = CTensor(2,2,3,1)
        CTensor(2,3,3,3) = CTensor(3,3,2,3)
        CTensor(1,2,2,2) = CTensor(2,2,1,2)
        CTensor(3,1,3,3) = CTensor(3,3,3,1)
        CTensor(1,2,3,3) = CTensor(3,3,1,2)
        CTensor(3,1,2,3) = CTensor(2,3,3,1)
        CTensor(1,2,2,3) = CTensor(2,3,1,2)
        CTensor(1,2,3,1) = CTensor(3,1,1,2)
        
        PITensor(1,1,1,1) = list(33)    ! PI11
        PITensor(2,2,2,2) = list(33)
        PITensor(3,3,3,3) = list(33)
        
        PITensor(1,1,2,2) = list(34)    ! PI12
        PITensor(1,1,3,3) = list(34)
        PITensor(2,2,3,3) = list(34)
        PITensor(2,2,1,1) = list(34)
        PITensor(3,3,1,1) = list(34)
        PITensor(3,3,2,2) = list(34)
        
        PITensor(2,3,2,3) = list(35)    ! PI44
        PITensor(3,1,3,1) = list(35)
        PITensor(3,1,3,1) = list(35)
        
        fpa    = (/ list(36), list(37), list(38) /)    ! fpa_x, fpa_y, fpa_z
        spa    = (/ list(39), list(40), list(41) /)    ! spa_x, spa_y, spa_z
        
        ! Get transformation matrix 
        t_matrix = trans_matrix(fpa, spa)
        
        ! Transform tensors according to crystal orientation
        call transformtensor(t_matrix, CTensor, CTensor_rot)
        call transformtensor(t_matrix, PITensor, PITensor_rot)        
        
        call xyz2lam(lambda, elem, xyzs, nod, vn)
        call fieldsc3D_simple(elem, lambda, u)
        call fieldsc3D_simple_aux(elem, lambda, u_aux, gu_aux)
        
        X_Psi    = u(1)
        X_phi_n  = u(2)
        X_phi_p  = u(3)
        
        ! Calculate hole/electron mobilities from given mechanical solution
        ! TODO: check if u_aux is valid
        
        ! Replicate nu from 3DSOLIDMECHANICS (needed for calculating the stress tensor)
          nu_mech(:,:,:,:) = 0._DPC
          nu_mech(1,1,1,1) = cmplx(CTensor_rot(1,1,1,1), 0._DP,DPC) ! C'_1111
          nu_mech(1,2,1,1) = cmplx(CTensor_rot(1,1,1,2), 0._DP,DPC) ! C'_1112
          nu_mech(1,3,1,1) = cmplx(CTensor_rot(1,1,3,1), 0._DP,DPC) ! C'_1131
          nu_mech(2,1,1,1) = cmplx(CTensor_rot(1,1,1,2), 0._DP,DPC) ! C'_1112
          nu_mech(2,2,1,1) = cmplx(CTensor_rot(1,2,1,2), 0._DP,DPC) ! C'_1212
          nu_mech(2,3,1,1) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu_mech(3,1,1,1) = cmplx(CTensor_rot(1,1,3,1), 0._DP,DPC) ! C'_1131
          nu_mech(3,2,1,1) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu_mech(3,3,1,1) = cmplx(CTensor_rot(3,1,3,1), 0._DP,DPC) ! C'_3131
          
          nu_mech(1,1,1,2) = cmplx(CTensor_rot(1,1,1,2), 0._DP,DPC) ! C'_1112
          nu_mech(1,2,1,2) = cmplx(CTensor_rot(1,1,2,2), 0._DP,DPC) ! C'_1122
          nu_mech(1,3,1,2) = cmplx(CTensor_rot(1,1,2,3), 0._DP,DPC) ! C'_1123
          nu_mech(2,1,1,2) = cmplx(CTensor_rot(1,2,1,2), 0._DP,DPC) ! C'_1212
          nu_mech(2,2,1,2) = cmplx(CTensor_rot(2,2,1,2), 0._DP,DPC) ! C'_2212
          nu_mech(2,3,1,2) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu_mech(3,1,1,2) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu_mech(3,2,1,2) = cmplx(CTensor_rot(2,2,3,1), 0._DP,DPC) ! C'_2231
          nu_mech(3,3,1,2) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          
          nu_mech(1,1,1,3) = cmplx(CTensor_rot(1,1,3,1), 0._DP,DPC) ! C'_1131
          nu_mech(1,2,1,3) = cmplx(CTensor_rot(1,1,2,3), 0._DP,DPC) ! C'_1123
          nu_mech(1,3,1,3) = cmplx(CTensor_rot(1,1,3,3), 0._DP,DPC) ! C'_1133
          nu_mech(2,1,1,3) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu_mech(2,2,1,3) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu_mech(2,3,1,3) = cmplx(CTensor_rot(3,3,1,2), 0._DP,DPC) ! C'_3312
          nu_mech(3,1,1,3) = cmplx(CTensor_rot(3,1,3,1), 0._DP,DPC) ! C'_3131
          nu_mech(3,2,1,3) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu_mech(3,3,1,3) = cmplx(CTensor_rot(3,3,3,1), 0._DP,DPC) ! C'_3331
          
          nu_mech(1,1,2,1) = cmplx(CTensor_rot(1,1,1,2), 0._DP,DPC) ! C'_1112
          nu_mech(1,2,2,1) = cmplx(CTensor_rot(1,2,1,2), 0._DP,DPC) ! C'_1212
          nu_mech(1,3,2,1) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu_mech(2,1,2,1) = cmplx(CTensor_rot(1,1,2,2), 0._DP,DPC) ! C'_1122
          nu_mech(2,2,2,1) = cmplx(CTensor_rot(2,2,1,2), 0._DP,DPC) ! C'_2212
          nu_mech(2,3,2,1) = cmplx(CTensor_rot(2,2,3,1), 0._DP,DPC) ! C'_2231
          nu_mech(3,1,2,1) = cmplx(CTensor_rot(1,1,2,3), 0._DP,DPC) ! C'_1123
          nu_mech(3,2,2,1) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu_mech(3,3,2,1) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
           
          nu_mech(1,1,2,2) = cmplx(CTensor_rot(1,2,1,2), 0._DP,DPC) ! C'_1212
          nu_mech(1,2,2,2) = cmplx(CTensor_rot(2,2,1,2), 0._DP,DPC) ! C'_2212
          nu_mech(1,3,2,2) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu_mech(2,1,2,2) = cmplx(CTensor_rot(2,2,1,2), 0._DP,DPC) ! C'_2212
          nu_mech(2,2,2,2) = cmplx(CTensor_rot(2,2,2,2), 0._DP,DPC) ! C'_2222
          nu_mech(2,3,2,2) = cmplx(CTensor_rot(2,2,2,3), 0._DP,DPC) ! C'_2223
          nu_mech(3,1,2,2) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu_mech(3,2,2,2) = cmplx(CTensor_rot(2,2,2,3), 0._DP,DPC) ! C'_2223
          nu_mech(3,3,2,2) = cmplx(CTensor_rot(2,3,2,3), 0._DP,DPC) ! C'_2323
           
          nu_mech(1,1,2,3) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu_mech(1,2,2,3) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu_mech(1,3,2,3) = cmplx(CTensor_rot(3,3,1,2), 0._DP,DPC) ! C'_3312
          nu_mech(2,1,2,3) = cmplx(CTensor_rot(2,2,3,1), 0._DP,DPC) ! C'_2231
          nu_mech(2,2,2,3) = cmplx(CTensor_rot(2,2,2,3), 0._DP,DPC) ! C'_2223
          nu_mech(2,3,2,3) = cmplx(CTensor_rot(2,2,3,3), 0._DP,DPC) ! C'_2233
          nu_mech(3,1,2,3) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu_mech(3,2,2,3) = cmplx(CTensor_rot(2,3,2,3), 0._DP,DPC) ! C'_2323
          nu_mech(3,3,2,3) = cmplx(CTensor_rot(3,3,2,3), 0._DP,DPC) ! C'_3323
          
          nu_mech(1,1,3,1) = cmplx(CTensor_rot(1,1,3,1), 0._DP,DPC) ! C'_1131
          nu_mech(1,2,3,1) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu_mech(1,3,3,1) = cmplx(CTensor_rot(3,1,3,1), 0._DP,DPC) ! C'_3131
          nu_mech(2,1,3,1) = cmplx(CTensor_rot(1,1,2,3), 0._DP,DPC) ! C'_1123
          nu_mech(2,2,3,1) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu_mech(2,3,3,1) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu_mech(3,1,3,1) = cmplx(CTensor_rot(1,1,3,3), 0._DP,DPC) ! C'_1133
          nu_mech(3,2,3,1) = cmplx(CTensor_rot(3,3,1,2), 0._DP,DPC) ! C'_3312
          nu_mech(3,3,3,1) = cmplx(CTensor_rot(3,3,3,1), 0._DP,DPC) ! C'_3331
          
          nu_mech(1,1,3,2) = cmplx(CTensor_rot(3,1,1,2), 0._DP,DPC) ! C'_3112
          nu_mech(1,2,3,2) = cmplx(CTensor_rot(2,2,3,1), 0._DP,DPC) ! C'_2231
          nu_mech(1,3,3,2) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu_mech(2,1,3,2) = cmplx(CTensor_rot(2,3,1,2), 0._DP,DPC) ! C'_2312
          nu_mech(2,2,3,2) = cmplx(CTensor_rot(2,2,2,3), 0._DP,DPC) ! C'_2223
          nu_mech(2,3,3,2) = cmplx(CTensor_rot(2,3,2,3), 0._DP,DPC) ! C'_2323
          nu_mech(3,1,3,2) = cmplx(CTensor_rot(3,3,1,2), 0._DP,DPC) ! C'_3312
          nu_mech(3,2,3,2) = cmplx(CTensor_rot(2,2,3,3), 0._DP,DPC) ! C'_2233
          nu_mech(3,3,3,2) = cmplx(CTensor_rot(3,3,2,3), 0._DP,DPC) ! C'_3323
           
          nu_mech(1,1,3,3) = cmplx(CTensor_rot(3,1,3,1), 0._DP,DPC) ! C'_3131
          nu_mech(1,2,3,3) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu_mech(1,3,3,3) = cmplx(CTensor_rot(3,3,3,1), 0._DP,DPC) ! C'_3331
          nu_mech(2,1,3,3) = cmplx(CTensor_rot(2,3,3,1), 0._DP,DPC) ! C'_2331
          nu_mech(2,2,3,3) = cmplx(CTensor_rot(2,3,2,3), 0._DP,DPC) ! C'_2323
          nu_mech(2,3,3,3) = cmplx(CTensor_rot(3,3,2,3), 0._DP,DPC) ! C'_3323
          nu_mech(3,1,3,3) = cmplx(CTensor_rot(3,3,3,1), 0._DP,DPC) ! C'_3331
          nu_mech(3,2,3,3) = cmplx(CTensor_rot(3,3,2,3), 0._DP,DPC) ! C'_3323
          nu_mech(3,3,3,3) = cmplx(CTensor_rot(3,3,3,3), 0._DP,DPC) ! C'_3333
        
        ! Compute stress tensor from auxliary solution (code based on nugradu calculation in fieldsc3D)
        do inat=1,nnat_aux
          do jnat=1,nnat_aux
            stress_temp = matmul(nu_mech(:,:,inat,jnat),gu_aux(:,jnat))
            stress_tensor(1,inat) = stress_tensor(1,inat) + stress_temp(1)
            stress_tensor(2,inat) = stress_tensor(2,inat) + stress_temp(2)
            stress_tensor(3,inat) = stress_tensor(3,inat) + stress_temp(3)
          end do
        end do
        
        ! Symmetry
        stress_tensor(2,1) = stress_tensor(1,2)
        stress_tensor(3,1) = stress_tensor(1,3)
        stress_tensor(3,2) = stress_tensor(2,3)
       
       ! Compute delta tensor ( [[delta]] = [[[[PI]]]] * [[sigma]] )
       ! The delta tensor contains the relative changes in resistivity due to the piezoresistive effect
        call multtensor(PITensor_rot, real(stress_tensor, DP), delta_tensor)
        
        ! Calculate new hole/electron mobilities
        mu_n_x = mu_n / (1._DP + delta_tensor(1,1))
        mu_n_y = mu_n / (1._DP + delta_tensor(2,2))
        mu_n_z = mu_n / (1._DP + delta_tensor(3,3))

        mu_p_x = mu_p / (1._DP + delta_tensor(1,1))
        mu_p_y = mu_p / (1._DP + delta_tensor(2,2))
        mu_p_z = mu_p / (1._DP + delta_tensor(3,3))
        
        Phi_Bn = X_Psi-X_phi_n
        Phi_Bp = X_phi_p-X_Psi

        V_T     = KBoltz*T_ref/Elch
        
        !X_N = n_i * exp((Phi_Bn)/V_T)
        !X_P = n_i * exp((Phi_Bp)/V_T)
        
        
        if ( abs(Phi_Bn).gt.11._DP ) then
          X_N = n_i * (1._DP + (Phi_Bn)/V_T + ((Phi_Bn)/V_T)**2/2._DP)
        else
          X_N = n_i * exp((Phi_Bn)/V_T)
        end if
        
        if ( abs(Phi_Bp).gt.11._DP ) then
          X_P = n_i * (1._DP + (Phi_Bp)/V_T + ((Phi_Bp)/V_T)**2/2._DP)
        else
          X_P = n_i * exp((Phi_Bp)/V_T)
        end if
 
        
        ! Calculate Generation/Recombination rate term
        !GenRate = 0._DP
        if (abs(X_phi_p-X_phi_n) .lt. 10*epsilon(X_phi_n)) then
           Genrate = 0._DP         
        else
           GenRate = (X_N*X_P - n_i**2) / (tau_n*(X_P+n_i)+tau_p*(X_N+n_i))
        end if
         
        ! nu is of the form nu(3,3,nnat,nnat)
        nu(:,:,:,:) = 0._DPC
        ! nu_11
        nu(1,1,1,1) = cmplx(eps0*epsr(1),0._DP,DPC)
        nu(2,2,1,1) = cmplx(eps0*epsr(2),0._DP,DPC)
        nu(3,3,1,1) = cmplx(eps0*epsr(3),0._DP,DPC)
        
        ! nu_22
        nu(1,1,2,2) = cmplx(mu_n_x * X_N,0._DP,DPC)
        nu(2,2,2,2) = cmplx(mu_n_y * X_N,0._DP,DPC)
        nu(3,3,2,2) = cmplx(mu_n_z * X_N,0._DP,DPC)
        
        ! nu_33 ! Continuity Eq multiplied by (-1)
        nu(1,1,3,3) = cmplx(mu_p_x * X_P,0._DP,DPC)
        nu(2,2,3,3) = cmplx(mu_p_y * X_P,0._DP,DPC)
        nu(3,3,3,3) = cmplx(mu_p_z * X_P,0._DP,DPC)

        if (any(real(nu(:,:,:,:)).lt.0)) then
          print*,nu
          print*,'Element', elem
        end if
        
        ! alpha is of the form alpha(nnat,nnat)
        alpha(:,:) = 0._DPC        
        beta(:,:,:) = 0._DPC
        gamma(:,:,:) = 0._DPC
        
        f1(1) = cmplx(Elch*(X_P - X_N + N_D - N_A), 0._DP,DPC)
        f1(2) = cmplx(GenRate, 0._DP,DPC)
        f1(3) = cmplx(-GenRate, 0._DP,DPC) ! Continuity Eq multiplied by (-1)
        
        f2(:,:) = 0._DPC          
      case default
      print*, '***** Unknown Physics Mode!'
      stop

      end select
!
      deallocate(list)
!
      return
      end subroutine pdecoeff3D_scalar
