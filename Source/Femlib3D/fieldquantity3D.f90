      subroutine fieldquantity3D_scalar(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use feminterface,   only: fetchmatparameters
      use feminterface3d, only: field3D, fieldsc3D, lam2xyz, pdecoeff3D
      use femtypes
      use globalvariables3d, only: x, nod, vn, pi, physics, nnat, eps0, Elch, dom, dommat, KBoltz
      use matconstants, only: numparam
      implicit none
      integer (I4B) :: elem
      real (DP) :: lambda(4), phi
      complex (DPC) :: zs
      character (len=*) :: fieldtype, unit, descriptor
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
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
!
!------------------------------------------------------------------------------
!
!  The routine gives the required field quantity with descriptor and unit
!  depending on the physics mode used for a point given in barycentric coordinates.
!  To obtain the time value of the field quantity (for plotting), the real part
!  of zs has got to be taken.
!
!------------------------------------------------------------------------------
!  Input:
!            elem        number of the element
!            fieldtype   field quantity to obtain
!            lambda      barycentric coordinates for point
!            phi         phase angle (in degrees) for harmonic quantities
!  Output:
!            zs          value of field quantity at lambda in elem
!            descriptor  description of field quantity (e.g. electric field strength)
!            unit        unit of field quantity
!            ok          =.false. if fieldtype is not available for the physics mode
!
!  local variables
      real(DP)                :: xyzt(3), tensor(3,3)
      real( DP), allocatable  :: list(:)
      integer (I4B)           :: matindex, nature
      complex (DPC)           :: zval(15,nnat), cfak, z1(3), z2(3), z3(3), z4(3)
      complex (DPC)           :: nu(3,3,nnat,nnat), alpha(nnat,nnat), beta(3,nnat,nnat)
      complex (DPC)           :: gamma(3,nnat,nnat), f1(nnat), f2(3,nnat)
      real (DP)               :: T_ref, n_i, mu_n, mu_p, V_T, X_n, X_p, L_0, N_A,N_D
      logical                 :: typ(5)
!------------------------------------------------------------------------------
!
!  initializations
      ok=.true.
      descriptor = ''
      unit = ''
!
!  get world coordinates
      call lam2xyz(lambda,elem,xyzt,nod,vn)
!  get pde coefficients
      if(associated(x)) call pdecoeff3D(elem,xyzt,nu,alpha,beta,f1,f2,gamma)
      
      select case (physics)
!
!------------------------------------------------------------------------------
!
!__
! SEMICONDUCTOR:
      case ('3DPOISSON')
          select case (fieldtype)
            case ('PHI')
              descriptor = 'Potential'
              unit = 'V'
              typ=(/.true., .false., .false., .false., .false./)
              call fieldsc3D(elem,lambda,typ,zval)
              cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
              nature=1
              zs = real(zval(1,nature)*cfak)
          
            
            case default
            ok=.false.
        
         end select
          
      case ('3DSEMICONDUCTOR')
           !  determine material index for given element in region
           matindex=dommat(dom(elem))
           !  read parameters from internal list and assign them to local variables
           allocate (list(numparam))
           call fetchmatparameters(list, matindex, xyzt)
           T_ref   = list(6)
           mu_n    = list(7)
           mu_p    = list(8)
           n_i     = list(9)        
           V_T     = KBoltz*T_ref/Elch
        
         select case (fieldtype)
          case ('PSI')
            descriptor = 'Potential'
            unit = 'V'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = real(zval(1,nature)*cfak)
            
          case ('PHI_N')
            descriptor = 'Electron Quasi-Fermi-Level'
            unit = 'V'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=2
            zs = real(zval(1,nature)*cfak)
            
          case ('PHI_P')
            descriptor = 'Hole Quasi-Fermi-Level'
            unit = 'V'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=3
            zs = real(zval(1,nature)*cfak)  
            
          case ('N')
            descriptor = 'Electron Density'
            unit = '1/m3'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            X_N = n_i*exp((zval(1,1)-zval(1,2))/V_T)

            zs = real(X_N*cfak)     
            
          case ('P')
            descriptor = 'Hole Density'
            unit = '1/m3'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
      
            X_P = n_i*exp((-zval(1,1)+zval(1,3))/V_T)
            
            zs = real(X_P*cfak)
            
            case ('N_A')
            descriptor = 'Dopant concentration'
            unit = '1/m3'     
            
            zs = real(list(4))
            
            case ('N_D')
            descriptor = 'Dopant concentration'
            unit = '1/m3'     
            
            zs = real(list(5))
            
          case ('N_D-N_A')
            descriptor = 'Dopant concentration'
            unit = '1/m3'     
            
            zs = real(list(5)-list(4))
            
          case ('E_FIELD')
            descriptor = 'E-Field Magnitude'
            unit = 'V/m'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            nature = 1
            z1(1) = real(-zval(4,nature)*cfak) ! x-component
            z1(2) = real(-zval(5,nature)*cfak) ! y-component
            z1(3) = real(-zval(6,nature)*cfak) ! z-component

            zs = sqrt(z1(1)**2+z1(2)**2+z1(3)**2) 
            
          case ('J_N_DRIFT')
            descriptor = 'Electron Drift Current Density Magnitude'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            X_N = n_i*exp((zval(1,1)-zval(1,2))/V_T)

            ! J_n_drift = q*n*mu_n*(-grad(Psi))
            z1(1) = real(Elch*X_N*mu_n*-zval(4,1)*cfak) ! x-component
            z1(2) = real(Elch*X_N*mu_n*-zval(5,1)*cfak) ! y-component
            z1(3) = real(Elch*X_N*mu_n*-zval(6,1)*cfak) ! z-component        

            zs = sqrt(z1(1)**2+z1(2)**2+z1(3)**2) 
          
          case ('J_P_DRIFT')
            descriptor = 'Hole Drift Current Density Magnitude'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            X_P = n_i*exp((-zval(1,1)+zval(1,3))/V_T)
            
            ! J_p_drift = q*p*mu_p*(-grad(Psi))
            z1(1) = real(Elch*X_P*mu_p*-zval(4,1)*cfak) ! x-component
            z1(2) = real(Elch*X_P*mu_p*-zval(5,1)*cfak) ! y-component
            z1(3) = real(Elch*X_P*mu_p*-zval(6,1)*cfak) ! z-component      

            zs = sqrt(z1(1)**2+z1(2)**2+z1(3)**2) 
          
          case ('J_N_DIFF')
            descriptor = 'Electron Diffusion Current Density Magnitude'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            X_N = n_i*exp((zval(1,1)-zval(1,2))/V_T)
            
            ! J_n_diff = q*mu_n*n*grad(Psi-phi_n)
            z1(1) = real(Elch*mu_n*X_N*(zval(4,1)-zval(4,2))*cfak) ! x-component
            z1(2) = real(Elch*mu_n*X_N*(zval(5,1)-zval(5,2))*cfak) ! y-component
            z1(3) = real(Elch*mu_n*X_N*(zval(6,1)-zval(6,2))*cfak) ! z-component 

            zs = sqrt(z1(1)**2+z1(2)**2+z1(3)**2) 
            
          case ('J_P_DIFF')
            descriptor = 'Hole Diffusion Current Density Magnitude'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            X_P = n_i*exp((-zval(1,1)+zval(1,3))/V_T)
            
            ! J_p_diff = q*mu_p*p*grad(Psi-phi_p)
            z1(1) = real(Elch*mu_p*X_P*(zval(4,1)-zval(4,3))*cfak) ! x-component
            z1(2) = real(Elch*mu_p*X_P*(zval(5,1)-zval(5,3))*cfak) ! y-component
            z1(3) = real(Elch*mu_p*X_P*(zval(6,1)-zval(6,3))*cfak) ! z-component      

            zs = sqrt(z1(1)**2+z1(2)**2+z1(3)**2) 
            
          case ('J_N')
            descriptor = 'Electron Current Density Magnitude'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            X_N = n_i*exp((zval(1,1)-zval(1,2))/V_T)
            
            ! J_n = -q*mu_n*n*grad(phi_n)
            z1(1) = real(-Elch*mu_n*X_N*zval(4,2)*cfak) ! x-component
            z1(2) = real(-Elch*mu_n*X_N*zval(5,2)*cfak) ! y-component
            z1(3) = real(-Elch*mu_n*X_N*zval(6,2)*cfak) ! z-component 
            
            zs = sqrt(z1(1)**2+z1(2)**2+z1(3)**2) 
            
          case ('J_P')
            descriptor = 'Hole Current Density Magnitude'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            X_P = n_i*exp((-zval(1,1)+zval(1,3))/V_T)
            
            ! J_p = -q*mu_p*p*grad(phi_p)
            z1(1) = real(-Elch*mu_p*X_P*zval(4,3)*cfak) ! x-component
            z1(2) = real(-Elch*mu_p*X_P*zval(5,3)*cfak) ! y-component
            z1(3) = real(-Elch*mu_p*X_P*zval(6,3)*cfak) ! z-component 
            
            zs = sqrt(z1(1)**2+z1(2)**2+z1(3)**2) 
            
          case ('J')
            descriptor = 'Total Current Density Magnitude'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            X_N = n_i*exp((zval(1,1)-zval(1,2))/V_T)
            X_P = n_i*exp((-zval(1,1)+zval(1,3))/V_T)
            
            ! Electron Total Current
            ! J_n = -q*mu_n*n*grad(phi_n)
            z1(1) = real(-Elch*mu_n*X_N*zval(4,2)*cfak) ! x-component
            z1(2) = real(-Elch*mu_n*X_N*zval(5,2)*cfak) ! y-component
            z1(3) = real(-Elch*mu_n*X_N*zval(6,2)*cfak) ! z-component  
          
            ! Hole Total Current
            ! J_p = -q*mu_p*p*grad(phi_p)
            z2(1) = real(-Elch*mu_p*X_P*zval(4,3)*cfak) ! x-component
            z2(2) = real(-Elch*mu_p*X_P*zval(5,3)*cfak) ! y-component
            z2(3) = real(-Elch*mu_p*X_P*zval(6,3)*cfak) ! z-component  

            ! Total Current
            z3 = z1 + z2
            
            zs = sqrt(z3(1)**2+z3(2)**2+z3(3)**2) 
            
            
          case default
          ok=.false.
        
        end select
        deallocate(list)
!__
! STRAIN / STRESS      
      case ('3DSOLIDMECHANICS')
        nature=3
        select case (fieldtype)
        case ('DX')
          descriptor = 'X-Displacement'
          unit = 'm'
          typ=(/.true., .false., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=1
          zs = real(zval(1,nature)*cfak)
        case ('DY')
          descriptor = 'Y-Displacement'
          unit = 'm'
          typ=(/.true., .false., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=2
          zs = real(zval(1,nature)*cfak)
        case ('DZ')
          descriptor = 'Z-Displacement'
          unit = 'm'
          typ=(/.true., .false., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          nature=3
          zs = real(zval(1,nature)*cfak)
        case ('VAN_MISES_STRESS','VON_MISES_STRESS')
          descriptor = 'von Mises stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stress tensor components
          tensor(1,1) = real(zval( 8,1)  *cfak)
          tensor(2,2) = real(zval( 9,2)  *cfak)
          tensor(3,3) = real(zval(10,3)  *cfak)
          tensor(2,3) = 0.5_DP*real((zval(9,3)+zval(10,2))  *cfak)
          tensor(1,3) = 0.5_DP*real((zval(8,3)+zval(10,1))  *cfak)
          tensor(1,2) = 0.5_DP*real((zval(8,2)+zval( 9,1))  *cfak)
          zs = sqrt( ( (tensor(1,1)-tensor(2,2))**2 +                   &
         &             (tensor(2,2)-tensor(3,3))**2 +                   &
         &             (tensor(3,3)-tensor(1,1))**2 +                   &
         &         6*(  tensor(1,2)**2 +                                &
         &              tensor(2,3)**2 +                                &
         &              tensor(1,3)**2  )  ) /2._DP )
        case ('SXX')
          descriptor = 'Sxx stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stress tensor components
          zs = real(zval( 8,1)  *cfak)
        case ('SYY')
          descriptor = 'Syy stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stress tensor components
          zs = real(zval( 9,2)  *cfak)
        case ('SZZ')
          descriptor = 'Szz stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stress tensor components
          zs = real(zval(10,3)  *cfak)
        case ('SYZ','SZY')
          descriptor = 'Syz stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stress tensor components
          zs = 0.5_DP*real((zval(9,3)+zval(10,2))  *cfak)
        case ('SXZ','SZX')
          descriptor = 'Sxz stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stress tensor components
          zs = 0.5_DP*real((zval(8,3)+zval(10,1))  *cfak)
        case ('SXY','SYX')
          descriptor = 'Sxy stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stress tensor components
          zs = 0.5_DP*real((zval(8,2)+zval( 9,1))  *cfak)
        case ('EXX')
          descriptor = 'Exx strain'
          unit = ''
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stain tensor components
          zs = real(zval( 4,1)  *cfak)
        case ('EYY')
          descriptor = 'Eyy strain'
          unit = ''
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stain tensor components
          zs = real(zval( 5,2)  *cfak)
        case ('EZZ')
          descriptor = 'Ezz strain'
          unit = ''
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stain tensor components
          zs = real(zval( 6,3)  *cfak)
        case ('EYZ','EZY')
          descriptor = 'Eyz strain'
          unit = ''
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stain tensor components
          zs = 0.5_DP*real((zval( 5,3)+zval( 6,2))  *cfak)
        case ('EXZ','EZX')
          descriptor = 'Exz strain'
          unit = ''
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stain tensor components
          zs = 0.5_DP*real((zval( 4,3)+zval( 6,1))  *cfak)
        case ('EXY','EYX')
          descriptor = 'Exy strain'
          unit = ''
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stain tensor components
          zs = 0.5_DP*real((zval( 4,2)+zval( 5,1))  *cfak)
!        case ('TRESCA STRESS')
!          descriptor = 'Tresca stress'
!          unit = 'N/m2'
!          typ = (/.false.,.true.,.false.,.false.,.false./)
!          call field(elem,lambda,typ,zval)
!          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!          call principlevalues(real(zval(7,1)*cfak), real(zval(8,2)*cfak),  &
!     &             0.5_DP*real((zval(7,2)+zval(8,1))*cfak),sigma1,sigma2)
!          zs = max(abs(sigma1),abs(sigma2),abs(sigma1-sigma2))
!          typ = .false.

        case default
          ok=.false.
        end select
        
        case ('3DELECTROSTATICS')
        nature=1
        select case (fieldtype)
          case ('U')
            descriptor = 'Voltage'
            unit = 'V'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = zval(1,nature)*cfak
          case ('Ex')
            descriptor = 'Electric Field Strength'
            unit = 'V'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(4,nature)*cfak
          case ('Ey')
            descriptor = 'Electric Field Strength'
            unit = 'V'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(5,nature)*cfak
          case ('Ez')
            descriptor = 'Electric Field Strength'
            unit = 'V'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(6,nature)*cfak
          case ('E')
            descriptor = 'electric field strength'
            unit = 'V/m'
            typ = (/.false.,.true.,.false.,.false.,.false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = sqrt(zval(4,nature)**2 + zval(5,nature)**2 + zval(6,nature)**2)
          case ('Dx')
            descriptor = 'Electric Flux Density'
            unit = 'AS/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(8,nature)*cfak*eps0
          case ('Dy')
            descriptor = 'Electric Flux Density'
            unit = 'AS/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(9,nature)*cfak*eps0
          case ('Dz')
            descriptor = 'Electric Flux Density'
            unit = 'AS/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(10,nature)*cfak*eps0
          case ('D')
            descriptor = 'Electric Flux Density'
            unit = 'AS/m2'
            typ = (/.false.,.true.,.false.,.false.,.false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = sqrt(zval(8,nature)**2 + zval(9,nature)**2 + zval(10,nature)**2)           
          case default
          ok=.false.
          end select
          
          case ('3DHEATTR')
        nature=1
        select case (fieldtype)
          case ('T')
            descriptor = 'Temperature'
            unit = 'K'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = zval(1,nature)*cfak
          
          case ('Qx')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(8,nature)*cfak
          case ('Qy')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(9,nature)*cfak
          case ('Qz')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(10,nature)*cfak
          case ('Q')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ = (/.false.,.true.,.false.,.false.,.false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = sqrt(zval(8,nature)**2 + zval(9,nature)**2 + zval(10,nature)**2)           
          case default
          ok=.false.
          end select

case ('3DTHERMOELECTRIC')
        select case (fieldtype)
        case ('T')
            nature=1
            descriptor = 'Temperature'
            unit = 'K'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            zs = zval(1,nature)*cfak
          
          case ('QX')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(8,nature)*cfak
          case ('QY')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(9,nature)*cfak
          case ('QZ')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(10,nature)*cfak
          case ('Q')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ = (/.false.,.true.,.false.,.false.,.false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = sqrt(zval(8,nature)**2 + zval(9,nature)**2 + zval(10,nature)**2)           
          
          case ('POTENTIAL')
            descriptor = 'Voltage'
            unit = 'V'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=2
            zs = zval(1,nature)*cfak
          
          case ('Jx')
            descriptor = 'CURRENT DENSITY IN X'
            unit = 'A/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=2
            zs = -zval(8,nature)*cfak

          case ('Jy')
            descriptor = 'CURRENT DENSITY IN Y'
            unit = 'A/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=2
            zs = -zval(9,nature)*cfak

          case ('Jz')
            descriptor = 'CURRENT DENSITY IN Z'
            unit = 'A/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=2
            zs = -zval(10,nature)*cfak
          case default
          ok=.false.
          end select
          
        case ('3DNONLINHEATTR')
        nature=1
        select case (fieldtype)
          case ('T')
            descriptor = 'Temperature'
            unit = 'K'
            typ=(/.true., .false., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = zval(1,nature)*cfak
          
          case ('QX')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(8,nature)*cfak
          case ('QY')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(9,nature)*cfak
          case ('QZ')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = -zval(10,nature)*cfak
          case ('Q')
            descriptor = 'HEAT FLUX'
            unit = 'W/m2'
            typ = (/.false.,.true.,.false.,.false.,.false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = sqrt(zval(8,nature)**2 + zval(9,nature)**2 + zval(10,nature)**2)           
          case default
          ok=.false.
          end select
!__
! Error:
      case default
        print *,"***** unknown physics mode:",physics
        stop
!__
!
      end select
      end subroutine fieldquantity3D_scalar



      subroutine fieldquantity3D_vector(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use feminterface3d, only: field3D, fieldsc3D, lam2xyz, pdecoeff3D
      use feminterface,   only: fetchmatparameters
      use femtypes
      use globalvariables3d, only: nod, vn, pi, physics, nnat, Elch, KBoltz, dom, dommat, eps0
      use matconstants, only: numparam
      implicit none
      integer (I4B) :: elem
      real (DP) :: lambda(4), phi
      complex (DPC) :: zs(:)
      character (len=*) :: fieldtype, unit, descriptor
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
!------------------------------------------------------------------------------
!
!  The routine gives the required field quantity with descriptor and unit
!  depending on the physics mode used for a point given in barycentric coordinates.
!  To obtain the time value of the field quantity (for plotting), the real part
!  of zs has got to be taken.
!
!------------------------------------------------------------------------------
!  Input:
!            elem        number of the element
!            fieldtype   field quantity to obtain
!            lambda      barycentric coordinates for point
!            phi         phase angle (in degrees) for harmonic quantities
!  Output:
!            zs          value of field quantity at lambda in elem
!            descriptor  description of field quantity (e.g. electric field strength)
!            unit        unit of field quantity
!            ok          =.false. if fieldtype is not available for the physics mode
!
!  local variables
      real(DP)                :: xyzt(3)
      real( DP), allocatable  :: list(:)
      integer (I4B)           :: matindex, nature
      complex (DPC)           :: zval(15,nnat), cfak, z1(3), z2(3), z3(3), z4(3)
      complex (DPC)           :: nu(3,3,nnat,nnat), alpha(nnat,nnat), beta(3,nnat,nnat)
      complex (DPC)           :: gamma(3,nnat,nnat), f1(nnat), f2(3,nnat)
      real (DP)               :: T_ref, n_i, mu_n, mu_p, V_T, X_n, X_p, L_0, N_A,N_D
      logical                 :: typ(5)
!------------------------------------------------------------------------------
!
!  initializations
      ok=.true.
      descriptor = ''
      unit = ''
!
!  get world coordinates
      call lam2xyz(lambda,elem,xyzt,nod,vn)
!  get pde coefficients
      call pdecoeff3D(elem,xyzt,nu,alpha,beta,f1,f2,gamma)

      select case (physics)
!
!------------------------------------------------------------------------------
!
          
      case ('3DPOISSON')
        nature=1
				select case (fieldtype)
        case ('E')
          descriptor = 'Electric Field Stength'
          unit = 'V/m'
          typ=(/.false., .true., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp( cmplx( 0._DP, phi*pi/180._DP, DPC ) )
          zs(1) = -zval(4,nature)*cfak
          zs(2) = -zval(5,nature)*cfak
          zs(3) = -zval(6,nature)*cfak

          case default
          ok=.false.
        
        end select
      
      case ('3DSEMICONDUCTOR')
          !  determine material index for given element in region
          matindex=dommat(dom(elem))
          !  read parameters from internal list and assign them to local variables
          allocate (list(numparam))
          call fetchmatparameters(list, matindex, xyzt)
          T_ref   = list(6)
          mu_n    = list(7)
          mu_p    = list(8)
          n_i     = list(9)
          
          V_T     = KBoltz*T_ref/Elch
        
         select case (fieldtype)           
          case ('E_FIELD')
            descriptor = 'E-Field'
            unit = 'V/m'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            nature = 1
            zs(1) = real(-zval(4,nature)*cfak) ! x-component
            zs(2) = real(-zval(5,nature)*cfak) ! y-component
            zs(3) = real(-zval(6,nature)*cfak) ! z-component

            
          case ('J_N_DRIFT')
            descriptor = 'Electron Drift Current Density'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))

            zs(1) = real(Elch*mu_n*zval(1,2)*-zval(4,1)*cfak) ! x-component
            zs(2) = real(Elch*mu_n*zval(1,2)*-zval(5,1)*cfak) ! y-component
            zs(3) = real(Elch*mu_n*zval(1,2)*-zval(6,1)*cfak) ! z-component        

          case ('J_P_DRIFT')
            descriptor = 'Hole Drift Current Density'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))

            zs(1) = real(Elch*mu_p*zval(1,3)*-zval(4,1)*cfak) ! x-component
            zs(2) = real(Elch*mu_p*zval(1,3)*-zval(5,1)*cfak) ! y-component
            zs(3) = real(Elch*mu_p*zval(1,3)*-zval(6,1)*cfak) ! z-component      
          
          case ('J_N_DIFF')
            descriptor = 'Electron Diffusion Current Density'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))

            zs(1) = real(Elch*mu_n*V_T*zval(4,2)*cfak) ! x-component
            zs(2) = real(Elch*mu_n*V_T*zval(5,2)*cfak) ! y-component
            zs(3) = real(Elch*mu_n*V_T*zval(6,2)*cfak) ! z-component 

          case ('J_P_DIFF')
            descriptor = 'Hole Diffusion Current Density'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))

            zs(1) = real(-Elch*mu_p*V_T*zval(4,3)*cfak) ! x-component
            zs(2) = real(-Elch*mu_p*V_T*zval(5,3)*cfak) ! y-component
            zs(3) = real(-Elch*mu_p*V_T*zval(6,3)*cfak) ! z-component      

          case ('J_N')
            descriptor = 'Electron Current Density'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            ! Electron Drift Current
            z1(1) = real(Elch*mu_n*zval(1,2)*-zval(4,1)*cfak) ! x-component
            z1(2) = real(Elch*mu_n*zval(1,2)*-zval(5,1)*cfak) ! y-component
            z1(3) = real(Elch*mu_n*zval(1,2)*-zval(6,1)*cfak) ! z-component  
          
            ! Electron Diffusion Current
            z2(1) = real(Elch*mu_n*V_T*zval(4,2)*cfak) ! x-component
            z2(2) = real(Elch*mu_n*V_T*zval(5,2)*cfak) ! y-component
            z2(3) = real(Elch*mu_n*V_T*zval(6,2)*cfak) ! z-component 
            
            ! Electron Total Current
            zs = z1 + z2
            
          case ('J_P')
            descriptor = 'Hole Current Density'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            ! Hole Drift Current
            z1(1) = real(Elch*mu_p*zval(1,3)*-zval(4,1)*cfak) ! x-component
            z1(2) = real(Elch*mu_p*zval(1,3)*-zval(5,1)*cfak) ! y-component
            z1(3) = real(Elch*mu_p*zval(1,3)*-zval(6,1)*cfak) ! z-component   
          
            ! Hole Diffusion Current
            z2(1) = real(-Elch*mu_p*V_T*zval(4,3)*cfak) ! x-component
            z2(2) = real(-Elch*mu_p*V_T*zval(5,3)*cfak) ! y-component
            z2(3) = real(-Elch*mu_p*V_T*zval(6,3)*cfak) ! z-component      
            
            ! Hole Total Current
            zs = z1 + z2
            
          case ('J')
            descriptor = 'Total Current Density'
            unit = 'A/m'
            typ=(/.true., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            
            ! Electron Drift Current
            z1(1) = real(Elch*mu_n*zval(1,2)*-zval(4,1)*cfak) ! x-component
            z1(2) = real(Elch*mu_n*zval(1,2)*-zval(5,1)*cfak) ! y-component
            z1(3) = real(Elch*mu_n*zval(1,2)*-zval(6,1)*cfak) ! z-component  
          
            ! Electron Diffusion Current
            z2(1) = real(Elch*mu_n*V_T*zval(4,2)*cfak) ! x-component
            z2(2) = real(Elch*mu_n*V_T*zval(5,2)*cfak) ! y-component
            z2(3) = real(Elch*mu_n*V_T*zval(6,2)*cfak) ! z-component 
            
            ! Electron Total Current
            z3 = z1 + z2
            
            ! Hole Drift Current
            z1(1) = real(Elch*mu_p*zval(1,3)*-zval(4,1)*cfak) ! x-component
            z1(2) = real(Elch*mu_p*zval(1,3)*-zval(5,1)*cfak) ! y-component
            z1(3) = real(Elch*mu_p*zval(1,3)*-zval(6,1)*cfak) ! z-component   
          
            ! Hole Diffusion Current
            z2(1) = real(-Elch*mu_p*V_T*zval(4,3)*cfak) ! x-component
            z2(2) = real(-Elch*mu_p*V_T*zval(5,3)*cfak) ! y-component
            z2(3) = real(-Elch*mu_p*V_T*zval(6,3)*cfak) ! z-component     
            
            ! Hole Total Current
            z4 = z1 + z2
            
            ! Total Current
            zs = z3 + z4
            
          case default
          ok=.false.
        
          end select
!__
! STRAIN / STRESS
      case ('3DSOLIDMECHANICS')
        nature=3
        select case (fieldtype)
        case ('DEL')
          descriptor = 'Displacement'
          unit = 'm'
          typ=(/.true., .false., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  nature 1   => x
          zs(1) = zval(1,1)*cfak
!  nature 2   => y
          zs(2) = zval(1,2)*cfak
!  nature 3   => z
          zs(3) = zval(1,3)*cfak
          case default
          ok=.false.
        end select
      
       case ('3DELECTROSTATICS')
        nature=1
        select case (fieldtype)
        case ('E')
          descriptor = 'Electric Field Stength'
          unit = 'V/m'
          typ=(/.false., .true., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = -zval(4,nature)*cfak
          zs(2) = -zval(5,nature)*cfak
          zs(3) = -zval(6,nature)*cfak
        case ('D')
          descriptor = 'Electric Flux Density'
          unit = 'AS/m2'
          typ=(/.false., .true., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = -zval(8,nature)*cfak*eps0
          zs(2) = -zval(9,nature)*cfak*eps0
          zs(3) = -zval(10,nature)*cfak*eps0
          case default
          ok=.false.
        end select
        
       case ('3DTHERMOELECTRIC')
!parameter( 1) =     angle      ,   0.        ,                              !  angle between the  x-axis and the first principle axis of the material
!parameter( 2) =     angle      ,   0.        ,                              !  angle between the  y-axis and the first principle axis of the material
!
!parameter( 3) =     kappa11    ,   0.02      ,       W/m.K                  !  thermal conductivity in the first principle axis
!parameter( 4) =     kappa22    ,   0.02      ,       W/m.K                  !  thermal conductivity in the second principle axis
!parameter( 5) =     kappa33    ,   0.02      ,       W/m.K                  !  thermal conductivity in the third principle axis
!
!parameter( 6) =     sigma11    ,   1.        ,       S/m                    !  electrical conductivity in the first principle axis
!parameter( 7) =     sigma22    ,   1.        ,       S/m                    !  electrical conductivity in the second principle axis
!parameter( 8) =     sigma33    ,   1.        ,       S/m                    !  electrical conductivity in the third principle axis
!
!parameter( 9) =      seeb11    ,   0.        ,       V/K                    !  thermoelectric coupling Seebeck coefficient in the first principle axis
!parameter(10) =      seeb22    ,   0.        ,       V/K                    !  thermoelectric coupling Seebeck coefficient in the second principle axis
!parameter(11) =      seeb33    ,   0.        ,       V/K                    !  thermoelectric coupling Seebeck coefficient in the third principle axis
         matindex = dommat(dom(elem))
         allocate(list(numparam))
         call fetchmatparameters(list,matindex,xyzt,elem)
         select case (fieldtype)
         case ('QTH')
          descriptor = 'Heatflux Density Thermally Induced '
          unit = 'W/m2'
          typ=(/.false., .true., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = -list(3)*zval(4,1)*cfak
          zs(2) = -list(4)*zval(5,1)*cfak
          zs(3) = -list(5)*zval(6,1)*cfak
        case ('QEL')
          descriptor = 'Heatflux Density Electrically Induced'
          unit = 'W/m2'
          typ=(/.true., .true., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = ((list(9)*zval(1,1)*list(6))*((-zval(4,2))-(list(9)*zval(4,1))))*cfak
          zs(2) = ((list(10)*zval(1,1)*list(7))*((-zval(5,2))-(list(10)*zval(5,1))))*cfak
          zs(3) = ((list(11)*zval(1,1)*list(8))*((-zval(6,2))-(list(11)*zval(6,1))))*cfak

        case ('QTOT')
          descriptor = 'Heatflux Density Total value'
          unit = 'W/m2'
          typ=(/.true., .true., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = (-list(3)*zval(4,1))*cfak+((list(9)*zval(1,1)*list(6))*(-zval(4,2)-(list(9)*zval(4,1))))*cfak
          zs(2) = (-list(4)*zval(5,1))*cfak+((list(10)*zval(1,1)*list(7))*(-zval(5,2)-(list(10)*zval(5,1))))*cfak
          zs(3) = (-list(5)*zval(6,1))*cfak+((list(11)*zval(1,1)*list(8))*(-zval(6,2)-(list(11)*zval(6,1))))*cfak
        case ('Q')
            descriptor = 'HEAT FLUX MAGNITUDE'
            unit = 'W/m2'
            typ = (/.false.,.true.,.false.,.false.,.false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            nature=1
            zs = sqrt(zval(8,nature)**2 + zval(9,nature)**2 + zval(10,nature)**2)   
        case ('E')
            nature=2
            descriptor = 'Electric Field'
            unit = 'V/m'
            typ=(/.false., .true., .false., .false., .false./)
            call fieldsc3D(elem,lambda,typ,zval)
            cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
            zs(1) = zval(4,nature)*cfak
            zs(2) = zval(5,nature)*cfak
            zs(3) = zval(6,nature)*cfak
           

          case default
          ok=.false.
        end select
        
        case ('3DHEATTR')
        nature=1
        select case (fieldtype)
        case ('Q')
          descriptor = 'Heat Flux'
          unit = 'W/m2'
          typ=(/.false., .true., .false., .false., .false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
          zs(1) = -zval(8,nature)*cfak
          zs(2) = -zval(9,nature)*cfak
          zs(3) = -zval(10,nature)*cfak
          case default
          ok=.false.
        end select
!__
! Error:
      case default
        print *,"***** unknown physics mode:",physics
        stop
!__
!
      end select
      end subroutine fieldquantity3D_vector



      subroutine fieldquantity3D_tensor(elem,fieldtype,lambda,phi,zs,descriptor,unit,ok)
      use feminterface3d, only: field3D, fieldsc3D, lam2xyz, pdecoeff3D
      use femtypes
      use globalvariables3d, only: nod, vn, pi, physics, nnat
      use matconstants, only: numparam
      implicit none
      integer (I4B) :: elem
      real (DP) :: lambda(4), phi
      complex (DPC) :: zs(:,:)
      character (len=*) :: fieldtype, unit, descriptor
      logical :: ok
      intent (in) :: elem, fieldtype, lambda, phi
      intent (out) :: descriptor, unit, zs, ok
!------------------------------------------------------------------------------
!
!  The routine gives the required field quantity with descriptor and unit
!  depending on the physics mode used for a point given in barycentric coordinates.
!  To obtain the time value of the field quantity (for plotting), the real part
!  of zs has got to be taken.
!
!------------------------------------------------------------------------------
!  Input:
!            elem        number of the element
!            fieldtype   field quantity to obtain
!            lambda      barycentric coordinates for point
!            phi         phase angle (in degrees) for harmonic quantities
!  Output:
!            zs          value of field quantity at lambda in elem
!            descriptor  description of field quantity (e.g. electric field strength)
!            unit        unit of field quantity
!            ok          =.false. if fieldtype is not available for the physics mode
!
!  local variables
      real(DP)              :: xyzt(3)
      integer (I4B)         :: nature
      complex (DPC)         :: zval(15,nnat), cfak
      complex (DPC)         :: nu(3,3,nnat,nnat), alpha(nnat,nnat), beta(3,nnat,nnat)
      complex (DPC)         :: gamma(3,nnat,nnat), f1(nnat), f2(3,nnat)
      logical               :: typ(5)
!------------------------------------------------------------------------------
!
!  initializations
      ok=.true.
      descriptor = ''
      unit = ''
!
!  get world coordinates
      call lam2xyz(lambda,elem,xyzt,nod,vn)
!  get pde coefficients
      call pdecoeff3D(elem,xyzt,nu,alpha,beta,f1,f2,gamma)
      select case (physics)
!
!------------------------------------------------------------------------------
!
!__
! STRAIN / STRESS
      case ('3DSOLIDMECHANICS')
        nature=3
        select case (fieldtype)
        case ('STRESS')
          descriptor = 'Stress'
          unit = 'N/m2'
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stress tensor components
!  Sxx
          zs(1,1) = real(zval( 8,1)  *cfak)
!  Syy
          zs(2,2) = real(zval( 9,2)  *cfak)
!  Szz
          zs(3,3) = real(zval(10,3)  *cfak)
!  Syz, Szy
          zs(2,3) = 0.5_DP*real((zval(9,3)+zval(10,2))  *cfak)
          zs(3,2) = zs(2,3)
!  Sxz, Szx
          zs(1,3) = 0.5_DP*real((zval(8,3)+zval(10,1))  *cfak)
          zs(3,1) = zs(1,3)
!  Sxy, Syx
          zs(1,2) = 0.5_DP*real((zval(8,2)+zval( 9,1))  *cfak)
          zs(2,1) = zs(1,2)
        case ('STRAIN')
          descriptor = 'Strain'
          unit = ''
          typ = (/.false.,.true.,.false.,.false.,.false./)
          call fieldsc3D(elem,lambda,typ,zval)
          cfak = exp(cmplx(0._DP,phi*pi/180._DP,DPC))
!  stress tensor components
!  Exx
          zs(1,1) = real(zval( 4,1)  *cfak)
!  Eyy
          zs(2,2) = real(zval( 5,2)  *cfak)
!  Ezz
          zs(3,3) = real(zval( 6,3)  *cfak)
!  Eyz, Ezy
          zs(2,3) = 0.5_DP*real((zval(5,3)+zval( 6,2))  *cfak)
          zs(3,2) = zs(2,3)
!  Exz, Ezx
          zs(1,3) = 0.5_DP*real((zval(4,3)+zval( 6,1))  *cfak)
          zs(3,1) = zs(1,3)
!  Exy, Eyx
          zs(1,2) = 0.5_DP*real((zval(4,2)+zval( 5,1))  *cfak)
          zs(2,1) = zs(1,2)
        case default
          ok=.false.
        end select
!__
! Error:
      case default
        print *,"***** unknown physics mode:",physics
        stop
!__
!
      end select
      end subroutine fieldquantity3D_tensor
