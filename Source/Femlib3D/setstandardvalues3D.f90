      subroutine setstandardvalues3D
      use feminterface, only: low2hi
      use femtypes
      use matconstants
      use globalvariables3D, only: physics
      implicit none
!
!    $Revision: 1.7 $
!    $Date: 2015/11/11 17:30:00 $
!    $Author: m_kasper $
!
!  set the default values for the material parameters 
!
      integer(I4B) i, j
      real(DP), allocatable :: svalue(:)
!
      select case(physics)
!
        case ('TEWAVE','TMWAVE')
        numparam = 19
!
        allocate(parameternames(numparam),svalue(numparam))
!  type of material
        parameternames(1)='pml'
        svalue(1)=0._DP
!  loss tangent for permeability
        parameternames(2)='tan_delta_mu'
        svalue(2)=0._DP
!  loss tangent for permittivity
        parameternames(3)='tan_delta_eps'
        svalue(3)=0._DP
!------------------------------------------------------------
! relative mu and epsilon for 3D
!------------------------------------------------------------
!  relative permeability in the first principle axis
        parameternames(4)='mur11'
        svalue(4)=1._DP
!  relative permeability in the second principle axis
        parameternames(5)='mur22'
        svalue(5)=1._DP
!  relative permeability in the third principle axis
        parameternames(6)='mur33'
        svalue(6)=1._DP
!  relative permittivity in the first principle axis
        parameternames(7)='epsr11'
        svalue(7)=1._DP
!  relative permittivity in the second principle axis
        parameternames(8)='epsr22'
        svalue(8)=1._DP
!  relative permittivity in the third principle axis
        parameternames(9)='epsr33'
        svalue(9)=1._DP
!------------------------------------------------------------
!  current density in the first principle axis
        parameternames(10)='J11'
        svalue(10)=0._DP
!  current density in the second principle axis
        parameternames(11)='J22'
        svalue(11)=0._DP
!  current density in the third principle axis
        parameternames(12)='J22'
        svalue(12)=0._DP
!  phase angle of the current density
        parameternames(13)='phase'
        svalue(13)=0._DP
!------------------------------------------------------------
!  direction of first principle axis 
        parameternames(14)='p1_v1'
        svalue(14)=1._DP
        parameternames(15)='p1_v2'
        svalue(15)=0._DP
        parameternames(16)='p1_v3'
        svalue(16)=0._DP
!  direction of second principle axis 
        parameternames(17)='p2_v1'
        svalue(17)=0._DP
        parameternames(18)='p2_v2'
        svalue(18)=1._DP
        parameternames(19)='p2_v3'
        svalue(19)=0._DP
!------------------------------------------------------------
      case ('3DHEATTR')
        numparam=14
!
        allocate(parameternames(numparam),svalue(numparam))
!  values for air!!!
!  thermal conductivity in the first principle axis
        parameternames(1)='lam11'
        svalue(1)=0.02_DP
!  thermal conductivity in the second principle axis
        parameternames(2)='lam22'
        svalue(2)=0.02_DP
!  thermal conductivity in the third principle axis
        parameternames(3)='lam33'
        svalue(3)=0.02_DP
!  angle between the x-axis and the first principle axis
        parameternames(4)='angle'
        svalue(4)=0._DP
!  convection coefficient on the top surface
        parameternames(5)='h_o'
        svalue(5)=0._DP
!  convection coefficient on the bottom surface
        parameternames(6)='h_u'
        svalue(6)=0._DP
!  thickness of the plate in m
        parameternames(7)='d'
        svalue(7)=0.001_DP
!  material density in kg/m3
        parameternames(8)='rho'
        svalue(8)=1.205_DP
!  specific heat capacitance
        parameternames(9)='c_p'
        svalue(9)=1000._DP
!  power density
        parameternames(10)='p'
        svalue(10)=0._DP
!  phase of power density
        parameternames(11)='phase'
        svalue(11)=0._DP
!  temperature of the surrounding fluid (in Kelvin)
        parameternames(12)='T_F'
        svalue(12)=293.15_DP
!  velocity in x-direction
        parameternames(13)='velx'
        svalue(13)=0._DP
!  velocity in y-direction
        parameternames(14)='vely'
        svalue(14)=0._DP
!----------------------------------------------------------------------------------------------------------------------
      case ('3DTHERMOELECTRIC')
        numparam=23
        allocate(parameternames(numparam),svalue(numparam))
!  values for air
!  thermal conductivity in the first principle axis
        parameternames(1)='lam11'
        svalue(1)=0.02_DP
!  thermal conductivity in the second principle axis
        parameternames(2)='lam22'
        svalue(2)=0.02_DP
!  thermal conductivity in the third principle axis
        parameternames(3)='lam33'
        svalue(3)=0.02_DP
!  angle between the x-axis and the first principle axis
        parameternames(4)='angle'
        svalue(4)=0._DP
!  material density in kg/m3
        parameternames(5)='rho'
        svalue(5)=1.205_DP
!  specific heat capacitance
        parameternames(6)='c_p'
        svalue(6)=1000._DP
!  power density
        parameternames(7)='p'
        svalue(7)=0._DP
!  phase of power density
        parameternames(8)='phase'
        svalue(8)=0._DP
!  temperature of the surrounding fluid (in Kelvin)
        parameternames(9)='T_F'
        svalue(9)=293.15_DP
!  velocity in x-direction
        parameternames(10)='velx'
        svalue(10)=0._DP
!  velocity in y-direction
        parameternames(11)='vely'
        svalue(11)=0._DP
!  velocity in z-direction
        parameternames(12)='velz'
        svalue(12)=0._DP
!  electric conductivity in the first principle axis
        parameternames(13)='kappa11'
        svalue(13)=1._DP
!  electric conductivity in the second principle axis
        parameternames(14)='kappa22'
        svalue(14)=1._DP
!  electric conductivity in the third principle axis
        parameternames(15)='kappa33'
        svalue(15)=1._DP
!  angle between the x-axis and the third principle axis
        parameternames(16)='angle1'
        svalue(16)=0._DP
!  current source in A/m^3
        parameternames(17)='source'
        svalue(17)=0._DP
!  external current density in A/m^2
        parameternames(18)='J_x'
        svalue(18)=0._DP
        parameternames(19)='J_y'
        svalue(19)=0._DP
        parameternames(20)='J_z'
        svalue(20)=0._DP
!  thermoelectric coupling
        parameternames(21)='Peltier'
        svalue(21)=0._DP
!  seebeck        
        parameternames(22)='Seebeck'
        svalue(22)=0._DP
! Reference temperature in K
        parameternames(23)='Tref'
        svalue(23)=1._DP
!-----------------------------------------------------------------------------------------------------------------------
      case ('3DSTATCURRENT')
        numparam=3
        allocate(parameternames(numparam),svalue(numparam))
        parameternames(1)='kappa11'
        svalue(1)=1._DP
        parameternames(2)='kappa22'
        svalue(2)=1._DP
        parameternames(3)='kappa33'
       svalue(3)=1._DP
!----------------------------------------------------------------------------------------------------------------------   
      case ('3DELECTROSTATICS')
        numparam=4
        allocate(parameternames(numparam),svalue(numparam))
        parameternames(1)='epsr11'
        svalue(1)=1._DP
        parameternames(2)='epsr22'
        svalue(2)=1._DP
        parameternames(3)='epsr33'
        svalue(3)=1._DP
        parameternames(4)='rho'
        svalue(4)=0._DP
!       
!----------------------------------------------------------------------------------------------------------------------   
      case ('3DSTOKESFLOW')
        numparam=2
        allocate(parameternames(numparam),svalue(numparam))
        parameternames(1)='mu'
        svalue(1)=1.73E-5_DP !(value for AIR)
        parameternames(2)='rho'
        svalue(2)=1.229_DP   !(value for AIR)
!       
      case default
!  this one is only implemented for calculating w/o a physics mode
        numparam=1
        allocate(parameternames(numparam),svalue(numparam))
        parameternames(1)='avalue'
        svalue(1)=1._DP
      end select
!
!  convert parameter names to UPPERCASE
      do i=1,numparam
        call low2hi(parameternames(i),len_trim(parameternames(i)))
      end do
!  allocate the array of parameters and copy the standard values to it
      do i=0,maxmat
        allocate(param(i)%d(numparam))
        do j=1,numparam
          param(i)%d(j)=svalue(j)
        end do
      end do
      deallocate(svalue)
      return
      end subroutine setstandardvalues3D
