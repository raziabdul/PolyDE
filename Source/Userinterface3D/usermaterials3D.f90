   subroutine usermaterials3D(matname,found,xyzs,names,values,numnames, elem)
      use femtypes
      use varmat_handler,        only : getvarmatprop
      use matconstants
      use feminterface,          only : number2string, putmatparam
      implicit none
      real (DP),            optional :: xyzs(:), values(:)
      integer (I4B),        optional :: numnames
      integer (I4B),        optional :: elem
      character (len=*)              :: matname
      character (len=*),    optional :: names(:)
      logical                        :: found
      intent (in)                    :: matname, xyzs, elem
      intent (out)                   :: names, values, found, numnames
!
!    $Revision: 1.4 $
!    $Date: 2014/07/03 13:23:03 $
!    $Author: m_kasper $
!
!
!  User programming interface for material parameters which cannot be obtained 
!  by use of material files (inhomogenous, spatial or parameter dependent materials)
!
!  Input:
!        matname     the actual material, for which parameters should be returned
!        xyzs        (world-)coordinates for which the material values are to be returned
!        elem        element number
!  Output:
!        found       = .true. if the material is defined in this file
!        names       list of names of material parameter to be given
!        values      list of values for the corresponding materialparameters 
!        numnames    number of name/ values pairs beeing returned
!
!  Notes:
!  -  in the case of a inhomogenous (spatial dependency) materials use matvar=.true.
!  -  material names (matname) of usermaterials are case sensitive 
!     (unlike filenames of materialfiles in Windows) 
!  -  the possible material parameters are listed in the documetation (PDEcoeff) 
!     and defined in the subroutine setstandardvalues  
!  -  the list of possible material parameters depends on the physics mode
!  -  parameter names are not case sensitve (use upper case)
!  -  if a parameter is not defined below, the standard value is used
!  -  you should avoid to give more than the required parameters 
!  
!
!  local variables:
!
      real (DP)                      :: value
      character(len=50)              :: paramname
      logical                        :: set
      select case (matname)
!######################################  add materials here below #####################################

      case ('my-material')
        found=.true.
        if (.not.present(names)) return
!  permeability
        names(1)='velx'
        values(1)=-sqrt(xyzs(1)**2+xyzs(2)**2)*sin(atan2(xyzs(2),xyzs(1)))
        names(2)='vely'
        values(2)=sqrt(xyzs(1)**2+xyzs(2)**2)*cos(atan2(xyzs(2),xyzs(1)))
        names(3)='kappa'
        values(3)=25.e6
!
        numnames=3
      case ('heat-material')
        found=.true.
        if (.not.present(names)) return
!  unit conductivity
        names(1)='kx'
        values(1)=1._DP
        names(2)='ky'
        values(2)=1._DP
        names(3)='kz'
        values(3)=1._DP
!
        numnames=3
        
      case ('Pdoped_nSi')
          found = .true.
          if (.not.present(names)) return
          ! eps_r
          names(1)='epsr11'
          values(1)= 11.7_DP
          names(2)='epsr22'
          values(2)= 11.7_DP
          names(3)='epsr33'
          values(3)= 11.7_DP
          if (present(xyzs)) then ! Material position present
          ! N_A
            names(4) ='N_A'
            paramname = 'N_A'
            call getvarmatprop(matname,paramname,xyzs,value,set)
            if (set) then
              values(4) = value
            else
              values(4) = 2.E21_DP
            end if
              
            names(5) ='N_D'
            paramname = 'N_D'
            call getvarmatprop(matname,paramname,xyzs,value,set)
            if (set) then
              values(5) = value
            else
              values(5) = 1.E21_DP
            end if

          else ! No material position specified, use standard values
            ! N_A
            names(4)='N_A'
            values(4)= 2.E21_DP ! Basic doping of the substrate
            ! N_D
            names(5) ='N_D'
            values(5) = 1.E21_DP ! Basic doping of the substrate
          end if
          ! T_ref
          names(6) ='T_ref'
          values(6) = 300_DP
          ! mu_n
          names(7) ='mu_n'
          values(7) = .1400_DP
          ! mu_p
          names(8) ='mu_p'
          values(8) = .0450_DP
          ! n_i
          names(9) ='n_i'
          values(9) = 1.0E16_DP
          ! tau_n
          names(10) ='tau_n'
          values(10) = 1.E-11_DP
          ! tau_p
          names(11) ='tau_p'
          values(11) = 1.E-11_DP
          
          if (present(numnames)) numnames = 11
          
! __ CIS-pSi __ 
      case ('CIS-pSi-Resistor')
          found = .true.
          if (.not.present(names)) return
          ! eps_r
          names(1)='epsr11'
          values(1)= 11.7_DP
          names(2)='epsr22'
          values(2)= 11.7_DP
          names(3)='epsr33'
          values(3)= 11.7_DP
          ! N_A
          names(4)='N_A'
          if (present(xyzs)) then
            paramname = 'N_A'
            call getvarmatprop(matname,paramname,xyzs,value,set)
            values(4)= value
            if (.not.set) then
              values(4)= 0._DP ! If no specific value given, set N_D = 0 |=> No Doping
            end if
          else
            values(4)= 0._DP ! If no position given, set N_D = 0 |=> No Doping
          end if
          ! N_D
          names(5) ='N_D'
          values(5) = 1.E21_DP ! Basic doping of the substrate
          ! T_ref
          names(6) ='T_ref'
          values(6) = 300_DP
          ! mu_n
          names(7) ='mu_n'
          values(7) = .1400_DP
          ! mu_p
          names(8) ='mu_p'
          values(8) = .0450_DP
          ! n_i
          names(9) ='n_i'
          values(9) = 1.0E16_DP
          ! tau_n
          names(10) ='tau_n'
          values(10) = 1.E-11_DP
          ! tau_p
          names(11) ='tau_p'
          values(11) = 1.E-11_DP
          
          if (present(numnames)) numnames = 11
          
      case ('CIS-pSi-Resistor2')
          found = .true.
          if (.not.present(names)) return
          ! eps_r
          names(1)='epsr11'
          values(1)= 11.7_DP
          names(2)='epsr22'
          values(2)= 11.7_DP
          names(3)='epsr33'
          values(3)= 11.7_DP
          ! N_A
          names(4)='N_A'
          if (present(xyzs)) then
            paramname = 'N_A'
            call getvarmatprop(matname,paramname,xyzs,value,set)
            values(4)= value
            if (.not.set) then
              values(4)= 0._DP ! If no specific value given, set N_D = 0 |=> No Doping
            end if
          else
            values(4)= 0._DP ! If no position given, set N_D = 0 |=> No Doping
          end if
          ! N_D
          names(5) ='N_D'
          values(5) = 1.E21_DP ! Basic doping of the substrate
          ! T_ref
          names(6) ='T_ref'
          values(6) = 300_DP
          ! mu_n
          names(7) ='mu_n'
          values(7) = .1400_DP
          ! mu_p
          names(8) ='mu_p'
          values(8) = .0480_DP
          ! n_i
          names(9) ='n_i'
          values(9) = 1.0E16_DP
          ! tau_n
          names(10) ='tau_n'
          values(10) = 1.E-11_DP
          ! tau_p
          names(11) ='tau_p'
          values(11) = 1.E-11_DP
          ! C11
          names(12) = 'C11'
          values(12) = 169.E9_DP
          ! C12
          names(13) = 'C12'
          values(13) = 63.9E9_DP
          ! C13
          names(14) = 'C11'
          values(14) = 63.9E9_DP
          ! C14
          names(15) = 'C14'
          values(15) = 0._DP
          ! C15
          names(16) = 'C15'
          values(16) = 0._DP
          ! C16
          names(17) = 'C16'
          values(17) = 0._DP
          ! C22
          names(18) = 'C22'
          values(18) = 169.E9_DP
          ! C23
          names(19) = 'C23'
          values(19) = 63.9E9_DP
          ! C24
          names(20) = 'C24'
          values(20) = 0._DP
          ! C25
          names(21) = 'C25'
          values(21) = 0._DP
          ! C26
          names(22) = 'C26'
          values(22) = 0._DP
          ! C33
          names(23) = 'C33'
          values(23) = 169.E9_DP
          ! C34
          names(24) = 'C34'
          values(24) = 0._DP
          ! C35
          names(25) = 'C35'
          values(25) = 0._DP
          ! C36
          names(26) = 'C36'
          values(26) = 0._DP
          ! C44
          names(27) = 'C44'
          values(27) = 79.55E9_DP
          ! C45
          names(28) = 'C45'
          values(28) = 0._DP
          ! C46
          names(29) = 'C45'
          values(29) = 0._DP
          ! C55
          names(30) = 'C55'
          values(30) = 79.55E9_DP
          ! C56
          names(31) = 'C56'
          values(31) = 0._DP
          ! C66
          names(32) = 'C66'
          values(32) = 79.55E9_DP
          ! PI11
          names(33) = 'PI11'
          values(33) = 6.6E-11_DP
          ! PI12
          names(34) = 'PI12'
          values(34) = -1.1E-11_DP
          ! PI44
          names(35) = 'PI44'
          values(35) = 138.1E-11_DP
          ! fpa_x
          names(36) = 'fpa_x'
          values(36) = 0._DP
          ! fpa_y
          names(37) = 'fpa_y'
          values(37) = 0._DP
          ! fpa_z
          names(38) = 'fpa_z'
          values(38) = 1._DP
          ! spa_x
          names(39) = 'spa_x'
          values(39) = 1._DP
          ! spa_y
          names(40) = 'spa_y'
          values(40) = 1._DP
          ! spa_z
          names(41) = 'spa_z'
          values(41) = 0._DP
        
          if (present(numnames)) numnames = 41
          
      case default
!  if the material is not defined in the file
        found=.false.
        return
      end select
!
      return
      end subroutine usermaterials3D
