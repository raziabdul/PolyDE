      subroutine initialize3D()
      use femtypes
      use feminterface,       only: zeit, getsetting, readsetting, setprocesspriority, physicsinfo
      !use feminterface3D,     only: 
      use globalvariables3D,  only: c0, eltype, nnat, omega, physics, pi, whatsystem
      implicit none
!
!    $Revision: 1.6 $
!    $Date: 2015/04/01 10:50:31 $
!    $Author: juryanatzki $
!
!  perform initial settings
!
      integer (I4B) ::  priority, n
      real (DP) :: freq, lam0
!      real (DP) :: xmin, xmax, ymin, ymax, h
!
!      print*,'hihi'
!  read in system info (Windows or Linux) and store it under the global variable
!  whatsystem
      call getsetting('WHATSYSTEM',whatsystem)
!
      call zeit(' ')
!
      priority=1
      call setprocesspriority(priority)
!
!  read physics mode from FEMsettings.txt
      call getsetting('PHYSICS_MODE',physics)
      call physicsinfo(physics,nnat,3)
!      print*,'in initialize3D, physics = ',physics
!
!  call settings for wavelength and frequency
      call getsetting('WAVELENGTH',lam0)
      call getsetting('FREQUENCY',freq)
!  if one value is 0, it is not used. frequency is stronger than wavelength
      if (freq .ge. 0._DP)  then
        omega = 2._DP*pi*freq
      else if (lam0 .gt. 0._DP) then
        omega = (c0 / lam0)*2._DP*pi
      else
        omega = 0._DP
      end if
!
!  get element type 
      call getsetting('ELEMENT_TYPE', eltype)
      call low2hi(eltype, 16)

      return
      end subroutine initialize3D
