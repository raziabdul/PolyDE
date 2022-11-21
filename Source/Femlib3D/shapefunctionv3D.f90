      subroutine shapefunctionv3D(lambda, vert, polylo, polyhi,   &
     &  nff, curl, xsi, cxsi, errcode)
      use feminterface3D, only: shapev3D, curlshapev3D, getgl
      use femtypes
      use globalvariables3D, only: polymaxv
      implicit none
      integer (I4B) :: polylo, polyhi, nff, errcode
      real (DP) :: vert(3,4), lambda(4)
      real (DP) :: xsi(:,:), cxsi(:,:)
      logical curl
      intent (in) :: lambda, vert, polylo, polyhi, nff, curl
      intent (out) :: xsi, cxsi, errcode
!
!    $Revision: 1.15 $
!    $Date: 2014/08/22 10:49:46 $
!    $Author: m_kasper $
!
!  evaluate the shape function and its curl at a point given in natural 
!  coordinates
!
!-----------------------------------------------------------------------
!
!  The barycentric coordinates are defined in the tetrahedron of specific 
!  orientation--See "tet_numbering.pdf" and "Barycentric Coordinates.doc"
!  
!-----------------------------------------------------------------------
!  
!  Input:
!            lambda            natural (barycentric) coordinates of the point 
!            v1, v2, v3, v4    coordinates of the tetrahedron vertices
!            polylo      lowest polynomial degree to consider                        
!            polyhi      highest polynomial degree                                   
!            nff         number of formfunctions                                     
!            curl        if .true. the curl will be computed                     
!  Output:                                                                           
!            xsi         vector of shapefunction at lambda                           
!            cxsi        curl of the shapefunction at xi, yi, zi                     
!                        where cxsi(i,:) representing the i-component, i=x,y,z 
!  Errorcodes:                                                                       
!            errcode     2001 = polylo is greater than polyhi                        
!                        2002 = polynomial degree cannot be higher than polymax      
!                        4001 = unavailable element type (currently only Nedelec)
!
!  global variables used: 
!            polymaxv
!         
!  local variables
      real(DP) :: gl(3,4)
!
!  Check values of polylo and polyhi
      if ( (polylo .gt. polymaxv) .or. (polyhi .gt. polymaxv) ) then
         errcode = 2002
         return
      end if
      if ( polylo .gt. polyhi ) then
        errcode = 2001
        return
      end if
!
!  Calculation of the gradients 
!
      call getgl(vert, gl)
!
      if (curl) then
!  Calculate curl vector for each component. The last argument is 
!  to identify the requested component
!
! TO DO : make it call only once and return all components
!
! x-component
        call curlshapev3D(lambda, gl, polylo, polyhi, cxsi(1,:), 1)
! y-component
        call curlshapev3D(lambda, gl, polylo, polyhi, cxsi(2,:), 2)
! z-component
        call curlshapev3D(lambda, gl, polylo, polyhi, cxsi(3,:), 3)
      end if
!
!  Calculate the shape function vector for each component
!
! TO DO : make it call only once and return all components
!
!   x-component
      call shapev3D(lambda, gl(1,:), polylo, polyhi, xsi(1,:))
!   y-component
      call shapev3D(lambda, gl(2,:), polylo, polyhi, xsi(2,:))
!   z-component
      call shapev3D(lambda, gl(3,:), polylo, polyhi, xsi(3,:))
!
      end subroutine shapefunctionv3D
!
!
!
      subroutine shapev3D(l, gl, polylo, polyhi, vec, eltype )
      Use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi
      real (DP)     :: l(4), gl(4)
      real (DP)     :: vec(:)
      character (len = *) :: eltype
      intent (in)   :: l, polylo, polyhi, eltype
      intent (out)  :: vec
!  
!  evaluate the shape functions for the corresponding polynomial degree
!  
!  Input:
!            l        coordinates to be calculated 
!            gl       component of gradient of triangle coordinates (x,y, or z) 
!            polylo   lowest polynomial degree to consider
!            polyhi   highest polynomial degree 
!  Output:
!            vec      array of shape functions at coordinates. 
!                     contains either x,y, or z component depending
!                     on the caller function
!
!  The polylo and polyhi variables are used to select the shapefunctions. If polylo .eq.
!  polyhi, only the new shape functions of the corresponding polynomial degree p are cal-
!  culated. That means: if polylo .eq. polyhi = 7 ==> xsi(29) to xsi(36) are calculated.
!  If polylo .neq. polyhi, shape functions starting from polylo to polyhi are calculated.
!  That means: if polylo = 2 and polyhi = 7 ==> xsi(4) up to xsi(36) are calculated.
!
!--------------------------------------------------------------------------
! LATEST SHAPE FUNCTIONS:
!
!  Schoeberl Functions 
!--------------------------------------------------------------------------
!
!
!  local variables:
      integer (I4B) k, p
      real (DP) L1, L2, L3, L4, gL1, gL2, gL3, gL4
!
!  intermediate variables to avoid arrays of l and gl
      L1 = l(1)
      L2 = l(2)
      L3 = l(3)
      L4 = l(4)
      gL1 = gl(1)
      gL2 = gl(2)
      gL3 = gl(3)
      gL4 = gl(4)
      
      p = polylo

! k depends on the formulations of shape functions used for the 3D
! vector bases. This is Nedelec type
      if (polylo .gt. 1) then
        k = ( p-1 )*( p+1 )*( p+2 )/2
      else
        k = 0
      end if
!
!-------------------------------------------------------------------!
!  BEGIN LOOP TO CALCULATE SHAPE FUNCTIONS FOR POLYNOMIAL DEGREE p  !
!-------------------------------------------------------------------!
      do
        if (p .gt. polyhi) exit    ! stopping criteria
!  Begin "case" control structure to select the appropriate shape functions for p
        poly: select case (p)
	case (1)   !  shape functions for p = 1
		vec(1-k)=L1*gL2-L2*gL1

		vec(2-k)=L2*gL3-L3*gL2

		vec(3-k)=L1*gL3-L3*gL1

		vec(4-k)=L1*gl4-L4*gL1

		vec(5-k)=L2*gl4-L4*gL2

		vec(6-k)=L3*gl4-L4*gL3

	case (2)   !  shape functions for p = 2
		vec(7-k)=L2*gL1+L1*gL2

		vec(8-k)=L3*gL2+L2*gL3

		vec(9-k)=L3*gL1+L1*gL3

		vec(10-k)=L4*gL1+L1*gl4

		vec(11-k)=L4*gL2+L2*gl4

		vec(12-k)=L4*gL3+L3*gl4

		vec(13-k)=gL2*L4*L3+(L4*gL3-L3*gl4)*L2

		vec(14-k)=-(L2*gL3-L3*gL2)*L4

		vec(15-k)=gL1*L4*L3+(L4*gL3-L3*gl4)*L1

		vec(16-k)=-(L1*gL3-L3*gL1)*L4

		vec(17-k)=L2*gL1*L4+(L4*gL2-L2*gl4)*L1

		vec(18-k)=-(L1*gL2-L2*gL1)*L4

		vec(19-k)=gL1*L3*L2+(L3*gL2-L2*gL3)*L1

		vec(20-k)=-(L1*gL2-L2*gL1)*L3

	case (3)   !  shape functions for p = 3
		vec(21-k)=-L2**2*gL1+(2._DP*(gL1-gL2)*L2+L1*gL2)*L1

		vec(22-k)=-L3**2*gL2+(2._DP*(gL2-gL3)*L3+L2*gL3)*L2

		vec(23-k)=-L3**2*gL1+(2._DP*(gL1-gL3)*L3+L1*gL3)*L1

		vec(24-k)=-L4**2*gL1+(2._DP*(gL1-gl4)*L4+L1*gl4)*L1

		vec(25-k)=-L4**2*gL2+(2._DP*(gL2-gl4)*L4+L2*gl4)*L2

		vec(26-k)=-L4**2*gL3+(2._DP*(gL3-gl4)*L4+L3*gl4)*L3

		vec(27-k)=gL2*L4*L3+(L4*gL3+L3*gl4)*L2

		vec(28-k)=-L4*L2*(L2-L4)*gL3-L4*L3*(L3-L4)*gL2+L3*L2*(L3-2._DP*L4+L2)*gl4

		vec(29-k)=-gL2*L4*L3**2+((2._DP*(gL2-gL3)*L4+L3*gl4)*L3+(L4*gL3-L3*gl4)*L2)*L2

		vec(30-k)=(L2*gL3-L3*gL2)*L4*(-L4+L2+L3)

		vec(31-k)=gL1*L4*L3+(L4*gL3+L3*gl4)*L1

		vec(32-k)=-L4*L3*(L3-L4)*gL1-L4*L1*(L1-L4)*gL3+L3*L1*(L3-2._DP*L4+L1)*gl4

		vec(33-k)=-gL1*L4*L3**2+((2._DP*(gL1-gL3)*L4+L3*gl4)*L3+(L4*gL3-L3*gl4)*L1)*L1

		vec(34-k)=(L1*gL3-L3*gL1)*L4*(-L4+L1+L3)

		vec(35-k)=L2*gL1*L4+(L4*gL2+L2*gl4)*L1

		vec(36-k)=(-L4*(L1-L4)*gL2+L2*(-2._DP*L4+L2+L1)*gl4)*L1-L4*(L2-L4)*L2*gL1

		vec(37-k)=-L2**2*gL1*L4+((2._DP*(gL1-gL2)*L4+L2*gl4)*L2+(L4*gL2-L2*gl4)*L1)*L1

		vec(38-k)=(L1*gL2-L2*gL1)*L4*(-L4+L1+L2)

		vec(39-k)=gL1*L3*L2+(L3*gL2+L2*gL3)*L1

		vec(40-k)=(L2*(-2._DP*L3+L2+L1)*gL3-L3*(L1-L3)*gL2)*L1-L3*(L2-L3)*L2*gL1

		vec(41-k)=-gL1*L3*L2**2+(L2*(2._DP*(gL1-gL2)*L3+L2*gL3)+(L3*gL2-L2*gL3)*L1)*L1

		vec(42-k)=(L1*gL2-L2*gL1)*L3*(-L3+L1+L2)

		vec(43-k)=gL1*L4*L3*L2+(gL2*L4*L3+(L3*gl4-L4*gL3)*L2)*L1

		vec(44-k)=gL1*L4*L3*L2+(gL2*L4*L3+(-L4*gL3-L3*gl4)*L2)*L1

		vec(45-k)=-(L1*gL2-L2*gL1)*L3*L4

	case (4)   !  shape functions for p = 4
		vec(46-k)=L2**3*gL1+((3._DP*gL2-8._DP*gL1)*L2**2+((3._DP*gL1-8._DP*gL2)*L2+L1*gL2)*L1)*L1

		vec(47-k)=L3**3*gL2+((3._DP*gL3-8._DP*gL2)*L3**2+((3._DP*gL2-8._DP*gL3)*L3+L2*gL3)*L2)*L2

		vec(48-k)=L3**3*gL1+(-(8._DP*gL1-3._DP*gL3)*L3**2+((3._DP*gL1-8._DP*gL3)*L3+L1*gL3)*L1)*L1

		vec(49-k)=L4**3*gL1+(-(8._DP*gL1-3._DP*gl4)*L4**2+((-8._DP*gl4+3._DP*gL1)*L4+L1*gl4)*L1)*L1

		vec(50-k)=L4**3*gL2+((3._DP*gl4-8._DP*gL2)*L4**2+((-8._DP*gl4+3._DP*gL2)*L4+L2*gl4)*L2)*L2

		vec(51-k)=L4**3*gL3+((3._DP*gl4-8._DP*gL3)*L4**2+((-8._DP*gl4+3._DP*gL3)*L4+L3*gl4)*L3)*L3

		vec(52-k)=-L4*L2*(L2-L4)*gL3-2._DP*(-gl4+gL2+gL3)*L4*L3*L2-L3*L2*(L3+L2)*gl4-L4*L3*(L3-L4)*gL2

		vec(53-k)=-gL2*L4*L3**2+((2._DP*(gL2-gL3)*L4-L3*gl4)*L3+(L4*gL3+L3*gl4)*L2)*L2

		vec(54-k)=(-gL2+8._DP*gl4)*L4*L3*L2**2+(8._DP*gl4-gL3)*L4*L3**2*L2+L4*L2*(L4**2+(-4._DP*L4+L2)*L2)*gL3-L3*L2*(4._DP*L4**2+(-2._DP+2._DP&
		*L1)*L4+(L1-1._DP)**2)*gl4+L4*L3*(L4**2+(-4._DP*L4+L3)*L3)*gL2

		vec(55-k)=L3**2*L4*(L3-L4)*gL2+((2._DP*(gL2-gL3)*L4**2+(-L4*(3._DP*gL2+gL1-gl4)-L3*gl4)*L3)*L3+(L4*L3*(gL1+3._DP*gL3-gl4)+L4**2*gL3+(L3&
		*gl4-L4*gL3)*L2)*L2)*L2

		vec(56-k)=L3**3*gL2*L4+(((3._DP*gL3-8._DP*gL2)*L4-L3*gl4)*L3**2+(((3._DP*gL2-8._DP*gL3)*L4+4._DP*L3*gl4)*L3+(L4*gL3-L3*gl4)*L2)*L2)*L2

		vec(57-k)=-(L2*gL3-L3*gL2)*L4*(6._DP*L4**2+(6._DP*L1-6._DP)*L4+(L1-1._DP)**2)

		vec(58-k)=-L4*L3*(L3-L4)*gL1+(L4**2*gL3+(2._DP*L4*(-gL1-gL3+gl4)-L3*gl4)*L3+(-L4*gL3-L3*gl4)*L1)*L1

		vec(59-k)=-gL1*L4*L3**2+((2._DP*(gL1-gL3)*L4-L3*gl4)*L3+(L4*gL3+L3*gl4)*L1)*L1

		vec(60-k)=L4*L3**2*L1*(8._DP*gl4-gL3)+L4*L1*(L4**2+(-4._DP*L4+L1)*L1)*gL3-L1**2*(-8._DP*gl4+gL1)*L4*L3-L3*L1*(4._DP*L4**2+(-2._DP+2._DP&
		*L2)*L4+(-1._DP+L2)**2)*gl4+L4**3*gL1*L3-4._DP*L4**2*L3**2*gL1+L3**3*gL1*L4

		vec(61-k)=L4*gL1*(L3-L4)*L3**2+((2._DP*(gL1-gL3)*L4**2+(-L4*(gL2+3._DP*gL1-gl4)-L3*gl4)*L3)*L3+(L4**2*gL3+L4*L3*(3._DP*gL3+gL2-gl4)+(L3&
		*gl4-L4*gL3)*L1)*L1)*L1

		vec(62-k)=L3**3*gL1*L4+((-(8._DP*gL1-3._DP*gL3)*L4-L3*gl4)*L3**2+(((3._DP*gL1-8._DP*gL3)*L4+4._DP*L3*gl4)*L3+(L4*gL3-L3*gl4)*L1)*L1)*L1

		vec(63-k)=-(L1*gL3-L3*gL1)*L4*(6._DP*L4**2-6._DP*L4+1._DP+(6._DP*L4-2._DP+L2)*L2)

		vec(64-k)=(-L4*(2._DP*L2+L1-L4)*gL2-L2*(-2._DP*L4+L2+L1)*gl4)*L1-L4*L2*(L2-L4+2._DP*L1)*gL1

		vec(65-k)=-L2**2*gL1*L4+((2._DP*(gL1-gL2)*L4-L2*gl4)*L2+(L4*gL2+L2*gl4)*L1)*L1

		vec(66-k)=(((-3._DP*L4**2+(8._DP*L4-L2)*L2)*L2-L2*(L2-9._DP*L4-L3+1._DP)*L1)*gl4+L4*(-L2**2+L4**2+(-4._DP*L4+L1)*L1)*gL2)*L1-L4*L2*(-L4&
		**2+(4._DP*L4-L2)*L2+L1**2)*gL1

		vec(67-k)=(L2*(L1-L2)*(-2._DP*L4+L2+L1)*gl4-L4*((2._DP*L4-L2)*L2+(-L4-2._DP*L2+L1)*L1)*gL2)*L1-L4*L2*((L4-L2)*L2+(L2-3._DP*L4-L3+1._DP)&
		*L1)*gL1

		vec(68-k)=L2**3*gL1*L4+((-L2*gl4+(3._DP*gL2-8._DP*gL1)*L4)*L2**2+(L2*(4._DP*L2*gl4+(3._DP*gL1-8._DP*gL2)*L4)+(L4*gL2-L2*gl4)*L1)*L1)*L1

		vec(69-k)=-(L1*gL2-L2*gL1)*L4*(6._DP*L4**2-6._DP*L4+1._DP+(6._DP*L4-2._DP+L3)*L3)

		vec(70-k)=-L3*(L2-L3)*L2*gL1+L1*(L2*(L3*(4._DP*gL3+2._DP*gl4)-L2*gL3)+(-L3*gL2-L2*gL3)*L1+L3**2*gL2)

		vec(71-k)=-gL1*L3*L2**2+(L2*(2._DP*(gL1-gL2)*L3-L2*gL3)+(L3*gL2+L2*gL3)*L1)*L1

		vec(72-k)=(((-3._DP*L3**2+(8._DP*L3-L2)*L2)*L2-L2*(L2-9._DP*L3-L4+1._DP)*L1)*gL3+L3*(L3**2-L2**2+(-4._DP*L3+L1)*L1)*gL2)*L1-L3*L2*(-L3*&
		*2+(4._DP*L3-L2)*L2+L1**2)*gL1

		vec(73-k)=L3*gL1*(L2-L3)*L2**2+((2._DP*(gL1-gL2)*L3**2+L2*(L3*(-3._DP*gL1+gL3-gl4)-L2*gL3))*L2+(L3*L2*(3._DP*gL2-gL3+gl4)+(L2*gL3-L3*gL&
		2)*L1+L3**2*gL2)*L1)*L1

		vec(74-k)=L2**3*gL1*L3+(((3._DP*gL2-8._DP*gL1)*L3-L2*gL3)*L2**2+(((3._DP*gL1-8._DP*gL2)*L3+4._DP*L2*gL3)*L2+(L3*gL2-L2*gL3)*L1)*L1)*L1

		vec(75-k)=-(L1*gL2-L2*gL1)*L3*(1._DP+L4**2-2._DP*L4+6._DP*(-L2-L1)*L3)

		vec(76-k)=gL1*L4*L3*L2+(gL2*L4*L3+(L4*gL3+L3*gl4)*L2)*L1

		vec(77-k)=(L3*(2._DP*L4**2-L4)*gL2+((-2._DP*L4**2+L4)*gL3+L3*(L4*(4._DP*gl4+gL1)-gl4))*L2)*L1-L4*L3*L2*(1._DP+L1-2._DP*L4)*gL1

		vec(78-k)=L1*(L3*(2._DP*L4**2-L4)*gL2+((-2._DP*L4**2+L4)*gL3+L3*(L4*(-4._DP*gl4-gL1)+gl4))*L2)-L3*L4*(-L4+L2+L3)*L2*gL1

		vec(79-k)=-gL1*L4*L3*(L2-L3)*L2+(gL2*L4*L3**2+((-2._DP*L4*gL3+L3*gl4)*L3+(L4*gL3-L3*gl4)*L2)*L2+(-gL2*L4*L3+(L4*gL3-L3*gl4)*L2)*L1)*L1

		vec(80-k)=-gL1*L4*L3*(L2-L3)*L2+(gL2*L4*L3**2+((-2._DP*L4*gL3-L3*gl4)*L3+(L4*gL3+L3*gl4)*L2)*L2+(-gL2*L4*L3+(L4*gL3+L3*gl4)*L2)*L1)*L1

		vec(81-k)=-gL1*L4*L3*L2**2+((2._DP*L3*(gL1-gL2)*L4+(L4*gL3-L3*gl4)*L2)*L2+(gL2*L4*L3+(L3*gl4-L4*gL3)*L2)*L1)*L1

		vec(82-k)=-gL1*L4*L3*L2**2+((2._DP*L3*(gL1-gL2)*L4+(L4*gL3+L3*gl4)*L2)*L2+(gL2*L4*L3+(-L4*gL3-L3*gl4)*L2)*L1)*L1

		vec(83-k)=(L1*gL2-L2*gL1)*L3*(-2._DP*L4+1._DP)*L4

		vec(84-k)=(L1*gL2-L2*gL1)*L3*(-L3+L1+L2)*L4

	case (5)   !  shape functions for p = 5
		vec(85-k)=-L2**4*gL1+((18._DP*gL1-4._DP*gL2)*L2**3+(-27._DP*(gL1-gL2)*L2**2+((4._DP*gL1-18._DP*gL2)*L2+L1*gL2)*L1)*L1)*L1

		vec(86-k)=-L3**4*gL2+((18._DP*gL2-4._DP*gL3)*L3**3+(-27._DP*(gL2-gL3)*L3**2+((4._DP*gL2-18._DP*gL3)*L3+L2*gL3)*L2)*L2)*L2

		vec(87-k)=-L3**4*gL1+((18._DP*gL1-4._DP*gL3)*L3**3+(-27._DP*(gL1-gL3)*L3**2+L1*(L1*gL3+(4._DP*gL1-18._DP*gL3)*L3))*L1)*L1

		vec(88-k)=-L4**4*gL1+((18._DP*gL1-4._DP*gl4)*L4**3+(-27._DP*(gL1-gl4)*L4**2+L1*(L1*gl4+(4._DP*gL1-18._DP*gl4)*L4))*L1)*L1

		vec(89-k)=-L4**4*gL2+((18._DP*gL2-4._DP*gl4)*L4**3+(-27._DP*(gL2-gl4)*L4**2+((4._DP*gL2-18._DP*gl4)*L4+L2*gl4)*L2)*L2)*L2

		vec(90-k)=-L4**4*gL3+((18._DP*gL3-4._DP*gl4)*L4**3+(-27._DP*(gL3-gl4)*L4**2+((4._DP*gL3-18._DP*gl4)*L4+L3*gl4)*L3)*L3)*L3

		vec(91-k)=L4*L3*(L4**2+(-4._DP*L4+L3)*L3)*gL2+(L4**3*gL3+(L4**2*(11._DP*gl4+8._DP*gL1)+(L4*(-12._DP*gl4-4._DP*gL1-gL3)+L3*gl4)*L3)*L3+(&
		-4._DP*L4**2*gL3+(L4*(-11._DP*gl4-3._DP*gL1+gL3)+2._DP*L3*gl4)*L3+(L4*gL3+L3*gl4)*L2)*L2)*L2

		vec(92-k)=L3**2*L4*(L3-L4)*gL2+((2._DP*(gL2-gL3)*L4**2+((3._DP*gL3-2._DP*gl4)*L4+L3*gl4)*L3)*L3+(-(-2._DP*gl4+3._DP*gL2)*L4*L3+L4**2*gL&
		3+(-L4*gL3-L3*gl4)*L2)*L2)*L2

		vec(93-k)=L3**3*gL2*L4+(((3._DP*gL3-8._DP*gL2)*L4+L3*gl4)*L3**2+(((3._DP*gL2-8._DP*gL3)*L4-4._DP*L3*gl4)*L3+(L4*gL3+L3*gl4)*L2)*L2)*L2

		vec(94-k)=-L4*L3*(L3-L4)*(L4**2+(-8._DP*L4+L3)*L3)*gL2+(L4**4*gL3+(-4._DP*L4**3*gl4+(-9._DP*(gL3-3._DP*gl4)*L4**2+(2._DP*(-9._DP*gl4+gL&
		3)*L4+L3*gl4)*L3)*L3)*L3+(-9._DP*L4**3*gL3+(-9._DP*(gL2-3._DP*gl4)*L4**2+(3._DP*L4*(-gL1-13._DP*gl4)+3._DP*L3*gl4)*L3)*L3+(9._DP*L4**2*&
		gL3+(2._DP*(-9._DP*gl4+gL2)*L4+3._DP*L3*gl4)*L3+(L3*gl4-L4*gL3)*L2)*L2)*L2)*L2

		vec(95-k)=-L4*L3**2*(L4**2+(-4._DP*L4+L3)*L3)*gL2+((2._DP*(gL2-gL3)*L4**3+(-L4**2*(3._DP*gL1+11._DP*gL2-gL3)+(2._DP*(gL2-4._DP*gl4)*L4+&
		L3*gl4)*L3)*L3)*L3+(L4**3*gL3+(L4**2*(3._DP*gL1-gL2+11._DP*gL3)+(3._DP*(gL2-gL3)*L4+L3*gl4)*L3)*L3+(L4*(-4._DP*L4-2._DP*L3+L2)*gL3-L3*(&
		-8._DP*L4+L3+L2)*gl4)*L2)*L2)*L2

		vec(96-k)=-L3**3*L4*(L3-L4)*gL2+((-(-3._DP*gL3+8._DP*gL2)*L4**2+(2._DP*L4*(5._DP*gL2+gL1)+L3*gl4)*L3)*L3**2+(((3._DP*gL2-8._DP*gL3)*L4*&
		*2+(L4*(7._DP*gl4-gL1)-3._DP*L3*gl4)*L3)*L3+(L4**2*gL3+(-2._DP*L4*(-5._DP*gL3-gL1)-3._DP*L3*gl4)*L3+(L3*gl4-L4*gL3)*L2)*L2)*L2)*L2

		vec(97-k)=-L4*L3**4*gL2+(((18._DP*gL2-4._DP*gL3)*L4+L3*gl4)*L3**3+((-27._DP*(gL2-gL3)*L4-9._DP*L3*gl4)*L3**2+(((4._DP*gL2-18._DP*gL3)*L&
		4+9._DP*L3*gl4)*L3+(L4*gL3-L3*gl4)*L2)*L2)*L2)*L2

		vec(98-k)=(L2*gL3-L3*gL2)*L4*(-L4+L2+L3)*(10._DP*L4**2+(10._DP*L1-10._DP)*L4+(L1-1._DP)**2)

		vec(99-k)=L4*gL1*(L4**2+(-4._DP*L4+L3)*L3)*L3+(L4**3*gL3+(L4**2*(8._DP*gL2+11._DP*gl4)+(L4*(-11._DP*gl4+gL1-3._DP*gL2)+L3*gl4)*L3)*L3+(&
		-4._DP*L4**2*gL3+(-L4*(12._DP*gl4+gL1+4._DP*gL2)+2._DP*L3*gl4)*L3+(L4*gL3+L3*gl4)*L1)*L1)*L1

		vec(100-k)=L4*gL1*(L3-L4)*L3**2+((2._DP*(gL1-gL3)*L4**2+((3._DP*gL3-2._DP*gl4)*L4+L3*gl4)*L3)*L3+(-(-2._DP*gl4+3._DP*gL1)*L4*L3+L4**2*g&
		L3+(-L4*gL3-L3*gl4)*L1)*L1)*L1

		vec(101-k)=L3**3*gL1*L4+((-(8._DP*gL1-3._DP*gL3)*L4+L3*gl4)*L3**2+(((3._DP*gL1-8._DP*gL3)*L4-4._DP*L3*gl4)*L3+(L4*gL3+L3*gl4)*L1)*L1)*L&
		1

		vec(102-k)=-L4*gL1*(L3-L4)*(L4**2+(-8._DP*L4+L3)*L3)*L3+(L4**4*gL3+(-4._DP*L4**3*gl4+(9._DP*L4**2*(3._DP*gl4-gL3)+(-2._DP*L4*(9._DP*gl4&
		-gL3)+L3*gl4)*L3)*L3)*L3+(-9._DP*L4**3*gL3+(-9._DP*(-3._DP*gl4+gL1)*L4**2+(-3._DP*L4*(gL2+13._DP*gl4)+3._DP*L3*gl4)*L3)*L3+(9._DP*L4**2&
		*gL3+(2._DP*(-9._DP*gl4+gL1)*L4+3._DP*L3*gl4)*L3+(L3*gl4-L4*gL3)*L1)*L1)*L1)*L1

		vec(103-k)=-L4*gL1*(L4**2+(-4._DP*L4+L3)*L3)*L3**2+((2._DP*(gL1-gL3)*L4**3+(-(12._DP*gL1+4._DP*gL2+gl4)*L4**2+(2._DP*(gL1-4._DP*gl4)*L4&
		+L3*gl4)*L3)*L3)*L3+(L4**3*gL3+((4._DP*gL2+12._DP*gL3+gl4)*L4**2+(3._DP*(gL1-gL3)*L4+L3*gl4)*L3)*L3+(L4*(-4._DP*L4-2._DP*L3+L1)*gL3-L3*&
		(-8._DP*L4+L3+L1)*gl4)*L1)*L1)*L1

		vec(104-k)=-L4*gL1*(L3-L4)*L3**3+((-(8._DP*gL1-3._DP*gL3)*L4**2+(2._DP*L4*(5._DP*gL1+gL2)+L3*gl4)*L3)*L3**2+(((3._DP*gL1-8._DP*gL3)*L4*&
		*2+(-L4*(gL2-7._DP*gl4)-3._DP*L3*gl4)*L3)*L3+(L4**2*gL3+(2._DP*L4*(gL2+5._DP*gL3)-3._DP*L3*gl4)*L3+(L3*gl4-L4*gL3)*L1)*L1)*L1)*L1

		vec(105-k)=-L4*L3**4*gL1+(((18._DP*gL1-4._DP*gL3)*L4+L3*gl4)*L3**3+((-27._DP*(gL1-gL3)*L4-9._DP*L3*gl4)*L3**2+L1*(L3*(9._DP*L3*gl4+(4._&
		DP*gL1-18._DP*gL3)*L4)+(L4*gL3-L3*gl4)*L1))*L1)*L1

		vec(106-k)=(L1*gL3-L3*gL1)*L4*(-L4+L1+L3)*(10._DP*L4**2-10._DP*L4+1._DP+(10._DP*L4-2._DP+L2)*L2)

		vec(107-k)=(((3._DP*L4**2+(-8._DP*L4+L2)*L2)*L2+L2*(L2-9._DP*L4-L3+1._DP)*L1)*gl4+L4*(L4**2+(-8._DP*L4+3._DP*L2)*L2+(4._DP*L2-4._DP*L4+&
		L1)*L1)*gL2)*L1+L4*L2*(L4**2+(-4._DP*L4+L2)*L2+(L2-11._DP*L4-3._DP*L3+3._DP)*L1)*gL1

		vec(108-k)=(-L2*(L1-L2)*(-2._DP*L4+L2+L1)*gl4-L4*((2._DP*L4-3._DP*L2)*L2+(L1-L4)*L1)*gL2)*L1-L4*L2*((L4-L2)*L2+(-2._DP*L4+3._DP*L1)*L1)&
		*gL1

		vec(109-k)=L2**3*gL1*L4+(L2**2*(L2*gl4+(3._DP*gL2-8._DP*gL1)*L4)+(L2*(-4._DP*L2*gl4+(3._DP*gL1-8._DP*gL2)*L4)+(L4*gL2+L2*gl4)*L1)*L1)*L&
		1

		vec(110-k)=(((-4._DP*L4**3+(27._DP*L4**2+(-18._DP*L4+L2)*L2)*L2)*L2+((27._DP*L4**2+(-36._DP*L4+3._DP*L2)*L2)*L2+L2*(3._DP*L2-18._DP*L4+&
		L1)*L1)*L1)*gl4+(L4**4+L4*(-9._DP*L4+2._DP*L2)*L2**2-L4*(9._DP*L4**2-3._DP*L2**2+(-9._DP*L4+L1)*L1)*L1)*gL2)*L1+L4*L2*(L4**3+(-9._DP*L4&
		**2+(9._DP*L4-L2)*L2)*L2+(-11._DP*L4+L2-2._DP*L3+2._DP)*L1**2)*gL1

		vec(111-k)=(-L2*(L1-L2)*(12._DP*L4**2-10._DP*L4+1._DP+(10._DP*L4-2._DP+L3)*L3)*gl4+(2._DP*L4**2*(-L4+2._DP*L2)*L2+L4*(L4**2+(8._DP*L4-3&
		._DP*L2)*L2+(-4._DP*L4-2._DP*L2+L1)*L1)*L1)*gL2)*L1+(-L4*(L4**2+(-4._DP*L4+L2)*L2)*L2**2+L4*L2*(2._DP*L4**2+(-8._DP*L4+2._DP*L2)*L2+(-4&
		._DP*L4+3._DP*L2)*L1)*L1)*gL1

		vec(112-k)=-L4*gL1*(L2-L4)*L2**3+(L1*(L1*((L2*gl4-L4*gL2)*L1+L2*(L4*(10._DP*gL2+2._DP*gL3)-3._DP*L2*gl4)+L4**2*gL2)+L2*(L2*(L4*(7._DP*g&
		l4-gL3)-3._DP*L2*gl4)+(3._DP*gL1-8._DP*gL2)*L4**2))+L2**2*(L2*(2._DP*L4*(5._DP*gL1+gL3)+L2*gl4)+(3._DP*gL2-8._DP*gL1)*L4**2))*L1

		vec(113-k)=-L4*L2**4*gL1+(L2**3*(L2*gl4+(18._DP*gL1-4._DP*gL2)*L4)+((-27._DP*(gL1-gL2)*L4-9._DP*L2*gl4)*L2**2+(L2*(9._DP*L2*gl4+(4._DP*&
		gL1-18._DP*gL2)*L4)+(L4*gL2-L2*gl4)*L1)*L1)*L1)*L1

		vec(114-k)=(L1*gL2-L2*gL1)*L4*(-L4+L1+L2)*(10._DP*L4**2-10._DP*L4+1._DP+(10._DP*L4-2._DP+L3)*L3)

		vec(115-k)=L3*gL1*(L3**2+(-4._DP*L3+L2)*L2)*L2+(L2*(L3**2*(11._DP*gL3+8._DP*gl4)+L2*(L3*(-gL2-12._DP*gL3-4._DP*gl4)+L2*gL3))+L1*(L2*(L3&
		*(gL2-11._DP*gL3-3._DP*gl4)+2._DP*L2*gL3)+(L3*gL2+L2*gL3)*L1-4._DP*L3**2*gL2)+L3**3*gL2)*L1

		vec(116-k)=L3*gL1*(L2-L3)*L2**2+((2._DP*(gL1-gL2)*L3**2+((-2._DP*gL3+3._DP*gL2)*L3+L2*gL3)*L2)*L2+(L3**2*gL2+(-3._DP*gL1+2._DP*gL3)*L3*&
		L2+(-L3*gL2-L2*gL3)*L1)*L1)*L1

		vec(117-k)=L2**3*gL1*L3+(((3._DP*gL2-8._DP*gL1)*L3+L2*gL3)*L2**2+(((3._DP*gL1-8._DP*gL2)*L3-4._DP*L2*gL3)*L2+(L3*gL2+L2*gL3)*L1)*L1)*L1

		vec(118-k)=-L3*gL1*(L2-L3)*(L3**2+(-8._DP*L3+L2)*L2)*L2+(L3**4*gL2+(-4._DP*L3**3*gL3+((-9._DP*gL2+27._DP*gL3)*L3**2+((-18._DP*gL3+2._DP&
		*gL2)*L3+L2*gL3)*L2)*L2)*L2+(-9._DP*L3**3*gL2+(L2*(3._DP*L3*(-13._DP*gL3-gl4)+3._DP*L2*gL3)+(-9._DP*gL1+27._DP*gL3)*L3**2)*L2+(9._DP*L3&
		**2*gL2+((-18._DP*gL3+2._DP*gL1)*L3+3._DP*L2*gL3)*L2+(L2*gL3-L3*gL2)*L1)*L1)*L1)*L1

		vec(119-k)=-L3*gL1*(L3**2+(-4._DP*L3+L2)*L2)*L2**2+((2._DP*(gL1-gL2)*L3**3+L2*(L3**2*(-gL3-12._DP*gL1-4._DP*gl4)+((-8._DP*gL3+2._DP*gL1&
		)*L3+L2*gL3)*L2))*L2+(L3**3*gL2+(L2*(3._DP*(gL1-gL2)*L3+L2*gL3)+L3**2*(gL3+12._DP*gL2+4._DP*gl4))*L2+(-L2*(L2-8._DP*L3+L1)*gL3+L3*(-4._&
		DP*L3-2._DP*L2+L1)*gL2)*L1)*L1)*L1

		vec(120-k)=-L3*gL1*(L2-L3)*L2**3+((L2*(2._DP*L3*(5._DP*gL1+gl4)+L2*gL3)+(3._DP*gL2-8._DP*gL1)*L3**2)*L2**2+(L2*((3._DP*gL1-8._DP*gL2)*L&
		3**2+L2*(L3*(7._DP*gL3-gl4)-3._DP*L2*gL3))+L1*(L2*(L3*(10._DP*gL2+2._DP*gl4)-3._DP*L2*gL3)+(L2*gL3-L3*gL2)*L1+L3**2*gL2))*L1)*L1

		vec(121-k)=-L3*L2**4*gL1+(((18._DP*gL1-4._DP*gL2)*L3+L2*gL3)*L2**3+((-27._DP*(gL1-gL2)*L3-9._DP*L2*gL3)*L2**2+(((4._DP*gL1-18._DP*gL2)*&
		L3+9._DP*L2*gL3)*L2+(L3*gL2-L2*gL3)*L1)*L1)*L1)*L1

		vec(122-k)=(L1*gL2-L2*gL1)*L3*(-L3+L1+L2)*(1._DP+L4**2-2._DP*L4+10._DP*(-L2-L1)*L3)

		vec(123-k)=L1*(L3*(2._DP*L4**2-L4)*gL2+((2._DP*L4**2-L4)*gL3+L3*(L4*(4._DP*gl4+gL1)-gl4))*L2)+gL1*L2*L4*L3*(2._DP*L4-1._DP-L1)

		vec(124-k)=-gL1*L4*L3*(L2-L3)*L2+(L2*(L3*(4._DP*L4*gL3+(2._DP*L4+L3)*gl4)+(-L4*gL3-L3*gl4)*L2)+(-gL2*L4*L3+(-L4*gL3-L3*gl4)*L2)*L1+gL2*&
		L4*L3**2)*L1

		vec(125-k)=-gL1*L4*L3*L2**2+((2._DP*L3*(gL1-gL2)*L4+(-L4*gL3-L3*gl4)*L2)*L2+(gL2*L4*L3+(L4*gL3+L3*gl4)*L2)*L1)*L1

		vec(126-k)=L4*L2*(L3*(-2._DP*gL2+(-14._DP+12._DP*L3)*gl4)-gL3)+L4**3*L2*(L3*(-12._DP*gL2-24._DP*gl4)-12._DP*gL3)+L4**2*L2**2*(-18._DP*L&
		3*gl4-6._DP*gL3)+L4**2*L2*(L3*(12._DP*gL2+(36._DP-18._DP*L3)*gl4)+7._DP*gL3)+L4*L2**2*(12._DP*L3*gl4+gL3)+6._DP*L2*L4**3*(L4+L2)*gL3+L3&
		*L2*(L1+L4)*gl4-6._DP*gL2*L4**4*L3+L4*(-L3**2+L3)*gL2+L4**2*L3*(6._DP*L3-7._DP)*gL2-6._DP*L4**3*L3*(-2._DP+L3)*gL2

		vec(127-k)=L4*(6._DP*L4**2-6._DP*L4+1._DP)*gL2*(L2+L1)*L3+(L4*(L4-1._DP)*(6._DP*L4**2-6._DP*L4+1._DP)*gL3+(-12._DP*(gL2-gl4)*L4**3+L4**&
		2*(-24._DP*gl4+12._DP*gL2)+L4*(-2._DP*gL2+12._DP*gl4)-gl4+L3*gl4*(-12._DP*L4+18._DP*L4**2+1._DP))*L3+(L3*gl4*(-12._DP*L4+18._DP*L4**2+1&
		._DP)+L4*(6._DP*L4**2-6._DP*L4+1._DP)*gL3)*L2)*L2

		vec(128-k)=-L4*L1*(-L3**2+(-2._DP*L3+L2)*L2+(L3+L2)*L1)*(-2._DP*L4+1._DP)*gL3+L3*L1*(-2._DP*L4**3+3._DP*L4**2-L4+(-4._DP*L4**2+2._DP*L4&
		)*L3+(2._DP*L4**2-4._DP*L4+1._DP+(8._DP*L4-2._DP)*L3)*L2)*gl4+L4*L3*(L1-L2)*(-L3+L1+L2)*(2._DP*L4-1._DP)*gL1

		vec(129-k)=(L3*L2*(-L3+L1+L2)*(3._DP*L4-1._DP)*gl4-L4*L2*(2._DP*L4**2-3._DP*L4+1._DP+(5._DP*L4-2._DP-2._DP*L3)*L3)*gL3+L4*L3*((L4-L3)*L&
		3-L2**2+(L1-L4)*L1)*gL2)*L1-L4*L3*L2*((L3-L4)*L3+(L4-L2)*L2+L1**2)*gL1

		vec(130-k)=gL1*L4*L3*(L3**2+(-4._DP*L3+L2)*L2)*L2+(L3**3*gL2*L4+((-3._DP*L4*gL3+L3*gl4)*L3**2+L2*(L3*(-4._DP*L3*gl4+(8._DP*gL3-gL2)*L4)&
		+(L3*gl4-L4*gL3)*L2))*L2+(-4._DP*gL2*L4*L3**2+(L3*(-4._DP*L3*gl4+(8._DP*gL3-gL1)*L4)+2._DP*(L3*gl4-L4*gL3)*L2)*L2+(gL2*L4*L3+(L3*gl4-L4&
		*gL3)*L2)*L1)*L1)*L1

		vec(131-k)=gL1*L4*L3*(L3**2+(-4._DP*L3+L2)*L2)*L2+(L3**3*gL2*L4+L2*((-3._DP*L4*gL3-L3*gl4)*L3**2+L2*(L3*(4._DP*L3*gl4+(8._DP*gL3-gL2)*L&
		4)+(-L4*gL3-L3*gl4)*L2))+(-4._DP*gL2*L4*L3**2+(L3*(4._DP*L3*gl4+(8._DP*gL3-gL1)*L4)-2._DP*(L4*gL3+L3*gl4)*L2)*L2+(gL2*L4*L3+(-L4*gL3-L3&
		*gl4)*L2)*L1)*L1)*L1

		vec(132-k)=gL1*L4*L3*(-L4+L2+L3)*L2**2+((-2._DP*L3*(gL1-gL2)*L4*(L3-L4)+L2*(L3*(L3*gl4+(-2._DP*gl4+3._DP*gL2)*L4)+(L3*gl4-L4*gL3)*L2+L4&
		**2*gL3))*L2+(-L4*L3*(L3-L4)*gL2+(L3*(-L3*gl4+(-3._DP*gL1+2._DP*gl4)*L4)-L4**2*gL3)*L2+(-gL2*L4*L3+(L4*gL3-L3*gl4)*L2)*L1)*L1)*L1

		vec(133-k)=gL1*L4*L3*(-L4+L2+L3)*L2**2+L1*(L2*(2._DP*L3*L4*(gL2-gL1)*(L3-L4)+L2*(L3*(L4*(3._DP*gL2+4._DP*gl4)-L3*gl4)+(-L4*gL3-L3*gl4)*&
		L2+L4**2*gL3))+L1*(L2*(L3*(-L4*(3._DP*gL1+4._DP*gl4)+L3*gl4)-L4**2*gL3)+(-gL2*L4*L3+(L4*gL3+L3*gl4)*L2)*L1-L4*L3*(L3-L4)*gL2))

		vec(134-k)=gL1*L4*L3*(L2-L3)*L2**2+((2._DP*L3**2*(gL1-gL2)*L4+L2*(L3*(L4*(-3._DP*gL1+gL3-gl4)-L3*gl4)+(L3*gl4-L4*gL3)*L2))*L2+(L3*L2*(L&
		3*gl4+(3._DP*gL2-gL3+gl4)*L4)+(-gL2*L4*L3+(L4*gL3-L3*gl4)*L2)*L1+gL2*L4*L3**2)*L1)*L1

		vec(135-k)=gL1*L4*L3*(L2-L3)*L2**2+((2._DP*L3**2*(gL1-gL2)*L4+L2*(L3*(L4*(-3._DP*gL1+gL3-gl4)+L3*gl4)+(-L4*gL3-L3*gl4)*L2))*L2+(L3*L2*(&
		-L3*gl4+(3._DP*gL2-gL3+gl4)*L4)+(-gL2*L4*L3+(L4*gL3+L3*gl4)*L2)*L1+gL2*L4*L3**2)*L1)*L1

		vec(136-k)=L4*L3*L2**3*gL1+(L2**2*(L3*(3._DP*gL2-8._DP*gL1)*L4+(L3*gl4-L4*gL3)*L2)+((L3*(3._DP*gL1-8._DP*gL2)*L4-4._DP*(L3*gl4-L4*gL3)*&
		L2)*L2+(gL2*L4*L3+(L3*gl4-L4*gL3)*L2)*L1)*L1)*L1

		vec(137-k)=L4*L3*L2**3*gL1+(L2**2*(L3*(3._DP*gL2-8._DP*gL1)*L4+(-L4*gL3-L3*gl4)*L2)+((L3*(3._DP*gL1-8._DP*gL2)*L4+4._DP*(L4*gL3+L3*gl4)&
		*L2)*L2+(gL2*L4*L3+(-L4*gL3-L3*gl4)*L2)*L1)*L1)*L1

		vec(138-k)=-(L1*gL2-L2*gL1)*L3*(6._DP*L4**2-6._DP*L4+1._DP)*L4

		vec(139-k)=(L1*gL2-L2*gL1)*L3*(-L3+L1+L2)*(2._DP*L4-1._DP)*L4

		vec(140-k)=-(L1*gL2-L2*gL1)*L3*(1._DP+L4**2-2._DP*L4+6._DP*(-L2-L1)*L3)*L4

        end select poly
!
        p = p + 1
      end do
!
      end subroutine shapev3D
!
!
!
      subroutine curlshapev3D( l, gl, polylo, polyhi, vec, axs )
      Use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi, axs
      real (DP)     :: l(4), gl(3,4)
      real (DP)     :: vec(:)
      intent (in)   :: l, gl, polylo, polyhi, axs
      intent (out)  :: vec
!  
!  evaluate the shape functions for the corresponding polynomial degree
!  
!  Input:
!            l        coordinates to be calculated
!            gl       gradient of volume coordinates where:
!                        row indices denote  x,y,z components
!                        column indices denote lambdas (1,2,3,4) 
!            polylo   lowest polynomial degree to consider
!            polyhi   highest polynomial degree
!            axs      requested component of the curl, 1, 2, 3 for x,y,z  
!  Output:
!            vec      array of either x,y, or z component of CURL of shape functions 
!                     at coordinates of a point of interest 
!                     
!
!  The polylo and polyhi variables are used to select the shapefunctions. If polylo .eq.
!  polyhi, only the new shape functions of the corresponding polynomial degree p are cal-
!  culated. That means: if polylo .eq. polyhi = 7 ==> xsi(29) to xsi(36) are calculated.
!  If polylo .neq. polyhi, shape functions starting from polylo to polyhi are calculated.
!  That means: if polylo = 2 and polyhi = 7 ==> xsi(4) up to xsi(36) are calculated.
!
!  local variables:
      integer (I4B) p, k
      real (DP) L1, L2, L3, L4
      real (DP) gL11, gL12, gL13, & 
&               gL21, gL22, gL23, &
&               gL31, gL32, gL33, &
&               gL41, gL42, gL43
      real (DP) chi12, chi13, chi14, chi23, chi24, chi34
!
!  intermediate variables to avoid l(3) and gl(3) arrays
      L1 = l(1)
      L2 = l(2)
      L3 = l(3)
      L4 = l(4)
!
! intermediate variables to avoid gl(3,4) arrays
! use explicit substitutions for gL: 
! gLmn = gL(m,n), m = Volume Coord, n = Cartesian: (x,y,z) = (1,2,3)
!  
  gL11 = gl(1,1)
  gL12 = gl(2,1)
  gL13 = gl(3,1)

  gL21 = gl(1,2)
  gL22 = gl(2,2)
  gL23 = gl(3,2)

  gL31 = gl(1,3)
  gL32 = gl(2,3)
  gL33 = gl(3,3)

  gL41 = gl(1,4)
  gL42 = gl(2,4)
  gL43 = gl(3,4)

!
! substitute chi: chi_ij = gLi x gLj
! chi in reality is a vector according to the right hand rule
! components of chi_ij:
!      chi_ij(x) = gLi(y)*gLj(z) - gLi(z)*gLj(y)
!      chi_ij(y) = gLi(z)*gLj(x) - gLi(x)*gLj(z)
!      chi_ij(z) = gLi(x)*gLj(y) - gLi(y)*gLj(x)
!
! however, here chi_ij is a scalar variable.
! substitution of chi_ij depends on the requested component given by axs
!
! TO DO: need to do either some permutations or use chi as vectors, 
! or it is very long

! NOTES: 1. The chi and curl implemented here are checked in maple for p=1.
!        2. The the stiffness matrix is rank-deficient:
!             Rank = 3 for p=1 and Rank = 11 for p=2. Is it serious?? 

!        To verify for p=2
      cartesian: select case (axs)
         case(1)  ! x axis
            chi12 = gL12*gL23 - gL13*gL22
            chi13 = gL12*gL33 - gL13*gL32
            chi14 = gL12*gL43 - gL13*gL42
            chi23 = gL22*gL33 - gL23*gL32
            chi24 = gL22*gL43 - gL23*gL42
            chi34 = gL32*gL43 - gL33*gL42
         case(2)  ! y axis
            chi12 = gL13*gL21 - gL11*gL23
            chi13 = gL13*gL31 - gL11*gL33
            chi14 = gL13*gL41 - gL11*gL43
            chi23 = gL23*gL31 - gL21*gL33
            chi24 = gL23*gL41 - gL21*gL43
            chi34 = gL33*gL41 - gL31*gL43
         case(3)  ! z axis
            chi12 = gL11*gL22 - gL12*gL21
            chi13 = gL11*gL32 - gL12*gL31
            chi14 = gL11*gL42 - gL12*gL41
            chi23 = gL21*gL32 - gL22*gL31
            chi24 = gL21*gL42 - gL22*gL41
            chi34 = gL31*gL42 - gL32*gL41
      end select cartesian

      p = polylo
! k depends on the formulations of shape functions used for the 3D
! vector bases.

! for now k is for Nedelec
      if (polylo .gt. 1) then
        k = ( p-1 )*( p+1 )*( p+2 )/2
      else
        k = 0
      end if
!
!-------------------------------------------------------------------------------!
!  BEGIN LOOP TO CALCULATE CURL OF SHAPE FUNCTIONS FOR POLYNOMIAL DEGREE p  !
!-------------------------------------------------------------------------------!
      do
        if (p .gt. polyhi) exit    ! stopping criteria
!  Begin "case" control structure to select the appropriate shape functions for p
        poly: select case (p)
	case (1)  ! curl of shape functions for p = 1
		vec(1-k)=2._DP*chi12

		vec(2-k)=2._DP*chi23

		vec(3-k)=2._DP*chi13

		vec(4-k)=2._DP*chi14

		vec(5-k)=2._DP*chi24

		vec(6-k)=2._DP*chi34

	case (2)   ! curl of shape functions for p = 2
		vec(7-k)=0._DP

		vec(8-k)=0._DP

		vec(9-k)=0._DP

		vec(10-k)=0._DP

		vec(11-k)=0._DP

		vec(12-k)=0._DP

		vec(13-k)=-2._DP*L2*chi34-2._DP*L3*chi24

		vec(14-k)=L2*chi34-L3*chi24-2._DP*L4*chi23

		vec(15-k)=-2._DP*L1*chi34-2._DP*L3*chi14

		vec(16-k)=L1*chi34-L3*chi14-2._DP*L4*chi13

		vec(17-k)=-2._DP*L1*chi24-2._DP*L2*chi14

		vec(18-k)=L1*chi24-L2*chi14-2._DP*L4*chi12

		vec(19-k)=-2._DP*L1*chi23-2._DP*L2*chi13

		vec(20-k)=L1*chi23-L2*chi13-2._DP*L3*chi12

	case (3)   ! curl of shape functions for p = 3
		vec(21-k)=0._DP

		vec(22-k)=0._DP

		vec(23-k)=0._DP

		vec(24-k)=0._DP

		vec(25-k)=0._DP

		vec(26-k)=0._DP

		vec(27-k)=0._DP

		vec(28-k)=L3*(-2._DP*L2*chi14+(2._DP*chi23-4._DP*chi24)*L4)+2._DP*L2*(L2*chi34+(-chi23-2._DP*chi34)*L4)+2._DP*L3**2*chi24

		vec(29-k)=2._DP*L3**2*chi24+(-4._DP*(chi24-chi34)*L3-2._DP*L2*chi34)*L2

		vec(30-k)=-2._DP*L4**2*chi23+((3._DP*chi23-2._DP*chi24)*L3+2._DP*L2*(0.15000000000000000000e1_DP*chi23+chi34))*L4+(L3+L2)*(-L2*chi34+L3&
		*chi24)

		vec(31-k)=0._DP

		vec(32-k)=L3*(-2._DP*L1*chi24+(2._DP*chi13-4._DP*chi14)*L4)+2._DP*L1*(L1*chi34+(-chi13-2._DP*chi34)*L4)+2._DP*L3**2*chi14

		vec(33-k)=2._DP*L3**2*chi14+(-4._DP*(chi14-chi34)*L3-2._DP*L1*chi34)*L1

		vec(34-k)=-2._DP*L4**2*chi13+((3._DP*chi13-2._DP*chi14)*L3+2._DP*(0.15000000000000000000e1_DP*chi13+chi34)*L1)*L4+(L1+L3)*(-L1*chi34+L3&
		*chi14)

		vec(35-k)=0._DP

		vec(36-k)=L2*(-4._DP*L4*(-0.50000000000000000000_DP*chi12+chi14)-2._DP*L1*chi34)+2._DP*L2**2*chi14+2._DP*L1*(L1*chi24+(-chi12-2._DP*chi&
		24)*L4)

		vec(37-k)=2._DP*L2**2*chi14+(-4._DP*(chi14-chi24)*L2-2._DP*L1*chi24)*L1

		vec(38-k)=-2._DP*L4**2*chi12+((3._DP*chi12-2._DP*chi14)*L2+2._DP*(0.15000000000000000000e1_DP*chi12+chi24)*L1)*L4+(L2+L1)*(-L1*chi24+L2&
		*chi14)

		vec(39-k)=0._DP

		vec(40-k)=L2*(-4._DP*L3*(-0.50000000000000000000_DP*chi12+chi13)+2._DP*L1*chi34)+2._DP*L2**2*chi13+2._DP*L1*(L1*chi23+(-chi12-2._DP*chi&
		23)*L3)

		vec(41-k)=2._DP*L2**2*chi13+(-4._DP*(chi13-chi23)*L2-2._DP*L1*chi23)*L1

		vec(42-k)=-2._DP*L3**2*chi12+((3._DP*chi12-2._DP*chi13)*L2+2._DP*(0.15000000000000000000e1_DP*chi12+chi23)*L1)*L3+(L2+L1)*(-L1*chi23+L2&
		*chi13)

		vec(43-k)=-2._DP*(L1*chi23+L2*chi13)*L4+2._DP*L1*L2*chi34

		vec(44-k)=-2._DP*(L1*chi23+L2*chi13)*L4-2._DP*L3*(L1*chi24+L2*chi14)

		vec(45-k)=(L1*chi23-L2*chi13-2._DP*L3*chi12)*L4+L3*(L1*chi24-L2*chi14)

	case (4)   ! curl of shape functions for p = 4
		vec(46-k)=0._DP

		vec(47-k)=0._DP

		vec(48-k)=0._DP

		vec(49-k)=0._DP

		vec(50-k)=0._DP

		vec(51-k)=0._DP

		vec(52-k)=0._DP

		vec(53-k)=0._DP

		vec(54-k)=L3*(L3*(-2._DP*L3*chi24+(-4._DP*chi23+0.16e2_DP*chi24)*L4)+(8._DP*chi23-6._DP*chi24)*L4**2)+L2*(L3*(2._DP*(chi14-chi24)*L3-0.&
		16e2_DP*chi14*L4)+L2*(2._DP*(chi14-chi34)*L3-2._DP*L2*chi34+(0.16e2_DP*chi34+4._DP*chi23)*L4)+(-6._DP*chi34-8._DP*chi23)*L4**2)

		vec(55-k)=((-2._DP*chi23+4._DP*chi24)*L4-2._DP*L3*chi24)*L3**2+((8._DP*L4*(-chi24+chi23+chi34)+(2._DP*chi24-4._DP*chi34)*L3)*L3+(2._DP*&
		(-chi23-2._DP*chi34)*L4+(-2._DP*chi34+4._DP*chi24)*L3+2._DP*L2*chi34)*L2)*L2

		vec(56-k)=-2._DP*L3**3*chi24+(0.16e2_DP*(chi24-0.37500000000000000000_DP*chi34)*L3**2+(-6._DP*(chi24-0.26666666666666666667e1_DP*chi34)&
		*L3-2._DP*L2*chi34)*L2)*L2

		vec(57-k)=-0.12e2_DP*L4**3*chi23+(-0.12e2_DP*L1*chi23+(-0.18e2_DP*chi24+6._DP*chi12)*L3+(0.18e2_DP*chi34-6._DP*chi13)*L2+0.12e2_DP*chi2&
		3)*L4**2-2._DP*(L1*chi23+(6._DP*chi24-chi12)*L3+(chi13-6._DP*chi34)*L2-chi23)*(L1-1._DP)*L4+(L1-1._DP)**2*(L2*chi34-L3*chi24)

		vec(58-k)=0._DP

		vec(59-k)=0._DP

		vec(60-k)=-L3**3*chi14+(L1*chi34-4._DP*L4*(-4._DP*chi14+chi13))*L3**2+(-2._DP*L1*chi24*(9._DP*L4-1._DP+L2)+(L1**2-7._DP*L4**2+(2._DP-2.&
		_DP*L2)*L4-(-1._DP+L2)**2)*chi14+8._DP*L4**2*chi13)*L3-((6._DP*L4**2+L3**2+2._DP*(-8._DP*L4+L3+L1)*L1)*chi34-4._DP*L4*(L1-2._DP*L4)*chi&
		13)*L1

		vec(61-k)=((-2._DP*chi13+4._DP*chi14)*L4-2._DP*L3*chi14)*L3**2+((8._DP*L4*(chi13+chi34-chi14)+(2._DP*chi14-4._DP*chi34)*L3)*L3+(2._DP*(&
		-chi13-2._DP*chi34)*L4+(-2._DP*chi34+4._DP*chi14)*L3+2._DP*L1*chi34)*L1)*L1

		vec(62-k)=-2._DP*L3**3*chi14+(0.16e2_DP*(chi14-0.37500000000000000000_DP*chi34)*L3**2+(-6._DP*(chi14-0.26666666666666666667e1_DP*chi34)&
		*L3-2._DP*L1*chi34)*L1)*L1

		vec(63-k)=-0.12e2_DP*L4**3*chi13+(-0.12e2_DP*L2*chi13+(-0.18e2_DP*chi14-6._DP*chi12)*L3+(-6._DP*chi23+0.18e2_DP*chi34)*L1+0.12e2_DP*chi&
		13)*L4**2-2._DP*(L2*chi13+(chi12+6._DP*chi14)*L3+(-6._DP*chi34+chi23)*L1-chi13)*(-1._DP+L2)*L4+(-1._DP+L2)**2*(L1*chi34-L3*chi14)

		vec(64-k)=0._DP

		vec(65-k)=0._DP

		vec(66-k)=-4._DP*L2*(-2._DP*(-0.75000000000000000000_DP*chi14+chi12)*L4**2+(L4*(chi12-4._DP*chi14)+0.50000000000000000000_DP*L2*chi14)*&
		L2)+(-8._DP*(0.75000000000000000000_DP*chi24+chi12)*L4**2+(-2._DP*(chi14-chi34)*L2-0.16e2_DP*L4*chi34)*L2+(4._DP*L4*chi12+(0.18e2_DP*L4&
		+2._DP*L3-2._DP)*chi24+2._DP*L2*chi34)*L1)*L1

		vec(67-k)=-2._DP*L2**2*(L2*chi14+L4*(chi12-2._DP*chi14))+(L2*(chi14-L3*chi14+(chi14-4._DP*chi24)*L2+(-9._DP*chi14+8._DP*chi12+8._DP*chi&
		24)*L4)+L1*((-4._DP*chi24+3._DP*chi14)*L2+(-6._DP*chi24-2._DP*chi12)*L4-2._DP*chi24*(-1._DP+L3)))*L1

		vec(68-k)=-2._DP*L2**3*chi14+(0.16e2_DP*(chi14-0.37500000000000000000_DP*chi24)*L2**2+(-6._DP*(chi14-0.26666666666666666667e1_DP*chi24)&
		*L2-2._DP*L1*chi24)*L1)*L1

		vec(69-k)=-0.12e2_DP*L4**3*chi12+(-0.12e2_DP*L3*chi12+(-6._DP*chi13-0.18e2_DP*chi14)*L2+(6._DP*chi23+0.18e2_DP*chi24)*L1+0.12e2_DP*chi1&
		2)*L4**2-2._DP*(-1._DP+L3)*(L3*chi12+(chi13+6._DP*chi14)*L2+(-6._DP*chi24-chi23)*L1-chi12)*L4+(-1._DP+L3)**2*(L1*chi24-L2*chi14)

		vec(70-k)=0._DP

		vec(71-k)=0._DP

		vec(72-k)=-4._DP*L2*(-2._DP*L3**2*(-0.75000000000000000000_DP*chi13+chi12)+(L3*(chi12-4._DP*chi13)+0.50000000000000000000_DP*L2*chi13)*&
		L2)+(-8._DP*(0.75000000000000000000_DP*chi23+chi12)*L3**2+(-2._DP*(chi34+chi13)*L2+0.16e2_DP*L3*chi34)*L2+((4._DP*chi12+0.16e2_DP*chi23&
		)*L3+(-4._DP*chi23-2._DP*chi13)*L2-2._DP*L1*chi23)*L1)*L1

		vec(73-k)=((4._DP*chi13-2._DP*chi12)*L3-2._DP*L2*chi13)*L2**2+((8._DP*(chi23+chi12-chi13)*L3+(2._DP*chi13-4._DP*chi23)*L2)*L2+(2._DP*(-&
		chi12-2._DP*chi23)*L3+(-2._DP*chi23+4._DP*chi13)*L2+2._DP*L1*chi23)*L1)*L1

		vec(74-k)=-2._DP*L2**3*chi13+(0.16e2_DP*(chi13-0.37500000000000000000_DP*chi23)*L2**2+(-6._DP*(chi13-0.26666666666666666667e1_DP*chi23)&
		*L2-2._DP*L1*chi23)*L1)*L1

		vec(75-k)=0.18e2_DP*chi12*(L2+L1)*L3**2+(-2._DP*chi12*(L4-1._DP)**2+2._DP*(6._DP*L2*chi13+(1._DP-L4)*chi14)*L2+(0.12e2_DP*(chi13-chi23)&
		*L2+2._DP*chi24*(L4-1._DP)-0.12e2_DP*L1*chi23)*L1)*L3+(L4-1._DP)**2*(L1*chi23-L2*chi13)

		vec(76-k)=0._DP

		vec(77-k)=-4._DP*(L1*chi23+L2*chi13)*L4**2+((8._DP*L1*chi34+2._DP*chi13)*L2+2._DP*L1*chi23)*L4-2._DP*L1*L2*chi34

		vec(78-k)=(-4._DP*L1*chi23-3._DP*L2*chi13+L3*chi12)*L4**2+((-1._DP+L3)*chi12*L3+(-8._DP*L3*chi14+(1._DP+L2)*chi13)*L2+((-8._DP*chi24+ch&
		i12)*L3+L2*chi13+2._DP*chi23)*L1)*L4+(2._DP*L1*chi24-L2*(-2._DP+L4)*chi14)*L3

		vec(79-k)=L2*(-2._DP*L1*(-3._DP*L3+L1-L4)*chi34+(2._DP*chi12-4._DP*chi13)*L4*L3)-2._DP*(L1*chi34-L4*chi13)*L2**2-2._DP*(-L1*chi23+L3*(c&
		hi12+2._DP*chi23))*L4*L1

		vec(80-k)=-2._DP*(L1*chi24+L2*chi14)*L3**2+L3*(L2*(-2._DP*L1*chi34+2._DP*(chi12-2._DP*chi13)*L4)-2._DP*L2**2*(chi24+chi34)-2._DP*L1*(-L&
		1*chi24+L4*(chi12+2._DP*chi23)))+2._DP*L4*(L2+L1)*(L1*chi23+L2*chi13)

		vec(81-k)=-2._DP*(L1*chi34-L4*chi13)*L2**2+4._DP*L1*(0.50000000000000000000_DP*L1*chi34+L4*(-chi13+chi23))*L2-2._DP*L1**2*L4*chi23

		vec(82-k)=2._DP*(L3*chi14+L4*chi13)*L2**2+4._DP*(L4*(-chi13+chi23)+L3*(-chi14+chi24))*L1*L2-2._DP*L1**2*(L3*chi24+L4*chi23)

		vec(83-k)=2._DP*(L1*chi23-L2*chi13-2._DP*L3*chi12)*L4**2+((4._DP*L1*chi24-4._DP*L2*chi14+2._DP*chi12)*L3-L1*chi23+L2*chi13)*L4+L3*(-L1*&
		chi24+L2*chi14)

		vec(84-k)=(L1*chi24-L2*chi14-2._DP*L4*chi12)*L3**2+(((3._DP*chi12-2._DP*chi13)*L2+2._DP*(0.15000000000000000000e1_DP*chi12+chi23)*L1)*L&
		4+(L2+L1)*(-L1*chi24+L2*chi14))*L3-L4*(L2+L1)*(L1*chi23-L2*chi13)

	case (5)   ! curl of shape functions for p = 5
		vec(85-k)=0._DP

		vec(86-k)=0._DP

		vec(87-k)=0._DP

		vec(88-k)=0._DP

		vec(89-k)=0._DP

		vec(90-k)=0._DP

		vec(91-k)=0._DP

		vec(92-k)=0._DP

		vec(93-k)=0._DP

		vec(94-k)=(0.18e2_DP*L4**3*(-0.44444444444444444444_DP*chi24+chi23)+((0.54e2_DP*chi24-0.36e2_DP*chi23)*L4**2+(6._DP*L4*(-6._DP*chi24+ch&
		i23)+2._DP*L3*chi24)*L3)*L3)*L3+(-0.18e2_DP*L4**3*(0.44444444444444444444_DP*chi34+chi23)+(-0.54e2_DP*L4**2*chi14+(6._DP*L4*(chi12-7._D&
		P*chi24+6._DP*chi14)+(2._DP*chi34+6._DP*chi24)*L3)*L3)*L3+(0.36e2_DP*L4**2*(chi23+0.15000000000000000000e1_DP*chi34)+(-0.36e2_DP*L4*(-c&
		hi14+chi34+0.16666666666666666667_DP*chi23)-6._DP*L3*chi14)*L3+((4._DP*chi34-2._DP*chi14)*L3-6._DP*L4*(6._DP*chi34+chi23)+2._DP*L2*chi3&
		4)*L2)*L2)*L2

		vec(95-k)=4._DP*L4*(L3+L2)*((-4._DP*chi24+chi23)*L3**2+((-8._DP*chi34-4._DP*chi23+8._DP*chi24)*L3+(chi23+4._DP*chi34)*L2)*L2)+L4**2*(L3&
		**2*(6._DP*chi24-8._DP*chi23)+((0.32e2_DP*chi23+0.12e2_DP*chi34-0.12e2_DP*chi24)*L3+(-6._DP*chi34-8._DP*chi23)*L2)*L2)-4._DP*((-0.50000&
		000000000000000_DP*L3**2+L3*L2)*chi24+L2*(-L3+0.50000000000000000000_DP*L2)*chi34)*(L3+L2)**2

		vec(96-k)=L3**3*(2._DP*L3*chi24+(2._DP*chi23-4._DP*chi24)*L4)+L2*(L3**2*((-0.44e2_DP*chi34-0.22e2_DP*chi23-0.32e2_DP*chi14)*L4+L3*(0.20&
		e2_DP*chi34+0.14e2_DP*chi14))+L2*(L3*(0.10e2_DP*L3*chi14+(0.12e2_DP*chi14+0.44e2_DP*chi34+0.22e2_DP*chi23)*L4)+L2*(2._DP*L2*chi34+(-0.2&
		0e2_DP*chi34-6._DP*chi14)*L3+(-2._DP*chi23-4._DP*chi34)*L4)))

		vec(97-k)=2._DP*L3**4*chi24+(-0.36e2_DP*(chi24-0.22222222222222222222_DP*chi34)*L3**3+(0.54e2_DP*(chi24-chi34)*L3**2+(-8._DP*(chi24-0.4&
		5000000000000000000e1_DP*chi34)*L3-2._DP*L2*chi34)*L2)*L2)*L2

		vec(98-k)=L4**3*((0.30e2_DP*chi34+0.60e2_DP*chi23)*L2+(0.60e2_DP*chi23-0.30e2_DP*chi24)*L3+0.20e2_DP*L4*chi23)+L4**2*(-2._DP*L4**2*chi2&
		3+((0.28e2_DP*chi24-0.36e2_DP*chi23)*L4+(-0.44e2_DP*chi23+0.48e2_DP*chi24)*L3)*L3+((-0.28e2_DP*chi34-0.36e2_DP*chi23)*L4+(0.48e2_DP*chi&
		24-0.88e2_DP*chi23-0.48e2_DP*chi34)*L3+(-0.44e2_DP*chi23-0.48e2_DP*chi34)*L2)*L2)-2._DP*(L1-1._DP)*L4*(L3*((-0.15000000000000000000e1_D&
		P*chi23+chi24)*(L1-1._DP)+L2*(9._DP*chi34-9._DP*chi24+2._DP*chi23))-L2*((0.15000000000000000000e1_DP*chi23+chi34)*(L1-1._DP)+L2*(-chi23&
		-9._DP*chi34))+L3**2*(chi23-9._DP*chi24))-0.20e2_DP*L4**4*chi23+(L1-1._DP)**2*(L3+L2)*(-L2*chi34+L3*chi24)

		vec(99-k)=0._DP

		vec(100-k)=0._DP

		vec(101-k)=0._DP

		vec(102-k)=(0.18e2_DP*(-0.44444444444444444444_DP*chi14+chi13)*L4**3+((0.54e2_DP*chi14-0.36e2_DP*chi13)*L4**2+(6._DP*L4*(chi13-6._DP*ch&
		i14)+2._DP*L3*chi14)*L3)*L3)*L3+(-0.18e2_DP*L4**3*(0.44444444444444444444_DP*chi34+chi13)+(-0.54e2_DP*L4**2*chi24+(-6._DP*L4*(chi12-6._&
		DP*chi24+7._DP*chi14)+(2._DP*chi34+6._DP*chi14)*L3)*L3)*L3+(0.36e2_DP*(0.15000000000000000000e1_DP*chi34+chi13)*L4**2+(-0.36e2_DP*L4*(-&
		chi24+chi34+0.16666666666666666667_DP*chi13)-6._DP*L3*chi24)*L3+((4._DP*chi34-2._DP*chi24)*L3-6._DP*L4*(6._DP*chi34+chi13)+2._DP*L1*chi&
		34)*L1)*L1)*L1

		vec(103-k)=4._DP*L4*(L1+L3)*((-4._DP*chi14+chi13)*L3**2+(-4._DP*(chi13-2._DP*chi14+2._DP*chi34)*L3+(4._DP*chi34+chi13)*L1)*L1)+L4**2*(L&
		3**2*(-8._DP*chi13+6._DP*chi14)+(0.48e2_DP*L3*(-0.25000000000000000000_DP*chi14+0.66666666666666666667_DP*chi13+0.25000000000000000000_&
		DP*chi34)+(-8._DP*chi13-6._DP*chi34)*L1)*L1)-4._DP*((-0.50000000000000000000_DP*L3**2+L3*L1)*chi14+L1*(-L3+0.50000000000000000000_DP*L1&
		)*chi34)*(L1+L3)**2

		vec(104-k)=L3**3*(2._DP*L3*chi14+(2._DP*chi13-4._DP*chi14)*L4)+L1*(L3**2*((6._DP*chi34-0.14e2_DP*chi14)*L3+(-0.44e2_DP*chi34-0.22e2_DP*&
		chi13-0.32e2_DP*chi24)*L4)+L1*(L3*(0.10e2_DP*L3*chi24+(0.32e2_DP*chi34-0.12e2_DP*chi14+0.22e2_DP*chi13)*L4)+L1*(2._DP*L1*chi34+(-0.14e2&
		_DP*chi34+6._DP*chi14)*L3+(-4._DP*chi34-2._DP*chi13)*L4)))

		vec(105-k)=2._DP*L3**4*chi14+(-0.36e2_DP*(chi14-0.22222222222222222222_DP*chi34)*L3**3+(0.54e2_DP*(chi14-chi34)*L3**2+(-8._DP*(chi14-0.&
		45000000000000000000e1_DP*chi34)*L3-2._DP*L1*chi34)*L1)*L1)*L1

		vec(106-k)=L4**3*(L1*(0.30e2_DP*chi34+0.40e2_DP*chi13)+L3*(0.40e2_DP*chi13-0.30e2_DP*chi14)-0.20e2_DP*chi13*(-1._DP+L2))+L4**2*(L3*(-0.&
		20e2_DP*L1*(chi13+chi34-chi14)+(2._DP-2._DP*L2)*(-0.16e2_DP*chi13+0.14e2_DP*chi14))+L1**2*(-0.20e2_DP*chi34-0.10e2_DP*chi13)+0.30e2_DP*&
		(-1._DP+L2)*L1*(0.93333333333333333333_DP*chi34+0.10666666666666666667e1_DP*chi13)+L3**2*(-0.10e2_DP*chi13+0.20e2_DP*chi14)-2._DP*chi13&
		*(-1._DP+L2)**2)+2._DP*(-1._DP+L2)*L4*(L3*((-chi14+0.15000000000000000000e1_DP*chi13)*(-1._DP+L2)+L1*(9._DP*chi14-9._DP*chi34-2._DP*chi&
		13))+L1*((0.15000000000000000000e1_DP*chi13+chi34)*(-1._DP+L2)+L1*(-9._DP*chi34-chi13))+L3**2*(-chi13+9._DP*chi14))-0.20e2_DP*L4**4*chi&
		13-(-1._DP+L2)**2*(L1+L3)*(L1*chi34-L3*chi14)

		vec(107-k)=0._DP

		vec(108-k)=0._DP

		vec(109-k)=0._DP

		vec(110-k)=2._DP*(9._DP*(-0.44444444444444444444_DP*chi14+chi12)*L4**3+(-0.18e2_DP*(chi12-0.15000000000000000000e1_DP*chi14)*L4**2+(3._&
		DP*L4*(chi12-6._DP*chi14)+L2*chi14)*L2)*L2)*L2+(-0.18e2_DP*(chi12+0.44444444444444444444_DP*chi24)*L4**3+(-0.54e2_DP*L4**2*chi34+(6._DP&
		*L4*(6._DP*chi24+chi12+0.12e2_DP*chi34)+(-4._DP*chi24-6._DP*chi34)*L2)*L2)*L2+(2._DP*((0.19e2_DP*L4-1._DP+L3)*chi12+0.27e2_DP*chi24*L4)&
		*L4+(-4._DP*L4*chi12+2._DP*(0.18e2_DP*L4-L1)*chi14+(0.72e2_DP*L4-6._DP*L2)*chi34)*L2+((6._DP*chi24+4._DP*chi14)*L2-4._DP*L4*(9._DP*chi2&
		4+chi12)+2._DP*L1*chi24)*L1)*L1)*L1

		vec(111-k)=L2**4*chi14+(-2._DP*chi14*L1+4._DP*L4*(chi12-2._DP*chi14))*L2**3+(L1**2*(-3._DP*chi34-6._DP*chi14)+L1*(L4*(-0.12e2_DP*chi12+&
		0.18e2_DP*chi34+0.24e2_DP*chi14)+2._DP*L3*chi34-2._DP*chi34)+(-8._DP*chi12+0.15e2_DP*chi14)*L4**2+0.10e2_DP*chi14*(-1._DP+L3)*L4+chi14*&
		(-1._DP+L3)**2)*L2**2+2._DP*L1*(L1**2*chi24+L1*(L4*(3._DP*chi34+0.12e2_DP*chi14-6._DP*chi12)-L3*chi34+chi34)+(-0.15e2_DP*chi14+0.16e2_D&
		P*chi12+0.15e2_DP*chi24)*L4**2-0.10e2_DP*(chi14-chi24)*(-1._DP+L3)*L4-(-1._DP+L3)**2*(chi14-chi24))*L2-L1**2*((-0.10e2_DP*L4+0.15e2_DP*&
		L4**2+1._DP+(0.10e2_DP*L4-2._DP+L3)*L3+(-8._DP*L4+L1)*L1)*chi24-4._DP*L4*(L1-2._DP*L4)*chi12)

		vec(112-k)=(4._DP*L4*(-chi14+0.50000000000000000000_DP*chi12)+2._DP*L2*chi14)*L2**3+(L2**2*((-6._DP*chi34-0.20e2_DP*chi14)*L2+(0.66e2_D&
		P*chi14+0.22e2_DP*chi13+0.12e2_DP*chi34)*L4)+((0.10e2_DP*L2*chi34+(-0.66e2_DP*chi14-0.22e2_DP*chi13-0.32e2_DP*chi34)*L4)*L2+(-2._DP*L4*&
		chi23+(-6._DP*L4+2._DP*L1-0.20e2_DP*L2)*chi24-6._DP*L2*chi34)*L1)*L1)*L1

		vec(113-k)=2._DP*L2**4*chi14+(-0.36e2_DP*(chi14-0.22222222222222222222_DP*chi24)*L2**3+(0.54e2_DP*(chi14-chi24)*L2**2+(-8._DP*(chi14-0.&
		45000000000000000000e1_DP*chi24)*L2-2._DP*L1*chi24)*L1)*L1)*L1

		vec(114-k)=L4**3*(L1*(0.30e2_DP*chi24+0.40e2_DP*chi12)+L2*(-0.30e2_DP*chi14+0.40e2_DP*chi12)-0.20e2_DP*chi12*(-1._DP+L3))+L4**2*(L2*(-0&
		.20e2_DP*L1*(chi24+chi12-chi14)+(-0.30e2_DP+0.30e2_DP*L3)*(-0.93333333333333333333_DP*chi14+0.10666666666666666667e1_DP*chi12))+L1**2*(&
		-0.20e2_DP*chi24-0.10e2_DP*chi12)+0.30e2_DP*(-1._DP+L3)*L1*(0.93333333333333333333_DP*chi24+0.10666666666666666667e1_DP*chi12)+L2**2*(0&
		.20e2_DP*chi14-0.10e2_DP*chi12)-2._DP*chi12*(-1._DP+L3)**2)-2._DP*(-1._DP+L3)*L4*(((-0.15000000000000000000e1_DP*chi12+chi14)*(-1._DP+L&
		3)+(-9._DP*chi14+chi12)*L2)*L2+(L2*(9._DP*chi24+2._DP*chi12-9._DP*chi14)+(-0.15000000000000000000e1_DP*chi12-chi24)*(-1._DP+L3)+L1*(9._&
		DP*chi24+chi12))*L1)-0.20e2_DP*L4**4*chi12-(-1._DP+L3)**2*(L2+L1)*(L1*chi24-L2*chi14)

		vec(115-k)=0._DP

		vec(116-k)=0._DP

		vec(117-k)=0._DP

		vec(118-k)=(0.18e2_DP*L3**3*(chi12-0.44444444444444444444_DP*chi13)+((-0.36e2_DP*chi12+0.54e2_DP*chi13)*L3**2+(6._DP*L3*(-6._DP*chi13+c&
		hi12)+2._DP*L2*chi13)*L2)*L2)*L2+(-0.18e2_DP*(0.44444444444444444444_DP*chi23+chi12)*L3**3+(0.54e2_DP*L3**2*chi34+(-6._DP*L3*(chi14+7._&
		DP*chi13+6._DP*chi34)+(2._DP*chi23+6._DP*chi13)*L2)*L2)*L2+(0.36e2_DP*(chi12+0.15000000000000000000e1_DP*chi23)*L3**2+(-0.36e2_DP*(chi1&
		3+0.16666666666666666667_DP*chi24+0.21666666666666666667e1_DP*chi23)*L3+6._DP*L2*chi34)*L2+((2._DP*chi13+6._DP*chi23)*L2-6._DP*L3*(6._D&
		P*chi23+chi12)+2._DP*L1*chi23)*L1)*L1)*L1

		vec(119-k)=4._DP*L3*(L2+L1)*((chi12-4._DP*chi13)*L2**2+L1*(L2*(8._DP*chi13-8._DP*chi23-4._DP*chi12)+(chi12+4._DP*chi23)*L1))+L3**2*(L2*&
		*2*(6._DP*chi13-8._DP*chi12)+(0.48e2_DP*L2*(-0.25000000000000000000_DP*chi13+0.66666666666666666667_DP*chi12+0.25000000000000000000_DP*&
		chi23)+(-8._DP*chi12-6._DP*chi23)*L1)*L1)-2._DP*(L2+L1)**2*(-L2**2*chi13+(2._DP*(chi13-chi23)*L2+L1*chi23)*L1)

		vec(120-k)=(4._DP*L3*(-chi13+0.50000000000000000000_DP*chi12)+2._DP*L2*chi13)*L2**3+(((6._DP*chi23-0.14e2_DP*chi13)*L2+(-0.44e2_DP*chi2&
		3+0.32e2_DP*chi34-0.22e2_DP*chi12)*L3)*L2**2+(L2*(-0.10e2_DP*L2*chi34+(0.22e2_DP*chi12+0.32e2_DP*chi23-0.12e2_DP*chi13)*L3)+(2._DP*(-ch&
		i12-2._DP*chi23)*L3+(-0.14e2_DP*chi23+6._DP*chi13)*L2+2._DP*L1*chi23)*L1)*L1)*L1

		vec(121-k)=2._DP*L2**4*chi13+(-0.36e2_DP*(chi13-0.22222222222222222222_DP*chi23)*L2**3+(0.54e2_DP*(chi13-chi23)*L2**2+(-8._DP*(chi13-0.&
		45000000000000000000e1_DP*chi23)*L2-2._DP*L1*chi23)*L1)*L1)*L1

		vec(122-k)=0.30e2_DP*chi12*(L2+L1)*L3**3+(-2._DP*chi12*(L4-1._DP)**2+(-0.70e2_DP*L2*chi12+(2._DP-2._DP*L4-0.30e2_DP*L2)*chi14)*L2+((-0.&
		80e2_DP*chi12-0.30e2_DP*chi23+0.30e2_DP*chi13)*L2+2._DP*chi24*(L4-1._DP)+(-0.40e2_DP*chi12-0.30e2_DP*chi23)*L1)*L1)*L3**2+((-2._DP*(chi&
		13-0.15000000000000000000e1_DP*chi12)*(L4-1._DP)**2-2._DP*(0.10e2_DP*L2*chi13+(1._DP-L4)*chi14)*L2)*L2+(2._DP*(L4-1._DP)**2*(0.15000000&
		000000000000e1_DP*chi12+chi23)+(2._DP*(chi14-chi24)*(L4-1._DP)+(-0.40e2_DP*chi13+0.20e2_DP*chi23)*L2)*L2+((0.40e2_DP*chi23-0.20e2_DP*ch&
		i13)*L2-2._DP*chi24*(L4-1._DP)+0.20e2_DP*L1*chi23)*L1)*L1)*L3-(L4-1._DP)**2*(L2+L1)*(L1*chi23-L2*chi13)

		vec(123-k)=0._DP

		vec(124-k)=0._DP

		vec(125-k)=0._DP

		vec(126-k)=chi34*L2*(2._DP*L1+(-2._DP-0.24e2_DP*L1)*L4+(0.12e2_DP+0.36e2_DP*L1)*L4**2-0.12e2_DP*L4**3)+2._DP*L4*(6._DP*L4**2-6._DP*L4+1&
		._DP)*(L2-L1)*chi23

		vec(127-k)=6._DP*L4**4*chi23+((6._DP*chi12+0.12e2_DP*chi24)*L3+(0.18e2_DP*chi23-0.12e2_DP*chi34)*L2-6._DP*chi23*(L1+2._DP))*L4**3+(0.18&
		e2_DP*L3**2*chi24+((0.54e2_DP*L2-0.24e2_DP-0.18e2_DP*L1)*chi24+0.36e2_DP*L2*chi34-6._DP*chi12)*L3+(-0.18e2_DP*chi23+0.12e2_DP*chi34)*L2&
		+chi23*(7._DP+6._DP*L1))*L4**2+(L2*(chi34-3._DP*chi13)+(L1+1._DP)*(-chi34+chi13)-0.12e2_DP*L3**2*chi24+((0.12e2_DP+0.12e2_DP*L1-0.36e2_&
		DP*L2)*chi24-0.24e2_DP*L2*chi34+chi12)*L3)*L4-L3*(-2._DP*L2*chi34+(L1-3._DP*L2-L3+1._DP)*chi24)

		vec(128-k)=L4**2*(L1*chi34*(-4._DP*L1-0.20e2_DP*L3-6._DP*L4+9._DP)+chi14*(2._DP*L4*L3+L3)+4._DP*chi13*(-1._DP+4._DP*L3+L4)*(L1-L2))+L4*&
		(((-1._DP+0.10e2_DP*L3+(0.24e2_DP*L3-6._DP-4._DP*L2)*L2+(-1._DP-4._DP*L3-4._DP*L2)*L1)*chi34+4._DP*L3**2*chi24)*L1-L3*(L4+(-1._DP+4._DP&
		*L4)*L2+(-3._DP+4._DP*L2+4._DP*L1)*L1)*chi14+2._DP*(L1-L2)*(L1-3._DP*L3+L2)*chi13)+L3**2*(-L1*chi24-L2*chi14)+L3*((-L3-L4)*(L1-L2)*chi1&
		4+L1*(L1-1._DP-6._DP*L2)*chi34)-2._DP*(L1*chi34+L3*chi14)*L4**3+L1*chi34*L2*(L2+L1+1._DP)

		vec(129-k)=L4*((-chi13+(2._DP*chi13+5._DP*(-chi14+chi13)*L3)*L3+(5._DP*L3*(-0.80000000000000000000_DP*chi12+chi14)-L2*chi13)*L2)*L2+(-c&
		hi23+(2._DP*chi23+5._DP*(-chi24+chi23)*L3)*L3+(-2._DP*chi34*(L3+3._DP)+(chi23+3._DP*chi34)*L2)*L2+(L2*(3._DP*chi34+chi13)+5._DP*(0.8000&
		0000000000000000_DP*chi12+chi24)*L3-L1*chi23)*L1)*L1)+((-chi34+L3*chi24)*L1-L3*chi14)*L2**2+L2*(L3**2*chi14*(L3+1._DP)+(-chi34*(-1._DP+&
		(-2._DP+2._DP*L3)*L3)+(L3*chi14-chi34)*L1)*L1)+L1*L3*(L1+1._DP+L3)*(L1-L3)*(chi34+chi14)-2._DP*(L1*chi23+L2*chi13)*L4**3+(L2**2*chi13+(&
		6._DP*L1*chi34+(2._DP*chi12-7._DP*chi13)*L3+3._DP*chi13)*L2+L1*(-2._DP*L3*chi12+(L1-7._DP*L3+3._DP)*chi23))*L4**2-L2**3*L3*chi14

		vec(130-k)=L2**2*(0.16e2_DP*(-0.25000000000000000000_DP*chi12+chi13)*L4*L3+(-2._DP*L4*chi13+(-0.20e2_DP*L3+4._DP*L1-2._DP*L4)*chi34)*L1&
		)+L2*(-6._DP*(chi13-0.13333333333333333333e1_DP*chi12)*L4*L3**2+(2._DP*L3*chi34*(8._DP*L4+7._DP*L3)+(2._DP*L4*chi13+(-4._DP*L4-0.20e2_D&
		P*L3+2._DP*L1)*chi34)*L1)*L1)+2._DP*(L1*chi34-L4*chi13)*L2**3-2._DP*L1*(3._DP*(chi23+0.13333333333333333333e1_DP*chi12)*L3**2+(-8._DP*(&
		0.25000000000000000000_DP*chi12+chi23)*L3+L1*chi23)*L1)*L4

		vec(131-k)=-2._DP*(L1*chi24+L2*chi14)*L3**3+((8._DP*L2*chi14+(-6._DP*chi13+8._DP*chi12)*L4)*L2+(-8._DP*L2*chi34-6._DP*(chi23+0.13333333&
		333333333333e1_DP*chi12)*L4+8._DP*L1*chi24)*L1)*L3**2-2._DP*(L2+L1)*((-8._DP*(-0.25000000000000000000_DP*chi12+chi13)*L4+L2*chi14)*L2+(&
		-L2*chi34-8._DP*L4*(0.25000000000000000000_DP*chi12+chi23)+L1*chi24)*L1)*L3-2._DP*L4*(L2+L1)**2*(L1*chi23+L2*chi13)

		vec(132-k)=2._DP*(L1*chi34-L4*chi13)*L2**3+(((-6._DP*chi23-4._DP*chi34)*L4+2._DP*L3*chi34)*L1-2._DP*L4*chi13*(L3-L4))*L2**2-2._DP*L1*(2&
		._DP*L4*(-chi13+chi23)*(L3-L4)+(-3._DP*L4*chi13+(-L2-3._DP*L4+1._DP)*chi34)*L1)*L2+2._DP*L1**2*L4*chi23*(-L4+L1+L3)

		vec(133-k)=-2._DP*(L3*chi14+L4*chi13)*L2**3+(2._DP*L4**2*chi13+(6._DP*(-0.33333333333333333333_DP*chi13+chi14)*L4-2._DP*L3*chi14)*L3+L1&
		*((-6._DP*chi24-2._DP*chi34)*L3+(-6._DP*chi23+2._DP*chi34)*L4))*L2**2+2._DP*L1*(((3._DP*chi13-chi34)*L4+L3*(3._DP*chi14+chi34))*L1-2._D&
		P*(chi13-chi23)*L4**2+6._DP*L3*(-chi14-0.33333333333333333333_DP*chi23+0.33333333333333333333_DP*chi13+chi24)*L4-2._DP*L3**2*(-chi14+ch&
		i24))*L2+2._DP*L1**2*(L4*(-L4+L1+L3)*chi23+L3*(L1-3._DP*L4+L3)*chi24)

		vec(134-k)=L2**2*(3._DP*L4*L3*(-0.33333333333333333333_DP*chi14-chi12+chi13)+L1*(L4*(-4._DP*chi34+6._DP*chi13)-6._DP*L3*chi34))-2._DP*L&
		1*L2*(4._DP*L4*L3*(3._DP*chi13+chi14-chi34)+(-3._DP*L4*chi13+(L4-3._DP*L3+L1)*chi34)*L1)+2._DP*L4*L1**2*((L1-3._DP*L3)*chi23-L3*chi24)+&
		2._DP*(L1*chi34-L4*chi13)*L2**3

		vec(135-k)=-2._DP*(L3*chi14+L4*chi13)*L2**3+(L1*(L3*(4._DP*chi34+6._DP*chi14)-L4*(4._DP*chi34-6._DP*chi13))+2._DP*(3._DP*L4*chi13+(L3+L&
		4)*chi14)*L3)*L2**2-L1*(-4._DP*L3*(2._DP*(-chi34+3._DP*chi23+chi24)*L4+L3*(-chi14+chi24))+L1*(L3*(6._DP*chi24+4._DP*chi34)-L4*(-6._DP*c&
		hi23+4._DP*chi34)))*L2+2._DP*L1**2*(L4*(L1-3._DP*L3)*chi23+L3*(L1-L3-L4)*chi24)

		vec(136-k)=-2._DP*L2**3*L4*chi13+((0.16e2_DP*(chi13-0.37500000000000000000_DP*chi23)*L4+2._DP*L2*chi34)*L2**2+((-6._DP*(chi13-0.2666666&
		6666666666667e1_DP*chi23)*L4-8._DP*L2*chi34)*L2+2._DP*(L2*chi34-L4*chi23)*L1)*L1)*L1

		vec(137-k)=-2._DP*(L3*chi14+L4*chi13)*L2**3+(0.16e2_DP*((chi13-0.37500000000000000000_DP*chi23)*L4+(chi14-0.37500000000000000000_DP*chi&
		24)*L3)*L2**2+(-6._DP*((chi13-0.26666666666666666667e1_DP*chi23)*L4+(chi14-0.26666666666666666667e1_DP*chi24)*L3)*L2-2._DP*(L3*chi24+L4&
		*chi23)*L1)*L1)*L1

		vec(138-k)=6._DP*(L1*chi23-L2*chi13-2._DP*L3*chi12)*L4**3+((-0.18e2_DP*L2*chi14+0.12e2_DP*chi12+0.18e2_DP*L1*chi24)*L3+6._DP*L2*chi13-6&
		._DP*L1*chi23)*L4**2+((-2._DP*chi12+0.12e2_DP*L2*chi14-0.12e2_DP*L1*chi24)*L3+L1*chi23-L2*chi13)*L4+L3*(L1*chi24-L2*chi14)

		vec(139-k)=(-4._DP*L3**2*chi12+((6._DP*chi12-4._DP*chi13)*L2+6._DP*L1*(0.66666666666666666667_DP*chi23+chi12))*L3-2._DP*(L2+L1)*(L1*chi&
		23-L2*chi13))*L4**2+((4._DP*L1*chi24-4._DP*L2*chi14+2._DP*chi12)*L3**2+(4._DP*L2**2*chi14+(4._DP*(chi14-chi24)*L1-3._DP*chi12+2._DP*chi&
		13)*L2-4._DP*L1*(0.12500000000000000000e1_DP*chi12+(-0.50000000000000000000_DP+L1)*chi24))*L3+(L2+L1)*(L1*chi23-L2*chi13))*L4+L3*(L1*ch&
		i24-L2*chi14)*(-L3+L1+L2)

		vec(140-k)=(L1*chi23-L2*chi13-2._DP*L3*chi12)*L4**3+((3._DP*L1*chi24+4._DP*chi12-3._DP*L2*chi14)*L3+2._DP*L2*chi13-2._DP*L1*chi23)*L4**&
		2+(0.18e2_DP*chi12*(L2+L1)*L3**2+(-2._DP*chi12+(4._DP*chi14+0.12e2_DP*L2*chi13)*L2+(0.12e2_DP*(chi13-chi23)*L2-4._DP*chi24-0.12e2_DP*L1&
		*chi23)*L1)*L3-L2*chi13+L1*chi23)*L4-6._DP*(-0.16666666666666666667_DP+(L2+L1)*L3)*L3*(L1*chi24-L2*chi14)

        end select poly
!
        p = p + 1
      end do
!
      end subroutine curlshapev3D
