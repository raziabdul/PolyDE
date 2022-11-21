      subroutine shapefunctionsc3D(lambda, vert, polylo, polyhi,   &
     &  nff, grad, xsi, gxsi, errcode)
      use feminterface3D, only: shapesc3D, gradshapesc3D, getgl
      use femtypes
      use globalvariables3D, only: polymaxsc
      implicit none
      integer (I4B) :: polylo, polyhi, nff, errcode
      real (DP) :: vert(3,4), lambda(4)
      real (DP) :: xsi(:), gxsi(:,:)
      logical grad
      intent (in) :: lambda, vert, polylo, polyhi, nff, grad
      intent (out) :: xsi, gxsi
      intent (inout) :: errcode
!
!    $Revision: 1.10 $
!    $Date: 2014/08/22 10:47:59 $
!    $Author: m_kasper $
!
!  evaluate the shape function and its grad at a point given in natural 
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
!            grad       if .true. the gradient will be computed                     
!  Output:                                                                           
!            xsi         vector of shapefunction at lambda                           
!            gxsi        gradient of the shapefunction at xi, yi, zi                     
!                        where gxsi(i,:) representing the i-component, i=x,y,z 
!  Errorcodes:                                                                       
!            errcode     2001 = polylo is greater than polyhi                        
!                        2002 = polynomial degree cannot be higher than polymax      
!                        4001 = unavailable element type (currently only Nedelec)
!
!  global variables used: 
!            polymaxsc
!         
!  local variables
      real(DP) :: gl(3,4)
!
!  Check values of polylo and polyhi
      if ( (polylo .gt. polymaxsc) .or. (polyhi .gt. polymaxsc) ) then
        errcode = 2002
        return
      end if
      if ( polylo .gt. polyhi ) then
        errcode = 2001
        return
      end if
!
      if (grad) then
!  Calculation of the gradients
        call getgl(vert, gl)
!  Calculate grad vector for each component. The last argument is 
!  to identify the requested component
!
! TO DO : make it call only once and return all components (may be nice but gives no speed-up)
!
! TO DO : gl(1,:), gl(2,:) creates a copy of a temporary array which may lead to lower performance
!         it would be better to have it as gl(:,1) (interchange indices such that memory ist continuous
!
! x-component
        call gradshapesc3D(lambda, gl(1,:), polylo, polyhi, gxsi(1,:))
! y-component
        call gradshapesc3D(lambda, gl(2,:), polylo, polyhi, gxsi(2,:))
! z-component
        call gradshapesc3D(lambda, gl(3,:), polylo, polyhi, gxsi(3,:))
      end if
!
!  Calculate the shape function vector
      call shapesc3D(lambda, polylo, polyhi, xsi(:))

      end subroutine shapefunctionsc3D
!
!
!
      subroutine shapesc3D(l, polylo, polyhi, vec)
      use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi
      real (DP)     :: l(4)
      real (DP)     :: vec(:)
      intent (in)   :: l, polylo, polyhi
      intent (out)  :: vec
!  
!  evaluate the shape functions for the corresponding polynomial degree
!  
!  Input:
!            l        coordinates to be calculated 
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
!
!  SIMPLE HIERARCHICAL FUNCTIONS
! 
!  The shape functions follow simple building rules:
!    vertex shape functions:   xsi = Li                                 first order only
!    edge shape functions:     xsi = Li*Lj 
!                                    * (Li-Lj)**u                       u=p-2
!    face shape functions:     xsi = Li*Lj*Lk 
!                                    * (Li-Lj)^u*(Lj-Lk)^v              u+v=p-3
!    volume shape functions:   xsi = Li*Lj*Lk*Ll 
!                                    * (Li-Lj)^u*(Lj-Lk)^v*(Lk-Ll)^w    u+v+w=p-4
!
!  Edges             1  |  2  |  3  |  4  |  5  |  6  |
!  indices(i,j)     1,2 | 2,3 | 1,3 | 1,4 | 2,4 | 3,4 |
!
!  Faces              1   |   2   |   3   |   4   |
!  indices(i,j,k)   2,3,4 | 1,3,4 | 1,2,4 | 1,2,3 | 
!--------------------------------------------------------------------------
!
!
!  local variables:
      integer (I4B) k, p
      real (DP) L1, L2, L3, L4
!
!  intermediate variables to avoid arrays of l and gl
      L1 = l(1)
      L2 = l(2)
      L3 = l(3)
      L4 = l(4)

      p = polylo

! k depends on the formulations of shape functions used.
! This is the standard C(0) type element (continuous)
      if (polylo .gt. 1) then
        k = ( p+1 )*( p+2 )*( p+3 )/2
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
        vec(  1-k) = L1
        vec(  2-k) = L2
        vec(  3-k) = L3
        vec(  4-k) = L4
    case (2)   !  shape functions for p = 2
!  edge DOF
        vec(  5-k) = L1*L2
        vec(  6-k) = L2*L3
        vec(  7-k) = L1*L3
        vec(  8-k) = L1*L4
        vec(  9-k) = L2*L4
        vec( 10-k) = L3*L4
    case (3)   !  shape functions for p = 3
!  edge DOF
        vec( 11-k) = L1 * L2 * (L1 - L2)
        vec( 12-k) = L2 * L3 * (L2 - L3)
        vec( 13-k) = L1 * L3 * (L1 - L3)
        vec( 14-k) = L1 * L4 * (L1 - L4)
        vec( 15-k) = L2 * L4 * (L2 - L4)
        vec( 16-k) = L3 * L4 * (L3 - L4)
!  face DOF
        vec( 17-k) = L2 * L4 * L3
        vec( 18-k) = L1 * L3 * L4
        vec( 19-k) = L1 * L4 * L2
        vec( 20-k) = L1 * L2 * L3
    case (4)   !  shape functions for p = 4
!  edge DOF
        vec( 21-k) = L1 * L2 * (L1 - L2)**2
        vec( 22-k) = L2 * L3 * (L2 - L3)**2
        vec( 23-k) = L1 * L3 * (L1 - L3)**2
        vec( 24-k) = L1 * L4 * (L1 - L4)**2
        vec( 25-k) = L2 * L4 * (L2 - L4)**2
        vec( 26-k) = L3 * L4 * (L3 - L4)**2
!  face DOF
        vec( 27-k) = L2 * L3 * L4 * (L2 - L3)
        vec( 28-k) = L2 * L3 * L4 * (L3 - L4)
        vec( 29-k) = L1 * L3 * L4 * (L1 - L3)
        vec( 30-k) = L1 * L3 * L4 * (L3 - L4)
        vec( 31-k) = L1 * L2 * L4 * (L1 - L2)
        vec( 32-k) = L1 * L2 * L4 * (L2 - L4)
        vec( 33-k) = L1 * L2 * L3 * (L1 - L2)
        vec( 34-k) = L1 * L2 * L3 * (L2 - L3)
!  volume DOF
        vec( 35-k) = L1 * L2 * L3 * L4

    case (5)   !  shape functions for p = 5
!  edge DOF
        vec( 36-k) = L1 * L2 * (L1 - L2)**3
        vec( 37-k) = L2 * L3 * (L2 - L3)**3
        vec( 38-k) = L1 * L3 * (L1 - L3)**3
        vec( 39-k) = L1 * L4 * (L1 - L4)**3
        vec( 40-k) = L2 * L4 * (L2 - L4)**3
        vec( 41-k) = L3 * L4 * (L3 - L4)**3
!  face DOF
        vec( 42-k) = L2 * L3 * L4 * (L2 - L3)**2
        vec( 45-k) = L1 * L3 * L4 * (L1 - L3)**2
        vec( 48-k) = L1 * L2 * L4 * (L1 - L2)**2
        vec( 51-k) = L1 * L2 * L3 * (L1 - L2)**2

        vec( 43-k) = L2 * L3 * L4 * (L2 - L3) * (L3 - L4)
        vec( 46-k) = L1 * L3 * L4 * (L1 - L3) * (L3 - L4)
        vec( 49-k) = L1 * L2 * L4 * (L1 - L2) * (L2 - L4)
        vec( 52-k) = L1 * L2 * L3 * (L1 - L2) * (L2 - L3)

        vec( 44-k) = L2 * L3 * L4 * (L3 - L4)**2
        vec( 47-k) = L1 * L3 * L4 * (L3 - L4)**2
        vec( 50-k) = L1 * L2 * L4 * (L2 - L4)**2
        vec( 53-k) = L1 * L2 * L3 * (L2 - L3)**2
!  volume DOF
        vec( 54-k) = L1 * L2 * L3 * L4 * (L1 - L2)
        vec( 55-k) = L1 * L2 * L3 * L4 * (L2 - L3)
        vec( 56-k) = L1 * L2 * L3 * L4 * (L3 - L4)

    case (6)   !  shape functions for p = 6
!  edge DOF
        vec( 57-k) = L1 * L2 * (L1 - L2)**4
        vec( 58-k) = L2 * L3 * (L2 - L3)**4
        vec( 59-k) = L1 * L3 * (L1 - L3)**4
        vec( 60-k) = L1 * L4 * (L1 - L4)**4
        vec( 61-k) = L2 * L4 * (L2 - L4)**4
        vec( 62-k) = L3 * L4 * (L3 - L4)**4
!  face DOF
        vec( 63-k) = L2 * L3 * L4 * (L2 - L3)**3
        vec( 67-k) = L1 * L3 * L4 * (L1 - L3)**3
        vec( 71-k) = L1 * L2 * L4 * (L1 - L2)**3
        vec( 75-k) = L1 * L2 * L3 * (L1 - L2)**3

        vec( 64-k) = L2 * L3 * L4 * (L2 - L3)**2 * (L3 - L4)
        vec( 68-k) = L1 * L3 * L4 * (L1 - L3)**2 * (L3 - L4)
        vec( 72-k) = L1 * L2 * L4 * (L1 - L2)**2 * (L2 - L4)
        vec( 76-k) = L1 * L2 * L3 * (L1 - L2)**2 * (L2 - L3)

        vec( 65-k) = L2 * L3 * L4 * (L2 - L3) * (L3 - L4)**2
        vec( 69-k) = L1 * L3 * L4 * (L1 - L3) * (L3 - L4)**2
        vec( 73-k) = L1 * L2 * L4 * (L1 - L2) * (L2 - L4)**2
        vec( 77-k) = L1 * L2 * L3 * (L1 - L2) * (L2 - L3)**2

        vec( 66-k) = L2 * L3 * L4 * (L3 - L4)**3
        vec( 70-k) = L1 * L3 * L4 * (L3 - L4)**3
        vec( 74-k) = L1 * L2 * L4 * (L2 - L4)**3
        vec( 78-k) = L1 * L2 * L3 * (L2 - L3)**3
!  volume DOF
        vec( 79-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2
        vec( 80-k) = L1 * L2 * L3 * L4 * (L2 - L3)**2
        vec( 81-k) = L1 * L2 * L3 * L4 * (L3 - L4)**2
        vec( 82-k) = L1 * L2 * L3 * L4 * (L1 - L2)*(L2 - L3)
        vec( 83-k) = L1 * L2 * L3 * L4 * (L1 - L2)*(L3 - L4)
        vec( 84-k) = L1 * L2 * L3 * L4 * (L2 - L3)*(L3 - L4)

    case (7)   !  shape functions for p = 7
!  edge DOF
        vec( 85-k) = L1 * L2 * (L1 - L2)**5
        vec( 86-k) = L2 * L3 * (L2 - L3)**5
        vec( 87-k) = L1 * L3 * (L1 - L3)**5
        vec( 88-k) = L1 * L4 * (L1 - L4)**5
        vec( 89-k) = L2 * L4 * (L2 - L4)**5
        vec( 90-k) = L3 * L4 * (L3 - L4)**5
!  face DOF
        vec( 91-k) = L2 * L3 * L4 * (L2 - L3)**4
        vec( 96-k) = L1 * L3 * L4 * (L1 - L3)**4
        vec(101-k) = L1 * L2 * L4 * (L1 - L2)**4
        vec(106-k) = L1 * L2 * L3 * (L1 - L2)**4

        vec( 92-k) = L2 * L3 * L4 * (L2 - L3)**3 * (L3 - L4)
        vec( 97-k) = L1 * L3 * L4 * (L1 - L3)**3 * (L3 - L4)
        vec(102-k) = L1 * L2 * L4 * (L1 - L2)**3 * (L2 - L4)
        vec(107-k) = L1 * L2 * L3 * (L1 - L2)**3 * (L2 - L3)

        vec( 93-k) = L2 * L3 * L4 * (L2 - L3)**2 * (L3 - L4)**2
        vec( 98-k) = L1 * L3 * L4 * (L1 - L3)**2 * (L3 - L4)**2
        vec(103-k) = L1 * L2 * L4 * (L1 - L2)**2 * (L2 - L4)**2
        vec(108-k) = L1 * L2 * L3 * (L1 - L2)**2 * (L2 - L3)**2

        vec( 94-k) = L2 * L3 * L4 * (L2 - L3) * (L3 - L4)**3
        vec( 99-k) = L1 * L3 * L4 * (L1 - L3) * (L3 - L4)**3
        vec(104-k) = L1 * L2 * L4 * (L1 - L2) * (L2 - L4)**3
        vec(109-k) = L1 * L2 * L3 * (L1 - L2) * (L2 - L3)**3

        vec( 95-k) = L2 * L3 * L4 * (L3 - L4)**4
        vec(100-k) = L1 * L3 * L4 * (L3 - L4)**4
        vec(105-k) = L1 * L2 * L4 * (L2 - L4)**4
        vec(110-k) = L1 * L2 * L3 * (L2 - L3)**4
!  volume DOF
        vec(111-k) = L1 * L2 * L3 * L4 * (L1 - L2)**3
        vec(112-k) = L1 * L2 * L3 * L4 * (L2 - L3)**3
        vec(113-k) = L1 * L2 * L3 * L4 * (L3 - L4)**3
        vec(114-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2 * (L2 - L3)
        vec(115-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2 * (L3 - L4)
        vec(116-k) = L1 * L2 * L3 * L4 * (L2 - L3)**2 * (L3 - L4)
        vec(117-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L2 - L3)**2
        vec(118-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L3 - L4)**2
        vec(119-k) = L1 * L2 * L3 * L4 * (L2 - L3) * (L3 - L4)**2
        vec(120-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L2 - L3) * (L3 - L4)

    case (8)   !  shape functions for p = 8
!  edge DOF
        vec(121-k) = L1 * L2 * (L1 - L2)**6
        vec(122-k) = L2 * L3 * (L2 - L3)**6
        vec(123-k) = L1 * L3 * (L1 - L3)**6
        vec(124-k) = L1 * L4 * (L1 - L4)**6
        vec(125-k) = L2 * L4 * (L2 - L4)**6
        vec(126-k) = L3 * L4 * (L3 - L4)**6
!  face DOF
        vec(127-k) = L2 * L3 * L4 * (L2 - L3)**5
        vec(133-k) = L1 * L3 * L4 * (L1 - L3)**5
        vec(139-k) = L1 * L2 * L4 * (L1 - L2)**5
        vec(145-k) = L1 * L2 * L3 * (L1 - L2)**5

        vec(128-k) = L2 * L3 * L4 * (L2 - L3)**4 * (L3 - L4)
        vec(134-k) = L1 * L3 * L4 * (L1 - L3)**4 * (L3 - L4)
        vec(140-k) = L1 * L2 * L4 * (L1 - L2)**4 * (L2 - L4)
        vec(146-k) = L1 * L2 * L3 * (L1 - L2)**4 * (L2 - L3)

        vec(129-k) = L2 * L3 * L4 * (L2 - L3)**3 * (L3 - L4)**2
        vec(135-k) = L1 * L3 * L4 * (L1 - L3)**3 * (L3 - L4)**2
        vec(141-k) = L1 * L2 * L4 * (L1 - L2)**3 * (L2 - L4)**2
        vec(147-k) = L1 * L2 * L3 * (L1 - L2)**3 * (L2 - L3)**2

        vec(130-k) = L2 * L3 * L4 * (L2 - L3)**2 * (L3 - L4)**3
        vec(136-k) = L1 * L3 * L4 * (L1 - L3)**2 * (L3 - L4)**3
        vec(142-k) = L1 * L2 * L4 * (L1 - L2)**2 * (L2 - L4)**3
        vec(148-k) = L1 * L2 * L3 * (L1 - L2)**2 * (L2 - L3)**3

        vec(131-k) = L2 * L3 * L4 * (L2 - L3) * (L3 - L4)**4
        vec(137-k) = L1 * L3 * L4 * (L1 - L3) * (L3 - L4)**4
        vec(143-k) = L1 * L2 * L4 * (L1 - L2) * (L2 - L4)**4
        vec(149-k) = L1 * L2 * L3 * (L1 - L2) * (L2 - L3)**4

        vec(132-k) = L2 * L3 * L4 * (L3 - L4)**5
        vec(138-k) = L1 * L3 * L4 * (L3 - L4)**5
        vec(144-k) = L1 * L2 * L4 * (L2 - L4)**5
        vec(150-k) = L1 * L2 * L3 * (L2 - L3)**5
!  volume DOF
        vec(151-k) = L1 * L2 * L3 * L4 * (L1 - L2)**4
        vec(152-k) = L1 * L2 * L3 * L4 * (L2 - L3)**4
        vec(153-k) = L1 * L2 * L3 * L4 * (L3 - L4)**4
        vec(154-k) = L1 * L2 * L3 * L4 * (L1 - L2)**3 * (L2 - L3)
        vec(155-k) = L1 * L2 * L3 * L4 * (L1 - L2)**3 * (L3 - L4)
        vec(156-k) = L1 * L2 * L3 * L4 * (L2 - L3)**3 * (L3 - L4)
        vec(157-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2 * (L2 - L3)**2
        vec(158-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2 * (L3 - L4)**2
        vec(159-k) = L1 * L2 * L3 * L4 * (L2 - L3)**2 * (L3 - L4)**2
        vec(160-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L2 - L3)**3
        vec(161-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L3 - L4)**3
        vec(162-k) = L1 * L2 * L3 * L4 * (L2 - L3) * (L3 - L4)**3
        vec(163-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2 * (L2 - L3) * (L3 - L4)
        vec(164-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L2 - L3)**2 * (L3 - L4)
        vec(165-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L2 - L3) * (L3 - L4)**2

    case (9)   !  shape functions for p = 9
!  edge DOF
        vec(166-k) = L1 * L2 * (L1 - L2)**7
        vec(167-k) = L2 * L3 * (L2 - L3)**7
        vec(168-k) = L1 * L3 * (L1 - L3)**7
        vec(169-k) = L1 * L4 * (L1 - L4)**7
        vec(170-k) = L2 * L4 * (L2 - L4)**7
        vec(171-k) = L3 * L4 * (L3 - L4)**7
!  face DOF
        vec(172-k) = L2 * L3 * L4 * (L2 - L3)**6
        vec(179-k) = L1 * L3 * L4 * (L1 - L3)**6
        vec(186-k) = L1 * L2 * L4 * (L1 - L2)**6
        vec(193-k) = L1 * L2 * L3 * (L1 - L2)**6

        vec(173-k) = L2 * L3 * L4 * (L2 - L3)**5 * (L3 - L4)
        vec(180-k) = L1 * L3 * L4 * (L1 - L3)**5 * (L3 - L4)
        vec(187-k) = L1 * L2 * L4 * (L1 - L2)**5 * (L2 - L4)
        vec(194-k) = L1 * L2 * L3 * (L1 - L2)**5 * (L2 - L3)

        vec(174-k) = L2 * L3 * L4 * (L2 - L3)**4 * (L3 - L4)**2
        vec(181-k) = L1 * L3 * L4 * (L1 - L3)**4 * (L3 - L4)**2
        vec(188-k) = L1 * L2 * L4 * (L1 - L2)**4 * (L2 - L4)**2
        vec(195-k) = L1 * L2 * L3 * (L1 - L2)**4 * (L2 - L3)**2

        vec(175-k) = L2 * L3 * L4 * (L2 - L3)**3 * (L3 - L4)**3
        vec(182-k) = L1 * L3 * L4 * (L1 - L3)**3 * (L3 - L4)**3
        vec(189-k) = L1 * L2 * L4 * (L1 - L2)**3 * (L2 - L4)**3
        vec(196-k) = L1 * L2 * L3 * (L1 - L2)**3 * (L2 - L3)**3

        vec(176-k) = L2 * L3 * L4 * (L2 - L3)**2 * (L3 - L4)**4
        vec(183-k) = L1 * L3 * L4 * (L1 - L3)**2 * (L3 - L4)**4
        vec(190-k) = L1 * L2 * L4 * (L1 - L2)**2 * (L2 - L4)**4
        vec(197-k) = L1 * L2 * L3 * (L1 - L2)**2 * (L2 - L3)**4

        vec(177-k) = L2 * L3 * L4 * (L2 - L3) * (L3 - L4)**5
        vec(184-k) = L1 * L3 * L4 * (L1 - L3) * (L3 - L4)**5
        vec(191-k) = L1 * L2 * L4 * (L1 - L2) * (L2 - L4)**5
        vec(198-k) = L1 * L2 * L3 * (L1 - L2) * (L2 - L3)**5

        vec(178-k) = L2 * L3 * L4 * (L3 - L4)**6
        vec(185-k) = L1 * L3 * L4 * (L3 - L4)**6
        vec(192-k) = L1 * L2 * L4 * (L2 - L4)**6
        vec(199-k) = L1 * L2 * L3 * (L2 - L3)**6
!  volume DOF
        vec(200-k) = L1 * L2 * L3 * L4 * (L1 - L2)**5
        vec(201-k) = L1 * L2 * L3 * L4 * (L2 - L3)**5
        vec(202-k) = L1 * L2 * L3 * L4 * (L3 - L4)**5
        vec(203-k) = L1 * L2 * L3 * L4 * (L1 - L2)**4 * (L2 - L3)
        vec(204-k) = L1 * L2 * L3 * L4 * (L1 - L2)**4 * (L3 - L4)
        vec(205-k) = L1 * L2 * L3 * L4 * (L2 - L3)**4 * (L3 - L4)
        vec(206-k) = L1 * L2 * L3 * L4 * (L1 - L2)**3 * (L2 - L3)**2
        vec(207-k) = L1 * L2 * L3 * L4 * (L1 - L2)**3 * (L3 - L4)**2
        vec(208-k) = L1 * L2 * L3 * L4 * (L2 - L3)**3 * (L3 - L4)**2
        vec(209-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2 * (L2 - L3)**3
        vec(210-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2 * (L3 - L4)**3
        vec(211-k) = L1 * L2 * L3 * L4 * (L2 - L3)**2 * (L3 - L4)**3
        vec(212-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L2 - L3)**4
        vec(213-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L3 - L4)**4
        vec(214-k) = L1 * L2 * L3 * L4 * (L2 - L3) * (L3 - L4)**4
        vec(215-k) = L1 * L2 * L3 * L4 * (L1 - L2)**3 * (L2 - L3) * (L3 - L4)
        vec(216-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L2 - L3)**3 * (L3 - L4)
        vec(217-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L2 - L3) * (L3 - L4)**3
        vec(218-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2 * (L2 - L3)**2 * (L3 - L4)
        vec(219-k) = L1 * L2 * L3 * L4 * (L1 - L2)**2 * (L2 - L3) * (L3 - L4)**2
        vec(220-k) = L1 * L2 * L3 * L4 * (L1 - L2) * (L2 - L3)**2 * (L3 - L4)**2
!
! General form of shape functions
!
!  Edge     u = p-2
!        vec(nnn-k) = L1 * L2 * (L1 - L2)**u
!        vec(nnn-k) = L2 * L3 * (L2 - L3)**u
!        vec(nnn-k) = L1 * L3 * (L1 - L3)**u
!        vec(nnn-k) = L1 * L4 * (L1 - L4)**u
!        vec(nnn-k) = L2 * L4 * (L2 - L4)**u
!        vec(nnn-k) = L3 * L4 * (L3 - L4)**u
!
!  Face     u + v = p-3
!        vec(nnn-k) = L2 * L3 * L4 * (L2 - L3)**u * (L3 - L4)**v
!        vec(nnn-k) = L1 * L3 * L4 * (L1 - L3)**u * (L3 - L4)**v
!        vec(nnn-k) = L1 * L2 * L4 * (L1 - L2)**u * (L2 - L4)**v
!        vec(nnn-k) = L1 * L2 * L3 * (L1 - L2)**u * (L2 - L3)**v
!
! Volume    u + v + w = p-4
!        vec(nnn-k) = L1 * L2 * L3 * L4 * (L1 - L2)**u * (L2 - L3)**v * (L3 - L4)**w
! 

        end select poly
!
        p = p + 1
      end do
!
      end subroutine shapesc3D
!
!
!
      subroutine gradshapesc3D( l, gl, polylo, polyhi, vec)
      use femtypes
      implicit none
      integer (I4B) :: polylo, polyhi
      real (DP)     :: l(4), gl(4)
      real (DP)     :: vec(:)
      intent (in)   :: l, gl, polylo, polyhi
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
!
!  Output:
!            vec      array of either x,y, or z component of GRAD of shape functions 
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
      real (DP) L1, L2, L3, L4, gL1, gL2, gL3, gL4
!
      p=polylo
!
!  intermediate variables to avoid l(3) and gl(3) arrays
      L1 = l(1)
      L2 = l(2)
      L3 = l(3)
      L4 = l(4)
!
! intermediate variables to avoid gl arrays
! use explicit substitutions for gL: 
!  
      gL1 = gl(1)
      gL2 = gl(2)
      gL3 = gl(3)
      gL4 = gl(4)

!
! k depends on the formulations of shape functions as above
      if (polylo .gt. 1) then
        k = ( p+1 )*( p+2 )*( p+3 )/2
      else
        k = 0
      end if
!
!-------------------------------------------------------------------------------!
!  BEGIN LOOP TO CALCULATE GRADIENT OF SHAPE FUNCTIONS FOR POLYNOMIAL DEGREE p  !
!-------------------------------------------------------------------------------!
      do
        if (p .gt. polyhi) exit    ! stopping criteria
!  Begin "case" control structure to select the appropriate shape functions for p
        poly: select case (p)
    case (1)  ! grad of shape functions for p = 1
        vec(  1-k) = gL1
        vec(  2-k) = gL2
        vec(  3-k) = gL3
        vec(  4-k) = gL4

    case (2)  ! grad of shape functions for p = 2
!  edge DOF
        vec(  5-k) = L1*gL2 + L2*gL1
        vec(  6-k) = L2*gL3 + L3*gL2
        vec(  7-k) = L1*gL3 + L3*gL1
        vec(  8-k) = L1*gL4 + L4*gL1
        vec(  9-k) = L2*gL4 + L4*gL2
        vec( 10-k) = L3*gL4 + L4*gL3

    case (3)  ! grad of shape functions for p = 3
!  edge DOF
        vec( 11-k) = (                                   &
     &               (L1 - L2) * ( gL1 * L2 + L1 * gL2 )                &
     &              + L1 * L2 * (gL1 - gL2) )
        vec( 12-k) = (                                   &
     &               (L2 - L3) * ( gL2 * L3 + L2 * gL3 )                &
     &              + L2 * L3 * (gL2 - gL3) )
        vec( 13-k) = (                                   &
     &               (L1 - L3) * ( gL1 * L3 + L1 * gL3 )                &
     &              + L1 * L3 * (gL1 - gL3) )
        vec( 14-k) = (                                   &
     &               (L1 - L4) * ( gL1 * L4 + L1 * gL4 )                &
     &              + L1 * L4 * (gL1 - gL4) )
        vec( 15-k) = (                                   &
     &               (L2 - L4) * ( gL2 * L4 + L2 * gL4 )                &
     &              + L2 * L4 * (gL2 - gL4) )
        vec( 16-k) = (                                   &
     &               (L3 - L4) * ( gL3 * L4 + L3 * gL4 )                &
     &              + L3 * L4 * (gL3 - gL4) )
!  face DOF
        vec(17-k) = gL2 * L4 * L3 + L2 * gL4 * L3 + L2 * L4 * gL3
        vec(18-k) = gL1 * L3 * L4 + L1 * gL3 * L4 + L1 * L3 * gL4
        vec(19-k) = gL1 * L4 * L2 + L1 * gL4 * L2 + L1 * L4 * gL2
        vec(20-k) = gL1 * L2 * L3 + L1 * gL2 * L3 + L1 * L2 * gL3

    case (4)  ! grad of shape functions for p = 4
!  edge DOF
        vec( 21-k) = (L1 - L2) * (                                   &
     &               (L1 - L2) * ( gL1 * L2 + L1 * gL2 )                &
     &              + L1 * L2 * 2 * (gL1 - gL2) )
        vec( 22-k) = (L2 - L3) * (                                   &
     &               (L2 - L3) * ( gL2 * L3 + L2 * gL3 )                &
     &              + L2 * L3 * 2 * (gL2 - gL3) )
        vec( 23-k) = (L1 - L3) * (                                   &
     &               (L1 - L3) * ( gL1 * L3 + L1 * gL3 )                &
     &              + L1 * L3 * 2 * (gL1 - gL3) )
        vec( 24-k) = (L1 - L4) * (                                   &
     &               (L1 - L4) * ( gL1 * L4 + L1 * gL4 )                &
     &              + L1 * L4 * 2 * (gL1 - gL4) )
        vec( 25-k) = (L2 - L4) * (                                   &
     &               (L2 - L4) * ( gL2 * L4 + L2 * gL4 )                &
     &              + L2 * L4 * 2 * (gL2 - gL4) )
        vec( 26-k) = (L3 - L4) * (                                   &
     &               (L3 - L4) * ( gL3 * L4 + L3 * gL4 )                &
     &              + L3 * L4 * 2 * (gL3 - gL4) )
!  face DOF
        vec( 27-k) = (                                          &
     &               (L2 - L3) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)  &
     &              + L2 * L3 * L4 * ((gL2 - gL3) ) )
        vec( 29-k) = (                                          &
     &               (L1 - L3) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)  &
     &              + L1 * L3 * L4 * ((gL1 - gL3) ) )
        vec( 31-k) = (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)  &
     &              + L1 * L2 * L4 * ((gL1 - gL2) ) )
        vec( 33-k) = (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)  &
     &              + L1 * L2 * L3 * ((gL1 - gL2) ) )

        vec( 28-k) = (                                          &
     &               (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)      &
     &              + L2 * L3 * L4 * ((gL3 - gL4) ) )
        vec( 30-k) = (                                          &
     &               (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)      &
     &              + L1 * L3 * L4 * ((gL3 - gL4) ) )
        vec( 32-k) = (                                          &
     &               (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)      &
     &              + L1 * L2 * L4 * ((gL2 - gL4) ) )
        vec( 34-k) = (                                          &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)      &
     &              + L1 * L2 * L3 * ((gL2 - gL3) ) )
!  volume DOF
        vec(35-k) = gL1 * L2 * L3 * L4 + L1 * gL2 * L3 * L4 + L1 * L2 * gL3 * L4 + L1 * L2 * L3 * gL4

    case (5)   !  grad of shape functions for p = 5
!  edge DOF
        vec( 36-k) = (L1 - L2)**2 * (                                   &
     &               (L1 - L2) * ( gL1 * L2 + L1 * gL2 )                &
     &              + L1 * L2 * 3 * (gL1 - gL2) )
        vec( 37-k) = (L2 - L3)**2 * (                                   &
     &               (L2 - L3) * ( gL2 * L3 + L2 * gL3 )                &
     &              + L2 * L3 * 3 * (gL2 - gL3) )
        vec( 38-k) = (L1 - L3)**2 * (                                   &
     &               (L1 - L3) * ( gL1 * L3 + L1 * gL3 )                &
     &              + L1 * L3 * 3 * (gL1 - gL3) )
        vec( 39-k) = (L1 - L4)**2 * (                                   &
     &               (L1 - L4) * ( gL1 * L4 + L1 * gL4 )                &
     &              + L1 * L4 * 3 * (gL1 - gL4) )
        vec( 40-k) = (L2 - L4)**2 * (                                   &
     &               (L2 - L4) * ( gL2 * L4 + L2 * gL4 )                &
     &              + L2 * L4 * 3 * (gL2 - gL4) )
        vec( 41-k) = (L3 - L4)**2 * (                                   &
     &               (L3 - L4) * ( gL3 * L4 + L3 * gL4 )                &
     &              + L3 * L4 * 3 * (gL3 - gL4) )
!  face DOF
        vec( 42-k) = (L2 - L3) * (                                          &
     &               (L2 - L3) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)  &
     &              + L2 * L3 * L4 * (2 * (gL2 - gL3) ) )
        vec( 45-k) = (L1 - L3) * (                                          &
     &               (L1 - L3) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)  &
     &              + L1 * L3 * L4 * (2 * (gL1 - gL3) ) )
        vec( 48-k) = (L1 - L2) * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)  &
     &              + L1 * L2 * L4 * (2 * (gL1 - gL2) ) )
        vec( 51-k) = (L1 - L2) * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)  &
     &              + L1 * L2 * L3 * (2 * (gL1 - gL2) ) )

        vec( 43-k) = (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * ((gL2 - gL3) * (L3 - L4) + (L2 - L3) * (gL3 - gL4) ) )
        vec( 46-k) = (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * ((gL1 - gL3) * (L3 - L4) + (L1 - L3) * (gL3 - gL4) ) )
        vec( 49-k) = (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * ((gL1 - gL2) * (L2 - L4) + (L1 - L2) * (gL2 - gL4) ) )
        vec( 52-k) = (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * ((gL1 - gL2) * (L2 - L3) + (L1 - L2) * (gL2 - gL3) ) )

        vec( 44-k) = (L3 - L4) * (                                          &
     &               (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)      &
     &              + L2 * L3 * L4 * (2 * (gL3 - gL4) ) )
        vec( 47-k) = (L3 - L4) * (                                          &
     &               (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)      &
     &              + L1 * L3 * L4 * (2 * (gL3 - gL4) ) )
        vec( 50-k) = (L2 - L4) * (                                          &
     &               (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)      &
     &              + L1 * L2 * L4 * (2 * (gL2 - gL4) ) )
        vec( 53-k) = (L2 - L3) * (                                          &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)      &
     &              + L1 * L2 * L3 * (2 * (gL2 - gL3) ) )
!  volume DOF
        vec( 54-k) = (                                                                 &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) ) )
        vec( 55-k) = (                                                                 &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * ((gL2 - gL3) ) )
        vec( 56-k) = (                                                                 &
     &               (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * ((gL3 - gL4)) )

    case (6)   !  grad of shape functions for p = 6
!  edge DOF
        vec( 57-k) = (L1 - L2)**3 * (                                   &
     &               (L1 - L2) * ( gL1 * L2 + L1 * gL2 )                &
     &              + L1 * L2 * 4 * (gL1 - gL2) )
        vec( 58-k) = (L2 - L3)**3 * (                                   &
     &               (L2 - L3) * ( gL2 * L3 + L2 * gL3 )                &
     &              + L2 * L3 * 4 * (gL2 - gL3) )
        vec( 59-k) = (L1 - L3)**3 * (                                   &
     &               (L1 - L3) * ( gL1 * L3 + L1 * gL3 )                &
     &              + L1 * L3 * 4 * (gL1 - gL3) )
        vec( 60-k) = (L1 - L4)**3 * (                                   &
     &               (L1 - L4) * ( gL1 * L4 + L1 * gL4 )                &
     &              + L1 * L4 * 4 * (gL1 - gL4) )
        vec( 61-k) = (L2 - L4)**3 * (                                   &
     &               (L2 - L4) * ( gL2 * L4 + L2 * gL4 )                &
     &              + L2 * L4 * 4 * (gL2 - gL4) )
        vec( 62-k) = (L3 - L4)**3 * (                                   &
     &               (L3 - L4) * ( gL3 * L4 + L3 * gL4 )                &
     &              + L3 * L4 * 4 * (gL3 - gL4) )
!  face DOF
        vec( 63-k) = (L2 - L3)**2 * (                                          &
     &               (L2 - L3) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)  &
     &              + L2 * L3 * L4 * (3 * (gL2 - gL3) ) )
        vec( 67-k) = (L1 - L3)**2 * (                                          &
     &               (L1 - L3) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)  &
     &              + L1 * L3 * L4 * (3 * (gL1 - gL3) ) )
        vec( 71-k) = (L1 - L2)**2 * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)  &
     &              + L1 * L2 * L4 * (3 * (gL1 - gL2) ) )
        vec( 75-k) = (L1 - L2)**2 * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)  &
     &              + L1 * L2 * L3 * (3 * (gL1 - gL2) ) )

        vec( 64-k) = (L2 - L3) * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (2 * (gL2 - gL3) * (L3 - L4) + (L2 - L3) * (gL3 - gL4) ) )
        vec( 68-k) = (L1 - L3) * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (2 * (gL1 - gL3) * (L3 - L4) + (L1 - L3) * (gL3 - gL4) ) )
        vec( 72-k) = (L1 - L2) * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (2 * (gL1 - gL2) * (L2 - L4) + (L1 - L2) * (gL2 - gL4) ) )
        vec( 76-k) = (L1 - L2) * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (2 * (gL1 - gL2) * (L2 - L3) + (L1 - L2) * (gL2 - gL3) ) )

        vec( 65-k) = (L3 - L4) * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * ( (gL2 - gL3) * (L3 - L4) + 2 * (L2 - L3) * (gL3 - gL4) ) )
        vec( 69-k) = (L3 - L4) * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * ( (gL1 - gL3) * (L3 - L4) + 2 * (L1 - L3) * (gL3 - gL4) ) )
        vec( 73-k) = (L2 - L4) * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * ( (gL1 - gL2) * (L2 - L4) + 2 * (L1 - L2) * (gL2 - gL4) ) )
        vec( 77-k) = (L2 - L3) * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * ( (gL1 - gL2) * (L2 - L3) + 2 * (L1 - L2) * (gL2 - gL3) ) )

        vec( 66-k) = (L3 - L4)**2 * (                                          &
     &               (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)      &
     &              + L2 * L3 * L4 * (3 * (gL3 - gL4) ) )
        vec( 70-k) = (L3 - L4)**2 * (                                          &
     &               (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)      &
     &              + L1 * L3 * L4 * (3 * (gL3 - gL4) ) )
        vec( 74-k) = (L2 - L4)**2 * (                                          &
     &               (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)      &
     &              + L1 * L2 * L4 * (3 * (gL2 - gL4) ) )
        vec( 78-k) = (L2 - L3)**2 * (                                          &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)      &
     &              + L1 * L2 * L3 * (3 * (gL2 - gL3) ) )
!  volume DOF
        vec( 79-k) = (L1 - L2) * (                                                                 &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) ) )
        vec( 80-k) = (L2 - L3) * (                                                                 &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (2 * (gL2 - gL3) ) )
        vec( 81-k) = (L3 - L4) * (                                                                 &
     &               (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (2 * (gL3 - gL4)) )
        vec( 82-k) = (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) + (L1 - L2) * (gL2 - gL3) ) )
        vec( 83-k) = (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L3 - L4) + (L1 - L2) * (gL3 - gL4)) )
        vec( 84-k) = (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL2 - gL3) * (L3 - L4) + (L2 - L3) * (gL3 - gL4)) )

    case (7)   !  grad of shape functions for p = 7
!  edge DOF
        vec( 85-k) = (L1 - L2)**4 * (                                   &
     &               (L1 - L2) * ( gL1 * L2 + L1 * gL2 )                &
     &              + L1 * L2 * 5 * (gL1 - gL2) )
        vec( 86-k) = (L2 - L3)**4 * (                                   &
     &               (L2 - L3) * ( gL2 * L3 + L2 * gL3 )                &
     &              + L2 * L3 * 5 * (gL2 - gL3) )
        vec( 87-k) = (L1 - L3)**4 * (                                   &
     &               (L1 - L3) * ( gL1 * L3 + L1 * gL3 )                &
     &              + L1 * L3 * 5 * (gL1 - gL3) )
        vec( 88-k) = (L1 - L4)**4 * (                                   &
     &               (L1 - L4) * ( gL1 * L4 + L1 * gL4 )                &
     &              + L1 * L4 * 5 * (gL1 - gL4) )
        vec( 89-k) = (L2 - L4)**4 * (                                   &
     &               (L2 - L4) * ( gL2 * L4 + L2 * gL4 )                &
     &              + L2 * L4 * 5 * (gL2 - gL4) )
        vec( 90-k) = (L3 - L4)**4 * (                                   &
     &               (L3 - L4) * ( gL3 * L4 + L3 * gL4 )                &
     &              + L3 * L4 * 5 * (gL3 - gL4) )
!  face DOF
        vec( 91-k) = (L2 - L3)**3 * (                                          &
     &               (L2 - L3) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)  &
     &              + L2 * L3 * L4 * (4 * (gL2 - gL3) ) )
        vec( 96-k) = (L1 - L3)**3 * (                                          &
     &               (L1 - L3) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)  &
     &              + L1 * L3 * L4 * (4 * (gL1 - gL3) ) )
        vec(101-k) = (L1 - L2)**3 * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)  &
     &              + L1 * L2 * L4 * (4 * (gL1 - gL2) ) )
        vec(106-k) = (L1 - L2)**3 * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)  &
     &              + L1 * L2 * L3 * (4 * (gL1 - gL2) ) )

        vec( 92-k) = (L2 - L3)**2 * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (3 * (gL2 - gL3) * (L3 - L4) + (L2 - L3) * (gL3 - gL4) ) )
        vec( 97-k) = (L1 - L3)**2 * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (3 * (gL1 - gL3) * (L3 - L4) + (L1 - L3) * (gL3 - gL4) ) )
        vec(102-k) = (L1 - L2)**2 * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (3 * (gL1 - gL2) * (L2 - L4) + (L1 - L2) * (gL2 - gL4) ) )
        vec(107-k) = (L1 - L2)**2 * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (3 * (gL1 - gL2) * (L2 - L3) + (L1 - L2) * (gL2 - gL3) ) )

        vec( 93-k) = (L2 - L3)*(L3 - L4) * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (2 * (gL2 - gL3) * (L3 - L4) + 2 * (L2 - L3) * (gL3 - gL4) ) )
        vec( 98-k) = (L1 - L3)*(L3 - L4) * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (2 * (gL1 - gL3) * (L3 - L4) + 2 * (L1 - L3) * (gL3 - gL4) ) )
        vec(103-k) = (L1 - L2)*(L2 - L4) * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (2 * (gL1 - gL2) * (L2 - L4) + 2 * (L1 - L2) * (gL2 - gL4) ) )
        vec(108-k) = (L1 - L2)*(L2 - L3) * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (2 * (gL1 - gL2) * (L2 - L3) + 2 * (L1 - L2) * (gL2 - gL3) ) )

        vec( 94-k) = (L3 - L4)**2 * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * ( (gL2 - gL3) * (L3 - L4) + 3 * (L2 - L3) * (gL3 - gL4) ) )
        vec( 99-k) = (L3 - L4)**2 * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * ( (gL1 - gL3) * (L3 - L4) + 3 * (L1 - L3) * (gL3 - gL4) ) )
        vec(104-k) = (L2 - L4)**2 * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * ( (gL1 - gL2) * (L2 - L4) + 3 * (L1 - L2) * (gL2 - gL4) ) )
        vec(109-k) = (L2 - L3)**2 * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * ( (gL1 - gL2) * (L2 - L3) + 3 * (L1 - L2) * (gL2 - gL3) ) )

        vec( 95-k) = (L3 - L4)**3 * (                                          &
     &               (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)      &
     &              + L2 * L3 * L4 * (4 * (gL3 - gL4) ) )
        vec(100-k) = (L3 - L4)**3 * (                                          &
     &               (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)      &
     &              + L1 * L3 * L4 * (4 * (gL3 - gL4) ) )
        vec(105-k) = (L2 - L4)**3 * (                                          &
     &               (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)      &
     &              + L1 * L2 * L4 * (4 * (gL2 - gL4) ) )
        vec(110-k) = (L2 - L3)**3 * (                                          &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)      &
     &              + L1 * L2 * L3 * (4 * (gL2 - gL3) ) )
!  volume DOF
        vec(111-k) = (L1 - L2)**2 * (                                                                 &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (3 * (gL1 - gL2) ) )
        vec(112-k) = (L2 - L3)**2 * (                                                                 &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (3 * (gL2 - gL3) ) )
        vec(113-k) = (L3 - L4)**2 * (                                                                 &
     &               (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (3 * (gL3 - gL4)) )
        vec(114-k) = (L1 - L2) * (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) * (L2 - L3) + (L1 - L2) * (gL2 - gL3) ) )
        vec(115-k) = (L1 - L2) * (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) * (L3 - L4) + (L1 - L2) * (gL3 - gL4)) )
        vec(116-k) = (L2 - L3) * (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (2 * (gL2 - gL3) * (L3 - L4) + (L2 - L3) * (gL3 - gL4)) )
        vec(117-k) = (L2 - L3) * (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) + 2 * (L1 - L2) * (gL2 - gL3) ) )
        vec(118-k) = (L3 - L4) * (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L3 - L4) + 2 * (L1 - L2) * (gL3 - gL4)) )
        vec(119-k) = (L3 - L4) * (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL2 - gL3) * (L3 - L4) + 2 * (L2 - L3) * (gL3 - gL4)) )
        vec(120-k) = (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                    &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) * (L3 - L4) + (L1 - L2) * (gL2 - gL3) * (L3 - L4) + (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )

    case (8)   !  grad of shape functions for p = 8
!  edge DOF
        vec(121-k) = (L1 - L2)**5 * (                                   &
     &               (L1 - L2) * ( gL1 * L2 + L1 * gL2 )                &
     &              + L1 * L2 * 6 * (gL1 - gL2) )
        vec(122-k) = (L2 - L3)**5 * (                                   &
     &               (L2 - L3) * ( gL2 * L3 + L2 * gL3 )                &
     &              + L2 * L3 * 6 * (gL2 - gL3) )
        vec(123-k) = (L1 - L3)**5 * (                                   &
     &               (L1 - L3) * ( gL1 * L3 + L1 * gL3 )                &
     &              + L1 * L3 * 6 * (gL1 - gL3) )
        vec(124-k) = (L1 - L4)**5 * (                                   &
     &               (L1 - L4) * ( gL1 * L4 + L1 * gL4 )                &
     &              + L1 * L4 * 6 * (gL1 - gL4) )
        vec(125-k) = (L2 - L4)**5 * (                                   &
     &               (L2 - L4) * ( gL2 * L4 + L2 * gL4 )                &
     &              + L2 * L4 * 6 * (gL2 - gL4) )
        vec(126-k) = (L3 - L4)**5 * (                                   &
     &               (L3 - L4) * ( gL3 * L4 + L3 * gL4 )                &
     &              + L3 * L4 * 6 * (gL3 - gL4) )
!  face DOF
        vec(127-k) = (L2 - L3)**4 * (                                          &
     &               (L2 - L3) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)  &
     &              + L2 * L3 * L4 * (5 * (gL2 - gL3) ) )
        vec(133-k) = (L1 - L3)**4 * (                                          &
     &               (L1 - L3) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)  &
     &              + L1 * L3 * L4 * (5 * (gL1 - gL3) ) )
        vec(139-k) = (L1 - L2)**4 * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)  &
     &              + L1 * L2 * L4 * (5 * (gL1 - gL2) ) )
        vec(145-k) = (L1 - L2)**4 * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)  &
     &              + L1 * L2 * L3 * (5 * (gL1 - gL2) ) )

        vec(128-k) = (L2 - L3)**3 * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (4 * (gL2 - gL3) * (L3 - L4) + (L2 - L3) * (gL3 - gL4) ) )
        vec(134-k) = (L1 - L3)**3 * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (4 * (gL1 - gL3) * (L3 - L4) + (L1 - L3) * (gL3 - gL4) ) )
        vec(140-k) = (L1 - L2)**3 * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (4 * (gL1 - gL2) * (L2 - L4) + (L1 - L2) * (gL2 - gL4) ) )
        vec(146-k) = (L1 - L2)**3 * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (4 * (gL1 - gL2) * (L2 - L3) + (L1 - L2) * (gL2 - gL3) ) )

        vec(129-k) = (L2 - L3)**2*(L3 - L4) * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (3 * (gL2 - gL3) * (L3 - L4) + 2 * (L2 - L3) * (gL3 - gL4) ) )
        vec(135-k) = (L1 - L3)**2*(L3 - L4) * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (3 * (gL1 - gL3) * (L3 - L4) + 2 * (L1 - L3) * (gL3 - gL4) ) )
        vec(141-k) = (L1 - L2)**2*(L2 - L4) * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (3 * (gL1 - gL2) * (L2 - L4) + 2 * (L1 - L2) * (gL2 - gL4) ) )
        vec(147-k) = (L1 - L2)**2*(L2 - L3) * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (3 * (gL1 - gL2) * (L2 - L3) + 2 * (L1 - L2) * (gL2 - gL3) ) )

        vec(130-k) = (L2 - L3)*(L3 - L4)**2 * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (2 * (gL2 - gL3) * (L3 - L4) + 3 * (L2 - L3) * (gL3 - gL4) ) )
        vec(136-k) = (L1 - L3)*(L3 - L4)**2 * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (2 * (gL1 - gL3) * (L3 - L4) + 3 * (L1 - L3) * (gL3 - gL4) ) )
        vec(142-k) = (L1 - L2)*(L2 - L4)**2 * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (2 * (gL1 - gL2) * (L2 - L4) + 3 * (L1 - L2) * (gL2 - gL4) ) )
        vec(148-k) = (L1 - L2)*(L2 - L3)**2 * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (2 * (gL1 - gL2) * (L2 - L3) + 3 * (L1 - L2) * (gL2 - gL3) ) )

        vec(131-k) = (L3 - L4)**3 * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * ( (gL2 - gL3) * (L3 - L4) + 4 * (L2 - L3) * (gL3 - gL4) ) )
        vec(137-k) = (L3 - L4)**3 * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * ( (gL1 - gL3) * (L3 - L4) + 4 * (L1 - L3) * (gL3 - gL4) ) )
        vec(143-k) = (L2 - L4)**3 * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * ( (gL1 - gL2) * (L2 - L4) + 4 * (L1 - L2) * (gL2 - gL4) ) )
        vec(149-k) = (L2 - L3)**3 * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * ( (gL1 - gL2) * (L2 - L3) + 4 * (L1 - L2) * (gL2 - gL3) ) )

        vec(132-k) = (L3 - L4)**4 * (                                          &
     &               (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)      &
     &              + L2 * L3 * L4 * (5 * (gL3 - gL4) ) )
        vec(138-k) = (L3 - L4)**4 * (                                          &
     &               (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)      &
     &              + L1 * L3 * L4 * (5 * (gL3 - gL4) ) )
        vec(144-k) = (L2 - L4)**4 * (                                          &
     &               (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)      &
     &              + L1 * L2 * L4 * (5 * (gL2 - gL4) ) )
        vec(150-k) = (L2 - L3)**4 * (                                          &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)      &
     &              + L1 * L2 * L3 * (5 * (gL2 - gL3) ) )
!  volume DOF
        vec(151-k) = (L1 - L2)**3 * (                                                                 &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (4 * (gL1 - gL2) ) )
        vec(152-k) = (L2 - L3)**3 * (                                                                 &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (4 * (gL2 - gL3) ) )
        vec(153-k) = (L3 - L4)**3 * (                                                                 &
     &               (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (4 * (gL3 - gL4)) )
     
        vec(154-k) = (L1 - L2)**2 * (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (3 * (gL1 - gL2) * (L2 - L3) + (L1 - L2) * (gL2 - gL3) ) )
        vec(155-k) = (L1 - L2)**2 * (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (3 * (gL1 - gL2) * (L3 - L4) + (L1 - L2) * (gL3 - gL4)) )
        vec(156-k) = (L2 - L3)**2 * (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (3 * (gL2 - gL3) * (L3 - L4) + (L2 - L3) * (gL3 - gL4)) )

        vec(157-k) = (L1 - L2)*(L2 - L3) * (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) * (L2 - L3) + 2 * (L1 - L2) * (gL2 - gL3) ) )
        vec(158-k) = (L1 - L2)*(L3 - L4) * (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) * (L3 - L4) + 2 * (L1 - L2) * (gL3 - gL4)) )
        vec(159-k) = (L2 - L3)*(L3 - L4) * (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (2 * (gL2 - gL3) * (L3 - L4) + 2 * (L2 - L3) * (gL3 - gL4)) )

        vec(160-k) = (L2 - L3)**2 * (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) + 3 * (L1 - L2) * (gL2 - gL3) ) )
        vec(161-k) = (L3 - L4)**2 * (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L3 - L4) + 3 * (L1 - L2) * (gL3 - gL4)) )
        vec(162-k) = (L3 - L4)**2 * (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL2 - gL3) * (L3 - L4) + 3 * (L2 - L3) * (gL3 - gL4)) )

        vec(163-k) = (L1 - L2) * (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) * (L2 - L3) * (L3 - L4) + (L1 - L2) * (gL2 - gL3) * (L3 - L4) + (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )
        vec(164-k) = (L2 - L3) * (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) * (L3 - L4) + 2 * (L1 - L2) * (gL2 - gL3) * (L3 - L4) + (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )
        vec(165-k) = (L3 - L4) * (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) * (L3 - L4) + (L1 - L2) * (gL2 - gL3) * (L3 - L4) + 2 * (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )

    case (9)   !  grad of shape functions for p = 9
!  edge DOF
        vec(166-k) = (L1 - L2)**6 * (                                   &
     &               (L1 - L2) * ( gL1 * L2 + L1 * gL2 )                &
     &              + L1 * L2 * 7 * (gL1 - gL2) )
        vec(167-k) = (L2 - L3)**6 * (                                   &
     &               (L2 - L3) * ( gL2 * L3 + L2 * gL3 )                &
     &              + L2 * L3 * 7 * (gL2 - gL3) )
        vec(168-k) = (L1 - L3)**6 * (                                   &
     &               (L1 - L3) * ( gL1 * L3 + L1 * gL3 )                &
     &              + L1 * L3 * 7 * (gL1 - gL3) )
        vec(169-k) = (L1 - L4)**6 * (                                   &
     &               (L1 - L4) * ( gL1 * L4 + L1 * gL4 )                &
     &              + L1 * L4 * 7 * (gL1 - gL4) )
        vec(170-k) = (L2 - L4)**6 * (                                   &
     &               (L2 - L4) * ( gL2 * L4 + L2 * gL4 )                &
     &              + L2 * L4 * 7 * (gL2 - gL4) )
        vec(171-k) = (L3 - L4)**6 * (                                   &
     &               (L3 - L4) * ( gL3 * L4 + L3 * gL4 )                &
     &              + L3 * L4 * 7 * (gL3 - gL4) )
!  face DOF
        vec(172-k) = (L2 - L3)**5 * (                                          &
     &               (L2 - L3) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)  &
     &              + L2 * L3 * L4 * (6 * (gL2 - gL3) ) )
        vec(179-k) = (L1 - L3)**5 * (                                          &
     &               (L1 - L3) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)  &
     &              + L1 * L3 * L4 * (6 * (gL1 - gL3) ) )
        vec(186-k) = (L1 - L2)**5 * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)  &
     &              + L1 * L2 * L4 * (6 * (gL1 - gL2) ) )
        vec(193-k) = (L1 - L2)**5 * (                                          &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)  &
     &              + L1 * L2 * L3 * (6 * (gL1 - gL2) ) )

        vec(173-k) = (L2 - L3)**4 * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (5 * (gL2 - gL3) * (L3 - L4) + (L2 - L3) * (gL3 - gL4) ) )
        vec(180-k) = (L1 - L3)**4 * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (5 * (gL1 - gL3) * (L3 - L4) + (L1 - L3) * (gL3 - gL4) ) )
        vec(187-k) = (L1 - L2)**4 * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (5 * (gL1 - gL2) * (L2 - L4) + (L1 - L2) * (gL2 - gL4) ) )
        vec(194-k) = (L1 - L2)**4 * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (5 * (gL1 - gL2) * (L2 - L3) + (L1 - L2) * (gL2 - gL3) ) )

        vec(174-k) = (L2 - L3)**3*(L3 - L4) * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (4 * (gL2 - gL3) * (L3 - L4) + 2 * (L2 - L3) * (gL3 - gL4) ) )
        vec(181-k) = (L1 - L3)**3*(L3 - L4) * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (4 * (gL1 - gL3) * (L3 - L4) + 2 * (L1 - L3) * (gL3 - gL4) ) )
        vec(188-k) = (L1 - L2)**3*(L2 - L4) * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (4 * (gL1 - gL2) * (L2 - L4) + 2 * (L1 - L2) * (gL2 - gL4) ) )
        vec(195-k) = (L1 - L2)**3*(L2 - L3) * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (4 * (gL1 - gL2) * (L2 - L3) + 2 * (L1 - L2) * (gL2 - gL3) ) )

        vec(175-k) = (L2 - L3)**2*(L3 - L4)**2 * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (3 * (gL2 - gL3) * (L3 - L4) + 3 * (L2 - L3) * (gL3 - gL4) ) )
        vec(182-k) = (L1 - L3)**2*(L3 - L4)**2 * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (3 * (gL1 - gL3) * (L3 - L4) + 3 * (L1 - L3) * (gL3 - gL4) ) )
        vec(189-k) = (L1 - L2)**2*(L2 - L4)**2 * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (3 * (gL1 - gL2) * (L2 - L4) + 3 * (L1 - L2) * (gL2 - gL4) ) )
        vec(196-k) = (L1 - L2)**2*(L2 - L3)**2 * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (3 * (gL1 - gL2) * (L2 - L3) + 3 * (L1 - L2) * (gL2 - gL3) ) )

        vec(176-k) = (L2 - L3)*(L3 - L4)**3 * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * (2 * (gL2 - gL3) * (L3 - L4) + 4 * (L2 - L3) * (gL3 - gL4) ) )
        vec(183-k) = (L1 - L3)*(L3 - L4)**3 * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * (2 * (gL1 - gL3) * (L3 - L4) + 4 * (L1 - L3) * (gL3 - gL4) ) )
        vec(190-k) = (L1 - L2)*(L2 - L4)**3 * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * (2 * (gL1 - gL2) * (L2 - L4) + 4 * (L1 - L2) * (gL2 - gL4) ) )
        vec(197-k) = (L1 - L2)*(L2 - L3)**3 * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * (2 * (gL1 - gL2) * (L2 - L3) + 4 * (L1 - L2) * (gL2 - gL3) ) )

        vec(177-k) = (L3 - L4)**4 * (                                          &
     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
     &              + L2 * L3 * L4 * ( (gL2 - gL3) * (L3 - L4) + 5 * (L2 - L3) * (gL3 - gL4) ) )
        vec(184-k) = (L3 - L4)**4 * (                                          &
     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
     &              + L1 * L3 * L4 * ( (gL1 - gL3) * (L3 - L4) + 5 * (L1 - L3) * (gL3 - gL4) ) )
        vec(191-k) = (L2 - L4)**4 * (                                          &
     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
     &              + L1 * L2 * L4 * ( (gL1 - gL2) * (L2 - L4) + 5 * (L1 - L2) * (gL2 - gL4) ) )
        vec(198-k) = (L2 - L3)**4 * (                                          &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
     &              + L1 * L2 * L3 * ( (gL1 - gL2) * (L2 - L3) + 5 * (L1 - L2) * (gL2 - gL3) ) )

        vec(178-k) = (L3 - L4)**5 * (                                          &
     &               (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)      &
     &              + L2 * L3 * L4 * (6 * (gL3 - gL4) ) )
        vec(185-k) = (L3 - L4)**5 * (                                          &
     &               (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)      &
     &              + L1 * L3 * L4 * (6 * (gL3 - gL4) ) )
        vec(192-k) = (L2 - L4)**5 * (                                          &
     &               (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)      &
     &              + L1 * L2 * L4 * (6 * (gL2 - gL4) ) )
        vec(199-k) = (L2 - L3)**5 * (                                          &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)      &
     &              + L1 * L2 * L3 * (6 * (gL2 - gL3) ) )
!  volume DOF
        vec(200-k) = (L1 - L2)**4 * (                                                                 &
     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (5 * (gL1 - gL2) ) )
        vec(201-k) = (L2 - L3)**4 * (                                                                 &
     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (5 * (gL2 - gL3) ) )
        vec(202-k) = (L3 - L4)**4 * (                                                                 &
     &               (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
     &              + L1 * L2 * L3 * L4 * (5 * (gL3 - gL4)) )

        vec(203-k) = (L1 - L2)**3 * (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (4 * (gL1 - gL2) * (L2 - L3) + (L1 - L2) * (gL2 - gL3) ) )
        vec(204-k) = (L1 - L2)**3 * (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (4 * (gL1 - gL2) * (L3 - L4) + (L1 - L2) * (gL3 - gL4)) )
        vec(205-k) = (L2 - L3)**3 * (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (4 * (gL2 - gL3) * (L3 - L4) + (L2 - L3) * (gL3 - gL4)) )

        vec(206-k) = (L1 - L2)**2*(L2 - L3) * (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (3 * (gL1 - gL2) * (L2 - L3) + 2 * (L1 - L2) * (gL2 - gL3) ) )
        vec(207-k) = (L1 - L2)**2*(L3 - L4) * (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (3 * (gL1 - gL2) * (L3 - L4) + 2 * (L1 - L2) * (gL3 - gL4)) )
        vec(208-k) = (L2 - L3)**2*(L3 - L4) * (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (3 * (gL2 - gL3) * (L3 - L4) + 2 * (L2 - L3) * (gL3 - gL4)) )

        vec(209-k) = (L1 - L2)*(L2 - L3)**2 * (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) * (L2 - L3) + 3 * (L1 - L2) * (gL2 - gL3) ) )
        vec(210-k) = (L1 - L2)*(L3 - L4)**2 * (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) * (L3 - L4) + 3 * (L1 - L2) * (gL3 - gL4)) )
        vec(211-k) = (L2 - L3)*(L3 - L4)**2 * (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * (2 * (gL2 - gL3) * (L3 - L4) + 3 * (L2 - L3) * (gL3 - gL4)) )

        vec(212-k) = (L2 - L3)**3 * (                                                       &
     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) + 4 * (L1 - L2) * (gL2 - gL3) ) )
        vec(213-k) = (L3 - L4)**3 * (                                                       &
     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L3 - L4) + 4 * (L1 - L2) * (gL3 - gL4)) )
        vec(214-k) = (L3 - L4)**3 * (                                                       &
     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
     &              + L1 * L2 * L3 * L4 * ((gL2 - gL3) * (L3 - L4) + 4 * (L2 - L3) * (gL3 - gL4)) )

        vec(215-k) = (L1 - L2)**2 * (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
     &              + L1 * L2 * L3 * L4 * (3 * (gL1 - gL2) * (L2 - L3) * (L3 - L4) + (L1 - L2) * (gL2 - gL3) * (L3 - L4) + (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )
        vec(216-k) = (L2 - L3)**2 * (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) * (L3 - L4) + 3 * (L1 - L2) * (gL2 - gL3) * (L3 - L4) + (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )
        vec(217-k) = (L3 - L4)**2 * (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) * (L3 - L4) + (L1 - L2) * (gL2 - gL3) * (L3 - L4) + 3 * (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )

        vec(218-k) = (L1 - L2)*(L2 - L3) * (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) * (L2 - L3) * (L3 - L4) + 2 * (L1 - L2) * (gL2 - gL3) * (L3 - L4) + (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )
        vec(219-k) = (L1 - L2)*(L3 - L4) * (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
     &              + L1 * L2 * L3 * L4 * (2 * (gL1 - gL2) * (L2 - L3) * (L3 - L4) + (L1 - L2) * (gL2 - gL3) * (L3 - L4) + 2 * (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )
        vec(220-k) = (L2 - L3)*(L3 - L4) * (                                                                                                 &
     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
     &              + L1 * L2 * L3 * L4 * ((gL1 - gL2) * (L2 - L3) * (L3 - L4) + 2 * (L1 - L2) * (gL2 - gL3) * (L3 - L4) + 2 * (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )

!
! General form of shape function gradients
!
!  Edge     u = p-2
!        vec(nnn-k) = (L1 - L2)**(u-1) * (
!                     (L1 - L2) * ( gL1 * L2 + L1 * gL2 )
!                    + L1 * L2 * u * (gL1 - gL2) )
!        vec(nnn-k) = (L2 - L3)**(u-1) * (
!                     (L2 - L3) * ( gL2 * L3 + L2 * gL3 )
!                    + L2 * L3 * u * (gL2 - gL3) )
!        vec(nnn-k) = (L1 - L3)**(u-1) * (
!                     (L1 - L3) * ( gL1 * L3 + L1 * gL3 )
!                    + L1 * L3 * u * (gL1 - gL3) )
!        vec(nnn-k) = (L1 - L4)**(u-1) * (
!                     (L1 - L4) * ( gL1 * L4 + L1 * gL4 )
!                    + L1 * L4 * u * (gL1 - gL4) )
!        vec(nnn-k) = (L2 - L4)**(u-1) * (
!                     (L2 - L4) * ( gL2 * L4 + L2 * gL4 )
!                    + L2 * L4 * u * (gL2 - gL4) )
!        vec(nnn-k) = (L3 - L4)**(u-1) * (
!                     (L3 - L4) * ( gL3 * L4 + L3 * gL4 )
!                    + L3 * L4 * u * (gL3 - gL4) )
!
!  Face     u + v = p-3
! Type      u = p-3, v=0
!        vec(nnn-k) = (L2 - L3)**(u-1) * (                                          &
!     &               (L2 - L3) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)      &
!     &              + L2 * L3 * L4 * (u * (gL2 - gL3) ) )
!        vec(nnn-k) = (L1 - L3)**(u-1) * (                                          &
!     &               (L1 - L3) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)      &
!     &              + L1 * L3 * L4 * (u * (gL1 - gL3) ) )
!        vec(nnn-k) = (L1 - L2)**(u-1) * (                                          &
!     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)      &
!     &              + L1 * L2 * L4 * (u * (gL1 - gL2) ) )
!        vec(nnn-k) = (L1 - L2)**(u-1) * (                                          &
!     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)      &
!     &              + L1 * L2 * L3 * (u * (gL1 - gL2) ) )
! Type      u>=1, v>=1
!        vec(nnn-k) = (L2 - L3)**(u-1)*(L3 - L4)**(v-1) * (                                          &
!     &               (L2 - L3) * (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)           &
!     &              + L2 * L3 * L4 * (u * (gL2 - gL3) * (L3 - L4) + v * (L2 - L3) * (gL3 - gL4) ) )
!        vec(nnn-k) = (L1 - L3)**(u-1)*(L3 - L4)**(v-1) * (                                          &
!     &               (L1 - L3) * (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)           &
!     &              + L1 * L3 * L4 * (u * (gL1 - gL3) * (L3 - L4) + v * (L1 - L3) * (gL3 - gL4) ) )
!        vec(nnn-k) = (L1 - L2)**(u-1)*(L2 - L4)**(v-1) * (                                          &
!     &               (L1 - L2) * (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)           &
!     &              + L1 * L2 * L4 * (u * (gL1 - gL2) * (L2 - L4) + v * (L1 - L2) * (gL2 - gL4) ) )
!        vec(nnn-k) = (L1 - L2)**(u-1)*(L2 - L3)**(v-1) * (                                          &
!     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)           &
!     &              + L1 * L2 * L3 * (u * (gL1 - gL2) * (L2 - L3) + v * (L1 - L2) * (gL2 - gL3) ) )
! Type      v = p-3, u=0
!        vec(nnn-k) = (L3 - L4)**(v-1) * (                                          &
!     &               (L3 - L4) * ((gL2 * L3 + L2 * gL3) * L4 + L2 * L3 * gL4)      &
!     &              + L2 * L3 * L4 * (v * (gL3 - gL4) ) )
!        vec(nnn-k) = (L3 - L4)**(v-1) * (                                          &
!     &               (L3 - L4) * ((gL1 * L3 + L1 * gL3) * L4 + L1 * L3 * gL4)      &
!     &              + L1 * L3 * L4 * (v * (gL3 - gL4) ) )
!        vec(nnn-k) = (L2 - L4)**(v-1) * (                                          &
!     &               (L2 - L4) * ((gL1 * L2 + L1 * gL2) * L4 + L1 * L2 * gL4)      &
!     &              + L1 * L2 * L4 * (v * (gL2 - gL4) ) )
!        vec(nnn-k) = (L2 - L3)**(v-1) * (                                          &
!     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 + L1 * L2 * gL3)      &
!     &              + L1 * L2 * L3 * (v * (gL2 - gL3) ) )
!
! Volume    u + v + w = p-4
!  type     u = p-4, w=0, v=0
!        vec(nnn-k) = (L1 - L2)**(u-1) * (                                                             &
!     &               (L1 - L2) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
!     &              + L1 * L2 * L3 * L4 * (u * (gL1 - gL2) ) )
!  type     v = p-4, w=0, u=0
!        vec(nnn-k) = (L2 - L3)**(v-1) * (                                                             &
!     &               (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
!     &              + L1 * L2 * L3 * L4 * (v * (gL2 - gL3) ) )
!  type     w = p-4, u=0, v=0
!        vec(nnn-k) = (L3 - L4)**(w-1) * (                                                             &
!     &               (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))  &
!     &              + L1 * L2 * L3 * L4 * (w * (gL3 - gL4)) )
!  type     u + v = p-4, w=0
!        vec(nnn-k) = (L1 - L2)**(u-1)*(L2 - L3)**(v-1) * (                                                       &
!     &               (L1 - L2) * (L2 - L3) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
!     &              + L1 * L2 * L3 * L4 * (u * (gL1 - gL2) * (L2 - L3) + v * (L1 - L2) * (gL2 - gL3) ) )
!  type     u + w = p-4, v=0
!        vec(nnn-k) = (L1 - L2)**(u-1)*(L3 - L4)**(w-1) * (                                                       &
!     &               (L1 - L2) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
!     &              + L1 * L2 * L3 * L4 * (u * (gL1 - gL2) * (L3 - L4) + w * (L1 - L2) * (gL3 - gL4)) )
!  type     v + w = p-4, u=0
!        vec(nnn-k) = (L2 - L3)**(v-1)*(L3 - L4)**(w-1) * (                                                       &
!     &               (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4)) &
!     &              + L1 * L2 * L3 * L4 * (v * (gL2 - gL3) * (L3 - L4) + w * (L2 - L3) * (gL3 - gL4)) )
!  type     u, v, w >= 1
!        vec(nnn-k) = (L1 - L2)**(u-1)*(L2 - L3)**(v-1)*(L3 - L4)**(w-1) * (                                                                                                 &
!     &               (L1 - L2) * (L2 - L3) * (L3 - L4) * ((gL1 * L2 + L1 * gL2) * L3 * L4 + L1 * L2 * (gL3 * L4 + L3 * gL4))                                                &
!     &              + L1 * L2 * L3 * L4 * (u * (gL1 - gL2) * (L2 - L3) * (L3 - L4) + v * (L1 - L2) * (gL2 - gL3) * (L3 - L4) + w * (L1 - L2) * (L2 - L3) * (gL3 - gL4)) )
!
!
         end select poly
     p = p + 1
    end do
!
      end subroutine gradshapesc3D
