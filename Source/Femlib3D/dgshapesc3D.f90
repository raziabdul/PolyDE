      subroutine dgshapesc3D(l, gl, nu, polylo, polyhi, vec, errcode)
      use femtypes
      use globalvariables3d, only: polymaxsc
      implicit none
      integer (I4B) polylo, polyhi, errcode
      real (DP) l(4), gl(3,4)
      complex (DPC) nu(:,:), vec(:)
      intent (in) :: l, gl, polylo, polyhi
      intent (out) :: vec, errcode
!-------------------------------------------------------------------------------
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
!    $Revision: 1.5 $
!    $Date: 2015/11/11 17:39:04 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!  Evaluate the Laplacian of the 3D scalar shape functions for the corresponding polynomial degree
!
!  NOTES: HAS NOT BEEN USED AS OF 3.11.11
!  Input:
!            l        lambda, i.e., natural (barycentric) coordinates of the point 
!            gl       gradient of volume coordinates where:
!                        row indices denote  x,y,z components
!                        column indices denote lambdas (1,2,3,4) 
!            polylo   lowest polynomial degree to consider
!            polyhi   highest polyaomial degree 
!  Output:
!            vec      array of DIV GRAD of shape functions 
!  Errorcodes:
!            errcode  2001 = polylo is greater than polyhi
!                     2002 = polynomial degree cannot be higher than polymax. 
!                            See globalvariables
!
!--------------------------------------------------------------------------
! LATEST SHAPE FUNCTIONS:
!
! Warburton, no simplify
! 
! 
!--------------------------------------------------------------------------
!  local variables
      integer (I4B) p,k
      real (DP) L1, L2, L3, L4
      real (DP) gL11, gL12, gL13, & 
&               gL21, gL22, gL23, &
&               gL31, gL32, gL33, &
&               gL41, gL42, gL43
      complex (DPC) vectemp(3)
      complex (DPC) n11, n12, n13, n14, &
&                    n22, n23, n24, &
&                         n33, n34, &
&                              n44 

! gLmn : gradient of lambda m in n coordinate
! n_ij: shortform for eta ij, scalar product of gLi * gLj
!
!  Check values of polylo and polyhi
      if ((polylo .gt. polymaxsc) .or. (polyhi .gt. polymaxsc)) then
        errcode = 2002
        return
      end if
      if (polylo .gt. polyhi) then
        errcode = 2001
        return
      end if
!
!  Calculation of the div of the gradients (Laplacian)

!  intermediate variables to avoid arrays of l and gl
      L1 = l(1)
      L2 = l(2)
      L3 = l(3)
      L4 = l(4)
!------------------------------------------------------------------------------
!  Substitute n in vec in terms of nu and grad of lambda's
!
!  n_ij = 1/2[ (nu*gl_i)*gl_j + (nu*gl_j)*gl_i ]
!
! expanding the right hand side,
!
! n_ij= nu_11*(delL_i/delx)(delL_j/delx) + nu_12*(delL_i/dely)(delL_j/delx) + nu_13*(delL_i/delz)(delL_j/delx) +  
!       nu_21*(delL_i/delx)(delL_j/dely) + nu_22*(delL_i/dely)(delL_j/dely) + nu_23*(delL_i/delz)(delL_j/dely) +  
!       nu_31*(delL_i/delx)(delL_j/delz) + nu_32*(delL_i/dely)(delL_j/delz) + nu_33*(delL_i/delz)(delL_j/delz)
!
! But here we use the intrinsic funcitons matmul and dot_product
!
! n_ij = transpose(gl_i) * nu * gl_j
!
!-----------------------------------------------------------------------------
!
! intermediate variables to avoid gl(3,4) arrays
! use explicit substitutions for gL: 
! gLmn = gL(m,n), m = Volume Coord, n = Cartesian: (x,y,z) = (1,2,3)
! WE MAY NOT NEED THIS, WILL VERIFY LATER
 
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
! substitute n: n_ij = grad lambda^T * nu * grad lambda

  vectemp = matmul(nu,gl(:,1))
  n11=dot_product(gl(:,1),vectemp)

  vectemp = matmul(nu,gl(:,2))
  n12=dot_product(gl(:,1),vectemp)
  n22=dot_product(gl(:,2),vectemp)

  vectemp = matmul(nu,gl(:,3))
  n13=dot_product(gl(:,1),vectemp)
  n23=dot_product(gl(:,2),vectemp)
  n33=dot_product(gl(:,3),vectemp)

  vectemp = matmul(nu,gl(:,4))
  n14=dot_product(gl(:,1),vectemp)
  n24=dot_product(gl(:,2),vectemp)
  n34=dot_product(gl(:,3),vectemp)
  n44=dot_product(gl(:,4),vectemp)



!!$  n11=nu(1,1)*gL11*gL11 + nu(1,2)*gL12*gL11 + nu(1,3)*gL13*gL11 +
!!$      nu(2,1)*gL11*gL12 + nu(2,2)*gL12*gL12 + nu(2,3)*gL13*gL13 +
!!$      nu(3,1)*gL11*gL13 + nu(3,2)*gL12*gL13 + nu(3,3)*gL13*gL13
!!$  n12=dot_product(gl(:,1),gl(:,2)) 
!!$  n13=dot_product(gl(:,1),gl(:,3))
!!$  n14=dot_product(gl(:,1),gl(:,4))
!!$  n22=dot_product(gl(:,2),gl(:,2))
!!$  n23=dot_product(gl(:,2),gl(:,3))
!!$  n24=dot_product(gl(:,2),gl(:,4))
!!$  n33=dot_product(gl(:,3),gl(:,3))
!!$  n34=dot_product(gl(:,3),gl(:,4))
!!$  n44=dot_product(gl(:,4),gl(:,4))

!  The polylo and polyhi variables are used to select the shape functions. If polylo .eq.
!  polyhi, only the new shape functions of the corresponding polynomial degree p are cal-
!  culated. That means: if polylo .eq. polyhi = 2 ==> xsi(5) to xsi(10) are calculated.
!  If polylo .neq. polyhi, shape functions starting from polylo to polyhi are calculated.
!  That means: if polylo = 2 and polyhi = 3 ==> xsi(5) up to xsi(20) are calculated.
!

      p = polylo

! k depends on the formulations of shape functions used.
! This is the standard C(0) type element (continuous)
      if (polylo .gt. 1) then
        k = ( polylo+1 )*( polylo+2 )*( polylo+3 )/6
      else
        k = 0
      end if

!
!-------------------------------------------------------------------------------!
!  BEGIN LOOP TO CALCULATE LAPLACIAN OF SHAPE FUNCTIONS FOR POLYNOMIAL DEGREE p !
!-------------------------------------------------------------------------------!
      do
        if (p .gt. polyhi) exit    ! stopping criteria
!  Begin "case" control structure to select the appropriate shape functions for p
        poly: select case (p)
!
!  TO Do the output does not agre to the difinition of shape functions
!
case (1)  ! Laplacian of shape functions for p = 1
        vec(  1-k) = 0_DP
        vec(  2-k) = 0_DP
        vec(  3-k) = 0_DP
        vec(  4-k) = 0_DP

case (2)  ! Laplacian of shape functions for p = 2
        vec(  5-k) = 2.0_DP*n12 
        vec(  6-k) = 2.0_DP*n23 
        vec(  7-k) = 2.0_DP*n13 
        vec(  8-k) = 2.0_DP*n14 
        vec(  9-k) = 2.0_DP*n24 
        vec( 10-k) = 2.0_DP*n34 

case (3)  ! Laplacian of shape functions for p = 3
!  edge DOF
        vec( 11-k) = n11 * ( 2*L2) + &
     &               n22 * (-2*L1) + &
     &               n12 * ( 2*(L1-L2))
        vec( 12-k) = n22 * ( 2*L3) + &
     &               n33 * (-2*L2) + &
     &               n23 * ( 2*(L2-L3))
        vec( 13-k) = n11 * ( 2*L3) + &
     &               n33 * (-2*L1) + &
     &               n13 * ( 2*(L1-L3))
        vec( 14-k) = n11 * ( 2*L4) + &
     &               n44 * (-2*L1) + &
     &               n14 * ( 2*(L1-L4))
        vec( 15-k) = n22 * ( 2*L4) + &
     &               n44 * (-2*L2) + &
     &               n24 * ( 2*(L2-L4))
        vec( 16-k) = n33 * ( 2*L4) + &
     &               n44 * (-2*L3) + &
     &               n34 * ( 2*(L3-L4))
!  face DOF
        vec( 17-k) = n24 * ( L3) + &
     &               n23 * ( L4) + &
     &               n34 * ( L2)
        vec( 18-k) = n13 * ( L4) + &
     &               n14 * ( L3) + &
     &               n34 * ( L1)
        vec( 19-k) = n14 * ( L2) + &
     &               n12 * ( L4) + &
     &               n24 * ( L1)
        vec( 20-k) = n12 * ( L3) + &
     &               n13 * ( L2) + &
     &               n23 * ( L1)

case (4)  ! Laplacian of shape functions for p = 4
        vec( 21-k) = n11 * (2*L2 * (3*L1 - 2*L2)) + &
     &               n22 * (2*L1 * (3*L2 - 2*L1)) + &
     &               n12 * (3*(L1-L2)**2 - 2*L1*L2)
        vec( 22-k) = n22 * (2*L3 * (3*L2 - 2*L3)) + &
     &               n33 * (2*L2 * (3*L3 - 2*L2)) + &
     &               n23 * (3*(L2-L3)**2 - 2*L2*L3)
        vec( 23-k) = n11 * (2*L3 * (3*L1 - 2*L3)) + &
     &               n33 * (2*L1 * (3*L3 - 2*L1)) + &
     &               n13 * (3*(L1-L3)**2 - 2*L1*L3)
        vec( 24-k) = n11 * (2*L4 * (3*L1 - 2*L4)) + &
     &               n44 * (2*L1 * (3*L4 - 2*L1)) + &
     &               n14 * (3*(L1-L4)**2 - 2*L1*L4)
        vec( 25-k) = n22 * (2*L4 * (3*L2 - 2*L4)) + &
     &               n44 * (2*L2 * (3*L4 - 2*L2)) + &
     &               n24 * (3*(L2-L4)**2 - 2*L2*L4)
        vec( 26-k) = n33 * (2*L4 * (3*L3 - 2*L4)) + &
     &               n44 * (2*L3 * (3*L4 - 2*L3)) + &
     &               n34 * (3*(L3-L4)**2 - 2*L3*L4)
!  face DOF
        vec( 27-k) = n22 * ( 2*L3*L4) +         &
     &               n33 * (-2*L2*L4) +         &
     &               n23 * (2*L4*(  L2-  L3)) + &
     &               n24 * (  L3*(2*L2-  L3)) + &
     &               n34 * (  L2*(  L2-2*L3))
        vec( 29-k) = n11 * ( 2*L3*L4) +         &
     &               n33 * (-2*L1*L4) +         &
     &               n13 * (2*L4*(  L1-  L3)) + &
     &               n14 * (  L3*(2*L1-  L3)) + &
     &               n34 * (  L1*(  L1-2*L3))
        vec( 31-k) = n11 * ( 2*L2*L4) +         &
     &               n22 * (-2*L1*L4) +         &
     &               n12 * (2*L4*(  L1-  L2)) + &
     &               n14 * (  L2*(2*L1-  L2)) + &
     &               n24 * (  L1*(  L1-2*L2))
        vec( 33-k) = n11 * ( 2*L2*L3) +         &
     &               n22 * (-2*L1*L3) +         &
     &               n12 * (2*L3*(  L1-  L2)) + &
     &               n13 * (  L2*(2*L1-  L2)) + &
     &               n23 * (  L1*(  L1-2*L2))
        vec( 28-k) = n33 * ( 2*L2*L4) +         &
     &               n44 * (-2*L2*L3) +         &
     &               n23 * (  L4*(2*L3-  L4)) + &
     &               n24 * (  L3*(  L3-2*L4)) + &
     &               n34 * (2*L2*(  L3-  L4))
        vec( 30-k) = n33 * ( 2*L1*L4) +         &
     &               n44 * (-2*L1*L3) +         &
     &               n13 * (  L4*(2*L3-  L4)) + &
     &               n14 * (  L3*(  L3-2*L4)) + &
     &               n34 * (2*L1*(  L3-  L4))
        vec( 32-k) = n22 * ( 2*L1*L4) +         &
     &               n44 * (-2*L1*L2) +         &
     &               n12 * (  L4*(2*L2-  L4)) + &
     &               n14 * (  L2*(  L2-2*L4)) + &
     &               n24 * (2*L1*(  L2-  L4))
        vec( 34-k) = n22 * ( 2*L1*L3) +         &
     &               n33 * (-2*L1*L2) +         &
     &               n12 * (  L3*(2*L2-  L3)) + &
     &               n13 * (  L2*(  L2-2*L3)) + &
     &               n23 * (2*L1*(  L2-  L3))
!  volume DOF
!        vec( 34-k) = ?

case (5)  ! Laplacian of shape functions for p = 5
        vec( 36-k) = (n11 * (3*L2 * (4*L1 - 2*L2)) + &
     &                n22 * (3*L1 * (4*L2 - 2*L1)) + &
     &                n12 * (4*(L1-L2)**2 - 6*L1*L2) ) * (L1 - L2)
        vec( 37-k) = (n22 * (3*L3 * (4*L2 - 2*L3)) + &
     &                n33 * (3*L2 * (4*L3 - 2*L2)) + &
     &                n23 * (4*(L2-L3)**2 - 6*L2*L3) ) * (L2 - L3)
        vec( 38-k) = (n11 * (3*L3 * (4*L1 - 2*L3)) + &
     &                n33 * (3*L1 * (4*L3 - 2*L1)) + &
     &                n13 * (4*(L1-L3)**2 - 6*L1*L3) ) * (L1 - L3)
        vec( 39-k) = (n11 * (3*L4 * (4*L1 - 2*L4)) + &
     &                n44 * (3*L1 * (4*L4 - 2*L1)) + &
     &                n14 * (4*(L1-L4)**2 - 6*L1*L4) ) * (L1 - L4)
        vec( 40-k) = (n22 * (3*L4 * (4*L2 - 2*L4)) + &
     &                n44 * (3*L2 * (4*L4 - 2*L2)) + &
     &                n24 * (4*(L2-L4)**2 - 6*L2*L4) ) * (L2 - L4)
        vec( 41-k) = (n33 * (3*L4 * (4*L3 - 2*L4)) + &
     &                n44 * (3*L3 * (4*L4 - 2*L3)) + &
     &                n34 * (4*(L3-L4)**2 - 6*L3*L4) ) * (L3 - L4)
!  face DOF
        vec( 42-k) = n22 * 2*L3*L4*( 3*L2-2*L3) +        &
     &               n33 * 2*L2*L4*( 3*L3-2*L2) +        &
     &               n23 * L4*(3*(L2-L3)**2 - 2*L2*L3) + &
     &               n24 * L3*(3*L2-  L3)*(L2-L3) +      &
     &               n34 * L2*(  L2-3*L3)*(L2-L3)
        vec( 45-k) = n11 * 2*L3*L4*( 3*L1-2*L3) +        &
     &               n33 * 2*L1*L4*( 3*L3-2*L1) +        &
     &               n13 * L4*(3*(L1-L3)**2 - 2*L1*L3) + &
     &               n14 * L3*(3*L1-  L3)*(L1-L3) +      &
     &               n34 * L1*(  L1-3*L3)*(L1-L3)
        vec( 48-k) = n11 * 2*L2*L4*( 3*L1-2*L2) +        &
     &               n22 * 2*L1*L4*( 3*L2-2*L1) +        &
     &               n12 * L4*(3*(L1-L2)**2 - 2*L1*L2) + &
     &               n14 * L2*(3*L1-  L2)*(L1-L2) +      &
     &               n24 * L1*(  L1-3*L2)*(L1-L2)
        vec( 51-k) = n11 * 2*L2*L3*( 3*L1-2*L2) +        &
     &               n22 * 2*L1*L3*( 3*L2-2*L1) +        &
     &               n12 * L3*(3*(L1-L2)**2 - 2*L1*L2) + &
     &               n13 * L2*(3*L1-  L2)*(L1-L2) +      &
     &               n23 * L1*(  L1-3*L2)*(L1-L2)
        vec( 43-k) = n22 * 2*L3*L4*(   L3-  L4) +                  &
     &               n33 * 2*L2*L4*(   L2-3*L3+  L4) +             &
     &               n44 * 2*L2*L3*(   L3-  L2) +                  &
     &               n23 * L4*(4*L2*L3 - 2*L4*(L2-L3) -3*L3**2 ) + &
     &               n24 * L3*(  L3-2*L4)*(2*L2-  L3) +            &
     &               n34 * L2*(4*L3*L4 + 2*L2*(L3-L4) -3*L3**2 )
        vec( 46-k) = n11 * 2*L3*L4*(   L3-  L4) +                  &
     &               n33 * 2*L1*L4*(   L1-3*L3+  L4) +             &
     &               n44 * 2*L1*L3*(   L3-  L1) +                  &
     &               n13 * L4*(4*L1*L3 - 2*L4*(L1-L3) -3*L3**2 ) + &
     &               n14 * L3*(  L3-2*L4)*(2*L1-  L3) +            &
     &               n34 * L1*(4*L3*L4 + 2*L1*(L3-L4) -3*L3**2 )
        vec( 49-k) = n11 * 2*L2*L4*(   L2-  L4) +                  &
     &               n22 * 2*L1*L4*(   L1-3*L2+  L4) +             &
     &               n44 * 2*L1*L2*(   L2-  L1) +                  &
     &               n12 * L4*(4*L1*L2 - 2*L4*(L1-L2) -3*L2**2 ) + &
     &               n14 * L3*(  L2-2*L4)*(2*L1-  L2) +            &
     &               n24 * L1*(4*L3*L4 + 2*L1*(L2-L4) -3*L2**2 )
        vec( 52-k) = n11 * 2*L2*L3*(   L2-  L3) +                  &
     &               n22 * 2*L1*L4*(   L1-3*L2+  L3) +             &
     &               n33 * 2*L1*L2*(   L2-  L1) +                  &
     &               n12 * L3*(4*L1*L2 - 2*L3*(L1-L2) -3*L2**2 ) + &
     &               n13 * L3*(  L2-2*L3)*(2*L1-  L2) +            &
     &               n23 * L1*(4*L3*L3 + 2*L1*(L2-L3) -3*L2**2 )
        vec( 44-k) = n33 * 2*L2*L4*( 3*L3-2*L4) +        &
     &               n44 * 2*L2*L3*( 3*L4-2*L3) +        &
     &               n23 * L4*(3*L3-  L4)*(L3-L4) +      &
     &               n24 * L3*(  L3-3*L4)*(L3-L4) +      &
     &               n34 * L2*(3*(L3-L4)**2 - 2*L3*L4)
        vec( 47-k) = n33 * 2*L1*L4*( 3*L3-2*L4) +        &
     &               n44 * 2*L1*L3*( 3*L4-2*L3) +        &
     &               n13 * L4*(3*L3-  L4)*(L3-L4) +      &
     &               n14 * L3*(  L3-3*L4)*(L3-L4) +      &
     &               n34 * L1*(3*(L3-L4)**2 - 2*L3*L4)
        vec( 50-k) = n22 * 2*L1*L4*( 3*L2-2*L4) +        &
     &               n44 * 2*L1*L2*( 3*L4-2*L2) +        &
     &               n12 * L4*(3*L2-  L4)*(L2-L4) +      &
     &               n14 * L2*(  L2-3*L4)*(L2-L4) +      &
     &               n24 * L1*(3*(L2-L4)**2 - 2*L2*L4)
        vec( 53-k) = n22 * 2*L1*L3*( 3*L2-2*L3) +        &
     &               n33 * 2*L1*L2*( 3*L3-2*L2) +        &
     &               n12 * L3*(3*L2-  L3)*(L2-L3) +      &
     &               n13 * L2*(  L2-3*L3)*(L2-L3) +      &
     &               n23 * L1*(3*(L2-L3)**2 - 2*L2*L3)
!  volume DOF
!        vec( 54-k) = ?
!        vec( 55-k) = ?
!        vec( 56-k) = ?


case (6)  ! Laplacian of shape functions for p = 6
        vec( 57-k) = (n11 * (4*L2 * (5*L1 - 2*L2)) + &
     &                n22 * (4*L1 * (5*L2 - 2*L1)) + &
     &                n12 * (5*(L1-L2)**2 -12*L1*L2) ) * (L1 - L2)**2
        vec( 58-k) = (n22 * (4*L3 * (5*L2 - 2*L3)) + &
     &                n33 * (4*L2 * (5*L3 - 2*L2)) + &
     &                n23 * (5*(L2-L3)**2 -12*L2*L3) ) * (L2 - L3)**2
        vec( 59-k) = (n11 * (4*L3 * (5*L1 - 2*L3)) + &
     &                n33 * (4*L1 * (5*L3 - 2*L1)) + &
     &                n13 * (5*(L1-L3)**2 -12*L1*L3) ) * (L1 - L3)**2
        vec( 60-k) = (n11 * (4*L4 * (5*L1 - 2*L4)) + &
     &                n44 * (4*L1 * (5*L4 - 2*L1)) + &
     &                n14 * (5*(L1-L4)**2 -12*L1*L4) ) * (L1 - L4)**2
        vec( 61-k) = (n22 * (4*L4 * (5*L2 - 2*L4)) + &
     &                n44 * (4*L2 * (5*L4 - 2*L2)) + &
     &                n24 * (5*(L2-L4)**2 -12*L2*L4) ) * (L2 - L4)**2
        vec( 62-k) = (n33 * (4*L4 * (5*L3 - 2*L4)) + &
     &                n44 * (4*L3 * (5*L4 - 2*L3)) + &
     &                n34 * (5*(L3-L4)**2 -12*L3*L4) ) * (L3 - L4)**2
!  face DOF
        vec( 63-k) = (n22 * 6*L3*L4*(2*L2-  L3) +                  &
     &                n33 * 6*L2*L4*(2*L3-  L2) +                  &
     &                n23 * 2*L4   *((L2-L3)**2-3*L2*L3)+          &
     &                n24 *   L3   *(  L2-  L3)*(4*L2-  L3) +      &
     &                n34 *   L2   *(  L2-  L3)*(  L2-4*L3) ) * (L2 - L3)
        vec( 67-k) = (n11 * 6*L3*L4*(2*L1-  L3) +                  &
     &                n33 * 6*L1*L4*(2*L3-  L1) +                  &
     &                n13 * 2*L4   *((L1-L3)**2-3*L1*L3)+          &
     &                n14 *   L3   *(  L1-  L3)*(4*L1-  L3) +      &
     &                n34 *   L1   *(  L1-  L3)*(  L1-4*L3) ) * (L1 - L3)
        vec( 71-k) = (n11 * 6*L2*L4*(2*L1-  L2) +                  &
     &                n22 * 6*L1*L4*(2*L2-  L1) +                  &
     &                n12 * 2*L4   *((L1-L2)**2-3*L1*L2)+          &
     &                n14 *   L2   *(  L1-  L2)*(4*L1-  L2) +      &
     &                n24 *   L1   *(  L1-  L2)*(  L1-4*L2) ) * (L1 - L2)
        vec( 75-k) = (n11 * 6*L2*L3*(2*L1-  L2) +                  &
     &                n22 * 6*L1*L3*(2*L2-  L1) +                  &
     &                n12 * 2*L3   *((L1-L2)**2-3*L1*L2)+          &
     &                n13 *   L2   *(  L1-  L2)*(4*L1-  L2) +      &
     &                n23 *   L1   *(  L1-  L2)*(  L1-4*L2) ) * (L1 - L2)
        vec( 64-k) = n22 * 2*L3*L4*(L3-L4)*(3*L2-2*L3) +                                                    &
     &               n33 * 2*L2*L4*(L2**2 - 6*L2*L3 + 6*L3**2 + 2*L2*L4 - 3*L3*L4) +                        &
     &               n44 * (-2)*L2*L3*(L2-L3)**2 +                                                          &
     &               n23 * L4*(6*L2**2*L3 - 3*L2**2*L4 - 12*L2*L3**2 + 8*L2*L3*L4 + 4*L3**3 - 3*L3**2*L4) + &
     &               n24 * L3*(L2-L3)*(L3-2*L4)*(3*L2-L3) +                                                 &
     &               n34 * 2*L2*(L2-L3)*(L2*(L3-L4) - 2*L3**2 + 3*L3*L4)
        vec( 68-k) = n11 * 2*L3*L4*(L3-L4)*(3*L1-2*L3) +                                                    &
     &               n33 * 2*L1*L4*(L1**2 - 6*L1*L3 + 6*L3**2 + 2*L1*L4 - 3*L3*L4) +                        &
     &               n44 * (-2)*L1*L3*(L1-L3)**2 +                                                          &
     &               n13 * L4*(6*L1**2*L3 - 3*L1**2*L4 - 12*L1*L3**2 + 8*L1*L3*L4 + 4*L3**3 - 3*L3**2*L4) + &
     &               n14 * L3*(L1-L3)*(L3-2*L4)*(3*L1-L3) +                                                 &
     &               n34 * 2*L1*(L1-L3)*(L1*(L3-L4) - 2*L3**2 + 3*L3*L4)
        vec( 72-k) = n11 * 2*L2*L4*(L2-L4)*(3*L1-2*L2) +                                                    &
     &               n22 * 2*L1*L4*(L1**2 - 6*L1*L2 + 6*L2**2 + 2*L1*L4 - 3*L2*L4) +                        &
     &               n44 * (-2)*L1*L2*(L1-L2)**2 +                                                          &
     &               n12 * L4*(6*L1**2*L2 - 3*L1**2*L4 - 12*L1*L2**2 + 8*L1*L2*L4 + 4*L2**3 - 3*L2**2*L4) + &
     &               n14 * L2*(L1-L2)*(L2-2*L4)*(3*L1-L2) +                                                 &
     &               n24 * 2*L1*(L1-L2)*(L1*(L2-L4) - 2*L2**2 + 3*L2*L4)
        vec( 76-k) = n11 * 2*L2*L3*(L2-L3)*(3*L1-2*L2) +                                                    &
     &               n22 * 2*L1*L3*(L1**2 - 6*L1*L2 + 6*L2**2 + 2*L1*L3 - 3*L2*L3) +                        &
     &               n33 * (-2)*L1*L2*(L1-L2)**2 +                                                          &
     &               n12 * L3*(6*L1**2*L2 - 3*L1**2*L3 - 12*L1*L2**2 + 8*L1*L2*L3 + 4*L2**3 - 3*L2**2*L3) + &
     &               n13 * L2*(L1-L2)*(L2-2*L3)*(3*L1-L2) +                                                 &
     &               n23 * 2*L1*(L1-L2)*(L1*(L2-L3) - 2*L2**2 + 3*L2*L3)
        vec( 65-k) = n22 * 2*L3*L4*(L3-L4)**2 +                                                             &
     &               n33 * 2*L2*L4*(3*L2*L3 - 2*L2*L4 - 6*L3**2 + 6*L3*L4 - L4**2) +                        &
     &               n44 * (-2)*L2*L3*(L2-L3)*(2*L3-3*L4) +                                                 &
     &               n23 * 2*L4*(L3-L4)*(3*L2*L3 - L2*L4 - 2*L3**2 + L3*L4) +                               &
     &               n24 * L3*(L3-L4)*(L3-3*L4)*(-L3+2*L2) +                                                &
     &               n34 * L2*(3*L2*L3**2 - 8*L2*L3*L4 + 3*L2*L4**2 - 4*L3**3 + 12*L3**2*L4 - 6*L3*L4*2)
        vec( 69-k) = n11 * 2*L3*L4*(L3-L4)**2 +                                                             &
     &               n33 * 2*L2*L4*(3*L1*L3 - 2*L1*L4 - 6*L3**2 + 6*L3*L4 - L4**2) +                        &
     &               n44 * (-2)*L1*L3*(L1-L3)*(2*L3-3*L4) +                                                 &
     &               n13 * 2*L4*(L3-L4)*(3*L1*L3 - L1*L4 - 2*L3**2 + L3*L4) +                               &
     &               n14 * L3*(L3-L4)*(L3-3*L4)*(-L3+2*L1) +                                                &
     &               n34 * L1*(3*L1*L3**2 - 8*L1*L3*L4 + 3*L1*L4**2 - 4*L3**3 + 12*L3**2*L4 - 6*L3*L4*2)
        vec( 73-k) = n11 * 2*L2*L4*(L2-L4)**2 +                                                             &
     &               n22 * 2*L2*L4*(3*L1*L2 - 2*L1*L4 - 6*L2**2 + 6*L2*L4 - L4**2) +                        &
     &               n44 * (-2)*L1*L2*(L1-L2)*(2*L2-3*L4) +                                                 &
     &               n12 * 2*L4*(L2-L4)*(3*L1*L2 - L1*L4 - 2*L2**2 + L2*L4) +                               &
     &               n14 * L2*(L2-L4)*(L2-3*L4)*(-L2+2*L1) +                                                &
     &               n24 * L1*(3*L1*L2**2 - 8*L1*L2*L4 + 3*L1*L4**2 - 4*L2**3 + 12*L2**2*L4 - 6*L2*L4*2)
        vec( 77-k) = n11 * 2*L2*L3*(L2-L3)**2 +                                                             &
     &               n22 * 2*L2*L3*(3*L1*L2 - 2*L1*L3 - 6*L2**2 + 6*L2*L3 - L3**2) +                        &
     &               n33 * (-2)*L1*L2*(L1-L2)*(2*L2-3*L3) +                                                 &
     &               n12 * 2*L3*(L2-L3)*(3*L1*L2 - L1*L3 - 2*L2**2 + L2*L3) +                               &
     &               n13 * L2*(L2-L3)*(L2-3*L3)*(-L2+2*L1) +                                                &
     &               n23 * L1*(3*L1*L2**2 - 8*L1*L2*L3 + 3*L1*L3**2 - 4*L2**3 + 12*L2**2*L3 - 6*L2*L3*2)
        vec( 66-k) = (n33 * 6*L2*L4*(2*L3-L4) +                   &
     &                n44 * 6*L2*L3*(2*L4-L3) +                   &
     &                n23 * L4*(L3-L4)**(4*L3-L4)+                &
     &                n24 * L3*(L3-L4)**(L3-4*L4) +               &
     &                n34 * 2*L2*(2*L3**2-7*L3*L4+2*L4**2) ) * (L3 - L4)
        vec( 70-k) = (n33 * 6*L1*L4*(2*L3-L4) +                   &
     &                n44 * 6*L1*L3*(2*L4-L3) +                   &
     &                n13 * L4*(L3-L4)**(4*L3-L4)+                &
     &                n14 * L3*(L3-L4)**(L3-4*L4) +               &
     &                n34 * 2*L1*(2*L3**2 - 7*L3*L4 + 2*L4**2) ) * (L3 - L4)
        vec( 74-k) = (n22 * 6*L1*L4*(2*L2-L4) +                   &
     &                n44 * 6*L1*L2*(2*L4-L2) +                   &
     &                n12 * L4*(L2-L4)**(4*L2-L4)+                &
     &                n14 * L3*(L2-L4)**(L2-4*L4) +               &
     &                n24 * 2*L1*(2*L2**2 - 7*L2*L4 + 2*L4**2) ) * (L2 - L4)
        vec( 78-k) = (n22 * 6*L1*L3*(2*L2-L3) +                   &
     &                n33 * 6*L1*L2*(2*L3-L2) +                   &
     &                n12 * L3*(L2-L3)**(4*L2-L3)+                &
     &                n13 * L3*(L2-L3)**(L2-4*L3) +               &
     &                n23 * 2*L1*(2*L2**2 - 7*L2*L3 + 2*L3**2) ) * (L2 - L3)

!  volume DOF
!        vec( 79-k) = ?
!        vec( 80-k) = ?
!        vec( 81-k) = ?
!        vec( 82-k) = ?
!        vec( 83-k) = ?
!        vec( 84-k) = ?


case (7)  ! Laplacian of shape functions for p = 7
        vec( 85-k) = (n11 * (5*L2 * (6*L1 - 2*L2)) + &
     &                n22 * (5*L1 * (6*L2 - 2*L1)) + &
     &                n12 * (6*(L1-L2)**2 -20*L1*L2) ) * (L1 - L2)**3
        vec( 86-k) = (n22 * (5*L3 * (6*L2 - 2*L3)) + &
     &                n33 * (5*L2 * (6*L3 - 2*L2)) + &
     &                n23 * (6*(L2-L3)**2 -20*L2*L3) ) * (L2 - L3)**3
        vec( 87-k) = (n11 * (5*L3 * (6*L1 - 2*L3)) + &
     &                n33 * (5*L1 * (6*L3 - 2*L1)) + &
     &                n13 * (6*(L1-L3)**2 -20*L1*L3) ) * (L1 - L3)**3
        vec( 88-k) = (n11 * (5*L4 * (6*L1 - 2*L4)) + &
     &                n44 * (5*L1 * (6*L4 - 2*L1)) + &
     &                n14 * (6*(L1-L4)**2 -20*L1*L4) ) * (L1 - L4)**3
        vec( 89-k) = (n22 * (5*L4 * (6*L2 - 2*L4)) + &
     &                n44 * (5*L2 * (6*L4 - 2*L2)) + &
     &                n24 * (6*(L2-L4)**2 -20*L2*L4) ) * (L2 - L4)**3
        vec( 90-k) = (n33 * (5*L4 * (6*L3 - 2*L4)) + &
     &                n44 * (5*L3 * (6*L4 - 2*L3)) + &
     &                n34 * (6*(L3-L4)**2 -20*L3*L4) ) * (L3 - L4)**3
!  face DOF

!  volume DOF


case (8)  ! Laplacian of shape functions for p = 8
        vec(121-k) = (n11 * (6*L2 * (7*L1 - 2*L2)) + &
     &                n22 * (6*L1 * (7*L2 - 2*L1)) + &
     &                n12 * (7*(L1-L2)**2 -30*L1*L2) ) * (L1 - L2)**4
        vec(122-k) = (n22 * (6*L3 * (7*L2 - 2*L3)) + &
     &                n33 * (6*L2 * (7*L3 - 2*L2)) + &
     &                n23 * (7*(L2-L3)**2 -30*L2*L3) ) * (L2 - L3)**4
        vec(123-k) = (n11 * (6*L3 * (7*L1 - 2*L3)) + &
     &                n33 * (6*L1 * (7*L3 - 2*L1)) + &
     &                n13 * (7*(L1-L3)**2 -30*L1*L3) ) * (L1 - L3)**4
        vec(124-k) = (n11 * (6*L4 * (7*L1 - 2*L4)) + &
     &                n44 * (6*L1 * (7*L4 - 2*L1)) + &
     &                n14 * (7*(L1-L4)**2 -30*L1*L4) ) * (L1 - L4)**4
        vec(125-k) = (n22 * (6*L4 * (7*L2 - 2*L4)) + &
     &                n44 * (6*L2 * (7*L4 - 2*L2)) + &
     &                n24 * (7*(L2-L4)**2 -30*L2*L4) ) * (L2 - L4)**4
        vec(126-k) = (n33 * (6*L4 * (7*L3 - 2*L4)) + &
     &                n44 * (6*L3 * (7*L4 - 2*L3)) + &
     &                n34 * (7*(L3-L4)**2 -30*L3*L4) ) * (L3 - L4)**4
!  face DOF

!  volume DOF


case (9)  ! Laplacian of shape functions for p = 9
        vec(166-k) = (n11 * (7*L2 * (8*L1 - 2*L2)) + &
     &                n22 * (7*L1 * (8*L2 - 2*L1)) + &
     &                n12 * (8*(L1-L2)**2 -42*L1*L2) ) * (L1 - L2)**5
        vec(167-k) = (n22 * (7*L3 * (8*L2 - 2*L3)) + &
     &                n33 * (7*L2 * (8*L3 - 2*L2)) + &
     &                n23 * (8*(L2-L3)**2 -42*L2*L3) ) * (L2 - L3)**5
        vec(168-k) = (n11 * (7*L3 * (8*L1 - 2*L3)) + &
     &                n33 * (7*L1 * (8*L3 - 2*L1)) + &
     &                n13 * (8*(L1-L3)**2 -42*L1*L3) ) * (L1 - L3)**5
        vec(169-k) = (n11 * (7*L4 * (8*L1 - 2*L4)) + &
     &                n44 * (7*L1 * (8*L4 - 2*L1)) + &
     &                n14 * (8*(L1-L4)**2 -42*L1*L4) ) * (L1 - L4)**5
        vec(170-k) = (n22 * (7*L4 * (8*L2 - 2*L4)) + &
     &                n44 * (7*L2 * (8*L4 - 2*L2)) + &
     &                n24 * (8*(L2-L4)**2 -42*L2*L4) ) * (L2 - L4)**5
        vec(171-k) = (n33 * (7*L4 * (8*L3 - 2*L4)) + &
     &                n44 * (7*L3 * (8*L4 - 2*L3)) + &
     &                n34 * (8*(L3-L4)**2 -42*L3*L4) ) * (L3 - L4)**5
!  face DOF

!  volume DOF


        end select poly
!
        p=p+1
      end do

!
    end subroutine dgshapesc3D
