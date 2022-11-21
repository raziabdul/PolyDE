      subroutine userbc3D( face, nature, bcnum, xyzs, pval, qval )
      use femtypes
      use feminterface,      only: fetchmatparameters
      use feminterface3D,    only: findelement3D
      use globalvariables3D, only: Elch, Kboltz, eps0, mu0, pi, pvalue, nod, vn, vv, numv, dommat, dom
      use matconstants
      implicit none
      integer (I4B) :: face, nature, bcnum
      real (DP) :: xyzs(3)
      complex(DPC) :: pval
      complex(DPC), optional :: qval(:)
      intent (in)  :: face, nature, bcnum, xyzs
      intent (out) :: pval, qval
!
!    $Revision: 1.9 $
!    $Date: 2014/02/27 15:06:55 $
!    $Author: m_kasper $
!
!  Boundary Condition for the case of scalar field
!
!  User programming interface for ANY Dirichlet and general of Neumann boundary conditions
!  Return values for coordinates xs, ys, zs.
!
!  Input:
!        bcnum          case selector from drawing layer text
!        xyzs           coordinates for determining Dirichlet
!
!  Output:
!        pval           scalar value for either the Dirichlet bc (u=p1) or General Neumann  p + sum( q(i)*u(i)., i=1..nnat)
!                       
!        qval(i)        vector for the General Neumann bc
!
!  local variables:
      integer (I4B)          :: elem
      real(DP)               :: T_ref, N_D, N_A, n_i
      real(DP),  allocatable :: list(:)
      logical                :: found
!__
! 
      select case (bcnum)
        case (11) !Semiconductor Mode, poisson equation, dirichlet bc
          if (nature.eq.1) then
            allocate(list(numparam))
            call findelement3D(xyzs,nod,vn,vv,numv,elem,found)
            call fetchmatparameters(list,dommat(dom(elem)),xyzs,elem)
          
            N_A   = list(4)
            N_D   = list(5)
            T_ref = list(6)
            n_i   = list(9)

            pval = pvalue(face,nature) + ((Kboltz*T_ref)/Elch)*asinh((N_D-N_A)/(2*n_i))
            
            deallocate(list)
          else
            pval = pvalue(face,nature)
          end if
        
        case (211)
          pval = pvalue(face,nature)
          qval = 0._DPC
        
      case default
!  if user bcnum is read, but no values is given, set x,y,z components values to 0
        print*,'Warning: userbc',bcnum,' is not defined, will use homogeneous case'
        pval = 0._DPC
        if ( present( qval ) ) qval = 0._DP
      end select
!
      return
      end subroutine userbc3D



      subroutine userbc3D_vec( bcnum, xyzs, pval, qval )
      use femtypes
      use globalvariables3D, only: eps0, mu0, omega, pi
      implicit none
      integer (I4B) :: bcnum
      real (DP) :: xyzs(3)
      complex(DPC) :: pval(3)
      complex(DPC), optional :: qval(3,3)
      intent (in) :: bcnum, xyzs
      intent (out) :: pval, qval
!
!  Boundary Condition for the case of vector valued field
!
!  User programming interface for ANY Dirichlet and general of Neumann boundary conditions
!  Return values for coordinates xs, ys, zs.
!
!  Input:
!        bcnum            case selector from drawing layer text
!        xyzs           coordinates for determining Dirichlet
!
!  Output:
!        pval(3)        vector values for either the Dirichlet bc (p1) or General Neumann
!                       (p2) with x, y, z components:  1=x, 2=y, 3=z
!        qval(3,3)      tensor values for the General Neumann bc.
!
!  local variables:
!      real (DP) s, startp(2), endp(2)
!      real (DP) dir(2), length, dirnormal(2)
      real (DP) k0, q0
!  Compute vacuum wavenumber
      k0 = omega*sqrt(mu0*eps0)
!  The factor in ABC for vacuum
      q0 = omega*sqrt(eps0/mu0)

!
      select case (bcnum)
      case (1)
!  half sine boundary functions in the x and y directions
! assume:
! -  1 x 1 face (unit square)
! -  the face to excite is in x-y plane at z=0
! -  the tangential field component is mathematically positive looking from 
!    OUTSIDE the cube
! 
         pval(1) =  sin( 1._DP*pi*xyzs(2) )
         pval(2) = -sin( 1._DP*pi*xyzs(1) )
! we don't care about the normal component of the face (points to +z) 
! so just let it 0
         pval(3) = 0._DP

      case (2)
!  full  sine boundary functions in the x and y directions
! assume:
! -  2 x 2 face (unit square)
! -  the face to excite is in x-y plane at z=0
! -  the tangential field component is mathematically positive looking from 
!    OUTSIDE the cube
! 
         pval(1) =  sin( 1._DP*pi*xyzs(2) )
         pval(2) = -sin( 1._DP*pi*xyzs(1) )
         pval(3) = 0._DP


      case default
!  if user bcnum is read, but no values is given, set x,y,z components values to 0
        print*,'Warning: userbc',bcnum,' is not defined, will use homogeneous case'
        pval = 0._DPC
        if ( present( qval ) ) qval = 0._DP
      end select
!
      return
      end subroutine userbc3D_vec