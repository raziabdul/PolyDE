      subroutine field3D(elem, lambda, u, curlu)
      use feminterface3D, only: shapefunctionv3D, pdecoeff3D, lam2xyz
      use feminterface3D, only: shapefunctionsc3D

      use globalvariables3D, only : eltype, nnat, nod, vgdof, vn, vp, x
      use femtypes
      implicit none
      integer (I4B) elem
      real (DP) lambda(4)
      complex (DPC) u(3,nnat), curlu(3,nnat)
      intent (in) ::  elem, lambda
      intent (out) ::  u, curlu
!------------------------------------------------------------------------------
!
!    $Revision: 1.6 $
!    $Date: 2014/05/26 10:50:59 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
!
!  Evaluate the field or potential and the derivative in one element
!    in the case of scalar elements:    the potential and the potential gradient
!    in the case of Nedelec elements:   The field vector and the curl of the field
!
!  Input:
!     elem      element number
!     lambda    vector of barycentric coordinates
!
!  Output:
!     u         vector field for Nedelec elements 
!               or u(1,:) is the potential for scalar elements 
!     curlu     curl of vector field for Nedelec elements 
!               or u(1,:) is the potential gradient for scalar elements 
!
!  local variables:
      integer (I4B) dof, i, errcode, inat
      integer (I4B), allocatable :: nff(:), vpelem(:)
      real (DP) :: vert(3,4)
      real (DP), allocatable :: xsis(:), xsi(:,:), cxsi(:,:)
!
!  vertex coordinates
      allocate(nff(nnat), vpelem(nnat))
      vert(1:3,1:4) = nod( 1:3, vn(1:4,elem) )
!
!
      u = 0._DPC
      curlu = 0._DPC

      do inat=1,nnat
        vpelem(inat) = vp(elem,inat)
        nff(inat) = size(vgdof(elem,inat)%d)

        if (eltype .eq. 'SCALAR') then
          allocate(xsis(nff(inat)), cxsi(3,nff(inat)))
          call shapefunctionsc3D(lambda, vert, 1, vpelem(inat),         &
     &                  nff(inat), .true., xsis, cxsi, errcode)
!  compute potential field for solution vector x and teh gradient
          do i=1,nff(inat)
            dof = vgdof(elem,inat)%d(i)
            if (dof .gt. 0 ) then
              u(1,inat) = u(1,inat) + xsis(i)*x(dof)
              curlu(1:3,inat) = curlu(1:3,inat) + cxsi(1:3,i)*x(dof)
            end if
          end do

          deallocate(xsis, cxsi)

        else

          allocate(xsi(3,nff(inat)), cxsi(3,nff(inat)))
!          call shapefunctionv3D(lambda, vert, 1, vpelem(inat),          &
!     &                  nff(inat), .true., xsi, cxsi, errcode)
          call shapefunctionv3D(lambda, vert, 1, vpelem(inat),          &
     &                  .true., xsi, cxsi, errcode)

     !  compute vector field for solution vector x and the curl
          do i=1,nff(inat)
            dof = vgdof(elem,inat)%d(i)
            if (dof .gt. 0 ) then
              u(1:3,inat) = u(1:3,inat) + xsi(1:3,i)*x(dof)
              curlu(1:3,inat) = curlu(1:3,inat) + cxsi(1:3,i)*x(dof)
            end if
          end do

          deallocate(xsi, cxsi)

        end if

      end do
 
      deallocate(nff, vpelem)

      return
      end subroutine field3D
