      subroutine getbcval( face, nature, xyzs, pval, qval )
      use feminterface3D, only: userbc3D
      use femtypes
      use globalvariables3d, only: bctype, pvalue, qvalue
      implicit none
      integer (I4B) :: face, nature
      real (DP) :: xyzs(3)
      complex(DPC) :: pval
      complex(DPC), optional :: qval(:)
      intent (in)   :: face, nature, xyzs
      intent (out)  :: pval, qval
!
!    $Revision: 1.5 $
!    $Date: 2014/08/22 11:15:05 $
!    $Author: m_kasper $
!
!  calculate values for boundary condition for given the BC number if 
!  edge or face and coordinates
!
!  Input:
!        face          face of this bc
!        nature        nature for which the bc is evaluated
!        xyzs          coordinates for determining dirichlet
!
!  Output:
!        pval           value for either the Dirichlet bc (p1) or General Neumann
!                       
!        qval(:)        values for the General Neumann bc.
!

!  If there is a homogenous Dirichlet bc or Neuman then assign values directly. 
!  Call userbc3d if bcnum is 1 to 99, 201-299. Print error in other cases.

!  zero bcnum refers to homogenous Dirichlet BC, 200 refers to homogenous Neumann
!  local variables
      integer (I4B) bcnum


      bcnum = bctype(face,nature)
      if ( bcnum .eq. 0) then
        pval = pvalue(face,nature)
      else if ( bcnum .eq. 200) then
        pval = pvalue(face,nature)
        qval = 0._DPC
!  for now, only use the own nature, the other are zero
        qval(nature) = qvalue(face,nature)
      else if (bcnum .ge.   1 .and. bcnum .lt. 100) then
        call userbc3D(face, nature, bcnum, xyzs, pval )
      else if (bcnum .ge. 201 .and. bcnum .lt. 300) then
        call userbc3D( face, nature,bcnum, xyzs, pval, qval )
      else
        print *, "No dirichlet boundary condition defined on this face."
      end if
!
      return
      end subroutine getbcval



      subroutine getbcval_vec( face, nature, xyzs, pval, qval )
      use feminterface3D, only: userbc3D_vec
      use femtypes
      use globalvariables3d, only: bctype, pvalue_vec, qvalue_vec
      implicit none
      integer (I4B) :: face, nature
      real (DP) :: xyzs(3)
      complex(DPC) :: pval(3)
      complex(DPC), optional :: qval(3,3)
      intent (in)   :: face, nature, xyzs
      intent (out)  :: pval, qval
!
!  calculate values for boundary condition for given the BC number if 
!  edge or face and coordinates
!
!  Input:
!        face          face of this bc
!        nature        nature for which the bc is evaluated
!        xyzs          coordinates for determining dirichlet
!
!  Output:
!        pval           value for either the Dirichlet bc (p1) or General Neumann
!                       
!        qval(:)        values for the General Neumann bc.
!

!  If there is a homogenous Dirichlet bc or Neuman then assign values directly. 
!  Call userbc3d if bcnum is 1 to 99, 201-299. Print error in other cases.

!  zero bcnum refers to homogenous Dirichlet BC, 200 refers to homogenous Neumann
!
!  local variables
      integer (I4B) bcnum, i

      bcnum = bctype(face,nature)
      if ( bcnum .eq. 0) then
        pval = pvalue_vec(:,face,nature)
      else if ( bcnum .eq. 200) then
        pval = pvalue_vec(:,face,nature)
        qval = 0._DPC
!  for now, only use the own nature, the other are zero
!  and we assume the tensor q only to be occupied in the diagonal
        do i=1,3
          qval(i,i) = qvalue_vec(i,face,nature)
        end do
      else if (bcnum .ge.   1 .and. bcnum .lt. 100) then
        call userbc3D_vec( bcnum, xyzs, pval )
      else if (bcnum .ge. 201 .and. bcnum .lt. 300) then
        call userbc3D_vec( bcnum, xyzs, pval, qval )
      else
        print *, "No dirichlet boundary condition defined on this face."
      end if
!
      return
      end subroutine getbcval_vec
