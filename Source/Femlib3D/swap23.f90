      subroutine swap23()
      use femtypes
      use feminterface3d, only: diangles, tetmeanratio, element_draw_vtk
      use globalvariables3D,   only: vn, nod, vv, dom, numv, pi
      implicit none
!
!------------------------------------------------------------------------------
!    $Revision: 1.1 $
!    $Date: 2014/07/28 11:06:00 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  swap - performs 2-3 swap
!
!  input:
!
!  output:
!
      integer (I4B), parameter :: criterionselector=2
!  criterionselector = 1        for GEOMPACK's  eta  quality criterion
!                    = 2        for minimum dihedral angle (relative to regular tetrahedron)
      integer (I4B) :: i, j, k, nb, mi(1), newel
!  fsn(1:3,f) is the node-set defined by local face 'f'
      integer (I4B),   parameter :: fsn(3,4) = reshape((/2,3,4,1,4,3,1,2,4,1,3,2/),(/3,4/))
!  eadjn(1:3,n) is the edge-set adjacent to local node 'n'
      integer (I4B),   parameter :: eadjn(3,4) = reshape((/1,3,4,1,2,5,2,3,6,4,5,6/),(/3,4/))
!  efencen(1:3,n) is the edge-set at the fence (not adjacent) to local node 'n'
      integer (I4B),   parameter :: efencen(3,4) = reshape((/2,5,6,3,4,6,1,4,5,1,2,3/),(/3,4/))
      real (DP) :: quali_i, quali_nb, quali_1, quali_2, quali_3, nodes(3,4), di_i(6), di_nb(6)
      logical :: done(numv*1.2)
!
      done = .false.
      do i = 1, numv
        if (done(i)) cycle
        done(i) = .true.
        do j = 1, 4
          nb = vv(j,i)
          if (nb .lt. i) cycle
          if (dom(i) .ne. dom(nb)) cycle
!  find local face number in neighbour
          do k = 1, 4
            if (vv(k,nb) .eq. i) exit
          end do
!  check whether the smallest dihedral angle of element i is at the cone of element i
!  if so a 2-3 swap cannot improve dihedral angles
          di_i = diangles(nod(1:3,vn(1:4,i )))
          mi=minloc(di_i)
          if (any(eadjn(:,j) .eq. mi(1))) cycle
!  check whether the smallest dihedral angle of element nb is at the cone of element nb
!  if so a 2-3 swap cannot improve dihedral angles
          di_nb = diangles(nod(1:3,vn(1:4,nb)))
          mi=minloc(di_nb)
          if (any(eadjn(:,k) .eq. mi(1))) cycle
!  we are only safe, if dihedral angles at the fence are smaller than Pi/2
          if (any(di_i(efencen(:,j)) .ge. pi/2)) cycle
          if (any(di_nb(efencen(:,k)) .ge. pi/2)) cycle
!  calculate the quality improvement
          if (criterionselector .eq. 1) then
            quali_i  = tetmeanratio(nod(1:3,vn(1:4,i )))
            quali_nb = tetmeanratio(nod(1:3,vn(1:4,nb)))
            nodes(1:3,1) = nod(1:3,vn(j,i))
            nodes(1:3,2) = nod(1:3,vn(k,nb))
            nodes(1:3,3) = nod(1:3,vn(fsn(2,j),i))
            nodes(1:3,4) = nod(1:3,vn(fsn(1,j),i))
            quali_1  = tetmeanratio( nodes )
            nodes(1:3,3) = nod(1:3,vn(fsn(3,j),i))
            nodes(1:3,4) = nod(1:3,vn(fsn(2,j),i))
            quali_2  = tetmeanratio( nodes )
            nodes(1:3,3) = nod(1:3,vn(fsn(1,j),i))
            nodes(1:3,4) = nod(1:3,vn(fsn(3,j),i))
            quali_3  = tetmeanratio( nodes )
          else
            quali_i  = minval(di_i)/1.23095941734077468213492917825_DP
            quali_nb = minval(di_nb)/1.23095941734077468213492917825_DP
            nodes(1:3,1) = nod(1:3,vn(j,i))
            nodes(1:3,2) = nod(1:3,vn(k,nb))
            nodes(1:3,3) = nod(1:3,vn(fsn(2,j),i))
            nodes(1:3,4) = nod(1:3,vn(fsn(1,j),i))
            quali_1  = minval(diangles( nodes ))/1.23095941734077468213492917825_DP
            nodes(1:3,3) = nod(1:3,vn(fsn(3,j),i))
            nodes(1:3,4) = nod(1:3,vn(fsn(2,j),i))
            quali_2  = minval(diangles( nodes ))/1.23095941734077468213492917825_DP
            nodes(1:3,3) = nod(1:3,vn(fsn(1,j),i))
            nodes(1:3,4) = nod(1:3,vn(fsn(3,j),i))
            quali_3  = minval(diangles( nodes ))/1.23095941734077468213492917825_DP
          end if
!
          if (min(quali_i, quali_nb) .lt. min(quali_1, quali_2, quali_3) ) then
print *,' should swap',i,nb,min(quali_i, quali_nb), min(quali_1, quali_2, quali_3)
print *,' from',quali_i, quali_nb
print *,' to',quali_1, quali_2, quali_3
print *,' angles',min(quali_i, quali_nb), min(quali_1, quali_2, quali_3)
!call element_draw_vtk((/i,nb/))
!  revisit the affected elements
done(vv(1:4,i)) = .false.
done(vv(1:4,nb)) = .false.
!done(vv(1:4,newel)) = .false.
          end if

        end do

      end do
pause
      return
      end subroutine




