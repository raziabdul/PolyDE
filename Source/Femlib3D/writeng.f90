      subroutine writeng(ok)
      use feminterface, only: getsetting
      use femtypes
      use globalvariables3D, only : sbc, dom, nnat, nod, numn, nums, numv, sn, vn
      implicit none
      logical ok
      intent (out) :: ok
!
!------------------------------------------------------------------------------
!    $Revision: 1.7 $
!    $Date: 2014/07/28 10:42:36 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  Write mesh data as neutral NETGEN mesh file. The file containing the mesh 
!  will be named basis.mesh. The information contained in this format is:
!
!  Neutral format:
!   1. Number of nodes (numn)
!   2. List of coordintes (nod)
!   3. Number of volume elements (numv)
!   4. domain (1 entry) and nodes (4 entries) for volume elements (dom, vn)
!   5. Number of surface elements (nums)
!   6. bc at surface (1 entry) and nodes (3 entries) for surface elements (bcs, sn)
!
!------------------------------------------------------------------------------
!  Input:
!    Regions
!      numdom      number of domains / regions
!      dommat      number for material in domain / region
!    Elements
!      numv        number of volume elements
!      nums        number of surface elements
!      sn          surface element -> nodes
!      vn          volume element -> nodes
!      dom         assignment of element to domain / region
!    Nodes
!      nod         world coordinates of nodes
!      numn        number of nodes
!      nsi         node - surface information
!                  NOT USED YET
!
!  Output:
!      ok          =.false. if an error occurs
!
!  local variables
      integer (I4B) :: i, unitid, ios
      real (DP) :: geoscale
      character (len=200) :: path
!
      external    ::    grglun
!
!  get project path from user's environment variables
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
!  open mesh.neutral file in project path for formatted reading (ASCII)
      open (unitid,file=path(1:len_trim(path))//'basis.mesh',            &
     &   form='formatted',position='REWIND',action='WRITE',iostat=ios)
      if (ios .ne. 0) then
        print*, '**** Error on opening file: ',path(1:len_trim(path))//'basis.mesh'
        print*, '**** IO Error No: ',ios
        ok=.false.
        return
      end if
      rewind unitid
!
!  write number of nodes
      write (unitid,*) numn
      print "(a,i8,a)","  writing ", numn," nodes"
!  write coordinates
      call getsetting('GEOMETRY_FACTOR',geoscale)
      if (geoscale .ne. 1._DP) then
        do i=1,numn
          write (unitid,*) nod(:,i)/ geoscale
        end do
      else
        do i=1,numn
          write (unitid,*) nod(:,i)
        end do
      end if

      do i=1,numn
        write (unitid,*) nod(:,i)
      end do
!
!  write number of volume elements
      write (unitid,*) numv
      print "(a,i8,a)","  writing ", numv," volume elements"
!
!  write domain and volume element - node map
      do i=1,numv
        write (unitid,*) dom(i),vn(:,i)
      end do
!
!  write number of surface elements
      write (unitid,*) nums
!
      print "(a,i8,a)","  writing ", nums," surface elements"
!
!  write domain and volume element - node map
      do i=1,nums
        write (unitid,*) sbc(i),sn(:,i)
      end do
!
      close (unitid)
      ok=.true.
      return
      end subroutine writeng
