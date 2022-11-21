      subroutine elmgeb(gbz,bzi,bzip,bzil,lzrb,layrb,matname)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: gbz, bzi(:), bzip(:), bzil(:), lzrb(:), layrb(:,:)
      character (len=*) :: matname(:)
      intent (in) :: lzrb, layrb, matname
      intent (inout) :: gbz, bzi, bzip, bzil
!
!    $Revision: 1.6 $
!    $Date: 2008/12/22 15:51:46 $
!    $Author: m_kasper $
!
!   eliminate regions which only have boundary branches ("holes")
!
!  In-/ Output:
!    gbz       number of regions
!    bzi       is a continous list of branches, including the inner branches. 
!              The sequence is the same as in the input file (netin.acd)
!    bzip      is a pointer (index) to the first branch in the list bzi of a 
!              region. For multiply connected regions inner boundary branches 
!              precede the outer boundary.  
!    bzil      for each region is a pointer (index) to the first 
!              (possibly inner) region in the in the array bzip
!  Input:
!    lzrb      for each branch the layer
!    layrb     type of boundary condition for branches on that layer (and for each nature)
!              0    = Dirichlet
!              1-99 = User Dirichlet
!              200  = general Neumann
!              300  = innner contour
!              400  = non-visble contour
!              1000 = comment layer
!              1001 = not recognized layer, treated as a comment layer
!
!  local variables
      integer (I4B) :: ggeb, i, j, k, kbndry, kinner, zanz, del, reg
      logical :: elim(gbz)
!
!  all regions which don't have a material beeing assigned yet are checked
!  whether they are surrounded only by boundary branches (Dirichlet Neumann or
!  general BC) if so these regions are considered to be holes and are eliminated
!
!  only one region so it cannot be a hole
      if (gbz.eq.1) return
!
!  loop over all regions and determine whether they are holes
      elim=.false.
      do i=1,gbz
!  only considere regions which do not have a material name
        if (matname(i) .ne. 'nothing') cycle
!  only considere regions which don't enclose further regions
        if (bzil(i+1)-bzil(i) .gt. 1) cycle
        kbndry=0
        kinner=0
        j=bzil(i)
        do k=bzip(j),bzip(j+1)-1
!  check for all branches of the region whether they have a boundary condition
          if (any( layrb(lzrb(abs(bzi(k))),:) .ge. 300)) then
            kinner=kinner+1
          else
            kbndry=kbndry+1
          end if
        end do
!
        if (kinner.ne.0) then
!  this region is not a hole, however it does not habe a materialname
          write(*,100) ' *WARNING: no material specified for region',i
100       format(a,i4)
          print*,'          the materialname "nothing" will be used'
          write(*,101) '           material has to be specified if',    &
     &                 ' regions are sperated by a branch'
101       format(a,a)
        else
!  this region is considered a hole, it is to be eliminated
          elim(i)=.true.
        end if
      end do
!
!  now eliminate regions
      ggeb=gbz
      del=0
      do i=1,gbz
        if (elim(i)) then
!  correct arrys bzi, bzip and bzil
!  copy the array eliminating the actual region
          reg=i-del
          j=bzil(reg)
          zanz=bzip(j+1)-bzip(j)
          do k=bzip(j+1),bzip(bzil(ggeb+1))-1
            bzi(k-zanz)=bzi(k)
          end do
          do k=bzil(reg),bzil(ggeb+1)
            bzip(k)= bzip(k+1)-zanz
          end do
          do k=reg,ggeb
            bzil(k)=bzil(k+1)-1
          end do
          ggeb=ggeb-1
          del=del+1
        end if
      end do
      gbz=ggeb
      return
      end subroutine elmgeb
