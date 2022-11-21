      subroutine testac(anzzwg,anzknt,xpoint,ypoint,zki,lzrb,layrb)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: anzzwg, anzknt, zki(:,:), layrb(:,:), lzrb(:)
      real (DP) :: xpoint(:), ypoint(:)
      intent (in) :: anzzwg, anzknt, xpoint, ypoint, zki, lzrb, layrb
!
!    $Revision: 1.8 $
!    $Date: 2014/07/01 14:39:43 $
!    $Author: m_kasper $
!
!  some topological test of mesh input
!
!  Input:
!    anzzwg    number of branches
!    anzknt    number of key-points
!    xpoint    x-coordinates of key-points
!    ypoint    x-coordinates of key-points
!    zki       key-points (geometry nodes) which determine the branches
!               zki(1,i)  start key-point of the branch
!               zki(2,i)  end key-point of the branch
!               zki(3,i)  third key-point of the branch, with
!                         /  = 0  :  straight line 
!               zki(3,i)  -  > 0  :  node located on an arc              
!                         \  < 0  :  midpoint of an arc
!    lzrb      for each branch the layer
!    layrb     type of boundary condition for branches on that layer (and for each nature)
!              0    = Dirichlet
!              1-99 = User Dirichlet
!              200  = general Neumann
!              300  = innner contour
!              400  = non-visble contour
!              1000 = comment layer
!              1001 = not recognized layer, treated as a comment layer
!  Output:
!                    only as text
!
!  local variables
      integer (I4B) :: i, nodes, ungrad, rndknt
      integer (I4B) :: pruef1(anzknt), pruef2(anzknt), pruef3(anzknt)
!
!  Test1:
!      simple test using the Euler characteristic
!      - number of key-point <= number of braches
!      - number of key-points with odd degree must be even
!
!  Test2:
!      test whether each key-point, which ist not a 'third point',
!      is linked to at lest two branches
!
!  Test3:
!      - at each key-point the number of adjacent boundary branches
!        (Dirichlet, Neumann, ..) has to be zero or two
!      - the number of boundary nodes must be larger than one
!
      if (anzknt.lt.2 .or. anzzwg.lt.2) then
        print*,'*ERROR: only ',anzknt,' key-point(s) and ',anzzwg,' branch(es)'
        return
      end if
!
      do i=1,anzknt
        pruef1(i)=0
        pruef2(i)=0
        pruef3(i)=0
      end do
      do i=1,anzzwg
        pruef1(zki(1,i))=pruef1(zki(1,i))+1
        pruef1(zki(2,i))=pruef1(zki(2,i))+1
!
        pruef2(zki(1,i))=pruef2(zki(1,i))+1
        pruef2(zki(2,i))=pruef2(zki(2,i))+1
        if (zki(3,i) .ne. 0) then
          pruef2(abs(zki(3,i)))=100
        end if
!
        if ( any( layrb(lzrb(i),:) .le. 299) .and.                      &
     &       all( layrb(lzrb(i),:) .ge. 0) ) then
!  this is a boundary branch
          pruef3(zki(1,i))=pruef3(zki(1,i))+1
          pruef3(zki(2,i))=pruef3(zki(2,i))+1
        end if
      end do
!
      nodes=0
      ungrad=0
      rndknt=0
      do i=1,anzknt
        if (pruef1(i) .eq. 0) then
          nodes=nodes+1
        else if (mod(pruef1(i),2).eq.1) then
          ungrad=ungrad+1
        end if
!
        if (pruef2(i) .le. 1) then
          write(*,1232) pruef2(i)
1232      format(1x,' *WARNING: key-point only has ',i2,' adjacent branches')
          print*,' coordinates:',xpoint(i),ypoint(i)
        end if
!
        if (pruef3(i) .gt. 0) then
          rndknt=rndknt+1
        end if
!
        if ( (pruef3(i) .eq. 1) .and. (pruef2(i) .gt. 1) ) then
          print*,'*WARNING: incorrect boundary condition?'
          write(*,1233) i
1233      format(1x,'key-point',i6,' only has one adjacent boundary branch')
          write (*,1235) '  coordinates:',xpoint(i),ypoint(i)
1235      format(a,2g23.13)
        end if
!
        if (pruef3(i) .gt. 2) then
          print*,'*WARNING: incorrect boundary condition?'
          write(*,1234) i,pruef3(i)
1234      format(1x,'key-point',i6,' has ',i2,' adjacent boundary branches')
          write (*,1235) '  coordinates:',xpoint(i),ypoint(i)
        end if
!
      end do
      if (rndknt .lt. 2) then
        print*,'*WARNING: incorrect boundary condition'
      end if
      if ( (nodes .gt. anzzwg) .or. (mod(ungrad,2) .eq. 1) ) then
        print*,'*ERROR: incorrect structure'
      end if
      print*
!
      return
      end subroutine testac
