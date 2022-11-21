      subroutine sortnt(anzknt,anzzwg,xpoint,ypoint,zki,zpz,lzrb)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: anzknt, anzzwg, zki(:,:), zpz(:), lzrb(:)
      real (DP) :: xpoint(:), ypoint(:)
      intent (in) :: anzknt, anzzwg
      intent (inout) :: xpoint, ypoint, zki, zpz, lzrb
!
!    $Revision: 1.7 $
!    $Date: 2008/12/22 15:43:11 $
!    $Author: m_kasper $
!
!  sort key-points and branches from left to right and from bottom to top
!
!  In-/ Output: 
!    anzknt    ermittelte Anzahl der Knoten
!    anzzwg    ermittelte Anzahl der Zweige
!    xpoint    x-koordinaten der Knoten
!    ypoint    y-koordinaten der Knoten
!    zki       Zweig-Knoten Information
!                zki(1,i)    Anfangsknoten der Zweige
!                zki(2,i)    Endknoten der Zweige
!                zki(3,i)    Dritter Knoten (Mittelpunkt) der Zweige
!    zpz       Zahl der Knoten auf dem Zweig
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
      integer (I4B) :: i, j, k, temp, maxki, maxkj, minki, minkj
      integer (I4B) :: neuknt(anzknt), pos(anzknt), neuzwg(anzzwg)
      integer (I4B) :: iwkzw(anzzwg)
      real (DP) :: wksp(anzknt)
!
!  sorting of key-points
!     1. by increasing x-coordinate
!     2. if equal by increasing y-coordinate
!
      neuknt(1)=1
      do i=2,anzknt
        k=0
        do j=i-1,1,-1
!  sort key-points by Straight Insertion
          if (( xpoint(neuknt(j)).lt.xpoint(i) ) .or.                   &
     &      ( xpoint(neuknt(j)).eq.xpoint(i) .and.                      &
     &      ypoint(neuknt(j)).lt.ypoint(i) ) ) then
            k=j
            exit
          end if
          neuknt(j+1)=neuknt(j)
        end do
        neuknt(k+1)=i
      end do
!  copying
      do i=1,anzknt
        wksp(i)=xpoint(i)
      end do
      do i=1,anzknt
        xpoint(i)=wksp(neuknt(i))
      end do
      do i=1,anzknt
        wksp(i)=ypoint(i)
      end do
      do i=1,anzknt
        ypoint(i)=wksp(neuknt(i))
      end do
      do i=1,anzknt
        pos(neuknt(i))=i
      end do
!  correct key-point index on branches
      do i=1,anzzwg
        if (zki(3,i).eq.0) then
          zki(1,i)=pos(zki(1,i))
          zki(2,i)=pos(zki(2,i))
!  for staight lines swap direction if needed
          if (zki(2,i).lt.zki(1,i)) then
            temp=zki(1,i)
            zki(1,i)=zki(2,i)
            zki(2,i)=temp
          end if
        else
          zki(1,i)=pos(zki(1,i))
          zki(2,i)=pos(zki(2,i))
!  midpoint or third point (assign the sign)
          zki(3,i)=sign(pos(abs(zki(3,i))),zki(3,i))
        end if
      end do
!  sorting of branches
      neuzwg(1)=1
      do i=2,anzzwg
        minki=min(zki(1,i),zki(2,i))
        maxki=max(zki(1,i),zki(2,i))
        k=0
        do j=i-1,1,-1
!  sorting of branches by Straight Insertion
          minkj=min(zki(1,neuzwg(j)),zki(2,neuzwg(j)))
          maxkj=max(zki(1,neuzwg(j)),zki(2,neuzwg(j)))
          if (( minkj.lt.minki ) .or.                                   &
     &      ( minkj.eq.minki .and.                                      &
     &      maxkj.lt.maxki ) ) then
            k=j
            exit
          end if
          neuzwg(j+1)=neuzwg(j)
        end do
        neuzwg(k+1)=i
      end do
!  copying branches
      do i=1,3
        do j=1,anzzwg
          iwkzw(j)=zki(i,j)
        end do
        do j=1,anzzwg
          zki(i,j)=iwkzw(neuzwg(j))
        end do
      end do
      do i=1,anzzwg
        iwkzw(i)=zpz(i)
      end do
      do i=1,anzzwg
        zpz(i)=iwkzw(neuzwg(i))
      end do
!  boundary conditions
      do i=1,anzzwg
        iwkzw(i)=lzrb(i)
      end do
      do i=1,anzzwg
        lzrb(i)=iwkzw(neuzwg(i))
      end do
      return
      end subroutine sortnt
