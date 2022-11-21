      subroutine putnod(x,y,accur,node,xpoint,ypoint,anzknt)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: x, y, accur, xpoint(:), ypoint(:)
      integer (I4B) :: node, anzknt
      intent (in) :: x, y, accur
      intent (out) :: node
      intent (inout) :: xpoint, ypoint, anzknt
!
!    $Revision: 1.5 $
!    $Date: 2008/12/22 15:12:38 $
!    $Author: m_kasper $
!
!  add a new key-point to the key-point list if it does not already exist
!
!  Input:
!             x,y      Koordinaten den (neuen) Knotens
!             accur    Genauigkeitsschranke (raster)
!  In-/ Output: 
!             xpoint   Liste der bereits exstierenden Knoten, wird aktualisiert
!             ypoint
!             anzknt   aktualisierte Anzahl der Knoten
!  Output:
!             node     Nummer des Knotens
      real (DP) :: delx, dely
      integer (I4B) :: i
!
      do i= 1,anzknt
        delx = x-xpoint(i)
        dely = y-ypoint(i)
        if (sqrt(delx*delx + dely*dely) .lt. accur ) then
!      print*,'Knoten existiert',accur
!  Knoten existiert bereits
          node=i
          return
        end if
      end do
!  Knoten einfuegen
      anzknt=anzknt+1
      node=anzknt
      xpoint(anzknt)=x
      ypoint(anzknt)=y
      return
      end subroutine putnod
