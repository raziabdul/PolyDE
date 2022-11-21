      subroutine setzwg(length,anzzwg,anzknt,zki,zpz,zpp)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: anzknt, anzzwg, zki(:,:), zpz(:)
      real (DP) :: length(:), zpp(:)
      intent (in) :: length, anzzwg, anzknt, zki
      intent (out) :: zpz, zpp
!
!    $Revision: 1.5 $
!    $Date: 2008/12/22 15:13:30 $
!    $Author: m_kasper $
!
!  determine the number of nodes and the distribution for all branches
!
!  Input:
!    length      Laenge der Zweige
!    anzknt      ermittelte Anzahl der Knoten
!    anzzwg      ermittelte Anzahl der Zweige
!    zki         Zweig-Knoten Information
!  Output:
!    zpz         Zahl der Punkte auf dem Zweig (incl. Anfangs- und Endpunkt)
!    zpp         Zweig-Punkt-Parameter (steuert Verteilung)
!                /  = 0  :  aequidistante Verteilung
!                -  > 0  :  Knoten am Zweiganfang verdichtet
!                \  < 0  :  Knoten am Zweigende verdichtet
!
!  local variables
      integer (I4B) :: i
      real (DP) :: lenmin, lenmax
      real (DP) :: width(anzknt)
!
      lenmin=huge(1._DP)
      lenmax=0._DP
      do i=1,anzzwg
!  kuerzeste und laengste Zweiglaenge bestimmen
        lenmin=min(length(i),lenmin)
        lenmax=max(length(i),lenmax)
      end do
      write(*,*) 'shortest branch: ',lenmax
      write(*,*) 'longest  branch: ',lenmin
!
!  Anzahl der Knoten auf den Zweigen und Knotenverschiebung bestimmen
      do i=1,anzknt
        width(i)=huge(1._DP)
      end do
!  Width gibt die gewuenschte Maschenweite in den Eingabeknoten an
      do i=1,anzzwg
!  Minimum aller angrenzenden Zweigelaengen
        width(zki(1,i))=min(width(zki(1,i)),length(i))
        width(zki(2,i))=min(width(zki(2,i)),length(i))
      end do
      do i=1,anzzwg
!     print*,'xxx',length(i),width(zki(1,i)),width(zki(2,i))
!  zpz   Anzahl der Knoten auf dem Zweig
!     print*,'#######',0.8*length(i)/(
!     *  (width(zki(1,i))+width(zki(2,i)))/2. ),
!     *   width(zki(1,i)),width(zki(2,i))
        zpz(i)=max(                                                     &
     &    int(0.8_DP*length(i)/                                         &
     &    ( (width(zki(1,i))+width(zki(2,i)))/2._DP )),                 &
     &    zpz(i), 3)
        zpz(i)=min(zpz(i),20)
!     print*,'#######',i,zpz(i),length(i)
        zpp(i)=0._DP
      end do
      return
      end subroutine setzwg
