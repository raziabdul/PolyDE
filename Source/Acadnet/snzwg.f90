      subroutine snzwg(xh,yh,zki,i,j,xp,yp,ls)
      use feminterface, only: inkrp, schnit, snkrgr, snkrkr
      use femtypes
      implicit none
      integer (I4B) :: zki(:,:), i, j
      real (DP) :: xh(:), yh(:), xp, yp
      logical :: ls
      intent (in) :: xh, yh, zki,i, j
      intent (out) :: xp, yp, ls
!
!    $Revision: 1.5 $
!    $Date: 2008/12/22 15:14:55 $
!    $Author: m_kasper $
!
!  test whether branches intersect
!
!    xh,yh   Koordinaten der Knoten
!    zki     Knotennummern der Zweige
!    i,j     Zu pruefende Zweige
!    xp,yp   Eventueller Schittpunkt
!    ls      =true , wenn es einen Schittpunkt gibt
!
!  local variables
      real (DP) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5
      real (DP) :: xa, ya, xb, yb, u, v, u1, u2, v1, v2, xp1, xp2
      real (DP) :: yp1, yp2, xpi, ypi, rausen, rho, flaec2
      integer (I4B) :: case, err
      logical :: ok
      real (DP), parameter :: eps=0.005_DP
!
      ls=.false.
      if ( (zki(3,i).eq.0) .and. (zki(3,j).eq.0) ) then
!  Fall 1 .........................   zwei Geraden
        x1=xh(zki(1,i))
        y1=yh(zki(1,i))
        x2=xh(zki(2,i))
        y2=yh(zki(2,i))
        x3=xh(zki(1,j))
        y3=yh(zki(1,j))
        x4=xh(zki(2,j))
        y4=yh(zki(2,j))
        case=1
!
      else if ( (zki(3,i).eq.0) .and. (zki(3,j).ne.0) ) then
!  Fall2..............................Gerade und Kreisbogen
        x1=xh(zki(1,i))
        y1=yh(zki(1,i))
        x2=xh(zki(2,i))
        y2=yh(zki(2,i))
        x3=xh(zki(1,j))
        y3=yh(zki(1,j))
        x4=xh(zki(2,j))
        y4=yh(zki(2,j))
        if (zki(3,j).lt.0) then
!  Kreismittelpunkt
          x0=xh(abs(zki(3,j)))
          y0=yh(abs(zki(3,j)))
        else
          x5=xh(abs(zki(3,j)))
          y5=yh(abs(zki(3,j)))
!  Kreismittelpunkt muss erst berechnet werden
          call inkrp(x3,y3,x5,y5,x4,y4,xpi,ypi,flaec2,rho,x0,y0,rausen)
        end if
        case=2
!
      else if ( (zki(3,i).ne.0) .and. (zki(3,j).eq.0) ) then
!  Fall2..............................Kreisbogen und Gerade
        x1=xh(zki(1,j))
        y1=yh(zki(1,j))
        x2=xh(zki(2,j))
        y2=yh(zki(2,j))
        x3=xh(zki(1,i))
        y3=yh(zki(1,i))
        x4=xh(zki(2,i))
        y4=yh(zki(2,i))
        if (zki(3,i).lt.0) then
!  Kreismittelpunkt
          x0=xh(abs(zki(3,i)))
          y0=yh(abs(zki(3,i)))
        else
          x5=xh(abs(zki(3,i)))
          y5=yh(abs(zki(3,i)))
!  Kreismittelpunkt muss erst berechnet werden
          call inkrp(x3,y3,x5,y5,x4,y4,xpi,ypi,flaec2,rho,x0,y0,rausen)
        end if
        case=2
!
      else
!  Fall3..............................zwei  Kreisboegen
        x1=xh(zki(1,i))
        y1=yh(zki(1,i))
        x2=xh(zki(2,i))
        y2=yh(zki(2,i))
        x3=xh(zki(1,j))
        y3=yh(zki(1,j))
        x4=xh(zki(2,j))
        y4=yh(zki(2,j))
        if (zki(3,i).lt.0) then
!  Kreismittelpunkt
          xa=xh(abs(zki(3,i)))
          ya=yh(abs(zki(3,i)))
        else
          x5=xh(abs(zki(3,i)))
          y5=yh(abs(zki(3,i)))
!  Kreismittelpunkt muss erst berechnet werden
          call inkrp(x1,y1,x5,y5,x2,y2,xpi,ypi,flaec2,rho,xa,ya,rausen)
        end if
        if (zki(3,j).lt.0) then
!  Kreismittelpunkt
          xb=xh(abs(zki(3,j)))
          yb=yh(abs(zki(3,j)))
        else
          x5=xh(abs(zki(3,j)))
          y5=yh(abs(zki(3,j)))
!  Kreismittelpunkt muss erst berechnet werden
          call inkrp(x3,y3,x5,y5,x4,y4,xpi,ypi,flaec2,rho,xb,yb,rausen)
        end if
        case=3
      end if
!
      select case (case)
      case (1)
!
!  Schnittpunkt fuer Geraden bestimmen
!
        call schnit(x1,y1,x2,y2,x3,y3,x4,y4,u,v,ok)
        if ((u.gt.eps).and.(u.lt.(1._DP-eps)).and.                      &
     &    (v.gt.eps).and.(v.lt.(1._DP-eps))) then
!  Schnittpunkt gefunden
          xp=x1+u*(x2-x1)
          yp=y1+u*(y2-y1)
          ls=.true.
        end if
!
      case(2)
!
!  Schittpunkt fuer Gerade und Kreibogen bestimmen
!
!     Print*,'Kreismittelpunkt',x0,y0
        call snkrgr(x1,y1,x2,y2,x3,y3,x4,y4,x0,y0,u1,u2,v1,v2,err)
!       print*,u1,u2,v1,v2
        if ((u1.gt.eps).and.(u1.lt.(1._DP-eps)).and.                    &
     &    (v1.gt.eps).and.(v1.lt.(1._DP-eps))) then
!  Schnittpunkt gefunden
          xp=x1+u1*(x2-x1)
          yp=y1+u1*(y2-y1)
          ls=.true.
        else if ((u2.gt.eps).and.(u2.lt.(1._DP-eps)).and.               &
     &    (v2.gt.eps).and.(v2.lt.(1._DP-eps))) then
!  Schnittpunkt gefunden
          xp=x1+u2*(x2-x1)
          yp=y1+u2*(y2-y1)
          ls=.true.
        end if
!
      case (3)
!  Schittpunkt fuer zwei Kreiboegen bestimmen
!
        call snkrkr(x1,y1,x2,y2,x3,y3,x4,y4,xa,ya,xb,yb,                &
     &    u1,u2,v1,v2,xp1,yp1,xp2,yp2,err)
        if ((u1.gt.eps).and.(u1.lt.(1._DP-eps)).and.                    &
     &    (v1.gt.eps).and.(v1.lt.(1._DP-eps))) then
!  Schnittpunkt gefunden
          xp=xp1
          yp=yp1
          ls=.true.
        else if ((u2.gt.eps).and.(u2.lt.(1._DP-eps)).and.               &
     &    (v2.gt.eps).and.(v2.lt.(1._DP-eps))) then
!  Schnittpunkt gefunden
          xp=xp2
          yp=yp2
          ls=.true.
        end if
      end select
      return
      end subroutine snzwg