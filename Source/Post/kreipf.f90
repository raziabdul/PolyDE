      subroutine kreipf(px,py,xp,yp,maxlen,vx,vy,acolor)
      use feminterface, only: pfeil2, pfeil
      use femtypes
      implicit none
      integer (I4B) :: px, py, acolor
      real (DP) :: maxlen
      real (DP) :: xp(:,:), yp(:,:), vx(:,:), vy(:,:)
      intent (in) :: px, py, xp, yp, maxlen, vx, vy, acolor

!    $Revision: 1.9 $
!    $Date: 2010/09/06 12:34:21 $
!    $Author: m_kasper $
!
!----------------------------------------------------------------------
!  Plots the arrow at the specified position
!----------------------------------------------------------------------
!
!  Input:
!     px, py    number of intended points in x and y direction
!     xp, yp    position of arrow
!     maxlen    maximum length of arrow
!     vx, vy    vector component in x and y direction
!     acolor    color of the arrow (not yet used)
!
!  local variables
      integer (I4B) :: i, j
      real (DP) :: lenx, leny, pfwink, pfkopf, pfbr, maxbr
      real (DP) :: vmax,  vlen(px,py)
! pi/180
      real (DP), parameter :: pid180=0.01745329251994329576923691_DP
!
      external :: pgdraw, pgmove, pgqci, pgsci, pgscr


      do i=1,px
        do j=1,py
          vlen(i,j) = sqrt(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
        end do
      end do
!
      vmax = maxval(vlen)
!

!  opening angle of the arrow head in degrees (between center-line and slope)
      pfwink = 20._DP
!  maximum lngth of the arrow head = 50% of the maximum arrow length
      pfkopf = 0.5_DP * maxlen
!  maximum arrow (body) width = 95%  of width of the arrow head
      maxbr = 0.95_DP * pfkopf*tan(pfwink*pid180)
!
      do  i=1,px
        do j=1,py
          pfbr = vlen(i,j)/vmax*maxbr
          if (pfbr .gt. maxbr/10._DP) then
!  Thicker arrow
            lenx = maxlen*vx(i,j)/vlen(i,j)
            leny = maxlen*vy(i,j)/vlen(i,j)
            call pfeil2(xp(i,j),yp(i,j),lenx,leny,pfwink,pfkopf,pfbr,.true.)
          else
            if (vlen(i,j)*maxlen/vmax*10._DP .gt. 1.2*pfkopf) then
!  Arrow with zero thickness
              lenx = maxlen*vx(i,j)/vmax*10._DP
              leny = maxlen*vy(i,j)/vmax*10._DP
              call pfeil(xp(i,j),yp(i,j),lenx,leny,pfwink,pfkopf,.true.,.true.)
            else if (vlen(i,j)*maxlen/vmax*100._DP .gt. 1.2*pfkopf) then
!  Only one line
              lenx = maxlen*vx(i,j)/vmax*10._DP
              leny = maxlen*vy(i,j)/vmax*10._DP
              call pfeil(xp(i,j),yp(i,j),lenx,leny,pfwink,pfkopf,.false.,.true.)
            else
!  Nothing
            end if
          end if
        end do
      end do
!
      return
      end subroutine kreipf
!
!
!
      subroutine pfeil2(x0,y0,lenx,leny,pfwink,pfkopf,pfbrt,mitte)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) :: x0, y0, lenx, leny, pfwink, pfkopf, pfbrt
      logical :: mitte
      intent (in) :: x0, y0, lenx, leny, pfwink, pfkopf, pfbrt, mitte
!  Draw an arrow with thick body
!                *
!         ********  *
!   foot  * Body      *   Head
!         ********  *
!                *
!
!  x0,y0        locus to draw the arrow
!  lenx,leny    x- and y-length of arrow (including head)
!  pfwink       opening angle (slope) of the arrow head
!  pfkopf       length of the head
!  pfbrt        half width of the arrow body
!  mitte        = .true   the arrow is centered at locus (x0,y0)
!               = .false. (x0,y0) is the foot of the arrow
!
!  local variables
      real (DP) :: gamma, xa, ya, xe, ye, wink, cospf, pfklen
! pi/180
      real (DP), parameter :: pid180=0.01745329251994329576923691_DP
!
      external :: pgdraw, pgmove, pgqci, pgsci, pgscr

      if (mitte) then
!  Starting point
        xa = x0-lenx/2._DP
        ya = y0-leny/2._DP
!  End point
        xe = x0+lenx/2._DP
        ye = y0+leny/2._DP
      else
!  Starting point
        xa = x0
        ya = y0
!  End point
        xe = x0+lenx
        ye = y0+leny
      end if
!
!  draw head
      wink = pid180*pfwink
      cospf = cos(wink)
!  length of limbs
      pfklen = pfkopf/cospf
      if (pfklen*sin(wink) .lt. pfbrt) then
        print*,'** arrow width is smaller than head width'
      end if
      gamma=atan2(leny,lenx)
!                3
!         1******2  *
!         *           4
!         0******6  *
!                5
      call pgmove(real(xa+pfbrt*sin(gamma)),                            &
     &            real(ya-pfbrt*cos(gamma)))
      call pgdraw(real(xa-pfbrt*sin(gamma)),                            &
     &            real(ya+pfbrt*cos(gamma)))
      call pgdraw(real(xe-pfkopf*cos(gamma)-pfbrt*sin(gamma)),          &
     &            real(ye-pfkopf*sin(gamma)+pfbrt*cos(gamma)))
      call pgdraw(real(xe-pfklen*cos(gamma-wink)),                      &
     &            real(ye-pfklen*sin(gamma-wink)))
      call pgdraw(real(xe),real(ye))
      call pgdraw(real(xe-pfklen*cos(gamma+wink)),                      &
     &            real(ye-pfklen*sin(gamma+wink)))
      call pgdraw(real(xe-pfkopf*cos(gamma)+pfbrt*sin(gamma)),          &
     &            real(ye-pfkopf*sin(gamma)-pfbrt*cos(gamma)))
      call pgdraw(real(xa+pfbrt*sin(gamma)),                            &
     &            real(ya-pfbrt*cos(gamma)))
      return
      end subroutine pfeil2
!
!
!
      subroutine pfeil(x0,y0,lenx,leny,pfwink,pfkopf,head,mitte)
      use feminterface, only :
      use femtypes
      implicit none
      real (DP) :: x0, y0, lenx, leny, pfwink, pfkopf
      logical :: head, mitte
      intent (in) :: x0, y0, lenx, leny, pfwink, pfkopf, head, mitte
!  Draw an arrow - the body is a line
!                *
!                *  *
!   foot  ********     *      Head
!                *  *
!                *
!
!  x0,y0        locus to draw the arrow
!  lenx,leny    x- and y-length of arrow (including head)
!  pfwink       opening angle (slope) of the arrow head
!  pfkopf       length of the head
!  head         = .false. if the arrow is drawn without head (just a line)
!  mitte        = .true   the arrow is centered at locus (x0,y0)
!               = .false. (x0,y0) is the foot of the arrow
!
!  local variables
      real (DP) :: gamma, xa, ya, xe, ye, wink, cospf, pfklen
! pi/180
      real (DP), parameter :: pid180=0.01745329251994329576923691_DP
!
      external :: pgdraw, pgmove, pgqci, pgsci, pgscr

       if (mitte) then
!  Starting point
        xa = x0-lenx/2._DP
        ya = y0-leny/2._DP
!  End point
        xe = x0+lenx/2._DP
        ye = y0+leny/2._DP
      else
!  Starting point
        xa = x0
        ya = y0
!  End point
        xe = x0+lenx
        ye = y0+leny
      end if
!
      if (head) then
!  draw head
        wink = pid180*pfwink
        cospf = cos(wink)
!  length of limbs
        pfklen = pfkopf/cospf
        gamma = atan2(leny,lenx)
!                4
!                *  *
!        0******1,5     3      Head
!                *  *
!                2
!  arrow
        call pgmove(real(xa),real(ya))
!  head
        call pgdraw(real(xe-pfkopf*cos(gamma)),                         &
     &              real(ye-pfkopf*sin(gamma)))
        call pgmove(real(xe-pfklen*cos(gamma-wink)),                    &
     &              real(ye-pfklen*sin(gamma-wink)))
        call pgdraw(real(xe),real(ye))
        call pgdraw(real(xe-pfklen*cos(gamma+wink)),                    &
     &              real(ye-pfklen*sin(gamma+wink)))
        call pgdraw(real(xe-pfklen*cos(gamma-wink)),                    &
     &              real(ye-pfklen*sin(gamma-wink)))
      else
!  no head
        call pgmove(real(xa),real(ya))
        call pgdraw(real(xe),real(ye))
      end if
      return
      end subroutine pfeil
!
!
!
      subroutine tensorgrid(px,py,xp,yp,maxlen,vx,vy,vxy,acolor)
      use feminterface, only: principlevalues
      use femtypes
      implicit none
      integer (I4B) :: px, py, acolor
      real (DP) :: maxlen
      real (DP) :: xp(:,:), yp(:,:), vx(:,:), vy(:,:), vxy(:,:)
      intent (in) :: px, py, xp, yp, maxlen, vx, vy, acolor
!----------------------------------------------------------------------
!  Plots a grid of tensors at the specified position
!----------------------------------------------------------------------
!
!  Input:
!     px, py         number of intended points in x and y direction
!     xp, yp         position of tensor
!     maxlen         maximum length of tensors
!     vx, vy, vxy    tensor components in x and y direction and shear component vxy
!     acolor         color of the tensor (not yet used)
!
!  local variables
      integer (I4B) :: i, j
      real (DP) :: xa, xe, xf, ya, ye, yf
      real (DP) :: vmax, lambda1(px,px), lambda2(px,px), angle(px,py)
! pi/180
      real (DP), parameter :: pid180=0.01745329251994329576923691_DP
!
      external :: pgdraw, pgmove, pgqci, pgsci, pgscr

!  get principle values and orientationn (angle) of first principle axis
      do i=1,px
        do j=1,py
          call principlevalues(vx(i,j),vy(i,j),vxy(i,j),lambda1(i,j),lambda2(i,j),angle(i,j))
        end do
      end do
!
      vmax = maxval(abs(lambda1))
!
      do  i=1,px
        do j=1,py
          if (abs(lambda1(i,j))/vmax .gt. 0.02) then
!  the tensor is representes as a rhombus (diamond) length and orientation of the
!  diagonals correspond to the orientation and magnitude in the principle axes of the tensor
!  Starting point
            xa = xp(i,j) - maxlen/vmax * lambda1(i,j) * cos(angle(i,j))/2._DP
            ya = yp(i,j) - maxlen/vmax * lambda1(i,j) * sin(angle(i,j))/2._DP
            call pgmove(real(xa),real(ya))
            xe = xp(i,j) + maxlen/vmax * lambda1(i,j) * cos(angle(i,j))/2._DP
            ye = yp(i,j) + maxlen/vmax * lambda1(i,j) * sin(angle(i,j))/2._DP
            call pgdraw(real(xe),real(ye))
            xf = xp(i,j) - maxlen/vmax * lambda2(i,j) * sin(angle(i,j))/2._DP
            yf = yp(i,j) + maxlen/vmax * lambda2(i,j) * cos(angle(i,j))/2._DP
            call pgdraw(real(xf),real(yf))
            call pgdraw(real(xa),real(ya))
            xf = xp(i,j) + maxlen/vmax * lambda2(i,j) * sin(angle(i,j))/2._DP
            yf = yp(i,j) - maxlen/vmax * lambda2(i,j) * cos(angle(i,j))/2._DP
            call pgdraw(real(xf),real(yf))
            call pgdraw(real(xe),real(ye))
          end if
        end do
      end do
!
      return
      end subroutine tensorgrid
!
!
!
      subroutine plottensor(px,py,xa,xe,ya,ye)
      use feminterface, only: elemnt, getpostsetting, fieldquantity, xy2lam
      use feminterface, only: tensorgrid, convertcolorstring
      use femtypes
      use globalvariables
      implicit none
      integer (I4B) px, py
      real (DP) xa, xe, ya, ye
      intent (in) :: px, py, xa, xe, ya, ye
!
!----------------------------------------------------------------------------
!  Fetches the position of the tensor and plots the tensor by calling tensorgrid
!----------------------------------------------------------------------------
!
!  Input:
!     px, py    number of intended points in x and y direction
!     xa, xe    plotting range in x direction (start and end point)
!     ya, ye    plotting range in y direction (start and end point)
!
!  local variables:
      integer (I4B) :: i, j, ielem, acolor, oldindex
      real (SP) :: r = 0.0, g = 0.0, b = 0.0
      real (SP) :: ocr, ocg, ocb
      real (DP) :: xp(px,py), yp(px,py), maxlen, phi, lambda(3)
      real (DP) :: vx(px,py), vy(px,py), vxy(px,py)
      complex (DPC) :: zten(2,2)
      logical :: ok, ok1
      character (len=100) :: fieldtype
      character (len=10) :: unit
      character (len=50) :: descriptor
!
      external :: pgdraw, pgmove, pgqci, pgqcr, pgsci, pgscr

!  get field quantity and phase from POSTsettings.txt
      call getpostsetting('FIELDTYPE',fieldtype)
      call getpostsetting('PHI',phi)
!
      do i=1,px
        do j=1,py
!  arrange points in x/y grid:   | - * - - * - - * - |
          xp(i,j) = xa + (real(i,DP)-0.5_DP)*(xe-xa)/real(px,DP)
          yp(i,j) = ya + (real(j,DP)-0.5_DP)*(ye-ya)/real(py,DP)
          call elemnt(xn,yn,e,n,xp(i,j),yp(i,j),ielem,en,ok)
          if (.not. ok) then
!  point not in an element. set values to 0
            vx(i,j) = 0._DP
            vy(i,j) = 0._DP
            vxy(i,j) = 0._DP
          else
!  fetch the components of the vector
            call xy2lam(xp(i,j),yp(i,j),ielem,lambda,xn,yn,e)
            call fieldquantity(ielem,fieldtype,lambda,phi,zten,descriptor,unit,ok1)
            if (.not.ok1) then
              print *,"***** fieldtype(tensor): ",trim(fieldtype)," unknown for ",trim(physics)
              stop
            end if
            vx(i,j) = real(zten(1,1),DP)
            vy(i,j) = real(zten(2,2),DP)
            vxy(i,j) = real(zten(1,2),DP)
          end if
        end do
      end do
!
!  determin maximum legth e.g. 120% of spacing
      maxlen = min((xe-xa)/px,(ye-ya)/py)*1.2_DP
      call convertcolorstring(r, g, b)
      call pgqci(oldindex)
      call pgqcr(21, ocr, ocg, ocb)
      call pgscr(21, r, g, b)
      call pgsci(21)
      acolor = 0
!
!  now draw the tensor plot
      call tensorgrid(px,py,xp,yp,maxlen,vx,vy,vxy,acolor)
!  reset color index
      call pgscr(21, ocr, ocg, ocb)
      call pgsci(oldindex)
      return
      end subroutine plottensor
