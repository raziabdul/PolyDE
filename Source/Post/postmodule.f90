      module postmodule
      implicit none
      public
      contains
!
!
!
      subroutine scale(palette,ncolor,fmin,fmax)
      use femtypes
      use feminterface, only:
      use globalvariables
      implicit none
      integer (I4B) :: ncolor
      integer (I4B) :: palette(:)
      real (DP) :: fmin,fmax
      intent (IN) :: fmin, fmax, palette, ncolor
!
!------------------------------------------------------------------------------
!
!    PolyDE -- a finite element simulation tool
!    Copyright (C) 2006 Institute for Micro Systems Technology,
!                       Hamburg University of Technology.
!
!    This file is part of PolyDE.
!
!    PolyDE is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by the Free Software Foundation; either version 2, or (at your
!    option) any later version.
!
!    PolyDE is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!    MA 02110-1301 USA.
!
!    $Revision: 1.18 $
!    $Date: 2014/07/01 14:44:45 $
!    $Author: m_kasper $
!
!------------------------------------------------------------------------------
!
!  Input:
!     fmin      minimum value for scale
!     fmax      maximum value for scale
!     ncolor    number of colors to be used
!     palette   color array in increasing order
!
!  local variables:
      integer (I4B) :: i, oldindex, rgb, zff
      real (SP) :: wxa, wxb, wya, wyb, xa,xb,ya,yb
      real (SP) :: vxa,vxb,vya,vyb
      real (SP) :: oldch, xch, ych
      real (SP) :: cr, cg, cb, ocr, ocg, ocb
      real (SP) :: xsmin, xsmax, ysmin, ysmax
      parameter (zff=z'ff')
!
      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgsfs, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgqcs, pgqwin, pgqvp, pgsvp, pgswin, pgwnad
      external :: pgrect, pgptxt 

! store the current world and viewport coords and the character height
      call pgqwin(wxa, wxb, wya, wyb)
      call pgqvp(0, xa, xb, ya, yb)
      call pgqch(oldch)
! determine the unit character height in ndc (normalized device coordinates) coords
      call pgsch(0.7)
      call pgqcs(0, xch, ych)
! use these to determine viewport coordinates for the wedge + annotation
      vxa = xa
      vxb = xb
      vyb = ya -(2.0*ych*oldch)
      vya = vyb -(0.4*3.0*ych*oldch)
! set the viewport for the color scale
      call pgsvp(vxa, vxb, vya, vyb)
! convert xmin,xmax,ymin,and ymax to real               
      xsmin = real(xmin)
      xsmax = real(xmax)
      ysmin = real(ymin)
      ysmax = real(ymax)
! set the window to be map to the viewport
      call pgswin(xsmin,xsmax,ysmin,ysmax)
! store the current color index
      call pgqci(oldindex)
! change the color index to 20, and store its color representation
      call pgsci(20)
      call pgqcr(20, ocr, ocg, ocb)
! set the fill type to solid
      call pgsfs(1)
! draw the rectangles and fill with colors
      do i = 1,abs(ncolor)
        rgb = palette(i)
        cr = real( iand(rgb           , zff), 8) /255.
        cg = real( iand(ishft(rgb, -8), zff), 8) /255.
        cb = real( iand(ishft(rgb,-16), zff), 8) /255.
        call pgscr(20, cr, cg, cb)
        call pgrect(xsmin+(i-1)*(xsmax-xsmin)/abs(ncolor),xsmin+i*(xsmax-xsmin)&
     &              /abs(ncolor), ysmin, ysmax)
      end do
! reset the color index
      call pgsci(oldindex)
! reset the color representation of color index 20
      call pgscr(20, ocr, ocg, ocb)
! set window
      call pgswin(real(fmin),real(fmax),0.0,1.0)
! draw labeled frame around viewport
      call pgqci(oldindex)
      call pgqcr(21,ocr,ocg,ocb)
      call pgscr(21,0.0,0.0,0.0)
      call pgsci(21)
      call pgbox('BCNST',0.0,0,'BC',0.0,0)
      call pgscr(21,ocr,ocg,ocb)
      call pgsci(oldindex)
! reset the viewport, window and character height
      call pgsvp(xa,xb,ya,yb)
      call pgswin(wxa,wxb,wya,wyb)
      call pgsch(oldch)
!
      end subroutine scale
!
!
!
      subroutine info
      use feminterface, only: getpostsetting, getsetting
      use femtypes
      use globalvariables
      implicit none
!
!  Routine that prints information about the plot on the right side of the graphics
!  window. Displayed information contains date and time, path, physics mode, 
!  frequency and error of the plot
!
!  Input:
!     res       squared residual for each element
!     sumref    value needed for computation of relative estimated error
!
!  local variables
      integer (I4B) :: a, b
      real (SP) :: vx1, vx2, vy1, vy2, wx1, wx2, wy1, wy2
      character (len=4)   :: information
      character (len=8)   :: dat
      character (len=10)  :: tim
      character (len=16)  :: displaydevice
      character (len=20)  :: tmp
      character (len=30)  :: short_path
      character (len=200) :: path
!
      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgptxt, pgqvp, pgqwin, pgswin, pgsvp, pgwnad
      external :: pglab, pgenv
      
!      
      call getpostsetting('INFO',information)
      call getpostsetting('DISPLAYDEVICE',displaydevice)
      displaydevice = displaydevice(1:len_trim(displaydevice))
      if (information .eq. 'NO') return
!  query viewport and window
      call pgqvp(0,vx1,vx2,vy1,vy2)
      call pgqwin(wx1,wx2,wy1,wy2)
      if ((displaydevice .eq. '/VCPS') .or. (displaydevice .eq. '/VPS')) then
        call pgsvp(0.0, 1.0, 0.0, 0.18)
      else
        call pgsvp(0.8, 1.0, 0.0,1.0)
      endif
      call pgwnad (0.0,1.0,0.0,1.0)
!  get time and date for information plot
      call date_and_time(dat,tim)
      call pgscf(1)
      call pgsch(0.7)
      call getsetting('PROJECTPATH',path)
!  now extracting last 2 folders from the path to short_path
      a = len(path)
      b = 0
      do
        if (b .eq. 3) exit
        if (whatsystem .eq. 'WINDOWS') then
          if((path(a:a) .ne. '\') .and. (b .lt. 3)) then
            a = a - 1
          else if((path(a:a) .eq. '\') .and. (b .lt. 3)) then
            a = a - 1
            b = b + 1
          end if
        else if (whatsystem .eq. 'LINUX') then
          if((path(a:a) .ne. '/') .and. (b .lt. 3)) then
            a = a - 1
          else if((path(a:a) .eq. '/') .and. (b .lt. 3)) then
            a = a - 1
            b = b + 1
          end if
        end if
      end do
      short_path = path((a+1):len(path))
      call pgwnad (0.0,1.0,0.0,1.0)
!  plot date and time
      call pgptxt(-1.0, 0.90, 0.0, 0.0, 'Date = '//dat(7:8)//'.'//dat(5:6)//'.'//dat(1:4))
      call pgptxt(-1.0, 0.80, 0.0, 0.0, 'Time = '//tim(1:2)//':'//tim(3:4))
!  plot project path and physics mode
      call pgptxt(-1.0, 0.70, 0.0, 0.0, 'Project Path = '//short_path)
      call pgptxt(-1.0, 0.60, 0.0, 0.0, 'Physics Mode ='//physics)
!  plot frequency in Hz
      write (tmp, '(a,g10.3)') '',omega/(2._DP*pi)
      call pgptxt(-1.0, 0.50, 0.0, 0.0, 'Frequency in Hz = '//tmp)
!  plot relative estimated error of the solution as
      write (tmp,'(a,g10.3)') '',100._DP*fem_accuracy
      call pgptxt(-1.0, 0.30, 0.0, 0.0, 'Relative Error in % ='//tmp)
!  reset viewport and window
      call pgsvp(vx1,vx2,vy1,vy2)
      call pgswin(wx1,wx2,wy1,wy2)
!
      end subroutine info
!
!
!
      subroutine label(fieldtype)
      use feminterface, only: fieldquantity, low2hi, convertcolorstring
      use femtypes
      use globalvariables
      implicit none
      character (len=*) :: fieldtype
      intent(in)  :: fieldtype
!
!  Subroutine to deal with title for plots. What is written in the label depends
!  on physics and selected quantity
!
! Input:
!     fieldtype      read from FIELDTYPE in POSTsettings.txt
!
! local variables:
      integer(I4B) :: oldindex
      real (SP) :: r = 0.0, g = 0.0, b = 0.0
      complex (DPC) :: zs
      logical :: ok
      character (len=10) :: unit
      character (len=50) :: descriptor

      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgqwin, pgsvp, pgwnad, pglab, pgenv


!------------------------------------------------------------------------------
!  retrieve descriptor and unit of fieldquantity
      call fieldquantity(1,fieldtype,(/1._DP,0._DP,0._DP/),0._DP,zs,descriptor,unit,ok)
      if (.not.ok) then
         print *,"***** fieldtype: ",trim(fieldtype)," unknown for ",trim(physics)
         stop
      end if
      call low2hi(descriptor,len(descriptor))
      call pgsch(1.0)
      call pgsci(oldindex)
      call convertcolorstring(r, g, b)
      call pgscr(21, r, g, b)
      call pgsci(21)
      if (unit .ne. '') then
        call pglab('', '', descriptor(1:len_trim(descriptor))//'; '//  &
     &             fieldtype(1:len_trim(fieldtype))//' in ['//unit(1:len_trim(unit))//']')
      else
        call pglab('', '', descriptor(1:len_trim(descriptor))//' '//  &
     &             fieldtype(1:len_trim(fieldtype)))
      end if
      call pgsci(oldindex)

!
      end subroutine label
!
!
!
      subroutine linegraph(fieldtype)
      use feminterface, only: elemnt, getpostsetting, xy2lam,           &
     &                        convertcolorstring, getsetting, fieldquantity
      use femtypes
      use globalvariables
      implicit none
      character (len=*) :: fieldtype
      intent (in) :: fieldtype
!
!  local variables
      integer (I4B) :: i, ierror=0, div, elem, linewidth, oldindex, unitid
      real (SP) :: r = 0.0, g = 0.0, b = 0.0
      real (SP), allocatable :: distvalsp(:), zvalsp(:)
      real (DP) :: dist, delta, lambda(3), phi
      real (DP) :: startp(2), endp(2), vec(2), unitvec(2)
      real (DP), allocatable :: distval(:), points(:,:), zval(:)
      complex (DPC) :: fieldval
      character (len=10) :: tmp, unit
      character (len=14) :: logscale, output
      character (len=50) :: descriptor
      character (len=200) :: path
      logical :: ok, ok1
!------------------------------------------------------------------------------
!
      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgqwin, pgsvp, pgwnad, pgenv

      call getpostsetting('STARTX',startp(1))
      call getpostsetting('STARTY',startp(2))
      call getpostsetting('ENDX',endp(1))
      call getpostsetting('ENDY',endp(2))
      call getpostsetting('DIVISION',div)
      call getpostsetting('PHI',phi)
      call getpostsetting('LOGSCALE',logscale)
      call getpostsetting('LINEWIDTH',linewidth)
      call getpostsetting('OUTPUT',output)
!  compute direction vector from start to end
      vec = endp - startp
!  compute length of vector and unitvector (length 1)
      dist = sqrt(vec(1)**2 + vec(2)**2)
      unitvec = vec/dist
!  compute delta depending on # of divisions on length
      delta = dist/div
!  allocate and initialize arrays
      allocate(distval(div+1), points(2,div+1), zval(div+1))
      distval = 0._DP
      points = 0._DP
      zval = 0._DP
!  loop over all points
      do i = 1, div+1
!  compute x,y coordinates for new point depending on # of division
        distval(i) = (i-1)*delta
        points(:,i)= startp + (i-1)*delta*unitvec
!  get element for point
        call elemnt(xn,yn,e,n,points(1,i),points(2,i),elem,en,ok)
        if (ok) then
          call xy2lam(points(1,i),points(2,i),elem,lambda,xn,yn,e)
          call fieldquantity(elem,fieldtype,lambda,phi,fieldval,descriptor,unit,ok1)
          if (.not.ok1) then
            print *,"***** fieldtype: ",trim(fieldtype)," unknown for ",trim(physics)
            stop
          end if
          zval(i) = real(fieldval,DP)
        else
          print *," **** Error in lineplot:"
          print "(a,g10.3,a,g10.3,a)"," point (",points(1,i),";",points(2,i),") not in an element"
          cycle
        end if
      end do
!  plot the graph
      allocate(distvalsp(div+1), zvalsp(div+1))
      distvalsp = 0._SP
      zvalsp = 0._SP
      distvalsp = real(distval,SP)
!  if logarithmic scale should be shown, then modify the values
      select case (logscale)
      case ('NO_LOG')
        zvalsp = real(zval,SP)
      case ('LOG')
        zvalsp = real(log10(zval),SP)
      case ('DB')
        zvalsp = real(20*log10(zval),SP)
      end select
!  writing to file
      if(output .eq. 'YES') then
      call getsetting('PROJECTPATH',path)
        call grglun(unitid)
!  open a file in the specified project path
        open (unitid,FILE=path(1:len_trim(path))//'linegraph.txt',STATUS='NEW',&
     &        FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
        if (ierror .ne. 0) then
          do i = 1,99
            write (tmp ,'(i2)') i
            open (unitid,FILE=path(1:len_trim(path))//'linegraph'//trim(adjustl(tmp))//'.txt',&
     &            STATUS='NEW',FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
            if (ierror .ne. 0) then
              cycle
            else
              print '(a)',"Writing output to "//path(1:len_trim(path))//'linegraph'//trim(adjustl(tmp))//'.txt'
              exit
            end if
          end do
        else
          print '(a)',"Writing output to "//path(1:len_trim(path))//'linegraph.txt'
        end if
!  write the data to the file
        do i = 1, size(distval)
          write (unitid,*) distval(i), zval(i)
        end do
        close (unitid)
      end if
!  set environment for pgplot
      if (minval(zvalsp) .ne. maxval(zvalsp)) then
        call pgenv(minval(distvalsp),maxval(distvalsp),minval(zvalsp),maxval(zvalsp),0,1)
      else
        call pgenv(minval(distvalsp),maxval(distvalsp),0.9_SP*maxval(zvalsp),1.1_SP*maxval(zvalsp),0,1)
      end if
!  change width of the line
      call pgslw(linewidth)
!  get color of line from POSTsettings.txt and set new color for line
      call convertcolorstring(r, g, b)
      call pgsci(oldindex)
      call pgscr(21, r, g, b)
      call pgsci(21)
      call pgline(div+1,distvalsp,zvalsp)
!
!  reset linewidth and color and deallocate arrays
      call pgsci(oldindex)
      call pgslw(1)
      deallocate(distval,distvalsp,points,zval,zvalsp)
      end subroutine linegraph
!
!
!
      subroutine plotarrow(px,py,xa,xe,ya,ye)
      use feminterface, only: elemnt, kreipf, getpostsetting, fieldquantity, xy2lam, convertcolorstring
      use femtypes
      use globalvariables
      implicit none
      integer (I4B) px, py
      real (DP) xa, xe, ya, ye
      intent (in) :: px, py, xa, xe, ya, ye
!
!----------------------------------------------------------------------------
!  Fetches the position of the arrow and plots the arrow by calling kreipf
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
      real (DP) :: vx(px,py), vy(px,py)
      complex (DPC) :: zvec(2)
      logical :: ok, ok1
      character (len=100) :: fieldtype
      character (len=10) :: unit
      character (len=50) :: descriptor
!
      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgqwin, pgsvp, pgwnad

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
          else
!  fetch the components of the vector
            call xy2lam(xp(i,j),yp(i,j),ielem,lambda,xn,yn,e)
            call fieldquantity(ielem,fieldtype,lambda,phi,zvec,descriptor,unit,ok1)
            if (.not.ok1) then
              print *,"***** fieldtype: ",trim(fieldtype)," unknown for ",trim(physics)
              stop
            end if
            vx(i,j) = real(zvec(1),DP)
            vy(i,j) = real(zvec(2),DP)
          end if
        end do
      end do
!
!  determin maximum legth of arrows e.g. 90% of spacing
      maxlen = min((xe-xa)/px,(ye-ya)/py)*0.90_DP
      call convertcolorstring(r, g, b)
      call pgqci(oldindex)
      call pgqcr(21, ocr, ocg, ocb)
      call pgscr(21, r, g, b)
      call pgsci(21)
      acolor = 0
!
!  now draw the arrow plot
      call kreipf(px,py,xp,yp,maxlen,vx,vy,acolor)
!  reset color index
      call pgscr(21, ocr, ocg, ocb)
      call pgsci(oldindex)
      return
      end subroutine plotarrow
!
!
!
      subroutine findrange(xa,xe,ya,ye,amin,amax,fieldtype)
      use feminterface, only: elemnt, fieldquantity, getpostsetting, xy2lam
      use femtypes
      use globalvariables
      implicit none
      real (DP) :: amin, amax, xa, xe, ya, ye, lambda(3), phi
      character (len=*) :: fieldtype
      intent (in) ::  xa, xe, ya, ye, fieldtype
      intent (out) ::  amin, amax
!----------------------------------------------------------------------
!  This subroutine finds the range for the plots
!  px : no. of points in x direction
!  py : no. of points in y direction
!  xa, xe :  structure range in x direction (start and end point)
!  ya, ye :  structure range in y direction (start and end point)
!  amin, amax : plotting range ; minimum and maximum
!----------------------------------------------------------------------
!
!  local variables
      integer (I4B) :: i, j, k, elem, level
      real(DP) :: zs, da
      real(DP) :: xs, ys, zsum, z2sum, sigma_z, zaverage
      complex (DPC) :: fieldval
      character (len=10) :: unit
      character (len=50) :: descriptor
      logical :: ok, ok1
!------------------------------------------------------------------------------
!
      external :: grglun, pgbeg, pgbox, pgcurs, pgdraw, pgend, pglabel, pgline, pgsch
      external :: pgqci, pgscf, pgsci, pgscr, pgslw, pgqch, pgqcr, pgmove
      external :: pgopen, pgpap, pgqwin, pgsvp, pgwnad

!  getting the value of PHI from POSTsettings
      call getpostsetting('PHI',phi)
!
!  k is the number of valid values
      amax = -huge(1._DP)
      amin =  huge(1._DP)
      k = 0
      zsum=0._DP
      z2sum=0._DP
      do level=1,10
        do i=0,2**(level-1)-1
          xs = (xe-xa)*real(1+2*i,DP)/real(2**level,DP) + xa
          do j=0,2**(level-1)-1
            ys = (ye-ya)*real(1+2*j,DP)/real(2**level,DP) + ya 
            call elemnt(xn,yn,e,n,xs,ys,elem,en,ok)
            if (ok) then
              call xy2lam(xs,ys,elem,lambda,xn,yn,e)
              call fieldquantity(elem,fieldtype,lambda,phi,fieldval,descriptor,unit,ok1)
              if (.not.ok1) then
                print *,"***** fieldtype: ",trim(fieldtype)," unknown for ",trim(physics)
                stop
              end if
              zs = real(fieldval,DP)
              zsum = zsum+zs
              z2sum = z2sum+zs**2
              amax = max(amax,zs)
              amin = min(amin,zs)
              k = k+1
            end if
          end do
        end do
        if (k .ge. 2000) exit
      end do
!
!  some statistics: average value and standard deviation
      zaverage=zsum/real(k,DP)
      sigma_z=sqrt((z2sum-zaverage*zsum)/real(k-1,DP))
!     
      amax = min(amax,zaverage+5._DP*sigma_z)
      amin = max(amin,zaverage-5._DP*sigma_z)
!
      da=amax-amin
!  max is 1% larger than maximum value
      amax = amax + 0.01_DP*da
!  min is 1% smaller than minimum value
      amin = amin - 0.01_DP*da
!
      return
      end subroutine findrange
!
!
!
      end module postmodule