      program post
      use feminterface, only: initialize, getsetting, getpostsetting, readpostsetting, &
                              zeit, zanfpp, meshinput, lin, struc, triplt, residual,   &
                              fpotentialp, convertcolorstring, angles, farbpalette,    &
                              farbe, elemnt, xy2lam, field, numberpoly, fieldquantity, &
                              putpostsetting, scatterparameters, total_loss, linfluid, &
                              linstokes, energyflow, currentflow, plottensor
      use femtypes
      use globalvariables
      use mpmodule
      use postmodule
      implicit none
!
!------------------------------------------------------------------------------
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
!    $Revision: 1.128 $
!    $Date: 2015/04/01 10:54:20 $
!    $Author: juryanatzki $
!------------------------------------------------------------------------------
!
!  Postprocessor for Polyde
!
      integer (I4B) :: errcode, open_check
      integer (I4B) :: bildnr, elem, i, ios, logunit, linewidth
      integer (I4B) :: nc=-16, nl=-16, oldindex
      integer (I4B) :: gridx, gridy, nature, succ, pgcurs
      integer (I4B), allocatable :: palette(:)
      real (SP) :: r = 0.0, g = 0.0, b = 0.0
      real (SP) :: margin, ocr, ocg, ocb
      real (SP) :: x_temp, y_temp
      real (SP) :: rxmin, rxmax, rymin, rymax
      real (DP) :: jdmesh, jdsolu
      real (DP) :: amin, amax, phi
      real (DP) :: h, w1, w2, w3, zcolmin, zcolmax
      real (DP) :: epsgl, resgl, resnl, error
      real (DP), allocatable :: sumres(:), sumref(:)
      real (DP) :: xplmin, xplmax, yplmin, yplmax, dch
      real (DP) :: startx, starty, endx, endy, lambda(3)
      real (DP) :: flux, swr, rc, tc, totloss  ! for scattering parameters and losses
      real (DP) :: current
      real (DP) :: fmaxi, fmini, fmaxl, fminl
      real (DP), allocatable :: wmin(:), res(:,:)
      complex (DPC) :: fieldval
      logical :: ok, gilt, eleinfo, grafik 
      logical :: colorfill, ext, hlin, int
      logical :: reell, istlin, konfom, zyl
      character (len=1)   :: ch
      character (len=30)   :: colorscale
      character (len=3)   :: pnumber, information
      character (len=10)  :: datatype, unit
      character (len=14)  :: displaydevice, drawingtype, logscale, region
      character (len=16)  :: estimator
      character (len=30)  :: callback, fieldtype
      character (len=50)  :: descriptor
      character (len=200) :: path

!------------------------------------------------------------------------------
!
      external :: grglun, pgbox, pgcurs, pgend, pglabel, pgline, pgsch, pgqci, pgsci, pgscr, pgslw, pgqcr
      external :: pgpap, pgqwin, pgsvp, pgwnad
      call initialize()
!
!  read path and open file post.log
      call getsetting('PROJECTPATH',path)
      call grglun(logunit)
      open (logunit,file=path(1:len_trim(path))//'post.log',            &
     &     form='FORMATTED',position='REWIND',iostat=ios)
!
! read the mesh
      call meshinput(ok,jdmesh,reell,gilt,istlin,konfom,zyl,epsgl,resgl,resnl,error)
      print *,'Dimension of structure'
      print "(a,g10.3,a,g10.3,a,g10.3,a,g10.3)                              &
     &    "," xmin: ",xmin," xmax: ",xmax," ymin: ",ymin," ymax: ",ymax
      call zeit(' Reading the Mesh')
      grafik=.false.
!
!  read in solution from file
      call getsetting('PHYSICS_MODE',physics)
      if (physics .eq.'FLUID') then
        call linfluid(ok, jdmesh, jdsolu, gilt, eleinfo, n, ndof, omega,&
     &         epsgl, ep, eg, xfluid, fem_accuracy)
      else if (physics .eq.'STOKES') then
        call linstokes(ok, jdmesh, jdsolu, gilt, eleinfo, n, ndof, omega,&
     &         epsgl, ep, eg, xfluid, fem_accuracy)
      else
        call lin(ok, jdmesh, jdsolu, gilt, eleinfo, n, ndof, omega,     &
     &         epsgl, nnat, ep, eg, x, fem_accuracy)
      end if
!
!
!
!  check for use of multiphysics and read data if yes
      if ((physics .eq. 'HEATTR') .or. (physics .eq. 'TEWAVE') .or. (physics .eq. 'TMWAVE')) then
        call initdata()
      end if
!
!
! 
!  print information about solution
      print "(a,g10.3/,a,g10.3)","frequency        : ",omega/(2._DP*pi)&
                              ,"angular frequency: ",omega
      if (omega .gt. 0) then
        print "(a,g10.3)","vacuum wavelength: ",(c0 / omega)*2._DP*pi
      end if
      print *,"Plotting ",n," elements with ",ndof," degrees of freedom."
      print "(a,g10.3)"," Relative estimated error in energy norm: ",100._DP*fem_accuracy
!
!
!  forever loop to draw plots
!
!  A drawing is generated in this loop. You can select what to draw in the file
!  POSTsettings.txt, which should reside in your project-folder. There are diffe-
!  rent commands to choose from:
!      ARROWPLOT      draws arrowplot for gradients, inplane-field and Poynting
!                     vector (last two only available for TE-/TMWAVE mode)
!      END            plotting / postprocessorwill be stopped
!      EXPORTDATA     write values of quantity fieldtype for given grid to file
!                     data.txt
!      FIELD          draws color plot of field distribution
!                     (red is positive, blue is negative, green is zero)
!                     (FIELDTYPE, DRAWINGTYPE, NUMCOLORS, AMIN, AMAX and PHI will
!                     have to be specified before the call once. Else the standard
!                     values will be chosen.)
!      LINEGRAPH      draws a graph of chosen field quantity along a line defined by
!                     STARTX, STARTY, ENDX and ENDY
!      LOSSES         compute total EM-losses in the domain and EM-losses along
!                     a given line.
!      MESH           draws the mesh
!      MESHQUALITY    draws the mesh's quality depending on the element's inner angles 
!                     (blue is good, red is bad)
!      MOUSE          mouse input for two points can be obtained. use mouse cursor and
!                     space key to select point under mouse cursor. first point is
!                     (startx,starty), second (endx,endy)
!      OPEN           opens a new graphics window in post
!                     (the extents will have to be specified first)
!      POINTVALUE     gives potential value of given point (startx/starty)
!      POLYORDER      draws the polynomial order of all elements
!                     (blue is 1, red is 10)
!      RESIDUAL       draws the (composed) residual of all elements
!                     (blue is low, red is the maximum residual)
!      SCALE          draws box around geometry with size of extends and tick marks.
!                     unit is meters (SI-units)
!      SCATTER        scatter parameters (SWR, transmission and reflection
!                     coefficient) along a line
!      STRUCT         draws the structure
!
!
!  initialize open_check for OPEN command (PS-files)
      open_check = 0
      do
        call readpostsetting(callback)
!
        select case(callback)
!------------------------------------------------------------------------------
!  drawing arrowplot
        case ('ARROWPLOT')
          if (.not.gilt) then
            print *,"**** Error No valid solution for plotting"
            cycle
          end if
          call getpostsetting('GRIDX',gridx)
          call getpostsetting('GRIDY',gridy)
          call plotarrow(gridx,gridy,xplmin,xplmax,yplmin,yplmax)
!------------------------------------------------------------------------------
!  compute electric current through port
        case ('CURRENT')
          if (physics .ne. 'STATCURRENT') then
            print *,"  *** Physics Mode selected is ",physics
            print *,"  *** no current computation available for ",physics
            cycle
          end if
!  cycle, if no valid solution data present
          if (.not.gilt) then
            print *,"**** Error No valid solution for plotting"
            cycle
          end if
!  get and check coordinates for losses along a line
          call getpostsetting('STARTX', startx)
          call getpostsetting('STARTY', starty)
          if ((startx .lt. xplmin) .or. (startx .gt. xplmax) .or. &
              (starty .lt. yplmin) .or. (starty .gt. yplmax)) then
            print *,"Point not in calculated domain"
            cycle
          end if
          call getpostsetting('ENDX', endx)
          call getpostsetting('ENDY', endy)
          if ((endx .lt. xplmin) .or. (endx .gt. xplmax) .or. &
              (endy .lt. yplmin) .or. (endy .gt. yplmax)) then
            print *,"Point not in calculated domain"
            cycle
          end if
!  set line width
          call getpostsetting('LINEWIDTH',linewidth)
          call pgslw(linewidth)
!  get color of line from POSTsettings.txt and set new color for line
          call convertcolorstring(r, g, b)
          call pgsci(oldindex)
          call pgscr(21, r, g, b)
          call pgsci(21)
          call pgline(2,(/real(startx,SP),real(endx,SP)/),(/real(starty,SP),real(endy,SP)/))
          call currentflow(startx,starty,endx,endy,current)
!  reset linewidth and color and deallocate arrays
          call pgsci(oldindex)
          call pgslw(1)
!------------------------------------------------------------------------------
!  leave the program
        case('END')
          print *,"END-command used. Plotting stopped."
          exit
!------------------------------------------------------------------------------
!  End Of File is reached
        case('EOF')
          print *,"End Of File reached"
          exit
!------------------------------------------------------------------------------
!  an error occured during reading the settings in readpostsettings.f90
        case('ERROR')
          exit
!------------------------------------------------------------------------------
!  write field values of given fieldtype to file data.txt
        case ('EXPORTDATA')
          if (.not.gilt) then
            print *,"**** Error No valid solution for calculating value"
            cycle
          end if
!
!  get fieldtype / -quantity, number of gridpoints in x- and y-direction and
!  phase from POSTsettings.txt
          call getpostsetting('DATATYPE',datatype)
          call getpostsetting('FIELDTYPE',fieldtype)
          call getpostsetting('PHI',phi)
          if (datatype .eq. 'GRID') then
            call getpostsetting('GRIDX',gridx)
            call getpostsetting('GRIDY',gridy)
          end if
!  export grid values to file gridvalues.txt
          call exportdata(datatype,fieldtype,phi,xplmin,xplmax,yplmin,yplmax,gridx,gridy)
!------------------------------------------------------------------------------
!  draw color plot of field distribution
        case ('FIELD')
          if (.not.gilt) then
            print *,"**** Error No valid solution for plotting"
            cycle
          end if
!  get what to plot as field distribution from POSTsettings.txt
          call getpostsetting('FIELDTYPE',fieldtype)
          call getpostsetting('NATURE',nature)
          if (nature .gt. nnat) then 
            print*,'***Error in POSTSETTINGS: NATURE is assigned to:',nature 
            print*,'***The PHYSICS MODE has only ',nnat,' natures'
            stop
          end if
!  get bounds for color scale and phase from plot
          call getpostsetting('NUMCOLORS',nc)
          call getpostsetting('NUMLINES',nl)
          if (nc .lt. 0 .or. nl .lt. 0) then
!  finde range (of color scale) within the actual display region
            call findrange(xplmin,xplmax,yplmin,yplmax,fmini,fmaxi,fieldtype)
            fmaxl = fmaxi
            fminl = fmini
          else
            call getpostsetting('AMIN',amin)
            call getpostsetting('AMAX',amax)
            if (amax .lt. amin) then
              print*,'*** minimum value is larger than maximum value (Amin/Amax)'
              fmaxl = amin
              fminl = amax
            else
              fmaxl = amax
              fminl = amin
            end if
          end if
          call getpostsetting('PHI',phi)
!  set the type of drawing:
!  FULLPLOT     contains false color plot with equipotentials
!  COLORPLOT    contains false color plot only
!  LINEPLOT     contains equipotentials only
          call getpostsetting('DRAWINGTYPE',drawingtype)
          drawingtype = drawingtype(1:len_trim(drawingtype))
          select case (drawingtype)
            case ('FULLPLOT')
              colorfill = .true.
              hlin      = .true.
              call convertcolorstring(r, g, b)
              call getpostsetting('LINEWIDTH',linewidth)
              call pgslw(linewidth)
              call getpostsetting('PALETTE',colorscale)
              allocate(palette(abs(nc)))
              call farbpalette(palette,abs(nc),colorscale)
!  draw the color scale
              call scale(palette,nc,fminl,fmaxl)
              deallocate(palette)
            case ('COLORPLOT')
              colorfill = .true.
              hlin      = .false.
              call getpostsetting('PALETTE',colorscale)
              allocate(palette(abs(nc)))
              call farbpalette(palette,abs(nc),colorscale)
!  draw the color scale
              call scale(palette,nc,fminl,fmaxl)
              deallocate(palette)
            case ('LINEPLOT')
              colorfill = .false.
              hlin      = .true.
              call convertcolorstring(r, g, b)
              call getpostsetting('LINEWIDTH',linewidth)
              call pgslw(linewidth)
          end select
!  do actual filling of triangles
          call fpotentialp(nc,nl,phi,fmaxl,fminl,colorfill,hlin,r,g,b,fieldtype)
          call zeit(' Color Plot of Field Distribution')
          colorfill = .false.
          hlin      = .false.
!------------------------------------------------------------------------------
        case ('INPUTDISPLAYSIZE')
          call getpostsetting('DISPLAYDEVICE',displaydevice)
          if ((displaydevice(2:2) .eq. 'W') .or. &
              (displaydevice(2:2) .eq. 'X')) then
            print "(a)","setting new values for XMIN,YMIN,XMAX and YMAX"
            print "(a)","waiting to acquire lower left corner with mouse ..."
            succ = pgcurs(x_temp,y_temp,ch)
            xplmin = real(x_temp,DP)
            yplmin = real(y_temp,DP)
            call putpostsetting('XMIN',xplmin)
            call putpostsetting('YMIN',yplmin)
            print "(a/,g10.3,g10.3)","coordinates for lower left corner: ",xplmin,yplmin
            print "(a)","waiting to acquire upper right corner with mouse ..."
            succ = pgcurs(x_temp,y_temp,ch)
            xplmax = real(x_temp,DP)
            yplmax = real(y_temp,DP)
            call putpostsetting('XMAX',xplmax)
            call putpostsetting('YMAX',yplmax)
            print "(a/,g10.3,g10.3)","coordinates for upper right corner: ",xplmax,yplmax
          end if
!------------------------------------------------------------------------------
        case('LABEL')
!  title of plot
          if (.not.gilt) then
            print *,"**** Error No valid solution for plotting"
            cycle
          end if
          call getpostsetting('FIELDTYPE',fieldtype)
          call label(fieldtype)
!------------------------------------------------------------------------------
        case('LINEGRAPH')
          call getpostsetting('FIELDTYPE',fieldtype)
!  draw the graph
          call linegraph(fieldtype)
!------------------------------------------------------------------------------
        case ('LOSSES')
          if ((physics .ne. 'TEWAVE') .and. (physics .ne. 'TMWAVE')) then
            print *,"  *** Physics Mode selected is ",physics
            print *,"  *** no EM-losses available for ",physics
            cycle
          end if
!  cycle, if no valid solution data present
          if (.not.gilt) then
            print *,"**** Error No valid solution for plotting"
            cycle
          end if
!  compute and display total losses in domain
          call getpostsetting('REGION', region)
          call total_loss(region,totloss)
          print *,"Total electromagnetic losses in domain: ",totloss
!  get and check coordinates for losses along a line
          call getpostsetting('STARTX', startx)
          call getpostsetting('STARTY', starty)
          if ((startx .lt. xplmin) .or. (startx .gt. xplmax) .or. &
              (starty .lt. yplmin) .or. (starty .gt. yplmax)) then
            print *,"Point not in calculated domain"
            cycle
          end if
          call getpostsetting('ENDX', endx)
          call getpostsetting('ENDY', endy)
          if ((endx .lt. xplmin) .or. (endx .gt. xplmax) .or. &
              (endy .lt. yplmin) .or. (endy .gt. yplmax)) then
            print *,"Point not in calculated domain"
            cycle
          end if
!  set line width
          call getpostsetting('LINEWIDTH',linewidth)
          call pgslw(linewidth)
!  get color of line from POSTsettings.txt and set new color for line
          call convertcolorstring(r, g, b)
          call pgsci(oldindex)
          call pgscr(21, r, g, b)
          call pgsci(21)
          call pgline(2,(/real(startx,SP),real(endx,SP)/),(/real(starty,SP),real(endy,SP)/))
          call energyflow(startx,starty,endx,endy,flux)
!  reset linewidth and color and deallocate arrays
          call pgsci(oldindex)
          call pgslw(1)
!------------------------------------------------------------------------------
!  draw mesh
        case ('MESH')
          call convertcolorstring(r, g, b)
          call getpostsetting('LINEWIDTH',linewidth)
          call pgslw(linewidth)
          call triplt(e,en,xn,yn,n,r,g,b)
!  draw polynomial order of elements as text in element's incircle
          call getpostsetting('PNUMBER',pnumber)
          if (pnumber .eq. 'YES') then
            call getpostsetting('NATURE',nature)
            if (nature .gt. nnat) then 
              print*,'***Error in POSTSETTINGS: NATURE is assigned to:',nature 
              print*,'***The PHYSICS MODE has only ',nnat,' natures'
              stop
            end if
            call numberpoly(dch,xplmin,xplmax,yplmin,yplmax,nature)
          end if
          call zeit(' Plotting the Mesh')
!------------------------------------------------------------------------------
!  draw the quality of the mesh
        case ('MESHQUALITY')
          call pglabel('','','Mesh Quality')
          call getpostsetting('LINEWIDTH',linewidth)
          call pgslw(linewidth)
          allocate(wmin(n))
          do i=1,n
            call angles(i,e,xn,yn,w1,w2,w3)
            wmin(i)=min(w1,w2,w3)
          end do 
          call getpostsetting('PALETTE',colorscale)
          nc = 10
          allocate (palette(nc))
          call farbpalette(palette,nc,colorscale)
          zcolmin = pi / 3.
          zcolmax = 0.
          call farbe(wmin,xn,yn,e,n,palette,nc,zcolmin,zcolmax)
          if (n < 20.e3) then
            r = 0.00
            call triplt(e,en,xn,yn,n,r,g,b)
          end if
          call scale(palette,nc,zcolmin,zcolmax)
          deallocate (palette, wmin)
          call zeit(' Plotting the Quality of the Mesh')
!------------------------------------------------------------------------------
        case ('MOUSE')
          call getpostsetting('DISPLAYDEVICE',displaydevice)    
          if ((displaydevice(2:2) .eq. 'W') .or.                        &
     &        (displaydevice(2:2) .eq. 'X')) then        
            print "(a)","waiting to acquire starting point with mouse ..."
            succ = pgcurs(x_temp,y_temp,ch)
            startx = real(x_temp,DP)
            starty = real(y_temp,DP)
            call putpostsetting('STARTX',startx)
            call putpostsetting('STARTY',starty)
            print "(a/,g10.3,g10.3)","coordinates for starting point: ",startx,starty
            print "(a)","waiting to acquire end point with mouse ..."
            succ = pgcurs(x_temp,y_temp,ch)
            endx = real(x_temp,DP)
            endy = real(y_temp,DP)
            call putpostsetting('ENDX',endx)
            call putpostsetting('ENDY',endy)
            print "(a/,g10.3,g10.3)","coordinates for end point: ",endx,endy
          end if
!------------------------------------------------------------------------------
! opens a new graphics window in post
        case ('OPEN')
          call getpostsetting('DISPLAYDEVICE',displaydevice)
          displaydevice = displaydevice(1:len_trim(displaydevice))
          call getpostsetting('XMIN',xplmin)
          call getpostsetting('XMAX',xplmax)
          call getpostsetting('YMIN',yplmin)
          call getpostsetting('YMAX',yplmax)
          call getpostsetting('INFO',information)
!  setting the paper size and view port with regards to the orientation
          if((open_check .gt. 0) .and.  ((displaydevice .eq. '/CPS')    &
     &        .or. (displaydevice .eq. '/PS')                           &
     &        .or. (displaydevice .eq. '/VCPS')                         &
     &        .or. (displaydevice .eq. '/VPS'))) then
            margin = 0.08
            if((displaydevice .eq. '/VCPS') .or. (displaydevice .eq. '/VPS')) then
              call pgpap(7.7,1.43)
              if(information .eq. 'YES') then
                call pgsvp(margin, 1.-margin, 0.2, 1.-margin)
              else
                call pgsvp(margin, 1.-margin-0.05, margin, 1.-margin)
              end if
            else
              call pgpap(11.0,0.7)
              if(information .eq. 'YES') then
                call pgsvp(margin, 0.8-0.02, margin, 1.-margin-0.03)
              else
                call pgsvp(margin, 1.-margin-0.05, margin+0.08, 1.-margin-0.03)
              end if
            end if 
            rxmin=real(xplmin)
            rxmax=real(xplmax)
            rymin=real(yplmin)
            rymax=real(yplmax)
            call pgwnad (rxmin,rxmax,rymin,rymax)
            call pgqwin(rxmin,rxmax,rymin,rymax)
            open_check = open_check + 1
!  open new window if no window was opened before
          else
!  compute default character heigth dch from extents
            dch = max(xplmax-xplmin,yplmax-yplmin)/40._DP
            if (xplmax .le. xplmin) then
              xplmin=xmin
              xplmax=xmax
            end if
            if (yplmax .le. yplmin) then
              yplmin=ymin
              yplmax=ymax
            end if
            call zanfpp(displaydevice,xplmin,xplmax,yplmin,yplmax,h,bildnr,information)
            grafik = .true.
            if ((displaydevice .eq. '/CPS') .or. (displaydevice .eq. '/PS') &
                .or. (displaydevice .eq. '/VCPS') .or. (displaydevice .eq. '/VPS')) then
              open_check = open_check + 1
            else
              open_check = 0
            end if
          endif
!------------------------------------------------------------------------------
!  display field value for given coordinates
        case ('POINTVALUE')
          if (.not.gilt) then
            print *,"**** Error No valid solution for calculating value"
            cycle
          end if
          call getpostsetting('FIELDTYPE',fieldtype)
          call getpostsetting('STARTX', startx)
          call getpostsetting('STARTY', starty)
          call getpostsetting('PHI',phi)
          if ((startx .lt. xmin) .or. (startx .gt. xmax) .or. &
              (starty .lt. ymin) .or. (starty .gt. ymax)) then
            print *,"Point is located outside limits"
            cycle
          end if
          call elemnt(xn,yn,e,n,startx,starty,elem,en,ok)
          if (.not.ok) then
            print *,"Point is not within the computational domain"
            cycle
          end if
!  get barycentric coordinates lambda of point in the old mesh
          call xy2lam(startx,starty,elem,lambda,xn,yn,e)
          call fieldquantity(elem,fieldtype,lambda,phi,fieldval,descriptor,unit,ok)
          if (ok) then
            print '(a,es10.3,a,es10.3,a)','field value at point (',startx,',',starty,') is:'
            print *,'real part     : ',real(fieldval,DP)
            print *,'imaginary part: ',aimag(fieldval)
            print *,'absolute value: ',abs(fieldval)
          else
            print *,"***** fieldtype: ",trim(fieldtype)," unknown for ",trim(physics)
          end if
!------------------------------------------------------------------------------
!  draw the polynomial order of the elements
        case ('POLYORDER')
          if (.not.gilt) then
            print *,"**** Error No valid solution for plotting"
            cycle
          end if
          call pglabel('','','Color Plot of Polynomial Order')
          call getpostsetting('NATURE',nature)
          if (nature .gt. nnat) then 
            print*,'***Error in POSTSETTINGS: NATURE is assigned to:',nature 
            print*,'***The PHYSICS MODE has only ',nnat,' natures'
            stop
          end if
          nc = polymax
          call convertcolorstring(r, g, b)
          call getpostsetting('PALETTE',colorscale)
          allocate (palette(nc))
          call farbpalette(palette,nc,colorscale)
          zcolmin = 0
          zcolmax = polymax
          call farbe(real(ep(:,nature),8),xn,yn,e,n,palette,nc,zcolmin,zcolmax)
!  draw polynomial order of elements as text in element's incircle
          call scale(palette,nc,zcolmin,zcolmax)
          deallocate(palette)
          call getpostsetting('PNUMBER',pnumber)
          if (pnumber .eq. 'YES') then
            call numberpoly(dch,xplmin,xplmax,yplmin,yplmax,nature)
          end if
          call zeit(' Color Plot of Polynomial Order for Elements')
!------------------------------------------------------------------------------
!  draw the residual of the solution
        case ('RESIDUAL')
          if (.not.gilt) then
            print *,"**** Error No valid solution for plotting"
            cycle
          end if
! compute error for solution
          call getsetting('ERROR_ESTIMATOR',estimator)
          estimator = estimator(1:len_trim(estimator))
          select case (estimator)
          case ('BOUNDARYRESIDUAL')
            ext = .true.
            int = .false.
          case ('INTERIORRESIDUAL')
            ext = .false.
            int = .true.
          case ('FULLRESIDUAL')
            ext = .true.
            int = .true.
          case default
            print "(a)","No specific error estimator chosen. Using defualt (full residual)."
            ext = .true.
            int = .true.
          end select
          allocate(res(n,nnat))
! temporary neglect for FLUID and probably STOKES
          if (physics .ne.'FLUID') then
          allocate(sumres(nnat),sumref(nnat))
          call residual(ext,int,.false.,errcode,res,sumres,sumref)
          deallocate(sumres,sumref)
          end if
          call getpostsetting('NATURE',nature)
          if (nature .gt. nnat) then 
            print*,'***Error in POSTSETTINGS: NATURE is assigned to:',nature 
            print*,'***The PHYSICS MODE has only ',nnat,' natures'
            stop
          end if
!  select how to calculate the residual
          call convertcolorstring(r, g, b)
          call getpostsetting('NUMCOLORS',nc)
          allocate (palette(abs(nc)))
          call getpostsetting('PALETTE',colorscale)
          call farbpalette(palette,abs(nc),colorscale)
          call getpostsetting('LOGSCALE',logscale)
          logscale = logscale(1:len_trim(logscale))
          select case (logscale)
            case ('NO_LOG')
              zcolmin = sum(sqrt(res(:,nature)))/n
              zcolmax = sqrt(maxval(res(:,nature)))
              print "(a)","plotting residual from mean to max value"
              print "(a,g10.3,a,g10.3)","mean value: ",zcolmin," max value: ",zcolmax
              call pglabel('','','Color Plot of Residual')
              call farbe(sqrt(res(:,nature)),xn,yn,e,n,palette,abs(nc),zcolmin,zcolmax)
              call scale(palette,nc,zcolmin,zcolmax)
            case ('LOG')
              zcolmin = log10(sum(sqrt(res(:,nature)))/n)
              zcolmax = log10(sqrt(maxval(res(:,nature))))
              print "(a)","plotting residual from mean to max value in logarithmic scale (base 10)"
              print "(a,g10.3,a,g10.3)","mean value: ",zcolmin," max value: ",zcolmax
              call pglabel('','','Color Plot of Residual (log10)')
              call farbe(log10(sqrt(res(:,nature))),xn,yn,e,n,palette,abs(nc),zcolmin,zcolmax)
              call scale(palette,nc,zcolmin,zcolmax)
          end select
          if (n < 1.e3) then
            call triplt(e,en,xn,yn,n,r,g,b)
          end if
          deallocate (palette,res)
          call zeit(' Color Plot of Residual for Elements')
!------------------------------------------------------------------------------
!  draw box with dimensions around structure and information
        case ('SCALE')
          if (grafik) then
            call convertcolorstring(r, g, b)
            call pgqci(oldindex)
            call pgqcr(21, ocr, ocg, ocb)
            call pgscr(21, r, g, b)
            call pgsci(21)
            call pgsch(0.7)
            call pgbox('BCTN',0.0,0,'BCTNV',0.0,0)
            call info
            call pgscr(21, ocr, ocg, ocb)
            call pgsci(oldindex)
            call pgsch(1.0)
          end if
!------------------------------------------------------------------------------
        case ('SCATTER')
          if ((physics .ne. 'TEWAVE') .and. (physics .ne. 'TMWAVE')) then
            print *,"  *** Physics Mode selected is ",physics
            print *,"  *** no scattering parameters available for ",physics
            cycle
          end if
          if (.not.gilt) then
            print *,"**** Error No valid solution for calculating value"
            cycle
          end if
          call getpostsetting('STARTX', startx)
          call getpostsetting('STARTY', starty)
          if ((startx .lt. xplmin) .or. (startx .gt. xplmax) .or. &
              (starty .lt. yplmin) .or. (starty .gt. yplmax)) then
            print *,"Point not in calculated domain"
            cycle
          end if
          call getpostsetting('ENDX', endx)
          call getpostsetting('ENDY', endy)
          if ((endx .lt. xplmin) .or. (endx .gt. xplmax) .or. &
              (endy .lt. yplmin) .or. (endy .gt. yplmax)) then
            print *,"Point not in calculated domain"
            cycle
          end if
          call scatterparameters(startx,starty,endx,endy,swr,rc,tc)
          print *,"Standing Wave Ratio:      ",swr
          print *,"Reflection Coefficient:   ",rc
          print *,"Transmission Coefficient: ",tc
!------------------------------------------------------------------------------
!  draw structure
        case ('STRUCT')
          call convertcolorstring(r, g, b)
          call getpostsetting('LINEWIDTH',linewidth)
          call pgslw(linewidth)
          call struc(gzz,zki,zrb,xbk,ybk,r,g,b)
          call zeit(' Plotting the Structure')
!------------------------------------------------------------------------------
!  drawing tensors on a grid
        case ('TENSORPLOT')
          if (.not.gilt) then
            print *,"**** Error No valid solution for plotting"
            cycle
          end if
          call getpostsetting('GRIDX',gridx)
          call getpostsetting('GRIDY',gridy)
          call plottensor(gridx,gridy,xplmin,xplmax,yplmin,yplmax)
!------------------------------------------------------------------------------
        case default
          print*,'*** unrecongized drawing command (ignored):  ', trim(callback)
        end select
      end do
!
      if (grafik) call pgend
      close(logunit)
!
      end program post
