      subroutine header(accur,xmin,xmax,ymin,ymax,unitfac)
      use feminterface, only: readky, getxyz, juld2cal, low2hi
      use femtypes
      implicit none
      real (DP) :: accur, xmax, xmin, ymax, ymin, unitfac
      intent (out) :: accur, xmax, xmin, ymax, ymin, unitfac
!
!    $Revision: 1.10 $
!    $Date: 2008/12/22 15:47:23 $
!    $Author: m_kasper $
!
!  read the header section of a dxf files 
!  and extract extends and the precision of the drawing 
!
!  Input:                 none (dxf file)
!  Output:   accur        precision (absolut)
!            xmin,ymin    lower left corner
!            xmax,ymax    upper right corner
!            unitfac      factor to convert to meters
!
!  known bugs:
!  - the AUTOCAD precision is not realy the accuracy of the drawing 
!    coodinate values, but is the number of digits in the display
!    of AUTOCAD drawing units. This means it is not really what we 
!    intent to have: the snap distance
!  - DXF Header variable $INSUNITS is used for evaluating unifrac.
!    However, this variable is not set appropriately in the DXF-file
!    and seems to be intended for units of bloc insertion.   
!    The dxf file does not contain other information about the units 
!    of the drawing.  
!
!  local variables
      integer (I4B) :: int, key, line, prec
      integer :: year, mon, day, hour, minute, second
      integer (I4B), parameter :: precmi=5
      real (DP) :: float, zmax, zmin, xsnap, ysnap, zsnap
      character (len=1000) :: string
      character (len=20) :: unit
      character (len=30) :: datstr
      character (len=3) :: dayofweek, month
!
      prec=precmi
      xmin=0._DP
      xmax=1._DP
      ymin=0._DP
      ymax=1._DP
      unitfac=1.e-3_DP
      unit='unitless'
!
100   call readky(key,string,float,int,line)
      if ((key.eq.0) .and. (string .eq. 'ENDSEC')) then
        write(*,222) xmin,ymin,unit
222     format(1x,'lower-left  drawing extends (x, y): (',2g10.3,') ',a)
        write(*,333) xmax,ymax,unit
333     format(1x,'upper-right drawing extends (x, y): (',2g10.3,') ',a)
        call low2hi(unit,len_trim(unit))
        write(*,444) unit(1:len_trim(unit)),unitfac
444     format(1x,'units are: ',a,' -> ',en11.3,' METERS')
!  estimate for the precision
        accur=min( (xmax-xmin) , (ymax-ymin) ) * 10._DP ** (-prec)
        return
      end if
      if (key.eq.9) then
        select case (string)
        case('$EXTMAX')
!  extends  - max
          call getxyz(xmax,ymax,zmax)
        case('$EXTMIN')
!  extends  - min
          call getxyz(xmin,ymin,zmin)
        case('$LUNITS')
        case('$LUPREC')
          call readky(key,string,float,int,line)
          if (key.eq.70) then
!  precision: the number of digits
            prec=max(int,precmi)
            write(*,111) prec,precmi
111         format(1x,'Precision is set to:',i3,' digits (minimum',i3,')')
          else
            print*,'*ERROR: syntax error'
          end if
        case('$SNAPUNIT')
          call getxyz(xsnap,ysnap,zsnap)
          print*,'Snap-Werte',xsnap,ysnap
        case('$INSUNITS')
          call readky(key,string,float,int,line)
          if (key.eq.70) then
            select case (int)
            case (0)
!  unitless
              unit='unitless'
              unitfac=1.e-3_DP
            case (1)
!  Inches
              unit='inches'
              unitfac=25.4e-3_DP
            case (2)
!  Feet
              unit='feet'
              unitfac=0.3048_DP
            case (3)
!  Miles
              unit='miles'
              unitfac=1609.344_DP
            case (4)
!  Millimeters
              unit='millimeters'
              unitfac=1.e-3_DP
            case (5)
!  Centimeters
              unit='centimeters'
              unitfac=1.e-2_DP
            case (6)
!  Meters
              unit='meters'
              unitfac=1._DP
            case (7)
!  Kilometers
              unit='kilometers'
              unitfac=1.e+3_DP
            case (8)
!  Microinches
              unit='microinches'
              unitfac=25.4e-9_DP
            case (9)
!  Mils
              unit='mils'
              unitfac=25.4e-6_DP
            case (10)
!  Yards
              unit='yards'
              unitfac=0.9144_DP
            case (11)
!  Angstroms
              unit='angstroms'
              unitfac=1.e-10_DP
            case (12)
!  Nanometers
              unit='nanometers'
              unitfac=1.e-9_DP
            case (13)
!  Microns
              unit='microns'
              unitfac=1.e-6_DP
            case (14)
!  Decimeters
              unit='decimeters'
              unitfac=1.e-1_DP
            case (15)
!  Decameters
              unit='decameters'
              unitfac=1.e1_DP
            case (16)
!  Hecometers
              unit='hectometers'
              unitfac=1.e2_DP
            case (17)
!  Gigameters
              unit='gigameters'
              unitfac=1.e9_DP
            case (18)
!  Astronomical units
              unit='astronomical units'
              unitfac=149597870691_DP
            case (19)
!  Light years
              unit='light yeard'
              unitfac=9460528405e6_DP
            case (20)
!  Parsecs
              unit='parsecs'
              unitfac=30856776e9_DP
            case default
            end select
          else
            print*,'*ERROR: syntax error'
          end if
        case ('$TDUPDATE')
          call readky(key,string,float,int,line)
          if (key.eq.40) then
            call juld2cal(float,year,mon,day,month,dayofweek,hour,minute,second)
            write(datstr,'(A3,'', '',I2.2,''. '',A3,''. '',I4,          &
     &        '' at '',I2,'':'',I2.2,'':''I2.2)')                       &
     &        dayofweek,day,month,year,hour,minute,second
            print*,'last update of the DXF-file:  ',datstr
          else
            print*,'*ERROR: syntax error'
          end if
        end select
      end if
!  ingnore all other header token
      goto 100
      end subroutine header
