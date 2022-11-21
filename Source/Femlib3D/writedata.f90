      subroutine writedata(keep_olddata,filename,plottype,plottitle,xunits,yunits,xdata,ydata)
      use feminterface, only: getsetting
      use femtypes
      implicit none
      integer   (I4B)     :: keep_olddata
      character (len=*)   :: filename,plottype,plottitle,xunits,yunits
      real      (DP)      :: xdata(:),ydata(:)
      intent    (in)      :: filename,plottype,plottitle,xunits,yunits,xdata,ydata, keep_olddata
!
!------------------------------------------------------------------------------
!    $Revision: 1.4 $
!    $Date: 2015/06/17 17:32:03 $
!    $Author: juryanatzki $
!------------------------------------------------------------------------------
!
! Header line Format for parsing in PolyDE: H#DataType;DataFileName;PlotType;PlotName
! Possible PlotTypes: XYPLOT, HISTOGRAM(, PLOTARRAY <- to be defined and implemented)
! Example:
!                               Preferred|     Array    | Data Descript.|   Name of the Plot
!                                Plottype|     Shape    |    OPTIONAL   |
!    WRITE(unitid,1313) "HEADER;HISTOGRAM#",numv,";",4,"#Elements;Angles#Solid Angles of all Mesh Elements"
!                                        |              |               |
!
!------------------------------------------------------------------------------
!  Input:
!       filename = Desired filename to which .dat will be added
!       plottype = Predefinition of the plottype
!      plottitle = Name of the Graph
!         xunits = Units of x-Axis
!         yunits = Units of y-Axis
!          xdata = Vector of x Data (should be of the size of ydata)
!          ydata = Vector of y Data (should be of the size of xdata)
!
!  local variables
      integer (I4B) i, unitid, ios
      character (len=200) :: path, plotsize
!
! Handle the case of different sizes of arrays xdata and ydata
    if (size(xdata).ne.size(ydata)) then
      print*, '**** Plotting Error: Not possible to handle different xdata and ydata vector sizes'
      return
    end if
! Get project path from user's environment variables
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)

!################# +       DATA WRITING       + #################
!################# +          START           + #################

1001    FORMAT (i)
1002    FORMAT (f)
1022    FORMAT (f,a1,f)
1003    FORMAT (a)
12222   FORMAT (f,f,f,f)
1013    FORMAT (i,a)
1313    FORMAT (A,I,A1,I1,A)

! Write MeshStats
    IF (keep_olddata.gt.0) THEN
          OPEN (unitid,file=path(1:LEN_TRIM(path))//filename(1:LEN_TRIM(filename)), &
    &     form='formatted',position='APPEND',action='WRITE',iostat=ios)
    ELSE
          OPEN (unitid,file=path(1:LEN_TRIM(path))//filename(1:LEN_TRIM(filename)), &
    &     form='formatted',position='REWIND',action='WRITE',iostat=ios)
    END IF

    IF (ios .ne. 0) THEN
        PRINT*, '**** Error on opening file: ',path(1:LEN_TRIM(path))//filename(1:LEN_TRIM(filename))
        PRINT*, '**** IO Error No: ',ios
        RETURN
    END IF
    
    WRITE(plotsize,"(I)") SIZE(xdata)
    ! Move text to the left
    plotsize = ADJUSTL(plotsize)
    WRITE(unitid,1003) "HEADER;"//plottype(1:LEN_TRIM(plottype))//"#"//plotsize(1:LEN_TRIM(plotsize))//";"//plotsize(1:LEN_TRIM(plotsize))//"#"//xunits(1:LEN_TRIM(xunits))//";"//yunits(1:LEN_TRIM(yunits))//"#"//plottitle(1:LEN_TRIM(plottitle))
    
    DO i = 1, SIZE(xdata)
      WRITE(unitid,1022) xdata(i), ";", ydata(i)
    END DO
    CLOSE(unitid)
    
!
!################# +          END OF          + #################
!################# +       DATA WRITING       + #################

      return
      end subroutine writedata
