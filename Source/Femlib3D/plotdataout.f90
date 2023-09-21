      subroutine plotdataout(filename, plottype, rowsDescription, columnsDescription, plottitle, dat, append)
      use feminterface,        only: getsetting
      use femtypes
      implicit none
      character (len=256) :: filename, plottype, rowsDescription, columnsDescription, plottitle
      real (DP)           :: dat(:,:)
      logical             :: append
      intent (in)         :: filename, plottype, rowsDescription, columnsDescription, plottitle, dat, append

!
!------------------------------------------------------------------------------
!    $Revision: 1.2 $
!    $Date: 2015/03/11 17:03:25 $
!    $Author: juryanatzki $
!------------------------------------------------------------------------------
!
!  Write plotting data down into file. The file contains:
! 
! Header line Format for parsing in PolyDE: H#DataType;DataFileName;PlotType;PlotName
! Possible DataTypes: 
! Possible PlotTypes: XYPLOT, HISTOGRAM, FIELDPLOT, PLOTARRAY
! Examples:
!                             Preferred  |   Array      | Data Descript.|   Name of the Plot
!                                Plottype|      Shape   |    OPTIONAL   |
!    WRITE(unitid,1313) "HEADER;HISTOGRAM#",numv,";",4,"#Elements;Angles#Solid Angles of all Mesh Elements"
!                                        |              |               |
!
! Now writing of the actual dataFiles follows:
!  HEADER
!  row1data;row1data;row1data;row1data;    ...
!  row2data;row2data;row2data;row2data;    ...
!  row3data;row3data;row3data;    ...
!
!------------------------------------------------------------------------------
!  local variables
      integer (I4B)        :: i, unitid, ios, irow, icol, nrow, ncol
      character   (len=10) :: dataRangeForm, formfor_plottypeForm, plottypeForm, formfor_rowsDescriptionForm, rowsDescriptionForm, &
     &                       formfor_columnsDescriptionForm, columnsDescriptionForm, formfor_plottitleForm, plottitleForm, tmpFilenumber
      character   (len=30) :: stringData
      character  (len=100) :: headerForm
      character  (len=200) :: path, lineString
      integer (I4B)        :: ierror
!
      external    ::    grglun


!  get project path from user's environment variables
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)

!  open a file in the specified project path
      open (unitid,FILE=path(1:len_trim(path))//trim(filename)//'.dat',STATUS='NEW',&
     &      FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
      if (ierror .ne. 0) then
        if (append) then
!  reopen file to append data
          open (unitid,file=path(1:len_trim(path))//trim(filename)//'.dat', &
     &          form='formatted',position='APPEND',action='WRITE',iostat=ios)
          if (ios .ne. 0) then
            print*, '**** Error on opening file: ',path(1:len_trim(path))//trim(filename)//'.dat'
            print*, '**** IO Error No: ',ios
            return
          end if
        else
!  Create a new file
          tmpFilenumber = ''
          do i = 1,99
            write (tmpFilenumber ,'(i2)') i
            open (unitid,FILE=path(1:len_trim(path))//trim(filename)//trim(adjustl(tmpFilenumber))//'.dat',&
     &            STATUS='NEW',FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
            if (ierror .ne. 0) then
              cycle
            else
              print '(a)',"Writing output to "//path(1:len_trim(path))//trim(filename)//trim(adjustl(tmpFilenumber))//'.dat'
              exit
            end if
          end do
        end if
      else
        print '(a)',"Writing output to "//path(1:len_trim(path))//trim(filename)//'.dat'
      end if

      nrow = size(dat,2)
      ncol = size(dat,1)


! Determine the form of the desctiption of rows and colomns
      write(dataRangeForm,'(A1,I1,A5,I1)') "I", ceiling(log10(real(nrow))), ",A1,I", ceiling(log10(real(ncol)))

! Determine the form of the plottype description
      write(formfor_plottypeForm,'(A5,I1,A1)') "(A1,I", ceiling(log10(real(len_trim(plottype)))), ")"
      write(plottypeForm,trim(formfor_plottypeForm)) "A", len_trim(plottype)

! Determine the form of dataDescription for rows
      write(formfor_rowsDescriptionForm,'(A5,I1,A1)') "(A1,I", ceiling(log10(real(len_trim(rowsDescription)))), ")"
      write(rowsDescriptionForm,trim(formfor_rowsDescriptionForm)) "A", len_trim(rowsDescription)

! Determine the form of dataDescription for columns
      write(formfor_columnsDescriptionForm,'(A5,I1,A1)') "(A1,I", ceiling(log10(real(len_trim(columnsDescription)))), ")"
      write(columnsDescriptionForm,trim(formfor_columnsDescriptionForm)) "A", len_trim(columnsDescription)

! Determine the form of the plottitle
      write(formfor_plottitleForm,'(A5,I1,A1)') "(A1,I", ceiling(log10(real(len_trim(plottitle)))), ")"
      write(plottitleForm,trim(formfor_plottitleForm)) "A", len_trim(plottitle)

      headerForm = "(A6,A1,"//trim(plottypeForm)//",A1,"//trim(dataRangeForm)//",A1,"//trim(rowsDescriptionForm)//",A1,"//trim(columnsDescriptionForm)//",A1,"//trim(plottitleForm)//")"

! Write Header
      print*,'** Generating Plot Data for ', trim(filename)
      write(unitid,trim(headerForm)) "HEADER",";",plottype(1:len_trim(plottype)),"#",nrow,";",ncol,"#",trim(rowsDescription),";",trim(columnsDescription),"#",trim(plottitle)
! Write Data
      do irow = 1, nrow
        lineString = ''
        do icol = 1, ncol
          write (stringData,"(F)") dat(icol,irow)
          lineString = trim(lineString)//trim(adjustl(stringData))
          if (icol .eq. ncol ) then
            exit
          end if
          lineString = trim(lineString)//";"
        end do
        write(unitid,"(A)") lineString(1:len_trim(lineString))
      end do

      close(unitid)

      return
end subroutine plotdataout
