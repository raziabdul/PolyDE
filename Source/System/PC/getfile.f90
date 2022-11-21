      subroutine getfile(file,path)
      use feminterface, only:
      use femtypes
!  following libraries are for Compaq Visual Fortran
      use DFLIB
      use DFWIN
!  following libraries are for Intel Visual Fortran
!      use ifqwin
!      use ifwin
      implicit none
      character (len=*) file, path
      intent (out) file, path
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
!    $Revision: 1.10 $
!    $Date: 2008/12/22 13:15:12 $
!    $Author: m_kasper $
!
! 
!  ask user to input a filename 
!
!  MS Windows Version
!  opens a dialog box for specifying the filename
!
!  local variables 
!
      character (len=300) filename
      type (T_OPENFILENAME) ofns                   ! OPENFILENAME Structure
      character(len=20*5) Filter
      character(len=60) title

      Filter     = 'DXF Files(*.DXF)' // char(0) //                     &
     &             '*.dxf'            // char(0) //                     &
     &             'All Files(*.*)'   // char(0) //                     &
     &             '*.*'              // char(0) // char(0)

      title = 'Acadnet Input'

      ofns%hwndOwner         = GETHWNDQQ(QWIN$FRAMEWINDOW)  ! Handle to the window that owns the dialog box
      ofns%hInstance         = NULL                         ! If the OFN_ENABLETEMPLATEHANDLE flag is set in the Flags member, hInstance is a handle to a memory object containing a dialog box template
      ofns%lpstrFilter       = LOC(Filter)                  ! Pointer to a buffer containing pairs of null-terminated filter strings. The last string in the buffer must be terminated by two NULL characters. 
      ofns%lpstrCustomFilter = NULL                         ! Pointer to a static buffer that contains a pair of null-terminated filter strings for preserving the filter pattern chosen by the user
      ofns%nMaxCustFilter    = NULL                         ! Specifies the size, in TCHARs, of the buffer identified by lpstrCustomFilter
      ofns%nFilterIndex      = 1                            ! Specifies the index of the currently selected filter in the File Types control
      ofns%lpstrFile         = LOC(filename)                ! Pointer to a buffer that contains a file name used to initialize the File Name edit control
      ofns%nMaxFile          = LEN(filename)                ! Specifies the size, in TCHARs, of the buffer pointed to by lpstrFile
      ofns%lpstrFileTitle    = NULL                         ! Pointer to a buffer that receives the file name and extension (without path information) of the selected file
      ofns%nMaxFileTitle     = NULL                         ! Specifies the size, in TCHARs, of the buffer pointed to by lpstrFileTitle
      ofns%lpstrInitialDir   = NULL                         ! Pointer to a NULL terminated string that can specify the initial directory
      ofns%lpstrTitle        = LOC(title)                   ! Pointer to a string to be placed in the title bar of the dialog box
      ofns%Flags             = NULL                         ! A set of bit flags you can use to initialize the dialog box
      ofns%nFileOffset       = NULL                         ! Specifies the zero-based offset, in TCHARs, from the beginning of the path to the file name in the string pointed to by lpstrFile
      ofns%nFileExtension    = NULL                         ! Specifies the zero-based offset, in TCHARs, from the beginning of the path to the file name extension in the string pointed to by lpstrFile
      ofns%lpstrDefExt       = NULL                         ! Pointer to a buffer that contains the default extension
      ofns%lCustData         = NULL                         ! Specifies application-defined data that the system passes to the hook procedure identified by the lpfnHook member
      ofns%lpfnHook          = NULL                         ! Pointer to a hook procedure. This member is ignored unless the Flags member includes the OFN_ENABLEHOOK flag
      ofns%lpTemplateName    = NULL                         ! Pointer to a null-terminated string that names a dialog template resource in the module identified by the hInstance member
      ofns%pvReserved        = NULL                         ! Reserved. Must be set to NULL
      ofns%dwReserved        = 0                            ! Reserved. Must be set to 0.
      ofns%FlagsEx           = 0                            ! A set of bit flags you can use to initialize the dialog box

      ofns%lStructSize       = sizeof (ofns)                ! Specifies the length, in bytes, of the structure. Windows 2000/XP: Use sizeof (OPENFILENAME) for this parameter

      if (getopenfilename(ofns)) then
        file=filename(index(filename,'\',.true.)+1:)
        path=filename(1:index(filename,'\',.true.))
        print '(a)',trim(path)//trim(file)
        print*
      else
        print*,'user abort'
        stop
      end if

      return
      end subroutine getfile 