      subroutine strtok(srce, delim, token)
      use feminterface, only:
      use femtypes
      implicit none
      character (len=*) :: token, srce, delim
      intent (in) :: srce, delim
      intent (out) :: token
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
!    $Revision: 1.9 $
!    $Date: 2014/07/15 13:04:51 $
!    $Author: m_kasper $
!
!     Tokenize a string in a similar manner to C.  The usage is
!     very similar to the C function.
!
!     Input:  Srce    =   Source string to tokenize. (see usage note)
!             Delim   =   Delimiter string.  Used to determine the
!                         beginning/end of each token in a string.
!
!     Output: Token  =    String Token.
!
!     Usage:  a) First Call Strtok with the string to tokenize as srce,
!                and the delimiter string used to tokenize scre is in delim
!             b) supsequent calls with srce=char(0)
!
!     CHARACTER*80 Token,SOURCE,DELIM,C
!          .
!          .      
!     SOURCE='This is a test. I hope that it Works! "eh" '
!     DELIM=' ,.;:"{}()!@#$%^&*'
!     call StrTok(SOURCE,DELIM,C)
!     PRINT *,C,' is the first token'
!     DO WHILE (C .ne. char(0))
!     PRINT *, C
!     call StrTok(char(0),DELIM,C)
!     ENDDO
!       
!
      character (len=1000) :: saves
      integer(I4B) :: start, length, begpos, endpos
      save :: start, saves
!
      if (srce(1:1) .ne. char(0)) then
        start = 1
        saves = srce
      endif
!
      begpos = start
      length = len_trim(saves)
5     continue
      if (begpos .le. length)  then 
        if (index(delim,saves(begpos:begpos)) .ne. 0) then
          begpos = begpos + 1
          goto 5
        end if
      endif
      if (begpos .gt. length) then
        token = char(0)
        return
      end if
      endpos = begpos
10    continue
      if (endpos .le. length) then 
        if (index(delim,saves(endpos:endpos)) .eq. 0) then
          endpos = endpos + 1
          goto 10
        end if
      endif
      token = saves(begpos:endpos-1)
      start=endpos
      return
      end subroutine strtok
!
!
!
      subroutine strtok1(srce, delim, token)
      use feminterface, only:
      use femtypes
      implicit none
      character (len=*) :: token, srce, delim
      intent (in) :: srce, delim
      intent (out) :: token
!
!     Tokenize a string in a similar manner to C.  The usage is
!     very similar to the C function.
!
!     Input:  Srce    =   Source string to tokenize. (see usage note)
!             Delim   =   Delimiter string.  Used to determine the
!                         beginning/end of each token in a string.
!
!     Output: Token  =    String Token.
!
!     Usage:  a) First Call Strtok with the string to tokenize as srce,
!                and the delimiter string used to tokenize scre is in delim
!             b) supsequent calls with srce=char(0)
!
!     CHARACTER*80 Token,SOURCE,DELIM,C
!          .
!          .      
!     SOURCE='This is a test. I hope that it Works! "eh" '
!     DELIM=' ,.;:"{}()!@#$%^&*'
!     call StrTok(SOURCE,DELIM,C)
!     PRINT *,C,' is the first token'
!     DO WHILE (C .ne. char(0))
!     PRINT *, C
!     call StrTok(char(0),DELIM,C)
!     ENDDO
!       
!
      character (len=1000) :: saves
      integer (I4B) :: start, begpos, endpos
      save :: start, saves
!
      if (srce(1:1) .ne. char(0)) then
        start = 1
        saves = srce
      endif
!
      begpos = start
      endpos=scan(saves(begpos:),delim)
      if (endpos .eq. 0) then 
        token = char(0)
      else 
        token = saves(begpos:endpos-1)
      end if
      start = endpos
!
      return
      end subroutine strtok1
!
!
!
      pure subroutine low2hi(string,length)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: length
      character (len=*) :: string
      intent (in) :: length
      intent(inout) :: string
!
!  Convert lower-case characters to UPPER-CASE
!
!  Input: 
!          length    upper limit of the text length up to which conversion 
!                    proceeds  
!  In-/ Output: 
!          string    text to be converted
!
!
      integer (I4B) :: i, cnumb
!
      do i=1,min(len_trim(string),length)
        cnumb=ichar(string(i:i))
        if ( cnumb.ge.97 .and. cnumb.le.122 ) then
          string(i:i)=char(cnumb-32)
        endif
      end do
      return
      end subroutine low2hi
!
!
!
      pure subroutine hi2low(string,length)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: length
      character (len=*) :: string
      intent (in) :: length
      intent(inout) :: string
!
!  Convert UPPER-CASE characters to lower-case
!
!  Input: 
!          length    upper limit of the text length up to which converion 
!                    proceeds  
!  In-/ Output: 
!          string    text to be converted
!
!
      integer (I4B) :: i,cnumb
!
      do i=1,min(len_trim(string),length)
        cnumb=ichar(string(i:i))
        if ( cnumb.ge.65 .and. cnumb.le.90 ) then
          string(i:i)=char(cnumb+32)
        endif
      end do
      return
      end subroutine hi2low
!
!
!
      subroutine string2number_SP(string,number,ok)
      use feminterface, only:
      use femtypes
      implicit none
      character(len=*) :: string
      real (SP) :: number
      logical, optional :: ok
      intent (in) :: string
      intent (out) :: number, ok
!  convert a string to a single precision number
      integer (I4B) :: ios
!
      read(string,*,iostat=ios) number
      if (ios .ne. 0) then
        if (present(ok)) then
          ok=.false.
        else
          print*,'Conversion Error in string2number (SP): ',string
        end if
      else if (present(ok)) then
        ok = .true.
      end if
      return
    end subroutine string2number_SP
!
!
!
    subroutine number2string_SP(number,string,ok)
      use femtypes
      implicit none
      character(len=*)  :: string
      real (SP)         :: number
      logical, optional :: ok
      intent (in)       :: number
      intent (out)      :: string, ok
!  convert a single precision number to a string
      integer (I4B)     :: ios
!
      write(string,*,iostat=ios) number
      if (ios .ne. 0) then
        if (present(ok)) then
          ok=.false.
        else
          print*,'Conversion Error in number2string (SP): ',number
        end if
      else if (present(ok)) then
        ok = .true.
      end if
      return
      end subroutine number2string_SP
!
!
!
      subroutine string2number_DP(string,number,ok)
      use feminterface, only:
      use femtypes
      implicit none
      character(len=*) :: string
      real (DP) :: number
      logical, optional :: ok
      intent (in) :: string
      intent (out) :: number, ok
!  convert a string to a double precision number
      integer (I4B) :: ios
!
      read(string,*,iostat=ios) number
      if (ios .ne. 0) then
        if (present(ok)) then
          ok=.false.
        else
          print*,'Conversion Error in string2number (DP): ',string
          print*,'Number Recognized (DP): ',number
          pause
        end if
      else if (present(ok)) then
        ok = .true.
      end if
      return
    end subroutine string2number_DP
!
!
!
    subroutine number2string_DP(number,string,ok)
      use femtypes
      implicit none
      character(len=*)  :: string
      real (DP)         :: number
      logical, optional :: ok
      intent (in)       :: number
      intent (out)      :: string, ok
!  convert a double precision number to a string
      integer (I4B)     :: ios
!
      write(string,*,iostat=ios) number
      if (ios .ne. 0) then
        if (present(ok)) then
          ok=.false.
        else
          print*,'Conversion Error in number2string (DP): ',number
        end if
      else if (present(ok)) then
        ok = .true.
      end if
      return
      end subroutine number2string_DP
!
!
!
      subroutine string2number_I(string,number,ok)
      use feminterface, only:
      use femtypes
      implicit none
      character(len=*) :: string
      integer (I4B) :: number
      logical, optional :: ok
      intent (in) :: string
      intent (out) :: number, ok
!  convert a string to a integer number
      integer (I4B) :: ios
!
      read(string,*,iostat=ios) number
      if (ios .ne. 0) then
        if (present(ok)) then
          ok=.false.
        else
          print*,'Conversion Error in string2number (I): ',string
        end if
      else if (present(ok)) then
        ok = .true.
      end if
      return
    end subroutine string2number_I
!
!
!
    subroutine number2string_I(number,string,ok)
      use femtypes
      implicit none
      character(len=*)  :: string
      integer (I4B)     :: number
      logical, optional :: ok
      intent (in)       :: number
      intent (out)      :: string, ok
!  convert a integer number to a string
      integer (I4B)     :: ios
!
      write(string,*,iostat=ios) number
      if (ios .ne. 0) then
        if (present(ok)) then
          ok=.false.
        else
          print*,'Conversion Error in number2string (I): ',number
        end if
      else if (present(ok)) then
        ok = .true.
      end if
      return
      end subroutine number2string_I


      subroutine tab2blank(string)
      use feminterface, only:
      use femtypes
      implicit none
      character(len=*) :: string
      intent (inout) :: string
!  convert all tabs (char(9)) to blanks
      integer(I4B) :: i
!
      do i=1,len_trim(string)
        if (ichar(string(i:i)) .eq. 9) then
          string(i:i)=' '
        endif
      enddo
      end subroutine

