      subroutine readky(key,string,float,int,line)
      use feminterface, only:
      use femtypes
      implicit none
      real (DP) :: float
      integer (I4B) :: line, int, key
      character (len=*) :: string
      intent (out) :: key, string, float, int, line
!
!    $Revision: 1.5 $
!    $Date: 2008/12/22 15:23:06 $
!    $Author: m_kasper $
!
!  read and return key (Group Code) and the corresponding value 
!  from the following line 
!  first line contains the key, depending on its value, the second line
!  is of type integer, floating-point or a character string
!
!  Output: 
!                 key       Schluesselwert
!                 string    String
!                 float     Real
!                 int       Integer
!                 line      Aktuelle Zeilennummer
!  Kommentare (key=999) werden nicht zurueckgegeben
!
!  local variables
      integer (I4B) :: iread, pos, lineo,i
      integer (I4B), parameter :: ndat=10
      integer (I4B) :: ibuf(ndat), keybuf(ndat)
      real (DP) :: fbuf(ndat)
      character (len=1000) :: chars(ndat)
      logical :: first,ende
      save :: pos, iread, lineo, ende
      save :: ibuf, keybuf, fbuf, chars, first
      data first /.true./
!
      if (first) then
        first=.false.
        lineo=0
        pos=0
        iread=0
        ende=.false.
      end if
!
      if (pos.eq.iread) then
        if (ende) goto 998
!
!  Lesen von ndat Datensaetzen vom Eingabefile (bis zum Ende des Files)
!         iread   Anzahl der (tatsaechlich) gelesenen Datensaetze
!
        iread=0
        do i=1,ndat
!  read one data set (i.e. one line)
101       read (unit=38,fmt=*,end=200,err=998) keybuf(i)
!
          select case (keybuf(i))
          case ( 0:9 )
!  key  = 0 .. 9  => string
            read (unit=38,fmt=333,end=200,err=998) chars(i)
            if (chars(i) .eq. 'EOF' ) then
!  ignoriere alles auf dem File nach dem Text EOF
              ende=.true.
              iread=i
              goto 200
            end if
333         format(a)
!
          case ( 10:59 )
!  key  = 10 .. 59  => float
            read (unit=38,fmt=*,end=200,err=998) fbuf(i)
!
          case ( 60:79 )
!  key  = 60 .. 79  => 16 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 90:99 )
!  key  = 90 .. 99  => 32 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 100, 102)
!  key  = 100 102  => string  (255 character)
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 105)
!  key  = 105     => string  representing hex value
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 110:149 )
!  key  = 110 .. 149  => float double precision
            read (unit=38,fmt=*,end=200,err=998) fbuf(i)
!
          case ( 170:179 )
!  key  = 170 .. 175  => 16 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 210:239 )
!  key  = 210 .. 239  => float double precision
            read (unit=38,fmt=*,end=200,err=998) fbuf(i)
!
          case ( 270:289 )
!  key  = 270 .. 279  => 16 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 290:299 )
!  key  = 290 .. 299  => boolean flag vaue
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 300:309)
!  key  = 300 .. 309  => arbitrary string
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 310:319)
!  key  = 310 .. 319  => string  representing hex value
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 320:329)
!  key  = 320 .. 329  => string  representing hex handle value
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 330:369 )
!  key  = 330 .. 369  => string  representing hex object ID
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 370:389 )
!  key  = 370 .. 389  => 16 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 390:399)
!  key  = 390 .. 399  => string  representing hex handle value
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 400:409 )
!  key  = 400 .. 409  => 16 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 410:419)
!  key  = 410 .. 419  => string
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 420:429 )
!  key  = 420 .. 429  => 32 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 430:439)
!  key  = 430 .. 439  => string
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 440:449 )
!  key  = 440 .. 449  => 32 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 450:459 )
!  key  = 1010 .. 1059  => float long
            read (unit=38,fmt=*,end=200,err=998) fbuf(i)
!
          case ( 460:469 )
!  key  = 460 .. 469  => float double precision
            read (unit=38,fmt=*,end=200,err=998) fbuf(i)
!
          case ( 470:479)
!  key  = 470 .. 479  => string
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 1000:1009 )
!  key  = 1000 .. 1009  => string
            read (unit=38,fmt=333,end=200,err=998) chars(i)
!
          case ( 1010:1059 )
!  key  = 1010 .. 1059  => float double precision
            read (unit=38,fmt=*,end=200,err=998) fbuf(i)
!
          case ( 1060:1070 )
!  key  = 1060 .. 1079  => 16 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 1071 )
!  key  = 1071  => 32 bit integer
            read (unit=38,fmt=*,end=200,err=998) ibuf(i)
!
          case ( 999 )
!  key  = 999  => Kommentar
            read (unit=38,fmt=*,end=200,err=998)
            goto 101
!
          case default 
            print*,'*ERROR: unknown key',keybuf(i),                         &
     &        ' read in line :',lineo+2*i-1
            read (unit=38,fmt=333,end=200,err=998) chars(i)
          end select
!
          iread=i
        end do
!
        goto 300
!
200     ende = .true.
300     pos=0
      end if
!
      pos=pos+1
      lineo = lineo + 2
      line=lineo
      key=keybuf(pos)
      string=chars(pos)
      int=ibuf(pos)
      float=fbuf(pos)
!
      return
998   print*,'*ERROR: read ERROR on the inut file at line ',lineo+iread+1
      return
!
!
      entry unread
!  stellt die letzte durch Readky gelesene zeile zurueck, so dass
!  diese beim naechsten Aufruf von Readky nochmals ausgegeben wird
      pos=pos-1
      lineo=lineo-2
      if (pos.lt.0) print*,'*ERROR in unread'
      return
      end subroutine readky
!
!
!