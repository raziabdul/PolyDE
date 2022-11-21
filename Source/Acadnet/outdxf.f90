      subroutine outdxf(xpoint,ypoint,zki,lzrb,laytxt,anzzwg,           &
     &  xgeb,ygeb,txtgeb,txthoc,anzgeb,lplayt)
      use feminterface, only: getsetting
      use femtypes
      implicit none
      real (DP) :: xpoint(:), ypoint(:), xgeb(:), ygeb(:), txthoc
      integer (I4B) :: lzrb(:), lplayt(:), zki(:,:), anzzwg, anzgeb
      integer (I4B) :: unitid
      character (len=*) :: txtgeb(:), laytxt(:)
      intent (in) :: xpoint, ypoint, zki, lzrb, laytxt, anzzwg
      intent (in) :: xgeb, ygeb, txtgeb, txthoc, anzgeb, lplayt
!
!    $Revision: 1.5 $
!    $Date: 2008/12/22 15:14:14 $
!    $Author: m_kasper $
!
!  output of the structure as a AutoCad drawing in DXF format to: acadnet.dxf 
!
!  Input:
!             xpoint    Liste der Knotenkoordinaten
!             ypoint
!             zki       Zweig-Knoten Information
!             lzrb      Randbedingungen der Zweige
!             anzzwg    Anzahl der Zweige
!             xgeb      Einfuegekoordinaten fuer Textstrings zur Kennzeichnung
!             ygeb      der Materialien der Gebiete
!             txtgeb    Texte der Gebiete
!             txthoc    Hoehe des Textes
!             anzgeb    Anzahl der Gebiete (Texte)
!
!             lplayt    Pointer auf layer fuer die Texte
!             layanz    Anzahl der gelesenen Layer
!             laytxt    Text-bezeichnung der Layer
!
!  local variables
      integer (I4B) :: i,lang
      real (DP) :: radius,rad1,rad2
      character (len=200) :: path
      real (DP), parameter ::  pid180=.017453292519943295769237_DP
!
!  Only the section ENTITIES is written
!
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
      open (unitid,file=path(1:len_trim(path))//'acadnet.dxf')
!
      write(unitid,2001) 0
!  Beginn einer SECTION
2001  format(i1)
2002  format(i2)
2003  format(i3)
2005  format(g16.10)
      write(unitid,1001)
1001  format('SECTION')
      write(unitid,2001) 2
!
!  Entities Section
!
      write(unitid,1002)
1002  format('ENTITIES')
!  Linien bzw. Kreisboegen ausgeben
      do i=1,anzzwg
        if (zki(3,i).eq.0) then
          write(unitid,2001) 0
!
!  Gerade
!
          write(unitid,1003)
1003      format ('LINE')
          write(unitid,2001) 8
! Layer
          lang=len_trim(laytxt(lzrb(i)))
          write(unitid,*) laytxt(lzrb(i))(1:lang)
          write(unitid,2002) 10
!  x-Koordinate des Anfangs - Punktes
          write(unitid,2005) xpoint(zki(1,i))
          write(unitid,2002) 20
!  y-Koordinate des Anfangs - Punktes
          write(unitid,2005) ypoint(zki(1,i))
          write(unitid,2002) 30
!  z-Koordinate des Anfangs -punktes
          write(unitid,2005) 0._DP
          write(unitid,2002) 11
!  x-Koordinate des End - Punktes
          write(unitid,2005) xpoint(zki(2,i))
          write(unitid,2002) 21
!  y-Koordinate des End - Punktes
          write(unitid,2005) ypoint(zki(2,i))
          write(unitid,2002) 31
!  z-Koordinate des End - Punktes
          write(unitid,2005) 0._DP
        else
          write(unitid,2001) 0
!
!  Kreisbogen
!
          write(unitid,1004)
1004      format ('ARC')
          write(unitid,2001) 8
!  Layer
          lang=len_trim(laytxt(lzrb(i)))
          write(unitid,*) laytxt(lzrb(i))(1:lang)
          write(unitid,2002) 10
!  x-Koordinate des Mittelpunkt
          write(unitid,2005) xpoint(abs(zki(3,i)))
          write(unitid,2002) 20
!  y-Koordinate des Mittelpunkt
          write(unitid,2005) ypoint(abs(zki(3,i)))
!      write(unitid,2002) 30
!      write(unitid,2005) 0.
!  Radius
          rad1 = sqrt( (xpoint(zki(1,i))-xpoint(abs(zki(3,i))))**2 +    &
     &      (ypoint(zki(1,i))-ypoint(abs(zki(3,i))))**2 )
          rad2 = sqrt( (xpoint(zki(2,i))-xpoint(abs(zki(3,i))))**2 +    &
     &      (ypoint(zki(2,i))-ypoint(abs(zki(3,i))))**2 )
          radius = ( rad1 + rad2 ) / 2._DP
          write(unitid,2002) 40
!  Radius
          write(unitid,2005) radius
!  Anfangswinkel
          write(unitid,2002) 50
!  Anfangswinkel in Grad
          write(unitid,2005)                                            &
     &      atan2(ypoint(zki(1,i))-ypoint(abs(zki(3,i))),               &
     &      xpoint(zki(1,i))-xpoint(abs(zki(3,i))) )/pid180
!  Endwinkel
          write(unitid,2002) 51
!  Endwinkel in Grad
          write(unitid,2005)                                            &
     &      atan2(ypoint(zki(2,i))-ypoint(abs(zki(3,i))) ,              &
     &      xpoint(zki(2,i))-xpoint(abs(zki(3,i))) )/pid180
        end if
      end do
!
      do i=1,anzgeb
!
!  Text
!
        write(unitid,2001) 0
        write(unitid,1006)
1006    format ('TEXT')
        write(unitid,2001) 8
!  Layer
        lang=len_trim(laytxt(lplayt(i)))
        write(unitid,*) laytxt(lplayt(i))(1:lang)
        write(unitid,2002) 10
!  x-Koordinate
        write(unitid,2005) xgeb(i)
        write(unitid,2002) 20
!  y-Koordinate
        write(unitid,2005) ygeb(i)
!     write(unitid,2002) 30
!  z-Koordinate
!     write(unitid,2005) 0.
        write(unitid,2002) 40
!  Texthoehe
        write(unitid,2005) txthoc
        write(unitid,2001) 1
!  Textwert (string)
        lang=len_trim(txtgeb(i))
        write(unitid,*) txtgeb(i)(1:lang)
      end do
!
      write (unitid,2001) 0
!  Ende der Section
      write (unitid,1008)
1008  format ('ENDSEC')
      write (unitid,2001) 0
!  Ende des Dxf-files
      write (unitid,1009)
1009  format ('EOF')
      close (unitid)
      return
      end subroutine outdxf
