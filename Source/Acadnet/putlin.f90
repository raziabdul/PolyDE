      subroutine putlin(x1,y1,x2,y2,r,phi1,phi2,gerad,accur,fakt1,      &
     &  layer,ok,xpoint,ypoint,zki,lzrb,zpz,length,anzzwg,anzknt)
      use feminterface, only: putnod, snzwg, reallocate
      use femtypes
      implicit none
      integer (I4B) :: layer, anzzwg, anzknt
      integer (I4B), pointer :: lzrb(:), zki(:,:), zpz(:)
      real (DP) :: x1, x2, y1, y2, r, phi1, phi2, accur, fakt1
      real (DP), pointer :: xpoint(:), ypoint(:), length(:)
      logical :: gerad, ok
      intent (in) :: x1, y1, x2, y2, r, phi1, phi2, gerad, layer, accur
      intent (in) :: fakt1
!      intent (out) :: length, ok ! org
      intent (inout) :: xpoint, ypoint, zki, lzrb, zpz, anzzwg, anzknt
      intent (inout) :: length, ok ! test

      !
!    $Revision: 1.6 $
!    $Date: 2008/12/22 15:21:39 $
!    $Author: m_kasper $
!
!  add a new branch legt einen neuen Zweig neu ab
!
!  Input:
!             x1,y1      Anfangspunkt   (bei Kreis Mittelpunkt)
!             x2,y2      Endpunkt
!             r          Radius (nur fuer Kreise)
!             phi1,phi2  anfangs und endwinkel (gradmass)
!             gerad      =.true. fuer  geraden
!                        =.false. fuer kreisboegen
!             layer      Layer der Linie
!             accur      Genauigkeit
!             fakt1      Punktdichte fuer Kreisboegen
!  In-/ Output: 
!             xpoint     Koordinaten-liste der bereits eingetragenen Punkte
!             ypoint           (wird aktualisiert)
!             zki        Zweig-Knoten Information
!             lzrb       Randbedingungen der Zweige
!             zpz        Anzagl der Knoten auf dem Zweig
!             anzzwg     Anzahl der Zweige
!             anzknt     Anzahl der Knoten
!  Output:
!             length     Zweiglaengen
!             ok        .true.
!                       .false. Wenn der Zweig nicht in die Liste
!                        aufgenommen wurde
!
!  local variables
      integer (I4B) :: n1, n2, n3, i, j
      real (DP) :: xa, ya, xe, ye, xp, yp, delphi, len
      real (DP) :: phi1b, phi2b
      real (DP), parameter ::  pid180=.017453292519943295769237_DP
      logical :: ls
!
      ok=.true.
      if (gerad) then
!
!  Teilprogramm fuer Geraden
!
        len=sqrt( (x1-x2)**2 + (y1-y2)**2 )
        if (len.lt.accur) then
          ok=.false.
          print*,'*WARNING length of segment is too short'
          print*,' from ( ',x1, ', ',y1,' )',' to ( ',x2, ', ',y2,' )'
          return
        end if
! reallocate arrays if necessary
        if (anzknt+2.gt.size(xpoint))                                   &
     &      xpoint=>reallocate(xpoint,2*size(xpoint)+2)
        if (anzknt+2.gt.size(ypoint))                                   &
     &      ypoint=>reallocate(ypoint,2*size(ypoint)+2)

        if (anzknt+2.gt.size(xpoint))                                   &
     &      xpoint=>reallocate(xpoint,2*size(xpoint))
        if (anzknt+2.gt.size(ypoint))                                   &
     &      ypoint=>reallocate(ypoint,2*size(ypoint))
        call putnod(x1,y1,accur,n1,xpoint,ypoint,anzknt)
        call putnod(x2,y2,accur,n2,xpoint,ypoint,anzknt)
        do i=1,anzzwg
          if ( zki(3,i) .eq. 0 ) then
            if ( (n1 .eq. zki(1,i) .and. n2 .eq.  zki(2,i) )  .or.      &
     &        (n2 .eq. zki(1,i) .and. n1 .eq.  zki(2,i) )  ) then
!             zweig existiert schon
              print*,'*WARNING: line segment already exists'
              print*,' from ( ',x1,',',y1,' )',' to ( ',x2,',',y2,' )'
              ok=.false.
              return
            end if
          end if
      end do
!  Zweig neu aufnehmen
!  reallocate arrays if neccessary
        if (anzzwg+1.gt.size(zki,2))                                    &
     &      zki=>reallocate(zki,size(zki,1),2*size(zki,2)+1)
        if (anzzwg+1.gt.size(lzrb))                                     &
     &    lzrb=>reallocate(lzrb,2*size(lzrb)+1)
        if (anzzwg+1.gt.size(length))                                   &
     &    length=>reallocate(length,2*size(length)+1)
        if (anzzwg+1.gt.size(zpz))                                      &
     &    zpz=>reallocate(zpz,2*size(zpz)+1)
        anzzwg=anzzwg+1
        zki(1,anzzwg)=n1
        zki(2,anzzwg)=n2
        zki(3,anzzwg)=0
        lzrb(anzzwg)=layer
        length(anzzwg)=len
        zpz(anzzwg) = 3
!       print*,'line: ',x1,y1,x2,y2,len
      else
!
!  Teilprogramm fuer Kreisboegen
!
        phi1b=phi1
        phi2b=phi2
        if (phi2b.lt.phi1b) phi2b=phi2b+360._DP
        phi1b=phi1b*pid180
        phi2b=phi2b*pid180
        delphi=phi2b-phi1b
        len=r * delphi
        if (len.lt.accur) then
          ok=.false.
          print*,'*WARNING: arc length is too short'
          print*,' angle form ',phi1,' to ',phi2,' radius = ',r
          return
        end if
!  reallocate arrays if necessary
        if (anzknt+3.gt.size(xpoint))                                   &
     &      xpoint=>reallocate(xpoint,2*size(xpoint)+3)
        if (anzknt+3.gt.size(ypoint))                                   &
     &      ypoint=>reallocate(ypoint,2*size(xpoint)+3)
!  Mittelpunkt
        call putnod(x1,y1,accur,n3,xpoint,ypoint,anzknt)
        n3=-n3
!  Anfang
        xa= x1 + r * cos(phi1b)
        ya= y1 + r * sin(phi1b)
        call putnod(xa,ya,accur,n1,xpoint,ypoint,anzknt)
!  Ende
        xe= x1 + r * cos(phi2b)
        ye= y1 + r * sin(phi2b)
        call putnod(xe,ye,accur,n2,xpoint,ypoint,anzknt)
!
        do i=1,anzzwg
          if ( zki(3,i) .ne. 0 ) then
            if ( (n1 .eq. zki(1,i)) .and. (n2 .eq. zki(2,i))            &
     &        .and. (n3 .eq. zki(3,i)) ) then
!  Zweig existiert schon
              print*,'*WARNING: arc segment already exists'
              print*,' angle from ',phi1,' to ',phi2,' radius = ',r
              ok=.false.
              return
            end if
          end if
        end do
!  Zweig neu aufnehmen
!  reallocate arrays if neccessary
        if (anzzwg+1.gt.size(zki,2))                                    &
     &      zki=>reallocate(zki,size(zki,1),2*size(zki,2)+1)
        if (anzzwg+1.gt.size(lzrb))                                     &
     &    lzrb=>reallocate(lzrb,2*size(lzrb)+1)
        if (anzzwg+1.gt.size(length))                                   &
     &    length=>reallocate(length,2*size(length)+1)
        if (anzzwg+1.gt.size(zpz))                                      &
     &    zpz=>reallocate(zpz,2*size(zpz)+1)
        anzzwg=anzzwg+1
        zki(1,anzzwg)=n1
        zki(2,anzzwg)=n2
        zki(3,anzzwg)=n3
        lzrb(anzzwg)=layer
        length(anzzwg)=len
!  Anzahl der Knoten auf einem Kreisbogen ist abhaengig vom
!  Winkel dieses Kreisbogens
        zpz(anzzwg) = int(delphi*fakt1*0.638_DP)
!  Garantiere mindestens 3 Knoten
        zpz(anzzwg) = max(3,zpz(anzzwg))
!      print*,'arc: ',x1,y1,r,len
      end if
!
!  testen auf Schnittpunkte
!
      do j=1,anzzwg-1
        call snzwg (xpoint,ypoint,zki,anzzwg,j,xp,yp,ls)
        if (ls) then
          write(*,666) anzzwg,j,xp,yp
666       format (' *ERROR: branches ',i6,' and ',i6,' intersect '      &
     &      ,'at the point x=',g13.6,'    y=',g13.6)
          ok=.false.
        end if
      end do
      return
      end subroutine putlin
