      subroutine gengeb(gbz,gzz,gkz,bzi,bzip,bzil,zki,xbk,ybk,area1)
      use feminterface, only: geng1, geng2
      use femtypes
      implicit none
      integer (I4B) :: gbz, gzz ,gkz, zki(:,:)
      integer (I4B), pointer :: bzi(:), bzip(:), bzil(:)
      real (DP) :: xbk(:), ybk(:)
      real (DP), pointer :: area1(:)
      intent (in) :: gzz, gkz, zki, xbk, ybk
!      intent (out) :: gbz, bzi, bzip, bzil, area1  ! org
      intent (inout) :: gbz, bzi, bzip, bzil, area1

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
!    $Revision: 1.6 $
!    $Date: 2008/12/22 13:09:36 $
!    $Author: m_kasper $
!
!-------------------------------------------------------------------------------
!
!  Generate the list of regions by connecting the branches
!
!  Input:
!    gzz        number of branches
!    gkz        number of key-points
!    maxgbz     maximum number of regions
!    zki        key-points (geometry nodes) which determine the branches
!               zki(1,i)  start key-point of the branch
!               zki(2,i)  end key-point of the branch
!               zki(3,i)  third key-point of the branch, with
!                         /  = 0  :  straight line 
!               zki(3,i)  -  > 0  :  node located on an arc              
!                         \  < 0  :  midpoint of an arc
!    xbk,ybk    coordinates of key-points
!  Output:
!    gbz        number of regions
!    bzi        is a continous list of branches, including the inner branches. 
!               The sequence is the same as in the input file (netin.acd)
!    bzip       is a pointer (index) to the first branch in the list bzi of a 
!               region. For multiply connected regions inner boundary branches 
!               precede the outer boundary.  
!    bzil       for each region is a pointer (index) to the first 
!               (possibly inner) region in the in the array bzip
!    area1      list with area of regions including the area of possibly inner
!               regions of multiply connected regions
!
!  local variables
      integer (I4B) :: i, j, count, negae, maxgbz
      integer (I4B), pointer :: plp(:), plist(:)
      real (DP), pointer :: area2(:)
      logical :: pos, neg
      parameter (maxgbz=100)
!
      do i=1,gzz
!  test for all branches whether start- end endpoint exist
        if ((zki(1,i).gt.gkz).or.(zki(2,i).gt.gkz)) then
          write (*,1001) i,zki(1,i),zki(2,i)
1001      format(' Error at branch ',i4,' node',i4,' or',i4,            &
     &      ' are not in the list of key-points')
          stop 10
        end if
        if (zki(1,i).eq.zki(2,i)) then
!  test whether start- end endpoint are identical
          write (*,1002) i,zki(1,i)
1002      format(' Error at branch',i4,' endpoints are identical nr.=',i4)
        end if
      end do
      write (*,8740)
8740  format(' determine regions')
      write (*,8743)
8743  format('  pass 1 : find simply connected regions')
!
      allocate (bzi(2*gzz), plist(2*gzz))
      allocate (area1(maxgbz), area2(maxgbz) )
      allocate (bzip(2*maxgbz+1), plp(2*maxgbz+1))
!
      gbz=0
      bzip(1)=1
      negae=0
!  plp      Point-List Pointer (compact storage)
!  plist    list of nodes of the actual region
      plp(1)=1
!  for all branches
      do i=1,gzz
        count=0
!  innere Zweige muessen einmal mehr (zweimal) durchlaufen werden
!  als auessere Zweige
        pos=.false.
        neg=.false.
        do j=1,bzip(gbz+1)-1
!  Fuer alle bereits generierten Gebiete:
          if (abs(bzi(j)).eq.i) then
!  Zweig schon einmal verwendet?
            count=count+1
!  Speichere die Richtung, in der er bereits verwendet wurde.
            if (bzi(j).eq.i) then
              pos=.true.
            else
              neg=.true.
            end if
          end if
        end do
        do j=1,plp(negae+1)-1
!  Fuer alle bereits generierten Gebiete:
          if (abs(plist(j)).eq.i) then
!  Zweig schon einmal verwendet?
            count=count+1
!  Speichere die Richtung, in der er bereits verwendet wurde.
            if (plist(j).eq.i) then
              pos=.true.
            else
              neg=.true.
            end if
          end if
        end do
        if (count.lt.2) then
!  Wenn der Zweig noch mindestens einmal durchlaufen werden muss,
!  dann rufe Unterprogramm zur Gebietsgenerierung mit Angabe
!  der noch zu durchlaufenden Richtung
          if ( .not. pos) then
            call geng1(i,+1,gbz,gzz,bzi,bzip,zki,xbk,ybk,area1,area2,   &
     &        negae,plp,plist)
          end if
          if ( .not. neg ) then
            call geng1(i,-1,gbz,gzz,bzi,bzip,zki,xbk,ybk,area1,area2,   &
     &        negae,plp,plist)
          end if
        end if
      end do
!
!  Setze jetzt die Startwerte fuer bzil
      allocate (bzil(gbz+1))
      bzil=(/ (i,i=1,gbz+1) /)
!
!  Wenn mehr als ein Gebiet mit negativer Flaeche gefunden wurde,
!  so liegt der Fall mehrfach zusammenhaengender Gebiete vor.
!  Rufe geng2, um diese Gebiete aufzuloesen.
!
      if (negae.gt.1) then
        write (*,8744)
8744    format('  pass 2 : find multiply connected regions')
        call geng2(gbz,bzi,bzip,bzil,                                   &
     &    zki,xbk,ybk,area1,area2,negae,plp,plist)
      end if
      deallocate (area2, plist, plp)
      return
      end subroutine gengeb
!
!
!
      subroutine geng1(startz,startd,gbz,gzz,bzi,bzip,zki,xbk,ybk,      &
     &  area1,area2,negae,plp,plist)
      use feminterface, only: inkrp, wink, reallocate
      use femtypes
      implicit none
      integer (I4B) :: startz, startd, gbz, gzz, zki(:,:), negae
      integer (I4B), pointer :: bzi(:), bzip(:), plp(:), plist(:)
      real (DP) :: xbk(:), ybk(:)
      real (DP), pointer :: area1(:), area2(:)
      intent (in) :: startz, startd, gzz, zki, xbk, ybk
!      intent (out) :: bzi, area1, area2, plist  ! org
      intent (inout) :: bzi, area1, area2, plist ! test
      intent (inout) :: gbz, bzip, negae, plp
!
!  generate a region for a specific starting branch
!
!  Input:
!    startz     branch for which we aim to find the region
!    startd     orientation of the starting branch
!    maxgbz     maximum number of regions
!    gzz        number of branches
!    zki        key-points (geometry nodes) which determine the branches
!               zki(1,i)  start key-point of the branch
!               zki(2,i)  end key-point of the branch
!               zki(3,i)  third key-point of the branch, with
!                         /  = 0  :  straight line 
!               zki(3,i)  -  > 0  :  node located on an arc              
!                         \  < 0  :  midpoint of an arc
!    xbk,ybk    coordinates of key-points
!  Output:
!    bzi        is a continous list of branches, including the inner branches. 
!               The sequence is the same as in the input file (netin.acd)
!    area1      list with area of regions including the area of possibly inner
!               regions of multiply connected regions
!    area2      area of regions with negative area
!    plist      list of nodes of the actual region
!  In-/ Output:
!    gbz        number of regions
!    bzip       is a pointer (index) to the first branch in the list bzi of a 
!               region. For multiply connected regions inner boundary branches 
!               precede the outer boundary.  
!    negae      number of regions with negative area
!    plp        Point-List Pointer (compact storage)
!
!  lcal variables
      real (DP) :: pi
      parameter (pi=3.14159265358979323846264338327950288419716_DP)
      integer (I4B) :: zaehl, neuzw, altzw, apunkt, lpunkt, vpunkt, i
      integer (I4B) :: maxgbz, dir, dir1
      real (DP) :: winkel, wink1, x, y, x1, y1, x2, y2, flaech, umfang
      real (DP) :: radius, xmitte, ymitte, flae2, w1, w2
      real (DP) :: xinnen, yinnen, rinnen
!
      if (gbz+1.lt.size(area1)) then
        gbz=gbz+1
      else
        maxgbz=2*size(area1)
        area1 => reallocate(area1,maxgbz)
        area2 => reallocate(area2,maxgbz)
        bzip  => reallocate(bzip,2*maxgbz+1)
        plp   => reallocate (plp,2*maxgbz+1)
!        write (*,1003)
!1003    format(' Error: too many regions')
      end if
!  Kompakte Abspeicherung:
!  Zaehl: Zeiger auf ersten Zweig des neuen Gebiets
      zaehl=bzip(gbz)
!  Setze ersten Zweig mit Richtung:
      bzi(zaehl)=startd*startz
      neuzw=startz
!  Merke den Anfangspunkt des ersten Zweigs
      if (startd.gt.0) then
        apunkt=zki(1,startz)
      else
        apunkt=zki(2,startz)
      end if
      dir=startd
!
3000  if (dir.gt.0) then
        vpunkt=zki(1,neuzw)
        lpunkt=zki(2,neuzw)
      else
        vpunkt=zki(2,neuzw)
        lpunkt=zki(1,neuzw)
      end if
!  Vpunkt: Vorletzter Punkt
!  Lpunkt: Letzter (aktueller) Punkt auf dem Rand
      altzw=neuzw
      neuzw=0
      winkel=-1._DP
      x=xbk(lpunkt)
      y=ybk(lpunkt)
!  Vermerke Koordinaten (x,y) des aktuellen Punkts
      if (zki(3,altzw).eq.0) then
        x1=xbk(vpunkt)
        y1=ybk(vpunkt)
      else
        if (zki(3,altzw).lt.0) then
          xmitte=xbk(-zki(3,altzw))
          ymitte=ybk(-zki(3,altzw))
        else
          call inkrp(xbk(vpunkt),ybk(vpunkt),x,y,xbk(zki(3,i)),         &
     &      ybk(zki(3,i)),                                              &
     &      xinnen,yinnen,flae2,rinnen,xmitte,ymitte,radius)
        end if
        y1=y+(x-xmitte)*dble(-dir)
        x1=x-(y-ymitte)*dble(-dir)
      end if
!  Vermerke Koordinaten (x1,y1), der vorangegangenen Punkts
      do i=1,gzz
!  Suche in der Zweigliste alle moeglichen Verzweigungen
        if ((i.ne.altzw).and.                                           &
     &    ((zki(1,i).eq.lpunkt).or.(zki(2,i).eq.lpunkt))) then
          if (zki(3,i).eq.0) then
!  Eine Gerade? Dann merke naechsten Punkt in (x2,y2)
            if (zki(1,i).eq.lpunkt) then
              x2=xbk(zki(2,i))
              y2=ybk(zki(2,i))
              dir1=1
            else
              x2=xbk(zki(1,i))
              y2=ybk(zki(1,i))
              dir1=-1
            end if
          else
            if (zki(3,i).lt.0) then
              xmitte=xbk(-zki(3,i))
              ymitte=ybk(-zki(3,i))
            else
              call inkrp(xbk(zki(1,i)),ybk(zki(1,i)),xbk(zki(2,i)),     &
     &          ybk(zki(2,i)),xbk(zki(3,i)),ybk(zki(3,i)),              &
     &          xinnen,yinnen,flae2,rinnen,xmitte,ymitte,radius)
            end if
!  Beim Kreisbogen werden die Koordinaten eines Punkts (x2,y2)
!  vermerkt, der auf der Tangente durch den aktuellen Punkt
!  liegt.
            if (zki(1,i).eq.lpunkt) then
              dir1=1
              y2=ybk(zki(1,i))+(xbk(zki(1,i))-xmitte)
              x2=xbk(zki(1,i))-(ybk(zki(1,i))-ymitte)
            else
              dir1=-1
              y2=ybk(zki(2,i))+(xbk(zki(2,i))-xmitte)*(-1._DP)
              x2=xbk(zki(2,i))-(ybk(zki(2,i))-ymitte)*(-1._DP)
            end if
          end if
!  ermittle den Winkel, unter dem die Verzweigung weitergeht
!  und verwende schliesslich den kleinsten Winkel, denn nur dann
!  ergibt sich ein Gebiet
          wink1=wink(x,y,x1,y1,x2,y2)
          if (wink1.gt.winkel) then
            winkel=wink1
            neuzw=i
            dir=dir1
          end if
        end if
      end do
      if (neuzw.eq.0) then
        write (*,2001)
2001    format(' Error: region cannot be closed')
        stop 10
      end if
!  Fuege den neuen ermittelten Zweig an die Liste fuer das Gebiet an.
      zaehl=zaehl+1
      if (zaehl.gt.2*gzz) then
!  Dieser Fall kann bei gueltigen erkannten Gebieten nicht
!  auftreten
        gbz=gbz-1
        return
      end if
      bzi(zaehl)=dir*neuzw
      if ((zki(1,neuzw).ne.apunkt).and.(zki(2,neuzw).ne.apunkt)) then
!  noch kein geschlossener Umlauf? Dann suche naechsten Zweig.
        goto 3000
      end if
      bzip(gbz+1)=zaehl+1
!  Das neue Gebiet ist fertig bestimmt.
!
      umfang=0._DP
      flaech=0._DP
!  Ermittle jetzt den Umfang und die Flaeche des neuen Gebiets
      do i=bzip(gbz),bzip(gbz+1)-1
        if (bzi(i).gt.0) then
          x1=xbk(zki(1,bzi(i)))
          y1=ybk(zki(1,bzi(i)))
          x2=xbk(zki(2,bzi(i)))
          y2=ybk(zki(2,bzi(i)))
        else
          x1=xbk(zki(2,-bzi(i)))
          y1=ybk(zki(2,-bzi(i)))
          x2=xbk(zki(1,-bzi(i)))
          y2=ybk(zki(1,-bzi(i)))
        end if
        if (zki(3,abs(bzi(i))).eq.0) then
          umfang=umfang+sqrt((x2-x1)**2+(y2-y1)**2)
!        flaech=flaech+0.25_DP*((y2-y1)*(x1+x2)-(x2-x1)*(y1+y2))
          flaech=flaech+0.5_DP*(x1*(y2-y1)-y1*(x2-x1))
        else
          if (zki(3,abs(bzi(i))).lt.0) then
            xmitte=xbk(-zki(3,abs(bzi(i))))
            ymitte=ybk(-zki(3,abs(bzi(i))))
            radius=sqrt((x1-xmitte)**2+(y1-ymitte)**2)
          else
            call inkrp(x1,y1,x2,y2,xbk(zki(3,abs(bzi(i)))),             &
     &        ybk(zki(3,abs(bzi(i)))),xinnen,yinnen,flae2,rinnen,       &
     &        xmitte,ymitte,radius)
          end if
          wink1=wink(xmitte,ymitte,x1,y1,x2,y2)
          if (bzi(i).lt.0) wink1=wink1-2._DP*pi
          umfang=umfang+radius*abs(wink1)
          w1=atan2(y1-ymitte,x1-xmitte)
          w2=atan2(y2-ymitte,x2-xmitte)
          flaech=flaech+                                                &
     &      (radius**2)/2._DP*wink1+                                    &
     &      radius/2._DP*(xmitte*(sin(w2)-sin(w1))-                     &
     &      ymitte*(cos(w2)-cos(w1)))
        end if
      end do
      if (flaech.lt.0) then
!  Wenn die Flaeche kleiner als Null ist, wurde auf dem Rand in der
!  falschen Richtung begonnen, oder es wurde eine innere Berandung
!  einer multiply connected region ermittelt. In diesem Fall
!  ist das Randpolygon zur weiteren Verwendung zu erhalten, so dass
!  es im Zwischenspeicher fuer den Pass 2 abgelegt wird.
!
        negae=negae+1
        do i=1,bzip(gbz+1)-bzip(gbz)
          plist(plp(negae)+i-1)=bzi(bzip(gbz)+i-1)
        end do
        plp(negae+1)=plp(negae)+bzip(gbz+1)-bzip(gbz)
        area2(negae)=flaech
!
!  Loesche die Daten wieder aus der echten Gebietsliste.
        gbz=gbz-1
      else
!  Ein neues gueltiges Gebiet wurde bestimmt.
        area1(gbz)=flaech
        write (*,4001) gbz,umfang,flaech
4001    format(' region',i5,'  perimeter:',g13.5,' area:',g13.5)
      end if
      return
      end subroutine geng1
!
!
!
      subroutine loesch(negae,merk,plp,plist,area2)
      use feminterface, only:
      use femtypes
      implicit none
      integer (I4B) :: negae, merk, plp(:), plist(:)
      real (DP) :: area2(:)
      intent (in) :: merk
      intent (inout) :: negae, plist, plp, area2
!
!  remove a region from lists plp, plist
!
!  Input:
!    merk       index of the region to be removed
!  In-/ Output:
!    negae      number of regions with negative area
!    plist      list of nodes of the actual region
!    plp        Point-List Pointer (compact storage)
!    area2      area of regions with negative area
!
!  local variables
      integer (I4B) :: i, merk1
!
      do i=1,plp(negae+1)-plp(merk+1)
        plist(plp(merk)+i-1)=plist(plp(merk+1)+i-1)
      end do
      do i=merk+1,negae
        merk1=plp(i+1)-plp(i)
        plp(i)=plp(i-1)+merk1
      end do
      do i=merk,negae-1
        area2(i)=area2(i+1)
      end do
      negae=negae-1
      return
      end subroutine loesch
!
!
!
      subroutine geng2(gbz,bzi,bzip,bzil,zki,xbk,ybk,area1,area2,     &
     &  negae,plp,plist)
      use feminterface, only: loesch, pinnen
      use femtypes
      implicit none
      integer (I4B) :: gbz
      integer (I4B) :: bzi(:), bzip(:), bzil(:), zki(:,:)
      integer (I4B) :: negae, plp(:), plist(:)
      real (DP) :: xbk(:), ybk(:), area1(:), area2(:)
      intent (in) :: gbz, zki, xbk, ybk, area1
      intent (out) :: bzi, bzip, bzil, area2, plist, plp
      intent (inout) :: negae
!
!  find multiply connected regions     
!
!  Input:
!    gbz        number of regions
!    zki        key-points (geometry nodes) which determine the branches
!               zki(1,i)  start key-point of the branch
!               zki(2,i)  end key-point of the branch
!               zki(3,i)  third key-point of the branch, with
!                         /  = 0  :  straight line 
!               zki(3,i)  -  > 0  :  node located on an arc              
!                         \  < 0  :  midpoint of an arc
!    xbk,ybk    coordinates of key-points
!    area1      list with area of regions including the area of possibly inner
!               regions of multiply connected regions
!  Output:
!    bzi        is a continous list of branches, including the inner branches. 
!               The sequence is the same as in the input file (netin.acd)
!    bzip       is a pointer (index) to the first branch in the list bzi of a 
!               region. For multiply connected regions inner boundary branches 
!               precede the outer boundary.  
!    bzil       for each region is a pointer (index) to the first 
!               (possibly inner) region in the in the array bzip
!    area2      area of regions with negative area
!    plist      list of nodes of the actual region
!    plp        Point-List Pointer (compact storage)
!  In-/ Output:
!    negae      number of regions with negative area
!
!  Suche zuerst die Zweigfolge, die die groesste negative
!  Flaeche umfasst. Es muss sich um den aeusseren Rand in
!  negativer Umlaufrichtung handeln. Er hat nichts mit
!  mehrfach zusammenhaengenden Gebieten zu tun und
!  muss nicht weiter betrachtet werden.
!
!  local variables 
      integer (I4B) :: i, j, merk, lrand, keypt, ptr, oben, nbr, start
      real (DP) :: wert
      real (DP) :: x, y
!
      merk=0
      wert=0._DP
      do i=1,negae
        if (abs(area2(i)).gt.wert) then
          merk=i
          wert=abs(area2(i))
        end if
      end do
!  Loesche den aeusseren Rand:
      call loesch(negae,merk,plp,plist,area2)
!
2000  continue
!
!  Teste alle existierenden Gebiete
      merk=0
      keypt=zki(1,abs(plist(plp(negae))))
      x=xbk(keypt)
      y=ybk(keypt)
      do i=1,gbz
!
!  Ein auesserer Rand muss eine groessere Flaeche umschliessen
!  als ein innerer Rand:
!
        if (area1(i).gt.abs(area2(negae))) then
!
!  Wenn bereits ein auesserer Rand gefunden ist,
!  koennen alle groesseren aeusseren Raender
!  sofort ausgeschlossen werden.
          if (merk.ne.0) then
            if (area1(i).gt.area1(merk)) cycle
          end if
!  Pruefe nur den aeusseren Rand eines jeden Gebiets
!  (den "letzten Rand")
!
          lrand=bzil(i+1)-1
          do j=bzip(lrand),bzip(lrand+1)-1
!  verhindere, dass ein Rand gegen sich selbst in umgedrehter
!  Richtung verglichen wird.
            if (abs(bzi(j)).eq.abs(plist(plp(negae)))) cycle
          end do
!
          if (pinnen(bzi,bzip(lrand),bzip(lrand+1)-1,                   &
     &      zki,xbk,ybk,x,y)) then
            merk=i
          end if
        end if
      end do
      if (merk.eq.0) then
        write (*,3500)
3500    format(                                                         &
     &    ' Error: can not resolve multiply connected region')
        stop 23
      else
!  Anzahl der einzufuegenden Zweige:
        nbr=plp(negae+1)-plp(negae)
!  An die Stelle Start soll eingefuegt werden
        ptr=bzil(merk+1)-1
        start=bzip(ptr)
!  Oben ist die Adresse des letzten gueltigen Zweiges
        oben=bzip(bzil(gbz+1))-1
!  Jetzt schaffe Platz fuer das Einfuegen der Liste
3600    continue
        bzi(oben+nbr)=bzi(oben)
        oben=oben-1
        if (oben.ge.start) goto 3600
!  Jetzt korrigiere die Pointer:
        oben=bzil(gbz+1)
3610    continue
        bzip(oben+1)=bzip(oben)+nbr
        oben=oben-1
        if (oben.ge.ptr) goto 3610
!
        do i=merk+1,gbz+1
          bzil(i)=bzil(i)+1
        end do
!
!  Nachdem alles vorbereitet ist, muss nun die Liste
!  des inneren Rands eingefuegt werden.
        do i=1,nbr
          bzi(start+i-1)=plist(plp(negae)+i-1)
        end do
!  Fertig
      end if
!
      negae=negae-1
      if (negae.gt.0) goto 2000
!
      return
      end subroutine geng2
