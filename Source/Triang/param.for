      integer kdim
      parameter (kdim=450000)
      integer maxgbz,maxgkz,maxgzz
      parameter (maxgbz=3000,maxgkz=6000,maxgzz=6000)
!     kdim      Anzahl der zulaessigen (Gesamt-) Netzknoten
!     mazgbz    Maximale Anzahl der Gebiete
!     maxgkz    Maximale Anzahl der (Eingabe-) Knoten
!     maxgzz    Maximale Anzahl der (Eingabe-) Zweige
!  Achtung: Aufgrund der kompakten Abspeicherung der
!  Gebiets-Zweig-Information muss maxgbz um 1 groesser sein
!  als die Zahl der zu verarbeitenden Gebiete!