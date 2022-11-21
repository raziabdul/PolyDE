      integer kdim
      parameter (kdim=450000)
      integer maxgbz,maxgkz,maxgzz
      parameter (maxgbz=30000,maxgkz=60000,maxgzz=30000)
!     kdim      Anzahl der zulaessigen (Gesamt-) Netzknoten
!     mazgbz    Maximale Anzahl der Gebiete
!     maxgkz    Maximale Anzahl der (Eingabe-) Knoten
!     maxgzz    Maximale Anzahl der (Eingabe-) Zweige
!  Achtung: Aufgrund der kompakten Abspeicherung der
!  Gebiets-Zweig-Information muss maxgbz um 1 groesser sein
!  als die Zahl der zu verarbeitenden Gebiete!