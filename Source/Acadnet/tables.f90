      subroutine tables(lytabl,layanz,laytxt,layrb,txtlen)
      use feminterface, only: readky, laytab
      use femtypes
      implicit none
      integer (I4B) :: layanz, txtlen
      integer (I4B), pointer :: layrb(:,:)
      character (len=txtlen), pointer :: laytxt(:)
      logical :: lytabl
      intent (in) :: txtlen
      intent (out) :: layanz, lytabl ! ORG layrb, laytxt
      intent (inout) :: layrb, laytxt
! Pointer layrb works with oneAPI if pointer is inout, not out 
!    $Revision: 1.6 $
!    $Date: 2008/12/22 15:44:48 $
!    $Author: m_kasper $
!
!  Lesen der Section Tables
!
!  Input:
!    txtlen    defined length for strings txtgeb
!  Output:  
!    layanz    Anzahl der gelesenen Layer
!    laytxt    Text-Bezeichnung der Layer
!    layrb     type of boundary condition for branches on that layer (and for each nature)
!              0    = Dirichlet
!              1-99 = User Dirichlet
!              200  = general Neumann
!              300  = innner contour
!              400  = non-visble contour
!              1000 = comment layer
!              1001 = not recognized layer, treated as a comment layer
!    lytabl    =.true. wenn eine Layertabelle gelesen wurde
!
!  local variables
      integer (I4B) :: int, key, line, i
      real (DP) :: float
      character (len=txtlen) :: string
!
      lytabl=.false.

100   call readky(key,string,float,int,line)
      if ((key.eq.0) .and. (string .eq. 'ENDSEC')) return
      if ((key.eq.0) .and. (string .eq. 'TABLE')) then
        call readky(key,string,float,int,line)
        if ((key.eq.2) .and. (string .eq. 'LAYER')) then
!  Layer Tabelle lesen
          call laytab(layanz,laytxt,layrb,txtlen)
          lytabl=.true.
        end if
      end if
      goto 100
      return
      end subroutine tables
