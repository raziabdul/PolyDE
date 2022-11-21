      subroutine writemshstats(ok)
      use feminterface, only: getsetting
      use feminterface3d, only: tetangles
      use femtypes
      use globalvariables3D, only : nod, numn, nums, numv, bcs, sn, dom, vn, nnat
      implicit none
      logical ok
      intent (out) :: ok
!
!------------------------------------------------------------------------------
!    $Revision: 1.2 $
!    $Date: 2014/02/11 16:24:37 $
!    $Author: juryanatzki $
!------------------------------------------------------------------------------
!
!  Write mesh data as neutral NETGEN mesh file. The file containing the mesh 
!  will be named basis.mesh. The information contained in this format is:
!
!  Neutral format:
!   1. Number of nodes (numn)
!   2. List of coordintes (nod)
!   3. Number of volume elements (numv)
!   4. domain (1 entry) and nodes (4 entries) for volume elements (dom, vn)
!   5. Number of surface elements (nums)
!   6. bc at surface (1 entry) and nodes (3 entries) for surface elements (bcs, sn)
!
!------------------------------------------------------------------------------
!  Input:
!    Regions
!      numdom      number of domains / regions
!      nummat      number for material in domain / region
!    Elements
!      numv        number of volume elements
!      nums        number of surface elements
!      sn          surface element -> nodes
!      vn          volume element -> nodes
!    Nodes
!      nod         world coordinates of nodes
!      numn        number of nodes
!      nsi         node - surface information
!                  NOT USED YET
!
!  Output:
!      ok          =.false. if an error occurs
!
!  local variables
      integer (I4B) i, unitid, ios
      real (DP) vert(3,4), tetangles_temp(4)
      character (len=200) :: path
      character (len=20) :: vnumber, filename, form
!
!  get project path from user's environment variables
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
!################# +      WRITE OUT MESH      + #################
!################# +          START           + #################
!
! open mesh.neutral file in created path for formatted reading (ASCII)
!
    form = ''
    WRITE(form,"(A5,I1,A1)") "(I", CEILING(LOG10(real(numv))), ")"
!
! Create new Directory for each Mesh

    vnumber = ''
    WRITE (vnumber,form) numv
    
    
      open (unitid,file=path(1:len_trim(path))//vnumber(1:len_trim(vnumber))//"basis_end.mesh",            &
     &   form='formatted',position='REWIND',action='WRITE',iostat=ios)
      if (ios .ne. 0) then
        print*, '**** Error on opening file: ',path(1:len_trim(path))//vnumber(1:len_trim(vnumber))//"basis_end.mesh"
        print*, '**** IO Error No: ',ios
        ok=.false.
        return
      end if
      rewind unitid
!
!  write number of nodes
      write (unitid,*) numn
      print "(a,i6,a)","  writing ", numn," nodes"
!
!  write coordinates
      do i=1,numn
        write (unitid,*) nod(:,i)
      end do
!
!  write number of volume elements
      write (unitid,*) numv
      print "(a,i6,a)","  writing ", numv," volume elements"
!
!  write domain and volume element - node map
      do i=1,numv
        write (unitid,*) dom(i),vn(:,i)
      end do
!
!  write number of surface elements
      write (unitid,*) nums
      print "(a,i6,a)","  writing ", nums," surface elements"
!  read domain and volume element - node map ------------------------------------ changed
!      do i=1,nums
!        read (unitid,*) bcs(i),sn(:,i)
!      end do
    !  do i=1,nums
    !   read(unitid,*) tempbcs, sn(:,i)
    !    do inat=1,nnat
    !      bcs(i,inat)=int(tempbcs/(1000**(nnat-inat)))
    !      tempbcs=tempbcs-(bcs(i,inat)*(1000**(nnat-inat)))
    !    end do
    !  end do
!--------------------------------------------------------------------------------change ends
!
!  write domain and volume element - node map
! QUICK-FIX for bcs(i,nnat)   TO BE corrected in future if writeng is in use
      do i=1,nums
        write (unitid,*) bcs(i,nnat),sn(:,i)
      end do
!
      close (unitid)
!
!################# +          END OF          + #################
!################# +      WRITE OUT MESH      + #################
!
!
!################# + MESH-QUALITY-INFORMATION + #################
!################# +          START           + #################
!
! Write into file "plotdata" - List of stored Plots:
!

! Header line Format for parsing in PolyDE: H#DataType;DataFileName;PlotType;PlotName
! Possible DataTypes: 
! Possible PlotTypes: XYPLOT, HISTOGRAM, FIELDPLOT, PLOTARRAY
! Examples:
!                             Preferred  |   Array      | Data Descript.|   Name of the Plot
!                                Plottype|      Shape   |    OPTIONAL   |
!    WRITE(unitid,1313) "HEADER;HISTOGRAM#",numv,";",4,"#Elements;Angles#Solid Angles of all Mesh Elements"
!                                        |              |               |
! Now writing of the actual dataFiles follows:

!form = ''
!    WRITE(form,"(A5,I1,A8)") "(A7,I", CEILING(LOG10(real(numv))), "A1,I1,A)"
!    
!1001    FORMAT (i)
!1002    FORMAT (f)
!1003    FORMAT (a)
!1004    FORMAT (f,f,f,f)
!1013    FORMAT (i,a)
!1313    FORMAT (A,I,A1,I1,A)
!
!! Write MeshStats
!    OPEN (unitid,file=path(1:LEN_TRIM(path))//"plotdata.dat", &
!    &     form='formatted',position='APPEND',action='WRITE',iostat=ios)
!    if (ios .ne. 0) then
!        print*, '**** Error on opening file: ',path(1:LEN_TRIM(path))//"plotdata.dat"
!        print*, '**** IO Error No: ',ios
!        ok=.false.
!        return
!    end if
!    
!    PRINT*,'Generating Plot Data for meshangles'
!    WRITE(unitid,1313) "HEADER;HISTOGRAM#",numv,";",4,"#Elements;Angles#Solid Angles of all Mesh Elements"
!    
!      DO i = 1, numv
!        vert(1:3,1:4) = nod( 1:3,vn(1:4,i) )
!        tetangles_temp = tetangles(vert)
!        WRITE(unitid,1004) tetangles_temp(1), tetangles_temp(2), tetangles_temp(3), tetangles_temp(4)
!      END DO
!    CLOSE(unitid)
    
!
!################# +          END OF          + #################
!################# + MESH-QUALITY-INFORMATION + #################






      ok=.true.
      return
      end subroutine writemshstats
