      subroutine matgeb(gbz,anzgeb,bzi,bzip,bzil,zki,area1,xpoint,      &
     &                  ypoint,xgeb,ygeb,txtgeb,matname,donetx)
      use feminterface, only: pinnen
      use femtypes
      implicit none
      integer (I4B) :: gbz, anzgeb, bzi(:), bzip(:), bzil(:)
      integer (I4B) :: zki(:,:), donetx(:)
      real (DP) :: xpoint(:), ypoint(:), xgeb(:), ygeb(:), area1(:)
      character (len=*) :: txtgeb(:), matname(:)
      intent (in) :: gbz, anzgeb, bzi, bzip, bzil, zki, area1
      intent (in) ::  xpoint, ypoint
      intent (out) :: matname
      intent (inout) :: txtgeb
!
!    $Revision: 1.5 $
!    $Date: 2008/12/22 15:52:48 $
!    $Author: m_kasper $
!
!  assign text (materials) to the regions
!
!  Input:
!    gbz       number of regions
!    anzgeb    number of texts read from the DXF-file (number of regions)
!    bzi       is a continous list of branches, including the inner branches. 
!              The sequence is the same as in the input file (netin.acd)
!    bzip      is a pointer (index) to the first branch in the list bzi of a 
!              region. For multiply connected regions inner boundary branches 
!              precede the outer boundary.  
!    bzil      for each region is a pointer (index) to the first 
!              (possibly inner) region in the in the array bzip
!    zki       key-points (geometry nodes) which determine the branches
!               zki(1,i)  start key-point of the branch
!               zki(2,i)  end key-point of the branch
!               zki(3,i)  third key-point of the branch, with
!                         /  = 0  :  straight line 
!               zki(3,i)  -  > 0  :  node located on an arc              
!                         \  < 0  :  midpoint of an arc
!    area1     list with area of regions including the area of possibly inner
!              regions of multiply connected regions
!    xpoint    x-coordinates of key-points
!    ypoint    x-coordinates of key-points
!    xgeb      alignment point of text
!    ygeb
!    txtgeb    the string of TEXT and MTEXT entities
!    donetx    a flag which is 1 if the text has been recognized as a boundary condition, otherwise 0
!
!  Output:
!    matname   materialnames of the regions
!
!  local variables
      integer (I4B) :: i, j, merk, lang, lang2, start, end, found(gbz)
      real (DP) :: x, y, minfl
!
      matname='nothing'
      found=0
!
!  search for each text the smallest region inside which the text is located
!
      do i=1,anzgeb
        if (donetx(i) .eq. 1) cycle
        txtgeb(i)=adjustl(txtgeb(i))
        lang=len_trim(txtgeb(i))
        x=xgeb(i)
        y=ygeb(i)
        merk=0
        minfl=huge(1._DP)
        do j=1,gbz
!  use the last region in the list bzi - this region includs the contained
!  inner regions of multiply connected regions
!
!     first and last branch number in the list bzi
          start=bzip( bzil(j+1) - 1 )
          end  =bzip( bzil(j+1) ) - 1
          if (pinnen(bzi,start,end,zki,xpoint,ypoint,x,y)) then
            if ( area1(j) .lt. minfl ) then
              merk =j
              minfl = area1(j)
            end if
          end if
        end do
!
        if (merk .eq. 0) then
          print*, '*WARNING: could not find a region for the text: "',  &
     &      txtgeb(i)(1:lang),'"'
          else
!
!  for the region: merk the text: txtgeb(i) was found
!
          if ( found(merk) .ne. 0 )  then
!  if a text (material) has been assigned previvously
            print*,'*WARNING: for the same region multiple strings had been found'
            print*,txtgeb(i)(1:lang)
            lang2=len_trim(txtgeb(found(merk)))
            print*,txtgeb(found(merk))(1:lang2)
          end if
          found(merk)=i
!
!  output the text as the material name
!
          matname(merk)=txtgeb(i)(1:lang)
        end if
      end do
      return
      end subroutine matgeb
