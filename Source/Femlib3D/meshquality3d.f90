      subroutine meshquality3d(vn,numv,nod)
      use femtypes
      use feminterface3d, only: tetmeanratio, diangles
      implicit none
      integer (I4B) :: vn(:,:), numv
      real (DP) :: nod(:,:)
      intent(in) :: vn, numv, nod
!
!------------------------------------------------------------------------------
!    $Revision: 1.8 $
!    $Date: 2014/08/22 10:43:07 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  meshquality3d - calculates the quality of all elemenet and display some statistics
!
!  input:
!  vn          element nodes
!  numv        number of elements
!  nod         node coordinates
!
!  output:     (ascii graphic)
!
      integer (I4B), parameter :: criterionselector=2
!  criterionselector = 1        for GEOMPACK's  eta  quality criterion
!                    = 2        for minimum dihedral angle (relative to regular tetrahedron)
      integer (I4B), parameter :: sections=25, height=20, leftmarg=8
      integer (I4B) :: i, j, occurence(0:sections)
      real (DP) :: quali(numv), min_q, max_q, averagequality, stddev, minquality, maxquality
      character(len=78) :: str
      character(2), parameter :: full_block = char(int(Z'E2'))//char(int(Z'96'))//char(int(Z'88'))
!                      horizontal     vertical       cross          left cross     right cross    fill           half-fill
!  extended Ascii (Codepage 437)
      character (1) :: ch=achar(196), cv=achar(179), cp=achar(197), cl=achar(195), cr=achar(180), cf=achar(219), cx=achar(220)
!  standard Ascii version
!      character (1) :: ch='-',        cv='|',        cp='+',        cl='+',        cr='+',        cf='X'      , cx='x'
!
!  equivalent Unicode code point, but not yet supported in Fortran
!      character (1) :: ch=achar(2500), cv=achar(2502), cp=achar(253C), cl=achar(251C), cr=achar(2524), cf=achar(2588), cx=achar(2584)                  
!
! UTF-8 encodings:
! codepage 437  unicode   unicode name                    UTF-8 hex
!   196         2500    'BOX DRAWINGS LIGHT HORIZONTAL'   0xE2 0x94 0x80
!   179
!   197
!   195
!   180
!   219         2588    'FULL BLOCK'                      0xE2 0x96 0x88
!   220
!   
      select case (criterionselector)

      case (1)

!$omp parallel do default(none) &
!$omp shared( nod,vn,quali,numv)     &
!$omp private( i)
        do i=1,numv
!  use GEOMPACK's  eta  quality criterion
          quali(i) = tetmeanratio(nod(1:3,vn(1:4,i)))
        end do
!$omp end parallel do

      case (2)

!$omp parallel do default(none) &
!$omp shared( nod,vn,quali,numv)     &
!$omp private( i)
        do i=1,numv
!  use minimum dihedral angle (relative to regular tetrahedron), 70.53° = 1.231 rad
          quali(i) = minval(diangles(nod(1:3,vn(1:4,i))))/1.23095941734077468213492917825_DP
        end do
!$omp end parallel do

      end select

      averagequality = sum(quali)/numv
      stddev = sqrt( dot_product(quali,quali)/numv - averagequality**2)
      minquality = minval(quali)
      maxquality = maxval(quali)
      
      min_q =  0._DP
      max_q =  1._DP
      
      occurence = 0
      do i = 1, numv
        j = nint((quali(i) - min_q)/(max_q - min_q)*sections)
        occurence(j) = occurence(j) + 1
      end do
      
      str(1:leftmarg)='                      '
      str(leftmarg-1:)='0.0'//repeat(ch,7)//'0.2'//repeat(ch,7)         &
     &               //'0.4'//repeat(ch,7)//'0.6'//repeat(ch,7)         &
     &               //'0.8'//repeat(ch,7)//'1.0'
      str(leftmarg+55:leftmarg+55+12)='Mesh Quality'
      print*,str

      do i= height, 1, -1
        select case (i)
        case (15)
          str(leftmarg- 5:leftmarg- 2)='15%'
          str(leftmarg+ 1:leftmarg+ 9)=repeat(ch,9)
          str(leftmarg+11:leftmarg+19)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+21:leftmarg+29)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+31:leftmarg+39)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+41:leftmarg+49)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg   :leftmarg   )=cl
          str(leftmarg+10:leftmarg+10)=cp
          str(leftmarg+20:leftmarg+20)=cp
          str(leftmarg+30:leftmarg+30)=cp
          str(leftmarg+40:leftmarg+40)=cp
          str(leftmarg+50:leftmarg+50)=cr
        case (10)
          str(leftmarg- 5:leftmarg- 2)='10%'
          str(leftmarg+ 1:leftmarg+ 9)=repeat(ch,9)
          str(leftmarg+11:leftmarg+19)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+21:leftmarg+29)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+31:leftmarg+39)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+41:leftmarg+49)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg   :leftmarg   )=cl
          str(leftmarg+10:leftmarg+10)=cp
          str(leftmarg+20:leftmarg+20)=cp
          str(leftmarg+30:leftmarg+30)=cp
          str(leftmarg+40:leftmarg+40)=cp
          str(leftmarg+50:leftmarg+50)=cr
        case (05)
          str(leftmarg- 5:leftmarg- 2)=' 5%'
          str(leftmarg+ 1:leftmarg+ 9)=repeat(ch,9)
          str(leftmarg+11:leftmarg+19)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+21:leftmarg+29)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+31:leftmarg+39)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg+41:leftmarg+49)=str(leftmarg+ 1:leftmarg+ 9)
          str(leftmarg   :leftmarg   )=cl
          str(leftmarg+10:leftmarg+10)=cp
          str(leftmarg+20:leftmarg+20)=cp
          str(leftmarg+30:leftmarg+30)=cp
          str(leftmarg+40:leftmarg+40)=cp
          str(leftmarg+50:leftmarg+50)=cr
        case default
          str(leftmarg:)='                                                  '
          str(leftmarg- 6:leftmarg- 3)='   '
          str(leftmarg- 1:leftmarg+ 2)=' '//cv//' '
          str(leftmarg+10:leftmarg+10)=cv
          str(leftmarg+20:leftmarg+20)=cv
          str(leftmarg+30:leftmarg+30)=cv
          str(leftmarg+40:leftmarg+40)=cv
          str(leftmarg+50:leftmarg+50)=cv
        end select

!  display some statistic at right margin
        select case (i)
        case (height)
          select case (criterionselector)
          case (1)
            str(leftmarg+55:leftmarg+55+12)='(eta quality)'
          case (2)
            str(leftmarg+55:leftmarg+55+12)='dihed. angles'
          end select
        case (height- 4)
          str(leftmarg+55:leftmarg+55+12)='aver. quality'
        case (height- 5)
          write(str(leftmarg+55:leftmarg+55+12),'(F9.7)') averagequality
        case (height- 8)
          str(leftmarg+55:leftmarg+55+12)='st. deviation'
        case (height- 9)
          write(str(leftmarg+55:leftmarg+55+12),'(F9.7)') stddev
        case (height-12)
          str(leftmarg+55:leftmarg+55+12)='min quality'
        case (height-13)
          write(str(leftmarg+55:leftmarg+55+12),'(F9.7)') minquality
        case (height-16)
          str(leftmarg+55:leftmarg+55+12)='max quality'
        case (height-17)
          write(str(leftmarg+55:leftmarg+55+12),'(F9.7)') maxquality
        case default
          str(leftmarg+55:leftmarg+55+12)='             '
        end select

        do j=0, sections
!  fill: lower bound is 75% of level, i.e. i-0.25
          if ( real(occurence(j))/real(numv)*100. .ge. (real(i)-0.25)) then 
            str(2*j+leftmarg : 2*j+leftmarg)=cf
!            str(2*j+leftmarg : 2*j+leftmarg)=full_block

            !  half-fill: lower bound is 25% of level, i.e. i-0.75
          else if ( real(occurence(j))/real(numv)*100. .ge. (real(i)-0.75)) then 
            str(2*j+leftmarg : 2*j+leftmarg)=cx
          end if
        end do
        print*,str
      end do

      str(1:leftmarg)='                      '
      str(leftmarg-1:)='0.0'//repeat(ch,7)//'0.2'//repeat(ch,7)         &
     &               //'0.4'//repeat(ch,7)//'0.6'//repeat(ch,7)         &
     &               //'0.8'//repeat(ch,7)//'1.0'
      print*,str
!      print*, occurence, sum(occurence)
      return
      end subroutine meshquality3d




      pure function diangles(vert)
      use femtypes
      use feminterface, only: cross_product
      use globalvariables3D, only: pi

      implicit none
      real   (DP) :: diangles(6), vert(3,4)
      intent (in) :: vert
!      intent (out) :: diangles
!
!  compute the dihedral angles of a tetrahedron
!
!  local variables
      real (DP) edg4(3), edg5(3), edg6(3)
      real (DP) fn1(3), fn2(3), fn3(3), fn4(3)

!  get edge vectors

!      edg1(1:3) = vert(1:3,2) - vert(1:3,1)
!      edg2(1:3) = vert(1:3,3) - vert(1:3,2)
!      edg3(1:3) = vert(1:3,3) - vert(1:3,1)
      edg4(1:3) = vert(1:3,4) - vert(1:3,1)
      edg5(1:3) = vert(1:3,4) - vert(1:3,2)
      edg6(1:3) = vert(1:3,4) - vert(1:3,3)

!  outward (or inward) faces normals

      fn1 = cross_product(edg5 , edg6)
      fn2 = cross_product(edg6 , edg4)
      fn3 = cross_product(edg4 , edg5)
      fn4 = -fn1 - fn2 - fn3

!  normalize
      fn1 = fn1 / sqrt(dot_product( fn1 , fn1 ))
      fn2 = fn2 / sqrt(dot_product( fn2 , fn2 ))
      fn3 = fn3 / sqrt(dot_product( fn3 , fn3 ))
      fn4 = fn4 / sqrt(dot_product( fn4 , fn4 ))
      
      diangles(1) = acos(dot_product( fn3 , fn4 ))
      diangles(2) = acos(dot_product( fn1 , fn4 ))
      diangles(3) = acos(dot_product( fn2 , fn4 ))
      diangles(4) = acos(dot_product( fn2 , fn3 ))
      diangles(5) = acos(dot_product( fn1 , fn3 ))
      diangles(6) = acos(dot_product( fn1 , fn2 ))

      diangles(1:6) = pi - diangles(1:6)

!  trihedral angles
!      triangles(1) = diangles(1) + diangles(3) + diangles(4)
!      triangles(2) = diangles(1) + diangles(2) + diangles(5)
!      triangles(3) = diangles(2) + diangles(3) + diangles(6)
!      triangles(4) = diangles(4) + diangles(5) + diangles(6)

      return
      end

