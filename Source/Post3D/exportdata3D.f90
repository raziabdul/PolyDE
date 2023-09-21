      subroutine exportdata3D(fieldtype,phi,origin,p1,p2,grid1,grid2)
      use feminterface, only: getsetting
      use feminterface3D, only: field3D, findelement3D, xyz2lam
      use femtypes
      use globalvariables3D, only: nod, numv, pi, vn, vv
      implicit none
      integer (I4B) :: grid1, grid2
      real (DP) :: phi, origin(3), p1(3), p2(3)
      character (len=10) :: fieldtype
      intent (in) :: fieldtype, phi, origin, p1, p2, grid1, grid2
!
!
!  Subroutine to export values for given fieldtype and phase as grid data on a
!  regular grid of grid1 * grid2 points to file data_<fieldtype><n>.txt, where
!  <n> is an integer. The file contains the fieldtype, origin of data and number
!  of grid points in first and second direction (from origin), p1 coords and
!  distances (val1), p2 coords and distances (val2) and matrix of grid values
!  (zval).
!
!  Comments in file additionally include date and time of solution.
!
!
!  Input:
!     fieldtype          field quantity to compute (e.g. POTENTIAL, EY ...)
!     phi                phase / time instant to compute
!     grid1, grid2       number of grid points in p1- and p2-direction
!                        (relative to origin)

!  Internal variables:
      integer (I4B) :: elem, i, ierror, j, unitid
      real (DP) :: dist1, dist2, vec1(3), vec2(3), unitvec1(3), unitvec2(3)
      real (DP) :: delta(2), lambda(4), point(3)
      real (DP), allocatable :: val1(:), val2(:), zval(:,:)
      character (len=8)  :: dat
      character (len=10) :: tim, tmp, tmp1
      character (len=23) :: fmtreal
      character (len=200):: path
      complex (DPC) :: u(3), curlu(3)
      logical :: found
!
      external :: grglun
!------------------------------------------------------------------------------
!  compute p1 and p2 direction
      vec1 = p1 - origin
      vec2 = p2 - origin
!  compute length of vector and unitvector (length 1)
      dist1 = sqrt(vec1(1)**2 + vec1(2)**2 + vec1(3)**2)
      dist2 = sqrt(vec2(1)**2 + vec2(2)**2 + vec2(3)**2)
      unitvec1 = vec1/dist1
      unitvec2 = vec2/dist2
!  check for orthogonal vectors
      if (dot_product(unitvec1,unitvec2) .ne. 0._DP) then
        print *,"*****Vectors spanning plane are not orthogonal."
        stop
      end if
!
!  compute distance between neighbouring grid points
      delta(1) = dist1/(grid1-1)
      delta(2) = dist2/(grid2-1)
!  allocate coordinate vectors and solution value matrix
      allocate(val1(grid1), val2(grid2), zval(grid1,grid2))
      val1 = (/((i-1)*delta(1) ,i=1, grid1)/)
      val2 = (/((i-1)*delta(2) ,i=1, grid2)/)
!  loop over all grid points to compute fieldvalue
      do i = 1, size(val1)
        do j = 1, size(val2)
!  get element for point
          point = origin + val1(i)*unitvec1 + val2(j)*unitvec2
          call findelement3D(point,nod,vn,vv,numv,elem,found)
!  if point is in an element, compute value
          if (found) then
            call xyz2lam(lambda, elem, point, nod, vn)
            call field3D(elem,lambda,u,curlu)
            select case(fieldtype)
            case ('UX')
              zval(i,j) = real(u(1)*exp(cmplx(0._DP,phi*pi/180._DP,DPC)),DP)
            case ('UY')
              zval(i,j) = real(u(2)*exp(cmplx(0._DP,phi*pi/180._DP,DPC)),DP)
            case ('UZ')
              zval(i,j) = real(u(3)*exp(cmplx(0._DP,phi*pi/180._DP,DPC)),DP)
            case ('CURLUX')
              zval(i,j) = real(curlu(1)*exp(cmplx(0._DP,phi*pi/180._DP,DPC)),DP)
            case ('CURLUY')
              zval(i,j) = real(curlu(2)*exp(cmplx(0._DP,phi*pi/180._DP,DPC)),DP)
            case ('CURLUZ')
              zval(i,j) = real(curlu(3)*exp(cmplx(0._DP,phi*pi/180._DP,DPC)),DP)
            case default
              print *,"*****Fieldtype ",trim(fieldtype)," not found!"
              stop
            end select
!  if point is not in an element, set value to 0
          else
            zval(i,j) = 0._DP
          end if
        end do
      end do
!
!  write data to file data<n>.txt
      call grglun(unitid)
      call getsetting('PROJECTPATH',path)
!  open a file in the specified project path
      open (unitid,FILE=path(1:len_trim(path))//'data_'//trim(fieldtype)//'.txt',STATUS='NEW',&
            FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
      if (ierror .ne. 0) then
        do i = 1,99
          write (tmp ,'(i2)') i
          open (unitid,FILE=path(1:len_trim(path))//'data_'//trim(fieldtype)//trim(adjustl(tmp))//'.txt',&
                STATUS='NEW',FORM='formatted',POSITION='REWIND',ACTION='WRITE',IOSTAT=ierror)
          if (ierror .ne. 0) then
            cycle
          else
            print '(a)',"Writing output to "//path(1:len_trim(path))//'data_'//trim(fieldtype)//trim(adjustl(tmp))//'.txt'
            exit
          end if
        end do
      else
        print '(a)',"Writing output to "//path(1:len_trim(path))//'data_'//trim(fieldtype)//'.txt'
      end if
!
!  write fieldtype of data contained in file
      write (unitid,'(a,/2x,a)') '% Fieldtype',fieldtype

!  write location of origin
      write (unitid,'(a)') '% Origin of data'
      write (unitid,'(3(1x,es15.8))') origin(1), origin(2), origin(3)

!  write location of p1 point
      write (unitid,'(a)') '% p1 point of data'
      write (unitid,'(3(1x,es15.8))') p1(1), p1(2), p1(3)

!  write location of p2 point
      write (unitid,'(a)') '% p2per point of data'
      write (unitid,'(3(1x,es15.8))') p2(1), p2(2), p2(3)

!  write number of grid points in p1- and p2-direction
      write (unitid,'(a)') '% Number of grid points in p1- and p2-direction'
      write (tmp ,'(i10)') size(val1)
      write (tmp1,'(i10)') size(val2)
      write (unitid,'(2x,a)') trim(adjustl(tmp))//' '//trim(adjustl(tmp1))

!  write header for grid of p1- and p2-distances
      write (unitid,'(a)') '% Grid of '//trim(adjustl(tmp)) // &
                             ' p1- and '//trim(adjustl(tmp1)) // &
                             ' p2-values'
!  write p1-values
      fmtreal = '('//trim(adjustl(tmp))//'(1x,es15.8)'//')'
      write (unitid,trim(fmtreal)) val1
!  write p2-values
      fmtreal = '('//trim(adjustl(tmp1))//'(1x,es15.8)'//')'
      write (unitid,trim(fmtreal)) val2
!  write date and time of solution
      call date_and_time(dat,tim)
      write (unitid,'(a)') '% Data ('//trim(fieldtype)//')' // &
                           ', Date = '//dat(7:8)//'.'//dat(5:6)//'.'//dat(1:4) // &
                           ', Time = '//tim(1:2)//':'//tim(3:4)
!  write field values as a grid1 * grid2 matrix
      write (tmp,'(i10)') size(zval(:,1))
      fmtreal = '('//trim(adjustl(tmp))//'(1x,es15.8)'//')'
      do i = 1, size(zval(1,:))
        write (unitid,trim(fmtreal)) zval(:,i)
      end do
!
      close (unitid)
      call grglun(unitid)
      deallocate (val1,val2,zval)
!
      end subroutine