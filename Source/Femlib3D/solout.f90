      subroutine solout(gilt, eleinfo, epsgl)
      use feminterface, only: getsetting
      use femtypes
      use globalvariables3D, only : nnat, numdof, numv, omega, resdof, vgdof, vp, x
      implicit none
      real (DP) :: epsgl
      logical :: gilt, eleinfo
      intent (in) :: gilt, eleinfo, epsgl
!
!------------------------------------------------------------------------------
!    $Revision: 1.6 $
!    $Date: 2015/04/01 16:56:08 $
!    $Author: juryanatzki $
!------------------------------------------------------------------------------
!
!  write the solution to the file solution
!
!  Input
!            eleinfo   =.true. if ep and eg should be written
!            epsgl     accuracy of the solution in the solution of the linear system
!            gilt      =.true. if the solution is valid
!
!  local variables
      integer (I4B) unitid, i, j, inat
      real (DP) jd
      character (len=200) path
!------------------------------------------------------------------------------
!
      call getsetting('PROJECTPATH',path)
      call grglun(unitid)
!
!  open the file loes in the specified path
      open (unitid,file=path(1:len_trim(path))//'solution', &
            form='unformatted',position='REWIND',action='WRITE')
!
!  write the data to the file
      write (unitid,err=1000) gilt, eleinfo, numv, numdof, omega, epsgl
!
      if (eleinfo) then
         write (unitid,err=1000) nnat
         do inat=1,nnat
           do j=1,numv
             write (unitid,err=1000) vp(j,inat)
!  number of Degree Of Freedom (dof) of this element
             do i=1,size(vgdof(j,inat)%d)
               write (unitid,err=1000) vgdof(j,inat)%d(i)
             end do
           end do
         end do
      end if
!
      if (gilt .and. associated(x)) then
         do j=1,numdof
            write (unitid,err=1000) x(j)
      !      print*,'x(',j,') = ', x(j)
         end do
         do j=1,numdof
            write (unitid,err=1000) 0._DP
            !write (unitid,err=1000) resdof(j)
         end do
      else
         do j=1,numdof
            write (unitid,err=1000) cmplx(0._DP,0._DP,DPC)
         end do
         do j=1,numdof
            write (unitid,err=1000) 0._DP
         end do
      end if
!
      close (unitid,err=1000)
      return
1000  print*,' an error occured while writing the solution vector'
      stop
      end subroutine solout