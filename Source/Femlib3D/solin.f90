      subroutine solin(ok, gilt, eleinfo, epsgl)
      use feminterface, only: destroyarrptr, getsetting
      use femtypes
      use globalvariables3D, only : eltype, nnat, numdof, numv, omega, resdof, vgdof, vp, x
      implicit none
      real (DP) :: epsgl
      logical :: ok, gilt, eleinfo
      intent (out) :: ok, gilt, eleinfo, epsgl
!
!------------------------------------------------------------------------------
!    $Revision: 1.7 $
!    $Date: 2014/07/15 13:00:41 $
!    $Author: m_kasper $
!------------------------------------------------------------------------------
!
!  read in the solution vector from file solution
!
!  Output:
!            ok        =.false. if an error occured
!            gilt      =.true. if the solution is valid
!            eleinfo   =.true. if ep and eg has been written
!            epsgl     accuracy of the solution in the solution of the linear system
!
!  local variables
!
      integer (I4B) :: unitid, i, j, neldof, numvstored, nfnat, inat
      integer (I4B) :: year, mon, day, hour, minute, second
      character(len=3) :: dayofweek, month
      character (len=30) :: datstr
      character (len=200) :: path
      logical :: deadpointer, exist
!------------------------------------------------------------------------------
!
      call getsetting('PROJECTPATH',path)
      inquire (file=path(1:len_trim(path))//'solution',exist=exist)
      if (.not. exist) then
        print*,'***** solution data file not present'
        print*,'***** call ignored'
        ok=.false.
        return
      end if
      call grglun(unitid)
!
!  open the file under the specified path
      open (unitid,file=path(1:len_trim(path))//'solution', &
            form='unformatted',status='old',position='REWIND',action='READ')
!
!  read in the data
      read (unitid,end=1000,err=1000) gilt, eleinfo, numvstored, numdof, omega, epsgl
      read (unitid,end=1000,err=1000) nfnat ! local variable
!
      if (numvstored .gt. 0) then
        if (numvstored .eq. numv) then
          if (nnat .eq. nfnat) then
            if (eleinfo) then
              if (associated(vp)) deallocate(vp)
              allocate(vp(numv,nnat))
              deadpointer=.true.
              if (associated(vgdof)) deadpointer = destroyarrptr(vgdof)
              if (deadpointer) allocate(vgdof(numv,nnat))
              do inat=1,nnat
                do j=1,numv
                  read (unitid,end=1000,err=1000) vp(j,inat)
!  number of Degree Of Freedom (dof) of this element (this is for Nedelec)
                  if ( eltype .eq. 'SCALAR') then
                    neldof = (vp(j,inat)+1)*(vp(j,inat)+2)*(vp(j,inat)+3)/6
                  elseif ( eltype .eq. 'COMPLETE') then
                    neldof = (vp(j,inat)+1)*(vp(j,inat)+2)*(vp(j,inat)+3)/2
! Default is Nedelec type element
                  else
                    neldof = vp(j,inat)*(vp(j,inat)+2)*(vp(j,inat)+3)/2
                  end if
!  allocate array of global dof of the element
                  allocate(vgdof(j,inat)%d(neldof))
                  do i=1,neldof
                  read (unitid,end=1000,err=1000) vgdof(j,inat)%d(i)
                  end do
                end do
              end do
            end if
!
            if (associated(x)) deallocate(x)
            if (associated(resdof)) deallocate(resdof)
            allocate(x(numdof))
            allocate(resdof(numdof))
!  read in the solution (stored for each DOF)
            print*,'numdof read : ', numdof
            do j=1,numdof
              read (unitid,end=1000,err=1000) x(j)
            end do
            do j=1,numdof
              read (unitid,end=1000,err=1000) resdof(j)
            end do
!
          else
            print*
            print*,'**** solution contains vector for a mesh of: ',numvstored,' elements'
            print*,'**** actual mesh size is: ',numv,' elements; ignoring the solution'
          end if
        end if
      end if
      close (unitid,err=1000)
!
      ok=.true.
      return
1000  print*,' an error occured while reading the solution file'
      stop
    end subroutine solin

!   
!
    subroutine solin_aux(ok, gilt, eleinfo, epsgl)
      use feminterface, only: destroyarrptr, getsetting
      use femtypes
      use globalvariables3D, only : eltype, nnat_aux, numdof_aux, omega, vgdof_aux, vp_aux, x_aux, numv
      implicit none
      real (DP) :: epsgl
      logical :: ok, gilt, eleinfo
      intent (out) :: ok, gilt, eleinfo, epsgl
!
!
!  read in an auxiliary solution vector from file solution_aux
!
!  Output:
!            ok        =.false. if an error occured
!            gilt      =.true. if the solution is valid
!            eleinfo   =.true. if ep and eg has been written
!            epsgl     accuracy of the solution in the solution of the linear system
!
!  local variables
!
      integer (I4B) :: unitid, i, j, neldof, numvstored, nfnat, inat
      integer (I4B) :: year, mon, day, hour, minute, second
      character(len=3) :: dayofweek, month
      character (len=30) :: datstr
      character (len=200) :: path, filename
      logical :: deadpointer, exist
!------------------------------------------------------------------------------
!
      call getsetting('PROJECTPATH',path)
      call getsetting('AUX_SOLUTION_FILE', filename)
      inquire (file=path(1:len_trim(path))//filename(1:len_trim(filename)),exist=exist)
      if (.not. exist) then
        print*,'***** solution data file not present'
        print*,'***** call ignored'
        ok=.false.
        return
      end if
      call grglun(unitid)
!
!  open the file under the specified path
      open (unitid,file=path(1:len_trim(path))//filename(1:len_trim(filename)), &
            form='unformatted',status='old',position='REWIND',action='READ')
!
!
      print*,'Reading auxiliary solution from file: ', filename(1:len_trim(filename))
!
!  read in the data
      read (unitid,end=1000,err=1000) gilt, eleinfo, numvstored, numdof_aux, omega, epsgl
      read (unitid,end=1000,err=1000) nfnat ! local variable
!
      if (numvstored .gt. 0) then
        if (numvstored .eq. numv) then
            nnat_aux = nfnat
            if (eleinfo) then
              if (associated(vp_aux)) deallocate(vp_aux)
              allocate(vp_aux(numvstored,nfnat))
              deadpointer=.true.
              if (associated(vgdof_aux)) deadpointer = destroyarrptr(vgdof_aux)
              if (deadpointer) allocate(vgdof_aux(numvstored,nfnat))
              do inat=1,nfnat
                do j=1,numvstored
                  read (unitid,end=1000,err=1000) vp_aux(j,inat)
!  number of Degree Of Freedom (dof) of this element (this is for Nedelec)
                  if ( eltype .eq. 'SCALAR') then
                    neldof = (vp_aux(j,inat)+1)*(vp_aux(j,inat)+2)*(vp_aux(j,inat)+3)/6
                  elseif ( eltype .eq. 'COMPLETE') then
                    neldof = (vp_aux(j,inat)+1)*(vp_aux(j,inat)+2)*(vp_aux(j,inat)+3)/2
! Default is Nedelec type element
                  else
                    neldof = vp_aux(j,inat)*(vp_aux(j,inat)+2)*(vp_aux(j,inat)+3)/2
                  end if
!  allocate array of global dof of the element
                  allocate(vgdof_aux(j,inat)%d(neldof))
                  do i=1,neldof
                    read (unitid,end=1000,err=1000) vgdof_aux(j,inat)%d(i)
                  end do
                end do
              end do
            end if
!
            if (associated(x_aux)) deallocate(x_aux)
            allocate(x_aux(numdof_aux))
!  read in the solution (stored for each DOF)
            print*,'numdof read : ', numdof_aux
            do j=1,numdof_aux
              read (unitid,end=1000,err=1000) x_aux(j)
            end do
!
          else
            print*
            print*,'**** Auxiliary solution contains vector for a mesh of: ',numvstored,' elements'
            print*,'**** Mesh size is: ',numv,' elements; ignoring the solution'
          end if
      end if
      close (unitid,err=1000)
!
      ok=.true.
      return
1000  print*,' an error occured while reading the solution file'
      stop
      end subroutine solin_aux