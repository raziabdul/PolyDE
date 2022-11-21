      subroutine umfsolver_z(a,b,x,n,ir,ic)
      use femtypes
      implicit none
      complex (DPC), pointer :: a(:)
      complex (DPC) b(:),x(:)
      integer (I4B) n
      integer (I4B), pointer :: ir(:),ic(:)
      intent (in) :: n
      intent (inout) :: b, x, ir, ic, a
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
!    $Revision: 1.16 $
!    $Date: 2014/03/07 14:05:34 $
!    $Author: m_kasper $
!
!    a       matrix (coo format)
!    b       right hand side
!    x       solution vector
!    n       order of matrix
!    ir      rox indices (coo format)
!    ic      column indices (coo format)
!
!  local variables
      integer lindex, job, ne, lvalue, d, dn, acopy, luilen
      integer dne, alen, wlen, axcopy
      integer keep(20), info(40), icntl(20)
      integer(8) :: memory, lvalue_l
      integer :: umfmemory_in_gb = 5
      integer, allocatable :: idex(:)
      complex(DPC), allocatable :: value(:), work (:)
      double precision  cntl(10), rinfo(20)

      external :: umz21i, umz2fa, umz2so, umd21i, umd2fa, umd2so
      external :: ums21i, ums2fa, ums2so
      
!  Initialize controls, and change default printing control
      call umz21i(keep, cntl, icntl)
!  Reset default printing control
      icntl(3) = 2
!  set non-zero pattern to be symmetric
      icntl(6) = 1
!  number of matrix entries
      ne=size(ir)
!
!  solve A x = b
!  if job=1, then a column-oriented form of the input matrix is preserved
!  otherwise, the input matrix is overwritten with its LU factors
      job=0
      if (job .eq. 1) then 
        acopy = ne+n+1
        axcopy= ne
      else 
        acopy = 0
        axcopy= 0
      end if
!  d = max (64, sqrt (n)) is the default
      d=max (64, int(sqrt (real(n))))
!  number of columns with more than d entries <=n
      dn=d
!
!  luilen should be keep(5) - keep(3) + 1 - acopy,
!  but keep(1-5) is only known after the call to umz2fa
      luilen=0
!  dne <= ne is the number of entries in such columns with more than d entries
      dne=ne/2+n
!
      alen=2*ne+11*n+11*dn+dne
!  wlen must be <=11*n + 3*dn + 8
      wlen=11*n+3*dn+8
!  lindex must be >3*ne+2*n*1
      lindex=max(3*ne+2*n+1, wlen+alen+acopy, wlen+luilen+acopy)
!
      memory = umfmemory_in_gb * int(2**30,8)           !  available memory in Bytes
! lvalue must be > 2*ne
      lvalue_l = max(11*ne, ne+axcopy, 13*lindex)       !  estimate of needed memory
      lvalue_l = min(lvalue_l, (memory-4*lindex)/16)    !  maximum range (memory limit) DPC=> 16Bytes
      lvalue_l = min(lvalue_l, 2**31-1)                 !  maximum range of 4 Bytes integer
      lvalue = lvalue_l
!
!  assign to value and idex
      allocate (idex(lindex))
!  row indices
      idex(   1:  ne) = ir(1:ne)
!  column indices
      idex(ne+1:2*ne) = ic(1:ne)
      deallocate (ir, ic)
      allocate (value(lvalue))
      value(1:ne) = a(1:ne)
      deallocate (a)

!  Initialize controls
      call umz21i(keep, cntl,icntl)
      icntl(6) = 1
      cntl(2) = 1.05

! Factorize A, Input matrix is not preserved.
      call umz2fa(n, ne, job, .false., lvalue, lindex, value, idex,     &
     &            keep, cntl, icntl, info, rinfo)
      if (info (1) .lt. 0) then
        print*,'Matrix is too large'
        if (info(1) .eq. -3) print*,'lindex too small'
        if (info(1) .eq. -4) print*,'lvalue too small'
        if (info(1) .eq. -5) print*,'lindex and lvalue too small'
        stop
      end if
!      print*,'memory needed in index ',info(19),' size of index',lindex 
!      print*,'memory needed in value ',info(21),' size of value',lvalue
!
! Solve Ax = b
      if (icntl(8) .eq. 0) then
        allocate( work(2*n) )
      else 
        allocate( work(4*n) )
      end if
!
      call umz2so(n, 0, .false., lvalue, lindex, value, idex,           &
     &            keep, b, x, work, cntl, icntl, info, rinfo)
!
      deallocate (work, value, idex )

      return
      end subroutine umfsolver_z



      subroutine umfsolver_d(a,b,x,n,ir,ic)
      use femtypes
      implicit none
      real (DP), pointer :: a(:)
      real (DP) :: b(:),x(:)
      integer (I4B) :: n
      integer (I4B), pointer :: ir(:),ic(:)
      intent (in) :: n
      intent (inout) :: b, x, ir, ic, a
!    a       matrix (coo format)
!    b       right hand side
!    x       solution vector
!    n       order of matrix
!    ir      rox indices (coo format)
!    ic      column indices (coo format)

!  local variables
      integer lindex, job, ne, lvalue, d, dn, acopy, luilen
      integer dne, alen, wlen, axcopy
      integer keep(20), info(40), icntl(20)
      integer(8) :: memory, lvalue_l
      integer :: umfmemory_in_gb = 5
      integer, allocatable :: idex(:)
      double precision, allocatable :: value(:), work (:)
      double precision  cntl(10), rinfo(20)
!
      external :: umz21i, umz2fa, umz2so, umd21i, umd2fa, umd2so
      external :: ums21i, ums2fa, ums2so

!  Initialize controls, and change default printing control
      call umd21i(keep,cntl,icntl)
!  Reset default printing control
      icntl(3) = 2
!  set non-zero pattern to be symmetric
      icntl(6) = 1
!  number of matrix entries
      ne=size(ir)
!
!  solve A x = b
!  if job=1, then a column-oriented form of the input matrix is preserved
!  otherwise, the input matrix is overwritten with its LU factors
      job=0
      if (job .eq. 1) then 
        acopy = ne+n+1
        axcopy= ne
      else 
        acopy = 0
        axcopy= 0
      end if
!  d = max (64, sqrt (n)) is the default
      d=max (64, int(sqrt (real(n))))
!  number of columns with more than d entries <=n
      dn=d
!
!  luilen should be keep(5) - keep(3) + 1 - acopy,
!  but keep(1-5) is only known after the call to umz2fa
      luilen=0
!  dne <= ne is the number of entries in such columns with more than d entries
      dne=ne/2+n
!
      alen=2*ne+11*n+11*dn+dne
!  wlen must be <=11*n + 3*dn + 8
      wlen=11*n+3*dn+8
!  lindex must be >3*ne+2*n*1
      lindex=max(3*ne+2*n+1, wlen+alen+acopy, wlen+luilen+acopy)
!
      memory = umfmemory_in_gb * int(2**30,8)           !  available memory in Bytes
! lvalue must be > 2*ne
      lvalue_l = max(11*ne, ne+axcopy, 13*lindex)       !  estimate of needed memory
      lvalue_l = min(lvalue_l, (memory-4*lindex)/8)     !  maximum range (memory limit)DP=> 8Bytes
      lvalue_l = min(lvalue_l, 2**31-1)                 !  maximum range of 4 Bytes integer
      lvalue = lvalue_l
!
!  assign to value and idex
      allocate (idex(lindex))
!  row indices
      idex(   1:  ne) = ir(1:ne)
!  column indices
      idex(ne+1:2*ne) = ic(1:ne)
      deallocate (ir, ic)
      allocate (value(lvalue))
      value(1:ne) = a(1:ne)
      deallocate (a)

!  Initialize controls
      call umd21i(keep, cntl,icntl)
      icntl(6) = 1
      cntl(2) = 1.05

! Factorize A, Input matrix is not preserved.
      call umd2fa(n,ne,job,.false.,lvalue,lindex,value,idex,keep,cntl,icntl,info,rinfo)
      if (info (1) .lt. 0) then
        print*,'Matrix is too large to solve with UMF'
        if (info(1) .eq. -3) print*,'lindex too small'
        if (info(1) .eq. -4) print*,'lvalue too small'
        if (info(1) .eq. -5) print*,'lindex and lvalue too small'
        stop
      end if
!      print*,'memory needed in index ',info(19),' size of index',lindex
!      print*,'memory needed in value ',info(21),' size of value',lvalue
!
! Solve Ax = b
      if (icntl(8) .eq. 0) then
        allocate( work(2*n) )
      else 
        allocate( work(4*n) )
      end if
!
      call umd2so(n, 0, .false., lvalue, lindex, value, idex, &
     &            keep, b, x, work, cntl, icntl, info, rinfo)
!
      deallocate (work, value, idex )

      return
      end subroutine umfsolver_d



      subroutine umfsolver_s(a,b,x,n,ir,ic)
      use femtypes
      implicit none
      real (SP), pointer :: a(:)
      real (SP) :: b(:),x(:)
      integer (I4B) :: n
      integer (I4B), pointer :: ir(:),ic(:)
      intent (in) :: n
      intent (inout) :: b, x, ir, ic, a
!    a       matrix (coo format)
!    b       right hand side
!    x       solution vector
!    n       order of matrix
!    ir      rox indices (coo format)
!    ic      column indices (coo format)

!  local variables
      integer lindex, job, ne, lvalue, d, dn, acopy, luilen
      integer dne, alen, wlen, axcopy
      integer keep(20), info(40), icntl(20)
      integer(8) :: memory, lvalue_l
      integer :: umfmemory_in_gb = 5
      integer, allocatable :: idex(:)
      real(SP), allocatable :: value(:), work (:)
      real(SP)  cntl(10), rinfo(20)
!
      external :: umz21i, umz2fa, umz2so, umd21i, umd2fa, umd2so
      external :: ums21i, ums2fa, ums2so

!  Initialize controls, and change default printing control
      call umd21i(keep,cntl,icntl)
!  Reset default printing control
      icntl(3) = 2
!  set non-zero pattern to be symmetric
      icntl(6) = 1
!  number of matrix entries
      ne=size(ir)
!
!  solve A x = b
!  if job=1, then a column-oriented form of the input matrix is preserved
!  otherwise, the input matrix is overwritten with its LU factors
      job=0
      if (job .eq. 1) then 
        acopy = ne+n+1
        axcopy= ne
      else 
        acopy = 0
        axcopy= 0
      end if
!  d = max (64, sqrt (n)) is the default
      d=max (64, int(sqrt (real(n))))
!  number of columns with more than d entries <=n
      dn=d
!
!  luilen should be keep(5) - keep(3) + 1 - acopy, 
!  but keep(1-5) is only known after the call to umz2fa 
      luilen=0
!  dne <= ne is the number of entries in such columns with more than d entries
      dne=ne/2+n
!
      alen=2*ne+11*n+11*dn+dne
!  wlen must be <=11*n + 3*dn + 8
      wlen=11*n+3*dn+8
!  lindex must be >3*ne+2*n*1
      lindex=max(3*ne+2*n+1, wlen+alen+acopy, wlen+luilen+acopy)
!
      memory = umfmemory_in_gb * int(2**30,8)           !  available memory in Bytes
! lvalue must be > 2*ne
      lvalue_l = max(11*ne, ne+axcopy, 13*lindex)       !  estimate of needed memory
      lvalue_l = min(lvalue_l, (memory-4*lindex)/4)     !  maximum range (memory limit)DP=> 4Bytes
      lvalue_l = min(lvalue_l, 2**31-1)                 !  maximum range of 4 Bytes integer
      lvalue = lvalue_l
!
!  assign to value and idex
      allocate (idex(lindex))
!  row indices
      idex(   1:  ne) = ir(1:ne)
!  column indices
      idex(ne+1:2*ne) = ic(1:ne)
      deallocate (ir, ic)
      allocate (value(lvalue))
      value(1:ne) = a(1:ne)
      deallocate (a)

!  Initialize controls
      call ums21i(keep, cntl,icntl)
      icntl(6) = 1
      cntl(2) = 1.05

! Factorize A, Input matrix is not preserved.
      call ums2fa(n,ne,job,.false.,lvalue,lindex,value,idex,keep,cntl,icntl,info,rinfo)
      if (info (1) .lt. 0) then
        print*,'Matrix is too large to solve with UMF'
        if (info(1) .eq. -3) print*,'lindex too small'
        if (info(1) .eq. -4) print*,'lvalue too small'
        if (info(1) .eq. -5) print*,'lindex and lvalue too small'
        stop
      end if
!      print*,'memory needed in index ',info(19),' size of index',lindex
!      print*,'memory needed in value ',info(21),' size of value',lvalue
!
! Solve Ax = b
      if (icntl(8) .eq. 0) then
        allocate( work(2*n) )
      else 
        allocate( work(4*n) )
      end if
!
      call ums2so(n, 0, .false., lvalue, lindex, value, idex, &
     &            keep, b, x, work, cntl, icntl, info, rinfo)
!
      deallocate (work, value, idex )

      return
      end subroutine umfsolver_s
