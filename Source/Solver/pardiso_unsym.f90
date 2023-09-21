      subroutine  pardiso_solver_DPC(a,b,x,n,ia,ja)
      use feminterface, only: qsortindex, getsetting
      use femtypes
      use mkl_pardiso
      use MKL_PARDISO_PRIVATE,      only: MKL_PARDISO_HANDLE
      implicit none
      complex (DPC), pointer :: a(:), b(:)
      complex (DPC) x(:)
      integer (I4B), pointer :: ia(:),ja(:)
      integer (I4B) n
      intent (in) :: n
      intent (out) :: x
!
!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
!
!  local variables:
!  Internal  solver  memory  pointer  for  64-bit  architectures
!  INTEGER*8  pt(64) -> KIND=8  !(work also for 32-bit)
!  Internal  solver  memory  pointer  for  32-bit  architectures
!  INTEGER*4  pt(64) -> KIND=4
      TYPE(MKL_PARDISO_HANDLE), allocatable :: PT(:)
! null addresses for unused input data
      INTEGER, pointer :: idum(:) => null()
      REAL(KIND=8), pointer :: ddum(:) => null()
      COMPLEX(DPC), pointer :: cdum(:) => null()

!  All  other  variables
      integer (I4B)                    :: solvermemory_in_gb
      real (DP)                        :: peak_memory

      INTEGER  maxfct,  mnum,  mtype,  phase,  nrhs,  error,  msglvl
      INTEGER  iparm(64)
      INTEGER  i, j

! For timers
    !real(KIND=8)  waltime1,  waltime2


      integer nnz_row_i, rowstart, rowend
      integer (I4B), allocatable :: indx(:)

      complex (DPC), allocatable :: as(:)
      integer (I4B), allocatable :: jas(:)
      real (DP), allocatable :: array(:)

! Some Pardiso paramenters
      DATA nrhs /1/,  maxfct /1/,  mnum /1/

! Set limit for peak memory for the solver
      print *, 'PARDISO with complex data type'
      call getsetting('LINSOLVER_PEAKMEM_GB', solvermemory_in_gb)

    ! revert the Polyde CSR to
    ! Pardiso standard CSR (ja indices in increasing order)
      allocate( as(SIZE(a)) )
      allocate( jas(SIZE(ja)) )

      print*,'ia size',size(ia)
      do i=1,size(ia)-1
        nnz_row_i = ia(i+1) - ia(i)
        ! take a row slice from ja 
        allocate(array(nnz_row_i))

        rowstart=ia(i)
        rowend=ia(i+1)-1
        ! set slice as negative to get indx in increasing order
        ! qsortindex takes a real argument
        array = real( -ja(rowstart:rowend), (DP) )
        !print*,'slice ', array

        allocate(indx(nnz_row_i))
        call qsortindex( array, indx, nnz_row_i)

        do j=rowstart,rowend
          ! sort a and ja through ordered indx
          as(j)  = a(indx(j-rowstart +1 ) + rowstart-1)
          jas(j) = ja(indx(j-rowstart +1 ) + rowstart-1)
        end do
        !print*,'ordered ', jas(rowstart:rowend)

        deallocate(indx)
        deallocate(array)
      end do

      print*,'Preparing PARDISO...'
    
! Clean handler and dummy addresses
      allocate( PT(64), idum(0), ddum(0), cdum(0))
      do i = 1, 64
        PT(i)%DUMMY = 0
      end do

!  Set  up  PARDISO  control  parameter
      do  i  =  1,  64
        iparm(i)  =  0
      end  do
      iparm(1)  =  1  !  no  solver  default
      iparm(2)  =  3  !  fill-in  reordering  from  METIS (OpenMP nested bissection)
      iparm(3)  =  1  !  numbers  of  processors, serial compilation
      iparm(4)  =  0  !  no  iterative-direct  algorithm
      iparm(5)  =  0  !  no  user  fill-in  reducing  permutation
      iparm(6)  =  0  !  =0  solution  on  the  first  n  components  of  x
      iparm(8)  =  9  !  numbers  of  iterative  refinement  steps
      iparm(10)  =  13  !  perturb  the  pivot  elements  with  1E-13
      iparm(11)  =  1  !  use  nonsymmetric  permutation  and  scaling  MPS
      iparm(14)  =  0  !  Output:  number  of  perturbed  pivots
      iparm(18)  =  -1  !  Output:  number  of  nonzeros  in  the  factor  LU
      iparm(19)  =  -1  !  Output:  Mflops  for  LU  factorization
      iparm(20)  =  0  !  Output:  Numbers  of  CG  Iterations
      error  =  0  !  initialize  error  flag
      msglvl  =  0  !  print  statistical  information
      !mtype  =  11  !  real  unsymmetric
      mtype  =  13  !  complex  unsymmetric

!  Reordering  and  Symbolic  Factorization,  This  step  also  allocates
!  all  memory  that  is  necessary  for  the  factorization
      phase  =  11  !  only  reordering  and  symbolic  factorization
      WRITE(*,*)  'Ready for reordering  and  symbolic  factorization ...  '
      CALL  pardiso  (PT,  maxfct,  mnum,  mtype,  phase,  n,  as,  ia, &
      &  jas, idum,  nrhs,  iparm,  msglvl,  cdum,  cdum,  error)
      WRITE(*,*)  'Reordering  completed  ...  '
      IF  (error  .NE.  0)  THEN
        WRITE(*,*)  'The  following symb. factor. ERROR  was  detected:  ',  error
        STOP
      END  IF
      WRITE(*,*)  'Number  of  nonzeros  in  factors  =  ',iparm(18)
      WRITE(*,*)  'Number  of  factorization  MFLOPS  =  ',iparm(19)

      peak_memory = real(max(iparm(15),iparm(16)+iparm(17)),DP)/(1024*1024)
      print *, '>>>>> Peak Memory', peak_memory
      if (peak_memory .gt. solvermemory_in_gb) then
        print *,' !!! Insufficient memory for the pardiso solver !!!'
      end if

!  Factorization.
    !WRITE(*,*)  'Ready for numerical factorization...  '
      phase  =  22  !  only  factorization
      CALL  pardiso  (pt,  maxfct,  mnum,  mtype,  phase,  n,  as,  ia,  &
      &  jas, idum,  nrhs,  iparm,  msglvl,  cdum,  cdum,  error)
      WRITE(*,*)  'Factorization  completed  ...  '
      IF  (error  .NE.  0)  THEN
        WRITE(*,*)  'The  following  num. factor. ERROR  was  detected:  ',  error
        STOP
      ENDIF

!  Back  substitution  and  iterative  refinement
    !WRITE(*,*)  'Reasy to to Back/Forward substitution (solve)...  '
      msglvl  =  0  !  print  statistical  information
      iparm(8)  =  2  !  max  numbers  of  iterative  refinement  steps
      phase  =  33  !  only  factorization
      CALL  pardiso  (pt,  maxfct,  mnum,  mtype,  phase,  n,  as,  ia, &
      &  jas, idum,  nrhs,  iparm,  msglvl,  b,  x,  error)
      IF  (error  .NE.  0)  THEN
        WRITE(*,*)  'The  following  solution ERROR  was  detected:  ',  error
        STOP
      ENDIF
      WRITE(*,*)  'Solve  completed  ...  '


!  Termination  and  release  of  memory
    !WRITE(*,*)  'Ready to terminate/release memory...  '
      phase  =  -1  !  release  internal  memory
      CALL  pardiso  (pt,  maxfct,  mnum,  mtype,  phase,  n,  cdum,  &
      &  idum,  idum, idum,  nrhs,  iparm,  msglvl,  cdum,  cdum,  error)
      IF  (error  .NE.  0)  THEN
        WRITE(*,*)  'The  following termination ERROR  was  detected:  ',  error
        STOP
      ENDIF

    !kill sorted matrix /a and column pointer /ja
      deallocate( as )
      deallocate( jas )

      deallocate( PT, idum, ddum, cdum)


      return
      end subroutine
!
!-------------------------------------------------------------------------------
!  
      subroutine  pardiso_solver_DP(a,b,x,n,ia,ja)
        use feminterface, only: qsortindex, getsetting, sort_asc_order
        use femtypes
        use mkl_pardiso
        use MKL_PARDISO_PRIVATE,      only: MKL_PARDISO_HANDLE
        implicit none
        real (DP), pointer :: a(:), b(:)
        real (DP) x(:)
        integer (I4B), pointer :: ia(:),ja(:)
        integer (I4B) n
        intent (in) :: n
        intent (out) :: x
  !
  !-------------------------------------------------------------------------------
  !
  !  local variables:
        TYPE(MKL_PARDISO_HANDLE), allocatable :: PT(:)
  ! null addresses for unused input data
        INTEGER, pointer :: idum(:) => null()
        REAL(KIND=8), pointer :: ddum(:) => null()
        real (DP), pointer :: cdum(:) => null()
  
  !  All  other  variables
        integer (I4B)                    :: solvermemory_in_gb
        real (DP)                        :: peak_memory
  
        INTEGER  maxfct,  mnum,  mtype,  phase,  nrhs,  error,  msglvl
        INTEGER  iparm(64)
        INTEGER  i, j
  
  ! For timers
      !real(KIND=8)  waltime1,  waltime2
  
  
        integer nnz_row_i, rowstart, rowend
        integer (I4B), allocatable :: indx(:)
    
  ! Some Pardiso paramenters
        DATA nrhs /1/,  maxfct /1/,  mnum /1/
  
  ! Set limit for peak memory for the solver
        print *, 'PARDISO with real data type'
        call getsetting('LINSOLVER_PEAKMEM_GB', solvermemory_in_gb)
  
      ! revert the Polyde CSR to
      ! Pardiso standard CSR (ja indices in increasing order)
        call sort_asc_order(a, ia, ja)
  
        print*,'Preparing PARDISO...'
      
  ! Clean handler and dummy addresses
        allocate( PT(64), idum(0), ddum(0), cdum(0))
        do i = 1, 64
          PT(i)%DUMMY = 0
        end do
  
  !  Set  up  PARDISO  control  parameter
        do  i  =  1,  64
          iparm(i)  =  0
        end  do
        iparm(1)  =  1  !  no  solver  default
        iparm(2)  =  3  !  fill-in  reordering  from  METIS (OpenMP nested bissection)
        iparm(3)  =  1  !  numbers  of  processors, serial compilation
        iparm(4)  =  0  !  no  iterative-direct  algorithm
        iparm(5)  =  0  !  no  user  fill-in  reducing  permutation
        iparm(6)  =  0  !  =0  solution  on  the  first  n  components  of  x
        iparm(8)  =  9  !  numbers  of  iterative  refinement  steps
        iparm(10)  =  13  !  perturb  the  pivot  elements  with  1E-13
        iparm(11)  =  1  !  use  nonsymmetric  permutation  and  scaling  MPS
        iparm(14)  =  0  !  Output:  number  of  perturbed  pivots
        iparm(18)  =  -1  !  Output:  number  of  nonzeros  in  the  factor  LU
        iparm(19)  =  -1  !  Output:  Mflops  for  LU  factorization
        iparm(20)  =  0  !  Output:  Numbers  of  CG  Iterations
        error  =  0  !  initialize  error  flag
        msglvl  =  0  !  print  statistical  information
        mtype  =  11  !  real  unsymmetric
  
  !  Reordering  and  Symbolic  Factorization,  This  step  also  allocates
  !  all  memory  that  is  necessary  for  the  factorization
        phase  =  11  !  only  reordering  and  symbolic  factorization
        WRITE(*,*)  'Ready for reordering  and  symbolic  factorization ...  '

        CALL  pardiso  (PT,  maxfct,  mnum,  mtype,  phase,  n,  a,  ia, &
        &  ja, idum,  nrhs,  iparm,  msglvl,  cdum,  cdum,  error)

        WRITE(*,*)  'Reordering  completed  ...  '
        IF  (error  .NE.  0)  THEN
          WRITE(*,*)  'The  following symb. factor. ERROR  was  detected:  ',  error
          STOP
        END  IF
        WRITE(*,*)  'Number  of  nonzeros  in  factors  =  ',iparm(18)
        WRITE(*,*)  'Number  of  factorization  MFLOPS  =  ',iparm(19)
  
        peak_memory = real(max(iparm(15),iparm(16)+iparm(17)),DP)/(1024*1024)
        print *, '>>>>> Peak Memory', peak_memory
        if (peak_memory .gt. solvermemory_in_gb) then
          print *,' !!! Insufficient memory for the pardiso solver !!!'
        end if
  
  !  Factorization.
      !WRITE(*,*)  'Ready for numerical factorization...  '
        phase  =  22  !  only  factorization
        CALL  pardiso  (pt,  maxfct,  mnum,  mtype,  phase,  n,  a,  ia,  &
        &  ja, idum,  nrhs,  iparm,  msglvl,  cdum,  cdum,  error)
        WRITE(*,*)  'Factorization  completed  ...  '
        IF  (error  .NE.  0)  THEN
          WRITE(*,*)  'The  following  num. factor. ERROR  was  detected:  ',  error
          STOP
        ENDIF
  
  !  Back  substitution  and  iterative  refinement
      !WRITE(*,*)  'Reasy to to Back/Forward substitution (solve)...  '
        msglvl  =  0  !  print  statistical  information
        iparm(8)  =  2  !  max  numbers  of  iterative  refinement  steps
        phase  =  33  !  Solve
        CALL  pardiso  (pt,  maxfct,  mnum,  mtype,  phase,  n,  a,  ia, &
        &  ja, idum,  nrhs,  iparm,  msglvl,  b,  x,  error)
        IF  (error  .NE.  0)  THEN
          WRITE(*,*)  'The  following  solution ERROR  was  detected:  ',  error
          STOP
        ENDIF
        WRITE(*,*)  'Solve  completed  ...  '
  
  
  !  Termination  and  release  of  memory
      !WRITE(*,*)  'Ready to terminate/release memory...  '
        phase  =  -1  !  release  internal  memory
        CALL  pardiso  (pt,  maxfct,  mnum,  mtype,  phase,  n,  cdum,  &
        &  idum,  idum, idum,  nrhs,  iparm,  msglvl,  cdum,  cdum,  error)
        IF  (error  .NE.  0)  THEN
          WRITE(*,*)  'The  following termination ERROR  was  detected:  ',  error
          STOP
        ENDIF
  
        deallocate( PT, idum, ddum, cdum)
  
  
        return
        end subroutine