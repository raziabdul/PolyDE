      subroutine lcsrzeroremover_ZZ(diag, lower, upper, lia, lja, symm, method, eps, rhs, x, eps_x)
      use femtypes
      use feminterface, only: reallocate
      implicit none
      complex(DPC), pointer :: diag(:), lower(:), upper(:)
      complex(DPC), pointer, optional :: x(:), rhs(:)
      real(DP), optional :: eps, eps_x
      integer(I4B), pointer :: lia(:), lja(:)
      integer (I4B) :: method
      logical :: symm
      intent(in) :: symm, method, eps_x, x
      intent(inout) :: diag, lower, upper, lia, lja, rhs


      integer (I4B) :: row, k, del, n, lnnz, start
      real (DP) :: delta 
print*,'function not implemented'
stop
      return
      end



      subroutine lcsrzeroremover_CC(diag, lower, upper, lia, lja, symm, method, eps, rhs, x, eps_x)
      use femtypes
      use feminterface, only: reallocate
      implicit none
      complex(SPC), pointer :: diag(:), lower(:), upper(:)
      complex(SPC), pointer, optional :: x(:), rhs(:)
      real(DP), optional :: eps, eps_x
      integer(I4B), pointer :: lia(:), lja(:)
      integer (I4B) :: method
      logical :: symm
      intent(in) :: symm, method, eps_x, x
      intent(inout) :: diag, lower, upper, lia, lja, rhs

      integer (I4B) :: row, k, del, n, lnnz, start
      real (DP) :: delta 
print*,'function not implemented'
stop

      return
      end



      subroutine lcsrzeroremover_DD(diag, lower, upper, lia, lja, symm, method, eps, rhs, x, eps_x)
      use femtypes
      use feminterface, only: reallocate, a_norm, f_norm
      implicit none
      real(DP), pointer :: diag(:), lower(:), upper(:)
      real(DP), pointer, optional :: x(:), rhs(:)
      real(DP), optional :: eps, eps_x
      integer(I4B), pointer :: lia(:), lja(:)
      integer (I4B) :: method
      logical :: symm
      intent(in) :: symm, method, eps_x, x
      intent(inout) :: diag, lower, upper, lia, lja, rhs
!
!  sparse matrix compactation by dropping small (off-diagonal) matrix entries 
!  if they are irrelevant for the solution
!
!    diag        diagonal of the matrix
!    lower       lower triangular matrix (compact storage)
!    upper       upper triangular matrix (compact storage)
!    b           right hand side vector
!    lia         compact storage information (LCSR format)
!    lja         compact storage information (LCSR format)
!    symm        = .true. if matrix is symmetric
!    method      method to decide whether a matrix entry can be dropped
!    eps         required accuracy (of solution)
!    x           solution vector            (optional)
!    eps_x       accuracy of approximate solution x provided

      integer (I4B) :: row, col, k, drops, n, lnnz, start
      real (DP) :: delta, alpha, x_anrm, b_nrm, x_nrm
      real (DP), allocatable :: ax(:)

      n = size(diag)                                             ! number of rows
      lnnz = size(lower)                                         ! number of non-zero entries in lower (or upper) matrix

alpha=1.0_DP

      select case(method)

      case (0)
!  x is not present
        delta = 1e2_DP * tiny(1._DP)
!  start dropping
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k)) .lt. delta .and.   &
     &          abs(upper(k)) .lt. delta ) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do

      case (1)
!  x is not present
        if (present(eps)) then
          delta = max(epsilon(1._DP),eps) / real(lnnz,DP)
        else
          delta = epsilon(1._DP) / real(lnnz,DP)
        end if
!  calculate row sum (l1-norm)
        allocate(ax(n))
        ax = abs(diag)
        do row = 2, n
          do k = lia(row-1)+1, lia(row)
            col = lja(k)
            ax(row) = ax(row) + abs(lower(k))
            ax(col) = ax(col) + abs(upper(k))
          end do
        end do
!  start dropping
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k)) .lt. delta*ax(row) .and.   &
     &          abs(upper(k)) .lt. delta*ax(col) ) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do
        deallocate(ax)

      case (2)
!  x is not present
        if (present(eps)) then
          delta = max(epsilon(1._DP),eps) / real(lnnz,DP)
        else
          delta = epsilon(1._DP) / real(lnnz,DP)
        end if
!
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k)*upper(k)/(diag(col)*diag(row))) .lt. delta**2) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do

      case (3)
!  rhs and x are present
        b_nrm = sqrt(dot_product(rhs,rhs))
        x_nrm = sqrt(dot_product(x,x))
        if (present(eps)) then
          delta = max(epsilon(1._DP),eps) / real(lnnz,DP)
        else
          delta = epsilon(1._DP) / real(lnnz,DP)
        end if
!
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k))*x_nrm*eps_x .lt. delta*b_nrm .and. &
     &          abs(upper(k))*x_nrm*eps_x .lt. delta*b_nrm) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha

              rhs(row) = rhs(row) - lower(k)*x(col) + abs(lower(k))*x(row)*alpha
              rhs(col) = rhs(col) - upper(k)*x(row) + abs(upper(k))*x(col)*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do

      case (4)
!  rhs and x are present
        x_anrm = a_norm(x,diag,lower,upper,lia,lja)
        eps = max(epsilon(1._DP),eps)

        delta = (eps_x + eps) / eps_x  * ((1._DP-eps_x)/(1._DP+eps_x))**2 / real(lnnz,DP) * x_anrm
!
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k)*x(col)*x(row)) .lt. delta .and. &
     &          abs(upper(k)*x(col)*x(row)) .lt. delta) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha

              rhs(row) = rhs(row) - lower(k)*x(col) + abs(lower(k))*x(row)*alpha
              rhs(col) = rhs(col) - upper(k)*x(row) + abs(upper(k))*x(col)*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do

      case (5)
!  rhs and x are present
        x_anrm = a_norm(x,diag,lower,upper,lia,lja)
        if (present(eps)) then
          delta = max(epsilon(1._DP),eps)**2 / real(lnnz,DP) * x_anrm /eps_x
        else
          delta = epsilon(1._DP) / real(lnnz,DP) * x_anrm/2._DP/(1+2*eps_x)
        !  delta = max(epsilon(1._DP),eps) / real(lnnz,DP) * x_anrm/2._DP/(1+2*eps_x)
        !else
        !  delta = epsilon(1._DP) / real(lnnz,DP) * x_anrm/2._DP/(1+2*eps_x)
        end if
!
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k)*x(col)*x(row)) .lt. delta .and. &
     &          abs(upper(k)*x(col)*x(row)) .lt. delta) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha

              rhs(row) = rhs(row) - lower(k)*x(col) + abs(lower(k))*x(row)*alpha
              rhs(col) = rhs(col) - upper(k)*x(row) + abs(upper(k))*x(col)*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do

      case (6)
!  rhs and x are present
        if (present(eps)) then
          delta = max(epsilon(1._DP),eps) / real(lnnz,DP)
        else
          delta = epsilon(1._DP) / real(lnnz,DP)
        end if
        delta=delta / eps_x
!  calculate row sum (l1-norm)
        allocate(ax(n))
        ax = abs(diag*x)
        do row = 2, n
          do k = lia(row-1)+1, lia(row)
            col = lja(k)
            ax(row) = ax(row) + abs(lower(k)*x(row))
            ax(col) = ax(col) + abs(upper(k)*x(col))
          end do
        end do

!
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k)*x(row)) .lt. delta*ax(row) .and. &
     &          abs(upper(k)*x(col)) .lt. delta*ax(col)) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha

              rhs(row) = rhs(row) - lower(k)*x(col) + abs(lower(k))*x(row)*alpha
              rhs(col) = rhs(col) - upper(k)*x(row) + abs(upper(k))*x(col)*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do
        deallocate(ax)


      end select

      lja => reallocate(lja,lnnz-drops)
      lower => reallocate(lower,lnnz-drops)
      if (symm) then 
        upper => lower
      else
        upper => reallocate(upper,lnnz-drops)
      end if

print '(a,f10.2,a)', 'removed', 100.*real(drops)/real(lnnz),' % matrix entries'
      return
      end 



      subroutine lcsrzeroremover_SS(diag, lower, upper, lia, lja, symm, method, eps, rhs, x, eps_x)
      use femtypes
      use feminterface, only: reallocate, a_norm, f_norm
      implicit none
      real(SP), pointer :: diag(:), lower(:), upper(:)
      real(SP), pointer, optional :: x(:), rhs(:)
      real(DP), optional :: eps, eps_x
      integer(I4B), pointer :: lia(:), lja(:)
      integer (I4B) :: method
      logical :: symm
      intent(in) :: symm, method
      intent(inout) :: diag, lower, upper, lia, lja

      integer (I4B) :: row, col, k, drops, n, lnnz, start
      real (SP) :: delta, alpha, x_anrm
      real (SP), allocatable :: ax(:)

      n = size(diag)
      lnnz = size(lower)                                    ! number of non-zero entries in lower (or upper) matrix

alpha=1.0_SP

      select case(method)

      case (1)
!  x is not present
        if (present(eps)) then
          delta = max(epsilon(1._SP),eps) / real(lnnz,SP)
        else
          delta = epsilon(1._SP) / real(lnnz,SP)
        end if
!  calculate li row sum
        allocate(ax(n))
        ax = abs(diag)
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            ax(row) = ax(row) + abs(lower(k))
            ax(col) = ax(col) + abs(upper(k))
          end do
        end do
!  start dropping
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k)) .lt. delta*ax(row) .and.   &
     &          abs(upper(k)) .lt. delta*ax(col) ) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do
        deallocate(ax)

      case (2)
!  x is not present
        if (present(eps)) then
          delta = max(epsilon(1._SP),eps) / real(lnnz,SP)
        else
          delta = epsilon(1._SP) / real(lnnz,SP)
        end if
!  start dropping
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k)*upper(k)/(diag(col)*diag(row))) .lt. delta**2) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do

      case (3)
!  x is present
        if (present(eps)) then
          delta = max(epsilon(1._SP),eps) / real(lnnz,SP) * sqrt(dot_product(rhs,rhs))
        else
          delta = epsilon(1._SP) / real(lnnz,SP) * sqrt(dot_product(rhs,rhs))
        end if

!
        drops = 0
        start = 1
        do row = 2, n
          do k = start, lia(row)
            col = lja(k)
            if (abs(lower(k)*x(col)) .lt. delta .and. &
     &          abs(upper(k)*x(row)) .lt. delta) then
              drops = drops + 1
              diag(row) = diag(row) + abs(lower(k))*alpha
              diag(col) = diag(col) + abs(upper(k))*alpha

              !rhs(row) = rhs(row) - lower(k)*x(col) + abs(lower(k))*x(row)*alpha
              !rhs(col) = rhs(col) - upper(k)*x(row) + abs(upper(k))*x(col)*alpha
            else
              lower(k-drops) = lower(k)
              if (.not.symm) upper(k-drops) = upper(k)
              lja(k-drops) = col
            end if
          end do
          start = lia(row)+1
          lia(row) = lia(row) - drops
        end do

      end select

      lja => reallocate(lja,lnnz-drops)
      lower => reallocate(lower,lnnz-drops)
      if (symm) then 
        upper => lower
      else
        upper => reallocate(upper,lnnz-drops)
      end if

print '(a,f10.2,a)', 'removed', 100.*real(drops)/real(lnnz),' % matrix entries'
      return
      end 



      function a_norm_DD(x,diag,lower,upper,lia,lja)
      use femtypes
      use feminterface, only:

! Calculate Energy Norm of vector x
      implicit none
      real(DP) :: diag(:), lower(:), upper(:)
      real(DP) :: x(:)
      real(DP) :: a_norm_DD
      integer(I4B) :: lia(:), lja(:)
      intent(in) :: x, diag, lower, upper, lia, lja

      real(DP), allocatable :: ax(:)
      integer (I4B) :: n, row, col, k
      real (DP) :: delta 

      n = size(x)
      allocate(ax(n))

      ax = diag * x
      do row = 2,n
        do k = lia(row-1)+1,lia(row)
          col = lja(k)
          ax(row) = ax(row) + lower(k)*x(col)
          ax(col) = ax(col) + upper(k)*x(row)
        end do
      end do
      a_norm_DD = dot_product(x,ax)
      deallocate(ax)

      return
      end



      function a_norm_SS(x,diag,lower,upper,lia,lja)
      use femtypes
      use feminterface, only:

! Calculate Energy Norm of vector x
      implicit none
      real(SP) :: diag(:), lower(:), upper(:)
      real(SP) :: x(:)
      real(SP) :: a_norm_SS
      integer(I4B) :: lia(:), lja(:)
      intent(in) :: x, diag, lower, upper, lia, lja

      real(SP), allocatable :: ax(:)
      integer (I4B) :: n, row, col, k
      real (SP) :: delta 

      n = size(x)
      allocate(ax(n))

      ax = diag * x
      do row = 2,n
        do k = lia(row-1)+1,lia(row)
          col = lja(k)
          ax(row) = ax(row) + lower(k)*x(col)
          ax(col) = ax(col) + upper(k)*x(row)
        end do
      end do
      a_norm_SS = dot_product(x,ax)
      deallocate(ax)

      return
      end


      function f_norm_DD(diag,lower,upper,lia,lja)
      use femtypes
      use feminterface, only:
!
! Calculate Frobenius Norm of Matrix
      implicit none
      real(SP) :: diag(:), lower(:), upper(:)
      real(SP) :: f_norm_DD
      integer(I4B) :: lia(:), lja(:)
      intent(in) :: diag, lower, upper, lia, lja

      f_norm_DD = sqrt(dot_product(lower,lower) + dot_product(upper,upper) + dot_product(diag,diag))

      return
      end



      function f_norm_SS(diag,lower,upper,lia,lja)
      use femtypes
      use feminterface, only:
!
! Calculate Frobenius Norm of Matrix
      implicit none
      real(SP) :: diag(:), lower(:), upper(:)
      real(SP) :: f_norm_SS
      integer(I4B) :: lia(:), lja(:)
      intent(in) :: diag, lower, upper, lia, lja

      f_norm_SS = sqrt(dot_product(lower,lower) + dot_product(upper,upper) + dot_product(diag,diag))

      return
      end
