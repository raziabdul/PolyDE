subroutine order3 ( i, j, k )
!
!******************************************************************************
!
!! ORDER3 reorders 3 integers into ascending order.
!
!
!  Purpose: 
!
!    Order I, J, K so that I <= J <= K.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input/output, integer I, J, K, on output are sorted into
!    nondecreasing order.
!
  implicit none
!
  integer i
  integer j
  integer k
  integer t
!
  if (j < i) then
    if (k < j) then
      call i_swap ( i, k )
    else if (k < i) then
      t = i
      i = j
      j = k
      k = t
    else
      call i_swap ( i, j )
    end if
  else
    if (k < i) then
      t = i
      i = k
      k = j
      j = t
    else if (k < j) then
      call i_swap ( j, k )
    end if
  end if

  return
end
