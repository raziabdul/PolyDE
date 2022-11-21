subroutine orderk ( k, ind )
!
!******************************************************************************
!
!! ORDERK reorders K elements of an array in nondecreasing order.
!
!
!  Purpose: 
!
!    Order K elements of array IND in nondecreasing order.
!    It is assume that K is small, say <= 15, so that insertion sort
!    is used. If K is larger, a faster sort such as heapsort should
!    be used.
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
!    Input, K - size of array IND.
!
!    Input/output, IND(1:K), an array, which is sorted on output.
!
  implicit none
!
  integer k
!
  integer i
  integer ind(k)
  integer j
  integer s
  integer t
!
  do i = 2, k

    t = ind(i)
    j = i

10  continue

     s = ind(j-1)

     if (t < s) then
       ind(j) = s
       j = j - 1
       if (j > 1) go to 10
     end if

     ind(j) = t

   end do

  return
end
