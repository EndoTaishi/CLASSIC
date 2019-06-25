!> \file
!! Utility to determine location of closest value in an array
!      * June 18/2019 - M. Fortier
!      *   - Significantly refactored the code, removing the GOTO statement
!      * MAY 12/2004 - L.SOLHEIM
!      *
!>     * GIVEN A MONOTONIC VECTOR V OF LENGTH N AND A VALUE X,
!!     * RETURN THE INDEX MVIDX SUCH THAT X IS BETWEEN
!!     * V(MVIDX) AND V(MVIDX+1).
!!     *
!!     * V MUST BE MONOTONIC, EITHER INCREASING OF DECREASING.
!!     * THERE IS NO CHECK ON WHETHER OR NOT THIS VECTOR IS
!!     * MONOTONIC.
!!     *
!!     * THIS FUNCTION RETURNS 1 OR N-1 IF X IS OUT OF RANGE.
!!     *
!!     * INPUT:
!!     *
!!     *   REAL    V(N) ...MONITONIC VECTOR (INCREASING OR DECREASING)
!!     *   INTEGER N    ...SIZE OF V
!!     *   REAL    X    ...SINGLE REAL VALUE
!!     *
!!     * OUTPUT:
!!     *
!!     *   V(MVIDX) .LE/>= X <=/>= V(MVIDX+1)
!!
integer function MVIDX(V,N,X)
  !-----------------------------------------------------------------------
  implicit none
  
  real :: X,V(N)
  integer :: N
  integer :: JL,JM,JU
  !-----------------------------------------------------------------------
  
  if (X == V(1)) then
    MVIDX = 1
    return
  else if (X == V(N)) then
    MVIDX = N - 1
    return
  end if
  JL = 1
  JU = N
  do while (JU - JL > 1)
    JM = (JU + JL) / 2
    if ((V(N) > V(1)).eqv.(X > V(JM))) then
      JL = JM
    else
      JU = JM
    end if
  end do
  MVIDX = JL
  return
end
