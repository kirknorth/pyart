!  Module: minimum_detection.f90


subroutine mds(power, fill_value, nr, ng, P, Q, R2, N)

   implicit none

   integer(kind=4), intent(in)                :: nr, ng
   real(kind=8), intent(in)                   :: fill_value
   real(kind=8), intent(in), dimension(nr,ng) :: power
   real(kind=8), intent(out), dimension(ng)   :: P, Q, R2, N

!  Define local variables ====================================================

   integer(kind=4)           :: g, r
   logical, dimension(nr,ng) :: power_mask
   real(kind=8)              :: P_tmp, Q_tmp, R2_tmp, N_tmp

!  ===========================================================================

!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nr, ng
   !f2py real(kind=8), intent(in)              :: fill_value
   !f2py real(kind=8), intent(in)              :: power
   !f2py real(kind=8), intent(out)             :: P, Q, R2, N

!  ===========================================================================

   P = fill_value
   Q = fill_value
   R2 = fill_value
   N = fill_value

   power_mask = power /= fill_value

   do g = 1, ng
      do r = 1, nr
!       Get the number of (possible) noisey gates (N) and then compute the
!       mean (P) and variance (Q) of the signal for the current iteration
        N_tmp = dble(count(power_mask(r:nr,g)))
        P_tmp = sum(power(r:nr,g), power_mask(r:nr,g)) / N_tmp
        Q_tmp = sum(power(r:nr,g)**2, power_mask(r:nr,g)) / N_tmp - P_tmp**2
        R2_tmp = P_tmp**2 / Q_tmp

!       If the white-noise criterion is met, then we have found our noise
!       floor for the current range and we can save the appropriate variables
!       and exit the current iteration
        if (R2_tmp > 1.d0) then
           P(g) = P_tmp
           Q(g) = Q_tmp
           R2(g) = power(r,g)
           N(g) = N_tmp
           exit
        endif

      enddo
   enddo

   return

end subroutine mds
