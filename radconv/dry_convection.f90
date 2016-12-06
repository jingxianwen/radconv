module dry_convection_mod

use constants_mod, only: grav, rdgas, rvgas, kappa, latent_heat, &
        cp_air, cp_vapor, cp_water

implicit none
private

public :: dry_convection

real(8), public, parameter :: d622 = rdgas/rvgas
real(8), public, parameter :: d378 = 1.-d622
real(8), public, parameter :: d608 = d378/d622

contains

subroutine dry_convection(tin, p_full, p_half, Tref)
                
real(8), intent(in), dimension(:) :: tin, p_full
real(8), intent(in), dimension(:) :: p_half
real(8), intent(out), dimension(:) :: Tref
real(8) :: T1, P1, P2, T2, DTref, deltap1, deltap2
integer :: n, k


   n = size(tin)
   Tref = tin *1.
   
   do k=n, 2, -1
          T1 = Tref(k)
          P1 = p_full(k)
          P2 = p_full(k-1)

          T2 = T1 *  (P2 / P1)**kappa
          if (Tref(k-1) .lt. T2) then
             deltap1 = p_half(k+1) - p_half(k)
             deltap2 = p_half(k) - p_half(k-1)
             DTref = T1 - T2

             T2 = (Tref(k) * deltap1 + Tref(k-1) * deltap2 &
                     - DTref*deltap2) / (deltap1 + deltap2)
             Tref(k-1) = T2
             Tref(k) = T2 + DTref
          end if
   end do

   do k=2,n 
          T1 = Tref(k)
          P1 = p_full(k)
          P2 = p_full(k-1)
          T2 = T1 *  (P2 / P1)**kappa
          if (Tref(k-1) .lt. T2) then
             deltap1 = p_half(k+1) - p_half(k)
             deltap2 = p_half(k) - p_half(k-1)
             DTref = T1 - T2
             T2 = (Tref(k) * deltap1 + Tref(k-1) * deltap2 &
                     - DTref*deltap2) / (deltap1 + deltap2)
             Tref(k-1) = T2
             Tref(k) = T2 + DTref
          end if
   end do
end subroutine dry_convection

end module dry_convection_mod