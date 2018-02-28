module radiation_mod

! ==================================================================================
! ==================================================================================

use constants_mod, only: stefan, cp_air, grav, pstd_mks
 
!==================================================================================
implicit none
private
!==================================================================================

!==================================================================================
public :: radiation_init, radiation_down, radiation_up, radiation_end
!==================================================================================

! module variables
real(8)    :: solar_constant  = 1360.0
real(8)    :: del_sol         = 1.4
real(8)    :: del_sw          = 0.0
real(8)    :: atm_abs         = 0.0
real(8)    :: sw_diff         = 0.0
real(8)    :: linear_tau      = 0.1     ! This is kappa * psurf / g for air (tau_1)
real(8)    :: albedo          = 0.3
real(8)    :: window          = 0.0     ! spectral window transparent to LW
real(8)    :: wv_exponent     = 4.0 
real(8)    :: solar_exponent  = 4.0 
real(8)    :: wv_tau          = 1.0     ! This is kappa for water vapor (opacity)
real(8)    :: reference_slp   = 1.e5    ! Reference pressure for solar absorption; I think this is in Pascals and needs to change to bars
logical    :: p_broad         = .false. ! Pressure broadening?
integer(8) :: NLAYER
integer(8) :: n

real(8), allocatable, dimension(:)  :: b, tdt_rad, entrop_rad, tdt_sw
real(8), allocatable, dimension(:)  :: up, down, net, solar_down, flux_rad, flux_sw
real(8), allocatable, dimension(:)  :: dtrans, lw_tau
real(8), allocatable, dimension(:)  :: dp_half
real(8), allocatable, dimension(:)  :: solar_tau
real(8) :: solar, lw_tau_0
real(8) :: b_surf
real(8) :: olr, swin
real(8) :: solar_tau_0, ss, p2
real(8), save                       :: pi, deg_to_rad, rad_to_deg

namelist/radiation_nml/ solar_constant, del_sol, &
           atm_abs, sw_diff, &
           linear_tau, wv_tau, del_sw, albedo, &
           window, wv_exponent, solar_exponent, reference_slp, p_broad

contains

subroutine output() 
  write(13, '(9999(f16.8))') ( flux_sw(n), n=1,NLAYER )
  write(14, '(9999(f16.8))') ( b(n), n=1,NLAYER )
  write(15, '(9999(f16.8))') ( tdt_rad(n), n=1,NLAYER )
  write(16, '(9999(f16.8))') ( tdt_sw(n), n=1,NLAYER )
  write(17, '(9999(f16.8))') ( dtrans(n), n=1,NLAYER )
  write(18, '(9999(f16.8))') ( lw_tau(n), n=1,NLAYER )
  write(19, '(9999(f16.8))') ( up(n), n=1,NLAYER )
  write(20, '(9999(f16.8))') ( down(n), n=1,NLAYER )
  write(21, '(9999(f16.8))') ( net(n), n=1,NLAYER )
  write(22, '(9999(f16.8))') ( solar_down(n), n=1,NLAYER )
  write(23, '(9999(f16.8))') ( flux_rad(n), n=1,NLAYER )
end subroutine

! ==================================================================================

subroutine radiation_init(num_levels)

  integer, intent(in)                     :: num_levels
  integer, dimension(3)                   :: half = (/1,2,4/)
  integer                                 :: ierr, io, unit, i

  ! read namelist
  open(unit=99, file='radiation.nml', ACTION="read")
  read  (99, nml=radiation_nml, end=100)
  100 close(99)
  
  NLAYER = num_levels
  pi    = 4.0*atan(1.)
  deg_to_rad = 2.*pi/360.
  rad_to_deg = 360.0/2./pi

  allocate (b                (num_levels))
  allocate (tdt_rad          (num_levels))
  allocate (tdt_sw           (num_levels))
  allocate (entrop_rad       (num_levels))
  allocate (dtrans           (num_levels))
  allocate (dp_half          (num_levels))
  allocate (up               (num_levels+1))
  allocate (down             (num_levels+1))
  allocate (net              (num_levels+1))
  allocate (solar_down       (num_levels+1))
  allocate (flux_rad         (num_levels+1))
  allocate (flux_sw          (num_levels+1))
  allocate (lw_tau           (num_levels+1))
  allocate (solar_tau        (num_levels+1))

  ! Output
  dtrans(:) = 0
  lw_tau(:) = 0
  up(:) = 0
  down(:) = 0
  net(:) = 0
  solar_down(:) = 0
  flux_rad(:) = 0
  flux_sw(:) = 0
  b(:) = 0
  tdt_rad(:) = 0
  tdt_sw(:) = 0
  call output()
  
  return
  
end subroutine radiation_init

! ==================================================================================

subroutine radiation_down (p_half, q, t, net_surf_sw_down, surf_lw_down)

  ! Begin the radiation calculation by computing downward fluxes.
  ! This part of the calculation does not depend on the surface temperature.
  real(8), intent(out) :: net_surf_sw_down
  real(8), intent(out) :: surf_lw_down
  real(8), intent(in) , dimension(:)  :: t, p_half, q
  integer                             :: i, j, k, n

  n = size(t,1)
  
  ! We are at the equator, so sin(lat) = 0
  ss = 0
  p2 = (1. - 3. * ss * ss) / 4.
  solar = 0.25 * solar_constant * (1.0 + del_sol * p2 + del_sw * ss)
  solar_tau_0 = (1.0 - sw_diff * ss * ss) * atm_abs

  ! Derivative of the pressure grid
  do k=1, n
     dp_half(k) = p_half(k + 1) - p_half(k)
  end do
  
  ! Longwave optical depth
  lw_tau = 0.
  do k = 1, n
     if (p_broad) then
        ! Air broadening, reference pressure = 1 bar
        lw_tau(k+1) = lw_tau(k) + q(k) * dp_half(k) * (p_half(k) / 1.e5)
     else
        lw_tau(k+1) = lw_tau(k) + q(k) * dp_half(k)
     end if
     solar_tau(k) = solar_tau_0 * (p_half(k) / reference_slp) ** solar_exponent
  end do
  ! NOTE: Removed arbitrary division by 10 in line below
  lw_tau = lw_tau * wv_tau / grav

  ! Add pressure broadening due to N2 here, eventually (quadratic tau)
  lw_tau(:) = lw_tau(:) + linear_tau * p_half(:) / p_half(n + 1)
  
  ! No radiation from spectral window
  b = (1.0 - window) * stefan * t * t * t * t

  do k = 1, n
    dtrans(k) = exp(-(lw_tau(k + 1) - lw_tau(k)))
  end do

  down(1) = 0.0
  do k = 1,n
    down(k+1) = down(k) * dtrans(k) + b(k) * (1.0 - dtrans(k))
  end do

  do k = 1,n+1
    solar_down(k) = solar * exp(-solar_tau(k))
  end do

  surf_lw_down     = down(n + 1)
  net_surf_sw_down = solar_down(n + 1) * (1. - albedo)
  swin = solar_down(1)
  
  return
end subroutine radiation_down

! ==================================================================================

subroutine radiation_up (p_half, t_surf, t, tdt, call_output)

  ! Now complete the radiation calculation by computing the upward and net fluxes.

  real(8), intent(in)                       :: t_surf
  real(8), intent(in) , dimension(:)        :: t, p_half
  real(8), intent(inout), dimension(:)      :: tdt
  logical, intent(in)                       :: call_output
  integer                                   :: j, k, n

  n = size(t,1)

  ! total flux from surface
  b_surf = stefan * t_surf * t_surf * t_surf * t_surf

  ! first deal with non-window upward flux
  up(n+1) = b_surf * (1.0 - window)
  do k = n,1,-1
    up(k) = up(k + 1) * dtrans(k) + b(k) * (1.0 - dtrans(k))
  end do

  ! add upward flux in spectral window
  do k = 1,n+1
   up(k) = up(k) + b_surf * window
  end do

  do k = 1,n+1
    net(k) = up(k) - down(k)
    flux_sw(k) = albedo * solar_down(n + 1) - solar_down(k)
    flux_rad(k) = net(k) + flux_sw(k)
  end do

  do k = 1,n
    tdt_rad(k) = (net(k+1) - net(k) - solar_down(k + 1) + solar_down(k))  
    tdt_sw(k) = (-solar_down(k + 1) + solar_down(k)) * grav / (cp_air * (p_half(k + 1) - p_half(k)))
    tdt(k) = tdt(k) + tdt_rad(k)
  end do
  
  ! Outgoing longwave
  olr = up(1)
  
  ! Output
  if (call_output) then
    call output()
  end if
  
  return
  
end subroutine radiation_up

! ==================================================================================
                                                                  
subroutine radiation_end
                                                                                                      
  deallocate (b, tdt_rad, tdt_sw, entrop_rad) 
  deallocate (up, down, net, solar_down, flux_rad, flux_sw)
  deallocate (dtrans)   
  deallocate (lw_tau, solar_tau)
  
end subroutine radiation_end

end module radiation_mod