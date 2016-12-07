! Radiative/convective equilibrium in a single lat/lon grid cell.
module radconv_mod

  ! Convection and radiation modules
  use qe_moist_convection_mod, only : moist_convection, qsat
  use simple_sat_vapor_pres_mod, only : escomp
  use radiation_mod, only : radiation_init, radiation_down, radiation_up, radiation_end
  use constants_mod, only : cp_air, cp_vapor, grav, stefan, rho0, cp_ocean, rdgas, kappa, hlv, rvgas
  use dry_convection_mod, only: dry_convection
  
  implicit none
  private
  public :: radconveq
  
  ! User options
  real(8) :: GRIDPOWER = 1.                               ! Pressure grid spacing exponent
  integer :: NLAYER = 100                                 ! Number of layers
  integer :: NITER = 1000                                 ! Number of iterations
  real(8) :: delta_t = 1                                  ! Time step (seconds)
  real(8) :: t_init = 280.                                ! Initial temperature (isothermal)
  real(8) :: q_init = 1.e-1                               ! Initial water mixing ratio (everywhere)
  real(8) :: ocean_depth = 0.1                            ! Slab ocean depth (m)
  real(8) :: surf_k = 1.                                  ! Newton's law coefficient for surface heating
  logical :: moistconv = .true.                           ! Do moist convection?
  integer :: thin = 1                                     ! Output thinning factor
  integer :: conv_iter = 20                               ! Convective iterations per step
  real(8) :: max_dt = 5.                                  ! Maximum temperature change per step
  real(8) :: evap_rate = 0.1                              ! Dimensionless surface evaporation rate
  real(8) :: tropopause = 0.5                             ! Initial guess at tropopause (bar)
  real(8) :: psurf = 1.                                   ! Initial surface pressure (bar)
  
  ! Convection params
  real(8), allocatable, dimension(:) :: t             ! Temperature profile
  real(8), allocatable, dimension(:) :: t_ad          ! Adiabat
  real(8), allocatable, dimension(:) :: p             ! Pressure profile
  real(8), allocatable, dimension(:) :: ph            ! Pressure profiles at half indices
  real(8), allocatable, dimension(:) :: q             ! Water mixing ratio profile
  real(8), allocatable, dimension(:) :: t_out         ! Local grid cell profiles
  real(8), allocatable, dimension(:) :: p_out
  real(8), allocatable, dimension(:) :: ph_out
  real(8), allocatable, dimension(:) :: q_out
  real(8), allocatable, dimension(:) :: rain              ! Rain profile
  real(8) :: convergence                                  ! Convergence metric
  real(8), allocatable, dimension(:) :: t_prev            ! Previous timestep arrays
  real(8), allocatable, dimension(:) :: q_prev
  real(8), allocatable, dimension(:) :: p_prev
  real(8) :: raining                                      ! Raining?
  real(8) :: Ep                                           ! ?
  real(8), allocatable, dimension(:) :: t_moistad         ! Moist adiabat profile
  real(8), allocatable, dimension(:) :: qs                ! Saturated mixing ratio          
  real(8) :: dqdt                                         ! Evaporation rate
  
  ! Radiation params
  real(8) :: net_surf_sw_down
  real(8) :: surf_lw_down
  real(8), allocatable, dimension(:) :: delta_f            ! Delta flux (for temperature tendency)
  real(8), allocatable, dimension(:) :: dt                 ! Temperature change
  real(8) :: tsdt                                          ! Surface temperature tendency
  real(8) :: surf_flux                                     ! Radiation imbalance at surface
  real(8) :: sens_heat                                     ! Sensible heat flux from surface
  real(8) :: ts                                            ! Surface temperature
  real(8) :: ps                                            ! Surface pressure
  real(8), allocatable, dimension(:) :: cp                 ! Heat capacity
  
  ! General stuff
  integer :: NOUT = 23                                    ! Total number of output files
  integer :: i, j, n
  logical :: call_output                                  ! Whether or not to output
  real(8) :: time = 0                                     ! Simulation time (seconds)
  
  namelist/radconv_nml/ NITER, NLAYER, delta_t, t_init, &
                        q_init, ocean_depth, surf_k, &
                        GRIDPOWER, tropopause, &
                        moistconv, thin, conv_iter, max_dt, &
                        evap_rate, psurf
                          
  contains
  
  subroutine output()
    write(1, '(9999(e16.8))') time
    write(2, '(9999(e16.8))') convergence
    write(3, '(9999(f16.8))') ( t(i), i=1,NLAYER )
    write(4, '(9999(e16.8))') ( q(i), i=1,NLAYER )
    write(5, '(9999(e16.8))') ( p(i), i=1,NLAYER )
    write(6, '(9999(e16.8))') ( rain(i), i=1,NLAYER )
    write(7, '(9999(e16.8))') ts 
    write(8, '(9999(e16.8))') ps 
    write(9, '(9999(e16.8))') ( t_ad(i), i=1,NLAYER )
    write(10, '(9999(e16.8))') ( ph(i), i=1,NLAYER )
    write(11, '(9999(e16.8))') ( t_moistad(i), i=1,NLAYER )
    write(12, '(9999(e16.8))') ( qs(i), i=1,NLAYER )
  end subroutine
  
  subroutine output_start()
    ! radconv.f90
    open(unit=1, file='output/time.dat', ACTION="write", STATUS="replace")
    open(unit=2, file='output/convergence.dat', ACTION="write", STATUS="replace")
    open(unit=3, file='output/t.dat', ACTION="write", STATUS="replace")
    open(unit=4, file='output/q.dat', ACTION="write", STATUS="replace")
    open(unit=5, file='output/p.dat', ACTION="write", STATUS="replace")  
    open(unit=6, file='output/rain.dat', ACTION="write", STATUS="replace")
    open(unit=7, file='output/ts.dat', ACTION="write", STATUS="replace")
    open(unit=8, file='output/ps.dat', ACTION="write", STATUS="replace")
    open(unit=9, file='output/t_ad.dat', ACTION="write", STATUS="replace")
    open(unit=10, file='output/ph.dat', ACTION="write", STATUS="replace")
    open(unit=11, file='output/t_moistad.dat', ACTION="write", STATUS="replace")
    open(unit=12, file='output/qsat.dat', ACTION="write", STATUS="replace")
    ! radiation.f90
    open(unit=13, file='output/flux_sw.dat', ACTION="write", STATUS="replace")    
    open(unit=14, file='output/b.dat', ACTION="write", STATUS="replace")    
    open(unit=15, file='output/tdt_rad.dat', ACTION="write", STATUS="replace")       
    open(unit=16, file='output/tdt_sw.dat', ACTION="write", STATUS="replace")    
    open(unit=17, file='output/dtrans.dat', ACTION="write", STATUS="replace") 
    open(unit=18, file='output/lw_tau.dat', ACTION="write", STATUS="replace")
    open(unit=19, file='output/up.dat', ACTION="write", STATUS="replace") 
    open(unit=20, file='output/down.dat', ACTION="write", STATUS="replace") 
    open(unit=21, file='output/net.dat', ACTION="write", STATUS="replace") 
    open(unit=22, file='output/solar_down.dat', ACTION="write", STATUS="replace") 
    open(unit=23, file='output/flux_rad.dat', ACTION="write", STATUS="replace")
  end subroutine
  
  subroutine output_end()
    do i=1,NOUT
      close(i)
    end do
  end subroutine
  
  subroutine radconveq()

    ! Read namelist
    open(unit=98, file='radconv.nml', ACTION="read")
    read  (98, nml=radconv_nml, end=100)
    100 close(98)

    ! Allocate arrays
    allocate (t (NLAYER))
    allocate (p (NLAYER))
    allocate (ph (NLAYER + 1))
    allocate (q (NLAYER))
    allocate (t_out (NLAYER))
    allocate (p_out (NLAYER))
    allocate (ph_out (NLAYER + 1))
    allocate (q_out (NLAYER))
    allocate (rain (NLAYER))
    allocate (delta_f (NLAYER))
    allocate (dt (NLAYER))
    allocate (t_prev (NLAYER))
    allocate (p_prev (NLAYER))
    allocate (q_prev (NLAYER))
    allocate (t_moistad (NLAYER))
    allocate (qs (NLAYER))
        
    ! Output files
    call output_start()
  
    ! Set up the pressure grids
    do i=1,NLAYER+1
      ph(i) = 1 - ((NLAYER + 1 - i) / real(NLAYER + 1)) ** GRIDPOWER
    end do
    do i=1,NLAYER
      p(i) = (ph(i + 1) - ph(i)) / log(ph(i + 1) / ph(i))
    end do
    
    ! Set up even temperature and moisture grids
    t = t_init
    do i = 1, NLAYER
      if (p(i) < tropopause) then
        q(i) = 0.
      else
        q(i) = q_init
      endif
    end do
     
    ts = t(NLAYER)
    ps = p(NLAYER)
    t_ad = t(NLAYER - 1) * (p / p(NLAYER - 1)) ** kappa

    ! Set up the convergence array
    convergence = 0.
    t_prev(:) = t
    q_prev(:) = q
    p_prev(:) = p
    
    ! Calculate dry adiabat
    t_ad = t(NLAYER - 1) * (p / p(NLAYER - 1)) ** kappa

    ! Calculate saturation
    do i=1, NLAYER
      call qsat(t(i), p(i) * 1.e5, qs(i))
    end do
    
    ! TODO: Calculate moist adiabat
    t_moistad = 0
    
    ! Log to file
    call output()
  
    ! Initialize radiation module
    call radiation_init(NLAYER)
  
    ! The main loop
    do n=1,NITER
      
      ! Are we outputting this iteration?
      call_output = ((mod(n, thin) .eq. 0) .or. (n .eq. NITER)) 
      
      ! Heat capacity
      cp = cp_air * (1 - q) + cp_vapor * q
        
      ! Downward/upward radiation --> temperature change
      delta_f = 0.
      call radiation_down (ph, q, t, net_surf_sw_down, surf_lw_down)
      call radiation_up (ph, ts, t, delta_f, call_output)
      dt = delta_f * grav / cp(:) / (ph(2:NLAYER+1) - ph(1:NLAYER)) * delta_t
      
      ! Limit the temperature change per step
      do i = 1, nlayer
        if (dt(i) > max_dt) then
          dt(i) = max_dt
        else if (dt(i) < -max_dt) then
          dt(i) = -max_dt
        end if
      end do      
      
      ! Re-compute at midpoint (as in Ray's Python module)
      delta_f = 0.
      call radiation_down (ph, q, t + 0.5 * dt, net_surf_sw_down, surf_lw_down)
      call radiation_up (ph, ts, t + 0.5 * dt, delta_f, call_output)
      dt = delta_f * grav / cp(:) / (ph(2:NLAYER+1) - ph(1:NLAYER)) * delta_t
      
      ! Limit the temperature change per step
      do i = 1, nlayer
        if (dt(i) > max_dt) then
          dt(i) = max_dt
        else if (dt(i) < -max_dt) then
          dt(i) = -max_dt
        end if
      end do  
      
      ! Finally, update the temperature
      t = t + dt

      if ((moistconv) .and. (conv_iter .gt. 0)) then

        ! Call moist convection module
        call moist_convection(t(:), q(:), 1.e5 * p(:), 1.e5 * ph(:), &
                              raining, t_out, q_out, p_out, ph_out, &
                              Ep, rain, conv_iter)
      
        ! Update arrays
        t(:) = t_out
        q(:) = q_out
        p(:) = p_out / 1.e5
        ph(:) = ph_out / 1.e5

      else if (conv_iter .gt. 0) then
        
        do i = 1, conv_iter
          ! Call dry convection module
          call dry_convection(t, p, ph, t_out)
          t(:) = t_out
        end do
  
      end if

      ! --- Adjust the surface temperature ---
      ! 1. Radiation imbalance. Positive if surface is heating
      surf_flux = net_surf_sw_down + surf_lw_down - stefan * ts * ts * ts * ts
      ! 2. Sensible heat flux. Positive if surface is heating
      sens_heat = surf_k * (t(NLAYER) - ts)
      ! 3. Convert from energy flux to temperature tendency
      tsdt = (surf_flux + sens_heat) / (rho0 * cp_ocean * ocean_depth)
      ! Update the surface temperature
      ts = ts + tsdt * delta_t
      ! Adjust bottom layer of atmosphere to conserve energy
      t(NLAYER) = t(NLAYER) - sens_heat * delta_t * grav / (ph(NLAYER + 1) - ph(NLAYER)) / cp(NLAYER)
      
      ! Surface pressure
      ps = p(NLAYER)
      
      ! Add moisture to atmosphere if necessary
      call qsat(T(NLAYER), p(NLAYER) * 1.e5, qs(NLAYER))
      if (q(NLAYER) .lt. qs(NLAYER)) then
        dqdt = evap_rate * (qs(NLAYER) - q(NLAYER))
        q(NLAYER) = q(NLAYER) + dqdt * delta_t
      end if
      
      ! Calculate convergence
      convergence = sum(((t - t_prev) / t_prev) ** 2) + &
                    sum(((p - p_prev) / p_prev) ** 2) + &
                    sum(((q - q_prev) / (q_prev + 1.e-15)) ** 2)
      convergence = log10(convergence)
      
      ! Update time
      time = time + delta_t
      
      ! Log profiles to file?
      if (call_output) then
      
        ! Calculate dry adiabat
        t_ad = t(NLAYER - 1) * (p / p(NLAYER - 1)) ** kappa
    
        ! Calculate saturation
        do i=1, NLAYER
          call qsat(t(i), p(i) * 1.e5, qs(i))
        end do

        ! TODO: Calculate moist adiabat
        t_moistad = 0
  
        call output()
      end if
      
      ! Update prev arrays
      t_prev(:) = t
      q_prev(:) = q
      p_prev(:) = p

    end do

    ! Deallocate stuff
    call radiation_end
    call output_end ()
    
  end subroutine
  
end module

! This is the user-facing routine
program radconv

  use radconv_mod, only : radconveq
  call radconveq()
  
end program