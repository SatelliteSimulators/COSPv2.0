! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! WHAT ARE the rights? Copyright (c) 2020 EUMETSAT?
! All rights reserved
!
! 
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!
! History
! June 2020:    Salomon Eliasson - initital version heavily based the MODIS simulator in COSPv2.0


! Notes on using the CLARA simulator: 
!  *) Required input also includes the optical thickness and cloud top pressure 
!     derived from the ISCCP simulator run with parameter top_height = 1.
!  *) Use the same tau pres as done by the MODIS simulator
!  *) Cloud particle sizes are specified as radii, measured in meters, though within the 
!     module we use units of microns. Where particle sizes are outside the bounds used in 
!     the CLARA retrieval libraries (parameters re_water_min, re_ice_min, etc.) the 
!     simulator returns missing values (re_fill)
!
! When error conditions are encountered this code calls the function complain_and_die, 
! supplied at the bottom of this module. Users probably want to replace this with 
! something more graceful. 
!
MODULE mod_effective_radius
!  USE MOD_COSP_CONFIG, only: R_UNDEF,modis_histTau,modis_histPres,numMODISTauBins,       &
!                             numMODISPresBins,numMODISReffIceBins,numMODISReffLiqBins,   &
!                             modis_histReffIce,modis_histReffLiq
  USE MOD_COSP_CONFIG, ONLY: R_UNDEF,modis_histTau,modis_histPres,numMODISTauBins,       &
                             numMODISPresBins
  USE COSP_KINDS,      ONLY: wp
  USE MOD_COSP_STATS,  ONLY: hist2D

  implicit none
  ! ##########################################################################
  ! Retrieval parameters
  integer, parameter :: &
       num_trial_res = 15              ! Increase to make the linear pseudo-retrieval of size more accurate
   
  real(wp) :: &
       ref_TauLimit,                  & ! How deep into the model cloud does Calipso see? 
       clara_opaque,                  & ! What do we consider an opaque cloud?
       phase_TauLimit,                & ! How deep into the cloud does the phase detection see?
       size_TauLimit,                 & ! Depth of the re retreivals
       phaseDiscrimination_Threshold, & ! What fraction of total extincton needs to be in a single
                                        ! category to make phase discrim. work? 
       re_fill,                       & !
       re_water_min,                  & ! Minimum effective radius (liquid)
       re_water_max,                  & ! Maximum effective radius (liquid)
       re_ice_min,                    & ! Minimum effective radius (ice)
       re_ice_max,                    & ! Minimum effective radius (ice)
       highCloudPressureLimit,        & ! High cloud pressure limit (Pa)
       lowCloudPressureLimit            ! Low cloud pressure limit (Pa)
  integer :: &
       phaseIsLiquid,                 & !
       phaseIsIce
       
  real(wp),dimension(num_trial_res) :: &
       trial_re_w, & ! Near-IR optical params vs size for retrieval scheme (liquid)
       trial_re_i    ! Near-IR optical params vs size for retrieval scheme (ice)
  real(wp),dimension(num_trial_res) :: &
       g_w,        & ! Assymettry parameter for size retrieval (liquid)
       g_i,        & ! Assymettry parameter for size retrieval (ice)
       w0_w,       & ! Single-scattering albedo for size retrieval (liquid)
       w0_i          ! Single-scattering albedo for size retrieval (ice)
  ! Algorithmic parameters
  real(wp),parameter :: &
     ice_density = 0.93_wp ! Liquid density is 1. 
      
contains
     
  ! ########################################################################################
  function weight_by_extinction(nLevels,tauIncrement, f, tauLimit) 
    ! INPUTS
    integer, intent(in)                    :: nLevels
    real(wp),intent(in),dimension(nLevels) :: tauIncrement, f
    real(wp),intent(in)                    :: tauLimit
    ! OUTPUTS
    real(wp)                               :: weight_by_extinction
    ! LOCAL VARIABLES
    real(wp)                               :: deltaX, totalTau, totalProduct
    integer                                :: i 
    
    ! Find the extinction-weighted value of f(tau), assuming constant f within each layer
    totalTau = 0._wp; totalProduct = 0._wp
    do i = 1, size(tauIncrement)
      if(totalTau + tauIncrement(i) > tauLimit) then 
        deltaX       = tauLimit - totalTau
        totalTau     = totalTau     + deltaX
        totalProduct = totalProduct + deltaX * f(i) 
      else
        totalTau     = totalTau     + tauIncrement(i) 
        totalProduct = totalProduct + tauIncrement(i) * f(i) 
      end if 
      if(totalTau >= tauLimit) exit
    end do 

    if (totalTau > 0._wp) then
       weight_by_extinction = totalProduct/totalTau
    else
       weight_by_extinction = 0._wp
    endif
    
  end function weight_by_extinction

  ! ########################################################################################
  pure function interpolate_to_min(x, y, yobs)
    ! INPUTS
    real(wp),intent(in),dimension(num_trial_res) :: x, y 
    real(wp),intent(in)                          :: yobs
    ! OUTPUTS
    real(wp)                                     :: interpolate_to_min
    ! LOCAL VARIABLES
    real(wp), dimension(num_trial_res)           :: diff
    integer                                      :: nPoints, minDiffLoc, lowerBound, upperBound
    
    ! Given a set of values of y as y(x), find the value of x that minimizes abs(y - yobs)
    !   y must be monotonic in x
 
    nPoints = size(y)
    diff(1:num_trial_res) = y(1:num_trial_res) - yobs
    minDiffLoc = minloc(abs(diff), dim = 1) 
    
    if(minDiffLoc == 1) then 
      lowerBound = minDiffLoc
      upperBound = minDiffLoc + 1
    else if(minDiffLoc == nPoints) then
      lowerBound = minDiffLoc - 1
      upperBound = minDiffLoc
    else
      if(diff(minDiffLoc-1) * diff(minDiffLoc) < 0) then
        lowerBound = minDiffLoc-1
        upperBound = minDiffLoc
      else 
        lowerBound = minDiffLoc
        upperBound = minDiffLoc + 1
      end if 
    end if 
    
    if(diff(lowerBound) * diff(upperBound) < 0) then     
      !
      ! Interpolate the root position linearly if we bracket the root
      !
      interpolate_to_min = x(upperBound) - & 
                           diff(upperBound) * (x(upperBound) - x(lowerBound)) / (diff(upperBound) - diff(lowerBound))
    else 
      interpolate_to_min = re_fill
    end if 

  end function interpolate_to_min

  ! ########################################################################################
  ! Optical properties
  ! ########################################################################################
  elemental function get_g_nir_old (phase, re)
    ! Polynomial fit for asummetry parameter g in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick

    ! INPUTS
    integer, intent(in) :: phase
    real(wp),intent(in) :: re
    ! OUTPUTS
    real(wp)            :: get_g_nir_old 
    ! LOCAL VARIABLES(parameters)
    real(wp), dimension(3), parameter :: &
         ice_coefficients         = (/ 0.7432,  4.5563e-3, -2.8697e-5 /), & 
         small_water_coefficients = (/ 0.8027, -1.0496e-2,  1.7071e-3 /), & 
         big_water_coefficients   = (/ 0.7931,  5.3087e-3, -7.4995e-5 /) 
   
    ! approx. fits from MODIS Collection 5 LUT scattering calculations
    if(phase == phaseIsLiquid) then
      if(re < 8.) then 
        get_g_nir_old = fit_to_quadratic(re, small_water_coefficients)
        if(re < re_water_min) get_g_nir_old = fit_to_quadratic(re_water_min, small_water_coefficients)
      else
        get_g_nir_old = fit_to_quadratic(re,   big_water_coefficients)
        if(re > re_water_max) get_g_nir_old = fit_to_quadratic(re_water_max, big_water_coefficients)
      end if 
    else
      get_g_nir_old = fit_to_quadratic(re, ice_coefficients)
      if(re < re_ice_min) get_g_nir_old = fit_to_quadratic(re_ice_min, ice_coefficients)
      if(re > re_ice_max) get_g_nir_old = fit_to_quadratic(re_ice_max, ice_coefficients)
    end if 
    
  end function get_g_nir_old

  ! ########################################################################################
  elemental function get_ssa_nir_old (phase, re)
    ! Polynomial fit for single scattering albedo in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick
    
    ! INPUTS
    integer, intent(in) :: phase
    real(wp),intent(in) :: re
    ! OUTPUTS
    real(wp)            :: get_ssa_nir_old
    ! LOCAL VARIABLES (parameters)
    real(wp), dimension(4), parameter :: ice_coefficients   = (/ 0.9994, -4.5199e-3, 3.9370e-5, -1.5235e-7 /)
    real(wp), dimension(3), parameter :: water_coefficients = (/ 1.0008, -2.5626e-3, 1.6024e-5 /) 
    
    ! approx. fits from MODIS Collection 5 LUT scattering calculations
    if(phase == phaseIsLiquid) then
       get_ssa_nir_old = fit_to_quadratic(re, water_coefficients)
       if(re < re_water_min) get_ssa_nir_old = fit_to_quadratic(re_water_min, water_coefficients)
       if(re > re_water_max) get_ssa_nir_old = fit_to_quadratic(re_water_max, water_coefficients)
    else
       get_ssa_nir_old = fit_to_cubic(re, ice_coefficients)
       if(re < re_ice_min) get_ssa_nir_old = fit_to_cubic(re_ice_min, ice_coefficients)
       if(re > re_ice_max) get_ssa_nir_old = fit_to_cubic(re_ice_max, ice_coefficients)
    end if
    
  end function get_ssa_nir_old
  
  elemental function get_g_nir (phase, re)
    !
    ! Polynomial fit for asummetry parameter g in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick
    !

    integer, intent(in) :: phase
    real(wp),    intent(in) :: re
    real(wp) :: get_g_nir 

    real(wp), dimension(3), parameter :: ice_coefficients         = (/ 0.7490, 6.5153e-3, -5.4136e-5 /), &
                                         small_water_coefficients = (/ 1.0364, -8.8800e-2, 7.0000e-3 /)
    real(wp), dimension(4), parameter :: big_water_coefficients   = (/ 0.6035, 2.8993e-2, -1.1051e-3, 1.5134e-5 /)

    ! approx. fits from MODIS Collection 6 LUT scattering calculations for 3.7 µm channel size retrievals
    if(phase == phaseIsLiquid) then 
       if(re < 7.) then
          get_g_nir = fit_to_quadratic(re, small_water_coefficients)
          if(re < re_water_min) get_g_nir = fit_to_quadratic(re_water_min, small_water_coefficients)
       else
          get_g_nir = fit_to_cubic(re, big_water_coefficients)
          if(re > re_water_max) get_g_nir = fit_to_cubic(re_water_max, big_water_coefficients)
       end if
    else
       get_g_nir = fit_to_quadratic(re, ice_coefficients)
      if(re < re_ice_min) get_g_nir = fit_to_quadratic(re_ice_min, ice_coefficients)
      if(re > re_ice_max) get_g_nir = fit_to_quadratic(re_ice_max, ice_coefficients)
    end if 
    
  end function get_g_nir

  ! --------------------------------------------
    elemental function get_ssa_nir (phase, re)
        integer, intent(in) :: phase
        real(wp),    intent(in) :: re
        real(wp)                :: get_ssa_nir
        !
        ! Polynomial fit for single scattering albedo in MODIS band 7 (near IR) as a function 
        !   of size for ice and water
        ! Fits from Steve Platnick
        !
        real(wp), dimension(4), parameter :: ice_coefficients   = (/ 0.9625, -1.8069e-2, 3.3281e-4,-2.2865e-6/)
        real(wp), dimension(3), parameter :: water_coefficients = (/ 1.0044, -1.1397e-2, 1.3300e-4 /)
        
        ! approx. fits from MODIS Collection 6 LUT scattering calculations
        if(phase == phaseIsLiquid) then
          get_ssa_nir = fit_to_quadratic(re, water_coefficients)
          if(re < re_water_min) get_ssa_nir = fit_to_quadratic(re_water_min, water_coefficients)
          if(re > re_water_max) get_ssa_nir = fit_to_quadratic(re_water_max, water_coefficients)
        else
          get_ssa_nir = fit_to_cubic(re, ice_coefficients)
          if(re < re_ice_min) get_ssa_nir = fit_to_cubic(re_ice_min, ice_coefficients)
          if(re > re_ice_max) get_ssa_nir = fit_to_cubic(re_ice_max, ice_coefficients)
        end if 

    end function get_ssa_nir

  

  ! ########################################################################################
  pure function fit_to_cubic(x, coefficients) 
    ! INPUTS
    real(wp),               intent(in) :: x
    real(wp), dimension(4), intent(in) :: coefficients
    ! OUTPUTS
    real(wp)                           :: fit_to_cubic  
    
    fit_to_cubic = coefficients(1) + x * (coefficients(2) + x * (coefficients(3) + x * coefficients(4)))
  end function fit_to_cubic
    
  ! ########################################################################################
  pure function fit_to_quadratic(x, coefficients) 
    ! INPUTS
    real(wp),               intent(in) :: x
    real(wp), dimension(3), intent(in) :: coefficients
    ! OUTPUTS
    real(wp)                           :: fit_to_quadratic
    
    fit_to_quadratic = coefficients(1) + x * (coefficients(2) + x * (coefficients(3)))
  end function fit_to_quadratic

  ! ########################################################################################
  ! Radiative transfer
  ! ########################################################################################
  pure function compute_toa_reflectace(nLevels,tau, g, w0)
    ! This wrapper reports reflectance only and strips out non-cloudy elements from the 
    ! calculation
    
    ! INPUTS
    integer,intent(in)                     :: nLevels
    real(wp),intent(in),dimension(nLevels) :: tau, g, w0
    ! OUTPUTS
    real(wp)                               :: compute_toa_reflectace
    ! LOCAL VARIABLES
    logical, dimension(nLevels)                   :: cloudMask
    integer, dimension(count(tau(1:nLevels) > 0)) :: cloudIndicies
    real(wp),dimension(count(tau(1:nLevels) > 0)) :: Refl,Trans
    real(wp)                                      :: Refl_tot, Trans_tot
    integer                                       :: i

    cloudMask(1:nLevels) = tau(1:nLevels) > 0. 
    cloudIndicies = pack((/ (i, i = 1, nLevels) /), mask = cloudMask) 
    do i = 1, size(cloudIndicies)
       call two_stream(tau(cloudIndicies(i)), g(cloudIndicies(i)), w0(cloudIndicies(i)), Refl(i), Trans(i))
    end do
    
    call adding_doubling(count(tau(1:nLevels) > 0),Refl(:), Trans(:), Refl_tot, Trans_tot)  
    
    compute_toa_reflectace = Refl_tot
    
  end function compute_toa_reflectace
 
  ! ########################################################################################
  pure subroutine two_stream(tauint, gint, w0int, ref, tra) 
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    ! INPUTS
    real(wp), intent(in)  :: tauint, gint, w0int
    ! OUTPUTS
    real(wp), intent(out) :: ref, tra
    ! LOCAL VARIABLES
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    integer, parameter :: beam = 2
    real(wp),parameter :: xmu = 0.866, minConservativeW0 = 0.9999999
    real(wp) :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
         rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den, th
    
    ! Compute reflectance and transmittance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    f   = gint**2
    tau = (1._wp - w0int * f) * tauint
    w0  = (1._wp - f) * w0int / (1._wp - w0int * f)
    g   = (gint - f) / (1._wp - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7._wp - w0* (4._wp + 3._wp * g)) / 4._wp
    gamma2 = -(1._wp - w0* (4._wp - 3._wp * g)) / 4._wp
    gamma3 =  (2._wp - 3._wp*g*xmu) / 4._wp
    gamma4 =   1._wp - gamma3

    if (w0int > minConservativeW0) then
      ! Conservative scattering
      if (beam == 1) then
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-exp(-tau/xmu)))
  
          ref = rh / (1._wp + gamma1 * tau)
          tra = 1._wp - ref       
      else if(beam == 2) then
          ref = gamma1*tau/(1._wp + gamma1*tau)
          tra = 1._wp - ref
      endif
    else
      ! Non-conservative scattering
      a1 = gamma1 * gamma4 + gamma2 * gamma3
      a2 = gamma1 * gamma3 + gamma2 * gamma4

      rk = sqrt(gamma1**2 - gamma2**2)
      
      r1 = (1._wp - rk * xmu) * (a2 + rk * gamma3)
      r2 = (1._wp + rk * xmu) * (a2 - rk * gamma3)
      r3 = 2._wp * rk *(gamma3 - a2 * xmu)
      r4 = (1._wp - (rk * xmu)**2) * (rk + gamma1)
      r5 = (1._wp - (rk * xmu)**2) * (rk - gamma1)
      
      t1 = (1._wp + rk * xmu) * (a1 + rk * gamma4)
      t2 = (1._wp - rk * xmu) * (a1 - rk * gamma4)
      t3 = 2._wp * rk * (gamma4 + a1 * xmu)
      t4 = r4
      t5 = r5

      beta = -r5 / r4         
  
      e1 = min(rk * tau, 500._wp) 
      e2 = min(tau / xmu, 500._wp) 
      
      if (beam == 1) then
         den = r4 * exp(e1) + r5 * exp(-e1)
         ref  = w0*(r1*exp(e1)-r2*exp(-e1)-r3*exp(-e2))/den
         den = t4 * exp(e1) + t5 * exp(-e1)
         th  = exp(-e2)
         tra = th-th*w0*(t1*exp(e1)-t2*exp(-e1)-t3*exp(e2))/den
      elseif (beam == 2) then
         ef1 = exp(-e1)
         ef2 = exp(-2*e1)
         ref = (gamma2*(1._wp-ef2))/((rk+gamma1)*(1._wp-beta*ef2))
         tra = (2._wp*rk*ef1)/((rk+gamma1)*(1._wp-beta*ef2))
      endif
    end if
  end subroutine two_stream

  ! ########################################################################################
  elemental function two_stream_reflectance(tauint, gint, w0int)
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    
    ! INPUTS
    real(wp), intent(in) :: tauint, gint, w0int
    ! OUTPUTS
    real(wp)             :: two_stream_reflectance
    ! LOCAL VARIABLES
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    integer, parameter :: beam = 2
    real(wp),parameter :: xmu = 0.866, minConservativeW0 = 0.9999999
    real(wp) :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
         rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den

    f   = gint**2
    tau = (1._wp - w0int * f) * tauint
    w0  = (1._wp - f) * w0int / (1._wp - w0int * f)
    g   = (gint - f) / (1._wp - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7._wp - w0* (4._wp + 3._wp * g)) / 4._wp
    gamma2 = -(1._wp - w0* (4._wp - 3._wp * g)) / 4._wp
    gamma3 =  (2._wp - 3._wp*g*xmu) / 4._wp
    gamma4 =   1._wp - gamma3

    if (w0int > minConservativeW0) then
      ! Conservative scattering
      if (beam == 1) then
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-exp(-tau/xmu)))
          two_stream_reflectance = rh / (1._wp + gamma1 * tau)
      elseif (beam == 2) then
          two_stream_reflectance = gamma1*tau/(1._wp + gamma1*tau)
      endif
        
    else    !

        ! Non-conservative scattering
         a1 = gamma1 * gamma4 + gamma2 * gamma3
         a2 = gamma1 * gamma3 + gamma2 * gamma4

         rk = sqrt(gamma1**2 - gamma2**2)
         
         r1 = (1._wp - rk * xmu) * (a2 + rk * gamma3)
         r2 = (1._wp + rk * xmu) * (a2 - rk * gamma3)
         r3 = 2._wp * rk *(gamma3 - a2 * xmu)
         r4 = (1._wp - (rk * xmu)**2) * (rk + gamma1)
         r5 = (1._wp - (rk * xmu)**2) * (rk - gamma1)
         
         t1 = (1._wp + rk * xmu) * (a1 + rk * gamma4)
         t2 = (1._wp - rk * xmu) * (a1 - rk * gamma4)
         t3 = 2._wp * rk * (gamma4 + a1 * xmu)
         t4 = r4
         t5 = r5

         beta = -r5 / r4         
         
         e1 = min(rk * tau, 500._wp) 
         e2 = min(tau / xmu, 500._wp) 
         
         if (beam == 1) then
           den = r4 * exp(e1) + r5 * exp(-e1)
           two_stream_reflectance  = w0*(r1*exp(e1)-r2*exp(-e1)-r3*exp(-e2))/den
         elseif (beam == 2) then
           ef1 = exp(-e1)
           ef2 = exp(-2*e1)
           two_stream_reflectance = (gamma2*(1._wp-ef2))/((rk+gamma1)*(1._wp-beta*ef2))
         endif
           
      end if
  end function two_stream_reflectance 

  ! ########################################################################################
  pure subroutine adding_doubling (npts,Refl, Tran, Refl_tot, Tran_tot)      
    ! Use adding/doubling formulas to compute total reflectance and transmittance from 
    ! layer values
    
    ! INPUTS
    integer,intent(in)                  :: npts
    real(wp),intent(in),dimension(npts) :: Refl,Tran
    ! OUTPUTS
    real(wp),intent(out)                :: Refl_tot, Tran_tot
    ! LOCAL VARIABLES
    integer :: i
    real(wp), dimension(npts) :: Refl_cumulative, Tran_cumulative
    
    Refl_cumulative(1) = Refl(1)
    Tran_cumulative(1) = Tran(1)    
    
    do i=2, npts
       ! place (add) previous combined layer(s) reflectance on top of layer i, w/black surface (or ignoring surface):
       Refl_cumulative(i) = Refl_cumulative(i-1) + Refl(i)*(Tran_cumulative(i-1)**2)/(1._wp - Refl_cumulative(i-1) * Refl(i))
       Tran_cumulative(i) = (Tran_cumulative(i-1)*Tran(i)) / (1._wp - Refl_cumulative(i-1) * Refl(i))
    end do
    
    Refl_tot = Refl_cumulative(size(Refl))
    Tran_tot = Tran_cumulative(size(Refl))
    
  end subroutine adding_doubling

end module mod_effective_radius
