! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
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
! History
! May 2015 - D. Swales - Original version
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_CLARA_INTERFACE
  USE COSP_KINDS,      ONLY: wp
  USE MOD_COSP_CONFIG, ONLY: R_UNDEF
  USE MOD_CLARA_SIM,   ONLY: num_trial_res,min_OpticalThickness,clara_opaque,&
                             refTauLimit,phase_TauLimit,size_TauLimit,re_fill,   &
                             phaseDiscrimination_Threshold,re_water_min,     &
                             re_water_max,re_ice_min,re_ice_max,               &
                             highCloudPressureLimit,lowCloudPressureLimit,phaseIsNone,   &
                             phaseIsLiquid,phaseIsIce,phaseIsUndetermined,trial_re_w,    &
                             trial_re_i,g_w,g_i,w0_w,w0_i, get_g_nir,get_ssa_nir
  IMPLICIT NONE

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  TYPE clara_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  TYPE CLARA_IN
     INTEGER,POINTER :: &
          Npoints,        & ! Number of horizontal gridpoints
          Ncolumns,       & ! Number of subcolumns
          Nlevels           ! Number of vertical levels
     INTEGER :: &
          Nsunlit           ! Number of sunlit lit pixels
     REAL(wp),ALLOCATABLE,DIMENSION(:) :: &
          sunlit,         & ! Sunlit scenes
          notSunlit         ! Dark scenes
     REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: &
          pres              ! Gridmean pressure at layer edges (Pa) 
     REAL(wp),POINTER ::  &
          tau(:,:,:),     & ! Subcolumn optical thickness @ 0.67 microns.
          liqFrac(:,:,:), & ! Liquid water fraction
          g(:,:,:),       & ! Subcolumn assymetry parameter  
          w0(:,:,:)         ! Subcolumn single-scattering albedo
  END TYPE CLARA_IN
CONTAINS
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROTUINE cosp_clara_init
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CLARA_INIT()
    integer :: i 
    
    ! Retrieval parameters
    min_OpticalThickness          = 0.225_wp   ! Minimum detectable optical thickness 
                                               ! i.e., global limit for day and night
                                               ! combined. This is in place until the
                                               ! variable limit based on atmospheric and
                                               ! surface conditions are in place (2906-2021)
    ref_TauLimit                  = 1._wp      ! How deep into the cloud does CO2 slicing 
                                               ! see? 
    phase_TauLimit                = 1._wp      ! How deep into the cloud does the phase 
                                               ! detection see?
    size_TauLimit                 = 2._wp      ! Depth of the re retreivals
    phaseDiscrimination_Threshold = 0.5_wp     ! What fraction of total extincton needs to 
                                               ! be in a single category to make phase 
                                               ! discrim. CLARA has no  'undefined'-category 
                                               ! (therefore 0.5 instead of 0.7 for MODIS)
    re_fill                       = -999._wp   ! Fill value
    re_water_min                  = 4._wp      ! Minimum effective radius (liquid)
    re_water_max                  = 30._wp     ! Maximum effective radius (liquid)
    re_ice_min                    = 5._wp      ! Minimum effective radius (ice)
    re_ice_max                    = 90._wp     ! Minimum effective radius (ice)
    highCloudPressureLimit        = 44000._wp  ! High cloud pressure limit (Pa)
    lowCloudPressureLimit         = 68000._wp  ! Low cloud pressure limit (Pa)
    phaseIsNone                   = 0          ! No cloud
    phaseIsLiquid                 = 1          ! Liquid cloud
    phaseIsIce                    = 2          ! Ice cloud

    ! Precompute near-IR optical params vs size for retrieval scheme    
    trial_re_w(1:num_trial_res) = re_water_min + (re_water_max - re_water_min) /         &
         (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)
    trial_re_i(1:num_trial_res) = re_ice_min   + (re_ice_max -   re_ice_min) /           &
         (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)

    ! Initialize estimates for size retrieval
    g_w(1:num_trial_res)  = get_g_nir(phaseIsLiquid,trial_re_w(1:num_trial_res))
    w0_w(1:num_trial_res) = get_ssa_nir(phaseIsLiquid,trial_re_w(1:num_trial_res))
    g_i(1:num_trial_res)  = get_g_nir(phaseIsIce,trial_re_i(1:num_trial_res))
    w0_i(1:num_trial_res) = get_ssa_nir(phaseIsIce,trial_re_i(1:num_trial_res))

  END SUBROUTINE COSP_CLARA_INIT
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE MOD_COSP_Modis_INTERFACE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_Modis_INTERFACE
