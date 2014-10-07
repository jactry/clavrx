! $Id: pixel_routines.f90 538 2014-09-15 20:30:48Z dbotambekov $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: cloud_cover_layers_module.f90 (src)
!       CLOUD_COVER_LAYERS (module)
!
! PURPOSE: this module houses routines for computing cloud cover or
!          fraction in each atmospheric layer including the total
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! Public routines used in this MODULE:
!
!--------------------------------------------------------------------------------------
module CLOUD_COVER_LAYERS
 use CONSTANTS
 use ALGORITHM_CONSTANTS
 use PIXEL_COMMON

 implicit none
 public:: COMPUTE_CLOUD_FRACTION_3x3,       &
          COMPUTE_CLOUD_COVER_LAYERS

  contains

!------------------------------------------------------------------------------
! compute cloud fraction over a 3x3 array using the Bayesian probability
!------------------------------------------------------------------------------
 subroutine COMPUTE_CLOUD_FRACTION_3x3(jmin,jmax)

  integer, intent(in):: jmin,jmax
  integer :: i
  integer :: j
  integer :: i1
  integer :: i2
  integer :: j1
  integer :: j2
  integer :: ii
  integer :: jj
  integer :: Ngood
  integer :: Ncloud
  integer,parameter :: n=1
  real:: Cloud_Fraction_Uncer_Temp

  !--- initialize
  Cloud_Fraction_3x3 = Missing_Value_Real4
  Cloud_Fraction_Uncer_3x3 = Missing_Value_Real4

  line_loop: DO j=jmin, jmax - jmin + 1

    j1 = max(jmin,j-n)
    j2 = min(jmax,j+n)

    element_loop: DO i = 1, Num_Pix

      i1 = max(1,i-n)
      i2 = min(Num_Pix,i+n)

      !--- check for a bad pixel pixel
      if (Bad_Pixel_Mask(i,j) == sym%YES) then
        cycle
      endif

      !--- compute cloud amount
      Ngood = count(Posterior_Cld_Probability(i1:i2,j1:j2) /= Missing_Value_Real4)        

      !--- see if there are any valid points
      if (Ngood == 0 ) then
         cycle
      endif

      !--- Total Cloud Fraction
      Ncloud = count(Posterior_Cld_Probability(i1:i2,j1:j2) >= 0.5)        
      Cloud_Fraction_3x3(i,j) = real(Ncloud)/real(Ngood)

      !--- compute the uncertainty of the total cloud fraction
      Cloud_Fraction_Uncer_Temp = 0.0
      Cloud_Fraction_Uncer_3x3(i,j) = 0.0
      do ii = i1,i2
         do jj = j1,j2  
           if (Posterior_Cld_Probability(ii,jj) >= 0.5) then
             Cloud_Fraction_Uncer_Temp = 1.0 - Posterior_Cld_Probability(ii,jj)
           else 
             Cloud_Fraction_Uncer_Temp = Posterior_Cld_Probability(ii,jj)
           endif
           Cloud_Fraction_Uncer_3x3(i,j) = Cloud_Fraction_Uncer_3x3(i,j) +   &
                                           Cloud_Fraction_Uncer_Temp

         enddo
      enddo
      Cloud_Fraction_Uncer_3x3(i,j) = Cloud_Fraction_Uncer_3x3(i,j)  / Ngood

    end do element_loop
 end do line_loop

 end subroutine COMPUTE_CLOUD_FRACTION_3x3
!==============================================================================
!
! Compute Cloud Cover Layers which is defined here as the fraction of
! high, middle and low-level cloud in a NxN array centered on each pixel
!
! input (from global arrays)
!  Cld_Layer_ACHA = -1 (N/A), 0 (clear), 1 (low cloud), 2 (mid cloud) and 3 (high cloud)
!
! output (through global arrays)
!  Cld_Layer = integer flag for each pixel that identifies layer
!  High_Cloud_Fraction_3x3 = cloud fraction (0-1) for high cloud over 3x3 array
!  Mid_Cloud_Fraction_3x3 = cloud fraction (0-1) for mid cloud over 3x3 array
!  Low_Cloud_Fraction_3x3 = cloud fraction (0-1) for low cloud over 3x3 array
!
!  Notes
!  1. cloud fractions have values limited by 1/N^2. If N=1, values are 1/9,2/9...
!  2. Cloud_Fraction_3x3 is the total cloud amount.  This is computed from the
!     mask only.  If arrays have pixels with failed ACHA, H+M+L /= Total.
!  3. H/M/L defined in ACHA but boundaries are 440 and 680 hPA (ISCCP)
!==============================================================================
 subroutine COMPUTE_CLOUD_COVER_LAYERS(jmin,jmax)

  integer, intent(in):: jmin,jmax
  integer :: i
  integer :: j
  integer :: i1
  integer :: i2
  integer :: j1
  integer :: j2
  integer :: ii
  integer :: jj
  integer :: Ngood
  integer :: Ncloud
  integer,parameter :: N=1

  !--------------------------------------------------------------------
  ! compute pixel-level cloud layer flag
  !--------------------------------------------------------------------
  !--- initialize
  Cld_Layer_ACHA = -1

  line_loop_layer: DO j=jmin, jmax - jmin + 1
    element_loop_layer: DO i = 1, Num_Pix
 
      !--- check for a bad pixel pixel
      if (Bad_Pixel_Mask(i,j) == sym%YES) then
        cycle
      endif

      !------- determine cloud layer based on pressure
      if (Pc_ACHA(i,j) /= MISSING_VALUE_REAL4) then
        if (Pc_ACHA(i,j) <= HIGH_CLOUD_MAX_PRESSURE_THRESH) then
           Cld_Layer_ACHA(i,j) = 3
        elseif (Pc_ACHA(i,j) < LOW_CLOUD_MIN_PRESSURE_THRESH) then
           Cld_Layer_ACHA(i,j) = 2
        else
           Cld_Layer_ACHA(i,j) = 1
        endif
      endif

      !--- set clear pixels to appropriate cloud layer
      if ((Cld_Mask(i,j) == sym%CLEAR) .or. (Cld_Mask(i,j) == sym%PROB_CLEAR)) then
           Cld_Layer_Acha(i,j) = 0
      endif

    end do element_loop_layer
 end do line_loop_layer


 !--------------------------------------------------------------------
 ! compute pixel-level cloud cover for each layer over 3x3 array
 !--------------------------------------------------------------------

 !--- initialize
 High_Cloud_Fraction_3x3 = Missing_Value_Real4
 Mid_Cloud_Fraction_3x3 = Missing_Value_Real4
 Low_Cloud_Fraction_3x3 = Missing_Value_Real4

 line_loop_cover: DO j=jmin, jmax - jmin + 1

    j1 = max(jmin,j-n)
    j2 = min(jmax,j+n)

    element_loop_cover: DO i = 1, Num_Pix

      i1 = max(1,i-N)
      i2 = min(Num_Pix,i+N)

      !--- check for a bad pixel pixel
      if (Bad_Pixel_Mask(i,j) == sym%YES) then
        cycle
      endif

      !------- determine cloud layer based on pressure
      Cld_Layer_ACHA(i,j) = -1
      if (Pc_ACHA(i,j) /= MISSING_VALUE_REAL4) then
        if (Pc_ACHA(i,j) <= 440.0) then
           Cld_Layer_ACHA(i,j) = 3
        elseif (Pc_ACHA(i,j) < 680.0) then
           Cld_Layer_ACHA(i,j) = 2
        else
           Cld_Layer_ACHA(i,j) = 1
        endif
      endif

      !--- set clear pixels to appropriate cloud layer
      if ((Cld_Mask(i,j) == sym%CLEAR) .or. (Cld_Mask(i,j) == sym%PROB_CLEAR)) then
           Cld_Layer_Acha(i,j) = 0
      endif

      !--- compute number valid (good) data points in nxn array
      Ngood = count(Cld_Layer_Acha(i1:i2,j1:j2) >= 0)

      !--- see if there are any valid points, if not skip this pixel
      if (Ngood == 0 ) then
         cycle
      endif

      !--- High Cloud Fraction
      Ncloud = count(Cld_Layer_Acha(i1:i2,j1:j2) == 3)
      High_Cloud_Fraction_3x3(i,j) = real(Ncloud)/real(Ngood)

      !--- Mid Cloud Fraction
      Ncloud = count(Cld_Layer_Acha(i1:i2,j1:j2) == 2)
      Mid_Cloud_Fraction_3x3(i,j) = real(Ncloud)/real(Ngood)

      !--- Low Cloud Fraction
      Ncloud = count(Cld_Layer_Acha(i1:i2,j1:j2) == 1)
      Low_Cloud_Fraction_3x3(i,j) = real(Ncloud)/real(Ngood)

    end do element_loop_cover
 end do line_loop_cover

 end subroutine COMPUTE_CLOUD_COVER_LAYERS

!-----------------------------------------------------------
! end of MODULE
!-----------------------------------------------------------
end module CLOUD_COVER_LAYERS
