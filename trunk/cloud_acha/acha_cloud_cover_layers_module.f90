! $Id:$
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
module ACHA_CLOUD_COVER_LAYERS

  use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, acha_output_struct,symbol_acha, &
           acha_input_struct, acha_rtm_nwp_struct

 implicit none
 public:: COMPUTE_CLOUD_COVER_LAYERS

 real, private, PARAMETER:: MISSING_VALUE_REAL4 = -999.0
 integer, private, PARAMETER:: MISSING_VALUE_INTEGER4 = -999
 type(symbol_acha), private :: symbol

 contains

!------------------------------------------------------------------------------
! compute cloud fraction over a 3x3 array using the Bayesian probability
!------------------------------------------------------------------------------
 subroutine COMPUTE_CLOUD_COVER_LAYERS(Input, Output)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(acha_input_struct), intent(inout) :: Input
  type(acha_output_struct), intent(inout) :: Output

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
  real:: Total_Cloud_Fraction_Uncer_Temp
  integer,parameter :: N=1
  real, parameter:: HIGH_CLOUD_MAX_PRESSURE_THRESH = 440.0
  real, parameter:: LOW_CLOUD_MIN_PRESSURE_THRESH = 680.0


  !-------------------------------------------------------------------------------
  ! Total Cloud Fraction and Its Uncertainty
  !-------------------------------------------------------------------------------

  !--- initialize
  Output%Total_Cloud_Fraction = MISSING_VALUE_REAL4
  Output%Total_Cloud_Fraction_Uncer = MISSING_VALUE_REAL4

  line_loop_total: DO j= 1, Input%Number_of_Lines

    j1 = max(1,j-N)
    j2 = min(Input%Number_of_Lines,j+N)

    element_loop_total: DO i = 1, Input%Number_of_Elements

      i1 = max(1,i-N)
      i2 = min(Input%Number_of_Elements,i+N)

      !--- check for a bad pixel pixel
      if (Input%Invalid_Data_Mask(i,j) == Symbol%YES) then
        cycle
      endif

      !--- compute cloud amount
      Ngood = count(Input%Cloud_Probability(i1:i2,j1:j2) /= MISSING_VALUE_REAL4)        

      !--- see if there are any valid points
      if (Ngood == 0 ) then
         cycle
      endif

      !--- Total Cloud Fraction
      Ncloud = count(Input%Cloud_Probability(i1:i2,j1:j2) >= 0.5)        
      Output%Total_Cloud_Fraction(i,j) = real(Ncloud)/real(Ngood)

      !--- compute the uncertainty of the total cloud fraction
      Total_Cloud_Fraction_Uncer_Temp = 0.0
      Output%Total_Cloud_Fraction_Uncer(i,j) = 0.0
      do ii = i1,i2
         do jj = j1,j2  
           if (Input%Cloud_Probability(ii,jj) >= 0.5) then
             Total_Cloud_Fraction_Uncer_Temp = 1.0 - Input%Cloud_Probability(ii,jj)
           else 
             Total_Cloud_Fraction_Uncer_Temp = Input%Cloud_Probability(ii,jj)
           endif
           Output%Total_Cloud_Fraction_Uncer(i,j) = Output%Total_Cloud_Fraction_Uncer(i,j) +   &
                                           Total_Cloud_Fraction_Uncer_Temp

         enddo
      enddo
      Output%Total_Cloud_Fraction_Uncer(i,j) = Output%Total_Cloud_Fraction_Uncer(i,j)  / Ngood

    end do element_loop_total
 end do line_loop_total

!==============================================================================
!
! Compute Cloud Cover Layers which is defined here as the fraction of
! high, middle and low-level cloud in a NxN array centered on each pixel
!
! input (from global arrays)
!  Output%Cloud_Layer = -1 (N/A), 0 (clear), 1 (low cloud), 2 (mid cloud) and 3 (high cloud)
!
! output (through global arrays)
!  Cld_Layer = integer flag for each pixel that identifies layer
!  Output%High_Cloud_Fraction = cloud fraction (0-1) for high cloud over 3x3 array
!  Output%Mid_Cloud_Fraction = cloud fraction (0-1) for mid cloud over 3x3 array
!  Output%Low_Cloud_Fraction = cloud fraction (0-1) for low cloud over 3x3 array
!
!  Notes
!  1. cloud fractions have values limited by 1/N^2. If N=1, values are 1/9,2/9...
!  2. Cloud_Fraction_3x3 is the total cloud amount.  This is computed from the
!     mask only.  If arrays have pixels with failed ACHA, H+M+L /= Total.
!  3. H/M/L defined in ACHA but boundaries are 440 and 680 hPA (ISCCP)
!==============================================================================

  !--------------------------------------------------------------------
  ! compute pixel-level cloud layer flag
  !--------------------------------------------------------------------

  !--- initialize
  Output%Cloud_Layer = -1

  line_loop_layer: DO j = 1, Input%Number_of_Lines
    element_loop_layer: DO i = 1, Input%Number_of_Elements
 
      !--- check for a bad pixel pixel
      if (Input%Invalid_Data_Mask(i,j) == Symbol%YES) then
        cycle
      endif

      !------- determine cloud layer based on pressure
      if (Output%Pc(i,j) /= MISSING_VALUE_REAL4) then
        if (Output%Pc(i,j) <= HIGH_CLOUD_MAX_PRESSURE_THRESH) then
           Output%Cloud_Layer(i,j) = 3
        elseif (Output%Pc(i,j) < LOW_CLOUD_MIN_PRESSURE_THRESH) then
           Output%Cloud_Layer(i,j) = 2
        else
           Output%Cloud_Layer(i,j) = 1
        endif
      endif

      !--- set clear pixels to appropriate cloud layer
      if ((Input%Cloud_Mask(i,j) == Symbol%CLEAR) .or. (Input%Cloud_Mask(i,j) == Symbol%PROB_CLEAR)) then
           Output%Cloud_Layer(i,j) = 0
      endif

    end do element_loop_layer
 end do line_loop_layer


 !--------------------------------------------------------------------
 ! compute pixel-level cloud cover for each layer over 3x3 array
 !--------------------------------------------------------------------

 !--- initialize
 Output%High_Cloud_Fraction = MISSING_VALUE_REAL4
 Output%Mid_Cloud_Fraction = MISSING_VALUE_REAL4
 Output%Low_Cloud_Fraction = MISSING_VALUE_REAL4

 line_loop_cover: DO j = 1, Input%Number_of_Lines

    j1 = max(1,j-N)
    j2 = min(Input%Number_of_Lines,j+N)

    element_loop_cover: DO i = 1, Input%Number_of_Elements

      i1 = max(1,i-N)
      i2 = min(Input%Number_of_Elements,i+N)

      !--- check for a bad pixel pixel
      if (Input%Invalid_Data_Mask(i,j) == Symbol%YES) then
        cycle
      endif

      !--- set clear pixels to appropriate cloud layer
      if ((Input%Cloud_Mask(i,j) == Symbol%CLEAR) .or. (Input%Cloud_Mask(i,j) == Symbol%PROB_CLEAR)) then
           Output%Cloud_Layer(i,j) = 0
      endif

      !--- compute number valid (good) data points in nxn array
      Ngood = count(Output%Cloud_Layer(i1:i2,j1:j2) >= 0)

      !--- see if there are any valid points, if not skip this pixel
      if (Ngood == 0 ) then
         cycle
      endif

      !--- High Cloud Fraction
      Ncloud = count(Output%Cloud_Layer(i1:i2,j1:j2) == 3)
      Output%High_Cloud_Fraction(i,j) = real(Ncloud)/real(Ngood)

      !--- Mid Cloud Fraction
      Ncloud = count(Output%Cloud_Layer(i1:i2,j1:j2) == 2)
      Output%Mid_Cloud_Fraction(i,j) = real(Ncloud)/real(Ngood)

      !--- Low Cloud Fraction
      Ncloud = count(Output%Cloud_Layer(i1:i2,j1:j2) == 1)
      Output%Low_Cloud_Fraction(i,j) = real(Ncloud)/real(Ngood)

    end do element_loop_cover
 end do line_loop_cover

 end subroutine COMPUTE_CLOUD_COVER_LAYERS

!-----------------------------------------------------------
! end of MODULE
!-----------------------------------------------------------
end module ACHA_CLOUD_COVER_LAYERS
