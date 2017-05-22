! $Id: acha_cloud_cover_layers_module.f90 1523 2016-02-22 16:30:27Z wstraka $
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
! Cloud Layer Definition
!
! Layered Cloud Layer Meanings
! 00000000 = 0 = Clear
! 00000001 = 1 = Low
! 00000010 = 2 = Mid
! 00000100 = 4 = High
! 00000101 = 5 = High and Low
! 00000110 = 6 = High and Mid
! 00000011 = 3 = Mid and Low
! 00000111 = 7 = High and Mid and Low
!
! to extract
! to see if low cloud, ibits(layer,0,1)
! to see if mid cloud, ibits(layer,1,1)
! to see if high cloud, ibits(layer,2,1)
!
! FL Layer Meanings
! 0 = Clear
! 1 = SFC - 5 kft
! 2 = 5 - 10 kft
! 3 = 10 - 18 kft
! 4 = 18 - 24 kft
! 5 = 24  TOA kft
! 6 = 5+4
! 7 = 5+4+3
! 8 = 5+4+3+2
! 9 = 5+4+3+2+1
! 10 = 4+3
! 11 = 4+3+2
! 12 = 4+3+2+1
! 13 = 3+2
! 14 = 3+2+1
! 15 = 2+1

!
!--------------------------------------------------------------------------------------
module CCL_MODULE

  use CCL_SERVICES_MOD, only : &
           real4, int1, int4, ccl_output_struct,ccl_symbol_struct, &
           ccl_input_struct, ccl_diag_struct

 implicit none
 public:: COMPUTE_CLOUD_COVER_LAYERS
 private:: COMPUTE_BOX_WIDTH

 type(ccl_symbol_struct), private :: symbol

 include 'ccl_parameters.inc'

 contains

!------------------------------------------------------------------------------
! compute cloud fraction over a 3x3 array using the Bayesian probability
!------------------------------------------------------------------------------
 subroutine COMPUTE_CLOUD_COVER_LAYERS(Input, Symbol_In, Output, Diag)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(ccl_input_struct), intent(inout) :: Input
  type(ccl_symbol_struct), intent(inout) :: Symbol_In
  type(ccl_output_struct), intent(inout) :: Output
  type(ccl_diag_struct), intent(inout), optional :: Diag
  integer, save:: Diag_Warning_Flag = 0


  integer:: Num_Elems
  integer:: Num_Lines
  integer :: i,j                    !pixel indices
  integer :: i1,i2,j1,j2            !pixel indices of ccl box
  integer :: i11,i22,j11,j22        !pixels inices of area which ccl result is valid
  integer :: ni, nj                 !size of ccl box
  integer:: N                       !width of ccl box
  integer:: M                       !width of ccl region
  integer :: Num_Good
  integer :: Num_Cloud
  integer:: Num_High, Num_Mid, Num_Low, Num_Clear, Num_All
  real:: Clear_Fraction

  real (kind=4), dimension(:,:), allocatable:: Pixel_Uncertainty
  integer (kind=1), dimension(:,:), allocatable:: pos, len      !ibits arguments


  !-------------------------------------------------------------------------------
  ! Total Cloud Fraction and Its Uncertainty
  !-------------------------------------------------------------------------------

  !--- copy input symbol to a module-wide variable
  symbol = Symbol_In

  !--- initialize
  Output%Total_Cloud_Fraction = MISSING_VALUE_REAL4
  Output%Total_Cloud_Fraction_Uncer = MISSING_VALUE_REAL4
  Output%High_Cloud_Fraction = MISSING_VALUE_REAL4
  Output%Mid_Cloud_Fraction = MISSING_VALUE_REAL4
  Output%Low_Cloud_Fraction = MISSING_VALUE_REAL4
  Output%Cloud_Layer = MISSING_VALUE_INTEGER1

  !--- initialize diagnostic output
  if (present(Diag) .and. Diag_Warning_Flag == 0) then
      print *, "CLAVR-x / CCL ===>  Diagnostic Output Turned On"
      Diag_Warning_Flag = 1
  endif
  if (present(Diag)) Diag%Array_1 = Missing_Value_Real4
  if (present(Diag)) Diag%Array_2 = Missing_Value_Real4
  if (present(Diag)) Diag%Array_3 = Missing_Value_Real4

  !--- Determine Box Width 
  call COMPUTE_BOX_WIDTH(Input%Sensor_Resolution_KM,CCL_BOX_WIDTH_KM, N)
  call COMPUTE_BOX_WIDTH(Input%Sensor_Resolution_KM,CCL_SPACING_KM, M)

!==============================================================================
!
! Compute Cloud Cover Layers which is defined here as the fraction of
! high, middle and low-level cloud in a NxN array centered on each pixel
!
! output (through global arrays)
!  Output%Cloud_Fraction = cloud fraction (0-1) for all clouds over array
!  Output%Cloud_Fraction_Uncertainty = cloud fraction uncertainty (0-1)
!  Output%High_Cloud_Fraction = cloud fraction (0-1) for high cloud over array
!  Output%Mid_Cloud_Fraction = cloud fraction (0-1) for mid cloud over array
!  Output%Low_Cloud_Fraction = cloud fraction (0-1) for low cloud over array
!
!  Notes
!  1. cloud fractions have values limited by 1/N^2. If N=1, values are 1/9,2/9...
!  2. Cloud_Fraction is the total cloud amount.  This is computed from the
!     mask only.  If arrays have pixels with failed ACHA, H+M+L /= Total.
!  3. H/M/L boundaries in ccl_parameters.inc
!==============================================================================

  !--------------------------------------------------------------------
  ! compute pixel-level cloud layer flag and H/M/L masks
  !--------------------------------------------------------------------
  Num_Elems = Input%Number_of_Elements
  Num_Lines = Input%Number_of_Lines
  allocate(Pixel_Uncertainty(Num_Elems, Num_Lines))

  !--- make cloud fraction pixel level uncertainty
  Pixel_Uncertainty = MISSING_VALUE_REAL4
  where(Input%Cloud_Probability >= 0.5)
     Pixel_Uncertainty = 1.0 - Input%Cloud_Probability
  endwhere
  where(Input%Cloud_Probability < 0.5)
     Pixel_Uncertainty = Input%Cloud_Probability
  endwhere

  !------- identify clear and H/M/L pixels
  Output%Cloud_Layer = 0
  where (Input%Pc /= MISSING_VALUE_REAL4 .and. Input%Pc <= HIGH_CLOUD_MAX_PRESSURE_THRESH)
        Output%Cloud_Layer = ibset(Output%Cloud_Layer,2) 
  endwhere
  where (Input%Pc /= MISSING_VALUE_REAL4 .and.  & 
         Input%Pc > HIGH_CLOUD_MAX_PRESSURE_THRESH .and. &
         Input%Pc < LOW_CLOUD_MIN_PRESSURE_THRESH) 
         Output%Cloud_Layer = ibset(Output%Cloud_Layer,1) 
  endwhere
  where (Input%Pc >= LOW_CLOUD_MIN_PRESSURE_THRESH)
        Output%Cloud_Layer = ibset(Output%Cloud_Layer,0) !0001
  endwhere

  !--- account for depth of clouds using base
  if (USE_BASE) then
    where (ibits(Output%Cloud_Layer,2,1) == 1 .and. Input%Pc_Base > HIGH_CLOUD_MAX_PRESSURE_THRESH)
         Output%Cloud_Layer = ibset(Output%Cloud_Layer,1)
    endwhere
    where (ibits(Output%Cloud_Layer,1,1) == 1 .and. Input%Pc_Base > LOW_CLOUD_MIN_PRESSURE_THRESH)
         Output%Cloud_Layer = ibset(Output%Cloud_Layer,0) 
    endwhere
  endif

  !--- add in lower clouds obscured by higher clouds
  if (USE_LOWER) then
    where (Input%Pc_Lower /= MISSING_VALUE_REAL4 .and. Input%Pc_Lower <= HIGH_CLOUD_MAX_PRESSURE_THRESH)
        Output%Cloud_Layer = ibset(Output%Cloud_Layer,2) 
    endwhere
    where (Input%Pc_Lower /= MISSING_VALUE_REAL4 .and.  & 
           Input%Pc_Lower > HIGH_CLOUD_MAX_PRESSURE_THRESH .and. &
           Input%Pc_Lower < LOW_CLOUD_MIN_PRESSURE_THRESH) 
        Output%Cloud_Layer = ibset(Output%Cloud_Layer,1) 
    endwhere
    where (Input%Pc_Lower >= LOW_CLOUD_MIN_PRESSURE_THRESH)
        Output%Cloud_Layer = ibset(Output%Cloud_Layer,0) !0001
    endwhere
  endif

  where (Input%Cloud_Mask == Symbol%CLEAR .or. Input%Cloud_Mask == Symbol%PROB_CLEAR)
         Output%Cloud_Layer = 0   !0000
  endwhere 

 !--------------------------------------------------------------------
 ! compute pixel-level cloud cover for each layer over the box
 ! N = ccl box size
 ! M = ccl result spacing
 !--------------------------------------------------------------------

 line_loop_cover: do j = 1, Num_Lines, 2*M+1

    j1 = max(1,j-N)
    j2 = min(Num_Lines,j+N)
    j11 = max(1,j-M)
    j22 = min(Num_Lines,j+M)

    element_loop_cover: do i = 1, Num_Elems, 2*M+1

      i1 = max(1,i-N)
      i2 = min(Num_Elems,i+N)
      i11 = max(1,i-M)
      i22 = min(Num_Elems,i+M)

      ni = i2 - i1 + 1 
      nj = i2 - i1 + 1 

      if (allocated(pos)) deallocate(pos)
      if (allocated(len)) deallocate(len)

      allocate(pos(ni,nj), len(ni,nj))
      len = 1
     
      !--- check for a bad pixel pixel
      if (Input%Invalid_Data_Mask(i,j) == Symbol%YES) cycle

      !-- High
      pos = 2
      Num_High = int(sum(real(ibits(Output%Cloud_Layer(i1:i2,j1:j2),pos,len))))

      !-- Mid
      pos = 1
      Num_Mid = int(sum(real(ibits(Output%Cloud_Layer(i1:i2,j1:j2),pos,len))))

      !-- Low
      pos = 0
      Num_Low = int(sum(real(ibits(Output%Cloud_Layer(i1:i2,j1:j2),pos,len))))


      Num_Clear = count(Input%Cloud_Probability(i1:i2,j1:j2) < 0.5 .and. &
                        Input%Cloud_Probability(i1:i2,j1:j2) /= MISSING_VALUE_REAL4)         
      Num_Cloud = count(Input%Cloud_Probability(i1:i2,j1:j2) >= 0.5)        
      Num_Good = count(Input%Cloud_Probability(i1:i2,j1:j2) /= MISSING_VALUE_REAL4)        
      Num_All = Num_Cloud + Num_Clear


      !print *, "test ", Num_High, Num_Mid, Num_Low, Num_Clear, Num_Cloud, Num_Good, Num_All

      !--- see if there are any valid mask points, if not skip this pixel
      if (Num_Good < COUNT_MIN_CCL) then
         cycle
      endif

      !--- Total Cloud Fraction
      Output%Total_Cloud_Fraction(i11:i22,j11:j22) = real(Num_Cloud)/real(Num_Good)

      !--- compute the uncertainty of the total cloud fraction
      Output%Total_Cloud_Fraction_Uncer(i11:i22,j11:j22) = sum(Pixel_Uncertainty(i1:i2,j1:j2))/real(Num_Good)

      !--- see if there are any valid CCL points, if not skip this pixel
      if (Num_All == 0 ) then
         cycle
      endif

      !--- High Cloud Fraction
      Output%High_Cloud_Fraction(i11:i22,j11:j22) = real(Num_High)/real(Num_All)

      !--- Mid Cloud Fraction
      Output%Mid_Cloud_Fraction(i11:i22,j11:j22) = real(Num_Mid)/real(Num_All)

      !--- Low Cloud Fraction - ignore emissivity calc in ECA
      Output%Low_Cloud_Fraction(i11:i22,j11:j22) = real(Num_Low)/real(Num_All)
      !print *, "Low Cloud Fraction = ", real(Num_Low)/real(Num_All)

    end do element_loop_cover
 end do line_loop_cover

 if (allocated(Pixel_Uncertainty)) deallocate(Pixel_Uncertainty)

 end subroutine COMPUTE_CLOUD_COVER_LAYERS
!----------------------------------------------------------------------
!--- determine cirrus box width
!---
!--- Sensor_Resolution_KM = the nominal resolution in kilometers
!--- Box_Width_KM = the width of the desired box in kilometers
!--- Box_Half_Width = the half width of the box in pixel-space
!----------------------------------------------------------------------
subroutine COMPUTE_BOX_WIDTH(Sensor_Resolution_KM,Box_Width_KM, &
                             Box_Half_Width)

   real, intent(in):: Sensor_Resolution_KM
   integer, intent(in):: Box_Width_KM
   integer, intent(out):: Box_Half_Width

   if (Sensor_Resolution_KM <= 0.0) then
       Box_Half_Width = 20
   else
       Box_Half_Width = int((Box_Width_KM / Sensor_Resolution_KM) / 2)
   endif

end subroutine COMPUTE_BOX_WIDTH

!-----------------------------------------------------------
! end of MODULE
!-----------------------------------------------------------
end module CCL_MODULE
