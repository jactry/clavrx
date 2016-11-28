! $Id$
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
           real4, int1, int4, acha_output_struct,acha_symbol_struct, &
           acha_input_struct, acha_diag_struct

 implicit none
 public:: COMPUTE_CLOUD_COVER_LAYERS
!public:: ASOS_CLOUD_CLASSIFICATION
 private:: COMPUTE_BOX_WIDTH

 type(acha_symbol_struct), private :: symbol

 include 'ccl_parameters.inc'

 contains

!------------------------------------------------------------------------------
! compute cloud fraction over a 3x3 array using the Bayesian probability
!------------------------------------------------------------------------------
 subroutine COMPUTE_CLOUD_COVER_LAYERS(Input, Symbol_In, Output, Diag)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(acha_input_struct), intent(inout) :: Input
  type(acha_symbol_struct), intent(inout) :: Symbol_In
  type(acha_output_struct), intent(inout) :: Output
  type(acha_diag_struct), intent(inout), optional :: Diag

  integer:: Num_Elems
  integer:: Num_Lines
  integer :: i,j                    !pixel indices
  integer :: i1,i2,j1,j2            !pixel indices of big box
  integer :: i11,i22,j11,j22        !pixel indices of skipped pixels
  integer:: N
  integer:: M
  integer :: Num_Good
  integer :: Num_Cloud
  integer:: Num_High, Num_Mid, Num_Low, Num_Clear, Num_All
  real:: Clear_Fraction

  integer (kind=int1), dimension(:,:), allocatable:: Mask_High
  integer (kind=int1), dimension(:,:), allocatable:: Mask_Mid
  integer (kind=int1), dimension(:,:), allocatable:: Mask_Low
  integer (kind=int1), dimension(:,:), allocatable:: Mask_Clear
  real (kind=real4), dimension(:,:), allocatable:: Pixel_Uncertainty


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
  Output%ASOS_Cloud_Code = MISSING_VALUE_INTEGER1
  Output%ASOS_Cloud_ECA = MISSING_VALUE_REAL4
  Output%ASOS_Cloud_Zmin = MISSING_VALUE_REAL4
  Output%ASOS_Cloud_Zmax = MISSING_VALUE_REAL4

  !--- initialize diagnostic output
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
  allocate(Mask_High(Num_Elems, Num_Lines)) 
  allocate(Mask_Mid(Num_Elems, Num_Lines)) 
  allocate(Mask_Low(Num_Elems, Num_Lines))
  allocate(Mask_Clear(Num_Elems, Num_Lines))
  allocate(Pixel_Uncertainty(Num_Elems, Num_Lines))

  !--- make cloud fraction pixel level uncertainty
  Pixel_Uncertainty = MISSING_VALUE_REAL4
!  where(Input%Cloud_Probability >= 0.5)
!     Pixel_Uncertainty = 1.0 - Input%Cloud_Probability
!  endwhere
!  where(Input%Cloud_Probability < 0.5)
!     Pixel_Uncertainty = Input%Cloud_Probability
!  endwhere

  do j = 1,Input%Number_of_Lines
      do i = 1, Input%Number_of_Elements
        if (Input%Cloud_Probability(i,j) < 0.5) THEN
          Pixel_Uncertainty(i,j) = Input%Cloud_Probability(i,j)
        endif
        
        if (Input%Cloud_Probability(i,j) >= 0.5) THEN
          Pixel_Uncertainty(i,j) = 1.0 - Input%Cloud_Probability(i,j)
        endif
       
      enddo
  enddo




  !------- identify clear and H/M/L pixels
  Mask_High = 0
  Mask_Mid = 0
  Mask_Low = 0
  Mask_Clear = 0
  where (Output%Pc /= MISSING_VALUE_REAL4 .and. Output%Pc <= HIGH_CLOUD_MAX_PRESSURE_THRESH)
        Mask_High = 1
  endwhere
  where (Output%Pc < LOW_CLOUD_MIN_PRESSURE_THRESH .and. Output%Pc > HIGH_CLOUD_MAX_PRESSURE_THRESH)
         Mask_Mid = 1
  endwhere
  where (Output%Pc >= LOW_CLOUD_MIN_PRESSURE_THRESH)
         Mask_Low = 1
  endwhere
  where (Input%Cloud_Mask == Symbol%CLEAR .or. Input%Cloud_Mask == Symbol%PROB_CLEAR)
         Mask_Clear = 1
  endwhere 

 !--------------------------------------------------------------------
 ! compute pixel-level cloud cover for each layer over the box
 !--------------------------------------------------------------------

 line_loop_cover: DO j = 1, Num_Lines, 2*M+1

    j1 = max(1,j-N)
    j2 = min(Num_Lines,j+N)
    j11 = max(1,j-M)
    j22 = min(Num_Lines,j+M)

    element_loop_cover: DO i = 1, Num_Elems, 2*M+1

      i1 = max(1,i-N)
      i2 = min(Num_Elems,i+N)
      i11 = max(1,i-M)
      i22 = min(Num_Elems,i+M)

      !--- check for a bad pixel pixel
      if (Input%Invalid_Data_Mask(i,j) == Symbol%YES) cycle

      !--- count all of the pixels in each layer
      Num_High = int(sum(real(Mask_High(i1:i2,j1:j2))))
      Num_Mid = int(sum(real(Mask_Mid(i1:i2,j1:j2))))
      Num_Low = int(sum(real(Mask_Low(i1:i2,j1:j2))))
      Num_Clear = int(sum(real(Mask_Clear(i1:i2,j1:j2))))
      Num_Cloud = count(Input%Cloud_Probability(i1:i2,j1:j2) >= 0.5)        
      Num_Good = count(Input%Cloud_Probability(i1:i2,j1:j2) /= MISSING_VALUE_REAL4)        
      Num_All = Num_High + Num_Mid + Num_Low + Num_Clear

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

    end do element_loop_cover
 end do line_loop_cover

 !-----------------------------------------------------------------------------------------------
 ! ASOS Calculations
 !
 ! ASOS_Cloud_Code Key:
 !
 ! 0  - CLEAR or LOW Only
 ! 11 - MID SCT 
 ! 12 - MID BKN
 ! 13 - MID OVC
 ! 21 - HIGH SCT 
 ! 22 - HIGH BKN
 ! 23 - HIGH OVC
 !-----------------------------------------------------------------------------------------------
 if (ASOS_FLAG) then
 line_loop_asos: DO j = 1, Input%Number_of_Lines, 2*M+1

    j1 = max(1,j-N)
    j2 = min(Input%Number_of_Lines,j+N)
    j11 = max(1,j-M)
    j22 = min(Input%Number_of_Lines,j+M)

    element_loop_asos: DO i = 1, Input%Number_of_Elements,2*M+1

      i1 = max(1,i-N)
      i2 = min(Input%Number_of_Elements,i+N)
      i11 = max(1,i-M)
      i22 = min(Input%Number_of_Elements,i+M)

      !--- check for a bad pixel pixel
      if (Input%Invalid_Data_Mask(i,j) == Symbol%YES) cycle

      !--- check for an ocean pixel
      if (Input%Surface_Type(i,j) == Symbol%WATER_SFC) cycle

      !--- check for pixels outside of CONUS
      if (Input%Latitude(i,j) > ASOS_LAT_MAX) cycle
      if (Input%Latitude(i,j) < ASOS_LAT_MIN) cycle
      if (Input%Longitude(i,j) < ASOS_LON_MIN) cycle
      if (Input%Latitude(i,j) > ASOS_LON_MAX) cycle


      !Clear Fraction
      Clear_Fraction = sum(Mask_Clear(i1:i2,j1:j2))/real(Num_All)

      !--- OVC
      if (Clear_Fraction < OVC_CLEAR_FRACTION_THRESH) then
         if (Output%High_Cloud_Fraction(i,j) > Output%Mid_Cloud_Fraction(i,j)) then
             Output%ASOS_Cloud_Code(i11:i22,j11:j22)= 23
         else
             Output%ASOS_Cloud_Code(i11:i22,j11:j22) = 13
         endif
      endif

      !--- BKN
      if (Clear_Fraction >= OVC_CLEAR_FRACTION_THRESH .and. Clear_Fraction < BKN_CLEAR_FRACTION_THRESH) then
         if (Output%Low_Cloud_Fraction(i,j) > LOW_CLOUD_FRACTION_THRESH) then
            if (Output%High_Cloud_Fraction(i,j) < 0.10) then
                 Output%ASOS_Cloud_Code(i11:i22,j11:j22) = 22
            else
                 Output%ASOS_Cloud_Code(i11:i22,j11:j22) = 12
            endif
         else
            if (Output%High_Cloud_Fraction(i,j) > Output%Mid_Cloud_Fraction(i,j)) then
                 Output%ASOS_Cloud_Code(i11:i22,j11:j22) = 22
            else
                 Output%ASOS_Cloud_Code(i11:i22,j11:j22) = 12
            endif
         endif
      endif

      !--- SCT
      if (Clear_Fraction > BKN_CLEAR_FRACTION_THRESH) then
          if (Output%High_Cloud_Fraction(i,j) > 0.0) then
              Output%ASOS_Cloud_Code(i11:i22,j11:j22) = 21
          elseif (Output%Mid_Cloud_Fraction(i,j) > 0.0) then
              Output%ASOS_Cloud_Code(i11:i22,j11:j22) = 11
          endif
      endif

      if (Output%ASOS_Cloud_Code(i,j) > 20) then
         Output%ASOS_Cloud_ECA(i11:i22,j11:j22) = sum(Mask_High(i1:i2,j1:j2)*Output%Ec(i1:i2,j1:j2))/real(Num_High)
         Output%ASOS_Cloud_Zmin(i11:i22,j11:j22) = minval(Mask_High(i1:i2,j1:j2)*Output%Zc(i1:i2,j1:j2))
         Output%ASOS_Cloud_Zmax(i11:i22,j11:j22) = maxval(Mask_High(i1:i2,j1:j2)*Output%Zc(i1:i2,j1:j2))
      elseif (Output%ASOS_Cloud_Code(i,j) > 10) then
         Output%ASOS_Cloud_ECA(i11:i22,j11:j22) = sum(Mask_Mid(i1:i2,j1:j2)*Output%Ec(i1:i2,j1:j2))/real(Num_Mid)
         Output%ASOS_Cloud_Zmin(i11:i22,j11:j22) = minval(Mask_Mid(i1:i2,j1:j2)*Output%Zc(i1:i2,j1:j2))
         Output%ASOS_Cloud_Zmax(i11:i22,j11:j22) = maxval(Mask_Mid(i1:i2,j1:j2)*Output%Zc(i1:i2,j1:j2))
      else
         Output%ASOS_Cloud_ECA(i11:i22,j11:j22) = MISSING_VALUE_REAL4
         Output%ASOS_Cloud_Zmin(i11:i22,j11:j22) = MISSING_VALUE_REAL4
         Output%ASOS_Cloud_Zmax(i11:i22,j11:j22) = MISSING_VALUE_REAL4
      endif

    end do element_loop_asos
 end do line_loop_asos
 endif

 if (allocated(Mask_High)) deallocate(Mask_High)
 if (allocated(Mask_Mid)) deallocate(Mask_Mid)
 if (allocated(Mask_Low)) deallocate(Mask_Low)
 if (allocated(Mask_Clear)) deallocate(Mask_Clear)
 if (allocated(Pixel_Uncertainty)) deallocate(Pixel_Uncertainty)

 end subroutine COMPUTE_CLOUD_COVER_LAYERS
!!-------------------------------------------------------------------------------------
!! ASOS_CLOUD_CLASSIFICATION
!!
!! Purpose: Use Cld_Mask, CTP and Ec to classify each pixel into the ASOS Cloud
!! Classes.
!!
!! Reference: The ATBD for the GOES Automated Surface Observing System (ASOS)
!!            Satellite Cloud Product.
!!
!!
!! The ASOS_Cloud_Class which is one of 10 values (1-10) as defined by
!! the table below. 
!!
!! The table below comes from ATBD.
!! Values in () are integer values used here. A value of 0 will be used for
!! pixels where this classification failed. Fill_Value = -128
!! 
!!                   Ec<0.33    0.33<Ec<0.65  0.65<Ec<0.95  0.95<ec
!!
!! CTP <= 400          b1(7)          c1(8)        d1(9)      e1(10)
!! 400<CTP<631         b2(3)          c2(4)        d2(5)      e2(6)
!! 631<CTP<950         -              -             -         f(2)
!! CTP > 950           a(1)           -             -          -
!!
!!-------------------------------------------------------------------------------------
!subroutine ASOS_CLOUD_CLASSIFICATION (Input, Symbol_In, Output)

! !===============================================================================
! !  Argument Declaration
! !==============================================================================

! type(acha_input_struct), intent(inout) :: Input
! type(acha_symbol_struct), intent(inout) :: Symbol_In
! type(acha_output_struct), intent(inout) :: Output

! integer :: i
! integer :: j
! real, parameter:: HIGH_CLOUD_MAX_PRESSURE_THRESH = 400.0
! real, parameter:: MID_CLOUD_MAX_PRESSURE_THRESH =  631.0
! real, parameter:: LOW_CLOUD_MAX_PRESSURE_THRESH =  950.0
! real, parameter, dimension(3):: Ec_Bounds = (/0.33,0.65,0.95/)

! !--- initialize
! Output%ASOS_Cloud_Class = MISSING_VALUE_INTEGER1

! line_loop: do j = 1, Input%Number_of_Lines
!   element_loop: do i = 1, Input%Number_of_Elements

!   !--- check for a bad pixel pixel
!   if (Input%Invalid_Data_Mask(i,j) == Symbol%YES) cycle

!   !--- check for failed height retrievals
!   if (Output%Pc(i,j) == MISSING_VALUE_REAL4 .or. &
!       Output%Ec(i,j) == MISSING_VALUE_REAL4) then
!         Output%ASOS_Cloud_Class(i,j) =  0
!   endif

!   !--- check for clear (a/1)
!   if (Input%Cloud_Type(i,j) == symbol%CLEAR_TYPE .or. &
!       Input%Cloud_Type(i,j) == symbol%PROB_CLEAR_TYPE .or. &
!       Output%Pc(i,j) > LOW_CLOUD_MAX_PRESSURE_THRESH) then

!       Output%ASOS_Cloud_Class(i,j) =  1
!       cycle

!   endif   

!   !--- check for low cloud (f/2)
!   if (Output%Pc(i,j) >= MID_CLOUD_MAX_PRESSURE_THRESH) then
!      Output%ASOS_Cloud_Class(i,j)  = 2
!   endif

!   !-- Do high and mid cloud classes
!   if (Output%PC(i,j) < HIGH_CLOUD_MAX_PRESSURE_THRESH) then
!      Output%ASOS_Cloud_Class(i,j)  = 7
!   else
!      Output%ASOS_Cloud_Class(i,j)  = 3
!   endif

!   !--- add no offset for thin clouds
!   if (Output%Ec(i,j) < Ec_Bounds(1))  cycle

!   !-- add an offset of 1 for moderate clouds
!   if (Output%Ec(i,j) >= Ec_Bounds(1) .and. Output%Ec(i,j) < Ec_Bounds(2))  then
!         Output%ASOS_Cloud_Class(i,j) = Output%ASOS_Cloud_Class(i,j) + 1
!         cycle
!   endif

!   !-- add an offset of 2 for thick clouds
!   if (Output%Ec(i,j) >= Ec_Bounds(2) .and. Output%Ec(i,j) < Ec_Bounds(3))  then
!         Output%ASOS_Cloud_Class(i,j) = Output%ASOS_Cloud_Class(i,j) + 2
!         cycle
!   endif

!   !-- add an offset of 3 for opaque clouds
!   if (Output%Ec(i,j) > Ec_Bounds(3))  then
!         Output%ASOS_Cloud_Class(i,j) = Output%ASOS_Cloud_Class(i,j) + 3
!         cycle
!   endif

!   end do element_loop
! end do line_loop


!end subroutine ASOS_CLOUD_CLASSIFICATION
!----------------------------------------------------------------------
!--- determine cirrus box width
!---
!--- Sensor_Resolution_KM = the nominal resolution in kilometers
!--- Box_Width_KM = the width of the desired box in kilometers
!--- Box_Half_Width = the half width of the box in pixel-space
!----------------------------------------------------------------------
subroutine COMPUTE_BOX_WIDTH(Sensor_Resolution_KM, &
                             Box_Width_KM, &
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
end module ACHA_CLOUD_COVER_LAYERS
