!$Id$
module UNIVERSAL_CLOUD_TYPE_MODULE 
!====================================================================
! Module Name: UNIVERSAL_CLOUD_TYPE_MODULE
!
! Function:  House routines needed for the universal_cloud_type routine
!
! Description: 
! This module houses the routines associated with the a cloud
! type routine that should run on all sensors.  This not the
! AWG cloud typing algorithm. This code is meant only for CLAVR-x
! and GSIP which are not yet able to run the AWG algorithms.
! 
! The goal of this algorithm is to make a phase/type that is 
! applicable for PATMOS-x climate processing of AVHRR/GOES/MODIS/VIIRS.
! It prioritizes consistency with common channels over
! absolute accuracy
!   
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: 
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
!
!  Step 1 - Setup Parameters Needed for Phasing
!  Step 2 - Compute Ice Probability
!  Step 3 - Compute Water Flag
!  Step 4 - Compute Cirrus Flag
!  Step 5 - Compute Overlap Flag
!  Step 6 - Apply LRC Logic
!  Step 7 - Deal with remaining unknown pixels
!  Step 8 - Derive Type from Phase with additional logic
!====================================================================
  use CONSTANTS
  use PIXEL_COMMON
  use NWP_COMMON
  use RT_UTILITIES
  use RTM_COMMON

  implicit none

  public::  UNIVERSAL_CLOUD_TYPE
  private:: COMPUTE_ICE_PROBABILITY_BASED_ON_TEMPERATURE
  private:: OPAQUE_CLOUD_HEIGHT_LOCAL
  private:: H2O_CLOUD_HEIGHT_LOCAL

  !--- module wide variables to make more efficient argument passing internally
  real (kind=real4), dimension(:), private, pointer:: T_Prof
  real (kind=real4), dimension(:), private, pointer:: Z_Prof
  real (kind=real4), dimension(:), private, allocatable:: P_Prof
  real (kind=real4), dimension(:), private, pointer:: Rad_11um_BB_Profile
  real (kind=real4), dimension(:), private, pointer:: Rad_H2O_BB_Profile

  !--- channel mapping for each system here
#ifdef ISCLAVRX
  include 'akh_cloud_type_clavrx_include_1.inc'
#endif
#ifdef ISGSIP
  include 'akh_cloud_type_gsip_include_1.inc'
#endif

  !--- thresholds
  real(kind=real4), parameter, private:: Ice_Temperature_Min = 243.0
  real(kind=real4), parameter, private:: Ice_Temperature_Max = 263.0
  real(kind=real4), parameter, private:: Rad_11um_Thresh = 2.0
  real(kind=real4), parameter, private:: Rad_H2O_Thresh = 0.25
  real, private, parameter:: Zclr_H2O_Moist_Thresh = 8.0
  real, private, parameter:: Bt_Ch27_Ch31_Covar_Cirrus_Moist_Thresh = 0.25
  real, private, parameter:: Bt_Ch27_Ch31_Covar_Cirrus_Thresh = 1.0 !0.5
  real, private, parameter:: Bt_Ch31_Btd_Ch31_Ch32_Covar_Cirrus_Thresh = -1.0
  real, private, parameter:: Beta_11um_85um_Ice_Thresh = 1.10
  real, private, parameter:: Beta_11um_12um_Overlap_Thresh = 0.95
  real, private, parameter:: Beta_11um_133um_Overlap_Thresh = 0.70
  real, private, parameter:: Emiss_67um_Cirrus_Thresh = 0.05
  real, private, parameter:: Emiss_133um_Cirrus_Thresh = 0.02
  real, private, parameter:: Fmft_Cold_Offset = 0.5 !K
  real, private, parameter:: Fmft_Cirrus_Thresh = 1.0 !K
  real, private, parameter:: Bt_11um_Std_Cirrus_Thresh = 4.0 !K

  CONTAINS
!====================================================================
!  record cvs version as a global variable for output to hdf
!====================================================================
#ifdef ISCLAVRX
subroutine SET_CLOUD_TYPE_VERSION()
   Cloud_Type_Version = "$Id$"
end subroutine SET_CLOUD_TYPE_VERSION
#endif

!====================================================================
! Subroutine Name: UNIVERSAL_CLOUD_TYPE
!
! Function: Compute cloud phase, cloud type and opaque and h2o cloud 
!           height/temperature/pressure
!
! Description: 
! This not the AWG cloud typing algorithm. 
! This code is meant only for CLAVR-x and GSIP which are not yet 
! able to run the AWG algorithms.
!   
! Dependencies: RTM data structure needed
!
! Restrictions:  11 um channel needs to be valid
!
! Reference: 
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
!====================================================================
subroutine UNIVERSAL_CLOUD_TYPE(Line_Start,Line_End)

  implicit none

  integer, intent(in):: Line_Start
  integer, intent(in):: Line_End
  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Num_Elem
  integer:: Num_Line
  integer:: Elem_Start
  integer:: Elem_End
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  integer:: Vza_Idx
  real:: Fmft
  integer:: i1,i2,j1,j2
  real:: T_Opa
  integer:: Elem_Lrc_Idx
  integer:: Line_Lrc_Idx
  integer (kind=int1), dimension(:,:), pointer:: Invalid_Pixel_Mask
  integer (kind=int1), dimension(:,:), pointer:: Cloud_Mask
  integer (kind=int1), dimension(:,:), pointer:: Cloud_Phase
  integer (kind=int1), dimension(:,:), pointer:: Cloud_Type
  integer (kind=int1), dimension(:), pointer:: Channel_On_Flag
  real (kind=real4), dimension(:,:), pointer:: Ref_16um
  real (kind=real4), dimension(:,:), pointer:: Ref_375um
  real (kind=real4), dimension(:,:), pointer:: Rad_11um
  real (kind=real4), dimension(:,:), pointer:: Rad_11um_Clear
  real (kind=real4), dimension(:,:), pointer:: Bt_H2O
  real (kind=real4), dimension(:,:), pointer:: Bt_H2O_Stddev
  real (kind=real4), dimension(:,:), pointer:: Rad_H2O
  real (kind=real4), dimension(:,:), pointer:: Rad_H2O_Clear
  real (kind=real4), dimension(:,:), pointer:: Bt_H2O_Clear
  real (kind=real4), dimension(:,:), pointer:: Covar_H2O_Window
  real (kind=real4), dimension(:,:), pointer:: Bt_11um
  real (kind=real4), dimension(:,:), pointer:: Bt_11um_Clear
  real (kind=real4), dimension(:,:), pointer:: Bt_11um_Max
  real (kind=real4), dimension(:,:), pointer:: Bt_11um_Min
  real (kind=real4), dimension(:,:), pointer:: Bt_11um_Mean
  real (kind=real4), dimension(:,:), pointer:: Bt_11um_Stddev
  real (kind=real4), dimension(:,:), pointer:: Bt_12um
  real (kind=real4), dimension(:,:), pointer:: Bt_12um_Clear
  real (kind=real4), dimension(:,:), pointer:: Tcld_Opa
  real (kind=real4), dimension(:,:), pointer:: Zcld_Opa
  real (kind=real4), dimension(:,:), pointer:: Pcld_Opa
  real (kind=real4), dimension(:,:), pointer:: Tcld_H2O
  real (kind=real4), dimension(:,:), pointer:: Zcld_H2O
  real (kind=real4), dimension(:,:), pointer:: Pcld_H2O
  real (kind=real4), dimension(:,:), pointer:: Zclr_H2O
  real (kind=real4), dimension(:,:), pointer:: Emiss_Window_Tropo
  real (kind=real4), dimension(:,:), pointer:: Emiss_H2O_Tropo
  real (kind=real4), dimension(:,:), pointer:: Beta_11um_85um_Tropo
  real (kind=real4), dimension(:,:), pointer:: Beta_11um_12um_Tropo
  real (kind=real4), dimension(:,:), pointer:: Beta_11um_133um_Tropo
  integer (kind=int1), pointer:: Sfc_Level
  integer (kind=int1), pointer:: Tropo_Level
  integer, dimension(:,:), pointer:: Vza_Idx_Rtm
  integer:: Number_Rtm_Levels

  integer, dimension(:,:), allocatable:: Cirrus_Flag
  integer, dimension(:,:), allocatable:: Water_Flag
  integer, dimension(:,:), allocatable:: Overlap_Flag
  real, dimension(:,:), allocatable:: Ice_Probability
  integer:: Fire_Flag

  integer, parameter:: n_sub = 5
! integer, dimension(n_sub,n_sub):: mask_sub, type_sub
  integer:: water_count
  integer:: ice_count

  integer(kind=int1) :: smoke_pixel
  integer(kind=int1) :: dust_pixel

  !--- begin executable code

  !------ set up pointers to global arrays 
#ifdef ISCLAVRX
  include 'akh_cloud_type_clavrx_include_2.inc'
#endif   
#ifdef ISGSIP
  include 'akh_cloud_type_gsip_include_2.inc'
#endif   

  !--- allocate
  Num_Line = size(Cloud_Type(0,:))

  allocate(Cirrus_Flag(Num_Elem,Num_Line))
  allocate(Water_Flag(Num_Elem,Num_Line))
  allocate(Overlap_Flag(Num_Elem,Num_Line))
  allocate(Ice_Probability(Num_Elem,Num_Line))

  !--- initialize
  Overlap_Flag = sym%NO
  Cirrus_Flag = sym%NO
  Water_Flag = sym%NO
  Ice_Probability = Missing_Value_Real4
  Fire_Flag = sym%NO
  Cloud_Phase = sym%UNKNOWN_PHASE

  !------------------------------------------------------------------
  ! Step #1: Check for non-cloud conditions and 
  !          compute parameters needed in future steps
  !------------------------------------------------------------------
  Line_Loop_Setup: do Line_Idx = Line_Start, Line_End
  Element_Loop_Setup: do Elem_Idx = 1, Num_Elem

     !--- skip bad pixels
     if (Invalid_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !--- save indices
     Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
     Vza_Idx = Vza_Idx_Rtm(Elem_Idx,Line_Idx)

     !--- point to RTM data structure components
#ifdef ISCLAVRX
     include 'akh_cloud_type_clavrx_include_3.inc'
#endif

#ifdef ISGSIP
     include 'akh_cloud_type_gsip_include_3.inc'
#endif

     !-------------------------------------------------------------
     ! Determine if a non-cloud type has been determined in the
     ! cloud mask, if so, set the type flag and exit
     !-------------------------------------------------------------
     Fire_Flag = BTEST(Cld_Test_Vector_Packed(2,Elem_Idx,Line_Idx), 7)
     if (Fire_Flag == sym%YES) then
        Cloud_Type(Elem_Idx,Line_Idx) = sym%FIRE_TYPE
        Cloud_Phase(Elem_Idx,Line_Idx) = sym%UNKNOWN_PHASE
        cycle
     endif

     !--- set clear to clear phase
     if (Cloud_Mask(Elem_Idx,Line_Idx) == sym%CLEAR) then
          Cloud_Phase(Elem_Idx,Line_Idx) = sym%CLEAR_PHASE
          cycle
     endif

     !--- set probably clear to clear phase
     if (Cloud_Mask(Elem_Idx,Line_Idx) == sym%PROB_CLEAR) then
          Cloud_Phase(Elem_Idx,Line_Idx) = sym%CLEAR_PHASE
          cycle
     endif

     !-------------------------------------------------------------
     ! Determine Opaque Cloud Temp and Height
     !-------------------------------------------------------------
     if (Channel_On_Flag(Chan_Idx_11um) == sym%NO) cycle

     call OPAQUE_CLOUD_HEIGHT_LOCAL(Rad_11um(Elem_Idx,Line_Idx), &
                            Tropo_Level, &
                            Sfc_Level, &
                            Tcld_Opa(Elem_Idx,Line_Idx), &
                            Zcld_Opa(Elem_Idx,Line_Idx), &
                            Pcld_Opa(Elem_Idx,Line_Idx))

     !------------------------------------------------------------
     ! Determine Water Vapor Height
     !------------------------------------------------------------
     if (Channel_On_Flag(Chan_Idx_67um) == sym%NO) cycle

     ! Determine Clear-sky Water Vapor Weighting Function Peak Height
     call H2O_CLEAR_SKY_PEAK_HEIGHT(Bt_H2O_Clear(Elem_Idx,Line_Idx), &
                                    Tropo_Level,  &
                                    Sfc_Level, &
                                    Zclr_H2O(Elem_Idx,Line_Idx))

     call H2O_CLOUD_HEIGHT_LOCAL(Rad_11um(Elem_Idx,Line_Idx),  &
                            Rad_11um_Clear(Elem_Idx,Line_Idx), &
                            Rad_H2O(Elem_Idx,Line_Idx), &
                            Rad_H2O_Clear(Elem_Idx,Line_Idx), &
                            Covar_H2O_Window(Elem_Idx,Line_Idx), &
                            Tropo_Level, &
                            Sfc_Level, &
                            Tcld_H2O(Elem_Idx,Line_Idx), &
                            Zcld_H2O(Elem_Idx,Line_Idx), &
                            Pcld_H2O(Elem_Idx,Line_Idx))

  enddo Element_Loop_Setup
  enddo Line_Loop_Setup

  !------------------------------------------------------------------------------
  ! Step #2: Compute Ice Probability based on best estimate of Cloud Temperature
  !------------------------------------------------------------------------------
  Line_Loop_Height: do Line_Idx = Line_Start, Line_End
  Element_Loop_Height: do Elem_Idx = 1, Num_Elem

     T_Opa = Tcld_H2O(Elem_Idx,Line_Idx)

     if (T_Opa == Missing_Value_Real4) T_Opa = Tcld_Opa(Elem_Idx,Line_Idx)

     if (T_Opa == Missing_Value_Real4) then 
        T_Opa = Bt_11um(Elem_Idx,Line_Idx)
     endif

     Ice_Probability(Elem_Idx,Line_Idx) = COMPUTE_ICE_PROBABILITY_BASED_ON_TEMPERATURE(T_Opa) 

  enddo Element_Loop_Height
  enddo Line_Loop_Height

  !----------------------------------------------------------------------------
  ! Step #3:  spectrally detect water clouds
  !
  !  These tests need to be conservative
  !----------------------------------------------------------------------------
  Water_Flag = sym%NO
  Line_Loop_WATER: do Line_Idx = Line_Start, Line_End
  Element_Loop_WATER: do Elem_Idx = 1, Num_Elem


   if (Cloud_Phase(Elem_Idx,Line_Idx) == sym%CLEAR_PHASE) cycle

   if (Tc_Opaque_Cloud(Elem_Idx,Line_Idx) > 240.0) then

   !---- 1.6 um Spectral Test for Water
   if (Channel_On_Flag(Chan_Idx_16um) == sym%YES) then
        if (Solzen(Elem_Idx,Line_Idx) < 80.0) then
          if (ch(6)%Ref_Toa_Clear(Elem_Idx,Line_Idx) < 20.0) then
          if (Ref_16um(Elem_Idx,Line_Idx) > 30.0) then
             Water_Flag(Elem_Idx,Line_Idx) = sym%YES
          endif
          endif
       endif
   endif

   !---- 3.75 um Spectral Test for Water - If triggered, ignore cirrus tests
   if (Channel_On_Flag(Chan_Idx_375um) == sym%YES) then
        if (Solzen(Elem_Idx,Line_Idx) < 80.0) then
         if (ch(20)%Sfc_Emiss(Elem_Idx,Line_Idx) > 0.90) then
          if (Ref_375um(Elem_Idx,Line_Idx) > 20.0) then
             Water_Flag(Elem_Idx,Line_Idx) = sym%YES
          endif
        elseif (Solzen(Elem_Idx,Line_Idx) > 90.0) then
          if (Ref_375um(Elem_Idx,Line_Idx) > 5.0) then
             Water_Flag(Elem_Idx,Line_Idx) = sym%YES
         endif
       endif
       endif
   endif

   endif

   if (Water_Flag(Elem_Idx,Line_Idx) == sym%YES) then
       Ice_Probability(Elem_Idx,Line_Idx) = 0.0
   endif

  enddo Element_Loop_WATER
  enddo Line_Loop_WATER

  !----------------------------------------------------------------------------
  ! Step #4:  detect cirrus clouds and assume they are ice
  !----------------------------------------------------------------------------
  Cirrus_Flag = sym%NO

  Line_Loop_CIRRUS: do Line_Idx = Line_Start, Line_End
  Element_Loop_CIRRUS: do Elem_Idx = 1, Num_Elem

   !---- don't detect cirrus if clear
   if (Cloud_Phase(Elem_Idx,Line_Idx) == sym%CLEAR_PHASE) cycle

   !---- don't detect cirrus if very high 11 um std deviation
   if (Bt_11um_Stddev(Elem_Idx,Line_Idx) > Bt_11um_Std_Cirrus_Thresh) cycle

     !--- split window
     if ((Channel_On_Flag(Chan_Idx_11um) == sym%YES) .and. &
        (Channel_On_Flag(Chan_Idx_12um) == sym%YES)) then

         Fmft  = (Bt_11um_Clear(Elem_Idx,Line_Idx) - Bt_12um_Clear(Elem_Idx,Line_Idx)) *  &
                 (Bt_11um(Elem_Idx,Line_Idx) - 260.0) / &
                 (Bt_11um_Clear(Elem_Idx,Line_Idx) - 260.0)

        if (Bt_11um_Clear(Elem_Idx,Line_Idx) <= 265.0) then
           Fmft = Fmft_Cold_Offset
        endif

        Fmft = (Bt_11um(Elem_Idx,Line_Idx) - Bt_12um(Elem_Idx,Line_Idx)) - Fmft

        if (Fmft > Fmft_Cirrus_Thresh) then
               Cirrus_Flag(Elem_Idx,Line_Idx) = sym%YES
        endif
     endif

      !--- 6.7 um covariance
      if (Channel_On_Flag(Chan_Idx_67um) == sym%YES) then

       if ((Covar_H2O_Window(Elem_Idx,Line_Idx) > 1.5) .and. Bt_Ch27_Max_3x3(Elem_Idx,Line_Idx) < 250.0) then
              Cirrus_Flag(Elem_Idx,Line_Idx) = sym%YES
       endif
 
      endif

     if (Cirrus_Flag(Elem_Idx,Line_Idx) == sym%YES) then
       Ice_Probability(Elem_Idx,Line_Idx) = 1.0
     endif

  enddo Element_Loop_CIRRUS
  enddo Line_Loop_CIRRUS

  !--------------------------------------------------------------------
  ! Step #5: Compute Overlap Flag
  !--------------------------------------------------------------------
  Overlap_Flag = sym%NO

  !--------------------------------------------------------------------
  ! Step #6: Compute Phase for Pixels that are LRCs  
  !--------------------------------------------------------------------
  Line_Loop_LRC: do Line_Idx = Line_Start, Line_End
  Element_Loop_LRC: do Elem_Idx = 1, Num_Elem

     if (Cloud_Phase(Elem_Idx,Line_Idx) == sym%CLEAR_PHASE) cycle

     Elem_Lrc_Idx = i_lrc(Elem_Idx,Line_Idx)
     Line_Lrc_Idx = j_lrc(Elem_Idx,Line_Idx)

     if (Elem_Lrc_Idx > 0 .and. Line_Lrc_Idx > 0) then
         if (Ice_Probability(Elem_Lrc_Idx,Line_Lrc_Idx) /= Missing_Value_Real4) then
            Ice_Probability(Elem_Idx,Line_Idx) = Ice_Probability(Elem_Lrc_Idx,Line_Lrc_Idx)
         endif
         if (Cirrus_Flag(ELem_Lrc_Idx,Line_Lrc_Idx) == sym%YES) then
           Ice_Probability(Elem_Idx,Line_Idx) = 1.0
         endif
     endif

     !--- phase based on ice cloud probability
     if (Ice_Probability(Elem_Idx,Line_Idx) > 0.50) then
          Cloud_Phase(Elem_Idx,Line_Idx) = sym%ICE_PHASE
     else
          Cloud_Phase(Elem_Idx,Line_Idx) = sym%WATER_PHASE
     endif

  enddo Element_Loop_LRC
  enddo Line_Loop_LRC

  !----------------------------------------------------------------------------
  ! Step #7:  Phase remaining unknown pixels by majority in area
  !----------------------------------------------------------------------------
  Line_Loop_LRC3: do Line_Idx = Line_Start, Line_End
  Element_Loop_LRC3: do Elem_Idx = 1, Num_Elem

    if (Cloud_Phase(Elem_Idx,Line_Idx) == sym%UNKNOWN_PHASE) then

        i1 = max(1,min(Num_Elem,Elem_Idx - n_sub))
        i2 = max(1,min(Num_Elem,Elem_Idx + n_sub))
        j1 = max(Line_Start,min(Line_End,Line_Idx - n_sub))
        j2 = max(Line_Start,min(Line_End,Line_Idx + n_sub))

        water_count = sum(Cloud_Phase(i1:i2,j1:j2),MASK=Cloud_Phase(i1:i2,j1:j2)==sym%WATER_PHASE)
        ice_count = sum(Cloud_Phase(i1:i2,j1:j2),MASK=Cloud_Phase(i1:i2,j1:j2)==sym%ICE_PHASE)

        if (ice_count > water_count) then
          Cloud_Phase(Elem_Idx,Line_Idx) = sym%ICE_PHASE
        else
          Cloud_Phase(Elem_Idx,Line_Idx) = sym%WATER_PHASE
        endif 

    endif

  enddo Element_Loop_LRC3
  enddo Line_Loop_LRC3

  !-----------------------------------------------------------------------------------
  ! set all water phase clouds with cold tops as supercooled water
  !-----------------------------------------------------------------------------------
  where(Cloud_Phase == sym%WATER_PHASE .and. Tcld_Opa < 273.0)
      Cloud_Phase = sym%SUPERCOOLED_PHASE
  endwhere

  !------------------------------------------------------------------
  ! Determine Type
  !------------------------------------------------------------------
  Line_Loop_Type: do Line_Idx = Line_Start, Line_End
  Element_Loop_Type: do Elem_Idx = 1, Num_Elem

     !--- save indices
     Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
     Vza_Idx = Vza_Idx_Rtm(Elem_Idx,Line_Idx)

     !--- skip bad pixels
     if (Invalid_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     Cloud_Type(Elem_Idx,Line_Idx) = sym%UNKNOWN_TYPE

     !-- check if indices are valid
     if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0) cycle

     if (Cloud_Mask(Elem_Idx,Line_Idx) == sym%CLEAR) then
          Cloud_Type(Elem_Idx,Line_Idx) = sym%CLEAR_TYPE
          cycle
     endif

     if (Cloud_Mask(Elem_Idx,Line_Idx) == sym%PROB_CLEAR) then
          Cloud_Type(Elem_Idx,Line_Idx) = sym%PROB_CLEAR_TYPE
          cycle
     endif

     !--- Look for smoke and dust if VIIRS 
     ! (function is in naive_bayesian_cloud_mask_module.f90)
     ! only for water phase clouds
     ! only confidently cloudy can come this far
     ! read mask bits and set type accordingly (0 or 1)
     if ( Viirs_Flag == sym%YES &
      .or. Iff_Viirs_Flag == sym%YES ) then

        ! check packed pixels if smoke & dust bits are set 
        smoke_pixel = ibits(Cld_Test_Vector_Packed(2,Elem_Idx,Line_Idx),4,1)
        if (smoke_pixel == 1) Cloud_Type(Elem_Idx,Line_Idx) = sym%SMOKE_TYPE

        dust_pixel = ibits(Cld_Test_Vector_Packed(2,Elem_Idx,Line_Idx),5,1)
        if (dust_pixel == 1) Cloud_Type(Elem_Idx,Line_Idx) = sym%DUST_TYPE

        ! cycle if pixel is dust o smoke
        if (smoke_pixel == 1 .or. dust_pixel == 1) cycle
     endif ! sensor check

     !--- supercooled
     if (Cloud_Phase(Elem_Idx,Line_Idx) == sym%SUPERCOOLED_PHASE) then
         Cloud_Type(Elem_Idx,Line_Idx) = sym%SUPERCOOLED_TYPE
         cycle
     endif

     !--- fog or water
     if (Cloud_Phase(Elem_Idx,Line_Idx) == sym%WATER_PHASE) then
      Cloud_Type(Elem_Idx,Line_Idx) = sym%WATER_TYPE
      if (Zcld_Opa(Elem_Idx,Line_Idx) < 1.0) then
          Cloud_Type(Elem_Idx,Line_Idx) = sym%FOG_TYPE
      endif
      cycle
     endif

     !---   ice types
     if (Cloud_Phase(Elem_Idx,Line_Idx) == sym%ICE_PHASE) then
        
         !--non-cirrus tests based on opacity
         if (Emiss_Window_Tropo(Elem_Idx,Line_Idx) > 0.95) then
              Cloud_Type(Elem_Idx,Line_Idx) = sym%OVERSHOOTING_TYPE
         else
              Cloud_Type(Elem_Idx,Line_Idx) = sym%OPAQUE_ICE_TYPE
         endif

         !--- if very cold, do not allow cirrus detection
         if (Bt_11um(Elem_Idx,Line_Idx) < 233.0) then
              Cloud_Type(Elem_Idx,Line_Idx) = sym%OPAQUE_ICE_TYPE
              cycle
         endif

         if (Emiss_Window_Tropo(Elem_Idx,Line_Idx) < 0.80) then

           Cloud_Type(Elem_Idx,Line_Idx) = sym%CIRRUS_TYPE

           !--- allow overlap where both 11/12 shows spectral inconsistency with single layer ice
           if ((Channel_On_Flag(Chan_Idx_11um) == sym%YES) .and. (Channel_On_Flag(Chan_Idx_12um) == sym%YES)) then
              if (Beta_11um_12um_Tropo(Elem_Idx,Line_Idx) < Beta_11um_12um_Overlap_Thresh) then
                Cloud_Type(Elem_Idx,Line_Idx) = sym%OVERLAP_TYPE
              endif
           endif

           !--- allow overlap where both 11/13.3 shows spectral inconsistency with single layer ice
           if ((Channel_On_Flag(Chan_Idx_11um) == sym%YES) .and. (Channel_On_Flag(Chan_Idx_133um) == sym%YES)) then
              if (Beta_11um_133um_Tropo(Elem_Idx,Line_Idx) < Beta_11um_133um_Overlap_Thresh) then
                Cloud_Type(Elem_Idx,Line_Idx) = sym%OVERLAP_TYPE
              endif
           endif

           !--- allow overlap where both cirrus and water clouds detected
           if (Cirrus_Flag(Elem_Idx,Line_Idx) == sym%YES .and. Water_Flag(Elem_Idx,Line_Idx) == sym%YES) then
                Cloud_Type(Elem_Idx,Line_Idx) = sym%OVERLAP_TYPE
           endif

         endif
     endif

   END DO Element_Loop_Type
END DO Line_Loop_Type

!-----------------------------------------------------------------------------
! deallocate allocated arrays
!-----------------------------------------------------------------------------
deallocate(Cirrus_Flag)
deallocate(Water_Flag)
deallocate(Overlap_Flag)
deallocate(Ice_Probability)

!-----------------------------------------------------------------------------
! nullify pointers
!-----------------------------------------------------------------------------
Invalid_Pixel_Mask => null()
Vza_Idx_Rtm => null()
Ref_16um => null()
Ref_375um => null()
Rad_11um => null()
Rad_11um_Clear => null()
Bt_11um => null()
Bt_11um_Mean => null()
Bt_11um_Min => null()
Bt_11um_Max => null()
Bt_11um_Stddev => null()
Bt_H2O => null()
Bt_H2O_Stddev => null()
Rad_H2O => null()
Rad_H2O_Clear => null()
Bt_H2O_Clear => null()
Covar_H2O_Window => null()
Tcld_H2O => null()
Pcld_H2O => null()
Zcld_H2O => null()
Zclr_H2O => null()
Tcld_Opa => null()
Pcld_Opa => null()
Zcld_Opa => null()
Cloud_Phase => null()
Cloud_Type => null()
Cloud_Mask => null()
Emiss_Window_Tropo => null()
Emiss_H2O_Tropo => null()
Beta_11um_85um_Tropo => null()
Beta_11um_12um_Tropo => null()
Beta_11um_133um_Tropo => null()
Channel_On_Flag => null()
T_Prof => null()
Z_Prof => null()
Rad_11um_BB_Profile => null()
Rad_H2O_BB_Profile => null()
deallocate(P_Prof)

end subroutine UNIVERSAL_CLOUD_TYPE
!====================================================================
! Function Name: Compute_Ice_Probability_Based_On_Temperature
!
! Function: Provide the probability that this pixel is ice 
!
! Description: 
!    Use the cloud temperature and an assumed relationship to 
!    determine the probability that the cloud is ice phase
!   
! Dependencies: 
!
! Restrictions: 
!
! Reference: 
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
!====================================================================
function COMPUTE_ICE_PROBABILITY_BASED_ON_TEMPERATURE(T_Opa) result(Ice_Prob)

real:: T_opa
real:: Ice_Prob

Ice_Prob = 1.0 - (T_opa-Ice_Temperature_Min)/(Ice_Temperature_Max - Ice_Temperature_Min)
Ice_Prob = min(1.0,Ice_Prob)
Ice_Prob = max(0.0,Ice_Prob)

end function COMPUTE_ICE_PROBABILITY_BASED_ON_TEMPERATURE

!====================================================================
! Function Name: OPAQUE_CLOUD_HEIGHT_LOCAL
!
! Function: estimate the cloud temperature/height/pressure
!
! Description: Use the 11um obs and assume the cloud is back and 
!           estimate height from 11 um BB cloud profile
!              
! Dependencies: 
!
! Restrictions: 
!
! Reference: 
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
!====================================================================
 subroutine OPAQUE_CLOUD_HEIGHT_LOCAL(Ch31_Rad, &
                                      Tropo_Level,  &
                                      Sfc_Level, &
                                      Tc_Opa, &
                                      Zc_Opa, &
                                      Pc_Opa)
   real, intent(in):: Ch31_Rad
   integer (kind=int1), intent(in):: Tropo_Level
   integer (kind=int1), intent(in):: Sfc_Level
   real, intent(out):: Pc_Opa
   real, intent(out):: Zc_Opa
   real, intent(out):: Tc_Opa

   integer:: Level_Idx
   integer:: Level_Idx_Start
   integer:: Level_Idx_End
   integer:: Level_Idx_Max_Valid_Cloud

   !--- initialize
   Pc_Opa =  Missing_Value_Real4
   Zc_Opa =  Missing_Value_Real4
   Tc_Opa =  Missing_Value_Real4

   !--- restrict levels to consider
   Level_Idx_Start = Tropo_Level 
   Level_Idx_End = Sfc_Level 

   !--- initialize levels
   Level_Idx_Max_Valid_Cloud = 0

   !--- check for stratospheric
   if (Ch31_Rad < Rad_11um_BB_Profile(Level_Idx_start)) then

        Level_Idx_Max_Valid_Cloud = Level_Idx_Start

   else

       level_loop: do Level_Idx = Level_Idx_Start, Level_Idx_End

         if (Ch31_Rad > Rad_11um_BB_Profile(Level_Idx)) then

           Level_Idx_Max_Valid_Cloud = Level_Idx

         endif

        end do level_loop

   endif

   !--- compute lowest pressure level with valid 11 micron emissivity
   if (Level_Idx_max_valid_Cloud > 0) then
        if (Level_Idx_Max_Valid_Cloud > 1) then
         Pc_Opa = P_prof(Level_Idx_Max_Valid_Cloud)
         Zc_Opa = Z_prof(Level_Idx_Max_Valid_Cloud)
         Tc_Opa = T_prof(Level_Idx_Max_Valid_Cloud)
        endif
   endif

 end subroutine OPAQUE_CLOUD_HEIGHT_LOCAL

!====================================================================
! Function Name: H2O_CLOUD_HEIGHT_LOCAL
!
! Function: estimate the cloud temperature/height/pressure
!
! Description: Use the 11um and 6.7um obs and the RTM cloud BB profiles
!              to perform h2o intercept on a pixel level. Filters
!              restrict this to high clouds only
!              
! Dependencies: 
!
! Restrictions: 
!
! Reference: 
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
!====================================================================
subroutine  H2O_CLOUD_HEIGHT_LOCAL(Rad_11um, &
                                   Rad_11um_Clear, &
                                   Rad_H2O, &
                                   Rad_H2O_Clear,  &
                                   Covar_H2O_Window, &
                                   Tropo_Level, &
                                   Sfc_Level, &
                                   Tc,  &
                                   Zc,  &
                                   Pc)

  real, intent(in):: Rad_11um
  real, intent(in):: Rad_H2O
  real, intent(in):: Rad_11um_Clear
  real, intent(in):: Rad_H2O_Clear
  real, intent(in):: Covar_H2O_Window
  integer (kind=int1), intent(in):: Sfc_Level
  integer (kind=int1), intent(in):: Tropo_Level
  real (kind=real4), intent(out) :: Tc
  real (kind=real4), intent(out) :: Pc
  real (kind=real4), intent(out) :: Zc

  real:: Rad_H2O_BB_Prediction
  real:: Slope
  real:: Intercept
  real:: Denominator
  integer:: ilev
  integer:: ilev_h2o

  !--- initialize
  Pc = Missing_Value_Real4
  Tc = Missing_Value_Real4
  Zc = Missing_Value_Real4

  !--- determine if a solution should be attempted
  if (Rad_11um_Clear - Rad_11um < Rad_11um_Thresh) then
      return 
  endif

  if (Rad_H2O_Clear - Rad_H2O < Rad_H2O_Thresh) then
      return
  endif

  if (Covar_H2O_Window /= Missing_Value_Real4 .and. Covar_H2O_Window < Bt_Ch27_Ch31_Covar_Cirrus_Thresh) then
      return 
  endif

 !--- attempt a solution

 !--- colder than tropo
 if (Rad_11um < Rad_11um_BB_Profile(Tropo_Level)) then

     ilev_h2o = Tropo_Level

 else   !if not, attempt solution

     !--- determine linear regress of h2o (y)  as a function of window (x)
      Denominator =  Rad_11um - Rad_11um_Clear

      if (Denominator < 0.0) then
             Slope = (Rad_H2O - Rad_H2O_Clear) / (Denominator)
             Intercept = Rad_H2O - Slope*Rad_11um
      else
            return 
      endif

      !--- brute force solution
      ilev_h2o = 0

      do ilev = Tropo_Level+1, Sfc_Level
          Rad_H2O_BB_Prediction = Slope*Rad_11um_BB_Profile(ilev) + Intercept

          if (Rad_H2O_BB_Prediction < 0) cycle

          if ((Rad_H2O_BB_Prediction > Rad_H2O_BB_Profile(ilev-1)) .and. & 
               (Rad_H2O_BB_Prediction <= Rad_H2O_BB_Profile(ilev))) then
               ilev_h2o = ilev
               exit
          endif

      enddo

 endif    !tropopause check

 !--- adjust back to full Rtm profile indices
 if (ilev_h2o > 0) then
       Pc = P_Prof(ilev_h2o)
       Tc = T_Prof(ilev_h2o)
       Zc = Z_Prof(ilev_h2o)
  endif

end subroutine H2O_CLOUD_HEIGHT_LOCAL
!----------------------------------------------------------------------
! Determine Clear-sky Water Vapor Weighting Function Peak Height
!
! Accomplish this by matching the compute 6.7 micron brightness
! temperature to the temperature profile
!----------------------------------------------------------------------
subroutine H2O_CLEAR_SKY_PEAK_HEIGHT(Bt_67um, &
                                    Tropo_Level,  &
                                    Sfc_Level, &
                                    Z_67um)

   real, intent(in):: Bt_67um
   integer (kind=int1), intent(in):: Tropo_Level
   integer (kind=int1), intent(in):: Sfc_Level
   real, intent(out):: Z_67um
   integer:: ilev


   !--- check for out-troposphere conditions
   if (Bt_67um == Missing_Value_Real4) then
          Z_67um = Missing_Value_Real4
          return
   endif
   if (Bt_67um < T_Prof(Tropo_Level)) then
          Z_67um = Z_Prof(Tropo_Level)
          return
   endif
   if (Bt_67um > T_Prof(Sfc_Level)) then
          Z_67um = Z_Prof(Sfc_Level)
          return
   endif

   do ilev = Tropo_Level+1, Sfc_Level
      if (Bt_67um <= T_Prof(ilev)) then
          Z_67um = Z_Prof(ilev)
          exit
      endif
   enddo

end subroutine H2O_CLEAR_SKY_PEAK_HEIGHT
!====================================================================
! Function Name: Covariance_LOCAL
!
! Function:
!    Compute the Covariance for two mxn arrays
!
! Description: Covariance = E(XY) - E(X)*E(Y)
!   
! Calling Sequence: BT_WV_BT_Window_Covar(Elem_Idx,Line_Idx) = Covariance( &
!                       sat%bt10(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                       sat%bt14(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                      Array_Width, Array_Hgt)
!   
!
! Inputs:
!   Array 1 - the first array (X)
!   Array 2 - the second array (Y)
!   Elem_size
!   Line_size
!
! Outputs: 
!   Covariance of X and Y
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Standard definition for the Covariance Computation
!
!====================================================================
FUNCTION Covariance_Local &
        (Array_One,Array_Two,Array_Width,Array_Hght,Invalid_Data_Mask) &
         RESULT(Covar_Array_One_Array_Two)
   real(kind=real4), intent(in), dimension(:,:):: Array_One
   real(kind=real4), intent(in), dimension(:,:):: Array_Two
   integer(kind=INT4), intent(in):: Array_Width
   integer(kind=INT4), intent(in):: Array_Hght
   integer(kind=INT1), intent(in), dimension(:,:):: Invalid_Data_Mask

   real(kind=real8):: Mean_Array_One
   real(kind=real8):: Mean_Array_Two
   real(kind=real8):: Mean_Array_One_x_Array_Two
   real(kind=real8):: Sum_Array_One
   real(kind=real8):: Sum_Array_Two
   real(kind=real8):: Sum_Array_One_x_Array_Two
   real(kind=real4):: Covar_Array_One_Array_Two

   !--- skip computation for pixel arrays with any missing data
   if (sum(Invalid_Data_Mask) > 0) then
      Covar_Array_One_Array_Two = Missing_Value_Real4
      return
   endif

   Sum_Array_One = sum(Array_One)
   Sum_Array_Two = sum(Array_Two)

   Mean_Array_One = Sum_Array_One / (Array_Width*Array_Hght)
   Mean_Array_Two = Sum_Array_Two / (Array_Width*Array_Hght)

   Sum_Array_One_x_Array_Two = sum(Array_One*Array_Two)
   Mean_Array_One_x_Array_Two = Sum_Array_One_x_Array_Two / (Array_Width*Array_Hght)
  
   Covar_Array_One_Array_Two  = Mean_Array_One_x_Array_Two - &
                                Mean_Array_One * Mean_Array_Two
  
 END FUNCTION Covariance_Local
!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module UNIVERSAL_CLOUD_TYPE_MODULE
