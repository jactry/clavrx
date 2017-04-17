!$Id: acha_clavrx_bridge_mod.f90 580 2014-10-08 03:38:52Z heidinger $
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Base Algorithm (ACBA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  the CLAVR-x processing system and the ACBA code.
!
! cld_height_acha    Zc_Acha
! cld_opd_dcomp      Tau_Dcomp
! cloud_water_path   Cwp
! surface_elevation
! ccl_nwp
! cloud_type
! land_class
! solar_zenith_angle
! latitude
! Fill_Value = -999.0
! QF_Fill = 1 
!
!------------------------------------------------------------------------------
module CLOUD_BASE_CLAVRX_BRIDGE

 use CLOUD_BASE
 use CLOUD_BASE_SERVICES
 use PIXEL_COMMON
   
 implicit none

 public :: CLOUD_BASE_BRIDGE

 private:: SET_DIAG, NULL_DIAG

 !--------------------------------------------------------------------
 ! define structures that will be arguments 
 !--------------------------------------------------------------------
 type(Symbol_acha), private :: Symbol
 type(acha_input_struct), private :: Input
 type(acha_output_struct), private :: Output
 type(acha_diag_struct), private :: Diag

 contains

!----------------------------------------------------------------------
! BRIDGE SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine CLOUD_BASE_BRIDGE()
 
   implicit none

   !---null pointers before filling them
   call NULL_INPUT_POINTERS()
   call NULL_OUTPUT_POINTERS()

   !-----------------------------------------------------------------------
   !---  CLAVR-x Bridge Section
   !-----------------------------------------------------------------------
   !--- initialize Input structure pointers

   !--- store integer values
   call SET_INPUT()

   !---- initalize Output structure
   call SET_OUTPUT()
  
   !----set Symbols to local values
   call SET_SYMBOL()

   !---- initialize diagnostic structure
   call SET_DIAG()

   !-----------------------------------------------------------------------
   !--- Call algorithm to make cloud geometrical boundaries
   !-----------------------------------------------------------------------
   call CLOUD_BASE_ALGORITHM(Input, Symbol, Output)
   !call CLOUD_BASE_ALGORITHM(Input, Symbol, Output, Diag)

   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT_POINTERS()
   call NULL_OUTPUT_POINTERS()

 end subroutine CLOUD_BASE_BRIDGE

 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding input
 !-----------------------------------------------------------------------------
 subroutine NULL_INPUT_POINTERS()
     Input%Invalid_Data_Mask =>  null()
     Input%Elem_Idx_Nwp =>   null()
     Input%Line_Idx_Nwp =>  null()
     Input%Elem_Idx_Opposite_Corner_NWP =>  null()
     Input%Line_Idx_Opposite_Corner_NWP =>  null()
     Input%Longitude_Interp_Weight_NWP =>  null()
     Input%Latitude_Interp_Weight_NWP =>  null()
     Input%Viewing_Zenith_Angle_Idx_Rtm =>  null()
     Input%Cosine_Zenith_Angle =>  null()
     Input%Sensor_Zenith_Angle =>  null()
     Input%Latitude =>  null()
     Input%Longitude =>  null()
     Input%Surface_Type =>  null()
     Input%Surface_Elevation =>  null()
     Input%Cloud_Mask =>  null()
     Input%Cloud_Probability => null()
     Input%Cloud_Type =>  null()
     Input%Elem_Idx_LRC_Input =>  null()
     Input%Line_Idx_LRC_Input =>   null()
     Input%Latitude_Pc =>  null()
     Input%Longitude_Pc =>  null()
     Input%Tc =>  null()
     Input%Ec =>  null()
     Input%Beta =>  null()
     Input%Pc =>  null()
     Input%Zc =>  null()
     Input%Tau =>  null()
     Input%Reff =>  null()
     Input%Tc_Uncertainty =>  null()
     Input%Ec_Uncertainty =>  null()
     Input%Beta_Uncertainty =>  null()
     Input%Pc_Uncertainty =>  null()
     Input%Zc_Uncertainty =>  null()
     Input%Lower_Cloud_Pressure =>  null()
     Input%Lower_Cloud_Temperature =>  null()
     Input%Lower_Cloud_Height =>  null()
     Input%Cdnc =>  null()
     Input%Hcld =>  null()
     Input%LCL =>  null()
     Input%CCL =>  null()
     Input%CWP => null()
     Input%CWP_nwp => null()
 end subroutine NULL_INPUT_POINTERS
 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding diagnostic
 !-----------------------------------------------------------------------------
 subroutine NULL_DIAG()
     Diag%Array_1 =>  null()
     Diag%Array_2 =>  null()
     Diag%Array_3 =>  null()
 end subroutine NULL_DIAG
 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding output
 !-----------------------------------------------------------------------------
 subroutine NULL_OUTPUT_POINTERS()
     Output%Zc_Top =>  null()
     Output%Zc_Base =>  null()
     Output%Pc_Top =>  null()
     Output%Pc_Base =>  null()
     Output%Zc_Base_Qf => null()
 end subroutine NULL_OUTPUT_POINTERS
 !-----------------------------------------------------------------------------
 ! Copy needed Symbol elements
 !-----------------------------------------------------------------------------
 subroutine SET_SYMBOL()
   Symbol%CLOUDY = sym%CLOUDY
   Symbol%PROB_CLOUDY = sym%PROB_CLOUDY
   Symbol%PROB_CLEAR = sym%PROB_CLEAR
   Symbol%CLEAR = sym%CLEAR

   Symbol%NO = sym%NO
   Symbol%YES = sym%YES

   Symbol%WATER_SFC = sym%WATER_SFC
   Symbol%EVERGREEN_NEEDLE_SFC = sym%EVERGREEN_NEEDLE_SFC
   Symbol%EVERGREEN_BROAD_SFC = sym%EVERGREEN_BROAD_SFC
   Symbol%DECIDUOUS_NEEDLE_SFC = sym%DECIDUOUS_NEEDLE_SFC
   Symbol%DECIDUOUS_BROAD_SFC = sym%DECIDUOUS_BROAD_SFC
   Symbol%MIXED_FORESTS_SFC = sym%MIXED_FORESTS_SFC
   Symbol%WOODLANDS_SFC = sym%WOODLANDS_SFC
   Symbol%WOODED_GRASS_SFC = sym%WOODED_GRASS_SFC
   Symbol%CLOSED_SHRUBS_SFC = sym%CLOSED_SHRUBS_SFC
   Symbol%OPEN_SHRUBS_SFC = sym%OPEN_SHRUBS_SFC
   Symbol%GRASSES_SFC = sym%GRASSES_SFC
   Symbol%CROPLANDS_SFC = sym%CROPLANDS_SFC
   Symbol%BARE_SFC = sym%BARE_SFC
   Symbol%URBAN_SFC = sym%URBAN_SFC

   Symbol%SHALLOW_OCEAN = sym%SHALLOW_OCEAN
   Symbol%LAND = sym%LAND
   Symbol%COASTLINE = sym%COASTLINE
   Symbol%SHALLOW_INLAND_WATER = sym%SHALLOW_INLAND_WATER
   Symbol%EPHEMERAL_WATER = sym%EPHEMERAL_WATER
   Symbol%DEEP_INLAND_WATER = sym%DEEP_INLAND_WATER
   Symbol%MODERATE_OCEAN = sym%MODERATE_OCEAN
   Symbol%DEEP_OCEAN = sym%DEEP_OCEAN

   Symbol%NO_SNOW = sym%NO_SNOW
   Symbol%SEA_ICE = sym%SEA_ICE
   Symbol%SNOW = sym%SNOW

   Symbol%CLEAR_TYPE = sym%CLEAR_TYPE
   Symbol%PROB_CLEAR_TYPE = sym%PROB_CLEAR_TYPE
   Symbol%FOG_TYPE = sym%FOG_TYPE
   Symbol%WATER_TYPE = sym%WATER_TYPE
   Symbol%SUPERCOOLED_TYPE = sym%SUPERCOOLED_TYPE
   Symbol%MIXED_TYPE = sym%MIXED_TYPE
   Symbol%OPAQUE_ICE_TYPE = sym%OPAQUE_ICE_TYPE
   Symbol%TICE_TYPE = sym%TICE_TYPE
   Symbol%CIRRUS_TYPE = sym%CIRRUS_TYPE
   Symbol%OVERLAP_TYPE = sym%OVERLAP_TYPE
   Symbol%OVERSHOOTING_TYPE = sym%OVERSHOOTING_TYPE
   Symbol%UNKNOWN_TYPE = sym%UNKNOWN_TYPE
   Symbol%DUST_TYPE = sym%DUST_TYPE
   Symbol%SMOKE_TYPE = sym%SMOKE_TYPE
   Symbol%FIRE_TYPE = sym%FIRE_TYPE

   Symbol%CLEAR_PHASE = sym%CLEAR_PHASE
   Symbol%WATER_PHASE = sym%WATER_PHASE
   Symbol%SUPERCOOLED_PHASE = sym%SUPERCOOLED_PHASE
   Symbol%MIXED_PHASE = sym%MIXED_PHASE
   Symbol%ICE_PHASE = sym%ICE_PHASE
   Symbol%UNKNOWN_PHASE = sym%UNKNOWN_PHASE
 end subroutine SET_SYMBOL

 subroutine SET_OUTPUT()
   Output%Zc_Top => ACHA%Zc_Top
   Output%Zc_Base => ACHA%Zc_Base
   Output%Pc_Top => ACHA%Pc_Top
   Output%Pc_Base => ACHA%Pc_Base
   Output%Zc_Base_Qf => ACHA%base_Quality_Flag
 end subroutine SET_OUTPUT

 subroutine SET_INPUT()

   Input%Number_of_Elements = Image%Number_Of_Elements
   Input%Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment
   Input%Num_Line_Max = Image%Number_Of_Lines_Per_Segment

   Input%Chan_Idx_67um = 27     !channel number for 6.7
   Input%Chan_Idx_85um = 29     !channel number for 8.5
   Input%Chan_Idx_11um = 31     !channel number for 11
   Input%Chan_Idx_12um = 32     !channel number for 12
   Input%Chan_Idx_133um = 33  !channel number for 13.3

   Input%Chan_On_67um = Sensor%Chan_On_Flag_Default(27)
   Input%Chan_On_85um = Sensor%Chan_On_Flag_Default(29)
   Input%Chan_On_11um = Sensor%Chan_On_Flag_Default(31)
   Input%Chan_On_12um = Sensor%Chan_On_Flag_Default(32)
   Input%Chan_On_133um = Sensor%Chan_On_Flag_Default(33)

   Input%Invalid_Data_Mask => Bad_Pixel_Mask
   Input%Elem_Idx_Nwp =>  I_Nwp
   Input%Line_Idx_Nwp => J_Nwp
   Input%Elem_Idx_Opposite_Corner_NWP => I_Nwp_x
   Input%Line_Idx_Opposite_Corner_NWP => J_Nwp_x
   Input%Longitude_Interp_Weight_NWP => Lon_Nwp_Fac
   Input%Latitude_Interp_Weight_NWP => Lat_Nwp_Fac
   Input%Viewing_Zenith_Angle_Idx_Rtm => Zen_Idx_Rtm

   Input%Cosine_Zenith_Angle => Geo%Coszen
   Input%Sensor_Zenith_Angle => Geo%Satzen
   Input%Latitude => Nav%Lat
   Input%Longitude => Nav%Lon

   Input%Surface_Type => Sfc%Sfc_Type

   Input%Surface_Elevation => Sfc%Zsfc
   Input%Cloud_Mask => CLDMASK%Cld_Mask
   Input%Cloud_Type => Cld_Type
   Input%Cloud_Probability => CLDMASK%Posterior_Cld_Probability

   Input%Elem_Idx_LRC_Input => I_LRC
   Input%Line_Idx_LRC_Input =>  J_LRC

   Input%Latitude_Pc => Nav%Lat_Pc
   Input%Longitude_Pc => Nav%Lon_Pc
   Input%Tc => ACHA%Tc
   Input%Ec => ACHA%Ec
   Input%Beta => ACHA%Beta
   Input%Pc => ACHA%Pc
   Input%Zc => ACHA%Zc
   Input%Tau => ACHA%Tau
   Input%Reff => ACHA%Reff
   Input%Tc_Uncertainty => ACHA%Tc_Uncertainty
   Input%Ec_Uncertainty => ACHA%Ec_Uncertainty
   Input%Beta_Uncertainty => ACHA%Beta_Uncertainty
   Input%Pc_Uncertainty => ACHA%Pc_Uncertainty
   Input%Zc_Uncertainty => ACHA%Zc_Uncertainty
   Input%Lower_Cloud_Pressure => ACHA%Lower_Pc
   Input%Lower_Cloud_Temperature => ACHA%Lower_Tc
   Input%Lower_Cloud_Height => ACHA%Lower_Zc
   Input%Cdnc => Cdnc_DCOMP
   Input%Hcld => Hcld_DCOMP
   Input%LCL => LCL_Height_Nwp_Pix
   Input%CCL => CCL_Height_Nwp_Pix
   Input%CWP => Cwp_Dcomp
   Input%CWP_nwp => Cwp_Nwp_Pix
 end subroutine SET_INPUT
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
 subroutine SET_DIAG
     Diag%Array_1 => Diag_Pix_Array_1 
     Diag%Array_2 => Diag_Pix_Array_2 
     Diag%Array_3 => Diag_Pix_Array_3 
 end subroutine SET_DIAG

end module CLOUD_BASE_CLAVRX_BRIDGE
