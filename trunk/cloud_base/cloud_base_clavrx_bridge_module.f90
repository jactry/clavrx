!$Id: acha_clavrx_bridge_mod.f90 580 2014-10-08 03:38:52Z heidinger $
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Base Algorithm (ACBA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  the CLAVR-x processing system and the ACBA code.
!
!------------------------------------------------------------------------------
module CLOUD_BASE_CLAVRX_BRIDGE

 use CLOUD_BASE
 use CLOUD_BASE_SERVICES
 use PIXEL_COMMON
   
 implicit none

 public :: CLOUD_BASE_BRIDGE

 !--------------------------------------------------------------------
 ! define structures that will be arguments 
 !--------------------------------------------------------------------
 type(Symbol_acha), private :: Symbol
 type(acha_input_struct), private :: Base_Input
 type(acha_output_struct), private :: Base_Output

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

   !-----------------------------------------------------------------------
   !--- Call algorithm to make ACHA optical and microphysical properties
   !-----------------------------------------------------------------------
   call CLOUD_BASE_ALGORITHM(Base_Input, &
                             Symbol, &
                             Base_Output)

   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT_POINTERS()
   call NULL_OUTPUT_POINTERS()

 end subroutine CLOUD_BASE_BRIDGE

 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding input to ACHA
 !-----------------------------------------------------------------------------
 subroutine NULL_INPUT_POINTERS()
     Base_Input%Invalid_Data_Mask =>  null()
     Base_Input%Elem_Idx_Nwp =>   null()
     Base_Input%Line_Idx_Nwp =>  null()
     Base_Input%Elem_Idx_Opposite_Corner_NWP =>  null()
     Base_Input%Line_Idx_Opposite_Corner_NWP =>  null()
     Base_Input%Longitude_Interp_Weight_NWP =>  null()
     Base_Input%Latitude_Interp_Weight_NWP =>  null()
     Base_Input%Viewing_Zenith_Angle_Idx_Rtm =>  null()
     Base_Input%Bt_67um =>  null()
     Base_Input%Bt_85um =>  null()
     Base_Input%Bt_11um =>  null()
     Base_Input%Bt_12um =>  null()
     Base_Input%Bt_133um =>  null()
     Base_Input%Rad_67um =>  null()
     Base_Input%Rad_11um =>  null()
     Base_Input%Cosine_Zenith_Angle =>  null()
     Base_Input%Sensor_Zenith_Angle =>  null()
     Base_Input%Sensor_Azimuth_Angle =>  null()
     Base_Input%Latitude =>  null()
     Base_Input%Longitude =>  null()
     Base_Input%Snow_Class =>  null()
     Base_Input%Surface_Type =>  null()
     Base_Input%Surface_Temperature => null()
     Base_Input%Surface_Air_Temperature =>  null()
     Base_Input%Tropopause_Temperature =>  null()
     Base_Input%Surface_Pressure =>  null()
     Base_Input%Surface_Elevation =>  null()
     Base_Input%Cloud_Mask =>  null()
     Base_Input%Cloud_Probability => null()
     Base_Input%Cloud_Type =>  null()
     Base_Input%Rad_Clear_67um =>  null()
     Base_Input%Rad_Clear_85um =>  null()
     Base_Input%Rad_Clear_11um =>  null()
     Base_Input%Rad_Clear_12um =>  null()
     Base_Input%Rad_Clear_133um =>  null()
     Base_Input%Surface_Emissivity_39um =>  null()
     Base_Input%Elem_Idx_LRC_Input =>  null()
     Base_Input%Line_Idx_LRC_Input =>   null()
     Base_Input%Tc_Cirrus_Sounder =>   null()
 end subroutine NULL_INPUT_POINTERS
 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding output to ACHA
 !-----------------------------------------------------------------------------
 subroutine NULL_OUTPUT_POINTERS()
     Base_Output%Latitude_Pc =>  null()
     Base_Output%Longitude_Pc =>  null()
     Base_Output%Tc =>  null()
     Base_Output%Ec =>  null()
     Base_Output%Beta =>  null()
     Base_Output%Pc =>  null()
     Base_Output%Zc =>  null()
     Base_Output%Tau =>  null()
     Base_Output%Reff =>  null()
     Base_Output%Tc_Uncertainty =>  null()
     Base_Output%Ec_Uncertainty =>  null()
     Base_Output%Beta_Uncertainty =>  null()
     Base_Output%Pc_Uncertainty =>  null()
     Base_Output%Zc_Uncertainty =>  null()
     Base_Output%Lower_Cloud_Pressure =>  null()
     Base_Output%Lower_Cloud_Temperature =>  null()
     Base_Output%Lower_Cloud_Height =>  null()
     Base_Output%Zc_Top =>  null()
     Base_Output%Zc_Base =>  null()
     Base_Output%Qf =>  null()
     Base_Output%OE_Qf =>  null()
     Base_Output%Packed_Qf =>  null()
     Base_Output%Packed_Meta_Data =>  null()
     Base_Output%Processing_Order  =>  null()
     Base_Output%Cost =>  null()
     Base_Output%Cloud_Layer =>  null()
     Base_Output%Total_Cloud_Fraction =>  null()
     Base_Output%Total_Cloud_Fraction_Uncer =>  null()
     Base_Output%High_Cloud_Fraction =>  null()
     Base_Output%Mid_Cloud_Fraction =>  null()
     Base_Output%Low_Cloud_Fraction =>  null()
 
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
   Base_Output%Latitude_Pc => Lat_Pc
   Base_Output%Longitude_Pc => Lon_Pc
   Base_Output%Tc => Tc_Acha
   Base_Output%Ec => Ec_Acha
   Base_Output%Beta => Beta_Acha
   Base_Output%Pc => Pc_Acha
   Base_Output%Zc => Zc_Acha
   Base_Output%Tau => Tau_Acha
   Base_Output%Reff => Reff_Acha
   Base_Output%Tc_Uncertainty => Tc_Acha_Uncertainty
   Base_Output%Ec_Uncertainty => Ec_Acha_Uncertainty
   Base_Output%Beta_Uncertainty => Beta_Acha_Uncertainty
   Base_Output%Pc_Uncertainty => Pc_Acha_Uncertainty
   Base_Output%Zc_Uncertainty => Zc_Acha_Uncertainty
   Base_Output%Lower_Cloud_Pressure => Pc_Lower_Cloud
   Base_Output%Lower_Cloud_Temperature => Tc_Lower_Cloud
   Base_Output%Lower_Cloud_Height => Zc_Lower_Cloud
   Base_Output%Zc_Top => Zc_Top_Acha
   Base_Output%Zc_Base => Zc_Base_Acha
   Base_Output%Qf => Acha_Quality_Flag
   Base_Output%OE_Qf => Acha_OE_Quality_Flags
   Base_Output%Packed_Qf => Acha_Packed_Quality_Flags
   Base_Output%Packed_Meta_Data => Acha_Packed_Meta_Data_Flags
   Base_Output%Processing_Order  => Acha_Processing_Order_Global
   Base_Output%Cost  => Cost_Acha
   Base_Output%Cloud_Layer  => Cld_Layer_Acha
   Base_Output%Total_Cloud_Fraction => Cloud_Fraction_3x3
   Base_Output%Total_Cloud_Fraction_Uncer => Cloud_Fraction_Uncer_3x3
   Base_Output%High_Cloud_Fraction => High_Cloud_Fraction_3x3
   Base_Output%Mid_Cloud_Fraction => Mid_Cloud_Fraction_3x3
   Base_Output%Low_Cloud_Fraction => Low_Cloud_Fraction_3x3
 end subroutine SET_OUTPUT

 subroutine SET_INPUT()
   Base_Input%Number_of_Elements = Num_Pix
   Base_Input%Number_of_Lines = Num_Scans_Read
   Base_Input%Number_of_Lines = Num_Scans_Per_Segment
   Base_Input%Num_Line_Max = Num_Scans_Per_Segment
   Base_Input%Process_Undetected_Cloud_Flag = Process_Undetected_Cloud_Flag
   Base_Input%Smooth_Nwp_Fields_Flag = Smooth_Nwp_Flag
   Base_Input%ACHA_Mode_Flag_In = ACHA_MODE
   Base_Input%Sensor_Resolution_KM = Sensor_Resolution_KM

   Base_Input%Chan_Idx_67um = 27     !channel number for 6.7
   Base_Input%Chan_Idx_85um = 29     !channel number for 8.5
   Base_Input%Chan_Idx_11um = 31     !channel number for 11
   Base_Input%Chan_Idx_12um = 32     !channel number for 12
   Base_Input%Chan_Idx_133um = 33  !channel number for 13.3

   Base_Input%Chan_On_67um = Chan_On_Flag_Default(27)
   Base_Input%Chan_On_85um = Chan_On_Flag_Default(29)
   Base_Input%Chan_On_11um = Chan_On_Flag_Default(31)
   Base_Input%Chan_On_12um = Chan_On_Flag_Default(32)
   Base_Input%Chan_On_133um = Chan_On_Flag_Default(33)

   Base_Input%Invalid_Data_Mask => Bad_Pixel_Mask
   Base_Input%Elem_Idx_Nwp =>  I_Nwp
   Base_Input%Line_Idx_Nwp => J_Nwp
   Base_Input%Elem_Idx_Opposite_Corner_NWP => I_Nwp_x
   Base_Input%Line_Idx_Opposite_Corner_NWP => J_Nwp_x
   Base_Input%Longitude_Interp_Weight_NWP => Lon_Nwp_Fac
   Base_Input%Latitude_Interp_Weight_NWP => Lat_Nwp_Fac
   Base_Input%Viewing_Zenith_Angle_Idx_Rtm => Zen_Idx_Rtm
   Base_Input%Bt_67um => ch(27)%Bt_Toa
   Base_Input%Bt_85um => ch(29)%Bt_Toa
   Base_Input%Bt_11um => ch(31)%Bt_Toa
   Base_Input%Bt_12um => ch(32)%Bt_Toa
   Base_Input%Bt_133um => ch(33)%Bt_Toa
   Base_Input%Covar_Bt_11um_67um => Covar_Ch27_Ch31_5x5

   Base_Input%Rad_67um => ch(27)%Rad_Toa
   Base_Input%Rad_11um => ch(31)%Rad_Toa
   Base_Input%Cosine_Zenith_Angle => Coszen
   Base_Input%Sensor_Zenith_Angle => Satzen
   Base_Input%Sensor_Azimuth_Angle => Sataz
   Base_Input%Latitude => Lat
   Base_Input%Longitude => Lon

   Base_Input%Snow_Class => Snow
   Base_Input%Surface_Type => Sfc_Type

   Base_Input%Surface_Temperature =>Tsfc_Nwp_Pix
   Base_Input%Surface_Air_Temperature => Tair_Nwp_Pix
   Base_Input%Tropopause_Temperature => Ttropo_Nwp_Pix
   Base_Input%Surface_Pressure => Psfc_Nwp_Pix

   Base_Input%Surface_Elevation => Zsfc
   Base_Input%Cloud_Mask => Cld_Mask
   Base_Input%Cloud_Type => Cld_Type
   Base_Input%Cloud_Probability => Posterior_Cld_Probability

   Base_Input%Rad_Clear_67um => ch(27)%Rad_Toa_Clear
   Base_Input%Rad_Clear_85um => ch(29)%Rad_Toa_Clear
   Base_Input%Rad_Clear_11um => ch(31)%Rad_Toa_Clear
   Base_Input%Rad_Clear_12um => ch(32)%Rad_Toa_Clear
   Base_Input%Rad_Clear_133um => ch(33)%Rad_Toa_Clear
   Base_Input%Surface_Emissivity_39um => ch(20)%Sfc_Emiss

   Base_Input%Elem_Idx_LRC_Input => I_LRC
   Base_Input%Line_Idx_LRC_Input =>  J_LRC
   Base_Input%Tc_Cirrus_Sounder =>  Tc_Cirrus_Co2
 end subroutine SET_INPUT

end module CLOUD_BASE_CLAVRX_BRIDGE
