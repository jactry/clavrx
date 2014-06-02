!$Id$
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Height Algorithm (ACHA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  processing systems and the ACHA code.
!
!------------------------------------------------------------------------------
module ACHA_CLAVRX_BRIDGE_MOD

 use AWG_CLOUD_HEIGHT
 use ACHA_SERVICES_MOD
   
 implicit none

 public :: AWG_CLOUD_HEIGHT_BRIDGE

 contains


!----------------------------------------------------------------------
! BEGINNING OF ACHA SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine AWG_CLOUD_HEIGHT_BRIDGE()
 
   implicit none

   !--------------------------------------------------------------------
   ! define structures that will be arguments to ACHA
   !--------------------------------------------------------------------
   type(symbol_acha) :: symbol
   type(acha_input_struct) :: Acha_Input
   type(acha_output_struct) :: Acha_Output

   !---null pointers before filling them
   call NULL_ACHA_POINTERS(Acha_Input, Acha_Output)

   !-----------------------------------------------------------------------
   !---  CLAVR-x Bridge Section
   !-----------------------------------------------------------------------
   !--- initialize Input structure pointers

   !--- store integer values
   Acha_Input%Number_of_Elements = Num_Pix
   Acha_Input%Number_of_Lines = Num_Scans_Read
   Acha_Input%Number_of_Lines = Num_Scans_Per_Segment
   Acha_Input%Num_Line_Max = Num_Scans_Per_Segment
   Acha_Input%Process_Undetected_Cloud_Flag = Process_Undetected_Cloud_Flag
   Acha_Input%Smooth_Nwp_Fields_Flag = Smooth_Nwp_Flag
   Acha_Input%ACHA_Mode_Flag_In = ACHA_MODE

   Acha_Input%Chan_Idx_67um = 27     !channel number for 6.7
   Acha_Input%Chan_Idx_85um = 29     !channel number for 8.5
   Acha_Input%Chan_Idx_11um = 31     !channel number for 11
   Acha_Input%Chan_Idx_12um = 32     !channel number for 12
   Acha_Input%Chan_Idx_133um = 33  !channel number for 13.3

   Acha_Input%Chan_On_67um = Chan_On_Flag_Default(27)
   Acha_Input%Chan_On_85um = Chan_On_Flag_Default(29)
   Acha_Input%Chan_On_11um = Chan_On_Flag_Default(31)
   Acha_Input%Chan_On_12um = Chan_On_Flag_Default(32)
   Acha_Input%Chan_On_133um = Chan_On_Flag_Default(33)

   Acha_Input%Invalid_Data_Mask => Bad_Pixel_Mask
   Acha_Input%Elem_Idx_Nwp =>  I_Nwp
   Acha_Input%Line_Idx_Nwp => J_Nwp
   Acha_Input%Elem_Idx_Opposite_Corner_NWP => I_Nwp_x
   Acha_Input%Line_Idx_Opposite_Corner_NWP => J_Nwp_x
   Acha_Input%Longitude_Interp_Weight_NWP => Lon_Nwp_Fac
   Acha_Input%Latitude_Interp_Weight_NWP => Lat_Nwp_Fac
   Acha_Input%Viewing_Zenith_Angle_Idx_Rtm => Zen_Idx_Rtm
   Acha_Input%Bt_67um => ch(27)%Bt_Toa
   Acha_Input%Bt_85um => ch(29)%Bt_Toa
   Acha_Input%Bt_11um => ch(31)%Bt_Toa
   Acha_Input%Bt_12um => ch(32)%Bt_Toa
   Acha_Input%Bt_133um => ch(33)%Bt_Toa
   Acha_Input%Covar_Bt_11um_67um => Covar_Ch27_Ch31_5x5

   Acha_Input%Rad_67um => ch(27)%Rad_Toa
   Acha_Input%Rad_11um => ch(31)%Rad_Toa
   Acha_Input%Cosine_Zenith_Angle => Coszen
   Acha_Input%Sensor_Zenith_Angle => Satzen
   Acha_Input%Sensor_Azimuth_Angle => Sataz
   Acha_Input%Latitude => Lat
   Acha_Input%Longitude => Lon

   Acha_Input%Snow_Class => Snow
   Acha_Input%Surface_Type => Sfc_Type

   Acha_Input%Surface_Temperature =>Tsfc_Nwp_Pix
   Acha_Input%Surface_Air_Temperature => Tair_Nwp_Pix
   Acha_Input%Tropopause_Temperature => Ttropo_Nwp_Pix
   Acha_Input%Surface_Pressure => Psfc_Nwp_Pix

   Acha_Input%Surface_Elevation => Zsfc
   Acha_Input%Cloud_Mask => Cld_Mask
   Acha_Input%Cloud_Type => Cld_Type

   Acha_Input%Rad_Clear_67um => ch(27)%Rad_Toa_Clear
   Acha_Input%Rad_Clear_85um => ch(29)%Rad_Toa_Clear
   Acha_Input%Rad_Clear_11um => ch(31)%Rad_Toa_Clear
   Acha_Input%Rad_Clear_12um => ch(32)%Rad_Toa_Clear
   Acha_Input%Rad_Clear_133um => ch(33)%Rad_Toa_Clear
   Acha_Input%Surface_Emissivity_39um => ch(20)%Sfc_Emiss

   Acha_Input%Elem_Idx_LRC_Input => I_LRC
   Acha_Input%Line_Idx_LRC_Input =>  J_LRC
     
   !---- initalize Output structure
   Acha_Output%Latitude_Pc => Lat_Pc
   Acha_Output%Longitude_Pc => Lon_Pc
   Acha_Output%Tc => Tc_Acha
   Acha_Output%Ec => Ec_Acha
   Acha_Output%Beta => Beta_Acha
   Acha_Output%Pc => Pc_Acha
   Acha_Output%Zc => Zc_Acha
   Acha_Output%Tau => Tau_Acha
   Acha_Output%Reff => Reff_Acha
   Acha_Output%Tc_Uncertainty => Tc_Acha_Uncertainty
   Acha_Output%Ec_Uncertainty => Ec_Acha_Uncertainty
   Acha_Output%Beta_Uncertainty => Beta_Acha_Uncertainty
   Acha_Output%Pc_Uncertainty => Pc_Acha_Uncertainty
   Acha_Output%Zc_Uncertainty => Zc_Acha_Uncertainty
   Acha_Output%Cloud_Layer => Cld_Layer_Acha
   Acha_Output%Lower_Cloud_Pressure => Pc_Lower_Cloud
   Acha_Output%Lower_Cloud_Temperature => Tc_Lower_Cloud
   Acha_Output%Lower_Cloud_Height => Zc_Lower_Cloud
   Acha_Output%Zc_Top => Zc_Top_Acha
   Acha_Output%Zc_Base => Zc_Base_Acha
   Acha_Output%Qf => Acha_Quality_Flag
   Acha_Output%OE_Qf => Acha_OE_Quality_Flags
   Acha_Output%Packed_Qf => Acha_Packed_Quality_Flags
   Acha_Output%Packed_Meta_Data => Acha_Packed_Meta_Data_Flags
   Acha_Output%Processing_Order  => Acha_Processing_Order_Global
  

   !----set symbols to local values
   symbol%CLOUDY = sym%CLOUDY
   symbol%PROB_CLOUDY = sym%PROB_CLOUDY
   symbol%PROB_CLEAR = sym%PROB_CLEAR
   symbol%CLEAR = sym%CLEAR

   symbol%NO = sym%NO
   symbol%YES = sym%YES

   symbol%WATER_SFC = sym%WATER_SFC
   symbol%EVERGREEN_NEEDLE_SFC = sym%EVERGREEN_NEEDLE_SFC
   symbol%EVERGREEN_BROAD_SFC = sym%EVERGREEN_BROAD_SFC
   symbol%DECIDUOUS_NEEDLE_SFC = sym%DECIDUOUS_NEEDLE_SFC
   symbol%DECIDUOUS_BROAD_SFC = sym%DECIDUOUS_BROAD_SFC
   symbol%MIXED_FORESTS_SFC = sym%MIXED_FORESTS_SFC
   symbol%WOODLANDS_SFC = sym%WOODLANDS_SFC
   symbol%WOODED_GRASS_SFC = sym%WOODED_GRASS_SFC
   symbol%CLOSED_SHRUBS_SFC = sym%CLOSED_SHRUBS_SFC
   symbol%OPEN_SHRUBS_SFC = sym%OPEN_SHRUBS_SFC
   symbol%GRASSES_SFC = sym%GRASSES_SFC
   symbol%CROPLANDS_SFC = sym%CROPLANDS_SFC
   symbol%BARE_SFC = sym%BARE_SFC
   symbol%URBAN_SFC = sym%URBAN_SFC

   symbol%SHALLOW_OCEAN = sym%SHALLOW_OCEAN
   symbol%LAND = sym%LAND
   symbol%COASTLINE = sym%COASTLINE
   symbol%SHALLOW_INLAND_WATER = sym%SHALLOW_INLAND_WATER
   symbol%EPHEMERAL_WATER = sym%EPHEMERAL_WATER
   symbol%DEEP_INLAND_WATER = sym%DEEP_INLAND_WATER
   symbol%MODERATE_OCEAN = sym%MODERATE_OCEAN
   symbol%DEEP_OCEAN = sym%DEEP_OCEAN

   symbol%NO_SNOW = sym%NO_SNOW
   symbol%SEA_ICE = sym%SEA_ICE
   symbol%SNOW = sym%SNOW

   symbol%CLEAR_TYPE = sym%CLEAR_TYPE
   symbol%PROB_CLEAR_TYPE = sym%PROB_CLEAR_TYPE
   symbol%FOG_TYPE = sym%FOG_TYPE
   symbol%WATER_TYPE = sym%WATER_TYPE
   symbol%SUPERCOOLED_TYPE = sym%SUPERCOOLED_TYPE
   symbol%MIXED_TYPE = sym%MIXED_TYPE
   symbol%OPAQUE_ICE_TYPE = sym%OPAQUE_ICE_TYPE
   symbol%TICE_TYPE = sym%TICE_TYPE
   symbol%CIRRUS_TYPE = sym%CIRRUS_TYPE
   symbol%OVERLAP_TYPE = sym%OVERLAP_TYPE
   symbol%OVERSHOOTING_TYPE = sym%OVERSHOOTING_TYPE
   symbol%UNKNOWN_TYPE = sym%UNKNOWN_TYPE
   symbol%DUST_TYPE = sym%DUST_TYPE
   symbol%SMOKE_TYPE = sym%SMOKE_TYPE
   symbol%FIRE_TYPE = sym%FIRE_TYPE

   symbol%CLEAR_PHASE = sym%CLEAR_PHASE
   symbol%WATER_PHASE = sym%WATER_PHASE
   symbol%SUPERCOOLED_PHASE = sym%SUPERCOOLED_PHASE
   symbol%MIXED_PHASE = sym%MIXED_PHASE
   symbol%ICE_PHASE = sym%ICE_PHASE
   symbol%UNKNOWN_PHASE = sym%UNKNOWN_PHASE

   !-----------------------------------------------------------------------
   !--- Call to AWG CLoud Height Algorithm (ACHA)
   !-----------------------------------------------------------------------
   call AWG_CLOUD_HEIGHT_ALGORITHM(Acha_Input, &
                                   symbol, &
                                   Acha_Output)
    
   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_ACHA_POINTERS(Acha_Input, Acha_Output)

   !-----------------------------------------------------------------------
   !---  read CVS Tag from ACHA and store in global variable for output
   !-----------------------------------------------------------------------
   call SET_ACHA_VERSION(Acha_Version)
                                   

 end subroutine AWG_CLOUD_HEIGHT_BRIDGE

 !-----------------------------------------------------------------------------
 !
 !-----------------------------------------------------------------------------
 subroutine NULL_ACHA_POINTERS(Acha_Input, Acha_Output)

     type(acha_input_struct), intent(inout) :: Acha_Input
     type(acha_output_struct), intent(inout) :: Acha_Output

     !--- null input pointers
     Acha_Input%Invalid_Data_Mask =>  NULL()
     Acha_Input%Elem_Idx_Nwp =>   NULL()
     Acha_Input%Line_Idx_Nwp =>  NULL()
     Acha_Input%Elem_Idx_Opposite_Corner_NWP =>  NULL()
     Acha_Input%Line_Idx_Opposite_Corner_NWP =>  NULL()
     Acha_Input%Longitude_Interp_Weight_NWP =>  NULL()
     Acha_Input%Latitude_Interp_Weight_NWP =>  NULL()
     Acha_Input%Viewing_Zenith_Angle_Idx_Rtm =>  NULL()
     Acha_Input%Bt_67um =>  NULL()
     Acha_Input%Bt_85um =>  NULL()
     Acha_Input%Bt_11um =>  NULL()
     Acha_Input%Bt_12um =>  NULL()
     Acha_Input%Bt_133um =>  NULL()
     Acha_Input%Rad_67um =>  NULL()
     Acha_Input%Rad_11um =>  NULL()
     Acha_Input%Cosine_Zenith_Angle =>  NULL()
     Acha_Input%Sensor_Zenith_Angle =>  NULL()
     Acha_Input%Sensor_Azimuth_Angle =>  NULL()
     Acha_Input%Latitude =>  NULL()
     Acha_Input%Longitude =>  NULL()
     Acha_Input%Snow_Class =>  NULL()
     Acha_Input%Surface_Type =>  NULL()
     Acha_Input%Surface_Temperature => NULL()
     Acha_Input%Surface_Air_Temperature =>  NULL()
     Acha_Input%Tropopause_Temperature =>  NULL()
     Acha_Input%Surface_Pressure =>  NULL()
     Acha_Input%Surface_Elevation =>  NULL()
     Acha_Input%Cloud_Mask =>  NULL()
     Acha_Input%Cloud_Type =>  NULL()
     Acha_Input%Rad_Clear_67um =>  NULL()
     Acha_Input%Rad_Clear_85um =>  NULL()
     Acha_Input%Rad_Clear_11um =>  NULL()
     Acha_Input%Rad_Clear_12um =>  NULL()
     Acha_Input%Rad_Clear_133um =>  NULL()
     Acha_Input%Surface_Emissivity_39um =>  NULL()
     Acha_Input%Elem_Idx_LRC_Input =>  NULL()
     Acha_Input%Line_Idx_LRC_Input =>   NULL()

     !--- null output pointers
     Acha_Output%Latitude_Pc =>  NULL()
     Acha_Output%Longitude_Pc =>  NULL()
     Acha_Output%Tc =>  NULL()
     Acha_Output%Ec =>  NULL()
     Acha_Output%Beta =>  NULL()
     Acha_Output%Pc =>  NULL()
     Acha_Output%Zc =>  NULL()
     Acha_Output%Tau =>  NULL()
     Acha_Output%Reff =>  NULL()
     Acha_Output%Tc_Uncertainty =>  NULL()
     Acha_Output%Ec_Uncertainty =>  NULL()
     Acha_Output%Beta_Uncertainty =>  NULL()
     Acha_Output%Pc_Uncertainty =>  NULL()
     Acha_Output%Zc_Uncertainty =>  NULL()
     Acha_Output%Cloud_Layer =>  NULL()
     Acha_Output%Lower_Cloud_Pressure =>  NULL()
     Acha_Output%Lower_Cloud_Temperature =>  NULL()
     Acha_Output%Lower_Cloud_Height =>  NULL()
     Acha_Output%Zc_Top =>  NULL()
     Acha_Output%Zc_Base =>  NULL()
     Acha_Output%Qf =>  NULL()
     Acha_Output%OE_Qf =>  NULL()
     Acha_Output%Packed_Qf =>  NULL()
     Acha_Output%Packed_Meta_Data =>  NULL()
     Acha_Output%Processing_Order  =>  NULL()
 
 end subroutine

end module ACHA_CLAVRX_BRIDGE_MOD
