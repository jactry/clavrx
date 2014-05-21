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
 public :: SET_ACHA_VERSION
   
   
 contains


 !====================================================================
 !  record cvs version as a global variable for output to hdf
 !
 !  THIS IS BROKEN- IT ONLY PASSES BRIDGE INFO - FIXME
 !
 !====================================================================
 subroutine SET_ACHA_VERSION()
   Acha_Version = "$Id$"
 end subroutine SET_ACHA_VERSION


!----------------------------------------------------------------------
! BEGINNING OF ACHA SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine AWG_CLOUD_HEIGHT_BRIDGE(Ctxt, Stat)

   implicit none
   type(AWG_CLX_CLOUD_HEIGHT_Ctxt) :: Ctxt
   integer(long) :: Stat
   integer(byte), pointer, dimension(:):: Chan_On_Flag_Default

   !--------------------------------------------------------------------
   ! define structures that will be arguments to ACHA
   !--------------------------------------------------------------------
   type(symbol_acha) :: symbol
   type(acha_input_struct) :: Acha_Input
   type(acha_output_struct) :: Acha_Output

   !---null pointers before filling them
   call NULL_ACHA_POINTERS(Acha_Input, Acha_Output)

   !-----------------------------------------------------------------------
   !---  NESDIS FRAMEWORK Bridge Section
   !-----------------------------------------------------------------------
     Acha_Input%Ctxt => Ctxt   

! initialize Input structure pointers
     Acha_Input%Number_of_Elements = Ctxt%SegmentInfo%Current_Column_Size
     Acha_Input%Number_of_Lines = Ctxt%SegmentInfo%Current_Row_Size

     Acha_Input%Num_Line_Max = Ctxt%SegmentInfo%Current_Row_Size
     Acha_Input%Process_Undetected_Cloud_Flag_Local = sym%NO
     Acha_Input%Smooth_Nwp_Fields_Flag = sym%YES


     CALL NFIP_CloudHeight_AchaMode(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Input%ACHA_Mode_Flag_In)
     

     !These need to be verified for AIT Framework.
     !Should be done via function call if possible
     Acha_Input%Chan_Idx_67um = 9     !channel number for 6.7
     Acha_Input%Chan_Idx_85um = 11     !channel number for 8.5
     Acha_Input%Chan_Idx_11um = 14     !channel number for 11
     Acha_Input%Chan_Idx_12um = 15     !channel number for 12
     Acha_Input%Chan_Idx_133um = 16  !channel number for 13.3

     CALL NFIA_CloudHeight_ChanOn(Ctxt%CLOUD_HEIGHT_Src1_T00, Chan_On_Flag_Default)
     Acha_Input%Chan_On_67um = Chan_On_Flag_Default(Acha_Input%Chan_Idx_67um)
     Acha_Input%Chan_On_85um = Chan_On_Flag_Default(Acha_Input%Chan_Idx_85um)
     Acha_Input%Chan_On_11um = Chan_On_Flag_Default(Acha_Input%Chan_Idx_11um)
     Acha_Input%Chan_On_12um = Chan_On_Flag_Default(Acha_Input%Chan_Idx_12um)
     Acha_Input%Chan_On_133um = Chan_On_Flag_Default(Acha_Input%Chan_Idx_133um)

     CALL NFIA_Sat_L1b_BadPixMsk(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14, Acha_Input%Invalid_Data_Mask)
     
     CALL NFIA_NWP_X_NWP(Ctxt%NWP_DATA_Src1_T00, Acha_Input%Elem_Idx_NWP)
     
     CALL NFIA_NWP_Y_NWP(Ctxt%NWP_DATA_Src1_T00, Acha_Input%Line_Idx_NWP)

!----- NEED TO ADD THESE ---------------     
     Acha_Input%Elem_Idx_Opposite_Corner_NWP => I_Nwp_x
     Acha_Input%Line_Idx_Opposite_Corner_NWP => J_Nwp_x
     Acha_Input%Longitude_Interp_Weight_NWP => Lon_Nwp_Fac
     Acha_Input%Latitude_Interp_Weight_NWP => Lat_Nwp_Fac
     Acha_Input%Viewing_Zenith_Angle_Idx_Rtm => Ivza_Rtm
!---------------------------------------

     if (Chan_On_Flag_Default(Chan_Idx_67um) == sym%YES) then
         CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI9,  Acha_Input%Bt_67um)
         CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI9,  Acha_Input%Rad_Clear_67um_Local)
     endif

      if (Chan_On_Flag_Default(Chan_Idx_85um) == sym%YES) then
         CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI11, Acha_Input%Bt_85um)
         CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI11, Acha_Input%Rad_Clear_85um_Local)
      endif

     if (Chan_On_Flag_Default(Chan_Idx_11um) == sym%YES) then
         CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14,  Acha_Input%Bt_11um)
         CALL NFIA_Sat_L1b_Rad(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14,  Acha_Input%Rad_11um)
         CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI14, Acha_Input%Rad_Clear_11um_Local)
     endif

      if (Chan_On_Flag_Default(Chan_Idx_12um) == sym%YES) then
         CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI15, Acha_Input%Bt_12um)
         CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI15, Acha_Input%Rad_Clear_12um_Local)
      endif
          
      if (Chan_On_Flag_Default(Chan_Idx_133um) == sym%YES) then
         CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI16, Acha_Input%Bt_133um)
         CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI16, Acha_Input%Rad_Clear_133um_Local)
      endif

     CALL NFIA_PseudoEmiss_Chn7(Ctxt%PSEUDO_EMISSIVITY_Src1_T00, Acha_Input%Surface_Emissivity_39um)

     CALL NFIA_Sat_Nav_COS_SatZen(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Acha_Input%Cosine_Zenith_Angle)
     
     CALL NFIA_Sat_Nav_SatZen(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Acha_Input%Sensor_Zenith_Angle)
     
!----- NEED TO ADD THIS TO FRAMEWORK ---------------     
     Acha_Input%Sensor_Azimuth_Angle => NULL()
!--------------------------------------     

     CALL NFIA_Sat_Nav_Lat(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Acha_Input%Latitude)
     
     CALL NFIA_Sat_Nav_Lon(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Acha_Input%Longitude)

     CALL NFIA_SnowMask_Mask(Ctxt%SNOW_MASK_Src1_T00, Acha_Input%Snow_Class)
     
     CALL NFIA_SfcType_Mask(Ctxt%SURFACE_TYPE_Src1_T00, Acha_Input%Surface_Type)

     !NEED 2D array call for framework - WCS3!
     CALL NFIP_NWP_TempSfc(Ctxt%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Acha_Input%Surface_Temperature)
     !------------------!
    
!----- NEED TO ADD THIS TO FRAMEWORK ---------------     
     Acha_Input%Surface_Air_Temperature => NULL()
!---------------------------------------

     !NEED 2D array call for framework - WCS3!
     CALL NFIP_NWP_TempTropo(Ctxt%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Acha_Input%Tropopause_Temperature)
     !------------------!

    CALL NFIA_NWP_PressSfc(Ctxt%NWP_DATA_Src1_T00, Acha_Input%Surface_Pressure)

     CALL NFIA_SfcElev_Elevation(Ctxt%SURFACE_ELEVATION_Src1_T00, Acha_Input%Surface_Elevation)
     
     CALL NFIA_CloudMask_Mask(Ctxt%CLOUD_MASK_Src1_T00, Acha_Input%Cloud_Mask_Local)
     
     CALL NFIA_CloudPhase_CldType(Ctxt%CLOUD_PHASE_Src1_T00, Acha_Input%Cloud_Type_Local)



!optional, if calculated. If null, then ACHA will calculate internal
!NOT IN FRAMEWORK. Not needed since ACHA will calculate
!     Acha_Input%Elem_Idx_LRC_Input => null()
!     Acha_Input%Line_Idx_LRC_Input =>  null()



!---- initalize Output structure

!----- NEED TO ADD THESE ---------------     
     Acha_Output%Latitude_Pc => Lat_Pc
     Acha_Output%Longitude_Pc => Lon_Pc
!---------------------------------------

     CALL NFIA_CloudHeight_CldTopTemp(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Tc)
     
     CALL NFIA_CloudHeight_CldTopEmss(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Ec)
     
     CALL NFIA_CloudHeight_CldBta1112(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Beta)
     
     CALL NFIA_CloudHeight_CldTopPres(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Pc)
     
     CALL NFIA_CloudHeight_CldTopHght(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Zc)
     
     CALL NFIA_CloudHeight_CldOptDpth(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Tau)
     
     CALL NFIA_CloudHeight_CldTopReff(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Reff)
     
     CALL NFIA_CloudHeight_TcError(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Tc_Uncertainty)
     
     CALL NFIA_CloudHeight_EcError(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Ec_Uncertainty)
     
     CALL NFIA_CloudHeight_B1112Error(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Beta_Uncertainty)
     
     CALL NFIA_CloudHeight_PcError(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Pc_Uncertainty)
     
     CALL NFIA_CloudHeight_ZcError(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Zc_Uncertainty)
     
     CALL NFIA_CloudHeight_CldLayer(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Cloud_Layer)
     
     CALL NFIA_CloudHeight_PcLowerCld(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Lower_Cloud_Pressure)
     
     CALL NFIA_CloudHeight_PcLowerCld(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Lower_Cloud_Pressure)
     
     CALL NFIA_CloudHeight_ZcLowerCld(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Lower_Cloud_Height)

!----- NEED TO ADD THESE ---------------
     Acha_Output%Zc_Top => Zc_Top_Acha
     Acha_Output%Zc_Base => Zc_Base_Acha
!---------------------------------------

     CALL NFIA_CloudHeight_CldHgtQF(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Qf)
     CALL NFIA_CloudHeight_Flag(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%OE_Qf)
     
!----- NEED TO ADD THESE ---------------     
     Acha_Output%Packed_Qf => Acha_Packed_Quality_Flags
     Acha_Output%Packed_Meta_Data => Acha_Packed_Meta_Data_Flags
!---------------------------------------

     CALL NFIA_CloudHeight_ProcOrder(Ctxt%CLOUD_HEIGHT_Src1_T00, Acha_Output%Processing_Order)


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
     Acha_Input%Cloud_Mask_Local =>  NULL()
     Acha_Input%Cloud_Type_Local =>  NULL()
     Acha_Input%Rad_Clear_67um_Local =>  NULL()
     Acha_Input%Rad_Clear_85um_Local =>  NULL()
     Acha_Input%Rad_Clear_11um_Local =>  NULL()
     Acha_Input%Rad_Clear_12um_Local =>  NULL()
     Acha_Input%Rad_Clear_133um_Local =>  NULL()
     Acha_Input%Surface_Emissivity_39um =>  NULL()
     Acha_Input%Elem_Idx_LRC_Input =>  NULL()
     Acha_Input%Line_Idx_LRC_Input =>   NULL()
#if defined (ISFWAIT)
     Acha_Input%Ctxt => NULL()   
#endif

                                   
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
