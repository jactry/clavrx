!$Id$
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Height Algorithm (ACHA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  processing systems and the ACHA code.
!
!------------------------------------------------------------------------------
module ACHA_GSIP_BRIDGE_MOD

 use AWG_CLOUD_HEIGHT
 use ACHA_SERVICES_MOD
   
 implicit none

 public :: AWG_CLOUD_HEIGHT_BRIDGE
 private:: SET_SYMBOL, SET_INPUT, SET_OUTPUT, NULL_INPUT, NULL_OUTPUT

 !--------------------------------------------------------------------
 ! define structures that will be arguments to ACHA
 !--------------------------------------------------------------------
 type(acha_symbol_struct), private :: Symbol
 type(acha_input_struct), private :: Input
 type(acha_output_struct), private :: Output

 contains

!----------------------------------------------------------------------
! ACHA BRIDGE SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine AWG_CLOUD_HEIGHT_BRIDGE()
 
   implicit none

   !---null pointers before filling them
   call NULL_INPUT()
   call NULL_OUTPUT()

   !-------------------------------------------
   !--- initialize structures
   !-------------------------------------------

   !--- store integer values
   call SET_INPUT()

   !---- initalize Output structure
   call SET_OUTPUT()
   print *, "Height out done"
  
   !----set symbols to local values
   call SET_SYMBOL()

   !-----------------------------------------------------------------------
   !--- Call to AWG CLoud Height Algorithm (ACHA)
   !-----------------------------------------------------------------------
   call AWG_CLOUD_HEIGHT_ALGORITHM(Input, &
                                    Symbol, &
                                    Output)

   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT()
   call NULL_OUTPUT()

 end subroutine AWG_CLOUD_HEIGHT_BRIDGE

 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding output to ACHA
 !-----------------------------------------------------------------------------
 subroutine NULL_INPUT()
     Input%Invalid_Data_Mask =>  null()
     Input%Bt_67um =>  null()
     Input%Bt_85um =>  null()
     Input%Bt_11um =>  null()
     Input%Bt_12um =>  null()
     Input%Bt_133um =>  null()
     Input%Rad_67um =>  null()
     Input%Rad_11um =>  null()
     Input%Cosine_Zenith_Angle =>  null()
     Input%Sensor_Zenith_Angle =>  null()
     Input%Sensor_Azimuth_Angle =>  null()
     Input%Surface_Temperature => null()
     Input%Surface_Air_Temperature =>  null()
     Input%Tropopause_Temperature =>  null()
     Input%Surface_Pressure =>  null()
     Input%Surface_Elevation =>  null()
     Input%Latitude =>  null()
     Input%Longitude =>  null()
     Input%Rad_Clear_67um =>  null()
     Input%Rad_Clear_85um =>  null()
     Input%Rad_Clear_11um =>  null()
     Input%Rad_Clear_12um =>  null()
     Input%Rad_Clear_133um =>  null()
     Input%Surface_Emissivity_39um =>  null()
     Input%Snow_Class =>  null()
     Input%Surface_Type =>  null()
     Input%Cloud_Mask =>  null()
     Input%Cloud_Probability => null()
     Input%Cloud_Type =>  null()
     Input%Elem_Idx_Nwp =>   null()
     Input%Line_Idx_Nwp =>  null()
     Input%Elem_Idx_Opposite_Corner_NWP =>  null()
     Input%Line_Idx_Opposite_Corner_NWP =>  null()
     Input%Viewing_Zenith_Angle_Idx_Rtm =>  null()
     Input%Latitude_Interp_Weight_NWP =>  null()
     Input%Longitude_Interp_Weight_NWP =>  null()
     Input%Elem_Idx_LRC_Input =>  null()
     Input%Line_Idx_LRC_Input =>   null()
     Input%Tc_Cirrus_Sounder =>   null()
 end subroutine NULL_INPUT
 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding output to ACHA
 !-----------------------------------------------------------------------------
 subroutine NULL_OUTPUT()
     Output%Latitude_Pc =>  null()
     Output%Longitude_Pc =>  null()
     Output%Tc =>  null()
     Output%Ec =>  null()
     Output%Beta =>  null()
     Output%Pc =>  null()
     Output%Zc =>  null()
     Output%Tau =>  null()
     Output%Reff =>  null()
     Output%Tc_Uncertainty =>  null()
     Output%Ec_Uncertainty =>  null()
     Output%Beta_Uncertainty =>  null()
     Output%Pc_Uncertainty =>  null()
     Output%Zc_Uncertainty =>  null()
     Output%Lower_Cloud_Pressure =>  null()
     Output%Lower_Cloud_Temperature =>  null()
     Output%Lower_Cloud_Height =>  null()
     Output%Cost =>  null()
     Output%Total_Cloud_Fraction =>  null()
     Output%Total_Cloud_Fraction_Uncer =>  null()
     Output%High_Cloud_Fraction =>  null()
     Output%Mid_Cloud_Fraction =>  null()
     Output%Low_Cloud_Fraction =>  null()
     Output%Cloud_Layer =>  null()
     Output%Qf =>  null()
     Output%OE_Qf =>  null()
     Output%Packed_Qf =>  null()
     Output%Packed_Meta_Data =>  null()
     Output%Processing_Order  =>  null()
     Output%Pc_Opaque =>  null()
     Output%Tc_Opaque =>  null()
     Output%Zc_Opaque =>  null()
     Output%Pc_H2O =>  null()
     Output%Tc_H2O =>  null()
     Output%Zc_H2O =>  null()
 end subroutine NULL_OUTPUT

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
   Output%Latitude_Pc => gsip_pix_prod%Lat_Pc
   Output%Longitude_Pc => gsip_pix_prod%Lon_Pc
   Output%Tc => gsip_pix_prod%ctt
   Output%Ec => gsip_pix_prod%cldemiss
   Output%Beta => gsip_pix_prod%r4_generic1
   Output%Pc => gsip_pix_prod%Cldp
   Output%Zc => gsip_pix_prod%cldz
   Output%Tau => gsip_pix_prod%cod_acha
   Output%Reff => gsip_pix_prod%r4_generic1
   Output%Tc_Uncertainty => gsip_pix_prod%Tc_error
   Output%Ec_Uncertainty => gsip_pix_prod%ec_error
   Output%Beta_Uncertainty => gsip_pix_prod%beta1112_error
   Output%Pc_Uncertainty =>  gsip_pix_prod%pc_error
   Output%Zc_Uncertainty => gsip_pix_prod%zc_error
   Output%Lower_Cloud_Pressure => gsip_pix_prod%Pc_Lower_Cloud
   Output%Lower_Cloud_Temperature => gsip_pix_prod%Tc_Lower_Cloud
   Output%Lower_Cloud_Height => gsip_pix_prod%Zc_Lower_Cloud
   Output%Qf => gsip_pix_prod%Cloud_Height_QF
   Output%OE_Qf => gsip_pix_prod%qcflg1
   Output%Packed_Qf => gsip_pix_prod%Acha_Packed_Quality_Flags
   Output%Packed_Meta_Data => gsip_pix_prod%Acha_Packed_Meta_Data_Flags
   Output%Processing_Order  => gsip_pix_prod%Processing_Order
   Output%Cost  => gsip_pix_prod%r4_generic2
   Output%Cloud_Layer  => gsip_pix_prod%cldlayer
   Output%Total_Cloud_Fraction => gsip_pix_prod%r4_generic2
   Output%Total_Cloud_Fraction_Uncer => gsip_pix_prod%r4_generic2
   Output%High_Cloud_Fraction => gsip_pix_prod%r4_generic2
   Output%Mid_Cloud_Fraction => gsip_pix_prod%r4_generic2
   Output%Low_Cloud_Fraction => gsip_pix_prod%r4_generic2
   Output%Pc_Opaque => gsip_pix_prod%r4_generic2
   Output%Tc_Opaque => gsip_pix_prod%r4_generic2
   Output%Zc_Opaque => gsip_pix_prod%r4_generic2
   Output%Pc_H2O => gsip_pix_prod%r4_generic3
   Output%Tc_H2O => gsip_pix_prod%r4_generic3
   Output%Zc_H2O => gsip_pix_prod%r4_generic3
 end subroutine SET_OUTPUT

 subroutine SET_INPUT()
   Input%Number_of_Elements = sat%nx
   Input%Number_of_Lines = Num_Scans_Per_Segment
   Input%Num_Line_Max = Num_Scans_Per_Segment
   Input%Process_Undetected_Cloud_Flag = sym%NO
   Input%Smooth_Nwp_Fields_Flag = Smooth_Nwp_Flag
   Input%ACHA_Mode_Flag_In = sat_info_gsip(1)%acha_mode
   !Sensor resolution
   Input%Sensor_Resolution_KM = WMO_Sensor_KM(sat_info_gsip(1)%WMO_Sc_Id)

   Input%Chan_Idx_67um = 9     !channel number for 6.7
   Input%Chan_Idx_85um = 11     !channel number for 8.5
   Input%Chan_Idx_11um = 14     !channel number for 11
   Input%Chan_Idx_12um = 15     !channel number for 12
   Input%Chan_Idx_133um = 16  !channel number for 13.3

   Input%Chan_On_67um = sat_info_gsip(1)%chanon(9)
   Input%Chan_On_85um = sat_info_gsip(1)%chanon(11)
   Input%Chan_On_11um = sat_info_gsip(1)%chanon(14)
   Input%Chan_On_12um = sat_info_gsip(1)%chanon(15)
   Input%Chan_On_133um = sat_info_gsip(1)%chanon(16)

   Input%Invalid_Data_Mask => bad_pix_mask(14,:,:)

   Input%Elem_Idx_Nwp =>  I_Nwp
   Input%Line_Idx_Nwp => J_Nwp
   Input%Elem_Idx_Opposite_Corner_NWP => I_Nwp_x
   Input%Line_Idx_Opposite_Corner_NWP => J_Nwp_x
   Input%Longitude_Interp_Weight_NWP => Lon_Nwp_Fac
   Input%Latitude_Interp_Weight_NWP => Lat_Nwp_Fac
   Input%Viewing_Zenith_Angle_Idx_Rtm => Ivza_Rtm
   
   if (Input%Chan_On_67um  == sym%YES) then
        Input%Bt_67um => bt9
        Input%Rad_67um => rad9
        Input%Rad_Clear_67um => Rad_Clear_Ch9_Rtm
  endif

   if (Input%Chan_On_85um  == sym%YES) then
        Input%Bt_85um => bt11
        Input%Rad_Clear_85um => Rad_Clear_Ch11_Rtm
   endif
   
   if (Input%Chan_On_11um  == sym%YES) then
        Input%Bt_11um => bt14
        Input%Rad_11um => rad14
        Input%Rad_Clear_11um => Rad_Clear_Ch14_Rtm
   endif
   
   if (Input%Chan_On_12um  == sym%YES) then
        Input%Bt_12um => bt15
        Input%Rad_Clear_12um => Rad_Clear_Ch15_Rtm
   endif
   
   if (Input%Chan_On_133um  == sym%YES) then
        Input%Bt_133um => bt16
        Input%Rad_Clear_133um => Rad_Clear_Ch16_Rtm        
   endif

   
   
   Input%Cosine_Zenith_Angle => Coszen
   Input%Sensor_Zenith_Angle => Satzen
   Input%Sensor_Azimuth_Angle => Sataz
   Input%Latitude => Lat
   Input%Longitude => Lon

   Input%Snow_Class => snow_mask
   Input%Surface_Type => Sfc_Type

   Input%Surface_Temperature =>Tsfc_Nwp_Pix
   Input%Surface_Air_Temperature => Tair_Nwp_Pix
   Input%Tropopause_Temperature => Ttropo_Nwp_Pix
   Input%Surface_Pressure => Psfc_Nwp_Pix

   Input%Surface_Elevation => Zsfc
   Input%Cloud_Mask =>  gsip_pix_prod%cldmask
   Input%Cloud_Type => gsip_pix_prod%Cldtype
   Input%Cloud_Probability => gsip_pix_prod%cldprob
   
   Input%Surface_Emissivity_39um => sfc_emiss_7

   Input%Elem_Idx_LRC_Input => I_LRC
   Input%Line_Idx_LRC_Input =>  J_LRC
   Input%Tc_Cirrus_Sounder =>  null()
 end subroutine SET_INPUT


 
 !-----------------------------------------------------------------------------
 !
 !-----------------------------------------------------------------------------
 function WMO_Sensor_KM(wmo_id) result (Sensor_KM)
      integer , intent(in) :: wmo_id
      real :: Sensor_KM

      select case (WMO_id)
      case(3)
	     Sensor_KM = 1.0
	   case(4)
	     Sensor_KM = 1.0
      case(55)
	     Sensor_KM = 3.0
      case(56)
	     Sensor_KM = 3.0
      case(57)
	     Sensor_KM = 3.0
      case(70)
	     Sensor_KM = 3.0
      case(171)
	     Sensor_KM = 4.0
	  case (172)
	     Sensor_KM = 13.0
	  case (200:209)
	     Sensor_KM = 1.0
	  case(223)
	     Sensor_KM = 1.0
      case(224)
	     Sensor_KM = 0.75
      case(252)
	     Sensor_KM = 4.0
      case(253)
	     Sensor_KM = 4.0
      case(254)
	     Sensor_KM = 4.0
      case(255)
	     Sensor_KM = 4.0
      case(256)
	     Sensor_KM = 4.0
      case(257)
	     Sensor_KM = 4.0
      case(258)
	     Sensor_KM = 4.0
      case(259)
	     Sensor_KM = 4.0
      case(152)
	     Sensor_KM = 2.0
      case(783)
	     Sensor_KM = 1.0
      case(784)
	     Sensor_KM = 1.0
	  case default
	     print*,'Please inform William that this sensor is missing in ACHA! ( william.straka@ssec.wisc.edu) wmo id: ', wmo_id 
		 stop	 
      end select


 end function WMO_Sensor_KM
 

end module ACHA_GSIP_BRIDGE_MOD
