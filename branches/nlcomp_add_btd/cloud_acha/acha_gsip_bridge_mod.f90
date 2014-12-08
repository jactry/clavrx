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
! public :: SET_ACHA_VERSION
   
   
 contains


 !====================================================================
 !  record cvs version as a global variable for output to hdf
 !
 !  THIS IS BROKEN- IT ONLY PASSES BRIDGE INFO - FIXME
 !
 !====================================================================
! subroutine SET_ACHA_VERSION()
!   Acha_Version = "$Id$"
! end subroutine SET_ACHA_VERSION


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
   !---  GSIP Bridge Section
   !-----------------------------------------------------------------------
   !--- initialize Input structure pointers

   !--- store integer values
   Acha_Input%Number_of_Elements = sat%nx
   Acha_Input%Number_of_Lines = Num_Scans_Per_Segment
   Acha_Input%Num_Line_Max = Num_Scans_Per_Segment
   Acha_Input%Process_Undetected_Cloud_Flag = sym%NO
   Acha_Input%Smooth_Nwp_Fields_Flag = Smooth_Nwp_Flag
   Acha_Input%ACHA_Mode_Flag_In = sat_info_gsip(1)%acha_mode
   !Sensor resolution
   Acha_Input%Sensor_Resolution_KM = WMO_Sensor_KM(sat_info_gsip(1)%WMO_Sc_Id)

   Acha_Input%Chan_Idx_67um = 9     !channel number for 6.7
   Acha_Input%Chan_Idx_85um = 11     !channel number for 8.5
   Acha_Input%Chan_Idx_11um = 14     !channel number for 11
   Acha_Input%Chan_Idx_12um = 15     !channel number for 12
   Acha_Input%Chan_Idx_133um = 16  !channel number for 13.3

   Acha_Input%Chan_On_67um = sat_info_gsip(1)%chanon(9)
   Acha_Input%Chan_On_85um = sat_info_gsip(1)%chanon(11)
   Acha_Input%Chan_On_11um = sat_info_gsip(1)%chanon(14)
   Acha_Input%Chan_On_12um = sat_info_gsip(1)%chanon(15)
   Acha_Input%Chan_On_133um = sat_info_gsip(1)%chanon(16)

   Acha_Input%Invalid_Data_Mask => bad_pix_mask(14,:,:)

   Acha_Input%Elem_Idx_Nwp =>  I_Nwp
   Acha_Input%Line_Idx_Nwp => J_Nwp
   Acha_Input%Elem_Idx_Opposite_Corner_NWP => I_Nwp_x
   Acha_Input%Line_Idx_Opposite_Corner_NWP => J_Nwp_x
   Acha_Input%Longitude_Interp_Weight_NWP => Lon_Nwp_Fac
   Acha_Input%Latitude_Interp_Weight_NWP => Lat_Nwp_Fac
   Acha_Input%Viewing_Zenith_Angle_Idx_Rtm => Ivza_Rtm

   if (Acha_Input%Chan_On_67um  == sym%YES) Acha_Input%Bt_67um => bt9
   if (Acha_Input%Chan_On_85um  == sym%YES) Acha_Input%Bt_85um => bt11
   if (Acha_Input%Chan_On_11um  == sym%YES) Acha_Input%Bt_11um => bt14
   if (Acha_Input%Chan_On_12um  == sym%YES) Acha_Input%Bt_12um => bt15
   if (Acha_Input%Chan_On_133um  == sym%YES) Acha_Input%Bt_133um => bt16

   if (Acha_Input%Chan_On_11um  == sym%YES) Acha_Input%Rad_11um => rad14
   Acha_Input%Cosine_Zenith_Angle => Coszen
   Acha_Input%Sensor_Zenith_Angle => Satzen
   Acha_Input%Sensor_Azimuth_Angle => Sataz
   Acha_Input%Latitude => Lat
   Acha_Input%Longitude => Lon

   Acha_Input%Snow_Class => snow_mask
   Acha_Input%Surface_Type => Sfc_Type

   Acha_Input%Surface_Temperature =>Tsfc_Nwp_Pix
   Acha_Input%Surface_Air_Temperature => Tair_Nwp_Pix
   Acha_Input%Tropopause_Temperature => Ttropo_Nwp_Pix
   Acha_Input%Surface_Pressure => Psfc_Nwp_Pix

   Acha_Input%Surface_Elevation => Zsfc
   Acha_Input%Cloud_Mask => gsip_pix_prod%cldmask
   Acha_Input%Cloud_Type => gsip_pix_prod%Cldtype

   Acha_Input%Rad_Clear_67um => Rad_Clear_Ch9_Rtm
   Acha_Input%Rad_Clear_85um => Rad_Clear_Ch11_Rtm
   Acha_Input%Rad_Clear_11um => Rad_Clear_Ch14_Rtm
   Acha_Input%Rad_Clear_12um => Rad_Clear_Ch15_Rtm
   Acha_Input%Rad_Clear_133um => Rad_Clear_Ch16_Rtm
   Acha_Input%Surface_Emissivity_39um => sfc_emiss_7

   Acha_Input%Elem_Idx_LRC_Input => null()
   Acha_Input%Line_Idx_LRC_Input =>  null()
     
   !---- initalize Output structure
   Acha_Output%Latitude_Pc => gsip_pix_prod%Lat_Pc
   Acha_Output%Longitude_Pc => gsip_pix_prod%Lon_Pc
   Acha_Output%Tc => gsip_pix_prod%ctt
   Acha_Output%Ec => gsip_pix_prod%cldemiss
   Acha_Output%Beta => gsip_pix_prod%r4_generic1
   Acha_Output%Pc => gsip_pix_prod%Cldp
   Acha_Output%Zc => gsip_pix_prod%cldz
   Acha_Output%Tau => gsip_pix_prod%cod_acha
   Acha_Output%Reff => gsip_pix_prod%r4_generic1
   Acha_Output%Tc_Uncertainty =>  gsip_pix_prod%Tc_error
   Acha_Output%Ec_Uncertainty => gsip_pix_prod%ec_error
   Acha_Output%Beta_Uncertainty => gsip_pix_prod%beta1112_error
   Acha_Output%Pc_Uncertainty => gsip_pix_prod%pc_error
   Acha_Output%Zc_Uncertainty => gsip_pix_prod%zc_error
   Acha_Output%Lower_Cloud_Pressure =>  gsip_pix_prod%Pc_Lower_Cloud
   Acha_Output%Lower_Cloud_Temperature => gsip_pix_prod%Tc_Lower_Cloud
   Acha_Output%Lower_Cloud_Height => gsip_pix_prod%Zc_Lower_Cloud
   Acha_Output%Zc_Top => gsip_pix_prod%Zc_Top_Acha
   Acha_Output%Zc_Base => gsip_pix_prod%Zc_Base_Acha
   Acha_Output%Qf => gsip_pix_prod%Cloud_Height_QF
   Acha_Output%OE_Qf => gsip_pix_prod%qcflg1
   Acha_Output%Packed_Qf => gsip_pix_prod%Acha_Packed_Quality_Flags
   Acha_Output%Packed_Meta_Data => gsip_pix_prod%Acha_Packed_Meta_Data_Flags
   Acha_Output%Processing_Order  => gsip_pix_prod%Processing_Order
   Acha_Output%Cost  => gsip_pix_prod%r4_generic2
  

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
     Acha_Output%Cost  => NULL()
 
 
 end subroutine

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
