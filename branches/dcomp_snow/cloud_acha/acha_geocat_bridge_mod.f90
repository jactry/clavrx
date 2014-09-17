!$Id: ACHA_GEOCAT_BRIDGE_MOD.f90 9 2014-01-31 08:19:35Z awalther $
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Height Algorithm (ACHA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  processing systems and the ACHA code.
!
!------------------------------------------------------------------------------
module ACHA_GEOCAT_BRIDGE_MOD

 use ACHA_SERVICES_MOD
 use AWG_CLOUD_HEIGHT
   
 implicit none

 public :: AWG_CLOUD_HEIGHT_BRIDGE
 private :: NULL_ACHA_POINTERS
 private :: WMO_Sensor_KM

 contains


!----------------------------------------------------------------------
! BEGINNING OF ACHA SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine AWG_CLOUD_HEIGHT_BRIDGE(Ialgo)
 
   implicit none
   INTEGER(KIND=Int4), INTENT(IN) :: Ialgo

   !--------------------------------------------------------------------
   ! define structures that will be arguments to ACHA
   !--------------------------------------------------------------------
   type(symbol_acha) :: symbol
   type(acha_input_struct) :: Acha_Input
   type(acha_output_struct) :: Acha_Output

   INTEGER(KIND=Int4) :: Num_lines, Num_Elems
   INTEGER(KIND=Int4) :: wmo_id

   !---null pointers before filling them
   call NULL_ACHA_POINTERS(Acha_Input, Acha_Output)

   !-----------------------------------------------------------------------
   !---  GEOCAT Bridge Section
   !-----------------------------------------------------------------------
   !--- initialize Input structure pointers

   !--- store integer values
   Acha_Input%Number_of_Elements = sat%nx
   Acha_Input%Number_of_Lines = sat%nscans_per_segment
   Acha_Input%Num_Line_Max = sat%nscans_per_segment
   Acha_Input%Process_Undetected_Cloud_Flag = sym%NO
   Acha_Input%Smooth_Nwp_Fields_Flag = sym%NO
   Acha_Input%ACHA_Mode_Flag_In = -1
   !Sensor resolution
   Acha_Input%Sensor_Resolution_KM = WMO_Sensor_KM(scinfo(sc_ind)%WMO_Sc_Id)
   
   Acha_Input%Chan_Idx_67um = 9     !channel number for 6.7
   Acha_Input%Chan_Idx_85um = 11     !channel number for 8.5
   Acha_Input%Chan_Idx_11um = 14     !channel number for 11
   Acha_Input%Chan_Idx_12um = 15     !channel number for 12
   Acha_Input%Chan_Idx_133um = 16  !channel number for 13.3

   Acha_Input%Chan_On_67um = out2(Ialgo)%ch_flg(9)
   Acha_Input%Chan_On_85um = out2(Ialgo)%ch_flg(11)
   Acha_Input%Chan_On_11um = out2(Ialgo)%ch_flg(14)
   Acha_Input%Chan_On_12um = out2(Ialgo)%ch_flg(15)
   Acha_Input%Chan_On_133um = out2(Ialgo)%ch_flg(16)

   Acha_Input%Invalid_Data_Mask => sat%bad_pixel_mask(14,:,:)
   Acha_Input%Elem_Idx_Nwp =>  sat%x_nwp
   Acha_Input%Line_Idx_Nwp => sat%y_nwp      
 
   Acha_Input%Viewing_Zenith_Angle_Idx_Rtm => sat%ivza
   Acha_Input%Bt_67um => sat%bt9
   Acha_Input%Bt_85um => sat%bt11
   Acha_Input%Bt_11um => sat%bt14
   Acha_Input%Bt_12um => sat%bt15
   Acha_Input%Bt_133um => sat%bt16

   Acha_Input%Rad_11um => sat%rad14
   Acha_Input%Cosine_Zenith_Angle => sat%cos_satzen
   Acha_Input%Sensor_Zenith_Angle => sat%satzen
   Acha_Input%Sensor_Azimuth_Angle => sat%sataz
   Acha_Input%Latitude => sat%lat
   Acha_Input%Longitude => sat%lon

   Acha_Input%Snow_Class => sat%snow_mask
   Acha_Input%Surface_Type => sat%sfc_type

   Acha_Input%Surface_Elevation => sat%Zsfc
   Acha_Input%Cloud_Mask => sat%cldmask
   Acha_Input%Cloud_Type => sat%Cldtype

   Acha_Input%Rad_Clear_67um => sat%rad_clr9
   Acha_Input%Rad_Clear_85um => sat%rad_clr11
   Acha_Input%Rad_Clear_11um => sat%rad_clr14
   Acha_Input%Rad_Clear_12um => sat%rad_clr15
   Acha_Input%Rad_Clear_133um => sat%rad_clr16
   Acha_Input%Surface_Emissivity_39um => sat%sfc_emiss7
   
   !GEOCAT is special, where we have to allocate the pixel level data
   ! instead of pointing to it

   Num_lines = Acha_Input%Number_of_Lines
   Num_Elems = Acha_Input%Number_of_Elements

   allocate (Acha_Input%Surface_Temperature(Num_Elems,Num_lines))
   allocate (Acha_Input%Surface_Air_Temperature(Num_Elems,Num_lines))
   allocate (Acha_Input%Tropopause_Temperature(Num_Elems,Num_lines))
   allocate (Acha_Input%Surface_Pressure(Num_Elems,Num_lines))
   allocate (Acha_Input%Elem_Idx_Opposite_Corner_NWP(Num_Elems,Num_lines))
   allocate (Acha_Input%Line_Idx_Opposite_Corner_NWP(Num_Elems,Num_lines))
   allocate (Acha_Input%Longitude_Interp_Weight_NWP(Num_Elems,Num_lines))
   allocate (Acha_Input%Latitude_Interp_Weight_NWP(Num_Elems,Num_lines))


   Acha_Input%Elem_Idx_Opposite_Corner_NWP = MISSING_VALUE_INT4
   Acha_Input%Line_Idx_Opposite_Corner_NWP = MISSING_VALUE_INT4
   Acha_Input%Longitude_Interp_Weight_NWP = MISSING_VALUE_REAL4
   Acha_Input%Latitude_Interp_Weight_NWP = MISSING_VALUE_REAL4
   


      
   CALL ACHA_NWP_Fill(ACHA_Input)
   
   !LRC will be done inside the science code
   
   Acha_Input%Elem_Idx_LRC_Input => null()
   Acha_Input%Line_Idx_LRC_Input =>  null()
     
   !---- initalize Output structure
   Acha_Output%Latitude_Pc => out2(Ialgo)%Lat_Pc
   Acha_Output%Longitude_Pc => out2(Ialgo)%Lon_Pc

   Acha_Output%Tc => out2(Ialgo)%Cldt
   Acha_Output%Ec => out2(Ialgo)%cldemiss
   Acha_Output%Beta => out2(Ialgo)%cldbeta1112
   Acha_Output%Pc =>  out2(Ialgo)%Cldp
   Acha_Output%Zc => out2(Ialgo)%cldz
   Acha_Output%Tau => out2(Ialgo)%cod_vis
   Acha_Output%Reff => out2(Ialgo)%cldreff
   Acha_Output%Tc_Uncertainty => out2(Ialgo)%Tc_error
   Acha_Output%Ec_Uncertainty => out2(Ialgo)%ec_error
   Acha_Output%Beta_Uncertainty => out2(Ialgo)%beta1112_error
   Acha_Output%Pc_Uncertainty => out2(Ialgo)%pc_error
   Acha_Output%Zc_Uncertainty => out2(Ialgo)%zc_error
   Acha_Output%Cloud_Layer => out2(Ialgo)%Cld_Layer
   Acha_Output%Lower_Cloud_Pressure => out2(Ialgo)%Pc_Lower_Cloud
   Acha_Output%Lower_Cloud_Temperature => out2(Ialgo)%Tc_Lower_Cloud
   Acha_Output%Lower_Cloud_Height => out2(Ialgo)%Zc_Lower_Cloud
  
   Acha_Output%Zc_Top => out2(Ialgo)%Zc_Top_Acha
   Acha_Output%Zc_Base => out2(Ialgo)%Zc_Base_Acha
  
   Acha_Output%Qf => out2(Ialgo)%Cloud_Height_QF
   Acha_Output%OE_Qf => out2(Ialgo)%qcflg1
   Acha_Output%Packed_Qf => out2(Ialgo)%Acha_Packed_Quality_Flags
   Acha_Output%Packed_Meta_Data => out2(Ialgo)%Acha_Packed_Meta_Data_Flags
   Acha_Output%Processing_Order  => out2(Ialgo)%Processing_Order
   Acha_Output%Cost  => out2(Ialgo)%Acha_Cost
  

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
   symbol%PROB_CLEAR_TYPE = sym%CLEAR_TYPE
   symbol%FOG_TYPE = sym%FOG_TYPE
   symbol%WATER_TYPE = sym%WATER_TYPE
   symbol%SUPERCOOLED_TYPE = sym%SUPERCOOLED_TYPE
   symbol%MIXED_TYPE = sym%MIXED_TYPE
   symbol%OPAQUE_ICE_TYPE = sym%TICE_TYPE
   symbol%TICE_TYPE = sym%TICE_TYPE
   symbol%CIRRUS_TYPE = sym%CIRRUS_TYPE
   symbol%OVERLAP_TYPE = sym%OVERLAP_TYPE
   symbol%OVERSHOOTING_TYPE = sym%OVERSHOOTING_TYPE
   symbol%UNKNOWN_TYPE = sym%UNKNOWN_TYPE
   symbol%DUST_TYPE = 11
   symbol%SMOKE_TYPE = 12
   symbol%FIRE_TYPE = 13

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
     
      ! deallocate NWP pixel vars
      If(allocated(Acha_Input%Surface_Temperature) ) deallocate (Acha_Input%Surface_Temperature)
      If(allocated(Acha_Input%Surface_Air_Temperature) ) deallocate (Acha_Input%Surface_Air_Temperature)
      If(allocated(Acha_Input%Tropopause_Temperature) ) deallocate (Acha_Input%Tropopause_Temperature)
      If(allocated(Acha_Input%Surface_Pressure) ) deallocate (Acha_Input%Surface_Pressure)

      If(allocated(Acha_Input%Elem_Idx_Opposite_Corner_NWP) ) deallocate (Acha_Input%Elem_Idx_Opposite_Corner_NWP)
      If(allocated(Acha_Input%Line_Idx_Opposite_Corner_NWP) ) deallocate (Acha_Input%Line_Idx_Opposite_Corner_NWP)
      If(allocated(Acha_Input%Longitude_Interp_Weight_NWP) ) deallocate (Acha_Input%Longitude_Interp_Weight_NWP)
      If(allocated(Acha_Input%Latitude_Interp_Weight_NWP) ) deallocate (Acha_Input%Latitude_Interp_Weight_NWP)
      
                                  
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
     Acha_Output%Cost  => NULL()
 
 
 end subroutine
 
 !-----------------------------------------------------------------------------
 !
 !-----------------------------------------------------------------------------
 function WMO_Sensor_KM(wmo_id) result (Sensor_KM)
      integer , intent(in) :: wmo_id
      real :: Sensor_KM

      select case (WMO_id)
      case(3)   ! Metop-01
            Sensor_KM = 1.0
      case(4)   ! Metop-02
            Sensor_KM = 1.0
      case(55)  ! Meteosat-08
             Sensor_KM = 3.0
      case(56)  ! Meteosat-09
             Sensor_KM = 3.0
      case(57)  ! Meteosat-10
             Sensor_KM = 3.0
      case(70)  ! Meteosat-11
             Sensor_KM = 3.0
      case(171) ! MTSAT-1R
             Sensor_KM = 4.0
      case(172) ! MTSAT-2
             Sensor_KM = 4.0
      case(200:209) ! NOAA-08 - NOAA-18
             Sensor_KM = 1.0
      case(223) ! NOAA-19
             Sensor_KM = 1.0
      case(224) ! NPP
             Sensor_KM = 0.75
      case(252) ! GOES-08
             Sensor_KM = 4.0
      case(253) ! GOES-09
             Sensor_KM = 4.0
      case(254) ! GOES-10
             Sensor_KM = 4.0
      case(255) ! GOES-11
             Sensor_KM = 4.0
      case(256) ! GOES-12
             Sensor_KM = 4.0
      case(257) ! GOES-13
             Sensor_KM = 4.0
      case(258) ! GOES-14
             Sensor_KM = 4.0
      case(259) ! GOES-15
             Sensor_KM = 4.0
      case(152) ! GMS-5
             Sensor_KM = 2.0
      case(783) ! MODIS Terra
             Sensor_KM = 1.0
      case(784) ! MODIS Aqua
             Sensor_KM = 1.0
      case default
             print*,'This sensor is missing a wmo id: ', wmo_id 
             stop
      end select


 end function WMO_Sensor_KM

end module ACHA_GEOCAT_BRIDGE_MOD
