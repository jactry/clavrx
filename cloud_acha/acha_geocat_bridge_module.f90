!$Id: ACHA_GEOCAT_BRIDGE_MOD.f90 9 2014-01-31 08:19:35Z awalther $
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Height Algorithm (ACHA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  processing systems and the ACHA code.
!
!------------------------------------------------------------------------------
module ACHA_GEOCAT_BRIDGE_MODULE

 use ACHA_SERVICES_MOD
 use AWG_CLOUD_HEIGHT
   
 implicit none

 public :: AWG_CLOUD_HEIGHT_BRIDGE
 private :: NULL_INPUT
 private :: NULL_OUTPUT
 private :: SET_INPUT
 private :: SET_OUTPUT
 private :: SET_SYMBOL
 private :: WMO_Sensor_KM

 !--------------------------------------------------------------------
 ! define structures that will be arguments to ACHA
 !--------------------------------------------------------------------
 type(symbol_acha) :: Symbol
 type(acha_input_struct) :: Input
 type(acha_output_struct) :: Output

 contains

!----------------------------------------------------------------------
! ACHA BRIDGE SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine AWG_CLOUD_HEIGHT_BRIDGE(Ialgo)
 
   implicit none
   integer(kind=int4), intent(in) :: Ialgo
   integer(kind=int4) :: Num_lines
   integer(kind=int4) :: Num_Elems
   integer(kind=int4) :: wmo_id

   !---null pointers before filling them
   call NULL_INPUT(Ialgo)
   call NULL_OUTPUT(Ialgo)

   !--- initialize input structure pointers
   call SET_INPUT(Ialgo)

   !--- initialize output structure pointers
   call SET_OUTPUT(Ialgo)

   !----set symbols to local values
   call SET_SYMBOL(Ialgo)

   ! GEOCAT is special, where we have to allocate the pixel level data
   ! instead of pointing to it

   Num_lines = Input%Number_of_Lines
   Num_Elems = Input%Number_of_Elements

   allocate (Input%Surface_Temperature(Num_Elems,Num_lines))
   allocate (Input%Surface_Air_Temperature(Num_Elems,Num_lines))
   allocate (Input%Tropopause_Temperature(Num_Elems,Num_lines))
   allocate (Input%Surface_Pressure(Num_Elems,Num_lines))
   allocate (Input%Elem_Idx_Opposite_Corner_NWP(Num_Elems,Num_lines))
   allocate (Input%Line_Idx_Opposite_Corner_NWP(Num_Elems,Num_lines))
   allocate (Input%Longitude_Interp_Weight_NWP(Num_Elems,Num_lines))
   allocate (Input%Latitude_Interp_Weight_NWP(Num_Elems,Num_lines))

   Input%Elem_Idx_Opposite_Corner_NWP = MISSING_VALUE_INT4
   Input%Line_Idx_Opposite_Corner_NWP = MISSING_VALUE_INT4
   Input%Longitude_Interp_Weight_NWP = MISSING_VALUE_REAL4
   Input%Latitude_Interp_Weight_NWP = MISSING_VALUE_REAL4

   CALL ACHA_NWP_Fill(Input)
   
   !LRC will be done inside the science code
   
   Input%Elem_Idx_LRC_Input => null()
   Input%Line_Idx_LRC_Input => null()
   Input%Tc_Cirrus_Sounder  => null()

   !-----------------------------------------------------------------------
   !--- Call to AWG CLoud Height Algorithm (ACHA)
   !-----------------------------------------------------------------------
   call AWG_CLOUD_HEIGHT_ALGORITHM(Input, &
                                   Symbol, &
                                   Output)
    
   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT(Ialgo)
   call NULL_OUTPUT(Ialgo)

 end subroutine AWG_CLOUD_HEIGHT_BRIDGE

 !-----------------------------------------------------------------------------
 !
 !-----------------------------------------------------------------------------
 subroutine NULL_INPUT(Ialgo)

     implicit none
     integer(kind=int4), intent(in) :: Ialgo

     Input%Invalid_Data_Mask =>  null()
     Input%Elem_Idx_Nwp =>   null()
     Input%Line_Idx_Nwp =>  null()
     Input%Viewing_Zenith_Angle_Idx_Rtm =>  null()

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
     Input%Latitude =>  null()
     Input%Longitude =>  null()
     Input%Snow_Class =>  null()
     Input%Surface_Type =>  null()
     Input%Surface_Elevation =>  null()
     Input%Cloud_Mask =>  null()
     Input%Cloud_Type =>  null()
     Input%Rad_Clear_67um =>  null()
     Input%Rad_Clear_85um =>  null()
     Input%Rad_Clear_11um =>  null()
     Input%Rad_Clear_12um =>  null()
     Input%Rad_Clear_133um =>  null()
     Input%Surface_Emissivity_39um =>  null()
     Input%Elem_Idx_LRC_Input =>  null()
     Input%Line_Idx_LRC_Input =>   null()
     Input%Tc_Cirrus_Sounder =>  null()

     ! deallocate NWP pixel vars
     If(allocated(Input%Surface_Temperature) ) deallocate (Input%Surface_Temperature)
     If(allocated(Input%Surface_Air_Temperature) ) deallocate (Input%Surface_Air_Temperature)
     If(allocated(Input%Tropopause_Temperature) ) deallocate (Input%Tropopause_Temperature)
     If(allocated(Input%Surface_Pressure) ) deallocate (Input%Surface_Pressure)

     If(allocated(Input%Elem_Idx_Opposite_Corner_NWP) ) deallocate (Input%Elem_Idx_Opposite_Corner_NWP)
     If(allocated(Input%Line_Idx_Opposite_Corner_NWP) ) deallocate (Input%Line_Idx_Opposite_Corner_NWP)
     If(allocated(Input%Longitude_Interp_Weight_NWP) ) deallocate (Input%Longitude_Interp_Weight_NWP)
     If(allocated(Input%Latitude_Interp_Weight_NWP) ) deallocate (Input%Latitude_Interp_Weight_NWP)

 end subroutine NULL_INPUT

 !-----------------------------------------------------------------------------
 !
 !-----------------------------------------------------------------------------
 subroutine NULL_OUTPUT(Ialgo)

     implicit none
     integer(kind=int4), intent(in) :: Ialgo

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
     Output%Zc_Top =>  null()
     Output%Zc_Base =>  null()
     Output%Qf =>  null()
     Output%OE_Qf =>  null()
     Output%Packed_Qf =>  null()
     Output%Packed_Meta_Data =>  null()
     Output%Processing_Order  =>  null()
     Output%Cost  => null()
     Output%Pc_Opaque =>  null()
     Output%Tc_Opaque =>  null()
     Output%Zc_Opaque =>  null()
     Output%Pc_H2O =>  null()
     Output%Tc_H2O =>  null()
     Output%Zc_H2O =>  null()
!---stw     Output%Cloud_Layer =>  null()
!---stw     Output%Total_Cloud_Fraction =>  null()
!---stw     Output%Total_Cloud_Fraction_Uncer =>  null()
!---stw     Output%High_Cloud_Fraction =>  null()
!---stw     Output%Mid_Cloud_Fraction =>  null()
!---stw     Output%Low_Cloud_Fraction =>  null()

 end subroutine NULL_OUTPUT

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

 !-----------------------------------------------------------------------------
 !
 !-----------------------------------------------------------------------------
 subroutine SET_INPUT(Ialgo)

     implicit none
     integer(kind=int4), intent(in) :: Ialgo

     Input%Number_of_Elements = sat%nx
     Input%Number_of_Lines = sat%nscans_per_segment
     Input%Num_Line_Max = sat%nscans_per_segment
     Input%Process_Undetected_Cloud_Flag = sym%NO
     Input%Smooth_Nwp_Fields_Flag = sym%NO
     Input%ACHA_Mode_Flag_In = -1
     !Sensor resolution
     Input%Sensor_Resolution_KM = WMO_Sensor_KM(scinfo(sc_ind)%WMO_Sc_Id)

     Input%Chan_Idx_67um = 9     !channel number for 6.7
     Input%Chan_Idx_85um = 11    !channel number for 8.5
     Input%Chan_Idx_11um = 14    !channel number for 11
     Input%Chan_Idx_12um = 15    !channel number for 12
     Input%Chan_Idx_133um = 16   !channel number for 13.3

     Input%Chan_On_67um = out2(Ialgo)%ch_flg(9)
     Input%Chan_On_85um = out2(Ialgo)%ch_flg(11)
     Input%Chan_On_11um = out2(Ialgo)%ch_flg(14)
     Input%Chan_On_12um = out2(Ialgo)%ch_flg(15)
     Input%Chan_On_133um = out2(Ialgo)%ch_flg(16)

     Input%Invalid_Data_Mask => sat%bad_pixel_mask(14,:,:)
     Input%Elem_Idx_Nwp =>  sat%x_nwp
     Input%Line_Idx_Nwp => sat%y_nwp

     Input%Viewing_Zenith_Angle_Idx_Rtm => sat%ivza
     Input%Bt_67um => sat%bt9
     Input%Bt_85um => sat%bt11
     Input%Bt_11um => sat%bt14
     Input%Bt_12um => sat%bt15
     Input%Bt_133um => sat%bt16

     Input%Rad_67um => sat%rad9
     Input%Rad_11um => sat%rad14
     Input%Cosine_Zenith_Angle => sat%cos_satzen
     Input%Sensor_Zenith_Angle => sat%satzen
     Input%Sensor_Azimuth_Angle => sat%sataz
     Input%Latitude => sat%lat
     Input%Longitude => sat%lon

     Input%Snow_Class => sat%snow_mask
     Input%Surface_Type => sat%sfc_type

     Input%Surface_Elevation => sat%Zsfc
     Input%Cloud_Mask => sat%cldmask
     Input%Cloud_Type => sat%Cldtype

     Input%Rad_Clear_67um => sat%rad_clr9
     Input%Rad_Clear_85um => sat%rad_clr11
     Input%Rad_Clear_11um => sat%rad_clr14
     Input%Rad_Clear_12um => sat%rad_clr15
     Input%Rad_Clear_133um => sat%rad_clr16
     Input%Surface_Emissivity_39um => sat%sfc_emiss7

 end subroutine SET_INPUT

 !-----------------------------------------------------------------------------
 !
 !-----------------------------------------------------------------------------
 subroutine SET_OUTPUT(Ialgo)

     implicit none
     integer(kind=int4), intent(in) :: Ialgo

     !--- Reinitialize output structure to maintain continuity with the
     !--- CLAVRx framework. Some Output elements are initialized inside
     !--- of awg_cloud_height.f90

     !---STW Debug
     out2(Ialgo)%cod_vis = MISSING_VALUE_REAL4
     out2(Ialgo)%cldreff = MISSING_VALUE_REAL4
     out2(Ialgo)%Tc_error = MISSING_VALUE_REAL4
     out2(Ialgo)%ec_error = MISSING_VALUE_REAL4
     out2(Ialgo)%beta1112_error = MISSING_VALUE_REAL4
     out2(Ialgo)%pc_error = MISSING_VALUE_REAL4
     out2(Ialgo)%zc_error = MISSING_VALUE_REAL4
     out2(Ialgo)%Pc_Lower_Cloud = MISSING_VALUE_REAL4
     out2(Ialgo)%Tc_Lower_Cloud = MISSING_VALUE_REAL4
     out2(Ialgo)%Zc_Lower_Cloud = MISSING_VALUE_REAL4
     out2(Ialgo)%Zc_Top_Acha = MISSING_VALUE_REAL4
     out2(Ialgo)%Zc_Base_Acha = MISSING_VALUE_REAL4
     out2(Ialgo)%Processing_Order = MISSING_VALUE_INT4
     out2(Ialgo)%Acha_Cost = MISSING_VALUE_REAL4
     out2(Ialgo)%Pc_Opaque_Cloud = MISSING_VALUE_REAL4
     out2(Ialgo)%Tc_Opaque_Cloud = MISSING_VALUE_REAL4
     out2(Ialgo)%Zc_Opaque_Cloud = MISSING_VALUE_REAL4
     out2(Ialgo)%Pc_H2O_Cloud = MISSING_VALUE_REAL4
     out2(Ialgo)%Tc_H2O_Cloud = MISSING_VALUE_REAL4
     out2(Ialgo)%Zc_H2O_Cloud = MISSING_VALUE_REAL4
     !---STW End Debug

     Output%Latitude_Pc => out2(Ialgo)%Lat_Pc
     Output%Longitude_Pc => out2(Ialgo)%Lon_Pc

     Output%Tc => out2(Ialgo)%Cldt
     Output%Ec => out2(Ialgo)%cldemiss
     Output%Beta => out2(Ialgo)%cldbeta1112
     Output%Pc =>  out2(Ialgo)%Cldp
     Output%Zc => out2(Ialgo)%cldz
     Output%Tau => out2(Ialgo)%cod_vis
     Output%Reff => out2(Ialgo)%cldreff
     Output%Tc_Uncertainty => out2(Ialgo)%Tc_error
     Output%Ec_Uncertainty => out2(Ialgo)%ec_error
     Output%Beta_Uncertainty => out2(Ialgo)%beta1112_error
     Output%Pc_Uncertainty => out2(Ialgo)%pc_error
     Output%Zc_Uncertainty => out2(Ialgo)%zc_error
     Output%Lower_Cloud_Pressure => out2(Ialgo)%Pc_Lower_Cloud
     Output%Lower_Cloud_Temperature => out2(Ialgo)%Tc_Lower_Cloud
     Output%Lower_Cloud_Height => out2(Ialgo)%Zc_Lower_Cloud

     Output%Zc_Top => out2(Ialgo)%Zc_Top_Acha
     Output%Zc_Base => out2(Ialgo)%Zc_Base_Acha

     Output%Qf => out2(Ialgo)%Cloud_Height_QF
     Output%OE_Qf => out2(Ialgo)%qcflg1
     Output%Packed_Qf => out2(Ialgo)%Acha_Packed_Quality_Flags
     Output%Packed_Meta_Data => out2(Ialgo)%Acha_Packed_Meta_Data_Flags
     Output%Processing_Order  => out2(Ialgo)%Processing_Order
     Output%Cost  => out2(Ialgo)%Acha_Cost
     Output%Pc_Opaque => out2(Ialgo)%Pc_Opaque_Cloud
     Output%Tc_Opaque => out2(Ialgo)%Tc_Opaque_Cloud
     Output%Zc_Opaque => out2(Ialgo)%Zc_Opaque_Cloud
     Output%Pc_H2O => out2(Ialgo)%Pc_H2O_Cloud
     Output%Tc_H2O => out2(Ialgo)%Tc_H2O_Cloud
     Output%Zc_H2O => out2(Ialgo)%Zc_H2O_Cloud

     !-------------------------------------------
     ! FIXME
     !----------------------------------------------
!---stw   Output%Cloud_Layer =>  null()
!---stw   Output%Total_Cloud_Fraction =>  null()
!---stw   Output%Total_Cloud_Fraction_Uncer =>  null()
!---stw   Output%High_Cloud_Fraction =>  null()
!---stw   Output%Mid_Cloud_Fraction =>  null()
!---stw   Output%Low_Cloud_Fraction =>  null()

 end subroutine SET_OUTPUT

 !-----------------------------------------------------------------------------
 !
 !-----------------------------------------------------------------------------
 subroutine SET_SYMBOL(Ialgo)

     implicit none
     integer(kind=int4), intent(in) :: Ialgo

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

 end subroutine SET_SYMBOL

end module ACHA_GEOCAT_BRIDGE_MODULE
