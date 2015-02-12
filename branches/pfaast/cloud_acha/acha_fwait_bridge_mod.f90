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
 private :: NULL_INPUT
 private :: NULL_OUTPUT
 private :: SET_INPUT
 private :: SET_OUTPUT
 private :: SET_SYMBOL
 private :: COVARIANCE_LOCAL
 private :: WMO_Sensor_KM
   
   
 contains

   REAL(SINGLE),  DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Covar_Ch27_Ch31_5x5

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

   integer(BYTE), dimension(:,:), POINTER  :: Cld_Mask
   integer(BYTE),dimension(:,:,:), POINTER  :: Cld_Test_Vector_Packed

   !--------------------------------------------------------------------
   ! define structures that will be arguments to ACHA
   !--------------------------------------------------------------------
   type(symbol_acha) :: symbol
   type(Input_struct) :: Input
   type(Output_struct) :: Output

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
  
   !----set symbols to local values
   call SET_SYMBOL()

   !-----------------------------------------------------------------------
   !--- Call to AWG CLoud Height Algorithm (ACHA)
   !-----------------------------------------------------------------------
   call AWG_CLOUD_HEIGHT_ALGORITHM(Input, &
                                    Symbol, &
                                    Output)

   !-----------------------------------------------------------------------
   !--- Call algorithm to make ACHA optical and microphysical properties
   !-----------------------------------------------------------------------
   call ACHA_COMP_ALGORITHM(Input, &
                            Symbol, &
                            Output)

   !-----------------------------------------------------------------------
   !--- Call to Geometrical Shadow Algorithm
   !-----------------------------------------------------------------------
   call CLOUD_SHADOW_RETR (  &
           Zc_Acha &
         , Solaz &
         , Solzen &
         , Lat &
         , Lon &
         , Lat_Pc &
         , Lon_Pc &
         , Shadow_Mask ) 

   !---- copy shadow result into cloud mask test bits
   
   !Set Cloud Mask local pointers
   CALL NFIA_CloudMask_Mask(Ctxt%CLOUD_MASK_Src1_T00, Cld_Mask)
   CALL NFIA_CloudMask_CldMaskPacked(Ctxt%CLOUD_MASK_Src1_T00, Cld_Test_Vector_Packed)

   
   where (Shadow_Mask == 1 .and. Cld_Mask == 0 )  
           Cld_Test_Vector_Packed ( 2 , :, : )  = ibset (Cld_Test_Vector_Packed ( 2 , :, : )  , 6 )
   end where
   
   Cld_Mask => null()
   NFIA_CloudMask_CldMaskPacked => null()

   !--- cloud cover layers
   call COMPUTE_CLOUD_COVER_LAYERS(Input,Symbol, Output)

   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT()
   call NULL_OUTPUT()                                 

 end subroutine AWG_CLOUD_HEIGHT_BRIDGE

 !-----------------------------------------------------------------------------
 !
 !-----------------------------------------------------------------------------


 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding input 
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
     Input%Covar_Bt_11um_67um =>  null()
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
     Input%Ctxt =>   null()
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
     
     if (allocated(Covar_Ch27_Ch31_5x5)) deallocate(Covar_Ch27_Ch31_5x5)
     
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

!--------------------------------------------------------
 subroutine SET_INPUT()
   type(AWG_CLX_CLOUD_HEIGHT_Ctxt) :: Ctxt
   integer(long) :: Stat
   integer :: WMO_Sc_Id, McIDAS_ID
   integer(byte), pointer, dimension(:):: Chan_On_Flag_Default
   !Variables for covariance calculation
   integer:: Num_Elem, Elem_Idx
   integer:: Num_Line, Line_Idx
   integer:: Num_Line_Max
   integer:: Elem_Idx_min
   integer:: Elem_Idx_max
   integer:: Elem_Idx_width
   integer:: Elem_Idx_segment_max
   integer:: Line_Idx_min
   integer:: Line_Idx_max
   integer:: Line_Idx_width
   integer:: Line_Idx_segment_max

   Input%Ctxt => Ctxt   

   Input%Number_of_Elements = Ctxt%SegmentInfo%Current_Column_Size
   Input%Number_of_Lines = Ctxt%SegmentInfo%Current_Row_Size
   Input%Num_Line_Max = Ctxt%SegmentInfo%Current_Row_Size
   Input%Process_Undetected_Cloud_Flag = sym%NO
   Input%Smooth_Nwp_Fields_Flag = sym%YES
   CALL NFIP_CloudHeight_AchaMode(Ctxt%CLOUD_HEIGHT_Src1_T00, Input%ACHA_Mode_Flag_In)
   
   ! McIDAS sensor ID
   CALL NFIP_Sat_Sat_ID(Ctxt%SATELLITE_DATA_Src1_T00, McIDAS_ID)
   CALL SatID_SSEC_to_WMO(McIDAS_ID, WMO_Sc_Id)
   
   Input%Sensor_Resolution_KM = WMO_Sensor_KM(WMO_Sc_Id)

   Input%Chan_Idx_67um = 9     !channel number for 6.7
   Input%Chan_Idx_85um = 11     !channel number for 8.5
   Input%Chan_Idx_11um = 14     !channel number for 11
   Input%Chan_Idx_12um = 15     !channel number for 12
   Input%Chan_Idx_133um = 16  !channel number for 13.3

    CALL NFIA_CloudHeight_ChanOn(Ctxt%CLOUD_HEIGHT_Src1_T00, Chan_On_Flag_Default)
     Input%Chan_On_67um = Chan_On_Flag_Default(Input%Chan_Idx_67um)
     Input%Chan_On_85um = Chan_On_Flag_Default(Input%Chan_Idx_85um)
     Input%Chan_On_11um = Chan_On_Flag_Default(Input%Chan_Idx_11um)
     Input%Chan_On_12um = Chan_On_Flag_Default(Input%Chan_Idx_12um)
     Input%Chan_On_133um = Chan_On_Flag_Default(Input%Chan_Idx_133um)

    CALL NFIA_Sat_L1b_BadPixMsk(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14, Input%Invalid_Data_Mask)
     
    CALL NFIA_NWP_X_NWP(Ctxt%NWP_DATA_Src1_T00, Input%Elem_Idx_NWP)
     
    CALL NFIA_NWP_Y_NWP(Ctxt%NWP_DATA_Src1_T00, Input%Line_Idx_NWP)

    CALL NFIA_NWP_X_NWP_Diag(Ctxt%NWP_DATA_Src1_T00, Input%Elem_Idx_Opposite_Corner_NWP)
    CALL NFIA_NWP_Y_NWP_Diag(Ctxt%NWP_DATA_Src1_T00, Input%Line_Idx_Opposite_Corner_NWP)
    CALL NFIA_NWP_LonWgtFac(Ctxt%NWP_DATA_Src1_T00, Input%Longitude_Interp_Weight_NWP)
    CALL NFIA_NWP_LatWgtFac(Ctxt%NWP_DATA_Src1_T00, Input%Latitude_Interp_Weight_NWP)

    CALL NFIA_RTM_ViewZenAngIndex(Ctxt%RTM_Src1_T00, Input%Viewing_Zenith_Angle_Idx_Rtm)   
   
     if (Chan_On_Flag_Default(Input%Chan_Idx_67um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI9,  Input%Bt_67um)
        CALL NFIA_Sat_L1b_Rad(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14,  Input%Rad_67um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI9,  Input%Rad_Clear_67um)
     endif

      if (Chan_On_Flag_Default(Input%Chan_Idx_85um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI11, Input%Bt_85um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI11, Input%Rad_Clear_85um)
      endif

     if (Chan_On_Flag_Default(Input%Chan_Idx_11um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14,  Input%Bt_11um)
        CALL NFIA_Sat_L1b_Rad(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14,  Input%Rad_11um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI14, Input%Rad_Clear_11um)
     endif

      if (Chan_On_Flag_Default(Input%Chan_Idx_12um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI15, Input%Bt_12um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI15, Input%Rad_Clear_12um)
      endif
          
      if (Chan_On_Flag_Default(Input%Chan_Idx_133um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI16, Input%Bt_133um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI16, Input%Rad_Clear_133um)
      endif
   
   ! ALLOCATE covariance array
   Num_Elem = Ctxt%SegmentInfo%Current_Column_Size
   Num_Line = Ctxt%SegmentInfo%Current_Row_Size
   
   
   ALLOCATE (Covar_Ch27_Ch31_5x5(Num_Elem,Num_Line))
   
   !COVARIANCE
   
    DO Elem_Idx = 1, Num_Elem
       DO Line_Idx = 1, Num_line

        !--- compute 5x5 arrays
        Elem_Idx_min = max(1,min(Elem_Idx - 2,Elem_Idx_segment_max))
        Elem_Idx_max = max(1,min(Elem_Idx + 2,Elem_Idx_segment_max))
        Line_Idx_min = max(1,min(Line_Idx - 2,Line_Idx_segment_max))
        Line_Idx_max = max(1,min(Line_Idx + 2,Line_Idx_segment_max))
        Line_Idx_width = Line_Idx_max - Line_Idx_min + 1
        Elem_Idx_width = Elem_Idx_max - Elem_Idx_min + 1

        IF ((Input%Chan_On_67um == sym%YES) .and. & 
            (Input%Chan_On_11um == sym%YES)) THEN

            Covar_Ch27_Ch31_5x5(Elem_Idx,Line_Idx) = Covariance_local(&
              Input%Bt_11um(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
              Input%Bt_67um(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               Elem_Idx_width, Line_Idx_width, &
               Input%Invalid_Data_Mask(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max))
        ENDIF

      ENDDO

    ENDDO


   
   
   Input%Covar_Bt_11um_67um => Covar_Ch27_Ch31_5x5

    CALL NFIA_Sat_Nav_COS_SatZen(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Cosine_Zenith_Angle)
     
    CALL NFIA_Sat_Nav_SatZen(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Sensor_Zenith_Angle)
     
     
    CALL NFIA_Sat_Nav_SatAzi(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Sensor_Azimuth_Angle)

    CALL NFIA_Sat_Nav_Lat(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Latitude)
     
    CALL NFIA_Sat_Nav_Lon(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Longitude)

    CALL NFIA_SnowMask_Mask(Ctxt%SNOW_MASK_Src1_T00, Input%Snow_Class)
     
    CALL NFIA_SfcType_Mask(Ctxt%SURFACE_TYPE_Src1_T00, Input%Surface_Type)

    CALL NFIP_NWP_TempSfc(Ctxt%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Input%Surface_Temperature)
    CALL NFIA_NWP_Temp2M(Ctxt%NWP_DATA_Src1_T00, Input%Surface_Air_Temperature)
    CALL NFIP_NWP_TempTropo(Ctxt%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Input%Tropopause_Temperature)
    CALL NFIA_NWP_PressSfc(Ctxt%NWP_DATA_Src1_T00, Input%Surface_Pressure)

    CALL NFIA_SfcElev_Elevation(Ctxt%SURFACE_ELEVATION_Src1_T00, Input%Surface_Elevation)
     
    CALL NFIA_CloudMask_Mask(Ctxt%CLOUD_MASK_Src1_T00, Input%Cloud_Mask)
     
    CALL NFIA_CloudPhase_CldType(Ctxt%CLOUD_PHASE_Src1_T00, Input%Cloud_Type)

   Input%Cloud_Probability => Posterior_Cld_Probability

   CALL NFIA_SfcEmis_SfcEmiss(Ctxt%SURFACE_EMISSIVITY_Src1_T00, CHN_ABI7, Input%Surface_Emissivity_39um)

   Input%Elem_Idx_LRC_Input => null()
   Input%Line_Idx_LRC_Input =>  null()
   Input%Tc_Cirrus_Sounder =>  null() ! need to fix in science code
 end subroutine SET_INPUT

!--------------------------------------------------------
 subroutine SET_OUTPUT()
   type(AWG_CLX_CLOUD_HEIGHT_Ctxt) :: Ctxt
   integer(long) :: Stat

    CALL NFIA_CloudHeight_Latitude_Pc(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Latitude_Pc)
    CALL NFIA_CloudHeight_Longitude_Pc(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Longitude_Pc)

    CALL NFIA_CloudHeight_CldTopTemp(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Tc)
     
    CALL NFIA_CloudHeight_CldTopEmss(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Ec)
     
    CALL NFIA_CloudHeight_CldBta1112(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Beta)
     
    CALL NFIA_CloudHeight_CldTopPres(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Pc)
     
    CALL NFIA_CloudHeight_CldTopHght(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Zc)
     
    CALL NFIA_CloudHeight_CldOptDpth(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Tau)
     
    CALL NFIA_CloudHeight_CldTopReff(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Reff)

    CALL NFIA_CloudHeight_TcError(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Tc_Uncertainty)
     
    CALL NFIA_CloudHeight_EcError(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Ec_Uncertainty)
     
    CALL NFIA_CloudHeight_B1112Error(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Beta_Uncertainty)
     
    CALL NFIA_CloudHeight_PcError(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Pc_Uncertainty)
     
    CALL NFIA_CloudHeight_ZcError(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Zc_Uncertainty)
     
    CALL NFIA_CloudHeight_TcLowerCld(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Lower_Cloud_Temperature)
          
    CALL NFIA_CloudHeight_PcLowerCld(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Lower_Cloud_Pressure)
     
    CALL NFIA_CloudHeight_ZcLowerCld(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Lower_Cloud_Height)

    CALL NFIA_CloudHeight_Zc_Top(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Zc_Top)
    CALL NFIA_CloudHeight_Zc_Base(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Zc_Base)

    CALL NFIA_CloudHeight_CldHgtQF(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Qf)
    CALL NFIA_CloudHeight_Flag(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%OE_Qf)
 
    CALL NFIA_CloudHeight_CldHgtQF(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Qf)

    CALL NFIA_CloudHeight_Flag(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%OE_Qf)
     
    CALL NFIA_CloudHeight_Packed_Qf(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Packed_Qf)
    CALL NFIA_CloudHeight_Packed_Mdata(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Packed_Meta_Data)
    CALL NFIA_CloudHeight_ProcOrder(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Processing_Order)

    CALL NFIA_CloudHeight_Cost(Ctxt%CLOUD_HEIGHT_Src1_T00, Output%Cost)
    
    CALL NFIA_CloudHeight_Cloud_Layer(Ctxt%CLOUD_HEIGHT_Src1_T00, &
                                      Output%Cloud_Layer)
    
   CALL NFIA_CloudHeight_Total_Cld_Frac(Ctxt%CLOUD_HEIGHT_Src1_T00,&
                                        Output%Total_Cloud_Fraction)
                                        
   CALL NFIA_CloudHeight_Total_Cld_Frac_Uncer(Ctxt%CLOUD_HEIGHT_Src1_T00, &
                                        Output%Total_Cloud_Fraction_Uncer)
                                        
   CALL NFIA_CloudHeight_High_Cld_Frac(Ctxt%CLOUD_HEIGHT_Src1_T00, &
                                        Output%High_Cloud_Fraction)
                                        
   CALL NFIA_CloudHeight_Mid_Cld_Frac(Ctxt%CLOUD_HEIGHT_Src1_T00, &
                                        Output%Mid_Cloud_Fraction)
                                        
   CALL NFIA_CloudHeight_Low_Cld_Frac(Ctxt%CLOUD_HEIGHT_Src1_T00, &
                                        Output%Low_Cloud_Fraction)
 end subroutine SET_OUTPUT


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
   function COVARIANCE_LOCAL &
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
   
   end function COVARIANCE_LOCAL


 !-----------------------------------------------------------------------------
 ! This function provides the sensor footprint
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
      case(224) ! S-NPP
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


end module ACHA_CLAVRX_BRIDGE_MOD
