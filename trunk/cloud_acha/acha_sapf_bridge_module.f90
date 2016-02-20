!$Id$
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Height Algorithm (ACHA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  processing systems and the ACHA code.
!
!------------------------------------------------------------------------------
module ACHA_CLAVRX_BRIDGE_MOD

 use ACHA_SERVICES_MOD
 use AWG_CLOUD_HEIGHT_ACHA
 use ACHA_CLOUD_COVER_LAYERS
 use ACHA_SHADOW
 use ACHA_COMP 
   
 implicit none

 public :: AWG_CLOUD_HEIGHT_BRIDGE
 private :: NULL_INPUT
 private :: NULL_OUTPUT
 private :: SET_INPUT
 private :: SET_OUTPUT
 private :: SET_SYMBOL
 private :: COVARIANCE_LOCAL
 private :: WMO_Sensor_KM
   
 REAL(SINGLE),  DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Covar_Ch27_Ch31_5x5
 REAL(SINGLE),  DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Dummy 
 integer(kind=int1),  DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Dummy1 

 type(NPP_VIIRS_CLD_HEIGHT_Ctxt), POINTER, PRIVATE :: Ctxt_ACHA

 integer, parameter, private :: N_Sc_Lut = 20



 type(acha_symbol_struct), PRIVATE  :: symbol
 type(acha_input_struct), PRIVATE  :: Input
 type(acha_output_struct), PRIVATE  :: Output
 
 contains


! Modes
! 0 - Use this mode to not call ACHA from the framework
! 1 - 11 um                          0           
! 2 - 11 + 6.7 um                    7
! 3 - 11 + 12 um                     1
! 4 - 11 + 13.3 um                   2
! 5 - 11 + 8.5 + 12 um               4
! 6 - 11 + 6.7 + 12 um               5
! 7 - 11 + 6.7 + 13.3 um             6
! 8 - 11 + 12 + 13.3 um              3

!----------------------------------------------------------------------
! BEGINNING OF ACHA SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine AWG_CLOUD_HEIGHT_BRIDGE(Ctxt, Stat)

   implicit none
   type(NPP_VIIRS_CLD_HEIGHT_Ctxt), target :: Ctxt
   integer(long) :: Stat
   integer:: Num_Elem, Num_Line

   integer(BYTE), dimension(:,:), POINTER  :: Cld_Mask
   integer(BYTE),dimension(:,:,:), POINTER  :: Cld_Test_Vector_Packed
   integer(BYTE), dimension(:,:), POINTER  :: Shadow_Mask
   REAL(SINGLE), DIMENSION(:,:), POINTER :: SolZen
   REAL(SINGLE), DIMENSION(:,:), POINTER :: SolAz
   
   REAL(SINGLE), DIMENSION(:,:), POINTER :: Conv_Cld_Prob
   REAL(SINGLE), DIMENSION(:,:), POINTER :: Supercooled_Cld_Prob
    
   Ctxt_ACHA => Ctxt

   !--------------------------------------------------------------------
   ! define structures that will be arguments to ACHA
   !--------------------------------------------------------------------

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
   
   call AWG_CLOUD_HEIGHT_ACHA_ALGORITHM(Input, &
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
   !ALLOCATE Shadow Mask
   Num_Elem = Ctxt%SegmentInfo%Current_Column_Size
   Num_Line = Ctxt%SegmentInfo%Current_Row_Size
   
   CALL NFIA_CloudHeight_Shadow_Mask(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, &
                                        Shadow_Mask)
   

   CALL NFIA_Sat_Nav_SolZen(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, SolZen)
   
   CALL NFIA_Sat_Nav_SolAzi(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, SolAz)
   
   call CLOUD_SHADOW_RETR (  &
           Output%Zc &
         , Solaz &
         , Solzen &
         , Input%Latitude &
         , Input%Longitude &
         , Output%Latitude_Pc &
         , Output%Longitude_Pc &
         , Shadow_Mask ) 
   Shadow_Mask => null()
    

   !--- cloud cover layers
   call COMPUTE_CLOUD_COVER_LAYERS(Input,Symbol, Output)  
   
   !--- Convective and supercooled cloud probability
   
   CALL NFIA_CloudHeight_SC_Cld_Prob(Ctxt%CLOUD_HEIGHT_Src1_T00, Supercooled_Cld_Prob)
   CALL NFIA_CloudHeight_Conv_Cld_Prob(Ctxt%CLOUD_HEIGHT_Src1_T00, Conv_Cld_Prob)

   CALL CONVECTIVE_CLOUD_PROBABILITY(Input, Output, Conv_Cld_Prob)

   CALL SUPERCOOLED_CLOUD_PROBABILITY(Input, Output, &
                                      Supercooled_Cld_Prob)
   
    
   Conv_Cld_Prob => null()
   Supercooled_Cld_Prob => null()

    

   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT()
   call NULL_OUTPUT()
    
   Ctxt_ACHA => null()                            
   
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


     if (allocated(Covar_Ch27_Ch31_5x5)) deallocate(Covar_Ch27_Ch31_5x5)
     
     !Until surface level stuff in framework, comment this out
     !Input%Surface_Temperature => null()
     !Input%Surface_Air_Temperature =>  null()
     !Input%Tropopause_Temperature =>  null()
     !Input%Surface_Pressure =>  null()
     
     !And then do this
   if (allocated(Input%Surface_Temperature)) deallocate(Input%Surface_Temperature)

   if (allocated(Input%Surface_Air_Temperature)) deallocate(Input%Surface_Air_Temperature)

   if (allocated(Input%Tropopause_Temperature)) deallocate(Input%Tropopause_Temperature)


   if (allocated(Input%Surface_Pressure)) deallocate(Input%Surface_Pressure)


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
     Output%Inversion_Flag  =>  null()

     Output%Pc_Opaque =>  null()

     !rchen  05/29/2015
     !Output%Tc_Opaque =>  null()
     !Output%Zc_Opaque =>  null()
     !Output%Pc_H2O =>  null()
     !Output%Tc_H2O =>  null()
     !Output%Zc_H2O =>  null()
!     if (allocated(Output%Pc_Opaque)) deallocate(Output%Pc_Opaque)
     if (allocated(Output%Tc_Opaque)) deallocate(Output%Tc_Opaque)
     if (allocated(Output%Zc_Opaque)) deallocate(Output%Zc_Opaque)
     if (allocated(Output%Pc_H2O)) deallocate(Output%Pc_H2O)
     if (allocated(Output%Tc_H2O)) deallocate(Output%Tc_H2O)
     if (allocated(Output%Zc_H2O)) deallocate(Output%Zc_H2O)
     
     
    !ASOS output - Nulled and turned off for now - WCS 20Feb, 2016

    if (allocated(Dummy)) deallocate(Dummy)
    if (allocated(Dummy1)) deallocate(Dummy1)
    Output%ASOS_Cloud_Code => null()
    Output%ASOS_Cloud_ECA =>  null()
    Output%ASOS_Cloud_Zmin =>  null()
    Output%ASOS_Cloud_Zmax =>  null()

     
     
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
   Symbol%PROB_CLEAR_TYPE = sym%CLEAR_TYPE !rchen use clear_type for prob_clear_type
   Symbol%FOG_TYPE = sym%FOG_TYPE
   Symbol%WATER_TYPE = sym%WATER_TYPE
   Symbol%SUPERCOOLED_TYPE = sym%SUPERCOOLED_TYPE
   Symbol%MIXED_TYPE = sym%MIXED_TYPE
   Symbol%OPAQUE_ICE_TYPE = sym%TICE_TYPE  !rchen use tice_type for OPAQUE_ICE_TYPE
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
   integer(long) :: Stat
   integer(KIND=INT4) :: WMO_Sc_Id, McIDAS_ID
!   integer(byte), pointer, dimension(:):: Chan_On_Flag_Default
   INTEGER(KIND=INT4), pointer, DIMENSION(:) :: ChnMap
   integer(byte), dimension(16):: Chan_On_Flag_Default
   integer(byte) :: ABI_Chan
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

   !Set up input structure
   Input%Ctxt => Ctxt_ACHA   

   Input%Number_of_Elements = Ctxt_ACHA%SegmentInfo%Current_Column_Size
   Input%Number_of_Lines = Ctxt_ACHA%SegmentInfo%Current_Row_Size
   Input%Num_Line_Max = Ctxt_ACHA%SegmentInfo%Current_Row_Size
   Input%Process_Undetected_Cloud_Flag = sym%NO
   Input%Smooth_Nwp_Fields_Flag = sym%YES
   CALL NFIP_CloudHeight_AchaMode(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Input%ACHA_Mode_Flag_In)
   
   ! McIDAS sensor ID
   CALL NFIP_Sat_Sat_ID(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, McIDAS_ID)
   
   CALL SatID_SSEC_to_WMO(McIDAS_ID, WMO_Sc_Id)
   
   Input%Sensor_Resolution_KM = WMO_Sensor_KM(WMO_Sc_Id)

   Input%Chan_Idx_67um = 9     !channel number for 6.7
   Input%Chan_Idx_85um = 11     !channel number for 8.5
   Input%Chan_Idx_11um = 14     !channel number for 11
   Input%Chan_Idx_12um = 15     !channel number for 12
   Input%Chan_Idx_133um = 16  !channel number for 13.3

   Chan_On_Flag_Default(:)= 0
   
   CALL NFIA_Sat_ChnMap_Flag(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, ChnMap)
     
    DO ABI_Chan = 1, 16
        IF (ChnMap(ABI_Chan) /= 0) THEN
            Chan_On_Flag_Default(ABI_Chan) = 1
        ENDIF   
     
    ENDDO
    
     Input%Chan_On_67um = Chan_On_Flag_Default(Input%Chan_Idx_67um)
     Input%Chan_On_85um = Chan_On_Flag_Default(Input%Chan_Idx_85um)
     Input%Chan_On_11um = Chan_On_Flag_Default(Input%Chan_Idx_11um)
     Input%Chan_On_12um = Chan_On_Flag_Default(Input%Chan_Idx_12um)
     Input%Chan_On_133um = Chan_On_Flag_Default(Input%Chan_Idx_133um)
     
    CALL NFIA_Sat_L1b_BadPixMsk(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14, Input%Invalid_Data_Mask)
     
    CALL NFIA_NWP_X_NWP(Ctxt_ACHA%NWP_DATA_Src1_T00, Input%Elem_Idx_NWP)
     
    CALL NFIA_NWP_Y_NWP(Ctxt_ACHA%NWP_DATA_Src1_T00, Input%Line_Idx_NWP)

    CALL NFIA_NWP_X_NWP_Diag(Ctxt_ACHA%NWP_DATA_Src1_T00, Input%Elem_Idx_Opposite_Corner_NWP)
    CALL NFIA_NWP_Y_NWP_Diag(Ctxt_ACHA%NWP_DATA_Src1_T00, Input%Line_Idx_Opposite_Corner_NWP)
    CALL NFIA_NWP_LonWgtFac(Ctxt_ACHA%NWP_DATA_Src1_T00, Input%Longitude_Interp_Weight_NWP)
    CALL NFIA_NWP_LatWgtFac(Ctxt_ACHA%NWP_DATA_Src1_T00, Input%Latitude_Interp_Weight_NWP)

    CALL NFIA_RTM_ViewZenAngIndex(Ctxt_ACHA%RTM_Src1_T00, Input%Viewing_Zenith_Angle_Idx_Rtm)   
   
     if (Chan_On_Flag_Default(Input%Chan_Idx_67um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI9,  Input%Bt_67um)
        CALL NFIA_Sat_L1b_Rad(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI9,  Input%Rad_67um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt_ACHA%RTM_Src1_T00, CHN_ABI9,  Input%Rad_Clear_67um)
     endif

      if (Chan_On_Flag_Default(Input%Chan_Idx_85um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI11, Input%Bt_85um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt_ACHA%RTM_Src1_T00, CHN_ABI11, Input%Rad_Clear_85um)
      endif

     if (Chan_On_Flag_Default(Input%Chan_Idx_11um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14,  Input%Bt_11um)
        CALL NFIA_Sat_L1b_Rad(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14,  Input%Rad_11um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt_ACHA%RTM_Src1_T00, CHN_ABI14, Input%Rad_Clear_11um)
     endif

      if (Chan_On_Flag_Default(Input%Chan_Idx_12um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI15, Input%Bt_12um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt_ACHA%RTM_Src1_T00, CHN_ABI15, Input%Rad_Clear_12um)
      endif
          
      if (Chan_On_Flag_Default(Input%Chan_Idx_133um) == sym%YES) then
        CALL NFIA_Sat_L1b_BrtTemp(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI16, Input%Bt_133um)
        CALL NFIA_RTM_Pixel_RadClr(Ctxt_ACHA%RTM_Src1_T00, CHN_ABI16, Input%Rad_Clear_133um)
      endif
   
   ! ALLOCATE covariance array
   Num_Elem = Ctxt_ACHA%SegmentInfo%Current_Column_Size
   Num_Line = Ctxt_ACHA%SegmentInfo%Current_Row_Size
   Elem_Idx_segment_max = Num_Elem
   Line_Idx_segment_max = Num_line
   
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

    CALL NFIA_Sat_Nav_COS_SatZen(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Cosine_Zenith_Angle)
     
    CALL NFIA_Sat_Nav_SatZen(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Sensor_Zenith_Angle)
     
     
    CALL NFIA_Sat_Nav_SatAzi(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Sensor_Azimuth_Angle)

    CALL NFIA_Sat_Nav_Lat(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Latitude)
     
    CALL NFIA_Sat_Nav_Lon(Ctxt_ACHA%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Longitude)

    CALL NFIA_SnowMask_Mask(Ctxt_ACHA%SNOW_MASK_Src1_T00, Input%Snow_Class)
     
    CALL NFIA_SfcType_Mask(Ctxt_ACHA%SURFACE_TYPE_Src1_T00, Input%Surface_Type)


!Until pixel level data is available, comment this out
!    CALL NFIP_NWP_TempSfc(Ctxt_ACHA%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Input%Surface_Temperature)
!    CALL NFIA_NWP_Temp2M(Ctxt_ACHA%NWP_DATA_Src1_T00, Input%Surface_Air_Temperature)
!    CALL NFIP_NWP_TempTropo(Ctxt_ACHA%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Input%Tropopause_Temperature)
!    CALL NFIA_NWP_PressSfc(Ctxt_ACHA%NWP_DATA_Src1_T00, Input%Surface_Pressure)
!
! And then do this

    allocate(Input%Surface_Temperature(Input%Number_of_Elements,Input%Number_of_Lines ))
    allocate(Input%Surface_Air_Temperature(Input%Number_of_Elements,Input%Number_of_Lines ))
    allocate(Input%Tropopause_Temperature(Input%Number_of_Elements,Input%Number_of_Lines ))
    allocate(Input%Surface_Pressure(Input%Number_of_Elements,Input%Number_of_Lines ))

      !Fill the pixel level NWP stuff for now
      CALL ACHA_NWP_Fill(Input)

    CALL NFIA_SfcElev_Elevation(Ctxt_ACHA%SURFACE_ELEVATION_Src1_T00, Input%Surface_Elevation)
     
    CALL NFIA_CloudMask_Mask(Ctxt_ACHA%CLOUD_MASK_Src1_T00, Input%Cloud_Mask)
     
    CALL NFIA_CloudPhase_CldType(Ctxt_ACHA%CLOUD_PHASE_Src1_T00, Input%Cloud_Type)
    
    
    CALL NFIA_CloudMask_CldProbability(Ctxt_ACHA%CLOUD_MASK_Src1_T00, Input%Cloud_Probability)

   CALL NFIA_SfcEmis_SfcEmiss(Ctxt_ACHA%SURFACE_EMISSIVITY_Src1_T00, CHN_ABI7, Input%Surface_Emissivity_39um)

   Input%Elem_Idx_LRC_Input => null()
   Input%Line_Idx_LRC_Input =>  null()
   Input%Tc_Cirrus_Sounder =>  null()

 end subroutine SET_INPUT

!--------------------------------------------------------
 subroutine SET_OUTPUT()
   integer(long) :: Stat
   integer:: Num_Elem, Num_Line


    CALL NFIA_CloudHeight_Latitude_Pc(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Latitude_Pc)
    CALL NFIA_CloudHeight_Longitude_Pc(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Longitude_Pc)

    CALL NFIA_CloudHeight_CldTopTemp(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Tc)
     
    CALL NFIA_CloudHeight_CldTopEmss(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Ec)
     
    CALL NFIA_CloudHeight_CldBta1112(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Beta)
     
    CALL NFIA_CloudHeight_CldTopPres(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Pc)
     
    CALL NFIA_CloudHeight_CldTopHght(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Zc)
     
    CALL NFIA_CloudHeight_CldOptDpth(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Tau)
     
    CALL NFIA_CloudHeight_CldTopReff(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Reff)

    CALL NFIA_CloudHeight_TcError(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Tc_Uncertainty)
     
    CALL NFIA_CloudHeight_EcError(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Ec_Uncertainty)
     
    CALL NFIA_CloudHeight_B1112Error(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Beta_Uncertainty)
     
    CALL NFIA_CloudHeight_PcError(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Pc_Uncertainty)
     
    CALL NFIA_CloudHeight_ZcError(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Zc_Uncertainty)
     
    CALL NFIA_CloudHeight_TcLowerCld(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Lower_Cloud_Temperature)
          
    CALL NFIA_CloudHeight_PcLowerCld(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Lower_Cloud_Pressure)
     
    CALL NFIA_CloudHeight_ZcLowerCld(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Lower_Cloud_Height)

    CALL NFIA_CloudHeight_CldHgtQF(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Qf)
    CALL NFIA_CloudHeight_Flag(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%OE_Qf)
 
    CALL NFIA_CloudHeight_CldHgtQF(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Qf)

    CALL NFIA_CloudHeight_Flag(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%OE_Qf)
     
    CALL NFIA_CloudHeight_Packed_Qf(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Packed_Qf)
    CALL NFIA_CloudHeight_Packed_Mdata(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Packed_Meta_Data)
    CALL NFIA_CloudHeight_ProcOrder(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Processing_Order)

    CALL NFIA_CloudHeight_Cost(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, Output%Cost)
    
    CALL NFIA_CloudHeight_CldLayer(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, &
                                      Output%Cloud_Layer)
        
   CALL NFIA_CloudHeight_Total_Cld_Frac(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00,&
                                        Output%Total_Cloud_Fraction)
                                        
   CALL NFIA_CloudHeight_TotCldFracUnc(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, &
                                        Output%Total_Cloud_Fraction_Uncer)
                                        
   CALL NFIA_CloudHeight_High_Cld_Frac(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, &
                                        Output%High_Cloud_Fraction)
                                        
   CALL NFIA_CloudHeight_Mid_Cld_Frac(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, &
                                        Output%Mid_Cloud_Fraction)
                                        
   CALL NFIA_CloudHeight_Low_Cld_Frac(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, &
                                        Output%Low_Cloud_Fraction)

   !Inversion flag
   CALL NFIA_CloudHeight_InverFlag(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, &
                                    Output%Inversion_Flag)
   
   
   ! ALLOCATE Dummy array
   Num_Elem = Ctxt_ACHA%SegmentInfo%Current_Column_Size
   Num_Line = Ctxt_ACHA%SegmentInfo%Current_Row_Size

   CALL NFIA_CloudHeight_Pc_Opaque(Ctxt_ACHA%CLOUD_HEIGHT_Src1_T00, &
                                    Output%Pc_Opaque)


!   ALLOCATE (Output%Pc_Opaque(Num_Elem,Num_Line))
   ALLOCATE (Output%Tc_Opaque(Num_Elem,Num_Line))
   ALLOCATE (Output%Zc_Opaque(Num_Elem,Num_Line))
   ALLOCATE (Output%Pc_H2O(Num_Elem,Num_Line))
   ALLOCATE (Output%Tc_H2O(Num_Elem,Num_Line))
   ALLOCATE (Output%Zc_H2O(Num_Elem,Num_Line))

   !ASOS output - Pointed to dummy arrays for now, since ASOS is turned off
   !              May need to put into XML at a later point pending 
   !              ASOS project - WCS 20 Feb, 2016
   
   ALLOCATE (Dummy(Num_Elem,Num_Line))
   ALLOCATE (Dummy1(Num_Elem,Num_Line))
   Output%ASOS_Cloud_Code => Dummy1
   Output%ASOS_Cloud_ECA =>  Dummy
   Output%ASOS_Cloud_Zmin =>  Dummy
   Output%ASOS_Cloud_Zmax =>  Dummy
   
   

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
      case(173) ! Himawari-8
             Sensor_KM = 2.0	     
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


!-------------------------------------------------------------------------------
! Convective Cloud Probability
!-------------------------------------------------------------------------------
subroutine CONVECTIVE_CLOUD_PROBABILITY(Input, Output, Conv_Cld_Prob)
  type(acha_input_struct), intent(inout) :: Input
  type(acha_output_struct), intent(inout) :: Output
  real(kind=real4), dimension(:,:), intent(out):: Conv_Cld_Prob
  real(kind=real4), dimension(:,:), pointer:: Emiss_11_Tropo

  real, parameter:: Btd_Thresh = -2.0
  real, parameter:: Etrop_Thresh_1 = 0.95
  real, parameter:: Etrop_Thresh_2 = 0.90
  real, parameter:: Tsfc_Thresh = 30.0
  
  !since CM already computes 11um Tropopause emissivity, just use it here
  Emiss_11_Tropo => null()
  CALL NFIA_CloudMask_Emiss11High(Ctxt_ACHA%CLOUD_MASK_Src1_T00, Emiss_11_Tropo)

  

  !--- initialize to 0 (no)
  Conv_Cld_Prob = 0.0

  !--- set bad data to missing
  where(Input%Invalid_Data_Mask ==sym%YES)
    Conv_Cld_Prob = Missing_Value_Real4
  endwhere

  !--- if only 11 micron, use a tight threshold
  if (Input%Chan_On_11um == sym%YES) then
    where(Emiss_11_Tropo >= Etrop_Thresh_1) 
      Conv_Cld_Prob = 1.0
    endwhere
  endif

  if (Input%Chan_On_67um == sym%YES .and. &
      Input%Chan_On_11um == sym%YES) then
    where(Emiss_11_Tropo >= Etrop_Thresh_2 .and. (Input%Bt_67um - Input%Bt_11um) >= Etrop_Thresh_2) 
     Conv_Cld_Prob = 1.0
    endwhere
  endif

  !--- limit false alarms over elevated terrain
  if (Input%Chan_On_11um == sym%YES) then
    where(Input%Surface_Temperature - Input%Bt_11um < Tsfc_Thresh)
       Conv_Cld_Prob = 0.0
    endwhere
  endif
  
  
  
  Emiss_11_Tropo => null()


end subroutine CONVECTIVE_CLOUD_PROBABILITY

!-------------------------------------------------------------------------------
! Supercooled Cloud Probability
!-------------------------------------------------------------------------------
subroutine SUPERCOOLED_CLOUD_PROBABILITY(Input, Output,&
                                         Supercooled_Cld_Prob)
  type(acha_input_struct), intent(inout) :: Input
  type(acha_output_struct), intent(inout) :: Output
  real(kind=real4), dimension(:,:), intent(out):: Supercooled_Cld_Prob

  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Number_Of_Elements
  integer:: Number_Of_Lines
  integer:: Tc_Idx
  
  !Supercooled LUT
  real(kind=real4), dimension(N_Sc_Lut):: &
  Sc_Tc_Lut = (/202.21,206.77,211.33,215.88,220.44,225.00,229.56,234.11,238.67,243.23, &
                247.79,252.34,256.90,261.46,266.01,270.57,275.13,279.69,284.24,288.80/)

  real(kind=real4), dimension(N_Sc_Lut):: &
  Sc_Prob_Lut = (/0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.002,0.051,0.243, &
                  0.470,0.658,0.800,0.886,0.888,0.870,0.142,0.004,0.004,0.004 /)

  
  

  !--- intialize local variables using global variables
  Number_Of_Elements = Input%Number_Of_Elements
  Number_Of_Lines = Input%Num_Line_Max

  !--- initialize to Missing
  Supercooled_Cld_Prob = Missing_Value_Real4

  !--- Loop through each pixel
  do Elem_Idx = 1, Number_Of_Elements 
     do Line_Idx = 1, Number_Of_Lines
        
       !--- filter bad
       if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
       !--- filter missing temps
       if (Output%Tc(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
       !--- filter with cloud type
       if (Input%Cloud_Type(Elem_Idx,Line_Idx) == Missing_Value_Int1) cycle

       if (Input%Cloud_Type(Elem_Idx,Line_Idx) /= sym%FOG_TYPE .and. &
           Input%Cloud_Type(Elem_Idx,Line_Idx) /= sym%WATER_TYPE .and. &
           Input%Cloud_Type(Elem_Idx,Line_Idx) /= sym%SUPERCOOLED_TYPE) then
           Supercooled_Cld_Prob(Elem_Idx,Line_Idx) = 0.0
           cycle
       endif

       Tc_Idx = minloc(abs(Output%Tc(Elem_Idx,Line_Idx) - Sc_Tc_Lut),1) 

       Tc_Idx = max(1,min(N_Sc_Lut,Tc_Idx))

       Supercooled_Cld_Prob(Elem_Idx,Line_Idx) = Sc_Prob_Lut(Tc_Idx) 

     enddo
  enddo

end subroutine SUPERCOOLED_CLOUD_PROBABILITY




end module ACHA_CLAVRX_BRIDGE_MOD
