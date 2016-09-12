!$Id$
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
module CLOUD_BASE_SAPF_BRIDGE

 use CLOUD_BASE
 use CLOUD_BASE_SERVICES
   
 implicit none

 public :: CLOUD_BASE_BRIDGE
 private :: NULL_INPUT_POINTERS
 private :: NULL_OUTPUT_POINTERS
 private :: SET_SYMBOL
 private :: SET_INPUT
 private :: SET_OUTPUT
 !private :: CLOUD_BASE_BRIDGE
 private :: CCL_LCL_CALC
 private :: VAPOR
 private :: VAPOR_ICE

 !--------------------------------------------------------------------
 ! define structures that will be arguments 
 !--------------------------------------------------------------------

 type(CLOUD_BASE_EN_Ctxt), POINTER, PRIVATE :: Ctxt_CBase

 type(Symbol_acha), private :: Symbol
 type(acha_input_struct), private :: Input
 type(acha_output_struct), private :: Output

 contains

!----------------------------------------------------------------------
! BRIDGE SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine CLOUD_BASE_BRIDGE(Ctxt, Stat)
 
   implicit none
   type(CLOUD_BASE_EN_Ctxt), target :: Ctxt
   integer(long) :: Stat


   Ctxt_CBase => Ctxt  

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
   !--- Call algorithm to make cloud geometrical boundaries
   !-----------------------------------------------------------------------
   call CLOUD_BASE_ALGORITHM(Input, &
                             Symbol, &
                             Output)

   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT_POINTERS()
   call NULL_OUTPUT_POINTERS()

   Ctxt_CBase => null()

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

     if (allocated(Input%CCL)) deallocate(Input%CCL)
     if (allocated(Input%LCL)) deallocate(Input%LCL)
     if (allocated(Input%CWP)) deallocate(Input%CWP)
     if (allocated(Input%CWP_nwp)) deallocate(Input%CWP_nwp)
     
     Input%Ctxt => null()

 end subroutine NULL_INPUT_POINTERS
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
   Symbol%PROB_CLEAR_TYPE = sym%CLEAR_TYPE
   Symbol%FOG_TYPE = sym%FOG_TYPE
   Symbol%WATER_TYPE = sym%WATER_TYPE
   Symbol%SUPERCOOLED_TYPE = sym%SUPERCOOLED_TYPE
   Symbol%MIXED_TYPE = sym%MIXED_TYPE
   Symbol%OPAQUE_ICE_TYPE = sym%TICE_TYPE
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



 subroutine SET_INPUT()

   integer(long) :: Stat
   INTEGER(KIND=INT4), pointer, DIMENSION(:) :: ChnMap
   integer(byte), dimension(16):: Chan_On_Flag_Default
   integer(byte) :: ABI_Chan

   !Variables for CCL/LCL/CWP calculation
   integer:: Num_Elem, Elem_Idx
   integer:: Num_Line, Line_Idx
   REAL(SINGLE), DIMENSION(:,:), POINTER :: lwp_dcomp, iwp_dcomp
   REAL(SINGLE) :: Tmpair_Nwp, Rhsfc_Nwp, lwp, iwp
   
   Num_line = Ctxt_CBase%SegmentInfo%Current_Row_Size
   Num_Elem = Ctxt_CBase%SegmentInfo%Current_Column_Size
   Input%Ctxt => Ctxt_CBase   

   Input%Number_of_Elements = Ctxt_CBase%SegmentInfo%Current_Column_Size
   Input%Number_of_Lines = Ctxt_CBase%SegmentInfo%Current_Row_Size
   Input%Num_Line_Max = Ctxt_CBase%SegmentInfo%Current_Row_Size
   
   Input%Chan_Idx_67um = 9     !channel number for 6.7
   Input%Chan_Idx_85um = 11     !channel number for 8.5
   Input%Chan_Idx_11um = 14     !channel number for 11
   Input%Chan_Idx_12um = 15     !channel number for 12
   Input%Chan_Idx_133um = 16  !channel number for 13.3

   Chan_On_Flag_Default(:)= 0
   
   CALL NFIA_Sat_ChnMap_Flag(Ctxt_CBase%SATELLITE_DATA_Src1_T00, ChnMap)
     
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

    CALL NFIA_Sat_L1b_BadPixMsk(Ctxt_CBase%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14, Input%Invalid_Data_Mask)
     
    CALL NFIA_NWP_X_NWP(Ctxt_CBase%NWP_DATA_Src1_T00, Input%Elem_Idx_NWP)
     
    CALL NFIA_NWP_Y_NWP(Ctxt_CBase%NWP_DATA_Src1_T00, Input%Line_Idx_NWP)

    CALL NFIA_NWP_X_NWP_Diag(Ctxt_CBase%NWP_DATA_Src1_T00, Input%Elem_Idx_Opposite_Corner_NWP)
    CALL NFIA_NWP_Y_NWP_Diag(Ctxt_CBase%NWP_DATA_Src1_T00, Input%Line_Idx_Opposite_Corner_NWP)
    CALL NFIA_NWP_LonWgtFac(Ctxt_CBase%NWP_DATA_Src1_T00, Input%Longitude_Interp_Weight_NWP)
    CALL NFIA_NWP_LatWgtFac(Ctxt_CBase%NWP_DATA_Src1_T00, Input%Latitude_Interp_Weight_NWP)

    CALL NFIA_RTM_ViewZenAngIndex(Ctxt_CBase%RTM_Src1_T00, Input%Viewing_Zenith_Angle_Idx_Rtm)   
   
    CALL NFIA_Sat_Nav_COS_SatZen(Ctxt_CBase%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Cosine_Zenith_Angle)
     
    CALL NFIA_Sat_Nav_SatZen(Ctxt_CBase%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Sensor_Zenith_Angle)

    CALL NFIA_Sat_Nav_Lat(Ctxt_CBase%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Latitude)     
    CALL NFIA_Sat_Nav_Lon(Ctxt_CBase%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Input%Longitude)

   CALL NFIA_SfcType_Mask(Ctxt_CBase%SURFACE_TYPE_Src1_T00, Input%Surface_Type)

    CALL NFIA_SfcElev_Elevation(Ctxt_CBase%SURFACE_ELEVATION_Src1_T00, Input%Surface_Elevation)
     
    CALL NFIA_CloudMask_Mask(Ctxt_CBase%CLOUD_MASK_Src1_T00, Input%Cloud_Mask)
     
    CALL NFIA_CloudPhase_CldType(Ctxt_CBase%CLOUD_PHASE_Src1_T00, Input%Cloud_Type)
    
    
    CALL NFIA_CloudMask_CldProbability(Ctxt_CBase%CLOUD_MASK_Src1_T00, Input%Cloud_Probability)
    
    
   CALL NFIA_CloudHeight_CldTopTemp(Ctxt_CBase%CLOUD_HEIGHT_Src1_T00, Input%Tc)
   CALL NFIA_CloudHeight_CldTopPres(Ctxt_CBase%CLOUD_HEIGHT_Src1_T00, Input%Pc)
   CALL NFIA_CloudHeight_CldTopHght(Ctxt_CBase%CLOUD_HEIGHT_Src1_T00, Input%Zc)
   CALL NFIA_CloudHeight_CldOptDpth(Ctxt_CBase%CLOUD_HEIGHT_Src1_T00, Input%Tau)
     

   !now for the vars that need to be allocated and determined
   
   ALLOCATE (Input%LCL(Num_Elem,Num_Line)) !LCL_Height_Nwp_Pix
   ALLOCATE (Input%CCL(Num_Elem,Num_Line)) !CCL_Height_Nwp_Pix
   ALLOCATE (Input%CWP(Num_Elem,Num_Line)) ! Cwp_Dcomp
   ALLOCATE (Input%CWP_nwp(Num_Elem,Num_Line)) !Cwp_Nwp_Pix
   
   !Now to do the determination of pixel level variables
   ! For CWP, we need LWP and IWP from DCOMP
   
   CALL NFIA_CloudMicro_CldLWP(   Ctxt_CBase%CLOUD_MICRO_Src1_T00, lwp_dcomp)
   CALL NFIA_CloudMicro_CldIWP(   Ctxt_CBase%CLOUD_MICRO_Src1_T00, iwp_dcomp)
   
   
    DO Elem_Idx = 1, Num_Elem
       DO Line_Idx = 1, Num_line
       
        !initalize variables
        Input%LCL(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Input%CCL(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Input%CWP(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Input%CWP_nwp(Elem_Idx,Line_Idx) = Missing_Value_Real4
        
        ! calculate DCOMP CWP
        IF ((lwp_dcomp(Elem_Idx,Line_Idx) /= Missing_Value_Real4) .OR. &
            (iwp_dcomp(Elem_Idx,Line_Idx) /= Missing_Value_Real4)) THEN
            
            !ensure values are 0 for missing pixels
            IF(lwp_dcomp(Elem_Idx,Line_Idx) .EQ. Missing_Value_Real4) lwp = 0
            IF(iwp_dcomp(Elem_Idx,Line_Idx) .EQ. Missing_Value_Real4) iwp = 0
            
            Input%CWP(Elem_Idx,Line_Idx) = iwp + lwp
        
        ENDIF
        
        
        !Fill NWP level values for CCL/LCL calculation and NWP CWP
        CALL NFIP_NWP_Temp2M(Ctxt_CBase%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Tmpair_Nwp)
        
        CALL NFIP_NWP_RH_2m(Ctxt_CBase%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Rhsfc_Nwp)
        

        !COMMENTED OUT UNTIL NWP_CWP IS INTEGERATED IN - WCS3
        
        !CALL NFIP_NWP_CWP(Ctxt_CBase%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Input%CWP_nwp(Elem_Idx,Line_Idx)) ! Need to put in SAPF

        !Determine CCL/LCL
        CALL  CCL_LCL_CALC(Tmpair_Nwp, Rhsfc_Nwp, Input%LCL(Elem_Idx,Line_Idx), &
                           Input%CCL(Elem_Idx,Line_Idx))
        
       ENDDO

    ENDDO
  
   
   
   !Variables not used, but in input structure
   Input%Elem_Idx_LRC_Input => null()
   Input%Line_Idx_LRC_Input =>  null()
   Input%Latitude_Pc => null()
   Input%Longitude_Pc => null()
   Input%Ec => null()
   Input%Beta => null()
   Input%Reff => null()
   Input%Tc_Uncertainty => null()
   Input%Ec_Uncertainty => null()
   Input%Beta_Uncertainty => null()
   Input%Pc_Uncertainty => null()
   Input%Zc_Uncertainty => null()
   Input%Lower_Cloud_Pressure => null()
   Input%Lower_Cloud_Temperature => null()
   Input%Lower_Cloud_Height => null()
   Input%Cdnc => null()
   Input%Hcld => null()

   
 end subroutine SET_INPUT


 subroutine SET_OUTPUT()

    CALL NFIA_CloudBase_CldTopHght(Ctxt_CBase%CLOUD_BASE_Src1_T00, Output%Zc_Top)
    CALL NFIA_CloudBase_CldBaseHght(Ctxt_CBase%CLOUD_BASE_Src1_T00, Output%Zc_Base)

    CALL NFIA_CloudBase_CldTopPres(Ctxt_CBase%CLOUD_BASE_Src1_T00, Output%Pc_Top)
     
    CALL NFIA_CloudBase_CldBasePres(Ctxt_CBase%CLOUD_BASE_Src1_T00, Output%Pc_Base)

    CALL NFIA_CloudBase_CldBaseQF(Ctxt_CBase%CLOUD_BASE_Src1_T00, Output%Zc_Base_Qf)
    
    !initialize DQF
    Output%Zc_Base_Qf = 1
    

 end subroutine SET_OUTPUT


 subroutine CCL_LCL_CALC(Temp, RH, LCL, CCL)
    REAL(SINGLE), INTENT(IN) :: Temp, RH   
    REAL(SINGLE), INTENT(OUT) :: LCL, CCL
    REAL(SINGLE) :: e, es, Td
    
    if (Temp > 180.0) then
        if (Temp > 253.0) then
            es = VAPOR(Temp)    !saturation vapor pressure wrt water hpa
        else
            es = VAPOR_ICE(Temp) !saturation vapor pressure wrt ice hpa
        endif
        e = es * RH / 100.0  !vapor pressure in hPa
        Td = 273.15 + 243.5 * alog(e / 6.112)  / (17.67 - alog(e/6.112))     !Dewpoint T in K
        LCL = 1000. * 0.125*(Temp - Td)  ! meters 
        CCL = 1000.* (Temp - Td)/4.4  ! meters 
    end if
    
 end subroutine CCL_LCL_CALC

!----------------------------------------------------------------
! functions to compute some needed water vapor parameters
!----------------------------------------------------------------
 function VAPOR(T) result(es)
                                                                     
!  T in Kelvin                                                          
!  es in mbar

  implicit none
  real, intent (in) :: T
  real :: es

   es = 6.112 * exp(17.67 * (T-273.16) / (T - 29.66))

  return 
end function VAPOR

!---- saturation vapor pressure for ice
function VAPOR_ICE(T) result(es)
   implicit none
   real, intent(in):: T
   real:: es
     es = 6.1078 * exp(21.8745584 * (T-273.16) / (T - 7.66))
  return 
end function VAPOR_ICE


end module CLOUD_BASE_SAPF_BRIDGE
