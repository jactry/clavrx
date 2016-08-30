MODULE Baseline_Cloud_Mask
!$Id: baseline_cloud_mask.f90,v 1.60 2012/11/27 00:03:49 ccalvert Exp $
!
! Module Name:
!   baseline_cloud_mask
!
! Function:
!   Contain the subroutines involved in the ABI cloud Mask
!
! Description:
!    This module holds the routines needed to compute the
!    ABI cloud Mask (ACM).  This is the upfront cloud Mask
!    that will be made available to all susequent ABI algorithms
!
! Reference:
!    The ABI Cloud Mask (ACM) Algorithm Theoretical Basis Document
!
! Dependencies:
!    Algorithm_Module Usage Module
!    Planck Module
!
! Calling Sequence:
!     Use baseline_cloud_mask
!
! Author: Andrew K. Heidinger, National Oceanic and Atmospheric Administration
! (c) This code is copyrighted by the author and all NOAA restrictions apply
!
! History:
!   2/2007 - Andrew Heidinger - Created based on CLAVR-x
!   7/2007 - Added Temporal IR Test
!   10/2008 - Removed ocean reflectance lookup tables
!   10/2008 - Added Dark Composite Chn2 Fields
!   10/2008 - Added MSG Gridded MODIS White Sky Chn2 Fields
!   10/2008 - Added mapped MODIS White Sky Chn2 and Chn5 Fields
!   12/2008 - Added SWIR and SST test space
!   1/2009 - Added Single Scattering Correction to Chn2 Reflectance
!
! Public Routines Within This Module
!   Baseline_Cloud_Mask_Main - main routine for cloud Mask
!
! Notes
!   - This cloud Mask will operate on more channels than is available on
!     GOES I-M, GOES NOP and SEVIRI.  Different channel combinations
!     can be used through channel selection IN the algorithm
!     definition module. However, do not specify a channel if needed
!     if it is not on the chosen platform.
!
!----------------------------------------------------------------------

!
!--- module use statements
!
use PIXEL_COMMON
use CONSTANTS
 use NWP_COMMON
 use RTM_COMMON
 use PLANCK
 use Message_handler
use calibration_constants, only: &
        sun_earth_distance  &
        , solar_ch20_nu
        
IMPLICIT NONE

PRIVATE

PUBLIC:: Baseline_Cloud_Mask_Main

PRIVATE:: Compute_Probably_Clear_Restoral
PRIVATE:: Compute_Probably_Cloudy
PRIVATE:: Clear_Chn2_Reflectance_Field
PRIVATE:: RUT_Routine
PRIVATE:: TUT_Routine
PRIVATE:: RTCT_Routine
PRIVATE:: ETROP_Routine
PRIVATE:: PFMFT_Routine
PRIVATE:: CIRH2O_Routine
PRIVATE:: RFMFT_Routine
PRIVATE:: TEMPIR_Routine
PRIVATE:: RGCT_Routine
PRIVATE:: RVCT_Routine
PRIVATE:: NIRREF_Chn5_Routine
PRIVATE:: NIRREF_Chn7_Routine
PRIVATE:: CIRREF_Routine
PRIVATE:: EMISS4_Routine
PRIVATE:: ULST_Routine
PRIVATE:: Set_Cmask_Thresholds

PRIVATE:: Compute_NWC
PRIVATE:: Compute_Clear_Sky_Scatter
PRIVATE:: Term_Refl_Norm

!
!--- include fixed thresholds via this include statement
!
INCLUDE 'baseline_cloud_mask_thresholds.inc'

CONTAINS

!====================================================================
! Baseline Cloud Mask from GOES-R AWG Cloud Application Team
!
! Principal Author: Andrew Heidinger
!
! Subroutine Name: Baseline_Cloud_Mask_Main
!
! Function:
!   Serve as the baseline cloud Mask for all ABI processing
!
! Description:
!   This routine derives a four-level cloud Mask for each pixel.
!   The 4 levels are clear, probably clear, probably cloudy and cloudy
!   The routine applies a series of cloud detection tests.
!   The Output includes all test decisions
!
! CALLing Sequence:
!   CALL Baseline_Cloud_Mask_Main(internal_algorithm_index)
!
! Inputs: All input passed through geocat structures (satellite, nwp, 
!         rtm and temporal)
!
! Outputs:
!      cloud_mask = out2(ialgo)%cldmask = 4 level cloud Mask
!      cloud_mask_packed = out2(ialgo)%cldmask_packed 
!                        =  test results packed into bytes
!      test_results = out2(ialgo)%qflg1 = actual test results (unpacked)
!      cloud_mask_qf = out2(ialgo)%qflg2 = quality flag for whole Mask
!      out2(ialgo)%r4_generic1 = generic REAL4 array used for diagnostic Output
!      out2(ialgo)%r4_generic2 = generic REAL4 array used for diagnostic Output
!      out2(ialgo)%r4_generic3 = generic REAL4 array used for diagnostic Output
!
! Dependencies:
!        GEOCAT satellite, rtm, nwp, and temporal structures must be
!        populated for this segment
!
! Restrictions:  TBD
!
!  History:
!   2/2007 - Andrew Heidinger - Created based on CLAVR-x
!   7/2007 - Added Temporal IR Test
!   1/2008 - Fixed bugs seen IN REAL-time processing
!   2/2008 - v3 delivered to AIT
!
! Local Variables:
! Elem_Idx - element number (the along scan pixel index)
! Line_Idx - line number (the across scan pixel index)
! Refl_Chn2 = 0.63 micron reflectance
! Rad_Chn7  = 3.9 micron radiance
! Refl_Chn4 =  1.38 micron reflectance
! Refl_Chn5 =  1.60 micron reflectance
! BT_Chn14  = 11.0 micron bt
! BT_Chn14_Clr = 11 micron clear bt
! Rad_Chn14  = 11.0 micron radiance
! Emiss_Chn7 = 3.9 micron emissivity
! Sfc_Emiss_Chn7_RTM = 3.9 micron surface emissivity
! X_NWP_Idx = nwp longitude cell
! Y_NWP_Idx = nwp latitude cell
! Rad_Chn7_Clr = clear 4 micron radiance
! Rad_Chn14_Clr = clear 11 micron radiance
! Cos_Sat_Zen  = cosine of satellite viewing zenith angle
! Cos_Sol_Zen  = cosine of solar viewing zenith angle
! Sfc_Idx_NWP = nwp level associated with surface
! Tropo_Idx_NWP = nwp level associated with tropopause
! View_Zen_Idx = viewing zenith angle bin
! rad14_bb_trop = BB 11 micron radiance at tropopause
! trans_atm_7 = BB 12 micron radiance at tropopause
! Sfc_Temp_Uni_NWP = 3x3 standard deviation of NWP surface temperature
! Refl_Chn2_Clear_Stddev_3x3 = 3x3 standard deviation of NWP surface temperature
! solar_7 = solar energy IN channel 7
! sed = sun earth distance factor
! glintzen = glint zenith angle
! landmask  = land Mask
! Coast_Type =  = coast Mask
! solarzen  = solar zenith angle
! desertmask = desert Mask
! snowmask  = snow Mask
!
!
! Elements of the Test Vector IN order as written OUT
!
!--- ancillary data flags                        Output location
! 1 - cloud Mask attempted                        byte 1 bit 1
! 2 - day                                         byte 1 bit 2
! 3 - terminator                                  byte 1 bit 3
! 4 - land                                        byte 1 bit 4
! 5 - coast                                       byte 1 bit 5
! 6 - glint                                       byte 1 bit 6
! 7 - desert                                      byte 1 bit 7
! 8 - snow                                        byte 1 bit 8
! 9 - cold surface                                byte 2 bit 1
!--- clear-sky spatial uniformity tests
! 10 - RUT                                        byte 2 bit 2
! 11 - TUT                                        byte 2 bit 3
!--- cloud tests using thermal emission
! 12 - RTCT                                       byte 2 bit 4
! 13 - ETROP                                      byte 2 bit 5
! 14 - PFMFT                                      byte 2 bit 6
! 15 - NFMFT                                      byte 2 bit 7
! 16 - RFMFT                                      byte 2 bit 8
! 17 - CIRH2O                                     byte 3 bit 1
! 18 - TEMPIR                                     byte 3 bit 2
! 19 - TERM_THERM_STAB                            byte 3 bit 3
!--- cloud tests using solar reflectance 
! 20 - RGCT                                       byte 3 bit 4
! 21 - RVCT                                       byte 3 bit 5
!--- cloud tests using swir solar reflectance 
! 22 - NIRREF                                     byte 3 bit 6
! 23 - CIRREF                                     byte 3 bit 7
!--- cloud tests using swir thermal emission
! 24 - EMISS4                                     byte 4 bit 8
! 25 - ULST                                       byte 4 bit 1
!--- restoral tests
! 26 - probably clear restoral                    byte 4 bit 2
! 27 - cloud edge (prob cloudy)                   byte 4 bit 3
!--- non-cloud tests 
! 28 - blank                                      byte 4 bit 4
! 29 - blank                                      byte 4 bit 5
! 30 - blank                                      byte 4 bit 6
! 31 - blank                                      byte 4 bit 7
! 32 - blank                                      byte 4 bit 8
!
!
! Note the following
!  - you must ensure that geocat array hold test_results IN big enough
!    this is done IN algorithm_mod.f90
!====================================================================
SUBROUTINE Baseline_Cloud_Mask_Main(Algo_Num)

   INTEGER(KIND=INT4), INTENT(IN) :: Algo_Num
   CHARACTER(len=*), PARAMETER:: Routine_Name =  &
                                 "baseline_cloud_mask: test_Baseline_Cloud_Mask_main"
   INTEGER(KIND=INT4) :: Elem_Idx
   INTEGER(KIND=INT4) :: Line_Idx
   INTEGER(KIND=INT4) :: Num_Elem
   INTEGER(KIND=INT4) :: Number_of_Lines_in_this_Segment
   INTEGER(KIND=INT4) :: Max_Num_Lines_per_Seg
   CHARACTER(LEN=10) :: Sat_Name
   CHARACTER(LEN=40) :: Algo_Name
   CHARACTER(LEN=100) :: Err_Message
   INTEGER(KIND=INT4) :: Chn_Num
   INTEGER(KIND=INT4) :: Elem_LRC_Idx
   INTEGER(KIND=INT4) :: Line_LRC_Idx
   INTEGER(KIND=INT4) :: Elem_NWC_Idx
   INTEGER(KIND=INT4) :: Line_NWC_Idx
   REAL(KIND=REAL4) :: Refl_Chn2
   REAL(KIND=REAL4) :: Refl_Chn4
   REAL(KIND=REAL4) :: Refl_Chn5
   REAL(KIND=REAL4) :: BT_Chn7
   REAL(KIND=REAL4) :: Rad_Chn7
   REAL(KIND=REAL4) :: Refl_Chn7
   REAL(KIND=REAL4) :: Emiss_Chn7
   REAL(KIND=REAL4) :: Emiss_Chn7_Clr
   REAL(KIND=REAL4) :: Rad_Chn7_Clr
   REAL(KIND=REAL4) :: BT_Chn9
   REAL(KIND=REAL4) :: BT_Chn10
   REAL(KIND=REAL4) :: BT_Chn14
   REAL(KIND=REAL4) :: BT_Chn14_Clr
   REAL(KIND=REAL4) :: BT_Chn15
   REAL(KIND=REAL4) :: BT_Chn15_Clr
   INTEGER(KIND=INT4) :: Alloc_Status
   INTEGER(KIND=INT4) :: Alloc_Status_Total
   INTEGER(KIND=INT4) :: Have_Prev_BT_Chn14_15min
   INTEGER(KIND=INT4) :: Have_Prev_BT_Chn14_Clr_15min
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: BT_WV_BT_Window_Corr
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: BT_WaterVapor_Mean_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: BT_WaterVapor_Max_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: BT_WaterVapor_Min_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: BT_WaterVapor_Stddev_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: BT_Chn14_Mean_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: BT_Chn14_Max_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: BT_Chn14_Min_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: BT_Chn14_Stddev_3x3
   INTEGER(KIND=INT4), DIMENSION(:,:), ALLOCATABLE:: X_NWC_Idx
   INTEGER(KIND=INT4), DIMENSION(:,:), ALLOCATABLE:: Y_NWC_Idx

   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Refl_Chn2_Mean_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Refl_Chn2_Max_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Refl_Chn2_Min_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Refl_Chn2_Stddev_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Sfc_Hgt_Mean_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Sfc_Hgt_Max_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Sfc_Hgt_Min_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Sfc_Hgt_Stddev_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Refl_Chn2_Clear_Mean_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Refl_Chn2_Clear_Max_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Refl_Chn2_Clear_Min_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Refl_Chn2_Clear_Stddev_3X3
   INTEGER(KIND=INT4):: Num_Pix
   INTEGER(KIND=INT4):: Snow_Mask
   INTEGER(KIND=INT4):: Desert_Mask
   INTEGER(KIND=INT4):: Coast_Type
   INTEGER(KIND=INT4):: Land_Type
   INTEGER(KIND=INT4):: Is_Day
   INTEGER(KIND=INT4):: Is_Terminator
   INTEGER(KIND=INT4):: Is_Land
   INTEGER(KIND=INT4):: Is_Coast
   INTEGER(KIND=INT4):: Is_Glint
   INTEGER(KIND=INT4):: Is_Desert
   INTEGER(KIND=INT4):: Is_Snow
   INTEGER(KIND=INT4):: Is_Valid_Pixel
   INTEGER(KIND=INT4), DIMENSION(16):: Is_Chn
   INTEGER(KIND=INT4):: Num_Tests
   INTEGER(KIND=INT4):: Error_Level
   REAL(KIND=REAL4):: Total_Precipitable_Water_NWP
   REAL(KIND=REAL4):: Total_Ozone_Path_NWP
   REAL(KIND=REAL4):: Sfc_Temp_Uni_NWP
   REAL(KIND=REAL4):: Sfc_Temp
   REAL(KIND=REAL4):: Sol_Zen
   REAL(KIND=REAL4):: Sfc_Hgt
   REAL(KIND=REAL4):: Scat_Zen
   REAL(KIND=REAL4):: Emiss_Tropo_Chn14_LRC
   REAL(KIND=REAL4):: BTD_Chn14_Chn15_NWC
   REAL(KIND=REAL4):: Emiss_Chn7_NWC
   REAL(KIND=REAL4):: BT_Chn14_NWC
   REAL(KIND=REAL4):: Planck_Emission_Chn7_Clr
   REAL(KIND=REAL4):: Chn7_Sol_Energy
   REAL(KIND=REAL4):: Sun_Earth_Dist
   REAL(KIND=REAL4):: Glint_Zen
   REAL(KIND=REAL4):: Sfc_Emiss_Chn7_RTM
   REAL(KIND=REAL4):: Atm_Trans_Chn7_RTM
   REAL(KIND=REAL4):: Atm_Solar_Trans_Chn7
   REAL(KIND=REAL4):: Solar_Rad_Chn7_Clr
   REAL(KIND=REAL4):: Solar_BT_Chn7_Clr
   REAL(KIND=REAL4):: Cos_Sat_Zen
   REAL(KIND=REAL4):: Cos_Scat_Zen
   REAL(KIND=REAL4):: Cos_Sol_Zen
   INTEGER(KIND=INT1):: Tropo_Idx_NWP
   INTEGER(KIND=INT1):: View_Zen_Idx
   INTEGER(KIND=INT1):: Sfc_Idx_NWP
   INTEGER(KIND=INT1):: Is_Coast_NWP
   INTEGER(KIND=INT4):: X_NWP_Idx
   INTEGER(KIND=INT4):: Y_NWP_Idx
   REAL(KIND=REAL4):: BT_Chn14_15min
   REAL(KIND=REAL4):: BT_Chn14_15min_Clr
   INTEGER:: Is_Cold_Surface !a flag that identifies cold surfaces
   INTEGER(KIND=INT4):: Uni_Land_Mask_Flag_Yes
   INTEGER(KIND=INT4):: Uni_Land_Mask_Flag_No
   REAL(KIND=REAL4):: Aerosol_Optical_Depth_Chn2
   REAL(KIND=REAL4):: NDSI
   
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Refl_Chn2_Clear

   REAL(KIND=REAL4):: RGCT_Threshold

   INTEGER:: Array_Right, Array_Left, Array_Top
   INTEGER:: Array_Bottom,Array_Width, Array_Hgt

   !--- quality flag 
   INTEGER(KIND=INT4):: IR_Test_Sum
   INTEGER(KIND=INT4):: VIS_Test_Sum
   INTEGER(KIND=INT4):: SWIR_SOLAR_Test_Sum
   INTEGER(KIND=INT4):: SWIR_THERMAL_Test_Sum
   INTEGER(KIND=INT4):: IR_Test_Mask
   INTEGER(KIND=INT4):: VIS_Test_Mask
   INTEGER(KIND=INT4):: SWIR_Test_Mask

   REAL(KIND=REAL4):: Air_Mass_Factor
   REAL(KIND=REAL4):: Transmission_Sing_Scat
   REAL(KIND=REAL4):: Refl_Sing_Scat

   !
   !--- cloud Mask arrays
   !
   INTEGER(KIND=INT4), DIMENSION(Total_Num_Tests):: Test_Bit_Depth

   !
   !--- local pointers
   !
   INTEGER(KIND=INT1), DIMENSION(:,:), POINTER :: Cloud_Mask
   INTEGER(KIND=INT1), DIMENSION(:,:), POINTER :: Cloud_Mask_Binary
   INTEGER(KIND=INT1), DIMENSION(:,:), POINTER :: Cloud_Mask_IR
   INTEGER(KIND=INT1), DIMENSION(:,:), POINTER :: Cloud_Mask_SST
   INTEGER(KIND=INT1), DIMENSION(:,:,:), POINTER :: Cloud_Mask_Packed
   INTEGER(KIND=INT1), DIMENSION(:,:,:), POINTER :: Test_Results
   INTEGER(KIND=INT1), DIMENSION(:,:), POINTER :: Cloud_Mask_QF
   INTEGER(KIND=INT1), DIMENSION(:,:), POINTER :: Cloud_Mask_Tmpy
   INTEGER(KIND=INT4), DIMENSION(:,:), POINTER :: X_LRC_Idx
   INTEGER(KIND=INT4), DIMENSION(:,:), POINTER :: Y_LRC_Idx
   INTEGER(KIND=INT1), DIMENSION(:,:), POINTER :: LRC_Mask
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER :: Emiss_Tropo_Chn14


   ! ------------- Variables needed for Term stability test
   REAL(KIND=REAL4) :: BT_Chn11
   REAL(KIND=REAL4) :: BT_Chn11_1Hr
   REAL(KIND=REAL4) :: BT_Chn14_1Hr
   REAL(KIND=REAL4) :: BT_Chn15_1Hr
   INTEGER(KIND=INT4):: Cmask_1Hr
   INTEGER(KIND=INT4) :: Have_Prev_BT_Chn11_1Hr
   INTEGER(KIND=INT4) :: Have_Prev_BT_Chn14_1Hr
   INTEGER(KIND=INT4) :: Have_Prev_BT_Chn15_1Hr
   INTEGER(KIND=INT4) :: Have_Prev_Cmask_1Hr

   !----------------------------------------------------------------------
   ! Executable Code
   !----------------------------------------------------------------------

   !--- store size of this segment into local variables
   Num_Elem = Image%Number_Of_Elements
   Number_of_Lines_in_this_Segment = Image%Number_Of_Lines_Read_This_Segment
   Max_Num_Lines_per_Seg = Image%Number_Of_Lines_Per_Segment


   !set the initial cloud mask algorithm name here (for usage in the temporal term test)
   Algo_Name = 'baseline_cmask_'

   !--- store name of sensor
   ! WCS - FIXME
   Sat_Name = trim(scinfo(sc_ind)%name)
   
   !-----------------------------------------------------------------
   ! Set sensor dependent thresholds
   !  This is due to a SRF issue on the 3.9 micron band
   !  Also compose sensor dependent algorithm name for temporal term test
   !-----------------------------------------------------------------
   
   CALL Set_Cmask_Thresholds(Sat_Name, Algo_Name)
   

   !-----------------------------------------------------------------
   ! Read IN background Chn2 Reflectance Field 
   !-----------------------------------------------------------------
   !CALL Clear_Chn2_Reflectance_Field(Num_Elem, &
   !                                  Max_Num_Lines_per_Seg, &
   !                                  Refl_Chn2_Clear)

   ! Data is already normalized in CLAVRx
   Refl_Chn2_Clear => ch(1)%Ref_Toa_Clear
   !-----------------------------------------------------------------
   ! Initialize masks
   !-----------------------------------------------------------------
   Cloud_Mask => Cld_Mask(i,j)
   Cloud_Mask_Binary => null()
   Cloud_Mask_IR => null()
   Cloud_Mask_SST => null()
   Cloud_Mask_Packed => Cld_Test_Vector_Packed(:,i,j)
   Test_Results =>
   Cloud_Mask_QF => 
   Cloud_Mask_Tmpy => 
   X_LRC_Idx => 
   Y_LRC_Idx => 
   LRC_Mask => 
   Emiss_Tropo_Chn14 => 
   
   !
   !--- set bit depths for packed Output (all here are 1 bit)
   !
   Test_Bit_Depth = 1

   !---------------------------------------------
   !--- initialize Output
   !---------------------------------------------

   !
   !---don't initialize to clear
   !
   Cloud_Mask = sym%PROB_CLEAR   
   Cloud_Mask_Binary = CLD_MASK_BIN_CLEAR  
   Cloud_Mask_IR = sym%PROB_CLEAR   
   Cloud_Mask_SST = sym%PROB_CLEAR   

   !
   !---This must be to zero. So the bit associated
   !
   Cloud_Mask_Packed = 0         

   !
   !with ivalid is initialized to zero
   !
   Test_Results = 0

   !Initialize quality flag to invalid
   Cloud_Mask_QF = INVALID_CMASK_RETREVAL

   !=======================================================================
   ! Compute spatial uniformity metrics
   !=======================================================================
   IF ((Sensor%Chan_On_Flag_Default(27) > 0) .OR. (Sensor%Chan_On_Flag_Default(28)  > 0)) THEN
      IF (Sensor%Chan_On_Flag_Default(28)  > 0) THEN     !Choose Chn10 if available
        CALL Compute_Spatial_Uniformity(1, 1, Space_Mask, ch(28)%Bt_Toa, BT_WaterVapor_Mean_3x3, &
                                        BT_WaterVapor_Max_3x3, BT_WaterVapor_Min_3x3, BT_WaterVapor_Stddev_3x3)
      ELSE                                        !If Chn10 not available, use Chn9
        CALL Compute_Spatial_Uniformity(1, 1, Space_Mask, ch(27)%Bt_Toa, BT_WaterVapor_Mean_3x3, &
                                        BT_WaterVapor_Max_3x3, BT_WaterVapor_Min_3x3, BT_WaterVapor_Stddev_3x3)
      ENDIF
   ENDIF

   CALL Compute_Spatial_Uniformity(1, 1, Space_Mask, ch(31)%Bt_Toa, BT_Chn14_Mean_3x3, &
                                   BT_Chn14_Max_3x3, BT_Chn14_Min_3x3, BT_Chn14_Stddev_3x3)
   CALL Compute_Spatial_Uniformity(1, 1, Space_Mask, ch(1)%Ref_Toa, Refl_Chn2_Mean_3X3, &
                                   Refl_Chn2_Max_3x3, Refl_Chn2_Min_3X3, Refl_Chn2_Stddev_3X3)
   CALL Compute_Spatial_Uniformity(1, 1, Space_Mask, Sfc%zsfc, Sfc_Hgt_Mean_3X3, &
                                   Sfc_Hgt_Max_3X3, Sfc_Hgt_Min_3X3, Sfc_Hgt_Stddev_3X3)
   CALL Compute_Spatial_Uniformity(1, 1, Space_Mask, Refl_Chn2_Clear, Refl_Chn2_Clear_Mean_3X3, &
                                   Refl_Chn2_Clear_Max_3x3, Refl_Chn2_Clear_Min_3X3, Refl_Chn2_Clear_Stddev_3X3)

   ! Temporal commented out due to CLAVR-x not having that capability

   !=======================================================================
   ! Load temporal information   (status = sym%FAILURE or sym%SUCCESS)
   !=======================================================================

   !--- previous 11 micron temp
!   Have_Prev_BT_Chn14_15min = Load_Temporal_Data(minus_i_15min,"",           &
!                                            temporal(minus_i_15min)%bt14)
!   Have_Prev_BT_Chn14_Clr_15min = Load_Temporal_Data(minus_i_15min,"",       &
!                                        temporal(minus_i_15min)%bt_clr14)

   !=======================================================================
   ! Load TERM_THERM_STAB temporal information   
   ! (status = sym%FAILURE or sym%SUCCESS)
   !=======================================================================
!   Have_Prev_BT_Chn11_1Hr = Load_Temporal_Data(minus_i_01hrs,"",           &
!                                            temporal(minus_i_01hrs)%bt11)
!   Have_Prev_BT_Chn14_1Hr = Load_Temporal_Data(minus_i_01hrs,"",           &
!                                            temporal(minus_i_01hrs)%bt14)
!   Have_Prev_BT_Chn15_1Hr = Load_Temporal_Data(minus_i_01hrs,"",           &
!                                            temporal(minus_i_01hrs)%bt15)
!   Have_Prev_Cmask_1Hr = Load_Temporal_Data(minus_i_01hrs,TRIM(Algo_Name), &
!                                            temporal(minus_i_01hrs)%cldmask)

   !=======================================================================
   ! Compute 11 micron emissivity at Tropopause
   !=======================================================================
   CALL Compute_Emiss_Tropo_Chn14(Emiss_Tropo_Chn14, &
                                  Number_of_Lines_in_this_Segment)

   !======================================================================
   ! compute local radiative center
   !======================================================================

   !--- determine which pixels to skip
   LRC_Mask = sym%NO
   
   WHERE(Emiss_Tropo_Chn14 /= Missing_Value_Real4)
        LRC_Mask = sym%YES
   ENDWHERE

   !--- initialize indices
   X_LRC_Idx = Missing_Value_INT4
   Y_LRC_Idx = Missing_Value_INT4

   !--- call routines to compute emiss lrc indices
   CALL Gradient2d(Emiss_Tropo_Chn14, &
                  Image%Number_Of_Elements, &
                  Image%Number_Of_Lines_Read_This_Segment, &
                  LRC_Mask, &
                  EMISS_TROPO_CHN14_GRADIENT_MIN, & 
                  EMISS_TROPO_CHN14_GRADIENT_MAX, & 
                  EMISS_TROPO_CHN14_GRADIENT_THRESH, & 
                  X_LRC_Idx, &
                  Y_LRC_Idx)

   !-------------------------------------------------------------------------------
   ! Ensure that pixels with values above the threshold are
   ! treated as LRC's.  This not the default behavior
   !-------------------------------------------------------------------------------
   DO Line_Idx=1, Number_of_Lines_in_this_Segment
      DO Elem_Idx=1, Num_Elem
         IF (Emiss_Tropo_Chn14(Elem_Idx,Line_Idx) >= EMISS_TROPO_CHN14_GRADIENT_THRESH ) THEN
          X_LRC_Idx(Elem_Idx,Line_Idx) = Elem_Idx
          Y_LRC_Idx(Elem_Idx,Line_Idx) = Line_Idx
         ENDIF
      ENDDO
   ENDDO

   !----------------------------------------------------------------------
   ! set flag to enforce uniformity IN land Mask for these computations
   !----------------------------------------------------------------------
   Uni_Land_Mask_Flag_Yes  = sym%YES
   Uni_Land_Mask_Flag_No  = sym%NO
   

   !----------------------------------------------------------------------
   ! compute NWC indices
   !----------------------------------------------------------------------
   Alloc_Status_Total = 0 
    
   ALLOCATE(X_NWC_Idx(Num_Elem,Max_Num_Lines_per_Seg),stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status 
    
   ALLOCATE(Y_NWC_Idx(Num_Elem,Max_Num_Lines_per_Seg),stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status 
   
   !======================================================================
   ! Check allocation of NWC arrays
   !======================================================================

   IF (Alloc_Status_Total /= 0) THEN
       WRITE (Err_Message, *) &
             'Error allocating NWC arrays'
             
             
        Error_Level = 2 ! AIT FATAL ERROR CODE
            
        CALL Display_Message("Baseline Cloud Mask", &
                 TRIM(Err_Message), &
                 Sym%FAILURE)
                 
        !  AIT Error Messaging
        !  CALL Error_Messaging (Routine_Name, Error_Message, Error_Level)           
        RETURN

   ENDIF

   !------------------------------------------------------------------------
   ! Call Routine to Compute Neighboring Warm Center
   !------------------------------------------------------------------------
   CALL Compute_NWC( &
                     ch(31)%Bt_Toa, &
                     NWC_PIXEL_RADIUS, &  !nbox
                     Uni_Land_Mask_Flag_Yes, &  !uni_land_mask_flag
                     Bad_Pixel_Mask(:,:), &
                     Sfc%Land, &
                     1, &
                     Num_Elem, &
                     1, &
                     Max_Num_Lines_per_Seg, &
                     X_NWC_Idx, & 
                     Y_NWC_Idx)

   !----------------------------------------------------------------------
   ! compute Correlation of Chn10 and Chn14
   !----------------------------------------------------------------------
   IF ((Sensor%Chan_On_Flag_Default(27) > 0) .OR. (Sensor%Chan_On_Flag_Default(28)  > 0) .AND. &
        (out2(Algo_Num)%ch_flg(14) > 0)) THEN
   
        Alloc_Status_Total = 0 
        

        ALLOCATE(BT_WV_BT_Window_Corr(Num_Elem,Max_Num_Lines_per_Seg), &
                 stat=Alloc_Status)
                 
        Alloc_Status_Total = Alloc_Status_Total + Alloc_Status 

        !======================================================================
        ! Check allocation of arrays
        !======================================================================
        
        IF (Alloc_Status_Total /= 0) THEN
           
            WRITE (Err_Message, *) &
             'Error allocating BT_WV_BT_Window_Corr arrays'
             
            Error_Level = 2 ! AIT FATAL ERROR CODE
            
            CALL Display_Message("Baseline Cloud Mask", &
                 TRIM(Err_Message), &
                 Sym%FAILURE)
                 
            !  AIT Error Messaging
            !  CALL Error_Messaging (Routine_Name, Error_Message, Error_Level)
            
            RETURN

        ENDIF
 
        Line_Loop_Corr: DO Line_Idx=1, Number_of_Lines_in_this_Segment
                Element_Loop_Corr: DO Elem_Idx=1, Num_Elem
                
                    Array_Right = max(1,min(Elem_Idx - 2,Num_Elem))
                    Array_Left = max(1,min(Elem_Idx + 2,Num_Elem))
                    Array_Top = max(1,min(Line_Idx - 2,Max_Num_Lines_per_Seg))
                    Array_Bottom = max(1,min(Line_Idx + 2,Max_Num_Lines_per_Seg))
                    Array_Width = Array_Left -Array_Right + 1
                    Array_Hgt = Array_Bottom -Array_Top + 1

                    IF (Space_Mask(Elem_Idx,Line_Idx) == sym%SPACE) THEN
                    
                        BT_WV_BT_Window_Corr(Elem_Idx,Line_Idx) = Missing_Value_Real4
                        
                        CYCLE
                        
                    ENDIF

                    
                    IF (Sensor%Chan_On_Flag_Default(28)  > 0) THEN
                         BT_WV_BT_Window_Corr(Elem_Idx,Line_Idx) = Pearson_Corr(&
                                                       ch(28)%Bt_Toa(Array_Right:Array_Left,Array_Top:Array_Bottom), &
                                                       ch(31)%Bt_Toa(Array_Right:Array_Left,Array_Top:Array_Bottom), &
                                                       Bad_Pixel_Mask(Array_Right:Array_Left,Array_Top:Array_Bottom), &
                                                       Bad_Pixel_Mask(Array_Right:Array_Left,Array_Top:Array_Bottom), &
                                                       Array_Width, Array_Hgt)
                    ELSE
                    
                         BT_WV_BT_Window_Corr(Elem_Idx,Line_Idx) = Pearson_Corr( &
                                                       ch(27)%Bt_Toa(Array_Right:Array_Left,Array_Top:Array_Bottom), &
                                                       ch(31)%Bt_Toa(Array_Right:Array_Left,Array_Top:Array_Bottom), &
                                                       Bad_Pixel_Mask(Array_Right:Array_Left,Array_Top:Array_Bottom), &
                                                       Bad_Pixel_Mask(Array_Right:Array_Left,Array_Top:Array_Bottom), &
                                                       Array_Width, Array_Hgt)
                    ENDIF
                    
                END DO Element_Loop_Corr
        END DO Line_Loop_Corr

   ENDIF
      
   !=======================================================================
   ! Loop over pixels apply cloud tests
   !=======================================================================

   Line_Loop_2: DO Line_Idx=1, Number_of_Lines_in_this_Segment

      Element_Loop_2: DO Elem_Idx=1, Num_Elem

         !
         !--- check for space pixel. QF already set to invalid
         !
         IF (Space_Mask(Elem_Idx,Line_Idx) == sym%SPACE) THEN
           
            CYCLE

         ENDIF

         !
         !--- check for pixel located beyond zone of operation. Set QF
         !
         IF (Geo%Solzen(Elem_Idx,Line_Idx) > SENSOR_ZEN_THRESH) THEN

            Cloud_Mask_QF(Elem_Idx,Line_Idx) = CMASK_OUTSIDE_SEN_ZEN_RANGE
            
            CYCLE

         ENDIF


         !
         !--- initialize channel flags
         !
         Is_Chn = sym%NO

         !
         !--- define aliases
         !  local name          global name                  
         !

         !
         !---nwp longitude cell
         !
         X_NWP_Idx =          I_Nwp(Elem_Idx,Line_Idx)         

         !
         !---nwp latitude cell
         !
         Y_NWP_Idx =          J_Nwp(Elem_Idx,Line_Idx)     

         !
         !--- LRC indices
         Elem_LRC_Idx = X_LRC_Idx(Elem_Idx,Line_Idx)
         Line_LRC_Idx = Y_LRC_Idx(Elem_Idx,Line_Idx)
         Emiss_Tropo_Chn14_LRC = Missing_Value_Real4
         IF (Elem_LRC_Idx > 0 .and. Line_LRC_Idx > 0) THEN
           Emiss_Tropo_Chn14_LRC = Emiss_Tropo_Chn14(Elem_LRC_Idx,Line_LRC_Idx)
         ENDIF

         !
         !--- NWC Indices
         Elem_NWC_Idx = X_NWC_Idx(Elem_Idx,Line_Idx)
         Line_NWC_Idx = Y_NWC_Idx(Elem_Idx,Line_Idx)

         !--- store BTD_Chn14_Chn15 at NWC
         BTD_Chn14_Chn15_NWC = Missing_Value_Real4
         IF (out2(Algo_Num)%ch_flg(14) > 0 .and. Sensor%Chan_On_Flag_Default(32) > 0) THEN
           IF (Elem_NWC_Idx > 0 .and. Line_NWC_Idx > 0) THEN
                BTD_Chn14_Chn15_NWC = ch(31)%Bt_Toa(Elem_NWC_Idx,Line_NWC_Idx) - &
                                      ch(32)%Bt_Toa(Elem_NWC_Idx,Line_NWC_Idx)
           ENDIF
         ENDIF

         !
         !---cosine of satellite viewing zenith angle
         !
         Cos_Sat_Zen  =       Geo%Coszen(Elem_Idx,Line_Idx)    

         !
         !---cosine of solar viewing zenith angle
         !
         Cos_Sol_Zen =     Geo%CosSolzen(Elem_Idx,Line_Idx)    

         !
         !---cosine of scattering zenith angle
         !
         Cos_Scat_Zen =     cos(Geo%Scatangle(Elem_Idx,Line_Idx)*dtor)

         !
         !---cosine of scattering zenith angle
         !
         Air_Mass_Factor =  1.0/Cos_Sol_Zen + 1.0/Cos_Sat_Zen

         !
         !---nwp level associated with surface
         !
         Sfc_Idx_NWP =          Rtm(X_NWP_Idx,Y_NWP_Idx)%Sfc_Level

         !
         !---nwp level associated with tropopause
         !
         Tropo_Idx_NWP =        Rtm(X_NWP_Idx,Y_NWP_Idx)%Tropo_Level 

         !
         !---viewing zenith angle bin
         !
         View_Zen_Idx =          Zen_Idx_Rtm(Elem_Idx,Line_Idx)        

         !
         !--- Surface Temperature and its 3x3 uniformity
         !
         Sfc_Temp_Uni_NWP =  Tmpsfc_uni_Nwp(X_NWP_Idx,Y_NWP_Idx)
     
         Sfc_Temp         =  Tmpsfc_Nwp(X_NWP_Idx,Y_NWP_Idx)

         !
         !--- Total Precipitable Water (g/m^2 or cm)
         !
         Total_Precipitable_Water_NWP =  Tpw_Nwp(X_NWP_Idx,Y_NWP_Idx)

         !
         !--- Total Ozone (Dobson Unit)
         !
         Total_Ozone_Path_NWP =  Ozone_Nwp(X_NWP_Idx,Y_NWP_Idx)

         !
         !---sun earth distance factor
         !
         Sun_Earth_Dist =           sat%sun_earth_distance         

         !
         !---glint zenith angle
         !
         Glint_Zen =      Sfc%Glint_Mask(Elem_Idx,Line_Idx)       

         !
         !---land type
         !
         Land_Type =      Sfc%Land(Elem_Idx,Line_Idx)    

         !
         !---coast type
         !
         Coast_Type =     Sfc%Coast_Mask(Elem_Idx,Line_Idx)    

         !
         !---Surface Height 
         !
         Sfc_Hgt =      Sfc%zsfc(Elem_Idx,Line_Idx)         

         !
         !---solar zenith angle
         !
         Sol_Zen =      Geo%Solzen(Elem_Idx,Line_Idx)         

         !
         !---solar scattering zenith angle
         !
         Scat_Zen =      Geo%Scatangle(Elem_Idx,Line_Idx)         

         !
         !---desert Mask
         !
         Desert_Mask =    Sfc%Desert_Mask(Elem_Idx,Line_Idx)    

         !
         !---snow Mask
         !
         Snow_Mask =      Sfc%Snow(Elem_Idx,Line_Idx)     

         !--------------------------------------------------------
         !--- local aliases for observations that are available
         !--------------------------------------------------------

         !--- Channel 2 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(1) > 0) THEN

            Is_Chn(2) = sym%YES

            !
            !---0.63 micron reflectance
            !
            Refl_Chn2 = ch(1)%Ref_Toa(Elem_Idx,Line_Idx)         

            !--- renormalize for improved terminator performance
            IF (Sol_Zen > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) THEN
              Refl_Chn2 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn2)
            ENDIF

         ENDIF

         !--- Channel 5 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(26) > 0) THEN

            Is_Chn(4) = sym%YES

            !
            !---1.38 micron reflectance
            !
            Refl_Chn4 =          ch(26)%Ref_Toa(Elem_Idx,Line_Idx)       

            !--- renormalize for improved terminator performance
            IF (Sol_Zen > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) THEN
              Refl_Chn4 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn4)
            ENDIF

         ENDIF

         !--- Channel 5 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(6) > 0) THEN

            Is_Chn(5) = sym%YES

            !
            !---1.60 micron reflectance
            !
            Refl_Chn5 =          ch(6)%Ref_Toa(Elem_Idx,Line_Idx)          

            !--- renormalize for improved terminator performance
            IF (Sol_Zen > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) THEN
              Refl_Chn5 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn5)
            ENDIF

         ENDIF

         !--- Channel 7 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(20) > 0) THEN

            Is_Chn(7) = sym%YES

            !
            !---3.9 micron brightness temperature
            !
            BT_Chn7  =         ch(20)%Bt_Toa(Elem_Idx,Line_Idx)        

            !
            !---3.9 micron radiance
            !
            Rad_Chn7  =         ch(20)%Rad_Toa(Elem_Idx,Line_Idx)        

            !
            !---3.9 um pseudo reflectance
            !
            Refl_Chn7 =        ch(20)%Ref_Toa(Elem_Idx,Line_Idx)      

            !--- renormalize for improved terminator performance
            IF (Sol_Zen > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) THEN
              Refl_Chn7 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn7)
            ENDIF

            !
            !---3.9 um emissivity
            !

            Emiss_Chn7 =        Ems_Ch20(Elem_Idx,Line_Idx)
                  
            Emiss_Chn7_NWC = Missing_Value_Real4
            IF (Elem_NWC_Idx > 0 .and. Line_NWC_Idx > 0) THEN
                Emiss_Chn7_NWC = Ems_Ch20(Elem_NWC_Idx,Line_NWC_Idx)
            ENDIF
           
            !
            !---3.9 um surface emissivity
            !
            Sfc_Emiss_Chn7_RTM =   ch(20)%Sfc_Emiss(Elem_Idx,Line_Idx)    

            !
            !---clear 4 micron radiance
            !
            Rad_Chn7_Clr =      ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx)     

            !
            !--- Atmospheric Transmission IN Channel 7
            !
            Atm_Trans_Chn7_RTM =   rtm(X_NWP_Idx,Y_NWP_Idx)%d(View_Zen_Idx)%ch(20)%trans_atm_clr7(Sfc_Idx_NWP) 

            !
            !---solar energy IN channel 7
            !
            Chn7_Sol_Energy =       scinfo(sc_ind)%S(7,sat%idet)    

         ENDIF


         !--- Channel 10 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(27) > 0) THEN

            Is_Chn(9) = sym%YES

            !
            !---6.7 micron bt
            !
            BT_Chn9  =         ch(27)%Bt_Toa(Elem_Idx,Line_Idx)

         ENDIF

         !--- Channel 10 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(28)  > 0) THEN

            Is_Chn(10) = sym%YES

            !
            !---7.3 micron bt
            !
            BT_Chn10  =         ch(28)%Bt_Toa(Elem_Idx,Line_Idx)

         ENDIF

         !--- Channel 11 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(29) > 0) THEN

            Is_Chn(11) = sym%YES

            !
            !---8.5 micron bt
            !
            BT_Chn11  =         ch(29)%Bt_Toa(Elem_Idx,Line_Idx)

         ENDIF


         !--- Channel 14 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(31) > 0) THEN

            Is_Chn(14) = sym%YES

            !
            !---11.0 micron bt
            !
            BT_Chn14  =         ch(31)%Bt_Toa(Elem_Idx,Line_Idx)         

            !
            !---11 um clear bt
            !
            BT_Chn14_Clr =      ch(31)%Bt_Toa_Clear(Elem_Idx,Line_Idx)     

            BT_Chn14_NWC = Missing_Value_Real4

            IF (Elem_NWC_Idx > 0 .and. Line_NWC_Idx > 0) THEN
                BT_Chn14_NWC = ch(31)%Bt_Toa(Elem_NWC_Idx,Line_NWC_Idx)
            ENDIF

         ENDIF

         !--- Channel 15 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(32) > 0) THEN

            Is_Chn(15) = sym%YES

            !
            !---12 micron bt
            !
            BT_Chn15  =  ch(32)%Bt_Toa(Elem_Idx,Line_Idx)         

            !
            !---12 um clear bt
            !
            BT_Chn15_Clr =  ch(32)%Bt_Toa_Clear(Elem_Idx,Line_Idx)     

         ENDIF

         !--------------------------------------------------------------------
         !--- compute a transmission term for solar +viewing path for ch7
         !--------------------------------------------------------------------
         IF (Is_Chn(7) == sym%YES) THEN

            Atm_Solar_Trans_Chn7 = 0.0

            IF ((Cos_Sol_Zen > 0.0).AND.(Cos_Sat_Zen > 0.0)) THEN

               Atm_Solar_Trans_Chn7 = Atm_Trans_Chn7_RTM ** (1.0 + Cos_Sat_Zen / Cos_Sol_Zen)

            ENDIF

         ENDIF

         
         !--------------------------------------------------------------------
         !---- alias temporal parameters
         !--------------------------------------------------------------------
         IF (Is_Chn(14) == sym%YES) THEN

            IF (Have_Prev_BT_Chn14_15min == sym%SUCCESS) THEN

               BT_Chn14_15min = temporal(minus_i_15min)%bt14%DATA(Elem_Idx,Line_Idx)

            ELSE

               BT_Chn14_15min = Missing_Value_Real4

            ENDIF

            IF (Have_Prev_BT_Chn14_Clr_15min == sym%SUCCESS) THEN
               BT_Chn14_15min_Clr = temporal(minus_i_15min)    &
                                     %bt_clr14%DATA(Elem_Idx,Line_Idx)

            ELSE

               BT_Chn14_15min_Clr = Missing_Value_Real4

            ENDIF

            !---- TERM_THERM_STAB parameters

            IF (Have_Prev_BT_Chn11_1Hr == sym%SUCCESS) THEN
               BT_Chn11_1Hr = temporal(minus_i_01hrs)    &
                                     %bt11%DATA(Elem_Idx,Line_Idx)
            
            ELSE

               BT_Chn11_1Hr = Missing_Value_Real4

            ENDIF


            IF (Have_Prev_BT_Chn14_1Hr == sym%SUCCESS) THEN
               BT_Chn14_1Hr = temporal(minus_i_01hrs)    &
                                     %bt14%DATA(Elem_Idx,Line_Idx)

            ELSE

               BT_Chn14_1Hr = Missing_Value_Real4

            ENDIF

            IF (Have_Prev_BT_Chn15_1Hr == sym%SUCCESS) THEN
               BT_Chn15_1Hr = temporal(minus_i_01hrs)    &
                                     %bt15%DATA(Elem_Idx,Line_Idx)

            ELSE

               BT_Chn15_1Hr = Missing_Value_Real4

            ENDIF


            IF (Have_Prev_Cmask_1Hr == sym%SUCCESS) THEN
               Cmask_1Hr = temporal(minus_i_01hrs)    &
                                     %cldmask%DATA(Elem_Idx,Line_Idx)

            ELSE

               Cmask_1Hr = Missing_Value_Real4

            ENDIF


         ENDIF

         !---------------------------------------------------------------------
         !--- define a coast flag based on nwp tsfc - this indicates
         !--- regions WHERE there is high spatial heterogeneity IN nwp and rtm
         !---------------------------------------------------------------------
         Is_Coast_NWP = sym%NO

         IF (Sfc_Temp_Uni_NWP > Sfc_Temp_Uni_NWP_Thresh) THEN

            Is_Coast_NWP = sym%YES

         ENDIF

         !---------------------------------------------------------------------
         !--- set flags that control how this pixel is processed
         !---------------------------------------------------------------------

         !
         !--- valid pixel - check for aLine_Loc bad channel that is used
         !
         Is_Valid_Pixel = sym%YES

         !
         !--- solar channels (7-16)
         !
         IF (Sol_Zen < Day_Sol_Zen_Thresh) THEN
           DO Chn_Num = 1, 6 

            IF ((out2(Algo_Num)%ch_flg(Chn_Num) > 0).AND.                  &
                (Bad_Pixel_Mask(Elem_Idx,Line_Idx)==sym%YES)) THEN

               Is_Chn(Chn_Num) = sym%NO

            ENDIF

           ENDDO

         ENDIF  

         !
         !--- thermal channels (7-16)
         !
         DO Chn_Num = 7, Nchan_Max

            IF ((out2(Algo_Num)%ch_flg(Chn_Num) > 0).AND.                  &
                (Bad_Pixel_Mask(Elem_Idx,Line_Idx)==sym%YES)) THEN

               Is_Chn(Chn_Num) = sym%NO

            ENDIF

         ENDDO

         !
         !--- Based on Channel Availability, determine validity of pixel
         !
         IF (Is_Chn(14) == sym%NO) THEN
               Is_Valid_Pixel = sym%NO
         ENDIF

         !
         !--- Also check for RTM calculations (use bt14_clr)
         !
         IF (out2(Algo_Num)%ch_flg(14) > 0) THEN
                 IF (BT_Chn14_Clr  < 200.0) then
                      Is_Valid_Pixel = sym%NO
                 ENDIF   
         ENDIF

         !
         !--- store valid Mask IN test_results location
         !
         Test_Results(1,Elem_Idx,Line_Idx) = Is_Valid_Pixel

         !
         !---- if a pixel is invalid, set QF and the skip to next one
         !
         IF (Test_Results(1,Elem_Idx,Line_Idx) == sym%NO) THEN

            Cloud_Mask_QF(Elem_Idx,Line_Idx) = INVALID_CMASK_BAD_CHN14
            
            CYCLE

         ENDIF


         !
         !--- day/night
         !
         Is_Day = sym%NO

         IF (Sol_Zen < Day_Sol_Zen_Thresh) THEN

            Is_Day = sym%YES

         ENDIF

         Test_Results(2,Elem_Idx,Line_Idx) = Is_Day

         !
         !--- terminiator
         !
         Is_Terminator = sym%NO
         if (Sol_Zen > 87.0 .And. Sol_Zen < 93.0) then
                 Is_Terminator = sym%YES
         endif

         Test_Results(3,Elem_Idx,Line_Idx) = Is_Terminator

         !
         !--- land/ocean
         !
         Is_Land = sym%NO

         IF (Land_Type == sym%LAND .or. Land_Type == sym%COASTLINE) THEN

            Is_Land = sym%YES

         ENDIF

         Test_Results(4,Elem_Idx,Line_Idx) = Is_Land

         !
         !--- coast
         !
         Is_Coast = sym%NO
 
         IF (Coast_Type /= sym%NO_COAST) THEN

            Is_Coast = sym%YES

         ENDIF

         Test_Results(5,Elem_Idx,Line_Idx) = Is_Coast

         !
         !--- glint
         !
         Is_Glint = sym%NO

         
         !-- Initial classification based on geometry
         IF ((Is_Day == sym%YES) .AND.  &
             (Is_Land == sym%NO) .AND.  &
             (Glint_Zen < Glint_Zen_Thresh))  THEN  

            Is_Glint = sym%YES

         ENDIF
         
         !-- restore cold pixels IN glint zone to be non-glint
         IF ((Is_Glint == sym%YES) .AND.              &
             ((BT_Chn14 < WATER_FREEZING_POINT) .OR.  &
             (BT_Chn14 < BT_Chn14_Clr - MAX_GLINT_CLR_OBS_BT_CHN14_DIFF))) THEN
             Is_Glint = sym%NO
         ENDIF
         !--- restore pixels with reflectance non-uni
         IF ((Is_Glint == sym%YES) .AND.              &
             (Refl_Chn2_Stddev_3X3(Elem_Idx,Line_Idx) >  &
                 Max_Glint_Clr_Rel_Refl2_Stddev_Thresh * &
                 Refl_Chn2_Mean_3X3(Elem_Idx,Line_Idx))) THEN
             Is_Glint = sym%NO
         ENDIF
   
        
         !--- store result into test vector
         Test_Results(6,Elem_Idx,Line_Idx) = Is_Glint

         !
         !--- desert
         !
         Is_Desert = sym%NO

         IF (Desert_Mask == sym%BRIGHT_DESERT) THEN

            Is_Desert = sym%YES

         ENDIF

         Test_Results(7,Elem_Idx,Line_Idx) = Is_Desert

         !
         !--- snow or ice
         !
         Is_Snow = sym%NO

         IF (Snow_Mask == sym%SNOW .or. Snow_Mask == sym%SEA_ICE) THEN

            Is_Snow = sym%YES

         ENDIF


         !
         !--- add some radiometric checks to prevent false snow
         !

         !--- only do this for nwp snow source - akh
         IF (BT_Chn14 > BT_Chn14_Snow_Thresh) THEN

            Is_Snow = sym%NO

         ENDIF

         Test_Results(8,Elem_Idx,Line_Idx) = Is_Snow

         !
         !--- Test to see if pixel is a cold surface
         !
     
         Is_Cold_Surface = sym%NO
     
         IF (Sfc_Temp < 265.0) THEN
            Is_Cold_Surface = sym%YES
         ENDIF 

         Test_Results(9,Elem_Idx,Line_Idx) = Is_Cold_Surface
         
         
         !--------------------------------------------------------------------
         !--- add single scattering rayleigh reflectance for chn2 clear estimate
         !--- do this only when a dark composite is not used.
         !--------------------------------------------------------------------
         Refl_Sing_Scat = 0.0
         
         ! WCS - clear sky reflectance in CLAVRX is already corrected.
         
!         IF(Sol_Zen < Day_Sol_Zen_Thresh) THEN

!            IF (Is_Land == sym%YES) THEN 
!                Aerosol_Optical_Depth_Chn2 = Aerosol_Optical_Depth_Chn2_Land
!            ELSE
!                Aerosol_Optical_Depth_Chn2 = Aerosol_Optical_Depth_Chn2_Ocean
!            ENDIF
!
!            CALL Compute_Clear_Sky_Scatter(Aerosol_Optical_Depth_Chn2, &
!                                           Aerosol_Single_Scatter_Albedo_Chn2, &
!                                           Aerosol_Asymmetry_Parameter, &
!                                           Rayleigh_Optical_Depth_Chn2, &
!                                           Sat_Name, &
!                                           Total_Precipitable_Water_NWP, &
!                                           Total_Ozone_Path_NWP, &
!                                           Scat_Zen, &
!                                           Cos_Sat_Zen, &
!                                           Cos_Sol_Zen, &
!                                           Refl_Chn2_Clear(Elem_Idx,Line_Idx)/100.0, &
!                                           Refl_Chn2_Clear(Elem_Idx,Line_Idx)/100.0, &
!                                           Transmission_Sing_Scat, &
!                                           Refl_Sing_Scat)

            !--- add it IN to the clear-sky estimate
 !           Refl_Chn2_Clear(Elem_Idx,Line_Idx) = Transmission_Sing_Scat * &
 !                   Refl_Chn2_Clear(Elem_Idx,Line_Idx) + Refl_Sing_Scat
                    
 !           Refl_Chn2_Clear_Max_3x3(Elem_Idx,Line_Idx) = Transmission_Sing_Scat * &
 !                   Refl_Chn2_Clear_Max_3x3(Elem_Idx,Line_Idx) + Refl_Sing_Scat
                    
 !           Refl_Chn2_Clear_Min_3x3(Elem_Idx,Line_Idx) = Transmission_Sing_Scat * &
  !                  Refl_Chn2_Clear_Min_3x3(Elem_Idx,Line_Idx) + Refl_Sing_Scat

            !--- renormalize for improved terminator performance
  !          IF (Sol_Zen > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) THEN
            
  !            Refl_Chn2_Clear(Elem_Idx,Line_Idx) = Term_Refl_Norm(Cos_Sol_Zen, &
  !                                  Refl_Chn2_Clear(Elem_Idx,Line_Idx))
                                    
   !           Refl_Chn2_Clear_Max_3x3(Elem_Idx,Line_Idx) = Term_Refl_Norm(Cos_Sol_Zen, &
   !                                 Refl_Chn2_Clear_Max_3x3(Elem_Idx,Line_Idx))
                                    
   !           Refl_Chn2_Clear_Min_3x3(Elem_Idx,Line_Idx) = Term_Refl_Norm(Cos_Sol_Zen, &
   !                                 Refl_Chn2_Clear_Min_3x3(Elem_Idx,Line_Idx))
                                    
   !         ENDIF

   !      ENDIF

         
         

         !=====================================================================
         ! Determine 0.65 micron clear sky background reflectance
         !=====================================================================
                  RGCT_Threshold = 45.0
                  IF (Is_Land == sym%YES .or. Is_Coast == sym%YES) THEN
                      RGCT_Threshold = 45.0
                      IF (Refl_Chn2_Clear(Elem_Idx, Line_Idx) /= Missing_Value_Real4) then

                           RGCT_Threshold = 10.0 + 1.2 * Refl_Chn2_Clear_Max_3x3(Elem_Idx, Line_Idx) + &
                                            Refl_Chn2_Clear_Stddev_3x3(Elem_Idx,Line_Idx)

                      ENDIF 
                  ENDIF
                  
                  IF (Is_Land == sym%NO) THEN
                      RGCT_Threshold = 99.0
                      IF (Is_Glint == sym%NO) THEN
                         IF (Refl_Chn2_Clear(Elem_Idx, Line_Idx) /= Missing_Value_Real4)  THEN
                                RGCT_Threshold = 5.0 + 1.2 * Refl_Chn2_Clear_Max_3x3(Elem_Idx, Line_Idx) 
                         ENDIF
                       ENDIF
                  ENDIF

         !=====================================================================
         ! Apply Clear Spatial Uniformity Tests
         !=====================================================================

         !
         !--- RUT 0.65 micron clear reflectance uniformity test
         !
         Num_Tests = 1 + Num_Ancil_Tests
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(2) == sym%YES) THEN

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RUT_Routine (&
                              Is_Land, &
                              Is_Snow, &
                              Is_Coast, &
                              Refl_Chn2_Clear(Elem_Idx,Line_Idx), &
                              Refl_Chn2_Stddev_3x3(Elem_Idx,Line_Idx), &
                              Sol_Zen)

          ENDIF 

         !
         !--- 11 micron clear uniformity test
         !
         Num_Tests = 2 + Num_Ancil_Tests
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(14) == sym%YES) THEN

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = TUT_Routine (&
                              Is_Land, &
                              Is_Coast, &
                              BT_Chn14_Stddev_3x3(Elem_Idx,Line_Idx), &
                              Sfc_Hgt_Stddev_3x3(Elem_Idx,Line_Idx))

         ENDIF

         !=====================================================================
         ! Apply Cloud Detection Tests
         !=====================================================================

         !--------------------------------------------------------------------
         ! infrared cloud tests
         !--------------------------------------------------------------------

         !--------------------------------------------------------------------
         !--- Relative Thermal Contrast Test - RTCT
         !--------------------------------------------------------------------
         Num_Tests = First_IR_Cld_Mask_Test + 0
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(14) == sym%YES) THEN

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RTCT_Routine(&
                              Is_Land,  &
                              Is_Coast,  &
                              Is_Snow, &
                              Is_Cold_Surface, &
                              BT_Chn14, &
                              BT_Chn14_Min_3x3(Elem_Idx,Line_Idx), &
                              BT_Chn14_Max_3x3(Elem_Idx,Line_Idx), &
                              Sfc_Hgt_Stddev_3x3(Elem_Idx,Line_Idx))
         ENDIF

         !--------------------------------------------------------------------
         !--- 11 micron Tropospheric Emissivity Test - ETROP
         !--------------------------------------------------------------------
         Num_Tests = First_IR_Cld_Mask_Test + 1
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(14) == sym%YES) THEN

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = ETROP_Routine( &
                               Is_Snow, &
                               Is_Land, &
                               Is_Coast, &
                               Is_Desert, &
                               Is_Cold_Surface, &
                               Land_Type, &
                               BT_Chn14, &
                               BT_Chn14_Clr, &
                               Emiss_Tropo_Chn14(Elem_Idx,Line_Idx), &
                               Emiss_Tropo_Chn14_LRC, &
                               BT_Chn14_Stddev_3x3(Elem_Idx,Line_Idx))
         ENDIF

         !--------------------------------------------------------------------
         !--- Positive 11-12 micron Test 
         !--------------------------------------------------------------------
         Num_Tests = First_IR_Cld_Mask_Test + 2
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF ((Is_Chn(14) == sym%YES) .and. (Is_Chn(15) == sym%YES)) THEN
            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = PFMFT_Routine( &
                      Is_Land, &
                      Is_Snow, &
                      Is_Cold_Surface, &
                      BT_Chn14, &
                      BT_Chn15, &
                      BT_Chn14_Clr, &
                      BT_Chn15_Clr, &
                      BT_Chn14_Stddev_3x3(Elem_Idx,Line_Idx)) 
         ENDIF

         !------------------------------------------------------------
         !--- NFMFT - Negative 11-12 micron Test
         !------------------------------------------------------------
         Num_Tests = First_IR_Cld_Mask_Test + 3
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF ((Is_Chn(14) == sym%YES) .and. (Is_Chn(15) == sym%YES)) THEN
             Test_Results(Num_Tests,Elem_Idx,Line_Idx) = NFMFT_Routine( &
                      Is_Land, &
                      Is_Snow, &
                      Is_Desert, &
                      BT_Chn14, &
                      BT_Chn15, &
                      BT_Chn14_Clr, &
                      BT_Chn15_Clr)
         ENDIF

         !---------------------------------------------------------------------
         !---  Relative FMFT
         !---------------------------------------------------------------------
         Num_Tests = First_IR_Cld_Mask_Test + 4
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF ((Is_Chn(14) == sym%YES) .and. (Is_Chn(15) == sym%YES)) THEN

           Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RFMFT_Routine( &
                      Is_Land, &
                      Is_Coast, &
                      BT_Chn14, &
                      BT_Chn15, &
                      BTD_Chn14_Chn15_NWC)
         ENDIF


         !---------------------------------------------------------------------
         !---  Cirrus Water Vapor Test (CIRH2O)
         !---------------------------------------------------------------------
         Num_Tests = First_IR_Cld_Mask_Test + 5
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (((Is_Chn(10) == sym%YES) .or. (Is_Chn(9) == sym%YES)) .and.  &
             (Is_Chn(14) == sym%YES)) THEN
             
            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = CIRH2O_Routine( &
                      Sfc_Hgt, &
                      BT_WV_BT_Window_Corr(Elem_Idx,Line_Idx), &
                      BT_WaterVapor_Stddev_3x3(Elem_Idx,Line_Idx), &
                      BT_Chn14_Stddev_3x3(Elem_Idx,Line_Idx), &
                      Total_Precipitable_Water_NWP, &
                      Cos_Sat_Zen)


         ENDIF

         !---------------------------------------------------------------------
         !---  Temporal IR test
         !---------------------------------------------------------------------
         Num_Tests = First_IR_Cld_Mask_Test + 6
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(14) == sym%YES) THEN

               Test_Results(Num_Tests,Elem_Idx,Line_Idx) = TEMPIR_Routine( &
                          BT_Chn14, &
                          BT_Chn14_Clr, &
                          BT_Chn14_15min, &
                          BT_Chn14_15min_Clr)
         ENDIF

         !---------------------------------------------------------------------
         !---  Terminator Temporal IR test
         !---------------------------------------------------------------------
         Num_Tests = First_IR_Cld_Mask_Test + 7
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(11) == sym%YES .AND. &
             Is_Chn(14) == sym%YES .AND. &
             Is_Chn(15) == sym%YES .AND. &
             BT_Chn11 /= Missing_Value_Real4 .AND. &
             BT_Chn14 /= Missing_Value_Real4 .AND. &
             BT_Chn15 /= Missing_Value_Real4 .AND. &
             Have_Prev_BT_Chn11_1Hr == sym%SUCCESS .AND. &
             Have_Prev_BT_Chn14_1Hr == sym%SUCCESS .AND. &
             Have_Prev_BT_Chn15_1Hr == sym%SUCCESS .AND. &
             Have_Prev_Cmask_1Hr == sym%SUCCESS .AND. &
             BT_Chn11_1Hr /= Missing_Value_Real4 .AND. &
             BT_Chn14_1Hr /= Missing_Value_Real4 .AND. &
             BT_Chn15_1Hr /= Missing_Value_Real4 .AND. &
             Cmask_1Hr /= Missing_Value_Real4) THEN

             Test_Results(Num_Tests,Elem_Idx,Line_Idx) = Term_Therm_Stab_Routine(&
                              Is_Land, &
                              Sol_Zen, &
                              BT_Chn11, &
                              BT_Chn11_1Hr, &
                              BT_Chn14, &
                              BT_Chn14_1Hr, &
                              BT_Chn15, &
                              BT_Chn15_1Hr, &
                              Cmask_1Hr)     
         ENDIF


         !=====================================================================
         ! Solar Reflectance Cloud Tests
         !=====================================================================

         !
         !--- Reflectance Gross Contrast Test RGCT
         !
         Num_Tests = First_Vis_Cld_Mask_Test + 0
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(2) == sym%YES) THEN

               Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RGCT_Routine( &
                            Is_Snow, &
                            Is_Glint, &
                            Sol_Zen, &
                            Refl_Chn2, &
                            RGCT_Threshold)

         ENDIF

         !----------------------------------------------------------------------
         !---  Relative Visible Contrast Test - RVCT
         !----------------------------------------------------------------------
         Num_Tests = First_Vis_Cld_Mask_Test + 1
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(2) == sym%YES) THEN

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RVCT_Routine( &
                         Is_Coast, &
                         Is_Snow, &
                         Is_Land, &
                         Refl_Chn2, &
                         Scat_Zen, &
                         Sol_Zen, &
                         Refl_Chn2_Clear_Stddev_3x3(Elem_Idx,Line_Idx), &
                         Refl_Chn2_Min_3x3(Elem_Idx,Line_Idx))
         ENDIF

         !=====================================================================
         ! SWIR Solar Cloud Tests
         !=====================================================================

         !---------------------------------------------------------------------
         !--- Near Infrared Reflectance Test (NIRREF)
         !---------------------------------------------------------------------
         Num_Tests = First_SWIR_Solar_Cld_Mask_Test + 0
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF ((Is_Chn(5) == sym%YES) .AND. (Is_Chn(2) == sym%YES)) THEN

            NDSI = (Refl_Chn2 - Refl_Chn5)/(Refl_Chn5 + Refl_Chn2)

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = NIRREF_Chn5_Routine(&
                         Is_Coast, &
                         Is_Snow, &
                         Sol_Zen, &
                         Sfc_Hgt, &
                         NDSI, &
                         Refl_Chn5)

         ENDIF

         IF ((Is_Chn(7) == sym%YES) .and. (Is_Chn(5) == sym%NO)) THEN

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = NIRREF_Chn7_Routine( &
                         Is_Snow, &
                         Sol_Zen, &
                         Sfc_Hgt, &
                         Refl_Chn7)

         ENDIF
 
         !---------------------------------------------------------------------
         !--- Cirrus Reflectance Test (CIRREF)
         !---------------------------------------------------------------------
         Num_Tests = First_SWIR_Solar_Cld_Mask_Test + 1
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(4) == sym%YES) THEN

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = CIRREF_Routine( &
                         Is_Snow, &
                         Refl_Chn4, &
                         Sfc_Hgt_Max_3x3(Elem_Idx,Line_Idx), &
                         Sol_Zen)

         ENDIF

         !=====================================================================
         ! SWIR Solar Thermal Tests
         !=====================================================================

         !---------------------------------------------------------------------
         !--- 3.9 micron emissivity test - EMISS_4
         !---------------------------------------------------------------------
         Num_Tests = First_SWIR_Thermal_Cld_Mask_Test + 0
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF ((Is_Chn(7) == sym%YES) .AND. &
             (Emiss_Chn7 /= Missing_Value_Real4) .AND. &
             (Rad_Chn7_Clr /= Missing_Value_Real4)) THEN

            !
            !--- make estimate of clear emissivity that includes solar
            !
            Planck_Emission_Chn7_Clr = planck_rad_fast(7,BT_Chn14_Clr)
            Solar_Rad_Chn7_Clr = Rad_Chn7_Clr

            IF (Is_Day == sym%YES .or. Is_Terminator == sym%YES) THEN
              Solar_Rad_Chn7_Clr = Rad_Chn7_Clr + (1.0 - Sfc_Emiss_Chn7_RTM)         &
                     * Atm_Solar_Trans_Chn7 * max(Cos_Sol_Zen,0.05) * (Chn7_Sol_Energy/PI)
            ENDIF

            Solar_BT_Chn7_Clr = planck_temp_fast(7,Solar_Rad_Chn7_Clr)

            Emiss_Chn7_Clr = Solar_Rad_Chn7_Clr / Planck_Emission_Chn7_Clr

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = EMISS4_Routine(Is_Glint,  &
                         Is_Land, &
                         Is_Snow, &
                         Is_Desert, &
                         BT_Chn14, &
                         Emiss_Chn7,  &
                         Emiss_Chn7_Clr, &
                         Sfc_Emiss_Chn7_RTM)
         ENDIF
         
         

         !---------------------------------------------------------------------
         !--- Uniform Low Stratus Test - ULST
         !---------------------------------------------------------------------
         Num_Tests = First_SWIR_Thermal_Cld_Mask_Test + 1
         Test_Results(Num_Tests,Elem_Idx,Line_Idx) = sym%NO

         IF (Is_Chn(7) == sym%YES) THEN

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = ULST_Routine(&
                       Is_Land, &
                       Is_Day, &
                       Is_Snow, &
                       Is_Cold_Surface, &
                       BT_Chn14,  &
                       Emiss_Chn7,  &
                       Emiss_Chn7_Clr, &
                       Emiss_Chn7_NWC, &
                       Sfc_Emiss_Chn7_RTM)

        ENDIF

        !======================================================================
        ! Combine each test into final 2 bit cloud Mask(0,1,2,3)
        !======================================================================

         !
         !--- full cloud Mask: vis + ir + swir tests
         !
         IF (SUM(Test_Results(First_Cld_Mask_Test:Last_Cld_Mask_Test,   &
            Elem_Idx,Line_Idx)) == Num_Cld_Mask_Tests*sym%NO) THEN  

            !
            !---check cld tests
            !
            IF (SUM(Test_Results(First_Clr_Uni_Test:Last_Clr_Uni_Test,  &
               Elem_Idx,Line_Idx)) == Num_Clr_Uni_Tests*sym%NO) THEN   

               !
               !---check uni tests
               !
               Cloud_Mask(Elem_Idx,Line_Idx) = sym%CLEAR

            ELSE

               Cloud_Mask(Elem_Idx,Line_Idx) = sym%PROB_CLEAR

            ENDIF

         ELSE

            Cloud_Mask(Elem_Idx,Line_Idx) = sym%CLOUDY

         ENDIF

         !
         !--- ir tests
         !
         IF ((SUM(Test_Results(First_IR_Cld_Mask_Test:Last_IR_Cld_Mask_Test,Elem_Idx,Line_Idx))  + &
              SUM(Test_Results(First_SWIR_Thermal_Cld_Mask_Test:Last_SWIR_Thermal_Cld_Mask_Test,Elem_Idx,Line_Idx)) )  &
              == (Num_IR_Cld_Mask_Tests + Num_SWIR_Thermal_Cld_Mask_Tests)*sym%NO) THEN  

         !
         !---check cld tests
         !

            IF (SUM(Test_Results(First_Clr_Uni_Test:Last_Clr_Uni_Test,      &
                Elem_Idx,Line_Idx)) == Num_Clr_Uni_Tests*sym%NO) THEN   

               !
               !---check uni tests
               !
               Cloud_Mask_IR(Elem_Idx,Line_Idx) = sym%CLEAR

            ELSE

               Cloud_Mask_IR(Elem_Idx,Line_Idx) = sym%PROB_CLEAR

            ENDIF

         ELSE

            Cloud_Mask_IR(Elem_Idx,Line_Idx) = sym%CLOUDY

         ENDIF

         !-----------------------------------------------------
         ! Assign a Quality Flag
         !----------------------------------------------------
         IR_Test_Sum = SUM(Test_Results(First_IR_Cld_Mask_Test:Last_IR_Cld_Mask_Test,Elem_Idx,Line_Idx))
         VIS_Test_Sum = SUM(Test_Results(First_VIS_Cld_Mask_Test:Last_VIS_Cld_Mask_Test,Elem_Idx,Line_Idx))
         SWIR_SOLAR_TEST_Sum = SUM(Test_Results(First_SWIR_SOLAR_Cld_Mask_Test: &
                                             Last_SWIR_SOLAR_Cld_Mask_Test,Elem_Idx,Line_Idx))
         SWIR_THERMAL_TEST_Sum = SUM(Test_Results(First_SWIR_THERMAL_Cld_Mask_Test: &
                                               Last_SWIR_THERMAL_Cld_Mask_Test,Elem_Idx,Line_Idx))

         Cloud_Mask_QF(Elem_Idx,Line_Idx) = IR_Test_Sum + VIS_Test_Sum +  &
                                            SWIR_Solar_Test_Sum + SWIR_Thermal_Test_Sum

         IR_Test_Mask = max(1,IR_Test_Sum)
         VIS_Test_Mask = max(1,VIS_Test_Sum)
         SWIR_Test_Mask = max(1,SWIR_THERMAL_Test_Sum + SWIR_SOLAR_Test_Sum)

         ! Set QF's based on available information
        
         IF ((Sensor%Chan_On_Flag_Default(20) > 0) .AND. (Is_Chn(7) == sym%NO)) THEN
             
             Cloud_Mask_QF(Elem_Idx,Line_Idx) = REDUCED_QUAL_BAD_CHN7
         
         ELSE IF ((Sensor%Chan_On_Flag_Default(1) > 0) .AND. &
                  (Is_Day == sym%YES) .AND. (Is_Chn(2) == sym%NO)) THEN
         
             Cloud_Mask_QF(Elem_Idx,Line_Idx) = REDUCED_QUAL_BAD_CHN2
             
          ELSE IF (((Is_Day == sym%YES) .AND. &
                    (((Sensor%Chan_On_Flag_Default(26) > 0) .AND. (Is_Chn(4) == sym%NO)) .OR. &
                    ((Sensor%Chan_On_Flag_Default(6) > 0) .AND. (Is_Chn(5) == sym%NO))))  .OR. &
                    ((Sensor%Chan_On_Flag_Default(28)  > 0) .AND. (Is_Chn(10) == sym%NO))  .OR. &             
                    ((Sensor%Chan_On_Flag_Default(29) > 0) .AND. (Is_Chn(11) == sym%NO))  .OR. &             
                    ((Sensor%Chan_On_Flag_Default(32) > 0) .AND. (Is_Chn(15) == sym%NO))) THEN

             Cloud_Mask_QF(Elem_Idx,Line_Idx) = REDUCED_QUAL_BAD_OTHER
 
         ELSE
         
             Cloud_Mask_QF(Elem_Idx,Line_Idx) = VALID_CMASK_RETRIEVAL
             
         ENDIF
                  
      END DO Element_Loop_2
   END DO Line_Loop_2
   
   !======================================================================
   ! Apply Restoral Tests Here
   !======================================================================

   !----------------------------------------------------------------------
   ! Probably Clear Restoral Test
   !----------------------------------------------------------------------

   Num_Tests = First_Res_Test + 0

   Num_Pix = 2 !pixel radius of window

   Cloud_Mask_Tmpy = Cloud_Mask

   CALL Compute_Probably_Clear_Restoral(Cloud_Mask_Tmpy,Cloud_Mask,   &
                                        Test_Results(1,:,:), &
                                        Test_Results(Num_Tests,:,:),Num_Pix)

   Cloud_Mask_Tmpy = 0


   !======================================================================
   ! Probably Cloudy Values
   !======================================================================
   Num_Tests = First_Res_Test + 1

   !
   !---pixel radius of window
   !

   Num_Pix = 1    

   Cloud_Mask_Tmpy = Cloud_Mask

   CALL Compute_Probably_Cloudy(Cloud_Mask_Tmpy,Cloud_Mask,  &
                                Test_Results(Num_Tests,:,:),Num_Pix)

   Cloud_Mask_Tmpy = 0


   !----------------------------------------------------------------------
   ! Probably Cloudy Restoral Test
   !----------------------------------------------------------------------

   Num_Tests = First_Res_Test + 2

   Num_Pix = 5 !pixel radius of window

   Cloud_Mask_Tmpy = Cloud_Mask

   Cloud_Mask_Tmpy = 0

   !======================================================================
   ! Pack test results into bytes for Output
   !======================================================================

   Line_Loop_3: DO Line_Idx=1, Image%Number_Of_Lines_Read_This_Segment

      Element_Loop_3: DO Elem_Idx=1, Num_Elem

         CALL PACK_BYTES(Test_Results(1:Total_Num_Tests,Elem_Idx,Line_Idx), &
                                   Test_Bit_Depth(1:Total_Num_Tests), &
                                   Cloud_Mask_Packed(1:5,Elem_Idx,Line_Idx))


      END DO Element_Loop_3

   END DO Line_Loop_3

   !======================================================================
   ! Make Binary Cloud Mask - WCS 03052012
   !======================================================================

   Line_Loop_4: DO Line_Idx=1, Image%Number_Of_Lines_Read_This_Segment

      Element_Loop_4: DO Elem_Idx=1, Num_Elem

         IF ((Cloud_Mask(Elem_Idx,Line_Idx) .EQ. sym%PROB_CLOUDY) .OR. &
             (Cloud_Mask(Elem_Idx,Line_Idx) .EQ. sym%CLOUDY)) THEN
             
             Cloud_Mask_Binary(Elem_Idx,Line_Idx) = CLD_MASK_BIN_CLOUD
             
         ENDIF


      END DO Element_Loop_4

   END DO Line_Loop_4

   !======================================================================
   ! nullify local pointers
   !======================================================================

   Cloud_Mask => null()
   Cloud_Mask_Binary => null()
   Cloud_Mask_IR => null()
   Cloud_Mask_Packed => null()
   Test_Results => null()
   Cloud_Mask_QF => null()
   Cloud_Mask_Tmpy => null()
   LRC_Mask => null()
   X_LRC_Idx => null()
   Y_LRC_Idx => null()
   Emiss_Tropo_Chn14 => null()
   !======================================================================
   ! deallocate
   !======================================================================

   CALL Destroy_Spatial_Uniformity(BT_WaterVapor_Mean_3x3, BT_WaterVapor_Max_3x3,  &
                                   BT_WaterVapor_Min_3x3, BT_WaterVapor_Stddev_3x3)
   CALL Destroy_Spatial_Uniformity(BT_Chn14_Mean_3x3, BT_Chn14_Max_3x3,  &
                                   BT_Chn14_Min_3x3, BT_Chn14_Stddev_3x3)
   CALL Destroy_Spatial_Uniformity(Refl_Chn2_Mean_3X3, Refl_Chn2_Max_3x3,  &
                                   Refl_Chn2_Min_3X3, Refl_Chn2_Stddev_3X3)
   CALL Destroy_Spatial_Uniformity(Sfc_Hgt_Mean_3X3, Sfc_Hgt_Max_3X3,  &
                                   Sfc_Hgt_Min_3X3, Sfc_Hgt_Stddev_3X3)

   !Initialize total deallocation status 
   Alloc_Status_Total = 0 
   
   IF (ALLOCATED(X_NWC_Idx)) DEALLOCATE(X_NWC_Idx,stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status
   
   IF (ALLOCATED(Y_NWC_Idx)) DEALLOCATE(Y_NWC_Idx,stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status
   
   IF (ALLOCATED(BT_WV_BT_Window_Corr)) DEALLOCATE(BT_WV_BT_Window_Corr,stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status

   !======================================================================
   ! deallocate 2nd darkest composite data array
   !======================================================================
  ! IF (ALLOCATED(Refl_Chn2_Clear)) DEALLOCATE(Refl_Chn2_Clear,stat=Alloc_Status)
   !Alloc_Status_Total = Alloc_Status_Total + Alloc_Status

   Refl_Chn2_Clear => null()
   
   IF (Alloc_Status /= 0) THEN

       WRITE (Err_Message, *) &
             'Error deallocating Chn32 Clear-sky Reflectance'
             
             
        Error_Level = 2 ! AIT FATAL ERROR CODE
            
        CALL Display_Message("Baseline Cloud Mask", &
                 TRIM(Err_Message), &
                 Sym%FAILURE)
                 
        !  AIT Error Messaging
        !  CALL Error_Messaging (Routine_Name, Error_Message, Error_Level)           
        RETURN      
      
   ENDIF


   !======================================================================
   ! Check deallocation of arrays
   !======================================================================

   IF (Alloc_Status_Total /= 0) THEN

        WRITE (Err_Message, *) &
            'Error deallocating arrays'
        Error_Level = 2 ! AIT FATAL ERROR CODE
            
        CALL Display_Message("Baseline Cloud Mask", &
                 TRIM(Err_Message), &
                 Sym%FAILURE)
                 
        !  AIT Error Messaging
        !  CALL Error_Messaging (Routine_Name, Error_Message, Error_Level)           
        RETURN      

   ENDIF
   

!
!---END SUBROUTINE Baseline_Cloud_Mask_Main
!


END SUBROUTINE Baseline_Cloud_Mask_Main



!====================================================================
! Subroutine Name: Compute_Probably_Clear_Restoral
!
! Function:
!    Routine to restore probably clear to clear
!
! Description:
!   This subroutine restores probably clear pixels to clear.
!
! Calling Sequence:
!   CALL Compute_Probably_Clear_Restoral(input, Output, Valid_Mask, Mask,Num_Pix)
!
! Inputs: 
!   input - the unrestored array
!   Valid_Mask - flag telling which elements have a valid cloud Mask
!   Num_Pix - radius of pixel window to search for a cloudy result
! 
! Outputs: 
!   Output - the restored array
!   Mask - flag telling which elements of input were restored
!
! Dependencies: None
!
! Restrictions:  None
!
!====================================================================
SUBROUTINE Compute_Probably_Clear_Restoral(input, Output, Valid_Mask, Mask,Num_Pix)
   INTEGER(KIND=INT1), INTENT(IN), DIMENSION(:,:):: input
   INTEGER(KIND=INT1), INTENT(IN), DIMENSION(:,:):: Valid_Mask
   INTEGER(KIND=INT1), INTENT(OUT), DIMENSION(:,:):: Output
   INTEGER(KIND=INT1), INTENT(OUT), DIMENSION(:,:):: Mask
   INTEGER, INTENT(IN):: Num_Pix
   INTEGER:: Elem_Idx  !index for pixel IN east-west direction
   INTEGER:: Line_Idx  !index for pixel IN north-south direction
   INTEGER:: Array_Right  !temporary index IN ielem-direction
   INTEGER:: Array_Left  !temporaty index IN ielem-direction
   INTEGER:: Array_Top  !temporary index IN iline-direction
   INTEGER:: Array_Bottom  !temporary index IN iline-direction

   !
   !--- initialize Output
   !
   Output = input
   Mask = sym%NO

   !
   !--- loop over scan lines IN segment
   !
   Line_Loop: DO Line_Idx=1, Image%Number_Of_Lines_Read_This_Segment

      !
      !--- determine y-dimensions of array to check
      !
      Array_Top = max(1,Line_Idx-Num_Pix)
      Array_Bottom = min(Image%Number_Of_Lines_Read_This_Segment,Line_Idx+Num_Pix)

      Element_Loop: DO Elem_Idx = 1, Image%Number_Of_Elements

         !
         !--- check to see if Mask is valid, if not cycle to next
         !
         IF (Valid_Mask(Elem_Idx,Line_Idx) == sym%NO) THEN
                 CYCLE
         ENDIF

         !
         !--- determine x-dimensions of array to check
         !
         Array_Right = max(1,Elem_Idx-Num_Pix)
         Array_Left = min(Image%Number_Of_Elements,Elem_Idx+Num_Pix)

         IF (input(Elem_Idx,Line_Idx) == sym%PROB_CLEAR) THEN


            IF ( (sym%CLEAR < sym%CLOUDY) .AND.  &
                 ( (MAXVAL(input(Array_Right:Array_Left,Array_Top:Array_Bottom)) /= sym%CLOUDY) .AND. &
                   (MAXVAL(input(Array_Right:Array_Left,Array_Top:Array_Bottom)) /= sym%PROB_CLOUDY))) THEN

               Mask(Elem_Idx,Line_Idx) = sym%YES
               Output(Elem_Idx,Line_Idx) = sym%CLEAR

            ENDIF

            IF ( (sym%CLEAR > sym%CLOUDY) .AND.  &
                 ( (MINVAL(input(Array_Right:Array_Left,Array_Top:Array_Bottom)) /= sym%CLOUDY) .AND. &
                   (MINVAL(input(Array_Right:Array_Left,Array_Top:Array_Bottom)) /= sym%PROB_CLOUDY))) THEN

               Mask(Elem_Idx,Line_Idx) = sym%YES
               Output(Elem_Idx,Line_Idx) = sym%CLEAR

            ENDIF

         ENDIF

      END DO Element_Loop

   END DO Line_Loop

!
! end SUBROUTINE Compute_Probably_Clear_Restoral
!
END SUBROUTINE Compute_Probably_Clear_Restoral

!====================================================================
! Subroutine Name: Compute_Probably_Cloudy
!
! Function:
!    Routine to compute probably cloudy pixels
!
! Description:
!   This subroutine computes which pixels are probably cloudy.
!
! Calling Sequence:
!   CALL Compute_Probably_Cloudy(input,Output,Mask,Num_Pix)
!
! Inputs: 
!   input - the unrestored array
!   Num_Pix - radius of pixel window to search for a cloudy result
! 
! Outputs: 
!   Output - the restored array
!   Mask - flag telling which elements of input were restored
!
! Dependencies: None
!
! Restrictions:  None
!
!====================================================================


SUBROUTINE  Compute_Probably_Cloudy(input,Output,Mask,Num_Pix)
   INTEGER(KIND=INT1), INTENT(IN), DIMENSION(:,:):: input
   INTEGER(KIND=INT1), INTENT(OUT), DIMENSION(:,:):: Output
   INTEGER(KIND=INT1), INTENT(OUT), DIMENSION(:,:):: Mask
   INTEGER, INTENT(IN):: Num_Pix
   INTEGER:: Elem_Idx
   INTEGER:: Line_Idx
   INTEGER:: Array_Right
   INTEGER:: Array_Left
   INTEGER:: Array_Top
   INTEGER:: Array_Bottom
   INTEGER:: Num_Cloudy
   INTEGER:: Num_Cloudy_Max



   !--------------------------------------------------------------------------
   ! Executable Code
   !--------------------------------------------------------------------------

   !
   !--- initialize Output
   !
   Output = input
   Mask = sym%NO

   !
   !--- loop over scan lines IN segment
   !
   Line_Loop: DO Line_Idx=1, Image%Number_Of_Lines_Read_This_Segment

      !
      !--- determine y-dimensions of array to check
      !
      Array_Top = max(1,Line_Idx-Num_Pix)
      Array_Bottom = min(Image%Number_Of_Lines_Read_This_Segment,Line_Idx+Num_Pix)

      Element_Loop: DO Elem_Idx=1, Image%Number_Of_Elements

         !
         !--- determine x-dimensions of array to check
         !
         Array_Right = max(1,Elem_Idx-Num_Pix)
         Array_Left = min(Image%Number_Of_Elements,Elem_Idx+Num_Pix)

         !
         !--- check to see if a cloudy pixel neighbors a non-cloudy pixel
         !
         IF (input(Elem_Idx,Line_Idx) == sym%CLOUDY) THEN

            Num_Cloudy = SUM(input(Array_Right:Array_Left,Array_Top:Array_Bottom))
            
            Num_Cloudy_Max = sym%CLOUDY * (Array_Left-Array_Right+1) * &
                             (Array_Bottom - Array_Top+1)

            !
            !--- if the value is less, set to probably cloudy
            !

            IF (Num_Cloudy /= Num_Cloudy_Max) THEN

               Mask(Elem_Idx,Line_Idx) = sym%YES
               Output(Elem_Idx,Line_Idx) = sym%PROB_CLOUDY

            ENDIF

         ENDIF

      END DO Element_Loop

   END DO Line_Loop
!
!end SUBROUTINE Compute_Probably_Cloudy
!
END SUBROUTINE Compute_Probably_Cloudy


!====================================================================
! Subroutine Name: Clear_Chn2_Reflectance_Field
!
! Function:
!    Routine to Read in Chn2 Reflectance Fields
!
! Description:
!   This subroutine reads in the channel 2 reflectance fields for a given 
!     segment
!
! Calling Sequence:
!   CALL Clear_Chn2_Reflectance_Field(Num_Elem, &
!                                     Max_Num_Lines_per_Seg, &
!                                     Refl_Chn2_Clear)
!
! Inputs: 
!  Num_Elem - Number of elements in this segment of data
!  Max_Num_Lines_per_Seg - Number of lines in this segment of data
! 
! Outputs: 
!  Refl_Chn2_Clear - Clear sky channel 2 albedo
!
! Dependencies: None
!
! Restrictions:  None
!
!====================================================================
! WCS - commented out this routine since CLAVR-x does this already

! SUBROUTINE Clear_Chn2_Reflectance_Field(Num_Elem, &
!                                         Max_Num_Lines_per_Seg, &
!                                         Refl_Chn2_Clear)
!
!   REAL(KIND=REAL4), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:):: Refl_Chn2_Clear
!   INTEGER(KIND=INT4),INTENT(IN) :: Num_Elem
!   INTEGER(KIND=INT4),INTENT(IN) :: Max_Num_Lines_per_Seg


   !--- If create one and populate background reflectance using default or global WS albedo data
!       IF (.NOT. ALLOCATED(Refl_Chn2_Clear)) ALLOCATE( &
!                           Refl_Chn2_Clear(Num_Elem,Max_Num_Lines_per_Seg))

!       Refl_Chn2_Clear = 5.0
!       WHERE(Sfc%Land == sym%LAND)
!          Refl_Chn2_Clear = 45.0
!       ENDWHERE 
!       !--- use MODIS white sky from mapped data if available
!       WHERE(sat%ws_albedo2 /= Missing_Value_Real4)
!          Refl_Chn2_Clear = sat%ws_albedo2
!       ENDWHERE 
!       !--- ensure consistent missing value
!       WHERE(Refl_Chn2_Clear <= 0.0)
!          Refl_Chn2_Clear = Missing_Value_Real4
!       ENDWHERE
!
! END SUBROUTINE Clear_Chn2_Reflectance_Field


!====================================================================
! Subroutine Name: Compute_Emiss_Tropo_Chn14
!
! Function:
!   Computes the 11 micron emissivity at Tropopause
!
! Description:
!   This subroutine computes 11 micron emissivity at Tropopause for a segment of
!   data
!
! Calling Sequence:
!   CALL Compute_Emiss_Tropo_Chn14(Emiss_Tropo_Chn14, &
!                                   Number_of_Lines_in_this_Segment)
!
! Inputs: 
!  Number_of_Lines_in_this_Segment - Number of lines in this segment of data
! 
! Outputs: 
!  Emiss_Tropo_Chn14 - Channel 14 emissivity at the tropopause
!
! Dependencies: None
!
! Restrictions:  None
!
!====================================================================
 SUBROUTINE Compute_Emiss_Tropo_Chn14(Emiss_Tropo_Chn14,&
                                      Number_of_Lines_in_this_Segment)
   INTEGER(KIND=INT4),INTENT(IN) :: Number_of_Lines_in_this_Segment
   REAL(KIND=REAL4), DIMENSION(:,:), INTENT(OUT):: Emiss_Tropo_Chn14
   INTEGER(KIND=INT1):: Tropo_Idx_NWP
   INTEGER(KIND=INT1):: View_Zen_Idx
   INTEGER:: X_NWP_Idx
   INTEGER:: Y_NWP_Idx
   INTEGER:: Elem_Idx
   INTEGER:: Line_Idx
   REAL(KIND=REAL4) :: Rad_Chn14
   REAL(KIND=REAL4) :: Clr_Rad_Chn14
   REAL(KIND=REAL4) :: Blkbdy_Tropo_Rad_Chn14

   !--- initialize
   Emiss_Tropo_Chn14 = Missing_Value_Real4

    Line_Loop: DO Line_Idx=1, Number_of_Lines_in_this_Segment
      Element_Loop: DO Elem_Idx = 1, Image%Number_Of_Elements

       IF (Space_Mask(Elem_Idx,Line_Idx) == sym%NO) THEN

            !
            !---nwp longitude cell
            !
            X_NWP_Idx =          I_Nwp(Elem_Idx,Line_Idx)         

            !
            !---nwp latitude cell
            !
            Y_NWP_Idx =          J_Nwp(Elem_Idx,Line_Idx)     

            !
            !---nwp level associated with tropopause
            !
            Tropo_Idx_NWP =        Rtm(X_NWP_Idx,Y_NWP_Idx)%Tropo_Level 

            !
            !---viewing zenith angle bin
            !
            View_Zen_Idx =          Zen_Idx_Rtm(Elem_Idx,Line_Idx)        

            !
            !---11 um radiance
            !
            Rad_Chn14  =        ch(31)%Rad_Toa(Elem_Idx,Line_Idx)

            !
            !---clear 11 micron radiance
            !
            Clr_Rad_Chn14 =     ch(31)%Rad_Toa_Clear(Elem_Idx,Line_Idx)

            !
            !---BB 11 um rad at tropopause
            !
            Blkbdy_Tropo_Rad_Chn14 = rtm(X_NWP_Idx,Y_NWP_Idx)%d(View_Zen_Idx)%ch(31)%cloud_prof14(Tropo_Idx_NWP)

            !
            !---Tropopause Emissivity
            !
            Emiss_Tropo_Chn14(Elem_Idx,Line_Idx) =  &
                  (Rad_Chn14 - Clr_Rad_Chn14) / (Blkbdy_Tropo_Rad_Chn14 - Clr_Rad_Chn14)

      END IF
    END DO Element_Loop
  END DO Line_Loop

 END SUBROUTINE Compute_Emiss_Tropo_Chn14

!====================================================================
! Subroutine Name: Compute_NWC
!
! Function:
!   Routine to compute spatial uniformity for every pixel in an array
! based on its nxn neighborhood
!
! Description:
!   This subroutine computes the nearest warm center for a given array of BT's 
!
! Calling Sequence:
!   CALL Compute_NWC( Input_Array, &
!                        Box_Size, &
!                        Uni_Land_Mask_Flag, &
!                        Bad_Mask, &
!                        Land_Mask,  &
!                        Elem_Start, &
!                        Elem_End, &
!                        Line_Start, &
!                        Line_End, &
!                        Loc_Max_Elem,  &
!                        Loc_Max_Line)
!
! Inputs: 
! Input_Array - the input array
! Box_Size - the size of the box (1=3x3, 2=5x5, ...)
! Bad_Mask - bad pixel mask array (only values with sym%YES will contribute)
! Land_Mask - Land mask array (only values with sym%YES will contribute)
! Elem_Start - starting x-index of array
! Elem_End - ending x-index of array
! Line_Start - starting y-index of array
! Line_End - ending y-index of array
! 
!
! Outputs: 
!  Loc_Max_Elem - Element index of the maximum value
!  Loc_Max_Line - Line index of the maximum value
!
! Dependencies: None
!
! Restrictions:  None
!
!====================================================================
SUBROUTINE Compute_NWC( Input_Array, &
                        Box_Size, &
                        Uni_Land_Mask_Flag, &
                        Bad_Mask, &
                        Land_Mask,  &
                        Elem_Start, &
                        Elem_End, &
                        Line_Start, &
                        Line_End, &
                        Loc_Max_Elem,  &
                        Loc_Max_Line)

  REAL(KIND=REAL4), DIMENSION(:,:), INTENT(IN):: Input_Array
  INTEGER(KIND=INT1), DIMENSION(:,:), INTENT(IN):: Bad_Mask
  INTEGER(KIND=INT1), DIMENSION(:,:), INTENT(IN):: Land_Mask
  REAL(KIND=REAL4):: Input_Array_NxN_Max
  INTEGER, INTENT(IN):: Uni_Land_Mask_Flag
  INTEGER, INTENT(IN):: Box_Size
  INTEGER, INTENT(IN):: Elem_Start
  INTEGER, INTENT(IN):: Elem_End
  INTEGER, INTENT(IN):: Line_Start
  INTEGER, INTENT(IN):: Line_End
  INTEGER, DIMENSION(:,:), INTENT(out):: Loc_Max_Elem
  INTEGER, DIMENSION(:,:), INTENT(out):: Loc_Max_Line
  INTEGER:: Elem_Idx 
  INTEGER:: Line_Idx
  INTEGER:: Elem_Idx_NxN_Right
  INTEGER:: Elem_Idx_NxN_Left
  INTEGER:: Line_Idx_NxN_Top
  INTEGER:: Line_Idx_NxN_Bottom
  INTEGER:: Elem_Idx_NxN
  INTEGER:: Line_Idx_NxN

  !--- initialize to missing
  Loc_Max_Elem = MISSING_VALUE_INT1
  Loc_Max_Line = MISSING_VALUE_INT1

  Line_Loop: DO Line_Idx = Line_Start, Line_End

  !--- set limits of NxN array in the j-direction
  Line_Idx_NxN_Top = max(Line_Start,Line_Idx-Box_Size)   !top index of local array
  Line_Idx_NxN_Bottom = min(Line_End,Line_Idx+Box_Size)     !bottom index of local array

  Element_Loop: DO Elem_Idx = Elem_Start, Elem_End

    !--- set limits of NxN array in the i-direction
    Elem_Idx_NxN_Right = max(Elem_Start,Elem_Idx - Box_Size)   !left index of local array
    Elem_Idx_NxN_Left = min(Elem_End,Elem_Idx + Box_Size)     !right index of local array

    !--- initialize
    Input_Array_NxN_Max = -1.0*huge(Input_Array_NxN_Max)

    !--- go through each element in NxN array
    Line_Loop_NxN: DO Line_Idx_NxN = Line_Idx_NxN_Top,Line_Idx_NxN_Bottom
      Element_Loop_NxN: DO Elem_Idx_NxN = Elem_Idx_NxN_Right,Elem_Idx_NxN_Left

        IF (Bad_Mask(Elem_Idx_NxN,Line_Idx_NxN) == sym%YES) THEN
          CYCLE
        ENDIF

        IF (Input_Array(Elem_Idx_NxN,Line_Idx_NxN) == MISSING_VALUE_REAL4) THEN
          CYCLE
        ENDIF
       
        IF ((Uni_Land_Mask_Flag == sym%YES) .and. &
             (Land_Mask(Elem_Idx,Line_Idx) /= Land_Mask(Elem_Idx_NxN,Line_Idx_Nxn))) THEN    
          CYCLE
        ENDIF

        IF (Input_Array(Elem_Idx_NxN,Line_Idx_Nxn) > Input_Array_NxN_Max) THEN
           Input_Array_NxN_Max = Input_Array(Elem_Idx_NxN,Line_Idx_Nxn)
           Loc_Max_Elem(Elem_Idx,Line_Idx) = Elem_Idx_NxN
           Loc_Max_Line(Elem_Idx,Line_Idx) = Line_Idx_NxN
        ENDIF

       END DO Element_Loop_NxN
    END DO Line_Loop_NxN

 END DO Element_Loop

END DO Line_Loop

END SUBROUTINE Compute_NWC

!--------------------------------------------------------------------------
! The subroutines and functions that contain the cloud mask tests go below here
!--------------------------------------------------------------------------
!====================================================================
! Subroutine Name: Set_Cmask_Thresholds
!
! Function:
!   Set Thresholds for tests that involve the 3.9 micron channel
!
! Description:
!   This Subroutine sets the thresholds for the Emiss_4 and ULST for a 
!   given satellite. This is needed because the 3.9 micron channel have
!   dIFferent spectral responce functions for various instruments.
!
! Calling Sequence:
!   CALL Set_Cmask_Thresholds(satellite_name)
!
! Inputs: Sat_name
!
! Outputs: None (thresholds in include are set)
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. In
!        this case, the thresholds for each satellite.
!
! Restrictions:  None
!
!====================================================================

 SUBROUTINE Set_Cmask_Thresholds(Sat_name, Algo_Name)
    CHARACTER(*), INTENT(IN) :: Sat_name
    CHARACTER(*), INTENT(inout) :: Algo_Name
    CHARACTER(len=1020) :: Algo_Name_Tmpy

    Algo_Name_Tmpy = TRIM(trim(Algo_Name)//"abi")
       
   ! Set thresholds for AVHRR  
   IF ((index(TRIM(Sat_name),'AVHRR')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_AVHRR
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_AVHRR
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_AVHRR
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_AVHRR
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_AVHRR
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_AVHRR
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_AVHRR

   ! Set thresholds for EOS TERRA MODIS
   ELSE IF ((index(TRIM(Sat_name),'Terra')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_MODIS
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_MODIS
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_MODIS
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_MODIS
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_MODIS
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_MODIS
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_MODIS

   ! Set thresholds for EOS AQUA MODIS
   ELSE IF ((index(TRIM(Sat_name),'Aqua')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_MODIS
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_MODIS
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_MODIS
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_MODIS
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_MODIS
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_MODIS
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_MODIS
         
   ! Set thresholds for MSG 
   ELSE IF ((index(TRIM(Sat_name),'Meteosat')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_SEVIRI
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_SEVIRI
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_SEVIRI
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_SEVIRI
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_SEVIRI
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_SEVIRI
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_SEVIRI
      Algo_Name_Tmpy = TRIM(trim(Algo_Name)//"seviri")

   ! Set thresholds for MTSAT
   ELSE IF ((index(TRIM(Sat_name),'MTSAT')) /= 0) THEN
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_MTSAT
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_MTSAT
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_MTSAT
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_MTSAT
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_MTSAT
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_MTSAT
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_MTSAT

   ! Set thresholds for GOES (I-P, ABI)
   ELSE
      ULST_Emiss_Chn7_DIFf_Thresh_Land = ULST_EMISS_CHN7_DIFF_THRESH_LAND_GOES
      ULST_Emiss_Chn7_DIFf_Thresh_Ocean = ULST_EMISS_CHN7_DIFF_THRESH_OCN_GOES
      ULST_Emiss_Chn7_DIFf_Thresh_Snow = ULST_EMISS_CHN7_DIFF_THRESH_SNOW_GOES
      EMISS4_Emiss_Chn7_Ocn_Thresh = EMISS4_EMISS_CHN7_OCN_THRESH_GOES
      EMISS4_Emiss_Chn7_Land_Thresh = EMISS4_EMISS_CHN7_LAND_THRESH_GOES
      EMISS4_Emiss_Chn7_Snow_Thresh = EMISS4_EMISS_CHN7_SNOW_THRESH_GOES
      EMISS4_Emiss_Chn7_Desert_Thresh = EMISS4_EMISS_CHN7_DESERT_THRESH_GOES
   ENDIF 
   
   
    Algo_Name = trim(Algo_Name_Tmpy)
   
 END SUBROUTINE Set_Cmask_Thresholds

!====================================================================
! Function Name: Compute_Clear_Sky_Scatter
!
! Function:
!   Computes the single scater and aerosol reflectance assuming that the
!   gas is mixed in with scattering
!
! Description:
   
!
! Calling Sequence:
!   Refl_Sing_Scat = Compute_Clear_Sky_Scatter(Aerosol_Optical_Depth_Chn2, &
!                                     Aerosol_Single_Scatter_Albedo_Chn2, &
!                                     Aerosol_Asymmetry_Parameter, &
!                                     Rayleigh_Optical_Depth_Chn2, &
!                                     Sat_Name, & 
!                                     TPW, &
!                                     Total_Ozone_Path_NWP, &
!                                     Scat_Zen, &
!                                     Cos_Sat_Zen, &
!                                     Cos_Sol_Zen, &
!                                     Sfc_Alb_View, &
!                                     Sfc_Alb_Sun, &
!                                     Transmission_Sing_Scat, &
!                                     Refl_Sing_Scat)
!
! Inputs:
!   Aerosol Optical depth at 0.64 micron
!   Aerosol Single Scatter Albedeo at 0.64 micron
!   Aerosol Asymmetry Parameter
!   Rayleigh Optical Depth at 0.64 micron
!   Satellite name
!   Total Precipitable Water
!   Total_Ozone_Path
!   scattering angle
!   cosine of the satellite zenith angle
!   cosine of the solar zenith angle
!   Clear sky Surface Albedo 
!   Clear Sky Surface Albedo
!   Single Scattering Transmission Coefficent
!
! Outputs: Single Scatter reflectance
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================


 SUBROUTINE Compute_Clear_Sky_Scatter(Aerosol_Optical_Depth_Chn2, &
                                     Aerosol_Single_Scatter_Albedo_Chn2, &
                                     Aerosol_Asymmetry_Parameter, &
                                     Rayleigh_Optical_Depth_Chn2, &
                                     Sat_Name, & 
                                     TPW, &
                                     Total_Ozone_Path_NWP, &
                                     Scat_Zen, &
                                     Cos_Sat_Zen, &
                                     Cos_Sol_Zen, &
                                     Sfc_Alb_View, &
                                     Sfc_Alb_Sun, &
                                     Transmission_Sing_Scat, &
                                     Refl_Sing_Scat)

   REAL, INTENT(IN):: Aerosol_Optical_Depth_Chn2
   REAL, INTENT(IN):: Aerosol_Single_Scatter_Albedo_Chn2
   REAL, INTENT(IN):: Aerosol_Asymmetry_Parameter
   REAL, INTENT(IN):: Rayleigh_Optical_Depth_Chn2
   CHARACTER (len=*), INTENT(IN):: Sat_Name 
   REAL, INTENT(IN):: TPW
   REAL, INTENT(IN):: Total_Ozone_Path_NWP
   REAL, INTENT(IN):: Scat_Zen
   REAL, INTENT(IN):: Cos_Sat_Zen
   REAL, INTENT(IN):: Cos_Sol_Zen
   REAL, INTENT(IN):: Sfc_Alb_View
   REAL, INTENT(IN):: Sfc_Alb_Sun
   REAL, INTENT(out):: Transmission_Sing_Scat
   REAL, INTENT(out):: Refl_Sing_Scat
   REAL:: Air_Mass_Factor
   REAL:: Aero_Phase_Funct
   REAL:: Ray_Phase_Funct
   REAL:: OD_Gas
   REAL:: OD_Total
   REAL:: OD_Scat_Total
   REAL:: OD_Iso_Total
   REAL:: Trans_Iso_Total_View
   REAL:: Trans_Iso_Total_Sun
   REAL:: OD_Iso_Scat_Total
   REAL:: Cos_Scat_Zen
   REAL:: Eff_Phase_Funct
   REAL:: Single_Scat_Alb
   REAL:: Refl_Sing_Scat_a
   REAL:: Refl_Sing_Scat_b
   REAL:: Refl_Sing_Scat_c
   REAL:: OD_ozone
   REAL:: OD_h2o
   REAL, DIMENSION(3):: OD_ozone_coef
   REAL, DIMENSION(3):: OD_h2o_coef

   !--- Set Default Gas Transmission Coefficients
   OD_ozone_coef = (/0.000566454,8.25224e-05,1.94007e-08/)
   OD_h2o_coef = (/  0.000044758, 0.00264790,-0.0000713698/)


   !--- Set Satellite SpecIFic Gas Transmission Coefficients
   IF ((index(TRIM(Sat_Name),'Meteosat')) /= 0) THEN
       OD_ozone_coef = (/0.000566454,8.25224e-05,1.94007e-08/)
       OD_h2o_coef = (/  0.000044758, 0.00264790,-0.0000713698/)
   ENDIF

   !--- compute cosine of scattering angle
   Cos_Scat_Zen = cos(Scat_Zen * dtor)

   !--- compute gaseous optical depth and transmission
   OD_h2o = OD_h2o_coef(1) + OD_h2o_coef(2)*TPW + OD_h2o_coef(3) * TPW**2
   
   OD_ozone = OD_ozone_coef(1) + OD_ozone_coef(2) * Total_Ozone_Path_NWP + &
              OD_ozone_coef(3) * Total_Ozone_Path_NWP**2

   OD_Gas = OD_h2o + OD_ozone

   !-- compute Rayleigh phase function
   Air_Mass_Factor = 1.0 / Cos_Sat_Zen + 1.0 / Cos_Sol_Zen
   
   Ray_Phase_Funct = 0.75 * (1.0 + Cos_Scat_Zen**2)

   !--- compute total transmission
   OD_Total = Aerosol_Optical_Depth_Chn2 + Rayleigh_Optical_Depth_Chn2 + OD_Gas
   
   Transmission_Sing_Scat = exp(-OD_Total * Air_Mass_Factor)

   OD_Iso_Total = (1.0 - Aerosol_Asymmetry_Parameter) * Aerosol_Optical_Depth_Chn2 + &
                   Rayleigh_Optical_Depth_Chn2 + OD_Gas
   
   Trans_Iso_Total_View = exp(-OD_Iso_Total / Cos_Sat_Zen)
   
   Trans_Iso_Total_Sun = exp(-OD_Iso_Total / Cos_Sol_Zen)

   !--- compute total scattering optical depth
   OD_Scat_Total = Aerosol_Single_Scatter_Albedo_Chn2 *&
                   Aerosol_Optical_Depth_Chn2 + Rayleigh_Optical_Depth_Chn2
   
   OD_Iso_Scat_Total = Aerosol_Single_Scatter_Albedo_Chn2 * (1.0 - Aerosol_Asymmetry_Parameter) * &
                    Aerosol_Optical_Depth_Chn2 + Rayleigh_Optical_Depth_Chn2

   !--- single scatter albedo
   Single_Scat_Alb = (Aerosol_Single_Scatter_Albedo_Chn2 * Aerosol_Optical_Depth_Chn2 + &
                     Rayleigh_Optical_Depth_Chn2) / ( OD_Total )

   !aerosol phase function (Henyey-Greenstein)
   Aero_Phase_Funct = (1.0 - Aerosol_Asymmetry_Parameter**2) / &
                    ( (1.0 + Aerosol_Asymmetry_Parameter**2 -  &
                       2.0 * Aerosol_Asymmetry_Parameter*Cos_Scat_Zen)**(1.5) )

   !--- compute effective phase function
   Eff_Phase_Funct = (Aerosol_Single_Scatter_Albedo_Chn2 * Aerosol_Optical_Depth_Chn2 * &
   Aero_Phase_Funct + Rayleigh_Optical_Depth_Chn2 * Ray_Phase_Funct) / (OD_Scat_Total)

   !--- compute single scatter reflectance (0-100%)
   Refl_Sing_Scat_a = Single_Scat_Alb * Eff_Phase_Funct / (4.0 * Air_Mass_Factor * &
                        Cos_Sat_Zen * Cos_Sol_Zen) * (1.0 - &
                        Transmission_Sing_Scat )

   Refl_Sing_Scat_b = (OD_Iso_Scat_Total / (2.0*Cos_Sol_Zen)) * &
                        Trans_Iso_Total_View * Sfc_Alb_View

   Refl_Sing_Scat_c = (OD_Iso_Scat_Total / (2.0 * Cos_Sat_Zen)) * &
                       Trans_Iso_Total_Sun * Sfc_Alb_Sun

   Refl_Sing_Scat = 100.0 * (Refl_Sing_Scat_a + Refl_Sing_Scat_b + Refl_Sing_Scat_c)

 END SUBROUTINE Compute_Clear_Sky_Scatter

!====================================================================
! Function Name: RUT_Routine
!
! Function:
!   Reflectance Uniformity Routine (RUT) Test
!
! Description:
!   Computes the Reflectance Uniformity Routine (RUT) test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RUT_Routine (&
!                              Is_Land, &
!                              Is_Snow, &
!                              Is_Coast, &
!                              Refl_Chn2_Clear(Elem_Idx,Line_Idx), &
!                              Refl_Chn2_Stddev_3x3(Elem_Idx,Line_Idx), &
!                              Sol_Zen)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel coast (YES/NO)
!   0.64 miron clear sky reflectance
!   Standard Deviation of the 0.64 micron reflectance over a 3x3 box
!   Solar zenith angle
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test (Test_Result)
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RUT_Routine(Is_Land, &
                     Is_Snow, &
                     Is_Coast, &
                     Refl_Chn2_Clear, &
                     Refl_Chn2_Stddev_3x3, &
                     Sol_Zen) &
                     RESULT(Test_Result)

       INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
       INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
       INTEGER(KIND=INT4), INTENT(IN) :: Is_Coast
       REAL(KIND=REAL4), INTENT(IN) :: Refl_Chn2_Clear
       REAL(KIND=REAL4), INTENT(IN) :: Refl_Chn2_Stddev_3x3
       REAL(KIND=REAL4), INTENT(IN) :: Sol_Zen
       INTEGER(KIND=int1) :: Test_Result

       REAL(KIND=REAL4) :: Test_Threshold

       Test_Result = sym%NO

       IF (Is_Snow == sym%NO .AND. Is_Coast == sym%NO) THEN

          IF (Sol_Zen < RUT_SOL_ZEN_THRESH) THEN

              !--- compute threshold
              IF (Is_Land == sym%YES) THEN
                 Test_Threshold = MAX(0.5,Refl_Chn2_Clear * REFL_CHN2_CLR_UNI_THRESH_LAND)
              ELSE
                 Test_Threshold = REFL_CHN2_CLR_UNI_THRESH_OCEAN
              ENDIF

              !--- apply test
              IF (Refl_Chn2_Stddev_3x3 > Test_Threshold) THEN
                     Test_Result = sym%YES
              ENDIF

           ENDIF

       ENDIF

       RETURN

 END FUNCTION RUT_Routine

!====================================================================
! Function Name: TUT_Routine
!
! Function:
!   Thermal Uniformity Routine (TUT) Test
!
! Description:
!   Computes the Thermal Uniformity Routine (TUT) test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = TUT_Routine (&
!                              Is_Land, &
!                              Is_Coast, &
!                              BT_Chn14_Stddev_3x3, &
!                              Sfc_Hgt_Stddev_3x3)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   Standard Deviation of the 11 micron BT over a 3x3 box
!   Standard Deviation of the surface height over a 3x3 box
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION TUT_Routine(Is_Land, &
                      Is_Coast, &
                      BT_Chn14_Stddev_3x3, &
                      Sfc_Hgt_Stddev_3x3)  &
                      RESULT (Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land 
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Coast
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Stddev_3x3
   REAL(KIND=REAL4), INTENT(IN) :: Sfc_Hgt_Stddev_3x3
   INTEGER(KIND=int1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold

   Test_Result = sym%NO

   IF (Is_Coast == sym%NO) THEN

         !
         !7K/km is the adiabatic lapse rate
         !
         Test_Threshold = 3.0 * 7.0*Sfc_Hgt_Stddev_3X3/1000.0  

         IF (Is_Land == sym%YES) THEN

      IF (BT_Chn14_Stddev_3x3 > BT_CHN14_CLR_UNI_THRESH_LAND  &
          +Test_Threshold) THEN

         Test_Result = sym%YES

      ENDIF

         ELSE

      IF (BT_Chn14_Stddev_3x3 > BT_CHN14_CLR_UNI_THRESH_OCN  &
           +Test_Threshold) THEN

         Test_Result = sym%YES

      ENDIF

         ENDIF

   ENDIF

   RETURN

 END FUNCTION TUT_Routine

!====================================================================
! Function Name: RTCT_Routine
!
! Function:
!   Relative Thermal Contrast Test
!
! Description:
!   Computes the Relative Thermal Contrast Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RTCT_Routine(&
!                              Is_Land,  &
!                              Is_Coast,  &
!                              Is_Snow, &
!                              Is_Cold_Surface, &
!                              BT_Chn14, &
!                              BT_Chn14_Min_3x3(Elem_Idx,Line_Idx), &
!                              BT_Chn14_Max_3x3(Elem_Idx,Line_Idx), &
!                              Sfc_Hgt_Stddev_3x3(Elem_Idx,Line_Idx))
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel cold surface (YES/NO)
!   11 micron BT 
!   Minimum of the 11 micron BT in a 3x3 box
!   Maximum of the 11 micron BT in a 3x3 box
!   Standard Deviation of the surface height over a 3x3 box
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RTCT_Routine(Is_Land, &
                      Is_Coast,  &
                      Is_Snow, &
                      Is_Cold_Surface, &
                      BT_Chn14, &
                      BT_Chn14_Min_3x3, &
                      BT_Chn14_Max_3x3, &
                      Sfc_Hgt_Stddev_3x3) &
                      RESULT (Test_Result)

   INTEGER (KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER (KIND=INT4), INTENT(IN) :: Is_Coast
   INTEGER (KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER (KIND=INT4), INTENT(IN) :: Is_Cold_Surface
   REAL (KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL (KIND=REAL4), INTENT(IN) :: BT_Chn14_Min_3x3
   REAL (KIND=REAL4), INTENT(IN) :: BT_Chn14_Max_3x3
   REAL (KIND=REAL4), INTENT(IN) :: Sfc_Hgt_Stddev_3x3
   INTEGER (KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold

   Test_Result = sym%NO

   IF (Is_Cold_Surface == sym%NO) THEN

       IF (Is_Land == sym%NO) THEN

            Test_Threshold = RTCT_OCN_THRESH

       ELSE

            Test_Threshold = RTCT_LAND_THRESH

       ENDIF

       !
       !7K/km is the adiabatic lapse rate
       !
       Test_Threshold = Test_Threshold + 3.0 * 7.0*Sfc_Hgt_Stddev_3X3/1000.0

       IF ((Is_Coast == sym%NO) .AND.  &
           (Is_Snow == sym%NO) .AND.  &
           (BT_Chn14_Min_3x3 <= 300.0) .AND.  &
           (BT_Chn14_Max_3x3 - BT_Chn14 > Test_Threshold)) THEN

           Test_Result = sym%YES

       ENDIF

   ENDIF

 END FUNCTION RTCT_Routine

!====================================================================
! Function Name: ETROP_Routine
!
! Function:
!   11 micron Tropospheric Emissivity Test
!
! Description:
!   Computes the 11 micron Tropospheric Emissivity Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  ETROP_Routine( &
!                               Is_Snow, &
!                               Is_Land, &
!                               Is_Coast, &
!                               Is_Desert, &
!                               Is_Cold_Surface, &
!                               Land_Type, &
!                               BT_Chn14, &
!                               BT_Chn14_Clr, &
!                               Emiss_Tropo_Chn14(Elem_Idx,Line_Idx), &
!                               Emiss_Tropo_Chn14_LRC, &
!                               BT_Chn14_Stddev_3x3(Elem_Idx,Line_Idx))
!
! Inputs:
!   Is pixel snow (YES/NO)
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   Is pixel desert (YES/NO)
!   Is pixel cold surface (YES/NO)
!   Is pixel land type
!   11 micron BT 
!   Clear sky 11 micron BT 
!   11 micron tropopause emissivity
!   11 micron tropopause emissivity at the Local radiative center
!   Standard Deviation of the 11 micron BT over a 3x3 box
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION ETROP_Routine(Is_Snow, &
                        Is_Land, &
                        Is_Coast, &
                        Is_Desert, &
                        Is_Cold_Surface, &
                        Land_Type, &
                        BT_Chn14, &
                        BT_Chn14_Clr, &
                        Emiss_Tropo_Chn14, &
                        Emiss_Tropo_Chn14_LRC, &
                        BT_Chn14_Stddev_3x3) &
                        RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Coast
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Desert
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Cold_Surface
   INTEGER(KIND=INT4), INTENT(IN) :: Land_Type
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Clr
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Tropo_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Tropo_Chn14_LRC
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Stddev_3x3
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Value
   REAL(KIND=REAL4) :: Test_Threshold
   REAL(KIND=REAL4) :: Test_LRC_Threshold

   Test_Result = sym%NO

   IF ((BT_Chn14 > 170.0) .AND.  &
       (BT_Chn14 < 315.0) .AND.  &
       (BT_Chn14_Clr > 240.0)) THEN

       !--- assign threshold
       Test_Threshold = EMISS_CHN14_TROPO_LAND_THRESH
       Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_LAND_THRESH
       
       IF (Is_Land == sym%NO .AND. Is_Coast == sym%NO) THEN 
         Test_Threshold = EMISS_CHN14_TROPO_OCN_THRESH
         Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_OCN_THRESH
       ENDIF
       IF (Is_Snow == sym%YES) THEN 
         Test_Threshold = EMISS_CHN14_TROPO_SNOW_THRESH
         Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_SNOW_THRESH
       ENDIF
       IF (Is_Desert == sym%YES) THEN 
         Test_Threshold = EMISS_CHN14_TROPO_DESERT_THRESH
         Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_DESERT_THRESH
       ENDIF
       IF (Is_Cold_Surface == sym%YES) THEN 
         Test_Threshold = EMISS_CHN14_TROPO_COLD_SURFACE_THRESH
         Test_LRC_Threshold = EMISS_CHN14_TROPO_LRC_COLD_SURFACE_THRESH
       ENDIF

       !--- select value to test
       Test_Value = Emiss_Tropo_Chn14

       !--- apply test
       IF (Test_Value > Test_Threshold) THEN
         Test_Result = sym%YES
       ENDIF

       !--- apply LRC portion of test
       IF (Emiss_Tropo_Chn14_LRC /= Missing_Value_REAL4) THEN

         !--- select value to test
         Test_Value = Emiss_Tropo_Chn14_LRC
        
         !--- apply test
         IF (Test_Value > Test_LRC_Threshold) THEN
          Test_Result = sym%YES
         ENDIF

       ENDIF



      !--- apply a tigher threshold for highly nonunIForm pixels
      Test_Threshold = 0.20
      IF ((Is_Coast == sym%NO) .AND. (BT_Chn14_Stddev_3x3 > 0.5)) THEN
         IF (Test_Value > Test_Threshold) THEN
             Test_Result = sym%YES
         ENDIF
      ENDIF 


      !--- perform a restoral near land where sst field is often errorneous
      !--- select value to test
      Test_Value = Emiss_Tropo_Chn14

      IF (Test_Result == sym%YES) THEN
         IF ((Test_Value < 0.20) .AND.  &
            (BT_Chn14_Stddev_3x3 < 1.0) .AND. &
            (Land_Type /= sym%LAND) .AND.  &
            (Land_Type /= sym%DEEP_OCEAN)) THEN

            Test_Result = sym%NO

         END IF 
      END IF 

   ENDIF

   RETURN

 END FUNCTION ETROP_Routine

!====================================================================
! Function Name: PFMFT_Routine
!
! Function:
!   Positive Split Window DIFferences Test (PFMFT) 
!
! Description:
!   Computes the Positive Split Window DIFferences Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = PFMFT_Routine( &
!                      Is_Land, &
!                      Is_Snow, &
!                      Is_Cold_Surface, &
!                      BT_Chn14, &
!                      BT_Chn15, &
!                      BT_Chn14_Clr, &
!                      BT_Chn15_Clr, &
!                      BT_Chn14_Stddev_3x3)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel cold surface (YES/NO)
!   11 micron BT 
!   12 micron BT 
!   Clear sky 11 micron BT 
!   Clear sky 12 micron BT 
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION PFMFT_Routine( &
                      Is_Land, &
                      Is_Snow, &
                      Is_Cold_Surface, &
                      BT_Chn14, &
                      BT_Chn15, &
                      BT_Chn14_Clr, &
                      BT_Chn15_Clr, &
                      BT_Chn14_Stddev_3x3) &
                      RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Cold_Surface
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Clr
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15_Clr
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Stddev_3x3
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold
   REAL(KIND=REAL4) :: Test_Value

   Test_Result = sym%NO

   IF((BT_Chn14 < BT_CHN14_MAX_FMFT_THRESH) .AND. &
      (BT_Chn14_Clr - BT_Chn15_Clr > BTDIFF_CHN14_CHN15_MIN_FMFT_THRESH) .AND. &
      (BT_Chn14_Stddev_3x3 > PFMFT_BT_CHN14_STDDEV_3x3_THRESH)) THEN

      IF ((BT_Chn14 > 270.0) .AND. (BT_Chn14_Clr > 270.0)) then
        Test_Value = (BT_Chn14 - BT_Chn15) -  &
                     (BT_Chn14_Clr - BT_Chn15_Clr) *(BT_Chn14 - 260.0) / &
                     (BT_Chn14_Clr - 260.0)
      ELSE
        Test_Value = (BT_Chn14 - BT_Chn15)
      ENDIF
 
 
      !-- set appropriate threshold
      IF (Is_Land == sym%NO) THEN
            Test_Threshold = PFMFT_OCEAN_THRESH
      ELSE
            Test_Threshold = PFMFT_LAND_THRESH
      ENDIF
      IF (Is_Snow == sym%YES) THEN
            Test_Threshold = PFMFT_SNOW_THRESH
      ENDIF
      IF (Is_Cold_Surface == sym%YES) THEN
            Test_Threshold = PFMFT_COLD_SURFACE_THRESH
      ENDIF

      !--- Apply Test
      IF (Test_Value > Test_Threshold) THEN
          Test_Result = sym%YES
      ENDIF
 
   ENDIF
 
   RETURN

 END FUNCTION PFMFT_Routine
 
!====================================================================
! Function Name: NFMFT_Routine
!
! Function:
!   Test for Negative Split Window DIFferences (NFMFT)
!
! Description:
!   Computes the Negative Split Window DIFferences Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = NFMFT_Routine( &
!                      Is_Land, &
!                      Is_Snow, &
!                      BT_Chn14, &
!                      BT_Chn15, &
!                      BT_Chn14_Clr, &
!                      BT_Chn15_Clr)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel desert (YES/NO)
!   11 micron BT 
!   12 micron BT 
!   Clear sky 11 micron BT 
!   Clear sky 12 micron BT 
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION NFMFT_Routine( &
                      Is_Land, &
                      Is_Snow, &
                      Is_Desert, &
                      BT_Chn14, &
                      BT_Chn15, &
                      BT_Chn14_Clr, &
                      BT_Chn15_Clr) &
                      RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Desert
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Clr
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15_Clr
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold
   REAL(KIND=REAL4) :: Test_Value

   Test_Result = sym%NO

   !--- skip this test for elevated values of BT_Chn14-BT_Chn15
   IF (BT_Chn14 - BT_Chn15 > NFMFT_BTD_CHN14_CHN15_MAX_THRESH) THEN
      RETURN
   ENDIF

   !-- set appropriate threshold
   Test_Threshold = NFMFT_LAND_THRESH
   IF (Is_Land == sym%NO) THEN
         Test_Threshold = NFMFT_OCEAN_THRESH
   ENDIF
   IF (Is_Desert == sym%YES) THEN
        Test_Threshold = NFMFT_DESERT_THRESH
   ENDIF
   IF (Is_Snow == sym%YES) THEN
        Test_Threshold = NFMFT_SNOW_THRESH
   ENDIF

   !--- Compute Value to Test
   Test_Value = (BT_Chn14 - BT_Chn15) -  &
                (BT_Chn14_Clr - BT_Chn15_Clr)

   !--- Apply Test
   IF (Test_Value < Test_Threshold) THEN
       Test_Result = sym%YES
   ENDIF

 END FUNCTION NFMFT_Routine

!====================================================================
! Function Name: RFMFT_Routine
!
! Function:
!   Relative Split Window DIFferences Test (RFMFT) 
!
! Description:
!   Computes the Relative Split Window DIFferences Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RFMFT_Routine( &
!                      Is_Land, &
!                      Is_Coast, &
!                      BT_Chn14, &
!                      BT_Chn15, &
!                      BTD_Chn14_Chn15_NWC)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   11 micron BT 
!   12 micron BT 
!   11 - 12 micron BT DIFf at the nearest warm center
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RFMFT_Routine( &
                      Is_Land, &
                      Is_Coast, &
                      BT_Chn14, &
                      BT_Chn15, &
                      BTD_Chn14_Chn15_NWC) &
                      RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Coast
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn15
   REAL(KIND=REAL4), INTENT(IN) :: BTD_Chn14_Chn15_NWC
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Value
   REAL(KIND=REAL4) :: Test_Threshold

   !--- initialize
   Test_Result = sym%NO

   !--- do not apply over hot land due to NWP surface temperature biases
   IF ((Is_Land == sym%YES) .AND. (BT_Chn14 > RFMFT_BT_CHN14_MAX_THRESH)) THEN
        Test_Result = sym%NO
        RETURN
   ENDIF

   IF (BT_Chn14 - BT_Chn15 < RFMFT_BTDIFF_CHN14_CHN15_MIN_THRESH) THEN
        Test_Result = sym%NO
        RETURN
   ENDIF

   IF (Is_Coast == sym%YES) THEN
        Test_Result = sym%NO
        RETURN
   ENDIF

   !--- test for high departures - cirrus

   !--- pick thresholds
   Test_Threshold = RFMFT_HI_LAND_THRESH
   IF (Is_Land == sym%NO) THEN
        Test_Threshold = RFMFT_HI_OCEAN_THRESH
   ENDIF

   !--- apply test
   Test_Value = abs((BT_Chn14 - BT_Chn15) - BTD_Chn14_Chn15_NWC)

   IF (Test_Value > Test_Threshold) THEN
        Test_Result = sym%YES
   ENDIF

 END FUNCTION RFMFT_Routine

!====================================================================
! Function Name: RFMFT_Routine
!
! Function:
!   Temporal IR Test (RFMFT) 
!
! Description:
!   Computes the Temporal IR Test (YES/NO), as described in the ABI cloud
!   mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = TEMPIR_Routine( &
!                          BT_Chn14, &
!                          BT_Chn14_Clr, &
!                          BT_Chn14_15min, &
!                          BT_Chn14_15min_Clr)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel coast (YES/NO)
!   11 micron BT 
!   Clear sky 11 micron BT 
!   11 micron BT from 15 minutes ago for same pixel
!   Clear sky 11 micron BT from 15 minutes ago for same pixel
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION  TEMPIR_Routine( &
                          BT_Chn14, &
                          BT_Chn14_Clr, &
                          BT_Chn14_15min, &
                          BT_Chn14_15min_Clr) &
                          RESULT(Test_Result)

   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_Clr
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_15min
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14_15min_Clr
   REAL(KIND=REAL4) :: Test_Threshold
   INTEGER(KIND=INT1) :: Test_Result

   !--- initialize output
   Test_Result = sym%NO

   IF ((BT_Chn14_15min /= Missing_Value_REAL4) .AND. &
       (BT_Chn14_15min_Clr /= Missing_Value_REAL4) .AND. &
       (BT_Chn14_15min < 330.0) .AND. &
       (BT_Chn14_15min_Clr < 330.0)) THEN

      !--- make threshold
      Test_Threshold = (BT_Chn14_15min_Clr - BT_Chn14_Clr) + &
                        TEMPIR_BT_CHN14_15MIN_TEMPORAL_OFFSET

      !--- apply test
      IF ((BT_Chn14_15min - BT_Chn14) > Test_Threshold)  THEN

        Test_Result = sym%YES
 
      ENDIF

   ENDIF

   RETURN

 END FUNCTION TEMPIR_Routine

!====================================================================
! Function Name: RGCT_Routine
!
! Function:
!   Reflectance Gross Constrast Test (RGCT)
!
! Description:
!   Computes the Reflectance Gross Constrast Test (RGCT) (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RGCT_Routine( &
!                       Is_Snow, &
!                       Is_Glint, &
!                       Sol_Zen, &
!                       Refl_Chn2, &
!                       RGCT_Threshold)
!
! Inputs:
!   Is pixel snow (YES/NO)
!   Is pixel glint (YES/NO)
!   Solar zenith angle
!   Air mass factor
!   0.64 miron reflectance
!   RGCT threshold for that pixel
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RGCT_Routine( &
                       Is_Snow, &
                       Is_Glint, &
                       Sol_Zen, &
                       Refl_Chn2, &
                       RGCT_Threshold) &
                       RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Glint
   REAL(KIND=REAL4), INTENT(IN) :: Sol_Zen
   REAL(KIND=REAL4), INTENT(IN) :: Refl_Chn2
   REAL(KIND=REAL4), INTENT(IN) :: RGCT_Threshold
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold


   !--- initialize
   Test_Result = sym%NO


   !--- apply test
   IF (Is_Snow == sym%NO .AND. Is_Glint == sym%NO) THEN

        IF (Sol_Zen < RGCT_Sol_Zen_Thresh) THEN 

            Test_Threshold = RGCT_Threshold

            IF (Refl_Chn2 > Test_Threshold) THEN

                Test_Result =  sym%YES
  
            ENDIF

        ENDIF

   ENDIF

   RETURN

 END FUNCTION RGCT_Routine

!====================================================================
! Function Name: RVCT_Routine
!
! Function:
!   Relative Visible Contrast Test (RVCT)
!
! Description:
!   Computes the Relative Visible Contrast Test (RGCT) (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
!
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = RVCT_Routine( &
!                        Is_Coast, &
!                        Is_Snow, &
!                        Is_Land, &
!                        Refl_Chn2, &
!                        Scat_Zen, &
!                        Sol_Zen, &
!                        Refl_Chn2_Clear_Stddev_3x3, &
!                        Refl_Chn2_Min_3x3)
!
! Inputs:
!   Is pixel coast (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel land(YES/NO)
!   Is pixel glint (YES/NO)
!   0.64 miron reflectance
!   Scattering zenith angle
!   Solar zenith angle
!   Standard Deviation of the clear sky 0.64 micron reflectance over a 3x3 box
!   Minima of the 0.64 micron reflectance over a 3x3 box
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION RVCT_Routine( &
                        Is_Coast, &
                        Is_Snow, &
                        Is_Land, &
                        Refl_Chn2, &
                        Scat_Zen, &
                        Sol_Zen, &
                        Refl_Chn2_Clear_Stddev_3x3, &
                        Refl_Chn2_Min_3x3) &
                        RESULT(Test_Result)

   INTEGER(KIND=INT4),INTENT(IN) :: Is_Coast
   INTEGER(KIND=INT4),INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4),INTENT(IN) :: Is_Land
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn2
   REAL(KIND=REAL4),INTENT(IN) :: Scat_Zen
   REAL(KIND=REAL4),INTENT(IN) :: Sol_Zen
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn2_Clear_Stddev_3x3
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn2_Min_3x3
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold

   !--- initialize
   Test_Result = sym%NO

   !--- derive threshold
   Test_Threshold = 999.0
   IF ( (Is_Coast == sym%NO) .AND. (Is_Snow == sym%NO) .AND.  &
        (Scat_Zen > RVCT_SCAT_ZEN_THRESH) .AND. &
        (Sol_Zen < RVCT_SOL_ZEN_THRESH)) THEN

      Test_Threshold = 5.0 

      IF (Is_Land == sym%YES) THEN
           Test_Threshold = Test_Threshold + 5.0 + 4.0*Refl_Chn2_Clear_Stddev_3x3 
      ENDIF

      IF (Is_Land == sym%NO .OR. Refl_Chn2_Clear_Stddev_3x3 <= 0.0) THEN
           Test_Threshold = 10.0
      ENDIF

      !--- apply test
      IF (Refl_Chn2 - Refl_Chn2_Min_3x3 > Test_Threshold) THEN
           Test_Result = sym%YES
      ENDIF

   ENDIF

   RETURN

 END FUNCTION RVCT_Routine

!====================================================================
! Function Name: NIRREF_Chn5_Routine
!
! Function:
!   1.6 micron (Ch4) Near-IR Reflectance Test
!
! Description:
!   Computes the 1.6 micron (Ch4) Near-IR Reflectance Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  NIRREF_Chn5_Routine( &
!                         Is_Coast, &
!                         Is_Snow, &
!                         Sol_Zen, &
!                         Sfc_Hgt, &
!                         NDSI, &
!                         Refl_Chn5)
!
! Inputs:
!   Is pixel coast (YES/NO)
!   Is pixel snow (YES/NO)
!   Solar zenith angle
!   Surface height
!   Normalized Difference Snow Index
!   1.5 micron reflection
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION NIRREF_Chn5_Routine(&
                         Is_Coast, &
                         Is_Snow, &
                         Sol_Zen, &
                         Sfc_Hgt, &
                         NDSI, &
                         Refl_Chn5) &
                         RESULT(Test_Result)

   INTEGER(KIND=INT4),INTENT(IN) :: Is_Coast
   INTEGER(KIND=INT4),INTENT(IN) :: Is_Snow
   REAL(KIND=REAL4),INTENT(IN) :: Sol_Zen
   REAL(KIND=REAL4),INTENT(IN) :: Sfc_Hgt
   REAL(KIND=REAL4),INTENT(IN) :: NDSI
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn5
   INTEGER(KIND=INT1) :: Test_Result

    !--- initialize
    Test_Result = sym%NO

    !--- apply test
    IF ( (Is_Coast == sym%NO) .AND.  &
         (Is_Snow == sym%YES) .AND.  &
         (Sol_Zen < NIRREF_CHN5_SOL_ZEN_THRESH) .AND.  &
         (NDSI < NIRREF_NDSI_THRESH_SNOW) .AND.  &
         (Sfc_Hgt < NIRREF_SFC_HGT_LIMIT)) THEN

         IF (Refl_Chn5 > NIRREF_CHN5_REFL_THRESH_SNOW) THEN

            Test_Result = sym%YES

         ENDIF

   ENDIF

   RETURN

 END FUNCTION NIRREF_Chn5_Routine
 
!====================================================================
! Function Name: NIRREF_Chn7_Routine
!
! Function:
!   3.9 micron (Ch7) Near-IR Reflectance Test
!
! Description:
!   Computes the 3.9 micron (Ch7) Near-IR Reflectance Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  NIRREF_Chn7_Routine( &
!                         Is_Snow, &
!                         Sol_Zen, &
!                         Sfc_Hgt, &
!                         Refl_Chn7)
!
! Inputs:
!   Is pixel snow (YES/NO)
!   Solar zenith angle
!   Surface height
!   3.9 micron reflectance
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION NIRREF_Chn7_Routine(&
                         Is_Snow, &
                         Sol_Zen, &
                         Sfc_Hgt, &
                         Refl_Chn7) &
                         RESULT(Test_Result)

   INTEGER(KIND=INT4),INTENT(IN) :: Is_Snow
   REAL(KIND=REAL4),INTENT(IN) :: Sol_Zen
   REAL(KIND=REAL4),INTENT(IN) :: Sfc_Hgt
   REAL(KIND=REAL4),INTENT(IN) :: Refl_Chn7
   INTEGER(KIND=INT1) :: Test_Result

    !--- initialize
    Test_Result = sym%NO

    !--- apply test
    IF ( (Is_Snow == sym%YES) .AND.  &
         (Sol_Zen < NIRREF_CHN7_SOL_ZEN_THRESH) .AND.  &
         (Sfc_Hgt < NIRREF_SFC_HGT_LIMIT)) THEN

         IF (Refl_Chn7 > NIRREF_CHN7_REFL_THRESH_SNOW) THEN

            Test_Result = sym%YES

         ENDIF

   ENDIF

   RETURN

 END FUNCTION NIRREF_Chn7_Routine

!====================================================================
! Function Name: CIRREF_Routine
!
! Function:
!   Cirrus Reflectance Test (1.38 micron), CIRREF
!
! Description:
!   Computes the Cirrus Reflectance Test (1.38 micron) (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  CIRREF_Routine( &
!                         Is_Snow, &
!                         Refl_Chn4, &
!                         Sfc_Hgt_Max_3x3, &
!                         CIRREF_Sfc_Hgt_Limit, &
!                         Sol_Zen, &
!                         Refl_Test_Sol_Zen_Thresh)!
! Inputs:
!   Is pixel snow (YES/NO)
!   1.38 micron reflectance
!   Minima of the surface height over a 3x3 box
!   CIRREF surface height limit threshold
!   Solar zenith angle
!   reflectance test solar zenith threshold
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION  CIRREF_Routine( &
                         Is_Snow, &
                         Refl_Chn4, &
                         Sfc_Hgt_Max_3x3, &
                         Sol_Zen) &
                         RESULT(Test_Result)

   INTEGER(KIND=INT4) :: Is_Snow
   REAL(KIND=REAL4), INTENT(IN) :: Refl_Chn4
   REAL(KIND=REAL4), INTENT(IN) :: Sfc_Hgt_Max_3x3
   REAL(KIND=REAL4), INTENT(IN) :: Sol_Zen
   INTEGER(KIND=INT1) :: Test_Result

   !--- initialize
   Test_Result = sym%NO

   !--- apply test
   IF ( (Is_Snow == sym%NO) .AND.             &
        (Sfc_Hgt_Max_3X3 < CIRREF_SFC_HGT_LIMIT) .AND. &
        (Sol_Zen < CIRREF_SOL_ZEN_THRESH)) THEN

      IF (Refl_Chn4 > CIRREF_THRESH) THEN

         Test_Result = sym%YES

      ENDIF

   ENDIF

   RETURN

 END FUNCTION CIRREF_Routine

!====================================================================
! Function Name: EMISS4_Routine
!
! Function:
!   4 micron (Ch7) Emissivity (EMISS4) Cloud Test 
!
! Description:
!   Computes the 4 micron (Ch7) Emissivity (EMISS4) Cloud Test (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  EMISS4_Routine(Is_Glint,  &
!                         Is_Land, &
!                         Is_Snow, &
!                         Is_Desert, &
!                         BT_Chn14, &
!                         Emiss_Chn7,  &
!                         Emiss_Chn7_Clr
! Inputs:
!   Is pixel glint (YES/NO)
!   Is pixel land (YES/NO)
!   Is pixel snow (YES/NO)
!   Is pixel Desert (YES/NO)
!   11 micron BT
!   3.9 micron emissvity
!   Clear sky 3.9 micron emissivity
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================
 FUNCTION EMISS4_Routine(Is_Glint,  &
                         Is_Land, &
                         Is_Snow, &
                         Is_Desert, &
                         BT_Chn14, &
                         Emiss_Chn7,  &
                         Emiss_Chn7_Clr, &
                         Sfc_Emiss_Chn7_RTM) &
                         RESULT (Test_Result)
  INTEGER(KIND=INT4), INTENT(IN):: Is_Glint
  INTEGER(KIND=INT4), INTENT(IN):: Is_Land
  INTEGER(KIND=INT4), INTENT(IN):: Is_Snow
  INTEGER(KIND=INT4), INTENT(IN):: Is_Desert
  REAL(KIND=REAL4), INTENT(IN):: BT_Chn14
  REAL(KIND=REAL4), INTENT(IN):: Emiss_Chn7
  REAL(KIND=REAL4), INTENT(IN):: Emiss_Chn7_Clr
  REAL(KIND=REAL4), INTENT(IN):: Sfc_Emiss_Chn7_RTM
  INTEGER(KIND=INT1):: Test_Result
  REAL(KIND=REAL4):: Test_Value
  REAL(KIND=REAL4):: Test_Threshold

  Test_Result = sym%NO

   !
   !---avoid glint
   !
   IF (Is_Glint == sym%YES) THEN   
      RETURN
   ENDIF

   !
   !---avoid warm pixels
   !
   IF (BT_Chn14 > EMISS4_BT_CHN14_MAX_THRESH) THEN   
      RETURN
   ENDIF

   !--- Determine threshold based on surface condition
   Test_Threshold = EMISS4_Emiss_Chn7_Ocn_Thresh
   IF (Is_Land == sym%Yes) THEN
       Test_Threshold = EMISS4_Emiss_Chn7_Land_Thresh
   ENDIF
   IF (Is_Desert == sym%Yes) THEN
       Test_Threshold = EMISS4_Emiss_Chn7_Desert_Thresh
   ENDIF
   IF (Is_Snow == sym%Yes) THEN
       Test_Threshold = EMISS4_Emiss_Chn7_Snow_Thresh
   ENDIF

   !--- Augment in presence of very low emissive surfaces
   IF (Sfc_Emiss_Chn7_RTM < EMISS4_SFC_EMISS_CHN7_THRESH) then
       Test_Threshold = Test_Threshold + 0.5
   ENDIF

   !--- Compute value to test
   Test_Value = (Emiss_Chn7 - Emiss_Chn7_Clr) / Emiss_Chn7_Clr

   !--- apply test
   IF (Test_Value > Test_Threshold) THEN

       Test_Result = sym%YES

   ENDIF


   RETURN

END FUNCTION EMISS4_Routine


!====================================================================
! Function Name: ULST_Routine
!
! Function:
!   UnIForm Low Stratus Test - ULST
!
! Description:
!   Computes the UnIForm Low Stratus Test (ULST) (YES/NO), as described in
!   the ABI cloud mask ATBD
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) =  ULST_Routine(&
!                       Is_Land, &
!                       Is_Day, &
!                       Is_Snow, &
!                       BT_Chn14,  &
!                       BT_Chn14_Stddev_3x3,  &
!                       Emiss_Chn7,  &
!                       Emiss_Chn7_Clr, &
!                       Emiss_Tropo_Chn14)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Is pixel day (YES/NO)
!   Is pixel snow (YES/NO)
!   11 micron BT
!   Standard Deviation of the 11 micron BT over a 3x3 box
!   3.9 micron emissvity
!   Clear sky 3.9 micron emissivity
!   11 micron tropopause emissvity
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
!====================================================================

 FUNCTION ULST_Routine(&
                       Is_Land, &
                       Is_Day, &
                       Is_Snow, &
                       Is_Cold_Surface, &
                       BT_Chn14,  &
                       Emiss_Chn7,  &
                       Emiss_Chn7_Clr, &
                       Emiss_Chn7_NWC, &
                       Sfc_Emiss_Chn7) &
                       RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN) :: Is_Land
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Day
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Snow
   INTEGER(KIND=INT4), INTENT(IN) :: Is_Cold_Surface
   REAL(KIND=REAL4), INTENT(IN) :: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Chn7
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Chn7_Clr
   REAL(KIND=REAL4), INTENT(IN) :: Emiss_Chn7_NWC
   REAL(KIND=REAL4), INTENT(IN) :: Sfc_Emiss_Chn7
   INTEGER(KIND=INT1) :: Test_Result
   REAL(KIND=REAL4) :: Test_Threshold 

   !--- initialize
   Test_Result = sym%NO

   !--- exclude day time data
   IF (Is_Day == sym%YES) THEN
        RETURN
   ENDIF

   !--- exclude day time data
   IF (Is_Cold_Surface == sym%YES) THEN
        RETURN
   ENDIF

   !--- exclude very warm pixels
   IF (BT_Chn14 > ULST_BT_CHN14_MAX_THRESH) THEN
        RETURN
   ENDIF

   !--- exclude pixels with large Emiss_Chn7
   IF (Emiss_Chn7 > ULST_EMISS_CHN7_MAX_THRESH) THEN
        RETURN
   ENDIF

   !--- exclude deserts
   IF (Sfc_Emiss_Chn7 < ULST_EMISS_CHN7_SFC_THRESH) THEN
        RETURN
   ENDIF

   !--- exclude values whose clear values are suspect
   IF ((Emiss_Chn7_Clr < ULST_EMISS_CHN7_CLR_MIN_THRESH) .OR. &
       (Emiss_Chn7_Clr > ULST_EMISS_CHN7_CLR_MAX_THRESH)) THEN
        RETURN
   ENDIF

!  !--- Apply Test for Difference with NWC
   Test_Threshold = ULST_EMISS_CH7_NWC_THRESH
   IF (( Emiss_Chn7_NWC > 0) .AND. &
       ( Emiss_Chn7_NWC - Emiss_Chn7 > Test_Threshold)) THEN
       Test_Result = sym%YES
   ENDIF

   !--- Set Threshold
   Test_Threshold = ULST_EMISS_CHN7_DIFF_THRESH_OCEAN
   IF (Is_Land == sym%YES) THEN
        Test_Threshold = ULST_EMISS_CHN7_DIFF_THRESH_LAND
   ENDIF
   IF (Is_Snow == sym%YES) THEN
        Test_Threshold = ULST_EMISS_CHN7_DIFF_THRESH_SNOW
   ENDIF

   !--- Apply Test for DIFferences Relative to Clear Sky
   IF (Emiss_Chn7_Clr - Emiss_Chn7 > Test_Threshold) THEN
        Test_Result = sym%YES
   ENDIF


   !--- Bulk Threshold
   IF (Emiss_Chn7 < 0.80) THEN
        Test_Result = sym%YES
   ENDIF

   RETURN

 END FUNCTION ULST_Routine

!====================================================================
! Function Name: TERM_THERM_STAB_Routine
!
! Function:
!    Maia terminator stability test
!
! Description:
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = Term_Therm_Stab_Routine(&
!            Is_Land, &
!            Sol_Zen, &
!            BT_Chn11, &
!            BT_Chn11_1Hr, &
!            BT_Chn14, &
!            BT_Chn14_1Hr, &
!            BT_Chn15, &
!            BT_Chn15_1Hr, &
!            Cmask_1Hr)
!
! Inputs:
!   Is pixel land (YES/NO)
!   Solar Zenith Angle
!   Current and 1 hour previous 8.5 micron BT 
!   Current and 1 hour previous 11 micron BT 
!   Current and 1 hour previous 12 micron BT 
!   1 hour previous Cloud Mask
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
! Reference: MAIA paper (need to insert)
!
!====================================================================
 FUNCTION Term_Therm_Stab_Routine(&
             Is_Land, &
             Sol_Zen, &
             BT_Chn11, &
             BT_Chn11_1Hr, &
             BT_Chn14, &
             BT_Chn14_1Hr, &
             BT_Chn15, &
             BT_Chn15_1Hr, &
             Cmask_1Hr) &
      RESULT(Test_Result)

   INTEGER(KIND=INT4), INTENT(IN):: Is_Land
   REAL(KIND=REAL4), INTENT(IN):: Sol_Zen
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn11
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn11_1Hr
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn14
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn14_1Hr
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn15
   REAL(KIND=REAL4), INTENT(IN):: BT_Chn15_1Hr
   INTEGER(KIND=INT4), INTENT(IN):: Cmask_1Hr
   INTEGER(KIND=INT1):: Test_Result
   REAL(KIND=REAL4):: BT14_DIFf 
   REAL(KIND=REAL4):: Test_Value 
   REAL(KIND=REAL4):: Threshold


   !--- initialize
   Test_Result = sym%NO

   !--- exclude day time data
   IF (Sol_Zen < TERM_THERM_STAB_SOLZEN_MIN_THRESH .OR.  &
       Sol_Zen > TERM_THERM_STAB_SOLZEN_MAX_THRESH) THEN
        RETURN
   ENDIF
   
   !--- exclude large B14_DIFF
   BT14_DIFf = ABS(BT_Chn14 - BT_Chn14_1Hr)
   IF (BT14_DIFf > TERM_THERM_STAB_BT14_DIFF_THRESH) THEN
        RETURN
   ENDIF
  
   IF (Is_Land == sym%YES) THEN
        Test_Value = ABS((BT_Chn14_1Hr - BT_Chn11_1Hr) - &
            (BT_Chn14 - BT_Chn11))
        Threshold = TERM_THERM_STAB_BTD_14_11_THRESH
   
   ELSE
        Test_Value = ABS((BT_Chn14_1Hr - BT_Chn15_1Hr) - &
            (BT_Chn14 - BT_Chn15))
        Threshold = TERM_THERM_STAB_BTD_14_15_THRESH
   ENDIF
   
   !--- Apply Test 
   IF (Test_Value < Threshold .and. Cmask_1Hr == sym%CLOUDY) THEN
        Test_Result = sym%YES
   ENDIF


   RETURN

 END FUNCTION Term_Therm_Stab_Routine


!====================================================================
! Function Name: CIRH2O_Routine
!
! Function:
!    Maia terminator stability test
!
! Description:
!   
! Calling Sequence:
!   Test_Results(Num_Tests,Elem_Idx,Line_Idx) = CIRH2O_Routine(&
!
! Inputs:
!
! Outputs: 
!   Function returns pass (sym%YES) or fail (sym%NO) result of the test via Test_Result
!
! Dependencies:
!        Cloud Mask threshold include file with various needed thresholds. 
!
! Restrictions:  None
!
! Reference: None
!
!====================================================================


 FUNCTION CIRH2O_Routine( &
                         Sfc_Hgt, &
                         WaterVapor_Window_Correlation, &
                         WaterVapor_Stddev, &
                         Window_Stddev, &
                         Total_Precipitable_Water, &
                         Cos_Sen_Zen)&
                         RESULT(Test_Result) 

   REAL(KIND=REAL4), INTENT(IN):: Sfc_Hgt
   REAL(KIND=REAL4), INTENT(IN):: WaterVapor_Window_Correlation
   REAL(KIND=REAL4), INTENT(IN):: WaterVapor_Stddev
   REAL(KIND=REAL4), INTENT(IN):: Window_Stddev
   REAL(KIND=REAL4), INTENT(IN):: Total_Precipitable_Water
   REAL(KIND=REAL4), INTENT(IN):: Cos_Sen_Zen
   INTEGER(KIND=INT1):: Test_Result

   Test_Result = sym%NO

   IF (Total_Precipitable_Water / Cos_Sen_Zen < CIRH2O_TPW_THRESH) THEN
      RETURN
   ENDIF

   IF (WaterVapor_Window_Correlation > CIRH2O_CORRELATION_THRESH .AND.  &
       WaterVapor_Stddev > CIRH2O_BT_CHN10_STDDEV_THRESH .AND. &
       Window_Stddev > CIRH2O_BT_CHN14_STDDEV_THRESH .AND. &
       Sfc_Hgt < CIRH2O_SFC_HGT_LIMIT) THEN
        Test_Result = sym%YES
   ENDIF


 END FUNCTION CIRH2O_Routine

!====================================================================
! Function Name: Pearson_Corr
!
! Function:
!    Compute the Pearson Correlation Coefficient for two mxn arrays
!
! Description: Pearson's product-moment coefficient
!   
! Calling Sequence: BT_WV_BT_Window_Corr(Elem_Idx,Line_Idx) = Pearson_Corr( &
!                       ch(28)%Bt_Toa(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                       ch(31)%Bt_Toa(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
!                       Bad_Pixel_Mask(Array_Right:Array_Left,Array_Top:Array_Bottom), &
!                       Bad_Pixel_Mask(Array_Right:Array_Left,Array_Top:Array_Bottom), &
!                      Array_Width, Array_Hgt)
!   
!
! Inputs:
!   Array 1
!   Array 2
!   Elem_size
!   Line_size
!
! Outputs: 
!   Pearson Correlation coefficent
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Standard definition for Pearson correlation
!
!====================================================================

FUNCTION Pearson_Corr(Array_One,Array_Two,Bad_Pixel_One, &
                      Bad_Pixel_Two, &
                      Array_Width,Array_Hght) RESULT(Pearson_Corr_Coeff)
   REAL(KIND=REAL4), INTENT(IN), DIMENSION(:,:):: Array_One
   REAL(KIND=REAL4), INTENT(IN), DIMENSION(:,:):: Array_Two
   INTEGER(KIND=INT1), INTENT(IN), DIMENSION(:,:):: Bad_Pixel_One
   INTEGER(KIND=INT1), INTENT(IN), DIMENSION(:,:):: Bad_Pixel_Two
   INTEGER(KIND=INT4), INTENT(IN):: Array_Width
   INTEGER(KIND=INT4), INTENT(IN):: Array_Hght
   REAL(KIND=REAL4), DIMENSION(Array_Width,Array_Hght):: Pearson_Corr_Term_1
   REAL(KIND=REAL4), DIMENSION(Array_Width,Array_Hght):: Pearson_Corr_Term_2
   REAL(KIND=REAL8):: Pearson_Corr_Top_Term_1
   REAL(KIND=REAL8):: Pearson_Corr_Top_Term_2
   REAL(KIND=REAL8):: Pearson_Corr_Bottom_Term_1
   REAL(KIND=REAL8):: Pearson_Corr_Bottom_Term_2
   REAL(KIND=REAL4):: Pearson_Corr_Coeff
   REAL(KIND=REAL8):: Mean_Array_One
   REAL(KIND=REAL8):: Mean_Array_Two
   REAL(KIND=REAL8):: Sum_Array_One
   REAL(KIND=REAL8):: Sum_Array_Two

   !--- skip computation for pixel arrays with any missing data
   IF (sum(Bad_Pixel_One) > 0 .OR. sum(Bad_Pixel_Two) > 0) THEN
      Pearson_Corr_Coeff = Missing_Value_Real4
      RETURN
   ENDIF


   Sum_Array_One = sum(Array_One)
   Sum_Array_Two = sum(Array_Two)

   Mean_Array_One = Sum_Array_One / (Array_Width*Array_Hght)
   Mean_Array_Two = Sum_Array_Two / (Array_Width*Array_Hght)

   Pearson_Corr_Term_1 = Array_One - Mean_Array_One
   Pearson_Corr_Term_2 = Array_Two - Mean_Array_Two

   Sum_Array_One = sum(Pearson_Corr_Term_1)
   Sum_Array_Two = sum(Pearson_Corr_Term_2)

   Mean_Array_One = 0.0
   Mean_Array_Two = 0.0

   Pearson_Corr_Top_Term_1 = sum(Pearson_Corr_Term_1*Pearson_Corr_Term_2)
   
   Pearson_Corr_Top_Term_2 = (Sum_Array_One*Sum_Array_Two) / (Array_Width*Array_Hght)
   
   Pearson_Corr_Bottom_Term_1 = sum(Pearson_Corr_Term_1**2) - &
                                ((Sum_Array_One)**2) / (Array_Width*Array_Hght)
                                 
   Pearson_Corr_Bottom_Term_2 = sum(Pearson_Corr_Term_2**2) - &
                                ((Sum_Array_Two)**2) / (Array_Width*Array_Hght)

   Pearson_Corr_Coeff = (Pearson_Corr_Top_Term_1 - Pearson_Corr_Top_Term_2) / &
                         sqrt(Pearson_Corr_Bottom_Term_1 * &
                              Pearson_Corr_Bottom_Term_2)
   
 END FUNCTION Pearson_Corr
 
!====================================================================
! Function Name: Term_Refl_Norm
!
! Function:
!    Renormalize reflectances to improve performance near the terminator 
! using the parameteization given by Li et. al. 2006
!
! Description: Renormalizes reflectances in the terminator region
!   
! Calling Sequence: Refl_Chn2 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn2)
!   
!
! Inputs:
!   Cosine of the Solar Zenith Angle
!   Channel 2 reflectance
!
! Outputs: 
!   Renormalized reflectance
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Li et. al. 2006
!
!====================================================================
 FUNCTION Term_Refl_Norm(Cos_Sol_Zen,Reflectance)  &
          RESULT(Reflectance_Normalized)

   REAL(KIND=REAL4), INTENT(IN):: Cos_Sol_Zen
   REAL(KIND=REAL4), INTENT(IN):: Reflectance
   REAL(KIND=REAL4):: Reflectance_Normalized
   REAL(KIND=REAL4):: Norm_Param

   Reflectance_Normalized = Reflectance * Cos_Sol_Zen

   Norm_Param = 24.35 / (2*Cos_Sol_Zen + sqrt(498.5225*(Cos_Sol_Zen**2) + 1) )

   Reflectance_Normalized = Reflectance_Normalized*Norm_Param 

 END FUNCTION Term_Refl_Norm


!--------------------------------------------------------------------------
! End ACM test Subroutines and functions 
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! End of the Module
!--------------------------------------------------------------------------
END MODULE Baseline_Cloud_Mask
