MODULE AWG_Baseline_Cloud_Mask
!$Id: awg_baseline_cloud_mask.f90,v 1.31 2011/04/18 16:14:36 wstraka Exp $
!
! Module Name:
!   AWG_Baseline_Cloud_Mask
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
! CALLing Sequence:
!     Use AWG_Baseline_Cloud_Mask
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
!   AWG_Baseline_Cloud_Mask_Main - main routine for cloud Mask
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
use calibration_constants, only: &
        sun_earth_distance  &
        , solar_ch20_nu


IMPLICIT NONE

PRIVATE

PUBLIC:: AWG_Baseline_Cloud_Mask_Main

PRIVATE:: Compute_Probably_Clear_Restoral
PRIVATE:: Compute_Probably_Cloudy
!PRIVATE:: Clear_Chn2_Reflectance_Field
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
INCLUDE 'awg_cloud_mask_tests_thresholds.inc'

CONTAINS
INCLUDE 'awg_cloud_mask_tests_code.inc'

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
! 29 - blank                                      byte 4 bit 4
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
SUBROUTINE AWG_Baseline_Cloud_Mask_Main()

   CHARACTER(len=*), PARAMETER:: Routine_Name =  &
                                 "AWG_Baseline_Cloud_Mask: test_AWG_Baseline_Cloud_Mask_main"
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
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: BT_WaterVapor_Stddev_3x3
   INTEGER(KIND=INT4), DIMENSION(:,:), ALLOCATABLE:: X_NWC_Idx
   INTEGER(KIND=INT4), DIMENSION(:,:), ALLOCATABLE:: Y_NWC_Idx

   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Refl_Chn2_Mean_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Refl_Chn2_Max_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Refl_Chn2_Min_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Refl_Chn2_Stddev_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Sfc_Hgt_Mean_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Sfc_Hgt_Max_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Sfc_Hgt_Min_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE:: Sfc_Hgt_Stddev_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Refl_Chn2_Clear_Mean_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Refl_Chn2_Clear_Max_3x3
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Refl_Chn2_Clear_Min_3X3
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER:: Refl_Chn2_Clear_Stddev_3X3
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
   INTEGER(KIND=INT1), DIMENSION(:,:,:), ALLOCATABLE, target:: Test_Results_Temp


   !
   !--- local pointers
   !
   INTEGER(KIND=INT1), DIMENSION(:,:), POINTER :: Cloud_Mask
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
   Max_Num_Lines_per_Seg = size(Nav%Lat,2)
   
   ! Need to have an array for the unpacked tests
   ALLOCATE(Test_Results_Temp(Num_Elem, Number_of_Lines_in_this_Segment, Total_Num_Tests))
   


   !set the initial cloud mask algorithm name here (for usage in the temporal term test)
   Algo_Name = 'awg_baseline_cmask_'

   !--- store name of sensor
   Sat_Name = trim(Sensor%Platform_Name)
   
   !-----------------------------------------------------------------
   ! Set sensor dependent thresholds
   !  This is due to a SRF issue on the 3.9 micron band
   !  Also compose sensor dependent algorithm name for temporal term test
   !-----------------------------------------------------------------
   
   CALL Set_Cmask_Thresholds(Sat_Name, Algo_Name)
   

   !-----------------------------------------------------------------
   ! Read IN background Chn2 Reflectance Field 
   ! Calculated by CLAVR-x
   !-----------------------------------------------------------------
                                     
    Refl_Chn2_Clear => ch(1)%Ref_Toa_Clear
    Refl_Chn2_Clear_Max_3x3 => Ref_Ch1_Clear_Max_3x3
    Refl_Chn2_Clear_Min_3x3 => Ref_Ch1_Clear_Min_3x3

   !-----------------------------------------------------------------
   ! Initialize masks
   !-----------------------------------------------------------------
   Cloud_Mask => Cld_Mask
   Cloud_Mask_Packed => Cld_Test_Vector_Packed
   Test_Results => Test_Results_Temp !will just null for now
   Cloud_Mask_QF => null() !will just null for now
   Cloud_Mask_Tmpy => One_Byte_Temp !need to find
   
   
   !inputs
   X_LRC_Idx => I_LRC
   Y_LRC_Idx => J_LRC
!   LRC_Mask => out2(Algo_Num)%LRC_Mask
   Emiss_Tropo_Chn14 => ch(31)%Emiss_Tropo
   
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
   ! Compute spatial uniformity metrics - Already done by CLAVRx.
   ! But set appropriate WV/IR 3x3
   !=======================================================================
   IF ((Sensor%Chan_On_Flag_Default(27) > 0) .OR. (Sensor%Chan_On_Flag_Default(29) > 0) .AND. &
        (Sensor%Chan_On_Flag_Default(31) > 0)) THEN
        
         IF (Sensor%Chan_On_Flag_Default(27) > 0) THEN
            BT_WaterVapor_Stddev_3x3 => Btd_Ch31_Ch27_Std_3x3
         ELSE
            BT_WaterVapor_Stddev_3x3 => Btd_Ch31_Ch29_Std_3x3
         
         ENDIF


   ENDIF

   !=======================================================================
   ! Load temporal information   (status = sym%FAILURE or sym%SUCCESS)
   !=======================================================================

   !--- previous 11 micron temp
   !
   
!   Have_Prev_BT_Chn14_15min = Load_Temporal_Data(minus_i_15min,"",           &
!                                            temporal(minus_i_15min)%bt14)
!   Have_Prev_BT_Chn14_Clr_15min = Load_Temporal_Data(minus_i_15min,"",       &
!                                        temporal(minus_i_15min)%bt_clr14)

   Have_Prev_BT_Chn14_15min = sym%FAILURE
   Have_Prev_BT_Chn14_Clr_15min = sym%FAILURE

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


    Have_Prev_BT_Chn11_1Hr = sym%FAILURE
    Have_Prev_BT_Chn14_1Hr = sym%FAILURE
    Have_Prev_BT_Chn15_1Hr = sym%FAILURE
    Have_Prev_Cmask_1Hr = sym%FAILURE



   !=======================================================================
   ! Compute 11 micron emissivity at Tropopause
   ! -------------- removed - CLAVR-x computes emiss11 already------------!
   !=======================================================================



   !======================================================================
   ! compute local radiative center
   !======================================================================
   
   !-------------- removed - CLAVR-x computes LRC already ------------------!

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
   ! WCS - must remain to match GS implmentation
   !----------------------------------------------------------------------
   IF ((Sensor%Chan_On_Flag_Default(27) > 0) .OR. (Sensor%Chan_On_Flag_Default(29) > 0) .AND. &
        (Sensor%Chan_On_Flag_Default(31) > 0)) THEN
   
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

                    IF (Sensor%Chan_On_Flag_Default(10) > 0) THEN
                         BT_WV_BT_Window_Corr(Elem_Idx,Line_Idx) = Pearson_Corr(&
                                                                    ch(29)%Bt_Toa(Array_Right:Array_Left,Array_Top:Array_Bottom), &
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
         IF (Geo%Satzen(Elem_Idx,Line_Idx) > SENSOR_ZEN_THRESH) THEN

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
         IF (Sensor%Chan_On_Flag_Default(14) > 0 .and. Sensor%Chan_On_Flag_Default(15) > 0) THEN
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
         Sfc_Temp_Uni_NWP =  tmpsfc_uni_nwp(Elem_Idx,Line_Idx)
     
         Sfc_Temp         =  Tsfc_Nwp_Pix(Elem_Idx,Line_Idx)

         !
         !--- Total Precipitable Water (g/m^2 or cm)
         !
         Total_Precipitable_Water_NWP =  Tpw_Nwp_Pix(Elem_Idx,Line_Idx)

         !
         !--- Total Ozone (Dobson Unit)
         !
         Total_Ozone_Path_NWP =  Ozone_Nwp_Pix(Elem_Idx,Line_Idx)

         !
         !---sun earth distance factor
         !
         Sun_Earth_Dist =           Sun_Earth_Distance

         !
         !---glint zenith angle
         !
         Glint_Zen =      Geo%Glintzen(Elem_Idx,Line_Idx)       

         !
         !---land type
         !
         Land_Type =      Sfc%Land_Mask(Elem_Idx,Line_Idx)    

         !
         !---coast type
         !
         Coast_Type =     Sfc%Coast_Mask(Elem_Idx,Line_Idx)    

         !
         !---Surface Height 
         !
         Sfc_Hgt =      Sfc%Zsfc(Elem_Idx,Line_Idx)         

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
         Desert_Mask =    Sfc%desert_mask(Elem_Idx,Line_Idx)    

         !
         !---snow Mask
         !
         Snow_Mask =      Sfc%snow(Elem_Idx,Line_Idx)     

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
            !--- already done in CLAVRx
            !IF (Sol_Zen > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) THEN
            !  Refl_Chn2 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn2)
            !ENDIF

         ENDIF

         !--- Channel 4 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(2) > 0) THEN

            Is_Chn(4) = sym%YES

            !
            !---1.38 micron reflectance
            !
            Refl_Chn4 =          ch(2)%Ref_Toa(Elem_Idx,Line_Idx)       

            !--- renormalize for improved terminator performance
            !--- already done in CLAVRx
            !IF (Sol_Zen > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) THEN
            !  Refl_Chn4 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn4)
            !ENDIF

         ENDIF

         !--- Channel 5 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(6) > 0) THEN

            Is_Chn(5) = sym%YES

            !
            !---1.60 micron reflectance
            !
            Refl_Chn5 =          ch(6)%Ref_Toa(Elem_Idx,Line_Idx)          

            !--- renormalize for improved terminator performance
            !--- already done in CLAVRx
            !IF (Sol_Zen > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) THEN
            !  Refl_Chn5 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn5)
            !ENDIF

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
            Rad_Chn7  =          ch(20)%Rad_Toa(Elem_Idx,Line_Idx)        

            !
            !---3.9 um pseudo reflectance
            !
            Refl_Chn7 =         ch(20)%Ref_Toa(Elem_Idx,Line_Idx)      

            !--- renormalize for improved terminator performance
            !--- already done in CLAVRx
            !IF (Sol_Zen > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) THEN
            !  Refl_Chn7 = Term_Refl_Norm(Cos_Sol_Zen,Refl_Chn7)
            !ENDIF

            !
            !---3.9 um emissivity
            !

            Emiss_Chn7 =        Rad_Ch20_ems(Elem_Idx,Line_Idx)
                  
            Emiss_Chn7_NWC = Missing_Value_Real4
            IF (Elem_NWC_Idx > 0 .and. Line_NWC_Idx > 0) THEN
                Emiss_Chn7_NWC = Rad_Ch20_ems(Elem_NWC_Idx,Line_NWC_Idx) !Need this
            ENDIF
           
            !
            !---3.9 um surface emissivity
            !
            Sfc_Emiss_Chn7_RTM =   ch(20)%sfc_emiss(Elem_Idx,Line_Idx)    

            !
            !---clear 4 micron radiance
            !
            Rad_Chn7_Clr =      ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx)     

            !
            !--- Atmospheric Transmission IN Channel 7
            !
            Atm_Trans_Chn7_RTM =   Rtm(X_NWP_Idx,Y_NWP_Idx)%d(View_Zen_Idx)%ch(20)%Trans_Atm_Profile(Sfc_Idx_NWP)
            
            !
            !---solar energy IN channel 7
            !
            Chn7_Sol_Energy =       solar_ch20_nu    

         ENDIF


         !--- Channel 9 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(27) > 0) THEN

            Is_Chn(9) = sym%YES

            !
            !---6.7 micron bt
            !
            BT_Chn9  =         ch(27)%Bt_Toa(Elem_Idx,Line_Idx)

         ENDIF

         !--- Channel 10 Aliases and Derived Parameters
         IF (Sensor%Chan_On_Flag_Default(29) > 0) THEN

            Is_Chn(10) = sym%YES

            !
            !---7.3 micron bt
            !
            BT_Chn10  =         ch(29)%Bt_Toa(Elem_Idx,Line_Idx)

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
         
         !--- WCS3 - REMOVED BECAUSE CLAVRX HAS NO TEMPORAL ------------------!

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

            IF ((Is_Chn(Chn_Num) > 0).AND.                  &
                (Bad_Pixel_Mask(Elem_Idx,Line_Idx)==sym%YES)) THEN

               Is_Chn(Chn_Num) = sym%NO

            ENDIF

           ENDDO

         ENDIF  

         !
         !--- thermal channels (7-16)
         !
         DO Chn_Num = 7, Nchan_Max

            IF ((Is_Chn(Chn_Num) > 0).AND.                  &
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
         IF (Sensor%Chan_On_Flag_Default(31) > 0) THEN
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
         
         IF(Sol_Zen < Day_Sol_Zen_Thresh) THEN

            !--- most of this is done behind the scenes of CLAVR-x, so commented out
            !--- and set to CLAVRx globals
!            Refl_Chn2_Clear(Elem_Idx,Line_Idx) = ch(1)%Ref_Toa_Clear(Elem_Idx,Line_Idx) 
                    
!            Refl_Chn2_Clear_Max_3x3(Elem_Idx,Line_Idx) = Ref_Ch1_Clear_Max_3x3(Elem_Idx,Line_Idx)
                    
!            Refl_Chn2_Clear_Min_3x3(Elem_Idx,Line_Idx) = Ref_Ch1_Clear_Min_3x3(Elem_Idx,Line_Idx)

            !--- renormalize for improved terminator performance
            !--- Already done by CLAVRX, so removing calls

         ENDIF

         
         

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
                              Bt_Ch31_Std_3x3(Elem_Idx,Line_Idx), &
                              Sfc_Hgt_Stddev_3x3(Elem_Idx,Line_Idx)) ! need - WCS

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
                              Bt_Ch31_Min_3x3(Elem_Idx,Line_Idx), &
                              Bt_Ch31_Max_3x3(Elem_Idx,Line_Idx), &
                              Sfc_Hgt_Stddev_3x3(Elem_Idx,Line_Idx)) ! need - WCS
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
                               Bt_Ch31_Std_3x3(Elem_Idx,Line_Idx))
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
                      Bt_Ch31_Std_3x3(Elem_Idx,Line_Idx)) 
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
                      Bt_Ch31_Std_3x3(Elem_Idx,Line_Idx), &
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
                         Sfc_Hgt_Max_3x3(Elem_Idx,Line_Idx), &  !need
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
            Planck_Emission_Chn7_Clr = planck_rad_fast(20,BT_Chn14_Clr)
            Solar_Rad_Chn7_Clr = Rad_Chn7_Clr

            IF (Is_Day == sym%YES .or. Is_Terminator == sym%YES) THEN
              Solar_Rad_Chn7_Clr = Rad_Chn7_Clr + (1.0 - Sfc_Emiss_Chn7_RTM)         &
                     * Atm_Solar_Trans_Chn7 * max(Cos_Sol_Zen,0.05) * (Chn7_Sol_Energy/PI)
            ENDIF

            Solar_BT_Chn7_Clr = planck_temp_fast(20,Solar_Rad_Chn7_Clr)

            Emiss_Chn7_Clr = Solar_Rad_Chn7_Clr / Planck_Emission_Chn7_Clr

            Test_Results(Num_Tests,Elem_Idx,Line_Idx) = EMISS4_Routine(Is_Glint,  &
                         Is_Coast, &
                         Is_Land, &
                         Is_Snow, &
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
        
         IF ((Sensor%Chan_On_Flag_Default(7) > 0) .AND. (Is_Chn(7) == sym%NO)) THEN
             
             Cloud_Mask_QF(Elem_Idx,Line_Idx) = REDUCED_QUAL_BAD_CHN7
         
         ELSE IF ((Sensor%Chan_On_Flag_Default(2) > 0) .AND. &
                  (Is_Day == sym%YES) .AND. (Is_Chn(2) == sym%NO)) THEN
         
             Cloud_Mask_QF(Elem_Idx,Line_Idx) = REDUCED_QUAL_BAD_CHN2
             
          ELSE IF (((Is_Day == sym%YES) .AND. &
                    (((Sensor%Chan_On_Flag_Default(4) > 0) .AND. (Is_Chn(4) == sym%NO)) .OR. &
                    ((Sensor%Chan_On_Flag_Default(5) > 0) .AND. (Is_Chn(5) == sym%NO))))  .OR. &
                    ((Sensor%Chan_On_Flag_Default(10) > 0) .AND. (Is_Chn(10) == sym%NO))  .OR. &             
                    ((Sensor%Chan_On_Flag_Default(11) > 0) .AND. (Is_Chn(11) == sym%NO))  .OR. &             
                    ((Sensor%Chan_On_Flag_Default(15) > 0) .AND. (Is_Chn(15) == sym%NO))) THEN

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
   ! nullify local pointers
   !======================================================================

   Cloud_Mask => null()
   Cloud_Mask_IR => null()
   Cloud_Mask_Packed => null()
   Test_Results => null()
   Cloud_Mask_QF => null()
   Cloud_Mask_Tmpy => null()
   LRC_Mask => null()
   X_LRC_Idx => null()
   Y_LRC_Idx => null()
   Emiss_Tropo_Chn14 => null()
   Refl_Chn2_Clear =>null()
   Refl_Chn2_Clear_Max_3x3 => null()
   Refl_Chn2_Clear_Min_3x3 => null()
   BT_WaterVapor_Stddev_3x3 =>null()


   !Initialize total deallocation status 
   Alloc_Status_Total = 0 
   
   IF (ALLOCATED(X_NWC_Idx)) DEALLOCATE(X_NWC_Idx,stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status
   
   IF (ALLOCATED(Y_NWC_Idx)) DEALLOCATE(Y_NWC_Idx,stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status
   
   IF (ALLOCATED(BT_WV_BT_Window_Corr)) DEALLOCATE(BT_WV_BT_Window_Corr,stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status
   
   IF (ALLOCATED(Test_Results_Temp)) DEALLOCATE(Test_Results_Temp,stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status   

   !======================================================================
   ! deallocate 2nd darkest composite data array
   !======================================================================
!   IF (ALLOCATED(Refl_Chn2_Clear)) DEALLOCATE(Refl_Chn2_Clear,stat=Alloc_Status)
   Alloc_Status_Total = Alloc_Status_Total + Alloc_Status

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


!
!---END SUBROUTINE AWG_Baseline_Cloud_Mask_main
!


END SUBROUTINE AWG_Baseline_Cloud_Mask_Main


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
! COMMENTED OUT SINCE CLAVR-X CALCULATES THIS
!
!====================================================================
! SUBROUTINE Clear_Chn2_Reflectance_Field(Num_Elem, &
!                                         Max_Num_Lines_per_Seg, &
!                                         Refl_Chn2_Clear)

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
       !--- use MODIS white sky from mapped data if available
!       WHERE(sat%ws_albedo2 /= Missing_Value_Real4)
!          Refl_Chn2_Clear = sat%ws_albedo2
!       ENDWHERE 
       !--- ensure consistent missing value
!       WHERE(Refl_Chn2_Clear <= 0.0)
!          Refl_Chn2_Clear = Missing_Value_Real4
!       ENDWHERE

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
            Clr_Rad_Chn14 =     ch(31)%rad_Toa_Clear(Elem_Idx,Line_Idx)

            !
            !---BB 11 um rad at tropopause
            !
            Blkbdy_Tropo_Rad_Chn14 = rtm(X_NWP_Idx,Y_NWP_Idx)%d(View_Zen_Idx)%ch(31)%Rad_BB_Cloud_Profile(Tropo_Idx_NWP)

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


!--------------------------------------------------------------------------
! End ACM test subroutines and functions 
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! End of the Module
!--------------------------------------------------------------------------
END MODULE AWG_Baseline_Cloud_Mask
