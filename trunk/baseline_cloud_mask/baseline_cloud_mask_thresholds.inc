!$Id: baseline_cloud_mask_thresholds.inc 1887 2016-11-24 20:03:09Z wstraka $
!----------------------------------------------------------------------
!--- parametes to control some internal computations (LRC and NWC)
!----------------------------------------------------------------------
 INTEGER, PRIVATE, PARAMETER :: NWC_PIXEL_RADIUS = 10
 REAL, PRIVATE, PARAMETER :: EMISS_TROPO_CHN14_GRADIENT_MIN = 0.0
 REAL, PRIVATE, PARAMETER :: EMISS_TROPO_CHN14_GRADIENT_MAX = 1.0
 REAL, PRIVATE, PARAMETER :: EMISS_TROPO_CHN14_GRADIENT_THRESH = 0.75

!----------------------------------------------------------------------
!--- determine positions of tests on test vector
!----------------------------------------------------------------------
 INTEGER, PRIVATE, PARAMETER :: FIRST_ANCIL_TEST = 1
 INTEGER, PRIVATE, PARAMETER :: NUM_ANCIL_TESTS = 9
 INTEGER, PRIVATE, PARAMETER :: LAST_ANCIL_TEST = FIRST_ANCIL_TEST + NUM_ANCIL_TESTS - 1

 INTEGER, PRIVATE, PARAMETER :: FIRST_CLR_UNI_TEST = NUM_ANCIL_TESTS + 1
 INTEGER, PRIVATE, PARAMETER :: NUM_CLR_UNI_TESTS = 2
 INTEGER, PRIVATE, PARAMETER :: LAST_CLR_UNI_TEST = FIRST_CLR_UNI_TEST + NUM_CLR_UNI_TESTS - 1

 INTEGER, PRIVATE, PARAMETER :: FIRST_IR_CLD_MASK_TEST = LAST_CLR_UNI_TEST + 1
 INTEGER, PRIVATE, PARAMETER :: NUM_IR_CLD_MASK_TESTS = 8
 INTEGER, PRIVATE, PARAMETER :: LAST_IR_CLD_MASK_TEST = FIRST_IR_CLD_MASK_TEST + NUM_IR_CLD_MASK_TESTS - 1

 INTEGER, PRIVATE, PARAMETER :: FIRST_VIS_CLD_MASK_TEST = LAST_IR_CLD_MASK_TEST + 1
 INTEGER, PRIVATE, PARAMETER :: NUM_VIS_CLD_MASK_TESTS = 2
 INTEGER, PRIVATE, PARAMETER :: LAST_VIS_CLD_MASK_TEST = FIRST_VIS_CLD_MASK_TEST + NUM_VIS_CLD_MASK_TESTS - 1

 INTEGER, PRIVATE, PARAMETER :: FIRST_SWIR_SOLAR_CLD_MASK_TEST = LAST_VIS_CLD_MASK_TEST + 1
 INTEGER, PRIVATE, PARAMETER :: NUM_SWIR_SOLAR_CLD_MASK_TESTS = 2
 INTEGER, PRIVATE, PARAMETER :: LAST_SWIR_SOLAR_CLD_MASK_TEST = &
                                FIRST_SWIR_SOLAR_CLD_MASK_TEST + NUM_SWIR_SOLAR_CLD_MASK_TESTS - 1

 INTEGER, PRIVATE, PARAMETER :: FIRST_SWIR_THERMAL_CLD_MASK_TEST = LAST_SWIR_SOLAR_CLD_MASK_TEST + 1
 INTEGER, PRIVATE, PARAMETER :: NUM_SWIR_THERMAL_CLD_MASK_TESTS = 2
 INTEGER, PRIVATE, PARAMETER :: LAST_SWIR_THERMAL_CLD_MASK_TEST = FIRST_SWIR_THERMAL_CLD_MASK_TEST + &
                                    NUM_SWIR_THERMAL_CLD_MASK_TESTS - 1

 INTEGER, PRIVATE, PARAMETER :: NUM_CLD_MASK_TESTS = NUM_VIS_CLD_MASK_TESTS + NUM_IR_CLD_MASK_TESTS +  &
                                             NUM_SWIR_SOLAR_CLD_MASK_TESTS + NUM_SWIR_THERMAL_CLD_MASK_TESTS

 INTEGER, PRIVATE, PARAMETER :: FIRST_CLD_MASK_TEST = FIRST_IR_CLD_MASK_TEST
 INTEGER, PRIVATE, PARAMETER :: LAST_CLD_MASK_TEST = LAST_SWIR_THERMAL_CLD_MASK_TEST
 
 INTEGER, PRIVATE, PARAMETER :: FIRST_RES_TEST = NUM_ANCIL_TESTS + NUM_CLR_UNI_TESTS + NUM_CLD_MASK_TESTS + 1
 INTEGER, PRIVATE, PARAMETER :: NUM_RES_TESTS = 3

 INTEGER, PRIVATE, PARAMETER :: NUM_NONCLD_TESTS = 4
 INTEGER, PRIVATE, PARAMETER :: FIRST_NONCLD_TEST = NUM_ANCIL_TESTS + NUM_CLR_UNI_TESTS +  &
                               NUM_CLD_MASK_TESTS + NUM_RES_TESTS + 1
 INTEGER, PRIVATE, PARAMETER :: LAST_NONCLD_TEST = FIRST_NONCLD_TEST - 1

 INTEGER, PRIVATE, PARAMETER :: TOTAL_NUM_TESTS = NUM_ANCIL_TESTS + NUM_CLR_UNI_TESTS +  &
                               NUM_CLD_MASK_TESTS + NUM_RES_TESTS + NUM_NONCLD_TESTS

!----------------------------------------------------------------------
!--- viewing geometry thresholds
!----------------------------------------------------------------------
 REAL, PRIVATE, PARAMETER:: SENSOR_ZEN_THRESH = 70.0
 REAL, PRIVATE, PARAMETER:: DAY_SOL_ZEN_THRESH = 87.0
 REAL, PRIVATE, PARAMETER:: TERMINATOR_MIN_SOL_ZEN_THRESH = 87.0
 REAL, PRIVATE, PARAMETER:: TERMINATOR_MAX_SOL_ZEN_THRESH = 93.0
 REAL, PRIVATE, PARAMETER:: RGCT_SOL_ZEN_THRESH = 80.0
 REAL, PRIVATE, PARAMETER:: NIRREF_CHN5_SOL_ZEN_THRESH = 80.0
 REAL, PRIVATE, PARAMETER:: NIRREF_CHN7_SOL_ZEN_THRESH = 80.0
 REAL, PRIVATE, PARAMETER:: CIRREF_SOL_ZEN_THRESH = 80.0
 REAL, PRIVATE, PARAMETER:: RUT_SOL_ZEN_THRESH = 80.0
 REAL, PRIVATE, PARAMETER:: RVCT_SOL_ZEN_THRESH = 83.0
 REAL, PRIVATE, PARAMETER:: RVCT_SCAT_ZEN_THRESH = 90.0
 REAL, PRIVATE, PARAMETER:: GLINT_ZEN_THRESH = 40.0
 REAL, PRIVATE, PARAMETER:: TERMINATOR_REFL_SOL_ZEN_THRESH = 60.0

 REAL, PRIVATE, PARAMETER:: WATER_FREEZING_POINT = 273.15
 REAL, PRIVATE, PARAMETER:: COLD_SFC_THRESH = 265.0
 REAL, PRIVATE, PARAMETER:: MAX_GLINT_CLR_OBS_BT_CHN14_DIFF = 5.0
 REAL, PRIVATE, PARAMETER:: MAX_GLINT_CLR_REL_REFL2_STDDEV_THRESH = 0.10
 REAL, PARAMETER, PRIVATE:: RAYLEIGH_OPTICAL_DEPTH_CHN2 = 0.05
 REAL, PARAMETER, PRIVATE:: RAYLEIGH_SINGLE_SCATTER_ALBEDO = 1.00

!--- Aerosol
 REAL, PARAMETER, PRIVATE:: AEROSOL_SINGLE_SCATTER_ALBEDO_CHN2 = 1.0
 REAL, PARAMETER, PRIVATE:: AEROSOL_OPTICAL_DEPTH_CHN2_LAND = 0.3
 REAL, PARAMETER, PRIVATE:: AEROSOL_OPTICAL_DEPTH_CHN2_OCEAN = 0.2
 REAL, PARAMETER, PRIVATE:: AEROSOL_ASYMMETRY_PARAMETER = 0.6

 !--- Quality Flags
 INTEGER, PRIVATE, PARAMETER :: VALID_CMASK_RETRIEVAL = 0
 INTEGER, PRIVATE, PARAMETER :: INVALID_CMASK_RETREVAL = 1
 INTEGER, PRIVATE, PARAMETER :: CMASK_OUTSIDE_SEN_ZEN_RANGE = 2
 INTEGER, PRIVATE, PARAMETER :: INVALID_CMASK_BAD_CHN14 = 3
 INTEGER, PRIVATE, PARAMETER :: REDUCED_QUAL_BAD_CHN7 = 4
 INTEGER, PRIVATE, PARAMETER :: REDUCED_QUAL_BAD_CHN2 = 5
 INTEGER, PRIVATE, PARAMETER :: REDUCED_QUAL_BAD_OTHER = 6

 !--- Binary CM values - WCS3 3/5/2012
 INTEGER, PRIVATE, PARAMETER :: CLD_MASK_BIN_CLEAR = 0
 INTEGER, PRIVATE, PARAMETER :: CLD_MASK_BIN_CLOUD = 1


!----------------------------------------------------------------------
! cloud detection thresholds
!----------------------------------------------------------------------
 REAL, PRIVATE, PARAMETER :: BT_CHN14_MAX_FMFT_THRESH = 310.0, &
                                         VALID_BT14_CLR_PIXEL_THRESH = 200.0, &
                                         BTDIFF_CHN14_CHN15_MIN_FMFT_THRESH = 0.0, &
                                         PFMFT_LAND_THRESH = 2.50, &
                                         PFMFT_OCEAN_THRESH = 0.80, &
                                         PFMFT_SNOW_THRESH = 1.0, &
                                         PFMFT_COLD_SURFACE_THRESH = 1.0, &
                                         PFMFT_BT_CHN14_STDDEV_3x3_THRESH = 0.30, &
                                         NFMFT_LAND_THRESH = -2.00, &
                                         NFMFT_OCEAN_THRESH = -1.00, &
                                         NFMFT_SNOW_THRESH = -5.00, &
                                         NFMFT_DESERT_THRESH = -5.00, &
                                         NFMFT_BTD_CHN14_CHN15_MAX_THRESH = 1.50, &
                                         BETA_TROPO_CHN14_CHN15_MIN_SWCIRRUS_THRESH = 0.95, &
                                         BETA_TROPO_CHN14_CHN15_MAX_SWCIRRUS_THRESH = 1.30, &
                                         BT_CHN14_MAX_SWCIRRUS_THRESH = 310.0, &
                                         EMISS_TROPO_CHN14_MAX_SWCIRRUS_THRESH = 0.95, &
                                         EMISS_TROPO_CHN14_MIN_SWCIRRUS_THRESH = 0.05, &
                                         BT_CHN14_MAX_TMFT_THRESH = 310.0, &
                                         BTDIFF_CHN7_CHN14_MIN_TMFT_THRESH = -10.0, &
                                         TMFT_LAND_THRESH =  5.0, & !10.0, &
                                         TMFT_OCEAN_THRESH = 3.0, & !5.0, &
                                         BT_CHN14_CLR_UNI_THRESH_LAND = 1.1, &
                                         BT_CHN14_CLR_UNI_THRESH_OCN = 0.60, &
                                         REFL_CHN2_CLR_UNI_THRESH_LAND = 0.20, &
                                         REFL_CHN2_CLR_UNI_THRESH_OCEAN = 1.0, &
                                         SFC_TEMP_UNI_NWP_THRESH = 4.0, &
                                         EMISS_CHN14_TROPO_OCN_THRESH = 0.10, &
                                         EMISS_CHN14_TROPO_LAND_THRESH = 0.30, &
                                         EMISS_CHN14_TROPO_SNOW_THRESH = 0.40, &
                                         EMISS_CHN14_TROPO_COLD_SURFACE_THRESH = 0.50, &
                                         EMISS_CHN14_TROPO_DESERT_THRESH = 0.40, &   
                                         EMISS_CHN14_TROPO_LRC_OCN_THRESH = 0.28, &
                                         EMISS_CHN14_TROPO_LRC_LAND_THRESH = 0.30, &
                                         EMISS_CHN14_TROPO_LRC_DESERT_THRESH = 0.40, &
                                         EMISS_CHN14_TROPO_LRC_SNOW_THRESH = 0.50, &
                                         EMISS_CHN14_TROPO_LRC_COLD_SURFACE_THRESH = 0.50, &    
                                         EMISS4_BT_CHN14_MAX_THRESH = 310.0, &
                                         EMISS4_SFC_EMISS_CHN7_THRESH = 0.90, &
                                         NIRREF_CHN5_REFL_THRESH_SNOW = 15.0, &      !NIRREF threshold on band 5 
                                         NIRREF_CHN7_REFL_THRESH_SNOW = 15.0, &      !NIRREF threshold on band 7 
                                         NIRREF_NDSI_THRESH_SNOW = 0.5, &            !NIRREF threshold on NDSI
                                         CIRREF_THRESH = 5.0, &                      !CIRREF threshold on band 4
                                         CIRREF_SFC_HGT_LIMIT = 2000.0, &     !maximum sfc alt (m) to apply CIRREF
                                         CIRH2O_SFC_HGT_LIMIT = 2000.0, &     !maximum sfc alt (m) to apply CIRH2O
                                         NIRREF_SFC_HGT_LIMIT = 1000.0, &     !maximum sfc alt (m) to apply NIRREF
                                         ULST_EMISS_CHN7_CLR_MIN_THRESH = 0.85, &
                                         ULST_EMISS_CHN7_CLR_MAX_THRESH = 1.25, &
                                         ULST_EMISS_CHN7_SFC_THRESH = 0.90, &
                                         ULST_EMISS_CHN7_MAX_THRESH = 0.95, &
                                         ULST_BT_CHN14_MAX_THRESH = 290.00, &
                                         ULST_EMISS_CH7_NWC_THRESH = 0.075, &
                                         RFMFT_HI_LAND_THRESH = 1.00, &
                                         RFMFT_HI_OCEAN_THRESH = 0.70, &
                                         RFMFT_BT_CHN14_MAX_THRESH=300.0, &
                                         RFMFT_BTDIFF_CHN14_CHN15_MIN_THRESH=1.0, &
                                         CIRH2O_CORRELATION_THRESH = 0.7, &
                                         CIRH2O_BT_CHN10_STDDEV_THRESH = 0.5, &
                                         CIRH2O_BT_CHN14_STDDEV_THRESH = 0.5, &
                                         CIRH2O_TPW_THRESH = 0.30, &
                                         TEMPIR_BT_CHN14_15MIN_TEMPORAL_OFFSET = 2.0, &
                                         TERM_THERM_STAB_SOLZEN_MIN_THRESH = 80.0, & 
                                         TERM_THERM_STAB_SOLZEN_MAX_THRESH = 93.0, & 
                                         TERM_THERM_STAB_BT14_DIFF_THRESH = 1.0, & 
                                         TERM_THERM_STAB_BTD_14_11_THRESH = 0.5, & 
                                         TERM_THERM_STAB_BTD_14_15_THRESH = 0.6 

 REAL, PRIVATE, PARAMETER :: RTCT_OCN_THRESH = 3.2   
 REAL, PRIVATE, PARAMETER :: RTCT_LAND_THRESH = 4.1


 !--- Sensor specific thresholds 
 REAL, PRIVATE, PARAMETER :: &
                                         EMISS4_EMISS_CHN7_OCN_THRESH_MODIS = 0.22, &
                                         EMISS4_EMISS_CHN7_LAND_THRESH_MODIS = 0.40, &
                                         EMISS4_EMISS_CHN7_SNOW_THRESH_MODIS = 0.40, &
                                         EMISS4_EMISS_CHN7_DESERT_THRESH_MODIS = 0.5, &
                                         
                                         EMISS4_EMISS_CHN7_OCN_THRESH_AVHRR = 0.22, &
                                         EMISS4_EMISS_CHN7_LAND_THRESH_AVHRR = 0.40, &
                                         EMISS4_EMISS_CHN7_SNOW_THRESH_AVHRR = 0.40, &
                                         EMISS4_EMISS_CHN7_DESERT_THRESH_AVHRR = 0.5, &
                                         
                                         EMISS4_EMISS_CHN7_OCN_THRESH_GOES = 0.1, &
                                         EMISS4_EMISS_CHN7_LAND_THRESH_GOES = 0.46, &
                                         EMISS4_EMISS_CHN7_SNOW_THRESH_GOES = 0.30, &
                                         EMISS4_EMISS_CHN7_DESERT_THRESH_GOES = 0.6, &                                         

                                         EMISS4_EMISS_CHN7_OCN_THRESH_SEVIRI = 0.1, &
                                         EMISS4_EMISS_CHN7_LAND_THRESH_SEVIRI = 0.46, &
                                         EMISS4_EMISS_CHN7_SNOW_THRESH_SEVIRI = 0.40, &
                                         EMISS4_EMISS_CHN7_DESERT_THRESH_SEVIRI = 0.6, &

                                         EMISS4_EMISS_CHN7_OCN_THRESH_MTSAT = 0.1, &
                                         EMISS4_EMISS_CHN7_LAND_THRESH_MTSAT = 0.46, &
                                         EMISS4_EMISS_CHN7_SNOW_THRESH_MTSAT = 0.40, &
                                         EMISS4_EMISS_CHN7_DESERT_THRESH_MTSAT = 0.6, &                                         

                                         ULST_EMISS_CHN7_DIFF_THRESH_LAND_GOES = 0.10, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_OCN_GOES = 0.12, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_SNOW_GOES = 0.12, &
                                         
                                         ULST_EMISS_CHN7_DIFF_THRESH_LAND_SEVIRI = 0.20, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_OCN_SEVIRI = 0.20, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_SNOW_SEVIRI = 0.12, &
                                         
                                         ULST_EMISS_CHN7_DIFF_THRESH_LAND_AVHRR = 0.12, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_OCN_AVHRR = 0.12, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_SNOW_AVHRR = 0.12, &

                                         ULST_EMISS_CHN7_DIFF_THRESH_LAND_MODIS = 0.12, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_OCN_MODIS = 0.12, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_SNOW_MODIS = 0.12, &

                                         ULST_EMISS_CHN7_DIFF_THRESH_LAND_MTSAT = 0.12, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_OCN_MTSAT = 0.12, &
                                         ULST_EMISS_CHN7_DIFF_THRESH_SNOW_MTSAT = 0.12


 REAL, PRIVATE, SAVE :: &
                                         ULST_Emiss_Chn7_Diff_Thresh_Land, &
                                         ULST_Emiss_Chn7_Diff_Thresh_Ocean, &
                                         ULST_Emiss_Chn7_Diff_Thresh_Snow, &
                                         EMISS4_Emiss_Chn7_Desert_Thresh, &
                                         EMISS4_Emiss_Chn7_Ocn_Thresh, &
                                         EMISS4_Emiss_Chn7_Land_Thresh, &
                                         EMISS4_Emiss_Chn7_Snow_Thresh 

!----------------------------------------------------------------------
! other thresholds
!----------------------------------------------------------------------
 REAL, PRIVATE, PARAMETER :: REFL_CHN2_SNOW_THRESH = 10.0    !minimum albedo to allow snow to exist
 REAL, PRIVATE, PARAMETER :: BT_CHN14_SNOW_THRESH = 277.0   !maximum temperature to allow snow to exist
