!$Id$
!----------------------------------------------------------------------
!
! Name:
!   awg_cld_hght_include_1.inc
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! Description:
!   Include file that declares a number of physical constants,
!   apriori estimates and uncertainty values
!
!----------------------------------------------------------------------

 !--- set starting index value
 integer(kind=INT4), parameter, private:: Element_Idx_Min = 1
 integer(kind=INT4), parameter, private:: Line_Idx_Min = 1
 integer(kind=INT4), parameter, private:: LRC_Meander_Flag = 1

 !--- number of meta data flags
 integer(kind=INT4), parameter :: NUM_META_DATA = 8

 !--- quality is determined by the ratio of estimated uncertainity relative
 !--- to that assume for the apriori (first guess).
 integer, private, parameter :: CTH_PARAM_FAILED_RETREVIAL = 0
 integer, private, parameter :: CTH_PARAM_LOW_QUALITY_RETREVIAL = 1  !retrived uncert. > 2/3 of apriori uncert.
 integer, private, parameter :: CTH_PARAM_2_3_APRIORI_RETREVIAL = 2  !retrieved uncert. < 2/3 of apriori uncert.
 integer, private, parameter :: CTH_PARAM_1_3_APRIORI_RETREVIAL = 3  !retrieved uncert. < 1/3 of apriori uncert.

 !--- LRC Computation Parameters
 integer, private, parameter :: Use_Lrc_Flag = 1  ! (0 = sym%NO, 1 = sym%YES)
 integer, private, parameter :: Grad_Flag_LRC = -1
 integer, private, parameter :: Max_LRC_Distance = 10
 real, private, parameter :: Min_LRC_Jump = 0.5
 real, private, parameter :: Max_LRC_Jump = 10.0
 real, private, parameter :: Min_Bt_11um_LRC = 220.0
 real, private, parameter :: Max_Bt_11um_LRC = 290.0

 !--- Cirrus Box Parameters
 integer, private, parameter :: USE_CIRRUS_FLAG = 1  ! (0 = sym%NO, 1 = sym%YES)
 real (kind=real4), parameter:: EMISSIVITY_MIN_TEMPERATURE_CIRRUS = 0.7
 integer (kind=int4), parameter:: COUNT_MIN_TEMPERATURE_CIRRUS = 5
 integer (kind=int4), parameter:: BOX_WIDTH_KM = 200

 !---------------------------------------------------------------------
 ! retrieval specific parameters
 !---------------------------------------------------------------------
 integer(kind=INT4), private, parameter :: Num_Levels_Rtm_Prof = 101

 !--- these parameters control the size of matrices
 integer(kind=INT4), parameter :: Num_Param = 3              !number of retrieved parameters

 !--- maximum number of iterations allowed
 integer(kind=INT4), parameter :: Iter_Idx_Max = 10              !maximum number of iterations

 !--- parameters that control the bottom-up estimation of Zc and Pc
 real(kind=real4), parameter, private:: ABS_LAPSE_RATE_UNCER_MIN = 0.0005 !K/m smallest lapse rate over height uncertainty calc
 real(kind=real4), parameter, private:: ABS_LAPSE_RATE_DP_DZ_UNCER_MIN = 0.00005 !K/hPa smallest lapse rate over height uncertainty calc

 real(kind=real4), parameter, private:: Min_P_inversion = 600.0 !min pressure to look for low-level inversion
 real(kind=real4), parameter, private:: Max_Delta_T_inversion = 15.0 !max temperature difference (surface -cloud) to look for low-level inversion
 real(kind=real4), parameter, private:: Delta_Psfc_Inversion = 50.0  !pressure above sfc to look for low-level inversion
 real(kind=real4), parameter, private:: Delta_T_Layer_Inversion = 0.0   !maximum temperature difference to define an inversion layer
 integer(kind=INT1), parameter, private:: ALLOW_STRATOSPHERE_SOLUTION_FLAG = 1  !allow interpolation in stratosphere

 !--- the parameters that provide the apriori values and their uncertainties
 real(kind=real4), parameter, private:: Tc_Ap_tropo_Offset_Cirrus = 15.0      !apriori Tc for opaque clouds

 real(kind=real4), parameter, private:: Tc_Ap_Uncer_Opaque  = 10.0            !apriori uncertainty of Tc for opaque clouds

!real(kind=real4), parameter, private:: Tc_Ap_Uncer_Cirrus = 40.0  !!           !apriori uncertainty of Tc for cirrus

 real(kind=real4), parameter, private:: Tc_Ap_Uncer_Cirrus = 20.0  !!           !apriori uncertainty of Tc for cirrus
                                                                    ! Originally it was 20 in the baseline version of ACHA

!real(kind=real4), parameter, private:: Ec_Ap_Uncer_Opaque = 0.2             !apriori uncertainty of ec for opaque clouds
 real(kind=real4), parameter, private:: Ec_Ap_Uncer_Opaque = 0.1             !apriori uncertainty of ec for opaque clouds
 real(kind=real4), parameter, private:: Ec_Ap_Uncer_Cirrus = 0.4

 real(kind=real4), parameter, private:: Beta_Ap_Water = 1.3
 real(kind=real4), parameter, private:: Beta_Ap_Uncer_Water = 0.2            !apriori uncertainty of  beta for ice
 real(kind=real4), parameter, private:: Beta_Ap_Ice = 1.06
 real(kind=real4), parameter, private:: Beta_Ap_Uncer_Ice = 0.2              !apriori uncertainty of  beta for water

 real(kind=real4), parameter, private:: Tau_Ap_Water_Phase= 3.0             !apriori estimate of tau for water phase clod
 real(kind=real4), parameter, private:: Tau_Ap_Fog_Type = 1.2               !apriori estimate of tau for fog cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Water_Type = 2.3             !apriori estimate of tau for water cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Mixed_Type = 2.3             !apriori estimate of tau for mixed cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Opaque_Ice_Type = 2.3        !apriori estimate of tau for opaque ice cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Cirrus_Type = 0.9            !apriori estimate of tau for cirrus cloud type
 real(kind=real4), parameter, private:: Tau_Ap_Overlap_Type = 2.0           !apriori estimate of tau for multilayer cloud type

 !--- specify calibration  errors
 real(kind=real4), parameter, private:: T11um_Cal_Uncer = 1.0
 real(kind=real4), parameter, private:: T11um_12um_Cal_Uncer = 0.5
!--> real(kind=real4), parameter, private:: T11um_12um_Cal_Uncer = 1.0  !BASELINE V5
 real(kind=real4), parameter, private:: T11um_133um_Cal_Uncer = 1.0
!--> real(kind=real4), parameter, private:: T11um_133um_Cal_Uncer = 2.0 !BASELINE V5
 real(kind=real4), parameter, private:: T11um_85um_Cal_Uncer = 0.5
 real(kind=real4), parameter, private:: T11um_67um_Cal_Uncer = 1.0

 !--- specify clear-sky model errors
 real(kind=real4), parameter, private:: T11um_Clr_Uncer_Water = 1.5
 real(kind=real4), parameter, private:: T11um_133um_Clr_Uncer_Water = 0.5
!--> real(kind=real4), parameter, private:: T11um_133um_Clr_Uncer_Water = 4.0 !BASELINE V5
 real(kind=real4), parameter, private:: T11um_12um_Clr_Uncer_Water = 0.5
!--> real(kind=real4), parameter, private:: T11um_12um_Clr_Uncer_Water = 1.0  !BASELINE V5
 real(kind=real4), parameter, private:: T11um_85um_Clr_Uncer_Water = 0.5
 real(kind=real4), parameter, private:: T11um_67um_Clr_Uncer_Water = 0.5

 real(kind=real4), parameter, private:: T11um_Clr_Uncer_Land = 5.0
 real(kind=real4), parameter, private:: T11um_133um_Clr_Uncer_Land = 1.0
!--> real(kind=real4), parameter, private:: T11um_133um_Clr_Uncer_Land = 4.0 !BASELINE V5
 real(kind=real4), parameter, private:: T11um_12um_Clr_Uncer_Land = 1.0
 real(kind=real4), parameter, private:: T11um_85um_Clr_Uncer_Land = 1.0
 real(kind=real4), parameter, private:: T11um_67um_Clr_Uncer_Land = 1.0

 !--- these parameters constrain the allowable solutions
 real(kind=real4), parameter, private:: Min_Allowable_Tc = 180.0  !K
 real(kind=real4), parameter, private:: Pc_Lower_Cloud_Offset = 200.0  !hpa

 !--- water cloud microphysical model
 INCLUDE 'acha_water_cloud_microphysical_model.inc'

 !--- ice cloud microphysical model
 INCLUDE 'acha_ice_cloud_microphysical_model.inc'

 !--- other options to explore sensitivity to microphysics
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_aggcol_b005.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_aggcol_b01.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_aggcol_b02.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_droxtal_b01.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_hbrose_b01.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_holcol_b01.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_plates_b01.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_laraggplates_b01.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_smaggplates_b01.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_solcol_b01.inc'
!INCLUDE 'acha_ice_cloud_microphysical_model_sevrough_sbrose_b01.inc'

 !--- limits on steps taken in retrieval
 real(kind=real4), parameter, private, DIMENSION(NUM_PARAM):: delta_x_max = (/20.0,0.1,0.1/)


 !--- parameters that allow for cloud geometrical thickness estimation
 real(kind=real4), parameter, private:: WATER_EXTINCTION = 25.00   ! 1/km

 real(kind=real4), parameter, private:: ICE_EXTINCTION1 = 1.71
 real(kind=real4), parameter, private:: ICE_EXTINCTION2 = 1.87
 real(kind=real4), parameter, private:: ICE_EXTINCTION3 = 2.24
 real(kind=real4), parameter, private:: ICE_EXTINCTION4 = 2.88
 real(kind=real4), parameter, private:: ICE_EXTINCTION5 = 4.74

 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION1 = 0.13
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION2 = 0.25
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION3 = 0.39
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION4 = 0.55
 real(kind=real4), parameter, private:: CIRRUS_EXTINCTION5 = 0.67

! real(kind=real4), parameter, private:: ICE_EXTINCTION = 5.00   ! 1/km
! real(kind=real4), parameter, private:: CIRRUS_EXTINCTION = 0.15 !0.50   ! 1/km

!------------------------------------------------------------------------
! Clear-sky Error Covariance Terms
! These come from two days of AQUA/MODIS (Jan and July 2009) with PFAAST
!------------------------------------------------------------------------

!--- Water
real(kind=real4),parameter,private:: Bt_67um_Mean_Water = 237.175
real(kind=real4),parameter,private:: Bt_85um_Mean_Water = 285.360
real(kind=real4),parameter,private:: Bt_11um_Mean_Water = 287.384
real(kind=real4),parameter,private:: Bt_12um_Mean_Water = 285.908
real(kind=real4),parameter,private:: Bt_133um_Mean_Water = 265.054
real(kind=real4),parameter,private:: Bt_67um_Mean_Land = 236.521
real(kind=real4),parameter,private:: Bt_85um_Mean_Land = 285.905
real(kind=real4),parameter,private:: Bt_11um_Mean_Land = 287.907
real(kind=real4),parameter,private:: Bt_12um_Mean_Land = 286.094
real(kind=real4),parameter,private:: Bt_133um_Mean_Land = 265.588
real(kind=real4),parameter,private:: Bt_67um_Mean_Desert = 205.415
real(kind=real4),parameter,private:: Bt_85um_Mean_Desert = 226.472
real(kind=real4),parameter,private:: Bt_11um_Mean_Desert = 228.628
real(kind=real4),parameter,private:: Bt_12um_Mean_Desert = 228.141
real(kind=real4),parameter,private:: Bt_133um_Mean_Desert = 218.694
real(kind=real4),parameter,private:: Bt_67um_Mean_Snow = 216.762
real(kind=real4),parameter,private:: Bt_85um_Mean_Snow = 247.841
real(kind=real4),parameter,private:: Bt_11um_Mean_Snow = 252.854
real(kind=real4),parameter,private:: Bt_12um_Mean_Snow = 252.213
real(kind=real4),parameter,private:: Bt_133um_Mean_Snow = 235.644
real(kind=real4),parameter,private:: Bt_67um_Mean_Arctic = 223.444
real(kind=real4),parameter,private:: Bt_85um_Mean_Arctic = 261.708
real(kind=real4),parameter,private:: Bt_11um_Mean_Arctic = 268.111
real(kind=real4),parameter,private:: Bt_12um_Mean_Arctic = 266.699
real(kind=real4),parameter,private:: Bt_133um_Mean_Arctic = 247.740
real(kind=real4),parameter,private:: Bt_67um_Mean_Antarctic = 216.762
real(kind=real4),parameter,private:: Bt_85um_Mean_Antarctic = 247.841
real(kind=real4),parameter,private:: Bt_11um_Mean_Antarctic = 252.854
real(kind=real4),parameter,private:: Bt_12um_Mean_Antarctic = 252.213
real(kind=real4),parameter,private:: Bt_133um_Mean_Antarctic = 235.644

real(kind=real4),parameter,private:: Bt_67um_Bt_67um_Covar_Water =  13.865
real(kind=real4),parameter,private:: Bt_85um_Bt_85um_Covar_Water =   2.656
real(kind=real4),parameter,private:: Bt_11um_Bt_11um_Covar_Water =   1.603
real(kind=real4),parameter,private:: Bt_12um_Bt_12um_Covar_Water =   1.830
real(kind=real4),parameter,private:: Bt_133um_Bt_133um_Covar_Water =   0.818
real(kind=real4),parameter,private:: Bt_11um_Bt_12um_Covar_Water =   1.679
real(kind=real4),parameter,private:: Bt_11um_Bt_133um_Covar_Water =   0.813
real(kind=real4),parameter,private:: Bt_11um_Bt_85um_Covar_Water =   1.450
real(kind=real4),parameter,private:: Bt_11um_Bt_67um_Covar_Water =   0.831
real(kind=real4),parameter,private:: Bt_85um_Bt_133um_Covar_Water =   0.732
real(kind=real4),parameter,private:: Bt_12um_Bt_133um_Covar_Water =   0.956
real(kind=real4),parameter,private:: Bt_12um_Bt_85um_Covar_Water =   1.493
real(kind=real4),parameter,private:: Bt_12um_Bt_67um_Covar_Water =   1.035
real(kind=real4),parameter,private:: Bt_67um_Bt_133um_Covar_Water =   1.129
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_12um_Covar_Water =  -0.076
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_67um_Covar_Water =   0.772
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_85um_Covar_Water =   0.154
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_133um_Covar_Water =   0.791
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_67um_Covar_Water =  13.806
real(kind=real4),parameter,private:: Btd_11um_85um_Btd_11um_85um_Covar_Water =   1.360
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_12um_Covar_Water =   0.075
real(kind=real4),parameter,private:: Btd_11um_133um_Btd_11um_133um_Covar_Water =   0.796
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_133um_Covar_Water =   0.067
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_85um_Covar_Water =  -0.033
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_67um_Covar_Water =   0.127
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_133um_Covar_Water =   1.088
real(kind=real4),parameter,private:: Bt_67um_Bt_67um_Covar_Land =  10.979
real(kind=real4),parameter,private:: Bt_85um_Bt_85um_Covar_Land =  29.028
real(kind=real4),parameter,private:: Bt_11um_Bt_11um_Covar_Land =  34.016
real(kind=real4),parameter,private:: Bt_12um_Bt_12um_Covar_Land =  31.408
real(kind=real4),parameter,private:: Bt_133um_Bt_133um_Covar_Land =   4.869
real(kind=real4),parameter,private:: Bt_11um_Bt_12um_Covar_Land =  32.499
real(kind=real4),parameter,private:: Bt_11um_Bt_133um_Covar_Land =  10.784
real(kind=real4),parameter,private:: Bt_11um_Bt_85um_Covar_Land =  30.984
real(kind=real4),parameter,private:: Bt_11um_Bt_67um_Covar_Land =   0.504
real(kind=real4),parameter,private:: Bt_85um_Bt_133um_Covar_Land =   9.772
real(kind=real4),parameter,private:: Bt_12um_Bt_133um_Covar_Land =  10.814
real(kind=real4),parameter,private:: Bt_12um_Bt_85um_Covar_Land =  29.487
real(kind=real4),parameter,private:: Bt_12um_Bt_67um_Covar_Land =   0.807
real(kind=real4),parameter,private:: Bt_67um_Bt_133um_Covar_Land =   1.207
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_12um_Covar_Land =   1.518
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_67um_Covar_Land =  33.514
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_85um_Covar_Land =   3.033
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_133um_Covar_Land =  23.233
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_67um_Covar_Land =  43.987
real(kind=real4),parameter,private:: Btd_11um_85um_Btd_11um_85um_Covar_Land =   1.078
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_12um_Covar_Land =   0.427
real(kind=real4),parameter,private:: Btd_11um_133um_Btd_11um_133um_Covar_Land =  17.318
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_133um_Covar_Land =   1.548
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_85um_Covar_Land =   0.021
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_67um_Covar_Land =   1.821
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_133um_Covar_Land =  23.936
real(kind=real4),parameter,private:: Bt_67um_Bt_67um_Covar_Desert =  12.948
real(kind=real4),parameter,private:: Bt_85um_Bt_85um_Covar_Desert =  44.966
real(kind=real4),parameter,private:: Bt_11um_Bt_11um_Covar_Desert =  53.745
real(kind=real4),parameter,private:: Bt_12um_Bt_12um_Covar_Desert =  50.583
real(kind=real4),parameter,private:: Bt_133um_Bt_133um_Covar_Desert =   6.422
real(kind=real4),parameter,private:: Bt_11um_Bt_12um_Covar_Desert =  51.851
real(kind=real4),parameter,private:: Bt_11um_Bt_133um_Covar_Desert =  15.636
real(kind=real4),parameter,private:: Bt_11um_Bt_85um_Covar_Desert =  44.917
real(kind=real4),parameter,private:: Bt_11um_Bt_67um_Covar_Desert =  -0.334
real(kind=real4),parameter,private:: Bt_85um_Bt_133um_Covar_Desert =  13.737
real(kind=real4),parameter,private:: Bt_12um_Bt_133um_Covar_Desert =  15.616
real(kind=real4),parameter,private:: Bt_12um_Bt_85um_Covar_Desert =  42.656
real(kind=real4),parameter,private:: Bt_12um_Bt_67um_Covar_Desert =   0.024
real(kind=real4),parameter,private:: Bt_67um_Bt_133um_Covar_Desert =   1.045
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_12um_Covar_Desert =   1.895
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_67um_Covar_Desert =  54.080
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_85um_Covar_Desert =   8.829
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_133um_Covar_Desert =  38.110
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_67um_Covar_Desert =  67.362
real(kind=real4),parameter,private:: Btd_11um_85um_Btd_11um_85um_Covar_Desert =   8.880
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_12um_Covar_Desert =   0.629
real(kind=real4),parameter,private:: Btd_11um_133um_Btd_11um_133um_Covar_Desert =  28.896
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_133um_Covar_Desert =   1.875
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_85um_Covar_Desert =  -0.365
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_67um_Covar_Desert =   2.254
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_133um_Covar_Desert =  39.489
real(kind=real4),parameter,private:: Bt_67um_Bt_67um_Covar_Snow =   7.450
real(kind=real4),parameter,private:: Bt_85um_Bt_85um_Covar_Snow =  33.423
real(kind=real4),parameter,private:: Bt_11um_Bt_11um_Covar_Snow =  32.312
real(kind=real4),parameter,private:: Bt_12um_Bt_12um_Covar_Snow =  32.714
real(kind=real4),parameter,private:: Bt_133um_Bt_133um_Covar_Snow =   9.407
real(kind=real4),parameter,private:: Bt_11um_Bt_12um_Covar_Snow =  32.484
real(kind=real4),parameter,private:: Bt_11um_Bt_133um_Covar_Snow =  16.957
real(kind=real4),parameter,private:: Bt_11um_Bt_85um_Covar_Snow =  31.103
real(kind=real4),parameter,private:: Bt_11um_Bt_67um_Covar_Snow =   1.016
real(kind=real4),parameter,private:: Bt_85um_Bt_133um_Covar_Snow =  16.397
real(kind=real4),parameter,private:: Bt_12um_Bt_133um_Covar_Snow =  17.089
real(kind=real4),parameter,private:: Bt_12um_Bt_85um_Covar_Snow =  31.286
real(kind=real4),parameter,private:: Bt_12um_Bt_67um_Covar_Snow =   1.070
real(kind=real4),parameter,private:: Bt_67um_Bt_133um_Covar_Snow =   1.234
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_12um_Covar_Snow =  -0.172
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_67um_Covar_Snow =  31.297
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_85um_Covar_Snow =   1.209
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_133um_Covar_Snow =  15.355
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_67um_Covar_Snow =  37.732
real(kind=real4),parameter,private:: Btd_11um_85um_Btd_11um_85um_Covar_Snow =   3.530
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_12um_Covar_Snow =   0.057
real(kind=real4),parameter,private:: Btd_11um_133um_Btd_11um_133um_Covar_Snow =   7.804
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_133um_Covar_Snow =  -0.041
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_85um_Covar_Snow =   0.011
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_67um_Covar_Snow =  -0.117
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_133um_Covar_Snow =  15.574
real(kind=real4),parameter,private:: Bt_67um_Bt_67um_Covar_Arctic =   2.916
real(kind=real4),parameter,private:: Bt_85um_Bt_85um_Covar_Arctic =  18.341
real(kind=real4),parameter,private:: Bt_11um_Bt_11um_Covar_Arctic =  20.321
real(kind=real4),parameter,private:: Bt_12um_Bt_12um_Covar_Arctic =  19.717
real(kind=real4),parameter,private:: Bt_133um_Bt_133um_Covar_Arctic =   4.462
real(kind=real4),parameter,private:: Bt_11um_Bt_12um_Covar_Arctic =  19.993
real(kind=real4),parameter,private:: Bt_11um_Bt_133um_Covar_Arctic =   9.062
real(kind=real4),parameter,private:: Bt_11um_Bt_85um_Covar_Arctic =  19.159
real(kind=real4),parameter,private:: Bt_11um_Bt_67um_Covar_Arctic =   0.328
real(kind=real4),parameter,private:: Bt_85um_Bt_133um_Covar_Arctic =   8.650
real(kind=real4),parameter,private:: Bt_12um_Bt_133um_Covar_Arctic =   8.967
real(kind=real4),parameter,private:: Bt_12um_Bt_85um_Covar_Arctic =  18.848
real(kind=real4),parameter,private:: Bt_12um_Bt_67um_Covar_Arctic =   0.346
real(kind=real4),parameter,private:: Bt_67um_Bt_133um_Covar_Arctic =   0.469
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_12um_Covar_Arctic =   0.328
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_67um_Covar_Arctic =  19.993
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_85um_Covar_Arctic =   1.162
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_133um_Covar_Arctic =  11.259
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_67um_Covar_Arctic =  22.581
real(kind=real4),parameter,private:: Btd_11um_85um_Btd_11um_85um_Covar_Arctic =   0.344
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_12um_Covar_Arctic =   0.052
real(kind=real4),parameter,private:: Btd_11um_133um_Btd_11um_133um_Covar_Arctic =   6.658
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_133um_Covar_Arctic =   0.233
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_85um_Covar_Arctic =   0.017
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_67um_Covar_Arctic =   0.346
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_133um_Covar_Arctic =  11.400
real(kind=real4),parameter,private:: Bt_67um_Bt_67um_Covar_Antarctic =   7.823
real(kind=real4),parameter,private:: Bt_85um_Bt_85um_Covar_Antarctic =  34.741
real(kind=real4),parameter,private:: Bt_11um_Bt_11um_Covar_Antarctic =  34.387
real(kind=real4),parameter,private:: Bt_12um_Bt_12um_Covar_Antarctic =  33.899
real(kind=real4),parameter,private:: Bt_133um_Bt_133um_Covar_Antarctic =  17.091
real(kind=real4),parameter,private:: Bt_11um_Bt_12um_Covar_Antarctic =  34.128
real(kind=real4),parameter,private:: Bt_11um_Bt_133um_Covar_Antarctic =  23.691
real(kind=real4),parameter,private:: Bt_11um_Bt_85um_Covar_Antarctic =  34.512
real(kind=real4),parameter,private:: Bt_11um_Bt_67um_Covar_Antarctic =  10.187
real(kind=real4),parameter,private:: Bt_85um_Bt_133um_Covar_Antarctic =  23.931
real(kind=real4),parameter,private:: Bt_12um_Bt_133um_Covar_Antarctic =  23.518
real(kind=real4),parameter,private:: Bt_12um_Bt_85um_Covar_Antarctic =  34.256
real(kind=real4),parameter,private:: Bt_12um_Bt_67um_Covar_Antarctic =  10.159
real(kind=real4),parameter,private:: Bt_67um_Bt_133um_Covar_Antarctic =   7.746
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_12um_Covar_Antarctic =   0.259
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_67um_Covar_Antarctic =  24.200
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_85um_Covar_Antarctic =  -0.125
real(kind=real4),parameter,private:: Bt_11um_Btd_11um_133um_Covar_Antarctic =  10.696
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_67um_Covar_Antarctic =  21.836
real(kind=real4),parameter,private:: Btd_11um_85um_Btd_11um_85um_Covar_Antarctic =   0.104
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_12um_Covar_Antarctic =   0.030
real(kind=real4),parameter,private:: Btd_11um_133um_Btd_11um_133um_Covar_Antarctic =   4.096
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_133um_Covar_Antarctic =   0.086
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_85um_Covar_Antarctic =   0.004
real(kind=real4),parameter,private:: Btd_11um_12um_Btd_11um_67um_Covar_Antarctic =   0.231
real(kind=real4),parameter,private:: Btd_11um_67um_Btd_11um_133um_Covar_Antarctic =   8.255
