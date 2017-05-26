!$Id: awg_cld_hght_clavrx_include_1.inc,v 1.2 2011/08/16 18:45:07 heidinger Exp $
!----------------------------------------------------------------------
!
! Name:
!   awg_cld_hght_alg.inc
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! Description:
!   Include file that declares a number of physical constants,
!   apriori estimates and uncertainty values
!
!----------------------------------------------------------------------

 !--- set starting index value
 INTEGER(KIND=INT4), PARAMETER, PRIVATE:: Element_Idx_Min = 1
 INTEGER(KIND=INT4), PARAMETER, PRIVATE:: Line_Idx_Min = 1
 INTEGER(KIND=INT4), PARAMETER, PRIVATE:: LRC_Meander_Flag = 1

 !---- set geocat flag
 INTEGER(KIND=INT4), PARAMETER :: Geocat_Flag = 0      !0 = not Geocat, 1 = Geocat

 !--- number of meta data flags
 INTEGER(KIND=INT4), PARAMETER :: NUM_META_DATA = 8

 !--- quality is determined by the ratio of estimated uncertainity relative
 !--- to that assume for the apriori (first guess).
 INTEGER, PRIVATE, PARAMETER :: CTH_PARAM_FAILED_RETREVIAL = 0
 INTEGER, PRIVATE, PARAMETER :: CTH_PARAM_LOW_QUALITY_RETREVIAL = 1  !retrived uncert. > 2/3 of apriori uncert.
 INTEGER, PRIVATE, PARAMETER :: CTH_PARAM_2_3_APRIORI_RETREVIAL = 2  !retrieved uncert. < 2/3 of apriori uncert.
 INTEGER, PRIVATE, PARAMETER :: CTH_PARAM_1_3_APRIORI_RETREVIAL = 3  !retrieved uncert. < 1/3 of apriori uncert.

 !--- LRC Computation Parameters
 INTEGER, PRIVATE, PARAMETER :: Use_Lrc_Flag = 1  ! (0 = sym%NO, 1 = sym%YES)
 INTEGER, PRIVATE, PARAMETER :: Grad_Flag_LRC = -1
 INTEGER, PRIVATE, PARAMETER :: Max_LRC_Distance = 10
 REAL, PRIVATE, PARAMETER :: Min_Bt_11um_LRC = 220.0
 REAL, PRIVATE, PARAMETER :: Max_Bt_11um_LRC = 290.0

 !---------------------------------------------------------------------
 ! PATMOS-x specific parameters
 !---------------------------------------------------------------------
 INTEGER(KIND=INT4), PRIVATE, PARAMETER :: Num_Levels_Rtm_Prof = 101
 INTEGER(KIND=INT4), PARAMETER :: Chan_Idx_67um = 27   !channel number for 6.7 micron
 INTEGER(KIND=INT4), PARAMETER :: Chan_Idx_73um = 28   !channel number for 7.3 micron
 INTEGER(KIND=INT4), PARAMETER :: Chan_Idx_85um = 29   !channel number for 8.5 micron
 INTEGER(KIND=INT4), PARAMETER :: Chan_Idx_11um = 31   !channel number for 11 micron
 INTEGER(KIND=INT4), PARAMETER :: Chan_Idx_12um = 32   !channel number for 12 micron
 INTEGER(KIND=INT4), PARAMETER :: Chan_Idx_133um = 33  !channel number for 13.3 micron

 !--- these parameters control the size of matrices
 INTEGER(KIND=INT4), PARAMETER :: Num_Param = 3              !number of retrieved parameters

 !--- maximum number of iterations allowed
 INTEGER(KIND=INT4), PARAMETER :: Iter_Max = 10              !maximum number of iterations

 !--- parameters that control the bottom-up estimation of Zc and Pc
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Lapse_Rate_Ocean = 9.8 !dry lapse rate over ocean for height K/km
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Min_P_inversion = 600.0 !min pressure to look for low-level inversion
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Delta_Psfc_Inversion = 50.0  !pressure above sfc to look for low-level inversion
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Delta_T_Layer_Inversion = 0.0   !maximum temperature difference to define an inversion layer

 !--- Bias correction values (sensor specific)
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: GAMMA_67um = 1.0     
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: GAMMA_73um = 1.0     
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: GAMMA_85um = 1.0     
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: GAMMA_11um = 1.0     
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: GAMMA_12um = 1.0     
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: GAMMA_133um = 1.0     

 !--- the parameters that provide the apriori values and their uncertainties
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tc_Ap_tropo_Offset_Cirrus = 15.0      !apriori Tc for opaque clouds
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tc_Ap_Uncer_Opaque  = 10.0            !apriori uncertainty of Tc for opaque clouds
 !REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tc_Ap_Uncer_Cirrus = 40.0  !!           !apriori uncertainty of Tc for cirrus
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tc_Ap_Uncer_Cirrus = 20.0  !!           !apriori uncertainty of Tc for cirrus
                                                                    ! Originally it was 20 in the baseline version of ACHA

 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Ec_Ap_Uncer_Opaque = 0.2             !apriori uncertainty of ec for opaque clouds
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Ec_Ap_Uncer_Cirrus = 0.4

 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Beta_Ap_Water = 1.3
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Beta_Ap_Uncer_Water = 0.2            !apriori uncertainty of  beta for ice
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Beta_Ap_Ice = 1.06
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Beta_Ap_Uncer_Ice = 0.2              !apriori uncertainty of  beta for water

 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tau_Ap_Fog_Type = 1.2               !apriori estimate of tau for fog cloud type
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tau_Ap_Water_Type = 2.3             !apriori estimate of tau for water cloud type
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tau_Ap_Mixed_Type = 2.3             !apriori estimate of tau for mixed cloud type
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tau_Ap_Opaque_Ice_Type = 2.3        !apriori estimate of tau for opaque ice cloud type
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tau_Ap_Cirrus_Type = 0.9            !apriori estimate of tau for cirrus cloud type
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Tau_Ap_Overlap_Type = 2.0           !apriori estimate of tau for multilayer cloud type

!--- specify calibration  errors
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_Cal_Uncer = 1.0
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_12um_Cal_Uncer = 0.5
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_12um_Cal_Uncer = 1.0  !BASELINE V5
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_133um_Cal_Uncer = 1.0
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_133um_Cal_Uncer = 2.0 !BASELINE V5
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_85um_Cal_Uncer = 0.5
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_73um_Cal_Uncer = 1.0
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_67um_Cal_Uncer = 1.0

!--- specify clear-sky model errors
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_Clr_Uncer_Water = 1.5
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_133um_Clr_Uncer_Water = 0.5
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_133um_Clr_Uncer_Water = 4.0 !BASELINE V5
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_12um_Clr_Uncer_Water = 0.5
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_12um_Clr_Uncer_Water = 1.0  !BASELINE V5
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_85um_Clr_Uncer_Water = 0.5
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_73um_Clr_Uncer_Water = 0.5
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_67um_Clr_Uncer_Water = 0.5

 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_Clr_Uncer_Land = 5.0
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_133um_Clr_Uncer_Land = 1.0
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_133um_Clr_Uncer_Land = 4.0 !BASELINE V5
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_12um_Clr_Uncer_Land = 1.0
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_85um_Clr_Uncer_Land = 1.0
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_73um_Clr_Uncer_Land = 1.0
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: T11um_67um_Clr_Uncer_Land = 1.0

!--- these parameters constrain the allowable solutions
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Min_Allowable_Tc = 160.0  !K
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: Pc_Lower_Cloud_Offset = 200.0  !hpa

!--- water clouds derived from Mie theory by A. Heidinger
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_11um_133um_FIT_WATER = -0.728113
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_11um_133um_FIT_WATER = 1.743389
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_11um_133um_FIT_WATER = -0.217447 !BASELINE V5
! REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_11um_133um_FIT_WATER = 1.250322  !BASELINE V5 
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_11um_85um_FIT_WATER = 0.930569
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_11um_85um_FIT_WATER = 0.048857
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_11um_73um_FIT_WATER = 0.532266
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_11um_73um_FIT_WATER = 0.438829
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_11um_67um_FIT_WATER = 0.268115
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_11um_67um_FIT_WATER = 0.702683

!--- these values are for ice clouds from Ping Yang's Database (10-100um)

 !- aggregates
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_11um_133um_FIT_ICE = -0.325909
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_11um_133um_FIT_ICE = 1.340083
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_11um_85um_FIT_ICE = 2.471422
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_11um_85um_FIT_ICE = -1.492580
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_11um_73um_FIT_ICE = 1.755726
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_11um_73um_FIT_ICE = -0.780324
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: A_BETA_11um_67um_FIT_ICE = 1.617212
 REAL(KIND=REAL4), PARAMETER, PRIVATE:: B_BETA_11um_67um_FIT_ICE = -0.646041


 REAL(KIND=REAL4), PARAMETER, PRIVATE, DIMENSION(NUM_PARAM):: delta_x_max = (/20.0,0.1,0.1/)

!------------------------------------------------------------------------
! Clear-sky Error Covariance Terms
! These come from two days of AQUA/MODIS (Jan and July 2009) with PFAAST
!------------------------------------------------------------------------

!--- Water
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Mean_Water = 237.175
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Mean_Water = 285.360
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Mean_Water = 287.384
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Mean_Water = 285.908
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Mean_Water = 265.054
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Mean_Land = 236.521
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Mean_Land = 285.905
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Mean_Land = 287.907
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Mean_Land = 286.094
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Mean_Land = 265.588
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Mean_Desert = 205.415
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Mean_Desert = 226.472
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Mean_Desert = 228.628
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Mean_Desert = 228.141
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Mean_Desert = 218.694
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Mean_Snow = 216.762
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Mean_Snow = 247.841
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Mean_Snow = 252.854
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Mean_Snow = 252.213
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Mean_Snow = 235.644
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Mean_Arctic = 223.444
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Mean_Arctic = 261.708
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Mean_Arctic = 268.111
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Mean_Arctic = 266.699
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Mean_Arctic = 247.740
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Mean_Antarctic = 216.762
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Mean_Antarctic = 247.841
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Mean_Antarctic = 252.854
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Mean_Antarctic = 252.213
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Mean_Antarctic = 235.644

REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_67um_Covar_Water =  13.865
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_73um_Covar_Water =   6.620
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_85um_Covar_Water =   2.656
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_11um_Covar_Water =   1.603
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_12um_Covar_Water =   1.830
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Bt_133um_Covar_Water =   0.818
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_12um_Covar_Water =   1.679
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_133um_Covar_Water =   0.813
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_85um_Covar_Water =   1.450
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_67um_Covar_Water =   0.831
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_73um_Covar_Water =   0.834
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_133um_Covar_Water =   0.732
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_133um_Covar_Water =   0.956
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_85um_Covar_Water =   1.493
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_67um_Covar_Water =   1.035
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_133um_Covar_Water =   1.129
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_133um_Covar_Water =   1.364
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_12um_Covar_Water =  -0.076
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_67um_Covar_Water =   0.772
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_73um_Covar_Water =   0.769
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_85um_Covar_Water =   0.154
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_133um_Covar_Water =   0.791
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_67um_Covar_Water =  13.806
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_73um_Covar_Water =   6.556
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_85um_Btd_11um_85um_Covar_Water =   1.360
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_12um_Covar_Water =   0.075
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_133um_Btd_11um_133um_Covar_Water =   0.796
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_133um_Covar_Water =   0.067
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_85um_Covar_Water =  -0.033
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_67um_Covar_Water =   0.127
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_73um_Covar_Water =   0.171
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_133um_Covar_Water =   1.088
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_133um_Covar_Water =   1.321
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_67um_Covar_Land =  10.979
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_73um_Covar_Land =   6.422
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_85um_Covar_Land =  29.028
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_11um_Covar_Land =  34.016
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_12um_Covar_Land =  31.408
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Bt_133um_Covar_Land =   4.869
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_12um_Covar_Land =  32.499
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_133um_Covar_Land =  10.784
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_85um_Covar_Land =  30.984
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_67um_Covar_Land =   0.504
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_73um_Covar_Land =   3.026
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_133um_Covar_Land =   9.772
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_133um_Covar_Land =  10.814
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_85um_Covar_Land =  29.487
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_67um_Covar_Land =   0.807
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_133um_Covar_Land =   1.207
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_133um_Covar_Land =   2.741
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_12um_Covar_Land =   1.518
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_67um_Covar_Land =  33.514
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_73um_Covar_Land =  30.991
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_85um_Covar_Land =   3.033
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_133um_Covar_Land =  23.233
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_67um_Covar_Land =  43.987
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_73um_Covar_Land =  34.386
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_85um_Btd_11um_85um_Covar_Land =   1.078
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_12um_Covar_Land =   0.427
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_133um_Btd_11um_133um_Covar_Land =  17.318
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_133um_Covar_Land =   1.548
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_85um_Covar_Land =   0.021
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_67um_Covar_Land =   1.821
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_73um_Covar_Land =   1.932
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_133um_Covar_Land =  23.936
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_133um_Covar_Land =  22.948
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_67um_Covar_Desert =  12.948
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_73um_Covar_Desert =   7.989
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_85um_Covar_Desert =  44.966
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_11um_Covar_Desert =  53.745
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_12um_Covar_Desert =  50.583
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Bt_133um_Covar_Desert =   6.422
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_12um_Covar_Desert =  51.851
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_133um_Covar_Desert =  15.636
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_85um_Covar_Desert =  44.917
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_67um_Covar_Desert =  -0.334
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_73um_Covar_Desert =   3.765
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_133um_Covar_Desert =  13.737
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_133um_Covar_Desert =  15.616
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_85um_Covar_Desert =  42.656
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_67um_Covar_Desert =   0.024
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_133um_Covar_Desert =   1.045
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_133um_Covar_Desert =   3.315
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_12um_Covar_Desert =   1.895
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_67um_Covar_Desert =  54.080
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_73um_Covar_Desert =  49.981
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_85um_Covar_Desert =   8.829
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_133um_Covar_Desert =  38.110
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_67um_Covar_Desert =  67.362
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_73um_Covar_Desert =  54.205
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_85um_Btd_11um_85um_Covar_Desert =   8.880
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_12um_Covar_Desert =   0.629
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_133um_Btd_11um_133um_Covar_Desert =  28.896
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_133um_Covar_Desert =   1.875
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_85um_Covar_Desert =  -0.365
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_67um_Covar_Desert =   2.254
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_73um_Covar_Desert =   2.427
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_133um_Covar_Desert =  39.489
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_133um_Covar_Desert =  37.660
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_67um_Covar_Snow =   7.450
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_73um_Covar_Snow =  10.772
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_85um_Covar_Snow =  33.423
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_11um_Covar_Snow =  32.312
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_12um_Covar_Snow =  32.714
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Bt_133um_Covar_Snow =   9.407
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_12um_Covar_Snow =  32.484
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_133um_Covar_Snow =  16.957
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_85um_Covar_Snow =  31.103
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_67um_Covar_Snow =   1.016
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_73um_Covar_Snow =  12.643
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_133um_Covar_Snow =  16.397
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_133um_Covar_Snow =  17.089
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_85um_Covar_Snow =  31.286
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_67um_Covar_Snow =   1.070
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_133um_Covar_Snow =   1.234
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_133um_Covar_Snow =   7.545
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_12um_Covar_Snow =  -0.172
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_67um_Covar_Snow =  31.297
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_73um_Covar_Snow =  19.669
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_85um_Covar_Snow =   1.209
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_133um_Covar_Snow =  15.355
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_67um_Covar_Snow =  37.732
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_73um_Covar_Snow =  17.798
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_85um_Btd_11um_85um_Covar_Snow =   3.530
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_12um_Covar_Snow =   0.057
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_133um_Btd_11um_133um_Covar_Snow =   7.804
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_133um_Covar_Snow =  -0.041
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_85um_Covar_Snow =   0.011
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_67um_Covar_Snow =  -0.117
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_73um_Covar_Snow =   0.079
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_133um_Covar_Snow =  15.574
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_133um_Covar_Snow =  10.257
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_67um_Covar_Arctic =   2.916
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_73um_Covar_Arctic =   2.322
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_85um_Covar_Arctic =  18.341
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_11um_Covar_Arctic =  20.321
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_12um_Covar_Arctic =  19.717
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Bt_133um_Covar_Arctic =   4.462
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_12um_Covar_Arctic =  19.993
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_133um_Covar_Arctic =   9.062
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_85um_Covar_Arctic =  19.159
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_67um_Covar_Arctic =   0.328
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_73um_Covar_Arctic =   3.610
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_133um_Covar_Arctic =   8.650
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_133um_Covar_Arctic =   8.967
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_85um_Covar_Arctic =  18.848
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_67um_Covar_Arctic =   0.346
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_133um_Covar_Arctic =   0.469
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_133um_Covar_Arctic =   2.156
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_12um_Covar_Arctic =   0.328
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_67um_Covar_Arctic =  19.993
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_73um_Covar_Arctic =  16.711
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_85um_Covar_Arctic =   1.162
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_133um_Covar_Arctic =  11.259
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_67um_Covar_Arctic =  22.581
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_73um_Covar_Arctic =  15.424
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_85um_Btd_11um_85um_Covar_Arctic =   0.344
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_12um_Covar_Arctic =   0.052
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_133um_Btd_11um_133um_Covar_Arctic =   6.658
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_133um_Covar_Arctic =   0.233
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_85um_Covar_Arctic =   0.017
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_67um_Covar_Arctic =   0.346
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_73um_Covar_Arctic =   0.319
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_133um_Covar_Arctic =  11.400
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_133um_Covar_Arctic =   9.805
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_67um_Covar_Antarctic =   7.823
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_73um_Covar_Antarctic =  13.391
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_85um_Covar_Antarctic =  34.741
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_11um_Covar_Antarctic =  34.387
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_12um_Covar_Antarctic =  33.899
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_133um_Bt_133um_Covar_Antarctic =  17.091
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_12um_Covar_Antarctic =  34.128
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_133um_Covar_Antarctic =  23.691
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_85um_Covar_Antarctic =  34.512
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_67um_Covar_Antarctic =  10.187
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Bt_73um_Covar_Antarctic =  20.088
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_85um_Bt_133um_Covar_Antarctic =  23.931
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_133um_Covar_Antarctic =  23.518
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_85um_Covar_Antarctic =  34.256
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_12um_Bt_67um_Covar_Antarctic =  10.159
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_67um_Bt_133um_Covar_Antarctic =   7.746
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_73um_Bt_133um_Covar_Antarctic =  14.648
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_12um_Covar_Antarctic =   0.259
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_67um_Covar_Antarctic =  24.200
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_73um_Covar_Antarctic =  14.299
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_85um_Covar_Antarctic =  -0.125
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Bt_11um_Btd_11um_133um_Covar_Antarctic =  10.696
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_67um_Covar_Antarctic =  21.836
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_73um_Covar_Antarctic =   7.602
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_85um_Btd_11um_85um_Covar_Antarctic =   0.104
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_12um_Covar_Antarctic =   0.030
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_133um_Btd_11um_133um_Covar_Antarctic =   4.096
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_133um_Covar_Antarctic =   0.086
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_85um_Covar_Antarctic =   0.004
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_67um_Covar_Antarctic =   0.231
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_12um_Btd_11um_73um_Covar_Antarctic =   0.128
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_67um_Btd_11um_133um_Covar_Antarctic =   8.255
REAL(KIND=REAL4),PARAMETER,PRIVATE:: Btd_11um_73um_Btd_11um_133um_Covar_Antarctic =   5.255
