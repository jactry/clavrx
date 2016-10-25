!$Id$
!----------------------------------------------------------------------
!
! Name:
!   acha_parameters.inc
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! Description:
!   Include file that declares a number of physical constants,
!   apriori estimates and uncertainty values
!
!----------------------------------------------------------------------

 !-----------------------------------------------------------
 ! include other include files to this file
 !-----------------------------------------------------------

 !--- water cloud microphysical model
 include 'acha_water_cloud_microphysical_model.inc'

 !--- ice cloud microphysical model
 include 'acha_ice_cloud_microphysical_model.inc'

 !--- clear-sky covariance terms
 include 'acha_clear_sky_covariances.inc'

 !------------------------------------------------------------
 ! ACHA parameters
 !------------------------------------------------------------
 integer(kind=INT4), parameter, private:: idiag_output = 0 ! (0=off, 1=on)
 real(kind=real4), parameter, private:: count_diag_max = 1000   ! maximum number of diagnostic files allowed

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
 integer, private, parameter :: USE_CIRRUS_FLAG = 0  ! (0 = sym%NO, 1 = sym%YES)
 real (kind=real4), parameter:: EMISSIVITY_MIN_CIRRUS = 0.7
 integer (kind=int4), parameter:: COUNT_MIN_CIRRUS = 3 !5
 integer (kind=int4), parameter:: CIRRUS_BOX_WIDTH_KM = 200 !200

 !---- Lower Cloud Parameters
 integer, private, parameter :: USE_LOWER_INTERP_FLAG = 1  ! (0 = sym%NO, 1 = sym%YES)
 integer (kind=int4), parameter:: COUNT_MIN_LOWER = 3 !5
 integer (kind=int4), parameter:: LOWER_BOX_WIDTH_KM = 200 !200
 real(kind=real4), parameter, private:: TC_LOWER_CLOUD_OFFSET = 10.0  !K

 !---------------------------------------------------------------------
 ! retrieval specific parameters
 !---------------------------------------------------------------------
 integer(kind=INT4), private, parameter :: Num_Levels_Rtm_Prof = 101

 !--- these parameters control the size of matrices
 integer(kind=INT4), parameter :: NUM_PARAM = 4              !number of retrieved parameters

 !--- maximum number of iterations allowed
 integer(kind=INT4), parameter :: Iter_Idx_Max = 10              !maximum number of iterations

 !--- limits on steps taken in retrieval
 real(kind=real4), parameter, private, DIMENSION(NUM_PARAM):: delta_x_max = (/20.0,0.1,0.1,5.0/)

 !--- parameters that control the conversion of Tc uncertainty to Pc and Zc uncertainty
 real(kind=real4), parameter, private:: ABS_LAPSE_RATE_DlnP_DZ_UNCER = 0.0001 !hPa/m
 real(kind=real4), parameter, private:: ABS_LAPSE_RATE_DT_DZ_UNCER = 0.007 !K/m

 !--- parameters that control the bottom-up estimation of Zc and Pc
 real(kind=real4), parameter, private:: MIN_P_INVERSION = 600.0 !min pressure to look for low-level inversion
 real(kind=real4), parameter, private:: MAX_DELTA_T_INVERSION = 20.0 !max temperature difference (surface -cloud) to look for low-level inversion
 real(kind=real4), parameter, private:: DELTA_PSFC_INVERSION = 50.0  !pressure above sfc to look for low-level inversion
 real(kind=real4), parameter, private:: DELTA_T_LAYER_INVERSION = 0.0   !maximum temperature difference to define an inversion layer
 integer(kind=INT1), parameter, private:: ALLOW_STRATOSPHERE_SOLUTION_FLAG = 1  !allow interpolation in stratosphere

 !--- the parameters that provide the apriori values and their uncertainties
!real(kind=real4), parameter, private:: Tc_Ap_Tropo_Offset_Cirrus = 15.0      !apriori Tc for opaque clouds

 real(kind=real4), parameter, private:: Tc_Ap_Uncer_Opaque  = 10.0            !apriori uncertainty of Tc for opaque clouds

!real(kind=real4), parameter, private:: Tc_Ap_Uncer_Cirrus = 20.0  !!           !apriori uncertainty of Tc for cirrus
                                                                    ! Originally it was 20 in the baseline version of ACHA

!real(kind=real4), parameter, private:: Ec_Ap_Uncer_Opaque = 0.2             !apriori uncertainty of ec for opaque clouds

 real(kind=real4), parameter, private:: Ts_Ap_Uncer = 1.0                    !apriori uncertainty of Tsfc
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

 !--- parameter that controls observed negative heights and pressures.
 real(kind=real4), parameter, private:: ZC_FLOOR = 75.0 !Set a bottom limit to Zc.

 !--- cirrus apriori for ice cloud temperature relative to tropopause
 !---- computed in 10 deg lat bands.  First bin is -90 to -80
 integer, parameter, private:: num_lat_cirrus_ap = 18
 real(kind=real4), dimension(num_lat_cirrus_ap), private, parameter:: &
    TC_CIRRUS_MEAN_LAT_VECTOR = (/  &
      17.1184, 16.8644, 14.4241, 15.7231, 14.6659, 17.5565, &
      27.6142, 17.9770, 13.5127, 16.2466, 22.4699, 24.6992, &
      18.7324, 18.3266, 15.2087, 14.8122, 18.0353, 19.9059/)
 real(kind=real4), dimension(num_lat_cirrus_ap), private, parameter:: &
    TC_CIRRUS_STDDEV_LAT_VECTOR = (/  &
      9.58205, 10.7957, 9.96308, 9.75967, 9.11487, 11.3290, &
      18.0040, 16.2742, 14.0061, 15.0945, 13.7971, 14.6018, &
      15.5322, 14.5822, 11.3727, 9.90263, 11.5416, 11.9146/)
