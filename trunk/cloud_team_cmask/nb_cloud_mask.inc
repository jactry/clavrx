 !--------------------------------------------------------------
 ! parameter used in naive_bayesian_cloud_mask
 !--------------------------------------------------------------

 integer, parameter, private:: int1 = selected_int_kind(1)
 integer, parameter, private:: int2 = selected_int_kind(3)
 integer, parameter, private:: int4 = selected_int_kind(8)
 integer, parameter, private:: int8 = selected_int_kind(10)
 integer, parameter, private:: real4 = selected_real_kind(6,37)
 integer, parameter, private:: real8 = selected_real_kind(15,307)

 !--- missing value
 real, parameter, private :: MISSING_VALUE_REAL4 = -999.0  
 real, parameter, private :: MISSING_VALUE_INT1 = -128          

 !--- numerical constants
 real, parameter, private:: pi = 3.14159265
 real, parameter, private:: dtor = pi/180.0

 !--- Set this to 1 to peform diagnostic calculations
 integer, parameter, private:: Do_By_Class_Flag = 1

 !--- Number of Tests
 integer, parameter, private:: NUMBER_OF_FLAGS = 36
 integer, parameter, private:: NUMBER_OF_FLAG_BYTES = 7
 integer, parameter, private:: NUMBER_OF_NONCLOUD_FLAGS = 18

 !--- Angular and Temperature Limit Thresholds
 real, parameter, private:: Reflectance_Gross_Solzen_Thresh = 80.0     !was 70.0
 real, parameter, private:: Reflectance_Spatial_Solzen_Thresh = 85.0
 real, parameter, private:: Reflectance_Gross_Airmass_Thresh = 5.0     
 real, parameter, private:: Reflectance_Gross_Lunzen_Thresh = 80.0
 real, parameter, private:: Radiance_Lunar_City_Thresh = 2.5e-08
 real, parameter, private:: Emiss_375um_Day_Solzen_Thresh = 85.0       !was 85.0
 real, parameter, private:: Emiss_375um_Night_Solzen_Thresh = 90.0     !was 80.0
 real, parameter, private:: Bt_11um_Cold_Scene_Thresh = 220.0    
 real, parameter, private:: Bt_375um_Cold_Scene_Thresh = 240.0    
 
 !--- EUMETCAST fire detection parameters
 real, parameter, private:: EumetCAST_Fire_Day_Solzen_Thresh = 70.0
 real, parameter, private:: EumetCAST_Fire_Night_Solzen_Thresh = 90.0
 real, private, parameter:: Bt_375um_Eumet_Fire_Day_Thresh = 310.0
 real, private, parameter:: Bt_Diff_Eumet_Fire_Day_Thresh = 8.0
 real, private, parameter:: Stddev_11um_Eumet_Fire_Day_Thresh = 1.0 
 real, private, parameter:: Stddev_375um_Eumet_Fire_Day_Thresh = 4.0
 real, private, parameter:: Bt_375um_Eumet_Fire_Night_Thresh = 290.0
 real, private, parameter:: Bt_Diff_Eumet_Fire_Night_Thresh = 0.0
 real, private, parameter:: Stddev_11um_Eumet_Fire_Night_Thresh = 1.0 
 real, private, parameter:: Stddev_375um_Eumet_Fire_Night_Thresh = 4.0
