 !--------------------------------------------------------------
 ! parameter used in naive_bayesian_cloud_mask
 !--------------------------------------------------------------

 !--- Set this to 1 to peform diagnostic calculations
 integer, parameter, private:: Do_By_Class_Flag = 1
 real, parameter, private:: Max_Cond_Ratio = 100.0 

 !--- Number of Tests
 integer, parameter, private:: NUMBER_OF_FLAGS = 36
 integer, parameter, private:: NUMBER_OF_FLAG_BYTES = 7
 integer, parameter, private:: NUMBER_OF_NONCLOUD_FLAGS = 18

 !--- Angular and Temperature Limit Thresholds
 real, parameter, private:: Reflectance_Gross_Solzen_Thresh = 80.0     !was 70.0
 real, parameter, private:: Reflectance_Spatial_Solzen_Thresh = 85.0
 real, parameter, private:: Reflectance_Gross_Airmass_Thresh = 100.0 !5.0   !turned off
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
