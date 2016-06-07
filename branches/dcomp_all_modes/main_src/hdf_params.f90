! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: hdf_params.f90 (src)
!       HDF_PARAMS (program)
!
! PURPOSE: This module contains routines used to read and write to the hdf
!          output files from CLAVR-X
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!
! COPYRIGHT
! (c) This code is copyrighted by the author and all NOAA restrictions apply
!
! Dependencies:  (The following are names of other CLAVR-x modules)
!  CONSTANTS
!  HDF
!  SCALING_PARAMETERS
!
! Calling Sequence:
!  use HDF_PARAMS
!
! Public Routines within this module
!  SCALE_VECTOR_I1_RANK1
!  SCALE_VECTOR_I1_RANK2
!  SCALE_VECTOR_I1_RANK3
!  SCALE_VECTOR_I2_RANK1
!  SCALE_VECTOR_I2_RANK2
!  SCALE_VECTOR_I2_RANK3
!  UNSCALE_VECTOR_I1_RANK1
!  WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES
!  READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES
!  READ_CLAVRX_HDF4_SDS_RANK1
!
!--------------------------------------------------------------------------------------
module HDF_PARAMS

use CONSTANTS, only: &
   int1 &
   , int2 &
   , int4 &
   , real4 &
   , missing_value_int1 &
   , missing_value_int2
   
use HDF, only: &
   SUCCEED &
   , DFNT_INT16 &
   , DFNT_CHAR8 &
   , DFNT_FLOAT32 &
   , DFNT_INT32 &
   , DFNT_INT8



use SCALING_PARAMETERS, only: &
    one_byte_min &
    , one_byte_max &
    , two_byte_min &
    , two_byte_max

implicit none
private
public:: SCALE_VECTOR_I1_RANK1, &
         SCALE_VECTOR_I1_RANK2, &
         SCALE_VECTOR_I1_RANK3, &
         SCALE_VECTOR_I2_RANK1, &
         SCALE_VECTOR_I2_RANK2, &
         SCALE_VECTOR_I2_RANK3, &
         UNSCALE_VECTOR_I1_RANK1, &
         WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES,  &
         READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES, &
         READ_CLAVRX_HDF4_SDS_RANK1


character(len=256), save, public:: renav_data_from


!--- scaling options
integer(kind=int1), parameter, public:: NO_SCALING = 0
integer(kind=int1), parameter, public:: LINEAR_SCALING = 1
integer(kind=int1), parameter, public:: LOG10_SCALING = 2
integer(kind=int1), parameter, public:: SQUARE_ROOT_SCALING = 3

contains

!-------------------------------------------------------------------------
! routine to write global attributes to clavrx files
!
!-------------------------------------------------------------------------
 subroutine WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES(hdf_file_id,data_type,file_name,file_1b, &
                           resolution_km, &
                           start_year,end_year,start_day,end_day,start_time,end_time,&
                           num_cells,num_cells_with_data,grid_format,dlat, &
                           therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                           modis_clr_alb_flag, nwp_opt, ch1_gain_low, ch1_gain_high, &
                           ch1_switch_count, ch1_dark_count, &
                           ch2_gain_low, ch2_gain_high, &
                           ch2_switch_count, ch2_dark_count, &
                           ch3a_gain_low, ch3a_gain_high, &
                           ch3a_switch_count, ch3a_dark_count, &
                           sun_earth_distance, &
                           c1, c2, a1_20, a2_20, nu_20, &
                           a1_31, a2_31, nu_31, a1_32, a2_32, nu_32, &
                           solar_Ch20_nu,timerr_seconds, &
                           acha_mode, dcomp_mode, wmo_sc_code, &
                           platform_name,sensor_name,dark_name,mask_name)
                           

 integer(kind=int4), intent(in):: hdf_file_id,num_cells,num_cells_with_data
 integer, intent(in):: therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                            modis_clr_alb_flag, nwp_opt

 integer(kind=int4), intent(in):: start_time
 integer(kind=int4), intent(in):: end_time
 integer(kind=int4), intent(in):: acha_mode
 integer(kind=int4), intent(in):: dcomp_mode
 integer(kind=int2), intent(in)::  start_year,end_year,start_day,end_day
 real(kind=real4), intent(in)::  resolution_km, dlat, ch1_gain_low, ch1_gain_high, &
                           ch2_gain_low, ch2_gain_high, &
                           ch3a_gain_low, ch3a_gain_high, &
                           ch1_switch_count, ch1_dark_count, &
                           ch2_switch_count, ch2_dark_count, &
                           ch3a_switch_count, ch3a_dark_count, &
                           sun_earth_distance, &
                           c1, c2, a1_20, a2_20, nu_20, &
                           a1_31, a2_31, nu_31, a1_32, a2_32, nu_32, &
                           solar_Ch20_nu,timerr_seconds
 integer(kind=int4):: wmo_sc_code

 character(len=*), intent(in):: file_name
 character(len=*), intent(in):: file_1b
 character(len=*), intent(in):: data_type
 character(len=*), intent(in):: grid_format
 character(len=*), intent(in):: platform_name
 character(len=*), intent(in):: sensor_name
 character(len=*), intent(in):: dark_name
 character(len=*), intent(in):: mask_name

 integer:: sfscatt, sfsnatt, hglibver
 integer:: major_v, minor_v, release
 integer(kind=int2):: istatus = 0
 character(len=80) :: hdf_ver
 character(len=36) :: machine

 character(len=6), parameter:: Conventions_String = "CF-1.6"
 character(len=38), parameter:: Metadata_Conventions_String = "CF-1.6, Unidata Dataset Discovery v1.0"
 character(len=42), parameter:: Standard_Name_Vocabulary_String = "CF Standard Name Table (v25, 05 July 2013)"
 character(len=13), parameter:: Naming_Authority_String = "gov.noaa.ncdc"
 character(len=37), parameter:: License_String = "No constraints on data access or use."


! Definition of strings used as HDF attributes.
!
! version history 4.1 - delivered to OSDPD in November 2006
! version history 4.2 - demonstrated on METOP
! version history 4.3 - included surface emissivity fields
! version history 4.4 - included lrc
! version history 5.0 - included modis white sky and ash protoype 
! version history 5.1 - rtm structures now 101 levels and 
!                       reorganized level-1b ingest to 
!                       read segment all at once prior to processing
!                       first version with working volcanic ash
! version history 5.2    bayesian cloud mask and DCOMP
! version history 6.0  - MODIS capability begin
! version history 6.5  - VIIRS capability begin

 character(len=36), parameter :: creator0 = "CLAVR-x + PATMOS-x ", &
                                 plang = "F90"

 character(len=36)  :: creator
 character(len=100):: Title_String
 character(len=100):: Calibration_String
 character(len=100):: Product_Version_String
 character(len=100):: Status_String
 character(len=100):: Institution_String
 character(len=100):: Program_String
 character(len=500):: Summary_String 
 character(len=200):: Variable_String
 character(len=500):: Keywords_String
 character(len=200):: Keywords_Vocabulary_String
 character(len=100):: Time_Coverage_Resolution_String
 character(len=100):: Metadata_Link_String
 character(len=100):: Spatial_Resolution_String

 include 'version.inc'

 !complete the creator string with the version number
 creator = trim(creator0)//trim(Product_Version_String)


 call getenv ("HOST",machine)
 if (len_trim(machine) == 0) call getenv ("HOSTNAME",machine)
 if (len_trim(machine) == 0) machine = "unknown"

!--- determine HDF library version
if (hglibver(major_v, minor_v, release, hdf_ver) /= SUCCEED) then
   print *, "could not determine HDF library version"
   stop 961
end if

!---- describe CLAVR-x
!istatus = sfscatt(hdf_file_id, "PROCESSING_NOTE", DFNT_CHAR8, len_trim(PROCESSINGSTRING), trim(PROCESSINGSTRING))+istatus
istatus = sfscatt(hdf_file_id, "HDF_LIB_VERSION", DFNT_CHAR8, len_trim(hdf_ver), trim(hdf_ver))+istatus
istatus = sfscatt(hdf_file_id, "MACHINE", DFNT_CHAR8, len_trim(machine), trim(machine))+istatus
istatus = sfscatt(hdf_file_id, "PROGLANG", DFNT_CHAR8, len_trim(plang), trim(plang))+istatus

!--- CF compliant global attributes required for NCDC delivery
istatus = sfscatt(hdf_file_id, "date_created", DFNT_CHAR8, len_trim(hdf_timestamp()), trim(hdf_timestamp()))+istatus
istatus = sfscatt(hdf_file_id, "product_version", DFNT_CHAR8,  &
                                len_trim(Product_Version_String), trim(Product_Version_String))+istatus
istatus = sfscatt(hdf_file_id,"summary", DFNT_CHAR8,len_trim(Summary_String),trim(Summary_String)) + istatus
istatus = sfscatt(hdf_file_id,"cdr_variable", DFNT_CHAR8,len_trim(Variable_String),trim(Variable_String)) + istatus
istatus = sfscatt(hdf_file_id,"institution", DFNT_CHAR8,len_trim(Institution_String),trim(Institution_String)) + istatus
istatus = sfscatt(hdf_file_id,"cdr_program", DFNT_CHAR8, len_trim(Program_String),trim(Program_String)) + istatus
istatus = sfscatt(hdf_file_id,"title", DFNT_CHAR8,len_trim(Title_String),trim(Title_String)) + istatus
istatus = sfscatt(hdf_file_id,"calibration_version", DFNT_CHAR8,len_trim(Calibration_String),trim(Calibration_String)) + istatus
istatus = sfscatt(hdf_file_id,"keywords", DFNT_CHAR8,len_trim(Keywords_String),trim(Keywords_String)) + istatus
istatus = sfscatt(hdf_file_id,"keywords_vocabulary", DFNT_CHAR8,len_trim(Keywords_Vocabulary_String),trim(Keywords_Vocabulary_String)) + istatus
istatus = sfscatt(hdf_file_id,"time_coverage_resolution", DFNT_CHAR8,len_trim(Time_Coverage_Resolution_String),trim(Time_Coverage_Resolution_String)) + istatus
istatus = sfscatt(hdf_file_id,"metadata_link", DFNT_CHAR8,len_trim(Metadata_Link_String),trim(Metadata_Link_String)) + istatus
istatus = sfscatt(hdf_file_id,"spatial_resolution", DFNT_CHAR8,len_trim(Spatial_Resolution_String),trim(Spatial_Resolution_String)) + istatus
istatus = sfscatt(hdf_file_id,"Conventions", DFNT_CHAR8,len_trim(Conventions_String),trim(Conventions_String)) + istatus
istatus = sfscatt(hdf_file_id,"title", DFNT_CHAR8,len_trim(Title_String),Title_String) + istatus
istatus = sfscatt(hdf_file_id,"Metadata_Conventions", DFNT_CHAR8,  &
                   len_trim(Metadata_Conventions_String),trim(Metadata_Conventions_String)) + istatus
istatus = sfscatt(hdf_file_id,"standard_name_vocabulary", DFNT_CHAR8,  &
                   len_trim(Standard_Name_Vocabulary_String),trim(Standard_Name_Vocabulary_String)) + istatus
istatus = sfscatt(hdf_file_id,"naming_authority", DFNT_CHAR8, &
                   len_trim(Naming_Authority_String),trim(Naming_Authority_String)) + istatus
istatus = sfscatt(hdf_file_id,"license", DFNT_CHAR8, len_trim(License_String),trim(License_String)) + istatus
istatus = sfscatt(hdf_file_id, "sensor", DFNT_CHAR8,len_trim(sensor_name),trim(sensor_name))+istatus
istatus = sfscatt(hdf_file_id, "platform", DFNT_CHAR8,len_trim(platform_name),trim(platform_name))+istatus

!--- describe the data
istatus = sfscatt(hdf_file_id, "FILENAME", DFNT_CHAR8, len_trim(file_name), trim(file_name))+istatus
istatus = sfscatt(hdf_file_id, "L1B", DFNT_CHAR8, len_trim(file_1b), trim(file_1b))+istatus
!istatus = sfsnatt(hdf_file_id, "RESOLUTION_KM", DFNT_FLOAT32, 1, resolution_km)+istatus
!istatus = sfsnatt(hdf_file_id, "spatial_resolution", DFNT_FLOAT32, 1, resolution_km)+istatus
istatus = sfsnatt(hdf_file_id, "START_YEAR", DFNT_INT16, 1, start_year)+istatus
istatus = sfsnatt(hdf_file_id, "START_DAY", DFNT_INT16, 1, start_day)+istatus
istatus = sfsnatt(hdf_file_id, "START_TIME", DFNT_FLOAT32, 1, start_time/3600000.0)+istatus
istatus = sfsnatt(hdf_file_id, "END_YEAR", DFNT_INT16, 1, end_year)+istatus
istatus = sfsnatt(hdf_file_id, "END_DAY", DFNT_INT16, 1, end_day)+istatus
istatus = sfsnatt(hdf_file_id, "END_TIME", DFNT_FLOAT32, 1, end_time/3600000.0)+istatus
istatus = sfsnatt(hdf_file_id, "ACHA_MODE", DFNT_INT32, 1, acha_mode)+istatus
istatus = sfsnatt(hdf_file_id, "DCOMP_MODE", DFNT_INT32, 1, dcomp_mode)+istatus
istatus = sfsnatt(hdf_file_id, "WMO_SATELLITE_CODE", DFNT_INT32, 1, wmo_sc_code)+istatus
istatus = sfscatt(hdf_file_id, "REFL_0_65UM_NOM_DARK_COMPOSITE_NAME", DFNT_CHAR8, &
                               len_trim(dark_name),trim(dark_name))+istatus
istatus = sfscatt(hdf_file_id, "NAIVE_BAYESIAN_CLOUD_MASK_NAME", DFNT_CHAR8, &
                               len_trim(mask_name),trim(mask_name))+istatus

!--- data type
istatus = sfscatt(hdf_file_id, "DATA_TYPE", DFNT_CHAR8,len_trim(data_type),trim(data_type))+istatus


!--- NCDC Attributes

!-- other global attributes for the level3 file.
     if (data_type(1:4) == "GRID") then
      istatus = sfsnatt(hdf_file_id, "NUM_CELLS_TOTAL", DFNT_INT32,1,num_cells)+istatus
      istatus = sfsnatt(hdf_file_id, "NUM_CELLS_WITH_DATA", DFNT_INT32,1,num_cells_with_data)+istatus
      istatus = sfsnatt(hdf_file_id, "GRIDCELL_RESOLUTION", DFNT_FLOAT32,1,dlat)+istatus
      istatus = sfscatt(hdf_file_id, "GRIDCELL_RESOLUTION_UNIT", DFNT_CHAR8,6,"degree")+istatus
      if (grid_format(1:10) == "EQUAL_AREA") then
       istatus = sfscatt(hdf_file_id, "GRIDCELL_FORMAT", DFNT_CHAR8,10,"EQUAL_AREA")+istatus
      else
       istatus = sfscatt(hdf_file_id, "GRIDCELL_FORMAT", DFNT_CHAR8,11,"EQUAL_ANGLE")+istatus
      endif
    endif

!---- processing flags
 istatus = sfscatt(hdf_file_id, "USE_1B_THERMAL_CALIBRATION_FLAG", DFNT_INT32,1,therm_cal_1b)+istatus
 istatus = sfscatt(hdf_file_id, "USE_1B_REFLECTANCE_CALIBRATION_FLAG", DFNT_INT32,1,Ref_cal_1b)+istatus
! istatus = sfscatt(hdf_file_id, "RENAVIGATION_DATA_FROM", DFNT_CHAR8, &
!     len_trim(renav_data_from), trim(renav_data_from))+istatus
 istatus = sfscatt(hdf_file_id, "RENAVIGATION_FLAG", DFNT_INT32,1,nav_opt)+istatus
 istatus = sfscatt(hdf_file_id, "USE_SST_ANALYSIS_FLAG", DFNT_INT32,1,use_sst_anal)+istatus
 istatus = sfsnatt(hdf_file_id, "NWP_OPT", DFNT_INT32,1,nwp_opt)+istatus
 istatus = sfsnatt(hdf_file_id, "MODIS_CLEAR_SKY_REFLECTANCE_FLAG", DFNT_INT32,1,modis_clr_alb_flag)+istatus

!-- reflectance channel calibration
istatus = sfsnatt(hdf_file_id, "CH1_GAIN_LOW", DFNT_FLOAT32,1,ch1_gain_low)+istatus
istatus = sfsnatt(hdf_file_id, "CH1_GAIN_HIGH", DFNT_FLOAT32,1,ch1_gain_high)+istatus
istatus = sfsnatt(hdf_file_id, "CH1_SWITCH_COUNT", DFNT_FLOAT32,1,ch1_switch_count)+istatus
istatus = sfsnatt(hdf_file_id, "CH1_DARK_COUNT", DFNT_FLOAT32,1,ch1_dark_count)+istatus
istatus = sfsnatt(hdf_file_id, "CH2_GAIN_LOW", DFNT_FLOAT32,1,ch2_gain_low)+istatus
istatus = sfsnatt(hdf_file_id, "CH2_GAIN_HIGH", DFNT_FLOAT32,1,ch2_gain_high)+istatus
istatus = sfsnatt(hdf_file_id, "CH2_SWITCH_COUNT", DFNT_FLOAT32,1,ch2_switch_count)+istatus
istatus = sfsnatt(hdf_file_id, "CH2_DARK_COUNT", DFNT_FLOAT32,1,ch2_dark_count)+istatus
istatus = sfsnatt(hdf_file_id, "CH3A_GAIN_LOW", DFNT_FLOAT32,1,ch3a_gain_low)+istatus
istatus = sfsnatt(hdf_file_id, "CH3A_GAIN_HIGH", DFNT_FLOAT32,1,ch3a_gain_high)+istatus
istatus = sfsnatt(hdf_file_id, "CH3A_SWITCH_COUNT", DFNT_FLOAT32,1,ch3a_switch_count)+istatus
istatus = sfsnatt(hdf_file_id, "CH3A_DARK_COUNT", DFNT_FLOAT32,1,ch3a_dark_count)+istatus
istatus = sfsnatt(hdf_file_id, "SUN_EARTH_DISTANCE", DFNT_FLOAT32,1,sun_earth_distance)+istatus

!--- thermal calibration constants
istatus = sfsnatt(hdf_file_id, "C1", DFNT_FLOAT32,1,c1)+istatus
istatus = sfsnatt(hdf_file_id, "C2", DFNT_FLOAT32,1,c2)+istatus
istatus = sfsnatt(hdf_file_id, "A_20", DFNT_FLOAT32,1,a1_20)+istatus
istatus = sfsnatt(hdf_file_id, "B_20", DFNT_FLOAT32,1,a2_20)+istatus
istatus = sfsnatt(hdf_file_id, "NU_20", DFNT_FLOAT32,1,nu_20)+istatus
istatus = sfsnatt(hdf_file_id, "A_31", DFNT_FLOAT32,1,a1_31)+istatus
istatus = sfsnatt(hdf_file_id, "B_31", DFNT_FLOAT32,1,a2_31)+istatus
istatus = sfsnatt(hdf_file_id, "NU_31", DFNT_FLOAT32,1,nu_31)+istatus
istatus = sfsnatt(hdf_file_id, "A_32", DFNT_FLOAT32,1,a1_32)+istatus
istatus = sfsnatt(hdf_file_id, "B_32", DFNT_FLOAT32,1,a2_32)+istatus
istatus = sfsnatt(hdf_file_id, "NU_32", DFNT_FLOAT32,1,nu_32)+istatus
istatus = sfsnatt(hdf_file_id, "SOLAR_20_NU", DFNT_FLOAT32,1,solar_Ch20_nu)+istatus
istatus = sfsnatt(hdf_file_id, "TIME_ERROR_SECONDS", DFNT_FLOAT32,1,timerr_seconds)+istatus

if (istatus /= 0) then
  print *, "error writing run-control flags as global HDF attributes"
  stop 963
endif

contains
function hdf_timestamp() result(string)
   character (len = 36) ::string
   character(len=10), dimension(3):: ctime

   call date_and_time(ctime(1), ctime(2), ctime(3))
   ! Timestamp string format in accordance with the ISO-8601 standard.
   string = ctime(1)(1:4)//"-"//ctime(1)(5:6)//"-"//ctime(1)(7:8)&
                //"T"//ctime(2)(1:2)//":"//ctime(2)(3:4)//":"//ctime(2)(5:6)&
                //ctime(3)(1:3)//":"//ctime(3)(4:5)
   return
end function hdf_timestamp
 end subroutine WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES

!------------------------------------------------------------------------------------------
! READ GLOBAL ATTRIBUTES FROM A CLAVR-x HDF FILE
!------------------------------------------------------------------------------------------
 subroutine READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES(hdf_file_id,data_type,file_name,file_1b, &
                           resolution_km, &
                           start_year,end_year,start_day,end_day,start_time,end_time,&
                           num_cells,num_cells_with_data,grid_format,grid_resolution, &
                           therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                           modis_clr_alb_flag, nwp_opt, ch1_gain_low, ch1_gain_high, &
                           ch1_switch_count, ch1_dark_count, &
                           ch2_gain_low, ch2_gain_high, &
                           ch2_switch_count, ch2_dark_count, &
                           ch3a_gain_low, ch3a_gain_high, &
                           ch3a_switch_count, ch3a_dark_count, &
                           sun_earth_distance, &
                           c1, c2, a_20, b_20, nu_20, &
                           a_31, b_31, nu_31, a_32, b_32, nu_32, &
                           solar_Ch20_nu, timerr_seconds, &
                           acha_mode, dcomp_mode, wmo_sc_code, platform_name, sensor_name, &
                           dark_name, mask_name, &
                           creator, plang, hdf_ver,hdf_timestamp )
                                                                                                                                                         
                                                                                                                                                         
     integer(kind=int4), intent(in):: hdf_file_id
     integer(kind=int4), intent(out):: num_cells,num_cells_with_data
     integer, intent(out):: therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                             modis_clr_alb_flag, nwp_opt
     integer(kind=int4), intent(out)::  start_time,end_time
     integer(kind=int4), intent(out)::  acha_mode, dcomp_mode, wmo_sc_code
     integer(kind=int2), intent(out)::  start_year,end_year,start_day,end_day
     real(kind=real4), intent(out)::  resolution_km 
     real(kind=real4), intent(out)::  grid_resolution, ch1_gain_low, ch1_gain_high, &
                           ch2_gain_low, ch2_gain_high, &
                           ch3a_gain_low, ch3a_gain_high, &
                           ch1_switch_count, ch1_dark_count, &
                           ch2_switch_count, ch2_dark_count, &
                           ch3a_switch_count, ch3a_dark_count, &
                           sun_earth_distance, &
                           c1, c2, a_20, b_20, nu_20, &
                           a_31, b_31, nu_31, a_32, b_32, nu_32, &
                           solar_Ch20_nu,timerr_seconds
                                                                                                                                                         
 character(len=*), intent(out):: file_name, file_1b, creator,  grid_format,  &
                                 plang, hdf_ver, hdf_timestamp, data_type, &
                                 platform_name, sensor_name, dark_name, mask_name
                                                                                                                                                         
integer(kind=int2):: istatus = 0
integer(kind=int4):: blank_int4
real(kind=real4):: blank_real4

!--- hdf  calls
integer:: sfrcatt, sfrnatt, sffattr


blank_int4 = 0
blank_real4 = 0.0

!--- attributes about the code used to make this data
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"PROCESSOR"), creator) + istatus
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"CREATED"), hdf_timestamp) + istatus
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"HDF_LIB_VERSION"), hdf_ver(1:48)) + istatus
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"PROGLANG"), plang) + istatus

!--- file names
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"FILENAME"), file_name) + istatus
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"L1B"), file_1b) + istatus

!--- temporal attributes
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "START_YEAR"), start_year) + istatus
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "START_DAY"), start_day) + istatus
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "START_TIME"), start_time) + istatus
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "END_YEAR"), end_year) + istatus
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "END_DAY"), end_day) + istatus
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "END_TIME"), end_time) + istatus
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "WMO_SATELLITE_CODE"),wmo_sc_code) + istatus
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "platform"),platform_name) + istatus           
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "sensor"),sensor_name) + istatus          

 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "REFL_0_65UM_NOM_DARK_COMPOSITE_NAME"), &
                                                     dark_name) + istatus
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "NAIVE_BAYESIAN_CLOUD_MASK_NAME"), &
                                                     mask_name) + istatus

 !--- algorithm options
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "ACHA_MODE"), acha_mode) + istatus
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "DCOMP_MODE"), dcomp_mode) + istatus

!--- data type
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id,"DATA_TYPE"), data_type) + istatus

!--- grid variables
if (data_type(1:4) == "GRID") then
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NUM_CELLS_WITH_DATA"), num_cells_with_data) + istatus
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NUM_CELLS_TOTAL"), num_cells) + istatus
 istatus = sfrcatt(hdf_file_id, sffattr(hdf_file_id, "GRIDCELL_FORMAT"), grid_format) + istatus
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "GRIDCELL_RESOLUTION"), grid_resolution) + istatus
 resolution_km = -999.0
else
 num_cells_with_data = blank_int4
 num_cells = blank_int4
 grid_format = "     " 
 grid_resolution = blank_real4
 resolution_km = -999.0 !if RESOLUTION_KM global attribute is not present, missing value
 istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "RESOLUTION_KM"), resolution_km) + istatus
endif

!--- flags
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "USE_1B_THERMAL_CALIBRATION_FLAG"), therm_cal_1b)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "USE_1B_REFLECTANCE_CALIBRATION_FLAG"), Ref_cal_1b)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "USE_SST_ANALYSIS_FLAG"), use_sst_anal)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NWP_FLAG"), nwp_opt)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "MODIS_CLEAR_SKY_REFLECTANCE_FLAG"), modis_clr_alb_flag)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "RENAVIGATION_FLAG"), nav_opt)

!--- calibration attributes
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "C1"), c1)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "C2"), c2)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "A_20"), a_20)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "B_20"), b_20)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NU_20"), nu_20)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "A_31"), a_31)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "B_31"), b_31)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NU_31"), nu_31)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "A_32"), a_32)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "B_32"), b_32)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "NU_32"), nu_32)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "SOLAR_20_NU"), solar_Ch20_nu)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "TIME_ERROR_SECONDS"), timerr_seconds)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "SUN_EARTH_DISTANCE"), sun_earth_distance)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH1_GAIN_LOW"), ch1_gain_low)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH1_GAIN_HIGH"), ch1_gain_high)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH1_SWITCH_COUNT"), ch1_switch_count)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH1_DARK_COUNT"), ch1_dark_count)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH2_GAIN_LOW"), ch2_gain_low)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH2_GAIN_HIGH"), ch2_gain_high)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH2_SWITCH_COUNT"), ch2_switch_count)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH2_DARK_COUNT"), ch2_dark_count)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH6_GAIN_LOW"), ch3a_gain_low)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH6_GAIN_HIGH"), ch3a_gain_high)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH6_SWITCH_COUNT"), ch3a_switch_count)
istatus = sfrnatt(hdf_file_id, sffattr(hdf_file_id, "CH6_DARK_COUNT"), ch3a_dark_count)

end subroutine READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES

!-----------------------------------------------------------------------------------------------------
! VECTOR SCALING ROUTINES
!
! iscaled = 0 = no scaling
!           1 = linear
!           2 = log10
!           3 = sqrt
!-----------------------------------------------------------------------------------------------------
 subroutine SCALE_VECTOR_I1_RANK1(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i1)
   real, dimension(:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int1), dimension(:),  intent(out):: temp_i1
   real, dimension(size(temp_r4,1)):: scratch_r4

   scratch_r4 = 0.0

!---- linear
    if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = one_byte_min + scratch_r4 * (one_byte_max - one_byte_min)
    endif
                                                                                                                                                          
!---- log10
    if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i1 = one_byte_min + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * (one_byte_max - one_byte_min)
     endif
!---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = one_byte_min + sqrt(scratch_r4) * (one_byte_max - one_byte_min)
    endif
                                                                                                                                                          
!--- set scaled missing values
    where (temp_r4 == unscaled_missing) 
         temp_i1 = missing_value_int1
    endwhere

 end subroutine SCALE_VECTOR_I1_RANK1

 subroutine SCALE_VECTOR_I1_RANK2(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i1)
   real, dimension(:,:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int1), dimension(:,:),  intent(out):: temp_i1
   real, dimension(size(temp_r4,1),size(temp_r4,2)):: scratch_r4

   scratch_r4 = 0.0
!---- linear
    if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = one_byte_min + scratch_r4 * (one_byte_max - one_byte_min)
    endif
!---- log10
    if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i1 = one_byte_min + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * (one_byte_max - one_byte_min)
     endif
!---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = one_byte_min + sqrt(scratch_r4) * (one_byte_max - one_byte_min)
    endif
!--- set scaled missing values
    where (temp_r4 == unscaled_missing)
         temp_i1 = missing_value_int1
    endwhere
 end subroutine SCALE_VECTOR_I1_RANK2

 subroutine SCALE_VECTOR_I1_RANK3(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i1)
   real, dimension(:,:,:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int1), dimension(:,:,:),  intent(out):: temp_i1
   real, dimension(size(temp_r4,1),size(temp_r4,2),size(temp_r4,3)):: scratch_r4

   scratch_r4 = 0.0
!---- linear
    if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = one_byte_min + scratch_r4 * (one_byte_max - one_byte_min)
    endif
!---- log10
    if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i1 = one_byte_min + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * (one_byte_max - one_byte_min)
     endif
!---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i1 = one_byte_min + sqrt(scratch_r4) * (one_byte_max - one_byte_min)
    endif
!--- set scaled missing values
    where (temp_r4 == unscaled_missing)
         temp_i1 = missing_value_int1
    endwhere
 end subroutine SCALE_VECTOR_I1_RANK3


 subroutine SCALE_VECTOR_I2_RANK1(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i2)
   real, dimension(:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int2), dimension(:), intent(out):: temp_i2
   real, dimension(size(temp_r4,1)):: scratch_r4

!---- linear
    if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = two_byte_min + scratch_r4 * (two_byte_max - two_byte_min)
    endif
!---- log10
    if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i2 = two_byte_min + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * (two_byte_max - two_byte_min)
    endif
!---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = two_byte_min + sqrt(scratch_r4) * (two_byte_max - two_byte_min)
    endif
!--- set scaled missing values
    where (temp_r4 == unscaled_missing)
         temp_i2 = missing_value_int2
    endwhere
 end subroutine SCALE_VECTOR_I2_RANK1

 subroutine SCALE_VECTOR_I2_RANK2(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i2)
   real, dimension(:,:), intent(in):: temp_r4
   integer, intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int2), dimension(:,:), intent(out):: temp_i2
   real, dimension(size(temp_r4,1),size(temp_r4,2)):: scratch_r4
        
!---- linear
    if (iscaled == 1) then
      
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       
       temp_i2 = two_byte_min + scratch_r4 * (two_byte_max - two_byte_min)
      
    endif
!---- log10
    if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i2 = two_byte_min + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * (two_byte_max - two_byte_min)
    endif
!---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = two_byte_min + sqrt(scratch_r4) * (two_byte_max - two_byte_min)
    endif
!--- set scaled missing values
    where (temp_r4 == unscaled_missing)
         temp_i2 = missing_value_int2
    endwhere
 end subroutine SCALE_VECTOR_I2_RANK2

 subroutine SCALE_VECTOR_I2_RANK3(temp_r4,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_i2)
   real, dimension(:,:,:), intent(in):: temp_r4
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   integer(kind=int2), dimension(:,:,:), intent(out):: temp_i2
   real, dimension(size(temp_r4,1),size(temp_r4,2),size(temp_r4,3)):: scratch_r4

!---- linear
    if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = two_byte_min + scratch_r4 * (two_byte_max - two_byte_min)
    endif
!---- log10
    if (iscaled == 2) then
      scratch_r4 = unscaled_missing
      where(temp_r4 > 0.0)
        scratch_r4 = log10(temp_r4)
      end where
      temp_i2 = two_byte_min + ((max(unscaled_min,min(unscaled_max,scratch_r4)) - unscaled_min)/ &
                    (unscaled_max - unscaled_min)) * (two_byte_max - two_byte_min)
    endif
!---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,(temp_r4 - unscaled_min)/(unscaled_max - unscaled_min)))
       temp_i2 = two_byte_min + sqrt(scratch_r4) * (two_byte_max - two_byte_min)
    endif
!--- set scaled missing values
    where (temp_r4 == unscaled_missing)
         temp_i2 = missing_value_int2
    endwhere
 end subroutine SCALE_VECTOR_I2_RANK3

!-----------------------------------------------------------------------------------------------------
! VECTOR UNSCALING ROUTINES
!-----------------------------------------------------------------------------------------------------
 subroutine UNSCALE_VECTOR_I1_RANK1(temp_i1,iscaled,unscaled_min,unscaled_max,unscaled_missing,temp_r4)
   integer(kind=int1), dimension(:), intent(in):: temp_i1
   integer(kind=int1), intent(in):: iscaled
   real, intent(in):: unscaled_min, unscaled_max, unscaled_missing
   real, dimension(:),  intent(out):: temp_r4
   real, dimension(size(temp_r4,1)):: scratch_r4
   integer (kind=int1):: scaled_min,scaled_max,scaled_missing

   scaled_min = one_byte_min
   scaled_max = one_byte_max
   scaled_missing = missing_value_int1

   scratch_r4 = 0.0
!---- linear
    if (iscaled == 1) then
       scratch_r4 = min(1.0,max(0.0,real(temp_i1 - scaled_min)/real(scaled_max - scaled_min)))
       temp_r4 = unscaled_min + scratch_r4 * (unscaled_max - unscaled_min)
    endif
                                                                                                                                                         
!---- log10
    if (iscaled == 2) then
       scratch_r4 = min(1.0,max(0.0,real(temp_i1 - scaled_min)/real(scaled_max - scaled_min)))
       temp_r4 = unscaled_min + scratch_r4 * (unscaled_max - unscaled_min)
       temp_r4 = 10**(temp_r4)
     endif

!---- square root
    if (iscaled == 3) then
       scratch_r4 = min(1.0,max(0.0,real(temp_i1 - scaled_min)/real(scaled_max - scaled_min)))
       temp_r4 = unscaled_min + (scratch_r4**2) * (unscaled_max - unscaled_min)
    endif
                                                                                                                                                         
                                                                                                                                                         
!--- set scaled missing values
    where (temp_i1 == scaled_missing)
         temp_r4 = unscaled_missing
    endwhere
                                                                                                                                                         
 end subroutine UNSCALE_VECTOR_I1_RANK1

!----------------------------------------------------------------------------------------------------
! HDF4 READ ROUTINES
! this routine reads in level3 sds's, it assumes
! 1. - you know the size of the array (num_cells_with_data)
!-----------------------------------------------------------------------------------------------------
 subroutine READ_CLAVRX_HDF4_SDS_RANK1(sd_id,sds_dim_input,sds_name,sds_data, &
                                       sds_data_type,scaled,sds_units, &
                                       unscaled_min,unscaled_max,unscaled_missing,&
                                       scaled_min,scaled_max,scaled_missing,&
                                       istatus)
  integer, intent(in):: sd_id,sds_dim_input
  character(len=*), intent(in):: sds_name
  real, intent(out):: unscaled_min, unscaled_max, unscaled_missing
  integer, intent(out):: scaled_min, scaled_max, scaled_missing
  real, intent(out), dimension(:):: sds_data
  integer(kind=int1), intent(out):: scaled
  character(len=*), intent(out):: sds_units
  integer, intent(out):: sds_data_type,istatus
  integer:: sds_id,sds_dim1
  integer:: num_attrs, sds_rank
  integer, dimension(1):: dimsizes
  character(72):: sds_name_temp

  real, dimension(size(sds_data)):: sds_data_temp

  integer, dimension(1):: sds_dims_1d, sds_start_1d, sds_stride_1d, sds_edges_1d


  integer(kind=int1), dimension(:), allocatable:: temp_i1
  integer(kind=int2), dimension(:), allocatable:: temp_i2
  integer(kind=int4), dimension(:), allocatable:: temp_i4
  real(kind=real4), dimension(:), allocatable:: temp_r4

  ! HDF function declarations
  integer:: sfselect,sfn2index,sfrnatt, sfrcatt, sffattr, sfrdata,sfendacc,sfginfo
 
  istatus = 0

!----------------------------------------------------------------------------
! open sds for reading
!----------------------------------------------------------------------------
sds_id = sfselect(sd_id, sfn2index(sd_id,sds_name))

istatus = sfginfo(sds_id, sds_name_temp, sds_rank, dimsizes, sds_data_type, num_attrs) + istatus
sds_dim1 = dimsizes(1)
sds_dims_1d = (/sds_dim1/)
sds_start_1d = (/ 0 /)      
sds_stride_1d = (/ 1 /)
sds_edges_1d = (/ sds_dim1 /)   
sds_units = " "

!--- read sds attributes
 istatus = sfrnatt(sds_id, sffattr(sds_id,"SCALED"), scaled)
 istatus = sfrcatt(sds_id, sffattr(sds_id,"UNITS"), sds_units)
 istatus = sfrnatt(sds_id, sffattr(sds_id,"RANGE_MISSING"), unscaled_missing)

!-- if scaled, read attributes that allow unscaling
 if (scaled > 0) then
   istatus = sfrnatt(sds_id, sffattr(sds_id,"SCALED_MISSING"), scaled_missing)
   istatus = sfrnatt(sds_id, sffattr(sds_id,"SCALED_MIN"), scaled_min)
   istatus = sfrnatt(sds_id, sffattr(sds_id,"SCALED_MAX"), scaled_max)
   istatus = sfrnatt(sds_id, sffattr(sds_id,"RANGE_MIN"), unscaled_min)
   istatus = sfrnatt(sds_id, sffattr(sds_id,"RANGE_MAX"), unscaled_max)
 endif

!-- check dimension against expectations
    if (sds_dim_input /= sds_dim1) then
      print *, "error, sds dimension differs from expectations, stopping", sds_dim_input, sds_dim1
      stop
     endif

!--- allocate arrays for holding data, read data and store in output array
if (sds_data_type == DFNT_INT8) then
    allocate(temp_i1(sds_dim1))
    istatus = sfrdata(sds_id, sds_start_1d, sds_stride_1d, sds_edges_1d, temp_i1) + istatus   
    sds_data = real(temp_i1)
elseif (sds_data_type == DFNT_INT16) then
    allocate(temp_i2(sds_dim1))
    istatus = sfrdata(sds_id, sds_start_1d, sds_stride_1d, sds_edges_1d, temp_i2) + istatus   
    sds_data = real(temp_i2)
elseif (sds_data_type == DFNT_INT32) then
    allocate(temp_i4(sds_dim1))
    istatus = sfrdata(sds_id, sds_start_1d, sds_stride_1d, sds_edges_1d, temp_i4) + istatus   
    sds_data = real(temp_i4)
elseif (sds_data_type == DFNT_FLOAT32) then
    allocate(temp_r4(sds_dim1))
    istatus = sfrdata(sds_id, sds_start_1d, sds_stride_1d, sds_edges_1d, temp_r4) + istatus   
    sds_data = temp_r4
else
    print *, "attempt to read unsupported data type, stopping"
    stop
endif

!---deallocate temp arrays
   if (allocated(temp_i1)) deallocate(temp_i1)
   if (allocated(temp_i2)) deallocate(temp_i2)
   if (allocated(temp_i4)) deallocate(temp_i4)
   if (allocated(temp_r4)) deallocate(temp_r4)

!--- close sds
   istatus = sfendacc(sds_id) + istatus

!--- unscale sds

if (scaled > 0) then 

    sds_data_temp = sds_data
!---- linear
    if (scaled == 1) then
       sds_data_temp = min(1.0,max(0.0,real(sds_data_temp - scaled_min)/real(scaled_max - scaled_min)))
       sds_data_temp = unscaled_min + sds_data_temp * (unscaled_max - unscaled_min)
    endif

!---- log10
    if (scaled == 2) then
       sds_data_temp = min(1.0,max(0.0,real(sds_data_temp - scaled_min)/real(scaled_max - scaled_min)))
       sds_data_temp = unscaled_min + sds_data_temp * (unscaled_max - unscaled_min)
       sds_data_temp = 10**(sds_data_temp)
     endif

!---- square root
    if (scaled == 3) then
       sds_data_temp = min(1.0,max(0.0,real(sds_data_temp - scaled_min)/real(scaled_max - scaled_min)))
       sds_data_temp = unscaled_min + (sds_data_temp**2) * (unscaled_max - unscaled_min)
    endif

!--- set scaled missing values
    where (sds_data == scaled_missing)
         sds_data_temp = unscaled_missing
    endwhere

    sds_data = sds_data_temp

 endif
 
end subroutine READ_CLAVRX_HDF4_SDS_RANK1



end module HDF_PARAMS
