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
!  UNSCALE_VECTOR_I1_RANK1
!  READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES
!  READ_CLAVRX_HDF4_SDS_RANK1
!
!--------------------------------------------------------------------------------------
module HDF_PARAMS

use CX_CONSTANTS_MOD, only: &
   int1 &
   , int2 &
   , int4 &
   , real4  
   
   implicit none
   include 'hdf.f90'   
private
public:: READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES, &
         READ_CLAVRX_HDF4_SDS_RANK1


character(len=256), save, public:: renav_data_from


!--- scaling options
integer(kind=int1), parameter, public:: NO_SCALING = 0
integer(kind=int1), parameter, public:: LINEAR_SCALING = 1
integer(kind=int1), parameter, public:: LOG10_SCALING = 2
integer(kind=int1), parameter, public:: SQUARE_ROOT_SCALING = 3

contains



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
