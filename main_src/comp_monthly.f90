! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: comp_monthly.f90 (src)
!       COMPILE_MONTHLY (program)
!
! PURPOSE: This code computes monthly files from daily files. 
!
! DESCRIPTION: For each parameter in the daily files, it computes
!              - the mean
!              - the standard deviation
!              - the maximum value
!              - the minimum value
!              - the number of daily values used 
!
!              This code takes input from a comp_monthly_input file that is 
!              composed of the following elements:
!              month  (1-12)
!              year (1981-200X)
!              node - asc, des, all
!              orbit - aft,mor,mid,all,two
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
! 
! REVISION HISTORY:
!   created August 2006
!
! NOTES:
!  1- the attributes for the first relevant file are used for the monthly files.
!     For files with multiple satellites, then calibration attributes are not
!     meaningful - maybe set to missing?  What does this mean for radiance?
!
! parameters
!  n_files_max - the maximum number of files that should be used to construct
!                a monthly average.  For 3 satellites with 2 files per day for
!                31 days - that gives about 200.
! min_count - minimum number of values required to make a monthly average
! satzen_thresh - maximum sensor zenith angle allowed in values used for monthly stats.
!
! monthly_type
!  1 - afternoon ascending     aft_asc
!  2 - afternoon descending    aft_des
!  3 - afternoon combined      aft
!  4 - morning ascending       mor_asc
!  5 - morning descending      mor_des
!  6 - morning combined        mor 
!  7 - mid-morning ascending       mid_asc
!  8 - mid-morning descending      mid_des
!  9 - mid-morning combined        mid 
! 10 - spare ascending       nxx_asc
! 11 - spare descending      nxx_des
! 12 - spare combined        nxx 
! 13 - all
! 14 - afternoon + morning   ampm
! 15 - morning + mid + afternoon  ammidpm
!
! Note, some hdf libraries limit the number of open files to be 32
! this version, does keep the files open
!--------------------------------------------------------------------------------------
program COMPILE_MONTHLY
 use CONSTANTS
 use HDF
 use HDF_PARAMS
 use NUMERICAL_ROUTINES
 use SCALING_PARAMETERS
 implicit none

!--- hard-wired parameters
 integer, parameter:: n_files_max = 200
 integer, parameter:: min_count_thresh = 1
 real(kind=real4):: satzen_thresh = 90.0
 integer, dimension(1:12,1981:2010):: aft_sat,mor_sat,mid_sat,spare_sat

!---
 character(len=100), dimension(n_files_max):: file_in
 integer, dimension(n_files_max):: file_mask
 integer:: sd_id
 integer:: sd_id_first
 integer:: file_id_first

 integer:: imonth, iyear,monthly_type,jday_start,jday_end,ileap, &
           iskip_stats,imean_write,istd_write,imin_write,imax_write,icount_write
 integer:: aft_sat_id, mor_sat_id, mid_sat_id,spare_sat_id
 integer:: sds_data_type_max, sds_data_type_min,sds_data_type_std,sds_data_type_count
 character(len=3):: aft_sat_string, mor_sat_string, mid_sat_string, spare_sat_string
 character(len=7):: aft_sat_dir_string,mor_sat_dir_string,mid_sat_dir_string,spare_sat_dir_string
 character(len=2):: temp_string
 real:: unscaled_min_std,unscaled_max_std
 character(len=7):: monthly_type_string
 real:: grid_resolution_input
 integer:: grid_format_flag_input,grid_format_flag_comp
 integer:: compress_flag

!----
 character(len=100)::  dir_in, dir_out, file_out
 integer:: ios
 integer(kind=int1):: asc_des_node  !(0=asc,1=des)
 integer(kind=int2):: jday

 character(len=4):: year_string
 character(len=3):: jday_string,month_string, &
                    grid_resolution_string
 character(len=1):: grid_format_string

!--- global attributes for composite
 integer(kind=int4):: sd_id_comp
 integer(kind=int4):: num_cells_comp
 integer:: therm_cal_1b_comp,Ref_cal_1b_comp,nav_flag_comp,use_sst_anal_comp, &
                        sst_anal_opt_comp, modis_clr_alb_flag_comp, nwp_flag_comp
 integer(kind=int4)::  start_time_comp,end_time_comp
 integer(kind=int4)::  acha_mode_comp, dcomp_mode_comp, wmo_sc_code_comp
 integer(kind=int2)::  start_year_comp,end_year_comp,start_day_comp,end_day_comp
 real(kind=real4)::  grid_resolution_comp, ch1_gain_low_comp, ch1_gain_high_comp, &
       ch2_gain_low_comp, ch2_gain_high_comp, &
       ch3a_gain_low_comp, ch3a_gain_high_comp, &
       ch1_switch_count_comp, ch1_dark_count_comp, &
       ch2_switch_count_comp, ch2_dark_count_comp, &
       ch3a_switch_count_comp, ch3a_dark_count_comp, &
       sun_earth_distance_comp, &
       c1_comp, c2_comp, a_3b_comp, b_3b_comp, nu_3b_comp, &
       a_31_comp, b_31_comp, nu_31_comp, a_32_comp, b_32_comp, nu_32_comp, &
       solar_3b_nu_comp,timerr_seconds_comp
       
 character(len=100):: file_name_comp, file_1b_comp, creator_comp, &
                    plang_comp, hdf_ver_comp, hdf_timestamp_comp, data_type_comp,grid_format_comp
 character(len=20):: platform_name_attribute_comp
 character(len=20):: platform_name_attribute_comp
 character(len=20):: sensor_name_attribute_comp
 character(len=120):: dark_name_comp
 character(len=120):: mask_name_comp

!--- sds attributes
  character(len=100):: sds_name
  real(kind=real4)::  resolution_km
  real:: unscaled_min, unscaled_max, unscaled_missing, unscaled_missing_count
  integer:: scaled_min, scaled_max, scaled_missing
  integer(kind=int1):: scaled,scaled_count
  character(len=100):: sds_units
  integer:: sds_data_type
  integer:: sds_rank, num_attrs
  integer, dimension(1):: dimsizes


!--- local
  integer:: ifile,n_files, i, istatus, sds_index, sds_id, istatus_write, &
            num_datasets, num_global_attrs,grid_format_flag
  real, dimension(:), allocatable:: z_daily,z_mean,z_std, &
                                    z_max,z_min,z_count,z_first
  real, dimension(:,:), allocatable:: satzen_daily

!--- hdf
  integer:: sfstart, sfend, sfginfo,sfselect, sffinfo, sfsnatt


!--- set compression flag
compress_flag = 1

!----
aft_sat = 0
mor_sat = 0
mid_sat = 0
spare_sat = 0
!                     1   2   3   4   5   6   7   8   9  10  11  12
aft_sat(:,1981) = (/  0,  0,  0,  0,  0,  0,  7,  7,  7,  7,  7,  7 /)
aft_sat(:,1982) = (/  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7 /)
aft_sat(:,1983) = (/  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7 /)
aft_sat(:,1984) = (/  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7 /)
aft_sat(:,1985) = (/  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9 /)
aft_sat(:,1986) = (/  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9 /)
aft_sat(:,1987) = (/  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9 /)
aft_sat(:,1988) = (/  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 11, 11 /)
aft_sat(:,1989) = (/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 /)
aft_sat(:,1990) = (/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 /)
aft_sat(:,1991) = (/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 /)
aft_sat(:,1992) = (/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 /)
aft_sat(:,1993) = (/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 /)
aft_sat(:,1994) = (/ 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 /)
aft_sat(:,1995) = (/ 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 /)
aft_sat(:,1996) = (/ 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 /)
aft_sat(:,1997) = (/ 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 /)
aft_sat(:,1998) = (/ 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 /)
aft_sat(:,1999) = (/ 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 /)
aft_sat(:,2000) = (/ 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 /)
aft_sat(:,2001) = (/ 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 /)
aft_sat(:,2002) = (/ 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 /)
aft_sat(:,2003) = (/ 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 /)
aft_sat(:,2004) = (/ 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 /)
aft_sat(:,2005) = (/ 16, 16, 16, 16, 16, 18, 18, 18, 18, 18, 18, 18 /)
aft_sat(:,2006) = (/ 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18 /)
aft_sat(:,2007) = (/ 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18 /)
aft_sat(:,2008) = (/ 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18 /)

!                     1   2   3   4   5   6   7   8   9  10  11  12
mor_sat(:,1981) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1982) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1983) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1984) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1985) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1986) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1987) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1988) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1989) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1990) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /)
mor_sat(:,1991) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0, 12, 12, 12 /)
mor_sat(:,1992) = (/ 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 /)
mor_sat(:,1993) = (/ 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 /)
mor_sat(:,1994) = (/ 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 /)
mor_sat(:,1995) = (/ 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 /)
mor_sat(:,1996) = (/ 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 /)
mor_sat(:,1997) = (/ 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 /)
mor_sat(:,1998) = (/ 12, 12, 12, 12, 12, 12, 12, 12, 12, 15, 15, 15 /)
mor_sat(:,1999) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)
mor_sat(:,2000) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)
mor_sat(:,2001) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)
mor_sat(:,2002) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)
mor_sat(:,2003) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)
mor_sat(:,2004) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)
mor_sat(:,2005) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)
mor_sat(:,2006) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)
mor_sat(:,2007) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)
mor_sat(:,2008) = (/ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 /)

!-----------------------------------------------------
! open input file
!----------------------------------------------------
open(unit=8,file="comp_monthly_input",status="old",action="read",position="rewind",iostat=ios)
if (ios /= 0) then
  print *, "error opening comp_monthly control file, ios = ", ios
  stop 1
else
  read(unit=8,fmt="(A)") dir_in
  read(unit=8,fmt="(A)") dir_out
  read(unit=8,fmt=*) imonth
  read(unit=8,fmt=*) iyear
  read(unit=8,fmt=*) monthly_type 
  read(unit=8,fmt=*) grid_resolution_input
  read(unit=8,fmt=*) grid_format_flag_input
endif
close(unit=8)

!-----------------------------------------------------------
! make a list of potential files that would comprise this file
!------------------------------------------------------------

!--- given the month, determine julian day range
   if (imonth == 1) then
      jday_start = 1
      jday_end = 31
      month_string = 'jan'
   elseif (imonth == 2) then
      jday_start = 32
      jday_end = 59
      month_string = 'feb'
   elseif (imonth == 3) then
      jday_start = 60 
      jday_end = 90
      month_string = 'mar'
   elseif (imonth == 4) then
      jday_start = 91
      jday_end = 120
      month_string = 'apr'
   elseif (imonth == 5) then
      jday_start = 121
      jday_end = 151
      month_string = 'may'
   elseif (imonth == 6) then
      jday_start = 152
      jday_end = 181
      month_string = 'jun'
   elseif (imonth == 7) then
      jday_start = 182
      jday_end = 212
      month_string = 'jul'
   elseif (imonth == 8) then
      jday_start = 213
      jday_end = 243
      month_string = 'aug'
   elseif (imonth == 9) then
      jday_start = 244
      jday_end = 273
      month_string = 'sep'
   elseif (imonth == 10) then
      jday_start = 274
      jday_end = 304
      month_string = 'oct'
   elseif (imonth == 11) then
      jday_start = 305
      jday_end = 334
      month_string = 'nov'
   elseif (imonth == 12) then
      jday_start = 335
      jday_end = 365
      month_string = 'dec'
   endif

!--- account for leap years
   ileap = 0
   ileap = leap_year_fct(iyear)

   if ((ileap == 1) .and. (imonth > 3)) then
           jday_start = jday_start + 1
           jday_end = jday_end + 1
   endif

   if ((ileap == 1) .and. (imonth == 2)) then
           jday_end = jday_end + 1
   endif


!--- determine which files to process for this month
    write (year_string,  '(I4.4)') iyear
    write (grid_format_string,   '(I1.1)') grid_format_flag_input

!--- make grid resolution string
 if (grid_resolution_input == 0.25) then
     grid_resolution_string = "025"
 elseif (grid_resolution_input == 0.5) then
     grid_resolution_string = "05"
 elseif (grid_resolution_input == 1.0) then
     grid_resolution_string = "1"
 elseif (grid_resolution_input == 2.5) then
     grid_resolution_string = "25"
 else
   print *, "invalid grid resolution, stopping: ", grid_resolution_input
   stop 3
 endif

!---- assign strings for filenames for each value of monthly_type
 if (monthly_type == 1) then
    monthly_type_string = "aft_asc"
 elseif (monthly_type == 2) then
    monthly_type_string = "aft_des"
 elseif (monthly_type == 3) then
    monthly_type_string = "aft"
 elseif (monthly_type == 4) then
    monthly_type_string = "mor_asc"
 elseif (monthly_type == 5) then
    monthly_type_string = "mor_des"
 elseif (monthly_type == 6) then
    monthly_type_string = "mor"
 elseif (monthly_type == 14) then
    monthly_type_string = "aft_mor"
 endif

!--- determine satellite id's that are active for this time
 aft_sat_id = aft_sat(imonth,iyear)
 mor_sat_id = mor_sat(imonth,iyear)
 mid_sat_id = mid_sat(imonth,iyear)
 spare_sat_id = spare_sat(imonth,iyear)

!--- make strings holding satellite names
    write (temp_string,  '(I2.2)') aft_sat_id
    aft_sat_string = "n"//temp_string
    aft_sat_dir_string = "noaa-"//temp_string
    write (temp_string,  '(I2.2)') mor_sat_id
    mor_sat_string = "n"//temp_string
    mor_sat_dir_string = "noaa-"//temp_string
    write (temp_string,  '(I2.2)') mid_sat_id
    mid_sat_string = "n"//temp_string
    mid_sat_dir_string = "noaa-"//temp_string
    write (temp_string,  '(I2.2)') spare_sat_id
    spare_sat_string = "n"//temp_string
    spare_sat_dir_string = "noaa-"//temp_string

   ifile = 0
   do jday = jday_start, jday_end

    write (jday_string,  '(I3.3)')jday 

!-- afternoon
    if (aft_sat_id > 0) then
!-- asc
      if (monthly_type == 1 .or. monthly_type == 3 .or. monthly_type == 13 .or.  &
          monthly_type == 14 .or. monthly_type == 15) then
          ifile = ifile + 1
          file_in(ifile) = "clavrx_"//trim(aft_sat_string)//"_asc_"// &
               trim(grid_resolution_string)//"_"//trim(grid_format_string)//"_"// &
               trim(year_string)//"_"//trim(jday_string)// &
               ".level3.hdf"
!         file_in(ifile) = aft_sat_dir_string//"/"//year_string//"/cell/"//trim(file_in(ifile))
      endif
!-- des
      if (monthly_type == 2 .or. monthly_type == 3 .or. monthly_type == 13 .or.  &
          monthly_type == 14 .or. monthly_type == 15) then
          ifile = ifile + 1
          file_in(ifile) = "clavrx_"//trim(aft_sat_string)//"_des_"// &
               trim(grid_resolution_string)//"_"//trim(grid_format_string)//"_"// &
               trim(year_string)//"_"//trim(jday_string)// &
               ".level3.hdf"
!         file_in(ifile) = aft_sat_dir_string//"/"//year_string//"/cell/"//trim(file_in(ifile))
      endif
     endif

!-- morning
    if (mor_sat_id > 0) then
!-- asc
      if (monthly_type == 4 .or. monthly_type == 6 .or. monthly_type == 13 .or.  &
          monthly_type == 14 .or. monthly_type == 15) then
          ifile = ifile + 1
          file_in(ifile) = "clavrx_"//trim(mor_sat_string)//"_asc_"// &
               trim(grid_resolution_string)//"_"//trim(grid_format_string)//"_"// &
               trim(year_string)//"_"//trim(jday_string)// &
               ".level3.hdf"
!         file_in(ifile) = mor_sat_dir_string//"/"//year_string//"/cell/"//trim(file_in(ifile))

      endif
!-- des
      if (monthly_type == 5 .or. monthly_type == 6 .or. monthly_type == 13 .or.  &
          monthly_type == 14 .or. monthly_type == 15) then
          ifile = ifile + 1
          file_in(ifile) = "clavrx_"//trim(mor_sat_string)//"_des_"// &
               trim(grid_resolution_string)//"_"//trim(grid_format_string)//"_"// &
               trim(year_string)//"_"//trim(jday_string)// &
               ".level3.hdf"
!         file_in(ifile) = mor_sat_dir_string//"/"//year_string//"/cell/"//trim(file_in(ifile))
      endif
     endif

   enddo
n_files = ifile
file_mask = 1


!--- test for which files are valid
  file_id_first = 0
  do ifile = 1, n_files
!--- heidinger convention

!--- open file
    sd_id = sfstart(trim(dir_in)//trim(file_in(ifile)), DFACC_READ)

!--- close file
    istatus = sfend(sd_id)   

!--- examine success
    if (sd_id == FAIL) then
      print *, "failed open for read on orbital gridded file, file = ",  &
                trim(dir_in)//trim(file_in(ifile))
      file_mask(ifile) = 0   !skip this file in further processing
    endif

!-- store first valid file id
    if ((file_mask(ifile) == 1) .and. (file_id_first == 0)) then
      file_id_first = ifile
    endif
  enddo

!---- check to if file found
  if (file_id_first == 0) then
     print *, "no valid files found, stopping "
     stop
  endif

!-------------------------------------------------------------------------------
!--- read global attributes for first valid file - these will serve as those
!--- for the final composite(except for time)
!-------------------------------------------------------------------------------
grid_format_comp=" "
sd_id_first = sfstart(trim(dir_in)//trim(file_in(file_id_first)), DFACC_READ)
call READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES(sd_id_first,data_type_comp,file_name_comp,file_1b_comp, &
                           resolution_km, &
                           start_year_comp,end_year_comp,start_day_comp,end_day_comp,start_time_comp,end_time_comp,&
                           num_cells_comp,num_cells_comp,grid_format_comp,grid_resolution_comp, &
                           therm_cal_1b_comp,Ref_cal_1b_comp,nav_flag_comp,use_sst_anal_comp,sst_anal_opt_comp, &
                           modis_clr_alb_flag_comp, nwp_flag_comp, ch1_gain_low_comp, ch1_gain_high_comp, &
                           ch1_switch_count_comp, ch1_dark_count_comp, &
                           ch2_gain_low_comp, ch2_gain_high_comp, &
                           ch2_switch_count_comp, ch2_dark_count_comp, &
                           ch3a_gain_low_comp, ch3a_gain_high_comp, &
                           ch3a_switch_count_comp, ch3a_dark_count_comp, &
                           sun_earth_distance_comp, &
                           c1_comp, c2_comp, a_3b_comp, b_3b_comp, nu_3b_comp, &
                           a_31_comp, b_31_comp, nu_31_comp, a_32_comp, b_32_comp, nu_32_comp, &
                           solar_3b_nu_comp, timerr_seconds_comp, &
                           acha_mode_comp, dcomp_mode_comp,  &
                           wmo_sc_code_comp, platform_name_attribute_comp,sensor_name_attribute_comp, &
                           dark_name_comp, mask_name_comp, &
                           creator_comp, plang_comp, hdf_ver_comp,hdf_timestamp_comp )

!--- do some checking for consistency with input
     if (grid_resolution_input /= grid_resolution_comp) then
        print *, "unexpected grid resolution, stopping ", grid_resolution_input, grid_resolution_comp
        stop
     endif

!--- determine number of sds's in these files (assume the same)
     istatus = sffinfo(sd_id,num_datasets,num_global_attrs)

!--- allocate space for stats
     allocate(z_daily(num_cells_comp),z_mean(num_cells_comp),z_std(num_cells_comp),z_max(num_cells_comp), &
              z_min(num_cells_comp), z_count(num_cells_comp),satzen_daily(num_cells_comp,n_files_max), &
              z_first(num_cells_comp))
     z_daily = 0.0
     z_mean = 0.0
     z_std = 0.0
     z_max = 0.0
     z_min = 0.0
     z_count = 0.0
     satzen_daily = 0.0
     z_first = 0.0

!--- convert grid_format into a flag for the output filename
    if (grid_format_comp(1:10) == "EQUAL_AREA") then
        grid_format_flag_comp = 0
    else
        grid_format_flag_comp = 1
    endif
     if (grid_format_flag_input /= grid_format_flag_comp) then
        print *, "unexpected format flag, stopping"
        stop
     endif

!--- make string needed for output filename
    write (year_string,  '(I4.4)') iyear
    write (jday_string,   '(I3.3)') jday
    write (grid_format_string,   '(I1.1)') grid_format_flag

!--- create time attributes for composite file
 start_year_comp = iyear
 end_year_comp = iyear
 start_day_comp = jday_start
 end_day_comp = jday_end
 start_time_comp = 0.0 * 3600000.0   !time in ms
 end_time_comp = 24.0 * 3600000.0    !time in ms

!--- construct output file name for this node
    file_out = "patmosx_"//trim(monthly_type_string)//"_"// &
               trim(grid_resolution_string)//"_"//trim(grid_format_string)//"_"// &
               trim(year_string)//"_"//trim(month_string)// &
               ".level3.hdf"


   print *, "Output written to = ", trim(file_out)

!--- open file for writing
    sd_id_comp = sfstart(trim(dir_out)//trim(file_out),DFACC_CREATE)

     if (sd_id_comp <  0) then
       print *, "composite HDF file creation failed. Exiting...", asc_des_node
       stop 38
     endif

!--- write global attributes to output file
 call WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES(sd_id_comp,"GRIDCELL",trim(file_out),trim(file_1b_comp), &
         resolution_km, &
         start_year_comp,end_year_comp,start_day_comp,end_day_comp,start_time_comp,end_time_comp,&
         num_cells_comp,num_cells_comp,grid_format_comp,grid_resolution_comp, &
         therm_cal_1b_comp,Ref_cal_1b_comp,nav_flag_comp,use_sst_anal_comp,sst_anal_opt_comp, &
         modis_clr_alb_flag_comp, nwp_flag_comp, ch1_gain_low_comp, ch1_gain_high_comp, &
         ch1_switch_count_comp, ch1_dark_count_comp, &
         ch2_gain_low_comp, ch2_gain_high_comp, &
         ch2_switch_count_comp, ch2_dark_count_comp, &
         ch3a_gain_low_comp, ch3a_gain_high_comp, &
         ch3a_switch_count_comp, ch3a_dark_count_comp, &
         sun_earth_distance_comp, &
         c1_comp, c2_comp, a_3b_comp, b_3b_comp, nu_3b_comp, &
         a_31_comp, b_31_comp, nu_31_comp, a_32_comp, b_32_comp, nu_32_comp, &
         solar_3b_nu_comp,timerr_seconds_comp, &
         acha_mode_comp, dcomp_mode_comp,wmo_sc_code_comp, &
         platform_name_attribute_comp,sensor_name_attribute_comp, &
         dark_name_comp, mask_name_comp)

!---- write additional attributes for a monthly file
istatus = sfsnatt(sd_id_comp, "JDAY_MIN", DFNT_FLOAT32, 1, real(jday_start,kind=real4))+istatus
istatus = sfsnatt(sd_id_comp, "JDAY_MAX", DFNT_FLOAT32, 1, real(jday_end,kind=real4))+istatus
istatus = sfsnatt(sd_id_comp, "JDAY_MAX", DFNT_FLOAT32, 1, real(jday_end,kind=real4))+istatus
istatus = sfsnatt(sd_id_comp, "AFTERNOON_SATELLITE_ID", DFNT_INT32, 1, aft_sat_id)+istatus
istatus = sfsnatt(sd_id_comp, "MORNING_SATELLITE_ID", DFNT_INT32, 1, mor_sat_id)+istatus
istatus = sfsnatt(sd_id_comp, "MID-MORNING_SATELLITE_ID", DFNT_INT32, 1, mid_sat_id)+istatus
istatus = sfsnatt(sd_id_comp, "SPARE_SATELLITE_ID", DFNT_INT32, 1, spare_sat_id)+istatus


!--- set scaling parameters for stats here
            sds_data_type_max = DFNT_INT8
            sds_data_type_min = DFNT_INT8
            sds_data_type_std = DFNT_INT8
            sds_data_type_count = DFNT_INT8
            scaled_count = 0
            unscaled_missing_count =  0.0

!--- loop through each parameter 

   istatus = 0
   sds_index_loop: do sds_index = 0,num_datasets-1


!-- determine name of this sds
         sds_id = sfselect(sd_id_first,sds_index)
         sds_name = " "
         istatus = sfginfo(sds_id, sds_name, sds_rank, dimsizes, sds_data_type, num_attrs) + istatus
         if (istatus /= 0) then
            print *, "could not find sds ", sds_index
            exit
         endif

   
!--- skip printing statistcs for sds's that are temporally uniform
         iskip_stats = 0
         if ((trim(sds_name) == "cell_latitude") .or. &
             (trim(sds_name) == "cell_longitude") .or. &
             (trim(sds_name) == "surface_type") .or. &
             (trim(sds_name) == "cell_index") .or. &
             (trim(sds_name) == "year") .or. &
             (trim(sds_name) == "spacecraft_id") .or. &
             (trim(sds_name) == "surface_elevation")) then
             iskip_stats = 1 
         endif

!--- control output of sds_stats

         imean_write = 1 
         istd_write = 1 
         icount_write = 1 
         imin_write = 1 
         imax_write = 1 

!--- skip sds's that have no meaning for multi-day files
         if ((trim(sds_name) == "asc_des_flag") .or. &
             (trim(sds_name) == "ch3a_on_flag") .or. &
             (trim(sds_name) == "orbit_number") .or. &
             (trim(sds_name) == "cld_type") .or. &
             (trim(sds_name) == "jday")) then
             imean_write = 0 * imean_write
             istd_write = 0 * istd_write
             icount_write = 0 * icount_write
             imin_write = 0 * imin_write 
             imax_write = 0  * imax_write
         endif

!--- skip sds's that have no meaning for multi-satellite files
         if (monthly_type > 12) then
            if (trim(sds_name) == "spacecraft_id") then
             imean_write = 0 * imean_write
             istd_write = 0 * istd_write
             icount_write = 0 * icount_write
             imin_write = 0 * imin_write 
             imax_write = 0  * imax_write
            endif
         endif

!--- skip sds's that are standard deviations (have std in sds_name)
        if (index(trim(sds_name),"std") > 0) then
             imean_write = 1 * imean_write
             istd_write = 0 * istd_write
             icount_write = 0 * icount_write
             imin_write = 0 * imin_write 
             imax_write = 0  * imax_write
        endif
       
!-- reset composite for this sds
      file_loop_2: do ifile = 1, n_files

!--- skip files that are not valid
         if (file_mask(ifile) == 0) then
            cycle
         endif

!--- output progress to standard output
         if (sds_index == 1) then
          print *, "compositing ", trim(file_in(ifile))
         endif

!--- for first valid file, initialize stats appropriately
         if (ifile == file_id_first) then
          z_mean = 0.0
          z_std = 0.0
          z_max = -1.0*huge(1.0)
          z_min = 1.0*huge(1.0)
          z_count = 0
         endif

!--- read in sensor zenith for first sds and store it
        if (sds_index == 1) then
         sd_id = sfstart(trim(dir_in)//trim(file_in(ifile)), DFACC_READ)
         call READ_CLAVRX_HDF4_SDS_RANK1(sd_id,num_cells_comp, &
                              trim("sensor_zenith"),satzen_daily(:,ifile), &
                              sds_data_type,scaled,sds_units, &
                              unscaled_min,unscaled_max,unscaled_missing,&
                              scaled_min,scaled_max,scaled_missing,&
                              istatus)
         istatus = sfend(sd_id)   
       endif

!--- read daily sds 
         sd_id = sfstart(trim(dir_in)//trim(file_in(ifile)), DFACC_READ)
         call READ_CLAVRX_HDF4_SDS_RANK1(sd_id,num_cells_comp, &
                              trim(sds_name),z_daily, &
                              sds_data_type,scaled,sds_units, &
                              unscaled_min,unscaled_max,unscaled_missing,&
                              scaled_min,scaled_max,scaled_missing,&
                              istatus)
         istatus = sfend(sd_id)   

!---- for first file, reset first time for each cell to missing
            if ((index(trim(sds_name),"utc_time") > 0).and.(ifile == file_id_first)) then
                 z_first = missing_value_real4
            endif

             do i = 1, num_cells_comp

!---  adjust time variables if straddling midnight
               if ((index(trim(sds_name),"utc_time") > 0).and. &
                   (z_daily(i) /= unscaled_missing)) then

!--- save first valid occurence
                   if (z_first(i) == unscaled_missing) then
                       z_first(i) = z_daily(i)
                   endif
                   if ((z_daily(i) - z_first(i) > 12.0).and.(z_first(i) /= unscaled_missing)) then
                       z_daily(i) = z_daily(i) - 24.0
                   endif
                   if ((z_daily(i) - z_first(i) < -12.0).and.(z_first(i) /= unscaled_missing)) then
                       z_daily(i) = z_daily(i) + 24.0
                   endif

               endif 

!--- update stats
               if (z_daily(i) /= unscaled_missing) then
                  if (satzen_daily(i,ifile) < satzen_thresh) then
                    z_count(i) = z_count(i) + 1
                    z_mean(i) = z_mean(i) + z_daily(i)
                    z_std(i) = z_std(i) + z_daily(i)**2
                    if (z_daily(i) > z_max(i)) then
                        z_max(i) = z_daily(i)
                    endif
                    if (z_daily(i) < z_min(i)) then
                        z_min(i) = z_daily(i)
                    endif
                 endif
               endif
              end do
      end do file_loop_2

!---- compute final stats and fill in missing where appropriate
      do i = 1, num_cells_comp
          if (z_count(i) >= min_count_thresh) then
                  z_mean(i) = z_mean(i) / z_count(i)
                  z_std(i) = sqrt(max(0.0,(z_std(i) / z_count(i) - z_mean(i)**2)))

!---  adjust time variables if necessary
                  if (index(trim(sds_name),"time") > 0) then
                    if (z_mean(i) > 24.0) then
                        z_mean(i) = z_mean(i) - 24.0
                    endif
                    if (z_mean(i) < 0.0) then
                        z_mean(i) = z_mean(i) + 24.0
                    endif
                  endif 
          else
            if (sds_data_type == DFNT_INT8) then
                  z_mean(i) = missing_value_real4
            endif
            if (sds_data_type == DFNT_INT16) then
                  z_mean(i) = missing_value_real4
            endif
            if (sds_data_type == DFNT_INT32) then
                  z_mean(i) = missing_value_real4
            endif
            if (sds_data_type == DFNT_FLOAT32) then
                  z_mean(i) = missing_value_real4
            endif
            z_count(i) = unscaled_missing_count
            z_max(i) = unscaled_missing
            z_min(i) = unscaled_missing
            z_std(i) = unscaled_missing

          endif
       enddo

!------------------------------------------------------
!--- write out monthly sds's to output file
!------------------------------------------------------

  if (iskip_stats == 0) then
!--- mean - use same scaling as in daily file
     if (imean_write == 1) then
         call WRITE_CLAVRX_HDF4_SDS(sd_id_comp,real(z_mean,kind=real4), &
                  trim(sds_name)//"_mean",&
                  sds_data_type,scaled,unscaled_min,unscaled_max,trim(sds_units), unscaled_missing, &
                  "grid_cell_number", compress_flag, istatus_write)
     endif

!--- max - use same scaling as in daily file (maybe switch to int1?)
      if (imax_write == 1) then
         call WRITE_CLAVRX_HDF4_SDS(sd_id_comp,real(z_max,kind=real4), &
                  trim(sds_name)//"_max",&
                  sds_data_type_max,scaled,unscaled_min,unscaled_max,trim(sds_units), unscaled_missing, &
                  "grid_cell_number", compress_flag, istatus_write)
      endif

!--- min - use same scaling as in daily file
      if (imin_write == 1) then
         call WRITE_CLAVRX_HDF4_SDS(sd_id_comp,real(z_min,kind=real4), &
                  trim(sds_name)//"_min",&
                  sds_data_type_min,scaled,unscaled_min,unscaled_max,trim(sds_units), unscaled_missing, &
                  "grid_cell_number", compress_flag, istatus_write)
      endif

!--- standard deviation - note change in scaling
      if (istd_write == 1) then
         unscaled_min_std = 0.0
         unscaled_max_std = maxval(z_std)
         call WRITE_CLAVRX_HDF4_SDS(sd_id_comp,real(z_std,kind=real4), &
                  trim(sds_name)//"_std",&
                  sds_data_type_std,scaled,unscaled_min_std,unscaled_max_std,trim(sds_units), unscaled_missing, &
                  "grid_cell_number", compress_flag, istatus_write)
     endif 

!--- number of values used in compilation - note change in scaling
     if (icount_write == 1) then
         call WRITE_CLAVRX_HDF4_SDS(sd_id_comp,real(z_count,kind=real4), &
                  trim(sds_name)//"_count",&
                  sds_data_type_count,scaled_count,unscaled_min,unscaled_max,"none", unscaled_missing, &
                  "grid_cell_number", compress_flag, istatus_write)
     endif

    endif

    if (iskip_stats == 1) then   !just write out mean (do not append _mean to sds_name)
         call WRITE_CLAVRX_HDF4_SDS(sd_id_comp,real(z_mean,kind=real4), &
                  trim(sds_name),&
                  sds_data_type,scaled,unscaled_min,unscaled_max,trim(sds_units), unscaled_missing, &
                  "grid_cell_number", compress_flag, istatus_write)
    endif

   end do sds_index_loop 

!--- close file
   istatus = sfend(sd_id_comp)


!--- deallocate memory
  deallocate(z_daily,z_mean,z_std,z_max,z_min,z_count,satzen_daily,z_first)


end program COMPILE_MONTHLY
