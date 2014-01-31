! $Id: comp_asc_des.f90,v 1.12.2.2 2014/01/26 04:48:32 heidinger Exp $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: comp_asc_des.f90 (src)
!       COMPILE_ASC_DES (program)
!
! PURPOSE: A main code generated one of the executables in the CLAVR-x
!          processing system
!
! DESCRIPTION: this code takes the level3 files created for each orbit
!              for one day from one satellite and writes separate level3 
!              files for the ascending and descending nodes. This code runs
!              after clavrxorb has processed the level-1b files for one day.
!              This program reads input (directories and filenames) from a file
!              called comp_asc_des_input.  Currently, there are no
!              command line arguments.  
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
! 
! (c) This code is copyrighted by the author and all NOAA restrictions apply
!
! Revision History:
!  created May 2004
!  August 2004 - Added error messaging (D. Donahue - STC)
!  June 2006 - Rewrote for version 4 CLAVR-x using standardized HDF read/write commands
!
! Format of Required Input File:
!    line 1: directory of orbital level3 files (input)
!    line 2: directory of daily level3 files (output)
!    line 3: year
!    line 4: julian day
!    line 5: satellite number (i.e. 18 = NOAA-18)
!    line 6: level-3 grid resolution (input files with different resolution are
!                                     skipped in the processing)
!    line 7: level-3 grid format.  0 = equal-area, 1 = equal-angle (input files
!                            with differnt format are skipped)
!    line 8+: orbital level3 file (one file per line)
!
! Reference:  Level3 file contents are given on
!             http://cimss.ssec.wisc.edu/patmosx
!
! Dependencies:  (The following are names of modules)
!    CONSTANTS
!    HDF
!    HDF_PARAMS
!    NUMERICAL_ROUTINES
!    SCALING_PARAMETERS
!
! Calling Sequence:  
!    comp_asc_des_level2 node nav_flag geo_flag
!--------------------------------------------------------------------------------------
program COMPILE_ASC_DES
 use CONSTANTS
 use HDF
 use HDF_PARAMS
 use NUMERICAL_ROUTINES
 use SCALING_PARAMETERS
 implicit none

!------
 integer, parameter:: n_files_max = 100
 character(len=100), dimension(n_files_max):: file_in
 integer, dimension(n_files_max):: sd_id

!----
 character(len=100)::  dir_in, dir_out, file_out, file_temp
 integer:: ios
 integer(kind=int1):: asc_des_node  !(0=asc,1=des)
 integer(kind=int2):: year
 integer(kind=int2):: jday
 integer(kind=int2):: wmo_id
 integer:: sd_id_first

 character(len=4):: year_string
 character(len=3):: jday_string,sat_string,node_string,grid_resolution_string
 character(len=1):: grid_format_string
 real(kind=real4):: grid_resolution_input
 integer(kind=int4):: grid_format_flag_input
 integer:: compress_flag

!--- global attributes for composite
 integer(kind=int4):: sd_id_comp
 integer(kind=int4):: num_cells_comp,num_cells_with_data_comp
 integer:: therm_cal_1b_comp,Ref_cal_1b_comp,nav_flag_comp,use_sst_anal_comp, &
                        sst_anal_opt_comp, modis_clr_alb_flag_comp, nwp_flag_comp
 integer(kind=int4)::  start_time_comp,end_time_comp
 integer(kind=int4)::  acha_mode_comp, dcomp_mode_comp, wmo_id_comp
 integer(kind=int2)::  start_year_comp,end_year_comp,start_day_comp,end_day_comp
 real(kind=real4)::  resolution_km
 real(kind=real4)::  grid_resolution_comp, ch1_gain_low_comp, ch1_gain_high_comp, &
       ch2_gain_low_comp, ch2_gain_high_comp, &
       ch3a_gain_low_comp, ch3a_gain_high_comp, &
       ch1_switch_count_comp, ch1_dark_count_comp, &
       ch2_switch_count_comp, ch2_dark_count_comp, &
       ch3a_switch_count_comp, ch3a_dark_count_comp, &
       sun_earth_distance_comp, &
       c1_comp, c2_comp, a_20_comp, b_20_comp, nu_20_comp, &
       a_31_comp, b_31_comp, nu_31_comp, a_32_comp, b_32_comp, nu_32_comp, &
       solar_Ch20_nu_comp,timerr_seconds_comp
                                                                                                                                                         
 character(len=100):: file_name_comp, file_1b_comp, creator_comp,   &
                      plang_comp, hdf_ver_comp, hdf_timestamp_comp, data_type_comp,grid_format_comp
 character(len=20):: sensor_name_attribute_comp
 character(len=120):: dark_name_comp
 character(len=120):: mask_name_comp

!--- global attributes for orbital files
 integer(kind=int4):: num_cells_orbit
 integer:: therm_cal_1b_orbit,Ref_cal_1b_orbit,nav_flag_orbit,use_sst_anal_orbit, &
                        sst_anal_opt_orbit, modis_clr_alb_flag_orbit, nwp_flag_orbit
 integer(kind=int4)::  start_time_orbit,end_time_orbit, wmo_id_temp_orbit
 integer(kind=int4)::  acha_mode_orbit,dcomp_mode_orbit
 integer(kind=int2)::  start_year_orbit,end_year_orbit,start_day_orbit,end_day_orbit
 real(kind=real4)::  grid_resolution_orbit, ch1_gain_low_orbit, ch1_gain_high_orbit, &
       ch2_gain_low_orbit, ch2_gain_high_orbit, &
       ch3a_gain_low_orbit, ch3a_gain_high_orbit, &
       ch1_switch_count_orbit, ch1_dark_count_orbit, &
       ch2_switch_count_orbit, ch2_dark_count_orbit, &
       ch3a_switch_count_orbit, ch3a_dark_count_orbit, &
       sun_earth_distance_orbit, &
       c1_orbit, c2_orbit, a_20_orbit, b_20_orbit, nu_20_orbit, &
       a_31_orbit, b_31_orbit, nu_31_orbit, a_32_orbit, b_32_orbit, nu_32_orbit, &
       solar_Ch20_nu_orbit, timerr_seconds_orbit
                                                                                                                                                         
 character(len=100):: file_name_orbit, file_1b_orbit, creator_orbit,  &
                    plang_orbit, hdf_ver_orbit, hdf_timestamp_orbit, data_type_orbit,grid_format_orbit
 character(len=20):: sensor_name_attribute_orbit
 character(len=120):: dark_name_orbit
 character(len=120):: mask_name_orbit

!--- sds attributes
  character(len=100):: sds_name
  real:: unscaled_min, unscaled_max, unscaled_missing
  integer:: scaled_min, scaled_max, scaled_missing
  integer(kind=int1):: scaled
  character(len=100):: sds_units
  integer:: sds_data_type
  integer:: sds_rank, num_attrs
  integer, dimension(1):: dimsizes


!--- local
  integer:: ifile,n_files, i, istatus, sds_index, sds_id, istatus_write, &
            num_datasets, num_global_attrs, iupdate,grid_format_flag
  integer, dimension(:), allocatable:: file_index_comp, file_mask
  integer, dimension(:), allocatable:: num_cells_with_data_orbit
  real, dimension(:), allocatable:: satzen_orbit, satzen_comp,z_comp,z_orbit,time_orbit,asc_des_orbit,wmo_id_orbit
  real, dimension(:,:), allocatable:: cell_index_orbit
  real:: time_thresh_1, time_thresh_2

!--- hdf
   integer:: sfstart, sfend, sfginfo,sfselect, sffinfo


!-- set time thresh
   time_thresh_1 = 2.0 !time in hours
   time_thresh_2 = 22.0 !time in hours
 
!-----------------------------------------------------
! open input file
!----------------------------------------------------
open(unit=8,file="comp_asc_des_input",status="old",action="read",position="rewind",iostat=ios)
if (ios /= 0) then
  print *, "error opening comp_asc_dec control file, ios = ", ios
  stop 1
else
  read(unit=8,fmt="(A)") dir_in
  read(unit=8,fmt="(A)") dir_out
  read(unit=8,fmt=*) year
  read(unit=8,fmt=*) jday
  read(unit=8,fmt=*) wmo_id
  read(unit=8,fmt=*) grid_resolution_input
  read(unit=8,fmt=*) grid_format_flag_input
endif


compress_flag = 1

!--- read in file names
  ifile = 0
  do i = 1, n_files_max
    read(unit=8,fmt=*,iostat=ios) file_temp
    if (ios /= 0) then
       exit
    endif
    ifile = ifile + 1
    file_in(ifile) = file_temp
  end do
  n_files = ifile


!--- open files for input
  do ifile = 1, n_files
    sd_id(ifile) = sfstart(trim(dir_in)//trim(file_in(ifile)), DFACC_READ)
    if (sd_id(ifile) == FAIL) then
      print *, "failed open for read on orbital gridded file, stopping"
      stop
    endif
  enddo

!--- read global attributes for first file - these will serve as those
!--- for the final composite(except for time)
grid_format_comp=" "
call READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES(sd_id(1),data_type_comp,file_name_comp,file_1b_comp, &
                           resolution_km, &
                           start_year_comp,end_year_comp,start_day_comp,end_day_comp,start_time_comp,end_time_comp,&
                           num_cells_comp,num_cells_with_data_comp,grid_format_comp,grid_resolution_comp, &
                           therm_cal_1b_comp,Ref_cal_1b_comp,nav_flag_comp,use_sst_anal_comp,sst_anal_opt_comp, &
                           modis_clr_alb_flag_comp, nwp_flag_comp, ch1_gain_low_comp, ch1_gain_high_comp, &
                           ch1_switch_count_comp, ch1_dark_count_comp, &
                           ch2_gain_low_comp, ch2_gain_high_comp, &
                           ch2_switch_count_comp, ch2_dark_count_comp, &
                           ch3a_gain_low_comp, ch3a_gain_high_comp, &
                           ch3a_switch_count_comp, ch3a_dark_count_comp, &
                           sun_earth_distance_comp, &
                           c1_comp, c2_comp, a_20_comp, b_20_comp, nu_20_comp, &
                           a_31_comp, b_31_comp, nu_31_comp, a_32_comp, b_32_comp, nu_32_comp, &
                           solar_Ch20_nu_comp, timerr_seconds_comp, &
                           acha_mode_comp, dcomp_mode_comp,  &
                           wmo_id_comp, sensor_name_attribute_comp, &
                           dark_name_comp, mask_name_comp, &
                           creator_comp, plang_comp, hdf_ver_comp,hdf_timestamp_comp )


!--- determine number of sds's in these files (assume the same)
     istatus = sffinfo(sd_id(1),num_datasets,num_global_attrs)

!--- allocate memory needed for composities
    allocate(file_index_comp(num_cells_comp),satzen_comp(num_cells_comp))
    file_index_comp = 0

    allocate(num_cells_with_data_orbit(n_files),file_mask(n_files))
    num_cells_with_data_orbit = 0
    file_mask = 0

!--- allocate array that holds cell indices for all orbits
     allocate(cell_index_orbit(num_cells_comp,n_files))  !this is saved and used again
     cell_index_orbit = 0.0

!--- allocate space for composited sds
     allocate(z_comp(num_cells_comp))
     z_comp = 0.0

!--- convert grid_format into a flag for the output filename
    if (grid_format_comp(1:10) == "EQUAL_AREA") then
        grid_format_flag = 0
    else
        grid_format_flag = 1
    endif


!--- check to see that file is consistent with expectations from input
    if (grid_format_flag /= grid_format_flag_input) then
       print *, "Unexpected grid format, stopping ", grid_format_flag_input, grid_format_flag
    endif
    if (grid_resolution_comp /= grid_resolution_input) then
       print *, "Unexpected grid resolution, stopping ", grid_resolution_input, grid_resolution_comp
    endif

!--- make string needed for output filename
    write (year_string,  '(I4.4)') year
    write (jday_string,   '(I3.3)') jday
    write (grid_format_string,   '(I1.1)') grid_format_flag


    select case(wmo_id)
      case(706)
        sat_string = "n06"
      case(707)
        sat_string = "n07"
      case(708)
        sat_string = "n05"
      case(200)
        sat_string = "n08"
      case(201)
        sat_string = "n09"
      case(202)
        sat_string = "n10"
      case(203)
        sat_string = "n11"
      case(204)
        sat_string = "n12"
      case(205)
        sat_string = "n14"
      case(206)
        sat_string = "n15"
      case(207)
        sat_string = "n16"
      case(208)
        sat_string = "n17"
      case(209)
        sat_string = "n18"
      case(223)
        sat_string = "n19"
      case default
        print *, "Unexpected WMO Spacecraft Id Number"
    end select

!--- create time attributes for composite file
 start_year_comp = year
 end_year_comp = year
 start_day_comp = jday
 end_day_comp = jday
 start_time_comp = 0.0 * 3600000.0   !time in ms
 end_time_comp = 24.0 * 3600000.0    !time in ms
 num_cells_with_data_comp = num_cells_comp
!file_1b_comp = ""

!--- make grid resolution string
 if (grid_resolution_comp == 0.25) then
     grid_resolution_string = "025"
 elseif (grid_resolution_comp == 0.5) then
     grid_resolution_string = "05"
 elseif (grid_resolution_comp == 1.0) then
     grid_resolution_string = "1"
 elseif (grid_resolution_comp == 2.5) then
     grid_resolution_string = "25"
 else
   print *, "invalid grid resolution, stopping: ", grid_resolution_comp
   stop 3
 endif

!---------------------------------------------------------------
! loop over asc / des
!---------------------------------------------------------------
asc_des_loop: do asc_des_node = 0, 1

 if (asc_des_node == 0) then
    node_string = "asc"
 else
    node_string = "des"
 endif


!--- initialize composite variables
     file_index_comp = 0
     cell_index_orbit = 0
     num_cells_with_data_orbit = 0
     satzen_comp = missing_value_real4 

!-----------------------------------------------------------------------------------------
!--- loop over files, reading satzen (the composite basis) and generate file index
!-----------------------------------------------------------------------------------------
     sd_id_first = 0
     do ifile = 1,n_files

!--- read global attributes of this file
             grid_format_orbit = " "
             call READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES(sd_id(ifile),data_type_orbit,file_name_orbit,file_1b_orbit, &
                           resolution_km, &
                           start_year_orbit,end_year_orbit,start_day_orbit,end_day_orbit,start_time_orbit,end_time_orbit,&
                           num_cells_orbit,num_cells_with_data_orbit(ifile),grid_format_orbit,grid_resolution_orbit, &
                           therm_cal_1b_orbit,Ref_cal_1b_orbit,nav_flag_orbit,use_sst_anal_orbit,sst_anal_opt_orbit, &
                           modis_clr_alb_flag_orbit,  &
                           nwp_flag_orbit, ch1_gain_low_orbit, ch1_gain_high_orbit, &
                           ch1_switch_count_orbit, ch1_dark_count_orbit, &
                           ch2_gain_low_orbit, ch2_gain_high_orbit, &
                           ch2_switch_count_orbit, ch2_dark_count_orbit, &
                           ch3a_gain_low_orbit, ch3a_gain_high_orbit, &
                           ch3a_switch_count_orbit, ch3a_dark_count_orbit, &
                           sun_earth_distance_orbit, &
                           c1_orbit, c2_orbit, a_20_orbit, b_20_orbit, nu_20_orbit, &
                           a_31_orbit, b_31_orbit, nu_31_orbit, a_32_orbit, b_32_orbit, nu_32_orbit, &
                           solar_Ch20_nu_orbit, timerr_seconds_orbit, &
                           acha_mode_orbit, dcomp_mode_orbit, &
                           wmo_id_temp_orbit, sensor_name_attribute_orbit, &
                           dark_name_orbit, mask_name_orbit, &
                           creator_orbit, plang_orbit, hdf_ver_orbit,hdf_timestamp_orbit )

!--- determine if this file fits our window
       file_mask(ifile) = 0

       if ((start_year_orbit == start_year_comp) .or. (start_year_orbit==end_year_comp) .or. &
           (end_year_orbit == start_year_comp) .or. (end_year_orbit==end_year_comp)) then

            if ((start_day_orbit == start_day_comp) .or. (start_day_orbit == end_day_orbit) .or. &
                (end_day_orbit == start_day_comp) .or. (end_day_orbit == end_day_orbit)) then
 
              file_mask(ifile) = 1

           endif

!-- store first valid file id
           if ((file_mask(ifile) == 1) .and. (sd_id_first == 0)) then
             sd_id_first = sd_id(ifile)
           endif

       endif

!--- check if file is of the right time
 
       if ((grid_resolution_comp /= grid_resolution_orbit) .or. & 
           (grid_format_comp(1:8) /= grid_format_orbit(1:8))) then
            file_mask(ifile) = 0
       endif

!--- if not in window, move to next orbit
       if (file_mask(ifile) == 0) then    !MAKE A FILE MASK
           cycle
       endif    

!--- allocate memory needed
          allocate(satzen_orbit(num_cells_with_data_orbit(ifile)))   !this is not saved
          allocate(time_orbit(num_cells_with_data_orbit(ifile)))   !this is not saved
          allocate(asc_des_orbit(num_cells_with_data_orbit(ifile)))   !this is not saved
          allocate(wmo_id_orbit(num_cells_with_data_orbit(ifile)))   !this is not saved
          satzen_orbit = 0.0
          time_orbit = 0.0
          asc_des_orbit = 0.0
          wmo_id_orbit = 0.0

!--- read cell_index
          call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile),"cell_index", &
                                  cell_index_orbit(1:num_cells_with_data_orbit(ifile),ifile), &
                                  sds_data_type,scaled,sds_units, &
                                  unscaled_min,unscaled_max,unscaled_missing,&
                                  scaled_min,scaled_max,scaled_missing,&
                                  istatus)

!--- read asc_des
          call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile),"asc_des_flag", &
                                  asc_des_orbit, &
                                  sds_data_type,scaled,sds_units, &
                                  unscaled_min,unscaled_max,unscaled_missing,&
                                  scaled_min,scaled_max,scaled_missing,&
                                  istatus)
!--- read time
          call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile),"utc_time", &
                                  time_orbit, &
                                  sds_data_type,scaled,sds_units, &
                                  unscaled_min,unscaled_max,unscaled_missing,&
                                  scaled_min,scaled_max,scaled_missing,&
                                  istatus)

!--- read spacecraft id
          call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile),"spacecraft_id", &
                                  wmo_id_orbit, &
                                  sds_data_type,scaled,sds_units, &
                                  unscaled_min,unscaled_max,unscaled_missing,&
                                  scaled_min,scaled_max,scaled_missing,&
                                  istatus)

!--- read satzen
            call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile),"sensor_zenith", &
                                  satzen_orbit, &
                                  sds_data_type,scaled,sds_units, &
                                  unscaled_min,unscaled_max,unscaled_missing,&
                                  scaled_min,scaled_max,scaled_missing,&
                                  istatus)

!--- make file index
            do i = 1, num_cells_with_data_orbit(ifile)

              iupdate = 1

!--- check for most nadir
                if ((satzen_orbit(i) >= satzen_comp(int(cell_index_orbit(i,ifile)))) .and. &
                    (satzen_comp(int(cell_index_orbit(i,ifile))) /= missing_value_real4)) then
                 iupdate = 0
                endif

!--- orbit straddles midnight, we want data before midnight and this is from next day
                if ((start_day_orbit < end_day_orbit) .and. (start_day_orbit == jday) .and. &
                  (time_orbit(i) < time_thresh_1)) then
                  iupdate = 0
                endif 

!--- orbit straddles midnight, we want data after midnight and this is from previous day
                if ((start_day_orbit < end_day_orbit) .and. (end_day_orbit == jday) .and. &
                  (time_orbit(i)) > time_thresh_2) then
                  iupdate = 0
                endif 

!--- wrong node
                if (asc_des_node /= asc_des_orbit(i)) then
                  iupdate = 0
                endif

!--- wrong satellite
              if (wmo_id /= wmo_id_orbit(i)) then
                  iupdate = 0
              endif

!--- choose this cell for composite if it met the tests
             if (iupdate == 1) then
                 file_index_comp(int(cell_index_orbit(i,ifile))) = ifile
                 satzen_comp(int(cell_index_orbit(i,ifile))) = satzen_orbit(i)
              endif

            enddo

!--- deallocate memory
     deallocate(satzen_orbit, time_orbit, asc_des_orbit, wmo_id_orbit)


     enddo


!--- construct output file name for this node
    file_out = "clavrx_"//trim(sat_string)//"_"//trim(node_string)//"_"// &
               trim(grid_resolution_string)//"_"//trim(grid_format_string)//"_"// &
               trim(year_string)//"_"//trim(jday_string)// &
               ".level3.hdf"

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
          num_cells_comp,num_cells_with_data_comp,grid_format_comp,grid_resolution_comp, &
          therm_cal_1b_comp,Ref_cal_1b_comp,nav_flag_comp,use_sst_anal_comp,sst_anal_opt_comp, &
          modis_clr_alb_flag_comp, nwp_flag_comp, ch1_gain_low_comp, ch1_gain_high_comp, &
          ch1_switch_count_comp, ch1_dark_count_comp, &
          ch2_gain_low_comp, ch2_gain_high_comp, &
          ch2_switch_count_comp, ch2_dark_count_comp, &
          ch3a_gain_low_comp, ch3a_gain_high_comp, &
          ch3a_switch_count_comp, ch3a_dark_count_comp, &
          sun_earth_distance_comp, &
          c1_comp, c2_comp, a_20_comp, b_20_comp, nu_20_comp, &
          a_31_comp, b_31_comp, nu_31_comp, a_32_comp, b_32_comp, nu_32_comp, &
          solar_Ch20_nu_comp,timerr_seconds_comp, &
          acha_mode_comp, dcomp_mode_comp,&
          wmo_id_comp, sensor_name_attribute_comp,dark_name_comp, mask_name_comp)


!--- loop through each parameter 

   istatus = 0
   sds_index_loop: do sds_index = 0,num_datasets-1


!-- loop over files
      file_loop_2: do ifile = 1, n_files

         if (file_mask(ifile) == 0) then
            cycle
         endif

         sds_id = sfselect(sd_id(ifile),sds_index)

         sds_name = " "
         istatus = sfginfo(sds_id, sds_name, sds_rank, dimsizes, sds_data_type, num_attrs) + istatus

         if (istatus /= 0) then
            print *, "could not find sds ", sds_index
            exit
         endif

!--- output progress to standard output
         if (sds_index == 1) then
          print *, "compositing ", trim(file_in(ifile)), " for node = ", node_string
         endif

!--- for first file, set comp to missing
!        if (sd_id(ifile) == sd_id_first) then
!         if (sds_data_type == DFNT_INT8) z_comp = missing_value_int1
!         if (sds_data_type == DFNT_INT16) z_comp = missing_value_int2
!         if (sds_data_type == DFNT_INT32) z_comp = missing_value_int4
!         if (sds_data_type == DFNT_FLOAT32) z_comp = missing_value_real4
!        endif
  

!--- allocate space for this orbit's sds
         allocate(z_orbit(dimsizes(1)))

!--- read sds for this orbit
         call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile), &
                              trim(sds_name),z_orbit, &
                              sds_data_type,scaled,sds_units, &
                              unscaled_min,unscaled_max,unscaled_missing,&
                              scaled_min,scaled_max,scaled_missing,&
                              istatus)

         if (sd_id(ifile) == sd_id_first) then
            z_comp = unscaled_missing
         endif

!--- populate composite
             do i = 1, num_cells_with_data_orbit(ifile)
              if (file_index_comp(int(cell_index_orbit(i,ifile))) == ifile) then
               z_comp(int(cell_index_orbit(i,ifile))) = z_orbit(i)
              endif
             enddo

!--- deallocate space for this orbit's sds
           deallocate(z_orbit)
         
      end do file_loop_2

!--- write out composited sds to output file
         call WRITE_CLAVRX_HDF4_SDS(sd_id_comp,real(z_comp,kind=real4), &
                  trim(sds_name),&
                  sds_data_type,scaled,unscaled_min,unscaled_max,trim(sds_units), unscaled_missing, &
                  "grid_cell_number", compress_flag, istatus_write)

   end do sds_index_loop 

!--- close file
   istatus = sfend(sd_id_comp)

end do asc_des_loop


!--- deallocate memory
  deallocate(file_index_comp,satzen_comp,num_cells_with_data_orbit,cell_index_orbit,z_comp)

end program COMPILE_ASC_DES 
