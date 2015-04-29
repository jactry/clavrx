! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: comp_time.f90 (src)
!       COMPILE_TIME (program)
!
! PURPOSE: this code takes the orbital gridcell files, combines them and writes
!          separate gridcell files for the ascending and descending nodes
!
! DESCRIPTION: the format of the files are the same.  The number of gridcells in the
!              orbital files is variable but the number of gridcells in the asc/des
!              files is the maximum number (full global coverage)
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
!   Created May 2004
!   August 2004 - Added error messaging (D. Donahue - STC)
!   June 2006 - Rewrote for version 4 CLAVR-x using standardized HDF read/write commands
!----------------------------------------------------------------------
program COMPILE_TIME
 use CONSTANTS
 use HDF
 use HDF_PARAMS
 use NUMERICAL_ROUTINES
 use SCALING_PARAMETERS
 implicit none

!------
 integer, parameter:: n_files_max = 100
 character(len=255), dimension(n_files_max):: file_in
 integer, dimension(n_files_max):: sd_id

!----
 character(len=355)::  dir_in, dir_out, file_out, file_temp
 integer:: ios
 integer(kind=int2):: jday, year, time
 real(kind=real4):: time_window, time_temp,time_diff, time_diff_comp

 character(len=4):: year_string
 character(len=3):: jday_string,grid_resolution_string
 character(len=2):: time_string
 character(len=1):: grid_format_string
 real(kind=real4):: grid_resolution_input
 integer(kind=int4):: grid_format_flag_input, compress_flag

!--- global attributes for composite
 integer(kind=int4):: sd_id_comp
 integer(kind=int4):: num_cells_comp,num_cells_with_data_comp
 integer:: therm_cal_1b_comp,Ref_cal_1b_comp,nav_opt_comp,use_sst_anal_comp, &
                         modis_clr_alb_flag_comp, nwp_flag_comp
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
                                                                                                                                                         
 character(len=255):: file_name_comp, file_1b_comp, creator_comp ,  &
                    plang_comp, hdf_ver_comp, hdf_timestamp_comp, data_type_comp,grid_format_comp
 character(len=20):: sensor_name_attribute_comp
 character(len=120):: dark_name_comp
 character(len=120):: mask_name_comp

!--- global attributes for orbital files
 integer(kind=int4):: num_cells_orbit
 integer:: therm_cal_1b_orbit,Ref_cal_1b_orbit,nav_opt_orbit,use_sst_anal_orbit, &
                         modis_clr_alb_flag_orbit, nwp_flag_orbit
 integer(kind=int4)::  start_time_orbit,end_time_orbit
 integer(kind=int4)::  acha_mode_orbit,dcomp_mode_orbit,wmo_sc_code_orbit
 integer(kind=int2)::  start_year_orbit,end_year_orbit,start_day_orbit,end_day_orbit
 real(kind=real4)::  resolution_km
 real(kind=real4)::  grid_resolution_orbit, ch1_gain_low_orbit, ch1_gain_high_orbit, &
       ch2_gain_low_orbit, ch2_gain_high_orbit, &
       ch3a_gain_low_orbit, ch3a_gain_high_orbit, &
       ch1_switch_count_orbit, ch1_dark_count_orbit, &
       ch2_switch_count_orbit, ch2_dark_count_orbit, &
       ch3a_switch_count_orbit, ch3a_dark_count_orbit, &
       sun_earth_distance_orbit, &
       c1_orbit, c2_orbit, a_3b_orbit, b_3b_orbit, nu_3b_orbit, &
       a_31_orbit, b_31_orbit, nu_31_orbit, a_32_orbit, b_32_orbit, nu_32_orbit, &
       solar_3b_nu_orbit,timerr_seconds_orbit
                                                                                                                                                         
 character(len=255):: file_name_orbit, file_1b_orbit, creator_orbit,  &
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
  real, dimension(:), allocatable:: time_comp,z_comp,z_orbit,time_orbit,jday_orbit,year_orbit
  real, dimension(:,:), allocatable:: cell_index_orbit

!--- hdf
   integer:: sfstart, sfend, sfginfo,sfselect, sffinfo

!

!-----------------------------------------------------
! open input file
!----------------------------------------------------
open(unit=8,file="comp_time_input",status="old",action="read",position="rewind",iostat=ios)
if (ios /= 0) then
  print *, "error opening comp_asc_dec control file, ios = ", ios
  stop 1
else
  read(unit=8,fmt="(A)") dir_in
  read(unit=8,fmt="(A)") dir_out
  read(unit=8,fmt=*) year
  read(unit=8,fmt=*) jday
  read(unit=8,fmt=*) time
  read(unit=8,fmt=*) time_window
  read(unit=8,fmt=*) grid_resolution_input
  read(unit=8,fmt=*) grid_format_flag_input
endif

!--- set compression flag (0 = none, 1 = gzip, 2 = szip)
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
call READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES(sd_id(1),data_type_comp,file_name_comp,file_1b_comp, &
                           resolution_km, &
                           start_year_comp,end_year_comp,start_day_comp,end_day_comp,start_time_comp,end_time_comp,&
                           num_cells_comp,num_cells_with_data_comp,grid_format_comp,grid_resolution_comp, &
                           therm_cal_1b_comp,Ref_cal_1b_comp,nav_opt_comp,use_sst_anal_comp, &
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
                           acha_mode_comp, dcomp_mode_comp, &
                           wmo_sc_code_comp, sensor_name_attribute_comp, &
                           dark_name_comp, mask_name_comp, &
                           creator_comp, plang_comp, hdf_ver_comp,hdf_timestamp_comp )

!--- determine number of sds's in these files (assume the same)
     istatus = sffinfo(sd_id(1),num_datasets,num_global_attrs)

!--- allocate memory needed for composities
     allocate(file_index_comp(num_cells_comp),time_comp(num_cells_comp))
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
    write (time_string,   '(I2.2)') time
    write (grid_format_string,   '(I1.1)') grid_format_flag

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

!--- initialize composite variables
     file_index_comp = 0
     cell_index_orbit = 0
     num_cells_with_data_orbit = 0
     time_comp = missing_value_real4 

!-----------------------------------------------------------------------------------------
!--- loop over files, reading time (the composite basis) and generate file index
!-----------------------------------------------------------------------------------------
     do ifile = 1,n_files

!--- read global attributes of this file
             call READ_CLAVRX_HDF_GLOBAL_ATTRIBUTES(sd_id(ifile),data_type_orbit,file_name_orbit,file_1b_orbit, &
                           resolution_km, &
                           start_year_orbit,end_year_orbit,start_day_orbit,end_day_orbit,start_time_orbit,end_time_orbit,&
                           num_cells_orbit,num_cells_with_data_orbit(ifile),grid_format_orbit,grid_resolution_orbit, &
                           therm_cal_1b_orbit,Ref_cal_1b_orbit,nav_opt_orbit,use_sst_anal_orbit, &
                           modis_clr_alb_flag_orbit,  &
                           nwp_flag_orbit, ch1_gain_low_orbit, ch1_gain_high_orbit, &
                           ch1_switch_count_orbit, ch1_dark_count_orbit, &
                           ch2_gain_low_orbit, ch2_gain_high_orbit, &
                           ch2_switch_count_orbit, ch2_dark_count_orbit, &
                           ch3a_gain_low_orbit, ch3a_gain_high_orbit, &
                           ch3a_switch_count_orbit, ch3a_dark_count_orbit, &
                           sun_earth_distance_orbit, &
                           c1_orbit, c2_orbit, a_3b_orbit, b_3b_orbit, nu_3b_orbit, &
                           a_31_orbit, b_31_orbit, nu_31_orbit, a_32_orbit, b_32_orbit, nu_32_orbit, &
                           solar_3b_nu_orbit, timerr_seconds_orbit, &
                           acha_mode_orbit, dcomp_mode_orbit, &
                           wmo_sc_code_orbit, sensor_name_attribute_orbit, &
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
       endif

!--- check if file is of the right time
       if ((grid_resolution_comp /= grid_resolution_orbit) .or. & 
           (grid_format_comp /= grid_format_orbit)) then
            file_mask(ifile) = 0
       endif

!--- if not in window, move to next orbit
       if (file_mask(ifile) == 0) then    !MAKE A FILE MASK
           cycle
       endif    

!--- allocate memory needed
          allocate(jday_orbit(num_cells_with_data_orbit(ifile)))   !this is not saved
          allocate(year_orbit(num_cells_with_data_orbit(ifile)))   !this is not saved
          allocate(time_orbit(num_cells_with_data_orbit(ifile)))   !this is not saved
          jday_orbit = 0.0
          year_orbit = 0.0
          time_orbit = 0.0

!--- read cell_index
          call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile),"cell_index", &
                                  cell_index_orbit(1:num_cells_with_data_orbit(ifile),ifile), &
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

!--- read day
          call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile),"jday", &
                                  jday_orbit, &
                                  sds_data_type,scaled,sds_units, &
                                  unscaled_min,unscaled_max,unscaled_missing,&
                                  scaled_min,scaled_max,scaled_missing,&
                                  istatus)
!--- read year
          call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile),"year", &
                                  year_orbit, &
                                  sds_data_type,scaled,sds_units, &
                                  unscaled_min,unscaled_max,unscaled_missing,&
                                  scaled_min,scaled_max,scaled_missing,&
                                  istatus)


!--- make file index
    do i = 1, num_cells_with_data_orbit(ifile)

      iupdate = 0
     
!---- compute time to account for jday_in being yesterday or tomorrow
      time_temp = time_orbit(i) + 24.0*(jday_orbit(i) - jday) +   &
                                365.0*24.0*(year_orbit(i) - year)

      time_diff = abs(time_temp - time)

!--- if composite missing, see if cell falls within window
      if (time_comp(int(cell_index_orbit(i,ifile))) == missing_value_real4) then

!--- if inside window and not missing, update
       if ((time_orbit(i) /= missing_value_real4) .and. &
          (time_diff < time_window)) then
          iupdate = 1
       endif

      else

!--- if composite not missing, see if closer than current composite value
        time_diff_comp = abs(time_comp(int(cell_index_orbit(i,ifile))) - time)
        if (time_diff < time_diff_comp) then
          iupdate = 1
         endif

      endif

!--- choose this cell for composite if it met the tests
              if (iupdate == 1) then
                 file_index_comp(int(cell_index_orbit(i,ifile))) = ifile
                 time_comp(int(cell_index_orbit(i,ifile))) = time_temp   !note, not time_orbit
              endif

      enddo

!--- deallocate memory
     deallocate(time_orbit, jday_orbit, year_orbit)


     enddo


!--- construct output file name for this node
    file_out = "clavrx_"//trim(time_string)//"z_"// &
               trim(grid_resolution_string)//"_"//trim(grid_format_string)//"_"// &
               trim(year_string)//"_"//trim(jday_string)// &
               ".level3.hdf"

!--- open file for writing
    sd_id_comp = sfstart(trim(dir_out)//trim(file_out),DFACC_CREATE)

     if (sd_id_comp <  0) then
       print *, "composite HDF file creation failed. Exiting..."
       stop 38
     endif

!--- write global attributes to output file
 call WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES(sd_id_comp,"GRIDCELL",trim(file_out),trim(file_1b_comp), &
         resolution_km, &
         start_year_comp,end_year_comp,start_day_comp,end_day_comp,start_time_comp,end_time_comp,&
         num_cells_comp,num_cells_with_data_comp,grid_format_comp,grid_resolution_comp, &
         therm_cal_1b_comp,Ref_cal_1b_comp,nav_opt_comp,use_sst_anal_comp, &
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
         acha_mode_comp, dcomp_mode_comp,&
         wmo_sc_code_comp, sensor_name_attribute_comp,dark_name_comp,mask_name_comp)


!--- loop through each parameter 

   istatus = 0
   sds_index_loop: do sds_index = 0,num_datasets-1

!-- reset composite for this sds
      z_comp = missing_value_real4

      file_loop_2: do ifile = 1, n_files

         if (file_mask(ifile) == 0) then
            cycle
         endif

         sds_id = sfselect(sd_id(ifile),sds_index)

         istatus = sfginfo(sds_id, sds_name, sds_rank, dimsizes, sds_data_type, num_attrs) + istatus

         if (istatus /= 0) then
            print *, "could not find sds ", sds_index
            exit
         endif

!--- output progress to standard output
         if (sds_index == 0) then
          print *, "compositing first sds for ", trim(file_in(ifile))
         endif

!--- allocate space for this orbit's sds

         allocate(z_orbit(dimsizes(1)))

!--- read sds for this orbit
         call READ_CLAVRX_HDF4_SDS_RANK1(sd_id(ifile),num_cells_with_data_orbit(ifile), &
                              trim(sds_name),z_orbit, &
                              sds_data_type,scaled,sds_units, &
                              unscaled_min,unscaled_max,unscaled_missing,&
                              scaled_min,scaled_max,scaled_missing,&
                              istatus)

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

!--- deallocate memory
  deallocate(file_index_comp,time_comp,num_cells_with_data_orbit,cell_index_orbit,z_comp)


end program COMPILE_TIME
