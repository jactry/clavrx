! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: user_options.f90 (src)
!       USER_OPTIONS (program)
!
! PURPOSE: CLAVR-x Module to house routines dealing with user options from the
!          default_options file or the command line
!
! DESCRIPTION: 
!             Routines in this module and their purpose:
!
!             READ_CLAVRXORB_DEFAULT_OPTIONS - Opens default file if necessary
!             QC_CLAVRXORB_OPTIONS - quality control options
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
!--------------------------------------------------------------------------------------
module USER_OPTIONS
 use PIXEL_COMMON
 use CONSTANTS
 use FILE_UTILITY

  implicit none

  public:: SETUP_USER_DEFINED_OPTIONS
  public:: CHECK_ALGORITHM_CHOICES

  private:: READ_CLAVRXORB_DEFAULT_OPTIONS, &
            HELPER, &
            QC_CLAVRXORB_OPTIONS

  character(24), parameter, private :: MOD_PROMPT = " USER_OPTIONS_ROUTINES: "
  integer :: Dcomp_Mode_User_Set
  contains
!------------------------------------------------------------------
!
!------------------------------------------------------------------
subroutine SETUP_USER_DEFINED_OPTIONS()

    call READ_CLAVRXORB_COMMANDLINE_OPTIONS()

    call QC_CLAVRXORB_OPTIONS()


end subroutine SETUP_USER_DEFINED_OPTIONS
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  subroutine CHECK_ALGORITHM_CHOICES()
     
   use  clavrx_message_module, only: mesg
 
 
 
     !------------------------------------------------------------------------
     !--- ACHA MODE Check
     !         (0=11; 1 = 11/12; 2 = 11/13.3; 
     !---       3 = 11,12,13; 4=8.5/11/12; 
     !---       5=11/12/6.7; 6=11/13.3/6.7; 7=11/6.7 )
     !------------------------------------------------------------------------

     if (Avhrr_Flag == sym%YES) then 
         if (Avhrr_1_Flag == sym%NO) then 
             if (Acha_Mode > 1) then
                     print *, EXE_PROMPT, &
                     "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT, &
                     "Changing to default for AVHRR"
                     Acha_Mode = Acha_Mode_Default_Avhrr 
             endif
         ELSE 
             if (Acha_Mode /= 0)  then 
                     print *, EXE_PROMPT, &
                     "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT, &
                     "Changing to default for AVHRR/1"
                     Acha_Mode = Acha_Mode_Default_Avhrr1 
             endif
         endif
     endif

print *, "Acha Mode = ", Acha_Mode
     if (Goes_Flag == sym%YES) then 

          if (Goes_Mop_Flag == sym%NO) then 

             if (Acha_Mode == 2 .or. Acha_Mode == 3 .or. &
                 Acha_Mode == 4 .or. Acha_Mode == 6)  then 
                     print *, EXE_PROMPT, &
                     "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT, &
                     "Changing to default for GOES-IL"
                     Acha_Mode = Acha_Mode_Default_Goes_IL
             endif

          ELSE 

             if (Acha_Mode == 1 .or. Acha_Mode == 3 .or.  &
                 Acha_Mode == 4 .or. Acha_Mode == 5)  then 
                     print *, EXE_PROMPT, &
                     "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT, &
                     "Changing to default for GOES-MP"
                     Acha_Mode = Acha_Mode_Default_Goes_MP
print *, "Acha Mode = ", Acha_Mode
             endif

          endif
     endif

     if (Viirs_Flag == sym%YES) then 

             if (Acha_Mode == 2 .or. Acha_Mode == 3 .or.  &
                 Acha_Mode == 5 .or. Acha_Mode == 6 .or. Acha_Mode == 7)  then 
                     print *, EXE_PROMPT, "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT,  "Changing to default for VIIRS"
                     Acha_Mode = Acha_Mode_Default_VIIRS
             endif
 
     endif

     if (Mtsat_Flag == sym%YES) then

             if (Acha_Mode == 2 .or. Acha_Mode == 3 .or. &
                 Acha_Mode == 4 .or. Acha_Mode == 6)  then
                     print *, EXE_PROMPT, "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT,  "Changing to default for MTSAT"
                     Acha_Mode = Acha_Mode_Default_MTSAT
             endif

     endif

     !--- check ACHA mode based on available channels
     if (Acha_Mode == 1 .and. &
         (Chan_On_Flag_Default(32)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 1 not possible with selected channels, ACHA Set to Mode 0'
         Acha_Mode = 0
     endif
     if (Acha_Mode == 2 .and. &
         (Chan_On_Flag_Default(33)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 2 not possible with selected channels, ACHA Set to Mode 0'
         Acha_Mode = 0
     endif
     if (Acha_Mode == 3 .and. &
         (Chan_On_Flag_Default(32)==sym%NO .or. Chan_On_Flag_Default(33)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 3 not possible with selected channels, ACHA Set to Mode 0'
         Acha_Mode = 0
     endif
     if (Acha_Mode == 4 .and. &
         (Chan_On_Flag_Default(29)==sym%NO .or. Chan_On_Flag_Default(32)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 4 not possible with selected channels, ACHA Set to Mode 0'
         Acha_Mode = 0
     endif
     if (Acha_Mode == 5 .and. &
         (Chan_On_Flag_Default(27)==sym%NO .or. Chan_On_Flag_Default(32)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 5 not possible with selected channels, ACHA Set to Mode 0'
         Acha_Mode = 0
     endif
     if (Acha_Mode == 6 .and. &
         (Chan_On_Flag_Default(27)==sym%NO .or. Chan_On_Flag_Default(33)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 6 not possible with selected channels, ACHA Set to Mode 0'
         Acha_Mode = 0
     endif
     if (Acha_Mode == 7 .and. &
         (Chan_On_Flag_Default(27)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 7 not possible with selected channels, ACHA Set to Mode 0'
         Acha_Mode = 0
     endif


     !-------------------------------------------------------------------------------------
     !--- DCOMP Mode Check
     !-------------------------------------------------------------------------------------

     !- dcomp mode 9 is Andys test code
     print*,'dcomp user set:   ===================> ', dcomp_mode_user_set
     dcomp_mode = dcomp_mode_user_set
   
     if (AVHRR_Flag == sym%YES .and. Dcomp_Mode /= 0  .and.Dcomp_Mode /= 9) then
        Dcomp_Mode = 3
        if (Sc_Id_WMO == 208) Dcomp_Mode = 1   !NOAA-17
        if (Sc_Id_WMO == 3) Dcomp_Mode = 1     !METOP-A
        if (Sc_Id_WMO == 4) Dcomp_Mode = 1     !METOP-B
        if (Sc_Id_WMO == 5) Dcomp_Mode = 1     !METOP-C
     endif
     if (GOES_Flag == sym%YES .and. Dcomp_Mode /= 0 .and.Dcomp_Mode /= 9) then
        Dcomp_Mode = 3
     endif
     if (Mtsat_Flag == sym%YES .and. Dcomp_Mode /= 0  .and.Dcomp_Mode /= 9) then
        Dcomp_Mode = 3
     endif
     if  ( Dcomp_Mode_user_set .ne. Dcomp_Mode) then
	   call mesg ( 'dcomp mode switched due to sensor  setting  ',color=91, level =-1)
     endif 
     print*,'dcomp mode used for this file: ' , dcomp_mode


     !--- check based on available channels
     if (Dcomp_Mode_User_Set == 1 .and. &
         (Chan_On_Flag_Default(1) == sym%NO .or. Chan_On_Flag_Default(6)==sym%NO)) then
         print *, EXE_PROMPT, 'DCOMP Mode 1 not possible with selected channels, DCOMP is now off'
     endif
     if (Dcomp_Mode_User_Set == 2 .and. &
         (Chan_On_Flag_Default(1) == sym%NO .or. Chan_On_Flag_Default(7)==sym%NO)) then
         print *, EXE_PROMPT, 'DCOMP Mode 2 not possible with selected channels, DCOMP is now off'
     endif
     if (Dcomp_Mode_User_Set == 3 .and. &
         (Chan_On_Flag_Default(1) == sym%NO .or. Chan_On_Flag_Default(20)==sym%NO)) then
         print *, EXE_PROMPT, 'DCOMP Mode 3 not possible with selected channels, DCOMP is now off'
     endif
    
     !--- Prevent level3 file creation for anything but AVHRR
     if (Level3_Flag == sym%YES .and. Avhrr_Flag == sym%NO) then
         Level3_Flag = sym%NO
         print *, EXE_PROMPT, "Level3 creation supported only for AVHRR, turning off"
     endif

  end subroutine CHECK_ALGORITHM_CHOICES
!------------------------------------------------------------------
! Read Parameters from AVHRR INPUT files and check them for errors
!
! parameters are pased in avhrr_pixel_common public memory
!------------------------------------------------------------------
  subroutine READ_CLAVRXORB_DEFAULT_OPTIONS(File_Default)
    character(len=*), intent(in):: File_Default
    integer::ios0
    integer::erstat
    integer:: Default_Lun
    
    print *, EXE_PROMPT, "DEFAULT FILE READ IN"
    print *, EXE_PROMPT, "Default file to be read in: ", trim(File_Default)

    Default_Lun = GET_LUN()

    open(unit=Default_Lun,file=trim(File_Default), &
          iostat = ios0, &
          action="read", &
          status="old", &
          position="rewind")

    erstat = 0
    if (ios0 /= 0) then    !avhrr_input_check
      erstat = 1
      print *, EXE_PROMPT, "error opening AVHRR default control file, ios0 = ", ios0
      stop 1
   else
      read(unit=Default_Lun,fmt=*) Ref_Cal_1b
      read(unit=Default_Lun,fmt=*) Therm_Cal_1b
      read(unit=Default_Lun,fmt=*) Bx_File_Flag
      read(unit=Default_Lun,fmt=*) Nav_Flag
      read(unit=Default_Lun,fmt=*) Nav_File_Flag
      read(unit=Default_Lun,fmt=*) Cmr_File_Flag
      read(unit=Default_Lun,fmt=*) Obs_File_Flag
      read(unit=Default_Lun,fmt=*) Geo_File_Flag
      read(unit=Default_Lun,fmt=*) Cld_File_Flag
      read(unit=Default_Lun,fmt=*) Sst_File_Flag
      read(unit=Default_Lun,fmt=*) Rtm_File_Flag
      read(unit=Default_Lun,fmt=*) Ash_File_Flag
      read(unit=Default_Lun,fmt=*) Level2_File_Flag
      read(unit=Default_Lun,fmt=*) Level3_Flag
      read(unit=Default_Lun,fmt=*) Cloud_Mask_Aux_Flag
      read(unit=Default_Lun,fmt=*) Cloud_Mask_Bayesian_Flag
      read(unit=Default_Lun,fmt=*) Blank_Flag
      read(unit=Default_Lun,fmt=*) Cld_Flag
      read(unit=Default_Lun,fmt=*) Aer_Flag
      read(unit=Default_Lun,fmt=*) Erb_Flag
      read(unit=Default_Lun,fmt=*) Ash_Flag
      read(unit=Default_Lun,fmt=*) Use_Sst_Anal_Default
      read(unit=Default_Lun,fmt=*) Sst_Anal_Opt
      read(unit=Default_Lun,fmt=*) Dlat
      read(unit=Default_Lun,fmt=*) Level3_Format
      read(unit=Default_Lun,fmt=*) Data_Comp_Flag
      read(unit=Default_Lun,fmt=*) Subset_Pixel_Hdf_Flag
      read(unit=Default_Lun,fmt=*) Nwp_Flag
      read(unit=Default_Lun,fmt=*) Rtm_Flag
      read(unit=Default_Lun,fmt=*) Modis_Clr_Alb_Flag
      read(unit=Default_Lun,fmt=*) Prob_Clear_Res_Flag
      read(unit=Default_Lun,fmt=*) Lrc_Flag
      read(unit=Default_Lun,fmt=*) Process_Undetected_Cloud_Flag
      read(unit=Default_Lun,fmt=*) Diag_Flag
      read(unit=Default_Lun,fmt=*) ASc_Flag_Diag
      read(unit=Default_Lun,fmt=*) Lat_Min_Limit
      read(unit=Default_Lun,fmt=*) Lat_Max_Limit
      read(unit=Default_Lun,fmt="(a)") Ancil_Data_Dir
      read(unit=Default_Lun,fmt="(a)") Gfs_Data_Dir
      read(unit=Default_Lun,fmt="(a)") Ncep_Data_Dir
      read(unit=Default_Lun,fmt="(a)") Cfsr_Data_Dir
      read(unit=Default_Lun,fmt="(a)") Oisst_Data_Dir
      read(unit=Default_Lun,fmt="(a)") Snow_Data_Dir
      read(unit=Default_Lun,fmt="(a)") GlobSnow_Data_Dir
      read(unit=Default_Lun,fmt="(a)") Dark_Comp_Data_Dir
      read(unit=Default_Lun,fmt="(a)") Temporary_Data_Dir
      read(unit=Default_Lun,fmt=*) Smooth_Nwp_Flag
      read(unit=Default_Lun,fmt=*) Use_Seebor
      read(unit=Default_Lun,fmt=*) Read_Hires_Sfc_Type
      read(unit=Default_Lun,fmt=*) Read_Land_Mask
      read(unit=Default_Lun,fmt=*) Read_Coast_Mask
      read(unit=Default_Lun,fmt=*) Read_Surface_Elevation
      read(unit=Default_Lun,fmt=*) Read_Volcano_Mask
      read(unit=Default_Lun,fmt=*) Solzen_min_limit, Solzen_max_limit
      read(unit=Default_Lun,fmt=*) Read_Snow_Mask
      read(unit=Default_Lun,fmt=*) Read_Dark_Comp
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(1:6)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(7:12)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(13:18)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(19:24)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(25:30)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(31:36)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(37:42)
      read(unit=Default_Lun,fmt=*) Dcomp_Mode_user_set
      read(unit=Default_Lun,fmt=*) Acha_Mode
      read(unit=Default_Lun,fmt="(a)") bayesian_cloud_mask_name
    endif
    close(unit=Default_Lun)

 end subroutine READ_CLAVRXORB_DEFAULT_OPTIONS
 
 subroutine READ_CLAVRXORB_COMMANDLINE_OPTIONS()
!------------------------------------------------------------------
! Command-line input option variables 
!------------------------------------------------------------------  
  
  character(len=30) :: fargv
  character(len=30) :: junk
  character(len=1) :: temp_string
  logical:: back
  character(len=300):: default_temp
  real :: int_temp
  integer :: fargc
  integer :: i
  integer :: Temp_Scans_Arg !--- temporary integer for number of scanlines
  integer*4 :: iargc

  temp_string = '.'  !--- tempoary string to search for in angle commandline

  !---- SET DEFAULT OPTIONS
  
  use_Default = sym%YES
  default_temp="./clavrxorb_default_options"
  file_list = "./clavrxorb_file_list"
  
  !--- set yes/no options to no as default
  Ref_Cal_1b = sym%NO
  therm_Cal_1b = sym%NO
  bx_file_Flag = sym%NO
  nav_Flag = 0
  prob_clear_res_Flag = 0
  lrc_Flag = sym%NO
  process_undetected_cloud_Flag = sym%NO
  nav_file_Flag = sym%NO 
  cmr_file_Flag = sym%NO
  cld_file_Flag = sym%NO
  sst_file_Flag = sym%NO
  obs_file_Flag = sym%NO
  geo_file_Flag = sym%NO
  rtm_file_Flag = sym%NO
  ash_file_Flag = sym%NO
  level2_file_Flag = sym%NO
  level3_Flag = sym%NO
  cld_Flag = sym%NO
  Aer_Flag = sym%NO
  Erb_Flag = sym%NO
  Ash_Flag = sym%NO
  Data_comp_Flag = 0
  Diag_Flag = 0
  Subset_pixel_hdf_Flag = 0
  use_seebor = 0
  Smooth_Nwp_Flag = 0
  read_hires_sfc_type = 0
  read_land_mask = 0
  read_coast_mask = 0
  read_surface_elevation = 0
  read_volcano_mask = 0
  Cloud_Mask_Aux_Flag = 0
  Cloud_Mask_Bayesian_Flag = 0
  use_sst_anal_Default = 0 ! Do not use OISST analisys
  sst_anal_opt = 0 !determine OISST file
  modis_clr_alb_Flag = 0 ! do not use clear-sky MODIS albedo maps
  read_snow_mask = 0
  output_scaled_reflectances = sym%NO !default is to output ref / cosSolzen
  Num_Scans_per_Segment = 240
  Temp_Scans_Arg = 0
 
  
  !--- Default is equal angle @ 0.5 degree resolution
  level3_format = 1
  dlat = 0.5
  
  !--- default solar zenith limits
  Solzen_Min_Limit= 0 
  Solzen_Max_Limit= 180.0
  
  !--- default latitude limits
  Lat_Min_Limit = -90.0
  Lat_Max_Limit =  90.0

  !--- default NWP and RTM options
  Nwp_Flag = 1     !--- use GFS
  Rtm_Flag = 1     !--- use PFAST
  
  !--- default ancillary data directory
  ancil_data_dir = './data/'
  gfs_data_dir = './data/gfs_archive/'
  ncep_data_dir = './data/ncep-reanalysis/'
  cfsr_data_dir = './data/cfsr/'
  oisst_data_dir = './data/oisst_archive/'
  snow_data_dir = './data/snow/'
 
  fargc = iargc()
  
  !--- first we will check to see if the default file is used
  !--- also check to see if the help file is to be displayed
  
  do i=1, fargc
        call getarg(i,fargv)
        if (trim(fargv) == "-help" .or. &
             trim(fargv) == "-h") then
          call HELPER()
          stop
        !No default file used
        
        else if  ( trim(fargv) == "-version" .or. &
             trim(fargv) == "-ver") then
          print*,'clavrx_5.4 trunk'
          print*,'$Header$'
          stop
        
        elseif (trim(fargv) == "-no_Default") then 
          use_Default = sym%NO
        !Different default file used
        elseif (trim(fargv) == "-default") then
          call getarg(i+1,default_temp)
          default_temp=trim(default_temp)
       endif
  enddo
  
  
  !---- If the default file is used, read it in first, then
  !---- check for other command line changes to the options
  

  if(Use_Default == sym%YES)  then
    ! call READ_CLAVRXORB_DEFAULT_OPTIONS(Default_Temp)
  endif
  if(Use_Default == sym%NO)  then
          print *, EXE_PROMPT, "Using standard defaults and command line options"
  endif 
 
  do i=1, fargc
        call getarg(i,fargv)
    
        !Change Ref_Cal_1b flag 
        if (trim(fargv) == "-Ref_Cal_1b") then 
          Ref_Cal_1b = sym%YES
        elseif (trim(fargv) == "-no_Ref_Cal_1b") then 
          Ref_Cal_1b = sym%NO
        
        !Change therm_Cal_1b flag
        elseif(trim(fargv) == "-therm_Cal_1b") then
          therm_Cal_1b = sym%YES
        elseif(trim(fargv) == "-no_Therm_Cal_1b") then
          therm_Cal_1b = sym%NO
        
        !Change Nav type
        elseif(trim(fargv) == "-l1bnav") then
          nav_Flag = 0
          
        !Change Nav output flag
        elseif(trim(fargv) == "-nav_file") then
          nav_file_Flag = sym%YES
        elseif(trim(fargv) == "-no_nav_file") then
          nav_file_Flag = sym%NO

        !Change RTM output flag
        elseif(trim(fargv) == "-rtm_file") then
          rtm_file_Flag = sym%YES
        elseif(trim(fargv) == "-no_rtm_file") then
          rtm_file_Flag = sym%NO

        !Change ash output flag
        elseif(trim(fargv) == "-ash_file") then
          ash_file_Flag = sym%YES
        elseif(trim(fargv) == "-no_ash_file") then
          ash_file_Flag = sym%NO

        !Change level2 output flag
        elseif(trim(fargv) == "-level2_file") then
          level2_file_Flag = sym%YES
        elseif(trim(fargv) == "-no_level2_file") then
          level2_file_Flag = sym%NO

        !Change cloud mask
        elseif(trim(fargv) == "-cloud_mask_aux") then
          Cloud_Mask_Aux_Flag = sym%YES
        elseif(trim(fargv) == "-cloud_mask_calc") then
          Cloud_Mask_Aux_Flag = sym%NO

        !Change data compression flag
        elseif(trim(fargv) == "-no_data_comp") then
          data_comp_Flag=0
        elseif(trim(fargv) == "-data_comp_gzip") then
          data_comp_Flag=1
        elseif(trim(fargv) == "-data_comp_szip") then
          data_comp_Flag=2

        !Subset pixel HDF
        elseif(trim(fargv) == "-subset_pixel_hdf") then
          subset_pixel_hdf_Flag = sym%YES
        elseif(trim(fargv) == "-no_subset_pixel_hdf") then
          subset_pixel_hdf_Flag = sym%NO

        !Change lat max/min for processing
        elseif(trim(fargv) == "-Lat_min_limit") then
           call getarg(i+1,junk)
           int_temp = scan(junk,temp_string, back) 
           if(int_temp > 0.0) read(junk,'(f6.3)') Lat_Min_Limit
           if(int_temp == 0.0) read(junk,'(f6.0)') Lat_Max_Limit
        elseif(trim(fargv) == "-lat_max_limit") then
          call getarg(i+1,junk)
          int_temp = scan(junk,temp_string, back) 
          if(int_temp > 0.0) read(junk,'(f6.3)') Lat_Max_Limit
          if(int_temp == 0.0) read(junk,'(f6.0)') Lat_Max_Limit

        !Change ancillary data directory
        elseif(trim(fargv) == "-ancil_data_dir") then
          call getarg(i+1,ancil_data_dir)
          ancil_data_dir=trim(ancil_data_dir)

        !Change GFS data directory
        elseif(trim(fargv) == "-gfs_data_dir") then
          call getarg(i+1,gfs_data_dir)
          gfs_data_dir=trim(gfs_data_dir)

        !Change NCEP reanalysis data directory
        elseif(trim(fargv) == "-ncep_data_dir") then
          call getarg(i+1,ncep_data_dir)
          ncep_data_dir=trim(ncep_data_dir)

        !Change CFSR reanalysis data directory
        elseif(trim(fargv) == "-cfsr_data_dir") then
          call getarg(i+1,cfsr_data_dir)
          cfsr_data_dir=trim(cfsr_data_dir)

        !Change OISST data directory
        elseif(trim(fargv) == "-oisst_data_dir") then
          call getarg(i+1,oisst_data_dir)
          oisst_data_dir=trim(oisst_data_dir)

        !Change SNOW data directory
        elseif(trim(fargv) == "-snow_data_dir") then
          call getarg(i+1,snow_data_dir)
          snow_data_dir=trim(snow_data_dir)

        !Smooth/not smooth NWP data
        elseif(trim(fargv) == "-smooth_nwp") then
          Smooth_Nwp_Flag = sym%YES
        elseif(trim(fargv) == "-no_smooth_nwp") then
          Smooth_Nwp_Flag = sym%NO

        !Use/not use SEEBOR emissivity
        elseif(trim(fargv) == "-use_seebor") then
          use_seebor = sym%YES
        elseif(trim(fargv) == "-no_seebor") then
          use_seebor = sym%NO

        !Change which surface type flag 
        elseif(trim(fargv) == "-high_res") then
          read_hires_sfc_type = 1
        elseif(trim(fargv) == "-low_res") then
          read_hires_sfc_type = 0

        !Read/not read land mask
        elseif(trim(fargv) == "-read_land_mask") then
          read_land_mask = sym%YES
        elseif(trim(fargv) == "-no_land_mask") then
          read_land_mask = sym%NO

        !Read/not read coast mask
        elseif(trim(fargv) == "-read_coast_mask") then
          read_coast_mask = sym%YES
        elseif(trim(fargv) == "-no_coast_mask") then
          read_coast_mask = sym%NO

        !Read/not read surface elevation
        elseif(trim(fargv) == "-no_surface_elevation") then
          read_surface_elevation = 0
        elseif(trim(fargv) == "-high_surface_elevation") then
          read_surface_elevation = 1
        elseif(trim(fargv) == "-low_surface_elevation") then
          read_surface_elevation = 2

        !Read/not read volcano mask
        elseif(trim(fargv) == "-read_volcano_mask") then
          read_volcano_mask = sym%YES
        elseif(trim(fargv) == "-no_volcano_mask") then
          read_volcano_mask = sym%NO

        !read snow mask
        elseif(trim(fargv) == "-snow") then
          read_snow_mask = sym%YES
        elseif(trim(fargv) == "-no_snow") then
          read_snow_mask = sym%NO

        !Output scaled reflectances or not
        elseif(trim(fargv) == "-output_scaled_ref") then
          output_scaled_reflectances = sym%YES
        elseif(trim(fargv) == "-no_output_scaled_ref") then
          output_scaled_reflectances = sym%NO

        ! change the number of scanlines per segment
        elseif(trim(fargv) == "-lines_per_seg") then
        call getarg(i+1,junk)
        read(junk,'(i4)') Temp_Scans_Arg
        if(Temp_Scans_Arg > 1) read(junk,'(i4)') Num_Scans_per_Segment

        !Change solar zenith angle limits
        elseif(trim(fargv) == "-solzen_min_limit") then
          call getarg(i+1,junk)
          int_temp = scan(junk,temp_string, back) 
          if(int_temp > 0.0) read(junk,'(f6.3)') Solzen_Min_Limit
          if(int_temp == 0.0) read(junk,'(f6.0)') Solzen_Min_Limit
        
       !Change dcomp mode
        elseif(trim(fargv) == "-dcomp_mode") then
          call getarg(i+1,junk)
          read(junk,'(i1)') Dcomp_Mode_user_set
          
       !Change dcomp mode
        elseif(trim(fargv) == "-Dcomp_Mode") then
          call getarg(i+1,junk)
          read(junk,'(i1)') Dcomp_Mode_user_set
          
    elseif(trim(fargv) == "-solzen_max_limit") then
       call getarg(i+1,junk)
       int_temp = scan(junk,temp_string, back) 
       if(int_temp > 0.0) read(junk,'(f6.3)') Solzen_Max_Limit
       if(int_temp == 0.0) read(junk,'(f6.0)') Solzen_Max_Limit
    elseif (trim(fargv) == "-filelist") then
       call getarg(i+1,file_list)
       file_list=trim(file_list)
    endif
  enddo

 end subroutine READ_CLAVRXORB_COMMANDLINE_OPTIONS

!-------------------------------------------------------------------------------
!--- QC options and modify as needed
!-------------------------------------------------------------------------------
 subroutine QC_CLAVRXORB_OPTIONS()

    integer:: erstat

    !---- Since the NWP controls everything, we first check if an NWP is being used
    !---- before anything else is checked.  If no nwp, we stop processing

    if ((Nwp_Flag < 0) .or. (Nwp_Flag > 4)) then
       print *,  EXE_PROMPT, "unrecognized value for Nwp_Flag: ", Nwp_Flag
       stop "6-Nwp_Flag"
    endif
    if (Nwp_Flag == 1) then
       print *,  EXE_PROMPT, "GFS data will be used"
    else if (Nwp_Flag == 2) then
       print *,  EXE_PROMPT, "NCEP Reanalysis data will be used"
    else if (Nwp_Flag == 3) then
       print *,  EXE_PROMPT, "NCEP Climate Forecast System Reanalysis data will be used"
    else if (Nwp_Flag == 4) then
       print *,  EXE_PROMPT, "GDAS Reanalysis data will be used"
    endif
    
    if (Nwp_Flag == 0) then
       print *,  EXE_PROMPT, "No choice made for NWP data, will not run algoritms or orbital level3 files"
       Cld_Flag = sym%NO
       Aer_Flag = sym%NO
       Erb_Flag = sym%NO
       Ash_Flag = sym%NO
       Cld_File_Flag = sym%NO
       Cmr_File_Flag = sym%NO
       Ash_File_Flag = sym%NO
       Sst_File_Flag = sym%NO
       Rtm_File_Flag = sym%NO
       Cloud_Mask_Bayesian_Flag = sym%NO
       Level3_Flag = sym%NO
       Cloud_Mask_Aux_Flag = sym%NO ! this is to determine if the lut's are being read in
    endif

    if (cloud_mask_bayesian_Flag == sym%YES) then
       print *, EXE_PROMPT, "Bayesian cloud mask will be generated"
    endif

    if (Ref_Cal_1b == sym%YES) then
       print *, EXE_PROMPT, "Reflectance Calibration within 1b will be used"
    endif

    if (therm_Cal_1b == sym%YES) then
       print *, EXE_PROMPT, "Thermal Calibration within 1b will be used"
    endif

    if (nav_Flag == 1) then
        print *,  EXE_PROMPT, "CLEVERNAV geolocation no longer supported, using REPOSNX"
        nav_Flag = 2
    endif

    if (nav_Flag == 2) then
          print *,  EXE_PROMPT, "REPOSNX geolocation adjustment done"
    endif

    if (nav_file_Flag == sym%YES) then
        print *,  EXE_PROMPT, "nav file will be created"
    endif

    if (cld_file_Flag == sym%YES) then
        print *,  EXE_PROMPT, "cld file will be created"
    endif

    if (cmr_file_Flag == sym%YES) then
        print *, EXE_PROMPT, "cmr file will be created"
    endif

    if (sst_file_Flag == sym%YES) then
        print *, EXE_PROMPT, "sst file will be created"
    endif

    if (ash_file_Flag == sym%YES) then
        print *, EXE_PROMPT, "vol ash file will be created"
    endif

    if (cld_Flag == sym%NO) then
        print *, EXE_PROMPT, "Cloud products will not be created"
    endif

    if (aer_Flag == sym%NO) then
        print *, EXE_PROMPT, "Aerosol products will not be created"
    endif

    if (erb_Flag == sym%NO) then
        print *, EXE_PROMPT, "Radiative Flux products will not be created"
    endif
      
    if (obs_file_Flag == sym%YES) then
        print *, EXE_PROMPT, "obs file will be created"
    endif

    if (geo_file_Flag == sym%YES) then
        print *, EXE_PROMPT, "geo file will be created"
    endif

    if (rtm_file_Flag == sym%YES) then
        print *, EXE_PROMPT, "rtm file will be created"
    endif

    if (level3_Flag == sym%YES) then
        print *, EXE_PROMPT, "level3 file will be created"
    endif

    if (Cloud_Mask_Aux_Flag == sym%YES) then
       print *,  EXE_PROMPT, "Cloud mask results will be read in from precomputed 1bx file"
    endif

    if (bx_file_Flag == 1) then
         if (nav_Flag == 1) then
            print *,  EXE_PROMPT, "Cannot create 1bx file and use external re-navigation"
            stop "3a"
         else
            print *, EXE_PROMPT, "1bx file will be created"
         endif
    endif


    if ((Cloud_Mask_Aux_Flag /= sym%NO_AUX_CLOUD_MASK) .and. (bx_file_Flag == sym%YES)) then
       erstat = 6
       print *,  EXE_PROMPT, "Cannot read cloud mask from both 1bx and write out a 1bx file, stopping"
       stop 6
    endif

    if (Rtm_Flag /=1) then
       print *,  EXE_PROMPT, "Only PFAST RTM implemented, stopping"
       stop
    endif


    print *, EXE_PROMPT, "Temporary Files will be written to ", trim(Temporary_Data_Dir)

 end subroutine QC_CLAVRXORB_OPTIONS
 
 
!--------------------------------------------------------------------------
! This subroutine outputs what each command line options
! are available to override the default file options or to set
! options different from the default options
!--------------------------------------------------------------------------
subroutine HELPER()
 
  print "(a,'help: option list')",EXE_PROMPT
  print *,"  This is a list of all of the command line options to override"
  print *,"     the file list options"
  print *," "

  print *,"  -default (default_temp)"
  print *, "  This option allows you to set which default file you want to use."
  print *, "  Initial default file is clavrx_Default_options"
  print *," "

  print *,"  -filelist (file_list)"
  print *, "  This option allows you to set which clavrxorb_file_list file you want to use."
  print *, "  Initial default file is the clavrxorb_file_list in the current directory"
  print *," "

  print *,"  -lines_per_seg (Num_Scans_per_Segment)"
  print *, "  specify the number of scanlines per segment (default is 240)"
  print *," "
 
  print *,"  -ref_Cal_1b"
  print *,"  Use the reflectance cal in level 1b file."
  print *," "

  print *,"  -no_ref_Cal_1b"
  print *,"  Do not us the reflectance cal in level 1b file."
  print *," "

  print *,"  -therm_Cal_1b"
  print *,"   Use the thermal calibration in the level 1b file"
  print *," "

  print *,"  -no_therm_Cal_1b"
  print *,"   Do not use the thermal calibration in the level 1b file"
  print *," "
  
  print *,"  -bx_file_Flag"
  print *, "  Fill in clavr-x bytes in the level 1b file"
  print *," "

  print *,"  -no_bx_file_Flag"
  print *, "  Do not fill in clavr-x bytes in the level 1b file"
  print *," "
  
  print *,"  -l1bnav"
  print *, "  Use the navigation data from level 1b file"
  print *," "

  print *,"  -clevernav"
  print *, "  Use the navigation data from a clevernav"
  print *," "

  print *,"  -nav_file"
  print *, "   Output a navigation file."
  print *," "
  print *,"  -no_nav_file"
  print *, "   Do not output a navigation file."
  print *," "

  print *,"  -cmr_file"
  print *, "   Output cmr file."
  print *," "
  print *,"  -no_cmr_file"
  print *, "   Do not output cmr file."
  print *," "
  
  print *,"  -obs_file"
  print *, "   Output obs file."
  print *," "
  print *,"  -no_obs_file"
  print *, "   Do not output obs file."
  print *," "
  
  print *,"  -geo_file"
  print *, "   Output geolocation file. "
  print *," "
  print *,"  -no_geo_file"
  print *, "   Output geolocation file. "
  print *," "

  print *,"  -cld_file"
  print *, "   Run and output cloud algorithms."
  print *," "
  print *,"  -no_cld_file"
  print *, "   Do not output cloud algorithms."
  print *," "

  print *,"  -sst_file"
  print *, "   Calculate and output SST data. "
  print *," "
  print *,"  -no_sst_file"
  print *, "   Do not output SST data. "
  print *," "

  print *,"  -ash_file"
  print *, "   Output VolAsh data file. "
  print *," "
  print *,"  -no_ash_file"
  print *, "   Do not output VolAsh data file. "
  print *," "

  print *,"  -rtm_file"
  print *, "   Output rtm data file. "
  print *," "
  print *,"  -no_rtm_file"
  print *, "   Do not output rtm data file. "
  print *," "

  print *,"  -level2_file"
  print *, "   Make Level-2 output."
  print *," "
  print *,"  -no_level2_file"
  print *, "   Don't make Level-2 output."
  print *," "

  print *,"  -level3_file"
  print *, "   Make level3 output."
  print *," "
  print *,"  -no_level3_file"
  print *, "   Don't make level3 output."
  print *," "
  
  print *,"  -cloud_mask_1b"
  print *, "   Read cloud mask from level 1b file."
  print *," "
  print *,"  -cloud_mask_calc"
  print *, "   recalculate cloud mask."
  print *," "
  
  print *,"  -cld_Flag"
  print *, "   Run cloud algorithms. "
  print *," "
  print *,"  -no_cld_Flag"
  print *, "   Don't run cloud algorithms. "
  print *," "
  
  print *,"  -aer_Flag"
  print *, "   Run aerosol algorithms. "
  print *," "
  print *,"  -no_aer_Flag"
  print *, "   Don't run aerosol algorithms. "
  print *," "
  
  print *,"  -erb_Flag"
  print *, "   Run radiative flux algorithms. "
  print *," "
  print *,"  -no_erb_Flag"
  print *, "   Don't run radiative flux algorithms. "
  print *," "

!  print *,"  -ash_Flag"
!  print *, "   Run Vol Ash algorithms. "
!  print *," "
!  print *,"  -no_ash_Flag"
!  print *, "   Don't run Vol Ash algorithms. "
!  print *," "

  print *,"  -use_sst_anal"
  print *, "   Use OISST analysis. "
  print *," "
  print *,"  -no_use_sst_anal"
  print *, "   Don't use OISST analysis. "
  print *," "
  
  print *,"  -use_oisst_cur"
  print *, "   Use the oisst.current OISST analysis."
  print *," "
  print *,"  -determine_oisst_file"
  print *, "   Determine which OISST analysis file to use."
  print *," "
  
  print *,"  -subset_pixel_hdf"
  print *, "  Subset HDF pixel data. " 
  print *," "
  print *,"  -no_subset_pixel_hdf"
  print *, "  Don't subset HDF pixel data. " 
  print *," "
  
  print *,"  -no_nwp"
  print *, "  Do not use NWP data. No algorithms will be processed. Also, no orbital data will be processed"
  print *," "

  print *,"  -gfs_nwp"
  print *, "  Use GFS data for NWP dataset. "
  print *," "
  
  print *,"  -ncep_nwp"
  print *, "  Use NCEP Reanalysis data"
  print *," "
  
  print *,"  -crtm"
  print *, "  Use CRTM. "
  print *, "  NOT AVAILABLE CURRENTLY. ONLY PFAST AVAILABLE" 
  print *," "
  
  print *,"  -pfast"
  print *, "  Use PFAST RTM. "
  print *," "
  
  print *,"  -use_clear_composite"
  print *,"    Use this option to use clear composite refs."
  print *," "
  print *,"  -no_clear_composite"
  print *,"    Use this option to NOT use clear composite refs."
  print *," "

  print *,"  -use_modis_clr_alb"
  print *,"    Use this option to use MODIS clear reflectance data."
  print *," "
  print *,"  -no_modis_clr_alb"
  print *,"    Use this option to NOT use MODIS clear reflectance data."
  print *," "
 
  print *,"  -prob_clear_res"
  print *, "  Restore probably clear pixels. "
  print *," "
  print *,"  -no_prob_clear_res"
  print *, "  Don't restore probably clear pixels. "
  print *," "

  print *,"  -lrc"
  print *, "  Include local radiative center calcs.  "
  print *," "
  print *,"  -no_lrc"
  print *, "  Don't include local radiative center calcs.  "
  print *," "
  
  print *,"  -crc"
  print *, "  Include cloud radiative center calcs.  "
  print *," "
  print *,"  -no_crc"
  print *, "  Don't include cloud radiative center calcs.  "
  print *," "

  print *,"  -lat_min_limit (Lat_Min_Limit)"
  print *, "  minimum latitude for processing(degrees)"
  print *," "
  print *,"  -lat_max_limit (Lat_Max_Limit)"
  print *, "  maximum latitude for processing (degrees)"
  print *," "

  print *,"  -ancil_data_dir (ancil_data_dir)"
  print *, "  change the location of the ancillary data directory"
  print *," "
  print *,"  -gfs_data_dir (gfs_data_dir)"
  print *, "  change the location of the gfs data directory"
  print *," "
  print *,"  -ncep_data_dir (ncep_data_dir)"
  print *, "  change the location of the ncep-reanalysis data directory"
  print *," "
  print *,"  -cfsr_data_dir (cfsr_data_dir)"
  print *, "  change the location of the cfsr data directory"
  print *," "
  print *,"  -oisst_data_dir (oisst_data_dir)"
  print *, "  change the location of the oisst data directory"
  print *," "
  print *,"  -snow_data_dir (snow_data_dir)"
  print *, "  change the location of the snow data directory"
  print *," "
  
  print *,"  -smooth_nwp"
  print *, "  Smooth the NWP fields. "
  print *," "
  print *,"  -no_smooth_nwp"
  print *, "  Don't smooth the NWP fields. "
  print *," "
  
  print *,"  -use_seebor"
  print *, "  Use Seebor emissivity dataset. "
  print *," "
  print *,"  -no_seebor"
  print *, "  Don't use Seebor emissivity dataset. "
  print *," "
  
  print *,"  -high_res"
  print *, "  Use high resolution surface type map. "
  print *," "
  print *,"  -low_res"
  print *, "  Use low resolution (8km) surface type map. "
  print *," "
  
  print *,"  -read_land_mask"
  print *, "  Read land mask in. "
  print *," "
  print *,"  -no_land_mask"
  print *, "  Don't read land mask in. "
  print *," "
  
  print *,"  -read_coast_mask"
  print *, "  Read coast mask in. "
  print *," "
  print *,"  -no_coast_mask"
  print *, "  Don't read coast mask in. "
  print *," "
  
  print *,"  -no_surface_elevation"
  print *, "  Don't read surface elevation map in. "
  print *," "
  print *,"  -high_surface_elevation"
  print *, "  Use 1km res GLOBE global elevation map. "
  print *," "
  print *,"  -low_surface_elevation"
  print *, "  Use 8km res GLOBE global elevation map. "
  print *," "

  
  print *,"  -read_volcano_mask"
  print *, "  Read volcano map in. "
  print *," "
  print *,"  -no_volcano_mask"
  print *, "  Don't read volcano map in. "
  print *," "
  
  print *,"  -Solzen_min_limit (Solzen_min_limit)"
  print *, "  Solar zenith angle minimum limit"
  print *," "
  
  print *,"  -Solzen_max_limit (Solzen_max_limit)"
  print *, "  Solar zenith angle maximum limit."
  print *," "  
  
  print *,"  -dcomp_mode"
  print *, "  dcomp mode 1,2 or 3 (1.6,2.2 or 3.8 micron)."
  print *," "  

  end subroutine HELPER

end module USER_OPTIONS
