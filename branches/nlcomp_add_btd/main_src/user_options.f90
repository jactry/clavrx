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
!             Public Routines in this module and their purpose:
!               SETUP_USER_DEFINED_OPTIONS:
!                    Reads user defined configuration
!           	      Called in process_clavrx.f90 outside any loop in the beginning
!
!                UPDATE_CONFIGURATION
!                     Updates configuarion for each file
!                     called in process_clavrx.f90 inside file loop
!                     Check algorithm modes and channel switches
!       
!
! AUTHORS:
!  	Andrew Heidinger, Andrew.Heidinger@noaa.gov
!   	Andi Walther
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
!
!  HISTORY:
!        16 Dec 2014: Switch to new option file  (AW)
!                       code cleaning
!
!         30 Dec 2014: Submitting to trunk

!
!   
!--------------------------------------------------------------------------------------
module USER_OPTIONS

   use PIXEL_COMMON, only: &
      chan_on_flag_default &
      , aer_flag &
      , ash_file_flag &
      , ash_flag &
      , goes_stride &
      
      , cld_flag &
      , cloud_mask_aux_flag &
      , cloud_mask_Bayesian_flag &
 
      , sasrab_flag &
      , nav_opt &
      , nwp_opt &
      , ref_cal_1b &
      , rtm_file_flag &
      , rtm_opt &
      , sst_file_flag &
      , temporary_data_dir &
      , therm_cal_1b &
      , ancil_data_dir &
      , cfsr_data_dir &
      , data_comp_flag &
  
      , file_list &
      , geo_file_flag &
      , gfs_data_dir &
      , lat_max_limit &
      , lat_min_limit &
      , level2_file_flag &
      , lrc_flag &
      , modis_clr_alb_flag &
      , ncep_data_dir &
      , num_scans_per_segment &
      , obs_file_flag &
      , oisst_data_dir &
      , output_scaled_reflectances &
      
      , process_undetected_cloud_flag &
      , read_coast_mask &
      , read_hires_sfc_type &
      , read_land_mask &
      , read_snow_mask &
      , read_surface_elevation &
      , read_volcano_mask &
      , smooth_nwp_flag &
      , snow_data_dir &
      , solzen_max_limit &
      , solzen_min_limit &
     
      , subset_pixel_hdf_flag &
      , use_default &
      , use_seebor &
      
      , acha_mode &
      , bayesian_cloud_mask_name &
      , Compress_Flag &
      , Nlcomp_Mode &
      , viirs_flag &
      , goes_mop_flag &
      , dcomp_mode &
      , avhrr_flag &
      , avhrr_1_flag &
      
      , goes_flag &
     
      , read_dark_comp &
      , mtsat_flag &
      , sc_id_wmo
      
      
   use CONSTANTS, only: &
      sym &
      , exe_prompt   
      
   use FILE_UTILITY, only: &
      get_lun
 
   use  clavrx_message_module, only: &
      mesg &
    , verb_lev

   implicit none
   private
   public  :: SETUP_USER_DEFINED_OPTIONS
   public  :: UPDATE_CONFIGURATION

   character(24), parameter, private :: MOD_PROMPT = " USER_OPTIONS_ROUTINES: "
   character ( len = 50 ) :: data_base_path
   integer :: Dcomp_Mode_User_Set
   integer :: Acha_Mode_User_set
   integer :: nlcomp_Mode_User_set
   integer :: mask_Mode_User_set
   integer :: expert_mode
   
   integer :: chan_on_flag_default_user_set ( 42)
   
   
   
   ! ---------------------------------------------------------------------------------
   ! Default Algorithm Modes - 
   ! ---------------------------------------------------------------------------------
   integer,parameter:: ACHA_Mode_Default_Avhrr = 3
   integer,parameter:: ACHA_Mode_Default_Avhrr1 = 1
   integer,parameter:: ACHA_Mode_Default_Goes_IL = 6
   integer,parameter:: ACHA_Mode_Default_Goes_MP = 7
   integer,parameter:: ACHA_Mode_Default_Goes_SNDR = 7
   integer,parameter:: ACHA_Mode_Default_COMS = 7
   integer,parameter:: ACHA_Mode_Default_VIIRS = 5
   integer,parameter:: ACHA_Mode_Default_MTSAT = 6
   integer,parameter:: ACHA_Mode_Default_SEVIRI = 8
   integer,parameter:: ACHA_Mode_Default_Modis = 8
   integer,parameter:: ACHA_Mode_Default_Fy2 = 7
   
contains

   ! ---------------------------------------------------------------------------------
   !  wrapper for initial clavrx option read
   ! ---------------------------------------------------------------------------------
   subroutine SETUP_USER_DEFINED_OPTIONS()
      call SET_DEFAULT_VALUES
      call DETERMINE_USER_CONFIG()
      call QC_CLAVRXORB_OPTIONS()

   end subroutine SETUP_USER_DEFINED_OPTIONS
   
   ! ---------------------------------------------------------------------------------
   !
   ! ---------------------------------------------------------------------------------
   subroutine set_default_values
           
      Aer_Flag = sym%YES
      Ash_Flag = sym%NO
      Data_comp_Flag = 0
      Subset_pixel_hdf_Flag = 0      
      modis_clr_alb_Flag = 1 ! do not use clear-sky MODIS albedo maps
      output_scaled_reflectances = sym%NO !default is to output ref / cosSolzen
      
      !--- default solar zenith limits
      Solzen_Min_Limit= 0 
      Solzen_Max_Limit= 180.0
       
      !  -- default what can be changed for expert mode
      Cloud_Mask_Bayesian_Flag = 1
      Dcomp_Mode_user_set = 3
      Acha_Mode_user_set = 1
      Nlcomp_Mode = 1            
      Level2_File_Flag = 1
      Rtm_File_Flag = 1
      Cld_Flag = 1
      num_scans_per_segment = 240
      Sasrab_Flag = 0
      Nwp_Opt = 1
      Rtm_Opt = 1 
      Compress_Flag = 1 
      Cloud_Mask_Aux_Flag = 0 
      bayesian_cloud_mask_name = 'default'
      Use_Seebor = 1 
      Read_Hires_Sfc_Type = 1 
      Read_Land_Mask = 1
      Read_Coast_Mask = 1
      Read_Surface_Elevation = 1
      Read_Volcano_Mask = 0
      Read_Snow_Mask = 1
      Read_Dark_Comp = 0
      Ref_Cal_1b = 1
      Therm_Cal_1b = 1    
      Nav_Opt = 0  
      goes_stride = 2
      Lrc_Flag = 1 
      Smooth_Nwp_Flag = 1  
      Process_Undetected_Cloud_Flag = 0
      Chan_On_Flag_Default_user_set(1:6) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_user_set(7:12) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_user_set(13:18) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_user_set(19:24) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_user_set(25:30) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_user_set(31:36) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_user_set(37:42) = [0,0,0,0,0,1]
      
   end subroutine
   

   !---------------------------------------------------------------------------------
   ! Read Parameters from AVHRR INPUT files and check them for errors
   !
   ! parameters are pased in avhrr_pixel_common public memory
   !---------------------------------------------------------------------------------
   subroutine READ_OPTION_FILE (File_Default)
      character(len=*), intent(in):: File_Default
      integer::ios0
      integer::erstat
      integer:: Default_Lun
            
      call mesg ("DEFAULT FILE READ IN",level = 5 )
      call mesg ("Default file to be read in: "//trim(File_Default),level = verb_lev % DEFAULT)

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
      end if
      
      read(unit=Default_Lun,fmt="(a)") Data_base_path
      read(unit=Default_Lun,fmt="(a)") Temporary_Data_Dir
      read(unit=Default_Lun,fmt=*) expert_mode
         
      if ( expert_mode  == 0 )  then
          close(unit=Default_Lun)
          return
      end if
      
      read(unit=Default_Lun,fmt=*) Cloud_Mask_Bayesian_Flag
      read(unit=Default_Lun,fmt=*) Dcomp_Mode_user_set
      read(unit=Default_Lun,fmt=*) Acha_Mode_user_set
      read(unit=Default_Lun,fmt=*) Nlcomp_Mode
         
      if ( expert_mode <= 1 )  then
          close(unit=Default_Lun)
          return
      end if
      
      read(unit=Default_Lun,fmt=*) Level2_File_Flag
      read(unit=Default_Lun,fmt=*) Rtm_File_Flag
      read(unit=Default_Lun,fmt=*) Cld_Flag
      read(unit=Default_Lun,fmt=*) num_scans_per_segment
      read(unit=Default_Lun,fmt=*) Sasrab_Flag
      read(unit=Default_Lun,fmt=*) Nwp_Opt
     
      read(unit=Default_Lun,fmt=*) Rtm_Opt
      read(unit=Default_Lun,fmt=*) Nav_Opt
      
      read(unit=Default_Lun,fmt=*) Compress_Flag
      read(unit=Default_Lun,fmt=*) Cloud_Mask_Aux_Flag

      read(unit=Default_Lun,fmt="(a)") bayesian_cloud_mask_name
      
      if ( expert_mode <= 2 )  then
          close(unit=Default_Lun)
          return
      end if 
      
      read(unit=Default_Lun,fmt=*) Use_Seebor
      read(unit=Default_Lun,fmt=*) Read_Hires_Sfc_Type
      read(unit=Default_Lun,fmt=*) Read_Land_Mask
      read(unit=Default_Lun,fmt=*) Read_Coast_Mask
      read(unit=Default_Lun,fmt=*) Read_Surface_Elevation
      read(unit=Default_Lun,fmt=*) Read_Volcano_Mask
      read(unit=Default_Lun,fmt=*) Read_Snow_Mask
      read(unit=Default_Lun,fmt=*) Read_Dark_Comp
         
      if ( expert_mode <= 3 ) then
          close(unit=Default_Lun)
          return
      end if
         
      read(unit=Default_Lun,fmt=*) Ref_Cal_1b
      read(unit=Default_Lun,fmt=*) Therm_Cal_1b
            
      if ( expert_mode <= 4 ) then
          close(unit=Default_Lun)
          return
      end if
      
      read(unit=Default_Lun,fmt=*) Lrc_Flag
      read(unit=Default_Lun,fmt=*) Smooth_Nwp_Flag
      read(unit=Default_Lun,fmt=*) Process_Undetected_Cloud_Flag
               
      if ( expert_mode <= 5 ) then
          close(unit=Default_Lun)
          return
      end if
      
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_user_set(1:6)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_user_set(7:12)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_user_set(13:18)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_user_set(19:24)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_user_set(25:30)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_user_set(31:36)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_user_set(37:42)
             
     
      close(unit=Default_Lun)
    

   end subroutine READ_OPTION_FILE
 
   subroutine DETERMINE_USER_CONFIG()
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
      integer :: iargc

      temp_string = '.'  !--- tempoary string to search for in angle commandline

      !---- SET DEFAULT OPTIONS
  
      use_Default = sym%YES
      default_temp="./clavrx_options"
      file_list = "./file_list"
  

  
      Temp_Scans_Arg = 0
      fargc = iargc()
  
      !--- first we will check to see if the default file is used
      !--- also check to see if the help file is to be displayed
  
      do i=1, fargc
         call getarg(i,fargv)
         if (trim(fargv) == "-help" .or. &
               trim(fargv) == "-h") then
            call HELPER()
            stop   
         else if  ( trim(fargv) == "-version" .or. &
             trim(fargv) == "-ver") then
            print*,&
            &'$Header$'
            stop
        
         else if (trim(fargv) == "-no_Default") then 
            use_Default = sym%NO
            !Different default file used
         else if (trim(fargv) == "-default") then
            call getarg(i+1,default_temp)
            default_temp=trim(default_temp)
         end if
      end do
  
  
      !---- If the default file is used, read it in first, then
      !---- check for other command line changes to the options
  

      if(Use_Default == sym%YES)  then
         call READ_OPTION_FILE (Default_Temp)
      end if
      
      if(Use_Default == sym%NO)  then
          print *, EXE_PROMPT, "Using standard defaults and command line options"
      end if 
      
      
      
 
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
          nav_opt = 0
          
        !Change RTM output flag
        elseif(trim(fargv) == "-rtm_file") then
          rtm_file_Flag = sym%YES
        elseif(trim(fargv) == "-no_rtm_file") then
          rtm_file_Flag = sym%NO

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

        !Smooth/not smooth NWP data
        elseif(trim(fargv) == "-smooth_nwp") then
          Smooth_Nwp_Flag = sym%YES
        elseif(trim(fargv) == "-no_smooth_nwp") then
          Smooth_Nwp_Flag = sym%NO

        !Read/not read volcano mask
        elseif(trim(fargv) == "-read_volcano_mask") then
          read_volcano_mask = sym%YES
        elseif(trim(fargv) == "-no_volcano_mask") then
          read_volcano_mask = sym%NO

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
  

      !--- default ancillary data directory
      ancil_data_dir = trim(data_base_path)//'/clavrx_ancil_data/'
      gfs_data_dir = trim(data_base_path)//'gfs/'
      ncep_data_dir = trim(data_base_path)//'/clavrx_ancil_data/ncep-reanalysis/'
      cfsr_data_dir = trim(data_base_path)//'/cfsr/'
      oisst_data_dir = trim(data_base_path)//'/clavrx_ancil_data/oisst/'
      snow_data_dir = trim(data_base_path)//'/snow/'
      
      
      chan_on_flag_default = chan_on_flag_default_user_set

   end subroutine DETERMINE_USER_CONFIG

   !-------------------------------------------------------------------------------
   !--- QC options and modify as needed
   !-------------------------------------------------------------------------------
   subroutine QC_CLAVRXORB_OPTIONS()

      integer:: erstat

      !---- Since the NWP controls everything, we first check if an NWP is being used
      !---- before anything else is checked.  If no nwp, we stop processing

      select case ( nwp_opt)
      case ( 0 )   
         print *,  EXE_PROMPT, "No choice made for NWP data, will not run algoritms or orbital level3 files"
         Cld_Flag = sym%NO
       
         Sasrab_Flag = sym%NO
         Rtm_File_Flag = sym%NO
         Cloud_Mask_Bayesian_Flag = sym%NO
         Cloud_Mask_Aux_Flag = sym%NO ! this is to determine if the lut's are being read in
      case ( 1 )
         call mesg ("GFS data will be used",level = verb_lev % DEFAULT)
      case ( 2 )
         call mesg ( "NCEP Reanalysis data will be used",level = verb_lev % DEFAULT)
      case ( 3 )
         call mesg ( "NCEP Climate Forecast System Reanalysis data will be used",level = verb_lev % DEFAULT)
      case ( 4 )
         call mesg ( "GDAS Reanalysis data will be used",level = verb_lev % DEFAULT)
      case default
         print *,  EXE_PROMPT, "unrecognized value for Nwp_Opt: ", Nwp_Opt
         stop "6-Nwp_Flag"
      
      end select
      
     
      if (cloud_mask_bayesian_Flag == sym%YES) then
         call mesg  ("Bayesian cloud mask will be generated")
      endif

      if (Ref_Cal_1b == sym%YES) then
         call mesg ("Reflectance Calibration within 1b will be used")
      endif

      if (therm_Cal_1b == sym%YES) then
         call mesg ("Thermal Calibration within 1b will be used")
      endif

      if (nav_Opt == 1) then
         call mesg ("CLEVERNAV geolocation no longer supported, using REPOSNX")
         nav_Opt = 2
       endif

      if (nav_opt == 2) then
         call mesg( "REPOSNX geolocation adjustment done")
      endif

      if (rtm_file_Flag == sym%YES) then
        call mesg( "rtm file will be created")
      endif

      if (Cloud_Mask_Aux_Flag == sym%YES) then
         print *,  EXE_PROMPT, "Cloud mask results will be read in from an aux file"
      endif

      if (Rtm_opt /=1) then
         print *,  EXE_PROMPT, "Only PFAST RTM implemented, switching to rtm ==1"
         rtm_opt = 1
      endif


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
      print *, "  Initial default file is clavrx_options"
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
      print *,"  AVHRR-ONLY: Do not us the reflectance cal in level 1b file."
      print *," "

      print *,"  -therm_Cal_1b"
      print *,"   AVHRR-ONLY: Use the thermal calibration in the level 1b file"
      print *," "

      print *,"  -no_therm_Cal_1b"
      print *,"   Do not use the thermal calibration in the level 1b file"
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
  
      print *,"  -cloud_mask_1b"
      print *, "   Read cloud mask from level 1b file."
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

      print *,"  -lat_min_limit (Lat_Min_Limit)"
      print *, "  minimum latitude for processing(degrees)"
      print *," "
      print *,"  -lat_max_limit (Lat_Max_Limit)"
      print *, "  maximum latitude for processing (degrees)"
      print *," "

      print *,"  -ancil_data_dir (ancil_data_dir)"
      print *, "  change the location of the ancillary data directory"
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
  
   ! -----------------------------------------------------------------
   ! -- wrapper for all updating tools for a new file
   !     called from PROCESS_CLAVRX inside file loop
   ! ---------------------------------------------------
   subroutine UPDATE_CONFIGURATION ( sensorname)
      character ( len =10) , intent(in) :: sensorname
      
      if ( expert_mode == 0 ) then
         acha_mode_user_set =  default_acha_mode ( sensorname )
         dcomp_mode_user_set = default_dcomp_mode ( sensorname )
      end if

      call CHECK_ALGORITHM_CHOICES(sensorname)
     
      call CHANNEL_SWITCH_ON (sensorname)

      if ( expert_mode < 3 .or. trim(bayesian_cloud_mask_name) == 'default') then
         bayesian_cloud_mask_name = default_nb_mask_classifier_file ( sensorname )
      end if
      
      call EXPERT_MODE_CHANNEL_ALGORITHM_CHECK ( sensorname ) 
      
      
      
   end subroutine UPDATE_CONFIGURATION
   
   !----------------------------------------------------------------------
   !  returns default acha mode 
   !----------------------------------------------------------------------
   integer function default_acha_mode ( sensorname )
      character ( len =10) , intent(in) :: sensorname
      
      
      
      select case ( trim(sensorname))
      
      case ( 'AVHRR')        
         default_acha_mode  = ACHA_Mode_Default_Avhrr
      case ( 'AVHRR_1')   
         default_acha_mode = ACHA_Mode_Default_Avhrr1   
      case ( 'GOES_MOP')      
         default_acha_mode  = ACHA_Mode_Default_Goes_MP
      case ( 'GOES')      
         default_acha_mode  = ACHA_Mode_Default_Goes_IL   
      case ( 'GOES_SNDR')
         default_acha_mode  = ACHA_Mode_Default_Goes_SNDR 
      case ( 'MTSAT')
          default_acha_mode  = ACHA_Mode_Default_MTSAT
      case ('SEVIRI')
         default_acha_mode  = ACHA_Mode_Default_SEVIRI
      case ('FY2')
         default_acha_mode  =  ACHA_Mode_Default_FY2 
      case ('VIIRS')
         default_acha_mode  =  ACHA_Mode_Default_VIIRS 
      case ('IFF_VIIRS')      
          default_acha_mode  = ACHA_Mode_Default_VIIRS      
      case ('IFF_AVHRR')      
         default_acha_mode  = ACHA_Mode_Default_Avhrr
      case ('COMS')
         default_acha_mode  = ACHA_Mode_Default_COMS
      case ('MODIS')
          default_acha_mode  = ACHA_Mode_Default_Modis 
      case ('MODIS_1KM')
          default_acha_mode  = ACHA_Mode_Default_Modis
      case default 
         print*,'sensor ',sensorname, ' is not set in check channels settings Inform andi.walther@ssec.wisc.edu'   
      end select
          
   end function default_acha_mode
   
   !-----------------------------------------------------------------
   !   returns default dcomp mode
   ! -----------------------------------------------------------------
   integer function default_dcomp_mode ( sensorname )
      character ( len =10) , intent(in) :: sensorname
   
      default_dcomp_mode = 3
      
      if (Sc_Id_WMO == 208) Dcomp_Mode = 1   !NOAA-17
      if (Sc_Id_WMO == 3) Dcomp_Mode = 1     !METOP-A
      if (Sc_Id_WMO == 4) Dcomp_Mode = 1     !METOP-B
      if (Sc_Id_WMO == 5) Dcomp_Mode = 1     !METOP-C
     
   end function default_dcomp_mode

!-----------------------------------------------------------------
!   returns default classifier name
!-----------------------------------------------------------------
   
   function default_nb_mask_classifier_file ( sensorname ) result (filename)
      character ( len =10) , intent(in) :: sensorname
      character ( len = 355 ) :: filename

      select case ( trim(sensorname))
      
      case ( 'AVHRR')        
         filename  = 'avhrr_default_nb_mask.nc'
      case ( 'AVHRR_1')   
         filename  = 'avhrr_default_nb_mask.nc'   
      case ( 'GOES_MOP')      
         filename  = 'goesnp_default_nb_mask.nc'
      case ( 'GOES')      
         filename  = 'goesim_default_nb_mask.nc'   
      case ( 'GOES_SNDR')
         filename  = 'goesnp_default_nb_mask.nc' 
      case ( 'MTSAT')
          filename  = 'mtsat_default_nb_mask.nc'
      case ('SEVIRI')
         filename  = 'avhrr_default_nb_mask.nc'
      case ('FY2')
         filename  = 'avhrr_default_nb_mask.nc' 
      case ('VIIRS')
         filename  = 'viirs_default_nb_mask.nc' 
      case ('IFF_VIIRS')      
          filename  = 'avhrr_default_nb_mask.nc'
      case ('IFF_AVHRR')      
        filename  = 'avhrr_default_nb_mask.nc'
      case ('COMS')
         filename  = 'avhrr_default_nb_mask.nc'
      case ('MODIS')
          filename  = 'modis_default_nb_mask.nc' 
      case ('MODIS_1KM')
          filename  = 'modis_default_nb_mask.nc'
      case default 
         print*,'sensor ',sensorname, ' is not set in check channels settings Inform andi.walther@ssec.wisc.edu'  
         stop 
      end select


     end function default_nb_mask_classifier_file

   !----------------------------------------------------------------------------
   !  check if algo mode set by user is possible
   !----------------------------------------------------------------------------
   
   subroutine CHECK_ALGORITHM_CHOICES(sensorname)
      character ( len =10) , intent(in) :: sensorname
      character ( len = 1 ) :: string_1
      
      integer :: possible_acha_modes ( 8 )
      integer :: possible_dcomp_modes ( 3)
 
      !------------------------------------------------------------------------
      !--- ACHA MODE Check
      !---      (       0 = off
      !---              1 = 11
      !----             2 = 11/6.7; 
      !---              3 = 11/12
      !----             4 = 11/13.4; 
      !---              5 = 11/12/8.5
      !---              6 = 11/12/6.7
      !---              7 = 11/13.3/6.7
      !---              8 = 11/12/13.3)
      !------------------------------------------------------------------------
      
      acha_mode = acha_mode_user_set     
      dcomp_mode = dcomp_mode_user_set
       
      possible_acha_modes = 0 
      possible_dcomp_modes = 0
         
      select case ( trim ( sensorname))
      
      case ( 'AVHRR')  
         possible_acha_modes(1:2)   = [1, 3]
         possible_dcomp_modes(1)    =  3 
      case ( 'AVHRR_1')   
         possible_acha_modes(1)     =  1
         possible_dcomp_modes(1)    =  3   
      case ( 'GOES_MOP')      
         possible_acha_modes(1:4)   =  [1, 2, 4, 7]
         possible_dcomp_modes(1)    =  3
      case ( 'GOES')      
         possible_acha_modes(1:4)   =  [1, 2, 3, 6] 
         possible_dcomp_modes(1)    =  3  
      case ( 'GOES_SNDR')
         possible_acha_modes(1:8)   =  [1, 2, 3, 4, 5, 6, 7, 8]  
         possible_dcomp_modes(1)    =  3
      case ( 'MTSAT')
          possible_acha_modes(1:4)  =  [4, 5, 7, 8]
          possible_dcomp_modes(1:2) =  [1,3]
      case ('SEVIRI')
         possible_acha_modes(1:8)  =   [1, 2, 3, 4, 5, 6, 7, 8]
         possible_dcomp_modes(1:2) =   [1, 3]
      case ('FY2')
         possible_acha_modes(1:2)  =   [1 , 2 ] 
         possible_dcomp_modes(1:2) =   [1 , 3 ]
      case ('VIIRS')
        possible_acha_modes(1:4)  =    [1, 3, 4, 5] 
        possible_dcomp_modes(1:3) =    [1, 2, 3]
         nlcomp_mode_user_set = 1  
      case ('IFF_VIIRS')      
         possible_acha_modes(1:4)  =   [1, 3, 5, 8]      
      case ('IFF_AVHRR')      
         possible_acha_modes(1:4)  =   [1, 3, 4, 8]
      case ('COMS')
         possible_acha_modes(1:3)  =   [1, 3, 6]
         possible_dcomp_modes(1:2) =   [1, 3]
      case ('MODIS')
         possible_acha_modes(1:8)  =   [1, 2, 3, 4, 5, 6, 7, 8] 
         possible_dcomp_modes(1:3) =   [1, 2, 3]
      case ('MODIS_1KM')
         possible_acha_modes(1:8)  =   [1, 2, 3, 4, 5, 6, 7, 8]
         possible_dcomp_modes(1:3) =   [1, 2, 3]
      case default 
         print*,'sensor ',sensorname, ' is not set in chaeck channels settings Inform andi.walther@ssec.wisc.edu'   
      end select
      
      if ( .not. ANY ( acha_mode_user_set == possible_acha_modes ) ) then
         acha_mode = default_acha_mode ( sensorname )
         print*,acha_mode_user_set,possible_acha_modes
         print*, 'User set ACHA mode not possible for '//trim(sensorname)//' switched to default '
      end if
 
      if ( .not. ANY ( dcomp_mode_user_set == possible_dcomp_modes ) ) then
         dcomp_mode = default_dcomp_mode ( sensorname )
         print*, 'User set DCOMP mode not possible for '//trim(sensorname)//' switched to default '
      end if

   end subroutine CHECK_ALGORITHM_CHOICES
   
   ! ----------------------------------------------------------------------
   !    returns all available sensors for this sensors
   ! ----------------------------------------------------------------------
   function existing_channels  (sensorname)  result( valid_channels )
      character ( len = 10 ) , intent(in) :: sensorname
      
      integer , target :: valid_channels ( 42) 
     
      
      valid_channels = -99
      select case ( trim(sensorname))
        
      case ( 'AVHRR')
         valid_channels (1:6) = [1,2,6,20,31,32]
      case ( 'GOES')      
         valid_channels (1:5) = [1,20,27, 31, 32]
      case ( 'GOES_MOP')      
         valid_channels (1:5) = [1,20,27, 31, 33]   
      case ( 'GOES_SNDR')
         valid_channels (1:18) = [1,20,21,23,24,25,30,31,32,33,34,35,36,37,38,39,40,41]      
      case ( 'MTSAT')
         valid_channels (1:5) = [1,20,27,31,32]  
      case ('SEVIRI')
         valid_channels (1:11) = [1,2,6,20,27,28,29,30,31,32,33]
      case ('FY2')
         valid_channels (1:5) = [1,20,27,31,32]    
      case ('VIIRS')
         valid_channels (1:22) = [1,2,3,4,5,6,7,8,9,15,20,22,26,29,31,32,37,38,39,40,41,42]       
      case ('IFF_VIIRS')      
         valid_channels (1:26) = [1,2,3,4,5,6,7,8,9,15,20,22,26,29,31,32,33,34,35,36,37,38,39,40,41,42]       
      case ('IFF_AVHRR')      
         valid_channels (1:19) = [1,2,6,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36]
      case ('COMS')
         valid_channels (1:5) = [1,20,27,31,32]
      case ('MODIS')
         valid_channels(1:12) = [1,2,6,7,8,20,26,27,29,31,32,33]   
      case ('MODIS_1KM')
         valid_channels(1:12) = [1,2,6,7,8,20,26,27,29,31,32,33]   
      case default 
         print*,'sensor ',sensorname, ' is not set in check channels settings Inform andi.walther@ssec.wisc.edu'   
      end select


   end function existing_channels
   
   !----------------------------------------------------------------------
   !   Channel settings
   !     will not be done for full-experts  ( expert mode 7 and higher)
   !
   !----------------------------------------------------------------------
   subroutine CHANNEL_SWITCH_ON (sensorname)
      character ( len =10) , intent(in) :: sensorname
      integer :: valid_channels ( 42)
      integer :: i
 
      ! expert can decide themselves
      if (expert_mode > 6 ) return
      
      valid_channels = existing_channels ( sensorname )
           
      Chan_On_Flag_Default =  0

      do i = 1, 42 
         if (valid_channels (i) < 0 ) cycle
         Chan_On_Flag_Default (valid_channels (i) ) = 1
      end do
   
   end subroutine CHANNEL_SWITCH_ON
   
   ! --------------------------------------------------------------------
   !  every incosistency between channel settings and algorithm mode 
   ! --------------------------------------------------------------------
   subroutine  EXPERT_MODE_CHANNEL_ALGORITHM_CHECK ( sensorname ) 
      character ( len =10) , intent(in) :: sensorname  
      
      integer :: valid_channels ( 42)
      integer :: i
      logical :: not_run_flag
      
      if ( expert_mode < 7 ) return

      chan_on_flag_default = chan_on_flag_default_user_set

      ! - turn off channels not available for this sensor
      
      valid_channels = existing_channels ( sensorname )
      
      do i = 1, 42 
         if ( any ( i == valid_channels )) cycle
         Chan_On_Flag_Default ( i ) = 0
      end do
      
      
       !--- check ACHA mode based on available channels
      if (Acha_Mode == 3 .and. &
         (Chan_On_Flag_Default(32)==sym%NO)) then
            not_run_flag = .true.
         
      endif
      if (Acha_Mode == 4 .and. &
         (Chan_On_Flag_Default(33)==sym%NO)) then
            not_run_flag = .true.
         
      endif
      if (Acha_Mode == 8 .and. &
         (Chan_On_Flag_Default(32)==sym%NO .or. Chan_On_Flag_Default(33)==sym%NO)) then
            not_run_flag = .true.
         
      endif
      if (Acha_Mode == 5 .and. &
         (Chan_On_Flag_Default(29)==sym%NO .or. Chan_On_Flag_Default(32)==sym%NO)) then
            not_run_flag = .true.
        
      endif
      if (Acha_Mode == 6 .and. &
         (Chan_On_Flag_Default(27)==sym%NO .or. Chan_On_Flag_Default(32)==sym%NO)) then
            not_run_flag = .true.
        
      endif
      if (Acha_Mode == 7 .and. &
         (Chan_On_Flag_Default(27)==sym%NO .or. Chan_On_Flag_Default(33)==sym%NO)) then
            not_run_flag = .true.
         
      endif
      if (Acha_Mode == 2 .and. &
         (Chan_On_Flag_Default(27)==sym%NO)) then
            not_run_flag = .true.
         
      endif
      
      
      if ( not_run_flag ) then
         print *, EXE_PROMPT, 'ACHA Mode ', acha_mode,' not possible with selected channels. ACHA and DCOMP  will not run.'
         Acha_Mode = 0
         Dcomp_mode = 0
      end if 

      !--- check based on available channels
      if (Dcomp_Mode == 1 .and. &
         (Chan_On_Flag_Default(1) == sym%NO .or. Chan_On_Flag_Default(6)==sym%NO)) then
         print *, EXE_PROMPT, 'DCOMP Mode 1 not possible with selected channels, DCOMP is now off'
      endif
      
      if (Dcomp_Mode == 2 .and. &
         (Chan_On_Flag_Default(1) == sym%NO .or. Chan_On_Flag_Default(7)==sym%NO)) then
         print *, EXE_PROMPT, 'DCOMP Mode 2 not possible with selected channels, DCOMP is now off'
      endif
      
      if (Dcomp_Mode == 3 .and. &
         (Chan_On_Flag_Default(1) == sym%NO .or. Chan_On_Flag_Default(20)==sym%NO)) then
         print *, EXE_PROMPT, 'DCOMP Mode 3 not possible with selected channels, DCOMP is now off'
      endif
  

   
   end subroutine EXPERT_MODE_CHANNEL_ALGORITHM_CHECK
   
   
  

end module USER_OPTIONS
