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
!       30 Dec 2014: Submitting to trunk
!
!   
!--------------------------------------------------------------------------------------
module USER_OPTIONS

   use PIXEL_COMMON, only: &
        Sensor &
      , Geo &
      , Nav &
      , Image  &
      , ACHA &
      , Aer_Flag &
      , Ash_File_Flag &
      , Ash_Flag &
      , Goes_Stride &
      , Cld_Flag &
      , Cloud_Mask_Aux_Flag &
      , Cloud_Mask_Bayesian_Flag &
      , Sasrab_Flag &
      , Nav_Opt &
      , Nwp_Opt &
      , Ref_Cal_1b &
      , Rtm_File_Flag &
      , Rtm_Opt &
      , Sst_File_Flag &
      , Temporary_Data_Dir &
      , Therm_Cal_1b &
      , Ancil_Data_Dir &
      , Cfsr_Data_Dir &
      , Gdas_Data_Dir &
      , File_List &
      , Geo_File_Flag &
      , Gfs_Data_Dir &
      , Level2_File_Flag &
      , Lrc_Flag &
      , Modis_Clr_Alb_Flag &
      , Ncep_Data_Dir &
      , Obs_File_Flag &
      , Oisst_Data_Dir &
      , Output_Scaled_Reflectances &
      , Process_Undetected_Cloud_Flag &
      , Read_Coast_Mask &
      , Read_Hires_Sfc_Type &
      , Read_Land_Mask &
      , Read_Snow_Mask &
      , Read_Surface_Elevation &
      , Read_Volcano_Mask &
      , Smooth_Nwp_Flag &
      , Snow_Data_Dir &
      , Use_Default &
      , Use_Seebor &
      , Bayesian_Cloud_Mask_Name &
      , Compress_Flag &
      , Nlcomp_Mode &
      , Dcomp_mode &
      , Avhrr_1_flag &
      , Read_Dark_Comp &
      , Globsnow_Data_Dir
      
      
   use CONSTANTS, only: &
      Sym &
      , Exe_Prompt &
      , Nchan_Clavrx
      
   use FILE_UTILITY, only: &
      Get_Lun
 
   use  CLAVRX_MESSAGE_MODULE, only: &
      Mesg &
    , Verb_Lev

   implicit none
   private
   public  :: SETUP_USER_DEFINED_OPTIONS
   public  :: UPDATE_CONFIGURATION

   character(24), parameter, private :: MOD_PROMPT = " USER_OPTIONS_ROUTINES: "
   character ( len = 250 ) :: Data_Base_Path
   integer :: Dcomp_Mode_User_Set
   integer :: Acha_Mode_User_Set
   integer :: Nlcomp_Mode_User_Set
   integer :: Mask_Mode_User_Set
   integer :: Expert_Mode
   
   integer :: Chan_On_Flag_Default_User_Set (Nchan_Clavrx)
   
   ! ---------------------------------------------------------------------------------
   ! Default Algorithm Modes - 
   ! ---------------------------------------------------------------------------------
   integer,parameter:: ACHA_Mode_Default_Avhrr = 3
   integer,parameter:: ACHA_Mode_Default_Avhrr1 = 1
   integer,parameter:: ACHA_Mode_Default_Goes_IL = 6
   integer,parameter:: ACHA_Mode_Default_Goes_MP = 7
   integer,parameter:: ACHA_Mode_Default_Goes_SNDR = 7
   integer,parameter:: ACHA_Mode_Default_COMS = 6
   integer,parameter:: ACHA_Mode_Default_VIIRS = 5
   integer,parameter:: ACHA_Mode_Default_MTSAT = 6
   integer,parameter:: ACHA_Mode_Default_SEVIRI = 8
   integer,parameter:: ACHA_Mode_Default_Modis = 8
   integer,parameter:: ACHA_Mode_Default_Fy2 = 7
   integer,parameter:: ACHA_Mode_Default_AHI = 1
   
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
   subroutine SET_DEFAULT_VALUES
           
      Aer_Flag = sym%YES
      Ash_Flag = sym%NO
      modis_clr_alb_Flag = 1 ! do not use clear-sky MODIS albedo maps
      output_scaled_reflectances = sym%NO !default is to output ref / cossolzen
      
      !--- default solar zenith limits
      Geo%Solzen_Min_Limit= 0 
      Geo%Solzen_Max_Limit= 180.0
      Geo%Satzen_Min_Limit= 0 
      Geo%Satzen_Max_Limit= 85.0
       
      !--- default what can be changed for expert mode
      Cloud_Mask_Bayesian_Flag = 1
      Dcomp_Mode_User_Set = 3
      Acha_Mode_User_Set = 1
      Nlcomp_Mode = 1            
      Level2_File_Flag = 1
      Rtm_File_Flag = 0
      Cld_Flag = 1
      Image%Number_Of_Lines_Per_Segment = 240
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
      Ref_Cal_1b = 0
      Therm_Cal_1b = 0    
      Nav_Opt = 0  
      goes_stride = 1
      Lrc_Flag = 1 
      Smooth_Nwp_Flag = 1  
      Process_Undetected_Cloud_Flag = 0
      Chan_On_Flag_Default_User_Set(1:6) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_User_Set(7:12) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_User_Set(13:18) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_User_Set(19:24) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_User_Set(25:30) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_User_Set(31:36) = [1,1,1,1,1,1]
      Chan_On_Flag_Default_User_Set(37:42) = [1,1,0,0,0,0]
      Chan_On_Flag_Default_User_Set(43:45) = [0,1,0]
      Nav%Lat_Max_Limit = 90.0
      Nav%Lat_Min_Limit = -90.0
      Nav%Lon_Max_Limit = 180.0
      Nav%Lon_Min_Limit = -180.0
      Geo%Satzen_Max_Limit = 85.0
      Geo%Satzen_Min_Limit = 0.0
      Geo%Solzen_Max_Limit = 180.0
      Geo%Solzen_Min_Limit = 0.0
      
   end subroutine SET_DEFAULT_VALUES
   

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
      integer:: GeoNav_Limit_Flag
            
      call MESG ("DEFAULT FILE READ IN",level = 5 )
      call MESG ("Default file to be read in: "//trim(File_Default),level = verb_lev % DEFAULT)

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
      read(unit=Default_Lun,fmt=*) Expert_Mode
         
      if ( Expert_Mode  == 0 )  then
          close(unit=Default_Lun)
          return
      end if
      
      read(unit=Default_Lun,fmt=*) Cloud_Mask_Bayesian_Flag
      read(unit=Default_Lun,fmt=*) Dcomp_Mode_User_Set
      read(unit=Default_Lun,fmt=*) Acha_Mode_User_Set
      read(unit=Default_Lun,fmt=*) Nlcomp_Mode
         
      if ( Expert_Mode <= 1 )  then
          close(unit=Default_Lun)
          return
      end if
      
      read(unit=Default_Lun,fmt=*) Level2_File_Flag
      read(unit=Default_Lun,fmt=*) Rtm_File_Flag
      read(unit=Default_Lun,fmt=*) Cld_Flag
      read(unit=Default_Lun,fmt=*) Image%Number_Of_Lines_Per_Segment
      read(unit=Default_Lun,fmt=*) Sasrab_Flag
      read(unit=Default_Lun,fmt=*) Nwp_Opt
     
      read(unit=Default_Lun,fmt=*) Rtm_Opt
      read(unit=Default_Lun,fmt=*) Nav_Opt
      
      read(unit=Default_Lun,fmt=*) Compress_Flag
      read(unit=Default_Lun,fmt=*) Cloud_Mask_Aux_Flag

      read(unit=Default_Lun,fmt="(a)") Bayesian_Cloud_Mask_Name
      
      if ( Expert_Mode <= 2 )  then
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
         
      if ( Expert_Mode <= 3 ) then
          close(unit=Default_Lun)
          return
      end if
         
      read(unit=Default_Lun,fmt=*) Ref_Cal_1b
      read(unit=Default_Lun,fmt=*) Therm_Cal_1b
            
      if ( Expert_Mode <= 4 ) then
          close(unit=Default_Lun)
          return
      end if
      
      read(unit=Default_Lun,fmt=*) Lrc_Flag
      read(unit=Default_Lun,fmt=*) Smooth_Nwp_Flag
      read(unit=Default_Lun,fmt=*) Process_Undetected_Cloud_Flag
               
      if ( Expert_Mode <= 5 ) then
          close(unit=Default_Lun)
          return
      end if

      ! --- Read lat, lon and sun angle high - low limits
      read(unit=Default_Lun,fmt=*) GeoNav_Limit_Flag

      if (GeoNav_Limit_Flag == 0) then
         Nav%Lat_Max_Limit = 90.0
         Nav%Lat_Min_Limit = -90.0
         Nav%Lon_Max_Limit = 180.0
         Nav%Lon_Min_Limit = -180.0
         Geo%Satzen_Max_Limit = 85.0
         Geo%Satzen_Min_Limit = 0.0
         Geo%Solzen_Max_Limit = 180.0
         Geo%Solzen_Min_Limit = 0.0
      else
         backspace(unit=Default_Lun)
         read(unit=Default_Lun,fmt=*) GeoNav_Limit_Flag, &
                                      Nav%Lat_Max_Limit, Nav%Lat_Min_Limit, &
                                      Nav%Lon_Max_Limit, Nav%Lon_Min_Limit, &
                                      Geo%Satzen_Max_Limit, Geo%Satzen_Min_Limit, &
                                      Geo%Solzen_Max_Limit, Geo%Solzen_Min_Limit
      endif

      !--- constrain values
      Nav%Lat_Max_Limit = min(Nav%Lat_Max_Limit, 90.0)
      Nav%Lat_Min_Limit = max(Nav%Lat_Min_Limit, -90.0)
      Nav%Lon_Max_Limit = min(Nav%Lon_Max_Limit, 180.0)
      Nav%Lon_Min_Limit = max(Nav%Lon_Min_Limit, -180.0)
      Geo%Satzen_Max_Limit = min(Geo%Satzen_Max_Limit, 85.0)
      Geo%Satzen_Min_Limit = max(Geo%Satzen_Min_Limit, 0.0)
      Geo%Solzen_Max_Limit = min(Geo%Solzen_Max_Limit, 180.0)
      Geo%Solzen_Min_Limit = max(Geo%Solzen_Min_Limit, 0.0)


      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(1:6)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(7:12)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(13:18)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(19:24)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(25:30)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(31:36)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(37:42)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default_User_Set(43:45)
             
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

      temp_string = '.'  !--- temporary string to search for in angle commandline

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
        elseif(trim(fargv) == "-no_output_comp") then
          Compress_Flag=0
        elseif(trim(fargv) == "-output_comp_gzip") then
          Compress_Flag=1
        elseif(trim(fargv) == "-output_comp_szip") then
          Compress_Flag=2

        !Change lat max/min for processing
        elseif(trim(fargv) == "-lat_min_limit") then
           call getarg(i+1,junk)
           int_temp = scan(junk,temp_string, back) 
           if(int_temp > 0.0) read(junk,'(f6.3)') Nav%Lat_Min_Limit
           if(int_temp == 0.0) read(junk,'(f6.0)') Nav%Lat_Max_Limit
        elseif(trim(fargv) == "-lat_max_limit") then
          call getarg(i+1,junk)
          int_temp = scan(junk,temp_string, back) 
          if(int_temp > 0.0) read(junk,'(f6.3)') Nav%Lat_Max_Limit
          if(int_temp == 0.0) read(junk,'(f6.0)') Nav%Lat_Max_Limit

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
        if(Temp_Scans_Arg > 1) read(junk,'(i4)') Image%Number_of_Lines_Per_Segment

        !Change solar zenith angle limits
        elseif(trim(fargv) == "-solzen_min_limit") then
          call getarg(i+1,junk)
          int_temp = scan(junk,temp_string, back) 
          if(int_temp > 0.0) read(junk,'(f6.3)') Geo%Solzen_Min_Limit
          if(int_temp == 0.0) read(junk,'(f6.0)') Geo%Solzen_Min_Limit
        
       !Change dcomp mode
        elseif(trim(fargv) == "-dcomp_mode") then
          call getarg(i+1,junk)
          read(junk,'(i1)') Dcomp_Mode_User_Set
          
          
         elseif(trim(fargv) == "-solzen_max_limit") then
            call getarg(i+1,junk)
            int_temp = scan(junk,temp_string, back) 
            if(int_temp > 0.0) read(junk,'(f6.3)') Geo%Solzen_Max_Limit
            if(int_temp == 0.0) read(junk,'(f6.0)') Geo%Solzen_Max_Limit
         elseif (trim(fargv) == "-filelist") then
            call getarg(i+1,file_list)
            file_list=trim(file_list)
         endif
      enddo
  

      !--- default ancillary data directory
      Ancil_Data_dir = trim(Data_Base_Path)
      Gfs_Data_Dir = trim(Data_Base_Path)//'/dynamic/gfs/'
      Ncep_Data_Dir = trim(Data_Base_Path)//'/dynamic/ncep-reanalysis/'
      Cfsr_Data_Dir = trim(Data_Base_Path)//'/dynamic/cfsr/'
      Gdas_Data_Dir = trim(Data_Base_Path)//'/dynamic/gdas/'
      Oisst_data_Dir = trim(Data_Base_Path)//'/dynamic/oisst/'
      Snow_Data_Dir = trim(Data_Base_Path)//'/dynamic/snow/hires/'
      Globsnow_Data_Dir = trim(Data_Base_Path)//'/dynamic/snow/globsnow/'
      
      
      Sensor%Chan_On_Flag_Default = Chan_On_Flag_Default_User_Set

   end subroutine DETERMINE_USER_CONFIG

   !-------------------------------------------------------------------------------
   !--- QC options and modify as needed
   !-------------------------------------------------------------------------------
   subroutine QC_CLAVRXORB_OPTIONS()

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
         call MESG ("GFS data will be used",level = verb_lev % DEFAULT)
      case ( 2 )
         call MESG ( "NCEP Reanalysis data will be used",level = verb_lev % DEFAULT)
      case ( 3 )
         call MESG ( "NCEP Climate Forecast System Reanalysis data will be used",level = verb_lev % DEFAULT)
      case ( 4 )
         call MESG ( "GDAS Reanalysis data will be used",level = verb_lev % DEFAULT)
      case default
         print *,  EXE_PROMPT, "unrecognized value for Nwp_Opt: ", Nwp_Opt
         stop "6-Nwp_Flag"
      
      end select
      
      if (Cloud_Mask_Bayesian_Flag == sym%YES) then
         call MESG  ("Bayesian cloud mask will be generated")
      endif

      if (Ref_Cal_1b == sym%YES) then
         call MESG ("Reflectance Calibration within 1b will be used")
      endif

      if (therm_Cal_1b == sym%YES) then
         call MESG ("Thermal Calibration within 1b will be used")
      endif

      if (nav_Opt == 1) then
         call MESG ("CLEVERNAV geolocation no longer supported, using REPOSNX")
         nav_Opt = 2
       endif

      if (nav_opt == 2) then
         call MESG( "REPOSNX geolocation adjustment done")
      endif

      if (rtm_file_Flag == sym%YES) then
        call MESG( "rtm file will be created")
      endif

      if (Cloud_Mask_Aux_Flag == sym%YES) then
         print *,  EXE_PROMPT, "Cloud mask results will be read in from an aux file"
      endif

      if (Rtm_opt /=1) then
         print *,  EXE_PROMPT, "Only PFAAST RTM implemented, switching to rtm ==1"
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

      print *,"  -lines_per_seg (Imager%Number_Of_Lines_Per_Segment)"
      print *, "  specify the number of lines per segment"
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
   subroutine UPDATE_CONFIGURATION (SensorName)
      character (len=*) , intent(in) :: SensorName
      
      if ( Expert_Mode == 0 ) then
         acha_mode_User_Set =  default_acha_mode ( SensorName )
         dcomp_mode_User_Set = default_dcomp_mode ( SensorName )
      end if

      call CHECK_ALGORITHM_CHOICES(SensorName)
     
      call CHANNEL_SWITCH_ON (SensorName)

      if ( Expert_Mode < 2 .or. trim(Bayesian_Cloud_Mask_Name) == 'default') then
         Bayesian_Cloud_Mask_Name = default_nb_mask_classifier_file ( SensorName )
      end if
      
      call EXPERT_MODE_CHANNEL_ALGORITHM_CHECK ( SensorName ) 


   end subroutine UPDATE_CONFIGURATION
   
   !----------------------------------------------------------------------
   !  returns default acha mode 
   !----------------------------------------------------------------------
   integer function default_acha_mode ( SensorName )
      character ( len =*) , intent(in) :: SensorName
      
      select case ( trim(SensorName))
      
      case ( 'AVHRR-2')        
         default_acha_mode  = ACHA_Mode_Default_Avhrr
      case ( 'AVHRR-3')        
         default_acha_mode  = ACHA_Mode_Default_Avhrr
      case ( 'AVHRR-1')   
         default_acha_mode = ACHA_Mode_Default_Avhrr1   
      case ( 'GOES-MP-IMAGER')      
         default_acha_mode  = ACHA_Mode_Default_Goes_MP
      case ( 'GOES-IL-IMAGER')      
         default_acha_mode  = ACHA_Mode_Default_Goes_IL   
      case ( 'GOES-IP-SOUNDER')
         default_acha_mode  = ACHA_Mode_Default_Goes_SNDR 
      case ( 'MTSAT-IMAGER')
          default_acha_mode  = ACHA_Mode_Default_MTSAT
      case ('SEVIRI')
         default_acha_mode  = ACHA_Mode_Default_SEVIRI
      case ('FY2-IMAGER')
         default_acha_mode  =  ACHA_Mode_Default_FY2 
      case ('VIIRS')
         default_acha_mode  =  ACHA_Mode_Default_VIIRS 
      case ('VIIRS-IFF')      
          default_acha_mode  = ACHA_Mode_Default_VIIRS      
      case ('AQUA-IFF')
          default_acha_mode  = ACHA_Mode_Default_Modis
      case ('AVHRR-IFF')      
         default_acha_mode  = ACHA_Mode_Default_Avhrr
      case ('COMS-IMAGER')
         default_acha_mode  = ACHA_Mode_Default_COMS
      case ('MODIS')
          default_acha_mode  = ACHA_Mode_Default_Modis 
      case ('MODIS-MAC')
          default_acha_mode  = ACHA_Mode_Default_Modis 
      case ('MODIS-CSPP')
          default_acha_mode  = ACHA_Mode_Default_Modis 
      case ('AHI')
          default_acha_mode  = ACHA_Mode_Default_AHI   
      case default 
         print*,'sensor ',SensorName, ' is not set in user_options.f90: check channels settings Inform andi.walther@ssec.wisc.edu'   
      end select
          
   end function default_acha_mode
   
   !-----------------------------------------------------------------
   !   returns default dcomp mode
   ! -----------------------------------------------------------------
   integer function default_dcomp_mode ( SensorName )
      character (len=*) , intent(in) :: SensorName
   
      default_dcomp_mode = 3
      
      if (Sensor%WMO_Id == 208) Dcomp_Mode = 1 !NOAA-17
      if (Sensor%WMO_Id == 3) Dcomp_Mode = 1   !METOP-A
      if (Sensor%WMO_Id == 4) Dcomp_Mode = 1   !METOP-B
      if (Sensor%WMO_Id == 5) Dcomp_Mode = 1   !METOP-C


!---- AKH - What about NOAA-16?
     
   end function default_dcomp_mode

!-----------------------------------------------------------------
!   returns default classifier name
!-----------------------------------------------------------------
   function default_nb_mask_classifier_file (SensorName) result (filename)
      character ( len = *) , intent(in) :: SensorName
      character ( len = 355 ) :: filename

      select case ( trim(SensorName))
      
      case ( 'AVHRR')        
         filename  = 'avhrr_default_nb_cloud_mask_lut.nc'
      case ( 'AVHRR-1')   
         filename  = 'avhrr_default_nb_cloud_mask_lut.nc'   
      case ( 'AVHRR-2')   
         filename  = 'avhrr_default_nb_cloud_mask_lut.nc'   
      case ( 'AVHRR-3')   
         filename  = 'avhrr_default_nb_cloud_mask_lut.nc'   
      case ( 'GOES-MP-IMAGER')      
         filename  = 'goesmp_default_nb_cloud_mask_lut.nc'
      case ( 'GOES-IL-IMAGER')      
         filename  = 'goesil_default_nb_cloud_mask_lut.nc'   
      case ( 'GOES-IP-SOUNDER')
         filename  = 'goesmp_default_nb_cloud_mask_lut.nc' 
      case ( 'MTSAT-IMAGER')
          filename  = 'mtsat_default_nb_cloud_mask_lut.nc'
      case ('SEVIRI')
         filename  = 'seviri_default_nb_cloud_mask_lut.nc'
      case ('FY2-IMAGER')
         filename  = 'avhrr_default_nb_cloud_mask_lut.nc' 
      case ('VIIRS')
         filename  = 'viirs_default_nb_cloud_mask_lut.nc' 
      case ('VIIRS-IFF')      
          filename  = 'viirs_default_nb_cloud_mask_lut.nc'
      case ('AQUA-IFF')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case ('AVHRR-IFF')      
        filename  = 'avhrr_default_nb_cloud_mask_lut.nc'
      case ('COMS-IMAGER')
         filename  = 'coms_default_nb_cloud_mask_lut.nc'
      case ('MODIS')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case ('MODIS-CSPP')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case ('MODIS-MAC')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'
      case ('AHI')
          filename  = 'modis_default_nb_cloud_mask_lut.nc'   
          
      case default 
         print*,'sensor ',SensorName, ' is not set in user_options.f90:  Inform andi.walther@ssec.wisc.edu'  
         stop 
      end select


     end function default_nb_mask_classifier_file

   !----------------------------------------------------------------------------
   !  check if algo mode set by user is possible
   !----------------------------------------------------------------------------
   subroutine CHECK_ALGORITHM_CHOICES(SensorName)
      character (len=*) , intent(in) :: SensorName
      
      integer :: possible_acha_modes ( 8 )
      integer :: possible_dcomp_modes ( 3 )
 
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
      
      acha % mode = acha_mode_User_Set     
      dcomp_mode = dcomp_mode_User_Set
       
      possible_acha_modes = 0 
      possible_dcomp_modes = 0
         
      select case ( trim ( SensorName))
      
      case ( 'AVHRR-3')  
         possible_acha_modes(1:2)   = [1, 3]
         possible_dcomp_modes(1)    =  3 
      case ( 'AVHRR-2')  
         possible_acha_modes(1:2)   = [1, 3]
         possible_dcomp_modes(1)    =  3 
      case ( 'AVHRR-1')   
         possible_acha_modes(1)     =  1
         possible_dcomp_modes(1)    =  3   
      case ( 'GOES-MP-IMAGER')      
         possible_acha_modes(1:3)   =  [1, 3, 7]
         possible_dcomp_modes(1)    =  3
      case ( 'GOES-IL-IMAGER')      
         possible_acha_modes(1:3)   =  [1,  3, 6] 
         possible_dcomp_modes(1)    =  3  
      case ( 'GOES-IP-SOUNDER')
         possible_acha_modes(1:8)   =  [1, 2, 3, 4, 5, 6, 7, 8]  
         possible_dcomp_modes(1)    =  3
      case ( 'MTSAT-IMAGER')
         possible_acha_modes(1:4)   =  [ 1, 2, 3 , 6 ]
         possible_dcomp_modes(1)    =  3
      case ('SEVIRI')
         possible_acha_modes(1:8)   =  [1, 2, 3, 4, 5, 6, 7, 8]
         possible_dcomp_modes(1:2)  =  [1, 3]
      case ('FY2-IMAGER')
         possible_acha_modes(1:2)   =  [1 , 2 ] 
         possible_dcomp_modes(1:1)  =  [3]
      case ('VIIRS')
         possible_acha_modes(1:3)   =  [1, 3, 5] 
         possible_dcomp_modes(1:3)  =  [1, 2, 3]
         nlcomp_mode_User_Set       =  1  
      case ('VIIRS-IFF')      
         possible_acha_modes(1:3)   =  [1, 3, 5]
         possible_dcomp_modes(1:3)  =  [1, 2, 3]
      case ('AQUA-IFF')
         possible_acha_modes(1:8)   =  [1, 2, 3, 4, 5, 6, 7, 8]
         possible_dcomp_modes(1:3)  =  [1, 2, 3]
      case ('AVHRR-IFF')
         possible_acha_modes(1:4)   =  [1, 3, 4, 8]
         possible_dcomp_modes(1)    =  3
      case ('COMS-IMAGER')
         possible_acha_modes(1:3)   =  [1, 3, 6]
         possible_dcomp_modes(1:1)  =  [3]
      case ('MODIS')
         possible_acha_modes(1:8)   =  [1, 2, 3, 4, 5, 6, 7, 8] 
         possible_dcomp_modes(1:3)  =  [1, 2, 3]
      case ('MODIS-CSPP')
         possible_acha_modes(1:8)   =  [1, 2, 3, 4, 5, 6, 7, 8]
         possible_dcomp_modes(1:3)  =  [1, 2, 3]
      case ('MODIS-MAC')
         possible_acha_modes(1:8)   =  [1, 2, 3, 4, 5, 6, 7, 8]
         possible_dcomp_modes(1:3)  =  [1, 2, 3]
      case ( 'AHI')
         possible_acha_modes(1:8)   =  [1, 2, 3, 4, 5, 6, 7, 8]
         possible_dcomp_modes(1:3)  =  [1, 2, 3]
      case default 
         print*,'sensor ',SensorName, ' is not set in check channels user_options settings Inform andi.walther@ssec.wisc.edu'   
      end select
      
      if ( acha_mode_user_set /= 0 .and. .not. ANY ( acha_mode_User_Set == possible_acha_modes ) ) then
         acha % mode = default_acha_mode ( SensorName )
         
         print*, 'User set ACHA mode not possible for '//trim(SensorName)//' switched to default ', default_acha_mode ( SensorName )
      end if
      
      if ( dcomp_mode_user_set /= 0 .and. .not. ANY ( dcomp_mode_User_Set == possible_dcomp_modes ) ) then
         dcomp_mode = default_dcomp_mode ( SensorName )
         print*, 'User set DCOMP mode not possible for '//trim(SensorName)//' switched to default ', default_dcomp_mode ( SensorName )
      end if

   end subroutine CHECK_ALGORITHM_CHOICES
   
   ! ----------------------------------------------------------------------
   !    returns all available sensors for this sensors
   ! ----------------------------------------------------------------------
   function Existing_Channels  (SensorName)  result( Valid_Channels )
      character (len = *) , intent(in) :: SensorName
      
      integer , target :: Valid_Channels (Nchan_Clavrx) 
      integer :: i 
      
      Valid_Channels = -99
      select case ( trim(SensorName))
        
      case ( 'AVHRR-1')
         Valid_Channels (1:4) = [1,2,20,31]
      case ( 'AVHRR-2')
         Valid_Channels (1:5) = [1,2,20,31,32]
      case ( 'AVHRR-3')
         Valid_Channels (1:6) = [1,2,6,20,31,32]
      case ( 'GOES-IL-IMAGER')      
         Valid_Channels (1:5) = [1,20,27,31,32]
      case ( 'GOES-MP-IMAGER')      
         Valid_Channels (1:5) = [1,20,27,31,33]   
      case ( 'GOES-IP-SOUNDER')
         Valid_Channels (1:18) = [1,20,21,23,24,25,30,31,32,33,34,35,36,37,38,39,40,41]      
      case ( 'MTSAT-IMAGER')
         Valid_Channels (1:5) = [1,20,27,31,32]  
      case ('SEVIRI')
         Valid_Channels (1:11) = [1,2,6,20,27,28,29,30,31,32,33]
      case ('FY2-IMAGER')
         Valid_Channels (1:5) = [1,20,27,31,32]    
      case ('VIIRS')
         Valid_Channels (1:22) = [1,2,3,4,5,6,7,8,9,15,20,22,26,29,31,32,39,40,41,42,43,44]
      case ('VIIRS-IFF')
         Valid_Channels (1:21) = [1,2,3,4,5,6,7,8,9,15,20,22,26,29,31,32,33,34,35,36,45]
      case ('AQUA-IFF')
         Valid_Channels (1:36) = [(i,i=1,36,1)]
      case ('AVHRR-IFF')
         Valid_Channels (1:20) = [1,2,6,20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36,45]
      case ('COMS-IMAGER')
         Valid_Channels (1:5) = [1,20,27,31,32]
      case ('MODIS')
         Valid_Channels(1:36) = [(i,i=1,36,1)]  
      case ('MODIS-MAC')
         Valid_Channels(1:36) = [(i,i=1,36,1)]  
      case ('MODIS-CSPP')
         Valid_Channels(1:36) = [(i,i=1,36,1)] 
      case ('AHI')
         Valid_Channels(1:16) = [1,2,3,4,6,7,20,27,28,29,30,31,32,33,37,38]    
      case default 
         print*,'sensor ',SensorName, ' is not set in check channels settings Inform andi.walther@ssec.wisc.edu'   
      end select

   end function Existing_Channels
   
   !----------------------------------------------------------------------
   !   Channel settings
   !     will not be done for full-experts  ( expert mode 7 and higher)
   !
   !----------------------------------------------------------------------
   subroutine CHANNEL_SWITCH_ON (SensorName)
      character (len=*) , intent(in) :: SensorName
      integer :: Valid_Channels (Nchan_Clavrx)
      integer :: i
 
      ! expert can decide themselves
      if (Expert_Mode > 6 ) return
      
      Valid_Channels = Existing_Channels ( SensorName )
          
      Sensor%Chan_On_Flag_Default =  0

      do i = 1, Nchan_Clavrx
         if (Valid_Channels (i) < 0 ) cycle
         Sensor%Chan_On_Flag_Default (Valid_Channels (i) ) = 1
      end do
   
   end subroutine CHANNEL_SWITCH_ON
   
   ! --------------------------------------------------------------------
   !  every incosistency between channel settings and algorithm mode 
   ! --------------------------------------------------------------------
   subroutine  EXPERT_MODE_CHANNEL_ALGORITHM_CHECK ( SensorName ) 
      character (len=*) , intent(in) :: SensorName  
      
      integer :: Valid_Channels (Nchan_Clavrx)
      integer :: i
      logical :: Not_Run_Flag
      
      if ( Expert_Mode < 6 ) return

      Sensor%Chan_On_Flag_Default = Chan_On_Flag_Default_User_Set

      ! - turn off channels not available for this sensor
      

      Valid_Channels = Existing_Channels ( SensorName )

      do i = 1, Nchan_Clavrx
         if ( any ( i == Valid_Channels )) cycle
         Sensor%Chan_On_Flag_Default ( i ) = 0
      end do

      !--- check ACHA mode based on available channels
      Not_Run_Flag = .false.
      if (ACHA%Mode > 0 .and. &
         (Sensor%Chan_On_Flag_Default(31)==sym%NO)) then
            Not_Run_Flag = .true.
      endif
      if (ACHA%Mode == 2 .and. &
         (Sensor%Chan_On_Flag_Default(27)==sym%NO)) then
            Not_Run_Flag = .true.
      endif
      if (ACHA%Mode == 3 .and. &
         (Sensor%Chan_On_Flag_Default(32)==sym%NO)) then
            Not_Run_Flag = .true.
      endif
      if (ACHA%Mode == 4 .and. &
         (Sensor%Chan_On_Flag_Default(33)==sym%NO)) then
            Not_Run_Flag = .true.
      endif
      if (ACHA%Mode == 5 .and. &
         (Sensor%Chan_On_Flag_Default(29)==sym%NO .or. Sensor%Chan_On_Flag_Default(32)==sym%NO)) then
            Not_Run_Flag = .true.
      endif
      if (ACHA%Mode == 6 .and. &
         (Sensor%Chan_On_Flag_Default(27)==sym%NO .or. Sensor%Chan_On_Flag_Default(32)==sym%NO)) then
            Not_Run_Flag = .true.
      endif
      if (ACHA%Mode == 7 .and. &
         (Sensor%Chan_On_Flag_Default(27)==sym%NO .or. Sensor%Chan_On_Flag_Default(33)==sym%NO)) then
            Not_Run_Flag = .true.
      endif

      if (ACHA%Mode == 8 .and. &
         (Sensor%Chan_On_Flag_Default(32)==sym%NO .or. Sensor%Chan_On_Flag_Default(33)==sym%NO)) then
            Not_Run_Flag = .true.
      endif
         
      if ( Not_Run_Flag ) then
         print *, EXE_PROMPT, 'ACHA Mode ', ACHA%Mode,' not possible with selected channels. ACHA and DCOMP  will not run.'
         ACHA%Mode = 0
         Dcomp_Mode = 0
      end if 

      !--- check based on available channels
      if (Dcomp_Mode == 1 .and. &
         (Sensor%Chan_On_Flag_Default(1) == sym%NO .or. Sensor%Chan_On_Flag_Default(6)==sym%NO)) then
         print *, EXE_PROMPT, 'DCOMP Mode 1 not possible with selected channels, DCOMP is now off'
      endif
      
      if (Dcomp_Mode == 2 .and. &
         (Sensor%Chan_On_Flag_Default(1) == sym%NO .or. Sensor%Chan_On_Flag_Default(7)==sym%NO)) then
         print *, EXE_PROMPT, 'DCOMP Mode 2 not possible with selected channels, DCOMP is now off'
      endif
      
      if (Dcomp_Mode == 3 .and. &
         (Sensor%Chan_On_Flag_Default(1) == sym%NO .or. Sensor%Chan_On_Flag_Default(20)==sym%NO)) then
         print *, EXE_PROMPT, 'DCOMP Mode 3 not possible with selected channels, DCOMP is now off'
      endif
  
   end subroutine EXPERT_MODE_CHANNEL_ALGORITHM_CHECK
   
end module USER_OPTIONS
