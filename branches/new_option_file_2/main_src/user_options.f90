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
!
!  HISTORY:
!        16 Dec 2014: Switch to new option file  (AW)
!                       code cleaning
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
   public :: SETUP_USER_DEFINED_OPTIONS
   public :: CHECK_ALGORITHM_CHOICES
   public :: CHECK_CHANNEL_SETTINGS

   character(24), parameter, private :: MOD_PROMPT = " USER_OPTIONS_ROUTINES: "
   character ( len = 50 ) :: data_base_path
   integer :: Dcomp_Mode_User_Set
   integer :: Acha_Mode_User_set
   integer :: expert_mode
   
   
   
   !---------------------------------------------------------------------------------
   ! Default Algorithm Modes - 
   !---------------------------------------------------------------------------------
   integer,parameter:: ACHA_Mode_Default_Avhrr = 3
   integer,parameter:: ACHA_Mode_Default_Avhrr1 = 1
   integer,parameter:: ACHA_Mode_Default_Goes_IL = 6
   integer,parameter:: ACHA_Mode_Default_Goes_MP = 7
   integer,parameter:: ACHA_Mode_Default_VIIRS = 5
   integer,parameter:: ACHA_Mode_Default_MTSAT = 6
   integer,parameter:: ACHA_Mode_Default_SEVIRI = 8
   integer,parameter:: ACHA_Mode_Default_Modis = 8
   
contains

   !------------------------------------------------------------------
   !
   !------------------------------------------------------------------
   subroutine SETUP_USER_DEFINED_OPTIONS()
      call SET_DEFAULT_VALUES
      call READ_CLAVRXORB_OPTIONS()
      call QC_CLAVRXORB_OPTIONS()

   end subroutine SETUP_USER_DEFINED_OPTIONS
   
   !
   !
   !
   subroutine set_default_values
   
   
            !--- set yes/no options to no as default
 
      
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
      bayesian_cloud_mask_name = 'viirs_default_bayes_mask.txt'
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
      Chan_On_Flag_Default(1:6) = [1,1,1,1,1,1]
      Chan_On_Flag_Default(7:12) = [1,1,1,1,1,1]
      Chan_On_Flag_Default(13:18) = [1,1,1,1,1,1]
      Chan_On_Flag_Default(19:24) = [1,1,1,1,1,1]
      Chan_On_Flag_Default(25:30) = [1,1,1,1,1,1]
      Chan_On_Flag_Default(31:36) = [1,1,1,1,1,1]
      Chan_On_Flag_Default(37:42) = [0,0,0,0,0,1]
   
   
   end subroutine
   
   
  
   !----------------------------------------------------------------------------
   !
   !----------------------------------------------------------------------------
   subroutine CHECK_ALGORITHM_CHOICES()

      character ( len = 1 ) :: string_1
 
 
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
      
      if (Avhrr_Flag == sym%YES) then 
         if (Avhrr_1_Flag == sym%NO) then 
             if (Acha_Mode /=1 .or. Acha_Mode /= 3) then
                     print *, EXE_PROMPT, &
                     "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT, &
                     "Changing to default for AVHRR ", Acha_Mode_Default_Avhrr
                     Acha_Mode = Acha_Mode_Default_Avhrr 
             endif
         ELSE 
             if (Acha_Mode /= 1)  then 
                     print *, EXE_PROMPT, &
                     "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT, &
                     "Changing to default for AVHRR/1 ",Acha_Mode_Default_Avhrr1
                     Acha_Mode = Acha_Mode_Default_Avhrr1 
             endif
         endif
         if ( expert_mode <= 1 ) bayesian_cloud_mask_name = 'avhrr_default_bayes_mask.txt'
      endif

      write (string_1,'(I1)') acha_mode
      call mesg ("Acha Mode = "//string_1)

      if (Goes_Flag == sym%YES) then 

          if (Goes_Mop_Flag == sym%NO) then 

             if (Acha_Mode == 4 .or. Acha_Mode == 8 .or. &
                 Acha_Mode == 5 .or. Acha_Mode == 7)  then 
                     print *, EXE_PROMPT, &
                     "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT, &
                     "Changing to default for GOES-IL ",Acha_Mode_Default_Goes_IL
                     Acha_Mode = Acha_Mode_Default_Goes_IL
             endif

          ELSE 

             if (Acha_Mode == 3 .or. Acha_Mode == 8 .or.  &
                 Acha_Mode == 5 .or. Acha_Mode == 6)  then 
                     print *, EXE_PROMPT, &
                     "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT, &
                     "Changing to default for GOES-MP"
                     Acha_Mode = Acha_Mode_Default_Goes_MP
                     
               write (string_1,'(I1)') acha_mode
               call mesg ("Acha Mode = "//string_1 )

             endif

          endif
           if ( expert_mode <= 1 ) bayesian_cloud_mask_name = 'goes_default_bayes_mask.txt' 
      endif

      if (Viirs_Flag == sym%YES) then 

             if (Acha_Mode == 4 .or. Acha_Mode == 8 .or.  &
                 Acha_Mode == 6 .or. Acha_Mode == 7 .or. Acha_Mode == 2)  then 
                     print *, EXE_PROMPT, "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT,  "Changing to default for VIIRS"
                     Acha_Mode = Acha_Mode_Default_VIIRS
             endif
         if ( expert_mode <= 1 ) bayesian_cloud_mask_name = 'viirs_default_bayes_mask.txt'
      endif

      if (Mtsat_Flag == sym%YES) then

             if (Acha_Mode == 4 .or. Acha_Mode == 8 .or. &
                 Acha_Mode == 5 .or. Acha_Mode == 7)  then
                     print *, EXE_PROMPT, "Acha_Mode incompatible with satellite observations"
                     print *, EXE_PROMPT,  "Changing to default for MTSAT"
                     Acha_Mode = Acha_Mode_Default_MTSAT
             endif
         if ( expert_mode <= 1 ) bayesian_cloud_mask_name = 'mtsat_default_bayes_mask.txt'       

      endif

      !--- check ACHA mode based on available channels
      if (Acha_Mode == 3 .and. &
         (Chan_On_Flag_Default(32)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 3 not possible with selected channels, ACHA Set to Mode 0'
         Acha_Mode = 0
      endif
      if (Acha_Mode == 4 .and. &
         (Chan_On_Flag_Default(33)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 4 not possible with selected channels. ACHA will not run.'
         Acha_Mode = 0
      endif
      if (Acha_Mode == 8 .and. &
         (Chan_On_Flag_Default(32)==sym%NO .or. Chan_On_Flag_Default(33)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 8 not possible with selected channels. ACHA will not run.'
         Acha_Mode = 0
      endif
      if (Acha_Mode == 5 .and. &
         (Chan_On_Flag_Default(29)==sym%NO .or. Chan_On_Flag_Default(32)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 5 not possible with selected channels. Acha will not run.'
         Acha_Mode = 0
      endif
      if (Acha_Mode == 6 .and. &
         (Chan_On_Flag_Default(27)==sym%NO .or. Chan_On_Flag_Default(32)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 6 not possible with selected channels. ACHA will not run.'
         Acha_Mode = 0
      endif
      if (Acha_Mode == 7 .and. &
         (Chan_On_Flag_Default(27)==sym%NO .or. Chan_On_Flag_Default(33)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 7 not possible with selected channels. ACHA will not run.'
         Acha_Mode = 0
      endif
      if (Acha_Mode == 2 .and. &
         (Chan_On_Flag_Default(27)==sym%NO)) then
         print *, EXE_PROMPT, 'ACHA Mode 2 not possible with selected channels. ACHA will not run.'
         Acha_Mode = 0
      endif


      !-------------------------------------------------------------------------------------
      !--- DCOMP Mode Check
      !-------------------------------------------------------------------------------------

      !- dcomp mode 9 is Andys test code
      write (string_1,'(I1)') dcomp_mode_user_set
      call mesg ('dcomp user set:   ===================>  '//string_1 ,level = 5 )
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
	      call mesg ( 'dcomp mode switched due to sensor  setting  ',color=91, level = 0 )
      endif 
     
      write (string_1,'(i1)') dcomp_mode
      call mesg ('dcomp mode used for this file: '//string_1)


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
    

   end subroutine CHECK_ALGORITHM_CHOICES
   
   
   !----------------------------------------------------------------------
   !   Channel settings
   !   
   !
   !----------------------------------------------------------------------
   subroutine CHECK_CHANNEL_SETTINGS(sensorname)
      character ( len =10) , intent(in) :: sensorname
      integer :: valid_channels ( 42)
      integer :: i
      
      ! expert can decide themselves
      if (expert_mode > 6 ) return
      
      valid_channels = -99
      
       Chan_On_Flag_Default =  0
      
    
      select case ( trim(sensorname))
      
      
      case ( 'AVHRR')
         valid_channels (1:6) = [1,2,6,20,31,32]
      case ( 'GOES')      
         valid_channels (1:6) = [1,20,27, 31, 32, 33]
         !if (goes_mop_flag == 1) valid_channels(5) = 33
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
         print*,'sensor ',sensorname, ' is not set in chaeck channels settings Inform andi.walther@ssec.wisc.edu'   
      end select
      
      do i = 1, 42 
         if (valid_channels (i) < 0 ) cycle
         Chan_On_Flag_Default (valid_channels (i)) = 1
      end do
   
   end subroutine CHECK_CHANNEL_SETTINGS
   
   
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
      
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(1:6)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(7:12)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(13:18)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(19:24)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(25:30)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(31:36)
      read(unit=Default_Lun,fmt=*) Chan_On_Flag_Default(37:42)
             
     
      close(unit=Default_Lun)
    

   end subroutine READ_CLAVRXORB_DEFAULT_OPTIONS
 
   subroutine READ_CLAVRXORB_OPTIONS()
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
            print*,'$Header$'
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
         call READ_CLAVRXORB_DEFAULT_OPTIONS(Default_Temp)
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

        !Read/not read land mask
        elseif(trim(fargv) == "-read_land_mask") then
          read_land_mask = sym%YES
        elseif(trim(fargv) == "-no_land_mask") then
          read_land_mask = sym%NO

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
  
  
  
   
  !--- default ancillary data directory
  ancil_data_dir = trim(data_base_path)//'/clavrx_ancil_data/'
  gfs_data_dir = trim(data_base_path)//'gfs/'
  ncep_data_dir = trim(data_base_path)//'/clavrx_ancil_data/ncep-reanalysis/'
  cfsr_data_dir = trim(data_base_path)//'/cfsr/'
  oisst_data_dir = trim(data_base_path)//'/clavrx_ancil_data/oisst/'
  snow_data_dir = './data/snow/'

 end subroutine READ_CLAVRXORB_OPTIONS

!-------------------------------------------------------------------------------
!--- QC options and modify as needed
!-------------------------------------------------------------------------------
 subroutine QC_CLAVRXORB_OPTIONS()

    integer:: erstat

    !---- Since the NWP controls everything, we first check if an NWP is being used
    !---- before anything else is checked.  If no nwp, we stop processing

    if ((Nwp_Opt < 0) .or. (Nwp_Opt > 4)) then
       print *,  EXE_PROMPT, "unrecognized value for Nwp_Opt: ", Nwp_Opt
       stop "6-Nwp_Flag"
    endif
    if (Nwp_Opt == 1) then
       call mesg ("GFS data will be used",level = verb_lev % DEFAULT)
    else if (Nwp_Opt == 2) then
       call mesg ( "NCEP Reanalysis data will be used",level = verb_lev % DEFAULT)
    else if (Nwp_Opt == 3) then
       call mesg ( "NCEP Climate Forecast System Reanalysis data will be used",level = verb_lev % DEFAULT)
    else if (Nwp_Opt == 4) then
       call mesg ( "GDAS Reanalysis data will be used",level = verb_lev % DEFAULT)
    endif
    
    if (Nwp_Opt == 0) then
       print *,  EXE_PROMPT, "No choice made for NWP data, will not run algoritms or orbital level3 files"
       Cld_Flag = sym%NO
       
       Sasrab_Flag = sym%NO
      
       
      
       Rtm_File_Flag = sym%NO
       Cloud_Mask_Bayesian_Flag = sym%NO
       Cloud_Mask_Aux_Flag = sym%NO ! this is to determine if the lut's are being read in
    endif

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

    if (cld_Flag == sym%NO) then
        print *, EXE_PROMPT, "Cloud products will not be created"
    endif

 

   

    if (rtm_file_Flag == sym%YES) then
        call mesg( "rtm file will be created")
    endif

   

    if (Cloud_Mask_Aux_Flag == sym%YES) then
       print *,  EXE_PROMPT, "Cloud mask results will be read in from an aux file"
    endif

 
    if (Rtm_opt /=1) then
       print *,  EXE_PROMPT, "Only PFAST RTM implemented, stopping"
       stop
    endif


    call mesg ("Temporary Files will be written to "//trim(Temporary_Data_Dir),level = verb_lev % VERBOSE )

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
 
  
  print *,"  -l1bnav"
  print *, "  Use the navigation data from level 1b file"
  print *," "

  print *,"  -clevernav"
  print *, "  Use the navigation data from a clevernav"
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
  print *,"  -cloud_mask_calc"
  print *, "   recalculate cloud mask."
  print *," "
  
  print *,"  -cld_Flag"
  print *, "   Run cloud algorithms. "
  print *," "
  print *,"  -no_cld_Flag"
  print *, "   Don't run cloud algorithms. "
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
