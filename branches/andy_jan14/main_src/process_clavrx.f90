! $Id$
  program PROCESS_CLAVRX
!-----------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: process_clavrx.f90 (src)
!       PROCESS_CLAVRX (program)
!       CLAVRXORB (executable)
!
! PURPOSE:
!
! DESCRIPTION: This code serves as the NESDIS operational cloud 
!      processing system (CLAVR-x). This code also serves as the cloud 
!      climate data generation system (PATMOS-x)
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
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
! This copyright pertains to all routines in the CLAVR-x system unless stated
!
! REVISION HISTORY:
! Version 5.0 - 2010 GEWEX submission
! Version 5.2 - Code Delivered to NCDC 
! Version 5.2.1 - MODIS Capability 
! Version 5.2.2 - GOES imager Capability, GlobSnow and CFSR capability
! Version 5.2.3 - Istvan Laszlo's Insolation (SASRAB) added
! Version 5.2.4 - Volcanic Ash added
! Version 5.2.5 - VIIRS M-band Added
! Version 5.2.6 - VIIRS DNB Support
! Version 5.2.7 - MTSAT, COMS Support
! Version 5.2.8 - GOES Sounder Support
! Version 5.3.0 - NCDC Delivery
!
!
! Basic Running Instruction
!  The input to this code is controlled through three mechanisms
!    1. command-line options  (type clavrxorb --help to see documentation)
!    2. a FILELIST - a list of level-1b files and directories (default name is
!      clavrxorb_File_list)
!    3. a OPTIONSLIST - a list of processing options (default is
!                     clavrxorb_Default_Options)
!
! Overview of capabilities.
!    CLAVRXORB can
!       - use AVHRR level-1b calibration or apply new calibration routines
!       - use AVHRR level-1b geolocation or apply new geolocation routines
!       - process NESDIS or AAPP AVHRR Level1b 
!       - process MYD021KM or MYD02SSH MODIS Level1b files
!       - process band-separated AREA files from GOES Imager/Sounder,  
!         SEVIRI, MTSAT-1R, MTSAT-2, COMS, FY-2
!       - generate pixel level cloud, aerosol and surface products
!       - write to a series of pixel-level hdf files
!       - write a level-3 file (gridded data for each orbit - AVHRR only)
!
! In general, CLAVRXORB uses global data arrays and structures to pass data
!
! Note, comments the begin with "Marker" refer to flowchart delivered to NCDC
!
! Web-page:  http://cimss.ssec.wisc.edu/clavr or 
!            http://cimss.ssec.wisc.edu/patmosx
!
! Channels 1 - 36 refer to MODIS or their analogs on other sensors
! Channels 37-42 are defined only for VIIRS
! Channel 37 - VIIRS I1 - 0.64 micron
! Channel 38 - VIIRS I2 - 0.865 micron
! Channel 39 - VIIRS I3 - 1.61 micron
! Channel 40 - VIIRS I4 - 3.74 micron
! Channel 41 - VIIRS I5 - 11.45 micron
! Channel 42 - VIIRS DNB - 0.7 micron
!
!-------------------------------------------------------------------------
 
!*****************************************************************************
! Marker: ACCESS MODULES 
!******************************************************************************
   use CONSTANTS
   use HDF
   use PIXEL_COMMON
   use PIXEL_ROUTINES
   use LEVEL2_ROUTINES
   use SURFACE_PROPERTIES
   use CLOUD_HEIGHT_ROUTINES
   use ACHA_CLAVRX_BRIDGE
   use CLOUD_BASE_CLAVRX_BRIDGE
   use DCOMP_CLAVRX_BRIDGE_MOD
   use NLCOMP_BRIDGE_MOD
   use AEROSOL_PROPERTIES
   use HDF_PARAMS
   
   use GLOBSNOW_READ_ROUTINES
   use GFS
   use NCEP_REANALYSIS
   use DCOMP_DERIVED_PRODUCTS_MODULE
   
   use RT_UTILITIES, only: &
        rtm_nvzen &
      , setup_pfaast &
      , setup_solar_rtm &
      , map_nwp_rtm  &
      , create_temp_nwp_vectors  &
      , destroy_temp_nwp_vectors &
      , get_pixel_nwp_rtm &
      , allocate_rtm &
      , deallocate_rtm &
      , deallocate_rtm_vars &
      , deallocate_rtm_cell  
      
   use RTM_COMMON,only: &
      nlevels_rtm
   
   use NWP_COMMON
   use SCALING_PARAMETERS
   use PLANCK
   use AVHRR_REPOSITION_ROUTINES
   use NB_CLOUD_MASK_CLAVRX_BRIDGE, only: &
       NB_CLOUD_MASK_BRIDGE
   use MODIS_MODULE
   use IFF_CLAVRX_BRIDGE
   use GOES_MODULE
   use LASZLO_INSOLATION
   use SEVIRI_MODULE
   use MTSAT_MODULE
   use COMS_MODULE
   use FY2_MODULE
   use SENSOR_MODULE
   
   use USER_OPTIONS,only:&
      SETUP_USER_DEFINED_OPTIONS &
      , update_configuration
      
   use CLAVRX_MESSAGE_MODULE, only: &
      mesg &
      , verb_lev
   use CLOUD_TYPE_BRIDGE_MODULE
   use SIMPLE_COD
 
   use dnb_retrievals_mod, only: &
      COMPUTE_LUNAR_REFLECTANCE
      
   use cr_config_mod 
   use date_tools_mod 
   
   use sfc_data_mod, only: &
      sfc_main_type 
      
  use nwp_data_mod , only: &
      nwp_main_type   
      
   implicit none
 
   !***********************************************************************
   ! Marker: DECLARE VARIABLES
   !***********************************************************************
   integer(kind=int4), parameter:: nmask_clavr = 5
   integer(kind=int4):: Nword_Clavr
   integer(kind=int4):: Nword_Clavr_Start
   integer(kind=int4):: Nrec_Avhrr_Header
   integer(kind=int4):: Ifile    
   integer(kind=int4):: File_List_Lun                 !logical unit number for File_list
   integer(kind=int4):: Segment_Number
   integer(kind=int4):: Skip_Processing_Flag
   integer(kind=int4):: Lat_Idx
   integer(kind=int4):: Lon_Idx
   integer(kind=int4):: Zen_Idx
   integer(kind=int4):: Elem_Idx  !generic pixel (along scan) index
   integer(kind=int4):: Line_Idx  !generic line (across scan) index
   integer(kind=int4):: Phase_Called_Flag
   logical:: Level1b_Exists
   integer, parameter:: num_Segment_Time_points=15
   real(kind=real4), dimension(num_Segment_Time_points):: Segment_Time_Point_Seconds
   real(kind=real4) :: Start_Time_Point_Hours
   real(kind=real4) :: End_Time_Point_Hours

   real(kind=real4) :: Segment_Time_Point_Seconds_temp
   real(kind=real4) :: Start_Time_Point_Hours_temp
   real(kind=real4) :: End_Time_Point_Hours_temp

   real(kind=real4) :: Total_Processing_Start_Time_Hours
   real(kind=real4) :: Total_Processing_End_Time_Hours
   real(kind=real4) :: Total_Processing_Time_seconds
   real(kind=real4) :: Orbital_Processing_Start_Time_Hours
   real(kind=real4) :: Orbital_Processing_End_Time_Hours
   real(kind=real4) :: Orbital_Processing_Time_Seconds
   character(len=180):: File_1b_Temp
   integer(kind=int4):: erstat
   real(kind=real4):: Time_Since_Launch
   integer(kind=int4):: err_reposnx_Flag
   integer(kind=int4), parameter:: One = 1
 
   integer(kind=int4):: ios
   integer(kind=int4):: File_Number
   integer(kind=int4):: Ierror_Level1b
   integer(kind=int4):: ierror_Nwp
     
   integer(kind=int4) :: ierror
   
  
   character(*), parameter :: PROGRAM_NAME = 'CLAVRXORB'

   integer, parameter:: LRC_Meander_Flag = 1
   integer, parameter:: Max_LRC_Distance = 10
   real, parameter:: Min_LRC_Jump = 0.0   !0.5
   real, parameter:: Max_LRC_Jump = 100.0 !10.0
   integer, parameter:: Missing_LRC_Value = -999
   integer, parameter:: Grad_Flag_LRC = -1
   real, parameter:: Min_Bt_11um_LRC = 220.0
   real, parameter:: Max_Bt_11um_LRC = 300.0
     
   ! GOES header structures
   TYPE (AREA_STRUCT) :: AREAstr
   TYPE (GVAR_NAV)    :: NAVstr
   
   logical :: dcomp_run
   
   character ( len = 30) :: string_30
   character ( len = 100) :: string_100
   
   type (conf_user_opt_type) :: config
   
   type ( date_type ) :: start_time_obj
   type ( date_type ) :: end_time_obj
   
   !------------- VIIRS variables --------------
   real(kind=real4), dimension(:,:), pointer :: lunar_ref
   real(kind=real4), dimension(:,:),allocatable :: lunar_placeholder
   
   type ( sfc_main_type ) :: sfc_obj
   type ( conf_user_opt_type) :: conf_obj
   type ( nwp_main_type ) :: nwp
   
  
   integer:: Chan_Idx , ii,jj
   
   

   !***********************************************************************
   ! Begin Executable Code
   !***********************************************************************
   call mesg ( '<-- Start of CLAVRXORB --> $Id$' &
      , level = verb_lev % MINIMAL , color = 43 )

   
   !----------------------------------------------------------------------------
   ! Initialize some flags
   !----------------------------------------------------------------------------
   Number_Of_Temporary_Files = 0
   Skip_Processing_Flag = sym%NO

   !----------------------------------------------------------------------------
   ! Determine time of the start of all processing
   !----------------------------------------------------------------------------
   Total_Processing_Start_Time_Hours = COMPUTE_TIME_HOURS()

   !------------------------------------------------------------------------------
   ! initialize previous date variables
   !-------------------------------------------------------------------------------
   Month_Prev = 0
   Sensor%WMO_Id_Previous = 0
   Start_Year_Prev = 0
   Start_Day_Prev = 0

   !*************************************************************************
   ! Marker: Read and Quality Check User Defined Options
   !*************************************************************************
   call config % set_config()
   call SETUP_USER_DEFINED_OPTIONS()
   ! - have to put  this toegther
   
   config % temp_path = trim(temporary_data_dir)
   print*,config % temp_path
   print*,'==============     files to process  ==========='
   do ii = 1, config % n_files
      print*,trim(config % file % infile(ii)), config % file % ETsensor (ii)
   end do   
   print*
   !-----------------------------------------------------------------------
   !--- set up surface radiative properties 
   !TODO    
   !-----------------------------------------------------------------------
   call SETUP_UMD_PROPS()
 
    !--------------------------------------------------------------------
    !--- setup clock corrections in memory
    !  -- this is avhrr only !!
    ! TODO
    !--------------------------------------------------------------------
    if (nav_opt== 2) then
      call SETUP_CLOCK_CORRECTIONS()
    endif

   !**********************************************************************
   ! Marker: Read file directories in FILELIST
   !**********************************************************************

   !--- print to screen which file list is used
   call mesg ( "CLAVR-x FILE LIST FILE USED: "//trim(File_list)  ) 

   Image%Level1b_Path = config % file % l1b_path
   Dir_Level2 = config % file % out_path
  
   !----------------------------------------------------------------------
   ! Marker: BEGIN LOOP OVER FILES
   !----------------------------------------------------------------------
   
   File_loop: do file_number = 1 , config % n_files
 
      !----------------------------------------------------------------------------
      ! Determine time of the start of the processing of this orbit
      !----------------------------------------------------------------------------
      Orbital_Processing_Start_Time_Hours = COMPUTE_TIME_HOURS()
      
      !--------------------------------------------------------------
      ! Determine if this level-1b file can be opended, if not skip
      !--------------------------------------------------------------
      file_1b_temp = config % file % infile (file_number)
      
      if (.not. file_test(config % file % infile_fullpath (file_number))) then
         print *, EXE_PROMPT, "ERROR: Level-1b file not found, skipping"
         cycle file_loop
      end if 

      !--------------------------------------------------------------
      !  Determine based on the L1b file's extension if the input file
      !  is gzipped. Announce this fact thru a boolean var.
      !--------------------------------------------------------------
      call DETERMINE_LEVEL1B_COMPRESSION(File_1b_Temp,L1b_Gzip,L1b_Bzip2)

      !------------------------------------------------------------------------
      ! Determine from which sensor this file comes from (MODIS,AVHRR or VIIRS)
      ! and populate sensor structure
      !------------------------------------------------------------------------
      call DETECT_SENSOR_FROM_FILE(AREAstr,NAVstr,Ierror)
      
      if (Ierror == sym%YES) then
         print *, EXE_PROMPT, "ERROR: Sensor could not be detected, skipping file "
         cycle file_loop
      end if

      !--- print to screen the file name
      call mesg (" " )
      call mesg ("<------------- Next Orbit ---------------> ",level = verb_lev % DEFAULT )

      !-----------------------------------------------------------------------
      ! knowing the file type, determine the expected number of elements and
      ! lines in the file
      ! for AVHRR, determine file type and some record lengths
      ! AVHRR Header is read here
      !-----------------------------------------------------------------------
      call  SET_FILE_DIMENSIONS(Image%Level1b_Full_Name,AREAstr,Nrec_Avhrr_Header,Nword_Clavr,Nword_Clavr_Start,Ierror) 
 
      if (Ierror == sym%YES) then
         print *, EXE_PROMPT, "ERROR: Could not set file dimentions, skipping file "
         cycle file_loop
      end if

      if (Image%Number_Of_Lines <= 0) then
         cycle file_loop    
      end if
  
      !-----------------------------------------------------------------------
      !--- Compute the time stamp for use in all generated HDF output files
      !  AW-2014-12-22 Why now? Why here?
      !-----------------------------------------------------------------------
      call HDF_TSTAMP()
   
      !-----------------------------------------------------------------------
      !--- set up pixel level arrays (size depends on sensor)
      !-----------------------------------------------------------------------
      !--- determine segment size here
      Line_Idx_Min_Segment = 1
      Line_Idx_Max_Segment = Image%Number_Of_Lines_Per_Segment

      !*************************************************************************
      ! Marker:  READ IN HEADER AND DETERMINE SOME CONSTANTS
      !*************************************************************************

      !----------------------------------------------------------------------
      ! Knowing the sensor, interogate files to start, end date and time
      !----------------------------------------------------------------------
      call SET_DATA_DATE_AND_TIME ( AREAstr)
      
      call start_time_obj % set_date_with_doy ( int(image % start_year) , int(image % start_doy) &
         , int(image % start_time / 1000./60./60.) &
         , int(modulo(image % start_time / 1000. / 60. , 60.)) &
         , int(modulo(image % start_time / 1000. , 60.)) ) 
      
       call end_time_obj % set_date_with_doy ( int(image % end_year) , int(image % end_doy) &
         , int(image % end_time / 1000./60./60.) &
         , int(modulo(image % end_time / 1000. / 60. , 60.)) &
         , int(modulo(image % end_time / 1000. , 60.)) )
         
      call start_time_obj % print_data(title='STARTTIME')
      call end_time_obj % print_data(title='STOPTIME')

      !----------------------------------------------------------------------
      ! Output sensor and image parameters to screen
      !----------------------------------------------------------------------
      call OUTPUT_SENSOR_TO_SCREEN()
      call OUTPUT_IMAGE_TO_SCREEN()
      call OUTPUT_PROCESSING_LIMITS_TO_SCREEN()
      
      !------------------------------------------------------------------
      ! Setup PFAAST (FAST IR RTM) for this particular sensor
      !AW   bring this later
      !------------------------------------------------------------------
      call SETUP_PFAAST(Sensor%WMO_Id)

      !------------------------------------------------------------------
      ! Setup Solar-channel RTM terms for this particular sensor
      !------------------------------------------------------------------
      call SETUP_SOLAR_RTM(Sensor%WMO_Id)

      !------------------------------------------------------------------
      ! update settings according sensor ( algo mode and channel settings 
      !    including turn-on and off)
      !------------------------------------------------------------------
      call UPDATE_CONFIGURATION (Sensor%Sensor_Name)
     
      !------------------------------------------------------------------
      ! Create pixel arrays which data for this segment
      !------------------------------------------------------------------
      call CREATE_PIXEL_ARRAYS()

      !------------------------------------------------------------------
      ! Read in Dark Sky Composite
      !------------------------------------------------------------------
      Dark_Composite_Name = "no_file"
      if (Read_Dark_Comp == sym%YES) then
         call DETERMINE_DARK_COMPOSITE_NAME(AREAstr)
      end if
   
      !*************************************************************************
      ! Marker:  READ IN SENSOR-SPECIFIC CONSTANTS
      !*************************************************************************

      !--- read in Instrument Constants from appropriate file
      call READ_INSTR_CONSTANTS()

      !--- read in Algorithm Constants from appropriate file
      call READ_ALGO_CONSTANTS()

   
      !--------------------------------------------------------------------
      ! read in sst analysis file
      !--------------------------------------------------------------------

      !Set to use default options. This will be turned off if there is NO daily OISST data
      Use_Sst_Anal = 0

 
     
 
      !*************************************************************************
      ! Marker:  READ IN NWP DATA
      !*************************************************************************
      call nwp % populate ( start_time_obj, end_time_obj ,  config  % ancil_path , nwp_opt  )
      
    
      
      !--- GFS
      if (Nwp_Opt == 1) then
         call READ_GFS_DATA(Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time, &
                           Image%End_Year, Image%End_Doy, Image%End_Time, Gfs_Data_Dir, ierror_Nwp) 
      end if

      !--- NCEP Reanalysis
      if (Nwp_Opt == 2) then
         call READ_NCEP_REANALYSIS_DATA(Image%Start_Year, Image%Start_Doy, Image%Start_Time,  &
                                        Image%End_Doy, Image%End_Time, Ncep_Data_Dir)
      end if

      !--- CFSR
      if (Nwp_Opt == 3) then
         call READ_GFS_DATA(Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time, &
                            Image%End_Year, Image%End_Doy, Image%End_Time, Cfsr_Data_Dir, ierror_Nwp) 
      endif

      !--- GDAS
      if (Nwp_Opt == 4) then
         call READ_GFS_DATA(Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time,  &
                            Image%End_Year, Image%End_Doy, Image%End_Time, Cfsr_Data_Dir, ierror_Nwp) 
      endif
     
      !---- if NWP is being read in, then proceeed in allocating RTM, NWP arrays
      if (Nwp_Opt /= 0) then
         
         !--- Quality control NWP fields
         call QC_NWP()
        
         !--- create temporary NWP vectors needed for RTM
         call CREATE_TEMP_NWP_VECTORS()

 
         !--- Compute mappings for NWP and RTM vertical coordinates
         call MAP_NWP_RTM(NLevels_Nwp, &
                          P_Std_Nwp, &
                          NLevels_Rtm, &
                          P_Std_Rtm)
         
         !--- allocate RTM structures
         call ALLOCATE_RTM(nlon_Nwp,nlat_Nwp)

      end if

     

      !*************************************************************************
      ! Marker: Populate or read in other lookup tables
      !*************************************************************************

      !--- planck and aerosol tables
      if (Sensor%WMO_Id /= Sensor%WMO_Id_Previous) then

         call  POPULATE_PLANCK_TABLES(a1_20,a2_20,nu_20, &
                                   a1_21,a2_21,nu_21, &
                                   a1_22,a2_22,nu_22, &
                                   a1_27,a2_27,nu_27, &
                                   a1_28,a2_28,nu_28, &
                                   a1_29,a2_29,nu_29, &
                                   a1_30,a2_30,nu_30, &
                                   a1_31,a2_31,nu_31, &
                                   a1_32,a2_32,nu_32, &
                                   a1_33,a2_33,nu_33, &
                                   a1_34,a2_34,nu_34, &
                                   a1_35,a2_35,nu_35, &
                                   a1_36,a2_36,nu_36, &
                                   a1_40,a2_40,nu_40, &
                                   a1_41,a2_41,nu_41)
 
         if (Aer_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR') > 0) then
            call READ_AER_CH123A_REF_LUTS(Ancil_Data_Dir,Sensor%WMO_Id)
         end if

         Sensor%WMO_Id_Previous = Sensor%WMO_Id
 
      end if
 
 

      !--- compute Sun-Earth distance
      Sun_Earth_Distance = 1.0 - 0.016729*cos(0.9856*(Image%Start_Doy-4.0)*dtor)
 
      !--- compute time since launch needed for reflectance calibration
      Time_Since_Launch = Image%Start_Year + Image%Start_Doy / 365.25  - launch_Date
 
      !--------------------------------------------------------------
      !-- Marker: Interpolate a clock correction for this orbit
      !--------------------------------------------------------------
      if ((trim(Sensor%Sensor_Name) == 'AVHRR-1') .or. &
          (trim(Sensor%Sensor_Name) == 'AVHRR-2') .or. &
          (trim(Sensor%Sensor_Name) == 'AVHRR-3')) then

         if(Nav_Opt == 2) then
            call INTERPOLATE_CLOCK_ERROR(Image%Start_Year, Image%Start_Time,  &
                                   Image%End_Year, Image%End_Time, &
                                   Sensor%WMO_Id,Nav%Timerr_Seconds)
           print *, EXE_PROMPT, "Clock Error = ", Nav%Timerr_Seconds
         endif  
      endif

      !*************************************************************************
      ! READ IN DATA RECORDS AND PERFORM QUALITY CHECKING
      !*************************************************************************

      !-------------------------------------------------------------------------------------------
      ! Marker: Begin loop over orbit segments
      !-------------------------------------------------------------------------------------------
      call mesg ( " ")
      call mesg ("Started Processing All Orbital Segments")


      !--- compute number of segments in this orbit 
      if (mod(Image%Number_Of_Lines,Image%Number_Of_Lines_Per_Segment) == 0) then
         Image%Number_Of_Segments = Image%Number_Of_Lines / Image%Number_Of_Lines_Per_Segment
      else
         Image%Number_Of_Segments = Image%Number_Of_Lines / Image%Number_Of_Lines_Per_Segment + 1
      endif

      !--- initialize time counters
      Segment_Time_Point_Seconds = 0.0
      Segment_Time_Point_Seconds_temp = 0.0

      Segment_loop: do Segment_Number = 1,Image%Number_Of_Segments
        
         !--- reset skip processing flag 
         Skip_Processing_Flag = sym%NO

         !--- reset Goes_Scan_Line_Flag
         Goes_Scan_Line_Flag = sym%NO

         !--- reset pixel arrays to missing for this segment
         call RESET_PIXEL_ARRAYS_TO_MISSING()

         !-----------------------------------------------------------------
         !---- Marker: Read level-1b data
         !-----------------------------------------------------------------
         Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

         call READ_LEVEL1B_DATA(Image%Level1b_Full_Name,Segment_Number, &
                                Time_Since_Launch,AREAstr,NAVstr,Nrec_Avhrr_Header,Ierror_Level1b)
                                
     
         if (Ierror_Level1b /= 0) then
            print *, EXE_PROMPT, "ERROR:  Error reading level1b, skipping this file"
            exit
         end if
         
         
          ! - read the sfc data for this segement                       
         conf_obj % ancil_path = trim(Ancil_Data_Dir)
         call sfc_obj % populate ( start_time_obj, nav % lat_1b, nav % lon_1b, config , nwp)
         
      
         !-------------------------------------------------------------------
         ! Modify Chan_On flags to account for channels read in
         !-------------------------------------------------------------------
         call SET_CHAN_ON_FLAG(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
         
         !-------------------------------------------------------------------
         ! Compute Lunar Reflectance
         !-------------------------------------------------------------------
         if (trim(Sensor%Sensor_Name) == 'VIIRS' .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then

            ! - check the angles if this is a good lunar scene
            ! - lun and solar zenith angle

            Lunar_Ref => ch(42)%Ref_Lunar_Toa
                  
            call COMPUTE_LUNAR_REFLECTANCE (ch(42)%Rad_Toa &
                     & , Geo%Solzen, Geo%Lunzen &
                     & , Image%Start_Year, month,day_of_month,Image%start_time &
                     & , Geo%moon_phase_angle  &
                     & , ancil_data_dir &
                     & , Lunar_placeholder)
            
            Lunar_Ref = Lunar_placeholder   
            Lunar_Ref=>null() 
                  
            Geo%Lunrelaz = abs ( Geo%Lunaz - Geo%Sataz )
            where ( Geo%Lunrelaz > 180. ) 
               Geo%Lunrelaz = 360.0 - Geo%Lunrelaz  
            end where
       
         end if
         
         End_Time_Point_Hours = COMPUTE_TIME_HOURS()

         !--- update time summation for level-1b processing
         Segment_Time_Point_Seconds(1) =  Segment_Time_Point_Seconds(1) + &
              60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

         !--- check to see that some data was read in
         if (Image%Number_Of_Lines_Read_This_Segment == 0) then
            print *, EXE_PROMPT, "WARNING: no scans read in, exiting segment processing loop"
            Skip_Processing_Flag = sym%YES
            !exit
         end if

         !---- go no further is no data is read in
         if (Skip_Processing_Flag == sym%NO) then    !skip_Processing_Flag
   
            !*******************************************************************
            ! Marker: Definition of pixel hdf files 
            !  if first segment, define pixel hdf files
            ! this has to be here because cal coefficients are only known
            ! after data has been read in and calibrated
            !*******************************************************************
            if (Segment_Number == 1) then

               !--- place algorithm cvs tags into global strings for output
               call SET_CLOUD_TYPE_VERSION()

               Num_Scans_Level2_Hdf = 0

               call DEFINE_HDF_FILE_STRUCTURES(Image%Number_Of_Lines, &
                              Dir_Rtm, &
                              Dir_Level2, &
                              Image%Level1b_Name, &
                              Rtm_File_Flag, &
                              Level2_File_Flag, &
                              c1,c2,a1_20,a2_20,nu_20, &
                              a1_31,a2_31,nu_31,a1_32,a2_32,nu_32,Solar_Ch20_Nu,&
                              Sun_Earth_Distance,Therm_Cal_1b, &
                              Ref_Cal_1b,Nav_Opt,Use_Sst_Anal, &
                              Modis_Clr_Alb_Flag,Nwp_Opt, &
                              Ch1_Gain_Low,Ch1_gain_High, &
                              Ch1_Switch_Count_Cal,Ch1_Dark_Count_Cal, &
                              Ch2_Gain_low,Ch2_Gain_High, &
                              Ch2_Switch_Count_Cal,Ch2_Dark_Count_Cal, &
                              Ch3a_Gain_low,Ch3a_Gain_High, &
                              Ch3a_Switch_Count_Cal,Ch3a_Dark_Count_Cal, &
                              Image%Start_Year,Image%End_Year,Image%Start_Doy,Image%End_Doy,&
                              Image%Start_Time,Image%End_Time)
            end if
   
            !----------------------------------------------------------------------
            ! check to see if this segment fits desired subset if subsetting
            !----------------------------------------------------------------------
            if (Subset_Pixel_Hdf_Flag == sym%YES) then

               !--- check to see if data fell within Solar zenith angle bounds
               if ((minval(Geo%Solzen(Image%Number_Of_Elements/2,:)) > Geo%Solzen_Max_Limit) .or.   &
                   & (maxval(Geo%Solzen(Image%Number_Of_Elements/2,:)) < Geo%Solzen_Min_Limit)) then
                  print *, "skipping seg for solzen"
                  cycle   Segment_loop
               end if

               !--- check to see if in latitude range (using diag values)
               if ((minval(Nav%Lat_1b(Image%Number_Of_Elements/2,:)) > Nav%Lat_Max_Limit) .or.   &
                  & (maxval(Nav%Lat_1b(Image%Number_Of_Elements/2,:)) < Nav%Lat_Min_Limit)) then
                  print *, "skipping seg for latitude"
                  cycle   Segment_loop 
               end if

            end if

            !*******************************************************************
            ! Marker: Recompute geolocation
            !*******************************************************************
            !--- use level 1b navigation
            Nav%Lat = Nav%Lat_1b
            Nav%Lon = Nav%Lon_1b

            !---  AVHRR Repositioning
            if (trim(Sensor%Sensor_Name) == 'AVHRR-1' .or. &
                trim(Sensor%Sensor_Name) == 'AVHRR-2' .or. &
                trim(Sensor%Sensor_Name) == 'AVHRR-3' ) then

              if (Nav_Opt == 2) then

               if ((Nav%Timerr_Seconds /= Missing_Value_Real4) .and. &
                   (Nav%Timerr_Seconds /= 0.0)) then 

                  call REPOSITION_FOR_CLOCK_ERROR(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment, &
                                  Nav%Timerr_Seconds,Err_Reposnx_Flag)
               else
                  Nav%Lat = Nav%Lat_1b
                  Nav%Lon = Nav%Lon_1b
               end if

             end if

            end if
   
            !--- compute time (local and utc) variables for this segment
            call CONVERT_TIME(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            if (Nwp_Opt /= 0) then

               !--- map each each into correct NWP cell
               call MAP_PIXEL_NWP(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)

               !--- compute needed NWP levels (sfc, tropo, inversion, ...)
               call COMPUTE_NWP_LEVELS_SEGMENT(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)

            end if
   
            !--- determine a pixel-level mask to exclude bad or unprocessible data
            call SET_BAD_PIXEL_MASK(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)

            !--- if the segment is 99% bad data, skip it
            if (Segment_Valid_Fraction < 0.01) then 
                  print *, EXE_PROMPT, "ERROR: Insufficient Data, skipping Segment ", Segment_Valid_Fraction
                  cycle Segment_loop 
            endif

            !*******************************************************************
            ! Marker: Interpolate ancillary data to each pixel in this segment
            !*******************************************************************
            Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

            !--- surface emissivity
            if (use_seebor == sym%YES) then

               

               do Chan_Idx = 20, 36
                  if (Chan_Idx == 26) cycle
                  if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
                      !TO_REMOVE_SFC
                     ch(Chan_Idx)%Sfc_Emiss = sfc_obj%emis(chan_idx) % data 
                     
                  endif
               enddo
            end if

            !--- mandatory fields - check for substitution of Bad_Pixel for space 

            !--- surface type
            !TO_REMOVE SFC
             sfc % sfc_type = sfc_obj % sfc_type % data  
              
            !--- surface elevation
            if (Read_Surface_Elevation /= 0) then
                 !TO_REMOVE_SFC 
               Sfc%Zsfc_Hires = real(sfc_obj % elevation % data ,kind=real4)
            end if

            !--- merge with nwp surface elevation
            if (Nwp_Opt /= 0) then
                call MERGE_NWP_HIRES_ZSFC(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
            endif

            !--- read coast mask
            if (Read_Coast_Mask == sym%YES) then
                !TO_REMOVE_SFC 
               sfc % coast = sfc_obj % coast_mask % data
            end if
            
         
            !--- read land mask
            if (Read_Land_Mask == sym%YES) then
                Sfc%Land = sfc_obj % land_class % data
            end if
  
            !--- modify land class with ndvi if available (helps mitigate navigation errors)
            call MODIFY_LAND_CLASS_WITH_NDVI(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            !--- read volcano mask
            if (Read_Volcano_Mask == sym%YES) then
               sfc%volcano_mask = sfc_obj % volcano % data
            end if

            !--- read Snow mask
            if (Read_Snow_Mask == sym%READ_SNOW_HIRES .and. sfc_obj % snow_class % is_set ) then
               Sfc%Snow_Hires = sfc_obj % snow_class % data
            end if
   
          
            !--- define binary land and coast masks (yes/no) from land and coast flags
            call COMPUTE_BINARY_LAND_COAST_MASKS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)   

            !---- ensure missing values for space scenes
            where (Space_Mask == sym%YES) 
               Sfc%Zsfc_Hires = Missing_Value_Real4
               Sfc%Coast = Missing_Value_Int1
               Sfc%Land = Missing_Value_Int1
               Sfc%Snow_Hires = Missing_Value_Int1
               Sfc%Snow_Glob = Missing_Value_Int1
               Sfc%Volcano_Mask = Missing_Value_Int1
               Sfc%Sfc_Type = Missing_Value_Int1
            end where

            !--- interpolate sst analyses to each pixel
            if (Use_Sst_Anal == 1) then
               sst_anal = sfc_obj % sst % data
               sst_anal_uni = sfc_obj % sst % stdv
               sst_anal_cice = sfc_obj % sea_ice_fraction % data
              
            end if

            !--- compute a coast mask relative to nwp data
            if (Nwp_Opt /= 0) then
               call COMPUTE_COAST_MASK_NWP(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
            end if

           
            
            !TO_REMOVE SFC
           
            do chan_idx =1,7 
               if (chan_idx == 3 .or. chan_idx ==4 ) cycle
               if (Sensor%Chan_On_Flag_Default(chan_idx) == sym%YES) then
                  if (  allocated ( ch(chan_idx)%Sfc_Ref_White_Sky )) deallocate ( ch(chan_idx)%Sfc_Ref_White_Sky )
                   allocate(ch(chan_idx)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                    
                  ch(chan_idx)%Sfc_Ref_White_Sky = sfc_obj%modis_w_sky(chan_idx) % data
               end if   
            end do
            
            
          
 
            !--- post process dark composite if one read in
            if (Read_Dark_Comp == sym%YES .and. Dark_Composite_Name /= "no_file") then
               call POST_PROCESS_GOES_DARK_COMPOSITE(Ref_Ch1_Dark_Composite, Sfc%Land)
            endif

            !--- check ancillary data
            call QUALITY_CONTROL_ANCILLARY_DATA(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            !-----------------------------------------------------------------------
            !--- Marker: Compute some fundamental pixel-level arrays
            !-----------------------------------------------------------------------

            !--- compute some common used pixel arrays 
            call COMPUTE_PIXEL_ARRAYS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            if (index(Sensor%Sensor_Name, 'VIIRS') > 0) then
               call COMPUTE_VIIRS_SST(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
            end if

            !--- normalize reflectances by the solar zenith angle and sun-earth distance
            call NORMALIZE_REFLECTANCES(Sun_Earth_Distance,Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
   
            !--- compute the channel 3b albedo arrays
            call CH3B_ALB(Sun_Earth_Distance,Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            !--- compute pixel level Snow map based on all ancillary data
            if (Nwp_Opt /= 0) then
               call COMPUTE_SNOW_FIELD(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
            end if

            !--- interpolate surface type field to each pixel in segment
            call GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)   

            !--- compute desert mask cloud detection
            Sfc%Desert_Mask =  DESERT_MASK_FOR_CLOUD_DETECTION(ch(20)%Sfc_Emiss, Nav%Lat, Sfc%Snow, Sfc%Sfc_Type)

            !--- compute city mask cloud detection
            if (Sensor%Chan_On_Flag_Default(42) == sym%YES) then
             Sfc%City_Mask =  CITY_MASK_FOR_CLOUD_DETECTION(ch(42)%Rad_Toa, Sfc%Sfc_Type)
            endif

            End_Time_Point_Hours = COMPUTE_TIME_HOURS()
                        

            !--- update time summation for level-1b processing
            Segment_Time_Point_Seconds(2) =  Segment_Time_Point_Seconds(2) + &
                60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            !*******************************************************************
            ! Marker: Compute nwp mapping and rtm values for each pixel in segment
            !*******************************************************************

            !--- needs an NWP to run.
            if (Nwp_Opt /= 0) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               !-- temporally interp skin temp for each segment (only ncep reanalysis)
               if (Nwp_Opt == 2) then
                  call TEMPORAL_INTERP_TMPSFC_NWP(Scan_Time(Line_Idx_Min_Segment), Scan_Time(Line_Idx_Min_Segment+Image%Number_Of_Lines_Read_This_Segment-1))
               end if

               !--- compute desired nwp parameters 
               call COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS()
               call COMPUTE_PIXEL_NWP_PARAMETERS(Smooth_Nwp_Flag)

               !--- compute a surface temperature from the NWP
               call MODIFY_TSFC_NWP_PIX(1,Image%Number_Of_Elements,Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               !--- compute pixel-level rtm parameters 
               Start_Time_Point_Hours_temp = COMPUTE_TIME_HOURS()
               call GET_PIXEL_NWP_RTM(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               End_Time_Point_Hours_temp = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds_temp =  Segment_Time_Point_Seconds_temp + &
                   & 60.0*60.0*(End_Time_Point_Hours_temp - Start_Time_Point_Hours_Temp)

               !--- apply atmospheric correction - needs rtm results
               call ATMOS_CORR(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               !--- compute surface products (Tsfc,Ndvi,Rsr ...)
               call SURFACE_REMOTE_SENSING(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(3) =  Segment_Time_Point_Seconds(3) + &
                  &  60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)
            end if

            !*******************************************************************
            ! Marker: Spatial metrics processing
            !*******************************************************************
            Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

            !--- populate needed spatial uniformity arrays
            call COMPUTE_SPATIAL_UNIFORMITY(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            !--- call spatial correlation routines
            call COMPUTE_SPATIAL_CORRELATION_ARRAYS()

            if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then

               !--- apply median filter to Bt_Ch20
               call COMPUTE_MEDIAN_SEGMENT(ch(20)%Bt_Toa,Bad_Pixel_Mask, One, &
                              1, Image%Number_Of_Elements,Line_Idx_Min_Segment,Line_Idx_Min_Segment+Image%Number_Of_Lines_Read_This_Segment-One, &
                              Bt_Ch20_Median_3x3,  &
                              Bt_Ch20_Std_Median_3x3)
 
               !--- apply median filter to Ems_Ch20
               call COMPUTE_MEDIAN_SEGMENT(Ems_Ch20,Bad_Pixel_Mask, One, &
                              1, Image%Number_Of_Elements,Line_Idx_Min_Segment,Line_Idx_Min_Segment+Image%Number_Of_Lines_Read_This_Segment-One, &
                              Ems_Ch20_Median_3x3,  &
                              Ems_Ch20_Std_Median_3x3)
            end if

            !--- populate local radiative center arrays
            if (LRC_Flag == sym%YES) then
               if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then

                  call LOCAL_LINEAR_RADIATIVE_CENTER(sym%YES,sym%NO, &
                                       LRC_Meander_Flag, &
                                       ch(31)%Bt_Toa(:,Line_Idx_Min_Segment:Image%Number_Of_Lines_Read_This_Segment-Line_Idx_Min_Segment+1), &
                                       1, &
                                       Image%Number_Of_Elements, &
                                       Line_Idx_Min_Segment,  &
                                       Image%Number_Of_Lines_Read_This_Segment, &
                                       Max_LRC_Distance,  &
                                       Min_LRC_Jump,  &
                                       Max_LRC_Jump,  &
                                       Grad_Flag_LRC,  &
                                       Missing_Value_Int4, &
                                       Bad_Pixel_Mask(:,Line_Idx_Min_Segment:Image%Number_Of_Lines_Read_This_Segment-Line_Idx_Min_Segment+1), &
                                       Min_Bt_11um_LRC,  &
                                       Max_Bt_11um_LRC, &
                                       I_LRC(:,Line_Idx_Min_Segment:Image%Number_Of_Lines_Read_This_Segment-Line_Idx_Min_Segment+1), &
                                       J_LRC(:,Line_Idx_Min_Segment:Image%Number_Of_Lines_Read_This_Segment-Line_Idx_Min_Segment+1))
               end if
            end if

            End_Time_Point_Hours = COMPUTE_TIME_HOURS()
            Segment_Time_Point_Seconds(4) =  Segment_Time_Point_Seconds(4) + &
               &  60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            !--- compute the glint mask
            call COMPUTE_GLINT()

            !*******************************************************************
            ! Marker: Generate pixel-level products
            !*******************************************************************

            !---- pixel level aerosol
            if (index(Sensor%Sensor_Name,'AVHRR') > 0 .and. Aer_Flag == sym%YES) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               call PIXEL_AER_RET_OCEAN(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(5) =  Segment_Time_Point_Seconds(5) + &
                    60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)
            end if

            !--- only apply cloud mask and type routines if nwp/rtm information available
            if (Cld_Flag == sym%YES .and. Nwp_Opt > 0) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               !--- simple cloud optical depth
!              call COMPUTE_SIMPLE_COD(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)               

               !--- cloud mask
               if (Cloud_Mask_Aux_Flag /= sym%USE_AUX_CLOUD_MASK) then
                  if (Cloud_Mask_Bayesian_Flag == sym%YES) then
                     call NB_CLOUD_MASK_BRIDGE(Segment_Number)
                  else
                     print *, "Only the Bayesian Cloud Mask is available, check selection"
                     stop 
                  end if
               end if

               if (Cloud_Mask_Aux_Flag == sym%USE_AUX_CLOUD_MASK .and. Cloud_Mask_Aux_Read_Flag == sym%YES) then
                  Cld_Mask = Cld_Mask_Aux
               end if

               !--- Compute Adjacent pixels Cloud Mask
               call ADJACENT_PIXEL_CLOUD_MASK(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               !--- if dark composite available for GOES-1km data, do this
               if (trim(Sensor%Sensor_Name) == 'GOES-IL-IMAGER' .or. trim(Sensor%Sensor_Name) == 'GOES-MP-IMAGER') then
                 if (Read_Dark_Comp == sym%YES .and. Dark_Composite_Name /= "no_file") then
                   call DARK_COMPOSITE_CLOUD_MASK(Cld_Mask)
                 endif
               end if

               !--- accumulate cloud mask metrics
               call COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS(Cloud_Mask_Count,Nonconfident_Cloud_Mask_Count)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(6) =  Segment_Time_Point_Seconds(6) + &
                     & 60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

               !--- cloud type 
               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               if (Cloud_Mask_Aux_Flag == sym%USE_AUX_CLOUD_MASK .and. trim(Sensor%Sensor_Name) == 'VIIRS') then

                  Cld_Type = Cld_Type_Aux
                  Cld_Phase = Cld_Phase_Aux
                  Phase_Called_Flag = sym%YES

               else  
                  call CLOUD_TYPE_BRIDGE()
                  Phase_Called_Flag = sym%YES
               end if

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
                  Segment_Time_Point_Seconds(7) =  Segment_Time_Point_Seconds(7) + &
                  & 60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            end if   !end of Cld_Flag check

            !--------------------------------------------------------------------
            !   Compute Cloud Properties (Height, Optical Depth, ...)
            !--------------------------------------------------------------------
            if (Cld_Flag == sym%YES .and. Nwp_Opt > 0) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
              

               !-------------------------------------------------------------------
               ! make co2 slicing height from sounder with using sounder/imager IFF
               !-------------------------------------------------------------------
               if (index(Sensor%Sensor_Name,'IFF') > 0) then

                   call CO2_SLICING_CLOUD_HEIGHT(Image%Number_Of_Elements,Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment, &
                                    P_Std_Rtm,Cld_Type, &
                                    Pc_Cirrus_Co2,Tc_Cirrus_Co2)
               endif

               if (ACHA%Mode == 0) then
                  call MODE_ZERO_CLOUD_HEIGHT(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               endif

               if (ACHA%Mode > 0) then 

                  !--- AWG CLoud Height Algorithm (ACHA) and associated products
                  call AWG_CLOUD_HEIGHT_BRIDGE()

                  !--- interpolate NWP wind profiles at cloud-top level
                  call COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

                  !--- interpolate ACHA cloud heights to flight level altitude.
                  call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

                  !--accumulate performance metrics
                  call COMPUTE_ACHA_PERFORMANCE_METRICS(ACHA%Processed_Count,ACHA%Valid_Count,ACHA%Success_Fraction)

               end if

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(8) =  Segment_Time_Point_Seconds(8) + &
                   &   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

    
               ! - lunar reflectance
               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               if (trim(Sensor%Sensor_Name) == 'VIIRS' .and. Sensor%Chan_On_Flag_Default(42) == sym % yes) then
                  if ( count (ch(42)%Ref_Lunar_Toa > 0) > 0 ) then
                     call awg_cloud_nlcomp_algorithm (  Iseg_In=Segment_Number) 
                  end if   
               end if   
          
               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(15) =  Segment_Time_Point_Seconds(15) + &
                     & 60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)  
     
               !--- cloud optical depth and effective radius from vis/nir approach
               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
    
               if (Dcomp_Mode > 0) then
        
                  call AWG_CLOUD_DCOMP_ALGORITHM( Iseg_In = Segment_Number , dcomp_run = dcomp_run)
                  call SET_DCOMP_VERSION()
                  if ( dcomp_run ) then
                     call COMPUTE_CLOUD_WATER_PATH(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
                     call COMPUTE_DCOMP_INSOLATION(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment,Sun_Earth_Distance)
                     call COMPUTE_ADIABATIC_CLOUD_PROPS(Line_Idx_Min_segment,Image%Number_Of_Lines_Read_This_Segment)
                     call COMPUTE_DCOMP_PERFORMANCE_METRICS(DCOMP_Processed_Count,DCOMP_Valid_Count)
                  end if
                  
               end if
    
               !--- compute precipation from optical properties
               if (Dcomp_Mode > 0 ) then
                  call COMPUTE_PRECIPITATION(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               end if

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(9) =  Segment_Time_Point_Seconds(9) + &
                  &  60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

               !--- CLoud Base Height Algorithm if ACHA was executed
               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               if (ACHA%Mode > 0) then 
                 call CLOUD_BASE_BRIDGE()
               endif
               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(10) =  Segment_Time_Point_Seconds(10) + &
                   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            end if

            !--- Non-cloud Detection
            if (Aer_Flag == sym%YES) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               !-->        call DUST_DETECTION_ALGORITHM(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               !-->        call SMOKE_DETECTION_ALGORITHM(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(11) =  Segment_Time_Point_Seconds(11) + &
                   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)
            end if

            !--- Volcanic ash - only for AVHRR/2 and AVHRR/3
            if (Ash_Flag > 0) then
               if (trim(Sensor%Sensor_Name) == 'AVHRR_2' .or. trim(Sensor%Sensor_Name) == 'AVHRR_3') then 
                  start_time_point_hours = COMPUTE_TIME_HOURS()
                  end_time_point_hours = COMPUTE_TIME_HOURS()
                  segment_time_point_seconds(11) =  segment_time_point_seconds(11) + &
                     60.0*60.0*(end_time_point_hours - start_time_point_hours)
               end if

            end if
    
            !--- radiative flux parameters
            Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
            
            if (AVHRR_1_Flag == sym%NO) then    !currently, no AVHRR/1 algorithm
               call COMPUTE_ERB(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
            end if

            !---  Run SASRAB
            if ( Sasrab_Flag == sym%YES) call INSOLATION(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            
            End_Time_Point_Hours = COMPUTE_TIME_HOURS()
            Segment_Time_Point_Seconds(12) =  Segment_Time_Point_Seconds(12) + &
                   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            if (Nwp_Opt > 0) then

               !--- assign clear sky quality flags
               call ASSIGN_CLEAR_SKY_QUALITY_FLAGS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            end if

            !--- generated cloud masked sst field
            if (Nwp_Opt > 0) then
                call COMPUTE_MASKED_SST(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
            end if

            !*******************************************************************
            ! Marker: Write to output files (pixel-level)
            !*******************************************************************
            Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

            call WRITE_PIXEL_HDF_RECORDS(Rtm_File_Flag,Level2_File_Flag)
            
            End_Time_Point_Hours = COMPUTE_TIME_HOURS()
            Segment_Time_Point_Seconds(13) =  Segment_Time_Point_Seconds(13) + &
                   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            !*************************************************************************
            ! Marker: RTM Structure Memory Deallocation
            !*************************************************************************

            !--- deallocate rtm profile arrays (only do if NWP is used)
            if (Nwp_Opt /= 0) then

               do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
                  do Elem_Idx = 1, Image%Number_Of_Elements
                     Lon_Idx = i_Nwp(Elem_Idx,Line_Idx)
                     Lat_Idx = j_Nwp(Elem_Idx,Line_Idx)
                     Zen_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)
                     !--- check to see if Lon_Idx and Lat_Idx are valid, if deallocate
                     if ((Lon_Idx >= 0) .and. (Lon_Idx <= Nlon_Nwp) .and. &
                         (Lat_Idx >= 0) .and. (Lat_Idx <= Nlat_Nwp) .and. &
                         (Zen_Idx >= 0) .and. (Zen_Idx <= Rtm_Nvzen)) then

                        if (Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx_Rtm(Elem_Idx,Line_Idx))%flag == sym%YES) then
                           call DEALLOCATE_RTM_VARS(Lon_Idx,Lat_Idx,Zen_Idx_Rtm(Elem_Idx,Line_Idx))
                        end if

                     end if
                  end do
               end do

               !--- deallocate rtm cells
               do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Image%Number_Of_Lines_Read_This_Segment - 1
                  do Elem_Idx = 1, Image%Number_Of_Elements

                     Lon_Idx = i_Nwp(Elem_Idx,Line_Idx)
                     Lat_Idx = j_Nwp(Elem_Idx,Line_Idx)

                     !--- check to see if Lon_Idx and Lat_Idx are valid, if deallocate
                     if ((Lon_Idx >= 0) .and. (Lon_Idx <= nlon_Nwp) .and. &
                           & (Lat_Idx >= 0) .and. (Lat_Idx <= nlat_Nwp)) then

                        if (Rtm(Lon_Idx,Lat_Idx)%Flag == sym%YES) then
                           call DEALLOCATE_RTM_CELL(Lon_Idx,Lat_Idx)
                        end if

                     end if
                  end do
               end do

            end if

            !--- screen output to mark progress through orbit
           write ( string_100, '(A22, I2, A16, I5 , 4X , I5)')  "processed segment #",Segment_Number," scanlines = ",  &
                        Scan_Number(Line_Idx_Min_Segment), Scan_Number(Image%Number_Of_Lines_Read_This_Segment)
           call mesg  ( string_100 )

         end if   !end Skip_Processing_Flag condition
        !*************************************************************************
        ! Marker: End of loop over orbital segments
        !*************************************************************************
         call sfc_obj % deallocate_all()
         call nwp % deallocate_sat_grid()
      end do Segment_loop

      call mesg ( "Finished Processing All Orbital Segments")
      call mesg ( " ")


      !*************************************************************************
      !   Marker: Close output pixel-level files
      !*************************************************************************

      !*************************************************************************
      ! Marker: Close non-static ancillary data files
      !*************************************************************************


      !*************************************************************************
      !Marker: Deallocate remaining arrays
      !*************************************************************************
      !--- main RTM structures and NWP arrays
      if (Nwp_Opt > 0) then
         call DESTROY_NWP_ARRAYS()
         call DESTROY_TEMP_NWP_VECTORS()
         call DEALLOCATE_RTM()
      end if

      !--- Deallocate memory from pixels arrays
      call DESTROY_PIXEL_ARRAYS()

      !--- remove files in temporary directory
      if (Number_Of_Temporary_Files > 0) then 
         do Ifile = 1, Number_Of_Temporary_Files
            print *, EXE_PROMPT, "Removing Temporary File: ", trim(Temporary_File_Name(Ifile))
            call system("rm "//trim(Temporary_Data_Dir)//trim(Temporary_File_Name(Ifile)))
         end do
      end if
      Number_Of_Temporary_Files = 0   !reset for next file

      !--- Determine time of the end of the processing of this orbit
      Orbital_Processing_End_Time_Hours = COMPUTE_TIME_HOURS()
      Orbital_Processing_Time_Seconds = 60.0*60.0*(Orbital_Processing_End_Time_Hours -  &
                                                Orbital_Processing_Start_Time_Hours)
      Orbital_Processing_Time_Minutes = Orbital_Processing_Time_Seconds/60.0

      !--- write algorithm attributes to level2
      call WRITE_ALGORITHM_ATTRIBUTES()

      !--- close pixel level hdf files
      call CLOSE_PIXEL_HDF_FILES(Rtm_File_Flag,Level2_File_Flag)

      !--- diagnostic screen output
      call mesg ("<----- Timing Results ----->")
      call mesg ("Time for Level-1b Processing (sec) = ", Segment_Time_Point_Seconds(1))
      call mesg ("Time for Ancil. Data Processing (sec) = ", Segment_Time_Point_Seconds(2))
      call mesg ("Time for RTM Processing (sec) = ", Segment_Time_Point_Seconds(3))
      call mesg ("Time for Spatial Processing (sec) = ", Segment_Time_Point_Seconds(4))
      call mesg ("Time for Aerosol Retrieval (sec) = ", Segment_Time_Point_Seconds(5))
      call mesg ("Time for Cloud Mask (sec) = ", Segment_Time_Point_Seconds(6))
      call mesg ("Time for Cloud Type (sec) = ", Segment_Time_Point_Seconds(7))
      call mesg ("Time for Cloud Height (sec) = ", Segment_Time_Point_Seconds(8))
      call mesg ("Time for Cloud Opt/Micro (sec) = ", Segment_Time_Point_Seconds(9))
      call mesg ("Time for NL-COMP (sec) = ", Segment_Time_Point_Seconds(15))
      call mesg ("Time for Cloud Base (sec) = ", Segment_Time_Point_Seconds(10))
      call mesg ("Time for Volcanic Ash (sec) = ", Segment_Time_Point_Seconds(11))
      call mesg ("Time for Earth Radiation Budget (sec) = ", Segment_Time_Point_Seconds(12))
      call mesg ("Time for Pixel-HDF Write (sec) = ", Segment_Time_Point_Seconds(13))
      !  call mesg ("Time for Grid-cell Compilation (sec) = ", Segment_Time_Point_Seconds(14)
      call mesg ("Total Time for Processing This Orbit (sec) = ", Orbital_Processing_Time_Seconds,level=verb_lev % MINIMAL)
      !  print *, EXE_PROMPT, "Temp Time for Processing This Orbit (sec) = ", Segment_Time_Point_Seconds_temp

      !--- add processing time to global attributes
  

      !*************************************************************************
      ! Marker: End loop over files
      !*************************************************************************
      call nwp % deallocate()
   end do File_Loop

   !*************************************************************************
   ! Marker: Close FILELIST
   !*************************************************************************
   close(unit=File_List_Lun)

   !----------------------------------------------------------------------
   ! DEALLOCATE MEMORY AND END PROGRAM
   !----------------------------------------------------------------------

   
   !*************************************************************************
   ! Marker: Final remaining memory deallocation
   !*************************************************************************

   !--- Determine time of the start of all processing
   Total_Processing_End_Time_Hours = COMPUTE_TIME_HOURS()
   Total_Processing_Time_seconds = 60.0*60.0*(Total_Processing_End_Time_Hours -  &
                                              Total_Processing_Start_Time_Hours)
    call mesg ( "Total Time for All Processing (sec) = ", Total_Processing_Time_seconds,level=verb_lev % MINIMAL)

   !---- print to screen that processing is done
   call mesg (  "<--------- End of CLAVRXORB ---------->",level=verb_lev % MINIMAL)

   stop
!--------- PFAAST error statements
 200   print *, EXE_PROMPT, "error reading PFAAST RTM coefficients, stopping"

!*************************************************************************
! Marker: End of code
!*************************************************************************
 end program PROCESS_CLAVRX
