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
! Version 2015a - Stable Release for CSPP
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
! Channel 37 = ABI Channel 8 = 6.2 micron
! Channel 38 = ABI Channel 13 = 10.4 micron
! Channels 39-44 are defined only for VIIRS
! Channel 39 - VIIRS I1 - 0.64 micron
! Channel 40 - VIIRS I2 - 0.865 micron
! Channel 41 - VIIRS I3 - 1.61 micron
! Channel 42 - VIIRS I4 - 3.74 micron
! Channel 43 - VIIRS I5 - 11.45 micron
! Channel 44 - VIIRS DNB - 0.7 micron
!
!-------------------------------------------------------------------------
 
!*****************************************************************************
! Marker: ACCESS MODULES 
!******************************************************************************
   use ACHA_CLAVRX_BRIDGE, only: &
      AWG_CLOUD_HEIGHT_BRIDGE
   
   use AWG_CLOUD_HEIGHT, only: &
      local_linear_radiative_center
      
   use AVHRR_PIXEL_AEROSOL, only: &
      pixel_aer_ret_ocean &
      , read_aer_ch123a_ref_luts
      
   use ASOS_CLAVRX_BRIDGE, only: &
      asos_bridge
   
   use AVHRR_REPOSITION_ROUTINES, only: &
      interpolate_clock_error &
      , reposition_for_clock_error &
      , setup_clock_corrections
   
   use BASELINE_CLOUD_MASK, only: &
      BASELINE_CLOUD_MASK_MAIN
   
   use calibration_constants, only: &
      Solar_Ch20_Nu &
      ,launch_Date &
      , sun_earth_distance 
   
   use CCL_CLAVRX_BRIDGE, only: &
      ccl_bridge
   
   use CLAVRX_MESSAGE_MODULE, only: &
      MESG &
      , VERB_LEV
   
   use CLAVRX_OLR_MODULE, only: &
      compute_olr &
      , setup_olr
   
   use CLAVRX_SST_MODULE, only: &
    compute_masked_sst &
    , compute_sst &
    , setup_sst
   
   use CLOUD_BASE_CLAVRX_BRIDGE, only: &
      cloud_base_bridge
   
   use CLOUD_TYPE_BRIDGE_MODULE, only: &
    cloud_type_bridge
   
   use CLOUD_HEIGHT_ROUTINES, only: &
      co2_slicing_cloud_height &
      , co2irw_cloud_height &
      , compute_altitude_from_pressure &
      , compute_cloud_top_level_nwp_wind_and_tpw &
      , compute_csbt_cloud_masks &
      , convective_cloud_probability &
      , CTP_MULTILAYER &
      , h2o_cloud_height &
      , make_cirrus_prior_temperature &
      , mode_zero_cloud_height &
      , modify_cloud_type_with_sounder  &
      , opaque_cloud_height &
      , opaque_transmission_height &
      , sounder_emissivity &
      , supercooled_cloud_probability
   
   use CX_CONSTANTS_MOD,only: &
      int1, int2,  int4, real4 &
      , DTOR &
      , Missing_Value_Int1 &
      , missing_value_int4 &
      , missing_value_real4 &
      , SOLAR_OBS_TYPE &
      , LUNAR_OBS_TYPE &
      , MIXED_OBS_TYPE &
      , THERMAL_OBS_TYPE &
      , sym &
      , EXE_PROMPT &
      , Nchan_Clavrx &
      , SFC_TYPE_SDS_NAME &
      , VOLCANO_MASK_SDS_NAME &
      , LAND_MASK_SDS_NAME &
      , COAST_MASK_SDS_NAME &
      , SURFACE_ELEV_SDS_NAME &
      , SNOW_MASK_SDS_NAME &
      , MODIS_ALB_2_13_SDS_NAME &
      , MODIS_ALB_1_64_SDS_NAME &
      , MODIS_ALB_1_24_SDS_NAME &
      , MODIS_ALB_0_66_SDS_NAME &
      , MODIS_ALB_0_86_SDS_NAME
   
   use cx_conf_mod , only:  &
      conf_main_type &
      , set_config
   
   use date_tools_mod, only: &
         leap_year_fct &
         , compute_month &
         , compute_day &
         , compute_time_hours
   
   
   use DCOMP_DERIVED_PRODUCTS_MODULE, only: &
      compute_adiabatic_cloud_props &
      , compute_cloud_water_path &
      , compute_dcomp_insolation &
      , compute_precipitation
   
   use DNB_RETRIEVALS_MOD, only: &
      COMPUTE_LUNAR_REFLECTANCE
   
   use DNCOMP_CLAVRX_BRIDGE_MOD, only: &
      AWG_CLOUD_DNCOMP_ALGORITHM &
      , awg_cloud_dncomp_algorithm_iband &
      , set_dcomp_version
   
   use file_tools, only: &
      file_test &
      , get_lun
   
   use GFS, only: &
    read_gfs_data
   
   use GLOBSNOW_READ_ROUTINES, only: &
      get_globsnow_filename &
      , get_pixel_globsnow_analysis &
      , read_globsnow_analysis_map
   
   use GOES_MODULE, only: &
         gvar_nav &
      , dark_composite_cloud_mask &
      , determine_dark_composite_name &
      , post_process_goes_dark_composite
      
   use CX_SSEC_AREAFILE_MOD,only: &
      area_header_type   
   
   use LAND_SFC_PROPERTIES, only: &
      land_grid_description &
      , get_snow_map_filename &
      , close_land_sfc_hdf &
      , open_land_sfc_hdf &
      , read_land_sfc_hdf
   
   use LASZLO_INSOLATION, only: &
      insolation
   
   use LEVEL2_ROUTINES, only: &
      write_pixel_hdf_records
   
   use NB_CLOUD_MASK_CLAVRX_BRIDGE, only: &
       NB_CLOUD_MASK_BRIDGE
   
   use NCEP_REANALYSIS, only: &
      read_ncep_reanalysis_data
   
   use NWP_COMMON   
   
   use NUMERICAL_TOOLS_MOD, only: &
           compute_median_segment
   
   use MODIS_MODULE, only:
   
    use OCA_MODULE, only: &
       READ_OCA
   
   use OISST_ANALYSIS, only: &
      GET_OISST_MAP_FILENAME &
      , get_pixel_sst_analysis &
      , read_oisst_analysis_map
   
   use PIXEL_COMMON
   
   
   use PIXEL_ROUTINES, only: &
        DESERT_MASK_FOR_CLOUD_DETECTION &
      , CITY_MASK_FOR_CLOUD_DETECTION &
      , ADJACENT_PIXEL_CLOUD_MASK &
      , ASSIGN_CLEAR_SKY_QUALITY_FLAGS &
      , ATMOS_CORR &
      , CH20_PSEUDO_REFLECTANCE &
      , compute_cloud_mask_performance_metrics &
      , compute_acha_performance_metrics &
      , compute_dcomp_performance_metrics &
      , compute_glint &
      , compute_pixel_arrays &
      , compute_snow_class &
      , compute_snow_class_nwp &
      , compute_snow_class_oisst &
      , compute_spatial_correlation_arrays &
      , COMPUTE_SPATIAL_UNIFORMITY &
      , convert_time &
      , determine_level1b_compression &
      , expand_space_mask_for_user_limits &
      , merge_nwp_hires_zsfc &
      , modify_land_class_with_ndvi &
      , normalize_reflectances &
      , quality_control_ancillary_data &
      , read_modis_white_sky_albedo &
      , set_bad_pixel_mask &
      , set_chan_on_flag &
      , set_solar_contamination_mask &
      , surface_remote_sensing
   
   use PLANCK, only: &
    populate_planck_tables &
    , populate_planck_tables_sounder
   

   use RT_UTILITIES, only: &
        RTM_NVZEN &
      , SETUP_SOLAR_RTM &
      , MAP_NWP_RTM &
      , CREATE_TEMP_NWP_VECTORS &
      , DESTROY_TEMP_NWP_VECTORS &
      , GET_PIXEL_NWP_RTM &
      , ALLOCATE_RTM &
      , DEALLOCATE_RTM &
      , DEALLOCATE_RTM_VARS &
      , DEALLOCATE_RTM_CELL 
   
   use RTM_COMMON, only: &
      NLEVELS_RTM &
      , P_Std_Rtm &
      , rtm
   
   use SENSOR_MODULE, only: &
    detect_sensor_from_file &
    , output_image_to_screen &
    , output_processing_limits_to_screen &
    , output_sensor_to_screen &
    , read_instr_constants &
    , read_level1b_data &
    , set_data_date_and_time &
    , set_file_dimensions
      
   use SFC_EMISS, only: &
      close_seebor_emiss &
      , open_seebor_emiss &
      , read_seebor_emiss
   
   use SIMPLE_COD, only: &
    compute_simple_lunar_cod &
    , compute_simple_solar_cod
      
   use SFC_PROP_UMD_MOD, only: &
    compute_binary_land_coast_masks &
    , get_pixel_sfc_emiss_from_sfc_type &
    , setup_umd_props
   
   use USER_OPTIONS, only: &
    setup_user_defined_options &
    , update_configuration
   
        
  
   implicit none
   
 
   !***********************************************************************
   ! Marker: DECLARE VARIABLES
   !***********************************************************************
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
   real(kind=real4) :: Total_Processing_Time_minutes
   real(kind=real4) :: Orbital_Processing_Start_Time_Hours
   real(kind=real4) :: Orbital_Processing_End_Time_Hours
   real(kind=real4) :: Orbital_Processing_Time_Seconds
   character(len=1020):: File_1b_Temp
   integer(kind=int4):: erstat
   real(kind=real4):: Time_Since_Launch
   integer(kind=int4):: err_reposnx_Flag
   integer(kind=int4), parameter:: One = 1
 
   integer(kind=int4):: ios
   integer(kind=int4):: File_Number
   integer(kind=int4):: Ierror_Level1b
   integer(kind=int4):: ierror_Nwp
   integer(kind=int4):: iperiod16   
   integer(kind=int4) :: ierror
   character(len=3):: Day_String
   character(len=1020):: Modis_White_Sky_0_66_Name
   character(len=1020):: Modis_White_Sky_0_86_Name
   character(len=1020):: Modis_White_Sky_1_24_Name
   character(len=1020):: Modis_White_Sky_1_64_Name
   character(len=1020):: Modis_White_Sky_2_13_Name
   character(len=1020):: Snow_Mask_File_Name
   character(len=1020):: oiSst_File_Name
   

   integer(kind=int4):: Emiss_File_Id = missing_value_int4
   integer(kind=int4):: Coast_Mask_Id = missing_value_int4
   integer(kind=int4):: Land_Mask_Id = missing_value_int4
   integer(kind=int4):: Sfc_Type_Id = missing_value_int4
   integer(kind=int4):: Volcano_Mask_Id = missing_value_int4
   integer(kind=int4):: Surface_Elev_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_0_66_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_0_86_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_1_24_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_1_64_Id = missing_value_int4
   integer(kind=int4):: Modis_Alb_2_13_Id = missing_value_int4
   integer(kind=int4):: Snow_Mask_Id = missing_value_int4
   type(Land_grid_Description) :: Coast_Mask_Str
   type(Land_grid_Description) :: Sfc_Type_Str
   type(Land_grid_Description) :: Land_Mask_Str
   type(Land_grid_Description) :: Volcano_Mask_Str
   type(Land_grid_Description) :: Surface_Elev_Str
   type(Land_grid_Description) :: Modis_Alb_0_66_Str
   type(Land_grid_Description) :: Modis_Alb_0_86_Str
   type(Land_grid_Description) :: Modis_Alb_1_24_Str
   type(Land_grid_Description) :: Modis_Alb_1_64_Str
   type(Land_grid_Description) :: Modis_Alb_2_13_Str
   type(Land_grid_Description) :: Snow_Mask_Str

   
   integer, parameter:: LRC_Meander_Flag = 1
   integer, parameter:: Max_LRC_Distance = 10
   real, parameter:: Min_LRC_Jump = 0.0   !0.5
   real, parameter:: Max_LRC_Jump = 100.0 !10.0
   
   integer, parameter:: Grad_Flag_LRC = -1
   real, parameter:: Min_Bt_11um_LRC = 220.0
   real, parameter:: Max_Bt_11um_LRC = 300.0
     
   ! GOES header structures
   TYPE (area_header_type) :: AREAstr
   TYPE (GVAR_NAV)    :: NAVstr
   
   logical :: dncomp_run
   
   character (len = 30) :: string_30
   character (len = 100) :: string_100
   
   !------------- VIIRS variables --------------
   real(kind=real4), dimension(:,:), pointer :: lunar_ref
   real(kind=real4), dimension(:,:),allocatable :: lunar_placeholder
   
   !--- mapping of modis channels to emissivity data-base (Emiss_Chan_Idx are ABI channels)
                                                            !20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38
   integer, dimension(20:45), parameter:: Emiss_Chan_Idx = (/ 7, 7, 7, 7, 7, 7, 0, 9,10,11,12,14,15,16,16,16,16, 8,13, &
                                                            !39,40,41,42,43,44,45
                                                              0, 0, 0, 7, 14, 0,16/)     !Check this
   integer:: Chan_Idx
   type (conf_main_type) :: conf
   
  
   
 !  interface 
 !  subroutine AWG_CLOUD_DNCOMP_ALGORITHM_IBAND (  iseg_in , infile, path, algorithm_started )  
    !--- input
 !     integer, intent(in),optional :: iseg_in
!		character(len=*) , intent(in) :: infile
!      character(len=*) , intent(in) :: path
!      ! - output 
!      logical , intent(out) :: algorithm_started
!   end subroutine AWG_CLOUD_DNCOMP_ALGORITHM_IBAND
!   end interface
   
   !      interface
   !    function open_land_sfc_hdf(data_dir, filename, grid_str) result(id)  
   !    use LAND_SFC_PROPERTIES, only: &
   !   land_grid_description
      
   !   CHARACTER(len=*), intent(in) :: data_dir, filename
   !   TYPE(land_grid_description), optional, intent(inout) :: grid_str
   !   end function open_land_sfc_hdf
   !   end interface  
   
  
   
   !***********************************************************************
   ! Begin Executable Code
   !***********************************************************************
   call set_config (conf)
   call conf % print_files
   call mesg ( '<----------  Start of CLAVRXORB ----------> $Id$' &
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

   call SETUP_USER_DEFINED_OPTIONS()


   !--- make directory for temporary files created during this run
   call system("mkdir "//trim(Temporary_Data_Dir))

   !*************************************************************************
   ! Marker: Open high spatial resolution ancillary data files
   !*************************************************************************
 
   call OPEN_STATIC_ANCIL_FILES()

   !-----------------------------------------------------------------------
   !--- set up surface radiative properties 
   !-----------------------------------------------------------------------
   call SETUP_UMD_PROPS()
 
   !--------------------------------------------------------------------
   !--- setup clock corrections in memory
   !--------------------------------------------------------------------
   if (nav_opt== 2) then
     call SETUP_CLOCK_CORRECTIONS()
   endif

   !**********************************************************************
   ! Marker: Read file directories in FILELIST
   !**********************************************************************

   !--- print to screen which file list is used
   call mesg ( "CLAVR-x FILE LIST FILE USED: "//trim(File_list)  ) 

   !--- open file containing list of level1b data to process
   File_list_lun = GET_LUN()
   open(unit=File_list_lun,file = trim(File_list),status="old",action="read",iostat=ios)
   if (ios /= 0) then
      write ( string_30, '(i8)') ios
      call mesg ("ERROR: Opening clavrxorb_File_list, iostat = "//trim(string_30) , level = verb_lev % ERROR) 
      stop
   end if

   !--- read directories from clavrxorb_input_Files
   read(unit=File_List_Lun,fmt="(a)") Image%Level1b_Path
   read(unit=File_List_Lun,fmt="(a)") Dir_Level2
  
   !--- reset file counter
   File_Number = 1

   !----------------------------------------------------------------------
   ! Marker: BEGIN LOOP OVER FILES
   !----------------------------------------------------------------------
   
   File_loop: do
 
      !----------------------------------------------------------------------
      ! Marker: READ IN CLAVRXORB_FILE_LIST AND SET FLAGS 
      !----------------------------------------------------------------------
      read(unit=File_list_lun,fmt="(a)",iostat=ios) File_1b_Temp
      if ( file_1b_temp == "") exit
      if (ios /= 0) then
         if (ios /= -1) then
            !-- non eof error
            erstat = 8
            call mesg ( "ERROR: Problem reading orbit names from control file" &
               , level = verb_lev % QUIET) 
            stop 8
         else
            !-- end of orbits
            if (File_Number == 1) then
               print *, EXE_PROMPT, "ERROR: No orbits to process, stopping"
               stop
            end if
            exit
         end if
      end if

 
      !----------------------------------------------------------------------------
      ! Determine time of the start of the processing of this orbit
      !----------------------------------------------------------------------------
      Orbital_Processing_Start_Time_Hours = COMPUTE_TIME_HOURS()

      !--------------------------------------------------------------
      ! Determine if this level-1b file can be opended, if not skip
      !--------------------------------------------------------------

      !**********************************************************************
      ! Marker: Prepare to read Level-1b file
      !**********************************************************************

      !-- see if level-1b file exists
      Level1b_Exists = file_test(trim(Image%Level1b_Path)//trim(File_1b_Temp))
      if (.not. Level1b_Exists ) then
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

      !-------------------------------------------------------
      ! reset record counters
      !-------------------------------------------------------
      File_Number = File_Number + 1

      !-----------------------------------------------------------------------
      ! knowing the file type, determine the expected number of elements and
      ! lines in the file
      ! for AVHRR, determine file type and some record lengths
      ! AVHRR Header is read here
      !-----------------------------------------------------------------------
      call  SET_FILE_DIMENSIONS(Image%Level1b_Full_Name &
               ,AREAstr &
               ,Nrec_Avhrr_Header &
               ,Ierror) 
 
      if (Ierror == sym%YES) then
         print *, EXE_PROMPT, "ERROR: Could not set file dimensions, skipping file "
         cycle file_loop
      end if

      if (Image%Number_Of_Lines <= 0) then
         print*,' File dimensions were not set correctly for this sensor ', sensor%sensor_name
         cycle file_loop    
      end if
  
    
      !-----------------------------------------------------------------------
      !--- set up pixel level arrays (size depends on sensor)
      !-----------------------------------------------------------------------
      !--- determine segment size here
      Line_Idx_Min_Segment = 1
      Line_Idx_Max_Segment = Image % Number_Of_Lines_Per_Segment

      !*************************************************************************
      ! Marker:  READ IN HEADER AND DETERMINE SOME CONSTANTS
      !*************************************************************************

      !----------------------------------------------------------------------
      ! Knowing the sensor, interogate files to start, end date and time
      !----------------------------------------------------------------------
      call SET_DATA_DATE_AND_TIME ( AREAstr)

      !----------------------------------------------------------------------
      ! Output sensor and image parameters to screen
      !----------------------------------------------------------------------
      call OUTPUT_SENSOR_TO_SCREEN()
      call OUTPUT_IMAGE_TO_SCREEN()
      call OUTPUT_PROCESSING_LIMITS_TO_SCREEN()
      
      !------------------------------------------------------------------
      ! Setup Solar-channel RTM terms for this particular sensor
      !------------------------------------------------------------------
      call SETUP_SOLAR_RTM(Sensor%WMO_Id)
      call SETUP_OLR()
      call SETUP_SST()

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
   
      !------------------------------------------------------------------
      ! compute the month and day of month
      !------------------------------------------------------------------

      ! Made ileap consistant with oisst call.
      ileap = int(leap_year_Fct(int(Image%Start_Year, kind=int4)),kind=int1)

      Month = int(COMPUTE_MONTH(int(Image%Start_Doy, kind=int4), int(ileap, kind=int4)),kind=int2)
      Day_Of_Month = int(COMPUTE_DAY(int(Image%Start_Doy, kind=int4), int(ileap, kind=int4)),kind=int2)


      !*************************************************************************
      ! Marker:  READ IN SENSOR-SPECIFIC CONSTANTS
      !*************************************************************************

      !--- read in Instrument Constants from appropriate file
      call READ_INSTR_CONSTANTS()


      !*************************************************************************
      ! Marker:  Open non-static high spatial resolution ancillary data
      !*************************************************************************

      !--------------------------------------------------------------------
      ! read in sst analysis file
      !--------------------------------------------------------------------

      !Set to use default options. This will be turned off if there is NO daily OISST data
      Oisst_File_Name= GET_OISST_MAP_FILENAME(Image%Start_Year,Image%Start_Doy, &
                                                 trim(OiSst_Data_Dir) )

      if (trim(OiSst_File_Name) == "no_file") then
          Use_Sst_Anal = .FALSE.
          call mesg ("WARNING: Could not find daily OISST file", level = verb_lev % WARNING )
      else
          call READ_OISST_ANALYSIS_MAP(OISST_File_Name)
          Use_Sst_Anal = .TRUE.
      end if
      
       !------- store previous day and year to prevent reading of same data for next orbit
      Start_Year_Prev = Image%Start_Year
      Start_Day_Prev = Image%Start_Doy  

      !--------------------------------------------------------------------
      !--- 5 km surface emiss
      !--------------------------------------------------------------------
      if (Use_Seebor == sym%YES) then
         print *, EXE_PROMPT, "Reading in SeeBor Data for month = ", Month
         call OPEN_SEEBOR_EMISS(trim(Ancil_Data_Dir)//"static/sfc_data", Month, Emiss_File_Id)
      end if

      !----------------------------------------------------------------------
      ! Open Modis White Sky Albedo Map appropriate for this day
      !----------------------------------------------------------------------
      if (Modis_Clr_Alb_Flag == sym%YES) then

         call MESG("Opening Modis clear albedo map")
         call OPEN_MODIS_WHITE_SKY_SFC_REFLECTANCE_FILES()

      endif
      !------------------------------------------------------------------
      ! Open Snow mask file
      !------------------------------------------------------------------
      call GET_SNOW_MASK()

      !*************************************************************************
      ! Marker:  READ IN NWP DATA
      !*************************************************************************

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
                            Image%End_Year, Image%End_Doy, Image%End_Time, Gdas_Data_Dir, ierror_Nwp) 
      endif
     
      !--- MERRA
      if (Nwp_Opt == 5) then
         call READ_GFS_DATA(Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time,  &
                            Image%End_Year, Image%End_Doy, Image%End_Time, Merra_Data_Dir, ierror_Nwp)
      endif

      !--- ERA INTERIM ANALYSIS
      if (Nwp_Opt == 6) then
         call READ_GFS_DATA(Nwp_Opt, Image%Start_Year, Image%Start_Doy, Image%Start_Time,  &
                            Image%End_Year, Image%End_Doy, Image%End_Time, Erai_Data_Dir, ierror_Nwp)
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

         call POPULATE_PLANCK_TABLES()

         !--- planck for 11 and 12um sounder ch
         if (trim(Sensor%Sensor_Name) == 'AVHRR-IFF' .or. &
             trim(Sensor%Sensor_Name) == 'VIIRS-IFF') then
            call POPULATE_PLANCK_TABLES_SOUNDER()
         endif
 
         if (Aer_Flag  .and. index(Sensor%Sensor_Name,'AVHRR') > 0) then
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

         if (trim(Sensor%Sensor_Name) == 'SEVIRI' .and. Cloud_Mask_Aux_Flag /= sym%NO_AUX_CLOUD_MASK) then
            call READ_OCA(trim(Image%Level1b_Path),trim(Image%Level1b_Name), &
                 Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment, &
                 Segment_Number,Image%Number_Of_Lines,Image%Number_Of_Lines_Read_This_Segment, &
                 Cost_Aux,Pc_Top1_Aux,Pc_Top2_Aux,Pc_Uncertainty1_Aux,Pc_Uncertainty2_Aux, &
                 Tau_Aux,Cld_Phase_Aux)
         endif

         !------------------------------------------------------------------
         ! Apply spatial limits
         !------------------------------------------------------------------
         call EXPAND_SPACE_MASK_FOR_USER_LIMITS(Space_Mask)
  
         !-------------------------------------------------------------------
         ! Modify Chan_On flags to account for channels read in
         !-------------------------------------------------------------------
         call SET_CHAN_ON_FLAG(Sensor%Chan_On_Flag_Default, Sensor%Chan_On_Flag_Per_Line)
        
         !-------------------------------------------------------------------
         ! Compute Lunar Reflectance
         !-------------------------------------------------------------------
         if ((trim(Sensor%Sensor_Name) == 'VIIRS' .or. trim(Sensor%Sensor_Name) == 'VIIRS-NASA') &
             .and. Sensor%Chan_On_Flag_Default(44)) then

           ! - check the angles if this is a good lunar scene
           ! - lun and solar zenith angle

           Lunar_Ref => Ch (44) % Ref_Lunar_Toa
                  
           call COMPUTE_LUNAR_REFLECTANCE (Ch (44) % Rad_Toa &
                      , Geo % Solzen, Geo%Lunzen &
                      , Image % Start_Year, Month, Day_Of_Month, Image % Start_Time &
                      , Geo % Moon_Phase_Angle  &
                      , Ancil_Data_Dir &
                      , Lunar_Placeholder)
            
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
            print *, EXE_PROMPT, "PROCESS CLAVRX WARNING: no scans read in, exiting segment processing loop"
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

               !--- compute desired nwp parameters 
               call COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS()

            end if
   
            !--- determine a pixel-level mask to exclude data solar contaminiation
            call SET_SOLAR_CONTAMINATION_MASK(Solar_Contamination_Mask)

            !--- determine a pixel-level mask to exclude bad or unprocessible data
            call SET_BAD_PIXEL_MASK(Bad_Pixel_Mask)

            !*******************************************************************
            ! Marker: Interpolate ancillary data to each pixel in this segment
            !*******************************************************************
            Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

            !--- map nwp parameters to pixel projection
            call COMPUTE_PIXEL_NWP_PARAMETERS(Smooth_Nwp_Flag)

            !--- map non-nwp ancillary data to the pixel projection
            call MAP_ANCIL_DATA_TO_PIXEL_GRID()

            !--- post process dark composite if one read in
            if (Read_Dark_Comp == sym%YES .and. Dark_Composite_Name /= "no_file") then
               call POST_PROCESS_GOES_DARK_COMPOSITE(Ref_Ch1_Dark_Composite)
            endif

            !--- check ancillary data and modify Bad_Pixel_Mask accordingly
            call QUALITY_CONTROL_ANCILLARY_DATA(Bad_Pixel_Mask)

            !-----------------------------------------------------------------------
            !--- Marker: Compute some fundamental pixel-level arrays
            !-----------------------------------------------------------------------
            !--- compute some common used pixel arrays 
            call COMPUTE_PIXEL_ARRAYS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            !--- normalize reflectances by the solar zenith angle and sun-earth distance
            call NORMALIZE_REFLECTANCES(Sun_Earth_Distance)
   
            !--- compute the channel 20 pseudo reflectance
            if (Sensor%Chan_On_Flag_Default(20)  .and.  Sensor%Chan_On_Flag_Default(31) ) then
              call CH20_PSEUDO_REFLECTANCE(Solar_Ch20_Nu,Geo%CosSolzen,ch(20)%Rad_Toa,ch(31)%Bt_Toa, &
                                           Sun_Earth_Distance,ch(20)%Ref_Toa,Ems_Ch20)
            endif

            !--- compute pixel level Snow map based on all ancillary data
            if (Nwp_Opt /= 0) then
               call COMPUTE_SNOW_CLASS(Sfc%Snow_NWP, Sfc%Snow_OISST, &
                                       Sfc%Snow_IMS,Sfc%Snow_GLOB, &
                                       Sfc%Land,Sfc%Snow)
            end if

            !--- SST if possible for this sensor (needs snow info for masking)
            call COMPUTE_SST()

            !--- interpolate surface type field to each pixel in segment
            call GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)   

            !--- compute desert mask cloud detection
            Sfc%Desert_Mask =  DESERT_MASK_FOR_CLOUD_DETECTION(ch(20)%Sfc_Emiss, Nav%Lat, Sfc%Snow, Sfc%Sfc_Type)

            !--- compute city mask cloud detection
            if (Sensor%Chan_On_Flag_Default(44) ) then
             Sfc%City_Mask =  CITY_MASK_FOR_CLOUD_DETECTION(ch(44)%Rad_Toa, Sfc%Sfc_Type)
            endif

            !--- update time summation for level-1b processing
            End_Time_Point_Hours = COMPUTE_TIME_HOURS()
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
                  call TEMPORAL_INTERP_TMPSFC_NWP(Scan_Time(Line_Idx_Min_Segment),  &
                                Scan_Time(Line_Idx_Min_Segment+Image%Number_Of_Lines_Read_This_Segment-1))
               end if

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

               !--- compute the channel 20 pseudo reflectance for clear-skies
               if (Sensor%Chan_On_Flag_Default(20)   .and.  Sensor%Chan_On_Flag_Default(31)) then
                  call CH20_PSEUDO_REFLECTANCE(Solar_Ch20_Nu,Geo%CosSolzen,ch(20)%Rad_Toa_Clear,ch(31)%Bt_Toa_Clear,Sun_Earth_Distance, &
                                               ch(20)%Ref_Toa_Clear,Ems_Ch20_Clear_Rtm)
 
               endif

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

            if (Sensor%Chan_On_Flag_Default(20) ) then

               !--- apply median filter to Bt_Ch20
               call COMPUTE_MEDIAN_SEGMENT(ch(20)%Bt_Toa,Bad_Pixel_Mask, One, &
                              1, Image%Number_Of_Elements,Line_Idx_Min_Segment, &
                              Line_Idx_Min_Segment+Image%Number_Of_Lines_Read_This_Segment-One, &
                              Bt_Ch20_Median_3x3,  &
                              Bt_Ch20_Std_Median_3x3)
 
               !--- apply median filter to Ems_Ch20
               call COMPUTE_MEDIAN_SEGMENT(Ems_Ch20,Bad_Pixel_Mask, One, &
                              1, Image%Number_Of_Elements,Line_Idx_Min_Segment, &
                              Line_Idx_Min_Segment+Image%Number_Of_Lines_Read_This_Segment-One, &
                              Ems_Ch20_Median_3x3,  &
                              Ems_Ch20_Std_Median_3x3)
            end if

            !--- populate local radiative center arrays
            if (LRC_Flag == sym%YES) then
               if (Sensor%Chan_On_Flag_Default(31) ) then

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

            !------------------------------------------------------------------------------
            ! compute glint masks for use in cloud detection
            !------------------------------------------------------------------------------

            if (Sensor%Chan_On_Flag_Default(31) ) then

             !--- solar glint mask
             if (Sensor%Chan_On_Flag_Default(1) ) then
               call COMPUTE_GLINT(Geo%Glintzen,ch(1)%Ref_Toa, Ref_Ch1_Std_3x3, Sfc%Glint_Mask)
             endif

             !--- lunar glint mask
             if (Sensor%Chan_On_Flag_Default(44) ) then
               call COMPUTE_GLINT(Geo%Glintzen_Lunar,ch(44)%Ref_Lunar_Toa,  &
                                  Ref_ChDNB_Lunar_Std_3x3, Sfc%Glint_Mask_Lunar)
             endif

            endif

            !*******************************************************************
            ! Marker: Generate pixel-level products
            !*******************************************************************

            !---- pixel level aerosol
            if (index(Sensor%Sensor_Name,'AVHRR') > 0 .and. Aer_Flag ) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               call PIXEL_AER_RET_OCEAN(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(5) =  Segment_Time_Point_Seconds(5) + &
                    60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)
            end if

            !--- only apply cloud mask and type routines if nwp/rtm information available
            if (Cld_Flag  .and. Nwp_Opt > 0) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               !--- simple cloud optical depth
               call COMPUTE_SIMPLE_SOLAR_COD(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)               
               call COMPUTE_SIMPLE_LUNAR_COD(Image%Number_Of_Elements,Image%Number_Of_Lines_Read_This_Segment)               

               call OPAQUE_CLOUD_HEIGHT()
               call CO2IRW_CLOUD_HEIGHT()
               call H2O_CLOUD_HEIGHT()
               call CO2_SLICING_CLOUD_HEIGHT(Image%Number_Of_Elements,Line_Idx_Min_Segment, &
                                             Image%Number_Of_Lines_Read_This_Segment, &
                                             P_Std_Rtm, &
                                             Pc_Co2,Tc_Co2,Zc_Co2)

               !--- cloud mask
               if (Cloud_Mask_Aux_Flag /= sym%USE_AUX_CLOUD_MASK) then
                  if (Cloud_Mask_Bayesian_Flag == sym%YES) then
                     call NB_CLOUD_MASK_BRIDGE(Segment_Number)
                  elseif (Cloud_Mask_Bayesian_Flag == sym%NO) then
                     call BASELINE_CLOUD_MASK_MAIN (Segment_Number)
                  else
                     print *, "Only the Bayesian and Baseline Cloud Masks are available, check selection"
                     stop 
                  end if
               end if

               if (Cloud_Mask_Aux_Flag == sym%USE_AUX_CLOUD_MASK .and. Cloud_Mask_Aux_Read_Flag == sym%YES) then
                  CLDMASK%Cld_Mask = Cld_Mask_Aux
               end if

               !--- Compute Adjacent pixels Cloud Mask
               call ADJACENT_PIXEL_CLOUD_MASK(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

               !--- if dark composite available for GOES-1km data, do this
               if (trim(Sensor%Sensor_Name) == 'GOES-IL-IMAGER' .or. trim(Sensor%Sensor_Name) == 'GOES-MP-IMAGER') then
                 if (Read_Dark_Comp == sym%YES .and. Dark_Composite_Name /= "no_file") then
                   call DARK_COMPOSITE_CLOUD_MASK(CLDMASK%Cld_Mask)
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
            if (Cld_Flag  .and. Nwp_Opt > 0) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               call CTP_MULTILAYER()

               Diag_Pix_Array_1 = Zc_H2O
               Diag_Pix_Array_2 = Zc_CO2IRW
               Diag_Pix_Array_3 = Ctp_Multilayer_Flag 

               !-------------------------------------------------------------------
               ! make co2 slicing height from sounder with using sounder/imager IFF
               !-------------------------------------------------------------------
               if (index(Sensor%Sensor_Name,'IFF') > 0) then

                  if (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then

                      call SOUNDER_EMISSIVITY()
                      call MODIFY_CLOUD_TYPE_WITH_SOUNDER (Cld_Temp_Sounder, Cld_Emiss_Sounder, Cld_Type)
                      call MAKE_CIRRUS_PRIOR_TEMPERATURE(Cld_Temp_Sounder, Cld_Press_Sounder, Cld_Emiss_Sounder,  &
                                                         Tc_Cirrus_Background, Zc_Cirrus_Background) 

                  else

                    call CO2_SLICING_CLOUD_HEIGHT(Image%Number_Of_Elements,Line_Idx_Min_Segment, &
                                     Image%Number_Of_Lines_Read_This_Segment, &
                                     P_Std_Rtm, &
                                     Pc_Co2,Tc_Co2,Zc_Co2)

                    call MODIFY_CLOUD_TYPE_WITH_SOUNDER (Tc_CO2, Ec_CO2, Cld_Type)

                    call MAKE_CIRRUS_PRIOR_TEMPERATURE(Tc_Co2, Pc_Co2, Ec_Co2,  &
                                                       Tc_Cirrus_Background, Zc_Cirrus_Background)

                  endif

               endif

               if (ACHA%Mode == 0) then
                  call MODE_ZERO_CLOUD_HEIGHT(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               endif

               if (ACHA%Mode > 0) then 

                  !--- AWG CLoud Height Algorithm (ACHA) and associated products
                  call AWG_CLOUD_HEIGHT_BRIDGE()
                  call ASOS_BRIDGE()

                  !--- interpolate NWP wind and tpw profiles at cloud-top level
                  call COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND_AND_TPW(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

                  !--- interpolate ACHA cloud heights to flight level altitude.
                  call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment,ACHA%Pc,ACHA%Alt)

                  !--accumulate performance metrics
                  call COMPUTE_ACHA_PERFORMANCE_METRICS(ACHA%Processed_Count,ACHA%Valid_Count,ACHA%Success_Fraction)

                  !-- make CSBT masks (Clear Sky Brightness Temperature)
                  call OPAQUE_TRANSMISSION_HEIGHT()
                  call COMPUTE_CSBT_CLOUD_MASKS()

               end if

               !--- Convective Cloud Probability
               call CONVECTIVE_CLOUD_PROBABILITY(Bad_Pixel_Mask,ch(31)%Bt_TOA,ch(27)%Bt_TOA,Ch(31)%Emiss_Tropo,Tsfc_Nwp_Pix,ACHA%Conv_Cld_Prob)
               call SUPERCOOLED_CLOUD_PROBABILITY(Bad_Pixel_Mask,Cld_Type,ACHA%Tc,ACHA%Supercooled_Cld_Prob)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(8) =  Segment_Time_Point_Seconds(8) + &
                   &   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)
    
               ! - lunar reflectance
               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               
               if ((trim(Sensor%Sensor_Name) == 'VIIRS' .or. trim(Sensor%Sensor_Name) == 'VIIRS-NASA') &
                       .and. Sensor%Chan_On_Flag_Default(44)  .and. Nlcomp_Mode > 0) then
                  if ( count (ch(44)%Ref_Lunar_Toa > 0) > 0 ) then
                     call AWG_CLOUD_DNCOMP_ALGORITHM( Iseg_In = Segment_Number , nlcomp_mode = .true. &
                        , algorithm_started = dncomp_run)
                  end if   
               end if   
          
               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(15) =  Segment_Time_Point_Seconds(15) + &
                     & 60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)  
     
               !--- cloud optical depth and effective radius from vis/nir approach
               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               if (Dcomp_Mode > 0) then
                  
                  call AWG_CLOUD_DNCOMP_ALGORITHM( Iseg_In = Segment_Number , algorithm_started = dncomp_run)
                  
                  if ( dcomp_mode .GE. 8 .and. trim(Sensor%Sensor_Name) == 'VIIRS' .and. Sensor%Chan_On_Flag_Default(39) ) then
                     call AWG_CLOUD_DNCOMP_ALGORITHM_IBAND ( infile = trim(File_1b_Temp), path= trim(dir_level2),algorithm_started = dncomp_run )                        
                    
                  end if 
                  
                  
                  call SET_DCOMP_VERSION()
                  if ( dncomp_run ) then
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
                 call CCL_BRIDGE()

                 !---Calculate flight level altitude of the base for AWC
                 call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment,ACHA%Pc_Base,ACHA%Base_Alt)

               endif
               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(10) =  Segment_Time_Point_Seconds(10) + &
                   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            end if

            !--- Non-cloud Detection
            if (Aer_Flag ) then

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
            

            !---  OLR
            call COMPUTE_OLR()

            !---  SASRAB
            if ( Sasrab_Flag ) call INSOLATION(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            End_Time_Point_Hours = COMPUTE_TIME_HOURS()
            Segment_Time_Point_Seconds(12) =  Segment_Time_Point_Seconds(12) + &
                   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            if (Nwp_Opt > 0) then

               !--- assign clear sky quality flags
               call ASSIGN_CLEAR_SKY_QUALITY_FLAGS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

            end if

            !--- generated cloud masked sst field
            call COMPUTE_MASKED_SST()

            !  endif   !end Skip_Processing_Flag condition

            !*******************************************************************
            ! Marker: Write to output files (pixel-level)
            !*******************************************************************
            Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
            
            
            ! = write to level2 files
            
           if ( Level2_File_Flag  == 1 ) then
               call WRITE_PIXEL_HDF_RECORDS(segment_number )
            end if
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

      end do Segment_loop

      call mesg ( "Finished Processing All Orbital Segments")
      call mesg ( " ")

   
      !*************************************************************************
      !   Marker: Close output pixel-level files
      !*************************************************************************

      !*************************************************************************
      ! Marker: Close non-static ancillary data files
      !*************************************************************************

      if (Emiss_File_Id > 0) call CLOSE_SEEBOR_EMISS(Emiss_File_Id)
      if (Modis_Alb_0_66_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_0_66_Id)
      if (Modis_Alb_0_86_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_0_86_Id)
      if (Modis_Alb_1_24_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_1_24_Id)
      if (Modis_Alb_1_64_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_1_64_Id)
      if (Modis_Alb_2_13_Id > 0) call CLOSE_LAND_SFC_HDF(Modis_Alb_2_13_Id)

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

      !--- remove files in temporary directory and then the directory
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
      call mesg ("Total Time for Processing This Orbit (sec) = ", Orbital_Processing_Time_Seconds,&
                                                           level=verb_lev % MINIMAL)
      !  print *, EXE_PROMPT, "Temp Time for Processing This Orbit (sec) = ", Segment_Time_Point_Seconds_temp

      !--- add processing time to global attributes
  

      !*************************************************************************
      ! Marker: End loop over files
      !*************************************************************************
 
   end do File_Loop

   !*************************************************************************
   ! Marker: Close FILELIST
   !*************************************************************************
   close(unit=File_List_Lun)

   !----------------------------------------------------------------------
   ! DEALLOCATE MEMORY AND END PROGRAM
   !----------------------------------------------------------------------

   !*************************************************************************
   !-- Marker: Close static ancillary data files
   !*************************************************************************

   !--- high resolution hdf ancillary data
   if (Sfc_Type_Id > 0) call CLOSE_LAND_SFC_HDF(Sfc_Type_Id)
   if (Coast_Mask_Id > 0)  call CLOSE_LAND_SFC_HDF(Coast_Mask_Id)
   if (Land_Mask_Id > 0)  call CLOSE_LAND_SFC_HDF(Land_Mask_Id)
   if (Volcano_Mask_Id > 0)  call CLOSE_LAND_SFC_HDF(Volcano_Mask_Id)
   if (Surface_Elev_Id > 0)  call CLOSE_LAND_SFC_HDF(Surface_Elev_Id)
   if (Snow_Mask_Id > 0)  call CLOSE_LAND_SFC_HDF(Snow_Mask_Id)

   !--- remove directory for temporary files
   call system("rmdir "//trim(Temporary_Data_Dir))

   !*************************************************************************
   ! Marker: Final remaining memory deallocation
   !*************************************************************************

   !--- Determine time of the start of all processing
   Total_Processing_End_Time_Hours = COMPUTE_TIME_HOURS()
   Total_Processing_Time_minutes = 60.0*(Total_Processing_End_Time_Hours -  &
                                              Total_Processing_Start_Time_Hours)
   call mesg ("Total Time for All Processing (minutes) = ", Total_Processing_Time_minutes, &
                                           level=verb_lev % MINIMAL)

   !---- print to screen that processing is done
   call mesg (  "<--------- End of CLAVRXORB ---------->",level=verb_lev % MINIMAL)

   stop
!--------- PFAAST error statements
 200   print *, EXE_PROMPT, "error reading PFAAST RTM coefficients, stopping"

!*************************************************************************
! Marker: End of main code
!*************************************************************************

contains

!*************************************************************************
! Marker: Begin of Subroutines
!*************************************************************************

!*************************************************************************
! open modis white sky reflectance files 
!*************************************************************************
   subroutine OPEN_MODIS_WHITE_SKY_SFC_REFLECTANCE_FILES()
   
 
      
    
         !--- determine 16 day period and its string value
         iperiod16 = 16 * ((Image%Start_Doy-1) / 16) + 1 
         write(Day_String,fmt="(i3.3)") iperiod16

         !--- Open Modis white sky 0.66 um albedo
         if (Sensor%Chan_On_Flag_Default(1) ) then
            Modis_White_Sky_0_66_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".0.659_x4.hdf"
            Modis_Alb_0_66_Str%sds_Name = MODIS_ALB_0_66_SDS_NAME
            Modis_Alb_0_66_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_0_66_Name), &
                                        grid_Str=Modis_Alb_0_66_Str)
         end if

         !--- Open Modis white sky 0.86 um albedo
         if (Sensor%Chan_On_Flag_Default(2) ) then
            Modis_White_Sky_0_86_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".0.858_x4.hdf"
            Modis_Alb_0_86_Str%sds_Name = MODIS_ALB_0_86_SDS_NAME
            Modis_Alb_0_86_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_0_86_Name), &
                                        grid_Str=Modis_Alb_0_86_Str)
         end if

         !--- Open Modis white sky 1.24 um albedo
         if (Sensor%Chan_On_Flag_Default(5) ) then
            Modis_White_Sky_1_24_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".1.24_x4.hdf"
            Modis_Alb_1_24_Str%sds_Name = MODIS_ALB_1_24_SDS_NAME
            Modis_Alb_1_24_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_1_24_Name), &
                                        grid_Str=Modis_Alb_1_24_Str)
         end if

         !--- Open Modis white sky 1.66 um albedo
         if (Sensor%Chan_On_Flag_Default(6)) then
            Modis_White_Sky_1_64_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".1.64_x4.hdf"
            Modis_Alb_1_64_Str%sds_Name = MODIS_ALB_1_64_SDS_NAME
            Modis_Alb_1_64_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_1_64_Name), &
                                        grid_Str=Modis_Alb_1_64_Str)
         end if                                    

         !--- Open Modis white sky 2.13 um albedo
         if (Sensor%Chan_On_Flag_Default(7) ) then
            Modis_White_Sky_2_13_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".2.13_x4.hdf"
         Modis_Alb_2_13_Str%sds_Name = MODIS_ALB_2_13_SDS_NAME
         Modis_Alb_2_13_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        trim(Modis_White_Sky_2_13_Name), &
                                        grid_Str=Modis_Alb_2_13_Str)
         end if                                    

  end subroutine OPEN_MODIS_WHITE_SKY_SFC_REFLECTANCE_FILES
  !****************************************************************************
  ! acquire the information for the snow classication from multiple sources
  !****************************************************************************
  subroutine GET_SNOW_MASK()

      if (Read_Snow_Mask == sym%READ_SNOW_HIRES) then
         Failed_IMS_Snow_Mask_Flag = sym%NO

         Snow_Mask_Str%sds_Name = SNOW_MASK_SDS_NAME
         Snow_Mask_File_Name = GET_SNOW_MAP_FILENAME(Image%Start_Year,Image%Start_Doy, &
                                                 trim(Snow_Data_Dir))
         if (trim(Snow_Mask_File_Name) == "no_file") then
            Failed_IMS_Snow_Mask_Flag = sym%YES
            call mesg ( "WARNING: Could not find Snow mask file ==> "// &
              Snow_Mask_File_Name, level = verb_lev % WARNING)
         else
            Snow_Mask_Id = OPEN_LAND_SFC_HDF(trim(Snow_Data_Dir), &
                                      Snow_Mask_File_Name, &
                                      grid_Str=Snow_Mask_Str)
            if (Snow_Mask_Id > 0) then
               Failed_IMS_Snow_Mask_Flag = sym%NO
               print *, EXE_PROMPT, "Snow mask file opened successfully "
            else
               Failed_IMS_Snow_Mask_Flag = sym%YES
               print *, EXE_PROMPT, "WARNING: Snow mask file open failed "
            end if
         end if
      end if


      if (Read_Snow_Mask == sym%READ_SNOW_GLOB) THEN 
         Failed_Glob_Snow_Mask_Flag = sym%NO

         Snow_Mask_File_Name = GET_GLOBSNOW_FILENAME(Image%Start_Year,Image%Start_Doy, &
                                                 trim(GlobSnow_Data_Dir))
         if (trim(Snow_Mask_File_Name) == "no_file") THEN
            Failed_Glob_Snow_Mask_Flag = sym%YES
            call mesg ( "WARNING:  Could not find GlobSnow mask file, using NWP ", level = verb_lev % WARNING )
         else
            CALL READ_GLOBSNOW_ANALYSIS_MAP(trim(GlobSnow_Data_Dir)//Snow_Mask_File_Name)
         end if
      endif

  end subroutine GET_SNOW_MASK
  !****************************************************************************
  ! open the hdf files that hold ancillary data
  !****************************************************************************
  subroutine OPEN_STATIC_ANCIL_FILES()
   !------------------------------------------------------------------------------
   !--- Read elevation data 
   !------------------------------------------------------------------------------
   if (Read_Surface_Elevation /= 0) then
      call mesg  ( "Opening surface elevation file", level = verb_lev % VERBOSE)
      Surface_Elev_Str%sds_Name = SURFACE_ELEV_SDS_NAME
     
      !--- read in which elevation type that is specified.
      if (Read_Surface_Elevation == 1) then
         Surface_Elev_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"static/sfc_data/", &
                                      "GLOBE_1km_digelev.hdf", &
                                      grid_Str=Surface_Elev_Str)
      else ! low resolution, Read_Surface_Elevation = 2
         Surface_Elev_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"static/sfc_data/", &
                                     "GLOBE_8km_digelev.hdf", &
                                      grid_Str=Surface_Elev_Str)
      end if

   endif

   !------------------------------------------------------------------
   ! Open coast mask file
   !------------------------------------------------------------------
   if (Read_Coast_Mask == sym%YES) then
      call mesg ( "Opening coast file", level = verb_lev % VERBOSE)
      Coast_Mask_Str%sds_Name = COAST_MASK_SDS_NAME
      Coast_Mask_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"static/sfc_data/", &
                                      "coast_mask_1km.hdf", &
                                      grid_Str=Coast_Mask_Str)
   end if

   !------------------------------------------------------------------
   ! Open land surface type file
   !------------------------------------------------------------------
   call mesg ( "Opening land surface type file", level = verb_lev % VERBOSE)
   Sfc_Type_Str%sds_Name = SFC_TYPE_SDS_NAME

   if (Read_Hires_sfc_type == sym%YES) then
      Sfc_Type_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                     "gl-latlong-1km-landcover.hdf", &
                                      grid_Str=Sfc_Type_Str)
   else
      Sfc_Type_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                     "gl-latlong-8km-landcover.hdf", &
                                      grid_Str=Sfc_Type_Str)
   end if

   !------------------------------------------------------------------
   ! Open land mask file
   !------------------------------------------------------------------
   if (Read_Land_Mask == sym%YES) then
      call mesg  ( "Opening land mask file" ,level = verb_lev % VERBOSE )
      Land_Mask_Str%sds_Name = LAND_MASK_SDS_NAME
      Land_Mask_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                    "lw_geo_2001001_v03m.hdf", &
                                     grid_Str=Land_Mask_Str)
   end if
 
   !------------------------------------------------------------------
   ! Open volcano mask file
   !------------------------------------------------------------------
   if (Read_Volcano_Mask == sym%YES) then
      call mesg  ( "Opening volcano mask file",level = verb_lev % VERBOSE )
      Volcano_Mask_Str%sds_Name = VOLCANO_MASK_SDS_NAME
      Volcano_Mask_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/static/sfc_data/", &
                                        "volcano_mask_1km.hdf", &
                                        grid_Str=Volcano_Mask_Str)
   end if

  end subroutine OPEN_STATIC_ANCIL_FILES
  !*******************************************************************************************************
  ! Read in various ancillary data and map it to the pixel grid for each segment 
  !*******************************************************************************************************
  subroutine  MAP_ANCIL_DATA_TO_PIXEL_GRID()

       !--- surface emissivity
       if (Use_Seebor == sym%YES) then

           !--- force channel 20 read used for desert definition
           call READ_SEEBOR_EMISS(Emiss_File_Id, Emiss_Chan_Idx(20), Nav%Lat, Nav%Lon, Space_Mask, ch(20)%Sfc_Emiss)

           do Chan_Idx = 21, Nchan_Clavrx
               if (ch(Chan_Idx)%Obs_Type /= MIXED_OBS_TYPE .and. &
                   ch(Chan_Idx)%Obs_Type /= THERMAL_OBS_TYPE) cycle

               if (Sensor%Chan_On_Flag_Default(Chan_Idx)) then
                 call READ_SEEBOR_EMISS(Emiss_File_Id, Emiss_Chan_Idx(Chan_Idx), Nav%Lat, Nav%Lon, Space_Mask, Ch(Chan_Idx)%Sfc_Emiss)
               endif

           enddo

        end if

        !--- mandatory fields - check for substitution of Bad_Pixel for space 
 
        !--- surface type
        call READ_LAND_SFC_HDF(Sfc_Type_Id, Sfc_Type_Str, Nav%Lat, Nav%Lon, Space_Mask, Sfc%Sfc_Type)

        !--- surface elevation
        if (Read_Surface_Elevation /= 0) then

               !--- read the high res data
               call READ_LAND_SFC_HDF(Surface_Elev_Id, Surface_Elev_Str, Nav%Lat, Nav%Lon, Space_Mask, Two_Byte_Temp)

               !---  convert to a real number
               Sfc%Zsfc_Hires = real(two_byte_temp,kind=real4)
               !--- values over water are missing, set to zero
               where(Sfc%Zsfc_Hires == Missing_Value_Real4)
                  Sfc%Zsfc_Hires = 0.0
               end where

        end if

        !--- merge with nwp surface elevation
        if (Nwp_Opt /= 0) then
                call MERGE_NWP_HIRES_ZSFC(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
        endif

        !--- read coast mask
        if (Read_Coast_Mask == sym%YES) then
             call READ_LAND_SFC_HDF(Coast_Mask_Id, Coast_Mask_Str, Nav%Lat, Nav%Lon, Space_Mask, Sfc%Coast)
        end if

        !--- read land mask
        if (Read_Land_Mask == sym%YES) then
            call READ_LAND_SFC_HDF(Land_Mask_Id, Land_Mask_Str, Nav%Lat, Nav%Lon, Space_Mask, Sfc%Land)
        end if

        !--- modify land class with ndvi if available (helps mitigate navigation errors)
        call MODIFY_LAND_CLASS_WITH_NDVI(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)

        !--- read volcano mask
        if (Read_Volcano_Mask == sym%YES) then
            call READ_LAND_SFC_HDF(Volcano_Mask_Id, Volcano_Mask_Str, Nav%Lat, Nav%Lon, Space_Mask, Sfc%Volcano_Mask)
        end if

        !--- read Snow mask
        if (Read_Snow_Mask == sym%READ_SNOW_HIRES .and. Failed_IMS_Snow_Mask_Flag == sym%NO) then
             call READ_LAND_SFC_HDF(Snow_Mask_Id, Snow_Mask_Str, Nav%Lat, Nav%Lon, Space_Mask, Sfc%Snow_IMS)
        end if
   
        if (Read_Snow_Mask == sym%READ_SNOW_GLOB .and. Failed_Glob_Snow_Mask_Flag == sym%NO ) then
            call GET_PIXEL_GLOBSNOW_ANALYSIS(Nav%Lat,Nav%Lon,Sfc%Land,Bad_Pixel_Mask,Sfc%Snow_Glob)
        end if
   
        !--- define binary land and coast masks (yes/no) from land and coast flags
        call COMPUTE_BINARY_LAND_COAST_MASKS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)   

        !--- interpolate sst analyses to each pixel
        if (Use_Sst_Anal) then
               call GET_PIXEL_SST_ANALYSIS(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
               call COMPUTE_SNOW_CLASS_OISST(SST_Anal_Cice,Sfc%Snow_OISST)
        end if

        !--- compute a coast mask relative to nwp data
        if (Nwp_Opt /= 0) then
            call COMPUTE_COAST_MASK_NWP(Line_Idx_Min_Segment,Image%Number_Of_Lines_Read_This_Segment)
        end if

        !--- compute a snow classification from NWP
        if (Nwp_Opt /= 0) then
                call COMPUTE_SNOW_CLASS_NWP(Weasd_NWP_Pix, Sea_Ice_Frac_NWP_Pix,Sfc%Snow_NWP)
        end if

        !---- ensure missing values for space scenes
        where (Space_Mask == sym%YES) 
               Sfc%Zsfc_Hires = Missing_Value_Real4
               Sfc%Coast = Missing_Value_Int1
               Sfc%Land = Missing_Value_Int1
               Sfc%Snow_IMS = Missing_Value_Int1
               Sfc%Snow_GLOB = Missing_Value_Int1
               Sfc%Snow_NWP = Missing_Value_Int1
               Sfc%Snow_OISST = Missing_Value_Int1
               Sfc%Volcano_Mask = Missing_Value_Int1
               Sfc%Sfc_Type = Missing_Value_Int1
        end where


        !--- interpolate white sky albedoes to each pixel in segment
        if (Modis_Clr_Alb_Flag == sym%YES) then

           if (Sensor%Chan_On_Flag_Default(1) ) then
                if (.not. allocated(ch(1)%Sfc_Ref_White_Sky)) then
                     allocate(ch(1)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(1)%Sfc_Ref_White_Sky = Missing_Value_Real4
                endif
                call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_0_66_Id, Modis_Alb_0_66_Str, &
                                        ch(1)%Sfc_Ref_White_Sky)
           endif
           if (Sensor%Chan_On_Flag_Default(2) ) then
                  if (.not. allocated(ch(2)%Sfc_Ref_White_Sky)) then
                     allocate(ch(2)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(2)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_0_86_Id, Modis_Alb_0_86_Str, &
                                          ch(2)%Sfc_Ref_White_Sky)
           endif
           if (Sensor%Chan_On_Flag_Default(5) ) then
                  if (.not. allocated(ch(5)%Sfc_Ref_White_Sky)) then
                     allocate(ch(5)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(5)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_1_24_Id, Modis_Alb_1_24_Str, &
                                          ch(5)%Sfc_Ref_White_Sky)
           endif
           if (Sensor%Chan_On_Flag_Default(6) ) then
                  if (.not. allocated(ch(6)%Sfc_Ref_White_Sky)) then
                     allocate(ch(6)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(6)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_1_64_Id, Modis_Alb_1_64_Str, &
                                          ch(6)%Sfc_Ref_White_Sky)
           endif
           if (Sensor%Chan_On_Flag_Default(7) ) then
                  if (.not. allocated(ch(7)%Sfc_Ref_White_Sky)) then
                     allocate(ch(7)%Sfc_Ref_White_Sky(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
                     ch(7)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_2_13_Id, Modis_Alb_2_13_Str, &
                                          ch(7)%Sfc_Ref_White_Sky)
           endif

        endif

    end subroutine  MAP_ANCIL_DATA_TO_PIXEL_GRID

!******************************************************************************
! End of Code
!******************************************************************************

 end program PROCESS_CLAVRX
