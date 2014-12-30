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
   use OISST_ANALYSIS
   use SURFACE_PROPERTIES
   use CLOUD_HEIGHT_ROUTINES
   use ACHA_CLAVRX_BRIDGE
   use CLOUD_BASE_CLAVRX_BRIDGE
   use DCOMP_CLAVRX_BRIDGE_MOD
   use NLCOMP_BRIDGE_MOD
   use AEROSOL_PROPERTIES
   use HDF_PARAMS
   use LAND_SFC_PROPERTIES
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
   use SFC_EMISS
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
   use USER_OPTIONS
   use CLAVRX_MESSAGE_MODULE, only: &
      mesg &
      , verb_lev
   use CLOUD_TYPE_BRIDGE_MODULE
   use SIMPLE_COD
 
   use dnb_retrievals_mod, only: &
      COMPUTE_LUNAR_REFLECTANCE
      
   implicit none
 
   !***********************************************************************
   ! Marker: DECLARE VARIABLES
   !***********************************************************************
   INTEGER(kind=int4), parameter:: nmask_clavr = 5
   INTEGER(kind=int4):: Nword_Clavr
   INTEGER(kind=int4):: Nword_Clavr_Start
   INTEGER(kind=int4):: Nrec_Avhrr_Header
   INTEGER(kind=int4):: Ifile    
   INTEGER(kind=int4):: File_List_Lun                 !logical unit number for File_list
   INTEGER(kind=int4):: Segment_Number
   INTEGER(kind=int4):: Skip_Processing_Flag
   INTEGER(kind=int4):: Lat_Idx
   INTEGER(kind=int4):: Lon_Idx
   INTEGER(kind=int4):: Zen_Idx
   INTEGER(kind=int4):: Elem_Idx  !generic pixel (along scan) index
   INTEGER(kind=int4):: Line_Idx  !generic line (across scan) index
   INTEGER(kind=int4):: Phase_Called_Flag
   logical:: Level1b_Exists
   INTEGER, parameter:: num_Segment_Time_points=15
   REAL(kind=REAL4), dimension(num_Segment_Time_points):: Segment_Time_Point_Seconds
   REAL(kind=REAL4) :: Start_Time_Point_Hours
   REAL(kind=REAL4) :: End_Time_Point_Hours

   REAL(kind=REAL4) :: Segment_Time_Point_Seconds_temp
   REAL(kind=REAL4) :: Start_Time_Point_Hours_temp
   REAL(kind=REAL4) :: End_Time_Point_Hours_temp

   REAL(kind=REAL4) :: Total_Processing_Start_Time_Hours
   REAL(kind=REAL4) :: Total_Processing_End_Time_Hours
   REAL(kind=REAL4) :: Total_Processing_Time_seconds
   REAL(kind=REAL4) :: Orbital_Processing_Start_Time_Hours
   REAL(kind=REAL4) :: Orbital_Processing_End_Time_Hours
   REAL(kind=REAL4) :: Orbital_Processing_Time_Seconds
   character(len=180):: File_1b_Temp
   character(len=180):: File_1b_Full
   INTEGER(kind=int4):: erstat
   REAL(kind=REAL4):: Time_Since_Launch
   INTEGER(kind=int4):: err_reposnx_Flag
   INTEGER(kind=int4), parameter:: One = 1
 
   INTEGER(kind=int4):: ios
   INTEGER(kind=int4):: File_Number
   INTEGER(kind=int4):: Ierror_Level1b
   INTEGER(kind=int4):: ierror_Nwp
   INTEGER(kind=int4):: iperiod16   
   INTEGER(kind=int4) :: ierror
   character(len=3):: Day_String
   character(len=100):: Modis_White_Sky_0_66_Name
   character(len=100):: Modis_White_Sky_0_86_Name
   character(len=100):: Modis_White_Sky_1_24_Name
   character(len=100):: Modis_White_Sky_1_64_Name
   character(len=100):: Modis_White_Sky_2_13_Name
   character(len=100):: Snow_Mask_File_Name
   character(len=100):: oiSst_File_Name
   character(*), parameter :: PROGRAM_NAME = 'CLAVRXORB'

   INTEGER(kind=int4):: Emiss_File_Id = missing_value_int4
   INTEGER(kind=int4):: Coast_Mask_Id = missing_value_int4
   INTEGER(kind=int4):: Land_Mask_Id = missing_value_int4
   INTEGER(kind=int4):: Sfc_Type_Id = missing_value_int4
   INTEGER(kind=int4):: Volcano_Mask_Id = missing_value_int4
   INTEGER(kind=int4):: Surface_Elev_Id = missing_value_int4
   INTEGER(kind=int4):: Modis_Alb_0_66_Id = missing_value_int4
   INTEGER(kind=int4):: Modis_Alb_0_86_Id = missing_value_int4
   INTEGER(kind=int4):: Modis_Alb_1_24_Id = missing_value_int4
   INTEGER(kind=int4):: Modis_Alb_1_64_Id = missing_value_int4
   INTEGER(kind=int4):: Modis_Alb_2_13_Id = missing_value_int4
   INTEGER(kind=int4):: Snow_Mask_Id = missing_value_int4
   TYPE(Land_grid_Description) :: Coast_Mask_Str
   TYPE(Land_grid_Description) :: Sfc_Type_Str
   TYPE(Land_grid_Description) :: Land_Mask_Str
   TYPE(Land_grid_Description) :: Volcano_Mask_Str
   TYPE(Land_grid_Description) :: Surface_Elev_Str
   TYPE(Land_grid_Description) :: Modis_Alb_0_66_Str
   TYPE(Land_grid_Description) :: Modis_Alb_0_86_Str
   TYPE(Land_grid_Description) :: Modis_Alb_1_24_Str
   TYPE(Land_grid_Description) :: Modis_Alb_1_64_Str
   TYPE(Land_grid_Description) :: Modis_Alb_2_13_Str
   TYPE(Land_grid_Description) :: Snow_Mask_Str

   
   INTEGER, PARAMETER:: LRC_Meander_Flag = 1
   INTEGER, PARAMETER:: Max_LRC_Distance = 10
   REAL, PARAMETER:: Min_LRC_Jump = 0.0   !0.5
   REAL, PARAMETER:: Max_LRC_Jump = 100.0 !10.0
   INTEGER, PARAMETER:: Missing_LRC_Value = -999
   INTEGER, PARAMETER:: Grad_Flag_LRC = -1
   REAL, PARAMETER:: Min_Bt_11um_LRC = 220.0
   REAL, PARAMETER:: Max_Bt_11um_LRC = 300.0
     
   ! GOES header structures
   TYPE (AREA_STRUCT) :: AREAstr
   TYPE (GVAR_NAV)    :: NAVstr
   
   logical :: dcomp_run
   
   character ( len = 30) :: string_30
   character ( len = 100) :: string_100
   
   
  
   !------------- VIIRS variables --------------
   real(kind=real4),    dimension(:,:), pointer :: lunar_ref
   real(kind=real4), dimension(:,:),allocatable :: lunar_placeholder
   
   
   
   !--- mapping of modis channels to emissivity data-base (Emiss_Chan_Idx are ABI channels)
                                                            !20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36
   integer, dimension(20:36), parameter:: Emiss_Chan_Idx = (/ 7, 7, 7, 7, 7, 7, 0,10,10,11,12,14,15,16,16,16,16/)
   integer:: Chan_Idx
   
   

   !***********************************************************************
   ! Begin Executable Code
   !***********************************************************************
   call mesg ( '<----------  Start of CLAVRXORB ----------> $Id$' , level = verb_lev % MINIMAL , color = 43 )

#ifdef HDF5LIBS
   call mesg ( "HDF5 USED" , level = verb_lev % VERBOSE )
#endif
   
   !----------------------------------------------------------------------------
   ! Initialize some flags
   !----------------------------------------------------------------------------
   Number_Of_Temporary_Files = 0
   Avhrr_Flag = sym%YES
   Goes_Flag = sym%NO
   Seviri_Flag = sym%NO
   Mtsat_Flag = sym%NO
   FY2_Flag = sym%NO
   COMS_Flag = sym%NO
   Modis_Flag = sym%NO
   Modis_5km_Flag = sym%NO
   Modis_1km_Flag = sym%NO
   Viirs_Flag = sym%NO
   Iff_Viirs_Flag = sym%NO
   Iff_Modis_Flag = sym%NO
  
 
   Skip_Processing_Flag = sym%NO

   !----------------------------------------------------------------------------
   ! Determine time of the start of all processing
   !----------------------------------------------------------------------------
   Total_Processing_Start_Time_Hours = COMPUTE_TIME_HOURS()

   !------------------------------------------------------------------------------
   ! initialize previous date variables
   !-------------------------------------------------------------------------------
   Month_Prev = 0
   Sc_Id_WMO_Prev = 0
   Start_Year_Prev = 0
   Start_Day_Prev = 0

   !*************************************************************************
   ! Marker: Read and Quality Check User Defined Options
   !*************************************************************************
  
   call SETUP_USER_DEFINED_OPTIONS()

   !*************************************************************************
   ! Marker: Open hgh spatial resolution ancillary data files
   !*************************************************************************
 
   !------------------------------------------------------------------------------
   !--- Read elevation data 
   !------------------------------------------------------------------------------
   if (Read_Surface_Elevation /= 0) then
      call mesg  ( "Opening surface elevation file", level = verb_lev % VERBOSE)
      Surface_Elev_Str%sds_Name = SURFACE_ELEV_SDS_NAME
     
      !--- read in which elevation type that is specified.
      if (Read_Surface_Elevation == 1) then
         Surface_Elev_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"sfc_data/", &
                                      "GLOBE_1km_digelev.hdf", &
                                      grid_Str=Surface_Elev_Str)
      else ! low resolution, Read_Surface_Elevation = 2
         Surface_Elev_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"sfc_data/", &
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
      Coast_Mask_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"sfc_data/", &
                                      "coast_mask_1km.hdf", &
                                      grid_Str=Coast_Mask_Str)
   end if

   !------------------------------------------------------------------
   ! Open land surface type file
   !------------------------------------------------------------------
   call mesg ( "Opening land surface type file", level = verb_lev % VERBOSE)
   Sfc_Type_Str%sds_Name = SFC_TYPE_SDS_NAME

   if (Read_Hires_sfc_type == sym%YES) then
      Sfc_Type_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/sfc_data/", &
                                     "gl-latlong-1km-landcover.hdf", &
                                      grid_Str=Sfc_Type_Str)
   else
      Sfc_Type_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/sfc_data/", &
                                     "gl-latlong-8km-landcover.hdf", &
                                      grid_Str=Sfc_Type_Str)
   end if

   !------------------------------------------------------------------
   ! Open land mask file
   !------------------------------------------------------------------
   if (Read_Land_Mask == sym%YES) then
      call mesg  ( "Opening land mask file" ,level = verb_lev % VERBOSE )
      Land_Mask_Str%sds_Name = LAND_MASK_SDS_NAME
      Land_Mask_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/sfc_data/", &
                                    "lw_geo_2001001_v03m.hdf", &
                                     grid_Str=Land_Mask_Str)
   end if
 
   !------------------------------------------------------------------
   ! Open volcano mask file
   !------------------------------------------------------------------
   if (Read_Volcano_Mask == sym%YES) then
      call mesg  ( "Opening volcano mask file",level = verb_lev % VERBOSE )
      Volcano_Mask_Str%sds_Name = VOLCANO_MASK_SDS_NAME
      Volcano_Mask_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/sfc_data/", &
                                        "volcano_mask_1km.hdf", &
                                        grid_Str=Volcano_Mask_Str)
    end if

    !-----------------------------------------------------------------------
    !--- set up surface radiative properties 
    !-----------------------------------------------------------------------
    call SETUP_UMD_PROPS()
 
    !--------------------------------------------------------------------
    !--- setup clock corrections in memory
    !--------------------------------------------------------------------
    if (Nav_Flag == 2) then
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
   read(unit=File_List_Lun,fmt="(a)") Dir_1b
   read(unit=File_List_Lun,fmt="(a)") Dir_1bx
   read(unit=File_List_Lun,fmt="(a)") Dir_Nav_in
   read(unit=File_List_Lun,fmt="(a)") Dir_Nav_Out
   read(unit=File_List_Lun,fmt="(a)") Dir_Cmr
   read(unit=File_List_Lun,fmt="(a)") Dir_Sst
   read(unit=File_List_Lun,fmt="(a)") Dir_Cld
   read(unit=File_List_Lun,fmt="(a)") Dir_Obs
   read(unit=File_List_Lun,fmt="(a)") Dir_Geo
   read(unit=File_List_Lun,fmt="(a)") Dir_Rtm
   read(unit=File_List_Lun,fmt="(a)") Dir_Ash
   read(unit=File_List_Lun,fmt="(a)") Dir_Level2
   read(unit=File_List_Lun,fmt="(a)") Dir_Level3

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
      Level1b_Exists = File_exists(trim(Dir_1b)//trim(File_1b_Temp))
      if (Level1b_Exists .eqv. .FALSE.) then
         print *, EXE_PROMPT, "ERROR: Level-1b file not found, skipping"
         cycle file_loop
      end if

      !--------------------------------------------------------------
      !  Determine based on the L1b file's extension if the input file
      !  is gzipped. Announce this fact thru a boolean var.
      !--------------------------------------------------------------
      call DETERMINE_LEVEL1B_COMPRESSION(File_1b_Temp,L1b_Gzip,L1b_Bzip2,File_1b_Full)

      !------------------------------------------------------------------------
      ! Determine from which sensor this file comes from (MODIS,AVHRR or VIIRS)
      !------------------------------------------------------------------------
      call DETECT_SENSOR_FROM_FILE(File_1b_Full,File_1b_Temp,AREAstr,NAVstr,Ierror)

      if (Ierror == sym%YES) then
         print *, EXE_PROMPT, "ERROR: Sensor could not be detected, skipping file "
         cycle file_loop
      end if

      !--- print to screen the file name
      call mesg (" " )
      call mesg ("<------------- Next Orbit ---------------> ",level = verb_lev % DEFAULT )
      call mesg ("file name = "//trim(File_1b) , level = verb_lev % MINIMAL )

      !-------------------------------------------------------
      ! reset record counters
      !-------------------------------------------------------
      File_Number = File_Number + 1


      !-----------------------------------------------------------------------
      ! knowing the file type, determine the expected number of elements and
      ! lines in the file
      ! for AVHRR, determine file type and some record lengths
      !-----------------------------------------------------------------------
      call  SET_FILE_DIMENSIONS(File_1b_Full,AREAstr,Nrec_Avhrr_Header,Nword_Clavr,Nword_Clavr_Start,Ierror) 
 
      if (Ierror == sym%YES) then
         print *, EXE_PROMPT, "ERROR: Could not set file dimentions, skipping file "
         cycle file_loop
      end if

      if (Num_Scans <= 0) then
         cycle file_loop    
      end if
  
      !-----------------------------------------------------------------------
      !--- Compute the time stamp for use in all generated HDF output files
      !-----------------------------------------------------------------------
      call HDF_TSTAMP()
   
      !-----------------------------------------------------------------------
      !--- set up pixel level arrays (size depends on sensor)
      !-----------------------------------------------------------------------
      !--- determine segment size here
      Line_Idx_Min_Segment = 1
      Line_Idx_Max_Segment = Num_Scans_Per_Segment

      !*************************************************************************
      ! Marker:  READ IN HEADER AND DETERMINE SOME CONSTANTS
      !*************************************************************************

      !----------------------------------------------------------------------
      ! Knowing the sensor, setup internal parameters needed for processing
      !----------------------------------------------------------------------
      call SET_SENSOR_CONSTANTS ( AREAstr)

      !------------------------------------------------------------------
      ! Turn off selected channels based on sensor
      !------------------------------------------------------------------
      call TURN_OFF_CHANNELS_BASED_ON_SENSOR(Avhrr_Flag,Avhrr_1_Flag, &
                                          Goes_Flag, Goes_Mop_Flag, &
                                          Goes_Sndr_Flag,Seviri_Flag,  &
                                          Mtsat_Flag, Viirs_Flag,  &
                                          Iff_Viirs_Flag, Iff_Avhrr_flag, &
                                          FY2_Flag, COMS_Flag)
   
      !------------------------------------------------------------------
      ! Setup PFAAST (FAST IR RTM) for this particular sensor
      !------------------------------------------------------------------
      call SETUP_PFAAST(Sc_Id_WMO)

      !------------------------------------------------------------------
      ! Setup Solar-channel RTM terms for this particular sensor
      !------------------------------------------------------------------
      call SETUP_SOLAR_RTM(Sc_Id_WMO)

      !------------------------------------------------------------------
      ! Create pixel arrays which data for this segment
      !------------------------------------------------------------------
      call CREATE_PIXEL_ARRAYS()

      !------------------------------------------------------------------
      ! Check to see if the channels available and selected allow
      ! for the generation of the algorithms
      !------------------------------------------------------------------
      call CHECK_ALGORITHM_CHOICES()

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
      ileap = leap_year_Fct(int(Start_Year, kind=int4))

      Month = COMPUTE_MONTH(int(Start_Day, kind=int4), int(ileap, kind=int4))
      Day_Of_Month = COMPUTE_DAY(int(Start_Day, kind=int4), int(ileap, kind=int4))

      !--- continue output to screen
      write ( string_30, '(I4,X,I3,X,F9.5)') Start_Year,Start_Day, &
             Start_Time/60.0/60.0/1000.0
      call mesg ("start year, day, time = "//string_30 ) 
      write ( string_30, '(I4,X,I3,X,F9.5)') End_Year,End_Day, &
             End_Time/60.0/60.0/1000.0
      call mesg ("End year, day, time = "//string_30 ) 
      
      

      !*************************************************************************
      ! Marker:  READ IN SENSOR-SPECIFIC CONSTANTS
      !*************************************************************************

      !--- read in Instrument Constants from appropriate file
      call READ_INSTR_CONSTANTS()

      !--- read in Algorithm Constants from appropriate file
      call READ_ALGO_CONSTANTS()

      !*************************************************************************
      ! Marker:  Open non-static high spatial resolution ancillary data
      !*************************************************************************

      !--------------------------------------------------------------------
      ! read in sst analysis file
      !--------------------------------------------------------------------

      !Set to use default options. This will be turned off if there is NO daily OISST data
      Use_Sst_Anal = Use_Sst_Anal_Default

      if (Use_Sst_Anal == sym%YES ) then  
         if (Sst_Anal_Opt < 2) then
            if(Start_Year /= Start_Year_Prev .or. Start_Day /= Start_Day_Prev) then

               OiSst_File_Name= GET_OISST_MAP_FILENAME(Start_Year,Start_Day, &
                                                 trim(OiSst_Data_Dir), &
                                                 Sst_Anal_Opt)

               if (trim(OiSst_File_Name) == "no_file") then
                  ! insert temp sst flag off here
                  Use_Sst_Anal = sym%NO
                  call mesg ("WARNING: Could not find daily OISST file", level = verb_lev % WARNING )
               else
                  call READ_OISST_ANALYSIS_MAP(oiSst_File_Name)
               end if
            endif
         else
            print *, EXE_PROMPT, "WARNING: NESDIS 100 km SST Option not yet available"
         end if
      end if

      !--------------------------------------------------------------------
      !--- 5 km surface emiss
      !--------------------------------------------------------------------
      if (use_seebor == sym%YES) then
         print *, EXE_PROMPT, "Reading in SeeBor Data for month = ", month
         call open_seebor_emiss(trim(Ancil_Data_Dir)//"sfc_data", month, Emiss_File_Id)
      end if

      !----------------------------------------------------------------------
      ! Read in Modis White Sky Albedo Map appropriate for this day
      !----------------------------------------------------------------------
      if (Modis_Clr_Alb_Flag == sym%YES) then

         call mesg ("Opening Modis clear albedo map")

         !--- determine 16 day period and its string value
         iperiod16 = 16 * ((Start_Day-1) / 16) + 1 
         write(Day_String,fmt="(i3.3)") iperiod16

         !--- Open Modis white sky 0.66 um albedo
         if (Chan_On_Flag_Default(1) == sym%YES) then
            Modis_White_Sky_0_66_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".0.659_x4.hdf"
            Modis_Alb_0_66_Str%sds_Name = MODIS_ALB_0_66_SDS_NAME
            Modis_Alb_0_66_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/sfc_data/", &
                                        trim(Modis_White_Sky_0_66_Name), &
                                        grid_Str=Modis_Alb_0_66_Str)
         end if

         !--- Open Modis white sky 0.86 um albedo
         if (Chan_On_Flag_Default(2) == sym%YES) then
            Modis_White_Sky_0_86_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".0.858_x4.hdf"
            Modis_Alb_0_86_Str%sds_Name = MODIS_ALB_0_86_SDS_NAME
            Modis_Alb_0_86_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/sfc_data/", &
                                        trim(Modis_White_Sky_0_86_Name), &
                                        grid_Str=Modis_Alb_0_86_Str)
         end if

         !--- Open Modis white sky 1.24 um albedo
         if (Chan_On_Flag_Default(5) == sym%YES) then
            Modis_White_Sky_1_24_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".1.24_x4.hdf"
            Modis_Alb_1_24_Str%sds_Name = MODIS_ALB_1_24_SDS_NAME
            Modis_Alb_1_24_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/sfc_data/", &
                                        trim(Modis_White_Sky_1_24_Name), &
                                        grid_Str=Modis_Alb_1_24_Str)
         end if

         !--- Open Modis white sky 1.66 um albedo
         if (Chan_On_Flag_Default(6) == sym%YES) then
            Modis_White_Sky_1_64_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".1.64_x4.hdf"
            Modis_Alb_1_64_Str%sds_Name = MODIS_ALB_1_64_SDS_NAME
            Modis_Alb_1_64_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/sfc_data/", &
                                        trim(Modis_White_Sky_1_64_Name), &
                                        grid_Str=Modis_Alb_1_64_Str)
         end if                                    

         !--- Open Modis white sky 2.13 um albedo
         if (Chan_On_Flag_Default(7) == sym%YES) then
            Modis_White_Sky_2_13_Name = "AlbMap.WS.c004.v2.0.00-04."// &
                                 Day_String//".2.13_x4.hdf"
         Modis_Alb_2_13_Str%sds_Name = MODIS_ALB_2_13_SDS_NAME
         Modis_Alb_2_13_Id = OPEN_LAND_SFC_HDF(trim(Ancil_Data_Dir)//"/sfc_data/", &
                                        trim(Modis_White_Sky_2_13_Name), &
                                        grid_Str=Modis_Alb_2_13_Str)
         end if                                    

      end if

      !------------------------------------------------------------------
      ! Open Snow mask file
      !------------------------------------------------------------------
      if (Read_Snow_Mask == sym%READ_SNOW_HIRES) then

         Failed_Hires_Snow_Mask_Flag = sym%NO

         Snow_Mask_Str%sds_Name = SNOW_MASK_SDS_NAME
         Snow_Mask_File_Name = GET_SNOW_MAP_FILENAME(Start_Year,Start_Day, &
                                                 trim(Snow_Data_Dir))
         if (trim(Snow_Mask_File_Name) == "no_file") then
            Failed_Hires_Snow_Mask_Flag = sym%YES
            call mesg ( "WARNING: Could not find Snow mask file ==> "// &
              Snow_Mask_File_Name, level = verb_lev % WARNING)
         else
            Snow_Mask_Id = OPEN_LAND_SFC_HDF(trim(Snow_Data_Dir), &
                                      Snow_Mask_File_Name, &
                                      grid_Str=Snow_Mask_Str)
            if (Snow_Mask_Id > 0) then
               Failed_Hires_Snow_Mask_Flag = sym%NO
               print *, EXE_PROMPT, "Snow mask file opened successfully "
            else
               Failed_Hires_Snow_Mask_Flag = sym%YES
               print *, EXE_PROMPT, "WARNING: Snow mask file open failed "
            end if
         end if
      end if


      IF (Read_Snow_Mask == sym%READ_SNOW_GLOB) THEN 
         Failed_Glob_Snow_Mask_Flag = sym%NO

         Snow_Mask_File_Name = GET_GLOBSNOW_FILENAME(Start_Year,Start_Day, &
                                                 trim(GlobSnow_Data_Dir))
         if (trim(Snow_Mask_File_Name) == "no_file") THEN
            Failed_Glob_Snow_Mask_Flag = sym%YES
            call mesg ( "WARNING:  Could not find GlobSnow mask file, using NWP ", level = verb_lev % WARNING )
         else
            CALL READ_GLOBSNOW_ANALYSIS_MAP(trim(GlobSnow_Data_Dir)//Snow_Mask_File_Name)
         end if

      END IF
    
      !*************************************************************************
      ! Marker:  READ IN NWP DATA
      !*************************************************************************

      !--- GFS
      if (Nwp_Flag == 1) then
         call READ_GFS_DATA(Nwp_Flag, Start_Year, Start_Day, Start_Time, End_Year, End_Day, End_Time, gfs_Data_Dir, ierror_Nwp) 
      end if

      !--- NCEP Reanalysis
      if (Nwp_Flag == 2) then
         call READ_NCEP_REANALYSIS_DATA(Start_Year, Start_Day, Start_Time, End_Day, End_Time, ncep_Data_Dir)
      end if

      !--- CFSR
      if (Nwp_Flag == 3) then
         call READ_GFS_DATA(Nwp_Flag, Start_Year, Start_Day, Start_Time, End_Year, End_Day, End_Time, cfsr_Data_Dir, ierror_Nwp) 
      endif

      !--- GDAS
      if (Nwp_Flag == 4) then
         call READ_GFS_DATA(Nwp_Flag, Start_Year, Start_Day, Start_Time, End_Year, End_Day, End_Time, cfsr_Data_Dir, ierror_Nwp) 
      endif
     
      !---- if NWP is being read in, then proceeed in allocating RTM, NWP arrays
      if (Nwp_Flag /= 0) then
         
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

      !------- store previous day and year to prevent reading of same data for next orbit
      Start_Year_Prev = Start_Year
      Start_Day_Prev = Start_Day

      !*************************************************************************
      ! Marker: Populate or read in other lookup tables
      !*************************************************************************

      !--- planck and aerosol tables
      if (Sc_Id_WMO /= Sc_Id_WMO_Prev) then

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
 
         if (Aer_Flag == sym%YES .and. Avhrr_Flag == sym%YES) then
            call READ_AER_CH123A_REF_LUTS(Ancil_Data_Dir,Sc_Id_WMO)
         end if

         Sc_Id_WMO_Prev = Sc_Id_WMO
 
      end if
 
      !--------------------------------------------------------------
      !--- name and open and write header to lbx file
      !--------------------------------------------------------------
      if (bx_File_Flag == sym%YES) then
         call WRITE_HEADER_1BX(File_1b)
      endif

      !--- compute Sun-Earth distance
      Sun_Earth_Distance = 1.0 - 0.016729*cos(0.9856*(Start_Day-4.0)*dtor)
 
      !--- compute time since launch needed for reflectance calibration
      Time_Since_Launch = Start_Year + Start_Day / 365.25  - launch_Date
 
      !--------------------------------------------------------------
      !-- Marker: Interpolate a clock correction for this orbit
      !--------------------------------------------------------------
      if (Avhrr_Flag == sym%YES .and. Nav_Flag == 2) then
         call INTERPOLATE_CLOCK_ERROR(Start_Year, Start_Time,  &
                                   End_Year, End_Time, &
                                   Sc_Id_WMO,timerr_seconds)
         print *, EXE_PROMPT, "Clock Error = ", timerr_seconds
      end if

      !*************************************************************************
      ! READ IN DATA RECORDS AND PERFORM QUALITY CHECKING
      !*************************************************************************

      !-------------------------------------------------------------------------------------------
      ! Marker: Begin loop over orbit segments
      !-------------------------------------------------------------------------------------------
      call mesg ( " ")
      call mesg ("Started Processing All Orbital Segments")


      !--- compute number of segments in this orbit 
      if (mod(Num_Scans,Num_Scans_Per_Segment) == 0) then
         Num_Segments = Num_Scans / Num_Scans_Per_Segment
      else
         Num_Segments = Num_Scans / Num_Scans_Per_Segment + 1
      end if

      !--- initialize time counters
      Segment_Time_Point_Seconds = 0.0
      Segment_Time_Point_Seconds_temp = 0.0

      Segment_loop: do Segment_Number = 1,Num_Segments
   
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

         call READ_LEVEL1B_DATA(File_1b_Full,Segment_Number,Time_Since_Launch,AREAstr,NAVstr,Nrec_Avhrr_Header,Ierror_Level1b)

         if (Ierror_Level1b /= 0) then
            print *, EXE_PROMPT, "ERROR:  Error reading level1b, skipping this file"
            exit
         end if

         !-------------------------------------------------------------------
         ! Modify Chan_On flags to account for channels read in
         !-------------------------------------------------------------------
         call SET_CHAN_ON_FLAG(Line_Idx_Min_Segment,Num_Scans_Read)
         
         !-------------------------------------------------------------------
         ! Compute Lunar Reflectance
         !-------------------------------------------------------------------
         if (Viirs_Flag == sym%YES .and. Chan_On_Flag_Default(42) == sym%YES) then

           ! - check the angles if this is a good lunar scene
           ! - lun and solar zenith angle

           Lunar_Ref => ch(42)%Ref_Lunar_Toa
                  
           call COMPUTE_LUNAR_REFLECTANCE (ch(42)%Rad_Toa &
                     & , Solzen, Lunzen &
                     & , start_year, month,day_of_month,start_time &
                     & , moon_phase_angle  &
                     & , ancil_data_dir &
                     & , Lunar_placeholder)
            
                  Lunar_Ref = Lunar_placeholder   
        
                  Lunar_Ref=>null() 
                  
                  Lunrelaz = abs ( Lunaz - Sataz )
                  where ( Lunrelaz > 180. ) 
                     Lunrelaz = 360.0 - Lunrelaz  
                  end where
       
            end if
         
         
         End_Time_Point_Hours = COMPUTE_TIME_HOURS()

         !--- update time summation for level-1b processing
         Segment_Time_Point_Seconds(1) =  Segment_Time_Point_Seconds(1) + &
              60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

         !--- check to see that some data was read in
         if (Num_Scans_Read == 0) then
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

               call DEFINE_HDF_FILE_STRUCTURES(Num_Scans, &
                              Dir_Rtm, &
                              Dir_Level2, &
                              File_1b, &
                              Rtm_File_Flag, &
                              Level2_File_Flag, &
                              c1,c2,a1_20,a2_20,nu_20, &
                              a1_31,a2_31,nu_31,a1_32,a2_32,nu_32,Solar_Ch20_Nu,&
                              Sun_Earth_Distance,Therm_Cal_1b, &
                              Ref_Cal_1b,Nav_Flag,Use_Sst_Anal, &
                              Sst_Anal_Opt, Modis_Clr_Alb_Flag,Nwp_Flag, &
                              Ch1_Gain_Low,Ch1_gain_High, &
                              Ch1_Switch_Count_Cal,Ch1_Dark_Count_Cal, &
                              Ch2_Gain_low,Ch2_Gain_High, &
                              Ch2_Switch_Count_Cal,Ch2_Dark_Count_Cal, &
                              Ch3a_Gain_low,Ch3a_Gain_High, &
                              Ch3a_Switch_Count_Cal,Ch3a_Dark_Count_Cal, &
                              Start_Year,End_Year,Start_Day,End_Day,&
                              Start_Time,End_Time)
            end if
   
            !----------------------------------------------------------------------
            ! check to see if this segment fits desired subset if subsetting
            !----------------------------------------------------------------------
            if (Subset_Pixel_Hdf_Flag == sym%YES) then

               !--- check to see if data fell within Solar zenith angle bounds
               if ((minval(Solzen(Num_Pix/2,:)) > Solzen_Max_Limit) .or.   &
                   & (maxval(Solzen(Num_Pix/2,:)) < Solzen_Min_Limit)) then
                  print *, "skipping seg for solzen"
                  cycle   Segment_loop
               end if

               !--- check to see if in latitude range (using diag values)
               if ((minval(Lat_1b(Num_Pix/2,:)) > Lat_Max_Limit) .or.   &
                  & (maxval(Lat_1b(Num_Pix/2,:)) < Lat_Min_Limit)) then
                  print *, "skipping seg for latitude"
                  cycle   Segment_loop 
               end if

            end if

            !*******************************************************************
            ! Marker: Recompute geolocation
            !*******************************************************************
            !--- use level 1b navigation
            if (Nav_Flag == 0 .or. Modis_Flag == sym%YES .or.  &
                &  Seviri_Flag == sym%YES .or. Mtsat_Flag == sym%YES .or. &
                &  Goes_Flag == sym%YES .or. Viirs_Flag == sym%YES .or. &
                &  FY2_Flag == sym%YES .or. COMS_Flag == sym%YES .or. &
                &  Goes_Sndr_Flag == sym%YES ) then
               
               Lat = Lat_1b
               Lon = Lon_1b
            end if

            !---  Repositioning
            if (Nav_Flag == 2 .and. Avhrr_Flag == sym%YES) then
               if ((timerr_seconds /= Missing_Value_Real4) .and. &
                      & (timerr_seconds /= 0.0)) then 
                  call REPOSITION_FOR_CLOCK_ERROR(Line_Idx_Min_Segment,Num_Scans_Read, &
                                  Timerr_Seconds,Err_Reposnx_Flag)
               else
                  Lat = Lat_1b
                  Lon = Lon_1b
               end if
            end if
   
            !--- compute time (local and utc) variables for this segment
            call CONVERT_TIME(Line_Idx_Min_Segment,Num_Scans_Read)

            if (Nwp_Flag /= 0) then

               !--- map each each into correct NWP cell
               call MAP_PIXEL_NWP(Num_Pix,Num_Scans_Read)

               !--- compute needed NWP levels (sfc, tropo, inversion, ...)
               call COMPUTE_NWP_LEVELS_SEGMENT(Num_Pix,Num_Scans_Read)

            end if

   
            !--- determine a pixel-level mask to exclude bad or unprocessible data
            call SET_BAD_PIXEL_MASK(Num_Pix,Num_Scans_Read)

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

               !--- force channel 20 read
               call READ_SEEBOR_EMISS(Emiss_File_Id, Emiss_Chan_Idx(20), lat, lon, Space_Mask, ch(20)%Sfc_Emiss)

               do Chan_Idx = 21, 36
                  if (Chan_Idx == 26) cycle
                  if (Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
                     call READ_SEEBOR_EMISS(Emiss_File_Id, Emiss_Chan_Idx(Chan_Idx), lat, lon, Space_Mask, ch(Chan_Idx)%Sfc_Emiss)
                  endif
               enddo
            end if

            !--- mandatory fields - check for substitution of Bad_Pixel for space 

            !--- surface type
            call READ_LAND_SFC_HDF(Sfc_Type_Id, Sfc_Type_Str, lat, lon, Space_Mask, Sfc_Type)

            !--- surface elevation
            if (Read_Surface_Elevation /= 0) then

               !--- read the high res data
               call READ_LAND_SFC_HDF(Surface_Elev_Id, Surface_Elev_Str, lat, lon, Space_Mask, two_byte_temp)

               !---  convert to a REAL number
               Zsfc_Hires = REAL(two_byte_temp,kind=REAL4)
               !--- values over water are missing, set to zero
               where(Zsfc_Hires == Missing_Value_Real4)
                  Zsfc_Hires = 0.0
               end where

            end if

            !--- merge with nwp surface elevation
            if (Nwp_Flag /= 0) then
                call MERGE_NWP_HIRES_ZSFC(Line_Idx_Min_Segment,Num_Scans_Read)
            endif

            !--- read coast mask
            if (Read_Coast_Mask == sym%YES) then
               call READ_LAND_SFC_HDF(Coast_Mask_Id, Coast_Mask_Str, Lat, Lon, Space_Mask, coast)
            end if

            !--- read land mask
            if (Read_Land_Mask == sym%YES) then
               call READ_LAND_SFC_HDF(Land_Mask_Id, Land_Mask_Str, Lat, Lon, Space_Mask, Land)
            end if

            !--- modify land class with ndvi if available (helps mitigate navigation errors)
            call MODIFY_LAND_CLASS_WITH_NDVI(Line_Idx_Min_Segment,Num_Scans_Read)

            !--- read volcano mask
            if (Read_Volcano_Mask == sym%YES) then
               call READ_LAND_SFC_HDF(Volcano_Mask_Id, Volcano_Mask_Str, lat, lon, Space_Mask, Volcano_Mask)
            end if

            !--- read Snow mask
            if (Read_Snow_Mask == sym%READ_SNOW_HIRES .and. Failed_Hires_Snow_Mask_Flag == sym%NO) then
               call READ_LAND_SFC_HDF(Snow_Mask_Id, Snow_Mask_Str, lat, lon, Space_Mask, Snow_Hires)
            end if
   
            if (Read_Snow_Mask == sym%READ_SNOW_GLOB .and. Failed_Glob_Snow_Mask_Flag == sym%NO ) then
               CALL GET_PIXEL_GLOBSNOW_ANALYSIS(lat,lon,Land,Bad_Pixel_Mask,Snow_Glob)
            end if
   
            !--- define binary land and coast masks (yes/no) from land and coast flags
            call COMPUTE_BINARY_LAND_COAST_MASKS(Line_Idx_Min_Segment,Num_Scans_Read)   

            !---- ensure missing values for space scenes
            where (Space_Mask == sym%YES) 
               Zsfc_Hires = Missing_Value_Real4
               Coast = Missing_Value_Int1
               Land = Missing_Value_Int1
               Snow_Hires = Missing_Value_Int1
               Snow_Glob = Missing_Value_Int1
               Volcano_Mask = Missing_Value_Int1
               Sfc_Type = Missing_Value_Int1
            end where

            !--- interpolate sst analyses to each pixel
            if (Use_Sst_Anal == 1) then
               if (Sst_Anal_Opt < 2) then
                  call GET_PIXEL_SST_ANALYSIS(Line_Idx_Min_Segment,Num_Scans_Read)
               end if
            end if

            !--- compute a coast mask relative to nwp data
            if (Nwp_Flag /= 0) then
               call COMPUTE_COAST_MASK_NWP(Line_Idx_Min_Segment,Num_Scans_Read)
            end if

            !--- interpolate white sky albedoes to each pixel in segment
            if (Modis_Clr_Alb_Flag == sym%YES) then
               if (Chan_On_Flag_Default(1) == sym%YES) then
                    if (.not. allocated(ch(1)%Sfc_Ref_White_Sky)) then
                     allocate(ch(1)%Sfc_Ref_White_Sky(Num_Pix,Num_Scans_Per_Segment))
                     ch(1)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_0_66_Id, Modis_Alb_0_66_Str, &
                                          ch(1)%Sfc_Ref_White_Sky)
               endif
               if (Chan_On_Flag_Default(2) == sym%YES) then
                  if (.not. allocated(ch(2)%Sfc_Ref_White_Sky)) then
                     allocate(ch(2)%Sfc_Ref_White_Sky(Num_Pix,Num_Scans_Per_Segment))
                     ch(2)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_0_86_Id, Modis_Alb_0_86_Str, &
                                          ch(2)%Sfc_Ref_White_Sky)
               endif
               if (Chan_On_Flag_Default(5) == sym%YES) then
                  if (.not. allocated(ch(5)%Sfc_Ref_White_Sky)) then
                     allocate(ch(5)%Sfc_Ref_White_Sky(Num_Pix,Num_Scans_Per_Segment))
                     ch(5)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_1_24_Id, Modis_Alb_1_24_Str, &
                                          ch(5)%Sfc_Ref_White_Sky)
               endif
               if (Chan_On_Flag_Default(6) == sym%YES) then
                  if (.not. allocated(ch(6)%Sfc_Ref_White_Sky)) then
                     allocate(ch(6)%Sfc_Ref_White_Sky(Num_Pix,Num_Scans_Per_Segment))
                     ch(6)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_1_64_Id, Modis_Alb_1_64_Str, &
                                          ch(6)%Sfc_Ref_White_Sky)
               endif
               if (Chan_On_Flag_Default(7) == sym%YES) then
                  if (.not. allocated(ch(7)%Sfc_Ref_White_Sky)) then
                     allocate(ch(7)%Sfc_Ref_White_Sky(Num_Pix,Num_Scans_Per_Segment))
                     ch(7)%Sfc_Ref_White_Sky = Missing_Value_Real4
                  endif
                  call READ_MODIS_WHITE_SKY_ALBEDO(Modis_Alb_2_13_Id, Modis_Alb_2_13_Str, &
                                          ch(7)%Sfc_Ref_White_Sky)
               endif
            endif
 
            !--- post process dark composite if one read in
            if (Read_Dark_Comp == sym%YES .and. Dark_Composite_Name /= "no_file") then
               call POST_PROCESS_GOES_DARK_COMPOSITE(Ref_Ch1_Dark_Composite, Land)
            endif

            !--- check ancillary data
            call QUALITY_CONTROL_ANCILLARY_DATA(Line_Idx_Min_Segment,Num_Scans_Read)

            !-----------------------------------------------------------------------
            !--- Marker: Compute some fundamental pixel-level arrays
            !-----------------------------------------------------------------------

            !--- compute some common used pixel arrays 
            call COMPUTE_PIXEL_ARRAYS(Line_Idx_Min_Segment,Num_Scans_Read)

            if (Viirs_Flag == sym%YES .or. Iff_Viirs_Flag == sym%YES) then
               call COMPUTE_VIIRS_SST(Line_Idx_Min_Segment,Num_Scans_Read)
            end if

            !--- if not viirs, normalize reflectances by the Solar zenith angle
            if (Viirs_Flag == sym%NO .and. Iff_Viirs_Flag == sym%NO .and. Iff_Avhrr_Flag == sym%NO) then
               call NORMALIZE_REFLECTANCES(Sun_Earth_Distance,Line_Idx_Min_Segment,Num_Scans_Read)
            end if
   
            !--- compute the channel 3b albedo arrays
            call CH3B_ALB(Sun_Earth_Distance,Line_Idx_Min_Segment,Num_Scans_Read)

            !--- compute pixel level Snow map based on all ancillary data
            if (Nwp_Flag /= 0) then
               call COMPUTE_SNOW_FIELD(Line_Idx_Min_Segment,Num_Scans_Read)
            end if

            !--- interpolate surface type field to each pixel in segment
            call GET_PIXEL_SFC_EMISS_FROM_SFC_TYPE(Line_Idx_Min_Segment,Num_Scans_Read)   

            !--- compute desert mask cloud detection
            Desert_Mask =  DESERT_MASK_FOR_CLOUD_DETECTION(ch(20)%Sfc_Emiss,Lat, Snow, Sfc_Type)

            !--- compute city mask cloud detection
            if (Chan_On_Flag_Default(42) == sym%YES) then
             City_Mask =  CITY_MASK_FOR_CLOUD_DETECTION(ch(42)%Rad_Toa, Sfc_Type)
            endif

            End_Time_Point_Hours = COMPUTE_TIME_HOURS()
                        

            !--- update time summation for level-1b processing
            Segment_Time_Point_Seconds(2) =  Segment_Time_Point_Seconds(2) + &
                60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            !*******************************************************************
            ! Marker: Compute nwp mapping and rtm values for each pixel in segment
            !*******************************************************************

            !--- needs an NWP to run.
            if (Nwp_Flag /= 0) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               !-- temporally interp skin temp for each segment (only ncep reanalysis)
               if (Nwp_Flag == 2) then
                  call TEMPORAL_INTERP_TMPSFC_NWP(Scan_Time(Line_Idx_Min_Segment), Scan_Time(Line_Idx_Min_Segment+Num_Scans_Read-1))
               end if

               !--- compute desired nwp parameters 
               call COMPUTE_SEGMENT_NWP_CLOUD_PARAMETERS()
               call COMPUTE_PIXEL_NWP_PARAMETERS(Smooth_Nwp_Flag)

               !--- compute a surface temperature from the NWP
               call MODIFY_TSFC_NWP_PIX(1,Num_Pix,Line_Idx_Min_Segment,Num_Scans_Read)

               !--- compute pixel-level rtm parameters 
               Start_Time_Point_Hours_temp = COMPUTE_TIME_HOURS()
               call GET_PIXEL_NWP_RTM(Line_Idx_Min_Segment,Num_Scans_Read)
               End_Time_Point_Hours_temp = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds_temp =  Segment_Time_Point_Seconds_temp + &
                   & 60.0*60.0*(End_Time_Point_Hours_temp - Start_Time_Point_Hours_temp)

               !--- apply atmospheric correction - needs rtm results
               call ATMOS_CORR(Line_Idx_Min_Segment,Num_Scans_Read)

               !--- compute surface products (Tsfc,Ndvi,Rsr ...)
               call SURFACE_REMOTE_SENSING(Line_Idx_Min_Segment,Num_Scans_Read)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(3) =  Segment_Time_Point_Seconds(3) + &
                  &  60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)
            end if

            !*******************************************************************
            ! Marker: Spatial metrics processing
            !*******************************************************************
            Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

            !--- populate needed spatial uniformity arrays
            call COMPUTE_SPATIAL_UNIFORMITY(Line_Idx_Min_Segment,Num_Scans_Read)

            !--- call spatial correlation routines
            call COMPUTE_SPATIAL_CORRELATION_ARRAYS()

            if (Chan_On_Flag_Default(20) == sym%YES) then

               !--- apply median filter to Bt_Ch20
               call COMPUTE_MEDIAN_SEGMENT(ch(20)%Bt_Toa,Bad_Pixel_Mask, One, &
                              1, Num_Pix,Line_Idx_Min_Segment,Line_Idx_Min_Segment+Num_Scans_Read-One, &
                              Bt_Ch20_Median_3x3,  &
                              Bt_Ch20_Std_Median_3x3)
 
               !--- apply median filter to Ems_Ch20
               call COMPUTE_MEDIAN_SEGMENT(Ems_Ch20,Bad_Pixel_Mask, One, &
                              1, Num_Pix,Line_Idx_Min_Segment,Line_Idx_Min_Segment+Num_Scans_Read-One, &
                              Ems_Ch20_Median_3x3,  &
                              Ems_Ch20_Std_Median_3x3)
            end if

            !--- populate local radiative center arrays
            if (LRC_Flag == sym%YES) then
               if (Chan_On_Flag_Default(31) == sym%YES) then

                  call LOCAL_LINEAR_RADIATIVE_CENTER(sym%YES,sym%NO, &
                                       LRC_Meander_Flag, &
                                       ch(31)%Bt_Toa(:,Line_Idx_Min_Segment:Num_Scans_Read-Line_Idx_Min_Segment+1), &
                                       1, &
                                       Num_Pix, &
                                       Line_Idx_Min_Segment,  &
                                       Num_Scans_Read, &
                                       Max_LRC_Distance,  &
                                       Min_LRC_Jump,  &
                                       Max_LRC_Jump,  &
                                       Grad_Flag_LRC,  &
                                       Missing_Value_Int4, &
                                       Bad_Pixel_Mask(:,Line_Idx_Min_Segment:Num_Scans_Read-Line_Idx_Min_Segment+1), &
                                       Min_Bt_11um_LRC,  &
                                       Max_Bt_11um_LRC, &
                                       I_LRC(:,Line_Idx_Min_Segment:Num_Scans_Read-Line_Idx_Min_Segment+1), &
                                       J_LRC(:,Line_Idx_Min_Segment:Num_Scans_Read-Line_Idx_Min_Segment+1))
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
            if (Avhrr_Flag == sym%YES .and. Aer_Flag == sym%YES) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               call PIXEL_AER_RET_OCEAN(Line_Idx_Min_Segment,Num_Scans_Read)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(5) =  Segment_Time_Point_Seconds(5) + &
                    60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)
            end if

            !--- only apply cloud mask and type routines if nwp/rtm information available
            if (Cld_Flag == sym%YES .and. Nwp_Flag > 0) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               !--- simple cloud optical depth
!              call COMPUTE_SIMPLE_COD(Num_Pix,Num_Scans_Read)               

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
               call ADJACENT_PIXEL_CLOUD_MASK(Line_Idx_Min_Segment,Num_Scans_Read)

               !--- if dark composite available for GOES-1km data, do this
               if (GOES_1km_Flag == sym%YES .and. Read_Dark_Comp == sym%YES .and. Dark_Composite_Name /= "no_file") then
                  call DARK_COMPOSITE_CLOUD_MASK(Cld_Mask)
               end if

               !--- accumulate cloud mask metrics
               call COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS(Cloud_Mask_Count,Nonconfident_Cloud_Mask_Count)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(6) =  Segment_Time_Point_Seconds(6) + &
                     & 60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

               !--- cloud type
               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               if (Cloud_Mask_Aux_Flag == sym%USE_AUX_CLOUD_MASK .and. Viirs_Flag == sym%YES) then

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
            if (Cld_Flag == sym%YES .and. Nwp_Flag > 0) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
              

               !-------------------------------------------------------------------
               ! make co2 slicing height from sounder with using sounder/imager IFF
               !-------------------------------------------------------------------
               if (IFF_VIIRS_FLAG == sym%YES .or. &
                   IFF_AVHRR_FLAG == sym%YES .or. &
                   IFF_MODIS_FLAG == sym%YES) then

                   call CO2_SLICING_CLOUD_HEIGHT(Num_Pix,Line_Idx_Min_Segment,Num_Scans_Read, &
                                    P_Std_Rtm,Cld_Type, &
                                    Pc_Cirrus_Co2,Tc_Cirrus_Co2)
               endif

               if (ACHA_Mode == 0) then
                  call MODE_ZERO_CLOUD_HEIGHT(Line_Idx_Min_Segment,Num_Scans_Read)
               endif

               if (ACHA_Mode > 0) then 

                  !--- AWG CLoud Height Algorithm (ACHA) and associated products
                  call AWG_CLOUD_HEIGHT_BRIDGE()

                  !--- interpolate NWP wind profiles at cloud-top level
                  call COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND(Line_Idx_Min_Segment,Num_Scans_Read)

                  !--- interpolate ACHA cloud heights to flight level altitude.
                  call COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min_Segment,Num_Scans_Read)

                  !--accumulate performance metrics
                  call COMPUTE_ACHA_PERFORMANCE_METRICS(Acha_Processed_Count,Acha_Valid_Count)

               end if

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(8) =  Segment_Time_Point_Seconds(8) + &
                   &   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

    
               ! - lunar reflectance
               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               if (viirs_flag == sym%yes .and. chan_on_flag_default(42) == sym % yes) then
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
                     call COMPUTE_CLOUD_WATER_PATH(Line_Idx_Min_Segment,Num_Scans_Read)
                     call COMPUTE_DCOMP_INSOLATION(Line_Idx_Min_Segment,Num_Scans_Read,Sun_Earth_Distance)
                     call  COMPUTE_ADIABATIC_CLOUD_PROPS(Line_Idx_Min_segment,Num_Scans_Read)
                     call COMPUTE_DCOMP_PERFORMANCE_METRICS(DCOMP_Processed_Count,DCOMP_Valid_Count)
                  end if
                  
               end if
    
               !--- compute precipation from optical properties
               if (Dcomp_Mode > 0 ) then
                  call COMPUTE_PRECIPITATION(Line_Idx_Min_Segment,Num_Scans_Read)
               end if

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(9) =  Segment_Time_Point_Seconds(9) + &
                  &  60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

               !--- CLoud Base Height Algorithm if ACHA was executed
               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               if (ACHA_Mode > 0) then 
                 call CLOUD_BASE_BRIDGE()
               endif
               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(10) =  Segment_Time_Point_Seconds(10) + &
                   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            end if

            !--- Non-cloud Detection
            if (Aer_Flag == sym%YES) then

               Start_Time_Point_Hours = COMPUTE_TIME_HOURS()

               !-->        call DUST_DETECTION_ALGORITHM(Line_Idx_Min_Segment,Num_Scans_Read)
               !-->        call SMOKE_DETECTION_ALGORITHM(Line_Idx_Min_Segment,Num_Scans_Read)

               End_Time_Point_Hours = COMPUTE_TIME_HOURS()
               Segment_Time_Point_Seconds(11) =  Segment_Time_Point_Seconds(11) + &
                   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)
            end if

            !--- Volcanic ash - only for AVHRR/2 and AVHRR/3
            if (Ash_Flag > 0) then
               if (avhrr_flag == sym%YES .and. avhrr_1_flag == sym%NO) then

                  start_time_point_hours = COMPUTE_TIME_HOURS()
                  end_time_point_hours = COMPUTE_TIME_HOURS()
                  segment_time_point_seconds(11) =  segment_time_point_seconds(11) + &
                     60.0*60.0*(end_time_point_hours - start_time_point_hours)
               end if

            end if
    
            !--- radiative flux parameters
            Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
            if (erb_Flag == sym%YES) then
               if (AVHRR_1_Flag == sym%NO) then    !currently, no AVHRR/1 algorithm
                  call COMPUTE_ERB(Line_Idx_Min_Segment,Num_Scans_Read)
               end if

               !--- don't run if not a geostationary satellite
               if (GOES_Flag == sym%YES .or. MTSAT_Flag == sym%YES .or. SEVIRI_Flag == sym%YES) then
                  call INSOLATION(Line_Idx_Min_Segment,num_scans_read)
               end if

            end if
            End_Time_Point_Hours = COMPUTE_TIME_HOURS()
            Segment_Time_Point_Seconds(12) =  Segment_Time_Point_Seconds(12) + &
                   60.0*60.0*(End_Time_Point_Hours - Start_Time_Point_Hours)

            if (Nwp_Flag > 0) then

               !--- assign clear sky quality flags
               call ASSIGN_CLEAR_SKY_QUALITY_FLAGS(Line_Idx_Min_Segment,Num_Scans_Read)

            end if

            !--- generated cloud masked sst field
            if (Nwp_Flag > 0) then
                call COMPUTE_MASKED_SST(Line_Idx_Min_Segment,Num_Scans_Read)
            end if

            !  endif   !end Skip_Processing_Flag condition

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
            if (Nwp_Flag /= 0) then

               do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
                  do Elem_Idx = 1, Num_Pix
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
               do Line_Idx = Line_Idx_Min_Segment, Line_Idx_Min_Segment + Num_Scans_Read - 1
                  do Elem_Idx = 1, Num_Pix

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
                        Scan_Number(Line_Idx_Min_Segment), Scan_Number(Num_Scans_Read)
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
      if (Nwp_Flag > 0) then
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
!--------- PFAST error statements
 200   print *, EXE_PROMPT, "error reading PFAST RTM coefficients, stopping"

!*************************************************************************
! Marker: End of code
!*************************************************************************
 end program PROCESS_CLAVRX
