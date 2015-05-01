! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: level2.f90 (src)
!       LEVEL2_ROUTINES (program)
!
! PURPOSE: Routines for creating, writing and closing pixel-level output files
!
! DESCRIPTION: 
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
!          5/09 - Created
!--------------------------------------------------------------------------------------
module LEVEL2_ROUTINES
   use CONSTANTS
   use PIXEL_COMMON
   use HDF
   use SCALING_PARAMETERS
   use HDF_PARAMS
   
   use clavrx_message_module

   implicit none
   private

   public:: DEFINE_HDF_FILE_STRUCTURES, &
            WRITE_PIXEL_HDF_RECORDS, &
            CLOSE_PIXEL_HDF_FILES, &
            WRITE_ALGORITHM_ATTRIBUTES

   private::DEFINE_PIXEL_2D_SDS, &
            DEFINE_PIXEL_3D_SDS

 !--- rtm indices
 integer, parameter, private:: Num_Rtm_Sds = 29
 integer, private, save:: Sd_Id_Rtm
 integer(kind=int4), dimension(Num_Rtm_Sds), save, private:: Sds_Id_Rtm

!----------------------------------------------------------------------
! the following variables taken from process_avhrr_clavr
!----------------------------------------------------------------------
!--- hdf specific variables
 integer(kind=int4),private:: Dim_Id

 !--- compression variables
 integer(kind=int4), dimension(2),private, save::  Comp_Prm
 integer(kind=int4), private, save::  Comp_Type

 integer(kind=int4),parameter,private:: Sds_Rank_1d = 1
 integer(kind=int4),dimension(1),private:: Sds_Dims_1d

 integer(kind=int4),parameter,private:: Sds_Rank_2d = 2
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Dims_2d
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Start_2d
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Stride_2d
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Edge_2d
 integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Chunk_Size_2d

 integer(kind=int4),parameter,private:: Sds_Rank_3d = 3
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Dims_3d
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Start_3d
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Stride_3d
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Edge_3d
 integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Chunk_Size_3d

 character(len=11), private, parameter:: MOD_PROMPT = "LEVEL2:"
 character(len=18), private, parameter:: coordinates_string = "longitude latitude"

 INCLUDE 'level2.inc'

 CONTAINS

!====================================================================
! SUBROUTINE Name: DEFINE_HDF_FILE_STRUCTURES
!
! Function:
!   Determines the structure of the level 2 files.
!
! Description:
!   This subroutine, given the inputs, determins and opens/creates the apppropriate
!       Level 2 (pixel level output) for various types of files. In addition
!       the HDF global attributes are put into the various HDF files 
!
!====================================================================
subroutine DEFINE_HDF_FILE_STRUCTURES(Num_Scans, &
   Dir_Rtm, &
   Dir_Level2, &
   File_1b, &
   Rtm_File_Flag, &
   Level2_File_Flag, &
   c1,c2,a1_20, &
   a2_20, &
   nu_20, &
   a1_31, &
   a2_31, &
   nu_31, &
   a1_32, &
   a2_32, &
   nu_32, &
   Solar_Ch20_nu, &
   Sun_Earth_Distance, &
   Therm_Cal_1b,  &
   Ref_Cal_1b, &
   nav_opt, &
   Use_Sst_anal, &
   Modis_clr_alb_flag, &
   Nwp_Opt, &
   Ch1_gain_low, &
   Ch1_gain_high, &
   Ch1_switch_count, &
   Ch1_dark_count,&
   Ch2_gain_low, &
   Ch2_gain_high, &
   Ch2_switch_count, &
   Ch2_dark_count, &
   Ch3a_Gain_low, &
   Ch3a_Gain_high, &
   Ch3a_Switch_Count, &
   Ch3a_Dark_Count,&
   Start_Year, &
   End_Year, &
   Start_Day, &
   End_Day, &
   Start_Time, &
   End_Time)


 character(len=*), intent(in):: Dir_Rtm
 character(len=*), intent(in):: Dir_Level2
 integer(kind=int4), intent(in):: Num_Scans
 character(len=*), intent(in):: File_1b
 integer, intent(in):: Rtm_File_Flag
 integer, intent(in):: Level2_File_Flag
 
 integer(kind=int4), intent(in) :: Modis_Clr_Alb_Flag
 integer(kind=int4), intent(in) :: Nwp_Opt
 integer, intent(in):: therm_cal_1b
 integer, intent(in):: nav_opt
 integer, intent(in):: use_sst_anal
 integer, intent(in):: Ref_cal_1b
 real(kind=real4), intent(in):: c1,c2,a1_20,a2_20,nu_20,a1_31,a2_31,nu_31,a1_32, &
                                a2_32,nu_32,solar_Ch20_nu,sun_earth_distance, &
                                Ch1_gain_low,Ch1_gain_high,Ch1_switch_count, &
                                Ch1_dark_count, Ch2_gain_low,Ch2_gain_high, &
                                Ch2_switch_count,Ch2_dark_count, Ch3a_gain_low,&
                                Ch3a_gain_high,Ch3a_switch_count, &
                                Ch3a_dark_count
 integer(kind=int4), intent(in):: Start_Time
 integer(kind=int4), intent(in):: End_Time
 integer(kind=int2), intent(in):: Start_Year
 integer(kind=int2), intent(in):: End_Year
 integer(kind=int2), intent(in):: Start_Day
 integer(kind=int2), intent(in):: End_Day
 character(len=4):: l1b_ext
 character(len=128):: File_1b_root
 character(len=128):: File_Rtm
 character(len=128):: File_Level2

 character(len=128):: Long_Name_Temp
 character(len=128):: Standard_Name_Temp
 character(len=128):: Sds_Name

 integer(kind=int4):: blank_int4
 character(len=1):: blank_char
 real:: blank_real4
 real(kind=real4):: Resolution_KM
 integer:: Istatus_Sum
 integer:: Istatus
 integer:: ipos
 integer:: ilen
 
 ! HDF function declarations
 integer:: sfstart
 integer:: sfcreate
 integer:: sfscatt
 integer:: sfsnatt
 integer:: erstat

  !--- begin executable code
  blank_int4 = 0
  blank_real4 = 0.0
  blank_char = " "

  !--- make a file name base for output
  File_1b_Root = trim (file_1b)

  !--- special processing for modis - remove hdf suffix
  l1b_ext = File_1b_Root(len_trim(File_1b_Root)-3:len_trim(File_1b_Root))
  if (trim(l1b_ext) == ".hdf") then
    File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-4)
  endif

  !--- special processing for viirs - remove hdf suffix - this hard coded for
  if (trim(Sensor%Sensor_Name) == 'VIIRS') then
    File_1b_Root = File_1b_Root(7:len_trim(File_1b_Root)-34)
  endif

  !--- special processing for ahi - remove B01.nc suffix - this hard coded for
  if (trim(Sensor%Sensor_Name) == 'AHI') then
    File_1b_Root = File_1b_Root(4:len_trim(File_1b_Root)-12)
  endif

  !--- special processing for IFF - remove hdf suffix - this hard coded for
! !PEATE files
! if (trim(Sensor%Sensor_Name) == 'AQUA-IFF' .or. trim(Sensor%Sensor_Name) == 'AQUA-IFF') then
!   File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-29)
! elseif (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then
!   File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-20)
! endif

  !--- do this for GOES names which are named goesxx_1_year_jday_hour.area
  if (trim(Sensor%Sensor_Name) == 'GOES-IL-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'GOES-MP-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'COMS-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'MTSAT-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'SEVIRI') then

    !-- remove area suffix
    l1b_ext = File_1b_Root(len_trim(File_1b_Root)-3:len_trim(File_1b_Root))
    if (trim(l1b_ext) == "area") then
     File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-5)
    endif
    !-- remove channel number
    if (trim(l1b_ext) == "area") then
        ipos = index(File_1b_Root, "_1_")
        ilen = len(File_1b_Root)
        File_1b_Root = File_1b_Root(1:ipos-1) // "_"//File_1b_Root(ipos+3:ilen)
    endif
  endif

  !--- add 1km tag for 1km GOES files
  if (index(Sensor%Sensor_Name,'GOES') > 0 .and. Sensor%Spatial_Resolution_Meters == 1000) then
    File_1b_Root = trim(File_1b_Root)//".1km"
  endif

  !--- add 'clavrx_' to the file name output
  File_1b_Root = 'clavrx_' // File_1b_Root

  !--- set Resolution_KM for global attribute
  Resolution_KM = Sensor%Spatial_Resolution_Meters / 1000.0

!-------------------------------------------------------------
! define chunking here
!-------------------------------------------------------------
     !--- 3d , first dim is sds dependent
     Sds_Dims_3d(2) = Image%Number_Of_Elements
     Sds_Dims_3d(3) = Num_Scans
     Sds_Chunk_Size_3d(2) = Image%Number_Of_Elements
     Sds_Chunk_Size_3d(3) = Num_Scans

     !-- dimension of 2d variables
     Sds_Dims_2d(1) = Image%Number_Of_Elements
     Sds_Dims_2d(2) = Image%Number_Of_Lines
     Sds_Chunk_Size_2d(1) = Image%Number_Of_Elements
     Sds_Chunk_Size_2d(2) = Image%Number_Of_Lines_Per_Segment

     !-- dimension of 1d variables
     Sds_Dims_1d(1) = Image%Number_Of_Lines

!-------------------------------------------------------------
! define compression here
!-------------------------------------------------------------
 Comp_Type = 0                  !no compression
 Comp_Prm(1) = 0
 Comp_Prm(2) = 0

 if (Compress_Flag == 1) then  !gzip compression
   Comp_Type = 4
   Comp_Prm(1) = 6
   Comp_Prm(2) = 0
 endif

 if (Compress_Flag == 2) then  !szip compression
   Comp_Type = 5
   Comp_Prm(1) = 32
   Comp_Prm(2) = 2
 endif

!---------------------------------------------------------
!-- define rtm file structure
!---------------------------------------------------------
 if (Rtm_File_Flag == sym%YES) then

     file_Rtm = trim(file_1b_root)//".rtm.hdf"
     call mesg ("creating RTM file "//trim(file_Rtm))

     Sd_Id_Rtm = sfstart(trim(dir_Level2)//trim(file_Rtm),DFACC_CREATE)

     if (Sd_Id_Rtm < 0) then
        print *, EXE_PROMPT, MOD_PROMPT, "Creation of RTM product file failed.  Exiting..."
        erstat = 74
        stop 74
     endif

     !--- write global attributes
     call WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES(Sd_Id_Rtm,"PIXEL",file_Rtm,file_1b_root, &
                           resolution_km, &
                           start_year,end_year,start_day,end_day,start_time,end_time,&
                           blank_int4,blank_int4,blank_char,blank_real4, &
                           therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                           modis_clr_alb_flag, nwp_opt, Ch1_gain_low, Ch1_gain_high, &
                           Ch1_switch_count, Ch1_dark_count, &
                           Ch2_gain_low, Ch2_gain_high, &
                           Ch2_switch_count, Ch2_dark_count, &
                           Ch3a_gain_low, Ch3a_gain_high, &
                           Ch3a_switch_count, Ch3a_dark_count, &
                           sun_earth_distance, &
                           c1, c2, a1_20, a2_20, nu_20, &
                           a1_31, a2_31, nu_31, a1_32, a2_32, nu_32, &
                           solar_Ch20_nu,Nav%Timerr_seconds, &
                           acha%mode, dcomp_mode,Sensor%WMO_Id, &
                           Sensor%Platform_Name, Sensor%Sensor_Name, &
                           Dark_Composite_Name,Bayesian_Cloud_Mask_Name)

     !-- reset status flag for error checking
     Istatus_Sum = 0

     !-- scan line
     Sds_Id_Rtm(1) = sfcreate(Sd_Id_Rtm,"scan_line_number",DFNT_INT32,Sds_Rank_1d,Sds_Dims_1d)
     Istatus_Sum = sfsnatt(Sds_Id_Rtm(1), "SCALED", DFNT_INT8, 1, sym%NO_SCALING) + Istatus_Sum
     Istatus_Sum = sfscatt(Sds_Id_Rtm(1), "units", DFNT_CHAR8, 4, "none") + Istatus_Sum
     Istatus_Sum = sfsnatt(Sds_Id_Rtm(1), "RANGE_MISSING", DFNT_FLOAT32, 1, Missing_Value_Real4) + Istatus_Sum

     !--- scan time
     Sds_Id_Rtm(2) = sfcreate(Sd_Id_Rtm,"scan_line_time",DFNT_FLOAT32,Sds_Rank_1d,Sds_Dims_1d)
     Istatus_Sum = sfsnatt(Sds_Id_Rtm(2), "SCALED", DFNT_INT8, 1, sym%NO_SCALING) + Istatus_Sum
     Istatus_Sum = sfscatt(Sds_Id_Rtm(2), "units", DFNT_CHAR8, 5, "hours") + Istatus_Sum
     Istatus_Sum = sfsnatt(Sds_Id_Rtm(2), "RANGE_MISSING", DFNT_FLOAT32, 1, Missing_Value_Real4) + Istatus_Sum

     !--- surface temperature
     call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(3),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "surface_temp_modeled", &
                              "surface_temperature_modeled", &
                              "modeled surface temperature", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Tsfc, Max_Tsfc, "K", Missing_Value_Real4, Istatus)
     Istatus_Sum = Istatus_Sum + Istatus

     !--- Bt_Ch20_Clear_Rtm
     if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(4),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_3_75um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_3_75_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 3.75 microns", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Bt20, Max_Bt20, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch20_clear_solar_Rtm
     if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(5),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_3_75um_nom_clear_sky_solar", &
                               "toa_brightness_temperature_assuming_clear_sky_with_solar_3.75_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "and including solar reflection at the nominal wavelength of 3.75 microns", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt20, Max_Bt20, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch31_Clear_Rtm
     if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(6),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_11_0um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_11_0_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 11.0 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt31, Max_Bt31, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch32_Clear_Rtm
     if (Sensor%Chan_On_Flag_Default(32) == sym%YES) then
        call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(7),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_12_0um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_12_0_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 12.0 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt32, Max_Bt32, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- etrop
     if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(8),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                                "emiss_11um_nom_tropopause", &
                                "emissivity_11_0_micron_nominal_tropopause", &
                                "emissivity at the tropopause modeled at the nominal wavelength of 11.0 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Etropo, Max_Etropo, "none", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ref_Ch1_Clear_Rtm
     if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(9),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "refl_0_65um_nom_clear_sky", &
                               "toa_bidirectional_reflectance_assuming_clear_sky_0_65_micron_nominal", &
                               "top of atmosphere bidirectional reflectance modeled assuming clear skies "// &
                               "at the nominal wavelength of 0.65 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ref_Ch2_Clear_Rtm
     if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(10),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "refl_0_86um_nom_clear_sky", &
                               "toa_bidirectional_reflectance_assuming_clear_sky_0_86_micron_nominal", &
                               "top of atmosphere bidirectional reflectance modeled assuming clear skies "// &
                               "at the nominal wavelength of 0.86 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Ref_Ch2, Max_Ref_Ch2, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch31_std_3x3
     call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(11),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_11_0um_nom_stddev_3x3", &
                              "brightness_temperature_11_0_micron_nominal_stddev_3x3", &
                              "3x3 pixel standard deviation of brightness temperature"// &
                              "at the nominal wavelength of 11.0 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt31_std, Max_Bt31_std, "K", Missing_Value_Real4, Istatus)
     Istatus_Sum = Istatus_Sum + Istatus

     !--- Bt_Ch31_Max_3x3
     if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(12),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                                "temp_11_0um_nom_max_3x3", &
                                "brightness_temperature_11_0_micron_nominal_max_3x3", &
                                "maximum 3x3 pixel brightness temperature"// &
                                "at the nominal wavelength of 11.0 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt31, Max_Bt31, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ref_Ch1_std_3x3
     if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(13),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                                "refl_0_65um_nom_stddev_3x3", &
                                "bidirectional_reflectance_0_65_micron_nom_stddev_3x3", &
                                "3x3 pixel standard deviation of reflectance"// &
                                "at the nominal wavelength of 0.65 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Ref_Ch1_std, Max_Ref_Ch1_std, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ref_Ch1_Min_3x3
     if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(14),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_min_3x3", &
                              "bidirectional_reflectance_0_65_micron_nom_min_3x3", &
                              "3x3 pixel minimum reflectance"// &
                              "at the nominal wavelength of 0.65 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
     endif

     !--- Btd_Ch31_Ch32_Bt_Ch31_Max_3x3
     if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(15),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "diff_ch31_ch32_Bt_ch31_max_3x3", &
                              "difference_11_minus_12_brightness_temperature_max_3x3", &
                              "3x3 pixel maximum of 11.0 micron minus 12.0 micron"// &
                              "brightness temperatures (nominal wavelengths)", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Btd_Ch31_Ch32, Max_Btd_Ch31_Ch32, "K", Missing_Value_Real4, Istatus)
     endif

     !--- Ems_Ch20
     if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(16),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                                "emiss_3_75um_nom", &
                                "top_of_atmosphere_emissivity_3_75_micron_nominal", &
                                "top of atmosphere emissivity at the nominal wavelength of 3.75 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Ems_Ch20, Max_Ems_Ch20, "none", Missing_Value_Real4, Istatus)
     endif

     !--- Ems_Ch20 clear
     if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(17),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                                "emiss_3_75um_nom_clear", &
                                "top_of_atmosphere_emissivity_3_75_micron_nominal_clear", &
                                "top of atmosphere clear sky emissivity estimate "// &
                                "at the nominal wavelength of 3.75 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Ems_Ch20, Max_Ems_Ch20, "none", Missing_Value_Real4, Istatus)
     endif

     !--- Ems_Ch20 median 3x3
     if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(18),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "emiss_3_75um_nom_median_3x3", &
                              "emissivity_3_75_micron_nominal_median_3x3", &
                              "3x3 pixel median top of atmosphere emissivity "// &
                              "at the nominal wavelength of 3.75 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ems_Ch20, Max_Ems_Ch20, "none", Missing_Value_Real4, Istatus)
     endif

     !--- Tc_opaque_cloud
     call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(19),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cld_temp_opaque", &
                              "cloud_top_temperature_opaque", &
                              "cloud top temperature assuming an opaque cloud", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Tc, Max_Tc, "K", Missing_Value_Real4, Istatus)

     !--- Bt_Ch20
     if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(20),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_3_75um_nom_median_3x3", &
                              "brightness_temperature_3.7_micron_nominal_median_3x3", &
                              "3x3 pixel median filtered brightness temperature "// &
                              "at the nominal wavelength of 3.75 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt20, Max_Bt20, "K", Missing_Value_Real4, Istatus)
     endif


     !--- Ch20 surface emissivity
     call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(21),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "surface_emiss_3_75um_nom", &
                              "surface_emissivity_3.7_micron_nominal", &
                              "surface emissivity at the nominal wavelength of 3.75 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_sfc_ems, Max_sfc_ems, "none", Missing_Value_Real4, Istatus)

     !--- oisst_uni
     call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(22),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "sst_background_uni_3x3", &
                              "sea_surface_skin_temperature_background_uni_3x3", &
                              "background sea surface skin temperature 3x3", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Sst_std, Max_Sst_std, "K", Missing_Value_Real4, Istatus)

     !--- Bt_Ch27_Clear_Rtm
     if (Sensor%Chan_On_Flag_Default(27) == sym%YES) then
        call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(23),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_6_7um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_6_7_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 6.7 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt27, Max_Bt27, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch28_Clear_Rtm
     if (Sensor%Chan_On_Flag_Default(28) == sym%YES) then
        call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(24),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_7_3um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_7_3_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 7.3 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt28, Max_Bt28, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch29_Clear_Rtm
     if (Sensor%Chan_On_Flag_Default(29) == sym%YES) then
        call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(25),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_8_5um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_8_5_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 8.5 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt29, Max_Bt29, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch33_Clear_Rtm
     if (Sensor%Chan_On_Flag_Default(33) == sym%YES) then
        call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(26),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_13_3um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_13_3_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 13.3 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt33, Max_Bt33, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Covariance of Ch27 and Ch31
     if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
          call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(27),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_11um_vs_67um_covar_5x5", &
                              "brightness_temperature_11_vs_67_micron_5x5_covariance", &
                              "5x5 pixel covariance of brightness temperatures at "// &
                              "11 microns vs 67 microns (nominal wavelengths)", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt_Covar, Max_Bt_Covar, "K^2", Missing_Value_Real4, Istatus)
     endif

     !--- Zc_opaque_cloud
     call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(28),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cld_height_opaque", &
                              "cloud_top_height_assuming_opaque", &
                              "cloud top height assuming an opaque cloud", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Zc, Max_Zc, "m", Missing_Value_Real4, Istatus)

     !--- Naive Bayesian Surface Type
     call DEFINE_PIXEL_2D_SDS(Sds_Id_Rtm(29),Sd_Id_Rtm,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "bayes_mask_sfc_type", &
                              "bayes_mask_sfc_type", &
                              "bayes_mask_sfc_type", &
                              DFNT_INT8, sym%NO_SCALING, &
                              0.0, 0.0, "none", Missing_Value_Real4, Istatus)

     !--- check for and report errors
     if (Istatus_Sum /= 0) then
       print *, EXE_PROMPT, MOD_PROMPT, "Error defining sds in rtm hdf file"
     endif

  endif


!---------------------------------------------------------
!-- define level2 file structure
!---------------------------------------------------------
  if (Level2_File_Flag == sym%YES) then
     File_Level2= trim(File_1b_Root)//".level2.hdf"
     call mesg (MOD_PROMPT//"creating level-2 file "//trim(file_Level2))
     Sd_Id_Level2 = sfstart(trim(Dir_Level2)//trim(file_Level2),DFACC_CREATE)
     if (Sd_Id_Level2 < 0) then
       erstat = 68
       print *, EXE_PROMPT, MOD_PROMPT, "Level-2 file creation failed. Exiting..."
       stop 68
     endif

    !--- write attributes
    call WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES(Sd_Id_Level2,"PIXEL",File_Level2,File_1b_Root, &
                           Resolution_KM, &
                           start_year,end_year,start_day,end_day,start_time,end_time,&
                           blank_int4,blank_int4,blank_char,blank_real4, &
                           therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                           modis_clr_alb_flag, nwp_opt, Ch1_gain_low, Ch1_gain_high, &
                           Ch1_switch_count, Ch1_dark_count, &
                           Ch2_gain_low, Ch2_gain_high, &
                           Ch2_switch_count, Ch2_dark_count, &
                           Ch3a_gain_low, Ch3a_gain_high, &
                           Ch3a_switch_count, Ch3a_dark_count, &
                           sun_earth_distance, &
                           c1, c2, a1_20, a2_20, nu_20, &
                           a1_31, a2_31, nu_31, a1_32, a2_32, nu_32, &
                           solar_Ch20_nu,Nav%Timerr_seconds, &
                           acha%mode, dcomp_mode,Sensor%WMO_Id, &
                           Sensor%Platform_Name, Sensor%Sensor_Name, &
                           Dark_Composite_Name,Bayesian_Cloud_Mask_Name)

     !-- reset status flag for error checking
     Istatus_Sum = 0

     !-- scan line variables

     !--- scan line number
     if (Sds_Num_Level2_Scanline_Flag == sym%YES) then
      Sds_Id_Level2(Sds_Num_Level2_Scanline) = sfcreate(Sd_Id_Level2,"scan_line_number",DFNT_INT32,Sds_Rank_1d,Sds_Dims_1d)
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Scanline), "SCALED", DFNT_INT8, 1, sym%NO_SCALING) + Istatus_Sum
      Istatus_Sum = sfscatt(Sds_Id_Level2(Sds_Num_Level2_Scanline), "units", DFNT_CHAR8, 4, "none") + Istatus_Sum
      Standard_Name_temp = "scan_line_number"
      Istatus_Sum = sfscatt(Sds_Id_Level2(Sds_Num_Level2_Scanline), "standard_name", DFNT_CHAR8, 13, "not specified") + Istatus_Sum
      Long_Name_temp = "scan line number"
      Istatus_Sum = sfscatt(Sds_Id_Level2(Sds_Num_Level2_Scanline), "long_name", DFNT_CHAR8,  &
                    len_trim(Long_Name_temp), trim(Long_Name_temp)) + Istatus_Sum
                    
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Scanline), "RANGE_MISSING", DFNT_FLOAT32, &
                    1,real(Missing_Value_Int4,kind=real4)) + Istatus_Sum
                   
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Scanline), "_FillValue", DFNT_INT32, &
                    1,Missing_Value_Int4) + Istatus_Sum
     endif

     !--- scanline time
     if (Sds_Num_Level2_Time_Flag == sym%YES) then
      Sds_Id_Level2(Sds_Num_Level2_Time) = sfcreate(Sd_Id_Level2,"scan_line_time",DFNT_FLOAT32,Sds_Rank_1d,Sds_Dims_1d)
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Time), "SCALED", DFNT_INT8, 1, sym%NO_SCALING) + Istatus_Sum
      Istatus_Sum = sfscatt(Sds_Id_Level2(Sds_Num_Level2_Time), "units", DFNT_CHAR8, 5, "hours_utc") + Istatus_Sum
      Istatus_Sum = sfscatt(Sds_Id_Level2(Sds_Num_Level2_Time), "standard_name", DFNT_CHAR8, 13, "not specified") + Istatus_Sum
      Long_Name_temp = "time for the scan line in fractional hours"
      Istatus_Sum = sfscatt(Sds_Id_Level2(Sds_Num_Level2_Time), "long_name", DFNT_CHAR8,  &
                   len_trim(Long_Name_temp), trim(Long_Name_temp)) + Istatus_Sum
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Time), "RANGE_MISSING", DFNT_FLOAT32,  &
                   1, Missing_Value_Real4) + Istatus_Sum
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Time), "_FillValue", DFNT_FLOAT32,  &
                   1, Missing_Value_Real4) + Istatus_Sum
     endif

     !--- bad scan flag
     if (Sds_Num_Level2_Bad_Scan_Flag == sym%YES) then
      Sds_Id_Level2(Sds_Num_Level2_Bad_Scan) = sfcreate(Sd_Id_Level2,"bad_scan_line_flag",DFNT_INT8,Sds_Rank_1d,Sds_Dims_1d)
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Bad_Scan), "SCALED", DFNT_INT8, 1, sym%NO_SCALING) + Istatus_Sum
      Istatus_Sum = sfscatt(Sds_Id_Level2(Sds_Num_Level2_Bad_Scan), "units", DFNT_CHAR8, 4, "none") + Istatus_Sum
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Bad_Scan), "_FillValue", DFNT_INT8,  &
                   1, Missing_Value_Int1) + Istatus_Sum
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Bad_Scan), "RANGE_MISSING", DFNT_FLOAT32,  &
                    1, Real(Missing_Value_Int1,kind=real4)) + Istatus_Sum
     endif

     !--- asc des flag
     if (Sds_Num_Level2_Asc_Flag_Flag == sym%YES) then
      Sds_Id_Level2(Sds_Num_Level2_Asc_Flag) = sfcreate(Sd_Id_Level2,"asc_des_flag",DFNT_INT8,Sds_Rank_1d,Sds_Dims_1d)
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Asc_Flag), "SCALED", DFNT_INT8, 1, sym%NO_SCALING) + Istatus_Sum
      Istatus_Sum = sfscatt(Sds_Id_Level2(Sds_Num_Level2_Asc_Flag), "units", DFNT_CHAR8, 4, "none") + Istatus_Sum
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Asc_Flag), "_FillValue", DFNT_INT8,  &
                   1, Missing_Value_Int1) + Istatus_Sum
      Istatus_Sum = sfsnatt(Sds_Id_Level2(Sds_Num_Level2_Asc_Flag), "RANGE_MISSING", DFNT_FLOAT32,  &
                    1, Real(Missing_Value_Int1,kind=real4)) + Istatus_Sum
     endif

     !--- Bad Pixel Mask
     if (Sds_Num_Level2_Bad_Pixel_Mask_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bad_Pixel_Mask),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "bad_pixel_mask", &
                               "not specified", &
                               "mask that distinguishes good(0) from bad(1) pixels ", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 1.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Gap Pixel Mask
     if (Sds_Num_Level2_Gap_Pixel_Mask_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Gap_Pixel_Mask),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "gap_pixel_mask", &
                               "not specified", &
                               "mask that distinguishes not in gap (0) from in-gap(1) pixels ", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 1.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Diagnostic 1
      if (Sds_Num_Level2_Diag1_Flag == sym%YES) then 
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Diag1),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "diagnostic_1", &
                              "not specified", &
                              "First diagnostic variable (contents will change)", &
                              DFNT_FLOAT32, sym%NO_SCALING, &
                              0.0,0.0, "unknown", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Diagnostic 2
      if (Sds_Num_Level2_Diag2_Flag == sym%YES) then 
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Diag2),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "diagnostic_2", &
                              "not specified", &
                              "Second diagnostic variable (contents will change)", &
                              DFNT_FLOAT32, sym%NO_SCALING, &
                              0.0,0.0, "unknown", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Diagnostic 3
      if (Sds_Num_Level2_Diag3_Flag == sym%YES) then 
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Diag3),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "diagnostic_3", &
                              "not specified", &
                              "third diagnostic variable (contents will change)", &
                              DFNT_FLOAT32, sym%NO_SCALING, &
                              0.0,0.0, "unknown", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  packed pixel quality meta-data
     if (Sds_Num_Level2_Meta_Data_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Meta_Data),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "packed_pixel_meta_data", &
                               "pixel_quality_flags_packed_into_one_byte", &
                               "order_and_depth: bad_pixel_mask(1),solar_contamination_mask(1),"//&
                               "ch6_on_pixel_mask(1),Bayes_Mask_Sfc_Type(3)", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", No_Attribute_Missing_Value, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- latitude
     if (Sds_Num_Level2_Lat_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Lat),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "latitude", &
                               "latitude", &
                               "latitude", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Lat, Max_Lat, "degrees_north", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- longitude
     if (Sds_Num_Level2_Lon_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Lon),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "longitude", &
                              "longitude", &
                              "longitude", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Lon, Max_Lon, "degrees_east", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- parallax corrected latitude
     if (Sds_Num_Level2_Latpc_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Latpc),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "latitude_pc", &
                               "latitude_parallax_corrected", &
                               "latitude_parallax_corrected_using_cloud_height", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Lat, Max_Lat, "degrees_north", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- parallax corrected longitude
     if (Sds_Num_Level2_Lonpc_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Lonpc),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "longitude_pc", &
                              "longitude_parallax_corrected", &
                              "longitude_parallax_corrected_using_cloud_height", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Lon, Max_Lon, "degrees_east", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- sensor zenith
     if (Sds_Num_Level2_Zen_Flag == sym%YES) then
      Sds_Name = "sensor_zenith_angle"
      if (NCDC_Attribute_Flag == 1) Sds_Name = "anchor_sensor_zenith"
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Zen),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              trim(Sds_Name), &
                              "sensor_zenith_angle", &
                              "sensor zenith for each pixel measured in degrees from nadir", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Zen, Max_Zen, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif 

     !--- solar zenith
     if (Sds_Num_Level2_Solzen_Flag == sym%YES) then
      Sds_Name = "solar_zenith_angle"
      if (NCDC_Attribute_Flag == 1) Sds_Name = "anchor_solar_zenith"
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Solzen),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              trim(Sds_Name), &
                              "solar_zenith_angle", &
                              "solar zenith for each pixel measured in degrees away from the sun (0=looking at sun)", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Solzen, Max_Solzen, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- relative azimuth
     if (Sds_Num_Level2_Relaz_Flag == sym%YES) then
      Sds_Name = "relative_azimuth_angle"
      if (NCDC_Attribute_Flag == 1) Sds_Name = "anchor_relative_azimuth"
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Relaz),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               trim(Sds_Name), &
                               "relative_sensor_azimuth_angle", &
                               "relative azimuth angle in degrees. 0 is the principal plane looking towards sun", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Relaz, Max_Relaz, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- solar azimuth
     if (Sds_Num_Level2_Solaz_Flag == sym%YES) then
      Sds_Name = "solar_azimuth_angle"
      if (NCDC_Attribute_Flag == 1) Sds_Name = "anchor_solar_azimuth"
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Solaz),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               trim(Sds_Name), &
                               "solar_azimuth_angle", &
                               "solar azimuth angle in degrees from north, pixel to sun, positive values are clockwise from north", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Solaz, Max_Solaz, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- sensor azimuth
     if (Sds_Num_Level2_Sataz_Flag == sym%YES) then
      Sds_Name = "sensor_azimuth_angle"
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Sataz),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               trim(Sds_Name), &
                               "sensor_azimuth_angle", &
                               "sensor azimuth angle in degrees from north, pixel to sensor, positive values are clockwise from north", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Sataz, Max_Sataz, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- glint zenith
     if (Sds_Num_Level2_Glintzen_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Glintzen),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "glint_zenith_angle", &
                              "glint_zenith_angle", &
                              "glint zenith for each pixel measured in degrees away from the specular image of sun", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Glintzen, Max_Glintzen, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- scattering zenith
     if (Sds_Num_Level2_Scatzen_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Scatzen),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "scattering_angle", &
                              "scattering_angle", &
                              "scattering angle for each pixel measured in degrees away from direction of forward scattering", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Scatang, Max_Scatang, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- lunar zenith
     if (Sds_Num_Level2_Lunzen_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Lunzen),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "lunar_zenith_angle", &
                              "lunar_zenith_angle", &
                              "lunar zenith for each pixel measured in degrees away from the moon (0=looking at sun)", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Solzen, Max_Solzen, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- lunar relative azimuth
     if (Sds_Num_Level2_LunRelaz_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_LunRelaz),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "lunar_relative_azimuth_angle", &
                               "lunar_relative_azimuth_angle", &
                               "relative azimuth angle in degrees. 0 is the principal plane looking towards moon", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Relaz, Max_Relaz, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- lunar azimuth
     if (Sds_Num_Level2_Lunaz_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Lunaz),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "lunar_azimuth_angle",  &
                               "lunar_azimuth_angle", &
                               "lunar azimuth angle in degrees from north", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Solaz, Max_Solaz, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !---  packed land cover
     if (Sds_Num_Level2_Packed_Land_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Packed_Land),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "packed_land_cover", &
                               "packed_land_cover", &
                               "land cover, snow and coast values packed into one byte, see patmos-x docs to unpack", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", No_Attribute_Missing_Value, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- glint mask
     if (Sds_Num_Level2_Glint_Mask_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Glint_Mask),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "glint_mask", &
                               "not specified", &
                               "glint mask (0=no) (1=yes)", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 1.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- coast mask
     if (Sds_Num_Level2_Coast_Mask_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Coast_Mask),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "coast_mask", &
                               "not specified", &
                               "coast mask (0=no) (1=yes)", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 1.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  surface type
     if (Sds_Num_Level2_Sfc_Type_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Sfc_Type),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "surface_type", &
                               "area_type", &
                               "UMD surface type: water=0,evergreen_needle=1,evergreen_broad=2,"// & 
                               "deciduous_needle=3,deciduous_broad=4,mixed_forest=5,woodlands=6,wooded_grass=7"// &
                               "closed_shrubs=8,open_shrubs=9,grasses=10,croplands=11,bare=12,urban=13",  &
                               DFNT_INT8, sym%NO_SCALING, Min_Sfc_Type, Max_Sfc_Type, &
                               "none", REAL(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  land classification
     if (Sds_Num_Level2_Land_Mask_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Land_Mask),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "land_class", &
                               "not specified", &
                               "land classes and values:shallow ocean=0,land=1,coastline=2,shallow inland water=3,"// & 
                               "ephemeral water=4,deep inland water=5,moderate ocean=6,deep ocean=7 ", &
                               DFNT_INT8, sym%NO_SCALING, Min_Land_Class, Max_Land_Class, &
                               "none", REAL(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  snow classification
     if (Sds_Num_Level2_Snow_Mask_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Snow_Mask),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "snow_class", &
                               "not specified", &
                               "snow classes and values:no snow/ice=1,sea_ice=2,snow=3",& 
                               DFNT_INT8, sym%NO_SCALING, Min_Snow_Class, Max_Snow_Class, &
                               "none", REAL(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  surface height
     if (Sds_Num_Level2_Zsfc_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Zsfc),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                              "surface_elevation", &
                              "surface_elevation", &
                              "surface elevation above mean sea level", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Zsfc, Max_Zsfc, "meters", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

!------------------------------------------------------------------------------------------------------------------------
!----- START OF OBSERVATIONS
!------------------------------------------------------------------------------------------------------------------------

     !--- Ch1 reflectance
     if (Sds_Num_Level2_Ch1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
      Sds_Name = "refl_0_65um_nom"
      if (NCDC_Attribute_Flag == 1) Sds_Name = "ch1_reflectance"
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch1),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              trim(Sds_Name), &
                              "toa_bidirectional_reflectance", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.65 microns - NOAA CDR", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch2 reflectance
     if (Sds_Num_Level2_Ch2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(2) == sym%YES) then
       Sds_Name = "refl_0_86um_nom"
       if (NCDC_Attribute_Flag == 1) Sds_Name = "ch2_reflectance"
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch2),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              trim(Sds_Name), &
                              "toa_bidirectional_reflectance", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.86 microns - NOAA CDR", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch2, Max_Ref_Ch2, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch3 reflectance
     if (Sds_Num_Level2_Ch3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(3) == sym%YES) then
       Sds_Name = "refl_0_47um_nom"
       if (NCDC_Attribute_Flag == 1) Sds_Name = "ch3_reflectance"
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch3),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              trim(Sds_Name), &
                              "not specified", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.47 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch3, Max_Ref_Ch3, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch4 reflectance
     if (Sds_Num_Level2_Ch4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(4) == sym%YES) then
       Sds_Name = "refl_0_55um_nom"
       if (NCDC_Attribute_Flag == 1) Sds_Name = "ch4_reflectance"
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch4),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              trim(Sds_Name), &
                              "not specified", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.55 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch4, Max_Ref_Ch4, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch5 reflectance
     if (Sds_Num_Level2_Ch5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(5) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch5),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_1_2um_nom", &
                              "not specified", &
                              "top of atmosphere reflectance at the nominal wavelength of 1.2 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch5, Max_Ref_Ch5, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch6 reflectance
     if (Sds_Num_Level2_Ch6_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) then
       Sds_Name = "refl_1_60um_nom"
       if (NCDC_Attribute_Flag == 1) Sds_Name = "ch3a_reflectance"
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch6),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              trim(Sds_Name), &
                              "toa_bidirectional_reflectance", &
                              "top of atmosphere reflectance at the nominal wavelength of 1.60 microns - NOAA CDR", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch6, Max_Ref_Ch6, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
      endif

     !--- Ch7 reflectance
     if (Sds_Num_Level2_Ch7_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(7) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch7),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_2_1um_nom", &
                              "toa_bidirectional_reflectance_2_1_micron_nominal", &
                              "top of atmosphere reflectance at the nominal wavelength of 2.1 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch7, Max_Ref_Ch7, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch8 reflectance
     if (Sds_Num_Level2_Ch8_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(8) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch8),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_41um_nom", &
                              "toa_bidirectional_reflectance_0_41_micron_nominal", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.41 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch8, Max_Ref_Ch8, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch9 reflectance
     if (Sds_Num_Level2_Ch9_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(9) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch9),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_44um_nom", &
                              "toa_bidirectional_reflectance_0_44_micron_nominal", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.44 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch9, Max_Ref_Ch9, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch17 reflectance
     if (Sds_Num_Level2_Ch17_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(17) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch17),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_95um_nom", &
                              "toa_bidirectional_reflectance_0_95_micron_nominal", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.95 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch17, Max_Ref_Ch17, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch18 reflectance
     if (Sds_Num_Level2_Ch18_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(18) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch18),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_93um_nom", &
                              "toa_bidirectional_reflectance_0_93_micron_nominal", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.93 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch18, Max_Ref_Ch18, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch19 reflectance
     if (Sds_Num_Level2_Ch19_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(19) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch19),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_94um_nom", &
                              "toa_bidirectional_reflectance_0_94_micron_nominal", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.94 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch19, Max_Ref_Ch19, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch26 reflectance
     if (Sds_Num_Level2_Ch26_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(26) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch26),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_1_38um_nom", &
                              "not specified", &
                              "top of atmosphere reflectance at the nominal wavelength of 1.38 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch26, Max_Ref_Ch26, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
      endif

     !--- ChDNB reflectance
     if (Sds_Num_Level2_ChDNB_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then 
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_ChDNB),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_sol_dnb_nom", &
                              "toa_bidirectional_reflectance_solar_dnb_nominal", &
                              "top of atmosphere reflectance solar dnb channel", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_ChDNB, Max_Ref_ChDNB, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
      endif
  
     !--- ChDNB reflectance lunar
     if (Sds_Num_Level2_ChDNB_lunar_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then 
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_ChDNB_lunar),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_lunar_dnb_nom", &
                              "toa_bidirectional_reflectance_lunar_dnb_nominal", &
                              "top of atmosphere reflectance lunar dnb channel", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_ChDNB_lunar, Max_Ref_ChDNB_lunar, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
      endif

     !--- Ch20 reflectance
     if (Sds_Num_Level2_Ch20_Ref_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch20_Ref),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_3_75um_nom", &
                              "toa_bidirectional_reflectance", &
                              "top of atmosphere reflectance at the nominal wavelength of 3.75 microns - NOAA CDR", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch20, Max_Ref_Ch20, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
      endif

     !--- Ch20 temperature
     if (Sds_Num_Level2_Ch20_Bt_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch20_Bt),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_3_75um_nom", &
                              "toa_brightness_temperature", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 3.75 microns - NOAA CDR", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt20, Max_Bt20, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch22 temperature
     if (Sds_Num_Level2_Ch22_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(22) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch22),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_3_9um_nom", &
                              "toa_brightness_temperature_3_9_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 3.9 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt22, Max_Bt22, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch27 temperature
     if (Sds_Num_Level2_Ch27_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch27),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_6_7um_nom", &
                              "toa_brightness_temperature_6_7_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 6.7 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt27, Max_Bt27, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Ch28 temperature
     if (Sds_Num_Level2_Ch28_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(28) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch28),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_7_3um_nom", &
                              "toa_brightness_temperature_7_3_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 7.3 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt28, Max_Bt28, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Ch29 temperature
     if (Sds_Num_Level2_Ch29_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(29) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch29),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_8_5um_nom", &
                              "toa_brightness_temperature_8_5_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 8.5 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt29, Max_Bt29, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Ch30 temperature
     if (Sds_Num_Level2_Ch30_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(30) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch30),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_9_7um_nom", &
                              "toa_brightness_temperature_9_7_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 9.7 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt30, Max_Bt30, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Ch31 temperature
     if (Sds_Num_Level2_Ch31_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch31),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_11_0um_nom", &
                              "toa_brightness_temperature", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 11.0 microns - NOAA CDR", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt31, Max_Bt31, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Ch32 temperature
     if (Sds_Num_Level2_Ch32_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch32),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_12_0um_nom", &
                              "toa_brightness_temperature", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 12.0 microns - NOAA CDR", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt32, Max_Bt32, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Ch33 temperature
     if (Sds_Num_Level2_Ch33_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(33) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch33),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_13_3um_nom", &
                              "toa_brightness_temperature_13_3_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 13.3 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt33, Max_Bt33, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Ch34 temperature
     if (Sds_Num_Level2_Ch34_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(34) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch34),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_13_6um_nom", &
                              "toa_brightness_temperature_13_6_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 13.6 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt34, Max_Bt34, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Ch35 temperature
     if (Sds_Num_Level2_Ch35_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(35) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch35),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_13_9um_nom", &
                              "toa_brightness_temperature_13_9_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 13.9 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt35, Max_Bt35, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Ch36 temperature
     if (Sds_Num_Level2_Ch36_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(36) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch36),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_14_2um_nom", &
                              "toa_brightness_temperature_14_2_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 14.2 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt36, Max_Bt36, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch37 temperature
     if (Sds_Num_Level2_Ch37_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(37) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch37),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_6_2um_nom", &
                              "toa_brightness_temperature_6_2_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 6.2 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt37, Max_Bt37, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch38 temperature
     if (Sds_Num_Level2_Ch38_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch38),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_10_4um_nom", &
                              "toa_brightness_temperature_10_4_micron_nominal", &
                              "top of atmosphere brightness temperature at the nominal wavelength of 10.4 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt38, Max_Bt38, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch1_Ref_std_3x3
     if (Sds_Num_Level2_Ch1_Std_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch1_std),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_stddev_3x3", &
                              "not specified", &
                              "standard deviation of the 0.63 micron reflectance computed over a 3x3 "//&
                              "pixel array", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch1_std, Max_Ref_Ch1_std, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
     !--- Bt_Ch31_std_3x3
     if (Sds_Num_Level2_Ch31_Std_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch31_std),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_11_0um_nom_stddev_3x3", &
                              "not specified", &
                              "standard deviation of the 11 micron brightness temperature "// &
                              "computed over a 3x3 pixel array", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Bt31_std, Max_Bt31_std, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

!---- Cloud Properties

     !--- cloud probability
     if (Sds_Num_Level2_Cldprob_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cldprob),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cloud_probability", &
                              "not specified", &
                              "probability of a pixel being cloudy from the Bayesian cloud mask", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_frac, Max_frac, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud mask
     if (Sds_Num_Level2_Cld_Mask_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cld_Mask),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cloud_mask", &
                               "not specified", &
                               "integer classification of the cloud mask including clear=0, probably-clear=1, "// &
                               "probably-cloudy=2, cloudy=3", &
                               DFNT_INT8, sym%NO_SCALING, Min_Cld_Mask, Max_Cld_Mask, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- adjacent pixel cloud mask
     if (Sds_Num_Level2_Adj_Pix_Cld_Mask_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Adj_Pix_Cld_Mask),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "adj_pix_cloud_mask", &
                               "adj_pix_cloud_mask", &
                               "integer classification of the adjacent pixel cloud mask including clear=0, "// &
                               "probably-clear=1, probably-cloudy=2, cloudy=3", &
                               DFNT_INT8, sym%NO_SCALING, Min_Cld_Mask, Max_Cld_Mask, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---- Packed Cloud Mask Test Vectors (First Byte / ACM Mask Only)
     if (Sds_Num_Level2_Cld_Tests_Flag == sym%YES) then

       Sds_Dims_3d(1) = Max_Num_Cld_Test_Bytes
       Sds_Chunk_Size_3d(1) = Max_Num_Cld_Test_Bytes

        call DEFINE_PIXEL_3D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cld_Tests),Sd_Id_Level2, &
                             Sds_Dims_3d, &
                             Sds_Chunk_Size_3d, &
                             "cloud_mask_test_packed_results", &
                             "cloud_mask_test_packed_results", &
                             "individual cloud mask packed test results", &
                             DFNT_INT8, sym%NO_SCALING, &
                             0.0, 0.0, "none", real(Missing_Value_Int1),  &
                             "byte_number", &
                             "pixel_elements_along_scan_direction", &
                             "scan_lines_along_track_direction", &
                             Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud mask from 1b
     if (Sds_Num_Level2_Cld_Mask_Aux_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cld_Mask_Aux),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cloud_mask_1b", &
                               "cloud_mask_1b", &
                               "integer classification of the cloud mask including clear=0, probably-clear=1, "// &
                               "probably-cloudy=2, cloudy=3 from mask read from level-1b file", &
                               DFNT_INT8, sym%NO_SCALING, Min_Cld_Mask, Max_Cld_Mask, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- bayesian mask surface type
     if (Sds_Num_Level2_Bayes_Sfc_Type_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bayes_Sfc_Type),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "bayes_mask_sfc_type", &
                               "bayes_mask_sfc_type", &
                               "integer classification of the surface type assumed in constructed the Bayesian cloud mask, "// &
                               "1=deep water,2=shallow ocean,3=land,4=snow,5=arctic,6=antarctic+greenland,7=desert", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 7.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

      
     !--- dust mask
     if (Sds_Num_Level2_Dust_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Dust),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "dust_mask", &
                               "not specified", &
                               "integer classification of the presence of dust (0=no,1=yes)", &
                               DFNT_INT8, sym%NO_SCALING, Min_Binary_Mask, Max_Binary_Mask, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- smoke mask
     if (Sds_Num_Level2_Smoke_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Smoke),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "smoke_mask", &
                               "not specified", &
                               "integer classification of the presence of smoke (0=no,1=yes)", &
                               DFNT_INT8, sym%NO_SCALING, Min_Binary_Mask, Max_Binary_Mask, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- shadow mask
     if (Sds_Num_Level2_Shadow_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Shadow),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "shadow_mask", &
                               "not specified", &
                               "integer classification of the presence of shadow (0=no,1=yes)", &
                               DFNT_INT8, sym%NO_SCALING, Min_Binary_Mask, Max_Binary_Mask, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- fire mask
     if (Sds_Num_Level2_Fire_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Fire),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "fire_mask", &
                               "not specified", &
                               "integer classification of the presence of fire (0=no,1=yes)", &
                               DFNT_INT8, sym%NO_SCALING, Min_Binary_Mask, Max_Binary_Mask, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud optical depth for Mask
     if (Cld_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1)==sym%YES .and. Sds_Num_Level2_Cod_Mask_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cod_Mask),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d,&
                               "cld_opd_mask", &
                               "atmosphere_optical_thickness_due_to_cloud_assuming_water_phase", &
                               "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                               "and water phase with 10 micron particle size determined for cloud mask use", &
                                DFNT_INT8, sym%LINEAR_SCALING, Min_Tau, Max_Tau, &
                               "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud type
     if (Sds_Num_Level2_Cld_Type_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cld_Type),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cloud_type", &
                               "not specified", &
                               "integer classification of the cloud type including clear "// &
                               "and aerosol type,0=clear,1=probably clear,2=fog,3=water,4=supercooled water,"//&
                               "5=mixed,6=opaque_ice,7=cirrus,8=overlapping,9=overshooting,10=unknown,11=dust,12=smoke", &
                               DFNT_INT8, sym%NO_SCALING, Min_Cld_Type, Max_Cld_Type, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud phase
     if (Sds_Num_Level2_Cld_Phase_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cld_Phase),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cloud_phase", &
                               "cloud_phase", &
                               "integer classification of the cloud phase including clear "// &
                               "and aerosol type,0=clear,1=water,2=supercooled water,3=mixed,4=ice,5=unknown", &
                               DFNT_INT8, sym%NO_SCALING, Min_Cld_Phase, Max_Cld_Phase, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud type 1b
     if (Sds_Num_Level2_Cld_Type_Aux_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cld_Type_Aux),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cloud_type_1b", &
                               "cloud_type_1b", &
                               "integer classification of the cloud type including clear "// &
                               "and aerosol type read from level-1b file", &
                               DFNT_INT8, sym%NO_SCALING, Min_Cld_Type, Max_Cld_Type, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud pressure from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ctp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cld_press_acha", &
                               "air_pressure_at_cloud_top", &
                               "cloud-top pressure computed using the AWG cloud height algorithm", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Pc, Max_Pc, &
                               "hPa", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud temperature from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctt_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ctt),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_temp_acha", &
                               "air_temperature_at_cloud_top", &
                               "cloud-top temperature computed using the AWG cloud height algorithm - NOAA CDR", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Tc, Max_Tc, &
                               "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud height from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cth),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_height_acha", &
                               "height_at_cloud_top", &
                               "cloud height computed using the AWG cloud height algorithm", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Zc, Max_Zc, &
                               "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud height from acha cloud top
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Top_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cth_Top),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_height_top_acha", &
                               "top_height_of_cloud", &
                               "estimate of actual cloud-top height computed using the AWG cloud height algorithm", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Zc, Max_Zc, &
                               "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud height from acha cloud base
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Base_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cth_Base),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_height_base_acha", &
                               "base_height_of_cloud", &
                               "estimate of actual cloud-base height computed using the AWG cloud height algorithm", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Zc, Max_Zc, &
                               "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud height of lower cloud from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Zc_Lower_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Zc_Lower),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_height_lower_acha", &
                               "height_of_lower_cloud", &
                               "estimate of actual cloud height of lower cloud computed using the AWG cloud height algorithm", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Zc, Max_Zc, &
                               "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud pressure of lower cloud from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Pc_Lower_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Pc_Lower),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_pressure_lower_acha", &
                               "pressure_of_lower_cloud", &
                               "estimate of actual cloud pressure of lower cloud computed using the AWG cloud height algorithm", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Pc, Max_Pc, &
                               "hPa", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud altitude from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Alt_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Alt),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_altitude_acha", &
                               "altitude_at_cloud_top", &
                               "cloud height altitude computed from pressure values", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Alt, Max_Alt, &
                               "feet", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  acha processing order
     if (Sds_Num_Level2_Acha_Order_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Acha_Order),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "acha_processing_order", &
                               "acha_processing_order", &
                               "integer classification of the order of processing with ACHA", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  acha inversion flag
     if (Sds_Num_Level2_Acha_Inver_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Acha_Inver),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "acha_inversion_flag", &
                               "acha_inversion_flag", &
                               "flag stating whether ACHA was processed assuming an inversion(1) or not(0)", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- cloud height from h2o
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_H2O_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cth_H2O),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_height_h2o", &
                               "height_at_cloud_top_from_h2o_intercept", &
                               "cloud-top height computed using the two-point h2o intercept", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Zc, Max_Zc, &
                               "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud height from opaque
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Opa_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cth_Opa),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_height_opaque", &
                               "height_at_cloud_top_assuming_opaque", &
                               "cloud-top height computed using assuming the cloud is opaque", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Zc, Max_Zc, &
                               "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud emissivity from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ec_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ec),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cld_emiss_acha", &
                              "convective_cloud_longwave_emissivity", &
                              "cloud emissivity at the nominal wavelength of 11 microns, "// &
                              "determined from the AWG cloud height algorithm - NOAA CDR", &
                              DFNT_INT8, sym%LINEAR_SCALING, Min_Ec, Max_Ec, &
                             "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- beta from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Beta_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Beta),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_beta_acha", &
                               "not specified", &
                               "cloud 11/12 micron beta value determined from the split-window method", &
                                DFNT_INT8, sym%LINEAR_SCALING, Min_beta, Max_beta, &
                               "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud-top height uncertainty from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Acha_Uncer_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cth_Acha_Uncer),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_height_uncer_acha", &
                               "not specified", &
                               "cloud height uncertainty computed using the AWG cloud height algorithm", &
                               DFNT_INT8, sym%LINEAR_SCALING, Min_Zc_Uncer, Max_Zc_Uncer, &
                               "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud-top temperature uncertainty from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctt_Acha_Uncer_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ctt_Acha_Uncer),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_temp_uncer_acha", &
                               "not specified", &
                               "cloud temperature uncertainty computed using the AWG cloud height algorithm", &
                               DFNT_INT8, sym%LINEAR_SCALING, Min_Tc_Uncer, Max_Tc_Uncer, &
                               "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- quality flag for ACHA Cloud-top Temperature
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Acha_Qf_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cth_Acha_Qf),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cld_temp_acha_qf", &
                               "not specified", &
                               "quality flag for cloud-top temperature from ACHA "// &
                               "not attempted=0, failed=1, low quality=2, high quality=3", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 3.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- quality flag for ACHA Cloud Emissivity
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ec_Acha_Qf_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ec_Acha_Qf),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cld_emiss_acha_qf", &
                               "not specified", &
                               "quality flag for 11.0 micron cloud emissivity from ACHA "// &
                               "not attempted=0, failed=1, low quality=2, high quality=3", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 3.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- quality flag for ACHA Cloud Beta
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Beta_Acha_Qf_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Beta_Acha_Qf),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cld_beta_acha_qf", &
                               "not specified", &
                               "quality flag for cloud 11.0/12.0 micron beta from ACHA "// &
                               "not attempted=0, failed=1, low quality=2, high quality=3", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 3.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud optical depth from ACHA
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Acha_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cod_Acha),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d,&
                               "cld_opd_acha", &
                               "atmosphere_optical_thickness_due_to_cloud", &
                               "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                               "determined from ACHA", &
                                DFNT_INT8, sym%LINEAR_SCALING, Min_Tau_Acha, Max_Tau_Acha, &
                               "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud effective particle radius from ACHA
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Acha_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ceps_Acha),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cld_reff_acha", &
                              "effective_radius_of_cloud_condensed_water_particles_at_cloud_top", &
                              "effective radius of cloud particles determined from ACHA; "//&
                              "see attributes for channels used", &
                              DFNT_INT8, sym%LINEAR_SCALING, Min_Reff, Max_Reff, &
                              "micron", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- packed quality flag for ACHA
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Acha_Quality_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Acha_Quality),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "acha_quality", &
                               "not specified", &
                               "quality flags for ACHA products "// &
                               "1:Processed (0=no,1=yes) "// &
                               "2:valid Tc retrieval (1=yes,0=no) "// &
                               "3:valid ec retrieval (1=yes,0=no) "// &
                               "4:valid beta retrieval (1=yes,0=no) "// &
                               "5:degraded Tc retrieval (1=yes,0=no) "// &
                               "6:degraded ec retrieval (1=yes,0=no) "// &
                               "7:degraded beta retrieval (1=yes,0=no) ", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", No_Attribute_Missing_Value, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- packed info flag for ACHA
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Acha_Info_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_ACHA_Info),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "acha_info", &
                               "not specified", &
                               "processing information for ACHA (0=no/1=yes) "// &
                               "1:Cloud Height Attempted "// &
                               "2:Bias Correction Employed "// &
                               "3:Ice Cloud Retrieval "//&
                               "4:Local Radiative Center Processing Used "//&
                               "5:Multi-layer Retrieval "//&
                               "6:Lower Cloud Interpolation Used "//&
                               "7:Boundary Layer Inversion Assumed ", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", No_Attribute_Missing_Value, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- acha O.E. cost
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Acha_Cost_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Acha_Cost),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cost_acha", &
                               "cost_acha", &
                               "final cost function value from acha optimal estimation", &
                                DFNT_INT8, sym%LINEAR_SCALING, Min_Acha_Cost, Max_Acha_Cost, &
                               "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- cloud optical depth from DCOMP
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cod),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d,&
                               "cld_opd_dcomp", &
                               "atmosphere_optical_thickness_due_to_cloud", &
                               "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                               "determined from DCOMP - NOAA CDR", &
                                DFNT_INT16, sym%LINEAR_SCALING, Min_tau, Max_tau, &
                               "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud effective particle radius from DCOMP
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ceps),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cld_reff_dcomp", &
                              "effective_radius_of_cloud_condensed_water_particles_at_cloud_top", &
                              "effective radius of cloud particles determined from DCOMP; "//&
                              "see attributes for channels used - NOAA CDR", &
                              DFNT_INT16, sym%LINEAR_SCALING, Min_Reff, Max_Reff, &
                              "micron", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud optical depth uncertainty from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Dcomp_Uncer_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cod_Dcomp_Uncer),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d,&
                              "cld_opd_dcomp_unc", &
                              "not specified", &
                              "uncertainty in the log10 cloud optical depth at the nominal wavelength "// &
                              "of 0.65 microns, determined from DCOMP; "// &
                              "see attributes for channels used", &
                              DFNT_INT16, sym%LINEAR_SCALING, Min_tau, Max_tau, &
                              "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud effective particle radius uncertainity from the dcomp approach
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Dcomp_Uncer_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ceps_Dcomp_Uncer),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_reff_dcomp_unc", &
                               "not specified", &
                               "uncertainty in the log10 effective radius of cloud particle determined from DCOMP; "// &
                               "see attributes for channels used", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Reff, Max_Reff, &
                               "micron", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif
      
     !--- quality flag for DCOMP Cloud Optical Depth
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Dcomp_Qf_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cod_Dcomp_Qf),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cld_opd_dcomp_qf", &
                               "not specified", &
                               "quality flag for cloud optical depth from DCOMP "// &
                               "not attempted=0, failed=1, low quality=2, high quality=3", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- quality flag for DCOMP Cloud Effective Radius
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Dcomp_Qf_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ceps_Dcomp_Qf),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "cld_reff_dcomp_qf", &
                               "not specified", &
                               "quality flag for cloud effective radius from DCOMP "// &
                               "not attempted=0, failed=1, low quality=2, high quality=3", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- packed quality flag for DCOMP
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Dcomp_Quality_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Dcomp_Quality),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "dcomp_quality", &
                               "not specified", &
                               "quality flags for DCOMP products "// &
                               "1:Processed (0=no,1=yes) "// &
                               "2:valid COD retrieval (0=yes,1=no) "// &
                               "3:valid REF retrieval (0=yes,1=no) "// &
                               "4:degraded COD retrieval (0=no,1=degraded) "// &
                               "5:degraded REF retrieval (0=no,1=degraded) "// &
                               "6:convergency (0=no,1=yes) "// &
                               "7:glint (0=no,1=yes) ", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", No_Attribute_Missing_Value, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- packed info flag for DCOMP 
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Dcomp_Info_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Dcomp_Info),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "dcomp_info", &
                               "not specified", &
                               "processing flags for DCOMP "// &
                               "1:info flag set ? (0=no,1=yes) "// &
                               "2:land/sea mask (0=land,1=sea) "// &
                               "3:day/night mask (0=Day,1=Night) "// &
                               "4:twilight (65-82 solar zenith) (0=no,1=yes) "// &
                               "5:snow (0=no,1= snow) "// &
                               "6:sea-ice (0=no,1=sea-ice) "// &
                               "7:phase (0=liquid,1=ice) "// &
                               "8:thick_cloud (0=no,1=yes) " // &
                               "9:thin_cloud (0=no,1=yes) ", &
                               DFNT_INT16, sym%NO_SCALING, 0.0, 0.0, &
                               "none", No_Attribute_Missing_Value, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- insolation from cloud properties
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_CldInsol_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_CldInsol),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "insolation_dcomp", &
                              "surface_downwelling_shortwave_flux_dcomp", &
                              "surface downwelling shortwave flux computed from the DCOMP cloud properties", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Insol, Max_Insol, "W m-2", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_CldInsol_Dif_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_CldInsol_Dif),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "insolation_diffuse_dcomp", &
                              "surface_downwelling_shortwave_flux_diffuse_dcomp", &
                              "diffuse component of the surface downwelling shortwave flux "// &
                              "computed from the DCOMP cloud properties", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Insol, Max_Insol, "W m-2", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud droplet number concentration
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cdnc_Dcomp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cdnc_Dcomp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cdnc_dcomp", &
                              "cloud_droplet_number_concentration", &
                              "cloud_droplet_number_concentration from DCOMP algorithm", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Cdnc, Max_Cdnc, "cm-3", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud geomterical thickness from DCOMP for adiabatic water clouds
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Hcld_Dcomp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Hcld_Dcomp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "hcld_dcomp", &
                              "geometrical_thickness_of_liquid_clouds", &
                              "geometrical_thickness_of_liquid_clouds based on DCOMP", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Hcld, Max_Hcld, "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- cloud optical depth from NLCOMP

     !--- cloud optical depth from NLCOMP
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Nlcomp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cod_Nlcomp),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d,&
                               "cld_opd_nlcomp", &
                               "atmosphere_optical_thickness_due_to_cloud", &
                               "cloud optical depth at the nominal wavelength of 0.65 microns, "//&
                               "determined from NLCOMP", &
                                DFNT_INT16, sym%LINEAR_SCALING, Min_tau, Max_tau, &
                               "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud effective particle radius from NLCOMP
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Nlcomp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ceps_Nlcomp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cld_reff_nlcomp", &
                              "effective_radius_of_cloud_particle", &
                              "effective radius of cloud particles determined from NLCOMP; "//&
                              "see attributes for channels used", &
                              DFNT_INT16, sym%LINEAR_SCALING, Min_Reff, Max_Reff, &
                              "micron", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud optical depth uncertainty from nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Nlcomp_Uncer_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cod_Nlcomp_Uncer),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d,&
                              "cld_opd_nlcomp_unc", &
                              "atmosphere_optical_thickness_due_to_cloud", &
                              "uncertainty in cloud optical depth at the nominal wavelength "// &
                              "of 0.65 microns, determined from NLCOMP", &
                              DFNT_INT16, sym%LINEAR_SCALING, Min_tau, Max_tau, &
                              "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud effective particle radius uncertainity from the nlcomp approach
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Nlcomp_Uncer_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ceps_Nlcomp_Uncer),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_reff_nlcomp_unc", &
                               "effective_radius_of_cloud_particle", &
                               "effective radius of cloud particle determined from NLCOMP; see attributes for channels used", &
                               DFNT_INT16, sym%LINEAR_SCALING, Min_Reff, Max_Reff, &
                               "micron", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif
      
     !--- packed quality flag for NLCOMP
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Nlcomp_Quality_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Nlcomp_Quality),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "nlcomp_quality", &
                               "nlcomp_quality", &
                               "quality flags for NLCOMP products "// &
                               " see documentation http://cimss.ssec.wisc.edu/clavr/ "// &
                               "1:Processed (0=no,1=yes) "// &
                               "2:valid COD retrieval (0=yes,1=no) "// &
                               "3:valid REF retrieval (0=yes,1=no) "// &
                               "4:degraded COD retrieval (0=no,1=degraded) "// &
                               "5:degraded REF retrieval (0=no,1=degraded) "// &
                               "6:convergency (0=no,1=yes) "// &
                               "7:glint (0=no,1=yes) ", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", No_Attribute_Missing_Value, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- packed info flag for NLCOMP 
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Nlcomp_Info_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Nlcomp_Info),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               "nlcomp_info", &
                               "nlcomp_info", &
                               "processing flags for NLCOMP "// &
                               "see http://cimss.ssec.wisc.edu/clavr/ "// &
                               "1: info flag set ? (0=no,1=yes) "// &
                               "2: land/sea mask (0=land,1=sea) "// &
                               "3: day/night mask (0=Day,1=Night) "// &
                               "4: twilight (65-82 solar zenith) (0=no,1=yes) "// &
                               "5: snow (0=no,1= snow) "// &
                               "6: sea-ice (0=no,1=sea-ice) "// &
                               "7: phase (0=liquid,1=ice) "// &
                               "8: thick_cloud (0=no,1=yes) " // &
                               "9: thin_cloud (0=no,1=yes) ", &
                               DFNT_INT16, sym%NO_SCALING, 0.0, 0.0, &
                               "none", No_Attribute_Missing_Value, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud albedo
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldalb_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cldalb),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cloud_albedo_0_65um_nom", &
                               "cloud_albedo", &
                               "cloud albedo at 0.65 microns nominal from DCOMP", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_albedo, Max_albedo, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud transmission
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldtrn_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cldtrn),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cloud_transmission_0_65um_nom", &
                               "not specified", &
                               "cloud transmission 0.65 microns nominal from DCOMP", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_transmission, Max_transmission, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif
 
     !--- cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldfrac_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_cldfrac),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cloud_fraction", &
                              "cloud_area_fraction", &
                              "cloud fraction computed over a 3x3 pixel array at the native resolution "// &
                              "centered on this pixel", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_frac, Max_frac, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- high cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_High_Cld_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_High_Cld),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "high_cloud_fraction", &
                              "high_cloud_area_fraction", &
                              "high cloud fraction computed over a 3x3 pixel array at the native resolution "// &
                              "centered on this pixel. High clouds have pressures less than 440 hPa.", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Frac, Max_Frac, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- mid cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Mid_Cld_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Mid_Cld),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "mid_cloud_fraction", &
                              "mid_cloud_area_fraction", &
                              "mid cloud fraction computed over a 3x3 pixel array at the native resolution "// &
                              "centered on this pixel. Mid clouds have pressures greater than 440 and "// &
                              "less than 680 hPa", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Frac, Max_Frac, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- low cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Low_Cld_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Low_Cld),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "low_cloud_fraction", &
                              "low_cloud_area_fraction", &
                              "low cloud fraction computed over a 3x3 pixel array at the native resolution "// &
                              "centered on this pixel. Low clouds have pressures greater than 680 hPa", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Frac, Max_Frac, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- cloud fraction uncertainity
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldfrac_Uncer_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cldfrac_Uncer),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cloud_fraction_uncertainty", &
                               "not specified", &
                               "cloud fraction uncertainty computed over a 3x3 array", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_frac, Max_frac, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


!---- Non-cloud Properties

     !--- etrop
     if (Sds_Num_Level2_Etrop_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Etrop),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "emissivity_11um_tropopause", &
                               "not specified", &
                               "emissivity at the nominal wavelength of 11 microns, "// &
                               "assuming the cloud was located at the Tropopause", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Etropo, Max_Etropo, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- ch 1 aerosol optical thickness
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot1_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_aot1),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d, &
                              "aot_0_65um_nom", &
                              "optical_thickness_of_atmosphere_layer_due_to_aerosol_0.65_micron nominal", &
                              "optical thickness of atmosphere layer due to aerosol "// &
                              "at the nominal wavelength of 0.65 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING,Min_aot, Max_aot, &
                              "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- ch 2 aerosol optical thickness
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot2_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Aot2),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d, &
                              "aot_0_86um_nom", &
                              "optical_thickness_of_atmosphere_layer_due_to_aerosol_0.86_micron nominal", &
                              "optical thickness of atmosphere layer due to aerosol "// &
                              "at the nominal wavelength of 0.86 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING,Min_aot, Max_aot, &
                              "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- ch 3a aerosol optical thickness
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot6_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Aot6),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d, &
                              "aot_1_6um_nom", &
                              "optical_thickness_of_atmosphere_layer_due_to_aerosol_1.6_micron nominal", &
                              "optical thickness of atmosphere layer due to aerosol "// &
                              "at the nominal wavelength of 1.6 microns", &
                              DFNT_INT16, sym%LINEAR_SCALING,Min_aot, Max_aot, &
                              "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- aerosol quality flag
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot_Qf_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Aot_Qf),Sd_Id_Level2,Sds_Dims_2d, Sds_Chunk_Size_2d, &
                              "aot_qf", &
                              "optical_thickness_of_atmosphere_layer_quality_flag", &
                              "quality flag for optical thickness of atmosphere layer", &
                              DFNT_INT8, sym%NO_SCALING,0.0, 0.0, &
                              "none", real(Missing_Value_Int1), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- outgoing longwave radiation
     if (Sds_Num_Level2_Olr_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Olr),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "olr", &
                               "toa_net_upward_longwave_flux", &
                               "top of atmosphere outgoing longwave radiation", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_olr, Max_olr, "W m-2", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- insolation
     if (Sasrab_Flag == sym%YES .and. Sds_Num_Level2_Insol_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Insol),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              !"test_product_1", &
                              !"test_product_1", &
                              !"test_product_1", &
                              "insolation_sasrab", &
                              "surface_downwelling_shortwave_flux_sasrab", &
                              "surface downwelling shortwave flux computed from the SASRAB routine", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Insol, Max_Insol, "W m-2", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- insolation diffuse
     if (Sasrab_Flag == sym%YES .and. Sds_Num_Level2_Insol_Dif_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Insol_Dif),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              !"test_product_2", &
                              !"test_product_2", &
                              !"test_product_2", &
                              "insolation_diffuse_sasrab", &
                              "surface_downwelling_shortwave_flux_diffuse_sasrab", &
                              "diffuse component of the surface downwelling shortwave flux "// &
                              "computed from the SASRAB routine", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Insol, Max_Insol, "W m-2", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- ndvi_sfc
     if (Sds_Num_Level2_Ndvi_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ndvi),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "ndvi_sfc", &
                               "normalized_difference_vegetation_index", &
                               "normalized difference vegetation index, atmospherically corrected", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Ndvi, Max_Ndvi, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- ndvi_sfc_from modis white sky
     if (Sds_Num_Level2_Ndvi_White_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ndvi_White),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "ndvi_sfc_white_sky_nwp", &
                               "normalized_difference_vegetation_index_at_surface_from_white_sky_reflectance", &
                               "normalized difference vegetation index, atmospherically corrected, modis white sky", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Ndvi, Max_Ndvi, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- surface temperature
     if (Sds_Num_Level2_Tsfc_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Tsfc),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "surface_temperature_retrieved", &
                               "surface_brightness_temperature", &
                               "surface temperature retrieved using "// &
                               "atmospherically corrected 11 micron radiance", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Tsfc, Max_Tsfc, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- ch31 atmospheric radiance
     if (Sds_Num_Level2_Ch31_Rad_Atm_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch31_Rad_Atm),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "atmos_rad_11_0um_nom_nwp", &
                               "atmospheric_radiance_11_0um_nom_nwp", &
                               "atmospheric emission 11 micron radiance at toa from nwp", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Ch31_Rad_Atm, Max_Ch31_Rad_Atm, "mW/m^2/cm^-1", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- ch31 atmospheric transmission


     if (Sds_Num_Level2_Ch31_Trans_Atm_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch31_Trans_Atm),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "atmos_trans_11_0um_nom_nwp", &
                               "atmospheric_transmission_11_0um_nom_nwp", &
                               "atmospheric tranmission 11 micron radiance at toa from nwp", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Trans, Max_Trans, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- ch31 atmospheric radiance downwelling surface
     if (Sds_Num_Level2_Ch31_Rad_Atm_Dwn_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch31_Rad_Atm_Dwn),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "atmos_rad_dwn_11_0um_nom_nwp", &
                               "atmospheric_radiance_downwelling_11_0um_nom_nwp", &
                               "atmospheric emission 11 micron radiance downward at sfc from nwp", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Ch31_Rad_Atm_Dwn, Max_Ch31_Rad_Atm_Dwn, "mW/m^2/cm^-1", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- ch31 surface emissivity
     if (Sds_Num_Level2_Ch31_Sfc_Emiss_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch31_Sfc_Emiss),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "sfc_emiss_11_0um_nom_nwp", &
                               "surface_emissivity_11_0um_nom_nwp", &
                               "surface emissivity at 11 micron from ancillary data", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Sfc_Ems, Max_Sfc_Ems, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- ch20 surface emissivity
     if (Sds_Num_Level2_Ch20_Sfc_Emiss_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch20_Sfc_Emiss),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "sfc_emiss_3_75um_nom_nwp", &
                               "surface_emissivity_3_75um_nom_nwp", &
                               "surface emissivity at 3.75 micron from ancillary data", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Sfc_Ems, Max_Sfc_Ems, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- surface air temperature
     if (Sds_Num_Level2_Tair_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Tair),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "surface_air_temperature_nwp", &
                               "surface_air_temperature_nwp", &
                               "surface air temperature from NWP ancillary data ", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Tsfc, Max_Tsfc, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- surface radiation temperature
     if (Sds_Num_Level2_Trad_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Trad),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "surface_radiation_temperature_retrieved", &
                               "surface_radiation_temperature_retrieved", &
                               "surface radiation temperature retrieved using "// &
                               "atmospherically corrected 11 micron radiance assuming a black surface", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Tsfc, Max_Tsfc, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- surface temperature background
     if (Sds_Num_Level2_Tsfc_Back_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Tsfc_Back),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "surface_temperature_background", &
                               "surface_temperature", &
                               "surface temperature assumed using "// &
                               "ancillary data sources", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Tsfc, Max_Tsfc, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- surface relative humidity
     if (Sds_Num_Level2_Rh_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Rh),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "surface_relative_humidity_nwp", &
                               "surface_relative_humidity", &
                               "near-surface relative humidity from NWP ancillary data ", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Rh, Max_Rh, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- surface pressure background
     if (Sds_Num_Level2_Psfc_Back_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Psfc_Back),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "surface_pressure_background", &
                               "surface_pressure_background", &
                               "surface pressure assumed using "// &
                               "ancillary data sources", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Psfc, Max_Psfc, "hPa", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- mean sea-level pressure background
     if (Sds_Num_Level2_Pmsl_Back_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Pmsl_Back),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "mean_sealevel_pressure_background", &
                               "mean_sealevel_pressure_background", &
                               "mean sealevel pressure assumed using "// &
                               "ancillary data sources", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Pmsl, Max_Pmsl, "hPa", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- K Index from NWP
     if (Sds_Num_Level2_Kindex_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Kindex),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "k_index_nwp", &
                               "k_index", &
                               "k index computed from "// &
                               "NWP ancillary data sources", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Kindex, Max_Kindex, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Cloud Water Path (NWP)
     if (Sds_Num_Level2_Cwp_Nwp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cwp_Nwp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cloud_water_path_nwp", &
                               "cloud_water_path", &
                               "cloud water path computed from NWP ancillary data sources", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Cwp, Max_Cwp, "g m-2", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Cloud Fraction (NWP)
     if (Sds_Num_Level2_Cfrac_Nwp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cfrac_Nwp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cloud_fraction_nwp", &
                               "cloud_fraction", &
                               "cloud fraction computed from NWP ancillary data sources", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Frac, Max_Frac, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Cloud Top Pressure (NWP)
     if (Sds_Num_Level2_Pc_Nwp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Pc_Nwp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "cld_press_nwp", &
                               "cld_press", &
                               "cloud-top pressure computed from NWP ancillary data sources", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Pc, Max_Pc, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Number of Cloud Layers (NWP)
     if (Sds_Num_Level2_Ncld_Nwp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ncld_Nwp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "number_cloud_layers_nwp", &
                               "number_cloud_layers", &
                               "number cloud layers in column from NWP ancillary data sources", &
                               DFNT_INT8, sym%NO_SCALING, &
                               0.0, 0.0, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---Cloud Type (NWP)
     if (Sds_Num_Level2_Cld_Type_Nwp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cld_Type_Nwp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cloud_type_nwp", &
                              "cloud_type", &
                              "cloud type from NWP ancillary data sources, see PATMOS-x documentation for key", &
                              DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                              "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif



     !--- tropopause temperature background
     if (Sds_Num_Level2_Temp_Tropo_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Temp_Tropo),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "tropopause_temperature_nwp", &
                               "tropopause_temperature", &
                               "tropopause temperature from NWP ancillary data", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Ttropo, Max_Ttropo, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  lifitng condensation level
     if (Sds_Num_Level2_LCL_Nwp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_LCL_Nwp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "lcl_nwp", &
                               "lifting_condensation_level_nwp", &
                               "lifting condensation level from NWP ancillary data", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Zc, Max_Zc, "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  convective condensation level
     if (Sds_Num_Level2_CCL_Nwp_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_CCL_Nwp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "ccl_nwp", &
                               "convective_condensation_level_nwp", &
                               "convective condensation level from NWP ancillary data", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Zc, Max_Zc, "m", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- remote sensing reflectance
     if (Sds_Num_Level2_Rsr_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Rsr),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "remote_sensing_reflectance", &
                               "remote_sensing_reflectance", &
                               "remote sensing reflectance (upward radiance/downward irradiance at surface)", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Rsr, Max_Rsr, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

!---- quality flags
     !--- first level2 quality flags
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Qf1_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Qf1),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "quality_flags_1", &
                               "quality_flags_1", &
                               "first set of packed quality flags, deprecated. Use *_qf variables.", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- second level2 quality flags
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Qf2_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Qf2),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "quality_flags_2", &
                               "quality_flags_2", &
                               "second set of packed quality flags, deprecated. Use *_qf variables.", &
                               DFNT_INT8, sym%NO_SCALING, 0.0, 0.0, &
                               "none", Real(Missing_Value_Int1,kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- channel 1 counts
     if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. Sds_Num_Level2_Ch1_Counts_Flag == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch1_Counts),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_counts", &
                              "not specified", &
                              "instrument counts for the nominal 0.65 micron channel", &
                              DFNT_INT16, sym%NO_SCALING, &
                              0.0, 0.0, "none", Real(Missing_Value_Int2,kind=real4), Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- channel 2 counts
     if (Sensor%Chan_On_Flag_Default(2) == sym%YES .and. Sds_Num_Level2_Ch2_Counts_Flag == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch2_Counts),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_86um_nom_counts", &
                              "not specified", &
                              "instrument counts for the nominal 0.86 micron channel", &
                              DFNT_INT16, sym%NO_SCALING, &
                              0.0, 0.0, "none", Real(Missing_Value_Int2,kind=real4), Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- channel 6 counts
     if (Sensor%Chan_On_Flag_Default(6) == sym%YES .and. Sds_Num_Level2_Ch6_Counts_Flag == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch6_Counts),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_1_60um_nom_counts", &
                              "not specified", &
                              "instrument counts for the nominal 1.6 micron channel", &
                              DFNT_INT16, sym%NO_SCALING, &
                              0.0, 0.0, "none", Real(Missing_Value_Int2,kind=real4), Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  total precipitable water
     if (Sds_Num_Level2_Tpw_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Tpw),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "total_precipitable_water_nwp", &
                              "atmosphere_mass_content_of_water_vapor", &
                              "total precipitable water from NWP ancillary data", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Tpw, Max_Tpw, "cm", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !---  total ozone
     if (Sds_Num_Level2_Ozone_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ozone),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "total_column_ozone_nwp", &
                              "total_column_ozone_nwp", &
                              "total ozone amount integrated over column from nwp ancillary data", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ozone, Max_Ozone, "DU", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- channel 20 reflectance atmospherically corrected
     if (Sds_Num_Level2_Ref_Ch20_Sfc_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch20_Sfc),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_3_75um_nom_atmos_corr", &
                              "toa_bidirectional_pseudo_reflectance_3_75_micron_atmos_corr", &
                              "observed pseudo-reflectance at the nominal wavelength of 3.75 microns, "// &
                              "atmospherically corrected, "// &
                              "expressed as a percentage using PATMOS-x calibration", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch20, Max_Ref_Ch20, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif 
     !--- channel 1 reflectance atmospherically corrected
     if (Sds_Num_Level2_Ref_Ch1_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch1_Sfc),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "refl_0_65um_nom_atmos_corr", &
                               "toa_bidirectional_pseudo_reflectance_0_65_micron_atmos_corr", &
                               "observed pseudo-reflectance at the nominal wavelength of 0.65 microns, "// &
                               "atmospherically corrected, "// &
                               "expressed as a percentage using PATMOS-x calibration", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- channel 2 reflectance atmospherically corrected
     if (Sds_Num_Level2_Ref_Ch2_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(2) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch2_Sfc),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "refl_0_86um_nom_atmos_corr", &
                               "toa_bidirectional_pseudo_reflectance_0_85_micron_atmos_corr", &
                               "observed pseudo-reflectance at the nominal wavelength of 0.85 microns, "// &
                               "atmospherically corrected, "// &
                               "expressed as a percentage using PATMOS-x calibration", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Ref_Ch2, Max_Ref_Ch2, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- channel 6 reflectance atmospherically corrected
     if (Sds_Num_Level2_Ref_Ch6_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch6_Sfc),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "refl_1_60um_nom_atmos_corr", &
                               "toa_bidirectional_pseudo_reflectance_1_60_micron_atmos_corr", &
                               "observed pseudo-reflectance at the nominal wavelength of 1.60 microns, "// &
                               "atmospherically corrected, "// &
                               "expressed as a percentage using PATMOS-x calibration", &
                               DFNT_INT8, sym%LINEAR_SCALING, &
                               Min_Ref_Ch6, Max_Ref_Ch6, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- sst with mask applied for AWIPS display
     if (Sds_Num_Level2_Sst_Masked_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Sst_Masked),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "pixel_sst_masked", &
                              "sea_surface_skin_temperature_masked", &
                              "sea surface skin temperature at the pixel with "// &
                              "land mask and cloud mask applied", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Sst, Max_Sst, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch1 reflectance unnormalized for AWIPS display
     if (Sds_Num_Level2_Ch1_Unnorm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch1_Unnorm),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_unnormalized", &
                              "toa_bidirectional_reflectance_0_65_micron_nominal_unormalized_to_solar_zenith", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.65 microns "// &
                              "unnormalized to the cosine of the solar zenith angle", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch2 reflectance unnormalized for AWIPS display
     if (Sds_Num_Level2_Ch2_Unnorm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(2) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch2_Unnorm),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_86um_nom_unnormalized", &
                              "toa_bidirectional_reflectance_0_86_micron_nominal_unormalized_to_solar_zenith", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.86 microns "// &
                              "unnormalized to the cosine of the solar zenith angle", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch2, Max_Ref_Ch2, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch6 reflectance unnormalized for AWIPS display
     if (Sds_Num_Level2_Ch6_Unnorm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch6_Unnorm),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_1_60um_nom_unnormalized", &
                              "toa_bidirectional_reflectance_1_60_micron_nominal_unormalized_to_solar_zenith", &
                              "top of atmosphere reflectance at the nominal wavelength of 1.60 microns "// &
                              "unnormalized to the cosine of the solar zenith angle", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch6, Max_Ref_Ch6, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ref_Ch1_Clear_Rtm
     if (Sds_Num_Level2_Ch1_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch1_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "refl_0_65um_nom_clear_sky", &
                               "toa_bidirectional_reflectance_assuming_clear_sky_0_65_micron_nominal", &
                               "top of atmosphere bidirectional reflectance modeled assuming clear skies "// &
                               "at the nominal wavelength of 0.65 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch20_Clear_Rtm
     if (Sds_Num_Level2_Ch20_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch20_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_3_75um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_3_75_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 11.0 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt20, Max_Bt20, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch27_Clear_Rtm
     if (Sds_Num_Level2_Ch27_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(27) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch27_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_6_7um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_6_7_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 6.7 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt27, Max_Bt27, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch28_Clear_Rtm
     if (Sds_Num_Level2_Ch28_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(28) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch28_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_7_3um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_7_3_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 7.3 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt28, Max_Bt28, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch29_Clear_Rtm
     if (Sds_Num_Level2_Ch29_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(29) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch29_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_8_5um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_8_5_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 8.5 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt29, Max_Bt29, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch30_Clear_Rtm
     if (Sds_Num_Level2_Ch30_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(30) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch30_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_9_7um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_9_7_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 9.7 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt30, Max_Bt30, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- Bt_Ch31_Clear_Rtm
     if (Sds_Num_Level2_Ch31_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch31_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_11_0um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 11.0 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt31, Max_Bt31, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif
 
     !--- Bt_Ch32_Clear_Rtm
     if (Sds_Num_Level2_Ch32_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(32) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch32_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_12_0um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_12_0_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 12.0 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt32, Max_Bt32, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch33_Clear_Rtm
     if (Sds_Num_Level2_Ch33_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(33) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch33_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_13_3um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_13_3_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 13.3 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt33, Max_Bt33, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch37_Clear_Rtm
     if (Sds_Num_Level2_Ch37_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(37) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch37_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_6_2um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_6_2_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 6.2 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt37, Max_Bt37, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Bt_Ch38_Clear_Rtm
     if (Sds_Num_Level2_Ch38_Clear_Flag == sym%YES .and.  Sensor%Chan_On_Flag_Default(38) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch38_Clear),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                               "temp_10_4um_nom_clear_sky", &
                               "toa_brightness_temperature_assuming_clear_sky_10_4_micron_nominal", &
                               "top of atmosphere brightness temperature modeled assuming clear skies "// &
                               "at the nominal wavelength of 10,4 microns", &
                                DFNT_INT16, sym%LINEAR_SCALING, &
                                Min_Bt38, Max_Bt38, "K", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ref_Ch1_Mean_3x3
     if (Sds_Num_Level2_Ch1_Mean_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch1_Mean),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_mean_3x3", &
                              "not specified", &
                              "mean of the 0.65 micron nominal reflectance computed over a 3x3 "//&
                              "pixel array", &
                               DFNT_INT16, sym%LINEAR_SCALING, &
                               Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
       Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- sst without cloud mask applied for AWIPS display
     if (Sds_Num_Level2_Sst_Unmasked_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Sst_Unmasked),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "pixel_sst_unmasked", &
                              "sea_surface_skin_temperature_unmasked", &
                              "sea surface skin temperature at the pixel with "// &
                              "land mask applied and cloud mask not applied", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Sst, Max_Sst, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- wind speed at 10m agl
     if (Sds_Num_Level2_Wnd_Spd_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Wnd_Spd),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "wind_speed_10m_nwp", &
                              "wind_speed_10m_above_ground", &
                              "wind speed from the NWP ancillary data at 10m above ground level", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Wndspd, Max_Wndspd, "m/s", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- wind direction at 10m agl
     if (Sds_Num_Level2_Wnd_Dir_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Wnd_Dir),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "wind_direction_10m_nwp", &
                              "wind_direction_10m_above_ground", &
                              "wind direction from the NWP ancillary data at 10m above ground level", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Wnddir, Max_Wnddir, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Ch1 Reflectance From Dark Composite
     if (Sds_Num_Level2_Ch1_Dark_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ch1_Dark),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_dark", &
                              "toa_bidirectional_reflectance_0_65_micron_nominal_dark_sky_composite", &
                              "top of atmosphere reflectance at the nominal wavelength of 0.65 microns "//&
                              "generated from a dark-sky compositing method", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Cloud Water Path from DCOMP
     if (Sds_Num_Level2_Cwp_Flag == sym%YES .and. Cld_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Cwp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "cloud_water_path", &
                              "not specified", &
                              "integrated total cloud water over whole column", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Cwp, Max_Cwp, "g m-2", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Rain Rate from DCOMP
     if (Sds_Num_Level2_Rain_Rate_Flag == sym%YES .and. Cld_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Rain_Rate),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "rain_rate", &
                              "rain_rate", &
                              "derived rain rate", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Rain_Rate, Max_Rain_Rate, "mm h-1", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif


     !--- Wind Speed at Cloud Top
     if (Sds_Num_Level2_Wnd_Spd_Cld_Top_Flag == sym%YES .and. Cld_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Wnd_Spd_Cld_Top),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "wind_speed_cloud_top_nwp", &
                              "wind_speed_cloud_top", &
                              "wind speed from the NWP ancillary data at cloud-top level", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Wndspd, Max_Wndspd, "m/s", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Wind Direction at Cloud Top
     if (Sds_Num_Level2_Wnd_Dir_Cld_Top_Flag == sym%YES .and. Cld_Flag == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Wnd_Dir_Cld_Top),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "wind_direction_cloud_top_nwp", &
                              "wind_direction_cloud_top", &
                              "wind direction from the NWP ancillary data at cloud-top level", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Wnddir, Max_Wnddir, "degrees", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Max of ChI1 at the M-band resolution
     if (Sds_Num_Level2_Ref_Max_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Max_ChI1),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_iband_max", &
                              "refl_0_65um_nom_iband_max", &
                              "refl_0_65um_nom_iband_max", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Min of ChI1 at the M-band resolution
     if (Sds_Num_Level2_Ref_Min_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Min_ChI1),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_iband_min", &
                              "refl_0_65um_nom_iband_min", &
                              "refl_0_65um_nom_iband_min", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Mean of ChI1 at the M-band resolution
     if (Sds_Num_Level2_Ref_Mean_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Mean_ChI1),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_iband_mean", &
                              "refl_0_65um_nom_iband_mean", &
                              "refl_0_65um_nom_iband_mean", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch1, Max_Ref_Ch1, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Uniformity Parameter of ChI1 at the M-band resolution
     if (Sds_Num_Level2_Ref_Uni_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Uni_ChI1),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_65um_nom_iband_uni", &
                              "refl_0_65um_nom_iband_uni", &
                              "refl_0_65um_nom_iband_uni", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Uni_Ch1, Max_Uni_Ch1, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Max of ChI2 at the M-band resolution
     if (Sds_Num_Level2_Ref_Max_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Max_ChI2),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_86um_nom_iband_max", &
                              "refl_0_86um_nom_iband_max", &
                              "refl_0_86um_nom_iband_max", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch2, Max_Ref_Ch2, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Min of ChI2 at the M-band resolution
     if (Sds_Num_Level2_Ref_Min_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Min_ChI2),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_86um_nom_iband_min", &
                              "refl_0_86um_nom_iband_min", &
                              "refl_0_86um_nom_iband_min", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch2, Max_Ref_Ch2, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Mean of ChI2 at the M-band resolution
     if (Sds_Num_Level2_Ref_Mean_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Mean_ChI2),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_86um_nom_iband_mean", &
                              "refl_0_86um_nom_iband_mean", &
                              "refl_0_86um_nom_iband_mean", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch2, Max_Ref_Ch2, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Uniformity Parameter of ChI2 at the M-band resolution
     if (Sds_Num_Level2_Ref_Uni_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Uni_ChI2),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_0_86um_nom_iband_uni", &
                              "refl_0_86um_nom_iband_uni", &
                              "refl_0_86um_nom_iband_uni", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Uni_Ch2, Max_Uni_Ch2, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Max of ChI3 at the M-band resolution
     if (Sds_Num_Level2_Ref_Max_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Max_ChI3),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_1_61um_nom_iband_max", &
                              "refl_1_61um_nom_iband_max", &
                              "refl_1_61um_nom_iband_max", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch6, Max_Ref_Ch6, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Min of ChI3 at the M-band resolution
     if (Sds_Num_Level2_Ref_Min_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Min_ChI3),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_1_61um_nom_iband_min", &
                              "refl_1_61um_nom_iband_min", &
                              "refl_1_61um_nom_iband_min", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch6, Max_Ref_Ch6, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Mean of ChI3 at the M-band resolution
     if (Sds_Num_Level2_Ref_Mean_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then                         
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Mean_ChI3),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_1_61um_nom_iband_mean", &
                              "refl_1_61um_nom_iband_mean", &
                              "refl_1_61um_nom_iband_mean", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Ref_Ch6, Max_Ref_Ch6, "%", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Uniformity Parameter of ChI3 at the M-band resolution                                                                      
     if (Sds_Num_Level2_Ref_Uni_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Ref_Uni_ChI3),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "refl_1_61um_nom_iband_uni", &
                              "refl_1_61um_nom_iband_uni", &
                              "refl_1_61um_nom_iband_uni", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Uni_Ch6, Max_Uni_Ch6, "none", Missing_Value_Real4, Istatus)                                        
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Max of ChI4 at the M-band resolution
     if (Sds_Num_Level2_Bt_Max_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt_Max_ChI4),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_3_74um_nom_iband_max", &
                              "temp_3_74um_nom_iband_max", &
                              "temp_3_74um_nom_iband_max", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Bt20, Max_Bt20, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Min of ChI4 at the M-band resolution
     if (Sds_Num_Level2_Bt_Min_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt_Min_ChI4),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_3_74um_nom_iband_min", &
                              "temp_3_74um_nom_iband_min", &
                              "temp_3_74um_nom_iband_min", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Bt20, Max_Bt20, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Mean of ChI4 at the M-band resolution
     if (Sds_Num_Level2_Bt_Mean_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then                          
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt_Mean_ChI4),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_3_74um_nom_iband_mean", &
                              "temp_3_74um_nom_iband_mean", &
                              "temp_3_74um_nom_iband_mean", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Bt20, Max_Bt20, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Uniformity Parameter of ChIs4 at the M-band resolution
     if (Sds_Num_Level2_Bt_Uni_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt_Uni_ChI4),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_3_74um_nom_iband_uni", &
                              "temp_3_74um_nom_iband_uni", &
                              "temp_3_74um_nom_iband_uni", &
                              DFNT_INT8, sym%LINEAR_SCALING, & 
                              Min_Uni_Ch5, Max_Uni_Ch5, "none", Missing_Value_Real4, Istatus)                                        
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Max of ChI5 at the M-band resolution
     if (Sds_Num_Level2_Bt_Max_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt_Max_ChI5),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_11_0um_nom_iband_max", &
                              "temp_11_0um_nom_iband_max", &
                              "temp_11_0um_nom_iband_max", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Bt31, Max_Bt31, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Min of ChI5 at the M-band resolution
     if (Sds_Num_Level2_Bt_Min_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt_Min_ChI5),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_11_0um_nom_iband_min", &
                              "temp_11_0um_nom_iband_min", &
                              "temp_11_0um_nom_iband_min", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Bt31, Max_Bt31, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Mean of ChI5 at the M-band resolution
     if (Sds_Num_Level2_Bt_Mean_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt_Mean_ChI5),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_11_0um_nom_iband_mean", &
                              "temp_11_0um_nom_iband_mean", &
                              "temp_11_0um_nom_iband_mean", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Bt31, Max_Bt31, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- Uniformity Parameter of ChI5 at the M-band resolution
     if (Sds_Num_Level2_Bt_Uni_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt_Uni_ChI5),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_11_0um_nom_iband_uni", &
                              "temp_11_0um_nom_iband_uni", &
                              "temp_11_0um_nom_iband_uni", &
                              DFNT_INT8, sym%LINEAR_SCALING, &
                              Min_Uni_Ch5, Max_Uni_Ch5, "none", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- 3.75 micron BT from Sounder
     if (Sds_Num_Level2_Bt375_Snd_Flag == sym%YES .and. index(Sensor%Sensor_Name,'IFF') > 0) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt375_Snd),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_3_75um_nom_sounder", &
                              "temp_3_75um_nom_sounder", &
                              "3.75 micron brightness temperature observed by the sounder", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt20, Max_Bt20, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- 11 micron BT from Sounder
     if (Sds_Num_Level2_Bt11_Snd_Flag == sym%YES .and. index(Sensor%Sensor_Name,'IFF') > 0) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt11_Snd),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_11_0um_nom_sounder", &
                              "temp_11_0um_nom_sounder", &
                              "11 micron brightness temperature observed by the sounder", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt31, Max_Bt31, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- 12 micron BT from Sounder
     if (Sds_Num_Level2_Bt12_Snd_Flag == sym%YES .and. index(Sensor%Sensor_Name,'IFF') > 0) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_Bt12_Snd),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "temp_12_0um_nom_sounder", &
                              "temp_12_0um_nom_sounder", &
                              "12 micron brightness temperature observed by the sounder", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Bt32, Max_Bt32, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- MJH AVHRR/HIRS IFF auxiliary products
     if (Sds_Num_Level2_HIRS_Cld_Temp_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_HIRS_Cld_Temp),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "hirs_cloud_temp", &
                              "hirs_cloud_temp", &
                              "Menzel's HIRS cloud top temperature from AVHRR-IFF", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Tc, Max_Tc, "K", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif
     if (Sds_Num_Level2_HIRS_Cld_Pres_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_HIRS_Cld_Pres),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "hirs_cloud_pres", &
                              "hirs_cloud_pres", &
                              "Menzel's HIRS cloud top pressure from AVHRR-IFF", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Pc, Max_Pc, "hPa", Missing_Value_Real4, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif
     if (Sds_Num_Level2_HIRS_Cld_Height_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_HIRS_Cld_Height),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "hirs_cloud_height", &
                              "hirs_cloud_height", &
                              "Menzel's HIRS cloud top height from AVHRR-IFF", &
                              DFNT_INT16, sym%LINEAR_SCALING, &
                              Min_Zc, Max_Zc, "m", Missing_Value_Real4, Istatus) ! following cld_height_acha 0 to 20,000m
      Istatus_Sum = Istatus_Sum + Istatus
     endif
     if (Sds_Num_Level2_HIRS_Mask_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_HIRS_Mask),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "hirs_mask", &
                              "hirs_mask", &
                              "1=there is real HIRS data covering this AVHRR pixel; 0=there is no (real) HIRS data covering this AVHRR pixel", &
                              DFNT_UINT8, sym%NO_SCALING, &
                              0., 1., "none", No_Attribute_Missing_Value, Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif
     if (Sds_Num_Level2_HIRS_ele_index_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_HIRS_ele_index),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "hirs_ele_index", &
                              "hirs_ele_index", &
                              "Element # in original HIRS granule for HIRS data at this AVHRR pixel", &
                              DFNT_INT16, sym%NO_SCALING, &
                              0., 1000., "none", Real(Missing_Value_Int2, kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif
     if (Sds_Num_Level2_HIRS_line_index_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
      call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(Sds_Num_Level2_HIRS_line_index),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              "hirs_line_index", &
                              "hirs_line_index", &
                              "Line # in original HIRS granule for HIRS data at this AVHRR pixel", &
                              DFNT_INT16, sym%NO_SCALING, &
                              0., 15000., "none", Real(Missing_Value_Int2, kind=real4), Istatus)
      Istatus_Sum = Istatus_Sum + Istatus
     endif

     !--- check for and report errors
     if (Istatus_Sum /= 0) then
       print *, EXE_PROMPT, MOD_PROMPT, "Error defining sds in level2 hdf file"
       stop
     endif
  endif

end subroutine DEFINE_HDF_FILE_STRUCTURES

!====================================================================
! SUBROUTINE Name: WRITE_PIXEL_HDF_RECORDS
!
! Function:
!   Writes output to appropriate Level 2 files
!
! Description:
!   This subroutine, given the flags that determine which files are 
!   being created, outputs the various pixel level outputs to the
!   appropriate level 2 files and appropriate SDSs for a given segment
!
!====================================================================
subroutine WRITE_PIXEL_HDF_RECORDS(Rtm_File_Flag,Level2_File_Flag)

 integer, intent(in):: Rtm_File_Flag
 integer, intent(in):: Level2_File_Flag
 integer:: Istatus
 integer:: Line_Idx

! HDF function declarations
 integer:: sfwdata

!-----------------------------------------------------------------------
! Get time of each scan line and convert to scale
!-----------------------------------------------------------------------
  where(Scan_Time == Missing_Value_Real4)
            Utc_Scan_Time_Hours = Missing_Value_Real4
  elsewhere
            Utc_Scan_Time_Hours = Scan_Time/60.0/60.0/1000.0
  endwhere

!--------------------------------------------------------------------------------
! determine start and edges for writing this segment's data to pixel hdf file
!-----------------------------------------------------------------------------
    Sds_Start_2d(1) = 0     !pixel dimension
    Sds_Start_2d(2) = Num_Scans_Level2_Hdf

    Sds_Stride_2d(1) = 1
    Sds_Stride_2d(2) = 1

    Sds_Edge_2d(1) = Image%Number_Of_Elements
    Sds_Edge_2d(2) = min(Image%Number_Of_Lines_Read_This_Segment,Image%Number_Of_Lines - Sds_Start_2d(2))

    if (Sds_Edge_2d(2) <= 0) then
      return
    endif

    !--- update Num_Scans_Level2_Hdf
    Num_Scans_Level2_Hdf = min(Image%Number_Of_Lines,Num_Scans_Level2_Hdf +  &
                               Image%Number_Of_Lines_Read_This_Segment)
!-------------------------------------------------------------------------
! write to rtm file
!-------------------------------------------------------------------------
   if (Rtm_File_Flag == sym%YES) then

   !--- reset status flag
   Istatus = 0

   !--- scan line number
   Istatus = sfwdata(Sds_Id_Rtm(1), Sds_Start_2d(2), Sds_Stride_2d(2), &
                        Sds_Edge_2d(2), scan_number(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

   !--- scan line time
   Istatus = sfwdata(Sds_Id_Rtm(2), Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                        Utc_Scan_Time_Hours(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

   !--- surface temperature
   call SCALE_VECTOR_I2_RANK2(Tsfc_Nwp_Pix,sym%LINEAR_SCALING,Min_Tsfc,Max_Tsfc,Missing_Value_Real4,Two_Byte_Temp)
   Istatus = sfwdata(Sds_Id_Rtm(3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                     Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

   !--- computed Ch20 temperature
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(20)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt20,Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- computed Ch20 temperature with solar
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Bt_Clear_Ch20_solar_Rtm,sym%LINEAR_SCALING,Min_Bt20,Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- computed Ch31 temperature
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(31)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt31,Max_Bt31,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(6), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- computed Ch32 temperature
   if (Sensor%Chan_On_Flag_Default(32) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(32)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt32,Max_Bt32,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(7), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- etrop
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(31)%Emiss_Tropo,sym%LINEAR_SCALING,Min_Etropo,Max_Etropo,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(8), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Ref_Ch1 clear
   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(1)%Ref_Toa_Clear,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(9), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                      Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Ref_Ch2 clear
   if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(2)%Ref_Toa_Clear,sym%LINEAR_SCALING,Min_Ref_Ch2,Max_Ref_Ch2,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(10), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Bt_Ch31_std_3x3
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Bt_Ch31_std_3x3,sym%LINEAR_SCALING,Min_Bt31_std,Max_Bt31_std,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(11), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Bt_Ch31_Max_3x3
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Bt_Ch31_Max_3x3,sym%LINEAR_SCALING,Min_Bt31,Max_Bt31,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(12), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Ref_Ch1_std_3x3
   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ref_Ch1_std_3x3,sym%LINEAR_SCALING,Min_Ref_Ch1_std,Max_Ref_Ch1_std,&
                                Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(13), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Ref_Ch1_Min_3x3
   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
    call SCALE_VECTOR_I2_RANK2(Ref_Ch1_Min_3x3,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
    Istatus = sfwdata(Sds_Id_Rtm(14), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                      Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Btd_Ch31_Ch32_Bt_Ch31_Max_3x3
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and.  Sensor%Chan_On_Flag_Default(32)==sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Btd_Ch31_Ch32_Bt_Ch31_Max_3x3,sym%LINEAR_SCALING, &
                                Min_Btd_Ch31_Ch32,Max_Btd_Ch31_Ch32,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(15), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Ems_Ch20
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ems_Ch20,sym%LINEAR_SCALING,Min_Ems_Ch20,Max_Ems_Ch20, &
                                Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(16), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Ems_Ch20 clear-sky
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ems_Ch20_clear_solar_Rtm,sym%LINEAR_SCALING,Min_Ems_Ch20,Max_Ems_Ch20, &
                                Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(17), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Ems_Ch20 median 3x3
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ems_Ch20_median_3x3,sym%LINEAR_SCALING,Min_Ems_Ch20,Max_Ems_Ch20, &
                                Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(18), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Tc_opaque_cloud
   call SCALE_VECTOR_I2_RANK2(Tc_Opaque_Cloud,sym%LINEAR_SCALING,Min_Tc,Max_Tc,Missing_Value_Real4,Two_Byte_Temp)
   Istatus = sfwdata(Sds_Id_Rtm(19), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                     Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

   !--- Ch20 temperature
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Bt_Ch20_median_3x3,sym%LINEAR_SCALING,Min_Bt20,Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(20), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Ch20 surface emissiviy
   call SCALE_VECTOR_I2_RANK2(ch(20)%Sfc_Emiss,sym%LINEAR_SCALING,Min_Sfc_Ems,Max_Sfc_Ems,Missing_Value_Real4,Two_Byte_Temp)
   Istatus = sfwdata(Sds_Id_Rtm(21), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                     Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

   !-- background sst uniformity
   call SCALE_VECTOR_I2_RANK2(sst_anal_uni,sym%LINEAR_SCALING,Min_Sst_std,Max_Sst_std,Missing_Value_Real4,Two_Byte_Temp)
   Istatus = sfwdata(Sds_Id_Rtm(22), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                     Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

   !--- computed Ch27 temperature
   if (Sensor%Chan_On_Flag_Default(27) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(27)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt27,Max_Bt27,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(23), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- computed Ch28 temperature
   if (Sensor%Chan_On_Flag_Default(28) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(28)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt28,Max_Bt28,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(24), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- computed Ch29 temperature
   if (Sensor%Chan_On_Flag_Default(29) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(29)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt29,Max_Bt29,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(25), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- computed Ch33 temperature
   if (Sensor%Chan_On_Flag_Default(33) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(33)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt33,Max_Bt33,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Rtm(26), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !-- 11um and 6.7um covariance
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
   call SCALE_VECTOR_I2_RANK2(Covar_Ch27_Ch31_5x5,sym%LINEAR_SCALING,Min_Bt_Covar,Max_Bt_Covar,Missing_Value_Real4,Two_Byte_Temp)
   Istatus = sfwdata(Sds_Id_Rtm(27), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                     Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
   endif

   !--- Zc_opaque_cloud
   call SCALE_VECTOR_I2_RANK2(Zc_Opaque_Cloud,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
   Istatus = sfwdata(Sds_Id_Rtm(28), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                     Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

   !---  bayes_sfc_mask
   Istatus = sfwdata(Sds_Id_Rtm(29), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                     Bayes_Mask_Sfc_Type_Global(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

   !--- check for and report errors
   if (Istatus /= 0) then
      print *, EXE_PROMPT, MOD_PROMPT, "Error writing to rtm file: ", Istatus
   endif

   endif

!-------------------------------------------------------------------------
! write to level2 file
!-------------------------------------------------------------------------
   if (Level2_File_Flag == sym%YES) then

      Istatus = 0
      
      !--- scan number
      if (Sds_Num_Level2_Scanline_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Scanline), Sds_Start_2d(2), Sds_Stride_2d(2),          &
                         Sds_Edge_2d(2), scan_number(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- scan time
      if (Sds_Num_Level2_Time_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Time), Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                         Utc_Scan_Time_Hours(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
 
      !--- bad scan flag
      if (Sds_Num_Level2_Bad_Scan_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bad_Scan), Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                         Bad_Scan_Flag(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- ascending/descending flag
      if (Sds_Num_Level2_Asc_Flag_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_asc_flag), Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                 Nav%Ascend(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Bad Pixel Mask
      if (Sds_Num_Level2_Bad_Pixel_Mask_Flag == sym%YES) then
       One_Byte_Temp = Bad_Pixel_Mask
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bad_Pixel_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Gap Pixel Mask
      if (Sds_Num_Level2_Gap_Pixel_Mask_Flag == sym%YES) then
       One_Byte_Temp = Gap_Pixel_Mask
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Gap_Pixel_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- diagnostic field #1
      if (Sds_Num_Level2_Diag1_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Diag1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                 Diag_Pix_Array_1(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- diagnostic field #2
      if (Sds_Num_Level2_Diag2_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Diag2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                 Diag_Pix_Array_2(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- diagnostic field #3
      if (Sds_Num_Level2_Diag3_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Diag3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                 Diag_Pix_Array_3(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif


      !--- packed pixel metadata
      if (Sds_Num_Level2_Meta_Data_Flag == sym%YES) then
       One_Byte_Temp = 0
       Temp_Mask = 0
       do Line_Idx = 1, Image%Number_Of_Lines_Per_Segment
         Temp_Mask(:,Line_Idx) = Sensor%Chan_On_Flag_Per_Line(6,Line_Idx)
       enddo
       One_Byte_Temp = ishft(Bayes_Mask_Sfc_Type_Global,3) + ishft(Temp_Mask,2)+ &
                       ishft(solar_contamination_mask,1) + bad_pixel_mask
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Meta_Data), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- latitude
      if (Sds_Num_Level2_Lat_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Nav%Lat,sym%LINEAR_SCALING,Min_Lat,Max_Lat,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Lat),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d, &
                  Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- longitude
      if (Sds_Num_Level2_Lon_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Nav%Lon,sym%LINEAR_SCALING,Min_Lon,Max_Lon,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_lon),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                   Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- latitude parallax_corrected
      if (Sds_Num_Level2_Latpc_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Nav%Lat_Pc,sym%LINEAR_SCALING,Min_Lat,Max_Lat,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Latpc),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d, &
                  Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- longitude parallax corrected
      if (Sds_Num_Level2_Lonpc_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Nav%Lon_Pc,sym%LINEAR_SCALING,Min_Lon,Max_Lon,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Lonpc),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                   Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- sensor zenith
      if (Sds_Num_Level2_Zen_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Satzen,sym%LINEAR_SCALING,Min_Zen,Max_Zen,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_zen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- solar zenith
      if (Sds_Num_Level2_Solzen_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Solzen,sym%LINEAR_SCALING,Min_Solzen,Max_Solzen,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Solzen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- relative azimuth 
      if (Sds_Num_Level2_Relaz_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Relaz,sym%LINEAR_SCALING,Min_Relaz,Max_Relaz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Relaz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- solar azimuth
      if (Sds_Num_Level2_Solaz_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Solaz,sym%LINEAR_SCALING,Min_Solaz,Max_Solaz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Solaz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- sensor/satellite azimuth
      if (Sds_Num_Level2_Sataz_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Sataz,sym%LINEAR_SCALING,Min_Sataz,Max_Sataz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Sataz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif


      !--- glint zenith
      if (Sds_Num_Level2_Glintzen_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Glintzen,sym%LINEAR_SCALING,Min_Glintzen,Max_Glintzen,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Glintzen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- scattering zenith
      if (Sds_Num_Level2_Scatzen_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Scatangle,sym%LINEAR_SCALING,Min_Scatang,Max_Scatang,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Scatzen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- lunar zenith
      if (Sds_Num_Level2_Lunzen_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Lunzen,sym%LINEAR_SCALING,Min_Solzen,Max_Solzen,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Lunzen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- lunar relative azimuth 
      if (Sds_Num_Level2_LunRelaz_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%LunRelaz,sym%LINEAR_SCALING,Min_Relaz,Max_Relaz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_LunRelaz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- lunar azimuth 
      if (Sds_Num_Level2_Lunaz_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Lunaz,sym%LINEAR_SCALING,Min_Solaz,Max_Solaz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Lunaz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- packed land cover (land,snow,coast masks)
      if (Sds_Num_Level2_Packed_Land_Flag == sym%YES) then
       One_Byte_Temp = 0
       One_Byte_Temp = ishft(Sfc%Land,5) + ishft(Sfc%Snow,3) + Sfc%Coast_Mask
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Packed_Land), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- glint mask
      if (Sds_Num_Level2_Glint_Mask_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Glint_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Glint_Mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- coast mask
      if (Sds_Num_Level2_Coast_Mask_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Coast_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Coast_Mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- surface type
      if (Sds_Num_Level2_Sfc_Type_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Sfc_Type), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Sfc_Type(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- land classification
      if (Sds_Num_Level2_Land_Mask_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Land_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Land(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Snow classification
      if (Sds_Num_Level2_Snow_Mask_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Snow_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Snow(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- surface elevation
      if (Sds_Num_Level2_Zsfc_Flag == sym%YES) then
          call SCALE_VECTOR_I2_RANK2(Sfc%Zsfc,sym%LINEAR_SCALING,Min_Zsfc,Max_Zsfc,Missing_Value_Real4,Two_Byte_Temp)
          Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Zsfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                    Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

!--------------------------------------------------------------------------------------------------
!--- observations
!--------------------------------------------------------------------------------------------------

      !-- Ch1 reflectance
      if (Sds_Num_Level2_Ch1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(ch(1)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

      endif

      !-- Ch2 reflectance
      if (Sds_Num_Level2_Ch2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(2) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(2)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch3 reflectance
      if (Sds_Num_Level2_Ch3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(3) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(3)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch3,Max_Ref_Ch3,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch4 reflectance
      if (Sds_Num_Level2_Ch4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(4) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(4)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch4,Max_Ref_Ch4,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch5 reflectance
      if (Sds_Num_Level2_Ch5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(5) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(5)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch5,Max_Ref_Ch5,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch6 reflectance
      if (Sds_Num_Level2_Ch6_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(6)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch6,Max_Ref_Ch6,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch6), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch7 reflectance
      if (Sds_Num_Level2_Ch7_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(7) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(7)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch7,Max_Ref_Ch7,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch7), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch8 reflectance
      if (Sds_Num_Level2_Ch8_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(8) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(8)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch8,Max_Ref_Ch8,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch8), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch9 reflectance
      if (Sds_Num_Level2_Ch9_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(9) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(9)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch9,Max_Ref_Ch9,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch9), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch17 reflectance
      if (Sds_Num_Level2_Ch17_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(17) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(17)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch17,Max_Ref_Ch17,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch17), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch18 reflectance
      if (Sds_Num_Level2_Ch18_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(18) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(18)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch18,Max_Ref_Ch18,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch18), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch19 reflectance
      if (Sds_Num_Level2_Ch19_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(19) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(19)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch19,Max_Ref_Ch19,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch19), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch26 reflectance
      if (Sds_Num_Level2_Ch26_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(26) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(26)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch26,Max_Ref_Ch26,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch26), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- chDNB reflectance
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES .and. Sds_Num_Level2_ChDNB_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(44)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_ChDNB,Max_Ref_ChDNB,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_ChDNB), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- chDNB reflectance lunar
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES .and. Sds_Num_Level2_ChDNB_lunar_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(44)%Ref_Lunar_Toa,sym%LINEAR_SCALING,Min_Ref_ChDNB_lunar,Max_Ref_ChDNB_lunar,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_ChDNB_lunar), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch20 reflectance
      if (Sds_Num_Level2_Ch20_Ref_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(20)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch20,Max_Ref_Ch20,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Ref), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
      !--- Ch20 temperature
      if (Sds_Num_Level2_Ch20_Bt_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(20)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt20,Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Bt), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
      !--- Ch22 temperature
      if (Sds_Num_Level2_Ch22_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(22) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(22)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt22,Max_Bt22,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch22), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
      !--- Ch27 temperature
      if (Sds_Num_Level2_Ch27_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(27)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt27,Max_Bt27,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch27), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
      !--- Ch28 temperature
      if (Sds_Num_Level2_Ch28_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(28) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(28)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt28,Max_Bt28,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch28), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch29 temperature
      if (Sds_Num_Level2_Ch29_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(29) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(29)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt29,Max_Bt29,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch29), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch30 temperature
      if (Sds_Num_Level2_Ch30_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(30) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(30)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt30,Max_Bt30,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch30), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch31 temperature
      if (Sds_Num_Level2_Ch31_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(31)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt31,Max_Bt31,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch32 temperature
      if (Sds_Num_Level2_Ch32_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(32)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt32,Max_Bt32,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch32), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch33 temperature
      if (Sds_Num_Level2_Ch33_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(33) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(33)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt33,Max_Bt33,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch33), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch34 temperature
      if (Sds_Num_Level2_Ch34_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(34) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(34)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt34,Max_Bt34,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch34), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch35 temperature
      if (Sds_Num_Level2_Ch35_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(35) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(35)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt35,Max_Bt35,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch35), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch36 temperature
      if (Sds_Num_Level2_Ch36_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(36) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(36)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt36,Max_Bt36,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch36), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch37 temperature
      if (Sds_Num_Level2_Ch37_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(37) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(37)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt37,Max_Bt37,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch37), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch38 temperature
      if (Sds_Num_Level2_Ch38_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(38)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt38,Max_Bt38,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch38), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch1_std_3x3
      if (Sds_Num_Level2_Ch1_Std_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Ref_Ch1_std_3x3,sym%LINEAR_SCALING,Min_Ref_Ch1_std,Max_Ref_Ch1_std, &
                                 Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Std), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Bt_Ch31_std_3x3
      if (Sds_Num_Level2_Ch31_Std_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Bt_Ch31_Std_3x3,sym%LINEAR_SCALING,Min_Bt31_std,Max_Bt31_std, &
                                 Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Std), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                      One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif


!-------- cloud properties

     !--- cloud probability
     if (Sds_Num_Level2_Cldprob_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(posterior_cld_probability, &
                                 sym%LINEAR_SCALING,Min_frac,Max_frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cldprob), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld mask
     if (Sds_Num_Level2_Cld_Mask_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Mask), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        cld_mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- adjacent pixel cloud mask
     if (Sds_Num_Level2_Adj_Pix_Cld_Mask_Flag == sym%YES) then
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Adj_Pix_Cld_Mask), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Adj_Pix_Cld_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld mask test vector (first byte - acm only)
     if (Sds_Num_Level2_Cld_Tests_Flag == sym%YES) then     

       Sds_Start_3d(1) = 0
       Sds_Start_3d(2) = 0
       Sds_Start_3d(3) = Sds_Start_2d(2)

       Sds_Stride_3d(1) = 1
       Sds_Stride_3d(2) = 1
       Sds_Stride_3d(3) = 1

       Sds_Edge_3d(1) = Max_Num_Cld_Test_Bytes
       Sds_Edge_3d(2) = Sds_Edge_2d(1)
       Sds_Edge_3d(3) = Sds_Edge_2d(2)

      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Tests),  &
                        Sds_Start_3d, &
                        Sds_Stride_3d, &
                        Sds_Edge_3d, &
                        Cld_Test_Vector_Packed(:,:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld mask 1b
     if (Sds_Num_Level2_Cld_Mask_Aux_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Mask_Aux), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Mask_Aux(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- bayes mask sfc type
     if (Sds_Num_Level2_Bayes_Sfc_Type_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bayes_Sfc_Type), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Bayes_Mask_Sfc_Type_Global(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- dust mask
     if (Sds_Num_Level2_Dust_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Dust), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Dust_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- smoke mask
     if (Sds_Num_Level2_Smoke_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Smoke), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Smoke_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- shadow mask
     if (Sds_Num_Level2_Shadow_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Shadow), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Shadow_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- fire mask
     if (Sds_Num_Level2_Fire_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Fire), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Fire_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth from mask
     if (Cld_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES .and. Sds_Num_Level2_Cod_Mask_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(ch(1)%Opd,sym%LINEAR_SCALING,Min_Tau,Max_Tau,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld type
     if (Sds_Num_Level2_Cld_Type_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Type), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Type(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld phase
     if (Sds_Num_Level2_Cld_Phase_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Phase), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Phase(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld type
     if (Sds_Num_Level2_Cld_Type_Aux_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Type_Aux), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Type_Aux(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld pressure
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Pc,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_ctp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld temperature
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctt_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Tc,sym%LINEAR_SCALING,Min_Tc,Max_Tc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctt), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Zc,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_cth), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld top height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Top_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Zc_Top,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Top), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld base height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Base_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Zc_Base,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Base), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- lower cld height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Zc_Lower_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Zc_Lower_Cloud,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Zc_Lower), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- lower cld pressure
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Pc_Lower_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Pc_Lower_Cloud,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Pc_Lower), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld altitude
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Alt_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(ACHA%Alt,sym%LINEAR_SCALING,Min_Alt,Max_Alt,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Alt), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- acha processing order
     if (Sds_Num_Level2_Acha_Order_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Order), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        ACHA%Processing_Order(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- acha inversion flag
     if (Sds_Num_Level2_Acha_Inver_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Inver), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        ACHA%Inversion_Flag(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld height from h2o
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_H2O_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Zc_H2O,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_H2O), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld height from opaque
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Opa_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Zc_Opaque_Cloud,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Opa), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld emissivity from split window
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ec_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Ec,sym%LINEAR_SCALING,Min_ec,Max_ec,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_ec), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- beta from split window
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Beta_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Beta,sym%LINEAR_SCALING,Min_beta,Max_beta,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_beta), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud temperature uncertainity from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctt_Acha_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Tc_Uncertainty,sym%LINEAR_SCALING,Min_Tc_Uncer, &
                                 Max_Tc_Uncer,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctt_Acha_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud height uncertainity from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Acha_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Zc_Uncertainty,sym%LINEAR_SCALING,Min_Zc_Uncer, &
                                 Max_Zc_Uncer,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Acha_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud temperature quality flag from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Acha_Qf_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%OE_Quality_Flags(1,:,:)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Acha_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud emissivity quality flag from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ec_Acha_Qf_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%OE_Quality_Flags(2,:,:)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ec_Acha_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud beta quality flag from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Beta_Acha_Qf_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%OE_Quality_Flags(3,:,:)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Beta_Acha_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Acha_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Tau,sym%LINEAR_SCALING,Min_Tau_Acha,Max_Tau_Acha,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Acha), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Acha_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Reff,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Acha), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- quality flag from Acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Acha_Quality_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%Packed_Quality_Flags
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Quality), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- info flag from Acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Acha_Info_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%Packed_Meta_Data_Flags
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Info), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Acha Cost
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Acha_Cost_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Cost,sym%LINEAR_SCALING,Min_Acha_Cost,Max_Acha_Cost,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Cost), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Tau_Dcomp,sym%LINEAR_SCALING,Min_tau,Max_tau,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Reff_Dcomp,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth cost from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Dcomp_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Tau_Dcomp_Cost,sym%LINEAR_SCALING,Min_tau,Max_tau,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Dcomp_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size cost from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Dcomp_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Reff_Dcomp_Cost,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Dcomp_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud optical depth quality flag from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Dcomp_Qf_Flag == sym%YES) then     
      One_Byte_Temp = Tau_Dcomp_Qf
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Dcomp_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud effective radius quality flag from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Dcomp_Qf_Flag == sym%YES) then     
      One_Byte_Temp = Reff_Dcomp_Qf
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Dcomp_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- quality flag from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Dcomp_Quality_Flag == sym%YES) then     
      One_Byte_Temp = Dcomp_quality_flag
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Dcomp_Quality), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- info flag  dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Dcomp_Info_Flag == sym%YES) then     
      Two_Byte_Temp = Dcomp_Info_Flag
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Dcomp_Info), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- insolation from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_CldInsol_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Insolation_Dcomp, &
                                 sym%LINEAR_SCALING,Min_Insol,Max_Insol,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_CldInsol), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- diffuse insolation from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_CldInsol_Dif_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Insolation_Dcomp_Diffuse, &
                                 sym%LINEAR_SCALING,Min_Insol,Max_Insol,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_CldInsol_Dif), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cdnc from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cdnc_Dcomp_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Cdnc_Dcomp, &
                                 sym%LINEAR_SCALING,Min_Cdnc,Max_Cdnc,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cdnc_Dcomp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud geo height from  dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Hcld_Dcomp_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Hcld_Dcomp, &
                                 sym%LINEAR_SCALING,Min_Hcld,Max_Hcld,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Hcld_Dcomp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth from nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Nlcomp_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Tau_Nlcomp,sym%LINEAR_SCALING,Min_tau,Max_tau,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Nlcomp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Nlcomp_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Reff_Nlcomp,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Nlcomp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth cost from nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Nlcomp_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Tau_Nlcomp_Cost,sym%LINEAR_SCALING,Min_Tau,Max_Tau,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Nlcomp_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size cost from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Nlcomp_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Reff_Nlcomp_Cost,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Nlcomp_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- quality flag from nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Nlcomp_Quality_Flag == sym%YES) then     
      One_Byte_Temp = Nlcomp_quality_flag
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Nlcomp_Quality), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- info flag from nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Nlcomp_Info_Flag == sym%YES) then     
      Two_Byte_Temp = Nlcomp_Info_Flag
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Nlcomp_Info), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud albedo
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldalb_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(cloud_063um_albedo, &
                                 sym%LINEAR_SCALING,Min_albedo,Max_albedo,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_cldalb), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud transmission
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldtrn_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(cloud_063um_transmission_solar, &
                                 sym%LINEAR_SCALING,Min_transmission,Max_transmission,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_cldtrn), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldfrac_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(cloud_fraction_3x3, &
                                 sym%LINEAR_SCALING,Min_frac,Max_frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_cldfrac), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud fraction uncertainty
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldfrac_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(Cloud_Fraction_Uncer_3x3, &
                                 sym%LINEAR_SCALING,Min_frac,Max_frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cldfrac_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- high cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_High_Cld_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(High_Cloud_Fraction_3x3, &
                                 sym%LINEAR_SCALING,Min_Frac,Max_Frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_High_Cld), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- mid cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Mid_Cld_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(Mid_Cloud_Fraction_3x3, &
                                 sym%LINEAR_SCALING,Min_Frac,Max_Frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Mid_Cld), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- low cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Low_Cld_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(Low_Cloud_Fraction_3x3, &
                                 sym%LINEAR_SCALING,Min_Frac,Max_Frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Low_Cld), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

!--- non-cloud props     
     !--- etrop
     if (Sds_Num_Level2_Etrop_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ch(31)%Emiss_Tropo,sym%LINEAR_SCALING,Min_Etropo,Max_Etropo,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Etrop), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- aot1
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot1_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(aot1,sym%LINEAR_SCALING,Min_aot,Max_aot,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_aot1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- aot2
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot2_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(aot2,sym%LINEAR_SCALING,Min_aot,Max_aot,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_aot2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- aot6
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot6_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(aot3a,sym%LINEAR_SCALING,Min_aot,Max_aot,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Aot6), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- aot qf
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot_Qf_Flag == sym%YES) then
        One_Byte_Temp = aot_qf
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Aot_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- olr
     if ( Sds_Num_Level2_Olr_Flag == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(olr, &
                                 sym%LINEAR_SCALING,Min_olr,Max_olr,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_olr), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- insolation
     if (Sasrab_Flag == sym%YES .and. Sds_Num_Level2_Insol_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Insolation_All_Sky, &
                                 sym%LINEAR_SCALING,Min_Insol,Max_Insol,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Insol), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- diffuse insolation
     if (Sasrab_Flag == sym%YES .and. Sds_Num_Level2_Insol_Dif_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Insolation_All_Sky_Diffuse, &
                                 sym%LINEAR_SCALING,Min_Insol,Max_Insol,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Insol_Dif), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- ndvi surface corrected
     if (Sds_Num_Level2_Ndvi_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(ndvi_sfc,sym%LINEAR_SCALING,Min_Ndvi,Max_Ndvi,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ndvi), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- ndvi surface corrected from modis white sky
     if (Sds_Num_Level2_Ndvi_White_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Ndvi_Sfc_White_Sky,sym%LINEAR_SCALING,Min_Ndvi,Max_Ndvi,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ndvi_White), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface temperature
     if (Sds_Num_Level2_Tsfc_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Tsfc_Retrieved,sym%LINEAR_SCALING,Min_Tsfc,Max_Tsfc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Tsfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface radiation temperature
     if (Sds_Num_Level2_Trad_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Trad_Retrieved,sym%LINEAR_SCALING,Min_Tsfc,Max_Tsfc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Trad), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface radiation temperature
     if (Sds_Num_Level2_Tair_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Tair_Nwp_Pix,sym%LINEAR_SCALING,Min_Tsfc,Max_Tsfc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Tair), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface temperature background
     if (Sds_Num_Level2_Tsfc_Back_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Tsfc_Nwp_Pix,sym%LINEAR_SCALING,Min_Tsfc,Max_Tsfc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Tsfc_Back), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface relative humidity background
     if (Sds_Num_Level2_Rh_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Rh_Nwp_Pix,sym%LINEAR_SCALING,Min_Rh,Max_Rh,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Rh), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface pressure background
     if (Sds_Num_Level2_Psfc_Back_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Psfc_Nwp_Pix,sym%LINEAR_SCALING,Min_Psfc,Max_Psfc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Psfc_Back), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- sea-level pressure background
     if (Sds_Num_Level2_Pmsl_Back_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Pmsl_Nwp_Pix,sym%LINEAR_SCALING,Min_Pmsl,Max_Pmsl,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Pmsl_Back), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- K index
     if (Sds_Num_Level2_Kindex_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(K_Index_Nwp_Pix,sym%LINEAR_SCALING,Min_Kindex,Max_Kindex,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Kindex), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Cwp from Nwp 
     if (Sds_Num_Level2_Cwp_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Cwp_Nwp_Pix,sym%LINEAR_SCALING,Min_Cwp,Max_Cwp,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cwp_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Cfrac from Nwp 
     if (Sds_Num_Level2_Cfrac_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Cfrac_Nwp_Pix,sym%LINEAR_SCALING,Min_Frac,Max_Frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cfrac_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Pc from Nwp 
     if (Sds_Num_Level2_Pc_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Pc_Nwp_Pix,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Pc_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Number of Cloud Layers from Nwp 
     if (Sds_Num_Level2_Ncld_Nwp_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ncld_Nwp), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Ncld_Layers_Nwp_Pix(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Cld Type from Nwp 
     if (Sds_Num_Level2_Cld_Type_Nwp_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Type_Nwp), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Type_Nwp_Pix(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- tropopause temperature nwp
     if (Sds_Num_Level2_Temp_Tropo_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Ttropo_Nwp_Pix,sym%LINEAR_SCALING,Min_Ttropo,Max_Ttropo,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Temp_Tropo), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- lifting condensation level nwp
     if (Sds_Num_Level2_LCL_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(LCL_Height_Nwp_Pix,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_LCL_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- convective condensation level nwp
     if (Sds_Num_Level2_CCL_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(CCL_Height_Nwp_Pix,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_CCL_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- remote sensing reflectance
     if (Sds_Num_Level2_Rsr_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(rsr,sym%LINEAR_SCALING,Min_Rsr,Max_Rsr,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Rsr), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

!--- quality flags

     !--- first level 2 packed quality flags
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Qf1_Flag == sym%YES) then
      One_Byte_Temp = 0
      One_Byte_Temp = ishft(ACHA%OE_Quality_Flags(1,:,:),6) + ishft(ACHA%OE_Quality_Flags(2,:,:),4) + &
                      ishft(ACHA%OE_Quality_Flags(3,:,:),2) + Tau_Dcomp_Qf
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Qf1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

     endif

     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Qf2_Flag == sym%YES) then
        !--- second level 2 packed quality flags
        One_Byte_Temp = 0
        One_Byte_Temp = ishft(Reff_Dcomp_Qf,6) + ishft(Aot_Qf,4) + ishft(Rsr_Qf,2) + cld_mask

        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Qf2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

    !--- ch1 counts
    if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. Sds_Num_Level2_Ch1_Counts_Flag == sym%YES) then
        Two_Byte_Temp = Ch1_Counts
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Counts), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch2 counts
    if (Sensor%Chan_On_Flag_Default(2) == sym%YES .and. Sds_Num_Level2_Ch2_Counts_Flag == sym%YES) then
        Two_Byte_Temp = Ch2_Counts
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch2_Counts), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch6 counts
    if (Sensor%Chan_On_Flag_Default(6) == sym%YES .and. Sds_Num_Level2_Ch6_Counts_Flag == sym%YES) then
        Two_Byte_Temp = Ch6_Counts
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch6_Counts), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- tpw
    if (Sds_Num_Level2_Tpw_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Tpw_Nwp_Pix,sym%LINEAR_SCALING,Min_Tpw,Max_Tpw,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Tpw), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ozone
    if (Sds_Num_Level2_Ozone_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Ozone_Nwp_Pix,sym%LINEAR_SCALING,Min_Ozone,Max_Ozone,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ozone), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch20 reflectance atmos corrected
    if (Sds_Num_Level2_Ref_Ch20_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(20)%Ref_Sfc,sym%LINEAR_SCALING,Min_Ref_Ch20,Max_Ref_Ch20,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch20_Sfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch1 reflectance atmos corrected
    if (Sds_Num_Level2_Ref_Ch1_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(1)%Ref_Sfc,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch1_Sfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch2 reflectance atmos corrected
    if (Sds_Num_Level2_Ref_Ch2_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(2) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(ch(2)%Ref_Sfc,sym%LINEAR_SCALING,Min_Ref_Ch2,Max_Ref_Ch2,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch2_Sfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch6 reflectance atmos corrected
    if (Sds_Num_Level2_Ref_Ch6_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(ch(6)%Ref_Sfc,sym%LINEAR_SCALING,Min_Ref_Ch6,Max_Ref_Ch6,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch6_Sfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- sst already masked
    if (Sds_Num_Level2_Sst_Masked_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Sst_Masked,sym%LINEAR_SCALING,Min_Sst,Max_Sst,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Sst_Masked), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                      One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif 

    !-- Ch1 reflectance
    if (Sds_Num_Level2_Ch1_Unnorm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(ch(1)%Ref_Toa_Unnorm,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Unnorm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
     
    !-- Ch2 reflectance
    if (Sds_Num_Level2_Ch2_Unnorm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(2) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(ch(2)%Ref_Toa_Unnorm,sym%LINEAR_SCALING,Min_Ref_Ch2,Max_Ref_Ch2,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch2_Unnorm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch6 reflectance
    if (Sds_Num_Level2_Ch6_Unnorm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(ch(6)%Ref_Toa_Unnorm,sym%LINEAR_SCALING,Min_Ref_Ch6,Max_Ref_Ch6,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch6_Unnorm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Ref_Ch1 clear
    if (Sds_Num_Level2_Ch1_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(1)%Ref_Toa_Clear,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Ch20 temperature clear
    if (Sds_Num_Level2_Ch20_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(20)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt20,Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch27 temperature clear
    if (Sds_Num_Level2_Ch27_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(27)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt27,Max_Bt27,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch27_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch28 temperature clear
    if (Sds_Num_Level2_Ch28_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(28) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(28)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt28,Max_Bt28,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch28_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch29 temperature clear
    if (Sds_Num_Level2_Ch29_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(29) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(29)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt29,Max_Bt29,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch29_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch30 temperature clear
    if (Sds_Num_Level2_Ch30_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(30) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(30)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt30,Max_Bt30,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch30_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch31 temperature clear
    if (Sds_Num_Level2_Ch31_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(31)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt31,Max_Bt31,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch32 temperature clear
    if (Sds_Num_Level2_Ch32_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(32)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt32,Max_Bt32,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch32_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch33 temperature clear
    if (Sds_Num_Level2_Ch33_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(33) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(33)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt33,Max_Bt33,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch33_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch37 temperature clear
    if (Sds_Num_Level2_Ch37_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(37) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(37)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt37,Max_Bt37,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch37_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch38 temperature clear
    if (Sds_Num_Level2_Ch38_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(38)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt38,Max_Bt38,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch38_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ref_Ch1 3x3 Mean
    if (Sds_Num_Level2_Ch1_Mean_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ref_Ch1_Mean_3x3,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Mean), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- sst unmasked by cloud
    if (Sds_Num_Level2_Sst_Unmasked_Flag == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Sst_Unmasked,sym%LINEAR_SCALING,Min_Sst,Max_Sst,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Sst_Unmasked), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- wind speed
    if (Sds_Num_Level2_Wnd_Spd_Flag == sym%YES) then
     call SCALE_VECTOR_I1_RANK2(Wnd_Spd_10m_Nwp_Pix,sym%LINEAR_SCALING,Min_Wndspd,Max_Wndspd,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Wnd_Spd), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- wind direction
    if (Sds_Num_Level2_Wnd_Dir_Flag == sym%YES) then
     call SCALE_VECTOR_I1_RANK2(Wnd_Dir_10m_Nwp_Pix,sym%LINEAR_SCALING,Min_Wnddir,Max_Wnddir,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Wnd_Dir), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- wind speed at cld top
    if (Sds_Num_Level2_Wnd_Spd_Cld_Top_Flag == sym%YES .and. Cld_Flag == sym%YES) then
     call SCALE_VECTOR_I1_RANK2(Wnd_Spd_Cld_Top_Nwp_Pix,sym%LINEAR_SCALING,Min_Wndspd,Max_Wndspd,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Wnd_Spd_Cld_Top), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- wind direction at cld top
    if (Sds_Num_Level2_Wnd_Dir_Cld_Top_Flag == sym%YES .and. Cld_Flag == sym%YES) then
     call SCALE_VECTOR_I1_RANK2(Wnd_Dir_Cld_Top_Nwp_Pix,sym%LINEAR_SCALING,Min_Wnddir,Max_Wnddir,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Wnd_Dir_Cld_Top), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch1 Dark-Sky Reflectance
    if (Sds_Num_Level2_Ch1_Dark_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Ref_Ch1_Dark_Composite,sym%LINEAR_SCALING,Min_Ref_Ch1, &
                                   Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Dark), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Cloud Water Path
    if (Sds_Num_Level2_Cwp_Flag == sym%YES .and. Cld_Flag == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Cwp_Dcomp,sym%LINEAR_SCALING,Min_Cwp, &
                                   Max_Cwp,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Rain Rate
    if (Sds_Num_Level2_Rain_Rate_Flag == sym%YES .and. Cld_Flag == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Rain_Rate_Dcomp,sym%LINEAR_SCALING,Min_Rain_Rate, &
                                   Max_Rain_Rate,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Rain_Rate), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I1 at M-band
    if (Sds_Num_Level2_Ref_Max_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Max_ChI1,sym%LINEAR_SCALING,Min_Ref_Ch1, &
                                   Max_Ref_Ch1,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Max_ChI1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Min of I1 at M-band
    if (Sds_Num_Level2_Ref_Min_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Min_ChI1,sym%LINEAR_SCALING,Min_Ref_Ch1, &
                                   Max_Ref_Ch1,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Min_ChI1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I1 at M-band
    if (Sds_Num_Level2_Ref_Mean_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Mean_ChI1,sym%LINEAR_SCALING,Min_Ref_Ch1, &
                                   Max_Ref_Ch1,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Mean_ChI1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I1 at M-band
    if (Sds_Num_Level2_Ref_Uni_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Uni_ChI1,sym%LINEAR_SCALING,Min_Uni_Ch1, &
                                   Max_Uni_Ch1,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Uni_ChI1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I2 at M-band
    if (Sds_Num_Level2_Ref_Max_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Max_ChI2,sym%LINEAR_SCALING,Min_Ref_Ch2, &
                                   Max_Ref_Ch2,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Max_ChI2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Min of I2 at M-band
    if (Sds_Num_Level2_Ref_Min_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Min_ChI2,sym%LINEAR_SCALING,Min_Ref_Ch2, &
                                   Max_Ref_Ch2,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Min_ChI2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I2 at M-band
    if (Sds_Num_Level2_Ref_Mean_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Mean_ChI2,sym%LINEAR_SCALING,Min_Ref_Ch2, &
                                   Max_Ref_Ch2,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Mean_ChI2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I2 at M-band
    if (Sds_Num_Level2_Ref_Uni_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Uni_ChI2,sym%LINEAR_SCALING,Min_Uni_Ch2, &
                                   Max_Uni_Ch2,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Uni_ChI2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I3 at M-band
    if (Sds_Num_Level2_Ref_Max_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Max_ChI3,sym%LINEAR_SCALING,Min_Ref_Ch6, &
                                   Max_Ref_Ch6,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Max_ChI3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Min of I3 at M-band
    if (Sds_Num_Level2_Ref_Min_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Min_ChI3,sym%LINEAR_SCALING,Min_Ref_Ch6, &
                                   Max_Ref_Ch6,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Min_ChI3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I3 at M-band
    if (Sds_Num_Level2_Ref_Mean_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Mean_ChI3,sym%LINEAR_SCALING,Min_Ref_Ch6, &
                                   Max_Ref_Ch6,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Mean_ChI3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I3 at M-band
    if (Sds_Num_Level2_Ref_Uni_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Uni_ChI3,sym%LINEAR_SCALING,Min_Uni_Ch6, &
                                   Max_Uni_Ch6,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Uni_ChI3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I4 at M-band
    if (Sds_Num_Level2_Bt_Max_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Max_ChI4,sym%LINEAR_SCALING,Min_Bt20, &
                                   Max_Bt20,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Max_ChI4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
                                                                                                                                     
    !--- Min of I4 at M-band
    if (Sds_Num_Level2_Bt_Min_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Min_ChI4,sym%LINEAR_SCALING,Min_Bt20, &
                                   Max_Bt20,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Min_ChI4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I4 at M-band
    if (Sds_Num_Level2_Bt_Mean_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Mean_ChI4,sym%LINEAR_SCALING,Min_Bt20, &
                                   Max_Bt20,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Mean_ChI4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I4 at M-band
    if (Sds_Num_Level2_Bt_Uni_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Uni_ChI4,sym%LINEAR_SCALING,Min_Uni_Ch5, &
                                   Max_Uni_Ch5,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Uni_ChI4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I5 at M-band
    if (Sds_Num_Level2_Bt_Max_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Max_ChI5,sym%LINEAR_SCALING,Min_Bt31, &
                                   Max_Bt31,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Max_ChI5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Min of I5 at M-band
    if (Sds_Num_Level2_Bt_Min_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Min_ChI5,sym%LINEAR_SCALING,Min_Bt31, &
                                   Max_Bt31,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Min_ChI5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I5 at M-band
    if (Sds_Num_Level2_Bt_Mean_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Mean_ChI5,sym%LINEAR_SCALING,Min_Bt31, &
                                   Max_Bt31,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Mean_ChI5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I5 at M-band
    if (Sds_Num_Level2_Bt_Uni_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Uni_ChI5,sym%LINEAR_SCALING,Min_Uni_Ch5, &
                                   Max_Uni_Ch5,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Uni_ChI5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif


    !--- ch31 atmospheric radiance
    if (Sds_Num_Level2_Ch31_Rad_Atm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(31)%Rad_Atm,sym%LINEAR_SCALING,Min_Ch31_Rad_Atm, &
                                   Max_Ch31_Rad_Atm,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Rad_Atm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch31 atmospheric transmission
    if (Sds_Num_Level2_Ch31_Trans_Atm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(31)%Trans_Atm,sym%LINEAR_SCALING,Min_Trans, &
                                   Max_Trans,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Trans_Atm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch31 downward atmospheric radiance
    if (Sds_Num_Level2_Ch31_Rad_Atm_Dwn_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(31)%Rad_Atm_Dwn_Sfc,sym%LINEAR_SCALING,Min_Ch31_Rad_Atm_Dwn, &
                                   Max_Ch31_Rad_Atm_Dwn,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Rad_Atm_Dwn), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch31 surface emissivity
    if (Sds_Num_Level2_Ch31_Sfc_Emiss_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(31)%Sfc_Emiss,sym%LINEAR_SCALING,Min_Sfc_Ems, &
                                   Max_Sfc_Ems,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Sfc_Emiss), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch20 surface emissivity
    if (Sds_Num_Level2_Ch20_Sfc_Emiss_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(20)%Sfc_Emiss,sym%LINEAR_SCALING,Min_Sfc_Ems, &
                                   Max_Sfc_Ems,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Sfc_Emiss), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- 3.75 micron BT for Sounder
    if (Sds_Num_Level2_Bt375_Snd_Flag == sym%YES .and. index(Sensor%Sensor_Name,'IFF') > 0) then
        call SCALE_VECTOR_I2_RANK2(Bt_375um_Sounder,sym%LINEAR_SCALING,Min_Bt20, &
                                   Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt375_Snd), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
    !--- 11 micron BT for Sounder
    if (Sds_Num_Level2_Bt11_Snd_Flag == sym%YES .and. index(Sensor%Sensor_Name,'IFF') > 0) then
        call SCALE_VECTOR_I2_RANK2(Bt_11um_Sounder,sym%LINEAR_SCALING,Min_Bt31, &
                                   Max_Bt31,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt11_Snd), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- 12 micron BT for Sounder
    if (Sds_Num_Level2_Bt12_Snd_Flag == sym%YES .and. index(Sensor%Sensor_Name,'IFF') > 0) then
        call SCALE_VECTOR_I2_RANK2(Bt_12um_Sounder,sym%LINEAR_SCALING,Min_Bt32, &
                                   Max_Bt32,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt12_Snd), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- MJH Menzel HIRS cloud temp for AVHRR/HIRS IFF
    if (Sds_Num_Level2_HIRS_Cld_Temp_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
        print *, 'Writing AVHRR/HIRS aux fields to level2'
        call SCALE_VECTOR_I2_RANK2(HIRS_Cld_Temp,sym%LINEAR_SCALING,Min_Tc, &
                                   Max_Tc,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_HIRS_Cld_Temp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
    !--- MJH Menzel HIRS cloud pressure for AVHRR/HIRS IFF
    if (Sds_Num_Level2_HIRS_Cld_Pres_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
        call SCALE_VECTOR_I2_RANK2(HIRS_Cld_Pres,sym%LINEAR_SCALING,Min_Pc, &
                                   Max_Pc,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_HIRS_Cld_Pres), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
    !--- MJH Menzel HIRS cloud heights for AVHRR/HIRS IFF
    if (Sds_Num_Level2_HIRS_Cld_Height_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
        call SCALE_VECTOR_I2_RANK2(HIRS_Cld_Height,sym%LINEAR_SCALING,Min_Zc, &
                                   Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_HIRS_Cld_Height), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
    !--- MJH HIRS mask
    if (Sds_Num_Level2_HIRS_Mask_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
        One_Byte_Temp = HIRS_Mask
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_HIRS_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
    !--- MJH HIRS ele collocation indices
    if (Sds_Num_Level2_HIRS_ele_index_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
        Two_Byte_Temp = HIRS_ele_index
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_HIRS_ele_index), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
    !--- MJH HIRS line collocation indices
    if (Sds_Num_Level2_HIRS_line_index_Flag == sym%YES .and. index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
        Two_Byte_Temp = HIRS_line_index
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_HIRS_line_index), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- check for and report errors
    if (Istatus /= 0) then
       print *, EXE_PROMPT, MOD_PROMPT, "Error writing to level2 file: ", Istatus
       stop
    endif
    
   endif

end subroutine WRITE_PIXEL_HDF_RECORDS


!====================================================================
! SUBROUTINE Name: CLOSE_PIXEL_HDF_FILES
!
! Function:
!   Closes appropriate Level 2 files
!
! Description:
!   This subroutine, given the flags that determine which files have 
!   been created, closes the open output files.
!
!====================================================================
subroutine CLOSE_PIXEL_HDF_FILES(Rtm_File_Flag,Level2_File_Flag)

 integer, intent(in):: Rtm_File_Flag
 integer, intent(in):: Level2_File_Flag

 integer:: Isds
 integer:: Istatus

! HDF function declarations
 integer:: sfsnatt
 integer:: sfendacc
 integer:: sfend


!------------------------------------------------------------------------
!--- close rtm file
!------------------------------------------------------------------------
  if (Rtm_File_Flag == sym%YES) then
    Istatus = 0
    Istatus = sfsnatt(Sd_Id_Rtm, "NUMBER_OF_ELEMENTS", DFNT_INT32,1,Image%Number_Of_Elements)+Istatus
    Istatus = sfsnatt(Sd_Id_Rtm, "NUMBER_OF_SCANS_LEVEL1B", DFNT_INT32,1,Image%Number_Of_Lines)+Istatus
    Istatus = sfsnatt(Sd_Id_Rtm, "NUMBER_OF_SCANS_LEVEL2", DFNT_INT32,1,Num_Scans_Level2_Hdf)+Istatus
    Istatus = sfsnatt(Sd_Id_Rtm, "PROCESSING_TIME_MINUTES", DFNT_FLOAT32,1,Orbital_Processing_Time_Minutes)+Istatus
    do Isds = 1, Num_Rtm_Sds
      if (sds_id_rtm(isds) /= 0 ) Istatus = sfendacc(Sds_Id_Rtm(Isds)) + Istatus
      sds_id_rtm(isds) = 0
    enddo
    Istatus = sfend(Sd_Id_Rtm) + Istatus
  endif

!------------------------------------------------------------------------
!--- close level2 file
!------------------------------------------------------------------------
  Istatus = 0
  Istatus = sfsnatt(Sd_Id_Level2, "NUMBER_OF_ELEMENTS", DFNT_INT32,1,Image%Number_Of_Elements)+Istatus
  Istatus = sfsnatt(Sd_Id_Level2, "NUMBER_OF_SCANS_LEVEL1B", DFNT_INT32,1,Image%Number_Of_Lines)+Istatus
  Istatus = sfsnatt(Sd_Id_Level2, "NUMBER_OF_SCANS_LEVEL2", DFNT_INT32,1,Num_Scans_Level2_Hdf)+Istatus
  Istatus = sfsnatt(Sd_Id_Level2, "PROCESSING_TIME_MINUTES", DFNT_FLOAT32,1,Orbital_Processing_Time_Minutes)+Istatus
  Istatus = sfsnatt(Sd_Id_Level2, "NONCONFIDENT_CLOUD_MASK_FRACTION", DFNT_FLOAT32,1,NONCONFIDENT_CLOUD_MASK_Fraction)+Istatus
  Istatus = sfsnatt(Sd_Id_Level2, "ACHA_SUCCESS_FRACTION", DFNT_FLOAT32,1,ACHA%Success_Fraction)+Istatus
  Istatus = sfsnatt(Sd_Id_Level2, "DCOMP_SUCCESS_FRACTION", DFNT_FLOAT32,1,DCOMP_Success_Fraction)+Istatus
  if (Level2_File_Flag == sym%YES) then
   do Isds = 1, Num_Level2_Sds
      
      if (sds_Id_Level2(Isds) /= 0 ) Istatus = sfendacc(Sds_Id_Level2(Isds)) + Istatus
      ! - do this if sds_id_level2 is saved from last file and not used here..
      ! - some compilers with ceratin flags save also without save attribute
      ! - 20141225-AW
      sds_Id_Level2(Isds) = 0
     
   enddo
   Istatus = sfend(Sd_Id_Level2) + Istatus
!--- errors are expected on close since all sds were not opened
!  if (Istatus /= 0) then
!     print *, EXE_PROMPT, MOD_PROMPT, "Error closing level2 hdf file"
!  endif
  endif

end subroutine CLOSE_PIXEL_HDF_FILES

!====================================================================
! SUBROUTINE Name: DEFINE_PIXEL_2D_SDS
!
! Function:
!   Defines a 2D SDS for the level 2 files
!
! Description:
!   This subroutine, given the inputs, creates the SDSs inside the various
!   files. 
!
!====================================================================
  subroutine DEFINE_PIXEL_2D_SDS(Sds_Id,     &
                                 Sd_Id,      &
                                 Sds_Dims,   &
                                 Sds_Chunk,  &
                                 Sds_Name,   &
                                 Sds_Standard_Name,   &
                                 Sds_Long_Name, &
                                 Sds_Type,   &
                                 scaled,     &
                                 Sds_Min,    &
                                 Sds_Max,    &
                                 Sds_Units,  &
                                 Sds_Missing,&
                                 Istatus)

    integer, intent(out):: Sds_Id
    integer, intent(in):: Sd_Id
    integer, dimension(2), intent(in):: Sds_Dims
    integer, dimension(2), intent(in):: Sds_Chunk
    integer, intent(in):: Sds_Type
    real, intent(in):: Sds_Min
    real, intent(in):: Sds_Max
    real, intent(in):: Sds_Missing
    integer(kind=int1), intent(in):: Scaled
    character(len=*), intent(in):: Sds_Name
    character(len=*), intent(in):: Sds_Standard_Name
    character(len=*), intent(in):: Sds_Long_Name
    character(len=*), intent(in):: Sds_Units
    integer, intent(out):: Istatus
    integer:: Dim_Id
    integer(kind=int4):: Scaled_Min
    integer(kind=int4):: Scaled_Max
    integer(kind=int4):: Scaled_Missing
    real(kind=real4):: Add_Offset
    real(kind=real4):: Scale_Factor

    integer:: sfcreate
    integer:: sfsnatt
    integer:: sfscatt
    integer:: sfdimid
    integer:: sfsdmname
    integer:: sfschnk
    
    real :: temp_vector_2  (2)
    integer (kind = int1)  :: temp_vector_2_int1  (2)  
    integer (kind = int2)  :: temp_vector_2_int2  (2) 
    
    integer (kind = int1), parameter  :: temp_indgen_2_int1  (2)  = [0,1]
    integer (kind = int1)  ,parameter :: temp_indgen_4_int1(4) = [0,1,2,3]
    integer (kind = int1)  , parameter :: temp_indgen_8_int1  (8) = [0,1,2,3,4,5,6,7]
    integer (kind = int1) ,parameter  :: temp_indgen_13_int1 (13) =  [0,1,2,3,4,5,6,7,8,9,10,11,12]
    integer (kind = int1)  ,parameter :: temp_indgen_14_int1 (14) =  [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    
    
    
    
    Istatus = 0 

    Sds_Id = sfcreate(Sd_Id,Sds_Name,Sds_Type,Sds_Rank_2d,Sds_Dims)

    !--- Attributes that go into all SDSs
    Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, scaled) + Istatus

    Istatus = sfscatt(Sds_Id, "units", DFNT_CHAR8, len_trim(Sds_Units), trim(Sds_Units)) + Istatus

    Istatus = sfscatt(Sds_Id, "standard_name", DFNT_CHAR8, len_trim(Sds_Standard_Name),  &
                               trim(Sds_Standard_Name)) + Istatus

    Istatus = sfscatt(Sds_Id, "long_name", DFNT_CHAR8, len_trim(Sds_Long_Name),  &
                                trim(Sds_Long_Name)) + Istatus

    if (Sds_Name /= "latitude" .and. Sds_Name /= "longitude") then 
       Istatus = sfscatt(Sds_Id, "coordinates", DFNT_CHAR8, len_trim(coordinates_string),  &
                          trim(coordinates_string)) + Istatus
    endif

    Dim_Id = sfdimid(Sds_Id, 0)
    Istatus = sfsdmname(Dim_Id,"pixel_elements_along_scan_direction") + Istatus

    Dim_Id = sfdimid(Sds_Id, 1)
    Istatus = sfsdmname(Dim_Id,"scan_lines_along_track_direction") + Istatus

    Istatus = sfschnk(Sds_Id,Sds_Chunk,Comp_Type,Comp_Prm)+Istatus

    !--- determine scaled ranges based on Sds_Type
    if (Sds_Type == DFNT_INT8) then
          Scaled_Min = One_Byte_Min
          Scaled_Max = One_Byte_Max
          Scaled_Missing = Missing_Value_Int1
    elseif (Sds_Type == DFNT_INT16) then
          Scaled_Min = Two_Byte_Min
          Scaled_Max = Two_Byte_Max
          Scaled_Missing = Missing_Value_Int2
    endif

    !--- write actual missng and actual rangel for all scaled variables
    if (Scaled /= sym%NO_SCALING) then
!     Istatus = sfsnatt(Sds_Id, "actual_missing", DFNT_FLOAT32, 1, Sds_Missing) + Istatus
      temp_vector_2 = [Sds_Min,Sds_Max]
     Istatus = sfsnatt(Sds_Id, "actual_range", DFNT_FLOAT32, 2, temp_vector_2) + Istatus
    endif

    !--- write valid_range for all scaled variables
    if (Scaled /= sym%NO_SCALING) then

     if (Sds_Type == DFNT_INT8) then
       temp_vector_2_int1 = int((/Scaled_Min,Scaled_Max/),kind=int1)
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2,  temp_vector_2_int1 ) + Istatus
     elseif (Sds_Type == DFNT_INT16) then
       temp_vector_2_int2 = int((/Scaled_Min,Scaled_Max/),kind=int2)
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2, temp_vector_2_int2) + Istatus
     elseif (Sds_Type == DFNT_FLOAT32) then
          temp_vector_2 = [Sds_Min,Sds_Max]
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2, temp_vector_2) + Istatus
     endif

    endif

    !--- write _FillValue for all except bitmasks (missing = -888)
    if (Sds_Missing /= -888.0) then

     !--- write Fill_Value
     if (Sds_Type == DFNT_INT8) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, int(Scaled_Missing,kind=int1)) + Istatus
     elseif (Sds_Type == DFNT_INT16) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, int(Scaled_Missing,kind=int2)) + Istatus
     elseif (Sds_Type == DFNT_FLOAT32) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, Sds_Missing) + Istatus
     endif

    endif

    !--- testing CF flag_meanings and flag_values attributes
    if (Sds_Name == "acha_quality") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 2, temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 143, "Processed "// &
                               " valid_Tc_retrieval "// &
                               " valid_ec_retrieval "// &
                               " valid_beta_retrieval "// &
                               " degraded_Tc_retrieval "// &
                               " degraded_ec_retrieval "// &
                               " degraded_beta_retrieval ") + Istatus
    end if

    if (Sds_Name == "acha_info") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 2, temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 196, "Cloud_Height_Attempted "// &
                               " Bias_Correction_Employed "// &
                               " Ice_Cloud_Retrieval "// &
                               " Local_Radiative_Center_Processing_Used "// &
                               " Multi_Layer_Retrieval "// &
                               " Lower_Cloud_Interpolation_Used "// &
                               " Boundary_Layer_Inversion_Assumed ") + Istatus
    end if

    if (Sds_Name == "dcomp_quality") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 2, temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 120, "Processed "// &
                               " valid_COD_retrieval "// &
                               " valid_REF_retrieval "// &
                               " degraded_COD_retrieval "// &
                               " degraded_REF_retrieval "// &
                               " convergency "// &
                               " glint ") + Istatus
    end if

    if (Sds_Name == "dcomp_info") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 2, temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 117, "info_flag_set "// &
                               " land_or_sea_mask "// &
                               " day_or_night mask "// &
                               " twilight_(65-82_solar_zenith) "// &
                               " snow "// &
                               " sea_ice "// &
                               " phase "// &
                               " thick_cloud ") + Istatus
    end if

    if (Sds_Name == "cld_opd_dcomp_qf") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cld_reff_dcomp_qf") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cld_temp_acha_qf") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cld_emiss_acha_qf") then
    
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cld_beta_acha_qf") then
      
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cloud_mask") then
    
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 47, "clear "// &
                               " probably_clear "// &
                               " probably_cloudy "// &
                               " cloudy ") + Istatus
    end if

    if (Sds_Name == "cloud_type") then
      
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 13, temp_indgen_13_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 120, "clear "// &
                               " probably_clear "// &
                               " fog "// &
                               " water "// &
                               " supercooled_water "// &
                               " mixed "// &
                               " opaque_ice "// &
                               " cirrus "// &
                               " overlapping "// &
                               " overshooting "// &
                               " unknown "// &
                               " dust "// &
                               " smoke ") + Istatus
    end if

    if (Sds_Name == "surface_type") then
      
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 14, temp_indgen_14_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 106, "water "// &
                               " evergreen_needle "// &
                               " evergreen_broad "// &
                               " deciduous_needle "// &
                               " deciduous_broad "// &
                               " mixed_forest "// &
                               " woodlands "// &
                               " wooded_grass "// &
                               " closed_shrubs "// &
                               " open_shrubs "// &
                               " grasses "// &
                               " croplands "// &
                               " bare "// &
                               " urban ") + Istatus
    end if

    if (Sds_Name == "land_class") then
      
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 8, temp_indgen_8_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 109, "ocean "// &
                               " land "// &
                               " coastline "// &
                               " shallow_inland_water "// &
                               " ephemeral_water "// &
                               " deep_inland_water "// &
                               " moderate_ocean "// &
                               " deep_ocean ") + Istatus
    end if

    if (Sds_Name == "snow_class") then
     
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 3,  temp_indgen_4_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 30, "no_snow_or_ice "// &
                               " sea_ice "// &
                               " snow ") + Istatus
    end if


    !--- Scaling Attributes that go into scaled SDSs
    if (Scaled > 0) then

     !--- determine offset and scale factor
     if (Sds_Min /= Sds_Max) then
       Scale_Factor = (Sds_Max - Sds_Min) / (Scaled_Max - Scaled_Min)
       Add_Offset =  Sds_Min - Scale_Factor*Scaled_Min
     else
       Scale_Factor = 1.0
       Add_Offset = 0.0
     endif

     !--- assume any one-byte integer is not to be scaled
     if (Sds_Missing == -128.0) then
        Scale_Factor = 1.0
        Add_Offset = 0.0
     endif

     !--- write remaining attributes
     if (Scaled == sym%LINEAR_SCALING) then
      Istatus = sfsnatt(Sds_Id, "add_offset", DFNT_FLOAT32, 1, Add_Offset) + Istatus
      Istatus = sfsnatt(Sds_Id, "scale_factor", DFNT_FLOAT32, 1, Scale_Factor) + Istatus
     endif 

     !--- write deprecated scaling attributes that supported non-linear scaling
     Istatus = sfsnatt(Sds_Id, "RANGE_MIN", DFNT_FLOAT32, 1, Sds_Min) + Istatus
     Istatus = sfsnatt(Sds_Id, "RANGE_MAX", DFNT_FLOAT32, 1, Sds_Max) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MISSING", DFNT_INT32, 1, Scaled_Missing) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MIN", DFNT_INT32, 1, Scaled_Min) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MAX", DFNT_INT32, 1, Scaled_Max) + Istatus


     if (Istatus /= 0) print *, "Error writing level2 2d sds named ", trim(Sds_Name)

   endif 

  end subroutine DEFINE_PIXEL_2D_SDS

!====================================================================
! SUBROUTINE Name: DEFINE_PIXEL_3D_SDS
!
! Function:
!   Defines a 3D SDS for the level 2 files
!
! Description:
!   This subroutine, given the inputs, creates the SDSs inside the various
!   files. 
!
!====================================================================
  subroutine DEFINE_PIXEL_3D_SDS(Sds_Id,     &
                                 Sd_Id,      &
                                 Sds_Dims,   &
                                 Sds_Chunk,  &
                                 Sds_Name,   &
                                 Sds_Standard_Name,   &
                                 Sds_Long_Name, &
                                 Sds_Type,   &
                                 scaled,     &
                                 Sds_Min,    &
                                 Sds_Max,    &
                                 Sds_Units,  &
                                 Sds_Missing,&
                                 Dim1_Name, &
                                 Dim2_Name, &
                                 Dim3_Name, &
                                 Istatus)

    integer, intent(out):: Sds_Id
    integer, intent(in):: Sd_Id
    integer, dimension(3), intent(in):: Sds_Dims
    integer, dimension(3), intent(in):: Sds_Chunk
    integer, intent(in):: Sds_Type
    real, intent(in):: Sds_Min
    real, intent(in):: Sds_Max
    real, intent(in):: Sds_Missing
    integer(kind=int1), intent(in):: scaled
    character(len=*), intent(in):: Sds_Name
    character(len=*), intent(in):: Sds_Standard_Name
    character(len=*), intent(in):: Sds_Long_Name
    character(len=*), intent(in):: Sds_Units
    character(len=*), intent(in):: Dim1_Name
    character(len=*), intent(in):: Dim2_Name
    character(len=*), intent(in):: Dim3_Name
    integer, intent(out):: Istatus
    integer:: Dim_Id
    integer(kind=int4):: Scaled_Min
    integer(kind=int4):: Scaled_Max
    integer(kind=int4):: Scaled_Missing
    real(kind=real4):: Add_Offset
    real(kind=real4):: Scale_Factor

    integer:: sfcreate
    integer:: sfsnatt
    integer:: sfscatt
    integer:: sfdimid
    integer:: sfsdmname
    integer:: sfschnk
    
    real :: temp_vector_2  (2)
    integer (kind = int1)  :: temp_vector_2_int1  (2)  
    integer (kind = int2)  :: temp_vector_2_int2  (2) 

    Istatus = 0 

    Sds_Id = sfcreate(Sd_Id,Sds_Name,Sds_Type,3,Sds_Dims)

    !--- Attributes that go into all SDSs
    Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, Scaled) + Istatus
    Istatus = sfscatt(Sds_Id, "units", DFNT_CHAR8, len_trim(Sds_Units), trim(Sds_Units)) + Istatus

    Istatus = sfscatt(Sds_Id, "standard_name", DFNT_CHAR8, len_trim(Sds_Standard_Name),  &
                               trim(Sds_Standard_Name)) + Istatus

    Istatus = sfscatt(Sds_Id, "long_name", DFNT_CHAR8, len_trim(Sds_Long_Name),  &
                                trim(Sds_Long_Name)) + Istatus
    Dim_Id = sfdimid(Sds_Id, 0)
    Istatus = sfsdmname(Dim_Id,trim(Dim1_Name)) + Istatus

    Dim_Id = sfdimid(Sds_Id, 1)
    Istatus = sfsdmname(Dim_Id,trim(Dim2_Name)) + Istatus

    Dim_Id = sfdimid(Sds_Id, 2)
    Istatus = sfsdmname(Dim_Id,trim(Dim3_Name)) + Istatus

    Istatus = sfschnk(Sds_Id,Sds_Chunk,Comp_Type,Comp_Prm)+Istatus

    !-- Range Missing written out regardless of Scaled Value
    Istatus = sfsnatt(Sds_Id, "RANGE_MISSING", DFNT_FLOAT32, 1, Sds_Missing) + Istatus
    temp_vector_2 = [Sds_Min,Sds_Max]
    Istatus = sfsnatt(Sds_Id, "actual_range", DFNT_FLOAT32, 2, temp_vector_2) + Istatus

    !--- determine scaled ranges based on Sds_Type
    if (Sds_Type == DFNT_INT8) then
          Scaled_Min = One_Byte_Min
          Scaled_Max = One_Byte_Max
          Scaled_Missing = Missing_Value_Int1
    elseif (Sds_Type == DFNT_INT16) then
          Scaled_Min = Two_Byte_Min
          Scaled_Max = Two_Byte_Max
          Scaled_Missing = Missing_Value_Int2
    endif
    
    !--- write valid_range and _FillValue for all except bitmasks (missing = -888)
    if (Sds_Missing /= -888.0) then
   
    !--- write valid range
    if (Sds_Type == DFNT_INT8) then
       temp_vector_2_int1 = int((/Scaled_Min,Scaled_Max/),kind=int1)
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2,  temp_vector_2_int1 ) + Istatus
     elseif (Sds_Type == DFNT_INT16) then
       temp_vector_2_int2 = int((/Scaled_Min,Scaled_Max/),kind=int2)
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2, temp_vector_2_int2) + Istatus
     elseif (Sds_Type == DFNT_FLOAT32) then
          temp_vector_2 = [Sds_Min,Sds_Max]
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2, temp_vector_2) + Istatus
     endif

     !--- write fill_value
     if (Sds_Type == DFNT_INT8) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, int(Scaled_Missing,kind=int1)) + Istatus
     elseif (Sds_Type == DFNT_INT16) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, int(Scaled_Missing,kind=int2)) + Istatus
     elseif (Sds_Type == DFNT_FLOAT32) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, Sds_Missing) + Istatus
     endif

    endif

    !--- Attributes that go into scaled SDSs
    if (Scaled > 0) then

     !--- determine offset and scale factor
     if (Sds_Min /= Sds_Max) then
       Scale_Factor = (Sds_Max - Sds_Min) / (Scaled_Max - Scaled_Min)
       Add_Offset =  Sds_Min - Scale_Factor*Scaled_Min
     else
       Scale_Factor = 1.0
       Add_Offset = 0.0
     endif
   
     !--- assume any one-byte integer is not to be scaled
     if (Sds_Missing == -128.0) then
        Scale_Factor = 1.0
        Add_Offset = 0.0
     endif

     !--- write remaining attributes
     Istatus = sfsnatt(Sds_Id, "RANGE_MIN", DFNT_FLOAT32, 1, Sds_Min) + Istatus
     Istatus = sfsnatt(Sds_Id, "RANGE_MAX", DFNT_FLOAT32, 1, Sds_Max) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MISSING", DFNT_INT32, 1, Scaled_Missing) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MIN", DFNT_INT32, 1, Scaled_Min) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MAX", DFNT_INT32, 1, Scaled_Max) + Istatus
     if (Scaled == sym%LINEAR_SCALING) then 
      Istatus = sfsnatt(Sds_Id, "add_offset", DFNT_FLOAT32, 1, Add_Offset) + Istatus
      Istatus = sfsnatt(Sds_Id, "scale_factor", DFNT_FLOAT32, 1, Scale_Factor) + Istatus
     endif

   endif

   if (Istatus /= 0) print *, "Error writing level2 3d sds named ", trim(Sds_Name)

  end subroutine DEFINE_PIXEL_3D_SDS

  !============================================================================
  !
  !============================================================================
  subroutine WRITE_ALGORITHM_ATTRIBUTES()

    integer:: sfscatt
    integer:: istatus
    
    istatus = 0
    istatus = sfscatt(Sd_Id_Level2, "CLOUD_MASK_VERSION", DFNT_CHAR8, len_trim(Cloud_Mask_Version), trim(Cloud_Mask_Version))+istatus
    istatus = sfscatt(Sd_Id_Level2, "CLOUD_MASK_THRESHOLDS_VERSION", DFNT_CHAR8,  &
                 len_trim(Cloud_Mask_Thresholds_Version), trim(Cloud_Mask_Thresholds_Version))+istatus 
    istatus = sfscatt(Sd_Id_Level2, "CLOUD_TYPE_VERSION", DFNT_CHAR8, len_trim(Cloud_Type_Version), trim(Cloud_Type_Version))+istatus 
    istatus = sfscatt(Sd_Id_Level2, "ACHA_VERSION", DFNT_CHAR8, len_trim(ACHA_Version), trim(ACHA_Version))+istatus 
    istatus = sfscatt(Sd_Id_Level2, "DCOMP_VERSION", DFNT_CHAR8, len_trim(dcomp_version), trim(dcomp_version))+istatus

  end subroutine WRITE_ALGORITHM_ATTRIBUTES

!============================================================================

end module LEVEL2_ROUTINES


