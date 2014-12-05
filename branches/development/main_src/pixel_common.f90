! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: pixel_common.f90 (src)
!       PIXEL_COMMON (program)
!
! PURPOSE: This module houses routines defining the global set of variables and
!          arrays used in CLAVR-x
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
! Public routines in this module:
!
! CREATE_PIXEL_ARRAYS - allocate memory for pixel level arrays
! DESTROY_PIXEL_ARRAYS - deallocate memory for pixel level arrays
! RESET_PIXEL_ARRAYS_TO_MISSING - set pixel arrays to missing
!
! File I/O: None
!
!  CLAVR-x uses MODIS channel numbers for all sensors
!
!  Modis Band   Avhrr Band   Abi Band   VIIRS   Wavelength
!     01            1          2         M5       0.659
!     02            2          3         M7       0.865
!     03            -          -         M3       0.470 
!     04            -          -         M4       0.555
!     05            -          -         M8       1.240
!     06            3a         5         M10      1.640
!     07            -          -         M11      2.130
!     08            -          -         M1       0.415
!     09            -          -         M2       0.443
!     10            -          -         -        0.490
!     11            -          -         -        0.531
!     12            -          -         -        0.565
!     13            -          -         -        0.653
!     14            -          -         -        0.681
!     15            -          -         M6       0.750
!     16            -          -         -        0.865
!     17            -          -         -        0.905
!     18            -          -         -        0.936
!     19            -          -         -        0.940
!     20            3b         7         M12      3.750
!     21            -          7         -        3.959
!     22            -          7         M13      3.959
!     23            -          -         -        4.050
!     24            -          -         -        4.465
!     25            -          -         -        4.515
!     26            -          -         M9       1.375
!     27            -          9         -        6.715
!     28            -         10         -        7.325
!     29            -         11         M14      8.550
!     30            -         12         -        9.730
!     31            4         14         M15     11.030
!     32            5         15         M16     12.020
!     33            -         16         -       13.335
!     34            -          -         -       13.635
!     35            -          -         -       13.935
!     36            -          -         -       14.235
!     37            -          -         I1       0.640
!     38            -          -         I2       0.865
!     39            -          -         I3       1.610
!     40            -          -         I4       3.740
!     41            -          -         I5      11.450
!     42            -          -         DNB      0.700
!
!  Description of variables in "ch" structure:
!
!  Rad_Toa = Observed Top of Atmosphere Radiance
!  Bt_Toa = Observed Top of Atmosphere Brightness Temperature
!  Rad_Toa_Clear = Simulated Top of Atmosphere Radiance under clear-sky
!  Rad_Atm = Simulated Radiance at Toa from atmospheric cloud-free emission
!  Trans_Atm = Simulated Transmission from Surface to Toa for cloud-free atmosphere
!              along viewing zenith angle path
!  Trans_Atm_Total = Simulated Transmission from Toa to Surface to Toa for cloud-free atmosphere
!                      along solar and viewing zenith angle path
!  Bt_Toa_Clear = Simulated Toa Brightness Temperature under clear-sky
!  Bt_Sfc = Observed Brightness Temperature adjusted as if measured at surface level
!  Rad_Sfc = Observed Radiance adjusted as if measured at surface level
!  Ref_Toa = Top of Atmosphere Reflectance
!  Ref_Toa_Unnorm = Top of Atmosphere Reflectance not normalized by cos(solzen)
!  Ref_Sfc = Observed Reflectance adjusted as if measured at surface level
!  Sfc_Ref_White_Sky - surface reflectance under diffuse illumination
!  Emiss_Tropo - emissity of cloud placed at Tropopause needed to match toa  radiance
!  Bt_Toa_Clear = Simulated Toa Reflectance under clear-sky
!  Unc = Uncertainty Flag (relevant only to MODIS)
!--------------------------------------------------------------------------------------
module PIXEL_COMMON

  use CONSTANTS
  implicit none
  private
  public:: CREATE_PIXEL_ARRAYS, &
           DESTROY_PIXEL_ARRAYS, &
           RESET_PIXEL_ARRAYS_TO_MISSING

  private:: CREATE_GEO_ARRAYS, RESET_GEO_ARRAYS, DESTROY_GEO_ARRAYS
  private:: CREATE_GEO_ANCHOR_ARRAYS, RESET_GEO_ANCHOR_ARRAYS, DESTROY_GEO_ANCHOR_ARRAYS
  private:: CREATE_NWP_PIX_ARRAYS, RESET_NWP_PIX_ARRAYS, DESTROY_NWP_PIX_ARRAYS
  private:: CREATE_REF_CHANNEL_ARRAYS, RESET_REF_CHANNEL_ARRAYS, DESTROY_REF_CHANNEL_ARRAYS
  private:: CREATE_THERM_CHANNEL_ARRAYS, RESET_THERM_CHANNEL_ARRAYS, DESTROY_THERM_CHANNEL_ARRAYS
  private:: CREATE_EXTRA_CHANNEL_ARRAYS, RESET_EXTRA_CHANNEL_ARRAYS, DESTROY_EXTRA_CHANNEL_ARRAYS
  private:: CREATE_LUNAR_ARRAYS, RESET_LUNAR_ARRAYS, DESTROY_LUNAR_ARRAYS
  private:: CREATE_BTD_ARRAYS, RESET_BTD_ARRAYS, DESTROY_BTD_ARRAYS
  private:: CREATE_SURFACE_ARRAYS, RESET_SURFACE_ARRAYS, DESTROY_SURFACE_ARRAYS
  private:: CREATE_ACHA_ARRAYS, RESET_ACHA_ARRAYS, DESTROY_ACHA_ARRAYS
  private:: CREATE_DCOMP_ARRAYS, RESET_DCOMP_ARRAYS, DESTROY_DCOMP_ARRAYS
  private:: CREATE_NLCOMP_ARRAYS, RESET_NLCOMP_ARRAYS, DESTROY_NLCOMP_ARRAYS
  private:: CREATE_SASRAB_ARRAYS, RESET_SASRAB_ARRAYS, DESTROY_SASRAB_ARRAYS
  private:: CREATE_OLR_ARRAYS, RESET_OLR_ARRAYS, DESTROY_OLR_ARRAYS
  private:: CREATE_AEROSOL_ARRAYS, RESET_AEROSOL_ARRAYS, DESTROY_AEROSOL_ARRAYS
  private:: CREATE_CLOUD_MASK_ARRAYS, RESET_CLOUD_MASK_ARRAYS, DESTROY_CLOUD_MASK_ARRAYS
  private:: CREATE_CLOUD_TYPE_ARRAYS, RESET_CLOUD_TYPE_ARRAYS, DESTROY_CLOUD_TYPE_ARRAYS
  private:: CREATE_DIAGNOSTIC_ARRAYS, RESET_DIAGNOSTIC_ARRAYS, DESTROY_DIAGNOSTIC_ARRAYS
  private:: CREATE_SFC_PROD_ARRAYS, RESET_SFC_PROD_ARRAYS, DESTROY_SFC_PROD_ARRAYS
  private:: CREATE_CLOUD_PROD_ARRAYS, RESET_CLOUD_PROD_ARRAYS, DESTROY_CLOUD_PROD_ARRAYS

  !--- arrays to keep track of files written to temporary directory, max 100 assumed
  integer, public, save:: Number_Of_Temporary_Files
  character(len=200),dimension(100), public, save:: Temporary_File_Name

  !---------------------------------------------------------------------------------
  ! CLAVR-x file list variables
  !---------------------------------------------------------------------------------
  type :: observations
    real, dimension(:,:), allocatable:: Rad_Toa
    real, dimension(:,:), allocatable:: Bt_Toa
    real, dimension(:,:), allocatable:: Rad_Toa_Clear
    real, dimension(:,:), allocatable:: Rad_Atm
    real, dimension(:,:), allocatable:: Trans_Atm
    real, dimension(:,:), allocatable:: Trans_Atm_Total
    real, dimension(:,:), allocatable:: Bt_Toa_Clear
    real, dimension(:,:), allocatable:: Bt_Sfc
    real, dimension(:,:), allocatable:: Rad_Sfc
    real, dimension(:,:), allocatable:: Ref_Toa
    real, dimension(:,:), allocatable:: Ref_Toa_Unnorm
    real, dimension(:,:), allocatable:: Ref_Sfc
    real, dimension(:,:), allocatable:: Ref_Toa_Clear
    real, dimension(:,:), allocatable:: Ref_Lunar_Toa
    real, dimension(:,:), allocatable:: Ref_Lunar_Toa_Clear
    real, dimension(:,:), allocatable:: Ref_Lunar_Sfc
    real, dimension(:,:), allocatable:: Sfc_Emiss
    real, dimension(:,:), allocatable:: Emiss_Tropo
    real, dimension(:,:), allocatable:: Sfc_Ref_White_Sky
    integer (kind=int1), dimension(:,:), allocatable:: Unc
  end type observations

  type(observations), dimension(42), public, target :: ch

  integer,public, save:: Bx_File_Flag
  integer,public, save:: Cmr_File_Flag
  integer,public, save:: Cloud_Mask_Aux_Flag
  integer,public, save:: Cloud_Mask_Aux_Read_Flag
  integer,public, save:: Cloud_Mask_Bayesian_Flag
  integer,public, save:: Level3_Flag
  integer,public, save:: Ref_cal_1b 
  integer,public, save:: Therm_cal_1b
  integer,public, save:: Nav_Flag       !0=level1b,1=clevernav,2=reposnx
  integer,public, save:: Nav_File_Flag  !yes, write out a navigation file
  integer,public, save:: Obs_File_Flag
  integer,public, save:: Geo_File_Flag
  integer,public, save:: Sst_File_Flag
  integer,public, save:: Cld_File_Flag
  integer,public, save:: Rtm_File_Flag
  integer,public, save:: Ash_File_Flag
  integer,public, save:: Level2_File_Flag
  integer,public, save:: Subset_Pixel_Hdf_Flag
  integer,public, save:: Use_Sst_Anal
  integer,public, save:: Use_Sst_Anal_Default
  integer,public, save:: L1b_Gzip
  integer,public, save:: L1b_Bzip2
  integer,public, save:: Diag_Flag
  integer,public, save:: Sst_Anal_Opt
  integer,public, save:: Data_Comp_Flag
  integer,public, save:: Use_seebor
  integer,public, save:: Read_Volcano_Mask
  integer,public, save:: Read_Land_Mask
  integer,public, save:: Read_Coast_Mask
  integer,public, save:: Read_Surface_Elevation
  integer,public, save:: Read_Hires_Sfc_Type
  integer,public, save:: Read_Snow_Mask
  integer,public, save:: Read_GlobSnow_Mask
  integer,public, save:: Read_Dark_Comp
  integer,public, save:: Machine_Byte_Ordering
  integer,public, save:: LRC_Flag  !local radiative center flag
  integer,public, save:: Process_Undetected_Cloud_Flag
  integer,public, save:: DCOMP_Mode
  integer,public, save:: ACHA_Mode
  integer,public, save:: Cld_Flag
  integer,public, save:: Blank_Flag
  integer,public, save:: Aer_Flag
  integer,public, save:: Ash_Flag
  integer,public, save:: Erb_Flag
  integer,public, save:: Nwp_Flag
  integer,public, save:: Modis_Clr_Alb_Flag
  integer,public, save:: Asc_Flag_Diag
  integer,public, save:: Rtm_Flag
  integer,public, save:: Prob_Clear_Res_Flag
  integer,public, save:: Smooth_Nwp_Flag

  !---------------------------------------------------------------------------------
  ! Flags Computed within CLAVR-x that describe the sensor data
  !---------------------------------------------------------------------------------
  integer,public, save:: AVHRR_GAC_Flag
  integer,public, save:: AVHRR_KLM_Flag
  integer,public, save:: AVHRR_AAPP_Flag
  integer,public, save:: Viirs_Flag
  integer,public, save:: Iff_Viirs_Flag
  integer,public, save:: Iff_Modis_Flag
  integer,public, save:: Seviri_Flag
  integer,public, save:: Mtsat_Flag
  integer,public, save:: FY2_Flag
  integer,public, save:: COMS_Flag
  integer,public, save:: Goes_Flag
  integer,public, save:: Goes_Mop_Flag
  integer,public, save:: Goes_Sndr_Flag
  integer,public, save:: Goes_1km_Flag
  integer,public, save:: Modis_Flag
  integer,public, save:: Modis_5km_Flag
  integer,public, save:: Modis_1km_Flag
  integer,public, save:: Modis_CSPP_Flag
  integer,public, save:: Modis_Aqua_Flag
  integer,public, save:: Modis_Aqua_Mac_Flag
  integer,public, save:: Modis_Terra_Flag
  integer,public, save:: Avhrr_Flag
  integer,public, save:: Avhrr_1_Flag

  !---------------------------------------------------------------------------------
  ! Internal Flags to communicate ancillary data information
  !---------------------------------------------------------------------------------
  integer,public, save:: Failed_Hires_Snow_Mask_Flag
  integer,public, save:: Failed_Glob_Snow_Mask_Flag
  integer,public, save:: Output_Scaled_Reflectances
  integer,public, save:: Ncdc_Level2_Flag


  !---------------------------------------------------------------------------------
  ! Default Algorithm Modes - (maybe move to user options)
  !---------------------------------------------------------------------------------
  integer,public,parameter:: ACHA_Mode_Default_Avhrr = 1
  integer,public,parameter:: ACHA_Mode_Default_Avhrr1 = 0
  integer,public,parameter:: ACHA_Mode_Default_Goes_IL = 5
  integer,public,parameter:: ACHA_Mode_Default_Goes_MP = 6
  integer,public,parameter:: ACHA_Mode_Default_Modis = 3
  integer,public,parameter:: ACHA_Mode_Default_VIIRS = 4
  integer,public,parameter:: ACHA_Mode_Default_MTSAT = 5

  !---------------------------------------------------------------------------------
  ! variables that are computed to serve as attributes in the output files
  !---------------------------------------------------------------------------------
  real(kind=real4), public, save:: Orbital_Processing_Time_Minutes
  real(kind=real4), public, save:: DCOMP_Success_Fraction
  real(kind=real4), public, save:: ACHA_Success_Fraction
  real(kind=real4), public, save:: ACHA_Processed_Count
  real(kind=real4), public, save:: ACHA_Valid_Count
  real(kind=real4), public, save:: DCOMP_Processed_Count
  real(kind=real4), public, save:: DCOMP_Valid_Count
  real(kind=real4), public, save:: Nonconfident_Cloud_Mask_Fraction
  real(kind=real4), public, save:: Nonconfident_Cloud_Mask_Count
  real(kind=real4), public, save:: Cloud_Mask_Count

  integer(kind=int4), public, save:: Byte_Swap_1b
  integer(kind=int4), public, save:: Lun_Level1b

  !---------------------------------------------------------------------------------
  ! CLAVR-x file list variables
  !---------------------------------------------------------------------------------
  character(len=128),public,save::File_1b
  character(len=128),public,save::File_1bx
  character(len=128),public,save::File_cmr
  character(len=128),public,save::File_sst
  character(len=128),public,save::File_1b_root
  character(len=128),public,save:: Ancil_Data_Dir
  character(len=128),public,save:: Gfs_Data_Dir
  character(len=128),public,save:: Ncep_Data_Dir
  character(len=128),public,save:: Cfsr_Data_Dir
  character(len=128),public,save:: Oisst_Data_Dir
  character(len=256),public,save:: Snow_Data_Dir
  character(len=256),public,save:: GlobSnow_Data_Dir
  character(len=128),public,save:: Dark_Comp_Data_Dir
  character(len=128),public,save:: Temporary_Data_Dir
  character(len=128),public,save:: Dir_1b
  character(len=128),public,save:: Dir_1bx
  character(len=128),public,save:: File_nav
  character(len=128),public,save:: Instr_Const_File
  character(len=128),public,save:: Algo_Const_File
  character(len=128),public,save:: Dir_cmr
  character(len=128),public,save:: Dir_sst
  character(len=128),public,save:: Dir_Level3
  character(len=128),public,save:: Dir_nav_in
  character(len=128),public,save:: Dir_nav_out
  character(len=128),public,save:: Dir_cld
  character(len=128),public,save:: Dir_obs
  character(len=128),public,save:: Dir_geo
  character(len=128),public,save:: Dir_Rtm
  character(len=128),public,save:: Dir_ash
  character(len=128),public,save:: Dir_Level2
  character(len=128),public,save:: Bayesian_Cloud_Mask_Name
  character(len=128),public,save:: Modis_Geo_Name
  character(len=128),public,save:: Modis_Cloud_Mask_Name
  character(len=128),public,save:: Dark_Composite_Name

  !----- IFF data files
  character(len=128),public,save:: IFF_File

  real(kind=real4), public, save:: Dlat
  real(kind=real4), public, save:: Lat_Min_Limit
  real(kind=real4), public, save:: Lat_Max_Limit
  real(kind=real4), public, save:: Solzen_Min_Limit
  real(kind=real4), public, save:: Solzen_Max_Limit
  real(kind=real4), public, save:: Timerr_Seconds

  !------------------------------------------------------------------
  !--- variables pertaining to scanline size
  !------------------------------------------------------------------
  integer, public, save:: Num_Scans_Per_Segment
  integer, public, save:: Num_Scans_Read
  integer, public, save:: Num_Segments
  integer, public, save:: Line_Idx_Min_Segment
  integer, public, save:: Line_Idx_Max_Segment
  integer(kind=int4), public, save:: Num_Pix
  integer(kind=int4), public, save:: l1b_Rec_Length
  integer(kind=int4), public, save:: Num_Anchors
  integer, public, save:: level3_format
  real(kind=real4), public, save:: dLat_hist2d

  !--- channel on/off flags
  integer(kind=int1), dimension(Nchan_Clavrx), public, target:: Chan_On_Flag_Default
  integer(kind=int1), dimension(:,:), allocatable, public:: Chan_On_Flag

  !------- pixel array declarations
  integer (kind=int1), dimension(:), allocatable, public:: Ascend
  integer (kind=int1), dimension(:), allocatable, public:: Bad_Scan_Flag
  integer (kind=int1), dimension(:), allocatable, public:: Ch3a_On_AVHRR

  integer (kind=int4), dimension(:), allocatable, public:: Scan_Number
  integer (kind=int4), dimension(:), allocatable, public:: Scan_Time
  integer (kind=int2), dimension(:), allocatable, public:: Scan_Day
  integer (kind=int2), dimension(:), allocatable, public:: Scan_Year
  real (kind=real4), dimension(:), allocatable, public:: Utc_Scan_Time_Hours
  real (kind=real4), dimension(:,:), allocatable, public:: Pixel_Local_Time_Hours
  real (kind=real4), dimension(:,:), allocatable, public:: Pixel_Time

  real (kind=real4), dimension(:,:), allocatable,save,public:: Solzen_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Lat_Anchor_1b
  real (kind=real4), dimension(:,:), allocatable,save,public:: Lon_Anchor_1b
  real (kind=real4), dimension(:,:), allocatable,save,public:: Satzen_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Relaz_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Solaz_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Sataz_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Glintzen_Anchor
  real (kind=real4), dimension(:,:), allocatable,save,public:: Scatangle_Anchor

  !---------------------------------------------------------------------------
  !--- parameters from avhrr Header
  !---------------------------------------------------------------------------
  integer(kind=int4), public, save:: Num_Scans
  integer(kind=int4), public, save:: Num_Scans_Level2_Hdf
  integer(kind=int2), public, save:: Sc_Id,Ver_1b,Data_Type, Num_Loc, Tip_Parity, Aux_Sync,  &
                                     Ramp_Auto_Cal,Start_Year_Prev, Start_Day_Prev, &
                                     Month,Month_Prev,Day_of_Month,Ileap
  character(len=6), public, save:: Sc_Id_Char
  character(len=7),public,save:: Proc_Block_Id
  integer(kind=int4), save, public:: Orbit_Number
  integer(kind=int4), save, public:: Start_Time
  integer(kind=int4), save, public:: End_Time
  integer(kind=int2), save, public:: Start_Year
  integer(kind=int2), save, public:: End_Year
  integer(kind=int2), save, public:: Start_Day
  integer(kind=int2), save, public:: End_Day

  !--- satellite and sensor descriptors
  integer(kind=int4), public, save:: Sc_Id_Internal
  integer(kind=int4), public, save:: Sc_Id_WMO
  integer(kind=int4), public, save:: Sc_Id_WMO_Prev
  character(len=20),public,save:: Sensor_Name_Attribute
  character(len=20),public,save:: Platform_Name_Attribute

  !--- instrument counts
  integer (kind=int2), dimension(:,:), allocatable, public,save:: Ch1_Counts
  integer (kind=int2), dimension(:,:), allocatable, public,save:: Ch2_Counts
  integer (kind=int2), dimension(:,:), allocatable, public,save:: Ch6_Counts
  real (kind=real4), dimension(:,:), allocatable, public,save:: Ch20_Counts_Filtered

  !--- calibrated observations
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChI2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChI3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_ChI4
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_ChI5
  real (kind=real4), dimension(:,:), allocatable, public, save,target:: Lunzen
  real (kind=real4), dimension(:,:), allocatable, public, save,target:: Lunaz
  real (kind=real4), dimension(:,:), allocatable, public, save,target:: LunRelaz
  real (kind=real4), public, save,target:: Moon_Phase_Angle
  real (kind=real4), public, save,target:: Moon_Illum_Frac

  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Mean_3x3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Max_3x3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Min_3x3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Std_3x3

  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Max_ChI5
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Min_ChI5
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Mean_ChI5
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Max_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Min_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Mean_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Uni_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Max_ChI2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Min_ChI2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Mean_ChI2



  real (kind=real4), dimension(:,:), allocatable, public, save:: Rad_Ch20_ems
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch31_LRC
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch31_Max_LRC
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch31_Std_LRC
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Solzen
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Lat
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Lon
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Lat_Pc
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Lon_Pc
  real (kind=real4), dimension(:,:), allocatable, public, save:: Lat_1b
  real (kind=real4), dimension(:,:), allocatable, public, save:: Lon_1b
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Glintzen
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Glintzen_Lunar
  real (kind=real4), dimension(:,:), allocatable, public, save,target:: Satzen
  real (kind=real4), dimension(:,:), allocatable, public, save,target:: Relaz
  real (kind=real4), dimension(:,:), allocatable, public, save,target:: Solaz
 
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Sataz
  real (kind=real4), dimension(:,:), allocatable, public, save:: Seczen
  real (kind=real4), dimension(:,:), allocatable, public, save:: Ems_Ch20
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zsfc
  real (kind=real4), dimension(:,:), allocatable, public, save:: Zsfc_Hires
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Coszen
  real (kind=real4), dimension(:,:), allocatable, public, save,target:: CosSolzen
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Scatangle
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Scatangle_Lunar
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Airmass

  real (kind=real4), dimension(:,:), allocatable, public, save:: Sst_Anal
  real (kind=real4), dimension(:,:), allocatable, public, save:: Sst_Anal_Err
  real (kind=real4), dimension(:,:), allocatable, public, save:: Sst_Anal_Uni
  real (kind=real4), dimension(:,:), allocatable, public, save:: Sst_Anal_Cice
  real (kind=real4), dimension(:,:), allocatable, public, save:: Tsfc_Retrieved
  real (kind=real4), dimension(:,:), allocatable, public, save:: Trad_Retrieved

  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Temp_Pix_Array
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Temp_Pix_Array_2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Diag_Pix_Array_1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Diag_Pix_Array_2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Diag_Pix_Array_3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Missing_Pixel_Array_Real4

  !uniformity metrics
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch31_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch31_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch31_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch31_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Ch1_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Ch1_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Ch1_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Ch1_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch20_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch20_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch20_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch20_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch27_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch27_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch27_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch27_Std_3x3

  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch27_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch27_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch27_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch27_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch29_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch29_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch29_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch29_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch32_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch32_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch32_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch32_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch33_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch33_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Btd_Ch31_Ch33_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch33_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ems_Ch20_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ems_Ch20_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ems_Ch20_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ems_Ch20_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Zsfc_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Zsfc_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Zsfc_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Zsfc_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Ems_Ch20_Median_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ems_Ch20_Std_Median_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Bt_Ch20_Median_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Bt_Ch20_Std_Median_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ref_Ch1_Sfc_White_Sky_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ref_Ch1_Sfc_White_Sky_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ref_Ch1_Sfc_White_Sky_Min_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ref_Ch1_Sfc_White_Sky_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch32_Bt_Ch31_Max_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Cloud_Fraction_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Cloud_Fraction_Uncer_3x3

  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Covar_Ch27_Ch31_5x5

  integer(kind=int4), dimension(:,:), allocatable, public, save:: Elem_Idx_Max_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Line_Idx_Max_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Elem_Idx_Min_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Line_Idx_Min_Bt_Ch31_3x3

  integer (kind=int1), dimension(:,:), allocatable, public:: Dcc_Mask !deep conv cld Mask

  real (kind=real4), dimension(:,:), allocatable, public:: Sst_Unmasked   !sst used in cld Mask
  real (kind=real4), dimension(:,:), allocatable, public:: Sst_Masked !sst where non-clear ocean is Masked
  real (kind=real4), dimension(:,:), allocatable, public:: Ndvi_Toa
  real (kind=real4), dimension(:,:), allocatable, public:: Ndsi_Toa
  real (kind=real4), dimension(:,:), allocatable, public, target:: Btd_Ch31_Ch32
  real (kind=real4), dimension(:,:), allocatable, public:: Btd_Ch20_Ch31
  real (kind=real4), dimension(:,:), allocatable, public:: Btd_Ch20_Ch32

  integer(kind=int1), dimension(:,:), allocatable, public, target:: Land
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Land_Modified
  integer(kind=int1), dimension(:,:), allocatable, public:: Coast
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Desert_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Sfc_Type
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Snow
  integer(kind=int1), dimension(:,:), allocatable, public:: Snow_Hires
  integer(kind=int1), dimension(:,:), allocatable, public:: Snow_Glob
  integer(kind=int1), dimension(:,:), allocatable, public:: Dust
  integer(kind=int1), dimension(:,:), allocatable, public:: Smoke
  integer(kind=int1), dimension(:,:), allocatable, public:: Fire
  integer(kind=int1), dimension(:,:), allocatable, public:: Solar_Contamination_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Bad_Pixel_Mask
  integer(kind=int1), dimension(:,:), allocatable, public:: Ch6_On_Pixel_Mask
  integer(kind=int1), dimension(:,:), allocatable, public:: Volcano_Mask
  integer(kind=int1), dimension(:,:), allocatable, public:: Coast_Mask_Nwp
  integer(kind=int1), dimension(:,:), allocatable, public:: Space_Mask
  integer(kind=int4), dimension(:,:), allocatable, public:: Sfc_Level_Rtm_Pixel
  real, public:: Segment_Valid_Fraction

  !--- viirs arrays
  integer(kind=int4), dimension(:,:), allocatable, public:: Gap_Pixel_Mask_Pattern
  integer(kind=int4), dimension(:,:), allocatable, public:: Gap_Line_Idx_Pattern
  integer(kind=int1), dimension(:,:), allocatable, public:: Gap_Pixel_Mask
  integer(kind=int4), dimension(:,:), allocatable, public:: Gap_Line_Idx
  integer(kind=int1), dimension(:,:), allocatable, public:: IFF_Gap_Mask


  !--- Mask arrays
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Land_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Coast_Mask 
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Glint_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Glint_Mask_Lunar
  integer(kind=int1), dimension(:,:), allocatable, public:: Bayes_Mask_Sfc_Type_Global

  !--- cloud Mask arrays
  integer (kind=int1), dimension(:,:,:), allocatable, public, save:: Cld_Test_Vector_Packed
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Cld_Mask_Aux
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Mask
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Adj_Pix_Cld_Mask
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Cld_Mask_Qf
  real (kind=int4),dimension(:,:),allocatable, public, save:: &
                                                       Posterior_Cld_Probability

  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Type
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Cld_Type_Aux
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Phase
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Cld_Phase_Aux
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Cirrus_Quality

  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Temp_Mask

!--- pixel level cloud props

  !-- ACHA cloud algorithm results
  integer(kind=int1), dimension(:,:), allocatable, public, target:: ACHA_Processing_Order_Global !changed to target
  integer(kind=int1), dimension(:,:), allocatable, public, target:: ACHA_Inversion_Flag_Global
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tc_ACHA
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ec_ACHA
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_ACHA
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_ACHA
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_Top_ACHA
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_Base_ACHA
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Beta_ACHA
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tau_ACHA
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Reff_ACHA
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tc_ACHA_Uncertainty
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ec_ACHA_Uncertainty
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Beta_ACHA_Uncertainty
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_ACHA_Uncertainty
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_ACHA_Uncertainty
  integer (kind=int1), dimension(:,:), allocatable, public, save, target:: ACHA_Quality_Flag
  integer (kind=int1), dimension(:,:), allocatable, public, save, target:: Cld_Layer_ACHA
  integer (kind=int1), dimension(:,:), allocatable, public, save, target:: Meta_Data_ACHA
  integer (kind=int1), dimension(:,:,:), allocatable, public, save, target:: ACHA_OE_Quality_Flags

     !--- Converted pressure to altitude heights.
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Alt_Acha

     !--- h2o heights
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tc_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zclr_H2O_Peak

     !-- DCOMP cloud algorithm results
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_ap
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: vis_Ref_fm
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Iwp_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Iwp_Tau_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Lwp_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Ice_Layer_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Water_Layer_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Scwater_Layer_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Rain_Rate_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: H_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: N_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_Cost
     real (kind=real4), dimension(:,:), allocatable, public, target, save:: Reff_DCOMP_Cost
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_Qf
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_Qf
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: DCOMP_Quality_Flag
     integer (kind=int2), dimension(:,:), allocatable, public,target, save:: DCOMP_Info_Flag
     integer (kind=int1), dimension(:,:), allocatable, public, save, target:: ACHA_Packed_Quality_Flags
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: ACHA_Packed_Meta_Data_Flags

     !-- Nlcomp cloud algorithm results
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_Nlcomp
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_Nlcomp
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_Nlcomp_Cost
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_Nlcomp_Cost
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: Nlcomp_Quality_Flag
     integer (kind=int2), dimension(:,:), allocatable, public,target, save:: Nlcomp_Info_Flag

     !--- pixel level aerosol props
     real (kind=real4), dimension(:,:), allocatable, public, save:: Aot1
     real (kind=real4), dimension(:,:), allocatable, public, save:: Aot2
     real (kind=real4), dimension(:,:), allocatable, public, save:: Aot3a
     integer (kind=int1), dimension(:,:), allocatable, public, save:: Aot_Qf

     !--- non-cloud properties
     integer (kind=int1), dimension(:,:), allocatable, public, save:: Ndvi_Qf
     integer (kind=int1), dimension(:,:), allocatable, public, save:: Tsfc_Qf

     !--- pixel level radiative flux props from DCOMP
     real (kind=real4), dimension(:,:), allocatable, public, save:: Olr
     integer (kind=int1), dimension(:,:), allocatable, public, save:: Olr_Qf
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_063um_Albedo
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_063um_Spherical_Albedo
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_063um_Transmission_View
     real (kind=real4), dimension(:,:), allocatable, public, save,target:: Cloud_063um_Transmission_Solar
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_DCOMP_Diffuse

     !--- SASRAB output
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_All_Sky
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_All_Sky_Diffuse
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_Clear_Sky
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_Cld_Opd
     real (kind=real4), dimension(:,:), allocatable, public, save:: Insolation_Aer_Opd
 
     !---- DCOMP params	 
     real (kind=real4), dimension(:,:), allocatable, public, save:: Cost_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public, save:: Error_Cov_Matrix_Cod
     real (kind=real4), dimension(:,:), allocatable, public, save:: Error_Cov_Matrix_Ref
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_2
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_3
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_4
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Wv1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Wv2
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Toc_Rfl1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Toc_Rfl2
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Virt_Alb1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Virt_Alb2
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Toc_Rfl_Unc1
     real (kind=real4), dimension(:,:), allocatable, public, save:: DCOMP_Diag_Toc_Rfl_Unc2

     INTEGER, public, save::DCOMP_VIS_CHN
     INTEGER, public, save::DCOMP_IR_CHN
     character(200), public, save:: DCOMP_version_string

     !--- pixel level radiative flux props
     real (kind=real4), dimension(:,:), allocatable, public, save:: Rsr
     integer (kind=int1), dimension(:,:), allocatable, public, save:: Rsr_Qf

     !-- pixel level atmospheric reflectance components
     real (kind=real4), dimension(:,:), allocatable, public, target:: Ref_Ch1_Clear_Min_3x3
     real (kind=real4), dimension(:,:), allocatable, public, target:: Ref_Ch1_Clear_Max_3x3
     real (kind=real4), dimension(:,:), allocatable, public, target:: Ref_Ch1_Clear_Mean_3x3
     real (kind=real4), dimension(:,:), allocatable, public, target:: Ref_Ch1_Clear_Std_3x3
     real (kind=real4), dimension(:,:), allocatable, public, target:: Ref_Ch1_Dark_Composite
     real (kind=real4), dimension(:,:), allocatable, public:: Ref_Ch20_LRC
     real (kind=real4), dimension(:,:), allocatable, public:: Ndvi_Sfc
     real (kind=real4), dimension(:,:), allocatable, public:: Bt_Ch31_Sfc
     real (kind=real4), dimension(:,:), allocatable, public:: Rad_Ch20_ems_Sfc
     real (kind=real4), dimension(:,:), allocatable, public:: Ems_Ch20_Sfc

  !--- integers values for output to files
  integer(kind=int1), dimension(:,:), public,save,allocatable, target:: One_Byte_Temp
  integer(kind=int2), dimension(:,:), public,save,allocatable, target:: Two_Byte_Temp

  !---volcanic ash parameters
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_Probability
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_Probability_IR
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_mass_Loading
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_height
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_pressure
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_temperature
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_Ir_optical_depth
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_vis_optical_depth
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_effective_radius
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_emiss
  integer(kind=int1), allocatable, dimension(:,:), public, save:: ash_temperature_Qf
  integer(kind=int1), allocatable, dimension(:,:), public, save:: ash_Emiss_Qf
  integer(kind=int1), allocatable, dimension(:,:), public, save:: ash_beta1112_Qf
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_Emiss_tot
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_beta1112
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_beta1112_tot
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_beta1112_Opaque
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_beta1112_tot_multi_Low
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_beta1112_tot_multi_mid
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_generic1
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_generic2
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_tpw
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_emissivity_Error
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_temperature_Error
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_height_Error
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_beta1112_Error
  integer(kind=int1), allocatable, dimension(:,:), public, save:: ash_Qf
  integer(kind=int1), allocatable, dimension(:,:), public, save:: ash_LRC_mask
  integer(kind=int4), allocatable, dimension(:,:), public, save:: ash_x_LRC
  integer(kind=int4), allocatable, dimension(:,:), public, save:: ash_y_LRC
  real(kind=real4), allocatable, dimension(:,:), public, save:: ash_ref_rat_nir_LRC
  integer(kind=int1), pointer, dimension(:,:), public, save:: ash_obj_mask
  integer(kind=int1), pointer, dimension(:,:), public, save:: ice_obj_mask
  integer(kind=int1), pointer, dimension(:,:), public, save:: hot_obj_mask
  real(kind=real4), allocatable, dimension(:), public, save:: pixel_area
  integer(kind=int4), public, save :: pixel_area_status

!--- nwp parameters
 integer, allocatable, dimension(:,:), public, save, target :: Zen_Idx_Rtm

 integer, allocatable, dimension(:,:), public, save, target :: I_Nwp
 integer, allocatable, dimension(:,:), public, save, target :: J_Nwp
 integer, allocatable, dimension(:,:), public, save, target :: I_Nwp_x
 integer, allocatable, dimension(:,:), public, save, target :: J_Nwp_x
 real, allocatable, dimension(:,:), public, save, target :: Lon_Nwp_Fac
 real, allocatable, dimension(:,:), public, save, target :: Lat_Nwp_Fac


!--- local radiative center
integer, allocatable, dimension(:,:), public, save, target :: i_LRC
integer, allocatable, dimension(:,:), public, save, target :: j_LRC

!--- rtm parameters
  real (kind=real4), dimension(:,:), allocatable, public, target:: Tsfc_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tair_Nwp_Pix !changed to target
  real (kind=real4), dimension(:,:), allocatable, public, save:: Rh_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Pmsl_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public, target:: Psfc_Nwp_Pix ! changed to target
  real (kind=real4), dimension(:,:), allocatable, public:: Tpw_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Ozone_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: K_Index_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Sc_Lwp_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Lwp_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Iwp_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Cwp_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Pc_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: LCL_Height_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: CCL_Height_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Cfrac_Nwp_Pix
  integer (kind=int1), dimension(:,:), allocatable, public:: Ncld_Layers_Nwp_Pix
  integer (kind=int1), dimension(:,:), allocatable, public:: Cld_Type_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Wnd_Spd_10m_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Wnd_Dir_10m_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Wnd_Spd_Cld_Top_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Wnd_Dir_Cld_Top_Nwp_Pix

  real (kind=real4), dimension(:,:), allocatable, public:: Trans_Atm_Ch20_Solar_Total_Rtm
  real (kind=real4), dimension(:,:), allocatable, public:: Trans_Atm_Ch20_Solar_Rtm
  real (kind=real4), dimension(:,:), allocatable, public:: Bt_Clear_Ch20_Solar_Rtm
  real (kind=real4), dimension(:,:), allocatable, public:: Rad_Clear_Ch20_Solar_Rtm
  real (kind=real4), dimension(:,:), allocatable, public:: Ems_Ch20_Clear_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ems_Ch20_Clear_Solar_Rtm
  real (kind=real4), dimension(:,:), allocatable, public:: Ems_Ch20_Clear_Solar_Sfc_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ttropo_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public, target:: Emiss_11um_Tropo_LRC
  real (kind=real4), dimension(:,:), allocatable, public, target:: Emiss_11um_Tropo_Nadir_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_12um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_85um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_67um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_133um_Tropo_Rtm


  real (kind=real4), dimension(:,:), allocatable, public:: Bt_Clear_Ch31_Rtm_unbiased
  real (kind=real4), dimension(:,:), allocatable, public:: Bt_Clear_Ch32_Rtm_unbiased
  real (kind=real4), dimension(:,:), allocatable, public:: Rad_Clear_Ch31_Rtm_unbiased
  real (kind=real4), dimension(:,:), allocatable, public:: Rad_Clear_Ch32_Rtm_unbiased


  real (kind=real4), dimension(:,:), allocatable, public, target:: Pc_Opaque_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Zc_Opaque_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Tc_Opaque_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Pc_Lower_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Zc_Lower_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Tc_Lower_Cloud

!--- modis white sky albedo maps
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ndvi_Sfc_White_Sky

!---- other static arrays carried by this module


  !--- Solar RTM Terms
  type, public :: solar_rtm_struct
      real, dimension(42,3):: Tau_H2O_Coef
      real, dimension(42):: Tau_Ray
      real, dimension(42):: Tau_O2
      real, dimension(42):: Tau_O3
      real, dimension(42):: Tau_CH4
      real, dimension(42):: Tau_CO2
      real, dimension(42):: Tau_Aer
      real, dimension(42):: Wo_Aer
      real, dimension(42):: G_Aer
  end type solar_rtm_struct
 
  type (solar_rtm_struct), public, save:: Solar_Rtm

!--- pixel arrays used in clavr_mod
  real(kind=real4), dimension(:,:), allocatable, public:: Pix_Data

!--- flags for using clavrxorb_Default_file
  integer ,public, save :: Use_Default

!--- clavrxorb_File_List filename
  character(len=300), public, save:: File_List

 contains

!----------------------------------------------------------------------------
! This routine allocate the memory for the pixel arrays 
!----------------------------------------------------------------------------
subroutine CREATE_PIXEL_ARRAYS()
  integer:: dim1 
  integer:: dim2 
  integer:: idx

  dim1 = Num_Pix
  dim2 = Num_Scans_Per_Segment

  do idx = 1,36
      if (Chan_On_Flag_Default(idx) == sym%YES) then
         allocate(Ch(idx)%Unc(dim1,dim2))
         if (idx <= 20 .or. idx == 26) then 
            allocate(Ch(idx)%Ref_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa_Unnorm(dim1,dim2))
            allocate(Ch(idx)%Ref_Sfc(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm_Total(dim1,dim2))
         endif
         if (idx >= 20 .and. idx /= 26) then 
            allocate(Ch(idx)%Rad_Toa(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa(dim1,dim2))
            allocate(Ch(idx)%Bt_Sfc(dim1,dim2))
            allocate(Ch(idx)%Rad_Sfc(dim1,dim2))
            allocate(Ch(idx)%Rad_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Rad_Atm(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm(dim1,dim2))
            allocate(Ch(idx)%Sfc_Emiss(dim1,dim2))
         endif
      endif
   enddo

   !--- DNB Variable
   idx = 42
   if (Chan_On_Flag_Default(idx) == sym%YES) then
      allocate(Ch(idx)%Rad_Toa(dim1,dim2))
      allocate(Ch(idx)%Ref_Toa(dim1,dim2))
      allocate(Ch(idx)%Ref_Lunar_Toa(dim1,dim2))
      allocate(Ch(idx)%Ref_Lunar_Toa_Clear(dim1,dim2))
      allocate(Ch(idx)%Ref_Lunar_Sfc(dim1,dim2))
      allocate(Ch(idx)%Trans_Atm_Total(dim1,dim2))
      allocate(Ch(idx)%Unc(dim1,dim2))
   endif

   allocate(Chan_On_Flag(Nchan_Clavrx,Num_Scans_Per_Segment))

   if ((Chan_On_Flag_Default(27) == sym%YES) .and.   &
       (Chan_On_Flag_Default(31) == sym%YES)) then
           allocate(Covar_Ch27_Ch31_5x5(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
   endif

   allocate(Bad_Pixel_Mask(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment), &
            Volcano_Mask(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment),&
            Space_Mask(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))

   allocate(Sfc_Level_Rtm_Pixel(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
   allocate(Solar_Contamination_Mask(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
    
   !--- VIIRS Arrays
   allocate(IFF_Gap_Mask(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
   allocate(Gap_Pixel_Mask(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
   allocate(Gap_Line_Idx(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))

   !--- 3x3 uni
   if (Chan_On_Flag_Default(27) == sym%YES) then
          allocate (Bt_Ch27_Mean_3x3(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
          allocate (Bt_Ch27_Max_3x3(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
          allocate (Bt_Ch27_Min_3x3(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
          allocate (Bt_Ch27_Std_3x3(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
   endif

   allocate(   &
          Dust(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment), &
          Smoke(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment), &
          Fire(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment), &
          Sst_Anal(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment), &
          Sst_Anal_Err(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment), &
          Sst_Anal_Cice(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment), &
          Sst_Anal_Uni(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))

  allocate(Beta_11um_12um_Tropo_Rtm(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
  allocate(Beta_11um_85um_Tropo_Rtm(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
  allocate(Beta_11um_67um_Tropo_Rtm(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
  allocate(Beta_11um_133um_Tropo_Rtm(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))

  allocate(Temp_Mask(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))


  allocate(Dcc_Mask(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))

  !--- geometry arrays
  call  CREATE_GEO_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  !--- anchor point arrays
  call  CREATE_GEO_ANCHOR_ARRAYS(Num_Anchors, Num_Scans_Per_Segment)
  !--- nwp fields interpolated to the pixel level
  call  CREATE_NWP_PIX_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_SURFACE_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  !--- ch1 arrays
  call  CREATE_REF_CHANNEL_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_THERM_CHANNEL_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_EXTRA_CHANNEL_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_LUNAR_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_BTD_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_ACHA_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_DCOMP_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_NLCOMP_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_SASRAB_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_OLR_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_AEROSOL_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_CLOUD_MASK_ARRAYS(Num_Pix, Num_Scans_Per_Segment, Max_Num_Cld_Test_Bytes)
  call  CREATE_CLOUD_TYPE_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_DIAGNOSTIC_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_SFC_PROD_ARRAYS(Num_Pix, Num_Scans_Per_Segment)
  call  CREATE_CLOUD_PROD_ARRAYS(Num_Pix, Num_Scans_Per_Segment)

  !--- pixel level parameters
   allocate(Zen_Idx_Rtm(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
   allocate(i_LRC(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))
   allocate(j_LRC(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))

   
  allocate(Ch3a_On_AVHRR(Num_Scans_Per_Segment), &
           Bad_Scan_Flag(Num_Scans_Per_Segment), &
           Scan_Number(Num_Scans_Per_Segment), &
           Scan_Time(Num_Scans_Per_Segment), &
           Scan_Day(Num_Scans_Per_Segment), &
           Scan_Year(Num_Scans_Per_Segment), &
           Utc_Scan_Time_Hours(Num_Scans_Per_Segment), &
           Pixel_Local_Time_Hours(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment), &
           Pixel_Time(Num_Pix,Line_Idx_Min_Segment:Line_Idx_Max_Segment))

  allocate(Pix_Data(nchan_Avhrr,Num_Pix))

   !--------------------------------------------------------------------------------
   ! Initialize variables that are not reset for each segment
   !--------------------------------------------------------------------------------

   !--- metrics - needed initialize counts to be zero but they accumulate through orbit
   ACHA_Processed_Count = 0
   ACHA_Valid_Count = 0
   ACHA_Success_Fraction = Missing_Value_Real4
   DCOMP_Processed_Count = 0
   DCOMP_Valid_Count = 0
   DCOMP_Success_Fraction = Missing_Value_Real4
   Nonconfident_Cloud_Mask_Fraction = Missing_Value_Real4
   Nonconfident_Cloud_Mask_Count = 0
   Cloud_Mask_Count = 0

   !--- other
   DCOMP_version_string = "not DCOMP source id  available"

end subroutine CREATE_PIXEL_ARRAYS

!------------------------------------------------------------------------------
subroutine DESTROY_PIXEL_ARRAYS()

  integer:: idx

  do idx = 1,42
      if (allocated(Ch(idx)%Rad_Toa)) deallocate(Ch(idx)%Rad_Toa)
      if (allocated(Ch(idx)%Bt_Toa)) deallocate(Ch(idx)%Bt_Toa)
      if (allocated(Ch(idx)%Rad_Toa_Clear)) deallocate(Ch(idx)%Rad_Toa_Clear)
      if (allocated(Ch(idx)%Bt_Toa_Clear)) deallocate(Ch(idx)%Bt_Toa_Clear)
      if (allocated(Ch(idx)%Rad_Sfc)) deallocate(Ch(idx)%Rad_Sfc)
      if (allocated(Ch(idx)%Bt_Sfc)) deallocate(Ch(idx)%Bt_Sfc)
      if (allocated(Ch(idx)%Rad_Atm)) deallocate(Ch(idx)%Rad_Atm)
      if (allocated(Ch(idx)%Trans_Atm)) deallocate(Ch(idx)%Trans_Atm)
      if (allocated(Ch(idx)%Trans_Atm_Total)) deallocate(Ch(idx)%Trans_Atm_Total)
      if (allocated(Ch(idx)%Ref_Toa)) deallocate(Ch(idx)%Ref_Toa)
      if (allocated(Ch(idx)%Ref_Toa_Unnorm)) deallocate(Ch(idx)%Ref_Toa_Unnorm)
      if (allocated(Ch(idx)%Ref_Toa_Clear)) deallocate(Ch(idx)%Ref_Toa_Clear)
      if (allocated(Ch(idx)%Ref_Sfc)) deallocate(Ch(idx)%Ref_Sfc)
      if (allocated(Ch(idx)%Ref_Lunar_Toa)) deallocate(Ch(idx)%Ref_Lunar_Toa)
      if (allocated(Ch(idx)%Ref_Lunar_Toa_Clear)) deallocate(Ch(idx)%Ref_Lunar_Toa_Clear)
      if (allocated(Ch(idx)%Ref_Lunar_Sfc)) deallocate(Ch(idx)%Ref_Lunar_Sfc)
      if (allocated(Ch(idx)%Sfc_Emiss)) deallocate(Ch(idx)%Sfc_Emiss)
      if (allocated(Ch(idx)%Sfc_Ref_White_Sky)) deallocate(Ch(idx)%Sfc_Ref_White_Sky)
      if (allocated(Ch(idx)%Emiss_Tropo)) deallocate(Ch(idx)%Emiss_Tropo)
      if (allocated(Ch(idx)%Unc)) deallocate(Ch(idx)%Unc)
  enddo

  if (allocated(Chan_On_Flag)) deallocate(Chan_On_Flag)

  if (allocated(Temp_Mask)) deallocate(Temp_Mask)

  deallocate(Ch3a_On_AVHRR, &
             Bad_Scan_Flag,  &
             Scan_Number, &
             Scan_Time,Scan_Day,Scan_Year,Utc_Scan_Time_Hours, &
             Pixel_Local_Time_Hours,Pixel_Time)

  if (allocated(Covar_Ch27_Ch31_5x5)) deallocate(Covar_Ch27_Ch31_5x5)

  deallocate(Pix_Data)


  deallocate(Bad_Pixel_Mask)

  call DESTROY_GEO_ARRAYS()
  call DESTROY_GEO_ANCHOR_ARRAYS()
  call DESTROY_NWP_PIX_ARRAYS()
  call DESTROY_SURFACE_ARRAYS()
  call DESTROY_REF_CHANNEL_ARRAYS()
  call DESTROY_THERM_CHANNEL_ARRAYS()
  call DESTROY_EXTRA_CHANNEL_ARRAYS()
  call DESTROY_LUNAR_ARRAYS()
  call DESTROY_BTD_ARRAYS()
  call DESTROY_ACHA_ARRAYS()
  call DESTROY_DCOMP_ARRAYS()
  call DESTROY_NLCOMP_ARRAYS()
  call DESTROY_SASRAB_ARRAYS()
  call DESTROY_OLR_ARRAYS()
  call DESTROY_AEROSOL_ARRAYS()
  call DESTROY_CLOUD_MASK_ARRAYS()
  call DESTROY_CLOUD_TYPE_ARRAYS()
  call DESTROY_DIAGNOSTIC_ARRAYS()
  call DESTROY_SFC_PROD_ARRAYS()
  call DESTROY_CLOUD_PROD_ARRAYS()

  deallocate(Volcano_Mask)
  deallocate(Space_Mask)
  deallocate(Sfc_Level_Rtm_Pixel)
  deallocate(Fire)
  deallocate(Dust)
  deallocate(Smoke)

  deallocate(Sst_Anal)
  deallocate(Sst_Anal_Err)
  deallocate(Sst_Anal_Cice)
  deallocate(Sst_Anal_Uni)

  deallocate(Solar_Contamination_Mask)
  deallocate(Dcc_Mask)

  !--- nwp and rtm indices
  if (allocated(Zen_Idx_Rtm)) deallocate(Zen_Idx_Rtm)

  !--- local radiative center indices
  if (allocated(i_LRC)) deallocate(i_LRC)
  if (allocated(j_LRC)) deallocate(j_LRC)

!--- pixel rtm parameters
  if (Chan_On_Flag_Default(27) == sym%YES) then
    if (allocated(Bt_Ch27_Mean_3x3)) deallocate(Bt_Ch27_Mean_3x3)
    if (allocated(Bt_Ch27_Max_3x3)) deallocate(Bt_Ch27_Max_3x3)
    if (allocated(Bt_Ch27_Min_3x3)) deallocate(Bt_Ch27_Min_3x3)
    if (allocated(Bt_Ch27_Std_3x3)) deallocate(Bt_Ch27_Std_3x3)
  endif
  if (Chan_On_Flag_Default(32) == sym%YES) then
    if (allocated(Rad_Clear_Ch32_Rtm_Unbiased)) deallocate(Rad_Clear_Ch32_Rtm_Unbiased)
    if (allocated(Bt_Clear_Ch32_Rtm_Unbiased)) deallocate(Bt_Clear_Ch32_Rtm_Unbiased)
  endif

  !--- ir cloud layer
  if (allocated(Beta_11um_12um_Tropo_Rtm)) deallocate(Beta_11um_12um_Tropo_Rtm)
  if (allocated(Beta_11um_67um_Tropo_Rtm)) deallocate(Beta_11um_67um_Tropo_Rtm)
  if (allocated(Beta_11um_85um_Tropo_Rtm)) deallocate(Beta_11um_85um_Tropo_Rtm)
  if (allocated(Beta_11um_133um_Tropo_Rtm)) deallocate(Beta_11um_133um_Tropo_Rtm)

  !--- VIIRS
  if (allocated(Gap_Pixel_Mask)) deallocate(Gap_Pixel_Mask)
  if (allocated(Gap_Line_Idx)) deallocate(Gap_Line_Idx)
  if (allocated(Gap_Pixel_Mask_Pattern)) deallocate(Gap_Pixel_Mask_Pattern)
  if (allocated(Gap_Line_Idx_Pattern)) deallocate(Gap_Line_Idx_Pattern)
  if (allocated(IFF_Gap_Mask)) deallocate(IFF_Gap_Mask)

end subroutine DESTROY_PIXEL_ARRAYS

!-----------------------------------------------------------------------------
!  reset pixel arrays to missing
!-----------------------------------------------------------------------------
subroutine RESET_PIXEL_ARRAYS_TO_MISSING()

      Bad_Scan_Flag = sym%NO        !not initialized to missing
      Bad_Pixel_Mask = sym%YES      !not initialized to missing
      Ch3a_On_AVHRR = Missing_Value_Int1

      call RESET_GEO_ARRAYS()
      call RESET_GEO_ANCHOR_ARRAYS()
      call RESET_NWP_PIX_ARRAYS()
      call RESET_SURFACE_ARRAYS()
      call RESET_REF_CHANNEL_ARRAYS()
      call RESET_THERM_CHANNEL_ARRAYS()
      call RESET_EXTRA_CHANNEL_ARRAYS()
      call RESET_LUNAR_ARRAYS()
      call RESET_BTD_ARRAYS()
      call RESET_ACHA_ARRAYS()
      call RESET_DCOMP_ARRAYS()
      call RESET_NLCOMP_ARRAYS()
      call RESET_SASRAB_ARRAYS()
      call RESET_OLR_ARRAYS()
      call RESET_AEROSOL_ARRAYS()
      call RESET_CLOUD_MASK_ARRAYS()
      call RESET_CLOUD_TYPE_ARRAYS()
      call RESET_DIAGNOSTIC_ARRAYS()
      call RESET_SFC_PROD_ARRAYS()
      call RESET_CLOUD_PROD_ARRAYS()

      if (Chan_On_Flag_Default(27) == sym%YES) THEN
        Bt_Ch27_Mean_3x3 = Missing_Value_Real4
        Bt_Ch27_Min_3x3 = Missing_Value_Real4
        Bt_Ch27_Max_3x3 = Missing_Value_Real4
        Bt_Ch27_Std_3x3 = Missing_Value_Real4
      endif

      Dcc_Mask = Missing_Value_Int1
      Sst_Anal = Missing_Value_Real4

      Zen_Idx_Rtm = Missing_Value_Int1

      if ((Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Chan_On_Flag_Default(32) == sym%YES)) then
        Beta_11um_12um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Chan_On_Flag_Default(27) == sym%YES)) then
        Beta_11um_67um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Chan_On_Flag_Default(29) == sym%YES)) then
        Beta_11um_85um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Chan_On_Flag_Default(33) == sym%YES)) then
        Beta_11um_133um_Tropo_Rtm = Missing_Value_Real4
      endif

      Space_Mask = Missing_Value_Int1
      Sfc_Level_Rtm_Pixel = Missing_Value_Int4
    
      Scan_Time = Missing_Value_Int4
      Utc_Scan_Time_Hours = Missing_Value_Real4
      Pixel_Local_Time_Hours = Missing_Value_Real4
      Pixel_Time = Missing_Value_Real4

      i_LRC = Missing_Value_Int4
      j_LRC = Missing_Value_Int4

      !--- note do not reset the patterns
      Gap_Pixel_Mask = sym%NO
      Gap_Line_Idx = Missing_Value_Int4
      IFF_Gap_Mask = sym%NO     

end subroutine RESET_PIXEL_ARRAYS_TO_MISSING
!-------------------------------------------------------------
! The following routines provide an internal structure to 
! to this module
!
! What follows are a pattern of 3 routines per data group
! create - allocates arrays
! reset - sets arrays to missing
! destroy - deallocate arrays
!-------------------------------------------------------------

!-------------------------------------------------------------
!
!-------------------------------------------------------------
subroutine CREATE_GEO_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(Ascend(dim2))
   allocate(Lat(dim1,dim2))
   allocate(Lon(dim1,dim2))
   allocate(Lat_1b(dim1,dim2))
   allocate(Lon_1b(dim1,dim2))
   allocate(Lat_Pc(dim1,dim2))
   allocate(Lon_Pc(dim1,dim2))
   allocate(Satzen(dim1,dim2))
   allocate(Solzen(dim1,dim2))
   allocate(Sataz(dim1,dim2))
   allocate(Solaz(dim1,dim2))
   allocate(Relaz(dim1,dim2))
   allocate(Coszen(dim1,dim2))
   allocate(Seczen(dim1,dim2))
   allocate(CosSolzen(dim1,dim2))
   allocate(Glintzen(dim1,dim2))
   allocate(Airmass(dim1,dim2))
   allocate(Scatangle(dim1,dim2))
end subroutine CREATE_GEO_ARRAYS
subroutine DESTROY_GEO_ARRAYS()
  deallocate(Ascend)
  deallocate(Lat)
  deallocate(Lon)
  deallocate(Lat_1b)
  deallocate(Lon_1b)
  deallocate(Lat_Pc)
  deallocate(Lon_Pc)
  deallocate(Satzen)
  deallocate(Solzen)
  deallocate(Sataz)
  deallocate(Solaz)
  deallocate(Relaz)
  deallocate(Coszen)
  deallocate(Seczen)
  deallocate(CosSolzen)
  deallocate(Glintzen)
  deallocate(Airmass)
  deallocate(Scatangle)
end subroutine DESTROY_GEO_ARRAYS
subroutine RESET_GEO_ARRAYS()
  Ascend = Missing_Value_Int1
  Lat = Missing_Value_Real4
  Lon = Missing_Value_Real4
  Lat_1b = Missing_Value_Real4
  Lon_1b = Missing_Value_Real4
  Lat_Pc = Missing_Value_Real4
  Lon_Pc = Missing_Value_Real4
  Satzen = Missing_Value_Real4
  Solzen = Missing_Value_Real4
  Coszen = Missing_Value_Real4
  Seczen = Missing_Value_Real4
  CosSolzen = Missing_Value_Real4
  Relaz = Missing_Value_Real4
  Solaz = Missing_Value_Real4
  Sataz = Missing_Value_Real4
  Glintzen = Missing_Value_Real4
  Scatangle = Missing_Value_Real4
  Airmass = Missing_Value_Real4
end subroutine RESET_GEO_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_GEO_ANCHOR_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(Lat_Anchor_1b(dim1,dim2))
   allocate(Lon_Anchor_1b(dim1,dim2))
   allocate(Solzen_Anchor(dim1,dim2))
   allocate(Satzen_Anchor(dim1,dim2))
   allocate(Scatangle_Anchor(dim1,dim2))
   allocate(Glintzen_Anchor(dim1,dim2))
   allocate(Relaz_Anchor(dim1,dim2))
   allocate(Solaz_Anchor(dim1,dim2))
   allocate(Sataz_Anchor(dim1,dim2))
end subroutine CREATE_GEO_ANCHOR_ARRAYS
subroutine RESET_GEO_ANCHOR_ARRAYS()
  Lat_Anchor_1b = Missing_Value_Real4
  Lon_Anchor_1b = Missing_Value_Real4
  Solzen_Anchor = Missing_Value_Real4
  Satzen_Anchor = Missing_Value_Real4
  Scatangle_Anchor = Missing_Value_Real4
  Glintzen_Anchor = Missing_Value_Real4
  Relaz_Anchor = Missing_Value_Real4
  Solaz_Anchor = Missing_Value_Real4
  Sataz_Anchor = Missing_Value_Real4
end subroutine RESET_GEO_ANCHOR_ARRAYS
subroutine DESTROY_GEO_ANCHOR_ARRAYS
  deallocate(Lat_Anchor_1b)
  deallocate(Lon_Anchor_1b)
  deallocate(Solzen_Anchor)
  deallocate(Satzen_Anchor)
  deallocate(Scatangle_Anchor)
  deallocate(Glintzen_Anchor)
  deallocate(Relaz_Anchor)
  deallocate(Solaz_Anchor)
  deallocate(Sataz_Anchor)
end subroutine DESTROY_GEO_ANCHOR_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_NWP_PIX_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(Tair_Nwp_Pix(dim1,dim2))
   allocate(Rh_Nwp_Pix(dim1,dim2))
   allocate(LCL_Height_Nwp_Pix(dim1,dim2))
   allocate(CCL_Height_Nwp_Pix(dim1,dim2))
   allocate(Tsfc_Nwp_Pix(dim1,dim2))
   allocate(Pmsl_Nwp_Pix(dim1,dim2))
   allocate(Psfc_Nwp_Pix(dim1,dim2))
   allocate(K_Index_Nwp_Pix(dim1,dim2))
   allocate(Sc_Lwp_Nwp_Pix(dim1,dim2))
   allocate(Lwp_Nwp_Pix(dim1,dim2))
   allocate(Iwp_Nwp_Pix(dim1,dim2))
   allocate(Cwp_Nwp_Pix(dim1,dim2))
   allocate(Pc_Nwp_Pix(dim1,dim2))
   allocate(Cfrac_Nwp_Pix(dim1,dim2))
   allocate(Ncld_Layers_Nwp_Pix(dim1,dim2))
   allocate(Cld_Type_Nwp_Pix(dim1,dim2))
   allocate(Tpw_Nwp_Pix(dim1,dim2))
   allocate(Ozone_Nwp_Pix(dim1,dim2))
   allocate(Ttropo_Nwp_Pix(dim1,dim2))
   allocate(Wnd_Spd_10m_Nwp_Pix(dim1,dim2))
   allocate(Wnd_Dir_10m_Nwp_Pix(dim1,dim2))
   allocate(Wnd_Spd_Cld_Top_Nwp_Pix(dim1,dim2))
   allocate(Wnd_Dir_Cld_Top_Nwp_Pix(dim1,dim2))
   allocate(I_Nwp(dim1,dim2))
   allocate(J_Nwp(dim1,dim2))
   allocate(I_Nwp_x(dim1,dim2))
   allocate(J_Nwp_x(dim1,dim2))
   allocate(Lon_Nwp_Fac(dim1,dim2))
   allocate(Lat_Nwp_Fac(dim1,dim2))
end subroutine CREATE_NWP_PIX_ARRAYS
subroutine RESET_NWP_PIX_ARRAYS()
   Tair_Nwp_Pix = Missing_Value_Real4
   Rh_Nwp_Pix = Missing_Value_Real4
   LCL_Height_Nwp_Pix = Missing_Value_Real4
   CCL_Height_Nwp_Pix = Missing_Value_Real4
   Tsfc_Nwp_Pix = Missing_Value_Real4
   Pmsl_Nwp_Pix = Missing_Value_Real4
   Psfc_Nwp_Pix = Missing_Value_Real4
   K_Index_Nwp_Pix = Missing_Value_Real4
   Sc_Lwp_Nwp_Pix = Missing_Value_Real4
   Lwp_Nwp_Pix = Missing_Value_Real4
   Iwp_Nwp_Pix = Missing_Value_Real4
   Cwp_Nwp_Pix = Missing_Value_Real4
   Pc_Nwp_Pix = Missing_Value_Real4
   Cfrac_Nwp_Pix = Missing_Value_Real4
   Ncld_Layers_Nwp_Pix = Missing_Value_Int1
   Cld_Type_Nwp_Pix = Missing_Value_Int1
   Tpw_Nwp_Pix = Missing_Value_Real4
   Ozone_Nwp_Pix = Missing_Value_Real4
   Ttropo_Nwp_Pix = Missing_Value_Real4
   Wnd_Spd_10m_Nwp_Pix = Missing_Value_Real4
   Wnd_Dir_10m_Nwp_Pix = Missing_Value_Real4
   Wnd_Spd_Cld_Top_Nwp_Pix = Missing_Value_Real4
   Wnd_Dir_Cld_Top_Nwp_Pix = Missing_Value_Real4
   I_Nwp = Missing_Value_Int4
   J_Nwp = Missing_Value_Int4
   I_Nwp_x = Missing_Value_Int4
   J_Nwp_x = Missing_Value_Int4
   Lon_Nwp_Fac = Missing_Value_Real4
   Lat_Nwp_Fac = Missing_Value_Real4
end subroutine RESET_NWP_PIX_ARRAYS
subroutine DESTROY_NWP_PIX_ARRAYS()
   deallocate(Tair_Nwp_Pix)
   deallocate(Rh_Nwp_Pix)
   deallocate(LCL_Height_Nwp_Pix)
   deallocate(CCL_Height_Nwp_Pix)
   deallocate(Tsfc_Nwp_Pix)
   deallocate(Pmsl_Nwp_Pix)
   deallocate(Psfc_Nwp_Pix)
   deallocate(K_Index_Nwp_Pix)
   deallocate(Sc_Lwp_Nwp_Pix)
   deallocate(Lwp_Nwp_Pix)
   deallocate(Iwp_Nwp_Pix)
   deallocate(Cwp_Nwp_Pix)
   deallocate(Pc_Nwp_Pix)
   deallocate(Cfrac_Nwp_Pix)
   deallocate(Ncld_Layers_Nwp_Pix)
   deallocate(Cld_Type_Nwp_Pix)
   deallocate(Tpw_Nwp_Pix)
   deallocate(Ozone_Nwp_Pix)
   deallocate(Ttropo_Nwp_Pix)
   deallocate(Wnd_Spd_10m_Nwp_Pix)
   deallocate(Wnd_Dir_10m_Nwp_Pix)
   deallocate(Wnd_Spd_Cld_Top_Nwp_Pix)
   deallocate(Wnd_Dir_Cld_Top_Nwp_Pix)
   deallocate(I_Nwp)
   deallocate(J_Nwp)
   deallocate(I_Nwp_x)
   deallocate(J_Nwp_x)
   deallocate(Lon_Nwp_Fac)
   deallocate(Lat_Nwp_Fac)
end subroutine DESTROY_NWP_PIX_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_REF_CHANNEL_ARRAYS(dim1,dim2)

   integer, intent(in):: dim1, dim2

   if (Chan_On_Flag_Default(1) == sym%YES) then
           allocate(Ch1_Counts(dim1,dim2))
           allocate(Ref_Ch1_Mean_3x3(dim1,dim2))
           allocate(Ref_Ch1_Max_3x3(dim1,dim2))
           allocate(Ref_Ch1_Min_3x3(dim1,dim2))
           allocate(Ref_Ch1_Std_3x3(dim1,dim2))
           allocate(Ref_Ch1_Clear_Mean_3x3(dim1,dim2))
           allocate(Ref_Ch1_Clear_Max_3x3(dim1,dim2))
           allocate(Ref_Ch1_Clear_Min_3x3(dim1,dim2))
           allocate(Ref_Ch1_Clear_Std_3x3(dim1,dim2))
           allocate(Ref_Ch1_Dark_Composite(dim1,dim2))
           allocate(Ref_Ch1_Sfc_White_Sky_Mean_3x3(dim1,dim2))
           allocate(Ref_Ch1_Sfc_White_Sky_Max_3x3(dim1,dim2))
           allocate(Ref_Ch1_Sfc_White_Sky_Min_3x3(dim1,dim2))
           allocate(Ref_Ch1_Sfc_White_Sky_Std_3x3(dim1,dim2))
   endif

   if (Chan_On_Flag_Default(2) == sym%YES) then
      allocate(Ch2_Counts(dim1,dim2))
   endif

   if (Chan_On_Flag_Default(6) == sym%YES) then
      allocate(Ch6_Counts(dim1,dim2))
   endif
   allocate(Ch6_On_Pixel_Mask(dim1,dim2))

end subroutine CREATE_REF_CHANNEL_ARRAYS

subroutine RESET_REF_CHANNEL_ARRAYS
   integer:: idx

   do idx = 1,42
      if (allocated(Ch(idx)%Rad_Toa)) Ch(idx)%Rad_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa)) Ch(idx)%Bt_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Rad_Toa_Clear)) Ch(idx)%Rad_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa_Clear)) Ch(idx)%Bt_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Sfc)) Ch(idx)%Bt_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Rad_Sfc)) Ch(idx)%Rad_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Rad_Atm)) Ch(idx)%Rad_Atm = Missing_Value_Real4
      if (allocated(Ch(idx)%Trans_Atm)) Ch(idx)%Trans_Atm = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa)) Ch(idx)%Ref_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa_Unnorm)) Ch(idx)%Ref_Toa_Unnorm = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa_Clear)) Ch(idx)%Ref_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Sfc)) Ch(idx)%Ref_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Toa)) Ch(idx)%Ref_Lunar_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Toa_Clear)) Ch(idx)%Ref_Lunar_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Sfc)) Ch(idx)%Ref_Lunar_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Sfc_Emiss)) Ch(idx)%Sfc_Emiss = Missing_Value_Real4
      if (allocated(Ch(idx)%Sfc_Ref_White_Sky)) Ch(idx)%Sfc_Ref_White_Sky = Missing_Value_Real4
      if (allocated(Ch(idx)%Emiss_Tropo)) Ch(idx)%Emiss_Tropo = Missing_Value_Real4
      if (allocated(Ch(idx)%Unc)) Ch(idx)%Unc = Missing_Value_Int1
   enddo

   if (Chan_On_Flag_Default(1) == sym%YES) then
      Ch1_Counts = Missing_Value_Int2
      Ref_Ch1_Mean_3x3 = Missing_Value_Real4
      Ref_Ch1_Min_3x3 = Missing_Value_Real4
      Ref_Ch1_Max_3x3 = Missing_Value_Real4
      Ref_Ch1_Std_3x3 = Missing_Value_Real4
      Ref_Ch1_Clear_Mean_3x3 = Missing_Value_Real4
      Ref_Ch1_Clear_Max_3x3 = Missing_Value_Real4
      Ref_Ch1_Clear_Min_3x3 = Missing_Value_Real4
      Ref_Ch1_Clear_Std_3x3 = Missing_Value_Real4
      Ref_Ch1_Dark_Composite = Missing_Value_Real4
      Ref_Ch1_Sfc_White_Sky_Mean_3x3 = Missing_Value_Real4
      Ref_Ch1_Sfc_White_Sky_Max_3x3 = Missing_Value_Real4
      Ref_Ch1_Sfc_White_Sky_Min_3x3 = Missing_Value_Real4
      Ref_Ch1_Sfc_White_Sky_Std_3x3 = Missing_Value_Real4
   endif

   if (Chan_On_Flag_Default(2) == sym%YES) then
      Ch2_Counts = Missing_Value_Int2
   endif

   if (Chan_On_Flag_Default(6) == sym%YES) then
      Ch6_Counts = Missing_Value_Int2
   endif
   Ch6_On_Pixel_Mask = Missing_Value_Int1

end subroutine RESET_REF_CHANNEL_ARRAYS
subroutine DESTROY_REF_CHANNEL_ARRAYS

  if (Chan_On_Flag_Default(1) == sym%YES) then
   if (allocated(Ch1_Counts)) deallocate (Ch1_Counts)
   if (allocated(Ref_Ch1_Mean_3x3)) deallocate (Ref_Ch1_Mean_3x3)
   if (allocated(Ref_Ch1_Max_3x3)) deallocate (Ref_Ch1_Max_3x3)
   if (allocated(Ref_Ch1_Min_3x3)) deallocate (Ref_Ch1_Min_3x3)
   if (allocated(Ref_Ch1_Std_3x3)) deallocate (Ref_Ch1_Std_3x3)
   if (allocated(Ref_Ch1_Clear_Mean_3x3)) deallocate (Ref_Ch1_Clear_Mean_3x3)
   if (allocated(Ref_Ch1_Clear_Max_3x3)) deallocate (Ref_Ch1_Clear_Max_3x3)
   if (allocated(Ref_Ch1_Clear_Min_3x3)) deallocate (Ref_Ch1_Clear_Min_3x3)
   if (allocated(Ref_Ch1_Clear_Std_3x3)) deallocate (Ref_Ch1_Clear_Std_3x3)
   if (allocated(Ref_Ch1_Dark_Composite)) deallocate (Ref_Ch1_Dark_Composite)
   if (allocated(Ref_Ch1_Sfc_White_Sky_Mean_3x3)) deallocate (Ref_Ch1_Sfc_White_Sky_Mean_3x3)
   if (allocated(Ref_Ch1_Sfc_White_Sky_Max_3x3)) deallocate (Ref_Ch1_Sfc_White_Sky_Max_3x3)
   if (allocated(Ref_Ch1_Sfc_White_Sky_Min_3x3)) deallocate (Ref_Ch1_Sfc_White_Sky_Min_3x3)
   if (allocated(Ref_Ch1_Sfc_White_Sky_Std_3x3)) deallocate (Ref_Ch1_Sfc_White_Sky_Std_3x3)
  endif

   if (Chan_On_Flag_Default(2) == sym%YES) then
      deallocate(Ch2_Counts)
   endif

   if (Chan_On_Flag_Default(6) == sym%YES) then
      deallocate(Ch6_Counts)
   endif
   deallocate(Ch6_On_Pixel_Mask)

end subroutine DESTROY_REF_CHANNEL_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_THERM_CHANNEL_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2

   if (Chan_On_Flag_Default(20) == sym%YES) then
       allocate(Rad_Ch20_ems(dim1,dim2))
       allocate(Ems_Ch20(dim1,dim2))
       allocate(Bt_Ch20_Mean_3x3(dim1,dim2))
       allocate(Bt_Ch20_Max_3x3(dim1,dim2))
       allocate(Bt_Ch20_Min_3x3(dim1,dim2))
       allocate(Bt_Ch20_Std_3x3(dim1,dim2))
       allocate(Ems_Ch20_Mean_3x3(dim1,dim2))
       allocate(Ems_Ch20_Max_3x3(dim1,dim2))
       allocate(Ems_Ch20_Min_3x3(dim1,dim2))
       allocate(Ems_Ch20_Std_3x3(dim1,dim2))
       allocate(Bt_Ch20_Median_3x3(dim1,dim2))
       allocate(Ems_Ch20_Median_3x3(dim1,dim2))
       allocate(Ems_Ch20_Std_Median_3x3(dim1,dim2))
       allocate(Bt_Ch20_Std_Median_3x3(dim1,dim2))
       allocate(Ch20_Counts_Filtered(dim1,dim2))
       allocate(Rad_Clear_Ch20_Solar_Rtm(dim1,dim2))
       allocate(Bt_Clear_Ch20_Solar_Rtm(dim1,dim2))
       allocate(Trans_Atm_Ch20_Solar_Rtm(dim1,dim2))
       allocate(Trans_Atm_Ch20_Solar_total_Rtm(dim1,dim2))
       allocate(Ems_Ch20_Clear_Rtm(dim1,dim2))
       allocate(Ems_Ch20_Clear_Solar_Rtm(dim1,dim2))
       allocate(Ems_Ch20_Clear_Solar_Sfc_Rtm(dim1,dim2))
       allocate(Rad_Ch20_ems_Sfc(dim1,dim2))
       allocate(Ems_Ch20_Sfc(dim1,dim2))
       allocate(Ref_Ch20_LRC(dim1,dim2))
   endif

   if (Chan_On_Flag_Default(31) == sym%YES) then
      allocate(Bt_Ch31_Mean_3x3(dim1,dim2))
      allocate(Bt_Ch31_Max_3x3(dim1,dim2))
      allocate(Bt_Ch31_Min_3x3(dim1,dim2))
      allocate(Bt_Ch31_Std_3x3(dim1,dim2))
      allocate(Elem_Idx_Max_Bt_Ch31_3x3(dim1,dim2))
      allocate(Line_Idx_Max_Bt_Ch31_3x3(dim1,dim2))
      allocate(Elem_Idx_Min_Bt_Ch31_3x3(dim1,dim2))
      allocate(Line_Idx_Min_Bt_Ch31_3x3(dim1,dim2))
      allocate(Emiss_11um_Tropo_Nadir_Rtm(dim1,dim2))
      allocate(Emiss_11um_Tropo_LRC(dim1,dim2))
      allocate(Rad_Clear_Ch31_Rtm_unbiased(dim1,dim2))
      allocate(Bt_Clear_Ch31_Rtm_unbiased(dim1,dim2))
      allocate(Bt_Ch31_LRC(dim1,dim2))
      allocate(Bt_Ch31_Max_LRC(dim1,dim2))
      allocate(Bt_Ch31_Std_LRC(dim1,dim2))
   endif

   if (Chan_On_Flag_Default(32) == sym%YES) then
     allocate(Rad_Clear_Ch32_Rtm_unbiased(dim1,dim2))
     allocate(Bt_Clear_Ch32_Rtm_unbiased(dim1,dim2))
   endif

end subroutine CREATE_THERM_CHANNEL_ARRAYS

subroutine RESET_THERM_CHANNEL_ARRAYS()

   if (Chan_On_Flag_Default(20) == sym%YES) then
       Rad_Ch20_ems = Missing_Value_Real4
       Ems_Ch20 = Missing_Value_Real4
       Bt_Ch20_Mean_3x3 = Missing_Value_Real4
       Bt_Ch20_Max_3x3 = Missing_Value_Real4
       Bt_Ch20_Min_3x3 = Missing_Value_Real4
       Bt_Ch20_Std_3x3 = Missing_Value_Real4
       Ems_Ch20_Mean_3x3 = Missing_Value_Real4
       Ems_Ch20_Max_3x3 = Missing_Value_Real4
       Ems_Ch20_Min_3x3 = Missing_Value_Real4
       Ems_Ch20_Std_3x3 = Missing_Value_Real4
       Bt_Ch20_Median_3x3 = Missing_Value_Real4
       Ems_Ch20_Median_3x3 = Missing_Value_Real4
       Ems_Ch20_Std_Median_3x3 = Missing_Value_Real4
       Bt_Ch20_Std_Median_3x3 = Missing_Value_Real4
       Ch20_Counts_Filtered = Missing_Value_Real4
       Rad_Clear_Ch20_Solar_Rtm = Missing_Value_Real4
       Bt_Clear_Ch20_Solar_Rtm = Missing_Value_Real4
       Trans_Atm_Ch20_Solar_Rtm = Missing_Value_Real4
       Trans_Atm_Ch20_Solar_total_Rtm = Missing_Value_Real4
       Ems_Ch20_Clear_Rtm = Missing_Value_Real4
       Ems_Ch20_Clear_Solar_Rtm = Missing_Value_Real4
       Ems_Ch20_Clear_Solar_Sfc_Rtm = Missing_Value_Real4
       Rad_Ch20_ems_Sfc = Missing_Value_Real4
       Ems_Ch20_Sfc = Missing_Value_Real4
       Ref_Ch20_LRC = Missing_Value_Real4
   endif

   if (Chan_On_Flag_Default(31) == sym%YES) then
      Bt_Ch31_Mean_3x3 = Missing_Value_Real4
      Bt_Ch31_Max_3x3 = Missing_Value_Real4
      Bt_Ch31_Min_3x3 = Missing_Value_Real4
      Bt_Ch31_Std_3x3 = Missing_Value_Real4
      Elem_Idx_Max_Bt_Ch31_3x3 = Missing_Value_Real4
      Line_Idx_Max_Bt_Ch31_3x3 = Missing_Value_Real4
      Elem_Idx_Min_Bt_Ch31_3x3 = Missing_Value_Real4
      Line_Idx_Min_Bt_Ch31_3x3 = Missing_Value_Real4
      Emiss_11um_Tropo_Nadir_Rtm = Missing_Value_Real4
      Emiss_11um_Tropo_LRC = Missing_Value_Real4
      Rad_Clear_Ch31_Rtm_unbiased = Missing_Value_Real4
      Bt_Clear_Ch31_Rtm_unbiased = Missing_Value_Real4
      Bt_Ch31_LRC = Missing_Value_Real4
      Bt_Ch31_Max_LRC = Missing_Value_Real4
      Bt_Ch31_Std_LRC = Missing_Value_Real4
   endif

   if (Chan_On_Flag_Default(32) == sym%YES) then
     Rad_Clear_Ch32_Rtm_unbiased = Missing_Value_Real4
     Bt_Clear_Ch32_Rtm_unbiased = Missing_Value_Real4
   endif

end subroutine RESET_THERM_CHANNEL_ARRAYS
subroutine DESTROY_THERM_CHANNEL_ARRAYS()

   if (Chan_On_Flag_Default(20) == sym%YES) then
       deallocate(Rad_Ch20_ems)
       deallocate(Ems_Ch20)
       deallocate(Bt_Ch20_Mean_3x3)
       deallocate(Bt_Ch20_Max_3x3)
       deallocate(Bt_Ch20_Min_3x3)
       deallocate(Bt_Ch20_Std_3x3)
       deallocate(Ems_Ch20_Mean_3x3)
       deallocate(Ems_Ch20_Max_3x3)
       deallocate(Ems_Ch20_Min_3x3)
       deallocate(Ems_Ch20_Std_3x3)
       deallocate(Bt_Ch20_Median_3x3)
       deallocate(Ems_Ch20_Median_3x3)
       deallocate(Ems_Ch20_Std_Median_3x3)
       deallocate(Bt_Ch20_Std_Median_3x3)
       deallocate(Ch20_Counts_Filtered)
       deallocate(Rad_Clear_Ch20_Solar_Rtm)
       deallocate(Bt_Clear_Ch20_Solar_Rtm)
       deallocate(Trans_Atm_Ch20_Solar_Rtm)
       deallocate(Trans_Atm_Ch20_Solar_total_Rtm)
       deallocate(Ems_Ch20_Clear_Rtm)
       deallocate(Ems_Ch20_Clear_Solar_Rtm)
       deallocate(Ems_Ch20_Clear_Solar_Sfc_Rtm)
       deallocate(Rad_Ch20_ems_Sfc)
       deallocate(Ems_Ch20_Sfc)
       deallocate(Ref_Ch20_LRC)
   endif
   if (Chan_On_Flag_Default(31) == sym%YES) then
      deallocate(Bt_Ch31_Mean_3x3)
      deallocate(Bt_Ch31_Max_3x3)
      deallocate(Bt_Ch31_Min_3x3)
      deallocate(Bt_Ch31_Std_3x3)
      deallocate(Elem_Idx_Max_Bt_Ch31_3x3)
      deallocate(Line_Idx_Max_Bt_Ch31_3x3)
      deallocate(Elem_Idx_Min_Bt_Ch31_3x3)
      deallocate(Line_Idx_Min_Bt_Ch31_3x3)
      deallocate(Emiss_11um_Tropo_Nadir_Rtm)
      deallocate(Emiss_11um_Tropo_LRC)
      deallocate(Rad_Clear_Ch31_Rtm_unbiased)
      deallocate(Bt_Clear_Ch31_Rtm_unbiased)
      deallocate(Bt_Ch31_LRC)
      deallocate(Bt_Ch31_Max_LRC)
      deallocate(Bt_Ch31_Std_LRC)
   endif
   if (Chan_On_Flag_Default(32) == sym%YES) then
     deallocate(Rad_Clear_Ch32_Rtm_unbiased)
     deallocate(Bt_Clear_Ch32_Rtm_unbiased)
  endif

end subroutine DESTROY_THERM_CHANNEL_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_EXTRA_CHANNEL_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Chan_On_Flag_Default(37) == sym%YES) then
           allocate(Ref_ChI1(2*dim1,2*dim2))
           allocate(Ref_Max_ChI1(dim1,dim2))
           allocate(Ref_Min_ChI1(dim1,dim2))
           allocate(Ref_Uni_ChI1(dim1,dim2))
           allocate(Ref_Mean_ChI1(dim1,dim2))
   endif
   if (Chan_On_Flag_Default(38) == sym%YES) then
           allocate(Ref_ChI2(2*dim1,2*dim2))
           allocate(Ref_Max_ChI2(dim1,dim2))
           allocate(Ref_Min_ChI2(dim1,dim2))
           allocate(Ref_Mean_ChI2(dim1,dim2))
   endif
   if (Chan_On_Flag_Default(39) == sym%YES) then
           allocate(Ref_ChI3(2*dim1,2*dim2))
   endif
   if (Chan_On_Flag_Default(40) == sym%YES) then
           allocate(Bt_ChI4(2*dim1,2*dim2))
   endif
   if (Chan_On_Flag_Default(41) == sym%YES) then
           allocate(Bt_ChI5(2*dim1,2*dim2))
           allocate(Bt_Max_ChI5(dim1,dim2))
           allocate(Bt_Min_ChI5(dim1,dim2))
           allocate(Bt_Mean_ChI5(dim1,dim2))
   endif
   if (Chan_On_Flag_Default(42) == sym%YES) then
           allocate(Ref_ChDNB_Lunar_Mean_3x3(dim1,dim2))
           allocate(Ref_ChDNB_Lunar_Max_3x3(dim1,dim2))
           allocate(Ref_ChDNB_Lunar_Min_3x3(dim1,dim2))
           allocate(Ref_ChDNB_Lunar_Std_3x3(dim1,dim2))
   endif
end subroutine CREATE_EXTRA_CHANNEL_ARRAYS

subroutine RESET_EXTRA_CHANNEL_ARRAYS()
      if (Chan_On_Flag_Default(37) == sym%YES) Ref_ChI1 = Missing_Value_Real4
      if (Chan_On_Flag_Default(37) == sym%YES) Ref_Max_ChI1 = Missing_Value_Real4
      if (Chan_On_Flag_Default(37) == sym%YES) Ref_Min_ChI1 = Missing_Value_Real4
      if (Chan_On_Flag_Default(37) == sym%YES) Ref_Uni_ChI1 = Missing_Value_Real4
      if (Chan_On_Flag_Default(37) == sym%YES) Ref_Mean_ChI1 = Missing_Value_Real4
      if (Chan_On_Flag_Default(38) == sym%YES) Ref_ChI2 = Missing_Value_Real4
      if (Chan_On_Flag_Default(38) == sym%YES) Ref_Max_ChI2 = Missing_Value_Real4
      if (Chan_On_Flag_Default(38) == sym%YES) Ref_Min_ChI2 = Missing_Value_Real4
      if (Chan_On_Flag_Default(38) == sym%YES) Ref_Mean_ChI2 = Missing_Value_Real4
      if (Chan_On_Flag_Default(39) == sym%YES) Ref_ChI3 = Missing_Value_Real4
      if (Chan_On_Flag_Default(40) == sym%YES) Bt_ChI4 = Missing_Value_Real4
      if (Chan_On_Flag_Default(41) == sym%YES) Bt_ChI5 = Missing_Value_Real4
      if (Chan_On_Flag_Default(41) == sym%YES) Bt_Max_ChI5 = Missing_Value_Real4
      if (Chan_On_Flag_Default(41) == sym%YES) Bt_Min_ChI5 = Missing_Value_Real4
      if (Chan_On_Flag_Default(41) == sym%YES) Bt_Mean_ChI5 = Missing_Value_Real4
      if (Chan_On_Flag_Default(42) == sym%YES) Ref_ChDNB_Lunar_Mean_3x3 = Missing_Value_Real4
      if (Chan_On_Flag_Default(42) == sym%YES) Ref_ChDNB_Lunar_Max_3x3 = Missing_Value_Real4
      if (Chan_On_Flag_Default(42) == sym%YES) Ref_ChDNB_Lunar_Min_3x3 = Missing_Value_Real4
      if (Chan_On_Flag_Default(42) == sym%YES) Ref_ChDNB_Lunar_Std_3x3 = Missing_Value_Real4
end subroutine RESET_EXTRA_CHANNEL_ARRAYS

subroutine DESTROY_EXTRA_CHANNEL_ARRAYS
  if (allocated(Ref_ChI1)) deallocate(Ref_ChI1)
  if (allocated(Ref_Max_ChI1)) deallocate(Ref_Max_ChI1)
  if (allocated(Ref_Min_ChI1)) deallocate(Ref_Min_ChI1)
  if (allocated(Ref_Uni_ChI1)) deallocate(Ref_Uni_ChI1)
  if (allocated(Ref_Mean_ChI1)) deallocate(Ref_Mean_ChI1)
  if (allocated(Ref_ChI2)) deallocate(Ref_ChI2)
  if (allocated(Ref_Max_ChI2)) deallocate(Ref_Max_ChI2)
  if (allocated(Ref_Min_ChI2)) deallocate(Ref_Min_ChI2)
  if (allocated(Ref_Mean_ChI2)) deallocate(Ref_Mean_ChI2)
  if (allocated(Ref_ChI3)) deallocate(Ref_ChI3)
  if (allocated(Bt_ChI4)) deallocate(Bt_ChI4)
  if (allocated(Bt_ChI5)) deallocate(Bt_ChI5)
  if (allocated(Bt_Max_ChI5)) deallocate(Bt_Max_ChI5)
  if (allocated(Bt_Min_ChI5)) deallocate(Bt_Min_ChI5)
  if (allocated(Bt_Mean_ChI5)) deallocate(Bt_Mean_ChI5)
  if (allocated(Ref_ChDNB_Lunar_Mean_3x3)) deallocate(Ref_ChDNB_Lunar_Mean_3x3)
  if (allocated(Ref_ChDNB_Lunar_Min_3x3)) deallocate(Ref_ChDNB_Lunar_Min_3x3)
  if (allocated(Ref_ChDNB_Lunar_Max_3x3)) deallocate(Ref_ChDNB_Lunar_Max_3x3)
  if (allocated(Ref_ChDNB_Lunar_Std_3x3)) deallocate(Ref_ChDNB_Lunar_Std_3x3)
end subroutine DESTROY_EXTRA_CHANNEL_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_LUNAR_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Chan_On_Flag_Default(42) == sym%YES) then
           allocate(Lunzen(dim1,dim2))
           allocate(Lunaz(dim1,dim2))
           allocate(LunRelaz(dim1,dim2))
           allocate(Scatangle_Lunar(dim1,dim2))
           allocate(Glintzen_Lunar(dim1,dim2))
   endif
end subroutine CREATE_LUNAR_ARRAYS
subroutine RESET_LUNAR_ARRAYS()
      if (Chan_On_Flag_Default(42) == sym%YES) Lunzen = Missing_Value_Real4
      if (Chan_On_Flag_Default(42) == sym%YES) Lunaz = Missing_Value_Real4
      if (Chan_On_Flag_Default(42) == sym%YES) LunRelaz = Missing_Value_Real4
      if (Chan_On_Flag_Default(42) == sym%YES) Moon_Phase_Angle = Missing_Value_Real4
      if (Chan_On_Flag_Default(42) == sym%YES) Moon_Illum_Frac = Missing_Value_Real4
      if (allocated(Scatangle_Lunar)) Scatangle_Lunar = Missing_Value_Real4
      if (allocated(Glintzen_Lunar)) Glintzen_Lunar = Missing_Value_Real4
end subroutine RESET_LUNAR_ARRAYS
subroutine DESTROY_LUNAR_ARRAYS()
  if (allocated(Lunzen))deallocate(Lunzen)
  if (allocated(Lunaz)) deallocate(Lunaz)
  if (allocated(LunRelaz)) deallocate(LunRelaz)
  if (allocated(Scatangle_Lunar)) deallocate(Scatangle_Lunar)
  if (allocated(Glintzen_Lunar)) deallocate(Glintzen_Lunar)
end subroutine DESTROY_LUNAR_ARRAYS

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_BTD_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Chan_On_Flag_Default(20) == sym%YES .and. Chan_On_Flag_Default(31) == sym%YES) then
      allocate(Btd_Ch20_Ch31(dim1,dim2))
   endif
   if (Chan_On_Flag_Default(20) == sym%YES .and. Chan_On_Flag_Default(32) == sym%YES) then
      allocate(Btd_Ch20_Ch32(dim1,dim2))
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(32) == sym%YES) then
      allocate(Btd_Ch31_Ch32(dim1,dim2))
      allocate(Btd_Ch31_Ch32_Mean_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch32_Max_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch32_Min_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch32_Std_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch32_Bt_Ch31_Max_3x3(dim1,dim2))
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(27) == sym%YES) then
      allocate(Btd_Ch31_Ch27_Mean_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch27_Max_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch27_Min_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch27_Std_3x3(dim1,dim2))
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(29) == sym%YES) then
      allocate(Btd_Ch31_Ch29_Mean_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch29_Max_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch29_Min_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch29_Std_3x3(dim1,dim2))
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(33) == sym%YES) then
      allocate(Btd_Ch31_Ch33_Mean_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch33_Max_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch33_Min_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch33_Std_3x3(dim1,dim2))
   endif
end subroutine CREATE_BTD_ARRAYS
subroutine RESET_BTD_ARRAYS()
   if (Chan_On_Flag_Default(20) == sym%YES .and. Chan_On_Flag_Default(31) == sym%YES) then
      Btd_Ch20_Ch31 = Missing_Value_Real4
   endif
   if (Chan_On_Flag_Default(20) == sym%YES .and. Chan_On_Flag_Default(32) == sym%YES) then
      Btd_Ch20_Ch32 = Missing_Value_Real4
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(32) == sym%YES) then
      Btd_Ch31_Ch32 = Missing_Value_Real4
      Btd_Ch31_Ch32_Mean_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch32_Max_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch32_Min_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch32_Std_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch32_Bt_Ch31_Max_3x3 = Missing_Value_Real4
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(27) == sym%YES) then
      Btd_Ch31_Ch27_Mean_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch27_Max_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch27_Min_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch27_Std_3x3 = Missing_Value_Real4
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(29) == sym%YES) then
      Btd_Ch31_Ch29_Mean_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch29_Max_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch29_Min_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch29_Std_3x3 = Missing_Value_Real4
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(33) == sym%YES) then
      Btd_Ch31_Ch33_Mean_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch33_Max_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch33_Min_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch33_Std_3x3 = Missing_Value_Real4
   endif
end subroutine RESET_BTD_ARRAYS
subroutine DESTROY_BTD_ARRAYS()
   if (Chan_On_Flag_Default(20) == sym%YES .and. Chan_On_Flag_Default(31) == sym%YES) then
      deallocate(Btd_Ch20_Ch31)
   endif
   if (Chan_On_Flag_Default(20) == sym%YES .and. Chan_On_Flag_Default(32) == sym%YES) then
      deallocate(Btd_Ch20_Ch32)
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(32) == sym%YES) then
      deallocate(Btd_Ch31_Ch32)
      deallocate(Btd_Ch31_Ch32_Mean_3x3)
      deallocate(Btd_Ch31_Ch32_Max_3x3)
      deallocate(Btd_Ch31_Ch32_Min_3x3)
      deallocate(Btd_Ch31_Ch32_Std_3x3)
      deallocate(Btd_Ch31_Ch32_Bt_Ch31_Max_3x3)
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(27) == sym%YES) then
      deallocate(Btd_Ch31_Ch27_Mean_3x3)
      deallocate(Btd_Ch31_Ch27_Max_3x3)
      deallocate(Btd_Ch31_Ch27_Min_3x3)
      deallocate(Btd_Ch31_Ch27_Std_3x3)
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(29) == sym%YES) then
      deallocate(Btd_Ch31_Ch29_Mean_3x3)
      deallocate(Btd_Ch31_Ch29_Max_3x3)
      deallocate(Btd_Ch31_Ch29_Min_3x3)
      deallocate(Btd_Ch31_Ch29_Std_3x3)
   endif
   if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(33) == sym%YES) then
      deallocate(Btd_Ch31_Ch33_Mean_3x3)
      deallocate(Btd_Ch31_Ch33_Max_3x3)
      deallocate(Btd_Ch31_Ch33_Min_3x3)
      deallocate(Btd_Ch31_Ch33_Std_3x3)
   endif
end subroutine DESTROY_BTD_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_SURFACE_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(Land(dim1,dim2))
   allocate(Land_Modified(dim1,dim2))
   allocate(Land_Mask(dim1,dim2))
   allocate(Coast(dim1,dim2))
   allocate(Coast_Mask(dim1,dim2))
   allocate(Coast_Mask_Nwp(dim1,dim2))
   allocate(Glint_Mask(dim1,dim2))
   allocate(Glint_Mask_Lunar(dim1,dim2))
   allocate(Desert_Mask(dim1,dim2))
   allocate(Snow_Hires(dim1,dim2))
   allocate(Snow_Glob(dim1,dim2))
   allocate(Snow(dim1,dim2))
   allocate(Sfc_Type(dim1,dim2))
   allocate(Zsfc(dim1,dim2))
   allocate(Zsfc_Hires(dim1,dim2))
   allocate(Zsfc_Mean_3x3(dim1,dim2))
   allocate(Zsfc_Max_3x3(dim1,dim2))
   allocate(Zsfc_Min_3x3(dim1,dim2))
   allocate(Zsfc_Std_3x3(dim1,dim2))
end subroutine CREATE_SURFACE_ARRAYS
subroutine RESET_SURFACE_ARRAYS
   Land = Missing_Value_Int1
   Land_Modified = Missing_Value_Int1
   Land_Mask = Missing_Value_Int1
   Coast = Missing_Value_Int1
   Coast_Mask = Missing_Value_Int1
   Coast_Mask_Nwp = Missing_Value_Int1
   Glint_Mask = Missing_Value_Int1
   Glint_Mask_Lunar = Missing_Value_Int1
   Desert_Mask = Missing_Value_Int1
   Snow_Hires = Missing_Value_Int1
   Snow_Glob = Missing_Value_Int1
   Snow = Missing_Value_Int1
   Sfc_Type = Missing_Value_Int1
   Zsfc = Missing_Value_Real4
   Zsfc_Hires = Missing_Value_Real4
   Zsfc_Mean_3x3 = Missing_Value_Real4
   Zsfc_Max_3x3 = Missing_Value_Real4
   Zsfc_Min_3x3 = Missing_Value_Real4
   Zsfc_Std_3x3 = Missing_Value_Real4
end subroutine RESET_SURFACE_ARRAYS
subroutine DESTROY_SURFACE_ARRAYS
   deallocate(Land)
   deallocate(Land_Modified)
   deallocate(Land_Mask)
   deallocate(Coast)
   deallocate(Coast_Mask)
   deallocate(Coast_Mask_Nwp)
   deallocate(Glint_Mask)
   deallocate(Glint_Mask_Lunar)
   deallocate(Desert_Mask)
   deallocate(Snow_Hires)
   deallocate(Snow_Glob)
   deallocate(Snow)
   deallocate(Sfc_Type)
   deallocate(Zsfc)
   deallocate(Zsfc_Hires)
   deallocate(Zsfc_Mean_3x3)
   deallocate(Zsfc_Max_3x3)
   deallocate(Zsfc_Min_3x3)
   deallocate(Zsfc_Std_3x3)
end subroutine DESTROY_SURFACE_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_ACHA_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES) then
      allocate(ACHA_Processing_Order_Global(dim1,dim2)) 
      allocate(ACHA_Inversion_Flag_Global(dim1,dim2)) 
      allocate(Tc_ACHA(dim1,dim2))
      allocate(Ec_ACHA(dim1,dim2))
      allocate(Beta_ACHA(dim1,dim2))
      allocate(Pc_ACHA(dim1,dim2))
      allocate(Zc_ACHA(dim1,dim2))
      allocate(Zc_Top_ACHA(dim1,dim2))
      allocate(Zc_Base_ACHA(dim1,dim2))
      allocate(Cld_Layer_ACHA(dim1,dim2))
      allocate(Tau_ACHA(dim1,dim2))
      allocate(Reff_ACHA(dim1,dim2))
      allocate(Tc_ACHA_Uncertainty(dim1,dim2))
      allocate(Ec_ACHA_Uncertainty(dim1,dim2))
      allocate(Beta_ACHA_Uncertainty(dim1,dim2))
      allocate(Zc_ACHA_Uncertainty(dim1,dim2))
      allocate(Pc_ACHA_Uncertainty(dim1,dim2))
      allocate(ACHA_Quality_Flag(dim1,dim2))
      allocate(ACHA_Packed_Meta_Data_Flags(dim1,dim2))
      allocate(ACHA_Packed_Quality_Flags(dim1,dim2))
      allocate(Meta_Data_ACHA(dim1,dim2))
      allocate(ACHA_OE_Quality_Flags(3,dim1,dim2))
      allocate(Alt_Acha(dim1,dim2))
   endif
end subroutine CREATE_ACHA_ARRAYS
subroutine RESET_ACHA_ARRAYS()
   if (Cld_Flag == sym%YES) then
      ACHA_Processing_Order_Global = Missing_Value_Int1
      ACHA_Inversion_Flag_Global = Missing_Value_Int1
      Tc_ACHA = Missing_Value_Real4
      Ec_ACHA = Missing_Value_Real4
      Beta_ACHA = Missing_Value_Real4
      Pc_ACHA = Missing_Value_Real4
      Zc_ACHA = Missing_Value_Real4
      Zc_Top_ACHA = Missing_Value_Real4
      Zc_Base_ACHA = Missing_Value_Real4
      Cld_Layer_ACHA = Missing_Value_Int1
      Tau_ACHA = Missing_Value_Real4
      Reff_ACHA = Missing_Value_Real4
      Tc_ACHA_Uncertainty = Missing_Value_Real4
      Ec_ACHA_Uncertainty = Missing_Value_Real4
      Beta_ACHA_Uncertainty = Missing_Value_Real4
      Zc_ACHA_Uncertainty = Missing_Value_Real4
      Pc_ACHA_Uncertainty = Missing_Value_Real4
      ACHA_Quality_Flag = 0
      ACHA_Packed_Meta_Data_Flags = 0
      ACHA_Packed_Quality_Flags = 0
      Meta_Data_ACHA = 0
      ACHA_OE_Quality_Flags = 0
      Alt_Acha = Missing_Value_Real4
   endif
end subroutine RESET_ACHA_ARRAYS
subroutine DESTROY_ACHA_ARRAYS()
   if (Cld_Flag == sym%YES) then
      deallocate(ACHA_Processing_Order_Global) 
      deallocate(ACHA_Inversion_Flag_Global) 
      deallocate(Tc_ACHA)
      deallocate(Ec_ACHA)
      deallocate(Beta_ACHA)
      deallocate(Pc_ACHA)
      deallocate(Zc_ACHA)
      deallocate(Zc_Top_ACHA)
      deallocate(Zc_Base_ACHA)
      deallocate(Cld_Layer_ACHA)
      deallocate(Tau_ACHA)
      deallocate(Reff_ACHA)
      deallocate(Tc_ACHA_Uncertainty)
      deallocate(Ec_ACHA_Uncertainty)
      deallocate(Beta_ACHA_Uncertainty)
      deallocate(Zc_ACHA_Uncertainty)
      deallocate(Pc_ACHA_Uncertainty)
      deallocate(ACHA_Quality_Flag)
      deallocate(ACHA_Packed_Meta_Data_Flags)
      deallocate(ACHA_Packed_Quality_Flags)
      deallocate(Meta_Data_ACHA)
      deallocate(ACHA_OE_Quality_Flags)
      deallocate(Alt_Acha)
   endif
end subroutine DESTROY_ACHA_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_DCOMP_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES) then
      allocate(Tau_DCOMP(dim1,dim2))
      allocate(Tau_DCOMP_Ap(dim1,dim2))
      allocate(Vis_Ref_Fm(dim1,dim2))
      allocate(Reff_DCOMP(dim1,dim2))
      allocate(Lwp_DCOMP(dim1,dim2))
      allocate(Iwp_DCOMP(dim1,dim2))
      allocate(Iwp_Tau_DCOMP(dim1,dim2))
      allocate(Cwp_DCOMP(dim1,dim2))
      allocate(Cwp_Ice_Layer_DCOMP(dim1,dim2))
      allocate(Cwp_Water_Layer_DCOMP(dim1,dim2))
      allocate(Cwp_Scwater_Layer_DCOMP(dim1,dim2))
      allocate(Rain_Rate_DCOMP(dim1,dim2))
      allocate(H_DCOMP(dim1,dim2))
      allocate(N_DCOMP(dim1,dim2))
      allocate(Tau_DCOMP_Cost(dim1,dim2))
      allocate(Reff_DCOMP_Cost(dim1,dim2))
      allocate(Tau_DCOMP_Qf(dim1,dim2))
      allocate(Reff_DCOMP_Qf(dim1,dim2))
      allocate(DCOMP_Quality_Flag(dim1,dim2))
      allocate(DCOMP_Info_Flag(dim1,dim2))
      allocate(Cloud_063um_Albedo(dim1,dim2))
      allocate(Cloud_063um_Spherical_Albedo(dim1,dim2))
      allocate(Cloud_063um_Transmission_View(dim1,dim2))
      allocate(Cloud_063um_Transmission_Solar(dim1,dim2))
      allocate(Insolation_DCOMP(dim1,dim2))
      allocate(Insolation_DCOMP_Diffuse(dim1,dim2))
      allocate(Cost_DCOMP(dim1,dim2))
      allocate(Error_Cov_Matrix_Cod(dim1,dim2))
      allocate(Error_Cov_Matrix_Ref(dim1,dim2))
      allocate(DCOMP_Diag_1(dim1,dim2))
      allocate(DCOMP_Diag_2(dim1,dim2))
      allocate(DCOMP_Diag_3(dim1,dim2))
      allocate(DCOMP_Diag_4(dim1,dim2))
      allocate(DCOMP_Diag_Wv1(dim1,dim2))
      allocate(DCOMP_Diag_Wv2(dim1,dim2))
      allocate(DCOMP_Diag_Toc_Rfl1(dim1,dim2))
      allocate(DCOMP_Diag_Toc_Rfl2(dim1,dim2))
      allocate(DCOMP_Diag_Virt_Alb1(dim1,dim2))
      allocate(DCOMP_Diag_Virt_Alb2(dim1,dim2))
      allocate(DCOMP_Diag_Toc_Rfl_Unc1(dim1,dim2))
      allocate(DCOMP_Diag_Toc_Rfl_Unc2(dim1,dim2))
   endif
end subroutine CREATE_DCOMP_ARRAYS
subroutine RESET_DCOMP_ARRAYS()
   if (Cld_Flag == sym%YES) then
      Tau_DCOMP = Missing_Value_Real4
      Tau_DCOMP_Ap = Missing_Value_Real4
      Vis_Ref_Fm = Missing_Value_Real4
      Reff_DCOMP = Missing_Value_Real4
      Lwp_DCOMP = Missing_Value_Real4
      Iwp_DCOMP = Missing_Value_Real4
      Iwp_Tau_DCOMP = Missing_Value_Real4
      Cwp_DCOMP = Missing_Value_Real4
      Cwp_Ice_Layer_DCOMP = Missing_Value_Real4
      Cwp_Water_Layer_DCOMP = Missing_Value_Real4
      Cwp_Scwater_Layer_DCOMP = Missing_Value_Real4
      Rain_Rate_DCOMP = Missing_Value_Real4
      H_DCOMP = Missing_Value_Real4
      N_DCOMP = Missing_Value_Real4
      Tau_DCOMP_Cost = Missing_Value_Real4
      Reff_DCOMP_Cost = Missing_Value_Real4
      Tau_DCOMP_Qf = Missing_Value_Int1
      Reff_DCOMP_Qf = Missing_Value_Int1
      DCOMP_Quality_Flag = 0
      DCOMP_Info_Flag = 0
      Cloud_063um_Albedo = Missing_Value_Real4
      Cloud_063um_Spherical_Albedo = Missing_Value_Real4
      Cloud_063um_Transmission_View = Missing_Value_Real4
      Cloud_063um_Transmission_Solar = Missing_Value_Real4
      Insolation_DCOMP = Missing_Value_Real4
      Insolation_DCOMP_Diffuse = Missing_Value_Real4
      Cost_DCOMP = Missing_Value_Real4
      Error_Cov_Matrix_Cod = Missing_Value_Real4
      Error_Cov_Matrix_Ref = Missing_Value_Real4
      DCOMP_Diag_1 = Missing_Value_Real4
      DCOMP_Diag_2 = Missing_Value_Real4
      DCOMP_Diag_3 = Missing_Value_Real4
      DCOMP_Diag_4 = Missing_Value_Real4
      DCOMP_Diag_Wv1 = Missing_Value_Real4
      DCOMP_Diag_Wv2 = Missing_Value_Real4
      DCOMP_Diag_Toc_Rfl1 = Missing_Value_Real4
      DCOMP_Diag_Toc_Rfl2 = Missing_Value_Real4
      DCOMP_Diag_Virt_Alb1 = Missing_Value_Real4
      DCOMP_Diag_Virt_Alb2 = Missing_Value_Real4
      DCOMP_Diag_Toc_Rfl_Unc1 = Missing_Value_Real4
      DCOMP_Diag_Toc_Rfl_Unc2 = Missing_Value_Real4
   endif
end subroutine RESET_DCOMP_ARRAYS
subroutine DESTROY_DCOMP_ARRAYS()
   if (Cld_Flag == sym%YES) then
      deallocate(Tau_DCOMP)
      deallocate(Tau_DCOMP_Ap)
      deallocate(Vis_Ref_Fm)
      deallocate(Reff_DCOMP)
      deallocate(Lwp_DCOMP)
      deallocate(Iwp_DCOMP)
      deallocate(Iwp_Tau_DCOMP)
      deallocate(Cwp_DCOMP)
      deallocate(Cwp_Ice_Layer_DCOMP)
      deallocate(Cwp_Water_Layer_DCOMP)
      deallocate(Cwp_Scwater_Layer_DCOMP)
      deallocate(Rain_Rate_DCOMP)
      deallocate(H_DCOMP)
      deallocate(N_DCOMP)
      deallocate(Tau_DCOMP_Cost)
      deallocate(Reff_DCOMP_Cost)
      deallocate(Tau_DCOMP_Qf)
      deallocate(Reff_DCOMP_Qf)
      deallocate(DCOMP_Quality_Flag)
      deallocate(DCOMP_Info_Flag)
      deallocate(Cloud_063um_Albedo)
      deallocate(Cloud_063um_Spherical_Albedo)
      deallocate(Cloud_063um_Transmission_View)
      deallocate(Cloud_063um_Transmission_Solar)
      deallocate(Insolation_DCOMP)
      deallocate(Insolation_DCOMP_Diffuse)
      deallocate(Cost_DCOMP)
      deallocate(Error_Cov_Matrix_Cod)
      deallocate(Error_Cov_Matrix_Ref)
      deallocate(DCOMP_Diag_1)
      deallocate(DCOMP_Diag_2)
      deallocate(DCOMP_Diag_3)
      deallocate(DCOMP_Diag_4)
      deallocate(DCOMP_Diag_Wv1)
      deallocate(DCOMP_Diag_Wv2)
      deallocate(DCOMP_Diag_Toc_Rfl1)
      deallocate(DCOMP_Diag_Toc_Rfl2)
      deallocate(DCOMP_Diag_Virt_Alb1)
      deallocate(DCOMP_Diag_Virt_Alb2)
      deallocate(DCOMP_Diag_Toc_Rfl_Unc1)
      deallocate(DCOMP_Diag_Toc_Rfl_Unc2)
   endif
end subroutine DESTROY_DCOMP_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_SASRAB_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Erb_Flag == sym%YES) then
      allocate(Insolation_All_Sky(dim1,dim2))
      allocate(Insolation_All_Sky_Diffuse(dim1,dim2))
      allocate(Insolation_Clear_Sky(dim1,dim2))
      allocate(Insolation_Cld_Opd(dim1,dim2))
      allocate(Insolation_Aer_Opd(dim1,dim2))
   endif
end subroutine CREATE_SASRAB_ARRAYS
subroutine RESET_SASRAB_ARRAYS()
   if (Erb_Flag == sym%YES) then
      Insolation_All_Sky = Missing_Value_Real4
      Insolation_All_Sky_Diffuse = Missing_Value_Real4
      Insolation_Clear_Sky = Missing_Value_Real4
      Insolation_Cld_Opd = Missing_Value_Real4
      Insolation_Aer_Opd = Missing_Value_Real4
   endif
end subroutine RESET_SASRAB_ARRAYS
subroutine DESTROY_SASRAB_ARRAYS()
   if (Erb_Flag == sym%YES) then
      deallocate(Insolation_All_Sky)
      deallocate(Insolation_All_Sky_Diffuse)
      deallocate(Insolation_Clear_Sky)
      deallocate(Insolation_Cld_Opd)
      deallocate(Insolation_Aer_Opd)
   endif
end subroutine DESTROY_SASRAB_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_OLR_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Erb_Flag == sym%YES) then
      allocate(Olr(dim1,dim2))
      allocate(Olr_Qf(dim1,dim2))
   endif
end subroutine CREATE_OLR_ARRAYS
subroutine RESET_OLR_ARRAYS
   if (Erb_Flag == sym%YES) then
      Olr = Missing_Value_Real4
      Olr_Qf = Missing_Value_Int1
   endif
end subroutine RESET_OLR_ARRAYS
subroutine DESTROY_OLR_ARRAYS
   if (Erb_Flag == sym%YES) then
      deallocate(Olr)
      deallocate(Olr_Qf)
   endif
end subroutine DESTROY_OLR_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_AEROSOL_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  allocate(Aot_Qf(dim1,dim2))
  if (Aer_Flag == sym%YES) then
     allocate(Aot1(dim1,dim2))
     allocate(Aot2(dim1,dim2))
     allocate(Aot3a(dim1,dim2))
  endif
end subroutine CREATE_AEROSOL_ARRAYS 
subroutine RESET_AEROSOL_ARRAYS()
   Aot_Qf = 0
   if (Aer_Flag == sym%YES) then
    Aot1 = Missing_Value_Real4
    Aot2 = Missing_Value_Real4
    Aot3a = Missing_Value_Real4
   endif
end subroutine RESET_AEROSOL_ARRAYS 
subroutine DESTROY_AEROSOL_ARRAYS()
  deallocate(Aot_Qf)
  if (Aer_Flag == sym%YES) then
     deallocate(Aot1)
     deallocate(Aot2)
     deallocate(Aot3a)
  endif
end subroutine DESTROY_AEROSOL_ARRAYS 
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_CLOUD_MASK_ARRAYS(dim1,dim2,dim3)
  integer, intent(in):: dim1, dim2, dim3
  allocate(Cld_Mask_Qf(dim1,dim2))
  if (Cld_Flag == sym%YES) then
     allocate(Cld_Mask_Aux(dim1,dim2))
     allocate(Cld_Test_Vector_Packed(dim3,dim1,dim2))
     allocate(Cld_Mask(dim1,dim2))
     allocate(Adj_Pix_Cld_Mask(dim1,dim2))
     allocate(Posterior_Cld_Probability(dim1,dim2))
     allocate(Bayes_Mask_Sfc_Type_Global(dim1,dim2))
  endif
end subroutine CREATE_CLOUD_MASK_ARRAYS
subroutine RESET_CLOUD_MASK_ARRAYS()
  Cld_Mask_Qf = Missing_Value_Int1
  if (Cld_Flag == sym%YES) then
     Cld_Mask = Missing_Value_Int1
     Cld_Mask_Aux = Missing_Value_Int1
     Adj_Pix_Cld_Mask = Missing_Value_Int1
     Posterior_Cld_Probability = Missing_Value_Real4
     Cld_Test_Vector_Packed = 0
     Bayes_Mask_Sfc_Type_Global = Missing_Value_Int1
  endif
end subroutine RESET_CLOUD_MASK_ARRAYS
subroutine DESTROY_CLOUD_MASK_ARRAYS()
  deallocate(Cld_Mask_Qf)
  if (Cld_Flag == sym%YES) then
     deallocate(Cld_Mask_Aux)
     deallocate(Cld_Test_Vector_Packed)
     deallocate(Cld_Mask)
     deallocate(Adj_Pix_Cld_Mask)
     deallocate(Posterior_Cld_Probability)
     deallocate(Bayes_Mask_Sfc_Type_Global)
  endif
end subroutine DESTROY_CLOUD_MASK_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_CLOUD_TYPE_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  if (Cld_Flag == sym%YES) then
     allocate(Cld_Type_Aux(dim1,dim2))
     allocate(Cld_Phase_Aux(dim1,dim2))
     allocate(Cld_Phase(dim1,dim2))
     allocate(Cld_Type(dim1,dim2))
     allocate(Cirrus_Quality(dim1,dim2))
  endif
end subroutine CREATE_CLOUD_TYPE_ARRAYS
subroutine RESET_CLOUD_TYPE_ARRAYS()
  if (Cld_Flag == sym%YES) then
      Cld_Phase = Missing_Value_Int1
      Cld_Type = Missing_Value_Int1
      Cld_Phase_Aux = Missing_Value_Int1
      Cld_Type_Aux = Missing_Value_Int1
      Cirrus_Quality = 1
  endif
end subroutine RESET_CLOUD_TYPE_ARRAYS
subroutine DESTROY_CLOUD_TYPE_ARRAYS
  if (Cld_Flag == sym%YES) then
     deallocate(Cld_Type_Aux)
     deallocate(Cld_Phase_Aux)
     deallocate(Cld_Phase)
     deallocate(Cld_Type)
     deallocate(Cirrus_Quality)
  endif
end subroutine DESTROY_CLOUD_TYPE_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_DIAGNOSTIC_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  allocate(Temp_Pix_Array(dim1,dim2))
  allocate(Temp_Pix_Array_2(dim1,dim2))
  allocate(Diag_Pix_Array_1(dim1,dim2))
  allocate(Diag_Pix_Array_2(dim1,dim2))
  allocate(Diag_Pix_Array_3(dim1,dim2))
  allocate(Missing_Pixel_Array_Real4(dim1,dim2))
  allocate(One_Byte_Temp(dim1,dim2))
  allocate(Two_Byte_Temp(dim1,dim2))
end subroutine CREATE_DIAGNOSTIC_ARRAYS
subroutine RESET_DIAGNOSTIC_ARRAYS()
  Temp_Pix_Array = Missing_Value_Real4
  Temp_Pix_Array_2 = Missing_Value_Real4
  Diag_Pix_Array_1 = Missing_Value_Real4
  Diag_Pix_Array_2 = Missing_Value_Real4
  Diag_Pix_Array_3 = Missing_Value_Real4
  Missing_Pixel_Array_Real4 = Missing_Value_Real4
  One_Byte_Temp = Missing_Value_Int1
  Two_Byte_Temp = Missing_Value_Int2
end subroutine RESET_DIAGNOSTIC_ARRAYS
subroutine DESTROY_DIAGNOSTIC_ARRAYS()
  deallocate(Temp_Pix_Array)
  deallocate(Temp_Pix_Array_2)
  deallocate(Diag_Pix_Array_1)
  deallocate(Diag_Pix_Array_2)
  deallocate(Diag_Pix_Array_3)
  deallocate(Missing_Pixel_Array_Real4)
  deallocate(One_Byte_Temp)
  deallocate(Two_Byte_Temp)
end subroutine DESTROY_DIAGNOSTIC_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_SFC_PROD_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  allocate(Tsfc_Retrieved(dim1,dim2))
  allocate(Trad_Retrieved(dim1,dim2))
  allocate(Tsfc_Qf(dim1,dim2))
  allocate(Ndsi_Toa(dim1,dim2))
  allocate(Ndvi_Toa(dim1,dim2))
  allocate(Ndvi_Qf(dim1,dim2))
  allocate(Ndvi_Sfc(dim1,dim2))
  allocate(Ndvi_Sfc_White_Sky(dim1,dim2))
  allocate(Rsr(dim1,dim2))
  allocate(Rsr_Qf(dim1,dim2))
  allocate(Sst_Unmasked(dim1,dim2))
  allocate(Sst_Masked(dim1,dim2))
end subroutine CREATE_SFC_PROD_ARRAYS
subroutine RESET_SFC_PROD_ARRAYS()
  Tsfc_Retrieved = Missing_Value_Real4
  Trad_Retrieved = Missing_Value_Real4
  Tsfc_Qf = Missing_Value_Int1
  Ndsi_Toa = Missing_Value_Real4
  Ndvi_Toa = Missing_Value_Real4
  Ndvi_Qf = Missing_Value_Int1
  Ndvi_Sfc = Missing_Value_Real4
  Ndvi_Sfc_White_Sky = Missing_Value_Real4
  Rsr = Missing_Value_Real4
  Rsr_Qf = Missing_Value_Int1
  Sst_Unmasked = Missing_Value_Real4
  Sst_Masked = Missing_Value_Real4
end subroutine RESET_SFC_PROD_ARRAYS
subroutine DESTROY_SFC_PROD_ARRAYS()
  if (allocated(Tsfc_Retrieved))  deallocate(Tsfc_Retrieved)
  if (allocated(Trad_Retrieved))  deallocate(Trad_Retrieved)
  if (allocated(Tsfc_Qf))  deallocate(Tsfc_Qf)
  if (allocated(Ndsi_Toa)) deallocate(Ndsi_Toa)
  if (allocated(Ndvi_Toa)) deallocate(Ndvi_Toa)
  if (allocated(Ndvi_Qf)) deallocate(Ndvi_Qf)
  if (allocated(Ndvi_Sfc)) deallocate(Ndvi_Sfc)
  if (allocated(Ndvi_Sfc_White_Sky)) deallocate(Ndvi_Sfc_White_Sky)
  if (allocated(Rsr)) deallocate(Rsr)
  if (allocated(Rsr_Qf)) deallocate(Rsr_Qf)
  if (allocated(Sst_Unmasked)) deallocate(Sst_Unmasked)
  if (allocated(Sst_Masked)) deallocate(Sst_Masked)
end subroutine DESTROY_SFC_PROD_ARRAYS
!-----------------------------------------------------------
! 
!-----------------------------------------------------------
subroutine CREATE_CLOUD_PROD_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  if (Cld_Flag == sym%YES) then
    allocate(Pc_Opaque_Cloud(dim1,dim2))
    allocate(Zc_Opaque_Cloud(dim1,dim2))
    allocate(Tc_Opaque_Cloud(dim1,dim2))
    allocate(Pc_Lower_Cloud(dim1,dim2))
    allocate(Zc_Lower_Cloud(dim1,dim2))
    allocate(Tc_Lower_Cloud(dim1,dim2))
    allocate(Pc_H2O(dim1,dim2))
    allocate(Tc_H2O(dim1,dim2))
    allocate(Zc_H2O(dim1,dim2))
    allocate(Zclr_H2O_Peak(dim1,dim2))
    allocate(Cloud_Fraction_3x3(dim1,dim2))
    allocate(Cloud_Fraction_Uncer_3x3(dim1,dim2))
  endif
end subroutine CREATE_CLOUD_PROD_ARRAYS
subroutine RESET_CLOUD_PROD_ARRAYS()
  if (Cld_Flag == sym%YES) then
     Pc_Opaque_Cloud = Missing_Value_Real4
     Zc_Opaque_Cloud = Missing_Value_Real4
     Tc_Opaque_Cloud = Missing_Value_Real4
     Pc_Lower_Cloud = Missing_Value_Real4
     Zc_Lower_Cloud = Missing_Value_Real4
     Tc_Lower_Cloud = Missing_Value_Real4
     Pc_H2O = Missing_Value_Real4
     Tc_H2O = Missing_Value_Real4
     Zc_H2O = Missing_Value_Real4
     Zclr_H2O_Peak = Missing_Value_Real4
     Cloud_Fraction_3x3 = Missing_Value_Real4
     Cloud_Fraction_Uncer_3x3 = Missing_Value_Real4
  endif
end subroutine RESET_CLOUD_PROD_ARRAYS
subroutine DESTROY_CLOUD_PROD_ARRAYS()
  if (Cld_Flag == sym%YES) then
     deallocate(Pc_Opaque_Cloud)
     deallocate(Zc_Opaque_Cloud)
     deallocate(Tc_Opaque_Cloud)
     deallocate(Pc_Lower_Cloud)
     deallocate(Zc_Lower_Cloud)
     deallocate(Tc_Lower_Cloud)
     deallocate(Pc_H2O)
     deallocate(Tc_H2O)
     deallocate(Zc_H2O)
     deallocate(Zclr_H2O_Peak)
     deallocate(Cloud_Fraction_3x3)
     deallocate(Cloud_Fraction_Uncer_3x3)
  endif
end subroutine DESTROY_CLOUD_PROD_ARRAYS
!-----------------------------------------------------------
! 
!-----------------------------------------------------------
subroutine CREATE_NLCOMP_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  if (Cld_Flag == sym%YES) then
     allocate(Tau_Nlcomp(dim1,dim2))
     allocate(Reff_Nlcomp(dim1,dim2))
     allocate(Nlcomp_Info_Flag(dim1,dim2))
     allocate(Nlcomp_Quality_Flag(dim1,dim2))
     allocate(Tau_Nlcomp_Cost(dim1,dim2))
     allocate(Reff_Nlcomp_Cost(dim1,dim2))
  endif
end subroutine CREATE_NLCOMP_ARRAYS
subroutine RESET_NLCOMP_ARRAYS()
  if (Cld_Flag == sym%YES) then
      Tau_Nlcomp = Missing_Value_Real4
      Reff_Nlcomp = Missing_Value_Real4
      Tau_Nlcomp_Cost = Missing_Value_Real4
      Reff_Nlcomp_Cost = Missing_Value_Real4
      Nlcomp_Quality_Flag = 0
      Nlcomp_Info_Flag = 0
  endif
end subroutine RESET_NLCOMP_ARRAYS
subroutine DESTROY_NLCOMP_ARRAYS()
  if (Cld_Flag == sym%YES) then
     deallocate(Tau_Nlcomp)
     deallocate(Reff_Nlcomp)
     deallocate(Nlcomp_Info_Flag)
     deallocate(Nlcomp_Quality_Flag)
     deallocate(Tau_Nlcomp_Cost)
     deallocate(Reff_Nlcomp_Cost)
  endif
end subroutine DESTROY_NLCOMP_ARRAYS
!-----------------------------------------------------------
! end of module
!-----------------------------------------------------------
end module PIXEL_COMMON