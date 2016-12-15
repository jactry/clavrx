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
!   CLAVR-x  Modis    Avhrr   ABI    AHI    VIIRS   Wavelength  Obs_Type
!     01       1       1      2      3      M5       0.659      solar
!     02       2       2      3      4      M7       0.865      solar
!     03       3       -      1      1      M3       0.470       solar
!     04       4       -      -      2      M4       0.555      solar
!     05       5       -      -      -      M8       1.240      solar
!     06       6       3a     5      5     M10       1.640      solar
!     07       7       -      6      6     M11       2.130      solar
!     08       8       -      -      -      M1       0.415      solar
!     09       9       -      -      -      M2       0.443      solar
!     10       10      -      -      -       -       0.490      solar
!     11       11      -      -      -       -       0.531      solar
!     12       12      -      -      -       -       0.565      solar
!     13       13      -      -      -       -       0.653      solar
!     14       14      -      -      -       -       0.681      solar
!     15       15      -      -      -      M6       0.750      solar
!     16       16      -      -      -       -       0.865      solar
!     17       17      -      -      -       -       0.905      solar
!     18       18      -      -      -       -       0.936      solar
!     19       19      -      -      -       -       0.940      solar
!     20       20      3b     7      7     M12       3.750      mixed
!     21       21      -      -      -       -       3.959      mixed
!     22       22      -      -      -     M13       3.959      mixed
!     23       23      -      -      -       -       4.050      therm
!     24       24      -      -      -       -       4.465      therm
!     25       25      -      -      -       -       4.515      therm
!     26       26      -      4      -      M9       1.375      solar
!     27       27      -      9      9       -       6.715      therm
!     28       28      -     10     10       -       7.325      therm
!     29       29      -     11     11     M14       8.550      therm
!     30       30      -     12     12       -       9.730      therm
!     31       31      4     14     14     M15      11.030      therm
!     32       32      5     15     15     M16      12.020      therm
!     33       32      -     16     16       -      13.335      therm
!     34       33      -      -      -       -      13.635      therm
!     35       34      -      -      -       -      13.935      therm
!     36       35      -      -      -       -      14.235      therm
!     37       36      -      8      8       -       6.200      therm
!     38       -       -     13     13       -      10.400      therm
!     39       -       -      -      -      I1       0.640      solar
!     40       -       -      -      -      I2       0.865      solar
!     41       -       -      -      -      I3       1.610      solar
!     42       -       -      -      -      I4       3.740      mixed
!     43       -       -      -      -      I5      11.450      therm
!     44       -       -      -      -     DNB       0.700      lunar
!     45*      -       -      -      -      -       13.335      therm
!
!     * = a pseudo 13.3 channel only AVHRR/HIRS and VIIRS/CRIS IFF
!
!  Description of variables in "ch" structure:
!
!  Rad_Toa = Observed Top of Atmosphere Radiance
!  Bt_Toa = Observed Top of Atmosphere Brightness Temperature
!  Rad_Toa_Clear = Simulated Top of Atmosphere Radiance under clear-sky
!  Rad_Atm = Simulated Radiance at Toa from atmospheric cloud-free emission
!  Rad_Atm_Dwn_Sfc = Simulated Downward Radiance at Sfc from atmospheric cloud-free emission
!  Trans_Atm = Simulated Transmission from Surface to Toa for cloud-free atmosphere
!              along viewing zenith angle path
!  Trans_Atm_Total = Simulated Transmission from Toa to Surface to Toa for cloud-free atmosphere
!                      along solar and viewing zenith angle path
!  Bt_Toa_Clear = Simulated Toa Brightness Temperature under clear-sky
!  Ref_Toa = Top of Atmosphere Reflectance
!  Ref_Toa_Unnorm = Top of Atmosphere Reflectance not normalized by cos(solzen)
!  Ref_Sfc = Observed Reflectance adjusted as if measured at surface level
!  Sfc_Ref_White_Sky - surface reflectance under diffuse illumination
!  Emiss_Tropo - emissity of cloud placed at Tropopause needed to match toa  radiance
!  Bt_Toa_Clear = Simulated Toa Reflectance under clear-sky
!  Unc = Uncertainty Flag (relevant only to MODIS)
!  Opd = Optical Depth
!  CSBT_Mask = Clear Sky Brighntess Temperature Mask
!  Opaque_Height = Maximum Height at which trans to space is 0.
!  Obs_Type = type of observation (solar,lunar,mixed,thermal).  Controls what is allocated
!--------------------------------------------------------------------------------------
module PIXEL_COMMON

  use CX_CONSTANTS_MOD, only: &
    real4 &
  , real8 &
  , int4 &
  , int2 &
  , int1 &
  , THERMAL_OBS_TYPE &
  , LUNAR_OBS_TYPE &
  , SOLAR_OBS_TYPE &
  , MIXED_OBS_TYPE &
  , Missing_Value_Int2 &
  , Missing_Value_Int4 &
  , Missing_Value_Int1 &
  , Missing_Value_Real4 &
  , sym &
  , nchan_clavrx &
  , Max_Num_Cld_Test_Bytes
  
  use CLAVRX_MESSAGE_MODULE, only: MESG, VERB_LEV
  implicit none
  private
  public:: CREATE_PIXEL_ARRAYS, &
           DESTROY_PIXEL_ARRAYS, &
           RESET_PIXEL_ARRAYS_TO_MISSING

  private:: CREATE_NAV_ARRAYS, RESET_NAV_ARRAYS, DESTROY_NAV_ARRAYS
  private:: CREATE_GEO_ARRAYS, RESET_GEO_ARRAYS, DESTROY_GEO_ARRAYS
  private:: CREATE_SENSOR_ARRAYS, RESET_SENSOR_ARRAYS, DESTROY_SENSOR_ARRAYS
  private:: CREATE_AVHRR_ANCHOR_ARRAYS, RESET_AVHRR_ANCHOR_ARRAYS, DESTROY_AVHRR_ANCHOR_ARRAYS
  private:: CREATE_NWP_PIX_ARRAYS, RESET_NWP_PIX_ARRAYS, DESTROY_NWP_PIX_ARRAYS
  private:: CREATE_REF_CHANNEL_ARRAYS, RESET_REF_CHANNEL_ARRAYS, DESTROY_REF_CHANNEL_ARRAYS
  private:: CREATE_THERM_CHANNEL_ARRAYS, RESET_THERM_CHANNEL_ARRAYS, DESTROY_THERM_CHANNEL_ARRAYS
  private:: CREATE_EXTRA_CHANNEL_ARRAYS, RESET_EXTRA_CHANNEL_ARRAYS, DESTROY_EXTRA_CHANNEL_ARRAYS
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
  character(len=1020),dimension(100), public, save:: Temporary_File_Name

  !---------------------------------------------------------------------------------
  ! CLAVR-x file list variables
  !---------------------------------------------------------------------------------
  type :: observations
    real, dimension(:,:), allocatable:: Rad_Toa
    real, dimension(:,:), allocatable:: Bt_Toa
    real, dimension(:,:), allocatable:: Rad_Toa_Clear
    real, dimension(:,:), allocatable:: Rad_Atm
    real, dimension(:,:), allocatable:: Rad_Atm_Dwn_Sfc
    real, dimension(:,:), allocatable:: Trans_Atm
    real, dimension(:,:), allocatable:: Trans_Atm_Total
    real, dimension(:,:), allocatable:: Bt_Toa_Clear
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
    real, dimension(:,:), allocatable:: Opd
    integer (kind=int1), dimension(:,:), allocatable:: Unc
    integer (kind=int1), dimension(:,:), allocatable:: CSBT_Mask
    real (kind=real4), dimension(:,:), allocatable:: Opaque_Height
    character(len=5) :: Obs_Type
  end type observations

  type :: geometry_definition
     real (kind=real4), dimension(:,:), allocatable:: Satzen
     real (kind=real4), dimension(:,:), allocatable:: Solzen
     real (kind=real4), dimension(:,:), allocatable:: Solaz
     real (kind=real4), dimension(:,:), allocatable:: Sataz
     real (kind=real4), dimension(:,:), allocatable:: Relaz
     real (kind=real4), dimension(:,:), allocatable:: Glintzen
     real (kind=real4), dimension(:,:), allocatable:: Seczen
     real (kind=real4), dimension(:,:), allocatable:: Coszen
     real (kind=real4), dimension(:,:), allocatable:: CosSolzen
     real (kind=real4), dimension(:,:), allocatable:: Scatangle
     real (kind=real4), dimension(:,:), allocatable:: Airmass
     real (kind=real4), dimension(:,:), allocatable:: Glintzen_Lunar
     real (kind=real4), dimension(:,:), allocatable:: Scatangle_Lunar
     real (kind=real4), dimension(:,:), allocatable:: Lunzen
     real (kind=real4), dimension(:,:), allocatable:: Lunaz
     real (kind=real4), dimension(:,:), allocatable:: LunRelaz
     double precision:: Moon_Phase_Angle
     real (kind=real4):: Moon_Illum_Frac
     real(kind=real4):: Solzen_Min_Limit
     real(kind=real4):: Solzen_Max_Limit
     real(kind=real4):: Satzen_Min_Limit
     real(kind=real4):: Satzen_Max_Limit
  end type geometry_definition

  type :: navigation_definition
   integer (kind=int1), dimension(:), allocatable:: Ascend
   integer (kind=int1), dimension(:,:), allocatable:: Sounder_Fov
   integer (kind=int1), dimension(:,:), allocatable:: Sounder_Fov_Mask
   integer (kind=int2), dimension(:,:), allocatable:: Sounder_X
   integer (kind=int2), dimension(:,:), allocatable:: Sounder_Y
   integer (kind=int2), dimension(:,:), allocatable:: Sounder_Fov_Segment_Idx
   real (kind=real4), dimension(:,:), allocatable:: Lat
   real (kind=real4), dimension(:,:), allocatable:: Lon
   real (kind=real4), dimension(:,:), allocatable:: Lat_Pc
   real (kind=real4), dimension(:,:), allocatable:: Lon_Pc
   real (kind=real4), dimension(:,:), allocatable:: Lat_1b
   real (kind=real4), dimension(:,:), allocatable:: Lon_1b
   real(kind=real4):: Lat_Min_Limit
   real(kind=real4):: Lat_Max_Limit
   real(kind=real4):: Lon_Min_Limit
   real(kind=real4):: Lon_Max_Limit
   real(kind=real4):: Timerr_Seconds
   logical :: lon_lat_limits_set
  end type navigation_definition

  type :: surface_definition
     integer(kind=int1), dimension(:,:), allocatable:: Land
     integer(kind=int1), dimension(:,:), allocatable:: Land_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Coast
     integer(kind=int1), dimension(:,:), allocatable:: Coast_Mask 
     integer(kind=int1), dimension(:,:), allocatable:: Coast_Mask_Nwp
     integer(kind=int1), dimension(:,:), allocatable:: Glint_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Glint_Mask_Lunar
     integer(kind=int1), dimension(:,:), allocatable:: Forward_Scatter_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Forward_Scatter_Mask_Lunar
     integer(kind=int1), dimension(:,:), allocatable:: Desert_Mask
     integer(kind=int1), dimension(:,:), allocatable:: City_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Volcano_Mask
     integer(kind=int1), dimension(:,:), allocatable:: Snow_OISST
     integer(kind=int1), dimension(:,:), allocatable:: Snow_NWP
     integer(kind=int1), dimension(:,:), allocatable:: Snow_IMS
     integer(kind=int1), dimension(:,:), allocatable:: Snow_GLOB
     integer(kind=int1), dimension(:,:), allocatable:: Snow
     integer(kind=int1), dimension(:,:), allocatable:: Sfc_Type
     real (kind=real4), dimension(:,:), allocatable:: Zsfc
     real (kind=real4), dimension(:,:), allocatable:: Zsfc_Hires
  end type surface_definition

  type :: sensor_definition
    character(len=32):: Sensor_Name
    integer(kind=int4):: Spatial_Resolution_Meters
    character(len=32):: Platform_Name
    integer(kind=int4):: WMO_Id
    integer(kind=int4):: WMO_Id_Previous
    character(len=1020):: Instr_Const_File
    character(len=1020):: Algo_Const_File
    real(kind=real8):: Geo_Sub_Satellite_Longitude
    real(kind=real8):: Geo_Sub_Satellite_Latitude
    logical, dimension(NCHAN_CLAVRX):: Chan_On_Flag_Default
    logical, dimension(:,:), allocatable:: Chan_On_Flag_Per_Line
  end type sensor_definition

  type :: image_definition
    character(len=1020):: Level1b_Name
    character(len=1020):: Level1b_Full_Name
    character(len=1020):: Level1b_Path
    integer(kind=int4):: Number_Of_Elements
    integer(kind=int4):: Number_Of_Lines
    integer(kind=int4):: Number_Of_Lines_Per_Segment
    integer(kind=int4):: Number_Of_Lines_Read_This_Segment
    integer(kind=int4):: Number_Of_Segments
    integer(kind=int4):: Segment_Number
    integer(kind=int2):: Start_Year
    integer(kind=int2):: Start_Doy
    integer(kind=int4):: Start_Time
    integer(kind=int2):: End_Year
    integer(kind=int2):: End_Doy
    integer(kind=int4):: End_Time
    character(len=1020) :: Auxiliary_Cloud_Mask_File_Name
    character(len=1020) :: Auxiliary_Geolocation_File_Name
  end type image_definition

  type :: acha_definition
    integer:: Mode
    real (kind=real4), dimension(:,:), allocatable:: Tc
    real (kind=real4), dimension(:,:), allocatable:: Ec
    real (kind=real4), dimension(:,:), allocatable:: Pc
    real (kind=real4), dimension(:,:), allocatable:: Zc
    real (kind=real4), dimension(:,:), allocatable:: Zc_Top
    real (kind=real4), dimension(:,:), allocatable:: Pc_Top
    real (kind=real4), dimension(:,:), allocatable:: Zc_Base
    real (kind=real4), dimension(:,:), allocatable:: Pc_Base
    real (kind=real4), dimension(:,:), allocatable:: Beta
    real (kind=real4), dimension(:,:), allocatable:: Tau
    real (kind=real4), dimension(:,:), allocatable:: Reff
    real (kind=real4), dimension(:,:), allocatable:: Tc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Ec_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Beta_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Zc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Pc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Lower_Tc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Lower_Zc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Lower_Pc_Uncertainty
    real (kind=real4), dimension(:,:), allocatable:: Alt
    real (kind=real4), dimension(:,:), allocatable:: Base_Alt
    real (kind=real4), dimension(:,:), allocatable:: Cost
    real (kind=real4), dimension(:,:), allocatable:: Lower_Pc
    real (kind=real4), dimension(:,:), allocatable:: Lower_Zc
    real (kind=real4), dimension(:,:), allocatable:: Lower_Tc
    integer(kind=int1), dimension(:,:), allocatable:: Processing_Order
    integer(kind=int1), dimension(:,:), allocatable:: Inversion_Flag
    integer (kind=int1), dimension(:,:), allocatable:: Quality_Flag
    integer (kind=int1), dimension(:,:), allocatable:: Meta_Data
    integer (kind=int1), dimension(:,:,:), allocatable:: OE_Quality_Flags
    integer (kind=int1), dimension(:,:), allocatable:: Packed_Quality_Flags
    integer (kind=int1), dimension(:,:), allocatable:: Packed_Meta_Data_Flags
    integer (kind=int1), dimension(:,:), allocatable:: base_Quality_Flag
    real (kind=real4), dimension(:,:), allocatable:: Conv_Cld_Prob
    real (kind=real4), dimension(:,:), allocatable:: Supercooled_Cld_Prob
    real(kind=real4):: Success_Fraction
    real(kind=real4):: Processed_Count
    real(kind=real4):: Valid_Count
    real (kind=real4), dimension(:,:), allocatable:: Ec_67um
    real (kind=real4), dimension(:,:), allocatable:: Ec_85um
    real (kind=real4), dimension(:,:), allocatable:: Ec_11um
    real (kind=real4), dimension(:,:), allocatable:: Ec_12um
    real (kind=real4), dimension(:,:), allocatable:: Ec_133um
  end type acha_definition

  type :: ccl_definition
    integer (kind=int1), dimension(:,:), allocatable:: Cld_Layer
    real(kind=real4), dimension(:,:), allocatable, public :: Cloud_Fraction
    real(kind=real4), dimension(:,:), allocatable, public :: Cloud_Fraction_Uncer
    real(kind=real4), dimension(:,:), allocatable, public :: High_Cloud_Fraction
    real(kind=real4), dimension(:,:), allocatable, public :: Mid_Cloud_Fraction
    real(kind=real4), dimension(:,:), allocatable, public :: Low_Cloud_Fraction
  end type ccl_definition

  type :: asos_definition
    integer (kind=int1), dimension(:,:), allocatable, public:: Code
    real(kind=real4), dimension(:,:), allocatable, public:: ECA
    real(kind=real4), dimension(:,:), allocatable, public:: Zmin
    real(kind=real4), dimension(:,:), allocatable, public:: Zmax
  end type asos_definition


  !---- declare structures using above types
  type(observations), dimension(Nchan_Clavrx), public, save, target :: Ch
  type(sensor_definition), public, save, target :: Sensor
  type(image_definition), public, save, target :: Image
  type(geometry_definition), public, save, target :: Geo
  type(navigation_definition), public, save, target :: Nav
  type(surface_definition), public, save, target :: Sfc
  type(acha_definition), public, save, target :: ACHA
  type(ccl_definition), public, save, target :: CCL
  type(asos_definition), public, save, target :: ASOS

  !---- declare other global variables
  integer,public, save:: Cmr_File_Flag
  integer,public, save:: Cloud_Mask_Aux_Flag
  integer,public, save:: Cloud_Mask_Aux_Read_Flag
  integer,public, save:: Cloud_Mask_Bayesian_Flag
  integer,public, save:: Ref_cal_1b 
  integer,public, save:: Therm_cal_1b
  integer,public, save:: Nav_Opt       !0=level1b,1=clevernav,2=reposnx
  integer,public, save:: Obs_File_Flag
  integer,public, save:: Geo_File_Flag
  integer,public, save:: Sst_File_Flag
  integer,public, save:: Cld_File_Flag
  integer,public, save:: Rtm_File_Flag
  integer,public, save:: Ash_File_Flag
  integer,public, save:: Level2_File_Flag
  integer,public, save:: Use_Sst_Anal
  logical,public, save:: Use_IR_Cloud_Type_Flag
  integer,public, save:: L1b_Gzip
  integer,public, save:: L1b_Bzip2
  
  integer,public, save:: Use_Seebor
  integer,public, save:: Read_Volcano_Mask
  integer,public, save:: Read_Land_Mask
  integer,public, save:: Read_Coast_Mask
  integer,public, save:: Read_Surface_Elevation
  integer,public, save:: Read_Hires_Sfc_Type
  integer,public, save:: Read_Snow_Mask
  integer,public, save:: Read_GLOBSnow_Mask
  integer,public, save:: Read_Dark_Comp
  integer,public, save:: Machine_Byte_Ordering
  integer,public, save:: LRC_Flag  !local radiative center flag
  integer,public, save:: Process_Undetected_Cloud_Flag
  integer,public, save:: DCOMP_Mode
  integer,public, save:: NLCOMP_Mode
  integer,public, save:: Mask_Mode
  integer,public, save:: Cld_Flag
  integer,public, save:: Blank_Flag
  integer,public, save:: Aer_Flag
  integer,public, save:: Ash_Flag
  integer,public, save:: Sasrab_Flag
  integer,public, save:: Nwp_Opt
  integer,public, save:: Modis_Clr_Alb_Flag

  integer,public, save:: Rtm_Opt
  integer,public, save:: Smooth_Nwp_Flag
  integer,public, save:: Compress_Flag
  !---------------------------------------------------------------------------------
  ! Flags Computed within CLAVR-x that describe the sensor data
  !---------------------------------------------------------------------------------
  integer,public, save:: AVHRR_GAC_Flag
  integer,public, save:: AVHRR_KLM_Flag
  integer,public, save:: AVHRR_AAPP_Flag
  integer,public, save:: AVHRR_1_Flag
  integer,public, save:: AVHRR_IFF_Flag
  integer,public, save:: Goes_Scan_Line_Flag
  

  !---------------------------------------------------------------------------------
  ! Internal Flags to communicate ancillary data information
  !---------------------------------------------------------------------------------
  integer,public, save:: Failed_IMS_Snow_Mask_Flag
  integer,public, save:: Failed_GLOB_Snow_Mask_Flag
  integer,public, save:: Output_Scaled_Reflectances
  integer,public, save:: Ncdc_Level2_Flag


  !---------------------------------------------------------------------------------
  ! variables that are computed to serve as attributes in the output files
  !---------------------------------------------------------------------------------
  real(kind=real4), public, save:: Orbital_Processing_Time_Minutes
  real(kind=real4), public, save:: DCOMP_Success_Fraction
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
  character(len=1020),public,save:: Ancil_Data_Dir
  character(len=1020),public,save:: CSV_File
  character(len=1020),public,save:: Gfs_Data_Dir
  character(len=1020),public,save:: Ncep_Data_Dir
  character(len=1020),public,save:: Cfsr_Data_Dir
  character(len=1020),public,save:: Merra_Data_Dir
  character(len=1020),public,save:: Gdas_Data_Dir
  character(len=1020),public,save:: Erai_Data_Dir
  character(len=1020),public,save:: Oisst_Data_Dir
  character(len=1020),public,save:: Snow_Data_Dir
  character(len=1020),public,save:: GLOBSnow_Data_Dir
  character(len=1020),public,save:: Dark_Comp_Data_Dir
  character(len=1020),public,save:: Temporary_Data_Dir
  character(len=1020),public,save:: File_Nav

  character(len=1020),public,save:: Dir_Rtm
  character(len=1020),public,save:: Dir_Level2
  character(len=1020),public,save:: Bayesian_Cloud_Mask_Name
  character(len=1020),public,save:: Dark_Composite_Name

  !----- IFF data files
  character(len=1020),public,save:: IFF_File

  !------------------------------------------------------------------
  !--- variables pertaining to scanline size
  !------------------------------------------------------------------
  integer, public, save:: Line_Idx_Min_Segment
  integer, public, save:: Line_Idx_Max_Segment 

  integer(kind=int4), public, save:: L1b_Rec_Length
  integer(kind=int4), public, save:: Num_Anchors
  integer , public, save :: Goes_Stride
 
  real(kind=real4), public, save:: dLat_hist2d

  !------- pixel array declarations
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
  integer(kind=int4), public, save:: Num_Scans_Level2_Hdf
  integer(kind=int2), public, save:: Sc_Id_Avhrr,AVHRR_Ver_1b, &
                                     AVHRR_Data_Type, Num_Loc, Tip_Parity, Aux_Sync,  &
                                     Ramp_Auto_Cal,Start_Year_Prev, Start_Day_Prev, &
                                     Month,Month_Prev,Day_of_Month,Ileap
  character(len=6), public, save:: Sc_Id_Char
  character(len=7),public,save:: Proc_Block_Id

  !--- instrument counts
  integer (kind=int2), dimension(:,:), allocatable, public,save:: Ch1_Counts
  integer (kind=int2), dimension(:,:), allocatable, public,save:: Ch2_Counts
  integer (kind=int2), dimension(:,:), allocatable, public,save:: Ch6_Counts
  real (kind=real4), dimension(:,:), allocatable, public,save:: Ch20_Counts_Filtered

  !--- sounder brightness temperatures  
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_375um_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_11um_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_12um_Sounder

  !--- MJH HIRS/AVHRR aux fields
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Cld_Temp_Sounder 
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Cld_Press_Sounder
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Cld_Height_Sounder 
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Cld_Emiss_Sounder 
   
  !--- calibrated observations
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChI2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChI3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_ChI4
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_ChI5

  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Mean_3x3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Max_3x3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Min_3x3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_ChDNB_Lunar_Std_3x3

  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Max_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Min_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Mean_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Uni_ChI1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Max_ChI2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Min_ChI2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Uni_ChI2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Mean_ChI2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Max_ChI3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Min_ChI3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Mean_ChI3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Ref_Uni_ChI3
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Max_ChI4
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Min_ChI4
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Uni_ChI4
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Mean_ChI4
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Max_ChI5
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Min_ChI5
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Uni_ChI5
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Mean_ChI5

! real (kind=real4), dimension(:,:), allocatable, public, save:: Rad_Ch20_Ems
  real (kind=real4), dimension(:,:), allocatable, public, save:: Ems_Ch20

  real (kind=real4), dimension(:,:), allocatable, public, save:: Sst_Anal
  real (kind=real4), dimension(:,:), allocatable, public, save:: Sst_Anal_Err
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Sst_Anal_Uni
  real (kind=real4), dimension(:,:), allocatable, public, save:: Sst_Anal_Cice
  real (kind=real4), dimension(:,:), allocatable, public, save:: Tsfc_Retrieved
  real (kind=real4), dimension(:,:), allocatable, public, save:: Trad_Retrieved

  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Temp_Pix_Array_1
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Temp_Pix_Array_2
  real (kind=real4), dimension(:,:), allocatable, public, save, target:: Temp_Pix_Array_3
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
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch20_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Bt_Ch27_Max_3x3

  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch27_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch29_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch32_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch33_Std_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Ems_Ch20_Median_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ems_Ch20_Std_Median_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Bt_Ch20_Median_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Bt_Ch20_Std_Median_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save:: Ref_Ch1_Sfc_White_Sky_Mean_3x3
  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Btd_Ch31_Ch32_Bt_Ch31_Max_3x3

  real(kind=real4), dimension(:,:), allocatable, public, save, target:: Covar_Ch27_Ch31_5x5

  integer(kind=int4), dimension(:,:), allocatable, public, save:: Elem_Idx_Max_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Line_Idx_Max_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Elem_Idx_Min_Bt_Ch31_3x3
  integer(kind=int4), dimension(:,:), allocatable, public, save:: Line_Idx_Min_Bt_Ch31_3x3

  real (kind=real4), dimension(:,:), allocatable, public:: Sst_Unmasked   !sst used in cld Mask
  real (kind=real4), dimension(:,:), allocatable, public:: Sst_Masked !sst where non-clear ocean is Masked
  real (kind=real4), dimension(:,:), allocatable, public:: Ndvi_Toa
  real (kind=real4), dimension(:,:), allocatable, public:: Ndsi_Toa
  real (kind=real4), dimension(:,:), allocatable, public:: Ndsi_Sfc
  real (kind=real4), dimension(:,:), allocatable, public, target:: Btd_Ch31_Ch32
  real (kind=real4), dimension(:,:), allocatable, public:: Btd_Ch20_Ch31
  real (kind=real4), dimension(:,:), allocatable, public:: Btd_Ch20_Ch32

  integer(kind=int1), dimension(:,:), allocatable, public, target:: Solar_Contamination_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Bad_Pixel_Mask
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
  integer(kind=int1), dimension(:,:), allocatable, public:: Bayes_Mask_Sfc_Type_Global
  integer(kind=int1), dimension(:,:), allocatable, public:: Shadow_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Dust_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Smoke_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Fire_Mask
  integer(kind=int1), dimension(:,:), allocatable, public, target:: Thin_Cirr_Mask

  !--- cloud Mask arrays
  integer (kind=int1), dimension(:,:,:), allocatable, public, save, target:: Cld_Test_Vector_Packed
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Mask
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Adj_Pix_Cld_Mask
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Cld_Mask_Qf
  real (kind=real4),dimension(:,:),allocatable, public, save, target:: &
                                                       Posterior_Cld_Probability
  real (kind=real4),dimension(:,:),allocatable, public, save, target:: &
                                                       Prior_Cld_Probability

  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Type
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Phase
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Type_IR
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Cld_Phase_IR
  integer (kind=int1),dimension(:,:),allocatable, public, save, target:: Ctp_Multilayer_Flag

  !--- Auxilliary variables
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Cld_Mask_Aux
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Cld_Type_Aux
  integer (kind=int1),dimension(:,:),allocatable, public, save:: Cld_Phase_Aux
  real (kind=real4),dimension(:,:),allocatable, public, save, target::Zc_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save:: Pc_Top1_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save:: Pc_Top2_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save:: Pc_Uncertainty1_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save:: Pc_Uncertainty2_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save:: Cost_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save:: Tau_Aux
  real (kind=real4), dimension(:,:), allocatable, public, save:: Reff_Aux

  !--- pixel level cloud props

     !--- h2o heights
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tc_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_H2O
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zclr_H2O_Peak
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Pc_CO2IRW
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Zc_CO2IRW
     real (kind=real4), dimension(:,:), allocatable, public, save, target:: Tc_CO2IRW

     !-- DCOMP cloud algorithm results
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_1
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_2
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_3
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_ap
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: vis_Ref_fm
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_1
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_2
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_3
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Iwp_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Iwp_Tau_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Lwp_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Ice_Layer_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Water_Layer_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cwp_Scwater_Layer_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Rain_Rate_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Hcld_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Cdnc_DCOMP
     real (kind=real4), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_Cost
     real (kind=real4), dimension(:,:), allocatable, public, target, save:: Reff_DCOMP_Cost
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: Tau_DCOMP_Qf
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: Reff_DCOMP_Qf
     integer (kind=int1), dimension(:,:), allocatable, public,target, save:: DCOMP_Quality_Flag
     integer (kind=int2), dimension(:,:), allocatable, public,target, save:: DCOMP_Info_Flag
     

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
     real (kind=real4), dimension(:,:), allocatable, public:: Ndvi_Sfc

  !--- scratch arrays
  integer(kind=int1), dimension(:,:), public,save,allocatable, target:: One_Byte_Temp
  integer(kind=int2), dimension(:,:), public,save,allocatable, target:: Two_Byte_Temp
  integer(kind=int1),dimension(:,:),allocatable, public, save, target:: Temp_Mask

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
  real (kind=real4), dimension(:,:), allocatable, public:: Weasd_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Sea_Ice_Frac_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Tpw_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Ozone_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: K_Index_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Sc_Lwp_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Lwp_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Iwp_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public, target:: Cwp_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Pc_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public, target:: LCL_Height_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public, target:: CCL_Height_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Cfrac_Nwp_Pix
  integer (kind=int1), dimension(:,:), allocatable, public:: Ncld_Layers_Nwp_Pix
  integer (kind=int1), dimension(:,:), allocatable, public:: Cld_Type_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Wnd_Spd_10m_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Wnd_Dir_10m_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Wnd_Spd_Cld_Top_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Wnd_Dir_Cld_Top_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Inversion_Base_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Inversion_Top_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Inversion_Strength_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public:: Tpw_Above_Cloud_Nwp_Pix

  real (kind=real4), dimension(:,:), allocatable, public:: Trans_Atm_Ch20_Solar_Total_Rtm
  real (kind=real4), dimension(:,:), allocatable, public:: Trans_Atm_Ch20_Solar_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ems_Ch20_Clear_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ttropo_Nwp_Pix
  real (kind=real4), dimension(:,:), allocatable, public, target:: Emiss_11um_Tropo_Nadir_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_12um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_85um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_67um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_133um_Tropo_Rtm
  real (kind=real4), dimension(:,:), allocatable, public, target:: Beta_11um_133fusum_Tropo_Rtm


  real (kind=real4), dimension(:,:), allocatable, public, target:: Pc_Opaque_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Zc_Opaque_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Tc_Opaque_Cloud
  real (kind=real4), dimension(:,:), allocatable, public, target:: Tc_Co2
  real (kind=real4), dimension(:,:), allocatable, public, target:: Pc_Co2
  real (kind=real4), dimension(:,:), allocatable, public, target:: Zc_Co2
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ec_Co2
  real (kind=real4), dimension(:,:), allocatable, public, target:: Tc_Cirrus_Background 
  real (kind=real4), dimension(:,:), allocatable, public, target:: Zc_Cirrus_Background 
!--- modis white sky albedo maps
  real (kind=real4), dimension(:,:), allocatable, public, target:: Ndvi_Sfc_White_Sky

!---- other static arrays carried by this module


  !--- Solar RTM Terms
  type, public :: solar_rtm_struct
      real, dimension(Nchan_Clavrx,3):: Tau_H2O_Coef
      real, dimension(Nchan_Clavrx):: Tau_Ray
      real, dimension(Nchan_Clavrx):: Tau_O2
      real, dimension(Nchan_Clavrx):: Tau_O3
      real, dimension(Nchan_Clavrx):: Tau_CH4
      real, dimension(Nchan_Clavrx):: Tau_CO2
      real, dimension(Nchan_Clavrx):: Tau_Aer
      real, dimension(Nchan_Clavrx):: Wo_Aer
      real, dimension(Nchan_Clavrx):: G_Aer
  end type solar_rtm_struct
 
  type (solar_rtm_struct), public, save:: Solar_Rtm

  !--- flags for using clavrxorb_Default_file
  integer ,public, save :: Use_Default

  !--- clavrxorb_File_List filename
  character(len=1020), public, save:: File_List

 contains

!----------------------------------------------------------------------------
! This routine allocate the memory for the pixel arrays 
!----------------------------------------------------------------------------
subroutine CREATE_PIXEL_ARRAYS()
  integer:: dim1 
  integer:: dim2 
  integer:: idx

  !--- allocate pixel arrays
  dim1 = Image%Number_Of_Elements
  dim2 = Image%Number_Of_Lines_Per_Segment

  !---- new
!----------------------------------------------------------------------
!  begin allocation of ch structure
!----------------------------------------------------------------------

  !--- set obs type for each channel
  Ch(1:19)%Obs_Type = SOLAR_OBS_TYPE
  Ch(20)%Obs_Type = MIXED_OBS_TYPE
  Ch(21)%Obs_Type = MIXED_OBS_TYPE
  Ch(22:25)%Obs_Type = THERMAL_OBS_TYPE
  Ch(26)%Obs_Type = SOLAR_OBS_TYPE
  Ch(27:38)%Obs_Type = THERMAL_OBS_TYPE
  Ch(39:41)%Obs_Type = SOLAR_OBS_TYPE
  Ch(42)%Obs_Type = MIXED_OBS_TYPE
  Ch(43)%Obs_Type = THERMAL_OBS_TYPE
  Ch(44)%Obs_Type = LUNAR_OBS_TYPE
  Ch(45)%Obs_Type = THERMAL_OBS_TYPE

  !--- loop through each that is on, allocate fields based on obs type 
  do idx = 1,Nchan_Clavrx
      if (Sensor%Chan_On_Flag_Default(idx) == sym%YES) then

        allocate(Ch(idx)%Unc(dim1,dim2))

        select case (ch(idx)%Obs_Type)

        case(SOLAR_OBS_TYPE)
            allocate(Ch(idx)%Ref_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa_Unnorm(dim1,dim2))
            allocate(Ch(idx)%Ref_Sfc(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm_Total(dim1,dim2))
            if (idx == 1) allocate(Ch(idx)%Opd(dim1,dim2))

        case(LUNAR_OBS_TYPE)
            allocate(Ch(idx)%Rad_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Lunar_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Lunar_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Ref_Lunar_Sfc(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm_Total(dim1,dim2))
            allocate(Ch(idx)%Opd(dim1,dim2))

        case(THERMAL_OBS_TYPE)
            allocate(Ch(idx)%Rad_Toa(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa(dim1,dim2))
            allocate(Ch(idx)%Rad_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Rad_Atm(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm(dim1,dim2))
            allocate(Ch(idx)%Sfc_Emiss(dim1,dim2))
            if (idx == 31) allocate(Ch(idx)%Rad_Atm_Dwn_Sfc(dim1,dim2))
            if (idx == 27 .or. idx == 29 .or. idx == 31 .or. idx == 32 .or.  &
                idx == 33 .or. idx == 45) allocate(Ch(idx)%Emiss_Tropo(dim1,dim2))
            if (idx >= 27 .and. idx <= 38)  allocate(Ch(idx)%CSBT_Mask(dim1,dim2))
            if (idx >= 27 .and. idx <= 38)  allocate(Ch(idx)%Opaque_Height(dim1,dim2))

        case(MIXED_OBS_TYPE)
            allocate(Ch(idx)%Ref_Toa(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa_Unnorm(dim1,dim2))
            allocate(Ch(idx)%Ref_Sfc(dim1,dim2))
            allocate(Ch(idx)%Ref_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm_Total(dim1,dim2))
            allocate(Ch(idx)%Rad_Toa(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa(dim1,dim2))
            allocate(Ch(idx)%Rad_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Bt_Toa_Clear(dim1,dim2))
            allocate(Ch(idx)%Rad_Atm(dim1,dim2))
            allocate(Ch(idx)%Trans_Atm(dim1,dim2))
            allocate(Ch(idx)%Sfc_Emiss(dim1,dim2))

        case default

           call MESG ("Arrays not allocated for unknown channel = ",idx) 

        end select

      endif

  enddo

  !--- force allocation of channel 20 surface emissivity for use in mask
  if (.not. allocated(Ch(20)%Sfc_Emiss)) allocate(Ch(20)%Sfc_Emiss(dim1,dim2))
  
!----------------------------------------------------------------------
!  end allocation of ch structure
!----------------------------------------------------------------------

   if ((Sensor%Chan_On_Flag_Default(27) == sym%YES) .and.   &
       (Sensor%Chan_On_Flag_Default(31) == sym%YES)) then
           allocate(Covar_Ch27_Ch31_5x5(dim1,dim2))
   endif

   allocate(Bad_Pixel_Mask(dim1,dim2), &
            Space_Mask(dim1,dim2))

   allocate(Sfc_Level_Rtm_Pixel(dim1,dim2))
   allocate(Solar_Contamination_Mask(dim1,dim2))
    
   !--- VIIRS Arrays
   allocate(IFF_Gap_Mask(dim1,dim2))
   allocate(Gap_Pixel_Mask(dim1,dim2))
   allocate(Gap_Line_Idx(dim1,dim2))

   !--- 3x3 uni
   if (Sensor%Chan_On_Flag_Default(27) == sym%YES) then
          allocate (Bt_Ch27_Max_3x3(dim1,dim2))
   endif

   allocate(   &
          Dust_Mask(dim1,dim2), &
          Smoke_Mask(dim1,dim2), &
          Fire_Mask(dim1,dim2), &
          Thin_Cirr_Mask(dim1,dim2), &
          Shadow_Mask(dim1,dim2), &
          Sst_Anal(dim1,dim2), &
          Sst_Anal_Err(dim1,dim2), &
          Sst_Anal_Cice(dim1,dim2), &
          Sst_Anal_Uni(dim1,dim2))

  allocate(Beta_11um_12um_Tropo_Rtm(dim1,dim2))
  allocate(Beta_11um_85um_Tropo_Rtm(dim1,dim2))
  allocate(Beta_11um_67um_Tropo_Rtm(dim1,dim2))
  allocate(Beta_11um_133um_Tropo_Rtm(dim1,dim2))
  allocate(Beta_11um_133fusum_Tropo_Rtm(dim1,dim2))

  allocate(Temp_Mask(dim1,dim2))

  !--- sensor arrays
  call  CREATE_SENSOR_ARRAYS(Nchan_Clavrx,dim2)
  !--- navigation arrays
  call  CREATE_NAV_ARRAYS(dim1, dim2)
  !--- geometry arrays
  call  CREATE_GEO_ARRAYS(dim1, dim2)
  !--- anchor point arrays
  call  CREATE_AVHRR_ANCHOR_ARRAYS(Num_Anchors, dim2)
  !--- nwp fields interpolated to the pixel level
  call  CREATE_NWP_PIX_ARRAYS(dim1, dim2)
  call  CREATE_SURFACE_ARRAYS(dim1, dim2)
  !--- ch1 arrays
  call  CREATE_REF_CHANNEL_ARRAYS(dim1, dim2)
  call  CREATE_THERM_CHANNEL_ARRAYS(dim1, dim2)
  call  CREATE_EXTRA_CHANNEL_ARRAYS(dim1, dim2)
! call  CREATE_LUNAR_ARRAYS(dim1, dim2)
  call  CREATE_BTD_ARRAYS(dim1, dim2)
  call  CREATE_ACHA_ARRAYS(dim1, dim2)
  call  CREATE_CCL_ARRAYS(dim1, dim2)
  call  CREATE_ASOS_ARRAYS(dim1, dim2)
  call  CREATE_DCOMP_ARRAYS(dim1, dim2)
  call  CREATE_NLCOMP_ARRAYS(dim1, dim2)
  call  CREATE_SASRAB_ARRAYS(dim1, dim2)
  call  CREATE_OLR_ARRAYS(dim1, dim2)
  call  CREATE_AEROSOL_ARRAYS(dim1, dim2)
  call  CREATE_CLOUD_MASK_ARRAYS(dim1, dim2, Max_Num_Cld_Test_Bytes)
  call  CREATE_CLOUD_TYPE_ARRAYS(dim1, dim2)
  call  CREATE_DIAGNOSTIC_ARRAYS(dim1, dim2)
  call  CREATE_SFC_PROD_ARRAYS(dim1, dim2)
  call  CREATE_CLOUD_PROD_ARRAYS(dim1, dim2)

  !--- pixel level parameters
   allocate(Zen_Idx_Rtm(dim1,dim2))
   allocate(i_LRC(dim1,dim2))
   allocate(j_LRC(dim1,dim2))

   
  allocate(Ch3a_On_AVHRR(dim2), &
           Bad_Scan_Flag(dim2), &
           Scan_Number(dim2), &
           Scan_Time(dim2), &
           Scan_Day(dim2), &
           Scan_Year(dim2), &
           Utc_Scan_Time_Hours(dim2), &
           Pixel_Local_Time_Hours(dim1,dim2), &
           Pixel_Time(dim1,dim2))

  !--------------------------------------------------------------------------------
  ! Initialize variables that are not reset for each segment
  !--------------------------------------------------------------------------------

  !--- metrics - needed initialize counts to be zero but they accumulate through orbit
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

  do idx = 1,Nchan_Clavrx
      if (allocated(Ch(idx)%Rad_Toa)) deallocate(Ch(idx)%Rad_Toa)
      if (allocated(Ch(idx)%Bt_Toa)) deallocate(Ch(idx)%Bt_Toa)
      if (allocated(Ch(idx)%Rad_Toa_Clear)) deallocate(Ch(idx)%Rad_Toa_Clear)
      if (allocated(Ch(idx)%Bt_Toa_Clear)) deallocate(Ch(idx)%Bt_Toa_Clear)
      if (allocated(Ch(idx)%Rad_Atm)) deallocate(Ch(idx)%Rad_Atm)
      if (allocated(Ch(idx)%Rad_Atm_Dwn_Sfc)) deallocate(Ch(idx)%Rad_Atm_Dwn_Sfc)
      if (allocated(Ch(idx)%Trans_Atm)) deallocate(Ch(idx)%Trans_Atm)
      if (allocated(Ch(idx)%Trans_Atm_Total)) deallocate(Ch(idx)%Trans_Atm_Total)
      if (allocated(Ch(idx)%Ref_Toa)) deallocate(Ch(idx)%Ref_Toa)
      if (allocated(Ch(idx)%Ref_Toa_Unnorm)) deallocate(Ch(idx)%Ref_Toa_Unnorm)
      if (allocated(Ch(idx)%Ref_Toa_Clear)) deallocate(Ch(idx)%Ref_Toa_Clear)
      if (allocated(Ch(idx)%Ref_Sfc)) deallocate(Ch(idx)%Ref_Sfc)
      if (allocated(Ch(idx)%Ref_Lunar_Toa)) deallocate(Ch(idx)%Ref_Lunar_Toa)
      if (allocated(Ch(idx)%Ref_Lunar_Toa_Clear)) deallocate(Ch(idx)%Ref_Lunar_Toa_Clear)
      if (allocated(Ch(idx)%Ref_Lunar_Sfc)) deallocate(Ch(idx)%Ref_Lunar_Sfc)
      if (allocated(Ch(idx)%Opd)) deallocate(Ch(idx)%Opd)
      if (allocated(Ch(idx)%Sfc_Emiss)) deallocate(Ch(idx)%Sfc_Emiss)
      if (allocated(Ch(idx)%Sfc_Ref_White_Sky)) deallocate(Ch(idx)%Sfc_Ref_White_Sky)
      if (allocated(Ch(idx)%Emiss_Tropo)) deallocate(Ch(idx)%Emiss_Tropo)
      if (allocated(Ch(idx)%Unc)) deallocate(Ch(idx)%Unc)
      if (allocated(Ch(idx)%CSBT_Mask)) deallocate(Ch(idx)%CSBT_Mask)
      if (allocated(Ch(idx)%Opaque_Height)) deallocate(Ch(idx)%Opaque_Height)
  enddo

  if (allocated(Temp_Mask)) deallocate(Temp_Mask)

  deallocate(Ch3a_On_AVHRR, &
             Bad_Scan_Flag,  &
             Scan_Number, &
             Scan_Time,Scan_Day,Scan_Year,Utc_Scan_Time_Hours, &
             Pixel_Local_Time_Hours,Pixel_Time)

  if (allocated(Covar_Ch27_Ch31_5x5)) deallocate(Covar_Ch27_Ch31_5x5)

  deallocate(Bad_Pixel_Mask)

  call DESTROY_SENSOR_ARRAYS()
  call DESTROY_NAV_ARRAYS()
  call DESTROY_GEO_ARRAYS()
  call DESTROY_AVHRR_ANCHOR_ARRAYS()
  call DESTROY_NWP_PIX_ARRAYS()
  call DESTROY_SURFACE_ARRAYS()
  call DESTROY_REF_CHANNEL_ARRAYS()
  call DESTROY_THERM_CHANNEL_ARRAYS()
  call DESTROY_EXTRA_CHANNEL_ARRAYS()
! call DESTROY_LUNAR_ARRAYS()
  call DESTROY_BTD_ARRAYS()
  call DESTROY_ACHA_ARRAYS()
  call DESTROY_CCL_ARRAYS()
  call DESTROY_ASOS_ARRAYS()
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

  deallocate(Space_Mask)
  deallocate(Sfc_Level_Rtm_Pixel)
  deallocate(Fire_Mask)
  deallocate(Thin_Cirr_Mask)
  deallocate(Dust_Mask)
  deallocate(Smoke_Mask)
  deallocate(Shadow_Mask)

  deallocate(Sst_Anal)
  deallocate(Sst_Anal_Err)
  deallocate(Sst_Anal_Cice)
  deallocate(Sst_Anal_Uni)

  deallocate(Solar_Contamination_Mask)
 

  !--- nwp and rtm indices
  if (allocated(Zen_Idx_Rtm)) deallocate(Zen_Idx_Rtm)

  !--- local radiative center indices
  if (allocated(i_LRC)) deallocate(i_LRC)
  if (allocated(j_LRC)) deallocate(j_LRC)

!--- pixel rtm parameters
  if (Sensor%Chan_On_Flag_Default(27) == sym%YES) then
    if (allocated(Bt_Ch27_Max_3x3)) deallocate(Bt_Ch27_Max_3x3)
  endif

  !--- ir cloud layer
  if (allocated(Beta_11um_12um_Tropo_Rtm)) deallocate(Beta_11um_12um_Tropo_Rtm)
  if (allocated(Beta_11um_67um_Tropo_Rtm)) deallocate(Beta_11um_67um_Tropo_Rtm)
  if (allocated(Beta_11um_85um_Tropo_Rtm)) deallocate(Beta_11um_85um_Tropo_Rtm)
  if (allocated(Beta_11um_133um_Tropo_Rtm)) deallocate(Beta_11um_133um_Tropo_Rtm)
  if (allocated(Beta_11um_133fusum_Tropo_Rtm)) deallocate(Beta_11um_133fusum_Tropo_Rtm)

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
      Bad_Pixel_Mask = sym%NO      !not initialized to missing
      Ch3a_On_AVHRR = Missing_Value_Int1

      call RESET_SENSOR_ARRAYS()
      call RESET_NAV_ARRAYS()
      call RESET_GEO_ARRAYS()
      call RESET_AVHRR_ANCHOR_ARRAYS()
      call RESET_NWP_PIX_ARRAYS()
      call RESET_SURFACE_ARRAYS()
      call RESET_REF_CHANNEL_ARRAYS()
      call RESET_THERM_CHANNEL_ARRAYS()
      call RESET_EXTRA_CHANNEL_ARRAYS()
      call RESET_BTD_ARRAYS()
      call RESET_ACHA_ARRAYS()
      call RESET_CCL_ARRAYS()
      call RESET_ASOS_ARRAYS()
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

      if (Sensor%Chan_On_Flag_Default(27) == sym%YES) THEN
        Bt_Ch27_Max_3x3 = Missing_Value_Real4
      endif

      
      Sst_Anal = Missing_Value_Real4
      Sst_Anal_Err = Missing_Value_Real4
      Sst_Anal_Cice = Missing_Value_Real4
      Sst_Anal_Uni = Missing_Value_Real4

      Zen_Idx_Rtm = Missing_Value_Int1

      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(32) == sym%YES)) then
        Beta_11um_12um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(27) == sym%YES)) then
        Beta_11um_67um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(29) == sym%YES)) then
        Beta_11um_85um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(33) == sym%YES)) then
        Beta_11um_133um_Tropo_Rtm = Missing_Value_Real4
      endif
      if ((Sensor%Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Sensor%Chan_On_Flag_Default(45) == sym%YES)) then
        Beta_11um_133fusum_Tropo_Rtm = Missing_Value_Real4
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
subroutine CREATE_SENSOR_ARRAYS(Nchan,dim2)
   integer, intent(in):: Nchan, dim2
   allocate(Sensor%Chan_On_Flag_Per_Line(Nchan,dim2))
end subroutine CREATE_SENSOR_ARRAYS
subroutine RESET_SENSOR_ARRAYS()
   if (allocated(Sensor%Chan_On_Flag_Per_Line)) Sensor%Chan_On_Flag_Per_Line = Missing_Value_Int1
end subroutine RESET_SENSOR_ARRAYS
subroutine DESTROY_SENSOR_ARRAYS()
   if (allocated(Sensor%Chan_On_Flag_Per_Line)) deallocate(Sensor%Chan_On_Flag_Per_Line)
end subroutine DESTROY_SENSOR_ARRAYS
!-------------------------------------------------------------
!
!-------------------------------------------------------------
subroutine CREATE_NAV_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(Nav%Ascend(dim2))
   allocate(Nav%Lat(dim1,dim2))
   allocate(Nav%Lon(dim1,dim2))
   allocate(Nav%Lat_1b(dim1,dim2))
   allocate(Nav%Lon_1b(dim1,dim2))
   allocate(Nav%Lat_Pc(dim1,dim2))
   allocate(Nav%Lon_Pc(dim1,dim2))
   if (index(Sensor%Sensor_Name,'IFF') > 0) then
     allocate(Nav%Sounder_Fov(dim1,dim2))
     allocate(Nav%Sounder_Fov_Mask(dim1,dim2))
     allocate(Nav%Sounder_Fov_Segment_Idx(dim1,dim2))
     allocate(Nav%Sounder_X(dim1,dim2))
     allocate(Nav%Sounder_Y(dim1,dim2))
   endif
end subroutine CREATE_NAV_ARRAYS
subroutine DESTROY_NAV_ARRAYS()
  if (allocated(Nav%Ascend)) deallocate(Nav%Ascend)
  if (allocated(Nav%Sounder_Fov)) deallocate(Nav%Sounder_Fov)
  if (allocated(Nav%Sounder_Fov_Mask)) deallocate(Nav%Sounder_Fov_Mask)
  if (allocated(Nav%Sounder_Fov_Segment_Idx)) deallocate(Nav%Sounder_Fov_Segment_Idx)
  if (allocated(Nav%Sounder_X)) deallocate(Nav%Sounder_X)
  if (allocated(Nav%Sounder_Y)) deallocate(Nav%Sounder_Y)
  if (allocated(Nav%Lat)) deallocate(Nav%Lat)
  if (allocated(Nav%Lon)) deallocate(Nav%Lon)
  if (allocated(Nav%Lat_1b)) deallocate(Nav%Lat_1b)
  if (allocated(Nav%Lon_1b)) deallocate(Nav%Lon_1b)
  if (allocated(Nav%Lat_Pc)) deallocate(Nav%Lat_Pc)
  if (allocated(Nav%Lon_Pc)) deallocate(Nav%Lon_Pc)
end subroutine
subroutine RESET_NAV_ARRAYS()
  if (allocated(Nav%Ascend)) Nav%Ascend = Missing_Value_Int1
  if (allocated(Nav%Sounder_Fov)) Nav%Sounder_Fov = Missing_Value_Int1
  if (allocated(Nav%Sounder_Fov_Mask)) Nav%Sounder_Fov_Mask = Missing_Value_Int1
  if (allocated(Nav%Sounder_Fov_Segment_Idx)) Nav%Sounder_Fov_Segment_Idx = Missing_Value_Int2
  if (allocated(Nav%Sounder_X)) Nav%Sounder_X = Missing_Value_Int2
  if (allocated(Nav%Sounder_Y)) Nav%Sounder_Y = Missing_Value_Int2
  if (allocated(Nav%Lat)) Nav%Lat = Missing_Value_Real4
  if (allocated(Nav%Lon)) Nav%Lon = Missing_Value_Real4
  if (allocated(Nav%Lat_1b)) Nav%Lat_1b = Missing_Value_Real4
  if (allocated(Nav%Lon_1b)) Nav%Lon_1b = Missing_Value_Real4
  if (allocated(Nav%Lat_Pc)) Nav%Lat_Pc = Missing_Value_Real4
  if (allocated(Nav%Lon_Pc)) Nav%Lon_Pc = Missing_Value_Real4
end subroutine RESET_NAV_ARRAYS
!------------------------------------------------------------------------------
!  routines to create, destroy and reset geo structure
!------------------------------------------------------------------------------
subroutine CREATE_GEO_ARRAYS(dim1,dim2)
    integer, intent(in):: dim1, dim2
    allocate (Geo%Satzen(dim1,dim2))
    allocate (Geo%Solzen(dim1,dim2))
    allocate (Geo%Sataz(dim1,dim2))
    allocate (Geo%Solaz(dim1,dim2))
    allocate (Geo%Relaz(dim1,dim2))
    allocate (Geo%Glintzen(dim1,dim2))
    allocate (Geo%Seczen(dim1,dim2))
    allocate (Geo%Coszen(dim1,dim2))
    allocate (Geo%Cossolzen(dim1,dim2))
    allocate (Geo%Scatangle(dim1,dim2))
    allocate (Geo%Airmass(dim1,dim2))
    if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
           allocate(Geo%Lunzen(dim1,dim2))
           allocate(Geo%Lunaz(dim1,dim2))
           allocate(Geo%LunRelaz(dim1,dim2))
           allocate(Geo%Scatangle_Lunar(dim1,dim2))
           allocate(Geo%Glintzen_Lunar(dim1,dim2))
    endif
end subroutine CREATE_GEO_ARRAYS
subroutine DESTROY_GEO_ARRAYS()
  if (allocated(Geo%Satzen)) deallocate(Geo%Satzen)
  if (allocated(Geo%Solzen)) deallocate(Geo%Solzen)
  if (allocated(Geo%Sataz)) deallocate(Geo%Sataz)
  if (allocated(Geo%Solaz)) deallocate(Geo%Solaz)
  if (allocated(Geo%Relaz)) deallocate(Geo%Relaz)
  if (allocated(Geo%Glintzen)) deallocate(Geo%Glintzen)
  if (allocated(Geo%Seczen)) deallocate(Geo%Seczen)
  if (allocated(Geo%Coszen)) deallocate(Geo%Coszen)
  if (allocated(Geo%Cossolzen)) deallocate(Geo%Cossolzen)
  if (allocated(Geo%Scatangle)) deallocate(Geo%Scatangle)
  if (allocated(Geo%Airmass)) deallocate(Geo%Airmass)
  if (allocated(Geo%Lunzen))deallocate(Geo%Lunzen)
  if (allocated(Geo%Lunaz)) deallocate(Geo%Lunaz)
  if (allocated(Geo%LunRelaz)) deallocate(Geo%LunRelaz)
  if (allocated(Geo%Scatangle_Lunar)) deallocate(Geo%Scatangle_Lunar)
  if (allocated(Geo%Glintzen_Lunar)) deallocate(Geo%Glintzen_Lunar)
end subroutine DESTROY_GEO_ARRAYS
subroutine RESET_GEO_ARRAYS()
  if (allocated(Geo%Satzen)) Geo%Satzen = Missing_Value_Real4
  if (allocated(Geo%Solzen)) Geo%Solzen = Missing_Value_Real4
  if (allocated(Geo%Sataz)) Geo%Sataz = Missing_Value_Real4
  if (allocated(Geo%Solaz)) Geo%Solaz = Missing_Value_Real4
  if (allocated(Geo%Relaz)) Geo%Relaz = Missing_Value_Real4
  if (allocated(Geo%Glintzen)) Geo%Glintzen = Missing_Value_Real4
  if (allocated(Geo%Seczen)) Geo%Seczen = Missing_Value_Real4
  if (allocated(Geo%Coszen)) Geo%Coszen = Missing_Value_Real4
  if (allocated(Geo%Cossolzen)) Geo%Cossolzen = Missing_Value_Real4
  if (allocated(Geo%Scatangle)) Geo%Scatangle = Missing_Value_Real4
  if (allocated(Geo%Airmass)) Geo%Airmass = Missing_Value_Real4
  if (allocated(Geo%Lunzen)) Geo%Lunzen = Missing_Value_Real4
  if (allocated(Geo%Lunaz)) Geo%Lunaz = Missing_Value_Real4
  if (allocated(Geo%LunRelaz)) Geo%LunRelaz = Missing_Value_Real4
  if (allocated(Geo%Scatangle_Lunar)) Geo%Scatangle_Lunar = Missing_Value_Real4
  if (allocated(Geo%Glintzen_Lunar)) Geo%Glintzen_Lunar = Missing_Value_Real4
  Geo%Moon_Phase_Angle = Missing_Value_Real4
  Geo%Moon_Illum_Frac = Missing_Value_Real4
end subroutine RESET_GEO_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_AVHRR_ANCHOR_ARRAYS(dim1,dim2)
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
end subroutine CREATE_AVHRR_ANCHOR_ARRAYS
subroutine RESET_AVHRR_ANCHOR_ARRAYS()
  Lat_Anchor_1b = Missing_Value_Real4
  Lon_Anchor_1b = Missing_Value_Real4
  Solzen_Anchor = Missing_Value_Real4
  Satzen_Anchor = Missing_Value_Real4
  Scatangle_Anchor = Missing_Value_Real4
  Glintzen_Anchor = Missing_Value_Real4
  Relaz_Anchor = Missing_Value_Real4
  Solaz_Anchor = Missing_Value_Real4
  Sataz_Anchor = Missing_Value_Real4
end subroutine RESET_AVHRR_ANCHOR_ARRAYS
subroutine DESTROY_AVHRR_ANCHOR_ARRAYS
  deallocate(Lat_Anchor_1b)
  deallocate(Lon_Anchor_1b)
  deallocate(Solzen_Anchor)
  deallocate(Satzen_Anchor)
  deallocate(Scatangle_Anchor)
  deallocate(Glintzen_Anchor)
  deallocate(Relaz_Anchor)
  deallocate(Solaz_Anchor)
  deallocate(Sataz_Anchor)
end subroutine DESTROY_AVHRR_ANCHOR_ARRAYS
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
   allocate(Weasd_Nwp_Pix(dim1,dim2))
   allocate(Sea_Ice_Frac_Nwp_Pix(dim1,dim2))
   allocate(Tpw_Nwp_Pix(dim1,dim2))
   allocate(Ozone_Nwp_Pix(dim1,dim2))
   allocate(Ttropo_Nwp_Pix(dim1,dim2))
   allocate(Wnd_Spd_10m_Nwp_Pix(dim1,dim2))
   allocate(Wnd_Dir_10m_Nwp_Pix(dim1,dim2))
   allocate(Wnd_Spd_Cld_Top_Nwp_Pix(dim1,dim2))
   allocate(Wnd_Dir_Cld_Top_Nwp_Pix(dim1,dim2))
   allocate(Inversion_Base_Nwp_Pix(dim1,dim2))
   allocate(Inversion_Top_Nwp_Pix(dim1,dim2))
   allocate(Inversion_Strength_Nwp_Pix(dim1,dim2))
   allocate(Tpw_Above_Cloud_Nwp_Pix(dim1,dim2))
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
   Sea_Ice_Frac_Nwp_Pix = Missing_Value_Real4
   Weasd_Nwp_Pix = Missing_Value_Real4
   Tpw_Nwp_Pix = Missing_Value_Real4
   Ozone_Nwp_Pix = Missing_Value_Real4
   Ttropo_Nwp_Pix = Missing_Value_Real4
   Wnd_Spd_10m_Nwp_Pix = Missing_Value_Real4
   Wnd_Dir_10m_Nwp_Pix = Missing_Value_Real4
   Wnd_Spd_Cld_Top_Nwp_Pix = Missing_Value_Real4
   Wnd_Dir_Cld_Top_Nwp_Pix = Missing_Value_Real4
   Inversion_Base_Nwp_Pix = Missing_Value_Real4
   Inversion_Top_Nwp_Pix = Missing_Value_Real4
   Inversion_Strength_Nwp_Pix = Missing_Value_Real4
   Tpw_Above_Cloud_Nwp_Pix = Missing_Value_Real4
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
   deallocate(Sea_Ice_Frac_Nwp_Pix)
   deallocate(Weasd_Nwp_Pix)
   deallocate(Tpw_Nwp_Pix)
   deallocate(Ozone_Nwp_Pix)
   deallocate(Ttropo_Nwp_Pix)
   deallocate(Wnd_Spd_10m_Nwp_Pix)
   deallocate(Wnd_Dir_10m_Nwp_Pix)
   deallocate(Wnd_Spd_Cld_Top_Nwp_Pix)
   deallocate(Wnd_Dir_Cld_Top_Nwp_Pix)
   deallocate(Inversion_Base_Nwp_Pix)
   deallocate(Inversion_Top_Nwp_Pix)
   deallocate(Inversion_Strength_Nwp_Pix)
   deallocate(Tpw_Above_Cloud_Nwp_Pix)
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

   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
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
   endif

   if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
      allocate(Ch2_Counts(dim1,dim2))
   endif

   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
      allocate(Ch6_Counts(dim1,dim2))
   endif

end subroutine CREATE_REF_CHANNEL_ARRAYS

subroutine RESET_REF_CHANNEL_ARRAYS
   integer:: idx

   do idx = 1,Nchan_Clavrx
      if (allocated(Ch(idx)%Rad_Toa)) Ch(idx)%Rad_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa)) Ch(idx)%Bt_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Rad_Toa_Clear)) Ch(idx)%Rad_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Bt_Toa_Clear)) Ch(idx)%Bt_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Rad_Atm)) Ch(idx)%Rad_Atm = Missing_Value_Real4
      if (allocated(Ch(idx)%Rad_Atm_Dwn_Sfc)) Ch(idx)%Rad_Atm_Dwn_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Trans_Atm)) Ch(idx)%Trans_Atm = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa)) Ch(idx)%Ref_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa_Unnorm)) Ch(idx)%Ref_Toa_Unnorm = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Toa_Clear)) Ch(idx)%Ref_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Sfc)) Ch(idx)%Ref_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Toa)) Ch(idx)%Ref_Lunar_Toa = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Toa_Clear)) Ch(idx)%Ref_Lunar_Toa_Clear = Missing_Value_Real4
      if (allocated(Ch(idx)%Ref_Lunar_Sfc)) Ch(idx)%Ref_Lunar_Sfc = Missing_Value_Real4
      if (allocated(Ch(idx)%Opd)) Ch(idx)%Opd = Missing_Value_Real4
      if (allocated(Ch(idx)%Sfc_Emiss)) Ch(idx)%Sfc_Emiss = Missing_Value_Real4
      if (allocated(Ch(idx)%Sfc_Ref_White_Sky)) Ch(idx)%Sfc_Ref_White_Sky = Missing_Value_Real4
      if (allocated(Ch(idx)%Emiss_Tropo)) Ch(idx)%Emiss_Tropo = Missing_Value_Real4
      if (allocated(Ch(idx)%Unc)) Ch(idx)%Unc = Missing_Value_Int1
      if (allocated(Ch(idx)%CSBT_Mask)) Ch(idx)%CSBT_Mask = Missing_Value_Int1
      if (allocated(Ch(idx)%Opaque_Height)) Ch(idx)%Opaque_Height = Missing_Value_Real4
   enddo

   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
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
   endif

   if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
      Ch2_Counts = Missing_Value_Int2
   endif

   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
      Ch6_Counts = Missing_Value_Int2
   endif

end subroutine RESET_REF_CHANNEL_ARRAYS
subroutine DESTROY_REF_CHANNEL_ARRAYS

  if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
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
  endif

   if (Sensor%Chan_On_Flag_Default(2) == sym%YES) then
      deallocate(Ch2_Counts)
   endif

   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
      deallocate(Ch6_Counts)
   endif

end subroutine DESTROY_REF_CHANNEL_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_THERM_CHANNEL_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2

   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
      !allocate(Rad_Ch20_Ems(dim1,dim2))
       allocate(Ems_Ch20(dim1,dim2))
       allocate(Bt_Ch20_Std_3x3(dim1,dim2))
       allocate(Bt_Ch20_Median_3x3(dim1,dim2))
       allocate(Ems_Ch20_Median_3x3(dim1,dim2))
       allocate(Ems_Ch20_Std_Median_3x3(dim1,dim2))
       allocate(Bt_Ch20_Std_Median_3x3(dim1,dim2))
       allocate(Ch20_Counts_Filtered(dim1,dim2))
       allocate(Trans_Atm_Ch20_Solar_Rtm(dim1,dim2))
       allocate(Trans_Atm_Ch20_Solar_Total_Rtm(dim1,dim2))
       allocate(Ems_Ch20_Clear_Rtm(dim1,dim2))
   endif

   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      allocate(Bt_Ch31_Mean_3x3(dim1,dim2))
      allocate(Bt_Ch31_Max_3x3(dim1,dim2))
      allocate(Bt_Ch31_Min_3x3(dim1,dim2))
      allocate(Bt_Ch31_Std_3x3(dim1,dim2))
      allocate(Elem_Idx_Max_Bt_Ch31_3x3(dim1,dim2))
      allocate(Line_Idx_Max_Bt_Ch31_3x3(dim1,dim2))
      allocate(Elem_Idx_Min_Bt_Ch31_3x3(dim1,dim2))
      allocate(Line_Idx_Min_Bt_Ch31_3x3(dim1,dim2))
      allocate(Emiss_11um_Tropo_Nadir_Rtm(dim1,dim2))
   endif

end subroutine CREATE_THERM_CHANNEL_ARRAYS

subroutine RESET_THERM_CHANNEL_ARRAYS()

   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
!      Rad_Ch20_Ems = Missing_Value_Real4
       Ems_Ch20 = Missing_Value_Real4
       Bt_Ch20_Std_3x3 = Missing_Value_Real4
       Bt_Ch20_Median_3x3 = Missing_Value_Real4
       Ems_Ch20_Median_3x3 = Missing_Value_Real4
       Ems_Ch20_Std_Median_3x3 = Missing_Value_Real4
       Bt_Ch20_Std_Median_3x3 = Missing_Value_Real4
       Ch20_Counts_Filtered = Missing_Value_Real4
       Trans_Atm_Ch20_Solar_Rtm = Missing_Value_Real4
       Trans_Atm_Ch20_Solar_Total_Rtm = Missing_Value_Real4
       Ems_Ch20_Clear_Rtm = Missing_Value_Real4
   endif

   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      Bt_Ch31_Mean_3x3 = Missing_Value_Real4
      Bt_Ch31_Max_3x3 = Missing_Value_Real4
      Bt_Ch31_Min_3x3 = Missing_Value_Real4
      Bt_Ch31_Std_3x3 = Missing_Value_Real4
      Elem_Idx_Max_Bt_Ch31_3x3 = Missing_Value_Real4
      Line_Idx_Max_Bt_Ch31_3x3 = Missing_Value_Real4
      Elem_Idx_Min_Bt_Ch31_3x3 = Missing_Value_Real4
      Line_Idx_Min_Bt_Ch31_3x3 = Missing_Value_Real4
      Emiss_11um_Tropo_Nadir_Rtm = Missing_Value_Real4
   endif

end subroutine RESET_THERM_CHANNEL_ARRAYS
subroutine DESTROY_THERM_CHANNEL_ARRAYS()

   if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
!      deallocate(Rad_Ch20_Ems)
       deallocate(Ems_Ch20)
       deallocate(Bt_Ch20_Std_3x3)
       deallocate(Bt_Ch20_Median_3x3)
       deallocate(Ems_Ch20_Median_3x3)
       deallocate(Ems_Ch20_Std_Median_3x3)
       deallocate(Bt_Ch20_Std_Median_3x3)
       deallocate(Ch20_Counts_Filtered)
       deallocate(Trans_Atm_Ch20_Solar_Rtm)
       deallocate(Trans_Atm_Ch20_Solar_Total_Rtm)
       deallocate(Ems_Ch20_Clear_Rtm)
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      deallocate(Bt_Ch31_Mean_3x3)
      deallocate(Bt_Ch31_Max_3x3)
      deallocate(Bt_Ch31_Min_3x3)
      deallocate(Bt_Ch31_Std_3x3)
      deallocate(Elem_Idx_Max_Bt_Ch31_3x3)
      deallocate(Line_Idx_Max_Bt_Ch31_3x3)
      deallocate(Elem_Idx_Min_Bt_Ch31_3x3)
      deallocate(Line_Idx_Min_Bt_Ch31_3x3)
      deallocate(Emiss_11um_Tropo_Nadir_Rtm)
   endif

end subroutine DESTROY_THERM_CHANNEL_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_EXTRA_CHANNEL_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (index(Sensor%Sensor_Name,'IFF') > 0) then
           allocate(Bt_375um_Sounder(dim1,dim2))
           allocate(Bt_11um_Sounder(dim1,dim2))
           allocate(Bt_12um_Sounder(dim1,dim2))
   endif
   if (index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
           ! MJH HIRS/AVHRR aux fields.
           ! Does this belong in CREATE_EXTRA_CHANNEL_ARRAYS??
           allocate(Cld_Temp_Sounder(dim1,dim2))
           allocate(Cld_Press_Sounder(dim1,dim2))
           allocate(Cld_Height_Sounder(dim1,dim2))
           allocate(Cld_Emiss_Sounder(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(39) == sym%YES) then
           allocate(Ref_ChI1(2*dim1,2*dim2))
           allocate(Ref_Max_ChI1(dim1,dim2))
           allocate(Ref_Min_ChI1(dim1,dim2))
           allocate(Ref_Uni_ChI1(dim1,dim2))
           allocate(Ref_Mean_ChI1(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(40) == sym%YES) then
           allocate(Ref_ChI2(2*dim1,2*dim2))
           allocate(Ref_Max_ChI2(dim1,dim2))
           allocate(Ref_Min_ChI2(dim1,dim2))
           allocate(Ref_Uni_ChI2(dim1,dim2))
           allocate(Ref_Mean_ChI2(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(41) == sym%YES) then
           allocate(Ref_ChI3(2*dim1,2*dim2))
           allocate(Ref_Max_ChI3(dim1,dim2))
           allocate(Ref_Min_ChI3(dim1,dim2))
           allocate(Ref_Uni_ChI3(dim1,dim2))
           allocate(Ref_Mean_ChI3(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(42) == sym%YES) then
           allocate(Bt_ChI4(2*dim1,2*dim2))
           allocate(Bt_Max_ChI4(dim1,dim2))
           allocate(Bt_Min_ChI4(dim1,dim2))
           allocate(Bt_Uni_ChI4(dim1,dim2))
           allocate(Bt_Mean_ChI4(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(43) == sym%YES) then
           allocate(Bt_ChI5(2*dim1,2*dim2))
           allocate(Bt_Max_ChI5(dim1,dim2))
           allocate(Bt_Min_ChI5(dim1,dim2))
           allocate(Bt_Uni_ChI5(dim1,dim2))
           allocate(Bt_Mean_ChI5(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
           allocate(Ref_ChDNB_Lunar_Mean_3x3(dim1,dim2))
           allocate(Ref_ChDNB_Lunar_Max_3x3(dim1,dim2))
           allocate(Ref_ChDNB_Lunar_Min_3x3(dim1,dim2))
           allocate(Ref_ChDNB_Lunar_Std_3x3(dim1,dim2))
   endif
end subroutine CREATE_EXTRA_CHANNEL_ARRAYS

subroutine RESET_EXTRA_CHANNEL_ARRAYS()
      if (Sensor%Chan_On_Flag_Default(39) == sym%YES) Ref_ChI1 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(39) == sym%YES) Ref_Max_ChI1 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(39) == sym%YES) Ref_Min_ChI1 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(39) == sym%YES) Ref_Uni_ChI1 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(39) == sym%YES) Ref_Mean_ChI1 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(40) == sym%YES) Ref_ChI2 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(40) == sym%YES) Ref_Max_ChI2 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(40) == sym%YES) Ref_Min_ChI2 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(40) == sym%YES) Ref_Uni_ChI2 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(40) == sym%YES) Ref_Mean_ChI2 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(41) == sym%YES) Ref_ChI3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(41) == sym%YES) Ref_Max_ChI3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(41) == sym%YES) Ref_Min_ChI3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(41) == sym%YES) Ref_Uni_ChI3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(41) == sym%YES) Ref_Mean_ChI3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(42) == sym%YES) Bt_ChI4 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(42) == sym%YES) Bt_Max_ChI4 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(42) == sym%YES) Bt_Min_ChI4 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(42) == sym%YES) Bt_Uni_ChI4 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(42) == sym%YES) Bt_Mean_ChI4 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(43) == sym%YES) Bt_ChI5 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(43) == sym%YES) Bt_Max_ChI5 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(43) == sym%YES) Bt_Min_ChI5 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(43) == sym%YES) Bt_Uni_ChI5 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(43) == sym%YES) Bt_Mean_ChI5 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) Ref_ChDNB_Lunar_Mean_3x3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) Ref_ChDNB_Lunar_Max_3x3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) Ref_ChDNB_Lunar_Min_3x3 = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) Ref_ChDNB_Lunar_Std_3x3 = Missing_Value_Real4
      if (index(Sensor%Sensor_Name,'IFF') > 0) then
          Bt_375um_Sounder = Missing_Value_Real4
          Bt_11um_Sounder = Missing_Value_Real4
          Bt_12um_Sounder = Missing_Value_Real4
      endif
      if (index(Sensor%Sensor_Name,'AVHRR-IFF') > 0) then
          Cld_Temp_Sounder = Missing_Value_Real4 ! MJH
          Cld_Press_Sounder = Missing_Value_Real4
          Cld_Height_Sounder = Missing_Value_Real4 
          Cld_Emiss_Sounder = Missing_Value_Real4 
      endif
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
  if (allocated(Ref_Uni_ChI2)) deallocate(Ref_Uni_ChI2)
  if (allocated(Ref_Mean_ChI2)) deallocate(Ref_Mean_ChI2)
  if (allocated(Ref_ChI3)) deallocate(Ref_ChI3)
  if (allocated(Ref_Max_ChI3)) deallocate(Ref_Max_ChI3)
  if (allocated(Ref_Min_ChI3)) deallocate(Ref_Min_ChI3)
  if (allocated(Ref_Uni_ChI3)) deallocate(Ref_Uni_ChI3)
  if (allocated(Ref_Mean_ChI3)) deallocate(Ref_Mean_ChI3)
  if (allocated(Bt_ChI4)) deallocate(Bt_ChI4)
  if (allocated(Bt_Max_ChI4)) deallocate(Bt_Max_ChI4)
  if (allocated(Bt_Min_ChI4)) deallocate(Bt_Min_ChI4)
  if (allocated(Bt_Uni_ChI4)) deallocate(Bt_Uni_ChI4)
  if (allocated(Bt_Mean_ChI4)) deallocate(Bt_Mean_ChI4)
  if (allocated(Bt_ChI5)) deallocate(Bt_ChI5)
  if (allocated(Bt_Max_ChI5)) deallocate(Bt_Max_ChI5)
  if (allocated(Bt_Min_ChI5)) deallocate(Bt_Min_ChI5)
  if (allocated(Bt_Uni_ChI5)) deallocate(Bt_Uni_ChI5)
  if (allocated(Bt_Mean_ChI5)) deallocate(Bt_Mean_ChI5)
  if (allocated(Ref_ChDNB_Lunar_Mean_3x3)) deallocate(Ref_ChDNB_Lunar_Mean_3x3)
  if (allocated(Ref_ChDNB_Lunar_Min_3x3)) deallocate(Ref_ChDNB_Lunar_Min_3x3)
  if (allocated(Ref_ChDNB_Lunar_Max_3x3)) deallocate(Ref_ChDNB_Lunar_Max_3x3)
  if (allocated(Ref_ChDNB_Lunar_Std_3x3)) deallocate(Ref_ChDNB_Lunar_Std_3x3)
  if (allocated(Bt_375um_Sounder)) deallocate(Bt_375um_Sounder)
  if (allocated(Bt_11um_Sounder)) deallocate(Bt_11um_Sounder)
  if (allocated(Bt_12um_Sounder)) deallocate(Bt_12um_Sounder)
  if (allocated(Cld_Temp_Sounder)) deallocate(Cld_Temp_Sounder) ! MJH
  if (allocated(Cld_Press_Sounder)) deallocate(Cld_Press_Sounder)
  if (allocated(Cld_Height_Sounder)) deallocate(Cld_Height_Sounder)
  if (allocated(Cld_Emiss_Sounder)) deallocate(Cld_Emiss_Sounder)
end subroutine DESTROY_EXTRA_CHANNEL_ARRAYS

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_BTD_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      allocate(Btd_Ch20_Ch31(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      allocate(Btd_Ch20_Ch32(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      allocate(Btd_Ch31_Ch32(dim1,dim2))
      allocate(Btd_Ch31_Ch32_Std_3x3(dim1,dim2))
      allocate(Btd_Ch31_Ch32_Bt_Ch31_Max_3x3(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
      allocate(Btd_Ch31_Ch27_Std_3x3(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(29) == sym%YES) then
      allocate(Btd_Ch31_Ch29_Std_3x3(dim1,dim2))
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(33) == sym%YES) then
      allocate(Btd_Ch31_Ch33_Std_3x3(dim1,dim2))
   endif
end subroutine CREATE_BTD_ARRAYS
subroutine RESET_BTD_ARRAYS()
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      Btd_Ch20_Ch31 = Missing_Value_Real4
   endif
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      Btd_Ch20_Ch32 = Missing_Value_Real4
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      Btd_Ch31_Ch32 = Missing_Value_Real4
      Btd_Ch31_Ch32_Std_3x3 = Missing_Value_Real4
      Btd_Ch31_Ch32_Bt_Ch31_Max_3x3 = Missing_Value_Real4
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
      Btd_Ch31_Ch27_Std_3x3 = Missing_Value_Real4
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(29) == sym%YES) then
      Btd_Ch31_Ch29_Std_3x3 = Missing_Value_Real4
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(33) == sym%YES) then
      Btd_Ch31_Ch33_Std_3x3 = Missing_Value_Real4
   endif
end subroutine RESET_BTD_ARRAYS
subroutine DESTROY_BTD_ARRAYS()
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
      deallocate(Btd_Ch20_Ch31)
   endif
   if (Sensor%Chan_On_Flag_Default(20) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      deallocate(Btd_Ch20_Ch32)
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      deallocate(Btd_Ch31_Ch32)
      deallocate(Btd_Ch31_Ch32_Std_3x3)
      deallocate(Btd_Ch31_Ch32_Bt_Ch31_Max_3x3)
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
      deallocate(Btd_Ch31_Ch27_Std_3x3)
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(29) == sym%YES) then
      deallocate(Btd_Ch31_Ch29_Std_3x3)
   endif
   if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. Sensor%Chan_On_Flag_Default(33) == sym%YES) then
      deallocate(Btd_Ch31_Ch33_Std_3x3)
   endif
end subroutine DESTROY_BTD_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_SURFACE_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   allocate(Sfc%Land(dim1,dim2))
   allocate(Sfc%Land_Mask(dim1,dim2))
   allocate(Sfc%Coast(dim1,dim2))
   allocate(Sfc%Coast_Mask(dim1,dim2))
   allocate(Sfc%Coast_Mask_Nwp(dim1,dim2))
   allocate(Sfc%Glint_Mask(dim1,dim2))
   allocate(Sfc%Glint_Mask_Lunar(dim1,dim2))
   allocate(Sfc%Forward_Scatter_Mask(dim1,dim2))
   allocate(Sfc%Forward_Scatter_Mask_Lunar(dim1,dim2))
   allocate(Sfc%Desert_Mask(dim1,dim2))
   allocate(Sfc%City_Mask(dim1,dim2))
   allocate(Sfc%Volcano_Mask(dim1,dim2))
   allocate(Sfc%Snow_NWP(dim1,dim2))
   allocate(Sfc%Snow_OISST(dim1,dim2))
   allocate(Sfc%Snow_IMS(dim1,dim2))
   allocate(Sfc%Snow_GLOB(dim1,dim2))
   allocate(Sfc%Snow(dim1,dim2))
   allocate(Sfc%Sfc_Type(dim1,dim2))
   allocate(Sfc%Zsfc(dim1,dim2))
   allocate(Sfc%Zsfc_Hires(dim1,dim2))
end subroutine CREATE_SURFACE_ARRAYS
subroutine RESET_SURFACE_ARRAYS
   Sfc%Land = Missing_Value_Int1
   Sfc%Land_Mask = Missing_Value_Int1
   Sfc%Coast = Missing_Value_Int1
   Sfc%Coast_Mask = Missing_Value_Int1
   Sfc%Coast_Mask_Nwp = Missing_Value_Int1
   Sfc%Glint_Mask = Missing_Value_Int1
   Sfc%Glint_Mask_Lunar = Missing_Value_Int1
   Sfc%Forward_Scatter_Mask = Missing_Value_Int1
   Sfc%Forward_Scatter_Mask_Lunar = Missing_Value_Int1
   Sfc%Desert_Mask = Missing_Value_Int1
   Sfc%City_Mask = Missing_Value_Int1
   Sfc%Volcano_Mask = Missing_Value_Int1
   Sfc%Snow_NWP = Missing_Value_Int1
   Sfc%Snow_OISST = Missing_Value_Int1
   Sfc%Snow_IMS = Missing_Value_Int1
   Sfc%Snow_GLOB = Missing_Value_Int1
   Sfc%Snow = Missing_Value_Int1
   Sfc%Sfc_Type = Missing_Value_Int1
   Sfc%Zsfc = Missing_Value_Real4
   Sfc%Zsfc_Hires = Missing_Value_Real4
end subroutine RESET_SURFACE_ARRAYS
subroutine DESTROY_SURFACE_ARRAYS
   deallocate(Sfc%Land)
   deallocate(Sfc%Land_Mask)
   deallocate(Sfc%Coast)
   deallocate(Sfc%Coast_Mask)
   deallocate(Sfc%Coast_Mask_Nwp)
   deallocate(Sfc%Glint_Mask)
   deallocate(Sfc%Glint_Mask_Lunar)
   deallocate(Sfc%Forward_Scatter_Mask)
   deallocate(Sfc%Forward_Scatter_Mask_Lunar)
   deallocate(Sfc%Desert_Mask)
   deallocate(Sfc%City_Mask)
   deallocate(Sfc%Volcano_Mask)
   deallocate(Sfc%Snow_NWP)
   deallocate(Sfc%Snow_OISST)
   deallocate(Sfc%Snow_IMS)
   deallocate(Sfc%Snow_GLOB)
   deallocate(Sfc%Snow)
   deallocate(Sfc%Sfc_Type)
   deallocate(Sfc%Zsfc)
   deallocate(Sfc%Zsfc_Hires)
end subroutine DESTROY_SURFACE_ARRAYS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_ACHA_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES) then

    allocate(ACHA%Tc(dim1,dim2)) 
    allocate(ACHA%Ec(dim1,dim2)) 
    allocate(ACHA%Pc(dim1,dim2)) 
    allocate(ACHA%Zc(dim1,dim2)) 
    allocate(ACHA%Zc_Top(dim1,dim2)) 
    allocate(ACHA%Pc_Top(dim1,dim2)) 
    allocate(Pc_Top1_Aux(dim1,dim2)) 
    allocate(Pc_Top2_Aux(dim1,dim2)) 
    allocate(ACHA%Zc_Base(dim1,dim2)) 
    allocate(ACHA%Pc_Base(dim1,dim2)) 
    allocate(ACHA%Beta(dim1,dim2)) 
    allocate(ACHA%Tau(dim1,dim2)) 
    allocate(ACHA%Reff(dim1,dim2)) 
    allocate(ACHA%Tc_Uncertainty(dim1,dim2)) 
    allocate(ACHA%Ec_Uncertainty(dim1,dim2)) 
    allocate(ACHA%Beta_Uncertainty(dim1,dim2)) 
    allocate(ACHA%Zc_Uncertainty(dim1,dim2)) 
    allocate(ACHA%Pc_Uncertainty(dim1,dim2)) 
    allocate(ACHA%Lower_Tc_Uncertainty(dim1,dim2)) 
    allocate(ACHA%Lower_Pc_Uncertainty(dim1,dim2)) 
    allocate(ACHA%Lower_Zc_Uncertainty(dim1,dim2)) 
    allocate(Pc_Uncertainty1_Aux(dim1,dim2)) 
    allocate(Pc_Uncertainty2_Aux(dim1,dim2)) 
    allocate(ACHA%Alt(dim1,dim2)) 
    allocate(ACHA%Base_Alt(dim1,dim2)) 
    allocate(ACHA%Cost(dim1,dim2)) 
    allocate(Cost_Aux(dim1,dim2)) 
    allocate(ACHA%Lower_Pc(dim1,dim2)) 
    allocate(ACHA%Lower_Zc(dim1,dim2)) 
    allocate(ACHA%Lower_Tc(dim1,dim2)) 
    allocate(ACHA%Processing_Order(dim1,dim2)) 
    allocate(ACHA%Inversion_Flag(dim1,dim2)) 
    allocate(ACHA%Quality_Flag(dim1,dim2)) 
    allocate(ACHA%Meta_Data(dim1,dim2)) 
    allocate(ACHA%OE_Quality_Flags(4,dim1,dim2)) 
    allocate(ACHA%Packed_Quality_Flags(dim1,dim2)) 
    allocate(ACHA%Packed_Meta_Data_Flags(dim1,dim2)) 
    allocate(ACHA%Conv_Cld_Prob(dim1,dim2)) 
    allocate(ACHA%Supercooled_Cld_Prob(dim1,dim2)) 
    allocate(ACHA%Base_Quality_Flag(dim1,dim2))
    allocate(ACHA%Ec_67um(dim1,dim2))
    allocate(ACHA%Ec_85um(dim1,dim2))
    allocate(ACHA%Ec_11um(dim1,dim2))
    allocate(ACHA%Ec_12um(dim1,dim2))
    allocate(ACHA%Ec_133um(dim1,dim2))
   endif

   !--- these accumulate through the whole image, do not reset with each segment
   ACHA%Processed_Count = 0
   ACHA%Valid_Count = 0
   ACHA%Success_Fraction = Missing_Value_Real4

end subroutine CREATE_ACHA_ARRAYS

subroutine RESET_ACHA_ARRAYS()

    ACHA%Tc = Missing_Value_Real4
    ACHA%Ec = Missing_Value_Real4
    ACHA%Pc = Missing_Value_Real4
    ACHA%Zc = Missing_Value_Real4
    ACHA%Zc_Top = Missing_Value_Real4
    ACHA%Pc_Top = Missing_Value_Real4
    Pc_Top1_Aux = Missing_Value_Real4
    Pc_Top2_Aux = Missing_Value_Real4
    ACHA%Zc_Base  = Missing_Value_Real4
    ACHA%Pc_Base  = Missing_Value_Real4
    ACHA%Beta = Missing_Value_Real4
    ACHA%Tau = Missing_Value_Real4
    ACHA%Reff = Missing_Value_Real4
    ACHA%Tc_Uncertainty = Missing_Value_Real4
    ACHA%Ec_Uncertainty = Missing_Value_Real4
    ACHA%Beta_Uncertainty = Missing_Value_Real4
    ACHA%Zc_Uncertainty = Missing_Value_Real4
    ACHA%Pc_Uncertainty = Missing_Value_Real4
    ACHA%Lower_Tc_Uncertainty = Missing_Value_Real4
    ACHA%Lower_Pc_Uncertainty = Missing_Value_Real4
    ACHA%Lower_Zc_Uncertainty = Missing_Value_Real4
    Pc_Uncertainty1_Aux = Missing_Value_Real4
    Pc_Uncertainty2_Aux = Missing_Value_Real4
    ACHA%Alt = Missing_Value_Real4
    ACHA%Base_Alt = Missing_Value_Real4
    ACHA%Cost = Missing_Value_Real4
    Cost_Aux = Missing_Value_Real4
    ACHA%Lower_Pc = Missing_Value_Real4
    ACHA%Lower_Zc = Missing_Value_Real4
    ACHA%Lower_Tc = Missing_Value_Real4
    ACHA%Processing_Order = Missing_Value_Int1
    ACHA%Inversion_Flag = Missing_Value_Int1
    ACHA%Quality_Flag = Missing_Value_Int1
    ACHA%Meta_Data = 0
    ACHA%OE_Quality_Flags = 0
    ACHA%Packed_Quality_Flags = 0
    ACHA%Packed_Meta_Data_Flags = 0
    ACHA%Conv_Cld_Prob = Missing_Value_Real4
    ACHA%Supercooled_Cld_Prob = Missing_Value_Real4
    ACHA%base_Quality_Flag = 1   ! Missing_Value_Int1
    ACHA%Ec_67um = Missing_Value_Real4
    ACHA%Ec_85um = Missing_Value_Real4
    ACHA%Ec_11um = Missing_Value_Real4
    ACHA%Ec_12um = Missing_Value_Real4
    ACHA%Ec_133um = Missing_Value_Real4

end subroutine RESET_ACHA_ARRAYS

subroutine DESTROY_ACHA_ARRAYS()

    deallocate(ACHA%Tc) 
    deallocate(ACHA%Ec) 
    deallocate(ACHA%Pc) 
    deallocate(ACHA%Zc) 
    deallocate(ACHA%Zc_Top) 
    deallocate(ACHA%Pc_Top) 
    deallocate(Pc_Top1_Aux) 
    deallocate(Pc_Top2_Aux) 
    deallocate(ACHA%Zc_Base) 
    deallocate(ACHA%Pc_Base) 
    deallocate(ACHA%Beta) 
    deallocate(ACHA%Tau) 
    deallocate(ACHA%Reff) 
    deallocate(ACHA%Tc_Uncertainty) 
    deallocate(ACHA%Ec_Uncertainty) 
    deallocate(ACHA%Beta_Uncertainty) 
    deallocate(ACHA%Zc_Uncertainty) 
    deallocate(ACHA%Pc_Uncertainty) 
    deallocate(ACHA%Lower_Tc_Uncertainty) 
    deallocate(ACHA%Lower_Zc_Uncertainty) 
    deallocate(ACHA%Lower_Pc_Uncertainty) 
    deallocate(Pc_Uncertainty1_Aux) 
    deallocate(Pc_Uncertainty2_Aux) 
    deallocate(ACHA%Alt) 
    deallocate(ACHA%Base_Alt) 
    deallocate(ACHA%Cost) 
    deallocate(Cost_Aux) 
    deallocate(ACHA%Lower_Pc) 
    deallocate(ACHA%Lower_Zc) 
    deallocate(ACHA%Lower_Tc) 
    deallocate(ACHA%Processing_Order) 
    deallocate(ACHA%Inversion_Flag) 
    deallocate(ACHA%Quality_Flag) 
    deallocate(ACHA%Meta_Data) 
    deallocate(ACHA%OE_Quality_Flags) 
    deallocate(ACHA%Packed_Quality_Flags) 
    deallocate(ACHA%Packed_Meta_Data_Flags) 
    deallocate(ACHA%Conv_Cld_Prob) 
    deallocate(ACHA%Supercooled_Cld_Prob) 
    deallocate(ACHA%base_Quality_Flag)
    deallocate(ACHA%Ec_67um) 
    deallocate(ACHA%Ec_85um) 
    deallocate(ACHA%Ec_11um) 
    deallocate(ACHA%Ec_12um) 
    deallocate(ACHA%Ec_133um) 

end subroutine DESTROY_ACHA_ARRAYS
!------------------------------------------------------------------------------
! Cloud Cover Layers (CCL) data structure routines
!------------------------------------------------------------------------------
subroutine CREATE_CCL_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES) then
    allocate (CCL%Cld_Layer(dim1,dim2))
    allocate (CCL%Cloud_Fraction(dim1,dim2))
    allocate (CCL%Cloud_Fraction_Uncer(dim1,dim2))
    allocate (CCL%High_Cloud_Fraction(dim1,dim2))
    allocate (CCL%Mid_Cloud_Fraction(dim1,dim2))
    allocate (CCL%Low_Cloud_Fraction(dim1,dim2))
   endif
end subroutine CREATE_CCL_ARRAYS
subroutine RESET_CCL_ARRAYS()
    if (allocated(CCL%Cld_Layer)) CCL%Cld_Layer = Missing_Value_Int1
    if (allocated(CCL%Cloud_Fraction)) CCL%Cloud_Fraction = Missing_Value_Real4
    if (allocated(CCL%Cloud_Fraction_Uncer)) CCL%Cloud_Fraction_Uncer = Missing_Value_Real4
    if (allocated(CCL%High_Cloud_Fraction)) CCL%High_Cloud_Fraction = Missing_Value_Real4
    if (allocated(CCL%Mid_Cloud_Fraction)) CCL%Mid_Cloud_Fraction = Missing_Value_Real4
    if (allocated(CCL%Low_Cloud_Fraction)) CCL%Low_Cloud_Fraction = Missing_Value_Real4
end subroutine RESET_CCL_ARRAYS
subroutine DESTROY_CCL_ARRAYS()
    if (allocated(CCL%Cld_Layer)) deallocate (CCL%Cld_Layer)
    if (allocated(CCL%Cloud_Fraction)) deallocate (CCL%Cloud_Fraction)
    if (allocated(CCL%Cloud_Fraction_Uncer)) deallocate (CCL%Cloud_Fraction_Uncer)
    if (allocated(CCL%High_Cloud_Fraction)) deallocate (CCL%High_Cloud_Fraction)
    if (allocated(CCL%Mid_Cloud_Fraction)) deallocate (CCL%Mid_Cloud_Fraction)
    if (allocated(CCL%Low_Cloud_Fraction)) deallocate (CCL%Low_Cloud_Fraction)
end subroutine DESTROY_CCL_ARRAYS
!------------------------------------------------------------------------------
! Automatic Surface Observing System (ASOS) data structure routines
!------------------------------------------------------------------------------
subroutine CREATE_ASOS_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES) then
    allocate(ASOS%Code(dim1,dim2))
    allocate(ASOS%ECA(dim1,dim2))
    allocate(ASOS%Zmax(dim1,dim2))
    allocate(ASOS%Zmin(dim1,dim2))
   endif
end subroutine CREATE_ASOS_ARRAYS
subroutine RESET_ASOS_ARRAYS()
   if (allocated(ASOS%Code)) ASOS%Code = Missing_Value_Int1
   if (allocated(ASOS%ECA)) ASOS%ECA = Missing_Value_Real4
   if (allocated(ASOS%Zmax)) ASOS%Zmax = Missing_Value_Real4
   if (allocated(ASOS%Zmin)) ASOS%Zmin = Missing_Value_Real4
end subroutine RESET_ASOS_ARRAYS
subroutine DESTROY_ASOS_ARRAYS()
   if (allocated(ASOS%Code)) deallocate(ASOS%Code)
   if (allocated(ASOS%ECA)) deallocate(ASOS%ECA)
   if (allocated(ASOS%Zmax)) deallocate(ASOS%Zmax)
   if (allocated(ASOS%Zmin)) deallocate(ASOS%Zmin)
end subroutine DESTROY_ASOS_ARRAYS

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
subroutine CREATE_DCOMP_ARRAYS(dim1,dim2)
   integer, intent(in):: dim1, dim2
   if (Cld_Flag == sym%YES) then
      allocate(Tau_DCOMP(dim1,dim2))
      allocate(Tau_DCOMP_1(dim1,dim2))
      allocate(Tau_DCOMP_2(dim1,dim2))
      allocate(Tau_DCOMP_3(dim1,dim2))
      allocate(Tau_Aux(dim1,dim2))
      allocate(Reff_Aux(dim1,dim2))
      allocate(Tau_DCOMP_Ap(dim1,dim2))
      allocate(Vis_Ref_Fm(dim1,dim2))
      allocate(Reff_DCOMP(dim1,dim2))
      allocate(Reff_DCOMP_1(dim1,dim2))
      allocate(Reff_DCOMP_2(dim1,dim2))
      allocate(Reff_DCOMP_3(dim1,dim2))
      allocate(Lwp_DCOMP(dim1,dim2))
      allocate(Iwp_DCOMP(dim1,dim2))
      allocate(Iwp_Tau_DCOMP(dim1,dim2))
      allocate(Cwp_DCOMP(dim1,dim2))
      allocate(Cwp_Ice_Layer_DCOMP(dim1,dim2))
      allocate(Cwp_Water_Layer_DCOMP(dim1,dim2))
      allocate(Cwp_Scwater_Layer_DCOMP(dim1,dim2))
      allocate(Rain_Rate_DCOMP(dim1,dim2))
      allocate(Hcld_DCOMP(dim1,dim2))
      allocate(Cdnc_DCOMP(dim1,dim2))
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
      Tau_DCOMP_1 = Missing_Value_Real4
      Tau_DCOMP_2 = Missing_Value_Real4
      Tau_DCOMP_3 = Missing_Value_Real4
      Tau_Aux = Missing_Value_Real4
      Reff_Aux = Missing_Value_Real4
      Tau_DCOMP_Ap = Missing_Value_Real4
      Vis_Ref_Fm = Missing_Value_Real4
      Reff_DCOMP = Missing_Value_Real4
      Reff_DCOMP_1 = Missing_Value_Real4
      Reff_DCOMP_2 = Missing_Value_Real4
      Reff_DCOMP_3 = Missing_Value_Real4
      Lwp_DCOMP = Missing_Value_Real4
      Iwp_DCOMP = Missing_Value_Real4
      Iwp_Tau_DCOMP = Missing_Value_Real4
      Cwp_DCOMP = Missing_Value_Real4
      Cwp_Ice_Layer_DCOMP = Missing_Value_Real4
      Cwp_Water_Layer_DCOMP = Missing_Value_Real4
      Cwp_Scwater_Layer_DCOMP = Missing_Value_Real4
      Rain_Rate_DCOMP = Missing_Value_Real4
      Hcld_DCOMP = Missing_Value_Real4
      Cdnc_DCOMP = Missing_Value_Real4
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
      deallocate(Tau_DCOMP_1)
      deallocate(Tau_DCOMP_2)
      deallocate(Tau_DCOMP_3)   
      deallocate(Tau_Aux)
      deallocate(Reff_Aux)
      deallocate(Tau_DCOMP_Ap)
      deallocate(Vis_Ref_Fm)
      deallocate(Reff_DCOMP)
      deallocate(Reff_DCOMP_1)
      deallocate(Reff_DCOMP_2)
      deallocate(Reff_DCOMP_3)
      deallocate(Lwp_DCOMP)
      deallocate(Iwp_DCOMP)
      deallocate(Iwp_Tau_DCOMP)
      deallocate(Cwp_DCOMP)
      deallocate(Cwp_Ice_Layer_DCOMP)
      deallocate(Cwp_Water_Layer_DCOMP)
      deallocate(Cwp_Scwater_Layer_DCOMP)
      deallocate(Rain_Rate_DCOMP)
      deallocate(Hcld_DCOMP)
      deallocate(Cdnc_DCOMP)
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
   if (Sasrab_Flag == sym%YES) then
      allocate(Insolation_All_Sky(dim1,dim2))
      allocate(Insolation_All_Sky_Diffuse(dim1,dim2))
      allocate(Insolation_Clear_Sky(dim1,dim2))
      allocate(Insolation_Cld_Opd(dim1,dim2))
      allocate(Insolation_Aer_Opd(dim1,dim2))
   endif
end subroutine CREATE_SASRAB_ARRAYS
subroutine RESET_SASRAB_ARRAYS()
   if (Sasrab_Flag == sym%YES) then
      Insolation_All_Sky = Missing_Value_Real4
      Insolation_All_Sky_Diffuse = Missing_Value_Real4
      Insolation_Clear_Sky = Missing_Value_Real4
      Insolation_Cld_Opd = Missing_Value_Real4
      Insolation_Aer_Opd = Missing_Value_Real4
   endif
end subroutine RESET_SASRAB_ARRAYS
subroutine DESTROY_SASRAB_ARRAYS()
   if (sasrab_Flag == sym%YES) then
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
   
      allocate(Olr(dim1,dim2))
   
end subroutine CREATE_OLR_ARRAYS
subroutine RESET_OLR_ARRAYS
   
      Olr = Missing_Value_Real4
   
end subroutine RESET_OLR_ARRAYS
subroutine DESTROY_OLR_ARRAYS
   
      deallocate(Olr)
   
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
     allocate(Prior_Cld_Probability(dim1,dim2))
     allocate(Bayes_Mask_Sfc_Type_Global(dim1,dim2))
  endif
  !Fix needed to get around issue when GFS turned off for AVHRR
  if (trim(Sensor%Sensor_Name) == 'AVHRR-1' .or. &
      trim(Sensor%Sensor_Name) == 'AVHRR-2' .or. &
      trim(Sensor%Sensor_Name) == 'AVHRR-3') then
     if (.not. allocated(Cld_Mask_Aux)) allocate(Cld_Mask_Aux(dim1,dim2))  
  endif
  
end subroutine CREATE_CLOUD_MASK_ARRAYS
subroutine RESET_CLOUD_MASK_ARRAYS()
  Cld_Mask_Qf = Missing_Value_Int1
  if (Cld_Flag == sym%YES) then
     Cld_Mask = Missing_Value_Int1
     Cld_Mask_Aux = Missing_Value_Int1
     Adj_Pix_Cld_Mask = Missing_Value_Int1
     Posterior_Cld_Probability = Missing_Value_Real4
     Prior_Cld_Probability = Missing_Value_Real4
     Cld_Test_Vector_Packed = 0
     Bayes_Mask_Sfc_Type_Global = Missing_Value_Int1
  endif
  !Fix needed to get around issue when GFS turned off for AVHRR
  if (trim(Sensor%Sensor_Name) == 'AVHRR-1' .or. &
      trim(Sensor%Sensor_Name) == 'AVHRR-2' .or. &
      trim(Sensor%Sensor_Name) == 'AVHRR-3') then
     Cld_Mask_Aux = Missing_Value_Int1
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
     deallocate(Prior_Cld_Probability)
     deallocate(Bayes_Mask_Sfc_Type_Global)
  endif
  !Fix needed to get around issue when GFS turned off for AVHRR
  if (trim(Sensor%Sensor_Name) == 'AVHRR-1' .or. &
      trim(Sensor%Sensor_Name) == 'AVHRR-2' .or. &
      trim(Sensor%Sensor_Name) == 'AVHRR-3') then
     if (allocated(Cld_Mask_Aux)) deallocate(Cld_Mask_Aux)  
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
     allocate(Cld_Phase_IR(dim1,dim2))
     allocate(Cld_Type_IR(dim1,dim2))
     allocate(Ctp_Multilayer_Flag(dim1,dim2))
     allocate(Zc_Aux(dim1,dim2))
  endif
end subroutine CREATE_CLOUD_TYPE_ARRAYS
subroutine RESET_CLOUD_TYPE_ARRAYS()
  if (Cld_Flag == sym%YES) then
      Cld_Phase = Missing_Value_Int1
      Cld_Type = Missing_Value_Int1
      Cld_Phase_IR = Missing_Value_Int1
      Cld_Type_IR = Missing_Value_Int1
      Cld_Phase_Aux = Missing_Value_Int1
      Cld_Type_Aux = Missing_Value_Int1
      Ctp_Multilayer_Flag = Missing_Value_Int1
      Zc_Aux = Missing_Value_Real4
  endif
end subroutine RESET_CLOUD_TYPE_ARRAYS
subroutine DESTROY_CLOUD_TYPE_ARRAYS
  if (Cld_Flag == sym%YES) then
     deallocate(Cld_Type_Aux)
     deallocate(Cld_Phase_Aux)
     deallocate(Cld_Phase)
     deallocate(Cld_Type)
     deallocate(Cld_Phase_IR)
     deallocate(Cld_Type_IR)
     deallocate(Ctp_Multilayer_Flag)
     deallocate(Zc_Aux)
  endif
end subroutine DESTROY_CLOUD_TYPE_ARRAYS
!-----------------------------------------------------------
!
!-----------------------------------------------------------
subroutine CREATE_DIAGNOSTIC_ARRAYS(dim1,dim2)
  integer, intent(in):: dim1, dim2
  allocate(Temp_Pix_Array_1(dim1,dim2))
  allocate(Temp_Pix_Array_2(dim1,dim2))
  allocate(Temp_Pix_Array_3(dim1,dim2))
  allocate(Diag_Pix_Array_1(dim1,dim2))
  allocate(Diag_Pix_Array_2(dim1,dim2))
  allocate(Diag_Pix_Array_3(dim1,dim2))
  allocate(Missing_Pixel_Array_Real4(dim1,dim2))
  allocate(One_Byte_Temp(dim1,dim2))
  allocate(Two_Byte_Temp(dim1,dim2))
end subroutine CREATE_DIAGNOSTIC_ARRAYS
subroutine RESET_DIAGNOSTIC_ARRAYS()
  Temp_Pix_Array_1 = Missing_Value_Real4
  Temp_Pix_Array_2 = Missing_Value_Real4
  Temp_Pix_Array_3 = Missing_Value_Real4
  Diag_Pix_Array_1 = Missing_Value_Real4
  Diag_Pix_Array_2 = Missing_Value_Real4
  Diag_Pix_Array_3 = Missing_Value_Real4
  Missing_Pixel_Array_Real4 = Missing_Value_Real4
  One_Byte_Temp = Missing_Value_Int1
  Two_Byte_Temp = Missing_Value_Int2
end subroutine RESET_DIAGNOSTIC_ARRAYS
subroutine DESTROY_DIAGNOSTIC_ARRAYS()
  deallocate(Temp_Pix_Array_1)
  deallocate(Temp_Pix_Array_2)
  deallocate(Temp_Pix_Array_3)
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
  allocate(Ndsi_Sfc(dim1,dim2))
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
  Ndsi_Sfc = Missing_Value_Real4
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
  if (allocated(Ndsi_Sfc)) deallocate(Ndsi_Sfc)
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
    allocate(Pc_H2O(dim1,dim2))
    allocate(Tc_H2O(dim1,dim2))
    allocate(Zc_H2O(dim1,dim2))
    allocate(Zclr_H2O_Peak(dim1,dim2))
    allocate(Zc_CO2IRW(dim1,dim2))
    allocate(Pc_CO2IRW(dim1,dim2))
    allocate(Tc_CO2IRW(dim1,dim2))
    allocate(Tc_Co2(dim1,dim2))
    allocate(Pc_Co2(dim1,dim2))
    allocate(Zc_Co2(dim1,dim2))
    allocate(Ec_Co2(dim1,dim2))
    allocate(Tc_Cirrus_Background(dim1,dim2))
    allocate(Zc_Cirrus_Background(dim1,dim2))
  endif
end subroutine CREATE_CLOUD_PROD_ARRAYS
subroutine RESET_CLOUD_PROD_ARRAYS()
  if (Cld_Flag == sym%YES) then
     Pc_Opaque_Cloud = Missing_Value_Real4
     Zc_Opaque_Cloud = Missing_Value_Real4
     Tc_Opaque_Cloud = Missing_Value_Real4
     Pc_H2O = Missing_Value_Real4
     Tc_H2O = Missing_Value_Real4
     Zc_H2O = Missing_Value_Real4
     Zclr_H2O_Peak = Missing_Value_Real4
     Pc_CO2IRW = Missing_Value_Real4
     Tc_CO2IRW = Missing_Value_Real4
     Zc_CO2IRW = Missing_Value_Real4
     Tc_Co2 = Missing_Value_Real4
     Pc_Co2 = Missing_Value_Real4
     Zc_Co2 = Missing_Value_Real4
     Ec_Co2 = Missing_Value_Real4
     Tc_Cirrus_Background = Missing_Value_Real4
     Zc_Cirrus_Background = Missing_Value_Real4
  endif
end subroutine RESET_CLOUD_PROD_ARRAYS
subroutine DESTROY_CLOUD_PROD_ARRAYS()
  if (Cld_Flag == sym%YES) then
     deallocate(Pc_Opaque_Cloud)
     deallocate(Zc_Opaque_Cloud)
     deallocate(Tc_Opaque_Cloud)
     deallocate(Pc_H2O)
     deallocate(Tc_H2O)
     deallocate(Zc_H2O)
     deallocate(Zclr_H2O_Peak)
     deallocate(Zc_CO2IRW)
     deallocate(Tc_CO2IRW)
     deallocate(Pc_CO2IRW)
     deallocate(Tc_Co2)
     deallocate(Pc_Co2)
     deallocate(Zc_Co2)
     deallocate(Ec_Co2)
     deallocate(Tc_Cirrus_Background)
     deallocate(Zc_Cirrus_Background)
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
