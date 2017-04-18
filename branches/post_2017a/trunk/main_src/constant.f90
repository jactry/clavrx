! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: constant.f90 (src)
!       CONSTANTS (program)
!
! PURPOSE: store and serve various constants for use in the CLAVR-x system
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE public
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED public USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! (c) This code is copyrighted by the author and all NOAA restrictions apply
! 
! Reference: CLAVR-x system description document
!
! Dependencies:  None
!
! Calling Sequece:
!   use CONSTANTS
!
! Public Routines within in this Module: None
!--------------------------------------------------------------------------------------
module CONSTANTS
  implicit none
  
  integer, parameter, public:: int1 = selected_int_kind(1)
  integer, parameter, public:: int2 = selected_int_kind(3)
  integer, parameter, public:: int4 = selected_int_kind(8)
  integer, parameter, public:: int8 = selected_int_kind(10)
  integer, parameter, public:: real4 = selected_real_kind(6,37)
  integer, parameter, public:: real8 = selected_real_kind(15,307)
  integer, parameter, public:: ipre = real4

  !--- Common Numerical Constants
  real (kind=real4), parameter, public:: g = 9.8
  real (kind=real4), parameter, public:: pi = 3.14159265
  real (kind=real4), parameter, public:: dtor = pi/180.0

  !--- Missing Values
  real (kind=real4), parameter, public:: Missing_Value_Real4 = -999.0
  real (kind=real8), parameter, public:: Missing_Value_Real8 = -999.0
  integer(kind=int1), parameter, public:: Missing_Value_Int1 = -128
  integer(kind=int2), parameter, public:: Missing_Value_Int2 = -32768
  integer(kind=int4), parameter, public:: Missing_Value_Int4 = -999
  real(kind=real4), parameter, public:: No_Attribute_Missing_Value = -888.0

  !--- other often used constants
  real(kind=real4), parameter, public:: Day_Solzen_Thresh_Mask = 85.0
  real (kind=real4), parameter, public:: TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH = 60.0

  !--- parameters used in retrievals
  real, parameter, public::  Cossolzen_Min_Solar =  0.1    !solar angle limit for daytime

  !--- channel numbers
  integer, parameter, public:: Nchan_Avhrr = 6
  integer, parameter, public:: Nchan_Clavrx = 45

  !--- maximum number of cloud mask tests - used to dimension arrays
  integer, parameter, public:: Max_Num_Cld_Tests = 32
  integer, parameter, public:: Max_Num_Cld_Test_Bytes = 7

  !--- error prompt id
  character(*), parameter, public :: EXE_PROMPT = "CLAVR-x>> "

  !--- useful constants
  integer(kind=int4), public, parameter :: MAX_STR_LEN = 256
  integer(kind=int4), public, parameter :: LARGE_HDF_NUMBER = 999999999

  !--- cvs strings to be written as attributes
  character(120), public :: ACHA_Version
  character(120), public :: DCOMP_Version
  character(120), public :: Cloud_Mask_Version
  character(120), public :: Cloud_Mask_Lut_Version
  character(120), public :: Cloud_Type_Version
  character(120), public :: Cloud_Type_IR_Version
  
  !--- define sds names in hdf files of relevant static ancillary data
  character(*), parameter, public :: SFC_TYPE_SDS_NAME = "surface_type"
  character(*), parameter, public :: COAST_MASK_SDS_NAME = "coast_mask"
  character(*), parameter, public :: VOLCANO_MASK_SDS_NAME = "volcano_mask"
  character(*), parameter, public :: LAND_MASK_SDS_NAME = "land_sea_mask"
  character(*), parameter, public :: SURFACE_ELEV_SDS_NAME = "surface_elevation"
  character(*), parameter, public ::  &
                     MODIS_ALB_0_66_SDS_NAME="Albedo_Map_0.659"
  character(*), parameter, public ::  &
                     MODIS_ALB_0_86_SDS_NAME="Albedo_Map_0.858"
  character(*), parameter, public ::  &
                     MODIS_ALB_1_24_SDS_NAME="Albedo_Map_1.24"
  character(*), parameter, public ::  &
                     MODIS_ALB_1_64_SDS_NAME="Albedo_Map_1.64"
  character(*), parameter, public ::  &
                     MODIS_ALB_2_13_SDS_NAME="Albedo_Map_2.13"
  character(*), parameter, public :: SNOW_MASK_SDS_NAME = "snow_ice_cover"

  !--- define a structure of symbols of clarity
  TYPE, public:: symbol_struct
    integer(kind=int1) :: CLOUDY = 3
    integer(kind=int1) :: PROB_CLOUDY = 2
    integer(kind=int1) :: PROB_CLEAR = 1
    integer(kind=int1) :: CLEAR = 0
    
    integer(kind=int1) :: CLEAR_BINARY = 0
    integer(kind=int1) :: CLOUDY_BINARY = 1

    integer(kind=int1) :: CLEAR_TYPE = 0
    integer(kind=int1) :: PROB_CLEAR_TYPE = 1
    integer(kind=int1) :: FOG_TYPE = 2
    integer(kind=int1) :: WATER_TYPE = 3
    integer(kind=int1) :: SUPERCOOLED_TYPE = 4
    integer(kind=int1) :: MIXED_TYPE = 5
    integer(kind=int1) :: OPAQUE_ICE_TYPE = 6
    integer(kind=int1) :: TICE_TYPE = 6
    integer(kind=int1) :: CIRRUS_TYPE = 7
    integer(kind=int1) :: OVERLAP_TYPE = 8
    integer(kind=int1) :: OVERSHOOTING_TYPE = 9
    integer(kind=int1) :: UNKNOWN_TYPE = 10
    integer(kind=int1) :: DUST_TYPE = 11
    integer(kind=int1) :: SMOKE_TYPE = 12
    integer(kind=int1) :: FIRE_TYPE = 13
                                                  
    integer(kind=int1) :: CLEAR_PHASE = 0
    integer(kind=int1) :: WATER_PHASE = 1
    integer(kind=int1) :: SUPERCOOLED_PHASE = 2
    integer(kind=int1) :: MIXED_PHASE = 3
    integer(kind=int1) :: ICE_PHASE = 4
    integer(kind=int1) :: UNKNOWN_PHASE = 5
                                                                                
    integer(kind=int1) :: NO_SPACE = 0
    integer(kind=int1) :: SPACE = 1
                                                                                
    integer(kind=int1) :: NO = 0
    integer(kind=int1) :: YES = 1

    integer(kind=int1) :: NO_AUX_CLOUD_MASK = 0
    integer(kind=int1) :: USE_AUX_CLOUD_MASK = 1
    integer(kind=int1) :: READ_BUT_DO_NOT_USE_AUX_CLOUD_MASK = 2
                                                                                
    !--- this apply to the sfc_type array
    integer(kind=int1) :: WATER_SFC = 0
    integer(kind=int1) :: EVERGREEN_NEEDLE_SFC = 1
    integer(kind=int1) :: EVERGREEN_BROAD_SFC = 2
    integer(kind=int1) :: DECIDUOUS_NEEDLE_SFC = 3
    integer(kind=int1) :: DECIDUOUS_BROAD_SFC = 4
    integer(kind=int1) :: MIXED_FORESTS_SFC = 5
    integer(kind=int1) :: WOODLANDS_SFC = 6
    integer(kind=int1) :: WOODED_GRASS_SFC = 7
    integer(kind=int1) :: CLOSED_SHRUBS_SFC = 8
    integer(kind=int1) :: OPEN_SHRUBS_SFC = 9
    integer(kind=int1) :: GRASSES_SFC = 10
    integer(kind=int1) :: CROPLANDS_SFC = 11
    integer(kind=int1) :: BARE_SFC = 12
    integer(kind=int1) :: URBAN_SFC = 13
                                                                                
    integer(kind=int1) :: NO_DESERT = 0
    integer(kind=int1) :: NIR_DESERT = 1
    integer(kind=int1) :: BRIGHT_DESERT = 2

    !--- this apply to the land flags (land array)
    integer(kind=int1) :: SHALLOW_OCEAN = 0
    integer(kind=int1) :: LAND = 1
    integer(kind=int1) :: COASTLINE = 2
    integer(kind=int1) :: SHALLOW_INLAND_WATER = 3
    integer(kind=int1) :: EPHEMERAL_WATER = 4
    integer(kind=int1) :: DEEP_INLAND_WATER = 5
    integer(kind=int1) :: MODERATE_OCEAN = 6
    integer(kind=int1) :: DEEP_OCEAN = 7

    integer(kind=int1) :: NO_VOLCANO = 0
    integer(kind=int1) :: CLOSE_VOLCANO = 1
    integer(kind=int1) :: VERY_CLOSE_VOLCANO = 2

    integer(kind=int1) :: NO_COAST = 0
    integer(kind=int1) :: COAST_1KM = 1
    integer(kind=int1) :: COAST_2KM = 2
    integer(kind=int1) :: COAST_3KM = 3
    integer(kind=int1) :: COAST_4KM = 4
    integer(kind=int1) :: COAST_5KM = 5
    integer(kind=int1) :: COAST_6KM = 6
    integer(kind=int1) :: COAST_7KM = 7
    integer(kind=int1) :: COAST_8KM = 8
    integer(kind=int1) :: COAST_9KM = 9
    integer(kind=int1) :: COAST_10KM = 10 

    integer(kind=int1) :: WATER_GEN = 0
    integer(kind=int1) :: COAST_GEN = 1
    integer(kind=int1) :: LAND_GEN  = 2
    integer(kind=int1) :: DESERT_GEN = 3
    integer(kind=int1) :: SNOW_GEN   = 4

    integer(kind=int1) :: NO_SNOW = 1
    integer(kind=int1) :: SEA_ICE = 2
    integer(kind=int1) :: SNOW = 3

    integer(kind=int1) :: READ_SNOW_NWP = 0
    integer(kind=int1) :: READ_SNOW_HIRES = 1
    integer(kind=int1) :: READ_SNOW_GLOB = 2

    integer(kind=int1) :: SUCCESS = 0
    integer(kind=int1) :: FAILURE = 1
    integer(kind=int1) :: INFORMATION = 2
    integer(kind=int1) :: WARNING = 3
    integer(kind=int1) :: EOF = 4
    integer(kind=int1) :: UNDEFINED = 5
    integer(kind=int1) :: EXISTS =  6
    integer(kind=int1) :: EXIT = 7

    integer(kind=int1) :: SNOW_NOT_AVAILABLE = 1
    integer(kind=int1) :: NWP_SNOW = 2
    integer(kind=int1) :: IMS_SNOW = 3

    integer(kind=int1) :: CONSTANT_SFC_EMISS = 1
    integer(kind=int1) :: TABLE_SFC_EMISS    = 2
    integer(kind=int1) :: SEEBOR_SFC_EMISS   =3

    integer(kind=int1) :: CONSTANT_SFC_ALB   = 1
    integer(kind=int1) :: TABLE_SFC_ALB      = 2
    integer(kind=int1) :: MODIS_SFC_ALB      = 3

    integer(kind=int4) :: LITTLE_ENDIAN      = 0
    integer(kind=int4) :: BIG_ENDIAN         = 1
    integer(kind=int4) :: SIGNED             = 1
    integer(kind=int4) :: UNSIGNED           = 0
    integer(kind=int4) :: SWAP               = 1
    integer(kind=int4) :: NOSWAP             = 0

    !--- scaling options
    integer(kind=int1) :: NO_SCALING = 0
    integer(kind=int1) :: LINEAR_SCALING = 1
    integer(kind=int1) :: LOG10_SCALING = 2
    integer(kind=int1) :: SQUARE_ROOT_SCALING = 3

    !General product quality flags - VOLCAT
    integer(kind=int4) :: INVALID_PRODUCT = 0
    integer(kind=int4) :: UNKNOWN_PRODUCT_QUALITY = 1
    integer(kind=int4) :: LOW_PRODUCT_QUALITY = 2
    integer(kind=int4) :: HIGH_PRODUCT_QUALITY_NOintERP = 3
    integer(kind=int4) :: HIGH_PRODUCT_QUALITY_intERP = 4
    integer(kind=int4) :: INVALID_FILLED = 5
    integer(kind=int4) :: INVALID_REPLACED = 6

  END TYPE symbol_struct

!ccm  character(len=4), parameter:: SOLAR_OBS_TYPE = "SOLAR"
  character(len=5), parameter:: SOLAR_OBS_TYPE = "SOLAR"
  character(len=5), parameter:: LUNAR_OBS_TYPE = "LUNAR"
  character(len=5), parameter:: MIXED_OBS_TYPE = "MIXED"
  character(len=5), parameter:: THERMAL_OBS_TYPE = "THERM"
     

  TYPE(symbol_struct), public, save :: sym
  
end module CONSTANTS
