!$Id: acha_clavrx_services_module.f90 1789 2016-09-28 22:20:51Z heidinger $
!------------------------------------------------------------------------------
!this module holds all the dependencies for ACHA for the various frameworks
!------------------------------------------------------------------------------
module ASOS_SERVICES_MOD

 use PLANCK
 use CONSTANTS
 use NWP_COMMON
 use RTM_COMMON
 use PIXEL_COMMON, only: &
       Ch, &
       Nav, &
       Geo, &
       Sensor, &
       Image, &
       ACHA, &
       ASOS, &
       Sfc, &
       Bad_Pixel_Mask, &
       Cld_Mask, &
       Cld_Type, &
       Posterior_Cld_Probability, &
       Diag_Pix_Array_1, &
       Diag_Pix_Array_2, &
       Diag_Pix_Array_3
 
 implicit none

 type, public :: asos_diag_struct
  real (kind=real4), dimension(:,:), pointer:: Array_1
  real (kind=real4), dimension(:,:), pointer:: Array_2
  real (kind=real4), dimension(:,:), pointer:: Array_3
 end type asos_diag_struct


 type, public :: asos_input_struct
 integer (kind=int4):: Number_of_Elements
 integer (kind=int4):: Number_Of_Lines
 real (kind=real4):: Sensor_Resolution_KM

 integer (kind=int1), dimension(:,:), pointer:: Invalid_Data_Mask
 real, dimension(:,:), pointer:: Latitude
 real, dimension(:,:), pointer:: Longitude
 integer (kind=int1),dimension(:,:), pointer:: Surface_Type
 integer (kind=int1),dimension(:,:), pointer:: Cloud_Mask
 real, dimension(:,:), pointer:: Cloud_Probability
 integer (kind=int1),dimension(:,:), pointer:: Cloud_Type
 real (kind=real4),dimension(:,:), pointer:: Pc
 real (kind=real4),dimension(:,:), pointer:: Zc
 real (kind=real4),dimension(:,:), pointer:: Ec
 
 end type asos_input_struct

!output structure
 type, public :: asos_output_struct
   integer(kind=int1), dimension(:,:), pointer:: ASOS_Cloud_Code
   real, dimension(:,:), pointer:: ASOS_Cloud_ECA
   real, dimension(:,:), pointer:: ASOS_Cloud_Zmin
   real, dimension(:,:), pointer:: ASOS_Cloud_Zmax
  end type asos_output_struct
  
!Symbol stucture

 type, public :: asos_symbol_struct
    integer(kind=int1) :: CLOUDY
    integer(kind=int1) :: PROB_CLOUDY
    integer(kind=int1) :: PROB_CLEAR
    integer(kind=int1) :: CLEAR

    integer(kind=int1) :: NO
    integer(kind=int1) :: YES

    integer(kind=int1) :: WATER_SFC
    integer(kind=int1) :: EVERGREEN_NEEDLE_SFC
    integer(kind=int1) :: EVERGREEN_BROAD_SFC
    integer(kind=int1) :: DECIDUOUS_NEEDLE_SFC
    integer(kind=int1) :: DECIDUOUS_BROAD_SFC
    integer(kind=int1) :: MIXED_FORESTS_SFC
    integer(kind=int1) :: WOODLANDS_SFC
    integer(kind=int1) :: WOODED_GRASS_SFC
    integer(kind=int1) :: CLOSED_SHRUBS_SFC
    integer(kind=int1) :: OPEN_SHRUBS_SFC
    integer(kind=int1) :: GRASSES_SFC
    integer(kind=int1) :: CROPLANDS_SFC
    integer(kind=int1) :: BARE_SFC
    integer(kind=int1) :: URBAN_SFC

    integer(kind=int1) :: SHALLOW_OCEAN
    integer(kind=int1) :: LAND
    integer(kind=int1) :: COASTLINE
    integer(kind=int1) :: SHALLOW_INLAND_WATER
    integer(kind=int1) :: EPHEMERAL_WATER
    integer(kind=int1) :: DEEP_INLAND_WATER
    integer(kind=int1) :: MODERATE_OCEAN
    integer(kind=int1) :: DEEP_OCEAN

    integer(kind=int1) :: NO_SNOW
    integer(kind=int1) :: SEA_ICE
    integer(kind=int1) :: SNOW

    integer(kind=int1) :: CLEAR_type
    integer(kind=int1) :: PROB_CLEAR_type
    integer(kind=int1) :: FOG_type
    integer(kind=int1) :: WATER_type
    integer(kind=int1) :: SUPERCOOLED_type
    integer(kind=int1) :: MIXED_type
    integer(kind=int1) :: OPAQUE_ICE_type
    integer(kind=int1) :: TICE_type
    integer(kind=int1) :: CIRRUS_type
    integer(kind=int1) :: OVERLAP_type
    integer(kind=int1) :: OVERSHOOTING_type
    integer(kind=int1) :: UNKNOWN_type
    integer(kind=int1) :: DUST_type
    integer(kind=int1) :: SMOKE_type
    integer(kind=int1) :: FIRE_type

    integer(kind=int1) :: CLEAR_PHASE
    integer(kind=int1) :: WATER_PHASE
    integer(kind=int1) :: SUPERCOOLED_PHASE
    integer(kind=int1) :: MIXED_PHASE
    integer(kind=int1) :: ICE_PHASE
    integer(kind=int1) :: UNKNOWN_PHASE
 end type asos_symbol_struct
 
 contains


end module ASOS_SERVICES_MOD
