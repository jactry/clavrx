!$Id: acha_clavrx_bridge_module.f90 1789 2016-09-28 22:20:51Z heidinger $
!------------------------------------------------------------------------------
!  NOAA AWG Cloud Height Algorithm (ACHA) Bridge Code
!
!  This module houses the routines that serve as a bridge between
!  the CLAVR-x processing system and the ACHA code.
!
!------------------------------------------------------------------------------
module ASOS_CLAVRX_BRIDGE

 use ASOS_SERVICES_MOD
 use ASOS_MODULE
 use CLAVRX_MESSAGE_MODULE, only: MESG
   
 implicit none

 public :: ASOS_BRIDGE
 private:: SET_SYMBOL, SET_INPUT, SET_OUTPUT, NULL_INPUT, NULL_OUTPUT
 private:: SET_DIAG, NULL_DIAG

 !--------------------------------------------------------------------
 ! define structures that will be arguments to ASOS
 !--------------------------------------------------------------------
 type(asos_symbol_struct), private :: Symbol
 type(asos_input_struct), private :: Input
 type(asos_output_struct), private :: Output
 type(asos_diag_struct), private :: Diag

 contains

!----------------------------------------------------------------------
! ASOS BRIDGE SUBROUTINE
!---------------------------------------------------------------------- 
 subroutine ASOS_BRIDGE()
 
   implicit none

   logical, save:: First_Call = .true.

   if (First_Call .eqv. .true.) then
       call MESG('ASOS starts ', color = 46)
   endif

   !---null pointers before filling them
   call NULL_INPUT()
   call NULL_OUTPUT()

   !-------------------------------------------
   !--- initialize structures
   !-------------------------------------------

   !--- store integer values
   call SET_INPUT()

   !---- initalize Output structure
   call SET_OUTPUT()
  
   !----set symbols to local values
   call SET_SYMBOL()

   !---- initialize diagnostic structure
   call SET_DIAG()

   !--- cloud cover layers
   call COMPUTE_ASOS(Input, Symbol, Output) !, Diag)

   !--- copy output into CLAVR-x variables
   ASOS%Code = Output%ASOS_Cloud_Code
   ASOS%ECA = Output%ASOS_Cloud_ECA
   ASOS%Zmin = Output%ASOS_Cloud_Zmin
   ASOS%Zmax = Output%ASOS_Cloud_Zmax

   !-----------------------------------------------------------------------
   !--- Null pointers after algorithm is finished
   !-----------------------------------------------------------------------
   call NULL_INPUT()
   call NULL_OUTPUT()

   First_Call = .false.

 end subroutine ASOS_BRIDGE

 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding input 
 !-----------------------------------------------------------------------------
 subroutine NULL_INPUT()
     Input%Invalid_Data_Mask =>  null()
     Input%Latitude =>  null()
     Input%Longitude =>  null()
     Input%Surface_Type =>  null()
     Input%Cloud_Mask =>  null()
     Input%Cloud_Probability => null()
     Input%Cloud_Type =>  null()
     Input%Pc =>  null()
     Input%Zc =>  null()
     Input%Ec =>  null()
 end subroutine NULL_INPUT
 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding input 
 !-----------------------------------------------------------------------------
 subroutine NULL_DIAG()
     Diag%Array_1 =>  null()
     Diag%Array_2 =>  null()
     Diag%Array_3 =>  null()
 end subroutine NULL_DIAG
 !-----------------------------------------------------------------------------
 ! Nullify the pointers holding output
 !-----------------------------------------------------------------------------
 subroutine NULL_OUTPUT()
   Output%ASOS_Cloud_Code => null()
   Output%ASOS_Cloud_ECA => null()
   Output%ASOS_Cloud_Zmin => null()
   Output%ASOS_Cloud_Zmax => null()
 end subroutine NULL_OUTPUT
 !-----------------------------------------------------------------------------
 ! Copy needed Symbol elements
 !-----------------------------------------------------------------------------
 subroutine SET_SYMBOL()
   Symbol%CLOUDY = sym%CLOUDY
   Symbol%PROB_CLOUDY = sym%PROB_CLOUDY
   Symbol%PROB_CLEAR = sym%PROB_CLEAR
   Symbol%CLEAR = sym%CLEAR

   Symbol%NO = sym%NO
   Symbol%YES = sym%YES

   Symbol%WATER_SFC = sym%WATER_SFC
   Symbol%EVERGREEN_NEEDLE_SFC = sym%EVERGREEN_NEEDLE_SFC
   Symbol%EVERGREEN_BROAD_SFC = sym%EVERGREEN_BROAD_SFC
   Symbol%DECIDUOUS_NEEDLE_SFC = sym%DECIDUOUS_NEEDLE_SFC
   Symbol%DECIDUOUS_BROAD_SFC = sym%DECIDUOUS_BROAD_SFC
   Symbol%MIXED_FORESTS_SFC = sym%MIXED_FORESTS_SFC
   Symbol%WOODLANDS_SFC = sym%WOODLANDS_SFC
   Symbol%WOODED_GRASS_SFC = sym%WOODED_GRASS_SFC
   Symbol%CLOSED_SHRUBS_SFC = sym%CLOSED_SHRUBS_SFC
   Symbol%OPEN_SHRUBS_SFC = sym%OPEN_SHRUBS_SFC
   Symbol%GRASSES_SFC = sym%GRASSES_SFC
   Symbol%CROPLANDS_SFC = sym%CROPLANDS_SFC
   Symbol%BARE_SFC = sym%BARE_SFC
   Symbol%URBAN_SFC = sym%URBAN_SFC

   Symbol%SHALLOW_OCEAN = sym%SHALLOW_OCEAN
   Symbol%LAND = sym%LAND
   Symbol%COASTLINE = sym%COASTLINE
   Symbol%SHALLOW_INLAND_WATER = sym%SHALLOW_INLAND_WATER
   Symbol%EPHEMERAL_WATER = sym%EPHEMERAL_WATER
   Symbol%DEEP_INLAND_WATER = sym%DEEP_INLAND_WATER
   Symbol%MODERATE_OCEAN = sym%MODERATE_OCEAN
   Symbol%DEEP_OCEAN = sym%DEEP_OCEAN

   Symbol%NO_SNOW = sym%NO_SNOW
   Symbol%SEA_ICE = sym%SEA_ICE
   Symbol%SNOW = sym%SNOW

   Symbol%CLEAR_TYPE = sym%CLEAR_TYPE
   Symbol%PROB_CLEAR_TYPE = sym%PROB_CLEAR_TYPE
   Symbol%FOG_TYPE = sym%FOG_TYPE
   Symbol%WATER_TYPE = sym%WATER_TYPE
   Symbol%SUPERCOOLED_TYPE = sym%SUPERCOOLED_TYPE
   Symbol%MIXED_TYPE = sym%MIXED_TYPE
   Symbol%OPAQUE_ICE_TYPE = sym%OPAQUE_ICE_TYPE
   Symbol%TICE_TYPE = sym%TICE_TYPE
   Symbol%CIRRUS_TYPE = sym%CIRRUS_TYPE
   Symbol%OVERLAP_TYPE = sym%OVERLAP_TYPE
   Symbol%OVERSHOOTING_TYPE = sym%OVERSHOOTING_TYPE
   Symbol%UNKNOWN_TYPE = sym%UNKNOWN_TYPE
   Symbol%DUST_TYPE = sym%DUST_TYPE
   Symbol%SMOKE_TYPE = sym%SMOKE_TYPE
   Symbol%FIRE_TYPE = sym%FIRE_TYPE

   Symbol%CLEAR_PHASE = sym%CLEAR_PHASE
   Symbol%WATER_PHASE = sym%WATER_PHASE
   Symbol%SUPERCOOLED_PHASE = sym%SUPERCOOLED_PHASE
   Symbol%MIXED_PHASE = sym%MIXED_PHASE
   Symbol%ICE_PHASE = sym%ICE_PHASE
   Symbol%UNKNOWN_PHASE = sym%UNKNOWN_PHASE
 end subroutine SET_SYMBOL

 subroutine SET_OUTPUT()
   Output%ASOS_Cloud_Code => ASOS%Code
   Output%ASOS_Cloud_ECA => ASOS%ECA
   Output%ASOS_Cloud_Zmin => ASOS%Zmin
   Output%ASOS_Cloud_Zmax => ASOS%Zmax
 end subroutine SET_OUTPUT
!--------------------------------------------------------
 subroutine SET_INPUT()

   Input%Number_of_Elements = Image%Number_Of_Elements
   Input%Number_of_Lines = Image%Number_Of_Lines_Per_Segment
   Input%Invalid_Data_Mask =>  Bad_Pixel_Mask
   Input%Latitude => Nav%Lat
   Input%Longitude => Nav%Lon
   Input%Surface_Type => Sfc%Sfc_Type
   Input%Cloud_Mask => Cld_Mask
   Input%Cloud_Probability => Posterior_Cld_Probability
   Input%Cloud_Type => Cld_Type
   Input%Pc => ACHA%Pc
   Input%Zc => ACHA%Zc
   Input%Ec => ACHA%Ec
 end subroutine SET_INPUT
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
 subroutine SET_DIAG
     Diag%Array_1 => Diag_Pix_Array_1 
     Diag%Array_2 => Diag_Pix_Array_2 
     Diag%Array_3 => Diag_Pix_Array_3 
 end subroutine SET_DIAG

end module ASOS_CLAVRX_BRIDGE
