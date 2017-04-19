!$Id: acha_clavrx_services_module.f90 1789 2016-09-28 22:20:51Z heidinger $
!------------------------------------------------------------------------------
!this module holds all the dependencies for CCL for the various frameworks
!------------------------------------------------------------------------------
module CCL_SERVICES_MOD

 use PLANCK
 use CX_CONSTANTS_MOD
 use NWP_COMMON
 use RTM_COMMON
 use PIXEL_COMMON, only: &
       Ch, &
       Nav, &
       Sensor, &
       Image, &
       ACHA, &
       CCL, &
       Sfc, &
       CLDMASK, &
       Bad_Pixel_Mask, &
       I_Lrc, &
       J_Lrc, &
       I_Nwp, &
       J_Nwp, &
       I_Nwp_x, &
       J_Nwp_x, &
       Lon_Nwp_Fac, &
       Lat_Nwp_Fac, &
       Zen_Idx_Rtm, &
       Smooth_Nwp_Flag, &
       Cld_Type, &
       Diag_Pix_Array_1, &
       Diag_Pix_Array_2, &
       Diag_Pix_Array_3
 
 implicit none

 public:: CCL_FETCH_PIXEL_NWP_RTM 

 integer(KIND=INT4), PRIVATE, PARAMETER :: Num_Levels_Rtm_Prof = 101

!ACHA input structure
! input structure

 type, public :: ccl_diag_struct
  real (kind=real4), dimension(:,:), pointer:: Array_1
  real (kind=real4), dimension(:,:), pointer:: Array_2
  real (kind=real4), dimension(:,:), pointer:: Array_3
 end type ccl_diag_struct


 type, public :: ccl_input_struct
 integer (kind=int4):: Number_of_Elements
 integer (kind=int4):: Number_Of_Lines
 real (kind=real4):: Sensor_Resolution_KM

 logical, dimension(:,:), pointer:: Invalid_Data_Mask
 real, dimension(:,:), pointer:: Latitude
 real, dimension(:,:), pointer:: Longitude
 integer (kind=int1),dimension(:,:), pointer:: Surface_Type
 integer (kind=int1),dimension(:,:), pointer:: Cloud_Mask
 real, dimension(:,:), pointer:: Cloud_Probability
 integer (kind=int1),dimension(:,:), pointer:: Cloud_Type
 real, dimension(:,:), pointer:: Pc
 real, dimension(:,:), pointer:: Pc_Base

 integer (kind=int4), dimension(:,:), pointer:: Elem_Idx_NWP
 integer (kind=int4), dimension(:,:), pointer:: Line_Idx_NWP
 integer (kind=int4), dimension(:,:), pointer:: Elem_Idx_Opposite_Corner_NWP
 integer (kind=int4), dimension(:,:), pointer:: Line_Idx_Opposite_Corner_NWP
 integer (kind=int4), dimension(:,:), pointer:: Viewing_Zenith_Angle_Idx_Rtm
 real (kind=real4), dimension(:,:), pointer:: Latitude_Interp_Weight_NWP
 real (kind=real4), dimension(:,:), pointer:: Longitude_Interp_Weight_NWP
 
 end type ccl_input_struct

 !---RTM and NWP pixel level structure
 type, public :: ccl_rtm_nwp_struct

   !-- Smooth NWP Fields flag
   integer:: Smooth_Nwp_Fields_Flag_Temp
   
   !-- NWP Levels
   integer:: Sfc_Level
   integer:: Tropo_Level

   !-- NWP profiles
   real, dimension(:), pointer :: T_Prof
   real, dimension(Num_Levels_Rtm_Prof) :: P_Prof
   real, dimension(:), pointer :: Z_Prof

   !-- NWP profiles used for spatial interpolation
   real, dimension(:), pointer :: T_Prof_1
   real, dimension(:), pointer :: T_Prof_2
   real, dimension(:), pointer :: T_Prof_3
   
   real, dimension(:), pointer :: Z_Prof_1
   real, dimension(:), pointer :: Z_Prof_2
   real, dimension(:), pointer :: Z_Prof_3
end type ccl_rtm_nwp_struct

!output structure
 type, public :: ccl_output_struct
   real, dimension(:,:), pointer:: Total_Cloud_Fraction
   real, dimension(:,:), pointer:: Total_Cloud_Fraction_Uncer
   real, dimension(:,:), pointer:: High_Cloud_Fraction
   real, dimension(:,:), pointer:: Mid_Cloud_Fraction
   real, dimension(:,:), pointer:: Low_Cloud_Fraction
   integer (kind=int1), dimension(:,:), pointer:: Cloud_Layer
  end type ccl_output_struct
  
!Symbol stucture

 type, public :: ccl_symbol_struct
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
 end type ccl_symbol_struct
 
 contains

!----------------------------------------------------------------------
! This subroutine gathers the necessary NWP and RTM profiles used for a given
! pixel for ACHA. 
!----------------------------------------------------------------------
 subroutine  CCL_FETCH_PIXEL_NWP_RTM(CCL_Input, symbol, &
                                     Elem_Idx, Line_Idx, CCL_RTM_NWP)
                                      
   type(ccl_input_struct), intent(inout) :: CCL_Input
   type(ccl_rtm_nwp_struct), intent(inout) :: CCL_RTM_NWP
   type(ccl_symbol_struct), intent(inout) :: symbol
   integer, intent(in) :: Elem_Idx
   integer, intent(in) :: Line_Idx
   integer:: Ivza
   integer:: Inwp
   integer:: Jnwp
   integer:: Inwp_x
   integer:: Jnwp_x
   real:: Inwp_Weight
   real:: Jnwp_Weight

   Inwp = CCL_Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx)
   Jnwp = CCL_Input%Line_Idx_Nwp(Elem_Idx,Line_Idx)
   
   Inwp_x = CCL_Input%Elem_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
   Jnwp_x = CCL_Input%Line_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
   
   Inwp_Weight = CCL_Input%Longitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
   Jnwp_Weight = CCL_Input%Latitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
   Ivza =  CCL_Input%Viewing_Zenith_Angle_Idx_Rtm(Elem_Idx,Line_Idx)


   !--- populate height and temperature profiles
   if (Inwp <= 0 .or. Jnwp <= 0) then
     print *, "bad nwp indices in awg"
   endif
   if (Allocated(Rtm(Inwp,Jnwp)%T_Prof) .eqv. .false.) then
      print *, "error, T_Prof not allocated"
   endif

   !initialize smooth NWP flag 
   CCL_Rtm_NWP%Smooth_Nwp_Fields_Flag_Temp = symbol%NO
    
   CCL_Rtm_NWP%Sfc_Level = Rtm(Inwp,Jnwp)%Sfc_Level
   CCL_Rtm_NWP%Tropo_Level = Rtm(Inwp,Jnwp)%Tropo_Level
   
   
   CCL_Rtm_NWP%Smooth_Nwp_Fields_Flag_Temp = symbol%NO
   
   !--- do various 101 level NWP Profiles
   CCL_Rtm_NWP%P_Prof = P_Std_Rtm

   CCL_Rtm_NWP%T_Prof => Rtm(Inwp,Jnwp)%T_Prof 
   CCL_Rtm_NWP%Z_Prof => Rtm(Inwp,Jnwp)%Z_Prof 

   !------------------------------------------------------
   ! Before smoothing profiles, ensure that all required
   ! rtm profiles are populated, if not, skip smoothing
   !------------------------------------------------------
   if ((Rtm(Inwp,Jnwp)%Flag == symbol%YES) .and. &
       (Rtm(Inwp_x,Jnwp)%Flag == symbol%YES) .and. &
       (Rtm(Inwp,Jnwp_x)%Flag == symbol%YES) .and. &
       (Rtm(Inwp_x,Jnwp_x)%Flag == symbol%YES)) then

        CCL_Rtm_NWP%Smooth_Nwp_Fields_Flag_Temp = symbol%YES
        
        CCL_Rtm_NWP%T_Prof_1 => Rtm(Inwp_x,Jnwp)%T_Prof 
        CCL_Rtm_NWP%T_Prof_2 => Rtm(Inwp,Jnwp_x)%T_Prof 
        CCL_Rtm_NWP%T_Prof_3 => Rtm(Inwp_x,Jnwp_x)%T_Prof 

        CCL_Rtm_NWP%Z_Prof_1 => Rtm(Inwp_x,Jnwp)%Z_Prof 
        CCL_Rtm_NWP%Z_Prof_2 => Rtm(Inwp,Jnwp_x)%Z_Prof 
        CCL_Rtm_NWP%Z_Prof_3 => Rtm(Inwp_x,Jnwp_x)%Z_Prof
        
   endif
   
 end subroutine CCL_FETCH_PIXEL_NWP_RTM


end module CCL_SERVICES_MOD
