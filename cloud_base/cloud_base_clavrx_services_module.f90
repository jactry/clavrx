!$Id: acha_services_clavrx_mod.f90 581 2014-10-08 03:39:08Z heidinger $
!------------------------------------------------------------------------------
!this module holds all the dependencies for ACHA for the various frameworks
!------------------------------------------------------------------------------
module CLOUD_BASE_SERVICES

 use PLANCK
 use CX_CONSTANTS_MOD
 use PIXEL_COMMON
 use NWP_COMMON
 use RTM_COMMON
 use NUMERICAL_TOOLS_MOD, only: LOCATE

 implicit none

 public:: FETCH_PIXEL_RTM_NWP 

 integer(KIND=INT4), PRIVATE, PARAMETER :: Num_Levels_Rtm_Prof = 101

!ACHA input structure
! input structure

 type, public :: acha_input_struct
 integer :: ACHA_Mode_Flag_In
 integer (kind=int4):: Number_of_Elements
 integer (kind=int4):: Number_Of_Lines
 integer (kind=int4):: Num_Line_Max

 !-- local pointers that point to global variables
 integer:: Chan_Idx_67um
 integer:: Chan_Idx_85um
 integer:: Chan_Idx_11um
 integer:: Chan_Idx_12um
 integer:: Chan_Idx_133um 
 integer:: Chan_On_67um
 integer:: Chan_On_85um
 integer:: Chan_On_11um
 integer:: Chan_On_12um
 integer:: Chan_On_133um

 integer (kind=int1), dimension(:,:), pointer:: Invalid_Data_Mask
 real, dimension(:,:), pointer:: Cosine_Zenith_Angle
 real, dimension(:,:), pointer:: Sensor_Zenith_Angle
 real, dimension(:,:), pointer:: Surface_Elevation
 real, dimension(:,:), pointer:: Latitude
 real, dimension(:,:), pointer:: Longitude
 integer (kind=int1),dimension(:,:), pointer:: Snow_Class
 integer (kind=int1),dimension(:,:), pointer:: Surface_Type
 integer (kind=int1),dimension(:,:), pointer:: Cloud_Mask
 real, dimension(:,:), pointer:: Cloud_Probability
 integer (kind=int1),dimension(:,:), pointer:: Cloud_Type
 integer (kind=int4), dimension(:,:), pointer:: Elem_Idx_NWP 
 integer (kind=int4), dimension(:,:), pointer:: Line_Idx_NWP 
 integer (kind=int4), dimension(:,:), pointer:: Elem_Idx_Opposite_Corner_NWP 
 integer (kind=int4), dimension(:,:), pointer:: Line_Idx_Opposite_Corner_NWP 
 integer (kind=int4), dimension(:,:), pointer:: Viewing_Zenith_Angle_Idx_Rtm
 real (kind=real4), dimension(:,:), pointer:: Latitude_Interp_Weight_NWP
 real (kind=real4), dimension(:,:), pointer:: Longitude_Interp_Weight_NWP
 integer(kind=int4), dimension(:,:), pointer :: Elem_Idx_LRC_Input
 integer(kind=int4), dimension(:,:), pointer :: Line_Idx_LRC_Input
 real, dimension(:,:), pointer:: Latitude_Pc                         !parallax corrected lat
 real, dimension(:,:), pointer:: Longitude_Pc                        !parallax corrected lon
 real, dimension(:,:), pointer:: Tc                                  !cloud temperature
 real, dimension(:,:), pointer:: Ec                                  !cloud emissivity
 real, dimension(:,:), pointer:: Beta                                !cloud beta
 real, dimension(:,:), pointer:: Pc                                  !cloud pressure
 real, dimension(:,:), pointer:: Zc                                  !cloud height
 real, dimension(:,:), pointer:: Tau                                 !cloud optical depth
 real, dimension(:,:), pointer:: Reff                                !cloud effective particle size
 real, dimension(:,:), pointer:: Tc_Uncertainty                      !uncertainty in cloud temperature
 real, dimension(:,:), pointer:: Ec_Uncertainty                      !uncertainty in cloud emissivity
 real, dimension(:,:), pointer:: Beta_Uncertainty                    !uncertainty in cloud beta
 real, dimension(:,:), pointer:: Pc_Uncertainty                      !uncertainty in cloud pressure
 real, dimension(:,:), pointer:: Zc_Uncertainty                      !uncertainty in cloud height
 real, dimension(:,:), pointer:: Lower_Cloud_Pressure                !pressure of lower cloud layer (if present)
 real, dimension(:,:), pointer:: Lower_Cloud_Temperature             !temperature of lower cloud layer (if present)
 real, dimension(:,:), pointer:: Lower_Cloud_Height                  !height of lower cloud layer (if present)
 real, dimension(:,:), pointer:: Cdnc                                !cloud droplet number concentration
 real, dimension(:,:), pointer:: Hcld                                !cloud geometrical thickness
 real, dimension(:,:), pointer:: LCL                                 !lifting condensation level
 real, dimension(:,:), pointer:: CCL                                 !convective condensation level
 real, dimension(:,:), pointer:: CWP                                 !cloud water path
 real, dimension(:,:), pointer:: CWP_nwp

 end type acha_input_struct

 !---RTM and NWP pixel level structure
 type, public :: acha_rtm_nwp_struct

   !-- Smooth NWP Fields flag
   integer:: Smooth_Nwp_Fields_Flag_Temp
   
   !-- NWP Levels
   integer:: Sfc_Level
   integer:: Tropo_Level

   !-- RTM profiles
   real, dimension(:), pointer :: Atm_Rad_Prof_67um
   real, dimension(:), pointer :: Atm_Rad_Prof_85um
   real, dimension(:), pointer :: Atm_Rad_Prof_11um
   real, dimension(:), pointer :: Atm_Rad_Prof_12um
   real, dimension(:), pointer :: Atm_Rad_Prof_133um
   real, dimension(:), pointer :: Atm_Trans_Prof_67um
   real, dimension(:), pointer :: Atm_Trans_Prof_85um
   real, dimension(:), pointer :: Atm_Trans_Prof_11um
   real, dimension(:), pointer :: Atm_Trans_Prof_12um
   real, dimension(:), pointer :: Atm_Trans_Prof_133um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_67um
   real, dimension(:), pointer :: Black_Body_Rad_Prof_11um

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
end type acha_rtm_nwp_struct

!output structure
 type, public :: acha_output_struct
   real, dimension(:,:), pointer:: Zc_Top
   real, dimension(:,:), pointer:: Zc_Base
   real, dimension(:,:), pointer:: Pc_Top
   real, dimension(:,:), pointer:: Pc_Base
   integer (kind=int1), dimension(:,:), pointer:: Zc_Base_Qf
 end type acha_output_struct

!Symbol stucture

 type, public :: symbol_acha
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
 end type symbol_acha
 
 contains

!----------------------------------------------------------------------
! This subroutine gathers the necessary NWP and RTM profiles used for a given
! pixel for ACHA. 
!----------------------------------------------------------------------
 subroutine  FETCH_PIXEL_RTM_NWP(Input, symbol, &
                                 Elem_Idx, Line_Idx, RTM_NWP)
                                      
   type(acha_input_struct), intent(inout) :: Input
   type(acha_rtm_nwp_struct), intent(inout) :: RTM_NWP
   type(symbol_acha), intent(inout) :: symbol
   integer, intent(in) :: Elem_Idx
   integer, intent(in) :: Line_Idx
   integer:: Ivza
   integer:: Inwp
   integer:: Jnwp
   integer:: Inwp_x
   integer:: Jnwp_x
   real:: Inwp_Weight
   real:: Jnwp_Weight

   Inwp = Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx)
   Jnwp = Input%Line_Idx_Nwp(Elem_Idx,Line_Idx)
   
   Inwp_x = Input%Elem_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
   Jnwp_x = Input%Line_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
   
   Inwp_Weight = Input%Longitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
   Jnwp_Weight = Input%Latitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
   Ivza =  Input%Viewing_Zenith_Angle_Idx_Rtm(Elem_Idx,Line_Idx)

   !--- populate height and temperature profiles
   if (Inwp <= 0 .or. Jnwp <= 0) then
     print *, "bad nwp indices in cloud base"
     return
   endif
   if (Allocated(Rtm(Inwp,Jnwp)%T_Prof) .eqv. .false.) then
      print *, "error, T_Prof not allocated"
   endif

   !initialize smooth NWP flag 
   RTM_NWP%Smooth_Nwp_Fields_Flag_Temp = symbol%NO
    
   RTM_NWP%Sfc_Level = Rtm(Inwp,Jnwp)%Sfc_Level
   RTM_NWP%Tropo_Level = Rtm(Inwp,Jnwp)%Tropo_Level
   
   RTM_NWP%Smooth_Nwp_Fields_Flag_Temp = symbol%NO
   
   !--- do various 101 level NWP Profiles
   RTM_NWP%P_Prof = P_Std_Rtm

   RTM_NWP%T_Prof => Rtm(Inwp,Jnwp)%T_Prof 
   RTM_NWP%Z_Prof => Rtm(Inwp,Jnwp)%Z_Prof 

   !------------------------------------------------------
   ! Before smoothing profiles, ensure that all required
   ! rtm profiles are populated, if not, skip smoothing
   !------------------------------------------------------
   if ((Rtm(Inwp,Jnwp)%Flag == symbol%YES) .and. &
       (Rtm(Inwp_x,Jnwp)%Flag == symbol%YES) .and. &
       (Rtm(Inwp,Jnwp_x)%Flag == symbol%YES) .and. &
       (Rtm(Inwp_x,Jnwp_x)%Flag == symbol%YES)) then

        RTM_NWP%Smooth_Nwp_Fields_Flag_Temp = symbol%YES
        
        RTM_NWP%T_Prof_1 => Rtm(Inwp_x,Jnwp)%T_Prof 
        RTM_NWP%T_Prof_2 => Rtm(Inwp,Jnwp_x)%T_Prof 
        RTM_NWP%T_Prof_3 => Rtm(Inwp_x,Jnwp_x)%T_Prof 

        RTM_NWP%Z_Prof_1 => Rtm(Inwp_x,Jnwp)%Z_Prof 
        RTM_NWP%Z_Prof_2 => Rtm(Inwp,Jnwp_x)%Z_Prof 
        RTM_NWP%Z_Prof_3 => Rtm(Inwp_x,Jnwp_x)%Z_Prof
        
   endif
   
   !---- RTM profiles
 
   !--- populate radiance and transmission profiles
   if (Input%Chan_On_67um == sym%YES) then
     RTM_NWP%Atm_Rad_Prof_67um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Rad_Atm_Profile
     RTM_NWP%Atm_Trans_Prof_67um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Trans_Atm_Profile
     RTM_NWP%Black_Body_Rad_Prof_67um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Rad_BB_Cloud_Profile
   endif

   if (Input%Chan_On_85um == sym%YES) then
     RTM_NWP%Atm_Rad_Prof_85um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(29)%Rad_Atm_Profile
     RTM_NWP%Atm_Trans_Prof_85um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(29)%Trans_Atm_Profile
   endif

   if (Input%Chan_On_11um == sym%YES) then
      RTM_NWP%Atm_Rad_Prof_11um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Rad_Atm_Profile
      RTM_NWP%Atm_Trans_Prof_11um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Trans_Atm_Profile
      RTM_NWP%Black_Body_Rad_Prof_11um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Rad_BB_Cloud_Profile
   endif
   
   if (Input%Chan_On_12um == sym%YES) then
      RTM_NWP%Atm_Rad_Prof_12um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(32)%Rad_Atm_Profile
      RTM_NWP%Atm_Trans_Prof_12um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(32)%Trans_Atm_Profile
   endif

   if (Input%Chan_On_133um == sym%YES) then
      RTM_NWP%Atm_Rad_Prof_133um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(33)%Rad_Atm_Profile
      RTM_NWP%Atm_Trans_Prof_133um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(33)%Trans_Atm_Profile
   endif
    
 end subroutine FETCH_PIXEL_RTM_NWP


end module CLOUD_BASE_SERVICES
