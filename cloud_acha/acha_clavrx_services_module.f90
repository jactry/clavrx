!$Id$
!------------------------------------------------------------------------------
!this module holds all the dependencies for ACHA for the various frameworks
!------------------------------------------------------------------------------
module ACHA_SERVICES_MOD

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
       Sfc, &
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
       Tair_Nwp_Pix, &
       Tsfc_Nwp_Pix, &
       Ttropo_Nwp_Pix, &
       Psfc_Nwp_Pix, &
       Smooth_Nwp_Flag, &
       Process_Undetected_Cloud_Flag, &
       Posterior_Cld_Probability, &
       Cld_Mask, &
       Cld_Type, &
       Cld_Test_Vector_Packed, &
       Shadow_Mask, &
       Cloud_Fraction_Uncer_3x3, &
       Tc_Cirrus_Co2, &
       Tc_Opaque_Cloud, &
       Pc_Opaque_Cloud, &
       Zc_Opaque_Cloud, &
       Tc_H2O, &
       Pc_H2O, &
       Cloud_Fraction_3x3, &
       Low_Cloud_Fraction_3x3, &
       Mid_Cloud_Fraction_3x3, &
       High_Cloud_Fraction_3x3, &
       Diag_Pix_Array_1, &
       Diag_Pix_Array_2, &
       Diag_Pix_Array_3
		 
	use NUMERICAL_ROUTINES, only: &
		invert_matrix &
		, locate
 implicit none

 public:: ACHA_FETCH_PIXEL_NWP_RTM 

 integer(KIND=INT4), PRIVATE, PARAMETER :: Num_Levels_Rtm_Prof = 101

!ACHA input structure
! input structure

 type, public :: acha_input_struct
 integer :: ACHA_Mode_Flag_In
 integer (kind=int4):: Number_of_Elements
 integer (kind=int4):: Number_Of_Lines
 integer (kind=int4):: Num_Line_Max
 integer (kind=int4):: Smooth_Nwp_Fields_Flag
 integer (kind=int4):: Process_Undetected_Cloud_Flag
 real (kind=real4):: Sensor_Resolution_KM

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
 real, dimension(:,:), pointer:: Bt_67um
 real, dimension(:,:), pointer:: Bt_85um
 real, dimension(:,:), pointer:: Bt_11um
 real, dimension(:,:), pointer:: Bt_12um
 real, dimension(:,:), pointer:: Bt_133um
 real, dimension(:,:), pointer:: Rad_67um
 real, dimension(:,:), pointer:: Rad_11um
 real, dimension(:,:), pointer:: Cosine_Zenith_Angle
 real, dimension(:,:), pointer:: Sensor_Zenith_Angle
 real, dimension(:,:), pointer:: Sensor_Azimuth_Angle
 real, dimension(:,:), pointer:: Surface_Temperature
 real, dimension(:,:), pointer:: Surface_Air_Temperature
 real, dimension(:,:), pointer:: Tropopause_Temperature
 real, dimension(:,:), pointer:: Surface_Pressure
 real, dimension(:,:), pointer:: Surface_Elevation
 real, dimension(:,:), pointer:: Latitude
 real, dimension(:,:), pointer:: Longitude
 real, dimension(:,:), pointer:: Rad_Clear_67um
 real, dimension(:,:), pointer:: Rad_Clear_85um
 real, dimension(:,:), pointer:: Rad_Clear_11um
 real, dimension(:,:), pointer:: Rad_Clear_12um
 real, dimension(:,:), pointer:: Rad_Clear_133um
 real, dimension(:,:), pointer:: Surface_Emissivity_39um 
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
 real (kind=real4), dimension(:,:), pointer:: Tc_Cirrus_Sounder
 
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
   real, dimension(:,:), pointer:: Latitude_Pc
   real, dimension(:,:), pointer:: Longitude_Pc
   real, dimension(:,:), pointer:: Tc
   real, dimension(:,:), pointer:: Ec
   real, dimension(:,:), pointer:: Beta
   real, dimension(:,:), pointer:: Pc
   real, dimension(:,:), pointer:: Zc
   real, dimension(:,:), pointer:: Tau
   real, dimension(:,:), pointer:: Reff
   real, dimension(:,:), pointer:: Tc_Uncertainty
   real, dimension(:,:), pointer:: Ec_Uncertainty
   real, dimension(:,:), pointer:: Beta_Uncertainty
   real, dimension(:,:), pointer:: Pc_Uncertainty
   real, dimension(:,:), pointer:: Zc_Uncertainty
   real, dimension(:,:), pointer:: Lower_Cloud_Pressure
   real, dimension(:,:), pointer:: Lower_Cloud_Temperature
   real, dimension(:,:), pointer:: Lower_Cloud_Height
   real, dimension(:,:), pointer:: Cost
   real, dimension(:,:), pointer:: Total_Cloud_Fraction
   real, dimension(:,:), pointer:: Total_Cloud_Fraction_Uncer
   real, dimension(:,:), pointer:: High_Cloud_Fraction
   real, dimension(:,:), pointer:: Mid_Cloud_Fraction
   real, dimension(:,:), pointer:: Low_Cloud_Fraction
   integer (kind=int1), dimension(:,:), pointer:: Cloud_Layer
   integer (kind=int1), dimension(:,:), pointer:: Qf
   integer (kind=int1), dimension(:,:,:), pointer:: OE_Qf
   integer (kind=int1), dimension(:,:), pointer :: Packed_Qf
   integer (kind=int1), dimension(:,:), pointer :: Packed_Meta_Data
   integer(kind=int1), dimension(:,:), pointer :: Processing_Order   
   integer(kind=int1), dimension(:,:), pointer :: Inversion_Flag
   real, dimension(:,:), pointer:: Pc_Opaque
   real, dimension(:,:), pointer:: Tc_Opaque
   real, dimension(:,:), pointer:: Zc_Opaque
   real, dimension(:,:), pointer:: Pc_H2O
   real, dimension(:,:), pointer:: Tc_H2O
   real, dimension(:,:), pointer:: Zc_H2O
  end type acha_output_struct
  
!Symbol stucture

 type, public :: acha_symbol_struct
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
 end type acha_symbol_struct
 
 contains

!----------------------------------------------------------------------
! This subroutine gathers the necessary NWP and RTM profiles used for a given
! pixel for ACHA. 
!----------------------------------------------------------------------
 subroutine  ACHA_FETCH_PIXEL_NWP_RTM(Acha_Input, symbol, &
                                      Elem_Idx, Line_Idx, Acha_RTM_NWP)
                                      
   type(acha_input_struct), intent(inout) :: Acha_Input
   type(acha_rtm_nwp_struct), intent(inout) :: Acha_RTM_NWP
   type(acha_symbol_struct), intent(inout) :: symbol
   integer, intent(in) :: Elem_Idx
   integer, intent(in) :: Line_Idx
   integer:: Ivza
   integer:: Inwp
   integer:: Jnwp
   integer:: Inwp_x
   integer:: Jnwp_x
   real:: Inwp_Weight
   real:: Jnwp_Weight

   Inwp = Acha_Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx)
   Jnwp = Acha_Input%Line_Idx_Nwp(Elem_Idx,Line_Idx)
   
   Inwp_x = Acha_Input%Elem_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
   Jnwp_x = Acha_Input%Line_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
   
   Inwp_Weight = Acha_Input%Longitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
   Jnwp_Weight = Acha_Input%Latitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
   Ivza =  Acha_Input%Viewing_Zenith_Angle_Idx_Rtm(Elem_Idx,Line_Idx)


   !--- populate height and temperature profiles
   if (Inwp <= 0 .or. Jnwp <= 0) then
     print *, "bad nwp indices in awg"
   endif
   if (Allocated(Rtm(Inwp,Jnwp)%T_Prof) .eqv. .false.) then
      print *, "error, T_Prof not allocated"
   endif

   !initialize smooth NWP flag 
   Acha_RTM_NWP%Smooth_Nwp_Fields_Flag_Temp = symbol%NO
    
   Acha_RTM_NWP%Sfc_Level = Rtm(Inwp,Jnwp)%Sfc_Level
   Acha_RTM_NWP%Tropo_Level = Rtm(Inwp,Jnwp)%Tropo_Level
   
   
   Acha_RTM_NWP%Smooth_Nwp_Fields_Flag_Temp = symbol%NO
   
   !--- do various 101 level NWP Profiles
   Acha_RTM_NWP%P_Prof = P_Std_Rtm

   Acha_RTM_NWP%T_Prof => Rtm(Inwp,Jnwp)%T_Prof 
   Acha_RTM_NWP%Z_Prof => Rtm(Inwp,Jnwp)%Z_Prof 

   !------------------------------------------------------
   ! Before smoothing profiles, ensure that all required
   ! rtm profiles are populated, if not, skip smoothing
   !------------------------------------------------------
   if ((Rtm(Inwp,Jnwp)%Flag == symbol%YES) .and. &
       (Rtm(Inwp_x,Jnwp)%Flag == symbol%YES) .and. &
       (Rtm(Inwp,Jnwp_x)%Flag == symbol%YES) .and. &
       (Rtm(Inwp_x,Jnwp_x)%Flag == symbol%YES)) then

        Acha_RTM_NWP%Smooth_Nwp_Fields_Flag_Temp = symbol%YES
        
        Acha_RTM_NWP%T_Prof_1 => Rtm(Inwp_x,Jnwp)%T_Prof 
        Acha_RTM_NWP%T_Prof_2 => Rtm(Inwp,Jnwp_x)%T_Prof 
        Acha_RTM_NWP%T_Prof_3 => Rtm(Inwp_x,Jnwp_x)%T_Prof 

        Acha_RTM_NWP%Z_Prof_1 => Rtm(Inwp_x,Jnwp)%Z_Prof 
        Acha_RTM_NWP%Z_Prof_2 => Rtm(Inwp,Jnwp_x)%Z_Prof 
        Acha_RTM_NWP%Z_Prof_3 => Rtm(Inwp_x,Jnwp_x)%Z_Prof
        
   endif
   
   !---- RTM profiles
 
   !--- populate radiance and transmission profiles
   if (Acha_Input%Chan_On_67um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_67um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Rad_Atm_Profile
     Acha_RTM_NWP%Atm_Trans_Prof_67um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Trans_Atm_Profile
     Acha_RTM_NWP%Black_Body_Rad_Prof_67um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(27)%Rad_BB_Cloud_Profile
   endif

   if (Acha_Input%Chan_On_85um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_85um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(29)%Rad_Atm_Profile
     Acha_RTM_NWP%Atm_Trans_Prof_85um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(29)%Trans_Atm_Profile
   endif

   if (Acha_Input%Chan_On_11um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_11um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Rad_Atm_Profile
      Acha_RTM_NWP%Atm_Trans_Prof_11um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Trans_Atm_Profile
      Acha_RTM_NWP%Black_Body_Rad_Prof_11um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(31)%Rad_BB_Cloud_Profile
   endif
   
   if (Acha_Input%Chan_On_12um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_12um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(32)%Rad_Atm_Profile
      Acha_RTM_NWP%Atm_Trans_Prof_12um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(32)%Trans_Atm_Profile
   endif

   if (Acha_Input%Chan_On_133um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_133um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(33)%Rad_Atm_Profile
      Acha_RTM_NWP%Atm_Trans_Prof_133um => Rtm(Inwp,Jnwp)%d(Ivza)%ch(33)%Trans_Atm_Profile
   endif
    
 end subroutine ACHA_FETCH_PIXEL_NWP_RTM


end module ACHA_SERVICES_MOD
