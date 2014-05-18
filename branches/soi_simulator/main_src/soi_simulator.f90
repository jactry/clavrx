!-----------------------------------------------------------------------
!  This module holds the radiative transfer quantities needed for
!  the algorithms
!
!
!-----------------------------------------------------------------------
module SOI_SIMULATOR

  use CONSTANTS
  use SOI
  use NUMERICAL_ROUTINES
  use RTM_COMMON
  use PIXEL_COMMON
  use NWP_COMMON
  use PLANCK

  implicit none
  private:: SETUP_SINGLE_SCATTERING_PROPS
  public:: COMPUTE_CLEAR_RAD, &
           COMPUTE_CLOUD_PROFILES_GFS,  &
           COMPUTE_CLOUD_RAD

  real, private, save , dimension(4)::  &
        Qe_log10_065_coef_water, wo_log10re_065_coef_water, g_log10re_065_coef_water, &
        Qe_log10_37_coef_water, wo_log10re_37_coef_water, g_log10re_37_coef_water, &
        Qe_log10_67_coef_water, wo_log10re_67_coef_water, g_log10re_67_coef_water, &
        Qe_log10_85_coef_water, wo_log10re_85_coef_water, g_log10re_85_coef_water, &
        Qe_log10_11_coef_water, wo_log10re_11_coef_water, g_log10re_11_coef_water, &
        Qe_log10_12_coef_water, wo_log10re_12_coef_water, g_log10re_12_coef_water, &
        Qe_log10_13_coef_water, wo_log10re_13_coef_water, g_log10re_13_coef_water
  
  real, private, save, dimension(4)::  &
        Qe_log10_065_coef_ice, wo_log10re_065_coef_ice, g_log10re_065_coef_ice, &
        Qe_log10_37_coef_ice, wo_log10re_37_coef_ice, g_log10re_37_coef_ice, &
        Qe_log10_67_coef_ice, wo_log10re_67_coef_ice, g_log10re_67_coef_ice, &
        Qe_log10_85_coef_ice, wo_log10re_85_coef_ice, g_log10re_85_coef_ice, &
        Qe_log10_11_coef_ice, wo_log10re_11_coef_ice, g_log10re_11_coef_ice, &
        Qe_log10_12_coef_ice, wo_log10re_12_coef_ice, g_log10re_12_coef_ice, &
        Qe_log10_13_coef_ice, wo_log10re_13_coef_ice, g_log10re_13_coef_ice

   integer, private, save:: Setup_Properties_Flag = 1

contains

!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
subroutine COMPUTE_CLOUDY_OBSERVATIONS_USING_SOI(Number_Elements,Number_Lines,Number_Levels)

   integer, intent(in):: Number_Elements
   integer, intent(in):: Number_Lines
   integer, intent(in):: Number_Levels
   integer:: Number_Levels_Rtm
   integer:: Number_Layers
   real, dimension(:), allocatable:: Cld_Opd_Profile, Qe_reference_Profile, Qe_Profile, &
                                     g_Profile, wo_Profile, Opd_Profile, Gas_Opd_Profile, B_Profile 
   integer, dimension(:), allocatable, save:: Lev_Idx_Rtm_Nwp
   real:: Bsfc, Bspace, mu_obs
   integer:: Elem_Idx,Line_Idx,Chan_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx,Zen_Idx,Lev_Idx,Lay_Idx
   integer:: Number_Levels_Rtm
   real:: Cld_Opd, Cld_wo, Rad_Toa
   real:: Surface_Emissivity

   Number_Layers = Number_Levels - 1
   Number_Levels_Rtm = size(P_Std_Rtm)

   !--- compute liquid and ice water layer profiles
   if (Setup_Properties_Flag) then
      call SETUP_SINGLE_SCATTERING_PROPS()
      allocate(Lev_Idx_Rtm_Nwp(Number_Levels_Rtm))
      call SETUP_NWP_TO_RTM_MAPPING(Press_Std_Rtm, Press_Std_Nwp, Lev_Idx_Rtm_Nwp)
      Setup_Properties_Flag = 0
   endif


   !---- allocate local profile vectors
   allocate(Cld_Opd_Profile(Number_Layers))
   allocate(Qe_reference_Profile, Qe_Profile, g_Profile, wo_Profile,  &
            Opd_Profile, Gas_Opd_Profile, source = Cld_Opd_Profile)
   allocate(B_Profile(Number_Levels))


   Element_Loop: do Elem_Idx = 1, Number_Elements
       Line_Loop: do Line_Idx = 1, Number_Lines

          if (Bad_Pixel_Mask(Elem_Idx,Line_Idx)) cycle

          Lon_Nwp_Idx = I_Nwp(Elem_Idx,Line_Idx) 
          Lat_Nwp_Idx = J_Nwp(Elem_Idx,Line_Idx) 
          Zen_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx) 

          if (Lon_Nwp_Idx < 1 .or. Lat_Nwp_Idx < 1)  cycle

          !filter clear?

          !--- loop over channels
          Chan_Loop: do Chan_Idx = 1, 36

            if (Chan_Idx /= 20 .and. & 
                Chan_Idx /= 27 .and. & 
                Chan_Idx /= 29 .and. & 
                Chan_Idx /= 31 .and. & 
                Chan_Idx /= 32 .and. & 
                Chan_Idx /= 33) then

                 cycle

             endif

             !--- convert trans to nwp levels
             Gas_Opd_Profile = COMPUTE_LAYER_GAS_OPD_PROFILE( &
                                   Rtm(Lon_Nwp_Idx,Lat_Nwp_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Profile)

             Surface_Emissivity = ch(Chan_Idx)%Sfc_Emiss(Elem_Idx,Line_Idx)
             Bsfc = PLANCK_RAD_FAST(Chan_Idx, Tsfc_Nwp_Pix(Elem_Idx,Line_Idx))
             Bspace = 0.0
             mu_obs = Coszen(Elem_Idx,Line_Idx)

             Level_Loop: do Lev_Idx = 1, Number_Levels
                B_Profile(Lev_Idx) = PLANCK_RAD_FAST(Chan_Idx, T_Prof_Nwp(Lev_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx))
                Trans_Gas(Lev_Idx) = Rtm(Lon_Nwp_Idx,Lat_Nwp_Idx)%d(Zen_Idx)%ch(Chan_Idx)% &
                                     Trans_Atm_Profile(Lev_Idx_Rtm_Nwp(Lev_Idx))
             enddo Level_Loop

             Layer_Loop: do Lay_Idx = 1, Number_Layers

               Gas_Opd_Profile(Lay_Idx) = -1.0*log(Trans_Gas(Lay_Idx+1)/Trans_Gas(Lay_Idx)) * mu_obs

               if (Cld_Phase_Prof_Nwp(Lay_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx) == sym%WATER_PHASE) then
                   call COMPUTE_CLOUD_OPT_PROP(sym%WATER_PHASE,Cld_Reff_Prof_Nwp(Lay_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx), &
                                Qe_ref_Profile(Lay_Idx),Qe_Profile(Lay_Idx),g_Profile(Lay_Idx),Cld_wo)
               elseif (Cld_Phase_Prof_Nwp(Lay_Idx) == sym%ICE_PHASE) then
                   call COMPUTE_CLOUD_OPT_PROP(sym%ICE_PHASE,Cld_Reff_Prof_Nwp(Lay_Idx,Lon_Nwp_Idx,Lat_Nwp_Idx), &
                                Qe_Reference_Profile(Lay_Idx),Qe_Profile(Lay_Idx),g_Profile(Lay_Idx),Cld_wo)
               else
                   Qe_Reference_Profile(Lay_Idx) = 1.0
                   Qe_Profile(Lay_Idx) = 0.0
                   g_Profile(Lay_Idx) = 0.0
                   Cld_wo = 0.0
               endif

               Cld_Opd = Cld_Opd_Profile(Lay_Idx) * Qe_Profile(Lay_Idx) / Qe_Reference_Profile(Lay_Idx)
               Opd_Profile(Lay_Idx) = Cld_Opd + Gas_Opd_Profile(Lay_Idx)

               if (Opd_Profile(Lay_Idx) > 0.00) then 
                wo_Profile(Lay_Idx) = Cld_wo*Cld_Opd / Opd_Profile(Lay_Idx)
               else
                wo_Profile(Lay_Idx) = Cld_wo
               endif
             enddo Layer_Loop 

             if (Cwp_Nwp_Pix(Elem_Idx,Line_Idx) > 1.0) then
                call FORWARD_MODEL_SOI(Number_Layers,Opd_Profile,wo_Profile,g_Profile,B_Profile, &
                             Surface_Emissivity,Bsfc, Bspace,mu_obs,Rad_Toa)

             else

                call FORWARD_MODEL_ABS_APPROX(Number_Layers,Opd_Profile,wo_Profile,g_Profile,B_Profile, &
                             Surface_Emissivity,Bsfc, Bspace,mu_obs,Rad_Toa)

             endif
        
             ch(Chan_Idx)%Rad_Toa(Elem_Idx,Line_Idx) = Rad_Toa 
             ch(Chan_Idx)%Bt_Toa(Elem_Idx, Line_Idx) = PLANCK_TEMP_FAST(Chan_Idx, ch(Chan_Idx)%Rad_Toa(Elem_Idx,Line_Idx)) 

          enddo Chan_Loop

       enddo Line_Loop
   enddo Element_Loop


end subroutine

!-----------------------------------------------------------------
subroutine ICE_SCATTERING_PROPERTIES(log10_re, Qe_vis,  &
                                     Qe_37, wo_37, g_37,&
                                     Qe_67, wo_67, g_67,&
                                     Qe_85, wo_85, g_85,&
                                     Qe_11, wo_11, g_11,&
                                     Qe_12, wo_12, g_12,&
                                     Qe_13, wo_13, g_13)

  real, intent(in):: log10_re
  real, intent(out):: Qe_vis
  real, intent(out):: Qe_37, wo_37, g_37
  real, intent(out):: Qe_67, wo_67, g_67
  real, intent(out):: Qe_85, wo_85, g_85
  real, intent(out):: Qe_11, wo_11, g_11
  real, intent(out):: Qe_12, wo_12, g_12
  real, intent(out):: Qe_13, wo_13, g_13


subroutine SETUP_SINGLE_SCATTERING_PROPERTIES()
  !--- ice - severly roughend agg columns b = 0.1
  Qe_log10_065_coef_water =   (/   2.2697,  -0.3151,    0.1360,  -0.0221 /)
  wo_log10re_065_coef_water = (/   1.0000,  -0.0000,   0.0000,  -0.0000 /)
  g_log10re_065_coef_water  = (/   0.7418,   0.0074,   0.0027,  -0.0018 /)
  Qe_log10re_37_coef_water = (/    2.2596,   0.5961,  -0.8401,   0.2461 /)
  wo_log10re_37_coef_water = (/    0.9381,   0.1148,  -0.3192,   0.0822 /)
  g_log10re_37_coef_water  = (/   1.0109,  -0.7069,   0.6398,  -0.1485 /)
  Qe_log10re_67_coef_water = (/  -1.3148,   8.8276,  -6.8874,   1.6612 /)
  wo_log10re_67_coef_water = (/   0.5835,   0.5432,  -0.6682,   0.1919 /)
  g_log10re_67_coef_water  = (/   0.7659,   0.1844,  -0.0288,  -0.0048 /)
  Qe_log10re_85_coef_water = (/  -2.5020,  10.6847,  -7.6003,   1.6969 /)
  wo_log10re_85_coef_water = (/   0.4144,   1.1855,  -1.1318,   0.2844 /)
  g_log10re_85_coef_water  = (/   0.7803,   0.0290,   0.0918,  -0.0274 /)
  Qe_log10re_11_coef_water = (/  -0.1574,   2.5725,  -0.9895,   0.1156 /)
  wo_log10re_11_coef_water = (/  -0.0500,   0.6949,  -0.2659,   0.0312 /)
  g_log10re_11_coef_water  = (/   0.4575,   0.9239,  -0.5672,   0.1183 /)
  Qe_log10re_12_coef_water = (/  -0.1559,   4.3089,  -2.8133,   0.6008 /)
  wo_log10re_12_coef_water = (/   0.0211,   0.7949,  -0.4311,   0.0825 /)
  g_log10re_12_coef_water  = (/   0.4196,   0.9945,  -0.6360,   0.1368 /)
  Qe_log10re_13_coef_water = (/  -0.6477,   6.4556,  -4.7423,   1.0951 /)
  wo_log10re_13_coef_water = (/   0.1532,   0.7362,  -0.4988,   0.1139 /)
  g_log10re_13_coef_water  = (/   0.3782,   1.0423,  -0.6597,   0.1408 /)

  !--- water - mie - spheres b = 0.1
  Qe_log10_065_coef_ice = (/   2.4692,  -0.6383,   0.3301,  -0.0616 /)
  wo_log10re_065_coef_ice = (/   1.0001,  -0.0002,   0.0003,  -0.0001 /)
  g_log10re_065_coef_ice  = (/   0.7283,   0.2926,  -0.2212,   0.0626 /)
  Qe_log10re_37_coef_ice = (/   4.2020,  -2.7306,   1.0465,  -0.1023 /)
  wo_log10re_37_coef_ice = (/   0.9653,   0.1125,  -0.2001,   0.0268 /)
  g_log10re_37_coef_ice  = (/   1.0217,  -0.7832,   0.7431,  -0.1864 /)
  Qe_log10re_67_coef_ice = (/  -3.0490,  15.1654, -12.4136,   3.0822 /)
  wo_log10re_67_coef_ice = (/   0.5288,   1.0839,  -1.2979,   0.3812 /)
  g_log10re_67_coef_ice  = (/   0.6315,   0.4040,  -0.1805,   0.0331 /)
  Qe_log10re_85_coef_ice = (/  -4.2799,  15.7967, -11.4845,   2.5978 /)
  wo_log10re_85_coef_ice = (/   0.0866,   2.1334,  -2.0143,   0.5305 /)
  g_log10re_85_coef_ice  = (/   0.3246,   1.1988,  -0.8124,   0.1903 /)
  Qe_log10re_11_coef_ice = (/  -0.8698,   2.8463,   0.2468,  -0.4790 /)
  wo_log10re_11_coef_ice = (/  -0.3134,   1.6573,  -1.0505,   0.2119 /)
  g_log10re_11_coef_ice  = (/  -0.0801,   2.0860,  -1.3774,   0.3018 /)
  Qe_log10re_12_coef_ice = (/  -0.4959,   3.1944,  -0.9197,  -0.0206 /)
  wo_log10re_12_coef_ice = (/  -0.1412,   0.7761,  -0.2803,   0.0273 /)
  g_log10re_12_coef_ice  = (/  -0.2389,   2.3538,  -1.5301,   0.3299 /)
  Qe_log10re_13_coef_ice = (/  -0.6141,   4.8231,  -2.6320,   0.4493 /)
  wo_log10re_13_coef_ice = (/  -0.1321,   0.8492,  -0.3691,   0.0546 /)
  g_log10re_13_coef_ice  = (/  -0.3462,   2.4821,  -1.5888,   0.3382 /)

end subroutine SETUP_SINGLE_SCATTERING_PROPERTIES

subroutine COMPUTE_LAYER_GAS_OPD_PROFILE(Trans_Profile_Rtm,Press_Profile_Rtm, Press_Profile_Nwp,Gas_Profile_Nwp)
   real, dimension(:), intent(in):: Trans_Profile_Rtm
   real, dimension(:), intent(in):: Press_Profile_Rtm
   real, dimension(:), intent(in):: Press_Profile_Nwp
   real, dimension(:), intent(in):: Gas_Opd_Profile_Nwp

   Num_Levels_Rtm = size(Press_Profile_Rtm)
   Num_Levels_Nwp = size(Press_Profile_Nwp)

   
end subroutine
!------------------------------------------------------------------------------
! Determine with RTM level is closest to an NWP level
! assume both NWP and RTM profiles start at TOA
!
!------------------------------------------------------------------------------
subroutine SETUP_NWP_TO_RTM_MAPPING(Pressure_Profile_Rtm, Pressure_Profile_Nwp, Lev_Idx_Rtm_Nwp)
   real, intent(in), dimension(:):: Pressure_Profile_Rtm
   real, intent(in), dimension(:):: Pressure_Profile_Nwp
   integer, intent(in), dimension(:):: Lev_Idx_Rtm_Nwp
   integer:: Lev_Idx_Rtm
   integer:: Lev_Idx_Nwp
   integer:: Number_Levels_Rtm
   integer:: Number_Levels_Nwp

   Number_Levels_Rtm = size(Pressure_Profile_Rtm)  
   Number_Levels_Nwp= size(Pressure_Profile_Nwp)  

  Nwp_Level_Loop: do Lev_Idx_Nwp = 1, Number_Levels_Nwp

    call LOCATE(Pressure_Profile_Rtm, Number_Levels_Rtm, Pressure_Profile_Nwp(Lev_Idx_Nwp),  &
                Lev_Idx_Rtm_Nwp(Lev_Idx_Nwp))

  enddo  Nwp_Level_Loop

end subroutine SETUP_NWP_TO_RTM_MAPPING

end module SOI_SIMULATOR
