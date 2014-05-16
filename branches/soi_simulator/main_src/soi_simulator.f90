!-----------------------------------------------------------------------
!  This module holds the radiative transfer quantities needed for
!  the algorithms
!
!
!-----------------------------------------------------------------------
module SOI_SIMULATOR

  use CONSTANTS
  use SOI
  use NWP_COMMON

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

   integer, private, save:: Setup_Properies_Flag = 1

contains

!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
subroutine COMPUTE_CLOUDY_OBSERVATIONS_USING_SOI()

   !--- compute liquid and ice water layer profiles
   if (Setup_Properties_Flag) then
      call SETUP_SINGLE_SCATTERING_PROPS()
      Setup_Properties_Flag = 0
   endif


   Element_Loop: do Elem_Idx = 1, Number_Elements
       Line_Loop: do Line_Idx = 1, Number_Lines

          if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) cycle

          Lon_Nwp_Idx = I_Nwp(Elem_Idx,Line_Idx) 
          Lat_Nwp_Idx = J_Nwp(Elem_Idx,Line_Idx) 
          Zen_Idx = Vza_Rtm_Idx(Elem_Idx,Line_Idx) 

          if (Lon_Nwp_Idx < 1 .or. Lat_Nwp_Idx < 1)  cycle

          !filter clear?

          !--- loop over channels
          Chan_Loop: do Chan_Idx = 1, 36

            if (Chan_Idx /= 20 .and. & 
                Chan_Idx /= 27 .and. & 
                Chan_Idx /= 29 .and. & 
                Chan_Idx /= 31 .and. & 
                Chan_Idx /= 32 .and. & 
                Chan_Idx /= 33)

                 cycle

             endif

             !--- convert trans to nwp levels
             Trans_Gas_Profile_Rtm = Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Profile

             Gas_Opd_Profile = COMPUTE_LAYER_GAS_OPD_PROFILE()

             Surface_Emissivity = ch(Chan_Idx)%Sfc_Emiss(Elem_Idx,Line_Idx)
             Bsfc = PLANCK_RAD_FAST(Chan_Idx, Tsfc_Nwp_Pix(Elem_Idx,Line_Idx))
             Bspace = 0.0
             mu_obs = Coszen(Elem_Idx,Line_Idx)

             Level_Loop: do Lev_Idx = 1, Number_Levels
                B_Profile = PLANCK_RAD_FAST(Chan_Idx, T_Profile(Elem_Idx,Line_Idx))
             enddo Level_Loop

             Layer_Loop: do Lay_Idx = 1, Number_Layers

               if (Cld_Phase_Prof_Nwp(Lay_Idx) = sym%WATER_PHASE) then
                   call COMPUTE_CLOUD_OPT_PROP(sym%WATER_PHASE,Cld_Reff_Prof_Nwp(Lay_Idx), &
                                Qe_ref_Profile(Lay_Idx),Qe_Profile(Lay_Idx),g_Profile(Lay_Idx),wo_Profile(Lay_Idx))
               elseif (Cld_Phase_Prof_Nwp(Lay_Idx) = sym%ICE_PHASE) then
                   call COMPUTE_CLOUD_OPT_PROP(sym%ICE_PHASE,Cld_Reff_Prof_Nwp(Lay_Idx), &
                                Qe_ref_Profile(Lay_Idx),Qe_Profile(Lay_Idx),g_Profile(Lay_Idx),wo_Profile(Lay_Idx))
               else
                   Qe_ref_Profile(Lay_Idx) = 1.0
                   Qe_Profile(Lay_Idx) = 0.0
                   g_Profile(Lay_Idx) = 0.0
                   wo_Profile(Lay_Idx) = 0.0
               endif

               Cld_Opd = Cld_Opd_Profile(Lay_Idx) * Qe_Profile(Lay_Idx) / Qe_ref_Profile(Lay_Idx)
               Opd_Profile(Lay_Idx) = Cld_Opd + Gas_Opd_Profile(Lay_Idx)

               if (Opd_Profile(Lay_Idx) > 0.00) then 
                wo_Profile(Lay_Idx) = wo_cld*Cld_Opd / Opd_Profile(Lay_Idx)
               else
                wo_Profile(Lay_Idx) = wo_cld
               endif
             enddo Layer_Loop 

             call FORWARD_MODEL_SOI(Number_Layers,Opd_Profile,wo_Profile,g_Profile,B_Profile, &
                             Surface_Emissivity,Bsfc, Bspace,mu_obs,rad_toa)
        
             ch(Chan_Idx)%Rad_Toa(Elem_Idx,Line_Idx) = rad_toa 

          enddo Chan_Loop


                

       enddo Line_Loop
   enddo Element_Loop

     !--- compute optical properties for this channel (wo,g,tau)

     !--- convert trans to nwp levels
     trans_gas_toa = Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Profile
     sfc_emissivity = ch(Chan_Idx)%Sfc_Emiss(Elem_Idx,Line_Idx)
     call COMPUTE_CLOUD_OPT_PROP(iphase,re_wat,Q_ref_wat,Q_wat,g_wat,wo_wat)
     call COMPUTE_CLOUD_OPT_PROP(iphase,re_ice,Q_ref_ice,Q_ice,g_ice,wo_ice)

     !--- make profiles of layer optical properties
     tau_wat_sum = 0.0
     tau_ice_sum = 0.0

     n = sfc_level(i,j)-1
     do k = 1,n
      tau_gas(k) =  -1.0*log(trans_gas_toa(k+1)/trans_gas_toa(k)) * mu_obs
      tau_wat = (Q_wat/Q_ref) * 1.5 * lwp_nwp(i,j,k) / re_wat
      tau_ice = (Q_ice/Q_ref) * 1.5 * iwp_nwp(i,j,k) / re_ice
      tau_eff(k) = tau_gas(k) + tau_wat + tau_ice
      wo_eff(k) = 0.0
      if (tau_eff(k) > 0.0) then
        wo_eff(k) = (wo_ice*tau_ice + wo_wat*tau_wat) / (tau_eff(k))
      endif
      asym_param_eff(k) = 0.0
      tau_scat_sum = wo_ice*tau_ice + wo_wat*tau_wat
      if (tau_scat_sum > 0.0) then
      asym_param_eff(k) = (g_ice*wo_ice*tau_ice + g_wat*wo_wat*tau_wat) / tau_scat_sum  
     endif

     tau_wat_sum = tau_wat_sum + tau_wat
     tau_ice_sum = tau_ice_sum + tau_ice

   enddo

!--- compute planck emission at surface and at each level

     !--- call forward model


   enddo

end subroutine

!---------------------------------------------------------------------------
! compute cloudy radiances using GFS profiles of cloud water/ice
! and assumed optical properties
!---------------------------------------------------------------------------
subroutine COMPUTE_CLOUD_RAD(i,j,i_sfc,j_sfc,c1,c2,nu_3b,a_3b,b_3b,nu_4,a_4,b_4,&
                             nu_5,a_5,b_5)

integer, intent(in):: i,j,i_sfc,j_sfc
real, intent(in):: c1,c2,nu_3b,a_3b,b_3b,nu_4,a_4,b_4,&
                   nu_5,a_5,b_5
real:: B_mean, T_mean, mu_obs, tau_ice, tau_wat, re_wat, re_ice, g_ice, g_wat, tau_scat_sum, &
       sfc_bb_emission, sfc_emissivity,nu,a,b,rad_toa,bt_toa,wo_ice,wo_wat,rad_layer,trans_toa, &
       Q_ref,Q_ice,Q_wat,tau_wat_sum,tau_ice_sum,space_bb_emission
integer:: k,n, ichan, lev

mu_obs = cos(satzen_nwp(i,j)*dtor)
re_ice = 30.0
re_wat = 10.0
Q_ref = 2.0

space_bb_emission = 0.0

channel_loop: do ichan = 3,5 

 if (ichan == 3) then
    wo_ice = 0.68
    wo_wat = 0.87
    g_ice = 0.90
    g_wat = 0.80
    trans_gas_toa = trans_atm_ch3b_nwp(i,j,:)
    sfc_emissivity = 0.98
    nu = nu_3b
    a = a_3b
    b = b_3b
    Q_ice = 2.2
    Q_wat = 2.3
 endif

 if (ichan == 4) then
    wo_ice = 0.48
    wo_wat = 0.5
    g_ice = 0.96
    g_wat = 0.92
    trans_gas_toa = trans_atm_ch4_nwp(i,j,:)
    sfc_emissivity = 0.98
    nu = nu_4
    a = a_4
    b = b_4
    Q_ice = 2.11
    Q_wat = 1.8
 endif

 if (ichan == 5) then
    wo_ice = 0.50
    wo_wat = 0.39
    g_ice = 0.92
    g_wat = 0.91
    trans_gas_toa = trans_atm_ch5_nwp(i,j,:)
    sfc_emissivity = 0.98
    nu = nu_5
    a = a_5
    b = b_5
    Q_ice = 2.3
    Q_wat = 1.8
 endif

!-------------------------------------------------
! make profiles of layer optical properties
!-------------------------------------------------
tau_wat_sum = 0.0
tau_ice_sum = 0.0

n = sfc_level(i,j)-1
do k = 1,n
   tau_gas(k) =  -1.0*log(trans_gas_toa(k+1)/trans_gas_toa(k)) * mu_obs
   tau_wat = (Q_wat/Q_ref) * 1.5 * lwp_nwp(i,j,k) / re_wat
   tau_ice = (Q_ice/Q_ref) * 1.5 * iwp_nwp(i,j,k) / re_ice
   tau_eff(k) = tau_gas(k) + tau_wat + tau_ice
   wo_eff(k) = 0.0
   if (tau_eff(k) > 0.0) then
    wo_eff(k) = (wo_ice*tau_ice + wo_wat*tau_wat) / (tau_eff(k))
   endif
   asym_param_eff(k) = 0.0
   tau_scat_sum = wo_ice*tau_ice + wo_wat*tau_wat
   if (tau_scat_sum > 0.0) then
    asym_param_eff(k) = (g_ice*wo_ice*tau_ice + g_wat*wo_wat*tau_wat) / tau_scat_sum  
   endif

   tau_wat_sum = tau_wat_sum + tau_wat
   tau_ice_sum = tau_ice_sum + tau_ice

enddo

!--- compute planck emission at surface and at each level
 sfc_bb_emission =  c1*(nu**3)/(exp((c2*nu)/ &
                        ((tmpsfc_nwp(i_sfc,j_sfc)-a)/b))-1.0)
do lev = 1,n+1
  planck_rad_level(lev) = c1*(nu**3)/(exp((c2*nu)/ &
                        ((level_temperature(i,j,lev)-a)/b))-1.0)
enddo

!------------------------------------------------------------------------
! call SOI
!-------------------------------------------------------------------------
 if (tot_lwp_nwp(i,j) + tot_iwp_nwp(i,j) > 1.0) then
   call FORWARD_MODEL_SOI(n,tau_eff(1:n),wo_eff(1:n),asym_param_eff(1:n),planck_rad_level(1:n+1), &
                        sfc_emissivity,sfc_bb_emission,space_bb_emission,mu_obs,rad_toa)

 else
   call FORWARD_MODEL_ABS_APPROX(n,tau_eff(1:n),wo_eff(1:n),asym_param_eff(1:n),planck_rad_level(1:n+1), &
                       sfc_emissivity,sfc_bb_emission,space_bb_emission, mu_obs,rad_toa)
 endif

bt_toa =  a + b * ((c2*nu) / log( 1.0 + (c1*(nu**3))/rad_toa))

if (ichan == 3) then
  ch3b_rad_cld_nwp(i,j) = rad_toa
  ch3b_bt_cld_nwp(i,j) = bt_toa
elseif (ichan == 4) then
  ch4_rad_cld_nwp(i,j) = rad_toa
  ch4_bt_cld_nwp(i,j) = bt_toa
elseif (ichan == 5) then
  ch5_rad_cld_nwp(i,j) = rad_toa
  ch5_bt_cld_nwp(i,j) = bt_toa
endif

end do channel_loop

end subroutine COMPUTE_CLOUD_RAD

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

end module SOI_SIMULATOR
