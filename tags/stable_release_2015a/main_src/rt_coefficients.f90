!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: rt_ceofficients.f90 (src)
!       RT_COEFFICIENTS (program)
!
! PURPOSE: This module holds the radiative transfer quantities needed for
!          the algorithms
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
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
!--------------------------------------------------------------------------------------
module RT_COEFFICIENTS

  use CONSTANTS
!  use GFS
  use SOI
  use NWP_COMMON

  implicit none
  private
  public:: COMPUTE_CLEAR_RAD, &
           COMPUTE_CLOUD_PROFILES_GFS, COMPUTE_CLOUD_RAD



! real, parameter, private:: g = 9.8, dtor = 0.017453293
  real, parameter, private:: g = 9.8
  
contains

!-------------------------------------------------------------------------
! this routine takes the transmission profile computed by the pCRTM
! and compute the needed radiance and transmission to space terms using 
! a simple non-scattering RT assumption
!-------------------------------------------------------------------------
subroutine COMPUTE_CLEAR_RAD(i,j,c1,c2,nu_3b,a_3b,b_3b,nu_4,a_4,b_4,&
                             nu_5,a_5,b_5)
integer, intent(in):: i,j
real, intent(in):: c1,c2,nu_3b,a_3b,b_3b,nu_4,a_4,b_4,&
                   nu_5,a_5,b_5
real:: B_mean, T_mean
integer:: k


!-----------------------------------------------------------
! compute atmospheric emission from a level to space
! levels go from 0 to nlevels
!------------------------------------------------------------
rad_atm_ch3b_nwp(i,j,:) = 0.0
rad_atm_ch4_nwp(i,j,:) = 0.0
rad_atm_ch5_nwp(i,j,:) = 0.0


do k = 2,sfc_level(i,j)

  T_mean = 0.5*(level_temperature(i,j,k-1) + level_temperature(i,j,k))

!-- ch3b
   B_mean = c1*(nu_3b**3)/(exp((c2*nu_3b)/((T_mean-a_3b)/b_3b))-1.0)
   rad_atm_ch3b_nwp(i,j,k) = rad_atm_ch3b_nwp(i,j,k-1) +  &
           (trans_atm_ch3b_nwp(i,j,k-1) - trans_atm_ch3b_nwp(i,j,k)) * B_mean

!-- ch4
   B_mean = c1*(nu_4**3)/(exp((c2*nu_4)/((T_mean-a_4)/b_4))-1.0)
   rad_atm_ch4_nwp(i,j,k) = rad_atm_ch4_nwp(i,j,k-1) +  &
           (trans_atm_ch4_nwp(i,j,k-1) - trans_atm_ch4_nwp(i,j,k)) * B_mean

!-- ch5
   B_mean = c1*(nu_5**3)/(exp((c2*nu_5)/((T_mean-a_5)/b_5))-1.0)
   rad_atm_ch5_nwp(i,j,k) = rad_atm_ch5_nwp(i,j,k-1) +  &
           (trans_atm_ch5_nwp(i,j,k-1) - trans_atm_ch5_nwp(i,j,k)) * B_mean

enddo

end subroutine COMPUTE_CLEAR_RAD


!----------------------------------------------------------------------------
! COMPUTE CLOUD PROFILES
!----------------------------------------------------------------------------
subroutine COMPUTE_CLOUD_PROFILES_GFS(i,j)

   integer, intent(in):: i,j
   integer:: ilev, ilay, n
   real:: ice_frac_top, ice_frac_bot,clwmr_ice_layer,  &
          clwmr_wat_layer,clwmr_min, lwp_sum, lwp_min

!---------------------------------------------------------------------
! convert clwmr profiles into optical depths
! iwp and lwp are in g/m^2
!---------------------------------------------------------------------
   lwp_nwp(i,j,:) = 0.0
   iwp_nwp(i,j,:) = 0.0
   n = sfc_level(i,j)-1
   do ilay = 1, n
      ice_frac_top  = min(1.0,max(0.0,(273.15 - level_temperature(i,j,ilay))/20.0))
      ice_frac_bot  = min(1.0,max(0.0,(273.15 - level_temperature(i,j,ilay+1))/20.0))
      clwmr_ice_layer = 0.5 * (ice_frac_top*level_clwmr(i,j,ilay) + ice_frac_bot*level_clwmr(i,j,ilay+1))
      clwmr_wat_layer = 0.5 * ((1.0-ice_frac_top)*level_clwmr(i,j,ilay) + (1.0-ice_frac_bot)*level_clwmr(i,j,ilay+1))
      iwp_nwp(i,j,ilay) = 1000.0*50.0 * clwmr_ice_layer * (level_pressure(i,j,ilay+1) - level_pressure(i,j,ilay))/g
      lwp_nwp(i,j,ilay) = 1000.0*50.0 * clwmr_wat_layer * (level_pressure(i,j,ilay+1) - level_pressure(i,j,ilay))/g
   enddo  

   tot_lwp_nwp(i,j) = sum(lwp_nwp(i,j,:))
   tot_iwp_nwp(i,j) = sum(iwp_nwp(i,j,:))

!----------------------------------------------------------------------
! compute number of distinct cloud layers in GFS profile
!----------------------------------------------------------------------
    ncld_layers_nwp(i,j) = 0
    if ((tot_iwp_nwp(i,j) > 2) .and. (tot_lwp_nwp(i,j) > 10)) then
     do ilay = 4, n 
       clwmr_min = 1.0e-05 * level_pressure(i,j,ilay)/1000.0    !
       if (level_clwmr(i,j,ilay) > clwmr_min .and. level_clwmr(i,j,ilay-1) < clwmr_min .and.  &
           level_clwmr(i,j,ilay-2) < clwmr_min) then 
           ncld_layers_nwp(i,j) = ncld_layers_nwp(i,j) + 1
       endif
     enddo
    endif

!   ihigh = 0
!   ilow = 0
!   if ((tot_iwp_nwp(i,j) > 2) .and. (tot_lwp_nwp > 2)) then
!   do ilay = 4, n 
!     iwp_toa = sum(iwp_nwp(i,j,1:ilay-1)
!     lwp_sfc = sum(lwp_nwp(i,j,ilay:n)
!     if ((iwp_toa == tot_iwp_nwp(i,j) .and.  &
!         (lwp_sfc == tot_lwp_nwp(i,j) .and.  &
!   enddo


!--------------------------------------------------------------------
! compute cloud temperature from GFS data
!--------------------------------------------------------------------
    lwp_sum = 0.0
    lwp_min = 2.0
    tc_nwp(i,j) = missing_value_real4
    do ilay = 2, n 
       if ((lwp_sum >= lwp_min) .or. &
           (lwp_nwp(i,j,ilay)+iwp_nwp(i,j,ilay) >= lwp_min)) then
           tc_nwp(i,j) = level_temperature(i,j,ilay)
           exit
       endif
       lwp_sum = lwp_sum + lwp_nwp(i,j,ilay) + iwp_nwp(i,j,ilay)
    enddo

end subroutine COMPUTE_CLOUD_PROFILES_GFS

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

end module RT_COEFFICIENTS
