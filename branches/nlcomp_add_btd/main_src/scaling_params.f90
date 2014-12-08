!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: scaling_params.f90 (src)
!       SCALING_PARAMETERS (program)
!
! PURPOSE: 
!         Module contains the scaling parameters for all scaling/unscaling performed
!         in CLAVR-x
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
module SCALING_PARAMETERS

 use CONSTANTS
 implicit none

 integer(kind=int4), parameter, public:: one_byte_max = 127, &   !(2**8)/2 -1
                                          one_byte_min = -127     !-(2**8)/2
 integer(kind=int4), parameter, public:: two_byte_max = 32767, & !(2**15)/2 - 1
                                          two_byte_min = -32767   !-(2**15)/2
  
 real, parameter, public:: Min_Ref_Ch1 = -2.0, Max_Ref_Ch1 = 120.0
 real, parameter, public:: Min_Uni_Ch1 =  0.0, Max_Uni_Ch1 = 10.0
 real, parameter, public:: Min_Ref_Ch2 = -2.0, Max_Ref_Ch2 = 120.0
 real, parameter, public:: Min_Uni_Ch2 =  0.0, Max_Uni_Ch2 = 10.0
 real, parameter, public:: Min_Ref_Ch3 = -2.0, Max_Ref_Ch3 = 120.0
 real, parameter, public:: Min_Ref_Ch4 = -2.0, Max_Ref_Ch4 = 120.0
 real, parameter, public:: Min_Ref_Ch5 = -2.0, Max_Ref_Ch5 = 120.0
 real, parameter, public:: Min_Uni_Ch5 =  0.0, Max_Uni_Ch5 = 10.0
 real, parameter, public:: Min_Ref_Ch6 = -2.0, Max_Ref_Ch6 = 120.0
 real, parameter, public:: Min_Ref_Ch7 = -2.0, Max_Ref_Ch7 = 120.0
 real, parameter, public:: Min_Ref_Ch8 = -2.0, Max_Ref_Ch8 = 120.0
 real, parameter, public:: Min_Ref_Ch9 = -2.0, Max_Ref_Ch9 = 120.0
 real, parameter, public:: Min_Ref_Ch17 = -2.0, Max_Ref_Ch17 = 120.0
 real, parameter, public:: Min_Ref_Ch18 = -2.0, Max_Ref_Ch18 = 120.0
 real, parameter, public:: Min_Ref_Ch19 = -2.0, Max_Ref_Ch19 = 120.0
 real, parameter, public:: Min_Ref_Ch20 = -20.0, Max_Ref_Ch20 = 80.0
 real, parameter, public:: Min_Ref_Ch26 = -2.0, Max_Ref_Ch26 = 60.0
 real, parameter, public:: Min_Ref_ChDNB = -2.0, Max_Ref_ChDNB = 150.0
 real, parameter, public:: Min_Ref_ChDNB_lunar = -2.0, Max_Ref_ChDNB_lunar = 150.0
 real, parameter, public:: Min_Rad_Ch20 = 0.0, Max_Rad_Ch20 = 4.0
 real, parameter, public:: Min_Rad_Ch31 = 0.0, Max_Rad_Ch31 = 200.0
 real, parameter, public:: Min_Rad_Ch32 = 0.0, Max_Rad_Ch32 = 200.0
 real, parameter, public:: Min_Bt20 = 180.0, Max_Bt20 = 340.0
 real, parameter, public:: Min_Bt21 = 180.0, Max_Bt21 = 340.0
 real, parameter, public:: Min_Bt22 = 180.0, Max_Bt22 = 340.0
 real, parameter, public:: Min_Bt23 = 180.0, Max_Bt23 = 340.0
 real, parameter, public:: Min_Bt24 = 180.0, Max_Bt24 = 340.0
 real, parameter, public:: Min_Bt25 = 180.0, Max_Bt25 = 340.0
 real, parameter, public:: Min_Bt27 = 180.0, Max_Bt27 = 340.0
 real, parameter, public:: Min_Bt28 = 180.0, Max_Bt28 = 340.0
 real, parameter, public:: Min_Bt29 = 180.0, Max_Bt29 = 340.0
 real, parameter, public:: Min_Bt30 = 180.0, Max_Bt30 = 340.0
 real, parameter, public:: Min_Bt31 = 180.0, Max_Bt31 = 340.0
 real, parameter, public:: Min_Bt32 = 180.0, Max_Bt32 = 340.0
 real, parameter, public:: Min_Bt33 = 180.0, Max_Bt33 = 340.0
 real, parameter, public:: Min_Bt34 = 180.0, Max_Bt34 = 340.0
 real, parameter, public:: Min_Bt35 = 180.0, Max_Bt35 = 340.0
 real, parameter, public:: Min_Bt36 = 180.0, Max_Bt36 = 340.0
 real, parameter, public:: Min_Btd_Ch31_Ch32 = -4.0, Max_Btd_Ch31_Ch32 = 12.0
 real, parameter, public:: Min_Btd_3b_4 = -20.0, Max_Btd_3b_4 = 100.0
 real, parameter, public:: Min_Ref_Ch1_std = 0.0, Max_Ref_Ch1_std = 20.0
 real, parameter, public:: Min_Ref_Ch2_std = 0.0, Max_Ref_Ch2_std = 20.0
 real, parameter, public:: Min_Ref_Ch6_std = 0.0, Max_Ref_Ch6_std = 20.0
 real, parameter, public:: Min_Ref_Ch20_std = 0.0, Max_Ref_Ch20_std = 20.0
 real, parameter, public:: Min_Rad_Ch20_std = 0.0, Max_Rad_Ch20_std = 4.0
 real, parameter, public:: Min_Rad_Ch31_std = 0.0, Max_Rad_Ch31_std = 20.0
 real, parameter, public:: Min_Rad_Ch32_std = 0.0, Max_Rad_Ch32_std = 20.0
 real, parameter, public:: Min_Bt20_std = 0.0, Max_Bt20_std = 20.0
 real, parameter, public:: Min_Bt31_std = 0.0, Max_Bt31_std = 20.0
 real, parameter, public:: Min_Bt32_std = 0.0, Max_Bt32_std = 20.0
 real, parameter, public:: Min_t45_std = 0.0, Max_t45_std = 12.0
 real, parameter, public:: Min_Psfc = 700.0, Max_Psfc = 1100.0
 real, parameter, public:: Min_Pmsl = 850.0, Max_Pmsl = 1100.0
 real, parameter, public:: Min_Kindex = -40.0, Max_Kindex = 80.0
 real, parameter, public:: Min_Tsfc = 220.0, Max_Tsfc = 340.0, Min_Tsfc_std = 0.0, Max_Tsfc_std = 20.0
 real, parameter, public:: Min_Sst = 265.0, Max_Sst = 315.0, Min_Sst_std = 0.0, Max_Sst_std = 20.0
 real, parameter, public:: Min_Lst = 220.0, Max_Lst = 340.0, Min_Lst_std = 0.0, Max_Lst_std = 20.0
 real, parameter, public:: Min_Ndvi = -0.5, Max_Ndvi = 1.0, Min_Ndvi_std = 0.0, Max_Ndvi_std = 1.0
 real, parameter, public:: Min_Ndsi = -0.5, Max_Ndsi = 1.0
 real, parameter, public:: Min_Zen = 0.0, Max_Zen = 90.0
 real, parameter, public:: Min_Relaz = 0.0, Max_Relaz = 180.0
 real, parameter, public:: Min_Solaz = -180.0, Max_Solaz = 180.0
 real, parameter, public:: Min_Solzen = 0.0, Max_Solzen = 180.0
 real, parameter, public:: Min_Sataz = -180.0, Max_Sataz = 180.0
 real, parameter, public:: Min_Glintzen = 0.0, Max_Glintzen = 180.0
 real, parameter, public:: Min_Scatang = 0.0, Max_Scatang = 180.0
 real, parameter, public:: Min_Zsfc = -500, Max_Zsfc = 10000
 real, parameter, public:: Min_Lat = -90.0, Max_Lat = 90.0
 real, parameter, public:: Min_Lon = -180, Max_Lon = 180.0
 real, parameter, public:: Min_Zc = 0.0, Max_Zc = 20000.0
 real, parameter, public:: Min_Zc_Uncer = 0.0, Max_Zc_Uncer = 10000.0
 real, parameter, public:: Min_Pc = 0.0, Max_Pc = 1100.0
 real, parameter, public:: Min_Tc = 160.0, Max_Tc = 320.0, Min_Tc_Std = 0.0, Max_Tc_Std = 40.0
 real, parameter, public:: Min_Tc_Uncer = 0.0, Max_Tc_Uncer = 100.0
 real, parameter, public:: Min_ec = 0.0, Max_ec = 1.0, Min_ec_std = 0.0, Max_ec_std = 1.0
 real, parameter, public:: Min_Beta = 0.0, Max_Beta = 2.0, Min_Beta_std = 0.0, Max_Beta_std = 2.0
 real, parameter, public:: Min_Tau = -0.2, Max_Tau = 160.0, Min_Tau_std = 0.0, Max_Tau_std = 100.0
 real, parameter, public:: Min_Tau_Acha = -0.2, Max_Tau_Acha = 8.0
 real, parameter, public:: Min_Acha_Cost = 0.0, Max_Acha_Cost = 18.0
 real, parameter, public:: Min_Reff = 0.0, Max_Reff = 160.0, Min_Reff_std = 0.0, Max_Reff_std = 100.0
 real, parameter, public:: Min_Hcld = 0.0, Max_Hcld = 4000.0
 real, parameter, public:: Min_Cdnc = 0.0, Max_Cdnc = 1000.0
 real, parameter, public:: Min_lwp = 0.0, Max_lwp = 2000.0, Min_lwp_std = 0.0, Max_lwp_std = 100.0
 real, parameter, public:: Min_iwp = 0.0, Max_iwp = 2000.0, Min_iwp_std = 0.0, Max_iwp_std = 100.0
 real, parameter, public:: Min_frac = 0.0, Max_frac = 1.0
 real, parameter, public:: Min_tpw = 0.0, Max_tpw = 10.0
 real, parameter, public:: Min_Ozone = 100.0, Max_Ozone = 550.0
 real, parameter, public:: Min_rh = 0.0, Max_rh = 110.0
 real, parameter, public:: Min_hght500 = 4500.0, Max_hght500 = 6500.0
 real, parameter, public:: Min_aot = -0.2, Max_aot = 5.0, Min_aot_std = 0.0, Max_aot_std = 1.0
 real, parameter, public:: Min_olr = 50.0, Max_olr = 350.0, Min_olr_std = 0.0, Max_olr_std = 100.0
 real, parameter, public:: Min_Insol = 0.0, Max_Insol = 1500.0
 real, parameter, public:: Min_weasd = 0.0, Max_weasd = 5500.0
 real, parameter, public:: Min_wndspd = 0.0, Max_wndspd = 50.0 
 real, parameter, public:: Min_wnddir = 0.0, Max_wnddir = 360.0 
 real, parameter, public:: Min_ems_Ch20 = 0.5, Max_ems_Ch20 = 3.0 
 real, parameter, public:: Min_etropo = -0.5, Max_etropo = 1.2 
 real, parameter, public:: Min_Ttropo = 160.0, Max_Ttropo = 260
 real, parameter, public:: Min_rsr = -2.0, Max_rsr = 10.0, Min_rsr_std = 0.0, Max_rsr_std = 2.0
 real, parameter, public:: Min_ash_mass = 0.0, Max_ash_mass = 100.0, Min_ash_mass_std = 0.0, Max_ash_mass_std = 100.0   !ton/km^2
 real, parameter, public:: Min_albedo = 0.0, Max_albedo = 1.0, Min_transmission = 0.0, Max_transmission = 1.0 
 real, parameter, public:: Min_dcomp_atmos_vis = 0.0,Max_dcomp_atmos_vis = 1.5, &
                           Min_dcomp_atmos_ir = 0.0, Max_dcomp_atmos_ir = 1.5 
 real, parameter, public:: Min_Sfc_Ems = 0.75, Max_Sfc_Ems = 1.0 
 real, parameter, public:: Min_Trans = 0.0, Max_Trans = 1.0 
 real, parameter, public:: Min_Ch31_Rad_Atm = 0.0, Max_Ch31_Rad_Atm = 100.0 
 real, parameter, public:: Min_Ch31_Rad_Atm_Dwn = 0.0, Max_Ch31_Rad_Atm_Dwn = 50.0 
 real, parameter, public:: Min_Bt_Covar = -10.0, Max_Bt_Covar = 10.0 
 real, parameter, public:: Min_Cwp = 0.0, Max_Cwp = 1200.0
 real, parameter, public:: Min_Rain_Rate = 0.0, Max_Rain_Rate = 32.0
 real, parameter, public:: Min_Binary_Mask = 0.0, Max_Binary_Mask = 1.0
 real, parameter, public:: Min_Cld_Mask = 0.0, Max_Cld_Mask = 3.0
 real, parameter, public:: Min_Cld_Type = 0.0, Max_Cld_Type = 13.0
 real, parameter, public:: Min_Sfc_Type = 0.0, Max_Sfc_Type = 13.0
 real, parameter, public:: Min_Land_Class = 0.0, Max_Land_Class = 7.0
 real, parameter, public:: Min_Snow_Class = 1.0, Max_Snow_Class = 3.0
 real, parameter, public:: Min_Coast_Class = 1.0, Max_Coast_Class = 10.0
 real, parameter, public:: Min_Cld_Phase = 0.0, Max_Cld_Phase = 5.0
 real, parameter, public:: Min_Alt = 0.0, Max_Alt = 100000.0

end module SCALING_PARAMETERS
