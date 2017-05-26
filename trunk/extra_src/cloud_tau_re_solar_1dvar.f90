! $Id: cloud_tau_re_solar_1dvar.f90,v 1.31 2013/06/08 18:47:18 wstraka Exp $
Module CLOUD_TAU_RE_SOLAR_MODULE
!======================================================================
!
! module name: CLOUD_TAU_RE_SOLAR_module
!
! Description:
!   Baseline Daytime Cloud optical and microphysical algorithm(s)
!
! author: AndRew Heidinger and the GOES-R AWG Cloud Team
!
! public routines:
!     CLOUD_TAU_RE_SOLAR - estimate tau and re from solar 
!                        reflectance only 
!                      - needs Tc for upstream algorithm
!     DEALLOCATE_MEMORY_FOR_MAIN_CLD_LUT_STRUCTURES
!
! private routines:
!     READ_CLD_REF_LUT_STRUC
!     READ_CLD_EMS_LUT_STRUC
!     INTERPOLATE_REFLECTANCE_VALUE - knowing refl table location, interp 5-dim ref tables
!     INTERPOLATE_EMISSIVITY_VALUE - knowing emiss table location, interp 3-dim ems tables
!     LOCATE_TABLE_POSITION
!     CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER - forward model for reflectance computation
!     CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER - forward model for thermal computation
!     ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE
!     DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES
!     ESTIMATE_TAU_SINGLE_CHANNEL_FIXED_REFF
!     ESTIMATE_ALBEDO_TRANSMISSION_SINGLE_LAYER
!
! version history:
!  04/28/07 - Created based on CLAVR-x algorithms
!
! All channels referenced to MODIS (1=0.65 um, 6 = 1.6 um, 
!                                   7 = 2.1, 20 = 3.75 um)
!
! Notes:
! 1. channel 7 is not included yet because not available on SEVIRI.
! 2. Error handling needs to be included.
!======================================================================
  use CONSTANTS
  use PIXEL_COMMON
  use PIXEL_ROUTINES
  use NUMERICAL_ROUTINES
  use SURFACE_PROPERTIES
  use NWP_COMMON
  use PLANCK
  use RT_UTILITIES
  use HDF

  Implicit None

!-----------------------------------------------------------------------
! declare public routines
!-----------------------------------------------------------------------
  public:: CLOUD_TAU_RE_SOLAR, &
           DEALLOCATE_MEMORY_FOR_MAIN_CLD_LUT_STRUCTURES

!-----------------------------------------------------------------------
! declare private routines
!-----------------------------------------------------------------------
  private:: READ_CLD_REF_LUT_STRUC, &
            READ_CLD_EMS_LUT_STRUC, &
            INTERPOLATE_REFLECTANCE_VALUE,&     !knowing refl table location, interp 5-dim ref tables
            INTERPOLATE_EMISSIVITY_VALUE,&      !knowing emiss table location, interp 3-dim ems tables
            LOCATE_TABLE_POSITION, &
            CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER, &  !forward model for reflectance computation
            CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER, &   !forward model for thermal computation
            ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE, &
            DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES, &
            ESTIMATE_TAU_SINGLE_CHANNEL_FIXED_REFF, &
            ESTIMATE_ALBEDO_TRANSMISSION_SINGLE_LAYER, &
            KNOWING_Z_COMPUTE_T_P_RTM

  !-----------------------------------------------------------------------
  ! define module-wide variables
  !-----------------------------------------------------------------------
  integer(kind=int1),parameter,private::WATER_PHASE = 1  !water index
  integer(kind=int1),parameter,private::MIXED_PHASE = 2  !mixed index
  integer(kind=int1),parameter,private::ICE_PHASE = 3    !ice index
  integer(kind=int1),parameter,private::Nchan_Max = 20   !max channels

  !--- include file
  include 'baseline_cloud_micro_day.inc'

 contains
!=======================================================================
! CLOUD_TAU_RE_SOLAR
! 
! Daytime Cloud Optical and Microphysical Retrieval Using Only Solar
! Reflectance
!
! imode = the mode in which this retrieval is being run and is 
!         determined by the channel selection.
!
! imode = 1 = ch1 and ch6
! imode = 2 = ch1 and ch7  (not coded)
! imode = 3 = ch1 and ch20 total
! imode = 4 = ch1 and ch20 solar only
!=======================================================================
 subroutine CLOUD_TAU_RE_SOLAR(iseg, nseg)

!--- AVHRR specific additions
   integer, intent(in):: iseg
   integer, intent(in):: nseg
   integer, parameter:: Algo_Idx = 1
   integer, parameter:: Nalgo = 1
   character(len=100):: lut_path
   integer, SAVE:: Sc_Id_Prev = 0
!--- end AVHRR specific additions

   integer:: itropo
   integer:: xnwp
   integer:: ynwp
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: ivza
   integer:: isfc
   integer:: imode

   !--- local variables that are aliases to elements of global structures
   real:: Refl_065um
   real:: Refl_160um
   real:: Refl_375um
   real:: Rad_375um
   real:: Bt_11um
   real(kind=real4):: zen
   real(kind=real4):: scat_zen
   real(kind=real4):: sol_zen
   real(kind=real4):: cos_zen
   real(kind=real4):: cos_sol_zen
   real(kind=real4):: Air_Mass
   real(kind=real4):: az
   real(kind=real4):: Ref_std  !reflectance variability
   real(kind=real4):: Ref_mean !reflectance variability
   real(kind=real4):: Ref_var  !reflectance variability
   integer(kind=int4):: cloudphase
   integer(kind=int4):: cloudtype
   integer(kind=int4):: cloudmask
   integer(kind=int4):: sfctype
   integer(kind=int4):: iparam
   real(kind=real4):: Tpw
   real(kind=real4):: Tsfc
   real(kind=real4):: Rad_375um_clr
   real(kind=real4):: T_cloud
   real(kind=real4):: Z_cloud
   real(kind=real4):: P_cloud
   integer(kind=int4):: undetected_cloud
   real(kind=real4):: R4_Dummy

   !--- local pointers to geocat public data structures
   real(kind=real4), dimension(:), pointer:: tlev
   real(kind=real4), dimension(:), pointer:: plev
   real(kind=real4), dimension(:), pointer:: zlev
   real(kind=real4), dimension(:), pointer:: Tpwlev

  !--- 1d-var retrieval arrays
  integer, parameter:: num_obs = 2
  integer, parameter:: num_param = 2
  integer, parameter:: iter_max = 20
  real, parameter:: delta_x_max = 0.2
  real, dimension(num_obs):: y
  real, dimension(num_obs):: f
  real, dimension(num_param):: x
  real, dimension(num_param):: x_ap
  real, dimension(num_param):: delta_x
  real, dimension(num_param):: delta_x_constrained
  real, dimension(num_param):: delta_x_dir
  real:: delta_x_distance
  real:: delta_x_distance_constrained
  real, dimension(num_obs,num_param):: K
  real, dimension(num_param,num_param):: Sa
  real, dimension(num_param,num_param):: Sa_inv
  real, dimension(num_param,num_param):: Sx
  real, dimension(num_param,num_param):: Sx_inv
  real, dimension(num_param,num_param):: E
  real, dimension(num_obs,num_obs):: Sy
  real, dimension(num_obs,num_obs):: Sy_inv

  !--- local variables
  real(kind=real4):: log10tau   !log10 of the optical depth (tau)
  real(kind=real4):: log10re    !log10 of the effective radius (re)
  real(kind=real4):: conv_test  !convergence metric (tested againist Conv_Crit)
  real(kind=real4):: r          !log10re interpolation weight
  real(kind=real4):: s          !log10tau interpolation weight (0-1)
  real(kind=real4):: t          !interpolation weight (0-1)
  real(kind=real4):: t2         !interpolation weight (0-1)
  real(kind=real4):: u          !interpolation weight (0-1)
  real(kind=real4):: v          !interpolation weight (0-1)
  integer(kind=int4):: ir       !interpolation index
  integer(kind=int4):: is       !interpolation index
  integer(kind=int4):: it       !interpolation index
  integer(kind=int4):: it2      !interpolation index
  integer(kind=int4):: iu       !interpolation index
  integer(kind=int4):: iv       !interpolation index
  real(kind=real4):: fm_err     !forward model error
  real(kind=real4):: Conv_Crit  !convergence criteria
  integer(kind=int4):: isingular  !flag noting a singular matrix
  integer(kind=int4):: ifail      !flag noting a failed retrieval
  integer(kind=int4):: iter       !iteration loop index
  integer(kind=int4):: ilev_rtm   !level index for rtm arrays
  integer(kind=int4):: ilev_nwp   !level index for nwp arrays
! real(kind=real4):: pc           !cloud top pressure
! real(kind=real4):: zc           !cloud top height
! real(kind=real4):: tc           !cloud top temperature
  real(kind=real4):: prof_Weight  !profile interpolation weight
  real(kind=real4):: Tpw_ac       !above-cloud Tpw_ac
  integer(kind=int1):: diag_output
  real(kind=real4):: Tau_gas      !gas optical depth
  real(kind=real4):: Tau_aer      !aerosol optical depth
  real(kind=real4):: Tau_ray      !rayleigh optical depth
  real(kind=real4):: Ref_Ss_065um     !single scatter reflectance
  real(kind=real4):: Ref_Ss_160um     !single scatter reflectance
  real(kind=real4):: Ref_Ss_375um     !single scatter reflectance
  real(kind=real4):: Ref_2_clear  !clear ch2 reflectance
  real(kind=real4):: cloud_albedo_ap_view
  real(kind=real4):: cloud_albedo_ap_sun

!--- local variables
  real(kind=real4):: Trans_Ac_065um
  real(kind=real4):: Trans_Ac_160um
  real(kind=real4):: Trans_Ac_375um
  real(kind=real4):: Trans_Ac_375um_zen
  real(kind=real4):: Rad_Ac_375um_Zen
  real(kind=real4):: Trans_Sfc_065um
  real(kind=real4):: Trans_Sfc_160um
  real(kind=real4):: Trans_Sfc_375um
  real(kind=real4):: Alb_Sfc_065um
  real(kind=real4):: Alb_Sfc_160um
  real(kind=real4):: Alb_Sfc_375um
  real(kind=real4):: Emiss_Sfc_375um
  real(kind=real4):: Solar_375um
  real(kind=real4):: sed
  real(kind=real4):: Rad_to_Ref_Fac_375um
  real(kind=real4):: Ref_Sol
  real(kind=real4):: dRef_Sol_dTau
  real(kind=real4):: dRef_Sol_dRe
  real(kind=real4):: Rad_Therm
  real(kind=real4):: Rad_Therm_Toc, Ref_Sol_Toc, f_toc
  real(kind=real4):: dRad_Therm_dTau
  real(kind=real4):: dRad_Therm_dRe
  real(kind=real4):: dB_dT
  real(kind=real4):: Bt
 
!-------------------------------------------------------------------------------
! Executable Code
!-------------------------------------------------------------------------------

!--- covergence criteria
   Conv_Crit = num_param/10.0

!--- set diagnostic output flag
!  diag_output = sym%YES
   diag_output = sym%NO

!--- initialize 1d-var arrays
   y = 0.0
   f = 0.0
   x = 0.0
   x_ap = 0.0
   delta_x = 0.0
   K = 0.0
   Sa = 0.0
   Sa_inv = 0.0
   Sy = 0.0
   Sy_inv = 0.0

!----------- make identity matrix
  E = 0.0
  do iparam = 1,num_param
    E(iparam,iparam) = 1.0
  end do

!--- initialize to missing
 Tau_Dcomp = Missing_Value_Real4
 Reff_Dcomp = Missing_Value_Real4
 Tau_Dcomp_Qf = 0
 Reff_Dcomp_Qf = 0
 Cloud_063um_albedo = Missing_Value_Real4
 Cloud_063um_transmission_view = Missing_Value_Real4
 Cloud_063um_transmission_solar = Missing_Value_Real4

!----------------------------------------------------------------------
! Check available channels and determine the mode of operation
!----------------------------------------------------------------------
if (Sensor%Chan_On_Flag_Default(1) == sym%NO) then
    print *, 'No solar_tau_reff calculations possible without Ref_065um_nom'
endif

!--------------------------------------------------------------------------
! determine which mode is on
! DCOMP_Mode = user selected mode
!--------------------------------------------------------------------------
 imode = 0
 if (DCOMP_Mode == 1 .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) imode = DCOMP_Mode
 if (DCOMP_Mode == 2 .and. Sensor%Chan_On_Flag_Default(7) == sym%YES) imode = DCOMP_Mode
 if (DCOMP_Mode == 3 .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) imode = DCOMP_Mode

 if (imode == 0 .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) imode = 3
 if (imode == 0 .and. Sensor%Chan_On_Flag_Default(7) == sym%YES) imode = 2
 if (imode == 0 .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) imode = 1

 if (imode == 0) then
    print *, 'No solar_tau_reff calculations possible ', Dcomp_Mode, Sensor%Chan_On_Flag_Default(6), Sensor%Chan_On_Flag_Default(7), Sensor%Chan_On_Flag_Default(20)
 endif

!if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
!   imode = 1    !use 1.6 micron
!else
!   imode = 3    !use 3.75 micron total reflectance
!   imode = 4    !use 3.75 micron pseudo reflectance
!endif

  !----------------------------------------------------------------------
  ! store path to ancillary data in a local variable
  !----------------------------------------------------------------------
  lut_path = TRIM(ancil_data_dir)//"/luts/cld/"

  !----------------------------------------------------------------------
  ! check to see if main structures need to be allocated
  !----------------------------------------------------------------------
  if (.NOT. allocated(Cld_Ref)) then
      call ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE(Nalgo)
  endif

  !----------------------------------------------------------------------
  !--- on the first segment, read tables
  !----------------------------------------------------------------------
  if (iseg == 1) then

  !--- check if proper tables are already read in for this sensor
  if ((Sc_Id_Prev == 0) .or. (Sc_Id_Prev /= Sensor%WMO_Id)) then

       !--- check if tables need to be deallocated first
       if (Sc_Id_Prev > 0) then
        call DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES(Algo_Idx)
       endif

       !--- Read channel 1 tables
       if (Sensor%Chan_On_Flag_Default(Chan_Idx_065um) == sym%YES) then
          call READ_TABLES(Algo_Idx,Chan_Idx_065um,lut_path)
       endif

       !--- Read 1.60 um channel 
       if (Sensor%Chan_On_Flag_Default(Chan_Idx_160um) == sym%YES) then
           call READ_TABLES(Algo_Idx,Chan_Idx_160um,lut_path)
       endif

       !--- Read 3.75 um channel 
       if (Sensor%Chan_On_Flag_Default(Chan_Idx_375um) == sym%YES) then
           call READ_TABLES(Algo_Idx,Chan_Idx_375um,lut_path)
       endif

   endif
  end if

  !----------------------------------------------------------------------
  ! check to see if there is any daytime data in this segment
  !----------------------------------------------------------------------
  sol_zen = minval(Geo%Solzen)
  if ((sol_zen >= Solzen_max).and.(sol_zen /= Missing_Value_Real4)) then
    return
  endif

  !======================================================================
  ! Loop over pixels in this segment
  !======================================================================
  line_loop_1: Do Line_Idx= 1, Image%Number_Of_Lines_Read_This_Segment
    element_loop_1: Do Elem_Idx= 1, Image%Number_Of_Elements
       
      !--- check for space pixel
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
        cycle
      end if

      if ((imode == 1) .and. (Sensor%Chan_On_Flag_Default(Chan_Idx_160um) == sym%NO)) then
        cycle
      endif
      if ((imode == 3) .and. (Sensor%Chan_On_Flag_Default(Chan_Idx_375um) == sym%NO)) then
        cycle
      endif

      !--- initilatize quality flag to be lowest quality
      Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 0
      Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 0

      undetected_cloud = sym%NO

      !--- define aliases
      Refl_065um = Missing_Value_Real4
      Refl_160um = Missing_Value_Real4
      Refl_375um = Missing_Value_Real4
      Rad_375um = Missing_Value_Real4
      if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
          Refl_065um = ch(1)%Ref_Toa(Elem_Idx,Line_Idx)          !0.6 micron reflectance
      endif
      if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
         Refl_160um = ch(6)%Ref_Toa(Elem_Idx,Line_Idx)          !1.6 micron reflectance
      endif
      if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
        Refl_375um = ch(20)%Ref_Toa(Elem_Idx,Line_Idx)         !3.9 micron reflectance
        Rad_375um = ch(20)%Rad_Toa(Elem_Idx,Line_Idx)          !3.9 micron reflectance
      endif
      Bt_11um  = ch(31)%Bt_Toa(Elem_Idx,Line_Idx)      !11.0 micron bt
      xnwp = i_nwp(Elem_Idx,Line_Idx)            !nwp longitude cell
      ynwp = j_nwp(Elem_Idx,Line_Idx)            !nwp latitude cell
      itropo = Tropo_Level_Nwp(xnwp,ynwp)  !nwp level associated with tropopause
      ivza = Zen_Idx_Rtm(Elem_Idx,Line_Idx)         !viewing zenith angle bin
      isfc = Sfc_Level_Nwp(xnwp,ynwp)      !first nwp level above surface
      zen = Geo%Satzen(Elem_Idx,Line_Idx)            !viewing zenith angle 
      scat_zen = Geo%scatangle(Elem_Idx,Line_Idx)    !scattering angle 
      sol_zen = Geo%Solzen(Elem_Idx,Line_Idx)        !solar zenith angle 
      cos_zen = Geo%coszen(Elem_Idx,Line_Idx)        !cosine viewing zenith angle 
      Air_Mass = Geo%Airmass(Elem_Idx,Line_Idx)      !air mass factor              !***** CHECK THIS
      cos_sol_zen = Geo%cosSolzen(Elem_Idx,Line_Idx) !cosine solar zenith angle 
      az = Geo%Relaz(Elem_Idx,Line_Idx)              !relative solar azimuth angle 
      T_cloud = ACHA%Tc(Elem_Idx,Line_Idx)       !cloud temperature
      Z_cloud = ACHA%Zc(Elem_Idx,Line_Idx)       !cloud height(check it is in m)


      !--- select cloud type to use
      cloudtype = Cld_Type(Elem_Idx,Line_Idx)

      cloudmask = CLDMASK%Cld_Mask(Elem_Idx,Line_Idx)    !cloud mask
      sfctype = Sfc%Sfc_Type(Elem_Idx,Line_Idx)     !surface type
      Ref_std = Ref_Ch1_Std_3x3(Elem_Idx,Line_Idx)   !measure of reflectance variability
      Ref_mean = Ref_Ch1_Mean_3x3(Elem_Idx,Line_Idx)   !measure of reflectance variability
      Ref_var = 0.0
      if (Ref_mean > 0.0) then
          Ref_var = Ref_std / Ref_mean
      endif
      Ref_var = min(1.0,Ref_var) 
      plev => P_std_nwp
      tlev => T_prof_nwp(:,xnwp,ynwp)
      zlev => Z_prof_nwp(:,xnwp,ynwp)
      Tpwlev => Tpw_prof_nwp(:,xnwp,ynwp)

      Tpw  = Tpw_Nwp_Pix(Elem_Idx,Line_Idx)                   !total column pw
      Tsfc  = Tsfc_Nwp_Pix(Elem_Idx,Line_Idx)                 !surface temperature
      Emiss_Sfc_375um = ch(20)%Sfc_Emiss(Elem_Idx,Line_Idx)     !3.9 micron surface emissivity
      Rad_375um_clr = ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx)  !clear 4 micron radiance
      Solar_375um = Solar_Ch20_Nu                       !solar energy in channel 7
      sed = sun_earth_distance                         !sun earth distance factor

      !--- compute factor for converting channel reflectance into a radiances
      Rad_to_Ref_Fac_375um = pi / cos_sol_zen /  (Solar_375um/sed**2) 

      !--- check for correct sensor and solar geometery
      if (zen >= zen_max) then
        cycle
      end if
      if (sol_zen >= Solzen_max) then
        cycle
      end if

      !--- check for clear pixels 
      if ((cloudmask == sym%CLEAR .or.  &
           cloudmask == sym%PROB_CLEAR .or.  &
           cloudtype == sym%CLEAR_TYPE .or. &
           cloudtype == sym%UNKNOWN_TYPE)) then

        if (process_undetected_cloud_flag == sym%YES) then

          !--- if not glint, assume a water cloud
          if (Sfc%Glint_Mask(Elem_Idx,Line_Idx)  == sym%NO) then
           cloudtype = sym%WATER_TYPE
           undetected_cloud = sym%YES
          else
           !--- assume clear-sky values
           Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 0
           Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 0
           Tau_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
           Reff_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
           Cloud_063um_albedo(Elem_Idx,Line_Idx) = Missing_Value_Real4
           Cloud_063um_transmission_solar(Elem_Idx,Line_Idx) = Missing_Value_Real4
           Cloud_063um_transmission_view(Elem_Idx,Line_Idx) = Missing_Value_Real4
           CYCLE
          endif

        else
           Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 0
           Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 0
           Tau_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
           Tau_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
           Cloud_063um_Albedo(Elem_Idx,Line_Idx) = 0.0
           Cloud_063um_Transmission_Solar(Elem_Idx,Line_Idx) = 1.0
           Cloud_063um_Transmission_View(Elem_Idx,Line_Idx) = 1.0
           CYCLE
        end if


      end if

      !--- check for valid cloud height, if not use obs brightness temp
      if (T_cloud == Missing_Value_Real4) then
        T_cloud = ch(31)%Bt_Toa(Elem_Idx,Line_Idx)
      end if

      !--- determine phase for this algorithm 
      cloudphase = ICE_PHASE
      if ((cloudtype == sym%FOG_TYPE) .or.  &
          (cloudtype == sym%WATER_TYPE) .or. & 
          (cloudtype == sym%SUPERCOOLED_TYPE)) then
          cloudphase = WATER_PHASE
      end if
      if ((cloudtype == sym%MIXED_TYPE)) then
          cloudphase = MIXED_PHASE
      endif

!----------------------------------------------------------------------
!--- Define the Observations
!----------------------------------------------------------------------

   !--- define y-vector - the observation vector
   if (imode == 1) then  ! 0.6 and 1.6 micron approach
      y(1) = Refl_065um/100.0
      y(2) = Refl_160um/100.0
   end if
   if (imode == 3) then  ! 0.6 and 3.9 micron based approach)
      y(1) = Refl_065um/100.0
      y(2) = Rad_375um * Rad_to_Ref_Fac_375um
   end if
   if (imode == 4) then  ! 0.6 and 3.9 micron pseudo reflectance based approach)
      y(1) = Refl_065um/100.0
      y(2) = Refl_375um/100.0
   end if

   !------------ measurement and forward model error
   Sy = 0.0
   Sy_inv = 0.0
   if (cloudphase == WATER_PHASE) then
     fm_err = 0.05      !fm model error
   else
     fm_err = 0.10
   end if

   Sy(1,1) = ((cal_err + fm_err + pp_err*Ref_var)*y(1) + 0.02)**2
   Sy(2,2) = ((cal_err + fm_err + pp_err*Ref_var)*y(2) + 0.02)**2

   if (imode == 3) then
     Sy(2,2) = ((cal_err + 1.0*fm_err + pp_err*Ref_var)*y(2) + 0.02)**2
   endif

   if (imode == 4) then
     Sy(2,2) = ((cal_err + 1.0*fm_err + pp_err*Ref_var)*y(2) + 0.02)**2
   endif


   !--- Invert Sy
   call INVERT_2x2(Sy,Sy_inv,isingular)
   if (isingular == sym%YES) then
     ifail = sym%YES
     Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 1
     Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 1
     Tau_Dcomp(Elem_Idx,Line_Idx) = X_Ap(1)
     Reff_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
     print *, "singular Sy"
     exit
   end if

   !--------------------------------------------------------------------------
   ! compute levels based on current estimate of cloud height
   !--------------------------------------------------------------------------

   !-- for nwp arrays
   call KNOWING_Z_COMPUTE_T_P_NWP(xnwp,ynwp,R4_Dummy,R4_Dummy,Z_cloud,ilev_nwp)  

   !-- for rtm arrays - note T_cloud and P_cloud are recomputed
   call KNOWING_Z_COMPUTE_T_P_RTM(xnwp,ynwp,P_cloud,R4_Dummy,Z_cloud,ilev_rtm)  

   !----------------------------------------------------------------------
   ! compute needed transmission terms
   !----------------------------------------------------------------------

     !--- determine amount of water vapor above the cloud
     prof_Weight = (Z_cloud - zlev(ilev_nwp)) /  (zlev(ilev_nwp+1)-zlev(ilev_nwp))

     Tpw_ac = Tpwlev(ilev_nwp) + &
              prof_Weight*(Tpwlev(ilev_nwp+1) - Tpwlev(ilev_nwp))

     !--- compute above cloud transmission for solar channels (h2o only)

     !-- above cloud ch1
     Tau_Aer =  Solar_Rtm%Tau_Aer(1) * (P_cloud / Psfc_Nwp(xnwp,ynwp))**2
     Tau_Ray =  Solar_Rtm%Tau_Ray(1) * (P_cloud / psfc_nwp(xnwp,ynwp))
     Tau_gas = Solar_Rtm%Tau_H2O_Coef(1,1) + Solar_Rtm%Tau_H2O_Coef(1,2)*Tpw_ac + Solar_Rtm%Tau_H2O_Coef(1,3)*(Tpw_ac**2)  + &
               Solar_Rtm%Tau_O3(1) + Solar_Rtm%Tau_O2(1) + Solar_Rtm%Tau_CH4(1) + Solar_Rtm%Tau_CO2(1)

     Trans_Ac_065um = exp(-1.0*Air_Mass*(Tau_Aer + Tau_Gas + Tau_Ray))
     Trans_Ac_065um = max(1.0e-06,min(1.0,Trans_Ac_065um))

     !-- total sky ch1
     Tau_Aer =  Solar_Rtm%Tau_Aer(1)
     Tau_Ray =  Solar_Rtm%Tau_Ray(1)
     Tau_Gas = Solar_Rtm%Tau_H2O_Coef(1,1) + Solar_Rtm%Tau_H2O_Coef(1,2)*Tpw + Solar_Rtm%Tau_H2O_Coef(1,3)*(Tpw**2)  + &
               Solar_Rtm%Tau_O3(1) + Solar_Rtm%Tau_O2(1) + Solar_Rtm%Tau_CH4(1) + Solar_Rtm%Tau_CO2(1)

     Trans_Sfc_065um = exp(-1.0*Air_Mass*(Tau_Aer + Tau_Gas + Tau_Ray))
     Trans_Sfc_065um = max(1.0e-06,min(1.0,Trans_Sfc_065um))


     !--- channel 6  - 1.6 micron

     !-- above cloud
     Tau_Aer =  Solar_Rtm%Tau_Aer(6) * (P_cloud / Psfc_Nwp(xnwp,ynwp))**2
     Tau_Ray =  Solar_Rtm%Tau_Ray(6) * (P_cloud / Psfc_Nwp(xnwp,ynwp))
     Tau_Gas = Solar_Rtm%Tau_H2O_Coef(6,1) + Solar_Rtm%Tau_H2O_Coef(6,2)*Tpw_Ac + Solar_Rtm%Tau_H2O_Coef(6,3)*(Tpw_Ac**2)  + &
               Solar_Rtm%Tau_O3(6) + Solar_Rtm%Tau_O2(6) + Solar_Rtm%Tau_CH4(6) + Solar_Rtm%Tau_CO2(6)

     Trans_Ac_160um = exp(-1.0*Air_Mass*(Tau_Aer + Tau_Gas + Tau_Ray))
     Trans_Ac_160um = max(1.0e-06,min(1.0,Trans_Ac_065um))

     !-- total sky
     Tau_Aer =  Solar_Rtm%Tau_Aer(6) 
     Tau_Ray =  Solar_Rtm%Tau_Ray(6)
     Tau_Gas = Solar_Rtm%Tau_H2O_Coef(6,1) + Solar_Rtm%Tau_H2O_Coef(6,2)*Tpw + Solar_Rtm%Tau_H2O_Coef(6,3)*(Tpw**2)  + &
               Solar_Rtm%Tau_O3(6) + Solar_Rtm%Tau_O2(6) + Solar_Rtm%Tau_CH4(6) + Solar_Rtm%Tau_CO2(6)

     Trans_Sfc_160um = exp(-1.0*Air_Mass*(Tau_Aer + Tau_Gas + Tau_Ray))
     Trans_Sfc_160um = max(1.0e-06,min(1.0,Trans_Sfc_065um))

     !--- channel 20 - 3.75/3.9 micron

     !--- use rtm structures
     Trans_Ac_375um = 1.0
     Trans_Sfc_375um = 1.0
     if ((imode == 3) .or. (imode == 4)) then
          Trans_Ac_375um_Zen =  rtm(xnwp,ynwp)%d(ivza)%ch(20)%Trans_Atm_Profile(ilev_rtm) +    &
               Prof_Weight*(rtm(xnwp,ynwp)%d(ivza)%ch(20)%Trans_Atm_Profile(ilev_rtm+1) -  &
                            rtm(xnwp,ynwp)%d(ivza)%ch(20)%Trans_Atm_Profile(ilev_rtm))
          Trans_Ac_375um = Trans_Ac_375um_Zen**(1.0/Cos_Sol_Zen)                    !now is the full path
          Trans_Sfc_375um = Trans_Atm_Ch20_Solar_Rtm(Elem_Idx,Line_Idx)

          Rad_Ac_375um_Zen =  rtm(xnwp,ynwp)%d(ivza)%ch(20)%Rad_Atm_Profile(ilev_rtm) +    &
               Prof_Weight*(rtm(xnwp,ynwp)%d(ivza)%ch(20)%Rad_Atm_Profile(ilev_rtm+1) -  &
                            rtm(xnwp,ynwp)%d(ivza)%ch(20)%Rad_Atm_Profile(ilev_rtm))
     end if

!----------------------------------------------------------------------
! set surface reflectance (0-1)
!----------------------------------------------------------------------
Alb_Sfc_065um = ch1_Sfc_alb_umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx))
Alb_Sfc_160um = ch6_Sfc_alb_umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx))
Alb_Sfc_375um = ch20_Sfc_alb_umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx))

!--- use value from clear composite if available over land
 if (Sfc%Sfc_Type(Elem_Idx,Line_Idx) /= sym%WATER_SFC) then

  if (Modis_Clr_Alb_Flag == sym%YES) then
   if (Sensor%Chan_On_Flag_Default(1) == sym%YES) then
    if (ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) > 0.0)  then
       Alb_Sfc_065um = ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)/100.0
    endif
   endif
   if (Sensor%Chan_On_Flag_Default(6) == sym%YES) then
    if (ch(6)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) > 0.0)  then
       Alb_Sfc_160um =ch(6)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)/100.0
    endif
   endif
  endif
 endif

!--- use sfc emissivitys for surface reflectance for 3.75 micron
 if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then
  if ((Emiss_Sfc_375um > 0.0) .and. (Sfc%Sfc_Type(Elem_Idx,Line_Idx) > 0))  then   
            Alb_Sfc_375um = 1.0 - Emiss_Sfc_375um
  endif  
 endif  

!------------------------------------------------------------------------------
!   estimate clear-sky toa ch1 reflectance
!------------------------------------------------------------------------------

    !--- compute rayleigh and aerosol contribution to 0.63 micron signal
     Tau_gas = Solar_Rtm%Tau_H2O_Coef(1,1) +  &
               Solar_Rtm%Tau_H2O_Coef(1,2)*Tpw_Ac +  &
               Solar_Rtm%Tau_H2O_Coef(1,3)*(Tpw_Ac**2)  + &
               Solar_Rtm%Tau_O3(1) +  &
               Solar_Rtm%Tau_O2(1) +  &
               Solar_Rtm%Tau_CH4(1) +  &
               Solar_Rtm%Tau_CO2(1)

     Tau_Gas = max(0.0,Tau_Gas)
     Tau_Aer =  Solar_Rtm%Tau_Aer(1) * (P_cloud / Psfc_Nwp(xnwp,ynwp))**2
     Tau_Ray =  Solar_Rtm%Tau_Ray(1) * (P_cloud / Psfc_Nwp(xnwp,ynwp))

     Cloud_Albedo_Ap_View = 0.0
     Cloud_Albedo_Ap_Sun = 0.0

     call COMPUTE_CLEAR_SKY_SCATTER(Tau_Aer, &
                                    Solar_Rtm%Wo_Aer(1), &
                                    Solar_Rtm%G_Aer(1), &
                                    Tau_Ray, &
                                    Tau_Gas, &
                                    Geo%Scatangle(Elem_Idx,Line_Idx), &
                                    Geo%Coszen(Elem_Idx,Line_Idx), &
                                    Geo%CosSolzen(Elem_Idx,Line_Idx),  &
                                    Cloud_Albedo_Ap_View, &
                                    Cloud_Albedo_Ap_Sun, &
                                    Ref_Ss_065um)

     Ref_Ss_065um = Ref_Ss_065um / 100.0   !convert from % to 0-1.0

     !------------------------------------------------------------------------------
     !--- Define the a priori
     !------------------------------------------------------------------------------

      !--- particle size apriori
      if (cloudphase == WATER_PHASE) then
        x_ap(2) = 1.0
      elseif (cloudphase == MIXED_PHASE) then
        x_ap(2) = 1.0
      elseif (cloudphase == ICE_PHASE) then
        x_ap(2) = 1.3
      end if

     !--- particle size apriori - fixed re retrieval

      x_ap(1) = 1.0   !can be anything here

      call LOCATE_TABLE_POSITION(Algo_Idx,Chan_Idx_065um,cloudphase,x_ap(1),x_ap(2),sol_zen,zen,az, &
                                 ir,is,it,it2,iu,iv,r,s,t,t2,u,v)

      call ESTIMATE_TAU_SINGLE_CHANNEL_FIXED_REFF(Algo_Idx,Chan_Idx_065um,cloudphase, &
                                             ir,it,it2,iu,iv, &
                                             Ref_Ss_065um,  &
                                             y(1), &
                                             Alb_Sfc_065um, &
                                             Trans_Ac_065um, &
                                             Trans_Sfc_065um, &
                                             x_ap(1), &
                                             cloud_albedo_ap_view, &
                                             cloud_albedo_ap_sun)

    !--- recompute clear sky scatter with cloud albedo terms
     call COMPUTE_CLEAR_SKY_SCATTER(Tau_aer, &
                                    Solar_Rtm%Wo_Aer(1), &
                                    Solar_Rtm%G_Aer(1), &
                                    Tau_ray, &
                                    Tau_gas, &
                                    Scat_Zen, &
                                    Cos_Zen, &
                                    Cos_Sol_Zen,  &
                                    Cloud_Albedo_Ap_View, &
                                    Cloud_Albedo_Ap_Sun, &
                                    Ref_Ss_065um)

      Ref_Ss_065um = Ref_Ss_065um / 100.0   !convert from % to 0-1.0

     !--- ignore single scattering reflectance in longer wavelength channels
     Ref_Ss_160um = 0.0
     Ref_Ss_375um = 0.0

     !--- compute final estimate of clear reflectance (0-1.0)
     Ref_2_clear = Ref_Ss_065um + Trans_Sfc_065um*Alb_Sfc_065um

      Tau_Dcomp_Ap(Elem_Idx,Line_Idx) = 10.0**x_ap(1)

      !---  a priori covariance matrix
      Sa = 0.0
      Sa(1,1) = (max(0.2,x_ap(1)))**2

      if (cloudphase == WATER_PHASE) then
        Sa(2,2) = (0.50)**2 
      else
        Sa(2,2) = (0.75)**2
      end if

      !--- Invert Sa
       call INVERT_2x2(Sa,Sa_inv,isingular)
       if (isingular == sym%YES) then
         ifail = sym%YES
         Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 1
         Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 1
         Tau_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
         Reff_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
         !print *, "singular matrix"
         EXIT
       end if

     !--- set first guess to a priori
       x = x_ap

if (diag_output == sym%YES) then
   print *, "<--- Begin New Retrieval for pixel = ", Elem_Idx,Line_Idx
   print *, "cloud type, phase = ", cloudtype, cloudphase
   print *, "angles = ", zen,sol_zen,az
   print *, "sfc ref = ", Alb_Sfc_065um,Alb_Sfc_160um,Alb_Sfc_375um
   print *, "Trans_Ac = ", Trans_Ac_065um,Trans_Ac_160um,Trans_Ac_375um
   print *, "y = ", y
   print *, "Sy = ", Sy
   print *, "x_ap = ",x_ap
   print *, "Sa = ",Sa
end if


!----------------------------------------------------------------------!
! call pixels darker than the clear reflectance missing
!----------------------------------------------------------------------
  if (y(1) < Ref_2_clear) then

      Tau_Dcomp(Elem_Idx,Line_Idx) = 0.0
      Cloud_063um_albedo(Elem_Idx,Line_Idx) = 0.0
      Cloud_063um_transmission_solar(Elem_Idx,Line_Idx) = 1.0
      Cloud_063um_transmission_view(Elem_Idx,Line_Idx) = 1.0
      cycle
  endif

!-----------------------------------------------------------------------
! begin retrieval loop
!-----------------------------------------------------------------------
iter = 0
conv_test = huge(conv_test)
ifail = sym%NO
retrieval_loop: do

 iter = iter + 1

 !--- check for convergence
 if ( (conv_test < Conv_Crit) .and. (iter > 1)) then
   if (diag_output == sym%YES) then
     print *,"converged after", iter, " iterations"
   end if
   ifail = sym%NO
   exit
 end if

 !--- check for convergence in terms of obs-fm
 if ( (iter > 2) .and. &
      (abs(y(1) - f(1)) < 0.5*cal_err*y(1)) .and. &
      (abs(y(2) - f(2)) < 0.5*cal_err*y(2))) then
   if (diag_output == sym%YES) then
     print *,"converged by obs-fm after", iter, " iterations"
   end if
   ifail = sym%NO
   exit
 end if

 !--- check for maximum number of iteration steps
 if (iter > iter_max) then
   if (diag_output == sym%YES) then
    print *,"failed convergence"
    STOP
   end if
   ifail = sym%YES
   Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 1
   Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 1
   Tau_Dcomp(Elem_Idx,Line_Idx) = 10.0 ** x_ap(1)
   Reff_Dcomp(Elem_Idx,Line_Idx) = 10.0 ** x_ap(2)
   exit
 end if


!--- some local variable for convience while iterating
log10tau = x(1)
log10re = x(2)

!--- interpolate within tables - assume the same for all channels and tables
 call LOCATE_TABLE_POSITION(Algo_Idx,Chan_Idx_065um,cloudphase,log10tau,log10re,sol_zen,zen,az, &
                            ir,is,it,it2,iu,iv,r,s,t,t2,u,v)

!--- call forward model for each channel used

    !--- 0.63 micron forward model
    call CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER(&
                         Algo_Idx,Chan_Idx_065um,cloudphase, &
                         r,s,t,t2,u,v,ir,is,it,it2,iu,iv, &
                         Trans_Ac_065um,Trans_Sfc_065um,Alb_Sfc_065um,Ref_Ss_065um, &
                         f(1),K(1,1),K(1,2))
   vis_Ref_fm(Elem_Idx,Line_Idx) = 100.0*f(1)    !store this as a diagnostic

   !--- 1.6 micron forward model
   if (imode == 1) then
    call CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER(&
                         Algo_Idx,Chan_Idx_160um,cloudphase, &
                         r,s,t,t2,u,v,ir,is,it,it2,iu,iv, &
                         Trans_Ac_160um,Trans_Sfc_160um,Alb_Sfc_160um, Ref_Ss_160um, &
                         f(2),K(2,1),K(2,2))
   endif

   !--- 3.75/3.9 micron forward model
   if ((imode == 3) .or. (imode == 4)) then

    !--- solar component
    call CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER(&
                         Algo_Idx,Chan_Idx_375um,cloudphase, &
                         r,s,t,t2,u,v,ir,is,it,it2,iu,iv, &
                         Trans_Ac_375um,Trans_Sfc_375um,Alb_Sfc_375um, Ref_Ss_375um, &
                         Ref_Sol, dRef_Sol_dtau, dRef_Sol_dRe)

    Ref_Sol_Toc = (Ref_Sol / Trans_Ac_375um_zen)

    !--- thermal component
    Rad_Therm = 0.0
    dRad_Therm_dtau = 0.0
    dRad_Therm_dRe = 0.0
    if (imode == 3) then
       call CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER(&
                         Algo_Idx,Chan_Idx_375um,cloudphase, &
                         r,s,u,ir,is,iu, &
                         Rad_375um_clr, T_cloud, Trans_Ac_375um_zen,Rad_Ac_375um_zen, &
                         Rad_Therm,dRad_Therm_dtau,dRad_Therm_dRe,Bt,dB_dT)

       Rad_Therm_Toc = (Rad_Therm - Rad_Ac_375um_zen) / Trans_Ac_375um_zen
    endif

    !--- combine solar and thermal
    f(2) = Ref_Sol + Rad_Therm * Rad_to_Ref_Fac_375um
    K(2,1) = dRef_Sol_dtau + dRad_Therm_dtau * Rad_to_Ref_Fac_375um
    K(2,2) = dRef_Sol_dRe + dRad_Therm_dRe * Rad_to_Ref_Fac_375um

    f_toc = Ref_Sol_Toc + Rad_Therm_Toc * Rad_to_Ref_Fac_375um

   endif

!--- Compute Sx
 Sx_inv = Sa_inv + matmul(transpose(K),matmul(Sy_inv,K)) !(Eq.102 Rodgers)
 call INVERT_2x2(Sx_inv,Sx,isingular)
 if (isingular == 1) then
   ifail = sym%YES
   Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 1
   Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 1
   Tau_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
   Reff_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
   exit
 end if

!--- compute next iteration step
 delta_x = matmul(Sx,(matmul(transpose(K),matmul(Sy_inv,(y-f))) +  &
                      matmul(Sa_inv,x_ap-x) ))

!--- check for convergence
  conv_test = abs(sum(delta_x*matmul(Sx_inv,delta_x)))

!--- control step size  (note change to preserve direction)
  delta_x_distance = sqrt( sum(delta_x**2))
  if (delta_x_distance > 0.0) then
   do iparam = 1,num_param
      delta_x_dir(iparam) = delta_x(iparam) / delta_x_distance
   ENDDO
   do iparam = 1,num_param
       delta_x_constrained(iparam) = sign( min(delta_x_max,abs(delta_x(iparam))) , delta_x(iparam) )
   end do
   delta_x_distance_constrained = sqrt( sum(delta_x_constrained**2))
   do iparam = 1,num_param
     delta_x(iparam) = delta_x_dir(iparam)*delta_x_distance_constrained
   end do
  endif
  

!--- update retrieved vector
  x = x + delta_x


if (diag_output  == sym%YES) then
  print *, "iter = ", iter
  print *, "f = ", f
  print *, "K = ", K
  print *, "Sx = ", Sx
  print *, "delta_x = ", delta_x
  print *, "new x = ", x
  print *, "conv test = ", conv_test
end if

end do retrieval_loop

!-----------------------------------------------------------------------
!--- compute albedo and transmission
!-----------------------------------------------------------------------
if (ifail == sym%NO) then

    !--- some local variable for convience while iterating
    log10tau = x(1)
    log10re = x(2)

    !--- interpolate within tables - assume the same for all channels and tables
    call LOCATE_TABLE_POSITION(Algo_Idx,Chan_Idx_065um,cloudphase,log10tau,log10re,sol_zen,zen,az, &
                            ir,is,it,it2,iu,iv,r,s,t,t2,u,v)

    call ESTIMATE_ALBEDO_TRANSMISSION_SINGLE_LAYER( &
                         Algo_Idx,cloudphase, &
                         r,s,t,t2,ir,is,it,it2, &
                         Cloud_063um_albedo(Elem_Idx,Line_Idx),  &
                         Cloud_063um_spherical_albedo(Elem_Idx,Line_Idx),  &
                         Cloud_063um_transmission_view(Elem_Idx,Line_Idx),  &
                         Cloud_063um_transmission_solar(Elem_Idx,Line_Idx))

 endif

!-----------------------------------------------------------------------
! end retrieval loop
!-----------------------------------------------------------------------

!--- delog for output and derive lwp and iwp
 if (ifail == sym%NO) then

  Tau_Dcomp(Elem_Idx,Line_Idx) = 10.0 ** x(1)
  Reff_Dcomp(Elem_Idx,Line_Idx) = 10.0 ** x(2)

 else  

   !--- for clear failed pixels, set to clear-sky values
   if ((Cloudmask == SYM%CLEAR) .or. (Cloudmask == SYM%PROB_CLEAR)) then
       Tau_Dcomp(Elem_Idx,Line_Idx) = 0.0
       Cloud_063um_Albedo(Elem_Idx,Line_Idx) = 0.0
       Cloud_063um_Transmission_Solar(Elem_Idx,Line_Idx) = 1.0
       Cloud_063um_Transmission_View(Elem_Idx,Line_Idx) = 1.0
   endif

 end if

 !------------------------------------------------------------------------------
 !--- determine quality flags for non-failed retrievals
 !------------------------------------------------------------------------------
 iparam = 1
 if (Sx(1,1) < Sa(1,1) / 4.0) then
  Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 3
 else
  Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 2
 endif

 if (Sx(2,2) < Sa(2,2) / 4.0) then
  Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 3
 else
   Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 2
 endif

 !--- set converged but undetected cloudy pixels as low quality
 if (undetected_cloud == sym%YES) then
       Tau_Dcomp_Qf(Elem_Idx,Line_Idx) = 2
       Reff_Dcomp_Qf(Elem_Idx,Line_Idx) = 2
 endif
 if (Tau_Dcomp_Qf(Elem_Idx,Line_Idx) == 1) then
     Tau_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
 endif
 if (Reff_Dcomp_Qf(Elem_Idx,Line_Idx) == 1) then
      Reff_Dcomp(Elem_Idx,Line_Idx) = Missing_Value_Real4
 endif


!-----------------------------------------------------------------------
! end loop over pixels in segment
!-----------------------------------------------------------------------
    end do element_loop_1
  end do line_loop_1

!======================================================================
! Nullify local pointers
!======================================================================
  if (iseg == nseg) then
   plev => null()
   tlev => null()
   zlev => null()
   Tpwlev => null()
  end if

!======================================================================
! Update previous Sc_Id name for orbit
!======================================================================
 !if (iseg == nseg) then
   Sc_Id_Prev = Sensor%WMO_Id
 !end if


 end subroutine CLOUD_TAU_RE_SOLAR

!======================================================================
! LOCATE POSITION IN TABLES
!
! This code finds position within the tables for the current conditions
!
! ir - is the index in re dimension such that value is between ir and ir+1
! is - is the index in tau dimension such that value is between is and is+1
! it - is the index in Solzen dimension such that value is between it and it+1
! it2 - is the index in zen / Solzen dimension such that value is between it2 and it2+1
! iu - is the index in zen dimension such that value is between iu and iu+1
! iv - is the index in Relaz dimension such that value is between iv and iv+1
!
! r - is the relative distance of value between ir and ir+1
!======================================================================
 subroutine LOCATE_TABLE_POSITION(Algo_Idx,Chan_Idx,Phase_Idx,tau,re,Solzen,zen,Relaz, &
                                   ir,is,it,it2,iu,iv,r,s,t,t2,u,v)

    integer(kind=int4), intent(in):: Algo_Idx,Chan_Idx
    integer(kind=int4), intent(in):: Phase_Idx
    real(kind=real4), intent(in):: tau,re,Solzen,zen,Relaz
    integer, intent(out):: ir,is,it,iu,iv,it2
    real(kind=real4), intent(out):: r,s,u,t,t2,v

!--- interpolate to find reflectance for this set of properties
   ir = max(1,min(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10re-1,          &
            1+int( (re - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(1)) / &
                (Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_log10re)  )))
   
   is = max(1,min(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau-1,        &
         1+int( (tau - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(1)) / &
                (Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_log10tau))))
 
   it = max(1,min(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nSolzen-1,       &
         1+int((Solzen - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(1)) / &
               (Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_Solzen))))

   it2 = max(1,min(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nSolzen-1,  &
          1+int((zen - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(1)) / &
                (Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_Solzen))))

   iu = max(1,min(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nzen-1, &
          1+int((zen - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(1)) / &
                (Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_zen))))

   !--- handle relative azimuth switching
   if (Relaz >= Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz_switch_pp) then

     iv = max(0,min(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nRelaz_in_pp-2,         &
              int((Relaz - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz_switch_pp) / &
               (Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_Relaz_in_pp))))  +        &
               Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nRelaz_out_pp + 1

   else

     iv = max(1,min(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nRelaz_out_pp,         &
              1+int((Relaz - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz(1)) / &
               (Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_Relaz_out_pp))))

   endif

!--- interpolation weights
   r = (re - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(ir)) / &
       ( Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(ir+1) - &
         Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(ir) )

   s = (tau - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(is)) / &
       ( Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(is+1) - &
         Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(is))

   t = (Solzen - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(it)) / &
       ( Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(it+1) - &
         Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(it))

   t2= (zen - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(it2)) / &
       ( Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(it2+1) - &
         Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(it2))

   u = (zen - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(iu)) / &
       ( Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(iu+1) - &
         Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(iu))

   v = (min(Relaz,360.0-Relaz) - Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz(iv)) / &
       ( Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz(iv+1) - &
         Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz(iv))

!  print *, "Relaz => ", Relaz, iv, Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nRelaz,   &
!                Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz_switch_pp
!  print *, "bound = ", Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz(iv:iv+1)
!  print *, "Relaz vector => ", Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz
!  print *, "v = ", v

 end subroutine LOCATE_TABLE_POSITION

!======================================================================
! interpolate within the reflectance lookup table once you
! determine the surrounding indices and weights
!
! r - weight used in re interpolation 
! s - weight used in tau interpolation 
! t - weight used in solar zenith interpolation 
! u - weight used in sensor zenith interpolation 
! v - weight used in relative azimuth interpolation 
!
! Ref_lut - subset of a reflectance table
! Tau_lut - subset of optical depth table
! Reff_lut - subset of effective radius table
!
! dimension of Ref_lut (r, s, t, u, v)
!
! temp_11 value at Tau_lut(is) and Reff_lut(ir)
! temp_21 value at Tau_lut(is+1) and Reff_lut(ir)
! temp_12 value at Tau_lut(is) and Reff_lut(ir+1)
! temp_22 value at Tau_lut(is+1) and Reff_lut(ir+1)
!
! note, dimensions of all arrays are 1:2
!======================================================================
 subroutine INTERPOLATE_REFLECTANCE_VALUE(r,s,t,u,v,t2, &
                Reff_lut,Tau_lut, &
                trn_lut,Trn2_lut,Sab_lut,Ref_lut, &
                ref,dRef_dtau,dRef_dRe, &
                trn,dTrn_dtau,dTrn_dRe, &
                Trn2,dTrn2_dtau,dTrn2_dRe, &
                 sab,dSab_dtau,dSab_dRe)
  real, intent(in):: r,s,t,u,v,t2
  real, dimension(2), intent(in):: Reff_lut, Tau_lut
  real, dimension(2,2,2), intent(in):: trn_lut,Trn2_lut
  real, dimension(2,2), intent(in):: Sab_lut
  real, dimension(2,2,2,2,2), intent(in):: Ref_lut
  real, intent(out):: ref, dRef_dtau, dRef_dRe, &
                      trn, dTrn_dtau, dTrn_dRe, &
                      Trn2, dTrn2_dtau, dTrn2_dRe, &
                      sab, dSab_dtau, dSab_dRe
  real:: temp_11, temp_12, temp_21, temp_22

   temp_11 = (1.0-t)*(1.0-u)*(1.0-v)*Ref_lut(1,1,1,1,1) + &
            (1.0-t)*(1.0-u)*(v)*Ref_lut(1,1,1,1,2) +   &
            (1.0-t)*(u)*(v)*Ref_lut(1,1,1,2,2) +     &
            (1.0-t)*(u)*(1.0-v)*Ref_lut(1,1,1,2,1) +   &
            (t)*(1.0-u)*(1.0-v)*Ref_lut(1,1,2,1,1) +   &
            (t)*(1.0-u)*(v)*Ref_lut(1,1,2,1,2) +     &
            (t)*(u)*(v)*Ref_lut(1,1,2,2,2) +       &
            (t)*(u)*(1.0-v)*Ref_lut(1,1,2,2,1)


   temp_21 = (1.0-t)*(1.0-u)*(1.0-v)*Ref_lut(1,2,1,1,1) + &
            (1.0-t)*(1.0-u)*(v)*Ref_lut(1,2,1,1,2) +   &
            (1.0-t)*(u)*(v)*Ref_lut(1,2,1,2,2) +     &
            (1.0-t)*(u)*(1.0-v)*Ref_lut(1,2,1,2,1) +   &
            (t)*(1.0-u)*(1.0-v)*Ref_lut(1,2,2,1,1) +   &
            (t)*(1.0-u)*(v)*Ref_lut(1,2,2,1,2) +     &
            (t)*(u)*(v)*Ref_lut(1,2,2,2,2) +       &
            (t)*(u)*(1.0-v)*Ref_lut(1,2,2,2,1)

   temp_12 = (1.0-t)*(1.0-u)*(1.0-v)*Ref_lut(2,1,1,1,1) + &
            (1.0-t)*(1.0-u)*(v)*Ref_lut(2,1,1,1,2) +   &
            (1.0-t)*(u)*(v)*Ref_lut(2,1,1,2,2) +     &
            (1.0-t)*(u)*(1.0-v)*Ref_lut(2,1,1,2,1) +   &
            (t)*(1.0-u)*(1.0-v)*Ref_lut(2,1,2,1,1) +   &
            (t)*(1.0-u)*(v)*Ref_lut(2,1,2,1,2) +     &
            (t)*(u)*(v)*Ref_lut(2,1,2,2,2) +       &
            (t)*(u)*(1.0-v)*Ref_lut(2,1,2,2,1)


   temp_22 = (1.0-t)*(1.0-u)*(1.0-v)*Ref_lut(2,2,1,1,1) + &
            (1.0-t)*(1.0-u)*(v)*Ref_lut(2,2,1,1,2) +   &
            (1.0-t)*(u)*(v)*Ref_lut(2,2,1,2,2) +     &
            (1.0-t)*(u)*(1.0-v)*Ref_lut(2,2,1,2,1) +   &
            (t)*(1.0-u)*(1.0-v)*Ref_lut(2,2,2,1,1) +   &
            (t)*(1.0-u)*(v)*Ref_lut(2,2,2,1,2) +     &
            (t)*(u)*(v)*Ref_lut(2,2,2,2,2) +       &
            (t)*(u)*(1.0-v)*Ref_lut(2,2,2,2,1)


   ref = (1.0-r)*(1.0-s)*temp_11 + (r)*(1.0-s)*temp_12 + &
         (1.0-r)*(s)*temp_21 + (r)*(s)*temp_22

   dRef_dtau = ((1.0-r)*(temp_21-temp_11) + (r)*(temp_22-temp_12)) /  &
                (Tau_lut(2) - Tau_lut(1))


   dRef_dRe = ((1.0-s)*(temp_12-temp_11) + (s)*(temp_22-temp_21)) /  &
                (Reff_lut(2) - Reff_lut(1))

!--- interp flux transmission at solar angle
   temp_11 = (1.0-t)*Trn_Lut(1,1,1) + (t)*Trn_Lut(1,1,2)
   temp_12 = (1.0-t)*Trn_Lut(2,1,1) + (t)*Trn_Lut(2,1,2)
   temp_21 = (1.0-t)*Trn_Lut(1,2,1) + (t)*Trn_Lut(1,2,2)
   temp_22 = (1.0-t)*Trn_Lut(2,2,1) + (t)*Trn_Lut(2,2,2)

   trn = (1.0-r)*(1.0-s)*temp_11 + (r)*(1.0-s)*temp_12 + &
         (1.0-r)*(s)*temp_21 + (r)*(s)*temp_22

   dTrn_dtau = ((1.0-r)*(temp_21-temp_11) + (r)*(temp_22-temp_12)) /  &
                (Tau_lut(2) - Tau_lut(1))

   dTrn_dRe = ((1.0-s)*(temp_12-temp_11) + (s)*(temp_22-temp_21)) /  &
                (Reff_lut(2) - Reff_lut(1))

!--- interp flux transmission at sensor angle
   temp_11 = (1.0-t2)*Trn2_lut(1,1,1) + (t2)*Trn2_lut(1,1,2)
   temp_12 = (1.0-t2)*Trn2_lut(2,1,1) + (t2)*Trn2_lut(2,1,2)
   temp_21 = (1.0-t2)*Trn2_lut(1,2,1) + (t2)*Trn2_lut(1,2,2)
   temp_22 = (1.0-t2)*Trn2_lut(2,2,1) + (t2)*Trn2_lut(2,2,2)

   Trn2 = (1.0-r)*(1.0-s)*temp_11 + (r)*(1.0-s)*temp_12 + &
         (1.0-r)*(s)*temp_21 + (r)*(s)*temp_22

   dTrn2_dtau = ((1.0-r)*(temp_21-temp_11) + (r)*(temp_22-temp_12)) /  &
                (Tau_lut(2) - Tau_lut(1))

   dTrn2_dRe = ((1.0-s)*(temp_12-temp_11) + (s)*(temp_22-temp_21)) /  &
                (Reff_lut(2) - Reff_lut(1))


!--- spherical albedo
  sab = (1.0-r)*(1.0-s)*Sab_lut(1,1) + (r)*(1.0-s)*Sab_lut(2,1) + &
        (1.0-r)*(s)*Sab_lut(1,2) + (r)*(s)*Sab_lut(2,2)

  dSab_dtau = ((1.0-r)*(Sab_lut(1,2)-Sab_lut(1,1)) +  &
               (r)*(Sab_lut(2,2)-Sab_lut(2,1))) /  &
                (Tau_lut(2) - Tau_lut(1))

  dSab_dRe = ((1.0-s)*(Sab_lut(2,1)-Sab_lut(1,1)) +  &
               (s)*(Sab_lut(2,2)-Sab_lut(1,2))) /  &
                (Reff_lut(2) - Reff_lut(1))


end subroutine INTERPOLATE_REFLECTANCE_VALUE

!=======================================================================
! interpolate within the emissivity lookup table once you
! determine the surrounding indices and weights
!
! r - weight used in re interpolation 
! s - weight used in tau interpolation 
! u - weight used in sensor zenith interpolation 
!  
! ems_lut - subset of a emissivity table
! trn_lut - subset of a emissivity table
! Tau_lut - subset of optical depth table
! Reff_lut - subset of effective radius table
!  
! dimension of ems_lut (r, s, u)
!
! temp_11 value at Tau_lut(is) and Reff_lut(ir)
! temp_21 value at Tau_lut(is+1) and Reff_lut(ir)
! temp_12 value at Tau_lut(is) and Reff_lut(ir+1)
! temp_22 value at Tau_lut(is+1) and Reff_lut(ir+1)
!
! note, dimensions of all arrays are 1:2
!======================================================================
 subroutine INTERPOLATE_EMISSIVITY_VALUE(r,s,u, &
                Reff_Lut,Tau_Lut, &
                Ems_Lut,Trn_Lut, &
                Ems,dEms_dtau,dEms_dRe, &
                Trn,dTrn_dtau,dTrn_dRe)
  real, intent(in):: r,s,u
  real, dimension(2), intent(in):: Reff_Lut, Tau_Lut
  real, dimension(2,2,2), intent(in):: Ems_Lut
  real, dimension(2,2,2), intent(in):: Trn_Lut
  real, intent(out):: Ems, dEms_dtau, dEms_dRe, &
                      Trn, dTrn_dtau, dTrn_dRe
  real:: temp_11, temp_12, temp_21, temp_22


! compute cloud emissivity
   temp_11 =  (1.0-u)*Ems_Lut(1,1,1) + (u)*Ems_Lut(1,1,2)
   temp_21 =  (1.0-u)*Ems_Lut(1,2,1) + (u)*Ems_Lut(1,2,2)
   temp_22 =  (1.0-u)*Ems_Lut(2,2,1) + (u)*Ems_Lut(2,2,2)
   temp_12 =  (1.0-u)*Ems_Lut(2,1,1) + (u)*Ems_Lut(2,1,2)

   ems = (1.0-r)*(1.0-s)*temp_11 + (r)*(1.0-s)*temp_12 + &
         (1.0-r)*(s)*temp_21 + (r)*(s)*temp_22

   ems = max(0.0,min(1.0,ems))  !needed?

   dEms_dtau =  ((1.0-r)*(temp_21-temp_11) + (r)*(temp_22-temp_12)) /  &
                (Tau_lut(2) - Tau_lut(1))

   dEms_dRe = ((1.0-s)*(temp_12-temp_11) + (s)*(temp_22-temp_21)) /  &
                (Reff_lut(2) - Reff_lut(1))

! compute cloud transmission
   temp_11 =  (1.0-u)*Trn_Lut(1,1,1) + (u)*Trn_Lut(1,1,2)
   temp_21 =  (1.0-u)*Trn_Lut(1,2,1) + (u)*Trn_Lut(1,2,2)
   temp_22 =  (1.0-u)*Trn_Lut(2,2,1) + (u)*Trn_Lut(2,2,2)
   temp_12 =  (1.0-u)*Trn_Lut(2,1,1) + (u)*Trn_Lut(2,1,2)

   trn = (1.0-r)*(1.0-s)*temp_11 + (r)*(1.0-s)*temp_12 + &
         (1.0-r)*(s)*temp_21 + (r)*(s)*temp_22

   trn = max(0.0,min(1.0,trn))   !needed?

   dTrn_dtau =  ((1.0-r)*(temp_21-temp_11) + (r)*(temp_22-temp_12)) /  &
                (Tau_lut(2) - Tau_lut(1))

   dTrn_dRe = ((1.0-s)*(temp_12-temp_11) + (s)*(temp_22-temp_21)) /  &
                (Reff_lut(2) - Reff_lut(1))

end subroutine INTERPOLATE_EMISSIVITY_VALUE


!======================================================================
! Forward Model for a Reflectance Channel
!
! this routine is appropriate for single layer cloud placed in an 
! absorbing atmosphere over a Lambertian Surface.  The accounting
! for the surface contribution is given by Eq 11 in 
! King, Michael, 1987: JAS, vol44, 1734-1751 
!======================================================================
subroutine CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER(&
                         Algo_Idx,Chan_Idx,Phase_Idx, &
                         r,s,t,t2,u,v,ir,is,it,it2,iu,iv, &
                         Trans_Ac,Trans_Sfc, Alb_Sfc, Ref_Ss, &
                         ref,dRef_dtau,dRef_dRe)

  integer, intent(in):: Algo_Idx,Chan_Idx,       &
                        ir,is,it,it2,iu,iv
  integer(kind=int4), intent(in):: Phase_Idx
  real, intent(in):: r,s,t,t2,u,v
  real, intent(in):: Trans_Ac
  real, intent(in):: Trans_Sfc
  real, intent(in):: Alb_Sfc
  real, intent(in):: Ref_Ss
  real, intent(out):: ref, dRef_dtau, dRef_dRe


  real:: trn, dTrn_dtau, dTrn_dRe
  real:: Trn2, dTrn2_dtau, dTrn2_dRe
  real:: sab,dSab_dtau,dSab_dRe,temp,Trans_Bc,Alb_Sfc_Temp

  real, dimension(2):: Tau_lut_sub
  real, dimension(2):: Reff_lut_sub
  real, dimension(2,2):: Sab_lut_sub
  real, dimension(2,2,2):: trn_lut_sub
  real, dimension(2,2,2):: Trn2_lut_sub
  real, dimension(2,2,2,2,2):: Ref_lut_sub

!----------------------------------------------------------------------
! construct sub-tables from appropriate lookup tables
!----------------------------------------------------------------------
  Reff_lut_sub = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(ir:ir+1)
  Tau_lut_sub = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(is:is+1)
  Ref_lut_sub = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ref(ir:ir+1,is:is+1,it:it+1,iu:iu+1,iv:iv+1)
  trn_lut_sub = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn(ir:ir+1,is:is+1,it:it+1)
  Trn2_lut_sub = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn(ir:ir+1,is:is+1,it2:it2+1)
  Sab_lut_sub = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%sph_alb(ir:ir+1,is:is+1)

!--- Interpolate the sub-tables
 call INTERPOLATE_REFLECTANCE_VALUE(r,s,t,u,v,t2, &
                Reff_lut_sub,Tau_lut_sub, &
                trn_lut_sub,Trn2_lut_sub,Sab_lut_sub,Ref_lut_sub, &
                ref,dRef_dtau,dRef_dRe, &
                trn,dTrn_dtau,dTrn_dRe, &
                Trn2,dTrn2_dtau,dTrn2_dRe, &
                sab,dSab_dtau,dSab_dRe)

!--- Estimate transmission from cloud to surface
 Trans_Bc = 0.0
 if (Trans_Ac > 0.0) then
  Trans_Bc = Trans_Sfc / Trans_Ac
 end if

 !--- Account for surface reflection and atmospheric transmission
 Alb_Sfc_Temp = Trans_Bc * Alb_Sfc
 temp = Alb_Sfc_Temp / (1.0 - Alb_Sfc_Temp * sab)
 ref = Trans_Ac * (ref + temp * trn*Trn2) + Ref_Ss

!--- this is wrong
! dRef_dtau = Trans_Ac * (dRef_dtau  + temp*trn*dTrn2_dtau + temp*Trn2*dTrn_dtau + &
!                            (temp*trn)**2*dSab_dtau)
!----

 dRef_dtau = Trans_Ac * (dRef_dtau  + &
                         temp*trn*dTrn2_dtau + &
                         temp*Trn2*dTrn_dtau + &
                         (trn*Trn2)*((Alb_Sfc_Temp**2)*dSab_dtau))/((1-Alb_Sfc_Temp*sab)**2)

 dRef_dRe = Trans_Ac *(dRef_dRe +  &
                       temp*trn*dTrn2_dRe +  &
                       temp*Trn2*dTrn_dRe + &
                       (trn*Trn2)*((Alb_Sfc_Temp**2)*dSab_dRe))/((1-Alb_Sfc_Temp*sab)**2)


end subroutine CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER

!======================================================================
! Forward Model for a Thermal Channel
!
! this routine is appropriate for single layer cloud placed in an 
! absorbing atmosphere over a Lambertian Surface. 
!======================================================================
subroutine CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER(&
                         Algo_Idx,Chan_Idx,Phase_Idx, &
                         r,s,u,ir,is,iu, &
                         Rad_clr,Tc,Trans_Ac,Rad_Ac,  &
                         rad,dRad_dtau,dRad_dRe,Bt,dB_dT)

  integer, intent(in):: Algo_Idx
  integer, intent(in):: Chan_Idx
  integer, intent(in):: ir
  integer, intent(in):: is
  integer, intent(in):: iu
  integer(kind=int4), intent(in):: Phase_Idx
  real, intent(in):: r,s,u
  real, intent(in):: Rad_clr
  real, intent(in):: Tc
  real, intent(in):: Trans_Ac
  real, intent(in):: Rad_Ac
  real, intent(out):: rad
  real, intent(out):: dRad_dtau
  real, intent(out):: dRad_dRe
  real, intent(out):: Bt
  real, intent(out):: dB_dT


  real:: ems
  real:: dEms_dtau
  real:: dEms_dRe
  real:: Trn
  real:: dTrn_dTau
  real:: dTrn_dRe
  real:: Bc

  real, dimension(2):: Tau_lut_sub
  real, dimension(2):: Reff_lut_sub
  real, dimension(2,2,2):: trn_lut_sub
  real, dimension(2,2,2):: ems_lut_sub

!----------------------------------------------------------------------
! construct sub-tables from appropriate lookup tables
!----------------------------------------------------------------------
  Reff_lut_sub = Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(ir:ir+1)
  Tau_lut_sub = Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(is:is+1)
  ems_lut_sub = Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%emiss(ir:ir+1,is:is+1,iu:iu+1)
  trn_lut_sub = Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trans(ir:ir+1,is:is+1,iu:iu+1)

!--- Interpolate the sub-tables
 call INTERPOLATE_EMISSIVITY_VALUE(r,s,u, &
                Reff_lut_sub,Tau_lut_sub, &
                ems_lut_sub,trn_lut_sub, &
                ems,dEms_dtau,dEms_dRe, &
                trn,dTrn_dtau,dTrn_dRe)

!--- compute planck emission at cloud temperature
 Bc  = planck_Rad_fast(Chan_Idx,Tc)

!--- compute toa radiance and its Jacobian
 Rad = Ems*Rad_Ac + Trans_Ac * Ems * Bc +  trn * Rad_clr


 dRad_dtau =  -Rad_Ac*dTrn_dtau + Trans_Ac*Bc*dEms_dtau + &
               Rad_clr*dTrn_dtau

 dRad_dRe = -Rad_Ac*dTrn_dRe + Trans_Ac*Bc*dEms_dRe + &
               Rad_clr*dTrn_dRe 

 Bt = planck_Temp_fast(Chan_Idx,rad,dB_dT = dB_dT)

end subroutine CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER

!======================================================================
! Estimate flux terms for a Reflectance Channel
!
!======================================================================
subroutine ESTIMATE_ALBEDO_TRANSMISSION_SINGLE_LAYER( &
                         Algo_Idx,Phase_Idx, &
                         r,s,t,t2,ir,is,it,it2, &
                         cloud_albedo,  &
                         cloud_spherical_albedo, &
                         cloud_transmission_view,  &
                         cloud_transmission_solar)

  integer, intent(in):: Algo_Idx,       &
                        ir,is,it,it2
  integer(kind=int4), intent(in):: Phase_Idx
  real, intent(in):: r,s,t,t2
  real, intent(out):: cloud_albedo
  real, intent(out):: cloud_spherical_albedo
  real, intent(out):: cloud_transmission_view
  real, intent(out):: cloud_transmission_solar

  real, dimension(2,2):: sph_alb_lut
  real, dimension(2,2,2):: alb_lut
  real, dimension(2,2,2):: trn_lut
  real, dimension(2,2,2):: Trn2_lut
  real:: temp_11
  real:: temp_12
  real:: temp_21
  real:: temp_22
  integer:: Chan_Idx

  !----------------------------------------------------------------------
  ! construct sub-tables from appropriate lookup tables
  !----------------------------------------------------------------------
  Chan_Idx = Chan_Idx_065um
  trn_lut = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn(ir:ir+1,is:is+1,it:it+1)
  Trn2_lut = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn(ir:ir+1,is:is+1,it2:it2+1)
  alb_lut = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%alb(ir:ir+1,is:is+1,it:it+1)
  sph_alb_lut = Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%sph_alb(ir:ir+1,is:is+1)

  !--- interp flux transmission at solar angle
  temp_11 = (1.0-t)*Trn_Lut(1,1,1) + (t)*Trn_Lut(1,1,2)
  temp_12 = (1.0-t)*Trn_Lut(2,1,1) + (t)*Trn_Lut(2,1,2)
  temp_21 = (1.0-t)*Trn_Lut(1,2,1) + (t)*Trn_Lut(1,2,2)
  temp_22 = (1.0-t)*Trn_Lut(2,2,1) + (t)*Trn_Lut(2,2,2)

  cloud_transmission_solar =  &
         (1.0-r)*(1.0-s)*temp_11 + (r)*(1.0-s)*temp_12 + &
         (1.0-r)*(s)*temp_21 + (r)*(s)*temp_22

  !--- interp flux transmission at view angle
  temp_11 = (1.0-t2)*Trn2_lut(1,1,1) + (t2)*Trn2_lut(1,1,2)
  temp_12 = (1.0-t2)*Trn2_lut(2,1,1) + (t2)*Trn2_lut(2,1,2)
  temp_21 = (1.0-t2)*Trn2_lut(1,2,1) + (t2)*Trn2_lut(1,2,2)
  temp_22 = (1.0-t2)*Trn2_lut(2,2,1) + (t2)*Trn2_lut(2,2,2)

  cloud_transmission_view =  &
         (1.0-r)*(1.0-s)*temp_11 + (r)*(1.0-s)*temp_12 + &
         (1.0-r)*(s)*temp_21 + (r)*(s)*temp_22

  !--- interp albedo at solar angle
  temp_11 = (1.0-t)*alb_lut(1,1,1) + (t)*alb_lut(1,1,2)
  temp_12 = (1.0-t)*alb_lut(2,1,1) + (t)*alb_lut(2,1,2)
  temp_21 = (1.0-t)*alb_lut(1,2,1) + (t)*alb_lut(1,2,2)
  temp_22 = (1.0-t)*alb_lut(2,2,1) + (t)*alb_lut(2,2,2)

  cloud_albedo =  &
         (1.0-r)*(1.0-s)*temp_11 + (r)*(1.0-s)*temp_12 + &
         (1.0-r)*(s)*temp_21 + (r)*(s)*temp_22

  !--- spherical albedo
  cloud_spherical_albedo  =  &
        (1.0-r)*(1.0-s)*sph_alb_lut(1,1) + (r)*(1.0-s)*sph_alb_lut(2,1) + &
        (1.0-r)*(s)*sph_alb_lut(1,2) + (r)*(s)*sph_alb_lut(2,2)

end subroutine ESTIMATE_ALBEDO_TRANSMISSION_SINGLE_LAYER

!=======================================================================
!=======================================================================
 subroutine ESTIMATE_TAU_SINGLE_CHANNEL_FIXED_REFF(Algo_Num,Chan_Idx,Phase_Idx, &
                                             ir,it,it2,iu,iv, &
                                             Ref_Ss,  &
                                             Ref_obs,  &
                                             Alb_Sfc, &
                                             Trans_Ac, &
                                             Trans_Sfc, &
                                             log10_Tau_apriori, &
                                             albedo_apriori_view, &
                                             albedo_apriori_sun)


  integer, intent(in):: Algo_Num
  integer, intent(in):: Chan_Idx
  integer, intent(in):: Phase_Idx
  integer, intent(in):: ir
  integer, intent(in):: it
  integer, intent(in):: it2
  integer, intent(in):: iu
  integer, intent(in):: iv
  real, intent(in):: Ref_Ss
  real, intent(in):: Ref_obs
  real, intent(in):: Alb_Sfc
  real, intent(in):: Trans_Ac
  real, intent(in):: Trans_Sfc
  real, intent(out):: log10_Tau_apriori
  real, intent(out):: albedo_apriori_view
  real, intent(out):: albedo_apriori_sun

  real:: Trans_Bc
  real:: Alb_Sfc_Temp
! integer, dimension(1):: Tau_Vector_index
  integer:: Tau_index, n_tau
  real:: dTau_dRef
  real:: dRef
  real :: Tau_interp_Weight

  real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: Tau_Vector
  real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: Ref_Cld_Vector
  real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: trn_Cld_Vector
  real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: Trn2_Cld_Vector
  real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: Alb_Cld_Vector
  real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: Alb2_Cld_Vector
  real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: Sab_Vector
  real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: Ref_Toa_Vector
  real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: Temp_Vector
! real, dimension(Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau):: abs_diff_Vector

!----------------------------------------------------------------------
! construct sub-tables from appropriate lookup tables
!----------------------------------------------------------------------
  Tau_Vector = Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau
  Ref_Cld_Vector = Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%ref(ir,:,it,iu,iv)
  trn_Cld_Vector = Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%trn(ir,:,it)
  Trn2_Cld_Vector = Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%trn(ir,:,it2)
  Alb_Cld_Vector = Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%alb(ir,:,it)
  Alb2_Cld_Vector = Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%alb(ir,:,it2)
  Sab_Vector = Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%sph_alb(ir,:)


   !--- add in single scattering
     Ref_Cld_Vector = Ref_Cld_Vector + Ref_Ss

   !--- Estimate transmission from cloud to surface
     Trans_Bc = 0.0
     if (Trans_Ac > 0.0) then
      Trans_Bc = Trans_Sfc / Trans_Ac
     end if

   !--- Account for surface reflection and atmospheric transmission
     Alb_Sfc_Temp = Trans_Bc * Alb_Sfc
     Temp_Vector = Alb_Sfc_Temp / (1.0 - Alb_Sfc_Temp * Sab_Vector)
     Ref_Toa_Vector = Trans_Ac * (Ref_Cld_Vector + Temp_Vector * trn_Cld_Vector*Trn2_Cld_Vector)

   !--- find closest tau to match Ref_obs
     N_Tau = Cld_Ref(Algo_Num)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau
     call LOCATE(Ref_Toa_Vector, n_tau, Ref_obs, Tau_index)
     Tau_index = min(n_tau-1,max(1,Tau_index))

     dRef = Ref_Toa_Vector(Tau_index+1)-Ref_Toa_Vector(Tau_index)
     dTau_dRef = 0.0
     if (dRef > 0) then
       dTau_dRef = (Tau_Vector(Tau_index+1) - Tau_Vector(Tau_index))/dRef
     endif
     log10_Tau_apriori = Tau_Vector(Tau_index) + dTau_dRef * (Ref_obs - Ref_Toa_Vector(Tau_index))


    !--- get albedoes for for future use in rayleigh correction
     Tau_interp_Weight = (log10_Tau_apriori - Tau_Vector(Tau_index)) /  &
                         (Tau_Vector(Tau_index+1)-Tau_Vector(Tau_index))

     albedo_apriori_view = Alb_Cld_Vector(Tau_index) +  Tau_interp_Weight* &
                           (Alb_Cld_Vector(Tau_index+1) - Alb_Cld_Vector(Tau_index))

     albedo_apriori_sun = Alb2_Cld_Vector(Tau_index) +  Tau_interp_Weight* &
                           (Alb2_Cld_Vector(Tau_index+1) - Alb2_Cld_Vector(Tau_index))

    !--- maybe faster - use minloc
    !abs_diff_Vector = abs(Ref_Toa_Vector - Ref_obs)
    !Tau_Vector_index = MINLOC(abs_diff_Vector)
    !log10_Tau_apriori = Tau_Vector(Tau_Vector_index(1))

 end subroutine ESTIMATE_TAU_SINGLE_CHANNEL_FIXED_REFF

!=======================================================================
! Generic Routine to read the tables for a channel
!=======================================================================
subroutine READ_TABLES(Algo_Idx,Chan_Idx,lut_path)

 integer(kind=int4), intent(in):: Algo_Idx,Chan_Idx
 character(len=*), intent(in):: lut_path


      !--- read in reflectance tables for all applicable channels
      if (Chan_Idx <= 20) then

         !--- set flag to indicate that this channel has been read in 
         Cld_Ref(Algo_Idx)%channel(Chan_Idx)%flag = sym%YES

          !--- set flag and read in water phase tables
          Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(WATER_PHASE)%flag = sym%YES
          call READ_CLD_REF_LUT_STRUC(lut_path,Algo_Idx,Chan_Idx,WATER_PHASE)

          !--- set flag and read in ice phase tables
          Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(ICE_PHASE)%flag = sym%YES
          call READ_CLD_REF_LUT_STRUC(lut_path,Algo_Idx,Chan_Idx,ICE_PHASE)

          !--- set flag and read in mixed phase tables
          Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(MIXED_PHASE)%flag = sym%YES
          call READ_CLD_REF_LUT_STRUC(lut_path,Algo_Idx,Chan_Idx,MIXED_PHASE)

      end if


      !--- read in emissivity tables for all applicable channels
      if (Chan_Idx >= 20) then

         !--- set flag to indicate that this channel has been read in 
         Cld_Ems(Algo_Idx)%channel(Chan_Idx)%flag = sym%YES

          !--- set flag and read in water phase tables
          Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(WATER_PHASE)%flag = sym%YES
          call READ_CLD_EMS_LUT_STRUC(lut_path,Algo_Idx,Chan_Idx,WATER_PHASE)

          !--- set flag and read in ice phase tables
          Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(ICE_PHASE)%flag = sym%YES
          call READ_CLD_EMS_LUT_STRUC(lut_path,Algo_Idx,Chan_Idx,ICE_PHASE)

          !--- set flag and read in mixed phase tables
          Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(MIXED_PHASE)%flag = sym%YES
          call READ_CLD_EMS_LUT_STRUC(lut_path,Algo_Idx,Chan_Idx,MIXED_PHASE)

      end if


end subroutine READ_TABLES
!=======================================================================
! Read Water Cloud Lookup Table based on structures
! 
!
! isat  = wmo satellite identIFication number
! Chan_Idx = ABI channel number
!
! sc_ind = spacecraft index
! scinfo(sc_ind)%id = wmo number
!
! This does ice and water
!=======================================================================
 subroutine READ_CLD_REF_LUT_STRUC(Ancil_Data_Path,Algo_Idx,Chan_Idx,Phase_Idx)

  integer, intent(in):: Algo_Idx,Chan_Idx
  integer(kind=int1), intent(in):: Phase_Idx
  character (len=*), intent(in):: Ancil_Data_Path
 
  character (len=128):: Lut_File
  integer:: astatus
  integer:: tstatus
  character (len=2):: Chan_String
  character (len=3):: Phase_String
  character (len=7):: Sat_String
  integer:: nlog10re
  integer:: nlog10tau
  integer:: nSolzen
  integer:: nzen
  integer:: nRelaz           !total number relative azimuths
  integer:: nRelaz_in_pp     !number of Relaz in the principal plane zone
  integer:: nRelaz_out_pp    !number of Relaz outside the principal plane zone
  real:: delta_Relaz_in_pp   !azimuth space in principal plane zone
  real:: delta_Relaz_out_pp  !azimuth space outside principal plane zone
  real:: Relaz_switch_pp     !azimuth defining principal plane zone

  !-- hdf
  integer:: sfstart, sfrcatt,sfn2index, sffattr, sfendacc,sfend,sfrdata,sfselect
  integer:: sd_id
  integer:: Sds_Id
  integer:: istatus

  integer, parameter:: sds_rank_1d = 1
  integer, dimension(sds_rank_1d):: Sds_Start_1d, Sds_Edge_1d, Sds_Stride_1d

  integer, parameter:: sds_rank_2d = 2
  integer, dimension(sds_rank_2d):: Sds_Start_2d, Sds_Edge_2d, Sds_Stride_2d, sds_dims_2d

  integer, parameter:: sds_rank_3d = 3
  integer, dimension(sds_rank_3d):: Sds_Start_3d, Sds_Edge_3d, Sds_Stride_3d, sds_dims_3d

  integer, parameter:: sds_rank_5d = 5
  integer, dimension(sds_rank_5d):: Sds_Start_5d, Sds_Edge_5d, Sds_Stride_5d, sds_dims_5d

!----------------------------------------------------------------------------------
! Executable Code
!----------------------------------------------------------------------------------

!--- construct appropriate file name - use actual satellite channel, not abi equivalent
  Phase_String = "wat"
  if (Phase_Idx == ICE_PHASE) then
    Phase_String = "ice"
  end if
  if (Phase_Idx == MIXED_PHASE) then
    Phase_String = "wat"                  !note assumed water tables apply to mixed phase clouds
  end if

  if (index(Sensor%Sensor_Name,'AVHRR') > 0) then
   if (Chan_Idx == 1) then
          Chan_String = "1"
   elseif (Chan_Idx == 6) then
          Chan_String = "3a"
   elseif (Chan_Idx == 20) then
          Chan_String = "3b"
   endif
  endif
  if (index(Sensor%Sensor_Name,'MODIS') > 0) then
   if (Chan_Idx == 1) then
          Chan_String = "1"
   elseif (Chan_Idx == 6) then
          Chan_String = "6"
   elseif (Chan_Idx == 20) then
          Chan_String = "20"
   endif
  endif
  if (index(Sensor%Sensor_Name,'GOES') > 0) then
   if (Chan_Idx == 1) then
          Chan_String = "1"
   elseif (Chan_Idx == 6) then
          Chan_String = "0"
   elseif (Chan_Idx == 20) then
          Chan_String = "2"
   endif
 endif

   select case (Sensor%WMO_Id)

      case(3)
        Sat_String  = "METOP-A"
      case(4)
        Sat_String  = "METOP-B"
      case(5)
        Sat_String  = "METOP-C"
      case(252)
        Sat_String  = "GOES-08"
      case(253)
        Sat_String  = "GOES-09"
      case(254)
        Sat_String  = "GOES-10"
      case(255)
        Sat_String  = "GOES-11"
      case(256)
        Sat_String  = "GOES-12"
      case(257)
        Sat_String  = "GOES-13"
      case(258)
        Sat_String  = "GOES-14"
      case(259)
        Sat_String  = "GOES-15"
      case(706)
        Sat_String  = "NOAA-06"
      case(223)
        Sat_String  = "NOAA-19"
      case(783) 
        Sat_String = "MODIS"
      case(784)  
        Sat_String = "MODIS"
      case default
         print*,'dcomp test code not working for sensor wmo: ',Sensor%WMO_Id
         print*, 'ask Andy.....'
         stop 99
   end select

  if (trim(Chan_String) == "0") return

  Lut_File = TRIM(Sat_String)//"_ch"//TRIM(Chan_String)//"_ref_lut_"//Phase_String//"_cld.hdf"


 !--- Open lookup table

 !-- hdf
 sd_id = sfstart(TRIM(Ancil_Data_Path)//TRIM(Lut_File), DFACC_READ) 

 if (sd_id < 0) then
     print *, EXE_PROMPT, "Error open cloud ref table ", trim(Ancil_Data_Path)//trim(Lut_File)
     stop 99
  endif

!---- Read in tables

!--- global attributes 
   istatus = 0

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"HEADER"),  &
                     Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%header)

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_EFF_RADII"),  &
                     nlog10re)

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_OPTICAL_DEPTHS"),  &
                     nlog10tau)

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_SOLAR_ZENITH_ANGLES"),  &
                     nSolzen) + istatus

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_SENSOR_ZENITH_ANGLES"),  &
                     nzen) + istatus

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_RELATIVE_AZIMUTH_ANGLES"),  &
                     nRelaz) + istatus

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_RELATIVE_AZIMUTH_ANGLES_IN_PP"),  &
                     nRelaz_in_pp) + istatus

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_RELATIVE_AZIMUTH_ANGLES_OUTSIDE_PP"),  &
                     nRelaz_out_pp) + istatus

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"AZIMUTH_SPACING_IN_PP"),  &
                     delta_Relaz_in_pp) + istatus

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"AZIMUTH_SPACING_OUTSIDE_PP"),  &
                     delta_Relaz_out_pp) + istatus

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"AZIMUTH_SWITCH_POINT"),  &
                     Relaz_switch_pp) + istatus

   if (istatus /= 0) then
      print "(a,'Error reading file attributes')",EXE_PROMPT
      stop 99
   endif

   !--- store local dimensions into data structures
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10re = nlog10re
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau = nlog10tau
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nSolzen = nSolzen
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nzen = nzen
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nRelaz = nRelaz
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nRelaz_in_pp = nRelaz_in_pp
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nRelaz_out_pp = nRelaz_out_pp
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_Relaz_out_pp = delta_Relaz_out_pp
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_Relaz_in_pp = delta_Relaz_in_pp
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz_switch_pp = Relaz_switch_pp

   !--- knowing the dimensions, allocate space
   tstatus = 0
   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen)) then
      allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(nzen), stat=astatus)
      tstatus = tstatus + astatus
   endif
   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen)) then
     allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(nSolzen), stat=astatus)
     tstatus = tstatus + astatus
   endif

   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz)) then
      allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz(nRelaz), stat=astatus)
      tstatus = tstatus + astatus
   endif

   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re)) then
      allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(nlog10re), stat=astatus)
      tstatus = tstatus + astatus
   endif

   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau)) then
      allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(nlog10tau), stat=astatus)
      tstatus = tstatus + astatus
   endif

   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio)) then
      allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio(nlog10re), stat=astatus)
      tstatus = tstatus + astatus
   endif

   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref)) then
     allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref(nlog10re), stat=astatus)
     tstatus = tstatus + astatus
   endif

   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%sph_alb)) then
     allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%sph_alb(nlog10re,nlog10tau), stat=astatus)
      tstatus = tstatus + astatus
   endif

   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn)) then
     allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn(nlog10re,nlog10tau,nSolzen), stat=astatus)
     tstatus = tstatus + astatus
   endif

   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%alb)) then
     allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%alb(nlog10re,nlog10tau,nSolzen), stat=astatus)
     tstatus = tstatus + astatus
   endif

   if (.NOT. allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ref)) then
     allocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ref(nlog10re,nlog10tau,nSolzen,nzen,nRelaz), stat=astatus)
     tstatus = tstatus + astatus
   endif

   if (tstatus /= 0) then
      print "(a,'Error allocating cloud reflectance lookup table arrays, stopping')",EXE_PROMPT
      STOP 99
   end if

   !--- initialize
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen  = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen  = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz  = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re  = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau  = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio  = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref  = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%sph_alb  = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%alb = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ref = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%rho = 0.0
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%lambda_ref = 0.0


!--- Read in values in the table
istatus = 0

!-- read in vectors
Sds_Start_1d = 0
Sds_Stride_1d = 1
Sds_Edge_1d = nzen

!-- sensor zenith angle
Sds_Edge_1d = nzen

Sds_Id = sfselect(sd_id,sfn2index(sd_id,"sensor_zenith_angle"))
istatus = sfrdata(Sds_Id,Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen) +  &
                  istatus
istatus = sfendacc(Sds_Id) + istatus

!-- solar sensor zenith angle
Sds_Edge_1d = nSolzen
Sds_Id = sfselect(sd_id,sfn2index(sd_id,"solar_zenith_angle"))
istatus = sfrdata(Sds_Id,Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen) +  &
                  istatus
istatus = sfendacc(Sds_Id) + istatus

!-- relative azimith angle
Sds_Edge_1d = nRelaz
Sds_Id = sfselect(sd_id,sfn2index(sd_id,"relative_azimuth_angle"))
istatus = sfrdata(Sds_Id,Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz) +  &
                  istatus
istatus = sfendacc(Sds_Id) + istatus

!-- log10 optical depth
Sds_Edge_1d = nlog10tau
Sds_Id = sfselect(sd_id,sfn2index(sd_id,"log10_optical_depth"))
istatus = sfrdata(Sds_Id,Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau) +  &
                  istatus
istatus = sfendacc(Sds_Id) + istatus

!-- log10 effective radius
Sds_Edge_1d = nlog10re
Sds_Id = sfselect(sd_id,sfn2index(sd_id,"log10_eff_radius"))
istatus = sfrdata(Sds_Id,Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re) +  &
                  istatus
istatus = sfendacc(Sds_Id) + istatus

!-- Qext ref
Sds_Edge_1d = nlog10re
Sds_Id = sfselect(sd_id,sfn2index(sd_id,"Qext_ref"))
istatus = sfrdata(Sds_Id,Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref) +  &
                  istatus
istatus = sfendacc(Sds_Id) + istatus

!-- ext_ratio
Sds_Edge_1d = nlog10re
Sds_Id = sfselect(sd_id,sfn2index(sd_id,"extinction_ratio"))
istatus = sfrdata(Sds_Id,Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio) +  &
                  istatus
istatus = sfendacc(Sds_Id) + istatus

!--- 2d arrays
Sds_Start_2d = 0
Sds_Stride_2d = 1

!--- spherical albedo
sds_dims_2d(1) = nlog10re
sds_dims_2d(2) = nlog10tau
Sds_Edge_2d(1) = nlog10re
Sds_Edge_2d(2) = nlog10tau
Sds_Id = sfselect(sd_id,sfn2index(sd_id,"spherical_albedo"))
istatus = sfrdata(Sds_Id,Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%sph_alb) +  &
                  istatus
istatus = sfendacc(Sds_Id) + istatus

!--- 3d arrays
Sds_Start_3d = 0
Sds_Stride_3d = 1

!--- albedo
sds_dims_3d(1) = nlog10re
sds_dims_3d(2) = nlog10tau
sds_dims_3d(3) = nSolzen
Sds_Edge_3d(1) = nlog10re
Sds_Edge_3d(2) = nlog10tau
    Sds_Edge_3d(3) = nSolzen
    Sds_Id = sfselect(sd_id,sfn2index(sd_id,"albedo"))
    istatus = sfrdata(Sds_Id,Sds_Start_3d, Sds_Stride_3d, Sds_Edge_3d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%alb) +  &
                  istatus
    istatus = sfendacc(Sds_Id) + istatus

    !-- transmission
    Sds_Id = sfselect(sd_id,sfn2index(sd_id,"transmission"))
    istatus = sfrdata(Sds_Id,Sds_Start_3d, Sds_Stride_3d, Sds_Edge_3d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn) +  &
                  istatus
    istatus = sfendacc(Sds_Id) + istatus

    !--- 5d arrays
    Sds_Start_5d = 0
    Sds_Stride_5d = 1

    !--- reflectance
    sds_dims_5d(1) = nlog10re
    sds_dims_5d(2) = nlog10tau
    sds_dims_5d(3) = nSolzen
    sds_dims_5d(4) = nzen
    sds_dims_5d(5) = nRelaz

    Sds_Edge_5d(1) = nlog10re
    Sds_Edge_5d(2) = nlog10tau
    Sds_Edge_5d(3) = nSolzen
    Sds_Edge_5d(4) = nzen
    Sds_Edge_5d(5) = nRelaz

    Sds_Id = sfselect(sd_id,sfn2index(sd_id,"reflectance"))
    istatus = sfrdata(Sds_Id,Sds_Start_5d, Sds_Stride_5d, Sds_Edge_5d, &
                  Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ref) +  &
                  istatus
    istatus = sfendacc(Sds_Id) + istatus

    if (istatus /= 0) then
      print "(a,'Error reading sds data from cloud lut files, stopping')",EXE_PROMPT
      STOP 99
    end if

    !--- Close the lookup table file
    istatus = sfend(sd_id) 

    if (istatus /= 0) then
      print "(a,'Error closing cloud lut files')",EXE_PROMPT
    end if

   !--- assuming everything is equally spaced, compute spacing here
   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_Solzen = & 
           Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(2) - & 
           Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen(1)

   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_zen = & 
           Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(2) - & 
           Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(1)

   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_log10tau = & 
           Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(2) - & 
           Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(1)

   Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_log10re = & 
           Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(2) - & 
           Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(1)
 
print *, "Lut successfully read in ", trim(Lut_File), Chan_Idx, Phase_Idx
!=======================================================================
! end of this routine
!=======================================================================
 end subroutine READ_CLD_REF_LUT_STRUC


!=======================================================================
! Read Cloud Emissivity Lookup Table based on structures
! 
!  
! isat  = wmo satellite identIFication number
! Chan_Idx = ABI channel number
!
! sc_ind = spacecraft index
! scinfo(sc_ind)%id = wmo number
!
!=======================================================================
 subroutine READ_CLD_EMS_LUT_STRUC(Ancil_Data_Path,Algo_Idx,Chan_Idx,Phase_Idx)

  integer, intent(in):: Algo_Idx,Chan_Idx
  integer(kind=int1), intent(in):: Phase_Idx
  character (len=*), intent(in):: Ancil_Data_Path

  character (len=128):: Lut_File
  integer:: astatus,tstatus
  character (len=2):: Chan_String
  character (len=3):: Phase_String
  character (len=7):: Sat_String
  integer:: nlog10re
  integer:: nlog10tau
  integer:: nzen

  !-- hdf
  integer:: sfstart, sfrcatt,sfn2index, sffattr, sfendacc,sfend,sfrdata,sfselect
  integer:: sd_id, Sds_Id, istatus

  integer, parameter:: sds_rank_1d = 1
  integer, dimension(sds_rank_1d):: Sds_Start_1d, Sds_Edge_1d, Sds_Stride_1d

  integer, parameter:: sds_rank_3d = 3
  integer, dimension(sds_rank_3d):: Sds_Start_3d, Sds_Edge_3d, Sds_Stride_3d, sds_dims_3d

!--- construct appropriate file name - use actual satellite channel, not abi equivalent
! write(Chan_String, '(I1.1)') scinfo(sc_ind)%ch_flg(Chan_Idx)

  Chan_String = "0"
  if (index(Sensor%Sensor_Name,'AVHRR') > 0) then
    if (Chan_Idx == 20) then
          Chan_String = "3b"
    endif
  endif
  if (index(Sensor%Sensor_Name,'MODIS') > 0) then
    if (Chan_Idx == 20) then
          Chan_String = "20"
    endif
  endif
  if (index(Sensor%Sensor_Name,'GOES') > 0) then
    if (Chan_Idx == 20) then
          Chan_String = "2"
    endif
  endif

  Phase_String = "wat"
  if (Phase_Idx == ICE_PHASE) then
    Phase_String = "ice"
  end if
  if (Phase_Idx == MIXED_PHASE) then
    Phase_String = "wat"                  !note assumed water tables apply to mixed phase clouds
  end if

  !--- make filenames
  select case (Sensor%WMO_Id)

      case(3)
        Sat_String  = "METOP-A"
      case(4)
        Sat_String  = "METOP-B"
      case(5)
        Sat_String  = "METOP-C"
      case(252)
        Sat_String  = "GOES-08"
      case(253)
        Sat_String  = "GOES-09"
      case(254)
        Sat_String  = "GOES-10"
      case(255)
        Sat_String  = "GOES-11"
      case(256)
        Sat_String  = "GOES-12"
      case(257)
        Sat_String  = "GOES-13"
      case(258)
        Sat_String  = "GOES-14"
      case(259)
        Sat_String  = "GOES-15"
      case(706)
        Sat_String  = "NOAA-06"
      case(223)
        Sat_String  = "NOAA-19"
      case(783)
        Sat_String = "MODIS"
      case(784)
        Sat_String = "MODIS"

  end select

  if (trim(Chan_String) == "0") return

  Lut_File = TRIM(Sat_String)//"_ch"//TRIM(Chan_String)//"_ems_lut_"//Phase_String//"_cld.hdf"

  !--- Open lookup table
  sd_id = sfstart(TRIM(Ancil_Data_Path)//TRIM(Lut_File),DFACC_READ)

  if (sd_id < 0) then
     print *, EXE_PROMPT, "Error open cloud ems table ", trim(Ancil_Data_Path)//trim(Lut_File)
     stop 99
  endif

   !--- global attributes 
   istatus = 0

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"HEADER"),  &
                     Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%header)

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_EFF_RADII"),  &
                     nlog10re)

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_OPTICAL_DEPTHS"),  &
                     nlog10tau)

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"NUMBER_SENSOR_ZENITH_ANGLES"),  &
                     nzen) + istatus

   !--- store local dimensions into data structures
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10re = nlog10re
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nlog10tau = nlog10tau
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%nzen = nzen

   tstatus = 0
   allocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(nzen), stat=astatus)
   tstatus = tstatus + astatus
   allocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(nlog10re), stat=astatus)
   tstatus = tstatus + astatus
   allocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(nlog10tau), stat=astatus)
   tstatus = tstatus + astatus
   allocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio(nlog10re), stat=astatus)
   tstatus = tstatus + astatus
   allocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref(nlog10re), stat=astatus)
   tstatus = tstatus + astatus
   allocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%emiss(nlog10re,nlog10tau,nzen), stat=astatus)
   tstatus = tstatus + astatus
   allocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trans(nlog10re,nlog10tau,nzen), stat=astatus)
   tstatus = tstatus + astatus

   if (tstatus /= 0) then
     print "(a,'Error allocating cloud emissivity lookup table arrays')",EXE_PROMPT
      STOP 99
   end if

   !--- initialize
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen  = 0.0
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re  = 0.0
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau  = 0.0
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio  = 0.0
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref  = 0.0
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%emiss = 0.0
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trans = 0.0
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%rho = 0.0
   Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%lambda_ref = 0.0

  !--- Read in values in the table

  !-- read in vectors
  Sds_Start_1d = 0
  Sds_Stride_1d = 1
  Sds_Edge_1d = nzen

  !-- sensor zenith angle
  Sds_Edge_1d = nzen
  Sds_Id = sfselect(sd_id,sfn2index(sd_id,"sensor_zenith_angle"))
  istatus = sfrdata(Sds_Id, Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen) +  &
                  istatus
  istatus = sfendacc(Sds_Id) + istatus

  !-- log10 optical depth
  Sds_Edge_1d = nlog10tau
  Sds_Id = sfselect(sd_id,sfn2index(sd_id,"log10_optical_depth"))
  istatus = sfrdata(Sds_Id, Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau) +  &
                  istatus
  istatus = sfendacc(Sds_Id) + istatus

  !-- log10 effective radius
  Sds_Edge_1d = nlog10re
  Sds_Id = sfselect(sd_id,sfn2index(sd_id,"log10_eff_radius"))
  istatus = sfrdata(Sds_Id, Sds_Start_1d, Sds_Stride_1d, Sds_Edge_1d, &
                  Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re) +  &
                  istatus
  istatus = sfendacc(Sds_Id) + istatus

  !--- 3d arrays
  Sds_Start_3d = 0
  Sds_Stride_3d = 1

  !--- emissivity
  sds_dims_3d(1) = nlog10re
  sds_dims_3d(2) = nlog10tau
  sds_dims_3d(3) = nzen
  Sds_Edge_3d(1) = nlog10re
  Sds_Edge_3d(2) = nlog10tau
  Sds_Edge_3d(3) = nzen

  Sds_Id = sfselect(sd_id,sfn2index(sd_id,"cloud_emissivity"))
  istatus = sfrdata(Sds_Id, Sds_Start_3d, Sds_Stride_3d, Sds_Edge_3d, &
                  Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%emiss) +  &
                  istatus
  istatus = sfendacc(Sds_Id) + istatus

  !--- transmission
  Sds_Id = sfselect(sd_id,sfn2index(sd_id,"cloud_transmission"))

  istatus = sfrdata(Sds_Id, Sds_Start_3d, Sds_Stride_3d, Sds_Edge_3d, &
                    Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trans) +  &
                    istatus

  istatus = sfendacc(Sds_Id) + istatus


  !--- Close the lookup table file
  istatus = sfend(sd_id) + istatus

  !--- check status of read
  if (istatus /= 0) then
     print *, "Error reading emissivity table for channel, phase = ",  &
              Chan_Idx, Phase_Idx
     STOP 99
  endif

  !--- assuming everything is equally spaced, compute spacing here
  Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_zen = &
           Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(2) - &
           Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen(1)

  Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_log10tau = &
           Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(2) - &
           Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau(1)

  Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%delta_log10re = &
           Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(2) - &
           Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re(1)


print *, "Lut successfully read in ", trim(Lut_File)
!=======================================================================
! end of this routine
!=======================================================================
 end subroutine READ_CLD_EMS_LUT_STRUC

!=======================================================================
!
!=======================================================================
 subroutine ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE(Nalgo)
   integer, intent(in):: Nalgo
   integer:: astatus

     !--- allocate main structures
     allocate(Cld_Ref(Nalgo),stat=astatus)
     if (astatus /= 0) then
       print "(a,'Not enough memory to allocate cld reflectance lut structure.')",EXE_PROMPT
       STOP
     end if
     Cld_Ref(:)%flag = sym%NO

     allocate(Cld_Ems(Nalgo),stat=astatus)
     if (astatus /= 0) then
       print "(a,'Not enough memory to allocate cld emissivity lut structure.')",EXE_PROMPT
       STOP
     end if

     Cld_Ems(:)%flag = sym%NO

 end subroutine ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE

 !======================================================================
 ! Deallocate threshold arrays on exit of last segment
 !======================================================================
 subroutine DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES(Algo_Idx)

    integer, intent(in):: Algo_Idx
    integer:: dstatus        !deallocate status 
    integer:: tstatus        !accumulated deallocate status
    integer:: Chan_Idx          !channel index
    integer:: Phase_Idx         !phase index

    print *, EXE_PROMPT, " Destroying Cloud LUTs" 

    !--------------------------------------------------------------
    !--- reflectance tables
    !--------------------------------------------------------------
    dstatus = 0
    tstatus = 0
     do Chan_Idx = 1, Nchan_Max
      do Phase_Idx = 1,Nphases
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Solzen,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Relaz,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%sph_alb)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%sph_alb,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trn,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%alb)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%alb,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ref)) then
         deallocate(Cld_Ref(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ref,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
      end do
     end do

     !--- check status of deallocating
     if (tstatus /= 0) then
      print "(a,'Error deallocating cloud reflectance lookup table arrays')",EXE_PROMPT
      STOP 99
     endif

    !--------------------------------------------------------------
    !--- emissivity tables
    !--------------------------------------------------------------
    dstatus = 0
    tstatus = 0
     do Chan_Idx = 1, Nchan_Max
      do Phase_Idx = 1,Nphases

       if (allocated(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen)) then
         deallocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%zen,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau)) then
         deallocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10tau,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re)) then
         deallocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%log10re,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio)) then
         deallocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%ext_ratio,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref)) then
         deallocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%Qext_ref,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trans)) then
         deallocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%trans,stat=dstatus)
         tstatus = tstatus + dstatus
       endif
       if (allocated(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%emiss)) then
         deallocate(Cld_Ems(Algo_Idx)%channel(Chan_Idx)%phase(Phase_Idx)%emiss,stat=dstatus)
         tstatus = tstatus + dstatus
       endif

      enddo

     enddo

     !--- check status of deallocating
     if (tstatus /= 0) then
      print "(a,'Error deallocating cloud reflectance lookup table arrays')",EXE_PROMPT
      STOP 99
     endif

   end subroutine DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES

!----------------------------------------------------------------------
!--- routine to deallocate main cloud lut structure
!----------------------------------------------------------------------
   subroutine DEALLOCATE_MEMORY_FOR_MAIN_CLD_LUT_STRUCTURES()
    integer:: dstatus      !deallocate status
    integer:: tstatus      !accumulated status

    !--- initialize status flags
    dstatus = 0
    tstatus = 0

    !--- reflectance
    deallocate(Cld_Ref,stat=dstatus)
    tstatus = tstatus + dstatus
    if (tstatus /= 0) then
      print "(a,'Error deallocating cloud reflectance lookup table arrays')",EXE_PROMPT
      STOP 99
    endif

    !--- emissivity
    deallocate(Cld_Ems,stat=dstatus)
    tstatus = tstatus + dstatus
    if (tstatus /= 0) then
      print "(a,'Error deallocating cloud emissivity lookup table arrays')",EXE_PROMPT
      STOP 99
    endif

 end subroutine DEALLOCATE_MEMORY_FOR_MAIN_CLD_LUT_STRUCTURES

!====================================================================
! SUBROUTINE Name: KNOWING_Z_COMPUTE_T_P_RTM
!
! Function:
!  derive pressure and temperature from a profile knowing height
!
!====================================================================
 subroutine KNOWING_Z_COMPUTE_T_P_RTM(Lon_Idx,Lat_Idx,P,T,Z,ilev)

  !--- arguments
  integer, intent(in):: Lon_Idx
  integer, intent(in):: Lat_Idx
  real, intent(in):: Z
  real, intent(out):: T
  real, intent(out):: P
  integer, intent(out):: ilev

  !--- local variables
  real:: dp
  real:: dt
  real:: dz
  integer:: kstart
  integer:: kend
  integer:: NLevels_temp

  !--- compute region of atmosphere to look
  kstart = Rtm(Lon_Idx,Lat_Idx)%Tropo_Level
  kend = Rtm(Lon_Idx,Lat_Idx)%Sfc_Level
  NLevels_temp = kend - kstart + 1

  !--- interpolate pressure profile
  call LOCATE(Rtm(Lon_Idx,Lat_Idx)%Z_Prof(kstart:kend),NLevels_temp,Z,ilev)
  ilev = max(1,min(NLevels_temp-1,ilev))
  ilev = ilev + kstart - 1

  dp = P_Std_Rtm(ilev+1) - P_Std_Rtm(ilev)
  dt = Rtm(Lon_Idx,Lat_Idx)%T_Prof(ilev+1) - Rtm(Lon_Idx,Lat_Idx)%T_Prof(ilev)
  dz = Rtm(Lon_Idx,Lat_Idx)%Z_Prof(ilev+1) - Rtm(Lon_Idx,Lat_Idx)%Z_Prof(ilev)

  !--- perform interpolation
  if (dz /= 0.0) then
   T = Rtm(Lon_Idx,Lat_Idx)%T_Prof(ilev) + dt/dz * (Z - Rtm(Lon_Idx,Lat_Idx)%Z_Prof(ilev))
   P = P_Std_Rtm(ilev) + dp/dz * (Z - Rtm(Lon_Idx,Lat_Idx)%Z_Prof(ilev))
  else
   T = Rtm(Lon_Idx,Lat_Idx)%T_Prof(ilev)
   P = P_Std_Rtm(ilev)
  endif

 end subroutine KNOWING_Z_COMPUTE_T_P_RTM

end module CLOUD_TAU_RE_SOLAR_module
