! $Id: cloud_tau_re_solar_direct.f90,v 1.11 2013/06/08 18:47:18 wstraka Exp $

Module CLOUD_TAU_RE_SOLAR_MODULE
! Baseline Daytime Cloud optical and microphysical algorithm(s)
!
! author: Andrew Heidinger and the GOES-R AWG Cloud Team
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
!     LOCATE_TABLE_POSITION_ANGLES
!     CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER_NOINTERP - forward model for reflectance computation
!     CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER_NOINTERP - forward model for thermal computation
!     COMPUTE_SOLAR_TRANS_TOA
!     ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE
!     DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES
!
! version history:
!  04/28/07 - Created based on CLAVR-x algorithms
!
! All channels referenced to ABI (2=0.6 um, 5 = 1.6 um, 
!                                 6 = 2.3, 7 = 3.9 um)
!
! Notes:
! 1.  channel 6 is not included yet because not available on SEVIRI.
! 2. Error handling needs to be included.
! 3. GEOCAT hdf routines should be used.
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
            LOCATE_TABLE_POSITION_ANGLES, &
            CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER_NOINTERP, &  !forward model for reflectance computation
            CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER_NOINTERP, &   !forward model for thermal computation
            COMPUTE_SOLAR_TRANS_TOA, &
            ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE, &
            DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES

!-----------------------------------------------------------------------
! define module-wide variables
!-----------------------------------------------------------------------
  INTEGER(KIND=int1),parameter,private::WATER_PHASE = 1  !water index
  INTEGER(KIND=int1),parameter,private::MIXED_PHASE = 2 !mixed index
  INTEGER(KIND=int1),parameter,private::ICE_PHASE = 3   !ice index
  INTEGER(KIND=int1),parameter,private::nchan_max = 7   !max channels 

  REAL(KIND=REAL4), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE:: f
  REAL(KIND=REAL4), DIMENSION(:,:), ALLOCATABLE, PRIVATE:: cost


!--- include file
 INCLUDE 'baseline_cloud_micro_day.inc'

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
! imode = 1 = ch2 and ch5
! imode = 2 = ch2 and ch6  (not coded)
! imode = 3 = ch2 and ch7 total
! imode = 4 = ch2 and ch7 solar only
!=======================================================================
 SUBROUTINE CLOUD_TAU_RE_SOLAR(jmin, jmax, iseg, nseg)

!--- AVHRR specific additions
   INTEGER, INTENT(in):: jmin
   INTEGER, INTENT(in):: jmax
   INTEGER, INTENT(in):: iseg
   INTEGER, INTENT(in):: nseg
   INTEGER:: ichan
   INTEGER, parameter:: ialgo = 1
   INTEGER, parameter:: nalgo = 1
   CHARACTER(LEN=100):: lut_path
   INTEGER, SAVE:: sc_id_prev = 0
!--- end AVHRR specific additions
   INTEGER:: nlog10tau_max
   INTEGER:: nlog10re_max
   INTEGER:: nlog10tau
   INTEGER:: nlog10re
   INTEGER:: itau, itau_1, itau_2, ire, ire_1, ire_2, i, j, quadrant
   REAL (KIND=REAL4):: dtau, dre, dy1_dtau, dy2_dtau, dy1_dre, dy2_dre, &
                       log10tau_interp, log10re_interp, weight,weight_sum

   INTEGER:: itropo
   INTEGER:: xnwp
   INTEGER:: ynwp
   INTEGER:: ielem
   INTEGER:: iline
   INTEGER:: ivza
   INTEGER:: isfc
   INTEGER:: imode
   INTEGER:: ierror

   INTEGER(KIND=int4):: cloudphase_previous
   REAL(KIND=real4):: zen_previous
   REAL(KIND=real4):: sol_zen_previous
   REAL(KIND=real4):: az_previous
   REAL(KIND=real4):: reuse_angle_thresh
   LOGICAL:: reuse_forward_model

!--- local variables that are aliases to elements of global structures
   REAL:: ref2
   REAL:: ref5
   REAL:: ref7
   REAL:: rad7
   REAL:: bt14
   REAL(KIND=real4):: zen
   REAL(KIND=real4):: sol_zen
   REAL(KIND=real4):: cos_zen
   REAL(KIND=real4):: cos_sol_zen
   REAL(KIND=real4):: az
   REAL(KIND=real4):: Ref_std    !reflectance variability
   REAL(KIND=real4):: Ref_mean  !reflectance variability
   REAL(KIND=real4):: Ref_var  !reflectance variability
   INTEGER(KIND=int4):: cloudphase
   INTEGER(KIND=int4):: cloudtype
   INTEGER(KIND=int4):: cloudmask
   INTEGER(KIND=int4):: sfctype
   INTEGER(KIND=int4):: iparam
   REAL(KIND=real4):: tpw
   REAL(KIND=real4):: Tsfc
   REAL(KIND=real4):: rad7_clr
   REAL(KIND=real4):: T_cloud
   REAL(KIND=real4):: Z_cloud
   REAL(KIND=real4):: P_cloud
   REAL(KIND=real4):: cost_minimum
   REAL(KIND=real4):: cost_corner_minimum
   INTEGER(KIND=INT4), dimension(2):: cost_minimum_indices

!--- local pointers to geocat public data structures
   REAL(KIND=real4), DIMENSION(:), pointer:: tlev
   REAL(KIND=real4), DIMENSION(:), pointer:: plev
   REAL(KIND=real4), DIMENSION(:), pointer:: zlev
   REAL(KIND=real4), DIMENSION(:), pointer:: tpwlev

!--- 1d-var retrieval arrays
  INTEGER, parameter:: num_obs = 2
  REAL, DIMENSION(num_obs):: y

!--- local variables
  REAL(KIND=real4):: log10tau   !log10 of the optical depth (tau)
  REAL(KIND=real4):: log10re    !log10 of the effective radius (re)
  REAL(KIND=real4):: r          !log10re interpolation weight
  REAL(KIND=real4):: s          !log10tau interpolation weight (0-1)
  REAL(KIND=real4):: t          !interpolation weight (0-1)
  REAL(KIND=real4):: t2         !interpolation weight (0-1)
  REAL(KIND=real4):: u          !interpolation weight (0-1)
  REAL(KIND=real4):: v          !interpolation weight (0-1)
  INTEGER(KIND=int4):: ir       !interpolation index
  INTEGER(KIND=int4):: is       !interpolation index
  INTEGER(KIND=int4):: it       !interpolation index
  INTEGER(KIND=int4):: it2      !interpolation index
  INTEGER(KIND=int4):: iu       !interpolation index
  INTEGER(KIND=int4):: iv       !interpolation index
  REAL(KIND=real4):: fm_err     !forward model error
  INTEGER(KIND=int4):: isingular  !flag noting a singular matrix
  INTEGER(KIND=int4):: ifail      !flag noting a failed retrieval
  INTEGER(KIND=int4):: iter       !iteration loop index
  INTEGER(KIND=int4):: ilev       !level index
! REAL(KIND=real4):: pc           !cloud top pressure
! REAL(KIND=real4):: zc           !cloud top height
! REAL(KIND=real4):: tc           !cloud top temperature
  REAL(KIND=real4):: prof_weight  !profile interpolation weight
  REAL(KIND=real4):: tpw_ac       !above-cloud tpw_ac
  INTEGER(KIND=int1):: diag_output
  REAL(KIND=real4):: tau_gas      !gas optical depth
  REAL(KIND=real4):: tau_aer      !aerosol optical depth
  REAL(KIND=real4):: wo_aer       !aerosol single scatter albedo
  REAL(KIND=real4):: g_aer        !asymmetry parameter
  REAL(KIND=real4):: tau_ray      !rayleigh optical depth
  REAL(KIND=real4):: Ref_ss_2     !single scatter reflectance
  REAL(KIND=real4):: Ref_ss_5     !single scatter reflectance
  REAL(KIND=real4):: Ref_ss_7     !single scatter reflectance
  REAL(KIND=real4):: Ref_2_clear  !clear ch2 reflectance
  REAL(KIND=real4):: log10_tau_apriori  !apriori estimate

!--- spatial uniformity
  REAL(KIND=real4), DIMENSION(:,:), allocatable:: ref2_mean
  REAL(KIND=real4), DIMENSION(:,:), allocatable:: ref2_max
  REAL(KIND=real4), DIMENSION(:,:), allocatable:: ref2_min
  REAL(KIND=real4), DIMENSION(:,:), allocatable:: ref2_uni

!--- local variables
  REAL(KIND=real4):: trans_ac_2
  REAL(KIND=real4):: trans_ac_5
  REAL(KIND=real4):: trans_ac_7
  REAL(KIND=real4):: trans_ac_7_zen
  REAL(KIND=real4):: Rad_ac_7_zen
  REAL(KIND=real4):: trans_sfc_2
  REAL(KIND=real4):: trans_sfc_5
  REAL(KIND=real4):: trans_sfc_7
  REAL(KIND=real4):: alb_sfc_2
  REAL(KIND=real4):: alb_sfc_5
  REAL(KIND=real4):: alb_sfc_7
  REAL(KIND=real4):: emiss_sfc_7
  REAL(KIND=real4):: solar_7
  REAL(KIND=real4):: sed
  REAL(KIND=real4):: ch7_Rad_to_Ref_fac
  REAL(KIND=real4):: Ref_sol
  REAL(KIND=real4):: Rad_therm
  REAL(KIND=real4):: dB_dT
  REAL(KIND=real4):: bt
 
  !--- local variable used in lwp, iwp, N and H computation
  REAL(KIND=real4), PARAMETER:: rho_water = 1000.0    !kg/m^3
  REAL(KIND=real4), PARAMETER:: rho_ice = 915.0       !kg/m^3
  REAL(KIND=real4), PARAMETER:: Q_eff_sca  = 2.0
  REAL(KIND=real4), PARAMETER:: drop_dis_width  = 0.8
  REAL(KIND=real4):: condensation_rate
  REAL(KIND=real4):: lwp_temp
   

!-------------------------------------------------------------------------------
! Executable Code
!-------------------------------------------------------------------------------

!--- set diagnostic output flag
   diag_output = sym%NO

!--- initialize 1d-var arrays
   y = 0.0

!--- initialize to missing
 tau_vis = missing_value_real4
 re_vis = missing_value_real4
 lwp_vis = missing_value_real4
 iwp_vis = missing_value_real4
 iwp_tau_vis = missing_value_real4
 H_vis = missing_value_real4
 N_vis = missing_value_real4
 tau_vis_qf = missing_value_int1
 re_vis_qf = missing_value_int1

!--------------------------------------------------------------------------
! AVHRR: determine which mode is on
!--------------------------------------------------------------------------
 IF (maxval(ch3a_on) > 0) THEN
    imode = 1    !use 1.6 micron
 ELSE
    imode = 3    !use 3.75 micron total reflectance
!   imode = 4    !use 3.75 micron pseudo reflectance
 ENDIF

!----------------------------------------------------------------------
! store path to ancillary data in a local variable
!----------------------------------------------------------------------
  lut_path = TRIM(ancil_data_dir)//"/luts/cld/"

!----------------------------------------------------------------------
! check to see if main structures need to be allocated
!----------------------------------------------------------------------
  IF (.NOT. ALLOCATED(cld_ref)) THEN
      call ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE(nalgo)
  ENDIF

!----------------------------------------------------------------------
!--- on the first segment, read tables
!----------------------------------------------------------------------
  IF (iseg == 1) THEN


!--- check if proper tables are already read in
  IF ((sc_id_prev == 0) .OR. &
      (sc_id_prev /= Sc_Id_Internal)) THEN

!--- check if tables need to be deallocated first
       IF (sc_id_prev > 0) THEN
        call DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES(ialgo)
       ENDIF

!--- Read channel 2 tables
       CALL READ_TABLES(ialgo,2,lut_path)

!--- Read second channel depending on mode selected
       IF (imode == 1) THEN
           CALL READ_TABLES(ialgo,5,lut_path)
       END IF
       IF (imode == 2) THEN
           CALL READ_TABLES(ialgo,6,lut_path)
       END IF
       IF (imode == 3) THEN
           CALL READ_TABLES(ialgo,7,lut_path)
       END IF
       IF (imode == 4) THEN
           CALL READ_TABLES(ialgo,7,lut_path)
       END IF

   ENDIF
  END IF

!----------------------------------------------------------------------
! allocate memory for forward and cost results - use max size
!----------------------------------------------------------------------
ichan = 2
nlog10tau_max =  maxval(cld_ref(ialgo)%channel(ichan)%phase(:)%nlog10tau)
nlog10re_max =  maxval(cld_ref(ialgo)%channel(ichan)%phase(:)%nlog10re)
if (.NOT. ALLOCATED(f)) ALLOCATE(f(num_obs,nlog10tau_max, nlog10re_max))
if (.NOT. ALLOCATED(cost)) ALLOCATE(cost(nlog10tau_max, nlog10re_max))

!----------------------------------------------------------------------
! check to see if there is any daytime data in this segment
!----------------------------------------------------------------------
  sol_zen = minval(Solzen)
  IF ((sol_zen >= Solzen_max).and.(sol_zen /= missing_value_real4)) THEN
    return
  ENDIF

  sol_zen_previous = missing_value_real4
  zen_previous = missing_value_real4
  az_previous = missing_value_real4
  cloudphase_previous = -1

!======================================================================
! Loop over pixels in this segment
!======================================================================
  line_loop_1: Do iline= jmin, jmax+jmin-1
    element_loop_1: Do ielem=1, num_pix

!--- check for space pixel
      IF (bad_pixel_mask(ielem,iline) == sym%YES) THEN
        cycle
      END IF

!     IF(ielem /= num_pix/2) then
!       cycle
!     END IF

!--- initilatize quality flag to be lowest quality
        tau_vis_qf(ielem,iline) = 0
        re_vis_qf(ielem,iline) = 0

!--- define aliases
      ref2 = alb1(ielem,iline)         !0.6 micron reflectance
      ref5 = alb3a(ielem,iline)         !1.6 micron reflectance
      ref7 = alb3b(ielem,iline)         !3.9 micron reflectance
      rad7 = rad3b(ielem,iline)         !3.9 micron reflectance
      bt14  = bt4(ielem,iline)        !11.0 micron bt
      xnwp = i_nwp(ielem,iline)        !nwp longitude cell
      ynwp = j_nwp(ielem,iline)        !nwp latitude cell
      itropo = tropo_level_nwp(xnwp,ynwp)         !nwp level associated with tropopause
      ivza = ivza_rtm(ielem,iline)            !viewing zenith angle bin
      isfc = sfc_level_nwp(xnwp,ynwp)             !first nwp level above surface
      zen = satzen(ielem,iline)           !viewing zenith angle 
      sol_zen = Solzen(ielem,iline)        !solar zenith angle 
      cos_zen = coszen(ielem,iline)    !cosine viewing zenith angle 
      cos_sol_zen = cosSolzen(ielem,iline) !cosine solar zenith angle 
      az = Relaz(ielem,iline)             !relative solar azimuth angle 
      T_cloud = Tc_acha(ielem,iline)         !cloud temperature
      Z_cloud = Zc_acha(ielem,iline)         !cloud temperature

!--- select cloud type to use
      if (cloud_type_ir_flag == sym%YES) then
       cloudtype = cld_type_ir(ielem,iline)
      else
       cloudtype = cld_type(ielem,iline)
      endif

      cloudmask = cld_mask(ielem,iline)    !cloud mask
      sfctype = sfc_type(ielem,iline)     !surface type
      Ref_std = alb1_std_3x3(ielem,iline)   !measure of reflectance variability
      Ref_mean = alb1_mean_3x3(ielem,iline)   !measure of reflectance variability
      Ref_var = 0.0
      if (Ref_mean > 0.0) then
          Ref_var = Ref_std / Ref_mean
      endif
      Ref_var = min(1.0,Ref_var) 
      plev => P_std_nwp
      tlev => T_prof_nwp(:,xnwp,ynwp)
      zlev => Z_prof_nwp(:,xnwp,ynwp)
      tpwlev => tpw_prof_nwp(:,xnwp,ynwp)
      tpw  = tpw_nwp(xnwp,ynwp)                       !total column pw
      Tsfc  = tsfc_rtm(ielem,iline)               !surface temperature
      emiss_sfc_7 = Sfc_Emiss_3b(ielem,iline)     !3.9 micron surface emissivity
      rad7_clr = Rad_Clear_Ch3b_rtm(ielem,iline)  !clear 4 micron radiance
      solar_7 = solar_3b_nu                       !solar energy in channel 7
      sed = sun_earth_distance                !sun earth distance factor

!--- compute factor for converting channel reflectance into a radiances
      ch7_Rad_to_Ref_fac = pi / cos_sol_zen /  (solar_7/sed**2) 

!--- check for correct sensor and solar geometery
      IF (zen >= zen_max) THEN
        cycle
      END IF
      IF (sol_zen >= Solzen_max) THEN
        cycle
      END IF

!--- check for clear pixels 
      IF (cloudmask == sym%CLEAR .or.  &
          cloudmask == sym%PROB_CLEAR .or.  &
          cloudmask == sym%PROB_CLOUDY .or.  &
          cloudtype == sym%CLEAR_TYPE) THEN
        cycle
      END IF

!--- check for valid cloud height
      IF (T_cloud == missing_value_real4) THEN
        cycle
      END IF

!--- check for correct phase
      IF (cloudtype == sym%UNKNOWN_TYPE) THEN
        cycle
      END IF

!--- determine phase for this algorithm 
      cloudphase = ICE_PHASE
      IF ((cloudtype == sym%FOG_TYPE) .or.  &
          (cloudtype == sym%WATER_TYPE) .or. & 
          (cloudtype == sym%SUPERCOOLED_TYPE)) THEN
          cloudphase = WATER_PHASE
      END IF
      IF ((cloudtype == sym%MIXED_TYPE)) then
          cloudphase = MIXED_PHASE
      ENDIF

!----------------------------------------------------------------------
!--- Define the Observations
!----------------------------------------------------------------------

!--- define y-vector - the observation vector
   IF (imode == 1) THEN  ! 0.6 and 1.6 micron approach
      y(1) = ref2/100.0
      y(2) = ref5/100.0
   END IF
   IF (imode == 3) THEN  ! 0.6 and 3.9 micron based approach)
      y(1) = ref2/100.0
      y(2) = rad7 * ch7_Rad_to_Ref_fac
   END IF
   IF (imode == 4) THEN  ! 0.6 and 3.9 micron pseudo reflectance based approach)
      y(1) = ref2/100.0
      y(2) = ref7/100.0
   END IF

!----------------------------------------------------------------------
! compute needed transmission terms
!----------------------------------------------------------------------

     !--- determine amount of water vapor above the cloud
     call KNOWING_Z_COMPUTE_T_P_NWP(xnwp,ynwp,P_cloud,T_cloud,Z_cloud,ilev)  

     prof_weight = (Z_cloud - zlev(ilev)) /  (zlev(ilev+1)-zlev(ilev))

     tpw_ac = tpwlev(ilev) + &
              prof_weight*(tpwlev(ilev+1) - tpwlev(ilev))

     !--- compute above cloud transmission for solar channels
     trans_ac_2 = 1.0
     trans_sfc_2 = 1.0
     trans_ac_2 =  COMPUTE_SOLAR_TRANS_TOA(tpw_ac,cos_zen,cos_sol_zen,h2o_tau_coef(Sc_Id_Internal,1,:))
     trans_sfc_2 =  COMPUTE_SOLAR_TRANS_TOA(tpw,cos_zen,cos_sol_zen,h2o_tau_coef(Sc_Id_Internal,1,:))

     trans_ac_5 = 1.0
     trans_sfc_5 = 1.0
     IF (imode == 1) THEN 
          trans_ac_5 =  COMPUTE_SOLAR_TRANS_TOA(tpw_ac,cos_zen,cos_sol_zen,h2o_tau_coef(Sc_Id_Internal,3,:))
          trans_sfc_5 =  COMPUTE_SOLAR_TRANS_TOA(tpw,cos_zen,cos_sol_zen,h2o_tau_coef(Sc_Id_Internal,3,:))
     END IF

     !--- use rtm structure for channel 7 - 3.9
     trans_ac_7 = 1.0
     trans_sfc_7 = 1.0
     IF (imode == 3 .or. imode == 4) THEN
          trans_ac_7_zen =  rtm(xnwp,ynwp)%d(ivza)%Trans_Atm_Clr3b(ilev) +    &
               prof_weight*(rtm(xnwp,ynwp)%d(ivza)%Trans_Atm_Clr3b(ilev+1) -  &
                            rtm(xnwp,ynwp)%d(ivza)%Trans_Atm_Clr3b(ilev))
          trans_ac_7 = trans_ac_7_zen**(1.0/cos_sol_zen)                    !now is the full path
          trans_sfc_7 = Trans_Atm_Ch3b_solar_rtm(ielem,iline)

          Rad_ac_7_zen =  rtm(xnwp,ynwp)%d(ivza)%Rad_Atm_Clr3b(ilev) +    &
               prof_weight*(rtm(xnwp,ynwp)%d(ivza)%Rad_Atm_Clr3b(ilev+1) -  &
                            rtm(xnwp,ynwp)%d(ivza)%Rad_Atm_Clr3b(ilev))
     END IF

!----------------------------------------------------------------------
! set surface reflectance (0-1)
!----------------------------------------------------------------------
alb_sfc_2 = ch1_sfc_alb_umd(sfc_type(ielem,iline))
alb_sfc_5 = ch3a_sfc_alb_umd(sfc_type(ielem,iline))
alb_sfc_7 = ch3b_sfc_alb_umd(sfc_type(ielem,iline))

!--- use value from clear composite if available over land
 IF (sfc_type(ielem,iline) /= sym%WATER_SFC) THEN
  IF (modis_clr_alb_flag == sym%YES) THEN
   IF (alb1_sfc_white_sky(ielem,iline) > 0.0)  THEN
       alb_sfc_2 = alb1_sfc_white_sky(ielem,iline)/100.0
   ENDIF
   IF (alb3a_sfc_white_sky(ielem,iline) > 0.0)  THEN
       alb_sfc_5 = alb3a_sfc_white_sky(ielem,iline)/100.0
   ENDIF
  ENDIF
 ENDIF

!--- use sfc emissivitys for surface reflectance for 3.75 micron
 IF ((emiss_sfc_7 > 0.0) .and. (sfc_type(ielem,iline) > 0))  THEN   
            alb_sfc_7 = 1.0 - emiss_sfc_7
 ENDIF  

!------------------------------------------------------------------------------
!   estimate clear-sky toa ch2 reflectance
!------------------------------------------------------------------------------

    !--- compute rayleigh and aerosol contribution to 0.63 micron signal
     tau_gas = -1.0*alog(trans_ac_2)
     if (tau_gas < 0.0) then 
        print *, "negative gas ", trans_ac_2, tau_gas
        stop
     endif

     tau_aer =  TAU_AEROSOL_BACKGROUND * (P_cloud / psfc_nwp(xnwp,ynwp))**2

     tau_ray =  tau_ray_corr(Sc_Id_Internal,1) * (P_cloud / psfc_nwp(xnwp,ynwp))

     CALL COMPUTE_CLEAR_SKY_SCATTER(tau_aer, &
                                    WO_AEROSOL_BACKGROUND, &
                                    G_AEROSOL_BACKGROUND, &
                                    tau_ray, &
                                    tau_gas, &
                                    scatangle(ielem,iline), &
                                    coszen(ielem,iline), &
                                    cosSolzen(ielem,iline),  &
                                    Ref_ss_2)

     Ref_ss_2 = Ref_ss_2 / 100.0   !convert from % to 0-1.0

     !--- ignore single scattering reflectance in longer wavelength channels
     Ref_ss_5 = 0.0
     Ref_ss_7 = 0.0

     !--- compute final estimate of clear reflectance (0-1.0)
     Ref_2_clear = Ref_ss_2 + trans_sfc_2*alb_sfc_2


!----------------------------------------------------------------------!
! call pixels darker than the clear reflectance missing
!----------------------------------------------------------------------
  IF (y(1) < Ref_2_clear) then
      cycle
  ENDIF



!---------------------------------------------------------------------
! check to see if you reuse the previous forward model computations
!---------------------------------------------------------------------

reuse_forward_model = .false.
reuse_angle_thresh = 5.0

if ( (cloudphase == cloudphase_previous) .and. &
     (abs(sol_zen-sol_zen_previous) < reuse_angle_thresh) .and. &
     (abs(zen-zen_previous) < reuse_angle_thresh) .and. &
     (abs(az-az_previous) < reuse_angle_thresh)) then
     reuse_forward_model = .true.
!    print *, "reusing forward model! ", sol_zen, sol_zen_previous
else
     reuse_forward_model = .false.
     sol_zen_previous = sol_zen
     zen_previous = zen
     az_previous = az
     cloudphase_previous = cloudphase
endif


IF (reuse_forward_model .eqv. .false.) THEN

!----------------------------------------------------------------------
! Loop over all tau and re and make a toa reflectance matrices
!----------------------------------------------------------------------
nlog10tau =  cld_ref(ialgo)%channel(2)%phase(cloudphase)%nlog10tau
nlog10re =  cld_ref(ialgo)%channel(2)%phase(cloudphase)%nlog10re

!--- find indices for this viewing geometry
 CALL LOCATE_TABLE_POSITION_ANGLES(ialgo,2,cloudphase,sol_zen,zen,az, &
                            it,it2,iu,iv)

!-------------------------------------------------------------------------------
do itau = 1, nlog10tau

 log10tau =  cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10tau(itau)

 do ire = 1, nlog10re

    log10re =  cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10tau(ire)

!--- ch2 forward model
    CALL CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER_NOINTERP(&
                         ialgo,2,cloudphase, &
                         ire,itau,it,it2,iu,iv, &
                         trans_ac_2,trans_sfc_2,alb_sfc_2,Ref_ss_2, &
                         f(1,itau,ire))

! print *, "vis ref = ", itau, log10tau, ire, log10re, f(1,itau,ire)
!--- ch5 forward model
   IF (imode == 1) THEN
    CALL CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER_NOINTERP(&
                         ialgo,5,cloudphase, &
                         ire,itau,it,it2,iu,iv, &
                         trans_ac_5,trans_sfc_5,alb_sfc_5, Ref_ss_5, &
                         f(2,itau,ire))
   ENDIF

!--- ch7 forward model
   IF (imode == 3 .or. imode == 4) THEN

!--- solar component
    CALL CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER_NOINTERP(&
                         ialgo,7,cloudphase, &
                         ire,itau,it,it2,iu,iv, &
                         trans_ac_7,trans_sfc_7,alb_sfc_7, Ref_ss_7, &
                         Ref_sol)

!--- thermal component
    Rad_therm = 0.0
    if (imode == 3) then
       CALL CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER_NOINTERP(&
                         ialgo,7,cloudphase, &
                         ire,itau,iu, &
                         rad7_clr, T_cloud, trans_ac_7_zen,Rad_ac_7_zen, &
                         Rad_therm,bt)
    ENDIF

!--- combine solar and thermal
    f(2,itau,ire) = Ref_sol + Rad_therm * ch7_Rad_to_Ref_fac

   ENDIF


!-----------------------------------------------------------------------
! End loop over all tau and re values
!-----------------------------------------------------------------------
  ENDDO 
 ENDDO

 ENDIF    !END CHECK FOR NEED TO RECOMPUTE TABLES

!-----------------------------------------------------------------------
! FIND OPTIMAL SOLUTION
!-----------------------------------------------------------------------
 !--- determine cost function
 cost(1:nlog10tau,1:nlog10re) =                        &
          (y(1) - f(1,1:nlog10tau,1:nlog10re))**2 +    &
          (y(2) - f(2,1:nlog10tau,1:nlog10re))**2

 !--- find location of cost minimum
 cost_minimum_indices = minloc(cost(1:nlog10tau,1:nlog10re))

!print *, "cost minimum indices = ", cost_minimum_indices
 !--- save values at mininum
 itau = cost_minimum_indices(1)
 ire = cost_minimum_indices(2)
 log10tau =  cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10tau(itau)
 log10re =  cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10re(ire)
 cost_minimum = cost(itau,ire)



!--- knowing optimal table indices, interpolate to find best 
itau_1 = max(1,itau-1)
itau_2 = min(nlog10tau,itau+1)
ire_1 = max(1,ire-1)
ire_2 = min(nlog10re,ire+1)

!-- find local corner
cost_corner_minimum = huge(1.00)
if (cost(itau_1,ire_1) <= cost_corner_minimum) then
  quadrant = 1
  cost_corner_minimum = cost(itau_1,ire_1)
endif
if (cost(itau_2,ire_1) <= cost_corner_minimum) then
  quadrant = 2
  cost_corner_minimum = cost(itau_2,ire_1)
endif
if (cost(itau_1,ire_2) <= cost_corner_minimum) then
  quadrant = 3
  cost_corner_minimum = cost(itau_1,ire_2)
endif
if (cost(itau_2,ire_2) <= cost_corner_minimum) then
  quadrant = 4
  cost_corner_minimum = cost(itau_2,ire_2)
endif

!--- define bounding indices of the quadrant
if (quadrant == 1) then !top left
  itau_1 = itau_1
  itau_2 = itau
  ire_1 = ire_1
  ire_2 = ire
endif
if (quadrant == 2) then !top right
  itau_1 = itau
  itau_2 = itau_2
  ire_1 = ire_1
  ire_2 = ire
endif
if (quadrant == 3) then   !bottom left
  itau_1 = itau_1
  itau_2 = itau
  ire_1 = ire
  ire_2 = ire_2
endif
if (quadrant == 4) then  !bottom right
  itau_1 = itau
  itau_2 = itau_2
  ire_1 = ire
  ire_2 = ire_2
endif
!--- find interpolated values
 weight_sum = 0.0
 log10tau_interp = 0.0
 log10re_interp = 0.0
 do i = itau_1,itau_2
   do j = ire_1,ire_2
      weight = (1.0/cost(i,j))
      weight_sum= weight_sum + weight
      log10tau_interp = log10tau_interp +  &
            cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10tau(i)*weight
      log10re_interp = log10re_interp + &
            cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10re(j)*weight
!           print *, i, j, cost(i,j),cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10tau(i), &
!                    cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10re(j), weight_sum
   enddo
 enddo
 log10tau_interp = log10tau_interp / weight_sum
 log10re_interp = log10re_interp / weight_sum

!print *, "optimal values = ", log10tau, log10re
!print *, "interp values = ", log10tau_interp, log10re_interp
!print *, "quadrant " ,quadrant

 !--- compute Kernel Matrix
 dtau = (cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10tau(itau_2) - &
         cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10tau(itau_1))

 dre = (cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10re(ire_2) - &
       cld_ref(ialgo)%channel(ichan)%phase(cloudphase)%log10re(ire_1))

 dy1_dtau = (f(1,itau_2,ire) - f(1,itau_1,ire))/dtau
 dy2_dtau = (f(2,itau_2,ire) - f(2,itau_1,ire))/dtau

 dy1_dre = (f(1,itau,ire_2) - f(1,itau,ire_1))/dre
 dy2_dre = (f(2,itau,ire_2) - f(2,itau,ire_1))/dre

!print *, "y = ", y
!print *, "optimal f = ", f(:,itau,ire)

 !--- convert back
 tau_vis(ielem,iline) = 10.0**log10tau_interp
 re_vis(ielem,iline) = 10.0**log10re_interp

!-----------------------------------------------------------------------
! END loop over pixels in segment
!-----------------------------------------------------------------------
    END DO element_loop_1
  END DO line_loop_1

!======================================================================
! Nullify local pointers
!======================================================================
  IF (iseg == nseg) THEN
!  tau => null()
!  re => null()
!  lwp => null()
!  iwp => null()
!  qf => null()
   plev => null()
   tlev => null()
   zlev => null()
   tpwlev => null()
  END IF

!======================================================================
! Update previous sc_id name for orbit
!======================================================================
  IF (iseg == nseg) THEN
   sc_id_prev = Sc_Id_Internal
  END IF


 END SUBROUTINE CLOUD_TAU_RE_SOLAR

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
 SUBROUTINE LOCATE_TABLE_POSITION_ANGLES(ialgo,ichan,iphase,Solzen,zen,Relaz, &
                                   it,it2,iu,iv)

    INTEGER(KIND=int4), INTENT(in):: ialgo,ichan
    INTEGER(KIND=int4), INTENT(in):: iphase
    REAL(KIND=real4), INTENT(in):: Solzen,zen,Relaz
    INTEGER, INTENT(out):: it,iu,iv,it2

!--- interpolate to find reflectance for this set of properties
   it = max(1,min(cld_ref(ialgo)%channel(ichan)%phase(iphase)%nSolzen-1,       &
         1+int((Solzen - cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen(1)) / &
               (cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_Solzen))))

   it2 = max(1,min(cld_ref(ialgo)%channel(ichan)%phase(iphase)%nSolzen-1,  &
          1+int((zen - cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen(1)) / &
                (cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_Solzen))))

   iu = max(1,min(cld_ref(ialgo)%channel(ichan)%phase(iphase)%nzen-1, &
          1+int((zen - cld_ref(ialgo)%channel(ichan)%phase(iphase)%zen(1)) / &
                (cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_zen))))

   !--- handle relative azimuth switching
   if (Relaz >= cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz_switch_pp) then

     iv = max(0,min(cld_ref(ialgo)%channel(ichan)%phase(iphase)%nRelaz_in_pp-2,         &
              int((Relaz - cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz_switch_pp) / &
               (cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_Relaz_in_pp))))  +        &
               cld_ref(ialgo)%channel(ichan)%phase(iphase)%nRelaz_out_pp + 1

   else

     iv = max(1,min(cld_ref(ialgo)%channel(ichan)%phase(iphase)%nRelaz_out_pp,         &
              1+int((Relaz - cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz(1)) / &
               (cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_Relaz_out_pp))))

   endif

 end SUBROUTINE LOCATE_TABLE_POSITION_ANGLES

!======================================================================
! Forward Model for a Reflectance Channel
!
! this routine is appropriate for single layer cloud placed in an 
! absorbing atmosphere over a Lambertian Surface.  The accounting
! for the surface contribution is given by Eq 11 in 
! King, Michael, 1987: JAS, vol44, 1734-1751 
!======================================================================
SUBROUTINE CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER_NOINTERP(&
                         ialgo,ichan,iphase, &
                         ir,is,it,it2,iu,iv, &
                         trans_ac,trans_sfc, alb_sfc, Ref_ss, &
                         ref)

  INTEGER, INTENT(in):: ialgo,ichan,       &
                        ir,is,it,it2,iu,iv
  INTEGER(KIND=int4), INTENT(in):: iphase
  REAL, INTENT(in):: trans_ac
  REAL, INTENT(in):: trans_sfc
  REAL, INTENT(in):: alb_sfc
  REAL, INTENT(in):: Ref_ss
  REAL, INTENT(out):: ref
  REAL:: trn
  REAL:: trn2
  REAL:: sab,temp,trans_bc,alb_sfc_temp

!----------------------------------------------------------------------
! construct sub-tables from appropriate lookup tables
!----------------------------------------------------------------------
  ref = cld_ref(ialgo)%channel(ichan)%phase(iphase)%ref(ir,is,it,iu,iv)
  trn = cld_ref(ialgo)%channel(ichan)%phase(iphase)%trn(ir,is,it)
  trn2 = cld_ref(ialgo)%channel(ichan)%phase(iphase)%trn(ir,is,it2)
  sab = cld_ref(ialgo)%channel(ichan)%phase(iphase)%sph_alb(ir,is)

!--- Estimate transmission from cloud to surface
 trans_bc = 0.0
 IF (trans_ac > 0.0) THEN
  trans_bc = trans_sfc / trans_ac
 END IF

!--- Account for surface reflection and atmospheric transmission
 alb_sfc_temp = trans_bc * alb_sfc
 temp = alb_sfc_temp / (1.0 - alb_sfc_temp * sab)
 ref = trans_ac * (ref + temp * trn*trn2) + Ref_ss

end SUBROUTINE CLOUD_FORWARD_REFLECTANCE_MODEL_SINGLE_LAYER_NOINTERP

!======================================================================
! Forward Model for a Thermal Channel
!
! this routine is appropriate for single layer cloud placed in an 
! absorbing atmosphere over a Lambertian Surface. 
!======================================================================
SUBROUTINE CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER_NOINTERP(&
                         ialgo,ichan,iphase, &
                         ir,is,iu, &
                         Rad_clr,Tc,trans_ac,Rad_ac,  &
                         rad,bt)

  INTEGER, INTENT(in):: ialgo
  INTEGER, INTENT(in):: ichan
  INTEGER, INTENT(in):: ir
  INTEGER, INTENT(in):: is
  INTEGER, INTENT(in):: iu
  INTEGER(KIND=int4), INTENT(in):: iphase
  INTEGER(KIND=int4):: ichan_avhrr
  REAL, INTENT(in):: Rad_clr
  REAL, INTENT(in):: Tc
  REAL, INTENT(in):: trans_ac
  REAL, INTENT(in):: Rad_ac
  REAL, INTENT(out):: rad
  REAL, INTENT(out):: bt


  REAL:: ems
  REAL:: trn
  REAL:: Bc

!---------------------------------------------------------------------
! translate channel to AVHRR
!---------------------------------------------------------------------
   ichan_avhrr = 0
   IF (ichan == 7) THEN 
       ichan_avhrr = 3
   ELSEIF (ichan == 15) THEN 
       ichan_avhrr = 4
   ELSEIF (ichan == 16) THEN 
       ichan_avhrr = 5
   ELSE
       PRINT *, EXE_PROMPT, "Invalid channel in cloud ems forward model"
   ENDIF

!----------------------------------------------------------------------
! construct sub-tables from appropriate lookup tables
!----------------------------------------------------------------------
  ems = cld_ems(ialgo)%channel(ichan)%phase(iphase)%emiss(ir,is,iu)
  trn = cld_ems(ialgo)%channel(ichan)%phase(iphase)%trans(ir,is,iu)

!--- compute planck emission at cloud temperature
 Bc  = planck_Rad_fast(ichan_avhrr,Tc)

!--- compute toa radiance and its Jacobian
 rad = ems*Rad_ac + trans_ac * ems * Bc +  trn * Rad_clr

 bt = planck_temp_fast(ichan_avhrr,rad)

end SUBROUTINE CLOUD_FORWARD_EMISSIVITY_MODEL_SINGLE_LAYER_NOINTERP

!=======================================================================
 Function COMPUTE_SOLAR_TRANS_TOA(tpw_ac,coszen,cosSolzen,coef) result(trans_solar_toa)

  REAL (KIND=real4), INTENT(in):: tpw_ac,coszen,cosSolzen
  REAL (KIND=real4), INTENT(in), DIMENSION(:):: coef
  REAL (KIND=real4):: trans_solar_toa
  REAL (KIND=real4):: tau
  REAL (KIND=real4):: airmass

  !--- compute airmass
  airmass = 1.0/coszen + 1.0/cosSolzen

  !--- compute optical depth
  tau = coef(1) + coef(2)*tpw_ac + coef(3)*(tpw_ac**2) 

  !--- constrain to be positive
  tau = max(0.0, tau)
 
  !--- compute transmission
  trans_solar_toa = exp ( -airmass * tau)

 END Function COMPUTE_SOLAR_TRANS_TOA


!=======================================================================
! Generic Routine to read the tables for a channel
!=======================================================================
SUBROUTINE READ_TABLES(ialgo,ichan,lut_path)

 INTEGER(KIND=int4), INTENT(in):: ialgo,ichan
 character(LEN=*), INTENT(in):: lut_path


!--- read in reflectance tables for all applicable channels
      IF (ichan <= 7) THEN

!--- set flag to indicate that this channel has been read in 
         cld_ref(ialgo)%channel(ichan)%flag = sym%YES

!--- set flag and read in water phase tables
          cld_ref(ialgo)%channel(ichan)%phase(WATER_PHASE)%flag = sym%YES
          CALL READ_CLD_REF_LUT_STRUC(lut_path,ialgo,ichan,WATER_PHASE)

!--- set flag and read in ice phase tables
          cld_ref(ialgo)%channel(ichan)%phase(ICE_PHASE)%flag = sym%YES
          CALL READ_CLD_REF_LUT_STRUC(lut_path,ialgo,ichan,ICE_PHASE)

!--- set flag and read in mixed phase tables
          cld_ref(ialgo)%channel(ichan)%phase(MIXED_PHASE)%flag = sym%YES
          CALL READ_CLD_REF_LUT_STRUC(lut_path,ialgo,ichan,MIXED_PHASE)

      END IF


!--- read in emissivity tables for all applicable channels
      IF (ichan >= 7) THEN

!--- set flag to indicate that this channel has been read in 
         cld_ems(ialgo)%channel(ichan)%flag = sym%YES

!--- set flag and read in water phase tables
          cld_ems(ialgo)%channel(ichan)%phase(WATER_PHASE)%flag = sym%YES
          CALL READ_CLD_EMS_LUT_STRUC(lut_path,ialgo,ichan,WATER_PHASE)

!--- set flag and read in ice phase tables
          cld_ems(ialgo)%channel(ichan)%phase(ICE_PHASE)%flag = sym%YES
          CALL READ_CLD_EMS_LUT_STRUC(lut_path,ialgo,ichan,ICE_PHASE)

!--- set flag and read in mixed phase tables
          cld_ems(ialgo)%channel(ichan)%phase(MIXED_PHASE)%flag = sym%YES
          CALL READ_CLD_EMS_LUT_STRUC(lut_path,ialgo,ichan,MIXED_PHASE)

      END IF


END SUBROUTINE READ_TABLES
!=======================================================================
! Read Water Cloud Lookup Table based on structures
! 
!
! isat  = wmo satellite identIFication number
! ichan = ABI channel number
!
! sc_ind = spacecraft index
! scinfo(sc_ind)%id = wmo number
!
! This does ice and water
!=======================================================================
 SUBROUTINE READ_CLD_REF_LUT_STRUC(ancil_data_path,ialgo,ichan,iphase)

  INTEGER, INTENT(in):: ialgo,ichan
  INTEGER(KIND=int1), INTENT(in):: iphase
  CHARACTER (LEN=*), INTENT(in):: ancil_data_path
 
  CHARACTER (LEN=128):: lut_file
  INTEGER:: astatus
  INTEGER:: tstatus
  CHARACTER (LEN=2):: chan_string
  CHARACTER (LEN=3):: phase_string
  CHARACTER (LEN=7):: sat_string
  INTEGER:: nlog10re
  INTEGER:: nlog10tau
  INTEGER:: nSolzen
  INTEGER:: nzen
  INTEGER:: nRelaz           !total number relative azimuths
  INTEGER:: nRelaz_in_pp     !number of Relaz in the principal plane zone
  INTEGER:: nRelaz_out_pp    !number of Relaz outside the principal plane zone
  REAL:: delta_Relaz_in_pp   !azimuth space in principal plane zone
  REAL:: delta_Relaz_out_pp  !azimuth space outside principal plane zone
  REAL:: Relaz_switch_pp     !azimuth defining principal plane zone

!-- hdf
INTEGER:: sfstart, sfrcatt,sfn2index, sffattr, sfendacc,sfend,sfrdata,sfselect
INTEGER:: sd_id
INTEGER:: sds_id
INTEGER:: istatus

INTEGER, PARAMETER:: sds_rank_1d = 1
INTEGER, DIMENSION(sds_rank_1d):: sds_start_1d, sds_edge_1d, sds_stride_1d

INTEGER, PARAMETER:: sds_rank_2d = 2
INTEGER, DIMENSION(sds_rank_2d):: sds_start_2d, sds_edge_2d, sds_stride_2d, sds_dims_2d

INTEGER, PARAMETER:: sds_rank_3d = 3
INTEGER, DIMENSION(sds_rank_3d):: sds_start_3d, sds_edge_3d, sds_stride_3d, sds_dims_3d

INTEGER, PARAMETER:: sds_rank_5d = 5
INTEGER, DIMENSION(sds_rank_5d):: sds_start_5d, sds_edge_5d, sds_stride_5d, sds_dims_5d

!----------------------------------------------------------------------------------
! Executable Code
!----------------------------------------------------------------------------------

!--- construct appropriate file name - use actual satellite channel, not abi equivalent
  phase_string = "wat"
  IF (iphase == ICE_PHASE) THEN
    phase_string = "ice"
  END IF
  IF (iphase == MIXED_PHASE) THEN
    phase_string = "wat"                  !note assumed water tables apply to mixed phase clouds
  END IF

  IF (ichan == 2) THEN
          chan_string = "1"
  ELSEIF (ichan == 5) THEN
          chan_string = "3a"
  ELSEIF (ichan == 7) THEN
          chan_string = "3b"
  ENDIF

!--- make filenames
  IF (phase_string /= "ice") THEN
    IF (Sc_Id_Internal < 4) THEN
        IF (Sc_Id_Internal == 2) THEN
          sat_string = "METOP-A"
        ENDIF
        IF (Sc_Id_Internal == 3) THEN
          sat_string = "METOP-B"
        ENDIF
    ENDIF
    IF (Sc_Id_Internal < 10 .and. Sc_Id_Internal > 4) THEN
         write(unit=sat_string,fmt="(a6,i1)") "NOAA-0",Sc_Id_Internal
    ENDIF
    IF (Sc_Id_Internal >= 10) THEN
         write(unit=sat_string,fmt="(a5,i2)")  "NOAA-",Sc_Id_Internal
    ENDIF

    lut_file = TRIM(sat_string)//"_ch"//TRIM(chan_string)//"_Ref_lut_"//phase_string//"_cld.hdf"

  ELSE

    sat_string = "AVHRR"
    lut_file = TRIM(sat_string)//"_ch"//TRIM(chan_string)//"_Ref_lut_mix1_"//phase_string//"_cld.hdf"

  ENDIF

!--- Open lookup table

!-- hdf
 sd_id = sfstart( &
                 TRIM(ancil_data_path)//TRIM(lut_file), &
                 DFACC_READ) 

 IF (sd_id < 0) THEN
     PRINT *, EXE_PROMPT, "Error open cloud ref table ", trim(ancil_data_path)//trim(lut_file)
     stop
  ENDIF

!---- Read in tables

!--- global attributes 
   istatus = 0

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"HEADER"),  &
                     cld_ref(ialgo)%channel(ichan)%phase(iphase)%header)

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

   IF (istatus /= 0) THEN
      PRINT "(a,'Error reading file attributes')",EXE_PROMPT
      stop
   ENDIF

!--- store local dimensions into data structures
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%nlog10re = nlog10re
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%nlog10tau = nlog10tau
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%nSolzen = nSolzen
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%nzen = nzen
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%nRelaz = nRelaz
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%nRelaz_in_pp = nRelaz_in_pp
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%nRelaz_out_pp = nRelaz_out_pp
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_Relaz_out_pp = delta_Relaz_out_pp
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_Relaz_in_pp = delta_Relaz_in_pp
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz_switch_pp = Relaz_switch_pp

!--- knowing the dimensions, allocate space
   tstatus = 0
   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%zen)) THEN
      ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%zen(nzen), stat=astatus)
      tstatus = tstatus + astatus
   ENDIF
   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen)) THEN
     ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen(nSolzen), stat=astatus)
     tstatus = tstatus + astatus
   ENDIF

   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz)) THEN
      ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz(nRelaz), stat=astatus)
      tstatus = tstatus + astatus
   ENDIF

   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10re)) THEN
      ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10re(nlog10re), stat=astatus)
      tstatus = tstatus + astatus
   ENDIF

   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10tau)) THEN
      ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10tau(nlog10tau), stat=astatus)
      tstatus = tstatus + astatus
   ENDIF

   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%ext_ratio)) THEN
      ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%ext_ratio(nlog10re), stat=astatus)
      tstatus = tstatus + astatus
   ENDIF

   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Qext_ref)) THEN
     ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Qext_ref(nlog10re), stat=astatus)
     tstatus = tstatus + astatus
   ENDIF

   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%sph_alb)) THEN
     ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%sph_alb(nlog10re,nlog10tau), stat=astatus)
      tstatus = tstatus + astatus
   ENDIF

   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%trn)) THEN
     ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%trn(nlog10re,nlog10tau,nSolzen), stat=astatus)
     tstatus = tstatus + astatus
   ENDIF

   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%alb)) THEN
     ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%alb(nlog10re,nlog10tau,nSolzen), stat=astatus)
     tstatus = tstatus + astatus
   ENDIF

   IF (.NOT. ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%ref)) THEN
     ALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%ref(nlog10re,nlog10tau,nSolzen,nzen,nRelaz), stat=astatus)
     tstatus = tstatus + astatus
   ENDIF

   IF (tstatus /= 0) THEN
      PRINT "(a,'Error allocating cloud reflectance lookup table arrays, stopping')",EXE_PROMPT
      STOP
   END IF

!--- initialize
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%zen  = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen  = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz  = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10re  = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10tau  = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%ext_ratio  = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%Qext_ref  = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%sph_alb  = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%trn = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%alb = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%ref = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%rho = 0.0
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%lambda_ref = 0.0


!--- Read in values in the table
istatus = 0

!-- read in vectors
sds_start_1d = 0
sds_stride_1d = 1
sds_edge_1d = nzen

!-- sensor zenith angle
sds_edge_1d = nzen

sds_id = sfselect(sd_id,sfn2index(sd_id,"sensor_zenith_angle"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%zen) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- solar sensor zenith angle
sds_edge_1d = nSolzen
sds_id = sfselect(sd_id,sfn2index(sd_id,"solar_zenith_angle"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- relative azimith angle
sds_edge_1d = nRelaz
sds_id = sfselect(sd_id,sfn2index(sd_id,"relative_azimuth_angle"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- log10 optical depth
sds_edge_1d = nlog10tau
sds_id = sfselect(sd_id,sfn2index(sd_id,"log10_optical_depth"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10tau) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- log10 effective radius
sds_edge_1d = nlog10re
sds_id = sfselect(sd_id,sfn2index(sd_id,"log10_eff_radius"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10re) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- Qext ref
sds_edge_1d = nlog10re
sds_id = sfselect(sd_id,sfn2index(sd_id,"Qext_ref"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%Qext_ref) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- ext_ratio
sds_edge_1d = nlog10re
sds_id = sfselect(sd_id,sfn2index(sd_id,"extinction_ratio"))
istatus = sfrdata(sds_id,sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%ext_ratio) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!--- 2d arrays
sds_start_2d = 0
sds_stride_2d = 1

!--- spherical albedo
sds_dims_2d(1) = nlog10re
sds_dims_2d(2) = nlog10tau
sds_edge_2d(1) = nlog10re
sds_edge_2d(2) = nlog10tau
sds_id = sfselect(sd_id,sfn2index(sd_id,"spherical_albedo"))
istatus = sfrdata(sds_id,sds_start_2d, sds_stride_2d, sds_edge_2d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%sph_alb) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!--- 3d arrays
sds_start_3d = 0
sds_stride_3d = 1

!--- albedo
sds_dims_3d(1) = nlog10re
sds_dims_3d(2) = nlog10tau
sds_dims_3d(3) = nSolzen
sds_edge_3d(1) = nlog10re
sds_edge_3d(2) = nlog10tau
sds_edge_3d(3) = nSolzen
sds_id = sfselect(sd_id,sfn2index(sd_id,"albedo"))
istatus = sfrdata(sds_id,sds_start_3d, sds_stride_3d, sds_edge_3d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%alb) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- transmission
sds_id = sfselect(sd_id,sfn2index(sd_id,"transmission"))
istatus = sfrdata(sds_id,sds_start_3d, sds_stride_3d, sds_edge_3d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%trn) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!--- 5d arrays
sds_start_5d = 0
sds_stride_5d = 1

!--- reflectance
sds_dims_5d(1) = nlog10re
sds_dims_5d(2) = nlog10tau
sds_dims_5d(3) = nSolzen
sds_dims_5d(4) = nzen
sds_dims_5d(5) = nRelaz

sds_edge_5d(1) = nlog10re
sds_edge_5d(2) = nlog10tau
sds_edge_5d(3) = nSolzen
sds_edge_5d(4) = nzen
sds_edge_5d(5) = nRelaz

sds_id = sfselect(sd_id,sfn2index(sd_id,"reflectance"))
istatus = sfrdata(sds_id,sds_start_5d, sds_stride_5d, sds_edge_5d, &
                  cld_ref(ialgo)%channel(ichan)%phase(iphase)%ref) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

    IF (istatus /= 0) THEN
      PRINT "(a,'Error reading sds data from cloud lut files, stopping')",EXE_PROMPT
      STOP
    END IF

!--- Close the lookup table file
    istatus = sfend(sd_id) 

    IF (istatus /= 0) THEN
      PRINT "(a,'Error closing cloud lut files')",EXE_PROMPT
    END IF

!--- assuming everything is equally spaced, compute spacing here
   cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_Solzen = & 
           cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen(2) - & 
           cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen(1)

   cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_zen = & 
           cld_ref(ialgo)%channel(ichan)%phase(iphase)%zen(2) - & 
           cld_ref(ialgo)%channel(ichan)%phase(iphase)%zen(1)

   cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_log10tau = & 
           cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10tau(2) - & 
           cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10tau(1)

   cld_ref(ialgo)%channel(ichan)%phase(iphase)%delta_log10re = & 
           cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10re(2) - & 
           cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10re(1)
 
!=======================================================================
! END of this routine
!=======================================================================
 END SUBROUTINE READ_CLD_REF_LUT_STRUC


!=======================================================================
! Read Cloud Emissivity Lookup Table based on structures
! 
!  
! isat  = wmo satellite identIFication number
! ichan = ABI channel number
!
! sc_ind = spacecraft index
! scinfo(sc_ind)%id = wmo number
!
!=======================================================================
 SUBROUTINE READ_CLD_EMS_LUT_STRUC(ancil_data_path,ialgo,ichan,iphase)

  INTEGER, INTENT(in):: ialgo,ichan
  INTEGER(KIND=int1), INTENT(in):: iphase
  CHARACTER (LEN=*), INTENT(in):: ancil_data_path

  CHARACTER (LEN=128):: lut_file
  INTEGER:: astatus,tstatus
  CHARACTER (LEN=2):: chan_string
  CHARACTER (LEN=3):: phase_string
  CHARACTER (LEN=7):: sat_string
  INTEGER:: nlog10re
  INTEGER:: nlog10tau
  INTEGER:: nzen

!-- hdf
  INTEGER:: sfstart, sfrcatt,sfn2index, sffattr, sfendacc,sfend,sfrdata,sfselect
  INTEGER:: sd_id, sds_id, istatus

  INTEGER, PARAMETER:: sds_rank_1d = 1
  INTEGER, DIMENSION(sds_rank_1d):: sds_start_1d, sds_edge_1d, sds_stride_1d

  INTEGER, PARAMETER:: sds_rank_3d = 3
  INTEGER, DIMENSION(sds_rank_3d):: sds_start_3d, sds_edge_3d, sds_stride_3d, sds_dims_3d

!--- construct appropriate file name - use actual satellite channel, not abi equivalent
! write(chan_string, '(I1.1)') scinfo(sc_ind)%ch_flg(ichan)

  IF (ichan == 2) THEN
          chan_string = "1"
  ELSEIF (ichan == 5) THEN
          chan_string = "3a"
  ELSEIF (ichan == 7) THEN
          chan_string = "3b"
  ENDIF

  phase_string = "wat"
  IF (iphase == ICE_PHASE) THEN
    phase_string = "ice"
  END IF
  IF (iphase == MIXED_PHASE) THEN
    phase_string = "wat"                  !note assumed water tables apply to mixed phase clouds
  END IF

!--- make filenames
  IF (phase_string /= "ice") THEN

    IF (Sc_Id_Internal < 4) THEN
        IF (Sc_Id_Internal == 2) THEN
          sat_string = "METOP-A"
        ENDIF
        IF (Sc_Id_Internal == 3) THEN
          sat_string = "METOP-B"
        ENDIF
    ENDIF
    IF (Sc_Id_Internal < 10 .and. Sc_Id_Internal > 4) THEN
         write(unit=sat_string,fmt="(a6,i1)") "NOAA-0",Sc_Id_Internal
    ENDIF
    IF (Sc_Id_Internal >= 10) THEN
         write(unit=sat_string,fmt="(a5,i2)")  "NOAA-",Sc_Id_Internal
    ENDIF

    lut_file = TRIM(sat_string)//"_ch"//TRIM(chan_string)//"_ems_lut_"//phase_string//"_cld.hdf"
  ELSE
    sat_string = "AVHRR"
    lut_file = TRIM(sat_string)//"_ch"//TRIM(chan_string)//"_ems_lut_mix1_"//phase_string//"_cld.hdf"
  ENDIF

!-- should include scinfo(sc_ind)%id - the wmo number
! lut_file = TRIM(scinfo(sc_ind)%name)//"_ch"//TRIM(chan_string)//"_ems_lut_"//phase_string//"_cld.hdf"

!--- Open lookup table
 sd_id = sfstart( &
                 TRIM(ancil_data_path)//TRIM(lut_file), &
                 DFACC_READ)

 IF (sd_id < 0) THEN
     PRINT *, EXE_PROMPT, "Error open cloud ems table ", trim(ancil_data_path)//trim(lut_file)
     stop
  ENDIF
!--- global attributes 
   istatus = 0

   istatus = sfrcatt(sd_id, &
                     sffattr(sd_id,"HEADER"),  &
                     cld_ems(ialgo)%channel(ichan)%phase(iphase)%header)

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
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%nlog10re = nlog10re
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%nlog10tau = nlog10tau
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%nzen = nzen

   tstatus = 0
   ALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%zen(nzen), stat=astatus)
   tstatus = tstatus + astatus
   ALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10re(nlog10re), stat=astatus)
   tstatus = tstatus + astatus
   ALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10tau(nlog10tau), stat=astatus)
   tstatus = tstatus + astatus
   ALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%ext_ratio(nlog10re), stat=astatus)
   tstatus = tstatus + astatus
   ALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%Qext_ref(nlog10re), stat=astatus)
   tstatus = tstatus + astatus
   ALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%emiss(nlog10re,nlog10tau,nzen), stat=astatus)
   tstatus = tstatus + astatus
   ALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%trans(nlog10re,nlog10tau,nzen), stat=astatus)
   tstatus = tstatus + astatus

    IF (tstatus /= 0) THEN
      PRINT "(a,'Error allocating cloud emissivity lookup table arrays')",EXE_PROMPT
      STOP
    END IF

!--- initialize
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%zen  = 0.0
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10re  = 0.0
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10tau  = 0.0
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%ext_ratio  = 0.0
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%Qext_ref  = 0.0
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%emiss = 0.0
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%trans = 0.0
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%rho = 0.0
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%lambda_ref = 0.0

!--- Read in values in the table

!-- read in vectors
sds_start_1d = 0
sds_stride_1d = 1
sds_edge_1d = nzen

!-- sensor zenith angle
sds_edge_1d = nzen
sds_id = sfselect(sd_id,sfn2index(sd_id,"sensor_zenith_angle"))
istatus = sfrdata(sds_id, sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ems(ialgo)%channel(ichan)%phase(iphase)%zen) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- log10 optical depth
sds_edge_1d = nlog10tau
sds_id = sfselect(sd_id,sfn2index(sd_id,"log10_optical_depth"))
istatus = sfrdata(sds_id, sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10tau) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- log10 effective radius
sds_edge_1d = nlog10re
sds_id = sfselect(sd_id,sfn2index(sd_id,"log10_eff_radius"))
istatus = sfrdata(sds_id, sds_start_1d, sds_stride_1d, sds_edge_1d, &
                  cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10re) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!--- 3d arrays
sds_start_3d = 0
sds_stride_3d = 1

!--- emissivity
sds_dims_3d(1) = nlog10re
sds_dims_3d(2) = nlog10tau
sds_dims_3d(3) = nzen
sds_edge_3d(1) = nlog10re
sds_edge_3d(2) = nlog10tau
sds_edge_3d(3) = nzen

sds_id = sfselect(sd_id,sfn2index(sd_id,"cloud_emissivity"))
istatus = sfrdata(sds_id, sds_start_3d, sds_stride_3d, sds_edge_3d, &
                  cld_ems(ialgo)%channel(ichan)%phase(iphase)%emiss) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus

!-- transmission
sds_id = sfselect(sd_id,sfn2index(sd_id,"cloud_transmission"))

istatus = sfrdata(sds_id, sds_start_3d, sds_stride_3d, sds_edge_3d, &
                  cld_ems(ialgo)%channel(ichan)%phase(iphase)%trans) +  &
                  istatus
istatus = sfendacc(sds_id) + istatus


!--- Close the lookup table file
  istatus = sfend(sd_id) + istatus

!--- check status of read
  IF (istatus /= 0) THEN
     PRINT *, "Error reading emissivity table for channel, phase = ",  &
              ichan, iphase
  ENDIF

!--- assuming everything is equally spaced, compute spacing here
   cld_ems(ialgo)%channel(ichan)%phase(iphase)%delta_zen = &
           cld_ems(ialgo)%channel(ichan)%phase(iphase)%zen(2) - &
           cld_ems(ialgo)%channel(ichan)%phase(iphase)%zen(1)

   cld_ems(ialgo)%channel(ichan)%phase(iphase)%delta_log10tau = &
           cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10tau(2) - &
           cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10tau(1)

   cld_ems(ialgo)%channel(ichan)%phase(iphase)%delta_log10re = &
           cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10re(2) - &
           cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10re(1)

!=======================================================================
! END of this routine
!=======================================================================
 END SUBROUTINE READ_CLD_EMS_LUT_STRUC

 SUBROUTINE ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE(nalgo)
   INTEGER, INTENT(in):: nalgo
   INTEGER:: astatus

!--- ALLOCATE main structures

     ALLOCATE(cld_ref(nalgo),stat=astatus)
     IF (astatus /= 0) THEN
       PRINT "(a,'Not enough memory to ALLOCATE cld reflectance lut structure.')",EXE_PROMPT
       STOP
     END IF
     cld_ref(:)%flag = sym%NO

     ALLOCATE(cld_ems(nalgo),stat=astatus)
     IF (astatus /= 0) THEN
       PRINT "(a,'Not enough memory to ALLOCATE cld emissivity lut structure.')",EXE_PROMPT
       STOP
     END IF

     cld_ems(:)%flag = sym%NO

 END SUBROUTINE ALLOCATE_MEMORY_FOR_MAIN_TABLE_STRUCTURE

!======================================================================
! Deallocate threshold arrays on exit of last segment
!======================================================================
  SUBROUTINE DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES(ialgo)

   INTEGER, INTENT(in):: ialgo
   INTEGER:: dstatus        !deallocate status 
   INTEGER:: tstatus        !accumulated deallocate status
   INTEGER:: ichan          !channel index
   INTEGER:: iphase         !phase index

    PRINT *, EXE_PROMPT, " Destroying Cloud LUTs" 

!--------------------------------------------------------------
!--- reflectance tables
!--------------------------------------------------------------
    dstatus = 0
    tstatus = 0
     DO ichan = 1, nchan_max
      DO iphase = 1,nphases
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%zen)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%zen,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Solzen,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Relaz,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10tau)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10tau,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10re)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%log10re,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%ext_ratio)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%ext_ratio,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Qext_ref)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%Qext_ref,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%sph_alb)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%sph_alb,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%trn)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%trn,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%alb)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%alb,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ref(ialgo)%channel(ichan)%phase(iphase)%ref)) THEN
         DEALLOCATE(cld_ref(ialgo)%channel(ichan)%phase(iphase)%ref,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
      END DO
     END DO

!--- check status of deallocating
     IF (tstatus /= 0) THEN
      PRINT "(a,'Error deallocating cloud reflectance lookup table arrays')",EXE_PROMPT
      STOP
     ENDIF

!--------------------------------------------------------------
!--- emissivity tables
!--------------------------------------------------------------
    dstatus = 0
    tstatus = 0
     DO ichan = 1, nchan_max
      DO iphase = 1,nphases

       IF (ALLOCATED(cld_ems(ialgo)%channel(ichan)%phase(iphase)%zen)) THEN
         DEALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%zen,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10tau)) THEN
         DEALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10tau,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10re)) THEN
         DEALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%log10re,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ems(ialgo)%channel(ichan)%phase(iphase)%ext_ratio)) THEN
         DEALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%ext_ratio,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ems(ialgo)%channel(ichan)%phase(iphase)%Qext_ref)) THEN
         DEALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%Qext_ref,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ems(ialgo)%channel(ichan)%phase(iphase)%trans)) THEN
         DEALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%trans,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF
       IF (ALLOCATED(cld_ems(ialgo)%channel(ichan)%phase(iphase)%emiss)) THEN
         DEALLOCATE(cld_ems(ialgo)%channel(ichan)%phase(iphase)%emiss,stat=dstatus)
         tstatus = tstatus + dstatus
       ENDIF

      END DO

     END DO

!--- check status of deallocating
     IF (tstatus /= 0) THEN
      PRINT "(a,'Error deallocating cloud reflectance lookup table arrays')",EXE_PROMPT
      STOP
     ENDIF

   END SUBROUTINE DEALLOCATE_MEMORY_FOR_TABLE_STRUCTURES

!----------------------------------------------------------------------
!--- routine to deallocate main cloud lut structure
!----------------------------------------------------------------------
   SUBROUTINE DEALLOCATE_MEMORY_FOR_MAIN_CLD_LUT_STRUCTURES()
    INTEGER:: dstatus      !deallocate status
    INTEGER:: tstatus      !accumulated status

!--- initialize status flags
    dstatus = 0
    tstatus = 0

!--- reflectance
    DEALLOCATE(cld_ref,stat=dstatus)
    tstatus = tstatus + dstatus
    IF (tstatus /= 0) THEN
      PRINT "(a,'Error deallocating cloud reflectance lookup table arrays')",EXE_PROMPT
      STOP
    ENDIF

!--- emissivity
    DEALLOCATE(cld_ems,stat=dstatus)
    tstatus = tstatus + dstatus
    IF (tstatus /= 0) THEN
      PRINT "(a,'Error deallocating cloud emissivity lookup table arrays')",EXE_PROMPT
      STOP
    ENDIF

 END SUBROUTINE DEALLOCATE_MEMORY_FOR_MAIN_CLD_LUT_STRUCTURES

END MODULE CLOUD_TAU_RE_SOLAR_MODULE
