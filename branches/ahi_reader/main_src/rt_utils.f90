!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: rt_utils.f90 (src)
!       RT_UTILITIES (program)
!
! PURPOSE: Perform needed Radiative Transfer Functions
!
! DESCRIPTION: CLAVR-x uses much radiative transfer.  This module stores and serves
!              the RT variables to other modules.  It serves as the interface
!              to external RT models (PFAAST) and holds convenient RT-specific
!              routines 
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
! Description of RTM Structure Members given in rtm_common.f90
!
! CLAVR-x has 45 channels.  
! Channels 1-36 are MODIS
! Channels 37-38 are ABI channels not on MODIS
! Channels 39-43 are the VIIRS I-bands
! Channel 44 is the VIIRS DNB
! Channel 45 is the IFF Pseudo MODIS ch 33
!
! Not all members of the RTM structure are populated for all channels.
! However, all members are allocated for any active cell
!
! Here is the current implementation
!
! There are 6 types of channels (FIX THIS)
! 1. MODIS IR-only channels = 21-36 excluding 26. 
! 2. Non-MODIS IR-only channels  37-38
! 3. MODIS (Solar + IR) channels = Channel 20
! 4. Supported Solar Channels = Channels 1,2,5,6,7 and 44
! 5. Unsupported Solar Channels = Channels 3,4,8-19,26,39-41
! 6. Unsupported IR Channels = Channels 42 and 43
!
! For each type described above, the following profiles are made:
! 1. Rad_Atm_Prof, Trans_Atm_Profile, Rad_BB_Cloud_Profile
! 2. All Profiles
! 3. Trans_Atm_Profile, Trans_Atm_Solar_Profile, Trans_Atm_Total_Profile
! 4. No Profiles
! 5. No Profiles
! 
! Viirs I-band support is limited (channels 39-43) and is being developed. No
! pixel-level toa clear-sky fields are generated for the I-band variables yet.
! 
!--------------------------------------------------------------------------------------
module RT_UTILITIES
 
   use CONSTANTS, only: &
      Real4 &
      , Int4 &
      , Int1 &
      , Sym &
      , Dtor &
      , Missing_Value_Real4  &
      , Exe_Prompt &
      , G &
      , PI &
      , SOLAR_OBS_TYPE &
      , LUNAR_OBS_TYPE &
      , MIXED_OBS_TYPE &
      , THERMAL_OBS_TYPE
      
   use NWP_COMMON, only : &
        Nlevels_Nwp &
      , T_Prof_Nwp &
      , P_Std_Nwp &
      , Rh_Prof_Nwp &
      , Tpw_Prof_Nwp &
      , Sfc_Level_Nwp &
      , Tmpair_Nwp &
      , Rhsfc_Nwp &
      , Psfc_Nwp &
      , Z_Prof_Nwp &
      , Ozone_Prof_Nwp &
      , Delta_T_Inversion &
      , P_Inversion_Min &     
      , P_Trop_Nwp &
      , P_Trop_Nwp 
      
   use PIXEL_COMMON, only: &
      Nwp_Opt &
      , Sensor &
      , Ch &
      , Geo &
      , Sfc &
      , Image &
      , Bad_Pixel_Mask &
      , Zen_Idx_Rtm &
      , I_Nwp &
      , J_Nwp &
      , I_Nwp_X &
      , J_Nwp_X &
      , Lon_Nwp_Fac &
      , Lat_Nwp_Fac &
      , Sfc_Level_Rtm_Pixel &
      , Solar_Rtm &
      , Space_Mask &
      , Use_Sst_Anal &
      , Sst_Anal &
      , Tsfc_Nwp_Pix &
      , Trans_Atm_Ch20_Solar_Rtm &
      , Trans_Atm_Ch20_Solar_Total_Rtm &
      , rad_clear_Ch20_Solar_Rtm &
      , Bt_Clear_Ch20_Solar_Rtm &
      , Ems_Ch20_Clear_Solar_Rtm &
      , Ems_Ch20_Clear_Rtm &
      , Ems_Ch20_Clear_Solar_Sfc_Rtm &
      , Beta_11um_12um_Tropo_Rtm &
      , Beta_11um_85um_Tropo_Rtm &
      , Beta_11um_67um_Tropo_Rtm &
      , Beta_11um_133um_Tropo_Rtm &
      , Pixel_Local_Time_Hours &
      , Ancil_Data_Dir
      
   use NUMERICAL_ROUTINES , only: &
         VAPOR &
       , COMPUTE_TIME_HOURS &
       , LOCATE
      
   use PLANCK, only: &
        PLANCK_RAD_FAST &
      , PLANCK_TEMP_FAST 

    use calibration_constants, only: &
       SOLAR_CH20_NU
      
   ! use SURFACE_PROPERTIES, only: 
   
   use RTM_COMMON, only: &
       Nlevels_Rtm &
      , Rtm &
      , P_Std_Rtm &
      , T_Std_Rtm &
      , Wvmr_Std_Rtm &
      , Ozmr_Std_Rtm 

    use CX_PFAAST_MOD, only: &
       COMPUTE_TRANSMISSION_PFAAST
   implicit none
   
   private:: EMISSIVITY, &
             BETA_RATIO
   
    
   public::  &
            SETUP_SOLAR_RTM, &
            MAP_NWP_RTM,  &
            CREATE_TEMP_NWP_VECTORS,  &
            DESTROY_TEMP_NWP_VECTORS, &      
            GET_PIXEL_NWP_RTM, &
            ALLOCATE_RTM, &
            DEALLOCATE_RTM, &
            DEALLOCATE_RTM_VARS, &
            DEALLOCATE_RTM_CELL
            
    integer, parameter, public:: Rtm_Nvzen = 50       

    real, parameter :: Co2_Ratio = 380.0 !in ppmv

    integer, parameter :: Chan_Idx_Min = 1
    integer, parameter :: Chan_Idx_Max = 45

    real(kind=real4),  save, dimension(Chan_Idx_Min:Chan_Idx_Max):: Gamma_Trans_Factor

    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Trans_Atm_Prof
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Trans_Atm_Solar_Prof
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Trans_Atm_Total_Prof
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Trans_Atm_Prof_Prev
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Trans_Atm_Solar_Prof_Prev
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Trans_Atm_Total_Prof_Prev
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Rad_Atm_Prof
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Rad_Atm_Dwn_Prof
    real, dimension(NLevels_Rtm,Chan_Idx_Min:Chan_Idx_Max),  save:: Rad_BB_Cloud_Prof

    integer, parameter:: Ilon_Stride = 0
    integer, parameter:: Ilat_Stride = 0
    integer, parameter::  Ivza_Stride = 0
    integer, dimension(NLevels_Rtm) :: k_Rtm_Nwp
    real, dimension(NLevels_Rtm) :: x_Rtm_Nwp
    real, dimension(NLevels_Rtm),  save::  &
                     T_Prof_Rtm,  &
                     Z_Prof_Rtm,  &
                     Wvmr_Prof_Rtm,  &
                     Ozmr_Prof_Rtm, &
                     Tpw_Prof_Rtm, &
                     Trans_Prof_Rtm

    integer, dimension(:), allocatable, save:: k_Nwp_Rtm
    real, dimension(:), allocatable, save:: x_Nwp_Rtm
    real, dimension(:), allocatable,  save:: Wvmr_Nwp

    !----- PFAAST Specific Settings (all private to this module)
    integer, dimension(Chan_Idx_Min:Chan_Idx_Max),  save:: Rtm_Chan_Idx
    character(len=20), dimension(Chan_Idx_Min:Chan_Idx_Max),  save:: Pfaast_Name
    integer,  save:: Sc_Id_Rtm
    character(len=20),  save:: Sc_Name_Rtm
   
    real, parameter::  Rtm_Vza_Binsize = 0.02

contains

   !====================================================================
   ! subroutine Name: CREATE_TEMP_NWP_VECTORS
   !
   ! Function:
   !   create and initialize NWP vectors used for RTM calcs
   !
   !====================================================================
   subroutine CREATE_TEMP_NWP_VECTORS()

      allocate( Wvmr_Nwp(NLevels_Nwp), &
           k_Nwp_Rtm(NLevels_Nwp), &
           x_Nwp_Rtm(NLevels_Nwp))

      Wvmr_Nwp = 0.0
      k_Nwp_Rtm = 0.0
      x_Nwp_Rtm = 0.0

   end subroutine CREATE_TEMP_NWP_VECTORS
   
   !====================================================================
   ! subroutine Name: DESTROY_TEMP_NWP_VECTORS
   !
   ! Function:
   !   destroy and initialize NWP vectors used for RTM calcs
   !
   !====================================================================
   subroutine DESTROY_TEMP_NWP_VECTORS()

      deallocate(Wvmr_Nwp,k_Nwp_Rtm, x_Nwp_Rtm)

   end subroutine DESTROY_TEMP_NWP_VECTORS

   !====================================================================
   ! subroutine Name: MAP_NWP_RTM
   !
   ! Function:
   !   develops the mapping of the NWP and RTM profile Levels
   !
   !====================================================================
   subroutine MAP_NWP_RTM(Num_Levels_Nwp_Profile, &
                        Press_Profile_Nwp, &
                        Num_Levels_Rtm_Profile, &
                        Press_Profile_Rtm)

      integer, intent(in):: Num_Levels_Nwp_Profile
      integer, intent(in):: Num_Levels_Rtm_Profile
      real, intent(in), dimension(:)::  Press_Profile_Nwp
      real, intent(in), dimension(:)::  Press_Profile_Rtm

      integer:: k
      integer:: k_temp
    
     
      !--- map NWP Levels to RTM Levels
      do k = 1, Num_Levels_Nwp_Profile
     
         !--- locate the nwp Level within the Rtm Levels 
         !--- P_Std_Nwp(k) should fall between Rtm Levels k_temp and k_temp +1
         call LOCATE(Press_Profile_Rtm, Num_Levels_Rtm_Profile, Press_Profile_Nwp(k), k_temp)
        
         k_Nwp_Rtm(k) = min(Num_Levels_Rtm_Profile-1,max(1,k_temp))

         !-- compute linear weighting factor
         x_Nwp_Rtm(k) = (Press_Profile_Nwp(k) - Press_Profile_Rtm(k_Nwp_Rtm(k))) / &
                      (Press_Profile_Rtm(k_Nwp_Rtm(k)+1) - Press_Profile_Rtm(k_Nwp_Rtm(k))) 
      end do

      !--- map RTM Levels to NWP Levels
      do k = 1, Num_Levels_Rtm_Profile
      
         !--- locate the Rtm Level within the Nwp Levels 
         !--- Press_Profile_Rtm(k) should fall between Nwp Levels k_temp and k_temp +1
         call LOCATE(Press_Profile_Nwp, Num_Levels_Nwp_Profile, Press_Profile_Rtm(k), k_temp)
         k_Rtm_Nwp(k) = min(Num_Levels_Nwp_Profile-1,max(1,k_temp))

         !-- compute linear weighting factor
         x_Rtm_Nwp(k) = (Press_Profile_Rtm(k) - Press_Profile_Nwp(k_Rtm_Nwp(k))) / &
                      (Press_Profile_Nwp(k_Rtm_Nwp(k)+1) - Press_Profile_Nwp(k_Rtm_Nwp(k)))
     
      end do 
      
   end subroutine MAP_NWP_RTM
   
   !====================================================================
   ! subroutine Name: CONVERT_ATMOS_PROF_NWP_RTM
   !
   ! Description:
   ! This routine interpolate the NWP profiles to profiles with the
   ! vertical spacing defined by the RTM model.  It operates on profiles
   ! stored in this module
   !
   ! INPUTS:
   !
   ! Highest_Level_Rtm_Nwp - highest Rtm Level that is below the highest nwp Level
   ! Lowest_Level_Rtm_Nwp - lowest Rtm Level that is above the lowest nwp Level
   ! Sfc_Level_Rtm - lowest Rtm Level above the surface
   ! P_near_Sfc_Nwp - lowest standard nwp Level above surface pressure
   !
   !====================================================================
   subroutine CONVERT_ATMOS_PROF_NWP_RTM(Num_Levels_NWP_Profile, &
                                       Surface_Level_Nwp, &
                                       Air_Temperature_Nwp, &
                                       Surface_Rh_Nwp, &
                                       Surface_Pressure_Nwp, &
                                       Press_Profile_Nwp, &
                                       T_Profile_Nwp, &
                                       Z_Profile_Nwp, &
                                       Wvmr_Profile_Nwp, &
                                       Ozmr_Profile_Nwp, &
                                       Num_Levels_Rtm_Profile, &
                                       Press_Profile_Rtm, &
                                       T_Profile_Rtm, &
                                       Z_Profile_Rtm, &
                                       Wvmr_Profile_Rtm, &
                                       Ozmr_Profile_Rtm, &
                                       T_Std_Profile_Rtm, &
                                       Wvmr_Std_Profile_Rtm, &
                                       Ozmr_Std_Profile_Rtm)
                                        

      integer, intent(in):: Num_Levels_Nwp_Profile
      integer(kind=int1), intent(in):: Surface_Level_Nwp 
      real, intent(in):: Air_Temperature_Nwp 
      real, intent(in):: Surface_Rh_Nwp 
      real, intent(in):: Surface_Pressure_Nwp 
      real, intent(in), dimension(:):: Press_Profile_Nwp 
      real, intent(in), dimension(:):: T_Profile_Nwp 
      real, intent(in), dimension(:):: Z_Profile_Nwp 
      real, intent(in), dimension(:):: Wvmr_Profile_Nwp 
      real, intent(in), dimension(:):: Ozmr_Profile_Nwp 
      integer, intent(in):: Num_Levels_Rtm_Profile
      real, intent(in), dimension(:):: Press_Profile_Rtm
      real, intent(out), dimension(:):: T_Profile_Rtm
      real, intent(out), dimension(:):: Z_Profile_Rtm
      real, intent(out), dimension(:):: Wvmr_Profile_Rtm
      real, intent(out), dimension(:):: Ozmr_Profile_Rtm
      real, intent(in), dimension(:):: T_Std_Profile_Rtm
      real, intent(in), dimension(:):: Wvmr_Std_Profile_Rtm
      real, intent(in), dimension(:):: Ozmr_Std_Profile_Rtm

      integer:: k
      integer:: Lowest_Level_Rtm_Nwp
      integer:: Highest_Level_Rtm_Nwp
      integer:: Sfc_Level_Rtm
      real:: dT_dP_near_Sfc
      real:: dWvmr_dP_near_Sfc
      real:: dZ_dP_near_Sfc
      real:: P_near_Sfc_Nwp
      real:: Wvmr_Sfc
      real:: es
      real:: e
      real:: T_Offset

      !--- initialize indices
      Lowest_Level_Rtm_Nwp = Num_Levels_Rtm_Profile
      Highest_Level_Rtm_Nwp = 1
      Sfc_Level_Rtm = Num_Levels_Rtm_Profile
      P_Near_Sfc_Nwp = Press_Profile_Nwp(Surface_Level_Nwp)
   
      !--- make Wvmr at sfc
      es = VAPOR(Air_Temperature_Nwp)
      e = Surface_Rh_Nwp * es / 100.0
      Wvmr_Sfc = 1000.0*0.622 * (e / (Surface_Pressure_Nwp - e))  !(g/kg)

      !--- determine some critical Levels in the Rtm profile
      do k = 1, Num_Levels_Rtm_Profile
         if (Press_Profile_Rtm(k) > Press_Profile_Nwp(1)) then
            Highest_Level_Rtm_Nwp = k
            exit
         endif
      enddo

      do k = Num_Levels_Rtm_Profile,1,-1
         if (Press_Profile_Rtm(k) < Surface_Pressure_Nwp) then
            Sfc_Level_Rtm = k
            exit
         endif
      enddo

      do k = Num_Levels_Rtm_Profile,1,-1
         if (Press_Profile_Rtm(k) < P_Near_Sfc_Nwp) then
            Lowest_Level_Rtm_Nwp = k
            exit
         endif
      enddo

      !--- compute T and Wvmr lapse rate near the surface
      dT_dP_near_Sfc = 0.0
      dZ_dP_near_Sfc = 0.0
      dWvmr_dP_near_Sfc = 0.0

      if (Surface_Pressure_Nwp /= Press_Profile_Nwp(Num_Levels_Nwp_Profile)) then
         dT_dP_near_Sfc =  &
             (T_Profile_Nwp(Surface_Level_Nwp) - Air_Temperature_Nwp)/ &
             (Press_Profile_Nwp(Surface_Level_Nwp) - Surface_Pressure_Nwp)
         dWvmr_dP_near_Sfc =  &
             (Wvmr_Profile_Nwp(Surface_Level_Nwp) - Wvmr_Sfc)/ &
             (Press_Profile_Nwp(Surface_Level_Nwp) - Surface_Pressure_Nwp)
      else
         dT_dP_near_Sfc =  &
             (T_Profile_Nwp(Surface_Level_Nwp-1) - Air_Temperature_Nwp)/ &
             (Press_Profile_Nwp(Surface_Level_Nwp-1) - Surface_Pressure_Nwp)
         dWvmr_dP_near_Sfc =  &
             (Wvmr_Profile_Nwp(Surface_Level_Nwp-1) - Wvmr_Sfc)/ &
             (Press_Profile_Nwp(Surface_Level_Nwp-1) - Surface_Pressure_Nwp)
      end if
   
      if (Press_Profile_Nwp(Surface_Level_Nwp-1) /= Press_Profile_Nwp(Surface_Level_Nwp)) then
         dZ_dP_near_Sfc =  &
             (Z_Profile_Nwp(Surface_Level_Nwp-1) - Z_Profile_Nwp(Surface_Level_Nwp))/ &
             (Press_Profile_Nwp(Surface_Level_Nwp-1) - Press_Profile_Nwp(Surface_Level_Nwp))
       
      end if

      ! --- compute temperature offset between standard and NWP profiles at top
      !--- this will be added to the standard profiles
      T_Offset = T_Profile_Nwp(1) - T_Std_Profile_Rtm(Highest_Level_Rtm_Nwp)

      !--- for Rtm Levels above the highest nwp Level, use Rtm standard
      !values
      do k = 1,Highest_Level_Rtm_Nwp-1
          Z_Profile_Rtm(k) = Z_Profile_Nwp(1)
          T_Profile_Rtm(k) = T_Std_Profile_Rtm(k) + T_Offset
          Wvmr_Profile_Rtm(k) = Wvmr_Std_Profile_Rtm(k)
          Ozmr_Profile_Rtm(k) = Ozmr_Std_Profile_Rtm(k)
      end do

      !--- Rtm Levels within standard nwp Levels above the surface
      do k = Highest_Level_Rtm_Nwp, Lowest_Level_Rtm_Nwp
          T_Profile_Rtm(k) = T_Profile_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (T_Profile_Nwp(k_Rtm_Nwp(k)+1) - T_Profile_Nwp(k_Rtm_Nwp(k)))
          Z_Profile_Rtm(k) = Z_Profile_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (Z_Profile_Nwp(k_Rtm_Nwp(k)+1) - Z_Profile_Nwp(k_Rtm_Nwp(k)))
          Wvmr_Profile_Rtm(k) = Wvmr_Profile_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (Wvmr_Profile_Nwp(k_Rtm_Nwp(k)+1) - Wvmr_Profile_Nwp(k_Rtm_Nwp(k)))

          Ozmr_Profile_Rtm(k) = Ozmr_Profile_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (Ozmr_Profile_Nwp(k_Rtm_Nwp(k)+1) - Ozmr_Profile_Nwp(k_Rtm_Nwp(k)))
      end do

      !--- Rtm Levels that are below the lowest nwp Level but above the surface
      do k = Lowest_Level_Rtm_Nwp+1, Sfc_Level_Rtm
           T_Profile_Rtm(k) = Air_Temperature_Nwp + dT_dP_near_Sfc * &
                      (Press_Profile_Rtm(k) - Surface_Pressure_Nwp)
           Wvmr_Profile_Rtm(k) = Wvmr_Sfc + dWvmr_dP_near_Sfc * &
                      (Press_Profile_Rtm(k) - Surface_Pressure_Nwp)
           Z_Profile_Rtm(k) = dZ_dP_near_Sfc * &
                      (Press_Profile_Rtm(k) - Surface_Pressure_Nwp)
           Ozmr_Profile_Rtm(k) = Ozmr_Profile_Nwp(Num_Levels_Nwp_Profile)
      end do

      !--- Rtm Levels below the surface
      do k = Sfc_Level_Rtm +1, Num_Levels_Rtm_Profile
         T_Profile_Rtm(k) = Air_Temperature_Nwp
         Z_Profile_Rtm(k) = dZ_dP_near_Sfc * (Press_Profile_Rtm(k) - Surface_Pressure_Nwp)
         Wvmr_Profile_Rtm(k) = Wvmr_Sfc
         Ozmr_Profile_Rtm(k) = Ozmr_Std_Profile_Rtm(k)
      end do

      !--- if using NCEP reanalysis which has no ozone profile, use default
      if (Nwp_Opt == 2) then
         Ozmr_Profile_Rtm = Ozmr_Std_Profile_Rtm
      end if

   end subroutine CONVERT_ATMOS_PROF_NWP_RTM


   !====================================================================
   ! subroutine Name: COMPUTE_CLEAR_RAD_PROFILES_RTM
   !
   ! Function:
   ! Computes clear sky radiance profiles
   !
   !====================================================================
   subroutine COMPUTE_CLEAR_RAD_PROFILES_RTM()

      integer:: Chan_Idx
      integer:: Lev_Idx
      real:: T_mean
      real:: B_mean
      real:: B_Level
      real:: Opd_Layer
      real:: Trans_Layer
      real:: Trans_Total


      !--- upwelling profiles
      Rad_Atm_Prof = Missing_Value_Real4
      Rad_BB_Cloud_Prof = Missing_Value_Real4

      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max

         if (Ch(Chan_Idx) %Obs_Type /= THERMAL_OBS_TYPE .and. &
             Ch(Chan_Idx)%Obs_Type /= MIXED_OBS_TYPE) cycle
         !if (Chan_Idx < 20) cycle
         !if (Chan_Idx == 26) cycle
         !if (Chan_Idx > 38) cycle
         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle

         Rad_Atm_Prof(1,Chan_Idx) = 0.0
         Rad_BB_Cloud_Prof(1,Chan_Idx) = 0.0

         do Lev_Idx = 2, NLevels_Rtm
    
            T_mean = 0.5*(T_Prof_Rtm(Lev_Idx-1) + T_Prof_Rtm(Lev_Idx))

            B_mean = PLANCK_RAD_FAST(Chan_Idx,T_mean)

            Rad_Atm_Prof(Lev_Idx,Chan_Idx) = Rad_Atm_Prof(Lev_Idx-1,Chan_Idx) +  &
              (Trans_Atm_Prof(Lev_Idx-1,Chan_Idx) - Trans_Atm_Prof(Lev_Idx,Chan_Idx)) * B_mean

            B_Level = PLANCK_RAD_FAST(Chan_Idx,T_Prof_Rtm(Lev_Idx))

            Rad_BB_Cloud_Prof(Lev_Idx,Chan_Idx) = Rad_Atm_Prof(Lev_Idx,Chan_Idx) +  &
                                          (Trans_Atm_Prof(Lev_Idx,Chan_Idx) * B_Level)

         end do
      end do

      !--- downwelling profiles
      Rad_Atm_Dwn_Prof = Missing_Value_Real4

      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max
         if (Chan_Idx /= 31) cycle
         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle

         Trans_Total = 1.0
         Rad_Atm_Dwn_Prof(1,Chan_Idx) = 0.0

         do Lev_Idx = 2, Nlevels_Rtm

            T_mean = 0.5*(T_Prof_Rtm(Lev_Idx) + T_Prof_Rtm(Lev_Idx-1))
            B_mean = PLANCK_RAD_FAST(Chan_Idx,T_mean)

            Opd_Layer = -1.0*log(Trans_Atm_Prof(Lev_Idx,Chan_Idx)/Trans_Atm_Prof(Lev_Idx-1,Chan_Idx)) 
            Opd_Layer = max(0.0,Opd_Layer)

            Trans_Layer = exp(-1.0*Opd_Layer)

            Rad_Atm_Dwn_Prof(Lev_Idx,Chan_Idx) = Trans_Total * Rad_Atm_Dwn_Prof(Lev_Idx-1,Chan_Idx) +  &
                                          (1.0-Trans_Layer) * B_mean

            Trans_Total = Trans_Total * Trans_Layer

         end do

      end do

   end subroutine COMPUTE_CLEAR_RAD_PROFILES_RTM

   !====================================================================
   ! subroutine Name: GET_PIXEL_NWP_RTM
   !
   ! Function:
   ! Calculates the PFAAST transmittance profiles for a given segment.
   !
   !====================================================================
   subroutine GET_PIXEL_NWP_RTM(Line_Idx_Min,Num_Lines)
      integer, intent(in):: Line_Idx_Min
      integer, intent(in):: Num_Lines
      integer:: Elem_Idx
      integer:: Line_Idx
      integer:: Sfc_Level_Idx
      integer:: Lon_Idx
      integer:: Lat_Idx
      integer:: Lon_Idx_x
      integer:: Lat_Idx_x
      integer:: Lon_Idx_Prev
      integer:: Lat_Idx_Prev
      integer:: Zen_Idx_Prev
      real:: Lat_x
      real:: Lon_x
      real:: Prof_Weight
      real:: Satzen_Mid_Bin
      integer:: Zen_Idx
      integer:: Error_Status
      integer:: Chan_Idx
      integer:: Chan_Idx_For_Pfaast

      real(kind=real4),save :: Segment_Time_Point_Seconds_Temp = 0
      real(kind=real4) :: Start_Time_Point_Hours_Temp
      real(kind=real4) :: End_Time_Point_Hours_Temp

      !--------------------------------------------------------------------------
      ! Compute Gamma_Trans_Factor for radiance bias
      !--------------------------------------------------------------------------
      call COMPUTE_GAMMA_FACTOR()
     

      Lon_Idx_Prev = 0
      Lat_Idx_Prev = 0
      Zen_Idx_Prev = 0

      !--- loop over pixels in segment
      line_loop: do Line_Idx = Line_Idx_Min, Num_Lines + Line_Idx_Min - 1
         element_loop: do Elem_Idx = 1, Image%Number_Of_Elements
                                                                       
            !--- check for bad scans
            if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
               cycle
            endif

            !--- check for space views
            if (Space_Mask(Elem_Idx,Line_Idx) == sym%YES) then
               cycle
            endif

            !--- compute viewing zenith bin for Rtm calculation
            Zen_Idx_Rtm(Elem_Idx,Line_Idx) =  &
              max(1,min(Rtm_Nvzen,ceiling(Geo%Coszen(Elem_Idx,Line_Idx)/Rtm_Vza_Binsize)))

            !--- store cell and angular indices
            Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
            Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
            Lon_Idx_x = I_Nwp_X(Elem_Idx,Line_Idx)
            Lat_Idx_x = J_Nwp_X(Elem_Idx,Line_Idx)
            Lon_x = Lon_Nwp_Fac(Elem_Idx,Line_Idx)
            Lat_x = Lat_Nwp_Fac(Elem_Idx,Line_Idx)
            Zen_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

            !--- if this is the first time this nwp cell is being processed do the following
            if (Rtm(Lon_Idx,Lat_Idx)%Flag == 0) then

               !--- allocate Rtm arrays
               call ALLOCATE_RTM_CELL(Lon_Idx,Lat_Idx,Rtm_Nvzen)

               !--- compute mixing ratio profile
               call COMPUTE_WVMR_PROFILE_NWP(P_Std_Nwp, &
                                     T_Prof_Nwp(:,Lon_Idx,Lat_Idx), &
                                     Rh_Prof_Nwp(:,Lon_Idx,Lat_Idx), &
                                     Wvmr_Nwp)

               !--- compute tpw profiles
               call COMPUTE_TPW_PROFILE_NWP(P_Std_Nwp, &
                                    Wvmr_Nwp,  &
                                    Tpw_Prof_Nwp(:,Lon_Idx,Lat_Idx))

               !--- convert the atmospheric profiles from nwp to Rtm pressure coords
               call CONVERT_ATMOS_PROF_NWP_RTM(NLevels_NWP, &
                                       Sfc_Level_Nwp(Lon_Idx,Lat_Idx), &
                                       Tmpair_Nwp(Lon_Idx,Lat_Idx), &
                                       Rhsfc_Nwp(Lon_Idx,Lat_Idx), &
                                       Psfc_Nwp(Lon_Idx,Lat_Idx), &
                                       P_Std_Nwp, &
                                       T_Prof_Nwp(:,Lon_Idx,Lat_Idx), &
                                       Z_Prof_Nwp(:,Lon_Idx,Lat_Idx), &
                                       Wvmr_Nwp, &
                                       Ozone_Prof_Nwp(:,Lon_Idx,Lat_Idx), &
                                       NLevels_Rtm, &
                                       P_Std_Rtm, &
                                       T_Prof_Rtm, &
                                       Z_Prof_Rtm, &
                                       Wvmr_Prof_Rtm, &
                                       Ozmr_Prof_Rtm, &
                                       T_Std_Rtm, &
                                       Wvmr_Std_Rtm, &
                                       Ozmr_Std_Rtm)

               !--- compute tpw profiles
               call COMPUTE_TPW_PROFILE_NWP(P_Std_Rtm, &
                                    Wvmr_Prof_Rtm,  &
                                    Tpw_Prof_Rtm)

               !--- store in Rtm structures
               Rtm(Lon_Idx,Lat_Idx)%T_Prof = T_Prof_Rtm
               Rtm(Lon_Idx,Lat_Idx)%Z_Prof = Z_Prof_Rtm
               Rtm(Lon_Idx,Lat_Idx)%Wvmr_Prof = Wvmr_Prof_Rtm
               Rtm(Lon_Idx,Lat_Idx)%Ozmr_Prof = Ozmr_Prof_Rtm
               Rtm(Lon_Idx,Lat_Idx)%Tpw_Prof = Tpw_Prof_Rtm

               !--- find key Levels (sfc, tropo, Inversions) for future processing
               call FIND_RTM_LEVELS(Lon_Idx,Lat_Idx)

               !--- set flag so not done again
               Rtm(Lon_Idx,Lat_Idx)%Flag = 1

            end if

            !--- determine if RTM needs to be called for this pixel, if not do the following
            if (Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Flag == 0) then   !check for existing RTM results

               call ALLOCATE_GLOBAL_RTM_STRUCTURE_ELEMENT(Lon_Idx,Lat_Idx,Zen_Idx)

               !--- determine the zenith angle for the RTM
               Satzen_Mid_Bin = acos((Zen_Idx-1)*Rtm_Vza_Binsize + Rtm_Vza_Binsize/2.0) / DTOR

               !--- store structures into RTM arguments
               T_Prof_Rtm = Rtm(Lon_Idx,Lat_Idx)%T_Prof
               Tpw_Prof_Rtm = Rtm(Lon_Idx,Lat_Idx)%Tpw_Prof
               Wvmr_Prof_Rtm = Rtm(Lon_Idx,Lat_Idx)%Wvmr_Prof
               Ozmr_Prof_Rtm = Rtm(Lon_Idx,Lat_Idx)%Ozmr_Prof

               !--------------------------------------------------------------
               ! Call Sensor Specific Fast IR RTM  
               !--------------------------------------------------------------

               !--- decide if this can be skipped
               if ((Lon_Idx_Prev == 0) .or. (abs(Lon_Idx-Lon_Idx_Prev) >= Ilon_Stride) .or. &
                  (Lat_Idx_Prev == 0) .or. (abs(Lat_Idx-Lat_Idx_Prev) >= Ilat_Stride) .or. &
                  (Zen_Idx_Prev == 0) .or. (abs(Zen_Idx - Zen_Idx_Prev) >= Ivza_Stride)) then  

                  Start_Time_Point_Hours_Temp = COMPUTE_TIME_HOURS()

                  !--------------------------------------------------------------
                  !  Call PFAAST Routines for IR channel Transmission Profiles
                  !--------------------------------------------------------------
                  do Chan_Idx = Chan_Idx_Min,Chan_Idx_Max

               !     if (Chan_Idx < 20) cycle
               !     if (Chan_Idx == 26) cycle
               !     if (Chan_Idx == 44) cycle

                     if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle
                     if (ch(Chan_Idx)%Obs_Type == SOLAR_OBS_TYPE) cycle
                     if (ch(Chan_Idx)%Obs_Type == LUNAR_OBS_TYPE) cycle
                     
                
                     Chan_Idx_For_Pfaast = Chan_Idx
                     if (Chan_Idx .eq. 45) Chan_Idx_For_Pfaast = 33
                     Sc_Name_Rtm = sensor_name_for_rtm(sensor%wmo_id,sensor%sensor_name, chan_idx)
                     
                     CALL COMPUTE_TRANSMISSION_PFAAST( &
                           trim(Ancil_Data_Dir) &
                        ,  T_Prof_rtm &
                        ,  Wvmr_Prof_Rtm &
                        ,  Ozmr_Prof_Rtm &
                        ,  Satzen_Mid_Bin &
                        ,  CO2_RATIO &
                        ,  Sc_Name_Rtm &
                        ,  Chan_Idx_For_Pfaast &
                        ,  Trans_Prof_Rtm &
                        ,  Use_Modis_Channel_Equivalent = .true.  ) 

                     !---- Copy the output to appropriate channel's tranmission vector
                     Trans_Atm_Prof(:,Chan_Idx) = Trans_Prof_Rtm

                  end do

                  !--------------------------------------------------------------
                  ! compute transmission profiles for solar channels
                  !--------------------------------------------------------------
                  do Chan_Idx = Chan_Idx_Min,Chan_Idx_Max

                     if (Chan_Idx >= 20 .and. Chan_Idx /= 26) cycle

                     call SOLAR_TRANS(Chan_Idx,Satzen_Mid_Bin,Error_Status) 

                     !---- Copy the output to appropriate channel's tranmission vector
                     Trans_Atm_Prof(:,Chan_Idx) = Trans_Prof_Rtm

                  end do
     
                  End_Time_Point_Hours_Temp = COMPUTE_TIME_HOURS()

                  Segment_Time_Point_Seconds_Temp =  Segment_Time_Point_Seconds_Temp + &
                    60.0*60.0*(End_Time_Point_Hours_Temp - Start_Time_Point_Hours_Temp)
    
                  !----------------------------------------------------------------------
                  ! Manipulate Transmission Profiles
                  !----------------------------------------------------------------------

                  !----------------------------------------------------------------------
                  ! Apply Gamma Scaling
                  !----------------------------------------------------------------------
                  do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max
                     if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES .and. Gamma_Trans_Factor(Chan_Idx)/=1.0) then
                        Trans_Atm_Prof(:,Chan_Idx)  = Trans_Atm_Prof(:,Chan_Idx) ** Gamma_Trans_Factor(Chan_Idx)
                     end if
                  end do

                  !----------------------------------------------------------------------
                  ! Compute Tranmissions along solar and total (two-way) paths
                  !----------------------------------------------------------------------
                  do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max
                     if (Chan_Idx > 21 .and. Chan_Idx /= 26) cycle

                     Trans_Atm_Total_Prof(:,Chan_Idx) = Trans_Atm_Prof(:,Chan_Idx) **  &
                                           (1.0/Geo%Cossolzen(Elem_Idx,Line_Idx))

                     Trans_Atm_Solar_Prof(:,Chan_Idx) = Trans_Atm_Prof(:,Chan_Idx) **  &
                                           (Geo%Coszen(Elem_Idx,Line_Idx)/Geo%Cossolzen(Elem_Idx,Line_Idx))

                  end do

                  Lon_Idx_Prev = Lon_Idx
                  Lat_Idx_Prev = Lat_Idx
                  Zen_Idx_Prev = Zen_Idx
                  Trans_Atm_Prof_Prev = Trans_Atm_Prof
                  Trans_Atm_Solar_Prof_Prev = Trans_Atm_Solar_Prof
                  Trans_Atm_Total_Prof_Prev = Trans_Atm_Total_Prof

               else

                  print *, "Skipped an RTM Call"
                  Trans_Atm_Prof = Trans_Atm_Prof_Prev
                  Trans_Atm_Solar_Prof = Trans_Atm_Solar_Prof_Prev
                  Trans_Atm_Total_Prof = Trans_Atm_Total_Prof_Prev

               end if    !skip rtm if statement

               !--- compute profiles of radiance (atm and bb cloud)
               call COMPUTE_CLEAR_RAD_PROFILES_RTM()

               !--- copy local rtm profiles back into global rtm structure
               call COPY_LOCAL_RTM_TO_GLOBAL_RTM_STRUCTURE(Lon_Idx,Lat_Idx,Zen_Idx)

               !---   set mask  to indicate this bin or this cell has been computed
               Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Flag = 1

            end if  !end check for previous RTM run

            !-----------------------------------------------------------------------------
            ! compute clear TOA radiance
            !-----------------------------------------------------------------------------  

            !---- compute the surface level for this pixel based on high-res elevation
            Sfc_Level_Rtm_Pixel(Elem_Idx,Line_Idx) = Rtm(Lon_Idx,Lat_Idx)%Sfc_Level

            if ((Sfc%Land(Elem_Idx,Line_Idx) == sym%LAND) .and. (Sfc%Zsfc(Elem_Idx,Line_Idx) /= Missing_Value_Real4)) then
               call LOCATE(Rtm(Lon_Idx,Lat_Idx)%Z_Prof,NLevels_Rtm,Sfc%Zsfc(Elem_Idx,Line_Idx),Sfc_Level_Rtm_Pixel(Elem_Idx,Line_Idx))
            endif
            Sfc_Level_Rtm_Pixel(Elem_Idx,Line_Idx) = max(1,min(NLevels_Rtm-1,Sfc_Level_Rtm_Pixel(Elem_Idx,Line_Idx)))

            !--- use this scalar for visual clarity
            Sfc_Level_Idx = Sfc_Level_Rtm_Pixel(Elem_Idx,Line_Idx)

            !--- if an sst analysis is available, use that
            if ((Sfc%Land_Mask(Elem_Idx,Line_Idx) == sym%NO) .and.  &
              (Sfc%Snow(Elem_Idx,Line_Idx) == sym%NO_SNOW) .and.  &
              (Use_Sst_Anal == sym%YES) .and.  &
              (Sst_Anal(Elem_Idx,Line_Idx) > 270.0 )) then
                 Tsfc_Nwp_Pix(Elem_Idx,Line_Idx) = Sst_Anal(Elem_Idx,Line_Idx)
            end if 

            !-- vertical interp weight
            Prof_Weight = (Sfc%Zsfc(Elem_Idx,Line_Idx) - Rtm(Lon_Idx,Lat_Idx)%Z_Prof(Sfc_Level_Idx)) / &
              (Rtm(Lon_Idx,Lat_Idx)%Z_Prof(Sfc_Level_Idx+1) - Rtm(Lon_Idx,Lat_Idx)%Z_Prof(Sfc_Level_Idx))

            !--- constrain - important when high res topo differs from low res nwp topo
            Prof_Weight = max(0.0,min(1.0,Prof_Weight))

            !--- call routine to compute radiative transfer terms such as Rad_Atm
            !--- Trans_Atm, and clear-sky radinace and brightness temperature
            !--- map global to local to allow for efficient looping
            call COMPUTE_CHANNEL_RT(Sfc_Level_Idx,Prof_Weight,Lon_Idx,Lat_Idx,Elem_Idx,Line_Idx,Zen_Idx)

            !--- compute Ch20 Emissivities
            call COMPUTE_CH20_EMISSIVITY(Elem_Idx,Line_Idx)

            !--- compute Emissivity at tropopause
            call COMPUTE_TROPOPAUSE_EMISSIVITIES(Elem_Idx,Line_Idx,Lon_Idx,Lat_Idx,Zen_Idx)

            !--- compute split-window beta ratio at tropopause
            call COMPUTE_BETA_RATIOES(Elem_Idx,Line_Idx)

         end do element_loop
      end do line_loop

      return

   end subroutine GET_PIXEL_NWP_RTM


   !====================================================================
   ! subroutine Name: FIND_RTM_LEVELS
   !
   ! Function:
   ! Finds various key Levels in the RTM (tropopause, etc), for each RTM gridcell
   !
   !====================================================================
   subroutine FIND_RTM_LEVELS(Lon_Idx,Lat_Idx)

      integer (kind=int4), intent(in):: Lon_Idx, Lat_Idx
      integer (kind=int4):: k

      !--- find surface Level - this is closest nwp Level above the actual surface
      Rtm(Lon_Idx,Lat_Idx)%Sfc_Level = NLevels_Rtm
      do k = NLevels_Rtm,1,-1
         if (P_Std_Rtm(k) < Psfc_Nwp(Lon_Idx,Lat_Idx)) then
            Rtm(Lon_Idx,Lat_Idx)%Sfc_Level = k
            exit
         end if
      end do

      !--------------------------------------------------------------------
      !--- find tropopause Level  based on tropopause pressure
      !--- tropopause is between tropopause_Level and tropopaue_Level + 1
      !--------------------------------------------------------------------
      do k = 1, Rtm(Lon_Idx,Lat_Idx)%Sfc_Level-1
         if ((P_Std_Rtm(k) <= P_trop_Nwp(Lon_Idx,Lat_Idx)) .and. &
            (P_Std_Rtm(k+1) > P_trop_Nwp(Lon_Idx,Lat_Idx))) then
            Rtm(Lon_Idx,Lat_Idx)%Tropo_Level = k
         end if
      end do

      !--- check if tropopause Level found
      if (Rtm(Lon_Idx,Lat_Idx)%Tropo_Level == 0) then
         print *, EXE_PROMPT, "Error, tropopause Level not found"
      end if         

      do k = 1, Rtm(Lon_Idx,Lat_Idx)%Sfc_Level-1
         if ((P_Std_Rtm(k) <= 850.0) .and. &
               (P_Std_Rtm(k+1) > 850.0)) then
               
            Rtm(Lon_Idx,Lat_Idx)%Level850 = k
         endif
      enddo
      
      do k = 1, Rtm(Lon_Idx,Lat_Idx)%Sfc_Level-1
         if ((P_Std_Rtm(k) <= 440.0) .and. &
            (P_Std_Rtm(k+1) > 440.0)) then
            
            Rtm(Lon_Idx,Lat_Idx)%Level440 = k
         endif
      enddo

      !---------------------------------------------------------------------
      ! find Inversion Level - highest Level Inversion below trop
      !---------------------------------------------------------------------
      Rtm(Lon_Idx,Lat_Idx)%Inversion_Level = 0
      if (Rtm(Lon_Idx,Lat_Idx)%Tropo_Level > 0 .and. Rtm(Lon_Idx,Lat_Idx)%Sfc_Level > 0) then
         do k = Rtm(Lon_Idx,Lat_Idx)%Tropo_Level,Rtm(Lon_Idx,Lat_Idx)%Sfc_Level-1
            if ((Rtm(Lon_Idx,Lat_Idx)%T_Prof(k) - Rtm(Lon_Idx,Lat_Idx)%T_Prof(k+1) > delta_t_Inversion) .and. &
                        (P_Std_Rtm(k) >= p_Inversion_min)) then
               Rtm(Lon_Idx,Lat_Idx)%Inversion_Level = k
               exit 
            endif
         enddo
      endif

   end subroutine FIND_RTM_LEVELS

   !====================================================================
   ! subroutine Name: COMPUTE_WVMR_PROFILE_NWP
   !
   ! Function:
   ! Computes the NWP WVMR profile
   !
   ! Input: Lon_Idx = Longitude Index for NWP cell
   !        Lat_Idx = Latitude Index for NWP cell
   !        Press_Profile = pressure profile (hPa)
   !        Temp_Profile = temperature profile (K)
   !        Rh_Profile = relative humidity profile (%)
   !
   ! Output: Wvmr_Profile = water vapor mixing ratio profile (g/kg)
   !         Tpw Profile = profile of Tpw from level to space (g/m^2)
   !
   !====================================================================
   subroutine COMPUTE_WVMR_PROFILE_NWP(Press_Profile, &
                                    Temp_Profile, &
                                    Rh_Profile, &
                                    Wvmr_Profile)

      real, intent(in), dimension(:):: Temp_Profile
      real, intent(in), dimension(:):: Press_Profile
      real, intent(in), dimension(:):: Rh_Profile
      real, intent(out), dimension(:):: Wvmr_Profile
      integer:: Lev_Idx
      real:: e
      real:: es

      integer:: Nlevels

      Nlevels = size(Press_Profile)

      !--- make Wvmr_Profile for use in RTM
      do Lev_Idx = 1, NLevels
         es = VAPOR(Temp_Profile(Lev_Idx))
         e = Rh_Profile(Lev_Idx) * es / 100.0
         Wvmr_Profile(Lev_Idx) = 1000.0*0.622 * (e / (Press_Profile(Lev_Idx) - e))  !(g/kg)
      end do

   end subroutine COMPUTE_WVMR_PROFILE_NWP

   !====================================================================
   ! subroutine Name: COMPUTE_TPW_PROFILE_NWP
   !
   ! Function:
   ! Computes the NWP TPW profile
   !
   ! Input: Lon_Idx = Longitude Index for NWP cell
   !        Lat_Idx = Latitude Index for NWP cell
   !        Press_Profile = pressure profile (hPa)
   !        Temp_Profile = temperature profile (K)
   !        Rh_Profile = relative humidity profile (%)
   !
   ! Output: Wvmr_Profile = water vapor mixing ratio profile (g/kg)
   !         Tpw Profile = profile of Tpw from level to space (g/m^2)
   !
   !====================================================================
   
   subroutine COMPUTE_TPW_PROFILE_NWP(Press_Profile, &
                                   Wvmr_Profile,  &
                                   Tpw_Profile)
   
      real, intent(in), dimension(:):: Press_Profile
      real, intent(in), dimension(:):: Wvmr_Profile
      real, intent(out), dimension(:):: Tpw_Profile
      
      integer:: Lay_Idx
      real :: w_mean
      real :: u_layer
      integer:: Nlevels

      Nlevels = size(Press_Profile)

      !--- make tpw profile for use in atmospheric correction
      Tpw_Profile(1) = 0.0
      do Lay_Idx = 1, NLevels-1   !layer index
         w_mean = 0.5*(Wvmr_Profile(Lay_Idx+1)+Wvmr_Profile(Lay_Idx)) / 1000.0  !(kg/kg)
         u_layer = (10.0/g)*(Press_Profile(Lay_Idx+1)-Press_Profile(Lay_Idx))*w_mean
         Tpw_Profile(Lay_Idx+1) = Tpw_Profile(Lay_Idx) + u_layer
      end do

   end subroutine COMPUTE_TPW_PROFILE_NWP

   !====================================================================
   !   subroutine Name: ALLOCATE_RTM
   !
   ! Function:
   ! Subroutine that allocates memory for nlon x nlat Rtm structure.
   !
   !====================================================================
   subroutine ALLOCATE_RTM(nx,ny)
      integer (kind=int4), intent(in) :: nx, ny

      integer :: Alloc_Status

      allocate(Rtm(nx,ny),stat=Alloc_Status)
      if (Alloc_Status /= 0) then
         print "(a,'Not enough memory to allocate Rtm_Params structure.')",EXE_PROMPT
         stop
      end if

   end subroutine ALLOCATE_RTM

   !====================================================================
   ! subroutine Name: DEALLOCATE_RTM
   !
   ! Function:
   !  Subroutine that deallocates memory used for nlon x nlat Rtm structure.
   !
   !====================================================================
   subroutine DEALLOCATE_RTM()

      integer :: Alloc_Status

      deallocate(Rtm,stat=Alloc_Status)
      if (Alloc_Status /= 0) then
         print "(a,'Error deallocating Rtm_Params structure.')",EXE_PROMPT
         stop
      end if

   end subroutine DEALLOCATE_RTM

   !====================================================================
   ! subroutine Name: ALLOCATE_RTM_CELL
   !
   ! Function:
   !  Subroutine to allocate memory for the RTM structure.
   !
   !====================================================================
   subroutine ALLOCATE_RTM_CELL(Lon_Idx,Lat_Idx,Nvza)
      integer (kind=int4), intent(in) :: Lon_Idx, Lat_Idx
      integer (kind=int4), intent(in) :: Nvza

      integer :: Alloc_Status

      allocate(Rtm(Lon_Idx, Lat_Idx)%d(Nvza),stat=Alloc_Status)
      allocate(Rtm(Lon_Idx, Lat_Idx)%T_Prof(NLevels_Rtm),stat=Alloc_Status)
      allocate(Rtm(Lon_Idx, Lat_Idx)%Z_Prof(NLevels_Rtm),stat=Alloc_Status)
      allocate(Rtm(Lon_Idx, Lat_Idx)%Wvmr_Prof(NLevels_Rtm),stat=Alloc_Status)
      allocate(Rtm(Lon_Idx, Lat_Idx)%Ozmr_Prof(NLevels_Rtm),stat=Alloc_Status)
      allocate(Rtm(Lon_Idx, Lat_Idx)%Tpw_Prof(NLevels_Rtm),stat=Alloc_Status)

      if (Alloc_Status /= 0) then
         print "(a,'Not enough memory to allocate Rtm_Params structure.')",EXE_PROMPT
         stop
      endif

      Rtm(Lon_Idx,Lat_Idx)%Flag = 1

   end subroutine ALLOCATE_RTM_CELL

   !====================================================================
   ! subroutine Name: DEALLOCATE_RTM_CELL
   !
   ! Function:
   !  Subroutine to deallocates memory for the RTM structure.
   !
   !====================================================================
   subroutine DEALLOCATE_RTM_CELL(Lon_Idx,Lat_Idx)
      integer (kind=int4), intent(in) :: Lon_Idx
      integer (kind=int4), intent(in) :: Lat_idx
      integer :: Alloc_Status
      
      deallocate(Rtm(Lon_Idx,Lat_Idx)%d,stat=Alloc_Status)
      deallocate(Rtm(Lon_Idx,Lat_Idx)%T_Prof,stat=Alloc_Status)
      deallocate(Rtm(Lon_Idx,Lat_Idx)%Z_Prof,stat=Alloc_Status)      
      deallocate(Rtm(Lon_Idx,Lat_Idx)%Wvmr_Prof,stat=Alloc_Status)
      deallocate(Rtm(Lon_Idx,Lat_Idx)%Ozmr_Prof,stat=Alloc_Status)
      deallocate(Rtm(Lon_Idx,Lat_Idx)%Tpw_Prof,stat=Alloc_Status)
      
      if (Alloc_Status /= 0) then
         print "(a,'Error deallocating Rtm_Params structure.')",EXE_PROMPT
         stop
      endif
  
      Rtm(Lon_Idx,Lat_Idx) % Flag = 0
      Rtm(Lon_Idx,Lat_Idx) % Sfc_Level = 0
      Rtm(Lon_Idx,Lat_Idx) % Inversion_Level = 0
      Rtm(Lon_Idx,Lat_Idx) % Tropo_Level = 0
      Rtm(Lon_Idx,Lat_Idx) % Level440 = 0
      Rtm(Lon_Idx,Lat_Idx) % Level850 = 0

   end subroutine DEALLOCATE_RTM_CELL

   !====================================================================
   ! subroutine Name: DEALLOCATE_RTM_VARS
   !
   ! Function:
   !  Subroutine to deallocate memory for the RTM profile arrays.
   !
   !====================================================================
   subroutine DEALLOCATE_RTM_VARS(Lon_Idx,Lat_Idx,Zen_Idx)
      integer (kind=int4), intent(in) :: Lon_Idx
      integer (kind=int4), intent(in) :: Lat_Idx
      integer (kind=int4), intent(in) :: Zen_Idx

      integer :: Chan_Idx
      integer :: Alloc_Status

      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max

         if (allocated(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Profile)) &
            deallocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Profile, stat=Alloc_Status)

         if (allocated(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Solar_Profile)) &
            deallocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Solar_Profile, stat=Alloc_Status)

         if (allocated(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Total_Profile)) &
            deallocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Total_Profile, stat=Alloc_Status)

         if (allocated(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Profile)) &
            deallocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Profile, stat=Alloc_Status)

         if (allocated(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Dwn_Profile)) &
            deallocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Dwn_Profile, stat=Alloc_Status)

         if (allocated(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_BB_Cloud_Profile)) &
            deallocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_BB_Cloud_Profile, stat=Alloc_Status)

         if (Alloc_Status /= 0) then
            print "(a,'Error destroying allocate channel Rtm profile.')",EXE_PROMPT
            stop
         endif

      enddo

      !--- set flag to zero (so we don't deallocate again)
      Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Flag = 0

   end subroutine DEALLOCATE_RTM_VARS

   !--------------------------------------------------------------
   ! Compute Gamma Factor for Radiance Bias Adjustment
   !
   ! based on satellite number, channel number and nwp source,
   ! set the value of Gamma_Trans_Factor to reduce RTM bias
   !
   ! nwp flag  (1=gfs,2=ncep reanalysis,3=cfsr)
   !--------------------------------------------------------------
   subroutine COMPUTE_GAMMA_FACTOR()

      !--- initialize to unity
      Gamma_Trans_Factor = 1.0

      !--- GOES-10
      if (Sensor%WMO_Id == 254) then
         if (Nwp_Opt == 3) then
            Gamma_Trans_Factor(20) = 1.25
            Gamma_Trans_Factor(27) = 0.79
            Gamma_Trans_Factor(31) = 1.15
            Gamma_Trans_Factor(31) = 1.05
         endif
       endif
       
      !--- GOES-11
      if (Sensor%WMO_Id == 255) then
         if (Nwp_opt == 3) then
            Gamma_Trans_Factor(20) = 1.35
            Gamma_Trans_Factor(27) = 0.74
            Gamma_Trans_Factor(31) = 1.05
            Gamma_Trans_Factor(32) = 1.05
         endif
      endif
      
      !--- GOES-12
      if (Sensor%WMO_Id == 256) then
         if (nwp_opt == 1 .or. nwp_opt == 3) then    !repeat of cfsr
            Gamma_Trans_Factor(20) = 1.45 
            Gamma_Trans_Factor(27) = 0.79 
            Gamma_Trans_Factor(31) = 1.15
            Gamma_Trans_Factor(33) = 1.15
         end if       
      end if
      
      !--- GOES-13
      if (Sensor%WMO_Id == 257) then
         if (nwp_opt == 1 .or. nwp_opt == 3) then
            Gamma_Trans_Factor(20) = 1.45
            Gamma_Trans_Factor(27) = 0.83
            Gamma_Trans_Factor(31) = 1.05
            Gamma_Trans_Factor(33) = 1.19
         end if      
      end if

      !--- GOES-15
      if (Sensor%WMO_Id == 259) then
         Gamma_Trans_Factor(20) = 1.55
         Gamma_Trans_Factor(27) = 0.728
         Gamma_Trans_Factor(31) = 1.05
         Gamma_Trans_Factor(33) = 1.075
      end if
    
      !--- MET-09
      if (Sensor%WMO_Id == 56) then
     
         Gamma_Trans_Factor(20) = 1.55
         Gamma_Trans_Factor(27) = 0.96
         Gamma_Trans_Factor(29) = 1.35
         Gamma_Trans_Factor(31) = 0.95
         Gamma_Trans_Factor(32) = 0.95
         Gamma_Trans_Factor(33) = 1.11
     
      endif
      
      !--- MTSAT-02
      if (Sensor%WMO_Id == 172) then
         if (nwp_opt == 1) then
            Gamma_Trans_Factor(20) = 1.25
            Gamma_Trans_Factor(27) = 0.87
            Gamma_Trans_Factor(31) = 1.15
            Gamma_Trans_Factor(32) = 1.15
         end if
      end if
      
      !--- COMS-1
      if (Sensor%WMO_Id == 810) then
         Gamma_Trans_Factor(20) = 1.25
         Gamma_Trans_Factor(27) = 0.936
         Gamma_Trans_Factor(31) = 1.05
         Gamma_Trans_Factor(32) = 1.075
      end if

   end subroutine COMPUTE_GAMMA_FACTOR


   !--------------------------------------------------------------------------------------------------
   !> subroutine NAME: sensor_name_for_rtm
   !!
   !! Description:
   !! Knowing the WMO Satellite Identification Number
   !!
   !--------------------------------------------------------------------------------------------------
   
   function sensor_name_for_rtm ( wmo_id, sensorname, chan_idx ) result ( sensor_name_rtm)
      
      integer, intent(in) :: wmo_id
      character (len =*) , intent(in) :: sensorname
      integer, intent(in) :: chan_idx
      character (len =20 ) ::  sensor_name_rtm
      integer :: i
    
      select case(WMO_Id)

      case(4) !METOP-A
         sensor_name_rtm = 'AVHRR-METOPA'

      case(3) !METOP-B
         sensor_name_rtm = 'AVHRR-METOPB'

      case(5) !METOP-C
        sensor_name_rtm = 'AVHRR-METOPC'

      case(55) !MSG-8
         sensor_name_rtm = 'SEVIRI-MSG08'
      
      case(56) !MSG-9
         sensor_name_rtm = 'SEVIRI-MSG09'
      
      case(57) !MSG-10
         sensor_name_rtm = 'SEVIRI-MSG10'
      
      case(171) !MTSAT-1R
         sensor_name_rtm = 'MTSAT-1'
      
      case(172) !MTSAT-2
         sensor_name_rtm = 'MTSAT-2'
         
      case(173) !AHI
         sensor_name_rtm = 'AHI'
               
      case(200) !NOAA-8
        sensor_name_rtm = 'AVHRR-NOAA08'

      case(201) !NOAA-9
        sensor_name_rtm = 'AVHRR-NOAA09'

      case(202) !NOAA-10
        sensor_name_rtm = 'AVHRR-NOAA10'

      case(203) !NOAA-11
        sensor_name_rtm = 'AVHRR-NOAA11'

      case(204) !NOAA-12
        sensor_name_rtm = 'AVHRR-NOAA12'

      case(205) !NOAA-14
        sensor_name_rtm = 'AVHRR-NOAA14'

      case(206) !NOAA-15
        sensor_name_rtm = 'AVHRR-NOAA15'

      case(207) !NOAA-16
        sensor_name_rtm = 'AVHRR-NOAA16'

      case(208) !NOAA-17
        sensor_name_rtm = 'AVHRR-NOAA17'

      case(209) !NOAA-18
        sensor_name_rtm = 'AVHRR-NOAA18'

      case(223) !NOAA-19
        sensor_name_rtm = 'AVHRR-NOAA19'

      case(224) !VIIRS - 
        sensor_name_rtm = 'VIIRS'

      case(252) !GOES-8
        sensor_name_rtm = 'GOES-8'

      case(253) !GOES-9
        sensor_name_rtm = 'GOES-9'

      case(254) !GOES-10
        sensor_name_rtm = 'GOES-10'

      case(255) !GOES-11
        sensor_name_rtm = 'GOES-11'

      case(256) !GOES-12
        sensor_name_rtm = 'GOES-12'

      case(257) !GOES-13
        sensor_name_rtm = 'GOES-13'

      case(258) !GOES-14
        sensor_name_rtm = 'GOES-14'

      case(259) !GOES-15
        sensor_name_rtm = 'GOES-15'

      case(706) !NOAA-6
        sensor_name_rtm = 'AVHRR-NOAA06'

      case(707) !NOAA-7
        sensor_name_rtm = 'AVHRR-NOAA07'

      case(708) !NOAA-5
        sensor_name_rtm = 'AVHRR-TIROSN'

      case(783) !MODIS 
          sensor_name_rtm = 'MODIS-TERRA'

      case(784) !MODIS 
         sensor_name_rtm = 'MODIS-AQUA'

      case(810) !COMS
         sensor_name_rtm ='FY2-1'
         
      case(514) !FY2D      
         sensor_name_rtm ='FY2-2'
         
      case(515) !FY2E          
         sensor_name_rtm ='FY2-3'
      
      case default
         print*,'sensor for WMO number not found in RT Utils  ', WMO_id  
         print*,'stopping ... Please fix this in rt_utils.F90'
         print*,' better tell andi.walther@ssec.wisc.edu'
         stop    
      end select
     
      
      
      if (trim ( Sensorname) == 'AVHRR-IFF') then
        
         !  sensor for channels 21:30 and 33:36 is HIRS
         if ( any ( chan_idx ==  [ (i,i=21,30,1) , 33,34,35,36] ) ) then
            ! - for this IFF sensor_name_rtm is initially set to AVHRR-<Satellite>
            select case(WMO_Id)
            
            case(4) !METOP-A
               sensor_name_rtm = 'HIRS-METOPA'

            case(3) !METOP-B
               sensor_name_rtm = 'HIRS-METOPB'

            case(5) !METOP-C
               sensor_name_rtm = 'HIRS-METOPC'
            
            case(200) !NOAA-8
               sensor_name_rtm = 'HIRS-NOAA08'

            case(201) !NOAA-9
               sensor_name_rtm = 'HIRS-NOAA09'

            case(202) !NOAA-10
               sensor_name_rtm = 'HIRS-NOAA10'

            case(203) !NOAA-11
               sensor_name_rtm = 'HIRS-NOAA11'

            case(204) !NOAA-12
               sensor_name_rtm = 'HIRS-NOAA12'

            case(205) !NOAA-14
               sensor_name_rtm = 'HIRS-NOAA14'

            case(206) !NOAA-15
               sensor_name_rtm = 'HIRS-NOAA15'

            case(207) !NOAA-16
               sensor_name_rtm = 'HIRS-NOAA16'

            case(208) !NOAA-17
               sensor_name_rtm = 'HIRS-NOAA17'

            case(209) !NOAA-18
               sensor_name_rtm = 'HIRS-NOAA18'

            case(223) !NOAA-19
               sensor_name_rtm = 'HIRS-NOAA19'
            
            case(706) !NOAA-6
               sensor_name_rtm = 'HIRS-NOAA06'

            case(707) !NOAA-7
               sensor_name_rtm = 'HIRS-NOAA07'

            case(708) !NOAA-5
               sensor_name_rtm = 'HIRS-TIROSN'
            
            case default
               print*,'sensor for WMO number not found in RT Utils for AVHRR-IFF  ', WMO_id  
               print*,'stopping ... Please fix this in rt_utils.F90'
               print*,' better tell andi.walther@ssec.wisc.edu'
               stop    
            end select
            
           
         end if   
      end if
   
      if (trim ( Sensorname) == 'VIIRS-IFF') then
         
         !  sensor for channels 27:28 and 33:36 is CRISP this is similar to MODIS-AQUA
         if ( any ( chan_idx ==  [27,28, 33,34,35,36,45] ) ) sensor_name_rtm   = 'MODIS-AQUA'
         
      end if
   
   
   end function sensor_name_for_rtm 
   

   !--------------------------------------------------------------------------------------------------
   !  set up values for the solar rtm calculations for this sensor
   !--------------------------------------------------------------------------------------------------
   subroutine SOLAR_TRANS(Chan_Idx,Zen_Ang,Error_Status)

      integer, intent(in):: Chan_Idx
      real, intent(in):: Zen_Ang
      integer, intent(out):: Error_Status

      real, dimension(3):: Tau_H2O_Coef
      real:: Tau_O2_Column
      real:: Tau_CO2_Column
      real:: Tau_CH4_Column
      real:: Tau_O3_Column

      real:: Tau_O2
      real:: Tau_CO2
      real:: Tau_CH4
      real:: Tau_O3
      real:: Tau_H2O

      real:: Tpw
      real:: mu
      real:: Tau_Gas

      integer:: Lev_Idx

      Error_Status = 1

      if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) return
!     if (Chan_Idx >= 20 .and. Chan_Idx /= 26 .and. Chan_Idx/= 44) return
      if (ch(Chan_Idx)%Obs_Type /= SOLAR_OBS_TYPE .or. ch(Chan_Idx)%Obs_Type /= LUNAR_OBS_TYPE) return

      Trans_Prof_Rtm = 1.0

      mu = cos(Zen_Ang*DTOR)
 
      !--- initialize
      Tau_H2O_Coef = Solar_Rtm%Tau_H2O_Coef(Chan_Idx,:)
      Tau_O2_Column = Solar_Rtm%Tau_O2(Chan_Idx)
      Tau_O3_Column = Solar_Rtm%Tau_O3(Chan_Idx)
      Tau_CO2_Column = Solar_Rtm%Tau_CO2(Chan_Idx)
      Tau_CH4_Column = Solar_Rtm%Tau_CH4(Chan_Idx)

      !--- loop through layers and fill profile
      do Lev_Idx = 1, Nlevels_Rtm
         Tpw = Tpw_Prof_Rtm(Lev_Idx)
         Tau_H2O = Tau_H2O_Coef(1) + Tau_H2O_Coef(2)*Tpw_Prof_Rtm(Lev_Idx) + Tau_H2O_Coef(3)*Tpw_Prof_Rtm(Lev_Idx)**2
         Tau_CO2 = Tau_Co2_Column * P_Std_Rtm(Lev_Idx) / P_Std_Rtm(Nlevels_Rtm)
         Tau_Ch4 = Tau_Ch4_Column * P_Std_Rtm(Lev_Idx) / P_Std_Rtm(Nlevels_Rtm)
         Tau_O2 = Tau_O2_Column * P_Std_Rtm(Lev_Idx) / P_Std_Rtm(Nlevels_Rtm)
         Tau_O3 = Tau_O3_Column * P_Std_Rtm(Lev_Idx) / P_Std_Rtm(Nlevels_Rtm)
         Tau_Gas = max(0.0,Tau_H2O + Tau_O3 + Tau_O2 + Tau_CO2 + Tau_CH4)
         Trans_Prof_Rtm(Lev_Idx) = exp(-Tau_Gas / mu)
      end do

      Error_Status = 0

   end subroutine SOLAR_TRANS
   
   !==============================================================================
   !
   !==============================================================================
   subroutine SETUP_SOLAR_RTM ( WMO_Id )

      integer, intent(in):: WMO_Id

      !------------------------------------------------------------------------
      ! create other arrays necessary for pixel processing
      ! note, n05-n14 are fill values as there is no ch3a
      !------------------------------------------------------------------------
      Solar_Rtm%Tau_H2O_Coef(1,:) = 0.0
      Solar_Rtm%Tau_Ray(1) = 0.0568  !ch1
      Solar_Rtm%Tau_O3(1) = 0.0342
      Solar_Rtm%Tau_O2(1) = 0.0007
      Solar_Rtm%Tau_CH4(1) = 0.00
      Solar_Rtm%Tau_CO2(1) = 0.00

      Solar_Rtm%Tau_H2O_Coef(2,:) = 0.0
      Solar_Rtm%Tau_Ray(2) = 0.0187  !ch2
      Solar_Rtm%Tau_O3(2) = 0.0008
      Solar_Rtm%Tau_O2(2) = 0.0137
      Solar_Rtm%Tau_CO2(2) = 0.00
      Solar_Rtm%Tau_CH4(2) = 0.00

      Solar_Rtm%Tau_H2O_Coef(5,:) = 0.0
      Solar_Rtm%Tau_Ray(5) = 0.00
      Solar_Rtm%Tau_O3(5) = 0.00
      Solar_Rtm%Tau_O2(5) = 0.00
      Solar_Rtm%Tau_CO2(5) = 0.0
      Solar_Rtm%Tau_CH4(5) = 0.0

      Solar_Rtm%Tau_H2O_Coef(6,:) = 0.0
      Solar_Rtm%Tau_Ray(6) = 0.0013
      Solar_Rtm%Tau_O3(6) = 0.00
      Solar_Rtm%Tau_O2(6) = 0.00
      Solar_Rtm%Tau_CO2(6) = 0.0161
      Solar_Rtm%Tau_CH4(6) = 0.0004

      Solar_Rtm%Tau_H2O_Coef(7,:) = 0.0
      Solar_Rtm%Tau_Ray(7) = 0.00
      Solar_Rtm%Tau_O3(7) = 0.00
      Solar_Rtm%Tau_O2(7) = 0.00
      Solar_Rtm%Tau_CO2(7) = 0.0
      Solar_Rtm%Tau_CH4(7) = 0.0

      Solar_Rtm%Tau_H2O_Coef(20,:) = 0.0
      Solar_Rtm%Tau_Ray(20) = 0.00
      Solar_Rtm%Tau_O3(20) = 0.00
      Solar_Rtm%Tau_O2(20) = 0.00
      Solar_Rtm%Tau_CO2(20) = 0.00
      Solar_Rtm%Tau_CH4(20) = 0.000

      Solar_Rtm%Tau_H2O_Coef(44,:) = 0.0
      Solar_Rtm%Tau_Ray(44) = 0.00
      Solar_Rtm%Tau_O3(44) = 0.00
      Solar_Rtm%Tau_O2(44) = 0.00
      Solar_Rtm%Tau_CO2(44) = 0.00
      Solar_Rtm%Tau_CH4(44) = 0.000

      !---- define channel aerosol properties
      Solar_Rtm%Tau_Aer(1) = 0.12
      Solar_Rtm%Wo_Aer(1) = 0.8
      Solar_Rtm%G_Aer(1) = 0.6
      Solar_Rtm%Tau_Aer(2) = 0.03
      Solar_Rtm%Wo_Aer(2) = 0.8
      Solar_Rtm%G_Aer(2) = 0.6
      Solar_Rtm%Tau_Aer(5) = 0.00
      Solar_Rtm%Wo_Aer(5) = 0.8
      Solar_Rtm%G_Aer(5) = 0.6
      Solar_Rtm%Tau_Aer(6) = 0.00
      Solar_Rtm%Wo_Aer(6) = 0.8
      Solar_Rtm%G_Aer(6) = 0.6
      Solar_Rtm%Tau_Aer(7) = 0.00
      Solar_Rtm%Wo_Aer(7) = 0.8
      Solar_Rtm%G_Aer(7) = 0.6
      Solar_Rtm%Tau_Aer(20) = 0.00
      Solar_Rtm%Wo_Aer(20) = 0.8
      Solar_Rtm%G_Aer(20) = 0.6
      Solar_Rtm%Tau_Aer(44) = 0.00
      Solar_Rtm%Wo_Aer(44) = 0.8
      Solar_Rtm%G_Aer(44) = 0.6


      select case(WMO_Id)

      case(3)      !Metop-A
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00005603, 0.00396121,-0.00012323/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02217304, 0.04638705,-0.00335894/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00003428, 0.00106910,-0.00001267/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00640370, 0.05971186,-0.00324950/)
        Solar_Rtm%Tau_Ray(1) = 0.05561
        Solar_Rtm%Tau_O3(1) =  0.02444
        Solar_Rtm%Tau_O2(1) =  0.00116
        Solar_Rtm%Tau_CO2(1) = 0.00000
        Solar_Rtm%Tau_CH4(1) = 0.00000
        Solar_Rtm%Tau_Ray(2) = 0.01769
        Solar_Rtm%Tau_O3(2) =  0.00102
        Solar_Rtm%Tau_O2(2) =  0.01420
        Solar_Rtm%Tau_CO2(2) = 0.00000
        Solar_Rtm%Tau_CH4(2) = 0.00000
        Solar_Rtm%Tau_Ray(6) = 0.00130
        Solar_Rtm%Tau_O3(6) =  0.00000
        Solar_Rtm%Tau_O2(6) =  0.00000
        Solar_Rtm%Tau_CO2(6) = 0.02145
        Solar_Rtm%Tau_CH4(6) = 0.00089


      case(4)      !Metop-B
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00007327, 0.00405108,-0.00012369/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02367473, 0.04800804,-0.00351956/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00002515, 0.00105920,-0.00001039/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00681261, 0.05585163,-0.00315752/)
        Solar_Rtm%Tau_Ray(1)= 0.05582
        Solar_Rtm%Tau_O3(1)=  0.02441
        Solar_Rtm%Tau_O2(1)=  0.00288
        Solar_Rtm%Tau_CO2(1)=-0.00000
        Solar_Rtm%Tau_CH4(1)=-0.00000
        Solar_Rtm%Tau_Ray(2)= 0.01775
        Solar_Rtm%Tau_O3(2)=  0.00104
        Solar_Rtm%Tau_O2(2)=  0.01974
        Solar_Rtm%Tau_CO2(2)= 0.00000
        Solar_Rtm%Tau_CH4(2)=-0.00000
        Solar_Rtm%Tau_Ray(6)= 0.00129
        Solar_Rtm%Tau_O3(6)= -0.00000
        Solar_Rtm%Tau_O2(6)= -0.00000
        Solar_Rtm%Tau_CO2(6)= 0.01928
        Solar_Rtm%Tau_CH4(6)= 0.00113 

      case(5)      !Metop-C - coped from METOP-B
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00007327, 0.00405108,-0.00012369/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02367473, 0.04800804,-0.00351956/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00002515, 0.00105920,-0.00001039/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00681261, 0.05585163,-0.00315752/)
        Solar_Rtm%Tau_Ray(1)= 0.05582
        Solar_Rtm%Tau_O3(1)=  0.02441
        Solar_Rtm%Tau_O2(1)=  0.00288
        Solar_Rtm%Tau_CO2(1)=-0.00000
        Solar_Rtm%Tau_CH4(1)=-0.00000
        Solar_Rtm%Tau_Ray(2)= 0.01775
        Solar_Rtm%Tau_O3(2)=  0.00104
        Solar_Rtm%Tau_O2(2)=  0.01974
        Solar_Rtm%Tau_CO2(2)= 0.00000
        Solar_Rtm%Tau_CH4(2)=-0.00000
        Solar_Rtm%Tau_Ray(6)= 0.00129
        Solar_Rtm%Tau_O3(6)= -0.00000
        Solar_Rtm%Tau_O2(6)= -0.00000
        Solar_Rtm%Tau_CO2(6)= 0.01928
        Solar_Rtm%Tau_CH4(6)= 0.00113 

      case(55)     !Meteosat-8
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00002914, 0.00265955, -0.00007179/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.00694695, 0.02886395, -0.00165606/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00008143, 0.00340042, -0.00010886/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00421946, 0.03202540, -0.00181524/)
        Solar_Rtm%Tau_Ray(1)= 0.05328
        Solar_Rtm%Tau_O3(1)=  0.02283
        Solar_Rtm%Tau_O2(1)=  0.00058
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(2)= 0.02046
        Solar_Rtm%Tau_O3(2)=  0.00127
        Solar_Rtm%Tau_O2(2)=  0.00032
        Solar_Rtm%Tau_CO2(2)= 0.00000
        Solar_Rtm%Tau_CH4(2)= 0.00000
        Solar_Rtm%Tau_Ray(6)= 0.00121
        Solar_Rtm%Tau_O3(6)= 0.00000
        Solar_Rtm%Tau_O2(6)= 0.00000
        Solar_Rtm%Tau_CO2(6)= 0.01303
        Solar_Rtm%Tau_CH4(6)= 0.00603
        Solar_Rtm%Tau_Ray(20) = 0.00000
        Solar_Rtm%Tau_O3(20) =  0.00012
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.09701
        Solar_Rtm%Tau_CH4(20) = 0.00842

      case(56)     !Meteosat-9
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00008122, 0.00268799,-0.00006496/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.00702844, 0.02808596,-0.00161884/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00003371, 0.00364672,-0.00012257/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00416429, 0.03168950,-0.00181495/)
        Solar_Rtm%Tau_Ray(1)= 0.05328     !These are for Meteosat-8
        Solar_Rtm%Tau_O3(1)=  0.02283
        Solar_Rtm%Tau_O2(1)=  0.00058
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(2)= 0.02046
        Solar_Rtm%Tau_O3(2)=  0.00127
        Solar_Rtm%Tau_O2(2)=  0.00032
        Solar_Rtm%Tau_CO2(2)= 0.00000
        Solar_Rtm%Tau_CH4(2)= 0.00000
        Solar_Rtm%Tau_Ray(6)= 0.00121
        Solar_Rtm%Tau_O3(6)= 0.00000
        Solar_Rtm%Tau_O2(6)= 0.00000
        Solar_Rtm%Tau_CO2(6)= 0.01303
        Solar_Rtm%Tau_CH4(6)= 0.00603
        Solar_Rtm%Tau_Ray(20) = 0.00000
        Solar_Rtm%Tau_O3(20) =  0.00012
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.09701
        Solar_Rtm%Tau_CH4(20) = 0.00842

      case(57)     !Meteosat-10
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00008122, 0.00268799,-0.00006496/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.00702844, 0.02808596,-0.00161884/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00003371, 0.00364672,-0.00012257/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00416429, 0.03168950,-0.00181495/)
        Solar_Rtm%Tau_Ray(1) = 0.05396
        Solar_Rtm%Tau_O3(1) =  0.02354
        Solar_Rtm%Tau_O2(1) =  0.00047
        Solar_Rtm%Tau_CO2(1) =-0.00000
        Solar_Rtm%Tau_CH4(1) =-0.00000
        Solar_Rtm%Tau_Ray(2) = 0.02057
        Solar_Rtm%Tau_O3(2) =  0.00129
        Solar_Rtm%Tau_O2(2) =  0.00046
        Solar_Rtm%Tau_CO2(2) =-0.00000
        Solar_Rtm%Tau_CH4(2) =-0.00000
        Solar_Rtm%Tau_Ray(6) = 0.00120
        Solar_Rtm%Tau_O3(6) = -0.00000
        Solar_Rtm%Tau_O2(6) = -0.00000
        Solar_Rtm%Tau_CO2(6) = 0.01191
        Solar_Rtm%Tau_CH4(6) = 0.00626
        Solar_Rtm%Tau_Ray(20) = 0.00000
        Solar_Rtm%Tau_O3(20) =  0.00013
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.09197
        Solar_Rtm%Tau_CH4(20) = 0.00843

      case(171)    !MTSAT-1R
        Solar_Rtm%Tau_H2O_Coef(1,:) = (/  0.000045, 0.002648, -0.000071/)  !MTSAT-1R FAKE
        Solar_Rtm%Tau_Ray(1)= 0.03695
        Solar_Rtm%Tau_O3(1)=  0.01055
        Solar_Rtm%Tau_O2(1)=  0.03034
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000

      case(172)    !MTSAT-2
        Solar_Rtm%Tau_H2O_Coef(1,:) = (/  0.000045, 0.002648, -0.000071/)  !MTSAT-2 FAKE
        Solar_Rtm%Tau_Ray(1)= 0.03695  !These are for MTSAT-1R
        Solar_Rtm%Tau_O3(1)=  0.01055
        Solar_Rtm%Tau_O2(1)=  0.03034
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000

      case(200)    !NOAA-8
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00134225, 0.00605508,-0.00028230/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.01975987, 0.04317701,-0.00308901/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00671906, 0.05136512,-0.00288289/)
        Solar_Rtm%Tau_Ray(1) = 0.05459
        Solar_Rtm%Tau_O3(1) =  0.02241
        Solar_Rtm%Tau_O2(1) =  0.00800
        Solar_Rtm%Tau_CO2(1) =-0.00000
        Solar_Rtm%Tau_CH4(1) =-0.00000
        Solar_Rtm%Tau_Ray(2) = 0.01891
        Solar_Rtm%Tau_O3(2) =  0.00129
        Solar_Rtm%Tau_O2(2) =  0.01530
        Solar_Rtm%Tau_CO2(2) = 0.00000
        Solar_Rtm%Tau_CH4(2) =-0.00000

      case(201)    !NOAA-9
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00035507, 0.00514254,-0.00018633/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02179904, 0.04657735,-0.00340790/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00725621, 0.06165216,-0.00337405/)
        Solar_Rtm%Tau_Ray(1) = 0.05579
        Solar_Rtm%Tau_O3(1) =  0.02227
        Solar_Rtm%Tau_O2(1) =  0.00544
        Solar_Rtm%Tau_CO2(1) = 0.00000
        Solar_Rtm%Tau_CH4(1) = 0.00000
        Solar_Rtm%Tau_Ray(2) = 0.01871
        Solar_Rtm%Tau_O3(2) =  0.00124
        Solar_Rtm%Tau_O2(2) =  0.01015
        Solar_Rtm%Tau_CO2(2) = 0.00000
        Solar_Rtm%Tau_CH4(2) = 0.00000

      case(202)    !NOAA-10
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00016548, 0.00408590,-0.00012677/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02135755, 0.04362566,-0.00313511/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00747077, 0.05748558,-0.00316209/)
        Solar_Rtm%Tau_Ray(1) = 0.05759
        Solar_Rtm%Tau_O3(1) =  0.02430
        Solar_Rtm%Tau_O2(1) =  0.00398
        Solar_Rtm%Tau_CO2(1) = 0.00000
        Solar_Rtm%Tau_CH4(1) = 0.00000
        Solar_Rtm%Tau_Ray(2) = 0.01832
        Solar_Rtm%Tau_O3(2) =  0.00115
        Solar_Rtm%Tau_O2(2) =  0.02603
        Solar_Rtm%Tau_CO2(2) = 0.00000
        Solar_Rtm%Tau_CH4(2) = 0.00000

      case(203)    !NOAA-11
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00034540, 0.00505752,-0.00018436/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02164948, 0.04401914,-0.00304622/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00773716, 0.05937098,-0.00317875/)
        Solar_Rtm%Tau_Ray(1) = 0.0570  !ch1
        Solar_Rtm%Tau_O3(1) = 0.0315
        Solar_Rtm%Tau_O2(1) = 0.0059
        Solar_Rtm%Tau_CH4(1) = 0.00
        Solar_Rtm%Tau_CO2(1) = 0.00
        Solar_Rtm%Tau_Ray(2) = 0.0200  !ch2
        Solar_Rtm%Tau_O3(2) = 0.0012
        Solar_Rtm%Tau_O2(2) = 0.0158
        Solar_Rtm%Tau_CO2(2) = 0.00
        Solar_Rtm%Tau_CH4(2) = 0.00

      case(204)    !NOAA-12
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00051982, 0.00545911,-0.00021035/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02174631, 0.04287286,-0.00300946/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00770047, 0.05084701,-0.00282219/)
        Solar_Rtm%Tau_Ray(1) = 0.0560  !ch1
        Solar_Rtm%Tau_O3(1)= 0.0312
        Solar_Rtm%Tau_O2(1) = 0.0057
        Solar_Rtm%Tau_CH4(1) = 0.00
        Solar_Rtm%Tau_CO2(1) = 0.00
        Solar_Rtm%Tau_Ray(2)= 0.0198  !ch2
        Solar_Rtm%Tau_O3(2) = 0.0011
        Solar_Rtm%Tau_O2(2) = 0.0122
        Solar_Rtm%Tau_CO2(2) = 0.00
        Solar_Rtm%Tau_CH4(2) = 0.00

      case(205)    !NOAA-14
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00042675, 0.00569601,-0.00021510/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02221496, 0.04598416,-0.00330379/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00721637, 0.05361031,-0.00299843/)
        Solar_Rtm%Tau_Ray(1) = 0.0555  !ch1
        Solar_Rtm%Tau_O3(1) = 0.0306
        Solar_Rtm%Tau_O2(1) = 0.0058
        Solar_Rtm%Tau_CH4(1) = 0.00
        Solar_Rtm%Tau_CO2(2) = 0.00
        Solar_Rtm%Tau_Ray(2) = 0.0190  !ch2
        Solar_Rtm%Tau_O3(2) = 0.0009
        Solar_Rtm%Tau_O2(2) = 0.0152
        Solar_Rtm%Tau_CO2(2) = 0.00
        Solar_Rtm%Tau_CH4(2) = 0.00

      case(206)    !NOAA-15
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00008182, 0.00393906,-0.00011248/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02240845, 0.04888891,-0.00372582/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00001497, 0.00106324,-0.00000750/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00787390, 0.06604162,-0.00378195/)
        Solar_Rtm%Tau_Ray(1) = 0.0565  !ch1
        Solar_Rtm%Tau_O3(1) = 0.0341
        Solar_Rtm%Tau_O2(1) = 0.0006
        Solar_Rtm%Tau_CH4(1) = 0.00
        Solar_Rtm%Tau_CO2(1) = 0.00
        Solar_Rtm%Tau_Ray(2) = 0.0190  !ch2
        Solar_Rtm%Tau_O3(2) = 0.0008
        Solar_Rtm%Tau_O2(2) = 0.0137
        Solar_Rtm%Tau_CO2(2) = 0.00
        Solar_Rtm%Tau_CH4(2) = 0.00
        Solar_Rtm%Tau_Ray(6) = 0.0013  !ch3a
        Solar_Rtm%Tau_O3(6) = 0.00
        Solar_Rtm%Tau_O2(6) = 0.00
        Solar_Rtm%Tau_CO2(6)= 0.0161
        Solar_Rtm%Tau_CH4(6) = 0.0006

   case(207)    !NOAA-16
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00006417, 0.00413998,-0.00013151/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02447558, 0.04929511,-0.00365359/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00001981, 0.00115096,-0.00001780/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00774936, 0.06153378,-0.00352855/)
        Solar_Rtm%Tau_Ray(1) = 0.0568  !ch1
        Solar_Rtm%Tau_O3(1) = 0.0342
        Solar_Rtm%Tau_O2(1) = 0.0007
        Solar_Rtm%Tau_CH4(1) = 0.00
        Solar_Rtm%Tau_CO2(1) = 0.00
        Solar_Rtm%Tau_Ray(2) = 0.0187  !ch2
        Solar_Rtm%Tau_O3(2) = 0.0008
        Solar_Rtm%Tau_O2(2) = 0.0137
        Solar_Rtm%Tau_CO2(2) = 0.00
        Solar_Rtm%Tau_CH4(2) = 0.00
        Solar_Rtm%Tau_Ray(6) = 0.0013  !ch3a
        Solar_Rtm%Tau_O3(6) = 0.00
        Solar_Rtm%Tau_O2(6) = 0.00
        Solar_Rtm%Tau_CO2(6) = 0.0161
        Solar_Rtm%Tau_CH4(6) = 0.0004

      case(208)    !NOAA-17
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00006082, 0.00376122,-0.00011058/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02492587, 0.04667360,-0.00332354/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00003421, 0.00107492,-0.00001075/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00777939, 0.05645507,-0.00307462/)
        Solar_Rtm%Tau_Ray(1) = 0.05499
        Solar_Rtm%Tau_O3(1) =  0.02394
        Solar_Rtm%Tau_O2(1) =  0.00272
        Solar_Rtm%Tau_CO2(1) = 0.00000
        Solar_Rtm%Tau_CH4(1) = 0.00000
        Solar_Rtm%Tau_Ray(2) = 0.01767
        Solar_Rtm%Tau_O3(2) =  0.00101
        Solar_Rtm%Tau_O2(2) =  0.01648
        Solar_Rtm%Tau_CO2(2) = 0.00000
        Solar_Rtm%Tau_CH4(2) = 0.00000
        Solar_Rtm%Tau_Ray(6) = 0.00129
        Solar_Rtm%Tau_O3(6) =  0.00000
        Solar_Rtm%Tau_O2(6) =  0.00000
        Solar_Rtm%Tau_CO2(6) = 0.01983
        Solar_Rtm%Tau_CH4(6) = 0.00093

      case(209)    !NOAA-18
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00007524, 0.00374030,-0.00010838/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02537023, 0.05129762,-0.00389744/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00002599, 0.00107186,-0.00000765/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00678576, 0.05510151,-0.00315387/)
        Solar_Rtm%Tau_Ray(1) = 0.05472
        Solar_Rtm%Tau_O3(1) =  0.02363
        Solar_Rtm%Tau_O2(1) =  0.00236
        Solar_Rtm%Tau_CO2(1) = 0.00000
        Solar_Rtm%Tau_CH4(1) = 0.00000
        Solar_Rtm%Tau_Ray(2) = 0.01727
        Solar_Rtm%Tau_O3(2) =  0.00096
        Solar_Rtm%Tau_O2(2) =  0.01568
        Solar_Rtm%Tau_CO2(2) = 0.00000
        Solar_Rtm%Tau_CH4(2) = 0.00000
        Solar_Rtm%Tau_Ray(6) = 0.00129
        Solar_Rtm%Tau_O3(6)=  0.00000
        Solar_Rtm%Tau_O2(6) =  0.00000
        Solar_Rtm%Tau_CO2(6) = 0.02001
        Solar_Rtm%Tau_CH4(6) = 0.00139

      case(223)    !NOAA-19
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00007145, 0.00377615,-0.00010227/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02050573, 0.04243041,-0.00304790/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00001813, 0.00103141,-0.00000540/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00707962, 0.05692835,-0.00310688/)
        Solar_Rtm%Tau_Ray(1) = 0.05436
        Solar_Rtm%Tau_O3(1) =  0.02321
        Solar_Rtm%Tau_O2(1) =  0.00089
        Solar_Rtm%Tau_CO2(1) = 0.00000
        Solar_Rtm%Tau_CH4(1) = 0.00000
        Solar_Rtm%Tau_Ray(2) = 0.01850
        Solar_Rtm%Tau_O3(2) =  0.00114
        Solar_Rtm%Tau_O2(2) =  0.01562
        Solar_Rtm%Tau_CO2(2) = 0.00000
        Solar_Rtm%Tau_CH4(2) = 0.00000
        Solar_Rtm%Tau_Ray(6) = 0.00128
        Solar_Rtm%Tau_O3(6) =  0.00000
        Solar_Rtm%Tau_O2(6) =  0.00000
        Solar_Rtm%Tau_CO2(6) = 0.01979
        Solar_Rtm%Tau_CH4(6) = 0.00109

      case(224)    !NPP-VIIRS
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00001172, 0.00053326,-0.00000720/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.00010465, 0.00210970,-0.00004976/)
        Solar_Rtm%Tau_H2O_Coef(5,:)  = (/  0.00002958, 0.00667158,-0.00023893/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00000041, 0.00102500, 0.00000419/)
        Solar_Rtm%Tau_H2O_Coef(7,:)  = (/ -0.00010325, 0.00118721,-0.00003250/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.01060415, 0.06901925,-0.00374645/)
        Solar_Rtm%Tau_H2O_Coef(44,:) = (/  0.00204467, 0.01079139,-0.00055673/)
        Solar_Rtm%Tau_Ray(1) = 0.04338
        Solar_Rtm%Tau_O3(1) =  0.01210
        Solar_Rtm%Tau_O2(1) =  0.00154
        Solar_Rtm%Tau_CO2(1) = 0.00000
        Solar_Rtm%Tau_CH4(1) = 0.00000
        Solar_Rtm%Tau_Ray(2) = 0.01582
        Solar_Rtm%Tau_O3(2) =  0.00062
        Solar_Rtm%Tau_O2(2) =  0.00003
        Solar_Rtm%Tau_CO2(2) = 0.00000
        Solar_Rtm%Tau_CH4(2) = 0.00000
        Solar_Rtm%Tau_Ray(5) = 0.00368
        Solar_Rtm%Tau_O3(5) =  0.00000
        Solar_Rtm%Tau_O2(5) =  0.01149
        Solar_Rtm%Tau_CO2(5) = 0.00046
        Solar_Rtm%Tau_CH4(5) = 0.00000
        Solar_Rtm%Tau_Ray(6) = 0.00131
        Solar_Rtm%Tau_O3(6) =  0.00000
        Solar_Rtm%Tau_O2(6) =  0.00000
        Solar_Rtm%Tau_CO2(6) = 0.02327
        Solar_Rtm%Tau_CH4(6) = 0.00093
        Solar_Rtm%Tau_Ray(7) = 0.00030
        Solar_Rtm%Tau_O3(7) =  0.00000
        Solar_Rtm%Tau_O2(7) =  0.00000
        Solar_Rtm%Tau_CO2(7) = 0.00000
        Solar_Rtm%Tau_CH4(7) = 0.05392
        Solar_Rtm%Tau_Ray(20) = 0.00001
        Solar_Rtm%Tau_O3(20) =  0.00069
        Solar_Rtm%Tau_O2(20) =  0.00000
        Solar_Rtm%Tau_CO2(20) = 0.00056
        Solar_Rtm%Tau_CH4(20) = 0.02148
        Solar_Rtm%Tau_Ray(44) = 0.04478
        Solar_Rtm%Tau_O3(44) =  0.01100
        Solar_Rtm%Tau_O2(44) =  0.02543
        Solar_Rtm%Tau_CO2(44) = 0.00000
        Solar_Rtm%Tau_CH4(44) = 0.00000


      case(252)    !GOES-8
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00092126, 0.00659908,-0.00029547/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00111326, 0.02145498,-0.00105303/)
        Solar_Rtm%Tau_Ray(1) = 0.06065
        Solar_Rtm%Tau_O3(1) = 0.01994
        Solar_Rtm%Tau_O2(1) = 0.01265
        Solar_Rtm%Tau_CO2(1) = 0.00000
        Solar_Rtm%Tau_CH4(1) = 0.00000
        Solar_Rtm%Tau_Ray(20) =-0.00000
        Solar_Rtm%Tau_O3(20) =  0.00000
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.00341
        Solar_Rtm%Tau_CH4(20) = 0.00682

      case(253)    !GOES-9
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00089897, 0.00678990,-0.00030832/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00045615, 0.01946014,-0.00092195/)
        Solar_Rtm%Tau_Ray(1)= 0.06006
        Solar_Rtm%Tau_O3(1)=  0.01938
        Solar_Rtm%Tau_O2(1)=  0.02214
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(20) =-0.00000
        Solar_Rtm%Tau_O3(20) =  0.00000
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.00358
        Solar_Rtm%Tau_CH4(20) = 0.00676

      case(254)    !GOES-10
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00109442, 0.00768902,-0.00036942/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00055447, 0.01943752,-0.00095680/)
        Solar_Rtm%Tau_Ray(1)= 0.05543
        Solar_Rtm%Tau_O3(1)=  0.01769
        Solar_Rtm%Tau_O2(1)=  0.02268
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(20) =-0.00000
        Solar_Rtm%Tau_O3(20) =  0.00000
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.00375
        Solar_Rtm%Tau_CH4(20) = 0.00656

      case(255)    !GOES-11
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00120161, 0.00811146,-0.00037514/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00088083, 0.02237515,-0.00110001/)
        Solar_Rtm%Tau_Ray(1)= 0.05372
        Solar_Rtm%Tau_O3(1)=  0.01672
        Solar_Rtm%Tau_O2(1)=  0.01615
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(20) =-0.00000
        Solar_Rtm%Tau_O3(20) =  0.00000
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.00344
        Solar_Rtm%Tau_CH4(20) = 0.00697

      case(256)    !GOES-12
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00109431, 0.00766287,-0.00036765/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00164272, 0.02421753,-0.00124759/)
        Solar_Rtm%Tau_Ray(1)= 0.05564
        Solar_Rtm%Tau_O3(1)=  0.01780
        Solar_Rtm%Tau_O2(1)=  0.01214
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(20) =-0.00000
        Solar_Rtm%Tau_O3(20) =  0.00000
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.00332
        Solar_Rtm%Tau_CH4(20) = 0.00705

      case(257)    !GOES-13
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00023748, 0.00478735,-0.00017657/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00236300, 0.02637727,-0.00142497/)
        Solar_Rtm%Tau_Ray(1)= 0.06124
        Solar_Rtm%Tau_O3(1)=  0.02180
        Solar_Rtm%Tau_O2(1)= 0.00593
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(20) = 0.00000
        Solar_Rtm%Tau_O3(20) =  0.00002
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.00920
        Solar_Rtm%Tau_CH4(20) = 0.00731

      case(258)    !GOES-14
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00026800, 0.00469394,-0.00016590/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00300622, 0.02934613,-0.00158632/)
        Solar_Rtm%Tau_Ray(1)= 0.06236
        Solar_Rtm%Tau_O3(1)=  0.02207
        Solar_Rtm%Tau_O2(1)=  0.00647
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(20) = 0.00000
        Solar_Rtm%Tau_O3(20) =  0.00003
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.00868
        Solar_Rtm%Tau_CH4(20) = 0.00790

      case(259)    !GOES-15
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00031354, 0.00469119,-0.00016307/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02535393, 0.05158247,-0.00390970/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00001013, 0.00103950,-0.00000234/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00266527, 0.02674780,-0.00147907/)
        Solar_Rtm%Tau_Ray(1)= 0.06229
        Solar_Rtm%Tau_O3(1)=  0.02202
        Solar_Rtm%Tau_O2(1)=  0.00665
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(20) = 0.00000
        Solar_Rtm%Tau_O3(20) =  0.00002
        Solar_Rtm%Tau_O2(20) = -0.00000
        Solar_Rtm%Tau_CO2(20) = 0.01666
        Solar_Rtm%Tau_CH4(20) = 0.00723

      case(706)    !NOAA-6
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00080858, 0.00517629,-0.00020678/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.01992662, 0.04281095,-0.00306524/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00001879, 0.00104313,-0.00000371/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00739157, 0.05705284,-0.00312311/)
        Solar_Rtm%Tau_Ray(1)= 0.05656
        Solar_Rtm%Tau_O3(1)=  0.02397
        Solar_Rtm%Tau_O2(1)=  0.00034
        Solar_Rtm%Tau_CO2(1)=-0.00000
        Solar_Rtm%Tau_CH4(1)=-0.00000
        Solar_Rtm%Tau_Ray(2)= 0.01858
        Solar_Rtm%Tau_O3(2)=  0.00122
        Solar_Rtm%Tau_O2(2)=  0.03535
        Solar_Rtm%Tau_CO2(2)= 0.00000
        Solar_Rtm%Tau_CH4(2)=-0.00000

      case(707)    !NOAA-7
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00019387, 0.00422749,-0.00014381/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02137234, 0.04687615,-0.00351913/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00001915, 0.00106146,-0.00000921/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00688883, 0.06105546,-0.00349372/)
        Solar_Rtm%Tau_Ray(1)= 0.05691
        Solar_Rtm%Tau_O3(1)=  0.02362
        Solar_Rtm%Tau_O2(1)=  0.00418
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(2)= 0.01824
        Solar_Rtm%Tau_O3(2)=  0.00116
        Solar_Rtm%Tau_O2(2)=  0.03971
        Solar_Rtm%Tau_CO2(2)= 0.00001
        Solar_Rtm%Tau_CH4(2)= 0.00000

      case(708)    !TIROS-N
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00320562, 0.01321214,-0.00075649/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.02192102, 0.04298591,-0.00310794/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00002771, 0.00106427,-0.00000733/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00688623, 0.05329971,-0.00298736/)
        Solar_Rtm%Tau_Ray(1)= 0.04112
        Solar_Rtm%Tau_O3(1)=  0.01015
        Solar_Rtm%Tau_O2(1)=  0.03489
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(2)= 0.01816
        Solar_Rtm%Tau_O3(2)=  0.00167
        Solar_Rtm%Tau_O2(2)=  0.04889
        Solar_Rtm%Tau_CO2(2)= 0.00000
        Solar_Rtm%Tau_CH4(2)= 0.00000

      case(783)    !TERRA-MODIS
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00009297, 0.00401504, -0.00011504/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.00058087, 0.00529998, -0.00023274/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00003207, 0.00097319, -0.00001085/)
        Solar_Rtm%Tau_H2O_Coef(7,:)  = (/  0.00155085, 0.01566729, -0.00073441/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00641049, 0.05349478, -0.00286482/)
        Solar_Rtm%Tau_Ray(1)= 0.05084
        Solar_Rtm%Tau_O3(1)=  0.02009
        Solar_Rtm%Tau_O2(1)=  0.00223
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(2)= 0.01624
        Solar_Rtm%Tau_O3(2)=  0.00066
        Solar_Rtm%Tau_O2(2)=  0.00004
        Solar_Rtm%Tau_CO2(2)= 0.00000
        Solar_Rtm%Tau_CH4(2)= 0.00000
        Solar_Rtm%Tau_Ray(6)= 0.00122
        Solar_Rtm%Tau_O3(6)=  0.00000
        Solar_Rtm%Tau_O2(6)=  0.00000
        Solar_Rtm%Tau_CO2(6)= 0.00588
        Solar_Rtm%Tau_CH4(6)= 0.00612

      case(784)    !AQUA-MODIS - copied from TERRA
        Solar_Rtm%Tau_H2O_Coef(1,:)  = (/  0.00009297, 0.00401504, -0.00011504/)
        Solar_Rtm%Tau_H2O_Coef(2,:)  = (/  0.00058087, 0.00529998, -0.00023274/)
        Solar_Rtm%Tau_H2O_Coef(6,:)  = (/ -0.00003207, 0.00097319, -0.00001085/)
        Solar_Rtm%Tau_H2O_Coef(7,:)  = (/  0.00155085, 0.01566729, -0.00073441/)
        Solar_Rtm%Tau_H2O_Coef(20,:) = (/  0.00641049, 0.05349478, -0.00286482/)
        Solar_Rtm%Tau_Ray(1)= 0.05084
        Solar_Rtm%Tau_O3(1)=  0.02009
        Solar_Rtm%Tau_O2(1)=  0.00223
        Solar_Rtm%Tau_CO2(1)= 0.00000
        Solar_Rtm%Tau_CH4(1)= 0.00000
        Solar_Rtm%Tau_Ray(2)= 0.01624
        Solar_Rtm%Tau_O3(2)=  0.00066
        Solar_Rtm%Tau_O2(2)=  0.00004
        Solar_Rtm%Tau_CO2(2)= 0.00000
        Solar_Rtm%Tau_CH4(2)= 0.00000
        Solar_Rtm%Tau_Ray(6)= 0.00122
        Solar_Rtm%Tau_O3(6)=  0.00000
        Solar_Rtm%Tau_O2(6)=  0.00000
        Solar_Rtm%Tau_CO2(6)= 0.00588
        Solar_Rtm%Tau_CH4(6)= 0.00612

      case(810)    !COMS
        Solar_Rtm%Tau_H2O_Coef(1,:) = (/  0.000045, 0.002648, -0.000071/)  !COMS FAKE
        Solar_Rtm%Tau_H2O_Coef(2,:) = (/  -0.000072, 0.003409, -0.000112/)
        Solar_Rtm%Tau_H2O_Coef(6,:) = (/ 0.00000000, 0.00000000,-0.00000000/)


      case default


      end select

   end subroutine SETUP_SOLAR_RTM
   
   !------------------------------------------------------------------------------
   ! Routine to compute some needed radiative transfer terms for the IR channels
   !
   ! Input:  Chan_Idx - number of the channel being used
   !         Sfc_Idx - level of the surface in the profiles
   !         Profile_Weight - interpolation weight for estimated the surface in the 
   !                          profiles
   !         Sfc_Emiss - emissivity of the surface for this channel
   !         Sfc_Temp - the temperature of the surface
   !         Rad_Atm_Profile - profile of radiance emitted from level to space
   !         Trans_Atm_Profile - profile of transmissio from level to space
   !
   ! Output: Rad_Atm - total radiance due to atmospheric emission
   !         Trans_Atm - total tranmission due to atmosphere
   !         Rad_Atm_Sfc - total radiance due both atmosphere and surface at TOA
   !         Bt_Atm_Sfc - Rad_Atm_Sfc expressed as a brightness temperature
   !------------------------------------------------------------------------------
   subroutine COMPUTE_CHANNEL_ATM_SFC_RAD_BT( &
                Chan_Idx, &
                Sfc_Idx, &
                Profile_Weight, &
                Sfc_Emiss, &
                Sfc_Temp, &
                Rad_Atm_Profile, &
                Trans_Atm_Profile, &
                Rad_Atm, &
                Trans_Atm, &
                Rad_Atm_Sfc, &
                Bt_Atm_Sfc)

      integer, intent(in):: Chan_Idx
      integer, intent(in):: Sfc_Idx
      real, intent(in):: Profile_Weight
      real, intent(in):: Sfc_Emiss
      real, intent(in):: Sfc_Temp
      real, intent(in), dimension(:):: Rad_Atm_Profile
      real, intent(in), dimension(:):: Trans_Atm_Profile
      real, intent(out):: Rad_Atm
      real, intent(out):: Trans_Atm
      real, intent(out):: Rad_Atm_Sfc
      real, intent(out):: Bt_Atm_Sfc
  
      real:: Sfc_Rad 

      Sfc_Rad = Sfc_Emiss * PLANCK_RAD_FAST(Chan_Idx,Sfc_Temp)

      Rad_Atm = Rad_Atm_Profile(Sfc_Idx) +  &
            (Rad_Atm_Profile(Sfc_Idx+1) - Rad_Atm_Profile(Sfc_Idx)) * Profile_Weight

      Trans_Atm = Trans_Atm_Profile(Sfc_Idx) +  &
              (Trans_Atm_Profile(Sfc_Idx+1) - Trans_Atm_Profile(Sfc_Idx)) * Profile_Weight

      Rad_Atm_Sfc = Rad_Atm + Trans_Atm * Sfc_Rad
    
      Bt_Atm_Sfc = PLANCK_TEMP_FAST(Chan_Idx,Rad_Atm_Sfc)

   end subroutine COMPUTE_CHANNEL_ATM_SFC_RAD_BT

   !------------------------------------------------------------------------------
   ! Routine to compute some needed radiative transfer terms for the IR channels
   !
   ! Input:  Sfc_Idx - level of the surface in the profiles
   !         Profile_Weight - interpolation weight for estimated the surface in the 
   !                          profiles
   !         Rad_Atm_Dwn_Profile - profile of radiance emitted from level to space
   !
   ! Output: Rad_Atm_Dwn_Sfc - total radiance due to atmospheric emission at surface
   !------------------------------------------------------------------------------
   subroutine COMPUTE_CHANNEL_ATM_DWN_SFC_RAD( &
                Sfc_Idx, &
                Profile_Weight, &
                Rad_Atm_Dwn_Profile, &
                Rad_Atm_Dwn_Sfc)

      integer, intent(in):: Sfc_Idx
      real, intent(in):: Profile_Weight
      real, intent(in), dimension(:):: Rad_Atm_Dwn_Profile
      real, intent(out):: Rad_Atm_Dwn_Sfc
  
      Rad_Atm_Dwn_Sfc = Rad_Atm_Dwn_Profile(Sfc_Idx) +  &
                    (Rad_Atm_Dwn_Profile(Sfc_Idx+1) - Rad_Atm_Dwn_Profile(Sfc_Idx)) * Profile_Weight

   end subroutine COMPUTE_CHANNEL_ATM_DWN_SFC_RAD

   !----------------------------------------------------------------------------------------
   ! This routine computes radiative transfer terms such as Rad_Atm
   ! Trans_Atm, and clear-sky radinace and brightness temperature
   !
   ! Input:  Sfc_Level_Idx - level just above the surface in the profiles
   !         Prof_Weight - interpolation weight for interpolating to surface level
   !         Lon_Idx = longitude index of NWP cell
   !         Lat_Idx = latitude index of NWP cell
   !         Zen_Idx = zenith angle index of RTM profile
   !
   ! Output:  (note this passed through global arrays) 
   !         Rad_Atm_ChX_Rtm = Radiance Emitted by Atmosphere in Channel X

   !----------------------------------------------------------------------------------------
   ! This routine computes radiative transfer terms such as Rad_Atm
   ! Trans_Atm, and clear-sky radinace and brightness temperature
   !
   ! Input:  Sfc_Level_Idx - level just above the surface in the profiles
   !         Prof_Weight - interpolation weight for interpolating to surface level
   !         Lon_Idx = longitude index of NWP cell
   !         Lat_Idx = latitude index of NWP cell
   !         Zen_Idx = zenith angle index of RTM profile
   !
   ! Output:  (note this passed through global arrays) 
   !         Rad_Atm_ChX_Rtm = Radiance Emitted by Atmosphere in Channel X
   !         Trans_Atm_ChX_Rtm = Transmission by Atmosphere in Channel X
   !         Rad_Clear_ChX_Rtm = Radiance at TOA for clear skies (atm + sfc)
   !         Bt_Clear_ChX_Rtm = BT at TOA for clear skies (atm + sfc)
   !----------------------------------------------------------------------------------------
   subroutine COMPUTE_CHANNEL_RT(Sfc_Level_Idx,Prof_Weight,Lon_Idx,Lat_Idx,Elem_Idx,Line_Idx,Zen_Idx)
    integer, intent(in):: Sfc_Level_Idx
      real, intent(in):: Prof_Weight
      integer, intent(in):: Lon_Idx
      integer, intent(in):: Lat_Idx
      integer, intent(in):: Elem_Idx
      integer, intent(in):: Line_Idx
      integer, intent(in):: Zen_Idx
      integer:: Chan_Idx
      real:: Sfc_Ref
      real:: Rad_Ch20_Temp

      !--------------------------------------------------------------
      ! Solar-Only channels, 1,2,6,7,DNB(44)
      !--------------------------------------------------------------
      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max
         if (Ch(Chan_Idx)%Obs_Type /= SOLAR_OBS_TYPE .and. &
             Ch(Chan_Idx)%Obs_Type /= LUNAR_OBS_TYPE) cycle
         
         select case (Chan_Idx)
         
         case (1,2,5,6,7,44)
            if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
               if (allocated(  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Total_Profile )) then
                  Ch(Chan_Idx)%Trans_Atm_Total(Elem_Idx,Line_Idx) = &
                      Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Ch(Chan_Idx)%Trans_Atm_Total_Profile(Sfc_Level_Idx) +  &
                     (Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Ch(Chan_Idx)%Trans_Atm_Total_Profile(Sfc_Level_Idx+1) -  &
                      Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Ch(Chan_Idx)%Trans_Atm_Total_Profile(Sfc_Level_Idx)) * Prof_Weight
               end if
            end if
         end select

      end do
      
      !--------------------------------------------------------------
      ! IR-only channels, 20-38 (except 26), 42, 43, 45
      !--------------------------------------------------------------

      !--- upwelling
      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max

         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle
    
         if (Ch(Chan_Idx)%Obs_Type /= THERMAL_OBS_TYPE .and. &
             Ch(Chan_Idx)%Obs_Type /= MIXED_OBS_TYPE ) cycle
         
         call COMPUTE_CHANNEL_ATM_SFC_RAD_BT( &
                Chan_Idx, &
                Sfc_Level_Idx, &
                Prof_Weight, &
                Ch(Chan_Idx)%Sfc_Emiss(Elem_Idx,Line_Idx), &
                Tsfc_Nwp_Pix(Elem_Idx,Line_Idx), &
                Rtm(Lon_Idx,Lat_Idx)%D(Zen_Idx)%Ch(Chan_Idx)%Rad_Atm_Profile, &
                Rtm(Lon_Idx,Lat_Idx)%D(Zen_Idx)%Ch(Chan_Idx)%Trans_Atm_Profile, &
                Ch(Chan_Idx)%Rad_Atm(Elem_Idx,Line_Idx), &
                Ch(Chan_Idx)%Trans_Atm(Elem_Idx,Line_Idx), &
                Ch(Chan_Idx)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                Ch(Chan_Idx)%Bt_Toa_Clear(Elem_Idx,Line_Idx))

      end do

      !--- downwelling (only channel 31) 
      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max
         if (Chan_Idx /= 31) cycle    
         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle

         call COMPUTE_CHANNEL_ATM_DWN_SFC_RAD( &
                Sfc_Level_Idx, &
                Prof_Weight, &
                Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%Ch(Chan_Idx)%Rad_Atm_Dwn_Profile, &
                Ch(Chan_Idx)%Rad_Atm_Dwn_Sfc(Elem_Idx,Line_Idx))
      end do

      !--------------------------------------------------------------
      !-- Ch20 Solar
      !--------------------------------------------------------------
      
      if (Sensor%Chan_On_Flag_Default(20) == sym%YES) then

         !--- add in solar component - does not account for glint
         Trans_Atm_Ch20_Solar_Rtm(Elem_Idx,Line_Idx) =  &
                  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Solar_Profile(Sfc_Level_Idx) + &
                 (Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Solar_Profile(Sfc_Level_Idx+1) - &
                  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Solar_Profile(Sfc_Level_Idx)) * &
                  Prof_Weight

         Trans_Atm_Ch20_Solar_Total_Rtm(Elem_Idx,Line_Idx) =  &
                  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Total_Profile(Sfc_Level_Idx) + &
                 (Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Total_Profile(Sfc_Level_Idx+1) - &
                  Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(20)%Trans_Atm_Total_Profile(Sfc_Level_Idx)) * &
                  Prof_Weight

         Rad_Clear_Ch20_Solar_Rtm(Elem_Idx,Line_Idx) = ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx)
         
         Bt_Clear_Ch20_Solar_Rtm(Elem_Idx,Line_Idx) = ch(20)%Bt_Toa_Clear(Elem_Idx,Line_Idx)

         if (Geo%Cossolzen(Elem_Idx,Line_Idx) >= 0.0) then

            Sfc_Ref = 1.0 - ch(20)%Sfc_Emiss(Elem_Idx,Line_Idx)

            Rad_Ch20_Temp = ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx) +  &
                   Trans_Atm_Ch20_Solar_Total_Rtm(Elem_Idx,Line_Idx) *  &
                   Sfc_Ref * (Geo%Cossolzen(Elem_Idx,Line_Idx)*Solar_Ch20_Nu / pi)

            Rad_Clear_Ch20_Solar_Rtm(Elem_Idx,Line_Idx) = Rad_Ch20_Temp
            Bt_Clear_Ch20_Solar_Rtm(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(20,Rad_Ch20_Temp)

         end if
      end if

   end subroutine COMPUTE_CHANNEL_RT
   
   !-------------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------------
   subroutine COMPUTE_CH20_EMISSIVITY(Elem_Idx,Line_Idx)
      integer, intent(in):: Elem_Idx
      integer, intent(in):: Line_Idx
      real:: Rad_Ch20_Temp
      real:: Ch20_Sfc_Rad

      if ((Sensor%Chan_On_Flag_Default(20) == sym%YES) .and. (Sensor%Chan_On_Flag_Default(31)==sym%YES)) then

         Ch20_Sfc_Rad = ch(20)%Sfc_Emiss(Elem_Idx,Line_Idx) * PLANCK_RAD_FAST(20,Tsfc_Nwp_Pix(Elem_Idx,Line_Idx))

         Rad_Ch20_Temp = PLANCK_RAD_FAST(20,ch(31)%Bt_Toa_Clear(Elem_Idx,Line_Idx))
         Ems_Ch20_Clear_Solar_Rtm(Elem_Idx,Line_Idx) = Rad_Clear_Ch20_Solar_Rtm(Elem_Idx,Line_Idx) / Rad_Ch20_Temp
         Ems_Ch20_Clear_Rtm(Elem_Idx,Line_Idx) = ch(20)%Rad_Toa_Clear(Elem_Idx,Line_Idx) / Rad_Ch20_Temp

         !--- atmos corrected
         Ems_Ch20_Clear_Solar_Sfc_Rtm(Elem_Idx,Line_Idx) = Ch20_Sfc_Rad +  &
                     Trans_Atm_Ch20_Solar_Rtm(Elem_Idx,Line_Idx) * &
                     (1.0 - ch(20)%Sfc_Emiss(Elem_Idx,Line_Idx))* &
                     Geo%Cossolzen(Elem_Idx,Line_Idx)*Solar_Ch20_Nu / pi
      end if

   end subroutine COMPUTE_CH20_EMISSIVITY
   
   !-------------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------------
   subroutine COMPUTE_TROPOPAUSE_EMISSIVITIES(Elem_Idx,Line_Idx,Lon_Idx,Lat_Idx,Zen_Idx)
      integer, intent(in):: Elem_Idx
      integer, intent(in):: Line_Idx
      integer, intent(in):: Lon_Idx
      integer, intent(in):: Lat_Idx
      integer, intent(in):: Zen_Idx
      integer:: Chan_Idx
      integer:: Lev_Bnd
      integer:: dim1
      integer:: dim2

      dim1 = Image%Number_Of_Elements
      dim2 = Image%Number_Of_Lines_Per_Segment

      Lev_Bnd = Rtm(Lon_Idx,Lat_Idx)%Tropo_Level  

      !--- check for missing tropopause level
      if (Rtm(Lon_Idx,Lat_Idx)%Tropo_Level == 0) then
         Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
         return
      end if

      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max
      
         select case (Chan_Idx)
         
         case(27,29,31,32,33) 
            if (Sensor%Chan_On_Flag_Default(Chan_Idx)==sym%YES) then

!              if (.not. allocated(ch(Chan_Idx)%Emiss_Tropo)) allocate(ch(Chan_Idx)%Emiss_Tropo(dim1,dim2))
               
               ch(Chan_Idx)%Emiss_Tropo(Elem_Idx,Line_Idx) =  &
                        EMISSIVITY(ch(Chan_Idx)%Rad_Toa(Elem_Idx,Line_Idx),  &
                        ch(Chan_Idx)%Rad_Toa_Clear(Elem_Idx,Line_Idx),  &
                        Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_BB_Cloud_Profile(Lev_Bnd))
            end if
         end select
      end do

   end subroutine COMPUTE_TROPOPAUSE_EMISSIVITIES
   
   !-------------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------------
   subroutine COMPUTE_BETA_RATIOES(Elem_Idx,Line_Idx)
      integer, intent(in):: Elem_Idx
      integer, intent(in):: Line_Idx

      !--- compute 11 and 12 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(32) == sym%YES) then

         Beta_11um_12um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(32)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx))
      end if
      
      !--- compute 11 and 8.5 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
          Sensor%Chan_On_Flag_Default(29) == sym%YES) then

         Beta_11um_85um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(29)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx))
      endif
  
      !--- compute 11 and 6.7 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(27) == sym%YES) then

         Beta_11um_67um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(27)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx))
      end if
      
      !--- compute 11 and 13.3 beta ratio at tropopause
      if (Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
         Sensor%Chan_On_Flag_Default(33) == sym%YES) then

         Beta_11um_133um_Tropo_Rtm(Elem_Idx,Line_Idx) = BETA_RATIO( &
                                        ch(33)%Emiss_Tropo(Elem_Idx,Line_Idx),  &
                                        ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx))
      endif

   end subroutine COMPUTE_BETA_RATIOES
   
   !====================================================================
   ! FUNCTION Name: BETA_RATIO
   !
   ! Function:
   !  Computes the beta ratio for two Emissivities. 
   !
   ! Input:  Emiss_top - emissivity in the numerator
   !         Emiss_bot - emissivity in the denominator
   !
   ! Output: Beta - the beta value from the two emissivities
   !
   !====================================================================
   function BETA_RATIO(Emiss_top, Emiss_bot) result(beta)
      real(kind=real4), intent(in) :: Emiss_top
      real(kind=real4), intent(in) :: Emiss_bot
      real(kind=real4) :: beta

      beta = Missing_Value_Real4

      if (Emiss_top > 0.0 .and. Emiss_top < 1.0 .and. &
            Emiss_bot > 0.0 .and. Emiss_bot < 1.0) then
            
         beta = alog(1.0 - Emiss_top)/alog(1.0 - Emiss_bot)
      end if

      return

   end function BETA_RATIO

   !====================================================================
   ! Function Name: EMISSIVITY
   !
   ! Function:
   !  Computes the  effective emissivity
   !
   ! Input:  Rad_Toa - channel radiance at top of atmosphere(toa)
   !         Rad_Clear_Tau - channel radiance at toa for clear skies
   !         Rad_Cloud_BB_Toa - channel radiance at TOA if cloud were a Black-Body
   !
   ! Output: Emiss - the effective cloud emissivity
   !
   !====================================================================
   function EMISSIVITY(Radiance_Toa, Radiance_Clear_Toa, Radiance_Cloud_BB_Toa) result(Emiss)
      real(kind=real4), intent(in) :: Radiance_Toa
      real(kind=real4), intent(in) :: Radiance_Clear_Toa
      real(kind=real4), intent(in) :: Radiance_Cloud_BB_Toa
      real(kind=real4) :: Emiss

      Emiss = Missing_Value_Real4

      if (Radiance_Cloud_BB_Toa /= Radiance_Clear_Toa) then
          Emiss = (Radiance_Toa - Radiance_Clear_Toa) / &
            (Radiance_Cloud_BB_Toa - Radiance_Clear_Toa) 
       end if

      return

   end function EMISSIVITY
   !====================================================================
   ! Function Name: NADIR_EMISSIVITY
   !
   ! Function: Convert a slant-path emissivity to a nadir value
   !   
   ! Input: Emiss - slant path emissivity
   !        Cos_Zen - cosine of viewing zenith angle
   !
   ! Output: Nadir_Emiss - emissivity if viewed nadir
   !
   !====================================================================
   function NADIR_EMISSIVITY(Emiss,Cos_Zen) result(Nadir_Emiss)
      real(kind=real4), intent(in) :: Emiss
      real(kind=real4), intent(in) :: Cos_Zen
      real(kind=real4) :: Nadir_Emiss
      real:: Tau_Temp

      Nadir_Emiss = Missing_Value_Real4

      if ((Emiss > 0.0) .and. (Emiss < 1.00)) then
         Tau_Temp = -Cos_Zen*alog(1.0 - Emiss)
         Nadir_Emiss = 1.0 - exp(-Tau_Temp)
      end if

      return

   end function NADIR_EMISSIVITY

   !===============================================================================
   !
   !===============================================================================
   subroutine COPY_LOCAL_RTM_TO_GLOBAL_RTM_STRUCTURE(Lon_Idx,Lat_Idx,Zen_Idx)

      integer, intent(in):: Lon_Idx
      integer, intent(in):: Lat_Idx
      integer, intent(in):: Zen_Idx
      integer:: Chan_Idx    

      do Chan_Idx = Chan_Idx_Min, Chan_Idx_Max
         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Total_Profile = Trans_Atm_Total_Prof(:,Chan_Idx)
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Solar_Profile = Trans_Atm_Solar_Prof(:,Chan_Idx)
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Profile = Trans_Atm_Prof(:,Chan_Idx)
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Profile = Rad_Atm_Prof(:,Chan_Idx)
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_BB_Cloud_Profile = Rad_BB_Cloud_Prof(:,Chan_Idx)
      end do

      Chan_Idx = 31
      if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
         Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Dwn_Profile = Rad_Atm_Dwn_Prof(:,Chan_Idx)
      end if

   end subroutine COPY_LOCAL_RTM_TO_GLOBAL_RTM_STRUCTURE
   
   !==============================================================================
   !
   !==============================================================================
   subroutine ALLOCATE_GLOBAL_RTM_STRUCTURE_ELEMENT(Lon_Idx,Lat_Idx,Zen_Idx)

      integer, intent(in):: Lon_Idx
      integer, intent(in):: Lat_Idx
      integer, intent(in):: Zen_Idx
      integer:: Alloc_Status
      integer:: Chan_Idx

      do Chan_Idx = Chan_Idx_Min,Chan_Idx_Max
         if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
            allocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Profile(NLevels_Rtm),stat=Alloc_Status)
            allocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Solar_Profile(NLevels_Rtm),stat=Alloc_Status)
            allocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Total_Profile(NLevels_Rtm),stat=Alloc_Status)
            allocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Profile(NLevels_Rtm),stat=Alloc_Status)
            allocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_BB_Cloud_Profile(NLevels_Rtm),stat=Alloc_Status)
         endif
      enddo

      Chan_Idx = 31
      if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%YES) then
         allocate(Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Rad_Atm_Dwn_Profile(NLevels_Rtm),stat=Alloc_Status)
      endif
   end subroutine ALLOCATE_GLOBAL_RTM_STRUCTURE_ELEMENT


!
!--- end of module
!
end module RT_UTILITIES
