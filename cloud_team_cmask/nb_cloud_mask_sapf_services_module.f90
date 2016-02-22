 !-----------------------------------------------------------------------------
 ! Input Structure
 !-----------------------------------------------------------------------------
 module NB_CLOUD_MASK_SERVICES

 use PCF_NPP_BAYES_CLOUD_MASK_Mod
 use CloudMask_Access_Mod
 use type_KINDS_AIT
 use Convert_Char
 use Error_Messaging_Module
 use GEOCAT_CONSTANTS
 use Fundamental_Constants_Geocat
 use Framework_Global_Variables_Module
 use Sat_Access_Mod
 use LandMask_Access_Mod
 use CoastMask_Access_Mod
 use DesertMask_Access_Mod
 use SfcElev_Access_Mod
 use SfcEmis_Access_Mod
 use SnowMask_Access_Mod
 use NWP_Access_Mod
 use RTM_Data_Access_Mod
 use CloudPhase_Access_Mod
 use Numerical_Routines
 use RTM_MODULE



 implicit none

 include 'nb_cloud_mask.inc'

 type, public :: et_cloudiness_class_type
      integer :: SPACE 
      integer :: MISSING
      integer :: CLOUDY
      integer :: PROB_CLOUDY
      integer :: PROB_CLEAR
      integer :: CLEAR
 end type

 type, public :: mask_input
    integer :: Num_Elem                                        !x-dimension of data arrays
    integer :: Num_Line                                        !number of lines of arrays with data
    integer :: Num_Line_Max                                    !y-dimension of data arrays
    integer(kind=int1) :: Invalid_Data_Mask                    !bad data mask (0=good,1=bad)
    integer :: Chan_On_041um                                   !flag if 0.41um channel on (0=no,1=yes)
    integer :: Chan_On_063um                                   !flag if 0.63um channel on (0=no,1=yes)
    integer :: Chan_On_086um                                   !flag if 0.86um channel on (0=no,1=yes)
    integer :: Chan_On_138um                                   !flag if 1.38um channel on (0=no,1=yes)
    integer :: Chan_On_160um                                   !flag if 1.60um channel on (0=no,1=yes)
    integer :: Chan_On_213um                                   !flag if 2.13um channel on (0=no,1=yes)
    integer :: Chan_On_375um                                   !flag if 3.75um channel on (0=no,1=yes)
    integer :: Chan_On_67um                                    !flag if 6.7um channel on (0=no,1=yes)
    integer :: Chan_On_85um                                    !flag if 8.5um channel on (0=no,1=yes)
    integer :: Chan_On_10um                                    !flag if 10.0um channel on (0=no,1=yes)
    integer :: Chan_On_11um                                    !flag if 11.0um channel on (0=no,1=yes)
    integer :: Chan_On_12um                                    !flag if 12.0um channel on (0=no,1=yes)
    integer :: Chan_On_I1_064um                                !flag if I1 0.64um channel on (0=no,1=yes)
    integer :: Chan_On_I4_374um                                !flag if I4 3.74um channel on (0=no,1=yes)
    integer :: Chan_On_I5_114um                                !flag if I5 11.4um channel on (0=no,1=yes)
    integer :: Chan_On_DNB                                     !flag if DNB channel on (0=no,1=yes)
    integer :: Use_Sounder_11um                                !flag for IFF files where both imager and sounder 11um are available    
    integer(kind=int1) :: Snow_Class                           !Snow Classification 
    integer(kind=int1) :: Land_Class                           !Land Classification
    integer(kind=int1) :: Oceanic_Glint_Mask                   !Mask of oceanic solar glint (0=no,1=yes)
    integer(kind=int1) :: Lunar_Oceanic_Glint_Mask             !Mask of oceanic lunar glint (0=no,1=yes)
    integer(kind=int1) :: Coastal_Mask                         !binary coast mask (0=no,1=yes)
    real(kind=real4) :: Solzen                                 !Solar zenith angle (degrees)
    real(kind=real4) :: Scatzen                                !Solar Scattering angle (degrees)
    real(kind=real4) :: Lunscatzen                             !Lunar Scattering angle (degrees)
    real(kind=real4) :: Senzen                                 !Sensor viewing zenith angle (degrees)
    real(kind=real4) :: Lunzen                                 !Lunar viewing zenith angle (degrees)
    real(kind=real4) :: Lat                                    !Latitude (degrees)
    real(kind=real4) :: Lon                                    !Longitude (degrees)
    real(kind=real4) :: Ref_041um                              !0.41 um toa reflectance (%)
    real(kind=real4) :: Ref_063um                              !0.63 um toa reflectance (%)
    real(kind=real4) :: Ref_063um_Clear                        !0.63 um toa reflectance for clear-sky (%)
    real(kind=real4) :: Ref_063um_Std                          !0.63 um toa reflectance 3x3 Std.  Dev. (%)
    real(kind=real4) :: Ref_063um_Min                          !Min 0.63 um toa reflectance over 3x3 (%)
    real(kind=real4) :: Ref_086um                              !0.86 um toa reflectance (%)
    real(kind=real4) :: Ref_138um                              !1.38 um toa reflectance (%)
    real(kind=real4) :: Ref_160um                              !1.60 um toa reflectance (%)
    real(kind=real4) :: Ref_160um_Clear                        !1.60 um toa reflectance for clear-sky (%)
    real(kind=real4) :: Ref_375um                              !3.75 um toa reflectance (%)
    real(kind=real4) :: Ref_375um_Clear                        !3.75 um toa reflectance for clear-sky (%)
    real(kind=real4) :: Ref_213um                              !2.13 um toa reflectance (%)
    real(kind=real4) :: Bt_375um                               !3.75 um toa brightness temp (K)
    real(kind=real4) :: Bt_375um_Std                           !3.75 um toa brightness temp 3x3 Std. Dev. (K)
    real(kind=real4) :: Emiss_375um                            !3.75 um pseudo toa emissivity
    real(kind=real4) :: Emiss_375um_Clear                      !3.75 um pseudo toa emissivity clear-sky
    real(kind=real4) :: Bt_67um                                !6.7 um toa brightness temperature (K)
    real(kind=real4) :: Bt_85um                                !8.5 um toa brightness temperature (K)
    real(kind=real4) :: Bt_10um                                !10 um toa brightness temperature (K)
    real(kind=real4) :: Bt_11um                                !11 um toa brightness temperature (K)
    real(kind=real4) :: Bt_11um_Sounder                        !11 um toa brightness temp from sounder (K) 
    real(kind=real4) :: Bt_11um_Std                            !11 um toa brightness temp 3x3 Std Dev (K)
    real(kind=real4) :: Bt_11um_Max                            !11 um toa brightness temp 3x3 Max (K)
    real(kind=real4) :: Bt_11um_Clear                          !11 um toa brightness temperature (K)
    real(kind=real4) :: Emiss_11um_Tropo                       !11 um tropo emiss
    real(kind=real4) :: Bt_12um                                !12 um toa brightness temperature (K)
    real(kind=real4) :: Bt_12um_Clear                          !12 um toa bright temp clear-sky (K)
    real(kind=real4) :: Ref_I1_064um_Std                       !2x2 std I1 064 um reflectance (%)
    real(kind=real4) :: Bt_I4_374um_Std                        !2x2 std I4 374 um brightness temp (K)
    real(kind=real4) :: Bt_I5_114um_Std                        !2x2 std I5 114 um brightness temp (K)
    real(kind=real4) :: Bt_11um_Bt_67um_Covar                  !covariance of 11 and 6.7 um bright temp.
    real(kind=real4) :: Sst_Anal_Uni                           !3x3 std of background sst field (K)
    real(kind=real4) :: Emiss_Sfc_375um                        !the surface emissivity at 3.75 um
    real(kind=real4) :: Rad_Lunar                              !Lunar toa radiance from DNB
    real(kind=real4) :: Ref_Lunar                              !Lunar reflectance from DNB (%)
    real(kind=real4) :: Ref_Lunar_Min                          !Min lunar reflectance over 3x3 (%)
    real(kind=real4) :: Ref_Lunar_Std                          !3x3 std dev of lunar ref from DNB (%)
    real(kind=real4) :: Ref_Lunar_Clear                        !Lunar reflectance for clear-skies (%)
    real(kind=real4) :: Zsfc                                   !surface altitude (km)
    integer :: Num_Segments                                    !number of segments in this data 
    integer(kind=int1) :: Solar_Contamination_Mask             !binary mask of solar contamination (0=no,1=yes)
    integer(kind=int1) :: Sfc_Type                             !surface type based on UMD classification

 end type mask_input 
 !-----------------------------------------------------------------------------
 ! Output Structure
 !-----------------------------------------------------------------------------
 type, public :: mask_output
    integer(kind=int1), dimension(NUMBER_OF_FLAG_BYTES) :: Cld_Flags_Packed         !array of packed results 
    integer(kind=int1) :: Cld_Mask_Bayes             !Derived 4-level cloud mask
    integer :: Cloud_Mask_Bayesian_Flag                                       !flag to tell if code should run
    real(kind=real4) :: Posterior_Cld_Probability    !posterior cloud probability (0-1)
    integer(kind=int1) :: Dust_Mask
    integer(kind=int1) :: Smoke_Mask
    integer(kind=int1) :: Fire_Mask
 end type mask_output

 !-----------------------------------------------------------------------------
 ! Diagnostic Output Structure
 !-----------------------------------------------------------------------------
 type, public :: diag_output
    real(kind=real4):: Array_1    !first diagnostic array
    real(kind=real4):: Array_2    !first diagnostic array
    real(kind=real4):: Array_3    !first diagnostic array
 end type diag_output 
 
 !-----------------------------------------------------------------------------
 ! Symbol Structure
 !-----------------------------------------------------------------------------
 type, public :: symbol_naive_bayesian
    integer(kind=int1) :: CLOUDY
    integer(kind=int1) :: PROB_CLOUDY
    integer(kind=int1) :: PROB_CLEAR
    integer(kind=int1) :: CLEAR
    integer(kind=int1) :: BINARY_CLR
    integer(kind=int1) :: BINARY_CLD
    integer(kind=int1) :: MISSING

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
 end type symbol_naive_bayesian


end module NB_CLOUD_MASK_SERVICES
