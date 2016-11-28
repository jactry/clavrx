 !-----------------------------------------------------------------------------
 ! Input Structure
 !-----------------------------------------------------------------------------
 module NB_CLOUD_MASK_SERVICES

   use CONSTANTS
   use PIXEL_COMMON
   use NUMERICAL_ROUTINES




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
    integer:: Num_Elem                                                 !x-dimension of data arrays
    integer:: Num_Line                                                 !number of lines of arrays with data
    integer:: Num_Line_Max                                             !y-dimension of data arrays
    integer(kind=int1), pointer:: Invalid_Data_Mask                    !bad data mask (0=good,1=bad)
    integer:: Chan_On_041um                                            !flag if 0.41um channel on (0=no,1=yes)
    integer:: Chan_On_063um                                            !flag if 0.63um channel on (0=no,1=yes)
    integer:: Chan_On_086um                                            !flag if 0.86um channel on (0=no,1=yes)
    integer:: Chan_On_138um                                            !flag if 1.38um channel on (0=no,1=yes)
    integer:: Chan_On_160um                                            !flag if 1.60um channel on (0=no,1=yes)
    integer:: Chan_On_213um                                            !flag if 2.13um channel on (0=no,1=yes)
    integer:: Chan_On_375um                                            !flag if 3.75um channel on (0=no,1=yes)
    integer:: Chan_On_67um                                             !flag if 6.7um channel on (0=no,1=yes)
    integer:: Chan_On_85um                                             !flag if 8.5um channel on (0=no,1=yes)
    integer:: Chan_On_11um                                             !flag if 11.0um channel on (0=no,1=yes)
    integer:: Chan_On_12um                                             !flag if 12.0um channel on (0=no,1=yes)
    integer:: Chan_On_I1_064um                                         !flag if I1 0.64um channel on (0=no,1=yes)
    integer:: Chan_On_I4_374um                                         !flag if I4 3.74um channel on (0=no,1=yes)
    integer:: Chan_On_I5_114um                                         !flag if I5 11.4um channel on (0=no,1=yes)
    integer :: Use_Sounder_11um                                !flag for IFF files where both imager and sounder 11um are available    
    integer:: Chan_On_DNB                                              !flag if DNB channel on (0=no,1=yes)
    integer(kind=int1), pointer:: Snow_Class                           !Snow Classification 
    integer(kind=int1), pointer:: Land_Class                           !Land Classification
    integer(kind=int1), pointer:: Oceanic_Glint_Mask                   !Mask of oceanic solar glint (0=no,1=yes)
    integer(kind=int1), pointer:: Lunar_Oceanic_Glint_Mask             !Mask of oceanic lunar glint (0=no,1=yes)
    integer(kind=int1), pointer:: Coastal_Mask                         !binary coast mask (0=no,1=yes)
    real(kind=real4), pointer:: Solzen                                 !Solar zenith angle (degrees)
    real(kind=real4), pointer:: Scatzen                                !Solar Scattering angle (degrees)
    real(kind=real4), pointer:: Lunscatzen                             !Lunar Scattering angle (degrees)
    real(kind=real4), pointer:: Senzen                                 !Sensor viewing zenith angle (degrees)
    real(kind=real4), pointer:: Lunzen                                 !Lunar viewing zenith angle (degrees)
    real(kind=real4), pointer:: Lat                                    !Latitude (degrees)
    real(kind=real4), pointer:: Lon                                    !Longitude (degrees)
    real(kind=real4), pointer:: Ref_041um                              !0.41 um toa reflectance (%)
    real(kind=real4), pointer:: Ref_063um                              !0.63 um toa reflectance (%)
    real(kind=real4), pointer:: Ref_063um_Clear                        !0.63 um toa reflectance for clear-sky (%)
    real(kind=real4), pointer:: Ref_063um_Std                          !0.63 um toa reflectance 3x3 Std.  Dev. (%)
    real(kind=real4), pointer:: Ref_063um_Min                          !Min 0.63 um toa reflectance over 3x3 (%)
    real(kind=real4), pointer:: Ref_086um                              !0.86 um toa reflectance (%)
    real(kind=real4), pointer:: Ref_138um                              !1.38 um toa reflectance (%)
    real(kind=real4), pointer:: Ref_160um                              !1.60 um toa reflectance (%)
    real(kind=real4), pointer:: Ref_160um_Clear                        !1.60 um toa reflectance for clear-sky (%)
    real(kind=real4), pointer:: Ref_375um                              !3.75 um toa reflectance (%)
    real(kind=real4), pointer:: Ref_375um_Clear                        !3.75 um toa reflectance for clear-sky (%)
    real(kind=real4), pointer:: Ref_213um                              !2.13 um toa reflectance (%)
    real(kind=real4), pointer:: Bt_375um                               !3.75 um toa brightness temp (K)
    real(kind=real4), pointer:: Bt_375um_Std                           !3.75 um toa brightness temp 3x3 Std. Dev. (K)
    real(kind=real4), pointer:: Emiss_375um                            !3.75 um pseudo toa emissivity
    real(kind=real4), pointer:: Emiss_375um_Clear                      !3.75 um pseudo toa emissivity clear-sky
    real(kind=real4), pointer:: Bt_67um                                !6.7 um toa brightness temperature (K)
    real(kind=real4), pointer:: Bt_85um                                !8.5 um toa brightness temperature (K)
    real(kind=real4), pointer:: Bt_11um                                !11 um toa brightness temperature (K)
    real(kind=real4), pointer:: Bt_11um_Std                            !11 um toa brightness temp 3x3 Std Dev (K)
    real(kind=real4), pointer:: Bt_11um_Max                            !11 um toa brightness temp 3x3 Max (K)
    real(kind=real4), pointer:: Bt_11um_Clear                          !11 um toa brightness temperature (K)
    real(kind=real4) :: Bt_11um_Sounder                        !11 um toa brightness temp from sounder (K) 
    real(kind=real4), pointer:: Emiss_11um_Tropo                       !11 um tropo emiss
    real(kind=real4), pointer:: Bt_12um                                !12 um toa brightness temperature (K)
    real(kind=real4), pointer:: Bt_12um_Clear                          !12 um toa bright temp clear-sky (K)
    real(kind=real4), pointer:: Ref_I1_064um_Std                       !2x2 std I1 064 um reflectance (%)
    real(kind=real4), pointer:: Bt_I4_374um_Std                        !2x2 std I4 374 um brightness temp (K)
    real(kind=real4), pointer:: Bt_I5_114um_Std                        !2x2 std I5 114 um brightness temp (K)
    real(kind=real4), pointer:: Bt_11um_Bt_67um_Covar                  !covariance of 11 and 6.7 um bright temp.
    real(kind=real4), pointer:: Sst_Anal_Uni                           !3x3 std of background sst field (K)
    real(kind=real4), pointer:: Emiss_Sfc_375um                        !the surface emissivity at 3.75 um
    real(kind=real4), pointer:: Rad_Lunar                              !Lunar toa radiance from DNB
    real(kind=real4), pointer:: Ref_Lunar                              !Lunar reflectance from DNB (%)
    real(kind=real4), pointer:: Ref_Lunar_Min                          !Min lunar reflectance over 3x3 (%)
    real(kind=real4), pointer:: Ref_Lunar_Std                          !3x3 std dev of lunar ref from DNB (%)
    real(kind=real4), pointer:: Ref_Lunar_Clear                        !Lunar reflectance for clear-skies (%)
    real(kind=real4), pointer:: Zsfc                                   !surface altitude (km)
    integer:: Num_Segments                                             !number of segments in this data 
    integer(kind=int1), pointer:: Solar_Contamination_Mask             !binary mask of solar contamination (0=no,1=yes)
    integer(kind=int1), pointer:: Sfc_Type                             !surface type based on UMD classification
 end type mask_input 
 !-----------------------------------------------------------------------------
 ! Output Structure
 !-----------------------------------------------------------------------------
 type, public :: mask_output
    integer(kind=int1), dimension(:), pointer:: Cld_Flags_Packed         !array of packed results 
    integer(kind=int1), pointer:: Cld_Mask_Bayes             !Derived 4-level cloud mask
    integer:: Cloud_Mask_Bayesian_Flag                                       !flag to tell if code should run
    real(kind=real4), pointer:: Posterior_Cld_Probability    !posterior cloud probability (0-1)
    integer(kind=int1), pointer:: Dust_Mask
    integer(kind=int1), pointer:: Smoke_Mask
    integer(kind=int1), pointer:: Fire_Mask
    integer(kind=int1), pointer:: Thin_Cirr_Mask
 end type mask_output

 !-----------------------------------------------------------------------------
 ! Diagnostic Output Structure
 !-----------------------------------------------------------------------------
 type, public :: diag_output
    real(kind=real4), pointer:: Array_1    !first diagnostic array
    real(kind=real4), pointer:: Array_2    !first diagnostic array
    real(kind=real4), pointer:: Array_3    !first diagnostic array
 end type diag_output 
 
 !-----------------------------------------------------------------------------
 ! Symbol Structure
 !-----------------------------------------------------------------------------
 type, public :: symbol_naive_bayesian
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
 end type symbol_naive_bayesian


end module NB_CLOUD_MASK_SERVICES
