!$Id$
!----------------------------------------------------------------------
! MODULE name: NB_CLOUD_MASK
! 
! Routines for the determination of the naive Bayesian cloud mask
!
! Authors: Andrew Heidinger, NOAA/NESDIS
!          Andi Walther, CIMSS
!          Denis Botambekov, CIMSS
!          William Straka, CIMSS
!
! DEPENDENCIES: Services_Module
!
! SIDE EFFECTS: None
!
! Cld_Flags Format
!
! Flag Number  bit-depth    bits      byte name
! 1            1            1         1    cloud mask attempted       
! 2            1            2         1    day_063     
! 3            1            3         1    day_063_spatial     
! 4            1            4         1    day_375
! 5            1            5         1    night_375
! 6            1            6         1    solar contamination 
! 7            1            7         1    coast
! 8            1            8         1    mountain
!----
! 9            1            9         2    forward scattering
! 10           1            10        2    cold scene 375
! 11           1            11        2    cold scene btd
! 12           1            12        2    glint
! 13           1            13        2    smoke detected
! 14           1            14        2    dust  detected
! 15           1            15        2    shadow detected
! 16           1            16        2    fire detected
!---
! 17           3            17-19     3    nbvcm surface type
! 18           1            20        3    thin cirrus detected (separate test)
!<-------------------- START OF CLOUD TESTS -------------------------->
! 19           2            21-22     3    T_11           (TGCT)
! 20           2            23-24     3    T_Max-T        (RTCT)
! ---
! 21           2            25-26     4    T_std          (TUT)
! 22           2            27-28     4    Emiss_Tropo
! 23           2            29-30     4    FMFT mask (Btd_11_12)
! 24           2            31-32     4    Btd_11_67 
!---
! 25           2            33-34     5    Bt_11_67_Covar 
! 26           2            35-36     5    Btd_11_85
! 27           2            37-38     5    Btd_375_11_All
! 28           2            39-40     5    Btd_375_11_Day
!---
! 29           2            41-42     6    Btd_375_11_Night
! 30           2            43-44     6    Spare
! 31           2            45-46     6    Ref_063_Day (RGCT)
! 32           2            47-48     6    Ref_std      (RUT)
!---
! 33           2            49-50     7    Ref_063_Min_3x3 (RVCT)
! 34           2            51-52     7    Ref_Ratio_Day   (RRCT)
! 35           2            53-54     7    Ref_138_Day     (CIRREF)
! 36           2            55-56     7    Ndsi_Day        (NIRREF)
!
!----------------------------------------------------------------------
module NB_CLOUD_MASK

 use NB_CLOUD_MASK_SERVICES
 use NB_CLOUD_MASK_LUT_MODULE

 implicit none

 private:: COMPUTE_BAYES_SFC_TYPE
 private:: split_window_test
 private:: reflectance_gross_contrast_test
 private:: relative_visible_contrast_test
 private:: reflectance_ratio_test
 private:: nir_reflectance_gross_contrast_test
 private:: emiss_375um_test
 private:: emiss_375um_day_test
 private:: emiss_375um_night_test
 private:: PACK_BITS_INTO_BYTES

 public:: NB_CLOUD_MASK_ALGORITHM
 public:: SET_CLOUD_MASK_VERSION
 public:: SET_CLOUD_MASK_THRESHOLDS_VERSION

 !--- set thresholds and algorithm specific constants
 include 'nb_cloud_mask.inc'

 ! --- set cloud mask probability thresholds based on sfc type
 real, dimension (7), parameter :: CONF_CLEAR_PROB_CLEAR_THRESH =  [0.01, 0.01, 0.01, 0.10, 0.10, 0.10, 0.10]
 real, dimension (7), parameter :: PROB_CLEAR_PROB_CLOUD_THRESH =  [0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50]
 real, dimension (7), parameter :: PROB_CLOUDY_CONF_CLOUD_THRESH = [0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90]

 !--- string to control on-screen prompts
 character(*), parameter, private :: EXE_PROMPT_CM = "Naive Bayesian Cloud Mask>> "

 !--------------------------------------------------------------------------------
 ! store table values as module-wide arrays
 !--------------------------------------------------------------------------------

 real, dimension(:), allocatable, private, save:: Cond_Ratio
 real, dimension(:), allocatable, private, save:: Posterior_Cld_Probability_By_Class
 real, dimension(:), allocatable, private, save:: Classifier_Value

 ! netCDF parameters
   integer, parameter, private :: sds_rank_1d = 1
   integer, dimension(sds_rank_1d), private :: sds_start_1d, sds_edge_1d, &
          sds_stride_1d

   integer, parameter, private :: sds_rank_2d = 2
   integer, dimension(sds_rank_2d), private :: sds_start_2d, sds_edge_2d, &
          sds_stride_2d, sds_dims_2d

   integer, parameter, private :: sds_rank_3d = 3
   integer, dimension(sds_rank_3d), private :: sds_start_3d, sds_edge_3d, &
          sds_stride_3d, sds_dims_3d

 contains
!====================================================================
!  record cvs version as a global variable for output to hdf
!====================================================================
 subroutine SET_CLOUD_MASK_VERSION(Cloud_Mask_Version)
   character(len=*), intent(out):: Cloud_Mask_Version
   Cloud_Mask_Version = "$Id$"
 end subroutine SET_CLOUD_MASK_VERSION

!====================================================================
!  pass threshold version to bridge
!  Cloud_Mask_Thresholds_Version is a module-wide variable
!  that is private
!====================================================================
 subroutine SET_CLOUD_MASK_THRESHOLDS_VERSION(Cloud_Mask_CVS_Tag)
   character(len=*), intent(out):: Cloud_Mask_CVS_Tag
   Cloud_Mask_CVS_Tag = Cloud_Mask_Thresholds_Version
 end subroutine SET_CLOUD_MASK_THRESHOLDS_VERSION

!====================================================================
! SUBROUTINE Name: CLOUD_MASK_NAIVE_BAYES
!
! Function:
!   CalcuLates the bayesian cloud mask. The bayesian cloud mask is
!   determined by utilizing the following surface types:
!
! Bayesian Surface Types
! 1 - Deep_Water
! 2 - Shallow_Water
! 3 - Unfrozen_Land
! 4 - Frozen_Land
! 5 - Arctic
! 6 - Antarctic
! 7 - Desert
!
!====================================================================
 subroutine NB_CLOUD_MASK_ALGORITHM( &
            Naive_Bayes_File_Name_Full_Path, &       !full path & name of nb classifer data file
            Symbol, &                                !local copy of sym structure
            Input,  &
            Output,  &
            Use_Prior_Table, &
            Diag)


   character (len=*), intent(in) :: Naive_Bayes_File_Name_Full_Path
   type(symbol_naive_bayesian), intent(in) :: Symbol
   type(mask_input), intent(in) :: Input
   type(mask_output), intent(out) :: Output
   logical, intent(in):: Use_Prior_Table
   type(diag_output), intent(out), Optional :: Diag

   !-- local pointers that point to global variables
   integer (kind=INT1), dimension(NUMBER_OF_FLAGS):: Cld_Flags
   integer (kind=INT1), dimension(NUMBER_OF_FLAGS):: Cld_Flag_Bit_Depth

   !internal variables
   integer:: Class_Idx
   integer:: Sfc_Idx
   integer:: Bin_Idx

   integer:: Oceanic_Glint_Flag
   integer:: Lunar_Oceanic_Glint_Flag
   integer:: Coastal_Flag
   integer:: Mountain_Flag
   integer:: Day_063_Flag
   integer:: Day_063_Spatial_Flag
   integer:: Night_Lunar_Flag
   integer:: Lunar_Spatial_Flag
   integer:: Day_375_Flag
   integer:: Night_375_Flag
   integer:: Forward_Scattering_Flag
   integer:: Solar_Contam_Flag
   integer:: Lunar_Forward_Scattering_Flag
   integer:: Cold_Scene_375um_Flag
   integer:: Cold_Scene_Flag
   integer:: Dry_Scene_Flag
   integer:: City_Flag

   real (kind=real4):: Airmass
   integer, parameter:: Spare_Value = 0
   
   real:: r
   real:: Prior_Yes_Temp

   !------------------------------------------------------------------------------------------
   !---  begin executable code
   !------------------------------------------------------------------------------------------
   
   !------------------------------------------------------------------------------------------
   !--- allocate memory for local arrays
   !------------------------------------------------------------------------------------------
   if (allocated(Cond_Ratio)) deallocate(Cond_Ratio)
   if (allocated(Classifier_Value)) deallocate(Classifier_Value)
   if (allocated(Posterior_Cld_Probability_By_Class)) deallocate(Posterior_Cld_Probability_By_Class)

   allocate(Cond_Ratio(N_Class)) 
   allocate(Posterior_Cld_Probability_By_Class(N_Class)) 
   allocate(Classifier_Value(N_Class)) 

   !--- Initialize flags
   Cld_Flags = 0
   Cld_Flag_Bit_Depth = 2     ! needs to be 2 for reserving space for missing tests

   !--- initialize diagnostic output
   if (present(Diag)) Diag%Array_1 = Missing_Value_Real4
   if (present(Diag)) Diag%Array_2 = Missing_Value_Real4
   if (present(Diag)) Diag%Array_3 = Missing_Value_Real4

   !------------------------------------------------------------------------------------------
   !  MAIN CODE
   !------------------------------------------------------------------------------------------

   !--- initialize output
   Output%Posterior_Cld_Probability = Missing_Value_Real4
   Output%Cld_Mask_Bayes = Missing_Value_Int1
   Output%Cld_Flags_Packed = 0

   Cld_Flags(1) = symbol%YES          ;    Cld_Flag_Bit_Depth(1) = 1
   !--- check for a bad pixel
   if (Input%Invalid_Data_Mask == symbol%YES) then
      Cld_Flags(1) = symbol%NO
   else

      !---  COMPUTE SURFACE TYPE
      Sfc_Idx = COMPUTE_BAYES_SFC_TYPE(Input%Land_Class, &
                                       Input%Coastal_Mask, &
                                       Input%Snow_Class, &
                                       Input%Sfc_Type, &
                                       Input%Lat, &
                                       Input%Lon, &
                                       Input%Sst_Anal_Uni, &
                                       Input%Emiss_Sfc_375um, &
                                       symbol)


      !----- compute prior
      if (Input%Prior /= MISSING_VALUE_REAL4 .and. Use_Prior_Table .eqv. .true.) then
         Prior_Yes_Temp = Input%Prior
      else
         Prior_Yes_Temp = Prior_Yes(Sfc_Idx)
      endif

      !--- check for a surface types without sufficient training for processing
      if (Skip_Sfc_Type_Flag(Sfc_Idx) == symbol%YES) then
          Cld_Flags(1) = symbol%NO
      else

          !--- set some flags to control processing
          Oceanic_Glint_Flag = Input%Oceanic_Glint_Mask
          Coastal_Flag = Input%Coastal_Mask
          Solar_Contam_Flag = Input%Solar_Contamination_Mask

          !--- compute airmass
          AirMass = 1.0/cos(Input%Solzen*dtor) + 1.0 / cos(Input%Senzen*dtor)

          !--- set day flag for 0.63 micron reflectance gross test
          if ((Input%Solzen > Reflectance_Gross_Solzen_Thresh) .or.  &
              (AirMass > Reflectance_Gross_Airmass_Thresh)) then
              Day_063_Flag = symbol%NO
          else
              Day_063_Flag = symbol%YES
          endif

          !--- set day flag for 0.63 micron reflectance spatial test
          if (Input%Solzen > Reflectance_Spatial_Solzen_Thresh) then
              Day_063_Spatial_Flag = symbol%NO
          else
              Day_063_Spatial_Flag = symbol%YES
          endif

          Lunar_Spatial_Flag = symbol%NO
          Night_Lunar_Flag = symbol%NO
          Lunar_Forward_Scattering_Flag = symbol%NO

          if (Input%Chan_On_DNB == symbol%YES) then
             Lunar_Oceanic_Glint_Flag = Input%Lunar_Oceanic_Glint_Mask
             Lunar_Spatial_Flag = symbol%YES
             if (Input%Lunzen > Reflectance_Spatial_Solzen_Thresh) then
               Lunar_Spatial_Flag = symbol%NO
             endif

             Night_Lunar_Flag = symbol%YES
             if (Input%Lunzen > Reflectance_Gross_Lunzen_Thresh) then
               Night_Lunar_Flag = symbol%NO
             endif

             Lunar_Forward_Scattering_Flag = symbol%NO
             if (Input%Lunscatzen < Forward_Scatter_Scatzen_Max_Thresh .and. &
                 Input%Lunzen < Forward_Scatter_Solzen_Max_Thresh) then
               Lunar_Forward_Scattering_Flag = symbol%YES
             endif

          endif

          if ((Input%Solzen > Emiss_375um_Day_Solzen_Thresh) .or. &
              (AirMass > Reflectance_Gross_Airmass_Thresh)) then
              Day_375_Flag = symbol%NO
          else
              Day_375_Flag = symbol%YES
          endif

          if (Input%Solzen < Emiss_375um_Night_Solzen_Thresh) then
              Night_375_Flag = symbol%NO
          else
              Night_375_Flag = symbol%YES
          endif

          if (Sfc_Idx /= 6 .and. Input%Zsfc > 2000.0) then
              Mountain_Flag = symbol%YES
          else
              Mountain_Flag = symbol%NO
          endif

!         if (Input%Scatzen < Forward_Scatter_Scatzen_Max_Thresh .and. &
!              Input%Solzen < Forward_Scatter_Solzen_Max_Thresh) then
! ---> DENIS TEST
          if (Input%Scat_Mask == 1) then
              Forward_Scattering_Flag = symbol%YES
          else
              Forward_Scattering_Flag = symbol%NO
          endif

          Cold_Scene_375um_Flag = symbol%NO
          if (Input%Chan_On_375um == symbol%YES) then
           if (Input%Bt_375um < Bt_375um_Cold_Scene_Thresh .and.  &
               Input%Bt_375um /= Missing_Value_Real4) then
              Cold_Scene_375um_Flag = symbol%YES
           endif
          endif

          Cold_Scene_Flag = symbol%NO
          if (Input%Sfc_Temp < Tsfc_Cold_Scene_Thresh) then
              Cold_Scene_Flag = symbol%YES
          endif

!         if (Input%Chan_On_11um == symbol%YES) then
!          if (Input%Bt_11um < Bt_11um_Cold_Scene_Thresh .and.  &
!              Input%Bt_11um /= Missing_Value_Real4) then
!             Cold_Scene_Btd_Flag = symbol%YES
!          endif
!         endif

          Dry_Scene_Flag = symbol%NO
          if (Input%Path_Tpw < Path_Tpw_Dry_Scene_Thresh) then
              Dry_Scene_Flag = symbol%YES
          endif

          !--- City Flag of DNB Lunar Tests
          City_Flag = symbol%NO 
          if (Input%Chan_On_DNB == symbol%YES) then
            if (Input%Rad_Lunar > Radiance_Lunar_City_Thresh) City_Flag = symbol%YES
          endif

          !----------------------------------------------------------------------------------
          !--- populate elements of Cld_Flags with processing flags
          !----------------------------------------------------------------------------------
          Cld_Flags(2) = Day_063_Flag          ;    Cld_Flag_Bit_Depth(2) = 1
          Cld_Flags(3) = Day_063_Spatial_Flag  ;    Cld_Flag_Bit_Depth(3) = 1
          Cld_Flags(4) = Day_375_Flag          ;    Cld_Flag_Bit_Depth(4) = 1
          Cld_Flags(5) = Night_375_Flag        ;    Cld_Flag_Bit_Depth(5) = 1
          Cld_Flags(6) = Solar_Contam_Flag     ;    Cld_Flag_Bit_Depth(6) = 1
          Cld_Flags(7) = Coastal_Flag          ;    Cld_Flag_Bit_Depth(7) = 1
          Cld_Flags(8) = Mountain_Flag         ;    Cld_Flag_Bit_Depth(8) = 1
          Cld_Flags(9) = Forward_Scattering_Flag;   Cld_Flag_Bit_Depth(9) = 1
          Cld_Flags(10) = Cold_Scene_375um_Flag ;   Cld_Flag_Bit_Depth(10) = 1
          Cld_Flags(11) = Cold_Scene_Flag      ;    Cld_Flag_Bit_Depth(11) = 1
          Cld_Flags(12) = Oceanic_Glint_Flag   ;    Cld_Flag_Bit_Depth(12) = 1
          Cld_Flags(13) = Spare_Value          ;    Cld_Flag_Bit_Depth(13) = 1
          Cld_Flags(14) = Spare_Value          ;    Cld_Flag_Bit_Depth(14) = 1
          Cld_Flags(15) = Spare_Value          ;    Cld_Flag_Bit_Depth(15) = 1
          Cld_Flags(16) = Spare_Value          ;    Cld_Flag_Bit_Depth(16) = 1
          Cld_Flags(17) = Sfc_Idx              ;    Cld_Flag_Bit_Depth(17) = 3
          Cld_Flags(18) = Spare_Value          ;    Cld_Flag_Bit_Depth(18) = 1


          !---------------------------------------------------------------------------------
          ! compute cloud probability for 1d NB classifiers
          !---------------------------------------------------------------------------------
          class_loop: do Class_Idx = 1, N_class

             Classifier_Value(Class_Idx) = Missing_Value_Real4
             Cond_Ratio(Class_Idx) =  1.0

             select case (trim(Classifier_Value_Name(Class_Idx,Sfc_Idx)))

                    case("T_11") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Bt_11um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um

                    case("T_Max-T") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Coastal_Flag == symbol%YES) cycle
                       if (Input%Bt_11um == Missing_Value_Real4) cycle
                       if (Input%Bt_11um_Max == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um_Max - Input%Bt_11um

                    case("T_Std") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Coastal_Flag == symbol%YES) cycle
                       if (Input%Bt_11um_Std == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um_Std

                    case("Emiss_Tropo") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Emiss_11um_Tropo == Missing_Value_Real4) cycle
                       if (Cold_Scene_Flag == symbol%YES) cycle
                       Classifier_Value(Class_Idx) = Input%Emiss_11um_Tropo

                    case("FMFT") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_12um == symbol%NO) cycle
                       if (Cold_Scene_Flag == symbol%YES) cycle
                       if (Input%Bt_11um == Missing_Value_Real4) cycle
                       if (Input%Bt_11um_Clear == Missing_Value_Real4) cycle
                       if (Input%Bt_12um == Missing_Value_Real4) cycle
                       if (Input%Bt_12um_Clear == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) =  &
                           split_window_test(Input%Bt_11um_Clear, Input%Bt_12um_Clear, &
                                             Input%Bt_11um, Input%Bt_12um)

                    case("Btd_11_67") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_67um == symbol%NO) cycle
                       if (Input%Bt_11um == Missing_Value_Real4) cycle
                       if (Input%Bt_67um == Missing_Value_Real4) cycle
                       if (Sfc_Idx == 1) cycle
                       if (Sfc_Idx == 2) cycle
                       if (Sfc_Idx == 3) cycle
                       if (Sfc_Idx == 7) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um - Input%Bt_67um
                       if (Input%Use_Sounder_11um == symbol%YES) then
                         Classifier_Value(Class_Idx) = Input%Bt_11um_Sounder - Input%Bt_67um
                       endif

                    case("Bt_11_67_Covar") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_67um == symbol%NO) cycle
                       if (Dry_Scene_Flag == symbol%YES) cycle
                       if (Input%Bt_11um_Bt_67um_Covar == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um_Bt_67um_Covar

                    case("Btd_11_85") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_85um == symbol%NO) cycle
                       if (Cold_Scene_Flag == symbol%YES) cycle
                       if (Input%Bt_11um == Missing_Value_Real4) cycle
                       if (Input%Bt_85um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um - Input%Bt_85um

                    case("Emiss_375") 
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       if (Input%Bt_375um == Missing_Value_Real4) cycle
                       if (Input%Emiss_375um == Missing_Value_Real4) cycle
                       if (Input%Emiss_375um_Clear == Missing_Value_Real4) cycle

                       Classifier_Value(Class_Idx) = emiss_375um_test( &
                                         Input%Emiss_375um,Input%Emiss_375um_Clear)

                    case("Emiss_375_Day") 
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       if (Day_375_Flag == symbol%NO) cycle
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       if (Input%Bt_375um == Missing_Value_Real4) cycle
                       if (Input%Emiss_375um == Missing_Value_Real4) cycle
                       if (Input%Emiss_375um_Clear == Missing_Value_Real4) cycle

                       Classifier_Value(Class_Idx) = emiss_375um_day_test( &
                                         Input%Emiss_375um,Input%Emiss_375um_Clear)

                    case("Emiss_375_Night") 
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                       if (Night_375_Flag == symbol%NO) cycle
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       if (Input%Bt_375um == Missing_Value_Real4) cycle
                       if (Input%Emiss_375um == Missing_Value_Real4) cycle
                       if (Input%Emiss_375um_Clear == Missing_Value_Real4) cycle

                       Classifier_Value(Class_Idx) = emiss_375um_night_test( &
                                        Input%Emiss_375um,Input%Emiss_375um_Clear)

                     case("Btd_375_11_All") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                      !if (All_375_Flag == symbol%NO) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       if (Forward_Scattering_Flag == symbol%YES) cycle
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       if (Cold_Scene_Flag == symbol%YES) cycle
                       if (Input%Bt_375um == Missing_Value_Real4) cycle
                       if (Input%Bt_11um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = (Input%Bt_375um - Input%Bt_11um) - &
                                                     (Input%Bt_375um_Clear - Input%Bt_11um_Clear)
                     case("Btd_375_11_Day") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                       if (Day_375_Flag == symbol%NO) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       if (Cold_Scene_Flag == symbol%YES) cycle
                       if (Forward_Scattering_Flag == symbol%YES) cycle
                       if (Input%Bt_375um == Missing_Value_Real4) cycle
                       if (Input%Bt_11um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = (Input%Bt_375um - Input%Bt_11um) - &
                                                     (Input%Bt_375um_Clear - Input%Bt_11um_Clear)

                     case("Btd_375_11_Night") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                       if (Night_375_Flag == symbol%NO) cycle
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       if (Cold_Scene_Flag == symbol%YES) cycle
                       if (Input%Bt_375um == Missing_Value_Real4) cycle
                       if (Input%Bt_11um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = (Input%Bt_375um - Input%Bt_11um) - &
                                                     (Input%Bt_375um_Clear - Input%Bt_11um_Clear)

                    case("Ref_063_Day")
                       if (Input%Solzen < 90.0) then
                         if (Input%Chan_On_063um == symbol%NO) cycle
                         if (Oceanic_Glint_Flag == symbol%YES) cycle
                         if (Forward_Scattering_Flag == symbol%YES) cycle
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Day_063_Flag == symbol%NO) cycle
                         if (Sfc_Idx == 4) cycle
                         if (Input%Ref_063um_Clear == Missing_Value_Real4) cycle
                         if (Input%Ref_063um == Missing_Value_Real4) cycle
                         Classifier_Value(Class_Idx) =  &
                             reflectance_gross_contrast_test(Input%Ref_063um_Clear, Input%Ref_063um)

                       else
                         if (Input%Chan_On_DNB == symbol%NO) cycle
                         if (Lunar_Oceanic_Glint_Flag == symbol%YES) cycle
                         if (Lunar_Forward_Scattering_Flag == symbol%YES) cycle
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Night_Lunar_Flag == symbol%NO) cycle
                         if (City_Flag == symbol%YES) cycle
                         if (Sfc_Idx == 4) cycle
                         if (Input%Ref_Lunar_Clear == Missing_Value_Real4) cycle
                         if (Input%Ref_Lunar == Missing_Value_Real4) cycle
                         if (Input%Ref_Lunar_Clear <= 0.) cycle
                         if (Input%Ref_Lunar <= 0.0) cycle
                         Classifier_Value(Class_Idx) =  &
                             reflectance_gross_contrast_test(Input%Ref_Lunar_Clear, &
                                                             Input%Ref_Lunar)
                       endif

                    case("Ref_Std")
                       if (Input%Solzen < 90.0) then
                         if (Input%Chan_On_063um == symbol%NO) cycle
                         if (Day_063_Spatial_Flag == symbol%NO) cycle 
!                         if (Forward_Scattering_Flag == symbol%YES) cycle ! ---> DENIS TEST
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Coastal_Flag == symbol%YES) cycle
                         if (Input%Ref_063um_Std == Missing_Value_Real4) cycle
                         Classifier_Value(Class_Idx) = Input%Ref_063um_Std
                       else
                         if (Input%Chan_On_DNB == symbol%NO) cycle
                         if (Lunar_Spatial_Flag == symbol%NO) cycle
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Coastal_Flag == symbol%YES) cycle
                         if (City_Flag == symbol%YES) cycle  
                         if (Input%Ref_Lunar_Std == Missing_Value_Real4) cycle
                         if (Input%Ref_Lunar_Std <= 0.0) cycle
                         Classifier_Value(Class_Idx) = Input%Ref_Lunar_Std
                       endif

                    case("Ref_063_Min_3x3_Day")
                       if (Input%Solzen < 90.0) then
                         if (Input%Chan_On_063um == symbol%NO) cycle
                         if (Day_063_Spatial_Flag == symbol%NO) cycle
!                         if (Forward_Scattering_Flag == symbol%YES) cycle ! ---> DENIS TEST
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Coastal_Flag == symbol%YES) cycle
                         if (Input%Ref_063um_Min == Missing_Value_Real4) cycle
                         if (Input%Ref_063um == Missing_Value_Real4) cycle
                         Classifier_Value(Class_Idx) = relative_visible_contrast_test( &
                                       Input%Ref_063um_Min, Input%Ref_063um)
                       else
                         if (Input%Chan_On_DNB == symbol%NO) cycle
                         if (Lunar_Spatial_Flag == symbol%NO) cycle
                         if (Night_Lunar_Flag == symbol%NO) cycle
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Coastal_Flag == symbol%YES) cycle
                         if (City_Flag == symbol%YES) cycle  
                         if (Input%Ref_Lunar == Missing_Value_Real4) cycle
                         if (Input%Ref_Lunar_Min == Missing_Value_Real4) cycle
                         if (Input%Ref_Lunar <= 0.0) cycle
                         if (Input%Ref_Lunar_Min <= 0.0) cycle
                         Classifier_Value(Class_Idx) = relative_visible_contrast_test( &
                                       Input%Ref_Lunar_Min, Input%Ref_Lunar)
                       endif

                    case("Ref_Ratio_Day")
                       if (Input%Chan_On_063um == symbol%NO) cycle
                       if (Input%Chan_On_086um == symbol%NO) cycle
                       if (Day_063_Flag == symbol%NO) cycle
!                       if (Forward_Scattering_Flag == symbol%YES) cycle ! ---> DENIS TEST
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       if (Input%Ref_063um == Missing_Value_Real4) cycle
                       if (Input%Ref_086um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = reflectance_ratio_test( &
                                         Input%Ref_063um,Input%Ref_086um)

                    case("Ref_138_Day")
                       if (Input%Chan_On_138um == symbol%NO) cycle
                       if (Forward_Scattering_Flag == symbol%YES) cycle
                       if (Day_063_Flag == symbol%NO) cycle
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Input%Ref_138um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = Input%Ref_138um

                    case("Ref_160_Day")
                       if (Input%Chan_On_160um == symbol%NO) cycle
                       if (Forward_Scattering_Flag == symbol%YES) cycle
                       if (Day_063_Flag == symbol%NO) cycle
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Input%Ref_160um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = Input%Ref_160um

                    case("Ref_375_Day")
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Forward_Scattering_Flag == symbol%YES) cycle
                       if (Day_063_Flag == symbol%NO) cycle
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Input%Ref_375um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = Input%Ref_375um

                    case("Ndsi_Day")
                       if (Input%Chan_On_063um == symbol%NO) cycle
                       if (Input%Chan_On_160um == symbol%NO) cycle
                       if (Forward_Scattering_Flag == symbol%YES) cycle
                       if (Day_063_Flag == symbol%NO) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       if (Input%Ref_160um == Missing_Value_Real4) cycle
                       if (Input%Ref_063um == Missing_Value_Real4) cycle
                       Classifier_Value(Class_Idx) = (Input%Ref_063um - Input%Ref_160um) /  &
                                                     (Input%Ref_063um + Input%Ref_160um)
                    
                     case default
                       print *, "Unknown Classifier Naive Bayesian Cloud Mask, returning"
                       print *, "Name = ",trim(Classifier_Value_Name(Class_Idx,1))
                       return
             end select

             !--- Turn off Classifiers if Chosen Metric is Missing
             if (Classifier_Value(Class_Idx) == Missing_Value_Real4) then
                cycle
             endif


             !--- interpolate class conditional values
             Bin_Idx = int ((Classifier_Value(Class_Idx) - Classifier_Bounds_Min(Class_Idx,Sfc_Idx))  /    &
                           Delta_Classifier_Bounds(Class_Idx,Sfc_Idx)) + 1
             Bin_Idx = max(1,min(N_bounds-1,Bin_Idx))

             Cond_Ratio(Class_Idx) = Class_Cond_Ratio(Bin_Idx,Class_Idx,Sfc_Idx)

!print *, "Cond Ratio T_Std = ", Class_Cond_Ratio(Bin_Idx,:,Sfc_Idx)
!stop

!            if (Class_Idx == 7) then
!               print *, Classifier_Value(Class_Idx), Classifier_Bounds_Min(Class_Idx,Sfc_Idx), Delta_Classifier_Bounds(Class_Idx,Sfc_Idx), Bin_Idx,  Class_Cond_Ratio(Bin_Idx,Class_Idx,Sfc_Idx)
!            endif

        enddo  class_loop 


        !------------------------------------------------------------------------------------------------------------
        !--- compute prosterior probabilites for each pixel
        !-----------------------------------------------------------------------------------------------------------

        r = product(Cond_Ratio)
        Output%Posterior_Cld_Probability =  1.0 / (1.0 + r/Prior_Yes_Temp - r)

        !--- constrain 
        if (r < 0.001) Output%Posterior_Cld_Probability = 1.0
        if (r > 99.0) Output%Posterior_Cld_Probability = 0.0

        !--- check for a missing prior
        if (Prior_Yes_Temp == Missing_Value_Real4) then 
          Output%Posterior_Cld_Probability = Missing_Value_Real4
          Posterior_Cld_Probability_By_Class = Missing_Value_Real4
        endif

        !------------------------------------------------------------------------------------------------------------
        !--- make a cloud mask
        !------------------------------------------------------------------------------------------------------------
        ! - based on type of srfc could be different thresholds
        if (Sfc_Idx > 0) then
           Output%Cld_Mask_Bayes = symbol%CLEAR
           if (Output%Posterior_Cld_Probability >= PROB_CLOUDY_CONF_CLOUD_THRESH(Sfc_Idx)) then
                   Output%Cld_Mask_Bayes = symbol%CLOUDY
           endif
           if ((Output%Posterior_Cld_Probability >= PROB_CLEAR_PROB_CLOUD_THRESH(Sfc_Idx)) .and. &
               (Output%Posterior_Cld_Probability < PROB_CLOUDY_CONF_CLOUD_THRESH(Sfc_Idx))) then
                   Output%Cld_Mask_Bayes = symbol%PROB_CLOUDY
           endif
           if ((Output%Posterior_Cld_Probability > CONF_CLEAR_PROB_CLEAR_THRESH(Sfc_Idx)) .and. &
               (Output%Posterior_Cld_Probability < PROB_CLEAR_PROB_CLOUD_THRESH(Sfc_Idx))) then
                   Output%Cld_Mask_Bayes = symbol%PROB_CLEAR
           endif
        endif

        !------------------------------------------------------------------------------------------------------------
        !--- compute probabilities for each class alone - used for flags - not
        !--- needed for mask or final probability - it should remain optional
        !------------------------------------------------------------------------------------------------------------
        if (Do_By_Class_Flag == symbol%YES) then

         do Class_Idx = 1, N_class
 
          r = Cond_Ratio(Class_Idx)
          Posterior_Cld_Probability_By_Class(Class_Idx) =  1.0 / (1.0 + r/Prior_Yes(Sfc_Idx) - r)

          !-- set cloud flags
          Cld_Flag_Bit_Depth(Class_To_Test_Idx(Class_Idx)) = 2
          Cld_Flags(Class_To_Test_Idx(Class_Idx)) = symbol%CLOUDY
          if (Posterior_Cld_Probability_By_Class(Class_Idx) <= CONF_CLEAR_PROB_CLEAR_THRESH(Sfc_Idx)) then
               Cld_Flags(Class_to_Test_Idx(Class_Idx)) = symbol%CLEAR

          elseif (Posterior_Cld_Probability_By_Class(Class_Idx) > CONF_CLEAR_PROB_CLEAR_THRESH(Sfc_Idx) .and. &
                  Posterior_Cld_Probability_By_Class(Class_Idx) <= PROB_CLEAR_PROB_CLOUD_THRESH(Sfc_Idx)) then
               Cld_Flags(Class_to_Test_Idx(Class_Idx)) = symbol%PROB_CLEAR

          elseif (Posterior_Cld_Probability_By_Class(Class_Idx) > PROB_CLEAR_PROB_CLOUD_THRESH(Sfc_Idx) .and. &
                  Posterior_Cld_Probability_By_Class(Class_Idx) <= PROB_CLOUDY_CONF_CLOUD_THRESH(Sfc_Idx)) then
               Cld_Flags(Class_to_Test_Idx(Class_Idx)) = symbol%PROB_CLOUDY
          endif


         !-----------------------------------------------------------------------
         !--- Diagnostic Output for debugging in CLAVR-x
         !-----------------------------------------------------------------------
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Ref_std") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Ref_063_Day") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Ref_063_Min_3x3_Day") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Ndsi_Day") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Btd_11_67") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Btd_375_11_Night") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "T_11") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Emiss_375") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Emiss_375_Day") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Emiss_375_Night") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Bt_11_67_Covar") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Emiss_Tropo") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "T_std") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "T_Max-T") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "FMFT") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Btd_375_11_Day") then
         !     if (present(Diag)) Diag%Array_1 = Classifier_Value(Class_Idx)
         !     if (present(Diag)) Diag%Array_2 = Posterior_Cld_Probability_By_Class(Class_Idx)
         !     if (present(Diag)) Diag%Array_3 = Cld_Flags(Class_To_Test_Idx(Class_Idx))
         ! endif

!        if (trim(Classifier_Value_Name(Class_Idx,1)) == "Btd_375_11_Night") Diag%Array_1 = Posterior_Cld_Probability_By_Class(Class_Idx)
!        if (trim(Classifier_Value_Name(Class_Idx,1)) == "Emiss_Tropo") Diag%Array_2 = Posterior_Cld_Probability_By_Class(Class_Idx)
!        if (trim(Classifier_Value_Name(Class_Idx,1)) == "FMFT") Diag%Array_3 = Posterior_Cld_Probability_By_Class(Class_Idx)

         enddo

         !-------------------------------------------------------------------
         ! Pack Bytes
         !-------------------------------------------------------------------
         call PACK_BITS_INTO_BYTES(Cld_Flags,Cld_Flag_Bit_Depth,Output%Cld_Flags_Packed(:))

        endif   !Do_By_Class_Loop
      endif   !If Skip_Sfc_Type_Flag
   endif   !If Invalid_Data_Mask

   !----- deallocate memory
   if (allocated(Cond_Ratio)) deallocate(Cond_Ratio)
   if (allocated(Classifier_Value)) deallocate(Classifier_Value)
   if (allocated(Posterior_Cld_Probability_By_Class)) deallocate(Posterior_Cld_Probability_By_Class)

 end subroutine NB_CLOUD_MASK_ALGORITHM

!====================================================================
! SUBROUTINE Name: COMPUTE_BAYES_SFC_TYPE
!
! Function:
!   Computes the bayesian surface type given the ancillary sfc data
!
!====================================================================
 function COMPUTE_BAYES_SFC_TYPE(Land_Temp, Coast_Temp, Snow_Temp, Sfc_Type_Temp, &
                                 Lat_Temp, Lon_Temp, Sst_Back_Uni_Temp,&
                                 Emiss_Sfc_375um_Temp, symbol) &
                                 result(Bayes_Mask_Sfc_Type_Temp)
   integer(kind=int1), intent(in):: Land_Temp
   integer(kind=int1), intent(in):: Coast_Temp
   integer(kind=int1), intent(in):: Snow_Temp
   integer(kind=int1), intent(in):: Sfc_Type_Temp
   real(kind=real4), intent(in):: Lat_Temp
   real(kind=real4), intent(in):: Lon_Temp
   real(kind=real4), intent(in):: Sst_Back_Uni_Temp
   real(kind=real4), intent(in):: Emiss_Sfc_375um_Temp
   TYPE(symbol_naive_bayesian), intent(in) :: symbol

   integer(kind=int4):: Bayes_Mask_Sfc_Type_Temp


   if (Land_Temp == symbol%LAND) then
           Bayes_Mask_Sfc_Type_Temp = 3
   else
           Bayes_Mask_Sfc_Type_Temp = 1
   endif

   !--- #2 - Shallow Ocean
   if ((Land_Temp == symbol%MODERATE_OCEAN) .or. &
       (Land_Temp == symbol%DEEP_INLAND_WATER) .or. &
       (Land_Temp == symbol%SHALLOW_INLAND_WATER) .or. &
       (Land_Temp == symbol%SHALLOW_OCEAN)) then
           Bayes_Mask_Sfc_Type_Temp = 2
   endif
   if ((Land_Temp /= symbol%LAND) .and. &
       (Sst_Back_Uni_Temp > 0.5)) then
           Bayes_Mask_Sfc_Type_Temp = 2
   endif

   !--- #3 Unfrozen_Land 
   if ((Land_Temp == symbol%LAND) .or. &
       (Land_Temp == symbol%COASTLINE) .or. &
       (Coast_Temp == symbol%YES) .or. &
       (Land_Temp == symbol%EPHEMERAL_WATER)) then

           Bayes_Mask_Sfc_Type_Temp = 3

   endif

   !--- #4 - Snow Covered Land
   if ((Lat_Temp > -60.0) .and. (Snow_Temp == symbol%SNOW)) then
           Bayes_Mask_Sfc_Type_Temp = 4
   endif

   !--- #5 - Arctic
   if ((Lat_Temp >= 0.0) .and. (Snow_Temp == symbol%SEA_ICE)) then
           Bayes_Mask_Sfc_Type_Temp = 5
   endif

   !--- #6 - Antarctic
   if ((Lat_Temp <= -60.0) .and. (Snow_Temp == symbol%SNOW)) then
           Bayes_Mask_Sfc_Type_Temp = 6
   endif
   if ((Lat_Temp <= 0.0) .and. (Snow_Temp == symbol%SEA_ICE)) then
           Bayes_Mask_Sfc_Type_Temp = 6
   endif
   if ((Lat_Temp >= 60.0) .and.  &
       (Lon_Temp > -75.0) .and. (Lon_Temp < -10.0) .and. &
       (Land_Temp == symbol%LAND .or. Land_Temp == symbol%COASTLINE) .and. &
       (Snow_Temp == symbol%SNOW)) then
           Bayes_Mask_Sfc_Type_Temp = 6
   endif

   !--- #7 - Desert
   if ( (Emiss_Sfc_375um_Temp < 0.90 ) .and. (abs(Lat_Temp) < 60.0) .and. &
       ((Sfc_Type_Temp == symbol%OPEN_SHRUBS_SFC) .or. (Sfc_Type_Temp == symbol%BARE_SFC)) ) then
           Bayes_Mask_Sfc_Type_Temp = 7
   endif
 
 
   if ( Bayes_Mask_Sfc_Type_Temp == 3 .and.  &
        Emiss_Sfc_375um_Temp < 0.93  .and.  &
        abs(Lat_Temp) < 60.0 .and. &
       ((Sfc_Type_Temp == symbol%OPEN_SHRUBS_SFC) .or.  &
        (Sfc_Type_Temp == symbol%CLOSED_SHRUBS_SFC) .or. &
        (Sfc_Type_Temp == symbol%GRASSES_SFC) .or.  &
        (Sfc_Type_Temp == symbol%BARE_SFC)) ) then
           Bayes_Mask_Sfc_Type_Temp = 7
   endif

  return 

 end function COMPUTE_BAYES_SFC_TYPE

!-----------------------------------------------------------------------------
! EUMETCAST Fire detection algorithm
!-----------------------------------------------------------------------------
  integer elemental function FIRE_TEST ( t11, t375, t11_std, t375_std, solzen)

     real, intent(in):: T11
     real, intent(in):: T375
     real, intent(in):: T11_Std
     real, intent(in):: T375_Std
     real, intent(in):: Solzen
     
     real :: Bt_375um_Eumet_Fire_Thresh
     real :: Bt_Diff_Eumet_Fire_Thresh
     real :: Stddev_11um_Eumet_Fire_Thresh 
     real :: Stddev_375um_Eumet_Fire_Thresh

     Fire_Test = 0

     if (T375_std /= Missing_Value_Real4 .and. &
         T11_std /= Missing_Value_Real4) then
         
         !Day
         if (Solzen < EumetCAST_Fire_Day_Solzen_Thresh) then 
            Bt_375um_Eumet_Fire_Thresh = Bt_375um_Eumet_Fire_day_Thresh
            Bt_Diff_Eumet_Fire_Thresh = Bt_Diff_Eumet_Fire_day_Thresh
            Stddev_11um_Eumet_Fire_Thresh = Stddev_11um_Eumet_Fire_Day_Thresh
            Stddev_375um_Eumet_Fire_Thresh = Stddev_375um_Eumet_Fire_Day_Thresh
         endif
         
         !Night
         if (solzen > EumetCAST_Fire_Night_Solzen_Thresh) then 
            Bt_375um_Eumet_Fire_Thresh = Bt_375um_Eumet_Fire_Night_Thresh
            Bt_Diff_Eumet_Fire_Thresh = Bt_Diff_Eumet_Fire_Night_Thresh
            Stddev_11um_Eumet_Fire_Thresh = Stddev_11um_Eumet_Fire_Night_Thresh
            Stddev_375um_Eumet_Fire_Thresh = Stddev_375um_Eumet_Fire_Night_Thresh
         endif
         
         if ((Solzen >= EumetCAST_Fire_Day_Solzen_Thresh) .and. &
             (Solzen <= EumetCAST_Fire_Night_Solzen_Thresh)) then
             
             !linear fit day -> night
             Bt_375um_Eumet_Fire_Thresh = ((-1.0)* Solzen) + 380.0
             Bt_Diff_Eumet_Fire_Thresh = ((-0.4)* Solzen) + 36.0
             
             !These two don't change, but 
             Stddev_11um_Eumet_Fire_Thresh = ((0.0)* solzen) + 1.0
             Stddev_375um_Eumet_Fire_Thresh = ((0.0)* solzen) + 4.0
         
         endif

       ! All of these conditions need to be met
       if ((T375 > bt_375um_Eumet_Fire_Thresh) .and. &
           ((T375 - t11) > Bt_Diff_Eumet_Fire_Thresh) .and. &
           (T375_Std > Stddev_375um_Eumet_Fire_Thresh) .and. &
           (T11_Std < Stddev_11um_Eumet_Fire_Thresh)) then
         Fire_Test = 1
       endif
         
     endif

  end function FIRE_TEST

  !-------------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------------
  real elemental function SPLIT_WINDOW_TEST ( t11_clear, t12_clear, t11, t12)

     real, intent(in):: t11_clear
     real, intent(in):: t12_clear 
     real, intent(in):: t11
     real, intent(in):: t12

     split_window_test  = (t11_clear - t12_clear) * (t11 - 260.0) / (t11_clear - 260.0)

     if (t11_clear <=265.0) split_window_test = 0.0

     split_window_test = (t11 - t12) - split_window_test

  end function SPLIT_WINDOW_TEST

  !-------------------------------------------------------------------------------------
  real elemental function REFLECTANCE_GROSS_CONTRAST_TEST (ref_clear,ref)
     real, intent(in):: ref_clear
     real, intent(in):: ref

     reflectance_gross_contrast_test = Missing_Value_Real4
     if (ref_clear /= Missing_Value_Real4) then
           reflectance_gross_contrast_test = ref - ref_clear
     endif

  end function REFLECTANCE_GROSS_CONTRAST_TEST

  !-------------------------------------------------------------------------------------
  real elemental function RELATIVE_VISIBLE_CONTRAST_TEST (ref_min,ref)
     real, intent(in):: ref_min
     real, intent(in):: ref

     relative_visible_contrast_test = Missing_Value_Real4
     if (ref_min /= Missing_Value_Real4) then
           relative_visible_contrast_test = ref - ref_min
     endif

  end function RELATIVE_VISIBLE_CONTRAST_TEST

  !-------------------------------------------------------------------------------------
  real elemental function REFLECTANCE_RATIO_TEST (ref_vis,ref_nir)
     real, intent(in):: ref_vis
     real, intent(in):: ref_nir

     reflectance_ratio_test = Missing_Value_Real4
     if (ref_vis /= Missing_Value_Real4 .and. ref_nir /= Missing_Value_Real4) then
           reflectance_ratio_test = ref_nir / ref_vis
     endif

  end function REFLECTANCE_RATIO_TEST

  !-------------------------------------------------------------------------------------
  real elemental function NIR_REFLECTANCE_GROSS_CONTRAST_TEST (ref_clear,ref)
     real, intent(in):: ref_clear
     real, intent(in):: ref

     nir_reflectance_gross_contrast_test = Missing_Value_Real4
     if (ref /= Missing_Value_Real4 .and. ref_clear /= Missing_Value_Real4) then
           nir_reflectance_gross_contrast_test = ref- ref_clear
     endif

  end function NIR_REFLECTANCE_GROSS_CONTRAST_TEST

  !-------------------------------------------------------------------------------------
  real elemental function EMISS_375UM_TEST (ems,ems_clear)
     real, intent(in):: ems_clear
     real, intent(in):: ems

     emiss_375um_test = Missing_Value_Real4
     if (ems /= Missing_Value_Real4 .and. ems_clear /= Missing_Value_Real4) then
           emiss_375um_test = (ems- ems_clear) / ems_clear
     endif

  end function EMISS_375UM_TEST

  !-------------------------------------------------------------------------------------
  real elemental function EMISS_375UM_DAY_TEST (ems,ems_clear)
     real, intent(in):: ems_clear
     real, intent(in):: ems

     emiss_375um_day_test = Missing_Value_Real4
     if (ems /= Missing_Value_Real4 .and. ems_clear /= Missing_Value_Real4) then
           emiss_375um_day_test = (ems- ems_clear) / ems_clear
     endif

  end function EMISS_375UM_DAY_TEST

  !-------------------------------------------------------------------------------------
  real elemental function EMISS_375UM_NIGHT_TEST (ems,ems_clear)
     real, intent(in):: ems_clear
     real, intent(in):: ems

     emiss_375um_night_test = Missing_Value_Real4
     if (ems /= Missing_Value_Real4 .and. ems_clear /= Missing_Value_Real4) then
           emiss_375um_night_test = ems
         ! emiss_375um_night_test = (ems- ems_clear) / ems_clear
     endif

  end function EMISS_375UM_NIGHT_TEST

!------------------------------------------------------------------------
! subroutine PACK_BITS_INTO_BYTES(input_bits,bit_depth,output_bytes)
!
! Routines to pack individual bytes into a single byte
!
! input:
! input_bits - vector of bytes to be packed into output_byte
! bit_start - vector of bit starting positions (1-7) for each input byte
! bit_depth - vector of bit depths (1-7) for each input byte (total can not exceed 8)
!
! output: 
! output_byte - byte variable that holds the bit values of the input_bits - can be i1 or i2
!
! local
!  n_in - number of elements of input vectors
!  i_in - index of input_bits (1-n_in)
!  n_out - number of elements of output vectors
!  i_out - index of output_bytes (1-n_out)
!
! Note:
! 1.  if the input byte has information in bits greater then bit depth - they are removed 
!
!
! Example, pack an input_byte wth  bit_start = 2 and bit depth 3 
!
! input byte
!           x x x
! _ _ _ _ _ _ _ _
!
! result of first ishft
! x x x
! _ _ _ _ _ _ _ _
!
! result of second ishft
!       x x x
! _ _ _ _ _ _ _ _

!
! Author: Andrew Heidinger
!
! Version History:  
! February 2006 - Created
!-----------------------------------------------------------------------------------

!--- This Version packs into one byte words
   subroutine PACK_BITS_INTO_BYTES (input_bits,bit_depth,output_bytes)
    integer(kind=int1), dimension(:), intent(in):: input_bits
    integer(kind=int1), dimension(:), intent(in):: bit_depth
    integer(kind=int1), dimension(:), intent(out):: output_bytes
    integer(kind=int1):: bit_start, bit_end, bit_offset
    integer(kind=int1):: temp_byte
    integer:: n_in,i_in,n_out,i_out
    integer, parameter:: word_bit_depth = 8

    !--- determine size of vectors
    n_in = size(input_bits)
    n_out = size(output_bytes)

    !--- reset output byte
    output_bytes = 0

    !--- initialize
    bit_offset = 0
    bit_start = 0
    bit_end = 0
    i_out = 1

    !--- loop through input bytes
    do i_in = 1, n_in

     !--- determine starting and ending bit locations
     bit_start = bit_offset + 1
     bit_end = bit_start + bit_depth(i_in) - 1

     !--- determine if this input byte will fit on current output byte, if not go to next
     if (bit_end > word_bit_depth) then
      i_out = i_out + 1
      bit_offset = 0
      bit_start = bit_offset + 1
      bit_end = bit_start + bit_depth(i_in) - 1
     endif

     !--- check for exceeding the space allowed for the packed bytes
     if (i_out > n_out) then
       print *, "ERROR: Insufficient space for bit packing ", i_out, n_out
       return
     endif

     !--- place input byte into correct position
     temp_byte =0
     temp_byte = ishft(input_bits(i_in),word_bit_depth-bit_depth(i_in))    !first ishft
     temp_byte = ishft(temp_byte,bit_end - word_bit_depth)                 !second ishft

     !--- modify output byte
     output_bytes(i_out) = output_bytes(i_out) + temp_byte

     !--- update bit offset
     bit_offset = bit_offset + bit_depth(i_in)

    enddo

  end subroutine  PACK_BITS_INTO_BYTES

!-----------------------------------------------------------------------------------

end module NB_CLOUD_MASK
