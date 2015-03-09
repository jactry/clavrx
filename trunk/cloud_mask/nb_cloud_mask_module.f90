! $Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/trunk/cloud_mask/naive_bayesian_cloud_mask_module.f90 629 2014-10-31 17:28:38Z awalther $
!
!----------------------------------------------------------------------
! MODULE name: NAIVE_BAYESIAN_CLOUD_MASK
! 
! Routines for the determination of the Naive Bayesian cloud mask
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! History:
!      2014/04/06:  new interface (AW)
!      2014/05/07:  added diagnostic type (Denis B)
!      2014/05/09:  added version (Denis B)
!
! DEPENDENCIES: None
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
! 12           1            12        2    blank
! 13           1            13        2    smoke detected
! 14           1            14        2    dust  detected
! 15           1            15        2    shadow detected
! 16           1            16        2    fire detected
!---
! 17           3            17-19     3    nbvcm surface type
! 18           1            20        3    blank   
!<-------------------- START OF CLOUD TESTS -------------------------->
! 19           2            21-22     3    T_11           (TGCT)
! 20           2            23-24     3    T_max-T        (RTCT)
! ---
! 21           2            25-26     4    T_Std          (TUT)
! 22           2            27-28     4    Emiss_Tropo
! 23           2            29-30     4    FMFT mask (Btd_11_12)
! 24           2            31-32     4    Btd_11_67 
!---
! 25           2            33-34     5    Bt_11_67_Covar 
! 26           2            35-36     5    Btd_11_85_Covar
! 27           2            37-38     5    Emiss_375
! 28           2            39-40     5    Emiss_375_Day
!---
! 29           2            41-42     6    Emiss_375_Night
! 30           2            43-44     6    Btd_375_11_Night
! 31           2            45-46     6    Ref_063_Day (RGCT)
! 32           2            47-48     6    Ref_std      (RUT)
!---
! 33           2            49-50     7    Ref_063_Min_3x3 (RVCT)
! 34           2            51-52     7    Ref_Ratio_Day   (RRCT)
! 35           2            53-54     7    Ref_138_Day     (CIRREF)
! 36           2            55-56     7    Ndsi_Day     (NIRREF)
!
!----------------------------------------------------------------------

module NB_CLOUD_MASK

 use NB_CLOUD_MASK_SERVICES
 use NETCDF

   implicit none
   
   type et_land_class_type
      integer :: FIRST = 0
      integer :: SHALLOW_OCEAN = 0
      integer :: LAND = 1
      integer :: COASTLINE = 2
      integer :: SHALLOW_INLAND_WATER = 3
      integer :: EPHEMERAL_WATER = 4
      integer :: DEEP_INLAND_WATER = 5
      integer :: MODERATE_OCEAN = 6
      integer :: DEEP_OCEAN = 7 
      integer :: LAST = 7     
   end type
   
   type ( et_land_class_type) , public :: ET_land_class 
   
   type et_snow_class_type
      integer :: FIRST = 1
      integer :: NO_SNOW = 1
      integer :: NO_SNOW_NOR_ICE = 1
      integer :: SEA_ICE = 2
      integer :: SNOW = 3
      integer :: LAST = 3   
   end type
   type ( et_snow_class_type) , public :: ET_snow_class 
   
   type et_cloudiness_class_type
      integer :: SPACE = 10
      integer :: MISSING = -128
      integer :: CLOUDY = 3
      integer :: PROB_CLOUDY = 2
      integer :: PROB_CLEAR = 1
      integer :: CLEAR = 0
   end type

   type ( et_cloudiness_class_type) , public :: ET_cloudiness_class 
      integer , parameter :: et_class_FIRST = 1
      integer , parameter :: et_class_T110 = 1
      integer , parameter :: et_class_TMAX_T = 2
      integer , parameter :: et_class_T_STD = 3
      integer , parameter :: et_class_E_TROP = 4
      integer , parameter :: et_class_FMFT = 5
      integer , parameter :: et_class_BTD_110_067 = 6
      integer , parameter :: et_class_BTD_110_067_COV = 7
      integer , parameter :: et_class_BTD_110_085 = 8
      integer , parameter :: et_class_E_037 = 9
      integer , parameter :: et_class_E_037_DAY = 10
      integer , parameter :: et_class_E_037_NGT = 11
      integer , parameter :: et_class_BTD_037_110_NGT = 12
      integer , parameter :: et_class_R_006_DAY = 13
      integer , parameter :: et_class_R_006_STD = 14
      integer , parameter :: et_class_R_006_MIN_3x3_DAY = 15
      integer , parameter :: et_class_R_RATIO_DAY = 16
      integer , parameter :: et_class_R_013_DAY = 17
      integer , parameter :: et_class_NDSI_DAY = 18
      integer , parameter :: et_class_LAST = 18
      
   type cloud_mask_geo_type
      real :: lat
      real :: lon
      real :: sol_zen
      real :: sat_zen
      real :: Lunar_zen
      real :: airmass
      real :: scat_angle
      real :: scat_angle_lunar
      logical :: lunar_glint_mask
      logical :: glint
      logical :: solar_conta
   end type cloud_mask_geo_type
   
   type cloud_mask_sfc_type
      integer :: land_class
      logical :: coast_mask
      integer :: snow_class
      integer :: sfc_type
      real :: dem
      real :: emis_ch20
      real :: sst_anal_uni
      logical :: is_city
   end type cloud_mask_sfc_type
   
   type cloud_mask_rtm_type 
      real :: bt_ch31_lrc
      real :: bt_ch31_3x3_max
      real :: bt_ch31_3x3_std
      real :: bt_ch20_3x3_std
      real :: emis_ch31_tropo
      real :: emis_ch32_tropo
      real :: emis_ch20_clear
      real :: bt_ch31_atm_sfc
      real :: bt_ch32_atm_sfc
      real :: bt_ch31_ch27_covar
      real :: ref_ch1_clear
      real :: ref_dnb_clear
   end type cloud_mask_rtm_type
   
   
   type cloud_mask_sat_viirs_iband_stats_type
      logical :: is_set
      real :: min
      real :: max
      real :: mean
      real :: std   
   end type cloud_mask_sat_viirs_iband_stats_type
   
   type cloud_mask_sat_viirs_iband_type
      logical :: is_set
      type ( cloud_mask_sat_viirs_iband_stats_type ) :: bt
      type ( cloud_mask_sat_viirs_iband_stats_type ) :: ref
   
   end type cloud_mask_sat_viirs_iband_type
   
   
   type cloud_mask_sat_type
      logical , dimension(44) :: chan_on
      real :: bt_ch20
      real :: bt_ch27
      real :: bt_ch29
      real :: bt_ch31
      real :: bt_ch32
      real :: ref_ch1
      real :: ref_ch2
      real :: ref_ch6
      real :: ref_ch7
      real :: ref_ch8
      real :: ref_ch26
      real :: emis_ch20_3x3_mean
      real :: ref_ch1_3x3_std
      real :: ref_ch1_3x3_min
      real :: ref_dnb_3x3_std
      real :: ref_dnb_3x3_min
      real :: ref_dnb_lunar
      type ( cloud_mask_sat_viirs_iband_type ) :: iband(5)
   end type cloud_mask_sat_type

   type cloud_mask_diagnostic
      real :: diagnostic_1
      real :: diagnostic_2
      real :: diagnostic_3
   end type cloud_mask_diagnostic

   type cloud_mask_version_type
      character (len = 355) :: cloud_mask_thresh_version_id
      character (len = 355) :: cloud_mask_version_id
   end type cloud_mask_version_type

   type cloud_mask_input_type
      character (len = 355) :: bayesian_mask_classifier      
      type (cloud_mask_geo_type) :: geo
      type (cloud_mask_sfc_type) :: sfc
      type (cloud_mask_rtm_type) :: rtm
      type (cloud_mask_sat_type) :: sat
   end type cloud_mask_input_type
         
   type bayes_coef_type
      character (len =120) :: cvs_version
      character (len =1000) :: file
      logical :: is_read = .false.
      integer :: n_class
      integer :: n_bounds
      integer :: N_sfc_bayes
      real, allocatable :: prior_no (:)
      real, allocatable :: prior_yes (:)
      real, allocatable :: Optimal_Posterior_Prob (:)
      real, allocatable :: Classifier_Bounds_Min(:,:)
      real, allocatable :: Classifier_Bounds_Max(:,:)
      real, allocatable :: Delta_Classifier_Bounds(:,:)
      integer, allocatable :: First_Valid_Classifier_Bounds(:,:)
      integer, allocatable :: Last_Valid_Classifier_Bounds(:,:)
      real, allocatable :: class_cond_no(:,:,:)
      real, allocatable :: class_cond_yes(:,:,:) 
      real, allocatable :: Class_Cond_Ratio(:,:,:) 
      real, allocatable :: cond_no(:,:,:)
      real, allocatable :: cond_yes(:,:,:)
      real, allocatable :: Cond_Ratio(:,:,:)
      logical, allocatable :: do_this_classifier (:,:)
      character (len=30), allocatable :: Classifier_Value_Name(:,:)
      integer, allocatable :: Classifier_Value_Name_enum(:)
      integer, allocatable :: flag_idx(:)
   
   contains
      procedure :: alloc => alloc_bayes_coef
      procedure :: dealloc => dealloc_bayes_coef
   end type bayes_coef_type
   
   type ( bayes_coef_type) , private , save :: bayes_coef
  
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
   !
   !
   !
   subroutine NB_CLOUD_MASK_ALGORITHM ( inp , erg , info_flags , diag , vers )
          
      implicit none
            
      type ( cloud_mask_input_type ) , intent ( in ) :: inp
      type ( cloud_mask_diagnostic ) , intent ( inout ) :: diag
      type ( cloud_mask_version_type ) , intent ( out ) :: vers
      real , intent ( out ) :: erg
      integer , intent ( out ) , optional :: info_flags ( 7 )
      
      integer :: sfc_type_number 
      integer :: class_idx, sfc_idx
      real :: Classifier_Value
      real, allocatable :: Cond_Yes(:)
      real, allocatable :: Cond_No(:)
      real, allocatable :: Cond_Ratio(:)
      
      integer :: bin_idx
      
      logical :: is_on_test
      logical :: is_mountain 
      logical :: is_day_375um
      logical :: is_day_063um
      logical :: is_day_063um_spatial_tests
      logical :: is_night_375um
      logical :: has_cold_btd
      logical :: is_cold_375um
      logical :: is_forward_scatter
      logical :: is_dust
      logical :: is_smoke
      logical :: is_cloud_shadow
      logical :: is_fire
      logical :: is_solar_contaminated
      logical :: use_lunar_refl_for_vis_tests
  
      real, parameter :: SOLZEN_DAY_THRESH = 85.0       !was 85.0
      real, parameter :: AIRMASS_THRESH = 5.0
      real, parameter :: SOLZEN_375UM_NIGHT_THRESH = 90.0
      real, parameter :: SOLZEN_375UM_DAY_THRESH = 85.0
      real, parameter :: SOLZEN_063UM_DAY_THRESH = 80.0
      real, parameter :: SOLZEN_063UM_DAY_THRESH_SPATIAL_TESTS = 85.0
      real, parameter :: BT_11UM_COLD_SCENE_THRESH = 220.0
      real, parameter :: BT_037UM_COLD_SCENE_THRESH = 240.0
      real, parameter :: MISSING_VALUE_REAL4 = -999.
      real, parameter :: Radiance_Lunar_City_Thresh = 2.5e-08
      real, parameter :: Scat_Angle_Lunar_Thresh = 80.0
      real, parameter :: Lunar_Zen_Thresh = 95.0
            
      integer :: pos_info_flag 
      integer :: idx_info_flag
      real :: class_contr
      real :: r
      
      ! ---    Executable  ----------------------------
      ! - read in classifer coeffiecients
      bayes_coef %  file =trim(inp  %  bayesian_mask_classifier ) 
      !if ( .not. bayes_coef % is_read) call READ_BAYES_COEFF ( ) 
      if ( .not. bayes_coef % is_read) call READ_BAYES_COEFF_NC ( ) 
      
      ! - set mask and thresholds version id
      vers % cloud_mask_thresh_version_id = bayes_coef % cvs_version 
      vers % cloud_mask_version_id = "$Id: naive_bayesian_cloud_mask_module.f90 1068 2015-03-06 23:31:56z dbotambekov $"

      ! - determine sfc type
      sfc_type_number =  BAYES_SFC_TYPE ( inp% geo % lat , inp % geo % lon &
          , inp % sfc % land_class , inp % sfc % coast_mask, inp % sfc % snow_class , inp % sfc % sfc_type &
          , inp % sfc % emis_ch20,  inp % sfc % sst_anal_uni )
          
         
      allocate ( Cond_Yes ( Bayes_Coef % N_Class ) )
      allocate ( Cond_No ( Bayes_Coef % N_Class ) )
      allocate ( Cond_Ratio ( Bayes_Coef % N_Class ) )

      Cond_Yes = 1.0
      Cond_No =  1.0
      Cond_Ratio =  1.0
          
      sfc_idx = sfc_type_number
           
      ! - several 0/1 flags
      is_mountain =  inp % sfc % dem  > 2000.0 &
                           .and. sfc_idx /= 6
                            
      use_lunar_refl_for_vis_tests = .false.
      if ( inp % sat % chan_on(44) ) then

         if ( inp % sat % ref_dnb_lunar >= 0. .and. &
             inp % geo %  scat_angle_lunar > Scat_Angle_Lunar_Thresh .and. &
             inp % geo % lunar_zen < Lunar_Zen_Thresh .and. &
            .not. inp % sfc % is_city  .and. &
            .not. is_mountain .and. &
            .not. inp % sfc % coast_mask .and. &
            .not. inp % sfc % snow_class == ET_snow_class % SNOW .and. &
            .not. inp % geo % lunar_glint_mask  )  then
             use_lunar_refl_for_vis_tests  = .true.      
         end if    
      end if                       
            
      has_cold_btd = .false.
      if ( inp % sat % chan_on(31) ) then
         if (  inp%sat %bt_ch31  < BT_11UM_COLD_SCENE_THRESH .and. &
                     inp%sat %bt_ch31 /= MISSING_VALUE_REAL4 ) then 
                  has_cold_btd = .true. 
         end if         
      end if
            
      is_cold_375um = .false.
      if ( inp % sat % chan_on(20) ) then
         if (  inp%sat %bt_ch20 < BT_037UM_COLD_SCENE_THRESH .and. &
                     inp%sat %bt_ch20 /= MISSING_VALUE_REAL4 ) then 
            is_cold_375um = .true. 
         end if         
      end if
            
      is_day_375um = .true.
      if (( inp % geo % sol_zen  >  SOLZEN_375UM_DAY_THRESH ) .or. &
                    (inp % geo %airmass > AIRMASS_THRESH  )) then           
         is_day_375um = .false.
      end if   
            
      is_night_375um = .true.
      if ( inp % geo % sol_zen  <  SOLZEN_375UM_NIGHT_THRESH )  then            
               is_night_375um = .false.
      end if   
               
      is_day_063um = .true.
      info_flags(1) = ibset ( info_flags(1) , 1 )
      if (( inp % geo % sol_zen >  SOLZEN_063UM_DAY_THRESH ) .or. &
                    (inp % geo %airmass > AIRMASS_THRESH  )) then           
         is_day_063um = .false.
         
      end if 
      
      is_day_063um_spatial_tests = .true.
      info_flags(1) = ibset ( info_flags(1) , 1 )
      if ( inp % geo % sol_zen >  SOLZEN_063UM_DAY_THRESH_SPATIAL_TESTS ) then           
         is_day_063um_spatial_tests = .false.
         
      end if
            
      is_forward_scatter = .false.
      if ( inp % geo %  scat_angle < 80. .and. inp % geo % sol_zen < 95.) then
         is_forward_scatter = .true.
      end if 
      
      ! --- check if DUST only for day time
      ! --- use VIIRS I1 band (ch37) if available 
      is_dust = .false.
      if ( inp % geo % sol_zen <= 85. ) then
         if ( inp % sat % chan_on (1) .and. inp % sat % chan_on (8) &
              .and. inp % sat % chan_on (37) ) then
            is_dust = DUST_DETECTION ( &
                 inp % sat % ref_ch8 &
               , inp % sat % ref_ch1 &
               , inp % sat % iband(1) % ref % std &
               , inp % geo % glint &
               , inp % sfc % land_class )
         else if ( inp % sat % chan_on (1) .and. inp % sat % chan_on (8) &
             .and. .not. inp % sat % chan_on (37) ) then
            is_dust = DUST_DETECTION ( &
                 inp % sat % ref_ch8 &
               , inp % sat % ref_ch1 &
               , inp % sat % ref_ch1_3x3_std &
               , inp % geo % glint &
               , inp % sfc % land_class ) 
         end if
      end if
      
      is_smoke = .false.
      ! --- check if SMOKE only for day time
      ! --- use VIIRS I1 band (ch37) if available
      is_smoke = .false.
      if ( inp % geo % sol_zen <= 85. ) then
         if ( inp % sat % chan_on (8) .and. inp % sat % chan_on (7) &
              .and. inp % sat % chan_on (37) ) then
            is_smoke = SMOKE_DETECTION ( &
                 inp % sat % ref_ch8 &
               , inp % sat % ref_ch7 &
               , inp % sat % iband(1) % ref % std &
               , inp % geo % sat_zen &
               , inp % geo % glint &
               , inp % sfc % land_class ) 
         else if ( inp % sat % chan_on (8) .and. inp % sat % chan_on (7) &
             .and. inp % sat % chan_on (1) .and. .not. inp % sat % chan_on (37)) then
            is_smoke = SMOKE_DETECTION ( &
                 inp % sat % ref_ch8 &
               , inp % sat % ref_ch7 &
               , inp % sat % ref_ch1_3x3_std &
               , inp % geo % sat_zen &
               , inp % geo % glint &
               , inp % sfc % land_class )
         end if
      end if
             
      ! --- check if CLOUD SHADOW
      is_cloud_shadow = .false.
      ! - TO ADD
      
      ! --- check if FIRE
      is_fire = .false.      
      if ( inp % sat % chan_on (31) &
           .and. inp % sat % chan_on (20) &
           .and. inp % sat % chan_on (41) &
           .and. inp % sat % chan_on (40) ) then
         is_fire = FIRE_DETECTION ( &
              inp % sat % bt_ch31 &
            , inp % sat % bt_ch20 &
            , inp % sat % iband(5) % ref % std &
            , inp % sat % iband(4) % ref % std &
            , inp % geo % sol_zen )
      else if ( inp % sat % chan_on (31) &
           .and. inp % sat % chan_on (20) &
           .and. .not. inp % sat % chan_on (41) &
           .and. .not. inp % sat % chan_on (40) ) then
         is_fire = FIRE_DETECTION ( &
              inp % sat % bt_ch31 &
            , inp % sat % bt_ch20 &
            , inp % rtm % bt_ch31_3x3_std &
            , inp % rtm % bt_ch20_3x3_std &
            , inp % geo % sol_zen )
      end if

      ! --- set some flags
      is_solar_contaminated = inp % geo % solar_conta
      
      info_flags = 0 
                                        info_flags(1) = ibset ( info_flags ( 1 ) , 0 )
      if ( is_day_063um )               info_flags(1) = ibset ( info_flags ( 1 ) , 1 )
      if ( is_day_063um_spatial_tests ) info_flags(1) = ibset ( info_flags ( 1 ) , 2 )
      if ( is_day_375um )               info_flags(1) = ibset ( info_flags ( 1 ) , 3 )
      if ( is_night_375um )             info_flags(1) = ibset ( info_flags ( 1 ) , 4 )
      if ( is_solar_contaminated )      info_flags(1) = ibset ( info_flags ( 1 ) , 5 )
      if ( inp % sfc % coast_mask )     info_flags(1) = ibset ( info_flags ( 1 ) , 6 )    
      if ( is_mountain )                info_flags(1) = ibset ( info_flags ( 1 ) , 7 )
      
      if ( is_forward_scatter )         info_flags(2) = ibset ( info_flags ( 2 ) , 0 )
      if ( is_cold_375um )              info_flags(2) = ibset ( info_flags ( 2 ) , 1 )
      if ( has_cold_btd )               info_flags(2) = ibset ( info_flags ( 2 ) , 2 )
      if ( is_smoke )                   info_flags(2) = ibset ( info_flags ( 2 ) , 4 )
      if ( is_dust )                    info_flags(2) = ibset ( info_flags ( 2 ) , 5 )
      if ( is_cloud_shadow )            info_flags(2) = ibset ( info_flags ( 2 ) , 6 )
      if ( is_fire )                    info_flags(2) = ibset ( info_flags ( 2 ) , 7 )
      
      ! -TODO : probably wrong
      info_flags(3) = ibits ( sfc_idx , 0 , 3 )
      
       ! diag % diagnostic_1 = -999.
       ! diag % diagnostic_2 = -999.
       ! diag % diagnostic_3 = -999.
                               
      ! - class loop 
      class_loop: do class_idx= 1,   bayes_coef % n_class
         
         ! - init
         is_on_test = .false.        
         classifier_value = -999.
         cond_yes (class_idx)  = 1.0
         cond_no  (class_idx)  = 1.0  
        
         select case ( bayes_coef % Classifier_Value_Name_enum (class_idx) )
          
         case ( et_class_T110 )
            
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( inp % sat % bt_ch31 < 0. ) cycle class_loop
            
            Classifier_Value = inp % sat % bt_ch31  
            
            is_on_test = .true.
            pos_info_flag = 4
            idx_info_flag = 3
                             
         case( et_class_TMAX_T )
            if ( .not. inp%sat % chan_on(31) ) cycle class_loop
            if ( inp % rtm % bt_ch31_3x3_max < 0. ) cycle class_loop
            if ( inp % sat % bt_ch31 < 0. ) cycle class_loop
            if ( is_mountain ) cycle class_loop
            if ( inp % sfc % coast_mask   ) cycle
            
            Classifier_Value = inp % rtm % bt_ch31_3x3_max  &
                              -  inp % sat % bt_ch31 
            is_on_test = .true.
            pos_info_flag = 6
            idx_info_flag = 3
           
         case ( et_class_T_STD )
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( inp % rtm % bt_ch31_3x3_std < 0. ) cycle class_loop
            if ( is_mountain ) cycle class_loop
            if ( inp % sfc % coast_mask ) cycle
            Classifier_Value = inp % rtm % bt_ch31_3x3_std 
            is_on_test = .true.
            pos_info_flag = 0
            idx_info_flag = 4
           
         case ( et_class_E_TROP )
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            
            Classifier_Value = inp % rtm % emis_ch31_tropo 
            
            is_on_test = .true.
            pos_info_flag = 2
            idx_info_flag = 4
           
         case  ( et_class_FMFT )
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( .not. inp % sat % chan_on(32) ) cycle class_loop
            if ( inp % rtm % bt_ch31_atm_sfc < 0. ) cycle class_loop
            if ( inp % rtm % bt_ch32_atm_sfc  < 0. ) cycle class_loop
            if ( inp % sat % bt_ch31  < 0.) cycle class_loop
            if ( inp % sat % bt_ch32  < 0. ) cycle class_loop
            
            if ( has_cold_btd ) cycle class_loop
           
            Classifier_Value = SPLIT_WINDOW_TEST ( inp % rtm % bt_ch31_atm_sfc  &
                                                 , inp % rtm % bt_ch32_atm_sfc  &
                                                 , inp % sat % bt_ch31 &
                                                 , inp % sat % bt_ch32 )
                                    
            is_on_test = .true. 
            pos_info_flag = 4
            idx_info_flag = 4
            
         case ( et_class_BTD_110_067 )
            if ( .not. inp % sat % chan_on(27) ) cycle class_loop
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( inp % sat % bt_ch31  < 0. ) cycle class_loop
            if ( inp % sat % bt_ch27  < 0. ) cycle class_loop            
            if ( has_cold_btd ) cycle class_loop
             
            Classifier_Value = inp % sat % bt_ch31 - inp % sat % bt_ch27
            
            is_on_test = .true.
            pos_info_flag = 6
            idx_info_flag = 4
        
         case ( et_class_BTD_110_067_COV )
            if ( .not. inp % sat % chan_on(27) ) cycle class_loop
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( inp % rtm % bt_ch31_ch27_covar < 0. ) cycle class_loop
            if ( has_cold_btd ) cycle class_loop
            Classifier_Value = inp % rtm % bt_ch31_ch27_covar
            is_on_test = .true.
            pos_info_flag = 0
            idx_info_flag = 5
           
         case ( et_class_BTD_110_085 )
            if ( .not. inp % sat % chan_on(29) ) cycle class_loop
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( inp % sat % bt_ch31  < 0. ) cycle class_loop
            if ( inp % sat % bt_ch29  < 0. ) cycle class_loop     
            if ( has_cold_btd ) cycle class_loop
          
            Classifier_Value = inp % sat % bt_ch31 - inp % sat % bt_ch29
            is_on_test = .true.
            pos_info_flag = 2
            idx_info_flag = 5
           
         case ( et_class_E_037 )
            
            if ( .not. inp %sat % chan_on(20) ) cycle class_loop                  
            if ( inp % geo % glint )  cycle class_loop                 
            if ( inp % sat % bt_ch20  <= 0 ) cycle class_loop
            if ( is_cold_375um ) cycle class_loop
            if ( inp % sat % emis_ch20_3x3_mean < 0.) cycle class_loop
            if ( inp % rtm % emis_ch20_clear  < 0.) cycle class_loop
           
            classifier_value = EMISS_375UM_TEST ( &
                                     inp % sat % emis_ch20_3x3_mean  &
                                   , inp % rtm % emis_ch20_clear  , is_night = .false.)
            is_on_test = .true.
            pos_info_flag = 4
            idx_info_flag = 5
           
         case ( et_class_E_037_DAY )
            if ( .not. inp % sat % chan_on(20) ) cycle class_loop
            if ( is_solar_contaminated) cycle class_loop 
            
            if ( inp % geo % glint  )  cycle class_loop
            
            if ( .not. is_day_375um ) cycle class_loop                  
            if ( is_cold_375um ) cycle class_loop
            if ( inp % sat % bt_ch20   <= 0. ) cycle class_loop
            if ( inp % sat % emis_ch20_3x3_mean < 0.) cycle class_loop
            if ( inp % rtm % emis_ch20_clear  < 0.) cycle class_loop
          
            classifier_value = EMISS_375UM_TEST ( &
                                     inp % sat % emis_ch20_3x3_mean  &
                                   , inp % rtm % emis_ch20_clear , is_night =.false.)
            is_on_test = .true.
            pos_info_flag = 6
            idx_info_flag = 5
           
         case ( et_class_E_037_NGT )
            if ( .not. inp % sat % chan_on(20) ) cycle
            if ( is_solar_contaminated ) cycle class_loop 
            if ( .not. is_night_375um ) cycle  
            if ( is_cold_375um ) cycle class_loop 
            if (  inp % sat % bt_ch20  <= 0 ) cycle class_loop
            if ( inp % sat % emis_ch20_3x3_mean < 0.) cycle class_loop
            if ( inp % rtm % emis_ch20_clear  < 0.) cycle class_loop
             
            classifier_value = EMISS_375UM_TEST ( &
                                     inp % sat % emis_ch20_3x3_mean  &
                                   , inp % rtm % emis_ch20_clear , is_night =.true.)
            
            is_on_test = .true.
            pos_info_flag = 0
            idx_info_flag = 6
              
        case ( et_class_BTD_037_110_NGT )
           if ( .not. inp % sat % chan_on(20) ) cycle
           if ( .not. inp % sat % chan_on(31) ) cycle
           if ( is_solar_contaminated ) cycle class_loop 
           if ( .not. is_night_375um ) cycle  
           if ( is_cold_375um ) cycle class_loop 
           if (  inp % sat % bt_ch20  <= 0. ) cycle class_loop
           if (  inp % sat % bt_ch31  <= 0. ) cycle class_loop
           
           classifier_value =  inp % sat % bt_ch20  -  inp % sat % bt_ch31
           
           is_on_test = .true.
           pos_info_flag = 2
           idx_info_flag = 6
              
        case ( et_class_R_006_DAY )
            ! - this solar test can be also applied for lunar
            !TODO - make own lunar visible coefficients
            !

            if ( use_lunar_refl_for_vis_tests ) then
               Classifier_Value = REFLECTANCE_GROSS_CONTRAST_TEST ( &
                                  inp % rtm % ref_dnb_clear &
                                , inp % sat % ref_dnb_lunar )
               is_on_test = .true.
               pos_info_flag = 4
               idx_info_flag = 6   
            else 
               if ( .not. inp % sat % chan_on(1) ) cycle
               if ( inp % geo % glint  )  cycle class_loop
               if ( is_forward_scatter )  cycle class_loop
               if ( is_mountain ) cycle class_loop
               if ( .not. is_day_063um ) cycle
               if ( inp % sfc % snow_class  == ET_snow_class % SNOW ) cycle
               
                            
               Classifier_Value = REFLECTANCE_GROSS_CONTRAST_TEST ( &
                                   inp % rtm % ref_ch1_clear &
                                 , inp % sat % ref_ch1 )
               is_on_test = .true.
               pos_info_flag = 4
               idx_info_flag = 6   
            end if
              
         case ( et_class_R_006_STD )
            
            if ( use_lunar_refl_for_vis_tests ) then
           
               Classifier_Value = inp % sat % ref_dnb_3x3_std 
           
               is_on_test = .true.
               pos_info_flag = 6
               idx_info_flag = 6 
            else
         
               if ( .not. inp % sat % chan_on(1) ) cycle class_loop
               if ( .not. is_day_063um_spatial_tests ) cycle
               if ( is_mountain  ) cycle class_loop
               if ( inp % sfc % coast_mask   ) cycle
               if ( .not. is_day_063um ) cycle class_loop
           
               Classifier_Value = inp % sat % ref_ch1_3x3_std 
           
               is_on_test = .true.
               pos_info_flag = 6
               idx_info_flag = 6 
            end if
            
         case ( et_class_R_006_MIN_3x3_DAY )
            if ( use_lunar_refl_for_vis_tests ) then
               Classifier_Value = RELATIVE_VISIBLE_CONTRAST_TEST ( &
                                    inp % sat % ref_dnb_3x3_min &
                                  , inp % sat % ref_dnb_lunar  ) 
                                         
               is_on_test = .true.
               pos_info_flag = 0
               idx_info_flag = 7 
            else
               if ( .not. inp % sat % chan_on(1) ) cycle class_loop
               if ( .not. is_day_063um_spatial_tests ) cycle
               if ( is_mountain  ) cycle class_loop
               if ( inp % sfc % coast_mask   ) cycle
          
               Classifier_Value = RELATIVE_VISIBLE_CONTRAST_TEST ( &
                                    inp % sat % ref_ch1_3x3_min &
                                  , inp % sat % ref_ch1  ) 
                                         
               is_on_test = .true.
               pos_info_flag = 0
               idx_info_flag = 7    
           end if
           
         case ( et_class_R_RATIO_DAY )
           if ( .not. inp %  sat % chan_on(1) ) cycle class_loop
           if ( .not. inp % sat % chan_on(2) ) cycle class_loop
           if ( .not. is_day_063um_spatial_tests ) cycle
           if ( is_mountain  ) cycle class_loop
           if ( inp % geo  % glint  )  cycle class_loop
          
           Classifier_Value = REFLECTANCE_RATIO_TEST (&
                                   inp % sat % ref_ch1 &
                                 , inp % sat % ref_ch2 ) 
           is_on_test = .true.
           pos_info_flag = 2
           idx_info_flag = 7  
              
         case ( et_class_R_013_DAY )
            if ( .not. inp %  sat % chan_on(26) ) cycle class_loop
            if ( is_forward_scatter )  cycle class_loop
            if ( .not. is_day_063um )  cycle class_loop
            if ( is_mountain ) cycle class_loop
            if ( inp % sat % ref_ch26  <= 0 ) cycle class_loop
           
            Classifier_Value = inp % sat % ref_ch26
            is_on_test = .true.   
            pos_info_flag = 4
            idx_info_flag = 7    
            
         case ( et_class_NDSI_DAY )
            if ( .not. inp % sat % chan_on(1) ) cycle class_loop
            if ( .not. inp % sat % chan_on(6) ) cycle class_loop
            if ( is_forward_scatter )  cycle class_loop
            if ( .not. is_day_063um ) cycle class_loop
            if ( inp % geo % glint )  cycle class_loop
            if ( inp % sat % ref_ch1 <= 0 ) cycle class_loop
            if ( inp % sat % ref_ch6 <= 0 ) cycle class_loop
           
            Classifier_Value = (inp % sat % ref_ch1 - inp % sat % ref_ch6) / &
                               (inp % sat % ref_ch1 + inp % sat % ref_ch6)
            is_on_test=.true.
            pos_info_flag = 6
            idx_info_flag = 7  
              
         case default
            print*,'unknown class ', bayes_coef % Classifier_Value_Name_enum (class_idx) &
                                   , bayes_coef % Classifier_Value_Name (class_idx,1)
         end select
       
       ! - AW 09/17/2014 
       ! - test can be switched off in bayesian file  
       if ( .not. bayes_coef % do_this_classifier (Class_Idx,Sfc_Idx) )  is_on_test = .false.
       
       ! --- all tests classifer
        if ( is_on_test ) then
            
            bin_idx = int (( classifier_value &
                 -  bayes_coef % Classifier_Bounds_Min(Class_Idx,Sfc_Idx))  /    &
                    bayes_coef % Delta_Classifier_Bounds(Class_Idx,Sfc_Idx)) + 1
             
            bin_idx =  max(1,min(bayes_coef % N_bounds-1,Bin_Idx)) 
            Cond_Yes (Class_Idx) = Bayes_Coef % Class_Cond_Yes ( Bin_Idx, Class_Idx, Sfc_Idx )
            Cond_No (Class_Idx)  = Bayes_Coef % Class_Cond_No ( Bin_Idx, Class_Idx, Sfc_Idx ) 
            Cond_Ratio (Class_Idx) = Bayes_Coef % Class_Cond_Ratio ( Bin_Idx, Class_Idx, Sfc_Idx )

            r = Cond_Ratio(Class_Idx)
            class_contr =  1.0 / (1.0 + r / Bayes_Coef % Prior_Yes(Sfc_Idx) - r)
!            class_contr =  bayes_coef % Prior_Yes(Sfc_Idx)* bayes_coef % class_cond_yes ( bin_idx, class_idx , sfc_idx ) / &
!                            ( bayes_coef % Prior_Yes(Sfc_Idx) * bayes_coef % class_cond_yes ( bin_idx, class_idx , sfc_idx )  &
!                            + bayes_coef % Prior_No(Sfc_Idx) * bayes_coef % class_cond_no ( bin_idx, class_idx , sfc_idx ) )
                            
            
            if ( class_contr > 0.1 .and. class_contr < 0.5)  info_flags ( idx_info_flag) = &
                            ibset ( info_flags ( idx_info_flag) , pos_info_flag )
            if ( class_contr >= 0.5 )  info_flags ( idx_info_flag) = ibset ( info_flags ( idx_info_flag) , pos_info_flag + 1)
            if ( class_contr > 0.9 ) then
                info_flags ( idx_info_flag ) = ibset ( info_flags ( idx_info_flag ) , pos_info_flag )
                info_flags ( idx_info_flag ) = ibset ( info_flags ( idx_info_flag ) , pos_info_flag + 1 )
            end if
           
         end if

         !-----------------------------------------------------------------------
         ! --- Diagnostic Output for debugging in CLAVR-x
         select case (  bayes_coef % Classifier_Value_Name_enum (class_idx))
           case ( et_class_T110 )
           case ( et_class_TMAX_T )
           case ( et_class_T_STD )
           case ( et_class_E_TROP )
           case ( et_class_FMFT )
           case ( et_class_BTD_110_067 )
           case ( et_class_BTD_110_067_COV )
           case ( et_class_BTD_110_085 )
           case ( et_class_E_037 )
           case ( et_class_E_037_DAY )
              !diag % diagnostic_1 =    classifier_value
              !diag % diagnostic_2 = class_contr
              !diag % diagnostic_3 = 1.  
           case ( et_class_E_037_NGT )
           case ( et_class_BTD_037_110_NGT )
           case ( et_class_R_006_DAY )
           case ( et_class_R_006_STD )
           case ( et_class_R_006_MIN_3x3_DAY )
           case ( et_class_R_RATIO_DAY )
           case ( et_class_R_013_DAY )
           case ( et_class_NDSI_DAY )
             
         end select
         !
                
      end do class_loop
                        
      ! - compute Posterior_Cld_Probability
           
      r = product(Cond_Ratio)
      erg =  1.0 / (1.0 + r / Bayes_Coef % Prior_Yes(Sfc_Idx) - r)
!     erg  = &
!           (bayes_coef % Prior_Yes(Sfc_Idx) * product(Cond_Yes)) / &
!           (bayes_coef % Prior_Yes(Sfc_Idx) * product(Cond_Yes) +  &        
!            bayes_coef % Prior_No(Sfc_Idx) * product(Cond_No))
                
            deallocate ( Cond_Yes )
            deallocate ( Cond_No )
            deallocate ( Cond_Ratio )
            
   end subroutine NB_CLOUD_MASK_ALGORITHM

   !-------------------------------------------------------------------------------------
   !   Reads the coefficients from file
   !-------------------------------------------------------------------------------------
   subroutine READ_BAYES_COEFF ( )

      use FILE_TOOLS, only : &
          GETLUN

      implicit none
      
      integer :: lun
      integer :: ios
     
      real :: time_diff_max
      integer :: n_class , n_bounds , n_sfc_bayes
      
      character (120) :: first_line
      character ( 72) :: header
      integer :: i_class
      integer :: i_sfc      ! loop variable
      integer :: i_sfc_file  ! number of surface in file
      
      integer :: int_dummy
      integer :: index_start
      integer :: index_end
      
      
      lun = GETLUN()
      open (unit=lun, file = trim(bayes_coef % file) , &
            action ="read", form="formatted", status="old", iostat=ios)
      if (ios/= 0) then
         print*, 'no bayesain cloud mask ', trim(bayes_coef % file)
         return
      end if
      
      ! - first two lines are header meta-data
      read ( unit=lun, fmt="(a120)") first_line 
      
      ! - in case we have an additional cvs id version line - very risky!
      !
      if ( index (first_line, '$Id:') > 0 ) then
         bayes_coef % cvs_version =  first_line
         read (unit=lun, fmt=*) header
      else
         bayes_coef % cvs_version= 'unknown'
         header = first_line
      end if
      
      read (unit=lun, fmt=*)
      read (unit=lun, fmt=*) time_diff_max, n_class, n_bounds, n_sfc_bayes
      bayes_coef % n_class = n_class
      bayes_coef % n_bounds = n_bounds

      call bayes_coef % alloc ( n_class, n_bounds, n_sfc_bayes )
      
      bayes_coef % do_this_classifier (:,:) = .true. 
      
      do i_sfc = 1 , n_sfc_bayes
         read (unit=lun, fmt=*, iostat=ios) i_sfc_file, Header
         read (unit=lun, fmt=*, iostat=ios) bayes_coef % Prior_Yes(i_sfc), bayes_coef % Prior_No(i_sfc)
         read (unit=lun, fmt=*, iostat=ios) bayes_coef % Optimal_Posterior_Prob(i_Sfc)
         
         do i_class = 1 , n_class
            
            read (unit=lun, fmt=*, iostat=ios)  &
                             int_dummy &
                          ,  index_start &
                          ,  index_end &
                          ,  bayes_coef %Classifier_Value_Name(i_class,i_sfc)
           
            
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'T_11' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_T110
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'T_max-T' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_TMAX_T
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'T_Std' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_T_STD
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Emiss_Tropo' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_E_TROP
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'FMFT' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_FMFT
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Btd_11_67' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_BTD_110_067
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Bt_11_67_Covar' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_BTD_110_067_COV
                          
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Btd_11_85' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_BTD_110_085
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Emiss_375' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_E_037
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Emiss_375_Day' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_E_037_DAY
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Emiss_375_Night' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_E_037_NGT
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Btd_375_11_Night' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_BTD_037_110_NGT
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Ref_063_Day' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_R_006_DAY
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Ref_std' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_R_006_STD
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Ref_063_Min_3x3_Day' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_R_006_MIN_3x3_DAY
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Ref_Ratio_Day' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_R_RATIO_DAY
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Ref_138_Day' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_R_013_DAY
               
            if ( trim(bayes_coef % Classifier_Value_Name(i_class,i_sfc)) == 'Ndsi_Day' ) &
                 bayes_coef % Classifier_Value_Name_enum(i_class) =  et_class_NDSI_DAY
                        
            read (unit=lun, fmt=*, iostat=ios)  &
                                 bayes_coef % Classifier_Bounds_Min(i_class,i_sfc)  &
                             ,   bayes_coef % Classifier_Bounds_Max(i_class,i_sfc)  &
                             ,   bayes_coef % Delta_Classifier_Bounds(i_class,i_sfc)
                                
             read (unit=lun, fmt=*, iostat=ios)  bayes_coef % Class_Cond_Yes (:,i_class,i_sfc) 
             read (unit=lun, fmt=*, iostat=ios)  bayes_coef % Class_Cond_No (:,i_class,i_sfc)
             
             ! --   AW 09/17/2014
             ! -- to switch off this class for this surface do set the start_index in the bayesian file to a 
             ! -- higher value than end_index
             ! -- There should be better ways to do this, but we have to adjust this first in the program (IDL) where we create 
             ! -- the bayesian files
             ! --  
            
             if ( index_start >= index_end) then
               bayes_coef % do_this_classifier (i_class,i_sfc) = .false.               
             end if
         end do
      end do
         
      bayes_coef % is_read = .true.
   
   end subroutine READ_BAYES_COEFF

!-------------------------------------------------------------------------------
! NetCDF LUT routines:
!-------------------------------------------------------------------------------

!====================================================================
! SUBROUTINE Name: READ_NAIVE_BAYES_NC
!
! Function:
!   Allocate and Read in the LUTs needed for Bayesian cloud mask tables
!
!====================================================================
 subroutine READ_BAYES_COEFF_NC ()

   implicit none

   !local variables
   integer:: Ncid
   integer:: Status
   integer:: i_class
   character(30):: Var_Name

   Bayes_Coef % Is_Read = .FALSE.

   print *,'Reading Bayesian Coefficients From -> ',trim(Bayes_Coef % file)
   Status = nf90_open(trim(Bayes_Coef % file), mode = nf90_nowrite, ncid = Ncid)
   if (Status /= nf90_noerr) then
      print *, 'ERROR: Bayesian Cloud Mask Classifier Open Failed '
      print *, 'Bayesian Cloud Mask Turned Off'
      return
   endif

   Status = nf90_get_att(Ncid, nf90_global, "data_file", bayes_coef % cvs_version)
   if (Status /= nf90_noerr) then
      print *, 'ERROR: Bayesian Cloud Mask Version Read Failed'
      return
   endif

   Status = nf90_get_att(Ncid, nf90_global, "n_class", Bayes_Coef % N_Class)
   Status = nf90_get_att(Ncid, nf90_global, "n_bounds_reg", Bayes_Coef % N_Bounds)
   Status = nf90_get_att(Ncid, nf90_global, "n_sfc_type", Bayes_Coef % N_Sfc_Bayes)

   !--- allocate
   call Bayes_Coef % Alloc ( Bayes_Coef % N_Class, Bayes_Coef % N_Bounds, &
                                Bayes_Coef % N_Sfc_Bayes )
   Bayes_Coef % Do_This_Classifier (:,:) = .true.

   !--- initialize
   Bayes_Coef % Prior_Yes = -999.
   Bayes_Coef % Prior_No = -999.
   Bayes_Coef % Optimal_Posterior_Prob = -999.
   Bayes_Coef % First_Valid_Classifier_Bounds = 0
   Bayes_Coef % Last_Valid_Classifier_Bounds = 0

   !Now read in to the 1D variables

   Var_name="prior_yes"
   call read_netcdf_1d_real ( Ncid, Bayes_Coef % N_Sfc_Bayes, Var_Name,Bayes_Coef % Prior_Yes )

   Var_name="prior_no"
   call read_netcdf_1d_real ( Ncid, Bayes_Coef % N_Sfc_Bayes, Var_Name,Bayes_Coef % Prior_No )

   Var_name="optimal_posterior_prob"
   call read_netcdf_1d_real ( Ncid, Bayes_Coef % N_Sfc_Bayes, Var_Name, &
                                    Bayes_Coef % Optimal_Posterior_Prob )

   !Now the 2D variables

   Sds_Start_2d = 1
   Sds_Edge_2d(1) = Bayes_Coef % N_Class
   Sds_Edge_2d(2) = Bayes_Coef % N_Sfc_Bayes

   Var_Name="bin_start"
   call read_netcdf_2d_real(Ncid, Sds_Start_2d, Sds_Edge_2d, &
                              Var_Name, Bayes_Coef % Classifier_Bounds_Min)

   Var_Name="bin_end" !real
   call read_netcdf_2d_real(Ncid, Sds_Start_2d, Sds_Edge_2d, &
                              Var_Name, Bayes_Coef % Classifier_Bounds_Max)

   Var_Name="delta_bin" !real
   call read_netcdf_2d_real(Ncid, Sds_Start_2d, Sds_Edge_2d, &
                              Var_Name, Bayes_Coef % Delta_Classifier_Bounds)

   Var_Name="first_valid_bounds" !integer
   call read_netcdf_2d_int(Ncid, Sds_Start_2d, Sds_Edge_2d, &
                              Var_Name, Bayes_Coef % First_Valid_Classifier_Bounds)

   Var_Name="last_valid_bounds" !integer
   call read_netcdf_2d_int(Ncid, Sds_Start_2d, Sds_Edge_2d, &
                              Var_Name, Bayes_Coef % Last_Valid_Classifier_Bounds)

   Var_Name="classifier_names" !character
   call read_netcdf_2d_char(Ncid, Sds_Start_2d, Sds_Edge_2d, &
                              Var_Name, Bayes_Coef % Classifier_Value_Name)

   !finally 3D variables
   Sds_Start_3d = 1
   Sds_Edge_3d(1) = Bayes_Coef % N_bounds-1
   Sds_Edge_3d(2) = Bayes_Coef % N_class
   Sds_Edge_3d(3) = Bayes_Coef % N_sfc_bayes

   Var_Name="class_cond_yes_reg" !real
   call read_netcdf_3d(ncid, Sds_Start_3d, Sds_edge_3d, &
                             Var_Name, Bayes_Coef % Class_Cond_Yes)

   Var_Name="class_cond_no_reg" !real
   call read_netcdf_3d(Ncid, Sds_Start_3d, Sds_Edge_3d, &
                             Var_Name, Bayes_Coef % Class_Cond_No)

   Var_Name="class_cond_ratio_reg" !real
   call read_netcdf_3d(Ncid, Sds_Start_3d, Sds_Edge_3d, &
                              Var_Name, Bayes_Coef % Class_Cond_Ratio)

   Status = nf90_close(Ncid)
 
   where (Bayes_Coef % First_Valid_Classifier_Bounds >= &
       Bayes_Coef % Last_Valid_Classifier_Bounds)
      Bayes_Coef % Do_This_Classifier = .false.
   endwhere

  ! --- set Classifier_Value_Name_enum
  do i_class = 1 , Bayes_Coef % N_Class

     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'T_11' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_T110
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'T_max-T' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_TMAX_T
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'T_Std' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_T_STD
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Emiss_Tropo' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_E_TROP
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'FMFT' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_FMFT
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Btd_11_67' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_BTD_110_067
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Bt_11_67_Covar' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_BTD_110_067_COV
                          
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Btd_11_85' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_BTD_110_085
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Emiss_375' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_E_037
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Emiss_375_Day' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_E_037_DAY
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Emiss_375_Night' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_E_037_NGT
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Btd_375_11_Night' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_BTD_037_110_NGT
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Ref_063_Day' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_R_006_DAY
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Ref_std' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_R_006_STD
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Ref_063_Min_3x3_Day' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_R_006_MIN_3x3_DAY
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Ref_Ratio_Day' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_R_RATIO_DAY
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Ref_138_Day' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_R_013_DAY
               
     if ( trim(bayes_coef % Classifier_Value_Name(i_class,1)) == 'Ndsi_Day' ) &
         bayes_coef % Classifier_Value_Name_enum(i_class) = et_class_NDSI_DAY
  enddo

   Bayes_Coef % Is_Read = .true.

 end subroutine READ_BAYES_COEFF_NC

   
   !-------------------------------------------------------------------------------------
   !  7 sfc bayesian surface  types 
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   function BAYES_SFC_TYPE ( lat , lon , land_class , coast , snow_class , sfc_type , emis_ch20 , sst_anal_uni)
     
       
      real  :: lat , lon , emis_ch20 , sst_anal_uni
      integer :: land_class , snow_class , sfc_type
      logical :: coast
      
      integer :: bayes_sfc_type
      integer , parameter :: CLOSED_SHRUBS_SFC = 8
      integer , parameter :: OPEN_SHRUBS_SFC = 9
      integer , parameter :: GRASSES_SFC = 10
      integer , parameter :: BARE_SFC = 12
     
      bayes_sfc_type = 0
     
      ! # deep ocean 
      if ( land_class == ET_land_class % DEEP_OCEAN ) then
         bayes_sfc_type = 1
      end if
      
      ! # 2 - shallow ocean/water
      if (  land_class  == ET_land_class % MODERATE_OCEAN .or. &
               land_class  == ET_land_class % DEEP_INLAND_WATER .or. &
               land_class  == ET_land_class % SHALLOW_INLAND_WATER .or. &
               land_class == ET_land_class % SHALLOW_OCEAN )  then
         bayes_sfc_type = 2
      end if   
     
      if ( land_class /= ET_land_class % LAND .and. sst_anal_uni > 0.5 ) then
         bayes_sfc_type = 2
      end if
     
      ! # 3 unfrozen land
      if (  land_class  == ET_land_class % LAND .or. &
               land_class  == ET_land_class % COASTLINE .or. &
               land_class  == ET_land_class % EPHEMERAL_WATER .or. &  
               coast  )  then
         bayes_sfc_type = 3
      end if
      

     
      ! -- # 4 snow covered land
      if ( (lat > - 60.) .and.  snow_class  == ET_snow_class % SNOW )  then
          bayes_sfc_type = 4
      end if
     
      ! -- # 5 Arctic
      if ( (lat >= 0.) .and.  snow_class  == ET_snow_class % SEA_ICE )  then
          bayes_sfc_type = 5
      end if
      
      ! -- # 6 Antarctic
      if ( (lat <= -60.) .and.  snow_class  == ET_snow_class % SNOW ) then
          bayes_sfc_type = 6
      end if
      
      if ( (lat <= -60.) .and.  snow_class  == ET_snow_class % SEA_ICE ) then
          bayes_sfc_type = 6
      end if
      
      ! -- #  6 greenland
      if ( lat >= 60.0 &
              .and. lon > -75. &
              .and. lon < -10.0  &    
              .and. (  land_class == ET_land_class % LAND &
              .or.  land_class  == ET_land_class % COASTLINE ) &
              .and.  snow_class  == ET_snow_class % SNOW ) then
               
          bayes_sfc_type = 6         
      end if     
         
      !-TODO
      
      ! -- # 7 Desert
      
      if ( emis_ch20 < 0.90 .and. abs ( lat ) < 60.0 &
         .and. ( sfc_type == OPEN_SHRUBS_SFC .or. sfc_type == BARE_SFC ) ) then
      
         bayes_sfc_type = 7
      end if   
      
      if ( bayes_sfc_type == 3 .and. &
             emis_ch20 < 0.93 .and. &
             abs ( lat ) < 60.0 .and. &
             ( sfc_type == OPEN_SHRUBS_SFC .or. &
             sfc_type == CLOSED_SHRUBS_SFC .or. &
             sfc_type == GRASSES_SFC .or. &
             sfc_type == BARE_SFC ) ) then
            
            bayes_sfc_type = 7
      end if
     
   end function BAYES_SFC_TYPE
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   subroutine ALLOC_BAYES_COEF ( this, n_class , n_bounds , n_sfc_bayes )
      class ( bayes_coef_type ) :: this
      integer :: n_class , n_bounds , n_sfc_bayes
      
      allocate (this % Prior_Yes(N_sfc_bayes))
      allocate (this % Prior_No(N_sfc_bayes))
      allocate (this % Optimal_Posterior_Prob(N_sfc_bayes))
      allocate (this % Classifier_Bounds_Min(N_class,N_sfc_bayes))
      allocate (this % Classifier_Bounds_Max(N_class,N_sfc_bayes))
      allocate (this % Delta_Classifier_Bounds(N_class,N_sfc_bayes))
      allocate (this % First_valid_Classifier_Bounds(N_class,N_sfc_bayes))
      allocate (this % Last_valid_Classifier_Bounds(N_class,N_sfc_bayes))
      allocate (this % do_this_classifier(N_class,N_sfc_bayes))
      allocate (this % Class_Cond_Yes(N_bounds-1,N_class,N_sfc_bayes))
      allocate (this % Class_Cond_No(N_bounds-1,N_class,N_sfc_bayes))
      allocate (this % Class_Cond_Ratio(N_bounds-1,N_class,N_sfc_bayes))
      allocate (this % Classifier_Value_Name(N_class,N_sfc_bayes))
      allocate (this % Classifier_Value_Name_enum(N_class))
      allocate (this % Flag_Idx(N_class))
            
   end subroutine ALLOC_BAYES_COEF
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   subroutine DEALLOC_BAYES_COEF ( this )
      class ( bayes_coef_type ) :: this
           
      deallocate (this % Prior_Yes)
      deallocate (this % Prior_No)
      deallocate (this % Optimal_Posterior_Prob)
      deallocate (this % do_this_classifier)
      deallocate (this % Classifier_Bounds_Min)
      deallocate (this % Classifier_Bounds_Max)
      deallocate (this % Delta_Classifier_Bounds)
      deallocate (this % First_valid_Classifier_Bounds)
      deallocate (this % Last_valid_Classifier_Bounds)
      deallocate (this % Class_Cond_Yes)
      deallocate (this % Class_Cond_No)
      deallocate (this % Class_Cond_Ratio)
      deallocate (this % Classifier_Value_Name)
      deallocate (this % Classifier_Value_Name_enum)
      deallocate (this % Flag_Idx)
            
   end subroutine DEALLOC_BAYES_COEF
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   real elemental function REFLECTANCE_RATIO_TEST ( ref_vis , ref_nir )
      real, intent(in):: ref_vis
      real, intent(in):: ref_nir
      real, parameter :: MISSING_VALUE_REAL4 = -999.
      
      reflectance_ratio_test = Missing_Value_Real4
      if (ref_vis /= Missing_Value_Real4 .and. ref_nir /= Missing_Value_Real4) then
           reflectance_ratio_test = ref_nir / ref_vis
      end if

   end function REFLECTANCE_RATIO_TEST
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------   
   real elemental function RELATIVE_VISIBLE_CONTRAST_TEST ( ref_min , ref )
      real, intent(in):: ref_min
      real, intent(in):: ref
      real, parameter :: MISSING_VALUE_REAL4 = -999.
      
      relative_visible_contrast_test = Missing_Value_Real4
      if (ref_min /= Missing_Value_Real4) then
           relative_visible_contrast_test = ref - ref_min
      end if

   end function RELATIVE_VISIBLE_CONTRAST_TEST
   
   !-------------------------------------------------------------------------------------
   !
   !------------------------------------------------------------------------------------- 
   real elemental function REFLECTANCE_GROSS_CONTRAST_TEST ( ref_clear , ref )
      real, intent(in):: ref_clear
      real, intent(in):: ref
      real, parameter :: MISSING_VALUE_REAL4 = -999.
      reflectance_gross_contrast_test = Missing_Value_Real4
      if (ref_clear /= Missing_Value_Real4) then
           reflectance_gross_contrast_test = ref - ref_clear
      endif

   end function REFLECTANCE_GROSS_CONTRAST_TEST
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   real elemental function SPLIT_WINDOW_TEST ( t11_clear, t12_clear, t11, t12)

      real, intent(in):: t11_clear
      real, intent(in):: t12_clear 
      real, intent(in):: t11
      real, intent(in):: t12

      split_window_test  = (t11_clear - t12_clear) * (t11 - 260.0) / (t11_clear - 260.0)

      if (t11_clear <= 265.0) split_window_test = 0.0

      split_window_test = (t11 - t12) - split_window_test

   end function SPLIT_WINDOW_TEST
  
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   real elemental function EMISS_375UM_TEST ( ems, ems_clear , is_night  )
      real, intent(in):: ems_clear
      real, intent(in):: ems
      logical ,intent(in) :: is_night
      
      real, parameter :: MISSING_VALUE_REAL4 = -999.

      emiss_375um_test = MISSING_VALUE_REAL4
      if (ems /= Missing_Value_Real4 .and. ems_clear /= Missing_Value_Real4) then
           emiss_375um_test = (ems- ems_clear) / ems_clear
      end if
      if ( is_night) emiss_375um_test = ems
  
  
   end function EMISS_375UM_TEST
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------   
   logical elemental function FIRE_DETECTION ( &
           bt_110 &
         , bt_037 &
         , bt_110_std &
         , bt_037_std &
         , sol_zen )
         
      real, intent(in):: bt_110
      real, intent(in):: bt_037
      real, intent(in):: bt_110_std
      real, intent(in):: bt_037_std
      real, intent(in):: sol_zen
      
      real, parameter :: EUMETCAST_FIRE_DAY_SOLZEN_THRESH = 70.0
      real, parameter :: EUMETCAST_FIRE_NIGHT_SOLZEN_THRESH = 90.0
      
      !---- EUMETCAST fire detection parameters
      real, parameter :: BT_375UM_EUMET_FIRE_DAY_THRESH = 310.0
      real, parameter :: BT_DIFF_EUMET_FIRE_DAY_THRESH = 8.0
      real, parameter :: STDDEV_11UM_EUMET_FIRE_DAY_THRESH = 1.0 
      real, parameter :: STDDEV_375UM_EUMET_FIRE_DAY_THRESH = 4.0
      
      real, parameter :: BT_375UM_EUMET_FIRE_NIGHT_THRESH = 290.0
      real, parameter :: BT_DIFF_EUMET_FIRE_NIGHT_THRESH = 0.0
      real, parameter :: STDDEV_11UM_EUMET_FIRE_NIGHT_THRESH = 1.0 
      real, parameter :: STDDEV_375UM_EUMET_FIRE_NIGHT_THRESH = 4.0
      
      real :: Bt_375um_Eumet_Fire_Thresh
      real :: Bt_Diff_Eumet_Fire_Thresh
      real :: Stddev_11um_Eumet_Fire_Thresh 
      real :: Stddev_375um_Eumet_Fire_Thresh
      
      integer :: part_of_day
      integer, parameter :: ET_part_of_day_DAY = 1
      integer, parameter :: ET_part_of_day_NIGHT = 2
      integer, parameter :: ET_part_of_day_TWILIGHT = 3

      fire_detection = .false.
                 
      ! - check if valid
      if ( bt_110_std < 0. .or. bt_037_std < 0.) then
         return
      end if 
      
      part_of_day = ET_part_of_day_TWILIGHT
      if ( sol_zen < EumetCAST_Fire_Day_Solzen_Thresh )  part_of_day = ET_part_of_day_DAY
      if ( sol_zen > EumetCAST_Fire_Night_Solzen_Thresh )  part_of_day = ET_part_of_day_NIGHT
      
      select case ( part_of_day )
      
      case( ET_part_of_day_DAY)
         Bt_375um_Eumet_Fire_Thresh       = BT_375UM_EUMET_FIRE_DAY_THRESH
         Bt_Diff_Eumet_Fire_Thresh        = BT_DIFF_EUMET_FIRE_DAY_THRESH
         Stddev_11um_Eumet_Fire_Thresh    = STDDEV_11UM_EUMET_FIRE_DAY_THRESH
         Stddev_375um_Eumet_Fire_Thresh   = STDDEV_375UM_EUMET_FIRE_DAY_THRESH
      
      case( ET_part_of_day_NIGHT)
         Bt_375um_Eumet_Fire_Thresh       = BT_375UM_EUMET_FIRE_NIGHT_THRESH
         Bt_Diff_Eumet_Fire_Thresh        = BT_DIFF_EUMET_FIRE_NIGHT_THRESH
         Stddev_11um_Eumet_Fire_Thresh    = STDDEV_11UM_EUMET_FIRE_NIGHT_THRESH
         Stddev_375um_Eumet_Fire_Thresh   = STDDEV_375UM_EUMET_FIRE_NIGHT_THRESH         
      
      case( ET_part_of_day_TWILIGHT)
         Bt_375um_Eumet_Fire_Thresh = ((-1.0)* sol_zen) + 380.0
         Bt_Diff_Eumet_Fire_Thresh = ((-0.4)* sol_zen) + 36.0
         Stddev_11um_Eumet_Fire_Thresh = ((0.0)* sol_zen) + 1.0
         Stddev_375um_Eumet_Fire_Thresh = ((0.0)* sol_zen) + 4.0
         
      end select
      
      ! All of these conditions need to be met
      if (  ( bt_037                > Bt_375um_Eumet_Fire_Thresh ) .and. &
            ( ( bt_037 - bt_110 )   > Bt_Diff_Eumet_Fire_Thresh) .and. &
            ( bt_037_std            > Stddev_375um_Eumet_Fire_Thresh) .and. &
            ( bt_110_std            < Stddev_11um_Eumet_Fire_Thresh) ) then
         Fire_detection = .true.
       endif
   
   
   end function FIRE_DETECTION
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   logical elemental function SMOKE_DETECTION ( &
           ref_004 &
         , ref_021 &
         , ref_006_std &
         , sat_zen &
         , is_glint &
         , land_class )
               
      real , intent(in) :: ref_004
      real , intent(in) :: ref_021
      real , intent(in) :: ref_006_std
      real , intent(in) :: sat_zen
      logical , intent(in) :: is_glint
      integer , intent(in) :: land_class  
      
      logical :: is_water_sfc      
     
      real, parameter :: pi = 3.14159265359
      real, parameter :: SMOKE_ST_DEV_LAND_THRESH = 0.1
      real, parameter :: SMOKE_ST_DEV_WATER_THRESH = 0.05
      real, parameter :: SMOKE_ST_DEV_LAND_GLINT_THRESH = 1.5
      real, parameter :: SMOKE_ST_DEV_WATER_GLINT_THRESH = 0.05
      real, parameter :: SMOKE_CAND_M11M1_REF_RATIO_THRESH = 0.25
      
      real :: ref_ratio
      real :: st_dev_thresh_cand 
      
      smoke_detection = .false.
            
      ! - check if valid
      if ( ref_004 < 0. .or. ref_021 < 0. .or. ref_006_std < 0.) then
         return
      end if 
      
      is_water_sfc = land_class == 0 .or. &
         & land_class >= 3 .and. land_class <= 7 
         
      ref_ratio = ref_021 / ref_004 
      
      if ( is_water_sfc .and. ref_ratio > ( SMOKE_CAND_M11M1_REF_RATIO_THRESH * cos (sat_zen*pi/180.0) ) ) return
      
      if ( is_glint .and. is_water_sfc )                     st_dev_thresh_cand = SMOKE_ST_DEV_WATER_GLINT_THRESH
      if ( is_glint .and. .not. ( is_water_sfc ) )           st_dev_thresh_cand = SMOKE_ST_DEV_LAND_GLINT_THRESH
      if ( .not. ( is_glint ) .and. is_water_sfc )           st_dev_thresh_cand = SMOKE_ST_DEV_WATER_THRESH
      if ( .not. ( is_glint ) .and. .not. ( is_water_sfc ) ) st_dev_thresh_cand = SMOKE_ST_DEV_LAND_THRESH
      
      if ( ref_006_std < st_dev_thresh_cand ) smoke_detection = .true.
      
   end function SMOKE_DETECTION
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   logical elemental function DUST_DETECTION ( &
              ref_004 &
            , ref_006 &
            , ref_006_std &
            , is_glint &
            , land_class )
            
      real , intent(in) :: ref_004
      real , intent(in) :: ref_006
      real , intent(in) :: ref_006_std
      logical , intent(in) :: is_glint
      integer , intent(in) :: land_class  
      
      logical :: is_water_sfc      
      
      real, parameter :: DUST_M1_REFL_THRESH = 0.8
      real, parameter :: DUST_CAND_M1M5_REFL_RATIO_THRESH = 0.25
      real, parameter :: DUST_ST_DEV_LAND_THRESH = 0.1
      real, parameter :: DUST_ST_DEV_WATER_THRESH = 0.05
      real, parameter :: DUST_ST_DEV_LAND_GLINT_THRESH = 0.7
      real, parameter :: DUST_ST_DEV_WATER_GLINT_THRESH = 0.05
      
      real :: st_dev_thresh_cand   
      
      dust_detection = .false.
      
      ! - check if valid
      if ( ref_004 < 0. .or. ref_006 < 0. .or. ref_006_std <= 0.) then
         return
      end if 
      
      ! -- exclude water 
      is_water_sfc = land_class == 0 .or. &
         & land_class >= 3 .and. land_class <= 7 
       
         
      if ( is_water_sfc .and. ( ref_004 >= DUST_M1_REFL_THRESH &
               .or. Ref_004 / Ref_006 >= DUST_CAND_M1M5_REFL_RATIO_THRESH )) then
         return            
      end if  
      
      ! adjust thresholds
      if ( is_glint .and. is_water_sfc )                 st_dev_thresh_cand = DUST_ST_DEV_WATER_GLINT_THRESH
      if ( is_glint .and. .not. (is_water_sfc) )         st_dev_thresh_cand = DUST_ST_DEV_LAND_GLINT_THRESH
      if ( .not. (is_glint) .and. is_water_sfc )         st_dev_thresh_cand = DUST_ST_DEV_WATER_THRESH
      if ( .not. (is_glint) .and. .not. (is_water_sfc) ) st_dev_thresh_cand = DUST_ST_DEV_LAND_THRESH
      
      if (  ref_006_std < st_dev_thresh_cand ) dust_detection = .true.
   
   end function DUST_DETECTION

   ! ----------------------------------------------------------
   ! Read in 1D arrays (used code from DCOMP reader
   ! ----------------------------------------------------------
   subroutine read_netcdf_1d_real (nc_file_id, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: var_dim
      character(30), intent(in) :: var_name
      real, intent(out), dimension(:) :: var_output

      integer :: nc_var_id
      integer :: status

      Sds_Start_1D = 1
      Sds_Stride_1D = 1
      Sds_Edge_1D = var_dim

      status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=Sds_Start_1D, count=Sds_Edge_1D)
      if (status /= nf90_noerr) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

   end subroutine read_netcdf_1d_real                                                                                              

   ! ----------------------------------------------------------
   ! Read in 1D arrays (used code from DCOMP reader
   ! ----------------------------------------------------------
   subroutine read_netcdf_1d_int (nc_file_id, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: var_dim
      character(len=*), intent(in) :: var_name
      integer, intent(out), dimension(:) :: var_output

      integer :: nc_var_id
      integer :: status

      Sds_Start_1D = 1
      Sds_Stride_1D = 1
      Sds_Edge_1D = var_dim

      status = nf90_inq_varid(nc_file_id,trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=Sds_Start_1D, count=Sds_Edge_1D)
      if (status /= nf90_noerr) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

   end subroutine read_netcdf_1d_int

   ! ----------------------------------------------------------
   ! Read in 2D arrays (used code from DCOMP reader)
   ! ----------------------------------------------------------
   subroutine read_netcdf_2d_real (nc_file_id, start_var, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: start_var(:)
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      real, intent(out), dimension(:,:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=start_var, count=var_dim)
      if ((status /= nf90_noerr)) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

   end subroutine read_netcdf_2d_real

   ! ----------------------------------------------------------
   ! Read in 2D arrays Integers
   ! ----------------------------------------------------------
   subroutine read_netcdf_2d_int (nc_file_id, start_var, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: start_var(:)
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      integer , intent(out), dimension(:,:) :: var_output

      integer :: nc_var_id
      integer :: status

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=start_var, count=var_dim)
      if ((status /= nf90_noerr)) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

   end subroutine read_netcdf_2d_int

   ! ----------------------------------------------------------
   ! Read in 2D arrays Characters
   ! ----------------------------------------------------------

   subroutine read_netcdf_2d_char (nc_file_id, start_var, var_dim, var_name, var_output)
        implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: start_var(:)
      integer, dimension(:), intent(in) :: var_dim

      character(len=*), intent(in) :: var_name
      character(len=30) , intent(out), dimension(:,:) :: var_output
      character(len=30), allocatable, dimension(:,:) :: var

      integer :: nc_var_id
      integer :: status, tmp1, tmp2, i, j
      integer, dimension(2) ::dimIDs

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !find dimentions
      status = nf90_inquire_variable(nc_file_id, nc_var_id, dimids = dimIDs)
      status = nf90_inquire_dimension(nc_file_id, dimIDs(1), len = tmp1)
      status = nf90_inquire_dimension(nc_file_id, dimIDs(2), len = tmp2)
      allocate (var(tmp1,tmp2))

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var, start=(/1,1/), count=(/tmp1,tmp2/) )
      if ((status /= nf90_noerr)) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

      !extract and save classifier names to the final array
      do i = 1, tmp2
        if ((var(i,1) .ge. 'a' .and. var(i,1) .le. 'z') &
        .or.(var(i,1) .ge. 'A' .and. var(i,1) .le. 'Z')) then
           var_output(i,:) = trim(var(i,1))
        endif
      enddo

      if (allocated(var)) deallocate (var)

   end subroutine read_netcdf_2d_char

  ! ----------------------------------------------------------
   ! Read in 3D arrays (used code from DCOMP reader
   ! ----------------------------------------------------------
   subroutine read_netcdf_3d (nc_file_id, start_var, var_dim, var_name, var_output)
         implicit none
      integer, intent(in) :: nc_file_id
      integer, intent(in) :: start_var(:)
      integer, dimension(:), intent(in) :: var_dim

      character(len=30), intent(in) :: var_name
      real, intent(out), dimension(:,:,:) :: var_output

      integer :: nc_var_id
      integer :: status = 0

      status = nf90_inq_varid(nc_file_id, trim(var_name), nc_var_id)
      if (status /= nf90_noerr) then
            print *, "Error: Unable to get variable id for ", trim(var_name)
            return
      endif

      !get Variable
      status = nf90_get_var(nc_file_id, nc_var_id, var_output, start=start_var, count=var_dim)
      if ((status /= nf90_noerr)) THEN
            print *,'Error: ',  trim(nf90_strerror(status)),'   ', trim(var_name)
            return
      ENDIF

   end subroutine read_netcdf_3d
   
!-------------------------------------------------------------------------------------   

end module NB_CLOUD_MASK

