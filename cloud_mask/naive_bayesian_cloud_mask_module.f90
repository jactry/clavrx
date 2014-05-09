! $Header$
!

module naive_bayesian_cloud_mask_module

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
      integer , parameter :: et_class_R_016_DAY = 18
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
      integer :: glint
      logical :: solar_conta
   end type cloud_mask_geo_type
   
   type cloud_mask_sfc_type
      integer :: land_class
      integer :: coast_mask
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
   
   type cloud_mask_sat_type
      logical , dimension(42) :: chan_on
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
   end type cloud_mask_sat_type

   type cloud_mask_diagnostic
      real :: diagnostic_1
      real :: diagnostic_2
      real :: diagnostic_3
   end type cloud_mask_diagnostic

   type cloud_mask_input_type
      character ( len =256) :: bayesian_mask_classifier      
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
      real, allocatable :: prior_no (:)
      real, allocatable :: prior_yes (:)
      real, allocatable :: Optimal_Posterior_Prob (:)
      
      real, allocatable :: Classifier_Bounds_Min(:,:)
      real, allocatable :: Classifier_Bounds_Max(:,:)
      real, allocatable :: Delta_Classifier_Bounds(:,:)
     
      real, allocatable :: class_cond_no(:,:,:)
      real, allocatable :: class_cond_yes(:,:,:) 
      real, allocatable :: cond_no(:,:,:)
      real, allocatable :: cond_yes(:,:,:)
      character (len=20), allocatable :: Classifier_Value_Name(:)
      integer, allocatable :: Classifier_Value_Name_enum(:)
      integer, allocatable :: flag_idx(:)
      
     
      
   
   contains
      procedure :: alloc => alloc_bayes_coef
      procedure :: dealloc => dealloc_bayes_coef
   end type bayes_coef_type
   
   type ( bayes_coef_type) , private , save :: bayes_coef
  
contains
   !
   !
   !
   subroutine cloud_mask_naive_bayes ( inp , erg , info_flags , diag )
          
      implicit none
            
      type ( cloud_mask_input_type ) , intent ( in ) :: inp
      type ( cloud_mask_diagnostic ) , intent ( inout ) :: diag
      real , intent ( out ) :: erg
      integer , intent ( out ) , optional :: info_flags ( 7 )
      
      integer :: sfc_type_number 
      integer :: class_idx, sfc_idx
      real :: Classifier_Value
      real, allocatable :: Cond_yes(:)
      real, allocatable :: Cond_no(:)
      
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
            
      integer :: pos_info_flag 
      integer :: idx_info_flag
      real :: class_contr
      
      ! ---    Executable  ----------------------------
      ! - read in classifer coeffiecients
      bayes_coef % file =trim(inp  % bayesian_mask_classifier ) 
      if ( .not. bayes_coef % is_read) call read_bayes_coeff ( ) 
      
      ! - determine sfc type
      
      sfc_type_number =  bayes_sfc_type ( inp% geo % lat , inp % geo % lat &
         & , inp % sfc % land_class , inp % sfc % coast_mask, inp % sfc % snow_class , inp % sfc % sfc_type &
         & , inp % sfc % emis_ch20,  inp % sfc % sst_anal_uni )
          
         
      allocate ( Cond_Yes ( bayes_coef % n_class ) )
      allocate ( Cond_No ( bayes_coef % n_class ) )
            
          
      sfc_idx = sfc_type_number
           
            ! - several 0/1 flags
      is_mountain =  inp % sfc % dem  > 2000.0 &
                            & .and. sfc_idx /= 6
                            
      use_lunar_refl_for_vis_tests = .false.
      if ( inp % sat % chan_on(42)) then
         if ( inp % sat % ref_dnb_lunar >= 0. .and. &
          ( inp % geo %  scat_angle_lunar  > 80. .or. inp % geo % lunar_zen  > 95. ) .and. &
          .not. is_mountain .and. &
          .not. inp % sfc % coast_mask .and. &
          .not. inp % sfc % snow_class  == ET_snow_class % SNOW )  then
            use_lunar_refl_for_vis_tests  = .true.      
         end if    
      end if                       
            
      has_cold_btd = .false.
      if ( inp % sat % chan_on(31) ) then
         if (  inp%sat %bt_ch31  < BT_11UM_COLD_SCENE_THRESH .and. &
                     & inp%sat %bt_ch31 /= MISSING_VALUE_REAL4 ) then 
                  has_cold_btd = .true. 
         end if         
      end if
            
      is_cold_375um = .false.
      if ( inp % sat % chan_on(20) ) then
         if (  inp%sat %bt_ch20 < BT_037UM_COLD_SCENE_THRESH .and. &
                     & inp%sat %bt_ch20 /= MISSING_VALUE_REAL4 ) then 
            is_cold_375um = .true. 
         end if         
      end if
            
      is_day_375um = .true.
      if (( inp % geo % sol_zen  >  SOLZEN_375UM_DAY_THRESH ) .or. &
                  &  (inp % geo %airmass > AIRMASS_THRESH  )) then           
         is_day_375um = .false.
      end if   
            
      is_night_375um = .true.
      if ( inp % geo % sol_zen  <  SOLZEN_375UM_NIGHT_THRESH )  then            
               is_night_375um = .false.
      end if   
               
      is_day_063um = .true.
      info_flags(1) = ibset ( info_flags(1) , 1 )
      if (( inp % geo % sol_zen >  SOLZEN_063UM_DAY_THRESH ) .or. &
                  &  (inp % geo %airmass > AIRMASS_THRESH  )) then           
         is_day_063um = .false.
         
      end if 
      
      is_day_063um_spatial_tests = .true.
      info_flags(1) = ibset ( info_flags(1) , 1 )
      if ( inp % geo % sol_zen >  SOLZEN_063UM_DAY_THRESH_SPATIAL_TESTS ) then           
         is_day_063um_spatial_tests = .false.
         
      end if
      
            
      is_forward_scatter = .false.
      if ( inp % geo %  scat_angle  < 80. .and. inp % geo % sol_zen  < 95.) then
         is_forward_scatter = .true.
      end if 
      
      is_smoke = .false.
      ! - TO ADD 
          
      is_dust = .false.
      if ( inp % sat % chan_on (1) .and. inp % sat % chan_on (8) ) then
         is_dust = dust_detection ( &
              inp % sat % ref_ch1 &
            , inp % sat % ref_ch8 &
            , inp % sat % ref_ch1_3x3_std &
            , inp % geo % glint &
            , inp % sfc % land_class ) 
      end if
      
      is_smoke = .false.
      if ( inp % sat % chan_on (1) .and. inp % sat % chan_on (7) ) then
         is_dust = smoke_detection ( &
              inp % sat % ref_ch1 &
            , inp % sat % ref_ch7 &
            , inp % sat % ref_ch1_3x3_std &
            , inp % geo % sat_zen &
            , inp % geo % glint &
            , inp % sfc % land_class ) 
      end if
            
      is_cloud_shadow = .false.
      ! - TO ADD
      
      is_fire = .false.      
      if ( inp % sat % chan_on (31) .and. inp % sat % chan_on (20) ) then
         is_fire = fire_detection ( &
              inp % sat % bt_ch31 &
            , inp % sat % bt_ch20 &
            , inp % rtm % bt_ch31_3x3_std &
            , inp % rtm % bt_ch20_3x3_std &
            , inp % geo % sol_zen )
      end if
      
      is_solar_contaminated = inp % geo % solar_conta
      
    
      
      info_flags = 0 
      
                                    info_flags(1) = ibset ( info_flags ( 1 ) , 0 )
      if ( is_day_063um)            info_flags(1) = ibset ( info_flags ( 1 ) , 1 )
      if ( is_day_063um_spatial_tests)  info_flags(1) = ibset ( info_flags ( 1 ) , 2 )
      if ( is_day_375um)            info_flags(1) = ibset ( info_flags ( 1 ) , 3 )
      if ( is_night_375um)          info_flags(1) = ibset ( info_flags ( 1 ) , 4 )
      if ( is_solar_contaminated)   info_flags(1) = ibset ( info_flags ( 1 ) , 5 )
      if ( inp % sfc % coast_mask)  info_flags(1) = ibset ( info_flags ( 1 ) , 6 )    
      if ( is_mountain )            info_flags(1) = ibset ( info_flags ( 1 ) , 7 )
      
      if ( is_forward_scatter )     info_flags(2) = ibset ( info_flags ( 2 ) , 0 )
      if ( is_cold_375um )          info_flags(2) = ibset ( info_flags ( 2 ) , 1 )
      if ( has_cold_btd )           info_flags(2) = ibset ( info_flags ( 2 ) , 2 )
      if ( is_smoke)                info_flags(2) = ibset ( info_flags ( 2 ) , 4 )
      if ( is_dust)                 info_flags(2) = ibset ( info_flags ( 2 ) , 5 )
      if ( is_cloud_shadow)         info_flags(2) = ibset ( info_flags ( 2 ) , 6 )
      if ( is_fire)                 info_flags(2) = ibset ( info_flags ( 2 ) , 7 )
      
      
      ! -TODO : probably wrong
      info_flags(3) = ibits ( sfc_idx , 0 , 3 )
      
                         
      ! - class loop 
      class_loop: do class_idx= 1,   bayes_coef % n_class
         
         ! - init
         is_on_test = .false.        
    
         cond_yes (class_idx)  = 1.0
         cond_no  (class_idx)  = 1.0  
        
         select case (  bayes_coef % Classifier_Value_Name_enum (class_idx))
          
         case( et_class_T110 )
            
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            Classifier_Value = inp % sat % bt_ch31  
            if (  inp%sfc %  land_class  == ET_land_class % DEEP_OCEAN ) then
               Classifier_Value = inp%rtm % bt_ch31_lrc 
            end if
            is_on_test = .true.
            pos_info_flag = 4
            idx_info_flag = 3
                             
         case( et_class_TMAX_T)
            if ( .not. inp%sat % chan_on(31) ) cycle class_loop
            if ( is_mountain  ) cycle class_loop
            if ( inp % sfc % coast_mask   ) cycle
            Classifier_Value = inp % rtm % bt_ch31_3x3_max  &
                               & -  inp % sat % bt_ch31 
            is_on_test = .true.
            pos_info_flag = 6
            idx_info_flag = 3
           
         case (et_class_T_STD)
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( is_mountain  ) cycle class_loop
            if ( inp % sfc % coast_mask ) cycle
            Classifier_Value = inp % rtm % bt_ch31_3x3_std 
            is_on_test = .true.
            pos_info_flag = 0
            idx_info_flag = 4
           
         case (et_class_E_TROP)
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
              
            Classifier_Value = inp % rtm % emis_ch31_tropo 
            
            is_on_test = .true.
            pos_info_flag = 2
            idx_info_flag = 4
           
         case  (et_class_FMFT)
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( .not. inp % sat % chan_on(32) ) cycle class_loop
            if ( has_cold_btd ) cycle class_loop
           
            Classifier_Value = split_window_test( inp % rtm % bt_ch31_atm_sfc  &
                                                  , inp % rtm % bt_ch32_atm_sfc  &
                                                  , inp % sat % bt_ch31 &
                                                  , inp % sat % bt_ch32 )
                                                 
                                    
            is_on_test = .true. 
            pos_info_flag = 4
            idx_info_flag = 4
            
         case (et_class_BTD_110_067)
            if ( .not. inp % sat % chan_on(27) ) cycle class_loop
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( has_cold_btd ) cycle class_loop
           
            Classifier_Value = inp % sat % bt_ch31  - inp % sat % bt_ch27
            
            is_on_test = .true.
            pos_info_flag = 6
            idx_info_flag = 4
        
         case ( et_class_BTD_110_067_COV)
            if ( .not. inp % sat % chan_on(27) ) cycle class_loop
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( has_cold_btd ) cycle class_loop
            Classifier_Value = inp % rtm % bt_ch31_ch27_covar
            is_on_test = .true.
            pos_info_flag = 0
            idx_info_flag = 5
           
         case( et_class_BTD_110_085)
            if ( .not. inp % sat % chan_on(29) ) cycle class_loop
            if ( .not. inp % sat % chan_on(31) ) cycle class_loop
            if ( has_cold_btd ) cycle class_loop
          
            Classifier_Value = inp % sat % bt_ch31  - inp % sat % bt_ch29
            is_on_test = .true.
            pos_info_flag = 2
            idx_info_flag = 5
           
         case( et_class_E_037 )
            if ( .not. inp %sat % chan_on(20) ) cycle class_loop                  
            if ( inp % geo % glint )  cycle class_loop                 
            if ( inp % sat % bt_ch20  <= 0 ) cycle class_loop
            if ( is_cold_375um ) cycle class_loop
           
            classifier_value = emiss_375um_test( &
                                  &  inp % sat % emis_ch20_3x3_mean  &
                                  & , inp % rtm % emis_ch20_clear  , is_night = .false.)
            is_on_test = .true.
            pos_info_flag = 4
            idx_info_flag = 5
           
         case( et_class_E_037_DAY)
            if ( .not. inp % sat % chan_on(20) ) cycle class_loop
            if ( is_solar_contaminated) cycle class_loop 
            if ( inp % geo % glint  )  cycle class_loop
            if ( .not. is_day_375um ) cycle class_loop                  
            if ( is_cold_375um ) cycle class_loop
            if ( inp % sat % bt_ch20   <= 0 ) cycle class_loop
          
            classifier_value = emiss_375um_test( &
                                  &   inp % sat % emis_ch20_3x3_mean  &
                                  & , inp % rtm % emis_ch20_clear , is_night =.false.)
            is_on_test = .true.
            pos_info_flag = 6
            idx_info_flag = 5
           
         case( et_class_E_037_NGT)
            if ( .not. inp % sat % chan_on(20) ) cycle
            if ( is_solar_contaminated) cycle class_loop 
            if ( .not. is_night_375um ) cycle  
            if (is_cold_375um ) cycle class_loop 
            if (  inp % sat % bt_ch20  <= 0 ) cycle class_loop
             
            classifier_value = emiss_375um_test( &
                                  & inp % sat % emis_ch20_3x3_mean  &
                                  & , inp % rtm % emis_ch20_clear , is_night =.true.)
            
            is_on_test = .true.
            pos_info_flag = 0
            idx_info_flag = 6
              
        case( et_class_BTD_037_110_NGT)
           if ( .not. inp % sat % chan_on(20) ) cycle
           if ( .not. inp % sat % chan_on(31) ) cycle
           if ( is_solar_contaminated) cycle class_loop 
           if ( .not. is_night_375um ) cycle  
           if (is_cold_375um ) cycle class_loop 
           if (  inp % sat % bt_ch20  <= 0 ) cycle class_loop
           
           classifier_value =  inp % sat % bt_ch20  -  inp % sat % bt_ch31
           
           is_on_test = .true.
           pos_info_flag = 2
           idx_info_flag = 6
              
        case( et_class_R_006_DAY)
        
        
            
            
            ! - this solar test can be also applied for lunar
            !TODO - make own lunar visible coefficients
            !
            if ( use_lunar_refl_for_vis_tests ) then
               Classifier_Value = reflectance_gross_contrast_test( &
                       &  inp % rtm %  ref_dnb_clear &
                       & ,   inp % sat % ref_dnb_lunar )
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
                            
               Classifier_Value = reflectance_gross_contrast_test( &
                       &  inp % rtm %  ref_ch1_clear &
                       & ,   inp % sat % ref_ch1 )
               is_on_test = .true.
               pos_info_flag = 4
               idx_info_flag = 6   
            
            end if
              
         case( et_class_R_006_STD)
            
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
            
         case( et_class_R_006_MIN_3x3_DAY )
         
            if ( use_lunar_refl_for_vis_tests ) then
               Classifier_Value = relative_visible_contrast_test ( &
                       &  inp % sat % ref_dnb_3x3_min &
                       & , inp % sat %ref_dnb_lunar  ) 
                                         
               is_on_test = .true.
               pos_info_flag = 0
               idx_info_flag = 7 
            else
                  
               if ( .not. inp % sat % chan_on(1) ) cycle class_loop
               if ( .not. is_day_063um_spatial_tests ) cycle
               if ( is_mountain  ) cycle class_loop
               if ( inp % sfc % coast_mask   ) cycle
          
               Classifier_Value = relative_visible_contrast_test ( &
                       &  inp % sat % ref_ch1_3x3_min &
                       & , inp % sat %ref_ch1  ) 
                                         
               is_on_test = .true.
               pos_info_flag = 0
               idx_info_flag = 7    
           end if
           
        case( et_class_R_RATIO_DAY)
           if ( .not. inp %  sat % chan_on(1) ) cycle class_loop
           if ( .not. inp % sat % chan_on(2) ) cycle class_loop
           if ( .not. is_day_063um_spatial_tests ) cycle
           if ( is_mountain  ) cycle class_loop
           if ( inp % geo  % glint  )  cycle class_loop
          
           Classifier_Value = reflectance_ratio_test (&
                                & inp % sat % ref_ch1 &
                                , inp % sat % ref_ch2 ) 
           is_on_test = .true.
           pos_info_flag = 2
           idx_info_flag = 7  
              
         case( et_class_R_013_DAY)
            if ( .not. inp %  sat % chan_on(26) ) cycle class_loop
            if ( is_forward_scatter )  cycle class_loop
            if ( .not. is_day_063um )  cycle class_loop
            if ( is_mountain  ) cycle class_loop
            if ( inp % sat % ref_ch26  <= 0 ) cycle class_loop
           
            Classifier_Value = inp % sat % ref_ch26
            is_on_test = .true.   
            pos_info_flag = 4
            idx_info_flag = 7    
            
         case ( et_class_R_016_Day)
            if ( .not. inp % sat % chan_on(6) ) cycle class_loop
            if ( is_forward_scatter )  cycle class_loop
            if ( .not. is_day_063um ) cycle class_loop
            if ( inp % geo % glint  )  cycle class_loop
            if ( inp % sat % ref_ch6 <= 0 ) cycle class_loop
           
            Classifier_Value = inp % sat % ref_ch6
            is_on_test=.true.
            pos_info_flag = 6
            idx_info_flag = 7  
              
         case default
            print*,'unknown class ', bayes_coef % Classifier_Value_Name_enum (class_idx) , &
                       &  bayes_coef % Classifier_Value_Name (class_idx)
         end select
       
       ! --- all tests classifer
        if ( is_on_test) then
            bin_idx = int (( classifier_value &
               &  -  bayes_coef % Classifier_Bounds_Min(Class_Idx,Sfc_Idx))  /    &
               &     bayes_coef % Delta_Classifier_Bounds(Class_Idx,Sfc_Idx)) + 1
             
            bin_idx =  max(1,min(bayes_coef % N_bounds-1,Bin_Idx)) 
            Cond_yes (class_idx) = bayes_coef % class_cond_yes ( bin_idx, class_idx , sfc_idx )
            Cond_no (class_idx)  = bayes_coef % class_cond_no ( bin_idx, class_idx , sfc_idx ) 
            
            class_contr =  bayes_coef % class_cond_yes ( bin_idx, class_idx , sfc_idx ) / &
                           & ( bayes_coef % class_cond_yes ( bin_idx, class_idx , sfc_idx )  &
                           & + bayes_coef % class_cond_no ( bin_idx, class_idx , sfc_idx ) )
                            
    
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
         !select case (  bayes_coef % Classifier_Value_Name_enum (class_idx))
         !  case ( et_class_T110 )
         !  case ( et_class_TMAX_T )
         !  case ( et_class_T_STD )
         !  case ( et_class_E_TROP )
         !  case ( et_class_FMFT )
         !  case ( et_class_BTD_110_067 )
         !  case ( et_class_BTD_110_067_COV )
         !  case ( et_class_BTD_110_085 )
         !  case ( et_class_E_037 )
         !  case ( et_class_E_037_DAY )
         !  case ( et_class_E_037_NGT )
         !  case ( et_class_BTD_037_110_NGT )
         !  case ( et_class_R_006_DAY )
         !  case ( et_class_R_006_STD )
         !  case ( et_class_R_006_MIN_3x3_DAY )
         !  case ( et_class_R_RATIO_DAY )
         !  case ( et_class_R_013_DAY )
         !  case ( et_class_R_016_Day )
         !     diag % diagnostic_1 = Classifier_Value  
         !     diag % diagnostic_2 = class_contr
         !     diag % diagnostic_3 = info_flags ( idx_info_flag )
         !end select
         !
                
      end do class_loop
                        
            ! - compute Posterior_Cld_Probability
           
            erg  = &
                  (bayes_coef % Prior_Yes(Sfc_Idx) * product(Cond_Yes)) / &
                  (bayes_coef % Prior_Yes(Sfc_Idx) * product(Cond_Yes) +  &        
                  bayes_coef % Prior_No(Sfc_Idx) * product(Cond_No))
                  
            deallocate ( Cond_Yes )
            deallocate ( Cond_No )
            
            
           
          
   end subroutine cloud_mask_naive_bayes
   

   !-------------------------------------------------------------------------------------
   !   Reads the coefficients from file
   !-------------------------------------------------------------------------------------
   subroutine read_bayes_coeff ( )
      use file_tools , only : getlun
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
      
      print*,'read coeffs ..'
      print*,trim(bayes_coef % file)
      
      lun = getlun()
      open ( unit=lun, file =  trim(bayes_coef % file) , &
            action ="read", form="formatted",status="old",iostat=ios)
      if (ios/= 0) then
         print*, 'no bayesain cloud mask ',trim(bayes_coef % file)
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
      
      read(unit=lun,fmt=*)
      read(unit=lun,fmt=*) time_diff_max, n_class, n_bounds, n_sfc_bayes
      bayes_coef % n_class = n_class
      bayes_coef % n_bounds = n_bounds

      call bayes_coef % alloc ( n_class, n_bounds, n_sfc_bayes )
      
      do i_sfc = 1 , n_sfc_bayes
         read (unit=lun,fmt=*,iostat=ios) i_sfc_file, Header
         read(unit=lun,fmt=*,iostat=ios) bayes_coef % Prior_Yes(i_sfc), bayes_coef % Prior_No(i_sfc)
         read(unit=lun,fmt=*,iostat=ios) bayes_coef % Optimal_Posterior_Prob(i_Sfc)
         
         do i_class = 1 , n_class
            
            read(unit=lun,fmt=*,iostat=ios)  &
                            int_dummy &
                          ,  int_dummy &
                          ,  int_dummy  &
                          ,  bayes_coef %Classifier_Value_Name(i_class)
           
            
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'T_11' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) = et_class_T110
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'T_max-T' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_TMAX_T
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'T_std' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_T_STD
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Emiss_tropo' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) = et_class_E_TROP
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'FMFT' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) = et_class_FMFT
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Btd_11_67' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_BTD_110_067
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Bt_11_67_Covar' ) &
               &  bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_BTD_110_067_COV
                          
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Btd_11_85' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_BTD_110_085
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Emiss_375' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_E_037
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Emiss_375_Day' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_E_037_DAY
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Emiss_375_Night' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_E_037_NGT
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Btd_375_11_Night' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_BTD_037_110_NGT
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Ref_063_Day' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_R_006_DAY
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Ref_std' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) = et_class_R_006_STD
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Ref_063_Min_3x3_Day' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) = et_class_R_006_MIN_3x3_DAY
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Ref_Ratio_Day' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) = et_class_R_RATIO_DAY
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Ref_138_Day' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_R_013_DAY
               
            if ( trim(bayes_coef %Classifier_Value_Name(i_class)) == 'Ref_160_Day' ) &
               & bayes_coef %Classifier_Value_Name_enum(i_class) =  et_class_R_016_DAY
                        
            read(unit=lun,fmt=*,iostat=ios)  &
                                bayes_coef % Classifier_Bounds_Min(i_class,i_sfc)  &
                             ,   bayes_coef % Classifier_Bounds_Max(i_class,i_sfc)  &
                             ,   bayes_coef % Delta_Classifier_Bounds(i_class,i_sfc)
                                
             read(unit=lun,fmt=*,iostat=ios)  bayes_coef %Class_Cond_Yes(:,i_class,i_sfc) 
             read(unit=lun,fmt=*,iostat=ios)  bayes_coef %Class_Cond_No(:,i_class,i_sfc)
             
             
             
         end do
      
      end do
      
         
      bayes_coef % is_read = .true.
	   	   
   
   end subroutine read_bayes_coeff
   
   !  7 sfc bayesian surface  types 
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   function bayes_sfc_type ( lat , lon, land_class , coast , snow_class , sfc_type , emis_ch20 , sst_anal_uni)
     
       
      real  :: lat ,lon, emis_ch20 , sst_anal_uni
      integer :: land_class,coast,snow_class , sfc_type
      
      
      integer :: bayes_sfc_type
      integer , parameter :: CLOSED_SHRUBS_SFC = 8
      integer , parameter :: OPEN_SHRUBS_SFC = 9
      integer , parameter :: GRASSES_SFC = 10
      integer , parameter :: BARE_SFC = 12
     
      bayes_sfc_type = 0
     
      ! # deep ocean 
      if ( land_class == ET_land_class % DEEP_OCEAN) then
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
               coast  == 1  )  then
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
      if (   lat >= 60.0 &
               & .and. lon > -75. &
               & .and. lon < -10.0  &    
               & .and. (  land_class == ET_land_class % LAND &
                .or.  land_class  == ET_land_class % COASTLINE ) &
               & .and.  snow_class  == ET_snow_class % SNOW ) then
               
          bayes_sfc_type = 6         
      end if     
         
      !-TODO
      
      ! -- # 7 Desert
      
      if ( emis_ch20 < 0.90 .and. abs ( lat ) < 60.0 &
         .and. ( sfc_type == OPEN_SHRUBS_SFC .or. sfc_type == BARE_SFC ) ) then
      
         bayes_sfc_type = 7
      end if   
      
      if ( bayes_sfc_type == 3 .and. &
         & emis_ch20 < 0.93 .and. &
         & abs ( lat ) < 60.0 .and. &
         & ( sfc_type == OPEN_SHRUBS_SFC .or. &
            & sfc_type == CLOSED_SHRUBS_SFC .or. &
            & sfc_type == GRASSES_SFC .or. &
            &  sfc_type == BARE_SFC ) &
            ) then
            
            bayes_sfc_type = 7
      end if
      
      
       
     
   end function bayes_sfc_type
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   subroutine alloc_bayes_coef ( this, n_class , n_bounds , n_sfc_bayes )
      class ( bayes_coef_type ) :: this
      integer :: n_class , n_bounds , n_sfc_bayes
      
      allocate(this%Prior_Yes(N_sfc_bayes))
      allocate(this%Prior_No(N_sfc_bayes))
      allocate(this%Optimal_Posterior_Prob(N_sfc_bayes))
      allocate(this%Classifier_Bounds_Min(N_class,N_sfc_bayes))
      allocate(this%Classifier_Bounds_Max(N_class,N_sfc_bayes))
      allocate(this%Delta_Classifier_Bounds(N_class,N_sfc_bayes))
      allocate(this%Class_Cond_Yes(N_bounds-1,N_class,N_sfc_bayes))
      allocate(this%Class_Cond_No(N_bounds-1,N_class,N_sfc_bayes))
      allocate(this%Classifier_Value_Name(N_class))
      allocate(this%Classifier_Value_Name_enum(N_class))
      allocate(this%Flag_Idx(N_class))
            
   end subroutine alloc_bayes_coef
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   subroutine dealloc_bayes_coef ( this )
      class ( bayes_coef_type ) :: this
           
      deallocate(this%Prior_Yes)
      deallocate(this%Prior_No)
      deallocate(this%Optimal_Posterior_Prob)
      
      deallocate(this%Classifier_Bounds_Min)
      deallocate(this%Classifier_Bounds_Max)
      deallocate(this%Delta_Classifier_Bounds)
      deallocate(this%Class_Cond_Yes)
      deallocate(this%Class_Cond_No)
      deallocate(this%Classifier_Value_Name)
      deallocate(this%Classifier_Value_Name_enum)
      deallocate(this%Flag_Idx)
           
   end subroutine dealloc_bayes_coef
   
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   real elemental function reflectance_ratio_test(ref_vis,ref_nir)
      real, intent(in):: ref_vis
      real, intent(in):: ref_nir
      real, parameter :: MISSING_VALUE_REAL4 = -999.
      
      reflectance_ratio_test = Missing_Value_Real4
      if (ref_vis /= Missing_Value_Real4 .and. ref_nir /= Missing_Value_Real4) then
           reflectance_ratio_test = ref_nir / ref_vis
      end if

   end function reflectance_ratio_test
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------   
   real elemental function relative_visible_contrast_test(ref_min,ref)
      real, intent(in):: ref_min
      real, intent(in):: ref
      real, parameter :: MISSING_VALUE_REAL4 = -999.
      
      relative_visible_contrast_test = Missing_Value_Real4
      if (ref_min /= Missing_Value_Real4) then
           relative_visible_contrast_test = ref - ref_min
      end if

   end function relative_visible_contrast_test
   
   
   !-------------------------------------------------------------------------------------
   !
   !------------------------------------------------------------------------------------- 
   real elemental function reflectance_gross_contrast_test(ref_clear,ref)
      real, intent(in):: ref_clear
      real, intent(in):: ref
      real, parameter :: MISSING_VALUE_REAL4 = -999.
      reflectance_gross_contrast_test = Missing_Value_Real4
      if (ref_clear /= Missing_Value_Real4) then
           reflectance_gross_contrast_test = ref - ref_clear
      endif

   end function reflectance_gross_contrast_test
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   real elemental function split_window_test( t11_clear, t12_clear, t11, t12)

      real, intent(in):: t11_clear
      real, intent(in):: t12_clear 
      real, intent(in):: t11
      real, intent(in):: t12

      split_window_test  = (t11_clear - t12_clear) * (t11 - 260.0) / (t11_clear - 260.0)

      if (t11_clear <= 265.0) split_window_test = 0.0

      split_window_test = (t11 - t12) - split_window_test

   end function split_window_test
  
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   real elemental function emiss_375um_test ( ems, ems_clear , is_night  )
      real, intent(in):: ems_clear
      real, intent(in):: ems
      logical ,intent(in) :: is_night
      
      real, parameter :: MISSING_VALUE_REAL4 = -999.

      emiss_375um_test = MISSING_VALUE_REAL4
      if (ems /= Missing_Value_Real4 .and. ems_clear /= Missing_Value_Real4) then
           emiss_375um_test = (ems- ems_clear) / ems_clear
      end if
      if ( is_night) emiss_375um_test = ems
  
  
   end function
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------   
   logical elemental function fire_detection ( &
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
            ( bt_110_std            < Stddev_11um_Eumet_Fire_Thresh)  &
               ) then
         Fire_detection = .true.
       endif
   
   
   end function fire_detection
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   logical elemental function  smoke_detection ( &
           ref_004 &
         , ref_021 &
         , ref_006_std &
         , sat_zen &
         , glint &
         , land_class )
               
      real , intent(in) :: ref_004
      real , intent(in) :: ref_021
      real , intent(in) :: ref_006_std
      real , intent(in) :: sat_zen
      integer , intent(in) :: glint
      integer , intent(in) :: land_class  
      
      logical :: is_water_sfc      
      logical :: is_glint
      
      
      real, parameter :: SMOKE_ST_DEV_LAND_THRESH = 0.2
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
         
      is_glint = glint == 1
      
      ref_ratio = ref_021 / ref_004 
      
      if ( is_water_sfc .and. ref_ratio > ( SMOKE_CAND_M11M1_REF_RATIO_THRESH * cosd (sat_zen) ) ) return
      
      
      if ( is_glint .and. is_water_sfc )                 st_dev_thresh_cand = SMOKE_ST_DEV_WATER_GLINT_THRESH
      if ( is_glint .and. .not. (is_water_sfc) )         st_dev_thresh_cand = SMOKE_ST_DEV_LAND_GLINT_THRESH
      if ( .not. ( is_glint) .and. is_water_sfc )        st_dev_thresh_cand = SMOKE_ST_DEV_WATER_THRESH
      if ( .not. (is_glint) .and. .not. (is_water_sfc) ) st_dev_thresh_cand = SMOKE_ST_DEV_LAND_THRESH
      
      if (  ref_006_std < st_dev_thresh_cand ) smoke_detection = .true.
      
   end function smoke_detection
   
   !-------------------------------------------------------------------------------------
   !
   !-------------------------------------------------------------------------------------
   logical elemental function dust_detection ( &
              ref_004 &
            , ref_006 &
            , ref_006_std &
            , glint &
            , land_class )
            
      real , intent(in) :: ref_004
      real , intent(in) :: ref_006
      real , intent(in) :: ref_006_std
      integer , intent(in) :: glint
      integer , intent(in) :: land_class  
      
      logical :: is_water_sfc      
      logical :: is_glint
              
      real, parameter :: dust_m1_refl_thresh = 0.1
      real, parameter :: dust_cand_m1m5_refl_ratio_thresh = 0.25

      real, parameter :: DUST_ST_DEV_LAND_THRESH = 0.2
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
      
      is_glint = glint == 1
            
      ! adjust thresholds
      
      if ( is_glint .and. is_water_sfc )                 st_dev_thresh_cand = DUST_ST_DEV_WATER_GLINT_THRESH
      if ( is_glint .and. .not. (is_water_sfc) )         st_dev_thresh_cand = DUST_ST_DEV_LAND_GLINT_THRESH
      if ( .not. ( is_glint) .and. is_water_sfc )        st_dev_thresh_cand = DUST_ST_DEV_WATER_THRESH
      if ( .not. (is_glint) .and. .not. (is_water_sfc) ) st_dev_thresh_cand = DUST_ST_DEV_LAND_THRESH
      
      if (  ref_006_std < st_dev_thresh_cand ) dust_detection = .true.
   
   
   end function dust_detection
   
   
   
end module naive_bayesian_cloud_mask_module
