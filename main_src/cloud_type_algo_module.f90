! $Header$
!
!
!  cloud type algorithm
!
!  HISTORY:
!     2014/04/30:  recoded from AK heidinger type code  (AW)
!
!   TODO : H20 profiling
!


module cloud_type_algo_module

   implicit none
   
   private
   public :: cloud_type_pixel
   public :: cloud_type_input_type
   public :: ET_cloud_type
   public :: set_cloud_phase
   private:: COMPUTE_ICE_PROBABILITY_BASED_ON_TEMPERATURE
  
    type  et_cloudtype_class_type
      integer :: FIRST = 0 
      integer :: CLEAR = 0
      integer :: PROB_CLEAR = 1
      integer :: FIRST_WATER = 2
      integer :: FOG = 2
      integer :: WATER = 3
      integer :: SUPERCOOLED = 4
      integer :: LAST_WATER = 4
      integer :: MIXED = 5
      integer :: FIRST_ICE = 6
      integer :: OPAQUE_ICE = 6
      integer :: TICE = 6
      integer :: CIRRUS = 7
      integer :: OVERLAP = 8
      integer :: OVERSHOOTING = 9
      integer :: LAST_ICE = 9
      integer :: UNKNOWN = 10
      integer :: DUST = 11
      integer :: SMOKE = 12
      integer :: FIRE = 13
      integer :: LAST = 13
   end type    
   type ( et_cloudtype_class_type ) :: ET_cloud_type
   
    type  et_cloudphase_class_type
      integer :: FIRST = 0 
      integer :: CLEAR = 0
      integer :: WATER = 1
      integer :: SUPERCOOLED = 2
      integer :: MIXED = 3
      integer :: ICE = 4
      integer :: UNKNOWN = 5
      integer :: LAST = 5
   end type    
   type ( et_cloudphase_class_type )  :: ET_cloud_phase  
   
   
   type cloud_type_sat_type
      logical , dimension(42) :: chan_on
      real :: ref_ch6
      real :: ref_ch20
      real :: rad_ch27
      real :: bt_ch27
      real :: rad_ch31 
      real :: bt_ch31 
      real :: bt_ch32   
   end type cloud_type_sat_type
   
   type cloud_type_rtm_type
      real, allocatable:: rad_ch31_bb_prof (:)
      real, allocatable:: rad_ch27_bb_prof (:)
      real, allocatable:: p_prof (:)
      real, allocatable:: z_prof (:)
      real, allocatable:: t_prof (:)
      integer :: tropo_lev
      integer :: sfc_lev
      real :: bt_ch31_3x3_max
      real :: bt_ch31_3x3_std
      real :: Covar_Ch27_Ch31_5x5
      real :: ref_ch6_clear
      real :: rad_ch27_atm_sfc
      real :: rad_ch31_atm_sfc
      real :: bt_ch31_atm_sfc
      real :: bt_ch32_atm_sfc
      real :: emiss_tropo_ch31     
      real :: Beta_11um_12um_Tropo
      real :: Beta_11um_133um_Tropo
       
   end type cloud_type_rtm_type
   
   type cloud_type_geo_type   
      real :: sol_zen
      real :: sat_zen
   end type cloud_type_geo_type
   
   type cloud_type_sfc_type   
      real :: emiss_ch20
   end type cloud_type_sfc_type   
   
   type cloud_type_input_type
      type ( cloud_type_sat_type) :: sat
      type ( cloud_type_rtm_type) :: rtm
      type ( cloud_type_geo_type) :: geo
      type ( cloud_type_sfc_type) :: sfc
   end type cloud_type_input_type
   
   real, parameter :: MISSING_REAL = -999.

contains
   ! -----------------------------------------------------------------------------
   ! - this computes cloud type without LRC correction
   ! - LRC correction can be done with theoptional force_ice keyword:
   ! - 
   ! -  LRC core is water and pixel is ice  ==>   correct the pixel to WATER
   ! -  LRC core is ice and pixel is water  ==> run this with force_ice = .true.
   ! ----------------------------------------------------------------------------- 
   
   subroutine cloud_type_pixel ( inp , ctype , ice_prob_out , force_ice  ) 
      implicit none
      type ( cloud_type_input_type) , intent(in) :: inp
      integer , intent ( out ) :: ctype 
      real, intent(out), optional :: ice_prob_out
      logical, intent(in), optional :: force_ice ! - this is convinient for LRC correction
      
      real :: ice_prob
      logical ::  force_ice_phase 
      logical :: is_cirrus
      logical :: is_water
      real :: t_opa
      real :: z_opa
      
      logical :: is_overlap
      logical :: experimental_on = .false.
      
      ! --   executable -----------------
      is_cirrus = .false.
      is_water = .false.
      
      force_ice_phase = .false.
      if (present (force_ice)) then
         if ( force_ice ) force_ice_phase = .true.          
      end if
      
      call get_ice_probabibility (inp &
         & , ice_prob , is_cirrus , is_water , t_opa , z_opa ) 
      
      ! - compute type from ice probablity phase discrimination
      if ( ice_prob > 0.5 .or. force_ice_phase ) then
         call determine_type_ice ( inp % Rtm % Emiss_tropo_ch31 &
            , inp % sat % bt_ch31 &
            , inp % rtm % Beta_11um_12um_Tropo &
            , inp % rtm % Beta_11um_133um_Tropo &
            , is_water, is_cirrus , ctype ) 
      else 
         call determine_type_water ( z_opa , t_opa , ctype ) 
      end if
      
      ! - optional output of ice probability
      if ( present ( ice_prob_out ) ) ice_prob_out = ice_prob
      
      ! - experimental overlap test
      
      if ( experimental_on ) then
      
         call overlap_test ( inp % sat % bt_ch31 &
                        , inp % sat % bt_ch32 &
                        , inp % sat % ref_ch6 & 
                        , inp % geo % sol_zen &
                        , inp % geo % sat_zen &
                        , is_overlap ) 
         
         print*,ctype,is_overlap,inp % sat % bt_ch31-inp % sat % bt_ch32,inp % sat % ref_ch6
      end if   
   end subroutine cloud_type_pixel
   
   ! ---------------------------------------------------------------
   !  returns type for ice phase pixels ( ice_probability gt 0.5)
   ! ---------------------------------------------------------------
   subroutine determine_type_ice ( emiss_tropo_11 &
      , bt_11 &
      , beta_11_12 &
      , beta_11_13 &
      , is_water &
      , is_cirrus &
      , ctype ) 
      
      implicit none
      
      real, intent(in) :: emiss_tropo_11
      real, intent(in) :: bt_11
      real, intent(in) :: beta_11_12
      real, intent(in) :: beta_11_13
      logical , intent(in) :: is_water
      logical , intent(in) :: is_cirrus
      integer, intent(out) :: ctype
      real, parameter :: BETA_11UM_12UM_OVERLAP_THRESH = 0.95
      real, parameter :: BETA_11UM_133UM_OVERLAP_THRESH = 0.70
      
           
      ctype = ET_cloud_type % OPAQUE_ICE
         
      if ( Emiss_tropo_11 > 0.95) then
         ctype = ET_cloud_type % OVERSHOOTING
         
      end if
         
      if ( Emiss_tropo_11 < 0.8) then
         ctype = ET_cloud_type % CIRRUS
         
         if ( beta_11_12 > 0. .and. beta_11_12 < BETA_11UM_12UM_OVERLAP_THRESH ) then
            ctype = ET_cloud_type % OVERLAP            
         end if 
         
         if ( beta_11_13 > 0. .and. beta_11_13 < BETA_11UM_133UM_OVERLAP_THRESH ) then
            ctype = ET_cloud_type % OVERLAP            
         end if   
            
         if ( is_cirrus .and. is_water ) ctype = ET_cloud_type % OVERLAP
         
         
         
            
      end if
         
      if (  bt_11 < 233.0 ) then
         ctype = ET_cloud_type % OPAQUE_ICE
      end if 
      
   
   end subroutine determine_type_ice
   
   ! ---------------------------------------------------------------
   !  returns type for water phase pixels ( ice_probability lt 0.5)
   ! ---------------------------------------------------------------
   subroutine determine_type_water ( z_opa, t_opa , ctype )
      real, intent (in ) :: z_opa
      real, intent ( in) :: t_opa
      integer, intent(out) :: ctype
   
      ctype = ET_cloud_type % WATER
         if ( t_opa < 273.0 .and. t_opa /= MISSING_REAL ) then
            ctype = ET_cloud_type % SUPERCOOLED
         else          
            if (z_opa < 1.0 .and. z_opa /= MISSING_REAL ) ctype = ET_cloud_type % FOG
         end if
   
   end subroutine determine_type_water
   
   
   ! -----------------------------------------------------------------------------------
   !  computes ice probability and cirrus water flag and opaque temperature and height
   ! ------------------------------------------------------------------------------------
   subroutine get_ice_probabibility (inp , ice_prob , is_cirrus , is_water , t_opa , z_opa ) 
      implicit none
      type ( cloud_type_input_type) , intent(in) :: inp
      real , intent(out) :: ice_prob
      logical, intent(out) :: is_cirrus
      logical, intent(out) :: is_water
      real, intent(out)    :: t_opa
      real, intent(out)    :: z_opa
     
      
      real, parameter:: ICE_TEMP_MIN = 243.0
      real, parameter:: ICE_TEMP_MAX = 263.0
      real, parameter:: FMFT_COLD_OFFSET = 0.5 !K
      real, parameter:: FMFT_CIRRUS_THRESH = 1.0 !K
      real, parameter :: BT_11UM_STD_CIRRUS_THRESH = 4.0 ! K
      
      real :: fmft
      real :: h2o_correct
         
      is_water = .false.
      is_cirrus = .false.
      

      !------------------------------------------------------------------------
      !  Determine the Cloud Height and Cloud Temperature for Typing
      !------------------------------------------------------------------------
      t_opa = MISSING_REAL
      z_opa = MISSING_REAL
      
      !---- if possible, use water vapor channel
      if ( inp % sat % chan_on ( 27) ) then
         call height_h2o_channel ( &
            inp % sat % rad_ch31 &
            , inp % rtm % rad_ch31_bb_prof &              
            , inp % rtm % rad_ch31_atm_sfc &
            , inp % sat % rad_ch27 &
            , inp % rtm % rad_ch27_bb_prof & 
            , inp % rtm % rad_ch27_atm_sfc &
            , inp % rtm % Covar_Ch27_Ch31_5x5 &
            , inp % rtm % tropo_lev &
            , inp % rtm % sfc_lev &
            , inp % rtm % t_prof &
            , inp % rtm % z_prof &
            , t_opa &
            , z_opa )
      end if
     
      !---- if no water vapor, use the atmos-corrected 11 micron
      if ( t_opa == MISSING_REAL) then
         call height_opaque ( &
            inp % sat % rad_ch31 &
            , inp % rtm % rad_ch31_bb_prof &
            , inp % rtm % tropo_lev &
            , inp % rtm % sfc_lev &
            , inp % rtm % t_prof &
            , inp % rtm % z_prof &
            , t_opa &
            , z_opa )
      end if
      
      !--- last resort, use raw 11 micron brightness temperature
      if ( t_opa == MISSING_REAL) t_opa = inp % sat % bt_ch31    
         

      !--- compute the ice probabilty based on our guess of cloud temperature
      ice_prob = COMPUTE_ICE_PROBABILITY_BASED_ON_TEMPERATURE(t_opa)
      
      
      !------------------------------------------------------------------------
      ! Modify ice_prob based on spectral tests for ice and water
      !------------------------------------------------------------------------

      !--- tests for water
      if ( t_opa > 240.0 ) then
                          
         ! 1.6 spectral test
         if ( inp % sat % chan_on ( 6 ) &
            .and. inp % geo % sol_zen < 80. &
            .and. inp % rtm % ref_ch6_clear < 20. &
            .and. inp % sat % ref_ch6 > 30. )  is_water = .true. 
         
         ! - 3.75 day   
         if ( ice_prob /= 0.0 &
            .and. inp % sat % chan_on ( 20 ) &
            .and. inp % geo % sol_zen < 80. &
            .and. inp % sfc % emiss_ch20 > 0.9 &
            .and. inp % sat % ref_ch20 > 20.0 ) is_water = .true.
         
         
         ! -3.75 night 
         if ( ice_prob /= 0.0 &
            .and. inp % sat % chan_on ( 20 ) &
            .and. inp % geo % sol_zen > 80. &          
            .and. inp % sat % ref_ch20 > 5.0 ) is_water = .true. 
         
          !--- modify ice_prob based on water tests
         if ( is_water) ice_prob = 0.0

      end if
      
      
      !--- tests for ice
      if ( inp % rtm % bt_ch31_3x3_std < BT_11UM_STD_CIRRUS_THRESH ) then
      
         ! - split window test
         if ( inp % sat % chan_on ( 31 ) &
            .and. inp % sat % chan_on ( 32 )) then
                     
            if ( inp % rtm % bt_ch31_atm_sfc <= 265.0 ) then
               h2o_correct = FMFT_COLD_OFFSET
            else 
               h2o_correct = (inp % rtm % bt_ch31_atm_sfc -inp % rtm %  bt_ch32_atm_sfc) &
                  & * (inp % sat % bt_ch31 -260.0) &
                  & / (inp % rtm % bt_ch31_atm_sfc -260.0)
            end if
            
            fmft = inp % sat % bt_ch31 - inp % sat % bt_ch32 - h2o_correct
            
            if ( fmft > FMFT_CIRRUS_THRESH )   is_cirrus = .true.            
         end if
         
         ! - 6.7 covariance test
         if ( inp % sat % chan_on ( 27 ) &
            .and. inp % sat % chan_on ( 31 ) ) then
            
            if ( inp%rtm%Covar_Ch27_Ch31_5x5 > 1.5 &
               .and. inp % rtm % bt_ch31_3x3_max < 250.0) is_cirrus = .true.         
         end if 
         
          !--- modify ice_prob based on cirrus tests
         if ( is_cirrus ) ice_prob = 1.0
      
      end if
      
      ! ++++++++++++++
      
          
   end subroutine get_ice_probabibility
   
   ! -----------------------------------------------
   !   computes Height parameters from NWP profiles
   !    Called in get_ice_probabibility
   ! ------------------------------------------------   
   subroutine height_h2o_channel ( &
           rad_11um &
         , rad_11um_rtm_prof &  
         , rad_11um_clear &
         , rad_h2o &
         , rad_h2o_rtm_prof &
         , rad_h2o_clear &
         , Covar_h2o_11 &
         , tropo_lev &
         , sfc_lev &
         , t_prof &
         , z_prof &
         , t_opa &
         , z_opa )
         
      implicit none      
      
      real, intent(in) :: rad_11um
      real, intent(in) :: rad_11um_clear
      real , intent(in) :: rad_11um_rtm_prof(:)
      real , intent(in) :: rad_h2o_rtm_prof(:)
      real, intent(in) :: rad_h2o
      real, intent(in) :: rad_h2o_clear
      real, intent(in) :: Covar_h2o_11
      integer, intent(in) :: tropo_lev
      integer, intent(in) :: sfc_lev
      real, intent(in) :: t_prof(:)
      real, intent(in) :: z_prof(:)   
      real, intent(out) :: t_opa
      real, intent(out) :: z_opa
      
      real, parameter :: RAD_11UM_THRESH = 2.0
      real, parameter :: RAD_67UM_THRESH = 0.25
      real, parameter :: BT_CH27_CH31_COVAR_CIRRUS_THRESH = 1.0
      
      integer :: cld_lev
      integer :: idx_lev
      
      real :: denominator
      real :: slope
      real :: intercept
      
      real :: rad_H2O_BB_prediction
      
      ! -------------------------------            
      t_opa = MISSING_REAL
      z_opa = MISSING_REAL
            
      ! some tests
      if ( rad_h2o < 0 ) return
      
      if ( rad_h2o_clear - rad_h2o < RAD_67UM_THRESH ) return
      
      if ( rad_11um_clear - rad_11um < RAD_11UM_THRESH ) return      
      if ( Covar_h2o_11 /= MISSING_REAL .and.  Covar_h2o_11 < BT_CH27_CH31_COVAR_CIRRUS_THRESH ) return
      
      ! - colder than tropopause
      if ( rad_11um < rad_11um_rtm_prof ( tropo_lev ) ) then
         cld_lev = tropo_lev
      else
         Denominator =  rad_11um - rad_11um_clear
         
         slope = (Rad_h2o - Rad_H2o_Clear) / (Denominator)
         intercept = Rad_h2o - Slope * Rad_11um
         
         cld_lev = 0
         
         do idx_lev = tropo_lev + 1, sfc_lev
            
            rad_H2O_BB_prediction = slope * rad_11um_rtm_prof (idx_lev) + Intercept
            
            if ( rad_H2O_BB_prediction < 0 ) cycle
           
            if ((rad_H2O_BB_Prediction > rad_h2o_rtm_prof(idx_lev-1)) .and. & 
               ( rad_H2O_BB_Prediction <= rad_h2o_rtm_prof(idx_lev))) then
               cld_lev = idx_lev
               exit
          endif
         
         end do
      
      end if 
      
      
      if (cld_lev > 0 ) then
         t_opa = t_prof( cld_lev )
         z_opa = z_prof( cld_lev )
      end if
      
   end subroutine height_h2o_channel   
   
   ! -----------------------------------------------
   !   computes Height parameters from NWP profiles
   !    Called in get_ice_probabibility
   ! ------------------------------------------------
   subroutine height_opaque ( &
        rad31 &
      , rad31_rtm_prof &
      , tropo_lev &
      , sfc_lev &
      , t_prof &
      , z_prof &
      , t_opa &
      , z_opa )
      
      implicit none
      
      real , intent(in) :: rad31
      real , intent(in) :: rad31_rtm_prof(:)
      integer , intent(in) :: tropo_lev
      integer , intent(in) :: sfc_lev
      real, intent(in) :: t_prof(:)
      real, intent(in) :: z_prof(:)
      real, intent(out) :: t_opa
      real, intent(out) :: z_opa
      integer :: cld_lev
      integer :: idx_lev
      
      t_opa = MISSING_REAL
      z_opa = MISSING_REAL
      cld_lev = 0 
      if ( rad31 < rad31_rtm_prof(tropo_lev) ) then
         cld_lev = tropo_lev
      else
         do idx_lev = tropo_lev , sfc_lev
            if ( rad31 > rad31_rtm_prof ( idx_lev ) ) cld_lev = idx_lev
         end do         
      end if
      
      if (cld_lev > 1 ) then
         t_opa = t_prof ( cld_lev)
         z_opa = z_prof (cld_lev )
      end if
        
   end subroutine height_opaque
   
   ! -----------------------------------------------------
   !   experimental code for additional overlap test
   !      this is just a placeholder now..
   !
   !      we need coeffcieints for vis and bt difference as
   !        function of solar and sensor angle
   !  
   !    05/01/2014 AW
   ! ---------------------------------------------------- 
   subroutine overlap_test ( bt_11 &
                        , bt_12 &
                        , ref_vis & 
                        , sol_zen &
                        ,  sat_zen &
                        , is_overlap )
                        
      real, intent(in) :: bt_11
      real, intent(in) :: bt_12
       real, intent(in) :: ref_vis
      real, intent(in) :: sol_zen
      real, intent(in) :: sat_zen
      
      logical, intent(out) :: is_overlap
      
      real :: btd_thresh
      real :: ref_vis_thresh
      
      is_overlap = .false.
     
       
      btd_thresh = 1. 
      ref_vis_thresh = 40.                 
                        
      if (( bt_11 - bt_12 ) > btd_thresh &
         .and.  ref_vis > ref_vis_thresh ) is_overlap = .true.                 
   end subroutine overlap_test
   
   !-----------------------------------------------------
   !
   !-----------------------------------------------------
    subroutine set_cloud_phase  ( ctype , cphase ) 
      integer, parameter :: int1 = selected_int_kind(1)
      integer(INT1) , intent(in) ::  ctype (:,:)
      integer(INT1) , intent(out) :: cphase (:,:)
      
      cphase = ET_cloud_phase % UNKNOWN
      
      where ( ctype ==  ET_cloud_type % CLEAR &
            .or. ctype ==  ET_cloud_type % PROB_CLEAR )
         cphase = ET_cloud_phase % CLEAR 
      end where
      
      where ( ctype >=  ET_cloud_type % FIRST_WATER &
            .and.  ctype >=  ET_cloud_type % LAST_WATER )
         cphase = ET_cloud_phase % WATER     
      end where
      
      where ( ctype ==  ET_cloud_type % SUPERCOOLED )
         cphase = ET_cloud_phase % SUPERCOOLED  
      end where
      
      where ( ctype ==  ET_cloud_type % MIXED )
         cphase = ET_cloud_phase % MIXED  
      end where
      
      where ( ctype >=  ET_cloud_type % FIRST_ICE &
            .and.  ctype >=  ET_cloud_type % LAST_ICE )
         cphase = ET_cloud_phase % ICE
     
      end where
      
      
   
   end subroutine set_cloud_phase
!====================================================================
! Function Name: Compute_Ice_Probability_Based_On_Temperature
!
! Function: Provide the probability that this pixel is ice
!
! Description:
!    Use the cloud temperature and an assumed relationship to
!    determine the probability that the cloud is ice phase
!
! Dependencies:
!
! Restrictions:
!
! Reference:
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
!====================================================================
function COMPUTE_ICE_PROBABILITY_BASED_ON_TEMPERATURE(T_Opa) result(Ice_Prob)

real:: T_opa
real:: Ice_Prob
real, parameter:: Ice_Temperature_Min = 243.0
real, parameter:: Ice_Temperature_Max = 263.0

Ice_Prob = 1.0 - (T_opa-Ice_Temperature_Min)/(Ice_Temperature_Max - Ice_Temperature_Min)
Ice_Prob = min(1.0,Ice_Prob)
Ice_Prob = max(0.0,Ice_Prob)

end function COMPUTE_ICE_PROBABILITY_BASED_ON_TEMPERATURE

end module cloud_type_algo_module
    
