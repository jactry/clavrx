! $Header:$
!
!
!
module cloud_phase_algo_module
   
  
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
      integer :: OPAQUE_ICE = 6
      integer :: TICE = 6
      integer :: CIRRUS = 7
      integer :: OVERLAP = 8
      integer :: OVERSHOOTING = 9
      integer :: UNKNOWN = 10
      integer :: DUST = 11
      integer :: SMOKE = 12
      integer :: FIRE = 13
      integer :: LAST = 13
   end type    
   type ( et_cloudtype_class_type ) :: ET_cloud_type
   
   type cloud_type_sat_type
      logical , dimension(42) :: chan_on
      real :: rad_ch31 
      real :: ref_ch6 
      real :: ref_ch20
      real :: bt_ch31 
      real :: bt_ch32   
   end type cloud_type_sat_type
   
   type cloud_type_rtm_type
      real, allocatable:: rad_ch31_bb_prof (:)
      real, allocatable:: p_prof (:)
      real, allocatable:: z_prof (:)
      real, allocatable:: t_prof (:)
      integer :: tropo_lev
      integer :: sfc_lev
      real :: bt_ch31_3x3_max
      real :: bt_ch31_3x3_std
      real :: Covar_Ch27_Ch31_5x5
      real :: ref_ch6_clear
      real :: bt_ch31_atm_sfc
      real :: bt_ch32_atm_sfc
      real :: emiss_tropo_ch31     
      real :: Beta_11um_12um_Tropo
      real :: Beta_11um_133um_Tropo
       
   end type cloud_type_rtm_type
   
   type cloud_type_geo_type   
      real :: sol_zen
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

contains

   ! - this computes cloud type without any LRC correction
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
      
      is_cirrus = .false.
      is_water = .false.
      
      force_ice_phase = .false.
      if (present (force_ice)) then
         if ( force_ice ) force_ice_phase = .true.          
      end if
      
      call get_ice_probabibility (inp , ice_prob , is_cirrus , is_water , t_opa , z_opa )
      
      if ( ice_prob > 0.5 .or. force_ice_phase ) then
         call determine_type_ice ( inp % Rtm % Emiss_tropo_ch31 &
            ,inp % sat % bt_ch31 &
            , inp % rtm % Beta_11um_12um_Tropo &
            , inp % rtm % Beta_11um_133um_Tropo &
            , is_water, is_cirrus , ctype ) 
      else 
         call determine_type_water ( z_opa , t_opa , ctype ) 
      end if
      
      if ( present ( ice_prob_out ) ) ice_prob_out = ice_prob
   
   end subroutine cloud_type_pixel
   
   !
   !
   !
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
   
   !
   !
   !
   subroutine determine_type_water ( z_opa, t_opa , ctype )
      real, intent (in ) :: z_opa
      real, intent ( in) :: t_opa
      integer, intent(out) :: ctype
   
      ctype = ET_cloud_type % WATER
         if ( t_opa < 273.0) then
            ctype = ET_cloud_type % SUPERCOOLED
         else          
            if (z_opa < 1.0 ) ctype = ET_cloud_type % FOG
         end if
   
   end subroutine determine_type_water
   
   
   !
   !
   subroutine get_ice_probabibility (inp , ice_prob , is_cirrus , is_water , t_opa , z_opa ) 
      implicit none
      type ( cloud_type_input_type) , intent(in) :: inp
      real , intent(out) :: ice_prob
      logical, intent(out) :: is_cirrus
      logical, intent(out) :: is_water
      real, intent(out)    :: t_opa
      real, intent(out)    :: z_opa
      
      integer :: phase
      integer :: ctype
      
      real, parameter:: ICE_TEMP_MIN = 243.0
      real, parameter:: ICE_TEMP_MAX = 263.0
      real, parameter:: FMFT_COLD_OFFSET = 0.5 !K
      real, parameter:: FMFT_CIRRUS_THRESH = 1.0 !K
      real, parameter :: BT_11UM_STD_CIRRUS_THRESH = 4.0 ! K
 
      
      
      real :: fmft
      real :: h2o_correct
     
     
      is_water = .false.
      is_cirrus = .false.
      
      ! - ice probability
      call height_opaque ( &
            inp % sat % rad_ch31 &
         , inp % rtm % rad_ch31_bb_prof &
         , inp % rtm % tropo_lev &
         , inp % rtm % sfc_lev &
         , inp % rtm % t_prof &
         , inp % rtm % z_prof &
         , t_opa &
         , z_opa ) 
         
      !AW TODO- if h2o channel is available you can do this more accuate
      
      ice_prob = 1.0 - (T_opa-ICE_TEMP_MIN)/(ICE_TEMP_MAX - ICE_TEMP_MIN) 
      Ice_Prob = min(1.0,Ice_Prob)
      Ice_Prob = max(0.0,Ice_Prob) 
      
      
      ! -check if water cloud
      
      
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
         
         if ( is_water) ice_prob = 0.0
      end if
      
      
      ! + test for cirrus +++++++
      
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
         ! -----------
         
         ! - 6.7 covariance
         if ( inp % sat % chan_on ( 27 ) &
            .and. inp % sat % chan_on ( 31 ) ) then
            
            if ( inp%rtm%Covar_Ch27_Ch31_5x5 > 1.5 &
               .and. inp % rtm % bt_ch31_3x3_max < 250.0) is_cirrus = .true.         
         end if 
         ! --------------
         
         if ( is_cirrus ) ice_prob = 1.0
      
      end if
      
      ! ++++++++++++++
      
     
      
      
   end subroutine get_ice_probabibility
   
   ! --------------------------------------------
   !
   ! --------------------------------------------
   subroutine height_opaque ( &
        rad31 &
      , rad31_rtm_prof &
      , tropo_lev &
      , sfc_lev &
      , t_prof &
      , z_prof &
      , t_opa &
      , z_opa )
      
      real , intent(in) :: rad31
      real , intent(in) :: rad31_rtm_prof(:)
      integer , intent(in) :: tropo_lev
      integer , intent(in) :: sfc_lev
      real, intent(in) :: t_prof(:)
      real, intent(in) :: z_prof(:)
      real, intent(out) :: t_opa
      real, intent(out) :: z_opa
      
      t_opa = -999.
      z_opa = -999.
      
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


end module cloud_phase_algo_module
    
