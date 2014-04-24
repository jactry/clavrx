! $Header:$
module cloud_phase_algo_module
   
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
   type ( et_cloudphase_class_type ) :: ET_cloud_phase
   
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
   
   
   
   type cloud_phase_sat_type
      logical , dimension(42) :: chan_on
      real :: rad_ch31 
      real :: ref_ch6 
      real :: ref_ch20
      real :: bt_ch31 
      real :: bt_ch32   
   end type cloud_phase_sat_type
   
   type cloud_phase_rtm_type
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
   end type cloud_phase_rtm_type
   
   type cloud_phase_geo_type   
      real :: sol_zen
   end type cloud_phase_geo_type
   
   type cloud_phase_sfc_type   
      real :: emiss_ch20
   end type cloud_phase_sfc_type   
   
   type cloud_phase_input_type
      type ( cloud_phase_sat_type) :: sat
      type ( cloud_phase_rtm_type) :: rtm
      type ( cloud_phase_geo_type) :: geo
      type ( cloud_phase_sfc_type) :: sfc
   end type cloud_phase_input_type

contains
   subroutine cloud_phase_algo (inp , phase , ctype) 
      implicit none
      type ( cloud_phase_input_type) , intent(in) :: inp
      integer, intent(out) :: phase
      integer, intent(out) :: ctype
      real :: ice_prob
      real, parameter:: ICE_TEMP_MIN = 243.0
      real, parameter:: ICE_TEMP_MAX = 263.0
      real, parameter:: FMFT_COLD_OFFSET = 0.5 !K
      real, parameter:: FMFT_CIRRUS_THRESH = 1.0 !K
      real, parameter :: BT_11UM_STD_CIRRUS_THRESH = 4.0 ! K
      
      real :: t_opa
      real :: z_opa
      real :: fmft
      real :: h2o_correct
      logical :: is_water_phase
      logical :: is_cirrus
      
      phase = ET_cloud_phase % UNKNOWN
      ctype = ET_cloud_type % UNKNOWN
      
      is_water_phase = .false.
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
            .and. inp % sat % ref_ch6 > 30. )  is_water_phase = .true. 
         
         ! - 3.75 day   
         if ( ice_prob /= 0.0 &
            .and. inp % sat % chan_on ( 20 ) &
            .and. inp % geo % sol_zen < 80. &
            .and. inp % sfc % emiss_ch20 > 0.9 &
            .and. inp % sat % ref_ch20 > 20.0 ) is_water_phase = .true.
         
         
         ! -3.75 night 
         if ( ice_prob /= 0.0 &
            .and. inp % sat % chan_on ( 20 ) &
            .and. inp % geo % sol_zen > 80. &          
            .and. inp % sat % ref_ch20 > 5.0 ) is_water_phase = .true. 
         
         if ( is_water_phase) ice_prob = 0.0
      end if
      
      
      ! - test for cirrus
      
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
         
         ! - 6.7 covariance
         if ( inp % sat % chan_on ( 27 ) &
            .and. inp % sat % chan_on ( 31 ) ) then
            
            if ( inp%rtm%Covar_Ch27_Ch31_5x5 > 1.5 &
               .and. inp % rtm % bt_ch31_3x3_max < 250.0) is_cirrus = .true.         
         end if 
         if ( is_cirrus ) ice_prob = 1.0
      
      end if
      
      if ( ice_prob > 0.5 ) then
         phase = ET_cloud_phase % ICE
         ctype = ET_cloud_type % OPAQUE_ICE
         
         if ( inp % Rtm % Emiss_tropo_ch31 > 0.95) then
            ctype = ET_cloud_type % OVERSHOOTING
         end if
         
         if ( inp % Rtm % Emiss_tropo_ch31 < 0.8) then
            ctype = ET_cloud_type % CIRRUS
            !TODO - beta stuff
            if ( is_cirrus .and. is_water_phase) ctype = ET_cloud_type % OVERLAP
            
         end if
         
         if ( inp % sat % bt_ch31 < 233.0 ) then
            ctype = ET_cloud_type % OPAQUE_ICE
         end if 
         
      else
         phase = ET_cloud_phase % WATER
         if ( t_opa < 273.0) then
            phase = ET_cloud_phase % SUPERCOOLED
            ctype = ET_cloud_type % SUPERCOOLED
         end if
         if (z_opa < 1.0 ) ctype = ET_cloud_type % FOG
      end if      
      
      
   end subroutine cloud_phase_algo
   
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
