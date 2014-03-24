! $Header$
!
!  HISTORY: 2014/01/12
!         : AW first verisob of NLCOMP for arrays
! 
!  TODO: Angles are wrong ( solar vs lunar
!
!
subroutine nlcomp_array ( input , output, debug_mode_user )
   use nlcomp_mod
   use nlcomp_interface_def_mod
   implicit none
   
   type (nlcomp_in_type) , intent(in) :: input
   type (nlcomp_out_type), intent(out) :: output
   integer , intent(in) , optional :: debug_mode_user
   
   ! parameters
   integer , parameter :: N_CHAN = 42
   real  :: ALBEDO_OCEAN (N_CHAN)
   real, parameter :: SAT_ZEN_MAX = 80.
   real, parameter :: SOL_ZEN_MIN = 90.
   
   real, parameter :: PI = 3.14159265359
   real, parameter :: MISSING_REAL4 = -999.
   
   ! - local logical arrays 
   real,  allocatable    :: air_mass_array( : , : )
   logical , allocatable :: is_cloud( : , : )
   logical , allocatable :: is_obs( : , : )
   logical , allocatable :: is_water_phase( : , : )
   logical , allocatable :: has_city_lights(:,:)
   
   
   ! - 
   type ( nlcomp_output_structure ) :: nlcomp_out
   
   integer :: debug_mode
   
   integer :: array_dim(2)
   integer :: dim_1 
   integer :: dim_2
   integer :: nr_lines
   integer :: nr_elem
   
   
   real ::calib_err ( N_CHAN )
   
   
   
    real :: cld_height
   real :: cld_press
   real :: cld_temp
   
   real :: rel_azi , lunar_rel_azi
   real :: sol_zen , sat_zen , lunar_zen
   integer :: line_idx , elem_idx
   
   integer :: chn_idx
   integer , parameter :: CHN_VIS = 42
   integer ,parameter  :: CHN_NIR = 20
   real(kind=real4) :: gas_coeff (3)
   
   real( kind = real4 ) :: trans_ozone ( N_CHAN )
   real( kind = real4 ) :: trans_unc_ozone ( N_CHAN )
   real( kind = real4 ) :: trans_rayleigh ( N_CHAN )
   real( kind = real4 ) :: trans_wvp ( N_CHAN )
   real( kind = real4 ) :: trans_unc_wvp ( N_CHAN )
   real( kind = real4 ) :: trans_total ( N_CHAN )
   
   real( kind = real4 ) :: assumed_tpw_error 
   real( kind = real4 ), parameter :: ozone_coeff (3)    = [ -0.000606266 , 9.77984e-05,-1.67962e-08 ]
   
   real ( kind = real4 ) :: refl_toc(40)
   real ( kind = real4 ) :: alb_sfc(40)
   real ( kind = real4 ) :: alb_unc_sfc(40)
   real ( kind = real4 ) :: rad_to_refl_factor
   real(kind=real4)   :: refl_toa = -999.
   
   real :: rad_clear_sky_toa_ch20 = -999.
   real :: rad_clear_sky_toc_ch20 = -999.
   
   real :: obs_vec(2) = [-999.,-999.]
   real :: obs_unc(2) = [-999.,-999.]
   real :: alb_vec(2) = [-999.,-999.]
   real :: alb_unc(2) = [-999.,-999.]
	real :: trans_vec(2) = [-999.,-999.]
   
  
   real , dimension(2) :: state_apriori
   
    ! -- nwp variables 
   real :: ozone_path_nwp
   
   ! - executable
   
   
   debug_mode = 1
   if ( present ( debug_mode_user)) debug_mode = debug_mode_user
   
   array_dim = shape ( input % sat % d )
   dim_1 = array_dim (1) 
   dim_2 = array_dim (2)
   nr_lines = array_dim(2)
   nr_elem = array_dim(1)
      
   ALBEDO_OCEAN (:) = 0.03
   calib_err ( : ) = 0.03
  
   allocate ( is_obs ( dim_1 , dim_2 ) &
                  ,  is_cloud ( dim_1 , dim_2 ) )
   allocate ( is_water_phase (  dim_1 , dim_2 ) )	
   allocate ( air_mass_array  ( dim_1 , dim_2 ) )
   allocate ( has_city_lights  ( dim_1 , dim_2 ) )
   
    
   ! - flag masks
   air_mass_array = 1.0 / cos (input % sat % d * pi / 180. ) &
                & + 1.0 / cos ( input % zen_lunar % d * pi / 180.)
   
  
   is_obs = input % is_valid % d   &
                       & .and. input % sat % d <= SAT_ZEN_MAX &
                       & .and. input % sol % d > SOL_ZEN_MIN &
                       & .and. input % zen_lunar % d < 80 &
                       & .and. input % refl (42) % d  >= 0. &
                       & .and. air_mass_array >= 2.
	
   has_city_lights = input % rad (42) % d > 1.E-06
   
    
   is_cloud =  is_obs &
                        & .and. ( input % cloud_mask % d == EM_cloud_mask % CLOUDY &
                        & .or. input % cloud_mask % d == EM_cloud_mask % PROB_CLOUDY ) &
                        & .and. input % cloud_temp % d > 10 
    
   is_water_phase = input % cloud_type % d == EM_cloud_type % FOG &
                        &  .or. input % cloud_type % d == EM_cloud_type % WATER &
                        &  .or. input % cloud_type % d == EM_cloud_type % SUPERCOOLED &
                        &  .or. input % cloud_type % d == EM_cloud_type % MIXED 
  
  !-allocation
   allocate ( output % cod % d         ( dim_1 , dim_2))
   allocate ( output % cps % d         ( dim_1 , dim_2))
   allocate ( output % cod_unc % d     ( dim_1 , dim_2))
   allocate ( output % ref_unc % d     ( dim_1 , dim_2))
   allocate ( output % cld_trn_sol % d ( dim_1 , dim_2))
   allocate ( output % cld_trn_obs % d ( dim_1 , dim_2))
   allocate ( output % cld_alb % d     ( dim_1 , dim_2))
   allocate ( output % cld_sph_alb % d ( dim_1 , dim_2))
   allocate ( output % info % d        ( dim_1 , dim_2))
   allocate ( output % quality % d     ( dim_1 , dim_2))
   allocate ( output % lwp % d         ( dim_1 , dim_2))
   allocate ( output % iwp % d         ( dim_1 , dim_2))
  
   
   output % cod % d           =  MISSING_REAL4 
   output % cps % d           =  MISSING_REAL4 
   output % cod_unc % d       =  MISSING_REAL4
   output % ref_unc % d       =  MISSING_REAL4 
   output % cld_trn_sol % d   =  MISSING_REAL4  
   output % cld_trn_obs % d   =  MISSING_REAL4   
   output % cld_alb % d       =  MISSING_REAL4  
   output % cld_sph_alb % d   =  MISSING_REAL4  

   
   line_loop: do line_idx = 1 , nr_lines
      elem_loop: do elem_idx = 1,   nr_elem
         
         
         if ( .not. is_cloud (elem_idx,line_idx)  ) cycle elem_loop
         
         ! - set local aliases
         cld_height     = input % cloud_hgt % d (elem_idx,line_idx)
         cld_press      = input % cloud_press % d (elem_idx,line_idx)
         cld_temp       = input % cloud_temp % d (elem_idx,line_idx)
         sol_zen        = input % sol % d (elem_idx,line_idx)
         lunar_zen      = input % zen_lunar % d (elem_idx,line_idx)
         sat_zen        = input % sat % d (elem_idx,line_idx)
         rel_azi        = input % azi % d (elem_idx,line_idx)
         lunar_rel_azi  = input % azi_lunar % d (elem_idx,line_idx)
         ozone_path_nwp = input % ozone_nwp % d (elem_idx,line_idx)
         
         
         ! - compute transmission 
              
         loop_chn: do chn_idx = 1 , 42
           
            if ( input % is_channel_on (chn_idx) .eqv. .false.) cycle  loop_chn
            
            trans_block : associate ( tpw_ac => input % tpw_ac % d (elem_idx,line_idx)  , &
                      & press_sfc => input % press_sfc  % d (elem_idx,line_idx) , &
                      & trans_chn20_ac_nadir => input % trans_ac_nadir(20) % d , & 
                      & refl_toa => input % refl (chn_idx)  % d (elem_idx, line_idx))
                     
               gas_coeff = input % gas_coeff ( chn_idx) % d  
               trans_ozone( chn_idx ) = 1.
               trans_rayleigh( chn_idx ) = 1.
            
               if ( chn_idx == CHN_VIS ) then
               
                  ! - use default ozone value for bad data
                  if (ozone_path_nwp < 100) ozone_path_nwp = 320
               
                  trans_ozone( chn_idx ) = exp ( -1. * ( ozone_coeff(1) &
                              & + ozone_coeff(2) *  ozone_path_nwp &
                              & + ozone_coeff(3) *  ozone_path_nwp**2))

                  trans_unc_ozone( chn_idx ) = trans_ozone(  chn_idx  ) -  exp ( -1. * ( ozone_coeff(1) &
                              & + ozone_coeff(2) *  (1.1 * ozone_path_nwp) &
                              & + ozone_coeff(3) * (1.1 * ozone_path_nwp)**2))
         
                  trans_ozone( chn_idx ) = min ( trans_ozone(  chn_idx  ) , 1. ) 
         
                  ! - rayleigh
                  trans_rayleigh (  chn_vis  ) = exp (-air_mass_array(elem_idx,line_idx)  &
                    &    * ( 0.044 *  (cld_press / press_sfc )) * 0.84)
            
               end if
              
               assumed_tpw_error = 1.2
                           
               trans_wvp ( chn_idx ) =  exp( - 1. * (gas_coeff(1) &
                          & + gas_coeff(2) * tpw_ac  &
                          & + gas_coeff(3) * ( tpw_ac ** 2 ) ) )
                                
               trans_unc_wvp ( chn_vis ) = abs(trans_wvp( chn_idx ) - exp ( -1. * (gas_coeff(1)   &
                           & + gas_coeff(2) * (assumed_tpw_error * tpw_ac) &
	                        & + gas_coeff(3) * ( ( assumed_tpw_error * tpw_ac ) **2 ) ) ) )       
                        
            end associate trans_block
           
            trans_total ( chn_idx ) = trans_rayleigh( chn_idx ) * trans_ozone( chn_idx ) * trans_wvp( chn_idx )
            trans_total ( chn_idx ) = trans_ozone( chn_idx ) * trans_wvp( chn_idx )
            
            
            refl_toc( chn_idx ) = refl_toa  * trans_total (chn_idx )
            
            alb_sfc( chn_idx ) =  ( input % alb_sfc ( chn_idx )  % d (elem_idx,line_idx) ) / 100.
            
            alb_sfc( chn_idx ) = max ( alb_sfc( chn_idx ) , ALBEDO_OCEAN (chn_idx) )
            
            alb_unc_sfc  (chn_idx) = 0.05
            
           
            
            if ( chn_idx == CHN_NIR ) then
               chn20_block : associate ( rad_toa => input % rad (chn_idx)  % d (elem_idx, line_idx) &
                        & , sun_earth_dist => input % sun_earth_dist &
                        & , solar_irradiance => input % solar_irradiance (chn_idx) )
               
                 
                  trans_total (chn_idx) = input % trans_ac_nadir ( chn_idx) % d  (elem_idx, line_idx)
                  rad_to_refl_factor = PI / cos ( sol_zen * PI / 180.) / ( solar_irradiance / input % sun_earth_dist ** 2 )
                  refl_toc( chn_idx ) = rad_toa * rad_to_refl_factor
                  
               end associate chn20_block
               
               rad_clear_sky_toc_ch20 = input % rad_clear_sky_toc ( 20) % d (elem_idx, line_idx) 
               rad_clear_sky_toa_ch20 = input % rad_clear_sky_toa ( 20) % d (elem_idx, line_idx)
               
            end if
             
         end do loop_chn
         
          ! - NIR
              
              
             
         obs_vec ( 1 ) = input % refl (CHN_VIS)  % d (elem_idx, line_idx) / 100.
         obs_unc ( 1 ) =   trans_unc_ozone ( CHN_VIS) +  trans_unc_wvp  ( CHN_VIS)  +calib_err (CHN_VIS)
         alb_vec ( 1 ) =  alb_sfc ( CHN_VIS)
        
         alb_unc ( 1) = 0.05
        
         trans_vec ( 1) = trans_total ( CHN_VIS )
        
         
         
         ! - nir channel
         obs_vec( 2 ) = input % rad ( CHN_NIR)  % d (elem_idx, line_idx)
         
         
         obs_unc( 2 ) = max ( trans_unc_wvp  ( CHN_NIR ) , 0.01 )  + calib_err (CHN_NIR)
         
         alb_vec( 2 ) = alb_sfc ( CHN_NIR)
         alb_unc( 2 ) = 0.05
         
         trans_vec ( 2) = trans_total ( CHN_NIR )
             
         ! - apriori
         state_apriori (1) = 0.7 * ( 100. * obs_vec(1) ) ** (0.9)
         state_apriori(1) = log10 ( max ( 0.1 , state_apriori(1) ) ) 
         state_apriori (2) = 1.3
         if  (is_water_phase ( elem_idx, line_idx) ) state_apriori(2) = 1.0
         
         call nlcomp_algorithm ( obs_vec , obs_unc  &
                & , alb_vec , alb_unc  , state_apriori , trans_vec  &
                & , lunar_zen , sat_zen &
                & , lunar_rel_azi , cld_temp , is_water_phase ( elem_idx, line_idx)  &
                & , rad_clear_sky_toc_ch20 , rad_clear_sky_toa_ch20  &
                & , nlcomp_out ,  ancil_path =  input % lut_path,  debug_in =  debug_mode  )
                
                
         output % cod % d (elem_idx,line_idx) = nlcomp_out % cod       
         output % cps % d (elem_idx,line_idx) = nlcomp_out % cps
         
         output % cod_unc % d ( elem_idx, line_idx) = nlcomp_out % codu
         output % ref_unc % d ( elem_idx, line_idx) = nlcomp_out % cpsu 
     
      end do elem_loop
   end do   line_loop 
   

end subroutine nlcomp_array
