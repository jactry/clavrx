! $Id$
!
!  HISTORY: 06/05/2014: changed filename for better naming convebtion
!
!           02/05/2014 : add AHI (AW)
subroutine dcomp_array_loop ( input, output , debug_mode_user)

   use dcomp_retrieval_mod
   
   use dcomp_interface_types_mod, only: &
      dcomp_in_type &
      , dcomp_out_type &
      , EM_cloud_type &
      , EM_cloud_mask &
      , EM_snow_class
   
   implicit none
   type (dcomp_in_type) , intent(in) :: input
   type (dcomp_out_type) :: output
   integer , intent(in) , optional :: debug_mode_user
   
   integer, parameter :: real4 = selected_real_kind(6,37)
   integer, parameter :: int4 = selected_int_kind(8)
   integer, parameter :: int1 = selected_int_kind(1)
   integer, parameter :: int2 = selected_int_kind(2)  
   integer :: nr_lines, nr_elem
   
  
    !- scalar local variables

   ! - number of possible channels in CLAVR-x
   integer , parameter :: N_CHAN = 45
  
   real :: refl_toa = -999.
   
   real ::calib_err ( N_CHAN ) 

   real :: obs_vec(2)   = [-999.,-999.]
   real :: obs_unc(2)   = [-999.,-999.]
   real :: alb_vec(2)   = [-999.,-999.]
   real :: alb_unc(2)   = [-999.,-999.]
	real :: trans_vec(2) = [-999.,-999.]
 
   real :: gas_coeff (3)
   
   real :: trans_ozone ( N_CHAN )
   real :: trans_unc_ozone ( N_CHAN )
   real :: trans_rayleigh ( N_CHAN )
   real :: trans_wvp ( N_CHAN )
   real :: trans_unc_wvp ( N_CHAN )
   real :: trans_total ( N_CHAN )
      
   integer , parameter :: CHN_VIS = 1
   integer  :: CHN_NIR
   integer  :: chn_idx

   real( kind = real4 ) :: assumed_tpw_error 
   real( kind = real4 ) :: ozone_coeff (3)  

   ! -- nwp variables 
   real( kind = real4 ) :: ozone_path_nwp
           
   real :: rad_clear_sky_toa_ch20 = -999.
   real :: rad_clear_sky_toc_ch20 = -999.
  
   real     , allocatable :: air_mass_array (:,:)
   logical  , allocatable :: cloud_array(:,:)
   logical  , allocatable :: obs_array(:,:)
   logical  , allocatable :: water_phase_array(:,:)
   
   integer ( kind = int2)  , allocatable :: info_flag ( :,:)
   integer ( kind = int2)  , allocatable :: quality_flag ( :,:)
   integer :: dim_1 
   integer :: dim_2
      
   integer  :: array_dim (2)
   real     :: state_apriori (2)
   
   
   integer :: debug_mode
      
   character (len = 20) :: sensor_name
      
   real :: cld_height
   real :: cld_press
   real :: cld_temp
   
   real :: rel_azi
   
   real ( kind = real4 ) :: refl_toc(N_CHAN)
   real ( kind = real4 ) :: alb_sfc(N_CHAN)
   real ( kind = real4 ) :: alb_unc_sfc(N_CHAN)
   real ( kind = real4 ) :: rad_to_refl_factor
   
   real, parameter :: SAT_ZEN_MAX = 75.
   real, parameter :: SOL_ZEN_MAX = 75.
   real, parameter :: PI = 4. * ATAN(1.)
   
   real ( kind = real4) :: ALBEDO_OCEAN (N_CHAN)
      
   type ( dcomp_output_structure ) :: dcomp_out
      
   real :: sol_zen 
   real :: sat_zen 
   
   integer :: line_idx 
   integer ::  elem_idx
   
   
   
   integer :: tried 
   integer :: success
   
   interface
      subroutine view2d(x,a,b,t)
         real,dimension(:,:) :: x
         real, intent(in), optional::a,b
         character ( len =*) , optional :: t
      end subroutine view2d
   end interface
   
   ! - executable ---------
   
   debug_mode = 4
   if ( present ( debug_mode_user)) debug_mode = debug_mode_user
  
   array_dim = shape ( input % sat % d )
   dim_1 = array_dim (1) 
   dim_2 = array_dim (2)
   nr_lines = array_dim(2)
   nr_elem = array_dim(1)
      
   ALBEDO_OCEAN (:) = 0.03
   calib_err ( : ) = 0.03
   
   
   allocate ( obs_array ( dim_1 , dim_2 ) &
                  ,  cloud_array ( dim_1 , dim_2 ) )
   allocate ( water_phase_array (  dim_1 , dim_2 ) )	
   allocate ( air_mass_array  ( dim_1 , dim_2 ) )		   
   
      
   air_mass_array = 1.0 / cos (input % sat % d * PI / 180. ) + 1.0 / cos ( input % sol % d * pi / 180.)
      
   obs_array = input % is_valid % d   &
                       & .and. input % sat % d <= sat_zen_max &
                       & .and. input % sol % d <= sol_zen_max &
                       & .and. input % refl (1) % d  >= 0. &
                       & .and. air_mass_array >= 2.
	 
   cloud_array =  obs_array &
                        & .and. ( input % cloud_mask % d == EM_cloud_mask % CLOUDY &
                        & .or. input % cloud_mask % d == EM_cloud_mask % PROB_CLOUDY ) &
                        & .and. input % cloud_temp % d > 10 
     
   water_phase_array = input % cloud_type % d == EM_cloud_type % FOG &
                        &  .or. input % cloud_type % d == EM_cloud_type % WATER &
                        &  .or. input % cloud_type % d == EM_cloud_type % SUPERCOOLED &
                        &  .or. input % cloud_type % d == EM_cloud_type % MIXED 
  

   !--- select sensor name based on wmo id number
   select case (input % sensor_wmo_id)
      case(3)
         sensor_name = 'METOP-B'
      case(4)
         sensor_name = 'METOP-A'	 
      case(55)
         sensor_name = 'Meteosat-8' 
      case(56)
         sensor_name = 'Meteosat-9' 
      case(57)
         sensor_name = 'Meteosat-10' 
      case(70)
         sensor_name = 'Meteosat-11' 
      case(171)
         sensor_name = 'MTSAT-1R'
      case (172)
         sensor_name = 'MTSAT-2'
      case (173)
         sensor_name = 'AHI'     
      case (200:204)
         write(sensor_name, "('NOAA-',i2.2)") input % sensor_wmo_id - 192 
      case (205:209)
         write(sensor_name, "('NOAA-',i2.2)") input % sensor_wmo_id - 191
      case(223)
         sensor_name = 'NOAA-19'	  
      case(224)
         sensor_name = 'VIIRS'
      case(252)
         sensor_name = 'GOES-08'
      case(253)
         sensor_name = 'GOES-09'
      case(254)
         sensor_name = 'GOES-10'
      case(255)
         sensor_name = 'GOES-11'
      case(256)
         sensor_name = 'GOES-12'
      case(257)
         sensor_name = 'GOES-13'
      case(258)
         sensor_name = 'GOES-14'
      case(259)
         sensor_name = 'GOES-15'
      case(705)
         sensor_name = 'NOAA-05'   
      case(706)
         sensor_name = 'NOAA-06'   
      case(707)
         sensor_name = 'NOAA-07'
      case(708)
         sensor_name = 'TIROS-N'   
      case(783)
         sensor_name = 'MODIS-TERRA'
      case(784)
         sensor_name = 'MODIS-AQUA'
      case(810)
         sensor_name = 'COMS-1'   
      case default
         print*,'please inform  andi.walther@ssec.wisc.edu wmo id: ', input % sensor_wmo_id 
         stop
   end select
  
   ozone_coeff  = [ -0.000606266 , 9.77984e-05,-1.67962e-08 ] 
   
   output = dcomp_out_type ( dim_1, dim_2 )

   where ( obs_array .and. .not. cloud_array )
      output % cld_trn_sol % d   =  1. 
      output % cld_trn_obs % d   =  1.  
      output % cld_alb % d       =  0.
      output % cld_sph_alb % d   =  0.   
   end where
    
   allocate ( info_flag ( dim_1, dim_2))
   allocate ( quality_flag ( dim_1, dim_2))
   
   ! - initialize
   quality_flag  = ibclr ( quality_flag , 0)
   quality_flag  = ibset ( quality_flag , 1)
   quality_flag  = ibset ( quality_flag , 2)
   quality_flag  = ibset ( quality_flag , 3)
   quality_flag  = ibset ( quality_flag , 4)
   quality_flag  = ibset ( quality_flag , 5)
   quality_flag  = ibclr ( quality_flag , 6)
   quality_flag  = ibclr ( quality_flag , 7)
   
   ! - initialize
   info_flag  = 0
      
   ! - check input options
   
   select case ( input % mode )
      case ( 1 ) 
         CHN_NIR = 6               
      case ( 2 )
         CHN_NIR = 7
      case ( 3)
         CHN_NIR = 20
      case default
        
   end select
    
   if ( input % is_channel_on (CHN_NIR) .eqv. .false.) then
      print*, 'dcomp NIR channel not set! ==> MODIS equaivalant channel: ', CHN_NIR
      print*, 'all dcomp results are set to missing values'
      return
   end if
   
   if ( input % is_channel_on (CHN_VIS) .eqv. .false.) then
      print*, 'dcomp VIS channel not set! ==> MODIS equaivalant channel: ', CHN_VIS
      print*, 'all dcomp results are set to missing values'
      return
   end if
   
   
   
   if ( debug_mode == 6 ) then
      call view2d ( input % cloud_press % d ,0., 1200., 'Cloud_press')
   end if
   
   line_loop: do line_idx = 1 , nr_lines
      elem_loop: do elem_idx = 1,   nr_elem
         
   
         if ( .not. cloud_array (elem_idx,line_idx)  ) cycle elem_loop
         
         ! - set aliases
         cld_height =  input % cloud_hgt % d (elem_idx,line_idx)
         cld_press  =  input % cloud_press % d (elem_idx,line_idx)
         cld_temp   =  input % cloud_temp % d (elem_idx,line_idx)
         sol_zen    =  input % sol % d (elem_idx,line_idx)
         sat_zen    =  input % sat % d (elem_idx,line_idx)
         rel_azi    =  input % azi % d (elem_idx,line_idx)
         ozone_path_nwp = input % ozone_nwp % d (elem_idx,line_idx)
                  
         ! - compute transmission 
              
         loop_chn: do chn_idx = 1 , 40
       
            if ( input % is_channel_on (chn_idx) .eqv. .false.) cycle  loop_chn
            
            trans_block : associate ( tpw_ac => input % tpw_ac % d (elem_idx,line_idx)  , &
                     
                     press_sfc => input % press_sfc  % d (elem_idx,line_idx) , &
                     trans_chn20_ac_nadir => input % trans_ac_nadir(20) % d , & 
                     refl_toa => input % refl (chn_idx)  % d (elem_idx, line_idx))
                     
               gas_coeff = input % gas_coeff ( chn_idx) % d  
               trans_ozone( chn_idx ) = 1.
               trans_rayleigh( chn_idx ) = 1.
            
               if ( chn_idx == CHN_VIS ) then
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
                                
               trans_unc_wvp ( chn_idx ) = abs(trans_wvp( chn_idx ) - exp ( -1. * (gas_coeff(1)   &
                           & + gas_coeff(2) * (assumed_tpw_error * tpw_ac) &
                           & + gas_coeff(3) * ( ( assumed_tpw_error * tpw_ac ) **2 ) ) ) )       
                        
            end associate trans_block
           
            trans_total ( chn_idx ) = trans_rayleigh( chn_idx ) * trans_ozone( chn_idx ) * trans_wvp( chn_idx )
            trans_total ( chn_idx ) = trans_ozone( chn_idx ) * trans_wvp( chn_idx )
            
            
            refl_toc( chn_idx ) = refl_toa  * trans_total (chn_idx )
            
            alb_sfc( chn_idx ) =  ( input % alb_sfc ( chn_idx )  % d (elem_idx,line_idx) ) / 100.
            
            alb_sfc( chn_idx ) = max ( alb_sfc( chn_idx ) , ALBEDO_OCEAN (chn_idx) )
            
            alb_unc_sfc  (chn_idx) = 0.05
                        
            if ( chn_idx == 20 ) then
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
       
         if ( input % mode == 3) then
            obs_vec( 2 ) = refl_toc ( 20 )
            obs_unc( 2 ) = obs_vec( 2 ) * 0.1 
         else
            obs_vec( 2 ) = input % refl ( CHN_NIR)  % d (elem_idx, line_idx) /100.
            obs_unc( 2 ) = max ( trans_unc_wvp  ( CHN_NIR ) , 0.01 )  + calib_err (CHN_NIR)
         end if
         
         
         alb_vec( 2 ) = alb_sfc ( CHN_NIR)
         alb_unc( 2 ) = 0.05
         
         trans_vec ( 2) = trans_total ( CHN_NIR )
              
         ! - apriori
         state_apriori (1) = 0.7 * ( 100. * obs_vec(1) ) ** (0.9)
         state_apriori(1) = log10 ( max ( 0.1 , state_apriori(1) ) ) 
         state_apriori (2) = 1.3
         if  (water_phase_array ( elem_idx, line_idx) ) state_apriori(2) = 1.0
         
                  
         call dcomp_algorithm ( &
                &   obs_vec &
                & , obs_unc  &
                & , alb_vec &
                & , alb_unc  &
                & , state_apriori &
                & , trans_vec  &
                & , sol_zen &
                & , sat_zen &
                & , rel_azi &
                & , cld_temp &
                & , water_phase_array ( elem_idx, line_idx)  &
                & , rad_clear_sky_toc_ch20 &
                & , rad_clear_sky_toa_ch20 &
                & , trim(sensor_name) &
                & , dcomp_out &
                & , input % mode &
                & , input % lut_path &
                & , debug_mode  )
        
        if ( debug_mode == 4  ) then
            print*,'=======================> input:',CHN_NIR
            print*,'Elem Line: ', elem_idx,line_idx
            print*,' Obs vector: ',obs_vec
            print*,' Obs uncert: ',obs_unc
            print*,' Surface Albedo: ', alb_vec
            print*,' Surface albedo unc: ',alb_unc
            print*,' Transmission: ',trans_vec
            print*, 'Angles: ',sol_zen,sat_zen,rel_azi
            print*, 'Cloud temp mask: ',cld_temp,water_phase_array( elem_idx, line_idx)
            print*, 'Ch20 rtm: ', rad_clear_sky_toc_ch20 , rad_clear_sky_toa_ch20
            print*, 'output: '
            print*, dcomp_out % cod, dcomp_out % cps
            print*
           
            print*,'==============================='
           
         end if
         
         output % cod % d (elem_idx,line_idx) = dcomp_out % cod
         output % cps % d (elem_idx, line_idx) = dcomp_out % cps 
     
         output % cod_unc % d ( elem_idx, line_idx) = dcomp_out % codu
         output % ref_unc % d ( elem_idx, line_idx) = dcomp_out % cpsu
         output % cld_trn_sol % d ( elem_idx, line_idx) = dcomp_out % cloud_trans_sol_vis  
         output % cld_trn_obs % d ( elem_idx, line_idx) = dcomp_out % cloud_trans_sat_vis   
         output % cld_alb % d ( elem_idx, line_idx)  = dcomp_out % cloud_alb_vis 
         output % cld_sph_alb % d ( elem_idx, line_idx)  = dcomp_out % cloud_sph_alb_vis 
         
         ! - DCOMP_QF_PROCESSION B0
         quality_flag (elem_idx,line_idx) = ibset ( quality_flag(elem_idx,line_idx) , 0)
         
         if ( dcomp_out % statusOK ) then
            ! - DCOMP_QF_COD_VALID B1
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 1)
            ! - DCOMP_QF_REF_VALID B2
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 2)
            ! - initial DCOMP_QF_COD_DEGRADED B3
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 3)
            ! - initial DCOMP_QF_REF_DEGRADED B4
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 4)
            ! --DCOMP_QF_REF_CONVERGENCY B5
            quality_flag (elem_idx,line_idx) = ibclr ( quality_flag(elem_idx,line_idx) , 5)
            
         end if                
      end do elem_loop
   end do   line_loop 

   ! DCOMP_INFO_LAND_SEA I0
   where (input % is_land % d )
      info_flag = ibset (  info_flag , 0 ) 
   end where 

  ! -DCOMP_INFO_DAY_NIGHT I1
   where (input % sol % d  > 82. )
	   info_flag = ibset ( info_flag , 1)
   end where
    
  ! - DCOMP_INFO_TWILIGHT I2
   where (input % sol % d > 65. .and. input % sol % d < 82.)
      info_flag = ibset ( info_flag , 2) 
   end where 
    
   ! - DCOMP_INFO_SNOW I3
   where ( input % snow_class % d == EM_snow_class % SNOW )
      info_flag = ibset ( info_flag , 3) 
   end where 
   
   ! - DCOMP_INFO_SEA_ICE I4
   where ( input % snow_class % d == EM_snow_class % SEA_ICE )
	   info_flag = ibset ( info_flag , 4) 
   end where
   
   ! -DCOMP_INFO_PHASE
   where ( .not. water_phase_array)
      info_flag = ibset ( info_flag, 5)
   end where
   
   ! -DCOMP_INFO_THICK_CLOUD
   where ( output % cod % d > 80 )
      info_flag = ibset ( info_flag, 6)
   end where
   
   ! -DCOMP_INFO_THIN_CLOUD
   where ( output % cod % d < 4 .and. output % cod % d > 0 )
      info_flag = ibset ( info_flag, 7)
   end where
 
   ! - DCOMP_QF_COD_DEGRADED B3
   where ( (  btest(info_flag,2) &
            .or. btest(info_flag,3) &
            .or. btest(info_flag,4) &
            .or. btest(info_flag,5)  &
            .or. btest(info_flag,6)) &
	        .and. btest(quality_flag,0)) 
         quality_flag = ibset ( quality_flag , 3 )
   end where
   
   ! - DCOMP_QF_REF_DEGRADED B4
   where (    ( btest(info_flag,2) &
             .or. btest(info_flag,3) &
             .or. btest(info_flag,4) &
             .or. btest(info_flag,5)  &
             .or. btest(info_flag,7)) &
	       &  .and. btest(quality_flag,0)) 
	        quality_flag = ibset ( quality_flag , 4 )
   end where
   
   ! - compute lwp
   where ( water_phase_array &              
	            .and. .not. btest ( quality_flag , 1 )  &
		         .and. .not. btest ( quality_flag , 2 ) )
          output % lwp % d =  output % cod % d * output % cps % d * 5.0 / 9.0
   end where
   
      ! - compute lwp
   where ( .not. water_phase_array &              
	            .and. .not. btest ( quality_flag , 1 )  &
		         .and. .not. btest ( quality_flag , 2 ) )
          output % iwp % d =  (output % cod % d  ** (1/0.84) ) / 0.065 
   end where
   
      
   output % quality % d = quality_flag
   output % info % d = info_flag

   ! compute successrate
   output % nr_obs = count ( obs_array)
   output % nr_clouds = count ( cloud_array)
   output % nr_success_cod = count (btest ( quality_flag , 1 ))
   output % nr_success_cps = count (btest ( quality_flag , 2 ))
   
   tried =  count (btest(quality_flag,0))
    
	output % successrate = 0.0	
      if ( tried > 0 ) then
         success = count(.not. btest(quality_flag,1) .and. .not. btest ( quality_flag,2) ) 
         output % successrate = success / tried
      end if
     
   deallocate ( obs_array &
	               ,  cloud_array  )
   deallocate ( water_phase_array )	
	  
   deallocate ( air_mass_array ) 	
   
   output % version = '$Id$'	  
   

end subroutine dcomp_array_loop



