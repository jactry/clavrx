! $Id: test_forward_computation.f90 5 2014-01-21 23:42:33Z awalther $

program test_forward_computation

   use dcomp_forward_mod
   implicit none
   real, dimension ( 2 ) :: obs_vec , alb_sfc , alb_sfc_u  , state_apr 
   real :: sol_zen, sat_zen , rel_azi ,cld_temp
   real ::  rad_abv_cld , trns_abv_cld , rad_sfc , rad_to_refl
    
   real, dimension ( 2 ) :: obs_fwd 
   real, dimension(2) :: trns_sol, trns_sat, sph_alb
   
   real, dimension(2,2) :: kernel
   real, dimension (2) :: state_vec
   real , dimension (2) :: air_trans_ac
   
   real  :: solar_rad_20
   real :: sun_earth_sistance 
   
   integer :: start_day
   
   real :: sun_earth_distance
   real :: solar_goes15
   real :: ew_20
   real :: kk
   integer :: jj
   logical :: wat_phase
   integer :: dcomp_mode
   dcomp_mode = 3
   wat_phase = 'true'
   
   obs_vec = [ 0.2211 , 0.433 ]
   
   alb_sfc = [ 0.15 , 0.15 ]
   alb_sfc_u = [ 0.01 , 0.01 ]
   state_vec =  [0.22 , 1. ] 
   air_trans_ac = [ 1.0 , 0.29 ]
   sol_zen  = 23.33
   sat_zen = 21.33
   rel_azi = 120.
   cld_temp = 276.
   
   rad_abv_cld    = 0.001
   trns_abv_cld   = 0.86
   rad_sfc        = 0.1
   
   start_day = 30
   sun_earth_distance = 1.0 - 0.016729 * cos ( 0.9856 * ( start_day - 4.0) * DTOR )
   solar_goes15 = 3.5111240
   ew_20 = 241.02064
   solar_rad_20 =  1000.* solar_goes15 / ew_20
   rad_to_refl    = PI / cos (sol_zen * DTOR )/ solar_rad_20 / sun_earth_distance ** 2
   
   
   kk = 254
 
   
   call  dcomp_forward_computation  ( &
                             state_vec  &
						   , sol_zen , sat_zen , rel_azi , wat_phase , 'MODIS-AQUA' , dcomp_mode &
						   , alb_sfc , air_trans_ac  &
						   , obs_fwd  &
						   , trns_sol &
						   , trns_sat &
						   , sph_alb &
						   , kernel , kk  &
						   , rad_abv_cld , rad_sfc  )  
        
   
   print*,'forward test=> ',obs_fwd
   
   
end program test_forward_computation
