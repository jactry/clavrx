! $Header$
subroutine nlcomp_start ( obs_vec &
                      & , obs_u &
                      & , alb_sfc &
                      & , alb_sfc_u &
                      & , state_apr &
                      & , air_trans_ac   &
                      & , sol_zen &
                      & , sat_zen &
                      & , rel_azi &
                      & , cld_temp &
                      & , water_phase &
							 & , rad_abv_cld &
                      & , rad_sfc   &
							 & , cod &
                      & , cps &
                      & , codu &
                      & , cpsu &
							 & , ancil_path &
                      & , debug_in ) 
         
                      
   use nlcomp_mod
   real, dimension ( 2 ) :: obs_vec, obs_u , alb_sfc , alb_sfc_u  
   real, dimension ( 2 ) :: state_apr 
   real, dimension ( 2 ) :: air_trans_ac 
   real, intent ( in ) :: sol_zen, sat_zen , rel_azi , cld_temp

   character ( len = 1024 ) , intent ( in ) :: ancil_path
   integer, intent(in), optional :: debug_in
  
   logical :: water_phase
   real, intent ( in ) ::  rad_abv_cld , rad_sfc 	  
   real, intent ( out ) :: cod , cps ,codu ,cpsu
 		
   call nlcomp_algorithm ( obs_vec , obs_u , alb_sfc , alb_sfc_u , state_apr , air_trans_ac &
                              & , sol_zen, sat_zen , rel_azi , cld_temp , water_phase &
							         & , rad_abv_cld ,  rad_sfc  &
							         & , cod , cps , codu , cpsu, debug_in = debug_in  , ancil_path =  ancil_path )  
							  
						  


end subroutine nlcomp_start
