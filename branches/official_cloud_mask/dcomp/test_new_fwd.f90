! $Id: test_new_fwd.f90 5 2014-01-21 23:42:33Z awalther $
program test_new_fwd

   use dcomp_fwd_mod
   integer , parameter :: n_chn = 4
   type ( pixel_vec ) :: pxl
   real, dimension ( 2 ) :: state_vec
   integer, dimension ( n_chn ) :: channels
   real, dimension ( n_chn ) :: alb_sfc , air_trans_ac
   real, dimension ( n_chn ) :: fwd_vec
   real , dimension ( n_chn )  :: cld_trns_sol
   real , dimension ( n_chn )  :: cld_trns_sat
   real , dimension ( n_chn )  :: cld_sph_alb
   real , dimension ( n_chn , 2 )  :: kernel
   real :: rad_toa , rad_toc	  
   character ( len = 10 ) :: sensor 
   
   pxl % sol_zen = 33.
   pxl % sat_zen = 43.
   pxl % rel_azi = 133.
   pxl % ctt = 228
   pxl % is_water_phase = .true.
   
   state_vec = [1.6, 0.56 ] 
   sensor = 'MODIS-AQUA'
   sensor='VIIRS'
   channels(:4) = [1 , 6,6 , 20 ] 
   alb_sfc(:4) = [ 0.1, 0.04 , 0.1 , 0.1 ]
   air_trans_ac(:4) = [ 0.9, 0.9 , 0.9, 0.9 ]
   
   rad_toa = 0.14
   rad_toc = 0.002
   
   call get_forward_model (   state_vec &
                            , pxl , sensor &
							, channels , alb_sfc &
							, air_trans_ac, fwd_vec &
							, cld_trns_sol &
							 , cld_trns_sat &
							 , cld_sph_alb &
							 , kernel &
							 , rad_toa &
							 , rad_toc )
							   

  print*,fwd_vec
  print*,kernel (:,1)
   print*,kernel (:,2)
  


end program test_new_fwd
