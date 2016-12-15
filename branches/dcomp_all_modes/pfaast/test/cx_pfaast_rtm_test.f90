! $Id$
program cx_pfaast_rtm_test
   use cx_pfaast_mod
   use cx_pfaast_constants_mod
   
   integer, parameter :: N_PROFILE = 101
   character(len = 200 ) :: ancil_data_path  
   real  :: temp (N_PROFILE)
   real  :: wvmr(N_PROFILE)
   real  :: ozmr(N_PROFILE)
   real  :: theta 
   character (len =40 ) :: sensor
   integer  :: kban_in 
   logical :: use_modis_channel_equivalent
   real  :: taut (N_PROFILE)  
   
   
   
   ancil_data_path = '/DATA/Ancil_data/clavrx_ancil_data/'
   temp = tstd
   wvmr = wstd
   ozmr = ostd
   theta = 62.
   sensor='VIIRS'
   kban_in = 32
   use_modis_channel_equivalent  = .true.
   
   
   
   call compute_transmission_pfaast ( &
       ancil_data_path &
       & ,temp &
       & ,wvmr &
       & ,ozmr & 
       & ,theta  &
       & ,sensor &
       & ,kban_in &
       & ,taut &
       & , use_modis_channel_equivalent )
   
       print*,taut

end program cx_pfaast_rtm_test
