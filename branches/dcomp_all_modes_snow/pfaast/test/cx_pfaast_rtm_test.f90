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
   character (len =40 ) :: sensor_list(3)
   integer :: ii, jj
   
   ancil_data_path = '/DATA/Ancil_data/clavrx_ancil_data/'
   temp = tstd
   wvmr = wstd
   ozmr = ostd
   theta = 0.
   sensor='VIIRS'
   sensor='MODIS-AQUA'
   sensor = 'GOES-16'
   kban_in = 32
   use_modis_channel_equivalent  = .true.
   
   sensor_list(1) =  'VIIRS     '
   sensor_list(2) =  'MODIS-AQUA'
   sensor_list(3) =  'GOES-16   '
   
   do ii =  20,33
      print*,'+++++++++++++++++++++++'
      print*,'channel: ',ii
      do jj = 1,3
   
         call compute_transmission_pfaast ( &
            ancil_data_path &
            & ,temp &
            & ,wvmr &
            & ,ozmr & 
            & ,theta  &
            & ,sensor_list(jj) &
            & ,ii &
            & ,taut &
            & , use_modis_channel_equivalent )
            
           
            write(*,'(2A,2x,f7.4)')  'Tottransm for ' &
               , trim(sensor_list(jj)), taut(101)
           
            
       
         end do
     end do

end program cx_pfaast_rtm_test
