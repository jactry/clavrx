! 
program cx_pfaast_test

use cx_pfaast_coef_mod


type(coef_type) :: coef
character(len=20) :: sensor
character(len=200) :: ancil_path
character (len =9) :: sensor_list(20)
integer :: kban 
integer :: kban_modis = 27 ! show results for modis channel 32

sensor_list(1) = 'AHI'
sensor_list(2) = 'MODIS-AQUA'
sensor_list(3) = 'MODIS-TERRA'
sensor_list(4) = 'VIIRS'
sensor_list(5) = 'GOES-12'
sensor_list(6) = 'GOES-13'
sensor_list(7) = 'GOES-14'
sensor_list(8) = 'GOES-15'
sensor_list(9) = 'GOES-16'

print*,'start..'
sensor='AHI'
!sensor='VIIRS'
!sensor='GOES-16'
sensor='MODIS-AQUA'
ancil_path='/DATA/Ancil_Data/clavrx_ancil_data/'


do i=1,9



   print*
   print*,trim(sensor_list(i))
   call coef % read_it (trim(sensor_list(i)),trim(ancil_path))

  
  
   print*,shape(coef % wvp_liquid)
  ! print*,coef % modis_channel_eqv
   
   kban = minloc ( abs ( coef % modis_channel_eqv - kban_modis ), 1)
   print*,'Native channel: ',kban
   print*,'-----------------------------------------'
   print*,'Max Value wvp liq : ',maxval(coef % wvp_liquid(:,:,kban) )
   print*,'wvp_liq(1,97):',coef%wvp_liquid(1,97,kban)
   print*,'Max Value dry : ',maxval(coef % dry(:,:,kban) )
   print*,'dry(1,97):',coef%dry(1,97,kban)
   
   
end do

end program
