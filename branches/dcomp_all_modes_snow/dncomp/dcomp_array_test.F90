! $Id$
!
!  PURPOSE: test program for array DCOMP 
!
!
program dcomp_array_test
   
   use dncomp_interface_def_mod, only: &
      dncomp_in_type &
      , dncomp_out_type
   
   
   type (dncomp_in_type)  :: inp
   type (dncomp_out_type) :: outp
   integer :: dim_1, dim_2
   integer, parameter :: N_CHN = 45
   logical :: chan_on ( N_CHN ) = .false.
   
   
   
   dim_1 = 2
   dim_2 = 2
   chan_on(:) = .true.
   inp = dncomp_in_type ( dim_1, dim_2, chan_on )
   
   inp % mode = 1
   
   inp % sensor_wmo_id = 224
   inp % lut_path = "/DATA/Ancil_Data/clavrx_ancil_data/static/luts/cld/"
   
   i =1 
   
   inp % refl(i) % d = reshape([20.,40.,50.,70.],(/2,2/))
  
   
   
   do i = 2, 19
      inp % refl(i) % d = 70.
      inp % alb_sfc ( i ) % d = 0.0      
   end do
   
   
   
   inp % sat % d = 22.
   inp % sol % d = 22.
   inp % azi % d = 22.
   
            ! - Cloud products
   inp % cloud_press % d = 890.
   inp % cloud_temp % d  = 266.
   inp % tau_acha % d    = 0.3
   inp % cloud_mask % d  = 3
   inp % is_valid % d = .TRUE.
   
   print*,'hallo array test'
   
   call dcomp_array_loop ( inp, outp)
   
   
   print*,outp % cod % d

end program dcomp_array_test
