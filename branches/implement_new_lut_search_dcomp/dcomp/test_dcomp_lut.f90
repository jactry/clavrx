program test_dcomp_lut
   use dcomp_data_pool_mod , only : &
       lut_obj , lut_output

   real:: sph,trn_sol,trn_sat,refl,ems
   real :: cod , cps
   integer :: idx_chn , idx_phase
   type ( lut_output) :: lut_data
   character ( len =10) :: sensor
   real :: sat_test

   call lut_obj % initialize( 'VIIRS')
   print*,'set angles'
   call lut_obj % set_angles ( 2.99,44.,33.)
   print*,'get_data'
   cod=1.2
   cps=1.2
   idx_chn = 20
   idx_phase = 1
   call lut_obj % get_data( idx_chn  &
                       & ,  idx_phase , cod, cps, lut_data)

   print*,lut_data%refl,lut_data%albsph,lut_data%trn_sol,lut_data%trn_sat,lut_data%ems, lut_data%dTrans_sat_dcod
   cod=0.5
   cps=1.2
   call lut_obj % get_data( idx_chn  &
                       & ,  idx_phase , cod, cps, lut_data)

   print*,lut_data%refl,lut_data%albsph,lut_data%trn_sol,lut_data%trn_sat,lut_data%ems
   
   call lut_obj % getProperty ( sensor = sensor , sat = sat_test )  
   print*,sensor,sat_test
end program test_dcomp_lut
