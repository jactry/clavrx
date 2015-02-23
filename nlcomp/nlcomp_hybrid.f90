module nlcomp_hybrid
   
   use dcomp_lut_mod , only : &
       lut_obj , lut_output

   type, public :: pixel_vec
      real :: sol_zen
	   real :: sat_zen
	   real :: rel_azi
      real :: lun_zen
      real :: lun_rel_azi
	   real :: ctt
      logical :: is_water_phase
   end type pixel_vec
    

contains
   subroutine vis_channel_cod ( &
        rfl_dnb &
      , pixel_vec &
      , cod
      )
      
   real, intent(in) :: rfl_dnb
   type ( pixel_vec ) , intent ( in) :: pixel
   
   real, intent(out) :: cod
   
   call lut_obj % initialize ( 'VIIRS' , ancil_path = lut_path)
   call lut_obj % set_angles ( pixel % sat_zen , pixel % lun_zen $
            , pixel % lun_rel_azi )   
      
      
   
   end subroutine vus_channel_cod


end module nlcomp_hybrid
