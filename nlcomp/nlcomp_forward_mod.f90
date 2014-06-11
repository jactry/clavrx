! $Id: dcomp_forward_mod.f90 385 2014-06-05 20:54:18Z awalther $
module nlcomp_forward_mod

   use dcomp_lut_mod , only : &
       lut_obj , lut_output
 
   
    use dcomp_science_tools_mod, only: &
      get_planck_radiance_39um  &
      , get_rad_refl_factor                
   
    use clavrx_planck_mod
   private
   public :: thick_cloud_cps
   
   
   integer, parameter, public:: real4 = selected_real_kind(6,37)
   
   real , parameter :: PI = 4.* ATAN(1.)
   real , parameter :: DTOR = PI / 180.
   
   type, public :: pixel_vec
      real :: sol_zen
	   real :: sat_zen
	   real :: rel_azi
      real :: lun_zen
      real :: lun_rel_azi
	   real :: ctt
      logical :: is_water_phase
   end type pixel_vec
    
   public :: nlcomp_forward_computation
   public :: compute_cld_albedo_vis
   
   real , save :: planck_rad20
   real , save :: rad_to_refl 
   
contains
   
   subroutine nlcomp_forward_computation ( &
                                    state_vec  &   ! -input
                                  , pixel &
                                  , sensor &
                                  
                                  , alb_sfc &
                                  , air_trans_ac &
                                  , fm_vec &         ! - output
                                  , cld_trns_sol &
                                  , cld_trns_sat &
                                  , cld_sph_alb &
                                  , kernel & 
                                  , rad_abv_cld &
                                  , rad_clear_toc  &
                                  , lut_path &
                                  , cld_albedo_vis  )
       
       
                          
   
      implicit none
      
      real , intent(in) :: state_vec (2)
      type ( pixel_vec ) , intent ( in ) :: pixel
      
      real,  intent ( in ) :: air_trans_ac ( : )   
      character ( len = * ) , intent (in ) :: sensor     
      real, intent(in) ::  alb_sfc (2)  
      real, optional, intent ( in ) ::   rad_abv_cld (42)
      real, optional, intent ( in ) ::   rad_clear_toc(42) 
      character ( len = 1024 ) , intent ( in ) , optional :: lut_path  
         
      real , intent ( out ) :: fm_vec ( 4 )
      real , intent ( out ) :: kernel( 4 , 2 )
      real , intent ( out ) :: cld_trns_sol( 2 ) 
      real , intent ( out ) :: cld_trns_sat( 2 )
      real , intent ( out ) :: cld_sph_alb( 2 )
      real , intent ( out ) , optional :: cld_albedo_vis
    
 
      integer :: phase_num
      
      real :: fm_nir_terr 
      real :: kernel_nir_terr_cod
      real :: kernel_nir_terr_cps         
      real :: trans_two_way
      real :: air_mass_two_way
      real :: air_mass_sat
      real :: trans_abv_cld
         
	   type ( lut_output) :: lut_data , lut_data31, lut_data32 , lut_data20

     
      
      integer :: MODIS_EQUIVALENT_39_CHN = 20
	   real :: Alb_Sfc_Term
      
      integer , parameter :: VIS_CHN = 42
      real :: bt31, bt32, bt20
      real :: rad31, rad32, rad20
      real :: planck_rad31, planck_rad32
      
      
     
	   ! - executable
      
      
      ! initialize lut object
      call lut_obj % initialize ( 'VIIRS' , ancil_path = lut_path)
      call lut_obj % set_angles ( pixel % sat_zen , pixel % lun_zen, pixel % lun_rel_azi )
      
      
      
      phase_num = 2
      if ( pixel % is_water_phase ) phase_num = 1
      
      
      air_mass_two_way = ( 1. / cos (pixel % lun_zen * DTOR ) + 1./ cos ( pixel % Sat_zen * DTOR ) )
      planck_rad20 = planck_tmp2rad ( pixel % ctt, trim(sensor), 20)
      planck_rad31 = planck_tmp2rad ( pixel % ctt, trim(sensor), 31)
      planck_rad32 = planck_tmp2rad ( pixel % ctt, trim(sensor), 32)
      
      
      ! forward element1
      
      call lut_obj % get_data ( 1 , phase_num , state_vec(1), state_vec(2) , lut_data)
     
      cld_trns_sol ( 1 ) = lut_data % trn_sol
      cld_trns_sat ( 1 ) = lut_data % trn_sat
      cld_sph_alb  ( 1 )  = lut_data % albsph
      !cld_albedo_vis  =  lut_data % alb
      
      Alb_Sfc_Term = max (0. , Alb_Sfc( 1 ) / (1.0 - Alb_Sfc(1 ) * lut_data % albsph )) 
      fm_vec(1) = lut_data% Refl + Alb_Sfc_Term * lut_data%Trn_sol * lut_data%Trn_sat 
       
      kernel(1,1) = lut_data%dRefl_dcod &
                     & +  Alb_Sfc_Term * lut_data%Trn_sol *  lut_data % dTrans_sat_dcod &
                     & + Alb_Sfc_Term * lut_data%Trn_sat  * lut_data%Dtrans_sol_Dcod &
                     & +  ((lut_data%Trn_sol * lut_data%Trn_sat  * Alb_Sfc( 1 ) * Alb_Sfc( 1 ) * lut_data%Dsph_alb_Dcod) &
                    /(( 1 - Alb_Sfc(1 ) * lut_data%albsph)**2))  
                                           
      kernel(1,2) = lut_data%dRefl_dcps &
                     & +  Alb_Sfc_Term * lut_data%Trn_sol *  lut_data % dTrans_sat_dcps &
                     & + Alb_Sfc_Term * lut_data%Trn_sat  * lut_data%Dtrans_sol_Dcps &
                     & +  ((lut_data%Trn_sol * lut_data%Trn_sat  * Alb_Sfc( 1 ) * Alb_Sfc( 1 ) * lut_data%Dsph_alb_Dcps) &
                    /(( 1 - Alb_Sfc(1 ) * lut_data%albsph)**2))  
                    
      trans_two_way = air_trans_ac( 1 ) ** air_mass_two_way   
      
      fm_vec(1) = fm_vec(1) * trans_two_way
	   kernel(1,1) = kernel(1,1) * trans_two_way
      kernel(1,2) = kernel(1,2) * trans_two_way
      
      ! - forward element2
     
      call lut_obj % get_data ( 20, phase_num , state_vec(1), state_vec(2) , lut_data20)
     
      

	   air_mass_sat =  1./ cos ( pixel % Sat_zen * dtor) 
		trans_abv_cld = air_trans_ac(2) ** air_mass_sat
                                   
      fm_nir_terr =   (rad_abv_cld(20) * lut_data20 % ems &
                  + trans_abv_cld * lut_data20 %  ems * planck_rad20 &
				  + lut_data20 % trn_ems * rad_clear_toc(20) )
            
      kernel_nir_terr_cod =  (rad_abv_cld(20) * lut_data20 %  Dems_Dcod & 
               + trans_abv_cld * Planck_Rad20 * lut_data20 % Dems_Dcod  &
               + rad_clear_toc(20) *lut_data20 % Dtrnems_Dcod)
             
      kernel_nir_terr_cps  =   (rad_abv_cld(20) * lut_data20 % Dems_Dcps & 
                + trans_abv_cld * Planck_Rad20 * lut_data20 % Dems_Dcps  &
                + rad_clear_toc(20)  *lut_data20 % Dtrnems_Dcps )
                       
      fm_vec(2) = fm_nir_terr
      kernel ( 2, 1) = kernel_nir_terr_cod 
      kernel ( 2, 2) = kernel_nir_terr_cps 
            
     
  
      call lut_obj % get_data ( 31, phase_num , state_vec(1), state_vec(2) , lut_data31)
      call lut_obj % get_data ( 32, phase_num , state_vec(1), state_vec(2) , lut_data32)  
      call lut_obj % get_data ( 20, phase_num , state_vec(1), state_vec(2) , lut_data20)
      rad20 = lut_data20 %  ems * planck_rad20 + lut_data20 % trn_ems * rad_clear_toc(20)            
      rad31 = lut_data31 %  ems * planck_rad31 + lut_data31 % trn_ems * rad_clear_toc(31)
      rad32 = lut_data32 %  ems * planck_rad32 + lut_data32 % trn_ems * rad_clear_toc(32)  
      bt20 =  planck_rad2tmp ( rad20, trim(sensor), 20)    
      bt31 =  planck_rad2tmp ( rad31, trim(sensor), 31)
      bt32 =  planck_rad2tmp ( rad32, trim(sensor), 32)  
      
      
       ! element 3
   
      fm_vec(3) = bt31 - bt32
      rad20_dcod = lut_data20 %  ems * planck_rad20 + lut_data20 % trn_ems * rad_clear_toc(20)
      bt20_dum =  planck_rad2tmp ( rad20 + , trim(sensor), 20)
      bt31_dum =  planck_rad2tmp ( rad31, trim(sensor), 31)
            
      kernel ( 3,1) = (fm_vec(3) - ( bt31 - bt32) ) / 0.01
      
      call lut_obj % get_data ( 31, phase_num , state_vec(1), state_vec(2) + 0.01 , lut_data31)
      call lut_obj % get_data ( 32, phase_num , state_vec(1), state_vec(2) + 0.01 , lut_data32)
      
      rad31 = lut_data31 %  ems * planck_rad31 + lut_data31 % trn_ems * rad_clear_toc(31)
      rad32 = lut_data32 %  ems * planck_rad32 + lut_data32 % trn_ems * rad_clear_toc(32)
    
      bt31 =  planck_rad2tmp ( rad31, trim(sensor), 31)
      bt32 =  planck_rad2tmp ( rad32, trim(sensor), 32) 
                         
      kernel ( 3,2) =   (fm_vec(3) - ( bt31 - bt32) ) / 0.01                   
                          
      
      ! element 4
      
      
      
     
      call lut_obj % get_data ( 31, phase_num , state_vec(1), state_vec(2) , lut_data31)
      
      rad31 = lut_data31 %  ems * planck_rad31 + lut_data31 % trn_ems * rad_clear_toc(31)
            
      bt31 =  planck_rad2tmp ( rad31, trim(sensor), 31)
      
   
      fm_vec(4) = bt20 - bt31
      
      call lut_obj % get_data ( 31, phase_num , state_vec(1)+0.01, state_vec(2) , lut_data31)
      call lut_obj % get_data ( 20, phase_num , state_vec(1)+0.01, state_vec(2) , lut_data20)
      rad31 = lut_data31 %  ems * planck_rad31 + lut_data31 % trn_ems * rad_clear_toc(31)
      rad20 = lut_data20 %  ems * planck_rad20 + lut_data20 % trn_ems * rad_clear_toc(20) 
      bt31 =  planck_rad2tmp ( rad31, trim(sensor), 31)
      bt20 =  planck_rad2tmp ( rad20, trim(sensor), 20)  
      
      kernel ( 4,1) = fm_vec(4) - ( bt20 - bt31 ) / 0.01 
      
      call lut_obj % get_data ( 31, phase_num , state_vec(1), state_vec(2)+0.01 , lut_data31)
      call lut_obj % get_data ( 20, phase_num , state_vec(1), state_vec(2)+0.01 , lut_data20)
      rad31 = lut_data31 %  ems * planck_rad31 + lut_data31 % trn_ems * rad_clear_toc(31)
      rad20 = lut_data20 %  ems * planck_rad20 + lut_data20 % trn_ems * rad_clear_toc(20) 
      bt31 =  planck_rad2tmp ( rad31, trim(sensor), 31)
      bt20 =  planck_rad2tmp ( rad20, trim(sensor), 20)  
                         
      kernel ( 4,2) = fm_vec(4) - ( bt20 - bt31 ) / 0.01 
      
           
        
	   
        
   end subroutine nlcomp_forward_computation
    
   !------------------------------------------------------------------------------
   !
   !------------------------------------------------------------------------------
   function thick_cloud_cps ( rfl_nir_obs,  pixel  ) result ( cps)
      implicit none
   
      real,  intent(in) :: rfl_nir_obs
      type ( pixel_vec) , intent(in) :: pixel
      
      
      
      
      integer :: channel_nir = 20
      real :: cps
      
      ! -- local
      
      integer  :: cps_pos          ! position in LUT
      real     :: cps_wgt ! position weights in LUT
      integer :: i
      real :: cps_vec_lut (9)   !- CPS vector in LUT      
      real :: rad (9) 
      real :: rfl_terr (9) 
     
      real :: Ems_lut(9)  ! - sub table
      real :: rfl_lut(9)
      real :: Rfl (9)						! - Reflectance at TOA under consideration of surface impact
      integer :: ii
      integer :: phase_num
      
      ! - executable
      
      ! - cps lut range [0.4,0.6,...,2.0]
      cps_vec_lut = (/(i,i=4,20,2)/)/10.
      
      phase_num = 2
      if ( pixel % is_water_phase ) phase_num = 1
      
      call lut_obj % thick_cloud_rfl ( channel_nir  , phase_num  , rfl_lut, ems_lut )
   
      rfl = rfl_lut  
           
         rad =  ems_lut * Planck_Rad20 
         rfl_terr = rad * Rad_To_Refl
         rfl = rfl +  rfl_terr
     
      
      cps_pos = -999
      cps = -999.
 
      do ii = 1, 8
         if ( rfl(ii)-rfl_nir_obs > 0 .and. rfl(ii+1)-rfl_nir_obs < 0 ) cps_pos = ii
      end do

      if (cps_pos .NE. -999) then
         cps_wgt = (rfl(cps_pos) - rfl_nir_obs)/(rfl(cps_pos) - rfl(cps_pos+1))
         cps = CPS_Vec_LUT(cps_pos) + (cps_wgt * 0.2)
      end if
 
      
      return					 
      
      
   end function thick_cloud_cps
 
       
end module nlcomp_forward_mod
