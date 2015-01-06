! $Header: https://svn.ssec.wisc.edu/repos/cloud_team_nlcomp/trunk/nlcomp_mod.f90 8 2014-01-31 08:14:58Z awalther $
module nlcomp_retrieval_mod
   private
   public :: nlcomp_algorithm
   
   type , public :: nlcomp_output_type
      logical :: statusOK
      real :: cod
      real :: cps
      real :: codu
      real :: cpsu
      real :: cloud_alb_vis
      real :: cloud_alb_vis_u
      real :: cloud_trans_sol_vis
      real :: cloud_trans_sat_vis
      real :: cloud_trans_sol_vis_u
      real :: cloud_trans_sat_vis_u
      real :: cloud_sph_alb_vis
      real :: cloud_sph_alb_vis_u
   end type nlcomp_output_type
   
   type :: nlcomp_input_chan_type
      real :: rad
      real :: rad_u
      real :: rfl
      real :: rfl_u
      real :: alb_sfc
      real :: alb_sfc_u
      real :: trans_air_abvcld
      real :: rad_abvcld_nwp
      real :: rad_sfc_nwp  
   end type nlcomp_input_chan_type
   
   type :: nlcomp_geo_type
      real :: sol_zen
      real :: sat_zen
      real :: rel_azi
      real :: lun_rel_azi
      real :: lun_zen
      real :: tsfc
   end type nlcomp_geo_type
   
   type :: nlcomp_prd_type
      real :: ctt
      logical :: cph
   end type nlcomp_prd_type
   
   type :: nlcomp_conf_type
      character (len = 1024 ) :: ancil_path
      integer :: debug_in   
   end type nlcomp_conf_type
   
   type :: nlcomp_state_type
      real :: a_priori ( 2 ) 
   end type nlcomp_state_type
   
   type, public :: nlcomp_input_type
      type ( nlcomp_input_chan_type ) :: chn ( 42)
      type ( nlcomp_geo_type ) :: geo
      type ( nlcomp_prd_type ) :: prd
      type ( nlcomp_conf_type) :: conf
      type ( nlcomp_state_type ) :: state
   end type nlcomp_input_type
   
   integer :: N_OBS = 4

contains

 subroutine nlcomp_algorithm ( inp &
                      & , nlcomp_out  ) 
                        			 
      use dcomp_math_tools_mod, only: &
         findinv , debug_mode
         
      use nlcomp_forward_mod, only : &
        nlcomp_forward_computation &
        , pixel_vec, thick_cloud_cps
        
      use clavrx_planck_mod 
      
      use nlcomp_hybrid_mod, only: &
         vis_channel_cod   &
         , cps_known_cod
   
      implicit none 
      
      type ( nlcomp_input_type) , intent(in) :: inp
      type (nlcomp_output_type) , intent ( out ) :: nlcomp_out
      
      real  :: obs_vec (N_OBS) 
      real  :: obs_u (N_OBS)
      
      real  :: alb_sfc (3)
      real  :: alb_sfc_u (3) 
      real  :: air_trans_ac (3) 
      real  :: state_apr (2)
      real  :: sol_zen, sat_zen , rel_azi , cld_temp
      logical :: cld_phase 
      real :: rad_abv_cld (42)
      real :: rad_clear_toc (42)
      
      
      real :: cod, cps, codu, cpsu
      integer :: debug_in
      character (len = 1024 )  :: ancil_path
      character (len = 1024 ) :: dcomp_ancil_path
     
      
      real, dimension ( 2, 2 ) :: S_a , S_a_inv
      real, dimension ( N_OBS, N_OBS ) :: S_y , S_y_inv
      real, dimension ( 2, 2 ) :: S_x , S_x_inv
      real, dimension ( N_OBS, N_OBS ) :: S_m  
      real, dimension ( 5, 5 ) :: S_b 
      real, dimension ( N_OBS, 5 ) :: kernel_b

      real :: obs_crl

      real, dimension ( 2 ) :: state_vec 
      real, dimension ( 2 ) :: delta_x
      real, dimension ( 4, 2) :: kernel
      real, dimension ( N_OBS ) :: obs_fwd 
      real, dimension ( 2 ) :: cld_trans_sol, cld_trans_sat, cld_sph_alb
	   
      character ( len = 5 ) :: sensor 
      integer :: iteration_idx
      integer :: errorflag 
      real :: conv_test
      real :: delta_dstnc
      real , parameter :: missing_real4_em = -999.
	  
      type ( pixel_vec ) :: pxl
	
      real :: max_step_size
     
      
      real :: bt_20
      real :: bt_31
      real :: bt_32

      ! - executable

     
      ! - okay lets start
      
      rad_abv_cld (20)  = inp % chn ( 20 ) % rad_abvcld_nwp
      rad_abv_cld (31)  = inp % chn ( 31 ) % rad_abvcld_nwp
      rad_abv_cld (32)  = inp % chn ( 32 ) % rad_abvcld_nwp
      
      rad_clear_toc(20) = inp % chn ( 20 ) % rad_sfc_nwp
      rad_clear_toc(31) = inp % chn ( 31 ) % rad_sfc_nwp
      rad_clear_toc(32) = inp % chn ( 32 ) % rad_sfc_nwp
      
      ! - first define the observation vector
      
      ! -  observation vector
      
      obs_vec ( 1 ) = inp % chn ( 42 ) % rfl
      obs_vec ( 2 ) = inp % chn ( 20 ) % rad
       
       
      
      pxl % sol_zen = inp % geo %sol_zen
      pxl % lun_zen = inp % geo %lun_zen
      pxl % sat_zen = inp % geo %sat_zen
      pxl % rel_azi = inp % geo %rel_azi
      pxl % lun_rel_azi = inp % geo %lun_rel_azi
      pxl % ctt = inp % prd %ctt
      pxl % is_water_phase = inp % prd %cph
      dcomp_ancil_path = trim ( inp % conf % ancil_path )
      alb_sfc ( 1) = inp % chn ( 42 ) % alb_sfc 
      

      bt_20 = planck_rad2tmp ( inp % chn ( 20 ) % rad , 'VIIRS' , 20 )
      bt_31 = planck_rad2tmp ( inp % chn ( 31 ) % rad , 'VIIRS' , 31 )
      bt_32 = planck_rad2tmp ( inp % chn ( 32 ) % rad , 'VIIRS' , 32 )
      
      obs_vec ( 3 ) = bt_31 - bt_32
      obs_vec ( 4 ) = bt_20 - bt_31
      
      print*,'bt s 20 31 32',bt_20,bt_31,bt_32, inp % geo % tsfc
      print*, 
      cod = -999.
      cps= -999.
      
      if ( obs_vec ( 1 ) > 1.0 ) return
      call vis_channel_cod (obs_vec ( 1 ),  pxl, alb_sfc ( 1) , dcomp_ancil_path , cod) 
      print*,'------  coputed cod ------- '
      call cps_known_cod ( obs_vec ,  cod , rad_clear_toc,  pxl , dcomp_ancil_path , cps ,inp % geo % tsfc )
      
     
      
      
      if ( cod > 0. ) then  
         nlcomp_out % cod = 10**cod
         nlcomp_out % cps = 10**cps
      end if
     
      
     
   end subroutine nlcomp_algorithm
end module nlcomp_retrieval_mod
