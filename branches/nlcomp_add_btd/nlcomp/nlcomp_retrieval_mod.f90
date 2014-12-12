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
      print*,inp % chn ( 20 ) % rad,inp % chn ( 31 ) % rad,inp % chn ( 32 ) % rad
      bt_20 = planck_rad2tmp ( inp % chn ( 20 ) % rad , 'VIIRS' , 20 )
      bt_31 = planck_rad2tmp ( inp % chn ( 31 ) % rad , 'VIIRS' , 31 )
      bt_32 = planck_rad2tmp ( inp % chn ( 32 ) % rad , 'VIIRS' , 32 )
      
      obs_vec ( 3 ) = bt_31 - bt_32
      obs_vec ( 4 ) = bt_20 - bt_31
      
     
      
      obs_u ( 1 ) = inp % chn ( 42 ) % rfl_u
      obs_u ( 2 ) = inp % chn ( 20 ) % rad_u
      obs_u ( 3 ) = obs_vec ( 3 ) * 0.4
      obs_u ( 4 ) = obs_vec ( 4 ) * 0.4
      if ( obs_vec(1) > 1.1 ) obs_u(1) = 100000.    
      ! - observation error cov
        
      obs_crl = 0.7  ! correlation between channel 42 and 20 
      S_m = 0.
      S_m (1,1) = ( max ( obs_u(1) * obs_vec(1) , 0.01 ) ) ** 2
      S_m (2,2) = ( max ( obs_u(2) * obs_vec(2) , 0.01 ) ) ** 2
      S_m (3,3) = ( 1. ) ** 2
      S_m (4,4) = ( 1. ) ** 2
      
      
      S_m (1,2) =  ( obs_u(2) * obs_vec(2) ) * (obs_u(1) * obs_vec(1) ) * obs_crl
      S_m (2,1) =  ( obs_u(2) * obs_vec(2) ) * (obs_u(1) * obs_vec(1) ) * obs_crl
      
      ! - a_priori 
      
      state_apr = inp % state % a_priori
      
      S_a = 0.0
      S_a(1,1)  =  inp % state % a_priori (1) ** 2 
      S_a(1,1)  = 0.8 ** 2
      S_a(2,2) =  0.9 ** 2    
 
      call findinv ( S_a , S_a_inv , 2 , errorflag)
      
      
      ! = forward model components vector TOADD
      S_b = 0.

      S_b(1,1) = ( alb_sfc_u(1) ) ** 2 
      S_b(2,2) = ( alb_sfc_u(2) ) ** 2
      S_b(3,3) =  1.
      S_b(4,4) =  1.
      S_b(5,5) =  1.

      pxl % sol_zen = inp % geo %sol_zen
      pxl % lun_zen = inp % geo %lun_zen
      pxl % sat_zen = inp % geo %sat_zen
      pxl % rel_azi = inp % geo %rel_azi
      pxl % lun_rel_azi = inp % geo %lun_rel_azi
      pxl % ctt = inp % prd %ctt
      pxl % is_water_phase = inp % prd %cph
       
      dcomp_ancil_path = trim ( inp % conf % ancil_path )
            
      air_trans_ac ( 1 ) = inp % chn ( 42 ) % trans_air_abvcld
      air_trans_ac ( 2 ) = inp % chn ( 20 ) % trans_air_abvcld
      
      alb_sfc ( 1) = inp % chn ( 42 ) % alb_sfc 
      alb_sfc ( 2) = inp % chn ( 20 ) % alb_sfc
      
      debug_mode = 1
      
      cod   = MISSING_REAL4_EM
      cps   = MISSING_REAL4_EM
      codu  = MISSING_REAL4_EM
      cpsu  = MISSING_REAL4_EM


      state_vec = state_apr
      

      iteration_idx = 0
         debug_mode = 6
      IF (debug_mode > 4 ) THEN
         PRINT *, "<--- Begin New Retrieval for pixel = "
         PRINT *, "cloud type, phase = ", cld_phase
         PRINT *, "angles = ", pxl%sat_zen,pxl%lun_zen,pxl%lun_rel_azi
         PRINT *, "sfc ref = ", inp % chn(42) % alb_sfc
         PRINT *, "Trans_Ac 42 = ", inp % chn(42) %trans_air_abvcld
         print *, "Trans_Ac 20 = ", inp % chn(20) %trans_air_abvcld
         PRINT *, "y = ", obs_vec
         PRINT *, "Sy = ", S_m
         PRINT *, "x_ap = ",inp % state % a_priori
         PRINT *, "Sa = ",S_a
      END IF
	   
	   Retrieval_Loop : do
         ! -- Update iteration counter
         iteration_Idx = iteration_Idx + 1
	   	
		   sensor ='VIIRS'
        
        ! Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
         call  nlcomp_forward_computation  ( &
               state_vec  &
            , pxl &
            , trim ( sensor ) &
            
            , alb_sfc  &
            , air_trans_ac &
            , obs_fwd &
            , cld_trans_sol &
            , cld_trans_sat &
            , cld_sph_alb &
            , kernel &
            , rad_abv_cld &
            , rad_clear_toc &
            , lut_path = dcomp_ancil_path   ) 
            
            
           
		   ! - define forward model vector
		   ! - first dimesnion : the two channels
		   ! - 1 sfc albedo vis ; 2 -  sfc albedo ir ; 3- rtm error in vis  4 - rtm error in nir 
		   ! -  5 - terrestiral part 
		 
         kernel_b = 0
         kernel_b ( 1, 1) = ( cld_trans_sol(1) * cld_trans_sat(1) ) / ( (1 - cld_sph_alb(1) * Alb_Sfc(1) ) ** 2. )
         kernel_b ( 1, 2) = 0.
         kernel_b ( 1, 3) = 0.04 
         kernel_b ( 1, 4) = 0.
         kernel_b ( 1, 5) = 0.
         
         kernel_b ( 2, 1) = 0.
         kernel_b ( 2, 2) = ( cld_trans_sol(2) * cld_trans_sat(2))/((1 - cld_sph_alb(2) * alb_Sfc(2)) **2. )
         kernel_b ( 2, 3) = 0.
         kernel_b ( 2, 4) = 0.02
         kernel_b ( 2, 5) = 0.05 * obs_vec(2)
         
         kernel_b ( 3, 1) = 0.
         kernel_b ( 3, 2) = ( cld_trans_sol(2) * cld_trans_sat(2))/((1 - cld_sph_alb(2) * alb_Sfc(2)) **2. )
         kernel_b ( 3, 3) = 0.
         kernel_b ( 3, 4) = 0.02
         kernel_b ( 3, 5) = 0.05 * obs_vec(2)
             
         ! - calculate observation error covariance 
         S_y = S_m + matmul (kernel_b, matmul (S_b, transpose (Kernel_B) ) )
         call findinv ( S_y , S_y_inv , 4 , errorflag)
 
		   !--compute Sx error covariance of solution x 
         S_x_inv = S_a_inv + matmul ( transpose ( Kernel ) , matmul ( S_y_inv , Kernel ) )
         
         call findinv ( S_x_inv, S_x , 2 , errorflag )
         
         ! - compute next iteration step
         delta_X = matmul( S_x , &
                &  (matmul( transpose ( kernel ), &
                &   matmul ( S_y_inv , ( obs_vec - obs_fwd ) ) ) +  &
                &   matmul ( S_a_inv , state_apr - state_vec ) ) )
       
		   ! - check for convergence
         Conv_Test = abs ( sum (delta_X * matmul ( S_x_inv , Delta_X ) ) )	
      
        debug_mode = 7
		
         if ( debug_mode > 4 ) then 
            
            print*
            print*,'iter = ',iteration_idx
            print*, 'f = ', obs_fwd
            print*, 'k = ',kernel
            print*, 'Sx = ',S_x
            print*,'delta x =', delta_x
            print*, ' new x =', state_vec + delta_x
            print*, 'conv test = ', conv_test
            print*,'obs_vec: ', obs_vec
            print*,'obs_fwd: ', obs_fwd
            print*,'state vec: ',state_vec
            print*
         end if
           
         max_step_size = 0.5
         delta_dstnc = sqrt ( sum ( delta_x ** 2 ) )
         
         if ( maxval ( abs(delta_x))   >  MAX_STEP_SIZE  ) then
            delta_x = delta_x * max_step_size / delta_dstnc 
         end if
        
         if ( debug_mode > 4 ) then
            print*,'real delta: ', delta_x
            print*,'new x = ', state_Vec + delta_X
         end if
         
         state_Vec = state_Vec + delta_x
        
         if ( conv_test < 0.2 ) then
	
            cod = 10 ** state_vec(1)
            cps = 10 ** state_vec(2)
            codu = S_x ( 1, 1 )
            cpsu = S_x ( 2, 2 )
		      if ( debug_mode > 4 ) then
               print*,' y:                   ', obs_vec(1:2)
               ! Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
               call  nlcomp_forward_computation  ( &
                     state_vec  , pxl &
						   , trim ( sensor ) , alb_sfc  &
						   , air_trans_ac &
						   , obs_fwd &
						   , cld_trans_sol &
						   , cld_trans_sat &
						   , cld_sph_alb &
						   , kernel , rad_abv_cld ,  rad_clear_toc, lut_path = dcomp_ancil_path   ) 
               print*,'fwd: ',obs_fwd
			  
		      end if
            exit retrieval_loop
		 
         end if
        
          if ( state_vec(1) > 2.0 .and. iteration_idx > 6 ) then
            state_vec(2) = thick_cloud_cps ( obs_vec(2) , pxl )  
             nlcomp_out % statusOK = .true.
             nlcomp_out % cod = 10**2.2
            nlcomp_out % codu = 1.0
             nlcomp_out % cps = 10 ** state_vec(2)
             nlcomp_out % cpsu = 1.0
            
            exit
         end if
		
         if ( iteration_idx > 20 ) then
            cod = -999.
            cps = -999.
            exit retrieval_loop	
         end if

      end do Retrieval_Loop
     
      nlcomp_out % cod = cod
      nlcomp_out % cps = cps
      
      
     
   end subroutine nlcomp_algorithm
end module nlcomp_retrieval_mod
