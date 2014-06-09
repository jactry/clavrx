! $Header: https://svn.ssec.wisc.edu/repos/cloud_team_nlcomp/trunk/nlcomp_mod.f90 8 2014-01-31 08:14:58Z awalther $
module nlcomp_retrieval_mod
   private
   public :: nlcomp_algorithm
   
   type , public :: nlcomp_output_structure
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
   end type nlcomp_output_structure

contains

 subroutine nlcomp_algorithm ( obs_vec &
                      & , obs_u &
                      & , alb_sfc &
                      & , alb_sfc_u  &
                      & , state_apr &
                      & , air_trans_ac  &
                      & , sol_zen &
                      & , sat_zen &
                      & , rel_azi &
                      & , cld_temp &
                      & , cld_phase &
                      & , rad_abv_cld &
                      & , rad_clear_toc   &
                      & , nlcomp_out &
                      & , debug_in  &
                      & , ancil_path ) 
                        			 
      use dcomp_math_tools_mod, only: &
         findinv , debug_mode
         
      use nlcomp_forward_mod, only : &
        nlcomp_forward_computation &
        , pixel_vec
   
      implicit none 
      
      real, dimension ( : ) :: obs_vec , obs_u , alb_sfc 
      real, dimension ( : ) :: alb_sfc_u , air_trans_ac  
      real, dimension ( 2 ) :: state_apr 
      real, intent ( in ) :: sol_zen, sat_zen , rel_azi , cld_temp
      logical, intent ( in ) :: cld_phase 
      real, intent ( in ) ::  rad_abv_cld , rad_clear_toc
      
      
      real :: cod, cps, codu, cpsu
      integer, intent(in), optional :: debug_in
      character (len = 1024 ), intent ( in ) , optional :: ancil_path
      character ( len =1024) :: dcomp_ancil_path
      type (nlcomp_output_structure) , intent ( out ) :: nlcomp_out
      
      real, dimension ( 2, 2 ) :: S_a , S_a_inv
      real, dimension ( 2, 2 ) :: S_y , S_y_inv
      real, dimension ( 2, 2 ) :: S_x , S_x_inv
      real, dimension ( 2, 2 ) :: S_m  
      real, dimension ( 5, 5 ) :: S_b 
      real, dimension ( 2, 5 ) :: kernel_b

      real :: obs_crl

      real, dimension ( 2 ) :: state_vec 
      real, dimension ( 2 ) :: delta_x
      real, dimension ( 2, 2) :: kernel
      real, dimension ( 2 ) :: obs_fwd 
      real, dimension ( 2 ) :: cld_trans_sol, cld_trans_sat, cld_sph_alb
	   
      character ( len = 5 ) :: sensor 
      integer :: iteration_idx
      integer :: errorflag 
      real :: conv_test
      real :: delta_dstnc
      real , parameter :: missing_real4_em = -999.
	  
      type ( pixel_vec ) :: pxl
	
      real :: max_step_size
      integer , dimension(2) :: channels

      ! -executable
      pxl % sol_zen = sol_zen
      pxl % sat_zen = sat_zen
      pxl % rel_azi = rel_azi
      pxl % ctt = cld_temp
      pxl % is_water_phase = cld_phase
     
      if ( present ( ancil_path )) then
         dcomp_ancil_path = trim ( ancil_path )
      else
         dcomp_ancil_path = '/data/Ancil_Data/clavrx_ancil_data/'
      end if
    
      debug_mode = 0
      if ( present ( debug_in )) debug_mode = debug_in

      cod = missing_real4_em
      cps = missing_real4_em
      codu = missing_real4_em
      cpsu = missing_real4_em

      S_a = 0.0
      S_a(1,1)  =  state_apr(1) ** 2 
      S_a(1,1)  = 0.8 ** 2
      S_a(2,2) =  0.9 ** 2    
 
      call findinv ( S_a , S_a_inv , 2 , errorflag)
 
      ! - observation error cov
      obs_crl = 0.7  
      S_m = 0.
      S_m (1,1) = ( max ( obs_u(1) * obs_vec(1) , 0.01 ) ) ** 2
      S_m (2,2) = ( max ( obs_u(2) * obs_vec(2) , 0.01 ) ) ** 2
      S_m (1,2) =  ( obs_u(2) * obs_vec(2) ) * (obs_u(1) * obs_vec(1) ) * obs_crl
      S_m (2,1) =  ( obs_u(2) * obs_vec(2) ) * (obs_u(1) * obs_vec(1) ) * obs_crl

      ! = forward model components vector
      S_b = 0.

      S_b(1,1) = ( alb_sfc_u(1) ) ** 2 
      S_b(2,2) = ( alb_sfc_u(2) ) ** 2
      S_b(3,3) =  1.
      S_b(4,4) =  1.
      S_b(5,5) =  1.

      state_vec = state_apr

      iteration_idx = 0

      IF (debug_mode > 4 ) THEN
         PRINT *, "<--- Begin New Retrieval for pixel = "
         PRINT *, "cloud type, phase = ", cld_phase
         PRINT *, "angles = ", sat_zen,sol_zen,rel_azi
         PRINT *, "sfc ref = ", alb_sfc(1:2)
         PRINT *, "Trans_Ac = ", air_trans_ac(1:2)
         PRINT *, "y = ", obs_vec(1:2)
         PRINT *, "Sy = ", S_m
         PRINT *, "x_ap = ",state_apr
         PRINT *, "Sa = ",S_a
      END IF
	   
	   Retrieval_Loop : do
         ! -- Update iteration counter
         iteration_Idx = iteration_Idx + 1
	   	
		   sensor ='VIIRS'
         channels = [1, 20   ]
        ! Start_Time_Point_Hours = COMPUTE_TIME_HOURS()
         call  nlcomp_forward_computation  ( &
            state_vec  , pxl &
            , trim ( sensor ) , channels , alb_sfc  &
            , air_trans_ac &
            , obs_fwd &
            , cld_trans_sol &
            , cld_trans_sat &
            , cld_sph_alb &
            , kernel , rad_abv_cld ,  rad_clear_toc, lut_path = dcomp_ancil_path   ) 

         !  time_acc = time_acc + compute_time_hours() - start_time_point_hours	
		 
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
         
         ! - calculate observation error covariance 
         S_y = S_m + matmul (kernel_b, matmul (S_b, transpose (Kernel_B) ) )
         call findinv ( S_y , S_y_inv , 2 , errorflag)
  
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
		
		
         if ( debug_mode > 4 ) then 
            print*
            print*,'iter = ',iteration_idx
            print*, 'f = ', obs_fwd
            print*, 'k = ',kernel
            print*, 'Sx = ',S_x
            print*,'delta x =', delta_x
            print*, ' new x =', state_vec + delta_x
            print*, 'conv test = ', conv_test
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
						   , trim ( sensor ) , channels , alb_sfc  &
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
