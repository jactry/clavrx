! $Id$
!
!  HISTORY: This files name was dcomp_mod.f90
!           changed filename for consistent convention
!
module dcomp_retrieval_mod

   private
   public :: dcomp_algorithm
   
   type , public :: dcomp_output_structure
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
   end type dcomp_output_structure

contains

   subroutine dcomp_algorithm ( &
           obs_vec &
         , obs_u  &
         , alb_sfc & 
         , alb_sfc_u  &
         , state_apr &
         , air_trans_ac   &
         , sol_zen &
         , sat_zen &
         , rel_azi &
         , cld_temp &
         , cld_phase &
         , snow_class &
         , rad_abv_cld &
         , rad_clear_toc  &
         , sensor &
         , output_str &
         , dcomp_mode &
         , ancil_path  &
         , debug_in ) 
   
                     
      use dcomp_math_tools_mod, only: &
         findinv , debug_mode
      
      use dcomp_forward_mod, only: &
         pixel_vec &
         , dcomp_forward_computation &
         , thick_cloud_cps
      
      implicit none 
      
      real, intent(in) :: obs_vec(:) 
      real, intent(in) :: obs_u (:) 
      real, intent(in) :: alb_sfc ( : )
      
      real, intent(in) :: alb_sfc_u (:)
      real, intent(in) :: air_trans_ac (:) 
      real, intent(in) :: state_apr ( 2 )
      real, intent(in) :: sol_zen
      real, intent(in) :: sat_zen
      real, intent(in) :: rel_azi
      real, intent(in) :: cld_temp
      
      logical, intent(in) :: cld_phase 
      integer(kind = 1), intent(in) :: snow_class
      real, intent(in) :: rad_abv_cld 
      real, intent(in) :: rad_clear_toc
      character ( len = * ) , intent ( in ) :: sensor
      integer , intent(in), optional :: dcomp_mode
      integer , intent(in), optional :: debug_in
      
      character (len = 1024 ), intent ( in ) , optional :: ancil_path
      
      type ( dcomp_output_structure ) , intent ( out ) :: output_str
      
      real :: cod , cps , codu , cpsu
      real :: S_a ( 2, 2 ), S_a_inv( 2, 2 )
      real :: S_y ( 2, 2 ), S_y_inv( 2, 2 )
      real :: S_x( 2, 2 ) , S_x_inv( 2, 2 )
      real :: S_m ( 2, 2 ) 
      real :: S_b (5,5) 
      real :: kernel_b(2,5)
      real :: obs_crl
      real :: state_vec (2)
      real :: delta_x (2)
      real :: kernel(2,2)
      real :: obs_fwd (2)
      real :: cld_trans_sol(2)
      real :: cld_trans_sat(2)
      real :: cld_sph_alb(2)
      
      integer :: iteration_idx
      integer :: errorflag 
      real :: conv_test
      real :: delta_dstnc
      real , parameter :: MISSING_REAL4_EM = -999.
      integer :: algo_mode
      type ( pixel_vec ) :: pxl
      
      real :: max_step_size
      integer , dimension(2) :: channels
      real :: cld_albedo_vis
      character ( len=1024) :: dcomp_ancil_path
       
      !     --   executable
      
      ! - initialize output
      output_str % statusOK = .false.
      output_str % cod  = -999.
      output_str % codu = -999.
      output_str % cps = -999.
      output_str % cpsu = -999.
      output_str % cloud_alb_vis = -999.
      output_str % cloud_alb_vis_u = -999.
      output_str % cloud_trans_sol_vis = -999.
      output_str % cloud_trans_sol_vis_u = -999.
      output_str % cloud_trans_sat_vis = -999.
      output_str % cloud_trans_sat_vis_u = -999.
      output_str % cloud_sph_alb_vis = -999.
      output_str % cloud_sph_alb_vis_u  = -999.
      
      ! -- observation
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
      
      algo_mode = 3
      if  ( present ( dcomp_mode )) algo_mode = dcomp_mode
      
      select case ( algo_mode )
         case ( 1 )
            channels = [ 1 , 6 ] 
         case(2)
            channels = [ 1 , 7 ]
         case ( 3 )
            channels = [ 1, 20 ]
         case ( 4)
            channels = [6,20]   
         case default
           print*,'this mode is not set stop'
           stop
      end select 
      
      debug_mode = 0
      if ( present ( debug_in )) debug_mode = debug_in
      
      cod  = missing_real4_em
      cps  = missing_real4_em
      codu = missing_real4_em
      cpsu = missing_real4_em
      
      S_a = 0.0
      S_a(1,1)  =  state_apr(1) ** 2 
      S_a(1,1)  = 0.8 ** 2
      S_a(2,2) =  0.65 ** 2    
      
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
      
      if (debug_mode > 4 ) then
         print *, "<--- begin new retrieval for pixel = "
         print *, "cloud type, phase = ", cld_phase
         print *, "angles = ", sat_zen,sol_zen,rel_azi
         print *, "sfc ref = ", alb_sfc(1:2)
         print *, "trans_ac = ", air_trans_ac(1:2)
         print *, "y = ", obs_vec(1:2)
         print *, "sy = ", s_m
         print *, "x_ap = ",state_apr
         print *, "sa = ",s_a
      end if

      Retrieval_Loop : do
         
         iteration_Idx = iteration_Idx + 1
         
         call  dcomp_forward_computation  ( &
            state_vec  , pxl &
            , trim ( sensor ) , channels , alb_sfc  &
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
         ! - first dimension : the two channels
         ! - 1 sfc albedo vis ; 2 -  sfc albedo ir ; 3- rtm error in vis  4 - rtm error in nir 
         ! -  5 - terrestrial part 

         kernel_b = 0
         kernel_b ( 1, 1) = ( cld_trans_sol(1) * cld_trans_sat(1) ) &
            & / ( (1 - cld_sph_alb(1) * Alb_Sfc(1) ) ** 2. )
         kernel_b ( 1, 2) = 0.
         kernel_b ( 1, 3) = 0.04 
         kernel_b ( 1, 4) = 0.
         kernel_b ( 1, 5) = 0.

         kernel_b ( 2, 1) = 0.
         kernel_b ( 2, 2) = ( cld_trans_sol(2) * cld_trans_sat(2)) &
            & / ( ( 1 - cld_sph_alb(2) * alb_Sfc(2)) ** 2. )
         kernel_b ( 2, 3) = 0.
         kernel_b ( 2, 4) = 0.02
         kernel_b ( 2, 5) = 0.05 * obs_vec(2)

         ! - calculate observation error covariance 
         S_y = S_m + matmul (kernel_b, matmul (S_b, transpose (Kernel_B) ) )
         call findinv ( S_y , S_y_inv , 2 , errorflag)

         !--compute Sx error covariance of solution x 
         S_x_inv = S_a_inv + matmul ( transpose ( Kernel ) , matmul ( S_y_inv , Kernel ) )
         call findinv ( S_x_inv, S_x , 2 , errorflag )

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
      
         ! - control step size 
         max_step_size = 0.5
         delta_dstnc = sqrt ( sum ( delta_x ** 2 ) )
         if ( maxval ( abs(delta_x))   >  max_step_size  ) then
            delta_x = delta_x * max_step_size / delta_dstnc 
         end if

         if ( debug_mode > 4 ) then
            print*,'real delta: ', delta_x
            print*,'new x = ', state_Vec + delta_X
         end if

         state_Vec = state_Vec + delta_x

         if ( conv_test < 0.08 ) then
   
            call  dcomp_forward_computation  ( &
              state_vec &
               , pxl &
               , trim ( sensor ) , channels , alb_sfc  &
               , air_trans_ac &
               , obs_fwd &
               , cld_trans_sol &
               , cld_trans_sat &
               , cld_sph_alb &
               , kernel &
               , rad_abv_cld &
               , rad_clear_toc &
               , lut_path = ancil_path &
               , cld_albedo_vis = cld_albedo_vis ) 
            
            output_str % statusOK = .true.
            output_str % cod  = 10 ** state_vec(1)
            output_str % codu = sqrt ( S_x ( 1, 1 ) )
            output_str % cps  = 10 ** state_vec(2)
            output_str % cpsu = sqrt ( S_x ( 2, 2 ) )
            
            output_str % cloud_alb_vis = cld_albedo_vis
            output_str % cloud_alb_vis_u = -999.
            output_str % cloud_trans_sol_vis = cld_trans_sol(1) 
            output_str % cloud_trans_sol_vis_u = -999.
            output_str % cloud_trans_sat_vis = cld_trans_sat(1)
            output_str % cloud_trans_sat_vis_u = -999.
            output_str % cloud_sph_alb_vis = cld_sph_alb(1)
            output_str % cloud_sph_alb_vis_u  = -999.
      
            if ( debug_mode > 4 ) then
               print*,' y:                   ', obs_vec(1:2)
               print*,'fwd of solution       ', obs_fwd
            end if
            
            if ( state_vec(1) > 2.2  ) then
             state_vec(2) = thick_cloud_cps ( obs_vec(2) , channels(2), pxl, dcomp_mode ) 
             output_str % statusOK = .true.
             output_str % cod = 10**2.2
             output_str % codu = 1.0
             output_str % cps = 10 ** state_vec(2)
             output_str % cpsu = 1.0
            
            exit
         end if
            
            
            
            exit retrieval_loop

         end if
         
         if ( state_vec(1) > 2.0 .and. iteration_idx > 6 ) then
            state_vec(2) = thick_cloud_cps ( obs_vec(2) , channels(2), pxl, dcomp_mode )  
            output_str % statusOK = .true.
            output_str % cod = 10**2.2
            output_str % codu = 1.0
            output_str % cps = 10 ** state_vec(2)
            output_str % cpsu = 1.0
            
            exit
         end if
 
         if ( iteration_idx > 20 )  exit retrieval_loop  
                  
      end do Retrieval_Loop 
      
   end subroutine dcomp_algorithm
   
end module dcomp_retrieval_mod
