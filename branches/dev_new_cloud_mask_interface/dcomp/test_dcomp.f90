! $Id$
!
!   several tests
!

program test_dcomp
   use dcomp_tools, only: &
      findinv
	  
  
   
   use dcomp_mod

   implicit none
   real , dimension (20) :: obs, obs_u , alb_sfc , alb_sfc_u , air_trans_ac
   real, dimension(2) ::  state_apr
   real :: rad_abv_cld , rad_sfc
   character( len = 2) :: color_string
   integer :: start_day
   real :: sol , sat , azi, cld_temp
   logical :: snow
   character (len = 100)  :: text
   
   character(len=20) :: sensor
   integer :: iflen , ier
   
   logical :: water_phase
   
   type ( dcomp_output_structure ) :: dcomp_results
   
   
   integer :: dcomp_mode
   integer :: debug_mode
   character ( len =20) :: host
   
   character ( len = 1024 ) :: ancil_path
   integer :: i , j
   real :: cod_erg(1000,1000)
   real :: ref_erg(1000,1000)
   
   call getenv("HOST",host)
   
   
    obs(1:2) =[ 0.5211 , 0.133]
    alb_sfc(1:2) =[0.12,0.12]
         
    cld_temp  = 276. 
    sol =  23.33 
    sat = 21. 
    azi = 80 
	 state_apr = [1.0,1.0] 
    
    sensor = 'GOES-15' 
    snow = .false. 
    water_phase = .true. 
    rad_abv_cld = 0.0002 
    rad_sfc =  0.14 
    air_trans_ac(1:2) = [0.9,0.8]
    
    debug_mode = 0
	obs_u(1:2) = [ 0.01 ,0.01]
   alb_sfc_u(1:2) = [0.01,0.01]
   
   dcomp_mode =  3 
  
  
   

   ancil_path = '/home/wstraka/geocat/data_algorithms/baseline_cloud_micro_day/version_1/'
   if ( host(1:4) == 'luna' ) ancil_path = '/DATA/Ancil_Data/clavrx_ancil_data/luts/cld/' 
   if ( host(1:4) == 'saga' ) ancil_path = '/data/Ancil_Data/clavrx_ancil_data/luts/cld/' 
   if ( host(1:4) == 'odin' ) ancil_path = '/data3/Ancil_Data/clavrx_ancil_data/luts/cld/' 


   do i=1,1000
      do j=1,1000
      obs(1) = i/1000.
      obs(2) = j/1000.
      call dcomp_algorithm ( obs , obs_u , alb_sfc , alb_sfc_u , state_apr , air_trans_ac &
                              & , sol, sat , azi , cld_temp , water_phase &
							  & , rad_abv_cld , rad_sfc , sensor &
							  & , dcomp_results , dcomp_mode = dcomp_mode &
							  & , debug_in = debug_mode , ancil_path = ancil_path ) 
                       
      
      print*, i, j, 'COD:' &
            ,dcomp_results%cod ,'REF: ',dcomp_results % cps 
      cod_erg(i,j) =    dcomp_results%cod
       ref_erg(i,j) =    dcomp_results%cps                 
   end do
   end do
   
   open(10,file='cod_result.dat')
   write(10,*) cod_erg
   close (10)
   open(10,file='ref_result.dat')
   write(10,*) ref_erg
   close (10)
   
   
end program test_dcomp
