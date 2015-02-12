! $Header: https://svn.ssec.wisc.edu/repos/cloud_team_nlcomp/trunk/nlcomp.f90 8 2014-01-31 08:14:58Z awalther $
program nlcomp_one_pixel_run

   use M_kracken
   
   use nlcomp_retrieval_mod, only: &
      nlcomp_algorithm	&
      ,  nlcomp_output_structure
   
   use  nlcomp_interface_def_mod  , only:  
		   
   implicit none
   real , dimension (20) :: obs, obs_u , alb_sfc , alb_sfc_u , air_trans_ac
   real, dimension(2) ::  state_apr
   real :: rad_abv_cld , rad_sfc
    character( len = 2) :: color_string
   integer :: start_day
   real :: sol_zen , sat_zen , rel_azi, cld_temp , lun_zen
   logical :: snow
   character (len = 100)  :: text
   real :: cod , cps , codu , cpsu
   character(len=20) :: sensor
  type ( nlcomp_output_structure ) :: nlcomp_results
   
   logical :: water_phase
   
   integer :: debug_mode
   character ( len =20) :: host
   
   character ( len = 1024 ) :: ancil_path

!  dummy input
! make cool input later
  
  
   call getenv("HOST",host)
   start_day = 30
    
   call kracken ("cmd" , "-obs1 0.5211 -obs2 0.133 -alb1 0.12 -alb2 0.12 &        
                & -ctt 276. -sol 123.33 -sat 21. -azi 80 -lun 53 &
					 & -apr1 1.0 -apr2 1.0  -snow ""#N#"" -wat ""#N#""  &
					 & -radabv 0.0002 -radsfc 0.14 -tvis 0.8 -dbg 0 &
					 & -obsu1 0.01 -obsu2 0.01 -albu1 0.01 -albu2 0.01  -tnr 0.8  -qq 0.8  ")
					 
    
	
   obs(1) = rget ("cmd_obs1")
   obs(2) = rget ("cmd_obs2")
   alb_sfc(1) = rget ("cmd_alb1")
   alb_sfc(2) = rget ("cmd_alb2")
   cld_temp = rget ("cmd_ctt")
   sol_zen  = rget ("cmd_sol")
   lun_zen = rget ("cmd_lun")
   sat_zen = rget ("cmd_sat")
   rel_azi = rget ("cmd_azi")
   rad_abv_cld = rget("cmd_radabv")
   
   rad_sfc = rget("cmd_radsfc")
   
   alb_sfc_u(1) = rget("cmd_albu1")
   alb_sfc_u(2) = rget("cmd_albu2")
   obs_u(1) = rget("cmd_obsu1")
   obs_u(2) = rget("cmd_obsu2")
   air_trans_ac(1) = rget("cmd_tvis")
   air_trans_ac(2) = rget("cmd_tnr")
 
   
   state_apr(1) = rget ("cmd_apr1")
   state_apr(2) = rget ("cmd_apr2")
   snow = lget("cmd_snow")
  
   water_phase = lget ("cmd_wat")
  
   debug_mode = iget("cmd_dbg")
   sensor = 'VIIRS' 
   !sensor = sensor(:iflen)	
	
   
   if ( debug_mode == 3 ) then 
      print*, 'NLCOMP stand alone processing'
	  print*,'to switch on this information set -dbg 3'
	  print*,'show usage set -dbg 2'
	  print*
	  print*,'input:   '
	  print*,'obs: ', obs(1:2)
	  print*,'sol, sat, rel_azi: ', sol_zen, sat_zen, rel_azi
	  print*,'alb_sfc: ',alb_sfc(1:2)
	  print*,'cloud top temperature: ', cld_temp,trim('K')
	  print*,'snow: ',snow
	  print*,'water phase: ', water_phase
     
      print*,'above cloud transmission: ', air_trans_ac(1:2)
      print*
   end if
   
   
   if ( debug_mode == 2 ) then 
      print*, 'DCOMP stand-alone processing'
	   print*,'to switch on this information set -dbg 2'
	   print*
	   print*,'usage: ...'
	   print*,' this shows the options with the default values if you do not set them: '
	   print*
	   print*,'./nlcomp -obs1 0.5211 -obs2 0.133 -alb1 0.12 -alb2 0.12'    
      print* ,' -ctt 276. -sol 23.33 -sat 21. -azi 80 '
	   print*,'-apr1 1.0 -apr2 1.0 -sen GOES-15'
	   print*,'-radabv 0.0002 -radsfc 0.14 -tvis 0.8 -dbg 0 '
	   print*,'-obsu1 0.01 -obsu2 0.01 -albu1 0.01 -albu2 0.01  -tnr 0.8  -qq 0.8  '
      print*,' boolean keywords :  -snow -wat'
	   print*
   end if
   
   ancil_path = '/data3/Ancil_Data/clavrx_ancil_data/static/luts/cld/'
   ancil_path = '/DATA/Ancil_Data/clavrx_ancil_data/static/luts/cld/' 
   if ( host(1:4) == 'luna' ) ancil_path = '/DATA/Ancil_Data/clavrx_ancil_data/static/luts/cld/' 
   if ( host(1:4) == 'saga' ) ancil_path = '/data/Ancil_Data/clavrx_ancil_data/static/luts/cld/' 
   
   call nlcomp_algorithm ( obs , obs_u , alb_sfc , alb_sfc_u , state_apr , air_trans_ac &
                              & , lun_zen, sat_zen , rel_azi , cld_temp , water_phase &
							  & , rad_abv_cld , rad_sfc &
							  & , nlcomp_results  &
							  & , debug_in = debug_mode , ancil_path = ancil_path )   ! - output
   write(color_string,'(I2)') 43
   text = '============= NLCOMP RESULTS==============='
   print*,achar(27)//'['//color_string//'m '//trim(text)//achar(27)//'[0m'
   print*,' ' , trim(sensor), ' NLCOMP output:  COD:',nlcomp_results%cod ,'REF: ', nlcomp_results%cps
	
  
    
end program nlcomp_one_pixel_run
