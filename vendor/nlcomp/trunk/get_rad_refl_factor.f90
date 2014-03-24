! $Header$
function get_rad_refl_factor ( sensor , sol_zen  ) result ( rad_to_refl )
   implicit none
   character (len =  * )  :: sensor
   real, intent ( in ) :: sol_zen
  
   real :: sun_earth_distance
   real, parameter :: PI = 3.14159265
   real, parameter :: DTOR = PI / 180.
   real :: rad_to_refl
   real :: solar, ew , solar_rad_20
   integer :: start_day
    
      start_day = 100
   
   
!- these are all sensor specifiv values and van be removed an computed in dcomp box!
   sun_earth_distance = 1.0 - 0.016729 * cos ( 0.9856 * ( start_day - 4.0) * DTOR )
   
   
  
   select case ( sensor )
   case ( 'GOES-11')
      solar = 2.2609792
	  ew = 155.07054   
   case ( 'GOES-12')
      solar = 4.204
	  ew = 270.43962  
   case ( 'GOES-13')
      solar = 3.2502353
	  ew = 222.80349
   case ( 'GOES-15')
      solar = 3.5111240
	  ew = 241.02064
   case ( 'MODIS-AQUA' , 'MODIS-TERRA')
       solar = 1.9553725
	   ew = 127.12892
   case ( 'VIIRS')
      solar = 2.2671891
	  ew = 140.86442
  
   case ('NOAA-07')
      solar =  4.5634907 
	  ew = 287.01810
   case ( 'NOAA-08')
      solar = 4.0794900
	  ew = 262.90717	  	   
   case ( 'NOAA-09')
      solar = 4.611
	  ew = 288.84289   
   case ('NOAA-10')
      solar =  4.2889941 
	  ew = 272.12618
   case ( 'NOAA-11')
      solar = 4.448
	  ew = 278.85792	  	   
   case ( 'NOAA-12')
      solar = 4.204
	  ew = 270.43962
   case ('NOAA-14')
      solar =  4.448 
	  ew = 284.69366
   case ( 'NOAA-15')
      solar = 3.781
	  ew = 236.53016    
   case ( 'NOAA-16')
      solar = 3.7372757
	  ew = 236.38144
   case ( 'NOAA-17')
      solar = 4.2443861
	  ew = 269.79606	  	  	   
   case ( 'NOAA-18')
      solar = 4.0620453
	  ew = 259.57508 
   case ( 'NOAA-19')
      solar = 4.1725582
	  ew = 265.07816   
   case ('MTSAT-2') 
      solar = 5.1190589
	  ew = 322.06623 
   case ('MTSAT-1R')
      solar = 4.3693553
	  ew = 280.93528
   case('METOP-A')
      solar = 4.642
	  ew = 291.06
   case('METOP-B')
      solar = 3.9731433
	  ew = 253.16404	 
   case ('Meteosat-8')
      solar = 5.3444818
	  ew = 365.59826
   case('Meteosat-9')
      solar = 5.4832609
	  ew = 374.42649
   case('Meteosat-10')
      solar = 5.4694165
	  ew = 374.16782	   
	      
   case default
      print*,'missing sensor calibrarytion'
	  print*, 'add to get_rad_refl_factor.f90'
      stop
   end select 	
   
   
   solar_rad_20 =  1000.* solar / ew
   rad_to_refl    = PI / cos (sol_zen * DTOR )/ solar_rad_20 / sun_earth_distance ** 2





end function get_rad_refl_factor
