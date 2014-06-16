! $Id$

module dcomp_science_tools_mod

   public :: get_rad_refl_factor
   public :: get_planck_radiance_39um

contains


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
      case ( 'GOES-08')
         solar = 2.1734486
         ew = 149.27568
      case ( 'GOES-09')
         solar = 2.2212146
         ew = 152.97360
      case ( 'GOES-10')
         solar = 2.2627592
         ew = 156.03971
      case ( 'GOES-11')
         solar = 2.2609792
         ew = 155.07054   
      case ( 'GOES-12')
         solar = 4.204
         ew = 270.43962  
      case ( 'GOES-13')
      solar = 3.2502353
      ew = 222.80349
   case ( 'GOES-14')
      solar = 3.4876949 
      ew = 237.43905 
   case ( 'GOES-15')
      solar = 3.5111240
      ew = 241.02064
   case ( 'MODIS-AQUA' , 'MODIS-TERRA')
       solar = 1.9553725
       ew = 127.12892
   case ( 'VIIRS')
      solar = 2.2671891
      ew = 140.86442
    case ('NOAA-05','TIROS-N')
      solar =  4.1604741 
      ew = 267.14124
   case ( 'NOAA-06')
      solar = 3.9980827
      ew = 254.00470
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
   case('GOES-16','ABI')
      solar = 3.5111240 ! faked from goes-15
      ew = 241.02064 ! faked
	case('COMS-1')
      solar = 4.8461549
      ew = 306.29122
         
   case default
      print*,'missing sensor calibration'
	  print*, 'add to get_rad_refl_factor.f90'
      stop
   end select 	
   
   
   solar_rad_20 =  1000.* solar / ew
   rad_to_refl    = PI / cos (sol_zen * DTOR )/ solar_rad_20 / sun_earth_distance ** 2





end function get_rad_refl_factor


   function get_planck_radiance_39um ( tmp , sensor ) result ( rad )
      implicit none
      real, intent ( in) :: tmp
      character (len =  * ) :: sensor
      !real, intent ( out), optional :: db_dt
   
   
   
      character(len = 20 ) ,save :: sensor_saved
      real :: rad
	  
      real :: db_dt_tmp
      integer, parameter:: nplanck = 161
      real, dimension(nplanck) :: B20
      real, parameter :: T_planck_min = 180.0  
      real, parameter :: delta_T_planck = 1.0
      real, dimension(nplanck) , save :: T_planck
      integer :: l , i
      
      real :: nu_20 , a1_20 , a2_20
      real , parameter :: c1 = 1.191062e-5
      real , parameter :: c2 = 1.4387863
      real :: c1_times_nu_20__3
      real :: c2_times_nu_20

      ! - 
   
      if ( sensor /= sensor_saved ) then
   
         select case ( sensor )
         case ('GOES-08')
            nu_20 = 2559.8724
            a1_20 = -0.63937510
            a2_20 = 1.0008597
         case ('GOES-09')
            nu_20 = 2556.0911
            a1_20 = -0.59821738
            a2_20 = 1.0007427
         case ('GOES-10')
            nu_20 = 2554.4746
            a1_20 = -0.63529295
            a2_20 = 1.0007974
         case ( 'GOES-11')
            nu_20 = 2562.0840
            a1_20 = -0.64474519
            a2_20 = 1.0007721
         case ( 'GOES-12')
            nu_20 = 2651.7708
            a1_20 = -1.9052739
            a2_20 = 1.003010 
         case ( 'GOES-13')
            nu_20 = 2563.5837
            a1_20 = -1.4757495
            a2_20 = 1.0021472
         case ( 'GOES-14')
            nu_20 = 2573.8799 
            a1_20 = -1.5595940 
            a2_20 =  1.0021890
         case ( 'GOES-15')
            nu_20 = 2562.1383
            a1_20 = -1.661627
            a2_20 = 1.0023207
         case ( 'MODIS-AQUA')
            nu_20 = 2642.3820
            a1_20 = -0.48854008
            a2_20 = 1.0005735
         case ( 'MODIS-TERRA')
            nu_20 = 2642.2448
	         a1_20 = -0.48752945
	         a2_20 = 1.0005738
         case ( 'VIIRS')
            nu_20 = 2708.3865
	         a1_20 = -0.59392036
	         a2_20 = 1.0006466	
         case ( 'NOAA-05','TIROS-N')
            nu_20 = 2655.7409
	         a1_20 = -1.6485446
	         a2_20 = 1.0020894      
         case ( 'NOAA-06')
            nu_20 = 2671.5433
	         a1_20 = -1.7667110
	         a2_20 = 1.0024428	      
         case ( 'NOAA-07')
            nu_20 = 2684.5233
	         a1_20 = -1.9488269
	         a2_20 = 1.0029260	      
         case ( 'NOAA-08')
            nu_20 = 2651.3776
	         a1_20 = -1.7764105
	         a2_20 = 1.0024260	  
         case ( 'NOAA-09')
            nu_20 = 2690.0451
	         a1_20 = -1.8832662
	         a2_20 = 1.0028978	      
         case ( 'NOAA-10')
            nu_20 = 2672.6164
	         a1_20 = -1.7986926
	         a2_20 = 1.0026426	      
         case ( 'NOAA-11')
            nu_20 = 2680.05
	         a1_20 = -1.738973
	         a2_20 = 1.003354 	  
         case ( 'NOAA-12')
            nu_20 = 2651.7708
	         a1_20 = -1.9052739
	         a2_20 = 1.003010	        
         case ( 'NOAA-14')
            nu_20 = 2654.25
	         a1_20 = -1.885330
	         a2_20 = 1.003839	
         case ( 'NOAA-15')
            nu_20 = 2695.9743
	         a1_20 = -1.624481
	         a2_20 = 1.001989		 
         case ( 'NOAA-16')
            nu_20 = 2681.2540
	         a1_20 = -1.6774586
	         a2_20 = 1.0017316 	 
         case ( 'NOAA-17')
            nu_20 = 2669.1414
	         a1_20 = -1.7002941
	         a2_20 = 1.0026724 	   
         case ( 'NOAA-18')
            nu_20 = 2660.6468
	         a1_20 = -1.7222650
	         a2_20 = 1.0028633	 
         case ( 'NOAA-19')
            nu_20 = 2670.2425
	         a1_20 = -1.6863857
	         a2_20 = 1.0025955		     
	      case ('MTSAT-1R')
	         nu_20 = 2647.9998
	         a1_20 = -2.455206
	         a2_20 = 1.0042972 
      case ('MTSAT-2')
	      nu_20 = 2680.1828
	      a1_20 = -2.3876343
	      a2_20 = 1.0042061   
      case ('METOP-A')
	      nu_20 = 2687.0392
	      a1_20 = -2.0653147
	      a2_20 = 1.0034418   
      case ('METOP-B')
	      nu_20 = 2664.3384
	      a1_20 = -1.7711318
	      a2_20 = 1.0029931   
      case ('Meteosat-8')
	      nu_20 = 2561.4547
	      a1_20 = -3.2692076
	      a2_20 = 1.0056489   
      case ('Meteosat-9')
	      nu_20 = 2562.2502
	      a1_20 = -3.2790754
	      a2_20 = 1.0059926   
      case ('Meteosat-10')
	      nu_20 = 2560.1576
	      a1_20 = -3.2146560
	      a2_20 = 1.0058230  
      case ('GOES-16','ABI')
	      nu_20 = 2562.1383 ! faked from goes-15
         a1_20 = -1.661627 ! faked from goes-15
         a2_20 = 1.0023207 ! faked from goes-15
	   case('COMS-1')
         nu_20 = 2675.0265
	      a1_20 = -2.2829416
	      a2_20 = 1.0037865   
    
	  
      case default
         print*,'missing sensor calibration for sensor ', sensor
	      print*, 'add to get_planck_radiance_39um'
            stop
      end select 	
	  
    
      c1_times_nu_20__3 =  c1 * nu_20 ** 3
      c2_times_nu_20 = c2 * nu_20
   
      do i = 1 , nplanck 
         t_planck(i) = T_planck_min + ( i - 1 ) * delta_T_planck
         B20(i) = c1_times_nu_20__3  / ( exp ( ( c2_times_nu_20 ) / &
              (( T_planck(i) - a1_20 ) / a2_20 ) ) - 1.0) 
      end do
      
      sensor_saved = sensor 
         
   end if
   
  
   
   
      		  

      l = ( tmp - T_planck_min ) / delta_T_planck
      l = max (1, min ( nplanck - 1 , l ) )
	  
      dB_dT_tmp = (B20(l+1)-B20(l))/(T_planck(l+1)-T_planck(l))
	  
      rad = B20(l) + (tmp - T_planck(l)) * (dB_dT_tmp)
  
  
   end function get_planck_radiance_39um




end module dcomp_science_tools_mod
