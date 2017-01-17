!$Id$

function sensorname_from_wmoid (id) result(sensor_name)

integer, intent(in) :: id
character(len=20) :: sensor_name

   !--- select sensor name based on wmo id number
   select case (id)
      case(3)
         sensor_name = 'METOP-B'
      case(4)
         sensor_name = 'METOP-A'	 
      case(55)
         sensor_name = 'Meteosat-8' 
      case(56)
         sensor_name = 'Meteosat-9' 
      case(57)
         sensor_name = 'Meteosat-10' 
      case(70)
         sensor_name = 'Meteosat-11' 
      case(171)
         sensor_name = 'MTSAT-1R'
      case (172)
         sensor_name = 'MTSAT-2'
      case (173)
         sensor_name = 'AHI'     
      case (200:204)
         write(sensor_name, "('NOAA-',i2.2)") id - 192 
      case (205:209)
         write(sensor_name, "('NOAA-',i2.2)") id - 191
      case(223)
         sensor_name = 'NOAA-19'	  
      case(224)
         sensor_name = 'VIIRS'
      case(252)
         sensor_name = 'GOES-08'
      case(253)
         sensor_name = 'GOES-09'
      case(254)
         sensor_name = 'GOES-10'
      case(255)
         sensor_name = 'GOES-11'
      case(256)
         sensor_name = 'GOES-12'
      case(257)
         sensor_name = 'GOES-13'
      case(258)
         sensor_name = 'GOES-14'
      case(259)
         sensor_name = 'GOES-15'
      case(270)
         sensor_name = 'GOES-16'   
      case(705)
         sensor_name = 'NOAA-05'   
      case(706)
         sensor_name = 'NOAA-06'   
      case(707)
         sensor_name = 'NOAA-07'
      case(708)
         sensor_name = 'TIROS-N'   
      case(783)
         sensor_name = 'MODIS-TERRA'
      case(784)
         sensor_name = 'MODIS-AQUA'
      case(810)
         sensor_name = 'COMS-1'   
      case default
         print*,'please inform  andi.walther@ssec.wisc.edu wmo id: ', id 
         stop
   end select
   
   

end function sensorname_from_wmoid
