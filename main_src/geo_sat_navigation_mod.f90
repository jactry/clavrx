! $Id$

module geo_sat_navigation_mod

   real ( kind = 8 ) , parameter :: PI = 4. * atan (1.)
   real ( kind = 8 ) , parameter :: DEG_TO_RAD = PI/180.
   real ( kind = 8 ) , parameter :: RAD_TO_DEG = 180./PI
   
   real ( kind = 8 ) , parameter :: R_POL = 6356.7523
   real ( kind = 8 ) , parameter :: R_EQ = 6378.1370
   real ( kind = 8 ) , parameter :: FLAT = 1.0/298.257222101
   real ( kind = 8 ) , parameter :: FLAT_P = 1/((1.- FLAT) * ( 1.- FLAT))
   
   real ( kind = 8 ) , parameter :: H_MSG = 42164.0
   real ( kind = 8 ) , parameter :: H_GOESR = 42164.16
   
   real ( kind = 8 ) , parameter :: MISSING_D = -999.000
   
contains
   
   subroutine fgf_to_earth ( xx, yy , cfac, coff, lfac, loff, sub_lon &
         , lon , lat )
         
      real(8) :: xx , yy
      real(8) :: cfac
      real(8) :: coff
      real(8) :: lfac
      real(8):: loff
      real(8) :: sub_lon   
      real(8) :: lat, lon
      real(8) :: lambda, theta
     
      
      call jmaidx_to_sat ( xx, yy, cfac, coff, lfac, loff, lambda , theta )
      
      
      
      call sat_to_earth ( lambda, theta, sub_lon, lon,lat )
      
   end subroutine fgf_to_earth 
   
   !
   !
   !
   subroutine earth_to_fgf ( lon, lat , cfac, coff, lfac, loff, sub_lon &
         , xx, yy )
         
      real(8) :: xx , yy
      real(8) :: cfac
      real(8) :: coff
      real(8) :: lfac
      real(8):: loff
      real(8) :: sub_lon   
      real(8) :: lat, lon
      real(8) :: lambda, theta
     
      
      call earth_to_sat(lon,lat , sub_lon, lambda, theta )
      
    
      
      call sat_to_jmaidx (lambda , theta , cfac, coff, lfac, loff,xx,yy) 
     
     
   end subroutine earth_to_fgf 
    
   !
   !
   subroutine jmaidx_to_sat ( xx, yy, cfac, coff, lfac, loff, lambda , theta )
      real(8) :: xx , yy
      real(8) :: cfac
      real(8) :: coff
      real(8) :: lfac
      real(8) :: loff
      real(8) :: lambda, theta
   
      real (8) , parameter :: FRNT_FACTOR =  2.0**16
      
      lambda = ((FRNT_FACTOR * ( xx- coff ))/ cfac ) * DEG_TO_RAD
      theta = ((FRNT_FACTOR * ( yy- loff ))/ lfac ) * DEG_TO_RAD 
      
     
   
   end subroutine jmaidx_to_sat
   
   !
   !
   !
   subroutine sat_to_jmaidx ( lambda, theta , cfac, coff, lfac, loff, xx , yy )
      real(8) :: xx , yy
      real(8) :: cfac
      real(8) :: coff
      real(8) :: lfac
      real(8) :: loff
      real(8), intent(in) :: lambda, theta
      
      real (8) , parameter :: FRNT_FACTOR =  2.0**(-16)
     
     
      xx = coff + ( lambda * FRNT_FACTOR * cfac * RAD_TO_DEG )
      yy = loff + ( theta * FRNT_FACTOR * lfac * RAD_TO_DEG )
      
      
   end subroutine sat_to_jmaidx
   
   !
   !
   !
   subroutine sat_to_earth ( lambda, theta, sub_lon, lon,lat )
      implicit none
      real(8) :: lambda, theta
      real(8) :: sub_lon 
      real(8) :: lat, lon
      
      real(8) :: sub_lon_radians
      real(8) :: x 
      real(8) :: y 
      real (8) :: h
      real(8) :: d, c1, c2, s_d, s_n,s_1,s_2,s_3,s_xy,geographic_lon,geographic_lat 
      
      sub_lon_radians = sub_lon * DEG_TO_RAD
      
     
      x = lambda
      y = theta
      h = H_MSG
      ! call geos_to_geos ( lambda, theta) 
      
      d = ( H * H - R_EQ * R_EQ)
      c1 = (h * cos(x) * cos(y)) * (h * cos(x) * cos(y))
      c2 = (cos(y) * cos(y) + FLAT_P * sin(y) * sin(y)) * d
      
      if ( c1 < c2 ) then
         lon = -999.00000
         lat = -999.0000
         return
      end if
      
      s_d = SQRT ( c1-c2)
      s_n = (h * cos(x) * cos(y) - s_d) / (cos(y) * cos(y) + FLAT_P * sin(y) * sin(y))
      s_1 = h -s_n * cos(x) * cos(y)
      s_2 = s_n * sin(x) * cos(y)
      s_3 = (-1) * s_n * sin(y)
      
      s_xy = SQRT (s_1 * s_1 + s_2 * s_2 )
      
      
      
      geographic_lon = atan ( s_2/s_1) + sub_lon_radians
      geographic_lat = atan(FLAT_P*(s_3/s_xy))
      
       
      lon= RAD_TO_DEG * geographic_lon
      lat= RAD_TO_DEG * geographic_lat
      
   
   end subroutine sat_to_earth
   
   subroutine earth_to_sat( lon,lat,sub_lon,lambda, theta )
      implicit none
      real(8) , intent(out) :: lambda, theta
      real(8) :: sub_lon 
      real(8) :: lat, lon
      
      real(8) :: sub_lon_radians
      real(8) :: x 
      real(8) :: y 
      real (8) :: h
      real(8) :: geocentric_lat 
      real(8) :: d, c1, c2, s_d, s_n,r_1,r_2,r_3,s_xy
      real (8) ::geographic_lon,geographic_lat 
      real (8) :: r_earth
      
      
      
      h = H_MSG
      
      geographic_lat = lat * DEG_TO_RAD
      geographic_lon = lon * DEG_TO_RAD
      
      sub_lon_radians = sub_lon * DEG_TO_RAD
      
      geocentric_lat = ATAN ((( R_POL ** 2)/(R_EQ **2 )) * TAN (geographic_lat))
      r_earth = R_POL / SQRT ( 1.0 - (( R_EQ ** 2 - R_POL ** 2) / (R_EQ * R_EQ )) &
         *  COS (geographic_lat) * COS(geographic_lat)   ) 
         
      r_1 = h - r_earth*cos(geocentric_lat)*cos(geographic_lon - sub_lon_radians)
      r_2 = -r_earth*cos(geocentric_lat)*sin(geographic_lon - sub_lon_radians)
      r_3 =r_earth*sin(geocentric_lat)
     
      lambda = -999.
      theta  = -999.
      
   
      lambda = atan (-r_2/r_1)
       theta = asin(-r_3/sqrt(r_1*r_1 + r_2*r_2 + r_3*r_3))
     
   end subroutine earth_to_sat
   
   
   
   subroutine geos_to_geos ( lambda, theta )
      real(8) :: lambda, theta
      real(8) :: lambda_loc, theta_loc
      
      lambda_loc = atan( tan(lambda)/cos(theta) )
      theta_loc = asin( sin(theta)*cos(lambda) )
      
      lambda = lambda_loc
      theta = theta_loc
      
   end subroutine 
   
   
end module geo_sat_navigation_mod
