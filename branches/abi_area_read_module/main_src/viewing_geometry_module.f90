! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: viewing_geometry_module.f90 (src)
!       VIEWING_GEOMETRY_MODULE (program)
!
! PURPOSE: provide functions for viewing angle computations
!
! DESCRIPTION: 
!              History:  This module was created to house the routines used to 
!                        compute the angles that define the sensor and solar 
!                        viewing geomety.
!             
!                        The first version is simply a collection of prexisting
!                        routines from other modules.  Future versions will move
!                        towards more robust, efficient and elegant code
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! Angle Conventions used in CLAVR-x
!
! Sensor Zenith Angle:
!       Defined as zero at nadir and positive in both direction away from
!       nadir.  This is not the scan angle.  It accounts for the earth's
!       curvature. Ranges from 0 to 90 degrees.
!
! Solar Zenith Angle:
!       Defined so zero degrees means the sun is straightaove the pixel.
!       Daytime ranges from 0 to 90.0  Full Range is 0 to 180.0
!
! Solar Azimuth Angle:
!       Represents the angle from the pixel to location on the earth where
!       the solar zenith angle is zero (sub-sun point).
!       If the sub-sun point is due south, the value is 180.0 or -180.0
!       If the sub-sun point is due north, the value is 0.0 
!       If the sub-sun point is due west, the value is -90.0 
!       If the sub-sun point is due east, the value is 90.0 
!       The range is from -180 to 180 degrees
!       This definition is the same as employed by MODIS.   
!
! Sensor Azimuth Angle:
!       Represents the angle from the pixel to location on the earth where
!       sensor zenith angle is 0 (the sub-sat point) 
!       If the sub-sat point is due south, the value is 180.0 or -180.0
!       If the sub-sat point is due north, the value is 0.0 
!       If the sub-sat point is due west, the value is -90.0 
!       If the sub-sat point is due east, the value is 90.0 
!       The range is from -180 to 180 degrees
!       This definition is the same as employed by MODIS.   
!
! Relative Azimuth Angle:
!       Represents the relative difference in solar and sensor azimuth
!       angles.  Defined so that values less than 90 degrees are for
!       pixels that are between the sub-sat and the sun.  Pixels with
!       values greater than 90 are the opposite side of the sun.
!       Glint always occours with values less than 90.0. This
!       definition has always been used in CLAVR-x
!       Range is 0 to 180.0
!
! Scattering Angle:
!       The angle between the sun, the pixel and sensor.  A value of
!       180.0 would imply pure backscatter.  A value of zero would
!       imply pure forward scatter (not seen in satellite remote sensing).
!       Values less than 90 imply forward scattering.
!       Range is 0 to 180.0
!
! Glint Zenith Angle:
!       The angle defined between the ray from the sensor to the pixel
!       and the specular reflection angle from the pixel to the sun.
!       A value of zero mean your are viewing the sun image if the
!       surface were a mirror. Values less than 40 imply glint
!       impacted observations over ocean.  Range 0 to 180.0  
!--------------------------------------------------------------------------------------
 module VIEWING_GEOMETRY_MODULE

  implicit none
  
  public:: sensor_zenith, &
           great_circle_angle, & 
           sensor_zenith_avhrr_anchor, & 
           relative_azimuth_avhrr_anchor, & 
           glint_angle, &
           relative_azimuth, &
           scattering_angle, &
           sensor_azimuth, &
           POSSOL, &
           COMPUTE_SENSOR_ZENITH_GEO    !not needed
           
  real, parameter, private :: R_EARTH = 6378.2064   !km
  real, parameter, private :: PI = 3.141592653589793
  real, parameter, private :: DTOR = PI / 180.0

  contains
!-------------------------------------------------
! compute greate circle angular distance
! this is used in AVHRR angle routines
!-------------------------------------------------
   real function great_circle_angle(satlon,satlat,pixlon,pixlat)

   real, intent(in):: satlon
   real, intent(in):: satlat
   real, intent(in):: pixlon
   real, intent(in):: pixlat
   real:: cosgeo

   cosgeo = cos(pixlat*DTOR)*cos(satlat*DTOR)*cos((pixlon-satlon)*DTOR) + &
            sin(pixlat*DTOR)*sin(satlat*DTOR)
   cosgeo = max(-1.0,min(1.0,cosgeo))
   great_circle_angle = acos(cosgeo) / DTOR

  end function great_circle_angle
!-------------------------------------------------
! sensor_zenith_avhrr_anchor
!
! input:
! anchor_index = index of anchor point (1-51)
! geox = great circle distance angle in degrees
! Reference: CLAVR-1
!
! Note, application of the senzor_zenith function
! to AVHRR does not work outside of tropics.
! it would be nice to have a generic routine
! 
!
!--------------------------------------------------
   real function sensor_zenith_avhrr_anchor(geox,anchor_index)

   real, intent(in):: geox
   integer, intent(in):: anchor_index
   integer, parameter:: subsat_anchor_index = 26
   real, parameter:: delta_scan_angle = 2.16512
   sensor_zenith_avhrr_anchor = geox +  &
                                abs(subsat_anchor_index - anchor_index)*delta_scan_angle

  end function sensor_zenith_avhrr_anchor
!-------------------------------------------------
! relative_azimuth_avhrr_anchor
!
! input:
! geox = great circle distance angle in degrees
! solzen_subsat = solar zenith angle at sub-sat point
! solzen_pix = solar zenith angle at pixel 
!
! Reference: CLAVR-1
!
!---------------------------------------------------
  real function relative_azimuth_avhrr_anchor(geox,solzen_subsat,solzen_pix)

   real, intent(in):: geox
   real, intent(in):: solzen_subsat
   real, intent(in):: solzen_pix
   real:: cosgeo
   real:: cossolzen_pix
   real:: numor
   real:: denom
   real:: psix

   cosgeo = cos(geox*DTOR)
   cosgeo = max(-1.0,min(1.0,cosgeo))
   cossolzen_pix = cos(solzen_pix*DTOR)
   numor= cosgeo*cossolzen_pix - cos(solzen_subsat*DTOR)
   denom = sin(geox*DTOR)* sqrt(1.0-cossolzen_pix**2)

   if (denom == 0.0) then
         psix = 1.0
     else
         psix = numor/denom
   endif

   psix = max(-1.0,min(1.0,psix))
   relative_azimuth_avhrr_anchor = acos(psix)/DTOR

  end function relative_azimuth_avhrr_anchor
!-------------------------------------------------
!
! input:
! h = height above earth of satellite orbit (km)
! satlon = subsatellite longitude
! satlat = subsatellite latitude
! pixlon = satellite pixel longitude
! pixlat = satellite pixel latitude
!
!--------------------------------------------------
   real function sensor_zenith(h,satlon,satlat,pixlon,pixlat)

   real, intent(in):: h
   real, intent(in):: satlon
   real, intent(in):: satlat
   real, intent(in):: pixlon
   real, intent(in):: pixlat
   real:: r
   real:: xlon
   real:: xlat
   real:: beta

   xlon = (pixlon - satlon)*DTOR
   xlat = (pixlat - satlat)*DTOR
   beta = acos( cos(xlat) * cos(xlon) )
   r = R_EARTH
   sensor_zenith = (r+h) * sin(beta) / sqrt( (r*r) + (r+h)*(r+h) - 2.0*r*(r+h)*cos(beta) )
   sensor_zenith = max(-1.0, min(1.0,sensor_zenith))
   sensor_zenith = asin(sensor_zenith) / DTOR

  end function sensor_zenith
  !-------------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------------
  real elemental function relative_azimuth ( sol_az ,sen_az )

  real, intent(in)::  sol_az
  real, intent(in):: sen_az

  relative_azimuth = abs(sol_az - sen_az)
  if (relative_azimuth > 180.0) then
       relative_azimuth = 360.0 - relative_azimuth
  endif
  relative_azimuth = 180.0 - relative_azimuth

  end function relative_azimuth
  !------------------------------------------------------------------------------------
  ! Glint angle  (the angle difference between direct "specular" reflection off
  ! the surface and actual reflection toward the satellite.) 
  !------------------------------------------------------------------------------------
  real elemental function glint_angle ( sol_zen , sen_zen , rel_az  )
     real, intent(in) :: sol_zen
     real, intent(in) :: sen_zen
     real, intent(in) :: rel_az

     glint_angle = cos ( sol_zen * DTOR ) * cos ( sen_zen * DTOR ) +    &
                   sin ( sol_zen * DTOR ) * sin ( sen_zen * DTOR ) * cos ( rel_az * DTOR )

     glint_angle = max(-1.0 , min( glint_angle ,1.0 ) )

     glint_angle = acos(glint_angle) / DTOR

  end function glint_angle
  !-------------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------------
  real elemental function scattering_angle(sol_zen, sen_zen, rel_az)

   real, intent(in):: sol_zen
   real, intent(in):: sen_zen
   real, intent(in):: rel_az

   scattering_angle = cos(sol_zen*DTOR)*cos(sen_zen*DTOR) -    &
                      sin(sol_zen*DTOR)*sin(sen_zen*DTOR)*cos(rel_az*DTOR)

   scattering_angle = -1.0*scattering_angle

   scattering_angle = max(-1.0,min(scattering_angle,1.0))

   scattering_angle = acos(scattering_angle) / DTOR

  end function scattering_angle
  !-------------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------------
  real elemental function sensor_azimuth(satlon,satlat,pixlon,pixlat) 

   real, intent(in):: satlon
   real, intent(in):: satlat
   real, intent(in):: pixlon
   real, intent(in):: pixlat

   real:: xlon
   real:: xlat
   real:: beta
   real:: sine_beta

   xlon = (pixlon - satlon)*DTOR 
   xlat = (pixlat - satlat)*DTOR 

   beta = acos( cos(xlat) * cos(xlon) )
   sine_beta = sin(beta)
   if (abs(sine_beta) > epsilon(sine_beta)) then 
     sensor_azimuth = sin(xlon) / sine_beta
     sensor_azimuth = min(1.0, max(-1.0,sensor_azimuth))
     sensor_azimuth = asin(sensor_azimuth) / DTOR
   else
     sensor_azimuth = 0.0
   endif
   if (xlat < 0.0) then
       sensor_azimuth = 180.0 - sensor_azimuth
   endif

   if (sensor_azimuth < 0.0) then
       sensor_azimuth = sensor_azimuth + 360.0
   endif

   !--- MODIS Convection (same as solar azimuth)
   sensor_azimuth = sensor_azimuth - 180.0

  end function sensor_azimuth
!---------------------------------------------------------------------
! subroutine POSSOL(jday,tu,xlon,xlat,asol,phis)
! Description: compute solar angles
!
! input:
!         jday - julian day
!         tu -  time of day - fractional hours
!         lon - latitude in degrees
!         lat - latitude in degrees
!      
! output:
!         asol - solar zenith angle 
!         phi - solar azimuth angle
!
!
! This routine comes from the 6S RTM http://6s.ltdri.org/
!----------------------------------------------------------------------
   subroutine POSSOL(jday,tu,xlon,xlat,asol,phis)
   implicit none
   integer, intent(in):: jday 
   real, intent(in):: tu,xlon,xlat
   real, intent(out):: asol, phis
   real:: tsm, xlo, xla, xj, a1, a2, et, tsv, ah, a3, delta, amuzero, elev, &
         az, caz, azim, pi2
! real, parameter:: pi = 3.14159265,DTOR = pi / 180.0
!
! solar position (zenithal angle asol, azimuthal angle phi
! in degrees
! jday is the number of the day in the month
!
!
! mean solar time
!
      tsm = tu + xlon/15.0
      xlo = xlon*DTOR
      xla = xlat*DTOR
      xj = real(jday)
!
! time equation (mn.dec)
!
      a1 = (1.00554*xj-6.28306)*DTOR
      a2 = (1.93946*xj+23.35089)*DTOR
      et = -7.67825*sin(a1) - 10.09176*sin(a2)
!
! true solar time
!
      tsv = tsm + et/60.0
      tsv = (tsv-12.0)
!
! hour angle
!
      ah = tsv*15.0*DTOR
!
! solar declination (in radian)
!
      a3 = (0.9683*xj-78.00878)*DTOR
      delta = 23.4856*sin(a3)*DTOR
!
! elevation, azimuth
!
      amuzero = sin(xla)*sin(delta)+cos(xla)*cos(delta)*cos(ah)
      elev=asin(amuzero)
      az = cos(delta)*sin(ah)/cos(elev)
      caz=(-cos(xla)*sin(delta)+sin(xla)*cos(delta)*cos(ah))/cos(elev)
!
      if (az >= 1.0) then
        azim = asin(1.0)
      elseif (az <= -1.0) then
        azim = asin(-1.0)
      else
        azim = asin(az)
      endif
!
      if (caz <= 0.0) then
        azim = pi - azim
      endif

      if ((caz > 0.0) .and. (az <= 0.0)) then
          azim = 2 * pi + azim
      endif
      azim = azim + pi
      pi2 = 2 * pi
      if (azim> pi2) then
          azim = azim - pi2
      endif
!
! conversion in degrees
!
      
     
      elev = elev / DTOR
      asol = 90.0 - elev
      phis = azim / DTOR


!
! conversion to MODIS definition
!

      if (phis > 180.0) phis = phis - 360.0

 end subroutine POSSOL
!-------------------------------------------------
!
!--------------------------------------------------
function COMPUTE_SENSOR_ZENITH_GEO(satlon,satlat,pixlon,pixlat) result(sen_zen)

   real, intent(in)::  satlon
   real, intent(in):: satlat
   real, intent(in)::  pixlon
   real, intent(in):: pixlat
   real:: sen_zen
   real:: xlon
   real:: xlat
   real:: beta

   xlon = (pixlon - satlon)*DTOR 
   xlat = (pixlat - satlat)*DTOR 

   beta = acos( cos(xlat) * cos(xlon) )

!  sen_zen = asin(max(-1.0_real8, min(1.0_real8, &
!            42164.0* sin(beta)/ sqrt(1.808e09 - 5.3725e08*cos(beta)))))

   sen_zen = asin(max(-1.0, min(1.0,42164.0* sin(beta)/ sqrt(1.808e09 - 5.3725e08*cos(beta)))))

   sen_zen = sen_zen / DTOR

end function COMPUTE_SENSOR_ZENITH_GEO

   !------------------------------------------------------------------
   ! input
   !   lon0 : scalar of reference longitude value in degree
   !   lat0 : scalar of reference latitude value in degree
   !   lon:    longitude value in degree ( any dimension)
   !   lat :   latitude value in degree ( any dimension)
   !
   ! output:
   !    distance on Earth in km
   !     
   !  assumption is a pure spherical Earth
   !  reference: http://www.movable-type.co.uk/scripts/gis-faq-5.1.html
   !
   !------------------------------------------------------------------
   elemental real function great_circle_distance ( lon0 ,lat0, lon ,lat )
      real, intent(in) :: lon0 , lat0
      real, intent(in) :: lat , lon
      
      real :: dlon , dlat
      real, parameter :: radius_earth = 6367.
      real :: a , c
      
      real :: lon0_rad, lat0_rad
      real :: lon_rad , lat_rad
      
      lon0_rad = lon0 * DTOR
      lon_rad = lon * DTOR
      lat0_rad = lat0 * DTOR
      lat_rad = lat * DTOR
      
      dlon = lon0_rad - lon_rad
      dlat = lat0_rad - lat_rad
      
      a = ( sin ( dlat / 2)) ** 2 + cos( lat_rad ) * (sin (dlon/2)) ** 2  
      c = 2 * asin (  sqrt(a))
      great_circle_distance = radius_earth * c
  
  
   end function great_circle_distance
  

end module VIEWING_GEOMETRY_MODULE


! this should go in a test dir
!program test_it
!use viewing_geometry_module

!print*,great_circle_distance( 1.,0.,[0.,1.],[0.,2])
!end program
