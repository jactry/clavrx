!$Id: cgms_nav.f90,v 1.4.2.2 2014/01/26 04:48:31 heidinger Exp $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: cgms_nav.f90 (src)
!       CGMS_NAV (program)
!
! PURPOSE: These routines do the navigation for all satellites with a fixed grid
!
! DESCRIPTION: Based on LRIT/HRIT Global Specification 
!              (CGMS 03, Issue 2.6, 12.08.1999)
!              This includes MTSAT and SEVIRI
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
!--------------------------------------------------------------------------------------
MODULE CGMS_NAV

  REAL(KIND(0.0d0)), PARAMETER, PRIVATE :: PI=3.14159265359d0

  REAL(KIND(0.0d0)), PARAMETER, PRIVATE ::  SAT_HEIGHT= 42164.0d0  ! distance from Earth centre to satellite    
  REAL(KIND(0.0d0)), PARAMETER, PRIVATE ::  R_EQ = 6378.169d0      ! radius from Earth centre to equator
  REAL(KIND(0.0d0)), PARAMETER, PRIVATE ::  R_POL= 6356.5838d0     !
  
  CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! subroutine pixcoord2geocoord_cgms
!  
!  PURPOSE: 
!  return the geograhic latitude and longitude of an Fixed image   
!    for a given pair of latitude/longitude.             
!    (based on the formulas given in Ref. [1])                
!                                                        
!                                                        
!  DEPENDENCIES:                                         
!    none                                                 
!                                                        
!                                                        
!  REFERENCE:                                            
!  [1] LRIT/HRIT Global Specification                     
!      (CGMS 03, Issue 2.6, 12.08.1999)                  
!                                                        
!  INPUT:                                                
!    row   (int) row-value of the pixel
!    colum (int) columb-value of the pixel
!    ccoff (int) coefficient of the scalling function    
!                      (see page 28, Ref [1])                                  
!    lloff (int) coefficient of the scalling function    
!                      (see page 28, Ref [1])                                  
!    lfac (int) coefficient of the scalling function    
!                      (see page 28, Ref [1])                                  
!    cfac (int) coefficient of the scalling function    
!                      (see page 28, Ref [1])                                  
!    deg_space (int) whether or not coefficents are in deg space
!    sub_lon (double) satellite subpoint in degree space
!                                                        
!  OUTPUT:                                               
!    latitude (double) geographic Latitude of the wanted pixel [Degrees]
!    longitude (double) geographic Longitude of the wanted pixel [Degrees]
!                                                        
!                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE pixcoord2geocoord_cgms( column,  row,  ccoff,  lloff, lfac, cfac, deg_space, sub_lon, latitude, longitude)

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: column, row, ccoff, lloff, lfac, cfac
  INTEGER, INTENT (IN) :: deg_space
  REAL(KIND(0.0d0)) , INTENT (IN) :: sub_lon
  REAL(KIND(0.0d0)) , INTENT (OUT) :: latitude, longitude

  
  REAL(KIND(0.0d0)) :: s1, s2, s3, sn, sd, sxy
  REAL(KIND(0.0d0)) :: x, y
  REAL(KIND(0.0d0)) :: longi, lati, sub_loni


  REAL(KIND(0.0d0)) :: sa

  INTEGER ::  c, l


  c=column
  l=row
  !Since we are going to radian space, convert sublon to radians
  sub_loni = sub_lon *PI / 180.0d0
  
  !  calculate viewing angle of the satellite by use of the equation 
  !  on page 28, Ref [1].
  
  x = (2.0d0**16.0d0 * ( DBLE(c) - DBLE(ccoff) )) / DBLE(CFAC)
  y = (2.0d0**16.0d0 * ( DBLE(l) - DBLE(lloff) )) / DBLE(LFAC)
  
  !MTSAT provides the linear coefficents in degree space, so
  !This simply 
  IF (deg_space == 1) then 
    x = x * (PI/180.0d0)
    y = y * (PI/180.0d0)
  ENDIF
  
  !  now calculate the inverse projection using equations on page 25, Ref. [1]  

  !  first check for visibility, whether the pixel is located on the Earth 
  !  surface or in space. 
  !  To do this calculate the argument to sqrt of "sd", which is named "sa". 
  !  If it is negative then the sqrt will return NaN and the pixel will be 
  !  located in space, otherwise all is fine and the pixel is located on the 
  !  Earth surface.

  sa =  (SAT_HEIGHT * cos(x) * cos(y) )**2 - (cos(y)*cos(y) + 1.006803d0 * sin(y)*sin(y)) * 1737121856.0d0

  ! take care if the pixel is in space, that an error code will be returned
  if ( sa .LE. 0.0 ) then
     latitude = -999.0
     longitude = -999.0
     return 
  end if
  
  ! now calculate the rest of the formulas using eq. on page 25 Ref [1]

  sd = sqrt( (SAT_HEIGHT * cos(x) * cos(y) )**2 - (cos(y)*cos(y) + 1.006803d0 * sin(y)*sin(y)) * 1737121856.0d0 )
  sn = (SAT_HEIGHT * cos(x) * cos(y) - sd) / ( cos(y)*cos(y) + 1.006803d0 * sin(y)*sin(y) ) 
  
  s1 = SAT_HEIGHT - sn * cos(x) * cos(y)
  s2 = sn * sin(x) * cos(y)
  s3 = -sn * sin(y)

  sxy = sqrt( s1*s1 + s2*s2 )

  ! using the previous calculations now the inverse projection can be
  ! calculated, which means calculating the lat./long. from the pixel
  ! row and column by equations on page 25, Ref [1].

  longi = atan(s2/s1) + sub_loni
  lati  = atan((1.006803d0*s3)/sxy)

  ! convert from radians into degrees
  latitude = (lati*(180.0d0/PI)) 
  longitude = longi*(180.0d0/PI)

END SUBROUTINE pixcoord2geocoord_cgms


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine geocoord2pixcoord_cgms                                    
!                                                       
! PURPOSE:                                              
!   return the pixel column and line of an fixed grid image 
!   for a given pair of geographic latitude/longitude.                   
!   (based on the formulas given in Ref. [1])                
!                                                       
!                                                       
! DEPENDENCIES:                                         
!   none                                       
!                                                       
!                                                       
! REFERENCE:                                            
! [1] LRIT/HRIT Global Specification                     
!     (CGMS 03, Issue 2.6, 12.08.1999)                  
!     for the parameters used in the program.
!                                                       
!                                                       
! MODIFICATION HISTORY:
!   Based on EUMETSAT code Version 1.01
!    Copyright(c) EUMETSAT 2005, 2009
!   Modified for all fixed grids by WCS3, 2011
!                                                       
!                                                       
! INPUT:                                                
!   latitude  (double) geographic Latitude of a point [Degrees] 
!   longitude (double) geographic Longitude of a point [Degrees]
!   ccoff (int)        coefficient of the scalling function   
!                      (see page 28, Ref [1])                                  
!   lloff (int)        coefficient of the scalling function   
!                      (see page 28, Ref [1])                 
!    lfac (int) coefficient of the scalling function    
!                      (see page 28, Ref [1])                                  
!    cfac (int) coefficient of the scalling function    
!                      (see page 28, Ref [1])                                  
!    deg_space (int) whether or not coefficents are in deg space
!    sub_lon (double) satellite subpoint in degree space
!                                                       
! OUTPUT:                                               
!   row    (int)       row-value of the pixel
!   column (int)       column-value of the pixel
!                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE geocoord2pixcoord_cgms( latitude,  longitude,  ccoff,  lloff, lfac, cfac, deg_space, sub_lon, column, row)

  IMPLICIT NONE


  REAL(KIND(0.0d0)), INTENT (IN) :: latitude, longitude
  REAL(KIND(0.0d0)), INTENT (IN) :: sub_lon
  INTEGER, INTENT (IN)  :: ccoff, lloff, lfac, cfac
  INTEGER, INTENT (IN)  :: deg_space
  INTEGER, INTENT (OUT) :: column, row
  

  INTEGER :: ccc=0, lll=0

  REAL(KIND(0.0d0)) :: lati, longi, subloni
  REAL(KIND(0.0d0)) :: c_lat
  REAL(KIND(0.0d0)) :: lat
  REAL(KIND(0.0d0)) :: lon
  REAL(KIND(0.0d0)) :: r1, r2, r3, rn, re, rl
  REAL(KIND(0.0d0)) :: xx, yy
  REAL(KIND(0.0d0)) :: cc, ll
  REAL(KIND(0.0d0)) :: dotprod


  lati= latitude
  longi= longitude
  
  ! Convert for international date line
  IF(longi .gt. 180.0d0) longi= longi - 360.0d0

  ! check if the values are sane, otherwise return error value
  if (lati .LT. -90.0d0 .OR. lati .GT. 90.0d0 .OR. longi .LT. -180.0d0 .OR. longi .GT. 180.0d0 ) then
     row = -999
     column = -999
     return
  end if

  ! convert them to radians 
  lat = lati*PI / 180.0d0
  lon = longi *PI / 180.0d0
  !Since we are going to radian space, convert sublon to radians
  subloni = sub_lon *PI / 180.0d0


  ! calculate the geocentric latitude from the       
  ! geographic one using equations on page 24, Ref. [1] 

  c_lat = atan ( (0.993243d0*(sin(lat)/cos(lat)) ))
      

  ! using c_lat calculate the length from the Earth 
  ! centre to the surface of the Earth ellipsoid    
  ! equations on page 23, Ref [1]                      
  
  re = R_POL / sqrt( (1.0d0 - 0.00675701d0 * cos(c_lat) * cos(c_lat) ) )


  ! calculate the forward projection using equations on page 24, Ref. [1]

  rl = re
  r1 = SAT_HEIGHT - rl * cos(c_lat) * cos(lon - subloni)
  r2 = - rl *  cos(c_lat) * sin(lon - subloni)
  r3 = rl * sin(c_lat)
  rn = sqrt( r1*r1 + r2*r2 +r3*r3 )

  ! check for visibility, whether the point on the Earth given by the
  ! latitude/longitude pair is visible from the satellte or not. This 
  ! is given by the dot product between the vectors of:
  ! 1) the point to the spacecraft,
  ! 2) the point to the centre of the Earth.
  ! If the dot product is positive the point is visible otherwise it
  ! is invisible.

  dotprod = r1*(rl * cos(c_lat) * cos(lon - subloni)) - r2*r2 - r3*r3*((r_EQ/R_POL)**2)

  if (dotprod .LE. 0.0d0 ) then
     column = -999
     row = -999
     return
  end if

  ! the forward projection is x and y 

  xx = atan( (-r2/r1) )
  yy = asin( (-r3/rn) )
  
  !If the coefficents are in degrees space (ex. MTSAT) then we have to convert

  IF (deg_space == 1) THEN
      xx = xx * (180.0d0/PI)
      yy = yy * (180.0d0/PI)
  ENDIF


  ! convert to pixel column and row using the scaling functions on 
  ! page 28, Ref. [1]. And finding nearest integer value for them.
  
  !convert to degree space back from radian space, because coeffs are
  !in degree space

  cc = DBLE(ccoff) + xx *  2.0d0**(-16.0d0) * DBLE(cfac)
  ll = DBLE(lloff) + yy *  2.0d0**(-16.0d0) * DBLE(lfac)


  ccc=nint(cc)
  lll=nint(ll)		

  column = ccc
  row = lll
     

END SUBROUTINE geocoord2pixcoord_cgms


END MODULE CGMS_NAV
