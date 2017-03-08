!  $Id$
!   this is for all basics scientific meteoerological etc routines
!
module cx_science_tools_mod
   
   
   use CX_CONSTANTS_MOD, only: &
      real4 , Missing_Value_Real4
      
   implicit none
   private
   real (kind=real4), parameter, public:: GLINT_ZEN_THRESH = 40.0
   real (kind=real4), parameter, public:: STEFAN_BOLTZMANN_CONSTANT = 5.670e-08  !W/m^2/K^4
   real (kind=real4), parameter, public:: SOLAR_CONSTANT = 1360.0   !W/m^2

  !--- default modis white sky values for Snow and water
   real, parameter, public:: REF_SFC_WHITE_SKY_WATER = 5.0 
   
   public :: VAPOR
   public :: VAPOR_ICE
   public :: WIND_SPEED
   public :: WIND_DIRECTION 

contains

   !---------------------------------------------------------------------
   ! SUBPROGRAM:  W3FC05        EARTH U,V WIND COMPONENTS TO DIR AND SPD
   !   PRGMMR: CHASE            ORG: NMC421      DATE:88-10-26
   !
   ! ABSTRACT: GIVEN THE TRUE (EARTH ORIENTED) WIND COMPONENTS
   !   COMPUTE THE WIND DIRECTION AND SPEED.
   !   INPUT WINDS AT THE POLE ARE ASSUMED TO FOLLOW THE WMO
   !   CONVENTIONS, WITH THE OUTPUT DIRECTION COMPUTED IN ACCORDANCE
   !   WITH WMO STANDARDS FOR REPORTING WINDS AT THE POLE.
   !   (SEE OFFICE NOTE 241 FOR WMO DEFINITION.)
   !
   ! PROGRAM HISTORY LOG:
   !   81-12-30  STACKPOLE, JOHN
   !   88-10-19  CHASE, P.   ALLOW OUTPUT VALUES TO OVERLAY INPUT
   !   89-01-21  R.E.JONES   CONVERT TO MICROSOFT FORTRAN 4.10
   !   90-06-11  R.E.JONES   CONVERT TO SUN FORTRAN 1.3
   !   91-03-30  R.E.JONES   SiliconGraphics FORTRAN
   !
   ! USAGE:    call W3FC05 (U, V, DIR, SPD)
   !
   !   INPUT ARGUMENT LIST:
   !     U        - real*4 EARTH-ORIENTED U-COMPONENT
   !     V        - real*4 EARTH-ORIENTED V-COMPONENT
   !
   !   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
   !     DIR      - real*4 WIND DIRECTION, DEGREES.  VALUES WILL
   !                BE FROM 0 TO 360 INCLUSIVE.
   !     SPD      - real*4 WIND SPEED IN SAME UNITS AS INPUT
   !---------------------------------------------------------------------

subroutine WIND_SPEED_AND_DIRECTION(u,v,dir,spd)
  real, intent(in) :: u
  real, intent(in) :: v
  real, intent(out) :: spd
  real, intent(out) :: dir
 

  real, parameter:: SPDTST = 1.0e-10
  real, parameter:: RTOD = 57.2957795
  real, parameter:: dchalf = 180.0

  spd = Missing_Value_Real4
  dir = Missing_Value_Real4
  if (u == Missing_Value_Real4 .or. v == Missing_Value_Real4) then
     return
  endif 
  spd = sqrt(u * u + v * v)
  if (spd < SPDTST) THEN
        dir = 0.0
  else 
        dir = atan2(u,v) * RTOD + DCHALF
  endif

  return

end subroutine WIND_SPEED_AND_DIRECTION
!-------------------------------------------------------------------------------------
! elemental funcions wind_speed and wind_direction taken from W3FC03 from UCAR
!
! W3FC05 header documentation follows:
!
! ABSTRACT: GIVEN THE TRUE (EARTH ORIENTED) WIND COMPONENTS
!   COMPUTE THE WIND DIRECTION AND SPEED.
!   INPUT WINDS AT THE POLE ARE ASSUMED TO FOLLOW THE WMO
!   CONVENTIONS, WITH THE OUTPUT DIRECTION COMPUTED IN ACCORDANCE
!   WITH WMO STANDARDS FOR REPORTING WINDS AT THE POLE.
!   (SEE OFFICE NOTE 241 FOR WMO DEFINITION.)
!
! PROGRAM HISTORY LOG:
!   81-12-30  STACKPOLE, JOHN
!   88-10-19  CHASE, P.   ALLOW OUTPUT VALUES TO OVERLAY INPUT
!   89-01-21  R.E.JONES   CONVERT TO MICROSOFT FORTRAN 4.10
!   90-06-11  R.E.JONES   CONVERT TO SUN FORTRAN 1.3
!   91-03-30  R.E.JONES   SiliconGraphics FORTRAN
!
! USAGE:    CALL W3FC05 (U, V, DIR, SPD)
!
!   INPUT ARGUMENT LIST:
!     U        - real*4 EARTH-ORIENTED U-COMPONENT
!     V        - real*4 EARTH-ORIENTED V-COMPONENT
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!     DIR      - real*4 WIND DIRECTION, DEGREES.  VALUES WILL
!                BE FROM 0 TO 360 INCLUSIVE.
!     SPD      - real*4 WIND SPEED IN SAME UNITS AS INPUT
!-------------------------------------------------------------------------------------
real elemental function WIND_SPEED ( u ,v )

  real, intent(in)::  u
  real, intent(in):: v

  if (u == Missing_Value_Real4 .or. v == Missing_Value_Real4) then
     wind_speed = Missing_Value_Real4
  else
     wind_speed = sqrt(u * u + v * v)
  endif

end function wind_speed
!-------------------------------------------------------------------------------------
! taken from W3FC03 from UCAR
!-------------------------------------------------------------------------------------
real elemental function WIND_DIRECTION ( u ,v )

  real, intent(in)::  u
  real, intent(in):: v
  real, parameter:: rtod = 57.2957795
  real, parameter:: dchalf = 180.0
  real, parameter:: spdtst = 1.0e-10

  if (u == Missing_Value_Real4 .or. v == Missing_Value_Real4) then
     wind_direction = Missing_Value_Real4
  else
     if (abs(u) < spdtst .and. abs(v) < spdtst) then
      wind_direction =  0.0
     else
      wind_direction = atan2(u,v) * rtod + dchalf
     endif
  endif

end function wind_direction

!----------------------------------------------------------------
! functions to compute some needed water vapor parameters
!----------------------------------------------------------------
 function VAPOR(T) result(es)
                                                                     
!  T in Kelvin                                                          
!  es in mbar

  implicit none
  real, intent (in) :: T
  real :: es

   es = 6.112 * exp(17.67 * (T-273.16) / (T - 29.66))

  return 
end function VAPOR

!---- saturation vapor pressure for ice
function VAPOR_ICE(T) result(es)
   implicit none
   real, intent(in):: T
   real:: es
     es = 6.1078 * exp(21.8745584 * (T-273.16) / (T - 7.66))
  return 
end function VAPOR_ICE

end module cx_science_tools_mod
