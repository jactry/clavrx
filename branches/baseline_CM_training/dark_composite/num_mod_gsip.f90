!------------------------------------------------------------------
! This module stores numerical routines needed in GSIP-fd
!
!
! CVS info:
!
!  $Id: num_mod_gsip.f90,v 1.1.1.1 2012/05/18 20:18:40 heidinger Exp $
!
!  $Log: num_mod_gsip.f90,v $
!  Revision 1.1.1.1  2012/05/18 20:18:40  heidinger
!
!
!  Revision 1.1  2010/08/12 00:58:18  heidinger
!  new
!
!  Revision 1.7  2005/09/09 21:32:27  heidinger
!  thors working version wth pCRTM
!
!  Revision 1.4  2004/12/03 23:13:31  heidinger
!  added diag routines
!
!  Revision 1.2  2004/10/28 02:11:51  heidinger
!  added goes-9 capability
!
!  Revision 1.1  2004/08/23 16:09:02  heidinger
!  no message
!
!
!-------------------------------------------------------------------
module NUMERICAL_ROUTINES

  use CONSTANTS

  implicit none

  private
  public :: LOCATE, POSSOL, JULIAN, COMPUTE_SATELLITE_ANGLES, &
            ICNVRT,COMPUTE_MONTH, COMPUTE_DAY, &
            VAPOR, SATVAPLIQ_OLD,SATVAPLIQ, SATVAPLIQB,SATVAPICE_OLD,SATVAPICE

contains

!------------------------------------------------------------------------
! Define calibration parameters
!
!  sat_id = satellite identification
!           8:  GOES-8
!           9:  GOES-9
!          10: GOES-10
!          11: GOES-11
!          12: GOES-12
!
! source of calibration - Mike Weinreb's GOES Calibration Webpage
!                         http://www.oso.noaa.gov/goes/goes-calibration/
!-------------------------------------------------------------------------
  subroutine DEFINE_GSIP_CAL(sat_id, Fo_2, year_launch,jday_launch, ch1_slope_degrad, &
                             ch1_slope, a, b, nu, slope, intercept)

!   Arguments
    integer,                           intent(in)  :: sat_id
    integer,                           intent(out) :: jday_launch, year_launch
    real (kind=real4),                 intent(out) :: Fo_2, ch1_slope, ch1_slope_degrad
    real (kind=real4), dimension(:),   intent(out) :: slope, intercept
    real (kind=real4), dimension(:,:), intent(out) :: nu, a, b

!   Local variables
    integer :: nchan

    nchan = size(slope)

    slope = reshape((/0.0, 227.3889, 38.8383, 5.2285, 5.0273, 5.5297/), (/nchan/))
    intercept = reshape((/29.0, 68.2167, 29.1287, 15.6854, 15.3332, 16.5892/), (/nchan/))
 
    if (sat_id == 8) then
       Fo_2 = 14.56                                              ! mW/m^2/cm^{-1}
       ch1_slope_degrad = 365.25 * 0.0001688
       ch1_slope = 0.1264
       year_launch = 1994
       jday_launch = 103
       nu(:,1) = reshape((/0.0, 2556.71, 1481.91, 934.30, 837.06, 751.91 /), (/nchan/))
       nu(:,2) = reshape((/0.0, 2558.62, 1481.91, 935.38, 837.00, 751.91 /), (/nchan/))
       a(:,1) = reshape((/0.0, -0.578526, -0.593903, -0.322585, -0.422571, -0.253449 /), (/nchan/))
       a(:,2) = reshape((/0.0, -0.581853, -0.593903, -0.351889, -0.466954, -0.253449 /), (/nchan/))
       b(:,1) = reshape((/0.0, 1.001512, 1.001418, 1.001271, 1.001170, 1.000743 /), (/nchan/))
       b(:,2) = reshape((/0.0, 1.001532, 1.001418, 1.001293, 1.001257, 1.000743 /), (/nchan/))
    elseif (sat_id == 9) then
       Fo_2 = 14.52                                              ! mW/m^2/cm^{-1} 
       ch1_slope_degrad = 365.25 * 0.000                         ! UNCALIBRATED VIS
       ch1_slope = 0.1365    !MODIS TERRA 2004 COMP
       year_launch = 1995
       jday_launch = 143
       nu(:,1) = reshape((/0.0, 2555.18, 1481.82, 934.59, 834.02, 751.91 /), (/nchan/))
       nu(:,2) = reshape((/0.0, 2555.18, 1481.91, 934.28, 834.09, 751.91 /), (/nchan/))
       a(:,1) = reshape((/0.0, -0.579908, -0.493016, -0.384798, -0.302995, -0.253449 /), (/nchan/))
       a(:,2) = reshape((/0.0, -0.579908, -0.493016, -0.384798, -0.302995, -0.253449 /), (/nchan/))
       b(:,1) = reshape((/0.0, 1.000942, 1.001076, 1.001293, 1.000941, 1.000743 /), (/nchan/))
       b(:,2) = reshape((/0.0, 1.000942, 1.001076, 1.001272, 1.000948, 1.000743 /), (/nchan/))
    elseif (sat_id == 10) then
       Fo_2 = 14.50                                              ! mW/m^2/cm^{-1} 
       ch1_slope_degrad = 365.25 * 0.0001022
       ch1_slope = 0.1165
       year_launch = 1997
       jday_launch = 115
       nu(:,1) = reshape((/0.0, 2552.9845, 1486.2212, 936.10260, 830.88473, 751.91 /), (/nchan/))
       nu(:,2) = reshape((/0.0, 2552.9845, 1486.2212, 935.98981, 830.89691, 751.91 /), (/nchan/))
       a(:,1) = reshape((/0.0, -0.60584483, -0.61653805, -0.27128884, -0.26505411, -0.253449 /), (/nchan/))
       a(:,2) = reshape((/0.0, -0.60584483, -0.61653805, -0.27064036, -0.26056452, -0.253449 /), (/nchan/))
       b(:,1) = reshape((/0.0, 1.0011017, 1.0014011, 1.0009674, 1.0009087, 1.000743 /), (/nchan/))
       b(:,2) = reshape((/0.0, 1.0011017, 1.0014011, 1.0009687, 1.0008962, 1.000743 /), (/nchan/))
    elseif (sat_id == 11) then
       Fo_2 = 14.58                                            ! mW/m^2/cm^{-1} --- CHECK THIS
       ch1_slope_degrad = 365.25 * 0.000                       ! UNCALIBRATED VIS
       ch1_slope = 0.1165
       year_launch = 2000
       jday_launch = 123
       nu(:,1) = reshape((/0.0, 2562.07, 1481.53, 931.76, 833.67, 751.91 /), (/nchan/))
       nu(:,2) = reshape((/0.0, 2562.07, 1481.53, 931.76, 833.04, 751.91 /), (/nchan/))
       a(:,1) = reshape((/0.0, -0.644790, -0.543401, -0.306809, -0.333216, -0.253449 /), (/nchan/))
       a(:,2) = reshape((/0.0, -0.644790, -0.543401, -0.306809, -0.315110, -0.253449 /), (/nchan/))
       b(:,1) = reshape((/0.0, 1.000775, 1.001495, 1.001274, 1.001000, 1.000743 /), (/nchan/))
       b(:,2) = reshape((/0.0, 1.000775, 1.001495, 1.001274, 1.000967, 1.000743 /), (/nchan/))
    elseif (sat_id == 12) then
       Fo_2 = 14.61                                            ! mW/m^2/cm^{-1} --- CHECK THIS
       ch1_slope_degrad = 365.25 * 0.000                       ! UNCALIBRATED VIS 
       ch1_slope = 0.1165
       year_launch = 2001
       jday_launch = 204
       nu(:,1) = reshape((/0.0, 2562.45, 1536.43, 933.21, 837.06, 751.91 /), (/nchan/))
       nu(:,2) = reshape((/0.0, 2562.45, 1536.94, 933.21, 837.00, 751.91 /), (/nchan/))
       a(:,1) = reshape((/0.0, -0.650731, -0.4764728, -0.360331, -0.422571, -0.253449 /), (/nchan/))
       a(:,2) = reshape((/0.0, -0.650731, -0.4775517, -0.360331, -0.466954, -0.253449 /), (/nchan/))
       b(:,1) = reshape((/0.0, 1.001520, 1.012420, 1.001306, 1.001170, 1.000743 /), (/nchan/))
       b(:,2) = reshape((/0.0, 1.001520, 1.012403, 1.001306, 1.001257, 1.000743 /), (/nchan/))
    endif

  end subroutine DEFINE_GSIP_CAL

!-------------------------------------------------------------------------
! Numerical recipes bisection search - x will be between xx(j) and xx(j+1)
!--------------------------------------------------------------------------
  subroutine LOCATE(xx, n, x, j)

!   Arguments
    integer,                        intent(in)  :: n
    integer,                        intent(out) :: j
    real (kind=ipre),               intent(in)  :: x
    real (kind=ipre), dimension(:), intent(in)  :: xx

!   Local variables
    integer :: i, jl, jm, ju

    jl = 0
    ju = n + 1
    do i = 1, 2*n
       if (ju-jl <= 1) then
          exit
       endif
       jm = (ju + jl) / 2
       if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    if (x == xx(1)) then
       j=1
    else if (x == xx(n)) then
       j = n - 1
    else
       j = jl
    endif

  end subroutine LOCATE

!---------------------------------------------------------------------
! Compute solar angles
!
! input:
!         jday = julian day
!         tu   = time of day - fractional hours
!         lon  = latitude in degrees
!         lat  = latitude in degrees
!      
!  output:
!         asol = solar zenith angle in degrees
!         phis = solar azimuth angle in degrees
!
!----------------------------------------------------------------------
  subroutine POSSOL(jday, tu, xlon, xlat, asol, phis)

!   Arguments
    integer,                        intent(in)  :: jday 
    real (kind=ipre),               intent(in)  :: tu
    real (kind=ipre), dimension(:), intent(in)  :: xlon, xlat
    real (kind=ipre), dimension(:), intent(out) :: asol, phis

!   Local variables
    real (kind=ipre):: tsm, xlo, xla, xj, a1, a2, et, tsv, ah, a3, delta, amuzero, elev, &
         az, caz, azim, pi2

    integer :: i, n

    n = size(xlon)

    do i = 1, n

! jday is the number of the day in the month
!
!
!      mean solar time
       tsm = tu + xlon(i)/15.0
       xlo = xlon(i)*dtor
       xla = xlat(i)*dtor
       xj = real(jday)

!      time equation (mn.dec)
       a1 = (1.00554*xj - 6.28306)  * dtor
       a2 = (1.93946*xj + 23.35089) * dtor
       et = -7.67825*sin(a1) - 10.09176*sin(a2)

!      true solar time
       tsv = tsm + et/60.0
       tsv = tsv - 12.0

!      hour angle
       ah = tsv*15.0*dtor

!      solar declination (in radian)
       a3 = (0.9683*xj - 78.00878) * dtor
       delta = 23.4856*sin(a3)*dtor

!     elevation, azimuth
      amuzero = sin(xla)*sin(delta) + cos(xla)*cos(delta)*cos(ah)
      elev = asin(amuzero)
      az = cos(delta)*sin(ah)/cos(elev)
      caz = (-cos(xla)*sin(delta) + sin(xla)*cos(delta)*cos(ah)) / cos(elev)

      if (az >= 1.0) then
         azim = asin(1.0)
      elseif (az <= -1.0) then
         azim = asin(-1.0)
      else
         azim = asin(az)
      endif

      if (caz <= 0.0) then
         azim = pi - azim
      endif

      if ((caz > 0.0) .and. (az <= 0.0)) then
         azim = 2 * pi + azim
      endif
      azim = azim + pi
      pi2 = 2 * pi
      if (azim > pi2) then
         azim = azim - pi2
      endif

!     conversion in degrees
      elev = elev / dtor
      asol(i) = 90.0 - elev
!     phis(i) = azim / dtor - 180.0
      phis(i) = azim / dtor  !akh - try to get 0 - 360

   enddo

 end subroutine POSSOL

!-------------------------------------------------
! Compute julian day (1-365/366)
!
! input:
!         iday   = integer day
!         imonth = integer month
!         iyear  = integer year (2 or 4 digits)
!         
! output:
!         jday  = julian day
!--------------------------------------------------
 subroutine JULIAN(iday,imonth,iyear,jday)

!  arguments
   integer, intent(in) ::  iday, imonth, iyear
   integer, intent(out) :: jday

!  local variables
   integer ::  j
   integer, dimension(12) :: jmonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)

   jday = iday
   if (modulo(iyear, 4) == 0) then
      jmonth(2)=29
   endif

   do j = 1, imonth - 1
      jday = jday + jmonth(j)
   end do

 end subroutine JULIAN

!--------------------------------------------------------------
! Subroutine to make geostationary satellite azimuth field
!   
!     xlon = longitude of the location (positive for western hemisphere) 
!     xlat = latitude of the location  (positive for northern hemisphere) 
!     isat = 0 for GOES-E,  1 for GOES-W, 2 for GOES-9 Pacific
!
!     zenith  = satellite zenith view angle 
!     azimuth = satellite azimuth angle clockwise from north
!--------------------------------------------------------------
 subroutine COMPUTE_SATELLITE_ANGLES(isat, xlon, xlat, zenith, azimuth)

!  arguments
   integer,                         intent(in)  :: isat
   real (kind=real4), dimension(:), intent(in)  :: xlon, xlat
   real (kind=real4), dimension(:), intent(out) :: zenith, azimuth

!  Local variables
   real (kind=real4), parameter :: EASTPOS = -75.0, WESTPOS = -135.0, PACPOS = 155.0
   real (kind=real8):: satlon, satlat, lat, lon, beta
   integer :: i, n

   n = size(xlon) 

   if (isat == 0) then
      satlon = EASTPOS
   elseif (isat == 1) then
      satlon = WESTPOS
   elseif (isat == 2) then
      satlon = PACPOS
   endif
   satlat = 0.0


   do i = 1, n

      lon = (xlon(i) - satlon) * dtor   ! in radians
      lat = (xlat(i) - satlat) * dtor   ! in radians

      beta = acos( cos(lat) * cos(lon) )

!     zenith angle      
      zenith(i) = asin(max(-1.0_real8, min(1.0_real8, &
              42164.0* sin(beta)/ sqrt(1.808e09 - 5.3725e08*cos(beta)))))
      zenith(i) = zenith(i) / dtor

!     azimuth angle
      azimuth(i) = sin(lon) / sin(beta)
      azimuth(i) = min(1.0, max(-1.0,azimuth(i)))
      azimuth(i) = asin(azimuth(i))
      azimuth(i) = azimuth(i) / dtor
      if (lat < 0.0) then
         azimuth(i) = 180.0 - azimuth(i)
      endif
      if (azimuth(i) < 0.0) then
         azimuth(i) = azimuth(i) + 360.0
      endif

   enddo

 end subroutine COMPUTE_SATELLITE_ANGLES

!------------------------------------------------------------------
! a good routine to convert integers to characters and vice versa
!------------------------------------------------------------------
      subroutine ICNVRT(WAY,NUM,STRING,LENGTH,IERR)
!
!       FUNCTION:
!F
!F        This subroutine does an integer-to-character conversion
!F        or a characater-to-integer conversion depending on the
!F        integer WAY:
!F                If WAY = 0 then an integer-to-character conversion
!F                is done. If WAY .NE. 0 then a character-to-integer
!F                conversion is done.
!F
!       USAGE:
!U
!U        CALL ICNVRT(WAY,NUM,STRING)
!U             where WAY, NUM, STRING, and LENGTH are defined below.
!U
!U        Example: CALL ICNVRT(0,1000,STRING,LENGTH)
!U                 on return STRING = '1000' and
!U                 LENGTH = 4.
!U         
!       INPUTS:
!I
!I        WAY - integer::; Determines which way the conversion goes:
!I              if WAY = 0 then an integer-to-character conversion
!I                         is performed;
!I              if WAY.NE.0 then a character-to-integer conversion
!I                         is performed.
!I
!I         NUM - integer::; an input only if WAY = 0. NUM is the integer
!I               number to be converted to a character expression.
!I
!I         STRING - CHARACTER; an input only if WAY .NE. 0. STRING
!I                is the character expression to be converted to an
!I                integer value. It contain no decimal points or 
!I                non-numeric characters other than possibly a
!I                sign. If STRING contains  a '+' sign, it will be
!I                stripped of it on return.
!I
!       OUTPUTS:
!O
!O         NUM - integer::; contains the integer:: representation of 
!O                STRING.
!O
!O         STRING - CHARACTER; contains the CHARACTER representation of
!O                  NUM.
!O
!O         LENGTH - integer::; The length of STRING to the first blank.
!O                  The significant part of STRING can be accessed with
!O                  the declaration STRING(1:LENGTH).
!O
!O         IERR - integer:: variable giving return condition:
!O                IERR = 0 for normal return;
!O                IERR = 1 if NUM cannot be converted to STRING because
!O                       STRING is too short or STRING cannot be
!O                       converted to NUM because STRING is too long.
!O                IERR = 2 if STRING contained a non-numeric character
!O                       other than a leading sign or something went
!O                       wrong with an integer-to-character conversion.
!O
!       ALGORITHM:
!A
!A         Nothing noteworthy, except that this subroutine will work
!A          for strange character sets where the character '1' doesn't
!A          follow '0', etc.
!A
!       MACHINE DEPendENCIES: CM
!M          The parameter MAXINT (below) should be set to the
!M          number of digits that an integer:: data type can have
!M          not including leading signs. For VAX FORTRAN V4.4-177
!M          MAXINT = 10.
!M
!M          NOTE: Under VAX FORTRAN V4.4-177, the
!M          error condition IERR = 1 will never occur for an
!M          integer-to-character conversion if STRING
!M          is allocated at least 11 bytes (CHARACTER*11).
!M
!       HISTORY:
!H
!H      written by:             bobby bodenheimer
!H      date:                   september 1986
!H      current version:        1.0
!H      modifications:          none
!H
!       ROUTINES CALLED:
!C
!C          NONE.
!C
!----------------------------------------------------------------------
!       written for:    The CASCADE Project
!                       Oak Ridge National Laboratory
!                       U.S. Department of Energy
!                       contract number DE-AC05-840R21400
!                       subcontract number 37B-7685 S13
!                       organization:  The University of Tennessee
!----------------------------------------------------------------------
!       THIS SOFTWARE IS IN THE PUBLIC DOMAIN
!       NO RESTRICTIONS ON ITS USE ARE IMPLIED
!----------------------------------------------------------------------
!
! Global Variables.
!
     integer,intent(in):: WAY
      integer, intent(out)::  LENGTH, IERR
      integer, intent(inout):: NUM
      character(len=*), intent(inout):: STRING
!
!
! Local Variables
!
      integer::       I
      integer::       MNUM
      integer::       M
      logical::       NEG
!
      integer, parameter::MAXINT=10
!
      NEG = .FALSE.
      IERR = 0
!
!  Integer-to-character conversion.
!
      if (WAY == 0) then
         STRING = " "
         if (NUM < 0) then
            NEG = .TRUE.
            MNUM = -NUM
            LENGTH = INT(LOG10(REAL(MNUM))) + 1
         else if (NUM == 0) then
            MNUM = NUM
            LENGTH = 1
         else
            MNUM = NUM
            LENGTH = INT(LOG10(REAL(MNUM))) + 1
         end if
         if (LENGTH > LEN(STRING)) then
            IERR = 1
            return
         end if
ten:     do I=LENGTH,1,-1    
            M=INT(REAL(MNUM)/10**(I-1))
            if (M == 0) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "0"
            else if (M == 1) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "1"
            else if (M == 2) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "2"
            else if (M == 3) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "3"
            else if (M == 4) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "4"
            else if (M == 5) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "5"
            else if (M == 6) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "6"
            else if (M == 7) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "7"
            else if (M == 8) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "8"
            else if (M == 9) then
               STRING(LENGTH-I+1:LENGTH-I+1) = "9"
            else
               IERR = 2
               return
            end if
            MNUM = MNUM - M*10**(I-1)
         end do ten

         if (NEG .eqv. .true.) then
            STRING = "-"//STRING
            LENGTH = LENGTH + 1
         end if
!
!  Character-to-integer conversion.
!
      else
         if (STRING(1:1) == "-") then
            NEG = .TRUE.
            STRING = STRING(2:LEN(STRING))
         end if
         if (STRING(1:1) == "+") STRING = STRING(2:LEN(STRING))
         NUM = 0
         LENGTH = INDEX(STRING," ") - 1
         if (LENGTH > MAXINT) then
            IERR = 1
            return
         end if
twenty:  do I=LENGTH,1,-1
            if (STRING(LENGTH-I+1:LENGTH-I+1) == "0") then
               M = 0
            else if (STRING(LENGTH-I+1:LENGTH-I+1) == "1") then
               M = 1
            else if (STRING(LENGTH-I+1:LENGTH-I+1) == "2") then
               M = 2
            else if (STRING(LENGTH-I+1:LENGTH-I+1) == "3") then
               M = 3
            else if (STRING(LENGTH-I+1:LENGTH-I+1) == "4") then
               M = 4
            else if (STRING(LENGTH-I+1:LENGTH-I+1) == "5") then
               M = 5
            else if (STRING(LENGTH-I+1:LENGTH-I+1) == "6") then
               M = 6
            else if (STRING(LENGTH-I+1:LENGTH-I+1) == "7") then
               M = 7
            else if (STRING(LENGTH-I+1:LENGTH-I+1) == "8") then
               M = 8
            else if (STRING(LENGTH-I+1:LENGTH-I+1) == "9") then
               M = 9
            else
               IERR = 2
               return
            end if
            NUM = NUM + INT(10**(I-1))*M
         end do twenty

         if (NEG .eqv. .true.) then
            NUM = -NUM
            STRING = '-'//STRING
            LENGTH = LENGTH + 1
         end if
      end if
!
!  Last lines of ICNVRT
!
   return

end subroutine ICNVRT

!---------------------------------------------
! Compute the month
!---------------------------------------------
 function COMPUTE_MONTH(jday, ileap) result(month)

!  arguments
   integer, intent(in) :: ileap
   integer, intent(in) :: jday
   integer :: month

   month = 0
   if (jday < 32) then
      month = 1
   elseif (jday < 61+ileap) then
      month = 2
   elseif (jday < 91+ileap) then
      month = 3
   elseif (jday < 121+ileap) then
      month = 4
   elseif (jday < 152+ileap) then
      month = 5
   elseif (jday < 182+ileap) then
      month = 6
   elseif (jday < 213+ileap) then
      month = 7
   elseif (jday < 244+ileap) then
      month = 8
   elseif (jday < 274+ileap) then
      month = 9
   elseif (jday < 305+ileap) then
      month = 10
   elseif (jday < 335+ileap) then
      month = 11
   else
      month = 12
   endif

 end function COMPUTE_MONTH

!--------------------------------------------
! Compute the day
!---------------------------------------------
 function COMPUTE_DAY(jday, ileap) result(day)

!  arguments
   integer, intent(in) :: ileap
   integer, intent(in) :: jday
   integer :: day

!  Local variables
   integer :: month

   month = 0
   if (jday < 32) then
      month = 1
      day = jday
   elseif (jday < 60) then
      day = jday - 31
   elseif ((jday == 60).and.(ileap == 1)) then
      day = jday - 31
   elseif ((jday == 60).and.(ileap == 0)) then
      day = jday - 59
   elseif (jday < 91+ileap) then
      day = jday - (59 + ileap)
   elseif (jday < 121+ileap) then
      day = jday - (90 + ileap)
   elseif (jday < 152+ileap) then
      day = jday - (120 + ileap)
   elseif (jday < 182+ileap) then
      day = jday - (151 + ileap)
   elseif (jday < 213+ileap) then
      day = jday - (181 + ileap)
   elseif (jday < 244+ileap) then
      day = jday - (212 + ileap)
   elseif (jday < 274+ileap) then
      day = jday - (243 + ileap)
   elseif (jday < 305+ileap) then
      day = jday - (273 + ileap)
   elseif (jday < 335+ileap) then
      day = jday - (304 + ileap)
   else
      day = jday - (334 + ileap)
   endif

 end function COMPUTE_DAY
!--------------------------------------------
!----------------------------------------------------------------
! functions to compute some needed water vapor parameters
!----------------------------------------------------------------
 function VAPOR(T) result(es)

!  This function computes saturation H2O vapor pressure (over liquid)
!  using Liebe's approximation (corrected).
!
!  T in Kelvin
!  es in mbar
!                   PWR 4/8/92
  implicit none
  real, intent (in) :: T
  real :: es, TH

  TH = 300.0 / T
  es = 35.3 * exp(22.64 * (1.0 - TH)) * TH**5

  return
end function VAPOR


function SATVAPLIQ_OLD(T) result(es)

  implicit none
  real, intent (in) :: T  ! temperature (°C)
  real :: es              ! saturation vapor pressure (Pa)
                          !    with respect to liquid water

  if (T .lt. -273.0)  then
     es = - 999.0
  else
     es = 100.0*(10.**(-7.90298*(373.15/(T+273.15)-1.0)              &
        + 5.02808*log10(373.15/(T+273.15))                           &
        - 1.3816E-07*(10.**(11.334*(1.0-((T+273.15)/373.15)))-1.0)   &
        + 8.1328E-03*(10.**(-3.49149*((373.15/(T+273.15))-1.0))-1.0) &
        + log10(1013.246)))

  endif

  return
end function SATVAPLIQ_OLD


function SATVAPLIQ(T) result(es)
! Teten's Formula.

  implicit none
  real, intent (in) :: T  ! temperature (°C)
  real :: es              ! saturation vapor pressure (Pa)
                          !    with respect to liquid water

  if (T .lt. -273.0)  then
     es = -999.0
  else
     es = 610.78*exp( 17.2693882*T / ( 237.3 + T ))
  endif

  return
end function SATVAPLIQ

function SATVAPLIQB(T) result(es)
! Bolton's (1980) Formula.

  implicit none
  real, intent (in) :: T  ! temperature (°C)
  real :: es              ! saturation vapor pressure (Pa)
                          !    with respect to liquid water

  if (T .lt. -273.0)  then
     es = -999.0
  else
     es = 610.78 * exp( 17.2693882*T / ( 237.3 + T ))
  endif

  return
end function SATVAPLIQB


function SATVAPICE_OLD(T) result(es)
! Goff-Gratch formula (List 1951)

  implicit none
  real, intent (in) :: T  ! temperature (°C)
  real :: es              ! saturation vapor pressure (Pa)
                          !    with respect to liquid water
! real :: SATVAPLIQ

  if (T .lt. -273.0)  then
     es = -999.0
  elseif (T .ge. 0.0) then
     es = SATVAPLIQ(T)
  else
     es = 100.0*(10.0**(-9.09718*(273.15/(T+273.15)-1.0)   &
          - 3.56654*log10(273.15/(T+273.15))               &
          + 0.876793*(1.0-(T+273.15)/273.15)               &
          + log10(6.1071)))
  endif

  return
end function SATVAPICE_OLD



function SATVAPICE(T) result(es)
!  The following is an empirical fit to the Smithsonian tables, by
!  Kerry Emanuel ("Atmospheric Convection", K. Emanuel (1994), Oxford
!  Press, pg 117).  Kerry claims it accurate to 0.14% in the range
!  -80 <= T <= 0C.  I believe e is in mb and T in K, but you may want
!  to spot check a few values to confirm.
!
!       e# = saturation vapor pressure over ice, T = temperature
!
!       ln e# = 23.33086 - (6111.72784/T) + 0.15215 ln T
!
!
!  Dennis Boccippio
!  djboccip@bolt.msfc.nasa.gov
!  http://www-cmpo.mit.edu/~djboccip/Dennis.html

  implicit none
  real, intent (in) :: T  ! temperature (°C)
  real :: es              ! saturation vapor pressure (Pa)
                          !    with respect to liquid water
  real :: TK              ! temperature in Kelvin
! real :: SATVAPliq

  if (T .lt. -273.0)  then
     es = -999.0
  elseif (T .ge. 0.0) then
     es = SATVAPliq(T)
  else
     TK = T + 273.15
     es = 100.0 * exp(23.33086 - (6111.72784/TK) + 0.15215*log(TK))
  endif

  return
end function SATVAPICE


end module NUMERICAL_ROUTINES
