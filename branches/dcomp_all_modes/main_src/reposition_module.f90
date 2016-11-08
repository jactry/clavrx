! $Id$ 
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: reposition_module.f90 (src)
!       AVHRR_REPOSITION_ROUTINES (program)
!
! PURPOSE: this module houses the non-Nagle routines for repositioning the
!          AVHRR lat and lon values for time corrections
!
! DESCRIPTION:  
!             Note, the mjdn numbers vary in this module.  The values used to
!             record the clock errors are referenced to Wednesday November 17, 1858
!
!             The values used by Fred Nagles routines are referenced to
!             12 Z Januar, 1970
!
!             The offset between the two is 40,587.5
!
!             In INTERPOLATE_CLOCK_ERROR, I will use the standard definition
!             In REPOSITION_FOR_CLOCK_ERROR, I will use Nagle's definition
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
! Public routines used in this module:
! REPOSITION_FOR_CLOCK_ERROR - compute the time for each pixel from the scan value
!--------------------------------------------------------------------------------------
module AVHRR_REPOSITION_ROUTINES
 use CONSTANTS
 use PIXEL_COMMON
 use NUMERICAL_TOOLS_MOD
 
  use date_tools_mod, only: &
         leap_year_fct &
         , compute_month &
         , compute_day
 
 implicit none
 private
 public:: REPOSITION_FOR_CLOCK_ERROR, &
          SETUP_CLOCK_CORRECTIONS, &
          INTERPOLATE_CLOCK_ERROR

 integer, parameter, private:: num_n07 = 15
 real, dimension(num_n07):: start_mjdn_n07, &
                            end_mjdn_n07, &
                            start_clock_error_n07, &
                            end_clock_error_n07

 integer, parameter, private:: num_n09 = 74
 real, dimension(num_n09):: start_mjdn_n09, &
                            end_mjdn_n09, &
                            start_clock_error_n09, &
                            end_clock_error_n09

 integer, parameter, private:: num_n11 = 50
 real, dimension(num_n11):: start_mjdn_n11, &
                            end_mjdn_n11, &
                            start_clock_error_n11, &
                            end_clock_error_n11

 integer, parameter, private:: num_n12 = 14
 real, dimension(num_n12):: start_mjdn_n12, &
                            end_mjdn_n12, &
                            start_clock_error_n12, &
                            end_clock_error_n12

 integer, parameter, private:: num_n14 = 28
 real, dimension(num_n14):: start_mjdn_n14, &
                            end_mjdn_n14, &
                            start_clock_error_n14, &
                            end_clock_error_n14

 integer, parameter, private:: num_n15 = 23
 real, dimension(num_n15):: start_mjdn_n15, &
                            end_mjdn_n15, &
                            start_clock_error_n15, &
                            end_clock_error_n15

 integer, parameter, private:: num_n16 = 12
 real, dimension(num_n16):: start_mjdn_n16, &
                            end_mjdn_n16, &
                            start_clock_error_n16, &
                            end_clock_error_n16

 integer, parameter, private:: num_n17 = 8
 real, dimension(num_n17):: start_mjdn_n17, &
                            end_mjdn_n17, &
                            start_clock_error_n17, &
                            end_clock_error_n17

 integer, parameter, private:: num_n05 = 0
 integer, parameter, private:: num_n06 = 0
 integer, parameter, private:: num_n08 = 0
 integer, parameter, private:: num_n10 = 0
 integer, parameter, private:: num_n18 = 0
 integer, parameter, private:: num_n19 = 0
 integer, parameter, private:: num_m01 = 0
 integer, parameter, private:: num_m02 = 0
 integer, parameter, private:: num_m03 = 0

 real, parameter::  mjdn_offset = 40587.5

 contains
!--------------------------------------------------------------------------
! Setup the clock corrections in memory
!--------------------------------------------------------------------------
 subroutine SETUP_CLOCK_CORRECTIONS()

 !--- NOAA-07
  start_mjdn_n07 = (/ &
                   44783.0, 44973.0, 45023.0, 45151.0, &
                   45360.0, 45516.0, 45528.0, 45781.0, &
                   46042.0, 46066.0, 46173.0, 46202.0, &
                   46247.0, 46370.0, 46579.0 /)

   end_mjdn_n07 = (/ &
                   44973.0, 45023.0, 45151.0, 45360.0, &
                   45516.0, 45528.0, 45781.0, 46042.0, &
                   46066.0, 46173.0, 46202.0, 46247.0, &
                   46370.0, 46579.0, 46589.0/)

   start_clock_error_n07 = (/ &
                    0.10, 0.80,-0.50,-0.60, &
                   -0.10, 1.50,-0.70,-1.10, &
                   -0.90, 0.00, 0.20,-1.10, &
                   -0.00,-0.30,-1.30 /)

   end_clock_error_n07 = (/ &
                    0.80, 1.10, 0.10, 0.30, &
                    0.50, 1.60, 0.30,-0.10, &
                   -0.80, 0.40, 0.30,-1.00, &
                    0.40, 0.30,-1.30/)

 !--- NOAA-09
  start_mjdn_n09 = (/ &
  46444.8, 46481.8, 46620.8, 46739.3, &
  46844.3, 46935.3, 47012.3, 47089.4, &
  47161.0, 47222.1, 47285.1, 47348.0, &
  47411.1, 47467.0, 47523.1, 47550.4, &
  47621.2, 47670.1, 47719.2, 47769.1, &
  47817.1, 47865.1, 47892.1, 47950.1, &
  47999.1, 48041.1, 48083.2, 48125.2, &
  48167.1, 48212.0, 48257.1, 48279.1, &
  48321.1, 48363.1, 48399.0, 48433.1, &
  48469.2, 48526.2, 48545.2, 48587.2, &
  48629.3, 48671.2, 48712.2, 48755.3, &
  48804.6, 48832.5, 48868.0, 48909.3, &
  48951.0, 48993.3, 49035.3, 49077.0, &
  49124.1, 49154.0, 49169.6, 49230.1, &
  49280.1, 49330.1, 49378.1, 49413.3, &
  49455.1, 49490.1, 49534.1, 49562.1, &
  49588.1, 49630.1, 49658.1, 49700.0, &
  49737.0, 49770.0, 49805.1, 49833.1, &
  49868.1, 49919.1/)

  end_mjdn_n09 = (/ &
  46481.4, 46620.4, 46738.8, 46843.9, &
  46934.8, 47011.9, 47088.9, 47161.0, &
  47222.0, 47285.0, 47347.9, 47411.0, &
  47467.0, 47523.0, 47549.5, 47620.5, &
  47669.1, 47718.5, 47769.0, 47816.2, &
  47864.4, 47891.9, 47949.4, 47998.5, &
  48040.9, 48082.9, 48124.5, 48166.5, &
  48205.4, 48257.0, 48278.9, 48321.0, &
  48362.2, 48398.9, 48432.2, 48467.1, &
  48488.0, 48544.9, 48587.0, 48629.0, &
  48671.0, 48711.8, 48754.9, 48804.0, &
  48828.5, 48857.6, 48909.0, 48950.6, &
  48993.0, 49035.0, 49076.7, 49123.8, &
  49153.8, 49169.3, 49200.8, 49279.8, &
  49328.9, 49377.9, 49412.8, 49453.1, &
  49489.8, 49534.0, 49561.8, 49587.9, &
  49629.9, 49657.9, 49699.6, 49736.9, &
  49769.8, 49804.8, 49832.9, 49867.9, &
  49918.8, 49931.8/)

  start_clock_error_n09 = (/ &
 -0.53, 0.45, 0.46, 0.47, 0.40, 0.39, &
  0.48, 0.48, 0.21, 0.36, 0.39, 0.40, &
  0.32, 0.35, 0.37, 0.75, 0.31, 0.36, &
  0.46, 0.45, 0.45, 0.41, 0.69, 0.36, &
  0.25, 0.31, 0.36, 0.40, 0.42, 0.33, &
 -0.07, 0.52, 0.47, 0.39, 0.16, 0.31, &
  0.39,-0.24, 0.36, 0.44, 0.52, 0.58, &
  0.75, 0.70,-0.01, 0.45, 0.71, 0.68, &
  0.69, 0.63, 0.57, 0.52, 0.25, 0.62, &
  1.02, 1.22, 0.77, 0.83, 0.42, 0.53, &
  0.35, 0.48, 0.20, 0.07, 0.53, 0.30, &
  0.67, 0.42, 0.41, 0.50, 0.56, 0.89, &
  0.91, 0.74/)

  end_clock_error_n09 = (/ &
 -0.80,-0.76,-0.80,-0.84,-0.85,-0.75, &
 -0.77,-0.79,-0.89,-0.87,-0.87,-0.93, &
 -0.89,-0.88,-0.25,-0.91,-0.85,-0.79, &
 -0.80,-0.79,-0.81,-0.30,-0.87,-0.99, &
 -0.93,-0.85,-0.81,-0.79,-0.72,-1.03, &
 -0.71,-0.76,-0.82,-0.72,-0.88,-0.78, &
 -0.21,-0.86,-1.03,-0.95,-0.90,-0.84, &
 -0.77,-0.99,-0.87,-0.42,-0.78,-0.79, &
 -0.87,-0.91,-0.96,-1.26,-0.88, 0.00, &
 -0.18,-0.72,-1.15,-1.04,-0.95,-1.03, &
 -1.02,-1.28,-0.94,-1.01,-1.15,-0.82, &
 -1.06,-1.11,-0.97,-0.90,-0.62,-0.57, &
 -1.24, 0.21/)

!---- NOAA-11
  start_mjdn_n11 = (/ &
  47430.7, 47523.0, 47635.0, 47810.0, &
  47892.3, 47978.0, 48049.0, 48097.0, &
  48105.3, 48125.0, 48188.0, 48237.3, &
  48294.1, 48335.0, 48384.1, 48433.0, &
  48489.1, 48538.0, 48587.0, 48636.1, &
  48678.1, 48727.1, 48769.1, 48832.0, &
  48874.0, 48930.4, 48958.0, 49000.0, &
  49035.0, 49077.1, 49126.0, 49169.1, &
  49222.0, 49259.0, 49301.1, 49343.1, &
  49385.1, 49420.0, 49462.0, 49497.0, &
  49534.0, 49561.9, 49574.0, 49609.0, &
  49651.0, 49693.0, 49730.8, 49770.0, &
  49812.0, 49925.7/)

  end_mjdn_n11 = (/ &
  47523.0, 47634.9, 47809.9, 47892.0, &
  47978.0, 48049.0, 48091.8, 48101.0, &
  48125.0, 48188.0, 48236.9, 48294.0, &
  48334.9, 48384.0, 48432.9, 48489.0, &
  48537.9, 48587.0, 48636.0, 48678.0, &
  48727.0, 48769.0, 48832.0, 48874.0, &
  48929.9, 48958.0, 49000.0, 49034.9, &
  49076.9, 49125.9, 49169.0, 49222.0, &
  49258.9, 49301.0, 49343.0, 49385.0, &
  49420.0, 49461.9, 49496.9, 49534.0, &
  49561.7, 49574.0, 49608.9, 49650.9, &
  49692.9, 49729.7, 49770.0, 49812.0, &
  49925.7, 50086.7/)

  start_clock_error_n11 = (/ &
  0.10,-0.45,-1.21,-0.99,-0.81,-0.48, &
 -0.33, 0.68, 0.10,-0.47,-0.35,-0.44, &
 -0.38,-0.59,-0.64,-0.67,-0.54,-0.57, &
 -0.56,-0.53,-0.64,-0.55,-0.63,-0.26, &
 -0.32,-0.05,-0.42,-0.46,-0.67,-0.68, &
 -0.55,-0.50,-0.21,-0.33,-0.33,-0.28, &
 -0.24,-0.38,-0.34,-0.43,-1.51, 0.23, &
 -0.45,-0.59,-0.47,-0.38,-0.36,-0.33, &
 -0.23, 2.81/)

 end_clock_error_n11 = (/ &
  0.55, 0.55, 0.92, 0.19, 0.51, 0.67, &
  0.40, 0.76, 0.53, 0.66, 0.54, 0.63, &
  0.42, 0.36, 0.35, 0.47, 0.44, 0.45, &
  0.47, 0.36, 0.46, 0.38, 0.74, 0.70, &
  0.96, 0.59, 0.56, 0.36, 0.33, 0.47, &
  0.51, 0.78, 0.68, 0.70, 0.73, 0.77, &
  0.64, 0.69, 0.57, 0.50,-0.81, 0.54, &
  0.45, 0.53, 0.65, 0.60, 0.67, 0.77, &
  2.81, 7.36/)

!--- NOAA-12
  start_mjdn_n12 = (/ &
  48420.0, 48443.0, 48494.7, 48538.0, &
  48712.0, 48800.8, 48988.0, 49413.1, &
  49534.2, 49562.0, 49847.1, 50190.1, &
  50630.0, 50899.0/)

  end_mjdn_n12 = (/ &
  48443.0, 48494.6, 48538.0, 48711.9, &
  48800.5, 48987.6, 49413.0, 49534.1, &
  49562.0, 49847.0, 50190.0, 50629.9, &
  50899.0, 51178.0/)
   
  start_clock_error_n12 = (/ &
  0.02,-1.22,-0.99,-1.21,-1.16,-0.36, &
 -0.34,-0.42,-1.17,-0.09,-0.41, 0.07, &
  0.68, 0.07/)

  end_clock_error_n12 = (/ &
  0.03,-1.17,-1.21,-1.20,-1.14,-0.33, &
  0.57,-0.17,-1.11, 0.61, 0.45,-0.32, &
  0.07,-0.12/)

 !--- NOAA-14
  start_mjdn_n14 = (/ &
    49730.9, 49931.1, 50083.0, 50162.1, &
    50305.0, 50372.0, 50575.0, 50630.1, &
    50645.0, 50680.0, 50750.1, 50925.1, &
    50995.0, 51156.0, 51178.9, 51310.0, &
    51401.0, 51489.0, 51562.0, 51646.0, &
    51786.0, 51919.0, 51990.0, 52150.0, &
    52290.0, 52402.0, 52479.0, 52640.0/)

  end_mjdn_n14 = (/ &
    49931.0, 50083.0, 50162.0, 50305.0, &
    50371.9, 50574.9, 50630.0, 50645.0, &
    50680.0, 50750.0, 50925.0, 50994.9, &
    51156.0, 51178.9, 51309.9, 51401.0, &
    51489.0, 51562.0, 51646.0, 51786.0, &
    51919.0, 51990.0, 52150.0, 52290.0, &
    52402.0, 52479.0, 52640.0, 52981.0/)

   start_clock_error_n14 = (/ &
   -0.37,-0.93, 0.53,-0.43, 0.07,-0.93, &
    0.18, 1.63, 0.25, 0.03,-0.90, 0.04, &
   -0.88,-0.05,-0.34, 0.22,-0.52,-0.65, &
   -0.50,-1.32,-1.18,-1.10,-1.56,-1.22, &
   -1.22,-1.35,-0.72,-0.81/)

   end_clock_error_n14 = (/ &
    0.59, 0.04, 1.07, 0.59, 0.58, 0.69, &
    0.64, 1.75, 0.55, 0.63, 0.55, 0.62, &
    0.45, 0.14, 0.72, 0.97, 0.15,-0.02, &
    0.22,-0.18,-0.10,-0.56,-0.25,-0.21, &
   -0.38,-0.72,-0.81,-1.00/)

!---- NOAA-15
  start_mjdn_n15 = (/ &
  50977.0, 50995.0, 51121.0, 51197.0, &
  51268.0, 51436.0, 51593.0, 51730.0, &
  51946.0, 52066.0, 52129.0, 52193.0, &
  52311.0, 52416.0, 52511.0, 52585.0, &
  52640.0, 52771.0, 52790.0, 52801.0, &
  52964.0, 52964.0, 52975.0/)

  end_mjdn_n15 = (/ &
  50995.0, 51121.0, 51197.0, 51268.0, &
  51436.0, 51593.0, 51730.0, 51946.0, &
  52066.0, 52129.0, 52192.0, 52311.0, &
  52416.0, 52511.0, 52585.0, 52640.0, &
  52771.0, 52790.0, 52801.0, 52963.9, &
  52964.0, 52975.0, 53283.0/)

  start_clock_error_n15 = (/ &
  1.11,-0.25, 0.38, 0.05,-0.81, 0.19, &
  0.03,-1.10,-0.30,-0.33,-0.79,-1.25, &
 -1.26,-1.40,-0.65,-0.92,-0.93,-0.20, &
 -0.97,-0.86,-0.34,-0.10,-0.89/)

  end_clock_error_n15 = (/ &
  1.26, 0.88, 1.05, 0.66, 0.70, 1.51, &
  0.90, 0.69, 0.70, 0.27,-0.24,-0.26, &
 -0.40,-0.65,-0.64,-0.93,-0.97,-0.17, &
 -0.97,-0.89,-0.34,-0.10,-1.01/)

!--- NOAA-16
  start_mjdn_n16 = (/ &
  51807.0, 51809.0, 51898.0, 52129.0, &
  52220.0, 52302.0, 52416.0, 52514.0, &
  52640.0, 52763.0, 52801.0, 52981.0/)

  end_mjdn_n16 = (/ &
  51809.0, 51898.0, 52129.0, 52220.0, &
  52302.0, 52416.0, 52514.0, 52640.0, &
  52763.0, 52801.0, 52981.0, 53283.0/)

  start_clock_error_n16 = (/ &
 -2.71, 0.44,-0.38,-0.60,-1.26,-0.89, &
 -1.43,-0.97,-0.95,-0.87,-0.95,-0.58/)

  end_clock_error_n16 = (/ &
 -2.80, 0.61, 0.38,-0.27,-0.89,-0.43, &
 -0.97,-0.95,-0.87,-0.86,-0.87,-0.56/)

!--- NOAA-17
  start_mjdn_n17 = (/ &
  52451.0, 52476.0, 52619.0, 52726.0, &
  52801.0, 52936.0, 52997.0, 53060.0/)

  end_mjdn_n17 = (/ &
  52476.0, 52619.0, 52726.0, 52801.0, &
  52936.0, 52997.0, 53060.0, 53261.0/)

  start_clock_error_n17 = (/ &
 -0.91,-0.93,-0.93,-0.91,-0.95,-0.91, &
 -0.96,-0.87/)

  end_clock_error_n17 = (/ &
 -0.93,-0.83,-0.82,-0.85,-0.85,-0.87, &
 -0.99,-0.97/)

 end subroutine SETUP_CLOCK_CORRECTIONS

!--------------------------------------------------------------------------
! Interpolate the clock error
!--------------------------------------------------------------------------
 subroutine INTERPOLATE_CLOCK_ERROR(start_year, start_itime, &
                                   end_year, end_itime, &
                                   sat_number,clock_error)

   integer(kind=int2), intent(in) :: start_year, &
                                     end_year

   integer(kind=int4), intent(in) :: start_itime, &
                                     end_itime

   integer(kind=int4), intent(in) :: sat_number

   real(kind=real4), intent(out) :: clock_error

   real(kind=real4), dimension(:), allocatable:: start_mjdn,  &
                                                 end_mjdn, &
                                                 start_clock_error, &
                                                 end_clock_error

   integer:: ibin
   integer:: day, year, hour, minute,  &
             start_month, end_month, start_dom, end_dom,ileap
   integer:: num
   integer(kind=int4):: month, isecond
   real(kind=real4):: second
   real(kind=real4):: temp
   real(kind=real4):: orbit_start_mjdn
   real(kind=real4):: orbit_end_mjdn
   real(kind=real4):: orbit_mean_mjdn
   real(kind=real8):: DTMJDN

!-------------------------------------------------------------------------------


    num = 0
    if (sat_number == 5) then
          num = num_n05
    elseif (sat_number == 6) then
          num = num_n06
    elseif (sat_number == 7) then
          num = num_n07
    elseif (sat_number == 8) then
          num = num_n08
    elseif (sat_number == 9) then
          num = num_n09
    elseif (sat_number == 10) then
          num = num_n10
    elseif (sat_number == 11) then
          num = num_n11
    elseif (sat_number == 12) then
          num = num_n12
    elseif (sat_number == 14) then
          num = num_n14
    elseif (sat_number == 15) then
          num = num_n15
    elseif (sat_number == 16) then
          num = num_n16
    elseif (sat_number == 17) then
          num = num_n17
    elseif (sat_number == 18) then
          num = num_n18
    elseif (sat_number == 19) then
          num = num_n19
    elseif (sat_number == 2) then
          num = num_m02
    elseif (sat_number == 1) then
          num = num_m01
    elseif (sat_number == 3) then
          num = num_m03
    endif

    
    !--- check to see if there are any correction data for this sensor
    if (num > 0) then
        allocate(start_mjdn(num))
        allocate(end_mjdn(num))
        allocate(start_clock_error(num))
        allocate(end_clock_error(num))
     else      !if not, then set error to missing and return
        clock_error = missing_value_real4
        return
     endif

    !--- copy sensor specific data to local arrays
    if (sat_number == 7) then
         start_mjdn = start_mjdn_n07 
         end_mjdn = end_mjdn_n07 
         start_clock_error = start_clock_error_n07 
         end_clock_error = end_clock_error_n07 
    elseif (sat_number == 9) then
         start_mjdn = start_mjdn_n09
         end_mjdn = end_mjdn_n09
         start_clock_error = start_clock_error_n09
         end_clock_error = end_clock_error_n09
    elseif (sat_number == 11) then
         start_mjdn = start_mjdn_n11
         end_mjdn = end_mjdn_n11
         start_clock_error = start_clock_error_n11
         end_clock_error = end_clock_error_n11
    elseif (sat_number == 12) then
         start_mjdn = start_mjdn_n12
         end_mjdn = end_mjdn_n12
         start_clock_error = start_clock_error_n12
         end_clock_error = end_clock_error_n12
    elseif (sat_number == 14) then
         start_mjdn = start_mjdn_n14
         end_mjdn = end_mjdn_n14
         start_clock_error = start_clock_error_n14
         end_clock_error = end_clock_error_n14
    elseif (sat_number == 15) then
         start_mjdn = start_mjdn_n15
         end_mjdn = end_mjdn_n15
         start_clock_error = start_clock_error_n15
         end_clock_error = end_clock_error_n15
    elseif (sat_number == 16) then
         start_mjdn = start_mjdn_n16
         end_mjdn = end_mjdn_n16
         start_clock_error = start_clock_error_n16
         end_clock_error = end_clock_error_n16
    elseif (sat_number == 17) then
         start_mjdn = start_mjdn_n17
         end_mjdn = end_mjdn_n17
         start_clock_error = start_clock_error_n17
         end_clock_error = end_clock_error_n17
    endif

    !--- convert time and date
    ileap = 0
    ileap = leap_year_fct(int(start_year,kind=int4))
    start_month = COMPUTE_MONTH(int(Image%Start_Doy,kind=real4),ileap)
    end_month = COMPUTE_MONTH(int(Image%End_Doy,kind=int4),ileap)
    start_dom = COMPUTE_DAY(int(Image%Start_Doy,kind=int4),ileap)
    end_dom = COMPUTE_DAY(int(Image%End_Doy,kind=int4),ileap)

    !--- determine MJDN at start of orbit
    day = start_dom
    month = start_month
    year = start_year
    temp = start_itime / 60.0 / 60.0/ 1000.0
    hour = int(temp)
    temp = (temp - hour)*60.0
    minute = int(temp)
    second = (temp - minute)*60.0
    isecond=int(second)
    orbit_start_mjdn = &
         real(DTMJDN(year,month,day,hour,minute,isecond),kind=real4) + &
         mjdn_offset

    !--- determine MJDN at end of orbit
    day = end_dom
    month = end_month
    year = end_year
    temp = end_itime / 60.0 / 60.0/ 1000.0
    hour = int(temp)
    temp = (temp - hour)*60.0
    minute = int(temp)
    second = (temp - minute)*60.0
    isecond=int(second)
    orbit_end_mjdn =  &
         real(DTMJDN(year,month,day,hour,minute,isecond),kind=real4) + &
         mjdn_offset
    orbit_mean_mjdn = 0.5*(orbit_start_mjdn + orbit_end_mjdn)

    !---- check to see if time is within the timespan of the clock
    !---- correction data for this satellite
    if ( (orbit_mean_mjdn >= end_mjdn(num)) .or. &
          (orbit_mean_mjdn <= start_mjdn(1))) then

        clock_error = missing_value_real4

        !--- deallocate the arrays
        if (allocated(start_mjdn)) deallocate(start_mjdn)
        if (allocated(end_mjdn)) deallocate(end_mjdn)
        if (allocated(start_clock_error)) deallocate(start_clock_error)
        if (allocated(end_clocK_error)) deallocate(end_clock_error)

        return

     endif

    

   !------ interpolate start time

   !--- determine which time period bin to use
   do ibin = 1, num
     if ( (orbit_mean_mjdn >= start_mjdn(ibin)) .and. &
          (orbit_mean_mjdn <= end_mjdn(ibin))) then
          exit
     endif
   enddo

   !---- check for a successful bin search
   if ((ibin == num+1)) then
           clock_error = missing_value_real4

   else   !attempt interpolaton

       clock_error = start_clock_error(ibin) + &
            (orbit_mean_mjdn - start_mjdn(ibin)) * &
            (end_clock_error(ibin) - start_clock_error(ibin)) / &
            (end_mjdn(ibin) - start_mjdn(ibin))

   endif

   !--- deallocate the arrays
   if (allocated(start_mjdn)) deallocate(start_mjdn)
   if (allocated(end_mjdn)) deallocate(end_mjdn)
   if (allocated(start_clock_error)) deallocate(start_clock_error)
   if (allocated(end_clocK_error)) deallocate(end_clock_error)


 end subroutine INTERPOLATE_CLOCK_ERROR
!--------------------------------------------------------------------------
! Reposition for a clock error
!
! compute the time for each pixel from the scanline time
!
!--------------------------------------------------------------------------
subroutine REPOSITION_FOR_CLOCK_ERROR(j1,j2,timerr,error_flag)

   integer, intent(in):: j1,j2
   real (kind=real4), intent(in):: timerr
   integer, intent(out):: error_flag
   integer:: i,j
   real, parameter:: senzen_max = 55.3 !degrees
   real:: alpha
   real:: second
   real:: temp
   integer:: day, year, hour, minute,  &
             start_month, end_month, start_dom, end_dom,ileap
   integer (kind=int4):: Pixel_offset, Pixel_spacing, ipixel
   integer(kind=int4):: month, isecond
   real (kind=real4), dimension(Image%Number_Of_Elements,j2):: lat_temp,lon_temp
   real (kind=real8), dimension(Image%Number_Of_Elements,j2):: Pixel_time_mjdn
   integer, dimension(j2):: valid_scan_number
   integer, dimension(j2):: valid_scan_index
   integer:: num_valid_scans, jj

   real (kind=real8), dimension(j1:j2-j1+1):: scan_time_mjdn
   real (kind=real8):: time_between_scans

   double precision:: DTMJDN


!--- initialize lat and lon from level-1b values
    nav % lat = nav % lat_1b
    nav % lon = nav % lon_1b

    valid_scan_index = 0
    scan_time_mjdn = missing_value_real4

!--- a convenient parameter
    alpha = 2.0*senzen_max / 360.0 / Image%Number_Of_Elements

!--- compute date terms
    ileap = 0
    if ((Image%start_year == 1984) .or. &
        (Image%start_year == 1988) .or. &
        (Image%start_year == 1992) .or. &
        (Image%start_year == 1996) .or. &
        (Image%start_year == 2000) .or. &
        (Image%start_year == 2004) .or. &
        (Image%start_year == 2008) .or. &
        (Image%start_year == 2012)) ileap = 1

    start_month = COMPUTE_MONTH(int(Image%Start_Doy,kind=real4),ileap)
    end_month = COMPUTE_MONTH(int(Image%End_Doy,kind=int4),ileap)
    start_dom = COMPUTE_DAY(int(Image%Start_Doy,kind=int4),ileap)
    end_dom = COMPUTE_DAY(int(Image%End_Doy,kind=int4),ileap)
    
!--- convert scan line times to MJDN
 jj = 1

 valid_scan_loop: do j = j1,j1+j2-1
     
!--- check for bad scans
      if ((scan_time(j) == missing_value_int4) .or. &
          (bad_scan_flag(j) == sym%YES))  then
           cycle 
      endif

      valid_scan_index(jj) = j
      valid_scan_number(jj) = scan_number(j)

!--- determine date of scan
      if (scan_time(j) >= image%start_time) then
         day = start_dom
         month = start_month
         year = image%start_year
      else
         day = end_dom
         month = end_month
         year = image%end_year
      endif

!--- determine hour, minute and second
      temp = scan_time(j) / 60.0 / 60.0/ 1000.0

      hour = int(temp)
      temp = (temp - hour)*60.0
      minute = int(temp)
      second = (temp - minute)*60.0

!--- determine MJDN of this scan
      isecond=int(second)
      scan_time_mjdn(jj) = DTMJDN(year,month,day,hour,minute,isecond)

      jj = jj + 1

   end do valid_scan_loop

   num_valid_scans = jj - 1

!  print *, "num_valid_scans = ", num_valid_scans

!-- set parameters based on data resolution
Pixel_Offset = 0
Pixel_Spacing = 1

if (AVHRR_GAC_Flag == sym%YES) then
        Pixel_offset = 8
        Pixel_spacing = 5
endif

!--- compute pixel level time
Pixel_time_loop:   do jj = 1, num_valid_scans

   !--- determine mean time between lines
    if (jj == num_valid_scans) then
       time_between_scans = (scan_time_mjdn(jj) - scan_time_mjdn(jj-1)) /  &
                            (valid_scan_number(jj) - valid_scan_number(jj-1))
    else
       time_between_scans = (scan_time_mjdn(jj+1) - scan_time_mjdn(jj)) /  &
                            (valid_scan_number(jj+1) - valid_scan_number(jj))
    endif

   !--- interpolate time for each pixel in scan
   ipixel = Pixel_offset
   do i = 1, Image%Number_Of_Elements
      Pixel_time_mjdn(i,jj) = scan_time_mjdn(jj) + (ipixel-1) * alpha * time_between_scans
      ipixel = ipixel + Pixel_spacing
   enddo
   lat_temp(:,jj) = nav % lat(:,valid_scan_index(jj))
   lon_temp(:,jj) = nav % lon(:,valid_scan_index(jj))

  end do Pixel_time_loop

!--- do repositioning

        !--- assuming reposnx treats first array element as 1
        call REPOSN (Pixel_time_mjdn(:,1:num_valid_scans), &
                     lat_temp(:,1:num_valid_scans), &
                     lon_temp(:,1:num_valid_scans), &
                     timerr,Image%Number_Of_Elements,num_valid_scans)

        !--- initialize output lat, lon to missing
        nav % lat = missing_value_real4
        nav % lon = missing_value_real4

        !--- populate output lat, lon with valid scans
        do jj = 1, num_valid_scans
         nav % lat(:,valid_scan_index(jj)) = lat_temp(:,jj)
         nav % lon(:,valid_scan_index(jj)) = lon_temp(:,jj)
!         print *, "nadir valid = ", bad_scan_flag(valid_scan_index(jj)),  &
!                       lat(205,valid_scan_index(jj)), lat_temp(205,jj), &
!                       Pixel_time_mjdn(205,jj)
        enddo

!  do j = j1,j1+j2-1
!  if ((bad_scan_flag(j) == sym%YES) .or. &
!      ((lat_1b(2005,j) > 44.0).and.(lat_1b(2005,j) < 46.0))) then
!     print *, "nadir lat = ", bad_scan_flag(j), lat_1b(205,j), lat(205,j),  &
!       lat(205,j) - lat_1b(205,j)
!   endif
!  enddo

!--- assign error flag - ask Fred how to this.
     error_flag = sym%NO

end subroutine REPOSITION_FOR_CLOCK_ERROR

end module AVHRR_REPOSITION_ROUTINES
