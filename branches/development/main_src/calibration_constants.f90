!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: calibration_constants.f90 (src)
!       CALIBRATION_CONSTANTS (program)
!
! PURPOSE: This module serves as a common block for passing the 
!          instrument and calibration coefficients
!
! DESCRIPTION: 
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
! NOTES:
!   File I/O - none
!
!   Public routines: none
!
!   Private routines: none
!--------------------------------------------------------------------------------------
module CALIBRATION_CONSTANTS

  use CONSTANTS

  implicit none

!--- Planck Constants
  real (kind=real4), parameter, public:: c1 = 1.191062e-5, &
                                         c2 = 1.4387863

  real (kind=real4), public,save:: A1_20,A2_20,Nu_20
  real (kind=real4), public,save:: A1_21,A2_21,Nu_21
  real (kind=real4), public,save:: A1_22,A2_22,Nu_22
  real (kind=real4), public,save:: A1_23,A2_23,Nu_23
  real (kind=real4), public,save:: A1_24,A2_24,Nu_24
  real (kind=real4), public,save:: A1_25,A2_25,Nu_25
  real (kind=real4), public,save:: A1_27,A2_27,Nu_27
  real (kind=real4), public,save:: A1_28,A2_28,Nu_28
  real (kind=real4), public,save:: A1_29,A2_29,Nu_29
  real (kind=real4), public,save:: A1_30,A2_30,Nu_30
  real (kind=real4), public,save:: A1_31,A2_31,Nu_31
  real (kind=real4), public,save:: A1_32,A2_32,Nu_32
  real (kind=real4), public,save:: A1_33,A2_33,Nu_33
  real (kind=real4), public,save:: A1_34,A2_34,Nu_34
  real (kind=real4), public,save:: A1_35,A2_35,Nu_35
  real (kind=real4), public,save:: A1_36,A2_36,Nu_36
  real (kind=real4), public,save:: A1_40,A2_40,Nu_40
  real (kind=real4), public,save:: A1_41,A2_41,Nu_41

!----- variables read in from instrument constant files (planck stored in constants module)
 real (kind=real4), save, public:: B0_3b,B1_3b,B2_3b,B0_4,B1_4,B2_4,B0_5,B1_5,B2_5, &
                      space_Rad_3b,space_Rad_4,space_Rad_5
 real(kind=real4), dimension(0:4,4),public,save:: Prt_Coef
 real(kind=real4), dimension(4),public,save:: Prt_Weight
 real(kind=real4),save,public:: Ch1_Gain_Low_0,Ch1_Gain_High_0,Ch2_Gain_Low_0,&
                Ch2_Gain_High_0,Ch3a_Gain_Low_0,Ch3a_Gain_High_0, &
                Ch1_Degrad_Low_1,Ch2_Degrad_Low_1,Ch3a_Degrad_Low_1, &
                Ch1_Degrad_High_1,Ch2_Degrad_High_1,Ch3a_Degrad_High_1, &
                Ch1_Degrad_Low_2,Ch2_Degrad_Low_2,Ch3a_Degrad_Low_2, &
                Ch1_Degrad_High_2,Ch2_Degrad_High_2,Ch3a_Degrad_High_2, &
                Ch1_Dark_Count,Ch2_Dark_Count,Ch3a_Dark_Count, &
                Ch1_Dark_Count_Cal,Ch2_Dark_Count_Cal,Ch3a_Dark_Count_Cal, &
                Ch1_Switch_Count,Ch2_Switch_Count, Ch3a_Switch_Count, &
                Ch1_Gain_Low,Ch1_Gain_High,Ch2_Gain_Low,Ch2_Gain_High, &
                Ch3a_Gain_Low,Ch3a_Gain_High, Ref_Ch1_Switch,Ref_Ch2_Switch, &
                Ref_Ch6_Switch,Ch1_Switch_Count_Cal,Ch2_Switch_Count_Cal, &
                Ch3a_Switch_Count_Cal

 real(kind=real4),save,public:: sun_Earth_distance
 real(kind=real4),save,public:: launch_date


 real(kind=real4),save,public:: Goes_Ch2_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch3_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch4_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch5_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch6_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch7_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch8_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch9_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch10_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch12_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch13_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch14_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch16_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch17_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch18_Thermal_Intercept
 real(kind=real4),save,public:: Goes_Ch2_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch3_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch4_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch5_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch6_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch7_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch8_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch9_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch10_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch12_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch13_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch14_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch16_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch17_Thermal_Slope
 real(kind=real4),save,public:: Goes_Ch18_Thermal_Slope

 real(kind=real8),save,public:: Goes_Input_Time = 0
 real(kind=real8),save,public:: Goes_Epoch_Time = 0
 real(kind=real8),save,public:: Goes_Sub_Satellite_Longitude
 real(kind=real8),save,public:: Goes_Sub_Satellite_Latitude

 !--- MCSST values
 real(kind=real4),save,public:: B1_day_mask,B2_day_mask, &
                                B3_day_mask,B4_day_mask

 real(kind=real4),save,public:: Solar_Ch20,Solar_Ch20_Nu,Ew_Ch20
 character(len=7),save,public:: Sat_Name

end module CALIBRATION_CONSTANTS
