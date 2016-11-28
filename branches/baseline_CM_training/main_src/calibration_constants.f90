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

  use CONSTANTS, only: &
  real4 & 
  , real8 &
  , nchan_clavrx

  implicit none

!--- Planck Constants
  real (kind=real4), parameter, public:: c1 = 1.191062e-5, &
                                         c2 = 1.4387863

  real (kind=real4), dimension(20:Nchan_Clavrx), public, save:: Planck_A1 = 0.0
  real (kind=real4), dimension(20:Nchan_Clavrx), public, save:: Planck_A2 = 1.0
  real (kind=real4), dimension(20:Nchan_Clavrx), public, save:: Planck_Nu = 0.0
  real (kind=real4), public, save :: Planck_A1_375um_Sndr = 0.0
  real (kind=real4), public, save :: Planck_A2_375um_Sndr = 1.0
  real (kind=real4), public, save :: Planck_Nu_375um_Sndr = 0.0
  real (kind=real4), public, save :: Planck_A1_11um_Sndr = 0.0
  real (kind=real4), public, save :: Planck_A2_11um_Sndr = 1.0
  real (kind=real4), public, save :: Planck_Nu_11um_Sndr = 0.0
  real (kind=real4), public, save :: Planck_A1_12um_Sndr = 0.0
  real (kind=real4), public, save :: Planck_A2_12um_Sndr = 1.0
  real (kind=real4), public, save :: Planck_Nu_12um_Sndr = 0.0

!----- variables read in from instrument constant files (planck stored in constants module)
 real (kind=real4), save, public:: B0_3b = 0.0
 real (kind=real4), save, public:: B1_3b = 0.0
 real (kind=real4), save, public:: B2_3b = 0.0
 real (kind=real4), save, public:: B0_4 = 0.0
 real (kind=real4), save, public:: B1_4 = 0.0
 real (kind=real4), save, public:: B2_4 = 0.0
 real (kind=real4), save, public:: B0_5 = 0.0
 real (kind=real4), save, public:: B1_5 = 0.0
 real (kind=real4), save, public:: B2_5 = 0.0
 real (kind=real4), save, public:: Space_Rad_3b = 0.0
 real (kind=real4), save, public:: Space_Rad_4 = 0.0
 real (kind=real4), save, public:: Space_Rad_5 = 0.0

 real(kind=real4), dimension(0:4,4),public,save:: Prt_Coef = 0.0
 real(kind=real4), dimension(4),public,save:: Prt_Weight = 0.0

 real(kind=real4),save,public:: Ch1_Gain_Low = 0.0
 real(kind=real4),save,public:: Ch1_Gain_High = 0.0
 real(kind=real4),save,public:: Ch1_Gain_Low_0 = 0.0
 real(kind=real4),save,public:: Ch1_Gain_High_0 = 0.0
 real(kind=real4),save,public:: Ch1_Degrad_Low_1 = 0.0
 real(kind=real4),save,public:: Ch1_Degrad_High_1 = 0.0
 real(kind=real4),save,public:: Ch1_Degrad_Low_2 = 0.0
 real(kind=real4),save,public:: Ch1_Degrad_High_2 = 0.0
 real(kind=real4),save,public:: Ch1_Dark_Count = 0.0
 real(kind=real4),save,public:: Ch1_Dark_Count_Cal = 0.0
 real(kind=real4),save,public:: Ch1_Switch_Count = 0.0
 real(kind=real4),save,public:: Ch1_Switch_Count_Cal = 0.0
 real(kind=real4),save,public:: Ref_Ch1_Switch = 0.0

 real(kind=real4),save,public:: Ch2_Gain_Low = 0.0
 real(kind=real4),save,public:: Ch2_Gain_High = 0.0
 real(kind=real4),save,public:: Ch2_Gain_Low_0 = 0.0
 real(kind=real4),save,public:: Ch2_Gain_High_0 = 0.0
 real(kind=real4),save,public:: Ch2_Degrad_Low_1 = 0.0
 real(kind=real4),save,public:: Ch2_Degrad_High_1 = 0.0
 real(kind=real4),save,public:: Ch2_Degrad_Low_2 = 0.0
 real(kind=real4),save,public:: Ch2_Degrad_High_2 = 0.0
 real(kind=real4),save,public:: Ch2_Dark_Count = 0.0
 real(kind=real4),save,public:: Ch2_Dark_Count_Cal = 0.0
 real(kind=real4),save,public:: Ch2_Switch_Count = 0.0
 real(kind=real4),save,public:: Ch2_Switch_Count_Cal = 0.0
 real(kind=real4),save,public:: Ref_Ch2_Switch = 0.0

 real(kind=real4),save,public:: Ch3a_Gain_Low = 0.0
 real(kind=real4),save,public:: Ch3a_Gain_High = 0.0
 real(kind=real4),save,public:: Ch3a_Gain_Low_0 = 0.0
 real(kind=real4),save,public:: Ch3a_Gain_High_0 = 0.0
 real(kind=real4),save,public:: Ch3a_Degrad_Low_1 = 0.0
 real(kind=real4),save,public:: Ch3a_Degrad_High_1 = 0.0
 real(kind=real4),save,public:: Ch3a_Degrad_Low_2 = 0.0
 real(kind=real4),save,public:: Ch3a_Degrad_High_2 = 0.0
 real(kind=real4),save,public:: Ch3a_Dark_Count = 0.0
 real(kind=real4),save,public:: Ch3a_Dark_Count_Cal = 0.0
 real(kind=real4),save,public:: Ch3a_Switch_Count = 0.0
 real(kind=real4),save,public:: Ch3a_Switch_Count_Cal = 0.0
 real(kind=real4),save,public:: Ref_Ch6_Switch = 0.0

 real(kind=real4),save,public:: Sun_Earth_Distance = 0.0
 real(kind=real4),save,public:: Launch_Date = 0.0



 real(kind=real8),save,public:: Goes_Input_Time = 0
 real(kind=real8),save,public:: Goes_Epoch_Time = 0

 real(kind=real4),save,public:: Solar_Ch20 = 0.0
 real(kind=real4),save,public:: Solar_Ch20_Nu = 0.0
 real(kind=real4),save,public:: Ew_Ch20 = 0.0
 character(len=7),save,public:: Sat_Name = ' '

end module CALIBRATION_CONSTANTS
