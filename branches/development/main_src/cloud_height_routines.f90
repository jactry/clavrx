!$Id: cloud_height_routines.f90,v 1.9.2.2 2014/01/26 04:48:32 heidinger Exp $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: cloud_height_routines.f90 (src)
!       CLOUD_HEIGHT_ROUTINES (program)
!
! PURPOSE: This module houses the routines associated with...
!          Cloud Algorithms other than ACHA - AWG Cloud Height Algorithm
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
!--------------------------------------------------------------------------------------
module CLOUD_HEIGHT_ROUTINES
  use CONSTANTS
  use PIXEL_COMMON
  use NWP_COMMON
  use RTM_COMMON
  use NUMERICAL_ROUTINES

  implicit none

  PUBLIC::  OPAQUE_CLOUD_HEIGHT   
  PUBLIC::  COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND
  PUBLIC::  COMPUTE_ALTITUDE_FROM_PRESSURE


  !--- include parameters for each system here
  INTEGER(KIND=INT4), PRIVATE, PARAMETER :: Chan_Idx_375um = 20  !channel number for 3.75 micron
  INTEGER(KIND=INT4), PRIVATE, PARAMETER :: Chan_Idx_67um = 27   !channel number for 6.7 micron
  INTEGER(KIND=INT4), PRIVATE, PARAMETER :: Chan_Idx_73um = 28   !channel number for 7.3 micron
  INTEGER(KIND=INT4), PRIVATE, PARAMETER :: Chan_Idx_85um = 29   !channel number for 8.5 micron
  INTEGER(KIND=INT4), PRIVATE, PARAMETER :: Chan_Idx_11um = 31   !channel number for 11 micron
  INTEGER(KIND=INT4), PRIVATE, PARAMETER :: Chan_Idx_12um = 32   !channel number for 12 micron
  INTEGER(KIND=INT4), PRIVATE, PARAMETER :: Chan_Idx_133um = 33  !channel number for 13.3 micron
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: Beta_Ap_Water = 1.3
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: Beta_Ap_Ice = 1.06
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: Ice_Temperature_Min = 243.0
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: Ice_Temperature_Max = 263.0

  REAL(KIND=REAL4), PARAMETER, PRIVATE:: Rad_Window_Thresh = 2.0
  REAL(KIND=REAL4), PARAMETER, PRIVATE:: Rad_H2O_Thresh = 0.25
  REAL, PRIVATE, PARAMETER:: Bt_Ch27_Ch31_Covar_Cirrus_Thresh = 0.5
  REAL, PRIVATE, PARAMETER:: Bt_Ch31_Btd_Ch31_Ch32_Covar_Cirrus_Thresh = -1.0
  REAL, PRIVATE, PARAMETER:: Beta_11um_85um_Ice_Thresh = 1.10
  REAL, PRIVATE, PARAMETER:: Beta_11um_12um_Overlap_Thresh = 0.95
  REAL, PRIVATE, PARAMETER:: Beta_11um_133um_Overlap_Thresh = 0.70
  REAL, PRIVATE, PARAMETER:: Emiss_67um_Cirrus_Thresh = 0.05
  REAL, PRIVATE, PARAMETER:: Emiss_133um_Cirrus_Thresh = 0.02

  CONTAINS

!----------------------------------------------------------------------
! Cloud Height for AVHRR/1.  Use precomputed opaque solutions
! for AVHRR/1 sensors.  Set quality flag to degraded (1)
!----------------------------------------------------------------------
subroutine  OPAQUE_CLOUD_HEIGHT(Line_Idx_min,Num_Lines)
                              
  INTEGER, INTENT(IN):: Line_Idx_min
  INTEGER, INTENT(IN):: Num_Lines

  INTEGER:: Elem_Idx
  INTEGER:: Line_Idx
  INTEGER:: Number_Of_Elements

  !--- intialize local variables using global variables
  Number_Of_Elements = num_pix

  !--- initialize output
  Tc_Acha =  Missing_Value_Real4
  Ec_Acha =  Missing_Value_Real4
  Beta_Acha =  Missing_Value_Real4
  Pc_Acha =  Missing_Value_Real4
  Zc_Acha =  Missing_Value_Real4
  Acha_Quality_Flag = 0
  Acha_OE_Quality_Flags = 0
  Cld_Layer_Acha = 0

  !--------------------------------------------------------------------------
  ! loop over pixels in scanlines
  !--------------------------------------------------------------------------
  Line_loop: do Line_Idx = Line_Idx_min,Num_Lines + Line_Idx_min - 1

    Element_loop:   do Elem_Idx = 1, Number_Of_Elements

    !--- check for a bad pixel
    if (Bad_pixel_mask(Elem_Idx,Line_Idx) == sym%YES) then
          cycle
    endif

    if (Tc_Opaque_Cloud(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
       Pc_Acha(Elem_Idx,Line_Idx) = Pc_Opaque_Cloud(Elem_Idx,Line_Idx)
       Tc_Acha(Elem_Idx,Line_Idx) = Tc_Opaque_Cloud(Elem_Idx,Line_Idx)
       Zc_Acha(Elem_Idx,Line_Idx) = Zc_Opaque_Cloud(Elem_Idx,Line_Idx)
       Ec_Acha(Elem_Idx,Line_Idx) = 1.0
       Acha_Quality_Flag(Elem_Idx,Line_Idx) = 1

       !------- determine cloud layer based on pressure
       if (Pc_Acha(Elem_Idx,Line_Idx) <= 440.0) then
         Cld_Layer_Acha(Elem_Idx,Line_Idx) = 3
       elseif (Pc_Acha(Elem_Idx,Line_Idx) < 680.0) then
         Cld_Layer_Acha(Elem_Idx,Line_Idx) = 2
       else
         Cld_Layer_Acha(Elem_Idx,Line_Idx) = 1
       endif

       !----- set beta passed on temp
       if (Tc_Acha(Elem_Idx,Line_Idx) < 260.0) then
           Beta_Acha(Elem_Idx,Line_Idx) = Beta_Ap_Ice
       else
           Beta_Acha(Elem_Idx,Line_Idx) = Beta_Ap_Water
       endif

    endif

    end do Element_loop

  end do Line_loop
end subroutine OPAQUE_CLOUD_HEIGHT

!----------------------------------------------------------------------
! routine to interpolate nwp winds at cloud height
!----------------------------------------------------------------------
subroutine COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND(Line_Idx_Min,Num_Lines)

  INTEGER, INTENT(IN):: Line_Idx_Min
  INTEGER, INTENT(IN):: Num_Lines
  INTEGER:: Num_Elem
  INTEGER:: Elem_Idx
  INTEGER:: Line_Idx
  INTEGER:: Line_Start
  INTEGER:: Line_End
  INTEGER:: Inwp
  INTEGER:: Jnwp
  INTEGER:: Ivza
  INTEGER:: Level
  REAL:: Tc_Temp
  REAL:: Zc_Temp
  REAL:: U_Wnd_Temp
  REAL:: V_Wnd_Temp

  !--- save elements
  Line_Start = Line_Idx_Min 
  Line_End = Line_Start + Num_Lines - 1
  Num_Elem = Num_Pix      !Num_Pix is a global variable

  !--- initialize
  Wnd_Spd_Cld_Top_Nwp_Pix = Missing_Value_Real4
  Wnd_Dir_Cld_Top_Nwp_Pix = Missing_Value_Real4

  !----------------------------------------------------------
  ! loop through segment
  !----------------------------------------------------------
  Line_Loop_1: do Line_Idx = Line_Start, Line_End
  Element_Loop_1: do Elem_Idx = 1, Num_Elem

     !--- save indices
     Inwp = I_Nwp(Elem_Idx,Line_Idx)
     Jnwp = J_Nwp(Elem_Idx,Line_Idx)
     Ivza = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !-- check if indices are valid
     if (Inwp < 0 .or. Jnwp < 0) cycle

     !-- check if cld-top pressure is valid
     if (Pc_Acha(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- determine Level in profiles 
     !--- note, wind profiles are at the native nwp resolution
     call KNOWING_P_COMPUTE_T_Z_NWP(Inwp,Jnwp, &
                                    Pc_Acha(Elem_Idx,Line_Idx), &
                                    Tc_Temp,Zc_Temp,Level)

     !--- assign u and v winds
     U_Wnd_Temp = U_Wnd_Prof_Nwp(Level,Inwp,Jnwp) 
     V_Wnd_Temp = V_Wnd_Prof_Nwp(Level,Inwp,Jnwp) 

     !-- determine speed and direction
     Wnd_Spd_Cld_Top_Nwp_Pix(Elem_Idx,Line_Idx) = wind_speed(U_Wnd_Temp, V_Wnd_Temp)
     Wnd_Dir_Cld_Top_Nwp_Pix(Elem_Idx,Line_Idx) = wind_direction(U_Wnd_Temp, V_Wnd_Temp)

  END DO Element_Loop_1
  END DO Line_Loop_1

end subroutine COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND

!----------------------------------------------------------------------
! routine to interpolate pressure to flight level altitude.
!----------------------------------------------------------------------
subroutine COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min,Num_Lines)

!--- Based on Sarah Monette calculations for HS3.  Her calculation is based on:
!--- http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf.

  INTEGER, INTENT(IN):: Line_Idx_Min
  INTEGER, INTENT(IN):: Num_Lines
  INTEGER:: Num_Elem
  INTEGER:: Elem_Idx
  INTEGER:: Line_Idx
  INTEGER:: Line_Start
  INTEGER:: Line_End
  INTEGER:: Inwp
  INTEGER:: Jnwp
  REAL:: Pc_Temp
  REAL:: Tc_Temp
  REAL:: Alt_Temp

  !--- Constants from Sarah Monette
  REAL, PARAMETER :: PW1 = 227.9 ! hPa
  REAL, PARAMETER :: PW2 = 56.89 ! hPa
  REAL, PARAMETER :: PW3 = 11.01 ! hPa
  REAL, PARAMETER :: P0 = 1013.25 ! hPa
  REAL, PARAMETER :: LR_OVER_G = 0.190263
  REAL, PARAMETER :: Z0 = 145422.16 ! feet
  REAL, PARAMETER :: LN_1 = -20859.0
  REAL, PARAMETER :: LN_2 = 149255.0
  REAL, PARAMETER :: PN_4 = 0.000470034
  REAL, PARAMETER :: PN_3 = -0.364267
  REAL, PARAMETER :: PN_2 = 47.5627
  REAL, PARAMETER :: PN_1 = -2647.45
  REAL, PARAMETER :: PN_0 = 123842.0

  !--- save elements
  Line_Start = Line_Idx_Min
  Line_End = Line_Start + Num_Lines - 1
  Num_Elem = Num_Pix      !Num_Pix is a global variable

  !--- initialize
  Alt_Acha = Missing_Value_Real4

  !----------------------------------------------------------
  ! loop through segment
  !----------------------------------------------------------
  Line_Loop_1: do Line_Idx = Line_Start, Line_End
  Element_Loop_1: do Elem_Idx = 1, Num_Elem

     !--- Initialize temporary value each time.
     Pc_Temp = Missing_Value_Real4
     Tc_Temp = Missing_Value_Real4
     Alt_Temp = Missing_Value_Real4

     !--- save indices
     Inwp = I_Nwp(Elem_Idx,Line_Idx)
     Jnwp = J_Nwp(Elem_Idx,Line_Idx)

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !--- check if indices are valid
     !--- stw May not need this.
     if (Inwp < 0 .or. Jnwp < 0) cycle

     !--- check if cld-top pressure is valid
     if (Pc_Acha(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- Place valid pressure in temp variable for readability in the
     !--- calculations below.
     Pc_Temp = Pc_Acha(Elem_Idx,Line_Idx)
     Tc_Temp = Tc_Acha(Elem_Idx,Line_Idx)

     !IF (Pc_Temp < PW3) THEN
     !  print*,"Pc_Temp : ", Pc_Temp, Tc_Temp
     !ENDIF

     !--- calculated altitude, in feet, from pressure.
     !--- 1st pivot point is directly from the pressure to
     !--- altitude from above reference.
     IF (Pc_Temp > PW1) THEN
       Alt_Temp = (1.0 - (Pc_Temp/P0)**LR_OVER_G) * Z0
!stw       print*,"Pc_Temp/Alt_Temp : ", Pc_Temp, Alt_Temp
     ENDIF

     !--- 2nd pivot point was modeled best with a natural log
     !--- fit.  From Sarah Monette.
     IF (Pc_Temp <= PW1 .AND. Pc_Temp >= PW2) THEN
       Alt_Temp = LN_1 * LOG(Pc_Temp) + LN_2
!stw       print*,"Pc_Temp/Alt_Temp : ", Pc_Temp, Alt_Temp
     ENDIF

     !--- 3rd pivot point. Modeled best with a polynomial
     !--- fit from Sarah Monette.
     IF (Pc_Temp < PW2 .AND. Pc_Temp >= PW3) THEN
       Alt_Temp = (PN_4*Pc_Temp**4) + (PN_3*Pc_Temp**3) + (PN_2*Pc_Temp**2) + &
                  (PN_1*Pc_Temp) + PN_0
!stw       print*,"Pc_Temp/Alt_Temp : ", Pc_Temp, Alt_Temp
     ENDIF

     IF (Pc_Temp < PW3) THEN
       Alt_Temp = Missing_Value_Real4
!stw       print*,"Pc_Temp/Alt_Temp : ", Pc_Temp, Alt_Temp
     ENDIF

     !--- Assign final altitude, in feet, to the level2 array.
     Alt_Acha(Elem_Idx,Line_Idx) = Alt_Temp

  END DO Element_Loop_1
  END DO Line_Loop_1

end subroutine COMPUTE_ALTITUDE_FROM_PRESSURE

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module CLOUD_HEIGHT_ROUTINES
