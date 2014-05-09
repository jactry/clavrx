!$Id$
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

  PRIVATE::  OPAQUE_CLOUD_HEIGHT   
  PRIVATE::  H2O_CLOUD_HEIGHT   
  PUBLIC::  COMPUTE_OPAQUE_CLOUD_HEIGHT   
  PUBLIC::  COMPUTE_H2O_CLOUD_HEIGHT   
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
  real(kind=real4), parameter, private:: Rad_11um_Thresh = 2.0
  real(kind=real4), parameter, private:: Rad_H2O_Thresh = 0.25
  real, private, parameter:: Bt_Ch27_Ch31_Covar_Cirrus_Moist_Thresh = 0.25
  real, private, parameter:: Bt_Ch27_Ch31_Covar_Cirrus_Thresh = 1.0 !0.5

  CONTAINS

!----------------------------------------------------------------------
! Compute opaque cloud height
!----------------------------------------------------------------------
subroutine  COMPUTE_OPAQUE_CLOUD_HEIGHT(Line_Idx_min,Num_Lines)
                              
  INTEGER, INTENT(IN):: Line_Idx_min
  INTEGER, INTENT(IN):: Num_Lines

  INTEGER:: Elem_Idx
  INTEGER:: Line_Idx
  INTEGER:: Number_Of_Elements
  INTEGER:: Nwp_Lon_Idx
  INTEGER:: Nwp_Lat_Idx
  INTEGER:: Vza_Idx

  !--- intialize local variables using global variables
  Number_Of_Elements = num_pix

  !--- initialize output
  Tc_Opaque_Cloud = Missing_Value_Real4
  Pc_Opaque_Cloud = Missing_Value_Real4
  Zc_Opaque_Cloud = Missing_Value_Real4

  !--------------------------------------------------------------------------
  ! loop over pixels in scanlines
  !--------------------------------------------------------------------------
  Line_loop: do Line_Idx = Line_Idx_min,Num_Lines + Line_Idx_min - 1
    Element_loop:   do Elem_Idx = 1, Number_Of_Elements

    !--- check for a bad pixel
    if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
          cycle
    endif

    !--- save indices
    Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
    Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
    Vza_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

    if (Nwp_Lon_Idx <= 0 .or. Nwp_Lat_Idx <= 0) cycle

    if (Chan_On_Flag_Default(31) == sym%YES) then

          call OPAQUE_CLOUD_HEIGHT(ch(31)%Rad_Toa(Elem_Idx,Line_Idx), &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(31)%Rad_BB_Cloud_Profile, &
                    P_Std_Rtm, &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Z_prof, &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%T_prof, &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level, &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level, &
                    Tc_Opaque_Cloud(Elem_Idx,Line_Idx),  &
                    Zc_Opaque_Cloud(Elem_Idx,Line_Idx),  &
                    Pc_Opaque_Cloud(Elem_Idx,Line_Idx))
    endif    

    end do Element_loop
  end do Line_loop
end subroutine COMPUTE_OPAQUE_CLOUD_HEIGHT

!----------------------------------------------------------------------
! Compute h2o cloud height
!----------------------------------------------------------------------
subroutine  COMPUTE_H2O_CLOUD_HEIGHT(Line_Idx_min,Num_Lines)
                              
  INTEGER, INTENT(IN):: Line_Idx_min
  INTEGER, INTENT(IN):: Num_Lines

  INTEGER:: Elem_Idx
  INTEGER:: Line_Idx
  INTEGER:: Number_Of_Elements
  INTEGER:: Nwp_Lon_Idx
  INTEGER:: Nwp_Lat_Idx
  INTEGER:: Vza_Idx

  !--- intialize local variables using global variables
  Number_Of_Elements = num_pix

  !--- initialize output
  Tc_H2O = Missing_Value_Real4
  Pc_H2O = Missing_Value_Real4
  Zc_H2O = Missing_Value_Real4

  !--------------------------------------------------------------------------
  ! loop over pixels in scanlines
  !--------------------------------------------------------------------------
  Line_loop: do Line_Idx = Line_Idx_min,Num_Lines + Line_Idx_min - 1
    Element_loop:   do Elem_Idx = 1, Number_Of_Elements

    !--- check for a bad pixel
    if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
          cycle
    endif

    !--- save indices
    Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
    Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
    Vza_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

    if (Nwp_Lon_Idx <= 0 .or. Nwp_Lat_Idx <= 0) cycle


    if (Chan_On_Flag_Default(27) == sym%YES .and. &
        Chan_On_Flag_Default(31) == sym%YES) then

          call H2O_CLOUD_HEIGHT(ch(31)%Rad_Toa(Elem_Idx,Line_Idx), &
                    ch(31)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(31)%Rad_BB_Cloud_Profile, &
                    ch(27)%Rad_Toa(Elem_Idx,Line_Idx), &
                    ch(27)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Idx)%ch(27)%Rad_BB_Cloud_Profile, &
                    P_Std_Rtm, &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Z_prof, &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%T_prof, &
                    Covar_Ch27_Ch31_5x5(Elem_Idx,Line_Idx), &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level, &
                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level, &
                    Tc_H2o(Elem_Idx,Line_Idx),  &
                    Zc_H2o(Elem_Idx,Line_Idx),  &
                    Pc_H2o(Elem_Idx,Line_Idx))
    endif

    end do Element_loop
  end do Line_loop

end subroutine COMPUTE_H2O_CLOUD_HEIGHT
!----------------------------------------------------------------------
! Compute the ACHA values without calling ACHA (mode = 0)
! 
! This needs to have the opaque and h2o cloud heights generated
!----------------------------------------------------------------------
subroutine  MODE_ZERO_CLOUD_HEIGHT(Line_Idx_min,Num_Lines)

  INTEGER, INTENT(IN):: Line_Idx_min
  INTEGER, INTENT(IN):: Num_Lines

  INTEGER:: Elem_Idx
  INTEGER:: Line_Idx
  INTEGER:: Number_Of_Elements
  INTEGER:: Nwp_Lon_Idx
  INTEGER:: Nwp_Lat_Idx
  INTEGER:: Vza_Idx

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
    if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
          cycle
    endif

    !--- combine heights into a coherent product
    if (Tc_H2O(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
       Pc_Acha(Elem_Idx,Line_Idx) = Pc_H2O(Elem_Idx,Line_Idx)
       Tc_Acha(Elem_Idx,Line_Idx) = Tc_H2O(Elem_Idx,Line_Idx)
       Zc_Acha(Elem_Idx,Line_Idx) = Zc_H2O(Elem_Idx,Line_Idx)
       Ec_Acha(Elem_Idx,Line_Idx) = 1.0
       Acha_Quality_Flag(Elem_Idx,Line_Idx) = 1

    elseif (Tc_Opaque_Cloud(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
       Pc_Acha(Elem_Idx,Line_Idx) = Pc_Opaque_Cloud(Elem_Idx,Line_Idx)
       Tc_Acha(Elem_Idx,Line_Idx) = Tc_Opaque_Cloud(Elem_Idx,Line_Idx)
       Zc_Acha(Elem_Idx,Line_Idx) = Zc_Opaque_Cloud(Elem_Idx,Line_Idx)
       Ec_Acha(Elem_Idx,Line_Idx) = 1.0
       Acha_Quality_Flag(Elem_Idx,Line_Idx) = 1
    else
       cycle
    endif

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

    end do Element_loop

  end do Line_loop
end subroutine MODE_ZERO_CLOUD_HEIGHT

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
  INTEGER:: Nwp_Lon_Idx
  INTEGER:: Nwp_Lat_Idx
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
     Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
     Ivza = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !-- check if indices are valid
     if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0) cycle

     !-- check if cld-top pressure is valid
     if (Pc_Acha(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- determine Level in profiles 
     !--- note, wind profiles are at the native nwp resolution
     call KNOWING_P_COMPUTE_T_Z_NWP(Nwp_Lon_Idx,Nwp_Lat_Idx, &
                                    Pc_Acha(Elem_Idx,Line_Idx), &
                                    Tc_Temp,Zc_Temp,Level)

     !--- assign u and v winds
     U_Wnd_Temp = U_Wnd_Prof_Nwp(Level,Nwp_Lon_Idx,Nwp_Lat_Idx) 
     V_Wnd_Temp = V_Wnd_Prof_Nwp(Level,Nwp_Lon_Idx,Nwp_Lat_Idx) 

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
  INTEGER:: Nwp_Lon_Idx
  INTEGER:: Nwp_Lat_Idx
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
     Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !--- check if indices are valid
     !--- stw May not need this.
     if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0) cycle

     !--- check if cld-top pressure is valid
     if (Pc_Acha(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- Place valid pressure in temp variable for readability in the
     !--- calculations below.
     Pc_Temp = Pc_Acha(Elem_Idx,Line_Idx)
     Tc_Temp = Tc_Acha(Elem_Idx,Line_Idx)

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

!====================================================================
! Function Name: OPAQUE_CLOUD_HEIGHT
!
! Function: estimate the cloud temperature/height/pressure
!
! Description: Use the 11um obs and assume the cloud is back and 
!           estimate height from 11 um BB cloud profile
!              
! Dependencies: 
!
! Restrictions: 
!
! Reference: 
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
!====================================================================
 subroutine OPAQUE_CLOUD_HEIGHT(Ch31_Rad, &
                                Ch31_Rad_BB_Prof, &
                                P_Prof, &
                                Z_Prof, &
                                T_Prof, &
                                Tropo_Level,  &
                                Sfc_Level, &
                                Tc_Opa, &
                                Zc_Opa, &
                                Pc_Opa)
   real, intent(in):: Ch31_Rad
   real, intent(in), dimension(:):: Ch31_Rad_BB_Prof
   real, intent(in), dimension(:):: P_Prof
   real, intent(in), dimension(:):: Z_Prof
   real, intent(in), dimension(:):: T_Prof
   integer (kind=int1), intent(in):: Tropo_Level
   integer (kind=int1), intent(in):: Sfc_Level
   real, intent(out):: Pc_Opa
   real, intent(out):: Zc_Opa
   real, intent(out):: Tc_Opa

   integer:: Level_Idx
   integer:: Level_Idx_Start
   integer:: Level_Idx_End
   integer:: Level_Idx_Max_Valid_Cloud

   !--- initialize
   Pc_Opa =  Missing_Value_Real4
   Zc_Opa =  Missing_Value_Real4
   Tc_Opa =  Missing_Value_Real4

   !--- restrict levels to consider
   Level_Idx_Start = Tropo_Level 
   Level_Idx_End = Sfc_Level 

   !--- initialize levels
   Level_Idx_Max_Valid_Cloud = 0

   !--- check for stratospheric
   if (Ch31_Rad < Ch31_Rad_BB_Prof(Level_Idx_start)) then

        Level_Idx_Max_Valid_Cloud = Level_Idx_Start

   else

       level_loop: do Level_Idx = Level_Idx_Start, Level_Idx_End

         if (Ch31_Rad > Ch31_Rad_BB_Prof(Level_Idx)) then

           Level_Idx_Max_Valid_Cloud = Level_Idx

         endif

        end do level_loop

   endif

   !--- compute lowest pressure level with valid 11 micron emissivity
   if (Level_Idx_max_valid_Cloud > 0) then
        if (Level_Idx_Max_Valid_Cloud > 1) then
         Pc_Opa = P_prof(Level_Idx_Max_Valid_Cloud)
         Zc_Opa = Z_prof(Level_Idx_Max_Valid_Cloud)
         Tc_Opa = T_prof(Level_Idx_Max_Valid_Cloud)
        endif
   endif

 end subroutine OPAQUE_CLOUD_HEIGHT

!====================================================================
! Function Name: H2O_CLOUD_HEIGHT
!
! Function: estimate the cloud temperature/height/pressure
!
! Description: Use the 11um and 6.7um obs and the RTM cloud BB profiles
!              to perform h2o intercept on a pixel level. Filters
!              restrict this to high clouds only
!              
! Dependencies: 
!
! Restrictions: 
!
! Reference: 
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
!====================================================================
subroutine  H2O_CLOUD_HEIGHT(Rad_11um, &
                             Rad_11um_Clear, &
                             Rad_11um_BB_Profile, &
                             Rad_H2O, &
                             Rad_H2O_Clear,  &
                             Rad_H2O_BB_Profile, &
                             P_Prof, &
                             Z_Prof, &
                             T_Prof, &
                             Covar_H2O_Window, &
                             Tropo_Level, &
                             Sfc_Level, &
                             Tc,  &
                             Zc,  &
                             Pc)

  real, intent(in):: Rad_11um
  real, intent(in):: Rad_H2O
  real, intent(in), dimension(:):: Rad_11um_BB_Profile
  real, intent(in), dimension(:):: Rad_H2O_BB_Profile
  real, intent(in), dimension(:):: P_Prof
  real, intent(in), dimension(:):: T_Prof
  real, intent(in), dimension(:):: Z_Prof
  real, intent(in):: Rad_11um_Clear
  real, intent(in):: Rad_H2O_Clear
  real, intent(in):: Covar_H2O_Window
  integer (kind=int1), intent(in):: Sfc_Level
  integer (kind=int1), intent(in):: Tropo_Level
  real (kind=real4), intent(out) :: Tc
  real (kind=real4), intent(out) :: Pc
  real (kind=real4), intent(out) :: Zc

  real:: Rad_H2O_BB_Prediction
  real:: Slope
  real:: Intercept
  real:: Denominator
  integer:: ilev
  integer:: ilev_h2o

  !--- initialize
  Pc = Missing_Value_Real4
  Tc = Missing_Value_Real4
  Zc = Missing_Value_Real4

  !--- determine if a solution should be attempted
  if (Rad_11um_Clear - Rad_11um < Rad_11um_Thresh) then
      return 
  endif

  if (Rad_H2O_Clear - Rad_H2O < Rad_H2O_Thresh) then
      return
  endif

  if (Covar_H2O_Window /= Missing_Value_Real4 .and. Covar_H2O_Window < Bt_Ch27_Ch31_Covar_Cirrus_Thresh) then
      return 
  endif

 !--- attempt a solution

 !--- colder than tropo
 if (Rad_11um < Rad_11um_BB_Profile(Tropo_Level)) then

     ilev_h2o = Tropo_Level

 else   !if not, attempt solution

     !--- determine linear regress of h2o (y)  as a function of window (x)
      Denominator =  Rad_11um - Rad_11um_Clear

      if (Denominator < 0.0) then
             Slope = (Rad_H2O - Rad_H2O_Clear) / (Denominator)
             Intercept = Rad_H2O - Slope*Rad_11um
      else
            return 
      endif

      !--- brute force solution
      ilev_h2o = 0

      do ilev = Tropo_Level+1, Sfc_Level
          Rad_H2O_BB_Prediction = Slope*Rad_11um_BB_Profile(ilev) + Intercept

          if (Rad_H2O_BB_Prediction < 0) cycle

          if ((Rad_H2O_BB_Prediction > Rad_H2O_BB_Profile(ilev-1)) .and. & 
               (Rad_H2O_BB_Prediction <= Rad_H2O_BB_Profile(ilev))) then
               ilev_h2o = ilev
               exit
          endif

      enddo

 endif    !tropopause check

 !--- adjust back to full Rtm profile indices
 if (ilev_h2o > 0) then
       Pc = P_Prof(ilev_h2o)
       Tc = T_Prof(ilev_h2o)
       Zc = Z_Prof(ilev_h2o)
  endif

end subroutine H2O_CLOUD_HEIGHT

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module CLOUD_HEIGHT_ROUTINES
