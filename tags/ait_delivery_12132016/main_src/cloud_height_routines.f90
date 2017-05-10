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
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PRIVATE
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PRIVATE USE. THEY ARE
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

  public::  COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND_AND_TPW
  public::  COMPUTE_ALTITUDE_FROM_PRESSURE
  public::  OPAQUE_CLOUD_HEIGHT
  public::  H2O_CLOUD_HEIGHT
  public::  CO2IRW_CLOUD_HEIGHT
  public::  CO2_SLICING_CLOUD_HEIGHT
  public::  CO2_SLICING_CLOUD_HEIGHT_NEW
  public::  OPAQUE_TRANSMISSION_HEIGHT
  public::  CONVECTIVE_CLOUD_PROBABILITY
  public::  SUPERCOOLED_CLOUD_PROBABILITY
  public::  MODIFY_CLOUD_TYPE_WITH_SOUNDER
  public::  COMPUTE_SOUNDER_MASK_SEGMENT
  public::  SOUNDER_EMISSIVITY
  public::  INTERCEPT_CLOUD_PRESSURE
  public::  SLICING_CLOUD_PRESSURE
  public::  CTP_MULTILAYER

  !--- include parameters for each system here
  integer(kind=int4), private, parameter :: Chan_Idx_375um = 20  !channel number for 3.75 micron
  integer(kind=int4), private, parameter :: Chan_Idx_67um = 27   !channel number for 6.7 micron
  integer(kind=int4), private, parameter :: Chan_Idx_73um = 28   !channel number for 7.3 micron
  integer(kind=int4), private, parameter :: Chan_Idx_85um = 29   !channel number for 8.5 micron
  integer(kind=int4), private, parameter :: Chan_Idx_11um = 31   !channel number for 11 micron
  integer(kind=int4), private, parameter :: Chan_Idx_12um = 32   !channel number for 12 micron
  integer(kind=int4), private, parameter :: Chan_Idx_133um = 33  !channel number for 13.3 micron
  real(kind=real4), parameter, private:: Beta_Ap_Water = 1.3
  real(kind=real4), parameter, private:: Beta_Ap_Ice = 1.06
  real(kind=real4), parameter, private:: Ice_Temperature_Min = 243.0
  real(kind=real4), parameter, private:: Ice_Temperature_Max = 263.0
  real(kind=real4), parameter, private:: Rad_11um_Thresh = 2.0
  real(kind=real4), parameter, private:: Rad_H2O_Thresh = 0.25
  real, private, parameter:: Bt_Ch27_Ch31_Covar_Cirrus_Moist_Thresh = 0.25
  real, private, parameter:: Bt_Ch27_Ch31_Covar_Cirrus_Thresh = 1.0 !0.5


  integer, private, parameter:: N_Sc_Lut = 20
  real(kind=real4), dimension(N_Sc_Lut), private, parameter:: &
  Sc_Tc_Lut = (/202.21,206.77,211.33,215.88,220.44,225.00,229.56,234.11,238.67,243.23, &
                247.79,252.34,256.90,261.46,266.01,270.57,275.13,279.69,284.24,288.80/)

  real(kind=real4), dimension(N_Sc_Lut), private, parameter:: &
  Sc_Prob_Lut = (/0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.002,0.051,0.243, &
                  0.470,0.658,0.800,0.886,0.888,0.870,0.142,0.004,0.004,0.004 /)

  contains

!----------------------------------------------------------------------
! Compute the ACHA values without calling ACHA (mode = 0)
!----------------------------------------------------------------------
subroutine CONVECTIVE_CLOUD_PROBABILITY(Bad_Pixel_Mask,Bt11,Bt67,Emiss_11_Tropo,Tsfc,Conv_Cld_Prob)
  integer(kind=int1), dimension(:,:), intent(in):: Bad_Pixel_Mask
  real, dimension(:,:), intent(in):: Bt11, Bt67, Emiss_11_Tropo,Tsfc
  real(kind=real4), dimension(:,:), intent(out):: Conv_Cld_Prob

  real, parameter:: Btd_Thresh = -2.0
  real, parameter:: Etrop_Thresh_1 = 0.95
  real, parameter:: Etrop_Thresh_2 = 0.90
  real, parameter:: Tsfc_Thresh = 30.0

  !--- initialize to 0 (no)
  Conv_Cld_Prob = 0.0

  !--- set bad data to missing
  where(Bad_Pixel_Mask ==sym%YES)
    Conv_Cld_Prob = Missing_Value_Real4
  endwhere

  !--- if only 11 micron, use a tight threshold
  if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
    where(Emiss_11_Tropo >= Etrop_Thresh_1) 
      Conv_Cld_Prob = 1.0
    endwhere
  endif

  if (Sensor%Chan_On_Flag_Default(27) == sym%YES .and. &
      Sensor%Chan_On_Flag_Default(31) == sym%YES) then
    where(Emiss_11_Tropo >= Etrop_Thresh_2 .and. (Bt67 - Bt11) >= Etrop_Thresh_2) 
     Conv_Cld_Prob = 1.0
    endwhere
  endif

  !--- limit false alarms over elevated terrain
  if (Sensor%Chan_On_Flag_Default(31) == sym%YES) then
    where(Tsfc - Bt11 < Tsfc_Thresh)
       Conv_Cld_Prob = 0.0
    endwhere
  endif

end subroutine CONVECTIVE_CLOUD_PROBABILITY
!--------------------------------------------------------------------------------------------------
! Supercooled Cloud Probability
!--------------------------------------------------------------------------------------------------
subroutine SUPERCOOLED_CLOUD_PROBABILITY(Bad_Pixel_Mask,Cloud_Type,Cloud_Temperature,Supercooled_Cld_Prob)
  integer(kind=int1), dimension(:,:), intent(in):: Bad_Pixel_Mask
  integer(kind=int1), dimension(:,:), intent(in):: Cloud_Type
  real(kind=real4), dimension(:,:), intent(in):: Cloud_Temperature
  real(kind=real4), dimension(:,:), intent(out):: Supercooled_Cld_Prob

  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Number_Of_Elements
  integer:: Number_Of_Lines
  integer:: Tc_Idx

  !--- intialize local variables using global variables
  Number_Of_Elements = Image%Number_Of_Elements
  Number_Of_Lines = Image%Number_Of_Lines_Read_This_Segment

  !--- initialize to Missing
  Supercooled_Cld_Prob = Missing_Value_Real4

  !--- Loop through each pixel
  do Elem_Idx = 1, Number_Of_Elements 
     do Line_Idx = 1, Number_Of_Lines
        
       !--- filter bad
       if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
       !--- filter missing temps
       if (Cloud_Temperature(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
       !--- filter with cloud type
       if (Cloud_Type(Elem_Idx,Line_Idx) == Missing_Value_Int1) cycle

       if (Cloud_Type(Elem_Idx,Line_Idx) /= sym%FOG_TYPE .and. &
           Cloud_Type(Elem_Idx,Line_Idx) /= sym%WATER_TYPE .and. &
           Cloud_Type(Elem_Idx,Line_Idx) /= sym%SUPERCOOLED_TYPE) then
           Supercooled_Cld_Prob(Elem_Idx,Line_Idx) = 0.0
           cycle
       endif

       Tc_Idx = minloc(abs(Cloud_Temperature(Elem_Idx,Line_Idx) - Sc_Tc_Lut),1) 

       Tc_Idx = max(1,min(N_Sc_Lut,Tc_Idx))

       Supercooled_Cld_Prob(Elem_Idx,Line_Idx) = Sc_Prob_Lut(Tc_Idx) 

     enddo
  enddo

end subroutine SUPERCOOLED_CLOUD_PROBABILITY

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine MODIFY_CLOUD_TYPE_WITH_SOUNDER (Tc_CO2, Ec_CO2, Cloud_Type)
   real(kind=real4), intent(in), dimension(:,:):: Tc_CO2
   real(kind=real4), intent(in), dimension(:,:):: Ec_CO2
   integer(kind=int1), intent(inout), dimension(:,:):: Cloud_Type

   real(kind=real4), parameter:: Ec_Thresh = 0.3
   real(kind=real4), parameter:: Tc_Thresh = 240.0 !K
   integer:: Elem_Idx,Line_Idx,Num_Elem, Num_Lines
   integer:: Elem_Lrc_Idx,Line_Lrc_Idx

   Num_Elem = Image%Number_Of_Elements
   Num_Lines = Image%Number_Of_Lines_Read_This_Segment

   where((Cloud_Type == sym%FOG_TYPE .or.  &
         Cloud_Type == sym%WATER_TYPE .or. &
         Cloud_Type == sym%SUPERCOOLED_TYPE) .and. &
         Ec_CO2 > Ec_Thresh .and. &
         Tc_CO2 /= Missing_Value_Real4 .and. &
         Tc_CO2 < Tc_Thresh)

     Cloud_Type = 99

   endwhere

   do Line_Idx = 1, Num_Lines
      do Elem_Idx = 1, Num_Elem
         Elem_Lrc_Idx = i_lrc(Elem_Idx,Line_Idx)
         Line_Lrc_Idx = j_lrc(Elem_Idx,Line_Idx)
         if (Elem_Lrc_Idx > 0 .and. Line_Lrc_Idx > 0) then
            if (Cloud_Type(Elem_Lrc_Idx,Line_Lrc_Idx) == 99) then
               Cloud_Type(Elem_Idx,Line_Idx) = 99
            endif
         endif
      enddo
   enddo

   where(Cloud_Type == 99)
    Cloud_Type = sym%OVERLAP_TYPE
   endwhere


end subroutine MODIFY_CLOUD_TYPE_WITH_SOUNDER

!----------------------------------------------------------------------
! Compute the ACHA values without calling ACHA (mode = 0)
! 
! This needs to have the opaque and h2o cloud heights generated
!----------------------------------------------------------------------
subroutine  MODE_ZERO_CLOUD_HEIGHT(Line_Idx_min,Num_Lines)

  integer, intent(in):: Line_Idx_min
  integer, intent(in):: Num_Lines

  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Number_Of_Elements

  !--- intialize local variables using global variables
  Number_Of_Elements = Image%Number_Of_Elements

  !--- initialize output
  ACHA%Tc =  Missing_Value_Real4
  ACHA%Ec =  Missing_Value_Real4
  ACHA%Beta =  Missing_Value_Real4
  ACHA%Pc =  Missing_Value_Real4
  ACHA%Zc =  Missing_Value_Real4
  ACHA%Quality_Flag = 0
  ACHA%OE_Quality_Flags = 0

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
       ACHA%Pc(Elem_Idx,Line_Idx) = Pc_H2O(Elem_Idx,Line_Idx)
       ACHA%Tc(Elem_Idx,Line_Idx) = Tc_H2O(Elem_Idx,Line_Idx)
       ACHA%Zc(Elem_Idx,Line_Idx) = Zc_H2O(Elem_Idx,Line_Idx)
       ACHA%Ec(Elem_Idx,Line_Idx) = 1.0
       ACHA%Quality_Flag(Elem_Idx,Line_Idx) = 1

    elseif (Tc_Opaque_Cloud(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
       ACHA%Pc(Elem_Idx,Line_Idx) = Pc_Opaque_Cloud(Elem_Idx,Line_Idx)
       ACHA%Tc(Elem_Idx,Line_Idx) = Tc_Opaque_Cloud(Elem_Idx,Line_Idx)
       ACHA%Zc(Elem_Idx,Line_Idx) = Zc_Opaque_Cloud(Elem_Idx,Line_Idx)
       ACHA%Ec(Elem_Idx,Line_Idx) = 1.0
       ACHA%Quality_Flag(Elem_Idx,Line_Idx) = 1
    else
       cycle
    endif

    !----- set beta passed on temp
    if (ACHA%Tc(Elem_Idx,Line_Idx) < 260.0) then
        ACHA%Beta(Elem_Idx,Line_Idx) = Beta_Ap_Ice
    else
        ACHA%Beta(Elem_Idx,Line_Idx) = Beta_Ap_Water
    endif

    end do Element_loop

  end do Line_loop
end subroutine MODE_ZERO_CLOUD_HEIGHT

!----------------------------------------------------------------------
! routine to interpolate nwp winds at cloud height
!----------------------------------------------------------------------
subroutine COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND_AND_TPW(Line_Idx_Min,Num_Lines)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  integer:: Num_Elem
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Line_Start
  integer:: Line_End
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  integer:: Nwp_Lon_Idx_x
  integer:: Nwp_Lat_Idx_x
  integer:: Ivza
  integer:: Level
  real:: Tc_Temp
  real:: Zc_Temp
  real:: U_Wnd_Temp
  real:: V_Wnd_Temp
  real:: dp
  real:: dtpw
  real:: Tpw_1
  real:: Tpw_2

  !--- save elements
  Line_Start = Line_Idx_Min 
  Line_End = Line_Start + Num_Lines - 1
  Num_Elem = Image%Number_Of_Elements      !Image%Number_Of_Elements is a global variable

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
     Nwp_Lon_Idx_x = I_Nwp_x(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx_x = J_Nwp_x(Elem_Idx,Line_Idx)
     Ivza = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !-- check if indices are valid
     if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0) cycle

     !-- check if cld-top pressure is valid
     if (ACHA%Pc(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- determine Level in profiles 
     !--- note, tpw and wind profiles are at the native nwp resolution
     call KNOWING_P_COMPUTE_T_Z_NWP(Nwp_Lon_Idx,Nwp_Lat_Idx, &
                                    ACHA%Pc(Elem_Idx,Line_Idx), &
                                    Tc_Temp,Zc_Temp,Level)

     !--- assign u and v winds
     U_Wnd_Temp = U_Wnd_Prof_Nwp(Level,Nwp_Lon_Idx,Nwp_Lat_Idx) 
     V_Wnd_Temp = V_Wnd_Prof_Nwp(Level,Nwp_Lon_Idx,Nwp_Lat_Idx) 

     !-- determine speed and direction
     Wnd_Spd_Cld_Top_Nwp_Pix(Elem_Idx,Line_Idx) = wind_speed(U_Wnd_Temp, V_Wnd_Temp)
     Wnd_Dir_Cld_Top_Nwp_Pix(Elem_Idx,Line_Idx) = wind_direction(U_Wnd_Temp, V_Wnd_Temp)

     !-- determine above cloud tpw (interpolate in x, y)
     Tpw_1 = INTERPOLATE_NWP(Tpw_Prof_Nwp(Level,Nwp_Lon_Idx,Nwp_Lat_Idx),     &
                             Tpw_Prof_Nwp(Level,Nwp_Lon_Idx_x,Nwp_Lat_Idx),   & 
                             Tpw_Prof_Nwp(Level,Nwp_Lon_Idx,Nwp_Lat_Idx_x),   & 
                             Tpw_Prof_Nwp(Level,Nwp_Lon_Idx_x,Nwp_Lat_Idx_x), & 
                             Lon_Nwp_Fac(Elem_Idx,Line_Idx),                  &
                             Lat_Nwp_Fac(Elem_Idx,Line_Idx))

     Tpw_2 = INTERPOLATE_NWP(Tpw_Prof_Nwp(Level+1,Nwp_Lon_Idx,Nwp_Lat_Idx),     &
                             Tpw_Prof_Nwp(Level+1,Nwp_Lon_Idx_x,Nwp_Lat_Idx),   & 
                             Tpw_Prof_Nwp(Level+1,Nwp_Lon_Idx,Nwp_Lat_Idx_x),   & 
                             Tpw_Prof_Nwp(Level+1,Nwp_Lon_Idx_x,Nwp_Lat_Idx_x), & 
                             Lon_Nwp_Fac(Elem_Idx,Line_Idx),                    &
                             Lat_Nwp_Fac(Elem_Idx,Line_Idx))


     !--- perform z interpolation
     dp = P_Std_Nwp(Level+1) - P_Std_Nwp(Level)

     dtpw = Tpw_2 - Tpw_1

     if (dp /= 0.0) then
       Tpw_Above_Cloud_Nwp_Pix(Elem_Idx,Line_Idx) = Tpw_1 + dtpw/dp * (ACHA%Pc(Elem_Idx,Line_Idx) - P_Std_Nwp(Level))
     else
       Tpw_Above_Cloud_Nwp_Pix(Elem_Idx,Line_Idx) = Tpw_1
     endif


  enddo Element_Loop_1
  enddo Line_Loop_1

end subroutine COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND_AND_TPW

!----------------------------------------------------------------------
! routine to interpolate pressure to flight level altitude.
!----------------------------------------------------------------------
subroutine COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min,Num_Lines,Pc_In,Alt_Out)

!--- Based on Sarah Monette calculations for HS3.  Her calculation is based on:
!--- http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf.

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  real, intent(inout), dimension(:,:) :: Alt_Out
  real, intent(in), dimension(:,:) :: Pc_In
  integer:: Num_Elem
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Line_Start
  integer:: Line_End
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  real:: Pc_Temp
  real:: Alt_Temp

  !--- Constants from Sarah Monette
  real, parameter :: PW1 = 227.9 ! hPa
  real, parameter :: PW2 = 56.89 ! hPa
  real, parameter :: PW3 = 11.01 ! hPa
  real, parameter :: P0 = 1013.25 ! hPa
  real, parameter :: LR_OVER_G = 0.190263
  real, parameter :: Z0 = 145422.16 ! feet
  real, parameter :: LN_1 = -20859.0
  real, parameter :: LN_2 = 149255.0
  real, parameter :: PN_4 = 0.000470034
  real, parameter :: PN_3 = -0.364267
  real, parameter :: PN_2 = 47.5627
  real, parameter :: PN_1 = -2647.45
  real, parameter :: PN_0 = 123842.0

  !--- save elements
  Line_Start = Line_Idx_Min
  Line_End = Line_Start + Num_Lines - 1
  Num_Elem = Image%Number_Of_Elements      !Image%Number_Of_Elements is a global variable

  !--- initialize
  Alt_Out = Missing_Value_Real4

  !----------------------------------------------------------
  ! loop through segment
  !----------------------------------------------------------
  Line_Loop_1: do Line_Idx = Line_Start, Line_End
  Element_Loop_1: do Elem_Idx = 1, Num_Elem

     !--- Initialize temporary value each time.
     Pc_Temp = Missing_Value_Real4
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
     if (Pc_In(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- Place valid pressure in temp variable for readability in the
     !--- calculations below.
     Pc_Temp = Pc_In(Elem_Idx,Line_Idx)

     !--- calculated altitude, in feet, from pressure.
     !--- 1st pivot point is directly from the pressure to
     !--- altitude from above reference.
     if (Pc_Temp > PW1) then
       Alt_Temp = (1.0 - (Pc_Temp/P0)**LR_OVER_G) * Z0
     endif

     !--- 2nd pivot point was modeled best with a natural log
     !--- fit.  From Sarah Monette.
     if (Pc_Temp <= PW1 .AND. Pc_Temp >= PW2) then
       Alt_Temp = LN_1 * LOG(Pc_Temp) + LN_2
     endif

     !--- 3rd pivot point. Modeled best with a polynomial
     !--- fit from Sarah Monette.
     if (Pc_Temp < PW2 .AND. Pc_Temp >= PW3) then
       Alt_Temp = (PN_4*Pc_Temp**4) + (PN_3*Pc_Temp**3) + (PN_2*Pc_Temp**2) + &
                  (PN_1*Pc_Temp) + PN_0
     endif

     if (Pc_Temp < PW3) then
       Alt_Temp = Missing_Value_Real4
     endif

     !--- Assign final altitude, in feet, to the level2 array.
     Alt_Out(Elem_Idx,Line_Idx) = Alt_Temp

  enddo Element_Loop_1
  enddo Line_Loop_1

end subroutine COMPUTE_ALTITUDE_FROM_PRESSURE

subroutine OPAQUE_CLOUD_HEIGHT() 

   integer:: Line_Idx
   integer:: Elem_Idx
   integer:: Num_Elements
   integer:: Num_Lines
   integer:: Nwp_Lon_Idx
   integer:: Nwp_Lat_Idx
   integer:: Vza_Rtm_Idx
   integer:: Tropo_Level_Idx
   integer:: Sfc_Level_Idx
   integer:: Level_Idx
   integer:: Level_Idx_Start
   integer:: Level_Idx_End

   Num_Elements = Image%Number_Of_Elements
   Num_Lines = Image%Number_Of_Lines_Per_Segment

   if (Sensor%Chan_On_Flag_Default(31)==sym%NO) return

   do Elem_Idx = 1, Num_Elements
     do Line_Idx = 1, Num_Lines

      !--- skip bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

      !--- indice aliases
      Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
      Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
      Vza_Rtm_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

      !-- check if indices are valid
      if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0 .or. Vza_Rtm_Idx < 0) cycle

      Tropo_Level_Idx = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
      Sfc_Level_Idx =   rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level
      Vza_Rtm_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

      !--- restrict levels to consider
      Level_Idx_Start = Tropo_Level_Idx
      Level_Idx_End = Sfc_Level_Idx

      !--- loop through levels
      Pc_Opaque_Cloud(Elem_Idx,Line_Idx) = P_Std_Rtm(Level_Idx_End)
      Zc_Opaque_Cloud(Elem_Idx,Line_Idx) = Rtm(Nwp_Lon_Idx, Nwp_Lat_Idx)%Z_Prof(Level_Idx_End)
      Tc_Opaque_Cloud(Elem_Idx,Line_Idx) = Rtm(Nwp_Lon_Idx, Nwp_Lat_Idx)%T_Prof(Level_Idx_End)

      level_loop: do Level_Idx = Level_Idx_Start, Level_Idx_End
        if (rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(31)%Rad_BB_Cloud_Profile(Level_Idx) >  &
            ch(31)%Rad_Toa(Elem_Idx,Line_Idx)) then
            Pc_Opaque_Cloud(Elem_Idx,Line_Idx) = P_Std_Rtm(Level_Idx - 1)
            Zc_Opaque_Cloud(Elem_Idx,Line_Idx) = Rtm(Nwp_Lon_Idx, Nwp_Lat_Idx)%Z_Prof(Level_Idx - 1)
            Tc_Opaque_Cloud(Elem_Idx,Line_Idx) = Rtm(Nwp_Lon_Idx, Nwp_Lat_Idx)%T_Prof(Level_Idx - 1)
           exit
       endif
      enddo Level_Loop

      enddo
   enddo

end subroutine OPAQUE_CLOUD_HEIGHT
!----------------------------------------------------------------------
! water vapor intercept
!----------------------------------------------------------------------
subroutine H2O_CLOUD_HEIGHT()

   integer:: Line_Idx
   integer:: Elem_Idx
   integer:: Num_Elements
   integer:: Num_Lines
   integer:: Nwp_Lon_Idx
   integer:: Nwp_Lat_Idx
   integer:: Vza_Rtm_Idx
   integer:: Tropo_Level_Idx
   integer:: Sfc_Level_Idx
   integer:: Level_Idx
   integer:: Level_Idx_Temp

   Num_Elements = Image%Number_Of_Elements
   Num_Lines = Image%Number_Of_Lines_Per_Segment

   if (Sensor%Chan_On_Flag_Default(31)==sym%NO) return
   if (Sensor%Chan_On_Flag_Default(27)==sym%NO) return

   do Elem_Idx = 1, Num_Elements
     do Line_Idx = 1, Num_Lines

      !--- skip bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

!     if (Cld_Type(Elem_Idx,Line_Idx) /= sym%CIRRUS_TYPE .and. & 
!         Cld_Type(Elem_Idx,Line_Idx) /= sym%OVERLAP_TYPE .and. & 
!         Cld_Type(Elem_Idx,Line_Idx) /= sym%OPAQUE_ICE_TYPE) then
!         cycle
!     endif


      !--- indice aliases
      Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
      Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
      Vza_Rtm_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

      !-- check if indices are valid
      if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0 .or. Vza_Rtm_Idx < 0) cycle

      Tropo_Level_Idx = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
      Sfc_Level_Idx =   rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level


      !--- call routine
      call INTERCEPT_CLOUD_PRESSURE(ch(31)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(31)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(31)%Rad_BB_Cloud_Profile, &
                               ch(27)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(27)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(27)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               P_Std_Rtm, &
                               Pc_H2O(Elem_Idx,Line_Idx), &
                               Level_Idx)      

      !--- compute temperature and height
      if (Pc_H2O(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
      if (Level_Idx == Missing_Value_Int4) cycle

      call KNOWING_P_COMPUTE_T_Z_NWP(Nwp_Lon_Idx,Nwp_Lat_Idx, &
                                     Pc_H2O(Elem_Idx,Line_Idx), &
                                     Tc_H2O(Elem_Idx,Line_Idx), &
                                     Zc_H2O(Elem_Idx,Line_Idx), &
                                     Level_Idx_Temp)
     

      enddo
   enddo

end subroutine H2O_CLOUD_HEIGHT
!----------------------------------------------------------------------
! water vapor intercept
!----------------------------------------------------------------------
subroutine CO2IRW_CLOUD_HEIGHT()

   integer:: Line_Idx
   integer:: Elem_Idx
   integer:: Num_Elements
   integer:: Num_Lines
   integer:: Nwp_Lon_Idx
   integer:: Nwp_Lat_Idx
   integer:: Vza_Rtm_Idx
   integer:: Tropo_Level_Idx
   integer:: Sfc_Level_Idx
   integer:: Level_Idx
   integer:: Level_Idx_Temp
   real, parameter:: Beta_Target = 1.1

   Num_Elements = Image%Number_Of_Elements
   Num_Lines = Image%Number_Of_Lines_Per_Segment

   if (Sensor%Chan_On_Flag_Default(27)==sym%NO) return
   if (Sensor%Chan_On_Flag_Default(31)==sym%NO) return
   if (Sensor%Chan_On_Flag_Default(33)==sym%NO) return

   do Elem_Idx = 1, Num_Elements
     do Line_Idx = 1, Num_Lines

      !--- skip bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
!     if (Cld_Type(Elem_Idx,Line_Idx) /= sym%CIRRUS_TYPE .and. & 
!         Cld_Type(Elem_Idx,Line_Idx) /= sym%OVERLAP_TYPE .and. & 
!         Cld_Type(Elem_Idx,Line_Idx) /= sym%OPAQUE_ICE_TYPE) then
!         cycle
!     endif


      !--- indice aliases
      Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
      Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
      Vza_Rtm_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

      !-- check if indices are valid
      if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0 .or. Vza_Rtm_Idx < 0) cycle

      Tropo_Level_Idx = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
      Sfc_Level_Idx =   rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level


      !--- call routine
      call SLICING_CLOUD_PRESSURE(ch(31)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(31)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(31)%Rad_BB_Cloud_Profile, &
                               ch(33)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(33)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(33)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               P_Std_Rtm, &
                               Beta_Target, &
                               Pc_CO2IRW(Elem_Idx,Line_Idx), &
                               Level_Idx)      
     
      !--- compute temperature and height
      if (Pc_Co2IRW(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
      if (Level_Idx == Missing_Value_Int4) cycle

      call KNOWING_P_COMPUTE_T_Z_NWP(Nwp_Lon_Idx,Nwp_Lat_Idx, &
                                     Pc_Co2IRW(Elem_Idx,Line_Idx), &
                                     Tc_Co2IRW(Elem_Idx,Line_Idx), &
                                     Zc_Co2IRW(Elem_Idx,Line_Idx), &
                                     Level_Idx_Temp)

      enddo
   enddo

end subroutine CO2IRW_CLOUD_HEIGHT

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine CTP_MULTILAYER()

   integer:: Line_Idx
   integer:: Elem_Idx
   integer:: Num_Elements
   integer:: Num_Lines

   Num_Elements = Image%Number_Of_Elements
   Num_Lines = Image%Number_Of_Lines_Per_Segment

   Ctp_Multilayer_Flag = Missing_Value_Int1

   !--- check if channels are on
   if (Sensor%Chan_On_Flag_Default(31)==sym%NO) return
   if (Sensor%Chan_On_Flag_Default(27)==sym%NO) return
   if (Sensor%Chan_On_Flag_Default(33)==sym%NO) return

 
   !--- loop through pixels in the segment
   do Elem_Idx = 1, Num_Elements
     do Line_Idx = 1, Num_Lines

      !--- skip bad pixels
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

      !---- filter on cloud type
      if (Cld_Type(Elem_Idx,Line_Idx) /= sym%CIRRUS_TYPE .and. & 
          Cld_Type(Elem_Idx,Line_Idx) /= sym%OVERLAP_TYPE .and. & 
          Cld_Type(Elem_Idx,Line_Idx) /= sym%OPAQUE_ICE_TYPE) then
          cycle
      endif


      !--- skip pixels without valid input
      if (Pc_H2O(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
      if (Pc_CO2IRW(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
      if (Ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle


      Ctp_Multilayer_Flag(Elem_Idx,Line_Idx) = sym%No

      !--- IRWCO2
      if (Ch(31)%Emiss_Tropo(Elem_Idx,Line_Idx) < 0.80) then
         if (abs(Pc_H2O(Elem_Idx,Line_Idx) - Pc_CO2IRW(Elem_Idx,Line_Idx)) > 50.0) then
           Ctp_Multilayer_Flag(Elem_Idx,Line_Idx) = sym%YES
         endif
      endif

      enddo
   enddo


end subroutine CTP_MULTILAYER

!----------------------------------------------------------------------
! Compute CO2 Slicing using Ch33, 34, 35 and ch 36
!----------------------------------------------------------------------
subroutine CO2_SLICING_CLOUD_HEIGHT(Num_Elem,Line_Idx_min,Num_Lines, &
                                    Pressure_Profile, &
                                    Pc_Co2,Tc_Co2,Zc_Co2)
  integer, intent(in):: Num_Elem
  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  real, intent(in), dimension(:):: Pressure_Profile
  real, intent(out), dimension(:,:):: Pc_Co2
  real, intent(out), dimension(:,:):: Tc_Co2
  real, intent(out), dimension(:,:):: Zc_Co2
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Line_Start
  integer:: Line_End
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  integer:: Vza_Rtm_Idx
  real:: Pc_33_34, Pc_34_35, Pc_35_36
  integer:: Tropo_Level_Idx
  integer:: Sfc_Level_Idx
  real:: Beta_Target
  integer (kind=int4), parameter:: COUNT_MIN_TEMPERATURE_CIRRUS = 2
  integer (kind=int4), parameter:: BOX_WIDTH_KM = 300
  real (kind=real4), parameter:: SOUNDER_RESOLUTION_KM = 20.0
  real (kind=real4), parameter:: PC_CIRRUS_MAX_THRESH = 440.0
  real (kind=real4), parameter:: EC_CIRRUS_MIN_THRESH = 0.2
  integer:: Lev_Idx_Temp
  integer,allocatable,dimension(:,:):: Pc_Lev_Idx

  real, dimension(3):: Pc_Temp

  allocate(Pc_Lev_Idx(Num_Elem,Num_Lines))

  Line_Start = Line_Idx_Min
  Line_End = Line_Start + Num_Lines - 1

  !--- intialize output
  Pc_Co2 = Missing_Value_Real4
  Tc_Co2 = Missing_Value_Real4
  Zc_Co2 = Missing_Value_Real4

  !---- check that all co2 channels are available
  if (Sensor%Chan_On_Flag_Default(33) == sym%NO .or. &
      Sensor%Chan_On_Flag_Default(34) == sym%NO .or. &
      Sensor%Chan_On_Flag_Default(35) == sym%NO .or. &
      Sensor%Chan_On_Flag_Default(36) == sym%NO) then
     return
  endif

  Line_Loop: do Line_Idx = Line_Start, Line_End
  Element_Loop: do Elem_Idx = 1, Num_Elem

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !--- skip data without sounder data
     if (ch(33)%Rad_Toa(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- indice aliases
     Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
     Vza_Rtm_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

     !-- check if indices are valid
     if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0 .or. Vza_Rtm_Idx < 0) cycle

     Tropo_Level_Idx = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
     Sfc_Level_Idx =   rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level

     !--- compute cloud top pressure using each channel pair
     Beta_Target = 1.0
     call SLICING_CLOUD_PRESSURE(ch(35)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(35)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(35)%Rad_BB_Cloud_Profile, &
                               ch(36)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(36)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(36)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               Pressure_Profile, &
                               Beta_Target, &
                               Pc_35_36, &
                               Pc_Lev_Idx(Elem_Idx,Line_Idx))
     Pc_Temp(1) = Pc_35_36
     if (Pc_35_36 /= Missing_Value_Real4) then
         Pc_Co2(Elem_Idx,Line_Idx) = Pc_35_36
         cycle
     endif

     call SLICING_CLOUD_PRESSURE(ch(34)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(34)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(34)%Rad_BB_Cloud_Profile, &
                               ch(35)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(35)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(35)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               Pressure_Profile, &
                               Beta_Target, &
                               Pc_34_35, &
                               Pc_Lev_Idx(Elem_Idx,Line_Idx))

     Pc_Temp(2) = Pc_34_35
     if (Pc_34_35 /= Missing_Value_Real4) then
         Pc_Co2(Elem_Idx,Line_Idx) = Pc_34_35
         cycle
     endif

     call SLICING_CLOUD_PRESSURE(ch(33)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(33)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(33)%Rad_BB_Cloud_Profile, &
                               ch(34)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(34)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(34)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               Pressure_Profile, &
                               Beta_Target, &
                               Pc_33_34, &
                               Pc_Lev_Idx(Elem_Idx,Line_Idx))

     Pc_Temp(3) = Pc_33_34
     if (Pc_33_34 /= Missing_Value_Real4) then
         Pc_Co2(Elem_Idx,Line_Idx) = Pc_33_34
         cycle
     endif

!--- alternative - take mean
!    Count_Valid = count(Pc_Temp /= Missing_Value_Real4)
!    if (Count_Valid > 0) then
!       Pc_Co2(Elem_Idx,Line_Idx) = sum(Pc_Temp, mask = Pc_Temp /= Missing_Value_Real4) / Count_Valid
!    endif
    
  enddo Element_Loop
  enddo Line_Loop

  !-------------------------------------------------------------------------
  ! determine temperature
  !-------------------------------------------------------------------------
  Line_Loop_2: do Line_Idx = Line_Start, Line_End
  Element_Loop_2: do Elem_Idx = 1, Num_Elem

     if (Pc_Co2(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
     if (Pc_Lev_Idx(Elem_Idx,Line_Idx) == Missing_Value_Int4) cycle

     !--- compute temperature
     call KNOWING_P_COMPUTE_T_Z_NWP(Nwp_Lon_Idx,Nwp_Lat_Idx, &
                                    Pc_Co2(Elem_Idx,Line_Idx), &
                                    Tc_Co2(Elem_Idx,Line_Idx), &
                                    Zc_Co2(Elem_Idx,Line_Idx), &
                                    Lev_Idx_Temp)

     !-- compute emissivity
     Ec_Co2(Elem_Idx,Line_Idx) = EMISSIVITY(ch(33)%Rad_Toa(Elem_Idx,Line_Idx),  &
                                 ch(33)%Rad_Toa_Clear(Elem_Idx,Line_Idx),  &
                                 rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(33)%Rad_BB_Cloud_Profile(Pc_Lev_Idx(Elem_Idx,Line_Idx)))

  enddo Element_Loop_2
  enddo Line_Loop_2

  deallocate(Pc_Lev_Idx)

end subroutine  CO2_SLICING_CLOUD_HEIGHT

!----------------------------------------------------------------------
! Compute CO2 Slicing using Ch33, 34, 35 and ch 36
!----------------------------------------------------------------------
subroutine CO2_SLICING_CLOUD_HEIGHT_NEW(Num_Elem,Line_Idx_min,Num_Lines, &
                                    Pressure_Profile,Cloud_Mask, &
                                    Pc_Co2,Tc_Co2,Zc_Co2)
  integer, intent(in):: Num_Elem
  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  integer(kind=int1), intent(in), dimension(:,:):: Cloud_Mask
  real, intent(in), dimension(:):: Pressure_Profile
  real, intent(out), dimension(:,:):: Pc_Co2
  real, intent(out), dimension(:,:):: Tc_Co2
  real, intent(out), dimension(:,:):: Zc_Co2
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Line_Start
  integer:: Line_End
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  integer:: Vza_Rtm_Idx
  real:: Pc_33_34, Pc_34_35, Pc_35_36
  integer:: Tropo_Level_Idx
  integer:: Sfc_Level_Idx
  real:: Beta_Target
  integer (kind=int4), parameter:: COUNT_MIN_TEMPERATURE_CIRRUS = 2
  integer (kind=int4), parameter:: BOX_WIDTH_KM = 300
  real (kind=real4), parameter:: SOUNDER_RESOLUTION_KM = 20.0
  real (kind=real4), parameter:: PC_CIRRUS_MAX_THRESH = 440.0
  real (kind=real4), parameter:: EC_CIRRUS_MIN_THRESH = 0.2
  integer:: Lev_Idx_Temp
  integer,allocatable,dimension(:,:):: Pc_Lev_Idx

  real, dimension(3):: Pc_Temp
  integer:: Count_Valid
  logical, allocatable, dimension(:):: Sounder_Fov_Mask
  integer:: Number_Sounder_Fov
  allocate(Pc_Lev_Idx(Num_Elem,Num_Lines))

  Line_Start = Line_Idx_Min
  Line_End = Line_Start + Num_Lines - 1

  !--- intialize output
  Pc_Co2 = Missing_Value_Real4
  Tc_Co2 = Missing_Value_Real4
  Zc_Co2 = Missing_Value_Real4

  !---- check that all co2 channels are available
  if (Sensor%Chan_On_Flag_Default(33) == sym%NO .or. &
      Sensor%Chan_On_Flag_Default(34) == sym%NO .or. &
      Sensor%Chan_On_Flag_Default(35) == sym%NO .or. &
      Sensor%Chan_On_Flag_Default(36) == sym%NO) then
     return
  endif

  !------- make a mask
  Number_Sounder_Fov = maxval(Nav%Sounder_Fov_Segment_Idx)
  allocate(Sounder_Fov_Mask(Number_Sounder_Fov))
  Sounder_Fov_Mask = .false.

  Line_Loop: do Line_Idx = Line_Start, Line_End
  Element_Loop: do Elem_Idx = 1, Num_Elem

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !--- skip data without sounder data
     if (ch(33)%Rad_Toa(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- indice aliases
     Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
     Vza_Rtm_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

     !-- check if indices are valid
     if (Nwp_Lon_Idx < 0 .or. Nwp_Lat_Idx < 0 .or. Vza_Rtm_Idx < 0) cycle

     Tropo_Level_Idx = rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Tropo_Level
     Sfc_Level_Idx =   rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%Sfc_Level

     !--- only do this for appropriate cloud types
     if (Cloud_Mask(Elem_Idx,Line_Idx) == sym%CLEAR_TYPE) cycle
     if (Cloud_Mask(Elem_Idx,Line_Idx) == sym%PROB_CLEAR_TYPE) cycle
     if (Nav%Sounder_Fov_Segment_Idx(Elem_Idx,Line_Idx) <= 0) cycle

     !--- if already done, move on
     if (Sounder_Fov_Mask(Nav%Sounder_Fov_Segment_Idx(Elem_Idx,Line_Idx)) .eqv. .true.) cycle

     !--- compute cloud top pressure using each channel pair
     Beta_Target = 1.0

     call SLICING_CLOUD_PRESSURE(ch(35)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(35)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(35)%Rad_BB_Cloud_Profile, &
                               ch(36)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(36)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(36)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               Pressure_Profile, &
                               Beta_Target, &
                               Pc_35_36, &
                               Pc_Lev_Idx(Elem_Idx,Line_Idx))
     Pc_Temp(1) = Pc_35_36
     if (Pc_35_36 /= Missing_Value_Real4) then
         Pc_Co2(Elem_Idx,Line_Idx) = Pc_35_36
         cycle
     endif

     call SLICING_CLOUD_PRESSURE(ch(34)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(34)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(34)%Rad_BB_Cloud_Profile, &
                               ch(35)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(35)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(35)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               Pressure_Profile, &
                               Beta_Target, &
                               Pc_34_35, &
                               Pc_Lev_Idx(Elem_Idx,Line_Idx))

     Pc_Temp(2) = Pc_34_35
     if (Pc_34_35 /= Missing_Value_Real4) then
         Pc_Co2(Elem_Idx,Line_Idx) = Pc_34_35
         cycle
     endif

     call SLICING_CLOUD_PRESSURE(ch(33)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(33)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(33)%Rad_BB_Cloud_Profile, &
                               ch(34)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(34)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(34)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               Pressure_Profile, &
                               Beta_Target, &
                               Pc_33_34, &
                               Pc_Lev_Idx(Elem_Idx,Line_Idx))

     Pc_Temp(3) = Pc_33_34
     if (Pc_33_34 /= Missing_Value_Real4) then
         Pc_Co2(Elem_Idx,Line_Idx) = Pc_33_34
         cycle
     endif

     Count_Valid = count(Pc_Temp /= Missing_Value_Real4)
     if (Count_Valid > 0) then
        Pc_Co2(Elem_Idx,Line_Idx) = sum(Pc_Temp, mask = Pc_Temp /= Missing_Value_Real4) / Count_Valid
     endif



     !--- set this fov as being processed
     Sounder_Fov_Mask(Nav%Sounder_Fov_Segment_Idx(Elem_Idx,Line_Idx)) = .true.

     !--- spread out information over all imager pixels within the sounder fov
     where(Nav%Sounder_Fov_Segment_Idx == Nav%Sounder_Fov_Segment_Idx(Elem_Idx,Line_Idx))
         Pc_Co2 = Pc_Co2(Elem_Idx,Line_Idx)
     end where
    
  enddo Element_Loop
  enddo Line_Loop


  deallocate(Sounder_Fov_Mask)

  !-------------------------------------------------------------------------
  ! determine temperature
  !-------------------------------------------------------------------------
  Line_Loop_2: do Line_Idx = Line_Start, Line_End
  Element_Loop_2: do Elem_Idx = 1, Num_Elem

     if (Pc_Co2(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- compute temperature
     call KNOWING_P_COMPUTE_T_Z_NWP(Nwp_Lon_Idx,Nwp_Lat_Idx, &
                                    Pc_Co2(Elem_Idx,Line_Idx), &
                                    Tc_Co2(Elem_Idx,Line_Idx), &
                                    Zc_Co2(Elem_Idx,Line_Idx), &
                                    Lev_Idx_Temp)

     !-- compute emissivity
     Pc_Lev_Idx(Elem_Idx,Line_Idx) = Lev_Idx_Temp
     if (Pc_Lev_Idx(Elem_Idx,Line_Idx) /= MISSING_VALUE_INT2) then
        Ec_Co2(Elem_Idx,Line_Idx) = EMISSIVITY(ch(33)%Rad_Toa(Elem_Idx,Line_Idx),  &
                                    ch(33)%Rad_Toa_Clear(Elem_Idx,Line_Idx),  &
                                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(33)%Rad_BB_Cloud_Profile(Pc_Lev_Idx(Elem_Idx,Line_Idx)))
     endif
  enddo Element_Loop_2
  enddo Line_Loop_2

  deallocate(Pc_Lev_Idx)

end subroutine CO2_SLICING_CLOUD_HEIGHT_NEW


!-------------------------------------------------------------------------
! Make Emissivity from Sounder Values
!-------------------------------------------------------------------------
subroutine SOUNDER_EMISSIVITY()
 
  real:: R4_Dummy
  integer:: Num_Elem
  integer:: Num_Lines
  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Lev_Idx
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  integer:: Vza_Rtm_Idx

  if (Sensor%Chan_On_Flag_Default(31) == sym%NO) return
  if (.not. allocated(Cld_Press_Sounder)) return
  if (.not. allocated(Cld_Emiss_Sounder)) return

  Num_Elem = Image%Number_Of_Elements  
  Num_Lines = Image%Number_Of_Lines_Per_Segment

  !-------------------------------------------------------------------------
  ! loop through and determine emissivity
  !-------------------------------------------------------------------------
  Line_Loop: do Line_Idx = 1, Num_Lines
  Element_Loop: do Elem_Idx = 1, Num_Elem

     !--- skip data without sounder data
     if (ch(31)%Rad_Toa(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
     if (Cld_Press_Sounder(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- indice aliases
     Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
     Vza_Rtm_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

     if (Cld_Press_Sounder(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
     if (Vza_Rtm_Idx < 1) cycle
     if (Nwp_Lat_Idx < 1) cycle
     if (Nwp_Lon_Idx < 1) cycle

     !--- compute temperature
     call KNOWING_P_COMPUTE_T_Z_NWP(Nwp_Lon_Idx,Nwp_Lat_Idx, &
                                    Cld_Press_Sounder(Elem_Idx,Line_Idx), &
                                    R4_Dummy, &
                                    R4_Dummy, &
                                    Lev_Idx)

     !-- compute emissivity
     if (Lev_Idx /= MISSING_VALUE_INT2) then
        Cld_Emiss_Sounder(Elem_Idx,Line_Idx) = EMISSIVITY(ch(31)%Rad_Toa(Elem_Idx,Line_Idx),  &
                                    ch(31)%Rad_Toa_Clear(Elem_Idx,Line_Idx),  &
                                    rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(31)%Rad_BB_Cloud_Profile(Lev_Idx))
     endif
  enddo Element_Loop
  enddo Line_Loop

end subroutine SOUNDER_EMISSIVITY


!-------------------------------------------------------------------------
! spatially interpolate
!-------------------------------------------------------------------------
subroutine MAKE_CIRRUS_PRIOR_TEMPERATURE(Tc_Co2, Pc_Co2, Ec_Co2,  &
                                         Tc_Cirrus_Background, Zc_Cirrus_Background)

  real, intent(in), dimension(:,:):: Tc_Co2
  real, intent(in), dimension(:,:):: Pc_Co2
  real, intent(in), dimension(:,:):: Ec_Co2
  real, intent(out), dimension(:,:):: Tc_Cirrus_Background
  real, intent(out), dimension(:,:):: Zc_Cirrus_Background

  integer:: Num_Elem
  integer:: Num_Lines
  integer:: Elem_Idx
  integer:: Line_Idx

  integer (kind=int4), parameter:: COUNT_MIN_TEMPERATURE_CIRRUS = 2
  integer (kind=int4), parameter:: BOX_WIDTH_KM = 300
  real (kind=real4), parameter:: SOUNDER_RESOLUTION_KM = 20.0
  real (kind=real4), parameter:: PC_CIRRUS_MAX_THRESH = 440.0
  real (kind=real4), parameter:: EC_CIRRUS_MIN_THRESH = 0.2
  integer (kind=int4):: Box_Width
  real:: Count_Temporary, Sum_Temporary, Temperature_Temporary
  real:: P_Temporary, Z_Temporary
  integer:: i1, i2, j1, j2, Lev_Idx, Ierr
  integer:: Nwp_Lon_Idx, Nwp_Lat_Idx
  real:: Z_Interp_Weight

  logical, dimension(:,:), allocatable:: Mask
  real, dimension(:,:), allocatable:: Tc_Cirrus_Background_Temp
  real, dimension(:,:), allocatable:: Zc_Cirrus_Background_Temp

  Num_Elem = Image%Number_Of_Elements
  Num_Lines = Image%Number_Of_Lines_Read_This_Segment

  !--- compute size of averaging window
  call COMPUTE_BOX_WIDTH(SOUNDER_RESOLUTION_KM,Box_Width_KM,Box_Width)

  !--- allocate temporary arrays needed here
  allocate(Mask(Num_Elem,Num_Lines))
  allocate(Tc_Cirrus_Background_Temp(Num_Elem,Num_Lines))
  allocate(Zc_Cirrus_Background_Temp(Num_Elem,Num_Lines))

  Tc_Cirrus_Background_Temp = Missing_Value_Real4
  Zc_Cirrus_Background_Temp = Missing_Value_Real4

  Mask = .false.

  !--- create mask which identifies pixels to be included in analysis
  where(Pc_Co2 /= Missing_Value_Real4 .and. &
        Pc_Co2 < PC_CIRRUS_MAX_THRESH .and. &
        Ec_Co2 > EC_CIRRUS_MIN_THRESH)
      Mask = .true.
  end where

  Element_Loop: do Elem_Idx = 1, Num_Elem, Box_Width/2

      i1 = min(Num_Elem,max(1,Elem_Idx - Box_Width))
      i2 = min(Num_Elem,max(1,Elem_Idx + Box_Width))

      Line_Loop: do Line_Idx = 1, Num_Lines, Box_Width/2
          j1 = min(Num_Lines,max(1,Line_Idx - Box_Width))
          j2 = min(Num_Lines,max(1,Line_Idx + Box_Width))

          Count_Temporary = count(Mask(i1:i2,j1:j2))
          Sum_Temporary = sum(Tc_Co2(i1:i2,j1:j2),Mask(i1:i2,j1:j2))
          if (Count_Temporary > COUNT_MIN_TEMPERATURE_CIRRUS) then
              Temperature_Temporary = Sum_Temporary / Count_Temporary
          else
              Temperature_Temporary = Missing_Value_Real4
          endif

          !---- interpolate for Zc
          Nwp_Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
          if (Nwp_Lon_Idx < 0) cycle
          Nwp_Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
          if (Nwp_Lat_Idx < 0) cycle
          if (Temperature_Temporary /= Missing_Value_Real4) then
               call KNOWING_T_COMPUTE_P_Z_NWP(Nwp_Lon_Idx,Nwp_Lat_Idx,P_Temporary, &
                                      Temperature_Temporary,Z_Temporary, &
                                      Lev_Idx,Elem_Idx,Line_Idx,Ierr,Z_Interp_Weight)
          else
              Z_Temporary = Missing_Value_Real4
          endif

          if (Temperature_Temporary /= Missing_Value_Real4) then
             Tc_Cirrus_Background_Temp(i1:i2,j1:j2) = Temperature_Temporary
          endif

          if (Z_Temporary /= Missing_Value_Real4) then
             Zc_Cirrus_Background_Temp(i1:i2,j1:j2) = Z_Temporary
          endif

      enddo Line_Loop
   enddo Element_Loop

   Tc_Cirrus_Background = Tc_Cirrus_Background_Temp
   Zc_Cirrus_Background = Zc_Cirrus_Background_Temp

  if (allocated(Mask)) deallocate(Mask)
  if (allocated(Tc_Cirrus_Background_Temp)) deallocate(Tc_Cirrus_Background_Temp)
  if (allocated(Zc_Cirrus_Background_Temp)) deallocate(Zc_Cirrus_Background_Temp)

end subroutine MAKE_CIRRUS_PRIOR_TEMPERATURE

   !====================================================================
   ! FUNCTION Name: BETA_RATIO
   !
   ! Function:
   !  Computes the beta ratio for two Emissivities. 
   !
   ! Input:  Emiss_top - emissivity in the numerator
   !         Emiss_bot - emissivity in the denominator
   !
   ! Output: Beta - the beta value from the two emissivities
   !
   !  Traditional Beta 11,12 um = BETA_RATIO(11,12)
   !
   !====================================================================
   function BETA_RATIO(Emiss_Bot, Emiss_Top) result(Beta)
      real(kind=real4), intent(in) :: Emiss_Top
      real(kind=real4), intent(in) :: Emiss_Bot
      real(kind=real4) :: Beta

      Beta = Missing_Value_Real4

      if (Emiss_Top > 0.0 .and. Emiss_Top < 1.0 .and. &
          Emiss_Bot > 0.0 .and. Emiss_Bot < 1.0) then

         Beta = alog(1.0 - Emiss_Top)/alog(1.0 - Emiss_Bot)

      end if

      return

   end function BETA_RATIO

   !====================================================================
   ! Function Name: EMISSIVITY
   !
   ! Function:
   !  Computes the  effective emissivity
   !
   ! Input:  Rad_Toa - channel radiance at top of atmosphere(toa)
   !         Rad_Clear_Tau - channel radiance at toa for clear skies
   !         Rad_Cloud_BB_Toa - channel radiance at TOA if cloud were a
   !         Black-Body
   !
   ! Output: Emiss - the effective cloud emissivity
   !
   !====================================================================
   function EMISSIVITY(Radiance_Toa, Radiance_Clear_Toa, Radiance_Cloud_BB_Toa) result(Emiss)
      real(kind=real4), intent(in) :: Radiance_Toa
      real(kind=real4), intent(in) :: Radiance_Clear_Toa
      real(kind=real4), intent(in) :: Radiance_Cloud_BB_Toa
      real(kind=real4) :: Emiss

      Emiss = Missing_Value_Real4

      if (Radiance_Cloud_BB_Toa /= Radiance_Clear_Toa) then
          Emiss = (Radiance_Toa - Radiance_Clear_Toa) / &
            (Radiance_Cloud_BB_Toa - Radiance_Clear_Toa)
       end if

      return

   end function EMISSIVITY
!----------------------------------------------------------------------
!--- determine cirrus box width
!---
!--- Sensor_Resolution_KM = the nominal resolution in kilometers
!--- Box_Width_KM = the width of the desired box in kilometers
!--- Box_Half_Width = the half width of the box in pixel-space
!----------------------------------------------------------------------
subroutine COMPUTE_BOX_WIDTH(Sensor_Resolution_KM,Box_Width_KM, &
                             Box_Half_Width)

   real, intent(in):: Sensor_Resolution_KM
   integer, intent(in):: Box_Width_KM
   integer, intent(out):: Box_Half_Width

   if (Sensor_Resolution_KM <= 0.0) then
       Box_Half_Width = 20
   else
       Box_Half_Width = int((Box_Width_KM / Sensor_Resolution_KM) / 2)
   endif

end subroutine COMPUTE_BOX_WIDTH

!----------------------------------------------------------------------
! Compute height where level transmission to space is opaque
!----------------------------------------------------------------------
subroutine OPAQUE_TRANSMISSION_HEIGHT()
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Lon_Idx
  integer:: Lat_Idx
  integer:: Zen_Idx
  integer:: Chan_Idx
  integer:: Lev_Idx
  real, dimension(:), pointer:: Trans_Prof
  real, dimension(:), pointer:: Z_Prof
  real, parameter:: Trans_Limit = 0.01

  do Elem_Idx = 1, Image%Number_Of_Elements
     do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment
      
     !--- indice aliases
     Lon_Idx = I_Nwp(Elem_Idx,Line_Idx)
     Lat_Idx = J_Nwp(Elem_Idx,Line_Idx)
     Zen_Idx = Zen_Idx_Rtm(Elem_Idx,Line_Idx)

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     !-- check if indices are valid
     if (Lon_Idx < 0 .or. Lat_Idx < 0) cycle

     Z_Prof => Rtm(Lon_Idx,Lat_Idx)%Z_Prof

     do Chan_Idx = 27,38

         if ( sensor % Chan_On_Flag_Default(Chan_Idx) == sym % no ) cycle

         Trans_Prof => Rtm(Lon_Idx,Lat_Idx)%d(Zen_Idx)%ch(Chan_Idx)%Trans_Atm_Profile 

         Lev_Idx = minloc(abs(Trans_Prof-Trans_Limit),1)   
         ch(Chan_Idx)%Opaque_Height(Elem_Idx,Line_Idx) = max(Z_Prof(Lev_Idx),sfc%Zsfc(Elem_Idx,Line_Idx))

         Trans_Prof => null()

      enddo

      Z_Prof => null()
        
     enddo
  enddo

end subroutine OPAQUE_TRANSMISSION_HEIGHT

!----------------------------------------------------------------------
! Compute CSBT Cloud Masks
! must be called at end of cloud processing chain (After ACHA)
!----------------------------------------------------------------------
subroutine COMPUTE_CSBT_CLOUD_MASKS()

  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Chan_Idx
  real, parameter:: Ch31_Mask_Cld_Prob_Max = 0.1
  real, parameter:: Covar_Ch27_Ch31_Max = 1.0
   
! if ( sensor % Chan_On_Flag_Default(27) == sym % no ) return
   
  do Elem_Idx = 1, Image%Number_Of_Elements
     do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     if (Posterior_Cld_Probability(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     do Chan_Idx = 27, 38

         if ( sensor % Chan_On_Flag_Default(Chan_Idx) == sym % no ) cycle
         if (ch(Chan_Idx)%Bt_Toa(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

           !--- initialize to cloudy
           ch(Chan_Idx)%CSBT_Mask(Elem_Idx,Line_Idx) = 3
           ch(Chan_Idx)%CSBT_Mask(Elem_Idx,Line_Idx) = 3

           !--- if full mask is clear, set channel masks to clear
           if (Posterior_Cld_Probability(Elem_Idx,Line_Idx) <= 0.10) then
              ch(Chan_Idx)%CSBT_Mask(Elem_Idx,Line_Idx) = 0
              cycle
           endif

           if (Posterior_Cld_Probability(Elem_Idx,Line_Idx) <= 0.25) then
              ch(Chan_Idx)%CSBT_Mask(Elem_Idx,Line_Idx) = 1
              cycle
           endif

           !if (Covar_Ch27_Ch31_5x5(Elem_idx,Line_Idx) <= Covar_Ch27_Ch31_Max) then
              if (Acha%Zc(Elem_Idx,Line_Idx) <= ch(Chan_Idx)%Opaque_Height(Elem_Idx,Line_Idx)) then
                ch(Chan_Idx)%CSBT_Mask(Elem_Idx,Line_Idx) = 2
              endif
           !endif

       enddo

     enddo
  enddo

end subroutine COMPUTE_CSBT_CLOUD_MASKS
!----------------------------------------------------------------------
! compute a mask that uniquely identifies all sounder fovs and make
! a mask that tells if that sounder fov has been processed
!
! input: 
! Nav%Sounder_X = x index of each sounder fov grouping
! Nav%Sounder_Y = y index of each sounder fov grouping (this is relative
!                 to the whole file not the segment
! Nav%Sounder_Fov = sounder index within each grouping
!
! output:
! Nav%Sounder_Fov_Segment_Idx which is a pixel-level
! index where each sounder fov in this segment is numbered from 1 to
! the total number of sounder fovs.
!----------------------------------------------------------------------
subroutine  COMPUTE_SOUNDER_MASK_SEGMENT()
  integer:: X_Idx, Y_Idx, Fov_Idx
  integer:: Fov_Segment_Idx
  integer:: Num_Line
  integer:: Num_Elem
  integer:: Sounder_X_Min
  integer:: Sounder_X_Max
  integer:: Sounder_Y_Min
  integer:: Sounder_Y_Max
  integer:: Sounder_Fov_Min
  integer:: Sounder_Fov_Max
  logical, allocatable, dimension(:,:):: Sounder_Valid_Mask

  !--- make a logical mask for identifying pixels with sounder data
  !--- needs to logical to be an argument for maxval and minval
  Num_Elem = Image%Number_Of_Elements  
  Num_Line = Image%Number_Of_Lines_Per_Segment
  allocate(Sounder_Valid_Mask(Num_Elem,Num_Line))
  Sounder_Valid_Mask =  .false.
  where(Nav%Sounder_X > 0 .and. Nav%Sounder_Y > 0)
     Sounder_Valid_Mask = .true.
  end where

  !--- determine ranges in sounder indices
  Sounder_X_Min = minval(Nav%Sounder_X,Sounder_Valid_Mask)
  Sounder_X_Max = maxval(Nav%Sounder_X,Sounder_Valid_Mask)
  Sounder_Y_Min = minval(Nav%Sounder_Y,Sounder_Valid_Mask)
  Sounder_Y_Max = maxval(Nav%Sounder_Y,Sounder_Valid_Mask)
  Sounder_Fov_Min = minval(Nav%Sounder_Fov,Sounder_Valid_Mask)
  Sounder_Fov_Max = maxval(Nav%Sounder_Fov,Sounder_Valid_Mask)

  !--- for HIRS there are no valid values for Sounder_Fov, should
  !--- treat as if all values are 1
  if (Sounder_Fov_Min == 0) Sounder_Fov_Min = 1
  if (Sounder_Fov_Max == 0) Sounder_Fov_Max = 1

  !--- loop over sounder indices and construct soundef fov idx for this segment
  Fov_Segment_Idx = 0
  Nav%Sounder_Fov_Segment_Idx = 0
  do X_Idx = Sounder_X_Min, Sounder_X_Max
    do Y_Idx = Sounder_Y_Min, Sounder_Y_Max
      do Fov_Idx = Sounder_Fov_Min, Sounder_Fov_Max

        where(Nav%Sounder_X == X_Idx .and. &
              Nav%Sounder_Y == Y_Idx .and. &
              Nav%Sounder_Fov == Fov_Idx)

          Nav%Sounder_Fov_Segment_Idx = Fov_Segment_Idx + 1

        end where

        Fov_Segment_Idx = maxval(Nav%Sounder_Fov_Segment_Idx)

      enddo
    enddo
  enddo

  !-- ensure missing value for pixels without sounder data
  where(Nav%Sounder_Fov_Segment_Idx == 0)
          Nav%Sounder_Fov_Segment_Idx = MISSING_VALUE_INT2
  end where

  !--- destroy logical mask created above
  deallocate(Sounder_Valid_Mask)

end subroutine  COMPUTE_SOUNDER_MASK_SEGMENT
!====================================================================
! Function Name: INTERCEPT_CLOUD_PRESSURE
!
! Function: estimate the cloud temperature/height/pressure
!
! Description: Use the 11um and an absorbing obs and the RTM cloud BB profiles
!              to perform intercept cloud pressure on a pixel level. Filters
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
subroutine  INTERCEPT_CLOUD_PRESSURE(Rad_Window, &
                                     Rad_Window_Clear, &
                                     Rad_Window_BB_Profile, &
                                     Rad_Abs, &
                                     Rad_Abs_Clear,  &
                                     Rad_Abs_BB_Profile,  &
                                     Tropo_Level, &
                                     Sfc_Level, &
                                     P_Prof, &
                                     Pc,  &
                                     Solution_Level)

  real, intent(in):: Rad_Window
  real, intent(in):: Rad_Abs
  real, intent(in), dimension(:):: Rad_Window_BB_Profile
  real, intent(in):: Rad_Window_Clear
  real, intent(in):: Rad_Abs_Clear
  real, intent(in), dimension(:):: Rad_Abs_BB_Profile
  integer, intent(in):: Sfc_Level
  integer, intent(in):: Tropo_Level
  real, intent(in), dimension(:):: P_Prof
  real, intent(out) :: Pc
  integer:: Solution_Level

  real:: Rad_Abs_BB_Prediction
  real:: Slope
  real:: Intercept
  real:: Denominator
  integer:: ilev

  real:: Rad_Window_Thresh
  real:: Rad_Abs_Thresh
  real, parameter:: Clear_Diff_Thresh = 0.02

  !--- initialize
  Pc = MISSING_VALUE_REAL4

  !--- set thresholds to detect useable signal
  Rad_Window_Thresh = Clear_Diff_Thresh*Rad_Window_Clear
  Rad_Abs_Thresh = Clear_Diff_Thresh*Rad_Abs_Clear

  !--- determine if a solution should be attempted
  if (Rad_Window_Clear - Rad_Window < Rad_Window_Thresh) return
  if (Rad_Window == MISSING_VALUE_REAL4) return
  if (Rad_Window_Clear == MISSING_VALUE_REAL4) return

  if (Rad_Abs_Clear - Rad_Abs < Rad_Abs_Thresh) return
  if (Rad_Abs == MISSING_VALUE_REAL4) return
  if (Rad_Abs_Clear == MISSING_VALUE_REAL4) return

 !--- attempt a solution than colder than tropo
 if (Rad_Window < Rad_Window_BB_Profile(Tropo_Level)) then

     Solution_Level = Tropo_Level

 else   !if not, attempt solution

     !--- determine linear regress of abs (y)  as a function of window (x)
      Denominator =  Rad_Window - Rad_Window_Clear

      if (Denominator < 0.0) then
             Slope = (Rad_Abs - Rad_Abs_Clear) / (Denominator)
             Intercept = Rad_Abs - Slope*Rad_Window
      else
            return
      endif

      !--- brute force solution
      Solution_Level = 0

      do ilev = Tropo_Level+1, Sfc_Level
          Rad_Abs_BB_Prediction = Slope*Rad_Window_BB_Profile(ilev) + Intercept

          if (Rad_Abs_BB_Prediction < 0) cycle

          if ((Rad_Abs_BB_Prediction > Rad_Abs_BB_Profile(ilev-1)) .and. &
               (Rad_Abs_BB_Prediction <= Rad_Abs_BB_Profile(ilev))) then
               Solution_Level = ilev
               exit
          endif

      enddo

 endif    !tropopause check

 !--- adjust back to full Rtm profile indices
 if (Solution_Level > 0) then
       Pc = P_Prof(Solution_Level)
 endif

end subroutine INTERCEPT_CLOUD_PRESSURE
!==================================================================================================
! Slicing Cloud Pressure - Requires two absorping channels 
!==================================================================================================
subroutine SLICING_CLOUD_PRESSURE(Ch_X_Rad_Toa, &
                                Ch_X_Rad_Toa_Clear, &
                                Ch_X_Rad_BB_Cloud_Profile, &
                                Ch_Y_Rad_Toa, &
                                Ch_Y_Rad_Toa_Clear, &
                                Ch_Y_Rad_BB_Cloud_Profile, &
                                Tropo_Level_Idx, &
                                Sfc_Level_Idx, &
                                P_Prof, &
                                Beta_Target, &
                                Pc,  &
                                Solution_Level)

   real, intent(in):: Ch_X_Rad_Toa
   real, intent(in):: Ch_X_Rad_Toa_Clear
   real, intent(in), dimension(:):: Ch_X_Rad_BB_Cloud_Profile
   real, intent(in):: Ch_Y_Rad_Toa
   real, intent(in):: Ch_Y_Rad_Toa_Clear
   real, intent(in), dimension(:):: Ch_Y_Rad_BB_Cloud_Profile
   integer, intent(in):: Tropo_Level_Idx
   integer, intent(in):: Sfc_Level_Idx
   real, intent(in), dimension(:):: P_Prof
   real, intent(in):: Beta_Target
   real, intent(out):: Pc
   integer, intent(out):: Solution_Level

   integer:: Lev_Idx
   integer:: Lev_Idx_Start
   integer:: Lev_Idx_End
   integer:: Num_Levels
   real:: Beta_X_Y
   real:: Beta_X_Y_Prev
   real:: Ch_X_Emissivity
   real:: Ch_Y_Emissivity

   !--- set levels for solution
   Num_Levels = size(P_Prof)
   Lev_Idx_Start = Tropo_Level_Idx
   Lev_Idx_End = Sfc_Level_Idx

   !--- initialize
   Pc = Missing_Value_Real4
   Beta_X_Y_Prev = Missing_Value_Real4
   Solution_Level = Missing_Value_Int4

   !--- loop through levels
   do Lev_Idx = Lev_Idx_Start, Lev_Idx_End  

      !---  make emissivities
      Ch_X_Emissivity = EMISSIVITY(Ch_X_Rad_Toa, Ch_X_Rad_Toa_Clear, Ch_X_Rad_BB_Cloud_Profile(Lev_Idx))
      if (Ch_X_Emissivity <= 0.0 .or. Ch_X_Emissivity >= 1.0) cycle

      Ch_Y_Emissivity = EMISSIVITY(Ch_Y_Rad_Toa, Ch_Y_Rad_Toa_Clear, Ch_Y_Rad_BB_Cloud_Profile(Lev_Idx))
      if (Ch_Y_Emissivity <= 0.0 .or. Ch_Y_Emissivity >= 1.0) cycle

      !--- make beta ratio
      Beta_X_Y = BETA_RATIO(Ch_X_Emissivity, Ch_Y_Emissivity)

      if (Beta_X_Y == Missing_Value_Real4) cycle

      !--- check is target beta is achieved
      if (Beta_X_Y <= Beta_Target .and. Beta_X_Y_Prev > Beta_Target .and. Beta_X_Y_Prev /= Missing_Value_Real4) then
         Solution_Level = Lev_Idx
         exit
      endif

      if (Beta_X_Y >= Beta_Target .and. Beta_X_Y_Prev < Beta_Target .and. Beta_X_Y_Prev /= Missing_Value_Real4) then
         Solution_Level = Lev_Idx
         exit
      endif

      Beta_X_Y_Prev = Beta_X_Y
      
   enddo

   !--- set output pressure if solution found
   if (Solution_Level /= Missing_Value_Int4) then
      Pc = P_Prof(Solution_Level)
   endif

end subroutine SLICING_CLOUD_PRESSURE
!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module CLOUD_HEIGHT_ROUTINES
