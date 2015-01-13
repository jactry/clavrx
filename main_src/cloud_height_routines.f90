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

  public::  COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND
  public::  COMPUTE_ALTITUDE_FROM_PRESSURE
  public::  CO2_SLICING_CLOUD_HEIGHT

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

  contains

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
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  integer:: Vza_Idx

  !--- intialize local variables using global variables
  Number_Of_Elements = Image%Number_Of_Elements

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

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  integer:: Num_Elem
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Line_Start
  integer:: Line_End
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  integer:: Ivza
  integer:: Level
  real:: Tc_Temp
  real:: Zc_Temp
  real:: U_Wnd_Temp
  real:: V_Wnd_Temp

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

  enddo Element_Loop_1
  enddo Line_Loop_1

end subroutine COMPUTE_CLOUD_TOP_LEVEL_NWP_WIND

!----------------------------------------------------------------------
! routine to interpolate pressure to flight level altitude.
!----------------------------------------------------------------------
subroutine COMPUTE_ALTITUDE_FROM_PRESSURE(Line_Idx_Min,Num_Lines)

!--- Based on Sarah Monette calculations for HS3.  Her calculation is based on:
!--- http://psas.pdx.edu/RocketScience/PressureAltitude_Derived.pdf.

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  integer:: Num_Elem
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Line_Start
  integer:: Line_End
  integer:: Nwp_Lon_Idx
  integer:: Nwp_Lat_Idx
  real:: Pc_Temp
  real:: Tc_Temp
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
     Alt_Acha(Elem_Idx,Line_Idx) = Alt_Temp

  enddo Element_Loop_1
  enddo Line_Loop_1

end subroutine COMPUTE_ALTITUDE_FROM_PRESSURE

!----------------------------------------------------------------------
! Compute CO2 Slicing
!----------------------------------------------------------------------
subroutine CO2_SLICING_CLOUD_HEIGHT(Num_Elem,Line_Idx_min,Num_Lines, &
                                    Pressure_Profile,Cloud_Type, &
                                    Pc_Cirrus_Co2,Tc_Cirrus_Co2)
  integer, intent(in):: Num_Elem
  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  integer(kind=int1), intent(in), dimension(:,:):: Cloud_Type
  real, intent(in), dimension(:):: Pressure_Profile
  real, intent(out), dimension(:,:):: Pc_Cirrus_Co2, Tc_Cirrus_Co2
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
  integer (kind=int4):: Box_Width
  real:: Zc_Cirrus_Co2
  real:: Count_Temporary, Sum_Temporary, Temperature_Temporary
  integer:: i1, i2, j1, j2
  integer:: Lev_Idx_Temp
  integer:: Pc_Lev_Idx


  logical, dimension(:,:), allocatable:: Mask
  real, dimension(:,:), allocatable:: Tc_Cirrus_Co2_Temp

  Line_Start = Line_Idx_Min
  Line_End = Line_Start + Num_Lines - 1

  !--- intialize output
  Pc_Cirrus_Co2 = Missing_Value_Real4
  Tc_Cirrus_Co2 = Missing_Value_Real4

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

     !--- only do this for appropriate cloud types
     if (Cloud_Type(Elem_Idx,Line_Idx) /= sym%CIRRUS_TYPE .and.   & 
!        Cloud_Type(Elem_Idx,Line_Idx) /= sym%PROB_CLEAR_TYPE .and.  &
!        Cloud_Type(Elem_Idx,Line_Idx) /= sym%CLEAR_TYPE .and.  &
         Cloud_Type(Elem_Idx,Line_Idx) /= sym%OVERLAP_TYPE) then
          cycle
     endif

     !--- compute cloud top pressure using each channel pair
     Beta_Target = 1.0

     call COMPUTE_BETA_PROFILE(ch(35)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(35)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(35)%Rad_BB_Cloud_Profile, &
                               ch(36)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(36)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(36)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               Pressure_Profile, &
                               Beta_Target, &
                               Pc_35_36,Pc_Lev_Idx)

     if (Pc_35_36 /= Missing_Value_Real4) then
         Pc_Cirrus_Co2(Elem_Idx,Line_Idx) = Pc_35_36
         cycle
     endif

     call COMPUTE_BETA_PROFILE(ch(34)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(34)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(34)%Rad_BB_Cloud_Profile, &
                               ch(35)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(35)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(35)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               Pressure_Profile, &
                               Beta_Target, &
                               Pc_34_35,Pc_Lev_Idx)

     if (Pc_34_35 /= Missing_Value_Real4) then
         Pc_Cirrus_Co2(Elem_Idx,Line_Idx) = Pc_34_35
         cycle
     endif

     call COMPUTE_BETA_PROFILE(ch(33)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(33)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(33)%Rad_BB_Cloud_Profile, &
                               ch(34)%Rad_Toa(Elem_Idx,Line_Idx), &
                               ch(34)%Rad_Toa_Clear(Elem_Idx,Line_Idx), &
                               rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(34)%Rad_BB_Cloud_Profile, &
                               Tropo_Level_Idx, &
                               Sfc_Level_Idx, &
                               Pressure_Profile, &
                               Beta_Target, &
                               Pc_33_34,Pc_Lev_Idx)

     if (Pc_33_34 /= Missing_Value_Real4) then
         Pc_Cirrus_Co2(Elem_Idx,Line_Idx) = Pc_33_34
         cycle
     endif
    
  enddo Element_Loop
  enddo Line_Loop

  !-------------------------------------------------------------------------
  ! determine temperature
  !-------------------------------------------------------------------------
  Line_Loop_2: do Line_Idx = Line_Start, Line_End
  Element_Loop_2: do Elem_Idx = 1, Num_Elem

     if (Pc_Cirrus_Co2(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- compute temperature
     call KNOWING_P_COMPUTE_T_Z_NWP(Nwp_Lon_Idx,Nwp_Lat_Idx, &
                                    Pc_Cirrus_Co2(Elem_Idx,Line_Idx), &
                                    Tc_Cirrus_Co2(Elem_Idx,Line_Idx), &
                                    Zc_Cirrus_Co2,Lev_Idx_Temp)

     !-- compute emissivity
     Ec_Cirrus_Co2(Elem_Idx,Line_Idx) = EMISSIVITY(ch(33)%Rad_Toa(Elem_Idx,Line_Idx),  &
                                                   ch(33)%Rad_Toa_Clear(Elem_Idx,Line_Idx),  &
                             rtm(Nwp_Lon_Idx,Nwp_Lat_Idx)%d(Vza_Rtm_Idx)%ch(33)%Rad_BB_Cloud_Profile(Pc_Lev_Idx))

  enddo Element_Loop_2
  enddo Line_Loop_2

  !-------------------------------------------------------------------------
  ! spatially interpolate
  !-------------------------------------------------------------------------

  !--- compute size of averaging window
  call COMPUTE_BOX_WIDTH(SOUNDER_RESOLUTION_KM,Box_Width_KM,Box_Width)

  !--- allocate temporary arrays needed here
  allocate(Mask(Num_Elem,Num_Lines))
  allocate(Tc_Cirrus_Co2_Temp(Num_Elem,Num_Lines))

  Tc_Cirrus_Co2_Temp = Missing_Value_Real4

  Mask = .false.

  !--- create mask which identifies pixels to be included in analysis
  where(Pc_Cirrus_Co2 /= Missing_Value_Real4 .and. &
        Pc_Cirrus_Co2 < PC_CIRRUS_MAX_THRESH .and. &
        Ec_Cirrus_Co2 > EC_CIRRUS_MIN_THRESH)
      Mask = .true.
  end where

  Element_Loop_3: do Elem_Idx = 1, Num_Elem, Box_Width/2

      i1 = min(Num_Elem,max(1,Elem_Idx - Box_Width))
      i2 = min(Num_Elem,max(1,Elem_Idx + Box_Width))

      Line_Loop_3: do Line_Idx = Line_Start, Line_End, Box_Width/2
          j1 = min(Num_Lines,max(1,Line_Idx - Box_Width))
          j2 = min(Num_Lines,max(1,Line_Idx + Box_Width))

          Count_Temporary = count(Mask(i1:i2,j1:j2))
          Sum_Temporary = sum(Tc_Cirrus_Co2(i1:i2,j1:j2),Mask(i1:i2,j1:j2))
          if (Count_Temporary > COUNT_MIN_TEMPERATURE_CIRRUS) then
              Temperature_Temporary = Sum_Temporary / Count_Temporary
          else
              Temperature_Temporary = Missing_Value_Real4
          endif

          if (Temperature_Temporary /= Missing_Value_Real4) then
             Tc_Cirrus_Co2_Temp(i1:i2,j1:j2) = Temperature_Temporary
          endif

      enddo Line_Loop_3
   enddo Element_Loop_3

   Tc_Cirrus_Co2 = Tc_Cirrus_Co2_Temp

  if (allocated(Mask)) deallocate(Mask)
  if (allocated(Tc_Cirrus_Co2_Temp)) deallocate(Tc_Cirrus_Co2_Temp)

end subroutine  CO2_SLICING_CLOUD_HEIGHT
!==================================================================================================
!
!==================================================================================================
subroutine COMPUTE_BETA_PROFILE(Ch_X_Rad_Toa, &
                                Ch_X_Rad_Toa_Clear, &
                                Ch_X_Rad_BB_Cloud_Profile, &
                                Ch_Y_Rad_Toa, &
                                Ch_Y_Rad_Toa_Clear, &
                                Ch_Y_Rad_BB_Cloud_Profile, &
                                Tropo_Level_Idx, &
                                Sfc_Level_Idx, &
                                Pressure_Profile, &
                                Beta_Target, &
                                Pc,  &
                                Pc_Lev_Idx)

   real, intent(in):: Ch_X_Rad_Toa
   real, intent(in):: Ch_X_Rad_Toa_Clear
   real, intent(in), dimension(:):: Ch_X_Rad_BB_Cloud_Profile
   real, intent(in):: Ch_Y_Rad_Toa
   real, intent(in):: Ch_Y_Rad_Toa_Clear
   real, intent(in), dimension(:):: Ch_Y_Rad_BB_Cloud_Profile
   integer, intent(in):: Tropo_Level_Idx
   integer, intent(in):: Sfc_Level_Idx
   real, intent(in), dimension(:):: Pressure_Profile
   real, intent(in):: Beta_Target
   real, intent(out):: Pc
   integer, intent(out):: Pc_Lev_Idx

   integer:: Lev_Idx
   integer:: Lev_Idx_Start
   integer:: Lev_Idx_End
   integer:: Num_Levels
   real:: Beta_X_Y
   real:: Beta_X_Y_Prev
   real:: Ch_X_Emissivity
   real:: Ch_Y_Emissivity

   !--- set levels for solution
   Num_Levels = size(Pressure_Profile)
   Lev_Idx_Start = Tropo_Level_Idx
   Lev_Idx_End = Sfc_Level_Idx

   !--- initialize
   Pc = Missing_Value_Real4
   Beta_X_Y_Prev = Missing_Value_Real4
   Pc_Lev_Idx = Missing_Value_Int4

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
         Pc_Lev_Idx = Lev_Idx
         exit
      endif

      if (Beta_X_Y >= Beta_Target .and. Beta_X_Y_Prev < Beta_Target .and. Beta_X_Y_Prev /= Missing_Value_Real4) then
         Pc_Lev_Idx = Lev_Idx
         exit
      endif

      Beta_X_Y_Prev = Beta_X_Y
      
   enddo

   !--- set output pressure if solution found
   if (Pc_Lev_Idx /= Missing_Value_Int4) then
      Pc = Pressure_Profile(Pc_Lev_Idx)
   endif

end subroutine COMPUTE_BETA_PROFILE
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
   !====================================================================
   function BETA_RATIO(Emiss_top, Emiss_bot) result(beta)
      real(kind=real4), intent(in) :: Emiss_top
      real(kind=real4), intent(in) :: Emiss_bot
      real(kind=real4) :: beta

      beta = Missing_Value_Real4

      if (Emiss_top > 0.0 .and. Emiss_top < 1.0 .and. &
          Emiss_bot > 0.0 .and. Emiss_bot < 1.0) then

         beta = alog(1.0 - Emiss_top)/alog(1.0 - Emiss_bot)
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
! End of Module
!----------------------------------------------------------------------

end module CLOUD_HEIGHT_ROUTINES
