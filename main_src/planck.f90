!$Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: planck.f90 (src)
!       PLANCK (program)
!
! PURPOSE: this module holds the routine to do rapid Planck
!          computations using a table lookup approach. 
!
! DESCRIPTION: This has been shown to speed up CLAVR-x over using the
!              explicit planck function with expontentials
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
! REVISION HISTORY:
!   August 31, 2006 - Created
!
! ROUTINES:
! public::
!   POPULATE_PLANCK_TABLES
!   PLANCK_RAD_FAST
!   PLANCK_TEMP_FAST
!   PLANCK_RAD
!   PLANCK_TEMP
!
!  Modis Band   Avhrr Band   Abi Band   Wavelength  
!     01            1          2           0.63
!     02            2          3           0.86
!     06            3a         5           1.60
!     20            3b         7           3.75
!     21            -          7           3.90
!     22            -          7           3.90
!     23            -          -           3.95
!     24            -          -           4.05
!     25            -          -           4.45
!     27            -          9           6.70
!     28            -         10           7.30
!     29            -         11           8.50
!     30            -         12           9.70
!     31            4         14           11.00
!     32            5         15           12.00
!     33            -         16           13.30
!     34            -          -           13.30
!     35            -          -           13.30
!     36            -          -           13.30
!     40            -          -            3.75
!     41            -          -           11.45
!
!--------------------------------------------------------------------------------------
 module PLANCK
  use CX_CONSTANTS_MOD
  use CALIBRATION_CONSTANTS
  use NUMERICAL_TOOLS_MOD
  use PIXEL_COMMON, only: &
        Sensor, &
        Ch

  implicit none
  private
  public:: POPULATE_PLANCK_TABLES
  public:: POPULATE_PLANCK_TABLES_SOUNDER
  public:: PLANCK_RAD_FAST
  public:: PLANCK_TEMP_FAST
  public:: PLANCK_RAD
  public:: PLANCK_TEMP
  public:: COMPUTE_BT_ARRAY
  public:: COMPUTE_BT_ARRAY_SOUNDER
  public:: CONVERT_RADIANCE
  private:: PLANCK_TEMP_FAST_SOUNDER

!-- planck tables arrays
  integer, parameter, private:: Nplanck = 161
  real (kind=int4), parameter, private:: T_Planck_min = 180.0,  &
                                         delta_T_Planck = 1.0
  real(kind=int4), dimension(20:Nchan_Clavrx,Nplanck), save, private:: BB_Rad = Missing_Value_Real4
  real(kind=int4), dimension(Nplanck), save, private:: BB_Rad_375um_Sndr = Missing_Value_Real4
  real(kind=int4), dimension(Nplanck), save, private:: BB_Rad_11um_Sndr = Missing_Value_Real4
  real(kind=int4), dimension(Nplanck), save, private:: BB_Rad_12um_Sndr = Missing_Value_Real4
  real(kind=int4), dimension(Nplanck), save, private:: T_Planck = Missing_Value_Real4

  contains

!-------------------------------------------------------------------------
 subroutine COMPUTE_BT_ARRAY(bt,rad,ichan,missing)

 real(kind=real4), dimension(:,:), intent(in):: Rad
 integer(kind=int4):: ichan
 real(kind=real4):: missing
 real(kind=real4), dimension(:,:), intent(out):: Bt
 integer:: i, j
 integer:: nx, ny

  nx = size(Rad,1)
  ny = size(Rad,2)

  Bt = missing

  do i = 1, nx
    do j = 1, ny
       if (Rad(i,j) /= missing) then
          Bt(i,j) = PLANCK_TEMP_FAST(ichan,Rad(i,j))
       endif
    enddo
  enddo

 end subroutine COMPUTE_BT_ARRAY

!-------------------------------------------------------------------------
 subroutine COMPUTE_BT_ARRAY_SOUNDER(bt,rad,ichan,missing)

 real(kind=real4), dimension(:,:), intent(in):: Rad
 integer(kind=int4):: ichan
 real(kind=real4):: missing
 real(kind=real4), dimension(:,:), intent(out):: Bt
 integer:: i, j
 integer:: nx, ny

  nx = size(Rad,1)
  ny = size(Rad,2)

  Bt = missing

  do i = 1, nx
    do j = 1, ny
       if (Rad(i,j) /= missing) then
          Bt(i,j) = PLANCK_TEMP_FAST_SOUNDER(ichan,Rad(i,j))
       endif
    enddo
  enddo

 end subroutine COMPUTE_BT_ARRAY_SOUNDER

!-------------------------------------------------------------------------
! subroutine POPULATE_PLANCK_TABLES()
!
! compute planck function tables
!
!-------------------------------------------------------------------------
 subroutine POPULATE_PLANCK_TABLES()

  integer:: i, ichan

  do i = 1, Nplanck
    T_Planck(i) = T_Planck_min + (i-1)*delta_T_Planck
  enddo

  do ichan = 20,Nchan_Clavrx
    if (Ch(ichan)%Obs_Type /= THERMAL_OBS_TYPE .and. &
        Ch(ichan)%Obs_Type /= MIXED_OBS_TYPE) cycle
     if (Sensor%Chan_On_Flag_Default(ichan)==sym%YES) then
         BB_Rad(ichan,:) = c1*(Planck_Nu(ichan)**3)/ &
               (exp((c2*Planck_Nu(ichan))/((T_Planck-Planck_A1(ichan))/Planck_A2(ichan)))-1.0)
     endif
  enddo

  end subroutine POPULATE_PLANCK_TABLES

!-------------------------------------------------------------------------
! subroutine POPULATE_PLANCK_TABLES_SOUNDER()
!
! compute planck function tables for IFF HIRS and CrIS ch
!
!-------------------------------------------------------------------------
 subroutine POPULATE_PLANCK_TABLES_SOUNDER()

     !--- 3.75um only in HIRS, CrIS doesn't have it
     if (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then
        BB_Rad_375um_Sndr = c1*(Planck_Nu_375um_Sndr**3)/ &
              (exp((c2*Planck_Nu_375um_Sndr)/ &
              ((T_Planck-Planck_A1_375um_Sndr)/Planck_A2_375um_Sndr))-1.0)
     endif

     BB_Rad_11um_Sndr = c1*(Planck_Nu_11um_Sndr**3)/ &
            (exp((c2*Planck_Nu_11um_Sndr)/ &
            ((T_Planck-Planck_A1_11um_Sndr)/Planck_A2_11um_Sndr))-1.0)

     BB_Rad_12um_Sndr = c1*(Planck_Nu_12um_Sndr**3)/ &
            (exp((c2*Planck_Nu_12um_Sndr)/ &
            ((T_Planck-Planck_A1_12um_Sndr)/Planck_A2_12um_Sndr))-1.0)


  end subroutine POPULATE_PLANCK_TABLES_SOUNDER

!------------------------------------------------------------------
! function PLANCK_RAD_FAST(ichan, T, dB_dT) result(B)
!
! Subroutine to convert brightness temperature to radiance using a
! look-up table *Function that returns a scalar*.
!------------------------------------------------------------------
  function PLANCK_RAD_FAST(ichan, T, dB_dT) result(B)
    integer (kind=int4), intent(in) :: ichan
    real (kind=real4), intent(in) :: T 
    real (kind=real4), optional, intent(out) :: dB_dT
    real (kind=real4) :: B, dB_dT_tmp
    integer:: l

    !--- check for appropriate channel
    if (Ch(ichan)%Obs_Type /= THERMAL_OBS_TYPE .and. &
        Ch(ichan)%Obs_Type /= MIXED_OBS_TYPE) return
    if (Sensor%Chan_On_Flag_Default(ichan) == sym%NO) then 
      print *, "unsupported channel number ",ichan," in Fast Planck Rad Computation, stopping"
      stop
    endif

    !--- compute planck emission for cloud temperature
    l = (T - T_Planck_min)/delta_T_Planck
    l = max(1,min(Nplanck-1,l))

    dB_dT_tmp = (BB_Rad(ichan,l+1)-BB_Rad(ichan,l))/(T_Planck(l+1)-T_Planck(l))
    B = BB_Rad(ichan,l) + (T - T_Planck(l)) * (dB_dT_tmp)

    if (present(dB_dT)) dB_dT = dB_dT_tmp

    return

  end function PLANCK_RAD_FAST

!------------------------------------------------------------------------
! function PLANCK_TEMP_FAST(ichan, B, dB_dT) result(T)
!
! Subroutine to convert radiance (B) to brightness temperature(T) using a
! look-up table *Function that returns a scalar*.
!------------------------------------------------------------------
  function PLANCK_TEMP_FAST(ichan, B, dB_dT) result(T)
    integer (kind=int4), intent(in) :: ichan
    real (kind=real4), intent(in) :: B 
    real (kind=real4), optional, intent(out) :: dB_dT
    real (kind=real4) :: T
    real (kind=real4) :: dB_dT_Tmp
    integer:: l

    T = Missing_Value_Real4
    dB_dT_Tmp = Missing_Value_Real4

    !--- check for appropriate channel
    if (Ch(ichan)%Obs_Type /= THERMAL_OBS_TYPE .and. &
        Ch(ichan)%Obs_Type /= MIXED_OBS_TYPE) return

    if (Sensor%Chan_On_Flag_Default(ichan) == sym%NO) then
      print *, "unsupported channel number ",ichan," in Fast Planck Temp Computation, stopping"
      stop
    endif

    !---- compute brightness temperature
    call LOCATE(BB_Rad(ichan,:),Nplanck,B,l)
    l = max(1,min(Nplanck-1,l))
    dB_dT_tmp = (BB_Rad(ichan,l+1)-BB_Rad(ichan,l))/(T_Planck(l+1)-T_Planck(l))
    T = T_Planck(l) + (B - BB_Rad(ichan,l)) / (dB_dT_tmp)

    if (present(dB_dT)) dB_dT = dB_dT_tmp

    return

  end function PLANCK_TEMP_FAST

!------------------------------------------------------------------------
! function PLANCK_TEMP_FAST_SOUNDER(ichan, B, dB_dT) result(T)
!
! Note: This function only for IFF HIRS and CrIS 11 and 12um ch.
! Subroutine to convert radiance (B) to brightness temperature(T) using a
! look-up table *Function that returns a scalar*.
!------------------------------------------------------------------
  function PLANCK_TEMP_FAST_SOUNDER(ichan, B) result(T)
    integer (kind=int4), intent(in) :: ichan
    real (kind=real4), intent(in) :: B
    real(kind=int4), dimension(Nplanck) :: BB_Rad_Tmp
    real (kind=real4) :: T
    real (kind=real4) :: dB_dT_Tmp
    integer:: l

    T = Missing_Value_Real4
    dB_dT_Tmp = Missing_Value_Real4

    !---- depending on ch # make temporate BB_Rad
    if (ichan == 20) then
       BB_Rad_Tmp = BB_Rad_375um_Sndr
    elseif (ichan == 31) then
       BB_Rad_Tmp = BB_Rad_11um_Sndr
    elseif (ichan == 32) then
       BB_Rad_Tmp = BB_Rad_12um_Sndr
    else
       print *, "Unsupported channel number ",ichan," in Fast Planck Temp Sounder Computation"
       stop
    endif

    !---- compute brightness temperature
    call LOCATE(BB_Rad_Tmp,Nplanck,B,l)
    l = max(1,min(Nplanck-1,l))
    dB_dT_tmp = (BB_Rad_Tmp(l+1)-BB_Rad_Tmp(l))/(T_Planck(l+1)-T_Planck(l))
    T = T_Planck(l) + (B - BB_Rad_Tmp(l)) / (dB_dT_tmp)

    return

  end function PLANCK_TEMP_FAST_SOUNDER

!----------------------------------------------------------------------
! function PLANCK_RAD(ichan, T) result(B)
!----------------------------------------------------------------------

  function PLANCK_RAD(ichan, T) result(B)
    integer (kind=int4), intent(in) :: ichan
    real (kind=real4), intent(in) :: T 
    real (kind=real4) :: B

    !--- check for appropriate channel
    if (Ch(ichan)%Obs_Type /= THERMAL_OBS_TYPE .and. &
        Ch(ichan)%Obs_Type /= MIXED_OBS_TYPE) return
    if (Sensor%Chan_On_Flag_Default(ichan) == sym%NO) then
      print *, "unsupported channel number ",ichan," in Planck Rad Computation, stopping"
      stop
    endif

    B = c1*(Planck_Nu(ichan)**3)/ &
        (exp((c2*Planck_Nu(ichan))/ &
        ((T-Planck_A1(ichan))/Planck_A2(ichan)))-1.0)

    return

  end function PLANCK_RAD
  
!----------------------------------------------------------------------
! function PLANCK_TEMP(ichan, B) result(T) 
!----------------------------------------------------------------------
  function PLANCK_TEMP(ichan, B) result(T)
    integer (kind=int4), intent(in) :: ichan
    real (kind=real4), intent(in) :: B 
    real (kind=real4) :: T

    !--- check for appropriate channel
    if (Ch(ichan)%Obs_Type /= THERMAL_OBS_TYPE .and. &
        Ch(ichan)%Obs_Type /= MIXED_OBS_TYPE) return
    if (Sensor%Chan_On_Flag_Default(ichan) == sym%NO) then
      print *, "unsupported channel number ",ichan," in Planck Temp Computation, stopping"
      stop
    endif

    T = Planck_A1(ichan) + Planck_A2(ichan) * ((c2*Planck_Nu(ichan)) / log( 1.0 + (c1*(Planck_Nu(ichan)**3))/B))

    return

  end function PLANCK_TEMP

 !----------------------------------------------------------------------
 ! Convert radiances based on channel centroid wavenumber (nu)
 ! from NASA standad (W m-2 um-1 sr-1) 
 ! to NOAA standard (mW/cm^2/cm^-1/str)
 !
 !----------------------------------------------------------------------
 subroutine CONVERT_RADIANCE(Radiance,Nu,Missing_Value)
  real (kind=real4), dimension(:,:), intent(inout):: Radiance
  real (kind=real4), intent(in):: Nu
  real (kind=real4), intent(in):: Missing_Value
 
  where(Radiance /= Missing_Value)
       Radiance = Radiance * (((10000.0 / Nu )**2) / 10.0)
  end where

  return

 end subroutine CONVERT_RADIANCE

end module PLANCK
