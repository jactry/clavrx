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
  use CONSTANTS
  use CALIBRATION_CONSTANTS
  use NUMERICAL_ROUTINES

  implicit none
   
  public:: POPULATE_PLANCK_TABLES
  public:: PLANCK_RAD_FAST
  public:: PLANCK_TEMP_FAST
  public:: PLANCK_RAD
  public:: PLANCK_TEMP
  public:: COMPUTE_BT_ARRAY
  public:: CONVERT_RADIANCE

!-- planck tables arrays
  integer, parameter, private:: nplanck = 161
  real (kind=int4), parameter, private:: T_planck_min = 180.0,  &
                                         delta_T_planck = 1.0
  real(kind=int4), dimension(nplanck), save, private:: B20
  real(kind=int4), dimension(nplanck), save, private:: B21
  real(kind=int4), dimension(nplanck), save, private:: B22
  real(kind=int4), dimension(nplanck), save, private:: B23
  real(kind=int4), dimension(nplanck), save, private:: B24
  real(kind=int4), dimension(nplanck), save, private:: B25
  real(kind=int4), dimension(nplanck), save, private:: B27
  real(kind=int4), dimension(nplanck), save, private:: B28
  real(kind=int4), dimension(nplanck), save, private:: B29
  real(kind=int4), dimension(nplanck), save, private:: B30
  real(kind=int4), dimension(nplanck), save, private:: B31
  real(kind=int4), dimension(nplanck), save, private:: B32
  real(kind=int4), dimension(nplanck), save, private:: B33
  real(kind=int4), dimension(nplanck), save, private:: B34
  real(kind=int4), dimension(nplanck), save, private:: B35
  real(kind=int4), dimension(nplanck), save, private:: B36
  real(kind=int4), dimension(nplanck), save, private:: B37
  real(kind=int4), dimension(nplanck), save, private:: B38
  real(kind=int4), dimension(nplanck), save, private:: B42
  real(kind=int4), dimension(nplanck), save, private:: B43
  real(kind=int4), dimension(nplanck), save, private:: T_planck

  contains

!-------------------------------------------------------------------------
 subroutine COMPUTE_BT_ARRAY(bt,rad,ichan,missing)

 real(kind=real4), dimension(:,:), intent(in):: rad
 integer(kind=int4):: ichan
 real(kind=real4):: missing
 real(kind=real4), dimension(:,:), intent(out):: bt
 integer:: i, j
 integer:: nx, ny

  nx = size(rad,1)
  ny = size(rad,2)

  bt = missing

  do i = 1, nx
    do j = 1, ny
       if (rad(i,j) /= missing) then
          bt(i,j) = PLANCK_TEMP_FAST(ichan,rad(i,j))
       endif
    enddo
  enddo

 end subroutine
!-------------------------------------------------------------------------
! subroutine POPULATE_PLANCK_TABLES(a1_20,a2_20,nu_20,a1_31,a2_31,nu_31,a1_32,a2_32,nu_32)
!
! compute planck function tables
!
!-------------------------------------------------------------------------
 subroutine POPULATE_PLANCK_TABLES(a1_20,a2_20,nu_20, &
                                   a1_21,a2_21,nu_21, &
                                   a1_22,a2_22,nu_22, &
                                   a1_27,a2_27,nu_27, &
                                   a1_28,a2_28,nu_28, &
                                   a1_29,a2_29,nu_29, &
                                   a1_30,a2_30,nu_30, &
                                   a1_31,a2_31,nu_31, &
                                   a1_32,a2_32,nu_32, &
                                   a1_33,a2_33,nu_33, &
                                   a1_34,a2_34,nu_34, &
                                   a1_35,a2_35,nu_35, &
                                   a1_36,a2_36,nu_36, &
                                   a1_37,a2_37,nu_37, &
                                   a1_38,a2_38,nu_38, &
                                   a1_42,a2_42,nu_42, &
                                   a1_43,a2_43,nu_43)


  real, intent(in):: a1_20,a2_20,nu_20, &
                     a1_21,a2_21,nu_21, &
                     a1_22,a2_22,nu_22, &
                     a1_27,a2_27,nu_27, &
                     a1_28,a2_28,nu_28, &
                     a1_29,a2_29,nu_29, &
                     a1_30,a2_30,nu_30, &
                     a1_31,a2_31,nu_31, &
                     a1_32,a2_32,nu_32, &
                     a1_33,a2_33,nu_33, &
                     a1_34,a2_34,nu_34, &
                     a1_35,a2_35,nu_35, &
                     a1_36,a2_36,nu_36, &
                     a1_37,a2_37,nu_37, &
                     a1_38,a2_38,nu_38, &
                     a1_42,a2_42,nu_42, &
                     a1_43,a2_43,nu_43
  integer:: i

  do i = 1, nplanck
    T_planck(i) = T_planck_min + (i-1)*delta_T_planck
    B20(i) = c1*(nu_20**3)/(exp((c2*nu_20)/ &
              ((T_planck(i)-a1_20)/a2_20))-1.0)
    B21(i) = c1*(nu_21**3)/(exp((c2*nu_21)/ &
              ((T_planck(i)-a1_21)/a2_21))-1.0)
    B22(i) = c1*(nu_22**3)/(exp((c2*nu_22)/ &
              ((T_planck(i)-a1_22)/a2_22))-1.0)
    B23(i) = c1*(nu_23**3)/(exp((c2*nu_23)/ &
              ((T_planck(i)-a1_23)/a2_23))-1.0)
    B24(i) = c1*(nu_24**3)/(exp((c2*nu_24)/ &
              ((T_planck(i)-a1_24)/a2_24))-1.0)
    B25(i) = c1*(nu_25**3)/(exp((c2*nu_25)/ &
              ((T_planck(i)-a1_25)/a2_25))-1.0)
    B27(i) = c1*(nu_27**3)/(exp((c2*nu_27)/ &
              ((T_planck(i)-a1_27)/a2_27))-1.0)
    B28(i) = c1*(nu_28**3)/(exp((c2*nu_28)/ &
              ((T_planck(i)-a1_28)/a2_28))-1.0)
    B29(i) = c1*(nu_29**3)/(exp((c2*nu_29)/ &
              ((T_planck(i)-a1_29)/a2_29))-1.0)
    B30(i) = c1*(nu_30**3)/(exp((c2*nu_30)/ &
              ((T_planck(i)-a1_30)/a2_30))-1.0)
    B31(i) = c1*(nu_31**3)/(exp((c2*nu_31)/ &
              ((T_planck(i)-a1_31)/a2_31))-1.0)
    B32(i) = c1*(nu_32**3)/(exp((c2*nu_32)/ &
              ((T_planck(i)-a1_32)/a2_32))-1.0)
    B33(i) = c1*(nu_33**3)/(exp((c2*nu_33)/ &
              ((T_planck(i)-a1_33)/a2_33))-1.0)
    B34(i) = c1*(nu_34**3)/(exp((c2*nu_34)/ &
              ((T_planck(i)-a1_34)/a2_34))-1.0)
    B35(i) = c1*(nu_35**3)/(exp((c2*nu_35)/ &
              ((T_planck(i)-a1_35)/a2_35))-1.0)
    B36(i) = c1*(nu_36**3)/(exp((c2*nu_36)/ &
              ((T_planck(i)-a1_36)/a2_36))-1.0)
    B37(i) = c1*(nu_37**3)/(exp((c2*nu_37)/ &
              ((T_planck(i)-a1_37)/a2_37))-1.0)
    B38(i) = c1*(nu_38**3)/(exp((c2*nu_38)/ &
              ((T_planck(i)-a1_38)/a2_38))-1.0)
    B42(i) = c1*(nu_42**3)/(exp((c2*nu_42)/ &
              ((T_planck(i)-a1_42)/a2_42))-1.0)
    B43(i) = c1*(nu_43**3)/(exp((c2*nu_43)/ &
              ((T_planck(i)-a1_43)/a2_43))-1.0)
  enddo
  

  end subroutine POPULATE_PLANCK_TABLES

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

!--- compute planck emission for cloud temperature
     l = (T - T_planck_min)/delta_T_planck
     l = max(1,min(nplanck-1,l))

    if (ichan == 20) then 
      dB_dT_tmp = (B20(l+1)-B20(l))/(T_planck(l+1)-T_planck(l))
      B = B20(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 21) then 
      dB_dT_tmp = (B21(l+1)-B21(l))/(T_planck(l+1)-T_planck(l))
      B = B21(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 22) then 
      dB_dT_tmp = (B22(l+1)-B22(l))/(T_planck(l+1)-T_planck(l))
      B = B22(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 23) then 
      dB_dT_tmp = (B23(l+1)-B23(l))/(T_planck(l+1)-T_planck(l))
      B = B23(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 24) then 
      dB_dT_tmp = (B24(l+1)-B24(l))/(T_planck(l+1)-T_planck(l))
      B = B24(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 25) then 
      dB_dT_tmp = (B25(l+1)-B24(l))/(T_planck(l+1)-T_planck(l))
      B = B25(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 27) then 
      dB_dT_tmp = (B27(l+1)-B27(l))/(T_planck(l+1)-T_planck(l))
      B = B27(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 28) then 
      dB_dT_tmp = (B28(l+1)-B28(l))/(T_planck(l+1)-T_planck(l))
      B = B28(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 29) then 
      dB_dT_tmp = (B29(l+1)-B29(l))/(T_planck(l+1)-T_planck(l))
      B = B29(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 30) then 
      dB_dT_tmp = (B30(l+1)-B30(l))/(T_planck(l+1)-T_planck(l))
      B = B30(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 31) then 
      dB_dT_tmp = (B31(l+1)-B31(l))/(T_planck(l+1)-T_planck(l))
      B = B31(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 32) then 
      dB_dT_tmp = (B32(l+1)-B32(l))/(T_planck(l+1)-T_planck(l))
      B = B32(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 33) then 
      dB_dT_tmp = (B33(l+1)-B33(l))/(T_planck(l+1)-T_planck(l))
      B = B33(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 34) then 
      dB_dT_tmp = (B34(l+1)-B34(l))/(T_planck(l+1)-T_planck(l))
      B = B34(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 35) then 
      dB_dT_tmp = (B35(l+1)-B35(l))/(T_planck(l+1)-T_planck(l))
      B = B35(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 36) then 
      dB_dT_tmp = (B36(l+1)-B36(l))/(T_planck(l+1)-T_planck(l))
      B = B36(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 37) then 
      dB_dT_tmp = (B37(l+1)-B37(l))/(T_planck(l+1)-T_planck(l))
      B = B37(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 38) then 
      dB_dT_tmp = (B38(l+1)-B38(l))/(T_planck(l+1)-T_planck(l))
      B = B38(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 42) then 
      dB_dT_tmp = (B42(l+1)-B42(l))/(T_planck(l+1)-T_planck(l))
      B = B42(l) + (T - T_planck(l)) * (dB_dT_tmp)
    elseif (ichan == 43) then 
      dB_dT_tmp = (B43(l+1)-B43(l))/(T_planck(l+1)-T_planck(l))
      B = B43(l) + (T - T_planck(l)) * (dB_dT_tmp)
    else
      print *, "unsupported channel number in Planck Computation, stopping"
      stop
    endif

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
    real (kind=real4) :: T, dB_dT_tmp
    integer:: l

    if (ichan == 20) then 
      call locate(B20,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B20(l+1)-B20(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B20(l)) / (dB_dT_tmp)
    elseif (ichan == 21) then 
      call locate(B21,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B21(l+1)-B21(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B21(l)) / (dB_dT_tmp)
    elseif (ichan == 22) then 
      call locate(B22,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B22(l+1)-B22(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B22(l)) / (dB_dT_tmp)
    elseif (ichan == 23) then 
      call locate(B23,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B23(l+1)-B23(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B23(l)) / (dB_dT_tmp)
    elseif (ichan == 24) then 
      call locate(B24,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B24(l+1)-B24(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B24(l)) / (dB_dT_tmp)
    elseif (ichan == 25) then 
      call locate(B25,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B25(l+1)-B25(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B25(l)) / (dB_dT_tmp)
    elseif (ichan == 27) then 
      call locate(B27,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B27(l+1)-B27(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B27(l)) / (dB_dT_tmp)
    elseif (ichan == 28) then 
      call locate(B28,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B28(l+1)-B28(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B28(l)) / (dB_dT_tmp)
    elseif (ichan == 29) then 
      call locate(B29,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B29(l+1)-B29(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B29(l)) / (dB_dT_tmp)
    elseif (ichan == 30) then 
      call locate(B30,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B30(l+1)-B30(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B30(l)) / (dB_dT_tmp)
    elseif (ichan == 31) then 
      call locate(B31,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B31(l+1)-B31(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B31(l)) / (dB_dT_tmp)
    elseif (ichan == 32) then 
      call locate(B32,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B32(l+1)-B32(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B32(l)) / (dB_dT_tmp)
    elseif (ichan == 33) then 
      call locate(B33,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B33(l+1)-B33(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B33(l)) / (dB_dT_tmp)
    elseif (ichan == 34) then 
      call locate(B34,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B34(l+1)-B34(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B34(l)) / (dB_dT_tmp)
    elseif (ichan == 35) then 
      call locate(B35,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B35(l+1)-B35(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B35(l)) / (dB_dT_tmp)
    elseif (ichan == 36) then 
      call locate(B36,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B36(l+1)-B36(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B36(l)) / (dB_dT_tmp)
    elseif (ichan == 37) then 
      call locate(B37,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B37(l+1)-B37(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B37(l)) / (dB_dT_tmp)
    elseif (ichan == 38) then 
      call locate(B38,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B38(l+1)-B38(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B38(l)) / (dB_dT_tmp)
    elseif (ichan == 42) then 
      call locate(B42,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B42(l+1)-B42(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B42(l)) / (dB_dT_tmp)
    elseif (ichan == 43) then 
      call locate(B43,nplanck,B,l)
      l = max(1,min(nplanck-1,l))
      dB_dT_tmp = (B43(l+1)-B43(l))/(T_planck(l+1)-T_planck(l))
      T = T_planck(l) + (B - B43(l)) / (dB_dT_tmp)
    else
      print *, "unsupported channel number in Planck Computation, stopping"
      stop
    endif

    if (present(dB_dT)) dB_dT = dB_dT_tmp

    return

  end function PLANCK_TEMP_FAST

!----------------------------------------------------------------------
! function PLANCK_RAD(ichan, T) result(B)
!----------------------------------------------------------------------

  function PLANCK_RAD(ichan, T) result(B)
    integer (kind=int4), intent(in) :: ichan
    real (kind=real4), intent(in) :: T 
    real (kind=real4) :: B

    if (ichan == 20) then 
       B = c1*(nu_20**3)/(exp((c2*nu_20)/((T-a1_20)/a2_20))-1.0)
    elseif (ichan == 21) then 
       B = c1*(nu_21**3)/(exp((c2*nu_21)/((T-a1_21)/a2_21))-1.0)
    elseif (ichan == 22) then 
       B = c1*(nu_22**3)/(exp((c2*nu_22)/((T-a1_22)/a2_22))-1.0)
    elseif (ichan == 23) then 
       B = c1*(nu_23**3)/(exp((c2*nu_23)/((T-a1_23)/a2_23))-1.0)
    elseif (ichan == 24) then 
       B = c1*(nu_24**3)/(exp((c2*nu_24)/((T-a1_24)/a2_24))-1.0)
    elseif (ichan == 25) then 
       B = c1*(nu_25**3)/(exp((c2*nu_25)/((T-a1_25)/a2_25))-1.0)
    elseif (ichan == 27) then 
       B = c1*(nu_27**3)/(exp((c2*nu_27)/((T-a1_27)/a2_27))-1.0)
    elseif (ichan == 28) then 
       B = c1*(nu_28**3)/(exp((c2*nu_28)/((T-a1_28)/a2_28))-1.0)
    elseif (ichan == 29) then 
       B = c1*(nu_29**3)/(exp((c2*nu_29)/((T-a1_29)/a2_29))-1.0)
    elseif (ichan == 30) then 
       B = c1*(nu_30**3)/(exp((c2*nu_30)/((T-a1_30)/a2_30))-1.0)
    elseif (ichan == 31) then 
       B = c1*(nu_31**3)/(exp((c2*nu_31)/((T-a1_31)/a2_31))-1.0)
    elseif (ichan == 32) then 
       B = c1*(nu_32**3)/(exp((c2*nu_32)/((T-a1_32)/a2_32))-1.0)
    elseif (ichan == 33) then 
       B = c1*(nu_33**3)/(exp((c2*nu_33)/((T-a1_33)/a2_33))-1.0)
    elseif (ichan == 34) then 
       B = c1*(nu_34**3)/(exp((c2*nu_34)/((T-a1_34)/a2_34))-1.0)
    elseif (ichan == 35) then 
       B = c1*(nu_35**3)/(exp((c2*nu_35)/((T-a1_35)/a2_35))-1.0)
    elseif (ichan == 36) then 
       B = c1*(nu_36**3)/(exp((c2*nu_36)/((T-a1_36)/a2_36))-1.0)
    elseif (ichan == 37) then 
       B = c1*(nu_37**3)/(exp((c2*nu_37)/((T-a1_37)/a2_37))-1.0)
    elseif (ichan == 38) then 
       B = c1*(nu_38**3)/(exp((c2*nu_38)/((T-a1_38)/a2_38))-1.0)
    elseif (ichan == 42) then 
       B = c1*(nu_42**3)/(exp((c2*nu_42)/((T-a1_42)/a2_42))-1.0)
    elseif (ichan == 43) then 
       B = c1*(nu_43**3)/(exp((c2*nu_43)/((T-a1_43)/a2_43))-1.0)
    else
      print *, "unsupported channel number in PLANCK_RAD, stopping"
      stop
    endif

    return

  end function PLANCK_RAD
  
!----------------------------------------------------------------------
! function PLANCK_TEMP(ichan, B) result(T) 
!----------------------------------------------------------------------
  function PLANCK_TEMP(ichan, B) result(T)
    integer (kind=int4), intent(in) :: ichan
    real (kind=real4), intent(in) :: B 
    real (kind=real4) :: T

    if (ichan == 20) then 
       T = a1_20 + a2_20 * ((c2*nu_20) / log( 1.0 + (c1*(nu_20**3))/B))
    elseif (ichan == 21) then 
       T = a1_21 + a2_21 * ((c2*nu_21) / log( 1.0 + (c1*(nu_21**3))/B))
    elseif (ichan == 22) then 
       T = a1_22 + a2_22 * ((c2*nu_22) / log( 1.0 + (c1*(nu_22**3))/B))
    elseif (ichan == 23) then 
       T = a1_23 + a2_23 * ((c2*nu_23) / log( 1.0 + (c1*(nu_23**3))/B))
    elseif (ichan == 24) then 
       T = a1_24 + a2_24 * ((c2*nu_24) / log( 1.0 + (c1*(nu_24**3))/B))
    elseif (ichan == 25) then 
       T = a1_25 + a2_25 * ((c2*nu_25) / log( 1.0 + (c1*(nu_25**3))/B))
    elseif (ichan == 27) then 
       T = a1_27 + a2_27 * ((c2*nu_27) / log( 1.0 + (c1*(nu_27**3))/B))
    elseif (ichan == 28) then 
       T = a1_28 + a2_28 * ((c2*nu_28) / log( 1.0 + (c1*(nu_28**3))/B))
    elseif (ichan == 29) then 
       T = a1_29 + a2_29 * ((c2*nu_29) / log( 1.0 + (c1*(nu_29**3))/B))
    elseif (ichan == 30) then 
       T = a1_30 + a2_30 * ((c2*nu_30) / log( 1.0 + (c1*(nu_30**3))/B))
    elseif (ichan == 31) then 
       T = a1_31 + a2_31 * ((c2*nu_31) / log( 1.0 + (c1*(nu_31**3))/B))
    elseif (ichan == 32) then 
       T = a1_32 + a2_32 * ((c2*nu_32) / log( 1.0 + (c1*(nu_32**3))/B))
    elseif (ichan == 33) then 
       T = a1_33 + a2_33 * ((c2*nu_33) / log( 1.0 + (c1*(nu_33**3))/B))
    elseif (ichan == 34) then 
       T = a1_34 + a2_34 * ((c2*nu_34) / log( 1.0 + (c1*(nu_34**3))/B))
    elseif (ichan == 35) then 
       T = a1_35 + a2_35 * ((c2*nu_35) / log( 1.0 + (c1*(nu_35**3))/B))
    elseif (ichan == 36) then 
       T = a1_36 + a2_36 * ((c2*nu_36) / log( 1.0 + (c1*(nu_36**3))/B))
    elseif (ichan == 37) then 
       T = a1_37 + a2_37 * ((c2*nu_37) / log( 1.0 + (c1*(nu_37**3))/B))
    elseif (ichan == 38) then 
       T = a1_38 + a2_38 * ((c2*nu_38) / log( 1.0 + (c1*(nu_38**3))/B))
    elseif (ichan == 40) then 
       T = a1_42 + a2_42 * ((c2*nu_42) / log( 1.0 + (c1*(nu_42**3))/B))
    elseif (ichan == 41) then 
       T = a1_43 + a2_43 * ((c2*nu_43) / log( 1.0 + (c1*(nu_43**3))/B))
    else
      print *, "unsupported channel number in PLANCK_TEMP, stopping"
      stop
    endif

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
