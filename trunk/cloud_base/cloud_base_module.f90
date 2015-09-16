! $Id: awg_cloud_height.f90 583 2014-10-08 03:43:36Z heidinger $
module CLOUD_BASE
!---------------------------------------------------------------------
! This module houses the routines associated with...
!
! Author: 
!
! Reference:
!
!----------------------------------------------------------------------
  use CLOUD_BASE_SERVICES
  use ACHA_SERVICES_MOD, only : LOCATE

  implicit none

  public:: CLOUD_BASE_ALGORITHM

  private:: INTERPOLATE_PROFILE_LOCAL
  private:: NULL_PIX_POINTERS 
  private:: KNOWING_Z_COMPUTE_T_P 

  !--- include the non-system specific variables
  include 'cloud_base_parameters.inc'

  !--- interpoLated profiles
  real, private, dimension(Num_Levels_RTM_Prof) :: Temp_Prof_RTM
  real, private, dimension(Num_Levels_RTM_Prof) :: Press_Prof_RTM
  real, private, dimension(Num_Levels_RTM_Prof) :: Hght_Prof_RTM
  integer, private:: Sfc_Level_RTM
  integer, private:: Tropo_Level_RTM

  real, private, PARAMETER:: MISSING_VALUE_REAL = -999.0
  integer, private, PARAMETER:: MISSING_VALUE_INTEGER = -999

  contains 

!------------------------------------------------------------------------------
! Cloud BaseE Height Algorithm 
!
! Author: Andrew Heidinger, NOAA
!
! Assumptions
!
! Limitations
!
! NOTE.  This algorithm use the same input and output structures as 
!        the AWG_CLOUD_HEIGHT_ALGORITHM.
!        Do not overwrite elements of the Output structure expect those
!        generated here.
!
!      Output%Tau
!      Output%Ec
!      Output%Reff
!      Output%Zc_Top
!      Output%Zc_Base
!      Output%Pc_Top
!      Output%Pc_Base
!
!----------------------------------------------------------------------
! modification history
!
!------------------------------------------------------------------------------
  subroutine CLOUD_BASE_ALGORITHM(Input, Symbol, Output)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(symbol_acha), intent(inout) :: Symbol
  type(acha_input_struct), intent(inout) :: Input
  type(acha_output_struct), intent(inout) :: Output

  !===============================================================================
  !  Pixel level RTM structure
  !===============================================================================
 
  type(acha_rtm_nwp_struct) :: RTM_NWP

  !===============================================================================
  !  Local Variable Declaration
  !===============================================================================

  integer (kind=int4):: Smooth_Nwp_Fields_Flag_Temp
  integer:: Elem_Idx
  integer:: Line_Idx
  real:: Inwp_Weight
  real:: Jnwp_Weight

  integer:: Cloud_Type

  !--- scalar local variables
  real (kind=real4):: Cloud_Extinction
  real (kind=real4):: Cloud_Geometrical_Thickness
  real (kind=real4):: Cloud_Geometrical_Thickness_Top_Offset
  real (kind=real4):: Cloud_Geometrical_Thickness_eff
  real (kind=real4):: Zc_Top_Max
  real (kind=real4):: Zc_Base_Min
  real (kind=real4):: R4_Dummy

  integer (kind=int4):: Itemp
  integer (kind=int4):: Ilev

!-----------------------------------------------------------------------
! BEGIN EXECUTABLE CODE
!-----------------------------------------------------------------------

   !-------------------------------------------------------------------------
   ! Initialization
   !-------------------------------------------------------------------------
   Output%Zc_Top = MISSING_VALUE_REAL
   Output%Zc_Base = MISSING_VALUE_REAL
   Output%Pc_Top = MISSING_VALUE_REAL
   Output%Pc_Base = MISSING_VALUE_REAL

   !--------------------------------------------------------------------------
   ! loop over pixels in scanlines
   !--------------------------------------------------------------------------
   Line_loop: do Line_Idx = Line_Idx_min,Input%Number_of_Lines + Line_Idx_min - 1

    Element_Loop:   do Elem_Idx = 1, Input%Number_of_Elements

     if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) cycle
     if ((Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx) <= 0) .or. &
         (Input%Line_Idx_Nwp(Elem_Idx,Line_Idx) <= 0)) then
         cycle
    endif

    !---- null profile pointers each time - WCS3
    call NULL_PIX_POINTERS(Input, RTM_NWP)

    !--- for convenience, save nwp indices to local variables
    Inwp_Weight = Input%Longitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
    Jnwp_Weight = Input%Latitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
    Cloud_Type = Input%Cloud_Type(Elem_Idx,Line_Idx)

    !-----------------------------------------------------------------------
    ! include code to setup local profiles correctly 
    !-----------------------------------------------------------------------
    
    !Call Services module
    call FETCH_PIXEL_RTM_NWP(Input, Symbol, &
                             Elem_Idx,Line_Idx, RTM_NWP)
    
    Sfc_Level_RTM = RTM_NWP%Sfc_Level
    Tropo_Level_RTM = RTM_NWP%Tropo_Level
    
    Press_Prof_RTM =  RTM_NWP%P_Prof

    !do smoothing routines here - WCS3
    if (RTM_NWP%Smooth_Nwp_Fields_Flag_Temp == symbol%YES) then

       !--- height profile       
       Hght_Prof_RTM = INTERPOLATE_PROFILE_LOCAL( RTM_NWP%Z_Prof, &
                                            RTM_NWP%Z_Prof_1, &
                                            RTM_NWP%Z_Prof_2, &
                                            RTM_NWP%Z_Prof_3, &
                                            Inwp_Weight,Jnwp_Weight)

      !--- temperature profile
      Temp_Prof_RTM = INTERPOLATE_PROFILE_LOCAL( RTM_NWP%T_Prof, &
                                           RTM_NWP%T_Prof_1, &
                                           RTM_NWP%T_Prof_2, &
                                           RTM_NWP%T_Prof_3, &
                                           Inwp_Weight,Jnwp_Weight)

    
    else

       Hght_Prof_RTM = RTM_NWP%Z_Prof
       Temp_Prof_RTM = RTM_NWP%T_Prof
    
    endif

    !-----------------------------------------------------------------------------
    !--- Cloud Base and Top
    !---
    !--- Note 1. Extinction values are in km^(-1)
    !--- Note 2. All heights and thickness are converted to meters
    !-----------------------------------------------------------------------------
    if (Input%Zc(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL .and. &
        Input%Tau(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL) then

       Cloud_Extinction = WATER_EXTINCTION

       if (Cloud_Type == symbol%OPAQUE_ICE_TYPE .or. &
           Cloud_Type == symbol%OVERSHOOTING_TYPE) then
           Itemp = int(Input%Tc(Elem_Idx,Line_Idx))
           select case (Itemp)
            case (:199)    ; Cloud_Extinction = ICE_EXTINCTION1
            case (200:219) ; Cloud_Extinction = ICE_EXTINCTION2
            case (220:239) ; Cloud_Extinction = ICE_EXTINCTION3
            case (240:259) ; Cloud_Extinction = ICE_EXTINCTION4
            case (260:)    ; Cloud_Extinction = ICE_EXTINCTION5
           end select
       endif

       if (Cloud_Type == symbol%CIRRUS_TYPE .or. &
           Cloud_Type == symbol%OVERLAP_TYPE) then
           Itemp = int(Input%Tc(Elem_Idx,Line_Idx))
           select case (Itemp)
            case (:199)    ; Cloud_Extinction = CIRRUS_EXTINCTION1
            case (200:219) ; Cloud_Extinction = CIRRUS_EXTINCTION2
            case (220:239) ; Cloud_Extinction = CIRRUS_EXTINCTION3
            case (240:259) ; Cloud_Extinction = CIRRUS_EXTINCTION4
             case (260:)    ; Cloud_Extinction = CIRRUS_EXTINCTION5
          end select
       endif

       Cloud_Geometrical_Thickness = Input%Tau(Elem_Idx,Line_Idx) / Cloud_Extinction   !(km)
       Cloud_Geometrical_Thickness = Cloud_Geometrical_Thickness * 1000.0 !(m)

       if (Input%Tau(Elem_Idx,Line_Idx) < 2.0) then 
          Cloud_Geometrical_Thickness_Top_Offset = Cloud_Geometrical_Thickness/2.0    !(m)
       else
          Cloud_Geometrical_Thickness_Top_Offset = 1000.0 / Cloud_Extinction !(m)
       endif

       Zc_Top_Max = Hght_Prof_RTM(Tropo_Level_RTM)
       Zc_Base_Min = Hght_Prof_RTM(Sfc_Level_RTM)

       !-- new code
       if (Input%Zc(Elem_Idx,Line_Idx) > Zc_Top_Max .and. Cloud_Type == symbol%OVERSHOOTING_TYPE) then
          Output%Zc_Top(Elem_Idx,Line_Idx) = Input%Zc(Elem_Idx,Line_Idx)

          ! compute Pc_Top from Zc_Top
          call KNOWING_Z_COMPUTE_T_P(Output%Pc_Top(Elem_Idx,Line_Idx),R4_Dummy,Output%Zc_Top(Elem_Idx,Line_Idx),Ilev)

       else
          Output%Zc_Top(Elem_Idx,Line_Idx) = min(Zc_Top_Max, &
                                             Input%Zc(Elem_Idx,Line_Idx) +  &
                                             Cloud_Geometrical_Thickness_Top_Offset)

          ! compute Pc_Top from Zc_Top
          call KNOWING_Z_COMPUTE_T_P(Output%Pc_Top(Elem_Idx,Line_Idx),R4_Dummy,Output%Zc_Top(Elem_Idx,Line_Idx),Ilev)

       endif


!  Compute cloud base 
     if (Cloud_Type == symbol%CIRRUS_TYPE .and. Input%Tau(Elem_Idx,Line_Idx) < 1.0 .and. Input%Zc(Elem_Idx,Line_Idx)  > 10000.0) then
       Output%Zc_Base(Elem_Idx,Line_Idx) = min(Input%Zc(Elem_Idx,Line_Idx), &
                                           max(Zc_Base_Min,  &
                                           Output%Zc_Top(Elem_Idx,Line_Idx) - Cloud_Geometrical_Thickness))
     else
       if (Input%Zc(Elem_Idx,Line_Idx) > 0 .and. (Input%CWP(Elem_Idx,Line_Idx)  > 0 .or. Input%CWP_nwp(Elem_Idx,Line_Idx) > 0)) then
         call CIRA_base_hgt(Input%Zc(Elem_Idx,Line_Idx),Input%CWP(Elem_Idx,Line_Idx), Input%CWP_NWP(Elem_Idx,Line_Idx) ,&
              Cloud_Type,Input%CCL(Elem_Idx,Line_Idx),Input%Surface_Elevation(Elem_Idx,Line_Idx), &
              Cloud_Geometrical_Thickness_eff,Output%Zc_Base(Elem_Idx,Line_Idx))
       endif
     endif 

       ! compute Pc_Base from Zc_Base
       call KNOWING_Z_COMPUTE_T_P(Output%Pc_Base(Elem_Idx,Line_Idx),R4_Dummy,Output%Zc_Base(Elem_Idx,Line_Idx),Ilev)

 endif


 !---- null profile pointers each time  - REALLY?
 CALL NULL_PIX_POINTERS(Input, RTM_NWP)

 end do Element_Loop

end do Line_Loop

end subroutine CLOUD_BASE_ALGORITHM

!----------------------------------------------------------------------------
! Function INTERPOLATE_PROFILE_LOCAL
!
! general interpoLation routine for profiles
!
! input:
! lonx - longitude weighting factor
! Latx = Latitude weighting factor
! z1 = data(ilon, iLat)
! z2 = data(ilonx,iLat)
! z3 = data(ilon,iLatx)
! z4 = data(ilonx,iLatx)
!
! output:
! z = interpoLated profile
!
!
!---------------------------------------------------------------------------
 function INTERPOLATE_PROFILE_LOCAL(z1,z2,z3,z4,lonx,Latx) result(z)

  real, dimension(:), intent(in):: z1
  real, dimension(:), intent(in):: z2
  real, dimension(:), intent(in):: z3
  real, dimension(:), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: Latx
  real, dimension(size(z1)):: z

  !--- linear inteprpoLation scheme
  z =  (1.0-lonx) * ((1.0-Latx) * z1 + (Latx)* z3) + &
           (lonx) * ((1.0-Latx) * z2 + (Latx)* z4)

 end function INTERPOLATE_PROFILE_LOCAL


!------------------------------------------------------------------------------
! Null Pixel Level Pointers 
!------------------------------------------------------------------------------
subroutine NULL_PIX_POINTERS(Input, RTM_NWP)

   type(acha_input_struct), intent(inout) :: Input
   type(acha_rtm_nwp_struct), intent(inout) :: RTM_NWP

   RTM_NWP%T_Prof => null()
   RTM_NWP%T_Prof_1 => null() 
   RTM_NWP%T_Prof_2 => null() 
   RTM_NWP%T_Prof_3 => null()

   RTM_NWP%Z_Prof => null() 
   RTM_NWP%Z_Prof_1 => null() 
   RTM_NWP%Z_Prof_2 => null() 
   RTM_NWP%Z_Prof_3 => null() 

   if (Input%Chan_On_67um == sym%YES) then
     RTM_NWP%Atm_Rad_Prof_67um =>  null()
     RTM_NWP%Atm_Trans_Prof_67um =>  null()
     RTM_NWP%Black_Body_Rad_Prof_67um => null()
   endif
   if (Input%Chan_On_85um == sym%YES) then
     RTM_NWP%Atm_Rad_Prof_85um =>  null()
     RTM_NWP%Atm_Trans_Prof_85um =>  null()
   endif
   if (Input%Chan_On_11um == sym%YES) then
     RTM_NWP%Atm_Rad_Prof_11um => null()
     RTM_NWP%Atm_Trans_Prof_11um => null()
     RTM_NWP%Black_Body_Rad_Prof_11um => null()
   endif
   if (Input%Chan_On_12um == sym%YES) then
     RTM_NWP%Atm_Rad_Prof_12um => null()
     RTM_NWP%Atm_Trans_Prof_12um => null()
   endif
   if (Input%Chan_On_133um == sym%YES) then
     RTM_NWP%Atm_Rad_Prof_133um => null()
     RTM_NWP%Atm_Trans_Prof_133um => null()
   endif
 
end subroutine NULL_PIX_POINTERS

!-----------------------------------------------------------------
! InterpoLate within profiles knowing Z to determine T and P
!-----------------------------------------------------------------
   subroutine KNOWING_Z_COMPUTE_T_P(P,T,Z,Ilev)

     real, intent(in):: Z
     real, intent(out):: T
     real, intent(out):: P
     integer, intent(out):: Ilev
     real:: dp
     real:: dt
     real:: dz

     !--- interpoLate pressure profile
     call LOCATE(Hght_Prof_RTM,Num_Levels_RTM_Prof,Z,Ilev)
     Ilev = max(1,min(Num_Levels_RTM_Prof-1,Ilev))

     dp = Press_Prof_RTM(Ilev+1) - Press_Prof_RTM(Ilev)
     dt = Temp_Prof_RTM(Ilev+1) - Temp_Prof_RTM(Ilev)
     dz = Hght_Prof_RTM(Ilev+1) - Hght_Prof_RTM(Ilev)

     !--- perform interpoLation
     if (dz /= 0.0) then
           T = Temp_Prof_RTM(Ilev) + dt/dz * (Z - Hght_Prof_RTM(Ilev))
           P = Press_Prof_RTM(Ilev) + dp/dz * (Z - Hght_Prof_RTM(Ilev))
     else
           T = Temp_Prof_RTM(Ilev)
           P = Press_Prof_RTM(Ilev)
     endif

   end subroutine KNOWING_Z_COMPUTE_T_P

!-----------------------------------------------------------------
! CIRA's base code, interpret from IDL codes 
!-----------------------------------------------------------------
subroutine CIRA_base_hgt(Zc,Cwp,Cwp_nwp,Cloud_Type,CCL,Surf_Elev,Cloud_Geometrical_Thickness,Zc_base)

  real(kind=real4), intent(in) :: Zc,Cwp,Cwp_nwp,CCL,Surf_Elev
  integer,intent(in) :: Cloud_Type
  real(kind=real4), intent(out) :: Cloud_Geometrical_Thickness,Zc_base

! local variables
  real(kind=real4) :: Zc_local, Cwp_local,Cwp_nwp_local
  integer :: ibin
  integer :: ibin_max
  integer :: icwp
  real :: zdelta
  real :: slope
  real :: yint
  integer, parameter :: nbin = 9, npara = 6, ncwp = 2
  real,dimension(nbin,npara,ncwp) :: regr_coeff
  real,parameter :: mincbh = 0.0, maxcbh = 20.0*1000
  integer, dimension(5) :: qf
  integer :: cbh_qf
  integer :: statusflag,status_below_sfc,status_high_top

!                      min cth       ;slope         y-int          r2      n    median CWP
  regr_coeff(1,:,1) = [0.00000,      2.25812,     0.405590,    0.0236532, 5921., 0.0710000]
  regr_coeff(1,:,2) = [0.00000,     0.997031,     0.516989,    0.0900793, 5881., 0.0710000]
  regr_coeff(2,:,1) = [2.00000,      6.10980,     0.664818,    0.0664282, 3624., 0.114000]
  regr_coeff(2,:,2) = [2.00000,     0.913021,      1.35698,    0.0735795, 3621., 0.114000]
  regr_coeff(3,:,1) = [4.00000,      11.5574,      1.22527,    0.0519277, 2340., 0.110000]
  regr_coeff(3,:,2) = [4.00000,      1.37922,      2.58661,    0.0695758, 2329., 0.110000]
  regr_coeff(4,:,1) = [6.00000,      14.5382,      1.70570,    0.0568334, 2535., 0.123000]
  regr_coeff(4,:,2) = [6.00000,      1.68711,      3.62280,    0.0501604, 2511., 0.123000]
  regr_coeff(5,:,1) = [8.00000,      9.09855,      2.14247,    0.0218789, 3588., 0.131000]
  regr_coeff(5,:,2) = [8.00000,      2.45953,      3.86957,    0.0727178, 3579., 0.131000]
  regr_coeff(6,:,1) = [10.0000,      13.5772,      1.86554,    0.0497041, 4249., 0.127000]
  regr_coeff(6,:,2) = [10.0000,      4.83087,      3.53141,     0.160008, 4218., 0.127000]
  regr_coeff(7,:,1) = [12.0000,      16.0793,      1.64965,    0.0695903, 3154., 0.115000]
  regr_coeff(7,:,2) = [12.0000,      5.05173,      3.98610,     0.180965, 3121., 0.115000]
  regr_coeff(8,:,1) = [14.0000,      14.6030,      2.00010,    0.0429476, 2744., 0.116000]
  regr_coeff(8,:,2) = [14.0000,      6.06439,      4.03301,     0.239837, 2717., 0.116000]
  regr_coeff(9,:,1) = [16.0000,      9.26580,      2.29640,    0.0113376, 1455., 0.0990000]
  regr_coeff(9,:,2) = [16.0000,      6.60431,      3.26442,     0.227116, 1449., 0.0990000]

  qf = [0,1,2,3,4]
  cbh_qf = qf(1)
  statusflag = 0
  status_below_sfc = -1
  status_high_top = -2

! start retrieval
       Zc_local = Zc/1000.
       CWP_local = CWP/1000.
       Cwp_nwp_local = Cwp_nwp/1000.

! force large cwp to cap at 1.2 kg/m2
       if (CWP_local > 1.2)    CWP_local = 1.2
       if (Zc_local > 20.0)    Zc_local = 20.0
       if (Cwp_nwp_local > 1.2) Cwp_nwp_local = 1.2

       if (Cwp_local < 0 .and. Cwp_nwp_local > 0) Cwp_local = Cwp_nwp_local


       zdelta = 2.0
       ibin = floor(Zc_local)/floor(zdelta)+1
       ibin_max = 9

       if (Zc_local > 18.0 .or. ibin > ibin_max)    ibin = ibin_max

       icwp = 1
       if (CWP_local > regr_coeff(ibin,6,1))  icwp = 2

       slope = regr_coeff(ibin,2,icwp)
       yint = regr_coeff(ibin,3,icwp)
       Cloud_Geometrical_Thickness = slope*CWP_local+yint
       Cloud_Geometrical_Thickness = Cloud_Geometrical_Thickness*1000.

       if (Cloud_Geometrical_Thickness > Zc_local*1000) statusflag = status_below_sfc

       Zc_Base = Zc_local*1000-Cloud_Geometrical_Thickness

! adjustment based on cloud type
! type 9 is overshooting
       if (Cloud_Type == 9 .and. CCL > 0.0 .and. CCL < Zc_Base) then
           if (Zc_Base < Zc_local*1000 .and. Zc_Base > 2000. .and. Cloud_Geometrical_Thickness > 9000.) then
               Zc_Base = CCL
           endif

           if (Zc_local < 10000./1000 .and. Cloud_Geometrical_Thickness > 8000.) then
               Zc_Base = CCL
           endif
       endif

! apply quality flag
       if (Zc_Base < Surf_Elev) then
           Zc_Base = Surf_Elev
           cbh_qf = qf(3)
       endif

       if (Zc_Base < mincbh .or. Zc_Base > maxcbh) then
           Zc_Base = MISSING_VALUE_REAL
           cbh_qf = qf(4)
       endif

       if (Zc_Base >= Zc_local*1000) then
           Zc_Base = MISSING_VALUE_REAL
           cbh_qf = qf(5)
       endif

end subroutine CIRA_base_hgt

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module CLOUD_BASE
