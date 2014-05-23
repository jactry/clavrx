! $Id$
module AWG_CLOUD_HEIGHT
!---------------------------------------------------------------------
! This module houses the routines associated with...
!
! ACHA - AWG Cloud Height Algorithm
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! Reference:
!
!  Heidinger, A.K., and M.J. Pavolonis, 2009: Gazing at Cirrus Clouds for 25 Years
!  through a Split Window. Part I: Methodology. J. Appl. Meteor. Climatol., 48,
!  1100-1116.
!
!  Heidinger, A. K.; Pavolonis, M. J.; Holz, R. E.; Baum, Bryan A. and Berthier,
!  S.. Using CALIPSO to explore the sensitivity to cirrus height in the infrared
!  observations from NPOESS/VIIRS and GOES-R/ABI. Journal of Geophysical
!  Research, Volume 115, 2010, Doi:10.1029/2009JD012152. 
!
! Meta Data Flags
! 1 - Cloud Height Attempted (0 = no / 1 = yes)
! 2 - Bias Correction Employed (0 = no / 1 = yes)
! 3 - Ice Cloud Retrieval (0 = no / 1 = yes)
! 4 - Local Radiatve Center Processing Used (0 = no / 1 = yes)
! 5 - Multi-layer Retrieval (0 = no / 1 = yes)
! 6 - Lower Cloud Interpolation Used (0 = no / 1 = yes)
! 7 - Boundary Layer Inversion Assumed  (0 = no / 1 = yes)
! 8 - NWP Profile Inversion Assumed (0 = no / 1 = yes)
!
! Packed Quality Flags
! 1 - Processed (0 = no / 1 = yes)
! 2 - Valid Tc Retrieval (0 = yes, 1 = no)
! 3 - Valid ec Retrieval (0 = yes, 1 = no)
! 4 - Valid beta Retrieval (0 = yes, 1 = no)
! 5 - degraded Tc Retrieval (0 = no, 1 = yes)
! 6 - degraded ec Retrieval (0 = no, 1 = yes)
! 7 - degraded beta Retrieval (0 = no, 1 = yes)
!
!----------------------------------------------------------------------
  use ACHA_SERVICES_MOD !acha_services_mod.f90 in akh_clavrx_src

  implicit none

  public::  AWG_CLOUD_HEIGHT_ALGORITHM
  public::  SET_ACHA_MODE
  public::  CHECK_ACHA_MODE
  public:: SET_ACHA_VERSION
  public:: LOCAL_LINEAR_RADIATIVE_CENTER

  private:: SPATIALLY_INTERPOLATE_LOWER_CLOUD_POSITION  
  private:: KNOWING_P_COMPUTE_T_Z
  private:: KNOWING_T_COMPUTE_P_Z
  private:: KNOWING_Z_COMPUTE_T_P
  private:: GENERIC_PROFILE_INTERPOLATION
  private:: OPTIMAL_ESTIMATION_ACHA
  private:: COMPUTE_FORWARD_MODEL_AND_KERNEL
  private:: COMPUTE_APRIORI_BASED_ON_TYPE
  private:: COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO
  private:: DETERMINE_SFC_TYPE_FORWARD_MODEL
  private:: SET_CLEAR_SKY_COVARIANCE_TERMS
  private:: COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE

  private:: DETERMINE_ACHA_MODE_BASED_ON_CHANNELS
  private:: INTERPOLATE_PROFILE_ACHA
  private:: INTERPOLATE_NWP_ACHA
  private:: DETERMINE_INVERSION_LEVEL
  private:: DETERMINE_OPAQUE_LEVEL
  private:: COMPUTE_REFERENCE_LEVEL_EMISSIVITY
  private:: COMPUTE_STANDARD_DEVIATION
  private:: COMPUTE_TAU_REFF_ACHA 
  private:: NULL_PIX_POINTERS 

  !--- include the non-system specific variables
  include 'awg_cld_hght_include_1.inc'

  !--- interpolated profiles
  real, private, dimension(Num_Levels_Rtm_Prof) :: Temp_Prof_Rtm
  real, private, dimension(Num_Levels_Rtm_Prof) :: Press_Prof_Rtm
  real, private, dimension(Num_Levels_Rtm_Prof) :: Hght_Prof_Rtm
  real, private, dimension(:), pointer :: Atm_Rad_Prof_67um_Rtm
  real, private, dimension(:), pointer :: Atm_Rad_Prof_85um_Rtm
  real, private, dimension(:), pointer :: Atm_Rad_Prof_11um_Rtm
  real, private, dimension(:), pointer :: Atm_Rad_Prof_12um_Rtm
  real, private, dimension(:), pointer :: Atm_Rad_Prof_133um_Rtm
  real, private, dimension(:), pointer :: Atm_Trans_Prof_67um_Rtm
  real, private, dimension(:), pointer :: Atm_Trans_Prof_85um_Rtm
  real, private, dimension(:), pointer :: Atm_Trans_Prof_11um_Rtm
  real, private, dimension(:), pointer :: Atm_Trans_Prof_12um_Rtm
  real, private, dimension(:), pointer :: Atm_Trans_Prof_133um_Rtm
  real, private, dimension(:), pointer :: Black_Body_Rad_Prof_11um_Rtm
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_67um_Rtm_Unbiased
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_85um_Rtm_Unbiased
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_11um_Rtm_Unbiased
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_12um_Rtm_Unbiased
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_133um_Rtm_Unbiased
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_67um_Rtm_Unbiased
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_85um_Rtm_Unbiased
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_11um_Rtm_Unbiased
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_12um_Rtm_Unbiased
  real, private, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_133um_Rtm_Unbiased
  integer, private:: Inver_Level_Rtm
  integer, private:: Sfc_Level_Rtm
  integer, private:: Tropo_Level_Rtm

  real, private:: Bt_67um_Mean
  real, private:: Bt_85um_Mean
  real, private:: Bt_11um_Mean
  real, private:: Bt_12um_Mean
  real, private:: Bt_133um_Mean

  real, private:: Bt_67um_Bt_67um_Covar
  real, private:: Bt_85um_Bt_85um_Covar
  real, private:: Bt_11um_Bt_11um_Covar
  real, private:: Bt_12um_Bt_12um_Covar
  real, private:: Bt_133um_Bt_133um_Covar

  real, private:: Bt_85um_Bt_133um_Covar
  real, private:: Bt_12um_Bt_133um_Covar
  real, private:: Bt_11um_Bt_67um_Covar
  real, private:: Bt_11um_Bt_85um_Covar
  real, private:: Bt_12um_Bt_85um_Covar
  real, private:: Bt_12um_Bt_67um_Covar
  real, private:: Bt_11um_Bt_12um_Covar
  real, private:: Bt_11um_Bt_133um_Covar
  real, private:: Bt_67um_Bt_133um_Covar
 
  real, private:: Bt_11um_Btd_11um_67um_Covar
  real, private:: Bt_11um_Btd_11um_85um_Covar
  real, private:: Bt_11um_Btd_11um_12um_Covar
  real, private:: Bt_11um_Btd_11um_133um_Covar

  real, private:: Btd_11um_67um_Btd_11um_67um_Covar
  real, private:: Btd_11um_85um_Btd_11um_85um_Covar
  real, private:: Btd_11um_12um_Btd_11um_12um_Covar
  real, private:: Btd_11um_133um_Btd_11um_133um_Covar
  real, private:: Btd_11um_12um_Btd_11um_133um_Covar
  real, private:: Btd_11um_12um_Btd_11um_67um_Covar
  real, private:: Btd_11um_12um_Btd_11um_85um_Covar
  real, private:: Btd_11um_67um_Btd_11um_133um_Covar

  real, private, PARAMETER:: Dt_Dz_Strato = -6500.0 !K/m
  real, private, PARAMETER:: Sensor_Zenith_Threshold = 70.0

  contains 

!------------------------------------------------------------------------------
! AWG Cloud Height Algorithm (ACHA)
!
! Author: Andrew Heidinger, NOAA
!
! Assumptions
!   1) No scattering
!   2) single layer cloud for cloud type /= 6
!   3) for overlap type, an opaque cloud 200 mb above the surface lies below 
!
! Limitations
!   1) sensitivity to Tc is low for thin clouds
!   2) little Emissivity sensitivity for low clouds
!
! input to the retrieval
!  y(1) = t4 - 11 micron brightness temperature
!  y(2) = t4 - t5 - the split window temperature
!
! the output of the retrieval
!  x(1) - the cloud temperature
!  x(2) - the 11 micron Emissivity at nadir
!  x(3) - the beta ratio for 11 and 12 microns
!
! This routine uses a 1d-var retrieval approach as outlined in Rodger (1976)
!
! input to the 1d-var approach
!  y - the vector of observations
!  x_Ap - the a apriori estimates of x
!  Sa - the error covariance matric of x_Ap
!  Sy - the error covariance of y (included calibration, forward model)
!
!  the Optimal Estimation Quality Flags are determined as follows
!  3 - estimated error < 1/3 a priori error
!  2 - estimated error < 2/3 a priori error
!  1 - any other converged retrieval
!  0 - a failed or unattempted retrieval
!
! the overall quality flag Description
! 0 - No retrieval attempted
! 1 - Retrieval attempted and failed
! 2 - Marginally Successful Retrieval
! 3 - Fully Successful Retrieval
!
!
! Meta Data
! 1 - Cloud Height Attempted (0 = no / 1 = yes)
! 2 - Bias Correction Employed (0 = no / 1 = yes)
! 3 - Ice Cloud Retrieval (0 = no / 1 = yes)
! 4 - Local Radiatve Center Processing Used (0 = no / 1 = yes)
! 5 - Multi-layer Retrieval (0 = no / 1 = yes)
! 6 - Lower Cloud Interpolation Used (0 = no / 1 = ! yes)
! 7 - Boundary Layer Inversion Assumed  (0 = ! no / 1 = yes)
!
! Processing Order Description
! pass 0 = treat all pixels as single layer ncc pixels (only done if Use_Lrc_Flag=no)
! pass 1 = non-multi-layer ncc pixels
! pass 2 = single layer water cloud pixels
! pass 3 = lrc multi-layer clouds
! pass 4 = all remaining clouds
!
!
!----------------------------------------------------------------------
! modification history
!
! July 2006 - Added beta as an element of x
! October 2006 - Added cloud lapse rate to make Tc more related to true 
!                cloud-top temperature
!
!
!------------------------------------------------------------------------------
  subroutine  AWG_CLOUD_HEIGHT_ALGORITHM(Input, symbol, Output)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(symbol_acha), intent(inout) :: symbol
  type(acha_input_struct), intent(inout) :: Input
  type(acha_output_struct), intent(inout) :: Output

  !===============================================================================
  !  Pixel level RTM structure
  !===============================================================================
 
  type(acha_rtm_nwp_struct) :: Acha_RTM_NWP

  !===============================================================================
  !  Local Variable Declaration
  !===============================================================================

  integer (kind=int4):: Smooth_Nwp_Fields_Flag_Temp
  integer :: ACHA_Mode_Flag
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Param_Idx
  integer:: i
  integer:: Fail_Flag
  integer:: Singular_Flag
  integer:: Ilev
  integer:: Inwp
  integer:: Jnwp
  integer:: Inwp_x
  integer:: Jnwp_x
  real:: Inwp_Weight
  real:: Jnwp_Weight
  integer:: Ivza
  integer:: ilrc
  integer:: jlrc
  integer:: Iter_Idx
  integer:: idiag_output
  integer:: ierror
  integer:: ipass
  integer:: ipass_min
  integer:: ipass_Max
  real:: Lapse_Rate
  real:: Lapse_Rate_dP_dZ
  real:: Delta_Cld_Temp_Sfc_Temp

  real:: conv_crit

  !--- ch27 variables
  real:: Rad_Ac_67um
  real:: Trans_Ac_67um
  real:: Rad_Clear_67um
  real:: Bc_67um
  real::  Bt_Clear_67um

  !--- ch29 variables
  real:: Rad_Ac_85um
  real:: Trans_Ac_85um
  real:: Rad_Clear_85um
  real:: Bc_85um
  real::  Bt_Clear_85um

  !--- ch31 variables
  real:: Rad_Ac_11um
  real:: Trans_Ac_11um
  real:: Rad_Clear_11um
  real:: Bc_11um
  real::  Bt_Clear_11um
  real::  Bt_11um_Lrc

  !--- ch32 variables
  real:: Rad_Ac_12um
  real:: Trans_Ac_12um
  real:: Rad_Clear_12um
  real:: Bc_12um
  real::  Bt_Clear_12um

  !--- ch33 variables
  real:: Rad_Ac_133um
  real:: Trans_Ac_133um
  real:: Rad_Clear_133um
  real:: Bc_133um
  real::  Bt_Clear_133um

  real:: a_Beta_11um_133um_fit
  real:: b_Beta_11um_133um_fit
  real:: a_Beta_11um_85um_fit
  real:: b_Beta_11um_85um_fit
  real:: a_Beta_11um_67um_fit
  real:: b_Beta_11um_67um_fit

  real:: Tc_Opaque_Level
  real:: Pc_Opaque_Level
  real:: Zc_Opaque_Level
  real:: Tsfc_Est
  real:: Tc_temp
  real:: Pc_temp
  real:: Zc_temp
  real:: Tc_Ap
  real:: Ec_Ap
  real:: Beta_Ap
  real:: Tc_Ap_Uncer
  real:: Ec_Ap_Uncer
  real:: Beta_Ap_Uncer
  real:: T11um_Clr_Uncer
  real:: T11um_12um_Clr_Uncer
  real:: T11um_133um_Clr_Uncer
  real:: T11um_85um_Clr_Uncer
  real:: T11um_67um_Clr_Uncer
  real:: R4_Dummy
  real:: Emiss_11um_Tropo
  integer:: Num_Obs
  integer:: Cloud_Type
  integer:: Cloud_Phase
  integer:: Undetected_Cloud
  integer:: Converged_Flag
  integer:: Sfc_Type_Forward_Model
  integer(kind=int1), dimension(NUM_META_DATA):: Meta_Data_Flags

  !--- 1d-var retrieval arrays
  real (kind=real4), allocatable, dimension(:):: y
  real (kind=real4), allocatable, dimension(:):: f
  real (kind=real4), allocatable, dimension(:):: x
  real (kind=real4), allocatable, dimension(:):: x_Ap
  real (kind=real4), allocatable, dimension(:):: Delta_X
  real (kind=real4), allocatable, dimension(:,:):: K
  real (kind=real4), allocatable, dimension(:,:):: Sa
  real (kind=real4), allocatable, dimension(:,:):: Sa_inv
  real (kind=real4), allocatable, dimension(:,:):: Sx
  real (kind=real4), allocatable, dimension(:,:):: Sx_inv
  real (kind=real4), allocatable, dimension(:,:):: E
  real (kind=real4), allocatable, dimension(:,:):: Sy
  real (kind=real4), allocatable, dimension(:,:):: Sy_inv
  real (kind=real4), allocatable, dimension(:):: y_variance
  real (kind=real4), allocatable, dimension(:):: Emiss_Vector
  real (kind=real4), allocatable, dimension(:,:):: AKM  !Averaging Kernel Matrix

  !--- local POINTERs to global arrays or data structures
  integer(kind=int4), allocatable, dimension(:,:):: Elem_Idx_LRC
  integer(kind=int4), allocatable, dimension(:,:):: Line_Idx_LRC
  integer(kind=int1), allocatable, dimension(:,:):: Skip_LRC_Mask
  real (kind=real4):: Bt_11um_Std
  real (kind=real4):: Btd_11um_67um_Std
  real (kind=real4):: Btd_11um_85um_Std
  real (kind=real4):: Btd_11um_12um_Std
  real (kind=real4):: Btd_11um_133um_Std

  !--- scalar local variables
  integer (kind=int4):: i1,i2,j1,j2
  real (kind=real4):: Cloud_Extinction
  real (kind=real4):: Cloud_Geometrical_Thickness
  real (kind=real4):: Cloud_Geometrical_Thickness_Top_Offset
  real (kind=real4):: Zc_Top_Max
  real (kind=real4):: Zc_Base_Min
  integer (kind=int4):: NWP_Profile_Inversion_Flag
  integer (kind=int4):: Itemp

!-----------------------------------------------------------------------
! BEGIN EXECUTABLE CODE
!-----------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! abort if no 11 um channel
  !----------------------------------------------------------------------------
  if (Input%Chan_On_11um == symbol%NO) THEN
     Output%Packed_Qf =  0
     return
  endif

  !---------------------------------------------------------------------------
  !-- Acha Mode set to  -1, determine based on channels
  !---------------------------------------------------------------------------
  ACHA_Mode_Flag = Input%ACHA_Mode_Flag_In
  if (ACHA_Mode_Flag < 0) then
    call DETERMINE_ACHA_MODE_BASED_ON_CHANNELS(symbol, &
                                      Acha_Mode_Flag, &
                                      Input%Chan_On_67um, &
                                      Input%Chan_On_85um,  &
                                      Input%Chan_On_11um,  &
                                      Input%Chan_On_12um,  &
                                      Input%Chan_On_133um)
  endif

  !--- determine number of channels
  select case(Acha_Mode_Flag)
     case(0)  !11 avhrr/1
       Num_Obs = 1
     case(1)  !11,12 avhrr/2/3
       Num_Obs = 2
     case(2)  !11,13.3 goes-nop
       Num_Obs = 2
     case(3)  !11,12,13.3 goes-r
       Num_Obs = 3
     case(4)  !11,12,8.5 viirs 
       Num_Obs = 3
     case(5)  !11,12,6.7
       Num_Obs = 3
     case(6)  !11,13.3,6.7
       Num_Obs = 3
     case(7)  !11,6.7
       Num_Obs = 2
  end select

  !--- allocate needed 2d arrays for processing this segment
  allocate(Elem_Idx_LRC(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Line_Idx_LRC(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Skip_LRC_Mask(Input%Number_of_Elements,Input%Number_of_Lines))

  !--- allocate 1D-VAR arrays based on number of channels
  allocate(y(Num_Obs))
  allocate(y_variance(Num_Obs))
  allocate(f(Num_Obs))
  allocate(x(Num_Param))
  allocate(x_Ap(Num_Param))
  allocate(Delta_X(Num_Param))
  allocate(K(Num_Obs,Num_Param))
  allocate(Sa(Num_Param,Num_Param))
  allocate(Sa_inv(Num_Param,Num_Param))
  allocate(Sx(Num_Param,Num_Param))
  allocate(Sx_inv(Num_Param,Num_Param))
  allocate(E(Num_Param,Num_Param))
  allocate(Sy(Num_Obs,Num_Obs))
  allocate(Sy_inv(Num_Obs,Num_Obs))
  allocate(Emiss_Vector(Num_Obs))
  allocate(AKM(Num_Param,Num_Param))

  !--- turn on diagnostic output
  idiag_output = symbol%NO

  !--- set convergence criterion
  conv_crit = Num_Param - 1.0

  !----------- make identity matrix
  E = 0.0
  do i = 1,Num_Param
    E(i,i) = 1.0
  enddo

  !--- initialize output
  Output%Tc =  Missing_Value_Real4
  Output%Ec =  Missing_Value_Real4
  Output%Beta =  Missing_Value_Real4
  Output%Pc =  Missing_Value_Real4
  Output%Zc =  Missing_Value_Real4
  Output%OE_Qf = 0
  Output%Qf = 0
  Output%Cloud_Layer = 0
  Meta_Data_Flags = 0
  Bc_67um = 0
  Bc_85um = 0
  Bc_11um = 0
  Bc_12um = 0
  
  !--------------------------------------------------------------------------
  ! spatial processing pixels
  ! compute local radiative centers using 11 um brightness temperature
  !---------------------------------------------------------------------------

  !--- construct a mask to select pixel for LRC computation
  Elem_Idx_LRC = Missing_Value_Int4
  Line_Idx_LRC = Missing_Value_Int4
  Skip_LRC_Mask = Input%Invalid_Data_Mask

  !--- call LRC routine
  if (Use_Lrc_Flag == symbol%YES) then

    if (associated(Input%Elem_Idx_LRC_Input) .and. &
        associated(Input%Line_Idx_LRC_Input)) then

      Elem_Idx_LRC = Input%Elem_Idx_LRC_Input
      Line_Idx_LRC = Input%Line_Idx_LRC_Input

    else

      where(Input%Cloud_Mask_Local == symbol%CLEAR .or. Input%Cloud_Mask_Local == symbol%PROB_CLEAR)
          Skip_LRC_Mask = symbol%YES
      endwhere

      call LOCAL_LINEAR_RADIATIVE_CENTER(symbol%YES,symbol%NO,&
                                         LRC_Meander_Flag, &
                                         Input%Bt_11um, &
                                         Element_Idx_Min, &
                                         Input%Number_of_Elements, & 
                                         Line_Idx_Min,  &
                                         Input%Number_of_Lines, & 
                                         Max_LRC_Distance,  &
                                         Min_LRC_Jump,  &
                                         Max_LRC_Jump,  &
                                         Grad_Flag_LRC,  &
                                         Missing_Value_Int4, &
                                         Skip_LRC_Mask, &
                                         Min_Bt_11um_Lrc,  &
                                         Max_Bt_11um_Lrc, &
                                         Elem_Idx_LRC,  &
                                         Line_Idx_LRC)

    endif
  endif

  !--------------------------------------------------------------------------
  ! determine processing order of pixels
  !--------------------------------------------------------------------------

  !-- intialize to missing is bad pixels and zero everywhere else
  Output%Processing_Order = Missing_Value_Int1
  where(Input%Invalid_Data_Mask == symbol%NO) 
    Output%Processing_Order =  0
  endwhere

  if (Use_Lrc_Flag == symbol%NO) then
      ipass_min = 0
      ipass_Max = 0
      Meta_Data_Flags(4) = symbol%NO
  else 
      ipass_min = 1
      ipass_Max = 4
      Meta_Data_Flags(4) = symbol%YES

      !--- loop through pixels, determine processing order
      Line_loop_0: do Line_Idx = Line_Idx_Min,Input%Number_of_Lines + Line_Idx_Min - 1
        Element_Loop_0:   do Elem_Idx = Element_Idx_Min, Input%Number_of_Elements + Element_Idx_Min - 1

        !--- skip data marked as bad
        if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) then
            cycle
        endif

        !--- skip data with viewing zenith angle greater than prescribed threshold
        if (Input%Sensor_Zenith_Angle(Elem_Idx,Line_Idx) > Sensor_Zenith_Threshold) then
            cycle
        endif

        ilrc = Elem_Idx_LRC(Elem_Idx,Line_Idx)
        jlrc = Line_Idx_LRC(Elem_Idx,Line_Idx)

        Cloud_Type = Input%Cloud_Type_Local(Elem_Idx,Line_Idx)

        !-- on pass 1, do single layer lrc's
        if ((Elem_Idx == ilrc) .and. (Line_Idx == jlrc) .and. &
            (Cloud_Type /= symbol%OVERLAP_TYPE)) then
             Output%Processing_Order(Elem_Idx,Line_Idx) = 1
             cycle
        endif

        !-- on pass 2, do non-lrc water clouds
        if (((Elem_Idx /= ilrc) .or. (Line_Idx /= jlrc)) .and. &
            (Cloud_Type == symbol%FOG_TYPE .or. &
            Cloud_Type == symbol%WATER_TYPE .or. &
            Cloud_Type == symbol%MIXED_TYPE .or. &
            Cloud_Type == symbol%SUPERCOOLED_TYPE)) then
            Output%Processing_Order(Elem_Idx,Line_Idx) = 2
            cycle
        endif

        !-- on pass 3, do lrc overlap clouds
        if ((Elem_Idx == ilrc) .and. (Line_Idx == jlrc) .and. &
            (Cloud_Type == symbol%OVERLAP_TYPE)) then
             Output%Processing_Order(Elem_Idx,Line_Idx) = 3
            cycle
        endif

        !--  on pass-4 do remaining
        if (Output%Processing_Order(Elem_Idx,Line_Idx) == 0) then
           Output%Processing_Order(Elem_Idx,Line_Idx) = 4
        endif

       end do Element_Loop_0
     end do Line_loop_0
   endif

  !--------------------------------------------------------------------------
  ! perform multiple passes if using lrc results
  !--------------------------------------------------------------------------

  pass_loop: do ipass = ipass_min, ipass_Max
  
   !--------------------------------------------------------------------------
   ! on the third pass, spatially interpolate water cloud temperature
   !--------------------------------------------------------------------------
   if ((ipass == 0) .or. (ipass == 3)) then
     call SPATIALLY_INTERPOLATE_LOWER_CLOUD_POSITION(symbol, &
                      ipass,Line_Idx_min,Input%Number_of_Lines, &
                      Input%Elem_Idx_Nwp,Input%Line_Idx_Nwp, &
                      Input%Invalid_Data_Mask, &
                      Input%Cloud_Type_Local, &
                      Input%Surface_Pressure, &
                      Output%Pc, &
                      Output%Lower_Cloud_Pressure, &
                      Output%Lower_Cloud_Temperature, &
                      Output%Lower_Cloud_Height)
   endif

   !--------------------------------------------------------------------------
   ! loop over pixels in scanlines
   !--------------------------------------------------------------------------
   Line_loop: do Line_Idx = Line_Idx_min,Input%Number_of_Lines + Line_Idx_min - 1

    Element_Loop:   do Elem_Idx = 1, Input%Number_of_Elements

    !---- null profile pointers each time - WCS3
    call NULL_PIX_POINTERS(Input, Acha_RTM_NWP)

    !--- check if pixel should be processd in this path
    if (ipass /= Output%Processing_Order(Elem_Idx,Line_Idx)) then
          cycle
    endif

    !--- check for a bad data or data outside of viewing angle range
    if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES .or. &
        Input%Sensor_Zenith_Angle(Elem_Idx,Line_Idx) > Sensor_Zenith_Threshold) then

          Output%Packed_Qf(Elem_Idx,Line_Idx) =  0
          Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0

    endif

    !--- for convenience, save nwp indices to local variables
    Inwp = Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx)
    Jnwp = Input%Line_Idx_Nwp(Elem_Idx,Line_Idx)
    Inwp_x = Input%Elem_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
    Jnwp_x = Input%Line_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
    Inwp_Weight = Input%Longitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
    Jnwp_Weight = Input%Latitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
    Ivza =  Input%Viewing_Zenith_Angle_Idx_Rtm(Elem_Idx,Line_Idx)
    ilrc = Elem_Idx_LRC(Elem_Idx,Line_Idx)
    jlrc = Line_Idx_LRC(Elem_Idx,Line_Idx)
    
    !--- Qc indices (check how Framework handles it)
    if (Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx) == Missing_Value_Int4 .or. &
        Input%Line_Idx_Nwp(Elem_Idx,Line_Idx) == Missing_Value_Int4 .or. &
        Input%Viewing_Zenith_Angle_Idx_Rtm(Elem_Idx,Line_Idx) == Missing_Value_Int4) then 
          Output%Packed_Qf(Elem_Idx,Line_Idx) =  0
          Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0
         cycle 
    endif

    !-----------------------------------------------------------------------
    ! include code to setup local profiles correctly 
    !-----------------------------------------------------------------------
    
    !Call framework services module
    call ACHA_FETCH_PIXEL_NWP_RTM(Input, symbol, &
                                  Elem_Idx,Line_Idx, Acha_RTM_NWP)
    
    
    Sfc_Level_RTM = Acha_RTM_NWP%Sfc_Level
    Tropo_Level_RTM = Acha_RTM_NWP%Tropo_Level
    
    Press_Prof_Rtm =  Acha_RTM_NWP%P_Prof

    !do smoothing routines here - WCS3
    if (Acha_RTM_NWP%Smooth_Nwp_Fields_Flag_Temp == symbol%YES) then

       !--- height profile       
       Hght_Prof_RTM = INTERPOLATE_PROFILE_ACHA( Acha_RTM_NWP%Z_Prof, &
                                                 Acha_RTM_NWP%Z_Prof_1, &
                                                 Acha_RTM_NWP%Z_Prof_2, &
                                                 Acha_RTM_NWP%Z_Prof_3, &
                                                 Inwp_Weight,Jnwp_Weight)

      !--- temperature profile
      Temp_Prof_RTM = INTERPOLATE_PROFILE_ACHA( Acha_RTM_NWP%T_Prof, &
                                                Acha_RTM_NWP%T_Prof_1, &
                                                Acha_RTM_NWP%T_Prof_2, &
                                                Acha_RTM_NWP%T_Prof_3, &
                                                Inwp_Weight,Jnwp_Weight)

    
    else

       Hght_Prof_RTM = Acha_RTM_NWP%Z_Prof
       Temp_Prof_RTM = Acha_RTM_NWP%T_Prof
    
    endif

   !--- populate radiance and transmission profiles
   if (Input%Chan_On_67um == sym%YES) then
     Atm_Rad_Prof_67um_RTM =>  Acha_RTM_NWP%Atm_Rad_Prof_67um_RTM
     Atm_Trans_Prof_67um_RTM =>  Acha_RTM_NWP%Atm_Trans_Prof_67um_RTM
   endif
   if (Input%Chan_On_85um == sym%YES) then
     Atm_Rad_Prof_85um_RTM =>  Acha_RTM_NWP%Atm_Rad_Prof_85um_RTM
     Atm_Trans_Prof_85um_RTM =>  Acha_RTM_NWP%Atm_Trans_Prof_85um_RTM
   endif
   if (Input%Chan_On_11um == sym%YES) then
      Atm_Rad_Prof_11um_RTM =>  Acha_RTM_NWP%Atm_Rad_Prof_11um_RTM
      Atm_Trans_Prof_11um_RTM =>  Acha_RTM_NWP%Atm_Trans_Prof_11um_RTM
      Black_Body_Rad_Prof_11um_RTM =>  Acha_RTM_NWP%Black_Body_Rad_Prof_11um_RTM
   endif
   if (Input%Chan_On_12um == sym%YES) then
      Atm_Rad_Prof_12um_RTM =>  Acha_RTM_NWP%Atm_Rad_Prof_12um_RTM
      Atm_Trans_Prof_12um_RTM =>  Acha_RTM_NWP%Atm_Trans_Prof_12um_RTM
   endif
   if (Input%Chan_On_133um == sym%YES) then
      Atm_Rad_Prof_133um_RTM =>  Acha_RTM_NWP%Atm_Rad_Prof_133um_RTM
      Atm_Trans_Prof_133um_RTM =>  Acha_RTM_NWP%Atm_Trans_Prof_133um_RTM
   endif

   !-----------------------------------------------------------------------
   !  find opaque levels
   !-----------------------------------------------------------------------
   CALL DETERMINE_OPAQUE_LEVEL(Input%Rad_11um(Elem_Idx,Line_Idx), &
                                Black_Body_Rad_Prof_11um_Rtm, &
                                Press_Prof_Rtm, &
                                Hght_Prof_Rtm, &
                                Temp_Prof_Rtm, &
                                Tropo_Level_Rtm, &
                                Sfc_Level_Rtm, &
                                Pc_Opaque_Level, &
                                Tc_Opaque_Level, &
                                Zc_Opaque_Level)

   !-------------------------------------------------------------------
   ! Apply Opaque Retrieval for Acha_Mode_Flag = 0, then cycle
   !-------------------------------------------------------------------
   if (Acha_Mode_Flag == 0) then
        if (((Input%Cloud_Mask_Local(Elem_Idx,Line_Idx) == symbol%CLEAR) .or.  &
            (Input%Cloud_Mask_Local(Elem_Idx,Line_Idx) == symbol%PROB_CLEAR)) .and. &
            (Input%Process_Undetected_Cloud_Flag_Local == symbol%NO)) then
          Output%Tc(Elem_Idx,Line_Idx) = Missing_Value_Real4
          Output%Pc(Elem_Idx,Line_Idx) = Missing_Value_Real4
          Output%Zc(Elem_Idx,Line_Idx) = Missing_Value_Real4
          Output%Ec(Elem_Idx,Line_Idx) = Missing_Value_Real4
          Output%Beta(Elem_Idx,Line_Idx) = Missing_Value_Real4
        else
          Output%Tc(Elem_Idx,Line_Idx) = Tc_Opaque_Level
          Output%Pc(Elem_Idx,Line_Idx) = Pc_Opaque_Level
          Output%Zc(Elem_Idx,Line_Idx) = Zc_Opaque_Level
          Output%Ec(Elem_Idx,Line_Idx) = 1.0
          Output%Beta(Elem_Idx,Line_Idx) = Missing_Value_Real4
        endif
        Output%Tc_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Output%Pc_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Output%Zc_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Output%Ec_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Output%Beta_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Output%Packed_Qf(Elem_Idx,Line_Idx) =  0
        Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0
        cycle
   endif


   !----------------------------------------------------------------------
   ! determine cloud phase from cloud type for convienience
   !----------------------------------------------------------------------
   Cloud_Type = Input%Cloud_Type_Local(Elem_Idx,Line_Idx)
   Cloud_Phase = symbol%UNKNOWN_PHASE

   if ( (Cloud_Type  == symbol%FOG_TYPE) .or. &
       (Cloud_Type  == symbol%WATER_TYPE) .or. &
       (Cloud_Type  == symbol%SUPERCOOLED_TYPE)) then
       Cloud_Phase = symbol%WATER_PHASE
   endif

   if ( (Cloud_Type  == symbol%CIRRUS_TYPE) .or. &
        (Cloud_Type  == symbol%OVERLAP_TYPE) .or. &
        (Cloud_Type  == symbol%TICE_TYPE)) then
        Cloud_Phase = symbol%ICE_PHASE
   endif


  !-----------------------------------------------------------------------
  !----- data quality check
  !-----------------------------------------------------------------------
  if ((Input%Bt_11um(Elem_Idx,Line_Idx) < 170.0) .or. &         !begin data check
      (Input%Bt_11um(Elem_Idx,Line_Idx) > 340.0) .or. &
      (Input%Surface_Temperature(Elem_Idx,Line_Idx) < 180.0) .or. &
      (Input%Surface_Temperature(Elem_Idx,Line_Idx) > 340.0) .or. &
      (Input%Tropopause_Temperature(Elem_Idx,Line_Idx) < 160.0) .or. &
      (Input%Tropopause_Temperature(Elem_Idx,Line_Idx) > 270.0)) then

     Output%Tc(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Output%Ec(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Output%Beta(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Output%Pc(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Output%Zc(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Output%Qf(Elem_Idx,Line_Idx) = 0
     Output%Cloud_Layer(Elem_Idx,Line_Idx) = 0

   else  !if passed data check then proceed with retrieval

     !---------------------------------------------------------------------
     ! select to do retrievals for all pixels or just cloudy ones
     !---------------------------------------------------------------------
     Undetected_Cloud = symbol%NO

     if ((Input%Cloud_Mask_Local(Elem_Idx,Line_Idx) == symbol%CLEAR) .or.  &
       (Input%Cloud_Mask_Local(Elem_Idx,Line_Idx) == symbol%PROB_CLEAR)) then

       if (Input%Process_Undetected_Cloud_Flag_Local == symbol%NO) then
               cycle
       else
              Undetected_Cloud = symbol%YES
       endif

     endif
  

   !----------------------------------------------------------------------
   !--- Set Meta Data Flags
   !----------------------------------------------------------------------
   Meta_Data_Flags(1) = symbol%YES
   if (Cloud_Phase == symbol%ICE_PHASE) then
       Meta_Data_Flags(3) = symbol%YES
   else
       Meta_Data_Flags(3) = symbol%NO
   endif
   if (Use_Lrc_Flag == symbol%YES) then
       Meta_Data_Flags(4) = symbol%YES
   else
       Meta_Data_Flags(4) = symbol%NO
   endif
   if (Cloud_Type == symbol%OVERLAP_TYPE) then
       Meta_Data_Flags(5) = symbol%YES
   else
       Meta_Data_Flags(5) = symbol%NO
   endif
   Meta_Data_Flags(6) = symbol%NO     !lower cloud interpolation
   Meta_Data_Flags(7) = symbol%NO     !low level inversion
   Meta_Data_Flags(8) = symbol%NO     !NWP profile inversion
   
   !-----------------------------------------------------------------------
   ! compute needed channel 3x3 standard deviations
   !-----------------------------------------------------------------------
   j1 = max(Line_Idx_Min, Line_Idx - 1)
   j2 = min(Input%Number_of_Lines, Line_Idx + 1)
   i1 = max(Element_Idx_Min, Elem_Idx - 1) 
   i2 = min(Input%Number_of_Elements, Elem_Idx + 1)

   Bt_11um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_11um(i1:i2,j1:j2),Input%Invalid_Data_Mask(i1:i2,j1:j2))

   if (Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 7) then
    Btd_11um_67um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_11um(i1:i2,j1:j2) -  Input%Bt_67um(i1:i2,j1:j2),&
                                                   Input%Invalid_Data_Mask(i1:i2,j1:j2))
   endif
   if (Acha_Mode_Flag == 4) then
    Btd_11um_85um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_11um(i1:i2,j1:j2) -  Input%Bt_85um(i1:i2,j1:j2), &
                                                  Input%Invalid_Data_Mask(i1:i2,j1:j2))
   endif

   if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag == 3 .or.  &
       Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 5) then
    Btd_11um_12um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_11um(i1:i2,j1:j2) -  Input%Bt_12um(i1:i2,j1:j2), &
                                                   Input%Invalid_Data_Mask(i1:i2,j1:j2))
   endif

   if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 6) then
    Btd_11um_133um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_11um(i1:i2,j1:j2) -  Input%Bt_133um(i1:i2,j1:j2), &
                                                    Input%Invalid_Data_Mask(i1:i2,j1:j2))
   endif
   
   !-----------------------------------------------------------------------
   ! assign values to y and x_Ap
   !----------------------------------------------------------------------

   !--- y - the observation vOutput%Ector
   select case(Acha_Mode_Flag)
     case(0)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
     case(1)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
     case(2)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_133um_Std**2 
     case(3)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
       y_variance(3) = Btd_11um_133um_Std**2 
     case(4)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_85um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
       y_variance(3) = Btd_11um_85um_Std**2 
     case(5)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_67um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(3) = Btd_11um_12um_Std**2 
       y_variance(3) = Btd_11um_67um_Std**2 
     case(6)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_67um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_133um_Std**2 
       y_variance(3) = Btd_11um_67um_Std**2 
     case(7)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_67um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_67um_Std**2 
     case DEFAULT
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
   end select

   if (idiag_output == symbol%YES) then
           print *, "ACHA_MODE = ", Acha_Mode_Flag
           print *, "y = ", y
   endif
   
   !-------------------------------------------------------------------
   ! Determine surface type for use in forward model
   ! 0 = Water
   ! 1 = Land
   ! 2 = Snow
   ! 3 = Desert
   ! 4 = Arctic
   ! 5 = Antarctic
   !-------------------------------------------------------------------
   call DETERMINE_SFC_TYPE_FORWARD_MODEL(symbol, &
                                         Input%Surface_Type(Elem_Idx,Line_Idx), &
                                         Input%Snow_Class (Elem_Idx,Line_Idx), &
                                         Input%Latitude(Elem_Idx,Line_Idx), &
                                         Input%Surface_Emissivity_39um(Elem_Idx,Line_Idx), &
                                         Sfc_Type_Forward_Model)

   !-------------------------------------------------------------------
   ! Based on fm surface type, set the clear-sky covariance terms
   !-------------------------------------------------------------------
   call SET_CLEAR_SKY_COVARIANCE_TERMS(Sfc_Type_Forward_Model)

   !-------------------------------------------------------------------
   ! These values are used in the baseline code.  In the latest
   ! code, the values from SET_CLEAR_SKY_COVARIANCE_TERMS are used
   !
   ! these are based on patmos-x clear data and are the 
   ! approx average of des+asc from August 2006 NOAA-18
   !-------------------------------------------------------------------
   if (Input%Surface_Type(Elem_Idx,Line_Idx) == symbol%WATER_SFC) then
     T11um_Clr_Uncer = T11um_Clr_Uncer_Water
     T11um_67um_Clr_Uncer = T11um_67um_Clr_Uncer_Water
     T11um_85um_Clr_Uncer = T11um_85um_Clr_Uncer_Water
     T11um_12um_Clr_Uncer = T11um_12um_Clr_Uncer_Water
     T11um_133um_Clr_Uncer = T11um_133um_Clr_Uncer_Water
   else
     T11um_Clr_Uncer = T11um_Clr_Uncer_land
     T11um_67um_Clr_Uncer = T11um_67um_Clr_Uncer_land
     T11um_85um_Clr_Uncer = T11um_85um_Clr_Uncer_land
     T11um_12um_Clr_Uncer = T11um_12um_Clr_Uncer_land
     T11um_133um_Clr_Uncer = T11um_133um_Clr_Uncer_land
   endif

   if (idiag_output == symbol%YES) then
           print *, "Clear-sky uncertainties = ", &
           T11um_Clr_Uncer,T11um_85um_Clr_Uncer,T11um_12um_Clr_Uncer, &
           T11um_133um_Clr_Uncer
   endif

  !--------------------------------------------------------------------
  ! pick a priori conditions
  !--------------------------------------------------------------------

  !--- logic for unmasked or untyped pixels (UndetOutput%Ected cloud)
  if (Undetected_Cloud == symbol%YES) then
         if (Tc_Opaque_Level < 260.0 .and.  &
             Tc_Opaque_Level /= Missing_Value_Real4) then
             Cloud_Type = symbol%CIRRUS_TYPE
         else
             Cloud_Type = symbol%FOG_TYPE
         endif
  endif

  !---- Compute 11um emissivity referenced to tropopause
  Emiss_11um_Tropo = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( &
                             Tropo_Level_Rtm, &
                             Input%Rad_11um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_11um_Local(Elem_Idx,Line_Idx))

  !---- select Output%Tc and Output%Ec apriori based on cloud type

  if ((ilrc /= Missing_Value_Int4) .and. &
      (jlrc /= Missing_Value_Int4)) then
           Bt_11um_Lrc =  Input%Bt_11um(ilrc,jlrc)
  else
           Bt_11um_Lrc = Missing_Value_Real4
  endif

  call COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO( &
                       symbol, &
                       Cloud_Phase, &
                       Emiss_11um_Tropo, &
                       Input%Tropopause_Temperature(Elem_Idx,Line_Idx), &
                        Input%Bt_11um(Elem_Idx,Line_Idx), &
                        Bt_11um_Lrc, &
                       Tc_Opaque_Level, &
                       Input%Cosine_Zenith_Angle(Elem_Idx,Line_Idx), &
                       Tc_Ap,Tc_Ap_Uncer, &
                       Ec_Ap,Ec_Ap_Uncer, &
                       Beta_Ap,Beta_Ap_Uncer)

   x_Ap(1) = Tc_Ap
   x_Ap(2) = Ec_Ap
   x_Ap(3) = Beta_Ap

   if (idiag_output == symbol%YES) then
           print *, "x_Ap = ", x_Ap
           print *, "x_Ap Uncer = ", Tc_Ap_Uncer,Ec_Ap_Uncer, Beta_Ap_Uncer
   endif

   !---- determine Output%Beta fit parameters based on phase (derived from type)
          
   !--- water phase clouds
   if (Cloud_Phase == symbol%WATER_PHASE) then
      a_Beta_11um_133um_fit = a_Beta_11um_133um_fit_Water
      b_Beta_11um_133um_fit = b_Beta_11um_133um_fit_Water
      a_Beta_11um_85um_fit = a_Beta_11um_85um_fit_Water
      b_Beta_11um_85um_fit = b_Beta_11um_85um_fit_Water
      a_Beta_11um_67um_fit = a_Beta_11um_67um_fit_Water
      b_Beta_11um_67um_fit = b_Beta_11um_67um_fit_Water
   else
   !--- ice phase clouds
      a_Beta_11um_133um_fit = a_Beta_11um_133um_fit_Ice
      b_Beta_11um_133um_fit = b_Beta_11um_133um_fit_Ice
      a_Beta_11um_85um_fit = a_Beta_11um_85um_fit_Ice
      b_Beta_11um_85um_fit = b_Beta_11um_85um_fit_Ice
      a_Beta_11um_67um_fit = a_Beta_11um_67um_fit_Ice
      b_Beta_11um_67um_fit = b_Beta_11um_67um_fit_Ice
  endif
  if (idiag_output == symbol%YES) then
        print *, "Output%Beta fit for 6.7 = ", a_Beta_11um_67um_fit, b_Beta_11um_67um_fit
        print *, "Output%Beta fit for 8.5 = ", a_Beta_11um_85um_fit, b_Beta_11um_85um_fit
        print *, "Output%Beta fit for 13.3 = ", a_Beta_11um_133um_fit, b_Beta_11um_133um_fit
  endif

  !--- now compute Sa
  Sa = 0.0
  Sa(1,1) = Tc_Ap_Uncer
  Sa(2,1) = 0.0 
  Sa(2,2) = Ec_Ap_Uncer
  Sa(1,2) = 0.0 
  Sa(3,3) = Beta_Ap_Uncer

  !--- modify a priori values based on lrc
  if ((ilrc /= Missing_Value_Int4) .and. &
      (jlrc /= Missing_Value_Int4)) then

        if ((Output%Tc(ilrc,jlrc) /= Missing_Value_Real4) .and. &
            (Output%Ec(ilrc,jlrc) > 0.00) .and. &
            (Output%Ec(ilrc,jlrc) <= 1.0)) then

          !-- use lrc value but weight uncertainty
          x_Ap(1) = Output%Tc(ilrc,jlrc)
          Sa(1,1) = 5.0 + (1.0-Output%Ec(ilrc,jlrc))*Tc_Ap_Uncer

       endif
  endif

  !--- square the individual elements to convert to variances (not a matmul)
  Sa = Sa**2

  !--- compute inverse of Sa matrix
  Singular_Flag =  INVERT_MATRIX(Sa, Sa_Inv, Num_Param)
  if (Singular_Flag == 1) then
    print *, "Cloud Height warning ==> Singular Sa in Split Window", &
           Elem_Idx,Line_Idx, Cloud_Type
    Fail_Flag = symbol%YES
    exit
  endif

!----------------------------------------------------------------------------
! compute clear-sky radiances
!-----------------------------------------------------------------------------

 !---------------------------------------------------------------------------
 !--- modify clear radiances to simulate that from an opaque cloud at 200 mb
 !--- above the surface when a multi-layer situation is suspOutput%Ected
 !---------------------------------------------------------------------------
 if (Cloud_Type == symbol%OVERLAP_TYPE .and.  &
    Output%Lower_Cloud_Pressure(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then

    Meta_Data_Flags(6) = symbol%YES

    Rad_Ac_11um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                            Hght_Prof_Rtm,Atm_Rad_Prof_11um_Rtm)

    Trans_Ac_11um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                            Hght_Prof_Rtm,Atm_Trans_Prof_11um_Rtm)

    Bc_11um = PLANCK_RAD_FAST(Input%Chan_Idx_11um,Output%Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
    Rad_Clear_11um = Rad_Ac_11um + Trans_Ac_11um*Bc_11um

  if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 5) then
     Rad_Ac_12um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Rad_Prof_12um_Rtm)

     Trans_Ac_12um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Trans_Prof_12um_Rtm)
     Bc_12um = PLANCK_RAD_FAST(Input%Chan_Idx_12um,Output%Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
     Rad_Clear_12um = Rad_Ac_12um + Trans_Ac_12um*Bc_12um
  endif

  if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 6) then
     Rad_Ac_133um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Rad_Prof_133um_Rtm)

     Trans_Ac_133um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Trans_Prof_133um_Rtm)

     Bc_133um = PLANCK_RAD_FAST(Input%Chan_Idx_133um,Output%Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
     Rad_Clear_133um = Rad_Ac_133um + Trans_Ac_133um*Bc_133um
  endif

  if (Acha_Mode_Flag == 4) then
     Rad_Ac_85um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Rad_Prof_85um_Rtm)

     Trans_Ac_85um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Trans_Prof_85um_Rtm)

     Bc_85um = PLANCK_RAD_FAST(Input%Chan_Idx_85um,Output%Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
     Rad_Clear_85um = Rad_Ac_85um + Trans_Ac_85um*Bc_85um
  endif

  if (Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 7) then
     Rad_Ac_67um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Rad_Prof_67um_Rtm)

     Trans_Ac_67um = GENERIC_PROFILE_INTERPOLATION(Output%Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Trans_Prof_67um_Rtm)

     Bc_67um = PLANCK_RAD_FAST(Input%Chan_Idx_67um,Output%Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
     Rad_Clear_67um = Rad_Ac_67um + Trans_Ac_67um*Bc_67um
  endif


  Tsfc_Est = Output%Lower_Cloud_Temperature(Elem_Idx,Line_Idx)

 else

 !---------------------------------------------------------------
 ! if not multi-layer, use existing clear sky radiances
 !---------------------------------------------------------------

  Tsfc_Est = Input%Surface_Temperature(Elem_Idx,Line_Idx)
  Rad_Clear_11um = Input%Rad_Clear_11um_Local(Elem_Idx,Line_Idx)
  if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag ==3 .or. Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 5) then
      Rad_Clear_12um = Input%Rad_Clear_12um_Local(Elem_Idx,Line_Idx)
  endif
  if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 6) then
      Rad_Clear_133um = Input%Rad_Clear_133um_Local(Elem_Idx,Line_Idx)
  endif
  if (Acha_Mode_Flag == 4) then
      Rad_Clear_85um = Input%Rad_Clear_85um_Local(Elem_Idx,Line_Idx)
  endif
  if (Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 7) then
      Rad_Clear_67um = Input%Rad_Clear_67um_Local(Elem_Idx,Line_Idx)
  endif

 endif

 if (idiag_output == symbol%YES) then
    print *, "Clear radiances = ", Rad_Clear_67um, Rad_Clear_85um, Rad_Clear_11um, Rad_Clear_12um, Rad_Clear_133um
 endif


!----------------------------------------------------------------
! Determine the level of the highest inversion (0=if none)
!----------------------------------------------------------------
Inver_Level_Rtm = DETERMINE_INVERSION_LEVEL(Tropo_Level_Rtm, Sfc_Level_Rtm, Input%Surface_Air_Temperature(Elem_Idx,Line_Idx))

!----------------------------------------------------------------
! if no inversion is present, check to see if clear radiance
! is warmer than observed radiance, if not then do not proceed
!----------------------------------------------------------------
if (( Inver_Level_Rtm == 0) .and. (Rad_Clear_11um < Input%Rad_11um(Elem_Idx,Line_Idx))) then

   !-- consider this pixel a failed retrieval
   Fail_Flag = symbol%YES

   Meta_Data_Flags(1) = symbol%NO

   !--- assign default values for this result
   Output%Tc(Elem_Idx,Line_Idx) = x_Ap(1)   !Missing_Value_Real4
   Output%Ec(Elem_Idx,Line_Idx) = x_Ap(2)   !Missing_Value_Real4
   Output%Beta(Elem_Idx,Line_Idx) = x_Ap(3) !Missing_Value_Real4
   Output%Tau(Elem_Idx,Line_Idx) = Missing_Value_Real4
   Output%Reff(Elem_Idx,Line_Idx) = Missing_Value_Real4
   Output%Qf(Elem_Idx,Line_Idx) = 1

   call KNOWING_T_COMPUTE_P_Z(symbol, &
                              Cloud_Type, &
                              Output%Pc(Elem_Idx,Line_Idx), &
                              Output%Tc(Elem_Idx,Line_Idx), &
                              Output%Zc(Elem_Idx,Line_Idx), &
                              Ilev,ierror,NWP_Profile_Inversion_Flag)
   if (ierror == 0) then
    Output%Zc(Elem_Idx,Line_Idx) = max(0.0,Output%Zc(Elem_Idx,Line_Idx))
    Output%Pc(Elem_Idx,Line_Idx) = min(Input%Surface_Pressure(Elem_Idx,Line_Idx), &
                                Output%Pc(Elem_Idx,Line_Idx))
   endif
   
   !-- go to next pixel
   cycle

endif

!-----------------------------------------------------------------
! start of retrieval loop
!-----------------------------------------------------------------
Iter_Idx = 0
Converged_Flag = symbol%NO
Fail_Flag = symbol%NO

!---- assign x to the first guess
x = x_Ap

Retrieval_Loop: do

 Iter_Idx = Iter_Idx + 1

  !---------------------------------------------------------------------
  ! estimate above cloud radiances and transmissions
  !---------------------------------------------------------------------
  Tc_temp = x(1)

  call KNOWING_T_COMPUTE_P_Z(symbol,Cloud_Type,Pc_temp,Tc_temp,Zc_temp,Ilev,ierror,NWP_Profile_Inversion_Flag)

  !--- compute above-cloud terms
  Rad_Ac_11um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_11um_Rtm)

  Trans_Ac_11um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_11um_Rtm)

  if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 5) then
     Rad_Ac_12um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_12um_Rtm)

     Trans_Ac_12um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_12um_Rtm)
  endif

  if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 6) then
    Rad_Ac_133um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_133um_Rtm)

    Trans_Ac_133um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_133um_Rtm)
  endif

  if (Acha_Mode_Flag == 4) then
    Rad_Ac_85um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_85um_Rtm)

    Trans_Ac_85um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_85um_Rtm)
  endif

  if (Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 7) then
    Rad_Ac_67um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_67um_Rtm)

    Trans_Ac_67um = GENERIC_PROFILE_INTERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_67um_Rtm)

  endif

  !--------------------------------------------------
  ! call abi_planck_temp(ichan,rad,bt)
  !--------------------------------------------------
   if (Input%Chan_On_67um == symbol%YES)  Bt_Clear_67um = PLANCK_TEMP_FAST(Input%Chan_Idx_67um, Rad_Clear_67um)
   if (Input%Chan_On_85um == symbol%YES)  Bt_Clear_85um = PLANCK_TEMP_FAST(Input%Chan_Idx_85um, Rad_Clear_85um)
   if (Input%Chan_On_11um == symbol%YES)  Bt_Clear_11um = PLANCK_TEMP_FAST(Input%Chan_Idx_11um, Rad_Clear_11um)
   if (Input%Chan_On_12um == symbol%YES)  Bt_Clear_12um = PLANCK_TEMP_FAST(Input%Chan_Idx_12um, Rad_Clear_12um)
   if (Input%Chan_On_133um == symbol%YES)  Bt_Clear_133um = PLANCK_TEMP_FAST(Input%Chan_Idx_133um, Rad_Clear_133um)
 
  !--------------------------------------------------
  ! call forward models
  !--------------------------------------------------
  if (idiag_output == symbol%YES) then
      print *, "Iter_Idx = ", Iter_Idx
  endif
  call COMPUTE_FORWARD_MODEL_AND_KERNEL(symbol,Acha_Mode_Flag,  &
           Input%Chan_On_67um,  &
           Input%Chan_On_85um, &
           Input%Chan_On_11um,  &
           Input%Chan_On_12um, &
           Input%Chan_On_133um, &
           Input%Chan_Idx_67um,  &
           Input%Chan_Idx_85um, &
           Input%Chan_Idx_11um,  &
           Input%Chan_Idx_12um, &
           Input%Chan_Idx_133um, &
           x, &
           Rad_Clear_67um, Rad_Ac_67um, Trans_Ac_67um, &
           Rad_Clear_85um, Rad_Ac_85um, Trans_Ac_85um, &
           Rad_Clear_11um, Rad_Ac_11um, Trans_Ac_11um, &
           Rad_Clear_12um, Rad_Ac_12um, Trans_Ac_12um, &
           Rad_Clear_133um, Rad_Ac_133um, Trans_Ac_133um, &
           a_Beta_11um_133um_fit, b_Beta_11um_133um_fit, &
           a_Beta_11um_85um_fit, b_Beta_11um_85um_fit, &
           a_Beta_11um_67um_fit, b_Beta_11um_67um_fit, &
           f, K,Emiss_Vector,idiag_output)

  !--------------------------------------------------
  ! compute the Sy convariance matrix
  !--------------------------------------------------

   Call COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE(  &
                                                 symbol, &
                                                 Emiss_Vector, &
                                                 Acha_Mode_Flag, &
                                                 y_variance, &
                                                 Sy) 


  if (idiag_output == symbol%YES) then
          print *, "Sa1 = ", Sa(1,1), Sa(1,2), Sa(1,3)
          print *, "Sa2 = ", Sa(2,1), Sa(2,2), Sa(2,3)
          print *, "Sa3 = ", Sa(3,1), Sa(3,2), Sa(3,3)
          print *, "Sy1 = ", Sy(1,:)
          print *, "Sy2 = ", Sy(2,:)
          print *, "shape of Sy = ", shape(Sy), Num_Obs
  endif
  !--------------------------------------------------
  ! call OE routine to advance the Iteration
  !--------------------------------------------------
  call OPTIMAL_ESTIMATION_ACHA(symbol,Iter_Idx,Iter_Idx_Max,Num_Param,Num_Obs, &
                         conv_crit,Delta_X_Max, &
                         y,f,x,x_Ap,K,Sy,Sa_inv, &
                         Sx,AKM,Delta_X,Converged_Flag,Fail_Flag, &
                         idiag_output)

! Diag_Pix_Array_1(Elem_Idx,Line_Idx) = AKM(1,1)
! Diag_Pix_Array_2(Elem_Idx,Line_Idx) = AKM(2,2)
! Diag_Pix_Array_3(Elem_Idx,Line_Idx) = AKM(3,3)

! Diag_Pix_Array_1(Elem_Idx,Line_Idx) = f(1) - y(1)
! Diag_Pix_Array_2(Elem_Idx,Line_Idx) = f(2) - y(2)
! Diag_Pix_Array_3(Elem_Idx,Line_Idx) = f(3) - y(3)

! Diag_Pix_Array_1(Elem_Idx,Line_Idx) = AKM(1,1)
! Diag_Pix_Array_2(Elem_Idx,Line_Idx) = f(1) - y(1)
! Diag_Pix_Array_3(Elem_Idx,Line_Idx) = x(1) - x_ap(1)

  !--- check for a failed Iteration
  if (Fail_Flag == symbol%YES) then
     exit
  endif

  !---------------------------------------------------------
  ! update retrieved vOutput%Ector
  !---------------------------------------------------------
  x = x + Delta_X

  if (idiag_output==symbol%YES) then  
          print *, "x = ", x
  endif

  !--------------------------------------------------------
  ! exit retrieval loop if converged
  !--------------------------------------------------------
  if (Converged_Flag == symbol%YES) then
       if (idiag_output==symbol%YES) then  
             print *, "convergence acheived ", Iter_Idx, x
             print *, '  '
       endif
       exit
  endif

  !-------------------------------------------------------
  ! constrain to reasonable values
  !-------------------------------------------------------
  x(1) = max(min_allowable_Tc,min(Tsfc_Est,x(1)))     !should we do this?
  x(2) = max(0.0,min(x(2),1.0))
  x(3) = max(0.8,min(x(3),1.8))

end do Retrieval_Loop

!=================================================================
! Begin Retrieval Post Processing
!=================================================================

!-----------------------------------------------------------------
! Successful Retrieval Post Processing
!-----------------------------------------------------------------
if (Fail_Flag == symbol%NO) then  !successful retrieval if statement

 !--- save retrieved vOutput%Ector into its output variables
 Output%Tc(Elem_Idx,Line_Idx) = x(1)
 Output%Beta(Elem_Idx,Line_Idx) = x(3)

 !--- save nadir adjusted emissivity and optical depth
 if (x(2) < 1.00) then

   call COMPUTE_TAU_REFF_ACHA(symbol, &
                              Output%Beta(Elem_Idx,Line_Idx), &
                              Input%Cosine_Zenith_Angle(Elem_Idx,Line_Idx), &
                              Cloud_Phase, &
                              x(2), &
                              Output%Ec(Elem_Idx,Line_Idx), &
                              Output%Tau(Elem_Idx,Line_Idx), &
                              Output%Reff(Elem_Idx,Line_Idx))

 else

   Output%Tau(Elem_Idx,Line_Idx) = 20.0
   Output%Ec(Elem_Idx,Line_Idx) = 1.0
   Output%Beta(Elem_Idx,Line_Idx) = 1.3
   Output%Reff(Elem_Idx,Line_Idx) = 20.0

 endif


 !--- save uncertainty estimates
 Output%Tc_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(1,1))
 Output%Ec_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(2,2))
 Output%Beta_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(3,3))

 !--- set quality flag for a successful retrieval
 Output%Qf(Elem_Idx,Line_Idx) = 3


 !---  for low clouds over water, force fixed lapse rate estimate of height
 Delta_Cld_Temp_Sfc_Temp = Input%Surface_Temperature(Elem_Idx,Line_Idx) - Output%Tc(Elem_Idx,Line_Idx)

 if (Input%Surface_Type(Elem_Idx,Line_Idx) == symbol%WATER_SFC .and. &
     Input%Snow_Class (Elem_Idx,Line_Idx) == symbol%NO_SNOW .and. &
     ((Delta_Cld_Temp_Sfc_Temp <  MAX_DELTA_T_INVERSION) .or. &
      (Cloud_Type == sym%WATER_TYPE) .or. &
      (Cloud_Type == sym%FOG_TYPE))) then

       !-- select lapse rate  (k/km)
       Lapse_Rate =  -0.061  + &
                     1.67*Delta_Cld_Temp_Sfc_Temp   &
                     -0.124*(Delta_Cld_Temp_Sfc_Temp**2)   &
                     +0.00343*(Delta_Cld_Temp_Sfc_Temp**3)

       Lapse_Rate = Lapse_Rate / 1000.0  !(K/m)

       !-- compute height
       Output%Zc(Elem_Idx,Line_Idx) = Delta_Cld_Temp_Sfc_Temp/Lapse_Rate + Input%Surface_Elevation(Elem_Idx,Line_Idx)

       !--- compute pressure
       call KNOWING_Z_COMPUTE_T_P(symbol,Output%Pc(Elem_Idx,Line_Idx),R4_Dummy,Output%Zc(Elem_Idx,Line_Idx),Ilev)

       !--- set meta data flag
       Meta_Data_Flags(7) = symbol%YES


 else   !the general top-down solution

        !--- Estimate height and pressure
        call KNOWING_T_COMPUTE_P_Z(symbol,Cloud_Type,Output%Pc(Elem_Idx,Line_Idx), &
                                  Output%Tc(Elem_Idx,Line_Idx), &
                                  Output%Zc(Elem_Idx,Line_Idx),&
                                  Ilev,ierror,NWP_Profile_Inversion_Flag)

        !--- check for NWP profile inversion and set meta data flag.
        if (NWP_Profile_Inversion_Flag == 1) then
          Meta_Data_Flags(8) = symbol%YES
        endif

 endif


 !--- compute height and pressure uncertainties 

 !--- compute lapse rate as dT/dZ
 Lapse_Rate =  (Temp_Prof_Rtm(Ilev+1) - Temp_Prof_Rtm(Ilev)) / &
               (Hght_Prof_Rtm(Ilev+1) - Hght_Prof_Rtm(Ilev))

 !--- compute lapse rate as dP/dZ
 Lapse_Rate_dP_dZ =  (Press_Prof_Rtm(Ilev+1) - Press_Prof_Rtm(Ilev)) / &
                     (Hght_Prof_Rtm(Ilev+1) - Hght_Prof_Rtm(Ilev))
 !-- Compute Height Uncertainty
 if (Lapse_Rate /= 0.0) THEN
      Output%Zc_Uncertainty(Elem_Idx,Line_Idx) = Output%Tc_Uncertainty(Elem_Idx,Line_Idx) / ABS(Lapse_Rate)
 endif

 !-- Compute Pressure Uncertainty
 if (Lapse_Rate /= 0.0 .and. Lapse_Rate_dP_dZ /= 0.0) THEN
     Output%Pc_Uncertainty(Elem_Idx,Line_Idx) = Output%Zc_Uncertainty(Elem_Idx,Line_Idx) * ABS(Lapse_Rate_dP_dZ)
 endif
 

 !--- quality flags of the retrieved parameters
 DO Param_Idx = 1,Num_Param    !loop over parameters
       if (Sx(Param_Idx,Param_Idx) < 0.111*Sa(Param_Idx,Param_Idx) ) THEN
            Output%OE_Qf(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_1_3_APRIORI_RETREVIAL
       elseif (Sx(Param_Idx,Param_Idx) < 0.444*Sa(Param_Idx,Param_Idx)) THEN
            Output%OE_Qf(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_2_3_APRIORI_RETREVIAL
       else
            Output%OE_Qf(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_LOW_QUALITY_RETREVIAL
       endif
 enddo
 
else   

!-----------------------------------------------------------------
! Failed Retrieval Post Processing
!-----------------------------------------------------------------

 !--- set output variables to apriori
 Output%Tc(Elem_Idx,Line_Idx) = x_Ap(1)   !Missing_Value_Real4
 Output%Ec(Elem_Idx,Line_Idx) = x_Ap(2)   !Missing_Value_Real4
 Output%Beta(Elem_Idx,Line_Idx) = x_Ap(3) !Missing_Value_Real4

 !--- set derived parameters to missing
 Output%Tau(Elem_Idx,Line_Idx) = Missing_Value_Real4
 Output%Reff(Elem_Idx,Line_Idx) = Missing_Value_Real4
 Output%Pc(Elem_Idx,Line_Idx) = Missing_Value_Real4
 Output%Zc(Elem_Idx,Line_Idx) = Missing_Value_Real4

 !--- set quality flags
 Output%OE_Qf(:,Elem_Idx,Line_Idx) = 0
 Output%Qf(Elem_Idx,Line_Idx) = 1

 !--- estimate height and pressure
 call KNOWING_T_COMPUTE_P_Z(symbol,Cloud_Type,Output%Pc(Elem_Idx,Line_Idx), &
                            Output%Tc(Elem_Idx,Line_Idx), &
                            Output%Zc(Elem_Idx,Line_Idx), &
                            Ilev,ierror,NWP_Profile_Inversion_Flag)

endif                              !end successful retrieval if statement

!------- determine cloud layer based on pressure
Output%Cloud_Layer(Elem_Idx,Line_Idx) = 0
if (Output%Pc(Elem_Idx,Line_Idx) <= 440.0) then
   Output%Cloud_Layer(Elem_Idx,Line_Idx) = 3
elseif (Output%Pc(Elem_Idx,Line_Idx) < 680.0) then
   Output%Cloud_Layer(Elem_Idx,Line_Idx) = 2
else
   Output%Cloud_Layer(Elem_Idx,Line_Idx) = 1
endif

!--- if retrieval done for an undetected pixel, label the Output%Qf
if (Undetected_Cloud == symbol%YES) then
 Output%Qf(Elem_Idx,Line_Idx) = 2
 Output%Cloud_Layer(Elem_Idx,Line_Idx) = 0
endif

!-----------------------------------------------------------------
! End Retrieval Post Processing
!-----------------------------------------------------------------

endif     ! ---------- end of data check
 
 
 !-----------------------------------------------------------------------------
 !--- Cloud Base and Top
 !---
 !--- Note 1. Extinction values are in km^(-1)
 !--- Note 2. All heights and thickness are converted to meters
 !-----------------------------------------------------------------------------
 if (Output%Zc(Elem_Idx,Line_Idx) /= Missing_Value_Real4 .and. &
     Output%Tau(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then

       Cloud_Extinction = WATER_EXTINCTION

       if(Cloud_Type == symbol%OPAQUE_ICE_TYPE .or. &
          Cloud_Type == symbol%OVERSHOOTING_TYPE) then
          Itemp = int(Output%Tc(Elem_Idx,Line_Idx))
          select case (Itemp)
            case (:199)    ; Cloud_Extinction = ICE_EXTINCTION1
            case (200:219) ; Cloud_Extinction = ICE_EXTINCTION2
            case (220:239) ; Cloud_Extinction = ICE_EXTINCTION3
            case (240:259) ; Cloud_Extinction = ICE_EXTINCTION4
            case (260:)    ; Cloud_Extinction = ICE_EXTINCTION5
          end select
       endif

       if(Cloud_Type == symbol%CIRRUS_TYPE .or. &
          Cloud_Type == symbol%OVERLAP_TYPE) then
          Itemp = int(Output%Tc(Elem_Idx,Line_Idx))
          select case (Itemp)
            case (:199)    ; Cloud_Extinction = CIRRUS_EXTINCTION1
            case (200:219) ; Cloud_Extinction = CIRRUS_EXTINCTION2
            case (220:239) ; Cloud_Extinction = CIRRUS_EXTINCTION3
            case (240:259) ; Cloud_Extinction = CIRRUS_EXTINCTION4
            case (260:)    ; Cloud_Extinction = CIRRUS_EXTINCTION5
          end select
       endif

       Cloud_Geometrical_Thickness = Output%Tau(Elem_Idx,Line_Idx) / Cloud_Extinction   !(km)
       Cloud_Geometrical_Thickness = Cloud_Geometrical_Thickness * 1000.0 !(m)

       if (Output%Tau(Elem_Idx,Line_Idx) < 2.0) then 
          Cloud_Geometrical_Thickness_Top_Offset = Cloud_Geometrical_Thickness/2.0    !(m)
       else
          Cloud_Geometrical_Thickness_Top_Offset = 1000.0 / Cloud_Extinction !(m)
       endif

       Zc_Top_Max = Hght_Prof_Rtm(Tropo_Level_Rtm)
       Zc_Base_Min = Hght_Prof_Rtm(Sfc_Level_Rtm)

       !-- new code
       if (Output%Zc(Elem_Idx,Line_Idx) > Zc_Top_Max .and. Cloud_Type == symbol%OVERSHOOTING_TYPE) then
          Output%Zc_Top(Elem_Idx,Line_Idx) = Output%Zc(Elem_Idx,Line_Idx)
       else
          Output%Zc_Top(Elem_Idx,Line_Idx) = min(Zc_Top_Max, &
                                             Output%Zc(Elem_Idx,Line_Idx) +  &
                                             Cloud_Geometrical_Thickness_Top_Offset)
       endif

       Output%Zc_Base(Elem_Idx,Line_Idx) = min(Output%Zc(Elem_Idx,Line_Idx), &
                                           max(Zc_Base_Min,  &
                                           Output%Zc_Top(Elem_Idx,Line_Idx) - Cloud_Geometrical_Thickness))

 endif

 !------------------------------------------------------------------------
 ! Pack Quality Flags Output
 !------------------------------------------------------------------------
 !--- bit1  
 Output%Packed_Qf(Elem_Idx,Line_Idx) =  1_int1

 !--- bit2
 if (Output%OE_Qf(1,Elem_Idx,Line_Idx)  /= CTH_PARAM_FAILED_RETREVIAL)  then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 2_int1
 endif

 !--- bit3
 if (Output%OE_Qf(2,Elem_Idx,Line_Idx)  /= CTH_PARAM_FAILED_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 4_int1
 endif

 !--- bit4
 if (Output%OE_Qf(3,Elem_Idx,Line_Idx)  /= CTH_PARAM_FAILED_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 8_int1
 endif

 !--- bit5
 if (Output%OE_Qf(1,Elem_Idx,Line_Idx)  == CTH_PARAM_LOW_QUALITY_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 16_int1
 endif

 !--- bit6
 if (Output%OE_Qf(2,Elem_Idx,Line_Idx)  == CTH_PARAM_LOW_QUALITY_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 32_int1
 endif

 !--- bit7
 if (Output%OE_Qf(3,Elem_Idx,Line_Idx)  == CTH_PARAM_LOW_QUALITY_RETREVIAL) then
     Output%Packed_Qf(Elem_Idx,Line_Idx) =    &
                     Output%Packed_Qf(Elem_Idx,Line_Idx) + 64_int1
 endif

 !------------------------------------------------------------------------
 ! Pack Meta Data for Output
 !------------------------------------------------------------------------
 do i = 1, 8
    Output%Packed_Meta_Data(Elem_Idx,Line_Idx) = Output%Packed_Meta_Data(Elem_Idx,Line_Idx) +  &
                                          (2**(i-1)) * Meta_Data_Flags(i)
 enddo


!---- null profile pointers each time - WCS3

   CALL NULL_PIX_POINTERS(Input, Acha_RTM_NWP)

 end do Element_Loop

end do Line_loop

end do pass_loop

!------------------------------------------------------------------------
! Apply Parallax Correction 
!------------------------------------------------------------------------
  call PARALLAX_ACHA(Output%Zc, Input%Surface_Elevation, &
                     Input%Latitude, Input%Longitude, &
                     Input%Sensor_Zenith_Angle, &
                     Input%Sensor_Azimuth_Angle, &
                     Output%Latitude_Pc,&
                     Output%Longitude_Pc) 

  !--- deallocate 2D arrays
  if (allocated(Elem_Idx_LRC)) deallocate(Elem_Idx_LRC)
  if (allocated(Line_Idx_LRC)) deallocate(Line_Idx_LRC)
  if (allocated(Skip_LRC_Mask)) deallocate(Skip_LRC_Mask)

  !--- deallocate 1D-VAR arrays
  deallocate(y)
  deallocate(y_variance)
  deallocate(f)
  deallocate(x)
  deallocate(x_Ap)
  deallocate(Delta_X)
  deallocate(K)
  deallocate(Sa)
  deallocate(Sa_inv)
  deallocate(Sx)
  deallocate(Sx_inv)
  deallocate(E)
  deallocate(Sy)
  deallocate(Sy_inv)
  deallocate(Emiss_Vector)
  deallocate(AKM)

end subroutine  AWG_CLOUD_HEIGHT_ALGORITHM

!-------------------------------------------------------------------------------
! Routine to spatially interpret water cloud temperature values to surrounding 
! pixels
!
! input:  interp_flag = 0 (do no interp, assume Zc=2km) /= 0 (do spatial interp)
!-------------------------------------------------------------------------------
  subroutine SPATIALLY_INTERPOLATE_LOWER_CLOUD_POSITION(symbol, &
                                                        Interp_Flag,Line_Idx_Min, &
                                                        Number_Of_Lines, &
                                                        Elem_Idx_Nwp,Line_Idx_Nwp, &
                                                        Invalid_Data_Mask, &
                                                        Cloud_Type_Local, &
                                                        Surface_Pressure, &
                                                        Cloud_Pressure, &
                                                        Lower_Cloud_Pressure, &
                                                        Lower_Cloud_Temperature, &
                                                        Lower_Cloud_Height)

      type(symbol_acha), intent(in) :: symbol
      integer, intent(in):: Interp_Flag    
      integer, intent(in):: Line_Idx_Min
      integer, intent(in):: Number_Of_Lines
      integer, intent(in), dimension(:,:):: Elem_Idx_Nwp
      integer, intent(in), dimension(:,:):: Line_Idx_Nwp
      integer(kind=int1), intent(in), dimension(:,:):: Invalid_Data_Mask
      integer(kind=int1), intent(in), dimension(:,:):: Cloud_Type_Local
      real, intent(in), dimension(:,:):: Surface_Pressure
      real, intent(in), dimension(:,:):: Cloud_Pressure
      real, intent(out), dimension(:,:):: Lower_Cloud_Pressure
      real, intent(out), dimension(:,:):: Lower_Cloud_Temperature
      real, intent(out), dimension(:,:):: Lower_Cloud_Height
      integer:: Line_Idx
      integer:: Elem_Idx
      integer:: Inwp
      integer:: Jnwp
      integer:: Inwp_x
      integer:: Jnwp_x
      integer:: Ilev
      integer:: Num_Elem
      integer:: num_line
      integer:: Line_start
      integer:: Line_end
      integer:: delem
      integer:: dline
      integer:: Elem_Idx_1
      integer:: Elem_Idx_2
      integer:: Line_Idx_1
      integer:: Line_Idx_2
      integer:: count_Valid
      integer, dimension(:,:), allocatable:: mask

      Num_Elem = size(Cloud_Type_Local,1)      !-----
      num_Line = size(Cloud_Type_Local,2)      !-----
      Line_start = Line_Idx_min
      Line_end = Number_Of_Lines + Line_Idx_min - 1

      if (Interp_Flag == 0) then 

       Line_loop_1: do Line_Idx = Line_start, Line_end 
         Element_Loop_1: do Elem_Idx = 1, Num_Elem 
           if (Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) then
               Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Missing_Value_Real4
               cycle
           endif

          Inwp = Elem_Idx_Nwp(Elem_Idx,Line_Idx)
          Jnwp = Line_Idx_Nwp(Elem_Idx,Line_Idx)

           if (Inwp == Missing_Value_Int4 .or. &
               Jnwp == Missing_Value_Int4) THEN
               
               cycle
               
           endif
          Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Surface_Pressure(Elem_Idx,Line_Idx) - Pc_Lower_Cloud_offset
        end do Element_Loop_1
       end do Line_loop_1


      else    !if Interp_Flag > 0 then do a spatial interpolation

        !--- set box width
        delem = 1
        dline = 1

        allocate(mask(Num_Elem,num_line))

         mask = 0

         where((Cloud_Type_Local == symbol%FOG_TYPE .or. &
             Cloud_Type_Local == symbol%WATER_TYPE .or. &
             Cloud_Type_Local == symbol%SUPERCOOLED_TYPE) .and. &
             Cloud_Pressure /= Missing_Value_Real4)
             mask = 1
         endwhere

         Line_loop_2: do Line_Idx = Line_start, Line_end 
             Element_Loop_2: do Elem_Idx = 1, Num_Elem 

             if (Cloud_Type_Local(Elem_Idx,Line_Idx)  /= symbol%OVERLAP_TYPE) then
                 cycle
             endif

             Elem_Idx_1 = min(Num_Elem,max(1,Elem_Idx-delem))
             Elem_Idx_2 = min(Num_Elem,max(1,Elem_Idx+delem))
             Line_Idx_1 = min(num_line,max(1,Line_Idx-dline))
             Line_Idx_2 = min(num_line,max(1,Line_Idx+dline))

             count_Valid = sum(mask(Elem_Idx_1:Elem_Idx_2,Line_Idx_1:Line_Idx_2))

             if (count_Valid > 0) then 
                 Lower_Cloud_Pressure(Elem_Idx,Line_Idx) =  &
                              sum(Cloud_Pressure(Elem_Idx_1:Elem_Idx_2,Line_Idx_1:Line_Idx_2)* & 
                                  mask(Elem_Idx_1:Elem_Idx_2,Line_Idx_1:Line_Idx_2)) / count_Valid
             endif
         
             end do Element_Loop_2
          end do Line_loop_2

          !--- fill missing values with default value
          Line_loop_3: do Line_Idx = Line_start, Line_end 
            Element_Loop_3: do Elem_Idx = 1, Num_Elem 
             if ((Lower_Cloud_Pressure(Elem_Idx,Line_Idx) == Missing_Value_Real4) .and. &
               (Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%NO)) then
               Inwp = Elem_Idx_Nwp(Elem_Idx,Line_Idx)
               Jnwp = Line_Idx_Nwp(Elem_Idx,Line_Idx)
               if (Inwp == Missing_Value_Int4 .or. &
                   Jnwp == Missing_Value_Int4) THEN
                    cycle
               endif 
               Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Surface_Pressure(Elem_Idx,Line_Idx) - Pc_Lower_Cloud_offset
             endif
            end do Element_Loop_3
           end do Line_loop_3

      deallocate(mask)

      endif   !end of Interp_Flag check

      !----------------------------------------------------------------
      !  Compute Height and Temperature
      !----------------------------------------------------------------
      Line_loop_4: do Line_Idx = Line_start, Line_end 
         Element_Loop_4: do Elem_Idx = 1, Num_Elem 

            !-- if a bad pixel, set to missing
            if (Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) then
               Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Missing_Value_Real4
               Lower_Cloud_Temperature(Elem_Idx,Line_Idx) = Missing_Value_Real4 
               Lower_Cloud_Height(Elem_Idx,Line_Idx) = Missing_Value_Real4 
               cycle
            endif

            !--- if not overlap, set to all missing 
            if (Cloud_Type_Local(Elem_Idx,Line_Idx) /= symbol%OVERLAP_TYPE .or. &
                Lower_Cloud_Pressure(Elem_Idx,Line_Idx) == Missing_Value_Real4) then
                Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Missing_Value_Real4
                Lower_Cloud_Height(Elem_Idx,Line_Idx) = Missing_Value_Real4
                Lower_Cloud_Temperature(Elem_Idx,Line_Idx) = Missing_Value_Real4
                cycle
            endif

            !--- compute T and Z from P
            Inwp = Elem_Idx_Nwp(Elem_Idx,Line_Idx)
            Jnwp = Line_Idx_Nwp(Elem_Idx,Line_Idx)
            call KNOWING_P_COMPUTE_T_Z( symbol, &
                                        Lower_Cloud_Pressure(Elem_Idx,Line_Idx), &
                                        Lower_Cloud_Temperature(Elem_Idx,Line_Idx), &
                                        Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                                        Ilev)
        end do Element_Loop_4
       end do Line_loop_4

   end subroutine SPATIALLY_INTERPOLATE_LOWER_CLOUD_POSITION

   !-----------------------------------------------------------------
   ! Interpolate within profiles knowing P to determine T and Z
   !-----------------------------------------------------------------
   subroutine KNOWING_P_COMPUTE_T_Z(symbol,P,T,Z,Ilev)

     type(symbol_acha), intent(in) :: symbol
     real, intent(in):: P
     real, intent(out):: T
     real, intent(out):: Z
     integer, intent(out):: Ilev
     real:: dp
     real:: dt
     real:: dz

     !--- interpolate pressure profile
     call LOCATE(Press_Prof_Rtm,Num_Levels_Rtm_Prof,P,Ilev)
     Ilev = max(1,min(Num_Levels_Rtm_Prof-1,Ilev))

     dp = Press_Prof_Rtm(Ilev+1) - Press_Prof_Rtm(Ilev)
     dt = Temp_Prof_Rtm(Ilev+1) - Temp_Prof_Rtm(Ilev)
     dz = Hght_Prof_Rtm(Ilev+1) - Hght_Prof_Rtm(Ilev)

     !--- perform interpolation
       if (dp /= 0.0) then
           T = Temp_Prof_Rtm(Ilev) + dt/dp * (P - Press_Prof_Rtm(Ilev))
           Z = Hght_Prof_Rtm(Ilev) + dz/dp * (P - Press_Prof_Rtm(Ilev))
       else
           T = Temp_Prof_Rtm(Ilev)
           Z = Hght_Prof_Rtm(Ilev)
       endif
   end subroutine KNOWING_P_COMPUTE_T_Z

   !-----------------------------------------------------------------
   ! Interpolate within profiles knowing Z to determine T and P
   !-----------------------------------------------------------------
   subroutine KNOWING_Z_COMPUTE_T_P(symbol,P,T,Z,Ilev)

     type(symbol_acha), intent(in) :: symbol
     real, intent(in):: Z
     real, intent(out):: T
     real, intent(out):: P
     integer, intent(out):: Ilev
     real:: dp
     real:: dt
     real:: dz

     !--- interpolate pressure profile
     call LOCATE(Hght_Prof_Rtm,Num_Levels_Rtm_Prof,Z,Ilev)
     Ilev = max(1,min(Num_Levels_Rtm_Prof-1,Ilev))

     dp = Press_Prof_Rtm(Ilev+1) - Press_Prof_Rtm(Ilev)
     dt = Temp_Prof_Rtm(Ilev+1) - Temp_Prof_Rtm(Ilev)
     dz = Hght_Prof_Rtm(Ilev+1) - Hght_Prof_Rtm(Ilev)

     !--- perform interpolation
     if (dz /= 0.0) then
           T = Temp_Prof_Rtm(Ilev) + dt/dz * (Z - Hght_Prof_Rtm(Ilev))
           P = Press_Prof_Rtm(Ilev) + dp/dz * (Z - Hght_Prof_Rtm(Ilev))
     else
           T = Temp_Prof_Rtm(Ilev)
           P = Press_Prof_Rtm(Ilev)
     endif

   end subroutine KNOWING_Z_COMPUTE_T_P

   !-----------------------------------------------------------------
   ! Interpolate within profiles knowing T to determine P and Z
   !-----------------------------------------------------------------
   subroutine KNOWING_T_COMPUTE_P_Z(symbol,Cloud_Type,P,T,Z,klev,ierr,Level_Within_Inversion_Flag)

     type(symbol_acha), intent(in) :: symbol
     integer, intent(in):: Cloud_Type
     real, intent(in):: T
     real, intent(out):: P
     real, intent(out):: Z
     integer, intent(out):: klev
     integer, intent(out):: ierr
     real:: dp
     real:: dt
     real:: dz
     integer:: kstart
     integer:: kend
     integer:: nlevels_temp
     integer, intent(out):: Level_Within_Inversion_Flag
     real:: R4_Dummy

     !--- initialization
     ierr = symbol%NO
     Z = Missing_Value_Real4
     P = Missing_Value_Real4

     !--- test for existence of a valid solution with troposphere
     kstart = Tropo_Level_Rtm
     kend = Sfc_Level_Rtm
     Nlevels_Temp = kend - kstart + 1

     !--- check to see if warmer than max, than assume at surface
     if (T > maxval(Temp_Prof_Rtm(kstart:kend))) then
         P = Press_Prof_Rtm(kend)
         Z = Hght_Prof_Rtm(kend)
         klev = kend - 1
         ierr = symbol%NO
         return
     endif

     !--- check to see if colder than min, than assume above tropopause
     !--- and either limit height to tropopause or extrapolate in stratosphere
     if ((T < minval(Temp_Prof_Rtm(kstart:kend))) .or. (klev < Tropo_Level_Rtm)) then
         if (ALLOW_STRATOSPHERE_SOLUTION_FLAG == 1 .and. Cloud_Type == symbol%OVERSHOOTING_TYPE) then
!--->      if (ALLOW_STRATOSPHERE_SOLUTION_FLAG == 1) then
           Z = Hght_Prof_Rtm(kstart) + (T - Temp_Prof_Rtm(kstart)) / Dt_Dz_Strato
           call KNOWING_Z_COMPUTE_T_P(symbol,P,R4_Dummy,Z,klev)
         else
           P = Press_Prof_Rtm(kstart)
           Z = Hght_Prof_Rtm(kstart)
           klev = kstart + 1
         endif
         ierr = symbol%NO
         return
     endif

     !--- if there is an inversion, look below first
     Level_Within_Inversion_Flag = 0
     if (Inver_Level_Rtm > 0 .and. Inver_Level_Rtm < Sfc_Level_Rtm) then
         kstart = Inver_Level_Rtm
         kend =  Sfc_Level_Rtm
         nlevels_temp = kend - kstart + 1
         call LOCATE(Temp_Prof_Rtm(kstart:kend),nlevels_temp,T,klev)
         if ((klev > 0) .and. (klev < nlevels_temp -1)) then
              klev = klev + kstart - 1
              Level_Within_Inversion_Flag = 1
         endif
      endif

    !--- if no solution within an inversion
    if (Level_Within_Inversion_Flag == 0) then
        kstart = Tropo_Level_Rtm
        kend = Sfc_Level_Rtm
        nlevels_temp = kend - kstart + 1
        call LOCATE(Temp_Prof_Rtm(kstart:kend),nlevels_temp,T,klev)
        klev = klev + kstart - 1
        klev = max(1,min(Num_Levels_Rtm_Prof-1,klev))
    endif

    !--- General Inversion
    dp = Press_Prof_Rtm(klev+1) - Press_Prof_Rtm(klev)
    dt = Temp_Prof_Rtm(klev+1) - Temp_Prof_Rtm(klev)
    dz = Hght_Prof_Rtm(klev+1) - Hght_Prof_Rtm(klev)

    if (dt /= 0.0) then
        P = Press_Prof_Rtm(klev) + dp/dt*(T-Temp_Prof_Rtm(klev))
        Z = Hght_Prof_Rtm(klev) + dz/dt*(T-Temp_Prof_Rtm(klev))
    else
        P = Press_Prof_Rtm(klev) 
        Z = Hght_Prof_Rtm(klev)
    endif

   end subroutine KNOWING_T_COMPUTE_P_Z

   !-----------------------------------------------------------------
   ! Interpolate within profiles knowing Z to determine above cloud
   ! radiative terms used in forward model
   !-----------------------------------------------------------------
   function GENERIC_PROFILE_INTERPOLATION(X_value,X_Profile,Y_Profile)  &
            result(Y_value) 

     real, intent(in):: X_value 
     real, dimension(:), intent(in):: X_Profile
     real, dimension(:), intent(in):: Y_Profile
     real:: Y_value

     integer:: Ilev
     real:: dx
     integer:: nlevels

     nlevels = size(X_Profile)

     !--- interpolate pressure profile
     call LOCATE(X_Profile,nlevels,X_value,Ilev)
     Ilev = max(1,min(nlevels-1,Ilev))

     dx = X_Profile(Ilev+1) - X_Profile(Ilev)

     !--- perform interpolation
     if (dx /= 0.0) then
        Y_value = Y_Profile(Ilev) +  &
                 (X_value - X_Profile(Ilev))  * &
                 (Y_Profile(Ilev+1) - Y_Profile(Ilev)) / dx
     else
          Y_value = Y_Profile(Ilev)
     endif
   end function GENERIC_PROFILE_INTERPOLATION

!------------------------------------------------------------------------
! subroutine to compute the Iteration in x due to optimal
! estimation
!
! The notation in this routine follows that of Clive Rodgers (1976,2000)
!
! input to this routine:
! Iter_Idx - the number of the current Iteration
! Iter_Idx_Max - the maximum number of Iterations allowed
! nx - the number of x values
! ny - the number of y values
! conv_crit - the convergence criteria
! y - the vector of observations
! f - the vector of observations predicted by the forward model
! x - the vector of retrieved parameters
! x_Ap - the vector of the apriori estimate of the retrieved parameters
! K - the Kernel Matrix
! Sy - the covariance matrix of y and f
! Sa_inv - the inverse of the covariance matrix of x_Ap
! Delta_X_Max - the maximum step allowed for each Delta_X value
!
! output of this routine:
! Sx - the covariance matrix of x 
! Delta_X - the increment in x for the next Iteration
! Converged_Flag - flag indicating if convergence was met (yes or no)
! Fail_Flag - flag indicating if this process failed (yes or no)
!
! local variables:
! Sx_inv - the inverse of Sx
! Delta_X_dir - the unit direction vectors for delta-x 
! Delta_X_distance - the total length in x-space of Delta_X
! Delta_X_constrained - the values of Delta_X after being constrained
!-----------------------------------------------------------------------
subroutine OPTIMAL_ESTIMATION_ACHA(symbol,Iter_Idx,Iter_Idx_Max,nx,ny, &
                              conv_crit,Delta_X_Max, &
                              y,f,x,x_Ap,K,Sy,Sa_inv, &
                              Sx,AKM,Delta_X,Converged_Flag,Fail_Flag, &
                              idiag_output)

  type(symbol_acha), intent(in) :: symbol
  integer, intent(in):: Iter_Idx
  integer, intent(in):: Iter_Idx_Max
  integer, intent(in):: idiag_output
  integer, intent(in):: ny
  integer, intent(in):: nx
  real(kind=real4), intent(in):: conv_crit
  real(kind=real4), dimension(:), intent(in):: Delta_X_Max
  real(kind=real4), dimension(:), intent(in):: y
  real(kind=real4), dimension(:), intent(in):: f
  real(kind=real4), dimension(:), intent(in):: x
  real(kind=real4), dimension(:), intent(in):: x_Ap
  real(kind=real4), dimension(:,:), intent(in):: K
  real(kind=real4), dimension(:,:), intent(in):: Sy
  real(kind=real4), dimension(:,:), intent(in):: Sa_inv
  real(kind=real4), dimension(:,:), intent(out):: Sx
  real(kind=real4), dimension(:,:), intent(out):: AKM
  real(kind=real4), dimension(:), intent(out):: Delta_X
  real(kind=real4), dimension(ny,ny):: Sy_inv
  real(kind=real4), dimension(nx,nx):: Sx_inv
  real(kind=real4), dimension(nx):: Delta_X_dir
  real(kind=real4), dimension(nx):: Delta_X_constrained
  real(kind=real4):: Delta_X_distance_constrained
  real(kind=real4):: Delta_X_distance
  integer, intent(out):: Fail_Flag
  integer, intent(out):: Converged_Flag
  integer:: Singular_Flag
  real:: Conv_Test
  integer:: ix
  integer:: m
  integer:: p

  m = size(Sy,1)
  p = size(Sx,1)

  Converged_Flag = symbol%NO
  Delta_X = Missing_Value_Real4
  Sx = Missing_Value_Real4

  Singular_Flag =  INVERT_MATRIX(Sy, Sy_Inv, m)
  if (Singular_Flag == symbol%YES) then
   print *, "Cloud Height warning ==> Singular Sy in ACHA "
   Fail_Flag = symbol%YES
   print *, "y = ", y
   print *, "Sy = ", Sy
  !stop
   return 
  endif

  !---- compute next step
  AKM = matmul(Transpose(K),matmul(Sy_inv,K))   !step saving
  Sx_inv = Sa_inv + AKM !(Eq.102 Rodgers)
  Singular_Flag =  INVERT_MATRIX(Sx_inv, Sx, p)
  if (Singular_Flag == symbol%YES) then
   print *, "Cloud Height warning ==> Singular Sx in ACHA "
   print *, "y = ", y
   print *, "f = ", f
   print *, "x_ap = ", x_ap
   print *, "x = ", x
   print *, "Sa_Inv = ", Sa_Inv
   print *, "Sy_Inv = ", Sy_Inv
   print *, "K = ", K
   Fail_Flag = symbol%YES
   return
  endif
  
  Delta_X = matmul(Sx,(matmul(Transpose(K),matmul(Sy_inv,(y-f))) +  &
                       matmul(Sa_inv,x_Ap-x) ))

  !--------------------------------------------------------------
  ! compute averaging kernel matrix (note partialy computed above)
  !--------------------------------------------------------------
  AKM = matmul(Sx,AKM) 

  !--------------------------------------------------------------
  ! check for convergence
  !--------------------------------------------------------------

  !--- compute convergence metric
  Conv_Test = abs(sum(Delta_X*matmul(Sx_inv,Delta_X)))

  !--- control step size  (note change to preserve direction)
  if (idiag_output == symbol%YES) then
          print *, "Sa_inv = ", Sa_inv
          print *, "Sy_inv = ", Sy_inv
          print *, "Sx_inv = ", Sx_inv
          print *, "Sx = ", Sx
          print *, "Delta_X = ", Delta_X
  endif
  Delta_X_distance = sqrt(sum(Delta_X**2))
  if (Delta_X_distance > 0.0) then
     DO ix = 1,nx
        Delta_X_dir(ix) = Delta_X(ix) / Delta_X_distance
     ENDDO
     DO ix = 1,nx
        Delta_X_constrained(ix) =  &
             sign(min(Delta_X_Max(ix),abs(Delta_X(ix))) , Delta_X(ix) )
     end DO
     Delta_X_distance_constrained = sqrt(sum(Delta_X_constrained**2))
     DO ix = 1,nx
        Delta_X(ix) = Delta_X_dir(ix)*Delta_X_distance_constrained
     end DO
  endif

  !--- check for non-traditional convergence
  if ((abs(Delta_X(1)) < 0.1) .and. (Iter_Idx > 1)) then
    Converged_Flag = symbol%YES
    Fail_Flag = symbol%NO
  endif

  !--- check for traditional convergence
  if (Conv_Test < conv_crit) then
      Converged_Flag = symbol%YES
      Fail_Flag = symbol%NO
  endif

  !--- check for exceeding allowed number of interactions
  if (Iter_Idx > Iter_Idx_Max) then
      Converged_Flag = symbol%NO
      Fail_Flag = symbol%YES
  endif

  end subroutine OPTIMAL_ESTIMATION_ACHA

  !----------------------------------------------------------------------
  !---
  !----------------------------------------------------------------------
  subroutine COMPUTE_Sy(symbol, &
                        Acha_Mode_Flag, &
                        Ec,              &
                        T11um_Cal_Uncer, &
                        T11um_12um_Cal_Uncer, &
                        T11um_133um_Cal_Uncer, &
                        T11um_85um_Cal_Uncer, &
                        T11um_67um_Cal_Uncer, &
                        T11um_Clr_Uncer, &
                        T11um_12um_Clr_Uncer, &
                        T11um_133um_Clr_Uncer, &
                        T11um_85um_Clr_Uncer, &
                        T11um_67um_Clr_Uncer, &
                        Bt_Ch31_Std, &
                        Btd_Ch31_Ch32_Std, &
                        Btd_Ch31_Ch33_Std, &
                        Btd_Ch31_Ch29_Std, &
                        Btd_Ch31_Ch27_Std, &
                        Sy)

  type(symbol_acha), intent(in) :: symbol
  integer(kind=int4),intent(in)::Acha_Mode_Flag
  real(kind=real4),intent(in)::Ec
  real(kind=real4),intent(in)::T11um_Cal_Uncer
  real(kind=real4),intent(in)::T11um_12um_Cal_Uncer
  real(kind=real4),intent(in)::T11um_133um_Cal_Uncer
  real(kind=real4),intent(in)::T11um_85um_Cal_Uncer
  real(kind=real4),intent(in)::T11um_67um_Cal_Uncer
  real(kind=real4),intent(in)::T11um_Clr_Uncer
  real(kind=real4),intent(in)::T11um_12um_Clr_Uncer
  real(kind=real4),intent(in)::T11um_133um_Clr_Uncer
  real(kind=real4),intent(in)::T11um_85um_Clr_Uncer
  real(kind=real4),intent(in)::T11um_67um_Clr_Uncer
  real(kind=real4),intent(in)::Bt_Ch31_Std
  real(kind=real4),intent(in)::Btd_Ch31_Ch32_Std
  real(kind=real4),intent(in)::Btd_Ch31_Ch33_Std
  real(kind=real4),intent(in)::Btd_Ch31_Ch29_Std
  real(kind=real4),intent(in)::Btd_Ch31_Ch27_Std
  real(kind=real4),dimension(:,:),intent(out)::Sy

  real(kind=real4):: T11um_Cld_Uncer
  real(kind=real4):: T11um_12um_Cld_Uncer
  real(kind=real4):: T11um_133um_Cld_Uncer
  real(kind=real4):: T11um_85um_Cld_Uncer
  real(kind=real4):: T11um_67um_Cld_Uncer

  !--- compute Sy 
  T11um_Cld_Uncer = max(1.0, Bt_Ch31_Std)
  T11um_12um_Cld_Uncer = max(0.5, Btd_Ch31_Ch32_Std)
  T11um_133um_Cld_Uncer = max(0.5, Btd_Ch31_Ch33_Std)
  T11um_85um_Cld_Uncer = max(0.5, Btd_Ch31_Ch29_Std)
  T11um_67um_Cld_Uncer = max(0.5, Btd_Ch31_Ch27_Std)

  Sy = 0.0
  Sy(1,1) = T11um_Cal_Uncer**2 + ((1.0-Ec)*T11um_Clr_Uncer)**2 + (T11um_Cld_Uncer)**2

  if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 4) then
    Sy(2,2) = T11um_12um_Cal_Uncer**2 + ((1.0-Ec)*T11um_12um_Clr_Uncer)**2 + (T11um_12um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 2) then
    Sy(2,2) = T11um_133um_Cal_Uncer**2 + ((1.0-Ec)*T11um_133um_Clr_Uncer)**2 + (T11um_133um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 6) then
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + ((1.0-Ec)*T11um_67um_Clr_Uncer)**2 + (T11um_67um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 6) then
    Sy(3,3) = T11um_133um_Cal_Uncer**2 + ((1.0-Ec)*T11um_133um_Clr_Uncer)**2 + (T11um_133um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 4) then
    Sy(3,3) = T11um_85um_Cal_Uncer**2 + ((1.0-Ec)*T11um_85um_Clr_Uncer)**2 + (T11um_85um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 7) then
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + ((1.0-Ec)*T11um_67um_Clr_Uncer)**2 + (T11um_67um_Cld_Uncer)**2
  endif


 end subroutine COMPUTE_Sy

 !---------------------------------------------------------------------
 !--- Compute the Forward Model Estimate (f) and its Kernel (df/dx)
 !---------------------------------------------------------------------
 subroutine COMPUTE_FORWARD_MODEL_AND_KERNEL( &
           symbol,Acha_Mode_Flag,  &
           Chan_On_67um, Chan_On_85um, Chan_On_11um, Chan_On_12um, Chan_On_133um, &
           Chan_Idx_67um, Chan_Idx_85um, Chan_Idx_11um, Chan_Idx_12um, Chan_Idx_133um, &
           x,    &
           Rad_Clear_67um, Rad_Ac_67um, Trans_Ac_67um,    &
           Rad_Clear_85um, Rad_Ac_85um, Trans_Ac_85um,    &
           Rad_Clear_11um, Rad_Ac_11um, Trans_Ac_11um,    &
           Rad_Clear_12um, Rad_Ac_12um, Trans_Ac_12um,    &
           Rad_Clear_133um, Rad_Ac_133um, Trans_Ac_133um, &
           a_Beta_11um_133um_fit, b_Beta_11um_133um_fit, &
           a_Beta_11um_85um_fit, b_Beta_11um_85um_fit,   &
           a_Beta_11um_67um_fit, b_Beta_11um_67um_fit,   &
           f,                                            & 
           K, &
           Emiss_Vector, &
           idiag_output)

  type(symbol_acha), intent(in) :: symbol
  integer(kind=int4), intent(in):: Acha_Mode_Flag
  integer, intent(in):: Chan_On_67um
  integer, intent(in):: Chan_On_85um
  integer, intent(in):: Chan_On_11um
  integer, intent(in):: Chan_On_12um
  integer, intent(in):: Chan_On_133um
  integer, intent(in):: Chan_Idx_67um
  integer, intent(in):: Chan_Idx_85um
  integer, intent(in):: Chan_Idx_11um
  integer, intent(in):: Chan_Idx_12um
  integer, intent(in):: Chan_Idx_133um
  real(kind=real4), dimension(:), intent(in):: x
  real(kind=real4), intent(in):: Rad_Clear_67um
  real(kind=real4), intent(in):: Rad_Ac_67um
  real(kind=real4), intent(in):: Trans_Ac_67um
  real(kind=real4), intent(in):: Rad_Clear_85um
  real(kind=real4), intent(in):: Rad_Ac_85um
  real(kind=real4), intent(in):: Trans_Ac_85um
  real(kind=real4), intent(in):: Rad_Clear_11um
  real(kind=real4), intent(in):: Rad_Ac_11um
  real(kind=real4), intent(in):: Trans_Ac_11um
  real(kind=real4), intent(in):: Rad_Clear_12um
  real(kind=real4), intent(in):: Rad_Ac_12um
  real(kind=real4), intent(in):: Trans_Ac_12um
  real(kind=real4), intent(in):: Rad_Clear_133um
  real(kind=real4), intent(in):: Rad_Ac_133um
  real(kind=real4), intent(in):: Trans_Ac_133um
  real(kind=real4), intent(in):: a_Beta_11um_133um_fit
  real(kind=real4), intent(in):: b_Beta_11um_133um_fit
  real(kind=real4), intent(in):: a_Beta_11um_85um_fit
  real(kind=real4), intent(in):: b_Beta_11um_85um_fit
  real(kind=real4), intent(in):: a_Beta_11um_67um_fit
  real(kind=real4), intent(in):: b_Beta_11um_67um_fit
  integer(kind=int4), intent(in):: idiag_output
  real(kind=real4), dimension(:), intent(out):: f 
  real(kind=real4), dimension(:,:), intent(out):: K
  real(kind=real4), dimension(:), intent(out):: Emiss_Vector

  real(kind=real4):: Tc
  real(kind=real4):: Bc_67um
  real(kind=real4):: Bc_85um
  real(kind=real4):: Bc_11um
  real(kind=real4):: Bc_12um
  real(kind=real4):: Bc_133um
  real(kind=real4):: Rad_67um
  real(kind=real4):: Rad_85um
  real(kind=real4):: Rad_11um
  real(kind=real4):: Rad_12um
  real(kind=real4):: Rad_133um
  real(kind=real4):: Emiss_67um
  real(kind=real4):: Emiss_85um
  real(kind=real4):: Emiss_11um
  real(kind=real4):: Emiss_12um
  real(kind=real4):: Emiss_133um
  real(kind=real4):: Trans_67um
  real(kind=real4):: Trans_85um
  real(kind=real4):: Trans_11um
  real(kind=real4):: Trans_12um
  real(kind=real4):: Trans_133um
  real(kind=real4):: Beta_11um_67um
  real(kind=real4):: Beta_11um_85um
  real(kind=real4):: Beta_11um_12um
  real(kind=real4):: Beta_11um_133um
  real(kind=real4):: dEmiss_67um_dEmiss_11um
  real(kind=real4):: dEmiss_85um_dEmiss_11um
  real(kind=real4):: dEmiss_12um_dEmiss_11um
  real(kind=real4):: dEmiss_133um_dEmiss_11um
  real(kind=real4):: dBeta_11um_67um_dBeta_11um_12um
  real(kind=real4):: dBeta_11um_85um_dBeta_11um_12um
  real(kind=real4):: dBeta_11um_133um_dBeta_11um_12um
  real(kind=real4):: dB_dT_67um
  real(kind=real4):: dB_dT_85um
  real(kind=real4):: dB_dT_11um
  real(kind=real4):: dB_dT_12um
  real(kind=real4):: dB_dT_133um
  real(kind=real4):: dB_dTc_67um
  real(kind=real4):: dB_dTc_85um
  real(kind=real4):: dB_dTc_11um
  real(kind=real4):: dB_dTc_12um
  real(kind=real4):: dB_dTc_133um


  !---  for notational convenience, rename elements of x to local variables
  Tc = x(1)
  Emiss_11um = min(x(2),0.999999)    !values must be below unity
  Beta_11um_12um = x(3)

  !--- compute planck Emission for cloud temperature
  if (Chan_On_67um == symbol%YES) Bc_67um = PLANCK_RAD_FAST( Chan_Idx_67um, Tc, dB_dT = dB_dTc_67um)
  if (Chan_On_85um == symbol%YES) Bc_85um = PLANCK_RAD_FAST( Chan_Idx_85um, Tc, dB_dT = dB_dTc_85um)
  if (Chan_On_11um == symbol%YES) Bc_11um = PLANCK_RAD_FAST( Chan_Idx_11um, Tc, dB_dT = dB_dTc_11um)
  if (Chan_On_12um == symbol%YES) Bc_12um = PLANCK_RAD_FAST( Chan_Idx_12um, Tc, dB_dT = dB_dTc_12um)
  if (Chan_On_133um == symbol%YES) Bc_133um = PLANCK_RAD_FAST( Chan_Idx_133um, Tc, dB_dT = dB_dTc_133um)

  if (idiag_output == symbol%YES) then
          print *, "Bc = ", Bc_67um, Bc_85um, Bc_11um, Bc_12um, Bc_133um
  endif

  !----- compute channel Emissivities

  !-- ch32
  dEmiss_12um_dEmiss_11um = Beta_11um_12um * (1.0-Emiss_11um)**(Beta_11um_12um - 1.0)
  Emiss_12um = 1.0 - (1.0-Emiss_11um)**Beta_11um_12um

  !--ch33
  Beta_11um_133um = a_Beta_11um_133um_fit + b_Beta_11um_133um_fit * Beta_11um_12um
  dBeta_11um_133um_dBeta_11um_12um = b_Beta_11um_133um_fit
  dEmiss_133um_dEmiss_11um = Beta_11um_133um * (1.0-Emiss_11um)**(Beta_11um_133um - 1.0)
  Emiss_133um = 1.0 - (1.0-Emiss_11um)**Beta_11um_133um

  !--ch29
  Beta_11um_85um = a_Beta_11um_85um_fit + b_Beta_11um_85um_fit * Beta_11um_12um
  dBeta_11um_85um_dBeta_11um_12um = b_Beta_11um_85um_fit
  dEmiss_85um_dEmiss_11um = Beta_11um_85um * (1.0-Emiss_11um)**(Beta_11um_85um - 1.0)
  Emiss_85um = 1.0 - (1.0-Emiss_11um)**Beta_11um_85um

  !--ch27
  Beta_11um_67um = a_Beta_11um_67um_fit + b_Beta_11um_67um_fit * Beta_11um_12um
  dBeta_11um_67um_dBeta_11um_12um = b_Beta_11um_67um_fit
  dEmiss_67um_dEmiss_11um = Beta_11um_67um * (1.0-Emiss_11um)**(Beta_11um_67um - 1.0)
  Emiss_67um = 1.0 - (1.0-Emiss_11um)**Beta_11um_67um

  if (idiag_output == symbol%YES) then
          print *, "Ec = ", Emiss_67um,Emiss_85um, Emiss_11um, Emiss_12um, Emiss_133um
  endif

 !--- define Transmission as complement of Emissivity
 Trans_67um = (1.0 - Emiss_67um)
 Trans_85um = (1.0 - Emiss_85um)
 Trans_11um = (1.0 - Emiss_11um)
 Trans_12um = (1.0 - Emiss_12um)
 Trans_133um = (1.0 - Emiss_133um)

 !--- forward model and Kernel for  11um channel
 Rad_11um = Emiss_11um*Rad_Ac_11um + Trans_Ac_11um * Emiss_11um * Bc_11um +  Trans_11um * Rad_Clear_11um
 f(1) = PLANCK_TEMP_FAST(Chan_Idx_11um,Rad_11um,dB_dT = dB_dT_11um)
 K(1,1) = (Trans_Ac_11um * Emiss_11um * dB_dTc_11um) / dB_dT_11um               !dT_11um / dT_c
 K(1,2) = (Rad_Ac_11um + Trans_Ac_11um*Bc_11um - Rad_Clear_11um)/dB_dT_11um     !dT_11um / dEmiss_c
 K(1,3) = 0.0                                                                   !dT_11um / dbeta

 !--- forward model for 11um - 12um
 if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 5) then
   Rad_12um = Emiss_12um*Rad_Ac_12um + Trans_Ac_12um * Emiss_12um * Bc_12um +  Trans_12um * Rad_Clear_12um
   f(2) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_12um,Rad_12um,dB_dT = dB_dT_12um)
   K(2,1) = K(1,1) - Trans_Ac_12um * Emiss_12um * dB_dTc_12um / dB_dT_12um   
   K(2,2) = K(1,2) -  (Rad_Ac_12um + Trans_Ac_12um*Bc_12um-Rad_Clear_12um)*&
                      (dEmiss_12um_dEmiss_11um)/dB_dT_12um  
   K(2,3) = (Rad_Ac_12um+Trans_Ac_12um*Bc_12um-Rad_Clear_12um)/ &
            dB_dT_12um*alog(1.0-Emiss_11um)*(1.0-Emiss_12um)    
 endif
 !--- forward model for 11um - 133um
 if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 6) then
   Rad_133um = Emiss_133um*Rad_Ac_133um + Trans_Ac_133um * Emiss_133um * Bc_133um +  Trans_133um * Rad_Clear_133um
   f(2) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_133um,Rad_133um,dB_dT = dB_dT_133um)
   K(2,1) = K(1,1) - Trans_Ac_133um*Emiss_133um*dB_dTc_133um/dB_dT_133um
   K(2,2) = K(1,2) - (Rad_Ac_133um + Trans_ac_133um*Bc_133um - Rad_Clear_133um)*&
                     (dEmiss_133um_dEmiss_11um)/dB_dT_133um
   K(2,3) = (Rad_Ac_133um + Trans_ac_133um*Bc_133um -Rad_Clear_133um)/ &
             dB_dT_133um * alog(1.0-Emiss_11um)*(1.0-Emiss_133um)
 endif
 if (Acha_Mode_Flag == 3) then
   Rad_133um = Emiss_133um*Rad_Ac_133um + Trans_Ac_133um * Emiss_133um * Bc_133um +  Trans_133um * Rad_Clear_133um
   f(3) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_133um,Rad_133um,dB_dT = dB_dT_133um)
   K(3,1) = K(1,1) - Trans_Ac_133um*Emiss_133um*dB_dTc_133um/dB_dT_133um
   K(3,2) = K(1,2) - (Rad_Ac_133um + Trans_ac_133um*Bc_133um - Rad_Clear_133um)*&
                     (dEmiss_133um_dEmiss_11um)/dB_dT_133um
   K(3,3) = (Rad_Ac_133um + Trans_ac_133um*Bc_133um -Rad_Clear_133um)/ &
             dB_dT_133um * alog(1.0-Emiss_11um)*(1.0-Emiss_133um)
 endif
 if (Acha_Mode_Flag == 4) then
   Rad_85um = Emiss_85um*Rad_Ac_85um + Trans_Ac_85um * Emiss_85um * Bc_85um +  Trans_85um * Rad_Clear_85um
   f(3) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_85um,Rad_85um,dB_dT = dB_dT_85um)
   K(3,1) = K(1,1) - Trans_Ac_85um*Emiss_85um*dB_dTc_85um/dB_dT_85um
   K(3,2) = K(1,2) - (Rad_Ac_85um + Trans_ac_85um*Bc_85um - Rad_Clear_85um)*&
                     (dEmiss_85um_dEmiss_11um)/dB_dT_85um
   K(3,3) = (Rad_Ac_85um + Trans_ac_85um*Bc_85um -Rad_Clear_85um)/ &
             dB_dT_85um * alog(1.0-Emiss_11um)*(1.0-Emiss_85um)
 endif
 if (Acha_Mode_Flag ==5 .or. Acha_Mode_Flag == 6) then
   Rad_67um = Emiss_67um*Rad_Ac_67um + Trans_Ac_67um * Emiss_67um * Bc_67um + Trans_67um * Rad_Clear_67um
   f(3) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_67um,Rad_67um,dB_dT = dB_dT_67um)
   K(3,1) = K(1,1) - Trans_Ac_67um*Emiss_67um*dB_dTc_67um/dB_dT_67um
   K(3,2) = K(1,2) - (Rad_Ac_67um + Trans_ac_67um*Bc_67um - Rad_Clear_67um)*&
                     (dEmiss_67um_dEmiss_11um)/dB_dT_67um
   K(3,3) = (Rad_Ac_67um + Trans_ac_67um*Bc_67um - Rad_Clear_67um)/ &
             dB_dT_67um * alog(1.0-Emiss_11um)*(1.0-Emiss_67um)
 endif
 if (Acha_Mode_Flag == 7) then
   Rad_67um = Emiss_67um*Rad_Ac_67um + Trans_Ac_67um * Emiss_67um * Bc_67um + Trans_67um * Rad_Clear_67um
   f(2) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_67um,Rad_67um,dB_dT = dB_dT_67um)
   K(2,1) = K(1,1) - Trans_Ac_67um*Emiss_67um*dB_dTc_67um/dB_dT_67um
   K(2,2) = K(1,2) - (Rad_Ac_67um + Trans_ac_67um*Bc_67um - Rad_Clear_67um)*&
                     (dEmiss_67um_dEmiss_11um)/dB_dT_67um
   K(2,3) = (Rad_Ac_67um + Trans_ac_67um*Bc_67um - Rad_Clear_67um)/ &
             dB_dT_67um * alog(1.0-Emiss_11um)*(1.0-Emiss_67um)
 endif

 !--- determine number of channels
  select case(Acha_Mode_Flag)
     case(0)  !avhrr, goes-im
       Emiss_Vector(1) = Emiss_11um
     case(1)  !avhrr, goes-im
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_12um
     case(2)  !goes-nop
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_133um
     case(3)  !goes-r
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_12um
       Emiss_Vector(3) = Emiss_133um
     case(4)  !viirs
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_12um
       Emiss_Vector(3) = Emiss_85um
     case(5)  !goes-im 3 chan
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_67um
       Emiss_Vector(3) = Emiss_12um
     case(6)  !goes-np 3 chan
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_67um
       Emiss_Vector(3) = Emiss_133um
     case(7)  !goes-np 3 chan
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_67um
  end select

 if (idiag_output == symbol%YES) then
         print *, "f = ", f
 endif

end subroutine COMPUTE_FORWARD_MODEL_AND_KERNEL

!----------------------------------------------------------------------
subroutine COMPUTE_APRIORI_BASED_ON_TYPE(symbol,Cloud_Type, &
                           Ttropo,T11um,mu,Tc_Opaque, &
                           Tc_Ap, Tc_Ap_Uncer, &
                           Ec_Ap, Ec_Ap_Uncer, &
                           Beta_Ap, Beta_Ap_Uncer)

  type(symbol_acha), intent(in) :: symbol
  integer, intent(in):: Cloud_Type
  real(kind=real4), intent(in):: Ttropo
  real(kind=real4), intent(in):: T11um
  real(kind=real4), intent(in):: mu
  real(kind=real4), intent(in):: Tc_Opaque
  real(kind=real4), intent(out):: Tc_Ap
  real(kind=real4), intent(out):: Ec_Ap
  real(kind=real4), intent(out):: Tc_Ap_Uncer
  real(kind=real4), intent(out):: Ec_Ap_Uncer
  real(kind=real4), intent(out):: Beta_Ap
  real(kind=real4), intent(out):: Beta_Ap_Uncer
  real(kind=real4):: Tc_Ap_Cirrus
  real(kind=real4):: Tc_Ap_Opaque
  real(kind=real4):: Tau_Ap

  !--- calipso values (not multiplier on uncer values)
  Tc_Ap_Cirrus = Ttropo + Tc_Ap_Tropo_OFFSET_Cirrus
  Tc_Ap_Opaque = Tc_Opaque
  if (Tc_Ap_Opaque == Missing_Value_Real4) then
     Tc_Ap_Opaque = T11um
  endif

  if (Cloud_Type == symbol%FOG_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_Fog_Type  
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = 2.0*Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == symbol%WATER_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_Water_Type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == symbol%SUPERCOOLED_TYPE .or. &
          Cloud_Type == symbol%MIXED_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_Mixed_Type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == symbol%TICE_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_Opaque_Ice_Type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == symbol%OVERSHOOTING_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_Opaque_Ice_Type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == symbol%CIRRUS_TYPE) then
     Tc_Ap = Tc_Ap_Cirrus
     Tau_Ap = Tau_Ap_Cirrus_Type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus
     Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus
  elseif (Cloud_Type == symbol%OVERLAP_TYPE) then
     Tc_Ap = Tc_Ap_Cirrus
     Tau_Ap = Tau_Ap_overlap_Type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus
     Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus
  endif

  !--- determine beta apriori and fit parameters based on 
  !--- phase (derived from type)
          
  !--- water phase clouds
  if ((Cloud_Type == symbol%FOG_TYPE) .or. &
     (Cloud_Type == symbol%WATER_TYPE) .or. &
      (Cloud_Type == symbol%SUPERCOOLED_TYPE)) THEN
      Beta_Ap = Beta_Ap_Water
      Beta_Ap_Uncer = Beta_Ap_Uncer_Water
  else
  !--- ice phase clouds
      Beta_Ap = Beta_Ap_Ice
      Beta_Ap_Uncer = Beta_Ap_Uncer_Ice
  endif

  !--- convert Tau_Ap to Emissivity
  Ec_Ap = 1.0 - exp(-Tau_Ap/mu)  !slow!

end subroutine COMPUTE_APRIORI_BASED_ON_TYPE

!----------------------------------------------------------------------
!---
!----------------------------------------------------------------------
subroutine COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO( &
                           symbol, &
                           Cloud_Phase, &
                           Emiss_11um_Tropo, &
                           Ttropo, &
                           T11um, &
                           T11um_Lrc, &
                           Tc_Opaque, &
                           Mu, &
                           Tc_Ap, &
                           Tc_Ap_Uncer, &
                           Ec_Ap, &
                           Ec_Ap_Uncer, &
                           Beta_Ap,  &
                           Beta_Ap_Uncer)

  type(symbol_acha), intent(in) :: symbol
  integer, intent(in):: Cloud_Phase
  real(kind=real4), intent(in):: Emiss_11um_Tropo
  real(kind=real4), intent(in):: Ttropo
  real(kind=real4), intent(in):: T11um
  real(kind=real4), intent(in):: T11um_Lrc
  real(kind=real4), intent(in):: Tc_Opaque
  real(kind=real4), intent(in):: Mu
  real(kind=real4), intent(out):: Tc_Ap
  real(kind=real4), intent(out):: Ec_Ap
  real(kind=real4), intent(out):: Beta_Ap
  real(kind=real4), intent(out):: Tc_Ap_Uncer
  real(kind=real4), intent(out):: Ec_Ap_Uncer
  real(kind=real4), intent(out):: Beta_Ap_Uncer

  real(kind=real4):: Tc_Ap_Cirrus
  real(kind=real4):: Tc_Ap_Opaque
  real(kind=real4):: Emiss_Weight
  real(kind=real4):: Emiss_Weight2

  !--- calipso values (not multiplier on uncer values)
  Tc_Ap_Cirrus = Ttropo + Tc_Ap_Tropo_OFFSET_Cirrus
  Tc_Ap_Opaque = Tc_Opaque
  if (T11um_Lrc /= Missing_Value_Real4) then
      Tc_Ap_Opaque = T11um_Lrc
  endif
  if (Tc_Ap_Opaque == Missing_Value_Real4) then
     Tc_Ap_Opaque = T11um
  endif

  if (Cloud_Phase /= symbol%ICE_PHASE) then
    Tc_Ap = Tc_Ap_Opaque
    Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
    Ec_Ap = 1.0 - exp(-1.0*Tau_Ap_Water_Phase/Mu)  !slow!
    Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
    Beta_Ap = Beta_Ap_Water
    Beta_Ap_Uncer = Beta_Ap_Uncer_Water
  endif

  if (Cloud_Phase == symbol%ICE_PHASE) then

    if (Emiss_11um_Tropo <= 0.0) then
            Emiss_Weight = 0.0
    elseif (Emiss_11um_Tropo > 1.0) then
            Emiss_Weight = 1.0
    else
            Emiss_Weight = Emiss_11um_Tropo
    endif

    Emiss_Weight2 = Emiss_Weight

    Tc_Ap = Emiss_Weight2*Tc_Ap_Opaque + &
            (1.0-Emiss_Weight2)*Tc_Ap_Cirrus

    Tc_Ap_Uncer = Emiss_Weight2*Tc_Ap_Uncer_Opaque + &
            (1.0-Emiss_Weight2)*Tc_Ap_Uncer_Cirrus

    Ec_Ap = min(0.99,max(0.1,Emiss_Weight)) 
    Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus

    Beta_Ap = Beta_Ap_Ice
    Beta_Ap_Uncer = Beta_Ap_Uncer_Ice

  endif

  end subroutine COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO

!-------------------------------------------------------------------
! Determine surface type for use in forward model
! 0 = Water
! 1 = Land
! 2 = Snow
! 3 = Desert
! 4 = Arctic
! 5 = Antarctic
!-------------------------------------------------------------------
subroutine DETERMINE_SFC_TYPE_FORWARD_MODEL(symbol, &
                                         Surface_Type, &
                                         Snow_Class, &
                                         Latitude, &
                                         Ch20_Surface_Emissivity, &
                                         Sfc_Type_Forward_Model)

type(symbol_acha), intent(in) :: symbol
integer(kind=int1), intent(in):: Surface_Type
integer(kind=int1), intent(in):: Snow_Class
real(kind=real4), intent(in):: Latitude
real(kind=real4), intent(in):: Ch20_Surface_Emissivity
integer(kind=int4), intent(out):: Sfc_Type_Forward_Model


if (Surface_Type == symbol%WATER_SFC) then
        Sfc_Type_Forward_Model = 0
else
        Sfc_Type_Forward_Model = 1   !Land
endif
if (Snow_Class == symbol%SNOW .and. &
    Latitude > -60.0) then
        Sfc_Type_Forward_Model = 2   !Snow
endif
if (Surface_Type /= symbol%WATER_SFC .and. &
    Snow_Class == symbol%NO_SNOW .and.  &
    Ch20_Surface_Emissivity > 0.90 .and.  &
    abs(Latitude) < 60.0) then
        Sfc_Type_Forward_Model = 3   !Desert
endif
if (Snow_Class == symbol%SEA_ICE .and. &
    Latitude > -60.0) then
        Sfc_Type_Forward_Model = 4   !Arctic
endif
if (Snow_Class /= symbol%NO_SNOW .and. Latitude < -60.0) then
        Sfc_Type_Forward_Model = 5   !Antartica
endif

end subroutine DETERMINE_SFC_TYPE_FORWARD_MODEL

!----------------------------------------------------------------------
! Compute Sy based on the clear-sky error covariance calculations.
! Using Andy's simpler expression
!
! This assumes that 
! Acha_Mode_Flag: 0=11um,1=11+12um,2=11+13.3um,3=11+12+13.3um,4=8.5+11+12um
!            6=11+6.7+13um
!
! Input:
! Emiss_Vector = a vector of emissivities in each channel. 
! Acha_Mode_Flag: 0=11um,1=11+12um,2=11+13.3um,3=11+12+13.3um,4=8.5+11+12um
! Sfc_Type_Forward_Model = the surface type used for covariance calcs 
! y_variance = the variance computed a 3x3 array for each element of y
!
! Output:
! Sy = error covariance matrix
!
!  Sy(i,i) = Cal_Err_y(i)^2 + Spatial_Variance_y(i) + (1-emiss_i)^2*y(i)_covar
!  Sy(i,j) = (1-emiss_i)*(1-emiss_j)*y(i)_y(x)_covar
!----------------------------------------------------------------------
subroutine COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE(   &
                                             symbol, &
                                             Emiss_Vector,  &
                                             Acha_Mode_Flag, &
                                             y_variance, &
                                             Sy) 


 type(symbol_acha), intent(in) :: symbol
 real(kind=real4), intent(in), dimension(:):: Emiss_Vector
 integer(kind=int4), intent(in):: Acha_Mode_Flag
 real(kind=real4), intent(in), dimension(:):: y_variance
 real(kind=real4), intent(out), dimension(:,:):: Sy
 real(kind=real4), dimension(size(y_variance)):: Sub_Pixel_Uncer

 real(kind=real4):: Emiss_11um
 
 real(kind=real4):: Trans2   

 Emiss_11um = min(1.0,max(0.0,Emiss_Vector(1)))

 Trans2 = (1.0 - Emiss_11um)**2  !cloud transmission squared

 !----------------------------------------------------------------
 !--- modify y_variance to represent a sub-pixel uncertainty 
 !--- assume that half of standard deviation is due to sub-pixel
 !--- heterogeneity and that this is a good estimate of the
 !--- forward model error due to sub-pixel heterogeneity
 !----------------------------------------------------------------
 Sub_Pixel_Uncer(1) = max(0.5,y_variance(1)/4.0)
 if (Acha_Mode_Flag > 0) then
    Sub_Pixel_Uncer(2) = max(0.25,y_variance(2)/4.0)
 endif
 if (Acha_Mode_Flag >= 3 .and. Acha_Mode_Flag /=7) then
    Sub_Pixel_Uncer(3) = max(0.25,y_variance(3)/4.0)
 endif

 if (Acha_Mode_Flag == 0) then
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
 endif
 if (Acha_Mode_Flag == 1) then
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_12um_Covar
    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_12um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_12um_Btd_11um_12um_Covar
 endif
 if (Acha_Mode_Flag == 2) then
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_133um_Covar
    Sy(2,1) = Trans2*Bt_11um_Btd_11um_133um_Covar
    Sy(2,2) = T11um_133um_Cal_Uncer**2 +  Sub_Pixel_Uncer(2)+ Trans2*Btd_11um_133um_Btd_11um_133um_Covar
 endif
 if (Acha_Mode_Flag == 3) then
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_12um_Covar
    Sy(1,3) = Trans2*Bt_11um_Btd_11um_133um_Covar
    Sy(2,2) = T11um_12um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + &
              Trans2*Btd_11um_12um_Btd_11um_12um_Covar
    Sy(2,1) = Sy(1,2)
    Sy(2,3) = Trans2*Btd_11um_12um_Btd_11um_133um_Covar
    Sy(3,3) = T11um_133um_Cal_Uncer**2 + Sub_Pixel_Uncer(3) + Trans2*Btd_11um_133um_Btd_11um_133um_Covar
    Sy(3,1) = Sy(1,3)
    Sy(3,2) = Sy(2,3)
 endif
 if (Acha_Mode_Flag == 4) then
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_12um_Covar
    Sy(1,3) = Trans2*Bt_11um_Btd_11um_85um_Covar

    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_12um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_12um_Btd_11um_12um_Covar
    Sy(2,3) = Trans2*Btd_11um_12um_Btd_11um_85um_Covar

    Sy(3,3) = T11um_85um_Cal_Uncer**2 + Sub_Pixel_Uncer(3) + Trans2*Btd_11um_85um_Btd_11um_85um_Covar
    Sy(3,1) = Sy(1,3)
    Sy(3,2) = Sy(2,3)
 endif

 if (Acha_Mode_Flag == 5) then
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_12um_Covar
    Sy(1,3) = Trans2*Bt_11um_Btd_11um_67um_Covar

    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_67um_Btd_11um_67um_Covar
    Sy(2,3) = Trans2*Btd_11um_12um_Btd_11um_67um_Covar

    Sy(3,3) = T11um_12um_Cal_Uncer**2 + Sub_Pixel_Uncer(3) + Trans2*Btd_11um_12um_Btd_11um_12um_Covar
    Sy(3,1) = Sy(1,3)
    Sy(3,2) = Sy(2,3)
 endif
 if (Acha_Mode_Flag == 6) then
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_133um_Covar
    Sy(1,3) = Trans2*Bt_11um_Btd_11um_67um_Covar

    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_67um_Btd_11um_67um_Covar
    Sy(2,3) = Trans2*Btd_11um_67um_Btd_11um_133um_Covar

    Sy(3,3) = T11um_133um_Cal_Uncer**2 + Sub_Pixel_Uncer(3) + Trans2*Btd_11um_133um_Btd_11um_133um_Covar
    Sy(3,1) = Sy(1,3)
    Sy(3,2) = Sy(2,3)
 endif
 if (Acha_Mode_Flag == 7) then
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_67um_Covar
    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_67um_Btd_11um_67um_Covar
 endif


end subroutine COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE

!----------------------------------------------------------------------
! Compute Sy based on the clear-sky error covariance calculations.
! This assumes that 
! Acha_Mode_Flag: 0=11um,1=11+12um,2=11+13.3um,3=11+12+13.3um,4=8.5+11+12um
!                 5=11+6.7+13.3um
!----------------------------------------------------------------------
subroutine SET_CLEAR_SKY_COVARIANCE_TERMS(Sfc_Type_Forward_Model)

  integer(kind=int4), intent(in):: Sfc_Type_Forward_Model

 !--- Water
 if (Sfc_Type_Forward_Model == 0) then
   Bt_67um_Mean = Bt_67um_Mean_Water
   Bt_85um_Mean = Bt_85um_Mean_Water
   Bt_11um_Mean = Bt_11um_Mean_Water
   Bt_12um_Mean = Bt_12um_Mean_Water
   Bt_133um_Mean = Bt_133um_Mean_Water
   
   Bt_67um_Bt_67um_Covar = Bt_67um_Bt_67um_Covar_Water
   Bt_85um_Bt_85um_Covar = Bt_85um_Bt_85um_Covar_Water
   Bt_11um_Bt_11um_Covar = Bt_11um_Bt_11um_Covar_Water
   Bt_12um_Bt_12um_Covar = Bt_12um_Bt_12um_Covar_Water
   Bt_133um_Bt_133um_Covar = Bt_133um_Bt_133um_Covar_Water

   Bt_85um_Bt_133um_Covar = Bt_85um_Bt_133um_Covar_Water
   Bt_12um_Bt_133um_Covar = Bt_12um_Bt_133um_Covar_Water
   Bt_11um_Bt_67um_Covar = Bt_11um_Bt_67um_Covar_Water
   Bt_11um_Bt_85um_Covar = Bt_11um_Bt_85um_Covar_Water
   Bt_12um_Bt_85um_Covar = Bt_12um_Bt_85um_Covar_Water
   Bt_12um_Bt_67um_Covar = Bt_12um_Bt_67um_Covar_Water
   Bt_11um_Bt_12um_Covar = Bt_11um_Bt_12um_Covar_Water
   Bt_11um_Bt_133um_Covar = Bt_11um_Bt_133um_Covar_Water
   Bt_67um_Bt_133um_Covar = Bt_67um_Bt_133um_Covar_Water

   Bt_11um_Btd_11um_67um_Covar = Bt_11um_Btd_11um_67um_Covar_Water
   Bt_11um_Btd_11um_85um_Covar = Bt_11um_Btd_11um_85um_Covar_Water
   Bt_11um_Btd_11um_12um_Covar = Bt_11um_Btd_11um_12um_Covar_Water
   Bt_11um_Btd_11um_133um_Covar = Bt_11um_Btd_11um_133um_Covar_Water

   Btd_11um_67um_Btd_11um_67um_Covar = Btd_11um_67um_Btd_11um_67um_Covar_Water
   Btd_11um_85um_Btd_11um_85um_Covar = Btd_11um_85um_Btd_11um_85um_Covar_Water
   Btd_11um_12um_Btd_11um_12um_Covar = Btd_11um_12um_Btd_11um_12um_Covar_Water
   Btd_11um_133um_Btd_11um_133um_Covar = Btd_11um_133um_Btd_11um_133um_Covar_Water
   Btd_11um_12um_Btd_11um_133um_Covar = Btd_11um_12um_Btd_11um_133um_Covar_Water
   Btd_11um_12um_Btd_11um_85um_Covar = Btd_11um_12um_Btd_11um_85um_Covar_Water
   Btd_11um_12um_Btd_11um_67um_Covar = Btd_11um_12um_Btd_11um_67um_Covar_Water
   Btd_11um_67um_Btd_11um_133um_Covar = Btd_11um_67um_Btd_11um_133um_Covar_Water
 endif
 !--- Land
 if (Sfc_Type_Forward_Model == 1) then
   Bt_67um_Mean = Bt_67um_Mean_Land
   Bt_85um_Mean = Bt_85um_Mean_Land
   Bt_11um_Mean = Bt_11um_Mean_Land
   Bt_12um_Mean = Bt_12um_Mean_Land
   Bt_133um_Mean = Bt_133um_Mean_Land
   
   Bt_67um_Bt_67um_Covar = Bt_67um_Bt_67um_Covar_Land
   Bt_85um_Bt_85um_Covar = Bt_85um_Bt_85um_Covar_Land
   Bt_11um_Bt_11um_Covar = Bt_11um_Bt_11um_Covar_Land
   Bt_12um_Bt_12um_Covar = Bt_12um_Bt_12um_Covar_Land
   Bt_133um_Bt_133um_Covar = Bt_133um_Bt_133um_Covar_Land

   Bt_85um_Bt_133um_Covar = Bt_85um_Bt_133um_Covar_Land
   Bt_12um_Bt_133um_Covar = Bt_12um_Bt_133um_Covar_Land
   Bt_11um_Bt_67um_Covar = Bt_11um_Bt_67um_Covar_Land
   Bt_11um_Bt_85um_Covar = Bt_11um_Bt_85um_Covar_Land
   Bt_12um_Bt_85um_Covar = Bt_12um_Bt_85um_Covar_Land
   Bt_12um_Bt_67um_Covar = Bt_12um_Bt_67um_Covar_Land
   Bt_11um_Bt_12um_Covar = Bt_11um_Bt_12um_Covar_Land
   Bt_11um_Bt_133um_Covar = Bt_11um_Bt_133um_Covar_Land
   Bt_67um_Bt_133um_Covar = Bt_67um_Bt_133um_Covar_Land

   Bt_11um_Btd_11um_67um_Covar = Bt_11um_Btd_11um_67um_Covar_Land
   Bt_11um_Btd_11um_85um_Covar = Bt_11um_Btd_11um_85um_Covar_Land
   Bt_11um_Btd_11um_12um_Covar = Bt_11um_Btd_11um_12um_Covar_Land
   Bt_11um_Btd_11um_133um_Covar = Bt_11um_Btd_11um_133um_Covar_Land

   Btd_11um_67um_Btd_11um_67um_Covar = Btd_11um_67um_Btd_11um_67um_Covar_Land
   Btd_11um_85um_Btd_11um_85um_Covar = Btd_11um_85um_Btd_11um_85um_Covar_Land
   Btd_11um_12um_Btd_11um_12um_Covar = Btd_11um_12um_Btd_11um_12um_Covar_Land
   Btd_11um_133um_Btd_11um_133um_Covar = Btd_11um_133um_Btd_11um_133um_Covar_Land
   Btd_11um_12um_Btd_11um_133um_Covar = Btd_11um_12um_Btd_11um_133um_Covar_Land
   Btd_11um_12um_Btd_11um_85um_Covar = Btd_11um_12um_Btd_11um_85um_Covar_Land
   Btd_11um_12um_Btd_11um_67um_Covar = Btd_11um_12um_Btd_11um_67um_Covar_Land
   Btd_11um_67um_Btd_11um_133um_Covar = Btd_11um_67um_Btd_11um_133um_Covar_Land
 endif

 !--- Snow
 if (Sfc_Type_Forward_Model == 2) then
   Bt_67um_Mean = Bt_67um_Mean_Snow
   Bt_85um_Mean = Bt_85um_Mean_Snow
   Bt_11um_Mean = Bt_11um_Mean_Snow
   Bt_12um_Mean = Bt_12um_Mean_Snow
   Bt_133um_Mean = Bt_133um_Mean_Snow
   
   Bt_67um_Bt_67um_Covar = Bt_67um_Bt_67um_Covar_Snow
   Bt_85um_Bt_85um_Covar = Bt_85um_Bt_85um_Covar_Snow
   Bt_11um_Bt_11um_Covar = Bt_11um_Bt_11um_Covar_Snow
   Bt_12um_Bt_12um_Covar = Bt_12um_Bt_12um_Covar_Snow
   Bt_133um_Bt_133um_Covar = Bt_133um_Bt_133um_Covar_Snow

   Bt_85um_Bt_133um_Covar = Bt_85um_Bt_133um_Covar_Snow
   Bt_12um_Bt_133um_Covar = Bt_12um_Bt_133um_Covar_Snow
   Bt_11um_Bt_67um_Covar = Bt_11um_Bt_67um_Covar_Snow
   Bt_11um_Bt_85um_Covar = Bt_11um_Bt_85um_Covar_Snow
   Bt_12um_Bt_85um_Covar = Bt_12um_Bt_85um_Covar_Snow
   Bt_12um_Bt_67um_Covar = Bt_12um_Bt_67um_Covar_Snow
   Bt_11um_Bt_12um_Covar = Bt_11um_Bt_12um_Covar_Snow
   Bt_11um_Bt_133um_Covar = Bt_11um_Bt_133um_Covar_Snow
   Bt_67um_Bt_133um_Covar = Bt_67um_Bt_133um_Covar_Snow

   Bt_11um_Btd_11um_67um_Covar = Bt_11um_Btd_11um_67um_Covar_Snow
   Bt_11um_Btd_11um_85um_Covar = Bt_11um_Btd_11um_85um_Covar_Snow
   Bt_11um_Btd_11um_12um_Covar = Bt_11um_Btd_11um_12um_Covar_Snow
   Bt_11um_Btd_11um_133um_Covar = Bt_11um_Btd_11um_133um_Covar_Snow

   Btd_11um_67um_Btd_11um_67um_Covar = Btd_11um_67um_Btd_11um_67um_Covar_Snow
   Btd_11um_85um_Btd_11um_85um_Covar = Btd_11um_85um_Btd_11um_85um_Covar_Snow
   Btd_11um_12um_Btd_11um_12um_Covar = Btd_11um_12um_Btd_11um_12um_Covar_Snow
   Btd_11um_133um_Btd_11um_133um_Covar = Btd_11um_133um_Btd_11um_133um_Covar_Snow
   Btd_11um_12um_Btd_11um_133um_Covar = Btd_11um_12um_Btd_11um_133um_Covar_Snow
   Btd_11um_12um_Btd_11um_85um_Covar = Btd_11um_12um_Btd_11um_85um_Covar_Snow
   Btd_11um_12um_Btd_11um_67um_Covar = Btd_11um_12um_Btd_11um_67um_Covar_Snow
   Btd_11um_67um_Btd_11um_133um_Covar = Btd_11um_67um_Btd_11um_133um_Covar_Snow
 endif

 !--- Desert
 if (Sfc_Type_Forward_Model == 3) then 
   Bt_67um_Mean = Bt_67um_Mean_Desert
   Bt_85um_Mean = Bt_85um_Mean_Desert
   Bt_11um_Mean = Bt_11um_Mean_Desert
   Bt_12um_Mean = Bt_12um_Mean_Desert
   Bt_133um_Mean = Bt_133um_Mean_Desert
   
   Bt_67um_Bt_67um_Covar = Bt_67um_Bt_67um_Covar_Desert
   Bt_85um_Bt_85um_Covar = Bt_85um_Bt_85um_Covar_Desert
   Bt_11um_Bt_11um_Covar = Bt_11um_Bt_11um_Covar_Desert
   Bt_12um_Bt_12um_Covar = Bt_12um_Bt_12um_Covar_Desert
   Bt_133um_Bt_133um_Covar = Bt_133um_Bt_133um_Covar_Desert

   Bt_85um_Bt_133um_Covar = Bt_85um_Bt_133um_Covar_Desert
   Bt_12um_Bt_133um_Covar = Bt_12um_Bt_133um_Covar_Desert
   Bt_11um_Bt_67um_Covar = Bt_11um_Bt_67um_Covar_Desert
   Bt_11um_Bt_85um_Covar = Bt_11um_Bt_85um_Covar_Desert
   Bt_12um_Bt_85um_Covar = Bt_12um_Bt_85um_Covar_Desert
   Bt_12um_Bt_67um_Covar = Bt_12um_Bt_67um_Covar_Desert
   Bt_11um_Bt_12um_Covar = Bt_11um_Bt_12um_Covar_Desert
   Bt_11um_Bt_133um_Covar = Bt_11um_Bt_133um_Covar_Desert
   Bt_67um_Bt_133um_Covar = Bt_67um_Bt_133um_Covar_Desert

   Bt_11um_Btd_11um_67um_Covar = Bt_11um_Btd_11um_67um_Covar_Desert
   Bt_11um_Btd_11um_85um_Covar = Bt_11um_Btd_11um_85um_Covar_Desert
   Bt_11um_Btd_11um_12um_Covar = Bt_11um_Btd_11um_12um_Covar_Desert
   Bt_11um_Btd_11um_133um_Covar = Bt_11um_Btd_11um_133um_Covar_Desert

   Btd_11um_67um_Btd_11um_67um_Covar = Btd_11um_67um_Btd_11um_67um_Covar_Desert
   Btd_11um_85um_Btd_11um_85um_Covar = Btd_11um_85um_Btd_11um_85um_Covar_Desert
   Btd_11um_12um_Btd_11um_12um_Covar = Btd_11um_12um_Btd_11um_12um_Covar_Desert
   Btd_11um_133um_Btd_11um_133um_Covar = Btd_11um_133um_Btd_11um_133um_Covar_Desert
   Btd_11um_12um_Btd_11um_133um_Covar = Btd_11um_12um_Btd_11um_133um_Covar_Desert
   Btd_11um_12um_Btd_11um_85um_Covar = Btd_11um_12um_Btd_11um_85um_Covar_Desert
   Btd_11um_12um_Btd_11um_67um_Covar = Btd_11um_12um_Btd_11um_67um_Covar_Desert
   Btd_11um_67um_Btd_11um_133um_Covar = Btd_11um_67um_Btd_11um_133um_Covar_Desert
 endif
 !--- Arctic
 if (Sfc_Type_Forward_Model == 4) then
   Bt_67um_Mean = Bt_67um_Mean_Arctic
   Bt_85um_Mean = Bt_85um_Mean_Arctic
   Bt_11um_Mean = Bt_11um_Mean_Arctic
   Bt_12um_Mean = Bt_12um_Mean_Arctic
   Bt_133um_Mean = Bt_133um_Mean_Arctic
   
   Bt_67um_Bt_67um_Covar = Bt_67um_Bt_67um_Covar_Arctic
   Bt_85um_Bt_85um_Covar = Bt_85um_Bt_85um_Covar_Arctic
   Bt_11um_Bt_11um_Covar = Bt_11um_Bt_11um_Covar_Arctic
   Bt_12um_Bt_12um_Covar = Bt_12um_Bt_12um_Covar_Arctic
   Bt_133um_Bt_133um_Covar = Bt_133um_Bt_133um_Covar_Arctic

   Bt_85um_Bt_133um_Covar = Bt_85um_Bt_133um_Covar_Arctic
   Bt_12um_Bt_133um_Covar = Bt_12um_Bt_133um_Covar_Arctic
   Bt_11um_Bt_67um_Covar = Bt_11um_Bt_67um_Covar_Arctic
   Bt_11um_Bt_85um_Covar = Bt_11um_Bt_85um_Covar_Arctic
   Bt_12um_Bt_85um_Covar = Bt_12um_Bt_85um_Covar_Arctic
   Bt_12um_Bt_67um_Covar = Bt_12um_Bt_67um_Covar_Arctic
   Bt_11um_Bt_12um_Covar = Bt_11um_Bt_12um_Covar_Arctic
   Bt_11um_Bt_133um_Covar = Bt_11um_Bt_133um_Covar_Arctic
   Bt_67um_Bt_133um_Covar = Bt_67um_Bt_133um_Covar_Arctic

   Bt_11um_Btd_11um_67um_Covar = Bt_11um_Btd_11um_67um_Covar_Arctic
   Bt_11um_Btd_11um_85um_Covar = Bt_11um_Btd_11um_85um_Covar_Arctic
   Bt_11um_Btd_11um_12um_Covar = Bt_11um_Btd_11um_12um_Covar_Arctic
   Bt_11um_Btd_11um_133um_Covar = Bt_11um_Btd_11um_133um_Covar_Arctic

   Btd_11um_67um_Btd_11um_67um_Covar = Btd_11um_67um_Btd_11um_67um_Covar_Arctic
   Btd_11um_85um_Btd_11um_85um_Covar = Btd_11um_85um_Btd_11um_85um_Covar_Arctic
   Btd_11um_12um_Btd_11um_12um_Covar = Btd_11um_12um_Btd_11um_12um_Covar_Arctic
   Btd_11um_133um_Btd_11um_133um_Covar = Btd_11um_133um_Btd_11um_133um_Covar_Arctic
   Btd_11um_12um_Btd_11um_133um_Covar = Btd_11um_12um_Btd_11um_133um_Covar_Arctic
   Btd_11um_12um_Btd_11um_85um_Covar = Btd_11um_12um_Btd_11um_85um_Covar_Arctic
   Btd_11um_12um_Btd_11um_67um_Covar = Btd_11um_12um_Btd_11um_67um_Covar_Arctic
   Btd_11um_67um_Btd_11um_133um_Covar = Btd_11um_67um_Btd_11um_133um_Covar_Arctic
 endif
 !--- Antarctic
 if (Sfc_Type_Forward_Model == 5) then
   Bt_67um_Mean = Bt_67um_Mean_Antarctic
   Bt_85um_Mean = Bt_85um_Mean_Antarctic
   Bt_11um_Mean = Bt_11um_Mean_Antarctic
   Bt_12um_Mean = Bt_12um_Mean_Antarctic
   Bt_133um_Mean = Bt_133um_Mean_Antarctic
   
   Bt_67um_Bt_67um_Covar = Bt_67um_Bt_67um_Covar_Antarctic
   Bt_85um_Bt_85um_Covar = Bt_85um_Bt_85um_Covar_Antarctic
   Bt_11um_Bt_11um_Covar = Bt_11um_Bt_11um_Covar_Antarctic
   Bt_12um_Bt_12um_Covar = Bt_12um_Bt_12um_Covar_Antarctic
   Bt_133um_Bt_133um_Covar = Bt_133um_Bt_133um_Covar_Antarctic

   Bt_85um_Bt_133um_Covar = Bt_85um_Bt_133um_Covar_Antarctic
   Bt_12um_Bt_133um_Covar = Bt_12um_Bt_133um_Covar_Antarctic
   Bt_11um_Bt_67um_Covar = Bt_11um_Bt_67um_Covar_Antarctic
   Bt_11um_Bt_85um_Covar = Bt_11um_Bt_85um_Covar_Antarctic
   Bt_12um_Bt_85um_Covar = Bt_12um_Bt_85um_Covar_Antarctic
   Bt_12um_Bt_67um_Covar = Bt_12um_Bt_67um_Covar_Antarctic
   Bt_11um_Bt_12um_Covar = Bt_11um_Bt_12um_Covar_Antarctic
   Bt_11um_Bt_133um_Covar = Bt_11um_Bt_133um_Covar_Antarctic
   Bt_67um_Bt_133um_Covar = Bt_67um_Bt_133um_Covar_Antarctic

   Bt_11um_Btd_11um_67um_Covar = Bt_11um_Btd_11um_67um_Covar_Antarctic
   Bt_11um_Btd_11um_85um_Covar = Bt_11um_Btd_11um_85um_Covar_Antarctic
   Bt_11um_Btd_11um_12um_Covar = Bt_11um_Btd_11um_12um_Covar_Antarctic
   Bt_11um_Btd_11um_133um_Covar = Bt_11um_Btd_11um_133um_Covar_Antarctic

   Btd_11um_67um_Btd_11um_67um_Covar = Btd_11um_67um_Btd_11um_67um_Covar_Antarctic
   Btd_11um_85um_Btd_11um_85um_Covar = Btd_11um_85um_Btd_11um_85um_Covar_Antarctic
   Btd_11um_12um_Btd_11um_12um_Covar = Btd_11um_12um_Btd_11um_12um_Covar_Antarctic
   Btd_11um_133um_Btd_11um_133um_Covar = Btd_11um_133um_Btd_11um_133um_Covar_Antarctic
   Btd_11um_12um_Btd_11um_133um_Covar = Btd_11um_12um_Btd_11um_133um_Covar_Antarctic
   Btd_11um_12um_Btd_11um_85um_Covar = Btd_11um_12um_Btd_11um_85um_Covar_Antarctic
   Btd_11um_12um_Btd_11um_67um_Covar = Btd_11um_12um_Btd_11um_67um_Covar_Antarctic
   Btd_11um_67um_Btd_11um_133um_Covar = Btd_11um_67um_Btd_11um_133um_Covar_Antarctic
 endif

end subroutine SET_ClEAR_SKY_COVARIANCE_TERMS

!----------------------------------------------------------------------
!  In GEOCAT, Acha_Mode_Flag can not be passed in, use this routine to 
!  determine it in based on the available channels
!----------------------------------------------------------------------
subroutine DETERMINE_ACHA_MODE_BASED_ON_CHANNELS(symbol, &
                                                 Acha_Mode_Flag, &
                                                 Chan_On_67um, &
                                                 Chan_On_85um, &
                                                 Chan_On_11um, &
                                                 Chan_On_12um, &
                                                 Chan_On_133um)

  type(symbol_acha), intent(in) :: symbol
  integer, intent(INOUT):: Acha_Mode_Flag
  integer, intent(in):: Chan_On_67um
  integer, intent(in):: Chan_On_85um
  integer, intent(in):: Chan_On_11um
  integer, intent(in):: Chan_On_12um
  integer, intent(in):: Chan_On_133um

  if (Acha_Mode_Flag == -1) then
     if (Chan_On_11um == symbol%YES .and. Chan_On_12um == symbol%YES) then
        if (Chan_On_133um == symbol%YES) then
            Acha_Mode_Flag = 3          ! 11/12/13.3 um
        elseif (Chan_On_85um == symbol%YES) then
            Acha_Mode_Flag = 4          ! 8.5/11/12 um
        elseif (Chan_On_67um == symbol%YES) then
            Acha_Mode_Flag = 5          ! 6.7/11/12 um
        else
            Acha_Mode_Flag = 1          ! 11/12 um
        endif
     endif
  endif
  if (Acha_Mode_Flag == -1) then
     if (Chan_On_12um == symbol%NO) then
        if (Chan_On_67um == symbol%YES .and. Chan_On_133um == symbol%YES) then
              Acha_Mode_Flag = 6       !6.7/11/13.3
        endif
        if (Chan_On_67um == symbol%NO .and. Chan_On_133um == symbol%YES) then
              Acha_Mode_Flag = 2       !11/13.3
        endif
     endif
  endif
end subroutine  DETERMINE_ACHA_MODE_BASED_ON_CHANNELS

!----------------------------------------------------------------------------
! Function INTERPOLATE_PROFILE_ACHA
!
! general interpolation routine for profiles
!
! input:
! lonx - longitude weighting factor
! latx = latitude weighting factor
! z1 = data(ilon, ilat)
! z2 = data(ilonx,ilat)
! z3 = data(ilon,ilatx)
! z4 = data(ilonx,ilatx)
!
! output:
! z = interpolated profile
!
!
!---------------------------------------------------------------------------
 function INTERPOLATE_PROFILE_ACHA(z1,z2,z3,z4,lonx,latx) result(z)

  real, dimension(:), intent(in):: z1
  real, dimension(:), intent(in):: z2
  real, dimension(:), intent(in):: z3
  real, dimension(:), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  real, dimension(size(z1)):: z

  !--- linear inteprpolation scheme
  z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
           (lonx) * ((1.0-latx) * z2 + (latx)* z4)

 end function INTERPOLATE_PROFILE_ACHA

!----------------------------------------------------------------------------
! Function INTERPOLATE_NWP_ACHA
!
! general interpolation routine for nwp fields
!
! description of arguments
! ilon, ilat - nwp indices of closest nwp cell
! ilonx,ilatx - nwp indices of nwp cells of diagnoal of bounding box
! lonx - longitude weighting factor
! latx = latitude weighting factor
! z1 = data(ilon, ilat)
! z2 = data(ilonx,ilat)
! z3 = data(ilon,ilatx)
! z4 = data(ilonx,ilatx)
!
! output:
! z = interpolated data point
!
!---------------------------------------------------------------------------
 function INTERPOLATE_NWP_ACHA(z1,z2,z3,z4,lonx,latx) result(z)

  real, intent(in):: z1
  real, intent(in):: z2
  real, intent(in):: z3
  real, intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  real:: z

  !--- linear inteprpolation scheme
  z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
           (lonx) * ((1.0-latx) * z2 + (latx)* z4)

 end function INTERPOLATE_NWP_ACHA

 !---------------------------------------------------------------------
 ! find Inversion Level - highest Level Inversion below trop
 !---------------------------------------------------------------------
 function DETERMINE_INVERSION_LEVEL(Tropo_Level, Sfc_Level, Sfc_Air_Temp) result(Inversion_Level)
   integer, intent(in):: Tropo_Level
   integer, intent(in):: Sfc_Level
   real, intent(in):: Sfc_Air_Temp
   integer:: Inversion_Level
   integer:: k

   Inversion_Level = 0

   !--- old
!  do k = Tropo_Level, Sfc_Level-1
!     if ((Temp_Prof_Rtm(k) - Temp_Prof_Rtm(k+1) > Delta_T_Layer_Inversion)
!     .and. &
!         (Press_Prof_Rtm(k) >= P_Inversion_Min)) then
!         Inversion_Level = k
!         exit
!     endif
!  enddo


   !---- new
   !---- 
   do k = Tropo_Level, Sfc_Level-1
      if (Temp_Prof_Rtm(k) > Sfc_Air_Temp) then
          Inversion_Level = minval(maxloc(Temp_Prof_Rtm(k:Sfc_Level-1)))
          Inversion_Level = Inversion_Level - 1 + k
          exit
      endif
   enddo

 end function DETERMINE_INVERSION_LEVEL
 !---------------------------------------------------------------------
 ! Find Opaque Cloud Level - highest Level Inversion below trop
 !---------------------------------------------------------------------
 subroutine DETERMINE_OPAQUE_LEVEL(Radiance_11um, &
                                   Black_Body_Rad_Prof_11um, &
                                   Press_Prof, &
                                   Height_Prof, &
                                   Temp_Prof, &
                                   Tropo_Level, &
                                   Sfc_Level, &
                                   Pc_Opaque_Level, &
                                   Tc_Opaque_Level, &
                                   Zc_Opaque_Level)
                             
   real(kind=real4), intent(in):: Radiance_11um
   real(kind=real4), intent(in), dimension(:):: Black_Body_Rad_Prof_11um
   real(kind=real4), intent(in), dimension(:):: Press_Prof
   real(kind=real4), intent(in), dimension(:):: Height_Prof
   real(kind=real4), intent(in), dimension(:):: Temp_Prof
   integer(kind=int4), intent(in):: Tropo_Level
   integer(kind=int4), intent(in):: Sfc_Level
   real(kind=real4), intent(out):: Pc_Opaque_Level
   real(kind=real4), intent(out):: Tc_Opaque_Level
   real(kind=real4), intent(out):: Zc_Opaque_Level
   integer:: Ilev
   integer:: Ilev_Start
   integer:: Ilev_End

   !--- initialize
   Pc_Opaque_Level =  Missing_Value_Real4
   Zc_Opaque_Level =  Missing_Value_Real4
   Tc_Opaque_Level =  Missing_Value_Real4

   !--- restrict levels to consider
   Ilev_Start = Tropo_Level 
   Ilev_End = Sfc_Level

   !--- loop through levels
   level_loop: do Ilev = Ilev_Start, Ilev_End
      Pc_Opaque_Level = Press_Prof(Ilev-1)
      Zc_Opaque_Level = Height_Prof(Ilev-1)
      Tc_Opaque_Level = Temp_Prof(Ilev-1)
     if (Black_Body_Rad_Prof_11um(Ilev) > Radiance_11um) then
         exit
     endif
   end do Level_Loop

 end subroutine DETERMINE_OPAQUE_LEVEL
 !----------------------------------------------------------------------
 !
 ! Compute the IR Emissivity at a Reference Level
 !
 ! Ref_Level refers to a level index in the profiles
 ! Toa_Radiance = top of atmosphere radiance
 ! Toa_Radiance_Clear = top of atmosphere radiance under clear-skies
 !
 ! 
 ! Black_Body_Rad_Prof_11um_Rtm - this is in memory
 !
 !----------------------------------------------------------------------
 function COMPUTE_REFERENCE_LEVEL_EMISSIVITY(Ref_Level,Toa_Radiance, &
                                             Toa_Radiance_Clear)     &
                                        result(Emissivity_Ref_Level)

    integer(kind=int4), intent(in):: Ref_Level
    real (kind=real4), intent(in):: Toa_Radiance  
    real (kind=real4), intent(in):: Toa_Radiance_Clear
    real (kind=real4):: Emissivity_Ref_Level

    Emissivity_Ref_Level = &
    (Toa_Radiance - Toa_Radiance_Clear) / &
    (Black_Body_Rad_Prof_11um_Rtm(Ref_Level) - Toa_Radiance_Clear)

 end function COMPUTE_REFERENCE_LEVEL_EMISSIVITY

 !----------------------------------------------------------------------
 !  Local Linear Radiative Center
 !----------------------------------------------------------------------
 subroutine LOCAL_LINEAR_RADIATIVE_CENTER(symbol_yes,symbol_no, &
                                          Meander_Flag, &
                                          Grid_Data, &
                                          Element_Start, Number_Of_Elements, & 
                                          Line_Start, Number_Of_Lines, & 
                                          Max_Grad_Distance, &
                                          Min_Grad_Value, &
                                          Max_Grad_Value, &
                                          Grad_Flag,  &
                                          Missing_LRC_Value, &
                                          Skip_LRC_Mask, &
                                          Min_Grid_Data_Valid, Max_Grid_Data_Valid, &
                                          Elem_Idx_LRC, Line_Idx_LRC)

  integer(kind=int1), intent(in) :: symbol_yes
  integer(kind=int1), intent(in) :: symbol_no
  integer, intent(in):: Meander_Flag
  real (kind=real4), intent(in), dimension(:,:) :: Grid_Data
  integer (kind=int4), intent(in):: Element_Start
  integer (kind=int4), intent(in):: Number_of_Elements
  integer (kind=int4), intent(in):: Line_Start
  integer (kind=int4), intent(in):: Number_of_Lines
  integer (kind=int4), intent(in):: Max_Grad_Distance
  real (kind=real4), intent(in):: Min_Grad_Value
  real (kind=real4), intent(in):: Max_Grad_Value
  integer (kind=int4), intent(in):: Grad_Flag
  integer (kind=int4), intent(in):: Missing_LRC_Value
  integer (kind=int1), intent(in), dimension(:,:):: Skip_LRC_Mask
  real (kind=real4), intent(in):: Min_Grid_Data_Valid
  real (kind=real4), intent(in):: Max_Grid_Data_Valid
  integer (kind=int4), intent(out), dimension(:,:):: Elem_Idx_LRC
  integer (kind=int4), intent(out), dimension(:,:):: Line_Idx_LRC
  real, dimension(3,3):: Grad_Array
  integer, dimension(2):: Grad_Indices
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Elem_Idx_Previous
  integer:: Line_Idx_Previous
  integer:: Elem_Idx_Next
  integer:: Line_Idx_Next
  real:: Grad_Temp
  integer:: Element_End
  integer:: Line_End
  integer:: ipoint
  integer:: Elem_Idx_dir
  integer:: Line_Idx_dir
 
  Element_End = Number_of_Elements + Element_Start - 1
  Line_End = Number_of_Lines + Line_Start - 1

  !--- initialize
  Elem_Idx_LRC = Missing_LRC_Value
  Line_Idx_LRC = Missing_LRC_Value

!----------------------------------------------------------------------
! loop through pixels in segment
!----------------------------------------------------------------------
Element_Loop:  do Elem_Idx = Element_Start+1, Element_End-1
Line_Loop:    do Line_Idx = Line_Start+1, Line_End-1

      !--- skip data due to mask
      if (Skip_LRC_Mask(Elem_Idx,Line_Idx) == symbol_YES) cycle

      !-- check for out of bounds data
      if (Grad_Flag ==  1 .and. Grid_Data(Elem_Idx,Line_Idx) < Min_Grid_Data_Valid) cycle
      if (Grad_Flag ==  -1 .and. Grid_Data(Elem_Idx,Line_Idx) > Max_Grid_Data_Valid) cycle

      !-- check for data that already meets LRC criteria
      if ((Grad_Flag ==  1 .and. Grid_Data(Elem_Idx,Line_Idx) > Max_Grid_Data_Valid) .or. &
          (Grad_Flag ==  -1 .and. Grid_Data(Elem_Idx,Line_Idx) < Min_Grid_Data_Valid)) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx
              cycle
      endif

      !--- initialize previous variables
      Elem_Idx_Previous = Elem_Idx
      Line_Idx_Previous = Line_Idx

      !---- go long gradient and check for a reversal or saturation
Gradient_Loop:    do ipoint = 1,Max_Grad_Distance

        !--- compute local gradient, find strongest gradient in 3x3 array and compute direction
        if (ipoint == 1 .or. Meander_Flag == symbol_YES) then

         !--- construct 3x3 array for analysis
         Grad_Array =  &
           Grid_Data(Elem_Idx_Previous-1:Elem_Idx_Previous+1,Line_Idx_Previous-1:Line_Idx_Previous+1) -  &
           Grid_Data(Elem_Idx_Previous,Line_Idx_Previous)

         !--- look for bad data
         if (minval(Grad_Array) == Missing_Value_Real4) exit 

         !--- compute local gradients, find strongest gradient
         if (Grad_Flag == 1) then
          Grad_Indices = maxloc(Grad_Array)
         else
          Grad_Indices = minloc(Grad_Array)
         endif 

         !--- compute direction
         Elem_Idx_Dir = Grad_Indices(1) - 2
         Line_Idx_Dir = Grad_Indices(2) - 2

         !--- check for pixels that are located at  minima/maxima
         if (Elem_Idx_Dir == 0 .and. Line_Idx_Dir == 0) then
           Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Previous
           Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Previous
           exit
         endif

         !--- on first step, only proceed if gradient magnitude exceeds a threshold
         if (ipoint == 1) then
            if (abs(Grad_Array(Grad_Indices(1),Grad_Indices(2))) < Min_Grad_Value) then
              exit
            endif
         endif

         !--- check for going up to steep of a gradient
         if (abs(Grad_Array(Grad_Indices(1),Grad_Indices(2))) > Max_Grad_Value) then
           exit
         endif

        endif

        !-- select next point on the path
        Elem_Idx_Next = Elem_Idx_Previous + Elem_Idx_Dir
        Line_Idx_Next = Line_Idx_Previous + Line_Idx_Dir

        !--- check for hitting segment boundaries
        if (Elem_Idx_Next == Element_Start .or. Elem_Idx_Next == Element_End .or. &
             Line_Idx_Next == Line_Start .or. Line_Idx_Next == Line_End) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Previous
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Previous
              exit
         endif

         !--- check for hitting bad data
         if (Skip_LRC_Mask(Elem_Idx_Next,Line_Idx_Next) == symbol_YES) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Previous
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Previous
              exit
         endif

         !--- check for sign reversal
         if (Meander_Flag == symbol_NO) then

          Grad_Temp = Grid_Data(Elem_Idx_Next,Line_Idx_Next) -  &
                      Grid_Data(Elem_Idx_Previous,Line_Idx_Previous)

          if (Grad_Flag * Grad_Temp < 0) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Previous
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Previous
              exit
          endif
         endif

         !--- check for saturation
         if (Grad_Flag == 1 .and. Grid_Data(Elem_Idx_Next,Line_Idx_Next) > Max_Grid_Data_Valid) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Next
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Next
              exit
         endif
         if (Grad_Flag == -1 .and. Grid_Data(Elem_Idx_Next,Line_Idx_Next) < Min_Grid_Data_Valid) then
              Elem_Idx_LRC(Elem_Idx,Line_Idx) = Elem_Idx_Next
              Line_Idx_LRC(Elem_Idx,Line_Idx) = Line_Idx_Next
              exit
         endif

         !--- store position
         Elem_Idx_Previous = Elem_Idx_Next
         Line_Idx_Previous = Line_Idx_Next

      enddo Gradient_Loop

    end do Line_Loop
  end do Element_Loop

end subroutine LOCAL_LINEAR_RADIATIVE_CENTER
!----------------------------------------------------------------------
! Local Routine for a Standard Deviation
!
! Data_Array - input array of real numbers
! Invalid_Mask = 0 for good pixels, 1 for invalid pixels
! Stddev = Standard Deviation for valid pixels in Data Array
!
! Num_Good = number of valid data point in array
!
! If Num_Good < 2, we do nothing
!----------------------------------------------------------------------
function COMPUTE_STANDARD_DEVIATION(Data_Array,Invalid_Mask) Result(Stddev_of_Array_r4)
   real(kind=real4), dimension(:,:), intent(in):: Data_Array
   integer(kind=int1), dimension(:,:), intent(in):: Invalid_Mask
   real:: Stddev_of_Array_r4
   real(kind=real8):: Stddev_of_Array_r8
   real(kind=real8):: Data_Sum
   real(kind=real8):: Data_Sum_Squared
   real(kind=real8):: Num_Good
   real(kind=real8):: temp

   Num_Good = real(sum(1 - Invalid_Mask))

   if (Num_Good == 0) then
    Stddev_of_Array_r8 = Missing_Value_Real4
   elseif (Num_Good == 1) then
    Stddev_of_Array_r8 = 0.0
   else
    Data_Sum = sum(Data_Array * (1.0 - Invalid_Mask))
    Data_Sum_Squared = sum((Data_Array*(1.0-Invalid_Mask))**2)
    temp = Data_Sum_Squared / Num_Good - (Data_Sum/Num_Good)**2
    if (temp > 0.0) then
      Stddev_of_Array_r8 = sqrt(temp)
    else
      Stddev_of_Array_r8 = 0.0
    endif
   endif 

    Stddev_of_Array_r4 = real(Stddev_of_Array_r8,kind=real4)

end function

!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
subroutine COMPUTE_TAU_REFF_ACHA(symbol,  &
                                 Beta, &
                                 Cosine_Zenith_Angle, &
                                 Cloud_Phase, &
                                 Ec_Slant, & 
                                 Ec, &
                                 Tau, &
                                 Reff)

   type(symbol_acha), intent(in) :: symbol
   real(kind=real4), intent(in):: Beta
   real(kind=real4), intent(in):: Cosine_Zenith_Angle
   real(kind=real4), intent(in):: Ec_Slant
   integer(kind=int4), intent(in):: Cloud_Phase
   real(kind=real4), intent(out):: Ec
   real(kind=real4), intent(out):: Tau
   real(kind=real4), intent(out):: Reff
   real(kind=real4):: Qe_vis
   real(kind=real4):: Qe_11um
   real(kind=real4):: wo_11um
   real(kind=real4):: g_11um
   real(kind=real4):: Tau_Abs_11um
   real(kind=real4):: Temp_R4
   real(kind=real4):: log10_Reff

   Tau = Missing_Value_Real4
   Reff = Missing_Value_Real4
   Ec = Missing_Value_Real4

   if (Cloud_Phase == symbol%ICE_PHASE) then
    Temp_R4 = A_Re_Beta_FIT_ICE +  &
              B_Re_Beta_FIT_ICE*Beta + &
              C_Re_Beta_FIT_ICE*Beta**2 + &
              D_Re_Beta_FIT_ICE*Beta**3
   else
    Temp_R4 = A_Re_Beta_FIT_WATER +  &
              B_Re_Beta_FIT_WATER*Beta + &
              C_Re_Beta_FIT_WATER*Beta**2 + &
              D_Re_Beta_FIT_WATER*Beta**3
   endif

   Reff = 10.0**(1.0/Temp_R4)

   if (Reff > 0.0) then
     log10_Reff = alog10(Reff)
   else
     return
   endif

   if (Cloud_Phase == symbol%ICE_PHASE) then

     Qe_Vis = A_Qe_065um_FIT_ICE +  &
              B_Qe_065um_FIT_ICE*log10_Reff + &
              C_Qe_065um_FIT_ICE*log10_Reff**2

     Qe_11um = A_Qe_11um_FIT_ICE +  &
               B_Qe_11um_FIT_ICE*log10_Reff + &
               C_Qe_11um_FIT_ICE*log10_Reff**2

     wo_11um = A_wo_11um_FIT_ICE +  &
               B_wo_11um_FIT_ICE*log10_Reff + &
               C_wo_11um_FIT_ICE*log10_Reff**2

     g_11um =  A_g_11um_FIT_ICE +  &
               B_g_11um_FIT_ICE*log10_Reff + &
               C_g_11um_FIT_ICE*log10_Reff**2

   else

     Qe_Vis = A_Qe_065um_FIT_WATER +  &
              B_Qe_065um_FIT_WATER*log10_Reff + &
              C_Qe_065um_FIT_WATER*log10_Reff**2

     Qe_11um = A_Qe_11um_FIT_WATER +  &
               B_Qe_11um_FIT_WATER*log10_Reff + &
               C_Qe_11um_FIT_WATER*log10_Reff**2

     wo_11um = A_wo_11um_FIT_WATER +  &
               B_wo_11um_FIT_WATER*log10_Reff + &
               C_wo_11um_FIT_WATER*log10_Reff**2

     g_11um =  A_g_11um_FIT_WATER +  &
               B_g_11um_FIT_WATER*log10_Reff + &
               C_g_11um_FIT_WATER*log10_Reff**2

   endif

   Tau_Abs_11um = -Cosine_Zenith_Angle*alog(1.0 - Ec_Slant) 

   Ec = 1.0 - exp(-Tau_Abs_11um)
   
   Tau = (Qe_vis / Qe_11um) * Tau_Abs_11um / (1.0 - wo_11um * g_11um)

   !--- set negative values to be missing - added by Y Li
   if (Tau < 0) then
      Tau = Missing_Value_Real4
      Reff= Missing_Value_Real4
   endif

end subroutine COMPUTE_TAU_REFF_ACHA 

!---------------------------------------------------------------------------
! Compute Parallax Correction
!
! This routine generates new lat and lon arrays that are parallax
! corrected based on the cloud height
!
! Input: Senzen - sensor viewing zenith angle (deg) 
!        Senaz  - sensor azimuth angle (deg)
!        Lat - uncorrected latitude (deg)
!        Lon  - uncorrected longitude (deg)
!        Zsfc  - surface elevation (m)
!        Zcld  - cloud height (m)
!
! Output
!       Lat_Pc - corrected latitude
!       Lon_Pc - corrected longitude
!
!---------------------------------------------------------------------------
subroutine PARALLAX_ACHA(Zcld,Zsfc,Lat,Lon,Senzen,Senaz,Lat_Pc,Lon_Pc) 
   real, intent(in), dimension(:,:):: Zcld
   real, intent(in), dimension(:,:):: Zsfc
   real, intent(in), dimension(:,:):: Lat
   real, intent(in), dimension(:,:):: Lon
   real, intent(in), dimension(:,:):: Senzen
   real, intent(in), dimension(:,:):: Senaz
   real, intent(out), dimension(:,:):: Lat_Pc
   real, intent(out), dimension(:,:):: Lon_Pc
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: Num_Elem
   integer:: Num_Line
   real:: Total_Displacement
   real:: Delta_Lon
   real:: Delta_Lat
   real:: Lon_Spacing_Per_m
   real,parameter:: Lat_Spacing_Per_m = 9.0909e-06   ! ( = 1.0/110000.0 m )

   Num_Elem = size(Zcld,1) 
   Num_Line = size(Zcld,2)

   !--- initialize output to standard values
   Lat_Pc = Lat
   Lon_Pc = Lon

   !--- loop over pixels in segment
   element_loop: do Elem_Idx = 1, Num_Elem
    line_loop: do Line_Idx = 1, Num_Line

     !--- check for valid data
     if (Zcld(Elem_Idx,Line_Idx) == Missing_Value_Real4 .or. &
         Senzen(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle

     !--- compute correction
     Total_Displacement = max(0.0,tan(Senzen(Elem_Idx,Line_Idx)*Dtor)* &
                                     (Zcld(Elem_Idx,Line_Idx) - Zsfc(Elem_Idx,Line_Idx)))

     Lon_Spacing_Per_m = Lat_Spacing_Per_m * cos(Lat(Elem_Idx,Line_Idx)*Dtor)

     Delta_Lon = sin(Senaz(Elem_Idx,Line_Idx)*Dtor)*Total_Displacement * Lon_Spacing_Per_m
     Delta_Lat = cos(Senaz(Elem_Idx,Line_Idx)*Dtor)*Total_Displacement * Lat_Spacing_Per_m

     !--- generate output positions
     Lat_Pc(Elem_Idx,Line_Idx) = Lat(Elem_Idx,Line_Idx) + Delta_Lat
     Lon_Pc(Elem_Idx,Line_Idx) = Lon(Elem_Idx,Line_Idx) + Delta_Lon

    enddo line_loop
   enddo element_loop

end subroutine PARALLAX_ACHA 


!----------------------------------------------------------------------
!---   ACHA MODE Check
!---
!--- Input:
!---        Chan_On_67um = channel on/off setting for 6.7 micron channel
!---        Chan_On_85um = channel on/off setting for 8.5 micron channel
!---        Chan_On_11um = channel on/off setting for 11 micron channel
!---        Chan_On_12um = channel on/off setting for 12 micron channel
!---        Chan_On_133um = channel on/off setting for 13.3 micron channel
!---
!---  Output:
!---        Acha_Mode_Output = a mode for ACHA that works for available 
!---                           channels. Will be set to 0-7 if success
!---                           -1 is no option is possible
!---
!--- Mode Definition
!---     0 = 11um only 
!---     1 = 11um and 12um 
!---     2 = 11um and 13.3um 
!---     3 = 11um, 12um and 13
!---     4 = 8.5um, 11um and 12um; 
!---     5 = 6.7um, 11um and 12um 
!---     6 = 6.7um, 11um and 13.3um
!---     7 = 6.7um and 11 um 
!----------------------------------------------------------------------
subroutine SET_ACHA_MODE(symbol, &
                         Chan_On_67um, &
                         Chan_On_85um, &
                         Chan_On_11um, &
                         Chan_On_12um, &
                         Chan_On_133um, &
                         Acha_Mode_Output)

   type(symbol_acha), intent(in) :: symbol
   integer, intent(in) :: Chan_On_67um
   integer, intent(in) :: Chan_On_85um
   integer, intent(in) :: Chan_On_11um
   integer, intent(in) :: Chan_On_12um
   integer, intent(in) :: Chan_On_133um
   integer, intent(out) :: Acha_Mode_Output

   if (Chan_On_11um == symbol%NO) then
      Acha_Mode_Output = -1
      return
   endif

   Acha_Mode_Output = 0

   if (Chan_On_85um == symbol%YES) then
      if (Chan_On_12um == symbol%YES) Acha_Mode_Output = 4
   endif

   if (Chan_On_67um == symbol%YES) then
      Acha_Mode_Output = 7
      if (Chan_On_12um == symbol%YES) Acha_Mode_Output = 5
      if (Chan_On_133um == symbol%YES) Acha_Mode_Output = 6
   endif

   if (Chan_On_133um == symbol%YES) then
      if (Chan_On_67um == symbol%NO) Acha_Mode_Output = 2
      if (Chan_On_12um == symbol%YES) Acha_Mode_Output = 3
   endif

end subroutine SET_ACHA_MODE

!----------------------------------------------------------------------
!--- check that the Acha_Mode is consistent with available channels
!--- if consistent, Acha_Mode_Error_Flag = 0, if not, flag = 1
!----------------------------------------------------------------------
subroutine CHECK_ACHA_MODE(symbol, &
                           Acha_Mode_Input, &
                           Chan_On_67um, &
                           Chan_On_85um, &
                           Chan_On_11um, &
                           Chan_On_12um, &
                           Chan_On_133um, &
                           Acha_Mode_Error_Flag)

   type(symbol_acha), intent(in) :: symbol
   integer, intent(in) :: Acha_Mode_Input
   integer, intent(in) :: Chan_On_67um
   integer, intent(in) :: Chan_On_85um
   integer, intent(in) :: Chan_On_11um
   integer, intent(in) :: Chan_On_12um
   integer, intent(in) :: Chan_On_133um
   integer, intent(out) :: Acha_Mode_Error_Flag

   Acha_Mode_Error_Flag = 0

   if (Chan_On_11um == symbol%NO) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_67um == symbol%NO) .and. &
       ((Acha_Mode_Input == 5) .or. &
        (Acha_Mode_Input == 6) .or. &
        (Acha_Mode_Input == 7))) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_85um == symbol%NO) .and. (Acha_Mode_Input == 4)) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_12um == symbol%NO) .and. &
       ((Acha_Mode_Input == 1) .or. &
        (Acha_Mode_Input == 4) .or. &
        (Acha_Mode_Input == 5))) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_133um == symbol%NO) .and. &
       ((Acha_Mode_Input == 2) .or. &
        (Acha_Mode_Input == 3) .or. &
        (Acha_Mode_Input == 6))) then
       Acha_Mode_Error_Flag = 1
       return
   endif

end subroutine CHECK_ACHA_MODE

!------------------------------------------------------------------------------
! Null Pixel Level Pointers 
!------------------------------------------------------------------------------
subroutine NULL_PIX_POINTERS(Input, Acha_RTM_NWP)

   type(acha_input_struct), intent(inout) :: Input
   type(acha_rtm_nwp_struct), intent(inout) :: Acha_RTM_NWP

   !---- null profile pointers each time
   if (Input%Chan_On_67um == sym%YES) then
     Atm_Rad_Prof_67um_RTM => NULL()
     Atm_Trans_Prof_67um_RTM => NULL()
   endif
   if (Input%Chan_On_85um == sym%YES) then
     Atm_Rad_Prof_85um_RTM => NULL()
     Atm_Trans_Prof_85um_RTM => NULL()
   endif
   if (Input%Chan_On_11um == sym%YES) then
      Atm_Rad_Prof_11um_RTM => NULL()
      Atm_Trans_Prof_11um_RTM => NULL()
      Black_Body_Rad_Prof_11um_RTM =>NULL()
   endif
   if (Input%Chan_On_12um == sym%YES) then
      Atm_Rad_Prof_12um_RTM => NULL()
      Atm_Trans_Prof_12um_RTM => NULL()
   endif
   if (Input%Chan_On_133um == sym%YES) then
      Atm_Rad_Prof_133um_RTM => NULL()
      Atm_Trans_Prof_133um_RTM => NULL()
   endif

   Acha_RTM_NWP%T_Prof => NULL()
   Acha_RTM_NWP%T_Prof_1 => NULL() 
   Acha_RTM_NWP%T_Prof_2 => NULL() 
   Acha_RTM_NWP%T_Prof_3 => NULL()

   Acha_RTM_NWP%Z_Prof => NULL() 
   Acha_RTM_NWP%Z_Prof_1 => NULL() 
   Acha_RTM_NWP%Z_Prof_2 => NULL() 
   Acha_RTM_NWP%Z_Prof_3 => NULL() 

   if (Input%Chan_On_67um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_67um_RTM =>  NULL()
     Acha_RTM_NWP%Atm_Trans_Prof_67um_RTM =>  NULL()
   endif
   if (Input%Chan_On_85um == sym%YES) then
     Acha_RTM_NWP%Atm_Rad_Prof_85um_RTM =>  NULL()
     Acha_RTM_NWP%Atm_Trans_Prof_85um_RTM =>  NULL()
   endif
     
   if (Input%Chan_On_11um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_11um_RTM => NULL()
      Acha_RTM_NWP%Atm_Trans_Prof_11um_RTM => NULL()
      Acha_RTM_NWP%Black_Body_Rad_Prof_11um_RTM => NULL()
   endif
   if (Input%Chan_On_12um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_12um_RTM => NULL()
      Acha_RTM_NWP%Atm_Trans_Prof_12um_RTM => NULL()
   endif
   if (Input%Chan_On_133um == sym%YES) then
      Acha_RTM_NWP%Atm_Rad_Prof_133um_RTM => NULL()
      Acha_RTM_NWP%Atm_Trans_Prof_133um_RTM=> NULL()
   endif
 

end subroutine NULL_PIX_POINTERS
!====================================================================
!  record svn version as a global variable for output to hdf
!====================================================================
subroutine SET_ACHA_VERSION(Acha_Version)
   character(len=*):: Acha_Version
   Acha_Version = "$Id$"
end subroutine SET_ACHA_VERSION

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module AWG_CLOUD_HEIGHT
