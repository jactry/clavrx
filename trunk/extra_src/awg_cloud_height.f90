! $Id:$
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
!
!
! This code is from akh_geocat_algorithms repository
! modified by Andrew Heidinger for use in CLAVR-x for GOES-16 testing
! here are the tags of the code.
! %Id: awg_cloud_height.f90,v 1.13 2011/12/16 21:38:24 wstraka Exp %
! CLAVR-x ACHA CVS tag: awg_cloud_height.f90,v 1.5 2011/10/19 22:00:57 heidinger
!----------------------------------------------------------------------
  use ALGORITHM_MODULE_USAGE

  IMPLICIT NONE

  PUBLIC::  AWG_CLOUD_HEIGHT_ALGORITHM
  PUBLIC::  SET_ACHA_VERSION

  PRIVATE:: SPATIALLY_intERPOLATE_LOWER_CLOUD_POSITION  
  PRIVATE:: KNOWING_P_COMPUTE_T_Z
  PRIVATE:: KNOWING_T_COMPUTE_P_Z
  PRIVATE:: KNOWING_Z_COMPUTE_T_P
  PRIVATE:: GENERIC_PROFILE_intERPOLATION
  PRIVATE:: OPTIMAL_ESTIMATION
  PRIVATE:: COMPUTE_FORWARD_MODEL_AND_KERNEL
  PRIVATE:: COMPUTE_APRIORI_BASED_ON_TYPE
  PRIVATE:: COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO
  PRIVATE:: DETERMINE_SFC_TYPE_FORWARD_MODEL
  PRIVATE:: SET_CLEAR_SKY_COVARIANCE_TERMS
  PRIVATE:: Compute_Sy
  PRIVATE:: Compute_Sy_Based_On_Clear_Sky_Covariance_Andy
  PRIVATE:: Compute_Sy_Based_On_Clear_Sky_Constant_Covariance_Chang
  PRIVATE:: Compute_Sy_Based_On_Clear_Sky_Variational_Covariance_Chang

  PRIVATE:: DETERMINE_ACHA_MODE_BASED_ON_CHANNELS
  PRIVATE:: intERPOLATE_PROFILE_ACHA
  PRIVATE:: intERPOLATE_NWP_ACHA
  PRIVATE:: DETERMINE_INVERSION_LEVEL
  PRIVATE:: DETERMINE_OPAQUE_LEVEL
  PRIVATE:: COMPUTE_REFERENCE_LEVEL_EMISSIVITY
  PUBLIC:: LOCAL_LINEAR_RADIATIVE_CENTER
  PRIVATE:: COMPUTE_STANDARD_DEVIATION

  !--- include parameters for each system here
  !INCLUDE 'awg_cld_hght_clavrx_include_1.inc'
  INCLUDE 'awg_cld_hght_geocat_include_1.inc'

  !--- interpolated profiles
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Temp_Prof_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Press_Prof_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Hght_Prof_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_67um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_73um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_85um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_11um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_12um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_133um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_67um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_73um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_85um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_11um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_12um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_133um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Black_Body_Rad_Prof_11um_Rtm
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_67um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_73um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_85um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_11um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_12um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Rad_Prof_133um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_67um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_73um_Rtm_Unbiased 
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_85um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_11um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_12um_Rtm_Unbiased
  real, PRIVATE, dimension(Num_Levels_Rtm_Prof):: Atm_Trans_Prof_133um_Rtm_Unbiased
  integer, PRIVATE:: Inver_Level_Rtm
  integer, PRIVATE:: Sfc_Level_Rtm
  integer, PRIVATE:: Tropo_Level_Rtm

  real, PRIVATE:: Bt_67um_Mean
  real, PRIVATE:: Bt_85um_Mean
  real, PRIVATE:: Bt_11um_Mean
  real, PRIVATE:: Bt_12um_Mean
  real, PRIVATE:: Bt_133um_Mean

  real, PRIVATE:: Bt_67um_Bt_67um_Covar
  real, PRIVATE:: Bt_85um_Bt_85um_Covar
  real, PRIVATE:: Bt_11um_Bt_11um_Covar
  real, PRIVATE:: Bt_12um_Bt_12um_Covar
  real, PRIVATE:: Bt_133um_Bt_133um_Covar

  real, PRIVATE:: Bt_85um_Bt_133um_Covar
  real, PRIVATE:: Bt_12um_Bt_133um_Covar
  real, PRIVATE:: Bt_11um_Bt_67um_Covar
  real, PRIVATE:: Bt_11um_Bt_85um_Covar
  real, PRIVATE:: Bt_12um_Bt_85um_Covar
  real, PRIVATE:: Bt_12um_Bt_67um_Covar
  real, PRIVATE:: Bt_11um_Bt_12um_Covar
  real, PRIVATE:: Bt_11um_Bt_133um_Covar
  real, PRIVATE:: Bt_67um_Bt_133um_Covar

  real, PRIVATE:: Bt_11um_Btd_11um_67um_Covar
  real, PRIVATE:: Bt_11um_Btd_11um_85um_Covar
  real, PRIVATE:: Bt_11um_Btd_11um_12um_Covar
  real, PRIVATE:: Bt_11um_Btd_11um_133um_Covar

  real, PRIVATE:: Btd_11um_67um_Btd_11um_67um_Covar
  real, PRIVATE:: Btd_11um_85um_Btd_11um_85um_Covar
  real, PRIVATE:: Btd_11um_12um_Btd_11um_12um_Covar
  real, PRIVATE:: Btd_11um_133um_Btd_11um_133um_Covar
  real, PRIVATE:: Btd_11um_12um_Btd_11um_133um_Covar
  real, PRIVATE:: Btd_11um_12um_Btd_11um_67um_Covar
  real, PRIVATE:: Btd_11um_12um_Btd_11um_85um_Covar
  real, PRIVATE:: Btd_11um_67um_Btd_11um_133um_Covar

  CONTAINS

!------------------------------------------------------------------------------
! AWG Cloud Temperature, Cloud Emissivity (11 micron) and Cloud Beta (11/12
! micron)
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
  subroutine SET_ACHA_VERSION()
    !ACHA_Version = "$Id: awg_cloud_height.f90,v 1.13 2011/12/16 21:38:24 wstraka Exp $" 
  end subroutine SET_ACHA_VERSION

!----------------------------------------------------------------------
! modification history
!
! July 2006 - Added beta as an element of x
! October 2006 - Added cloud lapse rate to make Tc more related to true 
!                cloud-top temperature
!
!
!------------------------------------------------------------------------------
  subroutine  AWG_CLOUD_HEIGHT_ALGORITHM()
                              
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Param_Idx
  integer:: i
  integer:: ifail
  integer:: isingular
  integer:: Ilev
  integer:: Inwp
  integer:: Jnwp
  integer:: Ivza
  integer:: ilrc
  integer:: jlrc
  integer:: iter
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
  real:: BT_CL_67um

  !--- ch29 variables
  real:: Rad_Ac_85um
  real:: Trans_Ac_85um
  real:: Rad_Clear_85um
  real:: Bc_85um
  real:: BT_CL_85um

  !--- ch31 variables
  real:: Rad_Ac_11um
  real:: Trans_Ac_11um
  real:: Rad_Clear_11um
  real:: Bc_11um
  real:: BT_CL_11um
  real:: Bt_11um_Lrc

  !--- ch32 variables
  real:: Rad_Ac_12um
  real:: Trans_Ac_12um
  real:: Rad_Clear_12um
  real:: Bc_12um
  real:: BT_CL_12um

  !--- ch33 variables
  real:: Rad_Ac_133um
  real:: Trans_Ac_133um
  real:: Rad_Clear_133um
  real:: Bc_133um
  real:: BT_CL_133um

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
  real:: r4_dummy
  real:: Emiss_11um_Tropo
  integer:: Num_Obs
  integer:: Cloud_Type
  integer:: Cloud_Phase
  integer:: Undetected_Cloud
  integer(kind=int1), dimension(:,:), POintER :: Processing_Order
  integer:: iconverged
  integer:: Sfc_Type_Forward_Model
  integer(kind=int1), dimension(NUM_META_DATA):: Meta_Data_Flags

  !--- 1d-var retrieval arrays
  real (kind=real4), allocatable, dimension(:):: y
  real (kind=real4), allocatable, dimension(:):: f
  real (kind=real4), allocatable, dimension(:):: x
  real (kind=real4), allocatable, dimension(:):: x_Ap
  real (kind=real4), allocatable, dimension(:):: delta_x
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

  !--- local POintERs to global arrays or data structures
  real (kind=real4), POintER, dimension(:,:):: Tc
  real (kind=real4), POintER, dimension(:,:):: Ec
  real (kind=real4), POintER, dimension(:,:):: Beta
  real (kind=real4), POintER, dimension(:,:):: Pc
  real (kind=real4), POintER, dimension(:,:):: Zc
  real (kind=real4), POintER, dimension(:,:):: Tau
  real (kind=real4), POintER, dimension(:,:):: Tc_Uncertainty
  real (kind=real4), POintER, dimension(:,:):: Ec_Uncertainty
  real (kind=real4), POintER, dimension(:,:):: Beta_Uncertainty
  real (kind=real4), POintER, dimension(:,:):: Zc_Uncertainty
  real (kind=real4), POintER, dimension(:,:):: Pc_Uncertainty
  integer (kind=int1), POintER, dimension(:,:):: Cloud_Layer
  integer (kind=int1), POintER, dimension(:,:):: Cloud_Type_Local
  integer (kind=int1), POintER, dimension(:,:):: Cloud_Mask_Local
  real (kind=real4), POintER, dimension(:,:):: Ch20_Surface_Emissivity
  real (kind=real4), POintER, dimension(:,:):: Surface_Elevation
  real (kind=real4), POintER, dimension(:,:):: Latitude
  integer (kind=int1), POintER, dimension(:,:):: Snow_Class
  integer (kind=int1), POintER, dimension(:,:):: Surface_Type
  integer(kind=int4), ALLOCATABLE, dimension(:,:):: Elem_Idx_LRC
  integer(kind=int4), ALLOCATABLE, dimension(:,:):: Line_Idx_LRC
  integer(kind=int1), ALLOCATABLE, dimension(:,:):: Skip_LRC_Mask
  integer (kind=int1), POintER, dimension(:):: Chan_On
  real (kind=real4), POintER, dimension(:,:):: Bt_67um
  real (kind=real4), POintER, dimension(:,:):: Bt_85um
  real (kind=real4), POintER, dimension(:,:):: Bt_11um
  real (kind=real4), POintER, dimension(:,:):: Bt_12um
  real (kind=real4), POintER, dimension(:,:):: Bt_133um
  real (kind=real4):: Bt_11um_Std
  real (kind=real4):: Btd_11um_67um_Std
  real (kind=real4):: Btd_11um_85um_Std
  real (kind=real4):: Btd_11um_12um_Std
  real (kind=real4):: Btd_11um_133um_Std
  real (kind=real4), POintER, dimension(:,:):: Rad_11um
  real (kind=real4), POintER, dimension(:,:):: Cosine_Zenith_Angle
  integer(kind=int4), POintER, dimension(:,:):: Elem_Idx_NWP
  integer(kind=int4), POintER, dimension(:,:):: Line_Idx_NWP 
  integer(kind=int4), POintER, dimension(:,:):: Elem_Idx_Opposite_Corner_NWP 
  integer(kind=int4), POintER, dimension(:,:):: Line_Idx_Opposite_Corner_NWP 
  integer(kind=int4), POintER, dimension(:,:):: Viewing_Zenith_Angle_Idx_Rtm
  real(kind=real4), POintER, dimension(:,:):: Latitude_Interp_Weight_NWP
  real(kind=real4), POintER, dimension(:,:):: Longitude_Interp_Weight_NWP
  integer (kind=int1), POintER, dimension(:,:):: Invalid_Data_Mask
  integer (kind=int1), POintER, dimension(:,:,:):: OE_Qf
  integer (kind=int1), POintER, dimension(:,:):: Qf
  real(kind=real4), POintER, dimension(:,:):: Lower_Cloud_Pressure
  real(kind=real4), POintER, dimension(:,:):: Lower_Cloud_Temperature
  real(kind=real4), POintER, dimension(:,:):: Lower_Cloud_Height
  real(kind=real4):: Surface_Temperature
  real(kind=real4), POintER, dimension(:,:):: Surface_Pressure
  real(kind=real4):: Tropopause_Temperature
  real(kind=real4), POintER, dimension(:,:):: Rad_Clear_67um_Local
  real(kind=real4), POintER, dimension(:,:):: Rad_Clear_85um_Local
  real(kind=real4), POintER, dimension(:,:):: Rad_Clear_11um_Local
  real(kind=real4), POintER, dimension(:,:):: Rad_Clear_12um_Local
  real(kind=real4), POintER, dimension(:,:):: Rad_Clear_133um_Local
  integer (kind=int1), POintER, dimension(:,:) :: Inversion_Present_Flag !WCS3
  integer(kind=int4), POintER, dimension(:,:):: Elem_Idx_LRC_CLAVRX
  integer(kind=int4), POintER, dimension(:,:):: Line_Idx_LRC_CLAVRX

  !--- scalar local variables
  integer (kind=int4):: Number_of_Elements
  integer (kind=int4):: Number_Of_Lines
  integer (kind=int4):: Smooth_Nwp_Fields_Flag
  integer (kind=int4):: Acha_Mode_Flag
  integer (kind=int4):: Process_Undetected_Cloud_Flag_Local
  integer (kind=int4):: i1,i2,j1,j2

!-----------------------------------------------------------------------
! BEGIN EXECUTABLE CODE
!-----------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! include code to set up profiles correctly
  !----------------------------------------------------------------------------
  INCLUDE 'awg_cld_hght_clavrx_include_2.inc'

  !---------------------------------------------------------------------------
  !-- Acha Mode set to  -1, determine based on channels
  !---------------------------------------------------------------------------
  if (Geocat_Flag == sym%YES) then
    call DETERMINE_ACHA_MODE_BASED_ON_CHANNELS(Acha_Mode_Flag,Chan_On)
  endif

  !--- determine number of channels
  SELECT CASE(Acha_Mode_Flag)
     CASE(0)  !avhrr-1
       Num_Obs = 1
     CASE(1)  !avhrr, goes-im
       Num_Obs = 2
     CASE(2)  !goes-nop
       Num_Obs = 2
     CASE(3)  !goes-r
       Num_Obs = 3
     CASE(4)  !viirs
       Num_Obs = 3
     CASE(5)  !test
       Num_Obs = 3
     CASE(6)  !test
       Num_Obs = 3
  END SELECT

  !--- allocate needed 2d arrays for processing this segment
  allocate(Elem_Idx_Lrc(Number_of_Elements,Number_of_Lines))
  allocate(Line_Idx_Lrc(Number_of_Elements,Number_of_Lines))
  allocate(Skip_LRC_Mask(Number_of_Elements,Number_of_Lines))

  !--- allocate 1D-VAR arrays based on number of channels
  allocate(y(Num_Obs))
  allocate(y_variance(Num_Obs))
  allocate(f(Num_Obs))
  allocate(x(num_param))
  allocate(x_Ap(num_param))
  allocate(delta_x(num_param))
  allocate(K(Num_Obs,num_param))
  allocate(Sa(num_param,num_param))
  allocate(Sa_inv(num_param,num_param))
  allocate(Sx(num_param,num_param))
  allocate(Sx_inv(num_param,num_param))
  allocate(E(num_param,num_param))
  allocate(Sy(Num_Obs,Num_Obs))
  allocate(Sy_inv(Num_Obs,Num_Obs))
  allocate(Emiss_Vector(Num_Obs))

  !--- turn on diagnostic output
  idiag_output = sym%NO

  !--- set convergence criterion
  conv_crit = num_param - 1.0

  !----------- make identity matrix
  E = 0.0
  do i = 1,num_param
    E(i,i) = 1.0
  enddo

  !--- initialize output
  Tc =  Missing_Value_Real4
  Ec =  Missing_Value_Real4
  Beta =  Missing_Value_Real4
  Pc =  Missing_Value_Real4
  Zc =  Missing_Value_Real4
  Elem_Idx_LRC = Missing_Value_Int4
  Line_Idx_LRC = Missing_Value_Int4
  Qf = 0
  Cloud_Layer = 0
  Meta_Data_Flags = 0
  !--------------------------------------------------------------------------
  ! spatial processing pixels
  ! compute local radiative centers using 11 um brightness temperature
  !---------------------------------------------------------------------------

  !--- construct a mask to select pixel for LRC computation
  Skip_LRC_Mask = Invalid_Data_Mask
  where(Cloud_Mask_Local == sym%CLEAR .or. Cloud_Mask_Local == sym%PROB_CLEAR)
          Skip_LRC_Mask = sym%YES
  endwhere

  !--- call LRC routine
  if (Geocat_Flag == sym%YES) then

    if (Use_Lrc_Flag == sym%YES) then

    call LOCAL_LINEAR_RADIATIVE_CENTER(LRC_Meander_Flag, &
                                       Bt_11um, &
                                       Element_Idx_Min, &
                                       Number_Of_Elements, & 
                                       Line_Idx_Min,  &
                                       Number_Of_Lines, & 
                                       Max_LRC_Distance,  &
                                       Grad_Flag_LRC,  &
                                       Missing_Value_Int4, &
                                       Skip_LRC_Mask, &
                                       Min_Bt_11um_LRC,  &
                                       Max_Bt_11um_LRC, &
                                       Elem_Idx_LRC,  &
                                       Line_Idx_LRC)
    endif

  else

    !-----------------------------------------------------------------------
    ! include code to access LRC indices in global memory in CLAVR-x
    !-----------------------------------------------------------------------
    Elem_Idx_LRC = Elem_Idx_LRC_CLAVRX
    Line_Idx_LRC = Line_Idx_LRC_CLAVRX

  endif

  !--------------------------------------------------------------------------
  ! determine processing order of pixels
  !--------------------------------------------------------------------------
  Processing_Order = 0

  if (Use_Lrc_Flag == sym%NO) then
      ipass_min = 0
      ipass_Max = 0
      Meta_Data_Flags(4) = sym%NO
  else 
      ipass_min = 1
      ipass_Max = 4
      Meta_Data_Flags(4) = sym%YES

      !--- loop through pixels, determine processing order
      Line_loop_0: do Line_Idx = Line_Idx_Min,Number_Of_Lines + Line_Idx_Min - 1
        Element_Loop_0:   do Elem_Idx = Element_Idx_Min, Number_of_Elements + Element_Idx_Min - 1

        if (Invalid_Data_Mask(Elem_Idx,Line_Idx) == sym%YES) then
            cycle
        endif

        ilrc = Elem_Idx_LRC(Elem_Idx,Line_Idx)
        jlrc = Line_Idx_LRC(Elem_Idx,Line_Idx)
        Cloud_Type = Cloud_Type_Local(Elem_Idx,Line_Idx)

        !-- on pass 1, do single layer lrc's
        if ((Elem_Idx == ilrc) .and. (Line_Idx == jlrc) .and. &
            (Cloud_Type /= sym%OVERLAP_TYPE)) then
             Processing_Order(Elem_Idx,Line_Idx) = 1
        endif

        !-- on pass 2, do non-lrc water clouds
        if ((Elem_Idx /= ilrc) .and. (Line_Idx /= jlrc) .and. &
            (Cloud_Type == sym%FOG_TYPE .or. &
            Cloud_Type == sym%WATER_TYPE .or. &
            Cloud_Type == sym%MIXED_TYPE .or. &
            Cloud_Type == sym%SUPERCOOLED_TYPE)) then
            Processing_Order(Elem_Idx,Line_Idx) = 2
        endif

        !-- on pass 3, do lrc overlap clouds
        if ((Elem_Idx == ilrc) .and. (Line_Idx == jlrc) .and. &
            (Cloud_Type == sym%OVERLAP_TYPE)) then
             Processing_Order(Elem_Idx,Line_Idx) = 3
        endif

        !--  on pass-4 do remaining
        if (Processing_Order(Elem_Idx,Line_Idx) == 0) then
           Processing_Order(Elem_Idx,Line_Idx) = 4
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
     call SPATIALLY_intERPOLATE_LOWER_CLOUD_POSITION(ipass,Line_Idx_min,Number_Of_Lines, &
                      Elem_Idx_NWP,Line_Idx_NWP, &
                      Invalid_Data_Mask, &
                      Cloud_Type_Local, &
                      Surface_Pressure, &
                      Pc, &
                      Lower_Cloud_Pressure, &
                      Lower_Cloud_Temperature, &
                      Lower_Cloud_Height)
   endif

   !--------------------------------------------------------------------------
   ! loop over pixels in scanlines
   !--------------------------------------------------------------------------
   Line_loop: do Line_Idx = Line_Idx_min,Number_Of_Lines + Line_Idx_min - 1

    Element_Loop:   do Elem_Idx = 1, Number_of_Elements

    !--- check if pixel should be processd in this path
    if (ipass /= Processing_Order(Elem_Idx,Line_Idx)) then
          cycle
    endif

    !--- check for a bad pixel
    if (Invalid_Data_Mask(Elem_Idx,Line_Idx) == sym%YES) then
          cycle
    endif

    !--- for convenience, save nwp indices to local variables
    Inwp = Elem_Idx_NWP(Elem_Idx,Line_Idx)
    Jnwp = Line_Idx_NWP(Elem_Idx,Line_Idx)
    Ivza =  Viewing_Zenith_Angle_Idx_Rtm(Elem_Idx,Line_Idx)
    ilrc = Elem_Idx_LRC(Elem_Idx,Line_Idx)
    jlrc = Line_Idx_LRC(Elem_Idx,Line_Idx)

    !--- Qc indices
    if (Inwp == Missing_Value_Int4 .or. &
        Jnwp == Missing_Value_Int4 .or. &
        Ivza == Missing_Value_Int4) then 
         cycle 
    endif

    !-----------------------------------------------------------------------
    ! include code to setup local profiles correctly 
    !-----------------------------------------------------------------------
    INCLUDE "awg_cld_hght_clavrx_include_3.inc"

    !-----------------------------------------------------------------------
    !  find opaque levels
    !-----------------------------------------------------------------------
    CALL DETERMINE_OPAQUE_LEVEL(Rad_11um(Elem_Idx,Line_Idx), &
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
        if (((Cloud_Mask_Local(Elem_Idx,Line_Idx) == sym%CLEAR) .or.  &
            (Cloud_Mask_Local(Elem_Idx,Line_Idx) == sym%PROB_CLEAR)) .and. &
            (Process_Undetected_Cloud_Flag_Local == sym%NO)) then
          Tc(Elem_Idx,Line_Idx) = Missing_Value_Real4
          Pc(Elem_Idx,Line_Idx) = Missing_Value_Real4
          Zc(Elem_Idx,Line_Idx) = Missing_Value_Real4
          Ec(Elem_Idx,Line_Idx) = Missing_Value_Real4
          Beta(Elem_Idx,Line_Idx) = Missing_Value_Real4
        else
          Tc(Elem_Idx,Line_Idx) = Tc_Opaque_Level
          Pc(Elem_Idx,Line_Idx) = Pc_Opaque_Level
          Zc(Elem_Idx,Line_Idx) = Zc_Opaque_Level
          Ec(Elem_Idx,Line_Idx) = 1.0
          Beta(Elem_Idx,Line_Idx) = Missing_Value_Real4
        endif
        Tc_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Pc_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Zc_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Ec_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        Beta_Uncertainty(Elem_Idx,Line_Idx) = Missing_Value_Real4
        cycle
   endif


   !----------------------------------------------------------------------
   ! determine cloud phase from cloud type for convienience
   !----------------------------------------------------------------------
   Cloud_Type = Cloud_Type_Local(Elem_Idx,Line_Idx)
   Cloud_Phase = sym%UNKNOWN_PHASE

   if ( (Cloud_Type  == sym%FOG_TYPE) .or. &
       (Cloud_Type  == sym%WATER_TYPE) .or. &
       (Cloud_Type  == sym%SUPERCOOLED_TYPE)) then
       Cloud_Phase = sym%WATER_PHASE
   endif

   if ( (Cloud_Type  == sym%CIRRUS_TYPE) .or. &
        (Cloud_Type  == sym%OVERLAP_TYPE) .or. &
        (Cloud_Type  == sym%TICE_TYPE)) then
        Cloud_Phase = sym%ICE_PHASE
   endif


  !-----------------------------------------------------------------------
  !----- data quality check
  !-----------------------------------------------------------------------
  if ((Bt_11um(Elem_Idx,Line_Idx) < 170.0) .or. &         !begin data check
     (Bt_11um(Elem_Idx,Line_Idx) > 340.0) .or. &
     (Surface_Temperature < 180.0) .or. &
     (Surface_Temperature > 340.0) .or. &
     (Tropopause_Temperature < 160.0) .or. &
     (Tropopause_Temperature > 270.0)) then

     Tc(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Ec(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Beta(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Pc(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Zc(Elem_Idx,Line_Idx) =  Missing_Value_Real4
     Qf(Elem_Idx,Line_Idx) = 0
     Cloud_Layer(Elem_Idx,Line_Idx) = 0

   else  !if passed data check then proceed with retrieval

   !---------------------------------------------------------------------
   ! select to do retrievals for all pixels or just cloudy ones
   !---------------------------------------------------------------------
   Undetected_Cloud = sym%NO

   if ((Cloud_Mask_Local(Elem_Idx,Line_Idx) == sym%CLEAR) .or.  &
       (Cloud_Mask_Local(Elem_Idx,Line_Idx) == sym%PROB_CLEAR)) then

       if (Process_Undetected_Cloud_Flag_Local == sym%NO) then
               cycle
       else
              Undetected_Cloud = sym%YES
       endif

   endif


   !----------------------------------------------------------------------
   !--- Set Meta Data Flags
   !----------------------------------------------------------------------
   Meta_Data_Flags(1) = sym%YES
   if (Cloud_Phase == sym%ICE_PHASE) then
       Meta_Data_Flags(3) = sym%YES
   else
       Meta_Data_Flags(3) = sym%NO
   endif
   if (Use_Lrc_Flag == sym%YES) then
       Meta_Data_Flags(4) = sym%YES
   else
       Meta_Data_Flags(4) = sym%NO
   endif
   if (Cloud_Type == sym%OVERLAP_TYPE) then
       Meta_Data_Flags(5) = sym%YES
   else
       Meta_Data_Flags(5) = sym%NO
   endif
   Meta_Data_Flags(6) = sym%NO     !lower cloud interpolation
   Meta_Data_Flags(7) = sym%NO     !low level inversion


   !-----------------------------------------------------------------------
   ! compute needed channel 3x3 standard deviations
   !-----------------------------------------------------------------------
   j1 = max(Line_Idx_Min, Line_Idx - 1)
   j2 = min(Number_of_Lines, Line_Idx + 1)
   i1 = max(Element_Idx_Min, Elem_Idx - 1) 
   i2 = min(Number_of_Elements, Elem_Idx + 1)

   Bt_11um_Std = COMPUTE_STANDARD_DEVIATION(Bt_11um(i1:i2,j1:j2),Invalid_Data_Mask(i1:i2,j1:j2))

   if (Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6) then
    Btd_11um_67um_Std = COMPUTE_STANDARD_DEVIATION(Bt_11um(i1:i2,j1:j2) - Bt_67um(i1:i2,j1:j2),&
                                                   Invalid_Data_Mask(i1:i2,j1:j2))
   endif
   if (Acha_Mode_Flag == 4) then
    Btd_11um_85um_Std = COMPUTE_STANDARD_DEVIATION(Bt_11um(i1:i2,j1:j2) - Bt_85um(i1:i2,j1:j2), &
                                                  Invalid_Data_Mask(i1:i2,j1:j2))
   endif

   if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag == 3 .or.  &
       Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 5) then
    Btd_11um_12um_Std = COMPUTE_STANDARD_DEVIATION(Bt_11um(i1:i2,j1:j2) - Bt_12um(i1:i2,j1:j2), &
                                                   Invalid_Data_Mask(i1:i2,j1:j2))
   endif

   if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 6) then
    Btd_11um_133um_Std = COMPUTE_STANDARD_DEVIATION(Bt_11um(i1:i2,j1:j2) - Bt_133um(i1:i2,j1:j2), &
                                                    Invalid_Data_Mask(i1:i2,j1:j2))
   endif

   !-----------------------------------------------------------------------
   ! assign values to y and x_Ap
   !----------------------------------------------------------------------

   !--- y - the observation vector
   SELECT CASE(Acha_Mode_Flag)
     CASE(1)
       y(1) = Bt_11um(Elem_Idx,Line_Idx)
       y(2) = Bt_11um(Elem_Idx,Line_Idx) - Bt_12um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
     CASE(2)
       y(1) = Bt_11um(Elem_Idx,Line_Idx)
       y(2) = Bt_11um(Elem_Idx,Line_Idx) - Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_11um_Std**2
       y_variance(2) = Btd_11um_133um_Std**2 
     CASE(3)
       y(1) = Bt_11um(Elem_Idx,Line_Idx)
       y(2) = Bt_11um(Elem_Idx,Line_Idx) - Bt_12um(Elem_Idx,Line_Idx)
       y(3) = Bt_11um(Elem_Idx,Line_Idx) - Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
       y_variance(3) = Btd_11um_133um_Std**2 
     CASE(4)
       y(1) = Bt_11um(Elem_Idx,Line_Idx)
       y(2) = Bt_11um(Elem_Idx,Line_Idx) - Bt_12um(Elem_Idx,Line_Idx)
       y(3) = Bt_11um(Elem_Idx,Line_Idx) - Bt_85um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
       y_variance(3) = Btd_11um_85um_Std**2 
     CASE(5)
       y(1) = Bt_11um(Elem_Idx,Line_Idx)
       y(2) = Bt_11um(Elem_Idx,Line_Idx) - Bt_12um(Elem_Idx,Line_Idx)
       y(3) = Bt_11um(Elem_Idx,Line_Idx) - Bt_67um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_11um_Std**2
       y_variance(3) = Btd_11um_12um_Std**2 
       y_variance(3) = Btd_11um_67um_Std**2 
     CASE(6)
       y(1) = Bt_11um(Elem_Idx,Line_Idx)
       y(2) = Bt_11um(Elem_Idx,Line_Idx) - Bt_133um(Elem_Idx,Line_Idx)
       y(3) = Bt_11um(Elem_Idx,Line_Idx) - Bt_67um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_11um_Std**2
       y_variance(2) = Btd_11um_133um_Std**2 
       y_variance(3) = Btd_11um_67um_Std**2 
     CASE DEFAULT
       y(1) = Bt_11um(Elem_Idx,Line_Idx)
       y(2) = Bt_11um(Elem_Idx,Line_Idx) - Bt_12um(Elem_Idx,Line_Idx)
       y_variance(1) = Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
   END SELECT

   if (idiag_output == sym%YES) then
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
   call DETERMINE_SFC_TYPE_FORWARD_MODEL(Surface_Type(Elem_Idx,Line_Idx), &
                                         Snow_Class(Elem_Idx,Line_Idx), &
                                         Latitude(Elem_Idx,Line_Idx), &
                                         Ch20_Surface_Emissivity(Elem_Idx,Line_Idx), &
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
   if (Surface_Type(Elem_Idx,Line_Idx) == sym%WATER_SFC) then
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

   if (idiag_output == sym%YES) then
           print *, "Clear-sky uncertainties = ", &
           T11um_Clr_Uncer,T11um_85um_Clr_Uncer,T11um_12um_Clr_Uncer, &
           T11um_133um_Clr_Uncer
   endif

  !--------------------------------------------------------------------
  ! pick a priori conditions
  !--------------------------------------------------------------------

  !--- logic for unmasked or untyped pixels (Undetected cloud)
  if (Undetected_Cloud == sym%YES) then
         if (Tc_Opaque_Level < 260.0 .and.  &
             Tc_Opaque_Level /= Missing_Value_Real4) then
             Cloud_Type = sym%CIRRUS_TYPE
         else
             Cloud_Type = sym%FOG_TYPE
         endif
  endif

  !---- select Tc and ec apriori based on cloud type
! call COMPUTE_APRIORI_BASED_ON_TYPE(Cloud_Type, &
!                      Ttropo_Rtm(Elem_Idx,Line_Idx), &
!                      Bt_11um(Elem_Idx,Line_Idx), &
!                      Cosine_Zenith_Angle(Elem_Idx,Line_Idx), &
!                      Tc_Opaque_cloud(Elem_Idx,Line_Idx), &
!                      Tc_Ap, Tc_Ap_Uncer, &
!                      Ec_Ap, Ec_Ap_Uncer, &
!                      Beta_Ap, Beta_Ap_Uncer)

  !---- Compute 11um emissivity referenced to tropopause
  Emiss_11um_Tropo = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( &
                             Tropo_Level_Rtm, &
                             Rad_11um(Elem_Idx,Line_Idx), &
                             Rad_Clear_11um_Local(Elem_Idx,Line_Idx))

  !---- select Tc and ec apriori based on cloud type

  if ((ilrc /= Missing_Value_Int4) .and. &
      (jlrc /= Missing_Value_Int4)) then
          Bt_11um_Lrc = Bt_11um(ilrc,jlrc)
  else
          Bt_11um_Lrc = Missing_Value_Real4
  endif

  call COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO( &
                       Cloud_Phase, &
                       Emiss_11um_Tropo, &
                       Tropopause_Temperature, &
                       Bt_11um(Elem_Idx,Line_Idx), &
                       Bt_11um_Lrc, &
                       Tc_Opaque_Level, &
                       Tc_Ap,Tc_Ap_Uncer, &
                       Ec_Ap,Ec_Ap_Uncer, &
                       Beta_Ap,Beta_Ap_Uncer)

 x_Ap(1) = Tc_Ap
 x_Ap(2) = Ec_Ap
 x_Ap(3) = Beta_Ap

   if (idiag_output == sym%YES) then
           print *, "x_Ap = ", x_Ap
           print *, "x_Ap Uncer = ", Tc_Ap_Uncer,Ec_Ap_Uncer, Beta_Ap_Uncer
   endif

   !---- determine beta fit parameters based on phase (derived from type)
          
   !--- water phase clouds
   IF (Cloud_Phase == sym%WATER_PHASE) then
      a_Beta_11um_133um_fit = a_Beta_11um_133um_fit_Water
      b_Beta_11um_133um_fit = b_Beta_11um_133um_fit_Water
      a_Beta_11um_85um_fit = a_Beta_11um_85um_fit_Water
      b_Beta_11um_85um_fit = b_Beta_11um_85um_fit_Water
      a_Beta_11um_67um_fit = a_Beta_11um_67um_fit_Water
      b_Beta_11um_67um_fit = b_Beta_11um_67um_fit_Water
   ELSE
   !--- ice phase clouds
      a_Beta_11um_133um_fit = a_Beta_11um_133um_fit_Ice
      b_Beta_11um_133um_fit = b_Beta_11um_133um_fit_Ice
      a_Beta_11um_85um_fit = a_Beta_11um_85um_fit_Ice
      b_Beta_11um_85um_fit = b_Beta_11um_85um_fit_Ice
      a_Beta_11um_67um_fit = a_Beta_11um_67um_fit_Ice
      b_Beta_11um_67um_fit = b_Beta_11um_67um_fit_Ice
  END IF
  if (idiag_output == sym%YES) then
        print *, "beta fit for 6.7 = ", a_Beta_11um_67um_fit, b_Beta_11um_67um_fit
        print *, "beta fit for 8.5 = ", a_Beta_11um_85um_fit, b_Beta_11um_85um_fit
        print *, "beta fit for 13.3 = ", a_Beta_11um_133um_fit, b_Beta_11um_133um_fit
  endif

  !--- now compute Sa
  Sa = 0.0
  Sa(1,1) = Tc_Ap_Uncer
  Sa(2,1) = 0.0 
  Sa(2,2) = Ec_Ap_Uncer
  Sa(1,2) = 0.0 
  Sa(3,3) = Beta_Ap_Uncer

  !--- modify a priori values based on lrc
! if ((ilrc /= Missing_Value_Int4) .and. &
!     (jlrc /= Missing_Value_Int4)) then

!       if ((Tc(ilrc,jlrc) /= Missing_Value_Real4) .and. &
!           (Ec(ilrc,jlrc) > 0.00) .and. &
!           (Ec(ilrc,jlrc) <= 1.0)) then

          !-- use lrc value but weight uncertainty
!         x_Ap(1) = Tc(ilrc,jlrc)
!         Sa(1,1) = 5.0 + (1.0-Ec(ilrc,jlrc))*Tc_Ap_Uncer

!      endif
! endif

  !--- square the individual elements to convert to variances (not a matmul)
  Sa = Sa**2

  !--- compute inverse of Sa matrix
  call INVERT_3x3(Sa,Sa_inv,isingular)
  if (isingular == 1) then
    print *, "Cloud Height warning ==> Singular Sa in Split Window", &
           Elem_Idx,Line_Idx, Cloud_Type
    ifail = sym%YES
    exit
  endif

!----------------------------------------------------------------------------
! compute clear-sky radiances
!-----------------------------------------------------------------------------

!---------------------------------------------------------------------------
!--- modify clear radiances to simulate that from an opaque cloud at 200 mb
!--- above the surface when a multi-layer situation is suspected
!---------------------------------------------------------------------------

 if (Cloud_Type == sym%OVERLAP_TYPE .and.  &
    Lower_Cloud_Pressure(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then

    Meta_Data_Flags(6) = sym%YES

    Rad_Ac_11um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                            Hght_Prof_Rtm,Atm_Rad_Prof_11um_Rtm)

    Trans_Ac_11um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                            Hght_Prof_Rtm,Atm_Trans_Prof_11um_Rtm)

    Bc_11um = PLANCK_RAD_FAST(Chan_Idx_11um,Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
    Rad_Clear_11um = Rad_Ac_11um + Trans_Ac_11um*Bc_11um

  if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 5) then
     Rad_Ac_12um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Rad_Prof_12um_Rtm)

     Trans_Ac_12um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Trans_Prof_12um_Rtm)
     Bc_12um = PLANCK_RAD_FAST(Chan_Idx_12um,Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
     Rad_Clear_12um = Rad_Ac_12um + Trans_Ac_12um*Bc_12um
  endif

  if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 6) then
     Rad_Ac_133um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Rad_Prof_133um_Rtm)

     Trans_Ac_133um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Trans_Prof_133um_Rtm)

     Bc_133um = PLANCK_RAD_FAST(Chan_Idx_133um,Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
     Rad_Clear_133um = Rad_Ac_133um + Trans_Ac_133um*Bc_133um
  endif

  if (Acha_Mode_Flag == 4) then
     Rad_Ac_85um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Rad_Prof_85um_Rtm)

     Trans_Ac_85um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Trans_Prof_85um_Rtm)

     Bc_85um = PLANCK_RAD_FAST(Chan_Idx_85um,Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
     Rad_Clear_85um = Rad_Ac_85um + Trans_Ac_85um*Bc_85um
  endif

  if (Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6) then
     Rad_Ac_67um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Rad_Prof_67um_Rtm)

     Trans_Ac_67um = GENERIC_PROFILE_intERPOLATION(Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                               Hght_Prof_Rtm,Atm_Trans_Prof_67um_Rtm)

     Bc_67um = PLANCK_RAD_FAST(Chan_Idx_67um,Lower_Cloud_Temperature(Elem_Idx,Line_Idx))
     Rad_Clear_67um = Rad_Ac_67um + Trans_Ac_67um*Bc_67um
  endif


  Tsfc_Est = Lower_Cloud_Temperature(Elem_Idx,Line_Idx)

 else

 !---------------------------------------------------------------
 ! if not multi-layer, use existing clear sky radiances
 !---------------------------------------------------------------

  Tsfc_Est = Surface_Temperature 
  Rad_Clear_11um = Rad_Clear_11um_Local(Elem_Idx,Line_Idx)
  if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag ==3 .or. Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 5) then
      Rad_Clear_12um = Rad_Clear_12um_Local(Elem_Idx,Line_Idx)
  endif
  if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 6) then
      Rad_Clear_133um = Rad_Clear_133um_Local(Elem_Idx,Line_Idx)
  endif
  if (Acha_Mode_Flag == 4) then
      Rad_Clear_85um = Rad_Clear_85um_Local(Elem_Idx,Line_Idx)
  endif
  if (Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6) then
      Rad_Clear_67um = Rad_Clear_67um_Local(Elem_Idx,Line_Idx)
  endif

 endif

 if (idiag_output == sym%YES) then
    print *, "Clear radiances = ", Rad_Clear_67um, Rad_Clear_85um, Rad_Clear_11um, Rad_Clear_12um, Rad_Clear_133um
 endif


!----------------------------------------------------------------
! Determine the level of the highest inversion (0=if none)
!----------------------------------------------------------------
Inver_Level_Rtm = DETERMINE_INVERSION_LEVEL(Tropo_Level_Rtm, Sfc_Level_Rtm)

!----------------------------------------------------------------
! if no inversion is present, check to see if clear radiance
! is warmer than observed radiance, if not then do not proceed
!----------------------------------------------------------------
if (( Inver_Level_Rtm == 0) .and. (Rad_Clear_11um < Rad_11um(Elem_Idx,Line_Idx))) then

   !-- consider this pixel a failed retrieval
   ifail = sym%YES

   Meta_Data_Flags(1) = sym%NO

   !--- assign default values for this result
   Tc(Elem_Idx,Line_Idx) = x_Ap(1)   !Missing_Value_Real4
   Ec(Elem_Idx,Line_Idx) = x_Ap(2)   !Missing_Value_Real4
   Beta(Elem_Idx,Line_Idx) = x_Ap(3) !Missing_Value_Real4
   Tau(Elem_Idx,Line_Idx) = Missing_Value_Real4
   Qf(Elem_Idx,Line_Idx) = 1

   call KNOWING_T_COMPUTE_P_Z(Pc(Elem_Idx,Line_Idx),Tc(Elem_Idx,Line_Idx),Zc(Elem_Idx,Line_Idx),Ilev,ierror)
   if (ierror == 0) then
    Zc(Elem_Idx,Line_Idx) = max(0.0,Zc(Elem_Idx,Line_Idx))
    Pc(Elem_Idx,Line_Idx) = min(Surface_Pressure(Inwp,Jnwp),Pc(Elem_Idx,Line_Idx))
   endif
   
   !-- go to next pixel
   cycle

endif

!-----------------------------------------------------------------
! start of retrieval loop
!-----------------------------------------------------------------
iter = 0
iconverged = sym%NO
ifail = sym%NO

!---- assign x to the first guess
x = x_Ap

Retrieval_Loop: do

 iter = iter + 1

  !---------------------------------------------------------------------
  ! estimate above cloud radiances and transmissions
  !---------------------------------------------------------------------
  Tc_temp = x(1)

  call KNOWING_T_COMPUTE_P_Z(Pc_temp,Tc_temp,Zc_temp,Ilev,ierror)

  !--- compute above-cloud terms
  Rad_Ac_11um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_11um_Rtm)

  Trans_Ac_11um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_11um_Rtm)

  if (Acha_Mode_Flag == 1 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 5) then
     Rad_Ac_12um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_12um_Rtm)

     Trans_Ac_12um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_12um_Rtm)
  endif

  if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 6) then
    Rad_Ac_133um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_133um_Rtm)

    Trans_Ac_133um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_133um_Rtm)
  endif

  if (Acha_Mode_Flag == 4) then
    Rad_Ac_85um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_85um_Rtm)

    Trans_Ac_85um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_85um_Rtm)
  endif

  if (Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6) then
    Rad_Ac_67um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Rad_Prof_67um_Rtm)

    Trans_Ac_67um = GENERIC_PROFILE_intERPOLATION(Zc_temp, &
                            Hght_Prof_Rtm,Atm_Trans_Prof_67um_Rtm)

  endif

  !--------------------------------------------------
  ! call abi_planck_temp(ichan,rad,bt)
  !--------------------------------------------------
   if (Chan_On(Chan_Idx_67um) == sym%YES) BT_CL_67um = PLANCK_TEMP_FAST(Chan_Idx_67um, Rad_Clear_67um)
   if (Chan_On(Chan_Idx_85um) == sym%YES) BT_CL_85um = PLANCK_TEMP_FAST(Chan_Idx_85um, Rad_Clear_85um)
   if (Chan_On(Chan_Idx_11um) == sym%YES) BT_CL_11um = PLANCK_TEMP_FAST(Chan_Idx_11um, Rad_Clear_11um)
   if (Chan_On(Chan_Idx_12um) == sym%YES) BT_CL_12um = PLANCK_TEMP_FAST(Chan_Idx_12um, Rad_Clear_12um)
   if (Chan_On(Chan_Idx_133um) == sym%YES) BT_CL_133um = PLANCK_TEMP_FAST(Chan_Idx_133um, Rad_Clear_133um)
 
  !--------------------------------------------------
  ! call forward models
  !--------------------------------------------------
  if (idiag_output == sym%YES) then
      print *, "Iter = ", iter
  endif
  call COMPUTE_FORWARD_MODEL_AND_KERNEL(Acha_Mode_Flag, x, &
           Rad_Clear_67um, Rad_Ac_67um, Trans_Ac_67um, &
           Rad_Clear_85um, Rad_Ac_85um, Trans_Ac_85um, &
           Rad_Clear_11um, Rad_Ac_11um, Trans_Ac_11um, &
           Rad_Clear_12um, Rad_Ac_12um, Trans_Ac_12um, &
           Rad_Clear_133um, Rad_Ac_133um, Trans_Ac_133um, &
           a_Beta_11um_133um_fit, b_Beta_11um_133um_fit, &
           a_Beta_11um_85um_fit, b_Beta_11um_85um_fit, &
           a_Beta_11um_67um_fit, b_Beta_11um_67um_fit, &
           f, K,Emiss_Vector,idiag_output)

  !--------------------------------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute the Sy convariance matrix
  !--------------------------------------------------

!----------------Uncomment for baseline---------------------

 ! call  COMPUTE_Sy(Acha_Mode_Flag, &
 !                 x(2), &
 !                 T11um_Cal_Uncer, &
 !                 T11um_12um_Cal_Uncer, &
 !                 T11um_133um_Cal_Uncer, &
 !                 T11um_85um_Cal_Uncer, &
 !                 T11um_67um_Cal_Uncer, &
 !                 T11um_Clr_Uncer, &
 !                 T11um_12um_Clr_Uncer, &
 !                 T11um_133um_Clr_Uncer, &
 !                 T11um_85um_Clr_Uncer, &
 !                 T11um_67um_Clr_Uncer, &
 !                 Bt_11um_Std, &
 !                 Btd_11um_12um_Std, &
 !                 Btd_11um_133um_Std, &
 !                 Btd_11um_85um_Std, &
 !                 Btd_11um_67um_Std, &
 !                 Sy)


!----------------Uncomment above for baseline---------------------



   Call Compute_Sy_Based_On_Clear_Sky_Covariance_Andy(  &
                                                 Emiss_Vector, &
                                                 Acha_Mode_Flag, &
                                                 y_variance, &
                                                 Sy) 

!    Call Compute_Sy_Based_On_Clear_Sky_Constant_Covariance_Chang(  &
!                                                  Emiss_Vector, &
!                                                  Acha_Mode_Flag, &
!                                                  y_variance, &
!                                                  Sy) 

!    Call Compute_Sy_Based_On_Clear_Sky_Variational_Covariance_Chang(  &
!                                                  Emiss_Vector, &
!                                                  Acha_Mode_Flag, &
!                                                  y_variance, &
!                                                  BT_CL_67um,&
!                                                  BT_CL_85um,&
!                                                  BT_CL_11um,&
!                                                  BT_CL_12um,&
!                                                  BT_CL_133um, &
!                                                  Sy) 

  if (idiag_output == sym%YES) then
          print *, "Sa1 = ", Sa(1,1), Sa(1,2), Sa(1,3)
          print *, "Sa2 = ", Sa(2,1), Sa(2,2), Sa(2,3)
          print *, "Sa3 = ", Sa(3,1), Sa(3,2), Sa(3,3)
          print *, "Sy1 = ", Sy(1,:)
          print *, "Sy2 = ", Sy(2,:)
          print *, "shape of Sy = ", shape(Sy), Num_Obs
  endif
  !--------------------------------------------------
  ! call OE routine to advance the iteration
  !--------------------------------------------------
  call OPTIMAL_ESTIMATION(iter,iter_Max,num_param,Num_Obs, &
                         conv_crit,delta_x_Max, &
                         y,f,x,x_Ap,K,Sy,Sa_inv, &
                         Sx,delta_x,iconverged,ifail, &
                         idiag_output)

  !--- check for a failed iteration
  if (ifail == sym%YES) then
     exit
  endif

  !---------------------------------------------------------
  ! update retrieved vector
  !---------------------------------------------------------
  x = x + delta_x

  if (idiag_output==sym%YES) then  
          print *, "x = ", x
  endif

  !--------------------------------------------------------
  ! exit retrieval loop if converged
  !--------------------------------------------------------
  if (iconverged == sym%YES) then
       if (idiag_output==sym%YES) then  
             print *, "convergence acheived ", iter, x
             print *, '  '
       endif
       exit
  endif

  !-------------------------------------------------------
  ! constrain to reasonable values
  !-------------------------------------------------------
  x(1) = max(min_allowable_Tc,min(Tsfc_Est,x(1)))     !should we do this?
  x(2) = max(0.0,min(x(2),1.0))
  x(3) = max(0.8,min(x(3),2.0))

end do Retrieval_Loop

!=================================================================
! Begin Retrieval Post Processing
!=================================================================

!-----------------------------------------------------------------
! Successful Retrieval Post Processing
!-----------------------------------------------------------------
if (ifail == sym%NO) then  !successful retrieval if statement

 !--- save retrieved vector into its output variables
 Tc(Elem_Idx,Line_Idx) = x(1)
 Beta(Elem_Idx,Line_Idx) = x(3)

 !--- save nadir adjusted emissivity and optical depth
 if (x(2) < 1.00) then
   !this is the 11um optical depth adjusted to nadir
   Tau(Elem_Idx,Line_Idx) = -Cosine_Zenith_Angle(Elem_Idx,Line_Idx)*alog(1.0 - x(2)) 
   Ec(Elem_Idx,Line_Idx) = 1.0 - exp(-Tau(Elem_Idx,Line_Idx))
   Tau(Elem_Idx,Line_Idx) = 2.0 * Tau(Elem_Idx,Line_Idx)    !this is now a "visible" optical depth (assuming 2.0 scale factor)
 else
   Tau(Elem_Idx,Line_Idx) = 20.0
   Ec(Elem_Idx,Line_Idx) = 1.0
   Beta(Elem_Idx,Line_Idx) = 1.3
 endif

 !--- save uncertainty estimates
 Tc_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(1,1))
 Ec_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(2,2))
 Beta_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(3,3))

 !--- set quality flag for a successful retrieval
 Qf(Elem_Idx,Line_Idx) = 3

 !--- Estimate height and pressure
 call KNOWING_T_COMPUTE_P_Z(Pc(Elem_Idx,Line_Idx),Tc(Elem_Idx,Line_Idx),Zc(Elem_Idx,Line_Idx),Ilev,ierror)

 !--- compute lapse rate as dT/dZ
 Lapse_Rate =  (Temp_Prof_Rtm(Ilev+1) - Temp_Prof_Rtm(Ilev)) / &
               (Hght_Prof_Rtm(Ilev+1) - Hght_Prof_Rtm(Ilev))

 !--- compute lapse rate as dP/dZ
 Lapse_Rate_dP_dZ =  (Press_Prof_Rtm(Ilev+1) - Press_Prof_Rtm(Ilev)) / &
                     (Hght_Prof_Rtm(Ilev+1) - Hght_Prof_Rtm(Ilev))

 !---  for low clouds over water, force fixed lapse rate estimate of height
 if (Surface_Type(Elem_Idx,Line_Idx) == sym%WATER_SFC .and. &
     Snow_Class(Elem_Idx,Line_Idx) == sym%NO_SNOW .and. &
     Pc(Elem_Idx,Line_Idx) > Min_P_Inversion) then

       Delta_Cld_Temp_Sfc_Temp = Surface_Temperature - Tc(Elem_Idx,Line_Idx)
       Zc(Elem_Idx,Line_Idx) = Delta_Cld_Temp_Sfc_Temp / &
                      LAPSE_RATE_OCEAN + &
                      Surface_Elevation(Elem_Idx,Line_Idx)/1000.0
       call KNOWING_Z_COMPUTE_T_P(Pc(Elem_Idx,Line_Idx),r4_dummy,Zc(Elem_Idx,Line_Idx),Ilev)
      !Lapse_Rate = LAPSE_RATE_OCEAN
       !-- Heidinger Test - based on CALIPSO analysis
       Lapse_Rate = -8.00 - 0.42*Delta_Cld_Temp_Sfc_Temp + 0.013*(Delta_Cld_Temp_Sfc_Temp)**2
       Meta_Data_Flags(7) = sym%YES

 endif

 !-- Compute Height Uncertainty
 IF (Lapse_Rate /= 0.0) THEN
      Zc_Uncertainty(Elem_Idx,Line_Idx) = Tc_Uncertainty(Elem_Idx,Line_Idx) / ABS(Lapse_Rate)
 ENDIF

 !-- Compute Pressure Uncertainty
 IF (Lapse_Rate /= 0.0 .and. Lapse_Rate_dP_dZ /= 0.0) THEN
      Pc_Uncertainty(Elem_Idx,Line_Idx) = Zc_Uncertainty(Elem_Idx,Line_Idx) * ABS(Lapse_Rate_dP_dZ)
 ENDIF

 !--- quality flags of the retrieved parameters
 DO Param_Idx = 1,Num_Param    !loop over parameters
       IF (Sx(Param_Idx,Param_Idx) < 0.111*Sa(Param_Idx,Param_Idx) ) THEN
            OE_Qf(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_1_3_APRIORI_RETREVIAL
       ELSEIF (Sx(Param_Idx,Param_Idx) < 0.444*Sa(Param_Idx,Param_Idx)) THEN
            OE_Qf(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_2_3_APRIORI_RETREVIAL
       ELSE
            OE_Qf(Param_Idx,Elem_Idx,Line_Idx) = CTH_PARAM_LOW_QUALITY_RETREVIAL
       END IF
 END DO

else   
!-----------------------------------------------------------------
! Failed Retrieval Post Processing
!-----------------------------------------------------------------

 !--- set output variables to apriori
 Tc(Elem_Idx,Line_Idx) = x_Ap(1)   !Missing_Value_Real4
 Ec(Elem_Idx,Line_Idx) = x_Ap(2)   !Missing_Value_Real4
 Beta(Elem_Idx,Line_Idx) = x_Ap(3) !Missing_Value_Real4

 !--- set derived parameters to missing
 Tau(Elem_Idx,Line_Idx) = Missing_Value_Real4
 Pc(Elem_Idx,Line_Idx) = Missing_Value_Real4
 Zc(Elem_Idx,Line_Idx) = Missing_Value_Real4

 !--- set quality flags
 OE_Qf(:,Elem_Idx,Line_Idx) = 0
 Qf(Elem_Idx,Line_Idx) = 1

 !--- estimate height and pressure
 call KNOWING_T_COMPUTE_P_Z(Pc(Elem_Idx,Line_Idx),Tc(Elem_Idx,Line_Idx),Zc(Elem_Idx,Line_Idx),Ilev,ierror)

endif                              !end successful retrieval if statement

!------- determine cloud layer based on pressure
Cloud_Layer(Elem_Idx,Line_Idx) = 0
if (Pc(Elem_Idx,Line_Idx) <= 440.0) then
   Cloud_Layer(Elem_Idx,Line_Idx) = 3
elseif (Pc(Elem_Idx,Line_Idx) < 680.0) then
   Cloud_Layer(Elem_Idx,Line_Idx) = 2
else
   Cloud_Layer(Elem_Idx,Line_Idx) = 1
endif

!--- if retrieval done for an undetected pixel, label the qf
if (Undetected_Cloud == sym%YES) then
 Qf(Elem_Idx,Line_Idx) = 2
 Cloud_Layer(Elem_Idx,Line_Idx) = 0
endif

!-----------------------------------------------------------------
! End Retrieval Post Processing
!-----------------------------------------------------------------

 endif     ! ---------- end of data check
 
 
 !------------ Set Inversion present flag with Metadata flag output
 
 !--------------- WCS3 --------------------
 if (Geocat_Flag == sym%YES) then
    Inversion_Present_Flag(Elem_Idx,Line_Idx) = Meta_Data_Flags(7)
 endif
 
 

 end do Element_Loop

end do Line_loop

end do pass_loop

Processing_Order => null()

  !--- deallocate 2D arrays
  deallocate(Elem_Idx_LRC)
  deallocate(Line_Idx_LRC)
  deallocate(Skip_LRC_Mask)

  !--- deallocate 1D-VAR arrays
  deallocate(y)
  deallocate(y_variance)
  deallocate(f)
  deallocate(x)
  deallocate(x_Ap)
  deallocate(delta_x)
  deallocate(K)
  deallocate(Sa)
  deallocate(Sa_inv)
  deallocate(Sx)
  deallocate(Sx_inv)
  deallocate(E)
  deallocate(Sy)
  deallocate(Sy_inv)
  deallocate(Emiss_Vector)

  !--------------------------------------------------------------------------
  ! Update Global Output Variables and destroy local POintERs and memory
  !--------------------------------------------------------------------------

  !----------------------------
  ! scale height to be in meters
  !----------------------------
 where (Zc /= Missing_Value_Real4)
    Zc = 1000.0* Zc
 endwhere

  !---- nullify local pointers
  Tc => null()
  Ec => null()
  Beta => null()
  Pc => null()
  Zc => null()
  Tau => null()
  Chan_On => null()
  Cloud_Layer => null()
  Cloud_Mask_Local => null()
  Cloud_Type_Local => null()
  Latitude => null()
  Snow_Class => null()
  Ch20_Surface_Emissivity => null()
  Surface_Elevation => null()
  Surface_Type => null()
  Elem_Idx_NWP =>  null()
  Line_Idx_NWP => null()
  Elem_Idx_Opposite_Corner_NWP => null()
  Line_Idx_Opposite_Corner_NWP => null()
  Latitude_Interp_Weight_NWP =>  null()
  Longitude_Interp_Weight_NWP => null()
  Viewing_Zenith_Angle_Idx_Rtm => null()
  Rad_11um => null()
  Bt_67um => null()
  Bt_85um => null()
  Bt_11um => null()
  Bt_12um => null()
  Bt_133um => null()
  Cosine_Zenith_Angle => null()
  Lower_Cloud_Pressure => null()
  Lower_Cloud_Temperature => null()
  Lower_Cloud_Height => null()
  Surface_Pressure => null()
  Rad_Clear_67um_Local => null()
  Rad_Clear_85um_Local => null()
  Rad_Clear_11um_Local => null()
  Rad_Clear_12um_Local => null()
  Rad_Clear_133um_Local => null()
  Invalid_Data_Mask => null()
  OE_Qf => null()
  Qf => null()
  Inversion_Present_Flag => null() !WCS3
  ! Elem_Idx_LRC_CLAVRX and Line_Idx_LRC_CLAVRX are only pointed to
  ! in CLAVR-x. Thus, nulling the pointers only needs to happen if this is done
  if (Geocat_Flag == sym%NO) then
      Elem_Idx_LRC_CLAVRX =>  null()
      Line_Idx_LRC_CLAVRX => null()
  endif

end subroutine  AWG_CLOUD_HEIGHT_ALGORITHM

!-------------------------------------------------------------------------------
! Routine to spatially interpret water cloud temperature values to surrounding 
! pixels
!
! input:  interp_flag = 0 (do no interp, assume Zc=2km) /= 0 (do spatial interp)
!-------------------------------------------------------------------------------
  subroutine SPATIALLY_intERPOLATE_LOWER_CLOUD_POSITION(Interp_Flag,Line_Idx_Min, &
                                                        Number_Of_Lines, &
                                                        Elem_Idx_Nwp,Line_Idx_Nwp, &
                                                        Invalid_Data_Mask, &
                                                        Cloud_Type_Local, &
                                                        Surface_Pressure, &
                                                        Cloud_Pressure, &
                                                        Lower_Cloud_Pressure, &
                                                        Lower_Cloud_Temperature, &
                                                        Lower_Cloud_Height)

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
      num_line = size(Cloud_Type_Local,2)      !-----
      Line_start = Line_Idx_min
      Line_end = Number_Of_Lines + Line_Idx_min - 1

      if (Interp_Flag == 0) then 

       Line_loop_1: do Line_Idx = Line_start, Line_end 
         Element_Loop_1: do Elem_Idx = 1, Num_Elem 
           if (Invalid_Data_Mask(Elem_Idx,Line_Idx) == sym%YES) then
               Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Missing_Value_Real4
               cycle
           endif
           IF (Inwp == Missing_Value_Int4 .or. &
               Jnwp == Missing_Value_Int4) THEN
               
               CYCLE
               
           ENDIF
          Inwp = Elem_Idx_Nwp(Elem_Idx,Line_Idx)
          Jnwp = Line_Idx_Nwp(Elem_Idx,Line_Idx)
          Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Surface_Pressure(Inwp,Jnwp) - Pc_Lower_Cloud_offset
        end do Element_Loop_1
       end do Line_loop_1


      else    !if Interp_Flag > 0 then do a spatial interpolation

        !--- set box width
        delem = 1
        dline = 1

        allocate(mask(Num_Elem,num_line))

         mask = 0

         where((Cloud_Type_Local == sym%FOG_TYPE .or. &
             Cloud_Type_Local == sym%WATER_TYPE .or. &
             Cloud_Type_Local == sym%SUPERCOOLED_TYPE) .and. &
             Cloud_Pressure /= Missing_Value_Real4)
             mask = 1
         endwhere

         Line_loop_2: do Line_Idx = Line_start, Line_end 
             Element_Loop_2: do Elem_Idx = 1, Num_Elem 

             if (Cloud_Type_Local(Elem_Idx,Line_Idx)  /= sym%OVERLAP_TYPE) then
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
               (Invalid_Data_Mask(Elem_Idx,Line_Idx) == sym%NO)) then
               Inwp = Elem_Idx_Nwp(Elem_Idx,Line_Idx)
               Jnwp = Line_Idx_Nwp(Elem_Idx,Line_Idx)
               IF (Inwp == Missing_Value_Int4 .or. &
                   Jnwp == Missing_Value_Int4) THEN
                    CYCLE
               ENDIF
               Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Surface_Pressure(Inwp,Jnwp) - Pc_Lower_Cloud_offset
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
            if (Invalid_Data_Mask(Elem_Idx,Line_Idx) == sym%YES) then
               Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Missing_Value_Real4
               Lower_Cloud_Temperature(Elem_Idx,Line_Idx) = Missing_Value_Real4 
               Lower_Cloud_Height(Elem_Idx,Line_Idx) = Missing_Value_Real4 
               cycle
            endif

            !--- if not overlap, set to all missing 
            if (Cloud_Type_Local(Elem_Idx,Line_Idx) /= sym%OVERLAP_TYPE .or. &
                Lower_Cloud_Pressure(Elem_Idx,Line_Idx) == Missing_Value_Real4) then
                Lower_Cloud_Pressure(Elem_Idx,Line_Idx) = Missing_Value_Real4
                Lower_Cloud_Height(Elem_Idx,Line_Idx) = Missing_Value_Real4
                Lower_Cloud_Temperature(Elem_Idx,Line_Idx) = Missing_Value_Real4
                cycle
            endif

            !--- compute T and Z from P
            Inwp = Elem_Idx_Nwp(Elem_Idx,Line_Idx)
            Jnwp = Line_Idx_Nwp(Elem_Idx,Line_Idx)
            call KNOWING_P_COMPUTE_T_Z( Lower_Cloud_Pressure(Elem_Idx,Line_Idx), &
                                        Lower_Cloud_Temperature(Elem_Idx,Line_Idx), &
                                        Lower_Cloud_Height(Elem_Idx,Line_Idx), &
                                        Ilev)
        end do Element_Loop_4
       end do Line_loop_4

   end subroutine SPATIALLY_intERPOLATE_LOWER_CLOUD_POSITION

   !-----------------------------------------------------------------
   ! Interpolate within profiles knowing P to determine T and Z
   !-----------------------------------------------------------------
   subroutine KNOWING_P_COMPUTE_T_Z(P,T,Z,Ilev)

     real, intent(in):: P
     real, intent(out):: T
     real, intent(out):: Z
     integer, intent(out):: Ilev
     real:: dp
     real:: dt
     real:: dz

     !--- interpolate pressure profile
    !call LOCATE(Press_Prof_Rtm,Num_Levels_Rtm_Prof,P,Ilev)
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
   subroutine KNOWING_Z_COMPUTE_T_P(P,T,Z,Ilev)

     real, intent(in):: Z
     real, intent(out):: T
     real, intent(out):: P
     integer, intent(out):: Ilev
     real:: dp
     real:: dt
     real:: dz
     integer:: kstart
     integer:: kend
     integer:: nlevels_temp

     !--- compute region of atmosphere to look
     kstart = Tropo_Level_Rtm
     kend = Sfc_Level_Rtm
     nlevels_temp = kend - kstart + 1

!    print *, "in KNOWING Z ", kstart, kend, nlevels_temp
!    print *, "in KNOWING Z ", Z
!    print *, "in KNOWING Z ", Hght_Prof_Rtm
     !--- interpolate pressure profile
     call LOCATE(Hght_Prof_Rtm(kstart:kend),nlevels_temp,Z,Ilev)
!    print *, "Ilev = ", Ilev
     Ilev = Ilev + kstart - 1
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
   subroutine KNOWING_T_COMPUTE_P_Z(P,T,Z,klev,ierr)
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
     integer:: Level_Within_Inversion_Flag

     !--- initialization
     ierr = sym%NO
     Z = Missing_Value_Real4
     P = Missing_Value_Real4

     !--- test for existence of a valid solution with troposphere
     kstart = Tropo_Level_Rtm
     kend = Sfc_Level_Rtm
     nlevels_temp = kend - kstart + 1

     !--- check to see if warmer than max, than assume at surface
     if (T > maxval(Temp_Prof_Rtm(kstart:kend))) then
         P = Press_Prof_Rtm(kend)
         Z = Hght_Prof_Rtm(kend)
         klev = kend - 1
         ierr = sym%NO
         return
     endif

     !--- check to see if colder than min, than assume tropopause
     if (T < minval(Temp_Prof_Rtm(kstart:kend))) then
         P = Press_Prof_Rtm(kstart)
         Z = Hght_Prof_Rtm(kstart)
         klev = kstart + 1
         ierr = sym%NO
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

    !--- if no solution within an inversion, look above
    if (Level_Within_Inversion_Flag == 0) then
        kstart = Tropo_Level_Rtm
        kend = Sfc_Level_Rtm
        nlevels_temp = kend - kstart + 1
        call LOCATE(Temp_Prof_Rtm(kstart:kend),nlevels_temp,T,klev)
        klev = klev + kstart - 1
        klev = max(1,min(Num_Levels_Rtm_Prof-1,klev))
    endif

    !-- if solution is above trop, set to trop values
    if (klev < Tropo_Level_Rtm) then
        P = Press_Prof_Rtm(Tropo_Level_Rtm)
        Z = Hght_Prof_Rtm(Tropo_Level_Rtm)

    else

        !--- determine derivatives
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

   !--- end above trop check
   endif

   end subroutine KNOWING_T_COMPUTE_P_Z

   !-----------------------------------------------------------------
   ! Interpolate within profiles knowing Z to determine above cloud
   ! radiative terms used in forward model
   !-----------------------------------------------------------------
   function GENERIC_PROFILE_intERPOLATION(X_value,X_profile,Y_profile)  &
            result(Y_value) 

     real, intent(in):: X_value 
     real, dimension(:), intent(in):: X_profile
     real, dimension(:), intent(in):: Y_profile
     real:: Y_value

     integer:: Ilev
     real:: dx
     integer:: nlevels


     nlevels = size(X_profile)

     !--- interpolate pressure profile
     call LOCATE(X_profile,nlevels,X_value,Ilev)
     Ilev = max(1,min(nlevels-1,Ilev))

     dx = X_profile(Ilev+1) - X_profile(Ilev)

     !--- perform interpolation
     if (dx /= 0.0) then
        Y_value = Y_profile(Ilev) +  &
                 (X_value - X_profile(Ilev))  * &
                 (Y_profile(Ilev+1) - Y_profile(Ilev)) / dx
     else
          Y_value = Y_profile(Ilev)
     endif
   end function GENERIC_PROFILE_intERPOLATION


!------------------------------------------------------------------------
! subroutine to compute the iteration in x due to optimal
! estimation
!
! The notation in this routine follows that of Clive Rodgers (1976,2000)
!
! input to this routine:
! iter - the number of the current iteration
! iter_Max - the maximum number of iterations allowed
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
! delta_x_Max - the maximum step allowed for each delta_x value
!
! output of this routine:
! Sx - the covariance matrix of x 
! delta_x - the increment in x for the next iteration
! iconverged - flag indicating if convergence was met (yes or no)
! ifail - flag indicating if this process failed (yes or no)
!
! local variables:
! Sx_inv - the inverse of Sx
! delta_x_dir - the unit direction vectors for delta-x 
! delta_x_distance - the total length in x-space of delta_x
! delta_x_constrained - the values of delta_x after being constrained
!-----------------------------------------------------------------------
subroutine OPTIMAL_ESTIMATION(iter,iter_Max,nx,ny, &
                              conv_crit,delta_x_Max, &
                              y,f,x,x_Ap,K,Sy,Sa_inv, &
                              Sx,delta_x,iconverged,ifail, &
                              idiag_output)

  integer, intent(in):: iter
  integer, intent(in):: iter_Max
  integer, intent(in):: idiag_output
  integer, intent(in):: ny
  integer, intent(in):: nx
  real(kind=real4), intent(in):: conv_crit
  real(kind=real4), dimension(:), intent(in):: delta_x_Max
  real(kind=real4), dimension(:), intent(in):: y
  real(kind=real4), dimension(:), intent(in):: f
  real(kind=real4), dimension(:), intent(in):: x
  real(kind=real4), dimension(:), intent(in):: x_Ap
  real(kind=real4), dimension(:,:), intent(in):: K
  real(kind=real4), dimension(:,:), intent(in):: Sy
  real(kind=real4), dimension(:,:), intent(in):: Sa_inv
  real(kind=real4), dimension(:,:), intent(out):: Sx
  real(kind=real4), dimension(:), intent(out):: delta_x
  real(kind=real4), dimension(ny,ny):: Sy_inv
  real(kind=real4), dimension(nx,nx):: Sx_inv
  real(kind=real4), dimension(nx):: delta_x_dir
  real(kind=real4), dimension(nx):: delta_x_constrained
  real(kind=real4):: delta_x_distance_constrained
  real(kind=real4):: delta_x_distance
  integer, intent(out):: ifail
  integer, intent(out):: iconverged
  integer:: isingular
  real:: conv_test
  integer:: ix

  iconverged = sym%NO
  delta_x = Missing_Value_Real4
  Sx = Missing_Value_Real4

  if (ny == 2) then
   call INVERT_2x2(Sy,Sy_inv,isingular)
  else
   call INVERT_3x3(Sy,Sy_inv,isingular)
  endif
  if (isingular == sym%YES) then
   print *, "Cloud Height warning ==> Singular Sy in ACHA "
   ifail = sym%YES
   print *, "y = ", y
   print *, "Sy = ", Sy
  !stop
   return 
  endif

  !---- compute next step
  Sx_inv = Sa_inv + matmul(Transpose(K),matmul(Sy_inv,K)) !(Eq.102 Rodgers)
  call INVERT_3x3(Sx_inv,Sx,isingular)
  if (isingular == sym%YES) then
   print *, "Cloud Height warning ==> Singular Sx in ACHA "
   ifail = sym%YES
   return
  endif
  
  delta_x = matmul(Sx,(matmul(Transpose(K),matmul(Sy_inv,(y-f))) +  &
                       matmul(Sa_inv,x_Ap-x) ))

  !--------------------------------------------------------------
  ! check for convergence
  !--------------------------------------------------------------

  !--- compute convergence metric
  conv_test = abs(sum(delta_x*matmul(Sx_inv,delta_x)))

  !--- control step size  (note change to preserve direction)
  if (idiag_output == sym%YES) then
          print *, "Sa_inv = ", Sa_inv
          print *, "Sy_inv = ", Sy_inv
          print *, "Sx_inv = ", Sx_inv
          print *, "Sx = ", Sx
          print *, "delta_x = ", delta_x
  endif
  delta_x_distance = sqrt(sum(delta_x**2))
  if (delta_x_distance > 0.0) then
     DO ix = 1,nx
        delta_x_dir(ix) = delta_x(ix) / delta_x_distance
     ENDDO
     DO ix = 1,nx
        delta_x_constrained(ix) =  &
             sign(min(delta_x_Max(ix),abs(delta_x(ix))) , delta_x(ix) )
     END DO
     delta_x_distance_constrained = sqrt(sum(delta_x_constrained**2))
     DO ix = 1,nx
        delta_x(ix) = delta_x_dir(ix)*delta_x_distance_constrained
     END DO
  endif
! print *, "constrained delta_x = ", delta_x

  !--- check for non-traditional convergence
  if ((abs(delta_x(1)) < 0.1) .and. (iter > 1)) then
    iconverged = sym%YES
    ifail = sym%NO
  endif

  !--- check for traditional convergence
  if (conv_test < conv_crit) then
      iconverged = sym%YES
      ifail = sym%NO
  endif

  !--- check for exceeding allowed number of interactions
  if (iter > iter_Max) then
      iconverged = sym%NO
      ifail = sym%YES
  endif

  end subroutine OPTIMAL_ESTIMATION

!----------------------------------------------------------------------
  subroutine COMPUTE_Sy(Acha_Mode_Flag, &
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


 end subroutine COMPUTE_Sy

 !---------------------------------------------------------------------
 subroutine COMPUTE_FORWARD_MODEL_AND_KERNEL(Acha_Mode_Flag, x,         &
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

  integer(kind=int4), intent(in):: Acha_Mode_Flag
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
  Bc_67um = PLANCK_RAD_FAST( Chan_Idx_67um, Tc, dB_dT = dB_dTc_67um)
  Bc_85um = PLANCK_RAD_FAST( Chan_Idx_85um, Tc, dB_dT = dB_dTc_85um)
  Bc_11um = PLANCK_RAD_FAST( Chan_Idx_11um, Tc, dB_dT = dB_dTc_11um)
  Bc_12um = PLANCK_RAD_FAST( Chan_Idx_12um, Tc, dB_dT = dB_dTc_12um)
  Bc_133um = PLANCK_RAD_FAST(Chan_Idx_133um, Tc, dB_dT = dB_dTc_133um)

  if (idiag_output == sym%YES) then
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

  if (idiag_output == sym%YES) then
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

 !--- determine number of channels
  SELECT CASE(Acha_Mode_Flag)
     CASE(1)  !avhrr, goes-im
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_12um
     CASE(2)  !goes-nop
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_133um
     CASE(3)  !goes-r
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_12um
       Emiss_Vector(3) = Emiss_133um
     CASE(4)  !viirs
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_12um
       Emiss_Vector(3) = Emiss_85um
     CASE(5)  !goes-im 3 chan
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_67um
       Emiss_Vector(3) = Emiss_12um
     CASE(6)  !goes-np 3 chan
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_67um
       Emiss_Vector(3) = Emiss_133um
  END SELECT

 if (idiag_output == sym%YES) then
         print *, "f = ", f
 endif

end subroutine COMPUTE_FORWARD_MODEL_AND_KERNEL

!----------------------------------------------------------------------
subroutine COMPUTE_APRIORI_BASED_ON_TYPE(Cloud_Type, &
                           Ttropo,T11um,mu,Tc_Opaque, &
                           Tc_Ap, Tc_Ap_Uncer, &
                           Ec_Ap, Ec_Ap_Uncer, &
                           Beta_Ap, Beta_Ap_Uncer)

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

  if (Cloud_Type == sym%FOG_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_fog_type  
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = 2.0*Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == sym%WATER_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_Water_type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == sym%SUPERCOOLED_TYPE .or. &
          Cloud_Type == sym%MIXED_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_mixed_type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == sym%TICE_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_Opaque_Ice_type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == sym%OVERSHOOTING_TYPE) then
     Tc_Ap = Tc_Ap_Opaque
     Tau_Ap = Tau_Ap_Opaque_Ice_type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
     Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
  elseif (Cloud_Type == sym%CIRRUS_TYPE) then
     Tc_Ap = Tc_Ap_Cirrus
     Tau_Ap = Tau_Ap_Cirrus_type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus
     Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus
  elseif (Cloud_Type == sym%OVERLAP_TYPE) then
     Tc_Ap = Tc_Ap_Cirrus
     Tau_Ap = Tau_Ap_overlap_type
     Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus
     Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus
  endif

  !--- determine beta apriori and fit parameters based on 
  !--- phase (derived from type)
          
   !--- water phase clouds
   IF ((Cloud_Type == sym%FOG_TYPE) .or. &
      (Cloud_Type == sym%WATER_TYPE) .or. &
      (Cloud_Type == sym%SUPERCOOLED_TYPE)) THEN
      Beta_Ap = Beta_Ap_Water
      Beta_Ap_Uncer = Beta_Ap_Uncer_Water
   ELSE
   !--- ice phase clouds
      Beta_Ap = Beta_Ap_Ice
      Beta_Ap_Uncer = Beta_Ap_Uncer_Ice
  END IF

  !--- convert Tau_Ap to Emissivity
  Ec_Ap = 1.0 - exp(-Tau_Ap/mu)  !slow!

  end subroutine COMPUTE_APRIORI_BASED_ON_TYPE

!----------------------------------------------------------------------
subroutine COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO( &
                           Cloud_Phase, Emiss_11um_Tropo, &
                           Ttropo,T11um,T11um_Lrc,Tc_Opaque, &
                           Tc_Ap, Tc_Ap_Uncer, &
                           Ec_Ap, Ec_Ap_Uncer, &
                           Beta_Ap, Beta_Ap_Uncer)

  integer, intent(in):: Cloud_Phase
  real(kind=real4), intent(in):: Emiss_11um_Tropo
  real(kind=real4), intent(in):: Ttropo
  real(kind=real4), intent(in):: T11um
  real(kind=real4), intent(in):: T11um_Lrc
  real(kind=real4), intent(in):: Tc_Opaque
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

  if (Cloud_Phase /= sym%ICE_PHASE) then
    Tc_Ap = Tc_Ap_Opaque
    Tc_Ap_Uncer = Tc_Ap_Uncer_Opaque
    Ec_Ap = 0.90          !------ parameter needed
    Ec_Ap_Uncer = Ec_Ap_Uncer_Opaque
    Beta_Ap = Beta_Ap_Water
    Beta_Ap_Uncer = Beta_Ap_Uncer_Water
  endif

  if (Cloud_Phase == sym%ICE_PHASE) then

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

!   Tc_Ap = Tc_Ap_Cirrus
!   Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus

    Ec_Ap = min(0.99,max(0.1,Emiss_Weight)) 
    Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus

    Beta_Ap = Beta_Ap_Ice
    Beta_Ap_Uncer = Beta_Ap_Uncer_Ice

!   print *, "ice apriori ", Tc_Ap, Emiss_11um_Tropo, Ec_Ap, Beta_Ap
!   print *, " --- ", Tc_Ap_Uncer, Ec_Ap_Uncer, Beta_Ap_Uncer

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
subroutine DETERMINE_SFC_TYPE_FORWARD_MODEL(Surface_Type, &
                                         Snow_Class, &
                                         Latitude, &
                                         Ch20_Surface_Emissivity, &
                                         Sfc_Type_Forward_Model)

integer(kind=int1), intent(in):: Surface_Type
integer(kind=int1), intent(in):: Snow_Class
real(kind=real4), intent(in):: Latitude
real(kind=real4), intent(in):: Ch20_Surface_Emissivity
integer(kind=int4), intent(out):: Sfc_Type_Forward_Model


if (Surface_Type == sym%WATER_SFC) then
        Sfc_Type_Forward_Model = 0
else
        Sfc_Type_Forward_Model = 1   !Land
endif
if (Snow_Class == sym%SNOW .and. &
    Latitude > -60.0) then
        Sfc_Type_Forward_Model = 2   !Snow
endif
if (Surface_Type /= sym%WATER_SFC .and. &
    Snow_Class == sym%NO_SNOW .and.  &
    Ch20_Surface_Emissivity > 0.90 .and.  &
    abs(Latitude) < 60.0) then
        Sfc_Type_Forward_Model = 3   !Desert
endif
if (Snow_Class == sym%SEA_ICE .and. &
    Latitude > -60.0) then
        Sfc_Type_Forward_Model = 4   !Arctic
endif
if (Snow_Class /= sym%NO_SNOW .and. Latitude < -60.0) then
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
!----------------------------------------------------------------------
subroutine Compute_Sy_Based_On_Clear_Sky_Covariance_Andy(   &
                                             Emiss_Vector,  &
                                             Acha_Mode_Flag, &
                                             y_variance, &
                                             Sy) 


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
 if (Acha_Mode_Flag >= 3) then
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
    Sy(1,3) = Trans2*Bt_11um_Btd_11um_85um_Covar
    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_67um_Btd_11um_67um_Covar
    Sy(2,3) = Trans2*Btd_11um_67um_Btd_11um_133um_Covar
    Sy(3,3) = T11um_12um_Cal_Uncer**2 + Sub_Pixel_Uncer(3) + Trans2*Btd_11um_12um_Btd_11um_12um_Covar
    Sy(3,1) = Sy(1,3)
    Sy(3,2) = Sy(2,3)
 endif
 if (Acha_Mode_Flag == 6) then
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_12um_Covar
    Sy(1,3) = Trans2*Bt_11um_Btd_11um_85um_Covar
    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_67um_Btd_11um_67um_Covar
    Sy(2,3) = Trans2*Btd_11um_67um_Btd_11um_133um_Covar
    Sy(3,3) = T11um_133um_Cal_Uncer**2 + Sub_Pixel_Uncer(3) + Trans2*Btd_11um_133um_Btd_11um_133um_Covar
    Sy(3,1) = Sy(1,3)
    Sy(3,2) = Sy(2,3)
 endif


end subroutine Compute_Sy_Based_On_Clear_Sky_Covariance_Andy

!----------------------------------------------------------------------
! Compute Sy based on the clear-sky constant error covariance calculations.
! Using Chang-Hwan's fuller expression
!
!
! This assumes that 
! Acha_Mode_Flag: 0=11um,1=11+12um,2=11+13.3um,3=11+12+13.3um,4=8.5+11+12um
!
! Input:
! Emiss_Vector = a vector of emissivities in each channel. 
! Acha_Mode_Flag: 0=11um,1=11+12um,2=11+13.3um,3=11+12+13.3um,4=8.5+11+12um
! Sfc_Type_Forward_Model = the surface type used for covariance calcs 
! y_variance = the variance computed a 3x3 array for each element of y
!
! Output:
! Sy = Constant error covariance matrix
!----------------------------------------------------------------------
subroutine Compute_Sy_Based_On_Clear_Sky_Constant_Covariance_Chang( &
                                             Emiss_Vector,  &
                                             Acha_Mode_Flag, &
                                             y_variance, &
                                             Sy) 


 real(kind=real4), intent(in), dimension(:):: Emiss_Vector
 integer(kind=int4), intent(in):: Acha_Mode_Flag
 real(kind=real4), intent(in), dimension(:):: y_variance
 real(kind=real4), intent(out), dimension(:,:):: Sy
 real(kind=real4):: Emiss_67um
 real(kind=real4):: Emiss_11um
 real(kind=real4):: Emiss_12um
 real(kind=real4):: Emiss_133um
 real(kind=real4):: Emiss_85um
 real(kind=real4):: Inst_Noise_variance_11um, &
                    Inst_Noise_variance_12um, &
                    Inst_Noise_variance_133um,&
                    Inst_Noise_variance_85um, &
                    Inst_Noise_variance_67um

Inst_Noise_variance_11um = (1.0-Bt_11um_Bt_12um_Covar/ &
                       sqrt(Bt_11um_Bt_11um_Covar*Bt_12um_Bt_12um_Covar))*Bt_11um_Bt_11um_Covar
Inst_Noise_variance_12um = (1.0-Bt_11um_Bt_12um_Covar/ &
                      sqrt(Bt_11um_Bt_11um_Covar*Bt_12um_Bt_12um_Covar))*Bt_12um_Bt_12um_Covar
Inst_Noise_variance_133um = (1.0-Bt_11um_Bt_133um_Covar/ &
                       sqrt(Bt_11um_Bt_11um_Covar*Bt_133um_Bt_133um_Covar))*Bt_133um_Bt_133um_Covar
Inst_Noise_variance_85um = (1.0-Bt_11um_Bt_85um_Covar/ &
                      sqrt(Bt_11um_Bt_11um_Covar*Bt_85um_Bt_85um_Covar))*Bt_85um_Bt_85um_Covar
Inst_Noise_variance_67um = (1.0-Bt_11um_Bt_67um_Covar/ &
                      sqrt(Bt_11um_Bt_11um_Covar*Bt_67um_Bt_67um_Covar))*Bt_67um_Bt_67um_Covar

 if (Acha_Mode_Flag == 0) then !NO MODE 0. Correct this - WCS3

    Emiss_11um = Emiss_Vector(1)
    Sy(1,1) = Inst_Noise_variance_11um  + y_variance(1) + &
              (1.0-Emiss_11um)**2 *(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)
 endif

 if (Acha_Mode_Flag == 1) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_12um = Emiss_Vector(2)

    Sy(1,1) = Inst_Noise_variance_11um + y_variance(1) + &
              ((1.0-Emiss_11um)**2) *(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_12um + y_variance(2) + &
              ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) + &
              ((1.0-Emiss_12um)**2)*(Bt_12um_Bt_12um_Covar-Inst_Noise_variance_12um) - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*(2.0*Bt_11um_Bt_12um_Covar)

    Sy(2,1) = Inst_Noise_variance_11um + &
              ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*Bt_11um_Bt_12um_Covar

    Sy(1,2) = Sy(2,1)
 endif

 if (Acha_Mode_Flag == 2) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_133um = Emiss_Vector(2)

    Sy(1,1) = Inst_Noise_variance_11um + y_variance(1) + &
              ((1.0-Emiss_11um)**2) * (Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_133um + y_variance(2) + &
              ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) + &
              ((1.0-Emiss_133um)**2)*(Bt_133um_Bt_133um_Covar-Inst_Noise_variance_133um) - &  
              (1.0-Emiss_11um)*(1.0-Emiss_133um) * (2.0*Bt_11um_Bt_133um_Covar)

    Sy(2,1) = Inst_Noise_variance_11um + &
	       ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) - &
               (1.0-Emiss_11um)*(1.0-Emiss_133um)*Bt_11um_Bt_133um_Covar

    Sy(1,2) = Sy(2,1)
 endif

 if (Acha_Mode_Flag == 3) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_12um = Emiss_Vector(2)
    Emiss_133um = Emiss_Vector(3)

    Sy(1,1) = Inst_Noise_variance_11um  + y_variance(1) + &
	         (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_12um  + y_variance(2) + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) + &
              (1.0-Emiss_12um)**2*(Bt_12um_Bt_12um_Covar-Inst_Noise_variance_12um) - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*2.0*(Bt_11um_Bt_12um_Covar)

    Sy(3,3) = Inst_Noise_variance_11um + Inst_Noise_variance_133um + y_variance(3) + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) + &
              (1.0-Emiss_133um)**2*(Bt_133um_Bt_133um_Covar-Inst_Noise_variance_133um) - &  
              (1.0-Emiss_11um)*(1.0-Emiss_133um)*2.0*(Bt_11um_Bt_133um_Covar)

    Sy(2,1) = Inst_Noise_variance_11um + &
	          (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*(Bt_11um_Bt_12um_Covar)

    Sy(3,1) = Inst_Noise_variance_11um + &
	           (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) - &
               (1.0-Emiss_11um)*(1.0-Emiss_133um)*(Bt_11um_Bt_133um_Covar)

    Sy(2,3) = Inst_Noise_variance_11um + & 
	    ((1.0-Emiss_11um)**2)* (Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) - &
            (1.0-Emiss_11um)*(1.0-Emiss_12um)* (Bt_11um_Bt_12um_Covar)  - &
            (1.0-Emiss_11um)*(1.0-Emiss_133um)* (Bt_11um_Bt_133um_Covar)  +  &
            (1.0-Emiss_12um)*(1.0-Emiss_133um)*(Bt_12um_Bt_133um_Covar) 

    Sy(3,2) = Sy(2,3)
    Sy(1,2) = Sy(2,1)
    Sy(1,3) = Sy(3,1)
 endif

 if (Acha_Mode_Flag == 4) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_12um = Emiss_Vector(2)
    Emiss_85um = Emiss_Vector(3)

    Sy(1,1) = Inst_Noise_variance_11um  + y_variance(1) + &
	       ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) 

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_12um + y_variance(2) + &
              ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  + &
              ((1.0-Emiss_12um)**2)*(Bt_12um_Bt_12um_Covar-Inst_Noise_variance_12um)  - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*2.0*(Bt_11um_Bt_12um_Covar) 

    Sy(3,3) = Inst_Noise_variance_11um + Inst_Noise_variance_85um + y_variance(3) + &
              ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  + &
              ((1.0-Emiss_85um)**2)*(Bt_85um_Bt_85um_Covar-Inst_Noise_variance_85um)  - &  
              (1.0-Emiss_11um)*(1.0-Emiss_85um)*2.0*(Bt_11um_Bt_85um_Covar) 

    Sy(2,1) = Inst_Noise_variance_11um + &
	      ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*(Bt_11um_Bt_12um_Covar) 

    Sy(3,1) = Inst_Noise_variance_11um + &
	       ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  - &
               (1.0-Emiss_11um)*(1.0-Emiss_85um)*(Bt_11um_Bt_85um_Covar) 

    Sy(2,3) = Inst_Noise_variance_11um + & 
	    ((1.0-Emiss_11um)**2) * (Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) - &
            (1.0-Emiss_11um)*(1.0-Emiss_12um)* (Bt_11um_Bt_12um_Covar)  - &
            (1.0-Emiss_11um)*(1.0-Emiss_85um)* (Bt_11um_Bt_85um_Covar) +  &
            (1.0-Emiss_12um)*(1.0-Emiss_85um)*(Bt_12um_Bt_85um_Covar) 

    Sy(3,2) = Sy(2,3)
    Sy(1,2) = Sy(2,1)
    Sy(1,3) = Sy(3,1)
 endif

 if (Acha_Mode_Flag == 5) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_12um = Emiss_Vector(2)
    Emiss_67um = Emiss_Vector(3)

    Sy(1,1) = Inst_Noise_variance_11um  + y_variance(1) + &
	         (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) 

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_12um  + y_variance(2) + &
              ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  + &
              ((1.0-Emiss_12um)**2)*(Bt_12um_Bt_12um_Covar-Inst_Noise_variance_12um)  - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*2.0*(Bt_11um_Bt_12um_Covar) 

    Sy(3,3) = Inst_Noise_variance_11um + Inst_Noise_variance_67um  + y_variance(3) + &
              ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  + &
              ((1.0-Emiss_67um)**2)*(Bt_67um_Bt_67um_Covar-Inst_Noise_variance_67um)  - &  
              (1.0-Emiss_11um)*(1.0-Emiss_67um)*2.0*(Bt_11um_Bt_67um_Covar) 

    Sy(2,1) = Inst_Noise_variance_11um + &
	      ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*(Bt_11um_Bt_12um_Covar) 

    Sy(3,1) = Inst_Noise_variance_11um + &
	       ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  - &
               (1.0-Emiss_11um)*(1.0-Emiss_67um)*(Bt_11um_Bt_67um_Covar) 

    Sy(2,3) = Inst_Noise_variance_11um + & 
	    ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  - &
            (1.0-Emiss_11um)*(1.0-Emiss_12um)* (Bt_11um_Bt_12um_Covar)  - &
            (1.0-Emiss_11um)*(1.0-Emiss_67um)* (Bt_11um_Bt_67um_Covar) +  &
            (1.0-Emiss_12um)*(1.0-Emiss_67um)*(Bt_12um_Bt_67um_Covar) 

    Sy(3,2) = Sy(2,3)
    Sy(1,2) = Sy(2,1)
    Sy(1,3) = Sy(3,1)
 endif

 if (Acha_Mode_Flag == 6) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_133um = Emiss_Vector(2)
    Emiss_67um = Emiss_Vector(3)

    Sy(1,1) = Inst_Noise_variance_11um  + y_variance(1) + &
	      ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) 

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_133um  + y_variance(2) + &
              ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) + &
              ((1.0-Emiss_133um)**2)*(Bt_133um_Bt_133um_Covar-Inst_Noise_variance_133um)  - &
              (1.0-Emiss_11um)*(1.0-Emiss_133um)*2.0*(Bt_11um_Bt_133um_Covar) 

    Sy(3,3) = Inst_Noise_variance_11um + Inst_Noise_variance_67um   + y_variance(3) + &
              ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  + &
              ((1.0-Emiss_67um)**2)*(Bt_67um_Bt_67um_Covar-Inst_Noise_variance_67um)  - &  
              (1.0-Emiss_11um)*(1.0-Emiss_67um)*2.0*(Bt_11um_Bt_67um_Covar) 

    Sy(2,1) = Inst_Noise_variance_11um + &
	      ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)  - &
              (1.0-Emiss_11um)*(1.0-Emiss_133um)*(Bt_11um_Bt_133um_Covar) 

    Sy(3,1) = Inst_Noise_variance_11um + &
	       ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) - &
               (1.0-Emiss_11um)*(1.0-Emiss_67um)*(Bt_11um_Bt_67um_Covar) 

    Sy(2,3) = Inst_Noise_variance_11um + & 
	        ((1.0-Emiss_11um)**2)*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um) - &
                (1.0-Emiss_11um)*(1.0-Emiss_133um)* (Bt_11um_Bt_133um_Covar) - &
                (1.0-Emiss_11um)*(1.0-Emiss_67um)* (Bt_11um_Bt_67um_Covar) +  &
                (1.0-Emiss_133um)*(1.0-Emiss_67um)*(Bt_67um_Bt_133um_Covar) 

    Sy(3,2) = Sy(2,3)
    Sy(1,2) = Sy(2,1)
    Sy(1,3) = Sy(3,1)
 endif

end subroutine Compute_Sy_Based_On_Clear_Sky_Constant_Covariance_Chang


!----------------------------------------------------------------------
! Compute Sy based on the clear-sky variational error covariance calculations.
! Using Chang-Hwan's fuller expression with constant errors
!
! This assumes that 
! Acha_Mode_Flag: 0=11um,1=11+12um,2=11+13.3um,3=11+12+13.3um,4=8.5+11+12um
!
! Input:
! Emiss_Vector = a vector of emissivities in each channel. 
! Acha_Mode_Flag: 0=11um,1=11+12um,2=11+13.3um,3=11+12+13.3um,4=8.5+11+12um
! Sfc_Type_Forward_Model = the surface type used for covariance calcs 
! y_variance = the variance computed a 3x3 array for each element of y
!
! Output:
! Sy = Variational error covariance matrix
!----------------------------------------------------------------------
subroutine Compute_Sy_Based_On_Clear_Sky_Variational_Covariance_Chang( &
                                             Emiss_Vector,  &
                                             Acha_Mode_Flag, &
                                             y_variance, &
                                             BT_CL_67um, &
                                             BT_CL_85um, &
                                             BT_CL_11um, &
                                             BT_CL_12um, &
                                             BT_CL_133um, &
                                             Sy)

 real(kind=real4), intent(in), dimension(:):: Emiss_Vector
 integer(kind=int4), intent(in):: Acha_Mode_Flag
 real(kind=real4), intent(in), dimension(:):: y_variance
 real(kind=real4), intent(in):: BT_CL_67um, BT_CL_85um,BT_CL_11um,BT_CL_12um,BT_CL_133um
 real(kind=real4), intent(out), dimension(:,:):: Sy
 real(kind=real4):: Emiss_67um
 real(kind=real4):: Emiss_11um
 real(kind=real4):: Emiss_12um
 real(kind=real4):: Emiss_133um
 real(kind=real4):: Emiss_85um
 real(kind=real4):: MOD_67um, &
                    MOD_85um, &
                    MOD_11um, &
                    MOD_12um, &
	            MOD_133um
 real(kind=real4):: MOD_67um_Mean, &
                    MOD_85um_Mean, &
	            MOD_11um_Mean, &
		    MOD_12um_Mean, &
		    MOD_133um_Mean
 real(kind=real4):: Inst_Noise_variance_11um, &
                    Inst_Noise_variance_12um, &
		    Inst_Noise_variance_133um,&
		    Inst_Noise_variance_85um, &
		    Inst_Noise_variance_67um

MOD_67um_Mean = Bt_67um_Mean - 200.0
MOD_85um_Mean = Bt_85um_Mean - 245.0
MOD_11um_Mean = Bt_11um_Mean - 245.0
MOD_12um_Mean = Bt_12um_Mean - 245.0
MOD_133um_Mean = Bt_133um_Mean - 223.0

MOD_67um = BT_CL_67um -200.0
MOD_85um = BT_CL_85um -245.0
MOD_11um = BT_CL_11um  - 245.0
MOD_12um = BT_CL_12um  - 245.0
MOD_133um = BT_CL_133um  - 223.0

Inst_Noise_variance_11um = (1.-Bt_11um_Bt_12um_Covar/ &
                  sqrt(Bt_11um_Bt_11um_Covar*Bt_12um_Bt_12um_Covar))*Bt_11um_Bt_11um_Covar
Inst_Noise_variance_12um = (1.-Bt_11um_Bt_12um_Covar/ &
                  sqrt(Bt_11um_Bt_11um_Covar*Bt_12um_Bt_12um_Covar))*Bt_12um_Bt_12um_Covar
Inst_Noise_variance_133um = (1.-Bt_11um_Bt_133um_Covar/ &
                  sqrt(Bt_11um_Bt_11um_Covar*Bt_133um_Bt_133um_Covar))*Bt_133um_Bt_133um_Covar
Inst_Noise_variance_85um = (1.-Bt_11um_Bt_85um_Covar/ &
                  sqrt(Bt_11um_Bt_11um_Covar*Bt_85um_Bt_85um_Covar))*Bt_85um_Bt_85um_Covar
Inst_Noise_variance_67um = (1.-Bt_11um_Bt_67um_Covar/ &
                  sqrt(Bt_11um_Bt_11um_Covar*Bt_67um_Bt_67um_Covar))*Bt_67um_Bt_67um_Covar

 if (Acha_Mode_Flag == 0) then !NO MODE 0. Correct this - WCS3

    Emiss_11um = Emiss_Vector(1)
    Sy(1,1) = Inst_Noise_variance_11um  +  &
              (1.0-Emiss_11um)**2 *(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2

 endif

 if (Acha_Mode_Flag == 1) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_12um = Emiss_Vector(2)

    Sy(1,1) = Inst_Noise_variance_11um +  y_variance(1) + &
              (1.0-Emiss_11um)**2 *(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_12um + y_variance(2) + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 + &
              (1.0-Emiss_12um)**2*(Bt_12um_Bt_12um_Covar-Inst_Noise_variance_12um)*(MOD_12um/MOD_12um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*(2.0*Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean)

    Sy(2,1) = Inst_Noise_variance_11um + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*Bt_11um_Bt_12um_Covar*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean)

    Sy(1,2) = Sy(2,1)
 endif

 if (Acha_Mode_Flag == 2) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_133um = Emiss_Vector(2)

    Sy(1,1) = Inst_Noise_variance_11um +  y_variance(1) + &
              (1.0-Emiss_11um)**2 * (Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_133um + y_variance(2) +  &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2+ &
              (1.0-Emiss_133um)**2*(Bt_133um_Bt_133um_Covar-Inst_Noise_variance_133um)*(MOD_133um/MOD_133um_Mean)**2 - &  
              (1.0-Emiss_11um)*(1.0-Emiss_133um)* (2.0*Bt_11um_Bt_133um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_133um/MOD_133um_Mean)

    Sy(2,1) = Inst_Noise_variance_11um + &
               (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
               (1.0-Emiss_11um)*(1.0-Emiss_133um)*Bt_11um_Bt_133um_Covar*(MOD_11um/MOD_11um_Mean)*(MOD_133um/MOD_133um_Mean)

    Sy(1,2) = Sy(2,1)
 endif

 if (Acha_Mode_Flag == 3) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_12um = Emiss_Vector(2)
    Emiss_133um = Emiss_Vector(3)

    Sy(1,1) = Inst_Noise_variance_11um  + y_variance(1) +  &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_12um + y_variance(2) + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 + &
              (1.0-Emiss_12um)**2*(Bt_12um_Bt_12um_Covar-Inst_Noise_variance_12um)*(MOD_12um/MOD_12um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*2.0*(Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean)

    Sy(3,3) = Inst_Noise_variance_11um + Inst_Noise_variance_133um + y_variance(3) +  &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 + &
              (1.0-Emiss_133um)**2*(Bt_133um_Bt_133um_Covar-Inst_Noise_variance_133um)*(MOD_133um/MOD_133um_Mean)**2 - &  
              (1.0-Emiss_11um)*(1.0-Emiss_133um)*2.0*(Bt_11um_Bt_133um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_133um/MOD_133um_Mean)

    Sy(2,1) = Inst_Noise_variance_11um + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*(Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean)

    Sy(3,1) = Inst_Noise_variance_11um + &
               (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
               (1.0-Emiss_11um)*(1.0-Emiss_133um)*(Bt_11um_Bt_133um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_133um/MOD_133um_Mean)

    Sy(2,3) = Inst_Noise_variance_11um + & 
            (1.0-Emiss_11um)**2* (Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
            (1.0-Emiss_11um)*(1.0-Emiss_12um)* (Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean) - &
            (1.0-Emiss_11um)*(1.0-Emiss_133um)* (Bt_11um_Bt_133um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_133um/MOD_133um_Mean)+  &
            (1.0-Emiss_12um)*(1.0-Emiss_133um)*(Bt_12um_Bt_133um_Covar)*(MOD_12um/MOD_12um_Mean)*(MOD_133um/MOD_133um_Mean)

!Sy(2,3) = 0
!Sy(2,1) = 0
!Sy(3,1) = 0

    Sy(3,2) = Sy(2,3)
    Sy(1,2) = Sy(2,1)
    Sy(1,3) = Sy(3,1)
 endif

 if (Acha_Mode_Flag == 4) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_12um = Emiss_Vector(2)
    Emiss_85um = Emiss_Vector(3)

    Sy(1,1) = Inst_Noise_variance_11um  + y_variance(1) + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_12um + y_variance(2) +  &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 + &
              (1.0-Emiss_12um)**2*(Bt_12um_Bt_12um_Covar-Inst_Noise_variance_12um)*(MOD_12um/MOD_12um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*2.0*(Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean)

    Sy(3,3) = Inst_Noise_variance_11um + Inst_Noise_variance_85um + y_variance(3) + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 + &
              (1.0-Emiss_85um)**2*(Bt_85um_Bt_85um_Covar-Inst_Noise_variance_85um)*(MOD_85um/MOD_85um_Mean)**2 - &  
              (1.0-Emiss_11um)*(1.0-Emiss_85um)*2.0*(Bt_11um_Bt_85um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_85um/MOD_85um_Mean)

    Sy(2,1) = Inst_Noise_variance_11um + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*(Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean)

    Sy(3,1) = Inst_Noise_variance_11um + &
               (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
               (1.0-Emiss_11um)*(1.0-Emiss_85um)*(Bt_11um_Bt_85um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_85um/MOD_85um_Mean)

    Sy(2,3) = Inst_Noise_variance_11um + & 
            (1.0-Emiss_11um)**2* (Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
            (1.0-Emiss_11um)*(1.0-Emiss_12um)* (Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean) - &
            (1.0-Emiss_11um)*(1.0-Emiss_85um)* (Bt_11um_Bt_85um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_85um/MOD_85um_Mean)+  &
            (1.0-Emiss_12um)*(1.0-Emiss_85um)*(Bt_12um_Bt_85um_Covar)*(MOD_12um/MOD_12um_Mean)*(MOD_85um/MOD_85um_Mean)

    Sy(3,2) = Sy(2,3)
    Sy(1,2) = Sy(2,1)
    Sy(1,3) = Sy(3,1)
 endif

 if (Acha_Mode_Flag == 5) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_12um = Emiss_Vector(2)
    Emiss_67um = Emiss_Vector(3)

    Sy(1,1) = Inst_Noise_variance_11um + y_variance(1) +  &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_12um + y_variance(2) + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 + &
              (1.0-Emiss_12um)**2*(Bt_12um_Bt_12um_Covar-Inst_Noise_variance_12um)*(MOD_12um/MOD_12um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*2.0*(Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean)

    Sy(3,3) = Inst_Noise_variance_11um + Inst_Noise_variance_67um + y_variance(3) +  &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 + &
              (1.0-Emiss_67um)**2*(Bt_67um_Bt_67um_Covar-Inst_Noise_variance_67um)*(MOD_67um/MOD_67um_Mean)**2 - &  
              (1.0-Emiss_11um)*(1.0-Emiss_67um)*2.0*(Bt_11um_Bt_67um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_67um/MOD_67um_Mean)

    Sy(2,1) = Inst_Noise_variance_11um + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_12um)*(Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean)

    Sy(3,1) = Inst_Noise_variance_11um + &
               (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
               (1.0-Emiss_11um)*(1.0-Emiss_67um)*(Bt_11um_Bt_67um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_67um/MOD_67um_Mean)

    Sy(2,3) = Inst_Noise_variance_11um + & 
            (1.0-Emiss_11um)**2* (Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
            (1.0-Emiss_11um)*(1.0-Emiss_12um)* (Bt_11um_Bt_12um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_12um/MOD_12um_Mean) - &
            (1.0-Emiss_11um)*(1.0-Emiss_67um)* (Bt_11um_Bt_67um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_67um/MOD_67um_Mean)+  &
            (1.0-Emiss_12um)*(1.0-Emiss_67um)*(Bt_12um_Bt_67um_Covar)*(MOD_12um/MOD_12um_Mean)*(MOD_67um/MOD_67um_Mean)

    Sy(3,2) = Sy(2,3)
    Sy(1,2) = Sy(2,1)
    Sy(1,3) = Sy(3,1)
 endif

 if (Acha_Mode_Flag == 6) then

    Emiss_11um = Emiss_Vector(1)
    Emiss_133um = Emiss_Vector(2)
    Emiss_67um = Emiss_Vector(3)

    Sy(1,1) = Inst_Noise_variance_11um  + y_variance(1) +  &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2

    Sy(2,2) = Inst_Noise_variance_11um + Inst_Noise_variance_133um + y_variance(2) + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 + &
              (1.0-Emiss_133um)**2*(Bt_133um_Bt_133um_Covar-Inst_Noise_variance_133um)*(MOD_133um/MOD_133um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_133um)*2.0*(Bt_11um_Bt_133um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_133um/MOD_133um_Mean)

    Sy(3,3) = Inst_Noise_variance_11um + Inst_Noise_variance_67um + y_variance(3) +  &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 + &
              (1.0-Emiss_67um)**2*(Bt_67um_Bt_67um_Covar-Inst_Noise_variance_67um)*(MOD_67um/MOD_67um_Mean)**2 - &  
              (1.0-Emiss_11um)*(1.0-Emiss_67um)*2.0*(Bt_11um_Bt_67um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_67um/MOD_67um_Mean)

    Sy(2,1) = Inst_Noise_variance_11um + &
              (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
              (1.0-Emiss_11um)*(1.0-Emiss_133um)*(Bt_11um_Bt_133um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_133um/MOD_133um_Mean)

    Sy(3,1) = Inst_Noise_variance_11um + &
               (1.0-Emiss_11um)**2*(Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
               (1.0-Emiss_11um)*(1.0-Emiss_67um)*(Bt_11um_Bt_67um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_67um/MOD_67um_Mean)

    Sy(2,3) = Inst_Noise_variance_11um + & 
            (1.0-Emiss_11um)**2* (Bt_11um_Bt_11um_Covar-Inst_Noise_variance_11um)*(MOD_11um/MOD_11um_Mean)**2 - &
            (1.0-Emiss_11um)*(1.0-Emiss_133um)* (Bt_11um_Bt_133um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_133um/MOD_133um_Mean) - &
            (1.0-Emiss_11um)*(1.0-Emiss_67um)* (Bt_11um_Bt_67um_Covar)*(MOD_11um/MOD_11um_Mean)*(MOD_67um/MOD_67um_Mean)+  &
            (1.0-Emiss_133um)*(1.0-Emiss_67um)*(Bt_67um_Bt_133um_Covar)*(MOD_133um/MOD_133um_Mean)*(MOD_67um/MOD_67um_Mean)

    Sy(3,2) = Sy(2,3)
    Sy(1,2) = Sy(2,1)
    Sy(1,3) = Sy(3,1)
 endif

end subroutine Compute_Sy_Based_On_Clear_Sky_Variational_Covariance_Chang

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
subroutine DETERMINE_ACHA_MODE_BASED_ON_CHANNELS(Acha_Mode_Flag,Chan_On)

  integer, intent(INout):: Acha_Mode_Flag
  integer(kind=int1), intent(in), dimension(:):: Chan_On

  if (Acha_Mode_Flag == -1) then
     if (Chan_On(Chan_Idx_11um)==sym%YES .and. Chan_On(Chan_Idx_12um)==sym%YES) then
        if (Chan_On(Chan_Idx_133um)==sym%YES) then
            Acha_Mode_Flag = 3          ! 11/12/13.3 um
        elseif (Chan_On(Chan_Idx_85um)==sym%YES) then
            Acha_Mode_Flag = 4          ! 8.5/11/12 um
        elseif (Chan_On(Chan_Idx_67um)==sym%YES) then
            Acha_Mode_Flag = 5          ! 6.7/11/12 um
        else
            Acha_Mode_Flag = 1          ! 11/12 um
        endif
     endif
  endif
  if (Acha_Mode_Flag == -1) then
     if (Chan_On(Chan_Idx_12um)==sym%NO) then
        if (Chan_On(Chan_Idx_67um)==sym%YES .and. Chan_On(Chan_Idx_133um) == sym%YES) then
              Acha_Mode_Flag = 6       !6.7/11/13.3
        endif
        if (Chan_On(Chan_Idx_67um)==sym%NO .and. Chan_On(Chan_Idx_133um) == sym%YES) then
              Acha_Mode_Flag = 2       !11/13.3
        endif
     endif
  endif
end subroutine  DETERMINE_ACHA_MODE_BASED_ON_CHANNELS

!----------------------------------------------------------------------------
! Function intERPOLATE_PROFILE_ACHA
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
 subroutine intERPOLATE_PROFILE_ACHA(z1,z2,z3,z4,lonx,latx,z)

  real, dimension(:), intent(in):: z1
  real, dimension(:), intent(in):: z2
  real, dimension(:), intent(in):: z3
  real, dimension(:), intent(in):: z4
  real, intent(in):: lonx
  real, intent(in):: latx
  real, dimension(:), intent(out):: z

  !--- linear inteprpolation scheme
  z =  (1.0-lonx) * ((1.0-latx) * z1 + (latx)* z3) + &
           (lonx) * ((1.0-latx) * z2 + (latx)* z4)

 end subroutine intERPOLATE_PROFILE_ACHA

!----------------------------------------------------------------------------
! Function intERPOLATE_NWP_ACHA
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
 function intERPOLATE_NWP_ACHA(z1,z2,z3,z4,lonx,latx) result(z)

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

 end function intERPOLATE_NWP_ACHA

 !---------------------------------------------------------------------
 ! find Inversion Level - highest Level Inversion below trop
 !---------------------------------------------------------------------
 function DETERMINE_INVERSION_LEVEL(Tropo_Level, Sfc_Level) result(Inversion_Level)
   integer, intent(in):: Tropo_Level
   integer, intent(in):: Sfc_Level
   integer:: Inversion_Level
   integer:: k

   Inversion_Level = 0
   do k = Tropo_Level, Sfc_Level-1
      if ((Temp_Prof_Rtm(k) - Temp_Prof_Rtm(k+1) > Delta_T_Layer_Inversion) .and. &
          (Press_Prof_Rtm(k) >= P_Inversion_Min)) then
          Inversion_Level = k
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
 subroutine LOCAL_LINEAR_RADIATIVE_CENTER(Meander_Flag, &
                                          Grid_Data, &
                                          Element_Start, Number_Of_Elements, & 
                                          Line_Start, Number_Of_Lines, & 
                                          Max_Grad_Distance, &
                                          Grad_Flag,  &
                                          Missing_LRC_Value, &
                                          Skip_LRC_Mask, &
                                          Min_Grid_Data_Valid, Max_Grid_Data_Valid, &
                                          ielem_LRC, iline_LRC)

  integer, intent(in):: Meander_Flag
  real (kind=real4), intent(in), dimension(:,:) :: Grid_Data
  integer (kind=int4), intent(in):: Element_Start
  integer (kind=int4), intent(in):: Number_of_Elements
  integer (kind=int4), intent(in):: Line_Start
  integer (kind=int4), intent(in):: Number_of_Lines
  integer (kind=int4), intent(in):: Max_Grad_Distance
  integer (kind=int4), intent(in):: Grad_Flag
  integer (kind=int4), intent(in):: Missing_LRC_Value
  integer (kind=int1), intent(in), dimension(:,:):: Skip_LRC_Mask
  real (kind=real4), intent(in):: Min_Grid_Data_Valid
  real (kind=real4), intent(in):: Max_Grid_Data_Valid
  integer (kind=int4), intent(out), dimension(:,:):: ielem_LRC
  integer (kind=int4), intent(out), dimension(:,:):: iline_LRC
  real, dimension(3,3):: Grad_Array
  integer, dimension(2):: Grad_Indices
  integer:: ielem
  integer:: iline
  integer:: ielem_Previous
  integer:: iline_Previous
  integer:: ielem_Next
  integer:: iline_Next
  real:: Grad_Temp
  integer:: Element_End
  integer:: Line_End
  integer:: ipoint
  integer:: ielem_dir
  integer:: iline_dir
 
  Element_End = Number_of_Elements + Element_Start - 1
  Line_End = Number_of_Lines + Line_Start - 1

  !--- initialize
  ielem_LRC = Missing_LRC_Value
  iline_LRC = Missing_LRC_Value

!----------------------------------------------------------------------
! loop through pixels in segment
!----------------------------------------------------------------------
Element_Loop:  do ielem = Element_Start+1, Element_End-1
Line_Loop:    do iline = Line_Start+1, Line_End-1

      !--- skip data due to mask
      if (Skip_LRC_Mask(ielem,iline) == sym%YES) cycle

      !-- check for out of bounds data
      if (Grad_Flag ==  1 .and. Grid_Data(ielem,iline) < Min_Grid_Data_Valid) cycle
      if (Grad_Flag ==  -1 .and. Grid_Data(ielem,iline) > Max_Grid_Data_Valid) cycle

      !-- check for data that already meets LRC criteria
      if ((Grad_Flag ==  1 .and. Grid_Data(ielem,iline) > Max_Grid_Data_Valid) .or. &
          (Grad_Flag ==  -1 .and. Grid_Data(ielem,iline) < Min_Grid_Data_Valid)) then
              ielem_LRC(ielem,iline) = ielem
              iline_LRC(ielem,iline) = iline
      endif

      !--- initialize previous variables
      ielem_Previous = ielem
      iline_Previous = iline

!     print *, "starting gradient path ", ielem, iline

      !---- go long gradient and check for a reversal or saturation
      do ipoint = 1,Max_Grad_Distance

!       print *, "ipoint = ", ipoint, ielem_Previous, iline_Previous

        !--- compute local gradient, find strongest gradient in 3x3 array and compute direction
        if (ipoint == 1 .or. Meander_Flag == sym%YES) then

         !--- construct 3x3 array for analysis
         Grad_Array =  &
           Grid_Data(ielem_Previous-1:ielem_Previous+1,iline_Previous-1:iline_Previous+1) -  &
           Grid_Data(ielem_Previous,iline_Previous)

         !--- look for bad data
         if (minval(Grad_Array) == Missing_Value_Real4) exit 

         !--- compute local gradients, find strongest gradient
         if (Grad_Flag == 1) then
          Grad_Indices = maxloc(Grad_Array)
         else
          Grad_Indices = minloc(Grad_Array)
         endif 

         !--- compute direction
         ielem_Dir = Grad_Indices(1)  - 2
         iline_Dir = Grad_Indices(2)  - 2

!        print *, "new direction = ", ielem_Dir, iline_Dir
         !--- check for pixels that are located at  minima/maxima
         if (ielem_Dir == 0 .and. iline_Dir == 0) then
           ielem_LRC(ielem,iline) = ielem_Previous
           iline_LRC(ielem,iline) = iline_Previous
!          print *, "valley found"
           exit
         endif

        endif

        !-- select next point on the path
        ielem_Next = ielem_Previous + ielem_Dir
        iline_Next = iline_Previous + iline_Dir

        !--- check for hitting segment boundaries
        if (ielem_Next == Element_Start .or. ielem_Next == Element_End .or. &
             iline_Next == Line_Start .or. iline_Next == Line_End) then
              ielem_LRC(ielem,iline) = ielem_Previous
              iline_LRC(ielem,iline) = iline_Previous
!             print *, "hit segment boundary"
              exit
         endif

         !--- check for hitting bad data
         if (Skip_LRC_Mask(ielem_Next,iline_Next) == sym%YES) then
              ielem_LRC(ielem,iline) = ielem_Previous
              iline_LRC(ielem,iline) = iline_Previous
!             print *, "bad data"
              exit
         endif

         !--- check for sign reversal
         if (Meander_Flag == sym%NO) then

          Grad_Temp = Grid_Data(ielem_Next,iline_Next) -  &
                      Grid_Data(ielem_Previous,iline_Previous)

          if (Grad_Flag * Grad_Temp < 0) then
              ielem_LRC(ielem,iline) = ielem_Previous
              iline_LRC(ielem,iline) = iline_Previous
!             print *, "sign reversal"
              exit
          endif
         endif

         !--- check for saturation
         if (Grad_Flag == 1 .and. Grid_Data(ielem_Next,iline_Next) > Max_Grid_Data_Valid) then
              ielem_LRC(ielem,iline) = ielem_Next
              iline_LRC(ielem,iline) = iline_Next
!             print *, "positive saturation"
              exit
         endif
         if (Grad_Flag == -1 .and. Grid_Data(ielem_Next,iline_Next) < Min_Grid_Data_Valid) then
              ielem_LRC(ielem,iline) = ielem_Next
              iline_LRC(ielem,iline) = iline_Next
!             print *, "negative saturation"
              exit
         endif

         !--- store position
         ielem_Previous = ielem_Next
         iline_Previous = iline_Next

      enddo

!   print *, "final point = ", ielem_LRC(ielem,iline), iline_LRC(ielem,iline)

!   if (ielem_LRC(ielem,iline) > 0) then
!   Pix_Array_1(ielem,iline) = Bt_Ch31(ielem_LRC(ielem,iline),iline_LRC(ielem,iline))
!   Pix_Array_2(ielem,iline) = sqrt(real(ielem_LRC(ielem,iline) - ielem)**2 + real(iline_LRC(ielem,iline) - iline)**2)
!   endif

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

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module AWG_CLOUD_HEIGHT
