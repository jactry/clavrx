!$Id$
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
! Modes
! 0 - Use this mode to not call ACHA from the framework
! 1 - 11 um                          0           
! 2 - 11 + 6.7 um                    7
! 3 - 11 + 12 um                     1
! 4 - 11 + 13.3 um                   2
! 5 - 11 + 8.5 + 12 um               4
! 6 - 11 + 6.7 + 12 um               5
! 7 - 11 + 6.7 + 13.3 um             6
! 8 - 11 + 12 + 13.3 um              3
! 9 - 11 + 12 + 13.3-pseudo um       -
!
! MULTI_LAYER_LOGIC_FLAG
! 0 - (baseline) just use the multilayer id in cloud type
! 1 - treat all multilayer like cirrus
! 2 - assume all cirrus are multilayer and let acha decide
!----------------------------------------------------------------------
!Changes needed to get into SAPF
!
! - Renamed AWG_CLOUD_HEIGHT_ALGORITHM to AWG_CLOUD_HEIGHT_ALGORITHM_ACHA
! - Renamed LOCAL_LINEAR_RADIATIVE_CENTER to LOCAL_LINEAR_RADIATIVE_CENTER_ACHA
! - Renamed module from AWG_CLOUD_HEIGHT to AWG_CLOUD_HEIGHT_ACHA
! - Had to redo Skip_LRC_Mask due to issues in Framework
!
! ** Note:  These changes are in the Framework repository only.
!
!----------------------------------------------------------------------
  use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, real8, dtor, &
           acha_output_struct,acha_symbol_struct, &
           acha_input_struct, acha_rtm_nwp_struct, &
           PLANCK_RAD_FAST, PLANCK_TEMP_FAST, &
           INVERT_MATRIX, ACHA_FETCH_PIXEL_NWP_RTM, &
           LOCATE, acha_diag_struct

  use ACHA_MICROPHYSICAL_MODULE

  implicit none

  public:: AWG_CLOUD_HEIGHT_ALGORITHM
  public:: CHECK_ACHA_MODE
  public:: SET_ACHA_VERSION
  public:: LOCAL_LINEAR_RADIATIVE_CENTER

  private:: COMPUTE_LOWER_CLOUD_TEMPERATURE
  private:: KNOWING_P_COMPUTE_T_Z
  private:: KNOWING_T_COMPUTE_P_Z
  private:: KNOWING_Z_COMPUTE_T_P
  private:: GENERIC_PROFILE_INTERPOLATION
  private:: OPTIMAL_ESTIMATION
  private:: COMPUTE_FORWARD_MODEL_AND_KERNEL
  private:: COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO
  private:: DETERMINE_SFC_TYPE_FORWARD_MODEL
  private:: SET_CLEAR_SKY_COVARIANCE_TERMS
  private:: COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE
  private:: COMPUTE_CIRRUS_APRIORI

  private:: DETERMINE_ACHA_MODE_BASED_ON_CHANNELS
  private:: INTERPOLATE_PROFILE_ACHA
  private:: DETERMINE_OPAQUE_CLOUD_HEIGHT
  private:: COMPUTE_REFERENCE_LEVEL_EMISSIVITY
  private:: COMPUTE_STANDARD_DEVIATION
  private:: NULL_PIX_POINTERS 
  private:: COMPUTE_TEMPERATURE_CIRRUS
  private:: COMPUTE_BOX_WIDTH
  private:: MEAN_SMOOTH
  private:: EMPIRICAL_LAPSE_RATE
  private:: OCEANIC_LAPSE_RATE_OLD

  !--- include the non-system specific variables
  include 'acha_parameters.inc'

  !--- interpolated profiles
  real, private, dimension(Num_Levels_RTM_Prof) :: Temp_Prof_RTM
  real, private, dimension(Num_Levels_RTM_Prof) :: Press_Prof_RTM
  real, private, dimension(Num_Levels_RTM_Prof) :: Hght_Prof_RTM
  integer, private, dimension(Num_Levels_RTM_Prof) :: Inver_Prof_RTM
  integer, private:: Inver_Top_Level_RTM
  integer, private:: Inver_Base_Level_RTM
  integer, private:: Sfc_Level_RTM
  integer, private:: Tropo_Level_RTM
  real, private:: Inver_Top_Height
  real, private:: Inver_Base_Height
  real, private:: Inver_Strength

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
  real, private, PARAMETER:: MISSING_VALUE_REAL4 = -999.0
  integer(kind=int1), private, PARAMETER:: MISSING_VALUE_INTEGER1 = -128
  integer(kind=int4), private, PARAMETER:: MISSING_VALUE_INTEGER4 = -999
  type(acha_symbol_struct), private :: symbol


  !-------------------------------------------------------------------------------------
  ! empirical lapse rate table data and metadata
  !-------------------------------------------------------------------------------------
  integer, private, parameter:: nts = 7
  integer, private, parameter:: ntcs = 9
  real, private, parameter:: ts_min = 270.0
  real, private, parameter:: dts = 5.0
  real, private, parameter:: tcs_min = -20.0
  real, private, parameter:: dtcs = 2.0

  real, private, dimension(nts,ntcs), parameter:: ocean_lapse_rate_table = reshape ((/ &
                          -7.3, -7.2, -7.3, -7.4, -7.4, -6.8, -6.2, &
                          -7.4, -7.3, -7.3, -7.4, -7.4, -7.0, -6.3, &
                          -7.5, -7.3, -7.3, -7.5, -7.6, -7.1, -6.5, &
                          -7.2, -7.1, -7.3, -7.5, -7.6, -7.2, -6.6, &
                          -6.9, -6.8, -7.1, -7.4, -7.5, -7.3, -7.0, &
                          -6.6, -6.6, -6.8, -7.0, -7.3, -7.4, -7.4, &
                          -6.7, -6.4, -6.4, -6.6, -7.0, -7.3, -7.6, &
                          -6.2, -5.8, -5.6, -5.8, -6.3, -6.8, -7.3, &
                          -5.8, -5.3, -5.0, -5.2, -5.9, -6.3, -6.8/), (/nts,ntcs/))

  real, private, dimension(nts,ntcs), parameter:: land_lapse_rate_table = reshape ((/ &
                           -5.2, -5.8, -6.2, -6.2, -6.4, -7.0, -7.7, &
                           -5.3, -5.8, -6.2, -6.3, -6.4, -7.1, -7.7, &
                           -5.2, -5.7, -6.0, -6.1, -6.4, -7.1, -7.7, &
                           -5.0, -5.4, -5.8, -5.9, -6.2, -6.9, -7.7, &
                           -5.0, -5.2, -5.5, -5.5, -5.8, -6.8, -7.8, &
                           -4.9, -5.0, -5.2, -4.9, -5.2, -6.2, -7.6, &
                           -4.7, -4.7, -4.8, -4.5, -4.8, -6.0, -7.5, &
                           -3.9, -4.0, -4.2, -3.9, -3.9, -5.3, -7.3, &
                           -3.3, -3.4, -3.7, -3.6, -3.5, -5.0, -7.3/), (/nts,ntcs/))

  ! surface emissivity
  real(kind=real4),private:: Emiss_Sfc_85um
  real(kind=real4),private:: Emiss_Sfc_11um
  real(kind=real4),private:: Emiss_Sfc_12um
  real(kind=real4),private:: Emiss_Sfc_133um
  real(kind=real4),private:: Emiss_Sfc_67um

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
! 6 - Lower Cloud InterpoLation Used (0 = no / 1 = ! yes)
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
! October 2006 - Added cloud lapse rate to make Tc more reLated to true 
!                cloud-top temperature
!
!
!------------------------------------------------------------------------------
  subroutine  AWG_CLOUD_HEIGHT_ALGORITHM(Input, Symbol_In, Output, Diag)

  !===============================================================================
  !  Argument Declaration
  !==============================================================================

  type(acha_input_struct), intent(inout) :: Input
  type(acha_symbol_struct), intent(in) :: Symbol_In
  type(acha_output_struct), intent(inout) :: Output
  type(acha_diag_struct), intent(inout), optional :: Diag

  !===============================================================================
  !  Pixel level RTM structure
  !===============================================================================
 
  type(acha_rtm_nwp_struct) :: ACHA_RTM_NWP

  !===============================================================================
  !  Local Variable Declaration
  !===============================================================================

  integer :: ACHA_Mode_Flag
  integer:: Elem_Idx
  integer:: Line_Idx
  integer:: Param_Idx
  integer:: i
  integer:: Singular_Flag
  integer:: Lev_Idx
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
  integer:: ierror
  integer:: Pass_Idx
  integer:: Pass_Idx_Min
  integer:: Pass_Idx_Max
  real:: Lapse_Rate
  real:: Abs_Lapse_Rate_dP_dZ
  real:: Delta_Cld_Temp_Sfc_Temp

  real:: Convergence_Criteria

  real:: Bs
  real:: Rad_Atm
  real:: Trans_Atm

  !--- ch27 variables
  real:: Rad_Ac_67um
  real:: Trans_Ac_67um
  real:: Trans_Bc_67um
  real:: Rad_Clear_67um

  !--- ch29 variables
  real:: Rad_Ac_85um
  real:: Trans_Ac_85um
  real:: Trans_Bc_85um
  real:: Rad_Clear_85um

  !--- ch31 variables
  real:: Rad_Ac_11um
  real:: Trans_Ac_11um
  real:: Trans_Bc_11um
  real:: Rad_Clear_11um

  !--- ch32 variables
  real:: Rad_Ac_12um
  real:: Trans_Ac_12um
  real:: Trans_Bc_12um
  real:: Rad_Clear_12um

  !--- ch33 variables
  real:: Rad_Ac_133um
  real:: Trans_Ac_133um
  real:: Trans_Bc_133um
  real:: Rad_Clear_133um

  real:: a_Beta_11um_133um_fit
  real:: b_Beta_11um_133um_fit
  real:: a_Beta_11um_85um_fit
  real:: b_Beta_11um_85um_fit
  real:: a_Beta_11um_67um_fit
  real:: b_Beta_11um_67um_fit

  real:: Tsfc_Est
  real:: Ts_temp
  real:: Tc_temp
  real:: Pc_temp
  real:: Ps_temp
  real:: Zc_Temp
  real:: Zs_Temp
  real:: Tc_Ap
  real:: Ec_Ap
  real:: Beta_Ap
  real:: Ts_Ap
  real:: Tc_Ap_Uncer
  real:: Ts_Ap_Uncer
  real:: Ec_Ap_Uncer
  real:: Tc_Ap_Imager
  real:: Tc_Ap_Sounder
  real:: Beta_Ap_Uncer
  real:: Sa_Tc_Imager
  real:: Sa_Tc_Sounder
  real:: T11um_Clr_Uncer
  real:: T11um_12um_Clr_Uncer
  real:: T11um_133um_Clr_Uncer
  real:: T11um_85um_Clr_Uncer
  real:: T11um_67um_Clr_Uncer
  real:: R4_Dummy
  real:: Emiss_11um_Tropo
  integer:: Num_Obs
  integer(kind=int1):: Cloud_Type
  integer:: Cloud_Phase
  integer:: Undetected_Cloud
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

  integer(kind=int4), dimension(:,:), allocatable:: Fail_Flag
  integer(kind=int4), dimension(:,:), allocatable:: Converged_Flag
  real (kind=real4), allocatable, dimension(:,:):: Temperature_Cirrus
  real (kind=real4), allocatable, dimension(:,:):: Temperature_Lower_Cloud_Apriori
  integer (kind=int4):: Box_Half_Width_Cirrus
  integer (kind=int4):: Box_Half_Width_Lower
  real (kind=real4), allocatable, dimension(:,:):: Pc_Opaque
  real (kind=real4), allocatable, dimension(:,:):: Zc_Opaque
  real (kind=real4), allocatable, dimension(:,:):: Tc_Opaque

  !--- local POINTERs to global arrays or data structures
  integer(kind=int4), allocatable, dimension(:,:):: Elem_Idx_LRC
  integer(kind=int4), allocatable, dimension(:,:):: Line_Idx_LRC
  integer(kind=int1), allocatable, dimension(:,:):: Skip_LRC_Mask
  real (kind=real4):: Bt_11um_Lrc
  real (kind=real4):: Bt_11um_Std
  real (kind=real4):: Btd_11um_67um_Std
  real (kind=real4):: Btd_11um_85um_Std
  real (kind=real4):: Btd_11um_12um_Std
  real (kind=real4):: Btd_11um_133um_Std

  !--- scalar local variables
  integer (kind=int4):: i1,i2,j1,j2
  integer (kind=int4):: NWP_Profile_Inversion_Flag
  logical:: Bad_Input_Flag

  integer:: count_diag
  integer:: lun_diag
  character(len=100):: file_name_diag

!-----------------------------------------------------------------------
! BEGIN EXECUTABLE CODE
!-----------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! copy symbol to a module-wide structure
   !--------------------------------------------------------------------
   symbol = Symbol_In 

  !----------------------------------------------------------------------------
  ! abort if no 11 um channel
  !----------------------------------------------------------------------------
  if (Input%Chan_On_11um == symbol%NO) THEN
     Output%Packed_Qf =  0
     return
  endif

  !--- initialize diagnostic output
  if (present(Diag)) Diag%Array_1 = Missing_Value_Real4
  if (present(Diag)) Diag%Array_2 = Missing_Value_Real4
  if (present(Diag)) Diag%Array_3 = Missing_Value_Real4
  
  !---------------------------------------------------------------------------
  !-- setup microphysical models
  !---------------------------------------------------------------------------
  call SETUP_ICE_MICROPHYSICAL_MODEL(Input%WMO_Id)

  !---------------------------------------------------------------------------
  !-- Acha Mode set to  -1, determine based on channels
  !---------------------------------------------------------------------------
  ACHA_Mode_Flag = Input%ACHA_Mode_Flag_In
  if (ACHA_Mode_Flag < 0) then
    call DETERMINE_ACHA_MODE_BASED_ON_CHANNELS( &
                                      Acha_Mode_Flag, &
                                      Input%Chan_On_67um, &
                                      Input%Chan_On_85um,  &
                                      Input%Chan_On_11um,  &
                                      Input%Chan_On_12um,  &
                                      Input%Chan_On_133um)
  endif

  !--- determine number of channels
  select case(Acha_Mode_Flag)
     case(1)  !11 avhrr/1
       Num_Obs = 1
     case(2)  !11,6.7
       Num_Obs = 2
     case(3)  !11,12 avhrr/2/3
       Num_Obs = 2
     case(4)  !11,13.3 goes-nop
       Num_Obs = 2
     case(5)  !11,12,8.5 viirs 
       Num_Obs = 3
     case(6)  !11,12,6.7
       Num_Obs = 3
     case(7)  !11,13.3,6.7
       Num_Obs = 3
     case(8)  !11,12,13.3 goes-r
       Num_Obs = 3
     case(9)  !11,12,13.3 goes-r with pseudo 13.3
       Num_Obs = 3
  end select

  !--- allocate needed 2d arrays for processing this segment
  allocate(Elem_Idx_LRC(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Line_Idx_LRC(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Skip_LRC_Mask(Input%Number_of_Elements,Input%Number_of_Lines))

  !--- allocate array for cirrus temperature
  allocate(Fail_Flag(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Converged_Flag(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Temperature_Cirrus(Input%Number_of_Elements,Input%Number_of_Lines))

  allocate(Pc_Opaque(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Tc_Opaque(Input%Number_of_Elements,Input%Number_of_Lines))
  allocate(Zc_Opaque(Input%Number_of_Elements,Input%Number_of_Lines))

  !--- allocate array to hold lowe cloud temp
  allocate(Temperature_Lower_Cloud_Apriori(Input%Number_of_Elements,Input%Number_of_Lines))

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

  !--- initialize diagnostic output
  count_diag = 0
  lun_diag = 0

  !--- set convergence criterion
  Convergence_Criteria = (Num_Param - 1.0) / 5.0

  !--- determine cirrus spatial interpolation box width
  call COMPUTE_BOX_WIDTH(Input%Sensor_Resolution_KM,CIRRUS_BOX_WIDTH_KM, Box_Half_Width_Cirrus)

  !--- determine lower cloud spatial interpolation box width
  call COMPUTE_BOX_WIDTH(Input%Sensor_Resolution_KM,LOWER_BOX_WIDTH_KM, Box_Half_Width_Lower)

  !----------- make identity matrix
  E = 0.0
  do i = 1,Num_Param
    E(i,i) = 1.0
  enddo

  !--- initialize output
  Output%Tc =  MISSING_VALUE_REAL4
  Output%Ec =  MISSING_VALUE_REAL4
  Output%Beta =  MISSING_VALUE_REAL4
  Output%Pc =  MISSING_VALUE_REAL4
  Output%Zc =  MISSING_VALUE_REAL4
  Output%OE_Qf = 0
  Output%Qf = 0
  Meta_Data_Flags = 0
  Output%Inversion_Flag = 0
  if (Input%Chan_On_67um == symbol%YES) Output%Ec_67um = MISSING_VALUE_REAL4
  if (Input%Chan_On_85um == symbol%YES) Output%Ec_85um = MISSING_VALUE_REAL4
  if (Input%Chan_On_11um == symbol%YES) Output%Ec_11um = MISSING_VALUE_REAL4
  if (Input%Chan_On_12um == symbol%YES) Output%Ec_12um = MISSING_VALUE_REAL4
  if (Input%Chan_On_133um == symbol%YES) Output%Ec_133um = MISSING_VALUE_REAL4
  
  !--------------------------------------------------------------------------
  ! spatial processing pixels
  ! compute local radiative centers using 11 um brightness temperature
  !---------------------------------------------------------------------------

  !--- construct a mask to select pixel for LRC computation
  Elem_Idx_LRC = MISSING_VALUE_INTEGER4
  Line_Idx_LRC = MISSING_VALUE_INTEGER4
  Skip_LRC_Mask = Input%Invalid_Data_Mask
  Temperature_Cirrus = MISSING_VALUE_REAL4
  Pc_Opaque = MISSING_VALUE_REAL4
  Tc_Opaque = MISSING_VALUE_REAL4
  Zc_Opaque = MISSING_VALUE_REAL4
  Temperature_Lower_Cloud_Apriori = MISSING_VALUE_REAL4

  !--- call LRC routine
  if (Use_Lrc_Flag == symbol%YES) then

    if (associated(Input%Elem_Idx_LRC_Input) .and. &
        associated(Input%Line_Idx_LRC_Input)) then

      Elem_Idx_LRC = Input%Elem_Idx_LRC_Input
      Line_Idx_LRC = Input%Line_Idx_LRC_Input

    else

      where(Input%Cloud_Mask == symbol%CLEAR .or. Input%Cloud_Mask == symbol%PROB_CLEAR)
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
                                         MISSING_VALUE_INTEGER4, &
                                         Skip_LRC_Mask, &
                                         Min_Bt_11um_Lrc,  &
                                         Max_Bt_11um_Lrc, &
                                         Elem_Idx_LRC,  &
                                         Line_Idx_LRC)

    endif
  endif

  !--------------------------------------------------------------------------
  ! Multi-Layer Logic Implemented via cloud type
  !-------------------------------------------------------------------------

  Output%Cloud_Type = Input%Cloud_Type
 
  if (MULTI_LAYER_LOGIC_FLAG == 1) then 
   where(Input%Cloud_Type == symbol%OVERLAP_TYPE)
     Output%Cloud_Type = symbol%CIRRUS_TYPE
   endwhere
  endif

  if (MULTI_LAYER_LOGIC_FLAG == 2) then 
   where(Input%Cloud_Type == symbol%CIRRUS_TYPE)
     Output%Cloud_Type = symbol%OVERLAP_TYPE
   endwhere
  endif

  !--------------------------------------------------------------------------
  ! determine processing order of pixels
  !--------------------------------------------------------------------------
  call COMPUTE_PROCESSING_ORDER(&
                                Input%Invalid_Data_Mask, Output%Cloud_Type,&
                                ELem_Idx_LRC,Line_Idx_LRC, &
                                Pass_Idx_Min,Pass_Idx_Max,USE_CIRRUS_FLAG, &
                                Output%Processing_Order) 

  !--------------------------------------------------------------------------
  ! Loop through pixels using the processing order
  !--------------------------------------------------------------------------

  pass_loop: do Pass_Idx = Pass_Idx_min, Pass_Idx_Max
  
   !--------------------------------------------------------------------------
   ! on the third pass, spatially interpolate water cloud temperature
   ! note, this first guess is stored in the Output Variable but it is
   ! over-written during the retrieval
   !--------------------------------------------------------------------------
   if ((Pass_Idx == 0) .or. (Pass_Idx == 3)) then

     call  COMPUTE_LOWER_CLOUD_TEMPERATURE(Output%Cloud_Type, &
                                        USE_LOWER_INTERP_FLAG, &
                                        Input%Surface_Temperature, &
                                        Output%Tc,&
                                        COUNT_MIN_LOWER,      &
                                        Box_Half_Width_Lower, &
                                        MISSING_VALUE_REAL4, &
                                        Temperature_Lower_Cloud_Apriori)
   endif

   !--------------------------------------------------------------------------
   ! loop over pixels in scanlines
   !--------------------------------------------------------------------------
   Line_loop: do Line_Idx = Line_Idx_min,Input%Number_of_Lines + Line_Idx_min - 1

    Element_Loop:   do Elem_Idx = 1, Input%Number_of_Elements

    !---- null profile pointers each time 
    call NULL_PIX_POINTERS(Input, ACHA_RTM_NWP)

    !--- check if pixel should be processd in this path
    if (USE_CIRRUS_FLAG == symbol%NO .or. Pass_Idx /= Pass_Idx_Max) then
      if (Pass_Idx /= Output%Processing_Order(Elem_Idx,Line_Idx)) then
          cycle
      endif
    endif

    !---------------------------------------------------------------
    ! Check to see if this pixel shoud be skipped
    !---------------------------------------------------------------
    Bad_Input_Flag = .false.

    if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) Bad_Input_Flag = .true.

    if (Input%Sensor_Zenith_Angle(Elem_Idx,Line_Idx) > Sensor_Zenith_Threshold) Bad_Input_Flag = .true.

    if (Input%Bt_11um(ELem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.

    if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 7) then
        if (Input%Bt_67um(ELem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (Acha_Mode_Flag == 5) then
        if (Input%Bt_85um(ELem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 8 .or. Acha_Mode_Flag == 9) then
        if (Input%Bt_12um(ELem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif
    if (Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 8 .or. Acha_Mode_Flag == 9) then
        if (Input%Bt_133um(ELem_Idx,Line_Idx) == MISSING_VALUE_REAL4) Bad_Input_Flag = .true.
    endif

    if (Bad_Input_Flag .eqv. .true.) then
          Output%Packed_Qf(Elem_Idx,Line_Idx) =  0
          Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0
          cycle
    endif

    !--- for convenience, save nwp indices to local variables
    Inwp = Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx)
    Jnwp = Input%Line_Idx_Nwp(Elem_Idx,Line_Idx)
    Inwp_x = Input%Elem_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
    Jnwp_x = Input%Line_Idx_Opposite_Corner_NWP(Elem_Idx,Line_Idx)
    Inwp_Weight = Input%Longitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
    Jnwp_Weight = Input%Latitude_Interp_Weight_NWP(Elem_Idx,Line_Idx)
    Ivza =  Input%Viewing_Zenith_Angle_Idx_RTM(Elem_Idx,Line_Idx)
    ilrc = Elem_Idx_LRC(Elem_Idx,Line_Idx)
    jlrc = Line_Idx_LRC(Elem_Idx,Line_Idx)
    Cloud_Type = Output%Cloud_Type(Elem_Idx,Line_Idx)
    
    !--- Qc indices
    if (Input%Elem_Idx_Nwp(Elem_Idx,Line_Idx) <= 0 .or. &
        Input%Line_Idx_Nwp(Elem_Idx,Line_Idx) <= 0 .or. &
        Input%Viewing_Zenith_Angle_Idx_RTM(Elem_Idx,Line_Idx) <= 0) then 
          Output%Packed_Qf(Elem_Idx,Line_Idx) =  0
          Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0
         cycle 
    endif

    !---  filter pixels for last pass for cirrus correction
    if (Pass_Idx == Pass_Idx_Max .and. USE_CIRRUS_FLAG == symbol%YES) then

        if (Output%Cloud_Type(Elem_Idx,Line_Idx) /= symbol%CIRRUS_TYPE .and. &
            Output%Cloud_Type(Elem_Idx,Line_Idx) /= symbol%OVERLAP_TYPE) then
             cycle
        endif

        !--- don't redo cirrus with valid lrc values
        if (ilrc > 0 .and. jlrc > 0) then 
            if (Output%Ec(ilrc,jlrc) > 0.7) cycle
        endif

    endif

    !-----------------------------------------------------------------------
    ! include code to setup local profiles correctly 
    !-----------------------------------------------------------------------
    
    !Call framework services module
    call ACHA_FETCH_PIXEL_NWP_RTM(Input, symbol, &
                                  Elem_Idx,Line_Idx, ACHA_RTM_NWP)
    
    
    Sfc_Level_RTM = ACHA_RTM_NWP%Sfc_Level
    Tropo_Level_RTM = ACHA_RTM_NWP%Tropo_Level
    
    Press_Prof_RTM =  ACHA_RTM_NWP%P_Prof

    !do smoothing routines here 
    if (ACHA_RTM_NWP%Smooth_Nwp_Fields_Flag_Temp == symbol%YES) then

       !--- height profile       
       Hght_Prof_RTM = INTERPOLATE_PROFILE_ACHA( ACHA_RTM_NWP%Z_Prof, &
                                                 ACHA_RTM_NWP%Z_Prof_1, &
                                                 ACHA_RTM_NWP%Z_Prof_2, &
                                                 ACHA_RTM_NWP%Z_Prof_3, &
                                                 Inwp_Weight,Jnwp_Weight)

      !--- temperature profile
      Temp_Prof_RTM = INTERPOLATE_PROFILE_ACHA( ACHA_RTM_NWP%T_Prof, &
                                                ACHA_RTM_NWP%T_Prof_1, &
                                                ACHA_RTM_NWP%T_Prof_2, &
                                                ACHA_RTM_NWP%T_Prof_3, &
                                                Inwp_Weight,Jnwp_Weight)
    else

       Hght_Prof_RTM = ACHA_RTM_NWP%Z_Prof
       Temp_Prof_RTM = ACHA_RTM_NWP%T_Prof
    
    endif

   !-----------------------------------------------------------------------
   !  find opaque cloud height
   !-----------------------------------------------------------------------
   call DETERMINE_OPAQUE_CLOUD_HEIGHT( &
                                Input%Rad_11um(Elem_Idx,Line_Idx), &
                                ACHA_RTM_NWP%Black_Body_Rad_Prof_11um, &
                                Press_Prof_RTM, &
                                Hght_Prof_RTM, &
                                Temp_Prof_RTM, &
                                Tropo_Level_RTM, &
                                Sfc_Level_RTM, &
                                Pc_Opaque(Elem_Idx,Line_Idx), &
                                Tc_Opaque(Elem_Idx,Line_Idx), &
                                Zc_Opaque(Elem_Idx,Line_Idx))

  !-------------------------------------------------------------------
  ! Apply Opaque Retrieval for Acha_Mode_Flag = 1, then cycle
  !-------------------------------------------------------------------
  if (Acha_Mode_Flag == 1) then
        if (((Input%Cloud_Mask(Elem_Idx,Line_Idx) == symbol%CLEAR) .or.  &
            (Input%Cloud_Mask(Elem_Idx,Line_Idx) == symbol%PROB_CLEAR)) .and. &
            (Input%Process_Undetected_Cloud_Flag == symbol%NO)) then
          Output%Tc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
          Output%Pc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
          Output%Zc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
          Output%Ec(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
          Output%Beta(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        else
          Output%Tc(Elem_Idx,Line_Idx) = Tc_Opaque(Elem_Idx,Line_Idx)
          Output%Pc(Elem_Idx,Line_Idx) = Pc_Opaque(Elem_Idx,Line_Idx)
          Output%Zc(Elem_Idx,Line_Idx) = Zc_Opaque(Elem_Idx,Line_Idx)
          Output%Ec(Elem_Idx,Line_Idx) = 1.0
          Output%Beta(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        endif
        Output%Tc_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Pc_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Zc_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Ec_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Beta_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
        Output%Packed_Qf(Elem_Idx,Line_Idx) =  0
        Output%Packed_Meta_Data(Elem_Idx,Line_Idx) =  0
        cycle
   endif


   !----------------------------------------------------------------------
   ! determine cloud phase from cloud type for convienience
   !----------------------------------------------------------------------
   Cloud_Phase = symbol%UNKNOWN_PHASE

   if ( (Cloud_Type  == symbol%FOG_TYPE) .or. &
       (Cloud_Type  == symbol%WATER_TYPE) .or. &
       (Cloud_Type  == symbol%SUPERCOOLED_TYPE)) then
       Cloud_Phase = symbol%WATER_PHASE
   endif

   if ( (Cloud_Type  == symbol%CIRRUS_TYPE) .or. &
        (Cloud_Type  == symbol%OVERLAP_TYPE) .or. &
        (Cloud_Type  == symbol%TICE_TYPE) .or. &
        (Cloud_Type  == symbol%OVERSHOOTING_TYPE)) then
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

     Output%Tc(Elem_Idx,Line_Idx) =  MISSING_VALUE_REAL4
     Output%Ec(Elem_Idx,Line_Idx) =  MISSING_VALUE_REAL4
     Output%Beta(Elem_Idx,Line_Idx) =  MISSING_VALUE_REAL4
     Output%Pc(Elem_Idx,Line_Idx) =  MISSING_VALUE_REAL4
     Output%Zc(Elem_Idx,Line_Idx) =  MISSING_VALUE_REAL4
     Output%Qf(Elem_Idx,Line_Idx) = 0

   else  !if passed data check then proceed with retrieval

     !---------------------------------------------------------------------
     ! select to do retrievals for all pixels or just cloudy ones
     !---------------------------------------------------------------------
     Undetected_Cloud = symbol%NO

     if ((Input%Cloud_Mask(Elem_Idx,Line_Idx) == symbol%CLEAR) .or.  &
         (Input%Cloud_Mask(Elem_Idx,Line_Idx) == symbol%PROB_CLEAR)) then

       if (Input%Process_Undetected_Cloud_Flag == symbol%NO) then
               cycle
       else
              Undetected_Cloud = symbol%YES
       endif

     endif


    !----------------------------------------------------------------------
    !  turn on diagnostic
    !----------------------------------------------------------------------
    lun_diag = 0
    if (idiag_output == symbol%YES .and. count_diag < count_diag_max) then
      print *, "ACHA Diagnostic File Created for this Pixel"
      count_diag = count_diag + 1
      write(file_name_diag,fmt="(I0.4)") count_diag
      file_name_diag = trim(file_name_diag)
      file_name_diag = 'acha_diag_output_'//trim(file_name_diag)//'.txt'
      lun_diag = get_lun_acha()
      open(unit=lun_diag,file=trim(file_name_diag),form="formatted",status='unknown',action='write') 
    endif

    if (lun_diag > 0) then
         write(unit=lun_diag,fmt=*) "Element, Line Indices (Relative to Segment) = ", Elem_Idx,Line_Idx
         write(unit=lun_diag,fmt=*) "Surface Elevation = ", Input%Surface_Elevation(Elem_Idx,Line_Idx)
         write(unit=lun_diag,fmt=*) "Latitude = ", Input%Latitude(Elem_Idx,Line_Idx)
         write(unit=lun_diag,fmt=*) "Longitude = ", Input%Longitude(Elem_Idx,Line_Idx)
         write(unit=lun_diag,fmt=*) "Zenith Angle = ", Input%Sensor_Zenith_Angle(Elem_Idx,Line_Idx)
         if (Input%Surface_Type(Elem_Idx,Line_Idx) == 0) write(unit=lun_diag,fmt=*) "WATER SURFACE"
         if (Input%Surface_Type(Elem_Idx,Line_Idx) > 0) write(unit=lun_diag,fmt=*) "LAND SURFACE"
         if (Input%SNOW_CLASS(Elem_Idx,Line_Idx) == 1) write(unit=lun_diag,fmt=*) "UNFROZEN SURFACE"
         if (Input%SNOW_CLASS(Elem_Idx,Line_Idx) > 1) write(unit=lun_diag,fmt=*) "FROZEN SURFACE"
         write(unit=lun_diag,fmt=*) "Cloud Type = ", Cloud_Type
         write(unit=lun_diag,fmt=*) "ACHA_MODE = ", ACHA_Mode_Flag
         write(unit=lun_diag,fmt=*) "Pass_Idx = ", Pass_Idx
         write(unit=lun_diag,fmt=*) "Clear 67 um Rad = ", Input%Rad_Clear_67um(Elem_Idx,Line_Idx)
         write(unit=lun_diag,fmt=*) "Clear 11 um Rad = ", Input%Rad_Clear_11um(Elem_Idx,Line_Idx)
         write(unit=lun_diag,fmt=*) "Clear 13.3 um Rad = ", Input%Rad_Clear_133um(Elem_Idx,Line_Idx)
         write(unit=lun_diag,fmt=*) "Clear 13.6 um Rad = ", Input%Rad_Clear_136um(Elem_Idx,Line_Idx)
         write(unit=lun_diag,fmt=*) "Clear 13.9 um Rad = ", Input%Rad_Clear_139um(Elem_Idx,Line_Idx)
         write(unit=lun_diag,fmt=*) "Clear 14.2 um Rad = ", Input%Rad_Clear_142um(Elem_Idx,Line_Idx)
         write(unit=lun_diag,fmt=*) "Obs 67 um Rad = ", PLANCK_RAD_FAST(Input%Chan_Idx_67um,Input%Bt_67um(Elem_Idx,Line_Idx))
         write(unit=lun_diag,fmt=*) "Obs 11 um Rad = ", PLANCK_RAD_FAST(Input%Chan_Idx_11um,Input%Bt_11um(Elem_Idx,Line_Idx))
         write(unit=lun_diag,fmt=*) "Obs 13.3 um Rad = ", PLANCK_RAD_FAST(Input%Chan_Idx_133um,Input%Bt_133um(Elem_Idx,Line_Idx))
         write(unit=lun_diag,fmt=*) "Obs 13.6 um Rad = ", PLANCK_RAD_FAST(Input%Chan_Idx_136um,Input%Bt_136um(Elem_Idx,Line_Idx))
         write(unit=lun_diag,fmt=*) "Obs 13.9 um Rad = ", PLANCK_RAD_FAST(Input%Chan_Idx_139um,Input%Bt_139um(Elem_Idx,Line_Idx))
         write(unit=lun_diag,fmt=*) "Obs 14.2 um Rad = ", PLANCK_RAD_FAST(Input%Chan_Idx_142um,Input%Bt_142um(Elem_Idx,Line_Idx))
         write(unit=lun_diag,fmt=*) "Profiles: idx, P, Z, T"
         do Lev_Idx = 1,Num_Levels_Rtm_Prof 
          write(unit=lun_diag,fmt="(I3,F8.2,F8.1,F8.3,6F8.2)") Lev_Idx,  &
                   Press_Prof_RTM(Lev_Idx), Hght_Prof_RTM(Lev_Idx), Temp_Prof_RTM(Lev_Idx), &
                   Acha_RTM_NWP%Black_Body_Rad_Prof_67um(Lev_Idx), Acha_RTM_NWP%Black_Body_Rad_Prof_11um(Lev_Idx), &
                   Acha_RTM_NWP%Black_Body_Rad_Prof_133um(Lev_Idx), &
                   Acha_RTM_NWP%Black_Body_Rad_Prof_136um(Lev_Idx), &
                   Acha_RTM_NWP%Black_Body_Rad_Prof_139um(Lev_Idx), &
                   Acha_RTM_NWP%Black_Body_Rad_Prof_142um(Lev_Idx)
         enddo
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
   Meta_Data_Flags(6) = symbol%NO     !lower cloud interpoLation
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

   if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 7) then
    Btd_11um_67um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_11um(i1:i2,j1:j2) -  Input%Bt_67um(i1:i2,j1:j2),&
                                                   Input%Invalid_Data_Mask(i1:i2,j1:j2))
   endif
   if (Acha_Mode_Flag == 5) then
    Btd_11um_85um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_11um(i1:i2,j1:j2) -  Input%Bt_85um(i1:i2,j1:j2), &
                                                    Input%Invalid_Data_Mask(i1:i2,j1:j2))
   endif

   if (Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 5 .or.  &
       Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 8 .or.  &
       Acha_Mode_Flag == 9) then
    Btd_11um_12um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_11um(i1:i2,j1:j2) -  Input%Bt_12um(i1:i2,j1:j2), &
                                                   Input%Invalid_Data_Mask(i1:i2,j1:j2))
   endif

   if (Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 7 .or. Acha_Mode_Flag == 8 .or. Acha_Mode_Flag == 9) then
    Btd_11um_133um_Std = COMPUTE_STANDARD_DEVIATION( Input%Bt_11um(i1:i2,j1:j2) -  Input%Bt_133um(i1:i2,j1:j2), &
                                                    Input%Invalid_Data_Mask(i1:i2,j1:j2))
   endif
   
   !-----------------------------------------------------------------------
   ! assign values to y and x_Ap
   !----------------------------------------------------------------------

   !--- y - the observation vOutput%Ector
   select case(Acha_Mode_Flag)
     case(1)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
     case(2)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_67um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_67um_Std**2 
     case(3)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
     case(4)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_133um_Std**2 
     case(5)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_85um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
       y_variance(3) = Btd_11um_85um_Std**2 
     case(6)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_67um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(3) = Btd_11um_12um_Std**2 
       y_variance(3) = Btd_11um_67um_Std**2 
     case(7)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_67um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_133um_Std**2 
       y_variance(3) = Btd_11um_67um_Std**2 
     case(8,9)
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y(3) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_133um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
       y_variance(3) = Btd_11um_133um_Std**2 
     case DEFAULT
       y(1) =  Input%Bt_11um(Elem_Idx,Line_Idx)
       y(2) =  Input%Bt_11um(Elem_Idx,Line_Idx) -  Input%Bt_12um(Elem_Idx,Line_Idx)
       y_variance(1) =  Bt_11um_Std**2
       y_variance(2) = Btd_11um_12um_Std**2 
   end select

   !-------------------------------------------------------------------
   ! Determine surface type for use in forward model
   ! 0 = Water
   ! 1 = Land
   ! 2 = Snow
   ! 3 = Desert
   ! 4 = Arctic
   ! 5 = Antarctic
   !-------------------------------------------------------------------
   call DETERMINE_SFC_TYPE_FORWARD_MODEL(Input%Surface_Type(Elem_Idx,Line_Idx), &
                                         Input%Snow_Class (Elem_Idx,Line_Idx), &
                                         Input%Latitude(Elem_Idx,Line_Idx), &
                                         Input%Surface_Emissivity_39um(Elem_Idx,Line_Idx), &
                                         Sfc_Type_Forward_Model)

   !-------------------------------------------------------------------
   ! Based on fm surface type, set the clear-sky covariance terms
   !-------------------------------------------------------------------
   call SET_CLEAR_SKY_COVARIANCE_TERMS(Sfc_Type_Forward_Model)

   !-------------------------------------------------------------------
   ! These values are used in the baseline code.  In the Latest
   ! code, the values from SET_CLEAR_SKY_COVARIANCE_TERMS are used
   !
   ! these are based on patmos-x clear data and are the 
   ! approx average of des+asc from August 2006 NOAA-18
   !-------------------------------------------------------------------
   !REMOVE
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

  !--------------------------------------------------------------------
  ! pick a priori conditions
  !--------------------------------------------------------------------

  !--- logic for unmasked or untyped pixels (UndetOutput%Ected cloud)
  if (Undetected_Cloud == symbol%YES) then
         if (Tc_Opaque(Elem_Idx,Line_Idx) < 260.0 .and.  &
             Tc_Opaque(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then
             Cloud_Type = symbol%CIRRUS_TYPE
         else
             Cloud_Type = symbol%FOG_TYPE
         endif
  endif

  !---- Compute 11um emissivity referenced to tropopause
  Emiss_11um_Tropo = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( &
                             Tropo_Level_RTM, &
                             Input%Rad_11um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_11um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_11um)

  !---- select Output%Tc and Output%Ec apriori based on cloud type

  if ((ilrc /= MISSING_VALUE_INTEGER4) .and. &
      (jlrc /= MISSING_VALUE_INTEGER4)) then
           Bt_11um_Lrc =  Input%Bt_11um(ilrc,jlrc)
  else
           Bt_11um_Lrc = MISSING_VALUE_REAL4
  endif

  call COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO( &
                       Cloud_Phase, &
                       Emiss_11um_Tropo, &
                       Input%Latitude(Elem_Idx,Line_Idx), &
                       Input%Tropopause_Temperature(Elem_Idx,Line_Idx), &
                       Input%Bt_11um(Elem_Idx,Line_Idx), &
                       Bt_11um_Lrc, &
                       Tc_Opaque(Elem_Idx,Line_Idx), &
                       Input%Cosine_Zenith_Angle(Elem_Idx,Line_Idx), &
                       Input%Snow_Class(Elem_Idx,Line_Idx), &
                       Input%Surface_Air_Temperature(Elem_Idx,Line_Idx), &
                       Tc_Ap,Tc_Ap_Uncer, &
                       Ec_Ap,Ec_Ap_Uncer, &
                       Beta_Ap,Beta_Ap_Uncer)

   if (lun_diag > 0) then 
     write(unit=lun_diag,fmt=*) "==========================================================="
     write(unit=lun_diag,fmt=*) "START OF OE RETRIEVAL"
     write(unit=lun_diag,fmt=*) "==========================================================="
     write(unit=lun_diag,fmt=*) "y = ", y
     write(unit=lun_diag,fmt=*) "initial x_ap = ", Tc_Ap, Ec_Ap, Beta_Ap, Ts_Ap
   endif 

  !------------------------------------------------------------------------
  ! Set Apriori to predetermined cirrus value if USE_CIRRUS_FLAG = Yes 
  !------------------------------------------------------------------------
  if (Pass_Idx == Pass_Idx_Max .and. USE_CIRRUS_FLAG == symbol%YES .and. &
      Temperature_Cirrus(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then
      Tc_Ap = Temperature_Cirrus(Elem_Idx,Line_Idx)
  endif

   if (lun_diag > 0) then 
     write(unit=lun_diag,fmt=*) "Cirrus Temperature = ", Temperature_Cirrus(Elem_Idx,Line_Idx)
   endif 

  !------------------------------------------------------------------------
  ! fill x_ap vector with a priori values  
  !------------------------------------------------------------------------
  Tsfc_Est = Input%Surface_Temperature(Elem_Idx,Line_Idx)

  if (Cloud_Type == symbol%OVERLAP_TYPE) then
    Ts_Ap = Temperature_Lower_Cloud_Apriori(Elem_Idx,Line_Idx)
    Ts_Ap_Uncer = Ts_Ap_Uncer_Lower_Cld
  else
    Ts_Ap = Tsfc_Est
    Ts_Ap_Uncer = Ts_Ap_Uncer_Sfc
  endif

  !------------------------------------------------------------------------
  ! fill x_ap vector with a priori values  
  !------------------------------------------------------------------------
   x_Ap(1) = Tc_Ap
   x_Ap(2) = Ec_Ap
   x_Ap(3) = Beta_Ap
   x_Ap(4) = Ts_Ap

   if (lun_diag > 0) then 
     write(unit=lun_diag,fmt=*) "final x_ap = ", x_ap
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

  !--- now compute Sa
  Sa = 0.0
  Sa(1,1) = Tc_Ap_Uncer
  Sa(2,2) = Ec_Ap_Uncer
  Sa(3,3) = Beta_Ap_Uncer
  Sa(4,4) = Ts_Ap_Uncer

   if (lun_diag > 0) then 
     write(unit=lun_diag,fmt=*) "initial Sa = ", Tc_Ap_Uncer,Ec_Ap_Uncer, Beta_Ap_Uncer, Ts_Ap_Uncer
   endif 

  !--- modify a priori values based on lrc
  if (Pass_Idx /= Pass_Idx_Max .or. USE_CIRRUS_FLAG == symbol%NO) then
    if ((ilrc /= MISSING_VALUE_INTEGER4) .and. &
        (jlrc /= MISSING_VALUE_INTEGER4)) then
         if ((Output%Tc(ilrc,jlrc) /= MISSING_VALUE_REAL4) .and. &
            (Output%Ec(ilrc,jlrc) > 0.00) .and. &
            (Output%Ec(ilrc,jlrc) <= 1.0)) then
          !-- use lrc value but weight uncertainty
          x_Ap(1) = Output%Tc(ilrc,jlrc)
          Sa(1,1) = 5.0 + (1.0-Output%Ec(ilrc,jlrc))*Tc_Ap_Uncer
        endif
    endif
  endif

  !--- square the individual elements to convert to variances (not a matmul)
  Sa = Sa**2

  if (lun_diag > 0) then 
    write(unit=lun_diag,fmt=*) "final Sa = ", Tc_Ap_Uncer,Ec_Ap_Uncer, Beta_Ap_Uncer, Ts_Ap_Uncer
  endif

  !------------------------------------------------------------------------
  ! If a sounder value is available for Tc apriori, combine it with 
  ! other value.  Do this for all passes and only cirrus or overlap
  !------------------------------------------------------------------------
  if (associated(Input%Tc_Cirrus_Sounder)) then
    if (Input%Tc_Cirrus_Sounder(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4 .and. &
      (Cloud_Type == symbol%CIRRUS_TYPE .or. Cloud_Type == symbol%OVERLAP_TYPE)) then

      Tc_Ap_Imager = x_Ap(1)  !K
      Sa_Tc_Imager = Sa(1,1)  !K^2
      Tc_Ap_Sounder = Input%Tc_Cirrus_Sounder(Elem_Idx,Line_Idx)  !K
      Sa_Tc_Sounder = 10.0**2    !K^2
      Sa(1,1) =   1.0/(1.0/Sa_Tc_Imager + 1.0/Sa_Tc_Sounder)
      x_Ap(1) = (Tc_Ap_Imager/Sa_Tc_Imager + Tc_Ap_Sounder/Sa_Tc_Sounder) *  Sa(1,1)

    endif
  endif

  if (lun_diag > 0) then
   write(unit=lun_diag,fmt=*) "x_ap after sounder  = ", x_ap
   write(unit=lun_diag,fmt=*) "S_a after sounder  = ", Sa(1,1), Sa(2,2), Sa(3,3), Sa(4,4)
  endif

  !--- compute inverse of Sa matrix
  Singular_Flag =  INVERT_MATRIX(Sa, Sa_Inv, Num_Param)
  if (Singular_Flag == 1) then
    print *, "Cloud Height warning ==> Singular Sa in ACHA", &
           Elem_Idx,Line_Idx, Cloud_Type
    Fail_Flag(Elem_Idx,Line_Idx) = symbol%YES
    exit
  endif

  !--------------------------------------------------
  ! assign surface emissivity for non-overlap type
  !--------------------------------------------------
  Emiss_Sfc_85um = 1.0
  Emiss_Sfc_11um = 1.0
  Emiss_Sfc_12um = 1.0
  Emiss_Sfc_133um = 1.0
  Emiss_Sfc_67um = 1.0
  if (Cloud_Type /= symbol%OVERLAP_TYPE) then
     if (Input%Chan_On_85um == symbol%YES) Emiss_Sfc_85um = Input%Surface_Emissivity_85um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_11um == symbol%YES) Emiss_Sfc_11um = Input%Surface_Emissivity_11um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_12um == symbol%YES) Emiss_Sfc_12um = Input%Surface_Emissivity_12um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_133um == symbol%YES) Emiss_Sfc_133um = Input%Surface_Emissivity_133um(Elem_Idx,Line_Idx)
     if (Input%Chan_On_67um == symbol%YES) Emiss_Sfc_67um = Input%Surface_Emissivity_67um(Elem_Idx,Line_Idx)
  endif

  
!----------------------------------------------------------------
! Determine the level of the highest inversion (0=if none)
!----------------------------------------------------------------
 call DETERMINE_INVERSION_CHARACTERISTICS(symbol%YES, &
                                          symbol%NO,  &
                                          Tropo_Level_RTM,              &
                                          Sfc_Level_RTM,                &
                                          Input%Surface_Air_Temperature(Elem_Idx,Line_Idx),&
                                          Input%Surface_Elevation(Elem_Idx,Line_Idx),      &
                                          Inver_Top_Level_RTM,    &
                                          Inver_Base_Level_RTM,   & 
                                          Inver_Top_Height,       &
                                          Inver_Base_Height,      & 
                                          Inver_Strength)

!-----------------------------------------------------------------
! start of retrieval loop
!-----------------------------------------------------------------
Iter_Idx = 0
Converged_Flag(Elem_Idx,Line_Idx) = symbol%NO
Fail_Flag(Elem_Idx,Line_Idx) = symbol%NO

!---- assign x to the first guess
x = x_Ap

Retrieval_Loop: do

 Iter_Idx = Iter_Idx + 1

  !---------------------------------------------------------------------
  ! estimate clear-sky radiative transfer terms used in forward model
  !---------------------------------------------------------------------
  Tc_Temp = x(1)
  Ts_Temp = x(4)

  call KNOWING_T_COMPUTE_P_Z(Cloud_Type,Pc_temp,Tc_temp,Zc_Temp,Lev_Idx,ierror,NWP_Profile_Inversion_Flag)
  call KNOWING_T_COMPUTE_P_Z(Cloud_Type,Ps_temp,Ts_temp,Zs_Temp,Lev_Idx,ierror,NWP_Profile_Inversion_Flag)

  !--- compute 11um radiative transfer terms
  Rad_Ac_11um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_11um)

  Trans_Ac_11um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_11um)

  Trans_Bc_11um = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_11um)

  if (Trans_Ac_11um > epsilon(Trans_Ac_11um)) then
     Trans_Bc_11um = Trans_Bc_11um / Trans_Ac_11um
  endif

  Rad_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_11um)
 
  Trans_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_11um)
  Bs = PLANCK_RAD_FAST(Input%Chan_Idx_11um,Ts_Temp)

  Rad_Clear_11um = Rad_Atm + Trans_Atm*Emiss_Sfc_11um*Bs   

  !--- compute 12um radiative transfer terms
  if (Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 8 .or. Acha_Mode_Flag == 9) then
     Rad_Ac_12um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_12um)

     Trans_Ac_12um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_12um)

     Trans_Bc_12um = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_12um) 

     if (Trans_Ac_12um > epsilon(Trans_Ac_12um)) then
       Trans_Bc_12um = Trans_Bc_12um / Trans_Ac_12um
     endif

     Rad_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                           Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_12um)
 
     Trans_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_12um)

     Bs = PLANCK_RAD_FAST(Input%Chan_Idx_12um,Ts_Temp)

     Rad_Clear_12um = Rad_Atm + Trans_Atm*Emiss_Sfc_12um*Bs

  endif

  !--- 13.3um clear radiative transfer terms
  if (Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 7 .or. Acha_Mode_Flag == 8 .or. Acha_Mode_Flag == 9) then
     Rad_Ac_133um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_133um)

     Trans_Ac_133um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_133um)

     Trans_Bc_133um = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_133um) 

     if (Trans_Ac_133um > epsilon(Trans_Ac_133um)) then
       Trans_Bc_133um = Trans_Bc_133um / Trans_Ac_133um
     endif

     Rad_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_133um)
 
     Trans_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_133um)

     Bs = PLANCK_RAD_FAST(Input%Chan_Idx_133um,Ts_Temp)

     Rad_Clear_133um = Rad_Atm + Trans_Atm*Emiss_Sfc_133um*Bs


  endif

  if (Acha_Mode_Flag == 5) then
     Rad_Ac_85um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_85um)

     Trans_Ac_85um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_85um)

     Trans_Bc_85um = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_85um)

     if (Trans_Ac_85um > epsilon(Trans_Ac_85um)) then
       Trans_Bc_85um = Trans_Bc_85um / Trans_Ac_85um
     endif

     Rad_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_85um)
 
     Trans_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_85um)

     Bs = PLANCK_RAD_FAST(Input%Chan_Idx_85um,Ts_Temp)

     Rad_Clear_85um = Rad_Atm + Trans_Atm*Emiss_Sfc_85um*Bs


  endif

  if (Acha_Mode_Flag == 2 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 7) then
     Rad_Ac_67um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_67um)

     Trans_Ac_67um = GENERIC_PROFILE_INTERPOLATION(Zc_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_67um)

     Trans_Bc_67um = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                            Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_67um)

     if (Trans_Ac_67um > epsilon(Trans_Ac_67um)) then
       Trans_Bc_67um = Trans_Bc_67um / Trans_Ac_67um
     endif

     Rad_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Rad_Prof_67um)
 
     Trans_Atm = GENERIC_PROFILE_INTERPOLATION(Zs_Temp, &
                             Hght_Prof_RTM,ACHA_RTM_NWP%Atm_Trans_Prof_67um)

     Bs = PLANCK_RAD_FAST(Input%Chan_Idx_67um,Ts_Temp)

     Rad_Clear_67um = Rad_Atm + Trans_Atm*Emiss_Sfc_67um*Bs


  endif

  !--------------------------------------------------
  ! call forward models
  !--------------------------------------------------
  call COMPUTE_FORWARD_MODEL_AND_KERNEL(Acha_Mode_Flag,  &
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
           Rad_Clear_67um, Rad_Ac_67um, Trans_Ac_67um, Trans_Bc_67um, &
           Rad_Clear_85um, Rad_Ac_85um, Trans_Ac_85um, Trans_Bc_85um, &
           Rad_Clear_11um, Rad_Ac_11um, Trans_Ac_11um, Trans_Bc_11um, &
           Rad_Clear_12um, Rad_Ac_12um, Trans_Ac_12um, Trans_Bc_12um, &
           Rad_Clear_133um, Rad_Ac_133um, Trans_Ac_133um, Trans_Bc_133um, &
           a_Beta_11um_133um_fit, b_Beta_11um_133um_fit, &
           a_Beta_11um_85um_fit, b_Beta_11um_85um_fit, &
           a_Beta_11um_67um_fit, b_Beta_11um_67um_fit, &
           f, K,Emiss_Vector)

  !--------------------------------------------------
  ! compute the Sy convariance matrix
  !--------------------------------------------------
  Call COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE( &
                                                 Emiss_Vector, &
                                                 Acha_Mode_Flag, &
                                                 y_variance, &
                                                 Sy) 

  !--------------------------------------------------
  ! call OE routine to advance the Iteration
  !--------------------------------------------------
  call OPTIMAL_ESTIMATION(Iter_Idx,Iter_Idx_Max,Num_Param,Num_Obs, &
                         Convergence_Criteria,Delta_X_Max, &
                         y,f,x,x_Ap,K,Sy,Sa_inv, &
                         Sx,AKM,Delta_x, &
                         Output%Cost(Elem_Idx,Line_Idx), &
                         Converged_Flag(Elem_Idx,Line_Idx),Fail_Flag(Elem_Idx,Line_Idx))

  !--- check for a failed Iteration
  if (Fail_Flag(Elem_Idx,Line_Idx) == symbol%YES) then
     exit
  endif

  !---------------------------------------------------------
  ! update retrieved Output%Vector
  !---------------------------------------------------------
  x = x + Delta_X

  if (lun_diag > 0) then
     write(unit=lun_diag,fmt=*) "Iter_Idx = ", Iter_Idx
     write(unit=lun_diag,fmt=*) "Kernel Matrix"
     write(unit=lun_diag,fmt=*) "K(1,:) = ", K(1,:)
     write(unit=lun_diag,fmt=*) "K(2,:) = ", K(2,:)
!    write(unit=lun_diag,fmt=*) "K(3,:) = ", K(3,:)
     write(unit=lun_diag,fmt=*) "f = ", f
     write(unit=lun_diag,fmt=*) "Emiss_Vector = ",Emiss_Vector
     write(unit=lun_diag,fmt=*) "Delta_x = ",Delta_x
     write(unit=lun_diag,fmt=*) "x = ",x
 endif


  !--------------------------------------------------------
  ! exit retrieval loop if converged
  !--------------------------------------------------------
  if (Converged_Flag(Elem_Idx,Line_Idx) == symbol%YES) then
       if (lun_diag > 0) then  
             write(unit=lun_diag,fmt=*) "convergence achieved ", Iter_Idx
       endif
       exit
  endif

  !-------------------------------------------------------
  ! constrain to reasonable values
  !-------------------------------------------------------
  x(1) = max(min_allowable_Tc,min(Tsfc_Est+5,x(1)))   
  x(2) = max(0.0,min(x(2),1.0))
  x(3) = max(0.8,min(x(3),1.8))
  x(4) = max(min_allowable_Tc,min(Tsfc_Est+10,x(4)))    

  if (lun_diag > 0) then
     write(unit=lun_diag,fmt=*) "constrained x = ",x
 endif

end do Retrieval_Loop

!=================================================================
! Begin Retrieval Post Processing
!=================================================================
 if (lun_diag > 0) then 
    write(unit=lun_diag,fmt=*) "final f = ", f
    write(unit=lun_diag,fmt=*) "final x = ", x
     write(unit=lun_diag,fmt=*) "==========================================================="
     write(unit=lun_diag,fmt=*) "END OF OE RETRIEVAL"
     write(unit=lun_diag,fmt=*) "==========================================================="
 endif 

!-----------------------------------------------------------------
! Successful Retrieval Post Processing
!-----------------------------------------------------------------
if (Fail_Flag(Elem_Idx,Line_Idx) == symbol%NO) then  !successful retrieval if statement

 !--- save retrievals into the output variables
 Output%Tc(Elem_Idx,Line_Idx) = x(1)
 Output%Ec(Elem_Idx,Line_Idx) = x(2)   !note, this is slant
 Output%Beta(Elem_Idx,Line_Idx) = x(3)

 if (Cloud_Type == symbol%OVERLAP_TYPE) then
    Output%Lower_Tc(Elem_Idx,Line_Idx) = x(4)
    call KNOWING_T_COMPUTE_P_Z(Cloud_Type,Output%Lower_Pc(Elem_Idx,Line_Idx), &
                               Output%Lower_Tc(Elem_Idx,Line_Idx), &
                               Output%Lower_Zc(Elem_Idx,Line_Idx), &
                               Lev_Idx,ierror,NWP_Profile_Inversion_Flag)
 endif


 !--- save uncertainty estimates
 Output%Tc_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(1,1))
 Output%Ec_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(2,2))
 Output%Beta_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(3,3))
 Output%Lower_Tc_Uncertainty(Elem_Idx,Line_Idx) = sqrt(Sx(4,4))   

 !-------------------------------------------------------------------------- 
 !--  If Lower Cloud is placed at surface - assume this single layer
 !--------------------------------------------------------------------------
 if (MULTI_LAYER_LOGIC_FLAG == 0 .or. MULTI_LAYER_LOGIC_FLAG == 2) then 
   if (Output%Lower_Zc(Elem_Idx,Line_Idx) < 1000.0 .and.  &
      Output%Cloud_Type(Elem_Idx,Line_Idx) == symbol%OVERLAP_TYPE) then
      Output%Lower_Zc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Output%Lower_Tc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Output%Lower_Pc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
      Output%Cloud_Type(Elem_Idx,Line_Idx) = symbol%CIRRUS_TYPE
      Output%Lower_Tc_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
   endif  
 endif

 !--- set quality flag for a successful retrieval
 Output%Qf(Elem_Idx,Line_Idx) = 3

 !--- Estimate height and pressure
 call KNOWING_T_COMPUTE_P_Z(Cloud_Type,Output%Pc(Elem_Idx,Line_Idx), &
                            Output%Tc(Elem_Idx,Line_Idx), &
                            Output%Zc(Elem_Idx,Line_Idx),&
                            Lev_Idx,ierror,NWP_Profile_Inversion_Flag)

 !--- check for NWP profile inversion and set meta data flag.
 if (NWP_Profile_Inversion_Flag == 1) then
     Meta_Data_Flags(8) = symbol%YES
 endif

 if (lun_diag > 0) then 
    write(unit=lun_diag,fmt=*) "Results after direct profile inversion "
    write(unit=lun_diag,fmt=*) "Cloud Temperature =",  Output%Tc(Elem_Idx,Line_Idx)
    write(unit=lun_diag,fmt=*) "Cloud Pressure =",  Output%Pc(Elem_Idx,Line_Idx)
    write(unit=lun_diag,fmt=*) "Cloud Height =",  Output%Zc(Elem_Idx,Line_Idx)
    write(unit=lun_diag,fmt=*) "NWP_Profile_Inversion_Flag = ", NWP_Profile_Inversion_Flag
 endif 
  
 !---  for low clouds over water, force fixed lapse rate estimate of height
 Delta_Cld_Temp_Sfc_Temp = Input%Surface_Temperature(Elem_Idx,Line_Idx) - Output%Tc(Elem_Idx,Line_Idx)
 Lapse_Rate = Missing_Value_Real4

 if (Input%Snow_Class (Elem_Idx,Line_Idx) == symbol%NO_SNOW .and. &
      ((Cloud_Type == symbol%WATER_TYPE) .or. &
       (Cloud_Type == symbol%FOG_TYPE) .or. & 
       (Cloud_Type == symbol%SUPERCOOLED_TYPE))) then

     if (Delta_Cld_Temp_Sfc_Temp <  MAX_DELTA_T_INVERSION) then

       !-- select lapse rate  (k/km)
       if (Input%Surface_Type(Elem_Idx,Line_Idx) == symbol%WATER_SFC) then 
           Lapse_Rate = EMPIRICAL_LAPSE_RATE(Input%Surface_Temperature(Elem_Idx,Line_Idx), &
                                       Output%Tc(Elem_Idx,Line_Idx), 0)
       else
           Lapse_Rate = EMPIRICAL_LAPSE_RATE(Input%Surface_Temperature(Elem_Idx,Line_Idx), &
                                       Output%Tc(Elem_Idx,Line_Idx), 1)
       endif

       !--- constrain lapse rate to be with -2 and -10 K/km
       Lapse_Rate = min(-2.0,max(-10.0,Lapse_Rate))

       !--- convert lapse rate to K/m
       Lapse_Rate = Lapse_Rate / 1000.0  !(K/m)

       !-- compute height
       Output%Zc(Elem_Idx,Line_Idx) = -1.0*Delta_Cld_Temp_Sfc_Temp/Lapse_Rate + Input%Surface_Elevation(Elem_Idx,Line_Idx)

       !--- Some negative cloud heights are observed because of bad height
       !--- NWP profiles.
       if (Output%Zc(Elem_Idx,Line_Idx) < 0) then
         Output%Zc(Elem_Idx,Line_Idx) = ZC_FLOOR
       endif

       !--- compute pressure
       call KNOWING_Z_COMPUTE_T_P(Output%Pc(Elem_Idx,Line_Idx),R4_Dummy,Output%Zc(Elem_Idx,Line_Idx),Lev_Idx)

       !--- set meta data flag
       Meta_Data_Flags(7) = symbol%YES
       Output%Inversion_Flag(Elem_Idx,Line_Idx) = 1

       endif
 endif

 if (lun_diag > 0) then 
    write(unit=lun_diag,fmt=*) "Results after oceanic inversion logic "
    write(unit=lun_diag,fmt=*) "Ocean Lapse Rate Input ",                            &
       MAX_DELTA_T_INVERSION, Input%Surface_Temperature(Elem_Idx,Line_Idx), & 
       Output%Tc(Elem_Idx,Line_Idx)
    write(unit=lun_diag,fmt=*) "Ocean Lapse Rate = ", Lapse_Rate
    write(unit=lun_diag,fmt=*) "Cloud Temperature =",  Output%Tc(Elem_Idx,Line_Idx)
    write(unit=lun_diag,fmt=*) "Cloud Pressure =",  Output%Pc(Elem_Idx,Line_Idx)
    write(unit=lun_diag,fmt=*) "Cloud Height =",  Output%Zc(Elem_Idx,Line_Idx)
    write(unit=lun_diag,fmt=*) "Meta_Data_Flags =",  Meta_Data_Flags
 endif 
 
 !-----------------------------------------------------------------------------
 !--- compute height and pressure uncertainties 
 !-----------------------------------------------------------------------------

 !-- Compute Height Uncertainty
 Output%Zc_Uncertainty(Elem_Idx,Line_Idx) = Output%Tc_Uncertainty(Elem_Idx,Line_Idx) /  &
                                            ABS_LAPSE_RATE_DT_DZ_UNCER

 if (Output%Lower_Tc_Uncertainty(Elem_Idx,Line_Idx) /= MISSING_VALUE_REAL4) then
      Output%Lower_Zc_Uncertainty(Elem_Idx,Line_Idx) = Output%Lower_Tc_Uncertainty(Elem_Idx,Line_Idx) /  &
                                                       ABS_LAPSE_RATE_DT_DZ_UNCER
 endif

 !-- Compute Pressure Uncertainty
 Output%Pc_Uncertainty(Elem_Idx,Line_Idx) = Output%Zc_Uncertainty(Elem_Idx,Line_Idx) *  &
                                            ABS_LAPSE_RATE_DlnP_DZ_UNCER * Output%Pc(Elem_Idx,Line_Idx)

 Output%LOWER_Pc_Uncertainty(Elem_Idx,Line_Idx) = Output%LOWER_Zc_Uncertainty(Elem_Idx,Line_Idx) *  &
                                            ABS_LAPSE_RATE_DlnP_DZ_UNCER * Output%LOWER_Pc(Elem_Idx,Line_Idx)

 !-----------------------------------------------------------------------------
 !--- quality flags of the retrieved parameters
 !-----------------------------------------------------------------------------
 do Param_Idx = 1,Num_Param    !loop over parameters
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
 Output%Tc(Elem_Idx,Line_Idx) = x_Ap(1)   !MISSING_VALUE_REAL4
 Output%Ec(Elem_Idx,Line_Idx) = x_Ap(2)   !MISSING_VALUE_REAL4
 Output%Beta(Elem_Idx,Line_Idx) = x_Ap(3) !MISSING_VALUE_REAL4

 if (Cloud_Type == symbol%OVERLAP_TYPE) then
    Output%Lower_Tc(Elem_Idx,Line_Idx) = x_Ap(4) !MISSING_VALUE_REAL4
    call KNOWING_T_COMPUTE_P_Z(Cloud_Type,Output%Lower_Pc(Elem_Idx,Line_Idx),&
                               Output%Lower_Tc(Elem_Idx,Line_Idx),&
                               Output%Lower_Zc(Elem_Idx,Line_Idx),&
                               Lev_Idx,ierror,NWP_Profile_Inversion_Flag)
 endif 

 !--- set derived parameters to missing
 Output%Tau(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
 Output%Reff(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
 Output%Pc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
 Output%Zc(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
 Output%Zc_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4
 Output%Pc_Uncertainty(Elem_Idx,Line_Idx) = MISSING_VALUE_REAL4

 !--- set quality flags
 Output%OE_Qf(:,Elem_Idx,Line_Idx) = 0
 Output%Qf(Elem_Idx,Line_Idx) = 1

 !--- estimate height and pressure
 call KNOWING_T_COMPUTE_P_Z(Cloud_Type,Output%Pc(Elem_Idx,Line_Idx), &
                            Output%Tc(Elem_Idx,Line_Idx), &
                            Output%Zc(Elem_Idx,Line_Idx), &
                            Lev_Idx,ierror,NWP_Profile_Inversion_Flag)

endif                              !end successful retrieval if statement

!--- if retrieval done for an undetected pixel, label the Output%Qf
if (Undetected_Cloud == symbol%YES) then
 Output%Qf(Elem_Idx,Line_Idx) = 2
endif

!-----------------------------------------------------------------
! End Retrieval Post Processing
!-----------------------------------------------------------------

endif     ! ---------- end of data check
 
 
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
                                                 int((2**(i-1)) * Meta_Data_Flags(i),kind=int1)
 enddo


 !-------------------------------------------------------------------------
 !--- spectral cloud emissivity
 !-------------------------------------------------------------------------
  call KNOWING_P_COMPUTE_T_Z(Output%Pc(Elem_Idx,Line_Idx),Tc_Temp,Zc_Temp,Lev_Idx)

  if (Input%Chan_On_67um == symbol%YES) then
        Output%Ec_67um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_67um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_67um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_67um)
  endif
  if (Input%Chan_On_85um == symbol%YES) then
        Output%Ec_85um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_85um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_85um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_85um)
  endif

  if (Input%Chan_On_11um == symbol%YES) then
        Output%Ec_11um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_11um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_11um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_11um)
  endif

  if (Input%Chan_On_12um == symbol%YES) then
        Output%Ec_12um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_12um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_12um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_12um)
  endif

  if (Input%Chan_On_133um == symbol%YES) then
        Output%Ec_133um(Elem_Idx,Line_Idx) = COMPUTE_REFERENCE_LEVEL_EMISSIVITY( Lev_Idx, &
                             Input%Rad_133um(Elem_Idx,Line_Idx), &
                             Input%Rad_Clear_133um(Elem_Idx,Line_Idx), &
                             ACHA_RTM_NWP%Black_Body_Rad_Prof_133um)
  endif


 !---- null profile pointers each time 
 call NULL_PIX_POINTERS(Input, ACHA_RTM_NWP)

 !---close diagnostic output
 if (lun_diag > 0) then
       close(unit=lun_diag)
 endif

 end do Element_Loop

end do Line_Loop

!---------------------------------------------------------------------------
! if selected, compute a background cirrus temperature and use for last pass
!---------------------------------------------------------------------------
if (USE_CIRRUS_FLAG == symbol%YES .and. Pass_Idx == Pass_Idx_Max - 1) then

    call COMPUTE_TEMPERATURE_CIRRUS( &
                 Output%Cloud_Type,   &
                 Output%Tc,          &
                 Output%Ec,          &
                 EMISSIVITY_MIN_CIRRUS, &
                 COUNT_MIN_CIRRUS,      &
                 Box_Half_Width_CIRRUS,      &
                 MISSING_VALUE_REAL4, &
                 Temperature_Cirrus)

endif

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

!------------------------------------------------------------------------
! clean-up and prepare for exit
!------------------------------------------------------------------------
  !--- deallocate 2D arrays
  if (allocated(Elem_Idx_LRC)) deallocate(Elem_Idx_LRC)
  if (allocated(Line_Idx_LRC)) deallocate(Line_Idx_LRC)
  if (allocated(Skip_LRC_Mask)) deallocate(Skip_LRC_Mask)
  if (allocated(Temperature_Cirrus)) deallocate(Temperature_Cirrus)
  if (allocated(Pc_Opaque)) deallocate(Pc_Opaque)
  if (allocated(Tc_Opaque)) deallocate(Tc_Opaque)
  if (allocated(Zc_Opaque)) deallocate(Zc_Opaque)
  if (allocated(Temperature_Lower_Cloud_Apriori)) deallocate(Temperature_Lower_Cloud_Apriori)
  if (allocated(Fail_Flag)) deallocate(Fail_Flag)
  if (allocated(Converged_Flag)) deallocate(Converged_Flag)

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

!-----------------------------------------------------------------
! InterpoLate within profiles knowing P to determine T and Z
!-----------------------------------------------------------------
subroutine KNOWING_P_COMPUTE_T_Z(P,T,Z,Lev_Idx)

     real, intent(in):: P
     real, intent(out):: T
     real, intent(out):: Z
     integer, intent(out):: Lev_Idx
     real:: dp
     real:: dt
     real:: dz

     !--- interpoLate pressure profile
     call LOCATE(Press_Prof_RTM,Num_Levels_RTM_Prof,P,Lev_Idx)
     Lev_Idx = max(1,min(Num_Levels_RTM_Prof-1,Lev_Idx))

     dp = Press_Prof_RTM(Lev_Idx+1) - Press_Prof_RTM(Lev_Idx)
     dt = Temp_Prof_RTM(Lev_Idx+1) - Temp_Prof_RTM(Lev_Idx)
     dz = Hght_Prof_RTM(Lev_Idx+1) - Hght_Prof_RTM(Lev_Idx)

     !--- perform interpoLation
       if (dp /= 0.0) then
           T = Temp_Prof_RTM(Lev_Idx) + dt/dp * (P - Press_Prof_RTM(Lev_Idx))
           Z = Hght_Prof_RTM(Lev_Idx) + dz/dp * (P - Press_Prof_RTM(Lev_Idx))
       else
           T = Temp_Prof_RTM(Lev_Idx)
           Z = Hght_Prof_RTM(Lev_Idx)
       endif

       !--- Some negative cloud heights are observed because  of bad height
       !--- NWP profiles.
       if (Z < 0) then
         Z = ZC_FLOOR
       endif

end subroutine KNOWING_P_COMPUTE_T_Z

!-----------------------------------------------------------------
! InterpoLate within profiles knowing Z to determine T and P
!-----------------------------------------------------------------
subroutine KNOWING_Z_COMPUTE_T_P(P,T,Z,Lev_Idx)

     real, intent(in):: Z
     real, intent(out):: T
     real, intent(out):: P
     integer, intent(out):: Lev_Idx
     real:: dp
     real:: dt
     real:: dz

     !--- interpoLate pressure profile
     call LOCATE(Hght_Prof_RTM,Num_Levels_RTM_Prof,Z,Lev_Idx)
     Lev_Idx = max(1,min(Num_Levels_RTM_Prof-1,Lev_Idx))

     dp = Press_Prof_RTM(Lev_Idx+1) - Press_Prof_RTM(Lev_Idx)
     dt = Temp_Prof_RTM(Lev_Idx+1) - Temp_Prof_RTM(Lev_Idx)
     dz = Hght_Prof_RTM(Lev_Idx+1) - Hght_Prof_RTM(Lev_Idx)

     !--- perform interpoLation
     if (dz /= 0.0) then
           T = Temp_Prof_RTM(Lev_Idx) + dt/dz * (Z - Hght_Prof_RTM(Lev_Idx))
           P = Press_Prof_RTM(Lev_Idx) + dp/dz * (Z - Hght_Prof_RTM(Lev_Idx))
     else
           T = Temp_Prof_RTM(Lev_Idx)
           P = Press_Prof_RTM(Lev_Idx)
     endif

end subroutine KNOWING_Z_COMPUTE_T_P

!-----------------------------------------------------------------
! InterpoLate within profiles knowing T to determine P and Z
!-----------------------------------------------------------------
subroutine KNOWING_T_COMPUTE_P_Z(Cloud_Type,P,T,Z,klev,ierr,Level_Within_Inversion_Flag)

     integer (kind=int1), intent(in):: Cloud_Type
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
     Z = MISSING_VALUE_REAL4
     P = MISSING_VALUE_REAL4
     klev = MISSING_VALUE_INTEGER4

     !--- test for existence of a valid solution with troposphere
     kstart = Tropo_Level_RTM
     kend = Sfc_Level_RTM
     Nlevels_Temp = kend - kstart + 1

     !--- check to see if warmer than max, than assume at surface
     if (T > maxval(Temp_Prof_RTM(kstart:kend))) then
         P = Press_Prof_RTM(kend)
         Z = Hght_Prof_RTM(kend)
         klev = kend - 1
         ierr = symbol%NO
         !--- Some negative cloud heights are observed because of bad height
         !--- NWP profiles.
         if (Z < 0) then
           Z = ZC_FLOOR
         endif
         return
     endif

     !--- check to see if colder than min, than assume above tropopause
     !--- and either limit height to tropopause or extrapoLate in stratosphere
     if (T < minval(Temp_Prof_RTM(kstart:kend))) then
         if (ALLOW_STRATOSPHERE_SOLUTION_FLAG == 1 .and. Cloud_Type == symbol%OVERSHOOTING_TYPE) then
           Z = Hght_Prof_RTM(kstart) + (T - Temp_Prof_RTM(kstart)) / Dt_Dz_Strato
           call KNOWING_Z_COMPUTE_T_P(P,R4_Dummy,Z,klev)
         else
           P = Press_Prof_RTM(kstart)
           Z = Hght_Prof_RTM(kstart)
           klev = kstart + 1
         endif
         ierr = symbol%NO
         return
     endif

     !--- if there is an inversion, look below first
     Level_Within_Inversion_Flag = 0
     if (Inver_Top_Level_RTM > 0 .and. Inver_Base_Level_RTM > 0) then
         kstart = Inver_Top_Level_RTM
         kend =  Inver_Base_Level_RTM
         nlevels_temp = kend - kstart + 1
         call LOCATE(Temp_Prof_RTM(kstart:kend),nlevels_temp,T,klev)
         if ((klev > 0) .and. (klev < nlevels_temp -1)) then
              klev = klev + kstart - 1
              Level_Within_Inversion_Flag = 1
         endif
      endif

    !--- if no solution within an inversion
    if (Level_Within_Inversion_Flag == 0) then
        kstart = Tropo_Level_RTM
        kend = Sfc_Level_RTM
        nlevels_temp = kend - kstart + 1
        call LOCATE(Temp_Prof_RTM(kstart:kend),nlevels_temp,T,klev)
        klev = klev + kstart - 1
        klev = max(1,min(Num_Levels_RTM_Prof-1,klev))
    endif

    !--- General Inversion
    dp = Press_Prof_RTM(klev+1) - Press_Prof_RTM(klev)
    dt = Temp_Prof_RTM(klev+1) - Temp_Prof_RTM(klev)
    dz = Hght_Prof_RTM(klev+1) - Hght_Prof_RTM(klev)

    if (dt /= 0.0) then
        P = Press_Prof_RTM(klev) + dp/dt*(T-Temp_Prof_RTM(klev))
        Z = Hght_Prof_RTM(klev) + dz/dt*(T-Temp_Prof_RTM(klev))
    else
        P = Press_Prof_RTM(klev) 
        Z = Hght_Prof_RTM(klev)
    endif

    !--- Some negative cloud heights are observed because of bad height
    !--- NWP profiles.
    if (Z < 0.0) then
      Z = ZC_FLOOR
    endif

end subroutine KNOWING_T_COMPUTE_P_Z

!-----------------------------------------------------------------
! InterpoLate within profiles knowing Z to determine above cloud
! radiative terms used in forward model
!-----------------------------------------------------------------
function GENERIC_PROFILE_INTERPOLATION(X_value,X_Profile,Y_Profile)  &
            result(Y_value) 

     real, intent(in):: X_value 
     real, dimension(:), intent(in):: X_Profile
     real, dimension(:), intent(in):: Y_Profile
     real:: Y_value

     integer:: Lev_Idx
     real:: dx
     integer:: nlevels

     nlevels = size(X_Profile)

     !--- interpoLate pressure profile
     call LOCATE(X_Profile,nlevels,X_value,Lev_Idx)
     Lev_Idx = max(1,min(nlevels-1,Lev_Idx))

     dx = X_Profile(Lev_Idx+1) - X_Profile(Lev_Idx)

     !--- perform interpoLation
     if (dx /= 0.0) then
        Y_value = Y_Profile(Lev_Idx) +  &
                 (X_value - X_Profile(Lev_Idx))  * &
                 (Y_Profile(Lev_Idx+1) - Y_Profile(Lev_Idx)) / dx
     else
          Y_value = Y_Profile(Lev_Idx)
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
! Convergence_Criteria - the convergence criteria
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
subroutine OPTIMAL_ESTIMATION(Iter_Idx,Iter_Idx_Max,nx,ny, &
                              Convergence_Criteria,Delta_X_Max, &
                              y,f,x,x_Ap,K,Sy,Sa_inv, &
                              Sx,AKM,Delta_x,Conv_Test,Converged_Flag,Fail_Flag)
       

  integer, intent(in):: Iter_Idx
  integer, intent(in):: Iter_Idx_Max
  integer, intent(in):: ny
  integer, intent(in):: nx
  real(kind=real4), intent(in):: Convergence_Criteria
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
  real(kind=real4), intent(out):: Conv_Test
  real(kind=real4), dimension(:), intent(out):: Delta_x
  real(kind=real4), dimension(ny,ny):: Sy_inv
  real(kind=real4), dimension(nx,nx):: Sx_inv
  real(kind=real4), dimension(nx):: Delta_x_dir
  real(kind=real4), dimension(nx):: Delta_x_constrained
  real(kind=real4):: Delta_X_distance_constrained
  real(kind=real4):: Delta_X_distance
  integer, intent(out):: Fail_Flag
  integer, intent(out):: Converged_Flag
  integer:: Singular_Flag
  integer:: ix
  integer:: m
  integer:: p

  m = size(Sy,1)
  p = size(Sx,1)

  Converged_Flag = symbol%NO
  Fail_Flag = symbol%NO
  Delta_X = MISSING_VALUE_REAL4
  Sx = MISSING_VALUE_REAL4

  Singular_Flag =  INVERT_MATRIX(Sy, Sy_Inv, m)
  if (Singular_Flag == symbol%YES) then
   print *, "Cloud Height warning ==> Singular Sy in ACHA "
   Fail_Flag = symbol%YES
   Converged_Flag = symbol%NO
   return 
  endif

  !---- compute next step
  AKM = matmul(Transpose(K),matmul(Sy_inv,K))   !step saving
  Sx_inv = Sa_inv + AKM !(Eq.102 Rodgers)
  Singular_Flag =  INVERT_MATRIX(Sx_inv, Sx, p)
  if (Singular_Flag == symbol%YES) then
   print *, "Cloud Height warning ==> Singular Sx in ACHA "
   Converged_Flag = symbol%NO
   Fail_Flag = symbol%YES
   return
  endif
  
  Delta_x = matmul(Sx,(matmul(Transpose(K),matmul(Sy_inv,(y-f))) +  &
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
  Delta_X_distance = sqrt(sum(Delta_X**2))
  if (Delta_X_distance > 0.0) then
     do ix = 1,nx
        Delta_X_dir(ix) = Delta_X(ix) / Delta_X_distance
     enddo
     do ix = 1,nx
        Delta_X_constrained(ix) =  &
             sign(min(Delta_X_Max(ix),abs(Delta_X(ix))) , Delta_X(ix) )
     enddo
     Delta_X_distance_constrained = sqrt(sum(Delta_X_constrained**2))
     do ix = 1,nx
        Delta_X(ix) = Delta_X_dir(ix)*Delta_X_distance_constrained
     enddo
  endif

  !--- check for non-traditional convergence
! if ((abs(Delta_X(1)) < 0.1) .and. (Iter_Idx > 1)) then
!   Converged_Flag = symbol%YES
!   Fail_Flag = symbol%NO
! endif

  !--- check for traditional convergence
  if (Conv_Test < Convergence_Criteria) then
      Converged_Flag = symbol%YES
      Fail_Flag = symbol%NO
  endif

  !--- check for exceeding allowed number of interactions
  if (Iter_Idx > Iter_Idx_Max) then
      Converged_Flag = symbol%NO
      Fail_Flag = symbol%YES
  endif

  end subroutine OPTIMAL_ESTIMATION

  !----------------------------------------------------------------------
  !---
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

  if (Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 8 .or. Acha_Mode_Flag == 9) then
    Sy(2,2) = T11um_12um_Cal_Uncer**2 + ((1.0-Ec)*T11um_12um_Clr_Uncer)**2 + (T11um_12um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 4) then
    Sy(2,2) = T11um_133um_Cal_Uncer**2 + ((1.0-Ec)*T11um_133um_Clr_Uncer)**2 + (T11um_133um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 7) then
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + ((1.0-Ec)*T11um_67um_Clr_Uncer)**2 + (T11um_67um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 7 .or. Acha_Mode_Flag == 8 .or. Acha_Mode_Flag == 9) then
    Sy(3,3) = T11um_133um_Cal_Uncer**2 + ((1.0-Ec)*T11um_133um_Clr_Uncer)**2 + (T11um_133um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 5) then
    Sy(3,3) = T11um_85um_Cal_Uncer**2 + ((1.0-Ec)*T11um_85um_Clr_Uncer)**2 + (T11um_85um_Cld_Uncer)**2
  endif
  if (Acha_Mode_Flag == 2) then
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + ((1.0-Ec)*T11um_67um_Clr_Uncer)**2 + (T11um_67um_Cld_Uncer)**2
  endif


 end subroutine COMPUTE_Sy

 !---------------------------------------------------------------------
 !--- Compute the Forward Model Estimate (f) and its Kernel (df/dx)
 !---------------------------------------------------------------------
 subroutine COMPUTE_FORWARD_MODEL_AND_KERNEL( &
           Acha_Mode_Flag,  &
           Chan_On_67um, Chan_On_85um, Chan_On_11um, Chan_On_12um, Chan_On_133um, &
           Chan_Idx_67um, Chan_Idx_85um, Chan_Idx_11um, Chan_Idx_12um, Chan_Idx_133um, &
           x,    &
           Rad_Clear_67um, Rad_Ac_67um, Trans_Ac_67um, Trans_Bc_67um,    &
           Rad_Clear_85um, Rad_Ac_85um, Trans_Ac_85um, Trans_Bc_85um,    &
           Rad_Clear_11um, Rad_Ac_11um, Trans_Ac_11um, Trans_Bc_11um,    &
           Rad_Clear_12um, Rad_Ac_12um, Trans_Ac_12um, Trans_Bc_12um,    &
           Rad_Clear_133um, Rad_Ac_133um, Trans_Ac_133um, Trans_Bc_133um,&
           a_Beta_11um_133um_fit, b_Beta_11um_133um_fit, &
           a_Beta_11um_85um_fit, b_Beta_11um_85um_fit,   &
           a_Beta_11um_67um_fit, b_Beta_11um_67um_fit,   &
           f,                                            & 
           K, &
           Emiss_Vector)

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
  real(kind=real4), intent(in):: Trans_Bc_67um
  real(kind=real4), intent(in):: Rad_Clear_85um
  real(kind=real4), intent(in):: Rad_Ac_85um
  real(kind=real4), intent(in):: Trans_Ac_85um
  real(kind=real4), intent(in):: Trans_Bc_85um
  real(kind=real4), intent(in):: Rad_Clear_11um
  real(kind=real4), intent(in):: Rad_Ac_11um
  real(kind=real4), intent(in):: Trans_Ac_11um
  real(kind=real4), intent(in):: Trans_Bc_11um
  real(kind=real4), intent(in):: Rad_Clear_12um
  real(kind=real4), intent(in):: Rad_Ac_12um
  real(kind=real4), intent(in):: Trans_Ac_12um
  real(kind=real4), intent(in):: Trans_Bc_12um
  real(kind=real4), intent(in):: Rad_Clear_133um
  real(kind=real4), intent(in):: Rad_Ac_133um
  real(kind=real4), intent(in):: Trans_Ac_133um
  real(kind=real4), intent(in):: Trans_Bc_133um
  real(kind=real4), intent(in):: a_Beta_11um_133um_fit
  real(kind=real4), intent(in):: b_Beta_11um_133um_fit
  real(kind=real4), intent(in):: a_Beta_11um_85um_fit
  real(kind=real4), intent(in):: b_Beta_11um_85um_fit
  real(kind=real4), intent(in):: a_Beta_11um_67um_fit
  real(kind=real4), intent(in):: b_Beta_11um_67um_fit
  real(kind=real4), dimension(:), intent(out):: f 
  real(kind=real4), dimension(:,:), intent(out):: K
  real(kind=real4), dimension(:), intent(out):: Emiss_Vector

  real(kind=real4):: Tc
  real(kind=real4):: Ts
  real(kind=real4):: Bc_67um
  real(kind=real4):: Bc_85um
  real(kind=real4):: Bc_11um
  real(kind=real4):: Bc_12um
  real(kind=real4):: Bc_133um
  real(kind=real4):: Bs_67um
  real(kind=real4):: Bs_85um
  real(kind=real4):: Bs_11um
  real(kind=real4):: Bs_12um
  real(kind=real4):: Bs_133um
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
  real(kind=real4):: dB_dTs_67um
  real(kind=real4):: dB_dTs_85um
  real(kind=real4):: dB_dTs_11um
  real(kind=real4):: dB_dTs_12um
  real(kind=real4):: dB_dTs_133um


  !---  for notational convenience, rename elements of x to local variables
  Tc = x(1)
  Emiss_11um = min(x(2),0.999999)    !values must be below unity
  Beta_11um_12um = x(3)
  Ts = x(4)

  !--- compute planck Emission for cloud temperature
  if (Chan_On_67um == symbol%YES) Bc_67um = PLANCK_RAD_FAST( Chan_Idx_67um, Tc, dB_dT = dB_dTc_67um)
  if (Chan_On_85um == symbol%YES) Bc_85um = PLANCK_RAD_FAST( Chan_Idx_85um, Tc, dB_dT = dB_dTc_85um)
  if (Chan_On_11um == symbol%YES) Bc_11um = PLANCK_RAD_FAST( Chan_Idx_11um, Tc, dB_dT = dB_dTc_11um)
  if (Chan_On_12um == symbol%YES) Bc_12um = PLANCK_RAD_FAST( Chan_Idx_12um, Tc, dB_dT = dB_dTc_12um)
  if (Chan_On_133um == symbol%YES) Bc_133um = PLANCK_RAD_FAST( Chan_Idx_133um, Tc, dB_dT = dB_dTc_133um)

  !--- compute planck Emission for surface temperature
  if (Chan_On_67um == symbol%YES) Bs_67um = PLANCK_RAD_FAST( Chan_Idx_67um, Ts, dB_dT = dB_dTs_67um)
  if (Chan_On_85um == symbol%YES) Bs_85um = PLANCK_RAD_FAST( Chan_Idx_85um, Ts, dB_dT = dB_dTs_85um)
  if (Chan_On_11um == symbol%YES) Bs_11um = PLANCK_RAD_FAST( Chan_Idx_11um, Ts, dB_dT = dB_dTs_11um)
  if (Chan_On_12um == symbol%YES) Bs_12um = PLANCK_RAD_FAST( Chan_Idx_12um, Ts, dB_dT = dB_dTs_12um)
  if (Chan_On_133um == symbol%YES) Bs_133um = PLANCK_RAD_FAST( Chan_Idx_133um, Ts, dB_dT = dB_dTs_133um)

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
 K(1,4) = ((1.0 - Emiss_11um) * Trans_Ac_11um * Trans_Bc_11um * dB_dTs_11um) / dB_dT_11um * Emiss_Sfc_11um    !dT_11um / dT_s

 !--- forward model for 11um - 12um
 if (Acha_Mode_Flag == 3 .or. Acha_Mode_Flag == 8 .or. Acha_Mode_Flag == 5 .or. Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 9) then
   Rad_12um = Emiss_12um*Rad_Ac_12um + Trans_Ac_12um * Emiss_12um * Bc_12um +  Trans_12um * Rad_Clear_12um
   f(2) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_12um,Rad_12um,dB_dT = dB_dT_12um)
   K(2,1) = K(1,1) - Trans_Ac_12um * Emiss_12um * dB_dTc_12um / dB_dT_12um   
   K(2,2) = K(1,2) -  (Rad_Ac_12um + Trans_Ac_12um*Bc_12um-Rad_Clear_12um)*&
                      (dEmiss_12um_dEmiss_11um)/dB_dT_12um  
   K(2,3) = (Rad_Ac_12um+Trans_Ac_12um*Bc_12um-Rad_Clear_12um)/ &
            dB_dT_12um*alog(1.0-Emiss_11um)*(1.0-Emiss_12um)    
   K(2,4) = K(1,4) -  (Trans_Ac_12um * Trans_12um * Trans_Bc_12um * dB_dTs_12um) / dB_dT_12um * Emiss_Sfc_12um
 endif
 !--- forward model for 11um - 133um
 if (Acha_Mode_Flag == 4 .or. Acha_Mode_Flag == 7) then
   Rad_133um = Emiss_133um*Rad_Ac_133um + Trans_Ac_133um * Emiss_133um * Bc_133um +  Trans_133um * Rad_Clear_133um
   f(2) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_133um,Rad_133um,dB_dT = dB_dT_133um)
   K(2,1) = K(1,1) - Trans_Ac_133um*Emiss_133um*dB_dTc_133um/dB_dT_133um
   K(2,2) = K(1,2) - (Rad_Ac_133um + Trans_ac_133um*Bc_133um - Rad_Clear_133um)*&
                     (dEmiss_133um_dEmiss_11um)/dB_dT_133um
   K(2,3) = (Rad_Ac_133um + Trans_ac_133um*Bc_133um -Rad_Clear_133um)/ &
             dB_dT_133um * alog(1.0-Emiss_11um)*(1.0-Emiss_133um)
   K(2,4) = K(1,4) -  (Trans_Ac_133um * Trans_133um * Trans_Bc_133um * dB_dTs_133um) / dB_dT_133um * Emiss_Sfc_133um
 endif
 if (Acha_Mode_Flag == 8 .or. Acha_Mode_Flag == 9) then
   Rad_133um = Emiss_133um*Rad_Ac_133um + Trans_Ac_133um * Emiss_133um * Bc_133um +  Trans_133um * Rad_Clear_133um
   f(3) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_133um,Rad_133um,dB_dT = dB_dT_133um)
   K(3,1) = K(1,1) - Trans_Ac_133um*Emiss_133um*dB_dTc_133um/dB_dT_133um
   K(3,2) = K(1,2) - (Rad_Ac_133um + Trans_ac_133um*Bc_133um - Rad_Clear_133um)*&
                     (dEmiss_133um_dEmiss_11um)/dB_dT_133um
   K(3,3) = (Rad_Ac_133um + Trans_ac_133um*Bc_133um -Rad_Clear_133um)/ &
             dB_dT_133um * alog(1.0-Emiss_11um)*(1.0-Emiss_133um)
   K(3,4) = K(1,4) -  (Trans_Ac_133um * Trans_133um * Trans_Bc_133um * dB_dTs_133um) / dB_dT_133um * Emiss_Sfc_133um
 endif
 if (Acha_Mode_Flag == 5) then
   Rad_85um = Emiss_85um*Rad_Ac_85um + Trans_Ac_85um * Emiss_85um * Bc_85um +  Trans_85um * Rad_Clear_85um
   f(3) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_85um,Rad_85um,dB_dT = dB_dT_85um)
   K(3,1) = K(1,1) - Trans_Ac_85um*Emiss_85um*dB_dTc_85um/dB_dT_85um
   K(3,2) = K(1,2) - (Rad_Ac_85um + Trans_ac_85um*Bc_85um - Rad_Clear_85um)*&
                     (dEmiss_85um_dEmiss_11um)/dB_dT_85um
   K(3,3) = (Rad_Ac_85um + Trans_ac_85um*Bc_85um -Rad_Clear_85um)/ &
             dB_dT_85um * alog(1.0-Emiss_11um)*(1.0-Emiss_85um)
   K(3,4) = K(1,4) -  (Trans_Ac_85um * Trans_85um * Trans_Bc_85um * dB_dTs_85um) / dB_dT_85um * Emiss_Sfc_85um
 endif
 if (Acha_Mode_Flag == 6 .or. Acha_Mode_Flag == 7) then
   Rad_67um = Emiss_67um*Rad_Ac_67um + Trans_Ac_67um * Emiss_67um * Bc_67um + Trans_67um * Rad_Clear_67um
   f(3) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_67um,Rad_67um,dB_dT = dB_dT_67um)
   K(3,1) = K(1,1) - Trans_Ac_67um*Emiss_67um*dB_dTc_67um/dB_dT_67um
   K(3,2) = K(1,2) - (Rad_Ac_67um + Trans_ac_67um*Bc_67um - Rad_Clear_67um)*&
                     (dEmiss_67um_dEmiss_11um)/dB_dT_67um
   K(3,3) = (Rad_Ac_67um + Trans_ac_67um*Bc_67um - Rad_Clear_67um)/ &
             dB_dT_67um * alog(1.0-Emiss_11um)*(1.0-Emiss_67um)
   K(3,4) = K(1,4) -  (Trans_Ac_67um * Trans_67um * Trans_Bc_67um * dB_dTs_67um) / dB_dT_67um * Emiss_Sfc_67um
 endif
 if (Acha_Mode_Flag == 2) then
   Rad_67um = Emiss_67um*Rad_Ac_67um + Trans_Ac_67um * Emiss_67um * Bc_67um + Trans_67um * Rad_Clear_67um
   f(2) = f(1) - PLANCK_TEMP_FAST(Chan_Idx_67um,Rad_67um,dB_dT = dB_dT_67um)
   K(2,1) = K(1,1) - Trans_Ac_67um*Emiss_67um*dB_dTc_67um/dB_dT_67um
   K(2,2) = K(1,2) - (Rad_Ac_67um + Trans_ac_67um*Bc_67um - Rad_Clear_67um)*&
                     (dEmiss_67um_dEmiss_11um)/dB_dT_67um
   K(2,3) = (Rad_Ac_67um + Trans_ac_67um*Bc_67um - Rad_Clear_67um)/ &
             dB_dT_67um * alog(1.0-Emiss_11um)*(1.0-Emiss_67um)
   K(2,4) = K(1,4) -  (Trans_Ac_67um * Trans_67um * Trans_Bc_67um * dB_dTs_67um) / dB_dT_67um * Emiss_Sfc_67um
 endif

 !--- determine number of channels
  select case(Acha_Mode_Flag)
     case(1)  !avhrr, goes-im
       Emiss_Vector(1) = Emiss_11um
     case(2)  !goes-np 3 chan
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_67um
     case(3)  !avhrr, goes-im
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_12um
     case(4)  !goes-nop
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_133um
     case(5)  !viirs
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_12um
       Emiss_Vector(3) = Emiss_85um
     case(6)  !goes-im 3 chan
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_67um
       Emiss_Vector(3) = Emiss_12um
     case(7)  !goes-np 3 chan
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_67um
       Emiss_Vector(3) = Emiss_133um
     case(8)  !goes-r
       Emiss_Vector(1) = Emiss_11um
       Emiss_Vector(2) = Emiss_12um
       Emiss_Vector(3) = Emiss_133um
  end select

end subroutine COMPUTE_FORWARD_MODEL_AND_KERNEL

!----------------------------------------------------------------------
!---
!----------------------------------------------------------------------
subroutine COMPUTE_APRIORI_BASED_ON_PHASE_ETROPO( &
                           Cloud_Phase, &
                           Emiss_11um_Tropo, &
                           Latitude, &
                           Ttropo, &
                           T11um, &
                           T11um_Lrc, &
                           Tc_Opaque, &
                           Mu, &
                           Snow_Flag, &
                           Tair, &
                           Tc_Ap, &
                           Tc_Ap_Uncer, &
                           Ec_Ap, &
                           Ec_Ap_Uncer, &
                           Beta_Ap,  &
                           Beta_Ap_Uncer)

  integer, intent(in):: Cloud_Phase
  real(kind=real4), intent(in):: Emiss_11um_Tropo
  real(kind=real4), intent(in):: Latitude
  real(kind=real4), intent(in):: Ttropo
  real(kind=real4), intent(in):: T11um
  real(kind=real4), intent(in):: T11um_Lrc
  real(kind=real4), intent(in):: Tc_Opaque
  real(kind=real4), intent(in):: Mu
  integer(kind=int1), intent(in):: Snow_Flag
  real(kind=real4), intent(in):: Tair
  real(kind=real4), intent(out):: Tc_Ap
  real(kind=real4), intent(out):: Ec_Ap
  real(kind=real4), intent(out):: Beta_Ap
  real(kind=real4), intent(out):: Tc_Ap_Uncer
  real(kind=real4), intent(out):: Ec_Ap_Uncer
  real(kind=real4), intent(out):: Beta_Ap_Uncer

  real(kind=real4):: Tc_Ap_Cirrus
  real(kind=real4):: Tc_Ap_Uncer_Cirrus
  real(kind=real4):: Tc_Ap_Opaque
  real(kind=real4):: Emiss_Weight
  real(kind=real4):: Emiss_Weight2

  !--- calipso values (not multiplier on uncer values)
  call compute_cirrus_apriori(Ttropo, Latitude, Tc_Ap_Cirrus, Tc_Ap_Uncer_Cirrus)

  Tc_Ap_Opaque = Tc_Opaque

  if (T11um_Lrc /= MISSING_VALUE_REAL4) then
      Tc_Ap_Opaque = T11um_Lrc
  endif
  if (Tc_Ap_Opaque == MISSING_VALUE_REAL4) then
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

!   if (Emiss_11um_Tropo <= 0.0) then
!           Emiss_Weight = 0.0
!   elseif (Emiss_11um_Tropo > 1.0) then
!           Emiss_Weight = 1.0
!   else
!           Emiss_Weight = Emiss_11um_Tropo
!   endif

!   Emiss_Weight2 = Emiss_Weight

!   Tc_Ap = Emiss_Weight2*Tc_Ap_Opaque + &
!           (1.0-Emiss_Weight2)*Tc_Ap_Cirrus

!   Tc_Ap_Uncer = Emiss_Weight2*Tc_Ap_Uncer_Opaque + &
!           (1.0-Emiss_Weight2)*Tc_Ap_Uncer_Cirrus

!   Ec_Ap = min(0.99,max(0.1,Emiss_Weight)) 
!   Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus

!   Beta_Ap = Beta_Ap_Ice
!   Beta_Ap_Uncer = Beta_Ap_Uncer_Ice

!==============
    Tc_Ap = min(Tc_Ap_Cirrus,Tc_Ap_Opaque)
    Tc_Ap_Uncer = Tc_Ap_Uncer_Cirrus
    Ec_Ap = min(0.99,max(0.1,Emiss_11um_Tropo)) 
    Ec_Ap_Uncer = Ec_Ap_Uncer_Cirrus
    Beta_Ap = Beta_Ap_Ice
    Beta_Ap_Uncer = Beta_Ap_Uncer_Ice

    if (Emiss_11um_Tropo > 0.90) then
      Tc_Ap = Tc_Ap_Opaque
      Tc_Ap_Uncer = Tc_Ap_Opaque
    endif



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
subroutine DETERMINE_SFC_TYPE_FORWARD_MODEL( &
                                         Surface_Type, &
                                         Snow_Class, &
                                         Latitude, &
                                         Ch20_Surface_Emissivity, &
                                         Sfc_Type_Forward_Model)

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
! Compute Sy based on the clear-sky error covariance calcuLations.
! Using Andy's simpler expression
!
! This assumes that 
! Acha_Mode_Flag: 1=11um,2=11+6.7um,3=11+12um,4=11+13.3um,5=8.5+11+12um
!                 6=11+6.7+12um,7=11+6.7+13.3um,8=11+12+13.3um
!
! Input:
! Emiss_Vector = a vector of emissivities in each channel. 
! Acha_Mode_Flag: 1=11um,2=11+6.7um,3=11+12um,4=11+13.3um,5=8.5+11+12um
!                 6=11+6.7+12um,7=11+6.7+13.3um,8=11+12+13.3um
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
 if (Acha_Mode_Flag >= 2) then
    Sub_Pixel_Uncer(2) = max(0.25,y_variance(2)/4.0)
 endif
 if (Acha_Mode_Flag >= 5) then
    Sub_Pixel_Uncer(3) = max(0.25,y_variance(3)/4.0)
 endif

 select case(Acha_Mode_Flag)
  case(1) 
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar

  case(2) 
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_67um_Covar
    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_67um_Btd_11um_67um_Covar

  case(3)
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_12um_Covar
    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_12um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_12um_Btd_11um_12um_Covar

  case(4) 
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_133um_Covar
    Sy(2,1) = Trans2*Bt_11um_Btd_11um_133um_Covar
    Sy(2,2) = T11um_133um_Cal_Uncer**2 +  Sub_Pixel_Uncer(2)+ Trans2*Btd_11um_133um_Btd_11um_133um_Covar

  case(5)
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_12um_Covar
    Sy(1,3) = Trans2*Bt_11um_Btd_11um_85um_Covar

    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_12um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_12um_Btd_11um_12um_Covar
    Sy(2,3) = Trans2*Btd_11um_12um_Btd_11um_85um_Covar

    Sy(3,3) = T11um_85um_Cal_Uncer**2 + Sub_Pixel_Uncer(3) + Trans2*Btd_11um_85um_Btd_11um_85um_Covar
    Sy(3,1) = Sy(1,3)
    Sy(3,2) = Sy(2,3)

  case(6)
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_12um_Covar
    Sy(1,3) = Trans2*Bt_11um_Btd_11um_67um_Covar

    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_67um_Btd_11um_67um_Covar
    Sy(2,3) = Trans2*Btd_11um_12um_Btd_11um_67um_Covar

    Sy(3,3) = T11um_12um_Cal_Uncer**2 + Sub_Pixel_Uncer(3) + Trans2*Btd_11um_12um_Btd_11um_12um_Covar
    Sy(3,1) = Sy(1,3)
    Sy(3,2) = Sy(2,3)

  case(7)
    Sy(1,1) = T11um_Cal_Uncer**2 + Sub_Pixel_Uncer(1) + Trans2*Bt_11um_Bt_11um_Covar
    Sy(1,2) = Trans2*Bt_11um_Btd_11um_133um_Covar
    Sy(1,3) = Trans2*Bt_11um_Btd_11um_67um_Covar

    Sy(2,1) = Sy(1,2)
    Sy(2,2) = T11um_67um_Cal_Uncer**2 + Sub_Pixel_Uncer(2) + Trans2*Btd_11um_67um_Btd_11um_67um_Covar
    Sy(2,3) = Trans2*Btd_11um_67um_Btd_11um_133um_Covar

    Sy(3,3) = T11um_133um_Cal_Uncer**2 + Sub_Pixel_Uncer(3) + Trans2*Btd_11um_133um_Btd_11um_133um_Covar
    Sy(3,1) = Sy(1,3)
    Sy(3,2) = Sy(2,3)

  case(8,9) 
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

  end select

end subroutine COMPUTE_SY_BASED_ON_CLEAR_SKY_COVARIANCE

!----------------------------------------------------------------------
! Compute Sy based on the clear-sky error covariance calcuLations.
! This assumes that 
! Acha_Mode_Flag: 1=11um,2=11+6.7um,3=11+12um,4=11+13.3um,5=8.5+11+12um
!                 6=11+6.7+12um,7=11+6.7+13.3um,8=11+12+13.3um,9-11+12+13.3um
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
subroutine DETERMINE_ACHA_MODE_BASED_ON_CHANNELS( &
                                                 Acha_Mode_Flag, &
                                                 Chan_On_67um, &
                                                 Chan_On_85um, &
                                                 Chan_On_11um, &
                                                 Chan_On_12um, &
                                                 Chan_On_133um)

  integer, intent(INOUT):: Acha_Mode_Flag
  integer, intent(in):: Chan_On_67um
  integer, intent(in):: Chan_On_85um
  integer, intent(in):: Chan_On_11um
  integer, intent(in):: Chan_On_12um
  integer, intent(in):: Chan_On_133um

  if (Acha_Mode_Flag == -1) then
     if (Chan_On_11um == symbol%YES .and. Chan_On_12um == symbol%YES) then
        if (Chan_On_133um == symbol%YES) then
            Acha_Mode_Flag = 8          ! 11/12/13.3 um
        elseif (Chan_On_85um == symbol%YES) then
            Acha_Mode_Flag = 5          ! 8.5/11/12 um
        elseif (Chan_On_67um == symbol%YES) then
            Acha_Mode_Flag = 6          ! 6.7/11/12 um
        else
            Acha_Mode_Flag = 3          ! 11/12 um
        endif
     endif
  endif

  if (Acha_Mode_Flag == -1) then
     if (Chan_On_12um == symbol%NO) then
        if (Chan_On_67um == symbol%YES .and. Chan_On_133um == symbol%YES) then
              Acha_Mode_Flag = 7       !6.7/11/13.3
        endif
        if (Chan_On_67um == symbol%NO .and. Chan_On_133um == symbol%YES) then
              Acha_Mode_Flag = 4       !11/13.3
        endif
        if (Chan_On_67um == symbol%YES .and. Chan_On_133um == symbol%NO) then
              Acha_Mode_Flag = 2       !11/6.7
        endif
     endif
  endif

  !--- if unsuccessful, resort to mode 1
  if (Acha_Mode_Flag == -1) then
     if (Chan_On_11um == symbol%YES) then
              Acha_Mode_Flag = 1       !11
     endif
  endif
   
   
end subroutine  DETERMINE_ACHA_MODE_BASED_ON_CHANNELS

!----------------------------------------------------------------------------
! Function INTERPOLATE_PROFILE_ACHA
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
 function INTERPOLATE_PROFILE_ACHA(z1,z2,z3,z4,lonx,Latx) result(z)

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

 end function INTERPOLATE_PROFILE_ACHA

 !-------------------------------------------------------------------------
 ! Input:
 ! Tropo_Level - level in RTM profiles of the Tropopause
 ! Sfc_Level - level in RTM profiles closest but above the surface
 ! Sfc_Air_Temp = air temperature at the surface level
 ! Sfc_Height = height of the surface level (m)
 ! Inversion_Top_Level -  level in RTM profiles closest to but below the top of
 !                         inversion
 ! Inversion_Base_Level -  level in RTM profiles closest to but above the base of
 !                         inversion
 ! Inversion_Strength - Temperature difference between Top and Base (K)
 ! Inversion_Base_Height -  Height of Inversion Base (m)
 ! Inversion_Top_Height -  Height of Inversion Top (m)
 !
 ! Input - via module-wide variables
 ! Press_Prof_RTM - pressure profile
 ! Hght_Prof_RTM - height profile
 ! Temp_Prof_RTM - temperature profile
 !
 ! Output - via module-wide variables
 ! Inver_Prof_RTM - level flags (0/1) if inversion present
 !--------------------------------------------------------------------------
 subroutine DETERMINE_INVERSION_CHARACTERISTICS(symbol_yes,               &
                                                symbol_no,                &
                                                Tropo_Level,              &
                                                Sfc_Level,                &
                                                Sfc_Air_Temp,             &
                                                Sfc_Height,               &
                                                Top_Lev_Idx,            &
                                                Base_Lev_Idx,           & 
                                                Inversion_Top_Height,     &
                                                Inversion_Base_Height,    & 
                                                Inversion_Strength)
   integer(kind=int1), intent(in) :: symbol_yes
   integer(kind=int1), intent(in) :: symbol_no
   integer, intent(in):: Tropo_Level
   integer, intent(in):: Sfc_Level
   real, intent(in):: Sfc_Air_Temp
   real, intent(in):: Sfc_Height
   real, intent(out):: Inversion_Top_Height
   real, intent(out):: Inversion_Base_Height
   integer, intent(out):: Top_Lev_Idx
   integer, intent(out):: Base_Lev_Idx
   real, intent(out):: Inversion_Strength
   integer:: k

   Inver_Prof_RTM = symbol_NO

   do k = Sfc_Level, Tropo_Level, -1

      if (Press_Prof_RTM(k) >= MIN_P_INVERSION) then

         if (Temp_Prof_RTM(k-1) - Temp_Prof_RTM(k) > DELTA_T_LAYER_INVERSION) then
            Inver_Prof_RTM(k-1:k) = symbol_YES
         endif

      endif
   enddo

   Top_Lev_Idx =  0
   do k = Tropo_Level,Sfc_Level,1
      if (Inver_Prof_RTM(k) == symbol_YES .and. Top_Lev_Idx == 0) then
         Top_Lev_Idx = k
         exit
      endif
   enddo

   Base_Lev_Idx = 0
   do k = Sfc_Level, Tropo_Level, -1
      if (Inver_Prof_RTM(k) == symbol_YES .and. Base_Lev_Idx == 0) then
         Base_Lev_Idx = k
         exit
      endif
   enddo

   Inversion_Strength  = Missing_Value_Real4
   Inversion_Base_Height = Missing_Value_Real4
   Inversion_Top_Height = Missing_Value_Real4

   !---- inversion top height (meters)
   if (Top_Lev_Idx /= 0)  Inversion_Top_Height = Hght_Prof_RTM(Top_Lev_Idx)

   !---- inversion base height (meters)
   if (Base_Lev_Idx /= 0) Inversion_Base_Height = Hght_Prof_RTM(Base_Lev_Idx)

   !--- assume inversion streches to surface if lowest level is the surface level
   if (Base_Lev_Idx == Sfc_Level .and.  Sfc_Height /= Missing_Value_Real4) then
    Inversion_Base_Height = Sfc_Height
   endif

   !--- inversion temperature strength
   if (Base_Lev_Idx /= 0 .and. Top_Lev_Idx /= 0) then
     Inversion_Strength = Temp_Prof_RTM(Top_Lev_Idx) - Temp_Prof_RTM(Base_Lev_Idx)

     !--- assume inversion streches to surface if lowest level is the surface level
     if (Base_Lev_Idx == Sfc_Level .and.  Sfc_Air_Temp /= Missing_Value_Real4) then
        Inversion_Strength = Temp_Prof_RTM(Top_Lev_Idx) - Sfc_Air_Temp
     endif
   endif

 end subroutine DETERMINE_INVERSION_CHARACTERISTICS


 !---------------------------------------------------------------------
 ! Find Opaque Cloud Level - highest Level Inversion below trop
 !---------------------------------------------------------------------
 subroutine DETERMINE_OPAQUE_CLOUD_HEIGHT( &
                                   Radiance_11um, &
                                   Black_Body_Rad_Prof_11um, &
                                   Press_Prof, &
                                   Height_Prof, &
                                   Temp_Prof, &
                                   Tropo_Level, &
                                   Sfc_Level, &
                                   Pc_Opaque, &
                                   Tc_Opaque, &
                                   Zc_Opaque)
                             
   real(kind=real4), intent(in):: Radiance_11um
   real(kind=real4), intent(in), dimension(:):: Black_Body_Rad_Prof_11um
   real(kind=real4), intent(in), dimension(:):: Press_Prof
   real(kind=real4), intent(in), dimension(:):: Height_Prof
   real(kind=real4), intent(in), dimension(:):: Temp_Prof
   integer(kind=int4), intent(in):: Tropo_Level
   integer(kind=int4), intent(in):: Sfc_Level
   real(kind=real4), intent(out):: Pc_Opaque
   real(kind=real4), intent(out):: Tc_Opaque
   real(kind=real4), intent(out):: Zc_Opaque
   integer:: Lev_Idx
   integer:: Lev_Idx_Start
   integer:: Lev_Idx_End

   !--- initialize
   Pc_Opaque =  MISSING_VALUE_REAL4
   Zc_Opaque =  MISSING_VALUE_REAL4
   Tc_Opaque =  MISSING_VALUE_REAL4

   !--- restrict levels to consider
   Lev_Idx_Start = Tropo_Level 
   Lev_Idx_End = Sfc_Level

   !--- loop through levels
   level_loop: do Lev_Idx = Lev_Idx_Start, Lev_Idx_End
      Pc_Opaque = Press_Prof(Lev_Idx-1)
      Zc_Opaque = Height_Prof(Lev_Idx-1)
      Tc_Opaque = Temp_Prof(Lev_Idx-1)
     if (Black_Body_Rad_Prof_11um(Lev_Idx) > Radiance_11um) then
         exit
     endif
   end do Level_Loop

   !--- Some negative cloud heights are observed because of bad height
   !--- NWP profiles.
   if (Zc_Opaque < 0) then
     Zc_Opaque = ZC_FLOOR
   endif

 end subroutine DETERMINE_OPAQUE_CLOUD_HEIGHT
 !----------------------------------------------------------------------
 !
 ! Compute the IR Emissivity at a Reference Level
 !
 ! Ref_Level refers to a level index in the profiles
 ! Toa_Radiance = top of atmosphere radiance
 ! Toa_Radiance_Clear = top of atmosphere radiance under clear-skies
 !
 ! 
 ! Black_Body_Rad_Prof_11um_RTM - this is in memory
 !
 !----------------------------------------------------------------------
 function COMPUTE_REFERENCE_LEVEL_EMISSIVITY(Ref_Level,Toa_Radiance, &
                                             Toa_Radiance_Clear, &
                                             Black_Body_Rad_Prof_11um) &
                                             result(Emissivity_Ref_Level)

    integer(kind=int4), intent(in):: Ref_Level
    real (kind=real4), intent(in):: Toa_Radiance  
    real (kind=real4), intent(in):: Toa_Radiance_Clear
    real (kind=real4), intent(in), dimension(:):: Black_Body_Rad_Prof_11um
    real (kind=real4):: Emissivity_Ref_Level

    Emissivity_Ref_Level = &
    (Toa_Radiance - Toa_Radiance_Clear) / &
    (Black_Body_Rad_Prof_11um(Ref_Level) - Toa_Radiance_Clear)

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
         if (minval(Grad_Array) == MISSING_VALUE_REAL4) exit 

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
    Stddev_of_Array_r8 = MISSING_VALUE_REAL4
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
! Compute Parallax Correction
!
! This routine generates new Lat and Lon arrays that are parallax
! corrected based on the cloud height
!
! Input: Senzen - sensor viewing zenith angle (deg) 
!        Senaz  - sensor azimuth angle (deg)
!        Lat - uncorrected Latitude (deg)
!        Lon  - uncorrected longitude (deg)
!        Zsfc  - surface elevation (m)
!        Zcld  - cloud height (m)
!
! Output
!       Lat_Pc - corrected Latitude
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
   real,parameter:: Lat_Spacing_Per_m = 8.9932e-06   ! ( = 1.0/111000.0 m )

   Num_Elem = size(Zcld,1) 
   Num_Line = size(Zcld,2)

   !--- initialize output to standard values
   Lat_Pc = Lat
   Lon_Pc = Lon

   !--- loop over pixels in segment
   element_loop: do Elem_Idx = 1, Num_Elem
    line_loop: do Line_Idx = 1, Num_Line

     !--- check for valid data
     if (Zcld(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4 .or. &
         Senzen(Elem_Idx,Line_Idx) == MISSING_VALUE_REAL4) cycle

     !--- compute correction
     Total_Displacement = max(0.0,tan(Senzen(Elem_Idx,Line_Idx)*Dtor)* &
                                     (Zcld(Elem_Idx,Line_Idx) - Zsfc(Elem_Idx,Line_Idx)))

     Lon_Spacing_Per_m = Lat_Spacing_Per_m / cos(Lat(Elem_Idx,Line_Idx)*Dtor)

     Delta_Lon = sin(Senaz(Elem_Idx,Line_Idx)*Dtor)*Total_Displacement * Lon_Spacing_Per_m
     Delta_Lat = cos(Senaz(Elem_Idx,Line_Idx)*Dtor)*Total_Displacement * Lat_Spacing_Per_m

     !--- generate output positions
     Lat_Pc(Elem_Idx,Line_Idx) = Lat(Elem_Idx,Line_Idx) + Delta_Lat
     Lon_Pc(Elem_Idx,Line_Idx) = Lon(Elem_Idx,Line_Idx) + Delta_Lon

    enddo line_loop
   enddo element_loop

end subroutine PARALLAX_ACHA 
!----------------------------------------------------------------------
!--- check that the Acha_Mode is consistent with available channels
!--- if consistent, Acha_Mode_Error_Flag = 0, if not, flag = 1
!----------------------------------------------------------------------
subroutine CHECK_ACHA_MODE( &
                           Acha_Mode_Input, &
                           Chan_On_67um, &
                           Chan_On_85um, &
                           Chan_On_11um, &
                           Chan_On_12um, &
                           Chan_On_133um, &
                           Acha_Mode_Error_Flag)

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
       ((Acha_Mode_Input == 2) .or. &
        (Acha_Mode_Input == 6) .or. &
        (Acha_Mode_Input == 7))) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_85um == symbol%NO) .and. (Acha_Mode_Input == 5)) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_12um == symbol%NO) .and. &
       ((Acha_Mode_Input == 3) .or. &
        (Acha_Mode_Input == 5) .or. &
        (Acha_Mode_Input == 6))) then
       Acha_Mode_Error_Flag = 1
       return
   endif

   if ((Chan_On_133um == symbol%NO) .and. &
       ((Acha_Mode_Input == 4) .or. &
        (Acha_Mode_Input == 7) .or. &
        (Acha_Mode_Input == 8))) then
       Acha_Mode_Error_Flag = 1
       return
   endif

end subroutine CHECK_ACHA_MODE

!------------------------------------------------------------------------------
! Null Pixel Level Pointers 
!------------------------------------------------------------------------------
subroutine NULL_PIX_POINTERS(Input, ACHA_RTM_NWP)

   type(acha_input_struct), intent(inout) :: Input
   type(acha_rtm_nwp_struct), intent(inout) :: ACHA_RTM_NWP

   ACHA_RTM_NWP%T_Prof => NULL()
   ACHA_RTM_NWP%T_Prof_1 => NULL() 
   ACHA_RTM_NWP%T_Prof_2 => NULL() 
   ACHA_RTM_NWP%T_Prof_3 => NULL()

   ACHA_RTM_NWP%Z_Prof => NULL() 
   ACHA_RTM_NWP%Z_Prof_1 => NULL() 
   ACHA_RTM_NWP%Z_Prof_2 => NULL() 
   ACHA_RTM_NWP%Z_Prof_3 => NULL() 

   if (Input%Chan_On_67um == symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_67um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_67um =>  NULL()
     ACHA_RTM_NWP%Black_Body_Rad_Prof_67um => NULL()
   endif
   if (Input%Chan_On_85um == symbol%YES) then
     ACHA_RTM_NWP%Atm_Rad_Prof_85um =>  NULL()
     ACHA_RTM_NWP%Atm_Trans_Prof_85um =>  NULL()
   endif
     
   if (Input%Chan_On_11um == symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_11um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_11um => NULL()
      ACHA_RTM_NWP%Black_Body_Rad_Prof_11um => NULL()
   endif
   if (Input%Chan_On_12um == symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_12um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_12um => NULL()
   endif
   if (Input%Chan_On_133um == symbol%YES) then
      ACHA_RTM_NWP%Atm_Rad_Prof_133um => NULL()
      ACHA_RTM_NWP%Atm_Trans_Prof_133um => NULL()
   endif
 

end subroutine NULL_PIX_POINTERS
!====================================================================
!  record svn version as a global variable for output to hdf
!====================================================================
subroutine SET_ACHA_VERSION(Acha_Version)
   character(len=*):: Acha_Version
   Acha_Version = "$Id$"
end subroutine SET_ACHA_VERSION
!====================================================================
! 
! Make a background field of cirrus temperature from appropriate
! retrievals and use as an apriori constraint
!
! Input
!   Cld_Type = standard cloud type values
!   Temperature_Cloud = Cloud-top Temperature 
!   Emissivity_Cloud = Cloud Emissvity
!   Emissivity_Thresh = threshold for determing source pixels
!   Count_Thresh = number of source pixels needed to make a target value
!   Box_Width = pixel dimension of averaging box
!   Missing = Missing value to be used
!
! Output
!   Temperature_Cirrus = cloud temperature of target pixels
!
! Local
!   Mask1 = mask of source pixels
!   Mask2 = mask of target pixels
!====================================================================
subroutine COMPUTE_TEMPERATURE_CIRRUS(Cld_Type, &
                                      Temperature_Cloud,&
                                      Emissivity_Cloud,&
                                      Emissivity_Thresh,&
                                      Count_Thresh, &
                                      Box_Width, &
                                      Missing, &
                                      Temperature_Cirrus)

   integer(kind=int1), intent(in), dimension(:,:):: Cld_Type
   real(kind=real4), intent(in), dimension(:,:):: Temperature_Cloud
   real(kind=real4), intent(in), dimension(:,:):: Emissivity_Cloud
   real(kind=int4), intent(in):: Emissivity_Thresh
   integer(kind=int4), intent(in):: Count_Thresh
   integer(kind=int4), intent(in):: Box_Width
   real(kind=int4), intent(in):: Missing
   real(kind=real4), intent(out), dimension(:,:):: Temperature_Cirrus
   integer(kind=int1), dimension(:,:), allocatable:: Mask1
   integer(kind=int1), dimension(:,:), allocatable:: Mask2

   integer:: Num_Elements
   integer:: Num_Lines

   Temperature_Cirrus = Missing

   Num_Elements = size(Temperature_Cirrus,1)
   Num_Lines = size(Temperature_Cirrus,2)

   allocate(Mask1(Num_Elements,Num_Lines))
   allocate(Mask2(Num_Elements,Num_Lines))

   !---- make source mask
   Mask1 = 0
   where( (Cld_Type == symbol%CIRRUS_TYPE .or. &
           Cld_Type == symbol%OPAQUE_ICE_TYPE .or.  &
           Cld_Type == symbol%OVERLAP_TYPE) .and.  &
           Temperature_Cloud /= Missing .and. &
           Emissivity_Cloud >= Emissivity_Thresh)
      Mask1 = 1
   end where

   !---- make target mask
   Mask2 = 0
   where( (Cld_Type == symbol%CIRRUS_TYPE .or. &
           Cld_Type == symbol%OPAQUE_ICE_TYPE) .and.  &
           Temperature_Cloud /= Missing .and. &
           Emissivity_Cloud < Emissivity_Thresh)
      Mask2 = 1
   end where


   call MEAN_SMOOTH(Mask1,Mask2,Missing,2,2,Count_Thresh,Box_Width,Num_Elements,Num_Lines, &
                    Temperature_Cloud,Temperature_Cirrus)

   !--------------------------------------
   deallocate(Mask1)
   deallocate(Mask2)

end subroutine COMPUTE_TEMPERATURE_CIRRUS
!-------------------------------------------------------------------------------
! Routine to spatially interpret water cloud temperature values to surrounding 
! pixels
!
! input:  interp_flag = 0 (do no interp, assume Zc=2km) /= 0 (do spatial interp)
!-------------------------------------------------------------------------------
subroutine COMPUTE_LOWER_CLOUD_TEMPERATURE(Cld_Type, &
                                        Interp_Flag, &
                                        Surface_Temperature, &
                                        Cloud_Temperature,&
                                        Count_Thresh, &
                                        Box_Width, &
                                        Missing, &
                                        Lower_Tc)

   integer(kind=int1), intent(in), dimension(:,:):: Cld_Type
   integer(kind=int4), intent(in):: Interp_Flag
   real, intent(in), dimension(:,:):: Surface_Temperature
   real(kind=real4), intent(in), dimension(:,:):: Cloud_Temperature
   integer(kind=int4), intent(in):: Count_Thresh
   integer(kind=int4), intent(in):: Box_Width
   real(kind=int4), intent(in):: Missing
   real(kind=real4), intent(out), dimension(:,:):: Lower_Tc

   integer(kind=int1), dimension(:,:), allocatable:: Mask1
   integer(kind=int1), dimension(:,:), allocatable:: Mask2
   integer:: Num_Elements
   integer:: Num_Lines

   !--- initialize output to missing
   Lower_Tc = Missing

   !--- grab size of these arrays
   Num_Elements = size(Cloud_Temperature,1)
   Num_Lines = size(Cloud_Temperature,2)

   !---- make output mask
   allocate(Mask2(Num_Elements,Num_Lines))
   Mask2 = 0
   where(Cld_Type == symbol%OVERLAP_TYPE)
          Mask2 = 1
   end where

   !--- set default to a static offset of surface pressure 
   where(Mask2 == 1 .and. Surface_Temperature /= Missing)
          Lower_Tc = Surface_Temperature - TC_LOWER_CLOUD_OFFSET
   end where

   !--- if no spatial interpolation is to be done, return
   if (Interp_Flag == symbol%NO) then
      deallocate(Mask2)
      return
   endif

   !---- make source mask
   allocate(Mask1(Num_Elements,Num_Lines))
   Mask1 = 0
   where( (Cld_Type == symbol%FOG_TYPE .or. &
              Cld_Type == symbol%WATER_TYPE .or.  &
              Cld_Type == symbol%SUPERCOOLED_TYPE) .and.  &
              Cloud_Temperature /= Missing)
             Mask1 = 1
   end where

   !--- call the spatial analysis routine
   call MEAN_SMOOTH(Mask1,Mask2,Missing,2,2,Count_Thresh,Box_Width,Num_Elements,Num_Lines, &
                    Cloud_Temperature,Lower_Tc)

   !--- deallocate memory
   deallocate(Mask1)
   deallocate(Mask2)

end subroutine COMPUTE_LOWER_CLOUD_TEMPERATURE

!--------------------------------------------------------------------------
! Determine processing order of pixels
!
! Processing Order Description
!
! pass 0 = Not Processed
! pass 1 = non-multi-layer ncc pixels
! pass 2 = single layer water cloud pixels
! pass 3 = lrc multi-layer clouds
! pass 4 = all remaining clouds
! pass 5 = if USE_CIRRUS_FLAG is set on, redo all thin cirrus using a priori
!          temperature from thicker cirrus.
!--------------------------------------------------------------------------
subroutine COMPUTE_PROCESSING_ORDER(Invalid_Data_Mask, &
                                    Cloud_Type, &
                                    ELem_Idx_LRC, Line_Idx_LRC, &  
                                    Pass_Idx_Min,Pass_Idx_Max, &
                                    USE_CIRRUS_FLAG, &
                                    Processing_Order) 
  
  integer(kind=int1), intent(in), dimension(:,:):: Invalid_Data_Mask
  integer(kind=int1), intent(in), dimension(:,:):: Cloud_Type
  integer(kind=int4), intent(in), dimension(:,:):: Elem_Idx_LRC, Line_Idx_LRC
  integer, intent(out):: Pass_Idx_Min, Pass_Idx_Max
  integer(kind=int4), intent(in):: USE_CIRRUS_FLAG
  integer(kind=int1), intent(out), dimension(:,:):: Processing_Order
  integer:: Number_of_Lines, Number_of_Elements
  integer:: Line_Idx, Elem_Idx
  integer:: ilrc, jlrc

  Number_of_Elements = size(Elem_Idx_LRC,1)
  Number_of_Lines = size(Elem_Idx_LRC,2)

  Processing_Order = MISSING_VALUE_INTEGER1
  where(Invalid_Data_Mask == symbol%NO) 
    Processing_Order =  0
  endwhere

  Pass_Idx_min = 1
  Pass_Idx_Max = 4

  if (USE_CIRRUS_FLAG == symbol%YES) Pass_Idx_Max = Pass_Idx_Max + 1

  !--- loop through pixels, determine processing order
  Line_Loop: do Line_Idx = 1, Number_of_Lines
     Element_Loop:   do Elem_Idx = 1, Number_of_Elements

        !--- skip data marked as bad
        if (Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) then
            cycle
        endif

        !--- skip data marked as bad
        if (Cloud_Type(Elem_Idx,Line_Idx) == symbol%CLEAR_TYPE .or. & 
            Cloud_Type(Elem_Idx,Line_Idx) == symbol%PROB_CLEAR_TYPE) then
            Processing_Order(Elem_Idx,Line_Idx) = 0
            cycle
        endif

        ilrc = Elem_Idx_LRC(Elem_Idx,Line_Idx)
        jlrc = Line_Idx_LRC(Elem_Idx,Line_Idx)

        !-- on pass 1, do single layer lrc's
        if ((Elem_Idx == ilrc) .and. (Line_Idx == jlrc) .and. &
            (Cloud_Type(Elem_Idx,Line_Idx) /= symbol%OVERLAP_TYPE)) then
             Processing_Order(Elem_Idx,Line_Idx) = 1
             cycle
        endif

        !-- on pass 2, do non-lrc water clouds
        if (((Elem_Idx /= ilrc) .or. (Line_Idx /= jlrc)) .and. &
            (Cloud_Type(Elem_Idx,Line_Idx) == symbol%FOG_TYPE .or. &
            Cloud_Type(Elem_Idx,Line_Idx) == symbol%WATER_TYPE .or. &
            Cloud_Type(Elem_Idx,Line_Idx) == symbol%MIXED_TYPE .or. &
            Cloud_Type(Elem_Idx,Line_Idx) == symbol%SUPERCOOLED_TYPE)) then
            Processing_Order(Elem_Idx,Line_Idx) = 2
            cycle
        endif

        !-- on pass 3, do lrc overlap clouds
        if ((Elem_Idx == ilrc) .and. (Line_Idx == jlrc) .and. &
            (Cloud_Type(Elem_Idx,Line_Idx) == symbol%OVERLAP_TYPE)) then
             Processing_Order(Elem_Idx,Line_Idx) = 3
            cycle
        endif

        !--  on pass-4 do remaining
        if (Processing_Order(Elem_Idx,Line_Idx) == 0) then
           Processing_Order(Elem_Idx,Line_Idx) = 4
        endif

       end do Element_Loop
     end do Line_Loop

end subroutine COMPUTE_PROCESSING_ORDER
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
  !-------------------------------------------------------------------------------------------------
  ! Smooth a field using a mean over an area 
  !
  ! Description
  !    Values of Z_in with Mask_In = 1 are used to populate pixels with Mask_Out = 1 
  !    Z_Out is computed as the mean of Z_In*Mask_In over a box whose size is 
  !    defined by N.  
  !
  ! Input
  !    Mask_In - binary mask of point used as the source of the smoothing
  !    Mask_Out - binary mask of points to have a result upon exit
  !    Missin = missing value used as fill for Z_Out
  !    Count_Thresh - number of source points to compue an output
  !    N - half-width of smoothing box (x and y)
  !    Num_Elements = size of array in x-direction
  !    Num_Lines = size of array in y-direction
  !    Z_In - source values
  !    Z_Out - output values
  !    di = number of pixels to skip in the i direction (0=none,1=every other ...)
  !    dj = number of pixels to skip in the j direction (0=none,1=every other ...)
  !
  !-------------------------------------------------------------------------------------------------
  subroutine MEAN_SMOOTH(Mask_In,Mask_Out,Missing,di,dj,Count_Thresh,N,Num_Elements, Num_Lines, Z_In,Z_Out)

     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_In
     integer (kind=int1), intent(in), dimension(:,:):: Mask_Out
     real (kind=real4), intent(in):: Missing
     integer (kind=int4), intent(in):: di
     integer (kind=int4), intent(in):: dj
     integer (kind=int4), intent(in):: Count_Thresh
     integer (kind=int4), intent(in):: N
     integer (kind=int4), intent(in):: Num_Elements
     integer (kind=int4), intent(in):: Num_lines
     real (kind=real4), intent(in), dimension(:,:), target:: Z_In
     real (kind=real4), intent(out), dimension(:,:):: Z_Out
     integer:: i
     integer:: j
     real:: Count_Temporary
     real, pointer, dimension(:,:):: Z_In_Sub
     integer (kind=int1), pointer, dimension(:,:):: Mask_In_Sub
     real (kind=real8):: Z_sum
     integer:: i1,i2,j1,j2

     do j = 1 + dj, Num_Lines-dj, dj + 1

        j1 = min(Num_Lines,max(1,j - N))
        j2 = min(Num_Lines,max(1,j + N))

        do i = 1 + N + di, Num_Elements - N - di, di + 1

          i1 = i - N
          i2 = i + N

          Z_Out(i,j) = Missing

          if (maxval(Mask_Out(i-di:i+di,j-dj:j+dj)) == 0) cycle
          !if (Mask_Out(i,j) == 0) cycle

          Mask_In_Sub => Mask_In(i1:i2,j1:j2)
          Z_In_Sub => Z_In(i1:i2,j1:j2)

          Count_Temporary = sum(real(Mask_In_Sub))
          if (Count_Temporary < Count_Thresh) cycle

          Z_Sum = sum(Z_In_Sub*Mask_In_Sub) / Count_Temporary

          Z_Out(i-di:i+di,j-dj:j+dj) = real(Z_Sum)

          Z_In_Sub => null() 
          Mask_In_Sub => null() 

      enddo
     enddo   

    !--- only values are missing where output mask is 0
    where(Mask_Out == 0)
      Z_Out = Missing
    endwhere

  end subroutine MEAN_SMOOTH


!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------
function EMPIRICAL_LAPSE_RATE(Tsfc, Tc, land_flag) result(lapse_rate)
  real, intent(in):: Tsfc
  real, intent(in):: Tc
  integer, intent(in):: land_flag  !(0=ocean,1=land)
  real:: Tcs
  real:: lapse_rate
  integer:: its, itcs

  Tcs = Tc - Tsfc
  its = int((Tsfc - ts_min) / dts) + 1
  its = max(1,min(nts,its))
  itcs = int((Tcs - tcs_min) / dtcs) + 1
  itcs = max(1,min(ntcs,itcs))

  if (land_flag == 0) then 
    lapse_rate = ocean_lapse_rate_table(its,itcs)
  else
    lapse_rate = land_lapse_rate_table(its,itcs)
  endif

end function EMPIRICAL_LAPSE_RATE

function OCEANIC_LAPSE_RATE_OLD(Tsfc, Tc) result(lapse_rate)
  real, intent(in):: Tsfc
  real, intent(in):: Tc
  real:: Tcs
  real:: lapse_rate

  Tcs = Tc - Tsfc

  lapse_rate =  -0.061 - 1.67*Tcs - 0.124*(Tcs**2) - 0.00343*(Tcs**3)
   

end function OCEANIC_LAPSE_RATE_OLD
!----------------------------------------------------------------------------
! estimate cirrus aprior temperature and uncertainty from a precomputed 
! latitude table (stored in acha_parameters.inc)
!----------------------------------------------------------------------------
subroutine compute_cirrus_apriori(t_tropo, latitude, tc_apriori, tc_apriori_uncer)
  real, intent(in):: t_tropo
  real, intent(in):: latitude
  real, intent(out):: tc_apriori
  real, intent(out):: tc_apriori_uncer

  integer:: lat_idx
  real, parameter:: lat_min = -90.0
  real, parameter:: delta_lat = -10.0

  lat_idx = int((latitude - lat_min) / delta_lat) + 1
  lat_idx = max(1,min(lat_idx, num_lat_cirrus_ap))
  
  tc_apriori = t_tropo + TC_CIRRUS_MEAN_LAT_VECTOR(lat_idx)
  tc_apriori_uncer = TC_CIRRUS_STDDEV_LAT_VECTOR(lat_idx)

end subroutine compute_cirrus_apriori

  FUNCTION get_lun_acha() RESULT( lun )


    ! -----------------
    ! Type declarations
    ! -----------------

    INTEGER :: lun
    LOGICAL :: file_open


    ! --------------------------------------------
    ! Initialise logical unit number and file_open
    ! --------------------------------------------

    lun = 9
    file_open = .TRUE.


    ! ------------------------------
    ! Start open loop for lun search
    ! ------------------------------

    lun_search: DO

      ! -- Increment logical unit number
      lun = lun + 1

      ! -- Check if file is open
      INQUIRE( lun, OPENED = file_open )

      ! -- Is this lun available?
      IF ( .NOT. file_open ) EXIT lun_search

    END DO lun_search

  END FUNCTION get_lun_acha

end module AWG_CLOUD_HEIGHT
