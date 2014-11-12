!$Id:$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: naive_bayesian_clavrx_bridge_module.f90 (src)
!       naive_bayesian_clavrx_bridge_module (program)
!
! PURPOSE: 
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
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
module NB_CLOUD_MASK_CLAVRX_BRIDGE

   ! -- MODULES USED
   use PIXEL_COMMON
   use NUMERICAL_ROUTINES
   use NB_CLOUD_MASK
   use CLOUD_MASK_ADDONS
   use NB_CLOUD_MASK_SERVICES

   implicit none

   public :: NB_CLOUD_MASK_BRIDGE

   private :: COVARIANCE_LOCAL
   private :: SET_SYMBOL
   private :: SET_INPUT
   private :: SET_OUTPUT
   private :: SET_DIAG
   private :: NULL_INPUT
   private :: NULL_OUTPUT
   private :: NULL_DIAG

   !--- define these structure as module wide
   type(mask_input), private :: Input   
   type(mask_output), private :: Output   
   type(diag_output), private :: Diag  
   type(symbol_naive_bayesian),private :: Symbol
   
contains
!----------------------------------------------------------------------
! Bridge Routine
!
! Note, the Diag argument is optional
!---------------------------------------------------------------------- 
 subroutine NB_CLOUD_MASK_BRIDGE(Segment_Number)
 
   implicit none

   integer, intent(in):: Segment_Number

   integer :: i, j
   character (len=355):: Ancil_Data_Path
   character (len=255):: Naive_Bayes_File_Name

   !---- set paths and mask classifier file name to their values in this framework
   Ancil_Data_Path = Ancil_Data_Dir
   Naive_Bayes_File_Name = Bayesian_Cloud_Mask_Name

   !--- set structure (symbol, input, output, diag)  elements to corresponding values in this framework
   call SET_SYMBOL()

   ! -----------    loop over pixels -----   
   line_loop: do i = 1, Num_Pix
      elem_loop: do  j = 1, Num_Scans_Read
         call SET_INPUT(i,j)
         call SET_OUTPUT(i,j)
         call SET_DIAG(i,j)

         !---call cloud mask routine
         call NB_CLOUD_MASK_ALGORITHM(Ancil_Data_Path,  &
                      Naive_Bayes_File_Name, &
                      Symbol,  &
                      Input, &
                      Output)
                     ! Diag)   !optional



         !--- call non-cloud detection routines (smoke, dust and fire)
         call NB_CLOUD_MASK_ADDONS_ALGORITHM(Symbol,  &
                                          Input, &
                                          Output) !, &
                                         ! Diag)   !optional

         !--- nullify pointers within these data structures
         call NULL_INPUT()
         call NULL_OUTPUT()
         call NULL_DIAG()
   
         !-----------------------------------------------------------------------
         ! CLAVR-x specific processing
         !-----------------------------------------------------------------------

         !--- unpack elements of the cloud test vector into clavr-x global arrays
         Bayes_Mask_Sfc_Type_Global(i,j) = ibits(Cld_Test_Vector_Packed(3,i,j),0,3)

      end do elem_loop
   end do line_loop

   !--- grab version tags for output as attributes in level2
   !--- only need to do this once, so do on first segment
   if (Segment_Number == 1) then
     call SET_CLOUD_MASK_VERSION(Cloud_Mask_Version)
     call SET_CLOUD_MASK_THRESHOLDS_VERSION(Cloud_Mask_Thresholds_Version)
   endif

   end subroutine NB_CLOUD_MASK_BRIDGE

   !====================================================================
   ! Function Name: Covariance_LOCAL
   !
   ! Function:
   !    Compute the Covariance for two mxn arrays
   !
   ! Description: Covariance = E(XY) - E(X)*E(Y)
   !   
   ! Calling Sequence: BT_WV_BT_Window_Covar(Elem_Idx,Line_Idx) = Covariance( &
   !                       sat%bt10(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
   !                       sat%bt14(Arr_Right:Arr_Left,Arr_Top:Arr_Bottom), &
   !                      Array_Width, Array_Hgt)
   !   
   !
   ! Inputs:
   !   Array 1 - the first array (X)
   !   Array 2 - the second array (Y)
   !   Elem_size
   !   Line_size
   !
   ! Outputs: 
   !   Covariance of X and Y
   !
   ! Dependencies:
   !        none
   !
   ! Restrictions:  None
   !
   ! Reference: Standard definition for the Covariance Computation
   !
   !====================================================================
   function COVARIANCE_LOCAL &
        (Array_One,Array_Two,Array_Width,Array_Hght,Invalid_Data_Mask) &
         RESULT(Covar_Array_One_Array_Two)

   real(kind=real4), intent(in), dimension(:,:):: Array_One
   real(kind=real4), intent(in), dimension(:,:):: Array_Two
   integer(kind=INT4), intent(in):: Array_Width
   integer(kind=INT4), intent(in):: Array_Hght
   integer(kind=INT1), intent(in), dimension(:,:):: Invalid_Data_Mask

   real(kind=real8):: Mean_Array_One
   real(kind=real8):: Mean_Array_Two
   real(kind=real8):: Mean_Array_One_x_Array_Two
   real(kind=real8):: Sum_Array_One
   real(kind=real8):: Sum_Array_Two
   real(kind=real8):: Sum_Array_One_x_Array_Two
   real(kind=real4):: Covar_Array_One_Array_Two

   !--- skip computation for pixel arrays with any missing data
   if (sum(Invalid_Data_Mask) > 0) then
      Covar_Array_One_Array_Two = Missing_Value_Real4
      return
   endif

   Sum_Array_One = sum(Array_One)
   Sum_Array_Two = sum(Array_Two)

   Mean_Array_One = Sum_Array_One / (Array_Width*Array_Hght)
   Mean_Array_Two = Sum_Array_Two / (Array_Width*Array_Hght)

   Sum_Array_One_x_Array_Two = sum(Array_One*Array_Two)
   Mean_Array_One_x_Array_Two = Sum_Array_One_x_Array_Two / (Array_Width*Array_Hght)
   
   Covar_Array_One_Array_Two  = Mean_Array_One_x_Array_Two - &
                                Mean_Array_One * Mean_Array_Two 
   
   end function COVARIANCE_LOCAL

   !============================================================================
   ! set symbols
   !============================================================================
   subroutine SET_SYMBOL()

      symbol%CLOUDY = sym%CLOUDY
      symbol%PROB_CLOUDY = sym%PROB_CLOUDY
      symbol%PROB_CLEAR = sym%PROB_CLEAR
      symbol%CLEAR = sym%CLEAR

      symbol%NO = sym%NO
      symbol%YES = sym%YES

      symbol%WATER_SFC = sym%WATER_SFC
      symbol%EVERGREEN_NEEDLE_SFC = sym%EVERGREEN_NEEDLE_SFC
      symbol%EVERGREEN_BROAD_SFC = sym%EVERGREEN_BROAD_SFC
      symbol%DECIDUOUS_NEEDLE_SFC = sym%DECIDUOUS_NEEDLE_SFC
      symbol%DECIDUOUS_BROAD_SFC = sym%DECIDUOUS_BROAD_SFC
      symbol%MIXED_FORESTS_SFC = sym%MIXED_FORESTS_SFC
      symbol%WOODLANDS_SFC = sym%WOODLANDS_SFC
      symbol%WOODED_GRASS_SFC = sym%WOODED_GRASS_SFC
      symbol%CLOSED_SHRUBS_SFC = sym%CLOSED_SHRUBS_SFC
      symbol%OPEN_SHRUBS_SFC = sym%OPEN_SHRUBS_SFC
      symbol%GRASSES_SFC = sym%GRASSES_SFC
      symbol%CROPLANDS_SFC = sym%CROPLANDS_SFC
      symbol%BARE_SFC = sym%BARE_SFC
      symbol%URBAN_SFC = sym%URBAN_SFC

      symbol%SHALLOW_OCEAN = sym%SHALLOW_OCEAN
      symbol%LAND = sym%LAND
      symbol%COASTLINE = sym%COASTLINE
      symbol%SHALLOW_INLAND_WATER = sym%SHALLOW_INLAND_WATER
      symbol%EPHEMERAL_WATER = sym%EPHEMERAL_WATER
      symbol%DEEP_INLAND_WATER = sym%DEEP_INLAND_WATER
      symbol%MODERATE_OCEAN = sym%MODERATE_OCEAN
      symbol%DEEP_OCEAN = sym%DEEP_OCEAN

      symbol%NO_SNOW = sym%NO_SNOW
      symbol%SEA_ICE = sym%SEA_ICE
      symbol%SNOW = sym%SNOW   
   end subroutine SET_SYMBOL

   !============================================================================
   ! set input pointers
   !============================================================================
   subroutine SET_INPUT(i,j)
      integer, intent (in) :: i, j

      Input%Num_Elem = Num_Pix
      Input%Num_Line = Num_Scans_Read
      Input%Num_Line_Max = Num_Scans_Per_Segment
      Input%Num_Segments = Num_Segments
      Input%Invalid_Data_Mask => Bad_Pixel_Mask(i,j)
      Input%Chan_On_041um = Chan_On_Flag_Default(8)
      Input%Chan_On_063um = Chan_On_Flag_Default(1)
      Input%Chan_On_086um = Chan_On_Flag_Default(2)
      Input%Chan_On_138um = Chan_On_Flag_Default(26)
      Input%Chan_On_160um = Chan_On_Flag_Default(6)
      Input%Chan_On_213um = Chan_On_Flag_Default(7)
      Input%Chan_On_375um = Chan_On_Flag_Default(20)
      Input%Chan_On_67um = Chan_On_Flag_Default(27)
      Input%Chan_On_85um = Chan_On_Flag_Default(29)
      Input%Chan_On_11um = Chan_On_Flag_Default(31)
      Input%Chan_On_12um = Chan_On_Flag_Default(32)
      Input%Chan_On_I1_064um = Chan_On_Flag_Default(37)
      Input%Chan_On_I4_374um = Chan_On_Flag_Default(40)
      Input%Chan_On_I5_114um = Chan_On_Flag_Default(41)
      Input%Chan_On_DNB = Chan_On_Flag_Default(42)
      Input%Snow_Class => Snow(i,j)
      Input%Land_Class => Land(i,j)
      Input%Oceanic_Glint_Mask => Glint_Mask(i,j)
      Input%Lunar_Oceanic_Glint_Mask => Glint_Mask_Lunar(i,j)
      Input%Coastal_Mask => Coast_Mask(i,j)
      Input%Solzen => Solzen(i,j)
      Input%Scatzen => Scatangle(i,j)
      Input%Lunscatzen => Scatangle_Lunar(i,j)
      Input%Senzen => Satzen(i,j)
      Input%Lunzen => Lunzen(i,j)
      Input%Lat => Lat(i,j)
      Input%Lon => Lon(i,j)
      Input%Ref_041um => ch(8)%Ref_Toa(i,j)
      Input%Ref_063um => ch(1)%Ref_Toa(i,j)
      Input%Ref_063um_Clear => ch(1)%Ref_Toa_Clear(i,j)
      Input%Ref_063um_Std => Ref_Ch1_Std_3x3(i,j)
      Input%Ref_063um_Min => Ref_Ch1_Min_3x3(i,j)
      Input%Ref_086um => ch(2)%Ref_Toa(i,j)
      Input%Ref_138um => ch(26)%Ref_Toa(i,j)
      Input%Ref_160um => ch(6)%Ref_Toa(i,j)
      Input%Ref_160um_Clear => ch(6)%Ref_Toa_Clear(i,j)
      Input%Ref_375um => ch(20)%Ref_Toa(i,j)
      Input%Ref_375um_Clear => ch(20)%Ref_Toa_Clear(i,j)
      Input%Ref_213um => ch(7)%Ref_Toa(i,j)
      Input%Bt_375um => ch(20)%Bt_Toa(i,j)
      Input%Bt_375um_Std => Bt_Ch20_Std_3x3(i,j)
      Input%Emiss_375um =>  Ems_Ch20_Median_3x3(i,j)
      Input%Emiss_375um_Clear => Ems_Ch20_Clear_Solar_Rtm(i,j)
      Input%Bt_67um => ch(27)%Bt_Toa (i,j)
      Input%Bt_85um => ch(29)%Bt_Toa(i,j)
      Input%Bt_11um => ch(31)%Bt_Toa(i,j)
      Input%Bt_11um_Std => Bt_Ch31_Std_3x3(i,j)
      Input%Bt_11um_Max => Bt_Ch31_Max_3x3(i,j)
      Input%Bt_11um_Clear => ch(31)%Bt_Toa_Clear(i,j)
      Input%Emiss_11um_Tropo => ch(31)%Emiss_Tropo(i,j)
      Input%Bt_12um => ch(32)%Bt_Toa(i,j)
      Input%Bt_12um_Clear => ch(32)%Bt_Toa_Clear(i,j)
      Input%Ref_I1_064um_Std => Ref_Uni_ChI1(i,j)
      Input%Bt_I4_374um_Std => Bt_Uni_ChI4(i,j)
      Input%Bt_I5_114um_Std => Bt_Uni_ChI5(i,j)
      Input%Bt_11um_Bt_67um_Covar => Covar_Ch27_Ch31_5x5(i,j)
      Input%Sst_Anal_Uni => Sst_Anal_Uni(i,j)
      Input%Emiss_Sfc_375um => ch(20)%Sfc_Emiss(i,j)
      Input%Rad_Lunar => ch(42)%Rad_Toa(i,j)
      Input%Ref_Lunar => ch(42)%Ref_Lunar_Toa(i,j)
      Input%Ref_Lunar_Min => Ref_ChDNB_Lunar_Min_3x3(i,j)
      Input%Ref_Lunar_Std => Ref_ChDNB_Lunar_Std_3x3(i,j)
      Input%Ref_Lunar_Clear => ch(42)%Ref_Lunar_Toa_Clear(i,j)
      Input%Zsfc => Zsfc(i,j)
      Input%Solar_Contamination_Mask => Solar_Contamination_Mask(i,j)
      Input%Sfc_Type => Sfc_Type(i,j)
   end subroutine SET_INPUT

   subroutine SET_OUTPUT(i,j)
      integer, intent (in) :: i, j

      Output%Cld_Flags_Packed => Cld_Test_Vector_Packed(:,i,j)
      Output%Cld_Mask_Bayes => Cld_Mask(i,j)
      Output%Posterior_Cld_Probability => Posterior_Cld_Probability(i,j)
      Output%Dust_Mask => Dust_Mask(i,j)
      Output%Smoke_Mask => Smoke_Mask(i,j)
      Output%Fire_Mask => Fire_Mask(i,j)
   end subroutine SET_OUTPUT

   subroutine SET_DIAG(i,j)
      integer, intent (in) :: i, j

      Diag%Array_1 => Diag_Pix_Array_1(i,j)
      Diag%Array_2 => Diag_Pix_Array_2(i,j)
      Diag%Array_3 => Diag_Pix_Array_3(i,j)
   end subroutine SET_DIAG

   !============================================================================
   ! nullify input pointers
   !============================================================================
   subroutine NULL_INPUT()
      Input%Invalid_Data_Mask => null()
      Input%Snow_Class => null()
      Input%Land_Class => null() 
      Input%Oceanic_Glint_Mask => null() 
      Input%Lunar_Oceanic_Glint_Mask => null()
      Input%Coastal_Mask => null()
      Input%Solzen => null()
      Input%Scatzen => null()
      Input%Lunscatzen => null()
      Input%Senzen => null()
      Input%Lunzen => null()
      Input%Lat => null()
      Input%Lon => null()
      Input%Ref_041um => null()
      Input%Ref_063um => null()
      Input%Ref_063um_Clear => null()
      Input%Ref_063um_Std => null()
      Input%Ref_063um_Min => null()
      Input%Ref_086um => null()
      Input%Ref_138um => null()
      Input%Ref_160um => null()
      Input%Ref_160um_Clear => null()
      Input%Ref_213um => null()
      Input%Bt_375um =>  null()
      Input%Bt_375um_Std => null()
      Input%Emiss_375um => null()
      Input%Emiss_375um_Clear => null()
      Input%Bt_67um => null()
      Input%Bt_85um => null()
      Input%Bt_11um => null()
      Input%Bt_11um_Std => null()
      Input%Bt_11um_Max => null()
      Input%Bt_11um_Clear => null()
      Input%Emiss_11um_Tropo => null()
      Input%Bt_12um => null()
      Input%Bt_12um_Clear => null()
      Input%Ref_I1_064um_Std => null()
      Input%Bt_I4_374um_Std => null()
      Input%Bt_I5_114um_Std => null()
      Input%Bt_11um_Bt_67um_Covar => null() 
      Input%Sst_Anal_Uni => null()
      Input%Emiss_Sfc_375um => null()
      Input%Rad_Lunar => null()
      Input%Ref_Lunar => null()
      Input%Ref_Lunar_Min => null()
      Input%Ref_Lunar_Std => null()
      Input%Ref_Lunar_Clear => null()
      Input%Zsfc => null()
      Input%Solar_Contamination_Mask => null()
      Input%Sfc_Type => null()
   end subroutine NULL_INPUT

   !============================================================================
   ! nullify output pointers
   !============================================================================
   subroutine NULL_OUTPUT()
    Output%Cld_Flags_Packed => null()
    Output%Cld_Mask_Bayes => null()
    Output%Posterior_Cld_Probability => null()
    Output%Dust_Mask => null()
    Output%Smoke_Mask => null()
    Output%Fire_Mask => null()
   end subroutine NULL_OUTPUT

   !============================================================================
   ! nullify diag pointers
   !============================================================================
   subroutine NULL_DIAG()
      Diag%Array_1 => null()
      Diag%Array_2 => null()
      Diag%Array_3 => null()
   end subroutine NULL_DIAG

   !============================================================================

end module NB_CLOUD_MASK_CLAVRX_BRIDGE
