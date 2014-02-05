!$Id: naive_bayesian_clavrx_bridge_module.f90 9 2014-01-31 08:19:35Z awalther $
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
module naive_bayesian_clavrx_bridge_module


   ! -- MODULES USED

   use CONSTANTS
   use PIXEL_COMMON
   use NUMERICAL_ROUTINES
   use NAIVE_BAYESIAN_CLOUD_MASK


   implicit none

   public :: AWG_CLOUD_BAYES_BRIDGE
   private :: covariance_local
   
   
contains
!----------------------------------------------------------------------
!
!---------------------------------------------------------------------- 
 subroutine AWG_CLOUD_BAYES_BRIDGE(Segment_Number)
 
   implicit none

   integer, intent(in):: Segment_Number
   character (len=120):: Ancil_Data_Path
   character (len=120):: Naive_Bayes_File_Name

   integer:: Num_Elem
   integer:: Num_Line
   integer:: Num_Line_Max
   type(symbol_naive_bayesian) :: symbol

   integer:: Chan_On_063um
   integer:: Chan_On_086um
   integer:: Chan_On_138um
   integer:: Chan_On_160um
   integer:: Chan_On_375um
   integer:: Chan_On_67um
   integer:: Chan_On_85um
   integer:: Chan_On_11um
   integer:: Chan_On_12um
   integer:: Chan_On_DNB
   
   
   !Initialize local pointers to global variables

   Ancil_Data_Path = Ancil_Data_Dir
   Naive_Bayes_File_Name = Bayesian_Cloud_Mask_Name
   Num_Elem = Num_Pix
   Num_Line = Num_Scans_Read
   Num_Line_Max = Num_Scans_Per_Segment
   
   
   !----set symbols to local values
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

   !------------------------------------------------------------------
   ! store channel mappings into flags sent through bridge
   ! clavrx uses the MODIS channel mapping
   !------------------------------------------------------------------
   Chan_On_063um = Chan_On_Flag_Default(1)
   Chan_On_086um = Chan_On_Flag_Default(2)
   Chan_On_138um = Chan_On_Flag_Default(26)
   Chan_On_160um = Chan_On_Flag_Default(6)
   Chan_On_375um = Chan_On_Flag_Default(20)
   Chan_On_67um = Chan_On_Flag_Default(27)
   Chan_On_85um = Chan_On_Flag_Default(29)
   Chan_On_11um = Chan_On_Flag_Default(31)
   Chan_On_12um = Chan_On_Flag_Default(32)
   Chan_On_DNB = Chan_On_Flag_Default(42)
   
   !Call Naive bayesian routine
   
   call CLOUD_MASK_NAIVE_BAYES(Ancil_Data_Path,  &
                               Naive_Bayes_File_Name, &
                               symbol,  &
                               Num_Elem,  &
                               Num_Line, &
                               Num_Line_Max, &
                               Bad_Pixel_Mask,  &
                               Cld_Test_Vector_Packed, &
                               Chan_On_063um,  &
                               Chan_On_086um,  &
                               Chan_On_138um,  &
                               Chan_On_160um,  &
                               Chan_On_375um,  &
                               Chan_On_67um,  &
                               Chan_On_85um,  &
                               Chan_On_11um,  &
                               Chan_On_12um,  &
                               Chan_On_DNB,  &
                               Snow,  &
                               Land, &
                               Glint_Mask,  &
                               Glint_Mask_Lunar,  &
                               Coast_Mask, &
                               Solzen,  &
                               Scatangle, &
                               Scatangle_Lunar, &
                               Satzen, &
                               Lunzen, &
                               Lat,  &
                               Lon, &
                               ch(1)%Ref_Toa, &
                               ch(1)%Ref_Toa_Clear, &
                               Ref_Ch1_Std_3x3,  &
                               Ref_Ch1_Min_3x3, &
                               ch(2)%Ref_Toa, &
                               ch(26)%Ref_Toa, &
                               ch(6)%Ref_Toa, &
                               ch(6)%Ref_Toa_Clear, &
                               ch(20)%Bt_Toa,  &
                               Bt_Ch20_Std_3x3,  &
                               Ems_Ch20_Median_3x3,  &    !needed?
                               Ems_Ch20_Clear_Solar_Rtm, &
                               ch(27)%Bt_Toa,  &
                               ch(29)%Bt_Toa, &
                               ch(31)%Bt_Toa,  &
                               Bt_Ch31_Std_3x3, &
                               Bt_Ch31_Max_3x3,  &
                               ch(31)%Bt_Toa_Clear,  &
                               Bt_Ch31_LRC,  &
                               ch(31)%Emiss_Tropo, &
                               Emiss_11um_Tropo_LRC,  &
                               ch(32)%Bt_Toa, &
                               ch(32)%Bt_Toa_Clear,  &
                               Covar_Ch27_Ch31_5x5, &
                               Sst_Anal_Uni, &
                               ch(20)%Sfc_Emiss, &
                               ch(42)%Rad_Toa, &
                               ch(42)%Ref_Lunar_Toa, &
                               Ref_ChDNB_Lunar_Min_3x3, &
                               Ref_ChDNB_Lunar_Std_3x3, &
                               ch(42)%Ref_Lunar_Toa_Clear, &
                               Zsfc,  &
                               Num_Segments, &
                               Solar_Contamination_Mask,  &
                               Sfc_Type,  &
                               Cld_Mask, &
                               Cloud_Mask_Bayesian_Flag, &      !REMOVE
                               Posterior_Cld_Probability, &
                               Diag_Pix_Array_1, &
                               Diag_Pix_Array_2, &
                               Diag_Pix_Array_3)

   
   !--- unpack elements of the cloud test vector into clavr-x global arrays
   Bayes_Mask_Sfc_Type_Global = ibits(Cld_Test_Vector_Packed(3,:,:),0,3)

   !--- grab version tags for output as attributes in level2
   !--- only need to do this once, so do on first segment
   if (Segment_Number == 1) then
     call SET_CLOUD_MASK_VERSION(Cloud_Mask_Version)
     call SET_CLOUD_MASK_THRESHOLDS_VERSION(Cloud_Mask_Thresholds_Version)
   endif

   end subroutine AWG_CLOUD_BAYES_BRIDGE

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
function covariance_local &
        (Array_One,Array_Two,Array_Width,Array_Hght,Invalid_Data_Mask) &
         RESULT(Covar_Array_One_Array_Two)

   real(kind=real4), intent(in), dimension(:,:):: Array_One
   real(kind=real4), intent(in), dimension(:,:):: Array_Two
   INTEGER(kind=INT4), intent(in):: Array_Width
   INTEGER(kind=INT4), intent(in):: Array_Hght
   INTEGER(kind=INT1), intent(in), dimension(:,:):: Invalid_Data_Mask

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
   
 end function covariance_local



end module naive_bayesian_clavrx_bridge_module
