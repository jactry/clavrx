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

! temp_11_0um_nom_sounder

   ! -- MODULES USED
   use PIXEL_COMMON, only: &
       Ch, &
       Geo, & 
       Nav, &
       Sfc, &
       Image, &
       Sensor, &
       Ancil_Data_Dir, &
       Bayesian_Cloud_Mask_Name, &
       Bayes_Mask_Sfc_Type_Global, &
       Bad_Pixel_Mask, &
       Cld_Mask, &
       Cld_Test_Vector_Packed, &
       Dust_Mask, &
       Fire_Mask, &
       Posterior_Cld_Probability, &
       Smoke_Mask, &
       SST_Anal_Uni, &
       Solar_Contamination_Mask, &
       Bt_11um_Sounder, &
       AVHRR_IFF_Flag, &
       Ref_Ch1_Std_3x3, &
       Ref_Ch1_Min_3x3, &
       Bt_Ch20_Std_3x3, &
       Ems_Ch20_Median_3x3, &
       Ems_Ch20_Clear_Rtm, &
       Covar_Ch27_Ch31_5x5, &
       Bt_Ch31_Std_3x3, &
       Bt_Ch31_Max_3x3, &
       Ref_ChDNB_Lunar_Std_3x3, &
       Ref_ChDNB_Lunar_Min_3x3, &
       Ref_Uni_ChI1, &
       Bt_Uni_ChI4, &
       Bt_Uni_ChI5, &
       Diag_Pix_Array_1, &
       Diag_Pix_Array_2, &
       Diag_Pix_Array_3

   use CONSTANTS, only: &
       sym, &
       CLOUD_MASK_VERSION, & 
       Cloud_Mask_Thresholds_Version
       
   use NB_CLOUD_MASK
   use CLOUD_MASK_ADDONS
   use NB_CLOUD_MASK_SERVICES
   use CLAVRX_MESSAGE_MODULE, only: MESG

   implicit none

   public :: NB_CLOUD_MASK_BRIDGE

   private :: COVARIANCE_LOCAL
   private :: SET_SYMBOL
   private :: SET_INPUT
   private :: SET_OUTPUT
   private :: SET_DIAG
   private :: MEDIAN_FILTER_MASK

   !--- define these structure as module wide
   type(mask_input), private :: Input   
   type(mask_output), private :: Output   
   type(diag_output), private :: Diag  
   type(symbol_naive_bayesian), private :: Symbol
   
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
   integer :: ELem_Idx_1, ELem_Idx_2, Line_Idx_1, Line_Idx_2
   integer, parameter :: N_box = 2
   character (len=555):: Naive_Bayes_File_Name_Full_Path

   logical, save:: First_Call = .true.
   integer, parameter:: Nmed = 2
   real(kind=real4):: Nmed_Total
   integer(kind=int1), dimension(:,:), allocatable:: I1_Temp_1
   integer(kind=int1), dimension(:,:), allocatable:: I1_Temp_2


   if (First_Call .eqv. .true.) then
       call MESG('NB Cloud Mask starts ', color = 46)
   endif

   !---- set paths and mask classifier file name to their values in this framework
   Naive_Bayes_File_Name_Full_Path = trim(Ancil_Data_Dir)// &
         "static/luts/nb_cloud_mask/"//trim(Bayesian_Cloud_Mask_Name)

   !--- set structure (symbol, input, output, diag)  
   !    elements to corresponding values in this framework
   call SET_SYMBOL()

   ! -----------    loop over pixels -----   
   line_loop: do i = 1, Image%Number_Of_Elements
      elem_loop: do  j = 1, Image%Number_Of_Lines_Read_This_Segment
         call SET_INPUT(i,j)

         !---call cloud mask routine
         call NB_CLOUD_MASK_ALGORITHM(  &
                      Naive_Bayes_File_Name_Full_Path, &
                      Symbol,  &
                      Input, &
                      Output)
                     ! Diag)   !optional

         !--- call non-cloud detection routines (smoke, dust and fire)
         call NB_CLOUD_MASK_ADDONS_ALGORITHM(Symbol,  &
                                          Input, &
                                          Output) !, &
                                          !Diag)   !optional

         call SET_OUTPUT(i,j)
         call SET_DIAG(i,j)
   
         !-----------------------------------------------------------------------
         ! CLAVR-x specific processing
         !-----------------------------------------------------------------------

         !--- unpack elements of the cloud test vector into clavr-x global arrays
         Bayes_Mask_Sfc_Type_Global(i,j) = ibits(Cld_Test_Vector_Packed(3,i,j),0,3)

      end do elem_loop
   end do line_loop
   !------------------------------------------------------------------------------
   ! Apply Median Filters
   !------------------------------------------------------------------------------
   Nmed_Total = (2.0*Nmed+1)**2
   allocate(I1_Temp_1(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
   allocate(I1_Temp_2(Image%Number_Of_Elements,Image%Number_Of_Lines_Per_Segment))
   I1_Temp_1 = Dust_Mask
   I1_Temp_2 = Smoke_Mask
   I1_Temp_1(Nmed:Image%Number_Of_Elements-Nmed,1+Nmed:Image%Number_Of_Lines_Read_This_Segment-Nmed) = 0
   I1_Temp_2(Nmed:Image%Number_Of_Elements-Nmed,1+Nmed:Image%Number_Of_Lines_Read_This_Segment-Nmed) = 0
   line_loop_median: do i = 1+Nmed, Image%Number_Of_Elements-Nmed
      elem_loop_median: do  j = 1+Nmed, Image%Number_Of_Lines_Read_This_Segment-Nmed
         I1_Temp_1(i,j) = nint(count(Dust_Mask(i-Nmed:i+Nmed,j-Nmed:j+Nmed) == sym%YES) / Nmed_Total)
         I1_Temp_2(i,j) = nint(count(Smoke_Mask(i-Nmed:i+Nmed,j-Nmed:j+Nmed) == sym%YES) / Nmed_Total)
      end do elem_loop_median
   end do line_loop_median
   Dust_Mask = I1_Temp_1
   Smoke_Mask = I1_Temp_2
   deallocate(I1_Temp_1)
   deallocate(I1_Temp_2)

   !-------------------------------------------------------------------------------------------------------
   !--- save dust, smoke, fire to the corresponding bit structure 
   !-------------------------------------------------------------------------------------------------------
   line_loop_pack: do i = 1, Image%Number_Of_Elements
      elem_loop_pack: do  j = 1, Image%Number_Of_Lines_Read_This_Segment
           if (Smoke_Mask(i,j) == 1) Cld_Test_Vector_Packed(2,i,j) = ibset(Cld_Test_Vector_Packed(2,i,j),4)
           if (Dust_Mask(i,j)  == 1) Cld_Test_Vector_Packed(2,i,j) = ibset(Cld_Test_Vector_Packed(2,i,j),5)
           if (Fire_Mask(i,j)  == 1) Cld_Test_Vector_Packed(2,i,j) = ibset(Cld_Test_Vector_Packed(2,i,j),7)
      end do elem_loop_pack
   end do line_loop_pack

   !------------------------------------------------------------------------------
   !-- CLAVR-x specific Processing
   !--- grab version tags for output as attributes in level2
   !--- only need to do this once, so do on first segment
   !------------------------------------------------------------------------------
   if (Segment_Number == 1) then
     call SET_CLOUD_MASK_VERSION(Cloud_Mask_Version)
     call SET_CLOUD_MASK_THRESHOLDS_VERSION(Cloud_Mask_Thresholds_Version)
   endif

   !-------------------------------------------------------------------------------
   ! on last segment, wipe out the lut from memory and reset is_read_flag to no
   !-------------------------------------------------------------------------------
   if (Segment_Number == Input%Num_Segments) then
       call RESET_NB_CLOUD_MASK_LUT()
   endif

   First_Call = .false.

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
   ! set input
   !============================================================================
   subroutine SET_INPUT(i,j)
      integer, intent (in) :: i, j

      Input%Num_Elem = Image%Number_Of_Elements
      Input%Num_Line = Image%Number_Of_Lines_Read_This_Segment
      Input%Num_Line_Max = Image%Number_Of_Lines_Per_Segment
      Input%Num_Segments = Image%Number_Of_Segments
      Input%Invalid_Data_Mask = Bad_Pixel_Mask(i,j)
      Input%Chan_On_041um = Sensor%Chan_On_Flag_Default(8)
      Input%Chan_On_063um = Sensor%Chan_On_Flag_Default(1)
      Input%Chan_On_086um = Sensor%Chan_On_Flag_Default(2)
      Input%Chan_On_138um = Sensor%Chan_On_Flag_Default(26)
      Input%Chan_On_160um = Sensor%Chan_On_Flag_Default(6)
      Input%Chan_On_213um = Sensor%Chan_On_Flag_Default(7)
      Input%Chan_On_375um = Sensor%Chan_On_Flag_Default(20)
      Input%Chan_On_67um = Sensor%Chan_On_Flag_Default(27)
      Input%Chan_On_85um = Sensor%Chan_On_Flag_Default(29)
      Input%Chan_On_10um = Sensor%Chan_On_Flag_Default(38)
      Input%Chan_On_11um = Sensor%Chan_On_Flag_Default(31)
      Input%Chan_On_12um = Sensor%Chan_On_Flag_Default(32)
      Input%Chan_On_I1_064um = Sensor%Chan_On_Flag_Default(39)
      Input%Chan_On_I4_374um = Sensor%Chan_On_Flag_Default(42)
      Input%Chan_On_I5_114um = Sensor%Chan_On_Flag_Default(43)
      Input%Chan_On_DNB = Sensor%Chan_On_Flag_Default(44)
      Input%Use_Sounder_11um = 0                        !note, off by default
      Input%Snow_Class = Sfc%Snow(i,j)
      Input%Land_Class = Sfc%Land(i,j)
      Input%Oceanic_Glint_Mask = Sfc%Glint_Mask(i,j)
      Input%Coastal_Mask = Sfc%Coast_Mask(i,j)
      Input%Solzen = Geo%Solzen(i,j)
      Input%Scatzen = Geo%Scatangle(i,j)
      Input%Senzen = Geo%Satzen(i,j)
      Input%Lat = Nav%Lat(i,j)
      Input%Lon = Nav%Lon(i,j)
      Input%Sst_Anal_Uni = Sst_Anal_Uni(i,j)
      Input%Emiss_Sfc_375um = ch(20)%Sfc_Emiss(i,j)    !note, this is on always
      Input%Zsfc = Sfc%Zsfc(i,j)
      Input%Solar_Contamination_Mask = Solar_Contamination_Mask(i,j)
      Input%Sfc_Type = Sfc%Sfc_Type(i,j)

      if (Input%Chan_On_041um == sym%YES)  then 
        Input%Ref_041um = ch(8)%Ref_Toa(i,j)
      endif
      if (Input%Chan_On_063um == sym%YES)  then 
        Input%Ref_063um = ch(1)%Ref_Toa(i,j)
        Input%Ref_063um_Clear = ch(1)%Ref_Toa_Clear(i,j)
        Input%Ref_063um_Std = Ref_Ch1_Std_3x3(i,j)
        Input%Ref_063um_Min = Ref_Ch1_Min_3x3(i,j)
      endif
      if (Input%Chan_On_086um == sym%YES)  then 
        Input%Ref_086um = ch(2)%Ref_Toa(i,j)
      endif
      if (Input%Chan_On_138um == sym%YES)  then 
        Input%Ref_138um = ch(26)%Ref_Toa(i,j)
      endif
      if (Input%Chan_On_160um == sym%YES)  then 
        Input%Ref_160um = ch(6)%Ref_Toa(i,j)
        Input%Ref_160um_Clear = ch(6)%Ref_Toa_Clear(i,j)
      endif
      if (Input%Chan_On_213um == sym%YES)  then 
        Input%Ref_213um = ch(7)%Ref_Toa(i,j)
      endif
      if (Input%Chan_On_375um == sym%YES)  then 
        Input%Ref_375um = ch(20)%Ref_Toa(i,j)
        Input%Ref_375um_Clear = ch(20)%Ref_Toa_Clear(i,j)
        Input%Bt_375um = ch(20)%Bt_Toa(i,j)
        Input%Bt_375um_Std = Bt_Ch20_Std_3x3(i,j)
        Input%Emiss_375um =  Ems_Ch20_Median_3x3(i,j)
        Input%Emiss_375um_Clear = Ems_Ch20_Clear_Rtm(i,j)
      endif
      if (Input%Chan_On_67um == sym%YES)  then 
        Input%Bt_67um = ch(27)%Bt_Toa (i,j)
        if (Input%Chan_On_11um == sym%YES)  then 
           Input%Bt_11um_Bt_67um_Covar = Covar_Ch27_Ch31_5x5(i,j)
        endif
      endif
      if (AVHRR_IFF_Flag == 1)  then
        Input%Use_Sounder_11um = sym%YES
        Input%Bt_11um_Sounder = Bt_11um_Sounder(i,j)
        Input%Bt_11um_Bt_67um_Covar = Missing_Value_Real4
      endif
      if (Input%Chan_On_85um == sym%YES)  then 
        Input%Bt_85um = ch(29)%Bt_Toa(i,j)
      endif
      if (Input%Chan_On_10um == sym%YES)  then                                                                                                      
        Input%Bt_10um = ch(38)%Bt_Toa(i,j)
      else
        Input%Bt_10um = Missing_Value_Real4
      endif
      if (Input%Chan_On_11um == sym%YES)  then 
        Input%Bt_11um = ch(31)%Bt_Toa(i,j)
        Input%Bt_11um_Std = Bt_Ch31_Std_3x3(i,j)
        Input%Bt_11um_Max = Bt_Ch31_Max_3x3(i,j)
        Input%Bt_11um_Clear = ch(31)%Bt_Toa_Clear(i,j)
        Input%Emiss_11um_Tropo = ch(31)%Emiss_Tropo(i,j)
      endif
      if (Input%Chan_On_12um == sym%YES)  then 
        Input%Bt_12um = ch(32)%Bt_Toa(i,j)
        Input%Bt_12um_Clear = ch(32)%Bt_Toa_Clear(i,j)
      endif
      if (Input%Chan_On_DNB == sym%YES)  then 
        Input%Lunscatzen = Geo%Scatangle_Lunar(i,j)
        Input%Lunar_Oceanic_Glint_Mask = Sfc%Glint_Mask_Lunar(i,j)
        Input%Rad_Lunar = ch(44)%Rad_Toa(i,j)
        Input%Ref_Lunar = ch(44)%Ref_Lunar_Toa(i,j)
        Input%Ref_Lunar_Min = Ref_ChDNB_Lunar_Min_3x3(i,j)
        Input%Ref_Lunar_Std = Ref_ChDNB_Lunar_Std_3x3(i,j)
        Input%Ref_Lunar_Clear = ch(44)%Ref_Lunar_Toa_Clear(i,j)
        Input%Lunzen = Geo%Lunzen(i,j)
      endif
      if (Input%Chan_On_I1_064um == sym%YES) then
         Input%Ref_I1_064um_Std = Ref_Uni_ChI1(i,j)
      endif
      if (Input%Chan_On_I4_374um == sym%YES) then
         Input%Bt_I4_374um_Std = Bt_Uni_ChI4(i,j)
      endif
      if (Input%Chan_On_I5_114um == sym%YES) then
         Input%Bt_I5_114um_Std = Bt_Uni_ChI5(i,j)
      endif
   end subroutine SET_INPUT

   subroutine SET_OUTPUT(i,j)
      integer, intent (in) :: i, j

      Cld_Test_Vector_Packed(:,i,j) = Output%Cld_Flags_Packed
      Cld_Mask(i,j) = Output%Cld_Mask_Bayes
      Posterior_Cld_Probability(i,j) = Output%Posterior_Cld_Probability
      Dust_Mask(i,j) = Output%Dust_Mask
      Smoke_Mask(i,j) = Output%Smoke_Mask
      Fire_Mask(i,j) = Output%Fire_Mask
   end subroutine SET_OUTPUT

   subroutine SET_DIAG(i,j)
      integer, intent (in) :: i, j

      Diag_Pix_Array_1(i,j) = Diag%Array_1
      Diag_Pix_Array_2(i,j) = Diag%Array_2
      Diag_Pix_Array_3(i,j) = Diag%Array_3
   end subroutine SET_DIAG

!==============================================================
! Median filter
!
! mask = 0 means use median, mask = 1 mean ignore
!==============================================================
subroutine MEDIAN_FILTER_MASK(z,mask,z_median)

! The purpose of this function is to find 
! median (emed), minimum (emin) and maximum (emax)
! for the array elem with nelem elements. 

  real, dimension(:,:), intent(in):: z
  real, intent(out):: z_median
  integer(kind=int1), dimension(:,:), intent(in):: mask
  integer:: i,j,k,nx,ny,nelem
  real, dimension(:), allocatable::x
  real:: u

  z_median = missing_value_real4

  nx = size(z,1)
  ny = size(z,2)

  nelem = nx * ny

  allocate(x(nelem))
  x = 0.0
  k = 0
  do i = 1, nx
    do j = 1, ny
      if (mask(i,j) == 0 .and. z(i,j) /= missing_value_int1) then
           k = k + 1
           x(k) = z(i,j)
      endif
   enddo
  enddo

  nelem = k

  if (nelem < 1) then
     if (allocated(x)) deallocate(x)
     return
  endif
  !--- sort the array into ascending order
  do i=1,nelem-1
   do j=i+1,nelem
    if(x(j)<x(i))then
     u=x(j)
     x(j)=x(i)
     x(i)=u
    end if
   end do
  end do

  !---- pick the median
  if(mod(nelem,2)==1)then
   i=nelem/2+1
   z_median=x(i)
  else
   i=nelem/2
   z_median=(x(i)+x(i+1))/2
   end if

  if (allocated(x)) deallocate(x)

end subroutine MEDIAN_FILTER_MASK

!============================================================================

end module NB_CLOUD_MASK_CLAVRX_BRIDGE
