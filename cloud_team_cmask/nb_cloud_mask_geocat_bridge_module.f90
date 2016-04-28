!$Id: 
!-------------------------------------------------------------------------------
!
! PURPOSE: Bridge for GEOCAT for NB Cloud mask. 
!
! IMPORTANT - WILL NOT WORK UNTIL GEOCAT HAS NETCDF CAPABILITY
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
!-------------------------------------------------------------------------------
module NB_CLOUD_MASK_GEOCAT_BRIDGE_MODULE


   use NB_CLOUD_MASK_GEOCAT_SERVICES_MODULE
   use NB_CLOUD_MASK

   implicit none

   public :: NB_CLOUD_MASK_BRIDGE

   private :: SET_SYMBOL
   private :: SET_INPUT
   private :: SET_OUTPUT
   private :: SET_DIAG
   private :: NULL_OUTPUT
   private :: NULL_DIAG

   !--- define these structure as module wide
   type(mask_input), private :: Input   
   type(mask_output), private :: Output   
   type(diag_output), private :: Diag  
   type(symbol_naive_bayesian),private :: Symbol
   
   !Make module wide variables
   character (len=120), TARGET, PRIVATE:: Ancil_Data_Path
   character (len=120), TARGET, PRIVATE:: Naive_Bayes_File_Name
   REAL(kind=real4), TARGET, PRIVATE :: Covar_Ch27_Ch31_5x5

   !Segment counter
   integer(kind=INT1), TARGET, PRIVATE:: Segment_Number_CM = 1

   !Make Iband and DNB flag
   integer (kind=INT1), DIMENSION(5), PRIVATE :: IBand_Flag
   integer (kind=INT1), PRIVATE :: DNB_Flag
   
   !allocatable
   integer (kind=INT1), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Solar_Contamination_Mask

   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ref_Ch1_Mean_3X3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ref_Ch1_Max_3x3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ref_Ch1_Min_3X3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ref_Ch1_Stddev_3X3

   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch31_Mean_3x3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch31_Max_3x3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch31_Min_3x3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch31_Stddev_3x3


   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch20_Stddev_3x3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch20_Mean_3x3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch20_Max_3x3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch20_Min_3x3


   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ems_Ch20_Median_3x3
   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ems_39_Med_3x3

   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ems_Ch20_Std_Median_3x3

   REAL(kind=real4), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Covar_67_11_5x5

   !Glint mask
   integer (kind=INT1), TARGET, PRIVATE :: Glint_Mask
      
   !Pointers
   INTEGER(kind=INT1), DIMENSION(:), POINTER, PRIVATE :: CHN_FLG
   REAL(KIND=REAL4), DIMENSION(:,:), POINTER :: Emiss_11um_Tropo_Rtm

!   INTEGER(LONG) :: Sfc_Idx_NWP
   REAL(kind=REAL4), PRIVATE :: Chn7_Sol_Energy

   
contains
!----------------------------------------------------------------------
! Bridge Routine
!
! Note, the Diag argument is optional
!---------------------------------------------------------------------- 
 subroutine NB_CLOUD_MASK_BRIDGE(Algo_Idx)
 
   implicit none

   INTEGER(KIND=INT4), INTENT(IN) :: Algo_Idx
   integer :: Line_Idx, Elem_Idx
   integer:: Num_Elem
   integer:: Num_Line
   integer:: Num_Line_Max
   integer:: Elem_Idx_min
   integer:: Elem_Idx_max
   integer:: Elem_Idx_width
   integer:: Elem_Idx_segment_max
   integer:: Line_Idx_min
   integer:: Line_Idx_max
   integer:: Line_Idx_width
   integer:: Line_Idx_segment_max
   integer:: VIIRS_375M_res_indx
   integer :: McIDAS_ID
   REAL(kind=real4) :: Glint_Zen_Thresh=40.0
   character (len=555):: Naive_Bayes_File_Name_Full_Path, Bayesian_Cloud_Mask_Name
   CHARACTER(LEN=1024):: Lut_Path



   !---- set paths and mask classifier file name to their values in this framework


  ! --   Set LUT file
  Lut_Path = TRIM(Out2(Algo_Idx)%Ancil_Subdir)//"/"

   !---- FIXME - need to figure out how to set bayesian name
   Naive_Bayes_File_Name_Full_Path = trim(Lut_Path)//trim(Bayesian_Cloud_Mask_Name)
         
   Num_Elem = sat%nx
   Num_Line = sat%ny
   Num_Line_Max = size(sat%lat,2)

   !allocate local arrays
   allocate(Ems_Ch20_Std_Median_3x3(num_elem,num_line))
   allocate(Ems_39_Med_3x3(num_elem,num_line))
   allocate(Covar_67_11_5x5(num_elem,num_line))

   !Solar Contamination
   allocate(Solar_Contamination_Mask(num_elem, num_line))
   !Only needed for AVHRR Counts
   Solar_Contamination_Mask = sym%NO
      
   !Point to channel flag
   CHN_FLG => out2(Algo_Idx)%ch_flg

   !--- I-Band Flag
   IBand_Flag(:) = sym%NO

   !--- DNB reflectance 
   DNB_Flag = sym%NO

   !--- set structure (symbol, input, output, diag)  elements to corresponding values in this framework
   call SET_SYMBOL()

   !Compute Spatial uniformity


   CALL Compute_Spatial_Uniformity(1, 1, sat%space_mask, sat%bt14, &
                                   Bt_Ch31_Mean_3x3, &
                                   Bt_Ch31_Max_3x3, Bt_Ch31_Min_3x3, &
                                   Bt_Ch31_Stddev_3x3)
   
   CALL Compute_Spatial_Uniformity(1, 1, sat%space_mask, sat%ref2, &
                                   Ref_Ch1_Mean_3X3, &
                                  Ref_Ch1_Max_3x3, Ref_Ch1_Min_3X3, &
                                   Ref_Ch1_Stddev_3X3)

   CALL Compute_Spatial_Uniformity(1, 1, sat%space_mask, sat%bt7, &
                                   Bt_Ch20_Mean_3x3, &
                                    Bt_Ch20_Max_3x3, Bt_Ch20_Min_3x3, &
                                   Bt_Ch20_Stddev_3x3)
        
    ! Calculate 3.9um emissivity median

   call COMPUTE_MEDIAN_SEGMENT(sat%ems7, sat%Bad_Pixel_Mask(14,:,:), 1, &
                              1, Num_Elem, 1, &
                              Num_line, &
                              Ems_39_Med_3x3,  &
                              Ems_Ch20_Std_Median_3x3)

   !=======================================================================
   ! Compute 11 micron emissivity at Tropopause
   !=======================================================================
   Emiss_11um_Tropo_Rtm => out2(Algo_Idx)%emiss11_high
   CALL Compute_Emiss_11um_Tropo_Rtm(Emiss_11um_Tropo_Rtm, &
                                  Num_Line)    
    
    ! -----------    loop over pixels -----   
   line_loop: do Line_Idx = 1, Num_Line
      elem_loop: do  Elem_Idx = 1, Num_Elem
      
        !-------------------------------------------------------------------
        ! Do space mask check here
        !-------------------------------------------------------------------

        IF (sat%space_mask(Elem_Idx,Line_Idx) == sym%YES) THEN
            out2(Algo_Idx)%cldmask(Elem_Idx,Line_Idx) = missing_value_int1
            CYCLE
        ENDIF
       



       
        !-------------------------------------------------------------------
        ! Do glint mask here
        !-------------------------------------------------------------------
      
        !--- initialize valid pixel to no
        Glint_Mask = sym%NO

        !--- skip land pixels
        IF ((sat%Land_Mask(Elem_Idx,Line_Idx) == sym%NO) .and. &
             sat%Snow_mask(Elem_Idx,Line_Idx) == sym%NO_SNOW) THEN

       !--- turn on in geometric glint cone and sufficient Ref_Ch1
            IF ((sat%glintzen(Elem_Idx,Line_Idx) < Glint_Zen_Thresh)) THEN

            !--- assume to be glint IF in geometric zone
                Glint_Mask = sym%YES

                IF (CHN_FLG(14) == sym%YES) THEN

                !--- exclude pixels colder than the freezing temperature
                    IF (sat%bt14(Elem_Idx,Line_Idx) < 273.15) THEN
                        Glint_Mask = sym%NO
                    ENDIF

            !--- exclude pixels colder than the surface
                    IF (sat%bt14(Elem_Idx,Line_Idx) < &
                        sat%bt_clr14(Elem_Idx,Line_Idx) - 5.0) THEN
                        Glint_Mask = sym%NO
                    ENDIF

                ENDIF

          !-turn off IF non-unIForm in reflectance
                IF (CHN_FLG(2) == sym%YES) THEN
                    IF (Ref_Ch1_Stddev_3X3(Elem_Idx,Line_Idx) > 1.0) THEN
                        Glint_Mask = sym%NO
                    ENDIF
                ENDIF

            ENDIF  !Glintzen check

        ENDIF    !land check



        !-------------------------------------------------------------------
        ! Do covariance here
        !-------------------------------------------------------------------
        Elem_Idx_min = max(1,min(Elem_Idx - 2,Num_Elem))
        Elem_Idx_max = max(1,min(Elem_Idx + 2,Num_Elem))
        Line_Idx_min = max(1,min(Line_Idx - 2,Num_line))
        Line_Idx_max = max(1,min(Line_Idx + 2,Num_line))
        Line_Idx_width = Line_Idx_max - Line_Idx_min + 1
        Elem_Idx_width = Elem_Idx_max - Elem_Idx_min + 1

        Covar_67_11_5x5(Elem_Idx,Line_Idx) = Missing_Value_Real4

         IF ((CHN_FLG(9) == sym%YES) .and. &
            ( CHN_FLG(14)== sym%YES)) THEN

            Covar_67_11_5x5(Elem_Idx,Line_Idx) = Covariance_local(&
              sat%bt14(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
              sat%bt9(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
              Elem_Idx_width, Line_Idx_width, &
               sat%Bad_Pixel_Mask(14,Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max))
        ENDIF
      

        
          
      
         ! Set inputs
         
         
         call SET_INPUT(Elem_Idx,Line_Idx)
         call SET_OUTPUT(Elem_Idx,Line_Idx, Algo_Idx)
         call SET_DIAG(Elem_Idx,Line_Idx)

         !---call cloud mask routine
!         call NB_CLOUD_MASK_ALGORITHM( &
!                      Naive_Bayes_File_Name_Full_Path, &
!                      Symbol,  &
!                      Input, &
!                      Output)

         !--- nullify pointers within these data structures
         call NULL_OUTPUT()
         call NULL_DIAG()
   
      end do elem_loop
   end do line_loop

!-------------------------------------------------------------------------------
! on last segment, wipe out the lut from memory and reset is_read_flag to no
!-------------------------------------------------------------------------------
!   if (Segment_Number_CM == Input%Num_Segments) then
!       call RESET_NB_CLOUD_MASK_LUT()
!   endif


   
   
   !Deallocate arrays
   deallocate(Solar_Contamination_Mask)
   deallocate(Ems_Ch20_Std_Median_3x3)

   CALL Destroy_Spatial_Uniformity(Ref_Ch1_Mean_3X3, Ref_Ch1_Max_3x3,  &
                                   Ref_Ch1_Min_3X3, Ref_Ch1_Stddev_3X3)
                                   
   CALL Destroy_Spatial_Uniformity(Bt_Ch31_Mean_3x3, Bt_Ch31_Max_3x3,  &
                                   Bt_Ch31_Min_3x3, Bt_Ch31_Stddev_3x3)
                                   
   CALL Destroy_Spatial_Uniformity(Bt_Ch20_Mean_3x3, Bt_Ch20_Max_3x3,  &
                                   Bt_Ch20_Min_3x3, Bt_Ch20_Stddev_3x3)
   
      
   !Pointers
   CHN_FLG => null()
   Emiss_11um_Tropo_Rtm =>null()
 
   
   !Increment segment number
   Segment_Number_CM = Segment_Number_CM +1

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

      Input%Num_Elem = sat%nx
      Input%Num_Line = sat%ny
      Input%Num_Line_Max = size(sat%lat,2)
      Input%Num_Segments = Sat%Iseg
      !------
      Input%Invalid_Data_Mask = sat%Bad_Pixel_Mask(14,i,j)
      Input%Chan_On_041um = CHN_FLG(1)
      Input%Chan_On_063um = CHN_FLG(2)
      Input%Chan_On_086um = CHN_FLG(3)
      Input%Chan_On_138um = CHN_FLG(4)
      Input%Chan_On_160um = CHN_FLG(5)
      Input%Chan_On_213um = CHN_FLG(6)
      Input%Chan_On_375um = CHN_FLG(7)
      Input%Chan_On_67um = CHN_FLG(9)
      Input%Chan_On_85um = CHN_FLG(11)
      Input%Chan_On_11um = CHN_FLG(14)
      Input%Chan_On_12um = CHN_FLG(15)
      Input%Chan_On_I1_064um = IBand_Flag(1)
      Input%Chan_On_I4_374um = IBand_Flag(4)
      Input%Chan_On_I5_114um = IBand_Flag(5)
      Input%Chan_On_DNB = DNB_Flag
      Input%Use_Sounder_11um = sym%NO

      Input%Bt_11um_Sounder = MISSING_VALUE_REAL4

      Input%Coastal_Mask = sat%Coast_Mask(i,j)
      Input%Snow_Class = sat%Snow_mask(i,j)
      Input%Land_Class = sat%Land_Mask(i,j)
      Input%Sfc_Type = sat%Sfc_Type(i,j)
      
      
      Input%Oceanic_Glint_Mask = Glint_Mask
      Input%Solzen = sat%SolZen(i,j)
      Input%Scatzen = sat%scatzen(i,j)
      Input%Senzen = sat%Satzen(i,j)
      Input%Lunzen = MISSING_VALUE_INT4 !No DNB in GEOCAT now
      Input%Lat = sat%Lat(i,j)
      Input%Lon = sat%Lon(i,j)
      
      IF (Input%Chan_On_063um == sym%YES) THEN
        Input%Ref_063um = sat%ref2(i,j)
        Input%Ref_063um_Clear = sat%ws_albedo2(i,j)
        Input%Ref_063um_Std = Ref_Ch1_Stddev_3X3(i,j)!FIXME
        Input%Ref_063um_Min = Ref_Ch1_Min_3X3(i,j)!FIXME     
      ENDIF
      
      IF (Input%Chan_On_086um == sym%YES) Input%Ref_086um = sat%ref3(i,j)
      IF (Input%Chan_On_138um == sym%YES) Input%Ref_138um = sat%ref4(i,j)

      IF (Input%Chan_On_160um == sym%YES) THEN
        Input%Ref_160um = sat%ref5(i,j)
        Input%Ref_160um_Clear = sat%ws_albedo5(i,j)
      ENDIF


      IF (Input%Chan_On_213um == sym%YES) Input%Ref_213um = sat%ref6(i,j)

      IF (Input%Chan_On_375um == sym%YES) THEN
        Input%Ref_375um = sat%ref7(i,j)
        Input%Ref_375um_Clear = MISSING_VALUE_REAL4 !Not filled or used for now
        Input%Bt_375um = sat%bt7(i,j)
        Input%Bt_375um_Std = Bt_Ch20_Stddev_3x3(i,j) 
        Input%Emiss_375um =  Ems_39_Med_3x3(i,j) 
        Input%Emiss_375um_Clear = sat%ems_sol_clr7(i,j)
        Input%Emiss_Sfc_375um = sat%sfc_emiss7(i,j)
      ENDIF
      
      IF (Input%Chan_On_67um == sym%YES) THEN
        Input%Bt_67um = sat%bt9(i,j)
        IF (Input%Chan_On_11um == sym%YES) Input%Bt_11um_Bt_67um_Covar = Covar_67_11_5x5(i,j) 
      
      ENDIF
      
      
      IF (Input%Chan_On_85um == sym%YES) Input%Bt_85um = sat%bt11(i,j)
      
      IF (Input%Chan_On_11um == sym%YES)  THEN
        Input%Bt_11um = sat%bt14(i,j)
        Input%Bt_11um_Std = Bt_Ch31_Stddev_3x3(i,j) 
        Input%Bt_11um_Max = Bt_Ch31_Max_3x3(i,j)
        Input%Bt_11um_Clear = sat%bt_clr14(i,j)
        Input%Emiss_11um_Tropo = Emiss_11um_Tropo_Rtm(i,j)
      ENDIF
      
      IF (Input%Chan_On_12um == sym%YES) THEN
        Input%Bt_12um = sat%bt15(i,j)
        Input%Bt_12um_Clear = sat%bt_clr15(i,j)
      ENDIF
     
      Input%Sst_Anal_Uni = sat%sst_clim_uni(i,j)
      Input%Zsfc = sat%Zsfc(i,j)
      Input%Solar_Contamination_Mask = Solar_Contamination_Mask(i,j)
   end subroutine SET_INPUT

   subroutine SET_OUTPUT(i,j, Algo_Idx)
      integer, intent (in) :: i, j
      INTEGER(KIND=INT4), INTENT(IN) :: Algo_Idx

      Output%Cld_Flags_Packed => out2(Algo_Idx)%cldmask_packed(:,i,j)
      Output%Cld_Mask_Bayes => out2(Algo_Idx)%cldmask(i,j)
      Output%Posterior_Cld_Probability => out2(Algo_Idx)%cld_probability(i,j)
      Output%Dust_Mask => null()
      Output%Smoke_Mask => null()
      Output%Fire_Mask => null()
   end subroutine SET_OUTPUT

   subroutine SET_DIAG(i,j)
      integer, intent (in) :: i, j

      Diag%Array_1 => null()
      Diag%Array_2 => null()
      Diag%Array_3 => null()
   end subroutine SET_DIAG

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


!====================================================================
! Subroutine Name: Compute_Emiss_11um_Tropo_Rtm
!
! Function:
!   Computes the 11 micron emissivity at Tropopause
!
! Description:
!   This subroutine computes 11 micron emissivity at Tropopause for a segment of
!   data
!
! Calling Sequence:
!   CALL Compute_Emiss_11um_Tropo_Rtm(Emiss_11um_Tropo_Rtm, &
!                                   Number_of_Lines_in_this_Segment)
!
! Inputs: 
!  Number_of_Lines_in_this_Segment - Number of lines in this segment of data
! 
! Outputs: 
!  Emiss_11um_Tropo_Rtm - Channel 14 emissivity at the tropopause
!
! Dependencies: None
!
! Restrictions:  None
!
!====================================================================
 SUBROUTINE Compute_Emiss_11um_Tropo_Rtm(Emiss_11um_Tropo_Rtm,&
                                      Number_of_Lines_in_this_Segment)
   INTEGER(KIND=INT4),INTENT(IN) :: Number_of_Lines_in_this_Segment
   REAL(KIND=REAL4), DIMENSION(:,:), INTENT(OUT):: Emiss_11um_Tropo_Rtm
   INTEGER(KIND=INT1):: Tropo_Idx_NWP
   INTEGER(KIND=INT1):: View_Zen_Idx
   INTEGER:: X_NWP_Idx
   INTEGER:: Y_NWP_Idx
   INTEGER:: Elem_Idx
   INTEGER:: Line_Idx
   REAL(KIND=REAL4) :: Rad_Chn14
   REAL(KIND=REAL4) :: Clr_Rad_Chn14
   REAL(KIND=REAL4) :: Blkbdy_Tropo_Rad_Chn14

   !--- initialize
   Emiss_11um_Tropo_Rtm = Missing_Value_Real4

    Line_Loop: DO Line_Idx=1, Number_of_Lines_in_this_Segment
      Element_Loop: DO Elem_Idx = 1, sat%nx

       IF (sat%space_mask(Elem_Idx,Line_Idx) == sym%NO) THEN

            !
            !---nwp longitude cell
            !
            X_NWP_Idx =          sat%x_nwp(Elem_Idx,Line_Idx)         

            !
            !---nwp latitude cell
            !
            Y_NWP_Idx =          sat%y_nwp(Elem_Idx,Line_Idx)     

            !
            !---nwp level associated with tropopause
            !
            Tropo_Idx_NWP =        nwp%dat(X_NWP_Idx,Y_NWP_Idx)%Tropo_Level 

            !
            !---viewing zenith angle bin
            !
            View_Zen_Idx =          sat%ivza(Elem_Idx,Line_Idx)        

            !
            !---11 um radiance
            !
            Rad_Chn14  =        sat%rad14(Elem_Idx,Line_Idx)

            !
            !---clear 11 micron radiance
            !
            Clr_Rad_Chn14 =     sat%rad_clr14(Elem_Idx,Line_Idx)

            !
            !---BB 11 um rad at tropopause
            !
            Blkbdy_Tropo_Rad_Chn14 = rtm(X_NWP_Idx,Y_NWP_Idx)%d(View_Zen_Idx)%cloud_prof14(Tropo_Idx_NWP)

            !
            !---Tropopause Emissivity
            !
            Emiss_11um_Tropo_Rtm(Elem_Idx,Line_Idx) =  &
                  (Rad_Chn14 - Clr_Rad_Chn14) / (Blkbdy_Tropo_Rad_Chn14 - Clr_Rad_Chn14)

      END IF
    END DO Element_Loop
  END DO Line_Loop

 END SUBROUTINE Compute_Emiss_11um_Tropo_Rtm




!---- MEDIAN Routine should go in num_mod

!==============================================================
! subroutine COMPUTE_MEDIAN(z,mask,z_median,z_mean,z_std_median)
!
! Median filter
!==============================================================
subroutine COMPUTE_MEDIAN(z,mask,z_median,z_mean,z_std_median)

! The purpose of this function is to find
! median (emed), minimum (emin) and maximum (emax)
! for the array elem with nelem elements.

 real, dimension(:,:), intent(in):: z
 real, intent(out):: z_median
 real, intent(out):: z_mean
 real, intent(out):: z_std_median
 integer(kind=int1), dimension(:,:), intent(in):: mask
 integer:: i,j,k,nx,ny,nelem
 real, dimension(:), allocatable::x
 real(kind=real4):: u

 z_median = missing_value_real4
 z_std_median = missing_value_real4
 z_mean = missing_value_real4

 nx = size(z,1)
 ny = size(z,2)

 nelem = nx * ny

 allocate(x(nelem))
 x = 0.0

 k = 0
 do i = 1, nx
   do j = 1, ny
      if (mask(i,j) == sym%NO .and. z(i,j) /= missing_value_real4) then
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

!--- compute standard deviation wrt median
  z_mean = sum(x(1:nelem))/nelem
  z_std_median = sqrt(sum((x(1:nelem) - z_median)**2) / nelem)


! if (z_std_median > 60.0) then
!         print *, "big std median ", z_std_median, nelem, x(1:nelem)
!         print *, "z_nxn = ", z
! endif

  if (allocated(x)) deallocate(x)

end subroutine COMPUTE_MEDIAN

!----------------------------------------------------------------------
! subroutine COMPUTE_MEDIAN_SEGMENT(z,mask,n,imin,imax,jmin,jmax,
!                                   z_median,z_std_median)
!
! Compute standard deviaion of an array wrt to the median
!----------------------------------------------------------------------
subroutine COMPUTE_MEDIAN_SEGMENT(z,mask,n,imin,imax,jmin,jmax, &
                                  z_median, &
                                  z_std_median)
  real(kind=real4), dimension(:,:), intent(in):: z
  integer(kind=int1), dimension(:,:), intent(in):: mask
  real(kind=real4), dimension(:,:), intent(out):: z_std_median
  real(kind=real4), dimension(:,:), intent(out):: z_median
! real(kind=real4), dimension(:,:), intent(out):: z_mean
  integer, intent(in):: n
  integer, intent(in):: imin
  integer, intent(in):: imax
  integer, intent(in):: jmin
  integer, intent(in):: jmax
  integer:: i
  integer:: j
  integer:: i1
  integer:: i2
  integer:: j1
  integer:: j2
  real(kind=real4) :: z_mean

  do i = imin, imax
    do j = jmin, jmax

     j1 = max(jmin,j-n)   !top index of local array
     j2 = min(jmax,j+n)   !bottom index of local array
     i1 = max(imin,i-n)   !left index of local array
     i2 = min(imax,i+n)   !right index of local array

!--- compute median
     call COMPUTE_MEDIAN(z(i1:i2,j1:j2),mask(i1:i2,j1:j2),z_median(i,j), &
                         z_mean,z_std_median(i,j))

     enddo
  enddo

end subroutine COMPUTE_MEDIAN_SEGMENT




end module NB_CLOUD_MASK_GEOCAT_BRIDGE_MODULE
