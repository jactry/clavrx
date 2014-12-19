!$Id:$
!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
module NB_CLOUD_MASK_SAPF_BRIDGE


   ! -- Framework specific modules
   USE PCF_NPP_BAYES_CLOUD_MASK_Mod
   USE TYPE_KINDS_AIT
   USE Convert_Char
   USE Error_Messaging_Module
   USE Framework_Global_Variables_Module
   USE LandMask_Access_Mod
   USE Sat_Access_Mod
   USE CoastMask_Access_Mod
   USE SnowMask_Access_Mod
   USE SfcEmis_Access_Mod
   USE DesertMask_Access_Mod
   USE CloudMask_Access_Mod
   USE SfcElev_Access_Mod
   USE SfcAlbd_Access_Mod
   USE SST_Access_Mod
   USE SfcType_Access_Mod
   USE PseudoEmiss_Access_Mod
   USE NWP_Access_Mod
   USE Numerical_Routines
   USE RTM_Data_Access_Mod
   USE RTM_Access_Mod
   !USE RTM_MODULE
   USE Planck_Module
   USE NF_PARM
   ! -- MODULES USED
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
   private :: COMPUTE_IBAND_STATISTICS

   !--- define these structure as module wide
   type(mask_input), private :: Input   
   type(mask_output), private :: Output   
   type(diag_output), private :: Diag  
   type(symbol_naive_bayesian),private :: Symbol
   
   !Make module wide variables
   character (len=120), TARGET, PRIVATE:: Ancil_Data_Path
   character (len=120), TARGET, PRIVATE:: Naive_Bayes_File_Name

   !Segment counter
   integer(kind=INT1), TARGET, PRIVATE:: Segment_Number_CM = 1

   !Make Iband and DNB flag
   integer (kind=INT1), DIMENSION(5), PRIVATE :: IBand_Flag
   integer (kind=INT1), PRIVATE :: DNB_Flag
   
   !allocatable
   integer (kind=INT1), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Solar_Contamination_Mask

   REAL (kind=SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Diag_Pix_Array_1
   REAL (kind=SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Diag_Pix_Array_2
   REAL (kind=SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Diag_Pix_Array_3

   integer (kind=INT1), TARGET, PRIVATE :: Glint_Mask
   REAL(SINGLE), TARGET, PRIVATE :: Covar_Ch27_Ch31_5x5
   
   
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ref_Ch1_Mean_3X3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ref_Ch1_Max_3x3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ref_Ch1_Min_3X3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ref_Ch1_Stddev_3X3

   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch31_Mean_3x3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch31_Max_3x3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch31_Min_3x3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch31_Stddev_3x3


   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ems_Ch20_Median_3x3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Ems_Ch20_Std_Median_3x3


   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch20_Stddev_3x3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch20_Mean_3x3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch20_Max_3x3
   REAL(SINGLE), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Bt_Ch20_Min_3x3


   !I-Band uniformity
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Dummy_IBand
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Ref_Uni_ChI1
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Bt_Uni_ChI4
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Bt_Uni_ChI5
   
   !Pointers
   INTEGER(LONG), DIMENSION(:), POINTER, PRIVATE :: CHN_FLG
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Latitude
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Longitude
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: SolZen
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: SatZen
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: ScatAng
   
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Emiss_11um_Tropo_Rtm
   INTEGER(BYTE), DIMENSION(:,:), POINTER, PRIVATE :: SpaceMask
   INTEGER(BYTE), DIMENSION(:,:), POINTER, PRIVATE :: CoastMask
   INTEGER(BYTE), DIMENSION(:,:), POINTER, PRIVATE :: LandMask
   INTEGER(BYTE), DIMENSION(:,:), POINTER, PRIVATE :: SnowMask
   INTEGER(BYTE), DIMENSION(:,:), POINTER, PRIVATE :: Surface_Type
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: SfcElev
   INTEGER(LONG), POINTER, PRIVATE, DIMENSION(:,:) :: Elem_Idx_LRC
   INTEGER(LONG), POINTER, PRIVATE, DIMENSION(:,:) :: Line_Idx_LRC

   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Ref_Ch5_Clear
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Ref_Ch2_Clear
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn14ClrBT
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn15ClrBT
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn7ClrRad
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn7Emiss
   INTEGER(BYTE), POINTER, PRIVATE, DIMENSION(:,:):: Bad_Pixel_Mask
   
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn11BT
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn14BT
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn15BT
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn7BT
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn9BT
   
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn7SfcEmiss

   !rchen define GlintAng
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: GlintZen
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: EmsCh7ClSlr
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Sst_Anal_Uni

!   INTEGER(LONG) :: Sfc_Idx_NWP
   REAL(SINGLE), PRIVATE :: Chn7_Sol_Energy

   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn1Refl
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn2Refl
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn3Refl
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn4Refl
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn5Refl
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn6Refl
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: Chn7Refl
   
   !I-Band Pointers
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: ChnI1Refl
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: ChnI4BT
   REAL(SINGLE), DIMENSION(:,:), POINTER, PRIVATE :: ChnI5BT
      
   !Outputs
   integer(BYTE), dimension(:,:), POINTER, PRIVATE  :: Cld_Mask
   real (SINGLE),dimension(:,:), POINTER, PRIVATE  :: Posterior_Cld_Probability
   integer(BYTE),dimension(:,:,:), POINTER, PRIVATE  :: Cld_Test_Vector_Packed
   integer(BYTE), dimension(:,:), POINTER, PRIVATE  :: Dust_Mask
   integer(BYTE), dimension(:,:), POINTER, PRIVATE  :: Smoke_Mask
   integer(BYTE), dimension(:,:), POINTER, PRIVATE  :: Fire_Mask
      
   
contains
!----------------------------------------------------------------------
! Bridge Routine
!
! Note, the Diag argument is optional
!---------------------------------------------------------------------- 
 subroutine NB_CLOUD_MASK_BRIDGE( Ctxt, Return_Status)
 
   implicit none

   TYPE(NPP_BAYES_CLOUD_MASK_Ctxt) :: Ctxt
   INTEGER(LONG) :: Return_Status

   integer :: Line_Idx, Elem_Idx
   integer:: Num_Elem
   integer:: Num_Line, Num_Scans_Read
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
   REAL(SINGLE) :: Glint_Zen_Thresh=40.0
   character (len=555):: Naive_Bayes_File_Name_Full_Path




   !---- set paths and mask classifier file name to their values in this framework
   Naive_Bayes_File_Name_Full_Path = Ctxt%CLOUD_MASK_Src1_T00%AncilPath
      
   Num_Elem = Ctxt%SegmentInfo%Current_Column_Size
   Num_Line = Ctxt%SegmentInfo%Current_Row_Size
   
   !allocate local arrays
      
   allocate(Diag_Pix_Array_1(num_elem,num_line))
   allocate(Diag_Pix_Array_2(num_elem,num_line))
   allocate(Diag_Pix_Array_3(num_elem,num_line))
   allocate(Emiss_11um_Tropo_Rtm(num_elem,num_line))
   allocate(Ems_Ch20_Median_3x3(num_elem,num_line))
   allocate(Ems_Ch20_Std_Median_3x3(num_elem,num_line))


   allocate(Dummy_IBand(num_elem,num_line))
   allocate(Ref_Uni_ChI1(num_elem,num_line))
   allocate(Bt_Uni_ChI4(num_elem,num_line))
   allocate(Bt_Uni_ChI5(num_elem,num_line))

   
   !Solar Contamination
   allocate(Solar_Contamination_Mask(num_elem, num_line))
   !Only needed for AVHRR Counts
   Solar_Contamination_Mask = sym%NO


   !Initialize pointers

   !Angles and geolocation
   CALL NFIA_Sat_ChnMap_CldMskChnUseFlg(Ctxt%SATELLITE_DATA_Src1_T00, CHN_FLG)
   CALL NFIA_Sat_Nav_ScatAng(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, ScatAng)
   CALL NFIA_Sat_Nav_SolZen(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, SolZen)
   CALL NFIA_Sat_Nav_Lat(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Latitude)
   CALL NFIA_Sat_Nav_Lon(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, Longitude)
   CALL NFIA_Sat_Nav_GlintAng(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, GlintZen)
   CALL NFIA_Sat_Nav_SatZen(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, SatZen)

   !Bad Pixel Mask
   CALL NFIA_Sat_L1b_BadPixMsk(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14, Bad_Pixel_Mask) !Bad_Pixel_Mask

   !Land/Coast/Surface Masks Surface Type
   CALL NFIA_Sat_Nav_SpaceMsk(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, SpaceMask)
   CALL NFIA_CoastMask_Mask(Ctxt%COAST_MASK_Src1_T00, CoastMask)
   CALL NFIA_LandMask_Mask(Ctxt%LAND_MASK_Src1_T00, LandMask)
   CALL NFIA_SnowMask_Mask(Ctxt%SNOW_MASK_Src1_T00, SnowMask)
   CALL NFIA_SfcType_Mask(Ctxt%SURFACE_TYPE_Src1_T00, Surface_Type)
   CALL NFIA_SfcElev_Elevation(Ctxt%SURFACE_ELEVATION_Src1_T00, SfcElev)

   
   !Brightness Temps
   CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI7, Chn7BT)   !3.9
   CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI9, Chn9BT)   !6.7
   CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI11, Chn11BT) !8.5
   CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14, Chn14BT) !11
   CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI15, Chn15BT) !12

   !reflectance percents
   CALL NFIA_Sat_L1b_ReflPrct(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI1, Chn1Refl) !0.47um
   CALL NFIA_Sat_L1b_ReflPrct(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI2, Chn2Refl) !0.64um
   CALL NFIA_Sat_L1b_ReflPrct(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI3, Chn3Refl) !0.86um
   CALL NFIA_Sat_L1b_ReflPrct(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI4, Chn4Refl) !1.3 um
   CALL NFIA_Sat_L1b_ReflPrct(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI5, Chn5Refl) !1.6um
   CALL NFIA_Sat_L1b_ReflPrct(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI6, Chn6Refl) !1.6um
   CALL NFIA_Sat_L1b_ReflPrct(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI7, Chn7Refl) !3.7um
   
   !RTM Clear sky BTs/Radiances
   CALL NFIA_RTM_Pixel_BtClr(Ctxt%RTM_Src1_T00, CHN_ABI14, Chn14ClrBT)
   CALL NFIA_RTM_Pixel_BtClr(Ctxt%RTM_Src1_T00, CHN_ABI15, Chn15ClrBT)
   CALL NFIA_RTM_Pixel_EmsCh7ClSlr(Ctxt%RTM_Src1_T00, EmsCh7ClSlr) 

   !Channel 7 emissivity and SEEBOR Emissivity
   CALL NFIA_PseudoEmiss_Chn7(Ctxt%PSEUDO_EMISSIVITY_Src1_T00, Chn7Emiss)
   CALL NFIA_SfcEmis_SfcEmiss(Ctxt%SURFACE_EMISSIVITY_Src1_T00, CHN_ABI7, Chn7SfcEmiss)
   CALL NFIP_Sat_Nav_SolEnergy(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI7, Chn7_Sol_Energy)
   
   !SST uniformity
   CALL NFIA_SST_SSTClimUnif(Ctxt%SST_Src1_T00, Sst_Anal_Uni)

   !MODIS White sky Albedo
   CALL NFIA_SfcAlbd_SfcAlbedo(Ctxt%SURFACE_ALBEDO_Src1_T00, 1, Ref_Ch2_Clear)
   CALL NFIA_SfcAlbd_SfcAlbedo(Ctxt%SURFACE_ALBEDO_Src1_T00, 3, Ref_Ch5_Clear)
   
   !need to calculate  Input%Ref_375um_Clear, ch(20)%Ref_Toa_Clear(i,j) 
   
   
   
   !Calculate 11um topopause emissivity
   CALL Compute_Emiss_Tropo_Chn14(Emiss_11um_Tropo_Rtm, Ctxt)


   !Need Cloud Mask outputs (FIXME)
   CALL NFIA_CloudMask_Mask(Ctxt%CLOUD_MASK_Src1_T00, Cld_Mask)
   CALL NFIA_CloudMask_CldProbability(Ctxt%CLOUD_MASK_Src1_T00, Posterior_Cld_Probability)
   CALL NFIA_CloudMask_CldMaskPacked(Ctxt%CLOUD_MASK_Src1_T00, Cld_Test_Vector_Packed)
   CALL NFIA_CloudMask_Dust_Mask(Ctxt%CLOUD_MASK_Src1_T00, Dust_Mask)
   CALL NFIA_CloudMask_Smoke_Mask(Ctxt%CLOUD_MASK_Src1_T00, Smoke_Mask)
   CALL NFIA_CloudMask_Fire_Mask(Ctxt%CLOUD_MASK_Src1_T00, Fire_Mask)
   
   
   !Compute Spatial uniformity

   
   CALL Compute_Spatial_Uniformity(1, 1, SpaceMask, Chn14BT, Bt_Ch31_Mean_3x3, &
                                   Bt_Ch31_Max_3x3, Bt_Ch31_Min_3x3, &
                                    Bt_Ch31_Stddev_3x3)
                                   
   CALL Compute_Spatial_Uniformity(1, 1, SpaceMask, Chn2Refl, Ref_Ch1_Mean_3X3, &
                                   Ref_Ch1_Max_3x3, Ref_Ch1_Min_3X3, &
                                   Ref_Ch1_Stddev_3X3)


   CALL Compute_Spatial_Uniformity(1, 1, SpaceMask, Chn7BT, Bt_Ch31_Mean_3x3, &
                                   Bt_Ch20_Max_3x3, Bt_Ch20_Min_3x3, &
                                   Bt_Ch20_Stddev_3x3)


   !--- apply median filter to 3.75um Emissivity
   !rchen replace Num_Scans_Read with Num_line
   call COMPUTE_MEDIAN_SEGMENT(Chn7Emiss, Bad_Pixel_Mask, 1, &
                              1, Num_Elem, 1, &
                              Num_line, &
                              Ems_Ch20_Median_3x3,  &
                              Ems_Ch20_Std_Median_3x3)
   

   !--- Calculate I-Band Uniformity
   IBand_Flag(:) = sym%NO
   
   ! McIDAS sensor ID
   CALL NFIP_Sat_Sat_ID(Ctxt%SATELLITE_DATA_Src1_T00, McIDAS_ID)
   
   !Use sensor ID to determin if I and DNB are available
   
   if (McIDAS_ID == NPP_VIIRS_MCIDAS_ID) then

!       VIIRS_375M_resolution_index = Ctxt%SegmentInfo%Res_Index(VIIRS_375M)

       IBand_Flag(1) = sym%YES
       IBand_Flag(4) = sym%YES
       IBand_Flag(5) = sym%YES

!       CALL NFIA_Sat_L1b_ReflPrct(Ctxt%SATELLITE_DATA_Src1_T00, &
!                                  VIIRS_375M_res_indx, CHN_VIIRS1, ChnI1Refl)
!       CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, &
!                                   VIIRS_375M_res_indx, CHN_VIIRS4, ChnI4BT)
!       CALL NFIA_Sat_L1b_BrtTemp(Ctxt%SATELLITE_DATA_Src1_T00, &
!                                   VIIRS_375M_res_indx, CHN_VIIRS5, ChnI5BT)

!       COMPUTE_IBAND_STATISTICS ( ChnI1Refl , Dummy_IBand, Dummy_IBand , &
!                                  Dummy_IBand, Ref_Uni_ChI1 ) 

!       COMPUTE_IBAND_STATISTICS ( ChnI4BT , Dummy_IBand, Dummy_IBand , &
!                                  Dummy_IBand, Bt_Uni_ChI4 ) 

!       COMPUTE_IBAND_STATISTICS ( ChnI5BT , Dummy_IBand, Dummy_IBand , &
!                                  Dummy_IBand, Bt_Uni_ChI5 ) 
   
   ENDIF
   IBand_Flag(:) = sym%NO

   !--- DNB reflectance 
   DNB_Flag = sym%NO



   !--- set structure (symbol, input, output, diag)  elements to corresponding values in this framework
   call SET_SYMBOL()

   ! -----------    loop over pixels -----   
   line_loop: do Line_Idx = 1, Num_Line
      elem_loop: do  Elem_Idx = 1, Num_Elem
      
      
        !-------------------------------------------------------------------
        ! Do glint mask here
        !-------------------------------------------------------------------
      
        !--- initialize valid pixel to no
        Glint_Mask = sym%NO


        !--- skip land pixels
        if ((LandMask(Elem_Idx,Line_Idx) == sym%NO) .and. &
          SnowMask(Elem_Idx,Line_Idx) == sym%NO_SNOW) then

       !--- turn on in geometric glint cone and sufficient Ref_Ch1
            if ((Glintzen(Elem_Idx,Line_Idx) < Glint_Zen_Thresh)) then

            !--- assume to be glint if in geometric zone
                Glint_Mask = sym%YES

                IF (CHN_FLG(14) == sym%YES) then

                !--- exclude pixels colder than the freezing temperature
                    IF (Chn14BT(Elem_Idx,Line_Idx) < 273.15) then
                        Glint_Mask = sym%NO
                    endif

            !--- exclude pixels colder than the surface
                    IF (Chn14BT(Elem_Idx,Line_Idx) < Chn14ClrBT(Elem_Idx,Line_Idx) - 5.0) then
                        Glint_Mask = sym%NO
                    endif

                endif

          !-turn off if non-uniform in reflectance
                IF (CHN_FLG(2) == sym%YES) then
                    !rchen
                    IF (Ref_Ch1_Stddev_3x3(Elem_Idx,Line_Idx) > 1.0) then
                        Glint_Mask = sym%NO
                    endif
                endif

            endif  !Glintzen check

        endif    !land check
      
      
        !-------------------------------------------------------------------
        ! Do covariance here
        !-------------------------------------------------------------------
        Elem_Idx_min = max(1,min(Elem_Idx - 2,Num_Elem))
        Elem_Idx_max = max(1,min(Elem_Idx + 2,Num_Elem))
        Line_Idx_min = max(1,min(Line_Idx - 2,Num_line))
        Line_Idx_max = max(1,min(Line_Idx + 2,Num_line))
        Line_Idx_width = Line_Idx_max - Line_Idx_min + 1
        Elem_Idx_width = Elem_Idx_max - Elem_Idx_min + 1
        
        Covar_Ch27_Ch31_5x5 = Missing_Value_Real4
      
         IF ((CHN_FLG(9) == sym%YES) .and. & 
            ( CHN_FLG(14)== sym%YES)) THEN

            Covar_Ch27_Ch31_5x5 = Covariance_local(&
              Chn14BT(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
              Chn9BT(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
              Elem_Idx_width, Line_Idx_width, &
               Bad_Pixel_Mask(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max))
        ENDIF
     
      
         ! Set inputs
         
         
         call SET_INPUT(Elem_Idx,Line_Idx,Ctxt)
         call SET_OUTPUT(Elem_Idx,Line_Idx)
         call SET_DIAG(Elem_Idx,Line_Idx)

         !---call cloud mask routine
         call NB_CLOUD_MASK_ALGORITHM( &
                      Naive_Bayes_File_Name_Full_Path, &
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
   
      end do elem_loop
   end do line_loop

!-------------------------------------------------------------------------------
! on last segment, wipe out the lut from memory and reset is_read_flag to no
!-------------------------------------------------------------------------------
   if (Segment_Number_CM == Input%Num_Segments) then
       call RESET_NB_CLOUD_MASK_LUT()
   endif


   
   
   !Deallocate arrays
   deallocate(Diag_Pix_Array_1)
   deallocate(Diag_Pix_Array_2)
   deallocate(Diag_Pix_Array_3)
   deallocate(Emiss_11um_Tropo_Rtm)
   deallocate(Ems_Ch20_Median_3x3)
   deallocate(Ems_Ch20_Std_Median_3x3)
   deallocate(Solar_Contamination_Mask)
   
   if (allocated(Ref_Ch1_Mean_3X3)) deallocate(Ref_Ch1_Mean_3X3)
   if (allocated(Ref_Ch1_Max_3x3)) deallocate(Ref_Ch1_Max_3x3)
   if (allocated(Ref_Ch1_Min_3X3)) deallocate(Ref_Ch1_Min_3X3)
   if (allocated(Ref_Ch1_Stddev_3X3)) deallocate(Ref_Ch1_Stddev_3X3)

   if (allocated(Bt_Ch31_Mean_3x3)) deallocate(Bt_Ch31_Mean_3x3)
   if (allocated(Bt_Ch31_Max_3x3)) deallocate(Bt_Ch31_Max_3x3)
   if (allocated(Bt_Ch31_Min_3x3)) deallocate(Bt_Ch31_Min_3x3)
   if (allocated(Bt_Ch31_Stddev_3x3)) deallocate(Bt_Ch31_Stddev_3x3)

   if (allocated(Bt_Ch20_Mean_3x3)) deallocate(Bt_Ch20_Mean_3x3)
   if (allocated(Bt_Ch20_Max_3x3)) deallocate(Bt_Ch20_Max_3x3)
   if (allocated(Bt_Ch20_Min_3x3)) deallocate(Bt_Ch20_Min_3x3)
   if (allocated(Bt_Ch20_Stddev_3x3)) deallocate(Bt_Ch20_Stddev_3x3)

      
   !Pointers
   CHN_FLG => null()
   Latitude => null()
   Longitude => null()
   SolZen => null()
   SatZen => null()
   ScatAng => null()  
   SpaceMask => null()
   CoastMask => null()
   LandMask => null()
   SnowMask => null()
   Surface_Type => null()
   SfcElev => null()
   Ref_Ch5_Clear => null()
   Ref_Ch2_Clear => null()
   Chn14ClrBT => null()
   Chn15ClrBT => null()
   Chn7ClrRad => null()
   Chn7Emiss => null()
   Bad_Pixel_Mask  => null()
   Chn11BT => null()
   Chn14BT => null()
   Chn15BT => null()
   Chn7BT => null()
   Chn9BT => null()
   Chn7SfcEmiss => null()
   GlintZen => null()
   EmsCh7ClSlr => null()
   Sst_Anal_Uni => null()
   Chn1Refl => null()
   Chn2Refl => null()
   Chn3Refl => null()
   Chn4Refl => null()
   Chn5Refl => null()
   Chn6Refl => null()
   Chn7Refl => null()
   Cld_Mask => null()
   Posterior_Cld_Probability => null()
   Cld_Test_Vector_Packed => null()
   Dust_Mask => null()
   Smoke_Mask => null()
   Fire_Mask => null()
   
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
   subroutine SET_INPUT(i,j, Ctxt)
      integer, intent (in) :: i, j
      TYPE(NPP_BAYES_CLOUD_MASK_Ctxt) :: Ctxt

      Input%Num_Elem = Ctxt%SegmentInfo%Current_Column_Size
      Input%Num_Line = Ctxt%SegmentInfo%Current_Row_Size
      Input%Num_Line_Max = Ctxt%SegmentInfo%Current_Row_Size
      Input%Num_Segments = Ctxt%SegmentInfo%Segment_total
      !------
      Input%Invalid_Data_Mask => Bad_Pixel_Mask(i,j)
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
      Input%Snow_Class => SnowMask(i,j)
      Input%Land_Class => LandMask(i,j)
      Input%Oceanic_Glint_Mask => Glint_Mask
      Input%Lunar_Oceanic_Glint_Mask => null() !No DNB in Framework now
      Input%Coastal_Mask => CoastMask(i,j)
      Input%Solzen => SolZen(i,j)
      Input%Scatzen => ScatAng(i,j)
      Input%Lunscatzen => null() !No DNB in Framework now
      Input%Senzen => Satzen(i,j)
      Input%Lunzen => null() !No DNB in Framework now
      Input%Lat => Latitude(i,j)
      Input%Lon => Longitude(i,j)
      Input%Ref_041um => Chn1Refl(i,j)
      Input%Ref_063um => Chn2Refl(i,j)
      Input%Ref_063um_Clear => Ref_Ch2_Clear(i,j)
      Input%Ref_063um_Std => Ref_Ch1_Stddev_3X3(i,j)
      Input%Ref_063um_Min => Ref_Ch1_Min_3x3(i,j)
      Input%Ref_086um => Chn3Refl(i,j)
      Input%Ref_138um => Chn4Refl(i,j)
      Input%Ref_160um => Chn5Refl(i,j)
      Input%Ref_160um_Clear => Ref_Ch5_Clear(i,j)
      Input%Ref_375um => Chn7Refl(i,j)
      Input%Ref_375um_Clear => null() !Not filled or used for now
      Input%Ref_213um => Chn6Refl(i,j)
      Input%Bt_375um => Chn7BT(i,j)
      Input%Bt_375um_Std => Bt_Ch20_Stddev_3x3(i,j)
      Input%Emiss_375um =>  Ems_Ch20_Median_3x3(i,j)
      Input%Emiss_375um_Clear => EmsCh7ClSlr(i,j)
      Input%Bt_67um => Chn9BT(i,j)
      Input%Bt_85um => Chn11BT(i,j)
      Input%Bt_11um => Chn14BT(i,j)
      Input%Bt_11um_Std => Bt_Ch31_Stddev_3x3(i,j)
      Input%Bt_11um_Max => Bt_Ch31_Max_3x3(i,j)
      Input%Bt_11um_Clear => Chn14ClrBT(i,j)
      Input%Emiss_11um_Tropo => Emiss_11um_Tropo_Rtm(i,j)
      Input%Bt_12um => Chn15BT(i,j)
      Input%Bt_12um_Clear => Chn15ClrBT(i,j)
      Input%Ref_I1_064um_Std => null()
      Input%Bt_I4_374um_Std => null()
      Input%Bt_I5_114um_Std => null()
      Input%Bt_11um_Bt_67um_Covar => Covar_Ch27_Ch31_5x5
      Input%Sst_Anal_Uni => Sst_Anal_Uni(i,j)
      Input%Emiss_Sfc_375um => Chn7SfcEmiss(i,j)
      Input%Rad_Lunar => null() !No DNB
      Input%Ref_Lunar => null() !No DNB
      Input%Ref_Lunar_Min => null() !No DNB
      Input%Ref_Lunar_Std => null() !No DNB
      Input%Ref_Lunar_Clear => null() !No DNB
      Input%Zsfc => SfcElev(i,j)
      Input%Solar_Contamination_Mask => Solar_Contamination_Mask(i,j)
      Input%Sfc_Type => Surface_Type(i,j)
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


!rchen copy Clear_Chn2_Reflectance_Field from baseline_cloud_mask.f90
!====================================================================
! Subroutine Name: Clear_Chn2_Reflectance_Field
!
! Function:
!    Routine to Read in Chn2 Reflectance Fields
!
! Description:
!   This subroutine reads in the channel 2 reflectance fields for a given
!     segment
!
! Calling Sequence:
!   CALL Clear_Chn2_Reflectance_Field(Refl_Chn2_Clear, Ctxt)
!
! Inputs:
!  Ctxt : Context pointer providing access to framework segment information
!
! Outputs:
!  Refl_Chn2_Clear - Clear sky channel 2 albedo
!
! Dependencies: None
!
! Restrictions:  None
!
!====================================================================
 SUBROUTINE Clear_Chn2_Reflectance_Field(Refl_Chn2_Clear, Ctxt)

   REAL(SINGLE), INTENT(INOUT), ALLOCATABLE, DIMENSION(:,:) :: Refl_Chn2_Clear
   INTEGER(LONG) :: Num_Elem
   INTEGER(LONG) :: Max_Num_Lines_per_Seg
   !tyu
   !TYPE(FW_Context), POINTER :: Ctxt
   TYPE(NPP_BAYES_CLOUD_MASK_Ctxt) :: Ctxt
   !tyu
   !=== INFO: generated declarations
   INTEGER(BYTE), DIMENSION(:,:), POINTER :: LandMask
   REAL(SINGLE), DIMENSION(:,:), POINTER :: WhiteSkyAlbChn2
   !=== INFO: generated initializations
   !REPLACE    CALL FWI_LandMask(Out=LandMask, Ctxt=Ctxt)
   CALL NFIA_LandMask_Mask(Ctxt%LAND_MASK_Src1_T00, LandMask)

   !REPLACE    CALL FWI_WhiteSkyAlb(Chn=CHN_ABI2, Out=WhiteSkyAlbChn2, Ctxt=Ctxt)
   !CALL NFIA_SfcAlbd_SfcAlbedo(Ctxt%SURFACE_ALBEDO_Src1_T00, CHN_ABI2, WhiteSkyAlbChn2)
   !tyu use real channel number : not the logical
   ! CHN_ABI2==>1 CHN_ABI3==>2 CHN_ABI5==>3
   CALL NFIA_SfcAlbd_SfcAlbedo(Ctxt%SURFACE_ALBEDO_Src1_T00, 1, WhiteSkyAlbChn2)

   !print *, "WhiteSkyAlbChn2"
   !print *,  WhiteSkyAlbChn2

   !=== INFO: end of generated code
   !REPLACE    CALL FWI_Seg( ElemCount=Num_Elem, LineCount=Max_Num_Lines_per_Seg, Ctxt=Ctxt )
   Num_Elem = Ctxt%SegmentInfo%Current_Column_Size
   !--RDR Use Allocate_Rows instead, Size diff in LandMask and Refl_Chn2_Clear causes problems
   ! with GNU on last segment
   !Max_Num_Lines_per_Seg = Ctxt%SegmentInfo%Current_Row_Size
   Max_Num_Lines_per_Seg = Ctxt%SegmentInfo%Allocate_Rows

   !--- If create one and populate background reflectance using default or global WS albedo data
       IF (.NOT. ALLOCATED(Refl_Chn2_Clear)) ALLOCATE( &
                           Refl_Chn2_Clear(Num_Elem,Max_Num_Lines_per_Seg))

       Refl_Chn2_Clear = 5.0
       WHERE(LandMask == TOPO_LAND)
          Refl_Chn2_Clear = 45.0
       ENDWHERE

       !--- use MODIS white sky from mapped data if available
       WHERE(WhiteSkyAlbChn2 /= Missing_Value_Real4)
          Refl_Chn2_Clear = WhiteSkyAlbChn2
       ENDWHERE

       !--- ensure consistent missing value
       WHERE(Refl_Chn2_Clear <= 0.0)
          Refl_Chn2_Clear = Missing_Value_Real4
       ENDWHERE

 END SUBROUTINE Clear_Chn2_Reflectance_Field

!-----------------------
! updated on 09/08/2010
!-----------------------
 SUBROUTINE Compute_Emiss_Tropo_Chn14(Emiss_Tropo_Chn14,Ctxt)
   INTEGER(LONG) :: Number_of_Lines_in_this_Segment
   REAL(SINGLE), DIMENSION(:,:), INTENT(OUT) :: Emiss_Tropo_Chn14
   !TYPE(FW_Context), POINTER :: Ctxt
   TYPE(NPP_BAYES_CLOUD_MASK_Ctxt) :: Ctxt
   INTEGER(LONG) :: Tropo_Idx_NWP
   !INTEGER(BYTE) :: View_Zen_Idx
   INTEGER :: X_NWP_Idx
   INTEGER :: Y_NWP_Idx
   INTEGER :: Elem_Idx
   INTEGER :: Line_Idx
   REAL(SINGLE) :: Rad_Chn14
   REAL(SINGLE) :: Clr_Rad_Chn14
   REAL(SINGLE) :: Blkbdy_Tropo_Rad_Chn14
   !=== INFO: generated declarations
   REAL(SINGLE), DIMENSION(:), POINTER :: BlackCldRadProfChn14
   REAL(SINGLE), DIMENSION(:,:), POINTER :: Chn14ClrRad
   INTEGER(LONG) :: Elems
   REAL(SINGLE), DIMENSION(:,:), POINTER :: Chn14Rad
   INTEGER(BYTE), DIMENSION(:,:), POINTER :: SpaceMask
   !=== INFO: generated initializations
   CALL NFIA_RTM_Pixel_RadClr(Ctxt%RTM_Src1_T00, CHN_ABI14, Chn14ClrRad)

   CALL NFIA_Sat_L1b_Rad(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, CHN_ABI14, Chn14Rad)

   CALL NFIA_Sat_Nav_SpaceMsk(Ctxt%SATELLITE_DATA_Src1_T00, COMMON_RESOLUTION, SpaceMask)

   Elems = Ctxt%SegmentInfo%Current_Column_Size
   Number_of_Lines_in_this_Segment = Ctxt%SegmentInfo%Current_Row_Size

   !=== INFO: end of generated code

   !--- initialize
   Emiss_Tropo_Chn14 = Missing_Value_Real4

    Line_Loop: DO Line_Idx=1, Number_of_Lines_in_this_Segment
      Element_Loop: DO Elem_Idx = 1, Elems

       IF (SpaceMask(Elem_Idx,Line_Idx) == NO) THEN

            !
            !---nwp level associated with tropopause
            CALL NFIP_NWP_TropoLevel(Ctxt%NWP_DATA_Src1_T00, Elem_Idx, Line_Idx, Tropo_Idx_NWP)

            !
            !---viewing zenith angle bin
            !
            ! View_Zen_Idx =          sat%ivza(Elem_Idx,Line_Idx)

            !
            !---11 um radiance
            !
            Rad_Chn14  =        Chn14Rad(Elem_Idx,Line_Idx)

            !
            !---clear 11 micron radiance
            !
            Clr_Rad_Chn14 =     Chn14ClrRad(Elem_Idx,Line_Idx)

            !
            !---BB 11 um rad at tropopause
            !

            CALL NFIA_RTM_Grid_CloudProf(Ctxt%RTM_Src1_T00, Elem_Idx, Line_Idx, CHN_ABI14, BlackCldRadProfChn14)

            Blkbdy_Tropo_Rad_Chn14 = BlackCldRadProfChn14(Tropo_Idx_NWP)
            !
            !---Tropopause Emissivity
            !
            Emiss_Tropo_Chn14(Elem_Idx,Line_Idx) =  &
                  (Rad_Chn14 - Clr_Rad_Chn14) / (Blkbdy_Tropo_Rad_Chn14 - Clr_Rad_Chn14)

      END IF
    END DO Element_Loop
  END DO Line_Loop

 END SUBROUTINE Compute_Emiss_Tropo_Chn14


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

   !----------------------------------------------------------------
   ! - iband has full file dimension of 6400 x1536
   ! - mband 3200 x 768
   !  - output of min_val ... is 3200 768
   !----------------------------------------------------------------
   subroutine COMPUTE_IBAND_STATISTICS ( iband_array , out_min_val, out_max_val , out_mean_val, out_std_val )
      implicit none
      real, dimension(:,:) , intent(in) :: iband_array
      real, dimension(:,:) , intent(out)  :: out_min_val, out_max_val , out_mean_val , out_std_val
      real, dimension(2,2) :: small_iband
      
      integer :: im , jm
      integer , dimension(2) ::  dim_m
      integer :: iband_x0, iband_x1 ,  iband_y0, iband_y1
     
      
      dim_m = shape ( out_min_val )
      
      out_min_val =-999.
      out_max_val =-999.
      out_mean_val =-999.
      
      do im = 1   , dim_m( 1)
        
         iband_x0 = ( im-1) * 2 +1
         iband_x1 = iband_x0 + 1
         do jm =1 , dim_m( 2 )
            iband_y0 = ( jm-1) * 2 +1
            iband_y1 = iband_y0 + 1
            small_iband = iband_array ( iband_x0 :  iband_x1 ,  iband_y0 :  iband_y1 )
            if ( minval ( small_iband ) > 0 ) then 
               out_min_val ( im, jm ) = minval ( small_iband )
               out_max_val ( im, jm ) = maxval ( small_iband )
               out_mean_val ( im, jm ) = sum ( small_iband ) / 4.0
               out_std_val ( im, jm ) = SQRT ( ( SUM ( (small_iband - out_mean_val(im, jm )) **2 ) ) / ( 4. - 1. ) ) 
            end if
         end do
      end do
   
   end subroutine COMPUTE_IBAND_STATISTICS


end module NB_CLOUD_MASK_SAPF_BRIDGE
