!$Id$
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
!
! Note, to use the diagnostic variables, do this
!   - set the Use_Diag flag to true
!   - turn on the Diag argument to the desirefd routine
!   - in the desired routine, set the diag variables to what you want
!   - when done, repeat this in reverse
!
!-------------------------------------------------------------------------------
module NB_CLOUD_MASK_GSIP_BRIDGE


   use NB_CLOUD_MASK_SERVICES
   use NB_CLOUD_MASK

   implicit none

   public :: NB_CLOUD_MASK_BRIDGE

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
   
   !Make module wide variables
   character (len=1020), TARGET, PRIVATE:: Ancil_Data_Path
   character (len=1020), TARGET, PRIVATE:: Naive_Bayes_File_Name

   !Segment counter
   integer(kind=INT1), TARGET, PRIVATE:: Segment_Number_CM = 1

   !Make Iband and DNB flag
   integer (kind=INT1), DIMENSION(5), PRIVATE :: IBand_Flag
   integer (kind=INT1), PRIVATE :: DNB_Flag
   
   !allocatable
   integer (kind=INT1), DIMENSION(:,:), ALLOCATABLE, TARGET, PRIVATE :: Solar_Contamination_Mask

   integer (kind=INT1), TARGET, PRIVATE :: Glint_Mask
      
   !Pointers
   INTEGER(kind=INT1), DIMENSION(:), POINTER, PRIVATE :: CHN_FLG

!   INTEGER(LONG) :: Sfc_Idx_NWP
   REAL(kind=REAL4), PRIVATE :: Chn7_Sol_Energy

   
contains
!----------------------------------------------------------------------
! Bridge Routine
!
! Note, the Diag argument is optional
!---------------------------------------------------------------------- 
 subroutine NB_CLOUD_MASK_BRIDGE( )
 
   implicit none

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
   character (len=1020):: Naive_Bayes_File_Name_Full_Path
   logical:: Use_Diag

   Use_Diag = .false.

   !---- set paths and mask classifier file name to their values in this framework
   Naive_Bayes_File_Name_Full_Path = trim(ancil_path)//"bayes/"//trim(Bayesian_Cloud_Mask_Name)
         
   Num_Elem = Num_Pix
   Num_Line = Num_Scans_Read
   Num_Line_Max = Num_Scans_Per_Segment

   !allocate local arrays
   !Solar Contamination
   allocate(Solar_Contamination_Mask(num_elem, num_line))
   !Only needed for AVHRR Counts
   Solar_Contamination_Mask = sym%NO
      
   !Point to channel flag
   CHN_FLG => sat_info_gsip(1)%chanon

   !--- I-Band Flag
   IBand_Flag(:) = sym%NO

   !--- DNB reflectance 
   DNB_Flag = sym%NO

   !--- set structure (symbol, input, output, diag)  elements to corresponding values in this framework
   call SET_SYMBOL()

   ! -----------    loop over pixels -----   
   line_loop: do Line_Idx = 1, Num_Line
      elem_loop: do  Elem_Idx = 1, Num_Elem
      
        !-------------------------------------------------------------------
        ! Do space mask check here
        !-------------------------------------------------------------------

        IF (Space_Mask(Elem_Idx,Line_Idx) == sym%YES) THEN
            gsip_pix_prod%cldmask(Elem_Idx,Line_Idx) = missing_value_int1
            CYCLE
        ENDIF
       
       
        !-------------------------------------------------------------------
        ! Do glint mask here
        !-------------------------------------------------------------------
      
        !--- initialize valid pixel to no
        Glint_Mask = sym%NO

        !--- skip land pixels
        IF ((Land_Mask(Elem_Idx,Line_Idx) == sym%NO) .and. &
          Snow_mask(Elem_Idx,Line_Idx) == sym%NO_SNOW) THEN

       !--- turn on in geometric glint cone and sufficient Ref_Ch1
            IF ((glintzen(Elem_Idx,Line_Idx) < Glint_Zen_Thresh)) THEN

            !--- assume to be glint IF in geometric zone
                Glint_Mask = sym%YES

                IF (CHN_FLG(14) == sym%YES) THEN

                !--- exclude pixels colder than the freezing temperature
                    IF (bt14(Elem_Idx,Line_Idx) < 273.15) THEN
                        Glint_Mask = sym%NO
                    ENDIF

            !--- exclude pixels colder than the surface
                    IF (bt14(Elem_Idx,Line_Idx) < Bt_Clear_Ch14_Rtm(Elem_Idx,Line_Idx) - 5.0) THEN
                        Glint_Mask = sym%NO
                    ENDIF

                ENDIF

          !-turn off IF non-unIForm in reflectance
                IF (CHN_FLG(2) == sym%YES) THEN
                    IF (Ref_Ch1_Std_3x3(Elem_Idx,Line_Idx) > 1.0) THEN
                        Glint_Mask = sym%NO
                    ENDIF
                ENDIF

            ENDIF  !Glintzen check

        ENDIF    !land check

         ! Set inputs
         call SET_INPUT(Elem_Idx,Line_Idx)

         !---call cloud mask routine
         call NB_CLOUD_MASK_ALGORITHM( &
                      Naive_Bayes_File_Name_Full_Path, &
                      Symbol,  &
                      Input, &
                      Output)

         call SET_OUTPUT(Elem_Idx,Line_Idx)
         call SET_DIAG(Elem_Idx,Line_Idx)

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
   deallocate(Solar_Contamination_Mask)
   
      
   !Pointers
   CHN_FLG => null()
 
   
   !Increment segment number
   Segment_Number_CM = Segment_Number_CM +1

   end subroutine NB_CLOUD_MASK_BRIDGE


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

      Input%Num_Elem = nx
      Input%Num_Line = num_scans_read
      Input%Num_Line_Max = Num_Scans_Per_Segment
      Input%Num_Segments = nseg
      !------
      Input%Invalid_Data_Mask => bad_pix_mask(14,i,j)
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

      Input%Coastal_Mask => Coast_Mask(i,j)
      Input%Snow_Class => Snow_mask(i,j)
      Input%Land_Class => Land_Mask(i,j)
      Input%Sfc_Type => Sfc_Type(i,j)
      
      
      Input%Oceanic_Glint_Mask => Glint_Mask
      Input%Solzen => SolZen(i,j)
      Input%Scatzen => scatzen(i,j)
      Input%Senzen => Satzen(i,j)
      Input%Lunzen => null() !No DNB in Framework now
      Input%Lat => Lat(i,j)
      Input%Lon => Lon(i,j)
      
      IF (Input%Chan_On_063um == sym%YES) THEN
        Input%Ref_063um => ref2(i,j)
        Input%Ref_063um_Clear => sat%ws_albedo2(i,j)
        Input%Ref_063um_Std => Ref_Ch1_Std_3x3(i,j)
        Input%Ref_063um_Min => Ref_Ch1_Min_3x3(i,j)      
      ENDIF
      
      IF (Input%Chan_On_086um == sym%YES) Input%Ref_086um => ref3(i,j)
      IF (Input%Chan_On_138um == sym%YES) Input%Ref_138um => ref4(i,j)

      IF (Input%Chan_On_160um == sym%YES) THEN
        Input%Ref_160um => ref5(i,j)
        Input%Ref_160um_Clear => sat%ws_albedo5(i,j)
      ENDIF


      IF (Input%Chan_On_213um == sym%YES) Input%Ref_213um => ref6(i,j)

      IF (Input%Chan_On_375um == sym%YES) THEN
        Input%Ref_375um => ref7(i,j)
        Input%Ref_375um_Clear => null() !Not filled or used for now
        Input%Bt_375um => bt7(i,j)
        Input%Bt_375um_Std => Bt_Ch20_Std_3x3(i,j)
        Input%Emiss_375um =>  Ems_Ch20_Median_3x3(i,j)
        Input%Emiss_375um_Clear => Ems_Ch7_Clear_Solar_Rtm(i,j)
        Input%Emiss_Sfc_375um => sfc_emiss_7(i,j)
      ENDIF
      
      IF (Input%Chan_On_67um == sym%YES) THEN
        Input%Bt_67um => bt9(i,j)
        IF (Input%Chan_On_11um == sym%YES) Input%Bt_11um_Bt_67um_Covar => Covar_Ch27_Ch31_5x5(i,j)
      
      ENDIF
      
      
      IF (Input%Chan_On_85um == sym%YES) Input%Bt_85um => bt11(i,j)
      
      IF (Input%Chan_On_11um == sym%YES)  THEN
        Input%Bt_11um => bt14(i,j)
        Input%Bt_11um_Std => Bt_Ch31_Std_3x3(i,j)
        Input%Bt_11um_Max => Bt_Ch31_Max_3x3(i,j)
        Input%Bt_11um_Clear => Bt_Clear_Ch14_Rtm(i,j)
        Input%Emiss_11um_Tropo => Emiss_11um_Tropo_Rtm(i,j)
      ENDIF
      
      IF (Input%Chan_On_12um == sym%YES) THEN
        Input%Bt_12um => bt15(i,j)
        Input%Bt_12um_Clear => Bt_Clear_Ch15_Rtm(i,j)
      ENDIF
     
      Input%Sst_Anal_Uni => Sst_Anal_Uni(i,j)
      Input%Zsfc => Zsfc(i,j)
      Input%Solar_Contamination_Mask => Solar_Contamination_Mask(i,j)
   end subroutine SET_INPUT

   subroutine SET_OUTPUT(i,j)
      integer, intent (in) :: i, j

      Output%Cld_Flags_Packed => gsip_pix_prod%cldmask_packed(:,i,j)
      Output%Cld_Mask_Bayes => gsip_pix_prod%cldmask(i,j)
      Output%Posterior_Cld_Probability => gsip_pix_prod%cldprob(i,j)
      Output%Dust_Mask => null()
      Output%Smoke_Mask => null()
      Output%Fire_Mask => null()
      Output%Thin_Cirr_Mask => null()
   end subroutine SET_OUTPUT

   subroutine SET_DIAG(i,j)
      integer, intent (in) :: i, j

      Diag%Array_1 => null()
      Diag%Array_2 => null()
      Diag%Array_3 => null()
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
    Output%Thin_Cirr_Mask => null()
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


end module NB_CLOUD_MASK_GSIP_BRIDGE
