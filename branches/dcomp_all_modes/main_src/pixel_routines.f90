! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: pixel_routines.f90 (src)
!       PIXEL_ROUTINES (program)
!
! PURPOSE: this module houses routines for computing some needed pixel-level arrays
!
! DESCRIPTION: 
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
! Public routines used in this MODULE:
! COMPUTE_PIXEL_ARRAYS - compute some commonly used arrays
! COMPUTE_TSFC - derive pixel-level surface temperature 
! ATMOS_CORR - perform atmospheric correction
! NORMALIZE_REFLECTANCES - divide reflectances by cosine solar zenith angle
! CH20_PSEUDO_REFLECTANCE - compute the channel 20 reflectance
! COMPUTE_SPATIAL_UNIFORMITY - compute metrics of radiance and reflectance
!                              spatial uniformity
! SPECTRAL_CORRECT_NDVI - apply a spectral correct to Ndvi to look like NOAA14
! ASSIGN_CLEAR_SKY_QUALITY_FLAGS - assign quality flags to clear-sky products
! CONVERT_TIME - compute a time in hours based on millisecond time in leveL1b
! COMPUTE_SNOW_CLASS - based on Snow information, make a Snow Classification
! COMPUTE_GLINT - derive a glint mask
!
! DETERMINE_LEVEL1B_COMPRESSION
!
!--------------------------------------------------------------------------------------
MODULE PIXEL_ROUTINES
   
   use CX_SCIENCE_TOOLS_MOD,only: &
      Glint_Zen_Thresh &
      , Ref_Sfc_White_Sky_Water
   
   use CX_CONSTANTS_MOD
      
   use LAND_SFC_PROPERTIES,only: &
      Land_grid_description &      
      , read_land_sfc_hdf
   
   use NUMERICAL_TOOLS_MOD,only: &
    COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES &
    , Covariance

   use NWP_COMMON,only: &
      Zsfc_nwp &
      , Bad_Nwp_Mask
   
   use PIXEL_COMMON   
      
   use PLANCK,only: &
      PLANCK_TEMP_FAST &
      , PLANCK_RAD_FAST &
      , PLANCK_TEMP_FAST
   
   use SFC_PROP_UMD_MOD,only: &
      Ch1_sfc_alb_umd &
      , Ch1_Snow_Sfc_Alb_Umd &
      , Ch2_Sfc_Alb_Umd &
      , Ch2_Snow_Sfc_Alb_Umd &
      , Ch6_Sfc_Alb_Umd &
      , Ch5_Snow_Sfc_Alb_Umd &
      , Ch6_Snow_Sfc_Alb_Umd &
      , Ch7_Snow_Sfc_Alb_Umd
   
   
   

 implicit none
 private
 public:: COMPUTE_PIXEL_ARRAYS, &
          SURFACE_REMOTE_SENSING,  &
          ATMOS_CORR,&
          NORMALIZE_REFLECTANCES,  &
          CH20_PSEUDO_REFLECTANCE,  &
          COMPUTE_SPATIAL_UNIFORMITY, &
          ASSIGN_CLEAR_SKY_QUALITY_FLAGS, &
          COMPUTE_SNOW_CLASS, &
          COMPUTE_SNOW_CLASS_NWP, &
          COMPUTE_SNOW_CLASS_OISST, &
          EXPAND_SPACE_MASK_FOR_USER_LIMITS, &
          SET_SOLAR_CONTAMINATION_MASK, &
          SET_BAD_PIXEL_MASK, &
          QUALITY_CONTROL_ANCILLARY_DATA,   &
          READ_MODIS_WHITE_SKY_ALBEDO,      &
          COMPUTE_CLEAR_SKY_SCATTER,        &
          COMPUTE_GLINT,                    &
          SET_CHAN_ON_FLAG,                 &
          COMPUTE_SPATIAL_CORRELATION_ARRAYS, &
          DETERMINE_LEVEL1B_COMPRESSION, &
          TERM_REFL_NORM, &
          MERGE_NWP_HIRES_ZSFC, &
          ADJACENT_PIXEL_CLOUD_MASK, &
          COMPUTE_ACHA_PERFORMANCE_METRICS, &
          COMPUTE_DCOMP_PERFORMANCE_METRICS, &
          COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS, &
          MODIFY_LAND_CLASS_WITH_NDVI, &
          DESERT_MASK_FOR_CLOUD_DETECTION, &
          CITY_MASK_FOR_CLOUD_DETECTION, &
          VIIRS_TO_MODIS

  private:: REMOTE_SENSING_REFLECTANCE, &
            NORMALIZED_DIFFERENCE_VEGETATION_INDEX, &
            NORMALIZED_DIFFERENCE_SNOW_INDEX, &
            COMPUTE_TSFC

  contains

   !----------------------------------------------------------------------
   ! set Chan_On_Flag for each to account for Ch3a/b switching on avhrr
   !
   ! this logic allows the default values to also be used to turn off
   ! channels
   !
   !  called by process_clavrx inside segment loop
   !    HISTORY: 2014/12/29 (AW); removed unused ch3a_on_avhrr for modis and goes
   !----------------------------------------------------------------------
   subroutine SET_CHAN_ON_FLAG(Chan_On_Flag_Default, Chan_On_Flag_Per_Line)

      logical, dimension(:), intent(in):: Chan_On_Flag_Default
      logical, dimension(:,:), intent(out):: Chan_On_Flag_Per_Line
      integer:: Number_of_Elements
      integer:: Number_of_Lines
      integer:: Line_Idx
      logical :: dcomp_first_valid_line_avhrr_set = .false. 
      
      Number_of_Elements = Image%Number_Of_Elements
      Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment
      
      line_loop: do Line_Idx = 1, Number_Of_Lines
          ! - change dcomp _mode according ch3a_on flag
          ! - this change of dcomp_mode is only possible once for one file
          ! - First daytime line determines dcomp mode for whole file
          ! - AW 02/13/2017
         if ( index(Sensor%Sensor_Name,'AVHRR') > 0 &
            .and. .not. dcomp_first_valid_line_avhrr_set  &
            .and. Geo%Solzen(1,line_idx) .lt. 82) then
            if (Ch3a_On_Avhrr(Line_Idx) == sym%YES) then
                dcomp_first_valid_line_avhrr_set = .true.
                dcomp_mode = 1
            end if
            if (Ch3a_On_Avhrr(Line_Idx) == sym%NO) then
               dcomp_first_valid_line_avhrr_set = .true.
                dcomp_mode = 3
            
            end if
         end if 
           
         ! - for all sensors : set chan_on_flag ( dimension [n_chn, n_lines] to default ) 
         Chan_On_Flag_Per_Line(:,Line_Idx) = Chan_On_Flag_Default   
         
         ! two exceptions
         if (trim(Sensor%Platform_Name)=='AQUA' .and. Chan_On_Flag_Default(6) ) then
            if (minval(ch(6)%Unc(:,Line_Idx)) >= 15) then
                 Chan_On_Flag_Per_Line(6,Line_Idx) = .false.
            end if  
         end if
         
         if (index(Sensor%Sensor_Name,'AVHRR') > 0) then
            if (Ch3a_On_Avhrr(Line_Idx) == sym%YES ) then
               Chan_On_Flag_Per_Line(6,Line_Idx) = Chan_On_Flag_Default(6)   
               Chan_On_Flag_Per_Line(20,Line_Idx) = .false.   
            end if
            if (Ch3a_On_Avhrr(Line_Idx) == sym%NO) then
               Chan_On_Flag_Per_Line(6,Line_Idx) = .false.  
               Chan_On_Flag_Per_Line(20,Line_Idx) = Chan_On_Flag_Default(20)   
            end if
         endif

      end do line_loop
   
   end subroutine SET_CHAN_ON_FLAG



!======================================================================
! Modify the space mask based the limits on lat, lon, satzen and solzen 
! Space mask was initially determined in the Level-1b Navigation
!======================================================================
subroutine EXPAND_SPACE_MASK_FOR_USER_LIMITS(Space_Mask)

   integer(kind=int1), dimension(:,:), intent(inout):: Space_Mask

   where(Nav%Lat_1b /= Missing_Value_Real4 .and. Nav%Lon_1b /= Missing_Value_Real4)
        Space_Mask = sym%NO
   elsewhere 
        Space_Mask = sym%YES
   end where

   where(isnan(Nav%Lat_1b) .or. isnan(Nav%Lon_1b))
        Space_Mask = sym%YES
   end where

   where(Nav%Lat_1b < Nav%Lat_Min_Limit .or. Nav%Lat_1b > Nav%Lat_Max_Limit)
        Space_Mask = sym%YES
   end where

   where(Nav%Lon_1b < Nav%Lon_Min_Limit .or. Nav%Lon_1b > Nav%Lon_Max_Limit)
        Space_Mask = sym%YES
   end where

   !--- Satzen limit
   where (Geo%Satzen > Geo%Satzen_Max_Limit .or. Geo%Satzen < Geo%Satzen_Min_Limit)
        Space_Mask = sym%YES
   end where

   !--- Solzen limit
   where (Geo%Solzen < Geo%Solzen_Min_Limit .or. Geo%Solzen > Geo%Solzen_Max_Limit .or. Geo%Satzen == Missing_Value_Real4)
        Space_Mask = sym%YES
   end where

end subroutine EXPAND_SPACE_MASK_FOR_USER_LIMITS
 
!======================================================================
! Check for solar contamination of what should be nighttime data
! 
! this is common in AVHRR and GOES.  Pixels with solar contamination
!
! are treated as bad pixel in the bad_pixel_mask
!
! Note that the Ch1_Counts used here are assumed have the dark count subtracted
! from them.
!
!======================================================================
subroutine SET_SOLAR_CONTAMINATION_MASK(Solar_Contamination_Mask) 

   logical, dimension(:,:), intent(out):: Solar_Contamination_Mask
   integer:: Number_of_Elements
   integer:: Number_of_Lines
   integer:: Elem_Idx
   integer:: Line_Idx
   integer, parameter:: Solar_Contamination_Thresh_AVHRR = 2
   integer, parameter:: Solar_Contamination_Thresh_GEO = 10

   Number_of_Elements = Image%Number_Of_Elements
   Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment

   !--- initialize 
   Solar_Contamination_Mask(:,1:Number_Of_Lines) = .FALSE.

   !---  loop through lines and elements
   line_loop: do Line_Idx = 1, Number_of_Lines
      element_loop: do Elem_Idx = 1, Number_of_Elements

        !--- check for solar contamination of nighttime data in AVHRR
        if (Sensor%Chan_On_Flag_Default(1) ) then

          if (index(Sensor%Sensor_Name,'AVHRR') > 0) then

            if ((Geo%Solzen(Elem_Idx,Line_Idx) > 90.0) .and. (Geo%Scatangle(Elem_Idx,Line_Idx) < 60.0)) then
              if (therm_cal_1b == sym%NO) then
                if (Ch1_Counts(Elem_Idx,Line_Idx) > Solar_Contamination_Thresh_AVHRR) then 
                   Solar_Contamination_Mask(Elem_Idx,Line_Idx) = .TRUE.
                endif
              else
                if (Ch1_Counts(Elem_Idx,Line_Idx) > Solar_Contamination_Thresh_AVHRR) then 
                   Solar_Contamination_Mask(Elem_Idx,Line_Idx) = .TRUE.
                endif 
              endif

            endif
  
          endif

          !--- check for solar contamination of nighttime data in GOES
          if (index(Sensor%Sensor_Name,'GOES') > 0) then
             if ((Geo%Solzen(Elem_Idx,Line_Idx) > 90.0) .and. (Geo%Scatangle(Elem_Idx,Line_Idx) < 60.0)) then
                if (Ch1_Counts(Elem_Idx,Line_Idx) > Solar_Contamination_Thresh_GEO) then
                   Solar_Contamination_Mask(Elem_Idx,Line_Idx) = .TRUE.
                endif 
             endif
          endif

        endif

        !--- check for solar contamination of nighttime data in GOES
        if (index(Sensor%Sensor_Name,'GOES') > 0) then
          if ((Geo%Solzen(Elem_Idx,Line_Idx) > 90.0) .and.  (Geo%Scatangle(Elem_Idx,Line_Idx) < 180.0)) then
            if (Ch1_Counts(Elem_Idx,Line_Idx)  > Solar_Contamination_Thresh_GEO) then
              Solar_Contamination_Mask(Elem_Idx,Line_Idx) = .TRUE.
            endif
          endif
        endif


        !until a more robust fix can be found, for now we will set all night
        ! pixels to have Solar_Contamination_Mask = sym%YES due to telescope
        ! contamination in 3.9um channel. Exists for FY-2C/D/E/G - WCS3
        
        if (index(Sensor%Sensor_Name,'FY2-IMAGER') > 0) then
          if ((Geo%Solzen(Elem_Idx,Line_Idx) > 90.0)) then
              Solar_Contamination_Mask(Elem_Idx,Line_Idx) = .TRUE.
          endif
        endif

      end do element_loop
   end do line_loop

end subroutine SET_SOLAR_CONTAMINATION_MASK

!======================================================================
! Check for bad pixels
!
! Apply tests to detect pixels that should not be processed.
!
!
! Bad_Pixel_Mask is meant to single mask that captures all reasons
! why a pixel should be skipped. This includes
!
! Space_Mask
! Solar_Contamination_Mask
! Missing NWP
! any one of channel tests
!
! Also, if many pixels in a line are bad, the whole line is set to bad
!
! Output:  Bad_Pixel_Mask 
! Input: taken from pixel common
!  
!======================================================================
	subroutine SET_BAD_PIXEL_MASK(Bad_Pixel_Mask)

   	logical, dimension(:,:), intent(out):: Bad_Pixel_Mask
   	integer:: Number_of_Elements
   	integer:: Number_of_Lines
   	integer:: Elem_Idx
   	integer:: Line_Idx
   	integer:: Number_Bad_Pixels
   	integer:: Number_Bad_Pixels_Thresh
   	integer:: Lon_Nwp_Idx
   	integer:: Lat_Nwp_Idx

   	Number_of_Elements = Image%Number_Of_Elements
   	Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment

   	Number_Bad_Pixels_Thresh = 0.9 * Image%Number_Of_Elements

		!----------------------------------------------------------------------
		!--- assign bad pixel mask based on scanline fatal flag
		!----------------------------------------------------------------------
   	!---- this is needed to ensure extra lines (beyond Num_Scans_Read) are bad
   	Bad_Pixel_Mask = .TRUE.

   	line_loop: do Line_Idx = 1, Number_of_Lines

      	!--- initialize
      	Bad_Pixel_Mask(:,Line_Idx) =.FALSE.

      	!--- check for a bad scan
      	if (Bad_Scan_Flag(Line_Idx) ) then

       		Bad_Pixel_Mask(:,Line_Idx) = .TRUE.

      	else

      		!--- if not a bad scan, check pixels on this scan
      		element_loop: do Elem_Idx = 1, Number_of_Elements

         		Bad_Pixel_Mask(Elem_Idx,Line_Idx) = Space_Mask(Elem_Idx,Line_Idx) == 1

         		!--- NaN checks on geolocation and geometry
         		if (isnan(Geo%Satzen(Elem_Idx,Line_Idx)) .or.  &
            		isnan(Geo%Solzen(Elem_Idx,Line_Idx))) then
            		Bad_Pixel_Mask(Elem_Idx,Line_Idx) = .TRUE.
         		endif

         		!--- Bad Relazimuth
         		if (Geo%Relaz(Elem_Idx,Line_Idx) == Missing_Value_Real4) then
            		Bad_Pixel_Mask(Elem_Idx,Line_Idx) = .TRUE.
         		endif

         		if ((Sensor%Chan_On_Flag_Default(31) ) .and. &
          			(Sensor%Chan_On_Flag_Default(32) )) then

          			!--- CALL any scan with a ridiculous pixel as bad 
          			!--- this is attempt data like NOAA-16 2004 023 where
          			!--- large fractions of scans are bad but not flagged as so
          			if (abs(ch(31)%Bt_Toa(Elem_Idx,Line_Idx) - ch(32)%Bt_Toa(Elem_Idx,Line_Idx)) > 20.0) then
               		Bad_Pixel_Mask(Elem_Idx,Line_Idx) = .TRUE.
          			endif

         		endif

         		!--- missing 11 um observations
         		if (Sensor%Chan_On_Flag_Default(31) ) then

         			if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < 150.0) then
            			Bad_Pixel_Mask(Elem_Idx,Line_Idx) = .TRUE.
         			endif

         			if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) > 350.0) then
            			Bad_Pixel_Mask(Elem_Idx,Line_Idx) = .TRUE.
         			endif

         			if (isnan(ch(31)%Bt_Toa(Elem_Idx,Line_Idx))) then
            			Bad_Pixel_Mask(Elem_Idx,Line_Idx) = .TRUE.
         			endif

         		endif

        

         		!--- check for solar zenith angle limits
         		if ((Geo%Solzen(Elem_Idx,Line_Idx) < Geo%Solzen_Min_Limit) .or. &
            		(Geo%Solzen(Elem_Idx,Line_Idx) > Geo%Solzen_Max_Limit)) then
             		Bad_Pixel_Mask(Elem_Idx,Line_Idx) = .TRUE.
         		endif

         		!--- NWP
         		if (Nwp_Opt /= 0) then
            		Lon_Nwp_Idx = i_Nwp(Elem_Idx,Line_Idx)
            		Lat_Nwp_Idx = j_Nwp(Elem_Idx,Line_Idx)
            		if (Lon_Nwp_Idx < 1 .or. Lat_Nwp_Idx < 1) then
                  	Bad_Pixel_Mask(Elem_Idx,Line_Idx) = .TRUE.
            		else
                  	if (Bad_Nwp_Mask(Lon_Nwp_Idx, Lat_Nwp_Idx) ==  sym%YES) Bad_Pixel_Mask(Elem_Idx,Line_Idx) = .TRUE.
            		endif
       	 		endif

      		end do element_loop

      	end if



		!TOCHECK -- is this still needed? (AW 01/05/2017)
      !---- consider any scanline with any solar contamination as a bad line
		!    if (maxval(Solar_Contamination_Mask(:,Line_Idx)) ) then
		!        Bad_Scan_Flag(Line_Idx) = sym%YES
		!        Bad_Pixel_Mask(:,Line_Idx) = sym%YES
		!    endif
     

   	end do line_loop

   	!-----------------------------------------------------------------------------------
   	! if the IDPS cloud mask is to be used for product generation, make sure that
   	! pixels within the gaps are considered bad
   	!-----------------------------------------------------------------------------------
   	if (Cloud_Mask_Aux_Flag == sym%USE_AUX_CLOUD_MASK .and. trim(Sensor%Sensor_Name) == 'VIIRS') then
      	where(Gap_Pixel_Mask == sym%YES )
         	Bad_Pixel_Mask = .TRUE.
      	end where
   	endif

   	!---------------------------------------------------------------------------------------
   	! Compute the fraction of the segment covered by valid data
   	!---------------------------------------------------------------------------------------
   	Segment_Valid_Fraction = 1.0 - count(Bad_Pixel_Mask(:,1:Number_of_Lines)) /  &
                                float(Number_of_Elements * Number_of_Lines)

	end subroutine SET_BAD_PIXEL_MASK
   
   !--------------------------------------------------------------------------
   !QUALITY_CONTROL_ANCILLARY_DATA
   !
   ! Apply some checks on ancillary data.  Call pixels bad when checks fail
   !
   ! Note, Bad_Pixel_Mask is modified but not created here. Do not initialize
   !  remove the loop (AW 07 March 2017)
   !--------------------------------------------------------------------------
   subroutine QUALITY_CONTROL_ANCILLARY_DATA (Bad_Pixel_Mask)
      logical, dimension(:,:), intent(inout):: Bad_Pixel_Mask
   
   
      where ( Sfc%Sfc_Type < 0. .OR. Sfc%Sfc_Type > 15)
         Bad_Pixel_Mask = .TRUE.
      end where
   

   end subroutine QUALITY_CONTROL_ANCILLARY_DATA




   !--------------------------------------------------------------------------
   ! COMPUTE_PIXEL_ARRAYS
   !
   ! compute quantities needed for other routines, (ie cloud mask)
   !
   ! input - none
   !
   ! output - only into shared memory
   !   seczen - secant of zenith angle
   !   Sst_Unmasked - an Sst used for cloud masking only
   !   Btd_Ch31_Ch32 - brightness temperature difference between Bt_Ch31 and Bt_Ch32 
   !   Btd_Ch20_Ch31 - brightness temperature difference between Bt_Ch20 and Bt_Ch31 
   !   Ref_ratio_Ch6_Ch1 - ratio of Ref_Ch6 / Ref_Ch1
   !   Ref_ratio_Ch20_Ch1 - ratio of Ref_Ch20 / Ref_Ch1
   !--------------------------------------------------------------------------
   subroutine COMPUTE_PIXEL_ARRAYS(j1,nj)

      integer, intent(in):: j1,nj
      integer(kind=int4):: j2

      j2 = j1 + nj - 1

      Geo%Coszen(:,j1:j2) = cos(Geo%Satzen(:,j1:j2)*dtor)

      Geo%Seczen(:,j1:j2) = 1.0 / Geo%Coszen(:,j1:j2)

      Geo%Cossolzen(:,j1:j2) = cos(Geo%Solzen(:,j1:j2)*dtor)

      Geo%Airmass(:,j1:j2) = 1.0
      
      where(Geo%Solzen(:,j1:j2) /= 0.0 .and. Geo%Coszen(:,j1:j2) /= 0.0) 
         Geo%Airmass(:,j1:j2) = 1.0/Geo%Cossolzen(:,j1:j2) + 1.0/Geo%Coszen(:,j1:j2)
      end where

      !--- other useful arrays 

      !--- channel 31 and channel 32 brightness temperature difference
      if ((Sensor%Chan_On_Flag_Default(31) ) .and. &
            (Sensor%Chan_On_Flag_Default(32) )) then
         Btd_Ch31_Ch32(:,j1:j2) = ch(31)%Bt_Toa(:,j1:j2) - ch(32)%Bt_Toa(:,j1:j2)
      end if

      !--- channel 20 and channel 31 brightness temperature difference
      if (Sensor%Chan_On_Flag_Default(20)  .and. Sensor%Chan_On_Flag_Default(31) ) then
         Btd_Ch20_Ch31 = ch(20)%Bt_Toa - ch(31)%Bt_Toa
         where(ch(20)%Bt_Toa == Missing_Value_Real4 .or. ch(31)%Bt_Toa == Missing_Value_Real4)
            Btd_Ch20_Ch31 = Missing_Value_Real4
         end where
      end if

      !--- channel 20 and channel 32 brightness temperature difference
      if (Sensor%Chan_On_Flag_Default(20)  .and. Sensor%Chan_On_Flag_Default(32) ) then
         Btd_Ch20_Ch32 = ch(20)%Bt_Toa - ch(32)%Bt_Toa
         where(ch(20)%Bt_Toa == Missing_Value_Real4 .or. ch(32)%Bt_Toa == Missing_Value_Real4)
            Btd_Ch20_Ch32 = Missing_Value_Real4
         end where
      end if

   end subroutine COMPUTE_PIXEL_ARRAYS

!-------------------------------------------------------------------------------
!--- populate the snow_class array based on all available sources of Snow data
!--
!--- Input:
!---  NWP_Wat_Eqv_Snow_Depth - water equivalent snow depth from nwp
!---  NWP_Sea_Ice_Frac - sea ice fracion from nwp
!---  SST_Sea_Ice_Frac - sea ice fracion from sst data source
!---  Snow_Class_IMS - high resolution snow class field (highest priority)
!---  Snow_Class_Global - ESA GlobSnow products (lower priority)
!---
!--- Output:
!---  Snow_Class_Final - final classificiation
!---
!--- Symbology:
!---  1 = sym%NO_SNOW
!---  2 = sym%SEA_ICE
!---  3 = sym%SNOW
!-------------------------------------------------------------------------------
 subroutine COMPUTE_SNOW_CLASS(Snow_Class_NWP, Snow_Class_OISST, Snow_Class_IMS, &
                               Snow_Class_Glob,Land_Class,Snow_Class_Final)
 
   integer(kind=int1), intent(in), dimension(:,:):: Snow_Class_NWP
   integer(kind=int1), intent(in), dimension(:,:):: Snow_Class_OISST
   integer(kind=int1), intent(in), dimension(:,:):: Snow_Class_IMS
   integer(kind=int1), intent(in), dimension(:,:):: Snow_Class_Glob
   integer(kind=int1), intent(in), dimension(:,:):: Land_Class
   integer(kind=int1), intent(out), dimension(:,:):: Snow_Class_Final
   integer(kind=int1):: Finished_Flag

   Snow_Class_Final = Missing_Value_Int1

   Finished_Flag = 0

   do while (Finished_Flag == 0)

      !--- High Res
      if (Read_Snow_Mask == sym%READ_SNOW_HIRES .and.     &
          Failed_IMS_Snow_Mask_flag == sym%NO) then
          Snow_Class_Final = Snow_Class_IMS
          Finished_Flag = 1
      endif

      !-- GlobSnow - does not work
      if (Read_Snow_Mask == sym%READ_SNOW_GLOB .and.   &
          Failed_Glob_Snow_Mask_Flag == sym%NO) then
          Snow_Class_Final = Snow_Class_Glob
          Finished_Flag = 1
      endif

      Snow_Class_Final = Snow_Class_Nwp

      !--- overwrite with oisst
      where(Snow_Class_OISST == sym%SEA_ICE)
        Snow_Class_Final = Snow_Class_OISST
      endwhere
      Finished_Flag = 1
       
   enddo

   !-- check for consistnecy of land and snow masks
   where(Snow_Class_Final == sym%SNOW .and. Land_Class /= sym%LAND)
             Snow_Class_Final = sym%SEA_ICE
   endwhere
   where(Snow_Class_Final == sym%SEA_ICE .and. Land_Class == sym%LAND)
             Snow_Class_Final = sym%SNOW
   endwhere

   !---- remove snow under certain conditions

   !-- can't be snow if warm
   if (Sensor%Chan_On_Flag_Default(31) ) then
    where (ch(31)%Bt_Toa > 277.0)
        Snow_Class_Final = sym%NO_SNOW
    endwhere
   endif

   !--- some day-specific tests
   if (Sensor%Chan_On_Flag_Default(1) ) then
    where (ch(1)%Ref_Toa < 10.0 .and. Geo%Solzen < 75.0)
        Snow_Class_Final = sym%NO_SNOW
    endwhere
   endif

 end subroutine COMPUTE_SNOW_CLASS
 !---------------------------------------------------------------------------------------
 ! Compute a Snow Classification from the NWP
 ! 
 ! threshold are empirically derived by comparing NWP fields to IMS (A.  Heidinger)
 !---------------------------------------------------------------------------------------
 subroutine COMPUTE_SNOW_CLASS_NWP(NWP_Wat_Eqv_Snow_Depth,NWP_Sea_Ice_Frac, Snow_Class_Nwp)

   real(kind=real4), intent(in), dimension(:,:):: NWP_Wat_Eqv_Snow_Depth,  &
                                                  NWP_Sea_Ice_Frac
   integer(kind=int1), intent(out), dimension(:,:):: Snow_Class_NWP

   !--- initialize all pixels as missing
   Snow_Class_NWP = Missing_Value_Int1

   !--- initialize pixel with valid data as no_snow
   where (NWP_Wat_Eqv_Snow_Depth >= 0.0 .or. NWP_Sea_Ice_Frac >=0.0) 
          Snow_Class_NWP = sym%NO_SNOW
   end where 

   !--- detect sea ice
   where (NWP_Sea_Ice_Frac > 0.5) 
          Snow_Class_NWP = sym%SEA_ICE
   end where

   !--- detect snow  (note snow can cover sea-ice so do this after sea-ice check)
   where (NWP_Wat_Eqv_Snow_Depth > 0.1) 
          Snow_Class_NWP = sym%SNOW
   end where 

 end subroutine COMPUTE_SNOW_CLASS_NWP
 !---------------------------------------------------------------------------------------
 ! Compute a Sea-Ice Classification from OISST Analysis
 ! Note =- OISST only provides Sea Ice, not Snow
 !---------------------------------------------------------------------------------------
 subroutine COMPUTE_SNOW_CLASS_OISST(SST_Sea_Ice_Frac, Snow_Class_OISST)

   real(kind=real4), intent(in), dimension(:,:):: SST_Sea_Ice_Frac
   integer(kind=int1), intent(out), dimension(:,:):: Snow_Class_OISST

   !--- initialize all to missing
   Snow_Class_OISST = Missing_Value_Int1

   !--- initialize valid values as no_snow
   where (SST_Sea_Ice_Frac >= 0.0) 
          Snow_Class_OISST = sym%NO_SNOW
   end where
   !--- detect sea ice
   where (SST_Sea_Ice_Frac > 0.5) 
          Snow_Class_OISST = sym%SEA_ICE
   end where

 end subroutine COMPUTE_SNOW_CLASS_OISST

!-------------------------------------------------------------------------------
!--- populate the snow_class array based on all available sources of Snow data
!--
!--- Input:
!---  NWP_Wat_Eqv_Snow_Depth - water equivalent snow depth from nwp
!---  NWP_Sea_Ice_Frac - sea ice fracion from nwp
!---  SST_Sea_Ice_Frac - sea ice fracion from sst data source
!---  Snow_Class_IMS - high resolution snow class field (highest priority)
!---  Snow_Class_Global - ESA GlobSnow products (lower priority)
!---
!--- Output:
!---  Snow_Class_Final - final classificiation
!---
!--- Symbology:
!---  1 = sym%NO_SNOW
!---  2 = sym%SEA_ICE
!---  3 = sym%SNOW
!-------------------------------------------------------------------------------
 subroutine COMPUTE_SNOW_FIELD(NWP_Wat_Eqv_Snow_Depth,NWP_Sea_Ice_Frac, SST_Sea_Ice_Frac,  &
                               Snow_Class_IMS,Snow_Class_Glob,Snow_Class_Final)

   real(kind=real4), intent(in), dimension(:,:):: NWP_Wat_Eqv_Snow_Depth,  &
                                                  NWP_Sea_Ice_Frac, SST_Sea_Ice_Frac
   integer(kind=int1), intent(in), dimension(:,:):: Snow_Class_IMS, Snow_Class_Glob
   integer(kind=int1), intent(out), dimension(:,:):: Snow_Class_Final
   integer(kind=int4):: Elem_Idx,Line_Idx,ihires

   Snow_Class_Final = Missing_Value_Int1

   line_loop: do Line_Idx = 1, Image%Number_Of_Lines_Read_This_Segment
     element_loop: do Elem_Idx= 1, Image%Number_Of_Elements

      !--- check for bad scans
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) ) then
       cycle
      endif

      !--- initialize valid pixels to NO_SNOW
      Snow_Class_Final(Elem_Idx,Line_Idx) = sym%NO_SNOW

      !--- if hires field is available, use it and ignore other sources
      ihires = sym%NO
      if (Read_Snow_Mask == sym%READ_SNOW_HIRES .and.     &
          Failed_IMS_Snow_Mask_flag == sym%NO) then
          Snow_Class_Final(Elem_Idx,Line_Idx) = Snow_Class_IMS(Elem_Idx,Line_Idx)
          ihires = sym%YES
      endif

      if (ihires == sym%NO) then

        !use nwp and/or Sst analysis
        if (Nwp_Opt > 0) then
            if (NWP_Wat_Eqv_Snow_Depth(Elem_Idx,Line_Idx) > 0.1) then  !this is Snow depth
                Snow_Class_Final(Elem_Idx,Line_Idx) = sym%SNOW
             endif
             if (.NOT. Use_Sst_Anal .and. NWP_Sea_Ice_Frac(Elem_Idx,Line_Idx) > 0.5) then
               Snow_Class_Final(Elem_Idx,Line_Idx) = sym%SEA_ICE
             endif
        endif

         !--- GlobSnow (Land Only)
         if (Read_Snow_Mask == sym%READ_SNOW_GLOB .and.   &
          Failed_Glob_Snow_Mask_Flag == sym%NO  .and.       &
          Sfc%Land(Elem_Idx,Line_Idx) == sym%LAND  .and.       &
          Nav%Lat(Elem_Idx,Line_Idx) >= 35.0 .and. &     !note, globSnow is NH only (35-85)
          Nav%Lat(Elem_Idx,Line_Idx) <= 85.0) then

            Snow_Class_Final(Elem_Idx,Line_Idx)  = Snow_Class_Glob(Elem_Idx,Line_Idx)
         endif  

         !-- correct for GlobSnow missing snow over Greenland and Arctic Islands
         !-- that have a barren surface type
         if ((Sfc%Sfc_Type(Elem_Idx,Line_Idx) == sym%Bare_Sfc .or. Sfc%Sfc_Type(Elem_Idx,Line_Idx) == sym%OPEN_SHRUBS_SFC) .and. &
             (abs(Nav%Lat(Elem_Idx,Line_Idx)) > 60.0) .and. &
             (NWP_Wat_Eqv_Snow_Depth(Elem_Idx,Line_Idx) > 0.1)) then
                Snow_Class_Final(Elem_Idx,Line_Idx) = sym%SNOW
         endif

         !-- correct for GlobSnow missing snow over mountains in v1.0
         !-- this will be fixed in future versions
         if ((Sfc%Zsfc(Elem_Idx,Line_Idx) > 1000.0) .and.  (NWP_Wat_Eqv_Snow_Depth(Elem_Idx,Line_Idx) > 0.001)) then
                Snow_Class_Final(Elem_Idx,Line_Idx) = sym%SNOW
         endif

         !--- pick up small lakes in regions marked snowy in glob snow 
         if ((Sfc%Land(Elem_Idx,Line_Idx) /= sym%LAND) .and. (Snow_Class_Glob(Elem_Idx,Line_Idx) == sym%SNOW .or. Snow_Class_Glob(Elem_Idx,Line_Idx) == sym%SEA_ICE)) then
          Snow_Class_Final(Elem_Idx,Line_Idx)  = sym%SEA_ICE
         endif


         !--- check for sea-ice from Sst analysis (OISST in this case)
         !--- and replace values from nwp
         !--- note, threshold 0.5 determined by matching IMS

         !--- under these conditions believe SST Analysis
         !--- and allow it to remove ice 
         if (use_Sst_anal  .and. &
             (Sfc%Land(Elem_Idx,Line_Idx) == sym%DEEP_OCEAN .or.     &
              Sfc%Land(Elem_Idx,Line_Idx) == sym%DEEP_INLAND_WATER .or. &
              Sfc%Land(Elem_Idx,Line_Idx) == sym%MODERATE_OCEAN .or. &
              Sfc%Land(Elem_Idx,Line_Idx) == sym%SHALLOW_OCEAN))  then

           if (Sst_Sea_Ice_Frac(Elem_Idx,Line_Idx) > 0.50) then     
             Snow_Class_Final(Elem_Idx,Line_Idx) = sym%SEA_ICE
           else
             Snow_Class_Final(Elem_Idx,Line_Idx) = sym%NO_SNOW
           endif

         endif

         !--- OISST does not capture sea-ice in shallow water near AntArctica
         if ((Sfc%Land(Elem_Idx,Line_Idx) == sym%SHALLOW_OCEAN) .and. &
             (Nav%Lat(Elem_Idx,Line_Idx) < -60.0) .and. &
             ((NWP_Sea_Ice_Frac(Elem_Idx,Line_Idx) > 0.50).or. &
              (NWP_Wat_Eqv_Snow_Depth(Elem_Idx,Line_Idx)>0.1))) then
                Snow_Class_Final(Elem_Idx,Line_Idx) = sym%SEA_ICE
         endif


         !--- allow sst analysis to detect any ice (but not remove it)
         !--- this is needed since Lake Erie is cover the OISST but
         !--- classified as shallow_inland_water
         if ((Sfc%Land(Elem_Idx,Line_Idx) /= sym%LAND) .and. (Sst_Sea_Ice_Frac(Elem_Idx,Line_Idx) > 0.50)) then     
             Snow_Class_Final(Elem_Idx,Line_Idx) = sym%SEA_ICE
         endif
          
         !-- final check for consistnecy of land and snow masks
         if ((Snow_Class_Final(Elem_Idx,Line_Idx) == sym%SNOW) .and. (Sfc%Land(Elem_Idx,Line_Idx) /= sym%LAND)) then
             Snow_Class_Final(Elem_Idx,Line_Idx) = sym%SEA_ICE
         endif
         if ((Snow_Class_Final(Elem_Idx,Line_Idx) == sym%SEA_ICE) .and. (Sfc%Land(Elem_Idx,Line_Idx) == sym%LAND)) then
             Snow_Class_Final(Elem_Idx,Line_Idx) = sym%SNOW
         endif
          
       endif

       !-------------------------------------------------------
       ! use some basic tests to modify Snow
       ! only do this when hires Snow field is unavailable
       !-------------------------------------------------------

       if (ihires == sym%NO) then

          !-- can't be snow if warm
         if (Sensor%Chan_On_Flag_Default(31) ) then
          if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) > 277.0) then
           Snow_Class_FInal(Elem_Idx,Line_Idx) = sym%NO_SNOW
          endif
         endif

         !--- some day-specific tests
         if (Geo%Solzen(Elem_Idx,Line_Idx) < 75.0) then  ! day check

         !--- conditions where Snow is not possible
           if (Sensor%Chan_On_Flag_Default(1) ) then
            if(ch(1)%Ref_Toa(Elem_Idx,Line_Idx) < 10.0) then
             Snow_Class_Final(Elem_Idx,Line_Idx) = sym%NO_SNOW
            endif
           endif

         endif           !day check

       endif            !hires Snow check

    end do element_loop
 end do   line_loop

 end subroutine COMPUTE_SNOW_FIELD

!------------------------------------------------------------------
! Compute a surface temperature by using the observed 
! 11 micron radiance, the rtm calculations and the surface emissivity
!
! July 2009 - made AVHRR/1 compliant
!
! Author: Andrew Heidinger
! 
!------------------------------------------------------------------
 subroutine COMPUTE_TSFC(jmin,jmax)

  integer, intent(in):: jmin
  integer, intent(in):: jmax

  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: xnwp
  integer:: ynwp

  real:: Rad11
  real:: Rad11_Atm
  real:: Trans11_Atm
  real:: Rad11_Atm_Dwn_Sfc
  real:: Emiss_Sfc11
  real:: Rad11_Sfc
  real:: B11_Sfc

  !------------------------------------------------------------------
  ! Loop over pixels and derive surface temp
  !------------------------------------------------------------------

  !--- initialize
  Tsfc_Retrieved = Missing_Value_Real4
  Trad_Retrieved = Missing_Value_Real4
  Tsfc_Qf = 0

  !--- if no ch31, abort
  if ( .NOT. Sensor%Chan_On_Flag_Default(31) ) then
     return 
  endif

  line_loop: do Line_Idx=jmin, jmax - jmin + 1
    element_loop: do Elem_Idx= 1, Image%Number_Of_Elements

      !--- check for a bad pixel pixel
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) ) then
        cycle
      endif

      !--- aliases for visual convenience
      Rad11 = ch(31)%Rad_Toa(Elem_Idx,Line_Idx)
      Xnwp = i_nwp(Elem_Idx,Line_Idx)                         !nwp latitude cell
      Ynwp = j_nwp(Elem_Idx,Line_Idx)                         !nwp longitude cell
      Rad11_Atm = ch(31)%Rad_Atm(Elem_Idx,Line_Idx)           !11 micron atmospheric radiance
      Trans11_Atm = ch(31)%Trans_Atm(Elem_Idx,Line_Idx)       !11 micron atmospheric transmittance
      Emiss_Sfc11 = ch(31)%Sfc_Emiss(Elem_Idx,Line_Idx)       !11 micron surface emissivity
      Rad11_Atm_Dwn_Sfc = ch(31)%Rad_Atm_Dwn_Sfc(Elem_Idx,Line_Idx)  !11 micron atmospheric radiance down at sfc

      !--- compute the radiance coming off the surface
      Rad11_Sfc = (Rad11 - Rad11_Atm) / Trans11_Atm - &
                  (1.0-Emiss_Sfc11)*Rad11_Atm_Dwn_Sfc

      !--- compute to a temperature
      Trad_Retrieved(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(31,Rad11_Sfc)

      !--- adjust for surface emissivity - this is now the black body emission at Tsfc
      B11_Sfc = Rad11_Sfc / Emiss_Sfc11

      !--- compute to a temperature
      Tsfc_Retrieved(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(31,B11_Sfc)

    end do element_loop
  end do line_loop

end subroutine COMPUTE_TSFC
!------------------------------------------------------------
! atmospheric correction
!
! note, optical depths taken from
!  Ignatov and Stowe, 2002, vol 59, JAS, pg 313-334
!  I used the MLW values
!-----------------------------------------------------------
subroutine ATMOS_CORR(Line_Idx_Min,Num_Lines)

   integer, intent(in):: Line_Idx_Min
   integer, intent(in):: Num_Lines
   integer:: Line_Idx
   integer:: Elem_Idx
   integer:: Line_Idx_Max
   integer:: Elem_Idx_Max
   integer:: Elem_Idx_Min
   integer:: Num_Elements
   integer:: Chan_Idx
   
   real:: Ref_ss
   real:: Albedo_View
   real:: Albedo_Sun
   real:: Trans_Total
   real:: Tau_Total
   real:: Tau_Ray
   real:: Tau_Gas
   real:: Tau_H2O
   real:: Tau_Aer
   real:: Wo_Aer
   real:: G_Aer

   real:: Source_Zen
   real:: Cos_Source_Zen
   real:: Airmass_Factor
   real:: Scattering_Angle

   Elem_Idx_Min = 1
   Num_Elements = Image%Number_Of_Elements
   Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
   Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  !---------- loop over scan lines
  line_loop: do Line_Idx = Line_Idx_Min, Line_Idx_Max

  element_loop: do Elem_Idx = Elem_Idx_Min, Elem_Idx_Max
      if ( line_idx .ne. 40 .or. elem_idx .ne. 70) cycle
     !--- check for bad individual pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) ) cycle

     channel_loop: do Chan_Idx = 1,Nchan_Clavrx

       if (.NOT. Sensor%Chan_On_Flag_Default(Chan_Idx)) cycle

       if (Chan_Idx /= 1 .and. Chan_Idx /= 2 .and. Chan_Idx /= 5 .and. &
           Chan_Idx /= 6 .and. Chan_Idx /= 7 .and. Chan_Idx /= 44) cycle

       !--- check for valid data
       if (Chan_Idx /= 44 .and. ch(Chan_Idx)%Ref_Toa(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
       if (Chan_Idx == 44) then
         if (ch(Chan_Idx)%Ref_Lunar_Toa(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
       endif

       !--- set source angle
       Source_Zen = Geo%SolZen(Elem_Idx,Line_Idx)
       Scattering_Angle = Geo%Scatangle(Elem_Idx,Line_Idx)
       if (Chan_Idx == 44) then
            Source_Zen = Geo%LunZen(Elem_Idx,Line_Idx)
            Scattering_Angle = Geo%Scatangle_Lunar(Elem_Idx,Line_Idx)
       endif

       !--- check for appropriate illumination
       if (Source_Zen >= 90.0) cycle

       Tau_H2O = Solar_Rtm%Tau_H2O_Coef(Chan_Idx,1) + Solar_Rtm%Tau_H2O_Coef(Chan_Idx,2)*Tpw_Nwp_Pix(Elem_Idx,Line_Idx) +  &
                   Solar_Rtm%Tau_H2O_Coef(Chan_Idx,3)*(Tpw_Nwp_Pix(Elem_Idx,Line_Idx)**2)
       Tau_Gas = max(0.0,Tau_H2O) + Solar_Rtm%Tau_O3(Chan_Idx) + Solar_Rtm%Tau_O2(Chan_Idx) &
               + Solar_Rtm%Tau_CO2(Chan_Idx) + Solar_Rtm%Tau_CH4(Chan_Idx)
       Tau_Aer = Solar_Rtm%Tau_Aer(Chan_Idx)
       Wo_Aer = Solar_Rtm%Wo_Aer(Chan_Idx)
       G_Aer = Solar_Rtm%G_Aer(Chan_Idx)
       Tau_Ray = Solar_Rtm%Tau_Ray(Chan_Idx)

       !------------------------------------------------------------------------
       ! select gas and surface reflectance parameters
       !------------------------------------------------------------------------
       select case (Chan_Idx)

       case(1)
         if (ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else 
              Albedo_View = Ch1_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif 
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch1_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(2)
         if (ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else 
              Albedo_View = Ch2_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif 
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch2_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(5)
         if (ch(5)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = ch(5)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else 
              Albedo_View = Ch6_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0 !Note there is no Ch5_Sfc_Alb_Umd
         endif 
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch5_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(6)
         if (ch(6)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = ch(6)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else 
              Albedo_View = Ch6_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif 
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch6_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(7)
         if (ch(7)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = ch(7)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else 
              Albedo_View = Ch6_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0   !Note there is no Ch7_Sfc_Alb_Umd
         endif 
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch7_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(44)  !DNB - use mean of ch1 and ch2 for sfc reflectance
         if (ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
            .and. ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = 0.5*(ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)+ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)) / 100.0
         else 
              Albedo_View = 0.5*(Ch1_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) + Ch2_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx))) / 100.0
         endif 
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = 0.5*(Ch1_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) &
               + Ch2_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx))) / 100.0
         endif
       end select

       !-- assume lambertian
       Albedo_Sun = Albedo_View

       Cos_Source_Zen = Missing_Value_Real4
       if (Source_Zen >= 0.0 .and. Source_Zen < 90.0) Cos_Source_Zen = cos(Source_Zen*DTOR)

       Airmass_Factor = 1.0 / Geo%CosZen(Elem_Idx,Line_Idx) + &
                        1.0 / Cos_Source_Zen

       !--- compute atmospheric scattering
       call COMPUTE_CLEAR_SKY_SCATTER(Tau_Aer, &
                                     Wo_Aer, &
                                     G_Aer, &
                                     Tau_Ray, &
                                     Tau_Gas, &
                                     Scattering_Angle, &
                                     Geo%Coszen(Elem_Idx,Line_Idx), &
                                     Cos_Source_Zen, &
                                     Albedo_View, &
                                     Albedo_Sun, &
                                     Ref_ss)

       !--- compute total transmission for combining terms
       Tau_Total = Tau_Aer + Tau_Ray + Tau_Gas
       Trans_Total = exp(-Tau_Total*Airmass_Factor)
       ! print*,'Routine ATNOS_CORR in PIXEL_ROTUINES: ',chan_idx,Tau_Aer , Tau_Ray , Tau_Gas, trans_total
       if (Chan_Idx /= 44) then

         !--- compute atmospherically corrected reflectance (at sfc level)s
         ch(Chan_Idx)%Ref_Sfc(Elem_Idx,Line_Idx) = (ch(Chan_Idx)%Ref_Toa(Elem_Idx,Line_Idx) - Ref_ss) / Trans_Total

         !--- compute top of clear-sky atmosphere reflectance
         ch(Chan_Idx)%Ref_Toa_Clear(Elem_Idx,Line_Idx) = Ref_ss + Trans_Total*100.0*Albedo_View

       else

         !--- compute atmospherically corrected reflectance (at sfc level)s
         ch(Chan_Idx)%Ref_Lunar_Sfc(Elem_Idx,Line_Idx) = (ch(Chan_Idx)%Ref_Lunar_Toa(Elem_Idx,Line_Idx) - Ref_ss) / Trans_Total

         !--- compute top of clear-sky atmosphere reflectance
         ch(Chan_Idx)%Ref_Lunar_Toa_Clear(Elem_Idx,Line_Idx) = Ref_ss + Trans_Total*100.0*Albedo_View

       endif
      
    end do  channel_loop

   end do element_loop

 end do line_loop

end subroutine ATMOS_CORR

!======================================================================
! Normalize the reflectances by the solar zenith angle cosine
!======================================================================
 	subroutine NORMALIZE_REFLECTANCES(Sun_Earth_Distance)
   	real(kind=real4), intent(in):: Sun_Earth_Distance
   	integer:: i,j, Chan_Idx
   	real:: Factor

   	! for these sensors, no correction is needed
   	if (trim(Sensor%Sensor_Name) == 'VIIRS') return 
   	if (trim(Sensor%Sensor_Name) == 'VIIRS-IFF') return 
   	if (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') return 

   	!--------------------------------------------------------------------
   	! loop through pixels and apply normalization factor
   	!--------------------------------------------------------------------
   	do j = 1, Image%Number_Of_Lines_Read_This_Segment

      	do i = 1, Image%Number_Of_Elements
			
				if (.NOT. bad_pixel_mask(i,j) .AND. geo % cossolzen(i,j) > 0.0 ) then
      

       			Factor = 1.0 / Geo%Cossolzen(i,j)

       			! for these sensors, a correction for sun earth distance is also needed
       			if ( (trim(Sensor%Sensor_Name) == 'AVHRR-1') .or. &
            		(trim(Sensor%Sensor_Name) == 'AVHRR-2') .or. &
            		(trim(Sensor%Sensor_Name) == 'AVHRR-3') .or. &
            		(trim(Sensor%Sensor_Name) == 'GOES-IL-IMAGER') .or. &
             		(trim(Sensor%Sensor_Name) == 'GOES-MP-IMAGER') .or. &
            		(trim(Sensor%Sensor_Name) == 'GOES-IP-SOUNDER') .or. &
            		(trim(Sensor%Sensor_Name) == 'SEVIRI') .or. &
            		(trim(Sensor%Sensor_Name) == 'MTSAT-IMAGER') .or. &
            		(trim(Sensor%Sensor_Name) == 'COMS-IMAGER') .or. &
            		(trim(Sensor%Sensor_Name) == 'FY2-IMAGER')) then

            			Factor = Factor * (Sun_Earth_Distance**2) 

       			endif

       			! apply correction
       			do Chan_Idx = 1,19
          			if (Sensor%Chan_On_Flag_Default(Chan_Idx) ) then
            			if (ch(Chan_Idx)%Ref_Toa(i,j) /= Missing_Value_Real4) then
             				ch(Chan_Idx)%Ref_Toa(i,j) = ch(Chan_Idx)%Ref_Toa(i,j) * Factor
            			endif
          			endif
       			enddo
				
       			if (Sensor%Chan_On_Flag_Default(26) ) then
          			if (ch(26)%Ref_Toa(i,j) /= Missing_Value_Real4) then
            			ch(26)%Ref_Toa(i,j) = ch(26)%Ref_Toa(i,j) * Factor
          			endif
       			endif

       			!  normalize by sun angle the dark sky composite
       			if ((Sensor%Chan_On_Flag_Default(1) )) then
         			if (Ref_Ch1_Dark_Composite(i,j) /= Missing_Value_Real4) then
            			Ref_Ch1_Dark_Composite(i,j) = Ref_Ch1_Dark_Composite(i,j) * Factor
         			endif
       			endif


       			!--- for avhrr, handle absense of ch6
       			if (index(Sensor%Sensor_Name,'AVHRR') > 0 .and. &
            		Sensor%Chan_On_Flag_Default(6) ) then
            		if (ch(6)%Ref_Toa(i,j) < 0) ch(6)%Ref_Toa(i,j) = Missing_Value_Real4
            		if (Ch3a_On_Avhrr(j) /= sym%YES) ch(6)%Ref_Toa(i,j) = Missing_Value_Real4
       			endif

       			!--- in terminator region, renormalize Channel 1 (maybe extend to all?)
       			if (Geo%Solzen(i,j) > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) then
          			if (Sensor%Chan_On_Flag_Default(1) )  &
               		ch(1)%Ref_Toa(i,j) = TERM_REFL_NORM(Geo%Cossolzen(i,j),ch(1)%Ref_Toa(i,j))
       			endif

       			!---- make unnormalized relflectances for AWIPS display ( can we retire this)
         		if (Sensor%Chan_On_Flag_Default(1) )  then
          			if (ch(1)%Ref_Toa(i,j) /= Missing_Value_Real4) then
            			ch(1)%Ref_Toa_Unnorm(i,j) = ch(1)%Ref_Toa(i,j) * Geo%Cossolzen(i,j)
          			endif
         		endif

         		if (Sensor%Chan_On_Flag_Default(2) ) then
          			if (ch(2)%Ref_Toa(i,j) /= Missing_Value_Real4) then
            			ch(2)%Ref_Toa_Unnorm(i,j) = ch(2)%Ref_Toa(i,j) * Geo%Cossolzen(i,j)
          			endif
         		endif

         		if (Sensor%Chan_On_Flag_Default(6) ) then
          			if (ch(6)%Ref_Toa(i,j) /= Missing_Value_Real4) then
            			ch(6)%Ref_Toa_Unnorm(i,j) = ch(6)%Ref_Toa(i,j) * Geo%Cossolzen(i,j)
          			endif
         		endif

      		else

       			!--- set to missing
       			do Chan_Idx = 1,19
          			if (Sensor%Chan_On_Flag_Default(Chan_Idx) ) ch(Chan_Idx)%Ref_Toa(i,j) = Missing_Value_Real4
       			enddo
       			if (Sensor%Chan_On_Flag_Default(26) ) ch(26)%Ref_Toa(i,j) = Missing_Value_Real4

       			!--- un-normalized refs for AWIPS
       			if (Sensor%Chan_On_Flag_Default(1) ) ch(1)%Ref_Toa_Unnorm(i,j) = Missing_Value_Real4
       			if (Sensor%Chan_On_Flag_Default(2) ) ch(2)%Ref_Toa_Unnorm(i,j) = Missing_Value_Real4
       			if (Sensor%Chan_On_Flag_Default(6) ) ch(6)%Ref_Toa_Unnorm(i,j) = Missing_Value_Real4
					
				end if


      	end do
   
   	end do

	end subroutine NORMALIZE_REFLECTANCES
!----------------------------------------------------------------------
! Compute Channel3b albedo
!
! input
!    Ch3a_on - 0 if Ch3b is on, 1 if Ch3a is on
!    Sun_Earth_Distance - sun-earth distance factor
!    Rad_Ch20 - vector of Ch3b radiance (mW/m^2/str/cm^1)
!    Bt_Ch31 - vector of ch4 brightness temp. (K)
!    Bt_Ch32 - vector of ch5 brightness temp. (K)
!    Solzen - vector of solar zenith angles (degrees)
!
!  output (passed through shared memory)
!    Ref_Ch20 - Ch3b albeDO computed using standard ch4 based method
!    Emiss_Ch20 - Ch3b emissivity computed using standrad Ch3 based method
!
!  internal
!    Rad_Ch20_ems - Ch3b emission radiance (mW/m^2/str/cm^-1)
!
! Note, Ref_Ch20_Sfc is computed in ATMOS_CORR
!
! Revision History
!   January 2003 - A. Heidinger
!
! --->    AW 10/20/2014
! ch(20) %ref_toa is the pseudo solar reflectance in 3.9 channels
!  Rad_obs = Rad_sol + ( 1 - R ) Rad_ch20_ems
!  Rad_obs = (R * F_0 * mu) / PI + ( 1 - R ) Rad_ch20_ems
!  == >   R = ( PI (Rad_obs - Rad_ch20_ems )) / ( F_o * mu - PI * Rad_ch20_ems )
!     see Kaufman and Remer IEEE 1994:
!   "Detection of  Forests Using Mid-IR Reflectance: An  Application for Aerosol Studies"
!   http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=297984
!
!-----------------------------------------------------------------------
subroutine CH20_PSEUDO_REFLECTANCE(Solar_Ch20_Nu,Cos_Solzen,Rad_Ch20,Bt_Ch31,Sun_Earth_Distance,Ref_Ch20,Emiss_Ch20)

  real(kind=real4), intent(in):: Solar_Ch20_Nu
  real(kind=real4), intent(in):: Sun_Earth_Distance
  real(kind=real4), dimension(:,:), intent(in):: Cos_Solzen
  real(kind=real4), dimension(:,:), intent(in):: Rad_Ch20
  real(kind=real4), dimension(:,:), intent(in):: Bt_Ch31
  real(kind=real4), dimension(:,:), intent(out):: Ref_Ch20
  real(kind=real4), dimension(:,:), intent(out):: Emiss_Ch20
  integer:: Elem_Idx, Line_Idx, Num_Elements, Num_Lines
  real :: Rad_Ch20_Ems
  real :: Solar_Irradiance

  Num_Elements = Image%Number_Of_Elements        !make local copy of a global variable
  Num_Lines = Image%Number_Of_Lines_Per_Segment  !make local copy of a global variable

  if (( .NOT. Sensor%Chan_On_Flag_Default(20) ) .or.  &
      ( .NOT. Sensor%Chan_On_Flag_Default(31) )) then   !start Ch3a_on check
        return
  endif

  !----------------------------------------------------------------------------
  !--- standard Ref_Ch20 computation
  !---------------------------------------------------------------------------
  do Line_Idx = 1,Num_Lines
      do Elem_Idx = 1, Num_Elements

        !--- check for bad scans
        if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) ) then
          cycle
        endif

        if (Bt_Ch31(Elem_Idx,Line_Idx) > 180.0) then
           Rad_Ch20_Ems = PLANCK_RAD_FAST(20,Bt_Ch31(Elem_Idx,Line_Idx))
           Emiss_Ch20(Elem_Idx,Line_Idx) = Rad_Ch20(Elem_Idx,Line_Idx) / Rad_Ch20_Ems
        else
           Rad_Ch20_Ems = Missing_Value_Real4
           Emiss_Ch20(Elem_Idx,Line_Idx) = Missing_Value_Real4
        endif
         
        if ((Rad_Ch20_Ems>0.0).and.(Rad_Ch20(Elem_Idx,Line_Idx)>0.0)) then
           Solar_Irradiance = max (0.0, (Solar_Ch20_Nu*Cos_Solzen(Elem_Idx,Line_Idx))/(Sun_Earth_Distance**2))
           Ref_Ch20(Elem_Idx,Line_Idx) = 100.0*pi*(Rad_Ch20(Elem_Idx,Line_Idx)-Rad_Ch20_Ems) /  &
                                                  (Solar_Irradiance - pi*Rad_Ch20_Ems)
        endif

        !--- constrain values
        if (Ref_Ch20(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Ref_Ch20(Elem_Idx,Line_Idx) = max(-50.0,min(100.0,Ref_Ch20(Elem_Idx,Line_Idx)))
        endif

      enddo

   enddo

end subroutine CH20_PSEUDO_REFLECTANCE
   ! - called from process_clavrx.f90
   subroutine COMPUTE_SPATIAL_CORRELATION_ARRAYS()

      integer:: Elem_Idx
      integer:: Elem_Idx_min
      integer:: Elem_Idx_max
      integer:: Elem_Idx_width
      integer:: Elem_Idx_segment_max
      integer:: Line_Idx
      integer:: Line_Idx_min
      integer:: Line_Idx_max
      integer:: Line_Idx_width
      integer:: Line_Idx_segment_max

      Elem_Idx_segment_max = Image%Number_Of_Elements
      Line_Idx_segment_max = Line_Idx_Min_Segment + Line_Idx_Max_Segment - 1
     
      do Elem_Idx = 1, Elem_Idx_segment_max
         do Line_Idx = 1, Line_Idx_segment_max
        
            !--- compute 5x5 arrays
            Elem_Idx_min = max(1,min(Elem_Idx - 2,Elem_Idx_segment_max))
            Elem_Idx_max = max(1,min(Elem_Idx + 2,Elem_Idx_segment_max))
            Line_Idx_min = max(1,min(Line_Idx - 2,Line_Idx_segment_max))
            Line_Idx_max = max(1,min(Line_Idx + 2,Line_Idx_segment_max))
            Line_Idx_width = Line_Idx_max - Line_Idx_min + 1
            Elem_Idx_width = Elem_Idx_max - Elem_Idx_min + 1
     
      
            if ((Sensor%Chan_On_Flag_Per_Line(27,Line_Idx) ) .and. & 
               (Sensor%Chan_On_Flag_Per_Line(31,Line_Idx) )) then
            
               Covar_Ch27_Ch31_5x5(Elem_Idx,Line_Idx) = Covariance(&
                  ch(31)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
                  ch(27)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
                  Elem_Idx_width, Line_Idx_width, &
                  Bad_Pixel_Mask(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max))
               
            end if

         end do

      end do
   
    
   end subroutine COMPUTE_SPATIAL_CORRELATION_ARRAYS


!------------------------------------------------------------------------
! Routine to assign clear-sky quality flags
!------------------------------------------------------------------------
  subroutine ASSIGN_CLEAR_SKY_QUALITY_FLAGS(jmin,jmax)

    integer, intent(in):: jmin,jmax
    integer:: i, j, i1, i2, j1, j2, n, max_Mask

!-- determine size of box for spatial filter
 n = 1
                                                                                                                                                
!--- initialize
 Tsfc_Qf = 0
 Ndvi_Qf = 0
 Rsr_Qf = 0


 do j = jmin, jmin+jmax-1
                                                                                                                                                
   !--- determine y-dimensions of array to check
   j1 = max(jmin,j-n)
   j2 = min(jmax,j+n)
                                                                                                                                                
    do i = 1, Image%Number_Of_Elements

      !--- check for bad scans
      if (Bad_Pixel_Mask(i,j) ) then
        cycle
      endif
                                                                                                                                                
      !--- determine x-dimensions of array to check
      i1 = max(1,i-n)
      i2 = min(Image%Number_Of_Elements,i+n)
        
      !--- initial cloud mask based
      if (CLDMASK%cld_Mask(i,j) == sym%CLEAR) then
         Tsfc_Qf(i,j) = 3
         Ndvi_Qf(i,j) = 3
         Rsr_Qf(i,j) = 3
      endif
      if (CLDMASK%cld_Mask(i,j) == sym%PROB_CLEAR) then
         Tsfc_Qf(i,j) = 2
         Ndvi_Qf(i,j) = 2
         Rsr_Qf(i,j) = 2
      endif
      if (CLDMASK%cld_Mask(i,j) == sym%PROB_CLOUDY) then
         Tsfc_Qf(i,j) = 1
         Ndvi_Qf(i,j) = 1
         Rsr_Qf(i,j) = 1
      endif
      if (CLDMASK%cld_Mask(i,j) == sym%CLOUDY) then
         Tsfc_Qf(i,j) = 0
         Ndvi_Qf(i,j) = 0
         Rsr_Qf(i,j) = 0
      endif

      !--- assign Ndvi over water to be low quality
      if (Sfc%Land_Mask(i,j) == sym%NO) then
        Ndvi_Qf(i,j) = 0
      endif

      !--- assign Ndvi at high angles to be low quality
      if (Geo%Solzen(i,j) > 75.0) then
       Ndvi_Qf(i,j) = min(1,int(Ndvi_Qf(i,j)))
      endif

      !--- modifcations of Rsr quality
      if (Sfc%Land_Mask(i,j) .EQ. 1 ) then    !ocean only
        Rsr_Qf(i,j) = 0
      endif
      if (Geo%Solzen(i,j) > 75.0) then    !sufficient light
        Rsr_Qf(i,j) = min(1,int(Rsr_Qf(i,j)))
      endif
      if (Geo%Glintzen(i,j) < Glint_Zen_Thresh) then   !outside glint
       Rsr_Qf(i,j) = min(1,int(Rsr_Qf(i,j)))
      endif
      if (Sfc%Snow(i,j) /= sym%NO_SNOW) then    !Snow
       Rsr_Qf(i,j) = min(1,int(Rsr_Qf(i,j)))
      endif


      !--- aot 
       if (aer_flag ) then

        Aot_Qf(i,j) = 0
        
        if ((CLDMASK%Cld_Mask(i,j) == sym%CLEAR) .and.  &
            (Sfc%Land_Mask(i,j) == sym%NO) .and.  &
            (Geo%Solzen(i,j) < 70.00) .and.  &
            (Sfc%Snow(i,j) == sym%NO_SNOW)) then
          Aot_Qf(i,j) = 1
          if (Geo%Glintzen(i,j) > Glint_Zen_Thresh) then
            if (Geo%Relaz(i,j) > 90.0) then
              Aot_Qf(i,j) = 3
            else
             Aot_Qf(i,j) = 2 - 1
           endif
         endif
        endif
                                                                                                                              
        !--- assign aerosol over Snow to be of low quality
        if (Sfc%Snow(i,j) /= sym%NO_SNOW) then
         Aot_Qf(i,j) = 0
        endif

        !--- assign high quality pixels around a cloudy results as qf = 2
        max_Mask = maxval(CLDMASK%cld_Mask(i1:i2,j1:j2))
      if (max_Mask >= 2) then
        if (Aot_Qf(i,j) == 3) then
          Aot_Qf(i,j) = 2
        endif
        if (Ndvi_Qf(i,j) == 3) then
          Ndvi_Qf(i,j) = 2
        endif
        if (Tsfc_Qf(i,j) == 3) then
          Tsfc_Qf(i,j) = 2
        endif
      endif

     !--- forcing the reporting of aerosol for this condition (A. Evan)
     if (Dust_Mask(i,j) .EQ. 1 ) then
         Aot_Qf(i,j) = 3
     endif

  endif ! end of Aerosol  QF


    end do
  end do

  end subroutine ASSIGN_CLEAR_SKY_QUALITY_FLAGS

!----------------------------------------------------------------------------
! compute spatial uniformity used for each pixel
!
!  note this can handle nx2 uniformity calculations.  If n /= 2, then make
!  sure this routine is called for each pixel in the scanline.  However,
!  because we still process two scanlines, we only do this once for each
!  scanline in the pair
!
! note, this routine attempts to return the local standard deviation of area
! surrounding the pixels.  This is done using the following approximation
! sigma = (max - min) / sqrt(n)  where n is the number of pixels that comprise
! the region.  For a 2x2, n=4, for 3x2, n=6 and so forth.
!
!  Now written assuming 3x3 uniformity except at edges of segments
!---------------------------------------------------------------------------
subroutine COMPUTE_SPATIAL_UNIFORMITY(jmin,jmax)
                                                                                                           
  integer, intent(in):: jmin,jmax
  integer:: nbox

  integer:: Nx
  integer:: Ny
  integer, dimension(:,:), allocatable:: Elem_Idx_Max
  integer, dimension(:,:), allocatable:: Line_Idx_Max
  integer, dimension(:,:), allocatable:: Elem_Idx_Min
  integer, dimension(:,:), allocatable:: Line_Idx_Min

  integer:: Uni_Land_Mask_Flag_Yes
  integer:: Uni_Land_Mask_Flag_No
  integer:: Elem_Idx
  integer:: Line_Idx
 
   !----------------------------------------------------------------------
   ! determine size of input arrays
   !----------------------------------------------------------------------
   Nx = size(Nav%Lat_1b,1)
   Ny = size(Nav%Lat_1b,2)

   !----------------------------------------------------------------------
   ! allocate memory for temporary arrays
   !----------------------------------------------------------------------
   allocate(Elem_Idx_Max(Nx,Ny))
   allocate(Line_Idx_Max(Nx,Ny))
   allocate(Elem_Idx_Min(Nx,Ny))
   allocate(Line_Idx_Min(Nx,Ny))

   Elem_Idx_Max = 0
   Line_Idx_Max = 0
   Elem_Idx_Min = 0
   Line_Idx_Min = 0

   !----------------------------------------------------------------------
   ! set flag to enforce uniformity in land mask for these computations
   !----------------------------------------------------------------------
   Uni_Land_Mask_Flag_Yes = sym%YES
   Uni_Land_Mask_Flag_No = sym%NO

   !----------------------------------------------------------------------
   ! compute 3x3 metric
   !----------------------------------------------------------------------
   nbox = 1   !1 gives a 3x3 box

   !--- Bt_Ch31 
   if (Sensor%Chan_On_Flag_Default(31) ) then

    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(31)%Bt_Toa,nbox, Uni_Land_Mask_Flag_No, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Bt_Ch31_Mean_3x3,Bt_Ch31_Min_3x3, &
                                       Bt_Ch31_Max_3x3,Bt_Ch31_Std_3x3, &
                                       Elem_Idx_Max,Line_Idx_Max,Elem_Idx_Min,Line_Idx_Min)
    !--- store indices
    Elem_Idx_Max_Bt_Ch31_3x3 = Elem_Idx_Max
    Line_Idx_Max_Bt_Ch31_3x3 = Line_Idx_Max
    Elem_Idx_Min_Bt_Ch31_3x3 = Elem_Idx_Min
    Line_Idx_Min_Bt_Ch31_3x3 = Line_Idx_Min

   endif

   !--- Ref_Ch1
   if (Sensor%Chan_On_Flag_Default(1) ) then
    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(1)%Ref_Toa,nbox,  &
                                       Uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Ref_Ch1_Mean_3x3,Ref_Ch1_Min_3x3, &
                                       Ref_Ch1_Max_3x3,Ref_Ch1_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif

   !--- Ref_ChDNB_Lunar
   if (Sensor%Chan_On_Flag_Default(44) ) then
    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(44)%Ref_Lunar_Toa,nbox,  &
                                       Uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Ref_ChDNB_Lunar_Mean_3x3,Ref_ChDNB_Lunar_Min_3x3, &
                                       Ref_ChDNB_Lunar_Max_3x3,Ref_ChDNB_Lunar_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif


   !--- Ref_Ch1 clear white sky from MODIS
   if (Sensor%Chan_On_Flag_Default(1) ) then
    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(1)%Sfc_Ref_White_Sky,nbox, &
                                       Uni_Land_Mask_Flag_No, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Ref_Ch1_Sfc_White_Sky_Mean_3x3, &
                                       Temp_Pix_Array_1, &
                                       Temp_Pix_Array_2, &
                                       Temp_Pix_Array_3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
    !--- fill in wholes in white sky albedo
    where(ch(1)%Sfc_Ref_White_Sky == Missing_Value_Real4 .and. &
           Ref_Ch1_Sfc_White_Sky_Mean_3x3 /= Missing_Value_Real4) 
           ch(1)%Sfc_Ref_White_Sky = Ref_Ch1_Sfc_White_Sky_Mean_3x3
    end where
   endif

   !--- Bt_Ch20
   if (Sensor%Chan_On_Flag_Default(20) ) then
    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(20)%Bt_Toa,nbox,  & 
                                       !uni_Land_Mask_flag_yes, &
                                       Uni_Land_Mask_Flag_No, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Temp_Pix_Array_1, Temp_Pix_Array_2, &
                                       Temp_Pix_Array_3,Bt_Ch20_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif

   !--- Btd_Ch31_Ch32
   if (Sensor%Chan_On_Flag_Default(31)  .and. Sensor%Chan_On_Flag_Default(32) ) then
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       Btd_Ch31_Ch32,nbox,  &
                                       !uni_Land_Mask_flag_yes, &
                                       Uni_Land_Mask_Flag_No, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Temp_Pix_Array_1, Temp_Pix_Array_2, &
                                       Temp_Pix_Array_3, Btd_Ch31_Ch32_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)

    !----- store Btd_Ch31_Ch32 at maximum Bt_Ch31 in surrounding pixels
    line_loop: do Line_Idx=jmin, jmax - jmin + 1
      element_loop: do Elem_Idx= 1, Image%Number_Of_Elements
       if ((Elem_Idx_Max_Bt_Ch31_3x3(Elem_Idx,Line_Idx) > 0) .and. &
           (Line_Idx_Max_Bt_Ch31_3x3(Elem_Idx,Line_Idx) > 0)) then
          Btd_Ch31_Ch32_Bt_Ch31_Max_3x3(Elem_Idx,Line_Idx) =  &
          Btd_Ch31_Ch32(Elem_Idx_Max_Bt_Ch31_3x3(Elem_Idx,Line_Idx),Line_Idx_Max_Bt_Ch31_3x3(Elem_Idx,Line_Idx)) 
       endif
      end do element_loop
    end do line_loop
   endif

   !--- Btd_Ch31_Ch33
   if (Sensor%Chan_On_Flag_Default(31)  .and. Sensor%Chan_On_Flag_Default(33) ) then
     Temp_Pix_Array_1 = ch(31)%Bt_Toa - ch(33)%Bt_Toa
     where(ch(31)%Bt_Toa == Missing_Value_Real4 .or. ch(33)%Bt_Toa == Missing_Value_Real4)
             Temp_Pix_Array_1 = Missing_Value_Real4
     end where
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       Temp_Pix_Array_1,nbox,  &
                                       Uni_Land_Mask_Flag_No, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Temp_Pix_Array_1, Temp_Pix_Array_2, &
                                       Temp_Pix_Array_3, Btd_Ch31_Ch33_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)

   endif
   !--- Btd_Ch31_Ch29
   if (Sensor%Chan_On_Flag_Default(31) .and. Sensor%Chan_On_Flag_Default(29)) then
     Temp_Pix_Array_1 = ch(31)%Bt_Toa - ch(29)%Bt_Toa
     where(ch(31)%Bt_Toa == Missing_Value_Real4 .or. ch(29)%Bt_Toa == Missing_Value_Real4)
             Temp_Pix_Array_1 = Missing_Value_Real4
     end where
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       Temp_Pix_Array_1,nbox,  &
                                       Uni_Land_Mask_Flag_No, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Temp_Pix_Array_1,Temp_Pix_Array_2, &
                                       Temp_Pix_Array_3,Btd_Ch31_Ch29_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif
   !--- Btd_Ch31_Ch27
   if (Sensor%Chan_On_Flag_Default(31) .and. Sensor%Chan_On_Flag_Default(27)) then
     Temp_Pix_Array_1 = ch(31)%Bt_Toa - ch(27)%Bt_Toa
     where(ch(31)%Bt_Toa == Missing_Value_Real4 .or. ch(27)%Bt_Toa == Missing_Value_Real4)
             Temp_Pix_Array_1 = Missing_Value_Real4
     end where
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       Temp_Pix_Array_1,nbox,  &
                                       Uni_Land_Mask_Flag_No, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Temp_Pix_Array_1, Temp_Pix_Array_2, &
                                       Temp_Pix_Array_3,  Btd_Ch31_Ch27_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)

   endif
  
   !--- Ref_Ch1_Clear
   if (Sensor%Chan_On_Flag_Default(1)) then
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(1)%Ref_Toa_Clear,nbox,  &
                                       Uni_Land_Mask_Flag_No, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Ref_Ch1_Clear_Mean_3x3,Ref_Ch1_Clear_Min_3x3, &
                                       Ref_Ch1_Clear_Max_3x3, Ref_Ch1_Clear_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif

   !--- Bt_Ch27
   if (Sensor%Chan_On_Flag_Default(27)) then
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(27)%Bt_Toa,nbox,  &
                                       Uni_Land_Mask_Flag_No, &
                                       Bad_Pixel_Mask, Sfc%Land_Mask, &
                                       1,Image%Number_Of_Elements,jmin,jmax, &
                                       Temp_Pix_Array_1,Temp_Pix_Array_2, &
                                       Bt_Ch27_Max_3x3,Temp_Pix_Array_3, &
                                       Elem_Idx_max,Line_Idx_max, Elem_Idx_min, Line_Idx_min)

  endif

  !--- deallocate memory allocated in this routine
  deallocate(Elem_Idx_max)
  deallocate(Line_Idx_max)
  deallocate(Elem_Idx_min)
  deallocate(Line_Idx_min)

end subroutine COMPUTE_SPATIAL_UNIFORMITY
!----------------------------------------------------------------------
!--- Compute a mask identifying presence of oceanic glint
!--- 
!--- input and output passed through global arrays
!----------------------------------------------------------------------
subroutine COMPUTE_GLINT(Source_GLintzen, Source_Ref_Toa, Source_Ref_Std_3x3, &
                         Source_Glint_Mask)

  real, dimension(:,:), intent(in):: Source_GlintZen
  real, dimension(:,:), intent(in):: Source_Ref_Toa
  real, dimension(:,:), intent(in):: Source_Ref_Std_3x3
  integer(kind=int1),  dimension(:,:), intent(out):: Source_Glint_Mask

  !--- define local variables
  integer:: Number_Of_Lines
  integer:: Number_Of_Elements
  integer:: Elem_Idx
  integer:: Line_Idx
  real:: Refl_Thresh


  !--- alias some global sizes into local values
  Number_Of_Lines = Image%Number_Of_Lines_Per_Segment
  Number_Of_Elements = Image%Number_Of_Elements

  Source_Glint_Mask = Missing_Value_Int1

     line_loop: do Line_Idx = 1, Number_Of_Lines

     element_loop: do Elem_Idx = 1, Number_Of_Elements

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) ) then
             cycle
     endif

     !--- initialize valid pixel to no
     Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO

     !--- skip land pixels
     if ((Sfc%Land_Mask(Elem_Idx,Line_Idx) == sym%NO) .and. &
          Sfc%Snow(Elem_Idx,Line_Idx) == sym%NO_SNOW) then

       !--- turn on in geometric glint cone and sufficient Ref_Ch1
       if ((Source_Glintzen(Elem_Idx,Line_Idx) < Glint_Zen_Thresh)) then

          !--- assume to be glint if in geometric zone
          Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%YES

          if (Sensor%Chan_On_Flag_Default(31) ) then

            !--- exclude pixels colder than the freezing temperature
            if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < 273.15) then
              Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
              cycle
            endif

            !--- exclude pixels colder than the surface
            if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < ch(31)%Bt_Toa_Clear(Elem_Idx,Line_Idx) - 5.0) then
              Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
              cycle
            endif

          endif

          !-turn off if non-uniform - but not near limb
          if (Geo%Satzen(Elem_Idx,Line_Idx) < 45.0) then 
           if (Sensor%Chan_On_Flag_Default(31) ) then
            if (Bt_Ch31_Std_3x3(Elem_Idx,Line_Idx) > 1.0) then
             Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
             cycle
            endif
           endif

           if (Sensor%Chan_On_Flag_Default(1) ) then
            if (Source_Ref_Std_3x3(Elem_Idx,Line_Idx) > 2.0) then
             Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
             cycle
            endif
           endif
          endif

          !-checks on the value of ch1
          if (Sensor%Chan_On_Flag_Default(1) ) then

            !-turn off if dark
            if (Source_Ref_Toa(Elem_Idx,Line_Idx) < 5.0) then
             Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
             cycle
            endif

            !-turn off if bright
            if (Source_Glintzen(Elem_Idx,Line_Idx) > 10.0 .and. &
                Source_Glintzen(Elem_Idx,Line_Idx) < 40.0) then

               Refl_Thresh = 25.0 - Source_Glintzen(Elem_Idx,Line_Idx)/3.0

               if (Source_Ref_Toa(Elem_Idx,Line_Idx) > Refl_Thresh) then
                  Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
                  cycle
               endif

            endif

!           if (Source_Glintzen(Elem_Idx,Line_Idx) > 20.0) then
!             if (Source_Ref_Toa(Elem_Idx,Line_Idx) > 15.0) then
!                 Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
!                 cycle
!             endif
!           endif

!           if (Source_Glintzen(Elem_Idx,Line_Idx) > 10.0) then
!             if (Source_Ref_Toa(Elem_Idx,Line_Idx) > 20.0) then
!                 Source_Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
!                 cycle
!             endif
!           endif

          endif

       endif  !Glintzen check

     endif    !land check

     enddo element_loop
   enddo line_loop


 end subroutine COMPUTE_GLINT

!----------------------------------------------------------------------
! Read MODIS white sky albedoes
!----------------------------------------------------------------------
subroutine READ_MODIS_WHITE_SKY_ALBEDO(modis_alb_id,modis_alb_str,Ref_Sfc_White_Sky)

    integer(kind=4), intent(in):: modis_alb_id
    TYPE(Land_grid_description), intent(in) :: modis_alb_str
    real(kind=real4), dimension(:,:), intent(out):: Ref_Sfc_White_Sky

    CALL READ_LAND_SFC_HDF(modis_alb_id, modis_alb_str, Nav%Lat, &
                          Nav%Lon, Space_Mask, Two_Byte_Temp)
    Ref_Sfc_White_Sky = 0.1*Two_Byte_Temp

    Ref_Sfc_White_Sky = 1.10*Ref_Sfc_White_Sky   !EMPIRICAL ADJUSTMENT

    where(Two_Byte_Temp == 32767)
             Ref_Sfc_White_Sky = Missing_Value_Real4
    endwhere

    !--- modify for water
    where(Sfc%Land_Mask == sym%NO)
            Ref_Sfc_White_Sky = Ref_Sfc_White_Sky_Water
    end where

end subroutine READ_MODIS_WHITE_SKY_ALBEDO

!------------------------------------------------------------
! compute the single scater and aerosol reflectance
! this assumes that the gas is mixed in with scattering
!-----------------------------------------------------------
subroutine COMPUTE_CLEAR_SKY_SCATTER(Tau_Aer, &
                                     wo_Aer, &
                                     g_Aer, &
                                     Tau_Ray, &
                                     Tau_gas, &
                                     scatangle, &
                                     Coszen, &
                                     Cossolzen, &
                                     Cloud_Albedo_View, &
                                     Cloud_Albedo_Sun, &
                                     Ref_ss)

   real, intent(in):: Tau_Aer
   real, intent(in):: wo_Aer
   real, intent(in):: g_Aer
   real, intent(in):: Tau_Ray
   real, intent(in):: Tau_gas
   real, intent(in):: scatangle
   real, intent(in):: Coszen
   real, intent(in):: Cossolzen
   real, intent(in):: Cloud_Albedo_View
   real, intent(in):: Cloud_Albedo_Sun
   real, intent(out):: Ref_ss

   real:: Airmass
   real:: P_Aer
   real:: P_Ray
   real:: Tau_Total
   real:: Tau_Scat_Total
   real:: Trans_Total
   real:: Tau_Iso_Total
   real:: Trans_Iso_Total_View
   real:: Trans_Iso_Total_Sun
   real:: Tau_Iso_Scat_Total
   real:: mu
   real:: Pf
   real:: wo
   real:: Ref_ss_a
   real:: Ref_ss_b
   real:: Ref_ss_c

      !--- compute cosine of scattering angle
      mu = cos(scatangle*dtor)

      !-- compute Rayleigh phase function
      Airmass = 1.0/Coszen + 1.0/Cossolzen
      P_Ray = 0.75*(1.0 + mu**2)

      !--- compute total transmission
      Tau_Total = Tau_Aer + Tau_Ray + Tau_gas
      Trans_Total = exp(-Tau_Total*Airmass)

      Tau_Iso_Total = (1.0-g_Aer)*Tau_Aer + Tau_Ray + Tau_gas
      Trans_Iso_Total_View = exp(-Tau_Iso_Total/Coszen)
      Trans_Iso_Total_Sun = exp(-Tau_Iso_Total/Cossolzen)

      !--- compute total scattering optical depth
      Tau_Scat_Total = wo_Aer*Tau_Aer + Tau_Ray
      Tau_Iso_Scat_Total = wo_Aer*(1.0-g_Aer)*Tau_Aer + Tau_Ray

      !--- single scatter albedo
      wo = (wo_Aer*Tau_Aer + Tau_Ray)/ ( Tau_Total )

      !aerosol phase function (Henyey-Greenstein)
      P_Aer = (1.0 - g_Aer**2)/( (1.0 + g_Aer**2 - 2.0*g_Aer*mu)**(1.5) )

      !--- compute effective phase function
      Pf = P_aer
      if (Tau_Scat_Total > 0.0) then
        Pf = (wo_Aer*Tau_Aer*P_Aer + Tau_Ray*P_Ray)/(Tau_Scat_Total)
      endif

      !--- compute single scatter reflectance (0-100%)
      Ref_ss_a = wo*Pf/(4.0*Airmass*Coszen*Cossolzen) * (1.0 - Trans_Total )

      Ref_ss_b = (Tau_Iso_Scat_Total / (2.0*Cossolzen)) *Trans_Iso_Total_View * Cloud_Albedo_View

      Ref_ss_c = (Tau_Iso_Scat_Total / (2.0*Coszen)) *Trans_Iso_Total_Sun * Cloud_Albedo_Sun

      Ref_ss = 100.0*(Ref_ss_a + Ref_ss_b + Ref_ss_c)

end subroutine COMPUTE_CLEAR_SKY_SCATTER



!==============================================================================
!
!==============================================================================
 subroutine DETERMINE_LEVEL1B_COMPRESSION(File_1b_Original,L1b_Gzip,L1b_Bzip2)
   character(len=*), intent(in):: File_1b_Original
   logical, intent(out):: L1b_Gzip
   logical, intent(out):: L1b_Bzip2
   character(len=1020):: System_String
   character(len=7):: L1b_ext


  !--- determine if the goes data is compressed
  L1b_ext = File_1b_Original(len_trim(File_1b_Original)-2: &
                             len_trim(File_1b_Original))

  !-- determine if gzipped
  if (trim(L1b_ext) == '.gz') then
     L1b_Gzip = .TRUE.
  else
     L1b_Gzip = .FALSE.
  endif

  !--- check if bzipped
  if (trim(L1b_ext) == 'bz2') then
     L1b_Bzip2 = .TRUE.
  else
     L1b_Bzip2 = .FALSE.
  endif

  !--- uncompress
  if (L1b_Gzip ) then
     Image%Level1b_Name = File_1b_Original(1:len(trim(File_1b_Original))-3)
     System_String = "gunzip -c "//trim(Image%Level1b_Path)//trim(File_1b_Original)// &
        " > "//trim(Temporary_Data_Dir)//trim(Image%Level1b_Name)
     call SYSTEM(System_String)

     Number_of_Temporary_Files = Number_of_Temporary_Files + 1
     Temporary_File_Name(Number_of_Temporary_Files) = trim(Image%Level1b_Name)

  elseif (L1b_Bzip2 ) then
     Image%Level1b_Name = File_1b_Original(1:len(trim(File_1b_Original))-4)
     System_String = "bunzip2 -c "//trim(Image%Level1b_Path)//trim(File_1b_Original)// &
        " > "//trim(Temporary_Data_Dir)//trim(Image%Level1b_Name)
     call SYSTEM(System_String)

     Number_of_Temporary_Files = Number_of_Temporary_Files + 1
     Temporary_File_Name(Number_of_Temporary_Files) = trim(Image%Level1b_Name)

  else
     Image%Level1b_Name = trim(File_1b_Original)
  endif

   !--- make a full file name
   if (L1b_Gzip  .or. L1b_bzip2 ) then
     Image%Level1b_Full_Name = trim(Temporary_Data_Dir)//trim(Image%Level1b_Name)
   else
    Image%Level1b_Full_Name = trim(Image%Level1b_Path)//trim(Image%Level1b_Name)
   endif

 end subroutine DETERMINE_LEVEL1B_COMPRESSION

!====================================================================
!
! Attempt to fix the land classification based on observed ndvi
!
! if the ndvi is high and the land class is not land, this pixel should be land
! if the ndvi is low and the land class is land, this pixel should be water
!
!====================================================================
subroutine MODIFY_LAND_CLASS_WITH_NDVI(Line_Idx_Min,Num_Lines)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines

  integer:: Line_Idx_Max
  integer:: Elem_Idx_Max
  integer:: Elem_Idx_Min
  integer:: Num_Elements
  integer:: Line_Idx
  integer:: Elem_Idx
  real:: ndvi_temp
  real, parameter:: Ndvi_Land_Threshold = 0.25
  real, parameter:: Ndvi_Water_Threshold = -0.25
  real, parameter:: Solzen_Threshold = 60.0
  real, parameter:: Ref_Ch2_Threshold = 60.0

  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  line_loop: do Line_Idx = Line_Idx_Min, Line_Idx_Max
    element_loop: do Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

    if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) ) cycle
    if (.NOT. Sensor%Chan_On_Flag_Default(1) ) cycle
    if (.NOT. Sensor%Chan_On_Flag_Default(2) ) cycle
    if (Geo%Solzen(Elem_Idx,Line_Idx) > Solzen_Threshold) cycle
    if (index(Sensor%Sensor_Name,'MODIS') > 0) cycle                         !modis ch2 saturates, need to modify for MODIS

    Ndvi_Temp = (ch(2)%Ref_Toa(Elem_Idx,Line_Idx) - ch(1)%Ref_Toa(Elem_idx,Line_Idx)) / &
                (ch(2)%Ref_Toa(Elem_Idx,Line_Idx) + ch(1)%Ref_Toa(Elem_idx,Line_Idx)) 

    if (Ndvi_Temp > Ndvi_Land_Threshold) then
      Sfc%Land(Elem_Idx,Line_Idx) = sym%LAND
    endif

    if (Ndvi_Temp < Ndvi_Water_Threshold .and. &
        Sfc%Land(Elem_Idx,Line_Idx) == sym%LAND) then
        Sfc%Land(Elem_Idx,Line_Idx) = sym%SHALLOW_INLAND_WATER
    endif

    enddo element_loop
  enddo line_loop

end subroutine MODIFY_LAND_CLASS_WITH_NDVI
!====================================================================
! Function Name: TERM_REFL_NORM
!
! Function:
!    Renormalize reflectances to improve performance near the terminator 
! using the parameteization given by Li and Shibata 2006
!
! Description: Renormalizes reflectances in the terminator region
!   
! Calling Sequence: Refl_Chn2 = TERM_REFL_NORM(Cos_Sol_Zen,Refl_Chn2)
!   
!
! Inputs:
!   Cosine of the Solar Zenith Angle
!   Channel 2 reflectance that is normalized by cosine of solar zenith
!
! Outputs: 
!   Renormalized reflectance
!
! Dependencies:
!        none
!
! Restrictions:  None
!
! Reference: Li and Shibata 2006 Eq.(6)
!    http://journals.ametsoc.org/doi/pdf/10.1175/JAS3682.1
!
!====================================================================
 FUNCTION TERM_REFL_NORM(Cos_Sol_Zen,Reflectance)  &
          RESULT(Reflectance_Normalized)

   real(kind=real4), intent(in):: Cos_Sol_Zen
   real(kind=real4), intent(in):: Reflectance
   real(kind=real4):: Reflectance_Normalized
   real(kind=real4):: Norm_Param

   Reflectance_Normalized = Reflectance * Cos_Sol_Zen

   Norm_Param = 24.35 / (2*Cos_Sol_Zen + sqrt(498.5225*(Cos_Sol_Zen**2) + 1) )

   Reflectance_Normalized = Reflectance_Normalized*Norm_Param

 end FUNCTION TERM_REFL_NORM


!====================================================================
! Routine Name: MERGE_NWP_HIRES_ZSFC
!
! Function:
! Merge the high resolution and low resolution surface elevation
! fields into one single field
!
! Inputs: 
!    Zsfc_Hires - passed via global arrays
!    Zsfc_Nwp - passed via global arrays
!
! Outputs: 
!    Zsfc - passed via global arrays
!====================================================================
subroutine MERGE_NWP_HIRES_ZSFC(Line_Idx_Min,Num_Lines)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  integer:: Num_Elements
  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Elem_Idx_Min
  integer:: Elem_Idx_Max
  integer:: Line_Idx_Max
  integer:: Lat_NWP_Idx
  integer:: Lon_NWP_Idx


  Elem_Idx_Min = 1
  Num_Elements = Image%Number_Of_Elements
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  line_loop: do Line_Idx = Line_Idx_Min, Line_Idx_Max
    element_loop: do Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

      !--- if no, geolocation, set to missing and go to next pixel
      if (Space_Mask(Elem_Idx,Line_Idx) .EQ. 1 ) then
          Sfc%Zsfc(Elem_Idx,Line_Idx) = Missing_Value_Real4
          cycle
      endif
 
      !--- if hires value available use it, or try NWP
      if (Sfc%Zsfc_Hires(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then 

           Sfc%Zsfc(Elem_Idx,Line_Idx) = Sfc%Zsfc_Hires(Elem_Idx,Line_Idx)

      !--- if this is ocean pixel, assume zero
      elseif (Sfc%Land(Elem_Idx,Line_Idx) == sym%SHALLOW_OCEAN .or. &
              Sfc%Land(Elem_Idx,Line_Idx) == sym%MODERATE_OCEAN .or. &
              Sfc%Land(Elem_Idx,Line_Idx) == sym%DEEP_OCEAN) then

              Sfc%Zsfc(Elem_Idx,Line_Idx) = 0.0     

      !--- try NWP
      else

         Lon_Nwp_Idx = I_Nwp(Elem_Idx,Line_Idx)
         Lat_Nwp_Idx = J_Nwp(Elem_Idx,Line_Idx)

         !--- if nwp not available, assume zero
         if (Lon_Nwp_Idx > 0 .and. Lat_Nwp_Idx > 0) then
           Sfc%Zsfc(Elem_Idx,Line_Idx) = Zsfc_Nwp(Lon_Nwp_Idx, Lat_Nwp_Idx)
         else
           Sfc%Zsfc(Elem_Idx,Line_Idx) = 0.0     
         endif

      endif
          
    enddo element_loop
  enddo line_loop

end subroutine MERGE_NWP_HIRES_ZSFC

!====================================================================
! Routine Name:
!
! Function:
! Compute maximum cloud mask number in adjusted 3x3 pixels
!
! Inputs:
!    Line_Start - minimum line
!    Number_of_Lines - lines to read
!
! Outputs:
!    Adj_Pix_Cld_Mask - passed via global arrays
!====================================================================
subroutine ADJACENT_PIXEL_CLOUD_MASK(Line_Start,Number_of_Lines)

  integer (kind=int4), intent(in):: Line_Start
  integer (kind=int4), intent(in):: Number_of_Lines
  integer:: Line_Idx
  integer:: Elem_Idx
  integer:: Number_of_Elements
  integer:: i1,i2,j1,j2

  Number_of_Elements = Image%Number_Of_Elements

  CLDMASK%Adj_Pix_Cld_Mask = Missing_Value_Int1

  line_loop: do Line_Idx = Line_Start, Number_of_Lines + Line_Start - 1

    j1 = max(1,Line_Idx - 1)
    j2 = min(Number_of_Lines,Line_Idx + 1)

    element_loop: do Elem_Idx = 1, Number_of_Elements

      i1 = max(1,Elem_Idx - 1)
      i2 = min(Number_of_Elements,Elem_Idx + 1)

      CLDMASK%Adj_Pix_Cld_Mask(Elem_Idx,Line_Idx) = maxval(CLDMASK%Cld_Mask(i1:i2,j1:j2))

    enddo element_loop
  enddo line_loop

end subroutine ADJACENT_PIXEL_CLOUD_MASK

!-----------------------------------------------------------
! Determine an ACHA Success fraction
!
! Success fraction will be defined as the number of points
! with a valid Tc versus the total number of points
! processed through ACHA.  This should not points where
! no retrieval was attempted (ie clear pixels)
!
! Note this arrays uses the global variable holding the
! packed ACHA quality flags
!
! ACHA Quality Flags
! 1 - Processed (1 = yes / 0 = no)
! 2 - Valid Tc Retrieval (1 = yes, 0 = no)
! 3 - Valid ec Retrieval (1 = yes, 0 = no)
! 4 - Valid beta Retrieval (1 = yes, 0 = no)
! 5 - degraded Tc Retrieval (1 = yes, 0 = no)
! 6 - degraded ec Retrieval (1 = yes, 0 = no)
! 7 - degraded beta Retrieval (1 = yes, 0 = no)! 
!-----------------------------------------------------------
subroutine COMPUTE_ACHA_PERFORMANCE_METRICS(Processed_Count,Valid_Count,Success_Fraction)

  real(kind=real4), intent(inout):: Processed_Count
  real(kind=real4), intent(inout):: Valid_Count
  real(kind=real4), intent(out):: Success_Fraction
  integer, parameter:: Count_Min = 10
  real:: Processed_Count_Segment
  real:: Valid_Count_Segment

  Processed_Count_Segment = count(btest(ACHA%Packed_Quality_Flags,0))

  Valid_Count_Segment = count((btest(int(ACHA%Packed_Quality_Flags),1)) .and.  &
                              (btest(int(ACHA%Packed_Quality_Flags),2)) .and. &
                              (btest(int(ACHA%Packed_Quality_Flags),3)))
  
  Processed_Count = Processed_Count + Processed_Count_Segment
  Valid_Count = Valid_Count + Valid_Count_Segment

  if (Processed_Count > Count_Min) then
    Success_Fraction = Valid_Count / Processed_Count
  else
    Success_Fraction = Missing_Value_Real4 
  endif

end subroutine COMPUTE_ACHA_PERFORMANCE_METRICS

!-----------------------------------------------------------
! Determine a DCOMP success fraction
!
! Success fraction will be defined as the number of points
! with a valid COD versus the total number of points
! processed through DCOMP.  This should not points where
! no retrieval was attempted (ie clear pixels)
!
! Note this arrays uses the global variable holding the
! packed DCOMP quality flags
!
!
! DCOMP Quality Flags
!"1:Processed (0=no,1=yes) "// &
!"2:valid COD retrieval (0=yes,1=no) "// &
!"3:valid REF retrieval (0=yes,1=no) "// &
!"4:degraded COD retrieval (0=no,1=degraded) "// &
!"5:degraded REF retrieval (0=no,1=degraded) "// &
!"6:convergency (1=no,0=yes) "// &
!"7:glint (0=no,1=yes) ", &
!-----------------------------------------------------------
subroutine COMPUTE_DCOMP_PERFORMANCE_METRICS(Dcomp_Processed_Count,Dcomp_Valid_Count)

  real(kind=real4), intent(inout):: Dcomp_Processed_Count
  real(kind=real4), intent(inout):: Dcomp_Valid_Count
  real:: Processed_Count_Segment
  real:: Valid_Count_Segment
  real, parameter:: Count_Min = 10.0


  Processed_Count_Segment = count(btest(Dcomp_Quality_Flag,0))
  Valid_Count_Segment = count((.not. btest(Dcomp_Quality_Flag,1)) .and. &
                              (.not. btest(Dcomp_Quality_Flag,2)) .and. &
                              btest(Dcomp_Quality_Flag,0) )

  
  Dcomp_Processed_Count = Dcomp_Processed_Count + Processed_Count_Segment
  Dcomp_Valid_Count = Dcomp_Valid_Count + Valid_Count_Segment

  if (DCOMP_Processed_Count > Count_Min) then
    DCOMP_Success_Fraction = DCOMP_Valid_Count / DCOMP_Processed_Count
  else
    DCOMP_Success_Fraction = Missing_Value_Real4 
  endif

end subroutine COMPUTE_DCOMP_PERFORMANCE_METRICS

   !-----------------------------------------------------------
   ! Determine a Fraction of pixels with a confident cloud mask
   !
   !-----------------------------------------------------------
   subroutine COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS(Cloud_Mask_Count,Nonconfident_Cloud_Mask_Count)
      
     
      integer:: Num_Elements
      integer:: Num_Lines
      real(kind=real4), intent(inout):: Cloud_Mask_Count
      real(kind=real4), intent(inout):: Nonconfident_Cloud_Mask_Count
      integer(kind=int1), dimension(:,:), allocatable:: Mask_local
      integer(kind=int1), dimension(:,:), allocatable:: Nonconfident_Mask_local
      integer, parameter:: Count_Min = 10
      real:: Count_Segment
      real:: Nonconfident_Count_Segment
  
      Num_Elements = Image%Number_Of_Elements  !make local copy of a global variable
      Num_Lines = Image%Number_Of_Lines_Per_Segment  !make local copy of a global variable

      allocate(Mask_local(Num_Elements,Num_Lines))
      allocate(Nonconfident_Mask_local(Num_Elements,Num_Lines))

      Mask_local = 0
      Nonconfident_Mask_local = 0

      where(CLDMASK%Cld_Mask == sym%CLEAR .or. CLDMASK%Cld_Mask == sym%PROB_CLEAR &
         .or. CLDMASK%Cld_Mask == sym%PROB_CLOUDY .or. CLDMASK%Cld_Mask == sym%Cloudy)
         Mask_local = 1
      end where

      where(CLDMASK%Cld_Mask == sym%PROB_CLEAR .or. CLDMASK%Cld_Mask == sym%PROB_CLOUDY)
         Nonconfident_Mask_local = 1
      end where

      Count_Segment = sum(real(Mask_local))
      if (Count_Segment < Count_Min) then
         return
      end if
  
      deallocate ( mask_local)

      Nonconfident_Count_Segment = sum(real(Nonconfident_Mask_local))
   
      deallocate (Nonconfident_Mask_local) 
   
      Cloud_Mask_Count = Cloud_Mask_Count + Count_Segment
      Nonconfident_Cloud_Mask_Count = Nonconfident_Cloud_Mask_Count + Nonconfident_Count_Segment

      if (Cloud_Mask_Count > Count_Min) then
         Nonconfident_Cloud_Mask_Fraction = Nonconfident_Cloud_Mask_Count / Cloud_Mask_Count
      else
         Nonconfident_Cloud_Mask_Fraction = Missing_Value_Real4 
      end if
  
   end subroutine COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS


!==============================================================================
!
! remote sensing reflectance - for ocean applications 
!
! Reference: (Smyth, Tyrrell and Tarrant, 2004 GRL, 31)
!
! note, 0.63 rayleigh optical passed to this from global memory
!==============================================================================
 real elemental function REMOTE_SENSING_REFLECTANCE ( &
                                             atmos_corrected_063_reflectance, &
                                             atmos_corrected_086_reflectance, &
                                             air_mass, &
                                             solar_zenith)
 

  real, intent(in):: atmos_corrected_063_reflectance
  real, intent(in):: atmos_corrected_086_reflectance
  real, intent(in):: air_mass
  real, intent(in):: solar_zenith
  real, parameter:: solar_zenith_max_threshold = 89.0

  if (solar_zenith < solar_zenith_max_threshold) then
     remote_sensing_reflectance = (atmos_corrected_063_reflectance -  &
                                   atmos_corrected_086_reflectance) /  &
                                   exp( -0.5*Solar_Rtm%Tau_Ray(1) * air_mass)
  else
     remote_sensing_reflectance = Missing_Value_Real4
  endif

 end function REMOTE_SENSING_REFLECTANCE

!==============================================================================
!
! normalized difference vegetation index - for land applications
!
! note, missing value passed from global memory
!==============================================================================
 real elemental function NORMALIZED_DIFFERENCE_VEGETATION_INDEX ( &
                                             atmos_corrected_063_reflectance, &
                                             atmos_corrected_086_reflectance, &
                                             solar_zenith)
 

  real, intent(in):: atmos_corrected_063_reflectance
  real, intent(in):: atmos_corrected_086_reflectance
  real, intent(in):: solar_zenith
  real, parameter:: solar_zenith_max_threshold = 89.0

  if (solar_zenith < solar_zenith_max_threshold) then
    normalized_difference_vegetation_index =  &
                   (atmos_corrected_086_reflectance - atmos_corrected_063_reflectance) /  &
                   (atmos_corrected_086_reflectance + atmos_corrected_063_reflectance) 
  else
    normalized_difference_vegetation_index = Missing_Value_Real4
  endif

 end function NORMALIZED_DIFFERENCE_VEGETATION_INDEX
!==============================================================================
!
! normalized difference snow index - for land applications
!
! note, missing value passed from global memory
!==============================================================================
 real elemental function NORMALIZED_DIFFERENCE_SNOW_INDEX ( &
                                             atmos_corrected_063_reflectance, &
                                             atmos_corrected_160_reflectance, &
                                             solar_zenith)
 

  real, intent(in):: atmos_corrected_063_reflectance
  real, intent(in):: atmos_corrected_160_reflectance
  real, intent(in):: solar_zenith
  real, parameter:: solar_zenith_max_threshold = 89.0

  if (solar_zenith < solar_zenith_max_threshold) then
    normalized_difference_snow_index =  &
                   (atmos_corrected_063_reflectance - atmos_corrected_160_reflectance) /  &
                   (atmos_corrected_063_reflectance + atmos_corrected_160_reflectance) 
  else
    normalized_difference_snow_index = Missing_Value_Real4
  endif

 end function NORMALIZED_DIFFERENCE_SNOW_INDEX
!==============================================================================
! A routine that call functions to populate some simple surface parameters
!==============================================================================
 subroutine SURFACE_REMOTE_SENSING(Line_Idx_Min,Line_Idx_Max)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Line_Idx_max

  if (Sensor%Chan_On_Flag_Default(31) ) then
     call COMPUTE_TSFC(Line_Idx_Min,Line_Idx_Max)
  endif

  if (Sensor%Chan_On_Flag_Default(1)  .and. &
      Sensor%Chan_On_Flag_Default(6) ) then

       Ndsi_Toa = NORMALIZED_DIFFERENCE_SNOW_INDEX(  &
                             ch(1)%Ref_Toa, &
                             ch(6)%Ref_Toa, &
                             Geo%Solzen)

       Ndsi_Sfc = NORMALIZED_DIFFERENCE_SNOW_INDEX(  &
                             ch(1)%Ref_Sfc, &
                             ch(6)%Ref_Sfc, &
                             Geo%Solzen)
  endif

  if (Sensor%Chan_On_Flag_Default(1)  .and. &
      Sensor%Chan_On_Flag_Default(2) ) then

       Ndvi_Toa = NORMALIZED_DIFFERENCE_VEGETATION_INDEX(  &
                             ch(1)%Ref_Toa, &
                             ch(2)%Ref_Toa, &
                             Geo%Solzen)

       Ndvi_Sfc = NORMALIZED_DIFFERENCE_VEGETATION_INDEX(  &
                             ch(1)%Ref_Sfc, &
                             ch(2)%Ref_Sfc, &
                             Geo%Solzen)

       Ndvi_Sfc_White_Sky = NORMALIZED_DIFFERENCE_VEGETATION_INDEX(  &
                             ch(1)%Sfc_Ref_White_Sky, &
                             ch(2)%Sfc_Ref_White_Sky, &
                             Geo%Solzen)

       Rsr = REMOTE_SENSING_REFLECTANCE( ch(1)%Ref_Sfc, &
                                         ch(2)%Ref_Sfc, &
                                         Geo%Airmass,       &
                                         Geo%Solzen)
  endif

 end subroutine SURFACE_REMOTE_SENSING
!==============================================================================
! COMPUTE_DESERT_MASK_FOR_CLOUD_DETECTION
!
! Purpose:  Feed the location of snow-free deserts to the cloud mask
!
!
!==============================================================================
integer(kind=int1) elemental function DESERT_MASK_FOR_CLOUD_DETECTION( &
                                       Emiss_Sfc_375um,                 &
                                       Lat,                             &
                                       Snow,                            &
                                       Surface_Type)

   real(kind=real4), intent(in):: Emiss_Sfc_375um
   real(kind=real4), intent(in):: Lat
   integer(kind=int1), intent(in):: Snow
   integer(kind=int1), intent(in):: Surface_Type

   Desert_Mask_For_Cloud_Detection = 0

   if ( Snow == sym%NO_SNOW .and.  &
        Surface_Type > 0 .and.         &
        Emiss_Sfc_375um < 0.93  .and.  &
        abs(Lat) < 60.0 .and. &
        ((Surface_Type == sym%OPEN_SHRUBS_SFC) .or.  &
         (Surface_Type == sym%CLOSED_SHRUBS_SFC) .or. &
         (Surface_Type == sym%GRASSES_SFC) .or.  &
         (Surface_Type == sym%BARE_SFC)) ) then

         Desert_Mask_For_Cloud_Detection = 1

   endif

end function DESERT_MASK_FOR_CLOUD_DETECTION
!==============================================================================
! COMPUTE_CITY_MASK_FOR_CLOUD_DETECTION
!
! Purpose:  Feed the location of cities to the cloud mask
!
!
!==============================================================================
integer(kind=int1) elemental function CITY_MASK_FOR_CLOUD_DETECTION( &
                                       Rad_Lunar,                     &
                                       Surface_Type)

   real(kind=real4), intent(in):: Rad_Lunar
   integer(kind=int1), intent(in):: Surface_Type

   real, parameter:: Radiance_Lunar_City_Thresh = 2.5e-08

   City_Mask_For_Cloud_Detection = 0

   !--- use surface type information
   if (Surface_Type == sym%URBAN_SFC) then
      City_Mask_For_Cloud_Detection = 0
   endif

   !----------------------------------------------------------------------------
   !--- if lunar radiance is available, assume large values are cities or other
   !--- surface surfaces of light that we treat as cities
   !--- note, need to check if lunar radiance is available.
   !----------------------------------------------------------------------------
   if (allocated(ch(44)%Rad_Toa)) then
     if (Rad_Lunar > Radiance_Lunar_City_Thresh) then
       City_Mask_For_Cloud_Detection = 1
     endif
   endif


end function CITY_MASK_FOR_CLOUD_DETECTION

!-----------------------------------------------------------------------------
! EUMETCAST Fire detection algorithm
!
!This implements the "Current Operational Algorithm" described in:
!TOWARDS AN IMPROVED ACTIVE FIRE MONITORING PRODUCT FOR MSG SATELLITES
!Sauli Joro, Olivier Samain, Ahmet Yildirim, Leo van de Berg, Hans Joachim Lutz
!EUMETSAT, Am Kavalleriesand 31, Darmstadt, Germany
!-----------------------------------------------------------------------------
  integer elemental function FIRE_TEST (T11,T375,T11_std,T375_std,Solzen)

     real, intent(in):: T11
     real, intent(in):: T375
     real, intent(in):: T11_Std
     real, intent(in):: T375_Std
     real, intent(in):: Solzen

     real :: Bt_375um_Eumet_Fire_Thresh
     real :: Bt_Diff_Eumet_Fire_Thresh
     real :: Stddev_11um_Eumet_Fire_Thresh
     real :: Stddev_375um_Eumet_Fire_Thresh

     !---- EUMETCAST fire detection parameters
     real, parameter :: EUMETCAST_FIRE_DAY_SOLZEN_THRESH = 70.0
     real, parameter :: EUMETCAST_FIRE_NIGHT_SOLZEN_THRESH = 90.0

     real, parameter :: BT_375UM_EUMET_FIRE_DAY_THRESH = 310.0
     real, parameter :: BT_DIFF_EUMET_FIRE_DAY_THRESH = 8.0
     real, parameter :: STDDEV_11UM_EUMET_FIRE_DAY_THRESH = 1.0
     real, parameter :: STDDEV_375UM_EUMET_FIRE_DAY_THRESH = 4.0

     real, parameter :: BT_375UM_EUMET_FIRE_NIGHT_THRESH = 290.0
     real, parameter :: BT_DIFF_EUMET_FIRE_NIGHT_THRESH = 0.0
     real, parameter :: STDDEV_11UM_EUMET_FIRE_NIGHT_THRESH = 1.0
     real, parameter :: STDDEV_375UM_EUMET_FIRE_NIGHT_THRESH = 4.0

     !--- initialize
     Fire_Test = 0

     
     !--- check if all needed data are non-missing
     if (T375 /= Missing_Value_Real4 .and. &
         T375_Std /= Missing_Value_Real4 .and. &
         T11 /= Missing_Value_Real4 .and. &
         T11_Std /= Missing_Value_Real4) then

         !Day
         if (Solzen < EumetCAST_Fire_Day_Solzen_Thresh) then
            Bt_375um_Eumet_Fire_Thresh = Bt_375um_Eumet_Fire_day_Thresh
            Bt_Diff_Eumet_Fire_Thresh = Bt_Diff_Eumet_Fire_day_Thresh
            Stddev_11um_Eumet_Fire_Thresh = Stddev_11um_Eumet_Fire_Day_Thresh
            Stddev_375um_Eumet_Fire_Thresh = Stddev_375um_Eumet_Fire_Day_Thresh
         endif

         !Night
         if (Solzen > EumetCAST_Fire_Night_Solzen_Thresh) then
            Bt_375um_Eumet_Fire_Thresh = Bt_375um_Eumet_Fire_Night_Thresh
            Bt_Diff_Eumet_Fire_Thresh = Bt_Diff_Eumet_Fire_Night_Thresh
            Stddev_11um_Eumet_Fire_Thresh = Stddev_11um_Eumet_Fire_Night_Thresh
            Stddev_375um_Eumet_Fire_Thresh = Stddev_375um_Eumet_Fire_Night_Thresh
         endif

         !Twilight
         if ((Solzen >= EumetCAST_Fire_Day_Solzen_Thresh) .and. &
             (Solzen <= EumetCAST_Fire_Night_Solzen_Thresh)) then

             !linear fit day -> night
             Bt_375um_Eumet_Fire_Thresh = ((-1.0)* Solzen) + 380.0
             Bt_Diff_Eumet_Fire_Thresh = ((-0.4)* Solzen) + 36.0

             !These two don't change, but 
             Stddev_11um_Eumet_Fire_Thresh = STDDEV_11UM_EUMET_FIRE_NIGHT_THRESH
             Stddev_375um_Eumet_Fire_Thresh = STDDEV_375UM_EUMET_FIRE_NIGHT_THRESH

         endif

       ! All of these conditions need to be met
       if ((T375 > Bt_375um_Eumet_Fire_Thresh) .and. &
           ((T375 - T11) > Bt_Diff_Eumet_Fire_Thresh) .and. &
           (T375_Std > Stddev_375um_Eumet_Fire_Thresh) .and. &
           (T11_Std < Stddev_11um_Eumet_Fire_Thresh)) then
         Fire_Test = 1
       endif

     endif

  end function FIRE_TEST
!-----------------------------------------------------------
! make VIIRS look like MODIS interms of spatial sampling
! 
! this is accomplished by averaging in the along scan direction
! appropriately to acheive a pixel size that grows with
! scan angle - as would be the case with MODIS or AVHRR
!
! this does not resample the data so there is over sampling
!
!-----------------------------------------------------------
subroutine VIIRS_TO_MODIS()

  real, dimension(3):: weight
  real, parameter, dimension(3):: weight_1_0 = (/0.00,1.00,0.00/)
  real, parameter, dimension(3):: weight_1_5 = (/0.25,1.00,0.25/)
  real, parameter, dimension(3):: weight_3_0 = (/1.00,1.00,1.00/)

  integer:: Elem_Idx, Line_Idx, Chan_Idx
  integer:: Number_of_Lines, Number_of_Elements
  integer:: i1, i2


  !--- set image size
  Number_of_Elements = Image%Number_Of_Elements
  Number_of_Lines = Image%Number_Of_Lines_Read_This_Segment

  !--- loop over channels
  do Chan_Idx = 1, NChan_Clavrx

    !--- check if channel is on
    if ( .NOT. Sensor%Chan_On_Flag_Default(Chan_Idx)) cycle

    !--- solar reflectances
    if (ch(Chan_Idx)%Obs_Type == SOLAR_OBS_TYPE .or. &
        ch(Chan_Idx)%Obs_Type == MIXED_OBS_TYPE) then

      do Line_Idx = 1, Number_of_Lines
         do Elem_Idx = 2, Number_of_Elements-1
  
           i1 = Elem_Idx - 1
           i2 = Elem_Idx + 1

           !--- pick weighting scheme
           weight = weight_1_5
           if (Elem_Idx <= 640 .or. Elem_Idx >= 2561) weight = weight_3_0
           if (Elem_Idx >= 1009 .and. Elem_Idx <= 2192) weight = weight_1_0
  
           !-- apply weighting
           Temp_Pix_Array_1(Elem_Idx,Line_Idx) = sum(weight*Ch(Chan_Idx)%Ref_Toa(i1:i2,Line_Idx)) / sum(weight)

         enddo
      enddo 

      !--- copy back
      ch(Chan_Idx)%Ref_Toa = Temp_Pix_Array_1

    endif

    !--- thermal
    if (ch(Chan_Idx)%Obs_Type == THERMAL_OBS_TYPE .or. &
        ch(Chan_Idx)%Obs_Type == MIXED_OBS_TYPE) then

      do Line_Idx = 1, Number_of_Lines
         do Elem_Idx = 2, Number_of_Elements-1

           i1 = Elem_Idx - 1
           i2 = Elem_Idx + 1
  
           !--- pick weighting scheme
           weight = weight_1_5
           if (Elem_Idx <= 640 .or. Elem_Idx >= 2561) weight = weight_3_0
           if (Elem_Idx >= 1009 .and. Elem_Idx <= 2192) weight = weight_1_0
  
           !-- apply weighting
           Temp_Pix_Array_1(Elem_Idx,Line_Idx) = sum(weight*Ch(Chan_Idx)%Rad_Toa(i1:i2,Line_Idx)) / sum(weight)
  
         enddo
      enddo 

      !--- copy back
      ch(Chan_Idx)%Rad_Toa = Temp_Pix_Array_1

      !--- compute BT
      do Line_Idx = 1, Number_of_Lines
         do Elem_Idx = 2, Number_of_Elements-1
            ch(Chan_Idx)%Bt_Toa(Elem_Idx,Line_Idx) = PLANCK_TEMP_FAST(Chan_Idx,ch(Chan_Idx)%Rad_Toa(Elem_Idx,Line_Idx))
         enddo
      enddo 

    endif

  enddo
  
end subroutine VIIRS_TO_MODIS


!-----------------------------------------------------------
! end of MODULE
!-----------------------------------------------------------
end MODULE PIXEL_ROUTINES
