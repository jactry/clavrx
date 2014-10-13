! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: pixel_routines.f90 (src)
!       PIXEL_ROUTINES (program)
!
! PURPOSE: this MODULE houses routines for computing some needed pixel-level arrays
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
! COMPUTE_ERB - compute outgoing longwave radiation (OLR)
! ATMOS_CORR - perform atmospheric correction
! NORMALIZE_REFLECTANCES - divide reflectances by cosine solar zenith angle
! CH3B_ALB - compute the channel 20 reflectance
! COMPUTE_SPATIAL_UNIFORMITY - compute metrics of radiance and reflectance
!                              spatial uniformity
! SPECTRAL_CORRECT_NDVI - apply a spectral correct to Ndvi to look like NOAA14
! ASSIGN_CLEAR_SKY_QUALITY_FLAGS - assign quality flags to clear-sky products
! CONVERT_TIME - compute a time in hours based on millisecond time in leveL1b
! COMPUTE_SNOW_FIELD - based on Snow information, make a Snow field.
! COMPUTE_GLINT - derive a glint mask
! COMPUTE_GLINT_LUNAR - derive a glint mask for lunar reflectance
!
! DETERMINE_LEVEL1B_COMPRESSION
!
!--------------------------------------------------------------------------------------
MODULE PIXEL_ROUTINES
 use CONSTANTS
 use ALGORITHM_CONSTANTS
 use PIXEL_COMMON
 use NUMERICAL_ROUTINES
 use NWP_COMMON
 use PLANCK
 use LAND_SFC_PROPERTIES
 use FILE_UTILITY
 use SURFACE_PROPERTIES

 implicit none
 public:: COMPUTE_PIXEL_ARRAYS, &
          SURFACE_REMOTE_SENSING,  &
          COMPUTE_ERB,  &
          ATMOS_CORR,&
          NORMALIZE_REFLECTANCES,  &
          CH3B_ALB,  &
          COMPUTE_SPATIAL_UNIFORMITY, &
          ASSIGN_CLEAR_SKY_QUALITY_FLAGS, &
          CONVERT_TIME, &
          COMPUTE_SNOW_FIELD, &
          SET_BAD_PIXEL_MASK, &
          QUALITY_CONTROL_ANCILLARY_DATA,   &
          READ_MODIS_WHITE_SKY_ALBEDO,      &
          COMPUTE_CLEAR_SKY_SCATTER,        &
          COMPUTE_MASKED_SST,               &
          COMPUTE_GLINT,                    &
          COMPUTE_GLINT_LUNAR,              &
          QC_MODIS,                         &
          SET_CHAN_ON_FLAG,                 &
          COMPUTE_SPATIAL_CORRELATION_ARRAYS, &
          TURN_OFF_CHANNELS_BASED_ON_SENSOR, &
          DETERMINE_LEVEL1B_COMPRESSION, &
          TERM_REFL_NORM, &
          MERGE_NWP_HIRES_ZSFC, &
          ADJACENT_PIXEL_CLOUD_MASK, &
          COMPUTE_VIIRS_SST, &
          COMPUTE_ACHA_PERFORMANCE_METRICS, &
          COMPUTE_DCOMP_PERFORMANCE_METRICS, &
          COMPUTE_CLOUD_MASK_PERFORMANCE_METRICS, &
          MODIFY_LAND_CLASS_WITH_NDVI, &
          DESERT_MASK_FOR_CLOUD_DETECTION, &
          CITY_MASK_FOR_CLOUD_DETECTION

  private:: REMOTE_SENSING_REFLECTANCE, &
            NORMALIZED_DIFFERENCE_VEGETATION_INDEX, &
            NORMALIZED_DIFFERENCE_SNOW_INDEX, &
            COMPUTE_TSFC

  contains
!----------------------------------------------------------------------------
! Routine to modify Chan_On_Default_Flag based on channels that are available
!----------------------------------------------------------------------------
subroutine TURN_OFF_CHANNELS_BASED_ON_SENSOR(Avhrr_Flag,Avhrr_1_Flag, &
                                              Goes_Flag, Goes_Mop_Flag, &
                                              Goes_Sndr_Flag, &
                                              Seviri_Flag, Mtsat_Flag, &
                                              Viirs_Flag, Iff_Viirs_Flag, &
                                              Iff_Avhrr_flag, &
                                              FY2_Flag, COMS_Flag)

   integer, intent(in):: Avhrr_Flag
   integer, intent(in):: Avhrr_1_Flag
   integer, intent(in):: Goes_Flag
   integer, intent(in):: Goes_Mop_Flag
   integer, intent(in):: Goes_Sndr_Flag
   integer, intent(in):: Seviri_Flag
   integer, intent(in):: Mtsat_Flag
   integer, intent(in):: Iff_Viirs_Flag
   integer, intent(in):: Iff_Avhrr_Flag
   integer, intent(in):: FY2_Flag
   integer, intent(in):: Viirs_Flag
   integer, intent(in):: COMS_Flag

   if (Avhrr_Flag == sym%YES) then
      Chan_On_Flag_Default(3:5) = sym%NO
      Chan_On_Flag_Default(7:19) = sym%NO
      Chan_On_Flag_Default(21:30) = sym%NO
      Chan_On_Flag_Default(33:36) = sym%NO
      if (Avhrr_1_Flag == sym%YES) then
            Chan_On_Flag_Default(32) = sym%NO
      endif
      Chan_On_Flag_Default(37:42) = sym%NO
   endif 

  !GOES
  if (Goes_Flag == sym%YES) then
       Chan_On_Flag_Default(2:19) = sym%NO
       Chan_On_Flag_Default(21:26) = sym%NO
       Chan_On_Flag_Default(28:30) = sym%NO
       Chan_On_Flag_Default(34:36) = sym%NO

       if (Goes_Mop_Flag == sym%YES) then
             Chan_On_Flag_Default(32) = sym%NO
       else
             Chan_On_Flag_Default(33) = sym%NO
       endif
       Chan_On_Flag_Default(37:42) = sym%NO
  endif

  !GOES Sounder
  if (Goes_Sndr_Flag == sym%YES) then
       Chan_On_Flag_Default(2:19) = sym%NO
       Chan_On_Flag_Default(22) = sym%NO
       Chan_On_Flag_Default(26) = sym%NO
       Chan_On_Flag_Default(29) = sym%NO
  endif

  !MTSAT
  if (Mtsat_Flag == sym%YES) then
       Chan_On_Flag_Default(2:19) = sym%NO
       Chan_On_Flag_Default(21:26) = sym%NO
       Chan_On_Flag_Default(28:30) = sym%NO
       Chan_On_Flag_Default(33:36) = sym%NO
       Chan_On_Flag_Default(37:42) = sym%NO
  endif

  ! SEVIRI
  !note 3.9 channel mapped to channel 20
  if (Seviri_Flag == sym%YES) then
       Chan_On_Flag_Default(3:5) = sym%NO
       Chan_On_Flag_Default(7:19) = sym%NO
       Chan_On_Flag_Default(21:26) = sym%NO
       Chan_On_Flag_Default(34:36) = sym%NO
       Chan_On_Flag_Default(37:42) = sym%NO
  endif

  !FY2-D/E
  if (FY2_Flag == sym%YES) then
       Chan_On_Flag_Default(2:19) = sym%NO
       Chan_On_Flag_Default(21:26) = sym%NO
       Chan_On_Flag_Default(28:30) = sym%NO
       Chan_On_Flag_Default(33:36) = sym%NO
       Chan_On_Flag_Default(37:42) = sym%NO
  ENDIF

  ! VIIRS  
  IF (Viirs_Flag == sym%YES) THEN
       Chan_On_Flag_Default(10:14) = sym%NO
       Chan_On_Flag_Default(16:19) = sym%NO
       Chan_On_Flag_Default(21) = sym%NO
       Chan_On_Flag_Default(23:25) = sym%NO
       Chan_On_Flag_Default(27:28) = sym%NO
       Chan_On_Flag_Default(30) = sym%NO
       Chan_On_Flag_Default(33:36) = sym%NO
  ENDIF

  ! VIIRS + CrIS in IFF
  ! note CrIS uses 33:36 channels
  IF (Iff_Viirs_Flag == sym%YES) THEN
       Chan_On_Flag_Default(10:14) = sym%NO
       Chan_On_Flag_Default(16:19) = sym%NO
       Chan_On_Flag_Default(21) = sym%NO
       Chan_On_Flag_Default(23:25) = sym%NO
       Chan_On_Flag_Default(27:28) = sym%NO
       Chan_On_Flag_Default(30) = sym%NO
  ENDIF

   ! AVHRR + HIRS in IFF
   ! note HIRS uses 21,23:25,27:30,33:36
   if (Iff_Avhrr_Flag == sym%YES) then
      Chan_On_Flag_Default(3:5) = sym%NO
      Chan_On_Flag_Default(7:19) = sym%NO
      Chan_On_Flag_Default(26) = sym%NO
      Chan_On_Flag_Default(37:42) = sym%NO
   endif
  
  !COMS
  if (COMS_Flag == sym%YES) then
       Chan_On_Flag_Default(2:19) = sym%NO
       Chan_On_Flag_Default(21:26) = sym%NO
       Chan_On_Flag_Default(28:30) = sym%NO
       Chan_On_Flag_Default(33:36) = sym%NO
       Chan_On_Flag_Default(37:42) = sym%NO
  endif
  
end subroutine TURN_OFF_CHANNELS_BASED_ON_SENSOR
!----------------------------------------------------------------------
! set Chan_On_Flag for each to account for Ch3a/b switching on avhrr
!
! this logic allows the default values to also be used to turn off
! channels
!----------------------------------------------------------------------
  subroutine SET_CHAN_ON_FLAG(jmin,nj)

     integer (kind=int4), intent(in):: jmin
     integer (kind=int4), intent(in):: nj
     integer:: Line_Idx

     Chan_On_Flag = sym%NO
     Ch6_On_Pixel_Mask = sym%NO

     line_loop: DO Line_Idx = jmin, nj - jmin + 1

       if (Modis_Flag == sym%YES) then

          Chan_On_Flag(:,Line_Idx) = Chan_On_Flag_Default   

          !--- Steve Platnick suggested to use UNC_Ch6 (Uncertainty of channel6)
          if (Modis_Aqua_Flag == sym%YES .OR. Modis_Aqua_Mac_Flag == sym%YES ) then
             if (Chan_On_Flag_Default(6) == sym%YES) then
               if (minval(ch(6)%Unc(:,Line_Idx)) >= 15) then
                 Chan_On_Flag(6,Line_Idx) = sym%NO 
               endif
             endif  
          endif

          !--- set ch3a_on flag for AVHRR algorithms
          Ch3a_On_Avhrr(Line_Idx) = -1
          if ((Chan_On_Flag_Default(20) == sym%NO) .and. &
              (Chan_On_Flag_Default(6) == sym%YES)) then
               Ch3a_On_Avhrr(Line_Idx) = 1
          endif
          if (Chan_On_Flag_Default(20) == sym%YES) then
               Ch3a_On_Avhrr(Line_Idx) = 0
          endif
       endif

       if (Avhrr_Flag == sym%YES) then
          Chan_On_Flag(:,Line_Idx) = sym%NO  
          Chan_On_Flag(:,Line_Idx) = Chan_On_Flag_Default
          if (Ch3a_On_Avhrr(Line_Idx) == sym%YES) then
             Chan_On_Flag(6,Line_Idx) = Chan_On_Flag_Default(6)   
             Chan_On_Flag(20,Line_Idx) = sym%NO   
          endif
          if (Ch3a_On_Avhrr(Line_Idx) == sym%NO) then
             Chan_On_Flag(6,Line_Idx) = sym%NO   
             Chan_On_Flag(20,Line_Idx) = Chan_On_Flag_Default(20)   
          endif
       endif

       if (Goes_Flag == sym%YES .or. Goes_Sndr_Flag == sym%YES .or. &
           Seviri_Flag == sym%YES .or. &
           Mtsat_Flag == sym%YES .or. Viirs_Flag == sym%YES ) then
          Chan_On_Flag(:,Line_Idx) = sym%NO         
          Chan_On_Flag(:,Line_Idx) = Chan_On_Flag_Default

          !--- set ch3a_on flag for AVHRR algorithms
          Ch3a_On_Avhrr(Line_Idx) = -1
          if ((Chan_On_Flag_Default(20) == sym%NO) .and. &
              (Chan_On_Flag_Default(6) == sym%YES)) then
               Ch3a_On_Avhrr(Line_Idx) = 1
          endif
          if (Chan_On_Flag_Default(20) == sym%YES) then
               Ch3a_On_Avhrr(Line_Idx) = 0
          endif
              
       endif

       !--- set 2d mask used for channel-6 (1.6 um)
       Ch6_On_Pixel_Mask(:,Line_Idx) = Chan_On_Flag(6,Line_Idx)

     end do line_loop

  end subroutine SET_CHAN_ON_FLAG

!----------------------------------------------------------------------
! rudimentary quality check of modis
!----------------------------------------------------------------------
subroutine QC_MODIS(jmin,nj)

  integer, intent(in):: jmin,nj
  integer:: Line_Idx

  Bad_Pixel_Mask = sym%NO

  line_loop: DO Line_Idx= jmin, nj- jmin + 1
     if (maxval(ch(31)%Rad_Toa(:,Line_Idx)) < 0.0) then
        Bad_Pixel_Mask(:,Line_Idx) = sym%YES
     endif
     if (maxval(lat_1b(:,Line_Idx)) < -100.0) then
        Bad_Pixel_Mask(:,Line_Idx) = sym%YES
     endif
  enddo line_loop

end subroutine QC_MODIS
 
!======================================================================
! Check for bad pixels
!======================================================================
subroutine SET_BAD_PIXEL_MASK(Number_of_Elements,Number_of_Lines)

   integer, intent(in):: Number_of_Elements
   integer, intent(in):: Number_of_Lines
   integer:: Elem_Idx
   integer:: Line_Idx
   integer:: Number_Bad_Pixels
   integer:: Number_Bad_Pixels_Thresh
   integer:: Lon_Nwp_Idx
   integer:: Lat_Nwp_Idx

   Number_Bad_Pixels_Thresh = 0.9 * Num_Pix

!----------------------------------------------------------------------
!--- assign bad pixel mask based on scanline fatal flag
!----------------------------------------------------------------------
   !---- this is needed to ensure extra lines (beyond Num_Scans_Read) are bad
   Bad_Pixel_Mask = sym%YES

   line_loop: DO Line_Idx = 1, Number_of_Lines

      !--- initialize
      Solar_Contamination_Mask(:,Line_Idx) = sym%NO
      Bad_Pixel_Mask(:,Line_Idx) = sym%NO
      Space_Mask(:,Line_Idx) = sym%NO_SPACE

     !--- check for a bad scan
     if (Bad_Scan_Flag(Line_Idx) == sym%YES) then

       Bad_Pixel_Mask(:,Line_Idx) = sym%YES
       Space_Mask(:,Line_Idx) = sym%SPACE

     else

      !--- if not a bad scan, check pixels on this scan
      element_loop: DO Elem_Idx = 1, Number_of_Elements

        !--- missing geolocation
        if ((Lat(Elem_Idx,Line_Idx) == Missing_Value_Real4) .or.  &
            (Lon(Elem_Idx,Line_Idx) == Missing_Value_Real4)) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
           Space_Mask(Elem_Idx,Line_Idx) = sym%SPACE
        endif

        !--- NaN checks on geolocation and geometry
        if (isnan(Lat(Elem_Idx,Line_Idx)) .or.  &
            isnan(Lon(Elem_Idx,Line_Idx)) .or.  &
            isnan(Satzen(Elem_Idx,Line_Idx)) .or.  &
            isnan(Solzen(Elem_Idx,Line_Idx))) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
           Space_Mask(Elem_Idx,Line_Idx) = sym%SPACE
        endif

        !--- Satzen limit
        if (Satzen(Elem_Idx,Line_Idx) >= Satzen_Thresh_Processing) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
        endif
 
        !--- Satzen limit
        if (Satzen(Elem_Idx,Line_Idx) == Missing_Value_Real4) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
        endif

        !--- Relaz limit
        if (Relaz(Elem_Idx,Line_Idx) == Missing_Value_Real4) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
        endif

        !--- Solzen limit
        if (Solzen(Elem_Idx,Line_Idx) == Missing_Value_Real4) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
        endif

        if ((Chan_On_Flag_Default(31) == sym%YES) .and. &
          (Chan_On_Flag_Default(32) == sym%YES)) then

          !--- CALL any scan with a ridiculous pixel as bad 
          !--- this is attempt data like NOAA-16 2004 023 where
          !--- large fractions of scans are bad but not flagged as so
          if (abs(ch(31)%Bt_Toa(Elem_Idx,Line_Idx) - ch(32)%Bt_Toa(Elem_Idx,Line_Idx)) > 20.0) then
              Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
          endif

        endif

        !--- space views for considered as bad pixels
        if (Space_Mask(Elem_Idx,Line_Idx) == sym%SPACE) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
        endif

        !--- missing 11 um observations
        if (Chan_On_Flag_Default(31) == sym%YES) then

         if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < 150.0) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
         endif

         if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) > 350.0) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
         endif

         if (isnan(ch(31)%Bt_Toa(Elem_Idx,Line_Idx))) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
         endif

        endif

        !--- AVHRR/3 Wedge Filter
        !--- for AVHRR/3, consider region where Ch3a is on but seeing night
        !--- as bad data.  In these regions there is no ch3b
!       if ((AVHRR_Flag == sym%YES) .and. &
!           ((Sc_Id_WMO <= 5) .or. (Sc_Id_WMO >= 206 .and. Sc_Id_WMO <= 223)) .and. &   !AVHRR/3 only
!           (Chan_On_Flag_Default(6) == sym%YES) .and. &
!           (Chan_On_Flag_Default(20) == sym%YES) .and. &
!           (Solzen(Elem_Idx,Line_Idx) > 90.0)) then
!           if (Chan_On_Flag_Default(20) == sym%YES) then
!            if (ch(20)%Bt_Toa(Elem_Idx,Line_Idx) < 0.0) then
!               Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
!            endif
!           endif
!       endif
        

        !--- check for solar zenith angle limits
        if ((Solzen(Elem_Idx,Line_Idx) < Solzen_Min_limit) .or. (Solzen(Elem_Idx,Line_Idx) > Solzen_Max_limit)) then
             Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
        endif

        !--- check for solar contamination of nighttime data in AVHRR
        if (Chan_On_Flag_Default(1) == sym%YES) then

          if (AVHRR_Flag == sym%YES) then

            if ((Solzen(Elem_Idx,Line_Idx) > 90.0) .and. (Scatangle(Elem_Idx,Line_Idx) < 60.0)) then
              if (therm_cal_1b == sym%NO) then
!               if (Ch1_Counts(Elem_Idx,Line_Idx) - Scan_Space_Counts_Avhrr(1,Line_Idx) > 2) then 
                if (Ch1_Counts(Elem_Idx,Line_Idx) - Ch1_Dark_Count > 2) then 
                   Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
                endif
             else
               if (Ch1_Counts(Elem_Idx,Line_Idx) - Ch1_Dark_Count > 2) then 
                Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
               endif 
             endif
            endif
  
          endif

          !--- check for solar contamination of nighttime data in GOES
          if (GOES_Flag == sym%YES) then
             if ((Solzen(Elem_Idx,Line_Idx) > 90.0) .and. (Scatangle(Elem_Idx,Line_Idx) < 60.0)) then
                if (Ch1_Counts(Elem_Idx,Line_Idx) - Ch1_Dark_Count > 2) then
                   Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
                endif 
             endif
          endif

        endif

        !--- CALL any bad pixel as being space (for ancil data interp)
        if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
           Space_Mask(Elem_Idx,Line_Idx) = sym%SPACE
        endif

        !--- check for solar contamination of nighttime data in GOES
        if (GOES_Flag == sym%YES) then
          if ((Solzen(Elem_Idx,Line_Idx) > 90.0) .and.  (Scatangle(Elem_Idx,Line_Idx) < 180.0)) then
            if (Ch1_Counts(Elem_Idx,Line_Idx) - Ch1_Dark_Count > 2) then
              Solar_Contamination_Mask(Elem_Idx,Line_Idx) = sym%YES
            endif
          endif
        endif

        !--- NWP
        if (Nwp_Flag /= 0) then
            Lon_Nwp_Idx = i_Nwp(Elem_Idx,Line_Idx)
            Lat_Nwp_Idx = j_Nwp(Elem_Idx,Line_Idx)
            if (Lon_Nwp_Idx < 1 .or. Lat_Nwp_Idx < 1) then
                 Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
            else
                 if (Bad_Nwp_Mask(Lon_Nwp_Idx, Lat_Nwp_Idx) == sym%YES) Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
            endif
        endif

      end do element_loop

     endif

     !------ if 90% of pixels on a line are bad, mark the whole scan line as bad (if not already)
     if (Bad_Scan_Flag(Line_Idx) == Missing_Value_Int1) then
         Number_Bad_Pixels = sum(Bad_Pixel_Mask(:,Line_Idx),Bad_Pixel_Mask(:,Line_Idx)==sym%YES)
         if (Number_Bad_Pixels > Number_Bad_Pixels_Thresh) then
           Bad_Scan_Flag(Line_Idx) = sym%YES
         else
           Bad_Scan_Flag(Line_Idx) = sym%NO
         endif
     endif


     !---- consider any scanline with any solar contamination as a bad line
!    if (maxval(Solar_Contamination_Mask(:,Line_Idx)) == sym%YES) then
!        Bad_Scan_Flag(Line_Idx) = sym%YES
!        Bad_Pixel_Mask(:,Line_Idx) = sym%YES
!    endif
     

   end do line_loop

  !-----------------------------------------------------------------------------------
  ! if the IDPS cloud mask is to be used for product generation, make sure that
  ! pixels within the gaps are considered bad
  !-----------------------------------------------------------------------------------
  if (Cloud_Mask_Aux_Flag == sym%USE_AUX_CLOUD_MASK .and. Viirs_Flag == sym%YES) then
      where(Gap_Pixel_Mask == sym%YES)
         Bad_Pixel_Mask = sym%YES
      endwhere
  endif

  !---------------------------------------------------------------------------------------
  ! Compute the fraction of the segment covered by valid data
  !---------------------------------------------------------------------------------------
  Segment_Valid_Fraction = 1.0 - sum(float(Bad_Pixel_Mask(:,1:Number_of_Lines))) /  &
                                float(Number_of_Elements * Number_of_Lines)

end subroutine SET_BAD_PIXEL_MASK
!--------------------------------------------------------------------------
!QUALITY_CONTROL_ANCILLARY_DATA
!
! Apply some checks on ancillary data.  Call pixels bad when checks fail
!--------------------------------------------------------------------------
subroutine QUALITY_CONTROL_ANCILLARY_DATA (j1,nj)
   integer, intent(in):: j1, nj
   integer:: j2, Elem_Idx,Line_Idx

   j2 = j1 + nj - 1

   DO Line_Idx = j1,j2+j1-1

      DO Elem_Idx = 1, Num_Pix

        !--- invalid sfc type observations
        if (Sfc_Type(Elem_Idx,Line_Idx) < 0 .or. Sfc_Type(Elem_Idx,Line_Idx) > 15) then
           Bad_Pixel_Mask(Elem_Idx,Line_Idx) = sym%YES
        endif


      enddo

   enddo

end subroutine QUALITY_CONTROL_ANCILLARY_DATA


!--------------------------------------------------------------------------
! CONVERT TIME
!
! compute the utc time for each scan in fractional hours and compute the
! local time for each pixel in fractional hours
!
!--------------------------------------------------------------------------
subroutine CONVERT_TIME(j1,j2)

   integer, intent(in):: j1,j2
   integer:: i,j

   !--- loop over scans
   do j = j1,j1+j2-1

    !--- loop over pixels
    do i = 1, Num_Pix

      !--- check for a bad pixel
      if (Bad_Pixel_Mask(i,j) == sym%YES) then
        cycle
      endif

      !--- convert ms time in leveL1b to utc time in hours
      Utc_Scan_Time_Hours(j) = Scan_Time(j) / 60.0 / 60.0/ 1000.0

      !--- compute local time based on utc and longitude
      Pixel_Local_Time_Hours(i,j) = Utc_Scan_Time_Hours(j) + lon(i,j) / 15.0

      !--- constrain local time to be between 0 and 24 
      if (Pixel_Local_Time_Hours(i,j) > 24.0) then 
          Pixel_Local_Time_Hours(i,j) = Pixel_Local_Time_Hours(i,j) - 24.0
      endif

      if (Pixel_Local_Time_Hours(i,j) < 0.0) then 
          Pixel_Local_Time_Hours(i,j) = Pixel_Local_Time_Hours(i,j) + 24.0
      endif

    enddo

   enddo

end subroutine CONVERT_TIME

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

   Coszen(:,j1:j2) = cos(Satzen(:,j1:j2)*dtor)

   Seczen(:,j1:j2) = 1.0 / Coszen(:,j1:j2)

   Cossolzen(:,j1:j2) = cos(Solzen(:,j1:j2)*dtor)

   Airmass(:,j1:j2) = 1.0
   where(Solzen(:,j1:j2) /= 0.0 .and. Coszen(:,j1:j2) /= 0.0) 
    Airmass(:,j1:j2) = 1.0/Cossolzen(:,j1:j2) + 1.0/Coszen(:,j1:j2)
   endwhere

!--- other useful arrays 

   !--- channel 31 and channel 32 brightness temperature difference
   if ((Chan_On_Flag_Default(31) == sym%YES) .and. &
       (Chan_On_Flag_Default(32) == sym%YES)) then
        Btd_Ch31_Ch32(:,j1:j2) = ch(31)%Bt_Toa(:,j1:j2) - ch(32)%Bt_Toa(:,j1:j2)
   endif

   !--- SST for cloud masking - MCSST
   ! --- edited by Denis Botambekov Oct,2012
   if (Viirs_Flag == sym%NO .and. Iff_Viirs_Flag == sym%NO) then
      if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(32) == sym%YES) then
         Sst_Unmasked(:,j1:j2) = b1_day_Mask * ch(31)%Bt_Toa(:,j1:j2) + b2_day_Mask * &
                    (Btd_Ch31_Ch32(:,j1:j2)) + b3_day_Mask * (Btd_Ch31_Ch32(:,j1:j2)) * &
                    (seczen(:,j1:j2)-1.0) + b4_day_Mask + 273.15
      else
         Sst_Unmasked(:,j1:j2) = Missing_Value_Real4
      endif
   endif

   !--- channel 20 and channel 31 brightness temperature difference
    if (Chan_On_Flag_Default(20) == sym%YES .and. Chan_On_Flag_Default(31) == sym%YES) then
     Btd_Ch20_Ch31 = ch(20)%Bt_Toa - ch(31)%Bt_Toa
     where(ch(20)%Bt_Toa == Missing_Value_Real4 .or. ch(31)%Bt_Toa == Missing_Value_Real4)
       Btd_Ch20_Ch31 = Missing_Value_Real4
     endwhere
    endif

   !--- channel 20 and channel 32 brightness temperature difference
    if (Chan_On_Flag_Default(20) == sym%YES .and. Chan_On_Flag_Default(32) == sym%YES) then
     Btd_Ch20_Ch32 = ch(20)%Bt_Toa - ch(32)%Bt_Toa
     where(ch(20)%Bt_Toa == Missing_Value_Real4 .or. ch(32)%Bt_Toa == Missing_Value_Real4)
       Btd_Ch20_Ch32 = Missing_Value_Real4
     endwhere
    endif

 end subroutine COMPUTE_PIXEL_ARRAYS

!-------------------------------------------------------------------------------
!--- populate the Snow array based on all available sources of Snow data
!-------------------------------------------------------------------------------
 subroutine COMPUTE_SNOW_FIELD(j1,nj)

   integer, intent(in):: j1,nj
   integer(kind=int4):: i,j, j2, inwp,jnwp, ihires

   j2 = j1 + nj - 1

   Snow = Missing_Value_Int1

j_loop:    DO j = j1,j2
                                                                                                                                         
  i_loop:    DO i = 1, Num_Pix

      Snow(i,j) = Missing_Value_Int1

      !--- check for bad scans
      if (Bad_Pixel_Mask(i,j) == sym%YES) then
       cycle
      endif

      !--- initialize valid pixels to NO_SNOW
      Snow(i,j) = sym%NO_SNOW

      !--- save nwp indices for convenience
      inwp = i_nwp(i,j)
      jnwp = j_nwp(i,j)
     
      !--- if hires available, use it and ignore other sources
      ihires = sym%NO
      if (read_Snow_Mask == sym%READ_SNOW_HIRES .and.     &
          Failed_Hires_Snow_Mask_flag == sym%NO) then
          Snow(i,j) = Snow_Hires(i,j)
          ihires = sym%YES
      endif

      if (ihires == sym%NO) then

        !use nwp and/or Sst analysis
        if (Nwp_Flag > 0) then
          if ((inwp > 0) .and. (jnwp > 0)) then
            if (Weasd_Nwp(inwp, jnwp) > 0.1) then  !this is Snow depth
                Snow(i,j) = sym%SNOW
             endif
             if (Use_Sst_Anal == sym%NO .and. Ice_Nwp(inwp, jnwp) > 0.5) then
               Snow(i,j) = sym%SEA_ICE
             endif
           endif
         endif

         !--- GlobSnow (Land Only)
         if (Read_Snow_Mask == sym%READ_SNOW_GLOB .and.   &
          Failed_Glob_Snow_Mask_Flag == sym%NO  .and.       &
          Land(i,j) == sym%LAND  .and.       &
          lat(i,j) >= 35.0 .and. &     !note, globSnow is NH only (35-85)
          lat(i,j) <= 85.0) then

            Snow(i,j)  = Snow_Glob(i,j)

         endif  


         !-- correct for GlobSnow missing snow over Greenland and Arctic Islands
         !-- that have a barren surface type
         if ((inwp > 0) .and. (jnwp > 0)) then
          if ((Sfc_Type(i,j) == sym%Bare_Sfc .or. Sfc_Type(i,j) == sym%OPEN_SHRUBS_SFC) .and. &
             (abs(Lat(i,j)) > 60.0) .and. &
             (Weasd_Nwp(inwp, jnwp) > 0.1)) then

                Snow(i,j) = sym%SNOW

          endif
         endif

         !-- correct for GlobSnow missing snow over mountains in v1.0
         !-- this will be fixed in future versions
         if ((inwp > 0) .and. (jnwp > 0)) then
          if ((Zsfc(i,j) > 1000.0) .and.  (Weasd_Nwp(inwp, jnwp) > 0.001)) then
                Snow(i,j) = sym%SNOW
          endif
         endif

         !--- pick up small lakes in regions marked snowy in glob snow 
         if ((Land(i,j) /= sym%LAND) .and. (Snow_Glob(i,j) == sym%SNOW .or. Snow_Glob(i,j) == sym%SEA_ICE)) then
          Snow(i,j)  = sym%SEA_ICE
         endif


         !--- check for sea-ice from Sst analysis (OISST in this case)
         !--- and replace values from nwp
         !--- note, threshold 0.5 determined by matching IMS

         !--- under these conditions believe SST Analysis
         !--- and allow it to remove ice 
         if (use_Sst_anal == sym%YES .and. &
             (Land(i,j) == sym%DEEP_OCEAN .or.     &
              Land(i,j) == sym%DEEP_INLAND_WATER .or. &
              Land(i,j) == sym%MODERATE_OCEAN .or. &
              Land(i,j) == sym%SHALLOW_OCEAN))  then

           if (Sst_Anal_Cice(i,j) > 0.50) then     
             Snow(i,j) = sym%SEA_ICE
           else
             Snow(i,j) = sym%NO_SNOW
           endif

         endif

         !--- OISST does not capture sea-ice in shallow water near AntArctica
         if ((inwp > 0) .and. (jnwp > 0)) then
          if ((Land(i,j) == sym%SHALLOW_OCEAN) .and. &
              (Lat(i,j) < -60.0) .and. &
              ((Ice_Nwp(inwp,jnwp) > 0.50).or.(Weasd_Nwp(inwp,jnwp)>0.1))) then
                Snow(i,j) = sym%SEA_ICE
          endif
         endif


         !--- allow sst analysis to detect any ice (but not remove it)
         !--- this is needed since Lake Erie is cover the OISST but
         !--- classified as shallow_inland_water
         if ((Land(i,j) /= sym%LAND) .and. (Sst_Anal_Cice(i,j) > 0.50)) then     
             Snow(i,j) = sym%SEA_ICE
         endif
          
         !-- final check for consistnecy of land and snow masks
         if ((Snow(i,j) == sym%SNOW) .and. (Land(i,j) /= sym%LAND)) then
             Snow(i,j) = sym%SEA_ICE
         endif
         if ((Snow(i,j) == sym%SEA_ICE) .and. (Land(i,j) == sym%LAND)) then
             Snow(i,j) = sym%SNOW
         endif

          
       endif

       !-------------------------------------------------------
       ! use some basic tests to modify Snow
       ! only DO this when hires Snow field is unavailable
       !-------------------------------------------------------

       if (ihires == sym%NO) then

          !-- can't be Snow if warm
         if (Chan_On_Flag_Default(31) == sym%YES) then
          if (ch(31)%Bt_Toa(i,j) > 277.0) then
           Snow(i,j) = sym%NO_SNOW
          endif
         endif

         !--- some day-specific tests
         if (Solzen(i,j) < 75.0) then  ! day check

         !--- conditions where Snow is not possible
           if (Chan_On_Flag_Default(1) == sym%YES) then
            if(ch(1)%Ref_Toa(i,j) < 10.0) then
             Snow(i,j) = sym%NO_SNOW
            endif
           endif

         endif           !day check

       endif            !hires Snow check

    end do i_loop
 end do   j_loop

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
  if (Chan_On_Flag_Default(31) == sym%NO) then
     return 
  endif

  line_loop: DO Line_Idx=jmin, jmax - jmin + 1
    element_loop: DO Elem_Idx= 1, Num_Pix

      !--- check for a bad pixel pixel
      if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
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

      !--- if water and a valid SST exists, use this as the surface temperature
      !--- or if no Sst is available, use the Tsfc computed above
      if ((Land_Mask(Elem_Idx,Line_Idx) == sym%NO) .and. &
          (Snow(Elem_Idx,Line_Idx) == sym%NO_SNOW)) then

          if (Sst_Unmasked(Elem_Idx,Line_Idx) == Missing_Value_Real4) then
            Sst_Unmasked(Elem_Idx,Line_Idx) = Tsfc_Retrieved(Elem_Idx,Line_Idx)
          else
            Tsfc_Retrieved(Elem_Idx,Line_Idx) = Sst_Unmasked(Elem_Idx,Line_Idx)
          endif 

      endif  

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
   Num_Elements = Num_Pix
   Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
   Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  !---------- loop over scan lines
  line_loop: DO Line_Idx = Line_Idx_Min, Line_Idx_Max

  element_loop: DO Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

     !--- check for bad individual pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     channel_loop: do Chan_Idx = 1,42

       if (Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle

       if (Chan_Idx /= 1 .and. Chan_Idx /= 2 .and. Chan_Idx /= 5 .and. &
           Chan_Idx /= 6 .and. Chan_Idx /= 7 .and. Chan_Idx /= 42) cycle

       !--- check for valid data
       if (Chan_Idx /= 42 .and. ch(Chan_Idx)%Ref_Toa(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
       if (Chan_Idx == 42) then
         if (ch(Chan_Idx)%Ref_Lunar_Toa(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
       endif

       !--- set source angle
       Source_Zen = SolZen(Elem_Idx,Line_Idx)
       Scattering_Angle = Scatangle(Elem_Idx,Line_Idx)
       if (Chan_Idx == 42) then
            Source_Zen = LunZen(Elem_Idx,Line_Idx)
            Scattering_Angle = Scatangle_Lunar(Elem_Idx,Line_Idx)
       endif

       !--- check for appropriate illumination
       if (Source_Zen >= 90.0) cycle

       Tau_H2O = Solar_Rtm%Tau_H2O_Coef(Chan_Idx,1) + Solar_Rtm%Tau_H2O_Coef(Chan_Idx,2)*Tpw_Nwp_Pix(Elem_Idx,Line_Idx) +  &
                   Solar_Rtm%Tau_H2O_Coef(Chan_Idx,3)*(Tpw_Nwp_Pix(Elem_Idx,Line_Idx)**2)
       Tau_Gas = max(0.0,Tau_H2O) + Solar_Rtm%Tau_O3(Chan_Idx) + Solar_Rtm%Tau_O2(Chan_Idx) + Solar_Rtm%Tau_CO2(Chan_Idx) + Solar_Rtm%Tau_CH4(Chan_Idx)
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
              Albedo_View = Ch1_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif 
         if (Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch1_Snow_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(2)
         if (ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else 
              Albedo_View = Ch2_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif 
         if (Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch2_Snow_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(5)
         if (ch(5)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = ch(5)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else 
              Albedo_View = Ch6_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0 !Note there is no Ch5_Sfc_Alb_Umd
         endif 
         if (Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch5_Snow_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(6)
         if (ch(6)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = ch(6)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else 
              Albedo_View = Ch6_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif 
         if (Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch6_Snow_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(7)
         if (ch(7)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = ch(7)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else 
              Albedo_View = Ch6_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0   !Note there is no Ch7_Sfc_Alb_Umd
         endif 
         if (Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = Ch7_Snow_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(42)  !DNB - use mean of ch1 and ch2 for sfc reflectance
         if (ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4 .and. ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Albedo_View = 0.5*(ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)+ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)) / 100.0
         else 
              Albedo_View = 0.5*(Ch1_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) + Ch2_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx))) / 100.0
         endif 
         if (Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then 
              Albedo_View = 0.5*(Ch1_Snow_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx)) + Ch2_Snow_Sfc_Alb_Umd(Sfc_Type(Elem_Idx,Line_Idx))) / 100.0
         endif
       end select

       !-- assume lambertian
       Albedo_Sun = Albedo_View

       Cos_Source_Zen = Missing_Value_Real4
       if (Source_Zen >= 0.0 .and. Source_Zen < 90.0) Cos_Source_Zen = cos(Source_Zen*DTOR)

       Airmass_Factor = 1.0 / CosZen(Elem_Idx,Line_Idx) + &
                        1.0 / Cos_Source_Zen

       !--- compute atmospheric scattering
       call COMPUTE_CLEAR_SKY_SCATTER(Tau_Aer, &
                                     Wo_Aer, &
                                     G_Aer, &
                                     Tau_Ray, &
                                     Tau_Gas, &
                                     Scattering_Angle, &
                                     Coszen(Elem_Idx,Line_Idx), &
                                     Cos_Source_Zen, &
                                     Albedo_View, &
                                     Albedo_Sun, &
                                     Ref_ss)

       !--- compute total transmission for combining terms
       Tau_Total = Tau_Aer + Tau_Ray + Tau_Gas
       Trans_Total = exp(-Tau_Total*Airmass_Factor)

       if (Chan_Idx /= 42) then

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

    !--- thermal atmospheric correction - use rtm parameters
   do Chan_Idx = 20, 36

      if (Chan_Idx == 26) cycle

      if (Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle

      ch(Chan_Idx)%Rad_Sfc(Elem_Idx,Line_Idx) = max(0.0,(ch(Chan_Idx)%Rad_Toa(Elem_Idx,Line_Idx) -  &
                                                         ch(Chan_Idx)%Rad_Atm(Elem_Idx,Line_Idx))/ &
                                                         ch(Chan_Idx)%Trans_Atm(Elem_Idx,Line_Idx))
   enddo

   !--- channel 20 reflectance correction
   if (Chan_On_Flag_Default(20) == sym%YES .and.   &
       Chan_On_Flag_Default(31) == sym%YES) then

       if (ch(31)%Bt_Sfc(Elem_Idx,Line_Idx) > 180.0) then
          Rad_Ch20_Ems_Sfc(Elem_Idx,Line_Idx) = PLANCK_RAD_FAST(20,ch(31)%Bt_Sfc(Elem_Idx,Line_Idx))
          Ems_Ch20_Sfc(Elem_Idx,Line_Idx) = ch(20)%Rad_Sfc(Elem_Idx,Line_Idx) / Rad_Ch20_Ems_Sfc(Elem_Idx,Line_Idx)
       else
         Rad_Ch20_Ems_Sfc(Elem_Idx,Line_Idx) = Missing_Value_Real4
         Ems_Ch20_Sfc(Elem_Idx,Line_Idx) = Missing_Value_Real4
       endif
       if ((Solzen(Elem_Idx,Line_Idx)<90.0).and.(Rad_Ch20_Ems_Sfc(Elem_Idx,Line_Idx)>0.0).and.(ch(20)%Rad_Sfc(Elem_Idx,Line_Idx)>0.0)) then
         ch(20)%Ref_Sfc(Elem_Idx,Line_Idx) = 100.0*pi*(ch(20)%Rad_Sfc(Elem_Idx,Line_Idx)-Rad_Ch20_Ems_Sfc(Elem_Idx,Line_Idx)) /  &
                  ((Solar_Ch20_Nu*Cossolzen(Elem_Idx,Line_Idx))/(Sun_Earth_Distance**2) - &
                  pi*Rad_Ch20_Ems_Sfc(Elem_Idx,Line_Idx) )
       endif

       !----- during the night, derive Ref_Ch20_Sfc from Ems_Ch20_Sfc
       if ((Solzen(Elem_Idx,Line_Idx)>=90.0).and.(Ems_Ch20_Sfc(Elem_Idx,Line_Idx)/=Missing_Value_Real4)) then
         ch(20)%Ref_Sfc(Elem_Idx,Line_Idx) = 100.0*(1.0-Ems_Ch20_Sfc(Elem_Idx,Line_Idx))
       endif

       !--- constrain values
       if (ch(20)%Ref_Sfc(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
        ch(20)%Ref_Sfc(Elem_Idx,Line_Idx) = max(-50.0,min(100.0,ch(20)%Ref_Sfc(Elem_Idx,Line_Idx)))
       endif

    endif


   end do element_loop

 end do line_loop

end subroutine ATMOS_CORR
!------------------------------------------------------------------
! derive the earth radiation budget parameters 
!------------------------------------------------------------------
 subroutine COMPUTE_ERB(Line_Idx_Min,Num_Lines)

  integer, intent(in):: Line_Idx_Min
  integer, intent(in):: Num_Lines
  integer:: Line_Idx_Max

  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  !--- OLR------------------------------------------------------------------
  ! this algorithm assumes that the window channels on avhrr are best used
  ! to derive the Olr component in the 8.5 to 12 micron window region.
  ! A seperate regression is used to predict the broadband Olr from the
  ! window-component of the Olr
  !--------------------------------------------------------------------------
   Olr(:,Line_Idx_Min:Line_Idx_Max) = Missing_Value_Real4
   Olr_Qf(:,Line_Idx_Min:Line_Idx_Max) = 0

   !--- check for required channels
   if (Chan_On_Flag_Default(31) == sym%NO) return
   if (Chan_On_Flag_Default(32) == sym%NO) return

   where(ch(31)%Bt_Toa(:,Line_Idx_Min:Line_Idx_Max) > 100.0 .and. &
         ch(31)%Bt_Toa(:,Line_Idx_Min:Line_Idx_Max) < 350.0 .and. &
         Bad_Pixel_Mask == sym%NO)

    !-- this makes the window eff flux temp (K) (win_Olr = sigma*win_temp^4)
    Olr(:,Line_Idx_Min:Line_Idx_Max) = win_0 + win_1*ch(32)%Bt_Toa(:,Line_Idx_Min:Line_Idx_Max) + &
                   win_2*ch(32)%Bt_Toa(:,Line_Idx_Min:Line_Idx_Max)*seczen(:,Line_Idx_Min:Line_Idx_Max) + &
                   win_3*(Btd_Ch31_Ch32(:,Line_Idx_Min:Line_Idx_Max))*seczen(:,Line_Idx_Min:Line_Idx_Max)

    !--- this make the Olr temp (K)
    Olr(:,Line_Idx_Min:Line_Idx_Max) = Olr_0 + Olr_1*Olr(:,Line_Idx_Min:Line_Idx_Max) + &
                       Olr_2*(Btd_Ch31_Ch32(:,Line_Idx_Min:Line_Idx_Max)) + &
                       Olr_3*(Btd_Ch31_Ch32(:,Line_Idx_Min:Line_Idx_Max))*seczen(:,Line_Idx_Min:Line_Idx_Max) + &
                       Olr_4*Olr(:,Line_Idx_Min:Line_Idx_Max)*seczen(:,Line_Idx_Min:Line_Idx_Max)

    !--- this make the Olr  flux (W/m^2)
    Olr(:,Line_Idx_Min:Line_Idx_Max) = stefan_boltzmann_constant * (Olr(:,Line_Idx_Min:Line_Idx_Max)**4)
    Olr_Qf(:,Line_Idx_Min:Line_Idx_Max) = 3
   endwhere


 end subroutine COMPUTE_ERB 

!======================================================================
! Normalize the reflectances by the solar zenith angle cosine
!======================================================================
 subroutine NORMALIZE_REFLECTANCES(Sun_Earth_Distance,j1,nj)
  real(kind=real4), intent(in):: Sun_Earth_Distance
  integer, intent(in):: j1
  integer, intent(in):: nj
  integer:: i,j,j2, Chan_Idx

  j2 = j1 + nj - 1

  DO j = j1,j2

     DO i = 1, Num_Pix

      if (Bad_Pixel_Mask(i,j) == sym%NO .and. Cossolzen(i,j) > 0.0) then

       do Chan_Idx = 1,19
          if (Chan_On_Flag_Default(Chan_Idx) == sym%YES) ch(Chan_Idx)%Ref_Toa(i,j) = ch(Chan_Idx)%Ref_Toa(i,j) / Cossolzen(i,j)
       enddo
       if (Chan_On_Flag_Default(26) == sym%YES) ch(26)%Ref_Toa(i,j) = ch(26)%Ref_Toa(i,j) / Cossolzen(i,j)

       if ((Chan_On_Flag_Default(1) == sym%YES)) then
         if (Ref_Ch1_Dark_Composite(i,j) /= Missing_Value_Real4) then
           Ref_Ch1_Dark_Composite(i,j) = (Ref_Ch1_Dark_Composite(i,j) / Cossolzen(i,j)) 
         endif
       endif

       if (Chan_On_Flag_Default(6) == sym%YES) then
           if (ch(6)%Ref_Toa(i,j) < 0) ch(6)%Ref_Toa(i,j) = Missing_Value_Real4
           if (AVHRR_Flag == sym%YES .and. Ch3a_On_Avhrr(j) /= sym%YES) ch(6)%Ref_Toa(i,j) = Missing_Value_Real4
       endif

       if (Avhrr_Flag == sym%YES ) then
        if (Chan_On_Flag_Default(1) == sym%YES) ch(1)%Ref_Toa(i,j) = ch(1)%Ref_Toa(i,j) * (Sun_Earth_Distance**2)
        if (Chan_On_Flag_Default(2) == sym%YES) ch(2)%Ref_Toa(i,j) = ch(2)%Ref_Toa(i,j) * (Sun_Earth_Distance**2)
        if (Chan_On_Flag_Default(6) == sym%YES .and. Ch3a_On_Avhrr(j)==sym%YES) then
              ch(6)%Ref_Toa(i,j) = ch(6)%Ref_Toa(i,j) * (Sun_Earth_Distance**2)
        endif
       endif

       if (Goes_Flag == sym%YES) then
        if (Chan_On_Flag_Default(1) == sym%YES) ch(1)%Ref_Toa(i,j) = ch(1)%Ref_Toa(i,j) * (Sun_Earth_Distance**2)
        if ((Chan_On_Flag_Default(1) == sym%YES)) then
           if (Ref_Ch1_Dark_Composite(i,j) /= Missing_Value_Real4) then
                   Ref_Ch1_Dark_Composite(i,j) = Ref_Ch1_Dark_Composite(i,j) * (Sun_Earth_Distance**2)
           endif
        endif
       endif


       !--- in terminator region, renormalize Channel 1 (maybe extend to all?)
       if (Solzen(i,j) > TERMINATOR_REFLECTANCE_SOL_ZEN_THRESH) then
          if (Chan_On_Flag_Default(1) == sym%YES)  &
              ch(1)%Ref_Toa(i,j) = TERM_REFL_NORM(Cossolzen(i,j),ch(1)%Ref_Toa(i,j))
       endif

       !---- make unnormalized relflectances for AWIPS display
        if (Chan_On_Flag_Default(1) == sym%YES)  then
          if (ch(1)%Ref_Toa(i,j) /= Missing_Value_Real4) then
            ch(1)%Ref_Toa_Unnorm(i,j) = ch(1)%Ref_Toa(i,j) * Cossolzen(i,j)
          endif
        endif

        if (Chan_On_Flag_Default(2) == sym%YES) then
          if (ch(2)%Ref_Toa(i,j) /= Missing_Value_Real4) then
            ch(2)%Ref_Toa_Unnorm(i,j) = ch(2)%Ref_Toa(i,j) * Cossolzen(i,j)
          endif
        endif

        if (Chan_On_Flag_Default(6) == sym%YES) then
          if (ch(6)%Ref_Toa(i,j) /= Missing_Value_Real4) then
            ch(6)%Ref_Toa_Unnorm(i,j) = ch(6)%Ref_Toa(i,j) * Cossolzen(i,j)
          endif
        endif

      else

       !--- set to missing
       do Chan_Idx = 1,19
          if (Chan_On_Flag_Default(Chan_Idx) == sym%YES) ch(Chan_Idx)%Ref_Toa(i,j) = Missing_Value_Real4
       enddo
       if (Chan_On_Flag_Default(26) == sym%YES) ch(26)%Ref_Toa(i,j) = Missing_Value_Real4

       !--- un-normalized refs for AWIPS
       if (Chan_On_Flag_Default(1) == sym%YES) ch(1)%Ref_Toa_Unnorm(i,j) = Missing_Value_Real4
       if (Chan_On_Flag_Default(2) == sym%YES) ch(2)%Ref_Toa_Unnorm(i,j) = Missing_Value_Real4
       if (Chan_On_Flag_Default(6) == sym%YES) ch(6)%Ref_Toa_Unnorm(i,j) = Missing_Value_Real4


      endif

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
!    Ems_Ch20 - Ch3b emissivity computed using standrad Ch3 based method
!
!  internal
!    Rad_Ch20_ems - Ch3b emission radiance (mW/m^2/str/cm^-1)
!
! Note, Ref_Ch20_Sfc is computed in ATMOS_CORR
!
! Revision History
!   January 2003 - A. Heidinger
!-----------------------------------------------------------------------
subroutine CH3B_ALB(Sun_Earth_Distance,j1,j2)
  real(kind=real4), intent(in):: Sun_Earth_Distance
  integer, intent(in):: j1, j2
  integer:: i,j

  DO j = j1,j1+j2-1

    !---- compute channel 3b albedo
    if ((Chan_On_Flag_Default(20) == sym%YES) .and.  &
        (Chan_On_Flag_Default(31) == sym%YES)) then   !start Ch3a_on check

      !----------------------------------------------------------------------------
      !--- standard Ref_Ch20 computation
      !---------------------------------------------------------------------------
      DO i = 1, Num_Pix


      !--- check for bad scans
        if (Bad_Pixel_Mask(i,j) == sym%YES) then
          cycle
        endif

        if (ch(31)%Bt_Toa(i,j) > 180.0) then
           Rad_Ch20_ems(i,j) = PLANCK_RAD_FAST(20,ch(31)%Bt_Toa(i,j))
           Ems_Ch20(i,j) = ch(20)%Rad_Toa(i,j) / Rad_Ch20_ems(i,j)
        else
           Rad_Ch20_ems(i,j) = Missing_Value_Real4
           Ems_Ch20(i,j) = Missing_Value_Real4
        endif

        if ((Solzen(i,j)<90.0).and.(Rad_Ch20_ems(i,j)>0.0).and.(ch(20)%Rad_Toa(i,j)>0.0)) then
           ch(20)%Ref_Toa(i,j) = 100.0*pi*(ch(20)%Rad_Toa(i,j)-Rad_Ch20_ems(i,j)) /  &
                        ((Solar_Ch20_nu*Cossolzen(i,j))/(Sun_Earth_Distance**2) - &
                          pi*Rad_Ch20_ems(i,j) )
        endif

         !---- at night, compute albeDO as 1 - emissivity
        if ((Solzen(i,j)>=90.0).and.(Ems_Ch20(i,j)/=Missing_Value_Real4)) then
          ch(20)%Ref_Toa(i,j) = 100.0*(1.0-Ems_Ch20(i,j))
        endif

        !--- constrain values
        if (ch(20)%Ref_Toa(i,j) /= Missing_Value_Real4) then
              ch(20)%Ref_Toa(i,j) = max(-50.0,min(100.0,ch(20)%Ref_Toa(i,j)))
        endif

      enddo

     endif   !--end of Ch3a_on check

   end do

!--- constrain values

end subroutine CH3B_ALB

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

    Elem_Idx_segment_max = Num_Pix
    Line_Idx_segment_max = Line_Idx_Min_Segment + Line_Idx_Max_Segment - 1

    DO Elem_Idx = 1, Elem_Idx_segment_max
       DO Line_Idx = 1, Line_Idx_segment_max

        !--- compute 5x5 arrays
        Elem_Idx_min = max(1,min(Elem_Idx - 2,Elem_Idx_segment_max))
        Elem_Idx_max = max(1,min(Elem_Idx + 2,Elem_Idx_segment_max))
        Line_Idx_min = max(1,min(Line_Idx - 2,Line_Idx_segment_max))
        Line_Idx_max = max(1,min(Line_Idx + 2,Line_Idx_segment_max))
        Line_Idx_width = Line_Idx_max - Line_Idx_min + 1
        Elem_Idx_width = Elem_Idx_max - Elem_Idx_min + 1

        if ((Chan_on_flag(27,Line_Idx) == sym%YES) .and. & 
            (Chan_on_flag(31,Line_Idx) == sym%YES)) then

            Covar_Ch27_Ch31_5x5(Elem_Idx,Line_Idx) = Covariance(&
               ch(31)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               ch(27)%Bt_Toa(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max), &
               Elem_Idx_width, Line_Idx_width, &
               Bad_Pixel_Mask(Elem_Idx_min:Elem_Idx_max,Line_Idx_min:Line_Idx_max))
        endif

      enddo

    enddo
    
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


 DO j = jmin, jmin+jmax-1
                                                                                                                                                
   !--- determine y-dimensions of array to check
   j1 = max(jmin,j-n)
   j2 = min(jmax,j+n)
                                                                                                                                                
    DO i = 1, Num_Pix

      !--- check for bad scans
      if (Bad_Pixel_Mask(i,j) == sym%YES) then
        cycle
      endif
                                                                                                                                                
      !--- determine x-dimensions of array to check
      i1 = max(1,i-n)
      i2 = min(Num_Pix,i+n)
        
      !--- initial cloud mask based
      if (cld_Mask(i,j) == sym%CLEAR) then
         Tsfc_Qf(i,j) = 3
         Ndvi_Qf(i,j) = 3
         Rsr_Qf(i,j) = 3
      endif
      if (cld_Mask(i,j) == sym%PROB_CLEAR) then
         Tsfc_Qf(i,j) = 2
         Ndvi_Qf(i,j) = 2
         Rsr_Qf(i,j) = 2
      endif
      if (cld_Mask(i,j) == sym%PROB_CLOUDY) then
         Tsfc_Qf(i,j) = 1
         Ndvi_Qf(i,j) = 1
         Rsr_Qf(i,j) = 1
      endif
      if (cld_Mask(i,j) == sym%CLOUDY) then
         Tsfc_Qf(i,j) = 0
         Ndvi_Qf(i,j) = 0
         Rsr_Qf(i,j) = 0
      endif

      !--- assign Ndvi over water to be low quality
      if (Land_Mask(i,j) == sym%NO) then
        Ndvi_Qf(i,j) = 0
      endif

      !--- assign Ndvi at high angles to be low quality
      if (Solzen(i,j) > 75.0) then
       Ndvi_Qf(i,j) = min(1,int(Ndvi_Qf(i,j)))
      endif

      !--- modifcations of Rsr quality
      if (Land_Mask(i,j) == sym%YES) then    !ocean only
        Rsr_Qf(i,j) = 0
      endif
      if (Solzen(i,j) > 75.0) then    !sufficient light
        Rsr_Qf(i,j) = min(1,int(Rsr_Qf(i,j)))
      endif
      if (Glintzen(i,j) < Glint_Zen_Thresh) then   !outside glint
       Rsr_Qf(i,j) = min(1,int(Rsr_Qf(i,j)))
      endif
      if (Snow(i,j) /= sym%NO_SNOW) then    !Snow
       Rsr_Qf(i,j) = min(1,int(Rsr_Qf(i,j)))
      endif


      !--- aot 
       if (aer_flag == sym%YES) then

        Aot_Qf(i,j) = 0
        
        if ((Cld_Mask(i,j) == sym%CLEAR) .and.  &
            (Land_Mask(i,j) == sym%NO) .and.  &
            (Solzen(i,j) < 70.00) .and.  &
            (Snow(i,j) == sym%NO_SNOW)) then
          Aot_Qf(i,j) = 1
          if (Glintzen(i,j) > Glint_Zen_Thresh) then
            if (Relaz(i,j) > 90.0) then
              Aot_Qf(i,j) = 3
            else
             Aot_Qf(i,j) = 2 - 1
           endif
         endif
        endif
                                                                                                                              
        !--- assign aerosol over Snow to be of low quality
        if (Snow(i,j) /= sym%NO_SNOW) then
         Aot_Qf(i,j) = 0
        endif

        !--- assign high quality pixels around a cloudy results as qf = 2
        max_Mask = maxval(cld_Mask(i1:i2,j1:j2))
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
     if (Dust_Mask(i,j) == sym%YES) then
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
!  because we still process two scanlines, we only DO this once for each
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
   Nx = size(Lat_1b,1)
   Ny = size(Lat_1b,2)

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
   if (Chan_On_Flag_Default(31) == sym%YES) then

    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(31)%Bt_Toa,nbox, uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
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
   if (Chan_On_Flag_Default(1) == sym%YES) then
    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(1)%Ref_Toa,nbox,  &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Ref_Ch1_Mean_3x3,Ref_Ch1_Min_3x3, &
                                       Ref_Ch1_Max_3x3,Ref_Ch1_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif

   !--- Ref_ChDNB_Lunar
   if (Chan_On_Flag_Default(42) == sym%YES) then
    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(42)%Ref_Lunar_Toa,nbox,  &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Ref_ChDNB_Lunar_Mean_3x3,Ref_ChDNB_Lunar_Min_3x3, &
                                       Ref_ChDNB_Lunar_Max_3x3,Ref_ChDNB_Lunar_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif


   !--- Ems_Ch20
   if (Chan_On_Flag_Default(20) == sym%YES) then
    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       Ems_Ch20,nbox,  &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Ems_Ch20_Mean_3x3,Ems_Ch20_Min_3x3, &
                                       Ems_Ch20_Max_3x3,Ems_Ch20_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif

   !--- Zsfc
   CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( & 
                                       Zsfc,nbox,  &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       zsfc_Mean_3x3,zsfc_Min_3x3, &
                                       zsfc_Max_3x3,zsfc_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)

   !--- Ref_Ch1 clear white sky from MODIS
   if (Chan_On_Flag_Default(1) == sym%YES) then
    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(1)%Sfc_Ref_White_Sky,nbox, &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Ref_Ch1_Sfc_White_Sky_Mean_3x3, &
                                       Ref_Ch1_Sfc_White_Sky_Min_3x3, &
                                       Ref_Ch1_Sfc_White_Sky_Max_3x3, &
                                       Ref_Ch1_Sfc_White_Sky_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
    !--- fill in wholes in white sky albedo
    where(ch(1)%Sfc_Ref_White_Sky == Missing_Value_Real4 .and. &
           Ref_Ch1_Sfc_White_Sky_Mean_3x3 /= Missing_Value_Real4) 
           ch(1)%Sfc_Ref_White_Sky = Ref_Ch1_Sfc_White_Sky_Mean_3x3
    end where
   endif

   !--- Bt_Ch20
   if (Chan_On_Flag_Default(20) == sym%YES) then
    CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(20)%Bt_Toa,nbox,  & 
                                       !uni_Land_Mask_flag_yes, &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Bt_Ch20_Mean_3x3,Bt_Ch20_Min_3x3, &
                                       Bt_Ch20_Max_3x3,Bt_Ch20_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif

   !--- Btd_Ch31_Ch32
   if (Chan_On_Flag_Default(31)==sym%YES .and. Chan_On_Flag_Default(32)==sym%YES) then
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       Btd_Ch31_Ch32,nbox,  &
                                       !uni_Land_Mask_flag_yes, &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Btd_Ch31_Ch32_Mean_3x3,Btd_Ch31_Ch32_Min_3x3, &
                                       Btd_Ch31_Ch32_Max_3x3,Btd_Ch31_Ch32_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)

    !----- store Btd_Ch31_Ch32 at maximum Bt_Ch31 in surrounding pixels
    line_loop: DO Line_Idx=jmin, jmax - jmin + 1
      element_loop: DO Elem_Idx= 1, Num_Pix
       if ((Elem_Idx_Max_Bt_Ch31_3x3(Elem_Idx,Line_Idx) > 0) .and. &
           (Line_Idx_Max_Bt_Ch31_3x3(Elem_Idx,Line_Idx) > 0)) then
          Btd_Ch31_Ch32_Bt_Ch31_Max_3x3(Elem_Idx,Line_Idx) =  &
          Btd_Ch31_Ch32(Elem_Idx_Max_Bt_Ch31_3x3(Elem_Idx,Line_Idx),Line_Idx_Max_Bt_Ch31_3x3(Elem_Idx,Line_Idx)) 
       endif
      end do element_loop
    end do line_loop
   endif

   !--- Btd_Ch31_Ch33
   if (Chan_On_Flag_Default(31)==sym%YES .and. Chan_On_Flag_Default(33)==sym%YES) then
     Temp_Pix_Array = ch(31)%Bt_Toa - ch(33)%Bt_Toa
     where(ch(31)%Bt_Toa == Missing_Value_Real4 .or. ch(33)%Bt_Toa == Missing_Value_Real4)
             Temp_Pix_Array = Missing_Value_Real4
     end where
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       Temp_Pix_Array,nbox,  &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Btd_Ch31_Ch33_Mean_3x3,Btd_Ch31_Ch33_Min_3x3, &
                                       Btd_Ch31_Ch33_Max_3x3,Btd_Ch31_Ch33_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)

   endif
   !--- Btd_Ch31_Ch29
   if (Chan_On_Flag_Default(31)==sym%YES .and. Chan_On_Flag_Default(29)==sym%YES) then
     Temp_Pix_Array = ch(31)%Bt_Toa - ch(29)%Bt_Toa
     where(ch(31)%Bt_Toa == Missing_Value_Real4 .or. ch(29)%Bt_Toa == Missing_Value_Real4)
             Temp_Pix_Array = Missing_Value_Real4
     end where
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       Temp_Pix_Array,nbox,  &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Btd_Ch31_Ch29_Mean_3x3,Btd_Ch31_Ch29_Min_3x3, &
                                       Btd_Ch31_Ch29_Max_3x3,Btd_Ch31_Ch29_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif
   !--- Btd_Ch31_Ch27
   if (Chan_On_Flag_Default(31)==sym%YES .and. Chan_On_Flag_Default(27)==sym%YES) then
     Temp_Pix_Array = ch(31)%Bt_Toa - ch(27)%Bt_Toa
     where(ch(31)%Bt_Toa == Missing_Value_Real4 .or. ch(27)%Bt_Toa == Missing_Value_Real4)
             Temp_Pix_Array = Missing_Value_Real4
     end where
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       Temp_Pix_Array,nbox,  &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Btd_Ch31_Ch27_Mean_3x3,Btd_Ch31_Ch27_Min_3x3, &
                                       Btd_Ch31_Ch27_Max_3x3,Btd_Ch31_Ch27_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)

   endif
  
   !--- Ref_Ch1_Clear
   if (Chan_On_Flag_Default(1)==sym%YES) then
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(1)%Ref_Toa_Clear,nbox,  &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Ref_Ch1_Clear_Mean_3x3,Ref_Ch1_Clear_Min_3x3, &
                                       Ref_Ch1_Clear_Max_3x3, Ref_Ch1_Clear_Std_3x3, &
                                       Elem_Idx_max, Line_Idx_max, Elem_Idx_min, Line_Idx_min)
   endif

   !--- Bt_Ch27
   if (Chan_On_Flag_Default(27)==sym%YES) then
     CALL COMPUTE_SPATIAL_UNIFORMITY_NxN_WITH_INDICES( &
                                       ch(27)%Bt_Toa,nbox,  &
                                       uni_Land_Mask_flag_no, &
                                       Bad_Pixel_Mask, Land_Mask, &
                                       1,Num_Pix,jmin,jmax, &
                                       Bt_Ch27_Mean_3x3,Bt_Ch27_Min_3x3, &
                                       Bt_Ch27_Max_3x3,Bt_Ch27_Std_3x3, &
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
subroutine COMPUTE_GLINT()

  !--- define local variables
  integer:: Number_Of_Lines
  integer:: Number_Of_Elements
  integer:: Elem_Idx
  integer:: Line_Idx


  !--- alias some global sizes into local values
  Number_Of_Lines = Num_Scans_Per_Segment
  Number_Of_Elements = Num_Pix

  Glint_Mask = Missing_Value_Int1

     line_loop: DO Line_Idx = 1, Number_Of_Lines

     element_loop: DO Elem_Idx = 1, Number_Of_Elements

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
             cycle
     endif

     !--- initialize valid pixel to no
     Glint_Mask(Elem_Idx,Line_Idx) = sym%NO

     !--- skip land pixels
     if ((Land_Mask(Elem_Idx,Line_Idx) == sym%NO) .and. &
          Snow(Elem_Idx,Line_Idx) == sym%NO_SNOW) then

       !--- turn on in geometric glint cone and sufficient Ref_Ch1
       if ((Glintzen(Elem_Idx,Line_Idx) < Glint_Zen_Thresh)) then

          !--- assume to be glint if in geometric zone
          Glint_Mask(Elem_Idx,Line_Idx) = sym%YES

          if (Chan_On_Flag_Default(31) == sym%YES) then

            !--- exclude pixels colder than the freezing temperature
            if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < 273.15) then
              Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
              cycle
            endif

            !--- exclude pixels colder than the surface
            if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < ch(31)%Bt_Toa_Clear(Elem_Idx,Line_Idx) - 5.0) then
              Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
              cycle
            endif

          endif

          !-turn off if non-uniform
          if (Chan_On_Flag_Default(31) == sym%YES) then
            if (Bt_Ch31_Std_3x3(Elem_Idx,Line_Idx) > 1.0) then
             Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
             cycle
            endif
          endif

          !-turn off if dark
          if (Chan_On_Flag_Default(1) == sym%YES) then
            if (ch(1)%Ref_Toa(Elem_Idx,Line_Idx) < 5.0) then
             Glint_Mask(Elem_Idx,Line_Idx) = sym%NO
             cycle
            endif
          endif

       endif  !Glintzen check

     endif    !land check

     enddo element_loop
   enddo line_loop


 end subroutine COMPUTE_GLINT
!----------------------------------------------------------------------
!--- Compute a mask identifying presence of oceanic glint in lunar ref
!--- 
!--- input and output passed through global arrays
!----------------------------------------------------------------------
subroutine COMPUTE_GLINT_LUNAR()

  !--- define local variables
  integer:: Number_Of_Lines
  integer:: Number_Of_Elements
  integer:: Elem_Idx
  integer:: Line_Idx


  !--- alias some global sizes into local values
  Number_Of_Lines = Num_Scans_Per_Segment
  Number_Of_Elements = Num_Pix

  Glint_Mask_Lunar = Missing_Value_Int1

     line_loop: DO Line_Idx = 1, Number_Of_Lines

     element_loop: DO Elem_Idx = 1, Number_Of_Elements

     !--- skip bad pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) then
             cycle
     endif

     !--- initialize valid pixel to no
     Glint_Mask_Lunar(Elem_Idx,Line_Idx) = sym%NO


     !--- skip land pixels
     if ((Land_Mask(Elem_Idx,Line_Idx) == sym%NO) .and. &
          Snow(Elem_Idx,Line_Idx) == sym%NO_SNOW) then

       !--- turn on in geometric glint cone and sufficient Ref_Ch1
       if ((Glintzen_Lunar(Elem_Idx,Line_Idx) < Glint_Zen_Thresh)) then

          !--- assume to be glint if in geometric zone
          Glint_Mask_Lunar(Elem_Idx,Line_Idx) = sym%YES

          if (Chan_On_Flag_Default(31) == sym%YES) then

            !--- exclude pixels colder than the freezing temperature
            if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < 273.15) then
              Glint_Mask_Lunar(Elem_Idx,Line_Idx) = sym%NO
              cycle
            endif

            !--- exclude pixels colder than the surface
            if (ch(31)%Bt_Toa(Elem_Idx,Line_Idx) < ch(31)%Bt_Toa_Clear(Elem_Idx,Line_Idx) - 5.0) then
              Glint_Mask_Lunar(Elem_Idx,Line_Idx) = sym%NO
              cycle
            endif

          endif

          !-turn off if non-uniform
          if (Chan_On_Flag_Default(31) == sym%YES) then
            if (Bt_Ch31_Std_3x3(Elem_Idx,Line_Idx) > 1.0) then
             Glint_Mask_Lunar(Elem_Idx,Line_Idx) = sym%NO
             cycle
            endif
          endif

          !-turn off if dark
          if (Chan_On_Flag_Default(42) == sym%YES) then
            if (ch(42)%Ref_Lunar_Toa(Elem_Idx,Line_Idx) < 5.0) then
             Glint_Mask_Lunar(Elem_Idx,Line_Idx) = sym%NO
             cycle
            endif
          endif

       endif  !Glintzen check

     endif    !land check

     enddo element_loop
   enddo line_loop


 end subroutine COMPUTE_GLINT_LUNAR

!----------------------------------------------------------------------
! Read MODIS white sky albedoes
!----------------------------------------------------------------------
subroutine READ_MODIS_WHITE_SKY_ALBEDO(modis_alb_id,modis_alb_str,Ref_Sfc_White_Sky)

    integer(kind=4), intent(in):: modis_alb_id
    TYPE(Land_grid_description), intent(in) :: modis_alb_str
    real(kind=real4), dimension(:,:), intent(out):: Ref_Sfc_White_Sky

    CALL READ_LAND_SFC_HDF(modis_alb_id, modis_alb_str, lat, &
                          lon, Space_Mask, Two_Byte_Temp)
    Ref_Sfc_White_Sky = 0.1*Two_Byte_Temp

    Ref_Sfc_White_Sky = 1.10*Ref_Sfc_White_Sky   !EMPIRICAL ADJUSTMENT

    where(Two_Byte_Temp == 32767)
             Ref_Sfc_White_Sky = Missing_Value_Real4
    endwhere

    !--- modify for water
    where(Land_Mask == sym%NO)
            Ref_Sfc_White_Sky = Ref_Sfc_White_Sky_Water
    endwhere

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


!===========================================================================
!-- Compute Cloud and Land Masked SST
!===========================================================================
subroutine COMPUTE_MASKED_SST(jmin,jmax)

  integer, intent(in):: jmin
  integer, intent(in):: jmax

  integer:: Elem_Idx
  integer:: Line_Idx

  Sst_Masked = Missing_Value_Real4

  line_loop: DO Line_Idx=jmin, jmax - jmin + 1
    element_loop: DO Elem_Idx= 1, Num_Pix

     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
 
     if (Land(Elem_Idx,Line_Idx) /= sym%LAND .and. Land(Elem_Idx,Line_Idx) /= sym%COASTLINE) then
        if (Cld_Mask(Elem_Idx,Line_Idx) == sym%CLEAR .or. Cld_Mask(Elem_Idx,Line_Idx) == sym%PROB_CLEAR) then
          Sst_Masked(Elem_Idx,Line_Idx) = Sst_Unmasked(Elem_Idx,Line_Idx)
        endif
     endif

    end do element_loop
  end do line_loop

 end subroutine COMPUTE_MASKED_SST

!==============================================================================
!
!==============================================================================
 subroutine DETERMINE_LEVEL1B_COMPRESSION(File_1b_Original,L1b_Gzip,L1b_bzip2,File_1b_Full)
   character(len=*), intent(in):: File_1b_Original
   integer(kind=int4), intent(out):: L1b_Gzip
   integer(kind=int4), intent(out):: L1b_bzip2
   character(len=*), intent(out):: File_1b_Full
   character(len=200):: System_String

   character(len=7):: L1b_ext


  !--- determine if the goes data is compressed
  L1b_ext = File_1b_Original(len_trim(File_1b_Original)-2: &
                             len_trim(File_1b_Original))

  !-- determine if gzipped
  if (trim(L1b_ext) == '.gz') then
     L1b_Gzip = sym%YES
  else
     L1b_Gzip = sym%NO
  endif

  !--- check if bzipped
  if (trim(L1b_ext) == 'bz2') then
     L1b_bzip2 = sym%YES
  else
     L1b_bzip2 = sym%NO
  endif

  !--- uncompress
  if (L1b_Gzip == sym%YES) then
     File_1b = File_1b_Original(1:len(trim(File_1b_Original))-3)
     System_String = "gunzip -c "//trim(dir_1b)//trim(File_1b_Original)// &
        " > "//trim(Temporary_Data_Dir)//trim(File_1b)
     CALL SYSTEM(System_String)

     Number_of_Temporary_Files = Number_of_Temporary_Files + 1
     Temporary_File_Name(Number_of_Temporary_Files) = trim(File_1b)

  elseif (L1b_bzip2 == sym%YES) then
     File_1b = File_1b_Original(1:len(trim(File_1b_Original))-4)
     System_String = "bunzip2 -c "//trim(dir_1b)//trim(File_1b_Original)// &
        " > "//trim(Temporary_Data_Dir)//trim(File_1b)
     CALL SYSTEM(System_String)

     Number_of_Temporary_Files = Number_of_Temporary_Files + 1
     Temporary_File_Name(Number_of_Temporary_Files) = trim(File_1b)

  else
     File_1b = trim(File_1b_Original)
  endif

   !--- make a full file name
   if (L1b_Gzip == sym%YES .or. L1b_bzip2 == sym%YES) then
     File_1b_Full = trim(Temporary_Data_Dir)//trim(File_1b)
   else
     File_1b_Full = trim(dir_1b)//trim(File_1b)
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
  Num_Elements = Num_Pix
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  line_loop: do Line_Idx = Line_Idx_Min, Line_Idx_Max
    element_loop: do Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

    if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle
    if (Chan_On_Flag_Default(1) == sym%NO) cycle
    if (Chan_On_Flag_Default(2) == sym%NO) cycle
    if (Solzen(Elem_Idx,Line_Idx) > Solzen_Threshold) cycle
    if (Modis_Flag == sym%YES) cycle                         !modis ch2 saturates, need to modify for MODIS

    Ndvi_Temp = (ch(2)%Ref_Toa(Elem_Idx,Line_Idx) - ch(1)%Ref_Toa(Elem_idx,Line_Idx)) / &
                (ch(2)%Ref_Toa(Elem_Idx,Line_Idx) + ch(1)%Ref_Toa(Elem_idx,Line_Idx)) 

    if (Ndvi_Temp > Ndvi_Land_Threshold) then
      Land(Elem_Idx,Line_Idx) = sym%LAND
    endif

    if (Ndvi_Temp < Ndvi_Water_Threshold .and. &
        Land(Elem_Idx,Line_Idx) == sym%LAND) then
        Land(Elem_Idx,Line_Idx) = sym%SHALLOW_INLAND_WATER
    endif

    enddo element_loop
  enddo line_loop

end subroutine MODIFY_LAND_CLASS_WITH_NDVI
!====================================================================
! Function Name: TERM_REFL_NORM
!
! Function:
!    Renormalize reflectances to improve performance near the terminator 
! using the parameteization given by Li et. al. 2006
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
! Reference: Li et. al. 2006
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
  Num_Elements = Num_Pix
  Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
  Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  line_loop: DO Line_Idx = Line_Idx_Min, Line_Idx_Max
    element_loop: DO Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

      !--- if no, geolocation, set to missing and go to next pixel
      if (Space_Mask(Elem_Idx,Line_Idx) == sym%YES) then
          Zsfc(Elem_Idx,Line_Idx) = Missing_Value_Real4
          cycle
      endif
 
      !--- if hires value available use it, or try NWP
      if (Zsfc_Hires(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then 

           Zsfc(Elem_Idx,Line_Idx) = Zsfc_Hires(Elem_Idx,Line_Idx)

      !--- if this is ocean pixel, assume zero
      elseif (Land(Elem_Idx,Line_Idx) == sym%SHALLOW_OCEAN .or. &
              Land(Elem_Idx,Line_Idx) == sym%MODERATE_OCEAN .or. &
              Land(Elem_Idx,Line_Idx) == sym%DEEP_OCEAN) then

              Zsfc(Elem_Idx,Line_Idx) = 0.0     

      !--- try NWP
      else

         Lon_Nwp_Idx = I_Nwp(Elem_Idx,Line_Idx)
         Lat_Nwp_Idx = J_Nwp(Elem_Idx,Line_Idx)

         !--- if nwp not available, assume zero
         if (Lon_Nwp_Idx > 0 .and. Lat_Nwp_Idx > 0) then
           Zsfc(Elem_Idx,Line_Idx) = Zsfc_Nwp(Lon_Nwp_Idx, Lat_Nwp_Idx)
         else
           Zsfc(Elem_Idx,Line_Idx) = 0.0     
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

  Number_of_Elements = Num_Pix

  Adj_Pix_Cld_Mask = Missing_Value_Int1

  line_loop: DO Line_Idx = Line_Start, Number_of_Lines + Line_Start - 1

    j1 = max(1,Line_Idx - 1)
    j2 = min(Number_of_Lines,Line_Idx + 1)

    element_loop: DO Elem_Idx = 1, Number_of_Elements

      i1 = max(1,Elem_Idx - 1)
      i2 = min(Number_of_Elements,Elem_Idx + 1)

      Adj_Pix_Cld_Mask(Elem_Idx,Line_Idx) = maxval(Cld_Mask(i1:i2,j1:j2))

    enddo element_loop
  enddo line_loop

end subroutine ADJACENT_PIXEL_CLOUD_MASK
!---------------------------------------------------------------------------------------
! COMPUTE_VIIRS_SST is used to compute Sst_Masked for VIIRS
! It purpouse is to replicate VIIRS SST team algorithm results
! called from process_clavrx.f90
!---------------------------------------------------------------------------------------
subroutine COMPUTE_VIIRS_SST(jmin,jmax)

  integer, intent(in):: jmin
  integer, intent(in):: jmax

  integer:: Elem_Idx
  integer:: Line_Idx
  real:: a0_day = 2.74438
  real:: a1_day = 0.996818
  real:: a2_day = 0.0720411
  real:: a3_day = 1.35385
  real:: a0_night = -5.31302
  real:: a1_night = 0.788114
  real:: a2_night = 1.04555
  real:: a3_night = -0.810677
  real:: a4_night = 0.145501
  real:: a5_night = 1.11931

  Sst_Unmasked = Missing_Value_Real4

  if (Chan_On_Flag_Default(20) == sym%NO) return
  if (Chan_On_Flag_Default(31) == sym%NO) return
  if (Chan_On_Flag_Default(32) == sym%NO) return

  line_loop: DO Line_Idx=jmin, jmax + jmin - 1
    element_loop: DO Elem_Idx= 1, Num_Pix

  ! Day time
     if (Land(Elem_Idx,Line_Idx) /= sym%LAND .and. Land(Elem_Idx,Line_Idx) /= sym%COASTLINE) then
        if (Solzen(Elem_Idx,Line_Idx) < 90.0 .and. sst_anal(Elem_Idx,Line_Idx) > 150.0) then
           if (Chan_On_Flag_Default(31) == sym%YES .and. Chan_On_Flag_Default(32) == sym%YES) then
              Sst_Unmasked(Elem_Idx,Line_Idx) = a0_day + a1_day * ch(31)%Bt_Toa(Elem_Idx,Line_Idx) + &
                   a2_day * Btd_Ch31_Ch32(Elem_Idx,Line_Idx) * (sst_anal(Elem_Idx,Line_Idx) - 273.15) + &
                   a3_day * Btd_Ch31_Ch32(Elem_Idx,Line_Idx) * (Seczen(Elem_Idx,Line_Idx) - 1)
           endif
        endif

  ! Night time
        if (Solzen(Elem_Idx,Line_Idx) .ge. 90.0) then
           if (Chan_On_Flag_Default(20) == sym%YES .and. Chan_On_Flag_Default(31) == sym%YES .and. &
                 Chan_On_Flag_Default(32) == sym%YES) then

                   Sst_Unmasked(Elem_Idx,Line_Idx) = a0_night + a1_night * ch(31)%Bt_Toa(Elem_Idx,Line_Idx) + &
                   a2_night * ch(20)%Bt_Toa(Elem_Idx,Line_Idx) + a3_night * ch(32)%Bt_Toa(Elem_Idx,Line_Idx) + &
                   a4_night * (ch(20)%Bt_Toa(Elem_Idx,Line_Idx) - ch(32)%Bt_Toa(Elem_Idx,Line_Idx)) * &
                   (Seczen(Elem_Idx,Line_Idx)-1) + a5_night * (Seczen(Elem_Idx,Line_Idx) - 1)
           endif 
        endif
     endif

    end do element_loop
  end do line_loop

end subroutine COMPUTE_VIIRS_SST

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
subroutine COMPUTE_ACHA_PERFORMANCE_METRICS(Acha_Processed_Count,Acha_Valid_Count)

  real(kind=real4), intent(inout):: ACHA_Processed_Count
  real(kind=real4), intent(inout):: ACHA_Valid_Count
  integer, parameter:: Count_Min = 10
  real:: Processed_Count_Segment
  real:: Valid_Count_Segment

  Processed_Count_Segment = count(btest(ACHA_Packed_Quality_Flags,0))

  Valid_Count_Segment = count((btest(int(ACHA_Packed_Quality_Flags),1)) .and.  &
                              (btest(int(ACHA_Packed_Quality_Flags),2)) .and. &
                              (btest(int(ACHA_Packed_Quality_Flags),3)))
  
  ACHA_Processed_Count = ACHA_Processed_Count + Processed_Count_Segment
  ACHA_Valid_Count = ACHA_Valid_Count + Valid_Count_Segment

  if (ACHA_Processed_Count > Count_Min) then
    ACHA_Success_Fraction = ACHA_Valid_Count / ACHA_Processed_Count
  else
    ACHA_Success_Fraction = Missing_Value_Real4 
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
  integer(kind=int1), dimension(:,:), allocatable:: Mask
  integer(kind=int1), dimension(:,:), allocatable:: Nonconfident_Mask
  integer, parameter:: Count_Min = 10
  real:: Count_Segment
  real:: Nonconfident_Count_Segment
  

  Num_Elements = Num_Pix  !make local copy of a global variable
  Num_Lines = Num_Scans_Per_Segment  !make local copy of a global variable

  allocate(Mask(Num_Elements,Num_Lines))
  allocate(Nonconfident_Mask(Num_Elements,Num_Lines))

  Mask = 0
  Nonconfident_Mask = 0

  where(Cld_Mask == sym%CLEAR .or. Cld_Mask == sym%PROB_CLEAR .or. Cld_Mask == sym%PROB_CLOUDY .or. Cld_Mask == sym%Cloudy)
     Mask = 1
  endwhere

  where(Cld_Mask == sym%PROB_CLEAR .or. Cld_Mask == sym%PROB_CLOUDY)
     Nonconfident_Mask = 1
  endwhere

  Count_Segment = sum(real(Mask))
  if (Count_Segment < Count_Min) then
      return
  endif

  Nonconfident_Count_Segment = sum(real(Nonconfident_Mask))

  Cloud_Mask_Count = Cloud_Mask_Count + Count_Segment
  Nonconfident_Cloud_Mask_Count = Nonconfident_Cloud_Mask_Count + Nonconfident_Count_Segment

  if (Cloud_Mask_Count > Count_Min) then
    Nonconfident_Cloud_Mask_Fraction = Nonconfident_Cloud_Mask_Count / Cloud_Mask_Count
  else
    Nonconfident_Cloud_Mask_Fraction = Missing_Value_Real4 
  endif
  

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

  if (Chan_On_Flag_Default(31) == sym%YES) then
     call COMPUTE_TSFC(Line_Idx_Min,Line_Idx_Max)
  endif

  if (Chan_On_Flag_Default(1) == sym%YES .and. &
      Chan_On_Flag_Default(6) == sym%YES) then

       Ndsi_Toa = NORMALIZED_DIFFERENCE_SNOW_INDEX(  &
                             ch(1)%Ref_Toa, &
                             ch(6)%Ref_Toa, &
                             Solzen)

       Ndsi_Sfc = NORMALIZED_DIFFERENCE_SNOW_INDEX(  &
                             ch(1)%Ref_Sfc, &
                             ch(6)%Ref_Sfc, &
                             Solzen)
  endif

  if (Chan_On_Flag_Default(1) == sym%YES .and. &
      Chan_On_Flag_Default(2) == sym%YES) then

       Ndvi_Toa = NORMALIZED_DIFFERENCE_VEGETATION_INDEX(  &
                             ch(1)%Ref_Toa, &
                             ch(2)%Ref_Toa, &
                             Solzen)

       Ndvi_Sfc = NORMALIZED_DIFFERENCE_VEGETATION_INDEX(  &
                             ch(1)%Ref_Sfc, &
                             ch(2)%Ref_Sfc, &
                             Solzen)

       Ndvi_Sfc_White_Sky = NORMALIZED_DIFFERENCE_VEGETATION_INDEX(  &
                             ch(1)%Sfc_Ref_White_Sky, &
                             ch(2)%Sfc_Ref_White_Sky, &
                             Solzen)

       Rsr = REMOTE_SENSING_REFLECTANCE( ch(1)%Ref_Sfc, &
                                         ch(2)%Ref_Sfc, &
                                         Airmass,       &
                                         Solzen)
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
   if (allocated(ch(42)%Rad_Toa)) then
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
! end of MODULE
!-----------------------------------------------------------
end MODULE PIXEL_ROUTINES
