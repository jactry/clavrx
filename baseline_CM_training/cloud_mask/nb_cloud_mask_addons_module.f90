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
! AUTHORS:  The PATMOS-x Team
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
module CLOUD_MASK_ADDONS

 use NB_CLOUD_MASK_SERVICES

 implicit none

 public:: NB_CLOUD_MASK_ADDONS_ALGORITHM

 private:: EUMETSAT_FIRE_TEST, &
           CLAVRX_SMOKE_OVER_WATER_TEST, &
           CLAVRX_SMOKE_OVER_LAND_TEST, &
           CLAVRX_DUST_TEST, &
           CLAVRX_THIN_CIRRUS_TEST

 !--- define structures
 include 'nb_cloud_mask.inc'

 contains
 !---------------------------------------------------------------------
 ! ubroutine that handle non-cloud tests in the enterprise mask
 ! These include smoke, dust and fire
 !---------------------------------------------------------------------
 subroutine NB_CLOUD_MASK_ADDONS_ALGORITHM( &
            Symbol, &                       !local copy of sym structure
            Input,  &
            Output,  &
            Diag)

  type(symbol_naive_bayesian), intent(in) :: Symbol
  type(mask_input), intent(in) :: Input
  type(mask_output), intent(out) :: Output
  type(diag_output), intent(out), Optional :: Diag

  real :: Min_Thresh
  real :: Max_Thresh

  !------------------------------------------------------------------------------------------
  !---  begin executable code
  !------------------------------------------------------------------------------------------

  !--- initialize to missing
  Output%Dust_Mask = MISSING_VALUE_INT1
  Output%Smoke_Mask = MISSING_VALUE_INT1
  Output%Fire_Mask = MISSING_VALUE_INT1
  Output%Thin_Cirr_Mask = MISSING_VALUE_INT1

   !--- initialize diagnostic output
   if (present(Diag)) Diag%Array_1 = Missing_Value_Real4
   if (present(Diag)) Diag%Array_2 = Missing_Value_Real4
   if (present(Diag)) Diag%Array_3 = Missing_Value_Real4

  !------------------------------------------------------------------------------------------
  !  MAIN LOOP
  !------------------------------------------------------------------------------------------

   !--- check for valid data
   if (Input%Invalid_Data_Mask == Symbol%NO) then

        !----------------------------------------------------------------------------
        ! CLAVR-x SMOKE TEST
        !
        ! Day and ice-free ocean only
        !
        !-------------------------------------------------------------------------
        if (Input%SNOW_Class == symbol%NO_SNOW .and. &
            (Input%Land_Class == symbol%DEEP_OCEAN .or. Input%Land_Class == symbol%MODERATE_OCEAN)) then

           Output%Smoke_Mask = 0 

           Output%Smoke_Mask = CLAVRX_SMOKE_OVER_WATER_TEST(Input%Ref_063um,Input%Ref_063um_Clear, &
                                                 Input%Ref_138um, &
                                                 Input%Ref_160um,Input%Ref_160um_Clear, &
                                                 Input%Ref_375um,Input%Ref_375um_Clear, &
                                                 Input%Bt_375um, &
                                                 Input%Bt_11um, Input%Bt_11um_Clear, &
                                                 Input%Bt_12um, Input%Bt_12um_Clear,&
                                                 Input%Chan_On_063um,Input%Chan_On_138um, &
                                                 Input%Chan_On_160um,Input%Chan_On_375um,Input%Chan_On_11um, &
                                                 Input%Chan_On_12um,Input%Emiss_11um_Tropo, &
                                                 Input%Ref_063um_Std, Input%Bt_11um_Std, &
                                                 Input%Solzen)

        endif

        if (Input%SNOW_Class == symbol%NO_SNOW .and. &
            (Input%Land_Class == symbol%LAND)) then

           Output%Smoke_Mask = 0 

           Output%Smoke_Mask = CLAVRX_SMOKE_OVER_LAND_TEST(Input%Ref_063um,Input%Ref_063um_Clear, &
                                                 Input%Ref_086um, &
                                                 Input%Ref_138um, &
                                                 Input%Ref_160um,Input%Ref_160um_Clear, &
                                                 Input%Ref_375um,Input%Ref_375um_Clear, &
                                                 Input%Bt_375um, &
                                                 Input%Bt_11um, Input%Bt_11um_Clear, &
                                                 Input%Bt_12um, Input%Bt_12um_Clear,&
                                                 Input%Chan_On_063um,Input%Chan_On_138um, &
                                                 Input%Chan_On_160um,Input%Chan_On_375um,Input%Chan_On_11um, &
                                                 Input%Chan_On_12um,Input%Emiss_11um_Tropo, &
                                                 Input%Ref_063um_Std, Input%Bt_11um_Std, &
                                                 Input%Solzen)

        endif

        !-------------------------------------------------------------------------
        !-- IR Dust algorithm
        !-------------------------------------------------------------------------
        if (Input%SNOW_Class == symbol%NO_SNOW .and. Input%Chan_On_12um == symbol%YES .and. &
            (Input%Land_Class == symbol%DEEP_OCEAN .or. Input%Land_Class == symbol%MODERATE_OCEAN)) then

              Output%Dust_Mask = 0

              Output%Dust_Mask = CLAVRX_DUST_TEST ( &
                    Input%Chan_On_85um,         &
                    Input%Chan_On_11um,         &
                    Input%Chan_On_12um,         &
                    Input%Bt_85um,            &
                    Input%Bt_11um,            &
                    Input%Bt_12um,            &
                    Input%Bt_11um_Clear,      &
                    Input%Bt_12um_Clear,      &
                    Input%Bt_11um_Std,        &
                    Input%Emiss_11um_Tropo)
        endif
 
      !-----------------------------------------------------------------------------
      ! EUMETSAT Fire Algorithm
      !-----------------------------------------------------------------------------
      if (Input%Chan_On_11um == symbol%YES .and.   &
          Input%Chan_On_375um == symbol%YES .and.  &
          Input%Land_Class == Symbol%Land .and.    &
          Input%Sfc_Type /= Symbol%BARE_SFC .and. &
          Input%Sfc_Type /= Symbol%OPEN_SHRUBS_SFC) then

          Output%Fire_Mask = 0

          Output%Fire_Mask = EUMETSAT_FIRE_TEST ( &
                   Input%Bt_11um &
                 , Input%Bt_375um &
                 , Input%Bt_11um_Std &
                 , Input%Bt_375um_Std &
                 , Input%Solzen )

      endif

      !-----------------------------------------------------------------------------
      ! VCM Thin Cirrus Algorithm (only for daytime)
      !-----------------------------------------------------------------------------
      if (Input%Chan_On_138um == symbol%YES .and.  &
          Input%Solzen <= Reflectance_Spatial_Solzen_Thresh) then
           
          Output%Thin_Cirr_Mask = 0

          ! --- get thresholds depending on the surface conditions
          Min_Thresh = Thin_Cirr_Min_Thresh
          Max_Thresh = Thin_Cirr_Max_Thresh

          ! desert
          if (Input%Sfc_Type == Symbol%BARE_SFC) then
            Min_Thresh = Desert_Thin_Cirr_Min_Thresh
            Max_Thresh = Desert_Thin_Cirr_Max_Thresh
          endif

          ! snow/ice
          if (Input%Snow_Class == Symbol%SNOW .or. &
              Input%Snow_Class == Symbol%SEA_ICE ) then
            Min_Thresh = Snow_Thin_Cirr_Min_Thresh
            Max_Thresh = Snow_Thin_Cirr_Max_Thresh
          endif

          Output%Thin_Cirr_Mask = CLAVRX_THIN_CIRRUS_TEST ( &
                   Input%Ref_138um,   &
                   Min_Thresh,        &
                   Max_Thresh )

      endif

    endif ! if Invalid_Data_Mask = NO

 end subroutine NB_CLOUD_MASK_ADDONS_ALGORITHM

!----------------------------------------------------------------------------
! CLAVR-x SMOKE TEST
!
! Daytime and Ice-Free Ocean Only.
!
! Coded: Andy Heidinger
!----------------------------------------------------------------------------
  integer elemental function CLAVRX_SMOKE_OVER_WATER_TEST(Refl_065,Refl_065_Clear,Refl_138, &
                                               Refl_160, Refl_160_Clear, &
                                               Refl_375, Refl_375_Clear, &
                                               Bt_375,&
                                               Bt_11, Bt_11_Clear, &
                                               Bt_12, Bt_12_Clear,&
                                               Chan_On_065,Chan_On_138,Chan_On_160,Chan_On_375,Chan_On_11, &
                                               Chan_On_12,Emiss_11_Tropo,Refl_065_Std,T11_Std,Solzen)
     real, intent(in):: Refl_065
     real, intent(in):: Refl_065_Clear
     real, intent(in):: Refl_138
     real, intent(in):: Refl_160
     real, intent(in):: Refl_160_Clear
     real, intent(in):: Refl_375
     real, intent(in):: Refl_375_Clear
     real, intent(in):: Bt_375
     real, intent(in):: Bt_11
     real, intent(in):: Bt_11_Clear
     real, intent(in):: Bt_12
     real, intent(in):: Bt_12_Clear
     real, intent(in):: Emiss_11_Tropo
     real, intent(in):: Refl_065_Std
     real, intent(in):: T11_Std
     integer, intent(in):: Chan_On_065
     integer, intent(in):: Chan_On_138
     integer, intent(in):: Chan_On_160
     integer, intent(in):: Chan_On_375
     integer, intent(in):: Chan_On_11
     integer, intent(in):: Chan_On_12
     real, intent(in):: Solzen

     integer:: IR_Flag
     integer:: VIS_Flag
     integer:: NIR_Flag
     integer:: NIR_IR_Flag
     integer:: SPLIT_WIN_Flag
     
     Clavrx_Smoke_Over_Water_Test = MISSING_VALUE_INT1

     if (Solzen < Solzen_Max_Smoke_Water_Thresh) then

           IR_Flag = 1
           VIS_Flag = 1
           NIR_Flag = 1
           NIR_IR_Flag = 1
           SPLIT_WIN_Flag = 1

           !--- IR test - smoke should be nearly invisible
           if (Chan_On_11 == 1) then
               if (Emiss_11_Tropo > Emiss_11_Tropo_Max_Smoke_Water_Thresh) IR_Flag = 0
               if (T11_Std > T11_Std_Max_Smoke_Water_Thresh) IR_Flag = 0
           endif

           !--- VIS test - smoke should be nearly invisible
           if (Chan_On_065 == 1) then
               if (Refl_065 - Refl_065_Clear > Refl_065_Max_Smoke_Water_Thresh .or.   &
                   Refl_065 - Refl_065_Clear < Refl_065_Min_Smoke_Water_Thresh) VIS_Flag = 0
               if (Refl_065_Std  > Refl_065_Std_Max_Smoke_Water_Thresh) VIS_Flag = 0
           endif

           !--- NIR Tests
           if (Chan_On_375 == 1) then
               if (Refl_375 - Refl_375_Clear > Refl_375_Max_Smoke_Water_Thresh) NIR_Flag = 0
               if (Refl_375 > Refl_375_Max_Smoke_Water_Thresh) NIR_Flag = 0
           endif

           if (Chan_On_160 == 1 .and. Chan_On_375 == 0) then
               if (Refl_160 > Refl_160_Max_Smoke_Water_Thresh) NIR_Flag = 0
           endif

           if (Chan_On_138 == 1) then
               if (Refl_138 > Refl_138_Max_Smoke_Water_Thresh) NIR_Flag = 0
           endif

           !--- NIR IR_Tests
           if (Chan_On_375 == 1 .and. Chan_On_11 == 1) then
               if ((Bt_375 - Bt_11) > Btd_4_11_Max_Smoke_Water_Thresh) NIR_IR_Flag = 0
           endif

           !--- SPLIT_WIN_Tests
           if (Chan_On_11 == 1 .and. Chan_On_12 == 1) then
             if (abs(split_window_test(Bt_11_Clear, Bt_12_Clear,Bt_11, Bt_12)) > 1) SPLIT_WIN_Flag = 0
           endif

           !--- combine into final answer
           Clavrx_Smoke_Over_Water_Test = IR_Flag * VIS_Flag * NIR_Flag * NIR_IR_Flag * SPLIT_WIN_Flag
   
    endif

  end function CLAVRX_SMOKE_OVER_WATER_TEST

!----------------------------------------------------------------------------
! CLAVR-x SMOKE TEST OVER LAND
!
! Daytime and Snow-Free Land Only.
!
! Coded: Andy Heidinger
!----------------------------------------------------------------------------
  integer elemental function CLAVRX_SMOKE_OVER_LAND_TEST(Refl_065,Refl_065_Clear, &
                                               Refl_086, Refl_138, &
                                               Refl_160, Refl_160_Clear, &
                                               Refl_375, Refl_375_Clear, &
                                               Bt_375,&
                                               Bt_11, Bt_11_Clear, &
                                               Bt_12, Bt_12_Clear,&
                                               Chan_On_065,Chan_On_138,Chan_On_160,Chan_On_375,Chan_On_11, &
                                               Chan_On_12,Emiss_11_Tropo,Refl_065_Std,T11_Std,Solzen)
     real, intent(in):: Refl_065
     real, intent(in):: Refl_065_Clear
     real, intent(in):: Refl_086
     real, intent(in):: Refl_138
     real, intent(in):: Refl_160
     real, intent(in):: Refl_160_Clear
     real, intent(in):: Refl_375
     real, intent(in):: Refl_375_Clear
     real, intent(in):: Bt_375
     real, intent(in):: Bt_11
     real, intent(in):: Bt_11_Clear
     real, intent(in):: Bt_12
     real, intent(in):: Bt_12_Clear
     real, intent(in):: Emiss_11_Tropo
     real, intent(in):: Refl_065_Std
     real, intent(in):: T11_Std
     integer, intent(in):: Chan_On_065
     integer, intent(in):: Chan_On_138
     integer, intent(in):: Chan_On_160
     integer, intent(in):: Chan_On_375
     integer, intent(in):: Chan_On_11
     integer, intent(in):: Chan_On_12
     real, intent(in):: Solzen
  

     real:: NIR_Smoke_Ratio
     integer:: IR_Flag
     integer:: VIS_Flag
     integer:: NIR_Flag
     integer:: NIR_IR_Flag
     integer:: SPLIT_WIN_Flag
     
     Clavrx_Smoke_Over_Land_Test = MISSING_VALUE_INT1

     if (Solzen < Solzen_Max_Smoke_Land_Thresh) then

           IR_Flag = 1
           VIS_Flag = 1
           NIR_Flag = 1
           NIR_IR_Flag = 1
           SPLIT_WIN_Flag = 1

           !--- IR test - smoke should be nearly invisible
           if (Chan_On_11 == 1) then
               if (Emiss_11_Tropo > Emiss_11_Tropo_Max_Smoke_Land_Thresh) IR_Flag = 0
               if (T11_Std > T11_Std_Max_Smoke_Land_Thresh) IR_Flag = 0
           endif

           !--- VIS test - smoke should be nearly invisible
           if (Chan_On_065 == 1) then
               if (Refl_065 - Refl_065_Clear > Refl_065_Max_Smoke_Land_Thresh .or.   &
                   Refl_065 - Refl_065_Clear < Refl_065_Min_Smoke_Land_Thresh) VIS_Flag = 0
               if (Refl_065_Std  > Refl_065_Std_Max_Smoke_Land_Thresh) VIS_Flag = 0
           endif

           !--- NIR Tests
           if (Chan_On_375 == 1) then
               if (Refl_375 - Refl_375_Clear > Refl_375_Max_Smoke_Land_Thresh) NIR_Flag = 0
               if (Refl_375 > Refl_375_Max_Smoke_Land_Thresh) NIR_Flag = 0
           endif

           !if (Chan_On_160 == 1 .and. Chan_On_375 == 0) then
           if (Chan_On_160 == 1) then
               NIR_Smoke_Ratio = (Refl_160 - Refl_065)/(Refl_086 - Refl_065)
               if (NIR_Smoke_Ratio > NIR_Smoke_Ratio_Max_Land_Thresh) NIR_Flag = 0
           endif

           if (Chan_On_138 == 1) then
               if (Refl_138 > Refl_138_Max_Smoke_Land_Thresh) NIR_Flag = 0
           endif

           !--- NIR IR_Tests
           if (Chan_On_375 == 1 .and. Chan_On_11 == 1) then
               if ((Bt_375 - Bt_11) > Btd_4_11_Max_Smoke_Land_Thresh) NIR_IR_Flag = 0
           endif

           !--- SPLIT_WIN_Tests
           if (Chan_On_11 == 1 .and. Chan_On_12 == 1) then
             if (abs(split_window_test(Bt_11_Clear, Bt_12_Clear,Bt_11, Bt_12)) > 1) SPLIT_WIN_Flag = 0
           endif

           !--- combine into final answer
           Clavrx_Smoke_Over_Land_Test = IR_Flag * VIS_Flag * NIR_Flag * NIR_IR_Flag * SPLIT_WIN_Flag
   
    endif

  end function CLAVRX_SMOKE_OVER_LAND_TEST

!----------------------------------------------------------------------------
! EUMETCAST Fire detection algorithm
!
! Reference: This implements the "Current Operational Algorithm" described in:
! TOWARDS AN IMPROVED ACTIVE FIRE MONITORING PRODUCT FOR MSG SATELLITES
! Sauli Joro, Olivier Samain, Ahmet Yildirim, Leo van de Berg, Hans Joachim Lutz
! EUMETSAT, Am Kavalleriesand 31, Darmstadt, Germany
!
! Coded by William Straka III
!
!-----------------------------------------------------------------------------
  integer elemental function EUMETSAT_FIRE_TEST(T11,T375,T11_Std,T375_Std,Solzen)

     real, intent(in):: T11
     real, intent(in):: T375
     real, intent(in):: T11_Std
     real, intent(in):: T375_Std
     real, intent(in):: Solzen

     real :: Bt_375um_Eumet_Fire_Thresh
     real :: Bt_Diff_Eumet_Fire_Thresh
     real :: Stddev_11um_Eumet_Fire_Thresh
     real :: Stddev_375um_Eumet_Fire_Thresh

     !--- initialize
     Eumetsat_Fire_Test = 0

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
         Eumetsat_Fire_Test = 1
       endif

     endif

  end function EUMETSAT_FIRE_TEST

   !---------------------------------------------------------------------------
   ! CLAVR-x IR DUST Algorithm
   !
   ! Coded: Andy Heidinger
   !---------------------------------------------------------------------------
   integer elemental function CLAVRX_DUST_TEST ( &
                    Chan_On_85,         &
                    Chan_On_11,         &
                    Chan_On_12,         &
                    Bt_85,            &
                    Bt_11,            &
                    Bt_12,            &
                    Bt_11_Clear,      &
                    Bt_12_Clear,      &
                    Bt_11_Std, &
                    Emiss_11_Tropo)

   integer , intent(in) :: Chan_On_85
   integer , intent(in) :: Chan_On_11
   integer , intent(in) :: Chan_On_12
   real , intent(in) :: Bt_85
   real , intent(in) :: Bt_11
   real , intent(in) :: Bt_12
   real , intent(in) :: Bt_11_Clear
   real , intent(in) :: Bt_12_Clear
   real , intent(in) :: Bt_11_Std
   real , intent(in) :: Emiss_11_Tropo
   real :: Btd_11_12
   real :: Btd_11_12_Metric
   integer:: Split_Win_Flag
   integer:: IR_Win_Flag
   integer:: IR_Win_Std_Flag
   integer:: IR_Win_85_Diff_Flag

   Split_Win_Flag = 1
   IR_Win_Flag = 1
   IR_Win_Std_Flag = 1
   IR_Win_85_Diff_Flag = 1

   !--- Split Window
   if (Chan_On_11 == 1 .and. Chan_On_12 == 1) then
      Btd_11_12 = (Bt_11 - Bt_12)
      Btd_11_12_Metric = split_window_test(Bt_11_Clear, Bt_12_Clear,Bt_11, Bt_12)
      if (Btd_11_12 > Btd_11_12_Max_Dust_Thresh) Split_Win_Flag = 0
      if (Btd_11_12_Metric > Btd_11_12_Metric_Max_Dust_Thresh) Split_Win_Flag = 0
      if ((Bt_11 - Bt_12) - (Bt_11_Clear - Bt_12_Clear) > Bt_11_12_Clear_Diff_Max_Dust_Thresh) Split_Win_Flag = 0
   endif

   !--- 8.5-11 should be moderately negative.  ice clouds are positive
   !--- water clouds are very negative
   if (Chan_On_11 == 1 .and. chan_On_12 == 1 .and. Chan_On_85 == 1) then
      if ((Bt_85 - Bt_11) > Btd_85_11_Max_Dust_Thresh) IR_Win_85_Diff_Flag = 0
      if ((Bt_85 - Bt_11) < Btd_85_11_Min_Dust_Thresh) IR_Win_85_Diff_Flag = 0
   endif

   !--- 11um variabilty
   if (Chan_On_11 == 1) then
      if (Bt_11_Std > Bt_11_Std_Max_Dust_Thresh) IR_Win_Std_Flag = 0
   endif

   !--- IR test - dust should low to moderate emissivity
   if (Chan_On_11 == 1) then
       if (Emiss_11_Tropo > Emiss_11_Tropo_Max_Dust_Thresh) IR_Win_Flag = 0
       if (Emiss_11_Tropo < Emiss_11_Tropo_Min_Dust_Thresh) IR_Win_Flag = 0
       if (Bt_11 - Bt_11_Clear < Bt_11_Clear_Diff_Min_Dust_Thresh) IR_Win_Flag = 0
   endif

   CLAVRX_DUST_TEST = Split_Win_Flag * IR_Win_Std_Flag * IR_Win_Flag * IR_Win_85_Diff_Flag

   end function CLAVRX_DUST_TEST


  !-------------------------------------------------------------------------------------
  !  CLAVR-x Split Window Test for Clear value for a given T11
  !-------------------------------------------------------------------------------------
  real elemental function SPLIT_WINDOW_TEST ( t11_clear, t12_clear, t11, t12)

     real, intent(in):: t11_clear
     real, intent(in):: t12_clear
     real, intent(in):: t11
     real, intent(in):: t12

     split_window_test  = (t11_clear - t12_clear) * (t11 - 260.0) / (t11_clear - 260.0) 

     if (t11_clear <=265.0) split_window_test = 0.0

     split_window_test = (t11 - t12) - split_window_test

  end function SPLIT_WINDOW_TEST

  !-------------------------------------------------------------------------------------
  !  CLAVR-x Thin Cirrus Test
  !-------------------------------------------------------------------------------------
  integer elemental function CLAVRX_THIN_CIRRUS_TEST ( r138, min_thresh, max_thresh )

     real, intent(in):: r138
     real, intent(in):: min_thresh
     real, intent(in):: max_thresh

     ! --- calculate test result
     if ( r138 >= min_thresh .and. &
          r138 <= max_thresh ) then
        clavrx_thin_cirrus_test = 1
     else
        clavrx_thin_cirrus_test = 0
     endif

  end function CLAVRX_THIN_CIRRUS_TEST

!-----------------------------------------------------------
! end of module
!-----------------------------------------------------------
end module CLOUD_MASK_ADDONS
