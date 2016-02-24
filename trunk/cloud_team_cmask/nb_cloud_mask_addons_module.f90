! $Id: pixel_routines.f90 538 2014-09-15 20:30:48Z dbotambekov $
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
! Public routines used in this MODULE:
!
!--------------------------------------------------------------------------------------
module CLOUD_MASK_ADDONS

 use NB_CLOUD_MASK_SERVICES

 implicit none

 public:: NB_CLOUD_MASK_ADDONS_ALGORITHM

 private:: EUMETSAT_FIRE_TEST, &
           CLAVRX_SMOKE_TEST, &
           CLAVRX_DUST_TEST, &
           IR_DUST_TEST2

!!--- define structures
 include 'nb_cloud_mask.inc'

 contains
 !---------------------------------------------------------------------
 !
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

  !------------------------------------------------------------------------------------------
  !---  begin executable code
  !------------------------------------------------------------------------------------------

  !--- initialize to missing
  Output%Dust_Mask = MISSING_VALUE_INT1
  Output%Smoke_Mask = MISSING_VALUE_INT1
  Output%Fire_Mask = MISSING_VALUE_INT1

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
        Output%Smoke_Mask = 0 ! MISSING_VALUE_INT1
        if (Input%SNOW_Class == symbol%NO_SNOW .and. &
            (Input%Land_Class == symbol%DEEP_OCEAN .or. Input%Land_Class == symbol%MODERATE_OCEAN)) then
           Output%Smoke_Mask = CLAVRX_SMOKE_TEST(Input%Ref_063um,Input%Ref_063um_Clear, &
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
        Output%Dust_Mask = 0
        if (Input%SNOW_Class == symbol%NO_SNOW .and. Input%Chan_On_12um == symbol%YES .and. &
            (Input%Land_Class == symbol%DEEP_OCEAN .or. Input%Land_Class == symbol%MODERATE_OCEAN)) then
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
        !-------------------------------------------------------------------------
        !-- Keiko Yamamoto Dust algorithm
        !-------------------------------------------------------------------------
        !Output%Dust_Mask = IR_DUST_TEST2 (     &
        !            Input%Land_Class,          &
        !            Input%Bt_85um,             &
        !            Input%Bt_10um,             &
        !            Input%Bt_11um,             &
        !            Input%Bt_12um,             &
        !            Input%Bt_11um_Std,         &
        !            Input%Senzen)
!
      !-----------------------------------------------------------------------------
      ! EUMETSAT Fire Algorithm
      !-----------------------------------------------------------------------------
      if (Input%Chan_On_11um == symbol%YES .and.   &
          Input%Chan_On_375um == symbol%YES .and.  &
          Input%Land_Class == Symbol%Land) then

          if (Input%Chan_On_I4_374um == symbol%YES .and. Input%Chan_On_I5_114um == symbol%YES ) then

             Output%Fire_Mask = EUMETSAT_FIRE_TEST ( &
                   Input%Bt_11um &
                 , Input%Bt_375um &
                 , Input%Bt_I5_114um_Std &
                 , Input%Bt_I4_374um_Std &
                 , Input%Solzen )

          elseif (Input%Chan_On_I4_374um == symbol%NO .and. Input%Chan_On_I5_114um == symbol%NO ) then

             Output%Fire_Mask = EUMETSAT_FIRE_TEST ( &
                   Input%Bt_11um &
                 , Input%Bt_375um &
                 , Input%Bt_11um_Std &
                 , Input%Bt_375um_Std &
                 , Input%Solzen )

          endif

      endif

    endif ! if Invalid_Data_Mask = NO

 end subroutine NB_CLOUD_MASK_ADDONS_ALGORITHM

!----------------------------------------------------------------------------
! CLAVR-x SMOKE TEST
!
! Daytime and Ice-Free Ocean
!----------------------------------------------------------------------------
  integer elemental function CLAVRX_SMOKE_TEST(Refl_065,Refl_065_Clear,Refl_138, &
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
     
     Clavrx_Smoke_Test = MISSING_VALUE_INT1

     if (Solzen < Solzen_Max_Smoke_Thresh) then

           IR_Flag = 1
           VIS_Flag = 1
           NIR_Flag = 1
           NIR_IR_Flag = 1
           SPLIT_WIN_Flag = 1

           !--- IR test - smoke should be nearly invisible
           if (Chan_On_11 == 1) then
               if (Emiss_11_Tropo > Emiss_11_Tropo_Max_Smoke_Thresh) IR_Flag = 0
               if (T11_Std > T11_Std_Max_Smoke_Thresh) IR_Flag = 0
           endif

           !--- VIS test - smoke should be nearly invisible
           if (Chan_On_065 == 1) then
               if (Refl_065 - Refl_065_Clear > Refl_065_Max_Smoke_Thresh .or.   &
                   Refl_065 - Refl_065_Clear < Refl_065_Min_Smoke_Thresh) VIS_Flag = 0
               if (Refl_065_Std  > Refl_065_Std_Max_Smoke_Thresh) VIS_Flag = 0
           endif

           !--- NIR Tests
           if (Chan_On_375 == 1) then
               if (Refl_375 - Refl_375_Clear > Refl_375_Max_Smoke_Thresh) NIR_Flag = 0
               if (Refl_375 > Refl_375_Max_Smoke_Thresh) NIR_Flag = 0
           endif

           if (Chan_On_160 == 1 .and. Chan_On_375 == 0) then
               if (Refl_160 > Refl_160_Max_Smoke_Thresh) NIR_Flag = 0
           endif

           if (Chan_On_138 == 1) then
               if (Refl_138 > Refl_138_Max_Smoke_Thresh) NIR_Flag = 0
           endif

           !--- NIR IR_Tests
           if (Chan_On_375 == 1 .and. Chan_On_11 == 1) then
               if ((Bt_11 - Bt_375) > Btd_4_11_Max_Smoke_Thresh) NIR_IR_Flag = 0
           endif

           !--- SPLIT_WIN_Tests
           if (Chan_On_11 == 1 .and. Chan_On_12 == 1) then
             if (abs(split_window_test(Bt_11_Clear, Bt_12_Clear,Bt_11, Bt_12)) > 1) SPLIT_WIN_Flag = 0
           endif

           !--- combine into final answer
           Clavrx_Smoke_Test = IR_Flag * VIS_Flag * NIR_Flag * NIR_IR_Flag * SPLIT_WIN_Flag
   
    endif

  end function CLAVRX_SMOKE_TEST

!----------------------------------------------------------------------------
! EUMETCAST Fire detection algorithm
!
! Reference: This implements the "Current Operational Algorithm" described in:
!TOWARDS AN IMPROVED ACTIVE FIRE MONITORING PRODUCT FOR MSG SATELLITES
!Sauli Joro, Olivier Samain, Ahmet Yildirim, Leo van de Berg, Hans Joachim Lutz
!EUMETSAT, Am Kavalleriesand 31, Darmstadt, Germany
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


   !---------------------------------------------------------------------------
   ! Function calculates Dust Mask using IR channels
   ! code is based on Keiko Yamamoto
   ! (Denis B. - February, 2016)
   !---------------------------------------------------------------------------
   integer elemental function IR_DUST_TEST2 ( &
                    Land_Class,       &
                    Bt85,             &
                    Bt10,             &
                    Bt11,             &
                    Bt12,             &
                    Bt11_Std,         &
                    Senzen)

   integer(kind=int1) , intent(in) :: Land_Class
   real , intent(in) :: Senzen
   real , intent(in) :: Bt85
   real , intent(in) :: Bt10
   real , intent(in) :: Bt11
   real , intent(in) :: Bt12
   real , intent(in) :: Bt11_Std
   real, parameter :: MISSING = -999.0
   real, parameter :: BT11_STD_THRESH = 1.0
   real, parameter :: BTD_12_11_THRESH = -0.5
   real, parameter :: BTD_11_85_MIN_THRESH = -1.5
   real, parameter :: BTD_11_85_MAX_THRESH = 1.0
   real, parameter :: BT85_THRESH = 243.0
   real, parameter :: R_LAND_THRESH = -0.1
   real, parameter :: G1_LAND_MIN_THRESH = -1.0
   real, parameter :: G1_LAND_MAX_THRESH = 3.5
   real, parameter :: G2_LAND_THRESH = -0.5
   real, parameter :: B1_LAND_THRESH = 243.0
   real, parameter :: B2_LAND_THRESH = 0.997
   real, parameter :: R_OCEAN_CLOUD_THRESH = 0.0
   real, parameter :: G1_OCEAN_CLOUD_THRESH = 1.5
   real, parameter :: G2_OCEAN_CLOUD_MIN_THRESH = -1.5
   real, parameter :: G2_OCEAN_CLOUD_MAX_THRESH = 0.8
   real, parameter :: B1_OCEAN_CLOUD_THRESH = 244.0
   real, parameter :: B2_OCEAN_CLOUD_THRESH = 1.0
   real, parameter :: R_OCEAN_DUST_MIN_THRESH = -1.2
   real, parameter :: R_OCEAN_DUST_MAX_THRESH = 0.0
   real, parameter :: G1_OCEAN_DUST_THRESH = 1.2
   real, parameter :: G2_OCEAN_DUST_MIN_THRESH = -0.5
   real, parameter :: G2_OCEAN_DUST_MAX_THRESH = 0.8
   real, parameter :: B1_OCEAN_DUST_THRESH = 244.0
   real, parameter :: B2_OCEAN_DUST_THRESH = 0.997
   real, parameter :: R_OCEAN_LAND_THRESH = 0.5
   real, parameter :: G_OCEAN_LAND_THRESH = 0.0
   real, parameter :: B_OCEAN_LAND_THRESH = 0.997
   real, parameter :: SENZEN_THRESH = 76.0
   real :: Btd_12_11
   real :: Btd_85_12
   real :: Btd_11_85
   real :: Bt85_11_Ratio
   real :: Bt_Ratio
   real :: MR, MG, MB

   ! --- initialize output
   Ir_Dust_Test2 = 0

   ! --- check if valid or Bt11_Std is too big
   !     or sensor zenith angle too big
   if (Bt11 ==  Missing_Value_Real4 .or. &
       Bt11_Std > BT11_STD_THRESH .or. &
       Senzen > SENZEN_THRESH) then
      return
   endif

   ! --- calculate needed variables
   Btd_12_11 = Bt12 - Bt11
   Btd_85_12 = Bt85 - Bt12
   Btd_11_85 = Bt11 - Bt85
   Bt85_11_Ratio = Bt85 / Bt11
   ! in case 10 micron is on
   if (Bt10 /= Missing_Value_Real4 .and. Btd_85_12 .ne. 0.0) then
      Bt_Ratio = (Bt10 - Bt11) / Btd_85_12
   endif

   ! --- 1st RGB check
   if (Btd_12_11 < BTD_12_11_THRESH .or. &
       Btd_11_85 < BTD_11_85_MIN_THRESH .or. &
       Btd_11_85 > BTD_11_85_MAX_THRESH .or. &
       Bt85 < BT85_THRESH) then
      return
   endif

   ! --- 2nd RGB check depends on surface
   ! Over Land Tests
   if (Land_Class .ge. 1 .and. Land_Class .le. 5) then
      MR = 1
      MG = 1
      MB = 1
      if (Btd_12_11 < R_LAND_THRESH) MR = 0
      ! in case 10 micron is on
      if (Bt10 /= Missing_Value_Real4) then
         if (Btd_11_85 > G1_LAND_MIN_THRESH .and. &
             Btd_11_85 < G1_LAND_MAX_THRESH .and. &
             Bt_Ratio < G2_LAND_THRESH) MG = 0
      ! in case 10 micron is off
      else
         if (Btd_11_85 > G1_LAND_MIN_THRESH .and. &
             Btd_11_85 < G1_LAND_MAX_THRESH) MG = 0
      endif
      if (Bt85 < B1_LAND_THRESH .and. &
          Bt85_11_Ratio > B2_LAND_THRESH) MB = 0
      ! decision
      if ((MR * MG * MB) /= 0) Ir_Dust_Test2 = 1
   endif

   ! Over Water Tests
   if (Land_Class == 0 .or. Land_Class > 5) then
      ! cloud like areas
      MR = 1
      MG = 1
      MB = 1
      if (Btd_12_11 < R_OCEAN_CLOUD_THRESH) MR = 0
      ! in case 10 micron is on
      if (Bt10 /= Missing_Value_Real4) then
         if (Btd_11_85 < G1_OCEAN_CLOUD_THRESH .and. &
             Bt_Ratio > G2_OCEAN_CLOUD_MIN_THRESH .and. &
             Bt_Ratio < G2_OCEAN_CLOUD_MAX_THRESH) MG = 0
      ! in case 10 micron is off
      else
         if (Btd_11_85 < G1_OCEAN_CLOUD_THRESH) MG = 0
      endif
      if (Bt85 < B1_OCEAN_CLOUD_THRESH .and. &
          Bt85_11_Ratio < B2_OCEAN_CLOUD_THRESH) MB = 0
      ! decision
      if (((MR + MG) * MB) /= 0) Ir_Dust_Test2 = 1

      ! dust like areas
      MR = 1
      MG = 1
      MB = 1
      if (Btd_12_11 > R_OCEAN_DUST_MIN_THRESH .and. &
          Btd_12_11 < R_OCEAN_DUST_MAX_THRESH) MR = 0
            ! in case 10 micron is on
      if (Bt10 .ne. Missing_Value_Real4) then
         if (Btd_11_85 < G1_OCEAN_DUST_THRESH .and. &
             Bt_Ratio > G2_OCEAN_DUST_MIN_THRESH .and. &
             Bt_Ratio < G2_OCEAN_DUST_MAX_THRESH) MG = 0
      ! in case 10 micron is off
      else
         if (Btd_11_85 < G1_OCEAN_DUST_THRESH) MG = 0
      endif
      if (Bt85 < B1_OCEAN_DUST_THRESH .and. &
          Bt85_11_Ratio > B2_OCEAN_DUST_THRESH) MB = 0
      ! decision
      if (((MR + MG) * MB) /= 0) Ir_Dust_Test2 = 1
   
      ! land slide like areas
      MR = 1
      MG = 1
      MB = 1
      if (Btd_11_85 > R_OCEAN_LAND_THRESH) MR = 0 
      if (Bt10 /= Missing_Value_Real4) then
         if (Bt_Ratio < G_OCEAN_LAND_THRESH) MG = 0
      ! in case 10 micron is off
      else
         MG = 0
      endif
      if (Bt85_11_Ratio > B_OCEAN_LAND_THRESH) MB = 0
      ! decision
      if ((MR + MG + MB) /= 0) Ir_Dust_Test2 = 1
   endif   
  
   end function IR_DUST_TEST2
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

!-----------------------------------------------------------
! end of module
!-----------------------------------------------------------
end module CLOUD_MASK_ADDONS
