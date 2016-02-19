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
           IR_DUST_TEST, &
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
        !-------------------------------------------------------------------------
        Output%Smoke_Mask = 0

        !-------------------------------------------------------------------------
        !-- IR Dust algorithm
        !-------------------------------------------------------------------------
        Output%Dust_Mask = IR_DUST_TEST2 (     &
                    Input%Land_Class,          &
                    Input%Bt_85um,             &
                    Input%Bt_10um,             &
                    Input%Bt_11um,             &
                    Input%Bt_12um,             &
                    Input%Bt_11um_Std,         &
                    Input%Senzen)

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
! Reference: Baum and Trepte, 1998
!
!----------------------------------------------------------------------------

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
   integer elemental function IR_DUST_TEST ( &
                    Solzen,             &
                    Emiss_375um,        &
                    Bt_11um,            &
                    Bt_12um,            &
                    Emiss_375um_Clear,  &
                    Bt_11um_Clear,      &
                    Bt_12um_Clear,      &
                    Bt_11um_Std)

   real , intent(in) :: Solzen
   real , intent(in) :: Emiss_375um
   real , intent(in) :: Bt_11um
   real , intent(in) :: Bt_12um
   real , intent(in) :: Emiss_375um_Clear
   real , intent(in) :: Bt_11um_Clear
   real , intent(in) :: Bt_12um_Clear
   real , intent(in) :: Bt_11um_Std
   real :: Btd_11_12
   real :: Btd_11_12_Metric
   real , parameter:: Btd_11_12_Metric_Max_Thresh = -1.0
   real , parameter:: Btd_11_12_Max_Thresh = -0.25 !-1.0
   real , parameter:: Bt_11_Std_Max_Thresh = 3.0

   Ir_Dust_Test = 0

   Btd_11_12 = (Bt_11um - Bt_12um)
   Btd_11_12_Metric = (Bt_11um_Clear - Bt_12um_Clear)*(Bt_11um-260.0) / (Bt_11um_Clear-260.0)
   if (Bt_11um < 260.0) Btd_11_12_Metric = 0.0
   Btd_11_12_Metric = Btd_11_12 - Btd_11_12_Metric

   if (Btd_11_12 < Btd_11_12_Max_Thresh .and. Btd_11_12_Metric < Btd_11_12_Metric_Max_Thresh) then
      Ir_Dust_Test = 1
   endif

   !--- turn off if btd is not too negative but it is spatially variable in 11um
   if (Bt_11um_Std > Bt_11_Std_Max_Thresh .and. Btd_11_12 > -1.5) then
      Ir_Dust_Test = 0
   endif

   if (Solzen > 90.0) then
     if (abs(Emiss_375um - Emiss_375um_Clear) > 0.1) then
      Ir_Dust_Test = 0
     endif
   endif 

   end function IR_DUST_TEST


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

!-----------------------------------------------------------
! end of module
!-----------------------------------------------------------
end module CLOUD_MASK_ADDONS
