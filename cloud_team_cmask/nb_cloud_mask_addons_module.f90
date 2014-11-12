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
           VCM_SMOKE_TEST, &
           VCM_DUST_TEST, &
           IR_DUST_TEST, &
           MEDIAN_FILTER

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

  !internal variables
!  integer :: Line_Start
!  integer :: Line_End
!  integer :: Elem_Idx
!  integer :: Line_Idx
!!  integer:: Num_Elem
!  integer:: Num_Line
!  integer:: Num_Line_Max
  integer:: i1, i2, j1, j2
  integer:: N_Median
  real, dimension(:,:), allocatable:: Median_Input, Median_Output
  integer(kind=int1), dimension(:,:), allocatable:: Median_Mask

  !------------------------------------------------------------------------------------------
  !---  begin executable code
  !------------------------------------------------------------------------------------------
!  Num_Elem = Input%Num_Elem
!  Num_Line = Input%Num_Line
!  Num_Line_Max = Input%Num_Line_Max

  !--- determine ending line number (Line_Start can differ from 1)
!  Line_Start = 1
!  Line_End = Line_Start + Num_Line - 1

  !--- initialize to missing
  Output%Dust_Mask = MISSING_VALUE_INT1
  Output%Smoke_Mask = MISSING_VALUE_INT1
  Output%Fire_Mask = MISSING_VALUE_INT1

  !------------------------------------------------------------------------------------------
  !  MAIN LOOP
  !------------------------------------------------------------------------------------------
!  line_loop:   do Line_Idx = Line_Start, Line_End
!   elem_loop:   do Elem_Idx = 1, Num_Elem

   !--- check for valid data
   if (Input%Invalid_Data_Mask == Symbol%NO) then

    !-------------------------------------------------------------------------
    ! Dust Detection 
    !-------------------------------------------------------------------------

        !-------------------------------------------------------------------------
        ! CLAVR-x IR Dust Detection 
        !-------------------------------------------------------------------------
!        if (Input%Chan_On_11um == symbol%YES  .and. Input%Chan_On_12um == symbol%YES) then
!
!             Output%Dust_Mask = IR_DUST_TEST( &
!                    Input%Solzen,  &
!                    Input%Emiss_375um,  &
!                    Input%Bt_11um,  &
!                    Input%Bt_12um,  &
!                    Input%Emiss_375um_Clear,  &
!                    Input%Bt_11um_Clear,  &
!                    Input%Bt_12um_Clear,  &
!                    Input%Bt_11um_Std)
!                
!
!           Diag%Array_1 = Input%Bt_11um - Input%Bt_12um
!
!          Diag%Array_2 = (Input%Bt_11um_Clear - Input%Bt_12um_Clear)* &
!                                           (Input%Bt_11um-260.0) / (Input%Bt_11um_Clear-260.0)
!          if (Input%Bt_11um < 260.0) Diag%Array_2 = 0.0
!          Diag%Array_3 = Diag%Array_1 - Diag%Array_2
!
!          Diag%Array_2 = (Input%Bt_11um_Clear - Input%Bt_12um_Clear)
!
!          !Diag%Array_2 = Input%Ref_375um / Input%Ref_063um
!
!           Diag%Array_2 = Input%Solzen
!
!           Diag%Array_3 = abs(Input%Emiss_375um - Input%Emiss_375um_Clear)
!
!        endif
   
     
       !--------------------------------------------------------------------------
       !  VCM Daytime Dust Detection
       !--------------------------------------------------------------------------
       if ( Input%Solzen <= 85.0 .and. &
            Input%Chan_On_063um == symbol%YES  .and. &
            Input%Chan_On_041um == symbol%YES) then


          if (Input%Chan_On_I1_064um == symbol%YES ) then

               Output%Dust_Mask = VCM_DUST_TEST ( &
                   Input%Ref_041um &
                 , Input%Ref_063um &
                 , Input%Ref_I1_064um_Std &
                 , Input%Oceanic_Glint_Mask &
                 , Input%Land_Class )

          elseif (Input%Chan_On_I1_064um == symbol%NO ) then

               Output%Dust_Mask = VCM_DUST_TEST ( &
                   Input%Ref_041um &
                 , Input%Ref_063um &
                 , Input%Ref_063um_Std &
                 , Input%Oceanic_Glint_Mask &
                 , Input%Land_Class )

          endif

      endif
       !---------------------------------------------------------------------
       ! VCM Daytime Smoke Detection
       !---------------------------------------------------------------------
       if ( Input%Solzen <= 85.0 .and. &
            Input%Chan_On_041um == symbol%YES .and.  &
            Input%Chan_On_213um == symbol%YES) then

            if (Input%Chan_On_I1_064um == symbol%YES ) then

                Output%Smoke_Mask = VCM_SMOKE_TEST ( &
                    Input%Ref_041um &
                  , Input%Ref_213um &
                  , Input%Ref_I1_064um_Std &
                  , Input%Senzen &
                  , Input%Oceanic_Glint_Mask &
                  , Input%Land_Class )

            else if (Input%Chan_On_063um == symbol%YES .and. &
                     Input%Chan_On_I1_064um == symbol%NO ) then

                Output%Smoke_Mask = VCM_SMOKE_TEST ( &
                    Input%Ref_041um &
                  , Input%Ref_213um &
                  , Input%Ref_063um_Std &
                  , Input%Senzen &
                  , Input%Oceanic_Glint_Mask &
                  , Input%Land_Class )

            endif

      endif
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
!   enddo elem_loop
!  enddo line_loop


!!!!!!!! Disabled until Andy's explanation (Denis B. 11/10/2014) !!!!!!!!!
  !----------------------------------------------------------------------------
  ! Apply Median Filter to Smoke and Dust
  !----------------------------------------------------------------------------
!   N_Median = 4
!   allocate(Median_Input(Num_Elem,Num_Line),Median_Mask(Num_Elem,Num_Line_Max),Median_Output(Num_Elem,Num_Line_Max))
!   Median_Input = real(Output%Dust_Mask)
!   Median_Mask = Input%Invalid_Data_Mask
!
!   do Line_Idx = Line_Start, Line_End
!     do Elem_Idx = 1, Num_Elem
!
!     j1 = max(Line_Start,Line_Idx-N_Median) !top index of local array
!     j2 = min(Line_End,Line_Idx+N_Median)   !bottom index of local array
!     i1 = max(1,Elem_Idx-N_Median)          !left index of local array
!     i2 = min(Num_Elem,Elem_Idx+N_Median)   !right index of local array
!
!     !--- compute median
!     call MEDIAN_FILTER(Median_Input(i1:i2,j1:j2),Median_Mask(i1:i2,j1:j2),Median_Output(Elem_Idx,Line_Idx))
!!    print *, "median filter test ", Elem_Idx, Line_Idx, i1, i2, j1,j2
!!    print *, "input = ", Median_Input(i1:i2,j1:j2)
!!    print *, "mask = ", Median_Mask(i1:i2,j1:j2)
!!    print *, "Output = ", Median_Output(Elem_Idx,Line_Idx)
!     enddo
!  enddo
!  Output%Dust_Mask = nint(Median_Output,kind=int1)
!
!  deallocate(Median_Input, Median_Mask, Median_Output)

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
  !----------------------------------------------------------------------------
  ! VIIRS VCM SMOKE TEST
  !
  !  Reference: ???
  !----------------------------------------------------------------------------
  integer elemental function VCM_SMOKE_TEST ( &
           Ref_004 &
         , Ref_021 &
         , Ref_006_Std &
         , sat_zen &
         , Is_glint &
         , Land_Class )

      real , intent(in) :: Ref_004
      real , intent(in) :: Ref_021
      real , intent(in) :: Ref_006_Std
      real , intent(in) :: sat_zen
      integer(kind=int1) , intent(in) :: Is_glint
      integer(kind=int1) , intent(in) :: Land_Class

      logical :: Is_water_sfc

      real, parameter :: pi = 3.14159265359
      real, parameter :: SMOKE_STD_DEV_LAND_THRESH = 0.1
      real, parameter :: SMOKE_STD_DEV_WATER_THRESH = 0.05
      real, parameter :: SMOKE_STD_DEV_LAND_GLINT_THRESH = 1.5
      real, parameter :: SMOKE_STD_DEV_WATER_GLINT_THRESH = 0.05
      real, parameter :: SMOKE_CAND_M11M1_REF_RATIO_THRESH = 0.25

      real :: Ref_ratio                                                                                                                                            
      real :: Std_Dev_thresh_cand

      Vcm_smoke_Test = 0

      ! - check if valid
      if ( Ref_004 < 0. .or. Ref_021 < 0. .or. Ref_006_Std < 0.) then
         return
      end if

      Is_water_sfc = Land_Class == 0 .or. &
         Land_Class >= 3 .and. Land_Class <= 7

      Ref_ratio = Ref_021 / Ref_004

      if ( Is_water_sfc .and. Ref_ratio > ( SMOKE_CAND_M11M1_REF_RATIO_THRESH &
                             * cos (sat_zen*pi/180.0) ) ) return

      if ( Is_glint == 1 .and. Is_water_sfc )  &
                Std_Dev_thresh_cand = SMOKE_STD_DEV_WATER_GLINT_THRESH
      if ( Is_glint == 1 .and. .not. ( Is_water_sfc ) ) &
                Std_Dev_thresh_cand = SMOKE_STD_DEV_LAND_GLINT_THRESH
      if ( Is_glint == 0 .and. Is_water_sfc ) &
                Std_Dev_thresh_cand = SMOKE_STD_DEV_WATER_THRESH
      if ( Is_glint == 0 .and. .not. ( Is_water_sfc ) ) &
                Std_Dev_thresh_cand = SMOKE_STD_DEV_LAND_THRESH

      if ( Ref_006_Std < Std_Dev_thresh_cand ) Vcm_smoke_Test = 1

  end function VCM_SMOKE_TEST

  !----------------------------------------------------------------------------
  ! VIIRS VCM DUST TEST
  !
  !  Reference: ???
  !----------------------------------------------------------------------------
  integer elemental function VCM_DUST_TEST ( &
              Ref_004 &
            , Ref_006 &
            , Ref_006_Std &
            , Is_glint &
            , Land_Class )

      real , intent(in) :: Ref_004
      real , intent(in) :: Ref_006
      real , intent(in) :: Ref_006_Std
      integer(kind=int1) , intent(in) :: Is_glint
      integer(kind=int1) , intent(in) :: Land_Class

      logical :: Is_water_sfc

      real, parameter :: DUST_M1_REFL_THRESH = 0.8
      real, parameter :: DUST_CAND_M1M5_REFL_RATIO_THRESH = 0.25
      real, parameter :: DUST_STD_DEV_LAND_THRESH = 0.1
      real, parameter :: DUST_STD_DEV_WATER_THRESH = 0.05
      real, parameter :: DUST_STD_DEV_LAND_GLINT_THRESH = 0.7
      real, parameter :: DUST_STD_DEV_WATER_GLINT_THRESH = 0.05

      real :: Std_Dev_thresh_cand

      Vcm_Dust_Test = 0

      ! - check if valid
      if ( Ref_004 < 0.0 .or. Ref_006 < 0.0 .or. Ref_006_Std <= 0.) then    
         return
      end if

      ! -- exclude water 
      Is_water_sfc = Land_Class == 0 .or. Land_Class >= 3 .and. Land_Class <= 7


      if ( Is_water_sfc .and. ( Ref_004 >= DUST_M1_REFL_THRESH &
        .or. Ref_004 / Ref_006 >= DUST_CAND_M1M5_REFL_RATIO_THRESH )) then
         return
      end if

      ! adjust thresholds
      if ( Is_Glint == 1 .and. Is_water_sfc ) &
           Std_Dev_thresh_cand = DUST_STD_DEV_WATER_GLINT_THRESH
      if ( Is_Glint == 1 .and. .not. (Is_water_sfc) ) &
           Std_Dev_thresh_cand = DUST_STD_DEV_LAND_GLINT_THRESH
      if ( Is_Glint == 0 .and. Is_water_sfc ) &
           Std_Dev_thresh_cand = DUST_STD_DEV_WATER_THRESH
      if ( Is_glint == 0 .and. .not. (Is_water_sfc) ) &
           Std_Dev_thresh_cand = DUST_STD_DEV_LAND_THRESH

      if ( Ref_006_Std < Std_Dev_thresh_cand ) Vcm_Dust_Test = 1

   end function VCM_DUST_TEST

   !---------------------------------------------------------------------------
   !
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

   !--- check for cloud
!  if (Solzen > 90.0) then
!    if (abs((Emiss_375um - Emiss_375um)) > 0.2) then
!     Ir_Dust_Test = 0
!    endif
!  endif 

   if (Solzen > 90.0) then
     if (abs(Emiss_375um - Emiss_375um_Clear) > 0.1) then
      Ir_Dust_Test = 0
     endif
   endif 

   end function IR_DUST_TEST
!==============================================================
! Median filter
!
! mask = 0 means use median, mask = 1 mean ignore
!==============================================================
subroutine MEDIAN_FILTER(z,mask,z_median)

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
      if (mask(i,j) == 0 .and. z(i,j) /= missing_value_real4) then
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

end subroutine MEDIAN_FILTER
!-----------------------------------------------------------
! end of module
!-----------------------------------------------------------
end module CLOUD_MASK_ADDONS
