!$Id:$
!----------------------------------------------------------------------
! MODULE name: NAIVE_BAYESIAN_CLOUD_MASK
! 
! Routines for the determination of the Naive Bayesian cloud mask
!
! Author: Andrew Heidinger, NOAA/NESDIS
!
! DEPENDENCIES:
!      CONSTANTS
!      PIXEL_COMMON
!      NUMERICAL_ROUTINES
!      FILE_UTILITY
!      CLOUD_MASK_ROUTINES
!
! SIDE EFFECTS: None
!
! Cld_Flags Format
!
! Flag Number  bit-depth    bits      byte name
! 1            1            1         1    cloud mask attempted       
! 2            1            2         1    day_063     
! 3            1            3         1    day_063_spatial     
! 4            1            4         1    day_375
! 5            1            5         1    night_375
! 6            1            6         1    solar contamination 
! 7            1            7         1    coast
! 8            1            8         1    mountain
!----
! 9            1            9         2    forward scattering
! 10           1            10        2    cold scene 375
! 11           1            11        2    cold scene btd
! 12           1            12        2    blank
! 13           1            13        2    smoke detected
! 14           1            14        2    dust  detected
! 15           1            15        2    shadow detected
! 16           1            16        2    fire detected
!---
! 17           3            17-19     3    nbvcm surface type
! 18           1            20        3    blank   
!<-------------------- START OF CLOUD TESTS -------------------------->
! 19           2            21-22     3    T_11           (TGCT)
! 20           2            23-24     3    T_max-T        (RTCT)
! ---
! 21           2            25-26     4    T_std          (TUT)
! 22           2            27-28     4    Emiss_tropo
! 23           2            29-30     4    FMFT mask (Btd_11_12)
! 24           2            31-32     4    Btd_11_67 
!---
! 25           2            33-34     5    Bt_11_67_Covar 
! 26           2            35-36     5    Btd_11_85_Covar
! 27           2            37-38     5    Emiss_375
! 28           2            39-40     5    Emiss_375_Day
!---
! 29           2            41-42     6    Emiss_375_Night
! 30           2            43-44     6    Btd_375_11_Night
! 31           2            45-46     6    Ref_063_Day (RGCT)
! 32           2            47-48     6    Ref_std      (RUT)
!---
! 33           2            49-50     7    Ref_063_Min_3x3 (RVCT)
! 34           2            51-52     7    Ref_Ratio_Day   (RRCT)
! 35           2            53-54     7    Ref_138_Day     (CIRREF)
! 36           2            55-56     7    Ref_160_Day     (NIRREF)
!
!----------------------------------------------------------------------
module NAIVE_BAYESIAN_CLOUD_MASK

 implicit none
 private
 private:: COMPUTE_BAYES_SFC_TYPE
 private:: READ_NAIVE_BAYES
 private:: SPLIT_WINDOW_TEST
 private:: reflectance_gross_contrast_test
 private:: relative_visible_contrast_test
 private:: reflectance_ratio_test
 private:: nir_reflectance_gross_contrast_test
 private:: emiss_375um_test
 private:: emiss_375um_day_test
 private:: emiss_375um_night_test
 private:: PACK_BITS_INTO_BYTES
 private:: fire_test
 public:: CLOUD_MASK_NAIVE_BAYES
 public:: SET_CLOUD_MASK_VERSION
 public:: SET_CLOUD_MASK_THRESHOLDS_VERSION

 !--- parameters
 integer, parameter, private:: int1 = selected_int_kind(1)
 integer, parameter, private:: int2 = selected_int_kind(3)
 integer, parameter, private:: int4 = selected_int_kind(8)
 integer, parameter, private:: int8 = selected_int_kind(10)
 integer, parameter, private:: real4 = selected_real_kind(6,37)
 integer, parameter, private:: real8 = selected_real_kind(15,307)
 integer, parameter, private:: ipre = real4

 !--- string to hold version id
 character(120), private, save :: Cloud_Mask_Thresholds_Version

 !--- string to control on-screen prompts
 character(*), parameter, private :: EXE_PROMPT_CM = "Naive Bayesian Cloud Mask>> "

 !--------------------------------------------------------------------------------
 ! store table values as module-wide arrays
 !--------------------------------------------------------------------------------
 integer, private, save:: N_class
 integer, private, save:: N_sfc_bayes 
 integer, private, save:: N_bounds
 real, dimension(:), allocatable, save:: Prior_Yes
 real, dimension(:), allocatable, save:: Prior_No
 real, dimension(:), allocatable, save:: Optimal_Posterior_Prob
 integer, dimension(:), allocatable, save:: Skip_Sfc_Type_Flag
 real, dimension(:,:), allocatable, private, save:: Classifier_Bounds_Min
 real, dimension(:,:), allocatable, private, save:: Classifier_Bounds_Max
 real, dimension(:,:), allocatable, private, save:: Delta_Classifier_Bounds
 integer, dimension(:,:), allocatable, private, save:: First_valid_Classifier_Bounds
 integer, dimension(:,:), allocatable, private, save:: Last_valid_Classifier_Bounds
 real, dimension(:,:,:), allocatable, private, save:: Class_Cond_Yes
 real, dimension(:,:,:), allocatable, private, save:: Class_Cond_No 
 character (len=20), dimension(:,:), allocatable, private, save:: Classifier_Value_Name
 integer, dimension(:), allocatable,private,save:: Class_To_Test_Idx

 real, dimension(:), allocatable, private, save:: Cond_Yes
 real, dimension(:), allocatable, private, save:: Cond_No
 real, dimension(:), allocatable, private, save:: Posterior_Cld_Probability_By_Class
 real, dimension(:), allocatable, private, save:: Classifier_Value

 !--- define strucutres
 include 'structure_definitions.inc'

 !--- set thresholds and algorithm specific constants
 include 'naive_bayesian_cloud_mask.inc'

 contains
!====================================================================
!  record cvs version as a global variable for output to hdf
!====================================================================
 subroutine SET_CLOUD_MASK_VERSION(Cloud_Mask_Version)
   character(len=*), intent(out):: Cloud_Mask_Version
   Cloud_Mask_Version = "$Id: naive_bayesian_cloud_mask_module.f90 86 2014-02-28 23:24:12Z heidinger $"
 end subroutine SET_CLOUD_MASK_VERSION

!====================================================================
!  pass threshold version to bridge
!  Cloud_Mask_Thresholds_Version is a module-wide variable
!  that is private
!====================================================================
 subroutine SET_CLOUD_MASK_THRESHOLDS_VERSION(Cloud_Mask_CVS_Tag)
   character(len=*), intent(out):: Cloud_Mask_CVS_Tag
   Cloud_Mask_CVS_Tag = Cloud_Mask_Thresholds_Version
 end subroutine SET_CLOUD_MASK_THRESHOLDS_VERSION

!====================================================================
! SUBROUTINE Name: READ_NAIVE_BAYES
!
! Function:
!   Allocate and Read in the LUTs needed for Bayesian cloud mask tables
!
!====================================================================
 subroutine READ_NAIVE_BAYES(Ancil_Data_Path, Naive_Bayes_File_Name, &
                             symbol, Cloud_Mask_Bayesian_Flag)

   character(len=*), intent(in):: Ancil_Data_Path
   character(len=*), intent(in):: Naive_Bayes_File_Name
   ! Need a method to flag things
   TYPE(symbol_naive_bayesian), intent(in) :: symbol
   integer, intent(out):: Cloud_Mask_Bayesian_Flag 
   
   !local variables
   integer:: lun
   integer:: Class_Idx
   integer:: Class_Idx_Temp
   integer:: Sfc_Idx
   integer:: Sfc_Idx_file
   integer:: ios
   integer:: ios_sum
   character(len=72):: Header
   real:: time_Diff_max
   integer*4 get_lun
   
   lun = GET_LUN()
   ios_sum = 0
   
   open(unit=lun,file=trim(Ancil_Data_Path)//"bayes/"//trim(Naive_Bayes_File_Name), &
        action="read",form="formatted",status="old",iostat=ios)
   if (ios /= 0) then
     print *, EXE_PROMPT_CM , 'ERROR: Bayesian Cloud Mask Classifier Open Failed '
     print *, EXE_PROMPT_CM , 'Bayesian Cloud Mask Turned Off'
     cloud_mask_bayesian_flag = symbol%NO
     return
   endif
   ios_sum = ios_sum + ios

   !--- skip first two lines which is a meta-data
   read (unit=lun,fmt="(a120)") Cloud_Mask_Thresholds_Version

   !--- check to see if this file has a cvs tag
   if (index(Cloud_Mask_Thresholds_Version,'$Id:') > 0) then
      print *, EXE_PROMPT_CM, "Cloud Mask Threshold Version: ", trim(Cloud_Mask_Thresholds_Version)
      read (unit=lun,fmt=*) Header
   else
      Header = Cloud_Mask_Thresholds_Version
      Cloud_Mask_Thresholds_Version = 'unknown'
   endif

   print *, EXE_PROMPT_CM, "Cloud Mask Bayesian File Header: ", trim(Header)

   read (unit=lun,fmt=*) 
   read (unit=lun,fmt=*) time_Diff_max, N_Class, N_bounds, N_sfc_bayes

   !--- allocate
   allocate(Prior_Yes(N_sfc_bayes))
   allocate(Prior_No(N_sfc_bayes))
   allocate(Optimal_Posterior_Prob(N_sfc_bayes))
   allocate(Skip_Sfc_Type_Flag(N_sfc_bayes))
   allocate(Classifier_Bounds_Min(N_class,N_sfc_bayes))
   allocate(Classifier_Bounds_Max(N_class,N_sfc_bayes))
   allocate(Delta_Classifier_Bounds(N_class,N_sfc_bayes))
   allocate(First_valid_Classifier_Bounds(N_class,N_sfc_bayes))
   allocate(Last_valid_Classifier_Bounds(N_class,N_sfc_bayes))
   allocate(Class_Cond_Yes(N_bounds-1,N_class,N_sfc_bayes))
   allocate(Class_Cond_No(N_bounds-1,N_class,N_sfc_bayes))
   allocate(Classifier_Value_Name(N_class,N_sfc_bayes))
   allocate(Class_To_Test_Idx(N_class))

   !--- initialize 
   Prior_Yes = Missing_Value_Real4
   Prior_No = Missing_Value_Real4
   Optimal_Posterior_Prob = Missing_Value_Real4
   Skip_Sfc_Type_Flag = 0
   First_valid_Classifier_Bounds = 0
   Last_valid_Classifier_Bounds = 0

   do Sfc_Idx = 1, N_sfc_bayes

   read (unit=lun,fmt=*,iostat=ios) Sfc_Idx_file, Header
   ios_sum = ios_sum + ios

   !--- check if file has the expected classifiers
   if (Sfc_Idx_file /= Sfc_Idx) then
         print *, EXE_PROMPT_CM," ERROR: Confused on Naive Bayesian Cloud Mask Classifiers"
         print *, EXE_PROMPT_CM , 'Bayesian Cloud Mask Turned Off'
         cloud_mask_bayesian_flag = symbol%NO
         return
   endif

   read(unit=lun,fmt=*,iostat=ios) Prior_Yes(Sfc_Idx), Prior_No(Sfc_Idx)
   ios_sum = ios_sum + ios

   read(unit=lun,fmt=*,iostat=ios) Optimal_Posterior_Prob(Sfc_Idx)
   ios_sum = ios_sum + ios


   !--- if data are missing, skip this surface type
   if (Optimal_Posterior_Prob(Sfc_Idx) == Missing_Value_Real4) then 
      Skip_Sfc_Type_Flag(Sfc_Idx) = symbol%YES 
   endif

   do Class_Idx = 1,N_class

            read(unit=lun,fmt=*,iostat=ios)  &
                                 Class_Idx_Temp, First_valid_Classifier_Bounds(Class_Idx,Sfc_Idx), &
                                 Last_valid_Classifier_Bounds(Class_Idx,Sfc_Idx), &
                                 Classifier_Value_Name(Class_Idx,Sfc_Idx)
            ios_sum = ios_sum + ios

            Classifier_Value_Name(Class_Idx,Sfc_Idx) = trim(Classifier_Value_Name(Class_Idx,Sfc_Idx))

            read(unit=lun,fmt=*,iostat=ios)  &
                                 Classifier_Bounds_Min(Class_Idx,Sfc_Idx), &
                                 Classifier_Bounds_Max(Class_Idx,Sfc_Idx), &
                                 Delta_Classifier_Bounds(Class_Idx,Sfc_Idx)
            ios_sum = ios_sum + ios

            read(unit=lun,fmt=*,iostat=ios) Class_Cond_Yes(:,Class_Idx,Sfc_Idx)
            ios_sum = ios_sum + ios

            read(unit=lun,fmt=*,iostat=ios) Class_Cond_No(:,Class_Idx,Sfc_Idx)
            ios_sum = ios_sum + ios

    end do

   end do 

   !--- close file
   close(unit=lun,iostat=ios)

   if (ios_sum == 0) then
           print *, EXE_PROMPT_CM, ' Bayesian Cloud Mask Data Read in Successfully'
   else
           print *, EXE_PROMPT_CM, 'ERROR:  Bayesian Cloud Mask Data Not Read in Successfully'
           print *, EXE_PROMPT_CM , 'Bayesian Cloud Mask Turned Off'
           Cloud_Mask_Bayesian_Flag = symbol%NO
   endif

 end subroutine READ_NAIVE_BAYES


!====================================================================
! SUBROUTINE Name: CLOUD_MASK_NAIVE_BAYES
!
! Function:
!   CalcuLates the bayesian cloud mask. The bayesian cloud mask is
!   determined by utilizing the following surface types:
!
! Bayesian Surface Types
! 1 - Deep_Water
! 2 - Shallow_Water
! 3 - Unfrozen_Land
! 4 - Frozen_Land
! 5 - Arctic
! 6 - Antarctic
! 7 - Desert
!
!====================================================================
 subroutine CLOUD_MASK_NAIVE_BAYES( &
            Ancil_Data_Path,  &             !path to ancillary data
            Naive_Bayes_File_Name, &        !name of nb classifer data file
            Symbol, &                       !local copy of sym structure
            Input,  &
            Output,  &
            Diag)


   character (len=120), intent(in) :: Ancil_Data_Path
   character (len=120), intent(in) :: Naive_Bayes_File_Name
   type(symbol_naive_bayesian), intent(in) :: Symbol
   type(mask_input), intent(in) :: Input
   type(mask_output), intent(out) :: Output
   type(diag_output), intent(out), Optional :: Diag

   !-- local pointers that point to global variables
   integer (kind=INT1), dimension(NUMBER_OF_FLAGS):: Cld_Flags
   integer (kind=INT1), dimension(NUMBER_OF_FLAGS):: Cld_Flag_Bit_Depth

   !internal variables
   integer :: Line_Start
   integer :: Line_End
   integer :: Elem_Idx
   integer :: Line_Idx
   integer:: Class_Idx
   integer:: Sfc_Idx
   integer:: Bin_Idx

   integer:: Num_Elem
   integer:: Num_Line
   integer:: Num_Line_Max

   integer:: Oceanic_Glint_Flag
   integer:: Lunar_Oceanic_Glint_Flag
   integer:: Coastal_Flag
   integer:: Mountain_Flag
   integer:: Day_063_Flag
   integer:: Day_063_Spatial_Flag
   integer:: Night_Lunar_Flag
   integer:: Lunar_Spatial_Flag
   integer:: Day_375_Flag
   integer:: Night_375_Flag
   integer:: Forward_Scattering_Flag
   integer:: Lunar_Forward_Scattering_Flag
   integer:: Cold_Scene_375um_Flag
   integer:: Cold_Scene_Btd_Flag
   integer:: Solar_Contam_Flag
   integer:: Smoke_Flag
   integer:: Fire_Flag
   integer:: Shadow_Flag
   integer:: Dust_Flag
   integer:: City_Flag

   integer, save:: Segment_Number_Local = 1
   real (kind=real4):: Airmass
   integer, parameter:: Spare_Value = 0
   
   !------------------------------------------------------------------------------------------
   !---  begin executable code
   !------------------------------------------------------------------------------------------
   Num_Elem = Input%Num_Elem
   Num_Line = Input%Num_Line
   Num_Line_Max = Input%Num_Line_Max
   
   !------------------------------------------------------------------------------------------
   !--- on first segment, read table
   !------------------------------------------------------------------------------------------
   if (Segment_Number_Local == 1) then
       call READ_NAIVE_BAYES(Ancil_Data_Path,Naive_Bayes_File_Name, &
                             symbol,Output%Cloud_Mask_Bayesian_Flag)

        !---set up Classifier to Test Mapping
        do Class_Idx = 1, N_Class

          select case (Classifier_Value_Name(Class_Idx,1))

                    case("T_11") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 1
                    case("T_max-T") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 2
                    case("T_std") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 3
                    case("Emiss_tropo") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 4
                    case("FMFT") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 5
                    case("Btd_11_67") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 6
                    case("Bt_11_67_Covar") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 7
                    case("Btd_11_85") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 8
                    case("Emiss_375") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 9
                    case("Emiss_375_Day") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 10
                    case("Emiss_375_Night") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 11
                    case("Btd_375_11_Night") 
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 12
                    case("Ref_063_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 13
                    case("Ref_std")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 14
                    case("Ref_063_Min_3x3_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 15
                    case("Ref_Ratio_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 16
                    case("Ref_138_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 17
                    case("Ref_160_Day")
                       Class_To_Test_Idx(Class_Idx) = NUMBER_OF_NONCLOUD_FLAGS + 18
                    case default
                       print *, "Unknown Classifier Naive Bayesian Cloud Mask, stopping"
           end select

         end do

   end if

   !--- determine ending line number (Line_Start can differ from 1)
   Line_Start = 1
   Line_End = Line_Start + Num_Line - 1

   !------------------------------------------------------------------------------------------
   !--- allocate memory for local arrays
   !------------------------------------------------------------------------------------------
   if (allocated(Cond_Yes)) deallocate(Cond_Yes)
   if (allocated(Cond_No)) deallocate(Cond_No)
   if (allocated(Classifier_Value)) deallocate(Classifier_Value)
   if (allocated(Posterior_Cld_Probability_By_Class)) deallocate(Posterior_Cld_Probability_By_Class)

   allocate(Cond_Yes(N_Class))
   allocate(Cond_No(N_Class)) 
   allocate(Posterior_Cld_Probability_By_Class(N_Class)) 
   allocate(Classifier_Value(N_Class)) 

   !--- Initialize flags
   Cld_Flags = 0
   Cld_Flag_Bit_Depth = 2     ! needs to be 2 for reserving space for missing tests

   !--- initialize diagnostic output
   if (present(Diag)) Diag%Array_1 = Missing_Value_Real4
   if (present(Diag)) Diag%Array_2 = Missing_Value_Real4
   if (present(Diag)) Diag%Array_3 = Missing_Value_Real4

   !------------------------------------------------------------------------------------------
   !  MAIN LOOP
   !------------------------------------------------------------------------------------------
   line_loop:   do Line_Idx = Line_Start, Line_End
   elem_loop:   do Elem_Idx = 1, Num_Elem

          Cld_Flags(1) = symbol%YES          ;    Cld_Flag_Bit_Depth(1) = 1
          !--- check for a bad pixel
          if (Input%Invalid_Data_Mask(Elem_Idx,Line_Idx) == symbol%YES) then
              Cld_Flags(1) = symbol%NO
              cycle
          endif

          !---  COMPUTE SURFACE TYPE
          Sfc_Idx = COMPUTE_BAYES_SFC_TYPE(Input%Land_Class(Elem_Idx,Line_Idx), &
                                           Input%Coastal_Mask(Elem_Idx,Line_Idx), &
                                           Input%Snow_Class(Elem_Idx,Line_Idx), &
                                           Input%Sfc_Type(Elem_Idx,Line_Idx), &
                                           Input%Lat(Elem_Idx,Line_Idx), &
                                           Input%Lon(Elem_Idx,Line_Idx), &
                                           Input%Sst_Anal_Uni(Elem_Idx,Line_Idx), &
                                           Input%Emiss_Sfc_375um(Elem_Idx,Line_Idx), &
                                           symbol)

          !--- check for a surface types without sufficient training for processing
          if (Skip_Sfc_Type_Flag(Sfc_Idx) == symbol%YES) then
              Cld_Flags(1) = symbol%NO
              cycle
          endif

          !--- set some flags to control processing
          Oceanic_Glint_Flag = Input%Oceanic_Glint_Mask(Elem_Idx,Line_Idx)
          Lunar_Oceanic_Glint_Flag = Input%Lunar_Oceanic_Glint_Mask(Elem_Idx,Line_Idx)
          Coastal_Flag = Input%Coastal_Mask(Elem_Idx,Line_Idx)
          Solar_Contam_Flag = Input%Solar_Contamination_Mask(Elem_Idx,Line_Idx)  

          !--- compute airmass
          AirMass = 1.0/cos(Input%Solzen(Elem_Idx,Line_Idx)*dtor) + 1.0 / cos(Input%Senzen(Elem_Idx,Line_Idx)*dtor)

          !--- set day flag for 0.63 micron reflectance gross test
          if ((Input%Solzen(Elem_Idx,Line_Idx) > Reflectance_Gross_Solzen_Thresh) .or.  &
              (AirMass > Reflectance_Gross_Airmass_Thresh)) then
              Day_063_Flag = symbol%NO
          else
              Day_063_Flag = symbol%YES
          endif

          !--- set day flag for 0.63 micron reflectance spatial test
          if (Input%Solzen(Elem_Idx,Line_Idx) > Reflectance_Spatial_Solzen_Thresh) then
              Day_063_Spatial_Flag = symbol%NO
          else
              Day_063_Spatial_Flag = symbol%YES
          endif

          Lunar_Spatial_Flag = symbol%NO
          Night_Lunar_Flag = symbol%NO
          Lunar_Forward_Scattering_Flag = symbol%NO

          if (Input%Chan_On_DNB == symbol%YES) then
             Lunar_Spatial_Flag = symbol%YES
             if (Input%Lunzen(Elem_Idx,Line_Idx) > Reflectance_Spatial_Solzen_Thresh) then
               Lunar_Spatial_Flag = symbol%NO
             endif

             Night_Lunar_Flag = symbol%YES
             if (Input%Lunzen(Elem_Idx,Line_Idx) > Reflectance_Gross_Lunzen_Thresh) then
               Night_Lunar_Flag = symbol%NO
             endif

             Lunar_Forward_Scattering_Flag = symbol%NO
             if (Input%Lunscatzen(Elem_Idx,Line_Idx) < 80.0 .and. Input%Lunzen(Elem_Idx,Line_Idx) < 95.0) then
               Lunar_Forward_Scattering_Flag = symbol%YES
             endif

          endif

          if ((Input%Solzen(Elem_Idx,Line_Idx) > Emiss_375um_Day_Solzen_Thresh) .or. &
              (AirMass > Reflectance_Gross_Airmass_Thresh)) then
              Day_375_Flag = symbol%NO
          else
              Day_375_Flag = symbol%YES
          endif

          if (Input%Solzen(Elem_Idx,Line_Idx) < Emiss_375um_Night_Solzen_Thresh) then
              Night_375_Flag = symbol%NO
          else
              Night_375_Flag = symbol%YES
          endif

          if (Sfc_Idx /= 6 .and. Input%Zsfc(Elem_Idx,Line_Idx) > 2000.0) then
              Mountain_Flag = symbol%YES
          else
              Mountain_Flag = symbol%NO
          endif

          if (Input%Scatzen(Elem_Idx,Line_Idx) < 80.0 .and. Input%Solzen(Elem_Idx,Line_Idx) < 95.0) then
              Forward_Scattering_Flag = symbol%YES
          else
              Forward_Scattering_Flag = symbol%NO
          endif

          Cold_Scene_375um_Flag = symbol%NO
          if (Input%Chan_On_375um == symbol%YES) then
           if (Input%Bt_375um(Elem_Idx,Line_Idx) < Bt_375um_Cold_Scene_Thresh .and.  &
               Input%Bt_375um(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Cold_Scene_375um_Flag = symbol%YES
           endif
          endif

          Cold_Scene_Btd_Flag = symbol%NO
          if (Input%Chan_On_11um == symbol%YES) then
           if (Input%Bt_11um(Elem_Idx,Line_Idx) < Bt_11um_Cold_Scene_Thresh .and.  &
               Input%Bt_11um(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              Cold_Scene_Btd_Flag = symbol%YES
           endif
          endif


          Smoke_Flag = symbol%NO
          Dust_Flag = symbol%NO
          Shadow_Flag = symbol%NO
          Fire_Flag = symbol%NO

          if (Input%Chan_On_11um == symbol%YES .and. Input%Chan_On_375um == symbol%YES) then
             Fire_Flag = fire_test(Input%Bt_11um(Elem_Idx,Line_Idx), &
                                   Input%Bt_375um(Elem_Idx,Line_Idx), &
                                   Input%Bt_11um_Std(Elem_Idx,Line_Idx), &
                                   Input%Bt_375um_Std(Elem_Idx,Line_Idx), &
                                   Input%Solzen(Elem_Idx,Line_Idx))
          endif

          !--- City Flag of DNB Lunar Tests
          City_Flag = symbol%NO 
          if (Input%Chan_On_DNB == symbol%YES) then
            if (Input%Rad_Lunar(Elem_Idx,Line_Idx) > Radiance_Lunar_City_Thresh) City_Flag = symbol%YES
          endif

          !----------------------------------------------------------------------------------
          !--- populate elements of Cld_Flags with processing flags
          !----------------------------------------------------------------------------------
          Cld_Flags(2) = Day_063_Flag          ;    Cld_Flag_Bit_Depth(2) = 1
          Cld_Flags(3) = Day_063_Spatial_Flag  ;    Cld_Flag_Bit_Depth(3) = 1
          Cld_Flags(4) = Day_375_Flag          ;    Cld_Flag_Bit_Depth(4) = 1
          Cld_Flags(5) = Night_375_Flag        ;    Cld_Flag_Bit_Depth(5) = 1
          Cld_Flags(6) = Solar_Contam_Flag     ;    Cld_Flag_Bit_Depth(6) = 1
          Cld_Flags(7) = Coastal_Flag          ;    Cld_Flag_Bit_Depth(7) = 1
          Cld_Flags(8) = Mountain_Flag         ;    Cld_Flag_Bit_Depth(8) = 1
          Cld_Flags(9) = Forward_Scattering_Flag;   Cld_Flag_Bit_Depth(9) = 1
          Cld_Flags(10) = Cold_Scene_375um_Flag ;   Cld_Flag_Bit_Depth(10) = 1
          Cld_Flags(11) = Cold_Scene_Btd_Flag  ;    Cld_Flag_Bit_Depth(11) = 1
          Cld_Flags(12) = Spare_Value          ;    Cld_Flag_Bit_Depth(12) = 1
          Cld_Flags(13) = Smoke_Flag           ;    Cld_Flag_Bit_Depth(13) = 1
          Cld_Flags(14) = Dust_Flag            ;    Cld_Flag_Bit_Depth(14) = 1
          Cld_Flags(15) = Shadow_Flag          ;    Cld_Flag_Bit_Depth(15) = 1
          Cld_Flags(16) = Fire_Flag            ;    Cld_Flag_Bit_Depth(16) = 1
          Cld_Flags(17) = Sfc_Idx              ;    Cld_Flag_Bit_Depth(17) = 3
          Cld_Flags(18) = Spare_Value          ;    Cld_Flag_Bit_Depth(18) = 1


          class_loop: do Class_Idx = 1, N_class

             Classifier_Value(Class_Idx) = Missing_Value_Real4
             Cond_Yes(Class_Idx) = 1.0 
             Cond_No(Class_Idx) =  1.0

             select case (Classifier_Value_Name(Class_Idx,Sfc_Idx))

                    case("T_11") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um(Elem_Idx,Line_Idx)

                    case("T_max-T") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Coastal_Flag == symbol%YES) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um_Max(Elem_Idx,Line_Idx) - Input%Bt_11um(Elem_Idx,Line_Idx)

                    case("T_std") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Coastal_Flag == symbol%YES) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um_Std(Elem_Idx,Line_Idx)

                    case("Emiss_tropo") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       Classifier_Value(Class_Idx) = Input%Emiss_11um_Tropo(Elem_Idx,Line_Idx)

                    case("FMFT") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_12um == symbol%NO) cycle
                       if (Cold_Scene_Btd_Flag == symbol%YES) cycle
                       Classifier_Value(Class_Idx) =  &
                           split_window_test(Input%Bt_11um_Clear(Elem_Idx,Line_Idx), Input%Bt_12um_Clear(Elem_Idx,Line_Idx), &
                                             Input%Bt_11um(Elem_Idx,Line_Idx), Input%Bt_12um(Elem_Idx,Line_Idx))

                    case("Btd_11_67") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_67um == symbol%NO) cycle
                       if (Cold_Scene_Btd_Flag == symbol%YES) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um(Elem_Idx,Line_Idx) - Input%Bt_67um(Elem_Idx,Line_Idx)

                    case("Bt_11_67_Covar") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_67um == symbol%NO) cycle
                       if (Cold_Scene_Btd_Flag == symbol%YES) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um_Bt_67um_Covar(Elem_Idx,Line_Idx)

                    case("Btd_11_85") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_85um == symbol%NO) cycle
                       if (Cold_Scene_Btd_Flag == symbol%YES) cycle
                       Classifier_Value(Class_Idx) = Input%Bt_11um(Elem_Idx,Line_Idx) - Input%Bt_85um(Elem_Idx,Line_Idx)

                    case("Emiss_375") 
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       if (Input%Chan_On_375um == symbol%YES) then
                          if (Input%Bt_375um(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
                       endif
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       Classifier_Value(Class_Idx) = emiss_375um_test( &
                                         Input%Emiss_375um(Elem_Idx,Line_Idx),Input%Emiss_375um_Clear(Elem_Idx,Line_Idx))

                    case("Emiss_375_Day") 
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       if (Day_375_Flag == symbol%NO) cycle
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       if (Input%Chan_On_375um == symbol%YES) then
                          if (Input%Bt_375um(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
                      endif
                      Classifier_Value(Class_Idx) = emiss_375um_day_test( &
                                         Input%Emiss_375um(Elem_Idx,Line_Idx),Input%Emiss_375um_Clear(Elem_Idx,Line_Idx))

                    case("Emiss_375_Night") 
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                       if (Night_375_Flag == symbol%NO) cycle
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       if (Input%Chan_On_375um == symbol%YES) then
                          if (Input%Bt_375um(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
                       endif
                       Classifier_Value(Class_Idx) = emiss_375um_night_test( &
                                        Input%Emiss_375um(Elem_Idx,Line_Idx),Input%Emiss_375um_Clear(Elem_Idx,Line_Idx))

                     case("Btd_375_11_Night") 
                       if (Input%Chan_On_11um == symbol%NO) cycle
                       if (Input%Chan_On_375um == symbol%NO) cycle
                       if (Solar_Contam_Flag == symbol%YES) cycle
                       if (Night_375_Flag == symbol%NO) cycle
                       if (Cold_Scene_375um_Flag == symbol%YES) cycle
                       if (Input%Chan_On_375um == symbol%YES) then
                         if (Input%Bt_375um(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
                       endif
                       Classifier_Value(Class_Idx) = Input%Bt_375um(Elem_Idx,Line_Idx) - &
                                                     Input%Bt_11um(Elem_Idx,Line_Idx)

                    case("Ref_063_Day")

                       if (Input%Solzen(Elem_Idx,Line_Idx) < 90.0) then

                         if (Input%Chan_On_063um == symbol%NO) cycle
                         if (Oceanic_Glint_Flag == symbol%YES) cycle
                         if (Forward_Scattering_Flag == symbol%YES) cycle
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Day_063_Flag == symbol%NO) cycle
                         if (Sfc_Idx == 4) cycle
                         Classifier_Value(Class_Idx) =  &
                             reflectance_gross_contrast_test(Input%Ref_063um_Clear(Elem_Idx,Line_Idx), &
                                                             Input%Ref_063um(Elem_Idx,Line_Idx))

                       else

                         if (Input%Chan_On_DNB == symbol%NO) cycle
                         if (Lunar_Oceanic_Glint_Flag == symbol%YES) cycle
                         if (Lunar_Forward_Scattering_Flag == symbol%YES) cycle
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Night_Lunar_Flag == symbol%NO) cycle
                         if (City_Flag == symbol%YES) cycle
                         if (Sfc_Idx == 4) cycle
                         Classifier_Value(Class_Idx) =  &
                             reflectance_gross_contrast_test(Input%Ref_Lunar_Clear(Elem_Idx,Line_Idx), &
                                                             Input%Ref_Lunar(Elem_Idx,Line_Idx))

                       endif

                    case("Ref_std")
                       if (Input%Solzen(Elem_Idx,Line_Idx) < 90.0) then
                         if (Input%Chan_On_063um == symbol%NO) cycle
                         if (Day_063_Spatial_Flag == symbol%NO) cycle 
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Coastal_Flag == symbol%YES) cycle
                         Classifier_Value(Class_Idx) = Input%Ref_063um_Std(Elem_Idx,Line_Idx)
                       else
                         if (Input%Chan_On_DNB == symbol%NO) cycle
                         if (Lunar_Spatial_Flag == symbol%NO) cycle
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Coastal_Flag == symbol%YES) cycle
                         if (City_Flag == symbol%YES) cycle  
                         Classifier_Value(Class_Idx) = Input%Ref_Lunar_Std(Elem_Idx,Line_Idx)
                       endif

                    case("Ref_063_Min_3x3_Day")
                       if (Input%Solzen(Elem_Idx,Line_Idx) < 90.0) then
                         if (Input%Chan_On_063um == symbol%NO) cycle
                         if (Day_063_Spatial_Flag == symbol%NO) cycle
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Coastal_Flag == symbol%YES) cycle
                         Classifier_Value(Class_Idx) = relative_visible_contrast_test( &
                                       Input%Ref_063um_Min(Elem_Idx,Line_Idx),Input%Ref_063um(Elem_Idx,Line_Idx))
                       else
                         if (Input%Chan_On_DNB == symbol%NO) cycle
                         if (Lunar_Spatial_Flag == symbol%NO) cycle
                         if (Night_Lunar_Flag == symbol%NO) cycle
                         if (Mountain_Flag == symbol%YES) cycle
                         if (Coastal_Flag == symbol%YES) cycle
                         if (City_Flag == symbol%YES) cycle  
                         Classifier_Value(Class_Idx) = relative_visible_contrast_test( &
                                       Input%Ref_Lunar_Min(Elem_Idx,Line_Idx),Input%Ref_Lunar(Elem_Idx,Line_Idx))
                       endif

                    case("Ref_Ratio_Day")
                       if (Input%Chan_On_063um == symbol%NO) cycle
                       if (Input%Chan_On_086um == symbol%NO) cycle
                       if (Day_063_Spatial_Flag == symbol%NO) cycle
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       Classifier_Value(Class_Idx) = reflectance_ratio_test( &
                                         Input%Ref_063um(Elem_Idx,Line_Idx),Input%Ref_086um(Elem_Idx,Line_Idx))

                    case("Ref_138_Day")
                       if (Input%Chan_On_138um == symbol%NO) cycle
                       if (Forward_Scattering_Flag == symbol%YES) cycle
                       if (Day_063_Flag == symbol%NO) cycle
                       if (Mountain_Flag == symbol%YES) cycle
                       if (Input%Chan_On_138um == symbol%YES) then
                         if (Input%Ref_138um(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
                       endif
                       Classifier_Value(Class_Idx) = Input%Ref_138um(Elem_Idx,Line_Idx)

                    case("Ref_160_Day")
                       if (Input%Chan_On_160um == symbol%NO) cycle
                       if (Forward_Scattering_Flag == symbol%YES) cycle
                       if (Day_063_Flag == symbol%NO) cycle
                       if (Oceanic_Glint_Flag == symbol%YES) cycle
                       if (Input%Chan_On_160um == symbol%YES) then
                         if (Input%Ref_160um(Elem_Idx,Line_Idx) == Missing_Value_Real4) cycle
                       endif
                       Classifier_Value(Class_Idx) = Input%Ref_160um(Elem_Idx,Line_Idx)
                    
                     case default
                       print *, "Unknown Classifier Naive Bayesian Cloud Mask, stopping"
 
             end select

             !--- Turn off Classifiers if Chosen Metric is Missing
             if (Classifier_Value(Class_Idx) == Missing_Value_Real4) then
                cycle
             endif


             !--- interpolate class conditional values
             Bin_Idx = int ((Classifier_Value(Class_Idx) - Classifier_Bounds_Min(Class_Idx,Sfc_Idx))  /    &
                           Delta_Classifier_Bounds(Class_Idx,Sfc_Idx)) + 1
             Bin_Idx = max(1,min(N_bounds-1,Bin_Idx))

             Cond_Yes(Class_Idx) = Class_Cond_Yes(Bin_Idx,Class_Idx,Sfc_Idx)
             Cond_No(Class_Idx) = Class_Cond_No(Bin_Idx,Class_Idx,Sfc_Idx)

        enddo  class_loop 


        !------------------------------------------------------------------------------------------------------------
        !--- compute prosterior probabilites for each pixel
        !-----------------------------------------------------------------------------------------------------------
        if (Prior_Yes(Sfc_Idx) == Missing_Value_Real4) then
          Output%Posterior_Cld_Probability(Elem_Idx,Line_Idx) = Missing_Value_Real4
          Posterior_Cld_Probability_By_Class = Missing_Value_Real4
        endif

        if (minval(Cond_Yes) == 0.0) then
          Output%Posterior_Cld_Probability(Elem_Idx,Line_Idx) = 0.0
        elseif (minval(Cond_No) == 0.0) then
          Output%Posterior_Cld_Probability(Elem_Idx,Line_Idx) = 1.0
        else
          Output%Posterior_Cld_Probability(Elem_Idx,Line_Idx) = &
                 (Prior_Yes(Sfc_Idx)*product(Cond_Yes)) / &
                 (Prior_Yes(Sfc_Idx)*product(Cond_Yes) +  &        
                 Prior_No(Sfc_Idx)*product(Cond_No))
        endif

        !------------------------------------------------------------------------------------------------------------
        !--- make a cloud mask
        !------------------------------------------------------------------------------------------------------------
        Output%Cld_Mask_Bayes(Elem_Idx,Line_Idx) = symbol%CLEAR
        if (Output%Posterior_Cld_Probability(Elem_Idx,Line_Idx) >= 0.9) then
                Output%Cld_Mask_Bayes(Elem_Idx,Line_Idx) = symbol%CLOUDY 
        endif
        if ((Output%Posterior_Cld_Probability(Elem_Idx,Line_Idx) >= 0.5) .and. &
            (Output%Posterior_Cld_Probability(Elem_Idx,Line_Idx) < 0.9)) then
                Output%Cld_Mask_Bayes(Elem_Idx,Line_Idx) = symbol%PROB_CLOUDY 
        endif
        if ((Output%Posterior_Cld_Probability(Elem_Idx,Line_Idx) > 0.1) .and. &
            (Output%Posterior_Cld_Probability(Elem_Idx,Line_Idx) < 0.5)) then
                Output%Cld_Mask_Bayes(Elem_Idx,Line_Idx) = symbol%PROB_CLEAR
        endif


        !------------------------------------------------------------------------------------------------------------
        !--- compute probabilities for each class alone - used for flags - not
        !--- needed for mask or final probability - it should remain optional
        !------------------------------------------------------------------------------------------------------------
        if (Do_By_Class_Flag == symbol%YES) then

         do Class_Idx = 1, N_class
 
          if (Cond_Yes(Class_Idx) == 0.0) then
              Posterior_Cld_Probability_By_Class(Class_Idx) = 0.0
          elseif (Cond_No(Class_Idx) == 0.0) then
              Posterior_Cld_Probability_By_Class(Class_Idx) = 1.0
          else
              Posterior_Cld_Probability_By_Class(Class_Idx) = &
                 (Prior_Yes(Sfc_Idx)*Cond_Yes(Class_Idx)) / &
                 (Prior_Yes(Sfc_Idx)*Cond_Yes(Class_Idx) +  &        
                  Prior_No(Sfc_Idx)*Cond_No(Class_Idx))

          endif

          !-- set cloud flags
          Cld_Flag_Bit_Depth(Class_To_Test_Idx(Class_Idx)) = 2
          Cld_Flags(Class_To_Test_Idx(Class_Idx)) = symbol%CLOUDY
          if (Posterior_Cld_Probability_By_Class(Class_Idx) <= 0.1) then
               Cld_Flags(Class_to_Test_Idx(Class_Idx)) = symbol%CLEAR

          elseif (Posterior_Cld_Probability_By_Class(Class_Idx) > 0.1 .and. &
                  Posterior_Cld_Probability_By_Class(Class_Idx) <= 0.5) then
               Cld_Flags(Class_to_Test_Idx(Class_Idx)) = symbol%PROB_CLEAR

          elseif (Posterior_Cld_Probability_By_Class(Class_Idx) > 0.5 .and. &
                  Posterior_Cld_Probability_By_Class(Class_Idx) <= 0.9) then
               Cld_Flags(Class_to_Test_Idx(Class_Idx)) = symbol%PROB_CLOUDY
          endif


         !-----------------------------------------------------------------------
         !--- Diagnostic Output for debugging in CLAVR-x
         !-----------------------------------------------------------------------
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Ref_std") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Ref_063_Day") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Ref_063_Min_3x3_Day") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Ref_160_Day") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Btd_11_67") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Btd_375_11_Night") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "T_11") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Emiss_375") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Emiss_375_Day") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Emiss_375_Night") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Bt_11_67_Covar") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Btd_11_85_Covar") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "Emiss_tropo") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "T_std") then
         !if (trim(Classifier_Value_Name(Class_Idx,1)) == "T_max-T") then
          if (trim(Classifier_Value_Name(Class_Idx,1)) == "FMFT") then
              if (present(Diag)) Diag%Array_1(Elem_Idx,Line_Idx) = Classifier_Value(Class_Idx)
              if (present(Diag)) Diag%Array_2(Elem_Idx,Line_Idx) = Posterior_Cld_Probability_By_Class(Class_Idx)
              if (present(Diag)) Diag%Array_3(Elem_Idx,Line_Idx) = Cld_Flags(Class_To_Test_Idx(Class_Idx))
          endif

         enddo

         !-------------------------------------------------------------------
         ! Pack Bytes
         !-------------------------------------------------------------------
         call PACK_BITS_INTO_BYTES(Cld_Flags,Cld_Flag_Bit_Depth,Output%Cld_Flags_Packed(:,Elem_Idx,Line_Idx))

        endif   !Do_By_Class_Loop


        
      enddo elem_loop
   enddo line_loop

   !--- on last segment, deallocate arrays
   if (Segment_Number_Local == Input%Num_Segments) then

       deallocate(Prior_Yes)
       deallocate(Prior_No)
       deallocate(Optimal_Posterior_Prob)
       deallocate(Skip_Sfc_Type_Flag)
       deallocate(Classifier_Bounds_Min)
       deallocate(Classifier_Bounds_Max)
       deallocate(Delta_Classifier_Bounds)
       deallocate(First_valid_Classifier_Bounds)
       deallocate(Last_valid_Classifier_Bounds)
       deallocate(Class_Cond_Yes)
       deallocate(Class_Cond_No)
       deallocate(Classifier_Value_Name)
       deallocate(Class_To_Test_Idx)

    endif

    !--- increment segment number
    if (Segment_Number_Local == Input%Num_Segments) then
       Segment_Number_Local = 1
    else
       Segment_Number_Local = Segment_Number_Local + 1
    endif

   !----- deallocate memory
   if (allocated(Cond_Yes)) deallocate(Cond_Yes)
   if (allocated(Cond_No)) deallocate(Cond_No)
   if (allocated(Classifier_Value)) deallocate(Classifier_Value)
   if (allocated(Posterior_Cld_Probability_By_Class)) deallocate(Posterior_Cld_Probability_By_Class)

 end subroutine CLOUD_MASK_NAIVE_BAYES

!====================================================================
! SUBROUTINE Name: COMPUTE_BAYES_SFC_TYPE
!
! Function:
!   Computes the bayesian surface type given the ancillary sfc data
!
!====================================================================
 function COMPUTE_BAYES_SFC_TYPE(Land_Temp, Coast_Temp, Snow_Temp, Sfc_Type_Temp, &
                                 Lat_Temp, Lon_Temp, Sst_Back_Uni_Temp,&
                                 Emiss_Sfc_375um_Temp, symbol) &
                                 result(Bayes_Mask_Sfc_Type_Temp)
   integer(kind=int1), intent(in):: Land_Temp
   integer(kind=int1), intent(in):: Coast_Temp
   integer(kind=int1), intent(in):: Snow_Temp
   integer(kind=int1), intent(in):: Sfc_Type_Temp
   real(kind=real4), intent(in):: Lat_Temp
   real(kind=real4), intent(in):: Lon_Temp
   real(kind=real4), intent(in):: Sst_Back_Uni_Temp
   real(kind=real4), intent(in):: Emiss_Sfc_375um_Temp
   TYPE(symbol_naive_bayesian), intent(in) :: symbol

   integer(kind=int4):: Bayes_Mask_Sfc_Type_Temp


   if (Land_Temp == symbol%LAND) then
           Bayes_Mask_Sfc_Type_Temp = 3
   else
           Bayes_Mask_Sfc_Type_Temp = 1
   endif

   !--- #2 - Shallow Ocean
   if ((Land_Temp == symbol%MODERATE_OCEAN) .or. &
       (Land_Temp == symbol%DEEP_INLAND_WATER) .or. &
       (Land_Temp == symbol%SHALLOW_INLAND_WATER) .or. &
       (Land_Temp == symbol%SHALLOW_OCEAN)) then
           Bayes_Mask_Sfc_Type_Temp = 2
   endif
   if ((Land_Temp /= symbol%LAND) .and. &
       (Sst_Back_Uni_Temp > 0.5)) then
           Bayes_Mask_Sfc_Type_Temp = 2
   endif

   !--- #3 Unfrozen_Land 
   if ((Land_Temp == symbol%LAND) .or. &
       (Land_Temp == symbol%COASTLINE) .or. &
       (Coast_Temp == symbol%YES) .or. &
       (Land_Temp == symbol%EPHEMERAL_WATER)) then

           Bayes_Mask_Sfc_Type_Temp = 3

   endif

   !--- #4 - Snow Covered Land
   if ((Lat_Temp > -60.0) .and. (Snow_Temp == symbol%SNOW)) then
           Bayes_Mask_Sfc_Type_Temp = 4
   endif

   !--- #5 - Arctic
   if ((Lat_Temp >= 0.0) .and. (Snow_Temp == symbol%SEA_ICE)) then
           Bayes_Mask_Sfc_Type_Temp = 5
   endif

   !--- #6 - Antarctic
   if ((Lat_Temp <= -60.0) .and. (Snow_Temp == symbol%SNOW)) then
           Bayes_Mask_Sfc_Type_Temp = 6
   endif
   if ((Lat_Temp <= 0.0) .and. (Snow_Temp == symbol%SEA_ICE)) then
           Bayes_Mask_Sfc_Type_Temp = 6
   endif
   if ((Lat_Temp >= 60.0) .and.  &
       (Lon_Temp > -75.0) .and. (Lon_Temp < -10.0) .and. &
       (Land_Temp == symbol%LAND .or. Land_Temp == symbol%COASTLINE) .and. &
       (Snow_Temp == symbol%SNOW)) then
           Bayes_Mask_Sfc_Type_Temp = 6
   endif

   !--- #7 - Desert
   if ( (Emiss_Sfc_375um_Temp < 0.90 ) .and. (abs(Lat_Temp) < 60.0) .and. &
       ((Sfc_Type_Temp == symbol%OPEN_SHRUBS_SFC) .or. (Sfc_Type_Temp == symbol%BARE_SFC)) ) then
           Bayes_Mask_Sfc_Type_Temp = 7
   endif
 
 
   if ( Bayes_Mask_Sfc_Type_Temp == 3 .and.  &
        Emiss_Sfc_375um_Temp < 0.93  .and.  &
        abs(Lat_Temp) < 60.0 .and. &
       ((Sfc_Type_Temp == symbol%OPEN_SHRUBS_SFC) .or.  &
        (Sfc_Type_Temp == symbol%CLOSED_SHRUBS_SFC) .or. &
        (Sfc_Type_Temp == symbol%GRASSES_SFC) .or.  &
        (Sfc_Type_Temp == symbol%BARE_SFC)) ) then
           Bayes_Mask_Sfc_Type_Temp = 7
   endif

  return 

 end function COMPUTE_BAYES_SFC_TYPE

!-----------------------------------------------------------------------------
! EUMETCAST Fire detection algorithm
!-----------------------------------------------------------------------------
  integer elemental function fire_test( t11, t375, t11_std, t375_std, solzen)

     real, intent(in):: T11
     real, intent(in):: T375
     real, intent(in):: T11_Std
     real, intent(in):: T375_Std
     real, intent(in):: Solzen
     
     real :: Bt_375um_Eumet_Fire_Thresh
     real :: Bt_Diff_Eumet_Fire_Thresh
     real :: Stddev_11um_Eumet_Fire_Thresh 
     real :: Stddev_375um_Eumet_Fire_Thresh

     Fire_Test = 0

     if (T375_std /= Missing_Value_Real4 .and. &
         T11_std /= Missing_Value_Real4) then
         
         !Day
         if (Solzen < EumetCAST_Fire_Day_Solzen_Thresh) then 
            Bt_375um_Eumet_Fire_Thresh = Bt_375um_Eumet_Fire_day_Thresh
            Bt_Diff_Eumet_Fire_Thresh = Bt_Diff_Eumet_Fire_day_Thresh
            Stddev_11um_Eumet_Fire_Thresh = Stddev_11um_Eumet_Fire_Day_Thresh
            Stddev_375um_Eumet_Fire_Thresh = Stddev_375um_Eumet_Fire_Day_Thresh
         endif
         
         !Night
         if (solzen > EumetCAST_Fire_Night_Solzen_Thresh) then 
            Bt_375um_Eumet_Fire_Thresh = Bt_375um_Eumet_Fire_Night_Thresh
            Bt_Diff_Eumet_Fire_Thresh = Bt_Diff_Eumet_Fire_Night_Thresh
            Stddev_11um_Eumet_Fire_Thresh = Stddev_11um_Eumet_Fire_Night_Thresh
            Stddev_375um_Eumet_Fire_Thresh = Stddev_375um_Eumet_Fire_Night_Thresh
         endif
         
         if ((Solzen >= EumetCAST_Fire_Day_Solzen_Thresh) .and. &
             (Solzen <= EumetCAST_Fire_Night_Solzen_Thresh)) then
             
             !linear fit day -> night
             Bt_375um_Eumet_Fire_Thresh = ((-1.0)* Solzen) + 380.0
             Bt_Diff_Eumet_Fire_Thresh = ((-0.4)* Solzen) + 36.0
             
             !These two don't change, but 
             Stddev_11um_Eumet_Fire_Thresh = ((0.0)* solzen) + 1.0
             Stddev_375um_Eumet_Fire_Thresh = ((0.0)* solzen) + 4.0
         
         endif

       ! All of these conditions need to be met
       if ((T375 > bt_375um_Eumet_Fire_Thresh) .and. &
           ((T375 - t11) > Bt_Diff_Eumet_Fire_Thresh) .and. &
           (T375_Std > Stddev_375um_Eumet_Fire_Thresh) .and. &
           (T11_Std < Stddev_11um_Eumet_Fire_Thresh)) then
         Fire_Test = 1
       endif
         
     endif

  end function fire_test

  !-------------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------------
  real elemental function split_window_test( t11_clear, t12_clear, t11, t12)

     real, intent(in):: t11_clear
     real, intent(in):: t12_clear 
     real, intent(in):: t11
     real, intent(in):: t12

     split_window_test  = (t11_clear - t12_clear) * (t11 - 260.0) / (t11_clear - 260.0)

     if (t11_clear <=265.0) split_window_test = 0.0

     split_window_test = (t11 - t12) - split_window_test

  end function split_window_test

  real elemental function reflectance_gross_contrast_test(ref_clear,ref)
     real, intent(in):: ref_clear
     real, intent(in):: ref

     reflectance_gross_contrast_test = Missing_Value_Real4
     if (ref_clear /= Missing_Value_Real4) then
           reflectance_gross_contrast_test = ref - ref_clear
     endif

  end function reflectance_gross_contrast_test

  real elemental function relative_visible_contrast_test(ref_min,ref)
     real, intent(in):: ref_min
     real, intent(in):: ref

     relative_visible_contrast_test = Missing_Value_Real4
     if (ref_min /= Missing_Value_Real4) then
           relative_visible_contrast_test = ref - ref_min
     endif

  end function relative_visible_contrast_test

  real elemental function reflectance_ratio_test(ref_vis,ref_nir)
     real, intent(in):: ref_vis
     real, intent(in):: ref_nir

     reflectance_ratio_test = Missing_Value_Real4
     if (ref_vis /= Missing_Value_Real4 .and. ref_nir /= Missing_Value_Real4) then
           reflectance_ratio_test = ref_nir / ref_vis
     endif

  end function reflectance_ratio_test

  real elemental function nir_reflectance_gross_contrast_test(ref_clear,ref)
     real, intent(in):: ref_clear
     real, intent(in):: ref

     nir_reflectance_gross_contrast_test = Missing_Value_Real4
     if (ref /= Missing_Value_Real4 .and. ref_clear /= Missing_Value_Real4) then
           nir_reflectance_gross_contrast_test = ref- ref_clear
     endif

  end function nir_reflectance_gross_contrast_test

  real elemental function emiss_375um_test(ems,ems_clear)
     real, intent(in):: ems_clear
     real, intent(in):: ems

     emiss_375um_test = Missing_Value_Real4
     if (ems /= Missing_Value_Real4 .and. ems_clear /= Missing_Value_Real4) then
           emiss_375um_test = (ems- ems_clear) / ems_clear
     endif

  end function emiss_375um_test

  real elemental function emiss_375um_day_test(ems,ems_clear)
     real, intent(in):: ems_clear
     real, intent(in):: ems

     emiss_375um_day_test = Missing_Value_Real4
     if (ems /= Missing_Value_Real4 .and. ems_clear /= Missing_Value_Real4) then
           emiss_375um_day_test = (ems- ems_clear) / ems_clear
     endif

  end function emiss_375um_day_test

  real elemental function emiss_375um_night_test(ems,ems_clear)
     real, intent(in):: ems_clear
     real, intent(in):: ems

     emiss_375um_night_test = Missing_Value_Real4
     if (ems /= Missing_Value_Real4 .and. ems_clear /= Missing_Value_Real4) then
           emiss_375um_night_test = ems
         ! emiss_375um_night_test = (ems- ems_clear) / ems_clear
     endif

  end function emiss_375um_night_test
!------------------------------------------------------------------------
! subroutine PACK_BITS_INTO_BYTES(input_bits,bit_depth,output_bytes)
!
! Routines to pack individual bytes into a single byte
!
! input:
! input_bits - vector of bytes to be packed into output_byte
! bit_start - vector of bit starting positions (1-7) for each input byte
! bit_depth - vector of bit depths (1-7) for each input byte (total can not exceed 8)
!
! output: 
! output_byte - byte variable that holds the bit values of the input_bits - can be i1 or i2
!
! local
!  n_in - number of elements of input vectors
!  i_in - index of input_bits (1-n_in)
!  n_out - number of elements of output vectors
!  i_out - index of output_bytes (1-n_out)
!
! Note:
! 1.  if the input byte has information in bits greater then bit depth - they are removed 
!
!
! Example, pack an input_byte wth  bit_start = 2 and bit depth 3 
!
! input byte
!           x x x
! _ _ _ _ _ _ _ _
!
! result of first ishft
! x x x
! _ _ _ _ _ _ _ _
!
! result of second ishft
!       x x x
! _ _ _ _ _ _ _ _

!
! Author: Andrew Heidinger
!
! Version History:  
! February 2006 - Created
!-----------------------------------------------------------------------------------

!--- This Version packs into one byte words
   subroutine PACK_BITS_INTO_BYTES(input_bits,bit_depth,output_bytes)
    integer(kind=int1), dimension(:), intent(in):: input_bits
    integer(kind=int1), dimension(:), intent(in):: bit_depth
    integer(kind=int1), dimension(:), intent(out):: output_bytes
    integer(kind=int1):: bit_start, bit_end, bit_offset
    integer(kind=int1):: temp_byte
    integer:: n_in,i_in,n_out,i_out
    integer, parameter:: word_bit_depth = 8

    !--- determine size of vectors
    n_in = size(input_bits)
    n_out = size(output_bytes)

    !--- reset output byte
    output_bytes = 0

    !--- initialize
    bit_offset = 0
    bit_start = 0
    bit_end = 0
    i_out = 1

    !--- loop through input bytes
    do i_in = 1, n_in

     !--- determine starting and ending bit locations
     bit_start = bit_offset + 1
     bit_end = bit_start + bit_depth(i_in) - 1

     !--- determine if this input byte will fit on current output byte, if not go to next
     if (bit_end > word_bit_depth) then
      i_out = i_out + 1
      bit_offset = 0
      bit_start = bit_offset + 1
      bit_end = bit_start + bit_depth(i_in) - 1
     endif

     !--- check for exceeding the space allowed for the packed bytes
     if (i_out > n_out) then
       print *, "ERROR: Insufficient space for bit packing ", i_out, n_out
       return
     endif

     !--- place input byte into correct position
     temp_byte =0
     temp_byte = ishft(input_bits(i_in),word_bit_depth-bit_depth(i_in))    !first ishft
     temp_byte = ishft(temp_byte,bit_end - word_bit_depth)                 !second ishft

     !--- modify output byte
     output_bytes(i_out) = output_bytes(i_out) + temp_byte

     !--- update bit offset
     bit_offset = bit_offset + bit_depth(i_in)

    enddo

  end subroutine  PACK_BITS_INTO_BYTES

end module
