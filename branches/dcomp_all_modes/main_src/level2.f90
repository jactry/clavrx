! $Id$
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3
!
! NAME: level2.f90 (src)
!       LEVEL2_ROUTINES (program)
!
! PURPOSE: Routines for creating, writing and closing pixel-level output files
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
! REVISION HISTORY:
!          5/09 - Created
!          May 2016 --   store product details in CSV file (AW)
!--------------------------------------------------------------------------------------
module LEVEL2_ROUTINES
   
   use CONSTANTS, only: &
      int4 &
      , int1 &
      , int2 &
      , sym &
      , real4 &
      , MISSING_VALUE_INT1 &
      , MISSING_VALUE_INT2 &
      , MISSING_VALUE_INT4 &
      , MISSING_VALUE_REAL4 &
      , EXE_PROMPT &
      , dcomp_version &
      , acha_version
      
   use PIXEL_COMMON, only: &
      sensor &
      , image &
      , compress_flag &
      , nav &
      , acha &
      , geo &
      , sfc &
      , ch &
      , dcomp_mode &
      , dark_composite_name &
      , bayesian_cloud_mask_name &
      , bad_scan_flag &
      , bad_pixel_mask &
      , utc_scan_time_hours &
      , scan_number &
      , utc_scan_time_hours &
      , scan_number &
      , gap_pixel_mask &
      , diag_pix_array_1 &
      , diag_pix_array_2 &
      , diag_pix_array_3 &
      , num_scans_level2_hdf &
      , dir_level2 &
      , rtm_file_flag &
      , therm_cal_1b &
      , ref_cal_1b &
      , nav_opt &
      , use_sst_anal &
      , modis_clr_alb_flag &
      , nwp_opt &
      , scan_time &
      , utc_scan_time_hours &
      , line_idx_min_segment &
      , ref_ch1_min_3x3 &
      , ref_ch1_std_3x3 &
      , posterior_cld_probability &
      , cld_mask &
      , adj_pix_cld_mask &
      , bayes_mask_sfc_type_global &
      , cld_type &
      , cld_phase &
      , zc_h2o &
      , zc_opaque_cloud &
      , tau_dcomp &
      , reff_dcomp &
      , tau_dcomp_cost &
      , reff_dcomp_cost &
      , tau_dcomp_qf &
      , reff_dcomp_qf &
      , dcomp_quality_flag &
      , dcomp_info_flag &
      , insolation_dcomp &
      , insolation_dcomp_diffuse &
      , tau_nlcomp &
      , reff_nlcomp &
      , tau_nlcomp_cost &
      , reff_nlcomp_cost &
      , nlcomp_quality_flag &
      , nlcomp_info_flag &
      , cloud_063um_albedo &
      , cloud_063um_transmission_solar &
      , cloud_fraction &
      , cloud_fraction_uncer &
      , ndvi_sfc &
      , Tsfc_Retrieved &
      , Tair_Nwp_Pix &
      , Tsfc_Nwp_Pix &
      , Rh_Nwp_Pix &
      , Psfc_Nwp_Pix &
      , Pmsl_Nwp_Pix &
      , K_Index_Nwp_Pix &
      , Cwp_Nwp_Pix &
      , Cfrac_Nwp_Pix &
      , Pc_Nwp_Pix &
      , Ncld_Layers_Nwp_Pix &
      , cld_type_nwp_pix &
      , Ttropo_Nwp_Pix &
      , LCL_Height_Nwp_Pix &
      , CCL_Height_Nwp_Pix &
      , ch1_counts &
      , ch2_counts &
      , ch6_counts &
      , tpw_nwp_pix &
      , Ref_Ch1_Mean_3x3 &
      , Sst_Unmasked &
      , Wnd_Spd_10m_Nwp_Pix &
      , Wnd_Dir_10m_Nwp_Pix &
      , Ref_Ch1_Dark_Composite &
      , cwp_dcomp &
      , rain_rate_dcomp &
      , Wnd_Spd_Cld_Top_Nwp_Pix &
      , wnd_dir_cld_top_nwp_pix &
      , orbital_processing_time_minutes &
      , nonconfident_cloud_mask_fraction &
      , dcomp_success_fraction
      
   use HDF, only: &
      DFACC_CREATE &
      , DFNT_INT8 &
      , DFNT_CHAR8 &
      , DFNT_FLOAT32 &
      , DFNT_INT32 &
      , DFNT_INT16
      
   use SCALING_PARAMETERS, only: &
      one_byte_min &
      , one_byte_max &
      , two_byte_min &
      , two_byte_max
      
   use HDF_PARAMS, only: &
      scale_vector_i2_rank2 &
      , write_clavrx_hdf_global_attributes
      
   use AVHRR_MODULE,only: &
     CH3A_GAIN_HIGH &
     , c1 &
     , c2 &
     , planck_a1 &
     , planck_a2 &
     , planck_nu &
     , solar_ch20_nu &
     , sun_earth_distance &
     , ch1_gain_low &
     , ch1_gain_high &
     , ch1_switch_count_cal &
     , ch1_dark_count_cal &
     , ch2_gain_low &
     , ch2_gain_high &
     , ch2_switch_count_cal &
     , ch2_dark_count_cal &
     , ch3a_gain_low &
     , ch3a_gain_high &
     , ch3a_switch_count_cal &
     , ch3a_dark_count_cal  &
     , cloud_mask_version   &  !- comes from anywhere else!
     , cloud_mask_thresholds_version &   !- comes from anywhere else!
     , acha_version      !- comes from anywhere else!
     
    
   use CLOUD_TYPE_BRIDGE_MODULE,only: &
      cloud_type_version &
      , SET_CLOUD_TYPE_VERSION
      
   use clavrx_message_module, only: &
       mesg
   
   use csv_mod, only: &
      csv_file_close_read &
      , csv_file_line_count &
      , csv_file_open_read &
      , csv_value_count
      
   use strings, only: &
      split

   implicit none
   private

   public:: WRITE_PIXEL_HDF_RECORDS, &
            CLOSE_PIXEL_HDF_FILES

   private::DEFINE_PIXEL_2D_SDS, &
            DEFINE_PIXEL_3D_SDS, &
            DEFINE_HDF_FILE_STRUCTURES, &
            WRITE_ALGORITHM_ATTRIBUTES, &
            file_root_from_l1b

 
   !----------------------------------------------------------------------
   ! the following variables taken from process_avhrr_clavr
   !----------------------------------------------------------------------
   !--- hdf specific variables
   integer(kind=int4),private:: Dim_Id

   !--- compression variables
   integer(kind=int4), dimension(2),private, save::  Comp_Prm
   integer(kind=int4), private, save::  Comp_Type

   integer(kind=int4),parameter,private:: Sds_Rank_1d = 1
   integer(kind=int4),dimension(1),private:: Sds_Dims_1d

   integer(kind=int4),parameter,private:: Sds_Rank_2d = 2
   integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Dims_2d
   integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Start_2d
   integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Stride_2d
   integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Edge_2d
   integer(kind=int4),dimension(Sds_Rank_2d),private:: Sds_Chunk_Size_2d

   integer(kind=int4),parameter,private:: Sds_Rank_3d = 3
   integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Dims_3d
   integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Start_3d
   integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Stride_3d
   integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Edge_3d
   integer(kind=int4),dimension(Sds_Rank_3d),private:: Sds_Chunk_Size_3d

   character(len=11), private, parameter:: MOD_PROMPT = "LEVEL2:"
   character(len=18), private, parameter:: coordinates_string = "longitude latitude"
 
 
    ! csv variables
   integer   ( kind = 4 ) csv_file_status
   integer   ( kind = 4 ) csv_file_unit
   integer   ( kind = 4 ) csv_record_status
   integer   ( kind = 4 ) i
   integer   ( kind = 4 ) line_num
   character ( len = 1000 ) record
   integer   ( kind = 4 ) value_count
   character (len = 300) :: before
   character ( len = 300) :: rec_arr ( 15 )
   integer :: j
   logical :: switch
   integer :: var_dim
   integer(kind = 1) :: scaling
   integer :: dtype
   real :: act_min, act_max
 

   integer :: Num_Level2_Sds 
   integer, private, save:: Sd_Id_Level2
   integer(kind=int4), allocatable, dimension(:):: Sds_Id_Level2
   integer(kind=int4), public, parameter:: NCDC_Attribute_Flag = 0
 
   logical :: file_is_open
 
   character ( len=200) ::csv_file_name

CONTAINS

   !====================================================================
   ! SUBROUTINE Name: DEFINE_HDF_FILE_STRUCTURES
   !
   ! Function:
   !   Determines the structure of the level 2 files.
   !
   ! Description:
   !   This subroutine, given the inputs, determins and opens/creates the apppropriate
   !       Level 2 (pixel level output) for various types of files. In addition
   !       the HDF global attributes are put into the various HDF files 
   !
   !====================================================================
   subroutine DEFINE_HDF_FILE_STRUCTURES(Num_Scans, &
      Dir_Level2, &
      File_1b, &
      Rtm_File_Flag, &
      Level2_File_Flag, &
      c1,c2,a1_20, &
      a2_20, &
      nu_20, &
      a1_31, &
      a2_31, &
      nu_31, &
      a1_32, &
      a2_32, &
      nu_32, &
      Solar_Ch20_nu, &
      Sun_Earth_Distance, &
      Therm_Cal_1b,  &
      Ref_Cal_1b, &
      nav_opt, &
      Use_Sst_anal, &
      Modis_clr_alb_flag, &
      Nwp_Opt, &
      Ch1_gain_low, &
      Ch1_gain_high, &
      Ch1_switch_count, &
      Ch1_dark_count,&
      Ch2_gain_low, &
      Ch2_gain_high, &
      Ch2_switch_count, &
      Ch2_dark_count, &
      Ch3a_Gain_low, &
      Ch3a_Gain_high, &
      Ch3a_Switch_Count, &
      Ch3a_Dark_Count,&
      Start_Year, &
      End_Year, &
      Start_Day, &
      End_Day, &
      Start_Time, &
      End_Time)


      character(len=*), intent(in):: Dir_Level2
      integer(kind=int4), intent(in):: Num_Scans
      character(len=*), intent(in):: File_1b
      integer, intent(in):: Rtm_File_Flag
      integer, intent(in):: Level2_File_Flag
 
      integer(kind=int4), intent(in) :: Modis_Clr_Alb_Flag
      integer(kind=int4), intent(in) :: Nwp_Opt
      integer, intent(in):: therm_cal_1b
      integer, intent(in):: nav_opt
      integer, intent(in):: use_sst_anal
      integer, intent(in):: Ref_cal_1b
      real(kind=real4), intent(in):: c1,c2,a1_20,a2_20,nu_20,a1_31,a2_31,nu_31,a1_32, &
                                a2_32,nu_32,solar_Ch20_nu,sun_earth_distance, &
                                Ch1_gain_low,Ch1_gain_high,Ch1_switch_count, &
                                Ch1_dark_count, Ch2_gain_low,Ch2_gain_high, &
                                Ch2_switch_count,Ch2_dark_count, Ch3a_gain_low,&
                                Ch3a_gain_high,Ch3a_switch_count, &
                                Ch3a_dark_count
      integer(kind=int4), intent(in):: Start_Time
      integer(kind=int4), intent(in):: End_Time
      integer(kind=int2), intent(in):: Start_Year
      integer(kind=int2), intent(in):: End_Year
      integer(kind=int2), intent(in):: Start_Day
      integer(kind=int2), intent(in):: End_Day
 
      character(len=1020):: File_1b_root
      character(len=1020):: File_Rtm
      character(len=1020):: File_Level2

      character(len=1020):: Long_Name_Temp
 
      character(len=128):: Sds_Name

      integer(kind=int4):: blank_int4
      character(len=1):: blank_char
      real:: blank_real4
      real(kind=real4):: Resolution_KM
      integer:: Istatus_Sum
      integer:: Istatus
      integer:: ipos
      integer:: ilen
 
      ! HDF function declarations
      integer:: sfstart
      integer:: sfcreate
      integer:: sfscatt
      integer:: sfsnatt
      integer:: erstat
 
      integer::k
   
   
      !--- begin executable code
      blank_int4 = 0
      blank_real4 = 0.0
      blank_char = " "
  
      File_1b_Root = file_root_from_l1b ( file_1b, Sensor%Sensor_Name,Sensor%Spatial_Resolution_Meters)

      !--- set Resolution_KM for global attribute
      Resolution_KM = Sensor%Spatial_Resolution_Meters / 1000.0

      !-------------------------------------------------------------
      ! define chunking here
      !-------------------------------------------------------------
      !--- 3d , first dim is sds dependent
      Sds_Dims_3d(2) = Image%Number_Of_Elements
      Sds_Dims_3d(3) = Num_Scans
      Sds_Chunk_Size_3d(2) = Image%Number_Of_Elements
      Sds_Chunk_Size_3d(3) = Num_Scans

      !-- dimension of 2d variables
      Sds_Dims_2d(1) = Image%Number_Of_Elements
      Sds_Dims_2d(2) = Image%Number_Of_Lines
      Sds_Chunk_Size_2d(1) = Image%Number_Of_Elements
      Sds_Chunk_Size_2d(2) = Image%Number_Of_Lines_Per_Segment

      !-- dimension of 1d variables
      Sds_Dims_1d(1) = Image%Number_Of_Lines

      !-------------------------------------------------------------
      ! define compression here
      !-------------------------------------------------------------
      Comp_Type = 0                  !no compression
      Comp_Prm(1) = 0
      Comp_Prm(2) = 0

      if (Compress_Flag == 1) then  !gzip compression
         Comp_Type = 4
         Comp_Prm(1) = 6
         Comp_Prm(2) = 0
      endif

      if (Compress_Flag == 2) then  !szip compression
         Comp_Type = 5
         Comp_Prm(1) = 32
         Comp_Prm(2) = 2
      endif


      !---------------------------------------------------------
      !-- define level2 file structure
      !---------------------------------------------------------
      if (Level2_File_Flag == sym%YES) then
         File_Level2= trim(File_1b_Root)//".level2.hdf"
         call mesg (MOD_PROMPT//"creating level-2 file "//trim(file_Level2))
         
         Sd_Id_Level2 = sfstart(trim(Dir_Level2)//trim(file_Level2),DFACC_CREATE)
         if (Sd_Id_Level2 < 0) then
            erstat = 68
            print *, EXE_PROMPT, MOD_PROMPT, "Level-2 file creation failed. Exiting..."
            stop 68
         end if

         !--- write attributes
         call WRITE_CLAVRX_HDF_GLOBAL_ATTRIBUTES(Sd_Id_Level2,"PIXEL",File_Level2,File_1b_Root, &
                           Resolution_KM, &
                           start_year,end_year,start_day,end_day,start_time,end_time,&
                           blank_int4,blank_int4,blank_char,blank_real4, &
                           therm_cal_1b,Ref_cal_1b,nav_opt,use_sst_anal, &
                           modis_clr_alb_flag, nwp_opt, Ch1_gain_low, Ch1_gain_high, &
                           Ch1_switch_count, Ch1_dark_count, &
                           Ch2_gain_low, Ch2_gain_high, &
                           Ch2_switch_count, Ch2_dark_count, &
                           Ch3a_gain_low, Ch3a_gain_high, &
                           Ch3a_switch_count, Ch3a_dark_count, &
                           sun_earth_distance, &
                           c1, c2, a1_20, a2_20, nu_20, &
                           a1_31, a2_31, nu_31, a1_32, a2_32, nu_32, &
                           solar_Ch20_nu,Nav%Timerr_seconds, &
                           acha%mode, dcomp_mode,Sensor%WMO_Id, &
                           Sensor%Platform_Name, Sensor%Sensor_Name, &
                           Dark_Composite_Name,Bayesian_Cloud_Mask_Name)

         !-- reset status flag for error checking
         Istatus_Sum = 0

      
         ! ---  rread in products from csv file
     
         csv_file_name='clavrx_level2_products.csv'
         call csv_file_line_count ( csv_file_name, line_num )
   
         call csv_file_open_read ( csv_file_name, csv_file_unit )
         Num_Level2_Sds = line_num
         allocate( Sds_Id_Level2(line_num))
      
      
         do i = 1, line_num
         
            read ( csv_file_unit, '(a)', iostat = csv_file_status ) record
            call csv_value_count ( record, csv_record_status, value_count )
            rec_arr = extract_single ( trim ( record ) )
        
            switch = trim(rec_arr(1))  .eq. "1"
         
            if ( i == 1) cycle
            read ( rec_arr(2), * ) var_dim
            read ( rec_arr(5), * ) dtype
            read ( rec_arr(6), * ) scaling
            read ( rec_arr(7), * ) act_min
            read ( rec_arr(8), * ) act_max
              
            if ( switch ) then
               select case (var_dim)
               case ( 1 )
                  select case (dtype)
                  case(1)
                  Sds_Id_Level2(i) = sfcreate(Sd_Id_Level2,trim(rec_arr(3)),DFNT_INT8,Sds_Rank_1d,Sds_Dims_1d)
                  Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "SCALED", DFNT_INT8, 1, scaling) + Istatus_Sum
                  Istatus_Sum = sfscatt(Sds_Id_Level2(i), "units", DFNT_CHAR8, 4, trim(rec_arr(14))) + Istatus_Sum
                  Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "_FillValue", DFNT_INT8,  &
                     1, Missing_Value_Int1) + Istatus_Sum
                  Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "RANGE_MISSING", DFNT_FLOAT32,  &
                     1, Real(Missing_Value_Int1,kind=real4)) + Istatus_Sum
               
                  case(3)
                  Sds_Id_Level2(i) = sfcreate(Sd_Id_Level2,"scan_line_number",DFNT_INT32,Sds_Rank_1d,Sds_Dims_1d)
                  Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "SCALED", DFNT_INT8, 1, sym%NO_SCALING) + Istatus_Sum
                  Istatus_Sum = sfscatt(Sds_Id_Level2(i), "units", DFNT_CHAR8, 4, "none") + Istatus_Sum
                  Istatus_Sum = sfscatt(Sds_Id_Level2(i), "standard_name", DFNT_CHAR8, 13, "not specified") + Istatus_Sum
                  Istatus_Sum = sfscatt(Sds_Id_Level2(i), "long_name", DFNT_CHAR8,  &
                     len_trim(Long_Name_temp), trim(Long_Name_temp)) + Istatus_Sum     
                  Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "RANGE_MISSING", DFNT_FLOAT32, &
                    1,real(Missing_Value_Int4,kind=real4)) + Istatus_Sum
                  Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "_FillValue", DFNT_INT32, &
                    1,Missing_Value_Int4) + Istatus_Sum
                    
                  case(4)
                  Sds_Id_Level2(i) = sfcreate(Sd_Id_Level2,trim(rec_arr(3)),DFNT_INT32,Sds_Rank_1d,Sds_Dims_1d)    
                  Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "SCALED", DFNT_INT8, 1, scaling) + Istatus_Sum  
               
                  Istatus_Sum = sfscatt(Sds_Id_Level2(i), "units", DFNT_CHAR8, 4,trim(rec_arr(14) )) + Istatus_Sum 
                  Istatus_Sum = sfscatt(Sds_Id_Level2(i), "standard_name", DFNT_CHAR8, 13, trim(rec_arr(12) )) + Istatus_Sum 
                  Istatus_Sum = sfscatt(Sds_Id_Level2(i), "long_name", DFNT_CHAR8,  &
                     len_trim(trim(rec_arr(14))),  trim(rec_arr(14)) ) + Istatus_Sum
                  Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "RANGE_MISSING", DFNT_FLOAT32, &
                     1,real(Missing_Value_Int4,kind=real4)) + Istatus_Sum           
                  Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "_FillValue", DFNT_INT32, &
                    1,Missing_Value_Int4) + Istatus_Sum  
                  end select  
               case(2)
                  select case (dtype)
                  
                  case(1)
                  call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(i),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               trim(rec_arr(3)), &
                               trim(rec_arr(12)), &
                               trim(rec_arr(15)), &
                               DFNT_INT8, scaling, act_min, act_max, &
                               trim(rec_arr(14) ), Real(Missing_Value_Int1,kind=real4), Istatus)
                  Istatus_Sum = Istatus_Sum + Istatus
                  
                  case(2)
                  call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(i),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              trim(rec_arr(3)), &
                               trim(rec_arr(12)), &
                               trim(rec_arr(15)), &
                               DFNT_INT16, scaling, act_min, act_max, &
                               trim(rec_arr(14) ), Missing_Value_Real4, Istatus)
                  Istatus_Sum = Istatus_Sum + Istatus
              
                  case(4)
                  call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(i),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d, &
                              trim(rec_arr(3)), &
                              trim(rec_arr(12)), &
                              trim(rec_arr(15)), &
                              DFNT_FLOAT32, sym%NO_SCALING, &
                              0.0,0.0, "unknown", Missing_Value_Real4, Istatus)
                  Istatus_Sum = Istatus_Sum + Istatus
                  end select
               
               end select
            end if
         end do
 
         call csv_file_close_read ( csv_file_name, csv_file_unit )


         !--- check for and report errors
         if (Istatus_Sum /= 0) then
            print *, EXE_PROMPT, MOD_PROMPT, "Error defining sds in level2 hdf file"
            stop
         end if
      end if
 
   end subroutine DEFINE_HDF_FILE_STRUCTURES
   
   !
   !
   !
   function extract_single ( record) result  (record_single)
      implicit none
      character ( len = * ) :: record
      character ( len = 300) :: record_single ( 15 )
      character (len = 1000 ) :: record_local
      character ( len =300) :: before
      integer :: i
      integer :: n_val
      
      record_local = trim (record)
      n_val = size ( record_single)
      record_single(:) = 'not_set' 
      do i = 1, n_val 
         call split ( record_local , ',', before)
         record_single ( i ) = trim(before)
      end do
      
   end function extract_single

   !====================================================================
   ! SUBROUTINE Name: WRITE_PIXEL_HDF_RECORDS
   !
   ! Function:
   !   Writes output to appropriate Level 2 files
   !
   ! Description:
   !   This subroutine, given the flags that determine which files are 
   !   being created, outputs the various pixel level outputs to the
   !   appropriate level 2 files and appropriate SDSs for a given segment
   !   this program is called in process_clavrx
   !
   !====================================================================
   subroutine WRITE_PIXEL_HDF_RECORDS(segment_number,Level2_File_Flag)
      integer, intent(in) :: segment_number
      integer, intent(in) :: Level2_File_Flag
   
      integer:: Istatus
      integer:: Line_Idx

      ! HDF function declarations
      integer:: sfwdata
      integer (kind=int1), allocatable :: data_dim1_dtype1(:)
      integer (kind=int4), allocatable :: data_dim1_dtype3(:)
      real(kind=real4), allocatable ::data_dim1_dtype4(:)
      integer (kind=int1), allocatable :: data_dim2_dtype1(:,:)
      real , allocatable :: data_dim2_dtype2(:,:)
      integer (kind=int4), allocatable :: data_dim2_dtype3(:,:)
      real(kind=real4), allocatable ::data_dim2_dtype4(:,:)
      character (len=40) :: name
      integer(kind=int2), dimension(:,:),allocatable :: Two_Byte_dummy
    
      if (Segment_Number == 1) then

         !--- place algorithm cvs tags into global strings for output
         call SET_CLOUD_TYPE_VERSION()

         Num_Scans_Level2_Hdf = 0

         call DEFINE_HDF_FILE_STRUCTURES(Image%Number_Of_Lines, &
                              Dir_Level2, &
                              Image%Level1b_Name, &
                              Rtm_File_Flag, &
                              Level2_File_Flag, &
                              c1,c2,planck_a1(20),planck_a2(20),planck_nu(20), &
                              planck_a1(31),planck_a2(31),planck_nu(31), &
                              planck_a1(32),planck_a2(32),planck_nu(32),Solar_Ch20_Nu,&
                              Sun_Earth_Distance,Therm_Cal_1b, &
                              Ref_Cal_1b,Nav_Opt,Use_Sst_Anal, &
                              Modis_Clr_Alb_Flag,Nwp_Opt, &
                              Ch1_Gain_Low,Ch1_gain_High, &
                              Ch1_Switch_Count_Cal,Ch1_Dark_Count_Cal, &
                              Ch2_Gain_low,Ch2_Gain_High, &
                              Ch2_Switch_Count_Cal,Ch2_Dark_Count_Cal, &
                              Ch3a_Gain_low,Ch3a_Gain_High, &
                              Ch3a_Switch_Count_Cal,Ch3a_Dark_Count_Cal, &
                              Image%Start_Year,Image%End_Year,Image%Start_Doy,Image%End_Doy,&
                              Image%Start_Time,Image%End_Time)
            end if
            
    
      
   !-----------------------------------------------------------------------
   ! Get time of each scan line and convert to scale
   !-----------------------------------------------------------------------
   where(Scan_Time == Missing_Value_Real4)
            Utc_Scan_Time_Hours = Missing_Value_Real4
   else where
            Utc_Scan_Time_Hours = Scan_Time/60.0/60.0/1000.0
   end where

   !--------------------------------------------------------------------------------
   ! determine start and edges for writing this segment's data to pixel hdf file
   !-----------------------------------------------------------------------------
   Sds_Start_2d(1) = 0     !pixel dimension
   Sds_Start_2d(2) = Num_Scans_Level2_Hdf

   Sds_Stride_2d(1) = 1
   Sds_Stride_2d(2) = 1

   Sds_Edge_2d(1) = Image%Number_Of_Elements
   Sds_Edge_2d(2) = min(Image%Number_Of_Lines_Read_This_Segment,Image%Number_Of_Lines - Sds_Start_2d(2))

   if (Sds_Edge_2d(2) <= 0) then
      return
   end if

   !--- update Num_Scans_Level2_Hdf
   Num_Scans_Level2_Hdf = min(Image%Number_Of_Lines,Num_Scans_Level2_Hdf +  &
                               Image%Number_Of_Lines_Read_This_Segment)
   
 
   
   
   !-------------------------------------------------------------------------
   ! write to level2 file
   !-------------------------------------------------------------------------
   if (Level2_File_Flag == sym%YES) then
         istatus = 0
            ! ---  re-read in products from csv file
      
      csv_file_name='clavrx_level2_products.csv'
      call csv_file_line_count ( csv_file_name, line_num )   
      call csv_file_open_read ( csv_file_name, csv_file_unit )
      
              allocate ( data_dim1_dtype1(sds_edge_2d(2)))
   
   allocate ( data_dim1_dtype3(sds_edge_2d(2)))
   allocate ( data_dim1_dtype4(sds_edge_2d(2)))
   allocate ( data_dim2_dtype1(sds_edge_2d(1),sds_edge_2d(2))) 
   allocate ( data_dim2_dtype2(sds_edge_2d(1),sds_edge_2d(2)))
   allocate ( data_dim2_dtype3(sds_edge_2d(1),sds_edge_2d(2)))
   allocate ( data_dim2_dtype4(sds_edge_2d(1),sds_edge_2d(2)))
   allocate ( two_byte_dummy(sds_edge_2d(1),sds_edge_2d(2)))   
         
      do i = 1, line_num
        
         read ( csv_file_unit, '(a)', iostat = csv_file_status ) record
         call csv_value_count ( record, csv_record_status, value_count )
         rec_arr = extract_single ( trim ( record ) )
        
         switch = trim(rec_arr(1))  .eq. "1"
          
         if ( i == 1) cycle
         read( rec_arr(2),  * ) var_dim
         name = trim( rec_arr(3) )
         read ( rec_arr(5), * ) dtype
         read ( rec_arr(6), * ) scaling
         read (rec_arr(7), * ) act_min
         read (rec_arr(8), * ) act_max
         
         if ( switch ) then
            include 'level2_assign.inc'
          
            select case (var_dim)
            case ( 1 )
               select case (dtype)
               case(1)
               Istatus = sfwdata(Sds_Id_Level2(i), Sds_Start_2d(2), Sds_Stride_2d(2),          &
                         Sds_Edge_2d(2), data_dim1_dtype1 ) + Istatus
               case(3)
                 
                 Istatus = sfwdata(Sds_Id_Level2(i), Sds_Start_2d(2), Sds_Stride_2d(2),          &
                         Sds_Edge_2d(2), data_dim1_dtype3 ) + Istatus
               case(4)
                 Istatus = sfwdata(Sds_Id_Level2(i), Sds_Start_2d(2), Sds_Stride_2d(2),          &
                         Sds_Edge_2d(2), data_dim1_dtype4 ) + Istatus
               
               end select  
            case(2)
               select case (dtype)
               case(1)
               
                  Istatus = sfwdata(Sds_Id_Level2(i), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                               data_dim2_dtype1) + Istatus
               case(2)
              
               call SCALE_VECTOR_I2_RANK2(data_dim2_dtype2,sym%LINEAR_SCALING,act_min,act_max,Missing_Value_Real4 &
                  ,Two_Byte_dummy)
                  Istatus = sfwdata(Sds_Id_Level2(i),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d, &
                  Two_Byte_Dummy ) + Istatus
                
               case(4)
                   Istatus = sfwdata(Sds_Id_Level2(i), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                       data_dim2_dtype4 ) + Istatus   
               end select   
            end select
         end if

      end do
      
         if ( allocated ( data_dim1_dtype1)) deallocate ( data_dim1_dtype1)
         if ( allocated ( data_dim1_dtype3)) deallocate ( data_dim1_dtype3)
         if ( allocated ( data_dim1_dtype4)) deallocate ( data_dim1_dtype4)
         if ( allocated ( data_dim2_dtype1)) deallocate ( data_dim2_dtype1)
         if ( allocated ( data_dim2_dtype2)) deallocate ( data_dim2_dtype2)
         if ( allocated ( data_dim2_dtype3)) deallocate ( data_dim2_dtype3)
         if ( allocated ( data_dim2_dtype4)) deallocate ( data_dim2_dtype4)
         if ( allocated ( two_byte_dummy)) deallocate (two_byte_dummy)

      call csv_file_close_read ( csv_file_name, csv_file_unit )
    !--- check for and report errors
    if (Istatus /= 0) then
       print *, EXE_PROMPT, MOD_PROMPT, "Error writing to level2 file: ", Istatus
       stop
    endif
    
   endif

   end subroutine WRITE_PIXEL_HDF_RECORDS


   !====================================================================
   ! SUBROUTINE Name: CLOSE_PIXEL_HDF_FILES
   !
   ! Function:
   !   Closes appropriate Level 2 files
   !
   ! Description:
   !   This subroutine, given the flags that determine which files have 
   !   been created, closes the open output files.
   !
   !====================================================================
   subroutine CLOSE_PIXEL_HDF_FILES(Level2_File_Flag)

 
      integer, intent(in):: Level2_File_Flag

      integer:: Isds
      integer:: Istatus

      ! HDF function declarations
      integer:: sfsnatt
      integer:: sfendacc
      integer:: sfend

      !--- write algorithm attributes to level2
      call WRITE_ALGORITHM_ATTRIBUTES()
      !------------------------------------------------------------------------
      !--- close level2 file
      !------------------------------------------------------------------------
      Istatus = 0
      Istatus = sfsnatt(Sd_Id_Level2, "NUMBER_OF_ELEMENTS", DFNT_INT32,1,Image%Number_Of_Elements)+Istatus
      Istatus = sfsnatt(Sd_Id_Level2, "NUMBER_OF_SCANS_LEVEL1B", DFNT_INT32,1,Image%Number_Of_Lines)+Istatus
      Istatus = sfsnatt(Sd_Id_Level2, "NUMBER_OF_SCANS_LEVEL2", DFNT_INT32,1,Num_Scans_Level2_Hdf)+Istatus
      Istatus = sfsnatt(Sd_Id_Level2, "PROCESSING_TIME_MINUTES", DFNT_FLOAT32,1,Orbital_Processing_Time_Minutes)+Istatus
      Istatus = sfsnatt(Sd_Id_Level2, "NONCONFIDENT_CLOUD_MASK_FRACTION", DFNT_FLOAT32,1,NONCONFIDENT_CLOUD_MASK_Fraction)+Istatus
      Istatus = sfsnatt(Sd_Id_Level2, "ACHA_SUCCESS_FRACTION", DFNT_FLOAT32,1,ACHA%Success_Fraction)+Istatus
      Istatus = sfsnatt(Sd_Id_Level2, "DCOMP_SUCCESS_FRACTION", DFNT_FLOAT32,1,DCOMP_Success_Fraction)+Istatus
    
      if (Level2_File_Flag == sym%YES) then
         do Isds = 1, Num_Level2_Sds
        
            if (sds_Id_Level2(Isds) /= 0 ) then
               Istatus = sfendacc(Sds_Id_Level2(Isds)) + Istatus
            end if
            sds_Id_Level2(Isds) = 0
     
         end do
      
         Istatus = sfend(Sd_Id_Level2) + Istatus

      end if

   end subroutine CLOSE_PIXEL_HDF_FILES

!====================================================================
! SUBROUTINE Name: DEFINE_PIXEL_2D_SDS
!
! Function:
!   Defines a 2D SDS for the level 2 files
!
! Description:
!   This subroutine, given the inputs, creates the SDSs inside the various
!   files. 
!
!====================================================================
  subroutine DEFINE_PIXEL_2D_SDS(Sds_Id,     &
                                 Sd_Id,      &
                                 Sds_Dims,   &
                                 Sds_Chunk,  &
                                 Sds_Name,   &
                                 Sds_Standard_Name,   &
                                 Sds_Long_Name, &
                                 Sds_Type,   &
                                 scaled,     &
                                 Sds_Min,    &
                                 Sds_Max,    &
                                 Sds_Units,  &
                                 Sds_Missing,&
                                 Istatus)

    integer, intent(out):: Sds_Id
    integer, intent(in):: Sd_Id
    integer, dimension(2), intent(in):: Sds_Dims
    integer, dimension(2), intent(in):: Sds_Chunk
    integer, intent(in):: Sds_Type
    real, intent(in):: Sds_Min
    real, intent(in):: Sds_Max
    real, intent(in):: Sds_Missing
    integer(kind=int1), intent(in):: Scaled
    character(len=*), intent(in):: Sds_Name
    character(len=*), intent(in):: Sds_Standard_Name
    character(len=*), intent(in):: Sds_Long_Name
    character(len=*), intent(in):: Sds_Units
    integer, intent(out):: Istatus
    integer:: Dim_Id
    integer(kind=int4):: Scaled_Min
    integer(kind=int4):: Scaled_Max
    integer(kind=int4):: Scaled_Missing
    real(kind=real4):: Add_Offset
    real(kind=real4):: Scale_Factor

    integer:: sfcreate
    integer:: sfsnatt
    integer:: sfscatt
    integer:: sfdimid
    integer:: sfsdmname
    integer:: sfschnk
    
    real :: temp_vector_2  (2)
    integer (kind = int1)  :: temp_vector_2_int1  (2)  
    integer (kind = int2)  :: temp_vector_2_int2  (2) 
    
    integer (kind = int1), parameter  :: temp_indgen_2_int1  (2)  = [0,1]
    integer (kind = int1)  ,parameter :: temp_indgen_4_int1(4) = [0,1,2,3]
    integer (kind = int1)  , parameter :: temp_indgen_8_int1  (8) = [0,1,2,3,4,5,6,7]
    integer (kind = int1) ,parameter  :: temp_indgen_13_int1 (13) =  [0,1,2,3,4,5,6,7,8,9,10,11,12]
    integer (kind = int1)  ,parameter :: temp_indgen_14_int1 (14) =  [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    
    
    
    
    Istatus = 0 

    Sds_Id = sfcreate(Sd_Id,Sds_Name,Sds_Type,Sds_Rank_2d,Sds_Dims)

    !--- Attributes that go into all SDSs
    Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, scaled) + Istatus

    Istatus = sfscatt(Sds_Id, "units", DFNT_CHAR8, len_trim(Sds_Units), trim(Sds_Units)) + Istatus

    Istatus = sfscatt(Sds_Id, "standard_name", DFNT_CHAR8, len_trim(Sds_Standard_Name),  &
                               trim(Sds_Standard_Name)) + Istatus

    Istatus = sfscatt(Sds_Id, "long_name", DFNT_CHAR8, len_trim(Sds_Long_Name),  &
                                trim(Sds_Long_Name)) + Istatus

    if (Sds_Name /= "latitude" .and. Sds_Name /= "longitude") then 
       Istatus = sfscatt(Sds_Id, "coordinates", DFNT_CHAR8, len_trim(coordinates_string),  &
                          trim(coordinates_string)) + Istatus
    endif

    Dim_Id = sfdimid(Sds_Id, 0)
    Istatus = sfsdmname(Dim_Id,"pixel_elements_along_scan_direction") + Istatus

    Dim_Id = sfdimid(Sds_Id, 1)
    Istatus = sfsdmname(Dim_Id,"scan_lines_along_track_direction") + Istatus

    Istatus = sfschnk(Sds_Id,Sds_Chunk,Comp_Type,Comp_Prm)+Istatus

    !--- determine scaled ranges based on Sds_Type
    if (Sds_Type == DFNT_INT8) then
          Scaled_Min = One_Byte_Min
          Scaled_Max = One_Byte_Max
          Scaled_Missing = Missing_Value_Int1
    elseif (Sds_Type == DFNT_INT16) then
          Scaled_Min = Two_Byte_Min
          Scaled_Max = Two_Byte_Max
          Scaled_Missing = Missing_Value_Int2
    endif

    !--- write actual missng and actual rangel for all scaled variables
    if (Scaled /= sym%NO_SCALING) then
!     Istatus = sfsnatt(Sds_Id, "actual_missing", DFNT_FLOAT32, 1, Sds_Missing) + Istatus
      temp_vector_2 = [Sds_Min,Sds_Max]
     Istatus = sfsnatt(Sds_Id, "actual_range", DFNT_FLOAT32, 2, temp_vector_2) + Istatus
    endif

    !--- write valid_range for all scaled variables
    if (Scaled /= sym%NO_SCALING) then

     if (Sds_Type == DFNT_INT8) then
       temp_vector_2_int1 = int((/Scaled_Min,Scaled_Max/),kind=int1)
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2,  temp_vector_2_int1 ) + Istatus
     elseif (Sds_Type == DFNT_INT16) then
       temp_vector_2_int2 = int((/Scaled_Min,Scaled_Max/),kind=int2)
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2, temp_vector_2_int2) + Istatus
     elseif (Sds_Type == DFNT_FLOAT32) then
          temp_vector_2 = [Sds_Min,Sds_Max]
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2, temp_vector_2) + Istatus
     endif

    endif

    !--- write _FillValue for all except bitmasks (missing = -888)
    if (Sds_Missing /= -888.0) then

     !--- write Fill_Value
     if (Sds_Type == DFNT_INT8) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, int(Scaled_Missing,kind=int1)) + Istatus
     elseif (Sds_Type == DFNT_INT16) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, int(Scaled_Missing,kind=int2)) + Istatus
     elseif (Sds_Type == DFNT_FLOAT32) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, Sds_Missing) + Istatus
     endif

    endif

    !--- testing CF flag_meanings and flag_values attributes
    if (Sds_Name == "acha_quality") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 2, temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 143, "Processed "// &
                               " valid_Tc_retrieval "// &
                               " valid_ec_retrieval "// &
                               " valid_beta_retrieval "// &
                               " degraded_Tc_retrieval "// &
                               " degraded_ec_retrieval "// &
                               " degraded_beta_retrieval ") + Istatus
    end if

    if (Sds_Name == "acha_info") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 2, temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 196, "Cloud_Height_Attempted "// &
                               " Bias_Correction_Employed "// &
                               " Ice_Cloud_Retrieval "// &
                               " Local_Radiative_Center_Processing_Used "// &
                               " Multi_Layer_Retrieval "// &
                               " Lower_Cloud_Interpolation_Used "// &
                               " Boundary_Layer_Inversion_Assumed ") + Istatus
    end if

    if (Sds_Name == "dcomp_quality") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 2, temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 120, "Processed "// &
                               " valid_COD_retrieval "// &
                               " valid_REF_retrieval "// &
                               " degraded_COD_retrieval "// &
                               " degraded_REF_retrieval "// &
                               " convergency "// &
                               " glint ") + Istatus
    end if

    if (Sds_Name == "dcomp_info") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 2, temp_indgen_2_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 117, "info_flag_set "// &
                               " land_or_sea_mask "// &
                               " day_or_night mask "// &
                               " twilight_(65-82_solar_zenith) "// &
                               " snow "// &
                               " sea_ice "// &
                               " phase "// &
                               " thick_cloud ") + Istatus
    end if

    if (Sds_Name == "cld_opd_dcomp_qf") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cld_reff_dcomp_qf") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cld_temp_acha_qf") then
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4, temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cld_emiss_acha_qf") then
    
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cld_beta_acha_qf") then
      
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
      Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 49, "not_attempted "// &
                               " failed "// &
                               " low_quality "// &
                               " high_quality ") + Istatus
    end if

    if (Sds_Name == "cloud_mask") then
    
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 4,  temp_indgen_4_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 47, "clear "// &
                               " probably_clear "// &
                               " probably_cloudy "// &
                               " cloudy ") + Istatus
    end if

    if (Sds_Name == "cloud_type") then
      
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 13, temp_indgen_13_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 120, "clear "// &
                               " probably_clear "// &
                               " fog "// &
                               " water "// &
                               " supercooled_water "// &
                               " mixed "// &
                               " opaque_ice "// &
                               " cirrus "// &
                               " overlapping "// &
                               " overshooting "// &
                               " unknown "// &
                               " dust "// &
                               " smoke ") + Istatus
    end if

    if (Sds_Name == "surface_type") then
      
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 14, temp_indgen_14_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 106, "water "// &
                               " evergreen_needle "// &
                               " evergreen_broad "// &
                               " deciduous_needle "// &
                               " deciduous_broad "// &
                               " mixed_forest "// &
                               " woodlands "// &
                               " wooded_grass "// &
                               " closed_shrubs "// &
                               " open_shrubs "// &
                               " grasses "// &
                               " croplands "// &
                               " bare "// &
                               " urban ") + Istatus
    end if

    if (Sds_Name == "land_class") then
      
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 8, temp_indgen_8_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 109, "ocean "// &
                               " land "// &
                               " coastline "// &
                               " shallow_inland_water "// &
                               " ephemeral_water "// &
                               " deep_inland_water "// &
                               " moderate_ocean "// &
                               " deep_ocean ") + Istatus
    end if

    if (Sds_Name == "snow_class") then
     
      Istatus = sfsnatt(Sds_Id, "flag_values", DFNT_INT8, 3,  temp_indgen_4_int1) + Istatus
            Istatus = sfscatt(Sds_Id, "flag_meanings", DFNT_CHAR8, 30, "no_snow_or_ice "// &
                               " sea_ice "// &
                               " snow ") + Istatus
    end if


    !--- Scaling Attributes that go into scaled SDSs
    if (Scaled > 0) then

     !--- determine offset and scale factor
     if (Sds_Min /= Sds_Max) then
       Scale_Factor = (Sds_Max - Sds_Min) / (Scaled_Max - Scaled_Min)
       Add_Offset =  Sds_Min - Scale_Factor*Scaled_Min
     else
       Scale_Factor = 1.0
       Add_Offset = 0.0
     endif

     !--- assume any one-byte integer is not to be scaled
     if (Sds_Missing == -128.0) then
        Scale_Factor = 1.0
        Add_Offset = 0.0
     endif

     !--- write remaining attributes
     if (Scaled == sym%LINEAR_SCALING) then
      Istatus = sfsnatt(Sds_Id, "add_offset", DFNT_FLOAT32, 1, Add_Offset) + Istatus
      Istatus = sfsnatt(Sds_Id, "scale_factor", DFNT_FLOAT32, 1, Scale_Factor) + Istatus
     endif 

     !--- write deprecated scaling attributes that supported non-linear scaling
     Istatus = sfsnatt(Sds_Id, "RANGE_MIN", DFNT_FLOAT32, 1, Sds_Min) + Istatus
     Istatus = sfsnatt(Sds_Id, "RANGE_MAX", DFNT_FLOAT32, 1, Sds_Max) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MISSING", DFNT_INT32, 1, Scaled_Missing) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MIN", DFNT_INT32, 1, Scaled_Min) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MAX", DFNT_INT32, 1, Scaled_Max) + Istatus


     if (Istatus /= 0) print *, "Error writing level2 2d sds named ", trim(Sds_Name)

   endif 

  end subroutine DEFINE_PIXEL_2D_SDS

!====================================================================
! SUBROUTINE Name: DEFINE_PIXEL_3D_SDS
!
! Function:
!   Defines a 3D SDS for the level 2 files
!
! Description:
!   This subroutine, given the inputs, creates the SDSs inside the various
!   files. 
!
!====================================================================
  subroutine DEFINE_PIXEL_3D_SDS(Sds_Id,     &
                                 Sd_Id,      &
                                 Sds_Dims,   &
                                 Sds_Chunk,  &
                                 Sds_Name,   &
                                 Sds_Standard_Name,   &
                                 Sds_Long_Name, &
                                 Sds_Type,   &
                                 scaled,     &
                                 Sds_Min,    &
                                 Sds_Max,    &
                                 Sds_Units,  &
                                 Sds_Missing,&
                                 Dim1_Name, &
                                 Dim2_Name, &
                                 Dim3_Name, &
                                 Istatus)

    integer, intent(out):: Sds_Id
    integer, intent(in):: Sd_Id
    integer, dimension(3), intent(in):: Sds_Dims
    integer, dimension(3), intent(in):: Sds_Chunk
    integer, intent(in):: Sds_Type
    real, intent(in):: Sds_Min
    real, intent(in):: Sds_Max
    real, intent(in):: Sds_Missing
    integer(kind=int1), intent(in):: scaled
    character(len=*), intent(in):: Sds_Name
    character(len=*), intent(in):: Sds_Standard_Name
    character(len=*), intent(in):: Sds_Long_Name
    character(len=*), intent(in):: Sds_Units
    character(len=*), intent(in):: Dim1_Name
    character(len=*), intent(in):: Dim2_Name
    character(len=*), intent(in):: Dim3_Name
    integer, intent(out):: Istatus
    integer:: Dim_Id
    integer(kind=int4):: Scaled_Min
    integer(kind=int4):: Scaled_Max
    integer(kind=int4):: Scaled_Missing
    real(kind=real4):: Add_Offset
    real(kind=real4):: Scale_Factor

    integer:: sfcreate
    integer:: sfsnatt
    integer:: sfscatt
    integer:: sfdimid
    integer:: sfsdmname
    integer:: sfschnk
    
    real :: temp_vector_2  (2)
    integer (kind = int1)  :: temp_vector_2_int1  (2)  
    integer (kind = int2)  :: temp_vector_2_int2  (2) 

    Istatus = 0 

    Sds_Id = sfcreate(Sd_Id,Sds_Name,Sds_Type,3,Sds_Dims)

    !--- Attributes that go into all SDSs
    Istatus = sfsnatt(Sds_Id, "SCALED", DFNT_INT8, 1, Scaled) + Istatus
    Istatus = sfscatt(Sds_Id, "units", DFNT_CHAR8, len_trim(Sds_Units), trim(Sds_Units)) + Istatus

    Istatus = sfscatt(Sds_Id, "standard_name", DFNT_CHAR8, len_trim(Sds_Standard_Name),  &
                               trim(Sds_Standard_Name)) + Istatus

    Istatus = sfscatt(Sds_Id, "long_name", DFNT_CHAR8, len_trim(Sds_Long_Name),  &
                                trim(Sds_Long_Name)) + Istatus
    Dim_Id = sfdimid(Sds_Id, 0)
    Istatus = sfsdmname(Dim_Id,trim(Dim1_Name)) + Istatus

    Dim_Id = sfdimid(Sds_Id, 1)
    Istatus = sfsdmname(Dim_Id,trim(Dim2_Name)) + Istatus

    Dim_Id = sfdimid(Sds_Id, 2)
    Istatus = sfsdmname(Dim_Id,trim(Dim3_Name)) + Istatus

    Istatus = sfschnk(Sds_Id,Sds_Chunk,Comp_Type,Comp_Prm)+Istatus

    !-- Range Missing written out regardless of Scaled Value
    Istatus = sfsnatt(Sds_Id, "RANGE_MISSING", DFNT_FLOAT32, 1, Sds_Missing) + Istatus
    temp_vector_2 = [Sds_Min,Sds_Max]
    Istatus = sfsnatt(Sds_Id, "actual_range", DFNT_FLOAT32, 2, temp_vector_2) + Istatus

    !--- determine scaled ranges based on Sds_Type
    if (Sds_Type == DFNT_INT8) then
          Scaled_Min = One_Byte_Min
          Scaled_Max = One_Byte_Max
          Scaled_Missing = Missing_Value_Int1
    elseif (Sds_Type == DFNT_INT16) then
          Scaled_Min = Two_Byte_Min
          Scaled_Max = Two_Byte_Max
          Scaled_Missing = Missing_Value_Int2
    endif
    
    !--- write valid_range and _FillValue for all except bitmasks (missing = -888)
    if (Sds_Missing /= -888.0) then
   
    !--- write valid range
    if (Sds_Type == DFNT_INT8) then
       temp_vector_2_int1 = int((/Scaled_Min,Scaled_Max/),kind=int1)
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2,  temp_vector_2_int1 ) + Istatus
     elseif (Sds_Type == DFNT_INT16) then
       temp_vector_2_int2 = int((/Scaled_Min,Scaled_Max/),kind=int2)
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2, temp_vector_2_int2) + Istatus
     elseif (Sds_Type == DFNT_FLOAT32) then
          temp_vector_2 = [Sds_Min,Sds_Max]
       Istatus = sfsnatt(Sds_Id, "valid_range", Sds_Type, 2, temp_vector_2) + Istatus
     endif

     !--- write fill_value
     if (Sds_Type == DFNT_INT8) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, int(Scaled_Missing,kind=int1)) + Istatus
     elseif (Sds_Type == DFNT_INT16) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, int(Scaled_Missing,kind=int2)) + Istatus
     elseif (Sds_Type == DFNT_FLOAT32) then
       Istatus = sfsnatt(Sds_Id, "_FillValue", Sds_Type, 1, Sds_Missing) + Istatus
     endif

    endif

    !--- Attributes that go into scaled SDSs
    if (Scaled > 0) then

     !--- determine offset and scale factor
     if (Sds_Min /= Sds_Max) then
       Scale_Factor = (Sds_Max - Sds_Min) / (Scaled_Max - Scaled_Min)
       Add_Offset =  Sds_Min - Scale_Factor*Scaled_Min
     else
       Scale_Factor = 1.0
       Add_Offset = 0.0
     endif
   
     !--- assume any one-byte integer is not to be scaled
     if (Sds_Missing == -128.0) then
        Scale_Factor = 1.0
        Add_Offset = 0.0
     endif

     !--- write remaining attributes
     Istatus = sfsnatt(Sds_Id, "RANGE_MIN", DFNT_FLOAT32, 1, Sds_Min) + Istatus
     Istatus = sfsnatt(Sds_Id, "RANGE_MAX", DFNT_FLOAT32, 1, Sds_Max) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MISSING", DFNT_INT32, 1, Scaled_Missing) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MIN", DFNT_INT32, 1, Scaled_Min) + Istatus
     Istatus = sfsnatt(Sds_Id, "SCALED_MAX", DFNT_INT32, 1, Scaled_Max) + Istatus
     if (Scaled == sym%LINEAR_SCALING) then 
      Istatus = sfsnatt(Sds_Id, "add_offset", DFNT_FLOAT32, 1, Add_Offset) + Istatus
      Istatus = sfsnatt(Sds_Id, "scale_factor", DFNT_FLOAT32, 1, Scale_Factor) + Istatus
     endif

   endif

   if (Istatus /= 0) print *, "Error writing level2 3d sds named ", trim(Sds_Name)

  end subroutine DEFINE_PIXEL_3D_SDS

  !============================================================================
  !
  !============================================================================
  subroutine WRITE_ALGORITHM_ATTRIBUTES()

    integer:: sfscatt
    integer:: istatus
    
    istatus = 0
    istatus = sfscatt(Sd_Id_Level2, "CLOUD_MASK_VERSION", DFNT_CHAR8, len_trim(Cloud_Mask_Version), trim(Cloud_Mask_Version))+istatus
    istatus = sfscatt(Sd_Id_Level2, "CLOUD_MASK_THRESHOLDS_VERSION", DFNT_CHAR8,  &
                 len_trim(Cloud_Mask_Thresholds_Version), trim(Cloud_Mask_Thresholds_Version))+istatus 
    istatus = sfscatt(Sd_Id_Level2, "CLOUD_TYPE_VERSION", DFNT_CHAR8, len_trim(Cloud_Type_Version), trim(Cloud_Type_Version))+istatus 
    istatus = sfscatt(Sd_Id_Level2, "ACHA_VERSION", DFNT_CHAR8, len_trim(ACHA_Version), trim(ACHA_Version))+istatus 
    istatus = sfscatt(Sd_Id_Level2, "DCOMP_VERSION", DFNT_CHAR8, len_trim(dcomp_version), trim(dcomp_version))+istatus

  end subroutine WRITE_ALGORITHM_ATTRIBUTES

!============================================================================

   function file_root_from_l1b (filename, sensor,spatial_resolution_meters) result (file_root)
      character (len = *), intent(in) :: filename
      character (len = *), intent(in) :: sensor
      integer, intent(in)::spatial_resolution_meters
      
      character (len = 200) :: file_root
      
      character(len=4):: l1b_ext
      integer:: ipos,ilen
      
      !--- make a file name base for output
      File_Root = trim (filename)
  

      !--- special processing for modis - remove hdf suffix
      l1b_ext = file_root(len_trim(file_root)-3:len_trim(file_root))
      if (trim(l1b_ext) == ".hdf") then
         file_root = file_root(1:len_trim(file_root)-4)
      endif

      !--- special processing for viirs - remove hdf suffix - this hard coded for
      if (trim(Sensor) == 'VIIRS') then
         file_root = file_root(7:len_trim(file_root)-34)
      endif

      !--- special processing for ahi - remove B01.nc suffix - this hard coded for
      if (trim(Sensor) == 'AHI') then
         file_root = file_root(4:len_trim(file_root)-12)
      endif

      !--- special processing for IFF - remove hdf suffix - this hard coded for
      ! !PEATE files
      ! if (trim(Sensor%Sensor_Name) == 'AQUA-IFF' .or. trim(Sensor%Sensor_Name) == 'AQUA-IFF') then
      !   file_root = file_root(1:len_trim(file_root)-29)
      ! elseif (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then
      !   file_root = file_root(1:len_trim(file_root)-20)
      ! endif

      !--- do this for GOES names which are named goesxx_1_year_jday_hour.area
      if (trim(Sensor) == 'GOES-IL-IMAGER' .or.  &
         trim(Sensor) == 'GOES-MP-IMAGER' .or.  &
         trim(Sensor) == 'COMS-IMAGER' .or.  &
         trim(Sensor) == 'MTSAT-IMAGER' .or.  &
         trim(Sensor) == 'FY2-IMAGER' .or.  &
         trim(Sensor) == 'SEVIRI') then

         !-- remove area suffix
         l1b_ext = file_root(len_trim(file_root)-3:len_trim(file_root))
         if (trim(l1b_ext) == "area") then
            file_root = file_root(1:len_trim(file_root)-5)
         endif
         !-- remove channel number
         if (trim(l1b_ext) == "area") then
            ipos = index(file_root, "_1_")
            ilen = len(file_root)
            file_root = file_root(1:ipos-1) // "_"//file_root(ipos+3:ilen)
         endif
      endif

      !--- add 1km tag for 1km GOES files
      if (index(Sensor,'GOES') > 0 .and. Spatial_Resolution_Meters == 1000) then
         file_root = trim(file_root)//".1km"
      endif

      !--- add 'clavrx_' to the file name output
      file_root = 'clavrx_' // file_root
      
   
   end function file_root_from_l1b


end module LEVEL2_ROUTINES


