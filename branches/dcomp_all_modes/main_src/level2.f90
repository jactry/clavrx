
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
!      One public Subroutine: WRITE_PIXEL_HDF_RECORDS
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
   
   use CX_CONSTANTS_MOD, only: &
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
      , acha_version &
      , cloud_mask_version &
      , cloud_mask_thresholds_version
   
   ! - many of those variables are in level2_assign.inc   
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
      , gap_pixel_mask &
      , diag_pix_array_1 &
      , diag_pix_array_2 &
      , diag_pix_array_3 &
      , dir_level2 &
      , rtm_file_flag &
      , therm_cal_1b &
      , ref_cal_1b &
      , nav_opt &
      , use_sst_anal &
      , modis_clr_alb_flag &
      , nwp_opt &
      , scan_time &
      , line_idx_min_segment &
      , ref_ch1_min_3x3 &
      , ref_ch1_std_3x3 &
      , cldmask &
      , cld_type &
      , cld_phase &
      , zc_h2o &
      , zc_opaque_cloud &
      , tau_dcomp &
      , tau_dcomp_1 &
      , tau_dcomp_2 &
      , tau_dcomp_3 &
      , reff_dcomp & 
      , reff_dcomp_1 &
      , reff_dcomp_2 &
      , reff_dcomp_3 &
      , tau_dcomp_cost &
      , reff_dcomp_cost &
      , tau_dcomp_qf &
      , reff_dcomp_qf &
      , dcomp_quality_flag &
      , dcomp_info_flag &
      , insolation_dcomp &
      , insolation_diffuse_dcomp &
      , tau_nlcomp &
      , reff_nlcomp &
      , tau_nlcomp_cost &
      , reff_nlcomp_cost &
      , nlcomp_quality_flag &
      , nlcomp_info_flag &
      , cloud_063um_albedo &
      , cloud_063um_transmission_solar &
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
      , dcomp_success_fraction&
      , cdnc_dcomp &
      , aot_qf &
      , Bt_Ch31_Max_3x3 &
      , Bt_Ch20_median_3x3 &
      , olr &
      , Ndvi_Sfc_White_Sky &
      , Sst_Anal_Uni &
      , ozone_nwp_pix &
      , dust_mask &
      , smoke_mask &
      , shadow_mask &
      , fire_mask &
      , cost_aux &
      , hcld_dcomp &
      , CCl &
      , ems_ch20 &
      , ems_ch20_clear_rtm &
      , Ems_Ch20_Median_3x3 &
      , Btd_Ch31_Ch32_Bt_Ch31_Max_3x3 &
      , cld_phase_aux &
      , CLDMASK &
      , cld_type_aux &
      , zc_aux &
      , zc_co2 &
      , pc_top1_aux &
      , pc_top2_aux &
      , tc_opaque_cloud &
      , pc_uncertainty1_aux &
      , pc_uncertainty2_aux &
      , tau_aux &
      , reff_aux &
      , Beta_11um_85um_Tropo_Rtm &
      , Beta_11um_12um_Tropo_Rtm &
      , aot1 &
      , csv_file &
      , Covar_Ch27_Ch31_5x5 &
      , Ref_Max_chi1 &
      , Ref_min_chi1 &
      , Ref_mean_chi1 &
      , Ref_uni_chi1  &
      , Ref_Max_chi2 &
      , Ref_min_chi2 &
      , Ref_mean_chi2 &
      , Ref_uni_chi2 &
      , Ref_Max_chi3 &
      , Ref_min_chi3 &
      , Ref_mean_chi3 &
      , Ref_uni_chi3 &
      , Bt_Max_chi4 &
      , Bt_min_chi4 &
      , Bt_mean_chi4 &
      , Bt_uni_chi4 &
      , Bt_Max_chi5 &
      , Bt_min_chi5 &
      , Bt_mean_chi5 &
      , Bt_uni_chi5 &
      , Cld_Type_Ir &
      , Cld_Phase_Ir &
      , Tpw_Above_Cloud_Nwp_Pix &
      , Inversion_Strength_Nwp_Pix &
      , Inversion_Base_Nwp_Pix &
      , Inversion_Top_Nwp_Pix &
		, Bt_Ch31_Std_3x3
      
    
   use CLOUD_TYPE_BRIDGE_MODULE,only: &
      cloud_type_version
      
   use clavrx_message_module, only: &
       mesg

   use cx_hdf_write_mod, only:  &
      hdf_file_open &
      , create_sds &
      , compress_sds &
      , write_sds &
      , add_att &
      , close_sds &
      , close_file
   
   use cx_prd_mod, only: &
      prd_dtype &
      , prd_individual_dtype
   
   use cx_string_tools_mod, only: &
    is_numeric
    
   use cx_scale_tools_mod, only: &
      scale_i2_rank2 &
      , scale_i1_rank2
   
   use cx_muri_clavrx_bridge_mod
   
   implicit none
   
   private
   include 'hdf.f90'
   

   public:: WRITE_PIXEL_HDF_RECORDS


 
   !----------------------------------------------------------------------
   ! the following variables taken from process_avhrr_clavr
   !----------------------------------------------------------------------
   !--- hdf specific variables
   integer(kind=int4) :: Dim_Id
   integer(kind=int4),parameter :: SDS_RANK_1D = 1
   integer(kind=int4),dimension(1) :: Sds_Dims_1d
   integer(kind=int4),parameter :: SDS_RANK_2D = 2
   integer(kind=int4),dimension(SDS_RANK_2D) :: Sds_Dims_2d
   integer(kind=int4),dimension(SDS_RANK_2D) :: Sds_Start_2d
   integer(kind=int4),dimension(SDS_RANK_2D) :: Sds_Stride_2d
   integer(kind=int4),dimension(SDS_RANK_2D) :: Sds_Edge_2d
   integer,dimension(SDS_RANK_2D) :: Sds_Chunk_Size_2d

   character(len=11), private, parameter:: MOD_PROMPT = "LEVEL2:"
   character(len=18), private, parameter:: COORDINATES_STRING = "longitude latitude" 
   integer(kind=int4),  save:: Num_Scans_Level2_Hdf
  
   type(prd_dtype), save, target :: prd
   type(prd_individual_dtype), pointer:: prd_i
   
   integer :: id_file
   
   integer(kind=int4), parameter :: two_byte_max = 32767, & !(2**15)/2 - 1
                                        two_byte_min = -32767   !-(2**15)/2
                                        
   integer(kind=int4), parameter :: one_byte_max = 127, & !(2**8)/2 - 1
                                        one_byte_min = -127   !-(2**8)/2
   

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
   subroutine DEFINE_HDF_FILE_STRUCTURES( &
      Dir_Level2, &
      File_1b )
      
      implicit none

      character(len=*), intent(in):: Dir_Level2
      
      character(len=*), intent(in):: File_1b
      character(len=1020):: File_1b_root
      character(len=1020):: File_Level2
     
      
      integer:: erstat, istatus
      
     
 
      integer :: ii
      real(kind=real4):: Add_Offset
      real(kind=real4):: Scale_Factor
      integer ::  idx_1, idx_2 , idx_3
      character :: substring_1
      character(len=2) :: substring_2
      logical :: change_switch
      integer :: ch_switch
                                        
      File_1b_Root = file_root_from_l1b ( file_1b, &
         Sensor%Sensor_Name,Sensor%Spatial_Resolution_Meters)

      !
      sds_dims_2d = (/ Image%Number_Of_Elements, Image%Number_Of_Lines /)
      Sds_Chunk_Size_2d(1) = Image%Number_Of_Elements
      Sds_Chunk_Size_2d(2) = Image%Number_Of_Lines_Per_Segment

      !-- dimension of 1d variables
      Sds_Dims_1d(1) = Image%Number_Of_Lines
      
      File_Level2= trim(File_1b_Root)//".level2.hdf"
      call mesg (MOD_PROMPT//"creating level-2 file "//trim(file_Level2))
         
      id_file = hdf_file_open(trim(Dir_Level2)//trim(file_Level2), create=.true.)
      
      
      
      if (id_file < 0) then
         erstat = 68
         print *, EXE_PROMPT, MOD_PROMPT, "Level-2 file creation failed. Exiting..."
         stop 68
      end if
                               
      do ii = 1, prd % num_products
         prd_i => prd % product(ii)
         ! --- switch off if channel is needed but not set
         ! --- 
         change_switch = .false.
         
         ! - clavrx name indicates if this is a channel sensitive ouptput
         idx_1 = index(prd_i % name_clavrx,'Ch')
         idx_2 = index(prd_i % name_clavrx,'ch(')
         idx_3 = index(prd_i % name_clavrx,'ch')
                 
         if ( idx_1 .NE. 0) then
            substring_1 =  prd_i % name_clavrx(idx_1+2:idx_1+2)
            substring_2 =  prd_i % name_clavrx(idx_1+2:idx_1+3)
            change_switch = .true.
         else if ( idx_2 .NE. 0) then
            substring_1 =  prd_i % name_clavrx(idx_2+3:idx_2+3)
            substring_2 =  prd_i % name_clavrx(idx_2+3:idx_2+4)
            change_switch = .true.
         else if ( idx_3 .NE. 0) then
            substring_1 =  prd_i % name_clavrx(idx_3+2:idx_3+2)
            substring_2 =  prd_i % name_clavrx(idx_3+2:idx_3+3)
            change_switch = .true.
         end if
         
         if (change_switch) then
             
            if (is_numeric(substring_2) .and. substring_2(2:) .ne. ')') then              
               read (substring_2,'(I2)') ch_switch
               if ( .NOT. sensor % chan_on_flag_default(ch_switch) ) prd_i % switch = .false.   
            else if (is_numeric(substring_1)) then             
               read (substring_1,'(I1)') ch_switch
               if ( .NOT. sensor % chan_on_flag_default(ch_switch) ) prd_i % switch = .false.            
            end if
            
         end if
         
         if (trim(prd_i % name_clavrx) .EQ. '"NOT_SET_YET"') prd_i % switch = .false.
         
         if (index(prd_i % name,'iband') .NE. 0 .and. trim(sensor % sensor_name) .NE. 'VIIRS') prd_i % switch = .false. 
         
         if ( prd_i % switch ) then
            
            select case ( prd_i % dim)
            case(1)
              prd_i % Sds_Id= create_sds (id_file, prd_i % name ,  Image%Number_Of_Lines, prd_i % dtype)
            case(2)
               prd_i % Sds_Id= create_sds (id_file, prd_i % name , sds_dims_2d , prd_i % dtype)
               istatus = compress_sds ( prd_i % Sds_Id,Compress_Flag, Sds_Chunk_Size_2d) 
            end select

            call add_att(  prd_i % Sds_Id, 'SCALED', prd_i % scaling)
            call add_att( prd_i % sds_id, 'unit', trim(prd_i % unit)) 
            call add_att( prd_i % sds_id, 'standard_name', trim(prd_i % standard_name))
            call add_att( prd_i % sds_id, 'long_name', trim(prd_i % long_name))  
            
            if ( prd_i % scaling .gt. 0 ) then
            
               select case ( prd_i % dtype ) 
               case(1)
                  scale_factor = (prd_i%act_max - prd_i%act_min )/(one_byte_max - one_byte_min)
                  add_offset = prd_i%act_min - scale_factor * one_byte_min
                  call add_att ( prd_i % sds_id, 'scaled_missing',Missing_Value_Int1 )
               case(2)
                  scale_factor = (prd_i%act_max - prd_i%act_min )/(two_byte_max - two_byte_min)
                  add_offset = prd_i%act_min - scale_factor * two_byte_min 
                  call add_att ( prd_i % sds_id, 'scaled_missing',Missing_Value_Int2 )  
               end select   
                          
               call add_att ( prd_i % sds_id, 'add_offset', add_offset)
               call add_att ( prd_i % sds_id, 'scale_factor',scale_factor)
               call add_att ( prd_i % sds_id, 'actual_min',prd_i%act_min )
               call add_att ( prd_i % sds_id, 'actual_max',prd_i%act_max )
            end if
         end if
      end do
 
   end subroutine DEFINE_HDF_FILE_STRUCTURES

   !====================================================================
   ! SUBROUTINE Name: WRITE_PIXEL_HDF_RECORDS
   !
   ! Function:
   !   Writes output to appropriate Level 2 files
   !
   !   call from process_clavrx
   !
   ! Description:
   !   This subroutine, given the flags that determine which files are 
   !   being created, outputs the various pixel level outputs to the
   !   appropriate level 2 files and appropriate SDSs for a given segment
   !   this program is called in process_clavrx
   !
   !====================================================================
   subroutine WRITE_PIXEL_HDF_RECORDS(segment_number)
      implicit none
      integer, intent(in) :: segment_number
      
      integer:: Istatus
    
      
      integer(kind=int1), allocatable :: data_dim1_dtype1(:)
      integer(kind=int4), allocatable :: data_dim1_dtype3(:)
      real(kind=real4), allocatable :: data_dim1_dtype4(:)
      integer(kind = int1), allocatable :: data_dim2_dtype1(:,:)
      integer( kind = int2), allocatable :: data_dim2_dtype2(:,:)
      integer(kind=int4), allocatable :: data_dim2_dtype3(:,:)
      real(kind=real4), allocatable :: data_dim2_dtype4(:,:)
      real(kind=real4), allocatable :: data_dim1_dtype_r4(:)
      real(kind=real4), allocatable :: data_dim2_dtype_r4(:,:)
     
      integer(kind=int2), allocatable :: Two_Byte_dummy(:,:)
      integer(kind=int1), allocatable :: One_Byte_dummy(:,:)
      integer :: ii
      character (len =50) :: name
      
      ! first segment
      if (Segment_Number == 1) then

         !--- place algorithm svn tags into global strings for output find better place for this..
        
         !- read csv file        
         if (.not. prd % is_set) call prd % read_products(csv_file)
         
         !- initialize  
         Num_Scans_Level2_Hdf = 0
         
         call DEFINE_HDF_FILE_STRUCTURES( &
                     Dir_Level2, &
                     Image%Level1b_Name)

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

      if (Sds_Edge_2d(2) .LE. 0) then
         return
      end if

      !--- update Num_Scans_Level2_Hdf
      Num_Scans_Level2_Hdf = min(Image%Number_Of_Lines,Num_Scans_Level2_Hdf +  &
                               Image%Number_Of_Lines_Read_This_Segment)

      !-------------------------------------------------------------------------
      ! write to level2 file
      !-------------------------------------------------------------------------
      
         istatus = 0
            ! ---  re-read in products from csv file
         
         !!! - all this data are used in level2_assign.inc file
         allocate ( data_dim1_dtype1(sds_edge_2d(2)))
         allocate ( data_dim1_dtype3(sds_edge_2d(2)))
         allocate ( data_dim1_dtype4(sds_edge_2d(2)))
         allocate ( data_dim2_dtype1(sds_edge_2d(1),sds_edge_2d(2))) 
         allocate ( data_dim2_dtype2(sds_edge_2d(1),sds_edge_2d(2)))
         allocate ( data_dim2_dtype3(sds_edge_2d(1),sds_edge_2d(2)))
         allocate ( data_dim2_dtype4(sds_edge_2d(1),sds_edge_2d(2)))
         allocate ( two_byte_dummy(sds_edge_2d(1),sds_edge_2d(2)))   
         allocate ( one_byte_dummy(sds_edge_2d(1),sds_edge_2d(2)))  
         allocate ( data_dim1_dtype_r4 (sds_edge_2d(2)))
         allocate ( data_dim2_dtype_r4 (sds_edge_2d(1),sds_edge_2d(2)))
        
        
       
         do ii = 1, prd % num_products
            prd_i => prd % product(ii)  
            name =  prd_i % name         
             
            if ( prd_i % switch ) then
               include 'level2_assign.inc'
               
               select case (prd_i % dim)
                  
               case(1)               
                  select case (prd_i % dtype)
                  case(1)
                  Istatus = write_sds ( prd_i % sds_id, Sds_Start_2d(2), Sds_Stride_2d(2),          &
                         Sds_Edge_2d(2), data_dim1_dtype1 ) + Istatus
                  
                  case(3)                       
                  Istatus = write_sds ( prd_i % sds_id,Sds_Start_2d(2), Sds_Stride_2d(2),          &
                         Sds_Edge_2d(2), data_dim1_dtype3 ) + Istatus
                  
                  case(4)
                  Istatus = write_sds ( prd_i % sds_id,Sds_Start_2d(2), Sds_Stride_2d(2),          &
                         Sds_Edge_2d(2), data_dim1_dtype4 ) + Istatus
               
                  end select  
               
               case(2)
                  select case (prd_i % dtype)
                  case(1)
                  if (prd_i % scaling == 1 ) then                     
                     call scale_i1_rank2(data_dim2_dtype_r4  &
                        ,prd_i % act_min,prd_i % act_max,Missing_Value_Real4 &
                        ,One_Byte_dummy)
                     Istatus = write_sds ( prd_i % sds_id,Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d, &
                        One_Byte_Dummy ) + Istatus 
                  else
                     
                     Istatus = write_sds ( prd_i % sds_id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                               data_dim2_dtype1) + Istatus
                  end if
                  
                  case(2)
                  if (prd_i % scaling == 1) then
                     call scale_i2_rank2(data_dim2_dtype_r4  &
                        ,prd_i % act_min,prd_i % act_max,Missing_Value_Real4 &
                        ,Two_Byte_dummy)
                     Istatus = write_sds ( prd_i % sds_id,Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d, &
                        Two_Byte_Dummy ) + Istatus
                  else 
                     Istatus = write_sds ( prd_i % sds_id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                               data_dim2_dtype2) + Istatus
                  end if 
                  
                  case(4)
                  !- this is only diagnostic variables, we don't scale..
                  Istatus = write_sds ( prd_i % sds_id, Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
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
         if ( allocated ( one_byte_dummy)) deallocate (one_byte_dummy)
                  
         !--- check for and report errors         
         if (Istatus /= 0) then
            print *, EXE_PROMPT, MOD_PROMPT, "Error writing to level2 file: ", Istatus
            stop
         endif
         
         ! final segment, write global attributes and closing.. bye bye
         if ( segment_number .EQ. Image%Number_Of_Segments) then 
           
            call add_global_attributes ( id_file) 
            do ii = 1, prd % num_products
                prd_i => prd % product(ii)
                if ( prd_i % switch ) then
                  call close_sds (  prd_i % sds_id)
                end if  
            end do
            call close_file (id_file)      
         end if
    
   end subroutine WRITE_PIXEL_HDF_RECORDS

   ! ----------------------------------------------------------------
   !
   ! ----------------------------------------------------------------
   subroutine add_global_attributes ( id_file)
       
       use pixel_common, only: &
         sensor &
         , image &
         , acha &
         , dcomp_mode
       
      integer, intent(in) :: id_file
      character(len=36) :: machine
      character(len=36) :: user
     
      call getenv ("HOST",machine)
      call getenv ("USER",user)
      
      call add_att (id_file, 'MACHINE',trim(machine))
      call add_att (id_file, 'date_created',hdf_timestamp())
      call add_att (id_file, 'user_created',trim(user))
      call add_att (id_file, 'institution','CIMSS')
      call add_att (id_file, 'sensor',trim(sensor% sensor_name))
      call add_att (id_file, 'satellite', trim(sensor%platform_name))
      call add_att (id_file,  'L1B', trim(Image%Level1b_Name))
      !call add_att (id_file,   'RESOLUTION_KM',
      !call add_att (id_file,   'spatial_resolution'
      call add_att (id_file,  'START_YEAR',image%start_year)
      call add_att (id_file,  'START_DAY_OF_YEAR',image%start_doy)
      call add_att (id_file,  'START_TIME_FRACTION_DAY',image%start_time/3600000.0)
      call add_att (id_file,  'END_YEAR',image%end_year)
      call add_att (id_file,  'END_DAY_OF_YEAR',image%end_doy)
      call add_att (id_file,  'END_TIME',image%end_time/3600000.0)
      call add_att (id_file,  'ACHA_MODE',acha%mode)
      call add_att (id_file,  'DCOMP_MODE',dcomp_mode)
      call add_att (id_file,  'WMO_SATELLITE_CODE',sensor%wmo_id)
      !call add_att (id_file,  'REFL_0_65UM_NOM_DARK_COMPOSITE_NAME'
      !call add_att (id_file,  'NAIVE_BAYESIAN_CLOUD_MASK_NAME', 
      
      call add_att ( id_file,'CLOUD_MASK_VERSION', trim(cloud_mask_version))
      call add_att ( id_file,'CLOUD_MASK_THRESHOLDS_VERSION',trim(CLOUD_MASK_THRESHOLDS_VERSION))
      call add_att ( id_file,'CLOUD_TYPE_VERSION',trim(CLOUD_TYPE_VERSION))
      call add_att ( id_file,'ACHA_VERSION',trim(ACHA_VERSION))
      call add_att ( id_file,'DCOMP_VERSION',trim(DCOMP_VERSION))
      
      call add_att ( id_file, 'NUMBER_OF_ELEMENTS', Image%Number_Of_Elements )
      call add_att ( id_file, 'NUMBER_OF_SCANS_LEVEL1B', Image%Number_Of_Lines) 
      call add_att ( id_file, 'NUMBER_OF_SCANS_LEVEL2', Num_Scans_Level2_Hdf)
    
    !  Istatus = sfsnatt(Sd_Id_Level2, "PROCESSING_TIME_MINUTES", DFNT_FLOAT32,1,Orbital_Processing_Time_Minutes)+Istatus
    !  Istatus = sfsnatt(Sd_Id_Level2, "NONCONFIDENT_CLOUD_MASK_FRACTION", DFNT_FLOAT32,1,NONCONFIDENT_CLOUD_MASK_Fraction)+Istatus
    !  Istatus = sfsnatt(Sd_Id_Level2, "ACHA_SUCCESS_FRACTION", DFNT_FLOAT32,1,ACHA%Success_Fraction)+Istatus
    !  Istatus = sfsnatt(Sd_Id_Level2, "DCOMP_SUCCESS_FRACTION", DFNT_FLOAT32,1,DCOMP_Success_Fraction)+Istatus
      
      
      contains
      function hdf_timestamp() result(string)
         character (len = 36) ::string
         character(len=10), dimension(3):: ctime

         call date_and_time(ctime(1), ctime(2), ctime(3))
         ! Timestamp string format in accordance with the ISO-8601 standard.
         string = ctime(1)(1:4)//"-"//ctime(1)(5:6)//"-"//ctime(1)(7:8)&
                //"T"//ctime(2)(1:2)//":"//ctime(2)(3:4)//":"//ctime(2)(5:6)&
                //ctime(3)(1:3)//":"//ctime(3)(4:5)
         return
      end function hdf_timestamp 
        
   end subroutine add_global_attributes
   
   !============================================================================
   !
   !   returns filename root for level-2 output from level-1b filename
   !
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
      if (trim(l1b_ext) == '.hdf') then
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
      if(trim(Sensor) == 'GOES-IL-IMAGER' .or.  &
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
