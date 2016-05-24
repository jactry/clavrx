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
!--------------------------------------------------------------------------------------
module LEVEL2_ROUTINES
   use CONSTANTS
   use PIXEL_COMMON
   use HDF
   use SCALING_PARAMETERS
   use HDF_PARAMS
   
   use clavrx_message_module
   
   use csv_mod
   use strings

   implicit none
   private

   public:: DEFINE_HDF_FILE_STRUCTURES, &
            WRITE_PIXEL_HDF_RECORDS, &
            CLOSE_PIXEL_HDF_FILES, &
            WRITE_ALGORITHM_ATTRIBUTES

   private::DEFINE_PIXEL_2D_SDS, &
            DEFINE_PIXEL_3D_SDS

 !--- rtm indices
 integer, parameter, private:: Num_Rtm_Sds = 28
 integer, private, save:: Sd_Id_Rtm
 integer(kind=int4), dimension(Num_Rtm_Sds), save, private:: Sds_Id_Rtm

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

 INCLUDE 'level2.inc'
 
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
 character(len=4):: l1b_ext
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
  integer   ( kind = 4 ) csv_file_status
  integer   ( kind = 4 ) csv_file_unit
  integer   ( kind = 4 ) csv_record_status
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) line_num
  character ( len = 400 ) record
  integer   ( kind = 4 ) value_count
  character (len =300) :: before
   character ( len = 300) :: rec_arr ( 13 )
  integer :: j
  logical :: switch
  integer :: var_dim
  integer(kind = 1) :: scaling
  integer :: dtype
   
   
  !--- begin executable code
  blank_int4 = 0
  blank_real4 = 0.0
  blank_char = " "

  !--- make a file name base for output
  File_1b_Root = trim (file_1b)
  
  
  ! ---
  csv_file_name='clavrx_level2_products.csv'
  call csv_file_line_count ( csv_file_name, line_num )
   print*,'CSV Line number: ',line_num
   
   
  

  !--- special processing for modis - remove hdf suffix
  l1b_ext = File_1b_Root(len_trim(File_1b_Root)-3:len_trim(File_1b_Root))
  if (trim(l1b_ext) == ".hdf") then
    File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-4)
  endif

  !--- special processing for viirs - remove hdf suffix - this hard coded for
  if (trim(Sensor%Sensor_Name) == 'VIIRS') then
    File_1b_Root = File_1b_Root(7:len_trim(File_1b_Root)-34)
  endif

  !--- special processing for ahi - remove B01.nc suffix - this hard coded for
  if (trim(Sensor%Sensor_Name) == 'AHI') then
    File_1b_Root = File_1b_Root(4:len_trim(File_1b_Root)-12)
  endif

  !--- special processing for IFF - remove hdf suffix - this hard coded for
! !PEATE files
! if (trim(Sensor%Sensor_Name) == 'AQUA-IFF' .or. trim(Sensor%Sensor_Name) == 'AQUA-IFF') then
!   File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-29)
! elseif (trim(Sensor%Sensor_Name) == 'AVHRR-IFF') then
!   File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-20)
! endif

  !--- do this for GOES names which are named goesxx_1_year_jday_hour.area
  if (trim(Sensor%Sensor_Name) == 'GOES-IL-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'GOES-MP-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'COMS-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'MTSAT-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'FY2-IMAGER' .or.  &
      trim(Sensor%Sensor_Name) == 'SEVIRI') then

    !-- remove area suffix
    l1b_ext = File_1b_Root(len_trim(File_1b_Root)-3:len_trim(File_1b_Root))
    if (trim(l1b_ext) == "area") then
     File_1b_Root = File_1b_Root(1:len_trim(File_1b_Root)-5)
    endif
    !-- remove channel number
    if (trim(l1b_ext) == "area") then
        ipos = index(File_1b_Root, "_1_")
        ilen = len(File_1b_Root)
        File_1b_Root = File_1b_Root(1:ipos-1) // "_"//File_1b_Root(ipos+3:ilen)
    endif
  endif

  !--- add 1km tag for 1km GOES files
  if (index(Sensor%Sensor_Name,'GOES') > 0 .and. Sensor%Spatial_Resolution_Meters == 1000) then
    File_1b_Root = trim(File_1b_Root)//".1km"
  endif

  !--- add 'clavrx_' to the file name output
  File_1b_Root = 'clavrx_' // File_1b_Root

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
     endif

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

     !-- scan line variables
    
     call csv_file_open_read ( csv_file_name, csv_file_unit )
      do i = 1, line_num
         
         read ( csv_file_unit, '(a)', iostat = csv_file_status ) record
         call csv_value_count ( record, csv_record_status, value_count )
         rec_arr = extract_single ( trim ( record ) )
        
         switch = trim(rec_arr(1))  .eq. "1"
         
         if ( i == 1) cycle
         read( rec_arr(2),  * ) var_dim
         read ( rec_arr(4), * ) dtype
         read ( rec_arr(5), * ) scaling
        
         print*,i,trim(rec_arr(3))
         
         if ( switch ) then
            select case (var_dim)
            case ( 1 )
               
               Sds_Id_Level2(i) = sfcreate(Sd_Id_Level2,trim(rec_arr(3)),DFNT_INT32,Sds_Rank_1d,Sds_Dims_1d)    
               Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "SCALED", DFNT_INT8, 1, scaling) + Istatus_Sum  
               print*,'1: ',istatus_sum
               Istatus_Sum = sfscatt(Sds_Id_Level2(i), "units", DFNT_CHAR8, 4,trim(rec_arr(13) )) + Istatus_Sum 
               Istatus_Sum = sfscatt(Sds_Id_Level2(i), "standard_name", DFNT_CHAR8, 13, trim(rec_arr(12) )) + Istatus_Sum 
               Istatus_Sum = sfscatt(Sds_Id_Level2(i), "long_name", DFNT_CHAR8,  &
                    len_trim(trim(rec_arr(14))),  trim(rec_arr(14)) ) + Istatus_Sum
               Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "RANGE_MISSING", DFNT_FLOAT32, &
                    1,real(Missing_Value_Int4,kind=real4)) + Istatus_Sum
                   
               Istatus_Sum = sfsnatt(Sds_Id_Level2(i), "_FillValue", DFNT_INT32, &
                    1,Missing_Value_Int4) + Istatus_Sum     
                 
            case(2)
               
               call DEFINE_PIXEL_2D_SDS(Sds_Id_Level2(i),Sd_Id_Level2,Sds_Dims_2d,Sds_Chunk_Size_2d,&
                               trim(rec_arr(3)), &
                               trim(rec_arr(12)), &
                               trim(rec_arr(14)), &
                               DFNT_INT8, scaling, 0.0, 1.0, &
                               trim(rec_arr(13) ), Real(Missing_Value_Int1,kind=real4), Istatus)
               Istatus_Sum = Istatus_Sum + Istatus
               
            end select
         end if
         
         print*,i,istatus_sum
        
        
     
      end do
     
   



     !--- check for and report errors
     if (Istatus_Sum /= 0) then
       print *, EXE_PROMPT, MOD_PROMPT, "Error defining sds in level2 hdf file"
       stop
     endif
  endif
  
  contains

   function extract_single ( record) result  (record_single)
      character ( len = * ) :: record
      character ( len = 300) :: record_single ( 13 )
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
  
  

end subroutine DEFINE_HDF_FILE_STRUCTURES

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
!
!====================================================================
subroutine WRITE_PIXEL_HDF_RECORDS(Rtm_File_Flag,Level2_File_Flag)

 integer, intent(in):: Rtm_File_Flag
 integer, intent(in):: Level2_File_Flag
 integer:: Istatus
 integer:: Line_Idx

! HDF function declarations
 integer:: sfwdata

!-----------------------------------------------------------------------
! Get time of each scan line and convert to scale
!-----------------------------------------------------------------------
  where(Scan_Time == Missing_Value_Real4)
            Utc_Scan_Time_Hours = Missing_Value_Real4
  elsewhere
            Utc_Scan_Time_Hours = Scan_Time/60.0/60.0/1000.0
  endwhere

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
    endif

    !--- update Num_Scans_Level2_Hdf
    Num_Scans_Level2_Hdf = min(Image%Number_Of_Lines,Num_Scans_Level2_Hdf +  &
                               Image%Number_Of_Lines_Read_This_Segment)

!-------------------------------------------------------------------------
! write to level2 file
!-------------------------------------------------------------------------
   if (Level2_File_Flag == sym%YES) then

      Istatus = 0
      
      !--- scan number
      if (Sds_Num_Level2_Scanline_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Scanline), Sds_Start_2d(2), Sds_Stride_2d(2),          &
                         Sds_Edge_2d(2), scan_number(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- scan time
      if (Sds_Num_Level2_Time_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Time), Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                         Utc_Scan_Time_Hours(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
 
      !--- bad scan flag
      if (Sds_Num_Level2_Bad_Scan_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bad_Scan), Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                         Bad_Scan_Flag(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- ascending/descending flag
      if (Sds_Num_Level2_Asc_Flag_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_asc_flag), Sds_Start_2d(2), Sds_Stride_2d(2), Sds_Edge_2d(2),  &
                 Nav%Ascend(Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Bad Pixel Mask
      if (Sds_Num_Level2_Bad_Pixel_Mask_Flag == sym%YES) then
       One_Byte_Temp = Bad_Pixel_Mask
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bad_Pixel_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Gap Pixel Mask
      if (Sds_Num_Level2_Gap_Pixel_Mask_Flag == sym%YES) then
       One_Byte_Temp = Gap_Pixel_Mask
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Gap_Pixel_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- diagnostic field #1
      if (Sds_Num_Level2_Diag1_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Diag1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                 Diag_Pix_Array_1(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- diagnostic field #2
      if (Sds_Num_Level2_Diag2_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Diag2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                 Diag_Pix_Array_2(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- diagnostic field #3
      if (Sds_Num_Level2_Diag3_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Diag3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d,  &
                 Diag_Pix_Array_3(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif


      !--- packed pixel metadata
      if (Sds_Num_Level2_Meta_Data_Flag == sym%YES) then
       One_Byte_Temp = 0
       Temp_Mask = 0
       do Line_Idx = 1, Image%Number_Of_Lines_Per_Segment
         Temp_Mask(:,Line_Idx) = Sensor%Chan_On_Flag_Per_Line(6,Line_Idx)
       enddo
       One_Byte_Temp = ishft(Bayes_Mask_Sfc_Type_Global,3) + ishft(Temp_Mask,2)+ &
                       ishft(solar_contamination_mask,1) + bad_pixel_mask
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Meta_Data), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- latitude
      if (Sds_Num_Level2_Lat_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Nav%Lat,sym%LINEAR_SCALING,Min_Lat,Max_Lat,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Lat),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d, &
                  Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- longitude
      if (Sds_Num_Level2_Lon_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Nav%Lon,sym%LINEAR_SCALING,Min_Lon,Max_Lon,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_lon),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                   Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- latitude parallax_corrected
      if (Sds_Num_Level2_Latpc_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Nav%Lat_Pc,sym%LINEAR_SCALING,Min_Lat,Max_Lat,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Latpc),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d, &
                  Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- longitude parallax corrected
      if (Sds_Num_Level2_Lonpc_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Nav%Lon_Pc,sym%LINEAR_SCALING,Min_Lon,Max_Lon,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Lonpc),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                   Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- sensor zenith
      if (Sds_Num_Level2_Zen_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Satzen,sym%LINEAR_SCALING,Min_Zen,Max_Zen,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_zen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- solar zenith
      if (Sds_Num_Level2_Solzen_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Solzen,sym%LINEAR_SCALING,Min_Solzen,Max_Solzen,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Solzen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- relative azimuth 
      if (Sds_Num_Level2_Relaz_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Relaz,sym%LINEAR_SCALING,Min_Relaz,Max_Relaz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Relaz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- solar azimuth
      if (Sds_Num_Level2_Solaz_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Solaz,sym%LINEAR_SCALING,Min_Solaz,Max_Solaz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Solaz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- sensor/satellite azimuth
      if (Sds_Num_Level2_Sataz_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Sataz,sym%LINEAR_SCALING,Min_Sataz,Max_Sataz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Sataz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif


      !--- glint zenith
      if (Sds_Num_Level2_Glintzen_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Glintzen,sym%LINEAR_SCALING,Min_Glintzen,Max_Glintzen,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Glintzen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- scattering zenith
      if (Sds_Num_Level2_Scatzen_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Scatangle,sym%LINEAR_SCALING,Min_Scatang,Max_Scatang,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Scatzen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- lunar zenith
      if (Sds_Num_Level2_Lunzen_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Lunzen,sym%LINEAR_SCALING,Min_Solzen,Max_Solzen,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Lunzen),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- lunar relative azimuth 
      if (Sds_Num_Level2_LunRelaz_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%LunRelaz,sym%LINEAR_SCALING,Min_Relaz,Max_Relaz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_LunRelaz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- lunar azimuth 
      if (Sds_Num_Level2_Lunaz_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(44) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Geo%Lunaz,sym%LINEAR_SCALING,Min_Solaz,Max_Solaz,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Lunaz),Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,  &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- packed land cover (land,snow,coast masks)
      if (Sds_Num_Level2_Packed_Land_Flag == sym%YES) then
       One_Byte_Temp = 0
       One_Byte_Temp = ishft(Sfc%Land,5) + ishft(Sfc%Snow,3) + Sfc%Coast_Mask
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Packed_Land), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- glint mask
      if (Sds_Num_Level2_Glint_Mask_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Glint_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Glint_Mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- coast mask
      if (Sds_Num_Level2_Coast_Mask_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Coast_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Coast_Mask(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- surface type
      if (Sds_Num_Level2_Sfc_Type_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Sfc_Type), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Sfc_Type(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- land classification
      if (Sds_Num_Level2_Land_Mask_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Land_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Land(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Snow classification
      if (Sds_Num_Level2_Snow_Mask_Flag == sym%YES) then
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Snow_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Sfc%Snow(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- surface elevation
      if (Sds_Num_Level2_Zsfc_Flag == sym%YES) then
          call SCALE_VECTOR_I2_RANK2(Sfc%Zsfc,sym%LINEAR_SCALING,Min_Zsfc,Max_Zsfc,Missing_Value_Real4,Two_Byte_Temp)
          Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Zsfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                    Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

!--------------------------------------------------------------------------------------------------
!--- observations
!--------------------------------------------------------------------------------------------------

      !-- Ch1 reflectance
      if (Sds_Num_Level2_Ch1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(ch(1)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

      endif

      !-- Ch2 reflectance
      if (Sds_Num_Level2_Ch2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(2) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(2)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch3 reflectance
      if (Sds_Num_Level2_Ch3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(3) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(3)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch3,Max_Ref_Ch3,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch4 reflectance
      if (Sds_Num_Level2_Ch4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(4) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(4)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch4,Max_Ref_Ch4,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch5 reflectance
      if (Sds_Num_Level2_Ch5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(5) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(5)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch5,Max_Ref_Ch5,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch6 reflectance
      if (Sds_Num_Level2_Ch6_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(6)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch6,Max_Ref_Ch6,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch6), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch7 reflectance
      if (Sds_Num_Level2_Ch7_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(7) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(7)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch7,Max_Ref_Ch7,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch7), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch8 reflectance
      if (Sds_Num_Level2_Ch8_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(8) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(8)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch8,Max_Ref_Ch8,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch8), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch9 reflectance
      if (Sds_Num_Level2_Ch9_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(9) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(9)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch9,Max_Ref_Ch9,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch9), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch17 reflectance
      if (Sds_Num_Level2_Ch17_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(17) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(17)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch17,Max_Ref_Ch17,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch17), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch18 reflectance
      if (Sds_Num_Level2_Ch18_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(18) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(18)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch18,Max_Ref_Ch18,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch18), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch19 reflectance
      if (Sds_Num_Level2_Ch19_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(19) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(19)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch19,Max_Ref_Ch19,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch19), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- ch26 reflectance
      if (Sds_Num_Level2_Ch26_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(26) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(26)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch26,Max_Ref_Ch26,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch26), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- chDNB reflectance
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES .and. Sds_Num_Level2_ChDNB_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(44)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_ChDNB,Max_Ref_ChDNB,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_ChDNB), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- chDNB reflectance lunar
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES .and. Sds_Num_Level2_ChDNB_lunar_Flag == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(44)%Ref_Lunar_Toa,sym%LINEAR_SCALING,Min_Ref_ChDNB_lunar,Max_Ref_ChDNB_lunar,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_ChDNB_lunar), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !-- Ch20 reflectance
      if (Sds_Num_Level2_Ch20_Ref_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(20)%Ref_Toa,sym%LINEAR_SCALING,Min_Ref_Ch20,Max_Ref_Ch20,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Ref), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
      !--- Ch20 temperature
      if (Sds_Num_Level2_Ch20_Bt_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(20)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt20,Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Bt), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
      !--- Ch22 temperature
      if (Sds_Num_Level2_Ch22_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(22) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(22)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt22,Max_Bt22,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch22), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
      !--- Ch27 temperature
      if (Sds_Num_Level2_Ch27_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(27)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt27,Max_Bt27,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch27), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif
      !--- Ch28 temperature
      if (Sds_Num_Level2_Ch28_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(28) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(28)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt28,Max_Bt28,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch28), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch29 temperature
      if (Sds_Num_Level2_Ch29_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(29) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(29)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt29,Max_Bt29,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch29), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch30 temperature
      if (Sds_Num_Level2_Ch30_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(30) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(30)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt30,Max_Bt30,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch30), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch31 temperature
      if (Sds_Num_Level2_Ch31_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(31)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt31,Max_Bt31,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch32 temperature
      if (Sds_Num_Level2_Ch32_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(32)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt32,Max_Bt32,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch32), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch33 temperature
      if (Sds_Num_Level2_Ch33_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(33) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(33)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt33,Max_Bt33,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch33), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch34 temperature
      if (Sds_Num_Level2_Ch34_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(34) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(34)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt34,Max_Bt34,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch34), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch35 temperature
      if (Sds_Num_Level2_Ch35_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(35) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(35)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt35,Max_Bt35,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch35), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch36 temperature
      if (Sds_Num_Level2_Ch36_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(36) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(36)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt36,Max_Bt36,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch36), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch37 temperature
      if (Sds_Num_Level2_Ch37_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(37) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(37)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt37,Max_Bt37,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch37), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch38 temperature
      if (Sds_Num_Level2_Ch38_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(38)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt38,Max_Bt38,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch38), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch45 temperature 13um Pseudo
      if (Sds_Num_Level2_Ch45_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(45) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Ch(45)%Bt_Toa,sym%LINEAR_SCALING,Min_Bt45,Max_Bt45,&
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch45), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch1_Min_3x3
      if (Sds_Num_Level2_Ch1_Min_3x3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Ref_Ch1_Min_3x3,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1, &
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Min_3x3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Ch1_Std_3x3
      if (Sds_Num_Level2_Ch1_Std_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Ref_Ch1_Std_3x3,sym%LINEAR_SCALING,Min_Ref_Ch1_std,Max_Ref_Ch1_std, &
                                 Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Std), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Bt_Ch31_std_3x3
      if (Sds_Num_Level2_Ch31_Std_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Bt_Ch31_Std_3x3,sym%LINEAR_SCALING,Min_Bt31_std,Max_Bt31_std, &
                                 Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Std), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                      One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif

      !--- Bt_Ch31_Max_3x3
      if (Sds_Num_Level2_Ch31_Max_3x3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(Bt_Ch31_Max_3x3,sym%LINEAR_SCALING,Min_Bt31,Max_Bt31, &
                                 Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Max_3x3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                      Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
      endif


!-------- cloud properties

     !--- cloud probability
     if (Sds_Num_Level2_Cldprob_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(posterior_cld_probability, &
                                 sym%LINEAR_SCALING,Min_frac,Max_frac,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cldprob), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld mask
     if (Sds_Num_Level2_Cld_Mask_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Mask), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        cld_mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- adjacent pixel cloud mask
     if (Sds_Num_Level2_Adj_Pix_Cld_Mask_Flag == sym%YES) then
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Adj_Pix_Cld_Mask), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Adj_Pix_Cld_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld mask test vector (first byte - acm only)
     if (Sds_Num_Level2_Cld_Tests_Flag == sym%YES) then     

       Sds_Start_3d(1) = 0
       Sds_Start_3d(2) = 0
       Sds_Start_3d(3) = Sds_Start_2d(2)

       Sds_Stride_3d(1) = 1
       Sds_Stride_3d(2) = 1
       Sds_Stride_3d(3) = 1

       Sds_Edge_3d(1) = Max_Num_Cld_Test_Bytes
       Sds_Edge_3d(2) = Sds_Edge_2d(1)
       Sds_Edge_3d(3) = Sds_Edge_2d(2)

      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Tests),  &
                        Sds_Start_3d, &
                        Sds_Stride_3d, &
                        Sds_Edge_3d, &
                        Cld_Test_Vector_Packed(:,:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- bayes mask sfc type
     if (Sds_Num_Level2_Bayes_Sfc_Type_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bayes_Sfc_Type), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Bayes_Mask_Sfc_Type_Global(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- dust mask
     if (Sds_Num_Level2_Dust_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Dust), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Dust_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- smoke mask
     if (Sds_Num_Level2_Smoke_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Smoke), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Smoke_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- shadow mask
     if (Sds_Num_Level2_Shadow_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Shadow), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Shadow_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- fire mask
     if (Sds_Num_Level2_Fire_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Fire), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Fire_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- csbt for ch27
     if (Sds_Num_Level2_Ch27_CSBT_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch27_CSBT), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        ch(27)%CSBT_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- csbt for ch28
     if (Sds_Num_Level2_Ch28_CSBT_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(28) == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch28_CSBT), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        ch(28)%CSBT_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- csbt for ch37
     if (Sds_Num_Level2_Ch37_CSBT_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(37) == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch37_CSBT), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        ch(37)%CSBT_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- csbt for ch33
     if (Sds_Num_Level2_Ch33_CSBT_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(33) == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch33_CSBT), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        ch(33)%CSBT_Mask(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth from mask
     if (Cld_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES .and. Sds_Num_Level2_Cod_Mask_Flag == sym%YES) then

      Temp_Pix_Array_1 = Ch(1)%Opd
      if (Sensor%Chan_On_Flag_Default(44) == sym%YES) then
        where(Geo%Solzen /= Missing_Value_Real4 .and. Geo%Solzen > 90.0)
            Temp_Pix_Array_1 = Ch(44)%Opd
        endwhere
      endif

      call SCALE_VECTOR_I1_RANK2(Temp_Pix_Array_1,sym%LINEAR_SCALING,Min_Tau,Max_Tau,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Mask), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

     endif

     !--- cld type
     if (Sds_Num_Level2_Cld_Type_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Type), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Type(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld phase
     if (Sds_Num_Level2_Cld_Phase_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Phase), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Phase(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld phase aux
     if (Sds_Num_Level2_Cld_Phase_Aux_Flag == sym%YES) then
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Phase_Aux), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Phase_Aux(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- auxiliary cld mask
     if (Sds_Num_Level2_Cld_Mask_Aux_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Mask_Aux), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Mask_Aux(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- auxiliary cld type
     if (Sds_Num_Level2_Cld_Type_Aux_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Type_Aux), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Type_Aux(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- auxiliary cld height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Aux_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Zc_Aux,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Aux), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- sndr cld height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Sndr_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Zc_Co2,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Sndr), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld pressure
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Pc,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_ctp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld temperature
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctt_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Tc,sym%LINEAR_SCALING,Min_Tc,Max_Tc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctt), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Zc,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_cth), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld top height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Top_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Zc_Top,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Top), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- estimated cld top pressure
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Top_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Pc_Top,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctp_Top), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- aux estimated cld top pressure layer1
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Top1_Aux_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Pc_Top1_Aux,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctp_Top1_Aux), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- aux estimated cld top pressure layer2
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Top2_Aux_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Pc_Top2_Aux,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctp_Top2_Aux), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld base height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Base_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Zc_Base,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Base), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- estimated cld base pressure
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Base_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Pc_Base,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctp_Base), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- lower cld height
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Zc_Lower_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Zc_Lower_Cloud,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Zc_Lower), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- lower cld pressure
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Pc_Lower_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Pc_Lower_Cloud,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Pc_Lower), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld altitude
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Alt_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(ACHA%Alt,sym%LINEAR_SCALING,Min_Alt,Max_Alt,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Alt), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld base altitude
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Base_Alt_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(ACHA%Base_Alt,sym%LINEAR_SCALING,Min_Alt,Max_Alt,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Base_Alt), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- acha processing order
     if (Sds_Num_Level2_Acha_Order_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Order), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        ACHA%Processing_Order(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- acha inversion flag
     if (Sds_Num_Level2_Acha_Inver_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Inver), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        ACHA%Inversion_Flag(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- acha cloud layer
     if (Sds_Num_Level2_Cld_Layer_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Layer), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        ACHA%Cld_Layer(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- acha convective cloud probability
     if (Sds_Num_Level2_Conv_Prob_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Conv_Cld_Prob,sym%LINEAR_SCALING,Min_Prob,Max_Prob,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Conv_Prob), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        One_Byte_Temp(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- acha supercooled cloud probability
     if (Sds_Num_Level2_Supercool_Prob_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Supercooled_Cld_Prob,sym%LINEAR_SCALING,Min_Prob,Max_Prob,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Supercool_Prob), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        One_Byte_Temp(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld height from h2o
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_H2O_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Zc_H2O,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_H2O), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld height from opaque
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Opa_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Zc_Opaque_Cloud,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Opa), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld temperature from opaque
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctt_Opa_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Tc_Opaque_Cloud,sym%LINEAR_SCALING,Min_Tc,Max_Tc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctt_Opa), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld emissivity from split window
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ec_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Ec,sym%LINEAR_SCALING,Min_ec,Max_ec,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_ec), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- beta from split window
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Beta_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Beta,sym%LINEAR_SCALING,Min_beta,Max_beta,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_beta), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud temperature uncertainity from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctt_Acha_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Tc_Uncertainty,sym%LINEAR_SCALING,Min_Tc_Uncer, &
                                 Max_Tc_Uncer,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctt_Acha_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud pressure uncertainity from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Acha_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Pc_Uncertainty,sym%LINEAR_SCALING,Min_Pc_Uncer, &
                                 Max_Pc_Uncer,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctp_Acha_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud pressure uncertainity layer1 from aux
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Aux_Uncer1_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Pc_Uncertainty1_Aux,sym%LINEAR_SCALING,Min_Pc_Uncer, &
                                 Max_Pc_Uncer,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctp_Aux_Uncer1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud pressure uncertainity layer2 from aux
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ctp_Aux_Uncer2_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Pc_Uncertainty2_Aux,sym%LINEAR_SCALING,Min_Pc_Uncer, &
                                 Max_Pc_Uncer,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ctp_Aux_Uncer2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud height uncertainity from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Acha_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Zc_Uncertainty,sym%LINEAR_SCALING,Min_Zc_Uncer, &
                                 Max_Zc_Uncer,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Acha_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud temperature quality flag from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cth_Acha_Qf_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%OE_Quality_Flags(1,:,:)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cth_Acha_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud emissivity quality flag from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ec_Acha_Qf_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%OE_Quality_Flags(2,:,:)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ec_Acha_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud beta quality flag from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Beta_Acha_Qf_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%OE_Quality_Flags(3,:,:)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Beta_Acha_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Acha_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Tau,sym%LINEAR_SCALING,Min_Tau_Acha,Max_Tau_Acha,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Acha), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size from acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Acha_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ACHA%Reff,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Acha), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- quality flag from Acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Acha_Quality_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%Packed_Quality_Flags
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Quality), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- info flag from Acha
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Acha_Info_Flag == sym%YES) then     
      One_Byte_Temp = ACHA%Packed_Meta_Data_Flags
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Info), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Acha Cost
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Acha_Cost_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(ACHA%Cost,sym%LINEAR_SCALING,Min_Acha_Cost,Max_Acha_Cost,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Acha_Cost), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Aux Cost
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Aux_Cost_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Cost_Aux,sym%LINEAR_SCALING,Min_Acha_Cost,Max_Acha_Cost,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Aux_Cost), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Tau_Dcomp,sym%LINEAR_SCALING,Min_tau,Max_tau,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth from aux
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Aux_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Tau_Aux,sym%LINEAR_SCALING,Min_tau,Max_tau,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Aux), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld reff from aux
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Aux_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Reff_Aux,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Aux), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Reff_Dcomp,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth cost from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Dcomp_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Tau_Dcomp_Cost,sym%LINEAR_SCALING,Min_tau,Max_tau,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Dcomp_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size cost from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Dcomp_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Reff_Dcomp_Cost,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Dcomp_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud optical depth quality flag from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Dcomp_Qf_Flag == sym%YES) then     
      One_Byte_Temp = Tau_Dcomp_Qf
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Dcomp_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud effective radius quality flag from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Dcomp_Qf_Flag == sym%YES) then     
      One_Byte_Temp = Reff_Dcomp_Qf
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Dcomp_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- quality flag from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Dcomp_Quality_Flag == sym%YES) then     
      One_Byte_Temp = Dcomp_quality_flag
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Dcomp_Quality), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- info flag  dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Dcomp_Info_Flag == sym%YES) then     
      Two_Byte_Temp = Dcomp_Info_Flag
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Dcomp_Info), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- insolation from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_CldInsol_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Insolation_Dcomp, &
                                 sym%LINEAR_SCALING,Min_Insol,Max_Insol,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_CldInsol), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- diffuse insolation from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_CldInsol_Dif_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Insolation_Dcomp_Diffuse, &
                                 sym%LINEAR_SCALING,Min_Insol,Max_Insol,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_CldInsol_Dif), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cdnc from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cdnc_Dcomp_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Cdnc_Dcomp, &
                                 sym%LINEAR_SCALING,Min_Cdnc,Max_Cdnc,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cdnc_Dcomp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud geo height from  dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Hcld_Dcomp_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Hcld_Dcomp, &
                                 sym%LINEAR_SCALING,Min_Hcld,Max_Hcld,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Hcld_Dcomp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth from nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Nlcomp_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Tau_Nlcomp,sym%LINEAR_SCALING,Min_tau,Max_tau,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Nlcomp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Nlcomp_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Reff_Nlcomp,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Nlcomp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld optical depth cost from nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cod_Nlcomp_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Tau_Nlcomp_Cost,sym%LINEAR_SCALING,Min_Tau,Max_Tau,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cod_Nlcomp_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cld particle size cost from dcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Ceps_Nlcomp_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I2_RANK2(Reff_Nlcomp_Cost,sym%LINEAR_SCALING,Min_Reff,Max_Reff,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ceps_Nlcomp_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- quality flag from nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Nlcomp_Quality_Flag == sym%YES) then     
      One_Byte_Temp = Nlcomp_quality_flag
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Nlcomp_Quality), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- info flag from nlcomp
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Nlcomp_Info_Flag == sym%YES) then     
      Two_Byte_Temp = Nlcomp_Info_Flag
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Nlcomp_Info), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud albedo
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldalb_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(cloud_063um_albedo, &
                                 sym%LINEAR_SCALING,Min_albedo,Max_albedo,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_cldalb), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud transmission
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldtrn_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(cloud_063um_transmission_solar, &
                                 sym%LINEAR_SCALING,Min_transmission,Max_transmission,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_cldtrn), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldfrac_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(Cloud_Fraction, &
                                 sym%LINEAR_SCALING,Min_frac,Max_frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_cldfrac), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- cloud fraction uncertainty
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Cldfrac_Uncer_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(Cloud_Fraction_Uncer, &
                                 sym%LINEAR_SCALING,Min_frac,Max_frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cldfrac_Uncer), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- high cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_High_Cld_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(High_Cloud_Fraction, &
                                 sym%LINEAR_SCALING,Min_Frac,Max_Frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_High_Cld), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- mid cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Mid_Cld_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(Mid_Cloud_Fraction, &
                                 sym%LINEAR_SCALING,Min_Frac,Max_Frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Mid_Cld), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- low cloud fraction
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Low_Cld_Flag == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(Low_Cloud_Fraction, &
                                 sym%LINEAR_SCALING,Min_Frac,Max_Frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Low_Cld), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

!--- non-cloud props     

     !--- emiss ch20
     if (Sds_Num_Level2_Ch20_Emiss_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Ems_Ch20,sym%LINEAR_SCALING,Min_Ems_Ch20,Max_Ems_Ch20,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Emiss), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- emiss ch20 clear
     if (Sds_Num_Level2_Ch20_Emiss_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Ems_Ch20_Clear_Rtm,sym%LINEAR_SCALING,Min_Ems_Ch20,Max_Ems_Ch20,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Emiss_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- emiss ch20 median 3x3
     if (Sds_Num_Level2_Ch20_Emiss_Median_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Ems_Ch20_Median_3x3,sym%LINEAR_SCALING,Min_Ems_Ch20,Max_Ems_Ch20,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Emiss_Median), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- bt ch20 median 3x3
     if (Sds_Num_Level2_Ch20_Bt_Median_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Bt_Ch20_Median_3x3,sym%LINEAR_SCALING,Min_Bt20,Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Bt_Median), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

    !--- Bt_11_67_Covar
    if (Sds_Num_Level2_Bt_11_67_Covar_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
        Sensor%Chan_On_Flag_Default(27) == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Covar_Ch27_Ch31_5x5,sym%LINEAR_SCALING,Min_Bt_Covar,Max_Bt_Covar,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_11_67_Covar), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Btd_Ch31_Ch32_Bt_Ch31_Max_3x3
    if (Sds_Num_Level2_Btd_Ch31_Ch32_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
        Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Btd_Ch31_Ch32_Bt_Ch31_Max_3x3,sym%LINEAR_SCALING,Min_Btd_Ch31_Ch32, &
                       Max_Btd_Ch31_Ch32,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Btd_Ch31_Ch32), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

     !--- etrop
     if (Sds_Num_Level2_Etrop_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then     
      call SCALE_VECTOR_I1_RANK2(ch(31)%Emiss_Tropo,sym%LINEAR_SCALING,Min_Etropo,Max_Etropo,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Etrop), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

    !--- Beta_11_67
    if (Sds_Num_Level2_Beta_11_67_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
        Sensor%Chan_On_Flag_Default(27) == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Beta_11um_67um_Tropo_Rtm,sym%LINEAR_SCALING,Min_Beta,Max_Beta,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Beta_11_67), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Beta_11_85
    if (Sds_Num_Level2_Beta_11_85_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
        Sensor%Chan_On_Flag_Default(29) == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Beta_11um_85um_Tropo_Rtm,sym%LINEAR_SCALING,Min_Beta,Max_Beta,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Beta_11_85), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Beta_11_12
    if (Sds_Num_Level2_Beta_11_12_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
        Sensor%Chan_On_Flag_Default(32) == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Beta_11um_12um_Tropo_Rtm,sym%LINEAR_SCALING,Min_Beta,Max_Beta,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Beta_11_12), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Beta_11_133
    if (Sds_Num_Level2_Beta_11_13_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES .and. &
        Sensor%Chan_On_Flag_Default(33) == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Beta_11um_133um_Tropo_Rtm,sym%LINEAR_SCALING,Min_Beta,Max_Beta,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Beta_11_13), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

     !--- aot1
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot1_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(aot1,sym%LINEAR_SCALING,Min_aot,Max_aot,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_aot1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- aot2
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot2_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(aot2,sym%LINEAR_SCALING,Min_aot,Max_aot,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_aot2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- aot6
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot6_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(aot3a,sym%LINEAR_SCALING,Min_aot,Max_aot,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Aot6), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- aot qf
     if (Aer_Flag == sym%YES .and. Sds_Num_Level2_Aot_Qf_Flag == sym%YES) then
        One_Byte_Temp = aot_qf
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Aot_Qf), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- olr
     if ( Sds_Num_Level2_Olr_Flag == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(olr, &
                                 sym%LINEAR_SCALING,Min_olr,Max_olr,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_olr), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- insolation
     if (Sasrab_Flag == sym%YES .and. Sds_Num_Level2_Insol_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Insolation_All_Sky, &
                                 sym%LINEAR_SCALING,Min_Insol,Max_Insol,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Insol), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- diffuse insolation
     if (Sasrab_Flag == sym%YES .and. Sds_Num_Level2_Insol_Dif_Flag == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Insolation_All_Sky_Diffuse, &
                                 sym%LINEAR_SCALING,Min_Insol,Max_Insol,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Insol_Dif), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- ndvi surface corrected
     if (Sds_Num_Level2_Ndvi_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(ndvi_sfc,sym%LINEAR_SCALING,Min_Ndvi,Max_Ndvi,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ndvi), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- ndvi surface corrected from modis white sky
     if (Sds_Num_Level2_Ndvi_White_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Ndvi_Sfc_White_Sky,sym%LINEAR_SCALING,Min_Ndvi,Max_Ndvi,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ndvi_White), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface temperature
     if (Sds_Num_Level2_Tsfc_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Tsfc_Retrieved,sym%LINEAR_SCALING,Min_Tsfc,Max_Tsfc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Tsfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface radiation temperature
     if (Sds_Num_Level2_Trad_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Trad_Retrieved,sym%LINEAR_SCALING,Min_Tsfc,Max_Tsfc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Trad), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface radiation temperature
     if (Sds_Num_Level2_Tair_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Tair_Nwp_Pix,sym%LINEAR_SCALING,Min_Tsfc,Max_Tsfc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Tair), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface temperature background
     if (Sds_Num_Level2_Tsfc_Back_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Tsfc_Nwp_Pix,sym%LINEAR_SCALING,Min_Tsfc,Max_Tsfc,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Tsfc_Back), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface temperature background uniformity 3x3
     if (Sds_Num_Level2_Tsfc_Uni_Back_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Sst_Anal_Uni,sym%LINEAR_SCALING,Min_Sst_std,Max_Sst_std,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Tsfc_Uni_Back), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface relative humidity background
     if (Sds_Num_Level2_Rh_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Rh_Nwp_Pix,sym%LINEAR_SCALING,Min_Rh,Max_Rh,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Rh), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- surface pressure background
     if (Sds_Num_Level2_Psfc_Back_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Psfc_Nwp_Pix,sym%LINEAR_SCALING,Min_Psfc,Max_Psfc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Psfc_Back), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- sea-level pressure background
     if (Sds_Num_Level2_Pmsl_Back_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Pmsl_Nwp_Pix,sym%LINEAR_SCALING,Min_Pmsl,Max_Pmsl,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Pmsl_Back), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- K index
     if (Sds_Num_Level2_Kindex_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(K_Index_Nwp_Pix,sym%LINEAR_SCALING,Min_Kindex,Max_Kindex,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Kindex), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Cwp from Nwp 
     if (Sds_Num_Level2_Cwp_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I2_RANK2(Cwp_Nwp_Pix,sym%LINEAR_SCALING,Min_Cwp,Max_Cwp,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cwp_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Cfrac from Nwp 
     if (Sds_Num_Level2_Cfrac_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Cfrac_Nwp_Pix,sym%LINEAR_SCALING,Min_Frac,Max_Frac,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cfrac_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Pc from Nwp 
     if (Sds_Num_Level2_Pc_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Pc_Nwp_Pix,sym%LINEAR_SCALING,Min_Pc,Max_Pc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Pc_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Number of Cloud Layers from Nwp 
     if (Sds_Num_Level2_Ncld_Nwp_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ncld_Nwp), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Ncld_Layers_Nwp_Pix(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- Cld Type from Nwp 
     if (Sds_Num_Level2_Cld_Type_Nwp_Flag == sym%YES) then     
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cld_Type_Nwp), Sds_Start_2d,Sds_Stride_2d,Sds_Edge_2d,     &
                        Cld_Type_Nwp_Pix(:,Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- tropopause temperature nwp
     if (Sds_Num_Level2_Temp_Tropo_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Ttropo_Nwp_Pix,sym%LINEAR_SCALING,Min_Ttropo,Max_Ttropo,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Temp_Tropo), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- lifting condensation level nwp
     if (Sds_Num_Level2_LCL_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(LCL_Height_Nwp_Pix,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_LCL_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- convective condensation level nwp
     if (Sds_Num_Level2_CCL_Nwp_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(CCL_Height_Nwp_Pix,sym%LINEAR_SCALING,Min_Zc,Max_Zc,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_CCL_Nwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

     !--- remote sensing reflectance
     if (Sds_Num_Level2_Rsr_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(rsr,sym%LINEAR_SCALING,Min_Rsr,Max_Rsr,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Rsr), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

!--- quality flags

     !--- first level 2 packed quality flags
     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Qf1_Flag == sym%YES) then
      One_Byte_Temp = 0
      One_Byte_Temp = ishft(ACHA%OE_Quality_Flags(1,:,:),6) + ishft(ACHA%OE_Quality_Flags(2,:,:),4) + &
                      ishft(ACHA%OE_Quality_Flags(3,:,:),2) + Tau_Dcomp_Qf
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Qf1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus

     endif

     if (Cld_Flag == sym%YES .and. Sds_Num_Level2_Qf2_Flag == sym%YES) then
        !--- second level 2 packed quality flags
        One_Byte_Temp = 0
        One_Byte_Temp = ishft(Reff_Dcomp_Qf,6) + ishft(Aot_Qf,4) + ishft(Rsr_Qf,2) + cld_mask

        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Qf2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
     endif

    !--- ch1 counts
    if (Sensor%Chan_On_Flag_Default(1) == sym%YES .and. Sds_Num_Level2_Ch1_Counts_Flag == sym%YES) then
        Two_Byte_Temp = Ch1_Counts
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Counts), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch2 counts
    if (Sensor%Chan_On_Flag_Default(2) == sym%YES .and. Sds_Num_Level2_Ch2_Counts_Flag == sym%YES) then
        Two_Byte_Temp = Ch2_Counts
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch2_Counts), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch6 counts
    if (Sensor%Chan_On_Flag_Default(6) == sym%YES .and. Sds_Num_Level2_Ch6_Counts_Flag == sym%YES) then
        Two_Byte_Temp = Ch6_Counts
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch6_Counts), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- tpw
    if (Sds_Num_Level2_Tpw_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Tpw_Nwp_Pix,sym%LINEAR_SCALING,Min_Tpw,Max_Tpw,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Tpw), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ozone
    if (Sds_Num_Level2_Ozone_Flag == sym%YES) then
      call SCALE_VECTOR_I1_RANK2(Ozone_Nwp_Pix,sym%LINEAR_SCALING,Min_Ozone,Max_Ozone,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ozone), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch20 reflectance atmos corrected
    if (Sds_Num_Level2_Ref_Ch20_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(20)%Ref_Sfc,sym%LINEAR_SCALING,Min_Ref_Ch20,Max_Ref_Ch20,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch20_Sfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch1 reflectance atmos corrected
    if (Sds_Num_Level2_Ref_Ch1_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
       call SCALE_VECTOR_I2_RANK2(ch(1)%Ref_Sfc,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch1_Sfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch2 reflectance atmos corrected
    if (Sds_Num_Level2_Ref_Ch2_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(2) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(ch(2)%Ref_Sfc,sym%LINEAR_SCALING,Min_Ref_Ch2,Max_Ref_Ch2,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch2_Sfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch6 reflectance atmos corrected
    if (Sds_Num_Level2_Ref_Ch6_Sfc_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(ch(6)%Ref_Sfc,sym%LINEAR_SCALING,Min_Ref_Ch6,Max_Ref_Ch6,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Ch6_Sfc), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                         One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- sst already masked
    if (Sds_Num_Level2_Sst_Masked_Flag == sym%YES) then
       call SCALE_VECTOR_I1_RANK2(Sst_Masked,sym%LINEAR_SCALING,Min_Sst,Max_Sst,Missing_Value_Real4,One_Byte_Temp)
       Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Sst_Masked), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                      One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif 

    !-- Ch1 reflectance
    if (Sds_Num_Level2_Ch1_Unnorm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(ch(1)%Ref_Toa_Unnorm,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Unnorm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
     
    !-- Ch2 reflectance
    if (Sds_Num_Level2_Ch2_Unnorm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(2) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(ch(2)%Ref_Toa_Unnorm,sym%LINEAR_SCALING,Min_Ref_Ch2,Max_Ref_Ch2,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch2_Unnorm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch6 reflectance
    if (Sds_Num_Level2_Ch6_Unnorm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(6) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(ch(6)%Ref_Toa_Unnorm,sym%LINEAR_SCALING,Min_Ref_Ch6,Max_Ref_Ch6,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch6_Unnorm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Ref_Ch1 clear
    if (Sds_Num_Level2_Ch1_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(1)%Ref_Toa_Clear,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Ch20 temperature clear
    if (Sds_Num_Level2_Ch20_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(ch(20)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt20,Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch27 temperature clear
    if (Sds_Num_Level2_Ch27_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(27) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(27)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt27,Max_Bt27,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch27_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch28 temperature clear
    if (Sds_Num_Level2_Ch28_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(28) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(28)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt28,Max_Bt28,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch28_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch29 temperature clear
    if (Sds_Num_Level2_Ch29_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(29) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(29)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt29,Max_Bt29,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch29_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch30 temperature clear
    if (Sds_Num_Level2_Ch30_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(30) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(30)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt30,Max_Bt30,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch30_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch31 temperature clear
    if (Sds_Num_Level2_Ch31_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(31)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt31,Max_Bt31,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch32 temperature clear
    if (Sds_Num_Level2_Ch32_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(32) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(32)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt32,Max_Bt32,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch32_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch33 temperature clear
    if (Sds_Num_Level2_Ch33_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(33) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(33)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt33,Max_Bt33,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch33_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch34 temperature clear
    if (Sds_Num_Level2_Ch34_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(34) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(34)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt34,Max_Bt34,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch34_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch35 temperature clear
    if (Sds_Num_Level2_Ch35_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(35) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(35)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt35,Max_Bt35,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch35_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch36 temperature clear
    if (Sds_Num_Level2_Ch36_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(36) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(36)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt36,Max_Bt36,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch36_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch37 temperature clear
    if (Sds_Num_Level2_Ch37_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(37) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(37)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt37,Max_Bt37,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch37_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch38 temperature clear
    if (Sds_Num_Level2_Ch38_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(38) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(38)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt38,Max_Bt38,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch38_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus 
    endif

    !--- Ch45 Pseudo temperature clear
    if (Sds_Num_Level2_Ch45_Clear_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(45) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ch(45)%Bt_Toa_Clear,sym%LINEAR_SCALING,Min_Bt45,Max_Bt45,Missing_Value_Real4,Two_Byte_Temp)
     Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch45_Clear), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Ref_Ch1 3x3 Mean
    if (Sds_Num_Level2_Ch1_Mean_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Ref_Ch1_Mean_3x3,sym%LINEAR_SCALING,Min_Ref_Ch1,Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Mean), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- sst unmasked by cloud
    if (Sds_Num_Level2_Sst_Unmasked_Flag == sym%YES) then
     call SCALE_VECTOR_I2_RANK2(Sst_Unmasked,sym%LINEAR_SCALING,Min_Sst,Max_Sst,Missing_Value_Real4,Two_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Sst_Unmasked), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- wind speed
    if (Sds_Num_Level2_Wnd_Spd_Flag == sym%YES) then
     call SCALE_VECTOR_I1_RANK2(Wnd_Spd_10m_Nwp_Pix,sym%LINEAR_SCALING,Min_Wndspd,Max_Wndspd,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Wnd_Spd), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- wind direction
    if (Sds_Num_Level2_Wnd_Dir_Flag == sym%YES) then
     call SCALE_VECTOR_I1_RANK2(Wnd_Dir_10m_Nwp_Pix,sym%LINEAR_SCALING,Min_Wnddir,Max_Wnddir,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Wnd_Dir), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- wind speed at cld top
    if (Sds_Num_Level2_Wnd_Spd_Cld_Top_Flag == sym%YES .and. Cld_Flag == sym%YES) then
     call SCALE_VECTOR_I1_RANK2(Wnd_Spd_Cld_Top_Nwp_Pix,sym%LINEAR_SCALING,Min_Wndspd,Max_Wndspd,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Wnd_Spd_Cld_Top), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- wind direction at cld top
    if (Sds_Num_Level2_Wnd_Dir_Cld_Top_Flag == sym%YES .and. Cld_Flag == sym%YES) then
     call SCALE_VECTOR_I1_RANK2(Wnd_Dir_Cld_Top_Nwp_Pix,sym%LINEAR_SCALING,Min_Wnddir,Max_Wnddir,Missing_Value_Real4,One_Byte_Temp)
      Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Wnd_Dir_Cld_Top), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                       One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Ch1 Dark-Sky Reflectance
    if (Sds_Num_Level2_Ch1_Dark_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(1) == sym%YES) then
        call SCALE_VECTOR_I2_RANK2(Ref_Ch1_Dark_Composite,sym%LINEAR_SCALING,Min_Ref_Ch1, &
                                   Max_Ref_Ch1,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch1_Dark), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                        Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Cloud Water Path
    if (Sds_Num_Level2_Cwp_Flag == sym%YES .and. Cld_Flag == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Cwp_Dcomp,sym%LINEAR_SCALING,Min_Cwp, &
                                   Max_Cwp,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Cwp), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !-- Rain Rate
    if (Sds_Num_Level2_Rain_Rate_Flag == sym%YES .and. Cld_Flag == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Rain_Rate_Dcomp,sym%LINEAR_SCALING,Min_Rain_Rate, &
                                   Max_Rain_Rate,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Rain_Rate), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I1 at M-band
    if (Sds_Num_Level2_Ref_Max_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Max_ChI1,sym%LINEAR_SCALING,Min_Ref_Ch1, &
                                   Max_Ref_Ch1,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Max_ChI1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Min of I1 at M-band
    if (Sds_Num_Level2_Ref_Min_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Min_ChI1,sym%LINEAR_SCALING,Min_Ref_Ch1, &
                                   Max_Ref_Ch1,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Min_ChI1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I1 at M-band
    if (Sds_Num_Level2_Ref_Mean_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Mean_ChI1,sym%LINEAR_SCALING,Min_Ref_Ch1, &
                                   Max_Ref_Ch1,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Mean_ChI1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I1 at M-band
    if (Sds_Num_Level2_Ref_Uni_ChI1_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(39) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Uni_ChI1,sym%LINEAR_SCALING,Min_Uni_Ch1, &
                                   Max_Uni_Ch1,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Uni_ChI1), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I2 at M-band
    if (Sds_Num_Level2_Ref_Max_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Max_ChI2,sym%LINEAR_SCALING,Min_Ref_Ch2, &
                                   Max_Ref_Ch2,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Max_ChI2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Min of I2 at M-band
    if (Sds_Num_Level2_Ref_Min_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Min_ChI2,sym%LINEAR_SCALING,Min_Ref_Ch2, &
                                   Max_Ref_Ch2,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Min_ChI2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I2 at M-band
    if (Sds_Num_Level2_Ref_Mean_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Mean_ChI2,sym%LINEAR_SCALING,Min_Ref_Ch2, &
                                   Max_Ref_Ch2,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Mean_ChI2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I2 at M-band
    if (Sds_Num_Level2_Ref_Uni_ChI2_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(40) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Uni_ChI2,sym%LINEAR_SCALING,Min_Uni_Ch2, &
                                   Max_Uni_Ch2,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Uni_ChI2), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I3 at M-band
    if (Sds_Num_Level2_Ref_Max_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Max_ChI3,sym%LINEAR_SCALING,Min_Ref_Ch6, &
                                   Max_Ref_Ch6,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Max_ChI3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Min of I3 at M-band
    if (Sds_Num_Level2_Ref_Min_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Min_ChI3,sym%LINEAR_SCALING,Min_Ref_Ch6, &
                                   Max_Ref_Ch6,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Min_ChI3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I3 at M-band
    if (Sds_Num_Level2_Ref_Mean_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Mean_ChI3,sym%LINEAR_SCALING,Min_Ref_Ch6, &
                                   Max_Ref_Ch6,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Mean_ChI3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I3 at M-band
    if (Sds_Num_Level2_Ref_Uni_ChI3_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(41) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Ref_Uni_ChI3,sym%LINEAR_SCALING,Min_Uni_Ch6, &
                                   Max_Uni_Ch6,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ref_Uni_ChI3), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I4 at M-band
    if (Sds_Num_Level2_Bt_Max_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Max_ChI4,sym%LINEAR_SCALING,Min_Bt20, &
                                   Max_Bt20,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Max_ChI4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
                                                                                                                                     
    !--- Min of I4 at M-band
    if (Sds_Num_Level2_Bt_Min_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Min_ChI4,sym%LINEAR_SCALING,Min_Bt20, &
                                   Max_Bt20,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Min_ChI4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I4 at M-band
    if (Sds_Num_Level2_Bt_Mean_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Mean_ChI4,sym%LINEAR_SCALING,Min_Bt20, &
                                   Max_Bt20,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Mean_ChI4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I4 at M-band
    if (Sds_Num_Level2_Bt_Uni_ChI4_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(42) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Uni_ChI4,sym%LINEAR_SCALING,Min_Uni_Ch5, &
                                   Max_Uni_Ch5,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Uni_ChI4), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Max of I5 at M-band
    if (Sds_Num_Level2_Bt_Max_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Max_ChI5,sym%LINEAR_SCALING,Min_Bt31, &
                                   Max_Bt31,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Max_ChI5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Min of I5 at M-band
    if (Sds_Num_Level2_Bt_Min_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Min_ChI5,sym%LINEAR_SCALING,Min_Bt31, &
                                   Max_Bt31,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Min_ChI5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Mean of I5 at M-band
    if (Sds_Num_Level2_Bt_Mean_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Mean_ChI5,sym%LINEAR_SCALING,Min_Bt31, &
                                   Max_Bt31,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Mean_ChI5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- Uni of I5 at M-band
    if (Sds_Num_Level2_Bt_Uni_ChI5_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(43) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(Bt_Uni_ChI5,sym%LINEAR_SCALING,Min_Uni_Ch5, &
                                   Max_Uni_Ch5,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt_Uni_ChI5), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif


    !--- ch31 atmospheric radiance
    if (Sds_Num_Level2_Ch31_Rad_Atm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(31)%Rad_Atm,sym%LINEAR_SCALING,Min_Ch31_Rad_Atm, &
                                   Max_Ch31_Rad_Atm,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Rad_Atm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch31 atmospheric transmission
    if (Sds_Num_Level2_Ch31_Trans_Atm_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(31)%Trans_Atm,sym%LINEAR_SCALING,Min_Trans, &
                                   Max_Trans,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Trans_Atm), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch31 downward atmospheric radiance
    if (Sds_Num_Level2_Ch31_Rad_Atm_Dwn_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(31)%Rad_Atm_Dwn_Sfc,sym%LINEAR_SCALING,Min_Ch31_Rad_Atm_Dwn, &
                                   Max_Ch31_Rad_Atm_Dwn,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Rad_Atm_Dwn), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch20 surface emissivity
    if (Sds_Num_Level2_Ch20_Sfc_Emiss_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(20) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(20)%Sfc_Emiss,sym%LINEAR_SCALING,Min_Sfc_Ems, &
                                   Max_Sfc_Ems,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch20_Sfc_Emiss), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- ch31 surface emissivity
    if (Sds_Num_Level2_Ch31_Sfc_Emiss_Flag == sym%YES .and. Sensor%Chan_On_Flag_Default(31) == sym%YES) then
        call SCALE_VECTOR_I1_RANK2(ch(31)%Sfc_Emiss,sym%LINEAR_SCALING,Min_Sfc_Ems, &
                                   Max_Sfc_Ems,Missing_Value_Real4,One_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Ch31_Sfc_Emiss), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          One_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- 3.75 micron BT for Sounder
    if (Sds_Num_Level2_Bt375_Snd_Flag == sym%YES .and. index(Sensor%Sensor_Name,'IFF') > 0) then
        call SCALE_VECTOR_I2_RANK2(Bt_375um_Sounder,sym%LINEAR_SCALING,Min_Bt20, &
                                   Max_Bt20,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt375_Snd), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif
    !--- 11 micron BT for Sounder
    if (Sds_Num_Level2_Bt11_Snd_Flag == sym%YES .and. index(Sensor%Sensor_Name,'IFF') > 0) then
        call SCALE_VECTOR_I2_RANK2(Bt_11um_Sounder,sym%LINEAR_SCALING,Min_Bt31, &
                                   Max_Bt31,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt11_Snd), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

    !--- 12 micron BT for Sounder
    if (Sds_Num_Level2_Bt12_Snd_Flag == sym%YES .and. index(Sensor%Sensor_Name,'IFF') > 0) then
        call SCALE_VECTOR_I2_RANK2(Bt_12um_Sounder,sym%LINEAR_SCALING,Min_Bt32, &
                                   Max_Bt32,Missing_Value_Real4,Two_Byte_Temp)
        Istatus = sfwdata(Sds_Id_Level2(Sds_Num_Level2_Bt12_Snd), Sds_Start_2d, Sds_Stride_2d, Sds_Edge_2d, &
                          Two_Byte_Temp(:, Line_Idx_Min_Segment:Sds_Edge_2d(2) + Line_Idx_Min_Segment - 1)) + Istatus
    endif

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
subroutine CLOSE_PIXEL_HDF_FILES(Rtm_File_Flag,Level2_File_Flag)

 integer, intent(in):: Rtm_File_Flag
 integer, intent(in):: Level2_File_Flag

 integer:: Isds
 integer:: Istatus

! HDF function declarations
 integer:: sfsnatt
 integer:: sfendacc
 integer:: sfend


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
      
      if (sds_Id_Level2(Isds) /= 0 ) Istatus = sfendacc(Sds_Id_Level2(Isds)) + Istatus
      ! - do this if sds_id_level2 is saved from last file and not used here..
      ! - some compilers with ceratin flags save also without save attribute
      ! - 20141225-AW
      sds_Id_Level2(Isds) = 0
     
   enddo
   Istatus = sfend(Sd_Id_Level2) + Istatus
!--- errors are expected on close since all sds were not opened
!  if (Istatus /= 0) then
!     print *, EXE_PROMPT, MOD_PROMPT, "Error closing level2 hdf file"
!  endif
  endif

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

end module LEVEL2_ROUTINES


