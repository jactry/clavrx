! $Id$
module AHI_CLAVRX_BRIDGE 
   
   use Pixel_Common , only : &
        Image &
      , Sensor &
      , Geo &
      , Nav &
      , Scan_Number &
      , Ancil_Data_Dir & 
      , Cloud_Mask_Aux_Flag &
      , Cloud_Mask_Aux_Read_Flag &
      , Cld_Mask_Aux &
      , Cld_Type_Aux &
      , Cld_Phase_Aux &
      , Scan_Time &
      , Gap_Pixel_Mask &
      , Ch 
   
   use CX_Constants_mod, only: &
      Int4 &
      ,Missing_Value_Real4
      
      
   use Calibration_Constants, only: &
            Planck_Nu
   
   use Planck, only: CONVERT_RADIANCE
   
   implicit none
   private
   public :: read_ahi_data
    
contains

   subroutine READ_AHI_DATA ( Segment_Number , File_Ch01 , Error_Out )
      use planck
      use cx_read_ahi_mod
      
      implicit none
      
      integer , intent(in) :: segment_number
      character(len=*), intent(in) :: file_ch01
      integer(kind=int4), intent(out) :: error_out
      
      integer :: modis_chn_list (16)
      real :: nu_list (7:16)
      type ( ahi_config_type )   :: ahi_c
      type ( ahi_data_out_type ) :: ahi_data
      integer :: y_start
      integer :: c_seg_lines
      integer, parameter :: NUM_CHN_AHI = 16
      integer, parameter :: SYM_YES = 1, SYM_NO =0
      integer :: i_chn
      logical :: is_solar_channel(16)
      
      integer :: modis_chn
      integer :: i_line
      integer :: offset_all(2), count_all(2)
     
      error_out = 0
      modis_chn_list = [ 3 , 4 , 1 , 2 , 6 , 7 , 20 , 37 ,  27 , 28,  &
                        29 , 30 , 38 , 31 , 32 , 33 ]
      
      ! List of solar constant values for radiance transformation
      nu_list(7:16) = [ Planck_Nu(20), Planck_Nu(37), Planck_Nu(27),  &
                        Planck_Nu(28), Planck_Nu(29), Planck_Nu(30),  &
                        Planck_Nu(38), Planck_Nu(31), Planck_Nu(32),  &
                        Planck_Nu(33)]
      
      is_solar_channel(7:16) = .false.
      is_solar_channel(1:6) = .true.
    
      ahi_c % file_base = file_ch01
      ahi_c % chan_on (:) = Sensor%Chan_On_Flag_Default ( modis_chn_list) == SYM_YES
    
      ahi_c % data_path = trim(Image%Level1b_Path)
      
      offset_all = [0,0]
      
      if ( nav % lon_lat_limits_set ) then
      
         ahi_c % lon_range =[Nav%Lon_Min_Limit,Nav%Lon_Max_Limit]
         ahi_c % lat_range =[Nav%Lat_Min_Limit,Nav%Lat_Max_Limit]
   
         call AHI_SEGMENT_INFORMATION_REGION( ahi_c , offset_all, count_all )
      end if
      
      y_start = offset_all(2) + ( segment_number -1 ) * Image%Number_Of_Lines_Per_Segment
      c_seg_lines = Image%Number_Of_Lines_Per_Segment
     
      if ( (c_seg_lines + y_start) > (Image%Number_Of_Lines+offset_all(2)) ) then
         c_seg_lines = Image%Number_Of_Lines - y_start + offset_all(2)  
      end if
 
      ahi_c % h5_offset = [offset_all(1),y_start]
      ahi_c % h5_count = [Image%Number_Of_Elements  , c_seg_lines] 
      
       
      call GET_AHI_DATA( ahi_c, ahi_data )
      
      if ( .not. ahi_data % success ) then
         print*,'AHI Read had errors '
         print*,'Please check messages above'
         print*,'stopping -- '
         print*,' In a later version of CLAVR-x the next file should be start now...'
         stop
      end if
     
      nav % lat_1b(:,1:c_seg_lines)    = ahi_data % geo % lat
      nav % lon_1b(:,1:c_seg_lines)    = ahi_data % geo % lon
      
      geo % sataz(:,1:c_seg_lines)        = ahi_data % geo % sataz
      geo % satzen(:,1:c_seg_lines)       = ahi_data % geo % satzen
      geo % solaz (:,1:c_seg_lines)       = ahi_data % geo % solaz 
      geo % solzen (:,1:c_seg_lines)      = ahi_data % geo % solzen 
      geo % relaz (:,1:c_seg_lines)       = ahi_data % geo % relaz
      geo % glintzen (:,1:c_seg_lines)    = ahi_data % geo % glintzen
      geo % scatangle (:,1:c_seg_lines)   = ahi_data % geo % scatangle
      ! nav % ascend (1:c_seg_lines)     = ahi_data % geo % ascend
   
      do i_chn = 1, NUM_CHN_AHI
         
         
         modis_chn = modis_chn_list (i_chn)
        
         if ( .not. ahi_data % chn ( i_chn ) % is_read ) then
            sensor % chan_on_flag_per_line (modis_chn ,1:c_seg_lines) = SYM_NO 
            cycle   
         end if
      
         if ( .not. ahi_c % chan_on(i_chn) ) cycle
      
         if ( is_solar_channel ( i_chn) ) then
           
            ch(modis_chn) % Ref_Toa ( : ,1:c_seg_lines)  =  ahi_data % chn (i_chn) % ref
         else
            if ( modis_chn > 38) then
               cycle
            end if
            
            ch(modis_chn) % Rad_Toa ( : ,1:c_seg_lines)  =  ahi_data % chn (i_chn) % rad
            call CONVERT_RADIANCE ( ch(modis_chn) % Rad_Toa ( : ,1:c_seg_lines) , nu_list(i_chn), -999. )
            call COMPUTE_BT_ARRAY ( ch(modis_chn)%bt_toa ( : ,1:c_seg_lines) &
               , ch(modis_chn)%rad_toa ( : ,1:c_seg_lines) &
                , modis_chn ,MISSING_VALUE_REAL4 )
 
         end if   
      
      end do
    
  
      Image%Number_Of_Lines_Read_This_Segment = c_seg_lines
      scan_number = [(i_line , i_line = y_start+ 1 , y_start+ Image%Number_Of_Lines_Per_Segment  , 1)]

      nav % ascend = 0 
      Cloud_Mask_Aux_Read_Flag = 0 
      
     
      
      ! - update time
      call AHI_DATA % TIME_START_OBJ % GET_DATE ( msec_of_day = Image%Start_Time  )
      call AHI_DATA % TIME_END_OBJ % GET_DATE ( msec_of_day = Image%End_Time  )     
         
      scan_time(1:c_seg_lines)   = Image%Start_Time + &
                                 ( scan_number(1:c_seg_lines) * (Image%End_Time - Image%Start_Time)) &
                                 / Image%Number_Of_Lines
      
      call AHI_DATA % DEALLOCATE_ALL 
 
 
   end subroutine READ_AHI_DATA 

end module AHI_CLAVRX_BRIDGE
