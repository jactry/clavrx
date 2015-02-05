! $Id:$
module ahi_clavrx_bridge
   
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
   
    use constants, only: &
      int4 

contains

   subroutine read_ahi_data ( segment_number ,  file_ch01 , error_out )
      use planck
      use cx_read_ahi_mod
      implicit none
      integer , intent(in) :: segment_number
      character(len=*), intent(in) :: file_ch01
      integer(kind=int4), intent(out) :: error_out
      
      integer :: modis_chn_list (16)
      type ( ahi_config_type ) :: ahi_c
      type ( ahi_data_out_type ) :: ahi_data
      integer :: y_start, c_seg_lines
      integer, parameter :: NUM_CHN_AHI = 16
      integer , parameter :: SYM_YES = 1, SYM_NO =0
      integer :: i_chn
      logical :: is_solar_channel(16)
      integer :: modis_chn
      integer :: i
     ! real, parameter :: MISSING_VALUE_REAL4 = -999.
       print*,'entering ahi read routines!!  WELCOME'
      print*,'AHI-TODO=> clean ahi_clavrx_bridge.f90 code'
      modis_chn_list = [ 3 , 4 , 1 , 2 , 6 , 7 , 20 , 43 ,  27 , 28  &
               , 29 , 30 , 44 , 31 , 32 , 33 ]
      
      is_solar_channel(7:16) = .false.
      is_solar_channel(1:6) = .true.
      print*,'reache the AHI bridge!!'
      
      ahi_c % file_base = file_ch01
   ahi_c % chan_on (:) = Sensor%Chan_On_Flag_Default ( modis_chn_list) == 1
    
   
   ahi_c % data_path = trim(Image%Level1b_Path)
   
   y_start = ( segment_number -1 ) * Image%Number_Of_Lines_Per_Segment
   c_seg_lines = min (  y_start + Image%Number_Of_Lines_Per_Segment -1 &
      , Image%Number_Of_Lines )  - y_start  + 1
   ahi_c % h5_offset = [0,y_start]
   ahi_c % h5_count = [Image%Number_Of_Elements  , c_seg_lines] 
   
   print*,'start reader'
   call get_ahi_data ( ahi_c, ahi_data )
   print*,'reader done'
   
   nav % lat_1b(:,1:c_seg_lines)    = ahi_data % geo % lat
   nav % lon_1b(:,1:c_seg_lines)    = ahi_data % geo % lon
   Print*, 'AHI:  =====> TODO Add scan_time in reader and bridge'
   !scan_time(1:c_seg_lines)   =ahi_data % geo % scan_time
   geo % sataz(:,1:c_seg_lines)     = ahi_data % geo % sataz
   geo % satzen(:,1:c_seg_lines)    = ahi_data % geo % satzen
   geo % solaz (:,1:c_seg_lines)    = ahi_data % geo % solaz 
   geo % solzen (:,1:c_seg_lines)   = ahi_data % geo % solzen 
   geo % relaz (:,1:c_seg_lines)   = ahi_data % geo % relaz
   geo % glintzen (:,1:c_seg_lines)   = ahi_data % geo % glintzen
    geo % scatangle (:,1:c_seg_lines)   = ahi_data % geo % scatangle
   print*,'AHI  =====> TODO add ascend, ( do we need this? What should it be? '
   !nav % ascend (1:c_seg_lines)     = ahi_data % geo % ascend
   
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
         
         if ( modis_chn > 36 ) then
            print*,'channel AHI, MODIS: ',i_chn,modis_chn
            print*,'AHI TODO: New Channels 43 and 44 not working for Planck and anything!!'
            cycle
         end if
        ch(modis_chn) % Rad_Toa ( : ,1:c_seg_lines)  =  ahi_data % chn (i_chn) % rad
       call compute_bt_array ( ch(modis_chn)%bt_toa ( : ,1:c_seg_lines) , ch(modis_chn)%rad_toa ( : ,1:c_seg_lines) &
                , modis_chn ,MISSING_VALUE_REAL4 )
        print*,'BT ', modis_chn, ch(modis_chn)%bt_toa (2750,190:199)
        
      end if   
      
   end do
    print*,'AHI:  =====> TODO: Add bt computations in ahi_clavrx_bridge.f90 !! Check if sensor is correctly set etc..'
   
   Image%Number_Of_Lines_Read_This_Segment = c_seg_lines
    do i = 1, Image%Number_Of_Lines_Per_Segment
         scan_number(i) = y_start + i 
      end do
      
    nav % ascend = 0 
    Cloud_Mask_Aux_Read_Flag = 0 
     call ahi_data % deallocate_all
  ! print*, '==> ',ahi_data % chn(1) % is_read
  ! print*, 'example REFLECTANCE [0-100] channel 7 (3.9um):',ahi_data % chn(7) % ref (2750,190:199)
  ! print*
  ! print*, 'example REFLECTANCE [0-100] channel 3:',ahi_data % chn(3) % ref (2750,190:199)
  ! print*
  ! print*, 'example RADIANCE:',ahi_data % chn(3) % rad (2750,190:199)
   
  ! print*,'example LON:',ahi_data % geo % lon (2750,190:199)
 
  !    stop
  
  
  print*,'leaving ahi read routines!!'
   end subroutine read_ahi_data

end module ahi_clavrx_bridge
