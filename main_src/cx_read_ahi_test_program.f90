program cx_read_ahi_test

   USE cx_read_ahi_mod, only: &
       ahi_config_type &
     , ahi_data_out_type &
     , GET_AHI_DATA &
     , AHI_SEGMENT_INFORMATION_REGION
   
   character (len =1024) :: file
   type ( ahi_config_type )   :: ahi_c
   type ( ahi_data_out_type ) :: ahi_data, ahi_data2
   integer :: i   
   real :: lon0,lon1
   real :: lat0,lat1
   integer :: offset_all(2), count_all(2)
   
   
      interface
      subroutine view2d(x,a,b,t)
         real,dimension(:,:) :: x
         real, intent(in), optional::a,b
         character ( len =*) , optional :: t
      end subroutine view2d
   end interface
   
      
   lon0=110
   lon1=130
   lat0=-30
   lat1=30   
      
   file='HS_H08_20150103_0300_B01_FLDK.nc'
   
   ahi_c % data_path = '/DATA/Satellite_Input/ahi/2015/'
   ahi_c % file_base = file
   ahi_c % chan_on(:) = .true.
   
   ahi_c % lon_range =[lon0,lon1]
   ahi_c % lat_range =[Lat0,lat1]
   
  !call AHI_SEGMENT_INFORMATION_REGION( ahi_c , offset_all, count_all )
   
  ! print*,offset_all,count_all
   
   ahi_c % h5_offset = [404,2004]
   ahi_c % h5_count = [100,200] 
   
   
  ! call GET_AHI_DATA( ahi_c, ahi_data ,only_nav=.true.)
   !call ahi_data % time_start_obj % print_data
   
   !do i=1,16
   
   !print*, ahi_data % chn(i) % is_read
   !print*, i, ahi_data % chn(1) % is_solar_channel
  ! print*, shape(ahi_data % chn(1) % ref)
   
   !print*,ahi_data % chn(1) % ref
   
   
   
   !AW okay works fine now with high resolution files
   
   offset_all = [1000,1000]
    count_all = [200,300]
   
   
   
   ahi_c % data_path = '/DATA/Satellite_Input/ahi/full_res/'
   ahi_c % file_base = 'HS_H08_20150617_0300_B01_FLDK.nc'
   ahi_c % chan_on(:) = .false.
   ahi_c % chan_on(5) = .true.
   ahi_c % h5_offset = offset_all
   ahi_c % h5_count = count_all
   call GET_AHI_DATA( ahi_c, ahi_data ) 
   call view2d(ahi_data % chn(5) % ref)
   ahi_c % chan_on(:) = .false.   
   ahi_c % chan_on(3) = .true.
   ahi_c % h5_offset = 4 * offset_all
   ahi_c % h5_count = 4 * count_all
   call GET_AHI_DATA( ahi_c, ahi_data2 )
   call view2d(ahi_data2 % chn(3) % ref)
   !print*, shape(ahi_data2 % chn(3) % ref)
  ! print*, shape(ahi_data % chn(5) % ref) 
   
   print*,'ssss'
   
   
   
   !end do

end program cx_read_ahi_test
