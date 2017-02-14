program cx_read_ahi_test

   USE cx_read_ahi_mod, only: &
       ahi_config_type &
     , ahi_data_out_type &
     , GET_AHI_DATA &
     , AHI_SEGMENT_INFORMATION_REGION &
     , get_var_dimension
   
  
   
   character (len =1024) :: file
   type ( ahi_config_type )   :: ahi_c
   type ( ahi_data_out_type ) :: ahi_data, ahi_data2
  
   real :: lon0,lon1
   real :: lat0,lat1
   integer :: offset_all(2), count_all(2)
   integer,pointer::dims(:)=>null()
   
      interface
      subroutine view2d(x,a,b,t)
         real,dimension(:,:) :: x
         real, intent(in), optional::a,b
         character ( len =*) , optional :: t
      end subroutine view2d
   end interface
   
   
         interface
      subroutine aw_image_idl(x,c,d,a,b,t)
         real,dimension(:,:) :: x
         real,dimension(:,:) :: c
         real,dimension(:,:) :: d
         real, intent(in), optional::a,b
         character ( len =*) , optional :: t
      end subroutine aw_image_idl
   end interface
   
      
   lon0=113
   lon1=120
   lat0=20
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
   
   offset_all = [1600,600]
    count_all = [200,200]
   
   
   
   ahi_c % data_path = '/DATA/Satellite_Input/ahi/full_res/'
   ahi_c % file_base = 'HS_H08_20150617_0300_B01_FLDK.nc'
   
   call AHI_SEGMENT_INFORMATION_REGION( ahi_c , offset_all, count_all )
   
   call  ahi_c % clean  
   
   ahi_c % chan_on(:) = .false.
   ahi_c % chan_on(3) = .true.
   ahi_c % chan_on(5) = .true.
   ahi_c % chan_on(1) = .true.
   ahi_c % chan_on(7) = .true.
   ahi_c % h5_offset = offset_all
   ahi_c % h5_count = count_all
   call GET_AHI_DATA( ahi_c, ahi_data ) 
    
  ! call aw_image_idl(ahi_data % chn(5) % ref, ahi_data%geo%lat, &
  !     ahi_data%geo%lon) 
  ! call sleep(4)     
  
 
   call aw_image_idl(ahi_data % chn(1) % ref, ahi_data%geo_124%lat, &
       ahi_data%geo_124%lon)
   call sleep(4)
  ! call aw_image_idl(ahi_data % chn(3) % ref, ahi_data%geo_3%lat, &
  !     ahi_data%geo_3%lon)
  !   call sleep(4)
   
       
    call  ahi_c % clean  
     ahi_c % data_path = '/DATA/Satellite_Input/ahi/NON_FLDK/ahi_seg_nc/'
   ahi_c % file_base = 'HS_H08_20160304_0330_B01_FLDK.nc'

   ahi_c % chan_on(:) = .false.
   ahi_c % chan_on(3) = .true.
   ahi_c % chan_on(5) = .true.
   ahi_c % chan_on(1) = .true.
   ahi_c % chan_on(7) = .true.
   ahi_c % h5_offset = offset_all
   ahi_c % h5_count = count_all
   
   call get_var_dimension ( ahi_c, 5, dims )
   
   call GET_AHI_DATA( ahi_c, ahi_data2 ) 
    
  
   
 !  call aw_image_idl(ahi_data2 % chn(5) % ref, ahi_data2%geo%lat, &
 !      ahi_data2%geo%lon)
   
 !  call   sleep(6)
 !  call aw_image_idl(ahi_data2 % chn(1) % ref, ahi_data2%geo%lat, &
 !      ahi_data2%geo%lon)
  ! call sleep(6)
   call aw_image_idl(ahi_data2 % chn(3) % ref, ahi_data2%geo%lat, &
       ahi_data2%geo%lon)
       
   call  ahi_data2 % deallocate_all
    call  ahi_c % clean 
     
     ahi_c % data_path = '/DATA/Satellite_Input/ahi/NON_FLDK/ahi_fd_nc/'
   ahi_c % file_base = 'HS_H08_20160304_0330_B01_FLDK.nc'

   ahi_c % chan_on(:) = .false.
   ahi_c % chan_on(3) = .true.
   ahi_c % chan_on(5) = .true.
   ahi_c % chan_on(1) = .true.
   ahi_c % chan_on(7) = .true.
   ahi_c % h5_offset = offset_all
   ahi_c % h5_count = count_all
   
   call get_var_dimension ( ahi_c, 5, dims )
   
   call GET_AHI_DATA( ahi_c, ahi_data2 ) 
   
  
    call sleep(6)
 !  call aw_image_idl(ahi_data2 % chn(5) % ref, ahi_data2%geo%lat, &
 !      ahi_data2%geo%lon)
   
 !  call   sleep(6)
 !  call aw_image_idl(ahi_data2 % chn(1) % ref, ahi_data2%geo%lat, &
 !      ahi_data2%geo%lon)
 !  call sleep(6)
   call aw_image_idl(ahi_data2 % chn(3) % ref, ahi_data2%geo%lat, &
       ahi_data2%geo%lon)
   
  
   
   !end do

end program cx_read_ahi_test
