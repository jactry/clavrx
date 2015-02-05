! $Id:$
module ahi_clavrx_bridge

    use constants, only: &
      int4 

contains

   subroutine read_ahi_data ( segment_number ,  file_ch01 , error_out )
   
      use cx_read_ahi_mod
   
      integer , intent(in) :: segment_number
      character(len=*), intent(in) :: file_ch01
      integer(kind=int4), intent(out) :: error_out
      
      integer :: modis_chn_list (16)
      type ( ahi_config_type ) :: ahi_c
      type ( ahi_data_out_type ) :: ahi_data
      
      modis_chn_list = [ 3 , 11 , 1 , 2 , 6 , 7 , 20 , 37 ,  27 , 29  &
               , 30 , 20 , 22 , 31 , 32 , 33 ]
      
      print*,'ereache the right bridge!!'
      
      ahi_c % file_base = 'HS_H08_20150125_0230_B01_FLDK.nc'
   ahi_c % chan_on (:) = .false.
   ahi_c % chan_on (5) = .true.
   ahi_c % chan_on (:) = .true.
   ahi_c % data_path = '/DATA/Satellite_Input/ahi/2015/'
   
   ahi_c % h5_offset = [0,0]
   ahi_c % h5_count = [5500,200] 
   
   print*,'start reader'
    call get_ahi_data ( ahi_c, ahi_data )
    print*,'reader done'
     print*, '==> ',ahi_data % chn(1) % is_read
   
   print*, ahi_data % chn(3) % ref (2750,190:199)
   print*, ahi_data % chn(3) % rad (2750,190:199)
   
   print*,ahi_data % geo % lon (2750,190:199)
      stop
   end subroutine read_ahi_data

end module ahi_clavrx_bridge
