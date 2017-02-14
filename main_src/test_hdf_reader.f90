program test
use cx_hdf_read_mod
use cx_data_io_tools_mod

character (len=200) :: file
integer :: nsds,natt
character ( len = MAXNCNAM), allocatable :: sds_name(:) 
character ( len = MAXNCNAM), allocatable :: att_name(:)
integer :: istatus
type(hdf_att), allocatable :: attrs(:)
integer :: i

character(len =MAXNCNAM), allocatable :: att_names(:)
   integer :: N_ATTS 
   
   
   
   



file='clavrx_npp_d20130101_t1210018_e1211259_b06116.level2.hdf'

!istatus =  hdf_get_finfo (file,nsds,sds_name,natt,att_name)
!print*,att_name


!print*
istatus = hdf_get_file_att (file,natt,attrs)

do i=1,natt
print*
   print*,i,trim(attrs(i)%name),attrs(i) % data % i2values,attrs(i) % data % i4values
enddo


N_ATTS = 6
      allocate ( att_names ( N_ATTS))
      att_names = (/ character(len=50) ::  'START_YEAR','START_DAY_OF_YEAR','DCOMP_VERSION','DCOMP_MODE','END_YEAR' &
                                          , 'END_DAY_OF_YEAR','END_TIME','WMO_SATELLITE_CODE' /)

!if ( hdf_get_file_att (file,natt,attrs,nattn = N_ATTS, att_name = att_names) < 0 ) then
!         print*,'error reading attributes'
!         stop
!      end if
!      print*,size(attrs)
!      print*,attrs(1) % data 


call copy_global_attributes ( file, 'test.hdf', exclude_list=['hallo','DCOMP_VERSION'])


end program test
