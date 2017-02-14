program cx_sds_test
   use cx_sds_io_mod,only: &
   cx_sds_finfo &
   , cx_sds_read &
   , cx_sds_read_raw
   
      use cx_sds_type_definitions_mod, only: &
      cx_sds_type &
      , cx_att_type &
      , cx_sds_data_type &
      , MAXNCNAM
   
   
   implicit none
   
   character(len=1024) :: file
   include 'cx_sds_constants.F90'

   type (cx_sds_type) :: pp
   integer :: ftype
   integer :: nsds
   integer :: natt                 
   character ( len = MAXNCNAM), allocatable :: sds_name(:)
   character ( len = MAXNCNAM), allocatable :: att_name(:)
   
   integer :: test
   
   type ( cx_sds_type), allocatable, target :: sds(:)
   type ( cx_sds_data_type), pointer :: pd
   real ,allocatable :: cwp(:)
   real ,allocatable :: cwp_2d(:,:)
   real :: wert(1)
  
  
   
   file = '/DATA/Satellite_Output/viirs/npp_d20130821_t0658182_e0659424_b09405.level2.hdf'
   
   test = cx_sds_finfo (file, ftype, nsds, sds_name, natt, att_name)
  
   test = cx_sds_read_raw ( file, 'cloud_water_path', sds)
   
   wert = sds(1) % get_att('RANGE_MAX')
   
  
   if ( cx_sds_read ( file, 'cloud_water_path', cwp_2d) < 0 ) goto 9999
   
   
   print*,cwp_2d
  stop 
    
   
   
   
   ! call sds(1).info
   !pd=>sds(1).data
   !call pd.info
   !print*,'ende'
   !allocate(cwp(3200*768))
   !call pd.transform_to_real(cwp)
    !  print*,cwp(1800:1900)
9999 print*,'end'      
end program cx_sds_test
