!  $Id:$
!   sfc data object class
!
!   HISTORY:
!     01/21/2015: AW created 
!
module cx_sfc_mod
   use sfc_tools, only : &
      read_hdf_global_attribute_float64 &
    , read_hdf_sds_dimensions &
    , read_hdf_sds
   
   implicit none
   
   private
   
   integer, parameter, public:: int1 = selected_int_kind(1)
   integer, parameter, public:: int2 = selected_int_kind(3)
   integer, parameter, public:: int4 = selected_int_kind(8)
   integer, parameter, public:: int8 = selected_int_kind(10)
   integer, parameter, public:: real4 = selected_real_kind(6,37)
   integer, parameter, public:: real8 = selected_real_kind(15,307)
   integer, parameter, public:: ipre = real4
    
   type, public :: sfc_file_metadata_type
      logical :: is_opened = .false.
      character ( len = 255 ) :: filename
      character(len=256) :: sds_name
      integer (kind=int4) :: num_lat
      integer (kind=int4) :: num_lon
      real (kind=real8) :: del_lat
      real (kind=real8) :: del_lon
      real (kind=real8) :: first_lat
      real (kind=real8) :: first_lon
      integer :: file_id
   end type sfc_file_metadata_type 
     
   type :: sfc_data_i1_type
      logical :: is_set
      type (sfc_file_metadata_type ) :: meta
      integer (kind=int1) , dimension(:,:), allocatable :: data
   end type sfc_data_i1_type
   
    type :: sfc_data_i2_type
      logical :: is_set
      type (sfc_file_metadata_type ) :: meta
      integer (int2) , dimension(:,:), allocatable :: data
   end type sfc_data_i2_type
   
   type sfc_data_r4_type
      logical :: do_use
      character (len=255) :: file 
      character ( len = 255) :: sds_name
      real, dimension(:,:) , allocatable :: data 
      real, dimension(:,:) , allocatable :: stdv
      type (sfc_file_metadata_type ) :: meta
   end type sfc_data_r4_type
   
      ! --- main sfc type
   type , public :: sfc_main_type
      integer :: dim_x
      integer :: dim_y 
      type ( sfc_data_i1_type ) :: land_class
      type ( sfc_data_i1_type ) :: sfc_type
      type ( sfc_data_i1_type ) :: coast_mask
      type ( sfc_data_i1_type ) :: snow_class
      type ( sfc_data_i1_type ) :: volcano
      type ( sfc_data_r4_type ) :: ndvi
      
      type ( sfc_data_r4_type ) :: sst
      type ( sfc_data_r4_type ) :: sea_ice_fraction
      type ( sfc_data_i2_type ) :: elevation
      type ( sfc_data_r4_type ) , dimension(7) :: modis_w_sky
      type ( sfc_data_r4_type ) , dimension(7) :: albedo_snow
      type ( sfc_data_r4_type ) , dimension(20:36) :: emis
      type ( sfc_data_r4_type )  :: z
      contains
         procedure :: populate
         procedure :: deallocate_all => deallocate_sfc
         procedure :: update_with_nwp
   end type sfc_main_type
      !        
   type, public :: sfc_config_type
      logical  :: chn_on ( 42)
   end type sfc_config_type
   

     
   type et_land_class_type
      integer :: FIRST = 0
      integer :: SHALLOW_OCEAN = 0
      integer :: LAND = 1
      integer :: COASTLINE = 2
      integer :: SHALLOW_INLAND_WATER = 3
      integer :: EPHEMERAL_WATER = 4
      integer :: DEEP_INLAND_WATER = 5
      integer :: MODERATE_OCEAN = 6
      integer :: DEEP_OCEAN = 7 
      integer :: LAST = 7     
   end type
   
   type ( et_land_class_type) , public :: ET_land_class 
   
   
   type et_snow_class_type
      integer :: FIRST = 1
      integer :: NO_SNOW = 1
      integer :: NO_SNOW_NOR_ICE = 1
      integer :: SEA_ICE = 2
      integer :: SNOW = 3
      integer :: LAST = 3   
   end type
   type ( et_snow_class_type) , public :: ET_snow_class   
contains  

   ! ================================================================
   !
   ! ================================================================
   subroutine populate ( this, date , lat , lon , conf,  nwp   )
      
      use date_tools_mod ,only: &
          date_type, period_16 &
          , first_doy_month
      use file_tools    
      use cx_nwp_mod , only: &
          nwp_main_type 
      use cr_config_mod
      use cr_oisst_mod
      
      ! use sfc_default_values_mod
      
      implicit none
      
      class ( sfc_main_type) :: this 
      type ( date_type ) , intent (in) :: date
      real, dimension(:,:) :: lat , lon
      type ( nwp_main_type ), intent(in) :: nwp
      type ( conf_user_opt_type ) , intent (in) :: conf
      
      character ( len =30), dimension(7) ::  w_sky_string
      character ( len =3) :: day_string
      character ( len =100) :: ancil_path
      real , parameter :: WHITE_SKY_EMPIRIC_FACTOR = 1.1
      character(len=4) :: year_str_emiss_files
      integer , dimension(20:36) :: chn_seebor
      integer :: dim_all(2)
      
      character ( len =200) ::snow_glob_filename
      character ( len =200) ::snow_hires_filename
      character ( len =200) ::oisst_filename
      
      integer :: i , nx , ny
      
      ! - executable
      
      
      dim_all = shape(lat)
      this % dim_x = dim_all(1)
      this % dim_y = dim_all(2)
    
      ! 1. land class     
      this % land_class % meta % filename = &
           &  trim(conf % ancil_path) // '/static/sfc_data/lw_geo_2001001_v03m.hdf'
      this % land_class % meta % sds_name = 'land_sea_mask'
      call open_file (this % land_class % meta )
      call read_meta (this % land_class % meta )
      call read_data_i1 ( this % land_class , lat , lon)
      call close_file ( this % land_class % meta )
    
      ! 2. sfc type  
      this % sfc_type % meta % filename = &
         & trim(conf % ancil_path) // '/static/sfc_data/gl-latlong-1km-landcover.hdf'
      this % sfc_type % meta % sds_name = 'surface_type'  
      call open_file (this % sfc_type % meta )
      call read_meta (this % sfc_type % meta )
      call read_data_i1 ( this % sfc_type , lat , lon)
      call close_file ( this % sfc_type % meta )
      
      ! 3. coastal mask
      this % coast_mask % meta % filename = &
         &  trim(conf % ancil_path) //'/static/sfc_data/coast_mask_1km.hdf'
      this % coast_mask % meta % sds_name = 'coast_mask'  
      call open_file (this % coast_mask % meta )
      call read_meta (this % coast_mask % meta )
      call read_data_i1 ( this % coast_mask , lat , lon)
      call close_file ( this % coast_mask % meta )
      
      ! 4. elevation
      this % elevation % meta % filename = &
         &  trim(conf % ancil_path) //'/static/sfc_data/GLOBE_1km_digelev.hdf'
      this % elevation % meta % sds_name = 'surface_elevation'  
      call open_file (this % elevation % meta )
      call read_meta (this % elevation % meta )    
      call read_data_i2 ( this % elevation , lat , lon) 
      call close_file ( this % elevation % meta )     
      this % elevation % is_set = .true.
      
      !5.volcanoe
      this % volcano % meta % filename = &
         &  trim(conf % ancil_path) //'/static/sfc_data/volcano_mask_1km.hdf'
      this % volcano % meta % sds_name = 'volcano_mask'  
      call open_file (this % volcano % meta )
      call read_meta (this % volcano % meta )  
      call read_data_i1 ( this % volcano , lat , lon) 
      call close_file ( this % volcano % meta )      
      this % volcano % is_set = .true.

      ! 6. snow class
      ! we use nwp data for snow
     
      !-TODO  differentiate snow options
      !- now only globSnow
      
      
      snow_hires_filename =  trim(conf % ancil_path) //'/dynamic/snow/hires/snow_map_4km_'//date % yymmdd//'.hdf' 
      snow_glob_filename  =  trim(conf % ancil_path) //'/dynamic/snow/globsnow/'//date%yyyy//'/GlobSnow_SWE_L3A_'//date%yyyymmdd//'_v1.0.hdf'
         
      if ( file_test (snow_hires_filename)) then  
         this % snow_class % meta % filename =  snow_hires_filename
         this % snow_class % meta % sds_name = 'snow_ice_cover'  
         call open_file (this % snow_class % meta )
         call read_meta (this % snow_class % meta )
         call read_data_i1 ( this % snow_class , lat , lon)
         this % snow_class % is_set = .true.
      else if (file_test (snow_glob_filename)) then
         this % snow_class % meta % filename = snow_glob_filename 
         
         print*,'glob snow data will be re-implemented soon '
         !this % snow_class = 
         this % snow_class % is_set = .false.
      else 
         
         this % snow_class % is_set = .false.
      end if
     
                
      ! 7. white sky albedo
      day_string =   period_16 ( date ) 
     
      ! - part of modis w sky name
      w_sky_string = 'none   '
      w_sky_string(1) = '0.659'
      w_sky_string(2) = '0.858' 
      w_sky_string(5) = '1.24' 
      w_sky_string(6) = '1.64'
      w_sky_string(7) = '2.13' 
     
      ancil_path = trim(conf % ancil_path) // '/static/sfc_data/'                                   
      do i = 1, 7 
         
         if ( trim(w_sky_string(i)) == 'none') cycle
         this %  modis_w_sky(i) % meta % filename =  trim(ancil_path)//'AlbMap.WS.c004.v2.0.00-04.' &
                   //trim(day_string)//'.'//trim(w_sky_string( i ))//'_x4.hdf'
         this %  modis_w_sky(i) % meta % sds_name = 'Albedo_Map_'//trim(w_sky_string( i ))
        
         call open_file (this % modis_w_sky(i) % meta )
         call read_meta (this % modis_w_sky(i) % meta )         
         call read_modis_white_sky ( this % modis_w_sky(i) , lat , lon)
         call close_file ( this % modis_w_sky(i)  % meta ) 
         
          !  Empirical correction of white sky albedo! 
          !  Personal communication: Jan Fokke Meiring 
          !  factor is  1.1 (AW 2013/05/29 )
         where (  this % modis_w_sky(i) % data /= -999.)
            this % modis_w_sky(i) % data &
                  & = WHITE_SKY_EMPIRIC_FACTOR * this % modis_w_sky(i) % data
         end where
          
         where ( this % land_class % data /= 1)
            this % modis_w_sky(i) % data = 5.
         end where
          
      end do
     
      
      ! 8. - seebor emis
      year_str_emiss_files = '2005'
     
      this %  emis (:) % meta % filename =  trim(ancil_path)//'global_emiss_intABI_' &
            & //trim(year_str_emiss_files)//trim(first_doy_month(date))//'.hdf'
      
      call open_file (this % emis(20) % meta )
      this % emis(21:) % meta %  file_id =  this % emis(20) % meta %  file_id 
      this % emis(21:) % meta %  is_opened = .true. 
      ! emis meta data hard coded
      this % emis(:) % meta % del_lat = 0.05
      this % emis(:) % meta % del_lon = 0.05
      this % emis(:) % meta % first_lat = 89.9750
      this % emis(:) % meta % first_lon = -179.975
      
      ! map modis channels to seebor numbers
      chn_seebor ( 20:36) = [ 7 , 7 , 7 , 7 , 7 , 7 , 0 , 10 , 10 , 11 &
                           , 0 , 14 , 15 , 16 , 16 , 16 , 16 ]
      
      do i = 20 , 36 
           
         if ( chn_seebor(i) == 0 ) cycle
         write(this %  emis(i) % meta % sds_name ,'(a,i0)') 'emiss', chn_seebor ( i )
         call read_emis ( this % emis(i) , lat , lon)
         
         where ( this % emis (i) % data < 0 )
            this % emis ( i) % data = 0.99 
         end where
         
         where ( lat == -999. )
            this % emis ( i) % data = -999. 
         end where
      
      end do  
      
      call close_file (this % emis(20) % meta )
      
      
      !9. SST from oisst
      oisst_filename = get_oisst_map_filename ( date,trim(conf % ancil_path) //'/dynamic/oisst/')
      this % sst % meta % filename = oisst_filename
      
      if (.not. oisst_filename == 'no_file') then
         nx = size (lat,1)
         ny = size(lat,2)
         if (.not. allocated (this % sst % data) ) allocate ( this % sst % data (nx,ny))
         if (.not. allocated (this % sst % stdv) ) allocate ( this % sst % stdv (nx,ny))
         if (.not. allocated (this % sea_ice_fraction % data) ) allocate ( this % sea_ice_fraction % data (nx,ny))
         
         call get_oisst_data (this % sst % meta % filename , conf % temp_path , lat ,lon , this % sst % data , this % sea_ice_fraction % data , this % sst % stdv)
         
      end if
  

   end subroutine
   
   ! ==========================================================
   !
   ! ==========================================================
   
   subroutine open_file ( file )
      implicit none
      type ( sfc_file_metadata_type ) :: file 
      logical :: file_exists
      integer :: fail
      integer :: dfacc_read
      integer :: sfstart
            
      DFACC_READ = 1
      FAIL = -1
       
      inquire(file = trim(file % filename ) , exist = file_exists)
     
      if (.not. file_exists) then
         print "(/,a,'Surface file, ',a,' does not exist.')",'ff', trim( file % filename )
         stop
      end if
     
      
      file % file_id = sfstart(trim(file % filename  ), DFACC_READ)
     
   end subroutine open_file

   ! ============================================================
   !
   ! ============================================================   
   subroutine close_file ( file )
      implicit none
      type ( sfc_file_metadata_type ) :: file 
      integer :: sfend
      integer :: istatus
      
      istatus = sfend ( file % file_id )
      
   
   end subroutine close_file
   
   ! ============================================================
   !
   ! ============================================================
   subroutine read_meta ( file )
      type (  sfc_file_metadata_type ) :: file
          
      real(kind=real8), parameter :: first_lat_default = -90.0_real8
      real(kind=real8), parameter :: first_lon_default = -180.0_real8
      real(kind=real8), parameter :: del_lat_default = 0.04_real8
      real(kind=real8), parameter :: del_lon_default = 0.04_real8
      
      real , parameter :: missing_value_real8 = -999.
       
      file%del_lat = read_hdf_global_attribute_float64(file%file_id, "dlat")
      file%del_lon = read_hdf_global_attribute_float64(file%file_id, "dlon")
      file%first_lat = read_hdf_global_attribute_float64(file%file_id, "first_lat")
      file%first_lon = read_hdf_global_attribute_float64(file%file_id, "first_lon")
   
      if (file%del_lat == missing_value_real8) file%del_lat = del_lat_default
      if (file%del_lon == missing_value_real8) file%del_lon = del_lon_default
      if (file%first_lat == missing_value_real8) file%first_lat = first_lat_default
      if (file%first_lon == missing_value_real8) file%first_lon = first_lon_default
  
      call read_hdf_sds_dimensions(file%file_id, file%sds_name, file%num_lat, file%num_lon) 
      
   end subroutine read_meta
   
   ! ===============================================================
   !    returns array with index in to read in buffer region
   !    buffer region is determined by lon/kat of satellite image
   ! ===============================================================
   
   subroutine index_in_buffer_array ( data_obj_meta, lat, lon, ilat_buffer , ilon_buffer,start_2d,stride_2d,edge_2d )
      type ( sfc_file_metadata_type ) :: data_obj_meta 
      real, dimension(:,:) :: lat , lon
      integer, dimension(2) :: start_2d, stride_2d, edge_2d
      
      real :: west_lon
      real :: east_lon
      real :: south_lat
      real :: north_lat
      integer :: nr_lat
      integer :: nr_lon
      integer :: idx_lon_first
      integer :: idx_lon_last
      integer :: idx_lat_first
      integer :: idx_lat_last
      integer  :: nx, ny
      
      integer ( int4 ) , allocatable , dimension(:,:) :: ilat , ilon
      integer ( int4 ) , allocatable , dimension(:,:) :: ilat_buffer , ilon_buffer
      integer :: temp
      
      north_lat = maxval ( lat , mask = lat >= -90.0 .and. lat <= 90.0 )
      south_lat = minval ( lat , mask = lat >= -90. .and. lat <= 90.0 )
      east_lon = maxval(lon, mask = lon >= -180.0 .and. lon <= 180.0)
      west_lon = minval(lon, mask = lon >= -180.0 .and. lon <= 180.0)
      
      idx_lon_first = int ( abs (west_lon - data_obj_meta%first_lon )  /  data_obj_meta % del_lon)
      idx_lon_last  = int ( abs (east_lon - data_obj_meta%first_lon )  /  data_obj_meta % del_lon)
      idx_lat_first = int ( abs (north_lat - data_obj_meta%first_lat ) /  data_obj_meta % del_lat)
      idx_lat_last  = int ( abs (south_lat - data_obj_meta%first_lat ) /  data_obj_meta % del_lat) 
      
      !- most files have north as first lat, if not we use this switch ( for volcano )
      if ( north_lat > data_obj_meta%first_lat ) then
         
         temp = idx_lat_first
         idx_lat_first = idx_lat_last
         idx_lat_last = temp
      end if
      
      nr_lat = idx_lat_last - idx_lat_first
      nr_lon = idx_lon_last - idx_lon_first
      
      nx = size ( lon,  1 )
      ny = size ( lon , 2 )
      
      allocate ( ilat (nx, ny ) , ilon (nx, ny ))
      allocate ( ilat_buffer (nx, ny ) , ilon_buffer (nx, ny ))
      
      ilat = int ( abs ( lat - data_obj_meta%first_lat ) / data_obj_meta % del_lat ) + 1
      ilon = int ( abs ( lon - data_obj_meta%first_lon ) / data_obj_meta % del_lon ) + 1
      
      ilat_buffer = ilat - idx_lat_first + 1
      ilon_buffer = ilon - idx_lon_first + 1
      
      
     
      deallocate ( ilon, ilat) 
      start_2d = [idx_lon_first , idx_lat_first]
      stride_2d = [ 1, 1 ]
      edge_2d = [ nr_lon  +1,  nr_lat +1 ]
      
   
   end subroutine index_in_buffer_array
   
   
   ! ====================================================================
   !
   ! ====================================================================
   subroutine read_data_i1 ( data_obj , lat , lon )
      type ( sfc_data_i1_type ) :: data_obj
      real, dimension(:,:) :: lat , lon
      integer, dimension(2) :: start_2d, stride_2d, edge_2d
      integer ( int1 ) , allocatable , dimension(:,:) :: buffer      
      integer ( int4 ) , allocatable , dimension(:,:) :: ilat_buffer , ilon_buffer
      integer  :: nx, ny
      integer :: i , j
      
      ! -executable
      call index_in_buffer_array ( data_obj%meta, lat, lon &             ! input
            , ilat_buffer , ilon_buffer , start_2d,stride_2d, edge_2d)   ! output
     
      call read_hdf_sds(data_obj%meta % file_id, trim(data_obj%meta%sds_name), start_2d, stride_2d, edge_2d, buffer )
    
      nx = size ( lon,  1 )
      ny = size ( lon , 2 )
     
      allocate ( data_obj % data ( nx , ny ) )
      
      data_obj % data  (:,:) = 0
     
      do i = 1  , nx 
         do j = 1 , ny
         
            if ( lat (i,j) == -999. .or. lon(i,j) == -999.   &
               .or. ilat_buffer (i,j) > edge_2d(2)  .or. ilon_buffer (i,j) > edge_2d(1) ) cycle
             
            data_obj % data  (i, j ) = buffer ( ilon_buffer ( i, j), ilat_buffer (i, j) )
         
         end do
      end do
   
      deallocate ( buffer )
      deallocate ( ilon_buffer,ilat_buffer)
      
   end subroutine read_data_i1
   
   ! ====================================================================
   !
   ! ====================================================================
   subroutine read_data_i2 ( data_obj , lat , lon )
      type ( sfc_data_i2_type ) :: data_obj
      real, dimension(:,:) :: lat , lon      
      integer, dimension(2) :: start_2d, stride_2d, edge_2d     
      integer ( int2 ) , allocatable , dimension(:,:) :: buffer     
      integer ( int4 ) , allocatable , dimension(:,:) :: ilat_buffer , ilon_buffer      
      integer  :: nx, ny
      integer :: i , j
      
      call index_in_buffer_array ( data_obj%meta, lat, lon &             ! input
            , ilat_buffer , ilon_buffer , start_2d,stride_2d, edge_2d)   ! output
 
      call read_hdf_sds(data_obj%meta % file_id, trim(data_obj%meta%sds_name), start_2d &
        , stride_2d, edge_2d, buffer )
     
      nx = size ( lon,  1 )
      ny = size ( lon , 2 )
      
         
      allocate ( data_obj % data ( nx , ny ) )
      data_obj % data ( :,: ) = 0
      do i = 1  , nx 
         do j = 1 , ny
            if ( lat (i,j) == -999.) cycle
            if ( lon(i,j) == -999. .or. ilat_buffer (i,j) > edge_2d(2)  .or. ilon_buffer (i,j) > edge_2d(1) ) cycle
               
            data_obj % data  (i, j ) =  buffer ( ilon_buffer ( i, j), ilat_buffer (i, j) ) 
         end do
      end do
      
      deallocate ( buffer )
      deallocate ( ilon_buffer,ilat_buffer )
      

   end subroutine read_data_i2
   
   ! ====================================================================
   !
   ! ====================================================================
   subroutine read_modis_white_sky ( data_obj , lat , lon )
      type ( sfc_data_r4_type ) :: data_obj
      real, dimension(:,:) :: lat , lon    
      integer, dimension(2) :: start_2d, stride_2d, edge_2d
      integer ( kind = int2)  , allocatable , dimension(:,:) :: buffer
      integer ( int4 ) , allocatable , dimension(:,:) :: ilat_buffer , ilon_buffer
      integer  :: nx, ny
      integer :: i , j
      
      
      call index_in_buffer_array ( data_obj%meta, lat, lon &             ! input
            , ilat_buffer , ilon_buffer , start_2d,stride_2d, edge_2d)   ! output
         
      call read_hdf_sds( &
            data_obj%meta % file_id &
         , trim(data_obj%meta%sds_name) &
         , start_2d &
         , stride_2d &
         , edge_2d &
         , buffer )
    
      nx = size ( lon , 1 )
      ny = size ( lon , 2 )
              
      allocate ( data_obj % data ( nx , ny ) )
     
      do i = 1  , nx 
         do j = 1 , ny
            ! -- factor computes white sky albedo in a range from [0,100]
           
            if ( lat (i,j) == -999. .or. lon(i,j) == -999.  .or. ilat_buffer (i,j) > edge_2d(2)  .or. ilon_buffer (i,j) > edge_2d(1)  ) cycle
         
            data_obj % data  (i, j ) = buffer ( ilon_buffer ( i, j), ilat_buffer (i, j) ) 
            if ( buffer ( ilon_buffer ( i, j), ilat_buffer (i, j) ) == 32767 ) data_obj % data  (i, j ) = -999.
      
         end do
      end do
      
 
      where (  data_obj % data  /= -999. )
         data_obj % data = 0.1 * data_obj % data
      end where
      
      deallocate ( buffer )
      deallocate (ilon_buffer,ilat_buffer)
   
   
   end subroutine read_modis_white_sky
   
   ! ====================================================================
   !
   ! ====================================================================  
   
   subroutine read_emis ( data_obj , lat , lon )
      type ( sfc_data_r4_type ) :: data_obj
      real, dimension(:,:) :: lat , lon
      integer(int4), dimension(2) :: start_2d, stride_2d, edge_2d
      real ( kind = real4)  , allocatable , dimension(:,:) :: buffer
      integer ( int4 ) , allocatable , dimension(:,:) :: ilat_buffer , ilon_buffer
      character(len =30) ::   name_scale , name_offset 
      integer  :: nx, ny
      integer :: i , j
       
      call index_in_buffer_array ( data_obj%meta, lat, lon &             ! input
            , ilat_buffer , ilon_buffer , start_2d,stride_2d, edge_2d)   ! output
   
      name_scale = 'scale_factor'
      name_offset = 'add_offset'
      
      call read_hdf_sds(data_obj%meta % file_id, trim(data_obj%meta%sds_name), start_2d, stride_2d, &
               &  edge_2d, name_scale, name_offset, buffer )
      
      nx = size ( lon, 1 )
      ny = size ( lon , 2)

      allocate ( data_obj % data ( nx , ny ) )
      
      do i = 1  , nx 
         do j = 1 , ny
            if ( lat (i,j) == -999. .or. lon(i,j) == -999.  .or. ilat_buffer (i,j) > edge_2d(2) .or. ilon_buffer (i,j) > edge_2d(1)) cycle
            ! -- factor computes white sky albedo in a range from [0,100]
            data_obj % data  (i, j ) =   buffer ( ilon_buffer ( i, j), ilat_buffer (i, j) )
         end do
      end do
 
      deallocate ( buffer )
      deallocate ( ilon_buffer,ilat_buffer)   
   
   end subroutine read_emis
 
   !============================== = = = = = = = = = = = = = = 
   !
   !============================================================
   
   subroutine update_with_nwp ( this , nwp, geo )
      use cx_nwp_mod, only:  nwp_main_type
      use cx_geo_mod, only:  geo_type
      implicit none   
      class ( sfc_main_type ) , intent ( in out ) :: this
      type ( nwp_main_type ) , intent ( in ) :: nwp
      type ( geo_type ) , intent ( in ) :: geo
      
      !- locals
      integer :: xnwp
      integer :: ynwp
      integer :: i , j
      
	
      
      if ( .not. allocated( this % z % data ) ) &
         & allocate ( this % z % data ( this % dim_x,  this % dim_y ) )
      
      if ( .not. this % snow_class % is_set ) then
         allocate ( this % snow_class % data ( this % dim_x,  this % dim_y ) )

         this % snow_class % data = et_snow_class % NO_SNOW_NOR_ICE
        
         do j = 1, this % dim_y 
            do i = 1, this % dim_x
           
               
					
               xnwp = geo % idx_nwp_x ( i , j)
               ynwp = geo % idx_nwp_y ( i , j)
            
			
               this % z % data ( i , j ) = nwp % zsfc ( xnwp , ynwp )
              
               if ( nwp % weasd ( xnwp , ynwp ) > 0.1 ) then
                  this % snow_class % data ( i, j) = et_snow_class % SNOW
               end if    
					
               if ( nwp % ice ( xnwp , ynwp ) > 0.5 ) then
                  this % snow_class % data ( i, j) = et_snow_class % SEA_ICE
               end if 
            end do
         end do   
         
     
      end if
            
      ! - update emis to default values for snow
      
      where ( this % snow_class % data == et_snow_class % SNOW )
         this % emis ( 20 ) % data = 0.984
         this % emis ( 29 ) % data = 0.979
         this % emis ( 31 ) % data = 0.979
         this % emis ( 32 ) % data = 0.977
      end where
		
	
         
   end subroutine update_with_nwp
   
   
   ! ============================== = = = = = = = = = = = = = = =======
   ! 
   ! ============================== = = = = = = = = = = = = = = =======
   
   subroutine deallocate_sfc ( this ) 
      class ( sfc_main_type ) :: this
      integer :: i_w_sky
      integer :: i_emis
       
      if ( allocated ( this % land_class % data ) ) deallocate  ( this % land_class % data )
      if ( allocated ( this % sfc_type % data ) ) deallocate  ( this % sfc_type % data )
      if ( allocated ( this % coast_mask % data ) ) deallocate  ( this % coast_mask % data )
      if ( allocated ( this % snow_class % data ) ) deallocate  ( this % snow_class % data )
      if ( allocated ( this % elevation % data ) ) deallocate  ( this % elevation % data )
      if ( allocated ( this % z % data ) ) deallocate  ( this % z % data )
      if ( allocated ( this % volcano % data )) deallocate ( this % volcano % data )
      if ( allocated ( this % sst % data )) deallocate ( this % sst % data )
      if ( allocated (this % sst % stdv) ) deallocate ( this % sst % stdv )
      if ( allocated (this % sea_ice_fraction % data) ) deallocate ( this % sea_ice_fraction % data )
      do i_w_sky = 1 , 7
         if ( allocated ( this % modis_w_sky (i_w_sky) % data )) &
            &  deallocate ( this % modis_w_sky (i_w_sky) % data )
      end do
      
      do i_emis = 20 , 36
         if ( allocated ( this % emis (i_emis) % data )) deallocate ( this % emis (i_emis) % data )
      end do 
      
   end subroutine deallocate_sfc
  
  
  
   
end module cx_sfc_mod



!program test_it
!!use sfc_data_mod!



